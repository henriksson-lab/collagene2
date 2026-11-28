//! Flat declaration of two PCR branches, reciprocal extension, then outer PCR.
//! This module compiles supplied oligos; it never derives material from a target.
use crate::{
    analysis::Primer,
    overlap_extension,
    verification::{Operation, OrientedSource, Source},
};
use collagene_model::{Alphabet, Error};
use serde::{Deserialize, Serialize};

#[derive(Clone, Debug, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct PrimerPair {
    pub forward: Primer,
    pub reverse: Primer,
    pub template_strand: crate::chemistry::Strand,
}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Request {
    pub left: PrimerPair,
    pub right: PrimerPair,
    pub outer: PrimerPair,
    pub extension: overlap_extension::Options,
}
fn bad(s: &str) -> Error {
    Error::Format(format!("overlap-extension PCR preparation: {s}"))
}
pub(crate) fn primer(p: &Primer) -> Result<(), Error> {
    if p.name.is_empty()
        || p.name.len() > 256
        || p.name.chars().any(char::is_control)
        || p.annealing.alphabet() != Alphabet::Dna
        || p.annealing.is_empty()
        || p.annealing.len() > 1000
        || p.tail.len() > 1000
        || !p
            .annealing
            .as_str()
            .bytes()
            .chain(p.tail.bytes())
            .all(|b| b"ACGT".contains(&b))
    {
        return Err(bad(
            "supply concrete DNA oligos with annealing1..1000, tail0..1000 and bounded names",
        ));
    }
    Ok(())
}
pub(crate) fn extension_shape(o: &overlap_extension::Options) -> Result<(), Error> {
    if o.min_paired_bases == 0
        || o.overlap < o.min_paired_bases
        || o.overlap > overlap_extension::MAX_SOURCE_BASES / 2
        || o.max_comparisons == 0
        || o.max_comparisons > overlap_extension::MAX_COMPARISONS
        || o.overlap > o.max_comparisons
        || o.max_candidates == 0
        || o.max_candidates > overlap_extension::MAX_CANDIDATES
    {
        return Err(bad("invalid selected terminal pairing or search bounds"));
    }
    Ok(())
}
impl Request {
    /// Borrowed shape only. Binding sites, actual homology and physical state
    /// remain the native verifier's work, not accepted preview evidence.
    pub fn validate_shape(&self) -> Result<(), Error> {
        for pair in [&self.left, &self.right, &self.outer] {
            primer(&pair.forward)?;
            primer(&pair.reverse)?;
        }
        if self.outer.template_strand != crate::chemistry::Strand::Top {
            return Err(bad(
                "outer PCR uses the newly synthesized extension's explicit Top frame",
            ));
        }
        extension_shape(&self.extension)
    }
    pub fn validate(&self) -> Result<(), Error> {
        self.validate_shape()
    }
    /// Conservative original-source work/material bounds; does no sequence
    /// search. Four stored source strands and both PCR branches are accounted.
    pub(crate) fn work(&self, source_bases: usize) -> Result<(usize, usize), Error> {
        self.validate_shape()?;
        let length = |pair: &PrimerPair| {
            source_bases
                .checked_add(pair.forward.tail.len())
                .and_then(|n| n.checked_add(pair.reverse.tail.len()))
                .ok_or(Error::TooLarge)
        };
        let a = length(&self.left)?;
        let b = length(&self.right)?;
        let sum = a.checked_add(b).ok_or(Error::TooLarge)?;
        if source_bases == 0
            || sum
                .checked_mul(2)
                .is_none_or(|n| n > overlap_extension::MAX_SOURCE_BASES)
            || self.extension.overlap >= a.min(b)
        {
            return Err(bad(
                "source/PCR upper bounds cannot support the declared reciprocal extension",
            ));
        }
        let extended = sum - self.extension.overlap;
        let pcr_work = |n: usize, pair: &PrimerPair| {
            n.checked_mul(pair.forward.annealing.len() + pair.reverse.annealing.len())
                .and_then(|n| n.checked_mul(4))
                .ok_or(Error::TooLarge)
        };
        let scan = sum
            .checked_mul(3)
            .and_then(|n| n.checked_add(self.extension.overlap))
            .ok_or(Error::TooLarge)?
            .min(self.extension.max_comparisons);
        let work = pcr_work(source_bases, &self.left)?
            .checked_add(pcr_work(source_bases, &self.right)?)
            .and_then(|n| n.checked_add(scan))
            .and_then(|n| n.checked_add(pcr_work(extended, &self.outer).ok()?))
            .ok_or(Error::TooLarge)?;
        // The source's initial reservation is handled by the shared caller.
        let extra = sum
            .checked_add(extended)
            .and_then(|n| n.checked_add(extended))
            .and_then(|n| n.checked_add(self.outer.forward.tail.len()))
            .and_then(|n| n.checked_add(self.outer.reverse.tail.len()))
            .ok_or(Error::TooLarge)?;
        Ok((work, extra))
    }
    pub(crate) fn append(
        &self,
        source: Source,
        steps: &mut Vec<Operation>,
    ) -> Result<Source, Error> {
        self.validate_shape()?;
        let first = steps.len();
        let right = first.checked_add(1).ok_or(Error::TooLarge)?;
        let extension = first.checked_add(2).ok_or(Error::TooLarge)?;
        let outer = first.checked_add(3).ok_or(Error::TooLarge)?;
        if outer >= 1000 {
            return Err(Error::TooLarge);
        }
        steps.push(Operation::Pcr {
            template: source.clone(),
            forward: self.left.forward.clone(),
            reverse: self.left.reverse.clone(),
            template_strand: self.left.template_strand,
        });
        steps.push(Operation::Pcr {
            template: source,
            forward: self.right.forward.clone(),
            reverse: self.right.reverse.clone(),
            template_strand: self.right.template_strand,
        });
        steps.push(Operation::OverlapExtension {
            left: OrientedSource {
                source: Source::Step(first),
                reverse: false,
            },
            right: OrientedSource {
                source: Source::Step(right),
                reverse: false,
            },
            options: self.extension.clone(),
        });
        steps.push(Operation::Pcr {
            template: Source::Step(extension),
            forward: self.outer.forward.clone(),
            reverse: self.outer.reverse.clone(),
            template_strand: self.outer.template_strand,
        });
        Ok(Source::Step(outer))
    }
}
