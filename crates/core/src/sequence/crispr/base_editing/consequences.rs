//! Optional ordinary-codon consequences on one explicitly assumed exact CDS.
//! Compound transcript paths and initiation-context effects are not inferred.
use super::*;
use crate::construct_qc::{AnnotationAssumption, annotation_revision};
use collagene_model::{Location, Sequence};
use std::collections::BTreeMap;

dto! {pub struct CdsRequest {pub annotation: AnnotationAssumption, pub genetic_code: u8}}
impl CdsRequest {
    pub(super) fn validate(&self, molecule: &Molecule) -> Result<(), Error> {
        if self.annotation.feature_index >= molecule.features.len()
            || !hex(&self.annotation.source_revision)
            || !hex(&self.annotation.annotation_revision)
            || !text(&self.annotation.rationale, 1024)
        {
            return Err(bad("invalid explicit CDS-role assumption"));
        }
        translation::CodeDescriptor::for_id(self.genetic_code)?;
        Ok(())
    }
}
dto! {
    pub struct CdsEvidence {
        pub feature_index: usize,
        pub genetic_code: u8,
        /// Means coordinate/register availability under the explicit assumption;
        /// it does not verify the annotation's biological role or every codon call.
        pub available: bool,
        #[serde(deserialize_with = "crate::encoding::required_option")]
        #[cfg_attr(feature = "json-schema", schemars(required))]
        pub code_revision: Option<String>,
        pub reason: String,
    }
}
dto! {pub enum ConsequenceKind {Synonymous, Missense, StopGained, StopLost, StopRetained, Unresolved}}
dto! {
    pub struct CodonConsequence {
        pub feature_index: usize,
        /// One-based ordinal in the explicitly assumed biological CDS frame.
        pub codon_number: usize,
        #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, usize, 3>")]
        pub reference_positions: Vec<usize>,
        pub before: String,
        pub after: String,
        pub before_amino_acid: char,
        pub after_amino_acid: char,
        pub kind: ConsequenceKind,
        pub at_cds_start: bool,
    }
}
struct Codon {
    index: usize,
    bases: [u8; 3],
    positions: [usize; 3],
    amino_acid: char,
}
pub(super) struct CodingContext {
    feature_index: usize,
    genetic_code: u8,
    reverse: bool,
    codons: Vec<Codon>,
}
fn interval(location: &Location) -> Option<(usize, usize, bool)> {
    match location {
        Location::Range {
            start,
            end,
            before: false,
            after: false,
        } => Some((*start, *end, false)),
        Location::Complement(inner) => match inner.as_ref() {
            Location::Range {
                start,
                end,
                before: false,
                after: false,
            } => Some((*start, *end, true)),
            _ => None,
        },
        _ => None,
    }
}
fn ordinary(sequence: &str, code: u8) -> Result<String, Error> {
    Ok(translation::translate(
        &Sequence::dna(sequence)?,
        &translation::TranslationOptions {
            genetic_code: code,
            ..Default::default()
        },
    )?
    .protein)
}
fn code_revision(code: u8) -> Result<String, Error> {
    let mut bases = String::with_capacity(192);
    for a in b"TCAG" {
        for b in b"TCAG" {
            for c in b"TCAG" {
                bases.extend([*a as char, *b as char, *c as char]);
            }
        }
    }
    hash(&(
        translation::CodeDescriptor::for_id(code)?,
        ordinary(&bases, code)?,
    ))
}
pub(super) fn prepare(
    molecule: &Molecule,
    request: Option<&CdsRequest>,
    sites: &[EditableSite],
) -> Result<(Option<CdsEvidence>, Option<CodingContext>), Error> {
    let Some(request) = request else {
        return Ok((None, None));
    };
    let index = request.annotation.feature_index;
    let feature = &molecule.features[index];
    let mut evidence = CdsEvidence {feature_index: index, genetic_code: request.genetic_code, available: false,
        code_revision: None,
        reason: "Codon consequences need a current explicit assumed CDS with exact single-interval biological coordinates, complete frame and matching genetic code. Compound, fuzzy, review-required and unsupported exception paths remain unresolved.".into()};
    let current_annotation = match annotation_revision(feature) {
        Ok(revision) => revision,
        Err(_) => {
            evidence.reason = "The selected annotation exceeds the canonical explicit-assumption identity bounds; ordinary codon consequences remain unresolved while finite window states can still be inspected.".into();
            return Ok((Some(evidence), None));
        }
    };
    if request.annotation.source_revision != molecule.revision()
        || request.annotation.annotation_revision != current_annotation
        || feature.needs_review
        || !feature.kind.eq_ignore_ascii_case("CDS")
    {
        return Ok((Some(evidence), None));
    }
    let Some((start, end, reverse)) = interval(&feature.location) else {
        return Ok((Some(evidence), None));
    };
    if end - start > 1_000_000 {
        evidence.reason = "The selected CDS exceeds the one-million-base consequence extraction bound; finite window states remain available and codon consequences are unresolved.".into();
        return Ok((Some(evidence), None));
    }
    let translated = match translation::translate_cds(molecule, feature) {
        Ok(translated) => translated,
        Err(error) => {
            evidence.reason = format!(
                "The canonical CDS translator cannot resolve this declared annotation: {error}"
            );
            return Ok((Some(evidence), None));
        }
    };
    if translated.codon_start != 1
        || translated.translation.trailing_bases != 0
        || translated.translation.options.genetic_code != request.genetic_code
        || feature.qualifiers.iter().any(|(key, _)| key == "partial")
    {
        return Ok((Some(evidence), None));
    }
    let extracted = feature.location.extract(&molecule.sequence)?;
    let mut indices = BTreeSet::new();
    for site in sites {
        if (start..end).contains(&site.reference_position) {
            let offset = if reverse {
                end - 1 - site.reference_position
            } else {
                site.reference_position - start
            };
            indices.insert(offset / 3);
        }
    }
    let mut codons = Vec::with_capacity(indices.len());
    let mut bundled = String::with_capacity(indices.len() * 3);
    for index in indices {
        let offset = index * 3;
        let bases = extracted
            .as_bytes()
            .get(offset..offset + 3)
            .ok_or_else(|| bad("incomplete ordinary codon at declared CDS boundary"))?;
        let bases = [bases[0], bases[1], bases[2]];
        bundled.extend(bases.map(char::from));
        let positions = [offset, offset + 1, offset + 2].map(|offset| {
            if reverse {
                end - 1 - offset
            } else {
                start + offset
            }
        });
        codons.push(Codon {
            index,
            bases,
            positions,
            amino_acid: 'X',
        });
    }
    if !codons.is_empty() {
        let amino = ordinary(&bundled, request.genetic_code)?;
        for (codon, amino) in codons.iter_mut().zip(amino.chars()) {
            codon.amino_acid = amino;
        }
    }
    evidence.available = true;
    evidence.code_revision = Some(code_revision(request.genetic_code)?);
    evidence.reason = "Ordinary codon calls use the canonical genetic-code owner on the explicitly assumed CDS strand/register. No annotation truth, actual initiation, protein function or expression is inferred; at_cds_start flags the ordinary first-codon comparison without initiating-M recoding.".into();
    Ok((
        Some(evidence),
        Some(CodingContext {
            feature_index: index,
            genetic_code: request.genetic_code,
            reverse,
            codons,
        }),
    ))
}
impl CodingContext {
    pub(super) fn assess(
        &self,
        selected: &[&EditableSite],
    ) -> Result<Vec<CodonConsequence>, Error> {
        let substitutions = selected
            .iter()
            .map(|site| (site.reference_position, site.reference_to as u8))
            .collect::<BTreeMap<_, _>>();
        let mut changed = Vec::new();
        let mut bundle = String::new();
        for codon in &self.codons {
            let mut after = codon.bases;
            for (offset, position) in codon.positions.iter().enumerate() {
                if let Some(base) = substitutions.get(position) {
                    after[offset] = if self.reverse {
                        collagene_model::complement(*base, Alphabet::Dna)
                    } else {
                        *base
                    };
                }
            }
            if after == codon.bases {
                continue;
            }
            bundle.extend(after.map(char::from));
            changed.push((codon, after));
        }
        if changed.is_empty() {
            return Ok(vec![]);
        }
        // One small translation per reference state, rather than one independent
        // decoder/cache allocation for each individual affected codon.
        let amino = ordinary(&bundle, self.genetic_code)?;
        changed
            .into_iter()
            .zip(amino.chars())
            .map(|((codon, after), after_amino_acid)| {
                let before_amino_acid = codon.amino_acid;
                let kind = if before_amino_acid == 'X' || after_amino_acid == 'X' {
                    ConsequenceKind::Unresolved
                } else if before_amino_acid == '*' && after_amino_acid == '*' {
                    ConsequenceKind::StopRetained
                } else if after_amino_acid == '*' {
                    ConsequenceKind::StopGained
                } else if before_amino_acid == '*' {
                    ConsequenceKind::StopLost
                } else if before_amino_acid == after_amino_acid {
                    ConsequenceKind::Synonymous
                } else {
                    ConsequenceKind::Missense
                };
                Ok(CodonConsequence {
                    feature_index: self.feature_index,
                    codon_number: codon.index + 1,
                    reference_positions: codon.positions.to_vec(),
                    before: String::from_utf8(codon.bases.to_vec())
                        .map_err(|_| bad("non-ASCII codon"))?,
                    after: String::from_utf8(after.to_vec()).map_err(|_| bad("non-ASCII codon"))?,
                    before_amino_acid,
                    after_amino_acid,
                    kind,
                    at_cds_start: codon.index == 0,
                })
            })
            .collect()
    }
}
pub(super) fn validate_report(report: &Report) -> Result<(), Error> {
    if report.cds.is_some() != report.request.cds.is_some() {
        return Err(bad("CDS consequence request/evidence shape differs"));
    }
    if let Some(cds) = &report.cds {
        translation::CodeDescriptor::for_id(cds.genetic_code)?;
        if report.request.cds.as_ref().is_none_or(|request| {
            request.annotation.feature_index != cds.feature_index
                || request.genetic_code != cds.genetic_code
        }) || !text(&cds.reason, 4096)
            || cds.available != cds.code_revision.is_some()
            || cds
                .code_revision
                .as_ref()
                .is_some_and(|revision| !hex(revision))
        {
            return Err(bad("invalid captured CDS consequence scope"));
        }
    }
    for state in &report.states {
        for codon in &state.codon_consequences {
            if codon.codon_number == 0
                || codon.reference_positions.len() != 3
                || codon
                    .reference_positions
                    .iter()
                    .any(|position| *position >= report.source_length)
                || codon.before.len() != 3
                || codon.after.len() != 3
                || Sequence::dna(&codon.before).is_err()
                || Sequence::dna(&codon.after).is_err()
                || !codon.before_amino_acid.is_ascii()
                || !codon.after_amino_acid.is_ascii()
                || !b"ACDEFGHIKLMNPQRSTVWYX*".contains(&(codon.before_amino_acid as u8))
                || !b"ACDEFGHIKLMNPQRSTVWYX*".contains(&(codon.after_amino_acid as u8))
                || report
                    .cds
                    .as_ref()
                    .is_none_or(|cds| !cds.available || cds.feature_index != codon.feature_index)
            {
                return Err(bad("invalid ordinary codon consequence evidence"));
            }
        }
    }
    Ok(())
}
