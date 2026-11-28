//! Explicit ideal hydrolytic 3′→5′ shortening of free duplex strand ends.
//! No enzyme identity, kinetic cut-length distribution, synthesis or ligation.
use crate::chemistry::{Duplex, Strand, StrandMapping, ThreePrime};
use collagene_model::{Error, Sequence, Topology};
use serde::{Deserialize, Serialize};
use std::io::Write;
pub const ENGINE_VERSION: &str = "three-prime-resection-1";
pub const MAX_SOURCE_BASES: usize = 2_000_000;
pub const MAX_REPORT_BYTES: usize = 16 << 20;
macro_rules! dto {($($x:tt)*)=>{#[derive(Debug,Clone,PartialEq,Eq,Serialize,Deserialize)]#[cfg_attr(feature="json-schema",derive(schemars::JsonSchema))]#[serde(deny_unknown_fields)]$($x)*};}
// Counts in the supplied oriented duplex. Both remove stored-strand suffixes:
// left shortens bottom3′; right shortens top3′. Original5′ ends never change.
dto! {pub struct Cuts {pub left:usize,pub right:usize}}
dto! {pub struct Removed {pub strand:Strand,pub start:usize,pub end:usize}}
dto! {pub struct Report {
    pub engine_version:String,pub source_geometry_revision:String,pub cuts:Cuts,
    pub product:Duplex,pub mappings:Vec<StrandMapping>,pub removed:Vec<Removed>,
    pub warnings:Vec<String>,
}}
fn bad(message: &str) -> Error {
    Error::Unsupported(format!("3′ resection: {message}"))
}
struct HashWriter {
    count: usize,
    hash: blake3::Hasher,
}
impl Write for HashWriter {
    fn write(&mut self, b: &[u8]) -> std::io::Result<usize> {
        self.count = self
            .count
            .checked_add(b.len())
            .filter(|n| *n <= MAX_REPORT_BYTES)
            .ok_or_else(|| std::io::Error::other("resection size limit"))?;
        self.hash.update(b);
        Ok(b.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}
pub(crate) fn revision(source: &impl Serialize) -> Result<String, Error> {
    let mut writer = HashWriter {
        count: 0,
        hash: blake3::Hasher::new(),
    };
    serde_json::to_writer(&mut writer, source).map_err(|_| Error::TooLarge)?;
    Ok(writer.hash.finalize().to_hex().to_string())
}
// Checks lengths before source traversal/allocation. Nicked substrates are
// excluded: their additional accessible3′ termini need a separate model.
pub fn preflight(source: &Duplex, cuts: &Cuts) -> Result<(), Error> {
    if source
        .top
        .len()
        .checked_add(source.bottom.len())
        .is_none_or(|n| n > MAX_SOURCE_BASES)
    {
        return Err(Error::TooLarge);
    }
    if source.topology != Topology::Linear || !source.nicks.is_empty() {
        return Err(bad(
            "requires nick-free linear duplex; circular/internal-nick processing is not declared",
        ));
    }
    if cuts.right >= source.top.len() || cuts.left >= source.bottom.len() {
        return Err(bad("resection cannot exhaust either strand"));
    }
    if (cuts.left > 0 && source.left.three_prime != ThreePrime::Hydroxyl)
        || (cuts.right > 0 && source.right.three_prime != ThreePrime::Hydroxyl)
    {
        return Err(bad(
            "positive cuts require a known accessible3′ hydroxyl; no unblocking is assumed",
        ));
    }
    let bottom_start = source
        .bottom_start
        .checked_add(cuts.left as i64)
        .ok_or(Error::TooLarge)?;
    let bottom_end = source
        .bottom_start
        .checked_add(source.bottom.len() as i64)
        .ok_or(Error::TooLarge)?;
    if bottom_start.max(0) >= bottom_end.min((source.top.len() - cuts.right) as i64) {
        return Err(bad("cuts remove the last paired core"));
    }
    source.validate()?;
    Ok(())
}
pub fn geometry(source: &Duplex, cuts: &Cuts) -> Result<Report, Error> {
    preflight(source, cuts)?;
    let top_end = source.top.len() - cuts.right;
    let bottom_end = source.bottom.len() - cuts.left;
    let mut product = Duplex {
        top: Sequence::dna(&source.top.as_str()[..top_end])?,
        bottom: Sequence::dna(&source.bottom.as_str()[..bottom_end])?,
        bottom_start: source.bottom_start + cuts.left as i64,
        topology: Topology::Linear,
        left: source.left,
        right: source.right,
        nicks: vec![],
    };
    if cuts.left > 0 {
        product.left.three_prime = ThreePrime::Hydroxyl;
    }
    if cuts.right > 0 {
        product.right.three_prime = ThreePrime::Hydroxyl;
    }
    product.validate()?;
    let mut removed = vec![];
    for (strand, start, end) in [
        (Strand::Top, top_end, source.top.len()),
        (Strand::Bottom, bottom_end, source.bottom.len()),
    ] {
        if start < end {
            removed.push(Removed { strand, start, end });
        }
    }
    let mappings = [(Strand::Top, top_end), (Strand::Bottom, bottom_end)]
        .into_iter()
        .map(|(strand, end)| StrandMapping {
            source_strand: strand,
            source_start: 0,
            source_end: end,
            product_strand: strand,
            product_start: 0,
            product_end: end,
            reverse: false,
            synthesized: false,
        })
        .collect();
    Ok(Report{engine_version:ENGINE_VERSION.into(),source_geometry_revision:revision(source)?,cuts:cuts.clone(),product,mappings,removed,warnings:vec!["Declared hydrolytic3′ shortening only. Cut lengths are explicit assumptions, not an enzyme/kit prediction. Original5′ phosphorylation is unchanged; no polymerase, ligase, host repair or blocked-end removal is assumed.".into()]})
}
// Fresh source-bound replay, never restoration of authority from a report.
pub fn replay(source: &Duplex, cuts: &Cuts, report: &Report) -> Result<(), Error> {
    if geometry(source, cuts)? != *report {
        return Err(bad("report differs from fresh source-bound resection"));
    }
    Ok(())
}

impl Report {
    /// Bounded internal documentary correspondence only. This neither resolves
    /// an original source nor reconstructs private material from imported data.
    pub fn validate_documentary(&self) -> Result<(), Error> {
        let p = &self.product;
        if self.engine_version != ENGINE_VERSION
            || self.source_geometry_revision.len() != 64
            || !self
                .source_geometry_revision
                .bytes()
                .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
        {
            return Err(bad("report engine/source geometry identity"));
        }
        let top = p
            .top
            .len()
            .checked_add(self.cuts.right)
            .ok_or(Error::TooLarge)?;
        let bottom = p
            .bottom
            .len()
            .checked_add(self.cuts.left)
            .ok_or(Error::TooLarge)?;
        if top.checked_add(bottom).is_none_or(|n| n > MAX_SOURCE_BASES)
            || self.warnings.len() > 1000
            || self.warnings.iter().any(|s| s.len() > 4096)
        {
            return Err(Error::TooLarge);
        }
        if p.topology != Topology::Linear
            || !p.nicks.is_empty()
            || (self.cuts.left > 0 && p.left.three_prime != ThreePrime::Hydroxyl)
            || (self.cuts.right > 0 && p.right.three_prime != ThreePrime::Hydroxyl)
        {
            return Err(bad("documentary product topology/nick/end-group mismatch"));
        }
        p.validate()?;
        if self.mappings.len() != 2
            || self.removed.len()
                != usize::from(self.cuts.right > 0) + usize::from(self.cuts.left > 0)
        {
            return Err(bad("copy/removal mapping count differs"));
        }
        for (i, (strand, end)) in [(Strand::Top, p.top.len()), (Strand::Bottom, p.bottom.len())]
            .into_iter()
            .enumerate()
        {
            let m = &self.mappings[i];
            if m.source_strand != strand
                || m.product_strand != strand
                || m.source_start != 0
                || m.source_end != end
                || m.product_start != 0
                || m.product_end != end
                || m.reverse
                || m.synthesized
            {
                return Err(bad("report claims non-copy strand ancestry"));
            }
        }
        let expected = [
            (Strand::Top, p.top.len(), top),
            (Strand::Bottom, p.bottom.len(), bottom),
        ]
        .into_iter()
        .filter(|(_, a, b)| a < b)
        .map(|(strand, start, end)| Removed { strand, start, end })
        .collect::<Vec<_>>();
        if self.removed != expected {
            return Err(bad("removed spans differ from explicit cuts"));
        }
        revision(self)?;
        Ok(())
    }
}
