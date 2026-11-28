//! Positional modification state on the actual stored 5′→3′ strands.
//! Explicit source-generation/synthesis assumptions, never imported metadata authority.
use crate::{
    chemistry::{
        self, Duplex, EndRepairMethod, EndTreatment, Strand,
        methylation::{self, MethylationState, Outcome, Sensitivity},
    },
    cloning::ligation::tailing::{self, TailingBase},
    sequence::analysis::{self, Enzyme},
};
use collagene_model::{Error, Molecule, Sequence, Topology};
use serde::{Deserialize, Serialize};
use std::io::Write;
pub const ENGINE_VERSION: &str = "modified-duplex-1";
pub const PROFILE_VERSION: &str = "positional-catalog-1";
const MAX_BASES: usize = 1_000_000;
const MAX_SPANS: usize = 16_384;
const MAX_BYTES: usize = 32 * 1024 * 1024;
macro_rules! dto{($($x:tt)*)=>{#[derive(Debug,Clone,PartialEq,Eq,Serialize,Deserialize)]#[cfg_attr(feature="json-schema",derive(schemars::JsonSchema))]#[serde(deny_unknown_fields)]$($x)*};}
macro_rules! en{($($x:tt)*)=>{#[derive(Debug,Clone,Copy,PartialEq,Eq,PartialOrd,Ord,Serialize,Deserialize)]#[cfg_attr(feature="json-schema",derive(schemars::JsonSchema))]$($x)*};}
en! {pub enum Chemical{N6A,C5C}}
en! {pub enum CallState{Unknown,Unmodified,Modified}}
dto! {pub struct ModificationSpan{pub chemical:Chemical,pub start:usize,pub end:usize,pub state:CallState}}
dto! {pub struct StrandModifications{pub n6a_default:CallState,pub c5c_default:CallState,pub spans:Vec<ModificationSpan>}}
dto! {pub struct GenerationAssumptions{pub dam:MethylationState,pub dcm:MethylationState,pub cpg:MethylationState}}
impl Default for GenerationAssumptions {
    fn default() -> Self {
        Self {
            dam: MethylationState::Unknown,
            dcm: MethylationState::Unknown,
            cpg: MethylationState::Unknown,
        }
    }
}
dto! {pub struct GenerationRecord{pub source_revision:String,pub assumptions:GenerationAssumptions}}
dto! {pub struct TrackedDuplex{pub duplex:Duplex,pub top_modifications:StrandModifications,pub bottom_modifications:StrandModifications,pub generation:Vec<GenerationRecord>}}
dto! {pub struct TrackedDigestFragment{pub duplex:TrackedDuplex,pub top_start:usize,pub bottom_start:i64}}
dto! {pub struct CheckedDigestOptions{pub selection:crate::cloning::restriction::catalog::CatalogSelection,pub max_comparisons:usize,pub max_sites:usize}}
/// Current checked-cleavage model defaults. These are authored limits, never a
/// clamp or conversion of an operator's explicit request.
pub fn checked_inspection_options() -> crate::cloning::restriction::inspection::InspectionOptions {
    crate::cloning::restriction::inspection::InspectionOptions {
        max_comparisons: 50_000_000,
        max_sites: 4096,
        max_fragments: 4097,
        max_report_bytes: 4 * 1024 * 1024,
    }
}
impl Default for CheckedDigestOptions {
    fn default() -> Self {
        Self {
            selection: checked_selection("EcoRI").expect("embedded checked enzyme"),
            max_comparisons: 20_000_000,
            max_sites: 4096,
        }
    }
}
dto! {pub struct CheckedSite{pub enzyme:String,pub recognition_start:usize,pub reverse:bool,pub top_cut:Option<usize>,pub bottom_cut:Option<usize>,pub outcome:Outcome,pub evidence:Vec<PositionEvidence>,pub reason:String}}
dto! {pub struct PositionEvidence{pub strand:Strand,pub index:usize,pub chemical:Chemical,pub state:CallState}}
dto! {pub struct CheckedDigestReport{pub engine_version:String,pub profile_version:String,pub profile_revision:String,pub source_revision:String,pub options:CheckedDigestOptions,pub complete:bool,pub products_determined:bool,pub comparisons:usize,pub sites:Vec<CheckedSite>,pub fragments:Vec<TrackedDigestFragment>,pub warnings:Vec<String>,pub limits:Vec<String>}}
dto! {
/// Bounded documentary cleavage evidence; no fragment DNA or private state.
pub struct CheckedCleavageEvidence{pub engine_version:String,pub profile_version:String,pub profile_revision:String,pub source_revision:String,pub options:CheckedDigestOptions,pub complete:bool,pub products_determined:bool,pub comparisons:usize,pub sites:Vec<CheckedSite>,pub warnings:Vec<String>,pub limits:Vec<String>}}
/// One enclosing-plan ledger shared by its preparation, source and recut scans.
/// Nested saved-origin replay remains a separate bounded native scope.
#[derive(Debug)]
pub struct CleavageBudget {
    limit: usize,
    used: usize,
}
impl CleavageBudget {
    pub fn new(limit: usize) -> Result<Self, Error> {
        if limit == 0 || limit > 100_000_000 {
            return Err(Error::TooLarge);
        }
        Ok(Self { limit, used: 0 })
    }
    pub fn remaining(&self) -> usize {
        self.limit - self.used
    }
    pub fn used(&self) -> usize {
        self.used
    }
    pub fn charge(&mut self, comparisons: usize) -> Result<(), Error> {
        self.used = self
            .used
            .checked_add(comparisons)
            .filter(|n| *n <= self.limit)
            .ok_or(Error::TooLarge)?;
        Ok(())
    }
}
/// Exact canonical identity, never an alias, isoschizomer or formulation lookup.
pub fn checked_selection(
    name: &str,
) -> Result<crate::cloning::restriction::catalog::CatalogSelection, Error> {
    let c = crate::cloning::restriction::catalog::builtin_catalog()?;
    let e = c
        .edition()
        .enzymes
        .iter()
        .find(|e| e.name == name)
        .ok_or_else(|| bad("exact embedded enzyme identity absent"))?;
    crate::cloning::restriction::catalog::select(&c, &[e.id.clone()])
}
dto! {pub struct ModificationSummary{pub engine_version:String,pub state_revision:String,pub top_length:usize,pub bottom_length:usize,pub known_modified:usize,pub unknown_eligible:usize,pub known_unmodified:usize,pub generation:Vec<GenerationRecord>}}
impl std::ops::Deref for TrackedDuplex {
    type Target = Duplex;
    fn deref(&self) -> &Duplex {
        &self.duplex
    }
}
fn bad(s: &str) -> Error {
    Error::Unsupported(format!("tracked duplex: {s}"))
}
struct Bound(usize);
impl Write for Bound {
    fn write(&mut self, b: &[u8]) -> std::io::Result<usize> {
        if b.len() > MAX_BYTES.saturating_sub(self.0) {
            return Err(std::io::Error::other("tracked state size"));
        }
        self.0 += b.len();
        Ok(b.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}
fn bounded(v: &impl Serialize) -> Result<(), Error> {
    serde_json::to_writer(Bound(0), v).map_err(|_| Error::TooLarge)
}
fn hash(v: &impl Serialize) -> Result<String, Error> {
    bounded(v)?;
    Ok(
        blake3::hash(&serde_json::to_vec(v).map_err(|_| Error::TooLarge)?)
            .to_hex()
            .to_string(),
    )
}
fn base(c: Chemical) -> u8 {
    match c {
        Chemical::N6A => b'A',
        Chemical::C5C => b'C',
    }
}
fn uniform(state: CallState) -> StrandModifications {
    StrandModifications {
        n6a_default: state,
        c5c_default: state,
        spans: vec![],
    }
}
impl StrandModifications {
    fn default_for(&self, c: Chemical) -> CallState {
        match c {
            Chemical::N6A => self.n6a_default,
            Chemical::C5C => self.c5c_default,
        }
    }
    fn validate_for(&self, seq: &Sequence) -> Result<(), Error> {
        if self.spans.len() > MAX_SPANS {
            return Err(Error::TooLarge);
        }
        let mut previous = None;
        for s in &self.spans {
            if s.start >= s.end
                || s.end > seq.len()
                || previous
                    .is_some_and(|(c, end)| s.chemical < c || (s.chemical == c && s.start < end))
            {
                return Err(bad(
                    "overlapping, unsorted or out-of-strand modification spans",
                ));
            }
            if s.state == CallState::Modified
                && !seq.as_str().as_bytes()[s.start..s.end]
                    .iter()
                    .all(|b| *b == base(s.chemical))
            {
                return Err(bad("modified span includes chemically incompatible base"));
            }
            previous = Some((s.chemical, s.end));
        }
        Ok(())
    }
    fn at(&self, c: Chemical, index: usize) -> CallState {
        let at = self
            .spans
            .partition_point(|s| (s.chemical, s.start) <= (c, index));
        if at > 0 {
            let s = &self.spans[at - 1];
            if s.chemical == c && index < s.end {
                return s.state;
            }
        }
        self.default_for(c)
    }
    fn intervals(&self, c: Chemical, start: usize, end: usize) -> Vec<(usize, usize, CallState)> {
        let mut out = vec![];
        let mut cursor = start;
        for s in self
            .spans
            .iter()
            .filter(|s| s.chemical == c && s.end > start && s.start < end)
        {
            let a = s.start.max(start);
            let b = s.end.min(end);
            if cursor < a {
                out.push((cursor, a, self.default_for(c)))
            }
            out.push((a, b, s.state));
            cursor = b;
        }
        if cursor < end {
            out.push((cursor, end, self.default_for(c)))
        }
        out
    }
}
fn append_span(
    out: &mut StrandModifications,
    c: Chemical,
    start: usize,
    end: usize,
    state: CallState,
    seq: &Sequence,
) -> Result<(), Error> {
    if start == end || state == CallState::Unknown {
        return Ok(());
    }
    if state == CallState::Modified {
        let bytes = seq.as_str().as_bytes();
        let mut at = start;
        while at < end {
            if bytes[at] != base(c) {
                at += 1;
                continue;
            }
            let begin = at;
            while at < end && bytes[at] == base(c) {
                at += 1
            }
            if out.spans.len() >= MAX_SPANS {
                return Err(Error::TooLarge);
            }
            out.spans.push(ModificationSpan {
                chemical: c,
                start: begin,
                end: at,
                state,
            });
        }
    } else {
        if out.spans.len() >= MAX_SPANS {
            return Err(Error::TooLarge);
        }
        out.spans.push(ModificationSpan {
            chemical: c,
            start,
            end,
            state,
        });
    }
    Ok(())
}
fn normalize(out: &mut StrandModifications, seq: &Sequence) -> Result<(), Error> {
    out.spans.sort_by_key(|s| (s.chemical, s.start));
    let mut merged: Vec<ModificationSpan> = vec![];
    for s in std::mem::take(&mut out.spans) {
        if let Some(last) = merged.last_mut() {
            if last.chemical == s.chemical && last.end == s.start && last.state == s.state {
                last.end = s.end;
                continue;
            }
        }
        merged.push(s);
    }
    out.spans = merged;
    for c in [Chemical::N6A, Chemical::C5C] {
        let spans = out
            .spans
            .iter()
            .filter(|s| s.chemical == c)
            .collect::<Vec<_>>();
        if spans.len() == 1 && spans[0].start == 0 && spans[0].end == seq.len() {
            let value = spans[0].state;
            match c {
                Chemical::N6A => out.n6a_default = value,
                Chemical::C5C => out.c5c_default = value,
            }
            out.spans.retain(|s| s.chemical != c);
        }
    }
    out.validate_for(seq)
}
struct Piece<'a> {
    source: &'a TrackedDuplex,
    strand: Strand,
    start: usize,
    len: usize,
    destination: Strand,
    offset: usize,
}

fn from_pieces(
    duplex: Duplex,
    pieces: &[Piece<'_>],
    synthesis: &[(Strand, usize, usize)],
    generation: Vec<GenerationRecord>,
) -> Result<TrackedDuplex, Error> {
    let mut out = TrackedDuplex {
        duplex,
        top_modifications: uniform(CallState::Unknown),
        bottom_modifications: uniform(CallState::Unknown),
        generation,
    };
    let mut covered = [vec![], vec![]];
    for p in pieces {
        let (seq, mods) = p.source.strand(p.strand);
        let target = match p.destination {
            Strand::Top => &out.duplex.top,
            Strand::Bottom => &out.duplex.bottom,
        };
        if p.start.checked_add(p.len).is_none_or(|end| end > seq.len())
            || p.offset
                .checked_add(p.len)
                .is_none_or(|end| end > target.len())
            || seq.as_str()[p.start..p.start + p.len] != target.as_str()[p.offset..p.offset + p.len]
        {
            return Err(bad(
                "source/product strand mapping differs from actual bases",
            ));
        }
        covered[usize::from(p.destination == Strand::Bottom)].push((p.offset, p.offset + p.len));
        let dest = match p.destination {
            Strand::Top => &mut out.top_modifications,
            Strand::Bottom => &mut out.bottom_modifications,
        };
        for c in [Chemical::N6A, Chemical::C5C] {
            for (a, b, state) in mods.intervals(c, p.start, p.start + p.len) {
                append_span(
                    dest,
                    c,
                    p.offset + a - p.start,
                    p.offset + b - p.start,
                    state,
                    target,
                )?;
            }
        }
    }
    for (s, a, b) in synthesis {
        let target = match s {
            Strand::Top => &out.duplex.top,
            Strand::Bottom => &out.duplex.bottom,
        };
        if a > b || *b > target.len() {
            return Err(bad("synthesized interval outside product strand"));
        }
        covered[usize::from(*s == Strand::Bottom)].push((*a, *b));
        let dest = match s {
            Strand::Top => &mut out.top_modifications,
            Strand::Bottom => &mut out.bottom_modifications,
        };
        for c in [Chemical::N6A, Chemical::C5C] {
            append_span(dest, c, *a, *b, CallState::Unmodified, target)?;
        }
    }
    for (s, mut intervals) in [Strand::Top, Strand::Bottom].into_iter().zip(covered) {
        intervals.sort_unstable();
        let mut end = 0;
        for (a, b) in intervals {
            if a != end || b < a {
                return Err(bad(
                    "overlapping or missing strand ancestry/synthesis intervals",
                ));
            }
            end = b;
        }
        if end != out.strand(s).0.len() {
            return Err(bad("incomplete strand ancestry"));
        }
    }
    normalize(&mut out.top_modifications, &out.duplex.top)?;
    normalize(&mut out.bottom_modifications, &out.duplex.bottom)?;
    out.validate()?;
    Ok(out)
}
fn provenance<'a>(
    sources: impl IntoIterator<Item = &'a TrackedDuplex>,
) -> Result<Vec<GenerationRecord>, Error> {
    let mut records = vec![];
    for source in sources {
        for record in &source.generation {
            if !records.contains(record) {
                if records.len() >= 1000 {
                    return Err(Error::TooLarge);
                }
                records.push(record.clone());
            }
        }
    }
    Ok(records)
}
fn duplex_preflight(duplex: &Duplex) -> Result<(), Error> {
    if duplex.nicks.len() > 10_000
        || duplex.top.len() > collagene_model::MAX_SEQUENCE_LEN
        || duplex.bottom.len() > collagene_model::MAX_SEQUENCE_LEN
    {
        return Err(Error::TooLarge);
    }
    bounded(duplex)?;
    duplex.validate()
}
impl TrackedDuplex {
    fn strand(&self, s: Strand) -> (&Sequence, &StrandModifications) {
        match s {
            Strand::Top => (&self.duplex.top, &self.top_modifications),
            Strand::Bottom => (&self.duplex.bottom, &self.bottom_modifications),
        }
    }
    /// Summarize actual positional states after validating this tracked operand.
    pub fn summary(&self) -> Result<ModificationSummary, Error> {
        self.validate()?;
        let mut counts = [0usize; 3];
        for s in [Strand::Top, Strand::Bottom] {
            let (seq, mods) = self.strand(s);
            for c in [Chemical::N6A, Chemical::C5C] {
                for (a, b, state) in mods.intervals(c, 0, seq.len()) {
                    let count = seq.as_str().as_bytes()[a..b]
                        .iter()
                        .filter(|b| **b == base(c))
                        .count();
                    counts[match state {
                        CallState::Unknown => 0,
                        CallState::Unmodified => 1,
                        CallState::Modified => 2,
                    }] += count;
                }
            }
        }
        Ok(ModificationSummary {
            engine_version: ENGINE_VERSION.into(),
            state_revision: hash(self)?,
            top_length: self.duplex.top.len(),
            bottom_length: self.duplex.bottom.len(),
            known_modified: counts[2],
            unknown_eligible: counts[0],
            known_unmodified: counts[1],
            generation: self.generation.clone(),
        })
    }
    pub fn checked_digest_preflight(&self, o: &CheckedDigestOptions) -> Result<(), Error> {
        o.validate()?;
        if self.duplex.top.len() > MAX_BASES || self.duplex.bottom.len() > MAX_BASES {
            return Err(Error::TooLarge);
        }
        self.validate()
    }
    fn index_at(&self, s: Strand, reference: usize) -> Option<usize> {
        match s {
            Strand::Top => (reference < self.duplex.top.len()).then_some(reference),
            Strand::Bottom => {
                let end = self.duplex.bottom_start + self.duplex.bottom.len() as i64;
                let index = end - 1 - reference as i64;
                if self.duplex.topology == Topology::Circular {
                    Some(index.rem_euclid(self.duplex.bottom.len() as i64) as usize)
                } else {
                    usize::try_from(index)
                        .ok()
                        .filter(|i| *i < self.duplex.bottom.len())
                }
            }
        }
    }
    pub fn validate(&self) -> Result<(), Error> {
        if self.duplex.top.len() > collagene_model::MAX_SEQUENCE_LEN
            || self.duplex.bottom.len() > collagene_model::MAX_SEQUENCE_LEN
            || self.top_modifications.spans.len() > MAX_SPANS
            || self.bottom_modifications.spans.len() > MAX_SPANS
            || self.generation.len() > 1000
        {
            return Err(Error::TooLarge);
        }
        bounded(self)?;
        self.duplex.validate()?;
        self.top_modifications.validate_for(&self.duplex.top)?;
        self.bottom_modifications
            .validate_for(&self.duplex.bottom)?;
        for g in &self.generation {
            if g.source_revision.len() != 64
                || !g
                    .source_revision
                    .bytes()
                    .all(|b| b.is_ascii_hexdigit() && !b.is_ascii_uppercase())
            {
                return Err(bad("invalid source-generation revision"));
            }
        }
        Ok(())
    }
    pub fn from_unknown_duplex(duplex: &Duplex) -> Result<Self, Error> {
        duplex_preflight(duplex)?;
        let out = Self {
            duplex: duplex.clone(),
            top_modifications: uniform(CallState::Unknown),
            bottom_modifications: uniform(CallState::Unknown),
            generation: vec![],
        };
        out.validate()?;
        Ok(out)
    }
    pub fn exact_pcr(sequence: &Sequence) -> Result<Self, Error> {
        Self::annealed_unmodified(&Duplex::exact_pcr(sequence)?)
    }
    /// Only for newly reconstructed standard-unmodified oligos, not imported state.
    pub fn annealed_unmodified(duplex: &Duplex) -> Result<Self, Error> {
        duplex_preflight(duplex)?;
        let out = Self {
            duplex: duplex.clone(),
            top_modifications: uniform(CallState::Unmodified),
            bottom_modifications: uniform(CallState::Unmodified),
            generation: vec![],
        };
        out.validate()?;
        Ok(out)
    }
    pub fn from_input(m: &Molecule, a: &GenerationAssumptions) -> Result<Self, Error> {
        let o = methylation::MethylationOptions::default();
        methylation::source_preflight(m, &o)?;
        let duplex = Duplex::from_molecule(m)?;
        let mut out = Self {
            duplex,
            top_modifications: uniform(CallState::Unmodified),
            bottom_modifications: uniform(CallState::Unmodified),
            generation: vec![GenerationRecord {
                source_revision: m.revision(),
                assumptions: a.clone(),
            }],
        };
        let n = m.sequence.len();
        for (pattern, c, offset, state) in [
            ("GATC", Chemical::N6A, 1, a.dam),
            ("CCAGG", Chemical::C5C, 1, a.dcm),
            ("CCTGG", Chemical::C5C, 1, a.dcm),
            ("CG", Chemical::C5C, 0, a.cpg),
        ] {
            let call = match state {
                MethylationState::Unmodified => CallState::Unmodified,
                MethylationState::FullyModified => CallState::Modified,
                _ => CallState::Unknown,
            };
            if call == CallState::Unmodified || pattern.len() > n {
                continue;
            }
            let limit = if m.topology == Topology::Circular {
                n
            } else {
                n - pattern.len() + 1
            };
            for start in 0..limit {
                if pattern
                    .bytes()
                    .enumerate()
                    .all(|(j, b)| m.sequence.as_str().as_bytes()[(start + j) % n] == b)
                {
                    let top = (start + offset) % n;
                    let bottom = n - 1 - (start + pattern.len() - 1 - offset) % n;
                    for (mods, index) in [
                        (&mut out.top_modifications, top),
                        (&mut out.bottom_modifications, bottom),
                    ] {
                        if mods.spans.len() >= MAX_SPANS {
                            return Err(Error::TooLarge);
                        }
                        mods.spans.push(ModificationSpan {
                            chemical: c,
                            start: index,
                            end: index + 1,
                            state: call,
                        });
                    }
                }
            }
        } // Unknown source targets override the unmodified default, so retain these spans.
        for mods in [&mut out.top_modifications, &mut out.bottom_modifications] {
            mods.spans.sort_by_key(|s| (s.chemical, s.start));
            let mut merged: Vec<ModificationSpan> = vec![];
            for s in std::mem::take(&mut mods.spans) {
                if let Some(last) = merged.last_mut() {
                    if last.chemical == s.chemical && last.end == s.start && last.state == s.state {
                        last.end = s.end;
                        continue;
                    }
                }
                merged.push(s);
            }
            mods.spans = merged;
        }
        out.validate()?;
        Ok(out)
    }
    pub fn state_at(
        &self,
        s: Strand,
        c: Chemical,
        index: usize,
    ) -> Result<Option<CallState>, Error> {
        let (seq, mods) = self.strand(s);
        if index >= seq.len() {
            return Err(bad("modification lookup outside strand"));
        }
        Ok((seq.as_str().as_bytes()[index] == base(c)).then(|| mods.at(c, index)))
    }
    pub fn reverse(&self) -> Result<Self, Error> {
        self.validate()?;
        let mut out = self.clone();
        out.duplex = self.duplex.reverse()?;
        std::mem::swap(&mut out.top_modifications, &mut out.bottom_modifications);
        out.validate()?;
        Ok(out)
    }
    pub fn treat(&self, treatment: EndTreatment) -> Result<Self, Error> {
        self.validate()?;
        let mut out = self.clone();
        out.duplex = self.duplex.treat(treatment)?;
        out.validate()?;
        Ok(out)
    }
    pub fn repair(&self, method: EndRepairMethod) -> Result<Self, Error> {
        self.validate()?;
        let report = chemistry::repair_ends(&self.duplex, method)?;
        let mut pieces = vec![];
        let mut synth = vec![];
        for m in &report.mappings {
            if m.synthesized {
                synth.push((m.product_strand, m.product_start, m.product_end));
            } else {
                if m.reverse {
                    return Err(bad("unsupported non-synthesis repair reversal"));
                }
                pieces.push(Piece {
                    source: self,
                    strand: m.source_strand,
                    start: m.source_start,
                    len: m.source_end - m.source_start,
                    destination: m.product_strand,
                    offset: m.product_start,
                });
            }
        }
        from_pieces(report.duplex, &pieces, &synth, self.generation.clone())
    }
    /// Explicitly denature/reanneal A.top with B.bottom and extend both paired
    /// 3′ ends. Only the selected strands retain their modification coverage;
    /// complementary synthesis uses the standard-unmodified assumption.
    ///
    /// The geometry is freshly reconstructed from these actual operands, never
    /// accepted from a serialized report or an intended reference sequence.
    /// Terminal-pair search completeness and uniqueness remain explicit in the
    /// returned report and must be assessed before claiming a determined product.
    pub fn overlap_extend(
        &self,
        other: &Self,
        options: &crate::cloning::overlap::extension::Options,
    ) -> Result<(Self, crate::cloning::overlap::extension::Geometry), Error> {
        self.validate()?;
        other.validate()?;
        let geometry =
            crate::cloning::overlap::extension::geometry(&self.duplex, &other.duplex, options)?;
        let pieces = [
            Piece {
                source: self,
                strand: Strand::Top,
                start: 0,
                len: self.duplex.top.len(),
                destination: Strand::Top,
                offset: 0,
            },
            Piece {
                source: other,
                strand: Strand::Bottom,
                start: 0,
                len: other.duplex.bottom.len(),
                destination: Strand::Bottom,
                offset: 0,
            },
        ];
        let annealed = from_pieces(
            geometry.annealed.clone(),
            &pieces,
            &[],
            provenance([self, other])?,
        )?;
        let extended = annealed.repair(EndRepairMethod::FillFivePrimeOverhangs)?;
        if extended.duplex != geometry.fill.duplex {
            return Err(bad(
                "tracked reciprocal extension differs from freshly computed native geometry",
            ));
        }
        Ok((extended, geometry))
    }
    pub fn tail(&self, base: TailingBase) -> Result<Self, Error> {
        self.validate()?;
        let report = tailing::tail_blunt_duplex(&self.duplex, base)?;
        let pieces = [
            Piece {
                source: self,
                strand: Strand::Top,
                start: 0,
                len: self.duplex.top.len(),
                destination: Strand::Top,
                offset: 0,
            },
            Piece {
                source: self,
                strand: Strand::Bottom,
                start: 0,
                len: self.duplex.bottom.len(),
                destination: Strand::Bottom,
                offset: 0,
            },
        ];
        let synth = report
            .additions
            .iter()
            .map(|a| (a.strand, a.product_start, a.product_end))
            .collect::<Vec<_>>();
        from_pieces(report.duplex, &pieces, &synth, self.generation.clone())
    }
}
fn slice_pieces<'a>(
    source: &'a TrackedDuplex,
    s: Strand,
    start: i64,
    len: usize,
    destination: Strand,
) -> Result<Vec<Piece<'a>>, Error> {
    let n = source.strand(s).0.len();
    if len > n {
        return Err(bad("fragment traverses a source strand more than once"));
    }
    let start = if source.duplex.topology == Topology::Circular {
        start.rem_euclid(n as i64) as usize
    } else {
        usize::try_from(start).map_err(|_| bad("fragment starts outside actual source strand"))?
    };
    if start >= n || source.duplex.topology == Topology::Linear && start + len > n {
        return Err(bad("fragment exceeds actual source strand"));
    }
    let first = len.min(n - start);
    let mut out = vec![Piece {
        source,
        strand: s,
        start,
        len: first,
        destination,
        offset: 0,
    }];
    if first < len {
        out.push(Piece {
            source,
            strand: s,
            start: 0,
            len: len - first,
            destination,
            offset: first,
        });
    }
    Ok(out)
}
fn map_fragment_validated(
    source: &TrackedDuplex,
    f: &chemistry::DigestFragment,
) -> Result<TrackedDigestFragment, Error> {
    f.duplex.validate()?;
    if f.top_start >= source.duplex.top.len()
        || f.bottom_start.checked_sub(f.top_start as i64) != Some(f.duplex.bottom_start)
    {
        return Err(bad(
            "fragment source boundaries contradict actual strand geometry",
        ));
    }
    if f.duplex.topology == Topology::Circular && f.duplex != source.duplex {
        return Err(bad("a circular digest survivor must be the intact source"));
    }
    let mut pieces = slice_pieces(
        source,
        Strand::Top,
        f.top_start as i64,
        f.duplex.top.len(),
        Strand::Top,
    )?;
    let bottom_end = source.duplex.bottom_start + source.duplex.bottom.len() as i64;
    let bottom_index = bottom_end
        .checked_sub(f.bottom_start)
        .and_then(|p| p.checked_sub(f.duplex.bottom.len() as i64))
        .ok_or_else(|| bad("extreme source bottom boundary"))?;
    pieces.extend(slice_pieces(
        source,
        Strand::Bottom,
        bottom_index,
        f.duplex.bottom.len(),
        Strand::Bottom,
    )?);
    let duplex = from_pieces(f.duplex.clone(), &pieces, &[], source.generation.clone())?;
    Ok(TrackedDigestFragment {
        duplex,
        top_start: f.top_start,
        bottom_start: f.bottom_start,
    })
}
pub fn map_digest_fragment(
    source: &TrackedDuplex,
    f: &chemistry::DigestFragment,
) -> Result<TrackedDigestFragment, Error> {
    source.validate()?;
    map_fragment_validated(source, f)
}
pub fn digest(
    source: &TrackedDuplex,
    enzymes: &[Enzyme],
) -> Result<Vec<TrackedDigestFragment>, Error> {
    source.validate()?;
    let fragments = chemistry::digest(&source.duplex, enzymes)?;
    let mut out = vec![];
    let mut bytes = 0usize;
    for f in &fragments {
        let tracked = map_fragment_validated(source, f)?;
        let size = serde_json::to_vec(&tracked)
            .map_err(|_| Error::TooLarge)?
            .len();
        bytes = bytes
            .checked_add(size)
            .filter(|n| *n <= MAX_BYTES)
            .ok_or(Error::TooLarge)?;
        out.push(tracked);
    }
    Ok(out)
}
pub fn ligate(
    sources: &[TrackedDuplex],
    circular: bool,
    allow_nicked: bool,
) -> Result<TrackedDuplex, Error> {
    if sources.is_empty() || sources.len() > 1000 {
        return Err(bad("ligation needs1..1000 sources"));
    }
    for s in sources {
        s.validate()?;
    }
    let ds = sources.iter().map(|s| s.duplex.clone()).collect::<Vec<_>>();
    let product = chemistry::ligate(&ds, circular, allow_nicked)?;
    let mut pieces = vec![];
    let mut top_offset = 0;
    let mut bottom_offset = 0;
    for source in sources {
        pieces.push(Piece {
            source,
            strand: Strand::Top,
            start: 0,
            len: source.duplex.top.len(),
            destination: Strand::Top,
            offset: top_offset,
        });
        top_offset += source.duplex.top.len();
    }
    for source in sources.iter().rev() {
        pieces.push(Piece {
            source,
            strand: Strand::Bottom,
            start: 0,
            len: source.duplex.bottom.len(),
            destination: Strand::Bottom,
            offset: bottom_offset,
        });
        bottom_offset += source.duplex.bottom.len();
    }
    if !circular {
        return from_pieces(product, &pieces, &[], provenance(sources)?);
    } // First form the actual concatenated linear bottom strand, then rotate its ledger
    // using exactly the alignment performed by chemistry::ligate.
    let linear = chemistry::ligate(&ds, false, allow_nicked)?;
    let tracked = from_pieces(linear, &pieces, &[], provenance(sources)?)?;
    let n = tracked.duplex.top.len();
    let shift = tracked.duplex.bottom_start.rem_euclid(n as i64) as usize;
    let mut aligned = vec![Piece {
        source: &tracked,
        strand: Strand::Top,
        start: 0,
        len: n,
        destination: Strand::Top,
        offset: 0,
    }];
    if shift == 0 {
        aligned.push(Piece {
            source: &tracked,
            strand: Strand::Bottom,
            start: 0,
            len: n,
            destination: Strand::Bottom,
            offset: 0,
        });
    } else {
        aligned.push(Piece {
            source: &tracked,
            strand: Strand::Bottom,
            start: shift,
            len: n - shift,
            destination: Strand::Bottom,
            offset: 0,
        });
        aligned.push(Piece {
            source: &tracked,
            strand: Strand::Bottom,
            start: 0,
            len: shift,
            destination: Strand::Bottom,
            offset: n - shift,
        });
    }
    from_pieces(product, &aligned, &[], tracked.generation.clone())
}

/// Reconstruct private homology-assembly state only from freshly supplied tracked
/// operands. Documentary mappings cannot be used as a material-import channel.
pub fn homology_assemble(
    sources: &[&TrackedDuplex],
    options: &crate::cloning::overlap::assembly::Options,
) -> Result<
    (
        Option<TrackedDuplex>,
        crate::cloning::overlap::assembly::Report,
    ),
    Error,
> {
    options.validate_shape(sources.len())?;
    crate::cloning::overlap::assembly::preflight_borrowed(
        &sources.iter().map(|s| &s.duplex).collect::<Vec<_>>(),
        options,
    )?;
    // Guard aggregate source material before the geometry clone. Source metadata
    // never authorizes state; callers must resolve actual tracked operands.
    let mut total = 0usize;
    for source in sources {
        total = total
            .checked_add(source.top.len())
            .and_then(|n| n.checked_add(source.bottom.len()))
            .filter(|n| *n <= crate::cloning::overlap::assembly::MAX_SOURCE_BASES)
            .ok_or(Error::TooLarge)?;
        source.validate()?;
    }
    let duplexes: Vec<_> = sources.iter().map(|s| s.duplex.clone()).collect();
    let mut report = crate::cloning::overlap::assembly::geometry(&duplexes, options)?;
    // Current Duplex::Nick records a boundary, not its 3′/5′ chemical groups.
    // The diagnostic report preserves those groups, but reducing blocked or
    // unknown junction chemistry to that private representation would be lossy.
    // Keep reference geometry available without granting reusable material.
    if report.product.as_ref().is_some_and(|p| !p.nicks.is_empty()) {
        report.incomplete = true;
        report.limits.push("Reusable material withheld: the current private nick representation cannot preserve the unsealed junction's blocked/unknown chemical groups; full nick chemistry remains unsupported.".into());
        report.validate_documentary()?;
        return Ok((None, report));
    }
    let Some(product) = &report.product else {
        return Ok((None, report));
    };
    let mut pieces = vec![];
    let mut synthesis = vec![];
    for sourced in &report.mappings {
        let m = &sourced.mapping;
        if m.synthesized {
            synthesis.push((m.product_strand, m.product_start, m.product_end));
        } else {
            pieces.push(Piece {
                source: sources[sourced.source],
                strand: m.source_strand,
                start: m.source_start,
                len: m.source_end - m.source_start,
                destination: m.product_strand,
                offset: m.product_start,
            });
        }
    }
    let tracked = from_pieces(
        product.clone(),
        &pieces,
        &synthesis,
        provenance(sources.iter().copied())?,
    )?;
    Ok((Some(tracked), report))
}
// Source-bound bridge: no serialized report is accepted as tracked material.
fn tracked_resection(
    source: &TrackedDuplex,
    report: &crate::cloning::overlap::resection::Report,
) -> Result<TrackedDuplex, Error> {
    let pieces = report
        .mappings
        .iter()
        .map(|m| Piece {
            source,
            strand: m.source_strand,
            start: m.source_start,
            len: m.source_end - m.source_start,
            destination: m.product_strand,
            offset: m.product_start,
        })
        .collect::<Vec<_>>();
    from_pieces(
        report.product.clone(),
        &pieces,
        &[],
        source.generation.clone(),
    )
}
impl TrackedDuplex {
    /// Fresh hydrolytic shortening of this actual tracked operand. No assumed
    /// phosphorylation, synthesis, blocked-end removal or joining is performed.
    pub fn resect_three_prime(
        &self,
        cuts: &crate::cloning::overlap::resection::Cuts,
    ) -> Result<(Self, crate::cloning::overlap::resection::Report), Error> {
        crate::cloning::overlap::resection::preflight(&self.duplex, cuts)?;
        self.validate()?;
        let report = crate::cloning::overlap::resection::geometry(&self.duplex, cuts)?;
        Ok((tracked_resection(self, &report)?, report))
    }
}
/// Retain each source's separate processed covalent strands and modifications.
/// The paired/gapped complex is documentary; this deliberately returns no
/// concatenated assembled TrackedDuplex or invented host-repair authority.
pub fn homology_anneal(
    sources: &[&TrackedDuplex],
    options: &crate::cloning::overlap::annealing::Options,
) -> Result<
    (
        Vec<TrackedDuplex>,
        crate::cloning::overlap::annealing::Report,
    ),
    Error,
> {
    options.validate_shape(sources.len())?;
    let duplexes = sources.iter().map(|s| &s.duplex).collect::<Vec<_>>();
    crate::cloning::overlap::annealing::preflight_borrowed(&duplexes, options)?;
    // Bound combined positional/generation metadata before copying any operand.
    crate::cloning::overlap::resection::revision(&sources)?;
    for source in sources {
        source.validate()?;
    }
    let report = crate::cloning::overlap::annealing::geometry_borrowed(&duplexes, options)?;
    let operands = sources
        .iter()
        .zip(&report.preparations)
        .map(|(source, preparation)| tracked_resection(source, preparation))
        .collect::<Result<Vec<_>, _>>()?;
    Ok((operands, report))
}

pub fn profile_revision() -> String {
    hash(&(PROFILE_VERSION,methylation::profile_revision(),"Exact canonical catalog identity/geometry; actual unmodified recognition permits arbitrary catalog cuts. NotI/SmaI full internal CpG and XbaI overlapping Dam rules; inherited unknown/modified positions require exact matching positional rules. Shared bounded recognition/context ledger; no formulation substitution.")).expect("fixed bounded profile")
}
impl CheckedDigestOptions {
    pub fn validate(&self) -> Result<(), Error> {
        self.selection.validate()?;
        if self.selection.enzymes.len() > 12
            || self.max_comparisons == 0
            || self.max_comparisons > 100_000_000
            || self.max_sites == 0
            || self.max_sites > 4096
        {
            return Err(bad("invalid checked enzyme/site/comparison bounds"));
        }
        Ok(())
    }
}

fn evidence_at(
    source: &TrackedDuplex,
    s: Strand,
    c: Chemical,
    reference: usize,
) -> Option<PositionEvidence> {
    let index = source.index_at(s, reference)?;
    source
        .state_at(s, c, index)
        .ok()
        .flatten()
        .map(|state| PositionEvidence {
            strand: s,
            index,
            chemical: c,
            state,
        })
}
fn site_outcome(
    source: &TrackedDuplex,
    site: &methylation::MethylationSite,
    profile: &methylation::EnzymeProfile,
) -> (Outcome, Vec<PositionEvidence>, String) {
    let n = source.duplex.top.len();
    let mut evidence = vec![];
    for j in 0..site.recognition_length {
        let at = (site.recognition_start + j) % n;
        for s in [Strand::Top, Strand::Bottom] {
            for c in [Chemical::N6A, Chemical::C5C] {
                if let Some(e) = evidence_at(source, s, c, at) {
                    evidence.push(e);
                }
            }
        }
    }
    let paired = source.duplex.topology == Topology::Circular
        || (site.recognition_start as i64 >= source.duplex.bottom_start.max(0)
            && (site.recognition_start + site.recognition_length) as i64
                <= (source.duplex.bottom_start + source.duplex.bottom.len() as i64).min(n as i64));
    if !paired {
        return (
            Outcome::Unknown,
            evidence,
            "Recognition sequence is not established as a complete paired-strand substrate.".into(),
        );
    }
    let inside = |position: usize| {
        (0..site.recognition_length).any(|j| (site.recognition_start + j) % n == position)
    };
    // Full internal CpG modification on both strands is the published blocking
    // profile. Partial/hemi subsets do not become a fully modified substrate.
    if matches!(site.enzyme.as_str(), "NotI" | "SmaI") {
        let cpg = site.effects[2]
            .overlaps
            .iter()
            .filter(|m| m.positions.iter().all(|p| inside(*p)))
            .collect::<Vec<_>>();
        if !cpg.is_empty()
            && cpg.iter().all(|m| {
                m.modified_bases.iter().all(|b| {
                    evidence_at(
                        source,
                        if b.reverse {
                            Strand::Bottom
                        } else {
                            Strand::Top
                        },
                        Chemical::C5C,
                        b.position,
                    )
                    .is_some_and(|e| e.state == CallState::Modified)
                })
            })
        {
            return(Outcome::Blocked,evidence,"All internal CpG targets on both actual strands are 5-methylcytosine; the NotI/SmaI full-CpG vendor profile blocks cleavage.".into());
        }
    }
    // For an intact overlapping Dam motif, vendor XbaI context guidance applies
    // to its actual two methyladenines. Broken contexts are NOT used to erase state.
    if site.enzyme == "XbaI" {
        for motif in &site.effects[0].overlaps {
            let pair = motif
                .modified_bases
                .iter()
                .filter_map(|b| {
                    evidence_at(
                        source,
                        if b.reverse {
                            Strand::Bottom
                        } else {
                            Strand::Top
                        },
                        Chemical::N6A,
                        b.position,
                    )
                })
                .collect::<Vec<_>>();
            if pair.len() == 2 && pair.iter().all(|e| e.state == CallState::Modified) {
                for e in pair {
                    if !evidence.contains(&e) {
                        evidence.push(e);
                    }
                }
                return(Outcome::Blocked,evidence,"Actual full Dam methyladenines at a published overlapping XbaI context block cleavage.".into());
            }
        }
    }
    let unresolved = evidence.iter().any(|e| {
        if e.state == CallState::Unmodified {
            return false;
        }
        let reference = match e.strand {
            Strand::Top => e.index,
            Strand::Bottom => {
                let r = source.duplex.bottom_start + source.duplex.bottom.len() as i64
                    - 1
                    - e.index as i64;
                if source.duplex.topology == Topology::Circular {
                    r.rem_euclid(n as i64) as usize
                } else {
                    r as usize
                }
            }
        };
        !site.effects.iter().enumerate().any(|(i, effect)| {
            let chemical = if i == 0 { Chemical::N6A } else { Chemical::C5C };
            chemical == e.chemical
                && profile.sensitivities[i] == Sensitivity::NotSensitive
                && effect.overlaps.iter().any(|m| {
                    m.certain
                        && m.modified_bases.iter().any(|b| {
                            b.position == reference && b.reverse == (e.strand == Strand::Bottom)
                        })
                })
        })
    });
    if unresolved {
        (Outcome::Unknown,evidence,"Inherited modified/unknown chemical positions lack a complete applicable positional rule. Missing or broken generation motifs do not establish unmodified bases; hemi/partial and some-combination effects remain unresolved.".into())
    } else {
        (Outcome::PermittedUnderProfile,evidence,"All relevant eligible chemical positions are known unmodified or covered by a published insensitive methylase context; reaction conditions remain assumptions.".into())
    }
}
pub fn checked_digest(
    source: &TrackedDuplex,
    o: &CheckedDigestOptions,
) -> Result<CheckedDigestReport, Error> {
    let mut budget = CleavageBudget::new(o.max_comparisons)?;
    checked_digest_with_budget(source, o, &mut budget)
}
/// Only fresh tracked operands enter this kernel. Catalog membership is checked
/// before scanning, and no identity inherits a rule merely by sharing geometry.
pub fn checked_digest_with_budget(
    source: &TrackedDuplex,
    o: &CheckedDigestOptions,
    budget: &mut CleavageBudget,
) -> Result<CheckedDigestReport, Error> {
    source.checked_digest_preflight(o)?;
    let catalog = crate::cloning::restriction::catalog::builtin_catalog()?;
    let enzymes = o.selection.resolved_enzymes(&catalog)?;
    let sequence = Molecule::dna(
        "checked source",
        source.duplex.top.as_str(),
        source.duplex.topology,
    )?;
    let initial = budget.used();
    let allowed = o.max_comparisons.min(budget.remaining());
    let mut r=CheckedDigestReport{engine_version:ENGINE_VERSION.into(),profile_version:PROFILE_VERSION.into(),profile_revision:profile_revision(),source_revision:hash(source)?,options:o.clone(),complete:true,products_determined:true,comparisons:0,sites:vec![],fragments:vec![],warnings:vec![
        "Actual positional N6-methyladenine/5-methylcytosine under declared generation and standard-unmodified synthesis assumptions; not measured modifications or experimental cleavage.".into(),
        "Exact catalog identity/geometry is retained. Unsupported identities do not inherit another enzyme/formulation's methylation rule; genuinely unmodified eligible positions need no sensitivity substitution.".into(),
        "Unknown/impaired sites withhold all determined fragments; blocked sites survive. Conditions, other modifications, kinetics and star activity remain outside this model.".into()],limits:vec![]};
    if allowed == 0 {
        r.complete = false;
        r.products_determined = false;
        r.limits.push("Shared cleavage comparison budget exhausted before scan; absent sites/products are not established.".into());
        r.validate()?;
        return Ok(r);
    }
    // Recognition scan supports all exact catalog geometries, including nonstarter
    // enzymes. It does not perform assumed-complete cleavage.
    let scan = crate::cloning::restriction::inspection::scan(
        "checked source",
        &sequence,
        &enzymes,
        &crate::cloning::restriction::inspection::InspectionOptions {
            max_comparisons: allowed as u64,
            max_sites: o.max_sites,
            max_fragments: o.max_sites + 1,
            max_report_bytes: 8 * 1024 * 1024,
        },
    )?;
    budget.charge(scan.comparisons as usize)?;
    if !scan.scan_complete {
        r.complete = false;
        r.products_determined = false;
        r.limits.push(
            "Bounded catalog recognition scan incomplete; missing sites are not excluded.".into(),
        );
    }
    // Existing positional rules are applicable only to the exact canonical
    // starter identity AND matching recognition/cut geometry. Context work is
    // charged to the same ledger, including unsuccessful scans.
    let starters = analysis::starter_enzymes();
    let supported = enzymes
        .iter()
        .filter(|e| {
            starters.iter().any(|p| {
                p.name == e.name
                    && p.recognition == e.recognition
                    && p.top_cut == e.top_cut
                    && p.bottom_cut == e.bottom_cut
            })
        })
        .map(|e| e.name.clone())
        .collect::<Vec<_>>();
    let published = methylation::profiles();
    for selected in &o.selection.enzymes {
        let note = if supported.contains(&selected.enzyme.name) {
            let p = published
                .iter()
                .find(|p| p.enzyme == selected.enzyme.name)
                .ok_or_else(|| bad("missing positional identity profile"))?;
            format!(
                "Exact positional model for catalog {}: {}. {}",
                selected.id,
                p.supplier_identity,
                p.notes.join(" ")
            )
        } else {
            format!(
                "Catalog {} has no exact applicable positional methylation profile here; only actually known unmodified relevant chemical positions permit cleavage. No same-cut/formulation substitution.",
                selected.id
            )
        };
        r.warnings.push(note);
    }
    // A proven absent site needs no methylase-context search. Likewise actual
    // known-unmodified recognition residues cannot require a sensitivity rule.
    // Reserve context work only for supported identities with an observed site
    // carrying a relevant modified/unknown actual position.
    let contextual_names = supported
        .iter()
        .filter(|name| {
            scan.sites.iter().any(|site| {
                if site.enzyme.as_str() != name.as_str() {
                    return false;
                }
                (0..site.recognition_length).any(|j| {
                    let at = (site.recognition_start + j) % source.duplex.top.len();
                    [Strand::Top, Strand::Bottom].into_iter().any(|strand| {
                        [Chemical::N6A, Chemical::C5C].into_iter().any(|chemical| {
                            evidence_at(source, strand, chemical, at)
                                .is_some_and(|e| e.state != CallState::Unmodified)
                        })
                    })
                })
            })
        })
        .cloned()
        .collect::<Vec<_>>();
    let remaining = allowed.saturating_sub(budget.used() - initial);
    let contextual = if !contextual_names.is_empty() && remaining > 0 {
        let c = methylation::assess(
            &sequence,
            &methylation::MethylationOptions {
                enzymes: contextual_names.clone(),
                dam: MethylationState::Unmodified,
                dcm: MethylationState::Unmodified,
                cpg: MethylationState::Unmodified,
                max_comparisons: remaining,
                max_sites: o.max_sites,
            },
        )?;
        budget.charge(c.comparisons)?;
        if !c.complete {
            r.complete = false;
            r.products_determined = false;
            r.limits.push(
                "Bounded positional context scan incomplete; no complete cleavage set established."
                    .into(),
            );
        }
        Some(c)
    } else {
        None
    };
    if remaining == 0 && !contextual_names.is_empty() {
        r.complete = false;
        r.products_determined = false;
        r.limits.push(
            "Shared cleavage comparison budget exhausted before positional context scan.".into(),
        );
    }
    let mut cuts = vec![];
    for site in &scan.sites {
        let mut evidence = vec![];
        for j in 0..site.recognition_length {
            let at = (site.recognition_start + j) % source.duplex.top.len();
            for strand in [Strand::Top, Strand::Bottom] {
                for chemical in [Chemical::N6A, Chemical::C5C] {
                    if let Some(e) = evidence_at(source, strand, chemical, at) {
                        evidence.push(e)
                    }
                }
            }
        }
        let n = source.duplex.top.len();
        let paired = source.duplex.topology == Topology::Circular
            || (site.recognition_start as i64 >= source.duplex.bottom_start.max(0)
                && (site.recognition_start + site.recognition_length) as i64
                    <= (source.duplex.bottom_start + source.duplex.bottom.len() as i64)
                        .min(n as i64));
        let (mut outcome, reason) = if !site.certain || !paired {
            (
                Outcome::Unknown,
                "Recognition is ambiguous or not an established complete paired-strand substrate."
                    .into(),
            )
        } else if let Some((c, known)) = contextual.as_ref().and_then(|c| {
            c.sites
                .iter()
                .find(|x| {
                    x.enzyme == site.enzyme
                        && x.recognition_start == site.recognition_start
                        && x.reverse == site.reverse
                })
                .map(|known| (c, known))
        }) {
            let profile = c
                .profiles
                .iter()
                .find(|p| p.enzyme == known.enzyme)
                .ok_or_else(|| bad("missing exact positional profile"))?;
            let (outcome, actual, reason) = site_outcome(source, known, profile);
            evidence = actual;
            (outcome, reason)
        } else if evidence.iter().all(|e| e.state == CallState::Unmodified) {
            (Outcome::PermittedUnderProfile,"Every relevant eligible chemical position is actually known unmodified; exact catalog cut geometry is used without borrowing a methylation profile.".into())
        } else {
            (Outcome::Unknown,"Modified/unknown actual positions lack an exact applicable positional rule; another formulation or same-cut enzyme is not substituted.".into())
        };
        if outcome == Outcome::PermittedUnderProfile
            && (site.top_cut.is_none() || site.bottom_cut.is_none())
        {
            outcome = Outcome::Unknown
        }
        if matches!(outcome, Outcome::Unknown | Outcome::Impaired) {
            r.products_determined = false
        }
        if outcome == Outcome::PermittedUnderProfile {
            let e = &enzymes[site.enzyme_index];
            cuts.push((
                analysis::CutSite {
                    enzyme: site.enzyme.clone(),
                    recognition_start: site.recognition_start,
                    reverse: site.reverse,
                    top_cut: site.top_cut.unwrap(),
                    bottom_cut: site.bottom_cut.unwrap(),
                    certain: true,
                },
                site.recognition_length,
                i64::from(e.bottom_cut) - i64::from(e.top_cut),
            ));
        }
        r.sites.push(CheckedSite {
            enzyme: site.enzyme.clone(),
            recognition_start: site.recognition_start,
            reverse: site.reverse,
            top_cut: site.top_cut,
            bottom_cut: site.bottom_cut,
            outcome,
            evidence,
            reason,
        });
    }
    r.comparisons = budget.used() - initial;
    if r.products_determined {
        match chemistry::digest_sites(&source.duplex, &cuts) {
            Ok(fragments) => {
                let mut bytes = serde_json::to_vec(&r).map_err(|_| Error::TooLarge)?.len();
                for f in &fragments {
                    let tracked = map_fragment_validated(source, f)?;
                    let size = serde_json::to_vec(&tracked)
                        .map_err(|_| Error::TooLarge)?
                        .len()
                        + 1;
                    if size > MAX_BYTES.saturating_sub(bytes) {
                        r.complete = false;
                        r.products_determined = false;
                        r.fragments.clear();
                        r.limits.push("Tracked fragment evidence exceeds bounded output; no complete product set returned.".into());
                        break;
                    }
                    bytes += size;
                    r.fragments.push(tracked);
                }
            }
            Err(Error::Unsupported(reason)) => {
                r.complete = false;
                r.products_determined = false;
                r.limits
                    .push(format!("Unsupported physical cut geometry: {reason}"));
            }
            Err(e) => return Err(e),
        }
    }
    r.validate()?;
    Ok(r)
}
impl CheckedDigestReport {
    /// Bounded internal shape checks; source authority requires fresh replay.
    pub fn validate(&self) -> Result<(), Error> {
        if self.engine_version.is_empty()
            || self.engine_version.len() > 128
            || self.profile_version.is_empty()
            || self.profile_version.len() > 128
            || self.profile_revision.len() != 64
            || self.source_revision.len() != 64
            || ![&self.profile_revision, &self.source_revision]
                .iter()
                .all(|s| {
                    s.bytes()
                        .all(|b| b.is_ascii_hexdigit() && !b.is_ascii_uppercase())
                })
            || self.options.selection.validate().is_err()
            || self.options.selection.enzymes.len() > 12
            || self.options.max_sites == 0
            || self.options.max_sites > 4096
            || self.options.max_comparisons == 0
            || self.options.max_comparisons > 100_000_000
            || self.sites.len() > self.options.max_sites
            || self.comparisons > self.options.max_comparisons
            || self.fragments.len() > self.options.max_sites + 1
            || self.limits.len() > 8
            || self.warnings.len() > 16
            || self
                .limits
                .iter()
                .chain(&self.warnings)
                .any(|s| s.len() > 4096)
            || self.complete != self.limits.is_empty()
            || (!self.complete || !self.products_determined) && !self.fragments.is_empty()
            || self.products_determined && self.fragments.is_empty()
        {
            return Err(bad("invalid bounded checked report shape"));
        }
        for site in &self.sites {
            if site.enzyme.len() > 64
                || !self
                    .options
                    .selection
                    .enzymes
                    .iter()
                    .any(|e| e.enzyme.name == site.enzyme)
                || site.evidence.len() > 256
                || site.reason.len() > 4096
                || self.products_determined
                    && matches!(site.outcome, Outcome::Unknown | Outcome::Impaired)
            {
                return Err(bad("invalid checked site evidence shape"));
            }
        }
        bounded(self)?;
        for fragment in &self.fragments {
            fragment.duplex.validate()?;
        }
        Ok(())
    }

    pub fn evidence(&self) -> CheckedCleavageEvidence {
        CheckedCleavageEvidence {
            engine_version: self.engine_version.clone(),
            profile_version: self.profile_version.clone(),
            profile_revision: self.profile_revision.clone(),
            source_revision: self.source_revision.clone(),
            options: self.options.clone(),
            complete: self.complete,
            products_determined: self.products_determined,
            comparisons: self.comparisons,
            sites: self.sites.clone(),
            warnings: self.warnings.clone(),
            limits: self.limits.clone(),
        }
    }
    pub fn validate_bindings_for(
        &self,
        source: &TrackedDuplex,
        o: &CheckedDigestOptions,
    ) -> Result<(), Error> {
        source.checked_digest_preflight(o)?;
        self.validate()?;
        if self.engine_version != ENGINE_VERSION
            || self.profile_version != PROFILE_VERSION
            || self.profile_revision != profile_revision()
            || self.source_revision != hash(source)?
            || self.options != *o
            || self.sites.len() > o.max_sites
            || self.comparisons > o.max_comparisons
            || self.fragments.len() > o.max_sites + 1
            || self.limits.len() > 8
            || self.warnings.len() > 16
            || (!self.complete || !self.products_determined) && !self.fragments.is_empty()
            || self.products_determined && (!self.complete || self.fragments.is_empty())
            || self.products_determined
                && self
                    .sites
                    .iter()
                    .any(|s| matches!(s.outcome, Outcome::Unknown | Outcome::Impaired))
        {
            return Err(bad(
                "checked report source/options/profile or shape mismatch",
            ));
        }
        for f in &self.fragments {
            f.duplex.validate()?;
        }
        Ok(())
    }
}
pub fn replay_checked(
    source: &TrackedDuplex,
    o: &CheckedDigestOptions,
    report: &CheckedDigestReport,
) -> Result<(), Error> {
    report.validate_bindings_for(source, o)?;
    if checked_digest(source, o)? != *report {
        return Err(bad("checked report differs from fresh positional replay"));
    }
    Ok(())
}

impl CheckedCleavageEvidence {
    pub fn validate(&self) -> Result<(), Error> {
        self.options.validate()?;
        if self.engine_version.is_empty()
            || self.engine_version.len() > 128
            || self.profile_version.is_empty()
            || self.profile_version.len() > 128
            || [&self.profile_revision, &self.source_revision]
                .iter()
                .any(|s| {
                    s.len() != 64
                        || !s
                            .bytes()
                            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
                })
            || self.comparisons > self.options.max_comparisons
            || self.sites.len() > self.options.max_sites
            || self.complete != self.limits.is_empty()
            || self.products_determined && !self.complete
            || self.limits.len() > 8
            || self.warnings.len() > 16
            || self
                .limits
                .iter()
                .chain(&self.warnings)
                .any(|s| s.len() > 4096)
        {
            return Err(bad("invalid checked cleavage documentary evidence"));
        }
        for site in &self.sites {
            if site.enzyme.len() > 64
                || !self
                    .options
                    .selection
                    .enzymes
                    .iter()
                    .any(|e| e.enzyme.name == site.enzyme)
                || site.evidence.len() > 256
                || site.reason.len() > 4096
                || self.products_determined
                    && matches!(site.outcome, Outcome::Unknown | Outcome::Impaired)
            {
                return Err(bad("invalid checked cleavage site evidence"));
            }
        }
        bounded(self)
    }
}
