//! Bounded guide/PAM geometry and selected-record Hamming off-target search.
//! No empirical efficiency score, external aligner, bulge model or physical cleavage products.
use crate::{Alphabet, Error, Location, Molecule, Sequence, Topology, analysis};
use serde::{Deserialize, Serialize};
use std::{
    collections::{BTreeMap, BTreeSet},
    io::Write,
};
pub const ENGINE_VERSION: &str = "crispr-1";
const MAX_BYTES: usize = 16 * 1024 * 1024;
const MAX_GUIDES: usize = 4096;
const MAX_REFERENCES: usize = 64;
const MAX_BASES: usize = 20_000_000;
macro_rules! dto{($($x:tt)*)=>{#[derive(Debug,Clone,PartialEq,Eq,Serialize,Deserialize)]#[cfg_attr(feature="json-schema",derive(schemars::JsonSchema))]#[serde(deny_unknown_fields)]$($x)*};}
macro_rules! en{($($x:tt)*)=>{#[derive(Debug,Clone,Copy,PartialEq,Eq,Serialize,Deserialize)]#[cfg_attr(feature="json-schema",derive(schemars::JsonSchema))]$($x)*};}
en! {pub enum PamSide{ThreePrime,FivePrime}}
dto! {pub struct NucleaseProfile{pub id:String,pub version:String,pub guide_length:usize,pub pam:String,pub pam_side:PamSide,pub nominal_top_offset:i32,pub nominal_bottom_offset:i32,pub custom:bool,pub sources:Vec<String>}}
dto! {pub struct GuideOptions{pub profile:NucleaseProfile,pub min_gc_bp:u16,pub max_gc_bp:u16,pub max_homopolymer:usize,pub max_candidates:usize,pub max_comparisons:usize}}
impl Default for GuideOptions {
    fn default() -> Self {
        Self {
            profile: spcas9(),
            min_gc_bp: 0,
            max_gc_bp: 10_000,
            max_homopolymer: 64,
            max_candidates: 4096,
            max_comparisons: 20_000_000,
        }
    }
}
dto! {pub struct GuideCandidate{pub id:String,pub guide_sequence:String,pub pam_sequence:String,pub guide_positions:Vec<usize>,pub pam_positions:Vec<usize>,pub reverse:bool,pub pam_certain:bool,pub concrete:bool,pub gc_bp:Option<u16>,pub longest_homopolymer:Option<usize>,pub passes_filters:bool,pub nominal_top_cut:Option<usize>,pub nominal_bottom_cut:Option<usize>,pub rejections:Vec<String>}}
dto! {pub struct GuideReport{pub engine_version:String,pub source_revision:String,pub options_revision:String,pub profile_revision:String,pub batch_revision:String,pub options:GuideOptions,pub source_length:usize,pub circular:bool,pub candidates:Vec<GuideCandidate>,pub complete:bool,pub comparisons:usize,pub limits:Vec<String>,pub warnings:Vec<String>}}
dto! {pub struct ReferenceInput{pub id:String,pub molecule:Molecule,pub build:Option<String>,pub accession:Option<String>}}
dto! {pub struct ReferenceRecord{pub id:String,pub revision:String,pub name:String,pub length:usize,pub circular:bool,pub build:Option<String>,pub accession:Option<String>}}
dto! {pub struct GuideQuery{pub id:String,pub sequence:String}}
dto! {pub struct IntendedTarget{pub reference_id:String,pub reference_revision:String,pub guide_positions:Vec<usize>,pub pam_positions:Vec<usize>,pub reverse:bool}}
dto! {pub struct OffTargetOptions{pub profile:NucleaseProfile,pub max_mismatches:usize,pub intended:Option<IntendedTarget>,pub reference_scope:String,pub max_comparisons:usize,pub max_index_entries:usize,pub max_hits:usize}}
impl Default for OffTargetOptions {
    fn default() -> Self {
        Self {
            profile: spcas9(),
            max_mismatches: 3,
            intended: None,
            reference_scope: "Explicitly selected reference records".into(),
            max_comparisons: 100_000_000,
            max_index_entries: 1_000_000,
            max_hits: 4096,
        }
    }
}
dto! {pub struct OffTargetHit{pub id:String,pub reference_index:usize,pub reference_guide:String,pub pam_sequence:String,pub guide_positions:Vec<usize>,pub pam_positions:Vec<usize>,pub reverse:bool,pub min_mismatches:usize,pub max_mismatches:usize,pub pam_certain:bool,pub definite:bool,pub is_intended:bool,pub nominal_top_cut:Option<usize>,pub nominal_bottom_cut:Option<usize>}}
dto! {pub struct OffTargetReport{pub engine_version:String,pub batch_revision:String,pub options_revision:String,pub profile_revision:String,pub query:GuideQuery,pub options:OffTargetOptions,pub references:Vec<ReferenceRecord>,pub hits:Vec<OffTargetHit>,pub complete:bool,pub exhaustive_for_selected_records:bool,pub comparisons:usize,pub index_entries:usize,pub index_used:bool,pub limits:Vec<String>,pub warnings:Vec<String>}}
pub fn spcas9() -> NucleaseProfile {
    NucleaseProfile {
        id: "spcas9".into(),
        version: "jinek-2012-geometry-v1".into(),
        guide_length: 20,
        pam: "NGG".into(),
        pam_side: PamSide::ThreePrime,
        nominal_top_offset: 17,
        nominal_bottom_offset: 17,
        custom: false,
        sources: vec!["https://pmc.ncbi.nlm.nih.gov/articles/PMC6286148/".into()],
    }
}
pub fn cas12a() -> NucleaseProfile {
    NucleaseProfile {
        id: "cas12a-tttv".into(),
        version: "declared-tttv23-nominal18-23-v1".into(),
        guide_length: 23,
        pam: "TTTV".into(),
        pam_side: PamSide::FivePrime,
        nominal_top_offset: 18,
        nominal_bottom_offset: 23,
        custom: false,
        sources: vec![
            "https://pmc.ncbi.nlm.nih.gov/articles/PMC4638220/".into(),
            "https://pmc.ncbi.nlm.nih.gov/articles/PMC7286691/".into(),
        ],
    }
}
fn bad(s: &str) -> Error {
    Error::Unsupported(format!("CRISPR: {s}"))
}
fn hex(s: &str) -> bool {
    s.len() == 64
        && s.bytes()
            .all(|b| b.is_ascii_hexdigit() && !b.is_ascii_uppercase())
}
fn label(s: &str, max: usize) -> bool {
    !s.is_empty() && s.len() <= max && !s.chars().any(char::is_control)
}
struct Bound(usize);
impl Write for Bound {
    fn write(&mut self, b: &[u8]) -> std::io::Result<usize> {
        if b.len() > MAX_BYTES.saturating_sub(self.0) {
            return Err(std::io::Error::other("bounded CRISPR source/report"));
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
impl NucleaseProfile {
    pub fn validate(&self) -> Result<(), Error> {
        if !label(&self.id, 128)
            || !label(&self.version, 128)
            || !(12..=64).contains(&self.guide_length)
            || self.pam.is_empty()
            || self.pam.len() > 16
            || Sequence::dna(&self.pam)?.as_str() != self.pam
            || !(-64..=self.guide_length as i32 + 64).contains(&self.nominal_top_offset)
            || !(-64..=self.guide_length as i32 + 64).contains(&self.nominal_bottom_offset)
            || self.sources.len() > 8
            || self
                .sources
                .iter()
                .any(|s| !label(s, 2048) || !s.starts_with("https://"))
        {
            return Err(bad("invalid declared nuclease profile"));
        }
        Ok(())
    }
    fn current(&self) -> Result<(), Error> {
        self.validate()?;
        if !self.custom && *self != spcas9() && *self != cas12a() {
            return Err(bad(
                "unknown or altered built-in profile; explicitly declare custom geometry",
            ));
        }
        if self.custom && ["spcas9", "cas12a-tttv"].contains(&self.id.as_str()) {
            return Err(bad("custom geometry must use a distinct declared identity"));
        }
        Ok(())
    }
}
impl GuideOptions {
    pub fn validate(&self) -> Result<(), Error> {
        self.profile.validate()?;
        if self.min_gc_bp > self.max_gc_bp
            || self.max_gc_bp > 10_000
            || self.max_homopolymer == 0
            || self.max_homopolymer > 64
            || self.max_candidates == 0
            || self.max_candidates > MAX_GUIDES
            || self.max_comparisons == 0
            || self.max_comparisons > 200_000_000
        {
            return Err(bad("invalid guide filters/resource bounds"));
        }
        Ok(())
    }
}
impl OffTargetOptions {
    pub fn validate(&self) -> Result<(), Error> {
        self.profile.validate()?;
        if self.max_mismatches > 6
            || !label(&self.reference_scope, 4096)
            || self.max_comparisons == 0
            || self.max_comparisons > 500_000_000
            || self.max_index_entries == 0
            || self.max_index_entries > 2_000_000
            || self.max_hits == 0
            || self.max_hits > MAX_GUIDES
        {
            return Err(bad("invalid Hamming/index/hit/work bounds"));
        }
        if let Some(i) = &self.intended {
            if !label(&i.reference_id, 256)
                || !hex(&i.reference_revision)
                || i.guide_positions.len() != self.profile.guide_length
                || i.pam_positions.len() != self.profile.pam.len()
            {
                return Err(bad("invalid exact intended locus"));
            }
        }
        Ok(())
    }
}
fn location(l: &Location, depth: usize, nodes: &mut usize) -> Result<(), Error> {
    *nodes = nodes.saturating_add(1);
    if depth > 32 || *nodes > 100_000 {
        return Err(Error::TooLarge);
    }
    match l {
        Location::Complement(x) => location(x, depth + 1, nodes)?,
        Location::Join(xs) | Location::Order(xs) => {
            for x in xs {
                location(x, depth + 1, nodes)?
            }
        }
        _ => (),
    }
    Ok(())
}
fn source_shape(m: &Molecule, nodes: &mut usize) -> Result<(), Error> {
    if m.sequence.alphabet() != Alphabet::Dna {
        return Err(bad("DNA source required"));
    }
    if m.features.len() > 4096 || m.metadata.len() > 4096 || m.name.len() > 4096 {
        return Err(Error::TooLarge);
    }
    *nodes = nodes.saturating_add(1 + m.features.len() + m.metadata.len());
    for f in &m.features {
        *nodes = nodes.saturating_add(f.qualifiers.len());
        if *nodes > 100_000 {
            return Err(Error::TooLarge);
        }
        location(&f.location, 0, nodes)?;
    }
    if *nodes > 100_000 {
        return Err(Error::TooLarge);
    }
    Ok(())
}
pub fn guide_preflight(m: &Molecule, o: &GuideOptions) -> Result<(), Error> {
    o.validate()?;
    o.profile.current()?;
    source_shape(m, &mut 0)?;
    bounded(m)?;
    m.validate()?;
    Ok(())
}
pub struct ReferenceBorrow<'a> {
    pub id: &'a str,
    pub molecule: &'a Molecule,
    pub build: Option<&'a str>,
    pub accession: Option<&'a str>,
}
fn query_validate(q: &GuideQuery, o: &OffTargetOptions) -> Result<(), Error> {
    o.validate()?;
    o.profile.current()?;
    if !label(&q.id, 256)
        || q.sequence.len() != o.profile.guide_length
        || !q.sequence.bytes().all(|b| b"ACGT".contains(&b))
    {
        return Err(bad(
            "query must be a named concrete A/C/G/T guide of declared length",
        ));
    }
    Ok(())
}
pub fn reference_batch_preflight(
    q: &GuideQuery,
    refs: &[ReferenceBorrow<'_>],
    o: &OffTargetOptions,
) -> Result<(), Error> {
    query_validate(q, o)?;
    if refs.is_empty() || refs.len() > MAX_REFERENCES {
        return Err(bad("select1..64 named immutable reference records"));
    }
    let mut ids = BTreeSet::new();
    let mut nodes = 0usize;
    let mut bases = 0usize;
    let mut bytes = Bound(0);
    for r in refs {
        if !label(r.id, 256)
            || !ids.insert(r.id)
            || r.build.is_some_and(|s| !label(s, 4096))
            || r.accession.is_some_and(|s| !label(s, 4096))
        {
            return Err(bad(
                "reference roles must be unique; provenance labels bounded and explicit",
            ));
        }
        source_shape(r.molecule, &mut nodes)?;
        bases = bases.saturating_add(r.molecule.sequence.len());
        if bases > MAX_BASES {
            return Err(Error::TooLarge);
        }
        serde_json::to_writer(&mut bytes, &(r.id, r.molecule, r.build, r.accession))
            .map_err(|_| Error::TooLarge)?;
    }
    for r in refs {
        r.molecule.validate()?;
    }
    if let Some(i) = &o.intended {
        let r = refs
            .iter()
            .find(|r| r.id == i.reference_id)
            .ok_or_else(|| bad("intended role is not in selected reference collection"))?;
        if r.molecule.revision() != i.reference_revision {
            return Err(bad("intended immutable reference revision is stale"));
        }
        let start = start_for_positions(r.molecule, &o.profile, &i.guide_positions, i.reverse)?;
        let (g, p) = positions(r.molecule, &o.profile, start, i.reverse);
        if g != i.guide_positions || p != i.pam_positions {
            return Err(bad(
                "intended locus does not match exact declared guide/PAM geometry",
            ));
        }
    }
    Ok(())
}
fn borrowed(refs: &[ReferenceInput]) -> Vec<ReferenceBorrow<'_>> {
    refs.iter()
        .map(|r| ReferenceBorrow {
            id: &r.id,
            molecule: &r.molecule,
            build: r.build.as_deref(),
            accession: r.accession.as_deref(),
        })
        .collect()
}
pub fn profile_revision(p: &NucleaseProfile) -> Result<String, Error> {
    p.validate()?;
    hash(p)
}
pub fn guide_revision(m: &Molecule, o: &GuideOptions) -> Result<String, Error> {
    guide_preflight(m, o)?;
    hash(&(ENGINE_VERSION, m.revision(), o))
}
pub fn off_target_revision(
    q: &GuideQuery,
    refs: &[ReferenceInput],
    o: &OffTargetOptions,
) -> Result<String, Error> {
    if refs.is_empty() || refs.len() > MAX_REFERENCES {
        return Err(bad("select1..64 references"));
    }
    reference_batch_preflight(q, &borrowed(refs), o)?;
    hash(&(
        ENGINE_VERSION,
        q,
        refs.iter()
            .map(|r| (&r.id, r.molecule.revision(), &r.build, &r.accession))
            .collect::<Vec<_>>(),
        o,
    ))
}
fn original(n: usize, at: usize, reverse: bool) -> usize {
    let at = at % n;
    if reverse { n - 1 - at } else { at }
}
fn call(m: &Molecule, at: usize, reverse: bool) -> u8 {
    let b = m.sequence.as_str().as_bytes()[original(m.sequence.len(), at, reverse)];
    if reverse {
        crate::model::complement(b, Alphabet::Dna)
    } else {
        b
    }
}
fn offsets(p: &NucleaseProfile) -> (usize, usize) {
    match p.pam_side {
        PamSide::ThreePrime => (0, p.guide_length),
        PamSide::FivePrime => (p.pam.len(), 0),
    }
}
fn positions(
    m: &Molecule,
    p: &NucleaseProfile,
    start: usize,
    reverse: bool,
) -> (Vec<usize>, Vec<usize>) {
    let (g, pam) = offsets(p);
    (
        (0..p.guide_length)
            .map(|i| original(m.sequence.len(), start + g + i, reverse))
            .collect(),
        (0..p.pam.len())
            .map(|i| original(m.sequence.len(), start + pam + i, reverse))
            .collect(),
    )
}
fn start_for_positions(
    m: &Molecule,
    p: &NucleaseProfile,
    g: &[usize],
    reverse: bool,
) -> Result<usize, Error> {
    if g.len() != p.guide_length
        || g.iter().any(|i| *i >= m.sequence.len())
        || p.guide_length + p.pam.len() > m.sequence.len()
    {
        return Err(bad("invalid original guide call indices"));
    }
    let oriented = if reverse {
        m.sequence.len() - 1 - g[0]
    } else {
        g[0]
    };
    let (go, _) = offsets(p);
    if m.topology == Topology::Circular {
        Ok((oriented + m.sequence.len() - go) % m.sequence.len())
    } else {
        let start = oriented
            .checked_sub(go)
            .ok_or_else(|| bad("guide/PAM starts before linear source"))?;
        if start + p.guide_length + p.pam.len() > m.sequence.len() {
            return Err(bad("guide/PAM extends beyond linear source"));
        }
        Ok(start)
    }
}
fn nominal(
    m: &Molecule,
    p: &NucleaseProfile,
    start: usize,
    reverse: bool,
) -> (Option<usize>, Option<usize>) {
    let (go, _) = offsets(p);
    let n = m.sequence.len() as i64;
    let boundary = |offset: i32| {
        let oriented = start as i64 + go as i64 + i64::from(offset);
        let reference = if reverse { n - oriented } else { oriented };
        if m.topology == Topology::Circular {
            Some(reference.rem_euclid(n) as usize)
        } else {
            usize::try_from(reference)
                .ok()
                .filter(|v| *v <= m.sequence.len())
        }
    };
    if reverse {
        (
            boundary(p.nominal_bottom_offset),
            boundary(p.nominal_top_offset),
        )
    } else {
        (
            boundary(p.nominal_top_offset),
            boundary(p.nominal_bottom_offset),
        )
    }
}
fn pam_match(
    m: &Molecule,
    p: &NucleaseProfile,
    start: usize,
    reverse: bool,
    work: &mut usize,
    max: usize,
) -> Result<Option<bool>, ()> {
    let (_, offset) = offsets(p);
    let mut certain = true;
    for (j, b) in p.pam.bytes().enumerate() {
        if *work >= max {
            return Err(());
        }
        *work += 1;
        let observed = analysis::base_mask(call(m, start + offset + j, reverse));
        let expected = analysis::base_mask(b);
        if observed & expected == 0 {
            return Ok(None);
        }
        certain &= observed & !expected == 0;
    }
    Ok(Some(certain))
}
// Internal bounded-window entry point. Caller preflights the profile/source once and
// supplies only valid oriented starts; this helper grants no verification authority.
pub(crate) fn candidate_at(
    m: &Molecule,
    o: &GuideOptions,
    start: usize,
    reverse: bool,
    revision: &str,
) -> Result<Option<GuideCandidate>, Error> {
    if start
        .checked_add(o.profile.guide_length + o.profile.pam.len())
        .is_none_or(|end| end > m.sequence.len())
    {
        return Err(bad("candidate window outside bounded linear source"));
    }
    make_candidate(m, o, start, reverse, revision)
}
fn make_candidate(
    m: &Molecule,
    o: &GuideOptions,
    start: usize,
    reverse: bool,
    revision: &str,
) -> Result<Option<GuideCandidate>, Error> {
    let p = &o.profile;
    let pam_certain = match pam_match(m, p, start, reverse, &mut 0, usize::MAX) {
        Ok(Some(c)) => c,
        _ => return Ok(None),
    };
    let (go, po) = offsets(p);
    let guide_sequence = (0..p.guide_length)
        .map(|j| call(m, start + go + j, reverse) as char)
        .collect::<String>();
    let pam_sequence = (0..p.pam.len())
        .map(|j| call(m, start + po + j, reverse) as char)
        .collect::<String>();
    let concrete = guide_sequence.bytes().all(|b| b"ACGT".contains(&b));
    let gc_count = guide_sequence.bytes().filter(|b| b"GC".contains(b)).count();
    let gc_bp = concrete.then(|| (gc_count * 10_000 / p.guide_length) as u16);
    let longest_homopolymer = concrete.then(|| {
        let mut last = 0;
        let mut run = 0;
        let mut longest = 0;
        for b in guide_sequence.bytes() {
            if b == last {
                run += 1
            } else {
                run = 1;
                last = b
            }
            longest = longest.max(run);
        }
        longest
    });
    let mut rejections = vec![];
    if !concrete {
        rejections
            .push("Guide contains ambiguous calls; only concrete guides are selectable.".into())
    }
    if !pam_certain {
        rejections.push("Only a possible PAM is established by the ambiguous source.".into())
    }
    if concrete
        && (gc_count * 10_000 < usize::from(o.min_gc_bp) * p.guide_length
            || gc_count * 10_000 > usize::from(o.max_gc_bp) * p.guide_length)
    {
        rejections.push("GC fraction is outside the explicitly declared filter.".into())
    }
    if longest_homopolymer.is_some_and(|n| n > o.max_homopolymer) {
        rejections.push("Homopolymer length exceeds the explicitly declared filter.".into())
    }
    let (g, pam) = positions(m, p, start, reverse);
    let (top, bottom) = nominal(m, p, start, reverse);
    let id = hash(&(revision, profile_revision(p)?, &g, &pam, reverse))?;
    Ok(Some(GuideCandidate {
        id,
        guide_sequence,
        pam_sequence,
        guide_positions: g,
        pam_positions: pam,
        reverse,
        pam_certain,
        concrete,
        gc_bp,
        longest_homopolymer,
        passes_filters: rejections.is_empty(),
        nominal_top_cut: top,
        nominal_bottom_cut: bottom,
        rejections,
    }))
}
fn warnings() -> Vec<String> {
    vec!["Sequence/PAM geometry only; nominal cut boundaries do not establish physical cleavage products, repair outcomes or experimental activity.".into(),"Cas12a TTTV/23nt with nominal18/23 offsets is an explicitly declared planning profile. Ortholog, sequence and time-dependent cleavage heterogeneity/gaps are not represented by these nominal boundaries.".into(),"GC/homopolymer thresholds are explicit filters, not an empirical guide-efficiency score. Chromatin, expression, polymorphisms, bulges and alternative PAMs outside the profile are not assessed.".into()]
}
pub fn design_guides(m: &Molecule, o: &GuideOptions) -> Result<GuideReport, Error> {
    let batch = guide_revision(m, o)?;
    let revision = m.revision();
    let mut r = GuideReport {
        engine_version: ENGINE_VERSION.into(),
        source_revision: revision.clone(),
        options_revision: hash(o)?,
        profile_revision: profile_revision(&o.profile)?,
        batch_revision: batch,
        options: o.clone(),
        source_length: m.sequence.len(),
        circular: m.topology == Topology::Circular,
        candidates: vec![],
        complete: true,
        comparisons: 0,
        limits: vec![],
        warnings: warnings(),
    };
    let total = o.profile.guide_length + o.profile.pam.len();
    let limit = if total > m.sequence.len() {
        0
    } else if r.circular {
        m.sequence.len()
    } else {
        m.sequence.len() - total + 1
    };
    'scan: for reverse in [false, true] {
        for start in 0..limit {
            let matched = match pam_match(
                m,
                &o.profile,
                start,
                reverse,
                &mut r.comparisons,
                o.max_comparisons,
            ) {
                Ok(m) => m,
                Err(()) => {
                    r.complete = false;
                    r.limits.push("Guide/PAM comparison budget reached; absent candidates are not established.".into());
                    break 'scan;
                }
            };
            if matched.is_none() {
                continue;
            }
            if r.candidates.len() >= o.max_candidates
                || o.profile.guide_length + o.profile.pam.len()
                    > o.max_comparisons.saturating_sub(r.comparisons)
            {
                r.complete = false;
                r.limits.push(
                    "Guide evidence/work budget reached; additional candidates may be missing."
                        .into(),
                );
                break 'scan;
            }
            r.comparisons += o.profile.guide_length + o.profile.pam.len();
            r.candidates
                .push(make_candidate(m, o, start, reverse, &revision)?.unwrap());
        }
    }
    r.validate()?;
    Ok(r)
}
pub fn validate_candidate(m: &Molecule, o: &GuideOptions, c: &GuideCandidate) -> Result<(), Error> {
    guide_preflight(m, o)?;
    let start = start_for_positions(m, &o.profile, &c.guide_positions, c.reverse)?;
    if make_candidate(m, o, start, c.reverse, &m.revision())?.as_ref() != Some(c)
        || !c.passes_filters
    {
        return Err(bad(
            "candidate differs from exact eligible source/profile/filter reconstruction",
        ));
    }
    Ok(())
}
fn records(refs: &[ReferenceInput]) -> Vec<ReferenceRecord> {
    refs.iter()
        .map(|r| ReferenceRecord {
            id: r.id.clone(),
            revision: r.molecule.revision(),
            name: r.molecule.name.clone(),
            length: r.molecule.sequence.len(),
            circular: r.molecule.topology == Topology::Circular,
            build: r.build.clone(),
            accession: r.accession.clone(),
        })
        .collect()
}
fn spend(work: &mut usize, max: usize) -> Result<(), ()> {
    if *work >= max {
        return Err(());
    }
    *work += 1;
    Ok(())
}
fn code(b: u8) -> u128 {
    match b {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        _ => 3,
    }
}
fn hit_at(
    q: &GuideQuery,
    r: &ReferenceInput,
    record: &ReferenceRecord,
    reference_index: usize,
    o: &OffTargetOptions,
    start: usize,
    reverse: bool,
    work: &mut usize,
) -> Result<Option<OffTargetHit>, ()> {
    let pam_certain = match pam_match(
        &r.molecule,
        &o.profile,
        start,
        reverse,
        work,
        o.max_comparisons,
    )? {
        Some(c) => c,
        None => return Ok(None),
    };
    let (go, po) = offsets(&o.profile);
    let mut minimum = 0;
    let mut maximum = 0;
    let mut sequence = String::new();
    for (j, b) in q.sequence.bytes().enumerate() {
        spend(work, o.max_comparisons)?;
        let observed = call(&r.molecule, start + go + j, reverse);
        let mask = analysis::base_mask(observed);
        let query = analysis::base_mask(b);
        minimum += usize::from(mask & query == 0);
        maximum += usize::from(mask != query);
        if minimum > o.max_mismatches {
            return Ok(None);
        }
        sequence.push(observed as char);
    }
    let pam_sequence = (0..o.profile.pam.len())
        .map(|j| call(&r.molecule, start + po + j, reverse) as char)
        .collect::<String>();
    let (g, p) = positions(&r.molecule, &o.profile, start, reverse);
    let is_intended = o.intended.as_ref().is_some_and(|i| {
        i.reference_id == record.id
            && i.reference_revision == record.revision
            && i.reverse == reverse
            && i.guide_positions == g
            && i.pam_positions == p
    });
    let (top, bottom) = nominal(&r.molecule, &o.profile, start, reverse);
    let id = hash(&(
        &q.sequence,
        &record.id,
        &record.revision,
        profile_revision(&o.profile).map_err(|_| ())?,
        &g,
        &p,
        reverse,
    ))
    .map_err(|_| ())?;
    Ok(Some(OffTargetHit {
        id,
        reference_index,
        reference_guide: sequence,
        pam_sequence,
        guide_positions: g,
        pam_positions: p,
        reverse,
        min_mismatches: minimum,
        max_mismatches: maximum,
        pam_certain,
        definite: pam_certain && maximum <= o.max_mismatches,
        is_intended,
        nominal_top_cut: top,
        nominal_bottom_cut: bottom,
    }))
}
fn limit(r: &mut OffTargetReport, message: &str) {
    r.complete = false;
    r.exhaustive_for_selected_records = false;
    if !r.limits.iter().any(|s| s == message) {
        r.limits.push(message.into());
    }
}
pub fn search_off_targets(
    q: &GuideQuery,
    refs: &[ReferenceInput],
    o: &OffTargetOptions,
) -> Result<OffTargetReport, Error> {
    let revision = off_target_revision(q, refs, o)?;
    let mut r = OffTargetReport {
        engine_version: ENGINE_VERSION.into(),
        batch_revision: revision,
        options_revision: hash(o)?,
        profile_revision: profile_revision(&o.profile)?,
        query: q.clone(),
        options: o.clone(),
        references: records(refs),
        hits: vec![],
        complete: true,
        exhaustive_for_selected_records: true,
        comparisons: 0,
        index_entries: 0,
        index_used: false,
        limits: vec![],
        warnings: warnings(),
    };
    r.warnings.push("Coverage is ONLY the exact selected immutable records, irrespective of build/accession/scope labels. No genome-wide absence or biological specificity is established. The search model is ungapped Hamming distance plus the declared PAM; bulges, structural variation and other PAMs are excluded.".into());
    r.warnings.push("Intended status tags only an exact named role/revision/original guide+PAM coordinate/orientation tuple. Equal DNA at another locus or in another selected record remains a separate hit.".into());
    let l = q.sequence.len();
    let total = l + o.profile.pam.len();
    let parts = o.max_mismatches + 1;
    let mut seeds: BTreeMap<(usize, u128), Vec<usize>> = BTreeMap::new();
    for i in 0..parts {
        let start = i * l / parts;
        let end = (i + 1) * l / parts;
        let v = q.sequence.as_bytes()[start..end]
            .iter()
            .fold(0u128, |v, b| (v << 2) | code(*b));
        seeds.entry((end - start, v)).or_default().push(start);
    }
    let lengths = seeds.keys().map(|(k, _)| *k).collect::<BTreeSet<_>>();
    let mut candidates: BTreeSet<(usize, bool, usize)> = BTreeSet::new();
    let mut fallback = vec![];
    'index: for (ri, reference) in refs.iter().enumerate() {
        let m = &reference.molecule;
        let n = m.sequence.len();
        if total > n {
            continue;
        }
        let mut ambiguous = false;
        for b in m.sequence.as_str().bytes() {
            if spend(&mut r.comparisons, o.max_comparisons).is_err() {
                limit(
                    &mut r,
                    "Reference classification/index/search comparison budget reached; unsearched records/positions remain.",
                );
                break 'index;
            }
            if !b"ACGT".contains(&b) {
                ambiguous = true;
                break;
            }
        }
        if ambiguous {
            fallback.push(ri);
            continue;
        }
        r.index_used = true;
        for reverse in [false, true] {
            for k in &lengths {
                let mask = if *k == 64 {
                    u128::MAX
                } else {
                    (1u128 << (2 * k)) - 1
                };
                let extent = if m.topology == Topology::Circular {
                    n + k - 1
                } else {
                    n
                };
                let mut rolling = 0u128;
                for at in 0..extent {
                    if spend(&mut r.comparisons, o.max_comparisons).is_err() {
                        limit(
                            &mut r,
                            "Reference classification/index/search comparison budget reached; unsearched records/positions remain.",
                        );
                        break 'index;
                    }
                    rolling = ((rolling << 2) | code(call(m, at, reverse))) & mask;
                    if at + 1 < *k {
                        continue;
                    }
                    let seed_start = at + 1 - k;
                    if seed_start >= n {
                        continue;
                    }
                    if let Some(offsets_list) = seeds.get(&(*k, rolling)) {
                        for query_start in offsets_list {
                            if r.index_entries >= o.max_index_entries {
                                limit(
                                    &mut r,
                                    "Seed index entry budget reached; additional candidate positions may be missing.",
                                );
                                break 'index;
                            }
                            r.index_entries += 1;
                            let (go, _) = offsets(&o.profile);
                            let subtract = query_start + go;
                            let start = if m.topology == Topology::Circular {
                                (seed_start + n - subtract) % n
                            } else {
                                let Some(start) = seed_start.checked_sub(subtract) else {
                                    continue;
                                };
                                if start + total > n {
                                    continue;
                                }
                                start
                            };
                            candidates.insert((ri, reverse, start));
                        }
                    }
                }
            }
        }
    }
    let mut output_bytes = serde_json::to_vec(&r).map_err(|_| Error::TooLarge)?.len();
    'verify: for (ri, reverse, start) in candidates {
        match hit_at(
            q,
            &refs[ri],
            &r.references[ri],
            ri,
            o,
            start,
            reverse,
            &mut r.comparisons,
        ) {
            Err(()) => {
                limit(
                    &mut r,
                    "Reference classification/index/search comparison budget reached; unsearched records/positions remain.",
                );
                break;
            }
            Ok(Some(hit)) => {
                let size = serde_json::to_vec(&hit).map_err(|_| Error::TooLarge)?.len() + 1;
                if r.hits.len() >= o.max_hits
                    || output_bytes.saturating_add(size) > MAX_BYTES - 16_384
                {
                    limit(
                        &mut r,
                        "Off-target hit/output budget reached; additional hits may be missing.",
                    );
                    break 'verify;
                }
                output_bytes += size;
                r.hits.push(hit);
            }
            Ok(None) => (),
        }
    }
    'fallback: for ri in fallback {
        let m = &refs[ri].molecule;
        let n = m.sequence.len();
        let count = if m.topology == Topology::Circular {
            n
        } else {
            n - total + 1
        };
        for reverse in [false, true] {
            for start in 0..count {
                match hit_at(
                    q,
                    &refs[ri],
                    &r.references[ri],
                    ri,
                    o,
                    start,
                    reverse,
                    &mut r.comparisons,
                ) {
                    Err(()) => {
                        limit(
                            &mut r,
                            "Reference classification/index/search comparison budget reached; unsearched records/positions remain.",
                        );
                        break 'fallback;
                    }
                    Ok(Some(hit)) => {
                        let size = serde_json::to_vec(&hit).map_err(|_| Error::TooLarge)?.len() + 1;
                        if r.hits.len() >= o.max_hits
                            || output_bytes.saturating_add(size) > MAX_BYTES - 16_384
                        {
                            limit(
                                &mut r,
                                "Off-target hit/output budget reached; additional hits may be missing.",
                            );
                            break 'fallback;
                        }
                        output_bytes += size;
                        r.hits.push(hit);
                    }
                    Ok(None) => (),
                }
            }
        }
    }
    r.hits.sort_by(|a, b| {
        (a.reference_index, a.reverse, &a.guide_positions).cmp(&(
            b.reference_index,
            b.reverse,
            &b.guide_positions,
        ))
    });
    r.validate()?;
    Ok(r)
}
fn geometry(
    g: &[usize],
    p: &[usize],
    reverse: bool,
    n: usize,
    circular: bool,
    profile: &NucleaseProfile,
) -> bool {
    if n == 0
        || g.len() != profile.guide_length
        || p.len() != profile.pam.len()
        || g.len() + p.len() > n
    {
        return false;
    }
    let combined = match profile.pam_side {
        PamSide::ThreePrime => g.iter().chain(p).copied().collect::<Vec<_>>(),
        PamSide::FivePrime => p.iter().chain(g).copied().collect(),
    };
    if combined.iter().any(|i| *i >= n) {
        return false;
    }
    combined.windows(2).all(|w| {
        if reverse {
            if w[0] == 0 {
                circular && w[1] == n - 1
            } else {
                w[1] + 1 == w[0]
            }
        } else if w[0] == n - 1 {
            circular && w[1] == 0
        } else {
            w[0] + 1 == w[1]
        }
    })
}
impl GuideReport {
    pub fn validate(&self) -> Result<(), Error> {
        self.options.validate()?;
        if !label(&self.engine_version, 128)
            || ![
                &self.source_revision,
                &self.options_revision,
                &self.profile_revision,
                &self.batch_revision,
            ]
            .iter()
            .all(|s| hex(s))
            || self.options_revision != hash(&self.options)?
            || self.profile_revision != profile_revision(&self.options.profile)?
            || self.source_length == 0
            || self.source_length > crate::model::MAX_SEQUENCE_LEN
            || self.candidates.len() > self.options.max_candidates
            || self.comparisons > self.options.max_comparisons
            || self.limits.len() > 8
            || self.warnings.len() > 16
            || self
                .limits
                .iter()
                .chain(&self.warnings)
                .any(|s| s.len() > 4096)
            || self.complete != self.limits.is_empty()
        {
            return Err(bad("invalid guide report bindings/bounds"));
        }
        let mut ids = BTreeSet::new();
        for c in &self.candidates {
            if !hex(&c.id)
                || !ids.insert(&c.id)
                || !geometry(
                    &c.guide_positions,
                    &c.pam_positions,
                    c.reverse,
                    self.source_length,
                    self.circular,
                    &self.options.profile,
                )
                || c.guide_sequence.len() != self.options.profile.guide_length
                || c.pam_sequence.len() != self.options.profile.pam.len()
                || !c
                    .guide_sequence
                    .bytes()
                    .chain(c.pam_sequence.bytes())
                    .all(|b| b"ACGTRYSWKMBDHVN".contains(&b))
                || c.rejections.len() > 8
                || c.rejections.iter().any(|s| s.len() > 4096)
                || c.gc_bp.is_some_and(|g| g > 10_000)
                || c.longest_homopolymer
                    .is_some_and(|n| n == 0 || n > self.options.profile.guide_length)
                || c.concrete != c.guide_sequence.bytes().all(|b| b"ACGT".contains(&b))
                || c.concrete != (c.gc_bp.is_some() && c.longest_homopolymer.is_some())
                || c.passes_filters != c.rejections.is_empty()
                || c.passes_filters && (!c.concrete || !c.pam_certain)
                || [c.nominal_top_cut, c.nominal_bottom_cut]
                    .into_iter()
                    .flatten()
                    .any(|i| {
                        if self.circular {
                            i >= self.source_length
                        } else {
                            i > self.source_length
                        }
                    })
            {
                return Err(bad("invalid guide candidate evidence shape"));
            }
        }
        bounded(self)
    }
    pub fn validate_bindings_for(&self, m: &Molecule, o: &GuideOptions) -> Result<(), Error> {
        guide_preflight(m, o)?;
        self.validate()?;
        if self.engine_version != ENGINE_VERSION
            || self.source_revision != m.revision()
            || self.source_length != m.sequence.len()
            || self.circular != (m.topology == Topology::Circular)
            || self.options != *o
            || self.batch_revision != guide_revision(m, o)?
        {
            return Err(bad(
                "guide report current source/options/engine binding mismatch",
            ));
        }
        Ok(())
    }
}
impl OffTargetReport {
    pub fn validate(&self) -> Result<(), Error> {
        self.options.validate()?;
        if !label(&self.engine_version, 128)
            || ![
                &self.batch_revision,
                &self.options_revision,
                &self.profile_revision,
            ]
            .iter()
            .all(|s| hex(s))
            || self.options_revision != hash(&self.options)?
            || self.profile_revision != profile_revision(&self.options.profile)?
            || !label(&self.query.id, 256)
            || self.query.sequence.len() != self.options.profile.guide_length
            || !self.query.sequence.bytes().all(|b| b"ACGT".contains(&b))
            || self.references.is_empty()
            || self.references.len() > MAX_REFERENCES
            || self.hits.len() > self.options.max_hits
            || self.comparisons > self.options.max_comparisons
            || self.index_entries > self.options.max_index_entries
            || self.limits.len() > 8
            || self.warnings.len() > 16
            || self
                .limits
                .iter()
                .chain(&self.warnings)
                .any(|s| s.len() > 4096)
            || self.complete != self.limits.is_empty()
            || self.complete != self.exhaustive_for_selected_records
        {
            return Err(bad("invalid selected-reference report bindings/bounds"));
        }
        let mut roles = BTreeSet::new();
        let mut bases = 0usize;
        for r in &self.references {
            if !label(&r.id, 256)
                || !roles.insert(&r.id)
                || !hex(&r.revision)
                || r.name.len() > 4096
                || r.length == 0
                || r.length > crate::model::MAX_SEQUENCE_LEN
                || r.build.as_ref().is_some_and(|s| !label(s, 4096))
                || r.accession.as_ref().is_some_and(|s| !label(s, 4096))
            {
                return Err(bad("invalid reference descriptor"));
            }
            bases = bases.saturating_add(r.length);
        }
        if bases > MAX_BASES {
            return Err(Error::TooLarge);
        }
        if let Some(i) = &self.options.intended {
            let r = self
                .references
                .iter()
                .find(|r| r.id == i.reference_id)
                .ok_or_else(|| bad("intended role absent from report descriptors"))?;
            if r.revision != i.reference_revision
                || !geometry(
                    &i.guide_positions,
                    &i.pam_positions,
                    i.reverse,
                    r.length,
                    r.circular,
                    &self.options.profile,
                )
            {
                return Err(bad("intended descriptor revision/geometry mismatch"));
            }
        }
        let mut ids = BTreeSet::new();
        for h in &self.hits {
            let r = self
                .references
                .get(h.reference_index)
                .ok_or_else(|| bad("off-target reference index outside collection"))?;
            let intended = self.options.intended.as_ref().is_some_and(|i| {
                i.reference_id == r.id
                    && i.reference_revision == r.revision
                    && i.guide_positions == h.guide_positions
                    && i.pam_positions == h.pam_positions
                    && i.reverse == h.reverse
            });
            if !hex(&h.id)
                || !ids.insert(&h.id)
                || !geometry(
                    &h.guide_positions,
                    &h.pam_positions,
                    h.reverse,
                    r.length,
                    r.circular,
                    &self.options.profile,
                )
                || h.reference_guide.len() != self.options.profile.guide_length
                || h.pam_sequence.len() != self.options.profile.pam.len()
                || !h
                    .reference_guide
                    .bytes()
                    .chain(h.pam_sequence.bytes())
                    .all(|b| b"ACGTRYSWKMBDHVN".contains(&b))
                || h.min_mismatches > self.options.max_mismatches
                || h.min_mismatches > h.max_mismatches
                || h.max_mismatches > self.options.profile.guide_length
                || h.definite != (h.pam_certain && h.max_mismatches <= self.options.max_mismatches)
                || h.is_intended != intended
                || [h.nominal_top_cut, h.nominal_bottom_cut]
                    .into_iter()
                    .flatten()
                    .any(|i| {
                        if r.circular {
                            i >= r.length
                        } else {
                            i > r.length
                        }
                    })
            {
                return Err(bad("invalid off-target hit evidence shape"));
            }
        }
        bounded(self)
    }
    pub fn validate_bindings_for(
        &self,
        q: &GuideQuery,
        refs: &[ReferenceInput],
        o: &OffTargetOptions,
    ) -> Result<(), Error> {
        let revision = off_target_revision(q, refs, o)?;
        self.validate()?;
        if self.engine_version != ENGINE_VERSION
            || self.query != *q
            || self.options != *o
            || self.references != records(refs)
            || self.batch_revision != revision
        {
            return Err(bad(
                "off-target report current query/source/options/engine binding mismatch",
            ));
        }
        Ok(())
    }
}
pub fn replay_guides(m: &Molecule, o: &GuideOptions, r: &GuideReport) -> Result<(), Error> {
    r.validate_bindings_for(m, o)?;
    if design_guides(m, o)? != *r {
        return Err(bad("guide evidence differs from fresh source replay"));
    }
    Ok(())
}
pub fn replay_off_targets(
    q: &GuideQuery,
    refs: &[ReferenceInput],
    o: &OffTargetOptions,
    r: &OffTargetReport,
) -> Result<(), Error> {
    r.validate_bindings_for(q, refs, o)?;
    if search_off_targets(q, refs, o)? != *r {
        return Err(bad(
            "off-target evidence differs from fresh selected-record replay",
        ));
    }
    Ok(())
}
