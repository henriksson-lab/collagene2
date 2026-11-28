//! Bounded sequence-context assessment under explicit source-generation assumptions.
//! No physical modified-base propagation, reaction prediction or supplier substitution.
use crate::analysis;
use collagene_model::{Alphabet, Error, Location, Molecule, Topology};
use serde::{Deserialize, Serialize};
use std::{collections::BTreeSet, io::Write};
pub const ENGINE_VERSION: &str = "methylation-1";
pub const PROFILE_VERSION: &str = "neb-starter-2026-11-24-v1";
const TABLE: &str =
    "https://www.neb.com/en/tools-and-resources/selection-charts/dam-dcm-and-cpg-methylation";
const XBA: &str = "https://www.neb.com/faqs/is-xbai-affected-by-methylation";
const BSA: &str = "https://www.neb.com/en-gb/products/r0535-bsai";
const MAX_BYTES: usize = 8 * 1024 * 1024;
macro_rules! dto{($($x:tt)*)=>{#[derive(Debug,Clone,PartialEq,Eq,Serialize,Deserialize)]#[cfg_attr(feature="json-schema",derive(schemars::JsonSchema))]#[serde(deny_unknown_fields)]$($x)*};}
macro_rules! enumeration{($($x:tt)*)=>{#[derive(Debug,Clone,Copy,PartialEq,Eq,Serialize,Deserialize)]#[cfg_attr(feature="json-schema",derive(schemars::JsonSchema))]$($x)*};}
enumeration! {pub enum MethylationState{Unmodified,FullyModified,Unknown,HemiOrMixed}}
enumeration! {pub enum Category{Dam,Dcm,CpG}}
enumeration! {pub enum Sensitivity{NotSensitive,Blocked,Impaired,BlockedByOverlapping,BlockedBySomeOverlaps,ImpairedBySomeOverlaps,Unsupported}}
enumeration! {pub enum Outcome{PermittedUnderProfile,Blocked,Impaired,Unknown}}
dto! {pub struct MethylationOptions{pub enzymes:Vec<String>,pub dam:MethylationState,pub dcm:MethylationState,pub cpg:MethylationState,pub max_comparisons:usize,pub max_sites:usize}}
impl Default for MethylationOptions {
    fn default() -> Self {
        Self {
            enzymes: analysis::starter_enzymes()
                .into_iter()
                .map(|e| e.name)
                .collect(),
            dam: MethylationState::Unknown,
            dcm: MethylationState::Unknown,
            cpg: MethylationState::Unknown,
            max_comparisons: 20_000_000,
            max_sites: 4096,
        }
    }
}
dto! {pub struct EnzymeProfile{pub enzyme:String,pub recognition:String,pub supplier_identity:String,pub sensitivities:[Sensitivity;3],pub sources:Vec<String>,pub notes:Vec<String>}}
dto! {pub struct ModifiedBase{pub position:usize,pub reverse:bool,pub base:char}}
dto! {pub struct MotifOverlap{pub motif:String,pub start:usize,pub positions:Vec<usize>,pub certain:bool,pub modified_bases:Vec<ModifiedBase>}}
dto! {pub struct SiteEffect{pub category:Category,pub state:MethylationState,pub sensitivity:Sensitivity,pub outcome:Outcome,pub overlaps:Vec<MotifOverlap>,pub reason:String}}
dto! {pub struct MethylationSite{pub index:usize,pub enzyme:String,pub recognition_start:usize,pub recognition_length:usize,pub reverse:bool,pub certain:bool,pub top_cut:Option<usize>,pub bottom_cut:Option<usize>,pub effects:Vec<SiteEffect>,pub outcome:Outcome}}
dto! {pub struct MethylationReport{pub engine_version:String,pub profile_version:String,pub profile_revision:String,pub source_revision:String,pub batch_revision:String,pub options_revision:String,pub options:MethylationOptions,pub source_length:usize,pub circular:bool,pub profiles:Vec<EnzymeProfile>,pub sites:Vec<MethylationSite>,pub complete:bool,pub comparisons:usize,pub limits:Vec<String>,pub warnings:Vec<String>}}
pub fn profiles() -> Vec<EnzymeProfile> {
    use Sensitivity::*;
    analysis::starter_enzymes().into_iter().map(|e|{let mut p=EnzymeProfile{enzyme:e.name.clone(),recognition:e.recognition,supplier_identity:format!("NEB {} (non-HF starter identity)",e.name),sensitivities:[NotSensitive;3],sources:vec![TABLE.into()],notes:vec![]};match e.name.as_str(){"EcoRI"=>p.sensitivities[2]=BlockedBySomeOverlaps,"EcoRV"=>p.sensitivities[2]=ImpairedBySomeOverlaps,"SmaI"|"NotI"=>p.sensitivities[2]=Blocked,"XbaI"=>{p.sensitivities[0]=BlockedByOverlapping;p.sources.push(XBA.into());},"BsaI"=>{p.supplier_identity="NEB BsaI R0535 (discontinued; not BsaI-HFv2)".into();p.sensitivities=[NotSensitive,Unsupported,BlockedBySomeOverlaps];p.sources=vec![BSA.into()];p.notes.push("Generic R0535 product sensitivity section reports some-overlap Dcm impairment, while its note reports blocking. This conflict is retained as unsupported; HFv2 data are not substituted.".into());},"BsmBI"=>{p.supplier_identity="Generic BsmBI starter identity; no formulation-specific profile verified".into();p.sensitivities=[Unsupported;3];p.notes.push("Current NEB table identifies BsmBI-v2. Its formulation is not silently substituted for generic BsmBI; sensitivity remains unsupported.".into());},_=>()}p}).collect()
}
fn err(s: &str) -> Error {
    Error::Unsupported(format!("methylation assessment: {s}"))
}
fn hash<T: Serialize>(v: &T) -> Result<String, Error> {
    Ok(
        blake3::hash(&serde_json::to_vec(v).map_err(|_| Error::TooLarge)?)
            .to_hex()
            .to_string(),
    )
}
fn hex(s: &str) -> bool {
    s.len() == 64
        && s.bytes()
            .all(|b| b.is_ascii_hexdigit() && !b.is_ascii_uppercase())
}
struct Bound(usize);
impl Write for Bound {
    fn write(&mut self, b: &[u8]) -> std::io::Result<usize> {
        if b.len() > MAX_BYTES.saturating_sub(self.0) {
            return Err(std::io::Error::other("bounded source/report"));
        }
        self.0 += b.len();
        Ok(b.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}
fn bounded<T: Serialize>(v: &T) -> Result<(), Error> {
    serde_json::to_writer(Bound(0), v).map_err(|_| Error::TooLarge)
}
fn location(l: &Location, depth: usize, n: &mut usize) -> Result<(), Error> {
    *n = n.saturating_add(1);
    if depth > 32 || *n > 100_000 {
        return Err(Error::TooLarge);
    }
    match l {
        Location::Complement(x) => location(x, depth + 1, n)?,
        Location::Join(xs) | Location::Order(xs) => {
            for x in xs {
                location(x, depth + 1, n)?
            }
        }
        _ => (),
    }
    Ok(())
}
impl MethylationOptions {
    pub fn validate(&self) -> Result<(), Error> {
        if self.enzymes.is_empty()
            || self.enzymes.len() > 12
            || self.max_sites == 0
            || self.max_sites > 4096
            || self.max_comparisons == 0
            || self.max_comparisons > 100_000_000
        {
            return Err(err("invalid enzyme count/site/work bounds"));
        }
        let known = profiles();
        let mut ids = BTreeSet::new();
        for name in &self.enzymes {
            if name.len() > 32
                || !ids.insert(name.to_ascii_lowercase())
                || !known.iter().any(|p| p.enzyme.eq_ignore_ascii_case(name))
            {
                return Err(err(
                    "enzyme IDs must be distinct members of the starter profile",
                ));
            }
        }
        Ok(())
    }
}
pub fn source_preflight(m: &Molecule, o: &MethylationOptions) -> Result<(), Error> {
    o.validate()?;
    if m.sequence.alphabet() != Alphabet::Dna {
        return Err(err("DNA required"));
    }
    if m.sequence.len() > 1_000_000
        || m.features.len() > 4096
        || m.metadata.len() > 4096
        || m.name.len() > 4096
    {
        return Err(Error::TooLarge);
    }
    let mut nodes = 1 + m.features.len() + m.metadata.len();
    for f in &m.features {
        nodes = nodes.saturating_add(f.qualifiers.len());
        if nodes > 100_000 {
            return Err(Error::TooLarge);
        }
        location(&f.location, 0, &mut nodes)?;
    }
    bounded(m)?;
    m.validate()?;
    Ok(())
}
pub fn options_revision(o: &MethylationOptions) -> Result<String, Error> {
    o.validate()?;
    hash(o)
}
pub fn profile_revision() -> String {
    hash(&(PROFILE_VERSION, profiles())).expect("fixed bounded profile")
}
pub fn batch_revision(m: &Molecule, o: &MethylationOptions) -> Result<String, Error> {
    source_preflight(m, o)?;
    hash(&(
        ENGINE_VERSION,
        profile_revision(),
        m.revision(),
        options_revision(o)?,
    ))
}
fn state(o: &MethylationOptions, c: Category) -> MethylationState {
    match c {
        Category::Dam => o.dam,
        Category::Dcm => o.dcm,
        Category::CpG => o.cpg,
    }
}
fn compare(
    bases: &[u8],
    start: usize,
    pattern: &[u8],
    circular: bool,
    work: &mut usize,
    max: usize,
) -> Result<Option<bool>, ()> {
    if pattern.len() > bases.len() || (!circular && start + pattern.len() > bases.len()) {
        return Ok(None);
    }
    let mut certain = true;
    for (j, p) in pattern.iter().enumerate() {
        if *work >= max {
            return Err(());
        }
        *work += 1;
        let b = analysis::base_mask(bases[(start + j) % bases.len()]);
        let mask = analysis::base_mask(*p);
        if b & mask == 0 {
            return Ok(None);
        }
        certain &= b & !mask == 0;
    }
    Ok(Some(certain))
}
fn context(
    m: &Molecule,
    start: usize,
    len: usize,
    c: Category,
    work: &mut usize,
    max: usize,
) -> Result<Vec<MotifOverlap>, ()> {
    let (pattern, modified) = match c {
        Category::Dam => ("GATC", 1),
        Category::Dcm => ("CCWGG", 1),
        Category::CpG => ("CG", 0),
    };
    let n = m.sequence.len();
    if pattern.len() > n {
        return Ok(vec![]);
    }
    let circle = m.topology == Topology::Circular;
    let mut starts = BTreeSet::new();
    for relative in -(pattern.len() as i64 - 1)..len as i64 {
        let at = start as i64 + relative;
        if circle {
            starts.insert(at.rem_euclid(n as i64) as usize);
        } else if at >= 0 && (at as usize) + pattern.len() <= n {
            starts.insert(at as usize);
        }
    }
    let mut out = vec![];
    for at in starts {
        if let Some(certain) = compare(
            m.sequence.as_str().as_bytes(),
            at,
            pattern.as_bytes(),
            circle,
            work,
            max,
        )? {
            let positions = (0..pattern.len()).map(|j| (at + j) % n).collect::<Vec<_>>();
            out.push(MotifOverlap {
                motif: pattern.into(),
                start: at,
                modified_bases: vec![
                    ModifiedBase {
                        position: positions[modified],
                        reverse: false,
                        base: if c == Category::Dam { 'A' } else { 'C' },
                    },
                    ModifiedBase {
                        position: positions[pattern.len() - 1 - modified],
                        reverse: true,
                        base: if c == Category::Dam { 'A' } else { 'C' },
                    },
                ],
                positions,
                certain,
            });
        }
    }
    Ok(out)
}
fn effect(
    c: Category,
    s: MethylationState,
    sensitivity: Sensitivity,
    overlaps: Vec<MotifOverlap>,
) -> SiteEffect {
    use Outcome::*;
    let (outcome, reason) = if s == MethylationState::HemiOrMixed {
        (
            Unknown,
            "Hemi/mixed modification is not modeled; no complete-cleavage claim.",
        )
    } else if s == MethylationState::Unmodified {
        (
            PermittedUnderProfile,
            "Explicitly unmodified for this modification category.",
        )
    } else if sensitivity == Sensitivity::Unsupported {
        (
            Unknown,
            "No supported non-conflicting formulation-specific sensitivity rule.",
        )
    } else if sensitivity == Sensitivity::NotSensitive {
        (
            PermittedUnderProfile,
            "Published profile reports no sensitivity to this category.",
        )
    } else if overlaps.is_empty() {
        (
            PermittedUnderProfile,
            "No overlapping methylase motif under the bounded current-source context model.",
        )
    } else if s == MethylationState::Unknown {
        (
            Unknown,
            "A relevant motif is possible but source modification state is unknown.",
        )
    } else if overlaps.iter().all(|x| !x.certain) {
        (
            Unknown,
            "Only possible motif overlaps exist in ambiguous sequence context.",
        )
    } else {
        match sensitivity {
            Sensitivity::Blocked | Sensitivity::BlockedByOverlapping => (
                Blocked,
                "Certain overlapping motif under explicit full modification assumption.",
            ),
            Sensitivity::Impaired => (
                Impaired,
                "Published impairment is not converted into complete cleavage or zero cleavage.",
            ),
            _ => (
                Unknown,
                "Published rule affects some overlap combinations; exact affected combinations are not established by this profile.",
            ),
        }
    };
    SiteEffect {
        category: c,
        state: s,
        sensitivity,
        outcome,
        overlaps,
        reason: reason.into(),
    }
}
fn aggregate(xs: impl Iterator<Item = Outcome>) -> Outcome {
    let mut outcome = Outcome::PermittedUnderProfile;
    for x in xs {
        if x == Outcome::Blocked {
            return x;
        }
        if x == Outcome::Unknown {
            outcome = x
        } else if x == Outcome::Impaired && outcome != Outcome::Unknown {
            outcome = x
        }
    }
    outcome
}
fn cuts(
    m: &Molecule,
    start: usize,
    len: usize,
    rev: bool,
    top: i32,
    bottom: i32,
) -> (Option<usize>, Option<usize>) {
    let n = m.sequence.len() as i64;
    let (a, b) = if rev {
        (
            start as i64 + len as i64 - i64::from(bottom),
            start as i64 + len as i64 - i64::from(top),
        )
    } else {
        (
            start as i64 + i64::from(top),
            start as i64 + i64::from(bottom),
        )
    };
    let cv = |p: i64| {
        if m.topology == Topology::Circular {
            Some(p.rem_euclid(n) as usize)
        } else if (0..=n).contains(&p) {
            Some(p as usize)
        } else {
            None
        }
    };
    (cv(a), cv(b))
}
pub fn assess(m: &Molecule, o: &MethylationOptions) -> Result<MethylationReport, Error> {
    let batch = batch_revision(m, o)?;
    let selected = profiles()
        .into_iter()
        .filter(|p| o.enzymes.iter().any(|s| s.eq_ignore_ascii_case(&p.enzyme)))
        .collect::<Vec<_>>();
    let mut r=MethylationReport{engine_version:ENGINE_VERSION.into(),profile_version:PROFILE_VERSION.into(),profile_revision:profile_revision(),source_revision:m.revision(),batch_revision:batch,options_revision:options_revision(o)?,options:o.clone(),source_length:m.sequence.len(),circular:m.topology==Topology::Circular,profiles:selected,sites:vec![],complete:true,comparisons:0,limits:vec![],warnings:vec!["Planning assessment under explicit source-generation assumptions, not measured modification evidence, physical modified-base propagation or a prediction of experimental cleavage.".into(),"Only Dam N6-adenine at GATC, Dcm 5-methylcytosine at CCWGG and CpG 5-methylcytosine are scoped. Other modifications, supplier conditions, accessibility, kinetics and star activity are not evaluated.".into(),"PermittedUnderProfile addresses modeled methylation only. Recognition certainty and available strand-cut substrate are separate evidence.".into(),"No complete digest products or modification propagation are inferred; old physical digest behavior is unchanged.".into()]};
    let enzymes = analysis::starter_enzymes();
    let mut evidence_bytes = 0usize;
    'scan: for p in r.profiles.clone() {
        let enzyme = enzymes.iter().find(|e| e.name == p.enzyme).unwrap();
        let rc = collagene_model::Sequence::dna(&p.recognition)?
            .reverse_complement()?
            .as_str()
            .as_bytes()
            .to_vec();
        let mut patterns = vec![(p.recognition.as_bytes(), false)];
        if rc != p.recognition.as_bytes() {
            patterns.push((&rc, true));
        }
        for (pattern, reverse) in patterns {
            let n = m.sequence.len();
            let limit = if pattern.len() > n {
                0
            } else if r.circular {
                n
            } else {
                n - pattern.len() + 1
            };
            for at in 0..limit {
                let matched = match compare(
                    m.sequence.as_str().as_bytes(),
                    at,
                    pattern,
                    r.circular,
                    &mut r.comparisons,
                    o.max_comparisons,
                ) {
                    Ok(x) => x,
                    Err(()) => {
                        r.complete = false;
                        r.limits.push("Recognition/context comparison budget reached; absent sites are not established.".into());
                        break 'scan;
                    }
                };
                let Some(certain) = matched else { continue };
                if r.sites.len() >= o.max_sites {
                    r.complete = false;
                    r.limits.push(
                        "Site evidence budget reached; additional sites may be missing.".into(),
                    );
                    break 'scan;
                }
                let mut effects = vec![];
                for (i, c) in [Category::Dam, Category::Dcm, Category::CpG]
                    .into_iter()
                    .enumerate()
                {
                    let overlaps = match context(
                        m,
                        at,
                        pattern.len(),
                        c,
                        &mut r.comparisons,
                        o.max_comparisons,
                    ) {
                        Ok(x) => x,
                        Err(()) => {
                            r.complete = false;
                            r.limits.push("Recognition/context comparison budget reached; partially examined site omitted.".into());
                            break 'scan;
                        }
                    };
                    effects.push(effect(c, state(o, c), p.sensitivities[i], overlaps));
                }
                let outcome = if certain {
                    aggregate(effects.iter().map(|e| e.outcome))
                } else {
                    Outcome::Unknown
                };
                let (top_cut, bottom_cut) = cuts(
                    m,
                    at,
                    pattern.len(),
                    reverse,
                    enzyme.top_cut,
                    enzyme.bottom_cut,
                );
                let site = MethylationSite {
                    index: r.sites.len(),
                    enzyme: p.enzyme.clone(),
                    recognition_start: at,
                    recognition_length: pattern.len(),
                    reverse,
                    certain,
                    top_cut,
                    bottom_cut,
                    effects,
                    outcome,
                };
                let size = serde_json::to_vec(&site)
                    .map_err(|_| Error::TooLarge)?
                    .len()
                    + 1;
                if evidence_bytes.saturating_add(size) > MAX_BYTES - 65_536 {
                    r.complete = false;
                    r.limits.push(
                        "Serialized evidence budget reached; additional sites may be missing."
                            .into(),
                    );
                    break 'scan;
                }
                evidence_bytes += size;
                r.sites.push(site);
            }
        }
    }
    r.validate()?;
    Ok(r)
}
impl MethylationReport {
    /// Bounded shape/provenance validation; imported evidence is never current science.
    pub fn validate(&self) -> Result<(), Error> {
        self.options.validate()?;
        if self.engine_version.is_empty()
            || self.engine_version.len() > 128
            || self.profile_version.is_empty()
            || self.profile_version.len() > 128
            || ![
                &self.profile_revision,
                &self.source_revision,
                &self.batch_revision,
                &self.options_revision,
            ]
            .iter()
            .all(|s| hex(s))
            || self.options_revision != options_revision(&self.options)?
            || self.source_length == 0
            || self.source_length > 1_000_000
            || self.profiles.len() != self.options.enzymes.len()
            || self.sites.len() > self.options.max_sites
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
            return Err(err("invalid bounded report header"));
        }
        let mut names = BTreeSet::new();
        for p in &self.profiles {
            if p.enzyme.len() > 32
                || !names.insert(&p.enzyme)
                || !self
                    .options
                    .enzymes
                    .iter()
                    .any(|s| s.eq_ignore_ascii_case(&p.enzyme))
                || p.recognition.len() > 16
                || p.recognition.is_empty()
                || !p.recognition.bytes().all(|b| b"ACGT".contains(&b))
                || p.supplier_identity.is_empty()
                || p.supplier_identity.len() > 4096
                || p.sources.len() > 8
                || p.sources.is_empty()
                || p.sources
                    .iter()
                    .any(|s| s.len() > 4096 || !s.starts_with("https://"))
                || p.notes.len() > 8
                || p.notes.iter().any(|s| s.len() > 4096)
            {
                return Err(err("invalid profile provenance"));
            }
        }
        let mut seen = BTreeSet::new();
        for (i, s) in self.sites.iter().enumerate() {
            let p = self
                .profiles
                .iter()
                .find(|p| p.enzyme == s.enzyme)
                .ok_or_else(|| err("site profile missing"))?;
            let l = s.recognition_length;
            if s.index != i
                || l != p.recognition.len()
                || l > self.source_length
                || s.recognition_start >= self.source_length
                || (!self.circular && s.recognition_start + l > self.source_length)
                || s.top_cut.is_some_and(|c| {
                    if self.circular {
                        c >= self.source_length
                    } else {
                        c > self.source_length
                    }
                })
                || s.bottom_cut.is_some_and(|c| {
                    if self.circular {
                        c >= self.source_length
                    } else {
                        c > self.source_length
                    }
                })
                || s.effects.len() != 3
                || !seen.insert((&s.enzyme, s.recognition_start, s.reverse))
            {
                return Err(err("invalid site geometry"));
            }
            for (j, e) in s.effects.iter().enumerate() {
                let category = [Category::Dam, Category::Dcm, Category::CpG][j];
                let (motif, offset) = match category {
                    Category::Dam => ("GATC", 1),
                    Category::Dcm => ("CCWGG", 1),
                    Category::CpG => ("CG", 0),
                };
                if e.category != category
                    || e.state != state(&self.options, category)
                    || e.sensitivity != p.sensitivities[j]
                    || e.reason.is_empty()
                    || e.reason.len() > 4096
                    || e.overlaps.len() > 32
                {
                    return Err(err("invalid category evidence"));
                }
                let mut starts = BTreeSet::new();
                for overlap in &e.overlaps {
                    if overlap.motif != motif
                        || overlap.positions.len() != motif.len()
                        || motif.len() > self.source_length
                        || overlap.start >= self.source_length
                        || (!self.circular && overlap.start + motif.len() > self.source_length)
                        || !starts.insert(overlap.start)
                        || overlap
                            .positions
                            .iter()
                            .enumerate()
                            .any(|(at, pos)| *pos != (overlap.start + at) % self.source_length)
                        || !overlap.positions.iter().any(|pos| {
                            (0..l).any(|at| *pos == (s.recognition_start + at) % self.source_length)
                        })
                    {
                        return Err(err("invalid motif coordinates"));
                    }
                    let base = if category == Category::Dam { 'A' } else { 'C' };
                    if overlap.modified_bases
                        != vec![
                            ModifiedBase {
                                position: overlap.positions[offset],
                                reverse: false,
                                base,
                            },
                            ModifiedBase {
                                position: overlap.positions[motif.len() - 1 - offset],
                                reverse: true,
                                base,
                            },
                        ]
                    {
                        return Err(err("invalid strand/base modification evidence"));
                    }
                }
                let expected = effect(category, e.state, e.sensitivity, e.overlaps.clone());
                if e.outcome != expected.outcome {
                    return Err(err(
                        "category outcome contradicts disclosed state/rule evidence",
                    ));
                }
            }
            let expected = if s.certain {
                aggregate(s.effects.iter().map(|e| e.outcome))
            } else {
                Outcome::Unknown
            };
            if s.outcome != expected {
                return Err(err("aggregate outcome contradicts site evidence"));
            }
        }
        bounded(self)
    }
    /// Current source/options/catalog binding only: does not rescan source sequence.
    pub fn validate_bindings_for(&self, m: &Molecule, o: &MethylationOptions) -> Result<(), Error> {
        source_preflight(m, o)?;
        self.validate()?;
        let selected = profiles()
            .into_iter()
            .filter(|p| o.enzymes.iter().any(|s| s.eq_ignore_ascii_case(&p.enzyme)))
            .collect::<Vec<_>>();
        if self.engine_version != ENGINE_VERSION
            || self.profile_version != PROFILE_VERSION
            || self.profile_revision != profile_revision()
            || self.profiles != selected
            || self.source_revision != m.revision()
            || self.source_length != m.sequence.len()
            || self.circular != (m.topology == Topology::Circular)
            || self.options != *o
            || self.options_revision != options_revision(o)?
            || self.batch_revision != batch_revision(m, o)?
        {
            return Err(err(
                "source/options/current engine/profile binding mismatch",
            ));
        }
        Ok(())
    }
}
/// Fresh exact computational replay, for auditing imported evidence.
pub fn replay(m: &Molecule, o: &MethylationOptions, r: &MethylationReport) -> Result<(), Error> {
    r.validate_bindings_for(m, o)?;
    if assess(m, o)? != *r {
        return Err(err("report differs from fresh source-context assessment"));
    }
    Ok(())
}
