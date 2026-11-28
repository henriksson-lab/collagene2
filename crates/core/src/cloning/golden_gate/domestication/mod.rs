//! Source/CDS-bound substitution intent, independent of physical preparation.
//! Native checks rescan the full DNA and every affected CDS. No imported report
//! constructs CheckedIntent; hosts additionally attest exact source events.
use crate::{
    cloning::{
        pcr::mutation_intent,
        restriction::catalog::{self as catalog_digest, CatalogSelection},
    },
    sequence::{analysis, codon_design, translation},
};
use collagene_model::{Alphabet, Error, Location, Molecule, Sequence, Topology};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};
pub const ENGINE_VERSION: &str = "source-cds-domestication-4";
pub const MAX_SITES: usize = 4096;
pub const MAX_FEATURES: usize = 256;
pub const MAX_PATCHES: usize = 256;
pub const MAX_WORK: u64 = 200_000_000;
const MAX_MAPPED: usize = 1_000_000;
macro_rules! dto {($($t:tt)*)=>{#[derive(Clone,Debug,Serialize,Deserialize)]#[cfg_attr(feature="json-schema",derive(schemars::JsonSchema))]#[serde(deny_unknown_fields)]$($t)*};}
dto! {#[doc="Zero-based, half-open original source coordinates; circular end may unroll at most one turn."] #[derive(PartialEq,Eq)]pub struct SourceRange{pub start:usize,pub end:usize}}
dto! {#[derive(PartialEq,Eq,PartialOrd,Ord)]pub struct SitePin{pub enzyme_id:String,pub start:usize,pub reverse:bool}}
dto! {#[derive(PartialEq,Eq)]pub struct FeaturePin{pub index:usize,pub revision:String}}
dto! {pub struct Limits{pub max_patches:usize,pub max_work:u64}}
impl Default for Limits {
    fn default() -> Self {
        Self {
            max_patches: 64,
            max_work: 50_000_000,
        }
    }
}
dto! {pub struct Request{
    pub source_revision:String,pub selection:CatalogSelection,pub scope:Vec<SourceRange>,
    pub retained_sites:Vec<SitePin>,pub coding_features:Vec<FeaturePin>,pub protected:Vec<SourceRange>,
    pub options:codon_design::DesignOptions,pub limits:Limits,
}}
dto! {#[derive(PartialEq,Eq)]pub struct Patch{#[doc="Zero-based original source position."] pub position:usize,pub expected:String,pub replacement:String}}
dto! {pub struct Selection{pub context_revision:String,pub candidate_id:String,pub patches:Vec<Patch>}}
dto! {pub struct Site{pub pin:SitePin,pub recognition:String,pub in_scope:bool,pub retained:bool}}
dto! {pub struct CdsAssessment{pub pin:FeaturePin,pub label:String,pub selected:bool,pub supported:bool,
pub reason:Option<String>,pub codon_start:Option<u8>,pub genetic_code:Option<u8>,pub has_terminal_stop:Option<bool>,pub protein:Option<String>}}
dto! {pub struct Inspection{pub context_revision:String,pub request:Request,pub sites:Vec<Site>,pub features:Vec<CdsAssessment>,pub blocking:Vec<String>,pub warnings:Vec<String>}}
dto! {pub enum ProposalStatus{CandidateFound,Blocked,SearchIncomplete,NoChangeRequired}}
dto! {pub struct FeatureSearch{pub feature_index:usize,pub status:codon_design::DesignStatus,pub search_complete:bool,pub work_used:u64,pub reserved_work:u64,pub minimum_changes_proven:bool}}
dto! {pub struct Proposal{pub context_revision:String,pub inspection:Inspection,pub candidate:Option<Selection>,pub check:Option<CheckReport>,pub status:ProposalStatus,pub search_complete:bool,pub searches:Vec<FeatureSearch>,pub warnings:Vec<String>}}
dto! {pub struct FeatureCheck{pub pin:FeaturePin,pub label:String,pub before_protein:String,pub after_protein:String,pub codon_start:u8,pub genetic_code:u8,pub has_terminal_stop:bool}}
dto! {pub struct CheckReport{pub context_revision:String,pub request:Request,pub selection:Selection,pub source_revision:String,pub target_revision:String,pub candidate_id:String,
pub sites_before:Vec<Site>,pub sites_after:Vec<Site>,pub features:Vec<FeatureCheck>,pub warnings:Vec<String>}}
/// Fresh private reference proof. It conveys no material/chemical authority.
pub struct CheckedIntent {
    source: Molecule,
    target: Molecule,
    report: CheckReport,
}
impl CheckedIntent {
    pub fn source(&self) -> &Molecule {
        &self.source
    }
    pub fn target(&self) -> &Molecule {
        &self.target
    }
    pub fn source_revision(&self) -> &str {
        &self.report.source_revision
    }
    pub fn request(&self) -> &Request {
        &self.report.request
    }
    pub fn patches(&self) -> &[Patch] {
        &self.report.selection.patches
    }
    pub fn report(&self) -> &CheckReport {
        &self.report
    }
    pub fn context_revision(&self) -> &str {
        &self.report.context_revision
    }
}
fn bad(s: &str) -> Error {
    Error::Unsupported(format!("domestication: {s}"))
}
fn hex(s: &str) -> bool {
    s.len() == 64
        && s.bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}
fn concrete(s: &str) -> bool {
    !s.is_empty() && s.bytes().all(|b| b"ACGT".contains(&b))
}
fn ranges(ranges: &[SourceRange], source: &Molecule) -> Result<(), Error> {
    let n = source.sequence.len();
    if ranges.len() > 256 {
        return Err(Error::TooLarge);
    }
    for r in ranges {
        if r.start >= n
            || r.end <= r.start
            || r.end - r.start > n
            || source.topology == Topology::Linear && r.end > n
        {
            return Err(bad("range exceeds original source or one circular turn"));
        }
    }
    Ok(())
}
fn contains(ranges: &[SourceRange], p: usize, n: usize) -> bool {
    ranges.iter().any(|r| {
        let q = if p < r.start { p + n } else { p };
        q >= r.start && q < r.end
    })
}
fn positions_of(sites: impl Iterator<Item = Site>, n: usize) -> BTreeSet<usize> {
    sites
        .flat_map(|s| (0..s.recognition.len()).map(move |i| (s.pin.start + i) % n))
        .collect()
}
fn request_shape(r: &Request) -> Result<(), Error> {
    if !hex(&r.source_revision)
        || r.scope.is_empty()
        || r.scope.len() > 256
        || r.protected.len() > 256
        || r.retained_sites.len() > MAX_SITES
        || r.coding_features.len() > MAX_FEATURES
        || r.limits.max_patches == 0
        || r.limits.max_patches > MAX_PATCHES
        || r.limits.max_work == 0
        || r.limits.max_work > MAX_WORK
        || r.options.max_work == 0
        || r.options.max_work > MAX_WORK
        || r.options.left_context.len() > 128
        || r.options.right_context.len() > 128
        || !r.options.left_context.is_empty()
        || !r.options.right_context.is_empty()
    {
        return Err(bad(
            "invalid bounded policy; source context is derived, never caller-substituted",
        ));
    }
    codon_design::options_revision(&r.options)?;
    r.selection.validate()?;
    if r.selection.enzymes.len() > 16 {
        return Err(Error::TooLarge);
    }
    let mut features = BTreeSet::new();
    for pin in &r.coding_features {
        if !hex(&pin.revision) || !features.insert(pin.index) {
            return Err(bad("invalid or duplicate selected feature pin"));
        }
    }
    let mut sites = BTreeSet::new();
    for pin in &r.retained_sites {
        if pin.enzyme_id.len() > 128 || !sites.insert(pin) {
            return Err(bad("invalid or duplicate retained site"));
        }
    }
    crate::encoding::bounded_json::hash(r)?;
    Ok(())
}
pub fn feature_pin(source: &Molecule, index: usize) -> Result<FeaturePin, Error> {
    let feature = source
        .features
        .get(index)
        .ok_or_else(|| bad("feature index is not retained"))?;
    Ok(FeaturePin {
        index,
        revision: crate::encoding::bounded_json::hash(feature)?,
    })
}
/// Borrowed guards before scans or source copies; host event ownership is separate.
pub fn preflight(source: &Molecule, r: &Request) -> Result<(), Error> {
    request_shape(r)?;
    collagene_model::bounds::model_preflight(&[source])?;
    source.validate()?;
    if source.sequence.alphabet() != Alphabet::Dna
        || !concrete(source.sequence.as_str())
        || source.revision() != r.source_revision
    {
        return Err(bad("exact concrete original DNA/source revision required"));
    }
    ranges(&r.scope, source)?;
    ranges(&r.protected, source)?;
    r.selection
        .validate_for(&catalog_digest::builtin_catalog()?)?;
    for pin in &r.coding_features {
        if feature_pin(source, pin.index)? != *pin
            || !source.features[pin.index].kind.eq_ignore_ascii_case("CDS")
        {
            return Err(bad("selected CDS/qualifier pin changed"));
        }
    }
    base_work(source, r)?;
    Ok(())
}
fn context(tag: &str, r: &Request) -> Result<String, Error> {
    request_shape(r)?;
    crate::encoding::bounded_json::hash(&(ENGINE_VERSION, tag, r))
}
pub fn inspect_revision(source: &Molecule, r: &Request) -> Result<String, Error> {
    preflight(source, r)?;
    context("inspect", r)
}
pub fn propose_revision(source: &Molecule, r: &Request) -> Result<String, Error> {
    preflight(source, r)?;
    context("propose", r)
}
fn patch_shape(r: &Request, patches: &[Patch]) -> Result<(), Error> {
    if patches.is_empty() || patches.len() > r.limits.max_patches {
        return Err(bad("choose a nonempty bounded literal substitution list"));
    }
    for p in patches {
        if p.expected.len() != 1
            || p.replacement.len() != 1
            || !concrete(&p.expected)
            || !concrete(&p.replacement)
            || p.expected == p.replacement
        {
            return Err(bad(
                "patches require different concrete single original/replacement bases",
            ));
        }
    }
    if patches.windows(2).any(|w| w[0].position >= w[1].position) {
        return Err(bad(
            "patches require strictly increasing unique source positions",
        ));
    }
    Ok(())
}
pub fn candidate_id(r: &Request, patches: &[Patch]) -> Result<String, Error> {
    request_shape(r)?;
    patch_shape(r, patches)?;
    crate::encoding::bounded_json::hash(&(ENGINE_VERSION, "candidate", r, patches))
}
pub fn check_revision(source: &Molecule, r: &Request, s: &Selection) -> Result<String, Error> {
    preflight(source, r)?;
    patch_shape(r, &s.patches)?;
    if s.context_revision != context("propose", r)?
        || s.candidate_id != candidate_id(r, &s.patches)?
    {
        return Err(bad("explicit candidate/source context changed"));
    }
    if s.patches
        .iter()
        .any(|p| p.position >= source.sequence.len())
    {
        return Err(bad("substitution lies outside original source"));
    }
    crate::encoding::bounded_json::hash(&(ENGINE_VERSION, "check", r, s))
}
fn scan(source: &Molecule, r: &Request) -> Result<Vec<Site>, Error> {
    let n = source.sequence.len();
    let mut sites = vec![];
    for e in &r.selection.enzymes {
        for hit in analysis::find_matches(source, &e.enzyme.recognition, true)? {
            if sites.len() >= MAX_SITES {
                return Err(Error::TooLarge);
            }
            if !hit.certain {
                return Err(bad(
                    "uncertain recognition occurrence cannot be domesticated as exact",
                ));
            }
            let pin = SitePin {
                enzyme_id: e.id.clone(),
                start: hit.start,
                reverse: hit.reverse,
            };
            sites.push(Site {
                in_scope: (0..hit.length).any(|i| contains(&r.scope, (hit.start + i) % n, n)),
                retained: r.retained_sites.contains(&pin),
                pin,
                recognition: e.enzyme.recognition.clone(),
            });
        }
    }
    sites.sort_by(|a, b| a.pin.cmp(&b.pin));
    Ok(sites)
}
fn uncertain_coverage(location: &Location) -> bool {
    match location {
        Location::Range { before, after, .. } => *before || *after,
        Location::Complement(inner) => uncertain_coverage(inner),
        Location::Join(parts) | Location::Order(parts) => parts.iter().any(uncertain_coverage),
        Location::Unresolved(_) => true,
    }
}
struct Mapping {
    pin: FeaturePin,
    label: String,
    positions: Vec<(usize, bool)>,
    translation: Option<translation::CdsTranslation>,
    reason: Option<String>,
    selected: bool,
    uncertain_coverage: bool,
}
fn append(
    location: &Location,
    reverse: bool,
    out: &mut Vec<(usize, bool)>,
    exact: &mut bool,
) -> Result<(), Error> {
    match location {
        Location::Range {
            start,
            end,
            before,
            after,
        } => {
            if *before || *after {
                *exact = false;
            }
            if out
                .len()
                .checked_add(end - start)
                .is_none_or(|n| n > MAX_MAPPED)
            {
                return Err(Error::TooLarge);
            }
            if reverse {
                out.extend((*start..*end).rev().map(|p| (p, true)));
            } else {
                out.extend((*start..*end).map(|p| (p, false)));
            }
        }
        Location::Complement(inner) => append(inner, !reverse, out, exact)?,
        Location::Join(parts) | Location::Order(parts) => {
            if matches!(location, Location::Order(_)) {
                *exact = false;
            }
            if reverse {
                for p in parts.iter().rev() {
                    append(p, true, out, exact)?;
                }
            } else {
                for p in parts {
                    append(p, false, out, exact)?;
                }
            }
        }
        Location::Unresolved(_) => *exact = false,
    }
    Ok(())
}
fn mappings(source: &Molecule, r: &Request) -> Result<Vec<Mapping>, Error> {
    let mut out = vec![];
    let mut total = 0usize;
    for (index, f) in source
        .features
        .iter()
        .enumerate()
        .filter(|(_, f)| f.kind.eq_ignore_ascii_case("CDS"))
    {
        if out.len() >= MAX_FEATURES {
            return Err(Error::TooLarge);
        }
        let mut positions = vec![];
        let mut exact = true;
        append(&f.location, false, &mut positions, &mut exact)?;
        total = total
            .checked_add(positions.len())
            .filter(|n| *n <= MAX_MAPPED)
            .ok_or(Error::TooLarge)?;
        if positions
            .iter()
            .map(|(p, _)| p)
            .collect::<BTreeSet<_>>()
            .len()
            != positions.len()
        {
            exact = false;
        }
        let selected = r.coding_features.iter().any(|p| p.index == index);
        let (translation, reason) = if !exact {
            (None,Some("Fuzzy, ordered, unresolved or repeated-position CDS is locked; no exact codon mapping is asserted.".into()))
        } else {
            match translation::translate_cds_validated(source, f, &r.source_revision) {
                Ok(t) => (Some(t), None),
                Err(e) => (None, Some(e.to_string())),
            }
        };
        out.push(Mapping {
            pin: feature_pin(source, index)?,
            label: f.label().chars().take(256).collect(),
            positions,
            translation,
            reason,
            selected,
            uncertain_coverage: uncertain_coverage(&f.location),
        });
    }
    Ok(out)
}
fn feature_assessments(maps: &[Mapping]) -> Vec<CdsAssessment> {
    maps.iter()
        .map(|m| CdsAssessment {
            pin: m.pin.clone(),
            label: m.label.clone(),
            selected: m.selected,
            supported: m.translation.is_some(),
            reason: m.reason.clone(),
            codon_start: m.translation.as_ref().map(|t| t.codon_start),
            genetic_code: m
                .translation
                .as_ref()
                .map(|t| t.translation.options.genetic_code),
            has_terminal_stop: m.translation.as_ref().map(|t| t.has_terminal_stop),
            protein: m.translation.as_ref().map(|t| t.protein.clone()),
        })
        .collect()
}
fn warnings() -> Vec<String> {
    vec!["Reference substitution intent only. Site removal and preserved modeled translation do not establish chemical preparation, mutation efficiency, gene function or experimental success.".into(),"All retained CDS annotations constrain affected edits. Unsupported/overlapping coding positions remain locked during proposals; explicit patches must independently preserve every affected supported CDS. Canonical annotation editing preserves supported unaffected claims; affected coordinate/translation/phase declarations remain bounded source evidence requiring review.".into()]
}
pub fn inspect(source: &Molecule, r: &Request) -> Result<Inspection, Error> {
    let context_revision = inspect_revision(source, r)?;
    let sites = scan(source, r)?;
    for pin in &r.retained_sites {
        if !sites.iter().any(|s| s.pin == *pin) {
            return Err(bad(
                "retained recognition pin is not an actual original occurrence",
            ));
        }
    }
    let maps = mappings(source, r)?;
    let blocking = maps
        .iter()
        .filter(|m| m.selected && m.translation.is_none())
        .map(|m| {
            format!(
                "{}: {}",
                m.label,
                m.reason.as_deref().unwrap_or("Unsupported CDS")
            )
        })
        .collect();
    let out = Inspection {
        context_revision,
        request: r.clone(),
        sites,
        features: feature_assessments(&maps),
        blocking,
        warnings: warnings(),
    };
    out.validate_binding(&out.context_revision)?;
    Ok(out)
}
fn affected(map: &Mapping, patches: &[Patch]) -> bool {
    if map.uncertain_coverage || map.positions.is_empty() && map.translation.is_none() {
        return true;
    }
    map.positions
        .iter()
        .any(|(p, _)| patches.iter().any(|edit| edit.position == *p))
}
fn patched(source: &Molecule, patches: &[Patch]) -> Result<Molecule, Error> {
    let mut target = source.clone();
    let mut bases = source.sequence.as_str().as_bytes().to_vec();
    for p in patches {
        bases[p.position] = p.replacement.as_bytes()[0];
    }
    target.sequence = Sequence::dna(std::str::from_utf8(&bases).map_err(|_| Error::TooLarge)?)?;
    Ok(target)
}
fn freeze_target(
    source: &Molecule,
    r: &Request,
    s: &Selection,
    letters: &Molecule,
    context: &str,
) -> Result<Molecule, Error> {
    let first = &s.patches[0];
    let interval = mutation_intent::Request {
        name: "Domestication intended reference".into(),
        edit_start: first.position,
        edit_end: first.position + 1,
        inserted: first.replacement.clone(),
    };
    // Per-base substitutions retain original coordinates and bounded citation
    // recipes. Never map a distant spanning interval as one replacement.
    let mut target = mutation_intent::build(source, &interval)?.target;
    for p in s.patches.iter().skip(1) {
        target = target.replace(p.position, p.position + 1, &p.replacement)?;
    }
    // Keep the canonical edit result. Replacing it with source annotations can
    // reactivate a codon/anticodon claim whose biological calls were changed.
    target
        .metadata
        .retain(|k, _| !k.starts_with("mutation_intent:"));
    for (key, value) in [
        ("domestication:engine", ENGINE_VERSION.to_owned()),
        ("domestication:source_revision", r.source_revision.clone()),
        ("domestication:context_revision", context.into()),
        ("domestication:candidate_id", s.candidate_id.clone()),
        (
            "domestication:patches",
            serde_json::to_string(&s.patches).map_err(|_| Error::TooLarge)?,
        ),
    ] {
        target.metadata.insert(key.into(), value);
    }
    if target.sequence != letters.sequence {
        return Err(bad("frozen source intent letters changed"));
    }
    target.validate()?;
    Ok(target)
}

pub fn check(source: &Molecule, r: &Request, s: &Selection) -> Result<CheckedIntent, Error> {
    let context_revision = check_revision(source, r, s)?;
    let before = scan(source, r)?;
    let n = source.sequence.len();
    for pin in &r.retained_sites {
        if !before.iter().any(|v| v.pin == *pin) {
            return Err(bad("retained recognition pin changed"));
        }
    }
    for p in &s.patches {
        if &source.sequence.as_str()[p.position..p.position + 1] != p.expected
            || contains(&r.protected, p.position, n)
        {
            return Err(bad("patch original base/protected position changed"));
        }
        for site in before.iter().filter(|v| v.retained || !v.in_scope) {
            if (0..site.recognition.len()).any(|i| (site.pin.start + i) % n == p.position) {
                return Err(bad(
                    "patch touches a retained/unscoped recognition occurrence",
                ));
            }
        }
    }
    let maps = mappings(source, r)?;
    let mapped = maps.iter().map(|m| m.positions.len() as u64).sum::<u64>();
    if base_work(source, r)?
        .checked_add(mapped.checked_mul(32).ok_or(Error::TooLarge)?)
        .is_none_or(|w| w > r.limits.max_work)
    {
        return Err(bad(
            "complete CDS/motif checking exceeds declared aggregate work",
        ));
    }
    let target_letters = patched(source, &s.patches)?;
    let after = scan(&target_letters, r)?;
    if after.iter().any(|site| site.in_scope && !site.retained) {
        return Err(bad("unwanted scoped recognition sites remain"));
    }
    if after.iter().any(|a| !before.iter().any(|b| b.pin == a.pin)) {
        return Err(bad(
            "substitutions create a new selected-enzyme site, including reverse strand/origin",
        ));
    }
    if before
        .iter()
        .filter(|s| s.retained || !s.in_scope)
        .any(|b| !after.iter().any(|a| a.pin == b.pin))
    {
        return Err(bad("wanted original recognition site was removed"));
    }
    // Additional explicitly declared motifs constrain full physical DNA, not
    // only a spliced CDS. The full scan is separate from codon search receipts.
    for motif in &r.options.forbidden_motifs {
        if !analysis::find_matches(&target_letters, motif, true)?.is_empty() {
            return Err(bad("explicit forbidden motif remains in full intended DNA"));
        }
    }
    let mut features = vec![];
    let checked_letters_revision = target_letters.revision();
    for map in maps.iter().filter(|m| affected(m, &s.patches)) {
        let original = map
            .translation
            .as_ref()
            .ok_or_else(|| bad("affected CDS has unsupported/uncertain translation semantics"))?;
        let feature = &source.features[map.pin.index];
        let translated = translation::translate_cds_validated(
            &target_letters,
            feature,
            &checked_letters_revision,
        )?;
        if original.protein != translated.protein
            || original.has_terminal_stop != translated.has_terminal_stop
            || original.translation.initiator_methionine_applied
                != translated.translation.initiator_methionine_applied
            || original.translation.trailing_bases != translated.translation.trailing_bases
            || original.translation.unresolved_codons > 0
            || translated.translation.unresolved_codons > 0
        {
            return Err(bad(
                "substitutions alter an affected CDS protein, initiator, terminal stop or partial-codon semantics",
            ));
        }
        // Unresolved leading/trailing partial codons are protected literally.
        let offset = usize::from(original.codon_start - 1);
        let complete = (map.positions.len().saturating_sub(offset)) / 3 * 3;
        if map.positions.iter().enumerate().any(|(i, (p, _))| {
            (i < offset || i >= offset + complete) && s.patches.iter().any(|v| v.position == *p)
        }) {
            return Err(bad("partial coding bases are locked"));
        }
        features.push(FeatureCheck {
            pin: map.pin.clone(),
            label: map.label.clone(),
            before_protein: original.protein.clone(),
            after_protein: translated.protein,
            codon_start: original.codon_start,
            genetic_code: original.translation.options.genetic_code,
            has_terminal_stop: original.has_terminal_stop,
        });
    }
    let target = freeze_target(source, r, s, &target_letters, &context_revision)?;
    let report = CheckReport {
        context_revision: context_revision.clone(),
        request: r.clone(),
        selection: s.clone(),
        source_revision: r.source_revision.clone(),
        target_revision: target.revision(),
        candidate_id: s.candidate_id.clone(),
        sites_before: before,
        sites_after: after,
        features,
        warnings: warnings(),
    };
    report.validate_binding(&context_revision)?;
    Ok(CheckedIntent {
        source: source.clone(),
        target,
        report,
    })
}
fn coding_text(source: &Molecule, positions: &[(usize, bool)]) -> String {
    positions
        .iter()
        .map(|(p, reverse)| {
            let b = source.sequence.as_str().as_bytes()[*p];
            if *reverse {
                collagene_model::complement(b, Alphabet::Dna) as char
            } else {
                b as char
            }
        })
        .collect()
}
fn runs(flags: &[bool]) -> Vec<codon_design::CodingRange> {
    let mut out = vec![];
    let mut i = 0;
    while i < flags.len() {
        if !flags[i] {
            i += 1;
            continue;
        }
        let start = i;
        while i < flags.len() && flags[i] {
            i += 1;
        }
        out.push(codon_design::CodingRange { start, end: i });
    }
    out
}
fn flank(source: &Molecule, position: usize, reverse: bool, left: bool, length: usize) -> String {
    let n = source.sequence.len();
    let mut out = vec![];
    for step in 1..=length {
        let delta = if left { -(step as i64) } else { step as i64 };
        let p = position as i64 + if reverse { -delta } else { delta };
        if source.topology == Topology::Linear && (p < 0 || p >= n as i64) {
            break;
        }
        let b = source.sequence.as_str().as_bytes()[p.rem_euclid(n as i64) as usize];
        out.push(if reverse {
            collagene_model::complement(b, Alphabet::Dna)
        } else {
            b
        });
    }
    if left {
        out.reverse();
    }
    String::from_utf8(out).unwrap_or_default()
}
pub fn propose(source: &Molecule, r: &Request) -> Result<Proposal, Error> {
    let context_revision = propose_revision(source, r)?;
    let inspection = inspect(source, r)?;
    let unwanted = inspection
        .sites
        .iter()
        .filter(|s| s.in_scope && !s.retained)
        .collect::<Vec<_>>();
    if unwanted.is_empty() {
        return Ok(Proposal {
            context_revision,
            inspection,
            candidate: None,
            check: None,
            status: ProposalStatus::NoChangeRequired,
            search_complete: true,
            searches: vec![],
            warnings: warnings(),
        });
    }
    let maps = mappings(source, r)?;
    let n = source.sequence.len();
    let mut frequency = BTreeMap::<usize, usize>::new();
    for m in &maps {
        for (p, _) in &m.positions {
            *frequency.entry(*p).or_default() += 1;
        }
    }
    let wanted_positions = positions_of(unwanted.iter().map(|s| (**s).clone()), n);
    let wanted_kept = positions_of(
        inspection
            .sites
            .iter()
            .filter(|s| s.retained || !s.in_scope)
            .cloned(),
        n,
    );
    let unknown_global = maps
        .iter()
        .any(|m| m.uncertain_coverage || m.translation.is_none() && m.positions.is_empty());
    let mut current = source.clone();
    let mut searches = vec![];
    let mut notes = warnings();
    // Reserve inspection and final checking, plus both constrained-search passes.
    // A failed/dropped independent checker must not return its budget to search.
    let mapped = maps.iter().map(|m| m.positions.len() as u64).sum::<u64>();
    let mut used = base_work(source, r)?
        .checked_mul(2)
        .and_then(|w| mapped.checked_mul(64).and_then(|m| w.checked_add(m)))
        .filter(|w| *w <= r.limits.max_work)
        .ok_or_else(|| {
            bad("inspection plus final independent checking exceeds aggregate work budget")
        })?;
    for m in maps.iter().filter(|m| m.selected) {
        let Some(t) = &m.translation else {
            notes.push(format!(
                "{} remains locked: {}",
                m.label,
                m.reason.as_deref().unwrap_or("Unsupported coding model")
            ));
            continue;
        };
        let offset = usize::from(t.codon_start - 1);
        let length = m.positions.len().saturating_sub(offset) / 3 * 3;
        if length == 0 || length > codon_design::MAX_RESIDUES * 3 {
            notes.push(format!(
                "{} has no bounded complete coding domain.",
                m.label
            ));
            continue;
        }
        let positions = &m.positions[offset..offset + length];
        let coding = coding_text(&current, positions);
        let original = Molecule::dna("source-bound CDS search", &coding, Topology::Linear)?;
        let mut locked = vec![false; length];
        let scoped = positions
            .iter()
            .map(|(p, _)| wanted_positions.contains(p))
            .collect::<Vec<_>>();
        for (i, (p, _)) in positions.iter().enumerate() {
            let codon_start = i / 3 * 3;
            let editable_codon = scoped[codon_start..codon_start + 3].iter().any(|v| *v);
            locked[i] = unknown_global
                || !editable_codon
                || contains(&r.protected, *p, n)
                || frequency.get(p).is_some_and(|v| *v > 1)
                || wanted_kept.contains(p);
            // Retain literal initiator semantics rather than let an internal-codon
            // synonym table reinterpret alternative recognized starts.
            if i < 3 {
                locked[i] = true;
            }
        }
        let mut options = r.options.clone();
        options.genetic_code = t.translation.options.genetic_code;
        let (first, rev) = positions[0];
        let (last, last_reverse) = positions[length - 1];
        options.left_context = flank(source, first, rev, true, 64);
        options.right_context = flank(source, last, last_reverse, false, 64);
        options.max_work = options
            .max_work
            .min(r.limits.max_work.saturating_sub(used) / 2);
        if options.max_work == 0 {
            notes.push("Aggregate source/search work budget exhausted; remaining CDS alternatives are unexamined.".into());
            break;
        }
        let policy = codon_design::ConstrainedPolicy {
            protected: runs(&locked),
            scoped_forbidden: r
                .selection
                .enzymes
                .iter()
                .map(|e| codon_design::ScopedMotif {
                    motif: e.enzyme.recognition.clone(),
                    ranges: runs(&scoped),
                })
                .filter(|s| !s.ranges.is_empty())
                .collect(),
            ..Default::default()
        };
        let reserved_work = options.max_work.checked_mul(2).ok_or(Error::TooLarge)?;
        used = used
            .checked_add(reserved_work)
            .filter(|v| *v <= r.limits.max_work)
            .ok_or(Error::TooLarge)?;
        let result = match codon_design::design_constrained(&original, &options, &policy) {
            Ok(v) => v,
            Err(e) => {
                notes.push(format!("{} search unavailable: {e}. Its full work reservation remains charged because actual expenditure is unknown.", m.label));
                searches.push(FeatureSearch {
                    feature_index: m.pin.index,
                    status: codon_design::DesignStatus::SearchIncomplete,
                    search_complete: false,
                    work_used: reserved_work,
                    reserved_work,
                    minimum_changes_proven: false,
                });
                continue;
            }
        };
        let actual_work = result.total_work_used;
        if actual_work > reserved_work {
            return Err(bad("constrained search exceeded its two-pass reservation"));
        }
        used = used
            .checked_sub(reserved_work)
            .and_then(|w| w.checked_add(actual_work))
            .filter(|w| *w <= r.limits.max_work)
            .ok_or(Error::TooLarge)?;
        searches.push(FeatureSearch {
            feature_index: m.pin.index,
            status: result.design.status,
            search_complete: result.design.search_complete,
            work_used: actual_work,
            reserved_work,
            minimum_changes_proven: result.minimum_changes_proven,
        });
        if let Some(candidate) = result.design.candidate {
            let mut bytes = current.sequence.as_str().as_bytes().to_vec();
            for ((p, reverse), b) in positions.iter().zip(candidate.sequence.as_str().bytes()) {
                bytes[*p] = if *reverse {
                    collagene_model::complement(b, Alphabet::Dna)
                } else {
                    b
                };
            }
            current.sequence =
                Sequence::dna(std::str::from_utf8(&bytes).map_err(|_| Error::TooLarge)?)?;
        }
    }
    let patches = source
        .sequence
        .as_str()
        .bytes()
        .zip(current.sequence.as_str().bytes())
        .enumerate()
        .filter_map(|(position, (a, b))| {
            (a != b).then(|| Patch {
                position,
                expected: (a as char).to_string(),
                replacement: (b as char).to_string(),
            })
        })
        .collect::<Vec<_>>();
    let mut candidate = None;
    let mut check_report = None;
    if !patches.is_empty() && patches.len() <= r.limits.max_patches {
        let selection = Selection {
            context_revision: context_revision.clone(),
            candidate_id: candidate_id(r, &patches)?,
            patches,
        };
        match check(source, r, &selection) {
            Ok(proof) => {
                check_report = Some(proof.report);
                candidate = Some(selection);
            }
            Err(e) => notes.push(format!(
                "Local coding search was not accepted by the independent full-source checker: {e}"
            )),
        }
    }
    if notes.len() > 16 {
        notes.truncate(15);
        notes.push("Further feature search diagnostics omitted from this bounded summary; inspect every source CDS assessment.".into());
    }
    let status = if candidate.is_some() {
        ProposalStatus::CandidateFound
    } else if searches.is_empty() {
        ProposalStatus::Blocked
    } else {
        ProposalStatus::SearchIncomplete
    };
    notes.push("Independent per-feature search conservatively locks overlapping/unselected CDS and original initiators. Joint overlapping-frame optimization and alternatives rejected by full-source rescanning are not exhausted; no global minimal-edit or impossibility claim is made. Compound-gap motif contexts may require explicit literal patches checked against every affected CDS.".into());
    Ok(Proposal {
        context_revision,
        inspection,
        candidate,
        check: check_report,
        status,
        search_complete: false,
        searches,
        warnings: notes,
    })
}
fn base_work(source: &Molecule, r: &Request) -> Result<u64, Error> {
    let weights = r
        .selection
        .enzymes
        .iter()
        .map(|e| e.enzyme.recognition.len() as u64 * 4)
        .chain(
            r.options
                .forbidden_motifs
                .iter()
                .map(|s| s.len() as u64 * 2),
        )
        .try_fold(0u64, |a, b| a.checked_add(b).ok_or(Error::TooLarge))?;
    (source.sequence.len() as u64)
        .checked_mul(weights)
        .filter(|w| *w <= r.limits.max_work)
        .ok_or(Error::TooLarge)
}
fn receipt_sites(sites: &[Site], r: &Request) -> Result<(), Error> {
    if sites.len() > MAX_SITES {
        return Err(Error::TooLarge);
    }
    let mut seen = BTreeSet::new();
    for site in sites {
        let enzyme = r
            .selection
            .enzymes
            .iter()
            .find(|e| e.id == site.pin.enzyme_id)
            .ok_or_else(|| bad("site enzyme not selected"))?;
        if site.recognition != enzyme.enzyme.recognition
            || site.retained != r.retained_sites.contains(&site.pin)
            || !seen.insert(&site.pin)
        {
            return Err(bad("site identity changed"));
        }
    }
    Ok(())
}
fn notes(values: &[String], max: usize) -> Result<(), Error> {
    if values.len() > max || values.iter().any(|s| s.len() > 4096) {
        return Err(Error::TooLarge);
    }
    Ok(())
}
impl Inspection {
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        if self.context_revision != expected
            || expected != context("inspect", &self.request)?
            || self.features.len() > MAX_FEATURES
        {
            return Err(bad("inspection binding changed"));
        }
        receipt_sites(&self.sites, &self.request)?;
        notes(&self.warnings, 16)?;
        notes(&self.blocking, MAX_FEATURES)?;
        let mut seen = BTreeSet::new();
        for f in &self.features {
            if !hex(&f.pin.revision)
                || !seen.insert(f.pin.index)
                || f.label.len() > 1024
                || f.protein.as_ref().is_some_and(|p| p.len() > MAX_MAPPED)
                || f.supported != f.protein.is_some()
                || f.reason.as_ref().is_some_and(|s| s.len() > 4096)
            {
                return Err(bad("CDS assessment shape changed"));
            }
        }
        crate::encoding::bounded_json::hash(self)?;
        Ok(())
    }
}
impl Proposal {
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        if self.context_revision != expected
            || expected != context("propose", &self.inspection.request)?
            || self.searches.len() > MAX_FEATURES
        {
            return Err(bad("proposal binding changed"));
        }
        self.inspection
            .validate_binding(&context("inspect", &self.inspection.request)?)?;
        notes(&self.warnings, 17)?;
        let mut charged = 0u64;
        let mut searched = BTreeSet::new();
        for search in &self.searches {
            charged = charged
                .checked_add(search.work_used)
                .filter(|w| *w <= self.inspection.request.limits.max_work)
                .ok_or(Error::TooLarge)?;
            if search.reserved_work > self.inspection.request.limits.max_work
                || search.work_used > search.reserved_work
                || !searched.insert(search.feature_index)
                || !self
                    .inspection
                    .request
                    .coding_features
                    .iter()
                    .any(|p| p.index == search.feature_index)
            {
                return Err(bad("search work/feature declaration changed"));
            }
        }
        match (&self.candidate, &self.check) {
            (Some(s), Some(c)) => {
                if !matches!(self.status, ProposalStatus::CandidateFound)
                    || s.context_revision != expected
                    || s.candidate_id != candidate_id(&self.inspection.request, &s.patches)?
                {
                    return Err(bad("candidate declaration changed"));
                }
                let token = crate::encoding::bounded_json::hash(&(
                    ENGINE_VERSION,
                    "check",
                    &self.inspection.request,
                    s,
                ))?;
                c.validate_binding(&token)?;
                if crate::encoding::bounded_json::hash(&c.selection)?
                    != crate::encoding::bounded_json::hash(s)?
                {
                    return Err(bad("candidate check belongs to another selection"));
                }
            }
            (None, None) => {
                if matches!(self.status, ProposalStatus::CandidateFound) {
                    return Err(bad("candidate status has no checked declaration"));
                }
            }
            _ => return Err(bad("incomplete candidate/check pair")),
        }
        crate::encoding::bounded_json::hash(self)?;
        Ok(())
    }
}
impl CheckReport {
    /// Read-only shape/binding checks. Fresh authority lives only in CheckedIntent.
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        request_shape(&self.request)?;
        patch_shape(&self.request, &self.selection.patches)?;
        if self.context_revision != expected
            || !hex(&self.target_revision)
            || self.source_revision != self.request.source_revision
            || self.candidate_id != self.selection.candidate_id
            || self.selection.context_revision != context("propose", &self.request)?
            || self.selection.candidate_id != candidate_id(&self.request, &self.selection.patches)?
            || expected
                != crate::encoding::bounded_json::hash(&(
                    ENGINE_VERSION,
                    "check",
                    &self.request,
                    &self.selection,
                ))?
            || self.features.len() > MAX_FEATURES
        {
            return Err(bad("check receipt binding changed"));
        }
        receipt_sites(&self.sites_before, &self.request)?;
        receipt_sites(&self.sites_after, &self.request)?;
        notes(&self.warnings, 16)?;
        if self.sites_after.iter().any(|s| s.in_scope && !s.retained)
            || self
                .sites_after
                .iter()
                .any(|s| !self.sites_before.iter().any(|b| b.pin == s.pin))
        {
            return Err(bad("checked site assertions changed"));
        }
        if self
            .sites_before
            .iter()
            .filter(|site| site.retained || !site.in_scope)
            .any(|before| !self.sites_after.iter().any(|after| after.pin == before.pin))
        {
            return Err(bad("retained site assertion changed"));
        }
        let mut seen = BTreeSet::new();
        for f in &self.features {
            if !hex(&f.pin.revision)
                || !seen.insert(f.pin.index)
                || f.label.len() > 1024
                || f.before_protein.len() > MAX_MAPPED
                || f.after_protein != f.before_protein
                || !(1..=3).contains(&f.codon_start)
            {
                return Err(bad("checked CDS assertions changed"));
            }
        }
        crate::encoding::bounded_json::hash(self)?;
        Ok(())
    }
}

pub mod preparation;
pub(crate) mod rounds;
