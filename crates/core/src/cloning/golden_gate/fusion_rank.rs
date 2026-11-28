//! Conditional published-profile ranking of retained source-bound fusion sets.
//! Ranking never changes a recipe, chooses a reaction winner or grants physical authority.
use crate::{
    cloning::{golden_gate::fusion_design as fusion, ligation::profiles},
    planning::plan::PlanInput,
};
use collagene_model::{Error, Molecule, Sequence};
use serde::{Deserialize, Serialize};
use std::{
    cmp::Ordering,
    collections::{BTreeMap, BTreeSet},
    io::{self, Write},
};
pub const ENGINE_VERSION: &str = "golden-gate-fusion-rank-1";
pub const MAX_REPORT_BYTES: usize = 32 << 20;
macro_rules! dto {($($t:tt)*)=>{#[derive(Clone,Debug,Serialize,Deserialize)]#[cfg_attr(feature="json-schema",derive(schemars::JsonSchema))]#[serde(deny_unknown_fields)]$($t)*};}
dto! {pub struct Request {
    pub design:fusion::Request,pub profile:profiles::ProfileRef,
    /// Explicit conditional assay assumption, not inferred from generic cleavage
    /// geometry/catalog names or a claim about actual supplied material.
    pub enzyme_formulation:String,
}}
dto! {#[derive(Copy,PartialEq,Eq)] pub enum ScoreClass {Positive,Zero,Unknown}}
dto! {pub struct PairEstimate {pub canonical_end:String,pub reverse_complement:String,pub estimate:profiles::Estimate}}
dto! {pub struct Score {pub class:ScoreClass,pub overall_fidelity:Option<f64>,pub log_overall_fidelity:Option<f64>,pub pairs:Vec<PairEstimate>}}
dto! {pub struct RankedCandidate {pub rank:usize,pub candidate_id:String,pub score:Score}}
dto! {pub struct Report {
    pub context_revision:String,pub request:Request,pub metadata:profiles::Metadata,
    /// Original source/intent context and candidate identities remain available
    /// for explicit materialization; profile ranking supplies no replacement DNA.
    pub search:fusion::Report,pub ranked:Vec<RankedCandidate>,pub limits:Vec<String>,pub warnings:Vec<String>,
}}
fn bad(message: &str) -> Error {
    Error::Unsupported(format!("Golden Gate fusion ranking: {message}"))
}
fn digest(value: &str) -> bool {
    value.len() == 64
        && value
            .bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}
fn rc(value: &str) -> Result<String, Error> {
    Ok(Sequence::dna(value)?.reverse_complement()?.as_str().into())
}
fn profile(request: &Request) -> Result<profiles::Summary, Error> {
    if request.enzyme_formulation.is_empty()
        || request.enzyme_formulation.len() > 128
        || request.enzyme_formulation.chars().any(char::is_control)
        || request.profile.profile_id.is_empty()
        || request.profile.profile_id.len() > 128
        || !digest(&request.profile.profile_revision)
    {
        return Err(bad(
            "choose an explicit bounded enzyme-formulation assumption",
        ));
    }
    let metadata = profiles::summaries()?
        .into_iter()
        .find(|m| {
            m.profile.profile_id == request.profile.profile_id
                && m.profile.profile_revision == request.profile.profile_revision
        })
        .ok_or_else(|| bad("published profile identity/revision is unavailable"))?;
    if metadata.profile.profile_id != request.profile.profile_id
        || metadata.profile.profile_revision != request.profile.profile_revision
    {
        return Err(bad("published profile identity or version changed"));
    }
    let geometry = request.design.base.enzyme.profile();
    let length =
        usize::try_from(geometry.bottom_cut - geometry.top_cut).map_err(|_| Error::TooLarge)?;
    if metadata.enzyme_formulation != request.enzyme_formulation {
        return Err(bad(
            "the explicit formulation assumption differs from the selected published profile; formulations are not silently substituted",
        ));
    }
    if metadata.geometry_family != request.design.base.enzyme || metadata.overhang_length != length
    {
        return Err(bad(
            "published profile recognition/cut geometry family or overhang length differs from the design",
        ));
    }
    Ok(metadata)
}
pub fn preflight(inputs: &[(&str, &Molecule)], request: &Request) -> Result<(), Error> {
    profile(request)?;
    fusion::preflight(inputs, &request.design)?;
    crate::encoding::bounded_json::hash(request)?;
    Ok(())
}
fn context(
    source_context: &str,
    request: &Request,
    metadata: &profiles::Summary,
) -> Result<String, Error> {
    crate::encoding::bounded_json::hash(&(ENGINE_VERSION, source_context, request, metadata))
}
pub fn inspect_revision(inputs: &[PlanInput], request: &Request) -> Result<String, Error> {
    let metadata = profile(request)?;
    context(
        &fusion::inspect_revision(inputs, &request.design)?,
        request,
        &metadata,
    )
}
fn pairs(ends: &[String], k: usize) -> Result<BTreeMap<String, String>, Error> {
    if ends.len() > crate::cloning::ligation::fidelity::MAX_ENDS {
        return Err(bad("bounded intended junction set required"));
    }
    let mut result = BTreeMap::new();
    for end in ends {
        if end.len() != k || !end.bytes().all(|b| b"ACGT".contains(&b)) {
            return Err(bad(
                "fusion set differs from the profile end length/alphabet",
            ));
        }
        let reverse = rc(end)?;
        let canonical = end.min(&reverse).clone();
        result.insert(canonical.clone(), rc(&canonical)?);
    }
    Ok(result)
}
fn summary(report: profiles::Report) -> Result<Score, Error> {
    let score = Score {
        class: if report.log_overall_fidelity.is_some() {
            ScoreClass::Positive
        } else if report.overall_fidelity == Some(0.0) {
            ScoreClass::Zero
        } else {
            ScoreClass::Unknown
        },
        overall_fidelity: report.overall_fidelity,
        log_overall_fidelity: report.log_overall_fidelity,
        pairs: report
            .pairs
            .into_iter()
            .map(|p| PairEstimate {
                canonical_end: p.canonical_end,
                reverse_complement: p.reverse_complement,
                estimate: p.estimate,
            })
            .collect(),
    };
    score.validate(&report.requested_ends, report.metadata.overhang_length)?;
    Ok(score)
}
impl Score {
    /// Compact mathematical consistency only. Original observation cells are not
    /// duplicated per candidate; replay native ranking to establish current scores.
    pub fn validate(&self, ends: &[String], k: usize) -> Result<(), Error> {
        let expected = pairs(ends, k)?;
        if self.pairs.len() != expected.len() {
            return Err(bad("pair score count changed"));
        }
        let mut product = 1.0f64;
        let mut logarithm = 0.0f64;
        let mut known = !ends.is_empty();
        let mut zero = false;
        for (pair, (canonical, reverse)) in self.pairs.iter().zip(expected) {
            if pair.canonical_end != canonical || pair.reverse_complement != reverse {
                return Err(bad("scored pair identity changed"));
            }
            match &pair.estimate {
                profiles::Estimate::Known {
                    correct_weight,
                    total_weight,
                    fidelity,
                } => {
                    if *total_weight == 0
                        || correct_weight > total_weight
                        || !fidelity.is_finite()
                        || *fidelity != *correct_weight as f64 / *total_weight as f64
                        || canonical == reverse
                    {
                        return Err(bad("invalid conditional pair weights/ratio"));
                    }
                    product *= fidelity;
                    if *correct_weight == 0 {
                        zero = true;
                    } else {
                        logarithm += fidelity.ln();
                    }
                }
                profiles::Estimate::Unknown { reason } => {
                    if (canonical == reverse)
                        != matches!(
                            reason,
                            crate::cloning::ligation::fidelity::UnknownReason::PalindromicEnd
                        )
                    {
                        return Err(bad("palindromic pair reason changed"));
                    }
                    known = false;
                }
            }
        }
        let overall = known.then_some(product);
        let log = (known && !zero).then_some(logarithm);
        let class = if log.is_some() {
            ScoreClass::Positive
        } else if overall == Some(0.0) {
            ScoreClass::Zero
        } else {
            ScoreClass::Unknown
        };
        if self.overall_fidelity != overall
            || self.log_overall_fidelity != log
            || self.class != class
            || log.is_some_and(|v| !v.is_finite())
        {
            return Err(bad("conditional overall/log score changed"));
        }
        Ok(())
    }
}
fn compare(a: &RankedCandidate, b: &RankedCandidate) -> Ordering {
    let category = |score: &Score| match score.class {
        ScoreClass::Positive => 0,
        ScoreClass::Zero => 1,
        ScoreClass::Unknown => 2,
    };
    category(&a.score)
        .cmp(&category(&b.score))
        .then_with(|| {
            if a.score.class == ScoreClass::Positive && b.score.class == ScoreClass::Positive {
                b.score
                    .log_overall_fidelity
                    .expect("validated positive score")
                    .total_cmp(
                        &a.score
                            .log_overall_fidelity
                            .expect("validated positive score"),
                    )
            } else {
                Ordering::Equal
            }
        })
        .then(a.candidate_id.cmp(&b.candidate_id))
}
pub fn inspect(inputs: &[PlanInput], request: &Request) -> Result<Report, Error> {
    let metadata = profile(request)?;
    // Exactly one bounded source search and one validated matrix index. Neither
    // is repeated for every retained set or inferred from an imported report.
    let search = fusion::inspect(inputs, &request.design)?;
    // XLSX decoding and matrix validation only occur in this native/worker path.
    // Revision/selection/wire callbacks above/below use fixed summaries alone.
    let profile = profiles::load(&request.profile)?;
    let prepared = profile.prepare()?;
    let mut ranked = Vec::with_capacity(search.candidates.len());
    for candidate in &search.candidates {
        let scope = fusion::candidate_end_scope(&search, candidate)?;
        ranked.push(RankedCandidate {
            rank: 0,
            candidate_id: candidate.id.clone(),
            score: summary(
                prepared.score_with_competitors(
                    &scope.intended_junctions,
                    &scope.physical_competitors,
                )?,
            )?,
        });
    }
    ranked.sort_by(compare);
    for (index, candidate) in ranked.iter_mut().enumerate() {
        candidate.rank = index + 1;
    }
    let mut limits = search.limits.clone();
    if ranked.iter().any(|c| c.score.class == ScoreClass::Unknown) {
        limits.push("Some retained sets have unknown estimates (missing counts, zero totals or palindromic ends); unknown is not zero and is ranked last.".into());
    }
    let context_revision = context(&search.context_revision, request, &metadata)?;
    let report=Report {context_revision:context_revision.clone(),request:request.clone(),metadata:profile.metadata().clone(),search,ranked,limits,
        warnings:vec![
            "Conditional published-frequency ranking only. Native recipe catalog identities and cleavage geometry remain unchanged; choosing a formulation hypothesis never converts generic BsaI to BsaI-HFv2 or substitutes Esp3I/BsmBI formulations in a protocol.".into(),
            "Actual assay ligase, buffer, cycling, substrate context, abundance and independence assumptions are not established for your material. Profile conditions/data/model revisions remain explicit; these values do not predict universal efficiency or experimental success.".into(),
            "Fixed linear outer cohesive ends are literal singleton competitors, without invented complementary ends or multiplied outer-end joining events. Only intended junction pairs contribute to the product; no intended joins gives an unknown/not-applicable estimate. Repeated ends do not establish concentration weighting.".into(),
            "Only retained mathematical search candidates are ranked; missing/pruned placements and unmodeled reaction paths remain outside this result. No global optimum, unique reaction winner, target-driven choice or automatic Apply is supplied.".into(),
            "Detailed orientation/partner cells are omitted from this compact ranking. Use explicit profile scoring for the same intended candidate junctions and separately declared literal outer-end competitors to inspect those counts; physical preparation still requires separate explicit materialization and verification.".into(),
        ]};
    report.validate_binding(&context_revision)?;
    Ok(report)
}
fn notes(values: &[String], max: usize) -> Result<(), Error> {
    if values.len() > max || values.iter().any(|s| s.len() > 4096) {
        return Err(Error::TooLarge);
    }
    Ok(())
}
fn metadata_matches(m: &profiles::Metadata, s: &profiles::Summary) -> bool {
    m.profile == s.profile
        && m.version == s.version
        && m.model == s.model
        && m.name == s.name
        && m.geometry_family == s.geometry_family
        && m.enzyme_formulation == s.enzyme_formulation
        && m.overhang_length == s.overhang_length
        && m.conditions == s.conditions
        && m.provenance == s.provenance
        && m.source_sha256 == s.source_sha256
        && m.cells == s.cells
}
impl Report {
    /// Read-only source/profile/candidate and compact-math binding. Imported data
    /// cannot restore a fresh ranking or authorize a physical preparation.
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        let summary = profile(&self.request)?;
        if !digest(expected)
            || self.context_revision != expected
            || !metadata_matches(&self.metadata, &summary)
            || !digest(&self.metadata.data_revision)
            || crate::encoding::bounded_json::hash(&self.request.design)?
                != crate::encoding::bounded_json::hash(&self.search.request)?
            || context(&self.search.context_revision, &self.request, &summary)? != expected
            || self.ranked.len() != self.search.candidates.len()
            || self.ranked.len() > fusion::MAX_CANDIDATES
        {
            return Err(bad("ranking source/profile/search binding changed"));
        }
        self.search
            .validate_binding(&self.search.context_revision)?;
        let mut seen = BTreeSet::new();
        for (index, ranked) in self.ranked.iter().enumerate() {
            let candidate = self
                .search
                .candidates
                .iter()
                .find(|c| c.id == ranked.candidate_id)
                .ok_or_else(|| bad("ranked candidate was not retained by original search"))?;
            if ranked.rank != index + 1 || !seen.insert(&ranked.candidate_id) {
                return Err(bad("rank or candidate identity changed"));
            }
            ranked
                .score
                .validate(&candidate.fusions, self.metadata.overhang_length)?;
        }
        if self
            .ranked
            .windows(2)
            .any(|w| compare(&w[0], &w[1]) == Ordering::Greater)
        {
            return Err(bad("ranking order changed"));
        }
        notes(&self.limits, 16)?;
        notes(&self.warnings, 16)?;
        serde_json::to_writer(&mut Size(0), self).map_err(|_| Error::TooLarge)?;
        Ok(())
    }
}
struct Size(usize);
impl Write for Size {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        self.0 = self
            .0
            .checked_add(bytes.len())
            .filter(|n| *n <= MAX_REPORT_BYTES)
            .ok_or_else(|| io::Error::other("fusion ranking report bound"))?;
        Ok(bytes.len())
    }
    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}
/// Explicit handoff only: preserve the original fusion source/intent context and
/// chosen stable ID. The caller supplies the current ranking context; no first-row
/// default and no profile/report receipt is accepted as chemical authority.
pub fn select(
    report: &Report,
    expected_context: &str,
    candidate_id: &str,
) -> Result<fusion::Selection, Error> {
    report.validate_binding(expected_context)?;
    if !report.ranked.iter().any(|c| c.candidate_id == candidate_id) {
        return Err(bad("explicit ranked candidate is not available"));
    }
    Ok(fusion::Selection {
        context_revision: report.search.context_revision.clone(),
        candidate_id: candidate_id.into(),
    })
}

/// Scope for a detailed conditional-score request, preserving the exact original
/// search candidate and actual outer-end orientation. No scientific replay.
pub fn score_scope(
    report: &Report,
    expected_context: &str,
    candidate_id: &str,
) -> Result<fusion::EndScope, Error> {
    report.validate_binding(expected_context)?;
    let candidate = report
        .search
        .candidates
        .iter()
        .find(|c| c.id == candidate_id)
        .ok_or_else(|| bad("explicit ranked candidate is unavailable"))?;
    fusion::candidate_end_scope(&report.search, candidate)
}
