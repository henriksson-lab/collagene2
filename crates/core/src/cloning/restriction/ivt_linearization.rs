//! Alternative single-enzyme IVT linearization hypotheses, never completed cuts.
//! Recognition, positional cleavage and fragment construction stay in their
//! existing owners. Only coding-coordinate projection and candidate policy live here.
use super::{catalog, conditions};
use crate::{
    chemistry::{
        DuplexSummary,
        methylation::Outcome,
        tracked_duplex::{
            self as tracked, CheckedCleavageEvidence, CheckedDigestOptions, CleavageBudget,
            GenerationAssumptions, ModificationSummary, TrackedDuplex,
        },
    },
    material::run::MaterialOrigin,
    planning::plan::{DesignPlan, Operation, PlanInput, Source},
    rna_chemistry::{self, Provenance, ivt},
    verification,
};
use collagene_model::{Error, Molecule, Topology};
use serde::{Deserialize, Serialize};
use std::io::Write;

pub const ENGINE: &str = "ivt-restriction-linearization-1";
pub const MAX_ENZYMES: usize = 64;
pub const MAX_REQUEST_BYTES: usize = 8 * 1024 * 1024;
pub const MAX_REPORT_BYTES: usize = 8 * 1024 * 1024;
pub const MAX_SITES_PER_ENZYME: usize = 64;
const PROMOTER_BASES: usize = ivt::profiles::T7_PROMOTER.len();
const STATEMENT: &str = "Alternative single-enzyme reference hypotheses under exact declared source generation or freshly replayed material, positional cleavage and documentary conditions. Eligible means only that the declared unique-cut/runoff/retention policy is satisfied. This calculation establishes no physical cleavage, experimental completion, isolation, yield, capping efficiency or automatic winner.";

#[derive(Debug, Clone, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub enum SourceInput {
    /// Uses the existing native original-input generation guard, never metadata recovery.
    DeclaredOriginal { generation: GenerationAssumptions },
    /// Existing closed native recipe is replayed; serialized tracked state is not accepted.
    Saved { origin: MaterialOrigin },
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Policy {
    #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, ivt::EndKind, 3>")]
    pub allowed_end_kinds: Vec<ivt::EndKind>,
    pub require_known_conditions: bool,
    pub require_resolved_ivt_checks: bool,
    pub require_complete_origin: bool,
}
impl Policy {
    fn validate(&self) -> Result<(), Error> {
        if self.allowed_end_kinds.is_empty() || self.allowed_end_kinds.len() > 3 {
            return Err(Error::Format(
                "explicit IVT end-kind policy required".into(),
            ));
        }
        for (i, kind) in self.allowed_end_kinds.iter().enumerate() {
            if self.allowed_end_kinds[..i].contains(kind) {
                return Err(Error::Format("duplicate IVT end kind".into()));
            }
        }
        Ok(())
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Limits {
    pub max_comparisons: usize,
    pub max_sites_per_enzyme: usize,
    pub max_report_bytes: usize,
}
impl Limits {
    fn validate(&self) -> Result<(), Error> {
        if !(1..=100_000_000).contains(&self.max_comparisons)
            || !(1..=MAX_SITES_PER_ENZYME).contains(&self.max_sites_per_enzyme)
            || !(1024..=MAX_REPORT_BYTES).contains(&self.max_report_bytes)
        {
            return Err(Error::TooLarge);
        }
        Ok(())
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Request {
    /// Base IVT declaration. linearization must be None: candidate geometry is derived here.
    pub ivt: ivt::Request,
    pub source: SourceInput,
    /// Alternatives, not an implicitly combined digestion cocktail.
    pub selection: catalog::CatalogSelection,
    #[serde(deserialize_with = "crate::encoding::required_option")]
    #[cfg_attr(feature = "json-schema", schemars(required))]
    pub conditions: Option<conditions::ReactionConditions>,
    pub policy: Policy,
    pub limits: Limits,
}
fn rna_error(error: rna_chemistry::Error) -> Error {
    match error {
        rna_chemistry::Error::Limit(_) => Error::TooLarge,
        other => Error::Format(other.to_string()),
    }
}
impl Request {
    pub fn validate(&self, source: &Molecule) -> Result<(), Error> {
        size(self, MAX_REQUEST_BYTES)?;
        self.ivt.validate_for(source).map_err(rna_error)?;
        self.policy.validate()?;
        self.limits.validate()?;
        let current = catalog::builtin_catalog()?;
        self.selection.validate_for(&current)?;
        if self.selection.enzymes.len() > MAX_ENZYMES || self.ivt.linearization.is_some() {
            return Err(Error::Format(
                "IVT picker requires bounded alternatives and no preselected cut geometry".into(),
            ));
        }
        if let Some(options) = &self.conditions {
            options.validate(&self.selection)?;
        }
        match &self.source {
            SourceInput::DeclaredOriginal { generation } => generation.validate()?,
            SourceInput::Saved { origin } => {
                let reference = match origin {
                    MaterialOrigin::Run(capsule) => {
                        crate::material::capsule::NativeOriginRef::Run(capsule)
                    }
                    MaterialOrigin::Derived(capsule) => {
                        crate::material::capsule::NativeOriginRef::Derived(capsule)
                    }
                };
                crate::material::capsule::native_origins_preflight([reference])?;
                if origin.event().revision != source.revision() {
                    return Err(Error::Format(
                        "selected IVT material source revision differs".into(),
                    ));
                }
            }
        }
        Ok(())
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum CandidateState {
    EligibleUnderDeclarations,
    Rejected,
    Unresolved,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum Severity {
    Reject,
    Unresolved,
    Advisory,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum IssueKind {
    IncompleteScan,
    UnresolvedCleavage,
    NoPermittedCut,
    MultiplePermittedCuts,
    UnavailableCutCoordinates,
    CutOffsetsOutsideSource,
    PromoterOrTranscriptNotRetained,
    RunoffMismatch,
    DisallowedEndKind,
    ConditionsNotEstablished,
    ConditionsBelowThreshold,
    IncompleteOrigin,
    IvtFailure,
    IvtUnresolved,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Issue {
    pub kind: IssueKind,
    pub severity: Severity,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct CutMapping {
    pub site_index: usize,
    pub top_bond: usize,
    pub bottom_bond: usize,
    /// Lifted offsets along coding direction. They are never normalized independently.
    pub coding_cut_offset: i64,
    pub template_cut_offset: i64,
    pub signed_template_minus_coding: i64,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Fragment {
    pub index: usize,
    pub source_top_start: usize,
    pub source_bottom_start: i64,
    pub geometry: DuplexSummary,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Candidate {
    pub enzyme_id: String,
    pub state: CandidateState,
    #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, Issue, 32>")]
    pub issues: Vec<Issue>,
    pub cleavage: CheckedCleavageEvidence,
    #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, Fragment, 65>")]
    pub fragments: Vec<Fragment>,
    #[serde(deserialize_with = "crate::encoding::required_option")]
    #[cfg_attr(feature = "json-schema", schemars(required))]
    pub conditions: Option<conditions::source::SourceReport>,
    #[serde(deserialize_with = "crate::encoding::required_option")]
    #[cfg_attr(feature = "json-schema", schemars(required))]
    pub mapping: Option<CutMapping>,
    #[serde(deserialize_with = "crate::encoding::required_option")]
    #[cfg_attr(feature = "json-schema", schemars(required))]
    pub retained_fragment: Option<usize>,
    #[serde(deserialize_with = "crate::encoding::required_option")]
    #[cfg_attr(feature = "json-schema", schemars(required))]
    pub ivt: Option<ivt::Report>,
}
fn state(issues: &[Issue]) -> CandidateState {
    if issues.iter().any(|i| i.severity == Severity::Reject) {
        CandidateState::Rejected
    } else if issues.iter().any(|i| i.severity == Severity::Unresolved) {
        CandidateState::Unresolved
    } else {
        CandidateState::EligibleUnderDeclarations
    }
}
fn issue(issues: &mut Vec<Issue>, kind: IssueKind, severity: Severity) {
    if !issues.iter().any(|i| i.kind == kind) {
        issues.push(Issue { kind, severity });
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Report {
    pub engine_version: String,
    pub context_revision: String,
    pub source_revision: String,
    pub source_topology: Topology,
    pub request: Request,
    pub source_modifications: ModificationSummary,
    pub source_history_incomplete: bool,
    pub comparisons: usize,
    /// All requested alternatives are represented, including budget-limited searches.
    #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, Candidate, 64>")]
    pub candidates: Vec<Candidate>,
    pub statement: String,
    pub report_revision: String,
}
pub fn revision(source: &Molecule, request: &Request) -> Result<String, Error> {
    request.validate(source)?;
    crate::encoding::bounded_json::hash(&(ENGINE, source.revision(), request))
}

/// Fresh source acquisition; no serialized TrackedDuplex enters the request.
fn material(source: &Molecule, input: &SourceInput) -> Result<(TrackedDuplex, bool), Error> {
    let (tracked, incomplete) = match input {
        SourceInput::Saved { origin } => {
            let fresh = origin.replay(source)?;
            (fresh.tracked().clone(), fresh.incomplete())
        }
        SourceInput::DeclaredOriginal { generation } => {
            let plan = DesignPlan {
                standard: None,
                name: "Declared original IVT template generation".into(),
                target_revision: source.revision(),
                inputs: vec![PlanInput {
                    id: "template".into(),
                    revision: source.revision(),
                    molecule: source.clone(),
                }],
                steps: vec![Operation::DeclareSourceModifications {
                    template: Source::Input("template".into()),
                    assumptions: generation.clone(),
                }],
            };
            let (report, mut states) = verification::engine::verify_with_tracked(&plan, source);
            let tracked = states.get_mut(0).and_then(Option::take).ok_or_else(|| {
                Error::Unsupported("original IVT source-generation declaration was refused by native source guards; reconstruct the actual originating material".into())
            })?;
            (
                tracked,
                report.status != verification::evidence::Status::ConsistentUnderAssumptions,
            )
        }
    };
    tracked.validate()?;
    if tracked.duplex.top != source.sequence
        || tracked.duplex.topology != source.topology
        || tracked.duplex.bottom_start != 0
        || tracked.duplex.bottom.len() != source.sequence.len()
        || !tracked.duplex.nicks.is_empty()
    {
        return Err(Error::Unsupported("IVT cut picker currently requires an intact fully paired source duplex without nicks or preexisting staggered ends".into()));
    }
    Ok((tracked, incomplete))
}

/// Project already established bonds; no recognition or cleavage is recalculated.
fn project(
    source_length: usize,
    topology: Topology,
    transcript: &ivt::Transcript,
    selected: &catalog::SelectedEnzyme,
    site: &tracked::CheckedSite,
    site_index: usize,
) -> Option<CutMapping> {
    let top = site.top_cut?;
    let bottom = site.bottom_cut?;
    // Documentary reports are untrusted too. Bound before signed coordinate
    // arithmetic, even though freshly inspected bonds already satisfy this.
    if source_length == 0
        || source_length > ivt::MAX_TEMPLATE_BASES
        || transcript.start >= source_length
        || top > source_length
        || bottom > source_length
    {
        return None;
    }
    let n = source_length as i64;
    let delta = i64::from(selected.enzyme.bottom_cut) - i64::from(selected.enzyme.top_cut);
    let observed = bottom as i64 - top as i64;
    if (topology == Topology::Circular && (observed - delta).rem_euclid(n) != 0)
        || (topology == Topology::Linear && observed != delta)
    {
        return None;
    }
    let coding = match transcript.coding_strand {
        ivt::CodingStrand::Forward => top as i64 - transcript.start as i64,
        ivt::CodingStrand::Reverse => transcript.start as i64 + 1 - bottom as i64,
    };
    let coding = if topology == Topology::Circular {
        coding.rem_euclid(n)
    } else {
        coding
    };
    Some(CutMapping {
        site_index,
        top_bond: top,
        bottom_bond: bottom,
        coding_cut_offset: coding,
        template_cut_offset: coding + delta,
        signed_template_minus_coding: delta,
    })
}

fn mapped_request(
    request: &Request,
    cut: &CutMapping,
    context: &str,
    enzyme_id: &str,
) -> Result<ivt::Request, Error> {
    let mut declaration = request.ivt.clone();
    declaration.linearization = Some(ivt::Linearization {
        coding_cut_offset: usize::try_from(cut.coding_cut_offset)
            .map_err(|_| Error::Format("negative coding cut offset".into()))?,
        template_cut_offset: usize::try_from(cut.template_cut_offset)
            .map_err(|_| Error::Format("negative template cut offset".into()))?,
        provenance: Provenance {
            statement: format!(
                "In-silico {enzyme_id} single-enzyme cut hypothesis under explicit source/positional profile; cutting is not established."
            ),
            reference: Some(format!(
                "collagene://ivt-linearization/{context}/{enzyme_id}"
            )),
        },
    });
    Ok(declaration)
}

fn contains_arc(
    source: &Molecule,
    container_start: i64,
    container_len: usize,
    start: i64,
    len: usize,
) -> bool {
    let distance = start - container_start;
    let distance = if source.topology == Topology::Circular {
        distance.rem_euclid(source.sequence.len() as i64)
    } else {
        distance
    };
    distance >= 0
        && (distance as usize)
            .checked_add(len)
            .is_some_and(|end| end <= container_len)
}
fn retains_template(
    source: &Molecule,
    transcript: &ivt::Transcript,
    fragment: &tracked::TrackedDigestFragment,
) -> bool {
    let (promoter, transcript_start) = match transcript.coding_strand {
        ivt::CodingStrand::Forward => (
            transcript.start as i64 - PROMOTER_BASES as i64,
            transcript.start as i64,
        ),
        ivt::CodingStrand::Reverse => (
            transcript.start as i64 + 1,
            transcript.start as i64 + 1 - transcript.length as i64,
        ),
    };
    let top = (fragment.top_start as i64, fragment.duplex.duplex.top.len());
    let bottom = (fragment.bottom_start, fragment.duplex.duplex.bottom.len());
    let template = if transcript.coding_strand == ivt::CodingStrand::Forward {
        bottom
    } else {
        top
    };
    contains_arc(source, top.0, top.1, promoter, PROMOTER_BASES)
        && contains_arc(source, bottom.0, bottom.1, promoter, PROMOTER_BASES)
        && contains_arc(
            source,
            template.0,
            template.1,
            transcript_start,
            transcript.length,
        )
}

fn one_conditions(
    options: &conditions::ReactionConditions,
    enzyme_id: &str,
) -> Result<conditions::ReactionConditions, Error> {
    let chosen = options
        .formulations
        .iter()
        .find(|f| f.enzyme_id == enzyme_id)
        .ok_or_else(|| Error::Format("missing explicit alternative formulation choice".into()))?;
    let mut one = options.clone();
    one.formulations = vec![chosen.clone()];
    Ok(one)
}

pub fn assess(source: &Molecule, request: &Request) -> Result<Report, Error> {
    let context = revision(source, request)?;
    let catalog = catalog::builtin_catalog()?;
    let (material, history_incomplete) = material(source, &request.source)?;
    let modifications = material.summary()?;
    let mut budget = CleavageBudget::new(request.limits.max_comparisons)?;
    // Reserve the request and source evidence before the first retained candidate.
    let mut charged = size(
        &(request, &modifications, STATEMENT),
        request.limits.max_report_bytes,
    )?;
    let mut candidates = Vec::with_capacity(request.selection.enzymes.len());
    for selected in &request.selection.enzymes {
        let selection = catalog::select(&catalog, &[selected.id.clone()])?;
        let checked = tracked::checked_digest_with_budget(
            &material,
            &CheckedDigestOptions {
                selection: selection.clone(),
                max_comparisons: request.limits.max_comparisons,
                max_sites: request.limits.max_sites_per_enzyme,
            },
            &mut budget,
        )?;
        let mut issues = Vec::new();
        if !checked.complete {
            issue(&mut issues, IssueKind::IncompleteScan, Severity::Unresolved);
        }
        if !checked.products_determined {
            issue(
                &mut issues,
                IssueKind::UnresolvedCleavage,
                Severity::Unresolved,
            );
        }
        let permitted: Vec<_> = checked
            .sites
            .iter()
            .enumerate()
            .filter(|(_, s)| s.outcome == Outcome::PermittedUnderProfile)
            .collect();
        if permitted.len() > 1 {
            issue(
                &mut issues,
                IssueKind::MultiplePermittedCuts,
                Severity::Reject,
            );
        }
        if permitted.is_empty() && checked.complete && checked.products_determined {
            issue(&mut issues, IssueKind::NoPermittedCut, Severity::Reject);
        }
        let condition_report = request
            .conditions
            .as_ref()
            .map(|options| {
                conditions::source::SourceReport::assess(
                    &catalog,
                    source,
                    &material,
                    &selection,
                    &one_conditions(options, &selected.id)?,
                )
            })
            .transpose()?;
        match condition_report
            .as_ref()
            .map(|c| c.report.documentary_status)
        {
            Some(conditions::EvidenceStatus::Known) => {}
            Some(
                conditions::EvidenceStatus::BelowThreshold
                | conditions::EvidenceStatus::TemperatureMismatch
                | conditions::EvidenceStatus::Conflict,
            ) => {
                issue(
                    &mut issues,
                    IssueKind::ConditionsBelowThreshold,
                    if request.policy.require_known_conditions {
                        Severity::Reject
                    } else {
                        Severity::Advisory
                    },
                );
            }
            _ => issue(
                &mut issues,
                IssueKind::ConditionsNotEstablished,
                if request.policy.require_known_conditions {
                    Severity::Unresolved
                } else {
                    Severity::Advisory
                },
            ),
        }
        if history_incomplete {
            issue(
                &mut issues,
                IssueKind::IncompleteOrigin,
                if request.policy.require_complete_origin {
                    Severity::Unresolved
                } else {
                    Severity::Advisory
                },
            );
        }
        let mut mapping = None;
        let mut retained_fragment = None;
        let mut ivt_report = None;
        if checked.complete && checked.products_determined && permitted.len() == 1 {
            let (index, site) = permitted[0];
            mapping = project(
                source.sequence.len(),
                source.topology,
                &request.ivt.transcript,
                selected,
                site,
                index,
            );
            if let Some(cut) = &mapping {
                let max_offset = match (source.topology, request.ivt.transcript.coding_strand) {
                    (Topology::Circular, _) => source.sequence.len(),
                    (_, ivt::CodingStrand::Forward) => {
                        source.sequence.len() - request.ivt.transcript.start
                    }
                    (_, ivt::CodingStrand::Reverse) => request.ivt.transcript.start + 1,
                } as i64;
                if cut.coding_cut_offset < 0
                    || cut.template_cut_offset < 0
                    || cut.coding_cut_offset > max_offset
                    || cut.template_cut_offset > max_offset
                {
                    issue(
                        &mut issues,
                        IssueKind::CutOffsetsOutsideSource,
                        Severity::Reject,
                    );
                } else {
                    let fragments: Vec<_> = checked
                        .fragments
                        .iter()
                        .enumerate()
                        .filter(|(_, f)| retains_template(source, &request.ivt.transcript, f))
                        .map(|(i, _)| i)
                        .collect();
                    if fragments.len() == 1 {
                        retained_fragment = Some(fragments[0]);
                    } else {
                        issue(
                            &mut issues,
                            IssueKind::PromoterOrTranscriptNotRetained,
                            Severity::Reject,
                        );
                    }
                    let declaration = mapped_request(request, cut, &context, &selected.id)?;
                    let report = ivt::assess(source, &declaration).map_err(rna_error)?;
                    if cut.template_cut_offset as usize != request.ivt.transcript.length {
                        issue(&mut issues, IssueKind::RunoffMismatch, Severity::Reject);
                    }
                    if !report
                        .linearization
                        .as_ref()
                        .is_some_and(|g| request.policy.allowed_end_kinds.contains(&g.kind))
                    {
                        issue(&mut issues, IssueKind::DisallowedEndKind, Severity::Reject);
                    }
                    if report
                        .findings
                        .iter()
                        .any(|f| f.verdict == ivt::Verdict::Fail)
                    {
                        issue(&mut issues, IssueKind::IvtFailure, Severity::Reject);
                    }
                    let promoter_known = report.findings.iter().any(|f| {
                        f.rule == ivt::Rule::ConsensusPromoter && f.verdict == ivt::Verdict::Pass
                    });
                    if !promoter_known
                        || (request.policy.require_resolved_ivt_checks
                            && report
                                .findings
                                .iter()
                                .any(|f| f.verdict == ivt::Verdict::Indeterminate))
                    {
                        issue(&mut issues, IssueKind::IvtUnresolved, Severity::Unresolved);
                    }
                    ivt_report = Some(report);
                }
            } else {
                issue(
                    &mut issues,
                    IssueKind::UnavailableCutCoordinates,
                    Severity::Unresolved,
                );
            }
        }
        let fragments = checked
            .fragments
            .iter()
            .enumerate()
            .map(|(index, fragment)| {
                Ok(Fragment {
                    index,
                    source_top_start: fragment.top_start,
                    source_bottom_start: fragment.bottom_start,
                    geometry: fragment.duplex.duplex.summary()?,
                })
            })
            .collect::<Result<Vec<_>, Error>>()?;
        let candidate = Candidate {
            enzyme_id: selected.id.clone(),
            state: state(&issues),
            issues,
            cleavage: checked.evidence(),
            fragments,
            conditions: condition_report,
            mapping,
            retained_fragment,
            ivt: ivt_report,
        };
        charged = charged
            .checked_add(size(
                &candidate,
                request.limits.max_report_bytes.saturating_sub(charged),
            )?)
            .ok_or(Error::TooLarge)?;
        candidates.push(candidate);
    }
    let mut report = Report {
        engine_version: ENGINE.into(),
        context_revision: context,
        source_revision: source.revision(),
        source_topology: source.topology,
        request: request.clone(),
        source_modifications: modifications,
        source_history_incomplete: history_incomplete,
        comparisons: budget.used(),
        candidates,
        statement: STATEMENT.into(),
        report_revision: String::new(),
    };
    report.report_revision = report.content_revision()?;
    size(&report, request.limits.max_report_bytes)?;
    report.validate_binding(&report.context_revision)?;
    Ok(report)
}

impl Report {
    fn content_revision(&self) -> Result<String, Error> {
        crate::encoding::bounded_json::hash(&(
            &self.engine_version,
            &self.context_revision,
            &self.source_revision,
            self.source_topology,
            &self.request,
            &self.source_modifications,
            self.source_history_incomplete,
            self.comparisons,
            &self.candidates,
            &self.statement,
        ))
    }
    /// Bounded documentary integrity only; fresh source claims require validate_for.
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        self.request.limits.validate()?;
        size(self, self.request.limits.max_report_bytes)?;
        self.request.policy.validate()?;
        let catalog = catalog::builtin_catalog()?;
        self.request.selection.validate_for(&catalog)?;
        let profile = ivt::profiles().map_err(rna_error)?.remove(0);
        if self.request.ivt.profile != profile.reference
            || self.request.ivt.linearization.is_some()
            || self.source_modifications.top_length == 0
            || self.source_modifications.top_length > ivt::MAX_TEMPLATE_BASES
            || self.source_modifications.top_length != self.source_modifications.bottom_length
        {
            return Err(Error::Format(
                "IVT source/profile documentary shape differs".into(),
            ));
        }
        if let Some(options) = &self.request.conditions {
            options.validate(&self.request.selection)?;
        }
        if self.engine_version != ENGINE
            || self.context_revision != expected
            || self.context_revision
                != crate::encoding::bounded_json::hash(&(
                    ENGINE,
                    &self.source_revision,
                    &self.request,
                ))?
            || self.source_revision != self.request.ivt.source_revision
            || self.report_revision != self.content_revision()?
            || self.statement != STATEMENT
            || self.candidates.len() != self.request.selection.enzymes.len()
            || self.comparisons > self.request.limits.max_comparisons
        {
            return Err(Error::Format(
                "IVT linearization report binding differs".into(),
            ));
        }
        let mut comparisons = 0usize;
        for (candidate, selected) in self.candidates.iter().zip(&self.request.selection.enzymes) {
            candidate.cleavage.validate()?;
            let one = catalog::select(&catalog, &[selected.id.clone()])?;
            if candidate.enzyme_id != selected.id
                || candidate.state != state(&candidate.issues)
                || candidate.issues.len() > 32
                || candidate.fragments.len() > MAX_SITES_PER_ENZYME + 1
                || candidate.cleavage.options.selection != one
                || candidate.cleavage.source_revision != self.source_modifications.state_revision
                || candidate.cleavage.profile_revision != tracked::profile_revision()
                || candidate.cleavage.engine_version != tracked::ENGINE_VERSION
                || candidate.cleavage.options.max_sites != self.request.limits.max_sites_per_enzyme
                || candidate.cleavage.options.max_comparisons != self.request.limits.max_comparisons
            {
                return Err(Error::Format("IVT alternative evidence differs".into()));
            }
            comparisons = comparisons
                .checked_add(candidate.cleavage.comparisons)
                .ok_or(Error::TooLarge)?;
            if let Some(conditions) = &candidate.conditions {
                conditions.validate_historical()?;
                let options = self
                    .request
                    .conditions
                    .as_ref()
                    .ok_or_else(|| Error::Format("unrequested IVT reaction choices".into()))?;
                if conditions.source_revision != self.source_revision
                    || conditions.material_revision != self.source_modifications.state_revision
                    || conditions.report.selection != one
                    || conditions.report.options != one_conditions(options, &selected.id)?
                {
                    return Err(Error::Format("IVT condition association differs".into()));
                }
            } else if self.request.conditions.is_some() {
                return Err(Error::Format("missing IVT conditions record".into()));
            }
            if let Some(report) = &candidate.ivt {
                report
                    .validate_binding(&report.context_revision)
                    .map_err(rna_error)?;
                let mapping = candidate
                    .mapping
                    .as_ref()
                    .ok_or_else(|| Error::Format("IVT declaration lacks cut mapping".into()))?;
                if report.request
                    != mapped_request(&self.request, mapping, &self.context_revision, &selected.id)?
                    || report.source_revision != self.source_revision
                    || report.source_topology != self.source_topology
                {
                    return Err(Error::Format("IVT candidate declaration differs".into()));
                }
            }
            if let Some(mapping) = &candidate.mapping {
                let site = candidate
                    .cleavage
                    .sites
                    .get(mapping.site_index)
                    .ok_or_else(|| Error::Format("IVT cut site index differs".into()))?;
                if site.outcome != Outcome::PermittedUnderProfile
                    || Some(mapping.clone())
                        != project(
                            self.source_modifications.top_length,
                            self.source_topology,
                            &self.request.ivt.transcript,
                            selected,
                            site,
                            mapping.site_index,
                        )
                {
                    return Err(Error::Format("IVT lifted cut mapping differs".into()));
                }
            }
            if candidate.state == CandidateState::EligibleUnderDeclarations
                && (!candidate.cleavage.complete
                    || !candidate.cleavage.products_determined
                    || candidate
                        .cleavage
                        .sites
                        .iter()
                        .filter(|s| s.outcome == Outcome::PermittedUnderProfile)
                        .count()
                        != 1
                    || candidate.mapping.is_none()
                    || candidate.ivt.is_none()
                    || candidate.retained_fragment.is_none())
            {
                return Err(Error::Format(
                    "eligible IVT candidate lacks its prerequisites".into(),
                ));
            }
            if candidate.state == CandidateState::EligibleUnderDeclarations {
                let report = candidate
                    .ivt
                    .as_ref()
                    .ok_or_else(|| Error::Format("missing IVT assessment".into()))?;
                let mapping = candidate
                    .mapping
                    .as_ref()
                    .ok_or_else(|| Error::Format("missing IVT cut mapping".into()))?;
                if mapping.template_cut_offset != self.request.ivt.transcript.length as i64
                    || !report
                        .linearization
                        .as_ref()
                        .is_some_and(|g| self.request.policy.allowed_end_kinds.contains(&g.kind))
                    || report
                        .findings
                        .iter()
                        .any(|f| f.verdict == ivt::Verdict::Fail)
                    || !report.findings.iter().any(|f| {
                        f.rule == ivt::Rule::ConsensusPromoter && f.verdict == ivt::Verdict::Pass
                    })
                    || (self.request.policy.require_resolved_ivt_checks
                        && report
                            .findings
                            .iter()
                            .any(|f| f.verdict == ivt::Verdict::Indeterminate))
                    || (self.request.policy.require_complete_origin
                        && self.source_history_incomplete)
                    || (self.request.policy.require_known_conditions
                        && !candidate
                            .conditions
                            .as_ref()
                            .is_some_and(|c| c.satisfies_request()))
                {
                    return Err(Error::Format(
                        "eligible IVT candidate contradicts its declared policy".into(),
                    ));
                }
            }
            if candidate
                .retained_fragment
                .is_some_and(|i| i >= candidate.fragments.len())
            {
                return Err(Error::Format("IVT retained fragment index differs".into()));
            }
        }
        if comparisons != self.comparisons {
            return Err(Error::Format("IVT comparison accounting differs".into()));
        }
        Ok(())
    }
    pub fn validate_for(&self, source: &Molecule, request: &Request) -> Result<(), Error> {
        self.validate_binding(&revision(source, request)?)?;
        let fresh = assess(source, request)?;
        if self.report_revision != fresh.report_revision {
            return Err(Error::Format("IVT linearization replay differs".into()));
        }
        Ok(())
    }
}

fn size(value: &impl Serialize, maximum: usize) -> Result<usize, Error> {
    struct Count {
        used: usize,
        maximum: usize,
    }
    impl Write for Count {
        fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
            self.used = self
                .used
                .checked_add(bytes.len())
                .filter(|n| *n <= self.maximum)
                .ok_or_else(|| std::io::Error::other("IVT candidate byte budget"))?;
            Ok(bytes.len())
        }
        fn flush(&mut self) -> std::io::Result<()> {
            Ok(())
        }
    }
    let mut count = Count { used: 0, maximum };
    serde_json::to_writer(&mut count, value).map_err(|_| Error::TooLarge)?;
    Ok(count.used)
}

#[cfg(test)]
mod tests;
