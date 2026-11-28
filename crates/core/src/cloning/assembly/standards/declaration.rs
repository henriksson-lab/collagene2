//! Advisory original-source declarations resolved into the existing standards checker.
//! No physical operation, implicit coordinate conversion or product gate lives here.
use super::{Report, Request, SourceBinding, SourcePart};
use crate::{
    encoding::bounded_json,
    planning::plan::{DesignPlan, PlanInput},
};
use collagene_assembly_standards::{assessment::Limits, profile::Profile};
use collagene_model::{Alphabet, Error};
use serde::{Deserialize, Serialize};
use std::{
    collections::{BTreeMap, BTreeSet},
    io::Write,
};

pub const DECLARED_ENGINE: &str = "declared-source-standard-review-1";
pub const MAX_DECLARATION_BYTES: usize = 128 * 1024;
pub const MAX_DECLARED_REVIEW_BYTES: usize = super::MAX_REPORT_BYTES + 64 * 1024;
const MAX_REQUEST_BYTES: usize = 16 * 1024 * 1024;
const MAX_PLAN_BYTES: usize = 32 * 1024 * 1024;
const SCOPE: &str = "Advisory format review of explicitly declared original-reference part intervals and body boundaries only. Native grammar, required/advisory site facts, additions, certainty and search limits are retained. This is not physical end validation, cleavage, assembly, coding-frame/function/expression or full-standard/vector compliance. Products, chemistry, independent target checks and global cloning status are unchanged by this advice. Saved-material/prepared-strand coordinates require a genuine native reference projection and are not inferred here.";
macro_rules! value {
    ($item:item) => {
        #[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
        #[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
        #[serde(deny_unknown_fields)]
        $item
    };
}
macro_rules! data {
    ($item:item) => {
        #[derive(Debug, Clone, Serialize, Deserialize)]
        #[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
        #[serde(deny_unknown_fields)]
        $item
    };
}
value! {
    pub struct ProfileRef {
        pub profile_id: String,
        pub profile_identity: String,
    }
}
value! {
    pub struct DeclaredPart {
        pub role: String,
        pub input_id: String,
        pub source_revision: String,
        pub start: usize,
        pub length: usize,
        pub reverse: bool,
        pub slot: String,
        pub body_start: usize,
        pub body_end: usize,
    }
}
value! {
    pub struct Declaration {
        pub profile: ProfileRef,
        #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, DeclaredPart, 128>")]
        pub parts: Vec<DeclaredPart>,
        pub limits: Limits,
    }
}
fn bad(message: &str) -> Error {
    Error::Format(format!("declared source standard: {message}"))
}
fn text(value: &str, maximum: usize) -> bool {
    !value.is_empty() && value.len() <= maximum && !value.chars().any(char::is_control)
}
fn hash(value: &str) -> bool {
    value.len() == 64
        && value
            .bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}
impl Declaration {
    /// Iterative fixed-size descriptors are admitted before any plan/source hash or clone.
    pub fn validate(&self) -> Result<(), Error> {
        if self.parts.is_empty() || self.parts.len() > 128 {
            return Err(Error::TooLarge);
        }
        bounded_json::preflight(self, MAX_DECLARATION_BYTES)?;
        self.limits.validate().map_err(super::lib_error)?;
        if !text(&self.profile.profile_id, 128) || !hash(&self.profile.profile_identity) {
            return Err(bad("invalid profile reference"));
        }
        let mut roles = BTreeSet::new();
        for part in &self.parts {
            if !text(&part.role, 96)
                || !roles.insert(part.role.as_str())
                || !text(&part.input_id, 128)
                || !hash(&part.source_revision)
                || !text(&part.slot, 96)
            {
                return Err(bad("invalid or duplicate role/input/profile declaration"));
            }
        }
        Ok(())
    }
    pub fn revision(&self) -> Result<String, Error> {
        self.validate()?;
        bounded_json::hash(&(DECLARED_ENGINE, self))
    }
}

/// One in-process allowance shared by native grammar reviews, including refusals.
/// Not serializable and never restored from retained documentary evidence.
#[derive(Debug)]
pub struct DeclarationBudget {
    limit: usize,
    charged: usize,
    observed: usize,
}
impl DeclarationBudget {
    pub fn new(limit: usize) -> Result<Self, Error> {
        if limit > 64_000_000 {
            return Err(Error::TooLarge);
        }
        Ok(Self {
            limit,
            charged: 0,
            observed: 0,
        })
    }
    /// Includes explicit reservations whose native failure returned no counter.
    pub fn charged(&self) -> usize {
        self.charged
    }
    /// Measured comparisons from successful, validated native reports only.
    pub fn observed_comparisons(&self) -> usize {
        self.observed
    }
    pub fn remaining(&self) -> usize {
        self.limit - self.charged
    }
    fn reserve(&mut self, amount: usize) -> Result<(), Error> {
        self.charged = self
            .charged
            .checked_add(amount)
            .filter(|n| *n <= self.limit)
            .ok_or(Error::TooLarge)?;
        Ok(())
    }
    fn settle(&mut self, reserved: usize, observed: usize) -> Result<(), Error> {
        if observed > reserved || reserved > self.charged {
            return Err(bad("native comparison ledger differs"));
        }
        let measured = self.observed.checked_add(observed).ok_or(Error::TooLarge)?;
        self.charged -= reserved - observed;
        self.observed = measured;
        Ok(())
    }
}
value! {
    pub enum Charge {
        NotRun,
        Observed { comparisons: usize },
        /// Native error returned without its internal counter. Retain the reserved
        /// allowance, explicitly not an assertion that all these comparisons ran.
        ReservedOnNativeRefusal { allowance: usize },
    }
}
impl Charge {
    pub fn charged(&self) -> usize {
        match self {
            Self::NotRun => 0,
            Self::Observed { comparisons } => *comparisons,
            Self::ReservedOnNativeRefusal { allowance } => *allowance,
        }
    }
}
value! {
    pub enum UnavailableKind {
        MissingInput,
        DuplicateInput,
        SourceRevisionMismatch,
        InvalidSource,
        InvalidInterval,
        InvalidBody,
        UnsupportedCoordinateScope,
        UnsupportedProfile,
        ProfileIdentityMismatch,
        ResourceLimit,
        SharedBudgetExhausted,
        NativeRefusal,
    }
}
value! {
    pub struct ReviewUnavailable {
        pub kind: UnavailableKind,
        pub part_index: Option<usize>,
        pub message: String,
    }
}
data! {
    pub enum DeclaredReviewState {
        Assessed(Report),
        Unavailable {
            #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, ReviewUnavailable, 128>")]
            reasons: Vec<ReviewUnavailable>,
        },
    }
}
data! {
    pub struct DeclaredReview {
        pub engine_version: String,
        pub plan_revision: String,
        pub declaration_revision: String,
        pub requested_limits: Limits,
        #[serde(deserialize_with = "crate::encoding::required_option")]
        #[cfg_attr(feature = "json-schema", schemars(required))]
        pub effective_limits: Option<Limits>,
        pub shared_remaining_before: usize,
        pub charge: Charge,
        pub state: DeclaredReviewState,
        pub scope: String,
        pub review_revision: String,
    }
}
fn unavailable(
    kind: UnavailableKind,
    part_index: Option<usize>,
    message: &str,
) -> DeclaredReviewState {
    DeclaredReviewState::Unavailable {
        reasons: vec![ReviewUnavailable {
            kind,
            part_index,
            message: message.into(),
        }],
    }
}
fn reason(
    reasons: &mut Vec<ReviewUnavailable>,
    kind: UnavailableKind,
    index: usize,
    message: &str,
) {
    reasons.push(ReviewUnavailable {
        kind,
        part_index: Some(index),
        message: message.into(),
    });
}
fn binding(part: &DeclaredPart, input: &PlanInput) -> SourceBinding {
    SourceBinding {
        id: part.role.clone(),
        revision: input.revision.clone(),
        start: part.start,
        length: part.length,
        reverse: part.reverse,
        source_length: input.molecule.sequence.len(),
        topology: input.molecule.topology,
    }
}
fn source_admission<'a>(
    plan: &'a DesignPlan,
    declaration: &Declaration,
) -> Result<Vec<&'a PlanInput>, DeclaredReviewState> {
    let mut inputs = BTreeMap::new();
    for input in &plan.inputs {
        if inputs.insert(input.id.as_str(), input).is_some() {
            return Err(unavailable(
                UnavailableKind::DuplicateInput,
                None,
                "Plan input identifiers are not unique; source format is unresolved.",
            ));
        }
    }
    let mut reasons = Vec::new();
    let mut resolved = Vec::with_capacity(declaration.parts.len());
    let mut bases = 0usize;
    for (index, part) in declaration.parts.iter().enumerate() {
        let Some(input) = inputs.get(part.input_id.as_str()).copied() else {
            reason(
                &mut reasons,
                UnavailableKind::MissingInput,
                index,
                "The exact declared input is unavailable.",
            );
            continue;
        };
        // One reason per part keeps the entire bounded set of missing roles visible.
        if input.revision != part.source_revision
            || input.molecule.revision() != part.source_revision
        {
            reason(
                &mut reasons,
                UnavailableKind::SourceRevisionMismatch,
                index,
                "The declared source revision differs from the exact plan input.",
            );
        } else if input.molecule.sequence.alphabet() != Alphabet::Dna
            || input.molecule.validate().is_err()
        {
            reason(
                &mut reasons,
                UnavailableKind::InvalidSource,
                index,
                "A valid DNA reference is required.",
            );
        } else if crate::chemistry::has_unmapped_physical_evidence(&input.molecule) {
            reason(
                &mut reasons,
                UnavailableKind::UnsupportedCoordinateScope,
                index,
                "Prepared/saved material requires a genuine native reference projection; stored strand coordinates are not inferred.",
            );
        } else if binding(part, input).validate().is_err() {
            reason(
                &mut reasons,
                UnavailableKind::InvalidInterval,
                index,
                "The declared original-reference interval is outside this source.",
            );
        } else if part.body_start >= part.body_end || part.body_end > part.length {
            reason(
                &mut reasons,
                UnavailableKind::InvalidBody,
                index,
                "Explicit oriented part-body boundaries must lie inside the selected interval.",
            );
        } else {
            bases = bases.checked_add(part.length).unwrap_or(usize::MAX);
            resolved.push(input);
        }
    }
    if !reasons.is_empty() {
        return Err(DeclaredReviewState::Unavailable { reasons });
    }
    if bases > declaration.limits.max_bases
        || declaration.parts.len() > declaration.limits.max_parts
    {
        return Err(unavailable(
            UnavailableKind::ResourceLimit,
            None,
            "Declared source parts exceed the explicit aggregate base/part bounds.",
        ));
    }
    // Measure the actual repeated request assets before any source clone; a
    // source used twice really is serialized twice by the existing Request owner.
    struct Remaining(usize);
    impl Write for Remaining {
        fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
            self.0 = self
                .0
                .checked_sub(bytes.len())
                .ok_or_else(|| std::io::Error::other("declared standard source request bound"))?;
            Ok(bytes.len())
        }
        fn flush(&mut self) -> std::io::Result<()> {
            Ok(())
        }
    }
    let mut count = Remaining(MAX_REQUEST_BYTES - MAX_DECLARATION_BYTES);
    for input in &resolved {
        if serde_json::to_writer(&mut count, &input.molecule).is_err() {
            return Err(unavailable(
                UnavailableKind::ResourceLimit,
                None,
                "Original source snapshots exceed the bounded native request representation.",
            ));
        }
    }
    Ok(resolved)
}

/// The optional declaration is advice within the one current plan, never a
/// physical operation or a second plan/protocol representation.
pub fn assess_declared(
    plan: &DesignPlan,
    budget: &mut DeclarationBudget,
) -> Result<Option<DeclaredReview>, Error> {
    let Some(declaration) = plan.standard.as_ref() else {
        return Ok(None);
    };
    declaration.validate()?;
    // Canonical enclosing-plan resource admission precedes local identity work.
    crate::material::capsule::plan_resource_preflight(plan)?;
    bounded_json::preflight(plan, MAX_PLAN_BYTES)?;
    let declaration_revision = declaration.revision()?;
    let mut review = DeclaredReview {
        engine_version: DECLARED_ENGINE.into(),
        plan_revision: plan.revision(),
        declaration_revision,
        requested_limits: declaration.limits.clone(),
        effective_limits: None,
        shared_remaining_before: budget.remaining(),
        charge: Charge::NotRun,
        state: unavailable(
            UnavailableKind::SharedBudgetExhausted,
            None,
            "Shared grammar-review allowance exhausted; no format/absence claim.",
        ),
        scope: SCOPE.into(),
        review_revision: String::new(),
    };
    if budget.remaining() != 0 {
        review.state = match evaluate(
            plan,
            declaration,
            budget,
            &mut review.effective_limits,
            &mut review.charge,
        ) {
            Ok(state) => state,
            Err(_) => {
                // An unexpected native retention refusal is advice, not a
                // physical failure. Keep any existing reservation charged.
                if review.charge == Charge::NotRun {
                    review.effective_limits = None;
                } else if let Some(limits) = &review.effective_limits {
                    review.charge = Charge::ReservedOnNativeRefusal {
                        allowance: limits.max_comparisons,
                    };
                }
                unavailable(
                    UnavailableKind::NativeRefusal,
                    None,
                    "Native advisory evaluation could not retain bounded evidence; no format conformity or absence claim.",
                )
            }
        };
    }
    // Native report payload is capped by super::MAX_REPORT_BYTES; this envelope
    // reserves an additional 64 KiB for its bounded pins/limits/charge/scope. If
    // late native association or retention fails, replace only advisory detail
    // with a compact typed refusal and retain the native reservation.
    let retained = review.identity().and_then(|revision| {
        review.review_revision = revision;
        review.validate_binding(&review.plan_revision, declaration)
    });
    if retained.is_err() {
        if let Some(limits) = &review.effective_limits {
            review.charge = Charge::ReservedOnNativeRefusal {
                allowance: limits.max_comparisons,
            };
        }
        review.state = unavailable(
            UnavailableKind::ResourceLimit,
            None,
            "Native advisory evidence could not be retained or associated within its bounds; no conformity claim.",
        );
        review.review_revision = review.identity()?;
        review.validate_binding(&review.plan_revision, declaration)?;
    }
    // Settle only AFTER final native association/hash/size validation. A late
    // refusal keeps the prior reservation even though no review is published.
    if let Charge::Observed { comparisons } = review.charge {
        let reserved = review
            .effective_limits
            .as_ref()
            .ok_or_else(|| bad("missing effective native allowance"))?
            .max_comparisons;
        budget.settle(reserved, comparisons)?;
    }
    Ok(Some(review))
}
fn evaluate(
    plan: &DesignPlan,
    declaration: &Declaration,
    budget: &mut DeclarationBudget,
    effective: &mut Option<Limits>,
    charge: &mut Charge,
) -> Result<DeclaredReviewState, Error> {
    let catalog = crate::cloning::restriction::catalog::builtin_catalog()?;
    let profile = match Profile::builtin(&declaration.profile.profile_id, &catalog) {
        Ok(profile) => profile,
        Err(_) => {
            return Ok(unavailable(
                UnavailableKind::UnsupportedProfile,
                None,
                "The declared installed profile is unavailable for the current catalogue.",
            ));
        }
    };
    if profile.identity() != declaration.profile.profile_identity {
        return Ok(unavailable(
            UnavailableKind::ProfileIdentityMismatch,
            None,
            "The exact installed profile/data identity changed; it was not repinned by name.",
        ));
    }
    let inputs = match source_admission(plan, declaration) {
        Ok(inputs) => inputs,
        Err(state) => return Ok(state),
    };
    let mut limits = declaration.limits.clone();
    limits.max_comparisons = limits.max_comparisons.min(budget.remaining());
    *effective = Some(limits.clone());
    let request = Request {
        profile_id: declaration.profile.profile_id.clone(),
        profile_identity: declaration.profile.profile_identity.clone(),
        parts: declaration
            .parts
            .iter()
            .zip(inputs)
            .map(|(part, input)| SourcePart {
                id: part.role.clone(),
                molecule: input.molecule.clone(),
                start: part.start,
                length: part.length,
                reverse: part.reverse,
                slot: part.slot.clone(),
                body_start: part.body_start,
                body_end: part.body_end,
            })
            .collect(),
        limits,
    };
    // The native checker shares one comparison allowance across all declared
    // parts. Reserve before calling so late native failure cannot reset it.
    let reserved = request.limits.max_comparisons;
    budget.reserve(reserved)?;
    *charge = Charge::ReservedOnNativeRefusal {
        allowance: reserved,
    };
    match super::assess(&request) {
        Ok(report) => {
            *charge = Charge::Observed {
                comparisons: report.assessment.comparisons,
            };
            Ok(DeclaredReviewState::Assessed(report))
        }
        Err(error) => Ok(unavailable(
            if matches!(error, Error::TooLarge) {
                UnavailableKind::ResourceLimit
            } else {
                UnavailableKind::NativeRefusal
            },
            None,
            &error.to_string().chars().take(1024).collect::<String>(),
        )),
    }
}
impl DeclaredReview {
    fn identity(&self) -> Result<String, Error> {
        bounded_json::preflight(
            &(
                &self.engine_version,
                &self.plan_revision,
                &self.declaration_revision,
                &self.requested_limits,
                &self.effective_limits,
                self.shared_remaining_before,
                &self.charge,
                &self.state,
                &self.scope,
            ),
            MAX_DECLARED_REVIEW_BYTES,
        )?;
        bounded_json::hash(&(
            &self.engine_version,
            &self.plan_revision,
            &self.declaration_revision,
            &self.requested_limits,
            &self.effective_limits,
            self.shared_remaining_before,
            &self.charge,
            &self.state,
            &self.scope,
        ))
    }
    /// Standalone current evidence shape/integrity. A plan-bound consumer must
    /// additionally call validate_binding with its exact declaration.
    pub fn validate_documentary(&self) -> Result<(), Error> {
        bounded_json::preflight(self, MAX_DECLARED_REVIEW_BYTES)?;
        self.requested_limits.validate().map_err(super::lib_error)?;
        if self.engine_version != DECLARED_ENGINE
            || !hash(&self.plan_revision)
            || !hash(&self.declaration_revision)
            || self.scope != SCOPE
            || self.review_revision != self.identity()?
            || self.shared_remaining_before > 64_000_000
            || self.charge.charged() > self.shared_remaining_before
        {
            return Err(bad("invalid standalone advisory evidence"));
        }
        if let Some(effective) = &self.effective_limits {
            effective.validate().map_err(super::lib_error)?;
            let mut expected = self.requested_limits.clone();
            expected.max_comparisons = expected.max_comparisons.min(self.shared_remaining_before);
            if effective != &expected || self.charge.charged() > effective.max_comparisons {
                return Err(bad("invalid effective advisory allowance"));
            }
        } else if self.charge != Charge::NotRun {
            return Err(bad("unexecuted advisory review has charges"));
        }
        match &self.state {
            DeclaredReviewState::Assessed(report) => {
                report.validate_binding(&report.context_revision)?;
                if Some(&report.assessment.limits) != self.effective_limits.as_ref()
                    || self.charge
                        != (Charge::Observed {
                            comparisons: report.assessment.comparisons,
                        })
                {
                    return Err(bad("native advisory charge/report disagreement"));
                }
            }
            DeclaredReviewState::Unavailable { reasons } => {
                if reasons.is_empty()
                    || reasons.len() > 128
                    || reasons.iter().any(|reason| {
                        !text(&reason.message, 4096)
                            || reason.part_index.is_some_and(|index| index >= 128)
                            || reason.kind == UnavailableKind::SharedBudgetExhausted
                                && self.shared_remaining_before != 0
                    })
                    || matches!(self.charge, Charge::Observed { .. })
                    || self.charge == Charge::NotRun && self.effective_limits.is_some()
                {
                    return Err(bad("invalid unavailable advisory evidence"));
                }
                if let Charge::ReservedOnNativeRefusal { allowance } = self.charge {
                    if self
                        .effective_limits
                        .as_ref()
                        .map(|limits| limits.max_comparisons)
                        != Some(allowance)
                    {
                        return Err(bad("invalid advisory refusal allowance"));
                    }
                }
            }
        }
        Ok(())
    }
    /// Bounded historical association/integrity, not current scientific authority.
    pub fn validate_binding(
        &self,
        expected_plan: &str,
        declaration: &Declaration,
    ) -> Result<(), Error> {
        declaration.validate()?;
        self.validate_documentary()?;
        if self.engine_version != DECLARED_ENGINE
            || !hash(expected_plan)
            || self.plan_revision != expected_plan
            || self.declaration_revision != declaration.revision()?
            || self.requested_limits != declaration.limits
            || self.scope != SCOPE
            || self.review_revision != self.identity()?
            || self.shared_remaining_before > 64_000_000
            || self.charge.charged() > self.shared_remaining_before
        {
            return Err(bad("review context/identity differs"));
        }
        if let Some(effective) = &self.effective_limits {
            effective.validate().map_err(super::lib_error)?;
            let mut expected = self.requested_limits.clone();
            expected.max_comparisons = expected.max_comparisons.min(self.shared_remaining_before);
            if effective != &expected || self.charge.charged() > effective.max_comparisons {
                return Err(bad("effective comparison allowance differs"));
            }
        } else if self.charge != Charge::NotRun {
            return Err(bad("unexecuted review has a native charge"));
        }
        match &self.state {
            DeclaredReviewState::Assessed(report) => {
                report.validate_binding(&report.context_revision)?;
                if report.profile.profile_id != declaration.profile.profile_id
                    || report.profile.profile_identity != declaration.profile.profile_identity
                    || Some(&report.assessment.limits) != self.effective_limits.as_ref()
                    || self.charge
                        != (Charge::Observed {
                            comparisons: report.assessment.comparisons,
                        })
                    || report.sources.len() != declaration.parts.len()
                {
                    return Err(bad("native report differs from declared review binding"));
                }
                for ((source, part), assessed) in report
                    .sources
                    .iter()
                    .zip(&declaration.parts)
                    .zip(&report.assessment.parts)
                {
                    if source.id != part.role
                        || source.revision != part.source_revision
                        || source.start != part.start
                        || source.length != part.length
                        || source.reverse != part.reverse
                        || assessed.slot != part.slot
                        || assessed.body_start != part.body_start
                        || assessed.body_end != part.body_end
                    {
                        return Err(bad("native source/body/slot correspondence differs"));
                    }
                }
            }
            DeclaredReviewState::Unavailable { reasons } => {
                if (self.charge == Charge::NotRun && self.effective_limits.is_some())
                    || reasons
                        .iter()
                        .any(|reason| reason.kind == UnavailableKind::SharedBudgetExhausted)
                        && self.shared_remaining_before != 0
                {
                    return Err(bad(
                        "unavailable review has contradictory execution evidence",
                    ));
                }
                if reasons.is_empty()
                    || reasons.len() > 128
                    || reasons.iter().any(|r| {
                        !text(&r.message, 4096)
                            || r.part_index.is_some_and(|i| i >= declaration.parts.len())
                    })
                    || matches!(self.charge, Charge::Observed { .. })
                {
                    return Err(bad("invalid unavailable-review evidence"));
                }
                if let Charge::ReservedOnNativeRefusal { allowance } = self.charge {
                    if self.effective_limits.as_ref().map(|l| l.max_comparisons) != Some(allowance)
                    {
                        return Err(bad("refusal reservation differs from effective allowance"));
                    }
                }
            }
        }
        Ok(())
    }
    /// Fresh checking is an execution/worker operation, never a GUI render hook.
    pub fn validate_for(&self, plan: &DesignPlan) -> Result<(), Error> {
        let declaration = plan
            .standard
            .as_ref()
            .ok_or_else(|| bad("plan has no declared standard"))?;
        declaration.validate()?;
        crate::material::capsule::plan_resource_preflight(plan)?;
        bounded_json::preflight(plan, MAX_PLAN_BYTES)?;
        self.validate_binding(&plan.revision(), declaration)?;
        let fresh = assess_declared(
            plan,
            &mut DeclarationBudget::new(self.shared_remaining_before)?,
        )?
        .ok_or_else(|| bad("declared review unavailable"))?;
        if fresh.review_revision != self.review_revision {
            return Err(bad("fresh declared review differs"));
        }
        Ok(())
    }
    pub fn map_interval(
        &self,
        part_index: usize,
        start: usize,
        end: usize,
    ) -> Result<Vec<(usize, usize)>, Error> {
        match &self.state {
            DeclaredReviewState::Assessed(report) => report
                .sources
                .get(part_index)
                .ok_or_else(|| bad("part index unavailable"))?
                .map_interval(start, end),
            DeclaredReviewState::Unavailable { .. } => Err(Error::Unsupported(
                "No resolved original-reference interval is available.".into(),
            )),
        }
    }
}
