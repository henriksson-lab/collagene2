//! Bound current-plan admission, independent branch scheduling and atomic commits.
use super::super::checks::{MAX_FINDING_BYTES, MAX_FINDINGS};
use super::operation::{MaterialBudget, StepAttempt, operation_rule};
use super::*;
use std::cell::Cell;

/// Resource admission deliberately does not validate current input identity or
/// operation semantics: those checks are accumulated in their own scopes below.
pub(super) fn admit(plan: &DesignPlan, target: &Molecule) -> Result<(), Error> {
    if plan.steps.len() > 1000 || plan.inputs.len() > 1000 || plan.name.len() > 1024 {
        return Err(Error::TooLarge);
    }
    let mut references = 0usize;
    for step in &plan.steps {
        let count = crate::planning::editing::reference_count(step);
        if count > 1000 {
            return Err(Error::TooLarge);
        }
        references = references
            .checked_add(count)
            .filter(|n| *n <= 100_000)
            .ok_or(Error::TooLarge)?;
        if step.sources().iter().any(|source| {
            matches!(source,
            Source::Input(id) if id.len() > 128 || id.chars().any(char::is_control))
        }) {
            return Err(Error::TooLarge);
        }
    }
    if plan
        .inputs
        .iter()
        .any(|input| input.id.len() > 128 || input.id.chars().any(char::is_control))
    {
        return Err(Error::TooLarge);
    }
    let mut models = plan
        .inputs
        .iter()
        .map(|input| &input.molecule)
        .collect::<Vec<_>>();
    models.push(target);
    collagene_model::bounds::model_preflight(&models)?;
    evidence_bytes(target, 32 * 1024 * 1024)?;
    crate::material::capsule::plan_resource_preflight(plan)?;
    ReferenceChecks::new(plan)?;
    for operation in &plan.steps {
        match operation {
            Operation::RecombineReference { request, .. } => {
                request.profile.profile.resource_preflight()?;
                crate::cloning::recombination::reference::limits_preflight(&request.limits)?;
                request.annotation_limits.validate()?;
            }
            Operation::MultiSiteRecombineReference { request, .. } => {
                request.reaction.resource_preflight()?;
                request.annotation_limits.validate()?;
            }
            Operation::TopoActions { request, .. } => {
                crate::cloning::topo::actions::resource_preflight(request)?;
            }
            Operation::Topo { request, .. } | Operation::CheckedTopo { request, .. } => {
                crate::cloning::topo::physical::steps_resource_preflight(
                    &request.steps,
                    &request.engine_limits,
                )
                .map_err(|error| {
                    Error::Unsupported(format!("TOPO transcript resources: {error}"))
                })?;
            }
            _ => (),
        }
    }
    Ok(())
}

/// Optional detail cannot consume the reserved per-step outcome/limit evidence.
/// Refusal is explicit; no overlong message is silently clipped into a verdict.
struct FindingBudget {
    bytes: usize,
    reserved_count: usize,
    reserved_bytes: usize,
    exhausted: bool,
}
impl FindingBudget {
    fn new(steps: usize) -> Self {
        Self {
            bytes: 2,
            reserved_count: steps + 4,
            reserved_bytes: (steps + 4) * 1024,
            exhausted: false,
        }
    }
    fn append(
        &mut self,
        report: &mut VerificationReport,
        finding: CheckFinding,
        mandatory: bool,
    ) -> bool {
        let count_limit = if mandatory {
            MAX_FINDINGS
        } else {
            MAX_FINDINGS - self.reserved_count
        };
        let bytes_limit = if mandatory {
            MAX_FINDING_BYTES
        } else {
            MAX_FINDING_BYTES - self.reserved_bytes
        };
        let size = finding
            .validate()
            .and_then(|()| evidence_bytes(&finding, bytes_limit));
        let Ok(size) = size else {
            self.exhausted = true;
            return false;
        };
        if report.findings.len() >= count_limit || self.bytes.saturating_add(size + 1) > bytes_limit
        {
            self.exhausted = true;
            return false;
        }
        self.bytes += size + 1;
        report.findings.push(finding);
        true
    }
}

pub(super) struct Execution<'a, 'scope> {
    report: VerificationReport,
    ancestry: Option<super::super::ancestry::TraceBuilder>,
    ancestry_error: Option<Error>,
    physical: Vec<Option<TrackedDuplex>>,
    reference_only: Vec<bool>,
    final_golden_gate: &'a mut Option<crate::cloning::golden_gate::GoldenGateReport>,
    scientific: &'a mut ReplayBudget<'scope>,
    finding_budget: FindingBudget,
    material_budget: MaterialBudget,
    limited: bool,
    retained_bytes: usize,
}
impl<'a, 'scope> Execution<'a, 'scope> {
    pub fn new(
        report: VerificationReport,
        final_golden_gate: &'a mut Option<crate::cloning::golden_gate::GoldenGateReport>,
        scientific: &'a mut ReplayBudget<'scope>,
        ancestry: Option<super::super::ancestry::TraceBuilder>,
        ancestry_error: Option<Error>,
    ) -> Self {
        let count = report.planned_steps;
        let retained_bytes = evidence_bytes(&report, 60 * 1024 * 1024).unwrap_or(60 * 1024 * 1024);
        Self {
            report,
            ancestry,
            ancestry_error,
            physical: vec![None; count],
            reference_only: vec![false; count],
            final_golden_gate,
            scientific,
            finding_budget: FindingBudget::new(count),
            material_budget: MaterialBudget::new(),
            limited: false,
            retained_bytes,
        }
    }
    fn finding(&mut self, finding: CheckFinding) {
        if !self.finding_budget.append(&mut self.report, finding, false) {
            self.limited = true;
        }
    }
    fn mandatory(&mut self, finding: CheckFinding) {
        // These fixed small messages fit the admission-time reserved headroom.
        let accepted = self.finding_budget.append(&mut self.report, finding, true);
        debug_assert!(accepted, "reserved mandatory verification evidence");
        self.limited |= !accepted;
    }
    fn step_finding(
        &mut self,
        index: usize,
        rule: RuleId,
        verdict: CheckVerdict,
        message: impl Into<String>,
    ) {
        self.finding(CheckFinding::new(
            rule,
            verdict,
            CheckScope::Step { step_index: index },
            message,
        ));
    }
    fn not_evaluated(&mut self, index: usize, reason: &str) {
        self.report.steps[index].state = StepState::NotEvaluated;
        self.mandatory(CheckFinding::new(
            RuleId::StepCompletion,
            CheckVerdict::Indeterminate,
            CheckScope::Step { step_index: index },
            reason,
        ));
    }
    pub fn reconstruct(&mut self, plan: &DesignPlan, target: &Molecule) {
        if plan.steps.is_empty() {
            self.finding(CheckFinding::new(
                RuleId::PlanShape,
                CheckVerdict::Fail,
                CheckScope::Plan,
                "A scientific cloning plan requires at least one operation.",
            ));
            return;
        }
        let target_bound = target.validate().is_ok() && plan.target_revision == target.revision();
        if !target_bound {
            self.finding(CheckFinding::new(RuleId::TargetBinding, CheckVerdict::Fail, CheckScope::Plan,
                "The independent target is invalid or its exact revision differs from the plan binding. Independent source and operation checks continue; target-dependent acceptance is unavailable."));
        }
        let mut counts = BTreeMap::<&str, usize>::new();
        for input in &plan.inputs {
            *counts.entry(&input.id).or_default() += 1;
        }
        let mut inputs = BTreeMap::new();
        let mut invalid_inputs = BTreeSet::new();
        for input in &plan.inputs {
            let scope = CheckScope::Input {
                input_id: input.id.clone(),
            };
            let mut valid = true;
            if input.id.trim().is_empty() || counts[input.id.as_str()] != 1 {
                valid = false;
                if invalid_inputs.insert(input.id.clone()) {
                    self.finding(CheckFinding::new(RuleId::InputIdentity, CheckVerdict::Fail, scope.clone(),
                        "Original input identity is empty or duplicated; no ambiguous entry is selected."));
                }
            }
            if let Err(error) = input.molecule.validate() {
                valid = false;
                self.finding(CheckFinding::new(
                    RuleId::InputShape,
                    CheckVerdict::Fail,
                    scope.clone(),
                    format!("Original input molecule is invalid: {error}"),
                ));
            }
            if input.revision != input.molecule.revision() {
                valid = false;
                self.finding(CheckFinding::new(
                    RuleId::InputRevision,
                    CheckVerdict::Fail,
                    scope,
                    "Original input revision differs from the current immutable molecule.",
                ));
            }
            if valid {
                inputs.insert(input.id.clone(), &input.molecule);
            } else {
                invalid_inputs.insert(input.id.clone());
            }
        }
        if let Some(trace) = &mut self.ancestry {
            for (id, molecule) in &inputs {
                if let Err(error) = trace.register_input(id, molecule) {
                    self.ancestry_error = Some(error);
                    break;
                }
            }
        }
        if self.ancestry_error.is_some() {
            self.ancestry = None;
        }
        let checks = match ReferenceChecks::new(plan) {
            Ok(checks) => checks,
            Err(_) => {
                self.limited = true;
                self.mandatory(CheckFinding::new(
                    RuleId::Admission,
                    CheckVerdict::Indeterminate,
                    CheckScope::Plan,
                    "The admitted reference catalogue/resource context could not be established.",
                ));
                for index in 0..plan.steps.len() {
                    self.not_evaluated(index, "The current reference/resource context is unavailable; no operation was evaluated.");
                }
                return;
            }
        };
        let plan_hash = self.report.plan_revision.clone();
        for (index, operation) in plan.steps.iter().enumerate() {
            if self.limited {
                self.not_evaluated(index, "Shared job admission, work or evidence budget was exhausted; this operation was not evaluated.");
                continue;
            }
            let finding_start = self.report.findings.len();
            if !target_bound && matches!(operation, Operation::CheckedTopo { .. }) {
                self.step_finding(
                    index,
                    RuleId::TopoTargetPolicy,
                    CheckVerdict::Fail,
                    "CheckedTopo requires the separately bound valid target; its target-dependent acceptance cannot run against an invalid or changed target.",
                );
            }
            // Check independent primer declarations even when a dependency/profile
            // fails. Do not simulate annealing or invent a usable substrate here.
            self.primer_admission(index, operation);
            let preflight_rule = Cell::new(operation_rule(operation));
            if let Err(error) = checks.validate_step(operation, &preflight_rule) {
                self.step_finding(
                    index,
                    preflight_rule.get(),
                    error_verdict(&error),
                    format!("Operation declaration was refused: {error}"),
                );
            }
            let mut blockers = Vec::new();
            let mut seen_sources = BTreeSet::new();
            for source in operation.sources() {
                let key = match source {
                    Source::Input(id) => (0u8, id.as_str(), 0usize),
                    Source::Step(producer) => (1u8, "", *producer),
                };
                if !seen_sources.insert(key) {
                    continue;
                }
                match source {
                    Source::Input(id) if invalid_inputs.contains(id) => {
                        blockers.push(source.clone());
                    }
                    Source::Input(id) if !inputs.contains_key(id) => {
                        self.finding(
                            CheckFinding::new(
                                RuleId::SourceDependency,
                                CheckVerdict::Fail,
                                CheckScope::Step { step_index: index },
                                "Operation references an unknown original input.",
                            )
                            .with_subjects(vec![source.clone()]),
                        );
                    }
                    Source::Step(producer) if *producer >= index => {
                        self.finding(CheckFinding::new(RuleId::SourceDependency, CheckVerdict::Fail,
                            CheckScope::Step { step_index: index }, "Step dependencies must refer strictly to earlier original plan indices.")
                            .with_subjects(vec![source.clone()]));
                    }
                    Source::Step(producer) if self.report.product(*producer).is_none() => {
                        blockers.push(source.clone());
                    }
                    _ => (),
                }
            }
            let failed = self.report.findings[finding_start..]
                .iter()
                .any(|f| f.verdict == CheckVerdict::Fail);
            let withheld = self.report.findings[finding_start..]
                .iter()
                .any(|f| f.verdict == CheckVerdict::Indeterminate);
            if self.limited {
                self.not_evaluated(index, "The shared check-evidence budget cannot retain this operation's complete findings; no product is committed.");
                continue;
            }
            if failed || withheld {
                self.report.steps[index].state = if failed {
                    StepState::Rejected
                } else {
                    StepState::Withheld
                };
                if !blockers.is_empty() {
                    self.finding(CheckFinding::new(RuleId::DependencyBlocked, CheckVerdict::Blocked,
                        CheckScope::Step { step_index: index }, "Independent declaration errors were reported; these unavailable source dependencies also prevent reconstruction.")
                        .with_subjects(blockers));
                    if self.limited {
                        self.not_evaluated(index, "The shared evidence budget cannot retain all direct dependency details alongside the independent errors; no product is asserted.");
                    }
                }
                continue;
            }
            if !blockers.is_empty() {
                let blocker_bytes =
                    evidence_bytes(&blockers, 60 * 1024 * 1024).unwrap_or(60 * 1024 * 1024);
                if self.retained_bytes.saturating_add(blocker_bytes) > 60 * 1024 * 1024 {
                    self.limited = true;
                    self.not_evaluated(index, "Direct blocker evidence exceeds the enclosing report budget; no complete blocking list or product is asserted.");
                    continue;
                }
                self.retained_bytes += blocker_bytes;
                self.report.steps[index].state = StepState::Blocked;
                self.report.steps[index].blocked_by = blockers;
                self.mandatory(CheckFinding::new(RuleId::DependencyBlocked, CheckVerdict::Blocked,
                    CheckScope::Step { step_index: index }, "Reconstruction is blocked by the explicit unavailable original input or earlier producing step(s) retained on this outcome."));
                continue;
            }
            if matches!(operation, Operation::Pcr { .. }) && self.scientific.pcr.remaining() == 0 {
                self.not_evaluated(index, "The shared native exact-PCR work budget is exhausted; no fresh binding search or product was evaluated.");
                continue;
            }
            let attempt = StepAttempt::run(
                plan,
                target,
                index,
                &plan_hash,
                &inputs,
                &self.report.steps,
                &self.physical,
                &self.reference_only,
                &mut self.material_budget,
                self.report.physical_processing_reports.len(),
                self.scientific,
                self.ancestry.as_mut(),
            );
            self.commit(index, operation, attempt);
            if self.report.steps[index].state != StepState::Reconstructed {
                if let Some(trace) = &mut self.ancestry {
                    if let Err(error) = trace.discard_step(index) {
                        self.ancestry_error = Some(error);
                    }
                }
            }
        }
        if self.finding_budget.exhausted {
            self.mandatory(CheckFinding::new(RuleId::CheckBudget, CheckVerdict::Indeterminate, CheckScope::Plan,
                "The shared named-check count, field or byte budget was reached. Detail was refused explicitly; remaining operations were not evaluated and no incomplete finding list establishes consistency."));
        }
        self.compare_target(plan, target, target_bound);
    }
    fn primer_admission(&mut self, index: usize, operation: &Operation) {
        let (forward, reverse, exact) = match operation {
            Operation::Pcr {
                forward, reverse, ..
            } => (forward, reverse, true),
            Operation::PcrWithMismatches {
                forward, reverse, ..
            } => (forward, reverse, false),
            _ => return,
        };
        for (role, primer) in [("Forward", forward), ("Reverse", reverse)] {
            let error: Option<(CheckVerdict, String)> = if primer.name.len() > 1024
                || primer.annealing.len() > 4096
                || primer.tail.len() > 4096
            {
                Some((
                    CheckVerdict::Indeterminate,
                    "primer exceeds the bounded 4096-base annealing/tail or 1024-byte name profile"
                        .into(),
                ))
            } else if primer.annealing.alphabet() != Alphabet::Dna {
                Some((
                    CheckVerdict::Fail,
                    "annealing sequence must use DNA alphabet".into(),
                ))
            } else if let Err(error) =
                analysis::Primer::new(&primer.name, primer.annealing.as_str(), &primer.tail)
            {
                Some((
                    CheckVerdict::Fail,
                    format!("primer declaration is invalid: {error}"),
                ))
            } else if exact
                && primer
                    .annealing
                    .as_str()
                    .bytes()
                    .any(|b| !analysis::base_mask(b).is_power_of_two())
            {
                Some((
                    CheckVerdict::Indeterminate,
                    "degenerate annealing is outside the exact-PCR profile".into(),
                ))
            } else {
                None
            };
            if let Some((verdict, message)) = error {
                self.step_finding(
                    index,
                    RuleId::PcrPrimer,
                    verdict,
                    format!("{role} {message}"),
                );
            }
        }
    }
    fn commit(&mut self, index: usize, operation: &Operation, attempt: StepAttempt) {
        let StepAttempt {
            result,
            rule,
            notes,
            incomplete,
            step_reports,
            physical_processing_reports,
            reference_step_reports,
            golden_gate,
            ancestry_error,
        } = attempt;
        if let Some(error) = ancestry_error {
            self.ancestry_error = Some(error);
            self.ancestry = None;
        }
        // The fresh final enumeration is a separate diagnostic receipt. A
        // truncated pool or unavailable selected candidate must remain visible
        // even when no product/private material can be committed. Consumers
        // inspect its own completeness flags; this receipt grants no chemistry.
        if index + 1 == self.report.planned_steps
            && matches!(
                operation,
                Operation::GoldenGate { .. } | Operation::SelectedGoldenGate { .. }
            )
        {
            *self.final_golden_gate = golden_gate;
        }
        let domain_rule = operation_rule(operation);
        for message in notes {
            self.step_finding(index, domain_rule, CheckVerdict::Informational, message);
            if self.limited {
                break;
            }
        }
        for message in incomplete {
            self.step_finding(index, domain_rule, CheckVerdict::Indeterminate, message);
            if self.limited {
                break;
            }
        }
        // Useful failed/withheld native evidence remains explicitly step-scoped;
        // only successful material is committed to private slots below.
        let staged_bytes = evidence_bytes(
            &(
                &step_reports,
                &physical_processing_reports,
                &reference_step_reports,
            ),
            60 * 1024 * 1024,
        )
        .and_then(|bytes| {
            let product_bytes = match &result {
                Ok((product, _, _)) => evidence_bytes(product, 60 * 1024 * 1024)?,
                Err(_) => 0,
            };
            bytes
                .checked_add(product_bytes)
                .and_then(|bytes| bytes.checked_add(128))
                .ok_or(Error::TooLarge)
        });
        let stage_fits = staged_bytes
            .as_ref()
            .is_ok_and(|bytes| self.retained_bytes.saturating_add(*bytes) <= 60 * 1024 * 1024)
            && self.report.physical_processing_reports.len() + physical_processing_reports.len()
                <= 64
            && self.report.reference_step_reports.len() + reference_step_reports.len() <= 64
            && self.report.step_reports.len() + step_reports.len() <= 1000;
        if !stage_fits {
            self.limited = true;
            self.report.steps[index].state = StepState::Withheld;
            self.mandatory(CheckFinding::new(RuleId::EvidenceBudget, CheckVerdict::Indeterminate,
                CheckScope::Step { step_index: index }, "Complete staged product/evidence cannot fit the enclosing report bounds. No product, tracked state or partial scientific evidence was committed; remaining expensive work is not evaluated."));
            return;
        }
        self.retained_bytes += staged_bytes.unwrap_or(0);
        self.report.step_reports.extend(step_reports);
        self.report
            .physical_processing_reports
            .extend(physical_processing_reports);
        self.report
            .reference_step_reports
            .extend(reference_step_reports);
        if self.limited {
            self.not_evaluated(index, "The shared check-evidence budget cannot retain this operation's full result; product and physical state are withheld.");
            return;
        }
        match result {
            Ok((product, tracked, reference_only)) => {
                self.report.steps[index].state = StepState::Reconstructed;
                self.report.steps[index].product = Some(product);
                self.physical[index] = tracked;
                self.reference_only[index] = reference_only;
                self.mandatory(CheckFinding::new(RuleId::StepCompletion, CheckVerdict::Pass,
                    CheckScope::Step { step_index: index }, "This operation's product was freshly reconstructed at its original step index; physical/reference-only limitations remain separate findings."));
            }
            Err(error) => {
                let verdict = error_verdict(&error);
                self.report.steps[index].state = if verdict == CheckVerdict::Fail {
                    StepState::Rejected
                } else {
                    StepState::Withheld
                };
                if matches!(rule, RuleId::ProductBudget | RuleId::EvidenceBudget) {
                    self.limited = true;
                }
                self.step_finding(
                    index,
                    rule,
                    verdict,
                    format!("Operation reconstruction was refused: {error}"),
                );
                if self.finding_budget.exhausted {
                    self.not_evaluated(index, "The operation error exceeded the shared named-check evidence bound; no product was committed.");
                }
            }
        }
    }
    fn compare_target(&mut self, plan: &DesignPlan, target: &Molecule, target_bound: bool) {
        let index = plan.steps.len() - 1;
        if !target_bound {
            self.mandatory(CheckFinding::new(
                RuleId::TargetAgreement,
                CheckVerdict::Blocked,
                CheckScope::Target { step_index: index },
                "The independent target binding failed. No product was compared against that target and no matching target or alternative target is inferred.",
            ));
            return;
        }
        if let Some(product) = self.report.final_product() {
            let matches = plan_target_matches(plan, product, target);
            self.mandatory(CheckFinding::new(RuleId::TargetAgreement,
                if matches { CheckVerdict::Pass } else { CheckVerdict::Fail },
                CheckScope::Target { step_index: index },
                if matches { "The designated final product matches the separately bound target under the current orientation/topology comparison policy. Independent branch errors and physical limitations remain separate checks." }
                else { "The designated final product differs from the independently bound target under the current sequence/orientation/topology policy." }));
        } else {
            self.mandatory(CheckFinding::new(RuleId::TargetAgreement, CheckVerdict::Blocked,
                CheckScope::Target { step_index: index }, "Target comparison is blocked because the designated final step has no established product; another successful branch is not substituted."));
        }
    }
    pub fn finish(
        mut self,
    ) -> (
        VerificationReport,
        Vec<Option<TrackedDuplex>>,
        Option<Result<super::super::ancestry::ReconstructionTrace, Error>>,
    ) {
        self.report.status = if self
            .report
            .findings
            .iter()
            .any(|f| f.verdict == CheckVerdict::Fail)
        {
            Status::Invalid
        } else if self.report.planned_steps > 0
            && self.report.reconstructed_count() == self.report.planned_steps
            && !self.report.findings.iter().any(|f| {
                matches!(
                    f.verdict,
                    CheckVerdict::Indeterminate | CheckVerdict::Blocked
                )
            })
            && self
                .report
                .findings
                .iter()
                .any(|f| f.rule == RuleId::TargetAgreement && f.verdict == CheckVerdict::Pass)
        {
            Status::ConsistentUnderAssumptions
        } else {
            Status::Incomplete
        };
        let trace = match (self.ancestry_error, self.ancestry) {
            (Some(error), _) => Some(Err(error)),
            (None, Some(trace)) => Some(Ok(trace.finish())),
            (None, None) => None,
        };
        (self.report, self.physical, trace)
    }
}
fn error_verdict(error: &Error) -> CheckVerdict {
    match error {
        Error::Unsupported(_) | Error::TooLarge => CheckVerdict::Indeterminate,
        _ => CheckVerdict::Fail,
    }
}
