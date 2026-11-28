//! Fresh deterministic reconstruction over the one current typed plan.
use super::checks::{CheckFinding, CheckScope, CheckVerdict, RuleId};
use super::{
    ENGINE_VERSION,
    evidence::{
        DerivedOriginEvidence, LigationOutcomeBinding, LigationOutcomeEvidence,
        MaterialStateEvidence, PhysicalProcessingEvidence, PhysicalProcessingStepReport,
        ReferenceStepEvidence, ReferenceStepReport, Status, StepOutcome, StepReport, StepState,
        VerificationReport, evidence_bytes,
    },
    pcr::{project_mutant_pcr, project_pcr_frame_annotations},
    projection::{project_into, project_into_with_collapse, region_map, reverse_reference_map},
    sources::{
        ReferenceChecks, declared_source_termini_preflight, recombination_inputs,
        retain_topo_source_incompleteness, validate_topo_shape,
    },
    terminal::{
        annealed_document_inner, carry_unmapped_evidence, repaired_document,
        simulate_oligo_annealing, tailed_document,
    },
};
use crate::{
    chemistry::{
        self, Duplex,
        tracked_duplex::{self as modified_duplex, CheckedDigestOptions, TrackedDuplex},
    },
    planning::plan::{DesignPlan, Operation, Source},
    sequence::analysis,
};
use collagene_model::annotations::MappingSegment;
use collagene_model::{Alphabet, Error, MAX_SEQUENCE_LEN, Molecule, Topology};
use std::collections::{BTreeMap, BTreeSet};
mod execution;
mod operation;
use super::replay_budget::ReplayBudget;
use execution::Execution;
pub fn verify(plan: &DesignPlan, target: &Molecule) -> VerificationReport {
    evaluate(plan, target).0
}

/// Fresh private material state from the same reconstruction as the report.
/// Only in-process recipe adapters can access it; serialized product metadata or
/// imported reports cannot supply this vector. Every slot retains its original
/// plan-step index; failed/reference-only steps supply no private material.
pub(crate) fn verify_with_physical(
    plan: &DesignPlan,
    target: &Molecule,
) -> (VerificationReport, Vec<Option<Duplex>>) {
    let (report, physical) = evaluate(plan, target);
    (
        report,
        physical
            .into_iter()
            .map(|state| state.map(|tracked| tracked.duplex))
            .collect(),
    )
}

/// Keep the full freshly reconstructed strand-local state for serial cloning.
/// This is deliberately private to the native core: a saved report, metadata or
/// deserialized tracked state cannot substitute for this reconstruction pass.
pub(crate) fn verify_with_tracked(
    plan: &DesignPlan,
    target: &Molecule,
) -> (VerificationReport, Vec<Option<TrackedDuplex>>) {
    evaluate(plan, target)
}

/// Capture the actual final Golden Gate enumeration from the same native pass,
/// including truncated/withheld pools without an established selected product.
/// This diagnostic receipt is separate from the private material slots and does
/// not confer product/chemical authority. Serialized reports cannot supply it.
pub(crate) fn verify_with_golden_gate_trace(
    plan: &DesignPlan,
    target: &Molecule,
) -> (
    VerificationReport,
    Vec<Option<Duplex>>,
    Option<crate::cloning::golden_gate::GoldenGateReport>,
) {
    let mut captured = None;
    let (report, physical) = evaluate_captured(
        plan,
        target,
        &mut captured,
        &mut modified_duplex::CleavageBudget::new(100_000_000)
            .expect("fixed native cleavage bound"),
    );
    (
        report,
        physical.into_iter().map(|s| s.map(|t| t.duplex)).collect(),
        captured,
    )
}

fn evaluate(
    plan: &DesignPlan,
    target: &Molecule,
) -> (VerificationReport, Vec<Option<TrackedDuplex>>) {
    evaluate_captured(
        plan,
        target,
        &mut None,
        &mut modified_duplex::CleavageBudget::new(100_000_000)
            .expect("fixed native cleavage bound"),
    )
}
/// Execute a native scientific scope even when no coordinate receipt is needed.
/// QC callers retain the same borrowed owner across nested source replays.
pub(crate) fn verify_with_tracked_scoped(
    plan: &DesignPlan,
    target: &Molecule,
    scientific: &mut ReplayBudget<'_>,
) -> (VerificationReport, Vec<Option<TrackedDuplex>>) {
    let (report, material, _) = evaluate_native(plan, target, &mut None, scientific, None);
    (report, material)
}
fn evaluate_captured(
    plan: &DesignPlan,
    target: &Molecule,
    final_golden_gate: &mut Option<crate::cloning::golden_gate::GoldenGateReport>,
    cleavage_budget: &mut modified_duplex::CleavageBudget,
) -> (VerificationReport, Vec<Option<TrackedDuplex>>) {
    let mut scientific = ReplayBudget::plan_local(cleavage_budget);
    let (report, material, _) =
        evaluate_native(plan, target, final_golden_gate, &mut scientific, None);
    (report, material)
}

/// One fresh execution, with private coordinate receipts produced only by its
/// actual native operations. Receipt refusal never changes ordinary outcomes.
pub(crate) fn verify_with_ancestry(
    plan: &DesignPlan,
    target: &Molecule,
) -> (
    VerificationReport,
    Result<super::ancestry::ReconstructionTrace, Error>,
) {
    let mut cleavage =
        modified_duplex::CleavageBudget::new(100_000_000).expect("fixed native cleavage bound");
    let mut scientific = ReplayBudget::plan_local(&mut cleavage);
    let (report, _, trace) = evaluate_native(
        plan,
        target,
        &mut None,
        &mut scientific,
        Some(
            super::ancestry::TraceBudget::new(super::ancestry::Limits::default())
                .expect("fixed trace limits"),
        ),
    );
    (report, trace.expect("requested native ancestry result"))
}

/// Nested retained-origin replay shares the enclosing actual trace budget.
/// The returned private receipts come from the same pass as tracked material.
pub(crate) fn verify_with_tracked_ancestry_budget(
    plan: &DesignPlan,
    target: &Molecule,
    scientific: &mut ReplayBudget<'_>,
    trace_budget: super::ancestry::TraceBudget,
) -> (
    VerificationReport,
    Vec<Option<TrackedDuplex>>,
    Result<super::ancestry::ReconstructionTrace, Error>,
) {
    let (report, material, trace) =
        evaluate_native(plan, target, &mut None, scientific, Some(trace_budget));
    (
        report,
        material,
        trace.expect("requested shared native ancestry result"),
    )
}

fn evaluate_native(
    plan: &DesignPlan,
    target: &Molecule,
    final_golden_gate: &mut Option<crate::cloning::golden_gate::GoldenGateReport>,
    scientific: &mut ReplayBudget<'_>,
    capture_ancestry: Option<super::ancestry::TraceBudget>,
) -> (
    VerificationReport,
    Vec<Option<TrackedDuplex>>,
    Option<Result<super::ancestry::ReconstructionTrace, Error>>,
) {
    let mut report = VerificationReport {
        standard_review: None,
        plan_revision: String::new(), target_revision: String::new(),
        engine_version: ENGINE_VERSION.into(), status: Status::Invalid,
        findings: vec![], planned_steps: 0, steps: vec![], exhaustive_alternative_search: false,
        step_reports: vec![], physical_processing_reports: vec![], reference_step_reports: vec![],
        assumptions: vec!["Exact annealing/overlap model; no experimental outcome is guaranteed.".into(),
            "RecombineReference is a bounded full-site-profile reciprocal reference exchange. All components are computed before an explicit selection; neither its selected DNA nor target agreement establishes chemical strand state, biological repair, efficacy or a selected experimental winner. Reference-only ancestry cannot supply physical PCR or annealing material.".into(),
            "TOPO freshly replays declared strands, pairing and the ordered transfer, flap-cleavage and nick-sealing transcript. Saved graphs retain all calls and consumed activations. Topo and CheckedTopo supply a restrictive reference view; only TopoActions can bridge a proven representable contour to native tracked material. No implicit repair, supplier-kit state, isolation or experimental success is asserted.".into(),
            "The typed Pcr operations share one enclosing-plan 200M actual-work ledger for exact binding, native coverage and copied product bases; failed searches retain charges. Other method engines, intentional-mismatch PCR and nested saved-origin replay retain their explicitly bounded native scopes.".into(),
            "Exact-PCR operations consider exact non-overlapping binding sites within one template turn, including competing oligo orientations and single-oligo products. Intentional-mismatch PCR uses its separately recorded mismatch and 3′ policy; ambiguous and overlapping-domain products remain incomplete previews.".into(),
            "Golden Gate enumerates a bounded one-copy fragment pool under explicit options; recut paths, carryover and alternative source histories remain evidence, not a claim of exhaustive kinetic mixture uniqueness.".into(),
            "Ordinary digest is an explicit complete-cleavage assumption. Golden Gate and CheckedDigest evaluate actual tracked positions under the current pinned sensitivity profile; reaction conditions, secondary structure and genome-wide specificity are not established.".into(),
            "Ordered overlap/ligation candidates are checked for consistency; alternative mixture products are not exhaustively searched.".into(),
            "Uncatalogued Digest assumes complete cleavage with 5′-phosphate/3′-OH ends. DeclareSourceModifications applies only to original inputs under a closed Dam/Dcm/CpG generation model; it is not a methylation measurement or evidence about other modification systems. Tracked descendants preserve actual strand states; unknown/impaired checked sites withhold fully determined products.".into(),
            "Exact PCR physical intermediates are ideal blunt duplexes with unmodified 5′-OH/3′-OH termini. Polymerase-specific terminal additions and modified primers are not inferred; single-nucleotide additions require an explicit Tail3Prime operation.".into(),
            "End repair performs only explicit 3′ fill-in or removal under its accessible-OH policy; 5′ chemistry and internal nicks are not repaired. Oligo annealing assumes independently declared linear ssDNA inputs and explicitly supplied terminal chemistry, searches exact antiparallel placements, and forms no covalent bonds.".into(),
            "Single-nucleotide 3′ tailing is an explicit ideal blunt/OH endpoint transform, not predicted polymerase activity. Method-specific overlap recommendations are advisory sequence checks, not proprietary reagent or experimental compatibility evidence.".into(),
            "Product annotations are projected from source coordinates; clipped, interrupted and uncertain features require review. Unmapped annotations remain source-graph evidence. Imported metadata cannot assert end chemistry.".into()],
    };
    let admission = execution::admit(plan, target);
    if let Err(error) = admission {
        let detail = error.to_string();
        let message = if detail.len() <= 3500 {
            format!("Request was refused before plan hashing or material reconstruction: {detail}")
        } else {
            "Request failed bounded admission before hashing or reconstruction; overlong diagnostic detail was refused explicitly.".into()
        };
        report.findings.push(CheckFinding::new(
            RuleId::Admission,
            CheckVerdict::Fail,
            CheckScope::Plan,
            message,
        ));
        report.status = Status::Invalid;
        return (
            report,
            vec![],
            capture_ancestry.map(|_| {
                Err(Error::Unsupported(
                    "Native ancestry request failed bounded plan admission".into(),
                ))
            }),
        );
    }
    match crate::cloning::assembly::standards::declaration::assess_declared(
        plan,
        &mut scientific.standards,
    ) {
        Ok(review) => report.standard_review = review,
        Err(_) => {
            report.findings.push(CheckFinding::new(
                RuleId::Admission,
                CheckVerdict::Fail,
                CheckScope::Plan,
                "Bounded advisory evidence envelope could not be constructed after resource admission; the request was refused before reconstruction.",
            ));
            return (
                report,
                vec![],
                capture_ancestry.map(|_| Err(Error::TooLarge)),
            );
        }
    }
    report.plan_revision = plan.revision();
    report.target_revision = target.revision();
    report.planned_steps = plan.steps.len();
    report.steps = (0..plan.steps.len())
        .map(|step_index| StepOutcome {
            step_index,
            state: StepState::NotEvaluated,
            product: None,
            blocked_by: vec![],
        })
        .collect();
    let (trace, trace_error) = match capture_ancestry
        .map(|budget| super::ancestry::TraceBuilder::with_budget(plan.steps.len(), budget))
    {
        Some(Ok(mut trace)) => match trace.bind_plan(&report.plan_revision) {
            Ok(()) => (Some(trace), None),
            Err(error) => (None, Some(error)),
        },
        Some(Err(error)) => (None, Some(error)),
        None => (None, None),
    };
    let mut execution = Execution::new(report, final_golden_gate, scientific, trace, trace_error);
    execution.reconstruct(plan, target);
    execution.finish()
}

/// Coordinate comparison only; this helper confers no physical or report authority.
/// Callers bind and bound models separately. Ordinary plans retain their oriented
/// comparison; final CheckedTopo uses only its explicitly declared allowances.
pub fn plan_target_matches(plan: &DesignPlan, product: &Molecule, target: &Molecule) -> bool {
    if product.topology != target.topology
        || product.sequence.alphabet() != target.sequence.alphabet()
    {
        return false;
    }
    let policy = match plan.steps.last() {
        Some(Operation::CheckedTopo { policy, .. }) => Some(policy),
        _ => None,
    };
    let same = |candidate: &str| {
        candidate == target.sequence.as_str()
            || target.topology == Topology::Circular
                && policy.is_none_or(|p| p.allow_circular_rotation)
                && analysis::rotation_offset(candidate, target.sequence.as_str()).is_some()
    };
    same(product.sequence.as_str())
        || policy.is_some_and(|p| p.allow_reverse_complement)
            && product
                .sequence
                .reverse_complement()
                .is_ok_and(|rc| same(rc.as_str()))
}

fn assemble_overlap(
    fragments: &[&Molecule],
    overlaps: &[usize],
    circular: bool,
    index: usize,
    findings: &mut Vec<String>,
    ancestry: &mut super::ancestry::native::NativeStep<'_>,
    source_tokens: &[Option<super::ancestry::SnapshotToken>],
    product_token: &mut Option<super::ancestry::SnapshotToken>,
) -> Result<Molecule, Error> {
    if fragments.is_empty()
        || fragments.len() > 1000
        || overlaps.len() != fragments.len() - 1 + usize::from(circular)
    {
        return Err(Error::Format(
            "ordered assembly requires one overlap per junction".into(),
        ));
    }
    for f in fragments {
        if f.topology != Topology::Linear || f.sequence.alphabet() != Alphabet::Dna {
            return Err(Error::Unsupported(
                "overlap inputs must be linear DNA fragments".into(),
            ));
        }
        if f.sequence
            .as_str()
            .bytes()
            .any(|b| analysis::base_mask(b).count_ones() != 1)
        {
            return Err(Error::Unsupported("ambiguous assembly sequence".into()));
        }
    }
    for (i, fragment) in fragments.iter().enumerate() {
        let incoming = if i > 0 {
            overlaps[i - 1]
        } else if circular {
            *overlaps.last().unwrap()
        } else {
            0
        };
        let outgoing = if i + 1 < fragments.len() {
            overlaps[i]
        } else if circular {
            *overlaps.last().unwrap()
        } else {
            0
        };
        if incoming.saturating_add(outgoing) >= fragment.sequence.len() {
            return Err(Error::Unsupported(format!(
                "fragment {} has overlapping or exhaustive overlap arms; a more detailed assembly model is required",
                i + 1
            )));
        }
    }
    let mut sequence = fragments[0].sequence.as_str().to_owned();
    for (i, fragment) in fragments.iter().enumerate().skip(1) {
        let overlap = overlaps[i - 1];
        let next = fragment.sequence.as_str();
        let previous = fragments[i - 1].sequence.as_str();
        if overlap == 0
            || overlap >= previous.len()
            || overlap >= next.len()
            || previous[previous.len() - overlap..] != next[..overlap]
        {
            return Err(Error::Format(format!(
                "step {}: invalid overlap at junction {i}",
                index + 1
            )));
        }
        if sequence.len() + next.len() - overlap > MAX_SEQUENCE_LEN {
            return Err(Error::TooLarge);
        }
        sequence.push_str(&next[overlap..]);
    }
    if circular {
        let overlap = overlaps[overlaps.len() - 1];
        if overlap == 0
            || overlap >= fragments[0].sequence.len()
            || overlap >= fragments.last().unwrap().sequence.len()
            || overlap >= sequence.len()
            || sequence[sequence.len() - overlap..] != sequence[..overlap]
        {
            return Err(Error::Format(
                "circular closure overlap does not match".into(),
            ));
        }
        sequence.truncate(sequence.len() - overlap);
    }
    let mut product = Molecule::dna(
        format!("Assembly {}", index + 1),
        &sequence,
        if circular {
            Topology::Circular
        } else {
            Topology::Linear
        },
    )?;
    let mut offset = 0usize;
    for (i, fragment) in fragments.iter().enumerate() {
        let length = fragment.sequence.len();
        let n = product.sequence.len();
        let destination = offset % n;
        let first = length.min(n - destination);
        let mut map = vec![MappingSegment {
            source_start: 0,
            source_end: first,
            destination_start: destination,
            reverse: false,
        }];
        if first < length {
            map.push(MappingSegment {
                source_start: first,
                source_end: length,
                destination_start: 0,
                reverse: false,
            });
        }
        let native_projection = project_into_with_collapse(
            fragment,
            &map,
            &mut product,
            findings,
            circular && fragments.len() == 1,
        )?;
        ancestry.project(
            source_tokens.get(i).copied().flatten(),
            fragment,
            product_token,
            &product,
            &map,
            &native_projection,
            super::ancestry::ProjectionKind::HomologyCollapse,
            &[],
        );
        offset += length;
        if i < overlaps.len() {
            offset -= overlaps[i];
        }
    }
    Ok(product)
}
