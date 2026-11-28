//! Literal source-format expectations, independent of preparation/product success.
use super::declaration::*;
use crate::{
    cloning::assembly::standards,
    planning::plan::{DesignPlan, PlanInput},
};
use collagene_assembly_standards::{
    assessment::{FindingKind, Limits, Outcome},
    profile::{MOCLO_2011_L0_ID, RFC10_ID, RFC1000_L0_ID},
};
use collagene_model::{Molecule, Topology};

fn plan(profile_id: &str, sources: &[(&str, &str)]) -> DesignPlan {
    let profile = standards::catalogue()
        .unwrap()
        .into_iter()
        .find(|p| p.profile_id == profile_id)
        .unwrap();
    let inputs: Vec<_> = sources
        .iter()
        .enumerate()
        .map(|(index, (_, sequence))| {
            let molecule =
                Molecule::dna(format!("literal input {index}"), sequence, Topology::Linear)
                    .unwrap();
            PlanInput {
                id: format!("input-{index}"),
                revision: molecule.revision(),
                molecule,
            }
        })
        .collect();
    let parts = inputs
        .iter()
        .zip(sources)
        .enumerate()
        .map(|(index, (input, (slot, _)))| DeclaredPart {
            role: format!("role-{index}"),
            input_id: input.id.clone(),
            source_revision: input.revision.clone(),
            start: 0,
            length: input.molecule.sequence.len(),
            reverse: false,
            slot: (*slot).into(),
            body_start: 0,
            body_end: input.molecule.sequence.len(),
        })
        .collect();
    DesignPlan {
        name: "advisory literal source standard".into(),
        target_revision: inputs[0].revision.clone(),
        inputs,
        steps: vec![],
        standard: Some(Declaration {
            profile: ProfileRef {
                profile_id: profile.profile_id,
                profile_identity: profile.profile_identity,
            },
            parts,
            limits: Limits::default(),
        }),
    }
}
fn review(plan: &DesignPlan) -> DeclaredReview {
    assess_declared(plan, &mut DeclarationBudget::new(64_000_000).unwrap())
        .unwrap()
        .unwrap()
}
fn native(review: &DeclaredReview) -> &standards::Report {
    match &review.state {
        DeclaredReviewState::Assessed(report) => report,
        other => panic!("expected native report: {other:?}"),
    }
}
fn unavailable(review: &DeclaredReview, kind: UnavailableKind) {
    match &review.state {
        DeclaredReviewState::Unavailable { reasons } => {
            assert!(reasons.iter().any(|reason| reason.kind == kind))
        }
        other => panic!("expected typed unavailable state: {other:?}"),
    }
}
#[test]
fn rfc10_requires_literal_formatted_additions_while_rfc1000_takes_declared_body_only() {
    const PREFIX: &str = "GAATTCGCGGCCGCTTCTAGAG";
    const BODY: &str = "ACACAC";
    const SUFFIX: &str = "TACTAGTAGCGGCCGCTGCAG";
    let formatted = format!("{PREFIX}{BODY}{SUFFIX}");
    let mut rfc10 = plan(RFC10_ID, &[("noncoding", &formatted)]);
    let part = &mut rfc10.standard.as_mut().unwrap().parts[0];
    part.body_start = PREFIX.len();
    part.body_end = PREFIX.len() + BODY.len();
    let result = review(&rfc10);
    assert_eq!(native(&result).assessment.outcome, Outcome::Conforms);
    assert_eq!(native(&result).assessment.parts[0].body_start, PREFIX.len());
    assert!(
        native(&review(&plan(RFC10_ID, &[("noncoding", BODY)])))
            .assessment
            .findings
            .iter()
            .any(|finding| matches!(finding.kind, FindingKind::Addition { .. }))
    );
    let body_only = review(&plan(RFC1000_L0_ID, &[("promoter", BODY)]));
    assert_eq!(native(&body_only).assessment.outcome, Outcome::Conforms);
    // The real RFC10 prefix/suffix are not silently stripped for a body profile.
    let mut incompatible = plan(RFC1000_L0_ID, &[("promoter", &formatted)]);
    incompatible.standard.as_mut().unwrap().parts[0].body_start = PREFIX.len();
    incompatible.standard.as_mut().unwrap().parts[0].body_end = PREFIX.len() + BODY.len();
    assert_eq!(
        native(&review(&incompatible)).assessment.outcome,
        Outcome::Violations
    );
}
#[test]
fn whole_moclo_four_and_five_module_paths_and_rfc1000_subpath_use_existing_order_grammar() {
    let four = plan(
        MOCLO_2011_L0_ID,
        &[
            ("promoter", "ACACAC"),
            ("five-prime-utr", "ACACAC"),
            ("cytosolic-cds", "ACACAC"),
            ("terminator", "ACACAC"),
        ],
    );
    let five = plan(
        MOCLO_2011_L0_ID,
        &[
            ("promoter", "ACACAC"),
            ("five-prime-utr", "ACACAC"),
            ("signal-peptide", "ACACAC"),
            ("cds-after-signal", "ACACAC"),
            ("terminator", "ACACAC"),
        ],
    );
    for source in [&four, &five] {
        let result = review(source);
        assert_eq!(native(&result).assessment.outcome, Outcome::Conforms);
        assert!(native(&result).successors.may_finish);
        assert!(result.scope.contains("not physical end validation"));
        assert!(result.scope.contains("full-standard/vector compliance"));
    }
    let mut wrong = four.clone();
    wrong.standard.as_mut().unwrap().parts.swap(1, 2);
    let result = review(&wrong);
    assert_eq!(native(&result).assessment.outcome, Outcome::Violations);
    assert!(
        native(&result)
            .assessment
            .findings
            .iter()
            .any(|f| matches!(f.kind, FindingKind::Grammar(_)))
    );
    // A complete four-module requirement must not be imposed on the installed
    // RFC1000 profile, whose declared nonempty contiguous subpaths are allowed.
    let subpath = plan(
        RFC1000_L0_ID,
        &[("five-prime-utr", "ACACAC"), ("cds", "ATGACACAC")],
    );
    assert_eq!(
        native(&review(&subpath)).assessment.outcome,
        Outcome::Conforms
    );
}
#[test]
fn reversed_circular_interval_and_real_native_finding_map_to_original_source_coordinates() {
    // Source[7..10]+source[0..3] is literal GGTCTC; reverse complement GAGACC.
    let mut source = plan(RFC1000_L0_ID, &[("promoter", "CTCAAAAGGT")]);
    source.inputs[0].molecule.topology = Topology::Circular;
    source.inputs[0].revision = source.inputs[0].molecule.revision();
    let input_revision = source.inputs[0].revision.clone();
    let part = &mut source.standard.as_mut().unwrap().parts[0];
    part.source_revision = input_revision;
    part.start = 7;
    part.length = 6;
    part.body_end = 6;
    part.reverse = true;
    let result = review(&source);
    let report = native(&result);
    assert_eq!(
        report.assessment.parts[0].sequence_blake3,
        blake3::hash(b"GAGACC").to_hex().to_string()
    );
    let finding = report
        .assessment
        .findings
        .iter()
        .find_map(|finding| match finding.kind {
            FindingKind::AvoidedSite {
                start,
                end,
                reverse,
                ..
            } => Some((start, end, reverse)),
            _ => None,
        })
        .unwrap();
    assert_eq!(finding, (0, 6, true));
    assert_eq!(
        result.map_interval(0, finding.0, finding.1).unwrap(),
        [(7, 10), (0, 3)]
    );
    assert_eq!(result.map_interval(0, 1, 5).unwrap(), [(8, 10), (0, 2)]);
}
#[test]
fn profile_source_body_and_missing_role_bindings_are_explicit_unavailable_not_new_geometry() {
    let source = plan(RFC1000_L0_ID, &[("promoter", "ACACAC")]);
    let mut wrong = source.clone();
    wrong.standard.as_mut().unwrap().profile.profile_identity = "0".repeat(64);
    unavailable(&review(&wrong), UnavailableKind::ProfileIdentityMismatch);
    let mut wrong = source.clone();
    wrong.standard.as_mut().unwrap().profile.profile_id = RFC10_ID.into();
    unavailable(&review(&wrong), UnavailableKind::ProfileIdentityMismatch);
    let mut wrong = source.clone();
    wrong.standard.as_mut().unwrap().profile.profile_id = "not-an-installed-standard".into();
    unavailable(&review(&wrong), UnavailableKind::UnsupportedProfile);
    let mut wrong = source.clone();
    wrong.standard.as_mut().unwrap().parts[0].source_revision = "0".repeat(64);
    unavailable(&review(&wrong), UnavailableKind::SourceRevisionMismatch);
    let mut wrong = source.clone();
    wrong.standard.as_mut().unwrap().parts[0].input_id = "missing-retained-role".into();
    unavailable(&review(&wrong), UnavailableKind::MissingInput);
    let mut wrong = source.clone();
    wrong.standard.as_mut().unwrap().parts[0].body_end = 7;
    unavailable(&review(&wrong), UnavailableKind::InvalidBody);
    let mut wrong = source.clone();
    wrong.standard.as_mut().unwrap().parts[0].start = 6;
    unavailable(&review(&wrong), UnavailableKind::InvalidInterval);
    let mut unknown_slot = source.clone();
    unknown_slot.standard.as_mut().unwrap().parts[0].slot = "undeclared-module-type".into();
    assert_eq!(
        native(&review(&unknown_slot)).assessment.outcome,
        Outcome::Violations
    );
}
#[test]
fn marker_evidence_does_not_become_reference_coordinates_and_advice_does_not_edit_inputs() {
    let mut source = plan(RFC1000_L0_ID, &[("promoter", "ACACAC")]);
    source.inputs[0].molecule.metadata.insert(
        "duplex_model".into(),
        "untrusted physical documentary state".into(),
    );
    source.inputs[0].revision = source.inputs[0].molecule.revision();
    source.standard.as_mut().unwrap().parts[0].source_revision = source.inputs[0].revision.clone();
    let before = serde_json::to_value(&source).unwrap();
    let result = review(&source);
    unavailable(&result, UnavailableKind::UnsupportedCoordinateScope);
    assert_eq!(serde_json::to_value(&source).unwrap(), before);
    assert!(result.map_interval(0, 0, 1).is_err());
}
#[test]
fn ambiguity_and_shared_exhaustion_never_become_source_format_passes() {
    let ambiguous = plan(RFC1000_L0_ID, &[("promoter", "NNNNNNNN")]);
    let result = review(&ambiguous);
    assert_ne!(native(&result).assessment.outcome, Outcome::Conforms);
    let source = plan(RFC1000_L0_ID, &[("promoter", "ACACACACACACACACACACACAC")]);
    let mut budget = DeclarationBudget::new(1).unwrap();
    let limited = assess_declared(&source, &mut budget).unwrap().unwrap();
    assert!(!native(&limited).assessment.complete);
    assert_eq!(budget.charged(), 1);
    assert_eq!(budget.observed_comparisons(), 1);
    assert_eq!(
        limited.effective_limits.as_ref().unwrap().max_comparisons,
        1
    );
    assert_eq!(
        limited.requested_limits.max_comparisons,
        Limits::default().max_comparisons
    );
    let exhausted = assess_declared(&source, &mut budget).unwrap().unwrap();
    unavailable(&exhausted, UnavailableKind::SharedBudgetExhausted);
    assert_eq!(exhausted.charge, Charge::NotRun);
    assert_eq!(budget.charged(), 1);
}
#[test]
fn late_native_refusal_reservation_is_distinct_from_measured_comparisons() {
    // 128 required source bindings alone cannot fit the native4KiB report.
    // Even after native findings are trimmed, metadata admission returns Err
    // without its consumed comparison counter. This is a real native refusal.
    let sources = vec![("promoter", "ACACAC"); 128];
    let mut source = plan(RFC1000_L0_ID, &sources);
    let limits = &mut source.standard.as_mut().unwrap().limits;
    limits.max_parts = 128;
    limits.max_report_bytes = 4096;
    limits.max_comparisons = 1000;
    let mut budget = DeclarationBudget::new(2000).unwrap();
    let result = assess_declared(&source, &mut budget).unwrap().unwrap();
    unavailable(&result, UnavailableKind::ResourceLimit);
    assert_eq!(
        result.charge,
        Charge::ReservedOnNativeRefusal { allowance: 1000 }
    );
    assert_eq!(budget.charged(), 1000);
    assert_eq!(budget.observed_comparisons(), 0);
    assert_eq!(budget.remaining(), 1000);
}
#[test]
fn current_declaration_association_and_integrity_are_checked_without_restoring_authority() {
    let source = plan(RFC1000_L0_ID, &[("promoter", "ACACAC")]);
    let result = review(&source);
    let declaration = source.standard.as_ref().unwrap();
    result
        .validate_binding(&source.revision(), declaration)
        .unwrap();
    result.validate_for(&source).unwrap();
    let mut changed = declaration.clone();
    changed.parts[0].reverse = true;
    assert!(
        result
            .validate_binding(&source.revision(), &changed)
            .is_err()
    );
    let mut tampered = result.clone();
    tampered.scope = "full standard and physical success verified".into();
    assert!(
        tampered
            .validate_binding(&source.revision(), declaration)
            .is_err()
    );
    let mut absent = source.clone();
    absent.standard = None;
    let mut budget = DeclarationBudget::new(100).unwrap();
    assert!(assess_declared(&absent, &mut budget).unwrap().is_none());
    assert_eq!(budget.charged(), 0);
}

#[test]
fn verifier_advice_does_not_change_products_or_scientific_status_and_is_required_nullable() {
    use crate::{
        planning::plan::{Operation, Source},
        verification::engine,
    };
    let mut declared = plan(RFC10_ID, &[("noncoding", "ACACAC")]);
    declared.steps = vec![Operation::Overlap {
        fragments: vec![Source::Input(declared.inputs[0].id.clone())],
        overlaps: vec![],
        circular: false,
    }];
    let target = declared.inputs[0].molecule.clone();
    let mut ordinary = declared.clone();
    ordinary.standard = None;
    let baseline = engine::verify(&ordinary, &target);
    let result = engine::verify(&declared, &target);
    assert!(baseline.final_product().is_some());
    assert_eq!(result.status, baseline.status);
    let ordinary_revision = ordinary.revision();
    let declared_revision = declared.revision();
    assert_ne!(ordinary_revision, declared_revision);
    let mut expected = baseline.final_product().unwrap().clone();
    assert_eq!(
        expected.metadata.get("plan_revision"),
        Some(&ordinary_revision)
    );
    assert_eq!(
        result
            .final_product()
            .unwrap()
            .metadata
            .get("plan_revision"),
        Some(&declared_revision)
    );
    // The changed declaration changes its genuine plan identity. Every other
    // product field, including DNA, annotations and physical evidence, is exact.
    expected
        .metadata
        .insert("plan_revision".into(), declared_revision);
    assert_eq!(result.final_product(), Some(&expected));
    let advisory = result.standard_review.as_ref().unwrap();
    advisory.validate_for(&declared).unwrap();
    assert_ne!(native(advisory).assessment.outcome, Outcome::Conforms);
    assert!(baseline.standard_review.is_none());
    result.validate_documentary().unwrap();
    result.validate_standard_for_plan(&declared).unwrap();
    assert!(result.validate_standard_for_plan(&ordinary).is_err());
    let mut missing_advice = result.clone();
    missing_advice.standard_review = None;
    assert!(
        missing_advice
            .validate_standard_for_plan(&declared)
            .is_err()
    );
    let mut unadmitted = declared.clone();
    unadmitted.name = "x".repeat(1025);
    let refused = engine::verify(&unadmitted, &target);
    assert!(refused.plan_revision.is_empty());
    assert!(refused.standard_review.is_none());
    refused.validate_standard_for_plan(&unadmitted).unwrap();

    let mut encoded = serde_json::to_value(&ordinary).unwrap();
    assert!(encoded.get("standard").unwrap().is_null());
    encoded.as_object_mut().unwrap().remove("standard");
    assert!(serde_json::from_value::<DesignPlan>(encoded).is_err());
    let mut encoded = serde_json::to_value(&baseline).unwrap();
    encoded.as_object_mut().unwrap().remove("standard_review");
    assert!(
        serde_json::from_value::<crate::verification::evidence::VerificationReport>(encoded)
            .is_err()
    );
}

#[test]
fn current_saved_run_retains_advice_and_refuses_old_record_contract() {
    use crate::{
        documents::workspace::Workspace,
        material::lineage,
        planning::plan::{Operation, Source},
    };
    let mut declared = plan(RFC1000_L0_ID, &[("promoter", "ACACAC")]);
    declared.steps = vec![Operation::Overlap {
        fragments: vec![Source::Input(declared.inputs[0].id.clone())],
        overlaps: vec![],
        circular: false,
    }];
    let target = Molecule::dna("Independent target", "ACACAC", Topology::Linear).unwrap();
    declared.target_revision = target.revision();
    let workspace = Workspace::from_documents(
        "advisory fixture",
        vec![declared.inputs[0].molecule.clone(), target],
        0,
    )
    .unwrap();
    let event = |index: usize| {
        let document = &workspace.documents[index];
        lineage::EventBinding {
            document_id: document.id.clone(),
            event_id: document.head_event.clone().unwrap(),
            revision: document.head.clone(),
        }
    };
    let commit = lineage::record_run(
        &workspace,
        &declared,
        event(1),
        vec![lineage::InputEventBinding {
            input_id: declared.inputs[0].id.clone(),
            event: event(0),
        }],
        true,
    )
    .unwrap();
    let record = lineage::load_run(&commit.workspace, &commit.receipt.run_id).unwrap();
    assert_eq!(record.schema_version, 2);
    let retained = record.assessment.standard_review.as_ref().unwrap();
    retained
        .validate_binding(
            &record.plan_revision,
            record.plan.standard.as_ref().unwrap(),
        )
        .unwrap();
    assert_eq!(native(retained).assessment.outcome, Outcome::Conforms);
    let bytes = lineage::save_run_record(&record).unwrap();
    let mut old: serde_json::Value = serde_json::from_slice(&bytes).unwrap();
    old["schema_version"] = serde_json::json!(1);
    assert!(lineage::load_run_record(&serde_json::to_vec(&old).unwrap()).is_err());
    let mut missing = record.clone();
    missing.assessment.standard_review = None;
    assert!(lineage::save_run_record(&missing).is_err());
}
