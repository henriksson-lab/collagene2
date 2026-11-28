//! Literal graph/slot oracle; no reconstructed-prefix assumptions or filtered slots.
use collagene_core::{
    chemistry::Strand,
    planning::plan::{DesignPlan, Operation, PlanInput, Source},
    sequence::analysis::Primer,
    verification::{
        checks::{CheckFinding, CheckScope, CheckVerdict, MAX_FINDINGS, RuleId},
        engine::verify,
        evidence::{Status, StepState, VerificationReport},
    },
};
use collagene_model::{Molecule, Topology};

const DNA: &str = "AGTCCGATGCTACGTAGGACTTCAGACCTGTA";
fn molecule() -> Molecule {
    Molecule::dna("independent source", DNA, Topology::Linear).unwrap()
}
fn pcr(template: Source) -> Operation {
    Operation::Pcr {
        template,
        forward: Primer::new("literal first eight", "AGTCCGAT", "").unwrap(),
        reverse: Primer::new("literal reverse last eight", "TACAGGTC", "").unwrap(),
        template_strand: Strand::Top,
    }
}
fn plan(steps: Vec<Operation>) -> (DesignPlan, Molecule) {
    let target = molecule();
    let source = molecule();
    (
        DesignPlan {
            standard: None,
            name: "independent branch fixture".into(),
            target_revision: target.revision(),
            inputs: vec![PlanInput {
                id: "stock".into(),
                revision: source.revision(),
                molecule: source,
            }],
            steps,
        },
        target,
    )
}
fn check(report: &VerificationReport, rule: RuleId, verdict: CheckVerdict, scope: CheckScope) {
    assert!(
        report
            .findings
            .iter()
            .any(|f| f.rule == rule && f.verdict == verdict && f.scope == scope),
        "missing {rule:?}/{verdict:?}/{scope:?}: {:?}",
        report.findings
    );
}
fn two_failures() -> (DesignPlan, Molecule) {
    let mut bad_pcr = pcr(Source::Input("stock".into()));
    let Operation::Pcr {
        forward, reverse, ..
    } = &mut bad_pcr
    else {
        unreachable!()
    };
    *forward = Primer::new("absent forward", "CCCCCCCC", "").unwrap();
    *reverse = Primer::new("absent reverse", "CCCCCCCC", "").unwrap();
    plan(vec![
        bad_pcr,
        Operation::Digest {
            template: Source::Input("stock".into()),
            enzymes: vec!["EcoRI".into()],
            fragment_index: 99,
        },
        pcr(Source::Input("stock".into())),
        pcr(Source::Step(0)),
        pcr(Source::Step(1)),
        pcr(Source::Step(2)),
    ])
}

#[test]
fn separate_failures_do_not_compact_later_products_or_erase_target_pass() {
    let (plan, target) = two_failures();
    let report = verify(&plan, &target);
    assert_eq!(report.status, Status::Invalid);
    assert_eq!(report.planned_steps, 6);
    assert_eq!(report.steps.len(), 6);
    assert_eq!(report.reconstructed_count(), 2);
    assert_eq!(
        report
            .steps
            .iter()
            .map(|s| s.step_index)
            .collect::<Vec<_>>(),
        vec![0, 1, 2, 3, 4, 5]
    );
    for index in [0, 1] {
        assert_eq!(report.steps[index].state, StepState::Rejected);
        assert!(report.product(index).is_none());
    }
    for index in [2, 5] {
        assert_eq!(report.steps[index].state, StepState::Reconstructed);
        let actual = report
            .product(index)
            .expect("independent branch retains actual index");
        assert_eq!(actual.sequence.as_str(), DNA);
        assert_eq!(actual.topology, Topology::Linear);
    }
    assert_eq!(report.steps[3].state, StepState::Blocked);
    assert_eq!(report.steps[3].blocked_by, vec![Source::Step(0)]);
    assert_eq!(report.steps[4].state, StepState::Blocked);
    assert_eq!(report.steps[4].blocked_by, vec![Source::Step(1)]);
    check(
        &report,
        RuleId::PcrBinding,
        CheckVerdict::Fail,
        CheckScope::Step { step_index: 0 },
    );
    check(
        &report,
        RuleId::DigestSelection,
        CheckVerdict::Fail,
        CheckScope::Step { step_index: 1 },
    );
    check(
        &report,
        RuleId::DependencyBlocked,
        CheckVerdict::Blocked,
        CheckScope::Step { step_index: 3 },
    );
    check(
        &report,
        RuleId::TargetAgreement,
        CheckVerdict::Pass,
        CheckScope::Target { step_index: 5 },
    );
    assert_eq!(report.final_product().unwrap().sequence.as_str(), DNA);
    report.validate_documentary().unwrap();
}

#[test]
fn wrong_target_binding_keeps_independent_checks_and_products_without_target_agreement() {
    let (mut plan, target) = two_failures();
    plan.target_revision = "0".repeat(64);
    let report = verify(&plan, &target);
    assert_eq!(report.status, Status::Invalid);
    assert_eq!(report.planned_steps, 6);
    assert_eq!(report.steps.len(), 6);
    assert_eq!(report.reconstructed_count(), 2);
    check(
        &report,
        RuleId::TargetBinding,
        CheckVerdict::Fail,
        CheckScope::Plan,
    );
    check(
        &report,
        RuleId::PcrBinding,
        CheckVerdict::Fail,
        CheckScope::Step { step_index: 0 },
    );
    check(
        &report,
        RuleId::DigestSelection,
        CheckVerdict::Fail,
        CheckScope::Step { step_index: 1 },
    );
    for index in [2, 5] {
        assert_eq!(report.steps[index].step_index, index);
        assert_eq!(report.steps[index].state, StepState::Reconstructed);
        let product = report
            .product(index)
            .expect("ordinary branch is independent of target authentication");
        assert_eq!(product.sequence.as_str(), DNA);
        assert_eq!(product.topology, Topology::Linear);
    }
    assert_eq!(report.steps[3].blocked_by, vec![Source::Step(0)]);
    assert_eq!(report.steps[4].blocked_by, vec![Source::Step(1)]);
    check(
        &report,
        RuleId::TargetAgreement,
        CheckVerdict::Blocked,
        CheckScope::Target { step_index: 5 },
    );
    assert!(
        !report.findings.iter().any(|finding| {
            finding.rule == RuleId::TargetAgreement && finding.verdict == CheckVerdict::Pass
        }),
        "matching letters cannot authenticate the independently supplied target"
    );
    report.validate_documentary().unwrap();
}

#[test]
fn stale_input_blocks_only_its_consumers_and_target_mismatch_remains_separate() {
    let (mut plan, mut target) = plan(vec![
        pcr(Source::Input("stale".into())),
        pcr(Source::Input("stock".into())),
    ]);
    let mut stale = plan.inputs[0].clone();
    stale.id = "stale".into();
    stale.revision = "0".repeat(64);
    plan.inputs.push(stale);
    target.sequence = collagene_model::Sequence::dna("TTTTTTTT").unwrap();
    plan.target_revision = target.revision();
    let report = verify(&plan, &target);
    assert_eq!(report.status, Status::Invalid);
    assert_eq!(report.steps[0].state, StepState::Blocked);
    assert_eq!(
        report.steps[0].blocked_by,
        vec![Source::Input("stale".into())]
    );
    assert!(report.product(0).is_none());
    assert_eq!(report.product(1).unwrap().sequence.as_str(), DNA);
    check(
        &report,
        RuleId::InputRevision,
        CheckVerdict::Fail,
        CheckScope::Input {
            input_id: "stale".into(),
        },
    );
    check(
        &report,
        RuleId::TargetAgreement,
        CheckVerdict::Fail,
        CheckScope::Target { step_index: 1 },
    );
    report.validate_documentary().unwrap();
}

#[test]
fn final_product_is_designated_last_slot_and_never_last_surviving_branch() {
    let (mut plan, target) = two_failures();
    plan.steps.truncate(4);
    let report = verify(&plan, &target);
    assert_eq!(report.status, Status::Invalid);
    assert_eq!(report.product(2).unwrap().sequence.as_str(), DNA);
    assert!(report.product(3).is_none());
    assert!(report.final_product().is_none());
}

#[test]
fn retained_slots_blockers_status_and_typed_checks_reject_documentary_tampering() {
    let (plan, target) = two_failures();
    let report = verify(&plan, &target);
    report.validate_documentary().unwrap();
    let mut missing = report.clone();
    missing.steps.pop();
    assert!(missing.validate_documentary().is_err());
    let mut duplicate = report.clone();
    duplicate.steps[2].step_index = 0;
    assert!(duplicate.validate_documentary().is_err());
    let mut ghost = report.clone();
    ghost.steps[0].product = report.product(2).cloned();
    assert!(ghost.validate_documentary().is_err());
    let mut fake_state = report.clone();
    fake_state.steps[2].state = StepState::Rejected;
    assert!(fake_state.validate_documentary().is_err());
    let mut future = report.clone();
    future.steps[3].blocked_by = vec![Source::Step(4)];
    assert!(future.validate_documentary().is_err());
    let mut successful_blocker = report.clone();
    successful_blocker.steps[3].blocked_by = vec![Source::Step(2)];
    assert!(successful_blocker.validate_documentary().is_err());
    let mut repeated = report.clone();
    repeated.steps[3].blocked_by = vec![Source::Step(0), Source::Step(0)];
    assert!(repeated.validate_documentary().is_err());
    let mut changed_status = report.clone();
    changed_status.status = Status::ConsistentUnderAssumptions;
    assert!(changed_status.validate_documentary().is_err());
    let mut absent_scope = report.clone();
    absent_scope.findings.push(CheckFinding::new(
        RuleId::PcrBinding,
        CheckVerdict::Pass,
        CheckScope::Step { step_index: 99 },
        "unassociated claim",
    ));
    assert!(absent_scope.validate_documentary().is_err());
    let mut unbounded = report;
    unbounded.findings = vec![
        CheckFinding::new(
            RuleId::PlanShape,
            CheckVerdict::Informational,
            CheckScope::Plan,
            "bounded fixture"
        );
        MAX_FINDINGS + 1
    ];
    assert!(unbounded.validate_documentary().is_err());
}

#[test]
fn named_predicates_require_their_own_rule_citation_and_bounded_details() {
    let finding = CheckFinding::new(
        RuleId::PcrBinding,
        CheckVerdict::Fail,
        CheckScope::Step { step_index: 0 },
        "No literal primer pair",
    );
    finding.validate().unwrap();
    let mut erased = finding.clone();
    erased.citations.clear();
    assert!(erased.validate().is_err());
    let mut relabeled = finding.clone();
    relabeled.rule = RuleId::TargetAgreement;
    assert!(relabeled.validate().is_err());
    let mut oversized = finding;
    oversized.message = "x".repeat(4097);
    assert!(oversized.validate().is_err());
}

#[test]
fn directly_owned_deep_finding_and_product_trees_are_refused_before_recursive_counting() {
    use collagene_model::{Error, Feature, Location};
    fn deep_location() -> Location {
        // Deliberately exceeds the canonical32-node depth, while this fixed100
        // depth remains safe to drop even if an assertion fails.
        let mut location = Location::range(0, 1);
        for _ in 0..100 {
            location = Location::Complement(Box::new(location));
        }
        location
    }
    let finding = CheckFinding::new(
        RuleId::AnnotationProjection,
        CheckVerdict::Informational,
        CheckScope::Plan,
        "Native caller supplied a deep standalone annotation coordinate",
    )
    .with_coordinates(vec![deep_location()]);
    assert!(matches!(finding.validate(), Err(Error::TooLarge)));

    let (plan, target) = plan(vec![pcr(Source::Input("stock".into()))]);
    let mut report = verify(&plan, &target);
    report.validate_documentary().unwrap();
    report.findings.push(finding);
    assert!(matches!(
        report.validate_documentary(),
        Err(Error::TooLarge)
    ));
    report.findings.pop();
    report.steps[0]
        .product
        .as_mut()
        .unwrap()
        .features
        .push(Feature {
            kind: "misc_feature".into(),
            location: deep_location(),
            qualifiers: vec![],
            needs_review: false,
        });
    assert!(matches!(
        report.validate_documentary(),
        Err(Error::TooLarge)
    ));
}

#[test]
fn failed_binding_searches_keep_one_shared_work_ledger() {
    use collagene_core::sequence::analysis::{
        ExactSearchBudget, TemplateView, simulate_pcr_geometry_view_with_budget,
    };
    let source = molecule();
    let forward = Primer::new("F", "CCCCCCCC", "").unwrap();
    let reverse = Primer::new("R", "CCCCCCCC", "").unwrap();
    let mut budget = ExactSearchBudget::new(1024).unwrap();
    let first = simulate_pcr_geometry_view_with_budget(
        TemplateView::reference(&source),
        &forward,
        &reverse,
        &mut budget,
    )
    .unwrap();
    assert!(first.is_empty());
    assert!(
        budget.used() > 0,
        "zero-hit searches still consume actual work"
    );
    let mut exhausted = false;
    for _ in 0..64 {
        let previous = budget.used();
        let result = simulate_pcr_geometry_view_with_budget(
            TemplateView::reference(&source),
            &forward,
            &reverse,
            &mut budget,
        );
        assert!(budget.used() >= previous);
        if result.is_err() {
            exhausted = true;
            break;
        }
        assert!(result.unwrap().is_empty());
    }
    assert!(
        exhausted,
        "repeated operations must not reset the shared ledger"
    );
    assert_eq!(budget.used(), 1024);
    assert_eq!(budget.remaining(), 0);
    assert!(
        simulate_pcr_geometry_view_with_budget(
            TemplateView::reference(&source),
            &forward,
            &reverse,
            &mut budget,
        )
        .is_err()
    );
    assert_eq!(budget.used(), 1024);
}
