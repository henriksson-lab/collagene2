use collagene_core::{
    chemistry::EndTreatment,
    cloning::{
        assembly::{
            diagnostics::{self as d, JunctionState, RowState, Stage},
            recipe::{self as a, FragmentSpec, Preparation},
        },
        ligation::recipe as l,
    },
    planning::plan::PlanInput,
    sequence::analysis::Primer,
    verification::evidence::Status,
};
use collagene_model::{Molecule, Topology};
fn input(id: &str, seq: &str) -> PlanInput {
    let molecule = Molecule::dna(id, seq, Topology::Linear).unwrap();
    PlanInput {
        id: id.into(),
        revision: molecule.revision(),
        molecule,
    }
}
fn target(seq: &str) -> Molecule {
    Molecule::dna("independent target", seq, Topology::Linear).unwrap()
}
fn pcr(role: &str, id: &str, f: &str, r: &str) -> FragmentSpec {
    FragmentSpec {
        generation: None,
        role: role.into(),
        input_id: id.into(),
        preparation: Preparation::ExactPcr {
            forward: Primer::new("literal F", f, "").unwrap(),
            reverse: Primer::new("literal R", r, "").unwrap(),
            template_strand: collagene_core::chemistry::Strand::Top,
        },
        reverse: false,
    }
}
fn three() -> (Vec<PlanInput>, Vec<FragmentSpec>, Molecule) {
    (
        vec![
            input("a", "ACGACCCCGGGGTAAC"),
            input("b", "TATGAGCTCCAGGTCA"),
            input("c", "GCTTACCGATGACCTA"),
        ],
        vec![
            pcr("first", "a", "ACGACCCC", "GTTACCCC"),
            pcr("middle", "b", "AAAAAAAA", "CCCCCCCC"),
            pcr("later", "c", "GCTTACCG", "TAGGTCAT"),
        ],
        target("ACGACCCCGGGGTAACTATGAGCTCCAGGTCAGCTTACCGATGACCTA"),
    )
}
#[test]
fn missing_middle_pcr_product_locates_that_row_and_keeps_later_steps_blocked_without_text_parsing()
{
    let (inputs, fragments, target) = three();
    let recipe = a::Recipe {
        name: "Missing middle PCR".into(),
        fragments,
        overlaps: vec![4, 4],
        circular: false,
        options: Default::default(),
    };
    let mut p = a::preview(&inputs, &target, &recipe).unwrap();
    assert_eq!(p.verification.products.len(), 1);
    let report = d::overlap(&p).unwrap();
    assert_eq!(
        report.rows.iter().map(|r| r.state).collect::<Vec<_>>(),
        vec![RowState::Reconstructed, RowState::Failed, RowState::Blocked]
    );
    let failure = report.failure.unwrap();
    assert_eq!(
        (failure.step, failure.row, failure.stage),
        (1, Some(1), Stage::Pcr)
    );
    assert_eq!(
        (report.rows[0].completed_steps, report.rows[0].total_steps),
        (1, 1)
    );
    assert!(
        report
            .junctions
            .iter()
            .all(|j| j.state == JunctionState::HomologyUnavailable)
    );
    p.verification.findings = vec!["Arbitrary presentation text claiming row99 and step100".into()];
    let located = d::overlap(&p).unwrap();
    assert_eq!(located.failure.unwrap().row, Some(1));
    assert_eq!(located.findings, p.verification.findings);
}
fn raw(mismatch: bool, max_comparisons: usize, scan: bool, wrong_target: bool) -> a::Preview {
    let inputs = vec![
        input("left", "AAAACCCC"),
        input("right", if mismatch { "GGGGTTTT" } else { "CCCCGGGG" }),
    ];
    let recipe = a::Recipe {
        name: "Explicit raw reference overlap".into(),
        fragments: vec![
            FragmentSpec {
                generation: None,
                role: "left".into(),
                input_id: "left".into(),
                preparation: Preparation::WholeLinear,
                reverse: false,
            },
            FragmentSpec {
                generation: None,
                role: "right".into(),
                input_id: "right".into(),
                preparation: Preparation::WholeLinear,
                reverse: false,
            },
        ],
        overlaps: vec![4],
        circular: false,
        options: collagene_core::cloning::overlap::profiles::ProfileOptions {
            scan_repeated_overlaps: scan,
            max_comparisons,
            ..Default::default()
        },
    };
    a::preview(
        &inputs,
        &target(if mismatch {
            "AAAACCCCGGGGTTTT"
        } else if wrong_target {
            "AAAACCCCGGGA"
        } else {
            "AAAACCCCGGGG"
        }),
        &recipe,
    )
    .unwrap()
}
#[test]
fn local_typed_profile_mismatch_keeps_both_prepared_references_and_does_not_guess_reaction_cause() {
    let p = raw(true, 10000, false, false);
    assert!(p.preparation_complete);
    assert_eq!(p.prepared.len(), 2);
    assert!(p.verification.products.is_empty());
    let report = d::overlap(&p).unwrap();
    assert!(
        report
            .rows
            .iter()
            .all(|r| r.state == RowState::InputReference)
    );
    let failure = report.failure.unwrap();
    assert_eq!(failure.stage, Stage::OverlapAssembly);
    assert_eq!(failure.row, None);
    assert_eq!(failure.upstream, None);
    assert_eq!(failure.downstream, None);
    assert_eq!(report.junctions[0].state, JunctionState::HomologyMismatch);
    assert_eq!(
        (report.junctions[0].upstream, report.junctions[0].downstream),
        (0, 1)
    );
}
#[test]
fn complete_products_do_not_turn_incomplete_evidence_or_wrong_target_into_local_preparation_failures()
 {
    let consistent = d::overlap(&raw(false, 10000, false, false)).unwrap();
    assert_eq!(consistent.status, Status::ConsistentUnderAssumptions);
    assert!(consistent.failure.is_none());
    assert_eq!(consistent.junctions[0].state, JunctionState::ExactHomology);
    let p = raw(false, 4, true, false);
    assert_eq!(p.verification.products.len(), p.compiled.plan.steps.len());
    let incomplete = d::overlap(&p).unwrap();
    assert_eq!(incomplete.status, Status::Incomplete);
    assert!(incomplete.failure.is_none());
    assert!(!incomplete.warnings.is_empty());
    let mismatch = d::overlap(&raw(false, 10000, false, true)).unwrap();
    assert_eq!(mismatch.status, Status::Invalid);
    assert_eq!(mismatch.failure.unwrap().stage, Stage::TargetComparison);
    assert_eq!(mismatch.junctions[0].state, JunctionState::ExactHomology);
}
#[test]
fn ligation_treatments_are_after_all_preparations_and_missing_assembly_never_names_an_invented_pair()
 {
    let (inputs, rows, t) = three();
    let recipe = l::Recipe {
        name: "Treatments after preparations".into(),
        fragments: rows
            .into_iter()
            .map(|source| l::Fragment {
                source,
                treatment: Some(EndTreatment::Phosphorylate),
            })
            .collect(),
        circular: false,
        allow_nicked: false,
    };
    let p = l::preview(
        &inputs,
        &t,
        &recipe,
        collagene_core::cloning::assembly::recipe::CleavagePolicy::AssumedComplete,
    )
    .unwrap();
    assert_eq!(p.verification.products.len(), 1);
    let report = d::ligation(&p).unwrap();
    assert_eq!(
        report.rows.iter().map(|r| r.state).collect::<Vec<_>>(),
        vec![RowState::Blocked, RowState::Failed, RowState::Blocked]
    );
    assert_eq!(
        (report.rows[0].completed_steps, report.rows[0].total_steps),
        (1, 2)
    );
    assert_eq!(
        (
            report.failure.as_ref().unwrap().step,
            report.failure.as_ref().unwrap().row
        ),
        (1, Some(1))
    );
    let inputs = inputs[..2].to_vec();
    let mut recipe = l::Recipe {
        name: "No implicit phosphorylation".into(),
        fragments: vec![
            l::Fragment {
                source: pcr("left", "a", "ACGACCCC", "GTTACCCC"),
                treatment: None,
            },
            l::Fragment {
                source: pcr("right", "b", "TATGAGCT", "TGACCTGG"),
                treatment: None,
            },
        ],
        circular: false,
        allow_nicked: false,
    };
    let p = l::preview(
        &inputs,
        &target("ACGACCCCGGGGTAACTATGAGCTCCAGGTCA"),
        &recipe,
        collagene_core::cloning::assembly::recipe::CleavagePolicy::AssumedComplete,
    )
    .unwrap();
    assert!(p.preparation_complete);
    let report = d::ligation(&p).unwrap();
    assert!(
        report
            .rows
            .iter()
            .all(|r| r.state == RowState::Reconstructed)
    );
    let f = report.failure.unwrap();
    assert_eq!(f.stage, Stage::LigationAssembly);
    assert_eq!((f.row, f.upstream, f.downstream), (None, None, None));
    assert_eq!(report.junctions[0].state, JunctionState::EndpointsAvailable);
    for row in &mut recipe.fragments {
        row.treatment = Some(EndTreatment::Phosphorylate);
    }
    let p = l::preview(
        &inputs,
        &target("ACGACCCCGGGGTAACTATGAGCTCCAGGTCA"),
        &recipe,
        collagene_core::cloning::assembly::recipe::CleavagePolicy::AssumedComplete,
    )
    .unwrap();
    assert_eq!(
        d::ligation(&p).unwrap().status,
        Status::ConsistentUnderAssumptions
    );
}
#[test]
fn forged_source_prepared_mapping_engine_and_context_shape_fail_closed() {
    let (inputs, fragments, t) = three();
    let recipe = a::Recipe {
        name: "Forged receipt guards".into(),
        fragments,
        overlaps: vec![4, 4],
        circular: false,
        options: Default::default(),
    };
    let p = a::preview(&inputs, &t, &recipe).unwrap();
    let mut stale = p.clone();
    stale.compiled.plan.inputs[0].revision = "0".repeat(64);
    assert!(d::overlap(&stale).is_err());
    let mut prepared = p.clone();
    prepared.prepared[0].sequence = collagene_model::Sequence::dna("AAAAAAAAAAAAAAAA").unwrap();
    assert!(d::overlap(&prepared).is_err());
    let mut engine = p.clone();
    engine.verification.engine_version = "archived engine with unknown prefix semantics".into();
    assert!(d::overlap(&engine).is_err());
    let mut context = p;
    context.context_revision = "stale".into();
    assert!(d::overlap(&context).is_err());
    let recipe = l::Recipe {
        name: "Forged native chain".into(),
        fragments: vec![
            l::Fragment {
                source: pcr("left", "a", "ACGACCCC", "GTTACCCC"),
                treatment: None,
            },
            l::Fragment {
                source: pcr("right", "b", "TATGAGCT", "TGACCTGG"),
                treatment: None,
            },
        ],
        circular: false,
        allow_nicked: false,
    };
    let mut p = l::preview(
        &inputs[..2],
        &target("ACGACCCCGGGGTAACTATGAGCTCCAGGTCA"),
        &recipe,
        collagene_core::cloning::assembly::recipe::CleavagePolicy::AssumedComplete,
    )
    .unwrap();
    p.compiled.rows[0].preparation.preparation_step = Some(1);
    assert!(d::ligation(&p).is_err());
}
