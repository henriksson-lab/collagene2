//! Whole cloning primers use literal edited source cores, not a target-selected product.
use collagene_core::{
    Molecule, Sequence, Topology,
    analysis::Primer,
    assembly_edit_intent::{self, Edit},
    assembly_edited_primer_design::{self as design, Request, RoleOrigin},
    assembly_intent, assembly_junction_edit,
    assembly_recipe::{Preparation, RunPreparation},
    assembly_substitution_intent::{self, Substitution},
    cloning_outputs,
    overlap_extension_primer_design::Options,
    overlap_primer_design::Fragment,
    overlap_profiles::ProfileOptions,
    run_lineage::{EventBinding, InputEventBinding},
    run_material::{self, Capsule, Selection},
    verification::{self, DesignPlan, Operation, PlanInput, Source, Status},
    workspace::Workspace,
};
const VECTOR: &str = "CACGTTCGTAGGACCA";
const ORIGINAL: &str = "AGTCCGATGCTACGTAGGACTTCAGACCTGTAGATCGTCA";
const EDITED: &str = "AGTCCGATGCTACGTATTTTCAGACCTGTAGATCGTCA";
const TARGET: &str = "CACGTTCGTAGGACCAAGTCCGATGCTACGTATTTTCAGACCTGTAGATCGTCA";
fn input(id: &str, text: &str) -> PlanInput {
    let molecule = Molecule::dna(id, text, Topology::Linear).unwrap();
    PlanInput {
        id: id.into(),
        revision: molecule.revision(),
        molecule,
    }
}
fn options() -> Options {
    let mut o = Options {
        homology_flank: 8,
        ..Default::default()
    };
    o.primer_options.min_length = 8;
    o.primer_options.max_length = 8;
    o.primer_options.min_gc_fraction = 0.;
    o.primer_options.max_gc_fraction = 1.;
    o.primer_options.max_tm_difference_c = 100.;
    o.primer_options.max_homopolymer = 200;
    o.primer_options.max_three_prime_complementarity = 200;
    o
}
fn fragment(role: &str, start: usize, end: usize, reverse: bool) -> Fragment {
    Fragment {
        role: role.into(),
        input_id: role.into(),
        start,
        end,
        reverse,
    }
}
fn request() -> Request {
    Request {
        intent: assembly_edit_intent::Request {
            base: assembly_substitution_intent::Request {
                base: assembly_intent::Request {
                    name: "Independent vector and literal edited insert".into(),
                    fragments: vec![
                        fragment("vector", 0, 16, false),
                        fragment("insert", 0, 40, false),
                    ],
                    circular: false,
                },
                junction_edits: vec![],
                substitutions: vec![],
            },
            edits: vec![Edit {
                role: "insert".into(),
                start: 16,
                end: 20,
                expected: "GGAC".into(),
                replacement: "TT".into(),
            }],
        },
        overlap: 8,
        options: options(),
        profile_options: ProfileOptions::default(),
        origins: vec![],
    }
}
fn inputs() -> Vec<PlanInput> {
    vec![input("vector", VECTOR), input("insert", ORIGINAL)]
}
fn check(report: &design::Report, expected: &str) {
    report.validate_binding(&report.context_revision).unwrap();
    assert_eq!(report.intent.target.sequence.as_str(), expected);
    assert_eq!(
        report.preview.verification.status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        report.preview.verification.findings
    );
    assert_eq!(
        report
            .preview
            .verification
            .products
            .last()
            .unwrap()
            .sequence
            .as_str(),
        expected
    );
}
#[test]
fn mixed_exact_and_edited_rows_prepare_actual_branches_before_whole_independent_target_check() {
    let inputs = inputs();
    let before = serde_json::to_value(&inputs).unwrap();
    let req = request();
    let p = design::design(&inputs, &req).unwrap();
    check(&p, TARGET);
    assert_eq!(
        p.context_revision,
        design::design_revision(&inputs, &req).unwrap()
    );
    assert_eq!(
        p.preview
            .verification
            .products
            .iter()
            .map(|m| m.sequence.len())
            .collect::<Vec<_>>(),
        vec![24, 26, 30, 38, 38, 54]
    );
    assert_eq!(p.rows[0].designs.len(), 1);
    assert_eq!(p.rows[1].designs.len(), 3);
    assert!(p.rows[0].edits.is_empty());
    assert_eq!(p.rows[1].edits.len(), 1);
    assert!(
        matches!(&p.recipe.fragments[0].preparation,Preparation::ExactPcr{reverse,..}if reverse.tail=="ATCGGACT")
    );
    assert!(
        matches!(&p.recipe.fragments[1].preparation,Preparation::OverlapExtensionPcr{request}if request.outer.forward.tail.is_empty()&&request.outer.reverse.tail.is_empty())
    );
    assert_eq!(p.tails, design::expected_tails(&inputs, &req).unwrap());
    assert_eq!(p.tails[0].homology, "AGTCCGAT");
    assert_eq!(p.tails[0].downstream_role.as_deref(), Some("insert"));
    assert_eq!(serde_json::to_value(&inputs).unwrap(), before);
}
#[test]
fn outgoing_homology_prefix_is_taken_after_the_downstream_literal_indel() {
    let mut req = request();
    req.overlap = 20;
    let tails = design::expected_tails(&inputs(), &req).unwrap();
    assert_eq!(tails[0].homology, "AGTCCGATGCTACGTATTTT");
    assert_eq!(tails[0].reverse_tail, "AAAATACGTAGCATCGGACT");
    let p = design::design(&inputs(), &req).unwrap();
    check(&p, TARGET);
    assert_eq!(p.preview.verification.products[0].sequence.len(), 36);
}
#[test]
fn source_orientation_is_preserved_for_reverse_ordinary_and_reverse_edited_contributions() {
    for (vector_reverse, insert_reverse) in [(true, false), (false, true), (true, true)] {
        let mut ins = inputs();
        let mut req = request();
        for (index, reverse) in [(0, vector_reverse), (1, insert_reverse)] {
            if reverse {
                ins[index].molecule.sequence =
                    ins[index].molecule.sequence.reverse_complement().unwrap();
                ins[index].revision = ins[index].molecule.revision();
                req.intent.base.base.fragments[index].reverse = true;
            }
        }
        let p = design::design(&ins, &req).unwrap();
        check(&p, TARGET);
        assert_eq!(p.recipe.fragments[0].reverse, vector_reverse);
        assert_eq!(p.recipe.fragments[1].reverse, insert_reverse);
        if vector_reverse {
            assert_eq!(p.tails[0].forward_tail, "ATCGGACT");
            assert!(p.tails[0].reverse_tail.is_empty());
        }
        if insert_reverse {
            assert_eq!(
                p.preview.prepared[1].sequence.as_str(),
                "TGACGATCTACAGGTCTGAAAATACGTAGCATCGGACT"
            );
        }
    }
}
#[test]
fn incoming_trims_junction_insert_and_snp_are_normalized_to_retained_core_coordinates() {
    let ins = vec![
        input("vector", &format!("GGGG{VECTOR}TTTT")),
        input("insert", &format!("CCCC{ORIGINAL}AAAA")),
    ];
    let mut req = request();
    req.intent.base.base.fragments[0].end = 24;
    req.intent.base.base.fragments[1].end = 48;
    req.intent.base.junction_edits = vec![assembly_junction_edit::Edit {
        junction: 0,
        remove_left: 4,
        remove_right: 4,
        inserted: "GC".into(),
    }];
    req.intent.edits[0].start = 20;
    req.intent.edits[0].end = 24;
    req.intent.base.substitutions = vec![Substitution {
        role: "insert".into(),
        offset: 27,
        expected: "A".into(),
        replacement: "G".into(),
    }];
    let p = design::design(&ins, &req).unwrap();
    check(
        &p,
        "GGGGCACGTTCGTAGGACCAGCAGTCCGATGCTACGTATTTTCGGACCTGTAGATCGTCAAAAA",
    );
    assert_eq!(p.rows[1].fragment.start, 4);
    assert_eq!(p.rows[1].fragment.end, 48);
    assert_eq!(
        p.rows[1]
            .edits
            .iter()
            .map(|e| (e.start, e.end))
            .collect::<Vec<_>>(),
        vec![(16, 20), (23, 24)]
    );
    assert_eq!(p.tails[0].outgoing_insert, "GC");
    assert_eq!(p.tails[0].reverse_tail, "ATCGGACTGC");
    assert_eq!(
        p.preview
            .verification
            .products
            .last()
            .unwrap()
            .sequence
            .len(),
        64
    );
}
#[test]
fn circular_closure_and_its_literal_insert_keep_first_retained_origin_and_uncertainty() {
    let mut req = request();
    req.intent.base.base.circular = true;
    req.intent.base.junction_edits = vec![assembly_junction_edit::Edit {
        junction: 1,
        remove_left: 0,
        remove_right: 0,
        inserted: "TA".into(),
    }];
    let p = design::design(&inputs(), &req).unwrap();
    check(
        &p,
        "CACGTTCGTAGGACCAAGTCCGATGCTACGTATTTTCAGACCTGTAGATCGTCATA",
    );
    assert_eq!(
        p.preview.verification.products.last().unwrap().topology,
        Topology::Circular
    );
    assert_eq!(p.tails[1].homology, "CACGTTCG");
    assert_eq!(p.tails[1].outgoing_insert, "TA");
    assert_eq!(p.tails[1].downstream_role.as_deref(), Some("vector"));
    req.options.max_comparisons = 18;
    let conditional = design::design(&inputs(), &req).unwrap();
    assert_eq!(conditional.preview.verification.status, Status::Incomplete);
    assert_eq!(
        conditional
            .preview
            .verification
            .products
            .last()
            .unwrap()
            .sequence,
        p.intent.target.sequence
    );
}
fn binding(w: &Workspace, i: usize) -> EventBinding {
    let d = &w.documents[i];
    EventBinding {
        document_id: d.id.clone(),
        event_id: d.head_event.clone().unwrap(),
        revision: d.head.clone(),
    }
}
fn saved() -> (Molecule, Capsule) {
    let src = input("original", ORIGINAL);
    let target =
        Molecule::dna("Independent original PCR stock", ORIGINAL, Topology::Linear).unwrap();
    let plan = DesignPlan {
        name: "Original scope unchanged".into(),
        target_revision: target.revision(),
        inputs: vec![src.clone()],
        steps: vec![Operation::Pcr {
            template: Source::Input("original".into()),
            forward: Primer::new("F", "AGTCCGAT", "").unwrap(),
            reverse: Primer::new("R", "TGACGATC", "").unwrap(),
        }],
    };
    let w = Workspace::from_documents("Origin", vec![src.molecule, target], 1).unwrap();
    let commit = cloning_outputs::finalize(
        &w,
        &plan,
        binding(&w, 1),
        vec![InputEventBinding {
            input_id: "original".into(),
            event: binding(&w, 0),
        }],
        &cloning_outputs::Options {
            name: "Saved original stock".into(),
            show_prepared: false,
            save_primers: false,
            allow_incomplete: false,
        },
    )
    .unwrap();
    let event = commit.receipt.construct.clone();
    let m = commit.workspace.revisions[&event.revision].clone();
    let origin = run_material::capsule(
        &commit.workspace,
        &Selection {
            run_id: commit.receipt.run.run_id,
            step_index: 0,
            event,
        },
    )
    .unwrap();
    (m, origin)
}
#[test]
fn per_role_closed_origin_replays_only_in_whole_scope_and_never_retargets_original_run() {
    let (source, origin) = saved();
    let immutable = serde_json::to_value(&origin).unwrap();
    let mut ins = inputs();
    ins[1].molecule = source;
    ins[1].revision = ins[1].molecule.revision();
    let mut req = request();
    req.origins = vec![RoleOrigin {
        role: "insert".into(),
        origin,
    }];
    let p = design::design(&ins, &req).unwrap();
    check(&p, TARGET);
    assert_eq!(p.rows[1].verification_status, Status::Incomplete);
    let Preparation::FromRun {
        origin,
        preparation: RunPreparation::OverlapExtensionPcr { .. },
    } = &p.recipe.fragments[1].preparation
    else {
        panic!("exact declared origin")
    };
    assert_eq!(serde_json::to_value(origin).unwrap(), immutable);
    assert_eq!(origin.target.sequence.as_str(), ORIGINAL);
    let mut wrong = ins;
    wrong[1].molecule.name.push('!');
    wrong[1].revision = wrong[1].molecule.revision();
    assert!(design::preflight(&wrong, &req).is_err());
}
#[test]
fn ambiguous_source_exact_pcr_cannot_be_selected_using_the_independent_whole_target() {
    let mut ins = inputs();
    ins[0] = input("vector", &format!("{VECTOR}{VECTOR}"));
    let p = design::design(&ins, &request()).unwrap();
    assert_eq!(p.intent.target.sequence.as_str(), TARGET);
    assert_eq!(p.preview.verification.status, Status::Incomplete);
    assert!(p.preview.verification.products.is_empty());
    assert!(!p.rows[0].search[0].unique_geometry);
    let mut wrong = Molecule::dna(
        "Wrong independent stock target",
        &format!("{VECTOR}{ORIGINAL}"),
        Topology::Linear,
    )
    .unwrap();
    wrong.features.clear();
    let mut plan = design::design(&inputs(), &request())
        .unwrap()
        .preview
        .compiled
        .plan;
    plan.target_revision = wrong.revision();
    let report = verification::verify(&plan, &wrong);
    assert_eq!(report.status, Status::Invalid);
    assert_eq!(report.products.last().unwrap().sequence.as_str(), TARGET);
}
#[test]
fn aggregate_work_explicit_tail_policy_stale_literals_and_receipt_tampering_are_guarded() {
    let mut req = request();
    req.options.outer_forward_tail = "AC".into();
    assert!(design::preflight(&inputs(), &req).is_err());
    req = request();
    req.options.primer_options.reverse_tail = "AC".into();
    assert!(design::preflight(&inputs(), &req).is_err());
    req = request();
    req.intent.edits[0].expected = "AAAA".into();
    assert!(design::preflight(&inputs(), &req).is_err());
    let stock = input("shared", &"ACGT".repeat(2500));
    let mut many = request();
    many.intent.base.base.fragments = vec![Fragment {
        role: "r0".into(),
        input_id: "shared".into(),
        start: 0,
        end: 10000,
        reverse: false,
    }];
    many.intent.edits = vec![Edit {
        role: "r0".into(),
        start: 5000,
        end: 5001,
        expected: "A".into(),
        replacement: "T".into(),
    }];
    assert!(design::preflight(&[stock.clone()], &many).is_ok());
    many.intent.base.base.fragments = (0..32)
        .map(|i| Fragment {
            role: format!("r{i}"),
            input_id: "shared".into(),
            start: 0,
            end: 10000,
            reverse: false,
        })
        .collect();
    many.intent.edits = (0..32)
        .map(|i| Edit {
            role: format!("r{i}"),
            start: 5000,
            end: 5001,
            expected: "A".into(),
            replacement: "T".into(),
        })
        .collect();
    assert!(design::preflight(&[stock], &many).is_err());
    let report = design::design(&inputs(), &request()).unwrap();
    let decoded: design::Report =
        serde_json::from_slice(&serde_json::to_vec(&report).unwrap()).unwrap();
    decoded.validate_binding(&report.context_revision).unwrap();
    for kind in 0..5 {
        let mut bad = report.clone();
        match kind {
            0 => bad.tails[0].homology = "AAAAAAAA".into(),
            1 => bad.rows[1].edits[0].replacement = "GG".into(),
            2 => bad.rows[0].search[0].selected_index = usize::MAX,
            3 => bad.intent.target.sequence = Sequence::dna(ORIGINAL).unwrap(),
            _ => bad.recipe.fragments[1].reverse = true,
        }
        assert!(
            bad.validate_binding(&report.context_revision).is_err(),
            "tamper{kind}"
        );
    }
}
