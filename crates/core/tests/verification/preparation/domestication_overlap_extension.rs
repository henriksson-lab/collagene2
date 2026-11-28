//! Original-source reciprocal PCR, independently checked literal substitutions,
//! and actual physical linear/circular products. No recoded template stock.
use collagene_core::{
    catalog_digest, codon_design, domestication as design, domestication_preparation as physical,
    overlap_extension_primer_design as oe, pdf, plan_oligos,
    run_lineage::{self, EventBinding, InputEventBinding},
    run_material::{self, MaterialOrigin},
    verification::{self, Operation, PlanInput, Source, Status},
    workspace::{self, Workspace},
};
use collagene_model::{Feature, Location, Molecule, Topology};
const ORIGINAL: &str = "AGTCCGATGCTACGTAGAATTCCAGACCTGTAGATCGTCA";
const EDITED: &str = "AGTCCGATGCTACGTAGACTTCCAGACCTGTAGATCGTCA";
fn options() -> oe::Options {
    let mut o = oe::Options {
        homology_flank: 8,
        min_paired_bases: 8,
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
fn fixture(topology: Topology) -> (PlanInput, physical::Request) {
    let mut source = Molecule::dna("Original EcoRI substrate", ORIGINAL, topology).unwrap();
    source.features.push(Feature {
        kind: "misc_feature".into(),
        location: Location::range(8, 12),
        qualifiers: vec![("label".into(), Some("Unchanged GCTA".into()))],
        needs_review: false,
    });
    let mut r = design::Request {
        source_revision: source.revision(),
        selection: catalog_digest::select(
            &catalog_digest::builtin_catalog().unwrap(),
            &["ecori".into()],
        )
        .unwrap(),
        scope: vec![design::SourceRange { start: 0, end: 40 }],
        retained_sites: vec![],
        coding_features: vec![],
        protected: vec![],
        options: codon_design::DesignOptions::default(),
        limits: design::Limits::default(),
    };
    r.source_revision = source.revision();
    let patches = vec![design::Patch {
        position: 18,
        expected: "A".into(),
        replacement: "C".into(),
    }];
    let selection = design::Selection {
        context_revision: design::propose_revision(&source, &r).unwrap(),
        candidate_id: design::candidate_id(&r, &patches).unwrap(),
        patches,
    };
    (
        PlanInput {
            id: physical::INPUT_ID.into(),
            revision: source.revision(),
            molecule: source,
        },
        physical::Request {
            origin: None,
            design: r,
            selection,
            route: physical::Route::OverlapExtension { options: options() },
        },
    )
}
fn repin(input: &mut PlanInput, r: &mut physical::Request) {
    input.revision = input.molecule.revision();
    r.design.source_revision = input.revision.clone();
    r.selection.context_revision = design::propose_revision(&input.molecule, &r.design).unwrap();
    r.selection.candidate_id = design::candidate_id(&r.design, &r.selection.patches).unwrap();
}
fn pin(w: &Workspace, index: usize) -> EventBinding {
    let d = &w.documents[index];
    EventBinding {
        document_id: d.id.clone(),
        event_id: d.head_event.clone().unwrap(),
        revision: d.head.clone(),
    }
}
#[test]
fn whole_linear_and_circular_sources_use_four_or_six_actual_steps_and_all_oligos() {
    for (topology, count) in [(Topology::Linear, 4), (Topology::Circular, 6)] {
        let (input, r) = fixture(topology);
        let original = input.molecule.clone();
        let checked = design::check(&original, &r.design, &r.selection).unwrap();
        assert_eq!(checked.target().sequence.as_str(), EDITED); // Before primer planning.
        let report = physical::materialize(&input, &r, None).unwrap();
        report.validate_binding(&report.context_revision).unwrap();
        assert_eq!(report.target.revision(), checked.target().revision());
        assert!(report.proposal_ready);
        assert_eq!(report.plan.steps.len(), count);
        assert_eq!(report.plan.inputs[0].molecule, original);
        assert_eq!(
            report.verification.status,
            Status::ConsistentUnderAssumptions,
            "{:?}",
            report.verification.findings
        );
        assert_eq!(report.verification.products.len(), count);
        let product = report.verification.products.last().unwrap();
        assert_eq!(product.sequence.as_str(), EDITED);
        assert_eq!(product.topology, topology);
        assert!(verification::plan_target_matches(
            &report.plan,
            product,
            &report.target
        ));
        assert!(
            matches!(&report.plan.steps[0],Operation::Pcr{template:Source::Input(id),..} if id==physical::INPUT_ID)
        );
        assert!(
            matches!(&report.plan.steps[1],Operation::Pcr{template:Source::Input(id),..} if id==physical::INPUT_ID)
        );
        assert!(matches!(
            &report.plan.steps[3],
            Operation::Pcr {
                template: Source::Step(2),
                ..
            }
        ));
        if count == 6 {
            assert!(matches!(
                &report.plan.steps[4],
                Operation::TreatEnds {
                    template: Source::Step(3),
                    treatment: collagene_core::chemistry::EndTreatment::Phosphorylate
                }
            ));
            assert!(
                matches!(&report.plan.steps[5],Operation::Ligate{fragments,circular:true,allow_nicked:false,..} if fragments.len()==1 && matches!(fragments[0].source,Source::Step(4)))
            );
        }
        let physical::MethodEvidence::OverlapExtension {
            preparation,
            designs,
            search,
            ..
        } = &report.method
        else {
            panic!("OE method evidence");
        };
        assert_eq!(report.request.selection.patches.len(), 1);
        assert_eq!(report.request.selection.patches[0].expected, "A");
        assert_eq!(report.request.selection.patches[0].replacement, "C");
        assert_eq!(designs.len(), 3);
        assert!(search.iter().all(|s| s.unique_geometry));
        assert!(
            preparation.outer.forward.tail.is_empty() && preparation.outer.reverse.tail.is_empty()
        );
        let orders = plan_oligos::build(&report.plan, &report.target, &Default::default()).unwrap();
        assert_eq!(orders.uses.len(), 6);
        assert_eq!(
            orders.oligos.records.len(),
            4,
            "shared outer endpoint oligos deduplicate exactly"
        );
        assert!(
            orders
                .uses
                .iter()
                .all(|u| [0, 1, 3].contains(&u.step_index))
        );
        let feature = product
            .features
            .iter()
            .find(|f| f.label() == "Unchanged GCTA")
            .unwrap();
        assert_eq!(feature.location.extract(&product.sequence).unwrap(), "GCTA");
        assert_eq!(input.molecule, original);
    }
}
#[test]
fn synonym_phase_protected_bases_and_wanted_sites_are_checked_before_preparation() {
    let (mut input, mut r) = fixture(Topology::Linear);
    input.molecule.features.push(Feature {
        kind: "CDS".into(),
        location: Location::range(16, 22),
        qualifiers: vec![
            ("codon_start".into(), Some("1".into())),
            ("label".into(), Some("Literal EF".into())),
        ],
        needs_review: false,
    });
    r.selection.patches[0].replacement = "G".into(); // GAA→GAG, unchanged Glu.
    repin(&mut input, &mut r);
    r.design.coding_features = vec![design::feature_pin(&input.molecule, 1).unwrap()];
    repin(&mut input, &mut r);
    let report = physical::materialize(&input, &r, None).unwrap();
    assert_eq!(report.check.features[0].before_protein, "EF");
    assert_eq!(report.check.features[0].after_protein, "EF");
    assert_eq!(&report.target.sequence.as_str()[16..22], "GAGTTC");
    let mut nonsyn = r.clone();
    nonsyn.selection.patches[0].replacement = "C".into();
    repin(&mut input, &mut nonsyn);
    assert!(physical::materialize(&input, &nonsyn, None).is_err());
    let mut protected = r.clone();
    protected.design.protected = vec![design::SourceRange { start: 18, end: 19 }];
    repin(&mut input, &mut protected);
    assert!(physical::materialize(&input, &protected, None).is_err());
    let mut wanted = r.clone();
    wanted.design.retained_sites = vec![design::SitePin {
        enzyme_id: "ecori".into(),
        start: 16,
        reverse: false,
    }];
    repin(&mut input, &mut wanted);
    assert!(physical::materialize(&input, &wanted, None).is_err());
}
#[test]
fn source_method_target_and_pair_tampering_do_not_become_checked_preparations() {
    let (input, r) = fixture(Topology::Linear);
    let report = physical::materialize(&input, &r, None).unwrap();
    let context = report.context_revision.clone();
    let mut changed = report.clone();
    changed.target = Molecule::dna("Unrelated target", ORIGINAL, Topology::Linear).unwrap();
    assert!(changed.validate_binding(&context).is_err());
    let mut changed = report.clone();
    if let physical::MethodEvidence::OverlapExtension { preparation, .. } = &mut changed.method {
        preparation.left.reverse.tail = "AAAA".into();
    }
    assert!(changed.validate_binding(&context).is_err());
    let mut changed = report.clone();
    changed.request.selection.patches[0].replacement = "T".into();
    assert!(changed.validate_binding(&context).is_err());
    let mut changed = report.clone();
    if let Operation::Pcr { template, .. } = &mut changed.plan.steps[1] {
        *template = Source::Step(0);
    }
    changed.verification.plan_revision = verification::plan_revision(&changed.plan);
    assert!(changed.validate_binding(&context).is_err());
    let mut changed = report.clone();
    if let physical::MethodEvidence::OverlapExtension { designs, .. } = &mut changed.method {
        designs[0].options.forward_tail = "TT".into();
    }
    assert!(changed.validate_binding(&context).is_err());
    let mut stale = input.clone();
    stale.revision = "a".repeat(64);
    assert!(physical::preflight(&stale, &r, None).is_err());
}
#[test]
fn saved_native_history_and_protocols_keep_original_source_and_independent_target() {
    for topology in [Topology::Linear, Topology::Circular] {
        let (input, r) = fixture(topology);
        let report = physical::materialize(&input, &r, None).unwrap();
        let w = Workspace::from_documents(
            "Original plus independent frozen intent",
            vec![input.molecule.clone(), report.target.clone()],
            1,
        )
        .unwrap();
        let heads = w
            .documents
            .iter()
            .map(|d| (d.id.clone(), d.head.clone(), d.head_event.clone()))
            .collect::<Vec<_>>();
        let html = verification::export_protocol_html(&report.plan, &report.target);
        assert!(html.contains("overlap") || html.contains("Overlap"));
        let pdf = pdf::export_protocol_pdf(&report.plan, &report.target).unwrap();
        assert!(pdf.bytes.starts_with(b"%PDF-"));
        assert_eq!(pdf.status, report.verification.status);
        let saved = run_lineage::record_run(
            &w,
            &report.plan,
            pin(&w, 1),
            vec![InputEventBinding {
                input_id: physical::INPUT_ID.into(),
                event: pin(&w, 0),
            }],
            false,
        )
        .unwrap();
        let w = workspace::load_workspace(&workspace::save_workspace(&saved.workspace).unwrap())
            .unwrap();
        for (id, head, event) in heads {
            let d = w.documents.iter().find(|d| d.id == id).unwrap();
            assert_eq!(d.head, head);
            assert_eq!(d.head_event, event);
        }
        let replay = run_lineage::replay_run(&w, &saved.receipt.run_id).unwrap();
        assert!(replay.products_match && replay.engine_matches && replay.status_matches);
        let p = saved.receipt.products.last().unwrap();
        let event = EventBinding {
            document_id: p.document_id.clone(),
            event_id: p.event_id.clone(),
            revision: p.revision.clone(),
        };
        let origin = MaterialOrigin::Run(
            run_material::capsule(
                &w,
                &run_material::Selection {
                    run_id: saved.receipt.run_id.clone(),
                    step_index: p.step_index,
                    event,
                },
            )
            .unwrap(),
        );
        let fresh = origin.replay(&w.revisions[&p.revision]).unwrap();
        assert_eq!(fresh.reference().sequence.as_str(), EDITED);
        assert!(fresh.tracked().duplex.nicks.is_empty());
        assert_eq!(
            collagene_core::modified_duplex::summary(fresh.tracked())
                .unwrap()
                .known_modified,
            0
        );
    }
}
#[test]
fn capped_terminal_search_remains_incomplete_without_target_fitting() {
    let (input, mut r) = fixture(Topology::Linear);
    let physical::Route::OverlapExtension { options } = &mut r.route else {
        unreachable!()
    };
    options.max_comparisons = 17;
    let report = physical::materialize(&input, &r, None).unwrap();
    assert_eq!(report.target.sequence.as_str(), EDITED);
    assert_eq!(
        report.verification.status,
        Status::Incomplete,
        "{:?}",
        report.verification.findings
    );
    assert_eq!(report.plan.inputs[0].molecule.sequence.as_str(), ORIGINAL);
    report.validate_binding(&report.context_revision).unwrap();
}
#[test]
fn competing_native_geometries_remain_an_explicit_nonready_diagnostic_proposal() {
    let (mut input, mut r) = fixture(Topology::Linear);
    input.molecule = Molecule::dna(
        "Repeated actual left endpoint",
        "AGTCCGATAGTCCGATGCTACGTAGAATTCCAGACCTGTAGATCGTCA",
        Topology::Linear,
    )
    .unwrap();
    r.design.scope[0].end = input.molecule.sequence.len();
    r.selection.patches[0].position = 26;
    repin(&mut input, &mut r);
    let report = physical::materialize(&input, &r, None).unwrap();
    assert!(!report.proposal_ready);
    assert!(
        report
            .proposal_findings
            .iter()
            .any(|s| s.contains("Diagnostic proposal only"))
    );
    let physical::MethodEvidence::OverlapExtension { search, .. } = &report.method else {
        panic!("OE evidence");
    };
    assert!(!search[0].unique_geometry);
    assert_ne!(
        report.verification.status,
        Status::ConsistentUnderAssumptions
    );
    assert!(report.verification.products.is_empty());
    assert_eq!(&report.target.sequence.as_str()[24..30], "GACTTC");
    report.validate_binding(&report.context_revision).unwrap();
    let mut forged = report.clone();
    forged.proposal_ready = true;
    assert!(forged.validate_binding(&report.context_revision).is_err());
}
#[test]
fn native_options_literals_markers_and_boundary_limits_fail_without_a_recoded_stock() {
    let (input, r) = fixture(Topology::Linear);
    for kind in 0..6 {
        let mut changed = r.clone();
        let physical::Route::OverlapExtension { options } = &mut changed.route else {
            unreachable!()
        };
        match kind {
            0 => options.outer_forward_tail = "AC".into(),
            1 => options.max_candidates = 0,
            2 => options.max_oligo_length = 8,
            3 => options.primer_options.max_pairs = 21,
            4 => options.primer_options.forward_tail = "AC".into(),
            _ => options.max_comparisons = 0,
        };
        assert!(
            physical::preflight(&input, &changed, None).is_err(),
            "invalid options {kind}"
        );
    }
    let mut wrong = r.clone();
    wrong.selection.patches[0].expected = "G".into();
    assert!(physical::materialize(&input, &wrong, None).is_err());
    let mut marked = input.clone();
    marked
        .molecule
        .metadata
        .insert("duplex_model".into(), "untrusted source summary".into());
    let mut r2 = r.clone();
    repin(&mut marked, &mut r2);
    assert!(physical::preflight(&marked, &r2, None).is_err());
    let mut boundary = r.clone();
    boundary.selection.patches[0] = design::Patch {
        position: 0,
        expected: "A".into(),
        replacement: "C".into(),
    };
    let mut same = input.clone();
    repin(&mut same, &mut boundary);
    assert!(physical::preflight(&same, &boundary, None).is_err());
}
