#[path = "../../common/material_library.rs"]
mod fixture;
use collagene_core::{
    cloning::assembly::intent::{ContributionPosition, Selector},
    documents::{history::DocumentHistory, workspace},
    material::lineage as run_lineage,
    planning::combinatorial as library,
    verification::{self, evidence::Status},
};
use collagene_model::{Molecule, Topology};
use fixture::*;
#[test]
fn actual_bottom_library_freezes_independent_target_and_releases_replayable_circle_after_restart() {
    let (w, d, row) = fixture();
    let before = workspace::save_workspace(&w).unwrap();
    let resolved = library::resolve(&d, &row, &w).unwrap();
    assert_eq!(resolved.target.sequence.as_str(), EXPECTED);
    let Some(library::IndependentIntent::Regions(intent)) = resolved.intent else {
        panic!("native actual-strand declaration")
    };
    assert_eq!(intent.fragments[0].original_sequence.as_str(), EXPECTED);
    assert_eq!(
        intent.fragments[0].source_position(0).unwrap(),
        ContributionPosition::Material {
            strand: collagene_core::chemistry::Strand::Bottom,
            offset: 0
        }
    );
    let mapping = &intent.fragments[0].top_mappings[0];
    assert_eq!(
        (
            mapping.source_start,
            mapping.source_end,
            mapping.destination_start
        ),
        (15, 21, 4)
    );
    let evaluated = library::evaluate(&d, &row, &w).unwrap();
    assert_eq!(
        evaluated.verification().status,
        Status::ConsistentUnderAssumptions
    );
    assert_eq!(
        evaluated
            .verification()
            .products
            .last()
            .unwrap()
            .sequence
            .as_str(),
        EXPECTED
    );
    assert_eq!(workspace::save_workspace(&w).unwrap(), before);
    let (committed, pin, run) = released();
    let restored =
        workspace::load_workspace(&workspace::save_workspace(&committed.workspace).unwrap())
            .unwrap();
    assert_eq!(
        restored.revisions[&pin.revision].sequence.as_str(),
        EXPECTED
    );
    let replay = run_lineage::replay_run(&restored, &run).unwrap();
    assert!(replay.products_match);
    assert_eq!(
        replay.current_report.status,
        Status::ConsistentUnderAssumptions
    );
    let retained = run_lineage::load_run(&restored, &run).unwrap();
    let html = verification::output::export_protocol_html(
        &retained.plan,
        &restored.revisions[&retained.target.revision],
    );
    assert!(html.contains(EXPECTED));
    assert!(
        collagene_core::output::pdf::export_protocol_pdf(
            &retained.plan,
            &restored.revisions[&retained.target.revision]
        )
        .unwrap()
        .bytes
        .starts_with(b"%PDF")
    );
}
#[test]
fn wrong_independent_target_state_source_alias_and_undeclared_or_different_origin_refuse() {
    let (mut w, d, mut row) = fixture();
    let wrong = Molecule::dna("Explicit wrong target", "AATTCGTAAA", Topology::Circular).unwrap();
    DocumentHistory::new("independent-wrong", wrong)
        .unwrap()
        .contribute(&mut w)
        .unwrap();
    row.target = library::TargetDeclaration::Existing(event(&w, w.documents.len() - 1));
    assert_eq!(
        library::evaluate(&d, &row, &w)
            .unwrap()
            .verification()
            .status,
        Status::Invalid
    );
    let (w, mut d, row) = fixture();
    let Selector::MaterialStrandRange {
        expected_state_revision,
        ..
    } = &mut d.bins[0].choices[0].selector
    else {
        unreachable!()
    };
    *expected_state_revision = "0".repeat(64);
    assert!(library::resolve(&d, &row, &w).is_err());
    let (w, mut d, row) = fixture();
    d.material_bindings.clear();
    assert!(library::resolve(&d, &row, &w).is_err());
    let (w, mut d, row) = fixture();
    d.bins[0].choices[0].source = event(&w, 2);
    assert!(
        library::resolve(&d, &row, &w).is_err(),
        "equal top letters are not a native event alias"
    );
    let (w, mut d, row) = fixture();
    d.material_bindings[0].origin = match d.material_bindings[0].origin.clone() {
        collagene_core::material::run::MaterialOrigin::Run(mut o) => {
            o.target.name.push('!');
            collagene_core::material::run::MaterialOrigin::Run(o)
        }
        _ => unreachable!(),
    };
    assert!(library::resolve(&d, &row, &w).is_err());
}

#[test]
fn actual_bottom_pcr_work_counts_both_seventeen_base_frame_and_real_coverage_scan() {
    use collagene_core::{
        chemistry::Strand,
        cloning::assembly::recipe::{Preparation, RunPreparation},
        material::run::MaterialOrigin,
        planning::combinatorial::batch,
        sequence::analysis::Primer,
    };
    let (w, pin, origin) = native_source("AGTCCGATGCTAC", "GTAGGTAGCATCGGACT", 13);
    let source = &w.revisions[&pin.revision];
    assert_eq!(source.sequence.len(), 13);
    let fresh = origin.replay(source).unwrap();
    assert_eq!(fresh.tracked().duplex.bottom.len(), 17);
    let state = fresh.tracked().summary().unwrap().state_revision;
    let MaterialOrigin::Run(capsule) = origin.clone() else {
        unreachable!()
    };
    let d = library::Design {
        id: "strand-cost".into(),
        name: "Real Bottom PCR cost".into(),
        circular: false,
        method: library::Method::Ligation {
            allow_nicked: false,
        },
        material_bindings: vec![library::LibraryMaterialBinding {
            event: pin.clone(),
            origin,
        }],
        bins: vec![library::Bin {
            id: "insert".into(),
            label: "Real native strand".into(),
            choices: vec![library::Choice {
                anneal_bottom: None,
                id: "bottom-seventeen".into(),
                label: "Bottom0..17".into(),
                source: pin,
                selector: Selector::MaterialStrandRange {
                    strand: Strand::Bottom,
                    start: 0,
                    end: 17,
                    expected_state_revision: state,
                },
                preparation_reverse: false,
                generation: None,
                treatment: None,
                preparation: Some(Preparation::FromRun {
                    origin: capsule,
                    preparation: RunPreparation::ExactPcr {
                        template_strand: Strand::Bottom,
                        forward: Primer::new("F", "GTAGGTAG", "").unwrap(),
                        reverse: Primer::new("R", "AGTCCGAT", "").unwrap(),
                    },
                }),
            }],
        }],
    };
    let row = library::SelectedConstruct {
        id: "bottom-pcr".into(),
        name: "Literal17 native PCR".into(),
        choices: vec![library::ChoiceSelection {
            bin_id: "insert".into(),
            choice_id: "bottom-seventeen".into(),
        }],
        target: library::TargetDeclaration::OriginalRegions {
            junction_edits: vec![],
            substitutions: vec![],
            edits: vec![],
        },
        overrides: None,
    };
    let resolved = library::resolve(&d, &row, &w).unwrap();
    assert_eq!(resolved.target.sequence.as_str(), "GTAGGTAGCATCGGACT");
    // Four binding orientations ×16 primer calls ×(17 actual frame +17
    // actual staggered coverage calls), plus13 original source bookkeeping.
    // The obsolete top-only estimate845 cannot pass this independently literal bound.
    assert_eq!(resolved.declared_work(), 2189);
    let mut session = batch::Session::begin(
        w,
        batch::Request {
            design: Some(d),
            rows: vec![batch::RowRequest::Ready(row)],
        },
    )
    .unwrap();
    let ctx = session.summary().context_revision;
    let row = session.advance(&ctx, 0).unwrap().row.unwrap();
    assert_eq!(row.status, Some(Status::ConsistentUnderAssumptions));
    assert_eq!(row.completed_steps, row.planned_steps);
}
