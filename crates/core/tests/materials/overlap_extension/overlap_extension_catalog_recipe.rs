//! Actual OE-PCR, exact catalog digest and optional repair remain one owned graph.
use collagene_core::{
    chemistry::EndRepairMethod,
    cloning::{
        assembly::{
            diagnostics as assembly_diagnostics,
            recipe::{self as a, FragmentSpec, Preparation, RunPreparation},
        },
        ligation::recipe as lig,
        overlap::{
            extension::{
                self as overlap_extension,
                pcr::{PrimerPair, Request},
            },
            profiles::ProfileOptions,
        },
        restriction::catalog as catalog_digest,
    },
    documents::workspace::Workspace,
    material::{
        lineage::{EventBinding, InputEventBinding},
        outputs as cloning_outputs,
        run::Selection,
    },
    planning::plan::{DesignPlan, Operation, PlanInput, Source},
    sequence::analysis::Primer,
    verification::evidence::Status,
};
use collagene_model::{Molecule, Topology};
const ORIGINAL: &str = "AGTCCGATGCTACGTAGGACTTCAGACCTGTAGATCGTCA";
const DIGEST: &str = "AATTCAGTCCGATGCTACGTATTTTCAGACCTGTAGATCGTCAG";
const FILLED: &str = "AATTCAGTCCGATGCTACGTATTTTCAGACCTGTAGATCGTCAGAATT";
fn input(m: Molecule) -> PlanInput {
    PlanInput {
        id: "original".into(),
        revision: m.revision(),
        molecule: m,
    }
}
fn original() -> Molecule {
    Molecule::dna("Original40", ORIGINAL, Topology::Linear).unwrap()
}
fn pair(f: &str, ft: &str, r: &str, rt: &str) -> PrimerPair {
    PrimerPair {
        forward: Primer::new("Forward", f, ft).unwrap(),
        reverse: Primer::new("Reverse", r, rt).unwrap(),
        template_strand: collagene_core::chemistry::Strand::Top,
    }
}
fn request() -> Request {
    Request {
        left: pair("AGTCCGAT", "", "TACGTAGC", "GGTCTGAAAA"),
        right: pair("TTCAGACC", "GCTACGTATT", "TGACGATC", ""),
        outer: pair("AGTCCGAT", "ATGCGAATTC", "TGACGATC", "ATGCGAATTC"),
        extension: overlap_extension::Options {
            overlap: 18,
            min_paired_bases: 8,
            ..Default::default()
        },
    }
}
fn selection() -> catalog_digest::CatalogSelection {
    let c = catalog_digest::builtin_catalog().unwrap();
    catalog_digest::select(&c, &["ecori".into()]).unwrap()
}
fn preparation(repair: Option<EndRepairMethod>) -> Preparation {
    Preparation::OverlapExtensionPcrCatalogDigest {
        conditions: None,
        request: request(),
        selection: selection(),
        fragment_index: 1,
        repair,
    }
}
fn fragment(preparation: Preparation, reverse: bool) -> FragmentSpec {
    FragmentSpec {
        generation: None,
        role: "edited-role".into(),
        input_id: "original".into(),
        preparation,
        reverse,
    }
}
fn linear(preparation: Preparation, reverse: bool) -> a::Recipe {
    a::Recipe {
        name: "Actual edited digest and repair".into(),
        fragments: vec![fragment(preparation, reverse)],
        overlaps: vec![],
        circular: false,
        options: ProfileOptions::default(),
    }
}
#[test]
fn outer_actual_pcr_is_cut_under_exact_catalog_pin_and_closes_real_sticky_ends() {
    let inputs = vec![input(original())];
    let before = serde_json::to_value(&inputs).unwrap();
    let target = Molecule::dna("Independent44bp circle", DIGEST, Topology::Circular).unwrap();
    let recipe = lig::Recipe {
        name: "OE to pinned EcoRI to actual ligation".into(),
        fragments: vec![lig::Fragment {
            source: fragment(preparation(None), false),
            treatment: None,
        }],
        circular: true,
        allow_nicked: false,
    };
    let p = lig::preview(
        &inputs,
        &target,
        &recipe,
        collagene_core::cloning::assembly::recipe::CleavagePolicy::AssumedComplete,
    )
    .unwrap();
    assert_eq!(
        p.verification.status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        p.verification.findings
    );
    assert_eq!(
        p.verification
            .steps
            .iter()
            .map(|step| step
                .product
                .as_ref()
                .expect("each declared step must reconstruct"))
            .map(|m| m.sequence.len())
            .collect::<Vec<_>>(),
        vec![26, 30, 38, 58, 44, 44]
    );
    assert!(matches!(&p.compiled.plan.steps[4],Operation::CatalogDigest{
 conditions: None,template:Source::Step(3),selection:pin,fragment_index:1}if pin==&selection()));
    assert_eq!(p.compiled.rows[0].preparation.preparation_step, Some(0));
    assert!(p.compiled.rows[0].preparation.repair_step.is_none());
    assert!(matches!(
        p.compiled.rows[0].preparation.prepared_source,
        Source::Step(4)
    ));
    assert_eq!(
        p.prepared[0]
            .as_ref()
            .expect("prepared row must reconstruct")
            .sequence
            .as_str(),
        DIGEST
    );
    assert!(
        p.product_physical
            .as_ref()
            .unwrap()
            .covalently_closed_in_model
    );
    assert_eq!(
        assembly_diagnostics::ligation(&p).unwrap().rows[0].total_steps,
        5
    );
    assert_eq!(serde_json::to_value(&inputs).unwrap(), before);
}
#[test]
fn optional_fill_in_uses_the_digest_strands_and_final_orientation_is_separate() {
    for (reverse, letters) in [
        (false, FILLED),
        (true, "AATTCTGACGATCTACAGGTCTGAAAATACGTAGCATCGGACTGAATT"),
    ] {
        let inputs = vec![input(original())];
        let target = Molecule::dna(
            "Independent48bp repaired fragment",
            letters,
            Topology::Linear,
        )
        .unwrap();
        let recipe = linear(
            preparation(Some(EndRepairMethod::FillFivePrimeOverhangs)),
            reverse,
        );
        let p = a::preview(&inputs, &target, &recipe).unwrap();
        assert_eq!(
            p.verification.status,
            Status::ConsistentUnderAssumptions,
            "{:?}",
            p.verification.findings
        );
        assert_eq!(p.compiled.rows[0].repair_step, Some(5));
        assert!(matches!(
            &p.compiled.plan.steps[5],
            Operation::RepairEnds {
                template: Source::Step(4),
                method: EndRepairMethod::FillFivePrimeOverhangs
            }
        ));
        assert_eq!(
            p.prepared[0]
                .as_ref()
                .expect("prepared row must reconstruct")
                .sequence
                .as_str(),
            FILLED
        );
        assert_eq!(
            p.verification
                .product(6)
                .expect("selected step must reconstruct")
                .sequence
                .as_str(),
            letters
        );
        assert_eq!(
            assembly_diagnostics::overlap(&p).unwrap().rows[0].total_steps,
            6
        );
    }
}
fn binding(w: &Workspace, index: usize) -> EventBinding {
    let d = &w.documents[index];
    EventBinding {
        document_id: d.id.clone(),
        event_id: d.head_event.clone().unwrap(),
        revision: d.head.clone(),
    }
}
#[test]
fn exact_saved_origin_precedes_all_seven_preparation_steps_without_retargeting() {
    let stock = original();
    let target =
        Molecule::dna("Independent original stock PCR", ORIGINAL, Topology::Linear).unwrap();
    let plan = DesignPlan {
        name: "Original PCR scope".into(),
        target_revision: target.revision(),
        inputs: vec![input(stock.clone())],
        steps: vec![Operation::Pcr {
            template: Source::Input("original".into()),
            forward: Primer::new("F", "AGTCCGAT", "").unwrap(),
            reverse: Primer::new("R", "TGACGATC", "").unwrap(),
            template_strand: collagene_core::chemistry::Strand::Top,
        }],
    };
    let w = Workspace::from_documents("Origin", vec![stock, target], 1).unwrap();
    let commit = cloning_outputs::finalize(
        &w,
        &plan,
        binding(&w, 1),
        vec![InputEventBinding {
            input_id: "original".into(),
            event: binding(&w, 0),
        }],
        &cloning_outputs::Options {
            name: "Original saved material".into(),
            show_prepared: false,
            save_primers: false,
            allow_incomplete: false,
        },
    )
    .unwrap();
    let event = commit.receipt.construct.clone();
    let capsule = collagene_core::material::capsule::capsule(
        &commit.workspace,
        &Selection {
            run_id: commit.receipt.run.run_id.clone(),
            step_index: 0,
            event: event.clone(),
        },
    )
    .unwrap();
    let pin = serde_json::to_value(&capsule).unwrap();
    let prep = Preparation::FromRun {
        origin: capsule,
        preparation: RunPreparation::OverlapExtensionPcrCatalogDigest {
            conditions: None,
            request: request(),
            selection: selection(),
            fragment_index: 1,
            repair: Some(EndRepairMethod::FillFivePrimeOverhangs),
        },
    };
    let target = Molecule::dna("Independent edited repaired48", FILLED, Topology::Linear).unwrap();
    let p = a::preview(
        &[input(commit.workspace.revisions[&event.revision].clone())],
        &target,
        &linear(prep, false),
    )
    .unwrap();
    assert_eq!(
        p.verification.status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        p.verification.findings
    );
    assert_eq!(p.compiled.assembly_step, 7);
    assert_eq!(p.compiled.rows[0].repair_step, Some(6));
    let Operation::ReconstructRunProduct { origin, .. } = &p.compiled.plan.steps[0] else {
        panic!("closed origin")
    };
    assert_eq!(serde_json::to_value(origin).unwrap(), pin);
    assert!(matches!(
        &p.compiled.plan.steps[1],
        Operation::Pcr {
            template: Source::Step(0),
            ..
        }
    ));
    assert!(matches!(
        &p.compiled.plan.steps[2],
        Operation::Pcr {
            template: Source::Step(0),
            ..
        }
    ));
    assert!(matches!(
        &p.compiled.plan.steps[5],
        Operation::CatalogDigest {
            template: Source::Step(4),
            ..
        }
    ));
    assert_eq!(
        assembly_diagnostics::overlap(&p).unwrap().rows[0].total_steps,
        7
    );
}
#[test]
fn stale_catalog_missing_fragments_and_broken_suffix_ownership_are_rejected() {
    let inputs = vec![input(original())];
    let target = Molecule::dna("Independent repaired", FILLED, Topology::Linear).unwrap();
    let mut prep = preparation(Some(EndRepairMethod::FillFivePrimeOverhangs));
    if let Preparation::OverlapExtensionPcrCatalogDigest { selection, .. } = &mut prep {
        selection.catalog_content_hash = "0".repeat(64);
    }
    assert!(a::preflight(&inputs, &target, &linear(prep, false)).is_err());
    let mut prep = preparation(None);
    if let Preparation::OverlapExtensionPcrCatalogDigest { fragment_index, .. } = &mut prep {
        *fragment_index = 999;
    }
    let failed = a::preview(&inputs, &target, &linear(prep, false)).unwrap();
    assert!(!failed.preparation_complete);
    assert_eq!(failed.verification.reconstructed_count(), 4);
    let p = a::preview(
        &inputs,
        &target,
        &linear(
            preparation(Some(EndRepairMethod::FillFivePrimeOverhangs)),
            false,
        ),
    )
    .unwrap();
    for kind in 0..4 {
        let mut bad = p.clone();
        match kind {
            0 => {
                if let Operation::CatalogDigest { template, .. } = &mut bad.compiled.plan.steps[4] {
                    *template = Source::Step(2);
                }
            }
            1 => {
                if let Operation::RepairEnds { template, .. } = &mut bad.compiled.plan.steps[5] {
                    *template = Source::Step(3);
                }
            }
            2 => bad.compiled.rows[0].repair_step = Some(4),
            _ => bad.compiled.rows[0].prepared_source = Source::Step(3),
        }
        bad.verification.plan_revision =
            collagene_core::planning::plan::DesignPlan::revision(&bad.compiled.plan);
        assert!(bad.validate().is_err(), "broken suffix{kind}");
    }
}
