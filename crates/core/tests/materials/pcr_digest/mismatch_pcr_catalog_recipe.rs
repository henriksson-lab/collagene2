//! Intentional mutation and real pinned cleavage, independent literal targets.
use collagene_core::{
    chemistry::{EndRepairMethod, FivePrime, Polarity},
    cloning::{
        assembly::recipe::{FragmentSpec, Preparation, RunPreparation},
        ligation::recipe as lig,
        pcr::mutagenesis::PcrOptions,
        restriction::catalog as catalog_digest,
    },
    documents::workspace::Workspace,
    material::{
        lineage::{EventBinding, InputEventBinding},
        outputs::{self as cloning_outputs, Options},
        run::Selection,
    },
    planning::plan::{DesignPlan, Operation, PlanInput, Source},
    sequence::analysis::Primer,
    verification::evidence::Status,
};
use collagene_model::{Molecule, Topology};
fn pin() -> catalog_digest::CatalogSelection {
    let c = catalog_digest::builtin_catalog().unwrap();
    let ids = ["EcoRI", "XhoI"]
        .iter()
        .map(|name| {
            c.edition()
                .enzymes
                .iter()
                .find(|e| e.name == *name)
                .unwrap()
                .id
                .clone()
        })
        .collect::<Vec<_>>();
    catalog_digest::select(&c, &ids).unwrap()
}
fn input(molecule: Molecule) -> PlanInput {
    PlanInput {
        id: "source".into(),
        revision: molecule.revision(),
        molecule,
    }
}
fn options() -> PcrOptions {
    let mut o = PcrOptions::default();
    o.binding.max_mismatches = 1;
    o.binding.min_three_prime_exact = 1;
    o
}
fn primers() -> (Primer, Primer) {
    (
        Primer::new("MutantF", "AGA", "CCGAATTC").unwrap(),
        Primer::new("R", "GGG", "GGCTCGAG").unwrap(),
    )
}
fn prep(repair: Option<EndRepairMethod>) -> Preparation {
    let (forward, reverse) = primers();
    Preparation::MismatchPcrCatalogDigest {
        conditions: None,
        forward,
        reverse,
        options: options(),
        selection: pin(),
        fragment_index: 1,
        repair,
    }
}
fn recipe(preparation: Preparation) -> lig::Recipe {
    lig::Recipe {
        name: "Mutant PCR and actual pinned digest".into(),
        fragments: vec![lig::Fragment {
            source: FragmentSpec {
                generation: None,
                role: "insert".into(),
                input_id: "source".into(),
                preparation,
                reverse: false,
            },
            treatment: None,
        }],
        circular: false,
        allow_nicked: false,
    }
}
fn target(dna: &str) -> Molecule {
    Molecule::dna(
        "Independent intended mutant fragment",
        dna,
        Topology::Linear,
    )
    .unwrap()
}
#[test]
fn ordinary_mismatch_amplicon_and_two_cuts_keep_real_outer_oh_and_internal_phosphate() {
    let inputs = vec![input(
        Molecule::dna("Original6bp", "AAACCC", Topology::Linear).unwrap(),
    )];
    let r = recipe(prep(None));
    let t = target("AATTCAGACCCC");
    let p = lig::preview(
        &inputs,
        &t,
        &r,
        collagene_core::cloning::assembly::recipe::CleavagePolicy::AssumedComplete,
    )
    .unwrap();
    assert_eq!(
        p.verification.status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        p.verification.findings
    );
    // F tail CCGAATTC + mutant domain AGA + original CCC + RC(GGCTCGAG).
    assert_eq!(
        p.verification
            .product(0)
            .expect("selected step must reconstruct")
            .sequence
            .as_str(),
        "CCGAATTCAGACCCCTCGAGCC"
    );
    assert_eq!(
        p.verification
            .product(1)
            .expect("selected step must reconstruct")
            .sequence
            .as_str(),
        "AATTCAGACCCC"
    );
    assert!(matches!(
        p.compiled.plan.steps[0],
        Operation::PcrWithMismatches { .. }
    ));
    assert!(matches!(
        p.compiled.plan.steps[1],
        Operation::CatalogDigest {
            template: Source::Step(0),
            fragment_index: 1,
            ..
        }
    ));
    let product = p.physical[0].as_ref().unwrap();
    let (left, right) = product.ends.as_ref().unwrap();
    assert_eq!((&*left.overhang, &*right.overhang), ("AATT", "TCGA"));
    assert_eq!(left.polarity, Polarity::FivePrime);
    assert_eq!(left.chemistry.five_prime, FivePrime::Phosphate);
    assert_eq!(right.chemistry.five_prime, FivePrime::Phosphate);
    let amplicon = serde_json::from_str::<collagene_core::chemistry::DuplexSummary>(
        &p.verification
            .product(0)
            .expect("selected step must reconstruct")
            .metadata["duplex_model"],
    )
    .unwrap();
    let (left, right) = amplicon.ends.unwrap();
    assert_eq!(left.chemistry.five_prime, FivePrime::Hydroxyl);
    assert_eq!(right.chemistry.five_prime, FivePrime::Hydroxyl);
    collagene_core::cloning::assembly::diagnostics::ligation(&p).unwrap();
    p.validate_binding(
        &lig::preview_revision(
            &inputs,
            &t,
            &r,
            collagene_core::cloning::assembly::recipe::CleavagePolicy::AssumedComplete,
        )
        .unwrap(),
    )
    .unwrap();
    let wrong = lig::preview(
        &inputs,
        &target("AATTCAAACCCC"),
        &r,
        collagene_core::cloning::assembly::recipe::CleavagePolicy::AssumedComplete,
    )
    .unwrap();
    assert_eq!(wrong.verification.status, Status::Invalid);
    assert!(wrong.product_physical.is_some());
}
fn binding(w: &Workspace, i: usize) -> EventBinding {
    let d = &w.documents[i];
    EventBinding {
        document_id: d.id.clone(),
        event_id: d.head_event.clone().unwrap(),
        revision: d.head.clone(),
    }
}
#[test]
fn saved_origin_chain_replays_before_mismatch_digest_and_explicit_repair() {
    let original = Molecule::dna("Original6bp", "AAACCC", Topology::Linear).unwrap();
    let original_target = target("AAACCC");
    let plan = DesignPlan {
        name: "Origin exactPCR".into(),
        inputs: vec![input(original.clone())],
        target_revision: original_target.revision(),
        steps: vec![Operation::Pcr {
            template: Source::Input("source".into()),
            forward: Primer::new("OriginF", "AAA", "").unwrap(),
            reverse: Primer::new("OriginR", "GGG", "").unwrap(),
            template_strand: collagene_core::chemistry::Strand::Top,
        }],
    };
    let w = Workspace::from_documents("Origin", vec![original, original_target], 1).unwrap();
    let commit = cloning_outputs::finalize(
        &w,
        &plan,
        binding(&w, 1),
        vec![InputEventBinding {
            input_id: "source".into(),
            event: binding(&w, 0),
        }],
        &Options {
            name: "Saved original PCR".into(),
            show_prepared: false,
            save_primers: false,
            allow_incomplete: false,
        },
    )
    .unwrap();
    let event = commit.receipt.construct;
    let origin = collagene_core::material::capsule::capsule(
        &commit.workspace,
        &Selection {
            run_id: commit.receipt.run.run_id,
            step_index: 0,
            event: event.clone(),
        },
    )
    .unwrap();
    let saved = commit.workspace.revisions[&event.revision].clone();
    let inputs = vec![input(saved)];
    let (forward, reverse) = primers();
    let r = recipe(Preparation::FromRun {
        origin,
        preparation: RunPreparation::MismatchPcrCatalogDigest {
            conditions: None,
            forward,
            reverse,
            options: options(),
            selection: pin(),
            fragment_index: 1,
            repair: Some(EndRepairMethod::Blunt),
        },
    });
    assert_eq!(
        r.fragments[0].source.preparation.catalog_selection(),
        Some(&pin())
    );
    // Filling the actual two 5′ extensions adds right TCGA to the top reference.
    let t = target("AATTCAGACCCCTCGA");
    let p = lig::preview(
        &inputs,
        &t,
        &r,
        collagene_core::cloning::assembly::recipe::CleavagePolicy::AssumedComplete,
    )
    .unwrap();
    assert_eq!(
        p.verification.status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        p.verification.findings
    );
    assert!(matches!(
        p.compiled.plan.steps[0],
        Operation::ReconstructRunProduct { .. }
    ));
    assert!(matches!(
        p.compiled.plan.steps[1],
        Operation::PcrWithMismatches {
            template: Source::Step(0),
            ..
        }
    ));
    assert!(matches!(
        p.compiled.plan.steps[2],
        Operation::CatalogDigest {
            template: Source::Step(1),
            ..
        }
    ));
    assert!(matches!(
        p.compiled.plan.steps[3],
        Operation::RepairEnds {
            template: Source::Step(2),
            ..
        }
    ));
    assert_eq!(p.compiled.rows[0].preparation.preparation_step, Some(0));
    assert_eq!(p.compiled.rows[0].preparation.repair_step, Some(3));
    let (left, right) = p.physical[0].as_ref().unwrap().ends.as_ref().unwrap();
    assert_eq!(left.polarity, Polarity::Blunt);
    assert_eq!(right.polarity, Polarity::Blunt);
    assert_eq!(left.chemistry.five_prime, FivePrime::Phosphate);
    assert_eq!(right.chemistry.five_prime, FivePrime::Phosphate);
    assert_eq!(
        collagene_core::cloning::assembly::diagnostics::ligation(&p)
            .unwrap()
            .rows[0]
            .total_steps,
        4
    );
    let mut bad = p.clone();
    if let Operation::CatalogDigest { template, .. } = &mut bad.compiled.plan.steps[2] {
        *template = Source::Step(0);
    }
    bad.verification.plan_revision =
        collagene_core::planning::plan::DesignPlan::revision(&bad.compiled.plan);
    assert!(bad.validate().is_err());
}
#[test]
fn incomplete_singleton_stays_conditionally_usable_and_never_becomes_verified_unique() {
    let inputs = vec![input(
        Molecule::dna("Repeated6bp domains", "AAACCCAAACCC", Topology::Linear).unwrap(),
    )];
    let mut preparation = prep(None);
    if let Preparation::MismatchPcrCatalogDigest { options, .. } = &mut preparation {
        options.max_products = 1;
    }
    let r = recipe(preparation);
    let p = lig::preview(
        &inputs,
        &target("AATTCAGACCCC"),
        &r,
        collagene_core::cloning::assembly::recipe::CleavagePolicy::AssumedComplete,
    )
    .unwrap();
    assert!(p.preparation_complete);
    assert!(p.product_physical.is_some());
    assert_eq!(
        p.verification.status,
        Status::Incomplete,
        "{:?}",
        p.verification.findings
    );
    assert!(
        p.verification
            .findings
            .iter()
            .any(|f| f.message.contains("incomplete"))
    );
    collagene_core::cloning::assembly::diagnostics::ligation(&p).unwrap();
    let wrong = lig::preview(
        &inputs,
        &target("AATTCAAACCCC"),
        &r,
        collagene_core::cloning::assembly::recipe::CleavagePolicy::AssumedComplete,
    )
    .unwrap();
    assert_eq!(wrong.verification.status, Status::Invalid);
}
#[test]
fn policy_pins_fragment_indices_and_secondary_operation_tampering_stay_bound() {
    let inputs = vec![input(
        Molecule::dna("Original6bp", "AAACCC", Topology::Linear).unwrap(),
    )];
    let r = recipe(prep(None));
    let t = target("AATTCAGACCCC");
    let p = lig::preview(
        &inputs,
        &t,
        &r,
        collagene_core::cloning::assembly::recipe::CleavagePolicy::AssumedComplete,
    )
    .unwrap();
    let mut changed = r.clone();
    if let Preparation::MismatchPcrCatalogDigest { options, .. } =
        &mut changed.fragments[0].source.preparation
    {
        options.binding.max_mismatches = 0;
    }
    assert_ne!(
        lig::preview_revision(
            &inputs,
            &t,
            &changed,
            collagene_core::cloning::assembly::recipe::CleavagePolicy::AssumedComplete
        )
        .unwrap(),
        p.context_revision
    );
    assert_ne!(
        lig::preview(
            &inputs,
            &t,
            &changed,
            collagene_core::cloning::assembly::recipe::CleavagePolicy::AssumedComplete
        )
        .unwrap()
        .verification
        .status,
        Status::ConsistentUnderAssumptions
    );
    let mut changed = r.clone();
    if let Preparation::MismatchPcrCatalogDigest { selection, .. } =
        &mut changed.fragments[0].source.preparation
    {
        selection.catalog_content_hash = "0".repeat(64);
    }
    assert!(
        lig::compile(
            &inputs,
            &t,
            &changed,
            collagene_core::cloning::assembly::recipe::CleavagePolicy::AssumedComplete
        )
        .is_err()
    );
    let mut changed = r.clone();
    if let Preparation::MismatchPcrCatalogDigest { fragment_index, .. } =
        &mut changed.fragments[0].source.preparation
    {
        *fragment_index = usize::MAX;
    }
    assert!(lig::source_preflight(&[("source", &inputs[0].molecule)], &t, &changed).is_err());
    let mut forged = p.clone();
    if let Operation::CatalogDigest { template, .. } = &mut forged.compiled.plan.steps[1] {
        *template = Source::Input("source".into());
    }
    forged.verification.plan_revision =
        collagene_core::planning::plan::DesignPlan::revision(&forged.compiled.plan);
    assert!(forged.validate().is_err());
    let decoded: lig::Preview = serde_json::from_slice(&serde_json::to_vec(&p).unwrap()).unwrap();
    decoded.validate_binding(&p.context_revision).unwrap();
}

#[test]
fn actual_bsa_i_mutant_fragment_passes_golden_gate_filter_and_saved_catalog_replay() {
    use collagene_core::cloning::golden_gate::GoldenGateEnzyme;
    use collagene_core::planning::combinatorial::{
        self as design, Bin, Choice, ChoiceSelection, Design, Method, SelectedConstruct,
        TargetDeclaration,
    };
    const ORIGINAL: &str = "GCTTACCGATGACCTAGTACGGCATTCGAAGC";
    const MUTANT: &str = "ACTTACCGATGACCTAGTACGGCATTCGAAGC";
    // Literal F11 + mutant32 + RC(R15). BsaI cuts off both recognition
    // sequences; ACTT is the top left5′ extension and TCAG the bottom right5′.
    const AMPLICON: &str = "ACGTGGTCTCAACTTACCGATGACCTAGTACGGCATTCGAAGCCTGATGAGACCACGT";
    let source = Molecule::dna("Original32bp BsaI source", ORIGINAL, Topology::Linear).unwrap();
    let intended = Molecule::dna("Independent mutant32bp", MUTANT, Topology::Linear).unwrap();
    let catalog = catalog_digest::builtin_catalog().unwrap();
    let select_one = |name: &str| {
        let id = catalog
            .edition()
            .enzymes
            .iter()
            .find(|e| e.name == name)
            .unwrap()
            .id
            .clone();
        catalog_digest::select(&catalog, &[id]).unwrap()
    };
    let selection = select_one("BsaI");
    let forward = Primer::new("BsaI mutantF", "ACTTACCG", "ACGTGGTCTCA").unwrap();
    let reverse = Primer::new("BsaI R", "GCTTCGAA", "ACGTGGTCTCATCAG").unwrap();
    let mut policy = options();
    policy.binding.min_three_prime_exact = 4;
    let preparation = Preparation::MismatchPcrCatalogDigest {
        conditions: None,
        forward: forward.clone(),
        reverse: reverse.clone(),
        options: policy.clone(),
        selection: selection.clone(),
        fragment_index: 1,
        repair: None,
    };
    assert_eq!(preparation.catalog_selection(), Some(&selection));
    let w = Workspace::from_documents(
        "Explicit single-fragment Golden Gate",
        vec![source.clone(), intended.clone()],
        1,
    )
    .unwrap();
    let mut declared = Design {
        material_bindings: vec![],
        id: "typeiis-library".into(),
        name: "Explicit mutant source preparation".into(),
        circular: false,
        bins: vec![Bin {
            id: "insert".into(),
            label: "Insert".into(),
            choices: vec![Choice {
                anneal_bottom: None,
                generation: None,
                id: "mutant".into(),
                label: "Explicit mismatch and BsaI".into(),
                source: binding(&w, 0),
                selector: collagene_core::cloning::assembly::intent::Selector::ReferenceRange {
                    start: 0,
                    end: 32,
                    reverse: false,
                },
                preparation_reverse: false,
                preparation: Some(preparation),
                treatment: None,
            }],
        }],
        method: Method::OrderedTypeIis {
            enzyme: GoldenGateEnzyme::BsaI,
            allow_nicked: false,
        },
    };
    let construct = SelectedConstruct {
        id: "checked-mutant".into(),
        name: "Checked mutant BsaI fragment".into(),
        choices: vec![ChoiceSelection {
            bin_id: "insert".into(),
            choice_id: "mutant".into(),
        }],
        target: TargetDeclaration::Existing(binding(&w, 1)),
        overrides: None,
    };
    let evaluated = design::evaluate(&declared, &construct, &w).unwrap();
    assert_eq!(
        evaluated.verification().status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        evaluated.verification().findings
    );
    assert_eq!(
        evaluated
            .verification()
            .product(0)
            .expect("selected step must reconstruct")
            .sequence
            .as_str(),
        AMPLICON
    );
    assert_eq!(
        evaluated
            .verification()
            .product(1)
            .expect("selected step must reconstruct")
            .sequence
            .as_str(),
        MUTANT
    );
    let design::Preview::Ligation(preview) = &evaluated.preview else {
        panic!("Golden Gate prepared fragments must retain ordinary actual ligation evidence")
    };
    let (left, right) = preview.physical[0].as_ref().unwrap().ends.as_ref().unwrap();
    assert_eq!((&*left.overhang, &*right.overhang), ("ACTT", "TCAG"));
    assert_eq!(left.chemistry.five_prime, FivePrime::Phosphate);
    assert_eq!(right.chemistry.five_prime, FivePrime::Phosphate);
    // Change ONLY the valid catalog pin: BbsI geometry cannot pass BsaI policy.
    if let Some(Preparation::MismatchPcrCatalogDigest { selection, .. }) =
        &mut declared.bins[0].choices[0].preparation
    {
        *selection = select_one("BbsI");
    }
    assert!(design::resolve(&declared, &construct, &w).is_err());
    // The existing saved CatalogDigest route remains eligible after the same
    // actual mismatch amplicon is finalized against its independent58bp target.
    let original_target =
        Molecule::dna("Independent full58bp amplicon", AMPLICON, Topology::Linear).unwrap();
    let origin_plan = DesignPlan {
        name: "Save actual BsaI-flanked mutant PCR".into(),
        target_revision: original_target.revision(),
        inputs: vec![input(source.clone())],
        steps: vec![Operation::PcrWithMismatches {
            template: Source::Input("source".into()),
            forward,
            reverse,
            options: policy,
        }],
    };
    let origin_w = Workspace::from_documents(
        "Original and independently specified amplicon",
        vec![source, original_target],
        1,
    )
    .unwrap();
    let committed = cloning_outputs::finalize(
        &origin_w,
        &origin_plan,
        binding(&origin_w, 1),
        vec![InputEventBinding {
            input_id: "source".into(),
            event: binding(&origin_w, 0),
        }],
        &Options {
            name: "Saved BsaI mutant amplicon".into(),
            show_prepared: false,
            save_primers: false,
            allow_incomplete: false,
        },
    )
    .unwrap();
    let event = committed.receipt.construct;
    let origin = collagene_core::material::capsule::capsule(
        &committed.workspace,
        &Selection {
            run_id: committed.receipt.run.run_id,
            step_index: 0,
            event: event.clone(),
        },
    )
    .unwrap();
    let saved = committed.workspace.revisions[&event.revision].clone();
    let existing = Preparation::FromRun {
        origin,
        preparation: RunPreparation::CatalogDigest {
            conditions: None,
            selection: selection.clone(),
            fragment_index: 1,
            repair: None,
        },
    };
    assert_eq!(existing.catalog_selection(), Some(&selection));
    let replayed = lig::preview(
        &[input(saved)],
        &intended,
        &recipe(existing),
        collagene_core::cloning::assembly::recipe::CleavagePolicy::AssumedComplete,
    )
    .unwrap();
    assert_eq!(
        replayed.verification.status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        replayed.verification.findings
    );
    assert_eq!(
        replayed
            .verification
            .final_product()
            .expect("final declared step must reconstruct")
            .sequence
            .as_str(),
        MUTANT
    );
    let (left, right) = replayed.physical[0]
        .as_ref()
        .unwrap()
        .ends
        .as_ref()
        .unwrap();
    assert_eq!((&*left.overhang, &*right.overhang), ("ACTT", "TCAG"));
}
