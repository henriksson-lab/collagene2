//! Literal target-free PCR → archive restart → real EcoRI fragments → independent
//! assembly target. References alone never manufacture private material state.
use collagene_core::{
    Feature, Location, Molecule, Topology,
    analysis::Primer,
    catalog_digest,
    chemistry::{Duplex, FivePrime},
    cloning_outputs,
    history::DocumentHistory,
    material_derivation as d, material_pcr, modified_duplex, run_digest_history,
    run_lineage::{self, EventBinding, InputEventBinding},
    run_material::{self, MaterialOrigin},
    verification::{self, DesignPlan, Operation, OrientedSource, PlanInput, Source, Status},
    workspace::{self, Workspace},
};
const STOCK: &str = "AGTCCGATGAATTCCGTAGGACTTCAGACCTGTAGATCGTCA";
const RIGHT: &str = "AATTCCGTAGGACTTCAGACCTGTAGATCGTCA";
fn event(w: &Workspace, index: usize) -> EventBinding {
    let doc = &w.documents[index];
    EventBinding {
        document_id: doc.id.clone(),
        event_id: doc.head_event.clone().unwrap(),
        revision: doc.head.clone(),
    }
}
fn source(w: &Workspace, pin: &EventBinding) -> Molecule {
    w.revisions[&pin.revision].clone()
}
fn enzyme() -> catalog_digest::CatalogSelection {
    let c = catalog_digest::builtin_catalog().unwrap();
    let e = c
        .edition()
        .enzymes
        .iter()
        .find(|e| e.name == "EcoRI")
        .unwrap();
    catalog_digest::select(&c, &[e.id.clone()]).unwrap()
}
fn primers() -> (Primer, Primer) {
    (
        Primer::new("Literal F", "AGTCCGAT", "").unwrap(),
        Primer::new("Literal R", "TGACGATC", "").unwrap(),
    )
}
fn pcr(
    w: &Workspace,
    pin: &EventBinding,
    origin: d::SourceOrigin,
    forward: Primer,
    reverse: Primer,
    name: &str,
) -> d::Request {
    let source = source(w, pin);
    let previews = match &origin {
        d::SourceOrigin::OriginalReference => {
            verification::simulate_pcr_documents(&source, &forward, &reverse).unwrap()
        }
        d::SourceOrigin::Saved(origin) => {
            material_pcr::simulate(&source, origin, &forward, &reverse).unwrap()
        }
    };
    assert_eq!(previews.len(), 1);
    d::Request {
        source: pin.clone(),
        origin,
        derivation: d::Derivation::ExactPcr {
            forward,
            reverse,
            selection: d::PcrSelection {
                selected_index: 0,
                expected_preview_revision: previews[0].molecule.revision(),
                geometry: previews[0].geometry.clone(),
            },
        },
        name: name.into(),
        save_primers: true,
    }
}
fn original() -> Workspace {
    let mut m = Molecule::dna("Literal original stock", STOCK, Topology::Linear).unwrap();
    m.features.push(Feature {
        kind: "misc_feature".into(),
        location: Location::range(18, 26),
        qualifiers: vec![("label".into(), Some("Source feature".into()))],
        needs_review: false,
    });
    Workspace::from_documents("No target until assembly", vec![m], 0).unwrap()
}
fn standalone() -> (Workspace, EventBinding) {
    let w = original();
    let (f, r) = primers();
    let request = pcr(
        &w,
        &event(&w, 0),
        d::SourceOrigin::OriginalReference,
        f,
        r,
        "Actual standalone PCR",
    );
    let before = workspace::save_workspace(&w).unwrap();
    let installed = d::prepare_install(&w, &request).unwrap();
    assert_eq!(workspace::save_workspace(&w).unwrap(), before);
    assert!(installed.receipt.primer_attachment.is_some());
    let product = installed.receipt.product;
    let w = workspace::load_workspace(&workspace::save_workspace(&installed.workspace).unwrap())
        .unwrap();
    (w, product)
}
#[test]
fn ordinary_pcr_restart_redigest_and_assembly_replay_form_one_complete_native_workflow() {
    let (w, pcr_pin) = standalone();
    assert_eq!(w.documents.len(), 2, "No fabricated independent target");
    assert!(run_lineage::list_runs(&w).unwrap().is_empty());
    let fresh = d::reconstruct(&w, &pcr_pin).unwrap();
    assert_eq!(fresh.reference().sequence.as_str(), STOCK);
    assert!(!fresh.incomplete());
    assert!(Duplex::from_molecule(fresh.reference()).is_err());
    assert_eq!(fresh.tracked().duplex.left.five_prime, FivePrime::Hydroxyl);
    assert_eq!(fresh.tracked().duplex.right.five_prime, FivePrime::Hydroxyl);
    let mods = modified_duplex::summary(fresh.tracked()).unwrap();
    assert_eq!(mods.known_modified, 0);
    assert_eq!(mods.unknown_eligible, 0);
    let pcr_capsule = d::capsule(&w, &pcr_pin).unwrap();
    let origin = d::SourceOrigin::Saved(MaterialOrigin::Derived(pcr_capsule));
    let preview = d::inspect_digest(
        "literal-derived-source",
        &source(&w, &pcr_pin),
        &origin,
        &enzyme(),
        &Default::default(),
    )
    .unwrap();
    assert!(preview.inspection.products_determined);
    assert_eq!(preview.inspection.fragments.len(), 2);
    assert_eq!(preview.fragment_modifications.len(), 2);
    let make = |index, name: &str| d::Request {
        source: pcr_pin.clone(),
        origin: origin.clone(),
        derivation: d::Derivation::CatalogDigest {
            selection: enzyme(),
            options: Default::default(),
            fragment_index: index,
        },
        name: name.into(),
        save_primers: false,
    };
    let left = d::prepare_install(&w, &make(0, "Actual left cut")).unwrap();
    let left_pin = left.receipt.product;
    let right = d::prepare_install(&left.workspace, &make(1, "Actual right cut")).unwrap();
    let right_pin = right.receipt.product;
    let mut w =
        workspace::load_workspace(&workspace::save_workspace(&right.workspace).unwrap()).unwrap();
    let right_fresh = d::reconstruct(&w, &right_pin).unwrap();
    assert_eq!(right_fresh.reference().sequence.as_str(), RIGHT);
    assert_eq!(right_fresh.tracked().duplex.bottom_start, 4);
    assert_eq!(
        right_fresh.reference().features[0].location.spans(),
        vec![(9, 17, false)]
    );
    let target = Molecule::dna("Independent assembly goal", STOCK, Topology::Linear).unwrap();
    DocumentHistory::new("independent-assembly-target", target.clone())
        .unwrap()
        .contribute(&mut w)
        .unwrap();
    let target_pin = event(&w, w.documents.len() - 1);
    let left_cap = d::capsule(&w, &left_pin).unwrap();
    let right_cap = d::capsule(&w, &right_pin).unwrap();
    let plan = DesignPlan {
        name: "Assemble freshly reconstructed sticky fragments".into(),
        target_revision: target.revision(),
        inputs: vec![
            PlanInput {
                id: "left".into(),
                revision: left_pin.revision.clone(),
                molecule: source(&w, &left_pin),
            },
            PlanInput {
                id: "right".into(),
                revision: right_pin.revision.clone(),
                molecule: source(&w, &right_pin),
            },
        ],
        steps: vec![
            Operation::ReconstructMaterial {
                template: Source::Input("left".into()),
                origin: left_cap,
            },
            Operation::ReconstructMaterial {
                template: Source::Input("right".into()),
                origin: right_cap,
            },
            Operation::Ligate {
                fragments: vec![
                    OrientedSource {
                        source: Source::Step(0),
                        reverse: false,
                    },
                    OrientedSource {
                        source: Source::Step(1),
                        reverse: false,
                    },
                ],
                circular: false,
                allow_nicked: false,
                options: Default::default(),
            },
        ],
    };
    let report = verification::verify(&plan, &target);
    assert_eq!(
        report.status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        report.findings
    );
    assert_eq!(report.products.last().unwrap().sequence.as_str(), STOCK);
    let committed = cloning_outputs::finalize(
        &w,
        &plan,
        target_pin,
        vec![
            InputEventBinding {
                input_id: "left".into(),
                event: left_pin,
            },
            InputEventBinding {
                input_id: "right".into(),
                event: right_pin,
            },
        ],
        &cloning_outputs::Options {
            name: "Final assembled product".into(),
            show_prepared: false,
            save_primers: false,
            allow_incomplete: false,
        },
    )
    .unwrap();
    let w = workspace::load_workspace(&workspace::save_workspace(&committed.workspace).unwrap())
        .unwrap();
    let selection = run_material::Selection {
        run_id: committed.receipt.run.run_id,
        step_index: 2,
        event: committed.receipt.construct,
    };
    let replay = run_material::reconstruct(&w, &selection).unwrap();
    assert_eq!(replay.reference().sequence.as_str(), STOCK);
    assert!(replay.tracked().duplex.nicks.is_empty());
}
#[test]
fn full_ordered_tails_and_source_run_and_legacy_digest_origins_are_replayed_not_renamed() {
    let w = original();
    let source_pin = event(&w, 0);
    let forward = Primer::new("Tail F", "AGTCCGAT", "GGC").unwrap();
    let reverse = Primer::new("Tail R", "TGACGATC", "ATG").unwrap();
    let request = pcr(
        &w,
        &source_pin,
        d::SourceOrigin::OriginalReference,
        forward.clone(),
        reverse.clone(),
        "Tailed actual PCR",
    );
    let installed = d::prepare_install(&w, &request).unwrap();
    assert_eq!(
        installed.workspace.revisions[&installed.receipt.product.revision]
            .sequence
            .as_str(),
        "GGCAGTCCGATGAATTCCGTAGGACTTCAGACCTGTAGATCGTCACAT"
    );
    let record =
        d::load_record(&installed.workspace.attachments[&installed.receipt.record_hash].bytes)
            .unwrap();
    let table = d::primer_table(&record).unwrap().unwrap();
    assert_eq!(
        table.records[0].raw_ordered_sequence.as_deref(),
        Some("GGCAGTCCGAT")
    );
    assert_eq!(
        table.records[1].raw_ordered_sequence.as_deref(),
        Some("ATGTGACGATC")
    );
    // Establish a conventional saved run independently, then derive target-free
    // PCR from that exact event and from its old target-free digest capsule.
    let mut w = original();
    let target = Molecule::dna("Declared conventional PCR goal", STOCK, Topology::Linear).unwrap();
    DocumentHistory::new("conventional-target", target.clone())
        .unwrap()
        .contribute(&mut w)
        .unwrap();
    let target_pin = event(&w, 1);
    let (f, r) = primers();
    let plan = DesignPlan {
        name: "Conventional original PCR".into(),
        target_revision: target.revision(),
        inputs: vec![PlanInput {
            id: "stock".into(),
            revision: source_pin.revision.clone(),
            molecule: source(&w, &source_pin),
        }],
        steps: vec![Operation::Pcr {
            template: Source::Input("stock".into()),
            forward: f.clone(),
            reverse: r.clone(),
        }],
    };
    let saved = cloning_outputs::finalize(
        &w,
        &plan,
        target_pin,
        vec![InputEventBinding {
            input_id: "stock".into(),
            event: source_pin,
        }],
        &cloning_outputs::Options {
            name: "Saved conventional PCR".into(),
            show_prepared: false,
            save_primers: false,
            allow_incomplete: false,
        },
    )
    .unwrap();
    let selection = run_material::Selection {
        run_id: saved.receipt.run.run_id,
        step_index: 0,
        event: saved.receipt.construct,
    };
    let origin = run_material::capsule(&saved.workspace, &selection).unwrap();
    let request = pcr(
        &saved.workspace,
        &selection.event,
        d::SourceOrigin::Saved(MaterialOrigin::Run(origin.clone())),
        f,
        r,
        "From exact saved run",
    );
    let after = d::prepare_install(&saved.workspace, &request).unwrap();
    assert_eq!(
        after.workspace.revisions[&after.receipt.product.revision]
            .sequence
            .as_str(),
        STOCK
    );
    let old = run_digest_history::prepare_install(
        &saved.workspace,
        &run_digest_history::Request {
            source: selection.event,
            origin,
            selection: enzyme(),
            options: Default::default(),
            fragment_index: 1,
        },
    )
    .unwrap();
    let pin = old.receipt.product;
    let cap = run_digest_history::capsule(&old.workspace, &pin).unwrap();
    let request = pcr(
        &old.workspace,
        &pin,
        d::SourceOrigin::Saved(MaterialOrigin::Digest(cap)),
        Primer::new("Inner F", "CCGTAGGA", "").unwrap(),
        Primer::new("Inner R", "TGACGATC", "").unwrap(),
        "From exact old digest",
    );
    let after = d::prepare_install(&old.workspace, &request).unwrap();
    assert_eq!(
        after.workspace.revisions[&after.receipt.product.revision]
            .sequence
            .as_str(),
        "CCGTAGGACTTCAGACCTGTAGATCGTCA"
    );
}
#[test]
fn geometry_source_event_and_opaque_record_forgery_are_rejected_without_mutating_input() {
    let w = original();
    let pin = event(&w, 0);
    let (f, r) = primers();
    let request = pcr(
        &w,
        &pin,
        d::SourceOrigin::OriginalReference,
        f,
        r,
        "Bound source PCR",
    );
    let before = workspace::save_workspace(&w).unwrap();
    let mut changed = request.clone();
    let d::Derivation::ExactPcr { selection, .. } = &mut changed.derivation else {
        unreachable!()
    };
    selection.expected_preview_revision = "0".repeat(64);
    assert!(d::prepare_install(&w, &changed).is_err());
    let mut changed = request.clone();
    let d::Derivation::ExactPcr { selection, .. } = &mut changed.derivation else {
        unreachable!()
    };
    selection.geometry.forward_start += 1;
    assert!(d::prepare_install(&w, &changed).is_err());
    assert_eq!(workspace::save_workspace(&w).unwrap(), before);
    let installed = d::prepare_install(&w, &request).unwrap();
    let product = installed.receipt.product;
    let mut cap = d::capsule(&installed.workspace, &product).unwrap();
    let parsed = d::load_record(&cap.record_bytes).unwrap();
    assert_eq!(parsed.product, product);
    let original = cap.record_bytes.clone();
    *cap.record_bytes.last_mut().unwrap() ^= 1;
    assert!(d::load_record(&cap.record_bytes).is_err());
    cap.record_bytes = original;
    cap.source_snapshot.name = "Tampered source".into();
    assert!(cap.validate_documentary().is_err());
    let reference = source(&installed.workspace, &product);
    let mut w = installed.workspace;
    DocumentHistory::new("equal-content-alias", reference)
        .unwrap()
        .contribute(&mut w)
        .unwrap();
    let alias = event(&w, w.documents.len() - 1);
    assert_eq!(alias.revision, product.revision);
    let cap = d::capsule(&w, &product).unwrap();
    assert!(d::validate_capsule_for_workspace(&w, &cap, &alias).is_err());
}
#[test]
fn saved_primer_table_is_atomic_presentation_and_omission_does_not_break_material_replay() {
    let (mut w, pin) = standalone();
    let record_hash = w
        .attachments
        .iter()
        .find(|(_, a)| a.kind == d::WORKSPACE_ATTACHMENT_KIND)
        .unwrap()
        .0
        .clone();
    let record = d::load_record(&w.attachments[&record_hash].bytes).unwrap();
    assert!(record.request.save_primers);
    assert!(d::primer_table(&record).unwrap().is_some());
    let keys = w
        .attachments
        .iter()
        .filter(|(_, a)| a.kind == "oligo-table-v1")
        .map(|(h, _)| h.clone())
        .collect::<Vec<_>>();
    assert_eq!(keys.len(), 1);
    for key in keys {
        w.attachments.remove(&key);
    }
    w.validate().unwrap();
    assert_eq!(
        d::reconstruct(&w, &pin)
            .unwrap()
            .reference()
            .sequence
            .as_str(),
        STOCK
    );
    let cap = d::capsule(&w, &pin).unwrap();
    let request = d::Request {
        source: pin.clone(),
        origin: d::SourceOrigin::Saved(MaterialOrigin::Derived(cap)),
        derivation: d::Derivation::CatalogDigest {
            selection: enzyme(),
            options: Default::default(),
            fragment_index: 1,
        },
        name: "No invented upstream primer order".into(),
        save_primers: true,
    };
    assert!(d::prepare_install(&w, &request).is_err());
}

#[test]
fn an_actual_reusable_incomplete_origin_stays_incomplete_after_new_unmodified_pcr() {
    use collagene_core::{homology_assembly, overlap_profiles};
    let a = Molecule::dna("Literal left source", "AGTCCGATGCTA", Topology::Linear).unwrap();
    let b = Molecule::dna(
        "Literal right source",
        "GCTATTCGAATTCCAGTC",
        Topology::Linear,
    )
    .unwrap();
    let target = Molecule::dna(
        "Independent source circle goal",
        "AGTCCGATGCTATTCGAATTCC",
        Topology::Circular,
    )
    .unwrap();
    let policy = homology_assembly::Options {
        overlaps: vec![4, 4],
        topology: Topology::Circular,
        min_paired_bases: 4,
        max_comparisons: 3 * (12 + 18) + 2 * 22 + 4 + 4,
        max_candidates: 100,
        ..Default::default()
    };
    let mut options = overlap_profiles::ProfileOptions::default();
    options.scan_repeated_overlaps = false;
    let plan = DesignPlan {
        name: "Actual sealed source with bounded alternatives".into(),
        target_revision: target.revision(),
        inputs: vec![
            PlanInput {
                id: "a".into(),
                revision: a.revision(),
                molecule: a.clone(),
            },
            PlanInput {
                id: "b".into(),
                revision: b.revision(),
                molecule: b.clone(),
            },
        ],
        steps: vec![
            Operation::Pcr {
                template: Source::Input("a".into()),
                forward: Primer::new("aF", "AGTC", "").unwrap(),
                reverse: Primer::new("aR", "TAGC", "").unwrap(),
            },
            Operation::Pcr {
                template: Source::Input("b".into()),
                forward: Primer::new("bF", "GCTA", "").unwrap(),
                reverse: Primer::new("bR", "GACT", "").unwrap(),
            },
            Operation::HomologyAssembly {
                fragments: vec![
                    OrientedSource {
                        source: Source::Step(0),
                        reverse: false,
                    },
                    OrientedSource {
                        source: Source::Step(1),
                        reverse: false,
                    },
                ],
                options,
                policy,
            },
        ],
    };
    assert_eq!(
        verification::verify(&plan, &target).status,
        Status::Incomplete
    );
    let w = Workspace::from_documents("Actual incomplete source", vec![a, b, target], 2).unwrap();
    let saved = cloning_outputs::finalize(
        &w,
        &plan,
        event(&w, 2),
        vec![
            InputEventBinding {
                input_id: "a".into(),
                event: event(&w, 0),
            },
            InputEventBinding {
                input_id: "b".into(),
                event: event(&w, 1),
            },
        ],
        &cloning_outputs::Options {
            name: "Reusable conditional circle".into(),
            show_prepared: false,
            save_primers: false,
            allow_incomplete: true,
        },
    )
    .unwrap();
    let selected = run_material::Selection {
        run_id: saved.receipt.run.run_id,
        step_index: 2,
        event: saved.receipt.construct,
    };
    let origin = run_material::capsule(&saved.workspace, &selected).unwrap();
    let request = pcr(
        &saved.workspace,
        &selected.event,
        d::SourceOrigin::Saved(MaterialOrigin::Run(origin)),
        Primer::new("Conditional F", "AGTCCGAT", "").unwrap(),
        Primer::new("Conditional R", "GGAATTC", "").unwrap(),
        "Conditional source PCR",
    );
    let prepared = d::prepare_install(&saved.workspace, &request).unwrap();
    assert!(prepared.receipt.incomplete);
    let product = prepared.receipt.product;
    let w = workspace::load_workspace(&workspace::save_workspace(&prepared.workspace).unwrap())
        .unwrap();
    let fresh = d::reconstruct(&w, &product).unwrap();
    assert!(fresh.incomplete());
    assert_eq!(
        fresh.reference().sequence.as_str(),
        "AGTCCGATGCTATTCGAATTCC"
    );
    assert!(fresh.findings_total() >= fresh.findings().len());
    let mods = modified_duplex::summary(fresh.tracked()).unwrap();
    assert_eq!(mods.known_modified, 0);
    assert_eq!(mods.unknown_eligible, 0);
}

#[test]
fn multiple_pcr_geometries_and_detached_chemical_summaries_never_become_implicit_material() {
    let repeated = Molecule::dna(
        "Repeated exact template",
        &format!("{STOCK}{STOCK}"),
        Topology::Linear,
    )
    .unwrap();
    let w = Workspace::from_documents("Multiple PCR products", vec![repeated.clone()], 0).unwrap();
    let (f, r) = primers();
    let previews = verification::simulate_pcr_documents(&repeated, &f, &r).unwrap();
    assert!(previews.len() > 1);
    let request = d::Request {
        source: event(&w, 0),
        origin: d::SourceOrigin::OriginalReference,
        derivation: d::Derivation::ExactPcr {
            forward: f,
            reverse: r,
            selection: d::PcrSelection {
                selected_index: 0,
                expected_preview_revision: previews[0].molecule.revision(),
                geometry: previews[0].geometry.clone(),
            },
        },
        name: "Explicit selection is not isolation".into(),
        save_primers: false,
    };
    let before = workspace::save_workspace(&w).unwrap();
    assert!(d::prepare_install(&w, &request).is_err());
    assert_eq!(workspace::save_workspace(&w).unwrap(), before);
    let (w, pin) = standalone();
    let cap = d::capsule(&w, &pin).unwrap();
    let (f, r) = primers();
    let mut request = pcr(
        &w,
        &pin,
        d::SourceOrigin::Saved(MaterialOrigin::Derived(cap)),
        f,
        r,
        "Exact material replay required",
    );
    request.origin = d::SourceOrigin::OriginalReference;
    assert!(
        d::prepare_install(&w, &request).is_err(),
        "Dropping the origin must not clear a physical source restriction"
    );
}
