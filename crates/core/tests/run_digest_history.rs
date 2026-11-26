//! Literal independent saved-PCR/circle cut products. A derived fragment is not
//! given a fabricated independent target merely to fit the cloning-run schema.
use collagene_core::{
    Feature, Location, Molecule, Topology,
    analysis::Primer,
    catalog_digest,
    chemistry::{Duplex, FivePrime},
    cloning_outputs, homology_assembly,
    methylation::MethylationState,
    modified_duplex::{self, GenerationAssumptions},
    run_digest_history as d,
    run_lineage::{EventBinding, InputEventBinding},
    run_material, transfer,
    verification::{self, DesignPlan, Operation, OrientedSource, PlanInput, Source, Status},
    workspace::{self, RevisionEvent, Workspace},
};
const STOCK: &str = "AGTCCGATGAATTCCGTAGGACTTCAGACCTGTAGATCGTCA";
const RIGHT_FRAGMENT: &str = "AATTCCGTAGGACTTCAGACCTGTAGATCGTCA";
fn event(w: &Workspace, index: usize) -> EventBinding {
    let doc = &w.documents[index];
    EventBinding {
        document_id: doc.id.clone(),
        event_id: doc.head_event.clone().unwrap(),
        revision: doc.head.clone(),
    }
}
fn enzyme() -> catalog_digest::CatalogSelection {
    let c = catalog_digest::builtin_catalog().unwrap();
    let id = c
        .edition()
        .enzymes
        .iter()
        .find(|e| e.name == "EcoRI")
        .unwrap()
        .id
        .clone();
    catalog_digest::select(&c, &[id]).unwrap()
}
fn input(id: &str, m: Molecule) -> PlanInput {
    PlanInput {
        id: id.into(),
        revision: m.revision(),
        molecule: m,
    }
}
fn finalize(
    w: Workspace,
    plan: DesignPlan,
    target: usize,
    sources: Vec<(&str, usize)>,
    incomplete: bool,
) -> (Workspace, d::Request) {
    let inputs = sources
        .into_iter()
        .map(|(id, index)| InputEventBinding {
            input_id: id.into(),
            event: event(&w, index),
        })
        .collect();
    let commit = cloning_outputs::finalize(
        &w,
        &plan,
        event(&w, target),
        inputs,
        &cloning_outputs::Options {
            name: "Exact saved native source".into(),
            show_prepared: false,
            save_primers: false,
            allow_incomplete: incomplete,
        },
    )
    .unwrap();
    let source = commit.receipt.construct.clone();
    let selection = run_material::Selection {
        run_id: commit.receipt.run.run_id.clone(),
        step_index: plan.steps.len() - 1,
        event: source.clone(),
    };
    let origin = run_material::capsule(&commit.workspace, &selection).unwrap();
    (
        commit.workspace,
        d::Request {
            source,
            origin,
            selection: enzyme(),
            options: Default::default(),
            fragment_index: 0,
        },
    )
}
fn pcr() -> (Workspace, d::Request) {
    let mut source = Molecule::dna("Original source stock", STOCK, Topology::Linear).unwrap();
    source.features.push(Feature {
        kind: "misc_feature".into(),
        location: Location::range(18, 26),
        qualifiers: vec![(
            "label".into(),
            Some("Retained actual source annotation".into()),
        )],
        needs_review: false,
    });
    let target = Molecule::dna("Independent full PCR target", STOCK, Topology::Linear).unwrap();
    let plan = DesignPlan {
        name: "Actual original PCR".into(),
        target_revision: target.revision(),
        inputs: vec![input("stock", source.clone())],
        steps: vec![Operation::Pcr {
            template: Source::Input("stock".into()),
            forward: Primer::new("F", "AGTCCGAT", "").unwrap(),
            reverse: Primer::new("R", "TGACGATC", "").unwrap(),
        }],
    };
    let w = Workspace::from_documents("Native original PCR", vec![source, target], 1).unwrap();
    let (w, mut request) = finalize(w, plan, 1, vec![("stock", 0)], false);
    request.fragment_index = 1;
    (w, request)
}
#[test]
fn real_saved_pcr_opens_an_independent_cut_fragment_without_target_bootstrap_and_replays_after_archive()
 {
    let (w, request) = pcr();
    let before = workspace::save_workspace(&w).unwrap();
    let context = d::installation_revision(&w, &request).unwrap();
    let prepared = d::prepare_install(&w, &request).unwrap();
    assert_eq!(prepared.receipt.context_revision, context);
    assert!(!prepared.receipt.origin_incomplete);
    prepared
        .receipt
        .validate_archive(&prepared.workspace)
        .unwrap();
    assert_eq!(workspace::save_workspace(&w).unwrap(), before);
    let product = &prepared.receipt.product;
    let reference = &prepared.workspace.revisions[&product.revision];
    assert_eq!(reference.sequence.as_str(), RIGHT_FRAGMENT);
    assert!(Duplex::from_molecule(reference).is_err());
    assert_eq!(reference.features[0].location.spans(), vec![(9, 17, false)]);
    assert_eq!(prepared.workspace.documents.len(), w.documents.len() + 1);
    assert_eq!(prepared.workspace.candidates.len(), w.candidates.len());
    assert_eq!(
        collagene_core::run_lineage::list_runs(&prepared.workspace)
            .unwrap()
            .len(),
        collagene_core::run_lineage::list_runs(&w).unwrap().len(),
        "No synthetic target or extra goal-checked cloning run was created"
    );
    let reopened =
        workspace::load_workspace(&workspace::save_workspace(&prepared.workspace).unwrap())
            .unwrap();
    let fresh = d::reconstruct(&reopened, product).unwrap();
    assert_eq!(fresh.reference().sequence.as_str(), RIGHT_FRAGMENT);
    assert_eq!(fresh.tracked().duplex.bottom_start, 4);
    assert_eq!(fresh.tracked().duplex.left.five_prime, FivePrime::Phosphate);
    assert_eq!(fresh.tracked().duplex.right.five_prime, FivePrime::Hydroxyl);
    assert!(fresh.tracked().duplex.nicks.is_empty());
    let summary = modified_duplex::summary(fresh.tracked()).unwrap();
    assert_eq!(summary.known_modified, 0);
    assert_eq!(summary.unknown_eligible, 0);
    let capsule = d::capsule(&reopened, product).unwrap();
    d::validate_capsule_for_workspace(&reopened, &capsule, product).unwrap();
    let received: d::Capsule =
        serde_json::from_slice(&serde_json::to_vec(&capsule).unwrap()).unwrap();
    assert_eq!(
        d::replay(&received, reference).unwrap().tracked(),
        fresh.tracked()
    );
    let record =
        d::load_record(&reopened.attachments[&prepared.receipt.record_hash].bytes).unwrap();
    assert!(
        d::source_bindings(&record)
            .unwrap()
            .contains(&request.source)
    );
    assert!(
        d::source_bindings(&record)
            .unwrap()
            .contains(&request.origin.target_event)
    );
    assert!(
        d::required_run_ids(&record)
            .unwrap()
            .contains(&request.origin.selection.run_id)
    );
    let installed = d::install(&w, &request).unwrap();
    installed.validate_binding(&context).unwrap();
}
#[test]
fn actual_physical_overlap_circle_digest_preserves_strands_and_inherited_search_limits() {
    let a = Molecule::dna("Actual left source", "AGTCCGATGCTA", Topology::Linear).unwrap();
    let b = Molecule::dna(
        "Actual right source",
        "GCTATTCGAATTCCAGTC",
        Topology::Linear,
    )
    .unwrap();
    let target = Molecule::dna(
        "Independent circle",
        "AGTCCGATGCTATTCGAATTCC",
        Topology::Circular,
    )
    .unwrap();
    let policy = homology_assembly::Options {
        overlaps: vec![4, 4],
        topology: Topology::Circular,
        min_paired_bases: 4,
        // Establish both original strand pairs, both declared homologies and
        // the actual product pair; cap only the subsequent alternative search.
        max_comparisons: 3 * (12 + 18) + 2 * 22 + 4 + 4,
        max_candidates: 100,
        ..Default::default()
    };
    let mut profile = collagene_core::overlap_profiles::ProfileOptions::default();
    profile.scan_repeated_overlaps = false;
    let plan = DesignPlan {
        name: "Two actual PCRs and physical closure".into(),
        target_revision: target.revision(),
        inputs: vec![input("left", a.clone()), input("right", b.clone())],
        steps: vec![
            Operation::Pcr {
                template: Source::Input("left".into()),
                forward: Primer::new("LF", "AGTC", "").unwrap(),
                reverse: Primer::new("LR", "TAGC", "").unwrap(),
            },
            Operation::Pcr {
                template: Source::Input("right".into()),
                forward: Primer::new("RF", "GCTA", "").unwrap(),
                reverse: Primer::new("RR", "GACT", "").unwrap(),
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
                options: profile,
                policy,
            },
        ],
    };
    let report = verification::verify(&plan, &target);
    assert_eq!(report.status, Status::Incomplete, "{:?}", report.findings);
    assert_eq!(
        report.products.len(),
        plan.steps.len(),
        "{:?}",
        report.findings
    );
    assert!(verification::plan_target_matches(
        &plan,
        report.products.last().unwrap(),
        &target
    ));
    let w =
        Workspace::from_documents("Actual saved physical circle", vec![a, b, target], 2).unwrap();
    let (w, request) = finalize(w, plan, 2, vec![("left", 0), ("right", 1)], true);
    let installed = d::prepare_install(&w, &request).unwrap();
    assert!(installed.receipt.origin_incomplete);
    assert!(installed.receipt.origin_findings_total >= installed.receipt.origin_findings.len());
    let fresh = d::reconstruct(&installed.workspace, &installed.receipt.product).unwrap();
    assert_eq!(
        fresh.reference().sequence.as_str(),
        "AATTCCAGTCCGATGCTATTCG"
    );
    assert_eq!(fresh.reference().topology, Topology::Linear);
    assert_eq!(fresh.tracked().duplex.top.len(), 22);
    assert_eq!(fresh.tracked().duplex.bottom.len(), 22);
    assert_eq!(fresh.tracked().duplex.bottom_start, 4);
    assert!(fresh.origin_incomplete());
    assert_eq!(fresh.tracked().duplex.left.five_prime, FivePrime::Phosphate);
    assert_eq!(
        fresh.tracked().duplex.right.five_prime,
        FivePrime::Phosphate
    );
}
#[test]
fn declared_original_modifications_survive_actual_full_circle_cleavage_and_capsule_replay() {
    let source = Molecule::dna(
        "Explicitly assumed original Dam stock",
        STOCK,
        Topology::Circular,
    )
    .unwrap();
    let target = Molecule::dna("Independent original circle", STOCK, Topology::Circular).unwrap();
    let plan = DesignPlan {
        name: "Closed source-generation assumption".into(),
        target_revision: target.revision(),
        inputs: vec![input("stock", source.clone())],
        steps: vec![Operation::DeclareSourceModifications {
            template: Source::Input("stock".into()),
            assumptions: GenerationAssumptions {
                dam: MethylationState::FullyModified,
                dcm: MethylationState::Unknown,
                cpg: MethylationState::Unknown,
            },
        }],
    };
    let w = Workspace::from_documents("Original declared source", vec![source, target], 1).unwrap();
    let (w, request) = finalize(w, plan, 1, vec![("stock", 0)], false);
    let original = run_material::reconstruct(&w, &request.origin.selection).unwrap();
    let expected = modified_duplex::summary(original.tracked()).unwrap();
    assert!(expected.known_modified > 0);
    let prepared = d::prepare_install(&w, &request).unwrap();
    let fresh = d::reconstruct(&prepared.workspace, &prepared.receipt.product).unwrap();
    let actual = modified_duplex::summary(fresh.tracked()).unwrap();
    assert_eq!(actual.known_modified, expected.known_modified);
    assert_eq!(actual.unknown_eligible, expected.unknown_eligible);
    assert_eq!(
        serde_json::to_value(&actual.generation).unwrap(),
        serde_json::to_value(&expected.generation).unwrap()
    );
    let capsule = d::capsule(&prepared.workspace, &prepared.receipt.product).unwrap();
    assert_eq!(
        d::replay(&capsule, fresh.reference()).unwrap().tracked(),
        fresh.tracked()
    );
}
#[test]
fn exact_source_origin_catalog_and_fragment_guards_refuse_changes_without_partial_install() {
    let (w, request) = pcr();
    let before = workspace::save_workspace(&w).unwrap();
    let mut bad = request.clone();
    bad.source = event(&w, 0);
    assert!(d::prepare_install(&w, &bad).is_err());
    let mut bad = request.clone();
    bad.fragment_index = 999;
    assert!(d::prepare_install(&w, &bad).is_err());
    let mut bad = request.clone();
    bad.options.max_comparisons = 1;
    assert!(d::prepare_install(&w, &bad).is_err());
    let mut bad = request.clone();
    bad.selection.catalog_content_hash = "00".repeat(32);
    assert!(d::prepare_install(&w, &bad).is_err());
    // A bounded, internally consistent recipe with the same DNA is still not
    // the exact retained originating run and cannot attest an installation.
    let mut bad = request.clone();
    let mut candidate = transfer::load_candidate(&bad.origin.candidate).unwrap();
    candidate.plan.name = "Different source recipe identity".into();
    bad.origin.candidate =
        transfer::save_candidate(&transfer::Candidate::new(candidate.plan)).unwrap();
    assert!(d::prepare_install(&w, &bad).is_err());
    let mut changed = w.clone();
    let mut histories = collagene_core::history::DocumentHistory::restore_all(&changed).unwrap();
    let source = histories
        .iter_mut()
        .find(|history| history.id == request.source.document_id)
        .unwrap();
    let mut renamed = source.molecule.clone();
    renamed.name = "A later saved source head".into();
    source
        .commit(renamed, "Name change after digest request")
        .unwrap();
    changed.documents.clear();
    changed.history.clear();
    for history in &histories {
        history.contribute(&mut changed).unwrap();
    }
    changed.validate().unwrap();
    assert!(d::prepare_install(&changed, &request).is_err());
    assert_eq!(workspace::save_workspace(&w).unwrap(), before);
    let installed = d::prepare_install(&w, &request).unwrap();
    let capsule = d::capsule(&installed.workspace, &installed.receipt.product).unwrap();
    let mut detached = capsule.source.clone();
    detached.metadata.remove("duplex_model");
    assert!(d::replay(&capsule, &detached).is_err());
}
#[test]
fn historical_science_record_remains_readable_but_cannot_restore_current_material_authority() {
    let (w, request) = pcr();
    let prepared = d::prepare_install(&w, &request).unwrap();
    let mut frame = prepared.workspace;
    let old = prepared.receipt.product;
    let mut record = d::load_record(
        &frame
            .attachments
            .remove(&prepared.receipt.record_hash)
            .unwrap()
            .bytes,
    )
    .unwrap();
    record.scientific_engine = "collagene-core/science-historical-fixture".into();
    record.context_revision = blake3::hash(
        &serde_json::to_vec(&(
            d::ENGINE_VERSION,
            &record.scientific_engine,
            &record.source_workspace_fingerprint,
            &record.request,
        ))
        .unwrap(),
    )
    .to_hex()
    .to_string();
    let mut product = frame.revisions.remove(&old.revision).unwrap();
    product.metadata.insert(
        "source:run_digest_history:context".into(),
        record.context_revision.clone(),
    );
    let revision = product.revision();
    let old_event = frame
        .history
        .iter()
        .find(|event| event.id == old.event_id)
        .unwrap()
        .clone();
    let event = RevisionEvent::new(
        &old.document_id,
        &revision,
        old_event.parents,
        old_event.operation,
    )
    .unwrap();
    frame.history.retain(|entry| entry.id != old.event_id);
    frame.history.push(event.clone());
    frame.revisions.insert(revision.clone(), product);
    let document = frame
        .documents
        .iter_mut()
        .find(|doc| doc.id == old.document_id)
        .unwrap();
    document.head = revision.clone();
    document.head_event = Some(event.id.clone());
    record.product = EventBinding {
        document_id: old.document_id,
        event_id: event.id,
        revision,
    };
    record.validate_documentary().unwrap();
    let bytes = serde_json::to_vec(&record).unwrap();
    let hash = blake3::hash(&bytes).to_hex().to_string();
    frame.attachments.insert(
        hash,
        workspace::WorkspaceAttachment {
            kind: d::WORKSPACE_ATTACHMENT_KIND.into(),
            name: "Historical target-free fragment".into(),
            media_type: "application/json".into(),
            bytes: bytes.clone(),
        },
    );
    d::load_record(&bytes).unwrap();
    d::validate_workspace_records(&frame).unwrap();
    let restored = workspace::load_workspace(&workspace::save_workspace(&frame).unwrap()).unwrap();
    assert!(d::reconstruct(&restored, &record.product).is_err());
    assert!(d::capsule(&restored, &record.product).is_err());
    let capsule = d::Capsule {
        schema_version: 1,
        source: restored.revisions[&record.request.source.revision].clone(),
        record,
    };
    capsule.validate_documentary().unwrap();
    assert!(
        d::replay(
            &capsule,
            &restored.revisions[&capsule.record.product.revision]
        )
        .is_err()
    );
}
