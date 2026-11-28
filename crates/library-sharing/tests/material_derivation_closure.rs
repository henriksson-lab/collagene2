//! Actual target-free PCR→digest followed by a separately verified physical run.
use collagene_core::{
    cloning::{overlap::assembly as homology_assembly, restriction::catalog as catalog_digest},
    documents::{
        transfer::Candidate,
        workspace::{self, RevisionEvent, Workspace, WorkspaceAttachment, WorkspaceDocument},
    },
    material::{
        derivation::{self as native, Derivation, PcrSelection, SourceOrigin},
        lineage::{EventBinding, InputEventBinding},
        outputs as cloning_outputs,
        run::{self as run_material, MaterialOrigin},
    },
    planning::plan::{DesignPlan, Operation, OrientedSource, PlanInput, Source},
    sequence::analysis::Primer,
    verification::{self, evidence::Status},
};
use collagene_library_sharing::{self as sharing, Selection};
use collagene_model::{Molecule, Topology};
const STOCK: &str = "AGTCCGATGAATTCCGTAGGACTTCAGACCTGTAGATCGTCA";
const FINAL: &str = "CACGTTCGAATTCCGTAGGACTTCAGACCTGTAGATCGTCA";
fn input(id: &str, molecule: Molecule) -> PlanInput {
    PlanInput {
        id: id.into(),
        revision: molecule.revision(),
        molecule,
    }
}
fn event(w: &Workspace, index: usize) -> EventBinding {
    let d = &w.documents[index];
    EventBinding {
        document_id: d.id.clone(),
        event_id: d.head_event.clone().unwrap(),
        revision: d.head.clone(),
    }
}
fn selection(e: &EventBinding) -> Selection {
    Selection {
        document_id: e.document_id.clone(),
        revision: e.revision.clone(),
        head_event: Some(e.event_id.clone()),
    }
}
struct Fixture {
    w: Workspace,
    final_event: EventBinding,
    final_run: run_material::Selection,
    digest: EventBinding,
    first_record: String,
    second_record: String,
    primer_table: String,
    private_head: String,
    source_id: String,
    original_source: String,
}
fn fixture() -> Fixture {
    let source = Molecule::dna("Original stock", STOCK, Topology::Linear).unwrap();
    let vector = Molecule::dna(
        "Independent vector stock",
        "CACGTTCGAATTCCGT",
        Topology::Linear,
    )
    .unwrap();
    let target = Molecule::dna(
        "Independent final assembly literal",
        FINAL,
        Topology::Linear,
    )
    .unwrap();
    let private = Molecule::dna(
        "UNRELATED PRIVATE DOCUMENT",
        "CCCCAAAAGGGGTTTT",
        Topology::Circular,
    )
    .unwrap();
    let w = Workspace::from_documents(
        "Private full local workspace",
        vec![source.clone(), vector.clone(), target.clone(), private],
        0,
    )
    .unwrap();
    let source_id = w.documents[0].id.clone();
    let original_source = source.revision();
    let forward = Primer::new("New original F", "AGTCCGAT", "").unwrap();
    let reverse = Primer::new("New original R", "TGACGATC", "").unwrap();
    let preview = verification::pcr::simulate_pcr_documents(
        &source,
        &forward,
        &reverse,
        collagene_core::chemistry::Strand::Top,
    )
    .unwrap()
    .remove(0);
    let first = native::prepare_install(
        &w,
        &native::Request {
            sources: vec![native::BoundSource {
                role: native::SourceRole::Template,
                event: event(&w, 0),
                origin: SourceOrigin::OriginalReference,
            }],
            derivation: Derivation::ExactPcr {
                forward,
                reverse,
                selection: PcrSelection {
                    selected_index: 0,
                    expected_preview_revision: preview.molecule.revision(),
                    geometry: preview.geometry,
                },
                template_strand: collagene_core::chemistry::Strand::Top,
            },
            name: "Actual target-free PCR".into(),
            save_primers: true,
        },
    )
    .unwrap();
    let first_record = first.receipt.record_hash.clone();
    let primer_table = first.receipt.primer_attachment.clone().unwrap();
    let c = catalog_digest::builtin_catalog().unwrap();
    let id = c
        .edition()
        .enzymes
        .iter()
        .find(|e| e.name == "EcoRI")
        .unwrap()
        .id
        .clone();
    let second = native::prepare_install(
        &first.workspace,
        &native::Request {
            sources: vec![native::BoundSource {
                role: native::SourceRole::Template,
                event: first.receipt.product.clone(),
                origin: SourceOrigin::Saved(MaterialOrigin::Derived(
                    native::capsule(&first.workspace, &first.receipt.product).unwrap(),
                )),
            }],
            derivation: Derivation::CatalogDigest {
                selection: catalog_digest::select(&c, &[id]).unwrap(),
                options: Default::default(),
                fragment_index: 1,
            },
            name: "Actual target-free EcoRI fragment".into(),
            save_primers: false,
        },
    )
    .unwrap();
    let second_record = second.receipt.record_hash.clone();
    let digest = second.receipt.product;
    let cut = second.workspace.revisions[&digest.revision].clone();
    assert_eq!(cut.sequence.as_str(), "AATTCCGTAGGACTTCAGACCTGTAGATCGTCA");
    let mut options = collagene_core::cloning::overlap::profiles::ProfileOptions::default();
    options.scan_repeated_overlaps = false;
    let plan = DesignPlan {
        name: "Saved derivation then physical assembly".into(),
        target_revision: target.revision(),
        inputs: vec![input("vector", vector), input("cut", cut)],
        steps: vec![
            Operation::Pcr {
                template: Source::Input("vector".into()),
                forward: Primer::new("Vector F", "CACGTTCG", "").unwrap(),
                reverse: Primer::new("Vector R", "ACGGAATT", "").unwrap(),
                template_strand: collagene_core::chemistry::Strand::Top,
            },
            Operation::ReconstructMaterial {
                template: Source::Input("cut".into()),
                origin: native::capsule(&second.workspace, &digest).unwrap(),
            },
            Operation::Pcr {
                template: Source::Step(1),
                forward: Primer::new("Insert F", "AATTCCGT", "").unwrap(),
                reverse: Primer::new("Insert R", "TGACGATC", "").unwrap(),
                template_strand: collagene_core::chemistry::Strand::Top,
            },
            Operation::HomologyAssembly {
                fragments: vec![
                    OrientedSource {
                        source: Source::Step(0),
                        reverse: false,
                    },
                    OrientedSource {
                        source: Source::Step(2),
                        reverse: false,
                    },
                ],
                options,
                policy: homology_assembly::Options {
                    overlaps: vec![8],
                    topology: Topology::Linear,
                    min_paired_bases: 8,
                    ..Default::default()
                },
            },
        ],
    };
    let report = verification::engine::verify(&plan, &target);
    assert_eq!(
        report.status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        report.findings
    );
    assert_eq!(report.products.len(), 4);
    assert!(verification::engine::plan_target_matches(
        &plan,
        report.products.last().unwrap(),
        &target
    ));
    let saved = cloning_outputs::finalize(
        &second.workspace,
        &plan,
        event(&second.workspace, 2),
        vec![
            InputEventBinding {
                input_id: "vector".into(),
                event: event(&second.workspace, 1),
            },
            InputEventBinding {
                input_id: "cut".into(),
                event: digest.clone(),
            },
        ],
        &Default::default(),
    )
    .unwrap();
    let final_event = saved.receipt.construct;
    let final_run = run_material::Selection {
        run_id: saved.receipt.run.run_id,
        step_index: 3,
        event: final_event.clone(),
    };
    let mut w = saved.workspace;
    // A later private current branch of the required source must stay local.
    let private = Molecule::dna(
        "Unpublished source branch",
        "GGGGCCCCAAAATTTT",
        Topology::Linear,
    )
    .unwrap();
    let private_head = private.revision();
    let previous = w.documents[0].head.clone();
    let previous_event = w.documents[0].head_event.clone().unwrap();
    let edited = RevisionEvent::new(
        &source_id,
        &private_head,
        vec![previous_event.clone()],
        "Private source edit after derivation",
    )
    .unwrap();
    w.revisions.insert(private_head.clone(), private);
    w.history.push(edited.clone());
    let doc = &mut w.documents[0];
    doc.undo.push(previous);
    doc.undo_events.push(previous_event);
    doc.head = private_head.clone();
    doc.head_event = Some(edited.id);
    let mut candidate = plan;
    candidate.name = "UNPUBLISHED CANDIDATE".into();
    w.candidates.push(Candidate::new(candidate));
    w.active_candidate = Some(w.candidates.last().unwrap().plan_revision.clone());
    let draft = WorkspaceAttachment {
        kind: collagene_core::cloning::assembly::library_history::ACTIVE_DRAFT_KIND.into(),
        name: "Private draft".into(),
        media_type: "application/json".into(),
        bytes: br#"{"private":"UNPUBLISHED DRAFT"}"#.to_vec(),
    };
    w.attachments
        .insert(blake3::hash(&draft.bytes).to_hex().to_string(), draft);
    w.validate().unwrap();
    Fixture {
        w,
        final_event,
        final_run,
        digest,
        first_record,
        second_record,
        primer_table,
        private_head,
        source_id,
        original_source,
    }
}
fn check_projected(f: &Fixture, projection: sharing::Projection) {
    let w = projection.workspace;
    assert_eq!(w.active_document, f.final_event.document_id);
    assert!(w.candidates.is_empty());
    assert!(!w.revisions.contains_key(&f.private_head));
    assert!(!w.documents.iter().any(|d| d.id == f.w.documents[3].id));
    assert_eq!(
        w.documents
            .iter()
            .find(|d| d.id == f.source_id)
            .unwrap()
            .head,
        f.original_source
    );
    assert!(w.attachments.contains_key(&f.first_record));
    assert!(w.attachments.contains_key(&f.second_record));
    assert!(!w.attachments.contains_key(&f.primer_table));
    assert!(
        !w.attachments.values().any(
            |a| a.kind == collagene_core::cloning::assembly::library_history::ACTIVE_DRAFT_KIND
        )
    );
    let restored = workspace::load_workspace(&workspace::save_workspace(&w).unwrap()).unwrap();
    let fresh = run_material::reconstruct(&restored, &f.final_run).unwrap();
    assert_eq!(fresh.reference().sequence.as_str(), FINAL);
    assert!(fresh.tracked().duplex.nicks.is_empty());
}
#[test]
fn public_and_private_selected_closures_keep_real_derivation_and_run_replay_without_private_heads()
{
    let f = fixture();
    let before = workspace::save_workspace(&f.w).unwrap();
    let selected = selection(&f.final_event);
    check_projected(
        &f,
        sharing::project_document_closure(&f.w, &selected).unwrap(),
    );
    check_projected(
        &f,
        sharing::project_private_document_closure(&f.w, &selected).unwrap(),
    );
    assert_eq!(workspace::save_workspace(&f.w).unwrap(), before);
}
#[test]
fn selecting_only_derived_fragment_keeps_typed_ancestors_without_unrelated_downstream_run() {
    let f = fixture();
    let projected = sharing::project_document_closure(&f.w, &selection(&f.digest)).unwrap();
    assert_eq!(projected.workspace.active_document, f.digest.document_id);
    assert!(
        projected
            .workspace
            .attachments
            .contains_key(&f.first_record)
    );
    assert!(
        projected
            .workspace
            .attachments
            .contains_key(&f.second_record)
    );
    assert!(
        collagene_core::material::lineage::list_runs(&projected.workspace)
            .unwrap()
            .is_empty()
    );
    let fresh = native::reconstruct(&projected.workspace, &f.digest).unwrap();
    assert_eq!(
        fresh.reference().sequence.as_str(),
        "AATTCCGTAGGACTTCAGACCTGTAGATCGTCA"
    );
    assert_eq!(fresh.tracked().duplex.bottom_start, 4);
}
#[test]
fn missing_owner_wrong_hash_and_exact_event_aliases_cannot_supply_derivation_authority() {
    let f = fixture();
    let selected = selection(&f.final_event);
    let mut missing = f.w.clone();
    missing.attachments.remove(&f.first_record);
    assert!(sharing::project_document_closure(&missing, &selected).is_err());
    let mut wrong = f.w.clone();
    let asset = wrong.attachments.remove(&f.second_record).unwrap();
    wrong.attachments.insert("00".repeat(32), asset);
    assert!(sharing::project_document_closure(&wrong, &selected).is_err());
    let mut w = f.w.clone();
    let alias_id = "equal-content-alias".to_string();
    let alias_event = RevisionEvent::new(
        &alias_id,
        &f.digest.revision,
        vec![],
        "Imported equal reference",
    )
    .unwrap();
    let alias = EventBinding {
        document_id: alias_id.clone(),
        event_id: alias_event.id.clone(),
        revision: f.digest.revision.clone(),
    };
    w.documents.push(WorkspaceDocument {
        id: alias_id,
        head: alias.revision.clone(),
        undo: vec![],
        redo: vec![],
        head_event: Some(alias.event_id.clone()),
        undo_events: vec![],
        redo_events: vec![],
    });
    w.history.push(alias_event);
    w.validate().unwrap();
    let capsule = native::capsule(&w, &f.digest).unwrap();
    assert!(native::validate_capsule_for_workspace(&w, &capsule, &alias).is_err());
    // Rebuild a valid event ID rather than relying on a generic hash mismatch.
    let orphan = RevisionEvent::new(
        &alias.document_id,
        &alias.revision,
        vec![],
        native::EVENT_OPERATION,
    )
    .unwrap();
    w.documents.last_mut().unwrap().head_event = Some(orphan.id.clone());
    let orphan_pin = EventBinding {
        document_id: alias.document_id,
        event_id: orphan.id.clone(),
        revision: alias.revision,
    };
    *w.history.last_mut().unwrap() = orphan;
    assert!(sharing::project_document_closure(&w, &selection(&orphan_pin)).is_err());
}
