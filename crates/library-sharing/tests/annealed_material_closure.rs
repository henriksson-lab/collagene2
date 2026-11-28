//! Selected adapter publication retains both exact oligo events, never private heads.
use collagene_core::{
    chemistry::{AnnealOptions, EndTreatment, FivePrime, Terminus, ThreePrime},
    documents::workspace::{self, RevisionEvent, Workspace},
    material::{
        derivation::{self as native, BoundSource, Derivation, SourceOrigin, SourceRole},
        lineage::EventBinding,
        run::MaterialOrigin,
    },
};
use collagene_library_sharing::{self as sharing, Selection};
use collagene_model::{Feature, Location, Molecule, Topology};
fn event(w: &Workspace, index: usize) -> EventBinding {
    let d = &w.documents[index];
    EventBinding {
        document_id: d.id.clone(),
        event_id: d.head_event.clone().unwrap(),
        revision: d.head.clone(),
    }
}
fn select(event: &EventBinding) -> Selection {
    Selection {
        document_id: event.document_id.clone(),
        revision: event.revision.clone(),
        head_event: Some(event.event_id.clone()),
    }
}
fn options(min_paired_bases: usize) -> AnnealOptions {
    AnnealOptions {
        min_paired_bases,
        top_chemistry: Terminus::PCR,
        bottom_chemistry: Terminus::PCR,
        ..Default::default()
    }
}
fn anneal(
    w: &Workspace,
    top: &EventBinding,
    bottom: &EventBinding,
    options: AnnealOptions,
) -> native::PreparedInstallation {
    let preview = native::annealing::preview_annealed_oligos(w, top, bottom, &options).unwrap();
    assert!(preview.preview.annealing.unique);
    let candidate = &preview.preview.annealing.candidates[0].id;
    native::prepare_install(
        w,
        &native::Request {
            sources: vec![
                BoundSource {
                    role: SourceRole::TopOligo,
                    event: top.clone(),
                    origin: SourceOrigin::OriginalReference,
                },
                BoundSource {
                    role: SourceRole::BottomOligo,
                    event: bottom.clone(),
                    origin: SourceOrigin::OriginalReference,
                },
            ],
            derivation: Derivation::AnnealedOligos {
                options,
                selection: preview.selection(candidate).unwrap(),
            },
            name: "Reusable declared annealed adapter".into(),
            save_primers: false,
        },
    )
    .unwrap()
}
fn fixture() -> (
    Workspace,
    EventBinding,
    EventBinding,
    EventBinding,
    String,
    String,
    String,
) {
    let mut top = Molecule::dna("Declared top oligo", "AATTCGTACT", Topology::Linear).unwrap();
    top.features.push(Feature {
        kind: "misc_feature".into(),
        location: Location::range(0, 4),
        qualifiers: vec![("label".into(), Some("Top 5prime overhang".into()))],
        needs_review: false,
    });
    let mut bottom = Molecule::dna("Declared bottom oligo", "AGTACG", Topology::Linear).unwrap();
    bottom.features.push(Feature {
        kind: "misc_feature".into(),
        location: Location::range(0, 2),
        qualifiers: vec![("label".into(), Some("Bottom source mark".into()))],
        needs_review: false,
    });
    let private = Molecule::dna("PRIVATE unrelated DNA", "GGGAAACCC", Topology::Circular).unwrap();
    let private_revision = private.revision();
    let w = Workspace::from_documents(
        "Two source adapter plus private lab document",
        vec![top, bottom, private],
        0,
    )
    .unwrap();
    let top = event(&w, 0);
    let bottom = event(&w, 1);
    let adapter = anneal(&w, &top, &bottom, options(6));
    let first_record = adapter.receipt.record_hash.clone();
    let source = adapter.receipt.product.clone();
    let treated = native::prepare_install(
        &adapter.workspace,
        &native::Request {
            sources: vec![BoundSource {
                role: SourceRole::Template,
                event: source.clone(),
                origin: SourceOrigin::Saved(MaterialOrigin::Derived(
                    native::capsule(&adapter.workspace, &source).unwrap(),
                )),
            }],
            derivation: Derivation::TerminalProcessing {
                action: native::TerminalAction::TreatEnds {
                    treatment: EndTreatment::Phosphorylate,
                },
            },
            name: "Phosphorylated saved adapter".into(),
            save_primers: false,
        },
    )
    .unwrap();
    let last_record = treated.receipt.record_hash;
    (
        treated.workspace,
        top,
        bottom,
        treated.receipt.product,
        first_record,
        last_record,
        private_revision,
    )
}
fn advance_private_head(w: &mut Workspace, index: usize, letters: &str) -> String {
    let old = event(w, index);
    let next = Molecule::dna("PRIVATE later oligo redesign", letters, Topology::Linear).unwrap();
    let revision = next.revision();
    let next_event = RevisionEvent::new(
        &old.document_id,
        &revision,
        vec![old.event_id.clone()],
        "Private later redesign after adapter retention",
    )
    .unwrap();
    w.revisions.insert(revision.clone(), next);
    w.history.push(next_event.clone());
    let d = &mut w.documents[index];
    d.undo.push(old.revision);
    d.undo_events.push(old.event_id);
    d.head = revision.clone();
    d.head_event = Some(next_event.id);
    revision
}
#[test]
fn selected_adapter_descendant_closes_both_original_roles_and_omits_later_private_heads() {
    let (mut w, top, bottom, product, first_record, last_record, private) = fixture();
    let new_top = advance_private_head(&mut w, 0, "GGGCCCAA");
    let new_bottom = advance_private_head(&mut w, 1, "TTTGGGCC");
    w.validate().unwrap();
    let before = workspace::save_workspace(&w).unwrap();
    let fresh = native::reconstruct(&w, &product).unwrap();
    assert_eq!(fresh.reference().sequence.as_str(), "AATTCGTACT");
    assert_eq!(fresh.tracked().duplex.bottom.as_str(), "AGTACG");
    assert_eq!(fresh.tracked().duplex.bottom_start, 4);
    assert_eq!(fresh.tracked().duplex.left.five_prime, FivePrime::Phosphate);
    assert_eq!(
        fresh.tracked().duplex.right.five_prime,
        FivePrime::Phosphate
    );
    assert_eq!(
        fresh.tracked().duplex.left.three_prime,
        ThreePrime::Hydroxyl
    );
    assert!(
        fresh
            .reference()
            .features
            .iter()
            .any(|f| f.label() == "Top 5prime overhang" && f.location == Location::range(0, 4))
    );
    assert!(
        fresh
            .reference()
            .features
            .iter()
            .any(|f| f.label() == "Bottom source mark"
                && f.location == Location::Complement(Box::new(Location::range(8, 10))))
    );
    for projected in [
        sharing::project_document_closure(&w, &select(&product)).unwrap(),
        sharing::project_private_document_closure(&w, &select(&product)).unwrap(),
    ] {
        assert!(projected.inventory.events.contains(&top.event_id));
        assert!(projected.inventory.events.contains(&bottom.event_id));
        assert!(projected.inventory.attachments.contains(&first_record));
        assert!(projected.inventory.attachments.contains(&last_record));
        assert!(!projected.workspace.revisions.contains_key(&private));
        assert!(!projected.workspace.revisions.contains_key(&new_top));
        assert!(!projected.workspace.revisions.contains_key(&new_bottom));
        for source in [&top, &bottom] {
            let restored = projected
                .workspace
                .documents
                .iter()
                .find(|d| d.id == source.document_id)
                .unwrap();
            assert_eq!(&restored.head, &source.revision);
            assert_eq!(restored.head_event.as_ref(), Some(&source.event_id));
        }
        let restarted =
            workspace::load_workspace(&workspace::save_workspace(&projected.workspace).unwrap())
                .unwrap();
        let replay = native::reconstruct(&restarted, &product).unwrap();
        assert_eq!(replay.tracked(), fresh.tracked());
        assert_eq!(replay.reference(), fresh.reference());
        let record = native::load_record(&restarted.attachments[&first_record].bytes).unwrap();
        assert_eq!(record.request.sources.len(), 2);
        assert_eq!(record.request.sources[0].role, SourceRole::TopOligo);
        assert_eq!(record.request.sources[1].role, SourceRole::BottomOligo);
        let event = restarted
            .history
            .iter()
            .find(|e| e.id == record.product.event_id)
            .unwrap();
        assert_eq!(
            event.parents,
            vec![top.event_id.clone(), bottom.event_id.clone()]
        );
    }
    assert_eq!(workspace::save_workspace(&w).unwrap(), before);
}
#[test]
fn same_exact_self_complementary_oligo_preserves_two_roles_and_one_unique_parent() {
    let w = Workspace::from_documents(
        "Self-complementary species, two explicit role uses",
        vec![
            Molecule::dna(
                "Self-complementary declared oligo",
                "AGTCGACT",
                Topology::Linear,
            )
            .unwrap(),
        ],
        0,
    )
    .unwrap();
    let source = event(&w, 0);
    let installed = anneal(&w, &source, &source, options(8));
    assert_eq!(installed.receipt.sources.len(), 2);
    assert_eq!(installed.receipt.sources[0].event, source);
    assert_eq!(installed.receipt.sources[1].event, source);
    let product = installed.receipt.product;
    let event = installed
        .workspace
        .history
        .iter()
        .find(|e| e.id == product.event_id)
        .unwrap();
    assert_eq!(event.parents, vec![source.event_id.clone()]);
    let projected =
        sharing::project_document_closure(&installed.workspace, &select(&product)).unwrap();
    assert_eq!(projected.workspace.documents.len(), 2);
    let capsule = native::capsule(&projected.workspace, &product).unwrap();
    assert_eq!(capsule.source_snapshots.len(), 2);
    assert_eq!(capsule.source_snapshots[0].role, SourceRole::TopOligo);
    assert_eq!(capsule.source_snapshots[1].role, SourceRole::BottomOligo);
    assert_eq!(
        capsule.source_snapshots[0].molecule,
        capsule.source_snapshots[1].molecule
    );
    let reference = &projected.workspace.revisions[&product.revision];
    assert_eq!(
        native::replay(&capsule, reference)
            .unwrap()
            .tracked()
            .duplex
            .bottom_start,
        0
    );
    let mut swapped = capsule.clone();
    swapped.source_snapshots.swap(0, 1);
    assert!(
        swapped.validate_documentary().is_err(),
        "Equal DNA does not erase strand roles"
    );
}
#[test]
fn selecting_unrelated_document_does_not_publish_adapter_records_or_either_private_oligo() {
    let (w, top, bottom, _product, first, last, private) = fixture();
    let private_event = event(&w, 2);
    assert_eq!(private_event.revision, private);
    let projected = sharing::project_document_closure(&w, &select(&private_event)).unwrap();
    assert_eq!(projected.workspace.documents.len(), 1);
    assert!(!projected.inventory.events.contains(&top.event_id));
    assert!(!projected.inventory.events.contains(&bottom.event_id));
    assert!(!projected.inventory.attachments.contains(&first));
    assert!(!projected.inventory.attachments.contains(&last));
    assert!(
        !projected
            .workspace
            .attachments
            .values()
            .any(|a| a.kind == native::WORKSPACE_ATTACHMENT_KIND)
    );
}
