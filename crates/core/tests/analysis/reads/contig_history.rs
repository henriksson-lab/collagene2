use collagene_core::{
    documents::{
        history::DocumentHistory,
        workspace::{self, Workspace},
    },
    material::lineage::{EventBinding, InputEventBinding},
    review::contigs::{self as contig, ContigOptions, ReadInput, history::*},
};
use collagene_model::{Molecule, Topology};
fn fixture() -> (Workspace, SavedAssembly, Vec<ReadInput>) {
    let sequence = "AACCGTAGTCTGACGATGCCTAGGTCATCA";
    let molecules = [&sequence[0..18], &sequence[9..27], &sequence[18..30]]
        .iter()
        .enumerate()
        .map(|(i, s)| Molecule::dna(format!("read {}", i + 1), s, Topology::Linear).unwrap())
        .collect::<Vec<_>>();
    let mut histories = molecules
        .iter()
        .enumerate()
        .map(|(i, m)| DocumentHistory::new(format!("source-{i}"), m.clone()).unwrap())
        .collect::<Vec<_>>();
    let bindings = histories
        .iter()
        .enumerate()
        .map(|(i, h)| InputEventBinding {
            input_id: format!("read_{}", i + 1),
            event: EventBinding {
                document_id: h.id.clone(),
                event_id: h.head_event().into(),
                revision: h.molecule.revision(),
            },
        })
        .collect::<Vec<_>>();
    let inputs = molecules
        .into_iter()
        .enumerate()
        .map(|(i, m)| ReadInput {
            id: format!("read_{}", i + 1),
            molecule: m,
        })
        .collect::<Vec<_>>();
    let mut options = ContigOptions::default();
    options.trim.min_retained_bases = 4;
    options.min_overlap = 6;
    options.both_strands = false;
    let report = contig::assemble(&inputs, &options).unwrap();
    let record = SavedAssembly::new(report, bindings).unwrap();
    histories[0]
        .commit(
            Molecule::dna("newer edited read head", "GATTACA", Topology::Linear).unwrap(),
            "Edit after assembly",
        )
        .unwrap();
    let mut w =
        Workspace::from_documents("assembly fixture", vec![histories[0].molecule.clone()], 0)
            .unwrap();
    w.documents.clear();
    w.revisions.clear();
    w.history.clear();
    for h in &histories {
        h.contribute(&mut w).unwrap();
    }
    w.active_document = histories[0].id.clone();
    (w, record, inputs)
}
#[test]
fn exact_historical_read_events_and_full_evidence_survive_native_archive() {
    let (mut w, record, original) = fixture();
    let before = w.documents.clone();
    let (hash, attachment) = make_attachment(&record).unwrap();
    w.attachments.insert(hash.clone(), attachment);
    w.validate().unwrap();
    let restored = workspace::load_workspace(&workspace::save_workspace(&w).unwrap()).unwrap();
    let summary = list(&restored).unwrap();
    assert_eq!(summary.len(), 1);
    assert_eq!(summary[0].attachment_hash, hash);
    assert_eq!(summary[0].read_count, 3);
    let saved = load_record(&restored.attachments[&hash].bytes).unwrap();
    let inputs = resolve_inputs(&restored, &saved).unwrap();
    assert_eq!(inputs, original);
    assert_eq!(
        serde_json::to_value(&restored.documents).unwrap(),
        serde_json::to_value(&before).unwrap()
    );
    assert_ne!(
        restored.head(&restored.active_document).unwrap().revision(),
        inputs[0].molecule.revision()
    );
    let replay = contig::assemble(&inputs, &saved.report.options).unwrap();
    assert_eq!(replay, saved.report);
}
#[test]
fn unknown_fields_bad_bindings_and_oversized_records_reject_without_state_installation() {
    let (w, record, _) = fixture();
    let before = workspace::save_workspace(&w).unwrap();
    let mut wrong = record.clone();
    wrong.inputs[0].event.document_id = wrong.inputs[1].event.document_id.clone();
    let (hash, a) = make_attachment(&wrong).unwrap();
    let mut next = w.clone();
    next.attachments.insert(hash, a);
    assert!(next.validate().is_err());
    assert_eq!(workspace::save_workspace(&w).unwrap(), before);
    let mut reversed = record.clone();
    reversed.inputs.reverse();
    assert!(reversed.validate().is_err());
    let mut raw = serde_json::to_value(&record).unwrap();
    raw["report"]["options"]["ignored_policy"] = true.into();
    assert!(load_record(&serde_json::to_vec(&raw).unwrap()).is_err());
    let mut raw = serde_json::to_value(&record).unwrap();
    raw["ignored_authority"] = "verified".into();
    assert!(load_record(&serde_json::to_vec(&raw).unwrap()).is_err());
    assert!(matches!(
        load_record(&vec![b' '; MAX_SAVED_BYTES + 1]),
        Err(collagene_model::Error::TooLarge)
    ));
}
#[test]
fn old_producer_engine_is_historical_evidence_and_does_not_pass_current_binding() {
    let (mut w, mut record, inputs) = fixture();
    record.report.engine_version = "historical-producer/0".into();
    for c in &mut record.report.contigs {
        if let Some(m) = &mut c.molecule {
            m.metadata
                .insert("contig:engine".into(), record.report.engine_version.clone());
        }
    }
    record.validate().unwrap();
    let (hash, a) = make_attachment(&record).unwrap();
    w.attachments.insert(hash.clone(), a);
    w.validate().unwrap();
    let saved = load_record(&w.attachments[&hash].bytes).unwrap();
    assert_eq!(saved.report.engine_version, "historical-producer/0");
    assert!(
        saved
            .report
            .validate_bindings_for(&inputs, &saved.report.options)
            .is_err()
    );
    let replay =
        contig::assemble(&resolve_inputs(&w, &saved).unwrap(), &saved.report.options).unwrap();
    assert_eq!(replay.engine_version, contig::CONTIG_ENGINE_VERSION);
    assert_ne!(replay.engine_version, saved.report.engine_version);
}
