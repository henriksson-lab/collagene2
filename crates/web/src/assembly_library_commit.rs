//! Pure host bridge: typed editable drafts plus documentary native recovery.
//! All scientific receipts remain invalidated by the existing proposal/application route.
use crate::{
    assembly_alternatives,
    assembly_draft_model::{self, Draft},
    assembly_library_topology::{self as topology, Proposal as TopologyProposal},
};
use collagene_core::{
    assembly_library_history::{self as history, Impact, Record},
    combinatorial_design,
    run_lineage::{self, EventBinding},
    workspace::Workspace,
};
pub struct Staged {
    pub workspace: Workspace,
    pub record_hash: String,
    pub after_hash: String,
    pub source_workspace_fingerprint: String,
    pub proposal_revision: String,
    pub draft: Draft,
}
/// Keep recovery and ordinary autosave on the same bounded raw-state codec.
pub fn raw_encode(draft: &Draft) -> Result<Vec<u8>, String> {
    assembly_draft_model::encode_raw(draft)
}
fn pin_key(pin: &EventBinding) -> (&str, &str, &str) {
    (&pin.document_id, &pin.event_id, &pin.revision)
}
/// The full host model has already been bounded/validated. Every original cached
/// molecule is compared to an actual retained native event, not merely equal DNA
/// elsewhere in the workspace. Include inactive/removed choices and working slots.
fn captures(frame: &Workspace, draft: &Draft) -> Result<Vec<EventBinding>, String> {
    let mut pins = Vec::new();
    for capture in assembly_alternatives::all_captures(draft) {
        let pin = &capture.source.binding;
        let exact = combinatorial_design::resolve_event(frame, pin).map_err(|e| e.to_string())?;
        let revision = exact.revision();
        if revision != capture.molecule.revision() || revision != pin.revision {
            return Err("An assembly capture differs from its exact retained source event.".into());
        }
        pins.push(pin.clone());
    }
    Ok(pins)
}
fn normalize(mut pins: Vec<EventBinding>) -> Result<Vec<EventBinding>, String> {
    if pins.len() > history::MAX_SOURCES {
        return Err("Assembly source union exceeds the recovery bound.".into());
    }
    pins.sort_by(|a, b| pin_key(a).cmp(&pin_key(b)));
    pins.dedup_by(|a, b| pin_key(a) == pin_key(b));
    Ok(pins)
}
fn impact(value: &topology::RowImpact) -> Impact {
    Impact {
        construct_id: value.construct_id.clone(),
        name: value.name.clone(),
        primer_review_bins: value.primer_review_bins.clone(),
        removed_preparation_bins: value.removed_preparation_bins.clone(),
        retained_internal_edits: value.retained_internal_edits,
        recovered_internal_edits: value.recovered_internal_edits,
        retained_junctions: value.retained_junctions,
        recovered_junctions: value.recovered_junctions,
        new_junctions: value.new_junctions.clone(),
        overlap_text_changed: value.overlap_text_changed,
        role_order_changed: value.role_order_changed,
        candidate_cleared: value.candidate_cleared,
        construct_removed: value.construct_removed,
        construct_added: value.construct_added,
    }
}
/// Stage both the raw replacement and its exact recovery in a private frame.
/// The UI must recheck its draft/epoch/fingerprint and persist this WHOLE frame in
/// one transaction before installing `draft`; an autosave is not this operation.
pub fn stage(
    frame: &Workspace,
    current: &Draft,
    proposal: &TopologyProposal,
    ack: &str,
    expected_source_fingerprint: &str,
) -> Result<Staged, String> {
    let fingerprint = run_lineage::workspace_fingerprint(frame).map_err(|e| e.to_string())?;
    if fingerprint != expected_source_fingerprint {
        return Err("Workspace changed after the library edit was proposed.".into());
    }
    let applied = topology::apply(current, proposal, ack)?;
    let before = raw_encode(&applied.recovery)?;
    let after = raw_encode(&applied.draft)?;
    let mut sources = captures(frame, &applied.recovery)?;
    sources.extend(captures(frame, &applied.draft)?);
    let sources = normalize(sources)?;
    let record = history::create(
        frame,
        before,
        &after,
        &applied.proposal_revision,
        sources,
        applied.impacts.iter().map(impact).collect(),
    )
    .map_err(|e| e.to_string())?;
    let commit = history::install(frame.clone(), record, after).map_err(|e| e.to_string())?;
    Ok(Staged {
        workspace: commit.workspace,
        record_hash: commit.record_hash,
        after_hash: commit.after_hash,
        source_workspace_fingerprint: commit.source_workspace_fingerprint,
        proposal_revision: applied.proposal_revision,
        draft: applied.draft,
    })
}
/// Historical data only. Full typed decoding and all-capture source attestation
/// happen independently of any native record claims before recovery is offered.
pub fn load_recovery(frame: &Workspace, record_hash: &str) -> Result<(Draft, Record), String> {
    let record = history::load(frame, record_hash).map_err(|e| e.to_string())?;
    let recovered = assembly_draft_model::decode(&record.before_draft)?;
    let pins = normalize(captures(frame, &recovered)?)?;
    if pins.iter().any(|pin| {
        !record
            .sources
            .iter()
            .any(|declared| pin_key(pin) == pin_key(declared))
    }) {
        return Err("Recovery record omits a captured original source binding.".into());
    }
    Ok((recovered, record))
}
/// Restoring is a new explicitly reviewed edit. Preserve the allocation highwater
/// across current, recovered and retained historical records; never reuse identities.
pub fn proposal_restore(
    frame: &Workspace,
    current: &Draft,
    record_hash: &str,
) -> Result<TopologyProposal, String> {
    current.validate()?;
    captures(frame, current)?;
    let (recovered, record) = load_recovery(frame, record_hash)?;
    let highwater =
        history::allocation_highwater(frame, &record.library_id).map_err(|e| e.to_string())?;
    topology::propose_restore(current, &recovered, highwater)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        assembly_draft_model::CapturedSource,
        cloning_methods_model::{Mode, OverlapPreparation},
        source_identity::SourceInfo,
    };
    use collagene_core::workspace::{RevisionEvent, WorkspaceAttachment};
    use collagene_model::{Molecule, Topology};
    use std::sync::Arc;
    fn fixture() -> (Workspace, Draft) {
        let mut frame = Workspace::from_documents(
            "host recovery",
            vec![
                Molecule::dna("backbone", "AACCGGTT", Topology::Linear).unwrap(),
                Molecule::dna("insert", "TTGGAACC", Topology::Linear).unwrap(),
                Molecule::dna("inactive alternative", "GGTTCCAA", Topology::Linear).unwrap(),
            ],
            0,
        )
        .unwrap();
        let mut draft = Draft::default();
        draft.circular = false;
        draft.methods.mode = Mode::Overlap;
        for (i, slot) in draft.slots.iter_mut().enumerate() {
            let document = &frame.documents[i];
            let molecule = frame.revisions[&document.head].clone();
            {
                let collagene_core::assembly_intent::Selector::ReferenceRange {
                    end: value, ..
                } = &mut slot.selector
                else {
                    panic!("fixture requires reference selector")
                };
                *value = molecule.sequence.len();
            }
            slot.captured = Some(Arc::new(CapturedSource {
                source: SourceInfo {
                    binding: EventBinding {
                        document_id: document.id.clone(),
                        event_id: document.head_event.clone().unwrap(),
                        revision: document.head.clone(),
                    },
                    name: molecule.name.clone(),
                    length: molecule.sequence.len(),
                    alphabet: molecule.sequence.alphabet(),
                    topology: molecule.topology,
                },
                molecule,
            }));
        }
        draft.sync_rows();
        draft.methods.overlaps = "4".into();
        for row in &mut draft.methods.rows {
            row.preparation = OverlapPreparation::WholeLinear;
            row.forward_tail = "AAA".into();
        }
        let mut draft = assembly_alternatives::enable(&draft, "host-library").unwrap();
        // Add a valid inactive source choice without changing the active working capture.
        let document = &frame.documents[2];
        let molecule = frame.revisions[&document.head].clone();
        let capture = Arc::new(CapturedSource {
            source: SourceInfo {
                binding: EventBinding {
                    document_id: document.id.clone(),
                    event_id: document.head_event.clone().unwrap(),
                    revision: document.head.clone(),
                },
                name: molecule.name.clone(),
                length: molecule.sequence.len(),
                alphabet: molecule.sequence.alphabet(),
                topology: molecule.topology,
            },
            molecule,
        });
        let library = draft.alternatives.as_mut().unwrap();
        let mut choice = library.bins[1].choices[0].clone();
        choice.id = format!("{}-choice-{}", library.id, library.next_id);
        library.next_id += 1;
        choice.slot.captured = Some(capture);
        library.bins[1].choices.push(choice);
        draft.validate().unwrap();
        let bytes = raw_encode(&draft).unwrap();
        let key = blake3::hash(&bytes).to_hex().to_string();
        frame.attachments.insert(
            key,
            WorkspaceAttachment {
                kind: history::ACTIVE_DRAFT_KIND.into(),
                name: "assembly-draft.json".into(),
                media_type: "application/json".into(),
                bytes,
            },
        );
        (frame, draft)
    }
    fn reorder(draft: &Draft) -> TopologyProposal {
        let ids = draft
            .alternatives
            .as_ref()
            .unwrap()
            .bins
            .iter()
            .rev()
            .map(|b| b.id.clone())
            .collect();
        topology::propose(draft, topology::Change::Reorder { bin_ids: ids }).unwrap()
    }
    #[test]
    fn removed_and_retained_internal_edit_counts_survive_native_history_archive() {
        let (frame, draft) = fixture();
        let draft = crate::assembly_edit_model::declare(&draft, 0, "3", "5", "TAG", None).unwrap();
        let draft = crate::assembly_substitution_model::declare(&draft, 1, "1", "A", None).unwrap();
        let removed = draft.alternatives.as_ref().unwrap().bins[0].id.clone();
        let proposal =
            topology::propose(&draft, topology::Change::Remove { bin_id: removed }).unwrap();
        let fp = run_lineage::workspace_fingerprint(&frame).unwrap();
        let staged = stage(&frame, &draft, &proposal, proposal.revision(), &fp).unwrap();
        let archive = collagene_core::workspace::save_workspace(&staged.workspace).unwrap();
        let reopened = collagene_core::workspace::load_workspace(&archive).unwrap();
        let (recovered, record) = load_recovery(&reopened, &staged.record_hash).unwrap();
        assert!(recovered == draft);
        assert_eq!(record.schema_version, history::SCHEMA_VERSION);
        assert_eq!(record.impacts[0].retained_internal_edits, 1);
        assert_eq!(record.impacts[0].recovered_internal_edits, 1);
        assert!(staged.draft.edits.is_empty());
        assert_eq!(staged.draft.substitutions.len(), 1);
        assert_eq!(
            reopened.attachments[&staged.record_hash].kind,
            history::KIND
        );
    }
    #[test]
    fn raw_unsaved_fields_and_inactive_captures_survive_crosslayer_commit() {
        let (frame, mut draft) = fixture();
        let frame_snapshot = serde_json::to_vec(&frame).unwrap();
        let saved = draft.alternatives.as_ref().unwrap().rows[0]
            .edits
            .methods
            .rows[0]
            .forward_tail
            .clone();
        draft.methods.rows[0].forward_tail = "CCCGGG".into();
        draft.methods.overlaps = "4,,unfinished".into();
        let before = raw_encode(&draft).unwrap();
        assert_ne!(before, assembly_draft_model::encode(&draft).unwrap());
        let proposal = reorder(&draft);
        let fp = run_lineage::workspace_fingerprint(&frame).unwrap();
        let staged = stage(&frame, &draft, &proposal, proposal.revision(), &fp).unwrap();
        assert_eq!(serde_json::to_vec(&frame).unwrap(), frame_snapshot);
        let (recovered, record) = load_recovery(&staged.workspace, &staged.record_hash).unwrap();
        assert_eq!(record.before_draft, before);
        assert_eq!(recovered.methods.overlaps, "4,,unfinished");
        assert_eq!(recovered.methods.rows[0].forward_tail, "CCCGGG");
        assert_eq!(
            recovered.alternatives.as_ref().unwrap().rows[0]
                .edits
                .methods
                .rows[0]
                .forward_tail,
            saved
        );
        assert_eq!(
            record.sources.len(),
            3,
            "inactive choice remains an attested original source"
        );
        assert_eq!(
            staged.workspace.attachments[&staged.after_hash].bytes,
            raw_encode(&staged.draft).unwrap()
        );
        assert_eq!(
            staged
                .workspace
                .attachments
                .values()
                .filter(|a| a.kind == history::ACTIVE_DRAFT_KIND)
                .count(),
            1
        );
    }
    #[test]
    fn changed_frame_draft_or_acknowledgement_never_mutates_live_inputs() {
        let (frame, draft) = fixture();
        let proposal = reorder(&draft);
        let fp = run_lineage::workspace_fingerprint(&frame).unwrap();
        let original = raw_encode(&draft).unwrap();
        let snapshot = serde_json::to_vec(&frame).unwrap();
        let mut changed = frame.clone();
        changed.name = "concurrent native edit".into();
        assert!(stage(&changed, &draft, &proposal, proposal.revision(), &fp).is_err());
        let mut changed = draft.clone();
        changed.name = "concurrent working edit".into();
        assert!(stage(&frame, &changed, &proposal, proposal.revision(), &fp).is_err());
        assert!(stage(&frame, &draft, &proposal, &"0".repeat(64), &fp).is_err());
        assert_eq!(raw_encode(&draft).unwrap(), original);
        assert_eq!(serde_json::to_vec(&frame).unwrap(), snapshot);
    }
    #[test]
    fn exact_native_membership_includes_inactive_choices_and_cached_snapshot() {
        let (frame, mut draft) = fixture();
        let library = draft.alternatives.as_mut().unwrap();
        let mut capture = library.bins[1].choices[1]
            .slot
            .captured
            .as_deref()
            .unwrap()
            .clone();
        capture.source.binding.document_id = "foreign-document".into();
        library.bins[1].choices[1].slot.captured = Some(Arc::new(capture));
        let proposal = reorder(&draft);
        let fp = run_lineage::workspace_fingerprint(&frame).unwrap();
        assert!(stage(&frame, &draft, &proposal, proposal.revision(), &fp).is_err());
        let (_, mut damaged) = fixture();
        let choice = &mut damaged.alternatives.as_mut().unwrap().bins[1].choices[1];
        let mut capture = choice.slot.captured.as_deref().unwrap().clone();
        capture.molecule =
            Molecule::dna("inactive alternative", "AACCAACC", Topology::Linear).unwrap();
        choice.slot.captured = Some(Arc::new(capture));
        assert!(
            raw_encode(&damaged).is_err(),
            "cached source bytes cannot diverge from their pinned molecular revision"
        );
        let (frame, draft) = fixture();
        let proposal = reorder(&draft);
        let fp = run_lineage::workspace_fingerprint(&frame).unwrap();
        let staged = stage(&frame, &draft, &proposal, proposal.revision(), &fp).unwrap();
        let mut record = history::load(&staged.workspace, &staged.record_hash).unwrap();
        record.sources.clear();
        let bytes = history::encode_record(&record).unwrap();
        let key = blake3::hash(&bytes).to_hex().to_string();
        let mut altered = staged.workspace;
        altered.attachments.insert(
            key.clone(),
            WorkspaceAttachment {
                kind: history::KIND.into(),
                name: "omitted-pins".into(),
                media_type: history::MEDIA_TYPE.into(),
                bytes,
            },
        );
        assert!(
            load_recovery(&altered, &key).is_err(),
            "native opaque validity cannot hide omitted host captures"
        );
    }
    #[test]
    fn historical_heads_remain_loadable_and_restore_is_a_new_monotonic_record() {
        let (frame, mut draft) = fixture();
        draft.methods.overlaps = "unfinished raw review".into();
        let proposal = reorder(&draft);
        let fp = run_lineage::workspace_fingerprint(&frame).unwrap();
        let staged = stage(&frame, &draft, &proposal, proposal.revision(), &fp).unwrap();
        let mut frame = staged.workspace;
        let old = frame.documents[0].clone();
        let molecule = Molecule::dna("changed later", "GGTTAACC", Topology::Linear).unwrap();
        let revision = molecule.revision();
        let event = RevisionEvent::new(
            old.id.clone(),
            revision.clone(),
            vec![old.head_event.clone().unwrap()],
            "later source edit",
        )
        .unwrap();
        frame.revisions.insert(revision.clone(), molecule);
        frame.history.push(event.clone());
        let doc = &mut frame.documents[0];
        doc.undo.push(doc.head.clone());
        doc.undo_events.push(doc.head_event.clone().unwrap());
        doc.head = revision;
        doc.head_event = Some(event.id);
        load_recovery(&frame, &staged.record_hash).unwrap();
        let mut current = staged.draft;
        current.alternatives.as_mut().unwrap().next_id += 50;
        let next = current.alternatives.as_ref().unwrap().next_id;
        let restore = proposal_restore(&frame, &current, &staged.record_hash).unwrap();
        assert!(restore.after().alternatives.as_ref().unwrap().next_id >= next);
        let fp = run_lineage::workspace_fingerprint(&frame).unwrap();
        let restored = stage(&frame, &current, &restore, restore.revision(), &fp).unwrap();
        assert_eq!(restored.draft.methods.overlaps, "unfinished raw review");
        assert_eq!(history::list(&restored.workspace).unwrap().len(), 2);
        assert!(
            history::allocation_highwater(&restored.workspace, "host-library").unwrap() >= next
        );
        assert_eq!(
            restored.workspace.documents[0].head, frame.documents[0].head,
            "restoration never rewinds native source heads"
        );
    }
}
