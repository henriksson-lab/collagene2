//! Narrow host installation guard for named-design attachment transactions.
//! This never replaces documents, history, candidates or the active working draft.
use crate::assembly_saved_designs_model::KIND;
use collagene_core::{
    run_lineage,
    workspace::{Workspace, WorkspaceAttachment},
};
use serde::Serialize;
use std::io::Write;

struct HashWriter {
    hash: blake3::Hasher,
    remaining: usize,
}
impl Write for HashWriter {
    fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
        if bytes.len() > self.remaining {
            return Err(std::io::Error::other("saved-design transaction size"));
        }
        self.remaining -= bytes.len();
        self.hash.update(bytes);
        Ok(bytes.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}
fn hash(value: &impl Serialize) -> Result<String, String> {
    let mut writer = HashWriter {
        hash: blake3::Hasher::new(),
        remaining: collagene_core::workspace::MAX_EXPANDED_BYTES,
    };
    serde_json::to_writer(&mut writer, value).map_err(|e| e.to_string())?;
    Ok(writer.hash.finalize().to_hex().to_string())
}
/// Content descriptors are sufficient only after native Workspace::validate:
/// revision keys bind complete molecules and candidate keys bind complete plans.
/// Do not expand or clone the rich models merely to compare unchanged state.
fn source_state(frame: &Workspace) -> Result<String, String> {
    hash(&(
        "saved-design-unchanged-workspace-v1",
        &frame.name,
        &frame.active_document,
        &frame.documents,
        frame.revisions.keys().collect::<Vec<_>>(),
        &frame.history,
        frame
            .candidates
            .iter()
            .map(|candidate| &candidate.plan_revision)
            .collect::<Vec<_>>(),
        &frame.active_candidate,
    ))
}
fn same_attachment(a: &WorkspaceAttachment, b: &WorkspaceAttachment) -> bool {
    a.kind == b.kind && a.name == b.name && a.media_type == b.media_type && a.bytes == b.bytes
}
/// Revalidate immediately before installing `after.attachments` only. A valid
/// named payload is editable state, never a scientific receipt. Deliberately do
/// not require hydration here: malformed/orphaned named attachments may be removed.
pub(crate) fn validate(
    before: &Workspace,
    after: &Workspace,
    current: &Workspace,
) -> Result<(), String> {
    let expected = run_lineage::workspace_fingerprint(before).map_err(|e| e.to_string())?;
    let actual = run_lineage::workspace_fingerprint(current).map_err(|e| e.to_string())?;
    if expected != actual {
        return Err(
            "The workspace changed before the saved-design update. Try the action again.".into(),
        );
    }
    after.validate().map_err(|e| e.to_string())?;
    if source_state(before)? != source_state(after)? {
        return Err("A saved-design update cannot change documents, history, candidates or workspace selection.".into());
    }
    for (key, value) in before
        .attachments
        .iter()
        .filter(|(_, value)| value.kind != KIND)
    {
        if !after
            .attachments
            .get(key)
            .is_some_and(|after| same_attachment(value, after))
        {
            return Err(
                "A saved-design update cannot remove or alter another workspace attachment.".into(),
            );
        }
    }
    for (key, value) in after
        .attachments
        .iter()
        .filter(|(_, value)| value.kind != KIND)
    {
        if !before
            .attachments
            .get(key)
            .is_some_and(|before| same_attachment(value, before))
        {
            return Err(
                "A saved-design update cannot add an unrelated workspace attachment.".into(),
            );
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        assembly_draft_model::{CapturedSource, Draft},
        assembly_saved_designs_model as designs,
        source_identity::SourceInfo,
    };
    use collagene_core::{
        Molecule, Topology,
        run_lineage::EventBinding,
        transfer::Candidate,
        verification::{DesignPlan, Operation, PlanInput, Source},
        workspace::RevisionEvent,
    };
    use std::sync::Arc;
    fn attachment(frame: &mut Workspace, kind: &str, name: &str, bytes: &[u8]) -> String {
        let key = blake3::hash(bytes).to_hex().to_string();
        frame.attachments.insert(
            key.clone(),
            WorkspaceAttachment {
                kind: kind.into(),
                name: name.into(),
                media_type: "application/json".into(),
                bytes: bytes.to_vec(),
            },
        );
        key
    }
    fn fixture() -> (Workspace, Draft) {
        let source = Molecule::dna("Original", "AACCGGTT", Topology::Linear).unwrap();
        let mut frame =
            Workspace::from_documents("Exact sources", vec![source.clone()], 0).unwrap();
        let doc = &frame.documents[0];
        let capture = Arc::new(CapturedSource {
            source: SourceInfo {
                binding: EventBinding {
                    document_id: doc.id.clone(),
                    event_id: doc.head_event.clone().unwrap(),
                    revision: doc.head.clone(),
                },
                name: source.name.clone(),
                length: source.sequence.len(),
                alphabet: source.sequence.alphabet(),
                topology: source.topology,
            },
            molecule: source.clone(),
        });
        let mut draft = Draft::default();
        for slot in &mut draft.slots {
            slot.captured = Some(capture.clone());
            slot.start = 0;
            slot.end = 8;
        }
        draft.sync_rows();
        draft.validate().unwrap();
        let raw = crate::assembly_draft_model::encode_raw(&draft).unwrap();
        attachment(
            &mut frame,
            "assembly-draft-v3",
            "Current working draft",
            &raw,
        );
        attachment(
            &mut frame,
            "user-note-v1",
            "Unrelated note",
            br#"{"keep":"untouched"}"#,
        );
        frame.validate().unwrap();
        (frame, draft)
    }
    fn assert_rejected(before: &Workspace, after: &Workspace, current: &Workspace) {
        let originals = [before, after, current].map(|w| serde_json::to_value(w).unwrap());
        assert!(validate(before, after, current).is_err());
        for (frame, value) in [before, after, current].into_iter().zip(originals) {
            assert_eq!(serde_json::to_value(frame).unwrap(), value);
        }
    }
    #[test]
    fn actual_named_save_rename_and_delete_preserve_autosave_foreign_data_and_sources() {
        let (frame, draft) = fixture();
        let original = serde_json::to_value(&frame).unwrap();
        let (saved, entry) = designs::save(&frame, "First design", &draft).unwrap();
        validate(&frame, &saved, &frame).unwrap();
        assert!(designs::load(&saved, &entry.id).unwrap() == draft);
        let renamed = designs::rename(&saved, &entry.id, "Renamed design").unwrap();
        validate(&saved, &renamed, &saved).unwrap();
        assert!(designs::load(&renamed, &entry.id).unwrap() == draft);
        let removed = designs::delete(&renamed, &entry.id).unwrap();
        validate(&renamed, &removed, &renamed).unwrap();
        assert_eq!(serde_json::to_value(&removed).unwrap(), original);
        assert_eq!(serde_json::to_value(&frame).unwrap(), original);
    }
    #[test]
    fn stale_parent_frame_cannot_discard_newer_working_draft_note_or_source_state() {
        let (frame, draft) = fixture();
        let (saved, _) = designs::save(&frame, "Bound save", &draft).unwrap();
        let mut current = frame.clone();
        attachment(
            &mut current,
            "user-note-v1",
            "Newer note",
            br#"{"newer":true}"#,
        );
        assert_rejected(&frame, &saved, &current);
        let mut current = frame.clone();
        let key = current
            .attachments
            .iter()
            .find(|(_, a)| a.kind == "assembly-draft-v3")
            .unwrap()
            .0
            .clone();
        current.attachments.get_mut(&key).unwrap().name = "Newer working draft label".into();
        current.validate().unwrap();
        assert_rejected(&frame, &saved, &current);
        let mut current = frame.clone();
        current.name = "Newer workspace name".into();
        assert_rejected(&frame, &saved, &current);
    }
    #[test]
    fn source_history_document_revision_and_candidate_changes_are_not_attachment_updates() {
        let (frame, draft) = fixture();
        let (saved, _) = designs::save(&frame, "Bound save", &draft).unwrap();
        let mut changed = saved.clone();
        changed.name = "Another workspace".into();
        assert_rejected(&frame, &changed, &frame);
        let mut changed = saved.clone();
        changed.documents[0].redo = vec![changed.documents[0].head.clone()];
        changed.documents[0].redo_events = vec![changed.documents[0].head_event.clone().unwrap()];
        changed.validate().unwrap();
        assert_rejected(&frame, &changed, &frame);
        let mut changed = saved.clone();
        let event = RevisionEvent::new(
            changed.documents[0].id.clone(),
            changed.documents[0].head.clone(),
            vec![changed.documents[0].head_event.clone().unwrap()],
            "New descriptive history",
        )
        .unwrap();
        changed.history.push(event);
        changed.validate().unwrap();
        assert_rejected(&frame, &changed, &frame);
        let mut changed = saved.clone();
        let extra = Molecule::dna("Unrelated new revision", "AATTCCGG", Topology::Linear).unwrap();
        changed.revisions.insert(extra.revision(), extra);
        changed.validate().unwrap();
        assert_rejected(&frame, &changed, &frame);
        let mut changed = saved;
        let source = changed.revisions[&changed.documents[0].head].clone();
        let candidate = Candidate::new(DesignPlan {
            name: "New unrelated candidate".into(),
            target_revision: source.revision(),
            inputs: vec![PlanInput {
                id: "source".into(),
                revision: source.revision(),
                molecule: source,
            }],
            steps: vec![Operation::TreatEnds {
                template: Source::Input("source".into()),
                treatment: collagene_core::chemistry::EndTreatment::Phosphorylate,
            }],
        });
        candidate.validate().unwrap();
        changed.active_candidate = Some(candidate.plan_revision.clone());
        changed.candidates.push(candidate);
        changed.validate().unwrap();
        assert_rejected(&frame, &changed, &frame);
    }
    #[test]
    fn unrelated_attachment_mutation_and_kind_relabeling_are_refused_but_bad_design_cleanup_is_narrow()
     {
        let (frame, draft) = fixture();
        let (saved, _) = designs::save(&frame, "Bound save", &draft).unwrap();
        let foreign = frame
            .attachments
            .iter()
            .find(|(_, a)| a.kind == "user-note-v1")
            .unwrap()
            .0
            .clone();
        let mut changed = saved.clone();
        changed.attachments.remove(&foreign);
        assert_rejected(&frame, &changed, &frame);
        let mut changed = saved.clone();
        changed.attachments.get_mut(&foreign).unwrap().name = "Different note".into();
        assert_rejected(&frame, &changed, &frame);
        let mut changed = saved.clone();
        changed.attachments.get_mut(&foreign).unwrap().kind = KIND.into();
        assert_rejected(&frame, &changed, &frame);
        let mut changed = saved;
        attachment(
            &mut changed,
            "user-note-v1",
            "Injected unrelated note",
            br#"{"injected":true}"#,
        );
        assert_rejected(&frame, &changed, &frame);
        let mut malformed = frame.clone();
        let key = attachment(
            &mut malformed,
            KIND,
            "Unreadable saved design",
            br#"{"unsupported":"payload"}"#,
        );
        malformed.validate().unwrap();
        let mut cleaned = malformed.clone();
        cleaned.attachments.remove(&key);
        validate(&malformed, &cleaned, &malformed).unwrap();
        assert_eq!(
            serde_json::to_value(&cleaned).unwrap(),
            serde_json::to_value(&frame).unwrap()
        );
    }
}
