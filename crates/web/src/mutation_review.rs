//! Install declared reference intent with exact computational ancestry. No physical
//! product or positive scientific result is installed by this operation.
use collagene_core::{
    history::DocumentHistory,
    mutation_intent::{self, Report},
    run_lineage::EventBinding,
    workspace::{Workspace, WorkspaceAttachment},
};

pub fn stage(
    frame: &Workspace,
    source: &EventBinding,
    report: &Report,
) -> Result<Workspace, String> {
    frame.validate().map_err(|e| e.to_string())?;
    let docs = DocumentHistory::restore_all(frame).map_err(|e| e.to_string())?;
    let original = docs
        .iter()
        .find(|doc| doc.id == source.document_id)
        .ok_or("The intended edit's original source document is unavailable.")?;
    let molecule = original
        .revision_molecule(&source.event_id)
        .map_err(|e| e.to_string())?;
    if molecule.revision() != source.revision || source.revision != report.source_revision {
        return Err("The intended edit differs from its exact original source event.".into());
    }
    let fresh = mutation_intent::build(&molecule, &report.request).map_err(|e| e.to_string())?;
    report
        .validate_binding(&fresh.context_revision)
        .map_err(|e| e.to_string())?;
    if fresh.target.revision() != report.target.revision()
        || fresh.deleted != report.deleted
        || fresh.origin_shift != report.origin_shift
    {
        return Err("The intended reference differs from the declared source edit.".into());
    }
    let revision = fresh.target.revision();
    let existing = docs.iter().find(|doc| {
        doc.head_revision() == revision
            && doc
                .summaries()
                .iter()
                .any(|event| event.current && event.parents == [source.event_id.clone()])
    });
    let mut staged = frame.clone();
    let target = if let Some(existing) = existing {
        existing.id.clone()
    } else {
        if docs.len() >= collagene_core::workspace::MAX_DOCUMENTS {
            return Err("The workspace document limit has been reached.".into());
        }
        let prefix = format!("mutation-intent-{}", &fresh.context_revision[..16]);
        let id = (1..=1001)
            .map(|n| format!("{prefix}-{n}"))
            .find(|id| !docs.iter().any(|doc| doc.id == *id))
            .ok_or("No intended-reference identity is available.")?;
        let derived = original.derive_event(&source.event_id, id,
            fresh.target.clone(),"Declared mutation intent from an exact source event; reference only, preparation and experimental outcome unverified")
            .map_err(|e|e.to_string())?;
        let id = derived.id.clone();
        derived.contribute(&mut staged).map_err(|e| e.to_string())?;
        id
    };
    let bytes =
        serde_json::to_vec(&serde_json::json!({"version":1,"source":source,"report":fresh}))
            .map_err(|e| e.to_string())?;
    if bytes.len() > collagene_core::workspace::MAX_ENTRY_BYTES {
        return Err("Intended-edit evidence exceeds the native attachment limit.".into());
    }
    let hash = blake3::hash(&bytes).to_hex().to_string();
    staged
        .attachments
        .entry(hash)
        .or_insert(WorkspaceAttachment {
            kind: "mutation-intent-v1".into(),
            name: "mutation-intent.json".into(),
            media_type: "application/json".into(),
            bytes,
        });
    staged.active_document = target;
    staged.validate().map_err(|e| e.to_string())?;
    Ok(staged)
}

#[cfg(test)]
mod tests {
    use super::*;
    use collagene_model::{Molecule, Topology};
    fn fixture() -> (Workspace, EventBinding, Report) {
        let source = Molecule::dna("original", "ACGTCAGT", Topology::Circular).unwrap();
        let frame = Workspace::from_documents("mutation", vec![source.clone()], 0).unwrap();
        let doc = &frame.documents[0];
        let event = EventBinding {
            document_id: doc.id.clone(),
            event_id: doc.head_event.clone().unwrap(),
            revision: doc.head.clone(),
        };
        let report = mutation_intent::build(
            &source,
            &mutation_intent::Request {
                name: "intended mutant".into(),
                edit_start: 2,
                edit_end: 4,
                inserted: "AA".into(),
            },
        )
        .unwrap();
        (frame, event, report)
    }
    #[test]
    fn exact_source_event_parent_and_original_are_preserved_and_repeated_creation_reuses_target() {
        let (frame, source, report) = fixture();
        let original = serde_json::to_value(&frame).unwrap();
        let staged = stage(&frame, &source, &report).unwrap();
        assert_eq!(staged.documents.len(), 2);
        assert_eq!(
            staged
                .head(&staged.active_document)
                .unwrap()
                .sequence
                .as_str(),
            "ACAACAGT"
        );
        let target = staged
            .documents
            .iter()
            .find(|doc| doc.id == staged.active_document)
            .unwrap();
        let event = staged
            .history
            .iter()
            .find(|e| Some(&e.id) == target.head_event.as_ref())
            .unwrap();
        assert_eq!(event.parents, [source.event_id.clone()]);
        assert_eq!(staged.documents[0].head, frame.documents[0].head);
        assert_eq!(
            staged.documents[0].head_event,
            frame.documents[0].head_event
        );
        assert_eq!(serde_json::to_value(&frame).unwrap(), original);
        assert_eq!(
            serde_json::to_value(stage(&staged, &source, &report).unwrap()).unwrap(),
            serde_json::to_value(&staged).unwrap()
        );
        assert_eq!(
            staged
                .attachments
                .values()
                .filter(|a| a.kind == "mutation-intent-v1")
                .count(),
            1
        );
    }
    #[test]
    fn forged_target_or_wrong_event_cannot_install_and_equal_dna_alias_keeps_distinct_ancestry() {
        let (frame, source, report) = fixture();
        let mut forged = report.clone();
        forged.target.sequence = collagene_model::Sequence::dna("AAAAAAAA").unwrap();
        assert!(stage(&frame, &source, &forged).is_err());
        let mut wrong = source.clone();
        wrong.event_id = format!("event-{}", "0".repeat(64));
        assert!(stage(&frame, &wrong, &report).is_err());
        let mut docs = DocumentHistory::restore_all(&frame).unwrap();
        let alias = docs[0]
            .fork_event(&source.event_id, "equal-dna-source")
            .unwrap();
        let other = EventBinding {
            document_id: alias.id.clone(),
            event_id: alias.head_event().into(),
            revision: alias.head_revision().into(),
        };
        docs.push(alias);
        let mut frame = frame.clone();
        docs[1].contribute(&mut frame).unwrap();
        frame.validate().unwrap();
        let first = stage(&frame, &source, &report).unwrap();
        let second = stage(&first, &other, &report).unwrap();
        assert_ne!(first.active_document, second.active_document);
        let target = second
            .documents
            .iter()
            .find(|d| d.id == second.active_document)
            .unwrap();
        assert_eq!(
            second
                .history
                .iter()
                .find(|e| Some(&e.id) == target.head_event.as_ref())
                .unwrap()
                .parents,
            [other.event_id]
        );
    }
}
