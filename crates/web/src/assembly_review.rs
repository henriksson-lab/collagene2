//! Atomic installation of an intended reference and its proposed preparation.
//! Exact source events are checked; no scientific success is restored or inferred.
use collagene_core::{
    assembly_intent::{self, Report},
    history::DocumentHistory,
    run_lineage::{EventBinding, InputEventBinding},
    transfer::Candidate,
    workspace::Workspace,
};
/// Independently declared reference intent, never a target copied from a simulator.
#[derive(Clone)]
pub enum IntendedReference {
    Regions(assembly_intent::Report),
    Junctions(collagene_core::assembly_junction_edit::Intent),
}
impl IntendedReference {
    pub fn target(&self) -> &collagene_core::Molecule {
        match self {
            Self::Regions(intent) => &intent.target,
            Self::Junctions(intent) => &intent.target,
        }
    }
    pub fn fragments(&self) -> &[assembly_intent::PlacedFragment] {
        match self {
            Self::Regions(intent) => &intent.fragments,
            Self::Junctions(intent) => &intent.fragments,
        }
    }
}
pub fn stage_any(
    frame: &Workspace,
    intent: &IntendedReference,
    sources: &[EventBinding],
    candidate: &Candidate,
) -> Result<Staged, String> {
    match intent {
        IntendedReference::Regions(intent) => stage(frame, intent, sources, candidate),
        IntendedReference::Junctions(intent) => stage_edited(frame, intent, sources, candidate),
    }
}
pub struct Staged {
    pub workspace: Workspace,
    pub target: EventBinding,
    pub sources: Vec<InputEventBinding>,
}
/// Stage an existing independently specified target and explicitly pinned inputs.
/// No target is constructed, source head rewound or scientific check restored.
/// Historical source events are valid only when supplied with their exact IDs;
/// equal content in another document is never used as a replacement.
pub fn stage_existing(
    frame: &Workspace,
    candidate: &Candidate,
    target: &EventBinding,
    sources: &[InputEventBinding],
) -> Result<Staged, String> {
    frame.validate().map_err(|e| e.to_string())?;
    candidate.validate().map_err(|e| e.to_string())?;
    if sources.len() != candidate.plan.inputs.len() || sources.len() > 128 {
        return Err("Pinned input count differs from the candidate output profile.".into());
    }
    let docs = DocumentHistory::restore_all(frame).map_err(|e| e.to_string())?;
    let intended = docs
        .iter()
        .find(|doc| doc.id == frame.active_document)
        .ok_or("No active independent target document is retained.")?;
    if intended.id != target.document_id
        || intended.head_event() != target.event_id
        || intended.head_revision() != target.revision
        || target.revision != candidate.plan.target_revision
    {
        return Err(
            "The exact active intended-target event changed; review the plan again.".into(),
        );
    }
    for (input, pinned) in candidate.plan.inputs.iter().zip(sources) {
        if input.id != pinned.input_id {
            return Err("Ordered pinned input identities differ from the candidate.".into());
        }
        let source = docs
            .iter()
            .find(|doc| doc.id == pinned.event.document_id)
            .ok_or("A pinned source document is unavailable.")?;
        let exact = source
            .revision_molecule(&pinned.event.event_id)
            .map_err(|e| e.to_string())?;
        if exact.revision() != pinned.event.revision
            || pinned.event.revision != input.revision
            || exact != input.molecule
        {
            return Err("An exact source event differs from its candidate input.".into());
        }
    }
    let mut staged = frame.clone();
    if !staged
        .candidates
        .iter()
        .any(|existing| existing.plan_revision == candidate.plan_revision)
    {
        staged.candidates.push(candidate.clone());
    }
    staged.active_candidate = Some(candidate.plan_revision.clone());
    staged.validate().map_err(|e| e.to_string())?;
    Ok(Staged {
        workspace: staged,
        target: target.clone(),
        sources: sources.to_vec(),
    })
}
pub fn retained_inputs(
    docs: &[DocumentHistory],
    sources: &[InputEventBinding],
    candidate: &Candidate,
) -> Result<Vec<collagene_core::Molecule>, String> {
    if sources.len() != candidate.plan.inputs.len() {
        return Err("Pinned input count differs from the candidate.".into());
    }
    candidate
        .plan
        .inputs
        .iter()
        .map(|input| {
            let pinned = sources
                .iter()
                .find(|source| source.input_id == input.id)
                .ok_or("A pinned input event is missing.")?;
            let doc = docs
                .iter()
                .find(|d| d.id == pinned.event.document_id)
                .ok_or("The pinned source document is unavailable.")?;
            let molecule = doc
                .revision_molecule(&pinned.event.event_id)
                .map_err(|e| e.to_string())?;
            if molecule.revision() != pinned.event.revision
                || pinned.event.revision != input.revision
            {
                return Err("Pinned source content differs from the preparation input.".into());
            }
            Ok(molecule)
        })
        .collect()
}
pub fn stage(
    frame: &Workspace,
    intent: &Report,
    sources: &[EventBinding],
    candidate: &Candidate,
) -> Result<Staged, String> {
    frame.validate().map_err(|e| e.to_string())?;
    candidate.validate().map_err(|e| e.to_string())?;
    if sources.len() != candidate.plan.inputs.len() || sources.is_empty() || sources.len() > 32 {
        return Err("Assembly source bindings do not match the candidate.".into());
    }
    let docs = DocumentHistory::restore_all(frame).map_err(|e| e.to_string())?;
    let mut parents = Vec::with_capacity(sources.len());
    for (binding, input) in sources.iter().zip(&candidate.plan.inputs) {
        let doc = docs
            .iter()
            .find(|d| d.id == binding.document_id)
            .ok_or("An assembly source document is unavailable.")?;
        let exact = doc
            .revision_molecule(&binding.event_id)
            .map_err(|e| e.to_string())?;
        if exact.revision() != binding.revision
            || binding.revision != input.revision
            || exact.revision() != input.molecule.revision()
        {
            return Err("A pinned assembly source does not match its preparation input.".into());
        }
        parents.push((doc, binding.event_id.as_str()));
    }
    let fresh = assembly_intent::build(&candidate.plan.inputs, &intent.request)
        .map_err(|e| e.to_string())?;
    intent
        .validate_binding(&fresh.context_revision)
        .map_err(|e| e.to_string())?;
    if intent.target.revision() != fresh.target.revision()
        || candidate.plan.target_revision != fresh.target.revision()
    {
        return Err(
            "The candidate does not target the independently declared source regions.".into(),
        );
    }
    let expected_target_revision = fresh.target.revision();
    let expected_parents = sources
        .iter()
        .map(|s| s.event_id.as_str())
        .collect::<std::collections::BTreeSet<_>>();
    let existing = docs.iter().find(|doc| {
        doc.head_revision() == expected_target_revision
            && doc
                .summaries()
                .iter()
                .find(|event| event.current)
                .is_some_and(|event| {
                    event
                        .parents
                        .iter()
                        .map(String::as_str)
                        .collect::<std::collections::BTreeSet<_>>()
                        == expected_parents
                })
    });
    let mut staged = frame.clone();
    let target = if let Some(existing) = existing {
        EventBinding {
            document_id: existing.id.clone(),
            event_id: existing.head_event().into(),
            revision: existing.head_revision().into(),
        }
    } else {
        if docs.len() >= 1000 {
            return Err("The workspace document limit has been reached.".into());
        }
        let prefix = format!("assembly-intent-{}", &fresh.context_revision[..16]);
        let id = (1..=1001)
            .map(|n| format!("{prefix}-{n}"))
            .find(|id| !docs.iter().any(|d| d.id == *id))
            .ok_or("No assembly document identity is available.")?;
        let derived=DocumentHistory::derive_from_events(&parents,id,fresh.target,"Declared assembly design intent from pinned source regions; expected reference only, preparation and experimental outcome unverified").map_err(|e|e.to_string())?;
        let target = EventBinding {
            document_id: derived.id.clone(),
            event_id: derived.head_event().into(),
            revision: derived.head_revision().into(),
        };
        derived.contribute(&mut staged).map_err(|e| e.to_string())?;
        target
    };
    staged.active_document = target.document_id.clone();
    if !staged
        .candidates
        .iter()
        .any(|c| c.plan_revision == candidate.plan_revision)
    {
        staged.candidates.push(candidate.clone());
    }
    staged.active_candidate = Some(candidate.plan_revision.clone());
    staged.validate().map_err(|e| e.to_string())?;
    Ok(Staged {
        workspace: staged,
        target,
        sources: sources
            .iter()
            .zip(&candidate.plan.inputs)
            .map(|(event, input)| InputEventBinding {
                input_id: input.id.clone(),
                event: event.clone(),
            })
            .collect(),
    })
}

fn stage_edited(
    frame: &Workspace,
    intent: &collagene_core::assembly_junction_edit::Intent,
    sources: &[EventBinding],
    candidate: &Candidate,
) -> Result<Staged, String> {
    if candidate.plan.inputs.is_empty() || candidate.plan.inputs.len() > 32 {
        return Err("Edited assembly inputs exceed the bounded staging profile.".into());
    }
    collagene_core::assembly_junction_edit::intent_preflight(
        &candidate
            .plan
            .inputs
            .iter()
            .map(|input| (input.id.as_str(), &input.molecule))
            .collect::<Vec<_>>(),
        &intent.base_request,
        &intent.edits,
    )
    .map_err(|e| e.to_string())?;
    frame.validate().map_err(|e| e.to_string())?;
    candidate.validate().map_err(|e| e.to_string())?;
    if sources.len() != candidate.plan.inputs.len() || sources.is_empty() || sources.len() > 32 {
        return Err("Assembly source bindings do not match the candidate.".into());
    }
    let docs = DocumentHistory::restore_all(frame).map_err(|e| e.to_string())?;
    let mut parents = Vec::with_capacity(sources.len());
    for (binding, input) in sources.iter().zip(&candidate.plan.inputs) {
        let doc = docs
            .iter()
            .find(|d| d.id == binding.document_id)
            .ok_or("An assembly source document is unavailable.")?;
        let exact = doc
            .revision_molecule(&binding.event_id)
            .map_err(|e| e.to_string())?;
        if exact.revision() != binding.revision
            || binding.revision != input.revision
            || exact.revision() != input.molecule.revision()
        {
            return Err("A pinned assembly source does not match its preparation input.".into());
        }
        parents.push((doc, binding.event_id.as_str()));
    }
    let fresh = collagene_core::assembly_junction_edit::build_intent(
        &candidate.plan.inputs,
        &intent.base_request,
        &intent.edits,
    )
    .map_err(|e| e.to_string())?;
    intent
        .validate_binding(&fresh.context_revision)
        .map_err(|e| e.to_string())?;
    if intent.target.revision() != fresh.target.revision()
        || candidate.plan.target_revision != fresh.target.revision()
    {
        return Err(
            "The candidate does not target the independently declared source regions and junction edits.".into(),
        );
    }
    let expected_target_revision = fresh.target.revision();
    let expected_parents = sources
        .iter()
        .map(|s| s.event_id.as_str())
        .collect::<std::collections::BTreeSet<_>>();
    let existing = docs.iter().find(|doc| {
        doc.head_revision() == expected_target_revision
            && doc
                .summaries()
                .iter()
                .find(|event| event.current)
                .is_some_and(|event| {
                    event
                        .parents
                        .iter()
                        .map(String::as_str)
                        .collect::<std::collections::BTreeSet<_>>()
                        == expected_parents
                })
    });
    let mut staged = frame.clone();
    let target = if let Some(existing) = existing {
        EventBinding {
            document_id: existing.id.clone(),
            event_id: existing.head_event().into(),
            revision: existing.head_revision().into(),
        }
    } else {
        if docs.len() >= 1000 {
            return Err("The workspace document limit has been reached.".into());
        }
        let prefix = format!("assembly-junction-intent-{}", &fresh.context_revision[..16]);
        let id = (1..=1001)
            .map(|n| format!("{prefix}-{n}"))
            .find(|id| !docs.iter().any(|d| d.id == *id))
            .ok_or("No assembly document identity is available.")?;
        let derived=DocumentHistory::derive_from_events(&parents,id,fresh.target,"Declared assembly design intent from pinned source regions and explicit junction edits; synthetic junction bases are design intent, not observed DNA or verified physical preparation").map_err(|e|e.to_string())?;
        let target = EventBinding {
            document_id: derived.id.clone(),
            event_id: derived.head_event().into(),
            revision: derived.head_revision().into(),
        };
        derived.contribute(&mut staged).map_err(|e| e.to_string())?;
        target
    };
    staged.active_document = target.document_id.clone();
    if !staged
        .candidates
        .iter()
        .any(|c| c.plan_revision == candidate.plan_revision)
    {
        staged.candidates.push(candidate.clone());
    }
    staged.active_candidate = Some(candidate.plan_revision.clone());
    staged.validate().map_err(|e| e.to_string())?;
    Ok(Staged {
        workspace: staged,
        target,
        sources: sources
            .iter()
            .zip(&candidate.plan.inputs)
            .map(|(event, input)| InputEventBinding {
                input_id: input.id.clone(),
                event: event.clone(),
            })
            .collect(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use collagene_core::{
        Molecule, Topology,
        overlap_primer_design::Fragment,
        verification::{DesignPlan, Operation, PlanInput, Source},
    };
    fn fixture() -> (Workspace, Report, Vec<EventBinding>, Candidate) {
        let mut frame = Workspace {
            name: "assembly workspace".into(),
            active_document: String::new(),
            documents: vec![],
            revisions: Default::default(),
            history: vec![],
            candidates: vec![],
            active_candidate: None,
            attachments: Default::default(),
        };
        let doc = DocumentHistory::new(
            "source",
            Molecule::dna("template", "ACGTACGT", Topology::Linear).unwrap(),
        )
        .unwrap();
        let binding = EventBinding {
            document_id: doc.id.clone(),
            event_id: doc.head_event().into(),
            revision: doc.head_revision().into(),
        };
        let inputs = vec![PlanInput {
            id: "source-1".into(),
            revision: binding.revision.clone(),
            molecule: doc.molecule.clone(),
        }];
        doc.contribute(&mut frame).unwrap();
        frame.active_document = doc.id.clone();
        let intent = assembly_intent::build(
            &inputs,
            &assembly_intent::Request {
                name: "intended prefix".into(),
                fragments: vec![Fragment {
                    role: "backbone".into(),
                    input_id: "source-1".into(),
                    start: 0,
                    end: 4,
                    reverse: false,
                }],
                circular: false,
            },
        )
        .unwrap();
        let candidate = Candidate::new(DesignPlan {
            name: "proposal".into(),
            target_revision: intent.target.revision(),
            inputs,
            steps: vec![Operation::Overlap {
                fragments: vec![Source::Input("source-1".into())],
                overlaps: vec![],
                circular: false,
            }],
        });
        (frame, intent, vec![binding], candidate)
    }
    fn existing_fixture() -> (Workspace, EventBinding, Vec<InputEventBinding>, Candidate) {
        let (frame, intent, bindings, candidate) = fixture();
        let staged = stage(&frame, &intent, &bindings, &candidate).unwrap();
        let mut frame = staged.workspace;
        frame.candidates.clear();
        frame.active_candidate = None;
        (frame, staged.target, staged.sources, candidate)
    }
    #[test]
    fn existing_target_outputs_stage_only_the_candidate_without_inventing_intent() {
        let (frame, target, sources, candidate) = existing_fixture();
        let before = serde_json::to_value(&frame).unwrap();
        let staged = stage_existing(&frame, &candidate, &target, &sources).unwrap();
        assert_eq!(staged.target, target);
        assert_eq!(
            serde_json::to_value(&staged.sources).unwrap(),
            serde_json::to_value(&sources).unwrap()
        );
        assert_eq!(staged.workspace.documents.len(), 2);
        assert_eq!(staged.workspace.active_document, frame.active_document);
        assert_eq!(
            staged.workspace.active_candidate.as_ref(),
            Some(&candidate.plan_revision)
        );
        assert_eq!(staged.workspace.candidates.len(), 1);
        assert!(staged.workspace.attachments.is_empty());
        let again = stage_existing(&staged.workspace, &candidate, &target, &sources).unwrap();
        assert_eq!(
            serde_json::to_value(&again.workspace).unwrap(),
            serde_json::to_value(&staged.workspace).unwrap()
        );
        let mut without_candidate = staged.workspace;
        without_candidate.candidates.clear();
        without_candidate.active_candidate = None;
        assert_eq!(serde_json::to_value(&without_candidate).unwrap(), before);
        assert_eq!(serde_json::to_value(&frame).unwrap(), before);
    }
    #[test]
    fn existing_target_outputs_use_original_source_events_without_rewinding_heads() {
        let (mut frame, target, sources, candidate) = existing_fixture();
        let mut docs = DocumentHistory::restore_all(&frame).unwrap();
        let source = docs
            .iter_mut()
            .find(|doc| doc.id == sources[0].event.document_id)
            .unwrap();
        source
            .commit(
                source.molecule.replace(0, 1, "T").unwrap(),
                "Changed original after inverse-PCR design",
            )
            .unwrap();
        let advanced = source.head_event().to_owned();
        frame.documents.clear();
        frame.history.clear();
        frame.revisions.clear();
        for doc in &docs {
            doc.contribute(&mut frame).unwrap();
        }
        let staged = stage_existing(&frame, &candidate, &target, &sources).unwrap();
        assert_eq!(staged.sources[0].event, sources[0].event);
        assert_eq!(
            staged
                .workspace
                .documents
                .iter()
                .find(|doc| doc.id == sources[0].event.document_id)
                .unwrap()
                .head_event
                .as_deref(),
            Some(advanced.as_str())
        );
        assert_eq!(
            staged.workspace.revisions[&sources[0].event.revision]
                .sequence
                .as_str(),
            "ACGTACGT"
        );
    }
    #[test]
    fn existing_target_outputs_reject_alias_substitution_and_same_dna_new_target_events() {
        let (mut frame, target, sources, candidate) = existing_fixture();
        let alias = DocumentHistory::new(
            "equal-source-alias",
            candidate.plan.inputs[0].molecule.clone(),
        )
        .unwrap();
        alias.contribute(&mut frame).unwrap();
        let before = serde_json::to_value(&frame).unwrap();
        let mut swapped = sources.clone();
        // Equal content in this other document cannot supply the original event.
        swapped[0].event.document_id = alias.id;
        assert!(stage_existing(&frame, &candidate, &target, &swapped).is_err());
        assert_eq!(serde_json::to_value(&frame).unwrap(), before);
        let mut docs = DocumentHistory::restore_all(&frame).unwrap();
        let intended = docs
            .iter_mut()
            .find(|doc| doc.id == target.document_id)
            .unwrap();
        let original = intended.molecule.clone();
        intended
            .commit(
                original.replace(0, 1, "T").unwrap(),
                "Temporary target edit",
            )
            .unwrap();
        intended
            .commit(original, "Reverted DNA in a distinct target event")
            .unwrap();
        assert_eq!(intended.head_revision(), target.revision);
        let next = EventBinding {
            document_id: intended.id.clone(),
            event_id: intended.head_event().into(),
            revision: intended.head_revision().into(),
        };
        assert_ne!(next.event_id, target.event_id);
        frame.documents.clear();
        frame.history.clear();
        frame.revisions.clear();
        for doc in &docs {
            doc.contribute(&mut frame).unwrap();
        }
        let before = serde_json::to_value(&frame).unwrap();
        assert!(stage_existing(&frame, &candidate, &target, &sources).is_err());
        assert_eq!(serde_json::to_value(&frame).unwrap(), before);
        assert!(stage_existing(&frame, &candidate, &next, &sources).is_ok());
    }
    #[test]
    fn existing_target_outputs_require_ordered_unique_input_ids_and_the_active_target() {
        let (mut frame, target, mut sources, candidate) = existing_fixture();
        let mut plan = candidate.plan;
        let mut second = plan.inputs[0].clone();
        second.id = "second-source-role".into();
        plan.inputs.push(second);
        sources.push(InputEventBinding {
            input_id: "second-source-role".into(),
            event: sources[0].event.clone(),
        });
        let candidate = Candidate::new(plan);
        assert!(stage_existing(&frame, &candidate, &target, &sources).is_ok());
        let mut incorrect = sources.clone();
        incorrect.swap(0, 1);
        assert!(stage_existing(&frame, &candidate, &target, &incorrect).is_err());
        incorrect = sources.clone();
        incorrect[1].input_id = incorrect[0].input_id.clone();
        assert!(stage_existing(&frame, &candidate, &target, &incorrect).is_err());
        assert!(stage_existing(&frame, &candidate, &target, &sources[..1]).is_err());
        frame.active_document = sources[0].event.document_id.clone();
        assert!(stage_existing(&frame, &candidate, &target, &sources).is_err());
    }
    #[test]
    fn target_and_candidate_install_atomically_without_a_positive_check() {
        let (frame, intent, bindings, candidate) = fixture();
        let staged = stage(&frame, &intent, &bindings, &candidate).unwrap();
        assert_eq!(frame.documents.len(), 1);
        assert!(frame.candidates.is_empty());
        assert_eq!(staged.workspace.documents.len(), 2);
        assert_eq!(
            staged.workspace.revisions[&staged.target.revision]
                .sequence
                .as_str(),
            "ACGT"
        );
        assert_eq!(staged.sources[0].event, bindings[0]);
        assert!(staged.workspace.attachments.is_empty());
        assert_eq!(
            staged
                .workspace
                .history
                .iter()
                .find(|e| e.id == staged.target.event_id)
                .unwrap()
                .parents,
            vec![bindings[0].event_id.clone()]
        );
    }
    #[test]
    fn a_pinned_historical_source_is_resolved_without_rewinding_its_head() {
        let (mut frame, intent, bindings, candidate) = fixture();
        let mut docs = DocumentHistory::restore_all(&frame).unwrap();
        let edited = docs[0].molecule.replace(0, 1, "T").unwrap();
        docs[0]
            .commit(edited, "Edited after selecting assembly region")
            .unwrap();
        let advanced = docs[0].head_revision().to_owned();
        frame.documents.clear();
        frame.revisions.clear();
        frame.history.clear();
        docs[0].contribute(&mut frame).unwrap();
        let inputs = vec![InputEventBinding {
            input_id: "source-1".into(),
            event: bindings[0].clone(),
        }];
        assert_eq!(
            retained_inputs(&docs, &inputs, &candidate).unwrap()[0]
                .sequence
                .as_str(),
            "ACGTACGT"
        );
        let staged = stage(&frame, &intent, &bindings, &candidate).unwrap();
        assert_eq!(staged.workspace.documents[0].head, advanced);
        assert_eq!(staged.sources[0].event.event_id, bindings[0].event_id);
    }
    #[test]
    fn repeated_review_reuses_target_only_for_the_same_source_events() {
        let (frame, intent, bindings, candidate) = fixture();
        let first = stage(&frame, &intent, &bindings, &candidate).unwrap();
        let second = stage(&first.workspace, &intent, &bindings, &candidate).unwrap();
        assert_eq!(first.target, second.target);
        assert_eq!(second.workspace.documents.len(), 2);
        let mut aliased = second.workspace;
        let original = DocumentHistory::restore_all(&aliased).unwrap()[0]
            .molecule
            .clone();
        let alias = DocumentHistory::new("distinct-source", original).unwrap();
        let source = EventBinding {
            document_id: alias.id.clone(),
            event_id: alias.head_event().into(),
            revision: alias.head_revision().into(),
        };
        alias.contribute(&mut aliased).unwrap();
        let distinct = stage(&aliased, &intent, &[source.clone()], &candidate).unwrap();
        assert_ne!(first.target, distinct.target);
        assert_eq!(
            distinct
                .workspace
                .history
                .iter()
                .find(|event| event.id == distinct.target.event_id)
                .unwrap()
                .parents,
            vec![source.event_id]
        );
    }
    #[test]
    fn stale_binding_or_forged_intended_dna_never_changes_workspace() {
        let (frame, mut intent, mut bindings, candidate) = fixture();
        bindings[0].revision = "0".repeat(64);
        assert!(stage(&frame, &intent, &bindings, &candidate).is_err());
        let (_, _, bindings, _) = fixture();
        intent.target.sequence = collagene_core::Sequence::dna("TTTT").unwrap();
        assert!(stage(&frame, &intent, &bindings, &candidate).is_err());
        assert_eq!(frame.documents.len(), 1);
    }
    fn edited_fixture() -> (Workspace, IntendedReference, Vec<EventBinding>, Candidate) {
        let (frame, plain, bindings, candidate) = fixture();
        let mut base = plain.request;
        base.name = "Explicit inserted and trimmed junction".into();
        base.fragments = vec![
            Fragment {
                role: "left".into(),
                input_id: "source-1".into(),
                start: 0,
                end: 4,
                reverse: false,
            },
            Fragment {
                role: "right".into(),
                input_id: "source-1".into(),
                start: 4,
                end: 8,
                reverse: false,
            },
        ];
        let edits = vec![collagene_core::assembly_junction_edit::Edit {
            junction: 0,
            remove_left: 1,
            remove_right: 1,
            inserted: "GG".into(),
        }];
        let edited = collagene_core::assembly_junction_edit::build_intent(
            &candidate.plan.inputs,
            &base,
            &edits,
        )
        .unwrap();
        assert_eq!(edited.target.sequence.as_str(), "ACGGGCGT");
        let mut plan = candidate.plan;
        plan.target_revision = edited.target.revision();
        (
            frame,
            IntendedReference::Junctions(edited),
            bindings,
            Candidate::new(plan),
        )
    }
    #[test]
    fn edited_junction_target_and_candidate_stage_with_exact_ancestry_and_retry_reuse() {
        let (frame, intent, bindings, candidate) = edited_fixture();
        let before = serde_json::to_value(&frame).unwrap();
        assert_eq!(intent.target().sequence.as_str(), "ACGGGCGT");
        assert_eq!(intent.fragments().len(), 2);
        let staged = stage_any(&frame, &intent, &bindings, &candidate).unwrap();
        assert_eq!(serde_json::to_value(&frame).unwrap(), before);
        assert_eq!(
            staged.workspace.revisions[&staged.target.revision]
                .sequence
                .as_str(),
            "ACGGGCGT"
        );
        assert_eq!(staged.workspace.documents[0].head, frame.documents[0].head);
        assert_eq!(
            staged
                .workspace
                .history
                .iter()
                .find(|e| e.id == staged.target.event_id)
                .unwrap()
                .parents,
            vec![bindings[0].event_id.clone()]
        );
        assert_eq!(
            staged.workspace.candidates[0].plan_revision,
            candidate.plan_revision
        );
        let again = stage_any(&staged.workspace, &intent, &bindings, &candidate).unwrap();
        assert_eq!(again.target, staged.target);
        assert_eq!(again.workspace.documents.len(), 2);
        // The original region path still dispatches to its unchanged staging contract.
        let (old, plain, sources, c) = fixture();
        let through = stage_any(
            &old,
            &IntendedReference::Regions(plain.clone()),
            &sources,
            &c,
        )
        .unwrap();
        let direct = stage(&old, &plain, &sources, &c).unwrap();
        assert_eq!(through.target, direct.target);
    }
    #[test]
    fn edited_equal_dna_aliases_do_not_reuse_a_target_from_other_source_events() {
        let (frame, intent, bindings, candidate) = edited_fixture();
        let first = stage_any(&frame, &intent, &bindings, &candidate).unwrap();
        let mut next = first.workspace;
        let original = next.revisions[&bindings[0].revision].clone();
        let alias = DocumentHistory::new("edited-intent-distinct-source", original).unwrap();
        let binding = EventBinding {
            document_id: alias.id.clone(),
            event_id: alias.head_event().into(),
            revision: alias.head_revision().into(),
        };
        alias.contribute(&mut next).unwrap();
        let other = stage_any(&next, &intent, &[binding.clone()], &candidate).unwrap();
        assert_ne!(first.target, other.target);
        assert_eq!(
            other
                .workspace
                .history
                .iter()
                .find(|e| e.id == other.target.event_id)
                .unwrap()
                .parents,
            vec![binding.event_id]
        );
    }
    #[test]
    fn edited_staging_rejects_forged_target_or_context_and_stale_source_events() {
        let (frame, intent, bindings, candidate) = edited_fixture();
        let before = serde_json::to_value(&frame).unwrap();
        let mut stale = bindings.clone();
        stale[0].event_id = "0".repeat(64);
        assert!(stage_any(&frame, &intent, &stale, &candidate).is_err());
        let IntendedReference::Junctions(original) = intent else {
            unreachable!()
        };
        let mut forged = original.clone();
        forged.target.sequence = collagene_core::Sequence::dna("TTTTTTTT").unwrap();
        assert!(
            stage_any(
                &frame,
                &IntendedReference::Junctions(forged),
                &bindings,
                &candidate
            )
            .is_err()
        );
        let mut forged = original;
        forged.context_revision = "0".repeat(64);
        assert!(
            stage_any(
                &frame,
                &IntendedReference::Junctions(forged),
                &bindings,
                &candidate
            )
            .is_err()
        );
        assert_eq!(serde_json::to_value(&frame).unwrap(), before);
    }
}
