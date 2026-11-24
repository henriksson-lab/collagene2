//! Atomic installation of an intended reference and its proposed preparation.
//! Exact source events are checked; no scientific success is restored or inferred.
use collagene_core::{
    assembly_intent::{self, Report},
    history::DocumentHistory,
    run_lineage::{EventBinding, InputEventBinding},
    transfer::Candidate,
    workspace::Workspace,
};
pub struct Staged {
    pub workspace: Workspace,
    pub target: EventBinding,
    pub sources: Vec<InputEventBinding>,
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
    if docs.len() >= 1000 {
        return Err("The workspace document limit has been reached.".into());
    }
    let prefix = format!("assembly-intent-{}", &fresh.context_revision[..16]);
    let id = (1..=1001)
        .map(|n| format!("{prefix}-{n}"))
        .find(|id| !docs.iter().any(|d| d.id == *id))
        .ok_or("No assembly document identity is available.")?;
    let derived=DocumentHistory::derive_from_events(&parents,id,fresh.target,
        "Declared assembly design intent from pinned source regions; expected reference only, preparation and experimental outcome unverified").map_err(|e|e.to_string())?;
    let target = EventBinding {
        document_id: derived.id.clone(),
        event_id: derived.head_event().into(),
        revision: derived.head_revision().into(),
    };
    let mut staged = frame.clone();
    derived.contribute(&mut staged).map_err(|e| e.to_string())?;
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
    fn stale_binding_or_forged_intended_dna_never_changes_workspace() {
        let (frame, mut intent, mut bindings, candidate) = fixture();
        bindings[0].revision = "0".repeat(64);
        assert!(stage(&frame, &intent, &bindings, &candidate).is_err());
        let (_, _, bindings, _) = fixture();
        intent.target.sequence = collagene_core::Sequence::dna("TTTT").unwrap();
        assert!(stage(&frame, &intent, &bindings, &candidate).is_err());
        assert_eq!(frame.documents.len(), 1);
    }
}
