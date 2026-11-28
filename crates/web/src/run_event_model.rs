//! Read-only exact-version navigation from a retained cloning run. Opening a
//! version never restores a scientific assessment or changes a document head.
use collagene_core::{
    history::DocumentHistory,
    run_lineage::{self, EventBinding},
    workspace::{self, Workspace},
};
use collagene_model::Molecule;

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Role {
    Target,
    Input { input_id: String },
    Product { step_index: usize },
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Request {
    pub run_id: String,
    pub binding: EventBinding,
    /// The host compares this with its live workspace epoch before acting.
    pub epoch: u64,
    pub role: Role,
}

#[derive(Clone, Debug)]
pub struct Snapshot {
    pub request: Request,
    pub molecule: Molecule,
    pub current_head: bool,
}

fn digest(value: &str) -> bool {
    value.len() == 64
        && value
            .bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}

fn identifier(value: &str) -> bool {
    !value.is_empty()
        && value.len() <= 128
        && value
            .bytes()
            .all(|b| b.is_ascii_alphanumeric() || matches!(b, b'_' | b'-' | b'.' | b':'))
}

/// Bound borrowed models before workspace validation hashes or serialized run
/// loading can traverse them. Final workspace validation retains the existing
/// aggregate archive/manifest/attachment budgets.
fn preflight(frame: &Workspace, request: &Request) -> Result<(), String> {
    if !digest(&request.run_id)
        || !identifier(&request.binding.document_id)
        || !digest(&request.binding.revision)
        || !request
            .binding
            .event_id
            .strip_prefix("event-")
            .is_some_and(digest)
        || matches!(&request.role, Role::Input { input_id } if input_id.trim().is_empty() || input_id.len()>128)
        || matches!(request.role, Role::Product { step_index } if step_index >= 1000)
    {
        return Err("The saved-run version request has invalid identities or bounds.".into());
    }
    if frame.documents.len() > workspace::MAX_DOCUMENTS
        || frame.revisions.len() > workspace::MAX_REVISIONS
        || frame.history.len() > workspace::MAX_HISTORY_EVENTS
        || frame.candidates.len() > workspace::MAX_CANDIDATES
        || frame.attachments.len() > workspace::MAX_ATTACHMENTS
        || frame
            .attachments
            .values()
            .any(|asset| asset.bytes.len() > workspace::MAX_ENTRY_BYTES)
    {
        return Err("The saved-run workspace exceeds its native bounds.".into());
    }
    // Raw letters and attachment bytes are a lower bound on the expanded
    // archive size; reject an impossible frame before hashing these contents.
    let mut minimum_bytes = 0usize;
    for bytes in frame
        .revisions
        .values()
        .map(|m| m.sequence.len())
        .chain(frame.attachments.values().map(|asset| asset.bytes.len()))
    {
        minimum_bytes = minimum_bytes
            .checked_add(bytes)
            .filter(|bytes| *bytes <= workspace::MAX_EXPANDED_BYTES)
            .ok_or("The retained workspace exceeds its expanded byte budget.")?;
    }
    for molecule in frame.revisions.values() {
        collagene_model::bounds::model_preflight(&[molecule]).map_err(|e| e.to_string())?;
    }
    for candidate in &frame.candidates {
        if candidate.plan.inputs.len() > 1000 || candidate.plan.steps.len() > 1000 {
            return Err("The retained candidate exceeds its native bounds.".into());
        }
        for input in &candidate.plan.inputs {
            collagene_model::bounds::model_preflight(&[&input.molecule])
                .map_err(|e| e.to_string())?;
        }
    }
    Ok(())
}

fn resolve_molecule<'a>(
    frame: &'a Workspace,
    request: &Request,
) -> Result<(&'a Molecule, bool), String> {
    preflight(frame, request)?;
    // load_run validates archive consistency and exact run cross-references. It
    // deliberately does not replay the scientific calculation.
    let record = run_lineage::load_run(frame, &request.run_id).map_err(|e| e.to_string())?;
    let member = match &request.role {
        Role::Target => record.target == request.binding,
        Role::Input { input_id } => record
            .inputs
            .iter()
            .any(|input| input.input_id == *input_id && input.event == request.binding),
        Role::Product { step_index } => record.products.iter().any(|product| {
            product.step_index == *step_index
                && product.document_id == request.binding.document_id
                && product.event_id == request.binding.event_id
                && product.revision == request.binding.revision
        }),
    };
    if !member {
        return Err(
            "The requested version is not the recorded source, target or product for this role."
                .into(),
        );
    }
    let document = frame
        .documents
        .iter()
        .find(|document| document.id == request.binding.document_id)
        .ok_or("The exact recorded document is no longer retained.")?;
    let event = frame
        .history
        .iter()
        .find(|event| event.id == request.binding.event_id)
        .ok_or("The exact recorded version is no longer retained.")?;
    if event.document_id != request.binding.document_id
        || event.revision != request.binding.revision
    {
        return Err("The recorded version differs from its original document binding.".into());
    }
    let molecule = frame
        .revisions
        .get(&event.revision)
        .ok_or("The recorded version's sequence is no longer retained.")?;
    let current_head = document.head_event.as_deref() == Some(request.binding.event_id.as_str())
        && document.head == request.binding.revision;
    Ok((molecule, current_head))
}

/// Clone only the selected, bounded immutable version for the read-only viewer.
/// Equal DNA in another event never substitutes for missing run membership.
pub fn resolve(frame: &Workspace, request: &Request) -> Result<Snapshot, String> {
    let (molecule, current_head) = resolve_molecule(frame, request)?;
    Ok(Snapshot {
        request: request.clone(),
        molecule: molecule.clone(),
        current_head,
    })
}

/// Explicitly fork the selected recorded version, retaining its exact event as
/// the new document's parent. No source head or saved assessment is changed.
pub fn fork(frame: &Workspace, request: &Request, new_id: &str) -> Result<Workspace, String> {
    resolve_molecule(frame, request)?;
    if !identifier(new_id)
        || frame.documents.iter().any(|document| document.id == new_id)
        || frame.documents.len() >= workspace::MAX_DOCUMENTS
        || frame.history.len() >= workspace::MAX_HISTORY_EVENTS
    {
        return Err("Choose an unused document identity within the workspace limits.".into());
    }
    let histories = DocumentHistory::restore_all(frame).map_err(|e| e.to_string())?;
    let source = histories
        .iter()
        .find(|document| document.id == request.binding.document_id)
        .ok_or("The recorded source document is unavailable.")?;
    let fork = source
        .fork_event(&request.binding.event_id, new_id)
        .map_err(|e| e.to_string())?;
    let mut next = frame.clone();
    fork.contribute(&mut next).map_err(|e| e.to_string())?;
    next.active_document = new_id.into();
    next.validate().map_err(|e| e.to_string())?;
    Ok(next)
}

#[cfg(test)]
mod tests {
    use super::*;
    use collagene_core::{
        analysis::Primer,
        run_lineage::InputEventBinding,
        verification::{DesignPlan, Operation, PlanInput, Source},
    };
    use collagene_model::Topology;

    const ORIGINAL: &str = "ACGACCCCGGGGTAAC";

    fn binding(frame: &Workspace, index: usize) -> EventBinding {
        let document = &frame.documents[index];
        EventBinding {
            document_id: document.id.clone(),
            event_id: document.head_event.clone().unwrap(),
            revision: document.head.clone(),
        }
    }

    fn fixture() -> (Workspace, Request) {
        let source = Molecule::dna("original source", ORIGINAL, Topology::Linear).unwrap();
        let target = Molecule::dna("independent target", ORIGINAL, Topology::Linear).unwrap();
        let frame = Workspace::from_documents(
            "saved run versions",
            vec![source.clone(), source.clone(), target.clone()],
            2,
        )
        .unwrap();
        let plan = DesignPlan {
            name: "Literal reference PCR".into(),
            target_revision: target.revision(),
            inputs: vec![PlanInput {
                id: "original rôle 1".into(),
                revision: source.revision(),
                molecule: source,
            }],
            steps: vec![Operation::Pcr {
                template_strand: collagene_core::chemistry::Strand::Top,
                template: Source::Input("original rôle 1".into()),
                forward: Primer::new("Forward", "ACGACCCC", "").unwrap(),
                reverse: Primer::new("Reverse", "GTTACCCC", "").unwrap(),
            }],
        };
        let original = binding(&frame, 0);
        let commit = run_lineage::record_run(
            &frame,
            &plan,
            binding(&frame, 2),
            vec![InputEventBinding {
                input_id: "original rôle 1".into(),
                event: original.clone(),
            }],
            false,
        )
        .unwrap();
        (
            commit.workspace,
            Request {
                run_id: commit.receipt.run_id,
                binding: original,
                epoch: 7,
                role: Role::Input {
                    input_id: "original rôle 1".into(),
                },
            },
        )
    }

    fn advance_and_revert_source(frame: &Workspace) -> Workspace {
        let mut histories = DocumentHistory::restore_all(frame).unwrap();
        let source = &mut histories[0];
        let original = source.molecule.clone();
        source
            .commit(
                Molecule::dna("later source", "TCGACCCCGGGGTAAC", Topology::Linear).unwrap(),
                "Later source edit",
            )
            .unwrap();
        source.commit(original, "Restore original letters").unwrap();
        let mut next = frame.clone();
        next.documents.clear();
        next.history.clear();
        for history in histories {
            history.contribute(&mut next).unwrap();
        }
        next.validate().unwrap();
        next
    }

    #[test]
    fn saved_roles_resolve_exact_versions_without_rewinding_a_later_head() {
        let (frame, source) = fixture();
        let current = resolve(&frame, &source).unwrap();
        assert!(current.current_head);
        assert_eq!(current.molecule.sequence.as_str(), ORIGINAL);
        let later = advance_and_revert_source(&frame);
        assert_eq!(later.documents[0].head, source.binding.revision);
        assert_ne!(
            later.documents[0].head_event.as_ref().unwrap(),
            &source.binding.event_id
        );
        let before = serde_json::to_value(&later).unwrap();
        let historical = resolve(&later, &source).unwrap();
        assert!(!historical.current_head);
        assert_eq!(historical.request, source);
        assert_eq!(historical.molecule.sequence.as_str(), ORIGINAL);
        assert_eq!(serde_json::to_value(&later).unwrap(), before);

        let record = run_lineage::load_run(&later, &source.run_id).unwrap();
        let target = Request {
            binding: record.target,
            role: Role::Target,
            ..source.clone()
        };
        assert!(resolve(&later, &target).unwrap().current_head);
        let product = &record.products[0];
        let product = Request {
            binding: EventBinding {
                document_id: product.document_id.clone(),
                event_id: product.event_id.clone(),
                revision: product.revision.clone(),
            },
            role: Role::Product { step_index: 0 },
            ..source
        };
        let snapshot = resolve(&later, &product).unwrap();
        assert!(snapshot.current_head);
        assert_eq!(snapshot.molecule.sequence.as_str(), ORIGINAL);
    }

    #[test]
    fn same_dna_aliases_new_events_wrong_roles_and_missing_records_are_refused() {
        let (frame, request) = fixture();
        let mut spoof = request.clone();
        spoof.binding = binding(&frame, 1);
        assert_eq!(spoof.binding.revision, request.binding.revision);
        assert!(resolve(&frame, &spoof).is_err());
        let later = advance_and_revert_source(&frame);
        spoof.binding = binding(&later, 0);
        assert_eq!(spoof.binding.revision, request.binding.revision);
        assert!(resolve(&later, &spoof).is_err());
        spoof = request.clone();
        spoof.role = Role::Target;
        assert!(resolve(&frame, &spoof).is_err());
        spoof.role = Role::Product { step_index: 0 };
        assert!(resolve(&frame, &spoof).is_err());
        spoof.role = Role::Input {
            input_id: "another-input".into(),
        };
        assert!(resolve(&frame, &spoof).is_err());
        spoof = request.clone();
        spoof.run_id = "0".repeat(64);
        assert!(resolve(&frame, &spoof).is_err());
        let mut missing = frame.clone();
        missing.attachments.clear();
        assert!(resolve(&missing, &request).is_err());
        let mut orphan = frame.clone();
        orphan
            .history
            .retain(|event| event.id != request.binding.event_id);
        assert!(resolve(&orphan, &request).is_err());
    }

    #[test]
    fn explicit_fork_retains_the_recorded_parent_and_all_original_heads_atomically() {
        let (frame, request) = fixture();
        let frame = advance_and_revert_source(&frame);
        let before = serde_json::to_value(&frame).unwrap();
        let next = fork(&frame, &request, "editable-source-copy").unwrap();
        assert_eq!(serde_json::to_value(&frame).unwrap(), before);
        assert_eq!(next.documents.len(), frame.documents.len() + 1);
        assert_eq!(
            serde_json::to_value(&next.documents[..frame.documents.len()]).unwrap(),
            serde_json::to_value(&frame.documents).unwrap()
        );
        assert_eq!(next.active_document, "editable-source-copy");
        let head = next.documents.last().unwrap().head_event.as_ref().unwrap();
        let event = next.history.iter().find(|event| event.id == *head).unwrap();
        assert_eq!(event.parents, vec![request.binding.event_id.clone()]);
        assert_eq!(event.revision, request.binding.revision);
        assert_eq!(
            next.head("editable-source-copy").unwrap().sequence.as_str(),
            ORIGINAL
        );
        assert_eq!(
            serde_json::to_value(&next.attachments).unwrap(),
            serde_json::to_value(&frame.attachments).unwrap()
        );
        assert_eq!(
            serde_json::to_value(&next.candidates).unwrap(),
            serde_json::to_value(&frame.candidates).unwrap()
        );
        assert_eq!(next.active_candidate, frame.active_candidate);
        assert!(fork(&next, &request, "editable-source-copy").is_err());
        assert!(fork(&frame, &request, &request.binding.document_id).is_err());
        assert!(fork(&frame, &request, "../invalid").is_err());
        assert_eq!(serde_json::to_value(&frame).unwrap(), before);
        assert!(!resolve(&next, &request).unwrap().current_head);
        // A fork is an editable descendant, not a replacement recorded source.
        let mut substitute = request;
        substitute.binding = EventBinding {
            document_id: "editable-source-copy".into(),
            event_id: head.clone(),
            revision: event.revision.clone(),
        };
        assert!(resolve(&next, &substitute).is_err());
    }
    #[test]
    fn unrelated_supported_candidate_does_not_block_saved_run_navigation() {
        let (mut frame, request) = fixture();
        let mut plan = run_lineage::load_run(&frame, &request.run_id).unwrap().plan;
        let input = plan.inputs[0].clone();
        plan.inputs.extend((1..129).map(|i| {
            let mut input = input.clone();
            input.id = format!("unused-{i}");
            input
        }));
        let proposal = collagene_core::transfer::Candidate::new(plan);
        proposal.validate().unwrap();
        frame.candidates.push(proposal);
        frame.validate().unwrap();
        assert_eq!(
            resolve(&frame, &request)
                .unwrap()
                .molecule
                .sequence
                .as_str(),
            ORIGINAL
        );
    }
}
