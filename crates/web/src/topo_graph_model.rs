//! A graph save binds an immutable native request to the complete local frame.
use collagene_core::{topo_graph_history as graph, workspace::Workspace};
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CaptureKind {
    InitialVector,
    ReviewedReaction,
    Reaction { pending_steps: usize },
}
pub struct Scope {
    pub context_revision: String,
    pub fingerprint: String,
    pub epoch: u64,
}
impl Scope {
    pub fn capture(
        frame: &Workspace,
        request: &graph::CaptureRequest,
        epoch: u64,
    ) -> Result<Self, String> {
        graph::preflight(frame, request).map_err(|e| e.to_string())?;
        Ok(Self {
            context_revision: graph::revision(frame, request).map_err(|e| e.to_string())?,
            fingerprint: collagene_core::run_lineage::workspace_fingerprint(frame)
                .map_err(|e| e.to_string())?,
            epoch,
        })
    }
    pub fn permit(&self, current: &Workspace, epoch: u64) -> Result<(), String> {
        if epoch != self.epoch
            || collagene_core::run_lineage::workspace_fingerprint(current)
                .map_err(|e| e.to_string())?
                != self.fingerprint
        {
            return Err("The workspace or graph source context changed before saving".into());
        }
        Ok(())
    }
}

pub fn initial_request(
    frame: &Workspace,
    draft: &crate::topo_preparation::Draft,
    name: String,
) -> Result<graph::CaptureRequest, String> {
    use crate::topo_preparation::{self, SourceMode};
    use collagene_core::run_lineage::EventBinding;
    if draft.source_mode != SourceMode::Original {
        return Err("An initial vector declaration uses the original-strand source mode".into());
    }
    let mut declared = draft.clone();
    let mut bindings = Vec::new();
    let mut models = Vec::new();
    for (index, (role, id)) in [
        ("vector_top", &draft.vector_top),
        ("vector_bottom", &draft.vector_bottom),
    ]
    .into_iter()
    .enumerate()
    {
        let doc = frame
            .documents
            .iter()
            .find(|doc| &doc.id == id)
            .ok_or("Choose both actual vector source documents")?;
        let model = frame.head(&doc.id).map_err(|e| e.to_string())?;
        if index == 0 {
            declared.vector_top_length = model.sequence.len();
        } else {
            declared.vector_bottom_length = model.sequence.len();
        }
        models.push(model);
        bindings.push(graph::SourceBinding {
            role: role.into(),
            event: EventBinding {
                document_id: doc.id.clone(),
                event_id: doc
                    .head_event
                    .clone()
                    .ok_or("The vector source has no exact native event")?,
                revision: doc.head.clone(),
            },
            origin: None,
        });
    }
    collagene_model::bounds::model_preflight(&models).map_err(|e| e.to_string())?;
    let request = graph::CaptureRequest {
        name,
        bindings,
        request: topo_preparation::initial_vector_request(&declared)?,
    };
    graph::preflight(frame, &request).map_err(|e| e.to_string())?;
    Ok(request)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::topo_preparation::{self as form, Draft, SourceMode};
    use collagene_core::{run_lineage::EventBinding, topo};
    use collagene_model::{Molecule, Topology};
    fn event(w: &Workspace, index: usize) -> EventBinding {
        let d = &w.documents[index];
        EventBinding {
            document_id: d.id.clone(),
            event_id: d.head_event.clone().unwrap(),
            revision: d.head.clone(),
        }
    }
    fn fixture() -> (Workspace, Draft) {
        let insert = Molecule::dna("Insert template", "CGAATG", Topology::Linear).unwrap();
        let w = Workspace::from_documents(
            "Graph source workflow",
            vec![
                Molecule::dna("Vector top", "AGGGGCATCCCTT", Topology::Linear).unwrap(),
                Molecule::dna("Vector bottom", "AGGGATGCCCCTT", Topology::Linear).unwrap(),
                insert.clone(),
                Molecule::dna(
                    "Independent target",
                    "AGGGGCATCCCTTCGAATGA",
                    Topology::Circular,
                )
                .unwrap(),
                insert,
            ],
            3,
        )
        .unwrap();
        let d=Draft{source_mode:SourceMode::Original,vector_top:w.documents[0].id.clone(),vector_bottom:w.documents[1].id.clone(),insert_template:w.documents[2].id.clone(),activation_confirmed:true,preparation_assumption:"Explicit supplied strand termini; native transfer history is a model, not measured stock".into(),forward_sequence:"CGA".into(),forward_annealing:3,reverse_sequence:"CAT".into(),reverse_annealing:3,add_a_tail:true,..Default::default()};
        (w, d)
    }
    fn chain_texts(report: &topo::TopoReport) -> Vec<String> {
        report
            .physical
            .product
            .chains
            .iter()
            .map(|chain| chain.sequence.as_str().to_owned())
            .collect()
    }
    #[test]
    fn native_initial_vector_one_transfer_archive_and_continuation_do_not_add_another_pcr() {
        let (frame, mut draft) = fixture();
        let heads = frame
            .documents
            .iter()
            .map(|d| (d.id.clone(), d.head.clone(), d.head_event.clone()))
            .collect::<Vec<_>>();
        // Initial capture is target-free and has only the two actual declared sources.
        let initial = initial_request(&frame, &draft, "Initial activated vector".into()).unwrap();
        assert_eq!(initial.bindings.len(), 2);
        let captured = graph::prepare_capture(&frame, &initial).unwrap();
        assert_eq!(captured.receipt.summary.source_bases, 26);
        assert_eq!(captured.receipt.summary.remaining_activations.len(), 2);
        assert!(captured.receipt.summary.released_activations.is_empty());
        let source =
            graph::documentary_capsule(&captured.workspace, &captured.receipt.selection).unwrap();
        draft.source_mode = SourceMode::SavedVector;
        draft.graph = Some(source.selection.clone());
        draft.vector_strands = ["0".into(), "1".into()];
        let raw_original = draft.clone();
        let seed = form::CaptureRequest {
            draft: draft.clone(),
            sources: vec![("insert_template".into(), event(&captured.workspace, 2))],
            target: event(&captured.workspace, 3),
            request: None,
            template_origin: None,
        };
        let material = form::capture(&captured.workspace, &seed).unwrap();
        assert_eq!(material.inputs.len(), 1);
        let orientations = topo::inspect_orientations(
            &material.inputs,
            &material.preparation.request,
            material.preparation.vector_strands,
            material.preparation.insert_strands,
            false,
        )
        .unwrap();
        let mut first = orientations.candidate_requests[0].clone();
        assert_eq!(first.steps.len(), 2);
        first.steps.truncate(1);
        let partial = graph::prepare_capture(
            &captured.workspace,
            &graph::CaptureRequest {
                name: "After first transfer".into(),
                bindings: vec![graph::SourceBinding {
                    role: "insert_template".into(),
                    event: event(&captured.workspace, 2),
                    origin: None,
                }],
                request: first,
            },
        )
        .unwrap();
        assert_eq!(partial.receipt.summary.source_bases, 40);
        assert_eq!(partial.receipt.summary.remaining_activations.len(), 1);
        assert_eq!(partial.receipt.summary.released_activations.len(), 1);
        let mut reopened = collagene_core::workspace::load_workspace(
            &collagene_core::workspace::save_workspace(&partial.workspace).unwrap(),
        )
        .unwrap();
        // New private source heads do not replace any graph's retained native events.
        let mut docs = collagene_core::history::DocumentHistory::restore_all(&reopened).unwrap();
        for doc in docs.iter_mut().take(3) {
            let mut changed = doc.molecule.clone();
            changed.name.push_str(" private");
            doc.commit(changed, "Private source branch").unwrap();
        }
        reopened.documents.clear();
        reopened.history.clear();
        for doc in &docs {
            doc.contribute(&mut reopened).unwrap();
        }
        draft.source_mode = SourceMode::ResumeGraph;
        draft.graph = Some(partial.receipt.selection.clone());
        draft.insert_strands = ["2".into(), "3".into()];
        draft.forward_sequence = "unfinished 1e".into();
        draft.reverse_annealing = 0;
        assert!(form::external_sources(&draft).is_empty());
        let resumed = form::capture(
            &reopened,
            &form::CaptureRequest {
                draft: draft.clone(),
                sources: vec![],
                target: event(&reopened, 3),
                request: None,
                template_origin: None,
            },
        )
        .unwrap();
        assert!(resumed.inputs.is_empty());
        assert!(resumed.template.is_none());
        assert_eq!(resumed.preparation.request.recipe.materials.len(), 1);
        assert!(matches!(
            &resumed.preparation.request.recipe.materials[0],
            topo::MaterialRecipe::SavedGraph { .. }
        ));
        let continued = topo::inspect_orientations(
            &resumed.inputs,
            &resumed.preparation.request,
            resumed.preparation.vector_strands,
            resumed.preparation.insert_strands,
            false,
        )
        .unwrap();
        assert_eq!(
            continued.candidate_requests[0].steps.len(),
            1,
            "Only the pending suffix belongs to the new reaction"
        );
        let final_graph =
            topo::simulate(&resumed.inputs, &continued.candidate_requests[0]).unwrap();
        assert_eq!(final_graph.physical.source_bases, 40);
        assert_eq!(final_graph.physical.product.released_activations.len(), 2);
        assert_eq!(final_graph.physical.product.nicks.len(), 2);
        let chains = chain_texts(&final_graph);
        assert!(chains.contains(&"AGGGGCATCCCTTCGAATGA".into()));
        assert!(chains.contains(&"AGGGATGCCCCTTCATTCGA".into()));
        // The inactive Original/Fresh-PCR controls remain literal, not silently repaired.
        assert_eq!(draft.vector_top, raw_original.vector_top);
        assert_eq!(
            draft.preparation_assumption,
            raw_original.preparation_assumption
        );
        assert_eq!(draft.forward_sequence, "unfinished 1e");
        for (id, head, head_event) in heads {
            let retained = partial
                .workspace
                .documents
                .iter()
                .find(|d| d.id == id)
                .unwrap();
            assert_eq!(retained.head, head);
            assert_eq!(retained.head_event, head_event);
        }
    }
    #[test]
    fn exact_event_graph_state_and_atomic_source_scope_refuse_stale_or_equal_dna_aliases() {
        let (frame, mut draft) = fixture();
        let request = initial_request(&frame, &draft, "Scoped vector".into()).unwrap();
        let scope = Scope::capture(&frame, &request, 7).unwrap();
        scope.permit(&frame, 7).unwrap();
        assert!(scope.permit(&frame, 8).is_err());
        let captured = graph::prepare_capture(&frame, &request).unwrap();
        assert!(
            scope.permit(&captured.workspace, 7).is_err(),
            "Another native frame must not overwrite the captured working snapshot"
        );
        draft.source_mode = SourceMode::SavedVector;
        draft.graph = Some(captured.receipt.selection.clone());
        draft.vector_strands = ["0".into(), "1".into()];
        let chosen = form::CaptureRequest {
            draft: draft.clone(),
            sources: vec![("insert_template".into(), event(&captured.workspace, 2))],
            target: event(&captured.workspace, 3),
            request: None,
            template_origin: None,
        };
        form::capture(&captured.workspace, &chosen).unwrap();
        let mut same_document = draft.clone();
        same_document.insert_template = captured.workspace.documents[3].id.clone();
        let explicit = form::capture(
            &captured.workspace,
            &form::CaptureRequest {
                draft: same_document,
                sources: vec![("insert_template".into(), event(&captured.workspace, 3))],
                target: event(&captured.workspace, 3),
                request: None,
                template_origin: None,
            },
        )
        .unwrap();
        assert_eq!(
            explicit.inputs[0].molecule.revision(),
            explicit.target.revision(),
            "Explicit source and expected-result roles may share a document without target-driven graph generation"
        );

        let alias = form::CaptureRequest {
            sources: vec![("insert_template".into(), event(&captured.workspace, 4))],
            ..chosen
        };
        assert_eq!(
            captured.workspace.documents[2].head,
            captured.workspace.documents[4].head
        );
        assert!(form::capture(&captured.workspace, &alias).is_err());
        draft.graph.as_mut().unwrap().state_revision = "0".repeat(64);
        assert!(
            form::capture(
                &captured.workspace,
                &form::CaptureRequest {
                    draft,
                    sources: vec![("insert_template".into(), event(&captured.workspace, 2))],
                    target: event(&captured.workspace, 3),
                    request: None,
                    template_origin: None
                }
            )
            .is_err()
        );
    }
}
