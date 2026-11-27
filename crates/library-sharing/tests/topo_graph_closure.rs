//! Only actual native graph consumers own graph dependency publication.
#[path = "../../core/tests/common/topo_graph.rs"]
mod literal;
use collagene_core::{
    cloning_outputs, topo as t, topo_graph_history as g,
    verification::{self, DesignPlan, Operation},
    workspace,
};
use collagene_library_sharing as sharing;
fn selected(e: &collagene_core::run_lineage::EventBinding) -> sharing::Selection {
    sharing::Selection {
        document_id: e.document_id.clone(),
        revision: e.revision.clone(),
        head_event: Some(e.event_id.clone()),
    }
}
#[test]
fn unrelated_plasmid_never_expands_to_saved_activated_graph_or_private_sources() {
    let first = literal::captured();
    let p = sharing::project_document_closure(
        &first.workspace,
        &selected(&literal::event(&first.workspace, 5)),
    )
    .unwrap();
    assert!(p.workspace.attachments.is_empty());
    assert_eq!(p.inventory.omitted.attachments, 1);
    assert_eq!(p.workspace.documents.len(), 1);
    assert!(
        !p.workspace
            .revisions
            .contains_key(&literal::event(&first.workspace, 0).revision)
    );
    assert!(
        first
            .workspace
            .attachments
            .contains_key(&first.receipt.record_hash)
    );
}
#[test]
fn graph_consuming_run_closes_both_records_and_exact_sources_then_replays() {
    let first = literal::captured();
    let second = g::prepare_capture(
        &first.workspace,
        &literal::partial(&first.workspace, &first.receipt.selection),
    )
    .unwrap();
    let request = literal::continuation(&second.workspace, &second.receipt.selection);
    let plan = DesignPlan {
        name: "Independent target check of actual graph continuation".into(),
        target_revision: literal::event(&second.workspace, 4).revision,
        inputs: vec![],
        steps: vec![Operation::CheckedTopo {
            sources: vec![],
            selection: t::Selection {
                batch_revision: t::batch_revision(&[], &request).unwrap(),
                component: 0,
                reference_seed: t::physical::Call {
                    strand: 0,
                    index: 0,
                },
            },
            request,
            policy: t::AcceptancePolicy {
                allow_nicks: true,
                ..Default::default()
            },
        }],
    };
    let captured = cloning_outputs::finalize(
        &second.workspace,
        &plan,
        literal::event(&second.workspace, 4),
        vec![],
        &cloning_outputs::Options {
            name: "Native graph checked reference".into(),
            show_prepared: false,
            save_primers: false,
            allow_incomplete: true,
        },
    )
    .unwrap();
    let p = sharing::project_document_closure(
        &captured.workspace,
        &selected(&captured.receipt.construct),
    )
    .unwrap();
    assert!(
        p.workspace
            .attachments
            .contains_key(&first.receipt.record_hash)
    );
    assert!(
        p.workspace
            .attachments
            .contains_key(&second.receipt.record_hash)
    );
    assert!(
        !p.workspace
            .documents
            .iter()
            .any(|d| d.id == captured.workspace.documents[5].id)
    );
    let reopened =
        workspace::load_workspace(&workspace::save_workspace(&p.workspace).unwrap()).unwrap();
    let inspected = g::inspect(&reopened, &second.receipt.selection).unwrap();
    assert_eq!(inspected.summary.source_bases, 40);
    let target = &reopened.revisions[&plan.target_revision];
    let replayed = verification::verify(&plan, target);
    assert!(verification::plan_target_matches(
        &plan,
        replayed.products.last().unwrap(),
        target
    ));
    assert_eq!(
        replayed.products.last().unwrap().sequence.as_str(),
        literal::FINAL
    );
    let mut missing = captured.workspace.clone();
    missing.attachments.remove(&first.receipt.record_hash);
    assert!(
        sharing::project_document_closure(&missing, &selected(&captured.receipt.construct))
            .is_err()
    );
}
