//! Current binary graph operations preserve native ownership and full graph scope.
#[path = "../../../common/topo_graph.rs"]
mod literal;
use collagene_core::{
    protocol::{
        self,
        request::Request,
        schema::{AnalysisJob as J, AnalysisResult as R},
    },
    topo as t, topo_graph_history as g, workspace,
};
fn run(revision: String, job: J) -> R {
    let bytes = protocol::request::Request::encode(&Request {
        id: 19,
        revision: revision.clone(),
        job,
    })
    .unwrap();
    assert_eq!(
        u16::from_le_bytes(bytes[4..6].try_into().unwrap()),
        protocol::schema::WORKER_SCHEMA_VERSION
    );
    let reply = protocol::response::Response::decode(
        &protocol::execution::process_request(&bytes).unwrap(),
    )
    .unwrap();
    assert_eq!(reply.revision, revision);
    reply.result
}
#[test]
fn native_capture_list_inspect_and_partial_resume_keep_all_calls_and_nicks() {
    let (w, r) = literal::fixture();
    let expected = g::revision(&w, &r).unwrap();
    let R::TopoGraphCaptured(first) = run(
        expected.clone(),
        J::CaptureTopoGraph {
            workspace: w.clone(),
            request: r,
        },
    ) else {
        panic!("capture")
    };
    first.validate_binding(&expected).unwrap();
    assert_eq!(first.receipt.summary.source_bases, 26);
    assert_eq!(first.receipt.summary.remaining_activations.len(), 2);
    assert_eq!(first.workspace.documents.len(), w.documents.len());
    assert_eq!(first.workspace.history.len(), w.history.len());
    let R::TopoGraphsListed(list) = run(
        g::list_revision(&first.workspace).unwrap(),
        J::ListTopoGraphs {
            workspace: first.workspace.clone(),
        },
    ) else {
        panic!("list")
    };
    assert_eq!(list.entries.len(), 1);
    let selection = first.receipt.selection.clone();
    let expected = g::inspect_revision(&first.workspace, &selection).unwrap();
    let R::TopoGraphInspected(inspected) = run(
        expected.clone(),
        J::InspectTopoGraph {
            workspace: first.workspace.clone(),
            selection: selection.clone(),
        },
    ) else {
        panic!("inspect")
    };
    inspected.validate_binding(&expected).unwrap();
    assert_eq!(inspected.summary.source_bases, 26);
    let second = g::prepare_capture(
        &first.workspace,
        &literal::partial(&first.workspace, &selection),
    )
    .unwrap();
    assert_eq!(second.receipt.summary.source_bases, 40);
    assert_eq!(second.receipt.summary.remaining_activations.len(), 1);
    assert_eq!(second.receipt.summary.released_activations.len(), 1);
    let reopened =
        workspace::load_workspace(&workspace::save_workspace(&second.workspace).unwrap()).unwrap();
    let request = literal::continuation(&reopened, &second.receipt.selection);
    let batch = t::batch_revision(&[], &request).unwrap();
    let selection = t::Selection {
        batch_revision: batch,
        component: 0,
        reference_seed: t::physical::Call {
            strand: 0,
            index: 0,
        },
    };
    let product = t::product(&[], &request, &selection).unwrap();
    assert_eq!(
        product.reference.as_ref().unwrap().sequence.as_str(),
        literal::FINAL
    );
    assert_eq!(product.projection.physical.source_bases, 40);
    assert_eq!(product.projection.physical.product.chains.len(), 2);
    assert_eq!(product.projection.physical.product.nicks.len(), 2);
    assert_eq!(
        product
            .projection
            .physical
            .product
            .released_activations
            .len(),
        2
    );
    assert_ne!(
        product.projection.physical.product.released_activations[0],
        product.projection.physical.product.released_activations[1]
    );
    let target = &reopened.revisions[&literal::event(&reopened, 4).revision];
    let mut policy = t::AcceptancePolicy::default();
    policy.allow_nicks = true;
    let assessment = t::check_target(&[], &request, &selection, target, &policy).unwrap();
    assert!(matches!(
        assessment.outcome,
        t::TargetOutcome::AcceptedUnderDeclaredGraph
    ));
    assert!(!matches!(
        t::check_target(
            &[],
            &request,
            &selection,
            target,
            &t::AcceptancePolicy::default()
        )
        .unwrap()
        .outcome,
        t::TargetOutcome::AcceptedUnderDeclaredGraph
    ));
}
#[test]
fn selector_state_record_and_workspace_ownership_cannot_be_substituted() {
    let first = literal::captured();
    let mut selected = first.receipt.selection.clone();
    selected.state_revision = "0".repeat(64);
    assert!(g::inspect(&first.workspace, &selected).is_err());
    let mut foreign = first.workspace.clone();
    foreign.attachments.remove(&first.receipt.record_hash);
    assert!(g::documentary_capsule(&foreign, &first.receipt.selection).is_err());
    let mut capsule = g::documentary_capsule(&first.workspace, &first.receipt.selection).unwrap();
    capsule.record_bytes[0] ^= 1;
    assert!(g::replay(&capsule).is_err());
    let mut sources = g::documentary_capsule(&first.workspace, &first.receipt.selection).unwrap();
    sources.source_snapshots[0].molecule.name = "Equal DNA, altered native source".into();
    assert!(g::replay(&sources).is_err());
}
