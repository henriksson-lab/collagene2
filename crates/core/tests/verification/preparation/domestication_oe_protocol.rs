#[path = "../../common/domestication_oe.rs"]
mod fixture;
use collagene_core::{
    cloning::golden_gate::domestication::{self, preparation as physical},
    planning::plan::Operation,
    protocol::{
        self,
        request::Request,
        schema::{AnalysisJob, AnalysisResult},
    },
    verification::{self, evidence::Status},
};
use collagene_model::Topology;
#[test]
fn whole_linear_and_circular_oe_materialize_the_frozen_literal_target_through_current_binary() {
    for topology in [Topology::Linear, Topology::Circular] {
        let (input, request) = fixture::fixture(topology);
        let source_before = input.molecule.clone();
        let intent =
            domestication::check(&input.molecule, &request.design, &request.selection).unwrap();
        assert_eq!(intent.target().sequence.as_str(), fixture::INTENDED);
        let target = intent.target().clone();
        let context = physical::revision(&input, &request, None).unwrap();
        let q = Request {
            id: 12,
            revision: context.clone(),
            job: AnalysisJob::PrepareDomestication {
                workspace: None,
                input: input.clone(),
                request: request.clone(),
            },
        };
        let encoded = protocol::request::Request::encode(&q).unwrap();
        assert_eq!(
            u16::from_le_bytes([encoded[4], encoded[5]]),
            protocol::schema::WORKER_SCHEMA_VERSION
        );
        let bytes = protocol::execution::process_request(&encoded).unwrap();
        let reply = protocol::response::Response::decode(&bytes).unwrap();
        assert_eq!(reply.id, 12);
        assert_eq!(reply.revision, context);
        let AnalysisResult::DomesticationPrepared(p) = reply.result else {
            panic!("Expected physical OE preparation")
        };
        p.validate_binding(&context).unwrap();
        assert!(p.proposal_ready, "{:?}", p.proposal_findings);
        assert_eq!(p.target, target);
        assert_eq!(p.plan.inputs[0].molecule, source_before);
        assert_eq!(
            p.plan.steps.len(),
            if topology == Topology::Circular { 6 } else { 4 }
        );
        assert!(
            matches!(&p.method,physical::MethodEvidence::OverlapExtension{designs,..} if designs.len()==3)
        );
        assert!(matches!(p.plan.steps[0], Operation::Pcr { .. }));
        assert!(matches!(p.plan.steps[1], Operation::Pcr { .. }));
        assert!(matches!(
            p.plan.steps[2],
            Operation::OverlapExtension { .. }
        ));
        assert!(matches!(p.plan.steps[3], Operation::Pcr { .. }));
        assert_eq!(
            p.verification.status,
            Status::ConsistentUnderAssumptions,
            "{:?}",
            p.verification.findings
        );
        assert_eq!(
            p.verification.products.last().unwrap().sequence.as_str(),
            fixture::INTENDED
        );
        assert!(verification::engine::plan_target_matches(
            &p.plan,
            p.verification.products.last().unwrap(),
            &target
        ));
        let wrong = verification::engine::verify(&p.plan, &source_before);
        assert_eq!(wrong.status, Status::Invalid);
        let mut changed = q.clone();
        if let AnalysisJob::PrepareDomestication {
            workspace: None,
            request,
            ..
        } = &mut changed.job
        {
            if let physical::Route::OverlapExtension { options } = &mut request.route {
                options.max_candidates = 99;
            }
        }
        assert!(
            protocol::execution::process_request(
                &protocol::request::Request::encode(&changed).unwrap()
            )
            .is_err()
        );
    }
}
