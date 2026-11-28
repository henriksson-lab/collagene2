#[path = "../../common/inverse_pcr_origin.rs"]
mod fixture;
use collagene_core::{
    cloning::pcr::mutagenesis,
    planning::plan::{Operation, Source},
    protocol::{
        self,
        request::Request,
        schema::{AnalysisJob, AnalysisResult},
    },
    verification::{self, evidence::Status},
};
#[test]
fn saved_run_and_derived_inverse_pcr_current_worker_bind_full_context_and_physical_prefix() {
    for derived in [false, true] {
        let (w, origin, _, source, target) = fixture::fixture(derived);
        let before = collagene_core::documents::workspace::save_workspace(&w).unwrap();
        let options = fixture::options();
        let context =
            mutagenesis::inverse_pcr_revision(&source, &target, Some(&origin), &options).unwrap();
        let request = Request {
            id: 44,
            revision: context.clone(),
            job: AnalysisJob::InversePcr {
                molecule: source.clone(),
                target: target.clone(),
                origin: Some(origin.clone()),
                options: options.clone(),
            },
        };
        let bytes = protocol::request::Request::encode(&request).unwrap();
        assert_eq!(
            u16::from_le_bytes([bytes[4], bytes[5]]),
            protocol::schema::WORKER_SCHEMA_VERSION
        );
        let reply = protocol::response::Response::decode(
            &protocol::execution::process_request(&bytes).unwrap(),
        )
        .unwrap();
        assert_eq!(reply.revision, context);
        let AnalysisResult::InversePcr(report) = reply.result else {
            panic!("Expected typed inversePCR report")
        };
        report.validate_binding(&context).unwrap();
        assert_eq!(report.plan.steps.len(), 4);
        assert!(matches!(
            report.plan.steps[0],
            Operation::ReconstructRunProduct { .. } | Operation::ReconstructMaterial { .. }
        ));
        assert!(matches!(
            report.plan.steps[1],
            Operation::Pcr {
                template: Source::Step(0),
                ..
            }
        ));
        assert_eq!(
            report.verification.status,
            Status::ConsistentUnderAssumptions,
            "{:?}",
            report.verification.findings
        );
        let product = report.verification.products.last().unwrap();
        // Explicit PCR-origin order: inserted GGA + retained source[5..32] + source[0..3].
        // The independent target keeps its separately declared coordinate origin.
        assert_eq!(
            product.sequence.as_str(),
            "GGAGATGCTACGTAGGACTTCAGACCTGTAAGT"
        );
        assert_eq!(target.sequence.as_str(), fixture::INTENDED);
        assert!(verification::engine::plan_target_matches(
            &report.plan,
            product,
            &target
        ));
        let mut changed = request.clone();
        if let AnalysisJob::InversePcr { options, .. } = &mut changed.job {
            options.inserted = "GGC".into();
        }
        assert!(
            protocol::execution::process_request(
                &protocol::request::Request::encode(&changed).unwrap()
            )
            .is_err()
        );
        let mut wrong = report.clone();
        wrong.context_revision = "0".repeat(64);
        assert!(wrong.validate_binding(&context).is_err());
        assert_eq!(
            collagene_core::documents::workspace::save_workspace(&w).unwrap(),
            before
        );
        let bad = verification::engine::verify(&report.plan, &source);
        assert_eq!(bad.status, Status::Invalid);
    }
}
#[test]
fn original_inverse_pcr_uses_explicit_null_origin_and_full_context_not_source_revision() {
    use collagene_model::{Molecule, Topology};
    let source = Molecule::dna("Original circle", fixture::STOCK, Topology::Circular).unwrap();
    let target =
        Molecule::dna("Independent intent", fixture::INTENDED, Topology::Circular).unwrap();
    let options = fixture::options();
    let revision = mutagenesis::inverse_pcr_revision(&source, &target, None, &options).unwrap();
    let request = Request {
        id: 1,
        revision,
        job: AnalysisJob::InversePcr {
            molecule: source.clone(),
            target,
            origin: None,
            options,
        },
    };
    let bytes = protocol::request::Request::encode(&request).unwrap();
    let reply = protocol::response::Response::decode(
        &protocol::execution::process_request(&bytes).unwrap(),
    )
    .unwrap();
    let AnalysisResult::InversePcr(report) = reply.result else {
        panic!("Expected inversePCR")
    };
    assert_eq!(report.plan.steps.len(), 3);
    assert_eq!(
        report.verification.status,
        Status::ConsistentUnderAssumptions
    );
    let mut wrong = request.clone();
    wrong.revision = source.revision();
    assert!(
        protocol::execution::process_request(&protocol::request::Request::encode(&wrong).unwrap())
            .is_err()
    );
    let mut json = serde_json::to_value(&request.job).unwrap();
    json["InversePcr"].as_object_mut().unwrap().remove("origin");
    assert!(serde_json::from_value::<AnalysisJob>(json).is_err());
}
