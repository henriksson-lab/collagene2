//! Unit-assay documentary evidence has its own bounded, exact request context.
use collagene_core::{
    cloning::restriction::{
        catalog as catalog_digest,
        conditions::{self, ChosenFormulation, EvidenceStatus, ReactionConditions},
    },
    protocol::{
        self,
        request::Request,
        response::Response,
        schema::{AnalysisJob, AnalysisResult},
    },
};
use collagene_model::Error;

fn request() -> Request {
    let catalog = catalog_digest::builtin_catalog().unwrap();
    let selection = catalog_digest::select(&catalog, &["ecori".into(), "bamhi".into()]).unwrap();
    let options = ReactionConditions {
        schema_version: 1,
        formulations: vec![
            ChosenFormulation {
                enzyme_id: "ecori".into(),
                formulation_id: Some("neb-r0101".into()),
            },
            ChosenFormulation {
                enzyme_id: "bamhi".into(),
                formulation_id: None,
            },
        ],
        buffer: "rCutSmart".into(),
        temperature_milli_c: 37000,
        incubation_minutes: 60,
        minimum_activity_percent: 75,
        activity_context: Some("NEB unit-assay relative activity".into()),
        temperature_context: Some("vendor recommended incubation".into()),
    };
    let revision = conditions::revision(&catalog, &selection, &options).unwrap();
    Request {
        id: 202,
        revision,
        job: AnalysisJob::AssessRestrictionConditions { selection, options },
    }
}
fn wire(magic: &[u8; 4], version: u16, value: &impl serde::Serialize) -> Vec<u8> {
    let payload = postcard::to_allocvec(value).unwrap();
    let mut bytes = magic.to_vec();
    bytes.extend(version.to_le_bytes());
    bytes.extend((payload.len() as u32).to_le_bytes());
    bytes.extend(blake3::hash(&payload).as_bytes());
    bytes.extend(payload);
    bytes
}
fn tag(bytes: &[u8]) -> u32 {
    let ((_, _, selector), _) =
        postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[42..]).unwrap();
    selector
}

#[test]
fn reaction_worker_uses_current_selectors_literal_context_and_retains_uncertainty() {
    let req = request();
    let AnalysisJob::AssessRestrictionConditions { selection, options } = &req.job else {
        unreachable!()
    };
    // Independent literal context; no molecule identity participates in this job.
    let expected = blake3::hash(
        &serde_json::to_vec(&(
            conditions::ENGINE_VERSION,
            &selection.catalog_edition,
            &selection.catalog_content_hash,
            selection,
            options,
        ))
        .unwrap(),
    )
    .to_hex()
    .to_string();
    assert_eq!(req.revision, expected);
    assert_eq!(protocol::schema::WORKER_SCHEMA_VERSION, 133);
    let input = req.encode().unwrap();
    assert_eq!(tag(&input), 202);
    assert_eq!(u16::from_le_bytes([input[4], input[5]]), 133);
    let output = protocol::execution::process_request(&input).unwrap();
    assert_eq!(tag(&output), 192);
    let response = Response::decode(&output).unwrap();
    assert_eq!((response.id, &response.revision), (req.id, &req.revision));
    let AnalysisResult::RestrictionConditions(report) = response.result else {
        panic!("expected conditions")
    };
    let c = catalog_digest::builtin_catalog().unwrap();
    report.validate_for(&c, selection, options).unwrap();
    assert_eq!(report.documentary_status, EvidenceStatus::BelowThreshold); // EcoRI 50% <75%
    assert_eq!(
        report.formulations[1].activity.status,
        EvidenceStatus::Unknown
    );
    assert!(report.formulations[1].formulation.is_none());
    let mut hf = request();
    let AnalysisJob::AssessRestrictionConditions { selection, options } = &mut hf.job else {
        unreachable!()
    };
    options.formulations[0].formulation_id = Some("neb-r3101".into());
    hf.revision = conditions::revision(&c, selection, options).unwrap();
    let response =
        Response::decode(&protocol::execution::process_request(&hf.encode().unwrap()).unwrap())
            .unwrap();
    let AnalysisResult::RestrictionConditions(report) = response.result else {
        panic!("expected HF conditions")
    };
    assert_eq!(
        report.formulations[0].activity.status,
        EvidenceStatus::Known
    ); // exact HF 100%
    assert_eq!(report.documentary_status, EvidenceStatus::Unknown); // unchosen BamHI remains unknown
}

#[test]
fn raw_worker_rejects_stale_context_old_schema_and_invalid_native_options() {
    let version = protocol::schema::WORKER_SCHEMA_VERSION;
    let mut stale = request();
    let AnalysisJob::AssessRestrictionConditions { options, .. } = &mut stale.job else {
        unreachable!()
    };
    options.incubation_minutes += 1;
    assert!(stale.encode().is_err());
    assert!(protocol::execution::process_request(&wire(b"CGQ\0", version, &stale)).is_err());
    let req = request();
    assert!(protocol::execution::process_request(&wire(b"CGQ\0", 128, &req)).is_err());
    let valid = protocol::execution::process_request(&req.encode().unwrap()).unwrap();
    let mut obsolete = valid;
    obsolete[4..6].copy_from_slice(&128u16.to_le_bytes());
    assert!(Response::decode(&obsolete).is_err());
    let mut forged = request();
    let AnalysisJob::AssessRestrictionConditions { selection, .. } = &mut forged.job else {
        unreachable!()
    };
    selection.enzymes[0].enzyme.top_cut += 1;
    assert!(protocol::execution::process_request(&wire(b"CGQ\0", version, &forged)).is_err());
    let mut large = request();
    let AnalysisJob::AssessRestrictionConditions { options, .. } = &mut large.job else {
        unreachable!()
    };
    options.formulations = vec![options.formulations[0].clone(); 65];
    assert!(large.encode().is_err());
    assert!(protocol::execution::process_request(&wire(b"CGQ\0", version, &large)).is_err());
    let mut oversized = request();
    let AnalysisJob::AssessRestrictionConditions { options, .. } = &mut oversized.job else {
        unreachable!()
    };
    options.buffer = "X".repeat(conditions::MAX_OPTIONS_BYTES + 1);
    assert!(oversized.encode().is_err());
    assert!(protocol::execution::process_request(&wire(b"CGQ\0", version, &oversized)).is_err());
}

#[test]
fn result_decoder_refuses_forged_binding_or_mutated_documentary_evidence() {
    let req = request();
    let response =
        Response::decode(&protocol::execution::process_request(&req.encode().unwrap()).unwrap())
            .unwrap();
    let mut forged = response.clone();
    forged.revision = "0".repeat(64);
    assert!(Response::decode(&forged.encode().unwrap()).is_err());
    let mut forged = response;
    let AnalysisResult::RestrictionConditions(report) = &mut forged.result else {
        unreachable!()
    };
    report.formulations[0].activity.status = EvidenceStatus::Known;
    assert!(Response::decode(&forged.encode().unwrap()).is_err());
}

#[test]
fn both_reaction_selectors_retain_the_ordinary_32_mib_frame_limit() {
    let maximum = collagene_core::formats::sequence::MAX_IMPORT_BYTES;
    assert_eq!(maximum, 32 * 1024 * 1024);
    let mut bytes = vec![0; 42 + maximum + 1];
    bytes[4..6].copy_from_slice(&protocol::schema::WORKER_SCHEMA_VERSION.to_le_bytes());
    bytes[6..10].copy_from_slice(&((maximum + 1) as u32).to_le_bytes());
    for (magic, selector, is_response) in [(b"CGQ\0", 202u32, false), (b"CGR\0", 192u32, true)] {
        bytes[..4].copy_from_slice(magic);
        let prefix = postcard::to_allocvec(&(202u32, "context", selector)).unwrap();
        bytes[42..42 + prefix.len()].copy_from_slice(&prefix);
        let result = if is_response {
            Response::decode(&bytes).map(|_| ())
        } else {
            protocol::execution::process_request(&bytes).map(|_| ())
        };
        assert!(matches!(result, Err(Error::TooLarge)));
    }
}
