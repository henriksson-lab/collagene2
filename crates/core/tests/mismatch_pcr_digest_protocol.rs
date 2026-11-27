//! Current mismatch-PCR catalog inspection binds source, policy and actual cut geometry.
use collagene_core::{
    Molecule, Topology,
    analysis::Primer,
    assembly_recipe::{FragmentSpec, Preparation, RunPreparation},
    catalog_digest, ligation_recipe,
    mismatch_pcr_digest_inspection::{self, Request as DigestRequest},
    mutagenesis::PcrOptions,
    pcr_digest_inspection,
    protocol::{self, AnalysisJob, AnalysisResult, Request},
    restriction_inspection::InspectionOptions,
    verification::{Operation, PlanInput, Status},
};
fn fixture() -> (Molecule, DigestRequest) {
    let source = Molecule::dna("Original", "ACGACCCCGGGGTAAC", Topology::Linear).unwrap();
    let c = catalog_digest::builtin_catalog().unwrap();
    let id = c
        .edition()
        .enzymes
        .iter()
        .find(|e| e.name == "EcoRI")
        .unwrap()
        .id
        .clone();
    let mut policy = PcrOptions::default();
    policy.binding.max_mismatches = 1;
    policy.binding.min_three_prime_exact = 2;
    (
        source,
        DigestRequest {
            base: pcr_digest_inspection::Request {
                source_role: "original".into(),
                forward: Primer::new("Mismatch F", "ACTACCCC", "CCGAATTC").unwrap(),
                reverse: Primer::new("Exact R", "GTTACCCC", "GGGAATTC").unwrap(),
                selection: catalog_digest::select(&c, &[id]).unwrap(),
                options: InspectionOptions::default(),
                template_strand: collagene_core::chemistry::Strand::Top,
            },
            pcr_options: policy,
        },
    )
}
fn wire(magic: &[u8; 4], v: u16, value: &impl serde::Serialize) -> Vec<u8> {
    let p = postcard::to_allocvec(value).unwrap();
    let mut b = magic.to_vec();
    b.extend(v.to_le_bytes());
    b.extend((p.len() as u32).to_le_bytes());
    b.extend(blake3::hash(&p).as_bytes());
    b.extend(p);
    b
}

#[test]
fn inspect162_154_binds_real_policy_amplicon_and_catalog_without_target_or_header_downgrade() {
    let (molecule, request) = fixture();
    let revision =
        mismatch_pcr_digest_inspection::inspection_revision(&molecule, None, &request, None)
            .unwrap();
    let r = Request {
        id: 162,
        revision: revision.clone(),
        job: AnalysisJob::InspectMismatchPcrDigest {
            molecule,
            origin: None,
            request,
            workspace: None,
        },
    };
    let bytes = protocol::encode_request(&r).unwrap();
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[42..]).unwrap();
    assert_eq!(tag, 157);
    assert!(protocol::process_request(&wire(b"CGQ\0", 77, &r)).is_err());
    let bytes = protocol::process_request(&bytes).unwrap();
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[42..]).unwrap();
    assert_eq!(tag, 149);
    let response = protocol::decode_response(&bytes).unwrap();
    let AnalysisResult::MismatchPcrDigestInspection(p) = &response.result else {
        panic!("inspection")
    };
    p.validate_binding(&revision).unwrap();
    assert!(!p.incomplete);
    assert_eq!(
        p.amplicon.as_ref().unwrap().sequence.as_str(),
        "CCGAATTCACTACCCCGGGGTAACGAATTCCC"
    );
    assert_eq!(
        p.inspection
            .as_ref()
            .unwrap()
            .fragments
            .iter()
            .map(|f| f.duplex.top_length)
            .collect::<Vec<_>>(),
        vec![3, 22, 7]
    );
    assert_eq!(p.request.pcr_options.binding.min_three_prime_exact, 2);
    assert!(protocol::decode_response(&wire(b"CGR\0", 77, &response)).is_err());
    let mut stale = r;
    stale.revision = "0".repeat(64);
    assert!(protocol::process_request(&protocol::encode_request(&stale).unwrap()).is_err());
    let mut stale = response;
    stale.revision = "0".repeat(64);
    assert!(
        protocol::decode_response(&wire(b"CGR\0", protocol::WORKER_SCHEMA_VERSION, &stale))
            .is_err()
    );
}
#[test]
fn mismatch_digest_preparations_roundtrip_current_native_evidence() {
    let (source, r) = fixture();
    let b = r.base;
    let prep = Preparation::MismatchPcrCatalogDigest {
        forward: b.forward.clone(),
        reverse: b.reverse.clone(),
        options: r.pcr_options.clone(),
        selection: b.selection.clone(),
        fragment_index: 1,
        repair: None,
    };
    assert_eq!(postcard::to_allocvec(&prep).unwrap()[0], 9);
    let run = RunPreparation::MismatchPcrCatalogDigest {
        forward: b.forward,
        reverse: b.reverse,
        options: r.pcr_options,
        selection: b.selection,
        fragment_index: 1,
        repair: None,
    };
    assert_eq!(postcard::to_allocvec(&run).unwrap()[0], 6);
    let inputs = vec![PlanInput {
        id: "source".into(),
        revision: source.revision(),
        molecule: source,
    }];
    let target = Molecule::dna(
        "Independent selected fragment",
        "AATTCACTACCCCGGGGTAACG",
        Topology::Linear,
    )
    .unwrap();
    let recipe = ligation_recipe::Recipe {
        name: "Mismatch PCR then pinned EcoRI".into(),
        fragments: vec![ligation_recipe::Fragment {
            source: FragmentSpec {
                generation: None,
                role: "insert".into(),
                input_id: "source".into(),
                preparation: prep,
                reverse: false,
            },
            treatment: None,
        }],
        circular: false,
        allow_nicked: false,
    };
    let request = Request {
        id: 2,
        revision: ligation_recipe::preview_revision(
            &inputs,
            &target,
            &recipe,
            collagene_core::assembly_recipe::CleavagePolicy::AssumedComplete,
        )
        .unwrap(),
        job: AnalysisJob::PreviewLigationRecipe {
            policy: collagene_core::assembly_recipe::CleavagePolicy::AssumedComplete,
            inputs,
            target: target.clone(),
            recipe,
        },
    };
    let e = protocol::process_request(&wire(b"CGQ\0", 77, &request)).unwrap_err();
    assert!(matches!(e, collagene_core::Error::Unsupported(_)));
    let response = protocol::decode_response(
        &protocol::process_request(&protocol::encode_request(&request).unwrap()).unwrap(),
    )
    .unwrap();
    let AnalysisResult::LigationRecipePreview(p) = &response.result else {
        panic!("preview")
    };
    assert_eq!(p.verification.status, Status::ConsistentUnderAssumptions);
    assert!(protocol::decode_response(&wire(b"CGR\0", 77, &response)).is_err());
    let r = Request {
        id: 3,
        revision: target.revision(),
        job: AnalysisJob::Verify {
            plan: p.compiled.plan.clone(),
            target,
        },
    };
    assert!(protocol::process_request(&wire(b"CGQ\0", 77, &r)).is_err());
    assert!(protocol::process_request(&wire(b"CGQ\0", 91, &r)).is_err());
    assert!(
        p.compiled
            .plan
            .steps
            .iter()
            .any(|o| matches!(o, Operation::PcrWithMismatches { .. }))
    );
}
