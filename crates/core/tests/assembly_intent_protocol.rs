use collagene_core::{
    Molecule, Topology,
    assembly_intent::{self, Request as IntentRequest},
    overlap_primer_design::Fragment,
    protocol::{self, AnalysisJob, AnalysisResult, Request},
    verification::PlanInput,
};
#[test]
fn declared_reference_transport_binds_source_and_regions_and_requires_its_new_header() {
    let molecule = Molecule::dna("original", "AACCGGTT", Topology::Linear).unwrap();
    let inputs = vec![PlanInput {
        id: "source-1".into(),
        revision: molecule.revision(),
        molecule,
    }];
    let intent = IntentRequest {
        name: "selected prefix".into(),
        fragments: vec![Fragment {
            role: "backbone".into(),
            input_id: "source-1".into(),
            start: 0,
            end: 4,
            reverse: false,
        }],
        circular: false,
    };
    let revision = assembly_intent::intent_revision(&inputs, &intent).unwrap();
    let request = Request {
        id: 998,
        revision: revision.clone(),
        job: AnalysisJob::BuildAssemblyIntent {
            inputs: inputs.clone(),
            request: intent.clone(),
        },
    };
    let bytes = protocol::encode_request(&request).unwrap();
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[42..]).unwrap();
    assert_eq!(tag, 98);
    let response = protocol::process_request(&bytes).unwrap();
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&response[42..]).unwrap();
    assert_eq!(tag, 90);
    let AnalysisResult::AssemblyIntent(report) =
        protocol::decode_response(&response).unwrap().result
    else {
        panic!("expected declared intent, not a verification result")
    };
    assert_eq!(report.target.sequence.as_str(), "AACC");
    assert_eq!(report.context_revision, revision);
    assert_eq!(report.fragments[0].output_start, 0);
    assert_eq!(report.fragments[0].output_end, 4);
    for old in [1u16, 35] {
        let mut downgraded = bytes.clone();
        downgraded[4..6].copy_from_slice(&old.to_le_bytes());
        assert!(protocol::process_request(&downgraded).is_err());
        let mut downgraded = response.clone();
        downgraded[4..6].copy_from_slice(&old.to_le_bytes());
        assert!(protocol::decode_response(&downgraded).is_err());
    }
    let mut changed = intent;
    changed.fragments[0].end = 5;
    let bytes = protocol::encode_request(&Request {
        id: 998,
        revision,
        job: AnalysisJob::BuildAssemblyIntent {
            inputs,
            request: changed,
        },
    })
    .unwrap();
    assert!(protocol::process_request(&bytes).is_err());
}
