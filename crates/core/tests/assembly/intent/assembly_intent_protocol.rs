use collagene_core::{
    cloning::assembly::intent::{
        self as assembly_intent, Contribution, Request as IntentRequest, Selector,
    },
    planning::plan::PlanInput,
    protocol::{
        self,
        request::Request,
        schema::{AnalysisJob, AnalysisResult},
    },
};
use collagene_model::{Molecule, Topology};
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
        fragments: vec![Contribution {
            role: "backbone".into(),
            input_id: "source-1".into(),
            selector: Selector::ReferenceRange {
                start: 0,
                end: 4,
                reverse: false,
            },
        }],
        circular: false,
        material_bindings: vec![],
    };
    let revision = assembly_intent::intent_revision(&inputs, &intent, None).unwrap();
    let request = Request {
        id: 998,
        revision: revision.clone(),
        job: AnalysisJob::BuildAssemblyIntent {
            inputs: inputs.clone(),
            request: intent.clone(),
            workspace: None,
        },
    };
    let bytes = protocol::request::Request::encode(&request).unwrap();
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[42..]).unwrap();
    assert_eq!(tag, 96);
    let response = protocol::execution::process_request(&bytes).unwrap();
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&response[42..]).unwrap();
    assert_eq!(tag, 88);
    let AnalysisResult::AssemblyIntent(report) = protocol::response::Response::decode(&response)
        .unwrap()
        .result
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
        assert!(protocol::execution::process_request(&downgraded).is_err());
        let mut downgraded = response.clone();
        downgraded[4..6].copy_from_slice(&old.to_le_bytes());
        assert!(protocol::response::Response::decode(&downgraded).is_err());
    }
    let mut changed = intent;
    if let Selector::ReferenceRange { end, .. } = &mut changed.fragments[0].selector {
        *end = 5;
    } else {
        panic!("reference selector required by fixture");
    }
    let bytes = protocol::request::Request::encode(&Request {
        id: 998,
        revision,
        job: AnalysisJob::BuildAssemblyIntent {
            inputs,
            request: changed,
            workspace: None,
        },
    })
    .unwrap();
    assert!(protocol::execution::process_request(&bytes).is_err());
}
