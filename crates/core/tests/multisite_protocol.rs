//! Current worker dispatch preserves source roles and independent reference intent.
use collagene_core::{
    Molecule, Topology,
    multisite_recombination::{self as adapter, AnnotationLimits, Request, Selection, multisite},
    protocol,
    recombination::reference::{Family, Input, Limits, Reaction, SiteVariant, SourceRecord},
};
use serde::Serialize;
fn core(id: u8) -> &'static str {
    match id {
        1 => "GTACAAA",
        2 => "GTACAAG",
        3 => "GTACAAC",
        4 => "GTACAAT",
        _ => panic!("fixture specificity"),
    }
}
fn site(family: Family, id: u8) -> String {
    match family {
        Family::L => format!("ACGG{}TCC", core(id)),
        Family::R => format!("TGGA{}GCT", core(id)),
        _ => panic!("fixture family"),
    }
}
fn exchange_fixture() -> (Vec<Input>, Request) {
    let mut profile=multisite::Profile{id:"synthetic-four-parent-LR".into(),edition:"independent-literal-v1".into(),sources:vec![],variants:vec![],limitations:vec!["Synthetic reference topology fixture, not measured biological specificity or a manufacturer profile".into()]};
    for id in 1..=4 {
        for (family, prefix) in [(Family::L, "L"), (Family::R, "R")] {
            let key = format!("synthetic-{prefix}{id}");
            let sequence = site(family, id);
            profile.sources.push(SourceRecord {
                id: key.clone(),
                url: "https://example.org/synthetic-reference-fixture".into(),
                raw_sha256: "1".repeat(64),
                vector_name: "Explicit synthetic bookkeeping fixture".into(),
                vector_length: sequence.len(),
                start_1based: 1,
                end_1based: sequence.len(),
                reverse_extraction: false,
                map_url: "https://example.org/synthetic-reference-fixture-map".into(),
                map_sha256: "2".repeat(64),
            });
            profile.variants.push(SiteVariant {
                id: key.clone(),
                family,
                subtype: id,
                sequence,
                overlap_start: 4,
                source_id: key,
            });
        }
    }
    let dna = [
        format!("{}CCCCA{}AAAAC", site(Family::L, 4), site(Family::R, 1)),
        format!("{}GGGGT{}TTTTC", site(Family::L, 1), site(Family::L, 2)),
        format!("{}AAAAG{}CCCCT", site(Family::R, 2), site(Family::L, 3)),
        format!("{}TTTTA{}GGGGA", site(Family::R, 4), site(Family::R, 3)),
    ];
    (
        dna.iter()
            .enumerate()
            .map(|(i, dna)| Input {
                role: format!("original-{i}"),
                molecule: Molecule::dna(
                    format!("Original circular parent{i}"),
                    dna,
                    Topology::Circular,
                )
                .unwrap(),
            })
            .collect(),
        Request {
            reaction: multisite::Request {
                profile,
                reaction: Reaction::LR,
                specificities: vec![4, 1, 2, 3],
                limits: Limits::default(),
            },
            annotation_limits: AnnotationLimits::default(),
        },
    )
}

const MAIN: &str = concat!(
    "GTACAATTCCCCCCATGGA",
    "GTACAAATCCGGGGTACGG",
    "GTACAAGGCTAAAAGACGG",
    "GTACAACGCTGGGGATGGA"
);
fn same_circle(actual: &str) -> bool {
    actual.len() == MAIN.len() && format!("{actual}{actual}").contains(MAIN)
}
fn wire(magic: &[u8; 4], version: u16, value: &impl Serialize) -> Vec<u8> {
    let payload = postcard::to_allocvec(value).unwrap();
    let mut out = magic.to_vec();
    out.extend(version.to_le_bytes());
    out.extend((payload.len() as u32).to_le_bytes());
    out.extend(blake3::hash(&payload).as_bytes());
    out.extend(payload);
    out
}
fn dispatch(request: &protocol::Request, tag: u32) -> protocol::Response {
    let bytes = protocol::encode_request(request).unwrap();
    let ((_, _, actual), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[42..]).unwrap();
    assert_eq!(actual, tag);
    protocol::decode_response(&protocol::process_request(&bytes).unwrap()).unwrap()
}
fn selected(inputs: &[Input], request: &Request) -> Selection {
    let batch = adapter::simulate(inputs, request).unwrap();
    batch
        .documents
        .iter()
        .find(|d| same_circle(d.molecule.sequence.as_str()))
        .expect("independent literal main76bp is a retained component")
        .selection
        .clone()
}
#[test]
fn simulate_job152_returns_source_bound_batch144_and_rejects_stale_sources() {
    let (inputs, request) = exchange_fixture();
    let revision = adapter::revision(&inputs, &request).unwrap();
    let request = protocol::Request {
        id: 152,
        revision: revision.clone(),
        job: protocol::AnalysisJob::SimulateMultiSiteRecombination {
            inputs: inputs.clone(),
            request,
        },
    };
    let response = dispatch(&request, 152);
    assert_eq!(response.id, 152);
    assert_eq!(response.revision, revision);
    let protocol::AnalysisResult::MultiSiteRecombination(batch) = response.result else {
        panic!("expected batch144")
    };
    batch.validate_binding(&revision).unwrap();
    assert_eq!(batch.sources.len(), 4);
    assert_eq!(batch.reference.candidates.len(), 1);
    assert_eq!(batch.documents.len(), 2);
    assert!(
        batch
            .documents
            .iter()
            .any(|d| same_circle(d.molecule.sequence.as_str()))
    );
    assert_eq!(
        batch.reference.input_revisions,
        inputs
            .iter()
            .map(|i| i.molecule.revision())
            .collect::<Vec<_>>()
    );
    let result = protocol::Response {
        id: 152,
        revision: revision.clone(),
        result: protocol::AnalysisResult::MultiSiteRecombination(batch),
    };
    let payload = wire(b"CGR\0", protocol::WORKER_SCHEMA_VERSION, &result);
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&payload[42..]).unwrap();
    assert_eq!(tag, 144);
    assert!(protocol::decode_response(&payload).is_ok());
    let mut stale = request;
    stale.revision = "0".repeat(64);
    assert!(protocol::process_request(&protocol::encode_request(&stale).unwrap()).is_err());
    let mut stale = result;
    stale.revision = "0".repeat(64);
    assert!(
        protocol::decode_response(&wire(b"CGR\0", protocol::WORKER_SCHEMA_VERSION, &stale))
            .is_err()
    );
}
#[test]
fn project_job153_result145_preserves_explicit_component_and_reference_only_boundary() {
    let (inputs, request) = exchange_fixture();
    let selection = selected(&inputs, &request);
    let revision = adapter::product_revision(&inputs, &request, &selection).unwrap();
    let req = protocol::Request {
        id: 153,
        revision: revision.clone(),
        job: protocol::AnalysisJob::ProjectMultiSiteRecombination {
            inputs,
            request,
            selection: selection.clone(),
        },
    };
    let response = dispatch(&req, 153);
    let protocol::AnalysisResult::MultiSiteRecombinationProduct(product) = response.result else {
        panic!("expected product145")
    };
    assert_eq!(product.selection, selection);
    assert_eq!(product.context_revision().unwrap(), revision);
    assert!(same_circle(product.molecule.sequence.as_str()));
    assert!(collagene_core::chemistry::has_unmapped_physical_evidence(
        &product.molecule
    ));
    let result = protocol::Response {
        id: 153,
        revision,
        result: protocol::AnalysisResult::MultiSiteRecombinationProduct(product),
    };
    let payload = wire(b"CGR\0", protocol::WORKER_SCHEMA_VERSION, &result);
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&payload[42..]).unwrap();
    assert_eq!(tag, 145);
    assert!(protocol::decode_response(&payload).is_ok());
    let mut stale = result;
    stale.revision = "0".repeat(64);
    assert!(
        protocol::decode_response(&wire(b"CGR\0", protocol::WORKER_SCHEMA_VERSION, &stale))
            .is_err()
    );
}
#[test]
fn check_job154_result146_binds_independent_target_and_never_accepts_wrong_target_revision() {
    let (inputs, request) = exchange_fixture();
    let selection = selected(&inputs, &request);
    let target = Molecule::dna("Independent literal76bp target", MAIN, Topology::Circular).unwrap();
    let revision = adapter::check_revision(&inputs, &request, &selection, &target).unwrap();
    let req = protocol::Request {
        id: 154,
        revision: revision.clone(),
        job: protocol::AnalysisJob::CheckMultiSiteRecombination {
            inputs: inputs.clone(),
            request: request.clone(),
            selection: selection.clone(),
            target: target.clone(),
        },
    };
    let response = dispatch(&req, 154);
    let protocol::AnalysisResult::MultiSiteRecombinationChecked(checked) = response.result else {
        panic!("expected checked146")
    };
    assert!(checked.assessment.target_matches);
    assert_eq!(checked.assessment.target_revision, target.revision());
    assert_eq!(checked.context_revision().unwrap(), revision);
    let response = protocol::Response {
        id: 154,
        revision: revision.clone(),
        result: protocol::AnalysisResult::MultiSiteRecombinationChecked(checked),
    };
    let payload = wire(b"CGR\0", protocol::WORKER_SCHEMA_VERSION, &response);
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&payload[42..]).unwrap();
    assert_eq!(tag, 146);
    assert!(protocol::decode_response(&payload).is_ok());
    let wrong = Molecule::dna(
        "Independent wrong target",
        &"A".repeat(76),
        Topology::Circular,
    )
    .unwrap();
    let wrong_revision = adapter::check_revision(&inputs, &request, &selection, &wrong).unwrap();
    assert_ne!(wrong_revision, revision);
    let mut wrong_req = protocol::Request {
        id: 155,
        revision: wrong_revision.clone(),
        job: protocol::AnalysisJob::CheckMultiSiteRecombination {
            inputs,
            request,
            selection,
            target: wrong,
        },
    };
    let response = dispatch(&wrong_req, 154);
    let protocol::AnalysisResult::MultiSiteRecombinationChecked(checked) = response.result else {
        panic!("expected checked wrong-target evidence")
    };
    assert!(!checked.assessment.target_matches);
    wrong_req.revision = revision.clone();
    assert!(protocol::process_request(&protocol::encode_request(&wrong_req).unwrap()).is_err());
    let changed = protocol::Response {
        id: 154,
        revision: wrong_revision,
        result: protocol::AnalysisResult::MultiSiteRecombinationChecked(checked),
    };
    assert!(
        protocol::decode_response(&wire(b"CGR\0", protocol::WORKER_SCHEMA_VERSION, &changed))
            .is_ok()
    );
    let mut changed = changed;
    changed.revision = revision;
    assert!(
        protocol::decode_response(&wire(b"CGR\0", protocol::WORKER_SCHEMA_VERSION, &changed))
            .is_err()
    );
}
#[test]
fn new_worker_jobs_and_results_cannot_be_relabelled_as_worker67() {
    let (inputs, request) = exchange_fixture();
    let selection = selected(&inputs, &request);
    let target = Molecule::dna("Independent literal target", MAIN, Topology::Circular).unwrap();
    let cases = [
        protocol::Request {
            id: 1,
            revision: adapter::revision(&inputs, &request).unwrap(),
            job: protocol::AnalysisJob::SimulateMultiSiteRecombination {
                inputs: inputs.clone(),
                request: request.clone(),
            },
        },
        protocol::Request {
            id: 2,
            revision: adapter::product_revision(&inputs, &request, &selection).unwrap(),
            job: protocol::AnalysisJob::ProjectMultiSiteRecombination {
                inputs: inputs.clone(),
                request: request.clone(),
                selection: selection.clone(),
            },
        },
        protocol::Request {
            id: 3,
            revision: adapter::check_revision(&inputs, &request, &selection, &target).unwrap(),
            job: protocol::AnalysisJob::CheckMultiSiteRecombination {
                inputs,
                request,
                selection,
                target,
            },
        },
    ];
    for request in cases {
        let current = protocol::encode_request(&request).unwrap();
        let response = protocol::process_request(&current).unwrap();
        let mut old_request = current;
        old_request[4..6].copy_from_slice(&67u16.to_le_bytes());
        assert!(protocol::process_request(&old_request).is_err());
        let mut old_result = response;
        old_result[4..6].copy_from_slice(&67u16.to_le_bytes());
        assert!(protocol::decode_response(&old_result).is_err());
    }
}
