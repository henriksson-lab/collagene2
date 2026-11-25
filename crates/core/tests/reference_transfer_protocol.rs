//! Binary copy/paste reference transfer retains annotations without chemical claims.
use collagene_core::{
    Error, Feature, Location, Molecule, Topology,
    protocol::{self, AnalysisJob, AnalysisResult, Request},
    reference_transfer::{self, Part, RequestCopy, RequestPaste, Selection},
};

fn fixture() -> (RequestCopy, Molecule) {
    let mut source = Molecule::dna("source", "AACGTTAC", Topology::Linear).unwrap();
    source.features.push(Feature {
        kind: "misc_feature".into(),
        location: Location::range(3, 5),
        qualifiers: vec![("label".into(), Some("copied source feature".into()))],
        needs_review: false,
    });
    let mut destination = Molecule::dna("destination", "TTAAAACC", Topology::Linear).unwrap();
    destination.features.push(Feature {
        kind: "misc_feature".into(),
        location: Location::range(4, 6),
        qualifiers: vec![("label".into(), Some("existing destination feature".into()))],
        needs_review: false,
    });
    (
        RequestCopy {
            source,
            selection: Selection {
                parts: vec![Part {
                    start: 2,
                    end: 6,
                    reverse: false,
                }],
                reverse: true,
            },
        },
        destination,
    )
}
fn copy_request(copy: RequestCopy) -> Request {
    Request {
        id: 53,
        revision: reference_transfer::copy_revision(&copy).unwrap(),
        job: AnalysisJob::CopyReferenceFragment { request: copy },
    }
}
fn paste_request(copy: &RequestCopy, destination: Molecule, fragment_revision: String) -> Request {
    let paste = RequestPaste {
        source: copy.source.clone(),
        destination,
        selection: copy.selection.clone(),
        fragment_revision,
        start: 2,
        end: 4,
    };
    Request {
        id: 54,
        revision: reference_transfer::paste_revision(&paste).unwrap(),
        job: AnalysisJob::PasteReferenceFragment { request: paste },
    }
}
fn round_trip(request: &Request) -> (Vec<u8>, protocol::Response) {
    let wire = protocol::process_request(&protocol::encode_request(request).unwrap()).unwrap();
    let response = protocol::decode_response(&wire).unwrap();
    assert_eq!(response.id, request.id);
    assert_eq!(response.revision, request.revision);
    (wire, response)
}
fn tag(wire: &[u8]) -> u32 {
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&wire[42..]).unwrap();
    tag
}

#[test]
fn typed_reverse_copy_and_paste_keep_source_and_destination_annotations() {
    let (copy, destination) = fixture();
    let source_before = copy.source.clone();
    let destination_before = destination.clone();
    let request = copy_request(copy.clone());
    let request_wire = protocol::encode_request(&request).unwrap();
    assert_eq!(tag(&request_wire), 120);
    let (reply_wire, reply) = round_trip(&request);
    assert_eq!(tag(&reply_wire), 112);
    let AnalysisResult::ReferenceFragment(fragment) = reply.result else {
        panic!("copy did not return a typed reference fragment");
    };
    fragment.validate_binding(&request.revision).unwrap();
    assert_eq!(fragment.molecule.sequence.as_str(), "AACG");
    assert_eq!(fragment.molecule.topology, Topology::Linear);
    let source_feature = fragment
        .molecule
        .features
        .iter()
        .find(|feature| feature.label() == "copied source feature")
        .unwrap();
    assert_eq!(
        source_feature.location,
        Location::Complement(Box::new(Location::range(1, 3)))
    );
    assert_eq!(
        source_feature
            .location
            .extract(&fragment.molecule.sequence)
            .unwrap(),
        "GT"
    );

    let request = paste_request(&copy, destination, fragment.context_revision);
    assert_eq!(tag(&protocol::encode_request(&request).unwrap()), 121);
    let (reply_wire, reply) = round_trip(&request);
    assert_eq!(tag(&reply_wire), 113);
    let AnalysisResult::ReferencePaste(report) = reply.result else {
        panic!("paste did not return a typed reference report");
    };
    report.validate_binding(&request.revision).unwrap();
    assert_eq!(report.product.sequence.as_str(), "TTAACGAACC");
    assert_eq!((report.inserted_start, report.inserted_end), (2, 6));
    let copied = report
        .product
        .features
        .iter()
        .find(|feature| feature.label() == "copied source feature")
        .unwrap();
    assert_eq!(
        copied.location,
        Location::Complement(Box::new(Location::range(3, 5)))
    );
    assert_eq!(
        copied.location.extract(&report.product.sequence).unwrap(),
        "GT"
    );
    let retained = report
        .product
        .features
        .iter()
        .find(|feature| feature.label() == "existing destination feature")
        .unwrap();
    assert_eq!(retained.location, Location::range(6, 8));
    assert!(!report.product.metadata.contains_key("verification_status"));
    assert_eq!(copy.source, source_before);
    if let AnalysisJob::PasteReferenceFragment { request } = request.job {
        assert_eq!(request.destination, destination_before);
    }
}

#[test]
fn schema53_selectors_are_rejected_before_old_payload_decode_and_bind_both_sources() {
    let (copy, destination) = fixture();
    let copied = reference_transfer::copy(&copy).unwrap();
    let requests = [
        copy_request(copy.clone()),
        paste_request(&copy, destination, copied.context_revision),
    ];
    for (request, expected_job, expected_result) in requests
        .iter()
        .zip([120, 121])
        .zip([112, 113])
        .map(|((r, j), v)| (r, j, v))
    {
        let wire = protocol::encode_request(request).unwrap();
        assert_eq!(tag(&wire), expected_job);
        let mut old = wire.clone();
        old[4..6].copy_from_slice(&52u16.to_le_bytes());
        assert!(matches!(
            protocol::process_request(&old),
            Err(Error::Unsupported(_))
        ));
        // Only the borrowed discriminator prefix exists. Its checksum is
        // deliberately stale: the new selector's minimum version comes first.
        let body =
            postcard::to_allocvec(&(request.id, request.revision.as_str(), expected_job)).unwrap();
        let mut prefix = wire[..42].to_vec();
        prefix[4..6].copy_from_slice(&52u16.to_le_bytes());
        prefix[6..10].copy_from_slice(&(body.len() as u32).to_le_bytes());
        prefix.extend_from_slice(&body);
        assert!(matches!(
            protocol::process_request(&prefix),
            Err(Error::Unsupported(_))
        ));
        let (response_wire, _) = round_trip(request);
        assert_eq!(tag(&response_wire), expected_result);
        let mut old_response = response_wire.clone();
        old_response[4..6].copy_from_slice(&52u16.to_le_bytes());
        assert!(matches!(
            protocol::decode_response(&old_response),
            Err(Error::Unsupported(_))
        ));
        let body = postcard::to_allocvec(&(request.id, request.revision.as_str(), expected_result))
            .unwrap();
        let mut prefix = response_wire[..42].to_vec();
        prefix[4..6].copy_from_slice(&52u16.to_le_bytes());
        prefix[6..10].copy_from_slice(&(body.len() as u32).to_le_bytes());
        prefix.extend_from_slice(&body);
        assert!(matches!(
            protocol::decode_response(&prefix),
            Err(Error::Unsupported(_))
        ));
        let mut stale = request.clone();
        stale.revision = "0".repeat(64);
        assert!(protocol::process_request(&protocol::encode_request(&stale).unwrap()).is_err());
    }
    let mut changed = requests[0].clone();
    if let AnalysisJob::CopyReferenceFragment { request } = &mut changed.job {
        request.source.name = "changed original source".into();
    }
    assert!(protocol::process_request(&protocol::encode_request(&changed).unwrap()).is_err());
    let mut changed = requests[1].clone();
    if let AnalysisJob::PasteReferenceFragment { request } = &mut changed.job {
        request.destination.name = "changed destination".into();
    }
    assert!(protocol::process_request(&protocol::encode_request(&changed).unwrap()).is_err());
}
