//! Appended worker contracts preserve explicit profiles, reciprocal alternatives and independent target requirements.
use collagene_core::{
    formats,
    protocol::{self, AnalysisJob as J, AnalysisResult as R, Request, Response},
    recombination::{self as a, reference},
};
fn fixture() -> (Vec<reference::Input>, a::Request, collagene_core::Molecule) {
    let m = formats::import_genbank(include_str!("../../../fixtures/recombination-bp-demo.gb"))
        .unwrap()
        .molecules;
    (
        vec![
            reference::Input {
                role: "insert".into(),
                molecule: m[0].clone(),
            },
            reference::Input {
                role: "donor".into(),
                molecule: m[1].clone(),
            },
        ],
        a::Request::new(
            a::ProfileSnapshot::manufacturer_context().unwrap(),
            reference::Reaction::BP,
        ),
        m[2].clone(),
    )
}
fn raw(magic: &[u8; 4], v: u16, p: &[u8]) -> Vec<u8> {
    let mut b = magic.to_vec();
    b.extend(v.to_le_bytes());
    b.extend((p.len() as u32).to_le_bytes());
    b.extend(blake3::hash(p).as_bytes());
    b.extend(p);
    b
}
fn run(revision: String, job: J) -> Response {
    let q = Request {
        id: 701,
        revision: revision.clone(),
        job,
    };
    let b = protocol::encode_request(&q).unwrap();
    assert_eq!(&b[4..6], &23u16.to_le_bytes());
    let r = protocol::decode_response(&protocol::process_request(&b).unwrap()).unwrap();
    assert_eq!(r.id, 701);
    assert_eq!(r.revision, revision);
    r
}
#[test]
fn inspect_check_select_reconstruct_independent_literal_fixture() {
    let (inputs, request, target) = fixture();
    let policy = reference::TargetPolicy {
        allow_reverse_complement: false,
    };
    let R::RecombinationInspection(report) = run(
        a::batch_revision(&inputs, &request).unwrap(),
        J::InspectRecombination {
            inputs: inputs.clone(),
            request: request.clone(),
        },
    )
    .result
    else {
        panic!()
    };
    report.validate_bindings_for(&inputs, &request).unwrap();
    assert!(report.reference.complete);
    assert_eq!(report.documents.len(), 2);
    assert_eq!(report.reference.candidates.len(), 1);
    let R::RecombinationTargetCheck(check) = run(
        a::check_revision(&inputs, &request, &target, &policy).unwrap(),
        J::CheckRecombinationTarget {
            inputs: inputs.clone(),
            request: request.clone(),
            target: target.clone(),
            policy,
        },
    )
    .result
    else {
        panic!()
    };
    assert_eq!(check.reference.matches.len(), 1);
    let matched = &check.reference.matches[0];
    assert!(matched.fully_specified_dna);
    let selection = a::Selection {
        batch_revision: report.batch_revision.clone(),
        candidate_id: matched.candidate,
        component_id: matched.product,
    };
    let R::RecombinationProduct(checked) = run(
        a::checked_product_revision(&inputs, &request, &selection, &target, &policy).unwrap(),
        J::SelectRecombinationComponent {
            inputs: inputs.clone(),
            request: request.clone(),
            selection: selection.clone(),
            target: target.clone(),
            policy,
        },
    )
    .result
    else {
        panic!()
    };
    checked
        .validate_bindings_for(&inputs, &request, &selection, &target, &policy)
        .unwrap();
    assert_eq!(checked.product.reciprocal_components.len(), 2);
    assert!(collagene_core::chemistry::Duplex::from_molecule(&checked.product.molecule).is_err());
    let wrong = collagene_core::Molecule::dna(
        "independent wrong",
        "ACGT",
        collagene_core::Topology::Circular,
    )
    .unwrap();
    let r = run(
        a::checked_product_revision(&inputs, &request, &selection, &wrong, &policy).unwrap(),
        J::SelectRecombinationComponent {
            inputs,
            request,
            selection,
            target: wrong,
            policy,
        },
    );
    assert!(matches!(r.result, R::Error(_)));
}
#[test]
fn old_headers_source_role_order_and_foreign_responses_are_rejected() {
    let (inputs, request, target) = fixture();
    let rev = a::batch_revision(&inputs, &request).unwrap();
    let q = Request {
        id: 3,
        revision: rev.clone(),
        job: J::InspectRecombination {
            inputs: inputs.clone(),
            request: request.clone(),
        },
    };
    let payload = postcard::to_allocvec(&q).unwrap();
    for v in 1..22 {
        assert!(
            protocol::process_request(&raw(b"CGQ\0", v, &payload)).is_err(),
            "version {v}"
        );
    }
    let report = a::simulate(&inputs, &request).unwrap();
    let r = Response {
        id: 3,
        revision: rev.clone(),
        result: R::RecombinationInspection(report),
    };
    let payload = postcard::to_allocvec(&r).unwrap();
    for v in 1..22 {
        assert!(protocol::decode_response(&raw(b"CGR\0", v, &payload)).is_err());
    }
    let mut swapped = inputs.clone();
    swapped.swap(0, 1);
    assert_ne!(a::batch_revision(&swapped, &request).unwrap(), rev);
    let mut renamed = inputs.clone();
    renamed[0].role = "equal DNA different role".into();
    assert!(
        protocol::process_request(
            &protocol::encode_request(&Request {
                id: 3,
                revision: rev,
                job: J::InspectRecombination {
                    inputs: renamed,
                    request: request.clone()
                }
            })
            .unwrap()
        )
        .is_err()
    );
    let mut r = r;
    r.revision = "0".repeat(64);
    assert!(
        protocol::decode_response(&raw(b"CGR\0", 22, &postcard::to_allocvec(&r).unwrap())).is_err()
    );
    let mut forged = request;
    forged.profile.revision = "0".repeat(64);
    assert!(
        a::check_revision(
            &inputs,
            &forged,
            &target,
            &reference::TargetPolicy {
                allow_reverse_complement: false
            }
        )
        .is_err()
    );
}
