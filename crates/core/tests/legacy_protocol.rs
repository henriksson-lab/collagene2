//! Independently defined historical layouts and hand-encoded Postcard fixtures.
//! Changing only today's envelope version does not prove backward compatibility.
use collagene_core::{
    analysis::Primer,
    protocol::*,
    verification::{DesignPlan, Operation, Source, Status, VerificationReport},
    *,
};
use serde::Serialize;

fn envelope(magic: &[u8; 4], version: u16, payload: &[u8]) -> Vec<u8> {
    let mut bytes = magic.to_vec();
    bytes.extend_from_slice(&version.to_le_bytes());
    bytes.extend_from_slice(&(payload.len() as u32).to_le_bytes());
    bytes.extend_from_slice(blake3::hash(payload).as_bytes());
    bytes.extend_from_slice(payload);
    bytes
}
fn wire(magic: &[u8; 4], version: u16, value: &impl Serialize) -> Vec<u8> {
    envelope(magic, version, &postcard::to_allocvec(value).unwrap())
}
fn molecule() -> Molecule {
    Molecule::dna("old product", "AACCGGTT", Topology::Circular).unwrap()
}

#[derive(Serialize)]
enum OldStatus {
    Invalid,
    ConsistentUnderAssumptions,
    Incomplete,
}
#[derive(Serialize)]
struct OldReport {
    plan_revision: String,
    target_revision: String,
    status: OldStatus,
    findings: Vec<String>,
    assumptions: Vec<String>,
    products: Vec<Molecule>,
    exhaustive_alternative_search: bool,
}
#[derive(Serialize)]
struct OldMatch {
    start: usize,
    length: usize,
    reverse: bool,
    certain: bool,
}
#[derive(Serialize)]
struct OldCut {
    enzyme: String,
    recognition_start: usize,
    reverse: bool,
    top_cut: usize,
    bottom_cut: usize,
    certain: bool,
}
#[derive(Serialize)]
struct OldPcr {
    sequence: Sequence,
    forward_start: usize,
    reverse_start: usize,
    crosses_origin: bool,
    forward_primer: String,
    reverse_primer: String,
}
#[derive(Serialize)]
enum OldResult {
    Search(Vec<OldMatch>),
    Restriction(Vec<OldCut>),
    Error(String),
    Pcr(Vec<OldPcr>),
    Verification(OldReport),
    DuplexDigest(Vec<collagene_core::chemistry::DigestFragment>),
}
#[derive(Serialize)]
struct OldResponse {
    id: u32,
    revision: String,
    result: OldResult,
}
fn old_response(result: OldResult) -> OldResponse {
    OldResponse {
        id: 301,
        revision: "source revision".into(),
        result,
    }
}

#[test]
fn hand_encoded_verification_without_engine_field_migrates_versions_one_and_two() {
    // id=17, revision="r", result discriminant4; report p,t,status1,
    // findings=["ok"], assumptions=["a"], products=[], exhaustive=false.
    let payload = [
        17, 1, b'r', 4, 1, b'p', 1, b't', 1, 1, 2, b'o', b'k', 1, 1, b'a', 0, 0,
    ];
    for version in [1, 2] {
        let decoded = decode_response(&envelope(b"CGR\0", version, &payload)).unwrap();
        assert_eq!(decoded.id, 17);
        assert_eq!(decoded.revision, "r");
        let AnalysisResult::Verification(report) = decoded.result else {
            panic!("historical discriminator changed")
        };
        assert_eq!(report.plan_revision, "p");
        assert_eq!(report.target_revision, "t");
        assert_eq!(report.status, Status::ConsistentUnderAssumptions);
        assert_eq!(report.findings, vec!["ok"]);
        assert_eq!(report.assumptions, vec!["a"]);
        assert!(report.products.is_empty());
        assert!(!report.exhaustive_alternative_search);
        assert_eq!(report.engine_version, LEGACY_ENGINE_VERSION);
        assert_ne!(
            report.engine_version,
            collagene_core::verification::ENGINE_VERSION
        );
    }
}
#[test]
fn historical_verification_preserves_every_field_and_product_evidence() {
    let mut product = molecule();
    product.features.push(Feature {
        kind: "misc_feature".into(),
        location: Location::Complement(Box::new(Location::Join(vec![
            Location::range(6, 8),
            Location::range(0, 2),
        ]))),
        qualifiers: vec![("note".into(), Some("original\nsource evidence".into()))],
        needs_review: true,
    });
    product
        .metadata
        .insert("legacy:unknown".into(), "retained".into());
    let expected_product = serde_json::to_value(&product).unwrap();
    for version in [1, 2] {
        for (old_status, status) in [
            (OldStatus::Invalid, Status::Invalid),
            (
                OldStatus::ConsistentUnderAssumptions,
                Status::ConsistentUnderAssumptions,
            ),
            (OldStatus::Incomplete, Status::Incomplete),
        ] {
            let old = old_response(OldResult::Verification(OldReport {
                plan_revision: "a".repeat(64),
                target_revision: product.revision(),
                status: old_status,
                findings: vec!["first finding".into(), "μ source".into()],
                assumptions: vec!["old assumption".into()],
                products: vec![product.clone()],
                exhaustive_alternative_search: true,
            }));
            let decoded = decode_response(&wire(b"CGR\0", version, &old)).unwrap();
            assert_eq!(decoded.id, 301);
            assert_eq!(decoded.revision, "source revision");
            let AnalysisResult::Verification(report) = decoded.result else {
                panic!()
            };
            assert_eq!(report.plan_revision, "a".repeat(64));
            assert_eq!(report.target_revision, product.revision());
            assert_eq!(report.status, status);
            assert_eq!(report.findings, vec!["first finding", "μ source"]);
            assert_eq!(report.assumptions, vec!["old assumption"]);
            assert!(report.exhaustive_alternative_search);
            assert_eq!(report.products.len(), 1);
            assert_eq!(
                serde_json::to_value(&report.products[0]).unwrap(),
                expected_product
            );
            assert_eq!(report.engine_version, LEGACY_ENGINE_VERSION);
        }
    }
}
#[test]
fn original_nonverification_response_layouts_keep_their_discriminants_and_fields() {
    for version in [1, 2] {
        let cases = [
            OldResult::Search(vec![OldMatch {
                start: 6,
                length: 4,
                reverse: true,
                certain: false,
            }]),
            OldResult::Restriction(vec![OldCut {
                enzyme: "manual".into(),
                recognition_start: 7,
                reverse: true,
                top_cut: 3,
                bottom_cut: 5,
                certain: true,
            }]),
            OldResult::Error("old scientific error".into()),
            OldResult::Pcr(vec![OldPcr {
                sequence: Sequence::dna("TTACGTAA").unwrap(),
                forward_start: 7,
                reverse_start: 2,
                crosses_origin: true,
                forward_primer: "left".into(),
                reverse_primer: "right".into(),
            }]),
        ];
        for old in cases {
            let report = decode_response(&wire(b"CGR\0", version, &old_response(old))).unwrap();
            match report.result {
                AnalysisResult::Search(hits) => {
                    let h = &hits[0];
                    assert_eq!(
                        (h.start, h.length, h.reverse, h.certain),
                        (6, 4, true, false)
                    );
                }
                AnalysisResult::Restriction(cuts) => {
                    let c = &cuts[0];
                    assert_eq!(
                        (
                            &*c.enzyme,
                            c.recognition_start,
                            c.reverse,
                            c.top_cut,
                            c.bottom_cut,
                            c.certain
                        ),
                        ("manual", 7, true, 3, 5, true)
                    );
                }
                AnalysisResult::Error(message) => assert_eq!(message, "old scientific error"),
                AnalysisResult::Pcr(products) => {
                    let p = &products[0];
                    assert_eq!(p.sequence.as_str(), "TTACGTAA");
                    assert_eq!(
                        (p.forward_start, p.reverse_start, p.crosses_origin),
                        (7, 2, true)
                    );
                    assert_eq!((&*p.forward_primer, &*p.reverse_primer), ("left", "right"));
                }
                _ => panic!("legacy enum order changed"),
            }
        }
    }
    let bytes = wire(b"CGR\0", 2, &old_response(OldResult::DuplexDigest(vec![])));
    assert!(
        matches!(decode_response(&bytes).unwrap().result,AnalysisResult::DuplexDigest(f)if f.is_empty())
    );
    let mut mislabeled = bytes;
    mislabeled[4..6].copy_from_slice(&1u16.to_le_bytes());
    assert!(decode_response(&mislabeled).is_err());
}
#[test]
fn version_three_report_with_recorded_engine_is_not_migrated_as_unrecorded() {
    let report = VerificationReport {
        plan_revision: "plan".into(),
        target_revision: "target".into(),
        engine_version: "historical-engine/science-3".into(),
        status: Status::Incomplete,
        findings: vec!["retained".into()],
        assumptions: vec![],
        products: vec![molecule()],
        exhaustive_alternative_search: false,
    };
    let bytes = wire(
        b"CGR\0",
        3,
        &Response {
            id: 19,
            revision: "revision".into(),
            result: AnalysisResult::Verification(report.clone()),
        },
    );
    let response = decode_response(&bytes).unwrap();
    let AnalysisResult::Verification(decoded) = response.result else {
        panic!()
    };
    assert_eq!(
        serde_json::to_value(decoded).unwrap(),
        serde_json::to_value(report).unwrap()
    );
    for wrong_version in [1u16, 2] {
        let mut mislabeled = bytes.clone();
        mislabeled[4..6].copy_from_slice(&wrong_version.to_le_bytes());
        assert!(decode_response(&mislabeled).is_err());
    }
    let old = old_response(OldResult::Verification(OldReport {
        plan_revision: "p".into(),
        target_revision: "t".into(),
        status: OldStatus::Invalid,
        findings: vec![],
        assumptions: vec![],
        products: vec![],
        exhaustive_alternative_search: false,
    }));
    assert!(decode_response(&wire(b"CGR\0", 3, &old)).is_err());
}

#[allow(dead_code)]
#[derive(Serialize)]
enum OldSource {
    Input(String),
    Step(usize),
}
#[allow(dead_code)]
#[derive(Serialize)]
enum OldOperation {
    Pcr {
        template: OldSource,
        forward: Primer,
        reverse: Primer,
    },
    Overlap {
        fragments: Vec<OldSource>,
        overlaps: Vec<usize>,
        circular: bool,
    },
}
#[derive(Serialize)]
struct OldInput {
    id: String,
    revision: String,
    molecule: Molecule,
}
#[derive(Serialize)]
struct OldPlan {
    name: String,
    target_revision: String,
    inputs: Vec<OldInput>,
    steps: Vec<OldOperation>,
}
#[allow(dead_code)]
#[derive(Serialize)]
enum OldJob {
    Search {
        molecule: Molecule,
        query: String,
    },
    Restriction {
        molecule: Molecule,
        enzyme: collagene_core::analysis::Enzyme,
    },
    Pcr {
        molecule: Molecule,
        forward: Primer,
        reverse: Primer,
    },
    Verify {
        target: Molecule,
        plan: OldPlan,
    },
}
#[derive(Serialize)]
struct OldRequest {
    id: u32,
    revision: String,
    job: OldJob,
}
#[test]
fn genuinely_old_verify_request_layout_reconstructs_without_new_operation_variants() {
    let left = Molecule::dna("left", "AAAACCCC", Topology::Linear).unwrap();
    let right = Molecule::dna("right", "CCCCGGGG", Topology::Linear).unwrap();
    let target = Molecule::dna("target", "AAAACCCCGGGG", Topology::Linear).unwrap();
    let old = OldRequest {
        id: 81,
        revision: target.revision(),
        job: OldJob::Verify {
            target: target.clone(),
            plan: OldPlan {
                name: "old overlap".into(),
                target_revision: target.revision(),
                inputs: vec![
                    OldInput {
                        id: "left".into(),
                        revision: left.revision(),
                        molecule: left,
                    },
                    OldInput {
                        id: "right".into(),
                        revision: right.revision(),
                        molecule: right,
                    },
                ],
                steps: vec![OldOperation::Overlap {
                    fragments: vec![
                        OldSource::Input("left".into()),
                        OldSource::Input("right".into()),
                    ],
                    overlaps: vec![4],
                    circular: false,
                }],
            },
        },
    };
    for version in [1, 2] {
        let bytes = wire(b"CGQ\0", version, &old);
        let reply = process_request(&bytes).unwrap();
        let decoded = decode_response(&reply).unwrap();
        assert_eq!(decoded.id, 81);
        let AnalysisResult::Verification(report) = decoded.result else {
            panic!()
        };
        assert_eq!(report.status, Status::ConsistentUnderAssumptions);
        assert_eq!(report.products[0].sequence, target.sequence);
        assert_eq!(
            report.engine_version,
            collagene_core::verification::ENGINE_VERSION
        );
    }
}
#[test]
fn schema_one_verify_rejects_each_new_physical_operation_before_execution() {
    let target = molecule();
    for operation in [
        Operation::Digest {
            template: Source::Input("missing".into()),
            enzymes: vec!["EcoRI".into()],
            fragment_index: 0,
        },
        Operation::TreatEnds {
            template: Source::Input("missing".into()),
            treatment: collagene_core::chemistry::EndTreatment::Phosphorylate,
        },
        Operation::Ligate {
            fragments: vec![],
            circular: true,
            allow_nicked: false,
        },
    ] {
        let request = Request {
            id: 9,
            revision: target.revision(),
            job: AnalysisJob::Verify {
                target: target.clone(),
                plan: DesignPlan {
                    name: "physical v2".into(),
                    target_revision: target.revision(),
                    inputs: vec![],
                    steps: vec![operation],
                },
            },
        };
        let bytes = wire(b"CGQ\0", 1, &request);
        assert!(
            matches!(process_request(&bytes),Err(Error::Unsupported(message))if message.contains("schema 2"))
        );
        assert!(process_request(&wire(b"CGQ\0", 2, &request)).is_ok());
    }
}
#[test]
fn every_truncated_legacy_payload_and_header_fails_without_panicking() {
    let payload = [
        17, 1, b'r', 4, 1, b'p', 1, b't', 1, 1, 2, b'o', b'k', 1, 1, b'a', 0, 0,
    ];
    for version in [1, 2] {
        let full = envelope(b"CGR\0", version, &payload);
        for length in 0..full.len() {
            assert!(
                std::panic::catch_unwind(|| decode_response(&full[..length]))
                    .unwrap()
                    .is_err()
            );
        }
        for length in 0..payload.len() {
            let short = envelope(b"CGR\0", version, &payload[..length]);
            assert!(
                std::panic::catch_unwind(|| decode_response(&short))
                    .unwrap()
                    .is_err()
            );
        }
        let mut trailing = payload.to_vec();
        trailing.push(0);
        assert!(decode_response(&envelope(b"CGR\0", version, &trailing)).is_err());
    }
    for length in 0..42 {
        assert!(
            std::panic::catch_unwind(|| process_request(&vec![0; length]))
                .unwrap()
                .is_err()
        );
    }
}
