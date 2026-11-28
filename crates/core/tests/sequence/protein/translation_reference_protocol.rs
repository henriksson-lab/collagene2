//! Source-bound binary translation and output requests; no external executables.
use collagene_core::{
    protocol::{
        self,
        request::Request,
        response::Response,
        schema::{AnalysisJob, AnalysisResult},
    },
    sequence::{
        translation::TranslationOptions,
        translation_reference::{self, Selection},
    },
};
use collagene_model::{Alphabet, Error, Feature, Location, Molecule, Topology};
fn fixture(cds: bool) -> translation_reference::Request {
    let mut source =
        Molecule::dna("original template", "ATGAAACCCGGGTTTTAA", Topology::Linear).unwrap();
    source.features.push(Feature {
        kind: "CDS".into(),
        location: Location::range(0, 18),
        qualifiers: vec![
            ("label".into(), Some("whole protein".into())),
            ("translation".into(), Some("MKPGF".into())),
        ],
        needs_review: false,
    });
    let selection = if cds {
        Selection::Cds { feature_index: 0 }
    } else {
        Selection::Region {
            location: Location::Join(vec![Location::range(0, 6), Location::range(12, 18)]),
            options: TranslationOptions::default(),
        }
    };
    translation_reference::Request { source, selection }
}
fn job(request: translation_reference::Request, output: bool) -> Request {
    let name = "translated protein";
    let revision = if output {
        translation_reference::output_revision(&request, name).unwrap()
    } else {
        translation_reference::inspect_revision(&request).unwrap()
    };
    Request {
        id: 54,
        revision,
        job: if output {
            AnalysisJob::PrepareTranslatedReference {
                request,
                name: name.into(),
            }
        } else {
            AnalysisJob::TranslateReference { request }
        },
    }
}
fn roundtrip(request: &Request) -> (Vec<u8>, Response) {
    let wire =
        protocol::execution::process_request(&protocol::request::Request::encode(request).unwrap())
            .unwrap();
    let response = protocol::response::Response::decode(&wire).unwrap();
    assert_eq!(response.id, request.id);
    assert_eq!(response.revision, request.revision);
    (wire, response)
}
fn tag(wire: &[u8]) -> u32 {
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&wire[42..]).unwrap();
    tag
}
fn response_wire(response: &Response) -> Vec<u8> {
    let body = postcard::to_allocvec(response).unwrap();
    let mut wire = b"CGR\0".to_vec();
    wire.extend_from_slice(&protocol::schema::WORKER_SCHEMA_VERSION.to_le_bytes());
    wire.extend_from_slice(&(body.len() as u32).to_le_bytes());
    wire.extend_from_slice(blake3::hash(&body).as_bytes());
    wire.extend(body);
    wire
}
#[test]
fn compound_region_and_cds_binary_outputs_have_independent_expected_proteins() {
    for cds in [false, true] {
        let input = fixture(cds);
        let original = input.source.clone();
        let expected = if cds { "MKPGF" } else { "MKF*" };
        let request = job(input.clone(), false);
        assert_eq!(
            tag(&protocol::request::Request::encode(&request).unwrap()),
            120
        );
        let (wire, response) = roundtrip(&request);
        assert_eq!(tag(&wire), 112);
        let AnalysisResult::TranslatedReference(report) = response.result else {
            panic!("missing translation report")
        };
        report.validate_binding(&request.revision).unwrap();
        report.validate_sources(&input).unwrap();
        assert_eq!(report.protein, expected);
        assert_eq!(report.source_revision, original.revision());
        assert_eq!(report.cds.is_some(), cds);
        let output = job(input.clone(), true);
        assert_ne!(output.revision, request.revision);
        assert_eq!(
            tag(&protocol::request::Request::encode(&output).unwrap()),
            121
        );
        let (wire, response) = roundtrip(&output);
        assert_eq!(tag(&wire), 113);
        let AnalysisResult::TranslatedOutputs(outputs) = response.result else {
            panic!("missing translated output")
        };
        outputs.validate_binding(&output.revision).unwrap();
        outputs
            .validate_sources(&input, "translated protein")
            .unwrap();
        assert_eq!(outputs.protein.sequence.as_str(), expected);
        assert_eq!(outputs.protein.sequence.alphabet(), Alphabet::Protein);
        assert_eq!(outputs.protein.topology, Topology::Linear);
        assert_eq!(outputs.protein.name, "translated protein");
        assert_eq!(input.source, original);
    }
}
#[test]
fn reverse_compound_location_uses_source_orientation() {
    let mut input = fixture(false);
    if let Selection::Region { location, .. } = &mut input.selection {
        *location = Location::Complement(Box::new(location.clone()));
    }
    let (_, response) = roundtrip(&job(input, false));
    let AnalysisResult::TranslatedReference(report) = response.result else {
        panic!("missing translation report")
    };
    // reverse-complement(ATGAAATTTTAA) = TTAAAATTTCAT.
    assert_eq!(report.protein, "LKFH");
}
#[test]
fn stale_source_options_output_name_and_context_fail_before_computation() {
    for output in [false, true] {
        let request = job(fixture(false), output);
        let mut stale = request.clone();
        stale.revision = "0".repeat(64);
        assert!(
            protocol::execution::process_request(
                &protocol::request::Request::encode(&stale).unwrap()
            )
            .is_err()
        );
        let mut source = request.clone();
        match &mut source.job {
            AnalysisJob::TranslateReference { request }
            | AnalysisJob::PrepareTranslatedReference { request, .. } => {
                request.source.name = "different exact source".into()
            }
            _ => unreachable!(),
        }
        assert!(
            protocol::execution::process_request(
                &protocol::request::Request::encode(&source).unwrap()
            )
            .is_err()
        );
        let mut options = request.clone();
        match &mut options.job {
            AnalysisJob::TranslateReference { request }
            | AnalysisJob::PrepareTranslatedReference { request, .. } => {
                if let Selection::Region { options, .. } = &mut request.selection {
                    options.frame = 1;
                }
            }
            _ => unreachable!(),
        }
        assert!(
            protocol::execution::process_request(
                &protocol::request::Request::encode(&options).unwrap()
            )
            .is_err()
        );
        if output {
            let mut changed = request.clone();
            if let AnalysisJob::PrepareTranslatedReference { name, .. } = &mut changed.job {
                *name = "another protein".into();
            }
            assert!(
                protocol::execution::process_request(
                    &protocol::request::Request::encode(&changed).unwrap()
                )
                .is_err()
            );
        }
    }
}
#[test]
fn schema54_discriminators_reject_schema53_body_prefixes_and_full_envelopes() {
    for (output, job_tag, result_tag) in [(false, 120, 112), (true, 121, 113)] {
        let request = job(fixture(false), output);
        let request_wire = protocol::request::Request::encode(&request).unwrap();
        let (reply_wire, _) = roundtrip(&request);
        for (wire, is_response, expected_tag) in [
            (&request_wire, false, job_tag),
            (&reply_wire, true, result_tag),
        ] {
            assert_eq!(tag(wire), expected_tag);
            let mut old = wire.clone();
            old[4..6].copy_from_slice(&53u16.to_le_bytes());
            let rejected = if is_response {
                protocol::response::Response::decode(&old).map(|_| ())
            } else {
                protocol::execution::process_request(&old).map(|_| ())
            };
            assert!(matches!(rejected, Err(Error::Unsupported(_))));
            // No new variant body exists; the new selector's minimum schema is
            // checked before payload decoding or the deliberately stale checksum.
            let body =
                postcard::to_allocvec(&(request.id, request.revision.as_str(), expected_tag))
                    .unwrap();
            let mut prefix = wire[..42].to_vec();
            prefix[4..6].copy_from_slice(&53u16.to_le_bytes());
            prefix[6..10].copy_from_slice(&(body.len() as u32).to_le_bytes());
            prefix.extend(body);
            let rejected = if is_response {
                protocol::response::Response::decode(&prefix).map(|_| ())
            } else {
                protocol::execution::process_request(&prefix).map(|_| ())
            };
            assert!(matches!(rejected, Err(Error::Unsupported(_))));
        }
    }
}
#[test]
fn translated_response_validation_rejects_wrong_context_and_inconsistent_native_content() {
    for output in [false, true] {
        let (_, response) = roundtrip(&job(fixture(false), output));
        let mut wrong = response.clone();
        wrong.revision = "0".repeat(64);
        assert!(protocol::response::Response::decode(&response_wire(&wrong)).is_err());
        let mut malformed = response.clone();
        match &mut malformed.result {
            AnalysisResult::TranslatedReference(report) => {
                report.translation.translated_codons = usize::MAX
            }
            AnalysisResult::TranslatedOutputs(outputs) => {
                outputs.protein = Molecule::dna("wrong alphabet", "ACG", Topology::Linear).unwrap()
            }
            _ => unreachable!(),
        }
        assert!(protocol::response::Response::decode(&response_wire(&malformed)).is_err());
    }
}
#[test]
fn ordinary_32_mib_caps_apply_to_both_jobs_and_results() {
    let maximum = collagene_core::formats::sequence::MAX_IMPORT_BYTES;
    // One reusable oversized envelope establishes that these selectors cannot
    // claim the separate 128 MiB workspace budget, before checksum/body decode.
    let mut oversized = vec![0; 42 + maximum + 1];
    oversized[4..6].copy_from_slice(&protocol::schema::WORKER_SCHEMA_VERSION.to_le_bytes());
    oversized[6..10].copy_from_slice(&((maximum + 1) as u32).to_le_bytes());
    for (magic, selector, response) in [
        (b"CGQ\0", 120, false),
        (b"CGQ\0", 121, false),
        (b"CGR\0", 112, true),
        (b"CGR\0", 113, true),
    ] {
        oversized[..4].copy_from_slice(magic);
        let body = postcard::to_allocvec(&(54u32, "context", selector as u32)).unwrap();
        oversized[42..42 + body.len()].copy_from_slice(&body);
        let rejected = if response {
            protocol::response::Response::decode(&oversized).map(|_| ())
        } else {
            protocol::execution::process_request(&oversized).map(|_| ())
        };
        assert!(matches!(rejected, Err(Error::TooLarge)));
    }
    drop(oversized);
    let mut input = fixture(false);
    input.source.original_genbank = Some("a".repeat(maximum + 1));
    // A caller-provided envelope revision cannot make oversized original source
    // evidence eligible for the workspace-only serialization exception.
    let request = Request {
        id: 54,
        revision: "0".repeat(64),
        job: AnalysisJob::TranslateReference { request: input },
    };
    assert!(matches!(
        protocol::request::Request::encode(&request),
        Err(Error::TooLarge)
    ));
}

#[test]
fn added_code_profiles_cross_current_binary_boundary_without_collapsing_dual_stops() {
    for (code, expected, possible) in [
        (4, "M**W", 0),
        (27, "MQQX", 1),
        (28, "MXXX", 3),
        (31, "MXXW", 2),
    ] {
        let source = Molecule::dna("literal new table", "ATGTAATAGTGA", Topology::Linear).unwrap();
        let input = translation_reference::Request {
            source: source.clone(),
            selection: Selection::Region {
                location: Location::range(0, 12),
                options: TranslationOptions {
                    genetic_code: code,
                    ..Default::default()
                },
            },
        };
        let request = job(input.clone(), false);
        let (_, response) = roundtrip(&request);
        let AnalysisResult::TranslatedReference(report) = response.result else {
            panic!("wrong result")
        };
        assert_eq!(report.protein, expected);
        assert_eq!(report.translation.possible_stop_codons, possible);
        report.validate_sources(&input).unwrap();
        assert_eq!(report.source_revision, source.revision());
    }
}
