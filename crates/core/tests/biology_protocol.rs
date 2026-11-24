//! Worker-boundary tests: independently specified biology survives binary dispatch,
//! preserves request identity, and reports computation errors without transport ambiguity.
use collagene_core::{
    alignment::{AlignmentMode, AlignmentOptions, DifferenceKind},
    protocol::{
        AnalysisJob, AnalysisResult, Request, Response, decode_response, encode_request,
        process_request,
    },
    translation::{InitiationPolicy, OrfOptions, TranslationOptions},
    *,
};
use serde::Serialize;

fn molecule(bases: &str) -> Molecule {
    Molecule::dna("worker fixture", bases, Topology::Linear).unwrap()
}
fn request(source: &Molecule, job: AnalysisJob) -> Request {
    Request {
        id: 0xfeed_beef,
        revision: source.revision(),
        job,
    }
}
fn run(request: &Request) -> Response {
    let bytes = encode_request(request).unwrap();
    assert_eq!(&bytes[..4], b"CGQ\0");
    assert_eq!(u16::from_le_bytes([bytes[4], bytes[5]]), 21);
    let response_bytes = process_request(&bytes).unwrap();
    assert_eq!(&response_bytes[..4], b"CGR\0");
    assert_eq!(
        u16::from_le_bytes([response_bytes[4], response_bytes[5]]),
        21
    );
    let response = decode_response(&response_bytes).unwrap();
    assert_eq!(response.id, request.id);
    assert_eq!(response.revision, request.revision);
    response
}

#[test]
fn translation_dispatch_preserves_genetic_code_frame_orientation_and_initiation() {
    // Reverse complement is AGTGGTGTAA; skip A and recode only initiating GTG.
    let source = molecule("TTACACCACT");
    let options = TranslationOptions {
        genetic_code: 11,
        frame: 1,
        reverse: true,
        initiation: InitiationPolicy::RequireRecognizedStart,
    };
    let response = run(&request(
        &source,
        AnalysisJob::Translate {
            molecule: source.clone(),
            options: options.clone(),
        },
    ));
    let AnalysisResult::Translation(report) = response.result else {
        panic!("wrong result variant")
    };
    assert_eq!(report.options, options);
    assert_eq!(report.protein, "MV*");
    assert!(report.initiator_methionine_applied);
    assert_eq!(report.trailing_bases, 0);
}

#[test]
fn cds_dispatch_extracts_compound_origin_feature_and_binds_annotation_revision() {
    let mut source = molecule("AAATAACCCATG");
    source.topology = Topology::Circular;
    let location = Location::Join(vec![Location::range(9, 12), Location::range(0, 6)]);
    source.features.push(Feature {
        kind: "CDS".into(),
        location: location.clone(),
        needs_review: false,
        qualifiers: vec![
            ("label".into(), Some("origin CDS".into())),
            ("transl_table".into(), Some("11".into())),
            ("codon_start".into(), Some("1".into())),
            ("translation".into(), Some("MK".into())),
        ],
    });
    let response = run(&request(
        &source,
        AnalysisJob::TranslateCds {
            molecule: source.clone(),
            feature_index: 0,
        },
    ));
    let AnalysisResult::CdsTranslation(report) = response.result else {
        panic!("wrong result variant")
    };
    assert_eq!(report.source_revision, source.revision());
    assert_eq!(report.location, location);
    assert_eq!(report.protein, "MK");
    assert_eq!(report.translation.options.genetic_code, 11);
    assert_eq!(report.codon_start, 1);
    assert_eq!(report.stored_translation_matches, Some(true));
    assert!(report.has_terminal_stop && report.complete_under_model);

    let mut changed = source.clone();
    changed.features[0].qualifiers[2].1 = Some("2".into());
    let stale = Request {
        id: 1,
        revision: source.revision(),
        job: AnalysisJob::TranslateCds {
            molecule: changed,
            feature_index: 0,
        },
    };
    assert!(
        matches!(process_request(&encode_request(&stale).unwrap()), Err(Error::Format(message))
        if message.contains("revision mismatch"))
    );
}

#[test]
fn orf_dispatch_preserves_circular_locations_and_explicit_limits() {
    let mut source = molecule("AAATAACCCATG");
    source.topology = Topology::Circular;
    let options = OrfOptions {
        min_amino_acids: 1,
        max_orfs: 7,
        both_strands: false,
        ..OrfOptions::default()
    };
    let response = run(&request(
        &source,
        AnalysisJob::FindOrfs {
            molecule: source.clone(),
            options: options.clone(),
        },
    ));
    let AnalysisResult::Orfs(report) = response.result else {
        panic!("wrong result variant")
    };
    assert_eq!(report.source_revision, source.revision());
    assert_eq!(report.options, options);
    assert!(!report.truncated);
    assert_eq!(report.orfs.len(), 1);
    let orf = &report.orfs[0];
    assert_eq!(orf.protein, "MK");
    assert_eq!(
        orf.location,
        Location::Join(vec![Location::range(9, 12), Location::range(0, 6)])
    );
    assert!(orf.crosses_origin && orf.has_start_codon && orf.has_stop_codon);
    assert!(!orf.reverse && !orf.partial_5prime && !orf.partial_3prime);
}

#[test]
fn alignment_dispatch_preserves_mode_scoring_and_difference_coordinates() {
    let source = molecule("GGACGTCC");
    let query = molecule("ACAT");
    let mut options = AlignmentOptions {
        mode: AlignmentMode::QueryInReference,
        ..Default::default()
    };
    options.scoring.match_score = 3;
    let response = run(&request(
        &source,
        AnalysisJob::Align {
            molecule: source.clone(),
            query: query.clone(),
            options,
        },
    ));
    let AnalysisResult::Alignment(report) = response.result else {
        panic!("wrong result variant")
    };
    assert_eq!(report.options, options);
    assert_eq!(report.score, 6); // Three matches at +3, one mismatch at -3.
    assert_eq!((report.reference_start, report.reference_end), (2, 6));
    assert_eq!((report.query_start, report.query_end), (0, 4));
    assert_eq!(report.reference_aligned, "ACGT");
    assert_eq!(report.query_aligned, "ACAT");
    assert_eq!(report.cigar, "2=1X1=");
    assert_eq!(report.differences.len(), 1);
    let difference = &report.differences[0];
    assert_eq!(difference.kind, DifferenceKind::Substitution);
    assert_eq!((difference.reference_start, difference.query_start), (4, 2));
    assert_eq!(
        (
            difference.reference_bases.as_str(),
            difference.query_bases.as_str()
        ),
        ("G", "A")
    );
}

#[test]
fn every_new_job_rejects_a_stale_source_before_analysis() {
    let source = molecule("ATGAAATAA");
    let jobs = [
        AnalysisJob::Translate {
            molecule: source.clone(),
            options: Default::default(),
        },
        AnalysisJob::TranslateCds {
            molecule: source.clone(),
            feature_index: 0,
        },
        AnalysisJob::FindOrfs {
            molecule: source.clone(),
            options: Default::default(),
        },
        AnalysisJob::Align {
            molecule: source.clone(),
            query: molecule("ATG"),
            options: Default::default(),
        },
    ];
    for job in jobs {
        let stale = Request {
            id: 17,
            revision: molecule("ATGCCCTAA").revision(),
            job,
        };
        assert!(
            matches!(process_request(&encode_request(&stale).unwrap()), Err(Error::Format(message))
            if message.contains("revision mismatch"))
        );
    }
}

#[test]
fn computation_failures_return_typed_error_responses_with_request_identity() {
    let source = molecule("ATGAAATAA");
    let mut invalid_query = molecule("ATG");
    invalid_query.features.push(Feature {
        kind: "misc_feature".into(),
        location: Location::range(0, 99),
        qualifiers: vec![],
        needs_review: false,
    });
    let jobs = [
        AnalysisJob::Translate {
            molecule: source.clone(),
            options: TranslationOptions {
                frame: 3,
                ..Default::default()
            },
        },
        AnalysisJob::TranslateCds {
            molecule: source.clone(),
            feature_index: usize::MAX,
        },
        AnalysisJob::FindOrfs {
            molecule: source.clone(),
            options: OrfOptions {
                max_orfs: 0,
                ..Default::default()
            },
        },
        AnalysisJob::Align {
            molecule: source.clone(),
            query: invalid_query,
            options: Default::default(),
        },
    ];
    for job in jobs {
        let response = run(&request(&source, job));
        let AnalysisResult::Error(message) = response.result else {
            panic!("expected analysis error")
        };
        assert!(!message.is_empty());
    }
}

// These legacy wire types deliberately do not use today's AnalysisJob/Result:
// appended enum variants must not change the discriminants of original variants.
#[derive(Serialize)]
enum LegacyJob {
    Search { molecule: Molecule, query: String },
}
#[derive(Serialize)]
struct LegacyRequest {
    id: u32,
    revision: String,
    job: LegacyJob,
}
#[derive(Serialize)]
enum LegacyResult {
    Search(Vec<collagene_core::analysis::Match>),
}
#[derive(Serialize)]
struct LegacyResponse {
    id: u32,
    revision: String,
    result: LegacyResult,
}

fn legacy_envelope(magic: &[u8; 4], version: u16, value: &impl Serialize) -> Vec<u8> {
    let payload = postcard::to_allocvec(value).unwrap();
    let mut bytes = magic.to_vec();
    bytes.extend_from_slice(&version.to_le_bytes());
    bytes.extend_from_slice(&(payload.len() as u32).to_le_bytes());
    bytes.extend_from_slice(blake3::hash(&payload).as_bytes());
    bytes.extend_from_slice(&payload);
    bytes
}

#[test]
fn original_worker_variants_remain_compatible_with_versions_one_and_two() {
    let source = molecule("AAACCC");
    for version in [1, 2] {
        let legacy = LegacyRequest {
            id: 72,
            revision: source.revision(),
            job: LegacyJob::Search {
                molecule: source.clone(),
                query: "AAA".into(),
            },
        };
        let reply = process_request(&legacy_envelope(b"CGQ\0", version, &legacy)).unwrap();
        let response = decode_response(&reply).unwrap();
        assert_eq!(
            (response.id, response.revision.as_str()),
            (72, source.revision().as_str())
        );
        let AnalysisResult::Search(matches) = response.result else {
            panic!("legacy variant changed")
        };
        assert_eq!(matches.len(), 1);
        assert_eq!((matches[0].start, matches[0].length), (0, 3));
        let old_reply = legacy_envelope(
            b"CGR\0",
            version,
            &LegacyResponse {
                id: 72,
                revision: source.revision(),
                result: LegacyResult::Search(matches),
            },
        );
        let decoded = decode_response(&old_reply).unwrap();
        assert_eq!(decoded.id, 72);
        assert!(matches!(decoded.result, AnalysisResult::Search(found) if found.len() == 1));
    }
}

#[test]
fn unsupported_worker_versions_and_corrupt_payloads_are_transport_errors() {
    let source = molecule("ATGAAATAA");
    let valid = encode_request(&request(
        &source,
        AnalysisJob::Translate {
            molecule: source.clone(),
            options: Default::default(),
        },
    ))
    .unwrap();
    let valid_response = process_request(&valid).unwrap();
    for version in [0u16, 22, u16::MAX] {
        let mut bytes = valid.clone();
        bytes[4..6].copy_from_slice(&version.to_le_bytes());
        assert!(matches!(
            process_request(&bytes),
            Err(Error::Unsupported(_))
        ));
        let mut response = valid_response.clone();
        response[4..6].copy_from_slice(&version.to_le_bytes());
        assert!(matches!(
            decode_response(&response),
            Err(Error::Unsupported(_))
        ));
    }
    let mut corrupt = valid;
    *corrupt.last_mut().unwrap() ^= 1;
    assert!(
        matches!(process_request(&corrupt), Err(Error::Format(message))
        if message.contains("checksum"))
    );
    let mut corrupt_response = valid_response;
    *corrupt_response.last_mut().unwrap() ^= 1;
    assert!(
        matches!(decode_response(&corrupt_response), Err(Error::Format(message))
        if message.contains("checksum"))
    );
}
