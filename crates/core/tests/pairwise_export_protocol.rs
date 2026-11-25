//! Portable pairwise reports use the ordinary binary worker envelope, not a
//! workspace exception or an external formatter.
use collagene_core::{
    Error, Feature, Location, Molecule, Topology,
    alignment::AlignmentOptions,
    pairwise_export::{self, Format},
    pairwise_history::Record,
    pairwise_review::{self, Options},
    pairwise_search,
    protocol::{self, AnalysisJob, AnalysisResult, Request, Response},
    run_lineage::EventBinding,
    workspace::Workspace,
};

fn record() -> Record {
    let reference = Molecule::dna("reference", "ACGTACGT", Topology::Linear).unwrap();
    let query = Molecule::dna("query", "ACCTACGT", Topology::Linear).unwrap();
    let frame =
        Workspace::from_documents("pairwise export", vec![reference.clone(), query.clone()], 0)
            .unwrap();
    let binding = |index: usize| {
        let doc = &frame.documents[index];
        EventBinding {
            document_id: doc.id.clone(),
            event_id: doc.head_event.clone().unwrap(),
            revision: doc.head.clone(),
        }
    };
    Record {
        schema_version: 1,
        reference: binding(0),
        query: binding(1),
        report: pairwise_review::compare(
            &reference,
            &query,
            &Options::Stored(AlignmentOptions::default()),
        )
        .unwrap(),
    }
}
fn request(format: Format) -> Request {
    let record = record();
    let reference_name = "Reference <retained>";
    let query_name = "Query & original";
    Request {
        id: 56,
        revision: pairwise_export::revision(&record, reference_name, query_name, false, format)
            .unwrap(),
        job: AnalysisJob::ExportPairwiseComparison {
            record,
            reference_name: reference_name.into(),
            query_name: query_name.into(),
            fresh: false,
            format,
        },
    }
}
fn selector(wire: &[u8]) -> u32 {
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&wire[42..]).unwrap();
    tag
}
fn roundtrip(request: &Request) -> (Vec<u8>, Response) {
    let wire = protocol::encode_request(request).unwrap();
    assert_eq!(selector(&wire), 125);
    let reply = protocol::process_request(&wire).unwrap();
    assert_eq!(selector(&reply), 117);
    let response = protocol::decode_response(&reply).unwrap();
    assert_eq!(response.id, request.id);
    assert_eq!(response.revision, request.revision);
    (reply, response)
}
fn response_wire(response: &Response) -> Vec<u8> {
    let body = postcard::to_allocvec(response).unwrap();
    let mut wire = b"CGR\0".to_vec();
    wire.extend_from_slice(&protocol::WORKER_SCHEMA_VERSION.to_le_bytes());
    wire.extend_from_slice(&(body.len() as u32).to_le_bytes());
    wire.extend_from_slice(blake3::hash(&body).as_bytes());
    wire.extend(body);
    wire
}

#[test]
fn all_formats_execute_inside_the_binary_worker_and_preserve_declared_rows() {
    for (format, extension, mime) in [
        (Format::AlignedFasta, "fasta", "text/plain;charset=utf-8"),
        (Format::Clustal, "aln", "text/plain;charset=utf-8"),
        (Format::Html, "html", "text/html;charset=utf-8"),
        (Format::Pdf, "pdf", "application/pdf"),
    ] {
        let (_, response) = roundtrip(&request(format));
        let AnalysisResult::ExportedPairwiseComparison(export) = response.result else {
            panic!("missing exported comparison")
        };
        export.validate().unwrap();
        assert_eq!(export.extension, extension);
        assert_eq!(export.mime, mime);
        if extension == "pdf" {
            assert!(export.bytes.starts_with(b"%PDF-"));
        } else {
            let text = String::from_utf8(export.bytes).unwrap();
            assert!(text.contains("ACGTACGT"));
            assert!(text.contains("ACCTACGT"));
            if extension == "html" {
                assert!(text.contains("&lt;retained&gt;"));
                assert!(!text.contains("Reference <retained>"));
            }
        }
    }
}

#[test]
fn format_names_freshness_and_exact_events_are_bound_before_export() {
    let original = request(Format::Html);
    for index in 0..5 {
        let mut changed = original.clone();
        if let AnalysisJob::ExportPairwiseComparison {
            reference_name,
            query_name,
            format,
            fresh,
            record,
        } = &mut changed.job
        {
            match index {
                0 => reference_name.push_str(" changed"),
                1 => query_name.push_str(" changed"),
                2 => *format = Format::Clustal,
                3 => *fresh = true,
                _ => record.reference.document_id.push_str("-alias"),
            }
        }
        let wire = protocol::encode_request(&changed).unwrap();
        assert!(protocol::process_request(&wire).is_err());
    }
    let mut wrong = original;
    wrong.revision = "0".repeat(64);
    assert!(protocol::process_request(&protocol::encode_request(&wrong).unwrap()).is_err());
}

#[test]
fn schema55_cannot_carry_new_selectors_even_when_the_variant_body_is_missing() {
    let request = request(Format::AlignedFasta);
    let original = protocol::encode_request(&request).unwrap();
    let (reply, _) = roundtrip(&request);
    for (wire, response, tag) in [(&original, false, 125u32), (&reply, true, 117u32)] {
        let mut old = wire.clone();
        old[4..6].copy_from_slice(&55u16.to_le_bytes());
        let rejected = if response {
            protocol::decode_response(&old).map(|_| ())
        } else {
            protocol::process_request(&old).map(|_| ())
        };
        assert!(matches!(rejected, Err(Error::Unsupported(_))));
        let body = postcard::to_allocvec(&(request.id, request.revision.as_str(), tag)).unwrap();
        let mut prefix = old[..42].to_vec();
        prefix[6..10].copy_from_slice(&(body.len() as u32).to_le_bytes());
        prefix.extend(body);
        let rejected = if response {
            protocol::decode_response(&prefix).map(|_| ())
        } else {
            protocol::process_request(&prefix).map(|_| ())
        };
        assert!(matches!(rejected, Err(Error::Unsupported(_))));
    }
}

#[test]
fn old_comparison_layout_remains_readable_under_schema55() {
    let record = record();
    let reference = Molecule::dna("reference", "ACGTACGT", Topology::Linear).unwrap();
    let query = Molecule::dna("query", "ACCTACGT", Topology::Linear).unwrap();
    let request = Request {
        id: 55,
        revision: record.report.context_revision,
        job: AnalysisJob::CompareSequences {
            reference,
            query,
            options: Options::Stored(AlignmentOptions::default()),
        },
    };
    let mut wire = protocol::encode_request(&request).unwrap();
    assert_eq!(selector(&wire), 119);
    wire[4..6].copy_from_slice(&55u16.to_le_bytes());
    let mut reply = protocol::process_request(&wire).unwrap();
    assert_eq!(selector(&reply), 111);
    reply[4..6].copy_from_slice(&55u16.to_le_bytes());
    assert!(matches!(
        protocol::decode_response(&reply).unwrap().result,
        AnalysisResult::ReviewedComparison(_)
    ));
}

#[test]
fn export_responses_reject_unsafe_metadata_and_keep_the_ordinary_32mib_budget() {
    let (_, response) = roundtrip(&request(Format::AlignedFasta));
    for index in 0..3 {
        let mut changed = response.clone();
        let AnalysisResult::ExportedPairwiseComparison(export) = &mut changed.result else {
            unreachable!()
        };
        match index {
            0 => export.mime = "text/javascript".into(),
            1 => export.extension = "../../script.js".into(),
            _ => export.warnings = vec!["a".repeat(65_537)],
        }
        assert!(protocol::decode_response(&response_wire(&changed)).is_err());
    }
    let mut malformed = request(Format::Html);
    if let AnalysisJob::ExportPairwiseComparison { reference_name, .. } = &mut malformed.job {
        *reference_name = "x".repeat(65_537);
    }
    assert!(protocol::encode_request(&malformed).is_err());

    let maximum = collagene_core::formats::MAX_IMPORT_BYTES;
    let mut wire = vec![0; 42 + maximum + 1];
    wire[4..6].copy_from_slice(&protocol::WORKER_SCHEMA_VERSION.to_le_bytes());
    wire[6..10].copy_from_slice(&((maximum + 1) as u32).to_le_bytes());
    for (magic, tag, response) in [
        (b"CGQ\0", 125u32, false),
        (b"CGR\0", 117u32, true),
        (b"CGQ\0", 126u32, false),
        (b"CGR\0", 118u32, true),
        (b"CGQ\0", 127u32, false),
        (b"CGR\0", 119u32, true),
    ] {
        wire[..4].copy_from_slice(magic);
        let body = postcard::to_allocvec(&(56u32, "export", tag)).unwrap();
        wire[42..42 + body.len()].copy_from_slice(&body);
        let rejected = if response {
            protocol::decode_response(&wire).map(|_| ())
        } else {
            protocol::process_request(&wire).map(|_| ())
        };
        assert!(matches!(rejected, Err(Error::TooLarge)));
    }
}

#[test]
fn aligned_row_search_runs_in_the_worker_with_native_coordinates_and_context_guards() {
    let report = record().report;
    let pattern = "ACCT";
    let options = pairwise_search::Options::default();
    let revision = pairwise_search::revision(&report, pattern, &options).unwrap();
    let request = Request {
        id: 57,
        revision,
        job: AnalysisJob::SearchPairwiseComparison {
            report,
            pattern: pattern.into(),
            options,
        },
    };
    let wire = protocol::encode_request(&request).unwrap();
    assert_eq!(selector(&wire), 126);
    let reply = protocol::process_request(&wire).unwrap();
    assert_eq!(selector(&reply), 118);
    let response = protocol::decode_response(&reply).unwrap();
    assert_eq!(response.revision, request.revision);
    let mut malformed = response.clone();
    let AnalysisResult::SearchedPairwiseComparison(matches) = &mut malformed.result else {
        panic!("missing searched comparison")
    };
    matches.hits[0].columns.1 = usize::MAX;
    assert!(protocol::decode_response(&response_wire(&malformed)).is_err());
    let AnalysisResult::SearchedPairwiseComparison(matches) = response.result else {
        panic!("missing searched comparison")
    };
    assert!(matches.complete);
    assert_eq!(matches.pattern, "ACCT");
    assert_eq!(matches.hits.len(), 1);
    let hit = &matches.hits[0];
    assert_eq!(hit.row, pairwise_search::Row::Query);
    assert_eq!(hit.columns, (0, 4));
    assert_eq!(hit.oriented, (0, 4));
    assert_eq!(hit.source_spans, vec![(0, 4)]);
    let mut stale = request.clone();
    if let AnalysisJob::SearchPairwiseComparison { pattern, .. } = &mut stale.job {
        *pattern = "ACGT".into();
    }
    assert!(protocol::process_request(&protocol::encode_request(&stale).unwrap()).is_err());
    for (mut old, response) in [(wire, false), (reply, true)] {
        old[4..6].copy_from_slice(&55u16.to_le_bytes());
        let rejected = if response {
            protocol::decode_response(&old).map(|_| ())
        } else {
            protocol::process_request(&old).map(|_| ())
        };
        assert!(matches!(rejected, Err(Error::Unsupported(_))));
    }
}

#[test]
fn annotation_context_worker_preserves_gap_columns_and_rejects_stale_or_old_envelopes() {
    let mut reference = Molecule::dna("reference", "AACCGGTT", Topology::Linear).unwrap();
    reference.features.push(Feature {
        kind: "misc_feature".into(),
        location: Location::range(2, 6),
        qualifiers: vec![("label".into(), Some("four source calls".into()))],
        needs_review: false,
    });
    let query = Molecule::dna("query", "AACCTGGTT", Topology::Linear).unwrap();
    let report = pairwise_review::compare(
        &reference,
        &query,
        &Options::Stored(AlignmentOptions::default()),
    )
    .unwrap();
    assert_eq!(report.comparison.alignment.reference_aligned, "AACC-GGTT");
    let revision =
        collagene_core::pairwise_context::revision(&report, &reference, &query, false).unwrap();
    let request = Request {
        id: 58,
        revision,
        job: AnalysisJob::PreparePairwiseContext {
            report: report.clone(),
            reference,
            query,
            translate: false,
        },
    };
    let wire = protocol::encode_request(&request).unwrap();
    assert_eq!(selector(&wire), 127);
    let reply = protocol::process_request(&wire).unwrap();
    assert_eq!(selector(&reply), 119);
    let response = protocol::decode_response(&reply).unwrap();
    assert_eq!(response.revision, request.revision);
    let AnalysisResult::PreparedPairwiseContext(context) = &response.result else {
        panic!("missing annotation context")
    };
    assert_eq!(context.context_revision, request.revision);
    let page = context.page(&report, 0).unwrap();
    assert_eq!(page.reference, "AACC-GGTT");
    assert_eq!(page.query, "AACCTGGTT");
    assert_eq!(
        page.reference_positions,
        vec![
            Some(0),
            Some(1),
            Some(2),
            Some(3),
            None,
            Some(4),
            Some(5),
            Some(6),
            Some(7)
        ]
    );
    assert_eq!(page.tracks.len(), 1);
    assert_eq!(page.tracks[0].annotation, "  >> >>  ");
    assert!(page.tracks[0].translation.is_none());
    let mut changed = request.clone();
    if let AnalysisJob::PreparePairwiseContext { translate, .. } = &mut changed.job {
        *translate = true;
    }
    assert!(protocol::process_request(&protocol::encode_request(&changed).unwrap()).is_err());
    let mut changed = response.clone();
    changed.revision = "0".repeat(64);
    assert!(protocol::decode_response(&response_wire(&changed)).is_err());
    for (mut old, response, tag) in [(wire, false, 127u32), (reply, true, 119u32)] {
        old[4..6].copy_from_slice(&55u16.to_le_bytes());
        let rejected = if response {
            protocol::decode_response(&old).map(|_| ())
        } else {
            protocol::process_request(&old).map(|_| ())
        };
        assert!(matches!(rejected, Err(Error::Unsupported(_))));
        let body = postcard::to_allocvec(&(request.id, request.revision.as_str(), tag)).unwrap();
        old.truncate(42);
        old[6..10].copy_from_slice(&(body.len() as u32).to_le_bytes());
        old.extend(body);
        let rejected = if response {
            protocol::decode_response(&old).map(|_| ())
        } else {
            protocol::process_request(&old).map(|_| ())
        };
        assert!(matches!(rejected, Err(Error::Unsupported(_))));
    }
}
