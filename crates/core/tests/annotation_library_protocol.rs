//! Revision-bound, versioned worker contracts for user feature and annotation evidence.
use collagene_core::{
    Feature, Location, Molecule, Topology,
    annotation_formats::{self, AnnotationFormat},
    feature_library::{self, FeatureLibrary, ScanOptions},
    protocol::{self, AnalysisJob, AnalysisResult, Request, Response},
};
fn source() -> Molecule {
    let mut molecule = Molecule::dna("source part", "CGATAC", Topology::Linear).unwrap();
    molecule.features.push(Feature {
        kind: "misc_feature".into(),
        location: Location::range(0, 6),
        qualifiers: vec![("label".into(), Some("Own source part".into()))],
        needs_review: false,
    });
    molecule
}
fn target() -> Molecule {
    Molecule::dna("independent target", "ATACG", Topology::Circular).unwrap()
}
fn library() -> FeatureLibrary {
    let mut result = FeatureLibrary::new("Own library", "1").unwrap();
    result.add_feature(&source(), 0).unwrap();
    result
}
fn linear_target() -> Molecule {
    Molecule::dna("reference", "TTTCGATACGGG", Topology::Linear).unwrap()
}
fn execute(molecule: &Molecule, job: AnalysisJob) -> Response {
    let request = Request {
        id: 823,
        revision: molecule.revision(),
        job,
    };
    let bytes = protocol::encode_request(&request).unwrap();
    assert_eq!(&bytes[4..6], &21u16.to_le_bytes());
    let response = protocol::decode_response(&protocol::process_request(&bytes).unwrap()).unwrap();
    assert_eq!(response.id, request.id);
    assert_eq!(response.revision, request.revision);
    response
}
#[test]
fn worker_scan_and_selected_apply_reconstruct_actual_source_bound_annotation() {
    let molecule = linear_target();
    let library = library();
    let AnalysisResult::FeatureLibraryScan(report) = execute(
        &molecule,
        AnalysisJob::ScanFeatureLibrary {
            molecule: molecule.clone(),
            library: library.clone(),
            options: ScanOptions::default(),
        },
    )
    .result
    else {
        panic!("not a scan")
    };
    assert_eq!(report.source_revision, molecule.revision());
    assert_eq!(report.library_revision, library.revision());
    assert!(report.complete);
    assert_eq!(report.hits.len(), 1);
    let hit = &report.hits[0];
    assert_eq!((hit.start, hit.length, hit.reverse), (3, 6, false));
    assert_eq!(
        (
            hit.concrete_matches,
            hit.compatible_ambiguous,
            hit.mismatches
        ),
        (6, 0, 0)
    );
    assert_eq!(
        (hit.identity_lower, hit.identity_upper, hit.coverage),
        (1., 1., 1.)
    );
    let selected = vec![hit.id.clone()];
    let AnalysisResult::AnnotatedDocument(annotated) = execute(
        &molecule,
        AnalysisJob::ApplyFeatureLibrary {
            molecule: molecule.clone(),
            library,
            report,
            selected,
        },
    )
    .result
    else {
        panic!("not applied")
    };
    assert_eq!(annotated.sequence, molecule.sequence);
    assert_eq!(annotated.features.len(), 1);
    assert_eq!(annotated.features[0].location.spans(), vec![(3, 9, false)]);
    assert!(annotated.features[0].needs_review);
    assert!(molecule.features.is_empty());
}
#[test]
fn selected_apply_rejects_tampered_report_and_old_reference_without_partial_edits() {
    let molecule = linear_target();
    let library = library();
    let mut report = feature_library::scan(&molecule, &library, &ScanOptions::default()).unwrap();
    let id = report.hits[0].id.clone();
    report.hits[0]
        .feature
        .qualifiers
        .push(("label".into(), Some("Forged label".into())));
    assert!(matches!(
        execute(
            &molecule,
            AnalysisJob::ApplyFeatureLibrary {
                molecule: molecule.clone(),
                library: library.clone(),
                report,
                selected: vec![id]
            }
        )
        .result,
        AnalysisResult::Error(_)
    ));
    let report = feature_library::scan(&molecule, &library, &ScanOptions::default()).unwrap();
    let changed = molecule.replace(0, 1, "A").unwrap();
    let selected = vec![report.hits[0].id.clone()];
    assert!(matches!(
        execute(
            &changed,
            AnalysisJob::ApplyFeatureLibrary {
                molecule: changed.clone(),
                library,
                report,
                selected
            }
        )
        .result,
        AnalysisResult::Error(_)
    ));
    assert!(molecule.features.is_empty());
}
#[test]
fn explicit_gff_mapping_import_apply_export_preserves_reverse_coordinates_and_raw_source() {
    let molecule = linear_target();
    let text =
        "##gff-version 3\nchr1\tfixture\tmisc_feature\t3\t6\t.\t-\t.\tID=f1;Name=source%3Blabel\n";
    let AnalysisResult::AnnotationImport(report) = execute(
        &molecule,
        AnalysisJob::ImportAnnotations {
            molecule: molecule.clone(),
            text: text.into(),
            format: AnnotationFormat::Gff3,
            seqid: "chr1".into(),
        },
    )
    .result
    else {
        panic!("not imported")
    };
    assert_eq!(report.reference_revision, molecule.revision());
    assert_eq!(report.original_bytes().unwrap(), text.as_bytes());
    assert_eq!(report.proposals.len(), 1);
    assert_eq!(
        report.proposals[0].feature.location.spans(),
        vec![(2, 6, true)]
    );
    let AnalysisResult::AnnotatedDocument(annotated) = execute(
        &molecule,
        AnalysisJob::ApplyAnnotations {
            molecule: molecule.clone(),
            report,
            selected: vec![0],
        },
    )
    .result
    else {
        panic!("not applied")
    };
    assert_eq!(annotated.sequence, molecule.sequence);
    assert_eq!(annotated.features.len(), 1);
    let AnalysisResult::AnnotationExport(export) = execute(
        &annotated,
        AnalysisJob::ExportAnnotations {
            molecule: annotated.clone(),
            format: AnnotationFormat::Gff3,
            seqid: "chr1".into(),
        },
    )
    .result
    else {
        panic!("not exported")
    };
    let reopened = annotation_formats::import_annotations(
        &molecule,
        &export.text,
        AnnotationFormat::Gff3,
        "chr1",
    )
    .unwrap();
    assert_eq!(
        reopened.proposals[0].feature.location,
        annotated.features[0].location
    );
}
#[test]
fn new_job_and_result_tags_require_seven_while_version_six_search_still_decodes() {
    let molecule = linear_target();
    let library = library();
    let scan = feature_library::scan(&molecule, &library, &ScanOptions::default()).unwrap();
    let annotation = annotation_formats::import_annotations(
        &molecule,
        "chr1\t0\t3\tname\t0\t+\n",
        AnnotationFormat::Bed,
        "chr1",
    )
    .unwrap();
    let jobs = vec![
        AnalysisJob::ScanFeatureLibrary {
            molecule: molecule.clone(),
            library: library.clone(),
            options: ScanOptions::default(),
        },
        AnalysisJob::ApplyFeatureLibrary {
            molecule: molecule.clone(),
            library,
            selected: vec![scan.hits[0].id.clone()],
            report: scan,
        },
        AnalysisJob::ImportAnnotations {
            molecule: molecule.clone(),
            text: "chr1\t0\t3\n".into(),
            format: AnnotationFormat::Bed,
            seqid: "chr1".into(),
        },
        AnalysisJob::ApplyAnnotations {
            molecule: molecule.clone(),
            report: annotation,
            selected: vec![0],
        },
        AnalysisJob::ExportAnnotations {
            molecule: molecule.clone(),
            format: AnnotationFormat::Bed,
            seqid: "chr1".into(),
        },
    ];
    for (i, job) in jobs.into_iter().enumerate() {
        let request = Request {
            id: 9,
            revision: molecule.revision(),
            job,
        };
        let bytes = protocol::encode_request(&request).unwrap();
        let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[42..]).unwrap();
        assert_eq!(tag, 22 + i as u32);
        let mut reply = protocol::process_request(&bytes).unwrap();
        reply[4..6].copy_from_slice(&6u16.to_le_bytes());
        assert!(protocol::decode_response(&reply).is_err());
        let mut old = bytes.clone();
        old[4..6].copy_from_slice(&6u16.to_le_bytes());
        assert!(protocol::process_request(&old).is_err());
        let mut stale = request;
        stale.revision = source().revision();
        assert!(protocol::process_request(&protocol::encode_request(&stale).unwrap()).is_err());
    }
    let mut legacy = protocol::encode_request(&Request {
        id: 12,
        revision: molecule.revision(),
        job: AnalysisJob::Search {
            molecule: molecule.clone(),
            query: "CGATAC".into(),
        },
    })
    .unwrap();
    legacy[4..6].copy_from_slice(&6u16.to_le_bytes());
    let response = protocol::decode_response(&protocol::process_request(&legacy).unwrap()).unwrap();
    assert!(matches!(response.result, AnalysisResult::Search(_)));
}
#[test]
fn large_new_job_header_is_rejected_before_its_payload_is_deserialized() {
    // A deliberately incomplete new-job payload: preflight must reject its declared
    // normal-message oversize independently of scientific/model deserialization.
    let mut payload = postcard::to_allocvec(&(9u32, "revision", 22u32)).unwrap();
    payload.resize(32 * 1024 * 1024 + 1, 0);
    let mut envelope = b"CGQ\0".to_vec();
    envelope.extend(7u16.to_le_bytes());
    envelope.extend((payload.len() as u32).to_le_bytes());
    envelope.extend(blake3::hash(&payload).as_bytes());
    envelope.extend(payload);
    assert!(matches!(
        protocol::process_request(&envelope),
        Err(collagene_core::Error::TooLarge)
    ));
    // A circular target shorter than the library is legal, and never creates a
    // second traversal of the template when no one-turn placement is possible.
    assert!(
        feature_library::scan(&target(), &library(), &ScanOptions::default())
            .unwrap()
            .hits
            .is_empty()
    );
}

#[test]
fn annotated_worker_result_cannot_smuggle_invalid_feature_coordinates() {
    let mut molecule = linear_target();
    molecule.features.push(Feature {
        kind: "misc_feature".into(),
        location: Location::range(0, 100),
        qualifiers: vec![],
        needs_review: true,
    });
    let response = Response {
        id: 1,
        revision: molecule.revision(),
        result: AnalysisResult::AnnotatedDocument(molecule),
    };
    let payload = postcard::to_allocvec(&response).unwrap();
    let mut envelope = b"CGR\0".to_vec();
    envelope.extend(7u16.to_le_bytes());
    envelope.extend((payload.len() as u32).to_le_bytes());
    envelope.extend(blake3::hash(&payload).as_bytes());
    envelope.extend(payload);
    assert!(protocol::decode_response(&envelope).is_err());
}
