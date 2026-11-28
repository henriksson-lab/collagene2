//! Revision-bound, versioned worker contracts for user feature and annotation evidence.
use collagene_core::{
    formats::annotations::{self as annotation_formats, AnnotationFormat},
    protocol::{
        self,
        request::Request,
        response::Response,
        schema::{AnalysisJob, AnalysisResult},
    },
    sequence::annotations::{
        library::{self as feature_library, FeatureLibrary, ScanOptions},
        merge::{self, Decision, Selection, Source},
    },
};
use collagene_model::{Feature, Location, Molecule, Topology};
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
fn merge_request(molecule: &Molecule, source: Source, ids: Vec<String>) -> merge::Request {
    merge::Request {
        schema_version: 1,
        target_revision: molecule.revision(),
        source,
        selections: ids
            .into_iter()
            .map(|proposal_id| Selection {
                proposal_id,
                decision: Decision::Append,
            })
            .collect(),
        limits: Default::default(),
    }
}
fn execute(molecule: &Molecule, job: AnalysisJob) -> Response {
    let request = Request {
        id: 823,
        revision: match &job {
            AnalysisJob::PreviewAnnotationMerge { molecule, request }
            | AnalysisJob::ApplyAnnotationMerge { molecule, request } => {
                merge::revision(molecule, request).unwrap()
            }
            _ => molecule.revision(),
        },
        job,
    };
    let bytes = protocol::request::Request::encode(&request).unwrap();
    assert_eq!(
        &bytes[4..6],
        &collagene_core::protocol::schema::WORKER_SCHEMA_VERSION.to_le_bytes()
    );
    let response = protocol::response::Response::decode(
        &protocol::execution::process_request(&bytes).unwrap(),
    )
    .unwrap();
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
    let AnalysisResult::AnnotationMergeApplied(applied) = execute(
        &molecule,
        AnalysisJob::ApplyAnnotationMerge {
            molecule: molecule.clone(),
            request: merge_request(&molecule, Source::Library { library, report }, selected),
        },
    )
    .result
    else {
        panic!("not applied")
    };
    let annotated = applied.molecule;
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
            AnalysisJob::ApplyAnnotationMerge {
                molecule: molecule.clone(),
                request: merge_request(
                    &molecule,
                    Source::Library {
                        library: library.clone(),
                        report
                    },
                    vec![id]
                ),
            }
        )
        .result,
        AnalysisResult::Error(_)
    ));
    let report = feature_library::scan(&molecule, &library, &ScanOptions::default()).unwrap();
    let changed = molecule.replace(0, 1, "A").unwrap();
    let selected = vec![report.hits[0].id.clone()];
    let original_context = merge::revision(
        &molecule,
        &merge_request(
            &molecule,
            Source::Library {
                library: library.clone(),
                report: report.clone(),
            },
            selected.clone(),
        ),
    )
    .unwrap();
    let stale = Request {
        id: 824,
        revision: original_context,
        job: AnalysisJob::ApplyAnnotationMerge {
            molecule: changed.clone(),
            request: merge_request(&changed, Source::Library { library, report }, selected),
        },
    };
    assert!(
        stale.encode().is_err(),
        "stale source pins refuse outbound publication"
    );
    assert!(changed.features.is_empty());
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
    let AnalysisResult::AnnotationMergeApplied(applied) = execute(
        &molecule,
        AnalysisJob::ApplyAnnotationMerge {
            molecule: molecule.clone(),
            request: merge_request(
                &molecule,
                Source::Dataset {
                    report: report.clone(),
                },
                vec![report.proposals[0].id.clone()],
            ),
        },
    )
    .result
    else {
        panic!("not applied")
    };
    let annotated = applied.molecule;
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
fn annotation_jobs_keep_current_selectors_and_reject_unsupported_versions() {
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
        AnalysisJob::PreviewAnnotationMerge {
            molecule: molecule.clone(),
            request: merge_request(
                &molecule,
                Source::Library {
                    library,
                    report: scan.clone(),
                },
                vec![scan.hits[0].id.clone()],
            ),
        },
        AnalysisJob::ImportAnnotations {
            molecule: molecule.clone(),
            text: "chr1\t0\t3\n".into(),
            format: AnnotationFormat::Bed,
            seqid: "chr1".into(),
        },
        AnalysisJob::ApplyAnnotationMerge {
            molecule: molecule.clone(),
            request: merge_request(
                &molecule,
                Source::Dataset {
                    report: annotation.clone(),
                },
                vec![annotation.proposals[0].id.clone()],
            ),
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
            revision: match &job {
                AnalysisJob::PreviewAnnotationMerge { molecule, request }
                | AnalysisJob::ApplyAnnotationMerge { molecule, request } => {
                    merge::revision(molecule, request).unwrap()
                }
                _ => molecule.revision(),
            },
            job,
        };
        let bytes = protocol::request::Request::encode(&request).unwrap();
        let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[42..]).unwrap();
        assert_eq!(tag, 22 + i as u32);
        let mut reply = protocol::execution::process_request(&bytes).unwrap();
        reply[4..6].copy_from_slice(&6u16.to_le_bytes());
        assert!(protocol::response::Response::decode(&reply).is_err());
        let mut old = bytes.clone();
        old[4..6].copy_from_slice(&6u16.to_le_bytes());
        assert!(protocol::execution::process_request(&old).is_err());
        let mut stale = request;
        stale.revision = source().revision();
        if matches!(
            &stale.job,
            AnalysisJob::PreviewAnnotationMerge { .. } | AnalysisJob::ApplyAnnotationMerge { .. }
        ) {
            assert!(stale.encode().is_err());
        } else {
            assert!(protocol::execution::process_request(&stale.encode().unwrap()).is_err());
        }
    }
}

#[test]
fn large_new_job_header_is_rejected_before_its_payload_is_deserialized() {
    // A deliberately incomplete new-job payload: preflight must reject its declared
    // normal-message oversize independently of scientific/model deserialization.
    let mut payload = postcard::to_allocvec(&(9u32, "revision", 22u32)).unwrap();
    payload.resize(32 * 1024 * 1024 + 1, 0);
    let mut envelope = b"CGQ\0".to_vec();
    envelope.extend(collagene_core::protocol::schema::WORKER_SCHEMA_VERSION.to_le_bytes());
    envelope.extend((payload.len() as u32).to_le_bytes());
    envelope.extend(blake3::hash(&payload).as_bytes());
    envelope.extend(payload);
    assert!(matches!(
        protocol::execution::process_request(&envelope),
        Err(collagene_model::Error::TooLarge)
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
    let molecule = linear_target();
    let library = library();
    let report = feature_library::scan(&molecule, &library, &ScanOptions::default()).unwrap();
    let request = merge_request(
        &molecule,
        Source::Library {
            library,
            report: report.clone(),
        },
        vec![report.hits[0].id.clone()],
    );
    let mut applied = merge::apply(&molecule, &request).unwrap();
    applied.molecule.features[0].location = Location::range(0, 100);
    let response = Response {
        id: 1,
        revision: merge::revision(&molecule, &request).unwrap(),
        result: AnalysisResult::AnnotationMergeApplied(applied),
    };
    let payload = postcard::to_allocvec(&response).unwrap();
    let mut envelope = b"CGR\0".to_vec();
    envelope.extend(collagene_core::protocol::schema::WORKER_SCHEMA_VERSION.to_le_bytes());
    envelope.extend((payload.len() as u32).to_le_bytes());
    envelope.extend(blake3::hash(&payload).as_bytes());
    envelope.extend(payload);
    assert!(protocol::response::Response::decode(&envelope).is_err());
}

#[test]
fn merge_context_replies_and_large_prefixes_keep_exact_current_binary_association() {
    use collagene_core::protocol::schema::WORKER_SCHEMA_VERSION;
    let molecule = linear_target();
    let library = library();
    let report = feature_library::scan(&molecule, &library, &ScanOptions::default()).unwrap();
    let request = merge_request(
        &molecule,
        Source::Library {
            library,
            report: report.clone(),
        },
        vec![report.hits[0].id.clone()],
    );
    let context = merge::revision(&molecule, &request).unwrap();
    let frame = |magic: &[u8; 4], payload: Vec<u8>| {
        let mut bytes = magic.to_vec();
        bytes.extend_from_slice(&WORKER_SCHEMA_VERSION.to_le_bytes());
        bytes.extend_from_slice(&(payload.len() as u32).to_le_bytes());
        bytes.extend_from_slice(blake3::hash(&payload).as_bytes());
        bytes.extend(payload);
        bytes
    };
    for (job, reply_selector) in [
        (
            AnalysisJob::PreviewAnnotationMerge {
                molecule: molecule.clone(),
                request: request.clone(),
            },
            205u32,
        ),
        (
            AnalysisJob::ApplyAnnotationMerge {
                molecule: molecule.clone(),
                request: request.clone(),
            },
            24u32,
        ),
    ] {
        let current = Request {
            id: 901,
            revision: context.clone(),
            job,
        };
        let bytes = current.encode().unwrap();
        let mut old = bytes.clone();
        old[4..6].copy_from_slice(&139u16.to_le_bytes());
        assert!(protocol::execution::process_request(&old).is_err());
        let output = protocol::execution::process_request(&bytes).unwrap();
        let ((_, _, tag), _) =
            postcard::take_from_bytes::<(u32, &str, u32)>(&output[42..]).unwrap();
        assert_eq!(tag, reply_selector);
        let mut reply = Response::decode(&output).unwrap();
        let mut old = output;
        old[4..6].copy_from_slice(&139u16.to_le_bytes());
        assert!(Response::decode(&old).is_err());
        match &mut reply.result {
            AnalysisResult::AnnotationMergePreview(preview) => {
                preview.report_revision = "0".repeat(64)
            }
            AnalysisResult::AnnotationMergeApplied(applied) => {
                applied.receipt.product_revision = "0".repeat(64)
            }
            _ => panic!("unexpected merge response"),
        }
        assert!(
            Response::decode(&frame(b"CGR\0", postcard::to_allocvec(&reply).unwrap())).is_err()
        );
        let mut stale = current;
        stale.revision = "0".repeat(64);
        assert!(stale.encode().is_err());
        assert!(
            protocol::execution::process_request(&frame(
                b"CGQ\0",
                postcard::to_allocvec(&stale).unwrap()
            ))
            .is_err()
        );
    }
    let payload_length = 32 * 1024 * 1024 + 1;
    for selector in [24u32, 205] {
        let mut prefix = b"CGR\0".to_vec();
        prefix.extend_from_slice(&WORKER_SCHEMA_VERSION.to_le_bytes());
        prefix.extend_from_slice(&(payload_length as u32).to_le_bytes());
        prefix.extend_from_slice(&[0; 32]);
        prefix.extend(postcard::to_allocvec(&(901u32, context.as_str(), selector)).unwrap());
        assert!(prefix.len() < 256);
        assert_eq!(
            Response::admitted_frame_limit(&prefix, 42 + payload_length).unwrap(),
            128 * 1024 * 1024
        );
        assert!(Response::decode(&prefix).is_err());
    }
}
