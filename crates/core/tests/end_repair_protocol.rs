use collagene_core::{
    chemistry::{AnnealOptions, EndRepairMethod},
    protocol::{self, AnalysisJob, AnalysisResult, Request, Response},
    transfer,
    verification::*,
    *,
};
fn source() -> Molecule {
    Molecule::dna("restriction source", "GAATTCAAAA", Topology::Circular).unwrap()
}
fn job(m: &Molecule, job: AnalysisJob) -> Vec<u8> {
    protocol::encode_request(&Request {
        id: 73,
        revision: m.revision(),
        job,
    })
    .unwrap()
}
fn response(bytes: &[u8]) -> Response {
    protocol::decode_response(&protocol::process_request(bytes).unwrap()).unwrap()
}
fn envelope(magic: &[u8; 4], version: u16, value: &impl serde::Serialize) -> Vec<u8> {
    let payload = postcard::to_allocvec(value).unwrap();
    let mut result = magic.to_vec();
    result.extend(version.to_le_bytes());
    result.extend((payload.len() as u32).to_le_bytes());
    result.extend(blake3::hash(&payload).as_bytes());
    result.extend(payload);
    result
}
#[test]
fn new_repair_and_anneal_jobs_bind_sources_and_preserve_explicit_evidence() {
    let m = source();
    let bytes = job(
        &m,
        AnalysisJob::RepairPreview {
            molecule: m.clone(),
            enzymes: vec!["EcoRI".into()],
            fragment_index: 0,
            method: EndRepairMethod::Blunt,
        },
    );
    assert_eq!(&bytes[4..6], &23u16.to_le_bytes());
    let result = response(&bytes);
    assert_eq!(result.id, 73);
    assert_eq!(result.revision, m.revision());
    let AnalysisResult::RepairPreview(preview) = result.result else {
        panic!("repair")
    };
    assert_eq!(preview.molecule.sequence.as_str(), "AATTCAAAAGAATT");
    assert_eq!(preview.source_revision, m.revision());
    assert_eq!(preview.repair.duplex.top.len(), 14);
    let top = Molecule::dna("top", "ACGTCG", Topology::Linear).unwrap();
    let bottom = Molecule::dna("bottom", "CGACGT", Topology::Linear).unwrap();
    let bytes = job(
        &top,
        AnalysisJob::AnnealPreview {
            molecule: top.clone(),
            bottom: bottom.clone(),
            options: AnnealOptions {
                min_paired_bases: 4,
                ..Default::default()
            },
        },
    );
    let result = response(&bytes);
    let AnalysisResult::AnnealPreview(preview) = result.result else {
        panic!("anneal")
    };
    assert_eq!(preview.top_revision, top.revision());
    assert_eq!(preview.bottom_revision, bottom.revision());
    assert!(preview.annealing.unique);
}
#[test]
fn appended_tags_and_new_results_cannot_claim_old_worker_versions() {
    let m = source();
    let jobs = vec![
        AnalysisJob::RepairPreview {
            molecule: m.clone(),
            enzymes: vec!["EcoRI".into()],
            fragment_index: 0,
            method: EndRepairMethod::Blunt,
        },
        AnalysisJob::AnnealPreview {
            molecule: Molecule::dna("top", "ACGTCG", Topology::Linear).unwrap(),
            bottom: Molecule::dna("bottom", "CGACGT", Topology::Linear).unwrap(),
            options: AnnealOptions {
                min_paired_bases: 4,
                ..Default::default()
            },
        },
    ];
    for (index, job_kind) in jobs.into_iter().enumerate() {
        let request = Request {
            id: 77,
            revision: match &job_kind {
                AnalysisJob::RepairPreview { molecule, .. }
                | AnalysisJob::AnnealPreview { molecule, .. } => molecule.revision(),
                _ => unreachable!(),
            },
            job: job_kind,
        };
        let encoded = protocol::encode_request(&request).unwrap();
        let ((_, _, tag), _) =
            postcard::take_from_bytes::<(u32, &str, u32)>(&encoded[42..]).unwrap();
        assert_eq!(tag, 27 + index as u32);
        let result = response(&encoded);
        for version in 1..8 {
            assert!(protocol::process_request(&envelope(b"CGQ\0", version, &request)).is_err());
            assert!(protocol::decode_response(&envelope(b"CGR\0", version, &result)).is_err());
        }
    }
}
#[test]
fn legacy_version_seven_jobs_still_work_and_stale_new_job_bindings_are_rejected() {
    let m = source();
    let request = Request {
        id: 41,
        revision: m.revision(),
        job: AnalysisJob::Search {
            molecule: m.clone(),
            query: "GAA".into(),
        },
    };
    assert!(matches!(
        response(&envelope(b"CGQ\0", 7, &request)).result,
        AnalysisResult::Search(_)
    ));
    let request = Request {
        id: 9,
        revision: "0".repeat(64),
        job: AnalysisJob::RepairPreview {
            molecule: m,
            enzymes: vec!["EcoRI".into()],
            fragment_index: 0,
            method: EndRepairMethod::Blunt,
        },
    };
    assert!(protocol::process_request(&protocol::encode_request(&request).unwrap()).is_err());
}
fn repair_plan() -> (DesignPlan, Molecule) {
    let source = source();
    let target = Molecule::dna("independent", "AATTCAAAAGAATT", Topology::Linear).unwrap();
    (
        DesignPlan {
            name: "repair".into(),
            target_revision: target.revision(),
            inputs: vec![PlanInput {
                id: "source".into(),
                revision: source.revision(),
                molecule: source,
            }],
            steps: vec![
                Operation::Digest {
                    template: Source::Input("source".into()),
                    enzymes: vec!["EcoRI".into()],
                    fragment_index: 0,
                },
                Operation::RepairEnds {
                    template: Source::Step(0),
                    method: EndRepairMethod::Blunt,
                },
            ],
        },
        target,
    )
}
#[test]
fn candidates_and_nested_plan_jobs_reject_new_operations_under_old_schema_labels() {
    let (plan, target) = repair_plan();
    let candidate = transfer::Candidate::new(plan.clone());
    let bytes = transfer::save_candidate(&candidate).unwrap();
    assert_eq!(&bytes[4..6], &7u16.to_le_bytes());
    assert_eq!(
        transfer::load_candidate(&bytes).unwrap().plan_revision,
        candidate.plan_revision
    );
    for version in 1..4 {
        assert!(transfer::load_candidate(&envelope(b"CGP\0", version, &candidate)).is_err());
    }
    for job in [
        AnalysisJob::Verify {
            plan: plan.clone(),
            target: target.clone(),
        },
        AnalysisJob::ExportProtocolPdf {
            plan,
            target: target.clone(),
        },
    ] {
        let request = Request {
            id: 92,
            revision: target.revision(),
            job,
        };
        for version in 1..8 {
            assert!(protocol::process_request(&envelope(b"CGQ\0", version, &request)).is_err());
        }
        assert!(protocol::process_request(&protocol::encode_request(&request).unwrap()).is_ok());
    }
}
#[test]
fn oversized_new_normal_jobs_are_rejected_before_scientific_deserialization() {
    let mut payload = postcard::to_allocvec(&(1u32, "revision", 27u32)).unwrap();
    payload.resize(32 * 1024 * 1024 + 1, 0);
    let mut bytes = b"CGQ\0".to_vec();
    bytes.extend(8u16.to_le_bytes());
    bytes.extend((payload.len() as u32).to_le_bytes());
    bytes.extend(blake3::hash(&payload).as_bytes());
    bytes.extend(payload);
    assert!(matches!(
        protocol::process_request(&bytes),
        Err(Error::TooLarge)
    ));
}
#[test]
fn malformed_repair_result_geometry_never_becomes_valid_preview_state() {
    let m = source();
    let bytes = job(
        &m,
        AnalysisJob::RepairPreview {
            molecule: m.clone(),
            enzymes: vec!["EcoRI".into()],
            fragment_index: 0,
            method: EndRepairMethod::Blunt,
        },
    );
    let mut output = response(&bytes);
    let AnalysisResult::RepairPreview(preview) = &mut output.result else {
        panic!("repair")
    };
    preview.repair.duplex.bottom_start = 10_000;
    assert!(protocol::decode_response(&envelope(b"CGR\0", 8, &output)).is_err());
}

#[test]
fn annealed_document_replays_report_in_worker_and_never_imports_physical_authority() {
    let top = Molecule::dna("top", "ACGTCG", Topology::Linear).unwrap();
    let bottom = Molecule::dna("bottom", "CGACGT", Topology::Linear).unwrap();
    let preview = simulate_oligo_annealing(
        &top,
        &bottom,
        &AnnealOptions {
            min_paired_bases: 4,
            ..Default::default()
        },
    )
    .unwrap();
    let request = Request {
        id: 12,
        revision: top.revision(),
        job: AnalysisJob::AnnealedDocument {
            molecule: top.clone(),
            bottom: bottom.clone(),
            preview: preview.clone(),
            candidate_index: 0,
        },
    };
    let encoded = protocol::encode_request(&request).unwrap();
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&encoded[42..]).unwrap();
    assert_eq!(tag, 29);
    let result = response(&encoded);
    let AnalysisResult::PhysicalPreviewDocument(molecule) = &result.result else {
        panic!("physical preview")
    };
    assert_eq!(molecule.sequence.as_str(), "ACGTCG");
    assert!(collagene_core::chemistry::Duplex::from_molecule(molecule).is_err());
    assert!(protocol::process_request(&envelope(b"CGQ\0", 7, &request)).is_err());
    assert!(protocol::decode_response(&envelope(b"CGR\0", 7, &result)).is_err());
    let mut stale = preview;
    stale.bottom_revision = "0".repeat(64);
    let bytes = job(
        &top,
        AnalysisJob::AnnealedDocument {
            molecule: top.clone(),
            bottom,
            preview: stale,
            candidate_index: 0,
        },
    );
    assert!(matches!(response(&bytes).result, AnalysisResult::Error(_)));
    let mut invalid = molecule.clone();
    invalid.features.push(Feature {
        kind: "misc_feature".into(),
        location: Location::range(0, 100),
        qualifiers: vec![],
        needs_review: true,
    });
    assert!(
        protocol::decode_response(&envelope(
            b"CGR\0",
            8,
            &Response {
                id: 12,
                revision: top.revision(),
                result: AnalysisResult::PhysicalPreviewDocument(invalid)
            }
        ))
        .is_err()
    );
}
