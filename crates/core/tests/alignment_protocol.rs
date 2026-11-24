//! Independently specified exchange wire bindings and legacy-version boundaries.
use collagene_core::{
    Alphabet,
    alignment_formats::{self as formats, AlignmentFormat},
    protocol::*,
};
const FASTA: &[u8] = b">alpha first\nA-CG\n>beta\nAT-G\n";
fn document() -> formats::AlignmentDocument {
    formats::import_alignment(FASTA, AlignmentFormat::AlignedFasta, Alphabet::Dna, "pairs").unwrap()
}
fn request(revision: String, job: AnalysisJob) -> Vec<u8> {
    encode_request(&Request {
        id: 91,
        revision,
        job,
    })
    .unwrap()
}
fn roundtrip(revision: String, job: AnalysisJob) -> Response {
    decode_response(&process_request(&request(revision, job)).unwrap()).unwrap()
}
fn encode_response(response: &Response) -> Result<Vec<u8>, collagene_core::Error> {
    let payload = postcard::to_allocvec(response).unwrap();
    let mut bytes = b"CGR\0".to_vec();
    bytes.extend(10u16.to_le_bytes());
    bytes.extend((payload.len() as u32).to_le_bytes());
    bytes.extend(blake3::hash(&payload).as_bytes());
    bytes.extend(payload);
    Ok(bytes)
}
#[test]
fn exact_import_and_native_load_have_distinct_input_evidence_bindings() {
    let hash = blake3::hash(FASTA).to_hex().to_string();
    let response = roundtrip(
        hash.clone(),
        AnalysisJob::ImportAlignment {
            bytes: FASTA.to_vec(),
            format: AlignmentFormat::AlignedFasta,
            alphabet: Alphabet::Dna,
            name: "pairs".into(),
        },
    );
    assert_eq!(response.id, 91);
    assert_eq!(response.revision, hash);
    let AnalysisResult::AlignmentDocument(d) = response.result else {
        panic!()
    };
    assert_eq!(d, document());
    assert_eq!(d.original_bytes().unwrap(), FASTA);
    let bytes = formats::save_alignment(&d).unwrap();
    let hash = blake3::hash(&bytes).to_hex().to_string();
    let response = roundtrip(hash.clone(), AnalysisJob::LoadAlignment { bytes });
    assert_eq!(response.revision, hash);
    let AnalysisResult::AlignmentDocument(loaded) = response.result else {
        panic!()
    };
    assert_eq!(loaded, d);
}
#[test]
fn row_extraction_has_literal_sequence_orientation_and_no_fabricated_annotation() {
    let d = document();
    let rev = alignment_revision(&d).unwrap();
    let rows = d.rows.iter().map(|r| r.id.clone()).collect();
    let response = roundtrip(
        rev.clone(),
        AnalysisJob::ExtractAlignmentRows {
            document: d,
            row_ids: rows,
        },
    );
    assert_eq!(response.revision, rev);
    let AnalysisResult::AlignmentRows(rows) = response.result else {
        panic!()
    };
    assert_eq!(rows[0].sequence.as_str(), "ACG");
    assert_eq!(rows[1].sequence.as_str(), "ATG");
    assert!(rows.iter().all(|r| r.features.is_empty()));
    assert_eq!(rows[0].metadata["alignment_current_blake3"], rev);
}
#[test]
fn export_and_save_bind_full_content_not_only_original_source() {
    let d = document();
    let rev = alignment_revision(&d).unwrap();
    let response = roundtrip(
        rev.clone(),
        AnalysisJob::ExportAlignment {
            document: d.clone(),
            format: AlignmentFormat::Stockholm,
        },
    );
    let AnalysisResult::AlignmentExport(export) = response.result else {
        panic!()
    };
    assert!(export.text.starts_with("# STOCKHOLM 1.0\n"));
    assert!(!export.losses.is_empty());
    let response = roundtrip(
        rev.clone(),
        AnalysisJob::SaveAlignment {
            document: d.clone(),
        },
    );
    let AnalysisResult::AlignmentSaved(bytes) = response.result else {
        panic!()
    };
    assert_eq!(formats::load_alignment(&bytes).unwrap(), d);
    let mut edited = d;
    edited.name = "renamed".into();
    assert_ne!(alignment_revision(&edited).unwrap(), rev);
    assert!(
        process_request(&request(
            rev,
            AnalysisJob::SaveAlignment { document: edited }
        ))
        .is_err()
    );
}
#[test]
fn incoming_byte_hash_mismatch_and_corrupt_inner_json_are_explicit() {
    assert!(
        process_request(&request(
            "wrong".into(),
            AnalysisJob::ImportAlignment {
                bytes: FASTA.to_vec(),
                format: AlignmentFormat::AlignedFasta,
                alphabet: Alphabet::Dna,
                name: "x".into()
            }
        ))
        .is_err()
    );
    let bytes = b"{}".to_vec();
    let hash = blake3::hash(&bytes).to_hex().to_string();
    let response = roundtrip(hash, AnalysisJob::LoadAlignment { bytes });
    assert!(matches!(response.result, AnalysisResult::Error(_)));
}
#[test]
fn every_new_request_and_result_rejects_pre_ten_labels_without_changing_old_tags() {
    let d = document();
    let rev = alignment_revision(&d).unwrap();
    let bytes = formats::save_alignment(&d).unwrap();
    let jobs = vec![
        AnalysisJob::ImportAlignment {
            bytes: FASTA.to_vec(),
            format: AlignmentFormat::AlignedFasta,
            alphabet: Alphabet::Dna,
            name: "pairs".into(),
        },
        AnalysisJob::ExportAlignment {
            document: d.clone(),
            format: AlignmentFormat::AlignedFasta,
        },
        AnalysisJob::ExtractAlignmentRows {
            document: d.clone(),
            row_ids: vec![],
        },
        AnalysisJob::SaveAlignment {
            document: d.clone(),
        },
        AnalysisJob::LoadAlignment {
            bytes: bytes.clone(),
        },
    ];
    for (index, job) in jobs.into_iter().enumerate() {
        let original = request(rev.clone(), job);
        let ((id, _, tag), _) =
            postcard::take_from_bytes::<(u32, &str, u32)>(&original[42..]).unwrap();
        assert_eq!((id, tag), (91, 34 + index as u32));
        assert_eq!(&original[4..6], &23u16.to_le_bytes());
        for old in 1..10u16 {
            let mut message = original.clone();
            message[4..6].copy_from_slice(&old.to_le_bytes());
            assert!(process_request(&message).is_err(), "tag{tag} old{old}");
        }
    }
    let results = vec![
        AnalysisResult::AlignmentDocument(d),
        AnalysisResult::AlignmentExport(formats::AlignmentExport {
            text: "x".into(),
            losses: vec![],
            renamings: Default::default(),
        }),
        AnalysisResult::AlignmentRows(vec![]),
        AnalysisResult::AlignmentSaved(bytes),
    ];
    for (index, result) in results.into_iter().enumerate() {
        let original = encode_response(&Response {
            id: 91,
            revision: rev.clone(),
            result,
        })
        .unwrap();
        let ((_, _, tag), _) =
            postcard::take_from_bytes::<(u32, &str, u32)>(&original[42..]).unwrap();
        assert_eq!(tag, 33 + index as u32);
        for old in 1..10u16 {
            let mut message = original.clone();
            message[4..6].copy_from_slice(&old.to_le_bytes());
            assert!(decode_response(&message).is_err(), "result{tag} old{old}");
        }
    }
}
#[test]
fn source_integrity_and_extraction_errors_cannot_be_hidden_in_typed_response() {
    let mut d = document();
    d.source.original_base64 = "not base64".into();
    let bytes = encode_response(&Response {
        id: 91,
        revision: "x".into(),
        result: AnalysisResult::AlignmentDocument(d),
    })
    .unwrap();
    assert!(decode_response(&bytes).is_err());
    let d = document();
    let rev = alignment_revision(&d).unwrap();
    let response = roundtrip(
        rev,
        AnalysisJob::ExtractAlignmentRows {
            document: d,
            row_ids: vec!["unknown".into()],
        },
    );
    assert!(matches!(response.result, AnalysisResult::Error(_)));
}
#[test]
fn protein_property_job_and_report_are_bound_and_cannot_claim_older_schema() {
    use collagene_core::{Molecule, Sequence, Topology, protein_properties::*};
    let mut molecule = Molecule::dna("protein", "ACW", Topology::Linear).unwrap();
    molecule.sequence = Sequence::new(Alphabet::Protein, "ACW").unwrap();
    let revision = molecule.revision();
    let encoded = request(
        revision.clone(),
        AnalysisJob::ProteinProperties {
            molecule: molecule.clone(),
            options: ProteinOptions::default(),
        },
    );
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&encoded[42..]).unwrap();
    assert_eq!(tag, 39);
    let response = decode_response(&process_request(&encoded).unwrap()).unwrap();
    let AnalysisResult::ProteinProperties(report) = response.result else {
        panic!()
    };
    assert_eq!(report.source_revision, revision);
    assert_eq!(report.residue_count, 3);
    for old in 1..10u16 {
        let mut bytes = encoded.clone();
        bytes[4..6].copy_from_slice(&old.to_le_bytes());
        assert!(process_request(&bytes).is_err());
    }
    let response = Response {
        id: 91,
        revision: revision.clone(),
        result: AnalysisResult::ProteinProperties(report),
    };
    let encoded = encode_response(&response).unwrap();
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&encoded[42..]).unwrap();
    assert_eq!(tag, 37);
    for old in 1..10u16 {
        let mut bytes = encoded.clone();
        bytes[4..6].copy_from_slice(&old.to_le_bytes());
        assert!(decode_response(&bytes).is_err());
    }
    let wrong = Response {
        revision: "stale".into(),
        ..response
    };
    assert!(decode_response(&encode_response(&wrong).unwrap()).is_err());
}
#[test]
fn saved_alignment_response_hash_must_match_exact_canonical_document_revision() {
    let d = document();
    let bytes = formats::save_alignment(&d).unwrap();
    let revision = alignment_revision(&d).unwrap();
    let response = Response {
        id: 91,
        revision,
        result: AnalysisResult::AlignmentSaved(bytes.clone()),
    };
    assert!(decode_response(&encode_response(&response).unwrap()).is_ok());
    let wrong = Response {
        revision: "other source".into(),
        ..response
    };
    assert!(decode_response(&encode_response(&wrong).unwrap()).is_err());
    // A second valid alignment cannot be substituted while keeping the old revision.
    let mut changed = d;
    changed.name = "different valid alignment".into();
    let substituted = Response {
        id: 91,
        revision: blake3::hash(&bytes).to_hex().to_string(),
        result: AnalysisResult::AlignmentSaved(formats::save_alignment(&changed).unwrap()),
    };
    assert!(decode_response(&encode_response(&substituted).unwrap()).is_err());
}
