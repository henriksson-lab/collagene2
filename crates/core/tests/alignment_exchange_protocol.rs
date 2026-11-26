use collagene_core::{
    Alphabet, Molecule, Topology,
    alignment_formats::{self as alignment, AlignmentFormat},
    protocol::{self, AnalysisJob, AnalysisResult, Request},
    workspace::{self, Workspace, WorkspaceAttachment},
};

const NEXUS: &[u8] = b"#NEXUS\nBEGIN DATA;\nDIMENSIONS NTAX=2 NCHAR=6;\nFORMAT DATATYPE=DNA GAP=-;\nMATRIX\n'alpha one' AC-GT-\nbeta A-TGT-\n;\nEND;\n";

fn version(mut bytes: Vec<u8>, value: u16) -> Vec<u8> {
    bytes[4..6].copy_from_slice(&value.to_le_bytes());
    bytes
}
fn request(job: AnalysisJob, revision: String) -> Vec<u8> {
    protocol::encode_request(&Request {
        id: 91,
        revision,
        job,
    })
    .unwrap()
}
fn document() -> alignment::AlignmentDocument {
    alignment::import_alignment(NEXUS, AlignmentFormat::Nexus, Alphabet::Dna, "source").unwrap()
}

#[test]
fn new_profiles_require_worker87_but_legacy_alignment_frames_keep_their_floor() {
    let bytes = request(
        AnalysisJob::ImportAlignment {
            bytes: NEXUS.to_vec(),
            format: AlignmentFormat::Nexus,
            alphabet: Alphabet::Dna,
            name: "source".into(),
        },
        blake3::hash(NEXUS).to_hex().to_string(),
    );
    assert!(protocol::process_request(&version(bytes.clone(), 86)).is_err());
    let result = protocol::process_request(&bytes).unwrap();
    assert!(protocol::decode_response(&version(result.clone(), 86)).is_err());
    let AnalysisResult::AlignmentDocument(document) =
        protocol::decode_response(&result).unwrap().result
    else {
        panic!("expected imported alignment")
    };
    assert_eq!(document.rows[0].identifier, "alpha one");
    assert_eq!(document.rows[0].aligned, "AC-GT-");
    let old = b">alpha\nAC-GT-\n>beta\nA-TGT-\n";
    let legacy = request(
        AnalysisJob::ImportAlignment {
            bytes: old.to_vec(),
            format: AlignmentFormat::AlignedFasta,
            alphabet: Alphabet::Dna,
            name: "legacy".into(),
        },
        blake3::hash(old).to_hex().to_string(),
    );
    let response = protocol::process_request(&version(legacy, 10)).unwrap();
    assert!(matches!(
        protocol::decode_response(&response).unwrap().result,
        AnalysisResult::AlignmentDocument(_)
    ));
}

#[test]
fn source_profiles_and_edited_space_labels_cannot_hide_in_old_save_export_or_extract_frames() {
    let document = document();
    let revision = protocol::alignment_revision(&document).unwrap();
    for job in [
        AnalysisJob::ExportAlignment {
            document: document.clone(),
            format: AlignmentFormat::Nexus,
        },
        AnalysisJob::SaveAlignment {
            document: document.clone(),
        },
        AnalysisJob::ExtractAlignmentRows {
            document: document.clone(),
            row_ids: vec![document.rows[0].id.clone()],
        },
    ] {
        let bytes = request(job, revision.clone());
        assert!(protocol::process_request(&version(bytes.clone(), 86)).is_err());
        let response = protocol::process_request(&bytes).unwrap();
        assert!(!matches!(
            protocol::decode_response(&response).unwrap().result,
            AnalysisResult::Error(_)
        ));
    }
    let mut edited = alignment::import_alignment(
        b">alpha\nACGT\n",
        AlignmentFormat::AlignedFasta,
        Alphabet::Dna,
        "legacy",
    )
    .unwrap();
    edited.rows[0].identifier = "alpha one".into();
    edited.validate().unwrap();
    assert!(!edited.source_is_current().unwrap());
    let bytes = request(
        AnalysisJob::SaveAlignment {
            document: edited.clone(),
        },
        protocol::alignment_revision(&edited).unwrap(),
    );
    assert!(protocol::process_request(&version(bytes.clone(), 86)).is_err());
    assert!(matches!(
        protocol::decode_response(&protocol::process_request(&bytes).unwrap())
            .unwrap()
            .result,
        AnalysisResult::AlignmentSaved(_)
    ));
    assert!(
        alignment::export_alignment(&edited, AlignmentFormat::AlignedFasta).is_err(),
        "unquoted FASTA cannot silently change a spaced identifier"
    );
}

#[test]
fn native_alignment_and_workspace_attachments_preserve_originals_and_enforce_profile_floor() {
    let document = document();
    let archive = alignment::save_alignment(&document).unwrap();
    let bytes = request(
        AnalysisJob::LoadAlignment {
            bytes: archive.clone(),
        },
        blake3::hash(&archive).to_hex().to_string(),
    );
    assert!(protocol::process_request(&version(bytes.clone(), 86)).is_err());
    let AnalysisResult::AlignmentDocument(reloaded) =
        protocol::decode_response(&protocol::process_request(&bytes).unwrap())
            .unwrap()
            .result
    else {
        panic!("expected native alignment")
    };
    assert_eq!(reloaded.original_bytes().unwrap(), NEXUS);
    let mut frame = Workspace::from_documents(
        "exchange",
        vec![Molecule::dna("plasmid", "ACGT", Topology::Linear).unwrap()],
        0,
    )
    .unwrap();
    frame.attachments.insert(
        blake3::hash(&archive).to_hex().to_string(),
        WorkspaceAttachment {
            kind: alignment::WORKSPACE_ATTACHMENT_KIND.into(),
            name: "alignment.cg2aln".into(),
            media_type: "application/json".into(),
            bytes: archive,
        },
    );
    let archive = workspace::save_workspace(&frame).unwrap();
    let bytes = request(
        AnalysisJob::LoadWorkspace {
            bytes: archive.clone(),
        },
        blake3::hash(&archive).to_hex().to_string(),
    );
    assert!(protocol::process_request(&version(bytes.clone(), 86)).is_err());
    assert!(matches!(
        protocol::decode_response(&protocol::process_request(&bytes).unwrap())
            .unwrap()
            .result,
        AnalysisResult::WorkspaceLoaded(_)
    ));
}
