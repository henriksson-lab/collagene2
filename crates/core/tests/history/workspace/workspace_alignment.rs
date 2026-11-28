//! Recognized alignment attachments are validated before any archive restoration.
use collagene_core::{
    alignment_formats::{self, AlignmentFormat, WORKSPACE_ATTACHMENT_KIND},
    workspace::*,
};
use collagene_model::{Alphabet, Molecule, Topology};
fn workspace(bytes: Vec<u8>) -> Workspace {
    let mut w = Workspace::from_documents(
        "workspace",
        vec![Molecule::dna("reference", "ACGT", Topology::Linear).unwrap()],
        0,
    )
    .unwrap();
    w.attachments.insert(
        blake3::hash(&bytes).to_hex().to_string(),
        WorkspaceAttachment {
            kind: WORKSPACE_ATTACHMENT_KIND.into(),
            name: "Imported alignment".into(),
            media_type: "application/json".into(),
            bytes,
        },
    );
    w
}
#[test]
fn exact_original_alignment_and_row_identity_survive_archive_exchange() {
    let original = b">one\nA-CG\n>two\nAT-G\n";
    let d = alignment_formats::import_alignment(
        original,
        AlignmentFormat::AlignedFasta,
        Alphabet::Dna,
        "aligned",
    )
    .unwrap();
    let bytes = alignment_formats::save_alignment(&d).unwrap();
    let w = workspace(bytes.clone());
    let archive = save_workspace(&w).unwrap();
    let loaded = load_workspace(&archive).unwrap();
    let attachment = loaded.attachments.values().next().unwrap();
    assert_eq!(attachment.bytes, bytes);
    let d2 = alignment_formats::load_alignment(&attachment.bytes).unwrap();
    assert_eq!(d2, d);
    assert_eq!(d2.original_bytes().unwrap(), original);
}
#[test]
fn recognized_typed_json_tampering_is_rejected_even_when_attachment_hash_matches() {
    for value in [
        serde_json::json!({}),
        serde_json::json!({"schema_version":99}),
    ] {
        let bytes = serde_json::to_vec(&value).unwrap();
        let w = workspace(bytes);
        assert!(w.validate().is_err());
        assert!(save_workspace(&w).is_err());
    }
    let original = b">one\nA-CG\n";
    let d = alignment_formats::import_alignment(
        original,
        AlignmentFormat::AlignedFasta,
        Alphabet::Dna,
        "aligned",
    )
    .unwrap();
    let mut value = serde_json::to_value(d).unwrap();
    value["source"]["original_base64"] = serde_json::json!("bad");
    let w = workspace(serde_json::to_vec(&value).unwrap());
    assert!(w.validate().is_err());
}
#[test]
fn opaque_unknown_attachments_remain_supported_but_known_alignment_media_type_is_explicit() {
    let mut w = workspace(b"{}".to_vec());
    w.attachments.values_mut().next().unwrap().kind = "future-format".into();
    assert!(w.validate().is_ok());
    let original = b">one\nACG\n";
    let d = alignment_formats::import_alignment(
        original,
        AlignmentFormat::AlignedFasta,
        Alphabet::Dna,
        "aligned",
    )
    .unwrap();
    let mut w = workspace(alignment_formats::save_alignment(&d).unwrap());
    w.attachments.values_mut().next().unwrap().media_type = "text/plain".into();
    assert!(w.validate().is_err());
}

#[test]
fn archive_load_rejects_invalid_recognized_alignment_even_with_consistent_hash_and_zip_crc() {
    use std::io::{Cursor, Read, Write};
    let mut w = workspace(b"{}".to_vec());
    w.attachments.values_mut().next().unwrap().kind = "opaque".into();
    let archive = save_workspace(&w).unwrap();
    let mut input = zip::ZipArchive::new(Cursor::new(archive)).unwrap();
    let mut output = zip::ZipWriter::new(Cursor::new(Vec::new()));
    for i in 0..input.len() {
        let mut entry = input.by_index(i).unwrap();
        let name = entry.name().to_string();
        let mut bytes = Vec::new();
        entry.read_to_end(&mut bytes).unwrap();
        if name == "manifest.json" {
            let mut value: serde_json::Value = serde_json::from_slice(&bytes).unwrap();
            for attachment in value["attachments"].as_object_mut().unwrap().values_mut() {
                attachment["kind"] = serde_json::json!(WORKSPACE_ATTACHMENT_KIND);
            }
            bytes = serde_json::to_vec(&value).unwrap();
        }
        output
            .start_file(
                name,
                zip::write::SimpleFileOptions::default()
                    .compression_method(zip::CompressionMethod::Stored),
            )
            .unwrap();
        output.write_all(&bytes).unwrap();
    }
    let archive = output.finish().unwrap().into_inner();
    assert!(load_workspace(&archive).is_err());
}
