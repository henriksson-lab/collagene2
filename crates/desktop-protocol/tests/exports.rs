use collagene_desktop_protocol as wire;
use wire::library::{Blob, Digest, Text};
fn request() -> wire::RequestEnvelope {
    wire::RequestEnvelope {
        request_id: 7,
        operation_id: wire::OperationId::new("op-00000000000000000000000000000007").unwrap(),
        command: wire::Command::ExportFile {
            suggested_name: Text::new("source.gb").unwrap(),
            bytes: Blob::new(b"ACGT\r\n\0\xff".to_vec()).unwrap(),
        },
    }
}
#[test]
fn appended_export_tags_and_small_frame_control_deferral_preserve_existing_protocol() {
    let req = request();
    let bytes = wire::encode_request(&req).unwrap();
    assert_eq!(wire::decode_request(&bytes).unwrap(), req);
    assert!(wire::export_request_hint(&bytes));
    assert!(wire::decode_control(&bytes).unwrap().is_none());
    let ((_, _, tag), _) =
        postcard::take_from_bytes::<(u64, &str, u32)>(&bytes[wire::HEADER_BYTES..]).unwrap();
    assert_eq!(tag, 11);
    let mut old = req.clone();
    old.command = wire::Command::Bootstrap;
    let bytes = wire::encode_request(&old).unwrap();
    assert!(!wire::export_request_hint(&bytes));
    let ((_, _, tag), _) =
        postcard::take_from_bytes::<(u64, &str, u32)>(&bytes[wire::HEADER_BYTES..]).unwrap();
    assert_eq!(tag, 0);
    assert!(wire::export_request_hint(b"DTQ\0").eq(&false));
}
#[test]
fn export_receipt_is_bound_to_exact_bytes_length_and_generated_portable_name() {
    let req = request();
    let wire::Command::ExportFile {
        suggested_name,
        bytes,
    } = &req.command
    else {
        unreachable!()
    };
    let reply = wire::ResponseEnvelope {
        request_id: req.request_id,
        operation_id: req.operation_id.clone(),
        outcome: wire::Outcome::Success(wire::Success::Exported {
            name: wire::export_file_name(&req.operation_id, suggested_name).unwrap(),
            hash: Digest::of(bytes),
            bytes: bytes.len() as u64,
        }),
    };
    let frame = wire::encode_response(&reply).unwrap();
    wire::decode_response_for(&frame, &req).unwrap();
    for mode in 0..3 {
        let mut wrong = reply.clone();
        let wire::Outcome::Success(wire::Success::Exported { name, hash, bytes }) =
            &mut wrong.outcome
        else {
            unreachable!()
        };
        match mode {
            0 => *hash = Digest::of(b"different"),
            1 => *bytes += 1,
            _ => *name = Text::new("arbitrary.gb").unwrap(),
        }
        assert!(wrong.validate_for(&req).is_err());
    }
    for invalid in ["", "../outside.gb", "C:\\export.gb", "x/y", "line\nbreak"] {
        assert!(wire::validate_export_name(invalid).is_err());
    }
    let unicode = "ä".repeat(128);
    let name = wire::export_file_name(&req.operation_id, &unicode).unwrap();
    assert!(name.len() < 256 && name.is_ascii());
    let long_name = format!("{}.cg2design", "construct".repeat(25));
    let name = wire::export_file_name(&req.operation_id, &long_name).unwrap();
    assert!(
        name.ends_with(".cg2design"),
        "long names must retain the exchange format extension"
    );
    assert!(name.len() < 256);
    assert!(Blob::<{ wire::MAX_EXPORT_BYTES }>::new(vec![0; wire::MAX_EXPORT_BYTES + 1]).is_err());
}
