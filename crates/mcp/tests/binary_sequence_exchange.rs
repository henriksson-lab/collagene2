//! Binary exchange through the real MCP SDK; hand-authored packet/trace sources.
use base64::{Engine as _, engine::general_purpose::STANDARD};
use collagene_mcp::CollageneServer;
use rmcp::{RoleClient, ServiceExt, model::CallToolRequestParams, service::RunningService};
use serde_json::{Value, json};
type Client = RunningService<RoleClient, ()>;
async fn call(c: &Client, name: &str, args: Value) -> Value {
    let r = c
        .call_tool(
            CallToolRequestParams::new(name.to_string())
                .with_arguments(args.as_object().unwrap().clone()),
        )
        .await
        .unwrap();
    assert_ne!(r.is_error, Some(true), "{name}: {:?}", r.content);
    r.structured_content.unwrap()
}
async fn reject(c: &Client, name: &str, args: Value) {
    let r = c
        .call_tool(
            CallToolRequestParams::new(name.to_string())
                .with_arguments(args.as_object().unwrap().clone()),
        )
        .await;
    assert!(r.is_err() || r.unwrap().is_error == Some(true));
}
async fn setup() -> (Client, tokio::task::JoinHandle<()>) {
    let (a, b) = tokio::io::duplex(65536);
    let server = tokio::spawn(async move {
        CollageneServer::default()
            .serve(a)
            .await
            .unwrap()
            .waiting()
            .await
            .unwrap();
    });
    (().serve(b).await.unwrap(), server)
}
fn document(result: &Value) -> Value {
    let d = &result["documents"][0];
    json!({"id":d["id"],"revision":d["revision"]})
}
async fn bytes(c: &Client, result: &Value) -> Vec<u8> {
    let chunk = call(
        c,
        "read_artifact",
        json!({"artifact_id":result["artifact_id"],"start":0,"max_bytes":16384}),
    )
    .await;
    assert_eq!(chunk["done"], true);
    STANDARD.decode(chunk["text"].as_str().unwrap()).unwrap()
}
#[derive(Clone)]
struct Field {
    tag: [u8; 4],
    number: u32,
    kind: u16,
    size: u16,
    count: u32,
    data: Vec<u8>,
}
fn field(tag: &[u8; 4], number: u32, kind: u16, size: u16, data: Vec<u8>) -> Field {
    Field {
        tag: *tag,
        number,
        kind,
        size,
        count: (data.len() / size as usize) as u32,
        data,
    }
}
fn shorts(values: &[i16]) -> Vec<u8> {
    values.iter().flat_map(|n| n.to_be_bytes()).collect()
}
fn abif_fixture(calls: &[u8], quality: bool) -> Vec<u8> {
    let positions: Vec<i16> = if calls.len() == 4 {
        vec![0, 2, 4, 7]
    } else {
        vec![0, 2, 4, 6, 7]
    };
    let mut fields = vec![
        field(b"FWO_", 1, 2, 1, b"TGCA".to_vec()),
        field(b"DATA", 9, 4, 2, shorts(&[-1, 0, 5, 7, 9, 11, 13, 15])),
        field(
            b"DATA",
            10,
            4,
            2,
            shorts(&[100, 101, 102, 103, 104, 105, 106, 107]),
        ),
        field(
            b"DATA",
            11,
            4,
            2,
            shorts(&[200, 201, 202, 203, 204, 205, 206, 207]),
        ),
        field(
            b"DATA",
            12,
            4,
            2,
            shorts(&[300, 301, 302, 303, 304, 305, 306, 307]),
        ),
        field(b"PBAS", 2, 2, 1, calls.to_vec()),
        field(b"PLOC", 2, 4, 2, shorts(&positions)),
        field(b"SMPL", 1, 18, 1, b"\x03syn".to_vec()),
        field(b"RevC", 1, 4, 2, shorts(&[0])),
    ];
    if quality {
        fields.push(field(
            b"PCON",
            2,
            2,
            1,
            (0..calls.len()).map(|i| (i * 10 + 10) as u8).collect(),
        ));
    }
    let mut bytes = vec![0; 128 + 28 * fields.len()];
    bytes[..4].copy_from_slice(b"ABIF");
    bytes[4..6].copy_from_slice(&101u16.to_be_bytes());
    bytes[6..10].copy_from_slice(b"tdir");
    bytes[10..14].copy_from_slice(&1u32.to_be_bytes());
    bytes[14..16].copy_from_slice(&1023u16.to_be_bytes());
    bytes[16..18].copy_from_slice(&28u16.to_be_bytes());
    bytes[22..26].copy_from_slice(&((fields.len() * 28) as u32).to_be_bytes());
    bytes[18..22].copy_from_slice(&(fields.len() as u32).to_be_bytes());
    bytes[26..30].copy_from_slice(&128u32.to_be_bytes());
    for (i, f) in fields.iter().enumerate() {
        let p = 128 + 28 * i;
        bytes[p..p + 4].copy_from_slice(&f.tag);
        bytes[p + 4..p + 8].copy_from_slice(&f.number.to_be_bytes());
        bytes[p + 8..p + 10].copy_from_slice(&f.kind.to_be_bytes());
        bytes[p + 10..p + 12].copy_from_slice(&f.size.to_be_bytes());
        bytes[p + 12..p + 16].copy_from_slice(&f.count.to_be_bytes());
        bytes[p + 16..p + 20].copy_from_slice(&(f.data.len() as u32).to_be_bytes());
        if f.data.len() <= 4 {
            bytes[p + 20..p + 20 + f.data.len()].copy_from_slice(&f.data);
        } else {
            let offset = bytes.len() as u32;
            bytes[p + 20..p + 24].copy_from_slice(&offset.to_be_bytes());
            bytes.extend_from_slice(&f.data);
        }
    }
    bytes
}
fn word(bytes: &mut [u8], at: usize, value: u32) {
    bytes[at..at + 4].copy_from_slice(&value.to_be_bytes());
}
fn scf_fixture(wide: bool, extras: bool, code: u32, calls: &[u8; 4]) -> Vec<u8> {
    let channels: [[u16; 8]; 4] = if wide {
        [
            [0, 65535, 0, 32768, 0, 1, 65534, 5],
            [20, 21, 22, 23, 24, 25, 26, 27],
            [1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000],
            [65535, 0, 65535, 0, 65535, 0, 65535, 0],
        ]
    } else {
        [
            [0, 255, 0, 128, 0, 1, 254, 5],
            [20, 21, 22, 23, 24, 25, 26, 27],
            [8, 8, 8, 8, 8, 8, 8, 8],
            [255, 0, 255, 0, 255, 0, 255, 0],
        ]
    };
    let precision = if wide { 2 } else { 1 };
    let base_offset = 128 + 32 * precision;
    let comment = if extras {
        b"NAME=actual\nMACH=retained instrument\n\0".as_slice()
    } else {
        b"NAME=actual\n\0".as_slice()
    };
    let comment_offset = base_offset + 48;
    let mut bytes = vec![0; comment_offset + comment.len()];
    bytes[..4].copy_from_slice(b".scf");
    word(&mut bytes, 4, 8);
    word(&mut bytes, 8, 128);
    word(&mut bytes, 12, 4);
    word(&mut bytes, 24, base_offset as u32);
    word(&mut bytes, 28, comment.len() as u32);
    word(&mut bytes, 32, comment_offset as u32);
    bytes[36..40].copy_from_slice(b"2.02");
    word(&mut bytes, 40, precision as u32);
    word(&mut bytes, 44, code);
    let mut at = 128;
    for sample in 0..8 {
        for channel in &channels {
            let value = channel[sample];
            if precision == 1 {
                bytes[at] = value as u8;
            } else {
                bytes[at..at + 2].copy_from_slice(&value.to_be_bytes());
            }
            at += precision;
        }
    }
    for (i, peak) in [0, 2, 4, 7].into_iter().enumerate() {
        let at = base_offset + i * 12;
        word(&mut bytes, at, peak);
        bytes[at + 4..at + 8].copy_from_slice(&[
            1 + i as u8,
            11 + i as u8,
            21 + i as u8,
            31 + i as u8,
        ]);
        bytes[at + 8] = calls[i];
    }
    bytes[comment_offset..].copy_from_slice(comment);
    if extras {
        let at = bytes.len();
        word(&mut bytes, 48, 4);
        word(&mut bytes, 52, at as u32);
        bytes.extend_from_slice(&[99, 88, 77, 66]);
        word(&mut bytes, 16, 2);
        bytes[56] = 9;
        bytes[base_offset + 9] = 7;
        bytes.push(42);
    }
    bytes
}

fn packet(kind: u8, data: &[u8]) -> Vec<u8> {
    let mut b = vec![kind];
    b.extend_from_slice(&(data.len() as u32).to_be_bytes());
    b.extend_from_slice(data);
    b
}
fn snapgene_fixture(protein: bool) -> Vec<u8> {
    let mut cookie = b"SnapGene".to_vec();
    cookie.extend_from_slice(&(if protein { 2u16 } else { 1 }).to_be_bytes());
    cookie.extend_from_slice(&15u16.to_be_bytes());
    cookie.extend_from_slice(&15u16.to_be_bytes());
    let mut bytes = packet(9, &cookie);
    bytes.extend(packet(
        if protein { 21 } else { 0 },
        if protein { b"\0MKWVTF" } else { b"\x03ACGTTG" },
    ));
    bytes
}
async fn imported(c: &Client, raw: &[u8], format: Value) -> Value {
    call(
        c,
        "import_sequence_file",
        json!({"data_base64":STANDARD.encode(raw),"format":format}),
    )
    .await
}
async fn molecule(c: &Client, d: &Value) -> collagene_core::Molecule {
    let out = call(c, "export_document", d.clone()).await;
    collagene_core::protocol::load_document(&bytes(c, &out).await).unwrap()
}
#[tokio::test]
async fn binary_sequence_profiles_keep_originals_and_export_reviewed_reference_files() {
    let (c, server) = setup().await;
    let tools = c.list_all_tools().await.unwrap();
    assert_eq!(tools.len(), collagene_mcp::CURRENT_TOOL_COUNT);
    for name in ["import_sequence_file", "export_xdna", "export_trace_scf"] {
        assert_eq!(
            tools.iter().find(|t| t.name == name).unwrap().input_schema["additionalProperties"],
            false
        );
    }
    for protein in [false, true] {
        let raw = snapgene_fixture(protein);
        let result = imported(&c, &raw, json!("SnapGene")).await;
        let d = document(&result);
        let m = molecule(&c, &d).await;
        assert_eq!(
            m.sequence.as_str(),
            if protein { "MKWVTF" } else { "ACGTTG" }
        );
        assert_eq!(
            m.sequence.alphabet(),
            if protein {
                collagene_core::Alphabet::Protein
            } else {
                collagene_core::Alphabet::Dna
            }
        );
        assert_eq!(
            collagene_core::snapgene::original_snapgene_bytes(&m).unwrap(),
            raw
        );
        let before = call(&c, "export_document", d.clone()).await;
        let review = call(&c, "export_xdna", json!({"document":d})).await;
        assert_eq!(review["state"], "review_required");
        assert!(review["artifact_id"].is_null());
        let output = call(
            &c,
            "export_xdna",
            json!({"document":d,"accept_omissions":true}),
        )
        .await;
        let emitted = bytes(&c, &output).await;
        assert_eq!(emitted[0], 0);
        assert_eq!(emitted[1], if protein { 4 } else { 1 });
        let roundtrip = imported(
            &c,
            &emitted,
            json!({"XdnaV0":{"name":"Explicit file name"}}),
        )
        .await;
        let r = molecule(&c, &document(&roundtrip)).await;
        assert_eq!(r.sequence, m.sequence);
        assert_eq!(r.topology, m.topology);
        assert_eq!(
            collagene_core::legacy_formats::original_legacy_bytes(&r).unwrap(),
            emitted
        );
        assert_eq!(call(&c, "export_document", d.clone()).await, before);
        let snapshot = bytes(&c, &before).await;
        let native = imported(&c, &snapshot, json!("NativeDocument")).await;
        assert_eq!(document(&native), d);
    }
    c.cancel().await.unwrap();
    server.await.unwrap();
}
#[tokio::test]
async fn scf_exports_preserve_actual_signal_accuracy_and_refuse_stale_calls_or_abif_scores() {
    let (c, server) = setup().await;
    let raw = scf_fixture(false, true, 2, b"ACRN");
    let result = imported(&c, &raw, json!("Scf")).await;
    let d = document(&result);
    let m = molecule(&c, &d).await;
    assert_eq!(
        collagene_core::traces::original_trace_bytes(&m).unwrap(),
        raw
    );
    let before = call(&c, "export_document", d.clone()).await;
    let review = call(
        &c,
        "export_trace_scf",
        json!({"document":d,"precision":"Smallest"}),
    )
    .await;
    assert_eq!(review["state"], "review_required");
    assert!(review["artifact_id"].is_null());
    assert!(!review["losses"].as_array().unwrap().is_empty());
    let out = call(
        &c,
        "export_trace_scf",
        json!({"document":d,"precision":"Smallest","accept_omissions":true}),
    )
    .await;
    assert_eq!(out["context_revision"], review["context_revision"]);
    let emitted = bytes(&c, &out).await;
    assert_eq!(&emitted[36..40], b"3.00");
    assert_eq!(&emitted[128..136], &[0, 255, 2, 127, 0, 129, 252, 10]);
    let original = collagene_core::traces::import_trace(&raw).unwrap();
    let regenerated = collagene_core::traces::import_trace(&emitted).unwrap();
    assert_eq!(original.channels, regenerated.channels);
    assert_eq!(original.qualities, regenerated.qualities);
    assert_eq!(original.peak_positions, regenerated.peak_positions);
    assert_eq!(original.called_symbols, regenerated.called_symbols);
    assert_eq!(call(&c, "export_document", d.clone()).await, before);
    let reversed = original
        .reverse_complement()
        .unwrap()
        .attached_molecule()
        .unwrap();
    let reversed = imported(
        &c,
        &collagene_core::protocol::save_document(&reversed).unwrap(),
        json!("NativeDocument"),
    )
    .await;
    let rd = document(&reversed);
    let out = call(
        &c,
        "export_trace_scf",
        json!({"document":rd,"precision":"SixteenBit","accept_omissions":true}),
    )
    .await;
    let actual = collagene_core::traces::import_trace(&bytes(&c, &out).await).unwrap();
    let expected = original.reverse_complement().unwrap();
    assert_eq!(actual.channels, expected.channels);
    assert_eq!(actual.qualities, expected.qualities);
    assert_eq!(actual.peak_positions, expected.peak_positions);
    assert_eq!(actual.sequence, expected.sequence);
    let edited = m.replace(0, 1, "T").unwrap();
    let edited = imported(
        &c,
        &collagene_core::protocol::save_document(&edited).unwrap(),
        json!("NativeDocument"),
    )
    .await;
    for ack in [false, true] {
        reject(
            &c,
            "export_trace_scf",
            json!({"document":document(&edited),"precision":"Smallest","accept_omissions":ack}),
        )
        .await;
    }
    let abif = abif_fixture(b"ACGT", true);
    let imported_abif = imported(&c, &abif, json!("Abif")).await;
    let ad = document(&imported_abif);
    let parsed = molecule(&c, &ad).await;
    assert_eq!(parsed.sequence.as_str(), "ACGT");
    assert_eq!(
        collagene_core::traces::original_trace_bytes(&parsed).unwrap(),
        abif
    );
    reject(
        &c,
        "export_trace_scf",
        json!({"document":ad,"precision":"Smallest","accept_omissions":true}),
    )
    .await;
    c.cancel().await.unwrap();
    server.await.unwrap();
}
#[tokio::test]
async fn format_mismatches_malformed_payloads_and_source_tampering_leave_sessions_unchanged() {
    let (c, server) = setup().await;
    let raw = snapgene_fixture(false);
    let imported_doc = imported(&c, &raw, json!("SnapGene")).await;
    let d = document(&imported_doc);
    let before = call(&c, "list_documents", json!({})).await;
    for (data, profile) in [
        (raw.clone(), json!("Scf")),
        (scf_fixture(false, false, 2, b"ACRN"), json!("Abif")),
        (b">text\nACGT\n".to_vec(), json!("SnapGene")),
        (raw.clone(), json!("NativeDocument")),
        (vec![0; 16], json!({"XdnaV0":{"name":"literal"}})),
    ] {
        reject(
            &c,
            "import_sequence_file",
            json!({"data_base64":STANDARD.encode(data),"format":profile}),
        )
        .await;
    }
    for name in ["", "bad\nname"] {
        reject(
            &c,
            "import_sequence_file",
            json!({"data_base64":STANDARD.encode(&raw),"format":{"XdnaV0":{"name":name}}}),
        )
        .await;
    }
    reject(
        &c,
        "import_sequence_file",
        json!({"data_base64":"!not-base64","format":"SnapGene"}),
    )
    .await;
    reject(
        &c,
        "import_sequence_file",
        json!({"data_base64":STANDARD.encode(&raw),"format":"Unknown"}),
    )
    .await;
    reject(
        &c,
        "import_sequence_file",
        json!({"data_base64":STANDARD.encode(&raw),"format":"SnapGene","path":"/tmp/untrusted"}),
    )
    .await;
    reject(
        &c,
        "export_trace_scf",
        json!({"document":d,"precision":"Smallest","accept_omissions":true}),
    )
    .await;
    reject(
        &c,
        "import_sequence_file",
        json!({"data_base64":"A".repeat((8*1024*1024usize).div_ceil(3)*4+1),"format":"SnapGene"}),
    )
    .await;
    let mut stale = d;
    stale["revision"] = json!("f".repeat(64));
    reject(
        &c,
        "export_xdna",
        json!({"document":stale,"accept_omissions":true}),
    )
    .await;
    assert_eq!(call(&c, "list_documents", json!({})).await, before);
    c.cancel().await.unwrap();
    server.await.unwrap();
}
