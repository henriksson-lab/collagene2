//! Real SDK routing: no direct tool-method shortcuts or mock quality/verification verdicts.
use base64::{Engine, engine::general_purpose::STANDARD};
use collagene_core::{
    Molecule, Topology,
    oligo_formats::{ColumnMapping, OligoImportOptions},
    trace_validation::TraceValidationOptions,
    transfer::Candidate,
    verification::{DesignPlan, Operation, PlanInput, Source},
    workspace::{self, RevisionEvent, Workspace, WorkspaceAttachment},
};
use collagene_mcp::CollageneServer;
use rmcp::{RoleClient, ServiceExt, model::CallToolRequestParams, service::RunningService};
use serde_json::{Value, json};

type Client = RunningService<RoleClient, ()>;
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
async fn call(client: &Client, name: &str, args: Value) -> Value {
    let response = client
        .call_tool(
            CallToolRequestParams::new(name.to_owned())
                .with_arguments(args.as_object().unwrap().clone()),
        )
        .await
        .unwrap();
    assert_ne!(
        response.is_error,
        Some(true),
        "{name}: {:?}",
        response.content
    );
    response.structured_content.unwrap()
}
async fn reject(client: &Client, name: &str, args: Value) {
    let result = client
        .call_tool(
            CallToolRequestParams::new(name.to_owned())
                .with_arguments(args.as_object().unwrap().clone()),
        )
        .await;
    assert!(
        result.is_err() || result.unwrap().is_error == Some(true),
        "accepted invalid {name}"
    );
}
async fn artifact(client: &Client, id: &Value) -> Vec<u8> {
    let mut offset = 0;
    let mut text = String::new();
    let mut encoding = String::new();
    loop {
        let chunk = call(
            client,
            "read_artifact",
            json!({"artifact_id":id,"start":offset,"max_bytes":16384}),
        )
        .await;
        encoding.clone_from(&chunk["encoding"].as_str().unwrap().to_string());
        text.push_str(chunk["text"].as_str().unwrap());
        let next = chunk["next_start"].as_u64().unwrap() as usize;
        if chunk["done"] == true {
            break;
        }
        assert!(next > offset);
        offset = next;
    }
    if encoding == "base64" {
        STANDARD.decode(text).unwrap()
    } else {
        text.into_bytes()
    }
}
async fn upload(client: &Client, bytes: &[u8]) -> Value {
    let begin = call(
        client,
        "begin_workspace_upload",
        json!({"expected_bytes":bytes.len(),"blake3":blake3::hash(bytes).to_hex().to_string()}),
    )
    .await;
    let id = &begin["upload_id"];
    for (i, chunk) in bytes.chunks(256 * 1024).enumerate() {
        let result = call(
            client,
            "append_workspace_upload",
            json!({"upload_id":id,"offset":i*256*1024,"data_base64":STANDARD.encode(chunk)}),
        )
        .await;
        assert_eq!(result["next_offset"], (i * 256 * 1024 + chunk.len()) as u64);
    }
    call(client, "finalize_workspace_upload", json!({"upload_id":id})).await
}
fn archive_fixture(documents: usize) -> Workspace {
    let mut all = (0..documents)
        .map(|i| Molecule::dna(format!("sequence {i}"), "CGATAC", Topology::Linear).unwrap())
        .collect::<Vec<_>>();
    all[0]
        .metadata
        .insert("opaque:source".into(), STANDARD.encode([0, 128, 255, 10]));
    let mut archive = Workspace::from_documents("History / candidates / evidence", all, 0).unwrap();
    let root = archive.documents[0].head.clone();
    let root_event = archive.documents[0].head_event.clone().unwrap();
    let changed = archive.revisions[&root].replace(2, 3, "C").unwrap();
    let changed_revision = changed.revision();
    archive.revisions.insert(changed_revision.clone(), changed);
    let event = RevisionEvent::new(
        "document-1",
        &changed_revision,
        vec![root_event.clone()],
        "Explicit base edit",
    )
    .unwrap();
    archive.documents[0].head = changed_revision.clone();
    archive.documents[0].undo = vec![root.clone()];
    archive.documents[0].head_event = Some(event.id.clone());
    archive.documents[0].undo_events = vec![root_event];
    archive.history.push(event);
    let candidate = Candidate::new(DesignPlan {
        name: "Unverified candidate, independently bound target".into(),
        target_revision: changed_revision,
        inputs: vec![PlanInput {
            id: "source".into(),
            revision: root.clone(),
            molecule: archive.revisions[&root].clone(),
        }],
        steps: vec![Operation::Overlap {
            fragments: vec![Source::Input("source".into())],
            overlaps: vec![],
            circular: false,
        }],
    });
    archive.active_candidate = Some(candidate.plan_revision.clone());
    archive.candidates.push(candidate);
    let attachment = WorkspaceAttachment {
        kind: "source-evidence".into(),
        name: "opaque.bin".into(),
        media_type: "application/octet-stream".into(),
        bytes: vec![0, 255, 128, 10, 1],
    };
    archive.attachments.insert(
        blake3::hash(&attachment.bytes).to_hex().to_string(),
        attachment,
    );
    archive.validate().unwrap();
    archive
}

#[tokio::test]
async fn native_workspace_preserves_1000_documents_exact_history_candidates_attachments_and_pagination()
 {
    let (client, server) = setup().await;
    let tools = client.list_all_tools().await.unwrap();
    for name in [
        "validate_reads",
        "import_fastq",
        "export_fastq",
        "import_oligos",
        "export_oligos",
        "import_workspace_archive",
        "export_workspace_archive",
        "begin_workspace_upload",
        "append_workspace_upload",
        "finalize_workspace_upload",
        "abort_workspace_upload",
        "list_workspaces",
    ] {
        let tool = tools.iter().find(|t| t.name == name).unwrap();
        assert_eq!(tool.input_schema["additionalProperties"], false, "{name}");
    }
    let original = archive_fixture(1000);
    let bytes = workspace::save_workspace(&original).unwrap();
    let imported = call(
        &client,
        "import_workspace_archive",
        json!({"data_base64":STANDARD.encode(&bytes)}),
    )
    .await;
    assert_eq!(imported["workspace"]["document_count"], 1000);
    assert_eq!(imported["workspace"]["candidates"][0]["verified"], false);
    let reference =
        json!({"id":imported["workspace"]["id"],"revision":imported["workspace"]["revision"]});
    let page = call(&client, "list_documents", json!({})).await;
    assert_eq!(page["documents"].as_array().unwrap().len(), 128);
    assert_eq!(page["total"], 1001);
    assert_eq!(page["next_start"], 128);
    let last = call(&client, "list_documents", json!({"start":896,"limit":128})).await;
    assert_eq!(last["documents"].as_array().unwrap().len(), 105);
    assert!(last["next_start"].is_null());
    reject(&client, "list_documents", json!({"limit":129})).await;
    let exported = call(&client, "export_workspace_archive", reference.clone()).await;
    let restored =
        workspace::load_workspace(&artifact(&client, &exported["artifact_id"]).await).unwrap();
    assert_eq!(
        serde_json::to_value(&restored).unwrap(),
        serde_json::to_value(&original).unwrap()
    );
    assert_eq!(exported["candidates_verified"], false);
    let duplicate = call(
        &client,
        "import_workspace_archive",
        json!({"data_base64":STANDARD.encode(&bytes)}),
    )
    .await;
    assert_eq!(duplicate["deduplicated"], true);
    assert_eq!(
        call(&client, "list_documents", json!({})).await["total"],
        1001
    );
    let mut stale = reference;
    stale["revision"] = json!("stale");
    reject(&client, "export_workspace_archive", stale).await;
    reject(
        &client,
        "verify_plan",
        json!({"plan_revision":original.candidates[0].plan_revision}),
    )
    .await;
    assert_eq!(
        call(&client, "list_workspaces", json!({})).await["workspaces"]
            .as_array()
            .unwrap()
            .len(),
        1
    );
    // Original planning snapshot import stays available alongside native archives.
    let plan = collagene_core::transfer::PlanningWorkspace {
        target_revision: original.documents[0].head.clone(),
        documents: vec![original.revisions[&original.documents[0].head].clone()],
    };
    let old = collagene_core::transfer::save_workspace(&plan).unwrap();
    assert_eq!(
        call(
            &client,
            "import_workspace",
            json!({"data_base64":STANDARD.encode(old)})
        )
        .await["target"]["revision"],
        plan.target_revision
    );
    client.cancel().await.unwrap();
    server.await.unwrap();
}

#[tokio::test]
async fn chunked_native_archive_larger_than_eight_mib_roundtrips_actual_binary_evidence() {
    let (client, server) = setup().await;
    let mut original = archive_fixture(1);
    // Deterministic xorshift data avoids a misleading tiny compressed test fixture.
    let mut state = 0x71a95a32u32;
    let mut data = Vec::with_capacity(9 * 1024 * 1024);
    for _ in 0..9 * 1024 * 1024 {
        state ^= state << 13;
        state ^= state >> 17;
        state ^= state << 5;
        data.push(state as u8);
    }
    let hash = blake3::hash(&data).to_hex().to_string();
    original.attachments.insert(
        hash.clone(),
        WorkspaceAttachment {
            kind: "source-evidence".into(),
            name: "large-original.bin".into(),
            media_type: "application/octet-stream".into(),
            bytes: data,
        },
    );
    let bytes = workspace::save_workspace(&original).unwrap();
    assert!(bytes.len() > 8 * 1024 * 1024);
    assert!(bytes.len() < 64 * 1024 * 1024);
    reject(
        &client,
        "import_workspace_archive",
        json!({"data_base64":STANDARD.encode(&bytes)}),
    )
    .await;
    let imported = upload(&client, &bytes).await;
    let exported = call(
        &client,
        "export_workspace_archive",
        json!({"id":imported["workspace"]["id"],"revision":imported["workspace"]["revision"]}),
    )
    .await;
    // Read non-four-byte-aligned base64 offsets as well as ordinary pagination.
    let a = call(
        &client,
        "read_artifact",
        json!({"artifact_id":exported["artifact_id"],"start":1,"max_bytes":7}),
    )
    .await;
    let b = call(
        &client,
        "read_artifact",
        json!({"artifact_id":exported["artifact_id"],"start":0,"max_bytes":8}),
    )
    .await;
    assert_eq!(
        a["text"].as_str().unwrap(),
        &b["text"].as_str().unwrap()[1..8]
    );
    let restored =
        workspace::load_workspace(&artifact(&client, &exported["artifact_id"]).await).unwrap();
    assert_eq!(
        restored.attachments[&hash].bytes,
        original.attachments[&hash].bytes
    );
    assert_eq!(restored.history.len(), 2);
    assert_eq!(
        restored.candidates[0].plan_revision,
        original.candidates[0].plan_revision
    );
    client.cancel().await.unwrap();
    server.await.unwrap();
}

#[tokio::test]
async fn upload_offsets_hashes_atomicity_abort_and_shared_reservations_are_enforced() {
    let (client, server) = setup().await;
    reject(
        &client,
        "begin_workspace_upload",
        json!({"expected_bytes":64*1024*1024+1,"blake3":"0".repeat(64)}),
    )
    .await;
    reject(
        &client,
        "begin_workspace_upload",
        json!({"expected_bytes":1,"blake3":"A".repeat(64)}),
    )
    .await;
    let large = call(
        &client,
        "begin_workspace_upload",
        json!({"expected_bytes":64*1024*1024,"blake3":"0".repeat(64)}),
    )
    .await;
    reject(
        &client,
        "import_document",
        json!({"text":">no partial import\nACGT"}),
    )
    .await;
    assert_eq!(call(&client, "list_documents", json!({})).await["total"], 0);
    call(
        &client,
        "abort_workspace_upload",
        json!({"upload_id":large["upload_id"]}),
    )
    .await;
    let bytes = workspace::save_workspace(&archive_fixture(1)).unwrap();
    let begin = call(
        &client,
        "begin_workspace_upload",
        json!({"expected_bytes":bytes.len(),"blake3":blake3::hash(&bytes).to_hex().to_string()}),
    )
    .await;
    let id = &begin["upload_id"];
    reject(
        &client,
        "append_workspace_upload",
        json!({"upload_id":id,"offset":1,"data_base64":STANDARD.encode(&bytes[..3])}),
    )
    .await;
    reject(
        &client,
        "append_workspace_upload",
        json!({"upload_id":id,"offset":0,"data_base64":""}),
    )
    .await;
    reject(
        &client,
        "finalize_workspace_upload",
        json!({"upload_id":id}),
    )
    .await;
    call(
        &client,
        "append_workspace_upload",
        json!({"upload_id":id,"offset":0,"data_base64":STANDARD.encode(&bytes[..3])}),
    )
    .await;
    reject(
        &client,
        "append_workspace_upload",
        json!({"upload_id":id,"offset":0,"data_base64":STANDARD.encode(&bytes[..3])}),
    )
    .await;
    reject(
        &client,
        "append_workspace_upload",
        json!({"upload_id":id,"offset":3,"data_base64":STANDARD.encode(vec![0;256*1024+1])}),
    )
    .await;
    call(
        &client,
        "append_workspace_upload",
        json!({"upload_id":id,"offset":3,"data_base64":STANDARD.encode(&bytes[3..])}),
    )
    .await;
    call(
        &client,
        "finalize_workspace_upload",
        json!({"upload_id":id}),
    )
    .await;
    reject(&client, "abort_workspace_upload", json!({"upload_id":id})).await;
    let before = call(&client, "list_documents", json!({})).await["total"].clone();
    let wrong = call(
        &client,
        "begin_workspace_upload",
        json!({"expected_bytes":3,"blake3":"0".repeat(64)}),
    )
    .await;
    call(
        &client,
        "append_workspace_upload",
        json!({"upload_id":wrong["upload_id"],"offset":0,"data_base64":STANDARD.encode(b"bad")}),
    )
    .await;
    reject(
        &client,
        "finalize_workspace_upload",
        json!({"upload_id":wrong["upload_id"]}),
    )
    .await;
    let malformed = call(
        &client,
        "begin_workspace_upload",
        json!({"expected_bytes":3,"blake3":blake3::hash(b"bad").to_hex().to_string()}),
    )
    .await;
    call(&client,"append_workspace_upload",json!({"upload_id":malformed["upload_id"],"offset":0,"data_base64":STANDARD.encode(b"bad")})).await;
    reject(
        &client,
        "finalize_workspace_upload",
        json!({"upload_id":malformed["upload_id"]}),
    )
    .await;
    assert_eq!(
        call(&client, "list_documents", json!({})).await["total"],
        before
    );
    client.cancel().await.unwrap();
    server.await.unwrap();
}

#[tokio::test]
async fn explicit_fastq_quality_validation_and_oligo_collection_evidence_are_revision_bound() {
    let (client, server) = setup().await;
    reject(
        &client,
        "import_fastq",
        json!({"text":"@r\nCGATAC\n+\nIIIIII\n"}),
    )
    .await;
    let imported = call(
        &client,
        "import_fastq",
        json!({"text":"@r\nCGATAC\n+\nIIIIII\n","encoding":"Phred33"}),
    )
    .await;
    let read = json!({"id":imported["documents"][0]["id"],"revision":imported["documents"][0]["revision"]});
    let reference = call(
        &client,
        "import_document",
        json!({"text":">independent\nTTTCGATACGGG\n"}),
    )
    .await;
    let reference = json!({"id":reference["documents"][0]["id"],"revision":reference["documents"][0]["revision"]});
    let mut options = TraceValidationOptions::default();
    options.trim.min_retained_bases = 1;
    let validation = call(
        &client,
        "validate_reads",
        json!({"reference":reference,"reads":[read],"options":options}),
    )
    .await;
    assert_eq!(validation["status"], "Incomplete");
    assert_eq!(validation["eligible_coverage"], 6);
    assert_eq!(validation["consensus_preview"], "NNNCGATACNNN");
    assert_eq!(validation["reads"][0]["quality_source"], "FastqPhred33");
    let report: Value =
        serde_json::from_slice(&artifact(&client, &validation["artifact_id"]).await).unwrap();
    assert_eq!(
        report["validation"]["reference_revision"],
        reference["revision"]
    );
    assert!(report["validation"]["reads"][0]["columns"][0]["stored_sample_index"].is_null());
    let exported = call(
        &client,
        "export_fastq",
        json!({"document":read,"encoding":"Phred64"}),
    )
    .await;
    assert_eq!(
        artifact(&client, &exported["artifact_id"]).await,
        b"@r\nCGATAC\n+\nhhhhhh\n"
    );
    let mut stale = read.clone();
    stale["revision"] = json!("stale");
    reject(
        &client,
        "validate_reads",
        json!({"reference":reference,"reads":[stale],"options":options}),
    )
    .await;
    reject(
        &client,
        "export_fastq",
        json!({"document":stale,"encoding":"Phred33"}),
    )
    .await;
    let text = "Name,Sequence,Modification,Amount,Unit\r\nplain,CGATAC,,10,nmol\r\nblocked,/5Biosg/ACGT,biotin,2,uM\r\n";
    let mapping = OligoImportOptions {
        columns: ColumnMapping {
            modification: Some(2),
            amount: Some(3),
            unit: Some(4),
            ..Default::default()
        },
        ..Default::default()
    };
    let oligos = call(
        &client,
        "import_oligos",
        json!({"text":text,"options":mapping,"selection":[0]}),
    )
    .await;
    assert_eq!(oligos["row_count"], 2);
    assert_eq!(oligos["documents"].as_array().unwrap().len(), 1);
    assert_eq!(oligos["rows"][1]["dna_supported"], false);
    let n = call(&client, "list_documents", json!({})).await["total"].clone();
    reject(
        &client,
        "import_oligos",
        json!({"text":text,"options":mapping,"selection":[0,1]}),
    )
    .await;
    assert_eq!(call(&client, "list_documents", json!({})).await["total"], n);
    let exported=call(&client,"export_oligos",json!({"collection":oligos["collection"],"delimiter":"Tab","mode":"MachineData","selection":[0,1]})).await;
    assert_eq!(
        artifact(&client, &exported["original_artifact_id"]).await,
        text.as_bytes()
    );
    let table = String::from_utf8(artifact(&client, &exported["table_artifact_id"]).await).unwrap();
    assert!(table.contains("/5Biosg/ACGT"));
    let package: Value =
        serde_json::from_slice(&artifact(&client, &exported["collection_artifact_id"]).await)
            .unwrap();
    assert_eq!(package["selection"], json!([0, 1]));
    assert_eq!(package["saved_selections"][0]["source_rows"], json!([0]));
    assert_eq!(
        package["table"]["records"][1]["raw_ordered_sequence"],
        "/5Biosg/ACGT"
    );
    assert!(
        !serde_json::to_string(&report)
            .unwrap()
            .contains("ConsistentUnderRecordedThresholds")
    );
    let mut stale = oligos["collection"].clone();
    stale["revision"] = json!("stale");
    reject(
        &client,
        "export_oligos",
        json!({"collection":stale,"delimiter":"Comma","mode":"SpreadsheetSafe"}),
    )
    .await;
    client.cancel().await.unwrap();
    server.await.unwrap();
}
