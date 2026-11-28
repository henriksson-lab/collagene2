//! Public JSON-RPC coverage of explicit paired-library declarations and outputs.
#[path = "../../../core/tests/common/paired_library.rs"]
mod fixture;
#[path = "../../../core/tests/common/pdf.rs"]
mod pdf_readback;

use base64::{Engine as _, engine::general_purpose::STANDARD};
use collagene_core::{
    documents::workspace,
    material::lineage,
    planning::{
        combinatorial::{batch, release},
        plan::Operation,
    },
};
use collagene_mcp::server::CollageneServer;
use rmcp::{RoleClient, ServiceExt, model::CallToolRequestParams, service::RunningService};
use serde_json::{Value, json};
use std::io::{Cursor, Read};

type Client = RunningService<RoleClient, ()>;

async fn call(client: &Client, name: &str, args: Value) -> Value {
    let result = client
        .call_tool(
            CallToolRequestParams::new(name.to_owned())
                .with_arguments(args.as_object().unwrap().clone()),
        )
        .await
        .unwrap();
    assert_ne!(result.is_error, Some(true), "{name}: {:?}", result.content);
    result.structured_content.unwrap()
}

async fn reject(client: &Client, name: &str, args: Value) {
    match client
        .call_tool(
            CallToolRequestParams::new(name.to_owned())
                .with_arguments(args.as_object().unwrap().clone()),
        )
        .await
    {
        Err(_) => {}
        Ok(result) => assert_eq!(
            result.is_error,
            Some(true),
            "{name} accepted an invalid declaration"
        ),
    }
}

async fn artifact(client: &Client, id: &Value) -> Vec<u8> {
    let mut start = 0;
    let mut text = String::new();
    let encoding = loop {
        let chunk = call(
            client,
            "read_artifact",
            json!({"artifact_id":id,"start":start,"max_bytes":16384}),
        )
        .await;
        text.push_str(chunk["text"].as_str().unwrap());
        if chunk["done"] == true {
            break chunk["encoding"].as_str().unwrap().to_owned();
        }
        let next = chunk["next_start"].as_u64().unwrap();
        assert!(next > start);
        start = next;
    };
    match encoding.as_str() {
        "base64" => STANDARD.decode(text).unwrap(),
        "utf-8" => text.into_bytes(),
        other => panic!("unexpected artifact encoding {other}"),
    }
}

fn workspace_ref(result: &Value) -> Value {
    json!({"id":result["workspace"]["id"],"revision":result["workspace"]["revision"]})
}

#[tokio::test]
async fn paired_library_public_sdk_keeps_selected_bottom_event_and_fresh_protocol_sources() {
    let (original, _, unselected_bottom) = fixture::fixture();
    let original_value = serde_json::to_value(&original).unwrap();
    let (hash, _) = original
        .attachments
        .iter()
        .find(|(_, a)| a.kind == batch::KIND)
        .unwrap();
    let retained = batch::load(&original, hash).unwrap();
    let design = serde_json::to_value(retained.request.design.as_ref().unwrap()).unwrap();
    let batch::RowRequest::Ready(construct) = &retained.request.rows[0] else {
        unreachable!()
    };
    let construct = serde_json::to_value(construct).unwrap();
    let matrix = serde_json::to_value(&retained.request).unwrap();
    assert!(design["bins"][0]["choices"][0]["anneal_bottom"].is_object());

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
    let client: Client = ().serve(b).await.unwrap();
    let imported = call(
        &client,
        "import_workspace_archive",
        json!({"data_base64":STANDARD.encode(workspace::save_workspace(&original).unwrap())}),
    )
    .await;
    let source_ref = workspace_ref(&imported);
    let before_documents = call(&client, "list_documents", json!({})).await;
    // Explicit null is a declared unpaired, unfinished alternative. An omitted
    // field is an unsupported encoding, rather than an inferred missing operand.
    let mut unfinished = design.clone();
    unfinished["bins"][0]["choices"][1]["anneal_bottom"] = Value::Null;
    unfinished["bins"][0]["choices"][1]["preparation"] = Value::Null;
    let page = call(
        &client,
        "page_combinatorial_design",
        json!({"workspace":source_ref,"design":unfinished,"first":0,"limit":2}),
    )
    .await;
    assert_eq!(page["verified"], false);
    let mut missing = design.clone();
    missing["bins"][0]["choices"][0]
        .as_object_mut()
        .unwrap()
        .remove("anneal_bottom");
    reject(
        &client,
        "inspect_combinatorial_construct",
        json!({"workspace":source_ref,"design":missing,"construct":construct}),
    )
    .await;
    let preview = call(
        &client,
        "inspect_combinatorial_construct",
        json!({"workspace":source_ref,"design":design,"construct":construct}),
    )
    .await;
    assert_eq!(preview["verified"], false);
    let evidence: Value =
        serde_json::from_slice(&artifact(&client, &preview["artifact_id"]).await).unwrap();
    let expected = vec![fixture::event(&original, 0), fixture::event(&original, 1)];
    let bindings: Vec<lineage::InputEventBinding> =
        serde_json::from_value(evidence["input_bindings"].clone()).unwrap();
    assert_eq!(
        bindings.iter().map(|i| i.event.clone()).collect::<Vec<_>>(),
        expected
    );
    assert_ne!(bindings[1].event.event_id, unselected_bottom);
    let anneal = &evidence["recipe"]["Ligation"]["recipe"]["fragments"][0]["source"]["preparation"]
        ["AnnealedOligos"];
    assert_eq!(anneal["candidate_index"], 0);
    assert_eq!(anneal["bottom_input_id"], bindings[1].input_id);
    assert_eq!(
        call(&client, "list_documents", json!({})).await,
        before_documents
    );

    let evaluated = call(
        &client,
        "evaluate_combinatorial_batch",
        json!({"workspace":source_ref,"matrix":matrix}),
    )
    .await;
    assert_eq!(evaluated["receipt"]["rows"][0]["input_count"], 2);
    let fresh: batch::Record =
        serde_json::from_slice(&artifact(&client, &evaluated["artifact_id"]).await).unwrap();
    assert_eq!(
        fresh.rows[0]
            .input_bindings
            .iter()
            .map(|i| i.event.clone())
            .collect::<Vec<_>>(),
        expected
    );
    let released=call(&client,"release_combinatorial_outputs",json!({"workspace":workspace_ref(&evaluated),
        "request":{"batch_hash":evaluated["receipt"]["attachment_hash"],"rows":[{"construct_id":"pair-A",
        "options":{"name":"Released paired AACC","show_prepared":false,"save_primers":false,"allow_incomplete":false}}]}})).await;
    assert_eq!(released["receipt"]["verified"], false);
    let manifest: release::Manifest =
        serde_json::from_slice(&artifact(&client, &released["artifact_id"]).await).unwrap();
    let native =
        workspace::load_workspace(&artifact(&client, &released["native_artifact_id"]).await)
            .unwrap();
    let row = &manifest.rows[0];
    assert_eq!(
        row.inputs
            .iter()
            .map(|i| i.event.clone())
            .collect::<Vec<_>>(),
        expected
    );
    assert_eq!(
        native.revisions[&row.output.construct.revision]
            .sequence
            .as_str(),
        "AACC"
    );
    let run = lineage::load_run(&native, &row.output.run.run_id).unwrap();
    assert_eq!(
        run.inputs
            .iter()
            .map(|i| i.event.clone())
            .collect::<Vec<_>>(),
        expected
    );
    assert!(run.plan.steps.iter().any(|op| matches!(
        op,
        Operation::AnnealOligos {
            candidate_index: 0,
            ..
        }
    )));
    for document in &original.documents {
        assert_eq!(
            native.head(&document.id).unwrap(),
            original.head(&document.id).unwrap()
        );
    }

    let outputs=call(&client,"export_combinatorial_outputs",json!({"workspace":workspace_ref(&released),"manifest_hash":released["receipt"]["manifest_hash"],
        "options":{"genbank":true,"fasta":false,"primer_csv":false,"primer_tsv":false,"protocol_html":true,"protocol_pdf":true,"construct_ids":["pair-A"]}})).await;
    let mut zip = zip::ZipArchive::new(Cursor::new(
        artifact(&client, &outputs["artifact_id"]).await,
    ))
    .unwrap();
    let (mut genbank, mut html, mut pdf) = (false, false, false);
    for index in 0..zip.len() {
        let mut entry = zip.by_index(index).unwrap();
        assert!(entry.enclosed_name().is_some());
        let path = entry.name().to_owned();
        let mut data = Vec::new();
        entry.read_to_end(&mut data).unwrap();
        if path.ends_with(".gb") {
            let decoded = collagene_core::formats::sequence::import_genbank(
                std::str::from_utf8(&data).unwrap(),
            )
            .unwrap();
            assert_eq!(decoded.molecules[0].sequence.as_str(), "AACC");
            genbank = true;
        } else if path.ends_with(".html") {
            let text = std::str::from_utf8(&data).unwrap();
            assert!(text.contains("Independent bottom"));
            assert!(text.contains("GGTT"));
            assert!(text.contains(&bindings[1].input_id));
            html = true;
        } else if path.ends_with(".pdf") {
            pdf_readback::check_xref(&data);
            let text = pdf_readback::extracted(&data)
                .split_whitespace()
                .collect::<Vec<_>>()
                .join(" ");
            assert!(text.contains("Independent bottom"));
            assert!(text.contains("Bottom ordered oligo"));
            assert!(text.contains("GGTT"));
            pdf = true;
        }
    }
    assert!(genbank && html && pdf);
    let source_export = call(&client, "export_workspace_archive", source_ref).await;
    let source_after =
        workspace::load_workspace(&artifact(&client, &source_export["artifact_id"]).await).unwrap();
    assert_eq!(serde_json::to_value(source_after).unwrap(), original_value);
    client.cancel().await.unwrap();
    server.await.unwrap();
}
