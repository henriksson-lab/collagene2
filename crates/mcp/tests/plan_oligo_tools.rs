//! Exact plan order output is not PCR or independent target verification.
use collagene_core::{analysis::Primer, mutagenesis::PcrOptions};
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
async fn call(c: &Client, name: &str, args: Value) -> Value {
    let result = c
        .call_tool(
            CallToolRequestParams::new(name.to_owned())
                .with_arguments(args.as_object().unwrap().clone()),
        )
        .await
        .unwrap();
    assert_ne!(result.is_error, Some(true), "{name}: {:?}", result.content);
    result.structured_content.unwrap()
}
async fn reject(c: &Client, name: &str, args: Value) {
    let result = c
        .call_tool(
            CallToolRequestParams::new(name.to_owned())
                .with_arguments(args.as_object().unwrap().clone()),
        )
        .await;
    assert!(result.is_err() || result.unwrap().is_error == Some(true));
}
async fn evidence(c: &Client, id: &Value) -> Value {
    let mut text = String::new();
    let mut start = 0;
    loop {
        let page = call(
            c,
            "read_artifact",
            json!({"artifact_id":id,"start":start,"max_bytes":16384}),
        )
        .await;
        text.push_str(page["text"].as_str().unwrap());
        if page["done"] == true {
            break;
        }
        start = page["next_start"].as_u64().unwrap();
    }
    serde_json::from_str(&text).unwrap()
}
fn named(imported: &Value, name: &str) -> Value {
    let row = imported["documents"]
        .as_array()
        .unwrap()
        .iter()
        .find(|d| d["name"] == name)
        .unwrap();
    json!({"id":row["id"],"revision":row["revision"]})
}
async fn text_artifact(client: &Client, id: &Value) -> String {
    let mut text = String::new();
    let mut start = 0;
    loop {
        let page = call(
            client,
            "read_artifact",
            json!({"artifact_id":id,"start":start,"max_bytes":16384}),
        )
        .await;
        text.push_str(page["text"].as_str().unwrap());
        if page["done"] == true {
            break;
        }
        start = page["next_start"].as_u64().unwrap();
    }
    text
}
#[tokio::test]
async fn actual_plan_oligos_keep_tails_duplicates_domain_uses_and_target_context_without_verifying()
{
    let (client, server) = setup().await;
    assert_eq!(
        client
            .list_tools(Default::default())
            .await
            .unwrap()
            .tools
            .len(),
        collagene_mcp::CURRENT_TOOL_COUNT
    );
    let imported=call(&client,"import_document",json!({"text":">source\nACGTACGTACGTACGT\n>independent_target\nTTTTCCCC\n>other_target\nAAAAGGGG\n"})).await;
    let source = named(&imported, "source");
    let target = named(&imported, "independent_target");
    let steps = json!([
        {"Pcr":{"template":{"Input":"s"},"forward":Primer::new("=forward","ACGT","GG").unwrap(),"reverse":Primer::new("same name","TTAA","").unwrap()}},
        {"PcrWithMismatches":{"template":{"Step":0},"forward":Primer::new("renamed","GACGT","G").unwrap(),"reverse":Primer::new("same name","CCGG","").unwrap(),"options":PcrOptions::default()}},
        {"TreatEnds":{"template":{"Step":1},"treatment":"Phosphorylate"}}
    ]);
    let proposal=call(&client,"propose_plan",json!({"name":"Actual declared orders","target":target,"inputs":[{"id":"s","document":source}],"steps":steps})).await;
    let args = json!({"plan_revision":proposal["plan_revision"],"format":"Csv"});
    let exported = call(&client, "export_plan_oligos", args.clone()).await;
    assert_eq!(exported["verified"], false);
    assert!(exported.get("verification_status").is_none());
    assert_eq!(exported["complete"], true);
    assert_eq!(exported["oligo_count"], 3);
    assert_eq!(exported["use_count"], 4);
    assert_eq!(exported["target_revision"], target["revision"]);
    let sheet = text_artifact(&client, &exported["table_artifact_id"]).await;
    assert!(sheet.contains("'=forward"));
    assert!(sheet.contains("GGACGT"));
    let report = evidence(&client, &exported["report_artifact_id"]).await;
    assert_eq!(report["uses"][0]["oligo_id"], report["uses"][2]["oligo_id"]);
    assert_ne!(report["uses"][1]["oligo_id"], report["uses"][3]["oligo_id"]);
    assert_eq!(report["uses"][0]["primer"]["tail"], "GG");
    assert_eq!(report["uses"][2]["primer"]["tail"], "G");
    assert!(report["uses"][2]["mismatch_options"].is_object());
    assert!(
        report["oligos"]["records"]
            .as_array()
            .unwrap()
            .iter()
            .all(|row| row["modification"].is_null() && row["quantity"].is_null())
    );
    let tsv=call(&client,"export_plan_oligos",json!({"plan_revision":proposal["plan_revision"],"format":"Tsv","options":{"deduplicate":false}})).await;
    assert_eq!(tsv["oligo_count"], 4);
    assert_ne!(tsv["context_revision"], exported["context_revision"]);
    assert!(
        text_artifact(&client, &tsv["table_artifact_id"])
            .await
            .contains('\t')
    );
    let other = named(&imported, "other_target");
    let second=call(&client,"propose_plan",json!({"name":"Same declared orders, other target","target":other,"inputs":[{"id":"s","document":source}],"steps":steps})).await;
    let second_export = call(
        &client,
        "export_plan_oligos",
        json!({"plan_revision":second["plan_revision"],"format":"Csv"}),
    )
    .await;
    assert_eq!(second_export["verified"], false);
    assert_eq!(second_export["target_revision"], other["revision"]);
    assert_ne!(
        second_export["context_revision"],
        exported["context_revision"]
    );
    assert_ne!(
        second_export["report_artifact_id"],
        exported["report_artifact_id"]
    );
    assert_eq!(
        second_export["table_artifact_id"],
        exported["table_artifact_id"]
    );
    assert_eq!(
        call(&client, "export_plan_oligos", args).await["report_artifact_id"],
        exported["report_artifact_id"]
    );
    reject(
        &client,
        "export_plan_oligos",
        json!({"plan_revision":"stale","format":"Csv"}),
    )
    .await;
    assert_eq!(call(&client, "list_documents", json!({})).await["total"], 3);
    drop(client);
    server.abort();
}
