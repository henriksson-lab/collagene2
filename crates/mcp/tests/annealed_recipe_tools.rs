//! Three exact originals produce two prepared fragments without replacing either oligo input.
use base64::{Engine, engine::general_purpose::STANDARD};
use collagene_core::{
    chemistry::{AnnealOptions, Terminus},
    verification::DesignPlan,
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
async fn artifact(c: &Client, id: &Value) -> String {
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
    text
}
async fn evidence(c: &Client, id: &Value) -> Value {
    serde_json::from_str(&artifact(c, id).await).unwrap()
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
#[tokio::test]
async fn pinned_backbone_and_two_actual_oligos_ligate_and_export_without_source_substitution() {
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
    // The independent top reference contains the 14bp opened backbone followed
    // by the 12bp top oligo. Eight paired insert bases leave literal AATT ends.
    let text=[("backbone","AAAAGAATTCTTTT","circular"),("top","AATTCCGATGCA","linear"),("bottom","AATTTGCATCGG","linear"),("target","AATTCTTTTAAAAGAATTCCGATGCA","circular"),("wrong","AATTCTTTTAAAAGAATTCCGATGCC","circular")].into_iter().map(|(n,d,t)|format!("LOCUS       {n:<20} {} bp DNA {t} SYN 01-JAN-2026\nFEATURES             Location/Qualifiers\nORIGIN\n        1 {}\n//\n",d.len(),d.to_ascii_lowercase())).collect::<String>();
    let imported = call(&client, "import_document", json!({"text":text})).await;
    let catalog = collagene_core::catalog_digest::builtin_catalog().unwrap();
    let id = catalog
        .edition()
        .enzymes
        .iter()
        .find(|e| e.name == "EcoRI")
        .unwrap()
        .id
        .clone();
    let selection = collagene_core::catalog_digest::select(&catalog, &[id]).unwrap();
    let mut options = AnnealOptions::default();
    options.min_paired_bases = 8;
    options.top_chemistry = Terminus::CLEAVED;
    options.bottom_chemistry = Terminus::CLEAVED;
    let recipe = json!({"name":"Pinned EcoRI backbone plus two original oligos","fragments":[
        {"source":{"generation":null,"role":"backbone","input_id":"backbone","preparation":{"CatalogDigest":{"selection":selection,"fragment_index":0,"repair":null}},"reverse":false},"treatment":null},
        {"source":{"generation":null,"role":"insert","input_id":"top","preparation":{"AnnealedOligos":{"bottom_input_id":"bottom","options":options,"candidate_index":0}},"reverse":false},"treatment":null}
    ],"circular":true,"allow_nicked":false});
    let inputs =
        ["backbone", "top", "bottom"].map(|id| json!({"id":id,"document":named(&imported,id)}));
    let args = json!({"inputs":inputs,"target":named(&imported,"target"),"recipe":recipe});
    let result = call(&client, "propose_ligation_assembly", args.clone()).await;
    assert_eq!(result["fragment_count"], 2);
    assert_eq!(
        result["verification"]["status"],
        "ConsistentUnderAssumptions"
    );
    assert_eq!(result["product_physical"]["top_length"], 26);
    assert_eq!(result["product_physical"]["nick_count"], 0);
    let full = evidence(&client, &result["artifact_id"]).await;
    let plan: DesignPlan = serde_json::from_value(full["compiled"]["plan"].clone()).unwrap();
    assert_eq!(plan.inputs.len(), 3);
    assert_eq!(plan.inputs[1].molecule.sequence.as_str(), "AATTCCGATGCA");
    assert_eq!(plan.inputs[2].molecule.sequence.as_str(), "AATTTGCATCGG");
    assert!(plan.steps.iter().any(|s|matches!(s,collagene_core::verification::Operation::AnnealOligos{top:collagene_core::verification::Source::Input(a),bottom:collagene_core::verification::Source::Input(b),..}if a=="top"&&b=="bottom")));
    let reference = json!({"plan_revision":result["plan_revision"]});
    assert_eq!(
        call(&client, "verify_plan", reference.clone()).await["status"],
        "ConsistentUnderAssumptions"
    );
    let html = call(&client, "export_protocol", reference.clone()).await;
    let text = artifact(&client, &html["artifact_id"]).await;
    assert!(text.contains("Anneal declared"));
    assert!(text.contains("Top ordered oligo"));
    assert!(text.contains("Bottom ordered oligo"));
    assert!(text.contains("AATTCCGATGCA"));
    assert!(text.contains("AATTTGCATCGG"));
    let pdf = call(&client, "export_protocol_pdf", reference).await;
    let bytes = STANDARD
        .decode(artifact(&client, &pdf["artifact_id"]).await)
        .unwrap();
    assert!(bytes.starts_with(b"%PDF-"));
    assert!(bytes.ends_with(b"%%EOF") || bytes.ends_with(b"%%EOF\n"));
    assert_eq!(bytes.len() as u64, pdf["bytes"].as_u64().unwrap());
    let mut stale = args.clone();
    stale["inputs"][2]["document"]["revision"] = json!("stale");
    reject(&client, "propose_ligation_assembly", stale).await;
    let mut candidate = args.clone();
    candidate["recipe"]["fragments"][1]["source"]["preparation"]["AnnealedOligos"]["candidate_index"] =
        json!(999999);
    reject(&client, "propose_ligation_assembly", candidate).await;
    let mut bounded = args.clone();
    bounded["recipe"]["fragments"][1]["source"]["preparation"]["AnnealedOligos"]["options"]["max_candidates"] =
        json!(1000000);
    reject(&client, "propose_ligation_assembly", bounded).await;
    let mut wrong = args;
    wrong["target"] = named(&imported, "wrong");
    assert_eq!(
        call(&client, "propose_ligation_assembly", wrong).await["verification"]["status"],
        "Invalid"
    );
    assert_eq!(call(&client, "list_documents", json!({})).await["total"], 5);
    drop(client);
    server.abort();
}
