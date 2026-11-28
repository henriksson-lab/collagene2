//! Real MCP transport: proposals remain advisory until original-source preparation.
use collagene_core::planning::plan::{DesignPlan, Operation, Source};
use collagene_mcp::server::CollageneServer;
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

fn row(role: &str, input: &str, right: &str) -> Value {
    json!({"role":role,"input_id":input,"start":0,"end":16,"reverse":false,"forward_length":8,"reverse_length":8,"right_overhang":right})
}
fn ligation_recipe(fragments: &Value) -> Value {
    json!({"name":"Two PCR-added Type IIS substrates","fragments":fragments.as_array().unwrap().iter().map(|f|json!({"source":f,"treatment":null})).collect::<Vec<_>>(),"circular":false,"allow_nicked":false})
}

#[tokio::test]
async fn site_bearing_proposals_preserve_originals_and_need_fresh_independent_ligation_check() {
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
    // Retained top strands are independently specified full16nt cores. The
    // first right fusion TATG is the next top's first4 bases, not target-derived.
    let imported=call(&client,"import_document",json!({"text":">left\nACGACCCCGGGGTAAC\n>right\nTATGAGCTCCAGGTCA\n>independent_target\nACGACCCCGGGGTAACTATGAGCTCCAGGTCA\n>wrong_target\nACGACCCCGGGGTAACTATGAGCTCCAGGTCT\n"})).await;
    let left = named(&imported, "left");
    let right = named(&imported, "right");
    let before_left = call(&client, "export_document", left.clone()).await;
    let before_right = call(&client, "export_document", right.clone()).await;
    let inputs = json!([{"id":"l","document":left},{"id":"r","document":right}]);
    let request = json!({"name":"Two explicit BsaI primer substrates","fragments":[row("left","l","TATG"),row("right","r","ACGT")],"enzyme":"BsaI","clamp":"ACGT","spacer":"A"});
    let args =
        json!({"inputs":inputs,"target":named(&imported,"independent_target"),"request":request});
    let designed = call(&client, "design_golden_gate_primers", args.clone()).await;
    assert_eq!(designed["state"], "advisory");
    assert_eq!(designed["verified"], false);
    assert!(designed.get("plan_revision").is_none());
    let fragments = designed["fragments"].as_array().unwrap();
    assert_eq!(fragments.len(), 2);
    let left_prep = &fragments[0]["preparation"]["PcrCatalogDigest"];
    let right_prep = &fragments[1]["preparation"]["PcrCatalogDigest"];
    assert_eq!(left_prep["forward"]["annealing"]["bases"], "ACGACCCC");
    assert_eq!(left_prep["forward"]["tail"], "ACGTGGTCTCA");
    assert_eq!(left_prep["reverse"]["annealing"]["bases"], "GTTACCCC");
    assert_eq!(left_prep["reverse"]["tail"], "ACGTGGTCTCACATA");
    assert_eq!(right_prep["forward"]["annealing"]["bases"], "TATGAGCT");
    assert_eq!(right_prep["reverse"]["annealing"]["bases"], "TGACCTGG");
    assert_eq!(right_prep["reverse"]["tail"], "ACGTGGTCTCAACGT");
    assert_eq!(left_prep["repair"], Value::Null);
    assert_eq!(right_prep["repair"], Value::Null);
    let report = evidence(&client, &designed["artifact_id"]).await;
    assert_eq!(report["inputs"], inputs);
    assert_eq!(report["target"], named(&imported, "independent_target"));
    assert_eq!(report["report"]["fragments"], designed["fragments"]);
    let repeat = call(&client, "design_golden_gate_primers", args.clone()).await;
    assert_eq!(repeat["context_revision"], designed["context_revision"]);
    assert_eq!(repeat["artifact_id"], designed["artifact_id"]);
    let proposal=call(&client,"propose_ligation_assembly",json!({"inputs":inputs,"target":named(&imported,"independent_target"),"recipe":ligation_recipe(&designed["fragments"])})).await;
    assert_eq!(proposal["verified"], false);
    assert_eq!(proposal["preparation_complete"], true);
    assert_eq!(
        proposal["verification"]["status"],
        "ConsistentUnderAssumptions"
    );
    assert_eq!(proposal["product_physical"]["top_length"], 32);
    assert_eq!(proposal["product_physical"]["nick_count"], 0);
    let full = evidence(&client, &proposal["artifact_id"]).await;
    let plan: DesignPlan = serde_json::from_value(full["compiled"]["plan"].clone()).unwrap();
    assert_eq!(plan.inputs.len(), 2);
    assert_eq!(plan.steps.len(), 5);
    assert_eq!(
        plan.inputs[0].molecule.sequence.as_str(),
        "ACGACCCCGGGGTAAC"
    );
    assert_eq!(
        plan.inputs[1].molecule.sequence.as_str(),
        "TATGAGCTCCAGGTCA"
    );
    assert_eq!(
        plan.inputs[0].revision,
        named(&imported, "left")["revision"].as_str().unwrap()
    );
    assert_eq!(
        plan.inputs[1].revision,
        named(&imported, "right")["revision"].as_str().unwrap()
    );
    assert!(matches!(&plan.steps[0],Operation::Pcr{template:Source::Input(id),..} if id=="l"));
    assert!(matches!(
        &plan.steps[1],
        Operation::CatalogDigest {
            template: Source::Step(0),
            ..
        }
    ));
    assert!(matches!(&plan.steps[2],Operation::Pcr{template:Source::Input(id),..} if id=="r"));
    assert!(matches!(
        &plan.steps[3],
        Operation::CatalogDigest {
            template: Source::Step(2),
            ..
        }
    ));
    assert!(matches!(&plan.steps[4], Operation::Ligate { .. }));
    assert_eq!(
        full["verification"]["products"].as_array().unwrap().len(),
        5
    );
    let mut wrong = args.clone();
    wrong["target"] = named(&imported, "wrong_target");
    let wrong_design = call(&client, "design_golden_gate_primers", wrong).await;
    assert_eq!(wrong_design["verified"], false);
    assert_eq!(wrong_design["fragments"], designed["fragments"]);
    assert_ne!(
        wrong_design["context_revision"],
        designed["context_revision"]
    );
    let mismatch=call(&client,"propose_ligation_assembly",json!({"inputs":inputs,"target":named(&imported,"wrong_target"),"recipe":ligation_recipe(&wrong_design["fragments"])})).await;
    assert_eq!(mismatch["verification"]["status"], "Invalid");
    let mut stale = args.clone();
    stale["inputs"][0]["document"]["revision"] = json!("0".repeat(64));
    reject(&client, "design_golden_gate_primers", stale).await;
    let mut out_of_bounds = args.clone();
    out_of_bounds["request"]["fragments"][0]["end"] = json!(17);
    reject(&client, "design_golden_gate_primers", out_of_bounds).await;
    let mut wrong_overhang = args;
    wrong_overhang["request"]["fragments"][0]["right_overhang"] = json!("AAA");
    reject(&client, "design_golden_gate_primers", wrong_overhang).await;
    assert_eq!(
        call(&client, "export_document", named(&imported, "left")).await["artifact_id"],
        before_left["artifact_id"]
    );
    assert_eq!(
        call(&client, "export_document", named(&imported, "right")).await["artifact_id"],
        before_right["artifact_id"]
    );
    assert_eq!(call(&client, "list_documents", json!({})).await["total"], 4);
    drop(client);
    server.abort();
}
