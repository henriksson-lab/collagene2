//! MCP previews physical strands explicitly and attests saved material before reuse.
use collagene_core::{
    Topology,
    assembly_recipe::{self, FragmentSpec, Preparation},
    chemistry::Terminus,
    homology_annealing::Options,
    three_prime_resection::Cuts,
};
use collagene_mcp::CollageneServer;
use rmcp::{RoleClient, ServiceExt, model::CallToolRequestParams, service::RunningService};
use serde_json::{Value, json};
type Client = RunningService<RoleClient, ()>;
const A: &str = "AGTCCGATGCTA";
const B: &str = "GCTATTCGAATTCAGTC";
const CIRCLE: &str = "AGTCCGATGCTATTCGAATTC";
const CUT_FILLED: &str = "AATTCAGTCCGATGCTATTCGAATT";
fn policy() -> Options {
    Options {
        overlaps: vec![4, 4],
        topology: Topology::Circular,
        cuts: vec![Cuts { left: 4, right: 4 }; 2],
        min_paired_bases: 4,
        max_comparisons: 100_000,
        max_candidates: 100,
    }
}
fn recipe() -> assembly_recipe::Recipe {
    assembly_recipe::Recipe {
        name: "Explicit physical circular assembly".into(),
        fragments: ["a", "b"]
            .into_iter()
            .map(|id| FragmentSpec {
                generation: None,
                role: id.into(),
                input_id: id.into(),
                preparation: Preparation::DeclaredDuplex {
                    left: Terminus::PCR,
                    right: Terminus::PCR,
                },
                reverse: false,
            })
            .collect(),
        overlaps: vec![4, 4],
        circular: true,
        options: Default::default(),
    }
}
async fn call(c: &Client, name: &str, args: Value) -> Value {
    let reply = c
        .call_tool(
            CallToolRequestParams::new(name.to_owned())
                .with_arguments(args.as_object().unwrap().clone()),
        )
        .await
        .unwrap();
    assert_ne!(reply.is_error, Some(true), "{name}: {:?}", reply.content);
    reply.structured_content.unwrap()
}
async fn reject(c: &Client, name: &str, args: Value) {
    let reply = c
        .call_tool(
            CallToolRequestParams::new(name.to_owned())
                .with_arguments(args.as_object().unwrap().clone()),
        )
        .await;
    assert!(reply.is_err() || reply.unwrap().is_error == Some(true));
}
async fn artifact(c: &Client, id: &Value) -> Value {
    let (mut start, mut text) = (0, String::new());
    loop {
        let page = call(
            c,
            "read_artifact",
            json!({"artifact_id":id,"start":start,"max_bytes":16384}),
        )
        .await;
        text.push_str(page["text"].as_str().unwrap());
        if page["done"] == true {
            return serde_json::from_str(&text).unwrap();
        }
        start = page["next_start"].as_u64().unwrap();
    }
}
fn named(imported: &Value, name: &str) -> Value {
    let row = imported["documents"]
        .as_array()
        .unwrap()
        .iter()
        .find(|row| row["name"] == name)
        .unwrap();
    json!({"id":row["id"],"revision":row["revision"]})
}
#[tokio::test]
async fn annealing_preview_retains_incomplete_strands_and_never_saves_or_grants_material() {
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
    let client = ().serve(b).await.unwrap();
    let tools = client.list_all_tools().await.unwrap();
    assert_eq!(tools.len(), collagene_mcp::CURRENT_TOOL_COUNT);
    assert!(
        tools
            .iter()
            .any(|tool| tool.name == "preview_annealing_recipe")
    );
    let text=[("a",A,"linear"),("b",B,"linear"),("target",CIRCLE,"circular"),("pcr_target",CIRCLE,"linear"),("cut_target",CUT_FILLED,"linear")].into_iter().map(|(name,dna,topology)|format!("LOCUS       {name:<20} {} bp DNA {topology} SYN 01-JAN-2026\nFEATURES             Location/Qualifiers\nORIGIN\n        1 {}\n//\n",dna.len(),dna.to_ascii_lowercase())).collect::<String>();
    let imported = call(&client, "import_document", json!({"text":text})).await;
    let inputs = json!([{"id":"a","document":named(&imported,"a")},{"id":"b","document":named(&imported,"b")}]);
    let target = named(&imported, "target");
    let args = json!({"inputs":inputs,"target":target,"recipe":recipe(),"policy":policy()});
    let before = call(&client, "list_documents", json!({})).await;
    let preview = call(&client, "preview_annealing_recipe", args.clone()).await;
    assert_eq!(preview["state"], "preview");
    assert_eq!(preview["candidate_saved"], false);
    assert_eq!(preview["verified"], false);
    assert_eq!(preview["status"], "Incomplete", "{preview}");
    assert_eq!(preview["product"]["length"], 21);
    assert_eq!(preview["physical"]["paired_diagnostic_available"], true);
    assert_eq!(preview["physical"]["reusable_material_available"], false);
    assert_eq!(preview["physical"]["incomplete"], true);
    assert_eq!(call(&client, "list_documents", json!({})).await, before);
    let mut inaccessible = recipe();
    for row in &mut inaccessible.fragments {
        row.preparation = Preparation::DeclaredDuplex {
            left: Terminus::UNKNOWN,
            right: Terminus::UNKNOWN,
        };
    }
    let unavailable = call(
        &client,
        "preview_annealing_recipe",
        json!({"inputs":inputs,"target":target,"recipe":inaccessible,"policy":policy()}),
    )
    .await;
    assert_eq!(unavailable["status"], "Incomplete");
    assert!(
        unavailable["findings"]
            .as_array()
            .unwrap()
            .iter()
            .any(|finding| finding
                .as_str()
                .is_some_and(|text| text.contains("3′") || text.contains("3-prime")))
    );
    assert!(unavailable["product"].is_null());
    assert!(unavailable["physical"].is_null());
    assert_eq!(unavailable["candidate_saved"], false);
    let evidence = artifact(&client, &preview["artifact_id"]).await;
    let native: assembly_recipe::Preview =
        serde_json::from_value(evidence["preview"].clone()).unwrap();
    native
        .validate_binding(preview["context_revision"].as_str().unwrap())
        .unwrap();
    assert_eq!(
        native
            .verification
            .products
            .last()
            .unwrap()
            .sequence
            .as_str(),
        CIRCLE
    );
    let proposed=call(&client,"propose_plan",json!({"name":"Keep explicitly previewed physical homology candidate","inputs":inputs,"target":target,"steps":native.compiled.plan.steps})).await;
    assert_eq!(proposed["verified"], false);
    let verified = call(
        &client,
        "verify_plan",
        json!({"plan_revision":proposed["plan_revision"]}),
    )
    .await;
    assert_eq!(verified["status"], "Incomplete");
    let mut changed = args;
    changed["policy"]["overlaps"][0] = json!(3);
    reject(&client, "preview_annealing_recipe", changed).await;
    drop(client);
    server.await.unwrap();
}
