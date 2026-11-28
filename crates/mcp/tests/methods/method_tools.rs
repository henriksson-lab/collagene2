//! SDK transports do not confer end chemistry or empirical profile compatibility.
use collagene_core::formats;
use collagene_mcp::CollageneServer;
use collagene_model::{Molecule, Topology};
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
    let r = c
        .call_tool(
            CallToolRequestParams::new(name.to_owned())
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
            CallToolRequestParams::new(name.to_owned())
                .with_arguments(args.as_object().unwrap().clone()),
        )
        .await;
    assert!(
        r.is_err() || r.unwrap().is_error == Some(true),
        "accepted{name}"
    );
}
async fn import(c: &Client, name: &str, s: &str, t: Topology) -> Value {
    let m = Molecule::dna(name, s, t).unwrap();
    let r = call(
        c,
        "import_document",
        json!({"text":formats::export_genbank(&m).unwrap()}),
    )
    .await;
    json!({"id":r["documents"][0]["id"],"revision":r["documents"][0]["revision"]})
}
async fn artifact(c: &Client, id: &Value) -> Value {
    let mut at = 0;
    let mut text = String::new();
    loop {
        let r = call(
            c,
            "read_artifact",
            json!({"artifact_id":id,"start":at,"max_bytes":16384}),
        )
        .await;
        text.push_str(r["text"].as_str().unwrap());
        if r["done"] == true {
            break;
        }
        at = r["next_start"].as_u64().unwrap();
    }
    serde_json::from_str(&text).unwrap()
}
#[tokio::test]
async fn tailing_reconstructs_digest_and_unique_pcr_with_explicit_strands_and_no_import_authority()
{
    let (c, s) = setup().await;
    let source = import(&c, "digest source", "GAATTCAAAA", Topology::Circular).await;
    let r=call(&c,"simulate_tailing",json!({"document":source,"recipe":{"Digest":{"enzymes":["EcoRI"],"fragment_index":0,"blunt_before":true}},"base":"T"})).await;
    assert_eq!(r["top_length"], 15);
    assert_eq!(r["bottom_start"], -1);
    assert_eq!(r["additions"].as_array().unwrap().len(), 2);
    let full = artifact(&c, &r["artifact_id"]).await;
    assert_eq!(
        full["tailing"]["tailing"]["duplex"]["top"]["bases"],
        "AATTCAAAAGAATTT"
    );
    assert_eq!(full["tailing"]["source_revision"], source["revision"]);
    reject(&c,"simulate_tailing",json!({"document":source,"recipe":{"Digest":{"enzymes":["EcoRI"],"fragment_index":0,"blunt_before":false}},"base":"T"})).await;
    let pcr = import(&c, "PCR source", "ACGACCCCGGGGTAAC", Topology::Linear).await;
    let primer = |name: &str, seq: &str| json!({"name":name,"annealing":{"alphabet":"Dna","bases":seq},"tail":""});
    let recipe = json!({"ExactPcr":{"template_strand":"Top","forward":primer("f","ACGA"),"reverse":primer("r","GTTA")}});
    let r = call(
        &c,
        "simulate_tailing",
        json!({"document":pcr,"recipe":recipe,"base":"A"}),
    )
    .await;
    let full = artifact(&c, &r["artifact_id"]).await;
    assert_eq!(
        full["tailing"]["tailing"]["duplex"]["left"]["five_prime"],
        "Hydroxyl"
    );
    let mut bad = recipe;
    bad["ExactPcr"]["forward"]["five_prime"] = json!("Phosphate");
    reject(
        &c,
        "simulate_tailing",
        json!({"document":pcr,"recipe":bad,"base":"A"}),
    )
    .await;
    c.cancel().await.unwrap();
    s.await.unwrap();
}
#[tokio::test]
async fn oriented_profile_recommendations_are_advisory_but_target_reconstruction_is_independent() {
    let (c, s) = setup().await;
    assert_eq!(
        c.list_tools(Default::default()).await.unwrap().tools.len(),
        collagene_mcp::CURRENT_TOOL_COUNT
    );
    let a = import(&c, "a", "AAAACC", Topology::Linear).await;
    let b = import(&c, "b", "AACCGG", Topology::Linear).await;
    let target = import(&c, "target", "AAAACCGGTT", Topology::Linear).await;
    let options = json!({"method":"HiFi","scan_repeated_overlaps":false,"max_comparisons":10000,"max_occurrences":100});
    let report=call(&c,"check_overlap_profile",json!({"documents":[a,b],"reverse":[false,true],"overlaps":[2],"circular":false,"options":options})).await;
    assert_eq!(report["geometry_valid"], true);
    assert_eq!(report["recommendation_compliant"], false);
    assert_eq!(report["exhaustive_alternative_search"], false);
    let full = artifact(&c, &report["artifact_id"]).await;
    assert_eq!(
        full["overlap_profile"]["orientations"],
        json!([false, true])
    );
    assert_eq!(
        full["overlap_profile"]["source_revisions"],
        json!([a["revision"], b["revision"]])
    );
    let step = json!({"ProfiledOverlap":{"fragments":[{"source":{"Input":"a"},"reverse":false},{"source":{"Input":"b"},"reverse":true}],"overlaps":[2],"circular":false,"options":options}});
    let p=call(&c,"propose_plan",json!({"name":"oriented profile","target":target,"inputs":[{"id":"a","document":a},{"id":"b","document":b}],"steps":[step]})).await;
    let r = call(
        &c,
        "verify_plan",
        json!({"plan_revision":p["plan_revision"]}),
    )
    .await;
    assert_eq!(r["status"], "ConsistentUnderAssumptions");
    let html = call(
        &c,
        "export_protocol",
        json!({"plan_revision":p["plan_revision"]}),
    )
    .await;
    assert!(html.to_string().contains("artifact"));
    let mut bad = step;
    bad["ProfiledOverlap"]["fragments"][0]["phosphorylate"] = json!(true);
    reject(&c,"propose_plan",json!({"name":"ignored chemistry","target":target,"inputs":[{"id":"a","document":a},{"id":"b","document":b}],"steps":[bad]})).await;
    c.cancel().await.unwrap();
    s.await.unwrap();
}
