//! Real SDK boundary, independently fixed synthetic scars, native event ownership and stale selection.
use base64::{Engine as _, engine::general_purpose::STANDARD};
use collagene_core::{Molecule, Topology, workspace};
use collagene_mcp::CollageneServer;
use rmcp::{RoleClient, ServiceExt, model::CallToolRequestParams, service::RunningService};
use serde_json::{Value, json};
type Client = RunningService<RoleClient, ()>;
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
    assert!(r.is_err() || r.unwrap().is_error == Some(true));
}
async fn artifact(c: &Client, id: &Value) -> Value {
    let (mut text, mut offset) = (String::new(), 0);
    loop {
        let p = call(
            c,
            "read_artifact",
            json!({"artifact_id":id,"start":offset,"max_bytes":8192}),
        )
        .await;
        text.push_str(p["text"].as_str().unwrap());
        if p["done"] == true {
            break;
        }
        offset = p["next_start"].as_u64().unwrap();
    }
    serde_json::from_str(&text).unwrap()
}
fn request() -> Value {
    let mut sources = vec![];
    let mut variants = vec![];
    for (i, core) in [
        (1, "GTACAAA"),
        (2, "GTACAAG"),
        (3, "GTACAAC"),
        (4, "GTACAAT"),
    ] {
        for (family, left, right) in [("L", "ACGG", "TCC"), ("R", "TGGA", "GCT")] {
            let id = format!("synthetic-{family}{i}");
            let sequence = format!("{left}{core}{right}");
            sources.push(json!({"id":id,"url":"https://example.org/independent-synthetic-fixture","raw_sha256":"1".repeat(64),"vector_name":"Explicit synthetic reference bookkeeping only","vector_length":sequence.len(),"start_1based":1,"end_1based":sequence.len(),"reverse_extraction":false,"map_url":"https://example.org/independent-synthetic-map","map_sha256":"2".repeat(64)}));
            variants.push(json!({"id":id,"family":family,"subtype":i,"sequence":sequence,"overlap_start":4,"source_id":id}));
        }
    }
    json!({"reaction":{"profile":{"id":"sdk-literal-four-parent","edition":"1","sources":sources,"variants":variants,"limitations":["Synthetic reference fixture; no claimed biological profile or measured specificity"]},"reaction":"LR","specificities":[4,1,2,3],"limits":collagene_core::recombination::reference::Limits::default()},"annotation_limits":collagene_core::recombination::AnnotationLimits::default()})
}
// These scars and source body strings are declared before calling any simulator.
const MAIN: &str = concat!(
    "GTACAAT", "TCC", "CCCCA", "TGGA", "GTACAAA", "TCC", "GGGGT", "ACGG", "GTACAAG", "GCT",
    "AAAAG", "ACGG", "GTACAAC", "GCT", "GGGGA", "TGGA"
);
const PARENTS: [&str; 4] = [
    "ACGGGTACAATTCCCCCCATGGAGTACAAAGCTAAAAC",
    "ACGGGTACAAATCCGGGGTACGGGTACAAGTCCTTTTC",
    "TGGAGTACAAGGCTAAAAGACGGGTACAACTCCCCCCT",
    "TGGAGTACAATGCTTTTTATGGAGTACAACGCTGGGGA",
];
#[tokio::test]
async fn sdk_multisite_exact_events_projection_independent_target_and_full_plan_evidence() {
    let (a, b) = tokio::io::duplex(65536);
    let task = tokio::spawn(async move {
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
    for name in [
        "inspect_multisite_recombination",
        "select_multisite_recombination_component",
        "check_multisite_recombination_target",
    ] {
        let t = tools.iter().find(|t| t.name == name).unwrap();
        assert_eq!(t.input_schema["additionalProperties"], false);
    }
    let mut molecules = PARENTS
        .iter()
        .enumerate()
        .map(|(i, s)| Molecule::dna(format!("parent{i}"), *s, Topology::Circular).unwrap())
        .collect::<Vec<_>>();
    molecules.push(
        Molecule::dna(
            "Independent source-body scar intent",
            MAIN,
            Topology::Circular,
        )
        .unwrap(),
    );
    molecules.push(
        Molecule::dna(
            "Independent wrong requirement",
            "ACGTTT",
            Topology::Circular,
        )
        .unwrap(),
    );
    let frame = workspace::Workspace::from_documents(
        "Immutable four-parent reference fixture",
        molecules.clone(),
        0,
    )
    .unwrap();
    let imported = call(
        &client,
        "import_workspace_archive",
        json!({"data_base64":STANDARD.encode(workspace::save_workspace(&frame).unwrap())}),
    )
    .await;
    let wr = json!({"id":imported["workspace"]["id"],"revision":imported["workspace"]["revision"]});
    let doc = |i: usize| {
        let rev = molecules[i].revision();
        json!({"id":rev,"revision":rev})
    };
    let mut sources = vec![];
    for i in 0..4 {
        let d = &frame.documents[i];
        let e = frame
            .history
            .iter()
            .find(|e| e.document_id == d.id)
            .unwrap();
        sources.push(json!({"role":format!("parent-{i}"),"document":doc(i),"native":{"workspace":wr,"event":{"document_id":d.id,"event_id":e.id,"revision":e.revision}}}));
    }
    let request = request();
    let r = call(
        &client,
        "inspect_multisite_recombination",
        json!({"sources":sources,"request":request,"limit":1}),
    )
    .await;
    assert_eq!(r["complete"], true);
    assert_eq!(r["candidate_count"], 1);
    let batch = artifact(&client, &r["artifact_id"]).await;
    assert_eq!(batch["sources"], json!(sources));
    assert_eq!(batch["batch"]["documents"].as_array().unwrap().len(), 2);
    let lengths = r["candidates"][0]["components"].as_array().unwrap();
    assert_eq!(
        lengths
            .iter()
            .map(|p| p["length"].as_u64().unwrap())
            .sum::<u64>(),
        molecules[..4]
            .iter()
            .map(|m| m.sequence.len() as u64)
            .sum::<u64>()
    );
    // Both reciprocal circles are76bp. Length is not a component identity.
    // Select the declared source0 MCS/cargo [14,19), rather than choosing a
    // component by target matching or trusting its generated product name.
    let products = batch["batch"]["reference"]["candidates"][0]["products"]
        .as_array()
        .unwrap();
    let main_components = products
        .iter()
        .enumerate()
        .filter(|(_, p)| {
            p["segments"].as_array().unwrap().iter().any(|segment| {
                segment["input"] == 0
                    && segment["source_start"].as_u64().unwrap() <= 14
                    && segment["source_end"].as_u64().unwrap() >= 19
                    && segment["reverse"] == false
            })
        })
        .map(|(id, _)| id)
        .collect::<Vec<_>>();
    assert_eq!(main_components.len(), 1);
    let component_id = main_components[0];
    let selection =
        json!({"batch_revision":r["batch_revision"],"candidate_id":0,"component_id":component_id});
    let checked = call(
        &client,
        "check_multisite_recombination_target",
        json!({"sources":sources,"request":request,"selection":selection,"target":doc(4)}),
    )
    .await;
    assert_eq!(checked["assessment"]["target_matches"], true);
    assert_eq!(checked["assessment"]["complete"], true);
    let reciprocal_selection = json!({"batch_revision":r["batch_revision"],"candidate_id":0,"component_id":1-component_id});
    let reciprocal_check = call(&client,"check_multisite_recombination_target",json!({"sources":sources,"request":request,"selection":reciprocal_selection,"target":doc(4)})).await;
    assert_eq!(
        reciprocal_check["assessment"]["target_matches"], false,
        "an equal-length reciprocal component must not be replaced by the target-matching component"
    );
    let before = call(&client, "list_documents", json!({})).await;
    let wrong = call(
        &client,
        "check_multisite_recombination_target",
        json!({"sources":sources,"request":request,"selection":selection,"target":doc(5)}),
    )
    .await;
    assert_eq!(wrong["assessment"]["target_matches"], false);
    assert_eq!(call(&client, "list_documents", json!({})).await, before);
    let selected = call(
        &client,
        "select_multisite_recombination_component",
        json!({"sources":sources,"request":request,"selection":selection}),
    )
    .await;
    assert_eq!(selected["document"]["length"], MAIN.len());
    assert_eq!(selected["verified"], false);
    let evidence = artifact(&client, &selected["artifact_id"]).await;
    assert_eq!(
        evidence["product"]["reciprocal_components"]
            .as_array()
            .unwrap()
            .len(),
        2
    );
    assert!(
        evidence["product"]["molecule"]["metadata"]
            .get("multisite_recombination:unverified_recipe")
            .is_some()
    );
    let stepsources = (0..4)
        .map(|i| json!({"role":format!("parent-{i}"),"source":{"Input":format!("input-{i}")}}))
        .collect::<Vec<_>>();
    let inputs = (0..4)
        .map(|i| json!({"id":format!("input-{i}"),"document":doc(i)}))
        .collect::<Vec<_>>();
    let proposed=call(&client,"propose_plan",json!({"name":"Explicit four-parent reference network","target":doc(4),"inputs":inputs,"steps":[{"MultiSiteRecombineReference":{"sources":stepsources,"request":request,"selection":selection}}]})).await;
    let verified = call(
        &client,
        "verify_plan",
        json!({"plan_revision":proposed["plan_revision"]}),
    )
    .await;
    assert_eq!(verified["reference_step_count"], 1);
    assert_eq!(verified["status"], "Incomplete");
    let full = artifact(&client, &verified["artifact_id"]).await;
    assert_eq!(full["reference_step_reports"].as_array().unwrap().len(), 1);
    assert!(
        full["reference_step_reports"][0]["evidence"]
            .get("MultiSiteRecombination")
            .is_some()
    );
    let mut forged = sources.clone();
    forged[0]["native"]["event"]["document_id"] = json!(frame.documents[1].id);
    reject(
        &client,
        "inspect_multisite_recombination",
        json!({"sources":forged,"request":request}),
    )
    .await;
    let mut stale = selection.clone();
    stale["batch_revision"] = json!("0".repeat(64));
    reject(
        &client,
        "select_multisite_recombination_component",
        json!({"sources":sources,"request":request,"selection":stale}),
    )
    .await;
    reject(
        &client,
        "inspect_multisite_recombination",
        json!({"sources":sources}),
    )
    .await;
    let mut changed = request.clone();
    changed["reaction"]["profile"]["edition"] = json!("another-edition");
    reject(
        &client,
        "select_multisite_recombination_component",
        json!({"sources":sources,"request":changed,"selection":selection}),
    )
    .await;
    client.cancel().await.unwrap();
    task.await.unwrap();
}
