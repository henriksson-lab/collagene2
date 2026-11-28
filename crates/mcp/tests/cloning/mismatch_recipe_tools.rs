//! Explicit mismatch policy crosses the actual MCP recipe, verifier and oligo handoff.
use collagene_core::{
    cloning::{ligation::recipe as ligation_recipe, pcr::mutagenesis::PcrOptions},
    planning::plan::Operation,
    sequence::analysis::Primer,
    verification::evidence::Status,
};
use collagene_mcp::server::CollageneServer;
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
async fn text(c: &Client, id: &Value) -> String {
    let (mut start, mut out) = (0, String::new());
    loop {
        let p = call(
            c,
            "read_artifact",
            json!({"artifact_id":id,"start":start,"max_bytes":16384}),
        )
        .await;
        out.push_str(p["text"].as_str().unwrap());
        if p["done"] == true {
            break;
        }
        start = p["next_start"].as_u64().unwrap();
    }
    out
}
fn named(imported: &Value, name: &str) -> Value {
    let d = imported["documents"]
        .as_array()
        .unwrap()
        .iter()
        .find(|d| d["name"] == name)
        .unwrap();
    json!({"id":d["id"],"revision":d["revision"]})
}
#[tokio::test]
async fn explicit_mismatch_recipe_keeps_policy_oligos_independent_target_and_immutable_source() {
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
    const ORIGINAL: &str = "AGTCCGATGCTACGTAGGACTTCAGACCTGTA";
    const MUTANT: &str = "AGACCGATGCTACGTAGGACTTCAGACCTGGA";
    let imported=call(&client,"import_document",json!({"text":format!(">original\n{ORIGINAL}\n>independent_mutant\n{MUTANT}\n>wrong_requirement\n{ORIGINAL}\n")})).await;
    let original = named(&imported, "original");
    let target = named(&imported, "independent_mutant");
    let mut policy = PcrOptions::default();
    policy.binding.max_mismatches = 1;
    policy.binding.min_three_prime_exact = 2;
    let forward = Primer::new("Forward mutation at source position3", "AGACCGAT", "").unwrap();
    let reverse = Primer::new("Reverse mutation at source position31", "TCCAGGTC", "").unwrap();
    let args = json!({"inputs":[{"id":"source","document":original}],"target":target,"recipe":{"standard": null,"name":"Explicit mismatch primers, independently required mutant","fragments":[{"source":{"generation":null,"role":"mutant","input_id":"source","preparation":{"MismatchPcr":{"forward":forward,"reverse":reverse,"options":policy}} ,"reverse":false},"treatment":null}],"circular":false,"allow_nicked":false}});
    let proposed = call(&client, "propose_ligation_assembly", args.clone()).await;
    assert_eq!(proposed["verified"], false);
    assert_eq!(
        proposed["verification"]["status"],
        "ConsistentUnderAssumptions"
    );
    let mut artifact: Value =
        serde_json::from_str(&text(&client, &proposed["artifact_id"]).await).unwrap();
    // MCP adds separately typed row diagnostics to the unchanged preview fields.
    let diagnostics: collagene_core::cloning::assembly::diagnostics::Report =
        serde_json::from_value(
            artifact
                .as_object_mut()
                .unwrap()
                .remove("diagnostics")
                .unwrap(),
        )
        .unwrap();
    let preview: ligation_recipe::Preview = serde_json::from_value(artifact).unwrap();
    assert_eq!(diagnostics.context_revision, preview.context_revision);
    preview.validate().unwrap();
    assert_eq!(
        preview.verification.status,
        Status::ConsistentUnderAssumptions
    );
    assert_eq!(
        preview.compiled.plan.inputs[0].molecule.sequence.as_str(),
        ORIGINAL
    );
    assert_eq!(
        preview
            .verification
            .final_product()
            .expect("final declared step must reconstruct")
            .sequence
            .as_str(),
        MUTANT
    );
    let Some(Operation::PcrWithMismatches {
        forward: actual_f,
        reverse: actual_r,
        options,
        ..
    }) = preview.compiled.plan.steps.first()
    else {
        panic!("actual mismatch-PCR operation retained")
    };
    assert_eq!(actual_f.ordered_sequence(), "AGACCGAT");
    assert_eq!(actual_r.ordered_sequence(), "TCCAGGTC");
    assert_eq!(options.binding.max_mismatches, 1);
    assert_eq!(options.binding.min_three_prime_exact, 2);
    assert_eq!(
        serde_json::to_value(options).unwrap(),
        serde_json::to_value(policy).unwrap()
    );
    let reference = json!({"plan_revision":proposed["plan_revision"]});
    assert_eq!(
        call(&client, "verify_plan", reference.clone()).await["status"],
        "ConsistentUnderAssumptions"
    );
    let protocol = call(&client, "export_protocol", reference.clone()).await;
    let html = text(&client, &protocol["artifact_id"]).await;
    assert!(html.contains("Intentional-mismatch PCR"));
    assert!(html.contains("AGACCGAT"));
    assert!(html.contains("TCCAGGTC"));
    assert!(html.contains("min_three_prime_exact"));
    let oligos = call(
        &client,
        "export_plan_oligos",
        json!({"plan_revision":proposed["plan_revision"],"format":"Csv"}),
    )
    .await;
    assert_eq!(oligos["oligo_count"], 2);
    let table = text(&client, &oligos["table_artifact_id"]).await;
    assert!(table.contains("AGACCGAT"));
    assert!(table.contains("TCCAGGTC"));
    let mut wrong = args;
    wrong["target"] = named(&imported, "wrong_requirement");
    let wrong = call(&client, "propose_ligation_assembly", wrong).await;
    assert_eq!(wrong["verification"]["status"], "Invalid");
    assert_eq!(
        call(
            &client,
            "verify_plan",
            json!({"plan_revision":wrong["plan_revision"]})
        )
        .await["status"],
        "Invalid"
    );
    let unchanged = call(
        &client,
        "read_sequence_region",
        json!({"document":original,"start":0,"end":32}),
    )
    .await;
    assert_eq!(unchanged["sequence"], ORIGINAL);
    assert_eq!(call(&client, "list_documents", json!({})).await["total"], 3);
    drop(client);
    server.abort();
}
