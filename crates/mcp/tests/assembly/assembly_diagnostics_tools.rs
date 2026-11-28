//! Additive typed locations distinguish preparation failure, blocked rows and target mismatch.
use collagene_core::sequence::analysis::Primer;
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

fn pcr(role: &str, input: &str, forward: &str, reverse: &str, tail: &str) -> Value {
    json!({"generation":null,"role":role,"input_id":input,"preparation":{"ExactPcr":{"template_strand":"Top","forward":Primer::new(format!("{role}-F"),forward,"").unwrap(),"reverse":Primer::new(format!("{role}-R"),reverse,tail).unwrap()}},"reverse":false})
}
fn profile() -> Value {
    json!({"method":"HiFi","scan_repeated_overlaps":false,"max_comparisons":100000,"max_occurrences":100})
}
async fn fixture(c: &Client) -> (Value, Value, Value) {
    let imported=call(c,"import_document",json!({"text":">A\nATGACCCCGGGGTTTCCAGA\n>B\nATCCGTTAGCGACCTTGGAA\n>C\nGCATAGTCGACTAGCTTGCA\n>Wrong\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n"})).await;
    let a = named(&imported, "A");
    let b = named(&imported, "B");
    let third = named(&imported, "C");
    let wrong = named(&imported, "Wrong");
    (
        json!([{"id":"a","document":a},{"id":"b","document":b}]),
        json!([{"id":"a","document":a},{"id":"b","document":b},{"id":"c","document":third}]),
        wrong,
    )
}
async fn matches_artifact(c: &Client, result: &Value) -> Value {
    let full = evidence(c, &result["artifact_id"]).await;
    assert_eq!(full["diagnostics"], result["diagnostics"]);
    assert_eq!(full["context_revision"], result["context_revision"]);
    assert!(full["compiled"].is_object());
    assert!(full["verification"].is_object());
    full
}
#[tokio::test]
async fn failed_middle_preparation_locates_the_executed_step_and_blocks_later_rows() {
    let (c, server) = setup().await;
    let (_, inputs, target) = fixture(&c).await;
    let rows = vec![
        pcr("first", "a", "ATGACCCC", "TCTGGAAA", ""),
        pcr("middle", "b", "AAAAAAAA", "CCCCCCCC", ""),
        pcr("last", "c", "GCATAGTC", "TGCAAGCT", ""),
    ];
    for ligation in [false, true] {
        let (tool, recipe, expected_step) = if ligation {
            (
                "propose_ligation_assembly",
                json!({"standard": null,"name":"failed middle physical preparation","fragments":rows.iter().map(|row|json!({"source":row,"treatment":"Phosphorylate"})).collect::<Vec<_>>(),"circular":false,"allow_nicked":false}),
                1,
            )
        } else {
            (
                "propose_overlap_assembly",
                json!({"standard": null,"name":"failed middle overlap preparation","fragments":rows,"overlaps":[15,15],"circular":false,"options":profile()}),
                1,
            )
        };
        let result = call(
            &c,
            tool,
            json!({"inputs":inputs,"target":target,"recipe":recipe}),
        )
        .await;
        assert_eq!(result["state"], "candidate");
        assert_eq!(result["verified"], false);
        assert_eq!(result["preparation_complete"], false);
        let d = &result["diagnostics"];
        assert_eq!(d["context_revision"], result["context_revision"]);
        assert_eq!(d["plan_revision"], result["plan_revision"]);
        assert_eq!(d["target_revision"], target["revision"]);
        // compile_preparations emits all PCR operations first; treatments follow.
        assert_eq!(d["rows"][0]["state"], "Reconstructed");
        assert_eq!(
            d["rows"][0]["completed_steps"],
            if ligation { 2 } else { 1 }
        );
        assert_eq!(d["rows"][0]["total_steps"], if ligation { 2 } else { 1 });
        assert_eq!(d["rows"][0]["role"], "first");
        assert_eq!(d["rows"][1]["state"], "Rejected");
        assert_eq!(d["rows"][1]["index"], 1);
        assert_eq!(d["rows"][1]["role"], "middle");
        assert_eq!(d["rows"][1]["input_id"], "b");
        assert_eq!(d["rows"][2]["state"], "Reconstructed");
        assert_eq!(
            d["rows"][2]["completed_steps"], d["rows"][2]["total_steps"],
            "independent later preparation must complete despite the rejected middle row"
        );
        assert_eq!(d["failures"][0]["row"], 1);
        assert_eq!(d["failures"][0]["step"], expected_step);
        assert_eq!(d["failures"][0]["stage"], "Pcr");
        assert!(
            d["failures"][0]["recovery"]
                .as_str()
                .is_some_and(|s| !s.is_empty())
        );
        let full = matches_artifact(&c, &result).await;
        assert_eq!(full["prepared"].as_array().unwrap().len(), 3);
        if ligation {
            assert_eq!(full["compiled"]["rows"][0]["treatment_step"], 3);
            assert!(full["compiled"]["plan"]["steps"][0].get("Pcr").is_some());
            assert!(full["compiled"]["plan"]["steps"][1].get("Pcr").is_some());
            assert!(
                full["compiled"]["plan"]["steps"][3]
                    .get("TreatEnds")
                    .is_some()
            );
        }
        assert_eq!(
            full["verification"]["steps"].as_array().unwrap().len(),
            full["compiled"]["plan"]["steps"].as_array().unwrap().len()
        );
        assert!(full["prepared"][1].is_null());
        assert!(full["prepared"][2].is_object());
        assert_eq!(
            full["verification"]["steps"][expected_step as usize]["state"],
            "Rejected"
        );
        assert!(full["verification"]["steps"][expected_step as usize]["product"].is_null());
        // Diagnostics do not change registry count, candidate identity or historical export routes.
        let repeated = call(
            &c,
            tool,
            json!({"inputs":inputs,"target":target,"recipe":recipe}),
        )
        .await;
        assert_eq!(repeated["artifact_id"], result["artifact_id"]);
        assert_eq!(repeated["plan_revision"], result["plan_revision"]);
    }
    assert_eq!(
        c.list_tools(Default::default()).await.unwrap().tools.len(),
        collagene_mcp::CURRENT_TOOL_COUNT
    );
    c.cancel().await.unwrap();
    server.await.unwrap();
}
#[tokio::test]
async fn final_target_mismatch_preserves_reconstructed_rows_and_real_junction_evidence() {
    let (c, server) = setup().await;
    let (inputs, _, target) = fixture(&c).await;
    for ligation in [false, true] {
        let left = pcr(
            "left",
            "a",
            "ATGACCCC",
            "TCTGGAAA",
            if ligation { "" } else { "AGGTCGCTAACGGAT" },
        );
        let right = pcr("right", "b", "ATCCGTTA", "TTCCAAGG", "");
        let (tool, recipe) = if ligation {
            (
                "propose_ligation_assembly",
                json!({"standard": null,"name":"wrong target is not failed PCR","fragments":[{"source":left,"treatment":"Phosphorylate"},{"source":right,"treatment":"Phosphorylate"}],"circular":false,"allow_nicked":false}),
            )
        } else {
            (
                "propose_overlap_assembly",
                json!({"standard": null,"name":"wrong target is not failed overlap preparation","fragments":[left,right],"overlaps":[15],"circular":false,"options":profile()}),
            )
        };
        let result = call(
            &c,
            tool,
            json!({"inputs":inputs,"target":target,"recipe":recipe}),
        )
        .await;
        assert_eq!(result["verified"], false);
        assert_eq!(result["preparation_complete"], true);
        assert_eq!(result["verification"]["status"], "Invalid");
        let d = &result["diagnostics"];
        assert_eq!(d["status"], "Invalid");
        assert_eq!(d["rows"][0]["state"], "Reconstructed");
        assert_eq!(d["rows"][1]["state"], "Reconstructed");
        assert_eq!(d["failures"][0]["stage"], "TargetComparison");
        assert!(d["failures"][0]["row"].is_null());
        assert_eq!(d["junctions"][0]["upstream"], 0);
        assert_eq!(d["junctions"][0]["downstream"], 1);
        assert_eq!(
            d["junctions"][0]["state"],
            if ligation {
                "EndpointsAvailable"
            } else {
                "ExactHomology"
            }
        );
        let full = matches_artifact(&c, &result).await;
        assert_eq!(full["prepared"].as_array().unwrap().len(), 2);
        assert_eq!(d["findings"], full["verification"]["findings"]);
        assert_eq!(d["plan_revision"], result["plan_revision"]);
    }
    c.cancel().await.unwrap();
    server.await.unwrap();
}
