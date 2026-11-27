use collagene_core::{
    Molecule, Topology,
    golden_gate::GoldenGateOptions,
    verification::{Operation, Source},
};
use collagene_mcp::CollageneServer;
use rmcp::{RoleClient, ServiceExt, model::CallToolRequestParams, service::RunningService};
use serde_json::{Value, json};
type Client = RunningService<RoleClient, ()>;
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
fn doc(v: &Value) -> Value {
    let d = &v["documents"][0];
    json!({"id":d["id"],"revision":d["revision"]})
}
async fn full(c: &Client, id: &Value) -> Value {
    let mut text = String::new();
    let mut start = 0;
    loop {
        let r = call(
            c,
            "read_artifact",
            json!({"artifact_id":id,"start":start,"max_bytes":16384}),
        )
        .await;
        text.push_str(r["text"].as_str().unwrap());
        if r["done"] == true {
            break;
        }
        start = r["next_start"].as_u64().unwrap();
    }
    serde_json::from_str(&text).unwrap()
}
#[tokio::test]
async fn actual_golden_gate_windows_reverse_bottom_coordinates_and_wrong_targets_share_the_existing_tool()
 {
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
    let c = ().serve(b).await.unwrap();
    assert_eq!(
        c.list_tools(Default::default()).await.unwrap().tools.len(),
        collagene_mcp::CURRENT_TOOL_COUNT
    );
    let source = Molecule::dna(
        "amplicon",
        "ATGCGGTCTCAATGACCCCGGGGTTTCCAGAATGATGAGACCGCAT",
        Topology::Linear,
    )
    .unwrap();
    let imported = call(
        &c,
        "import_document",
        json!({"text":collagene_core::formats::export_genbank(&source).unwrap()}),
    )
    .await;
    let source = doc(&imported);
    let before = call(&c, "export_document", source.clone()).await;
    for (reverse, wrong) in [(false, false), (true, false), (false, true)] {
        let target = Molecule::dna(
            format!("target-{reverse}-{wrong}"),
            if wrong {
                "AAAAAAAAAAAAAAAAAAAA"
            } else if reverse {
                "TCTGGAAACCCCGGGGTCAT"
            } else {
                "ATGACCCCGGGGTTTCCAGA"
            },
            Topology::Circular,
        )
        .unwrap();
        let imported = call(
            &c,
            "import_document",
            json!({"text":collagene_core::formats::export_genbank(&target).unwrap()}),
        )
        .await;
        let target = doc(&imported);
        let plan=call(&c,"propose_plan",json!({"name":"actual selected cut path","inputs":[{"id":"template","document":source}],"target":target,
            "steps":[Operation::DeclareSourceModifications{template:Source::Input("template".into()),assumptions:collagene_core::modified_duplex::GenerationAssumptions{dam:collagene_core::methylation::MethylationState::Unmodified,dcm:collagene_core::methylation::MethylationState::Unmodified,cpg:collagene_core::methylation::MethylationState::Unmodified}},Operation::GoldenGate{templates:vec![Source::Step(0)],options:GoldenGateOptions{max_fragments:1,..Default::default()},candidate_index:usize::from(reverse)}]})).await;
        let query = json!({"plan_revision":plan["plan_revision"],"options":collagene_core::assembly_junctions::Options{flank_bases:4,..Default::default()}});
        let inspected = call(&c, "inspect_assembly_junctions", query.clone()).await;
        assert_eq!(
            inspected["verification_status"],
            if wrong {
                "Invalid"
            } else {
                "ConsistentUnderAssumptions"
            }
        );
        assert_eq!(inspected["product_length"], 20);
        assert_eq!(inspected["junction_count"], 1);
        assert_eq!(inspected["complete"], true);
        let method = &inspected["method"]["GoldenGate"];
        assert_eq!(method["candidate_index"], usize::from(reverse));
        assert_eq!(method["path"][0]["top_start"], 11);
        assert_eq!(method["path"][0]["bottom_start"], 15);
        assert_eq!(method["path"][0]["template"], json!({"Step":0}));

        assert_eq!(method["path"][0]["reverse"], reverse);
        assert_eq!(inspected["junctions"][0]["boundary"], 0);
        assert_eq!(
            inspected["junctions"][0]["window"]["dna"],
            if reverse { "GGGGTCAT" } else { "CAGAATGA" }
        );
        let evidence = full(&c, &inspected["artifact_id"]).await;
        assert_eq!(evidence["method"], inspected["method"]);
        let declared: Molecule =
            serde_json::from_value(evidence["verification"]["products"][0].clone()).unwrap();
        assert_eq!(method["path"][0]["template_revision"], declared.revision());
        assert_eq!(
            evidence["verification"]["products"]
                .as_array()
                .unwrap()
                .len(),
            2
        );
        assert_eq!(evidence["context_revision"], inspected["context_revision"]);
        assert_eq!(
            call(&c, "inspect_assembly_junctions", query).await["artifact_id"],
            inspected["artifact_id"]
        );
    }
    assert_eq!(call(&c, "export_document", source).await, before);
    c.cancel().await.unwrap();
    server.await.unwrap();
}
