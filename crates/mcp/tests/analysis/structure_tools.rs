//! Actual SDK over memory pipes; no browser, HTTP, user files or scientific process.
use collagene_mcp::server::CollageneServer;
use rmcp::{RoleClient, ServiceExt, model::CallToolRequestParams, service::RunningService};
use serde_json::{Value, json};
type Client = RunningService<RoleClient, ()>;
const SOURCE: &str = "ATOM      1  CA AALA A  -2B     11.104  13.207   2.100  0.60 20.00           C  \r\nHETATM    2 CA    CA Z  17      12.000  13.000   2.000  1.00 10.00          CA2+\r\nCONECT    1    2\r\nEND\r\n";
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
    match c
        .call_tool(
            CallToolRequestParams::new(name.to_owned())
                .with_arguments(args.as_object().unwrap().clone()),
        )
        .await
    {
        Ok(r) => assert_eq!(r.is_error, Some(true)),
        Err(rmcp::service::ServiceError::McpError(e)) => assert_eq!(e.code.0, -32602),
        Err(e) => panic!("unexpected transport failure: {e}"),
    }
}
#[tokio::test]
async fn typed_import_reinspection_original_artifacts_and_source_roles_remain_documentary() {
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
    let tools = c.list_tools(Default::default()).await.unwrap();
    assert_eq!(tools.tools.len(), collagene_mcp::CURRENT_TOOL_COUNT);
    for name in ["import_structure", "inspect_structure"] {
        let tool = tools.tools.iter().find(|t| t.name == name).unwrap();
        let schema = serde_json::to_string(&tool.input_schema).unwrap();
        assert!(schema.contains("model") && schema.contains("hydrogens"));
    }
    let before = call(&c, "list_documents", json!({})).await["total"].clone();
    let imported = call(
        &c,
        "import_structure",
        json!({"text":SOURCE,"format":"Pdb"}),
    )
    .await;
    assert_eq!(imported["atom_count"], 2);
    assert_eq!(imported["models"], json!([1]));
    assert_eq!(imported["selected_atom_count"], 2);
    assert_eq!(imported["display_edge_count"], 1);
    assert_eq!(
        imported["source_hash"],
        blake3::hash(SOURCE.as_bytes()).to_hex().as_str()
    );
    let source = imported["source"].clone();
    assert_eq!(source["id"], source["revision"]);
    let mut original = String::new();
    let mut start = 0;
    loop {
        let part = call(
            &c,
            "read_artifact",
            json!({"artifact_id":source["id"],"start":start,"max_bytes":31}),
        )
        .await;
        assert_eq!(part["encoding"], "utf-8");
        original.push_str(part["text"].as_str().unwrap());
        if part["done"] == true {
            break;
        }
        start = part["next_start"].as_u64().unwrap();
    }
    assert_eq!(original.as_bytes(), SOURCE.as_bytes());
    let full = call(
        &c,
        "read_artifact",
        json!({"artifact_id":imported["artifact_id"],"start":0,"max_bytes":16384}),
    )
    .await;
    assert_eq!(full["done"], true);
    let report: collagene_core::protein::structure::Prepared =
        serde_json::from_str(full["text"].as_str().unwrap()).unwrap();
    report
        .validate_binding(imported["context_revision"].as_str().unwrap())
        .unwrap();
    assert_eq!(report.atoms[0].position, [11.104, 13.207, 2.1]);
    assert_eq!(report.scene.edges[0].atoms, [0, 1]);
    let inspected=call(&c,"inspect_structure",json!({"source":source,"selection":{"model":1,"chains":{"Author":"Z"},"alternates":"All","hydrogens":true}})).await;
    assert_eq!(inspected["selected_atom_count"], 1);
    assert_eq!(inspected["atom_preview"][0]["source_atom_index"], 1);
    assert_eq!(inspected["unresolved_connections"], 1);
    assert_eq!(inspected["source"], source);
    assert_ne!(inspected["context_revision"], imported["context_revision"]);
    let repeated = call(&c, "inspect_structure", json!({"source":source})).await;
    assert_eq!(repeated["artifact_id"], imported["artifact_id"]);
    let mut stale = source.clone();
    stale["revision"] = json!("f".repeat(64));
    reject(&c, "inspect_structure", json!({"source":stale})).await;
    let foreign = json!({"id":imported["artifact_id"],"revision":imported["artifact_id"]});
    reject(&c, "inspect_structure", json!({"source":foreign})).await;
    reject(&c,"inspect_structure",json!({"source":source,"selection":{"model":999,"chains":"All","alternates":"All","hydrogens":true}})).await;
    reject(
        &c,
        "import_structure",
        json!({"text":"ATOM      1","format":"Pdb"}),
    )
    .await;
    reject(
        &c,
        "import_structure",
        json!({"text":SOURCE,"format":"Mmcif"}),
    )
    .await;
    assert_eq!(call(&c, "list_documents", json!({})).await["total"], before);
    c.cancel().await.unwrap();
    server.await.unwrap();
}

#[tokio::test]
async fn mmcif_source_role_and_label_author_filters_keep_nonpolymer_identity() {
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
    let text = "data_authored\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.auth_asym_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nHETATM 7 O O HOH W . Water 4.0 5.0 6.0\n";
    let imported = call(
        &c,
        "import_structure",
        json!({"text":text,"format":"Mmcif"}),
    )
    .await;
    assert_eq!(imported["atom_count"], 1);
    assert_eq!(imported["atom_preview"][0]["atom"]["label"]["chain"], "W");
    assert_eq!(
        imported["atom_preview"][0]["atom"]["author"]["chain"],
        "Water"
    );
    assert_eq!(
        imported["atom_preview"][0]["atom"]["label"]["sequence"],
        Value::Null
    );
    let original = call(
        &c,
        "read_artifact",
        json!({"artifact_id":imported["source"]["id"],"start":0,"max_bytes":16384}),
    )
    .await;
    assert_eq!(original["mime_type"], "chemical/x-mmcif");
    assert_eq!(original["filename"], "structure.cif");
    assert_eq!(original["text"], text);
    let filtered=call(&c,"inspect_structure",json!({"source":imported["source"],"selection":{"model":1,"chains":{"Label":"W"},"alternates":"All","hydrogens":false}})).await;
    assert_eq!(filtered["selected_atom_count"], 1);
    assert_eq!(
        filtered["atom_preview"][0]["atom"]["position"],
        json!([4.0, 5.0, 6.0])
    );
    reject(&c,"inspect_structure",json!({"source":imported["source"],"selection":{"model":1,"chains":{"Label":"Water"},"alternates":"All","hydrogens":true}})).await;
    c.cancel().await.unwrap();
    server.await.unwrap();
}
