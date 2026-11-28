//! In-memory SDK regression coverage for documentary SBOL exchange and authority boundaries.
use base64::{Engine as _, engine::general_purpose::STANDARD};
use collagene_mcp::CollageneServer;
use collagene_model::{Molecule, Topology};
use rmcp::{RoleClient, ServiceExt, model::CallToolRequestParams, service::RunningService};
use serde_json::{Value, json};

type Client = RunningService<RoleClient, ()>;
const COMPONENT: &str = "https://example.org/sdk/component";
const SEQUENCE: &str = "https://example.org/sdk/sequence";
const HIERARCHY: &str = include_str!("../../../design-exchange/tests/fixtures/hierarchy.jsonld");

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
    assert!(
        result.is_err() || result.unwrap().is_error == Some(true),
        "{name} unexpectedly accepted"
    );
}

async fn artifact(c: &Client, id: &Value) -> Vec<u8> {
    let mut text = String::new();
    let mut start = 0;
    let encoding = loop {
        let page = call(
            c,
            "read_artifact",
            json!({"artifact_id":id,"start":start,"max_bytes":4096}),
        )
        .await;
        text.push_str(page["text"].as_str().unwrap());
        if page["done"] == true {
            break page["encoding"].as_str().unwrap().to_owned();
        }
        let next = page["next_start"].as_u64().unwrap();
        assert!(next > start);
        start = next;
    };
    match encoding.as_str() {
        "base64" => STANDARD.decode(text).unwrap(),
        "utf-8" => text.into_bytes(),
        other => panic!("unexpected artifact encoding: {other}"),
    }
}

fn graph(single_stranded: bool) -> String {
    let mut types = vec![
        json!({"@id":"https://identifiers.org/SBO:0000251"}),
        json!({"@id":"https://identifiers.org/SO:0000987"}),
    ];
    if single_stranded {
        types.push(json!({"@id":"https://identifiers.org/SO:0000984"}));
    }
    // Independent literal graph: this fixture does not use the exporter being tested.
    format!(
        "\n{}\n",
        json!({
            "@context":{"sbol":"http://sbols.org/v3#"},
            "@graph":[
                {"@id":COMPONENT,"@type":"sbol:Component",
                 "sbol:hasNamespace":{"@id":"https://example.org/sdk"},
                 "sbol:type":types,"sbol:hasSequence":{"@id":SEQUENCE}},
                {"@id":SEQUENCE,"@type":"sbol:Sequence",
                 "sbol:hasNamespace":{"@id":"https://example.org/sdk"},
                 "sbol:elements":"ACGTTGCA","sbol:encoding":{"@id":"https://identifiers.org/edam:format_1207"}}
            ]
        })
    )
}

fn reference(summary: &Value) -> Value {
    json!({"id":summary["id"],"revision":summary["revision"]})
}

async fn documents(c: &Client) -> Value {
    call(c, "list_documents", json!({"start":0,"limit":128})).await
}

#[tokio::test]
async fn retained_graph_exports_are_exact_and_source_roles_and_paging_are_enforced() {
    let (c, server) = setup().await;
    let original = graph(false);
    let imported = call(&c, "import_sbol", json!({"text":original})).await;
    let source = imported["source"].clone();
    assert_eq!(imported["component_count"], 1);
    assert_eq!(imported["object_count"], 2);
    assert_eq!(documents(&c).await["total"], 0);
    let exact = call(
        &c,
        "export_sbol_graph",
        json!({"source":source,"representation":"Original"}),
    )
    .await;
    assert_eq!(
        artifact(&c, &exact["artifact_id"]).await,
        original.as_bytes()
    );
    assert_eq!(
        exact["blake3"],
        blake3::hash(original.as_bytes()).to_hex().to_string()
    );
    let expanded = call(
        &c,
        "export_sbol_graph",
        json!({"source":source,"representation":"Expanded"}),
    )
    .await;
    let expanded_bytes = artifact(&c, &expanded["artifact_id"]).await;
    let expanded_graph: Value = serde_json::from_slice(&expanded_bytes).unwrap();
    let expanded_text = std::str::from_utf8(&expanded_bytes).unwrap();
    assert!(expanded_graph.is_array() || expanded_graph.is_object());
    assert!(expanded_text.contains("http://sbols.org/v3#Component"));
    assert!(expanded_text.contains(COMPONENT));
    assert!(expanded_text.contains("ACGTTGCA"));
    let inspected = call(
        &c,
        "inspect_sbol",
        json!({"source":source,"component_start":0,"max_components":1}),
    )
    .await;
    assert_eq!(inspected["components"][0]["component"], COMPONENT);
    for (start, maximum) in [(0, 0), (0, 51), (2, 1)] {
        reject(
            &c,
            "inspect_sbol",
            json!({"source":source,"component_start":start,"max_components":maximum}),
        )
        .await;
    }
    let mut stale = source.clone();
    stale["Imported"]["revision"] = json!("stale");
    reject(
        &c,
        "export_sbol_graph",
        json!({"source":stale,"representation":"Original"}),
    )
    .await;
    let foreign =
        json!({"Imported":{"id":imported["artifact_id"],"revision":imported["artifact_id"]}});
    reject(
        &c,
        "inspect_sbol",
        json!({"source":foreign,"component_start":0,"max_components":1}),
    )
    .await;
    assert_eq!(documents(&c).await["total"], 0);
    c.cancel().await.unwrap();
    server.await.unwrap();
}

#[tokio::test]
async fn hierarchy_projection_requires_review_and_roundtrips_only_reference_and_exact_source() {
    let (c, server) = setup().await;
    let imported = call(&c, "import_sbol", json!({"text":HIERARCHY})).await;
    let source = imported["source"].clone();
    let selection = json!({"component":"https://example.org/lab/plasmid", "sequence":"https://example.org/lab/plasmid_sequence", "topology":null});
    let review = call(
        &c,
        "project_sbol",
        json!({"source":source,"selection":selection}),
    )
    .await;
    assert_eq!(review["state"], "review_required");
    assert!(review["document"].is_null());
    assert!(review["workspace"].is_null());
    assert!(review["issue_count"].as_u64().unwrap() > 0);
    assert_eq!(documents(&c).await["total"], 0);
    let accepted = call(
        &c,
        "project_sbol",
        json!({"source":source,"selection":selection,"accept_omissions":true}),
    )
    .await;
    assert_eq!(accepted["receipt"]["state"], "created");
    assert_eq!(accepted["workspace"]["candidate_count"], 0);
    assert_eq!(documents(&c).await["total"], 1);
    let workspace = reference(&accepted["workspace"]);
    let archive = call(&c, "export_workspace_archive", workspace).await;
    assert_eq!(archive["candidates_verified"], false);
    let bytes = artifact(&c, &archive["artifact_id"]).await;
    let frame = collagene_core::workspace::load_workspace(&bytes).unwrap();
    assert_eq!(frame.documents.len(), 1);
    assert!(frame.candidates.is_empty());
    assert!(frame.active_candidate.is_none());
    assert_eq!(frame.history.len(), 1);
    assert_eq!(frame.history[0].operation, "Imported document");
    assert!(frame.history[0].parents.is_empty());
    assert_eq!(frame.attachments.len(), 1);
    let hash = imported["source_hash"].as_str().unwrap();
    let attached = &frame.attachments[hash];
    assert_eq!(
        attached.kind,
        collagene_library_sharing::SBOL_SOURCE_ATTACHMENT_KIND
    );
    assert_eq!(attached.bytes, HIERARCHY.as_bytes());
    let molecule = frame.head(&frame.active_document).unwrap();
    assert_eq!(molecule.sequence.as_str(), "AAAACCCCGGGGTTTT");
    assert_eq!(molecule.metadata["sbol:source_hash"], hash);
    assert!(molecule.metadata.keys().all(|key| key.starts_with("sbol:")));
    // A fresh SDK session demonstrates that graph access comes from the native
    // attachment, not the original session's imported-source membership set.
    let (fresh, fresh_server) = setup().await;
    let restored = call(
        &fresh,
        "import_workspace_archive",
        json!({"data_base64":STANDARD.encode(bytes)}),
    )
    .await;
    let restored_source =
        json!({"Workspace":{"workspace":reference(&restored["workspace"]),"source_hash":hash}});
    let exact = call(
        &fresh,
        "export_sbol_graph",
        json!({"source":restored_source,"representation":"Original"}),
    )
    .await;
    assert_eq!(
        artifact(&fresh, &exact["artifact_id"]).await,
        HIERARCHY.as_bytes()
    );
    assert_eq!(restored["workspace"]["candidate_count"], 0);
    fresh.cancel().await.unwrap();
    fresh_server.await.unwrap();
    c.cancel().await.unwrap();
    server.await.unwrap();
}

#[tokio::test]
async fn single_stranded_dna_cannot_be_projected_even_when_omissions_are_accepted() {
    let (c, server) = setup().await;
    let original = graph(true);
    let imported = call(&c, "import_sbol", json!({"text":original})).await;
    let source = imported["source"].clone();
    for accepted in [false, true] {
        reject(
            &c,
            "project_sbol",
            json!({"source":source,
            "selection":{"component":COMPONENT,"sequence":SEQUENCE,"topology":"Linear"},
            "accept_omissions":accepted}),
        )
        .await;
    }
    assert_eq!(documents(&c).await["total"], 0);
    let exact = call(
        &c,
        "export_sbol_graph",
        json!({"source":source,"representation":"Original"}),
    )
    .await;
    assert_eq!(
        artifact(&c, &exact["artifact_id"]).await,
        original.as_bytes()
    );
    c.cancel().await.unwrap();
    server.await.unwrap();
}

#[tokio::test]
async fn native_export_cannot_discard_physical_markers_or_mutate_its_source() {
    let (c, server) = setup().await;
    let plain = Molecule::dna("plain", "GATTACAC", Topology::Linear).unwrap();
    let plain_bytes = collagene_core::protocol::save_document(&plain).unwrap();
    let plain_import = call(
        &c,
        "import_sequence_file",
        json!({
            "data_base64":STANDARD.encode(&plain_bytes),"format":"NativeDocument"
        }),
    )
    .await;
    let plain_document = reference(&plain_import["documents"][0]);
    let ids = json!({"namespace":"https://example.org/export/",
        "component":"https://example.org/export/plain",
        "sequence":"https://example.org/export/plain_sequence"});
    let exported = call(
        &c,
        "export_sbol",
        json!({"document":plain_document,"ids":ids}),
    )
    .await;
    assert_eq!(exported["state"], "prepared");
    let graph_bytes = artifact(&c, &exported["artifact_id"]).await;
    let graph = collagene_design_exchange::import(&graph_bytes).unwrap();
    let projected = graph
        .project(
            &collagene_design_exchange::Selection {
                component: "https://example.org/export/plain".into(),
                sequence: "https://example.org/export/plain_sequence".into(),
                topology: None,
            },
            collagene_design_exchange::ProjectionPolicy::Strict,
        )
        .unwrap();
    assert_eq!(projected.molecule.sequence.as_str(), "GATTACAC");
    assert_eq!(projected.molecule.topology, Topology::Linear);
    let reimported = call(
        &c,
        "import_sbol",
        json!({
            "text":std::str::from_utf8(&graph_bytes).unwrap()
        }),
    )
    .await;
    assert_eq!(reimported["components"][0]["component"], ids["component"]);
    assert_eq!(
        reimported["components"][0]["sequences"][0]["iri"],
        ids["sequence"]
    );
    reject(&c, "export_sbol", json!({"document":plain_document,
        "ids":{"namespace":"https://example.org/export/","component":"https://example.org/export/same","sequence":"https://example.org/export/same"},
        "accept_omissions":true})).await;
    let plain_retained = call(&c, "export_document", plain_document).await;
    assert_eq!(
        artifact(&c, &plain_retained["artifact_id"]).await,
        plain_bytes
    );
    let mut molecule = Molecule::dna("guarded", "ACGTTGCA", Topology::Linear).unwrap();
    molecule.metadata.insert(
        "topo:physical_graph_required".into(),
        "retained physical evidence".into(),
    );
    let original = collagene_core::protocol::save_document(&molecule).unwrap();
    let imported = call(
        &c,
        "import_sequence_file",
        json!({"data_base64":STANDARD.encode(&original),"format":"NativeDocument"}),
    )
    .await;
    let document = reference(&imported["documents"][0]);
    for accepted in [false, true] {
        reject(&c, "export_sbol", json!({"document":document,
            "ids":{"namespace":"https://example.org/export","component":"https://example.org/export/component","sequence":"https://example.org/export/sequence"},
            "accept_omissions":accepted})).await;
    }
    let retained = call(&c, "export_document", document.clone()).await;
    assert_eq!(artifact(&c, &retained["artifact_id"]).await, original);
    assert_eq!(documents(&c).await["total"], 2);
    c.cancel().await.unwrap();
    server.await.unwrap();
}
