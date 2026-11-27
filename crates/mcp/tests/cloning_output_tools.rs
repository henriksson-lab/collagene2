//! Finalization retains exact native event history and source-bound order evidence.
use base64::{Engine, engine::general_purpose::STANDARD};
use collagene_core::{
    Molecule, Topology,
    analysis::Primer,
    chemistry::AnnealOptions,
    oligo_formats,
    workspace::{self, Workspace},
};
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
    assert_ne!(r.is_error, Some(true), "{name}:{:?}", r.content);
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
    let (mut at, mut text) = (0, String::new());
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
fn doc(m: &Molecule) -> Value {
    json!({"id":m.revision(),"revision":m.revision()})
}
fn binding(w: &Workspace, i: usize) -> Value {
    json!({"document_id":w.documents[i].id,"event_id":w.documents[i].head_event,"revision":w.documents[i].head})
}
fn workspace_ref(result: &Value) -> Value {
    json!({"id":result["workspace"]["id"],"revision":result["workspace"]["revision"]})
}
async fn native(c: &Client, id: &Value) -> Vec<u8> {
    let (mut at, mut text) = (0, String::new());
    loop {
        let page = call(
            c,
            "read_artifact",
            json!({"artifact_id":id,"start":at,"max_bytes":16384}),
        )
        .await;
        assert_eq!(page["encoding"], "base64");
        text.push_str(page["text"].as_str().unwrap());
        if page["done"] == true {
            break;
        }
        at = page["next_start"].as_u64().unwrap();
    }
    STANDARD.decode(text).unwrap()
}
async fn upload(c: &Client, w: &Workspace) -> Value {
    let bytes = workspace::save_workspace(w).unwrap();
    let begin = call(
        c,
        "begin_workspace_upload",
        json!({"expected_bytes":bytes.len(),"blake3":blake3::hash(&bytes).to_hex().to_string()}),
    )
    .await;
    for (offset, chunk) in bytes.chunks(4096).enumerate() {
        call(c,"append_workspace_upload",json!({"upload_id":begin["upload_id"],"offset":offset*4096,"data_base64":STANDARD.encode(chunk)})).await;
    }
    call(
        c,
        "finalize_workspace_upload",
        json!({"upload_id":begin["upload_id"]}),
    )
    .await
}
#[tokio::test]
async fn exact_event_finalization_saved_primers_native_history_and_replay_are_atomic_and_idempotent()
 {
    let (c, server) = setup().await;
    assert_eq!(
        c.list_tools(Default::default()).await.unwrap().tools.len(),
        collagene_mcp::CURRENT_TOOL_COUNT
    );
    let left = Molecule::dna("left", "ATGACCCCGGGGTTTCCAGA", Topology::Linear).unwrap();
    let right = Molecule::dna("right", "ATCCGTTAGCGACCTTGGAA", Topology::Linear).unwrap();
    let target = Molecule::dna(
        "independent target",
        "ATGACCCCGGGGTTTCCAGAATCCGTTAGCGACCTTGGAA",
        Topology::Linear,
    )
    .unwrap();
    let w = Workspace::from_documents(
        "Finalize fixture",
        vec![left.clone(), right.clone(), target.clone()],
        2,
    )
    .unwrap();
    let base_bytes = workspace::save_workspace(&w).unwrap();
    let imported = upload(&c, &w).await;
    let wr = workspace_ref(&imported);
    let recipe = json!({"name":"Actual prepared overlap","fragments":[{"generation":null,"role":"left","input_id":"l","preparation":{"ExactPcr":{"template_strand":"Top","forward":Primer::new("LF","ATGACCCC","").unwrap(),"reverse":Primer::new("LR","TCTGGAAA","AGGTCGCTAACGGAT").unwrap()}},"reverse":false},{"generation":null,"role":"right","input_id":"r","preparation":{"ExactPcr":{"template_strand":"Top","forward":Primer::new("RF","ATCCGTTA","").unwrap(),"reverse":Primer::new("RR","TTCCAAGG","").unwrap()}},"reverse":false}],"overlaps":[15],"circular":false,"options":{"method":"HiFi","scan_repeated_overlaps":false,"max_comparisons":100000,"max_occurrences":100}});
    let inputs = json!([{"id":"l","document":doc(&left)},{"id":"r","document":doc(&right)}]);
    let plan = call(
        &c,
        "propose_overlap_assembly",
        json!({"inputs":inputs,"target":doc(&target),"recipe":recipe}),
    )
    .await;
    assert_eq!(plan["verification"]["status"], "ConsistentUnderAssumptions");
    let args = json!({"workspace":wr,"plan_revision":plan["plan_revision"],"target_binding":binding(&w,2),"input_bindings":[{"input_id":"l","event":binding(&w,0)},{"input_id":"r","event":binding(&w,1)}],"options":{"name":"Reviewed construct","show_prepared":false,"save_primers":true,"allow_incomplete":false}});
    let before = call(&c, "list_workspaces", json!({})).await["workspaces"]
        .as_array()
        .unwrap()
        .len();
    let mut wrong = args.clone();
    wrong["target_binding"] = binding(&w, 0);
    reject(&c, "finalize_cloning_outputs", wrong).await;
    let mut stale = args.clone();
    stale["input_bindings"][0]["event"]["event_id"] = json!(format!("event-{}", "0".repeat(64)));
    reject(&c, "finalize_cloning_outputs", stale).await;
    let mut switched = args.clone();
    switched["input_bindings"][0]["event"] = binding(&w, 1);
    reject(&c, "finalize_cloning_outputs", switched).await;
    assert_eq!(
        call(&c, "list_workspaces", json!({})).await["workspaces"]
            .as_array()
            .unwrap()
            .len(),
        before
    );
    let blocked_fragments = recipe["fragments"]
        .as_array()
        .unwrap()
        .iter()
        .map(|source| json!({"source":source,"treatment":null}))
        .collect::<Vec<_>>();
    let blocked_plan=call(&c,"propose_ligation_assembly",json!({"inputs":inputs,"target":doc(&target),"recipe":{"name":"Unphosphorylated PCR ends cannot seal","fragments":blocked_fragments,"circular":false,"allow_nicked":false}})).await;
    assert_ne!(
        blocked_plan["verification"]["status"],
        "ConsistentUnderAssumptions"
    );
    let mut blocked_finalize = args.clone();
    blocked_finalize["plan_revision"] = blocked_plan["plan_revision"].clone();
    reject(&c, "finalize_cloning_outputs", blocked_finalize).await;
    assert_eq!(
        call(&c, "list_workspaces", json!({})).await["workspaces"]
            .as_array()
            .unwrap()
            .len(),
        before
    );
    let done = call(&c, "finalize_cloning_outputs", args.clone()).await;
    let next = workspace_ref(&done);
    assert_ne!(next, wr);
    assert_eq!(
        done["receipt"]["run"]["products"].as_array().unwrap().len(),
        3
    );
    assert_eq!(
        done["receipt"]["visible_products"]
            .as_array()
            .unwrap()
            .len(),
        1
    );
    assert_eq!(
        done["receipt"]["hidden_products"].as_array().unwrap().len(),
        2
    );
    let bytes = native(&c, &done["native_artifact_id"]).await;
    let decoded = workspace::load_workspace(&bytes).unwrap();
    assert_eq!(decoded.documents[0].head, w.documents[0].head);
    assert_eq!(decoded.documents[1].head, w.documents[1].head);
    assert_eq!(decoded.documents[2].head, w.documents[2].head);
    let construct = &done["receipt"]["construct"];
    assert_eq!(
        decoded.revisions[construct["revision"].as_str().unwrap()].name,
        "Reviewed construct"
    );
    let full = artifact(&c, &done["artifact_id"]).await;
    let commit: collagene_core::cloning_outputs::Commit = serde_json::from_value(full).unwrap();
    commit
        .validate_binding(done["receipt"]["context_revision"].as_str().unwrap())
        .unwrap();
    for product in &commit.receipt.run.products {
        assert!(
            decoded
                .history
                .iter()
                .any(|event| event.id == product.event_id && event.revision == product.revision)
        );
        assert!(decoded.revisions.contains_key(&product.revision));
    }
    let primer_hash = commit.receipt.primer_attachment.as_ref().unwrap();
    assert_eq!(decoded.attachments[primer_hash].kind, "oligo-table-v1");
    let oligos: oligo_formats::OligoTable =
        serde_json::from_slice(&decoded.attachments[primer_hash].bytes).unwrap();
    oligos.validate().unwrap();
    assert_eq!(oligos.records.len(), 4);
    let replay = call(
        &c,
        "replay_workspace_run",
        json!({"workspace":next,"run_id":commit.receipt.run.run_id}),
    )
    .await;
    assert_eq!(replay["products_match"], true);
    assert_eq!(replay["status_matches"], true);
    let again = call(&c, "finalize_cloning_outputs", args).await;
    assert_eq!(again["workspace"]["id"], done["workspace"]["id"]);
    assert_eq!(again["native_artifact_id"], done["native_artifact_id"]);
    assert_eq!(again["artifact_id"], done["artifact_id"]);
    assert_eq!(
        call(&c, "list_workspaces", json!({})).await["workspaces"]
            .as_array()
            .unwrap()
            .len(),
        before + 1
    );
    let reimported = call(
        &c,
        "import_workspace_archive",
        json!({"data_base64":STANDARD.encode(&bytes)}),
    )
    .await;
    assert_eq!(workspace_ref(&reimported), workspace_ref(&done));
    let original = call(&c, "export_workspace_archive", wr).await;
    assert_eq!(native(&c, &original["artifact_id"]).await, base_bytes);
    drop(c);
    server.abort();
}
#[tokio::test]
async fn incomplete_opaque_order_data_never_publishes_a_finalized_workspace() {
    let (c, server) = setup().await;
    let mapping = oligo_formats::OligoImportOptions {
        columns: oligo_formats::ColumnMapping {
            name: 0,
            ordered_sequence: Some(1),
            modification: Some(2),
            ..Default::default()
        },
        ..Default::default()
    };
    let table = oligo_formats::import_oligos(
        "Name,Sequence,Modification\nopaque,ACGT,unknown supplied modification\n",
        &mapping,
    )
    .unwrap();
    let record = &table.records[0];
    let mut top = Molecule::dna("opaque starting oligo", "ACGT", Topology::Linear).unwrap();
    top.metadata
        .insert("genbank.molecule_type".into(), "ss-DNA".into());
    top.metadata.insert(
        "oligo:original_record_base64".into(),
        record.original_record_base64.clone(),
    );
    top.metadata.insert(
        "oligo:original_record_blake3".into(),
        record.original_record_blake3.clone(),
    );
    top.metadata.insert(
        "oligo:source_mapping".into(),
        serde_json::to_string(&mapping).unwrap(),
    );
    top.metadata.insert(
        "oligo:source_fields".into(),
        serde_json::to_string(&record.original_fields).unwrap(),
    );
    let mut bottom = Molecule::dna("other starting oligo", "ACGT", Topology::Linear).unwrap();
    bottom
        .metadata
        .insert("genbank.molecule_type".into(), "ss-DNA".into());
    let target = Molecule::dna("independent anneal target", "ACGT", Topology::Linear).unwrap();
    let w = Workspace::from_documents(
        "Opaque order fixture",
        vec![top.clone(), bottom.clone(), target.clone()],
        2,
    )
    .unwrap();
    let imported = upload(&c, &w).await;
    let mut options = AnnealOptions::default();
    options.min_paired_bases = 4;
    let plan=call(&c,"propose_plan",json!({"name":"Declared anneal source orders","target":doc(&target),"inputs":[{"id":"top","document":doc(&top)},{"id":"bottom","document":doc(&bottom)}],"steps":[{"AnnealOligos":{"top":{"Input":"top"},"bottom":{"Input":"bottom"},"options":options,"candidate_index":0}}]})).await;
    let orders = call(
        &c,
        "export_plan_oligos",
        json!({"plan_revision":plan["plan_revision"],"format":"Csv"}),
    )
    .await;
    assert_eq!(orders["complete"], false);
    assert!(!orders["unavailable"].as_array().unwrap().is_empty());
    let before = call(&c, "list_workspaces", json!({})).await["workspaces"]
        .as_array()
        .unwrap()
        .len();
    reject(&c,"finalize_cloning_outputs",json!({"workspace":workspace_ref(&imported),"plan_revision":plan["plan_revision"],"target_binding":binding(&w,2),"input_bindings":[{"input_id":"top","event":binding(&w,0)},{"input_id":"bottom","event":binding(&w,1)}],"options":{"name":"Not publishable","show_prepared":false,"save_primers":true,"allow_incomplete":false}})).await;
    assert_eq!(
        call(&c, "list_workspaces", json!({})).await["workspaces"]
            .as_array()
            .unwrap()
            .len(),
        before
    );
    drop(c);
    server.abort();
}
