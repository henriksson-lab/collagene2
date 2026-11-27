#[path = "../../core/tests/common/domestication_oe.rs"]
mod fixture;
use base64::{Engine as _, engine::general_purpose::STANDARD};
use collagene_core::{
    Topology, verification,
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
async fn artifact(c: &Client, id: &Value) -> (Vec<u8>, String) {
    let mut text = String::new();
    let mut start = 0;
    let encoding;
    loop {
        let r = call(
            c,
            "read_artifact",
            json!({"artifact_id":id,"start":start,"max_bytes":16384}),
        )
        .await;
        text.push_str(r["text"].as_str().unwrap());
        if r["done"] == true {
            encoding = r["encoding"].as_str().unwrap().to_owned();
            break;
        }
        start = r["next_start"].as_u64().unwrap();
    }
    let bytes = if encoding == "base64" {
        STANDARD.decode(text).unwrap()
    } else {
        text.into_bytes()
    };
    assert_eq!(blake3::hash(&bytes).to_hex().as_str(), id.as_str().unwrap());
    (bytes, encoding)
}

/// Exercise the real SDK with a native archived product, not a metadata-only stock.
#[tokio::test]
async fn derived_stock_domestication_requires_exact_owner_and_exports_only_new_primers() {
    use collagene_core::{
        Feature, Location, Sequence,
        analysis::Primer,
        chemistry::Strand,
        domestication as design, domestication_preparation as physical,
        history::DocumentHistory,
        material_derivation as derivation,
        run_lineage::EventBinding,
        run_material::{self, MaterialOrigin},
    };
    let (mut original, mut request) = fixture::fixture(Topology::Linear);
    original.molecule.features.push(Feature {
        kind: "CDS".into(),
        location: Location::range(0, 12),
        qualifiers: vec![("gene".into(), Some("Retained seed coding feature".into()))],
        needs_review: false,
    });
    original.revision = original.molecule.revision();
    let w = Workspace::from_documents("Saved source owner", vec![original.molecule.clone()], 0)
        .unwrap();
    let doc = &w.documents[0];
    let seed = EventBinding {
        document_id: doc.id.clone(),
        event_id: doc.head_event.clone().unwrap(),
        revision: doc.head.clone(),
    };
    let forward = Primer::new("Upstream stock F", &fixture::ORIGINAL[..8], "").unwrap();
    let reverse = Primer::new(
        "Upstream stock R",
        Sequence::dna(&fixture::ORIGINAL[fixture::ORIGINAL.len() - 8..])
            .unwrap()
            .reverse_complement()
            .unwrap()
            .as_str(),
        "",
    )
    .unwrap();
    let source_origin = derivation::SourceOrigin::OriginalReference;
    let products = derivation::preview_pcr(
        &original.molecule,
        &source_origin,
        &forward,
        &reverse,
        Strand::Top,
    )
    .unwrap();
    assert_eq!(products.len(), 1);
    let stock = derivation::prepare_install(
        &w,
        &derivation::Request {
            name: "Actual PCR stock".into(),
            source: seed,
            origin: source_origin.clone(),
            save_primers: true,
            derivation: derivation::Derivation::ExactPcr {
                template_strand: Strand::Top,
                selection: derivation::PcrSelection {
                    selected_index: 0,
                    expected_preview_revision: products[0].molecule.revision(),
                    geometry: products[0].geometry.clone(),
                },
                forward,
                reverse,
            },
        },
    )
    .unwrap();
    assert!(stock.receipt.primer_attachment.is_some());
    let w =
        workspace::load_workspace(&workspace::save_workspace(&stock.workspace).unwrap()).unwrap();
    let event = stock.receipt.product;
    let origin = MaterialOrigin::Derived(derivation::capsule(&w, &event).unwrap());
    let actual = w.revisions.get(&event.revision).unwrap();
    assert_eq!(actual.sequence.as_str(), fixture::ORIGINAL);
    assert_eq!(actual.features.len(), 1);
    // Pins are rebuilt from the actual native product; no original feature/body alias.
    request.origin = Some(origin.clone());
    request.design.source_revision = actual.revision();
    request.design.coding_features = vec![design::feature_pin(actual, 0).unwrap()];
    request.selection.context_revision = design::propose_revision(actual, &request.design).unwrap();
    request.selection.candidate_id =
        design::candidate_id(&request.design, &request.selection.patches).unwrap();
    let checked = design::check(actual, &request.design, &request.selection).unwrap();
    assert_eq!(checked.target().sequence.as_str(), fixture::INTENDED);
    let archive = workspace::save_workspace(&w).unwrap();
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
    let imported = call(
        &c,
        "import_workspace_archive",
        json!({"data_base64":STANDARD.encode(&archive)}),
    )
    .await;
    let owner =
        json!({"id":imported["workspace"]["id"], "revision":imported["workspace"]["revision"]});
    let source = json!({"id":actual.revision(), "revision":actual.revision()});
    let before_source = call(&c, "export_document", source.clone()).await;
    let before_workspaces = call(&c, "list_workspaces", json!({})).await;
    let before_documents = call(&c, "list_documents", json!({})).await;
    reject(
        &c,
        "prepare_domestication",
        json!({"source":source,"request":request}),
    )
    .await;
    reject(
        &c,
        "prepare_domestication",
        json!({"source":source,"request":request,"origin_workspace":null}),
    )
    .await;
    assert_eq!(
        call(&c, "list_documents", json!({})).await,
        before_documents
    );
    assert_eq!(
        call(&c, "list_workspaces", json!({})).await,
        before_workspaces
    );
    let mut marker_only = request.clone();
    marker_only.origin = None;
    reject(
        &c,
        "prepare_domestication",
        json!({"source":source,"request":marker_only,"origin_workspace":owner}),
    )
    .await;
    assert_eq!(
        call(&c, "list_documents", json!({})).await,
        before_documents
    );
    let foreign_import = call(&c, "import_workspace_archive", json!({"data_base64":STANDARD.encode(workspace::save_workspace(&Workspace::from_documents("Foreign unrelated owner",vec![original.molecule.clone()],0).unwrap()).unwrap())})).await;
    let foreign = json!({"id":foreign_import["workspace"]["id"],"revision":foreign_import["workspace"]["revision"]});
    let after_foreign_import = call(&c, "list_documents", json!({})).await;
    reject(
        &c,
        "prepare_domestication",
        json!({"source":source,"request":request,"origin_workspace":foreign}),
    )
    .await;
    assert_eq!(
        call(&c, "list_documents", json!({})).await,
        after_foreign_import
    );
    // An equal complete molecule under a different event is not the retained origin.
    let mut aliases = w.clone();
    let histories = DocumentHistory::restore_all(&w).unwrap();
    let h = histories
        .iter()
        .find(|h| h.id == event.document_id)
        .unwrap();
    let alias = h
        .fork_event(&event.event_id, "different-equal-stock-event")
        .unwrap();
    let alias_event = EventBinding {
        document_id: alias.id.clone(),
        event_id: alias.head_event().into(),
        revision: alias.head_revision().into(),
    };
    alias.contribute(&mut aliases).unwrap();
    let imported_alias = call(
        &c,
        "import_workspace_archive",
        json!({"data_base64":STANDARD.encode(workspace::save_workspace(&aliases).unwrap())}),
    )
    .await;
    let alias_owner = json!({"id":imported_alias["workspace"]["id"],"revision":imported_alias["workspace"]["revision"]});
    let mut forged = request.clone();
    match forged.origin.as_mut().unwrap() {
        MaterialOrigin::Derived(cap) => cap.product = alias_event,
        _ => unreachable!(),
    }
    let before_rejections = call(&c, "list_documents", json!({})).await;
    let workspace_rejections = call(&c, "list_workspaces", json!({})).await;
    reject(
        &c,
        "prepare_domestication",
        json!({"source":source,"request":forged,"origin_workspace":alias_owner}),
    )
    .await;
    let seed_source = json!({"id":original.revision,"revision":original.revision});
    reject(
        &c,
        "prepare_domestication",
        json!({"source":seed_source,"request":request,"origin_workspace":owner}),
    )
    .await;
    assert_eq!(
        call(&c, "list_documents", json!({})).await,
        before_rejections
    );
    assert_eq!(
        call(&c, "list_workspaces", json!({})).await,
        workspace_rejections
    );
    assert_eq!(
        call(&c, "export_document", source.clone()).await,
        before_source
    );
    assert_eq!(before_documents["total"], after_foreign_import["total"]);
    let prepared = call(
        &c,
        "prepare_domestication",
        json!({"source":source,"request":request,"origin_workspace":owner}),
    )
    .await;
    assert_eq!(prepared["verified"], false);
    assert_eq!(prepared["proposal_ready"], true);
    assert_eq!(
        prepared["verification"]["status"],
        "ConsistentUnderAssumptions"
    );
    assert_eq!(prepared["operation_count"], 5);
    let (bytes, _) = artifact(&c, &prepared["artifact_id"]).await;
    let full: Value = serde_json::from_slice(&bytes).unwrap();
    let native: physical::Prepared = serde_json::from_value(full["prepared"].clone()).unwrap();
    assert_eq!(native.target.sequence.as_str(), fixture::INTENDED);
    assert_eq!(native.target.sequence, checked.target().sequence);
    assert!(
        matches!(&native.plan.steps[0],verification::Operation::ReconstructMaterial{template:verification::Source::Input(id),origin:cap} if id==physical::INPUT_ID && cap.product==event)
    );
    native
        .validate_binding(
            &physical::revision(
                &verification::PlanInput {
                    id: physical::INPUT_ID.into(),
                    revision: actual.revision(),
                    molecule: actual.clone(),
                },
                &request,
                Some(&w),
            )
            .unwrap(),
        )
        .unwrap();
    let oligos = call(
        &c,
        "export_plan_oligos",
        json!({"plan_revision":prepared["plan_revision"],"format":"Csv"}),
    )
    .await;
    assert_eq!(oligos["complete"], true);
    assert_eq!(oligos["use_count"], 6);
    let (csv, _) = artifact(&c, &oligos["table_artifact_id"]).await;
    let csv = String::from_utf8(csv).unwrap();
    assert!(!csv.contains("Upstream stock"));
    let wrong = call(&c,"propose_plan",json!({"name":"Wrong independent saved-stock target","target":source,"inputs":[{"id":physical::INPUT_ID,"document":source}],"steps":native.plan.steps})).await;
    assert_eq!(
        call(
            &c,
            "verify_plan",
            json!({"plan_revision":wrong["plan_revision"]})
        )
        .await["status"],
        "Invalid"
    );
    assert_eq!(call(&c, "export_document", source).await, before_source);
    assert_eq!(workspace::save_workspace(&w).unwrap(), archive);
    c.cancel().await.unwrap();
    server.await.unwrap();
}
