//! In-process SDK coverage: OE material is reconstructed without an assembly target.
use base64::{Engine as _, engine::general_purpose::STANDARD};
use collagene_core::{
    Molecule, Topology,
    analysis::Primer,
    catalog_digest,
    chemistry::Strand,
    cloning_outputs, material_derivation as native, overlap_extension,
    overlap_extension_pcr::{PrimerPair, Request as Preparation},
    run_lineage::{self, EventBinding, InputEventBinding},
    run_material::{self, MaterialSelection},
    verification::{DesignPlan, Operation, PlanInput, Source},
    workspace::{self, Workspace},
};
use collagene_mcp::{CollageneServer, MaterialDerivationOrigin as Origin};
use rmcp::{RoleClient, ServiceExt, model::CallToolRequestParams, service::RunningService};
use serde_json::{Value, json};
type Client = RunningService<RoleClient, ()>;
const ORIGINAL: &str = "AGTCCGATGCTACGTAGGACTTCAGACCTGTAGATCGTCA";
const AMPLICON: &str = "ATGCGAATTCAGTCCGATGCTACGTATTTTCAGACCTGTAGATCGTCAGAATTCGCAT";
const RIGHT_CUT: &str = "AATTCGCAT";
fn preparation() -> Preparation {
    let pair = |name: &str, f: &str, ft: &str, r: &str, rt: &str| PrimerPair {
        forward: Primer::new(format!("New {name} F"), f, ft).unwrap(),
        reverse: Primer::new(format!("New {name} R"), r, rt).unwrap(),
        template_strand: Strand::Top,
    };
    Preparation {
        left: pair("left", "AGTCCGAT", "", "TACGTAGC", "GGTCTGAAAA"),
        right: pair("right", "TTCAGACC", "GCTACGTATT", "TGACGATC", ""),
        outer: pair("outer", "AGTCCGAT", "ATGCGAATTC", "TGACGATC", "ATGCGAATTC"),
        extension: overlap_extension::Options {
            overlap: 18,
            min_paired_bases: 8,
            ..Default::default()
        },
    }
}
fn event(w: &Workspace, index: usize) -> EventBinding {
    let d = &w.documents[index];
    EventBinding {
        document_id: d.id.clone(),
        event_id: d.head_event.clone().unwrap(),
        revision: d.head.clone(),
    }
}
fn original() -> Workspace {
    Workspace::from_documents(
        "Source only; no intended assembly target",
        vec![Molecule::dna("Original reference", ORIGINAL, Topology::Linear).unwrap()],
        0,
    )
    .unwrap()
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
async fn call(c: &Client, name: &str, args: Value) -> Value {
    let response = c
        .call_tool(
            CallToolRequestParams::new(name.to_owned())
                .with_arguments(args.as_object().unwrap().clone()),
        )
        .await
        .unwrap();
    assert_ne!(
        response.is_error,
        Some(true),
        "{name}: {:?}",
        response.content
    );
    response.structured_content.unwrap()
}
async fn reject(c: &Client, args: Value) {
    match c
        .call_tool(
            CallToolRequestParams::new("material_derivation")
                .with_arguments(args.as_object().unwrap().clone()),
        )
        .await
    {
        Err(_) => {}
        Ok(response) => assert_eq!(
            response.is_error,
            Some(true),
            "Unexpected material: {response:?}"
        ),
    }
}
async fn artifact(c: &Client, id: &Value) -> (String, String) {
    let (mut at, mut text) = (0, String::new());
    loop {
        let page = call(
            c,
            "read_artifact",
            json!({"artifact_id":id,"start":at,"max_bytes":16384}),
        )
        .await;
        text.push_str(page["text"].as_str().unwrap());
        if page["done"] == true {
            return (text, page["encoding"].as_str().unwrap().into());
        }
        at = page["next_start"].as_u64().unwrap();
    }
}
async fn archive(c: &Client, installed: &Value) -> Workspace {
    let (text, encoding) = artifact(c, &installed["native_artifact_id"]).await;
    assert_eq!(encoding, "base64");
    workspace::load_workspace(&STANDARD.decode(text).unwrap()).unwrap()
}
async fn import(c: &Client, w: &Workspace) -> Value {
    let imported = call(
        c,
        "import_workspace_archive",
        json!({"data_base64":STANDARD.encode(workspace::save_workspace(w).unwrap())}),
    )
    .await;
    reference(&imported)
}
fn reference(reply: &Value) -> Value {
    json!({"id":reply["workspace"]["id"],"revision":reply["workspace"]["revision"]})
}
async fn preview(
    c: &Client,
    w: &Value,
    source: &EventBinding,
    origin: &Origin,
    preparation: &Preparation,
) -> native::OverlapExtensionPreview {
    let result = call(c, "material_derivation", json!({"workspace":w,"task":{"PreviewOverlapExtensionPcr":{"source":source,"origin":origin,"preparation":preparation}}})).await;
    let (text, encoding) = artifact(c, &result["artifact_id"]).await;
    assert_eq!(encoding, "utf-8");
    let evidence: Value = serde_json::from_str(&text).unwrap();
    let report: native::OverlapExtensionPreview =
        serde_json::from_value(evidence["preview"].clone()).unwrap();
    report
        .validate_binding(result["context_revision"].as_str().unwrap())
        .unwrap();
    report
}
fn install_request(
    source: &EventBinding,
    origin: &Origin,
    preparation: &Preparation,
    selection: &native::OverlapExtensionSelection,
) -> Value {
    json!({"source":source,"origin":origin,"derivation":{"OverlapExtensionPcr":{"preparation":preparation,"selection":selection}},"name":"Reusable edited OE material","save_primers":true})
}
async fn install(c: &Client, w: &Value, request: Value) -> Value {
    call(
        c,
        "material_derivation",
        json!({"workspace":w,"task":{"Install":{"request":request}}}),
    )
    .await
}
fn assert_material(w: &Workspace, receipt: &native::Receipt, incomplete: bool) {
    receipt.validate_archive(w).unwrap();
    let material = native::reconstruct(w, &receipt.product).unwrap();
    assert_eq!(material.reference().sequence.as_str(), AMPLICON);
    assert_eq!(material.incomplete(), incomplete);
    assert_eq!(material.tracked().duplex.bottom_start, 0);
    let modifications = collagene_core::modified_duplex::summary(material.tracked()).unwrap();
    assert_eq!(modifications.known_modified, 0);
    assert_eq!(modifications.unknown_eligible, 0);
    let table: collagene_core::oligo_formats::OligoTable =
        serde_json::from_slice(&w.attachments[receipt.primer_attachment.as_ref().unwrap()].bytes)
            .unwrap();
    assert_eq!(table.records.len(), 6);
    assert_eq!(
        table
            .records
            .iter()
            .map(|r| r.name.as_str())
            .collect::<Vec<_>>(),
        vec![
            "New left F",
            "New left R",
            "New right F",
            "New right R",
            "New outer F",
            "New outer R"
        ]
    );
    let capsule = native::capsule(w, &receipt.product).unwrap();
    let record = native::load_record(&capsule.record_bytes).unwrap();
    let ordered = native::primer_table(&record).unwrap().unwrap();
    assert_eq!(ordered.records.len(), 6);
    assert_eq!(
        ordered.records[4].raw_ordered_sequence.as_deref(),
        Some("ATGCGAATTCAGTCCGAT")
    );
    assert_eq!(
        ordered.records[5].raw_ordered_sequence.as_deref(),
        Some("ATGCGAATTCTGACGATC")
    );
}
#[tokio::test]
async fn original_oe_material_installs_six_new_primers_reopens_and_redigests_without_a_target() {
    let (c, server) = setup().await;
    let w = original();
    let input_bytes = workspace::save_workspace(&w).unwrap();
    let selected = event(&w, 0);
    let input = import(&c, &w).await;
    let before = call(&c, "list_documents", json!({})).await["total"].clone();
    let preparation = preparation();
    let report = preview(
        &c,
        &input,
        &selected,
        &Origin::OriginalReference,
        &preparation,
    )
    .await;
    assert!(report.reusable());
    assert!(!report.incomplete);
    assert_eq!(call(&c, "list_documents", json!({})).await["total"], before);
    let request = install_request(
        &selected,
        &Origin::OriginalReference,
        &preparation,
        &report.selection().unwrap(),
    );
    for field in ["expected_preview_revision", "expected_product_revision"] {
        let mut stale = request.clone();
        stale["derivation"]["OverlapExtensionPcr"]["selection"][field] = json!("0".repeat(64));
        reject(
            &c,
            json!({"workspace":input,"task":{"Install":{"request":stale}}}),
        )
        .await;
    }
    let mut changed = request.clone();
    changed["derivation"]["OverlapExtensionPcr"]["preparation"]["outer"]["forward"]["tail"] =
        json!("GG");
    reject(
        &c,
        json!({"workspace":input,"task":{"Install":{"request":changed}}}),
    )
    .await;
    let mut invented = request.clone();
    invented["derivation"]["OverlapExtensionPcr"]["selection"]["geometry"] =
        json!({"sequence":AMPLICON});
    reject(
        &c,
        json!({"workspace":input,"task":{"Install":{"request":invented}}}),
    )
    .await;
    assert_eq!(call(&c, "list_documents", json!({})).await["total"], before);
    let installed = install(&c, &input, request).await;
    assert!(
        installed["authority"]
            .as_str()
            .unwrap()
            .contains("No independent target")
    );
    assert_eq!(installed["source_workspace"], input);
    let receipt: native::Receipt = serde_json::from_value(installed["receipt"].clone()).unwrap();
    let reopened = archive(&c, &installed).await;
    assert_eq!(reopened.documents.len(), 2);
    assert!(run_lineage::list_runs(&reopened).unwrap().is_empty());
    assert_eq!(reopened.documents[0].head_event, w.documents[0].head_event);
    assert_material(&reopened, &receipt, false);
    // A fresh server can discover/replay the archived recipe, then digest its actual tails.
    c.cancel().await.unwrap();
    server.await.unwrap();
    let (c, server) = setup().await;
    let restarted = import(&c, &reopened).await;
    let origin = Origin::Saved(MaterialSelection::Derived(receipt.product.clone()));
    let inspect = call(&c, "inspect_material_origin", json!({"workspace":restarted,"event":receipt.product,"selection":MaterialSelection::Derived(receipt.product.clone())})).await;
    assert_eq!(inspect["origin_available"], true);
    let catalog = catalog_digest::builtin_catalog().unwrap();
    let choice = catalog_digest::select(&catalog, &["ecori".into()]).unwrap();
    let digest = call(&c, "material_derivation", json!({"workspace":restarted,"task":{"PreviewDigest":{"source":receipt.product,"origin":origin,"selection":choice,"options":collagene_core::restriction_inspection::InspectionOptions::default()}}})).await;
    assert_eq!(digest["products_determined"], true);
    assert_eq!(digest["fragment_count"], 3);
    let cut = install(&c, &restarted, json!({"source":receipt.product,"origin":origin,"derivation":{"CatalogDigest":{"selection":choice,"options":collagene_core::restriction_inspection::InspectionOptions::default(),"fragment_index":2}},"name":"Literal right EcoRI tail fragment","save_primers":false})).await;
    let cut_receipt: native::Receipt = serde_json::from_value(cut["receipt"].clone()).unwrap();
    let cut_w = archive(&c, &cut).await;
    let cut_material = native::reconstruct(&cut_w, &cut_receipt.product).unwrap();
    assert_eq!(cut_material.reference().sequence.as_str(), RIGHT_CUT);
    assert_eq!(cut_material.tracked().duplex.bottom_start, 4);
    assert!(cut_receipt.primer_attachment.is_none());
    assert_eq!(workspace::save_workspace(&w).unwrap(), input_bytes);
    c.cancel().await.unwrap();
    server.await.unwrap();
}
fn saved_source(run: bool) -> (Workspace, EventBinding, Origin) {
    let w = original();
    let source = w.revisions[&w.documents[0].head].clone();
    let f = Primer::new("Historical stock F", "AGTCCGAT", "").unwrap();
    let r = Primer::new("Historical stock R", "TGACGATC", "").unwrap();
    if !run {
        let product = native::preview_pcr(
            &source,
            &native::SourceOrigin::OriginalReference,
            &f,
            &r,
            Strand::Top,
        )
        .unwrap()
        .remove(0);
        let installed = native::prepare_install(
            &w,
            &native::Request {
                source: event(&w, 0),
                origin: native::SourceOrigin::OriginalReference,
                derivation: native::Derivation::ExactPcr {
                    forward: f,
                    reverse: r,
                    template_strand: Strand::Top,
                    selection: native::PcrSelection {
                        selected_index: 0,
                        expected_preview_revision: product.molecule.revision(),
                        geometry: product.geometry,
                    },
                },
                name: "Historical actual PCR stock".into(),
                save_primers: true,
            },
        )
        .unwrap();
        let event = installed.receipt.product;
        return (
            installed.workspace,
            event.clone(),
            Origin::Saved(MaterialSelection::Derived(event)),
        );
    }
    let target = Molecule::dna(
        "Historical independently supplied stock requirement",
        ORIGINAL,
        Topology::Linear,
    )
    .unwrap();
    let w = Workspace::from_documents(
        "Historical actual run",
        vec![source.clone(), target.clone()],
        0,
    )
    .unwrap();
    let plan = DesignPlan {
        name: "Historical exact stock PCR".into(),
        target_revision: target.revision(),
        inputs: vec![PlanInput {
            id: "stock".into(),
            revision: source.revision(),
            molecule: source,
        }],
        steps: vec![Operation::Pcr {
            template: Source::Input("stock".into()),
            forward: f,
            reverse: r,
            template_strand: Strand::Top,
        }],
    };
    let installed = cloning_outputs::finalize(
        &w,
        &plan,
        event(&w, 1),
        vec![InputEventBinding {
            input_id: "stock".into(),
            event: event(&w, 0),
        }],
        &cloning_outputs::Options {
            name: "Historical run product".into(),
            show_prepared: false,
            save_primers: true,
            allow_incomplete: false,
        },
    )
    .unwrap();
    let event = installed.receipt.construct;
    let origin = Origin::Saved(MaterialSelection::Run(run_material::Selection {
        run_id: installed.receipt.run.run_id,
        step_index: 0,
        event: event.clone(),
    }));
    (installed.workspace, event, origin)
}
#[tokio::test]
async fn saved_run_and_derived_sources_require_exact_native_owners_and_keep_new_orders_separate() {
    for run in [false, true] {
        let (c, server) = setup().await;
        let (w, source, origin) = saved_source(run);
        let original_bytes = workspace::save_workspace(&w).unwrap();
        let w = workspace::load_workspace(&original_bytes).unwrap();
        let input = import(&c, &w).await;
        let prep = preparation();
        let report = preview(&c, &input, &source, &origin, &prep).await;
        assert!(report.reusable());
        let before = call(&c, "list_documents", json!({})).await["total"].clone();
        let request = install_request(&source, &origin, &prep, &report.selection().unwrap());
        reject(&c,json!({"workspace":input,"task":{"PreviewOverlapExtensionPcr":{"source":source,"origin":"OriginalReference","preparation":prep}}})).await;
        let mut alias = source.clone();
        alias.event_id = event(&w, 0).event_id;
        reject(&c,json!({"workspace":input,"task":{"PreviewOverlapExtensionPcr":{"source":alias,"origin":origin,"preparation":prep}}})).await;
        let foreign = original();
        let foreign_ref = import(&c, &foreign).await;
        reject(
            &c,
            json!({"workspace":foreign_ref,"task":{"Install":{"request":request}}}),
        )
        .await;
        assert_eq!(call(&c, "list_documents", json!({})).await["total"], before);
        let installed = install(&c, &input, request).await;
        let receipt: native::Receipt =
            serde_json::from_value(installed["receipt"].clone()).unwrap();
        let reopened = archive(&c, &installed).await;
        assert_material(&reopened, &receipt, false);
        assert_eq!(
            run_lineage::list_runs(&reopened).unwrap().len(),
            usize::from(run)
        );
        assert_eq!(
            reopened
                .documents
                .iter()
                .find(|d| d.id == source.document_id)
                .unwrap()
                .head_event
                .as_ref(),
            Some(&source.event_id)
        );
        let exported = call(&c, "export_workspace_archive", input.clone()).await;
        let (text, encoding) = artifact(&c, &exported["artifact_id"]).await;
        assert_eq!(encoding, "base64");
        assert_eq!(STANDARD.decode(text).unwrap(), original_bytes);
        c.cancel().await.unwrap();
        server.await.unwrap();
    }
}
#[tokio::test]
async fn ambiguous_oe_products_are_withheld_and_capped_extension_stays_incomplete_after_reopen() {
    let (c, server) = setup().await;
    let repeated = Workspace::from_documents(
        "Repeated source",
        vec![
            Molecule::dna(
                "Duplicated",
                &format!("{ORIGINAL}{ORIGINAL}"),
                Topology::Linear,
            )
            .unwrap(),
        ],
        0,
    )
    .unwrap();
    let input = import(&c, &repeated).await;
    let before = call(&c, "list_documents", json!({})).await["total"].clone();
    let report = preview(
        &c,
        &input,
        &event(&repeated, 0),
        &Origin::OriginalReference,
        &preparation(),
    )
    .await;
    assert!(!report.reusable());
    assert!(report.selection().is_err());
    assert_eq!(call(&c, "list_documents", json!({})).await["total"], before);
    let w = original();
    let input = import(&c, &w).await;
    let source = event(&w, 0);
    let mut ambiguous = preparation();
    ambiguous.outer.forward = Primer::new("Ambiguous F", "T", "").unwrap();
    ambiguous.outer.reverse = Primer::new("Ambiguous R", "A", "").unwrap();
    let report = preview(&c, &input, &source, &Origin::OriginalReference, &ambiguous).await;
    assert!(!report.reusable());
    assert!(report.selection().is_err());
    let mut conditional = preparation();
    conditional.extension.max_comparisons = 18;
    let report = preview(
        &c,
        &input,
        &source,
        &Origin::OriginalReference,
        &conditional,
    )
    .await;
    assert!(report.reusable());
    assert!(report.incomplete);
    let installed = install(
        &c,
        &input,
        install_request(
            &source,
            &Origin::OriginalReference,
            &conditional,
            &report.selection().unwrap(),
        ),
    )
    .await;
    let receipt: native::Receipt = serde_json::from_value(installed["receipt"].clone()).unwrap();
    assert!(receipt.incomplete);
    assert_material(&archive(&c, &installed).await, &receipt, true);
    c.cancel().await.unwrap();
    server.await.unwrap();
}
