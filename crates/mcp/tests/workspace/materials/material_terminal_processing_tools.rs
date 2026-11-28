//! Saved terminal chemistry remains native material, with independent ligation intent.
use base64::{Engine as _, engine::general_purpose::STANDARD};
use collagene_core::{
    chemistry::{FivePrime, Polarity, Strand, ThreePrime},
    cloning::{
        assembly::recipe::{FragmentSpec, Preparation, RunPreparation},
        ligation::recipe as ligation_recipe,
        overlap::extension::{
            self as overlap_extension,
            pcr::{PrimerPair, Request as OePreparation},
        },
    },
    documents::workspace::{self, Workspace},
    material::{
        derivation as native,
        lineage::{self as run_lineage, EventBinding},
        run::{self as run_material, MaterialSelection},
    },
    planning::plan::DesignPlan,
    sequence::analysis::Primer,
    verification::evidence::PhysicalProcessingEvidence,
};
use collagene_mcp::server::CollageneServer;
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
fn binding(w: &Workspace, i: usize) -> EventBinding {
    let doc = &w.documents[i];
    EventBinding {
        document_id: doc.id.clone(),
        event_id: doc.head_event.clone().unwrap(),
        revision: doc.head.clone(),
    }
}
fn saved(source: &EventBinding) -> Value {
    json!({"Saved":MaterialSelection::Derived(source.clone())})
}
async fn pcr(c: &Client, w: &Value, source: &EventBinding, f: &str, r: &str) -> Value {
    let f = Primer::new("Historical stock F", f, "").unwrap();
    let r = Primer::new("Historical stock R", r, "").unwrap();
    let preview=call(c,"material_derivation",json!({"workspace":w,"task":{"PreviewPcr":{"source":source,"origin":"OriginalReference","forward":f,"reverse":r,"template_strand":"Top"}}})).await;
    assert_eq!(preview["count"], 1);
    call(c,"material_derivation",json!({"workspace":w,"task":{"Install":{"request":{"sources":[{"role":"template","event":source,"origin":"OriginalReference"}],"derivation":{"ExactPcr":{"forward":f,"reverse":r,"template_strand":"Top","selection":{"selected_index":0,"expected_preview_revision":preview["previews"][0]["preview_revision"]}}},"name":"Actual native PCR stock","save_primers":true}}}})).await
}
async fn terminal(c: &Client, w: &Value, source: &EventBinding, action: Value) -> Value {
    let preview=call(c,"material_derivation",json!({"workspace":w,"task":{"PreviewTerminalProcessing":{"source":source,"origin":saved(source),"action":action}}})).await;
    let (text, encoding) = artifact(c, &preview["artifact_id"]).await;
    assert_eq!(encoding, "utf-8");
    let evidence: Value = serde_json::from_str(&text).unwrap();
    let report: native::TerminalProcessingPreview =
        serde_json::from_value(evidence["preview"].clone()).unwrap();
    report
        .validate_binding(preview["context_revision"].as_str().unwrap())
        .unwrap();
    let installed=call(c,"material_derivation",json!({"workspace":w,"task":{"Install":{"request":{"sources":[{"role":"template","event":source,"origin":saved(source)}],"derivation":{"TerminalProcessing":{"action":action}},"name":"Explicit terminal preparation","save_primers":false}}}})).await;
    let receipt: native::Receipt = serde_json::from_value(installed["receipt"].clone()).unwrap();
    assert!(receipt.primer_attachment.is_none());
    let archive = archive(c, &installed).await;
    receipt.validate_archive(&archive).unwrap();
    let record = native::load_record(
        &native::capsule(&archive, &receipt.product)
            .unwrap()
            .record_bytes,
    )
    .unwrap();
    assert!(native::primer_table(&record).unwrap().is_none());
    installed
}
fn receipt_of(value: &Value) -> native::Receipt {
    serde_json::from_value(value["receipt"].clone()).unwrap()
}
fn ligation(w: &Workspace, v: &EventBinding, i: &EventBinding) -> ligation_recipe::Recipe {
    let row = |id: &str, event: &EventBinding| ligation_recipe::Fragment {
        source: FragmentSpec {
            generation: None,
            role: id.into(),
            input_id: id.into(),
            preparation: Preparation::FromDerived {
                origin: native::capsule(w, event).unwrap(),
                preparation: RunPreparation::WholeLinear,
            },
            reverse: false,
        },
        treatment: None,
    };
    ligation_recipe::Recipe {
        name: "Independently intended TA circle from retained native tails".into(),
        fragments: vec![row("vector", v), row("insert", i)],
        circular: true,
        allow_nicked: false,
    }
}
async fn propose(
    c: &Client,
    w: &Workspace,
    v: &EventBinding,
    i: &EventBinding,
    target: &Value,
) -> Value {
    call(c,"propose_ligation_assembly",json!({"inputs":[{"id":"vector","document":{"id":v.revision,"revision":v.revision}},{"id":"insert","document":{"id":i.revision,"revision":i.revision}}],"target":target,"recipe":ligation(w,v,i)})).await
}
#[tokio::test]
async fn saved_pcr_tails_and_phosphates_reopen_into_ordinary_ligation_with_extractable_alternatives()
 {
    let (c, server) = setup().await;
    let source = Workspace::from_documents(
        "Two references without a target",
        vec![
            Molecule::dna("Vector reference", "AAACCC", Topology::Linear).unwrap(),
            Molecule::dna("Insert reference", "GGGTTT", Topology::Linear).unwrap(),
        ],
        0,
    )
    .unwrap();
    let source_bytes = workspace::save_workspace(&source).unwrap();
    let original = import(&c, &source).await;
    let vector = pcr(&c, &original, &binding(&source, 0), "AAA", "GGG").await;
    let insert = pcr(&c, &reference(&vector), &binding(&source, 1), "GGG", "AAA").await;
    let vp = receipt_of(&vector).product;
    let ip = receipt_of(&insert).product;
    let vtail = terminal(
        &c,
        &reference(&insert),
        &vp,
        json!({"Tail3Prime":{"base":"T"}}),
    )
    .await;
    let itail = terminal(
        &c,
        &reference(&vtail),
        &ip,
        json!({"Tail3Prime":{"base":"A"}}),
    )
    .await;
    let vt = receipt_of(&vtail).product;
    let it = receipt_of(&itail).product;
    let tail_w = archive(&c, &itail).await;
    for (event, dna, base) in [(&vt, "AAACCCT", "T"), (&it, "GGGTTTA", "A")] {
        let fresh = native::reconstruct(&tail_w, event).unwrap();
        assert_eq!(fresh.reference().sequence.as_str(), dna);
        let summary = fresh.tracked().duplex.summary().unwrap();
        let (left, right) = summary.ends.unwrap();
        for end in [left, right] {
            assert_eq!(end.polarity, Polarity::ThreePrime);
            assert_eq!(end.overhang, base);
            assert_eq!(end.chemistry.five_prime, FivePrime::Hydroxyl);
            assert_eq!(end.chemistry.three_prime, ThreePrime::Hydroxyl);
        }
    }
    let target = Molecule::dna(
        "Independent literal TA target",
        "AAACCCTGGGTTTA",
        Topology::Circular,
    )
    .unwrap();
    let imported = call(
        &c,
        "import_document",
        json!({"text":collagene_core::formats::sequence::export_genbank(&target).unwrap()}),
    )
    .await;
    let row = &imported["documents"][0];
    let target_ref = json!({"id":row["id"],"revision":row["revision"]});
    let unsealed = propose(&c, &tail_w, &vt, &it, &target_ref).await;
    assert_eq!(
        unsealed["verification"]["status"], "Invalid",
        "tailing must never supply missing5′phosphates"
    );
    let before = call(&c, "list_documents", json!({})).await["total"].clone();
    let mut stale = reference(&itail);
    stale["revision"] = json!("0".repeat(64));
    reject(&c,json!({"workspace":stale,"task":{"PreviewTerminalProcessing":{"source":vt,"origin":saved(&vt),"action":{"TreatEnds":{"treatment":"Phosphorylate"}}}}})).await;
    reject(&c,json!({"workspace":original,"task":{"PreviewTerminalProcessing":{"source":vt,"origin":saved(&vt),"action":{"TreatEnds":{"treatment":"Phosphorylate"}}}}})).await;
    reject(&c,json!({"workspace":reference(&itail),"task":{"Install":{"request":{"sources":[{"role":"template","event":vt,"origin":saved(&vt)}],"derivation":{"TerminalProcessing":{"action":{"TreatEnds":{"treatment":"Phosphorylate"}}}},"name":"No oligos are synthesized","save_primers":true}}}})).await;
    assert_eq!(call(&c, "list_documents", json!({})).await["total"], before);
    let vphos = terminal(
        &c,
        &reference(&itail),
        &vt,
        json!({"TreatEnds":{"treatment":"Phosphorylate"}}),
    )
    .await;
    let iphos = terminal(
        &c,
        &reference(&vphos),
        &it,
        json!({"TreatEnds":{"treatment":"Phosphorylate"}}),
    )
    .await;
    let v = receipt_of(&vphos).product;
    let i = receipt_of(&iphos).product;
    let reopened = archive(&c, &iphos).await;
    assert!(run_lineage::list_runs(&reopened).unwrap().is_empty());
    assert_eq!(
        reopened.documents[0].head_event,
        source.documents[0].head_event
    );
    assert_eq!(reopened.documents[1].head, source.documents[1].head);
    c.cancel().await.unwrap();
    server.await.unwrap();
    let (c, server) = setup().await;
    let _ = import(&c, &reopened).await;
    let imported = call(
        &c,
        "import_document",
        json!({"text":collagene_core::formats::sequence::export_genbank(&target).unwrap()}),
    )
    .await;
    let row = &imported["documents"][0];
    let target_ref = json!({"id":row["id"],"revision":row["revision"]});
    let result = propose(&c, &reopened, &v, &i, &target_ref).await;
    assert_eq!(
        result["verification"]["status"],
        "ConsistentUnderAssumptions"
    );
    assert_eq!(result["product_physical"]["nick_count"], 0);
    let (text, encoding) = artifact(&c, &result["artifact_id"]).await;
    assert_eq!(encoding, "utf-8");
    let evidence: Value = serde_json::from_str(&text).unwrap();
    let report: collagene_core::verification::evidence::VerificationReport =
        serde_json::from_value(evidence["verification"].clone()).unwrap();
    let outcome = report
        .physical_processing_reports
        .iter()
        .find_map(|step| match &step.evidence {
            PhysicalProcessingEvidence::LigationOutcomes(e) => Some(&e.report),
            _ => None,
        })
        .unwrap();
    assert!(
        format!(
            "{}{}",
            outcome.intended.sequence.as_str(),
            outcome.intended.sequence.as_str()
        )
        .contains("AAACCCTGGGTTTA")
    );
    assert!(!outcome.alternatives.is_empty());
    assert!(outcome.alternatives.iter().any(|p| {
        format!("{}{}", p.sequence.as_str(), p.sequence.as_str()).contains("AAACCCTAAACCCA")
    }));
    let html = call(
        &c,
        "export_protocol",
        json!({"plan_revision":result["plan_revision"]}),
    )
    .await;
    let (text, encoding) = artifact(&c, &html["artifact_id"]).await;
    assert_eq!(encoding, "utf-8");
    assert!(text.contains("Full extractable outcome FASTA"));
    assert!(text.contains("Alternative"));
    let plan: DesignPlan = serde_json::from_value(evidence["compiled"]["plan"].clone()).unwrap();
    assert_eq!(
        plan.target_revision,
        target_ref["revision"].as_str().unwrap()
    );
    assert_eq!(plan.inputs.len(), 2);
    assert_eq!(workspace::save_workspace(&source).unwrap(), source_bytes);
    c.cancel().await.unwrap();
    server.await.unwrap();
}
fn conditional_oe_stock() -> (Workspace, EventBinding) {
    let source = Molecule::dna(
        "Original literal stock",
        "AGTCCGATGCTACGTAGGACTTCAGACCTGTAGATCGTCA",
        Topology::Linear,
    )
    .unwrap();
    let w =
        Workspace::from_documents("Conditional source-only OE recipe", vec![source], 0).unwrap();
    let pair = |f: &str, ft: &str, r: &str, rt: &str| PrimerPair {
        forward: Primer::new("Historical OE F", f, ft).unwrap(),
        reverse: Primer::new("Historical OE R", r, rt).unwrap(),
        template_strand: Strand::Top,
    };
    let prep = OePreparation {
        left: pair("AGTCCGAT", "", "TACGTAGC", "GGTCTGAAAA"),
        right: pair("TTCAGACC", "GCTACGTATT", "TGACGATC", ""),
        outer: pair("AGTCCGAT", "ATGCGAATTC", "TGACGATC", "ATGCGAATTC"),
        extension: overlap_extension::Options {
            overlap: 18,
            min_paired_bases: 8,
            max_comparisons: 18,
            ..Default::default()
        },
    };
    let origin = native::SourceOrigin::OriginalReference;
    let source = binding(&w, 0);
    let preview = native::preview_overlap_extension(&w, &source, &origin, &prep).unwrap();
    assert!(preview.incomplete);
    assert!(preview.reusable());
    let result = native::prepare_install(
        &w,
        &native::Request {
            sources: vec![collagene_core::material::derivation::BoundSource {
                role: collagene_core::material::derivation::SourceRole::Template,
                event: source,
                origin: origin,
            }],
            derivation: native::Derivation::OverlapExtensionPcr {
                preparation: prep,
                selection: preview.selection().unwrap(),
            },
            name: "Conditional freshly prepared actual stock".into(),
            save_primers: false,
        },
    )
    .unwrap();
    (result.workspace, result.receipt.product)
}
#[tokio::test]
async fn terminal_processing_preserves_conditional_origins_and_does_not_repair_unknown_ends() {
    let (c, server) = setup().await;
    let original = Workspace::from_documents(
        "Unknown terminal chemistry",
        vec![
            Molecule::dna(
                "Unmeasured original reference",
                "AGTCCGAT",
                Topology::Linear,
            )
            .unwrap(),
        ],
        0,
    )
    .unwrap();
    let original_ref = import(&c, &original).await;
    let source = binding(&original, 0);
    let before = call(&c, "list_documents", json!({})).await["total"].clone();
    reject(&c,json!({"workspace":original_ref,"task":{"PreviewTerminalProcessing":{"source":source,"origin":"OriginalReference","action":{"Tail3Prime":{"base":"A"}}}}})).await;
    assert_eq!(call(&c, "list_documents", json!({})).await["total"], before);
    let installed=call(&c,"material_derivation",json!({"workspace":original_ref,"task":{"Install":{"request":{"sources":[{"role":"template","event":source,"origin":"OriginalReference"}],"derivation":{"TerminalProcessing":{"action":{"TreatEnds":{"treatment":"Phosphorylate"}}}},"name":"Explicit5P only","save_primers":false}}}})).await;
    let receipt = receipt_of(&installed);
    assert!(receipt.primer_attachment.is_none());
    let original_phos = archive(&c, &installed).await;
    let fresh = native::reconstruct(&original_phos, &receipt.product).unwrap();
    assert_eq!(fresh.reference().sequence.as_str(), "AGTCCGAT");
    assert_eq!(fresh.tracked().duplex.left.five_prime, FivePrime::Phosphate);
    assert_eq!(
        fresh.tracked().duplex.right.five_prime,
        FivePrime::Phosphate
    );
    assert_eq!(fresh.tracked().duplex.left.three_prime, ThreePrime::Unknown);
    assert_eq!(
        fresh.tracked().duplex.right.three_prime,
        ThreePrime::Unknown
    );
    reject(&c,json!({"workspace":reference(&installed),"task":{"PreviewTerminalProcessing":{"source":receipt.product,"origin":saved(&receipt.product),"action":{"Tail3Prime":{"base":"A"}}}}})).await;
    let (conditional, source) = conditional_oe_stock();
    let original_bytes = workspace::save_workspace(&conditional).unwrap();
    let wref = import(&c, &conditional).await;
    let installed = terminal(
        &c,
        &wref,
        &source,
        json!({"TreatEnds":{"treatment":"Phosphorylate"}}),
    )
    .await;
    let receipt = receipt_of(&installed);
    assert!(receipt.incomplete);
    let next = archive(&c, &installed).await;
    let prior = native::reconstruct(&conditional, &source).unwrap();
    let fresh = native::reconstruct(&next, &receipt.product).unwrap();
    assert!(fresh.incomplete());
    assert_eq!(fresh.reference().sequence, prior.reference().sequence);
    assert_eq!(
        fresh.tracked().duplex.bottom_start,
        prior.tracked().duplex.bottom_start
    );
    assert_eq!(fresh.tracked().duplex.nicks, prior.tracked().duplex.nicks);
    assert_eq!(
        fresh.tracked().top_modifications,
        prior.tracked().top_modifications
    );
    assert_eq!(
        fresh.tracked().bottom_modifications,
        prior.tracked().bottom_modifications
    );
    assert_eq!(
        workspace::save_workspace(&conditional).unwrap(),
        original_bytes
    );
    c.cancel().await.unwrap();
    server.await.unwrap();
    let (c, server) = setup().await;
    let wref = import(&c, &next).await;
    let inspected=call(&c,"inspect_material_origin",json!({"workspace":wref,"event":receipt.product,"selection":MaterialSelection::Derived(receipt.product.clone())})).await;
    assert_eq!(inspected["origin_available"], true);
    let (text, encoding) = artifact(&c, &inspected["artifact_id"]).await;
    assert_eq!(encoding, "utf-8");
    let inspected: run_material::MaterialInspection = serde_json::from_str(&text).unwrap();
    let origin = inspected.selected.unwrap();
    let molecule = &next.revisions[&receipt.product.revision];
    assert!(origin.replay(molecule).unwrap().incomplete());
    c.cancel().await.unwrap();
    server.await.unwrap();
}
