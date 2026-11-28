//! AI inspects the actual OE amplicon and reuses its catalog choice in ordinary cloning.
use collagene_core::{
    chemistry::EndRepairMethod,
    cloning::{
        assembly::recipe::{FragmentSpec, Preparation},
        ligation::recipe as ligation_recipe,
        overlap::extension::{
            self as overlap_extension,
            digest::{Preview, Request},
            pcr::{PrimerPair, Request as Extension},
        },
        restriction::catalog as catalog_digest,
    },
    sequence::analysis::Primer,
};
use collagene_mcp::server::CollageneServer;
use rmcp::{RoleClient, ServiceExt, model::CallToolRequestParams, service::RunningService};
use serde_json::{Value, json};
type Client = RunningService<RoleClient, ()>;
const ORIGINAL: &str = "AGTCCGATGCTACGTAGGACTTCAGACCTGTAGATCGTCA";
const AMPLICON: &str = "ATGCGAATTCAGTCCGATGCTACGTATTTTCAGACCTGTAGATCGTCAGAATTCGCAT";
const TARGET: &str = "AATTCAGTCCGATGCTACGTATTTTCAGACCTGTAGATCGTCAGAATT";
fn pair(f: &str, ft: &str, r: &str, rt: &str) -> PrimerPair {
    PrimerPair {
        forward: Primer::new("F", f, ft).unwrap(),
        reverse: Primer::new("R", r, rt).unwrap(),
        template_strand: collagene_core::chemistry::Strand::Top,
    }
}
fn extension() -> Extension {
    Extension {
        left: pair("AGTCCGAT", "", "TACGTAGC", "GGTCTGAAAA"),
        right: pair("TTCAGACC", "GCTACGTATT", "TGACGATC", ""),
        outer: pair("AGTCCGAT", "ATGCGAATTC", "TGACGATC", "ATGCGAATTC"),
        extension: overlap_extension::Options {
            overlap: 18,
            min_paired_bases: 8,
            ..Default::default()
        },
    }
}
fn selection() -> catalog_digest::CatalogSelection {
    let catalog = catalog_digest::builtin_catalog().unwrap();
    let id = catalog
        .edition()
        .enzymes
        .iter()
        .find(|e| e.name == "EcoRI")
        .unwrap()
        .id
        .clone();
    catalog_digest::select(&catalog, &[id]).unwrap()
}
fn request() -> Request {
    Request {
        source_role: "insert".into(),
        preparation: extension(),
        selection: selection(),
        options: Default::default(),
    }
}
fn recipe(request: Extension) -> ligation_recipe::Recipe {
    ligation_recipe::Recipe {
        standard: None,
        name: "Physical OE, selected EcoRI fragment and explicit fill".into(),
        fragments: vec![ligation_recipe::Fragment {
            source: FragmentSpec {
                generation: None,
                role: "insert".into(),
                input_id: "source".into(),
                preparation: Preparation::OverlapExtensionPcrCatalogDigest {
                    request,
                    conditions: None,
                    selection: selection(),
                    fragment_index: 1,
                    repair: Some(EndRepairMethod::FillFivePrimeOverhangs),
                },
                reverse: false,
            },
            treatment: None,
        }],
        circular: false,
        allow_nicked: false,
    }
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
async fn artifact(c: &Client, id: &Value) -> Value {
    let (mut start, mut result) = (0, String::new());
    loop {
        let part = call(
            c,
            "read_artifact",
            json!({"artifact_id":id,"start":start,"max_bytes":16384}),
        )
        .await;
        result.push_str(part["text"].as_str().unwrap());
        if part["done"] == true {
            return serde_json::from_str(&result).unwrap();
        }
        start = part["next_start"].as_u64().unwrap();
    }
}
fn named(imported: &Value, name: &str) -> Value {
    let row = imported["documents"]
        .as_array()
        .unwrap()
        .iter()
        .find(|r| r["name"] == name)
        .unwrap();
    json!({"id":row["id"],"revision":row["revision"]})
}
#[tokio::test]
async fn actual_four_step_digest_and_compound_recipe_preserve_stock_and_conditional_evidence() {
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
    let tools = client.list_all_tools().await.unwrap();
    assert_eq!(tools.len(), collagene_mcp::CURRENT_TOOL_COUNT);
    assert!(
        tools
            .iter()
            .any(|tool| tool.name == "inspect_overlap_extension_digest")
    );
    let imported = call(
        &client,
        "import_document",
        json!({"text":format!(">original\n{ORIGINAL}\n>independent_filled_fragment\n{TARGET}\n")}),
    )
    .await;
    let original = named(&imported, "original");
    let target = named(&imported, "independent_filled_fragment");
    let before = call(&client, "list_documents", json!({})).await;
    let inspected = call(
        &client,
        "inspect_overlap_extension_digest",
        json!({"document":original,"origin":null,"workspace":null,"request":request(),"max_fragments":2}),
    )
    .await;
    assert_eq!(inspected["amplicon"]["length"], 58);
    assert_eq!(inspected["fragment_count"], 3);
    assert_eq!(inspected["fragments"].as_array().unwrap().len(), 2);
    assert_eq!(inspected["next_fragment_start"], 2);
    assert_eq!(inspected["fragments_done"], false);
    assert_eq!(inspected["selectable"], true);
    assert_eq!(inspected["incomplete"], false);
    assert_eq!(
        inspected["pcr_stages"]
            .as_array()
            .unwrap()
            .iter()
            .map(|s| s["product_length"].as_u64().unwrap())
            .collect::<Vec<_>>(),
        vec![26, 30, 58]
    );
    let preview: Preview =
        serde_json::from_value(artifact(&client, &inspected["artifact_id"]).await).unwrap();
    preview
        .validate_binding(inspected["context_revision"].as_str().unwrap())
        .unwrap();
    assert_eq!(
        preview.amplicon.as_ref().unwrap().sequence.as_str(),
        AMPLICON
    );
    assert_eq!(
        preview
            .inspection
            .as_ref()
            .unwrap()
            .fragments
            .iter()
            .map(|f| f.duplex.top_length)
            .collect::<Vec<_>>(),
        vec![5, 44, 9]
    );
    assert_eq!(
        call(&client, "list_documents", json!({})).await,
        before,
        "digest inspection must not install a target or a product document"
    );
    for oligo in 0..6 {
        let mut changed = request();
        let pair = match oligo / 2 {
            0 => &mut changed.preparation.left,
            1 => &mut changed.preparation.right,
            _ => &mut changed.preparation.outer,
        };
        let primer = if oligo % 2 == 0 {
            &mut pair.forward
        } else {
            &mut pair.reverse
        };
        primer.annealing = collagene_model::Sequence::dna("AAAAAAAA").unwrap();
        let rejected = call(
            &client,
            "inspect_overlap_extension_digest",
            json!({"document":original,"origin":null,"workspace":null,"request":changed}),
        )
        .await;
        assert_eq!(
            rejected["selectable"], false,
            "changed oligo {oligo}: {rejected}"
        );
        assert!(
            rejected["amplicon"].is_null(),
            "changed oligo {oligo}: {rejected}"
        );
        assert_ne!(rejected["context_revision"], inspected["context_revision"]);
    }
    let proposed=call(&client,"propose_ligation_assembly",json!({"inputs":[{"id":"source","document":original}],"target":target,"recipe":recipe(extension())})).await;
    assert_eq!(proposed["verified"], false);
    assert_eq!(
        proposed["verification"]["status"], "ConsistentUnderAssumptions",
        "{proposed}"
    );
    let evidence = artifact(&client, &proposed["artifact_id"]).await;
    let final_product: collagene_model::Molecule = serde_json::from_value(
        evidence["verification"]["steps"]
            .as_array()
            .unwrap()
            .last()
            .unwrap()["product"]
            .clone(),
    )
    .unwrap();
    assert_eq!(final_product.sequence.as_str(), TARGET);
    assert_eq!(
        call(
            &client,
            "verify_plan",
            json!({"plan_revision":proposed["plan_revision"]})
        )
        .await["status"],
        "ConsistentUnderAssumptions"
    );
    let wrong=call(&client,"propose_ligation_assembly",json!({"inputs":[{"id":"source","document":original}],"target":original,"recipe":recipe(extension())})).await;
    assert_eq!(wrong["verification"]["status"], "Invalid");
    let mut conditional = request();
    conditional.preparation.extension.max_comparisons = 18;
    let inspected = call(
        &client,
        "inspect_overlap_extension_digest",
        json!({"document":original,"origin":null,"workspace":null,"request":conditional}),
    )
    .await;
    assert_eq!(inspected["selectable"], true);
    assert_eq!(inspected["incomplete"], true);
    assert_eq!(inspected["extension"]["complete"], false);
    let mut conditional = extension();
    conditional.extension.max_comparisons = 18;
    let proposed=call(&client,"propose_ligation_assembly",json!({"inputs":[{"id":"source","document":original}],"target":target,"recipe":recipe(conditional)})).await;
    assert_eq!(
        proposed["verification"]["status"], "Incomplete",
        "{proposed}"
    );
    assert_eq!(
        call(
            &client,
            "read_sequence_region",
            json!({"document":original,"start":0,"end":40})
        )
        .await["sequence"],
        ORIGINAL
    );
    drop(client);
    server.await.unwrap();
}
