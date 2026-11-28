use collagene_core::{analysis::Primer, golden_gate::GoldenGateOptions};
use collagene_mcp::CollageneServer;
use collagene_model::{Molecule, Topology};
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
fn doc(v: &Value) -> Value {
    json!({"id":v["documents"][0]["id"],"revision":v["documents"][0]["revision"]})
}
async fn import(c: &Client, m: Molecule) -> Value {
    doc(&call(
        c,
        "import_document",
        json!({"text":collagene_core::formats::export_genbank(&m).unwrap()}),
    )
    .await)
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
async fn whole_reaction_tools_retain_vector_stuffer_amplicon_flanks_and_uncut_carryover() {
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
    // Independent construct precedes any preparation or reaction inspection.
    let target = import(
        &c,
        Molecule::dna(
            "desired",
            "ATGACCCCGGGGTTTCCAGAATGAGTTAGCGACCTTGGAA",
            Topology::Circular,
        )
        .unwrap(),
    )
    .await;
    let vector = import(
        &c,
        Molecule::dna(
            "whole vector",
            "ATGCGGTCTCAATGACCCCGGGGTTTCCAGAATGATGAGACCGCAT",
            Topology::Circular,
        )
        .unwrap(),
    )
    .await;
    let insert = import(
        &c,
        Molecule::dna("insert template", "ATGAGTTAGCGACCTTGGAA", Topology::Linear).unwrap(),
    )
    .await;
    let carry = import(
        &c,
        Molecule::dna("uncut carryover", "ACGTACGTACGT", Topology::Circular).unwrap(),
    )
    .await;
    let before = call(&c, "export_document", vector.clone()).await;
    let inputs = json!([{"id":"vector","document":vector},{"id":"insert","document":insert},{"id":"carry","document":carry}]);
    let request = json!({"name":"declared tube","rows":[{"generation":{"dam":"Unmodified","dcm":"Unmodified","cpg":"Unmodified"},"role":"vector","input_id":"vector","substrate":{"WholeSource":{"origin":null}}},{"generation":{"dam":"Unmodified","dcm":"Unmodified","cpg":"Unmodified"},"role":"amplicon","input_id":"insert","substrate":{"Prepared":{"preparation":{"ExactPcr":{"template_strand":"Top","forward":Primer::new("manualF","ATGAGTTA","ATGCGGTCTCA").unwrap(),"reverse":Primer::new("manualR","TTCCAAGG","ATGCGGTCTCATCAT").unwrap()}},"treatment":null}}},{"generation":{"dam":"Unmodified","dcm":"Unmodified","cpg":"Unmodified"},"role":"uncut","input_id":"carry","substrate":{"WholeSource":{"origin":null}}}],"options":GoldenGateOptions{max_fragments:3,..Default::default()}});
    let args = json!({"origin_workspaces":[],"inputs":inputs,"target":target,"request":request});
    let inspected = call(&c, "inspect_golden_gate_reaction", args.clone()).await;
    assert_eq!(inspected["state"], "advisory");
    assert_eq!(inspected["verified"], false);
    assert_eq!(inspected["carryover_sources"], json!([2]));
    assert_eq!(inspected["pool_available"], true);
    assert_eq!(inspected["products_determined"], true);
    assert_eq!(
        inspected["profile_revision"],
        collagene_core::modified_duplex::profile_revision()
    );
    let mut unknown = args.clone();
    unknown["request"]["rows"][0]["generation"] = Value::Null;
    let withheld = call(&c, "inspect_golden_gate_reaction", unknown).await;
    assert_eq!(withheld["products_determined"], false);
    assert_eq!(withheld["candidate_count"], 0);
    assert_eq!(withheld["enumeration_complete"], false);

    let evidence = full(&c, &inspected["artifact_id"]).await;
    assert_eq!(evidence["kind"], "golden-gate-reaction-v1");
    assert_eq!(
        evidence["report"]["pool"]["fragments"]
            .as_array()
            .unwrap()
            .len(),
        6
    );
    assert_eq!(
        evidence["report"]["probe_plan"]["steps"]
            .as_array()
            .unwrap()
            .len(),
        5
    );
    let a = inspected["candidates"]
        .as_array()
        .unwrap()
        .iter()
        .find(|a| a["target_matches"] == true && a["restriction_site_free"] == true)
        .unwrap();
    let selection = json!({"context_revision":inspected["context_revision"],"candidate_index":a["candidate_index"],"directed_sequence_key":a["directed_sequence_key"],"physical_sequence_key":a["physical_sequence_key"]});
    let selected = call(
        &c,
        "select_golden_gate_reaction_candidate",
        json!({"artifact_id":inspected["artifact_id"],"selection":selection}),
    )
    .await;
    assert_eq!(selected["verified"], false);
    assert_eq!(selected["artifact_id"], selected["candidate_artifact_id"]);
    assert_eq!(
        call(
            &c,
            "verify_plan",
            json!({"plan_revision":selected["plan_revision"]})
        )
        .await["status"],
        "ConsistentUnderAssumptions"
    );
    assert_eq!(
        call(&c, "inspect_golden_gate_reaction", args.clone()).await["artifact_id"],
        inspected["artifact_id"]
    );
    let mut capped = args;
    capped["request"]["options"]["max_states"] = 1.into();
    capped["request"]["options"]["min_fragments"] = 2.into();
    let partial = call(&c, "inspect_golden_gate_reaction", capped).await;
    assert_eq!(partial["enumeration_complete"], false);
    assert_eq!(partial["pool_available"], true);
    assert_eq!(partial["candidate_count"], 0);
    let mut bad = selection;
    bad["physical_sequence_key"] = "0".repeat(64).into();
    let result = c
        .call_tool(
            CallToolRequestParams::new("select_golden_gate_reaction_candidate").with_arguments(
                json!({"artifact_id":inspected["artifact_id"],"selection":bad})
                    .as_object()
                    .unwrap()
                    .clone(),
            ),
        )
        .await;
    assert!(result.is_err() || result.unwrap().is_error == Some(true));
    assert_eq!(call(&c, "export_document", vector).await, before);
    c.cancel().await.unwrap();
    server.await.unwrap();
}

#[tokio::test]
async fn saved_circular_whole_source_requires_exact_retained_workspace_and_cannot_be_redeclared() {
    use base64::{Engine as _, engine::general_purpose::STANDARD};
    use collagene_core::{
        chemistry::EndTreatment,
        cloning_outputs,
        run_lineage::{EventBinding, InputEventBinding},
        run_material::{self, MaterialOrigin},
        verification::{DesignPlan, Operation, OrientedSource, PlanInput, Source},
        workspace::{self, Workspace},
    };
    const STOCK: &str = "AGTCCGATGAATTCCGTAGGACTTCAGACCTGTAGATCGTCA";
    fn event(w: &Workspace, index: usize) -> EventBinding {
        let d = &w.documents[index];
        EventBinding {
            document_id: d.id.clone(),
            event_id: d.head_event.clone().unwrap(),
            revision: d.head.clone(),
        }
    }
    let stock = Molecule::dna("Original literal", STOCK, Topology::Linear).unwrap();
    let target = Molecule::dna("Independently declared circle", STOCK, Topology::Circular).unwrap();
    let w = Workspace::from_documents(
        "Actual PCR and ligation",
        vec![stock.clone(), target.clone()],
        1,
    )
    .unwrap();
    let plan = DesignPlan {
        name: "Explicit PCR, phosphate and circle closure".into(),
        target_revision: target.revision(),
        inputs: vec![PlanInput {
            id: "stock".into(),
            revision: stock.revision(),
            molecule: stock,
        }],
        steps: vec![
            Operation::Pcr {
                template: Source::Input("stock".into()),
                forward: Primer::new("F", "AGTCCGAT", "").unwrap(),
                reverse: Primer::new("R", "TGACGATC", "").unwrap(),
                template_strand: collagene_core::chemistry::Strand::Top,
            },
            Operation::TreatEnds {
                template: Source::Step(0),
                treatment: EndTreatment::Phosphorylate,
            },
            Operation::Ligate {
                fragments: vec![OrientedSource {
                    source: Source::Step(1),
                    reverse: false,
                }],
                circular: true,
                allow_nicked: false,
                options: Default::default(),
            },
        ],
    };
    let saved = cloning_outputs::finalize(
        &w,
        &plan,
        event(&w, 1),
        vec![InputEventBinding {
            input_id: "stock".into(),
            event: event(&w, 0),
        }],
        &cloning_outputs::Options {
            name: "Saved physical circle".into(),
            show_prepared: false,
            save_primers: false,
            allow_incomplete: false,
        },
    )
    .unwrap();
    let selection = run_material::Selection {
        run_id: saved.receipt.run.run_id.clone(),
        step_index: 2,
        event: saved.receipt.construct.clone(),
    };
    let origin = MaterialOrigin::Run(run_material::capsule(&saved.workspace, &selection).unwrap());
    let original_archive = workspace::save_workspace(&saved.workspace).unwrap();
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
        json!({"data_base64":STANDARD.encode(&original_archive)}),
    )
    .await;
    let workspace_ref =
        json!({"id":imported["workspace"]["id"],"revision":imported["workspace"]["revision"]});
    let wrong_target = import(
        &c,
        Molecule::dna(
            "Independent unrelated requirement",
            "AAAACCCC",
            Topology::Circular,
        )
        .unwrap(),
    )
    .await;
    let mut args = json!({"inputs":[{"id":"circle","document":{"id":selection.event.revision,"revision":selection.event.revision}}],"target":wrong_target,
        "origin_workspaces":[{"role":"whole","workspace":workspace_ref}],
        "request":{"name":"Saved circle remains a whole circular substrate","rows":[{"role":"whole","input_id":"circle","generation":null,"substrate":{"WholeSource":{"origin":origin}}}],"options":GoldenGateOptions{max_fragments:1,..Default::default()}}});
    let report = call(&c, "inspect_golden_gate_reaction", args.clone()).await;
    assert_eq!(report["pool_available"], true);
    assert_eq!(report["products_determined"], true);
    assert_eq!(report["carryover_sources"], json!([0]));
    assert_eq!(report["candidate_count"], 0);
    let evidence = full(&c, &report["artifact_id"]).await;
    assert_eq!(evidence["origin_workspaces"], args["origin_workspaces"]);
    assert_eq!(
        evidence["report"]["pool"]["source_cleavage"][0]["complete"],
        true
    );
    assert!(
        evidence["report"]["pool"]["source_cleavage"][0]["sites"]
            .as_array()
            .unwrap()
            .is_empty()
    );
    async fn refused(c: &Client, args: Value) {
        let result = c
            .call_tool(
                CallToolRequestParams::new("inspect_golden_gate_reaction")
                    .with_arguments(args.as_object().unwrap().clone()),
            )
            .await;
        assert!(result.is_err() || result.unwrap().is_error == Some(true));
    }
    let mut missing = args.clone();
    missing.as_object_mut().unwrap().remove("origin_workspaces");
    refused(&c, missing).await;
    let mut unattested = args.clone();
    unattested["origin_workspaces"] = json!([]);
    refused(&c, unattested).await;
    let unrelated = Workspace::from_documents(
        "Equal detached circle is not retained history",
        vec![saved.workspace.revisions[&selection.event.revision].clone()],
        0,
    )
    .unwrap();
    let unrelated = call(
        &c,
        "import_workspace_archive",
        json!({"data_base64":STANDARD.encode(workspace::save_workspace(&unrelated).unwrap())}),
    )
    .await;
    let mut foreign = args.clone();
    foreign["origin_workspaces"][0]["workspace"] =
        json!({"id":unrelated["workspace"]["id"],"revision":unrelated["workspace"]["revision"]});
    refused(&c, foreign).await;
    let mut reset = args.clone();
    reset["request"]["rows"][0]["generation"] =
        json!({"dam":"Unmodified","dcm":"Unmodified","cpg":"Unmodified"});
    refused(&c, reset).await;
    args["request"]["rows"][0]["substrate"]["WholeSource"]["origin"] = Value::Null;
    args["origin_workspaces"] = json!([]);
    let detached = call(&c, "inspect_golden_gate_reaction", args).await;
    assert_eq!(detached["pool_available"], false);
    assert_eq!(detached["candidate_count"], 0);
    let detached = full(&c, &detached["artifact_id"]).await;
    assert!(
        detached["report"]["verification"]["products"]
            .as_array()
            .unwrap()
            .is_empty()
    );
    assert_eq!(
        workspace::save_workspace(&saved.workspace).unwrap(),
        original_archive
    );
    c.cancel().await.unwrap();
    server.await.unwrap();
}
