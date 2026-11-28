//! MCP previews physical strands explicitly and attests saved material before reuse.
use base64::{Engine, engine::general_purpose::STANDARD};
use collagene_core::{
    chemistry::{EndRepairMethod, Terminus},
    cloning::{
        assembly::recipe::{self as assembly_recipe, FragmentSpec, Preparation, RunPreparation},
        overlap::assembly::{Options, OuterEndPolicy, Resection},
        restriction::catalog as catalog_digest,
    },
    documents::workspace::{self, Workspace},
    material::{
        lineage::{EventBinding, InputEventBinding},
        outputs as cloning_outputs,
        run::{self as run_material, Selection},
    },
    planning::plan::PlanInput,
    sequence::analysis::Primer,
};
use collagene_mcp::server::CollageneServer;
use collagene_model::{Molecule, Topology};
use rmcp::{RoleClient, ServiceExt, model::CallToolRequestParams, service::RunningService};
use serde_json::{Value, json};
type Client = RunningService<RoleClient, ()>;
const A: &str = "AGTCCGATGCTA";
const B: &str = "GCTATTCGAATTCAGTC";
const CIRCLE: &str = "AGTCCGATGCTATTCGAATTC";
const CUT_FILLED: &str = "AATTCAGTCCGATGCTATTCGAATT";
fn policy() -> Options {
    Options {
        overlaps: vec![4, 4],
        topology: Topology::Circular,
        resection: Resection::ExactHomology,
        outer_end_policy: OuterEndPolicy::JoiningEndsOnly,
        min_paired_bases: 4,
        max_comparisons: 100_000,
        max_candidates: 100,
    }
}
fn recipe() -> assembly_recipe::Recipe {
    assembly_recipe::Recipe {
        name: "Explicit physical circular assembly".into(),
        fragments: ["a", "b"]
            .into_iter()
            .map(|id| FragmentSpec {
                generation: None,
                role: id.into(),
                input_id: id.into(),
                preparation: Preparation::DeclaredDuplex {
                    left: Terminus::PCR,
                    right: Terminus::PCR,
                },
                reverse: false,
            })
            .collect(),
        overlaps: vec![4, 4],
        circular: true,
        options: Default::default(),
    }
}
async fn call(c: &Client, name: &str, args: Value) -> Value {
    let reply = c
        .call_tool(
            CallToolRequestParams::new(name.to_owned())
                .with_arguments(args.as_object().unwrap().clone()),
        )
        .await
        .unwrap();
    assert_ne!(reply.is_error, Some(true), "{name}: {:?}", reply.content);
    reply.structured_content.unwrap()
}
async fn reject(c: &Client, name: &str, args: Value) {
    let reply = c
        .call_tool(
            CallToolRequestParams::new(name.to_owned())
                .with_arguments(args.as_object().unwrap().clone()),
        )
        .await;
    assert!(reply.is_err() || reply.unwrap().is_error == Some(true));
}
async fn artifact(c: &Client, id: &Value) -> Value {
    let (mut start, mut text) = (0, String::new());
    loop {
        let page = call(
            c,
            "read_artifact",
            json!({"artifact_id":id,"start":start,"max_bytes":16384}),
        )
        .await;
        text.push_str(page["text"].as_str().unwrap());
        if page["done"] == true {
            return serde_json::from_str(&text).unwrap();
        }
        start = page["next_start"].as_u64().unwrap();
    }
}
fn named(imported: &Value, name: &str) -> Value {
    let row = imported["documents"]
        .as_array()
        .unwrap()
        .iter()
        .find(|row| row["name"] == name)
        .unwrap();
    json!({"id":row["id"],"revision":row["revision"]})
}
fn binding(w: &Workspace, index: usize) -> EventBinding {
    let d = &w.documents[index];
    EventBinding {
        document_id: d.id.clone(),
        event_id: d.head_event.clone().unwrap(),
        revision: d.head.clone(),
    }
}
fn recorded() -> cloning_outputs::Commit {
    let inputs = [("a", A), ("b", B)]
        .into_iter()
        .map(|(id, dna)| {
            let molecule = Molecule::dna(id, dna, Topology::Linear).unwrap();
            PlanInput {
                id: id.into(),
                revision: molecule.revision(),
                molecule,
            }
        })
        .collect::<Vec<_>>();
    let target = Molecule::dna("Original independent circle", CIRCLE, Topology::Circular).unwrap();
    let compiled =
        assembly_recipe::compile_homology(&inputs, &target, &recipe(), &policy()).unwrap();
    let w = Workspace::from_documents(
        "Saved physical homology source",
        vec![
            inputs[0].molecule.clone(),
            inputs[1].molecule.clone(),
            target,
        ],
        2,
    )
    .unwrap();
    cloning_outputs::finalize(
        &w,
        &compiled.plan,
        binding(&w, 2),
        vec![
            InputEventBinding {
                input_id: "a".into(),
                event: binding(&w, 0),
            },
            InputEventBinding {
                input_id: "b".into(),
                event: binding(&w, 1),
            },
        ],
        &cloning_outputs::Options {
            name: "Saved physical circle".into(),
            show_prepared: false,
            save_primers: false,
            allow_incomplete: false,
        },
    )
    .unwrap()
}
#[tokio::test]
async fn explicit_physical_recipe_and_attested_saved_product_support_actual_pcr_and_digest_reuse() {
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
            .any(|tool| tool.name == "preview_homology_recipe")
    );
    let text=[("a",A,"linear"),("b",B,"linear"),("target",CIRCLE,"circular"),("pcr_target",CIRCLE,"linear"),("cut_target",CUT_FILLED,"linear")].into_iter().map(|(name,dna,topology)|format!("LOCUS       {name:<20} {} bp DNA {topology} SYN 01-JAN-2026\nFEATURES             Location/Qualifiers\nORIGIN\n        1 {}\n//\n",dna.len(),dna.to_ascii_lowercase())).collect::<String>();
    let imported = call(&client, "import_document", json!({"text":text})).await;
    let inputs = json!([{"id":"a","document":named(&imported,"a")},{"id":"b","document":named(&imported,"b")}]);
    let target = named(&imported, "target");
    let args = json!({"inputs":inputs,"target":target,"recipe":recipe(),"policy":policy()});
    let before = call(&client, "list_documents", json!({})).await;
    let preview = call(&client, "preview_homology_recipe", args.clone()).await;
    assert_eq!(preview["state"], "preview");
    assert_eq!(preview["candidate_saved"], false);
    assert_eq!(preview["verified"], false);
    assert_eq!(preview["status"], "ConsistentUnderAssumptions", "{preview}");
    assert_eq!(preview["product"]["length"], 21);
    assert_eq!(preview["physical"]["paired_product_available"], true);
    assert_eq!(preview["physical"]["reusable_material_available"], true);
    assert_eq!(preview["physical"]["incomplete"], false);
    assert_eq!(call(&client, "list_documents", json!({})).await, before);
    // Paired reference strands must not be advertised as reusable material when
    // the declared original termini cannot form sealed junctions.
    let mut diagnostic_recipe = recipe();
    for row in &mut diagnostic_recipe.fragments {
        row.preparation = Preparation::DeclaredDuplex {
            left: Terminus::UNKNOWN,
            right: Terminus::UNKNOWN,
        };
    }
    let diagnostic = call(
        &client,
        "preview_homology_recipe",
        json!({"inputs":inputs,"target":target,"recipe":diagnostic_recipe,"policy":policy()}),
    )
    .await;
    assert_eq!(diagnostic["status"], "Incomplete");
    assert_eq!(diagnostic["physical"]["paired_product_available"], true);
    assert_eq!(diagnostic["physical"]["reusable_material_available"], false);
    assert_eq!(diagnostic["candidate_saved"], false);
    assert_eq!(call(&client, "list_documents", json!({})).await, before);
    let evidence = artifact(&client, &preview["artifact_id"]).await;
    let native: assembly_recipe::Preview =
        serde_json::from_value(evidence["preview"].clone()).unwrap();
    native
        .validate_binding(preview["context_revision"].as_str().unwrap())
        .unwrap();
    assert_eq!(
        native
            .verification
            .final_product()
            .expect("final declared step must reconstruct")
            .sequence
            .as_str(),
        CIRCLE
    );
    let proposed=call(&client,"propose_plan",json!({"name":"Keep explicitly previewed physical homology candidate","inputs":inputs,"target":target,"steps":native.compiled.plan.steps})).await;
    assert_eq!(proposed["verified"], false);
    let verified = call(
        &client,
        "verify_plan",
        json!({"plan_revision":proposed["plan_revision"]}),
    )
    .await;
    assert_eq!(verified["status"], "ConsistentUnderAssumptions");
    let mut changed = args;
    changed["policy"]["overlaps"][0] = json!(3);
    reject(&client, "preview_homology_recipe", changed).await;
    let commit = recorded();
    let selected = Selection {
        run_id: commit.receipt.run.run_id,
        step_index: 2,
        event: commit.receipt.construct,
    };
    let capsule = collagene_core::material::capsule::capsule(&commit.workspace, &selected).unwrap();
    let frame=call(&client,"import_workspace_archive",json!({"data_base64":STANDARD.encode(workspace::save_workspace(&commit.workspace).unwrap())})).await;
    let wr = json!({"id":frame["workspace"]["id"],"revision":frame["workspace"]["revision"]});
    let saved = json!({"id":selected.event.revision,"revision":selected.event.revision});
    let catalog = catalog_digest::builtin_catalog().unwrap();
    let id = catalog
        .edition()
        .enzymes
        .iter()
        .find(|e| e.name == "EcoRI")
        .unwrap()
        .id
        .clone();
    let pin = catalog_digest::select(&catalog, &[id]).unwrap();
    for (preparation, target_name, length) in [
        (
            RunPreparation::ExactPcr {
                forward: Primer::new("F", "AGTCCGAT", "").unwrap(),
                reverse: Primer::new("R", "GAATTCGA", "").unwrap(),
                template_strand: collagene_core::chemistry::Strand::Top,
            },
            "pcr_target",
            21,
        ),
        (
            RunPreparation::CatalogDigest {
                selection: pin,
                conditions: None,
                fragment_index: 0,
                repair: Some(EndRepairMethod::FillFivePrimeOverhangs),
            },
            "cut_target",
            25,
        ),
    ] {
        let recipe = assembly_recipe::Recipe {
            name: "Explicit saved physical material continuation".into(),
            fragments: vec![FragmentSpec {
                generation: None,
                role: "saved".into(),
                input_id: "saved".into(),
                preparation: Preparation::FromRun {
                    origin: capsule.clone(),
                    preparation,
                },
                reverse: false,
            }],
            overlaps: vec![],
            circular: false,
            options: Default::default(),
        };
        let policy = Options {
            overlaps: vec![],
            topology: Topology::Linear,
            ..policy()
        };
        let mut args = json!({"inputs":[{"id":"saved","document":saved}],"target":named(&imported,target_name),"recipe":recipe,"policy":policy});
        reject(&client, "preview_homology_recipe", args.clone()).await;
        args["origin_workspaces"] = json!([{"role":"saved","workspace":wr}]);
        let result = call(&client, "preview_homology_recipe", args).await;
        assert_eq!(result["status"], "ConsistentUnderAssumptions", "{result}");
        assert_eq!(result["product"]["length"], length);
        assert_eq!(result["physical"]["paired_product_available"], true);
    }
    assert_eq!(
        call(
            &client,
            "read_sequence_region",
            json!({"document":named(&imported,"a"),"start":0,"end":12})
        )
        .await["sequence"],
        A
    );
    assert_eq!(
        call(
            &client,
            "read_sequence_region",
            json!({"document":saved,"start":0,"end":21})
        )
        .await["sequence"],
        CIRCLE
    );
    drop(client);
    server.await.unwrap();
}
