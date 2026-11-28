//! In-process SDK acceptance of explicit TA/GC libraries and atomic release.
use base64::{Engine, engine::general_purpose::STANDARD};
use collagene_core::{
    chemistry::{EndTreatment, Strand},
    cloning::{
        assembly::{intent::Selector, recipe::Preparation},
        ligation::{
            tailing::TailingBase,
            tailing_recipe::{Action, Profile},
        },
    },
    documents::workspace::{self, Workspace},
    material::lineage::EventBinding,
    planning::combinatorial::{
        Bin, Choice, ChoiceSelection, Design, Method, SelectedConstruct, TailingAction,
        TargetDeclaration,
    },
    sequence::analysis::Primer,
};
use collagene_mcp::server::CollageneServer;
use collagene_model::{Molecule, Topology};
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
async fn reject(c: &Client, name: &str, args: Value) {
    let result = c
        .call_tool(
            CallToolRequestParams::new(name.to_owned())
                .with_arguments(args.as_object().unwrap().clone()),
        )
        .await;
    assert!(result.is_err() || result.unwrap().is_error == Some(true));
}
async fn setup() -> (Client, tokio::task::JoinHandle<()>) {
    let (a, b) = tokio::io::duplex(65_536);
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
async fn text_artifact(c: &Client, id: &Value) -> Value {
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
            break;
        }
        start = page["next_start"].as_u64().unwrap();
    }
    serde_json::from_str(&text).unwrap()
}
async fn native_artifact(c: &Client, id: &Value) -> Workspace {
    let (mut start, mut text) = (0, String::new());
    loop {
        let page = call(
            c,
            "read_artifact",
            json!({"artifact_id":id,"start":start,"max_bytes":16384}),
        )
        .await;
        assert_eq!(page["encoding"], "base64");
        text.push_str(page["text"].as_str().unwrap());
        if page["done"] == true {
            break;
        }
        start = page["next_start"].as_u64().unwrap();
    }
    workspace::load_workspace(&STANDARD.decode(text).unwrap()).unwrap()
}
fn binding(w: &Workspace, index: usize) -> EventBinding {
    let doc = &w.documents[index];
    EventBinding {
        document_id: doc.id.clone(),
        event_id: doc.head_event.clone().unwrap(),
        revision: doc.head.clone(),
    }
}
fn workspace_ref(value: &Value) -> Value {
    json!({"id":value["workspace"]["id"],"revision":value["workspace"]["revision"]})
}
async fn upload(c: &Client, w: &Workspace) -> Value {
    call(
        c,
        "import_workspace_archive",
        json!({"data_base64":STANDARD.encode(workspace::save_workspace(w).unwrap())}),
    )
    .await
}

fn fixture(profile: Profile) -> (Workspace, Design, SelectedConstruct) {
    let expected = match profile {
        Profile::Ta => "AAACCCTGGGTTTA",
        Profile::Gc => "AAACCCCGGGTTTG",
    };
    let w = Workspace::from_documents(
        "SDK native tailing library",
        vec![
            Molecule::dna("Vector stock", "AAACCC", Topology::Linear).unwrap(),
            Molecule::dna("Insert stock", "GGGTTT", Topology::Linear).unwrap(),
            Molecule::dna("Independent scarred target", expected, Topology::Circular).unwrap(),
            Molecule::dna("Independent wrong target", "AACCGG", Topology::Circular).unwrap(),
        ],
        0,
    )
    .unwrap();
    let (left, right) = match profile {
        Profile::Ta => (TailingBase::T, TailingBase::A),
        Profile::Gc => (TailingBase::C, TailingBase::G),
    };
    let d = Design {
        id: "sdk-tail-library".into(),
        name: "Literal complementary tailed circles".into(),
        circular: true,
        material_bindings: vec![],
        method: Method::Tailing {
            profile,
            actions: vec![
                TailingAction {
                    bin_id: "vector".into(),
                    action: Action::AddTail { base: left },
                },
                TailingAction {
                    bin_id: "insert".into(),
                    action: Action::AddTail { base: right },
                },
            ],
            allow_nicked: false,
        },
        bins: (0..2)
            .map(|i| Bin {
                id: if i == 0 { "vector" } else { "insert" }.into(),
                label: format!("Contribution {i}"),
                choices: vec![Choice {
                    anneal_bottom: None,
                    id: format!("choice-{i}"),
                    label: format!("Source {i}"),
                    source: binding(&w, i),
                    selector: Selector::ReferenceRange {
                        start: 0,
                        end: 6,
                        reverse: false,
                    },
                    preparation_reverse: false,
                    preparation: Some(Preparation::ExactPcr {
                        forward: Primer::new("F", if i == 0 { "AAA" } else { "GGG" }, "").unwrap(),
                        reverse: Primer::new("R", if i == 0 { "GGG" } else { "AAA" }, "").unwrap(),
                        template_strand: Strand::Top,
                    }),
                    treatment: Some(EndTreatment::Phosphorylate),
                    generation: None,
                }],
            })
            .collect(),
    };
    let c = SelectedConstruct {
        id: "correct".into(),
        name: "Independent tailed circle".into(),
        choices: vec![
            ChoiceSelection {
                bin_id: "vector".into(),
                choice_id: "choice-0".into(),
            },
            ChoiceSelection {
                bin_id: "insert".into(),
                choice_id: "choice-1".into(),
            },
        ],
        target: TargetDeclaration::Existing(binding(&w, 2)),
        overrides: None,
    };
    (w, d, c)
}
fn output_options() -> Value {
    json!({"name":"Released tailed circle","show_prepared":false,"save_primers":true,"allow_incomplete":false})
}
#[tokio::test]
async fn ta_gc_sdk_keeps_profile_and_target_verdicts_distinct_and_release_atomic() {
    let (c, server) = setup().await;
    for profile in [Profile::Ta, Profile::Gc] {
        let (w, d, row) = fixture(profile);
        let source_bytes = workspace::save_workspace(&w).unwrap();
        let imported = upload(&c, &w).await;
        let source = workspace_ref(&imported);
        let inspected = call(
            &c,
            "inspect_combinatorial_construct",
            json!({"workspace":source,"design":d,"construct":row}),
        )
        .await;
        let inspection = text_artifact(&c, &inspected["artifact_id"]).await;
        assert_eq!(
            inspection["recipe"]["Tailing"]["recipe"]["profile"],
            json!(profile)
        );
        assert_eq!(
            inspection["recipe"]["Tailing"]["recipe"]["fragments"]
                .as_array()
                .unwrap()
                .len(),
            2
        );
        assert_eq!(inspection["target_binding"], json!(binding(&w, 2)));
        assert_eq!(inspection["verified"], false);
        let mut missing = serde_json::to_value(&d).unwrap();
        missing["method"]["Tailing"]["actions"]
            .as_array_mut()
            .unwrap()
            .pop();
        reject(
            &c,
            "inspect_combinatorial_construct",
            json!({"workspace":source,"design":missing,"construct":row}),
        )
        .await;
        let mut wrong = row.clone();
        wrong.id = "wrong".into();
        wrong.target = TargetDeclaration::Existing(binding(&w, 3));
        let saved=call(&c,"evaluate_combinatorial_batch",json!({"workspace":source,"matrix":{"design":d,"rows":[{"Ready":row},{"Ready":wrong}]}})).await;
        let record = text_artifact(&c, &saved["artifact_id"]).await;
        assert_eq!(record["schema_version"], 3);
        assert_eq!(record["rows"][0]["status"], "ConsistentUnderAssumptions");
        assert_eq!(record["rows"][1]["status"], "Invalid");
        for i in 0..2 {
            assert_eq!(
                record["rows"][i]["tailing_profile"]["profile"],
                json!(profile)
            );
            assert_eq!(record["rows"][i]["tailing_profile"]["status"], "Matched");
            assert!(
                !record["rows"][i]["tailing_profile"]["findings"]
                    .as_array()
                    .unwrap()
                    .is_empty()
            );
        }
        let saved_ref = workspace_ref(&saved);
        let batch_hash = saved["receipt"]["attachment_hash"].clone();
        reject(
            &c,
            "release_combinatorial_outputs",
            json!({"workspace":saved_ref,"request":{"batch_hash":batch_hash,"rows":[
                {"construct_id":"correct","options":output_options()},
                {"construct_id":"wrong","options":output_options()},
            ]}}),
        )
        .await;
        let after_failed = call(
            &c,
            "list_combinatorial_releases",
            json!({"workspace":saved_ref,"limit":1}),
        )
        .await;
        assert_eq!(after_failed["total"], 0);
        let released = call(
            &c,
            "release_combinatorial_outputs",
            json!({"workspace":saved_ref,"request":{"batch_hash":batch_hash,"rows":[
                {"construct_id":"correct","options":output_options()},
            ]}}),
        )
        .await;
        let native = native_artifact(&c, &released["native_artifact_id"]).await;
        let manifest = text_artifact(&c, &released["artifact_id"]).await;
        let output = &manifest["rows"][0]["output"];
        assert!(output["primer_attachment"].is_string());
        let product = &native.revisions[output["construct"]["revision"].as_str().unwrap()];
        let expected = match profile {
            Profile::Ta => "AAACCCTGGGTTTA",
            Profile::Gc => "AAACCCCGGGTTTG",
        };
        assert_eq!(product.sequence.len(), expected.len());
        assert!(format!("{expected}{expected}").contains(product.sequence.as_str()));
        assert_eq!(product.topology, Topology::Circular);
        for doc in &w.documents {
            assert_eq!(native.head(&doc.id).unwrap(), w.head(&doc.id).unwrap());
        }
        let reopened = upload(&c, &native).await;
        let opened=call(&c,"open_combinatorial_release",json!({"workspace":workspace_ref(&reopened),"manifest_hash":released["receipt"]["manifest_hash"]})).await;
        assert_eq!(text_artifact(&c, &opened["artifact_id"]).await, manifest);
        assert_eq!(source_bytes, workspace::save_workspace(&w).unwrap());
    }
    c.cancel().await.unwrap();
    server.await.unwrap();
}
