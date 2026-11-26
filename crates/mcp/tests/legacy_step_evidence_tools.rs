//! A legacy TOPO graph-only failure remains downloadable without publishing invented DNA.
use collagene_core::{
    formats, topo as t,
    verification::{self, Operation, Source, Status, TopoSource, VerificationReport},
};
use collagene_mcp::CollageneServer;
use rmcp::{RoleClient, ServiceExt, model::CallToolRequestParams, service::RunningService};
use serde_json::{Value, json};
use t::physical::{self as p, projection};
type Client = RunningService<RoleClient, ()>;
async fn call(c: &Client, n: &str, a: Value) -> Value {
    let r = c
        .call_tool(
            CallToolRequestParams::new(n.to_owned()).with_arguments(a.as_object().unwrap().clone()),
        )
        .await
        .unwrap();
    assert_ne!(r.is_error, Some(true), "{n}: {:?}", r.content);
    r.structured_content.unwrap()
}
async fn artifact(c: &Client, id: &Value) -> Value {
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
            break;
        }
        at = page["next_start"].as_u64().unwrap();
    }
    serde_json::from_str(&text).unwrap()
}
fn reference(imported: &Value, name: &str) -> Value {
    let d = imported["documents"]
        .as_array()
        .unwrap()
        .iter()
        .find(|d| d["name"] == name)
        .unwrap();
    json!({"id":d["id"],"revision":d["revision"]})
}
fn strand(role: &str, length: usize, activated: bool) -> t::PreparedStrand {
    t::PreparedStrand {
        id: role.into(),
        segments: vec![t::SourceInterval {
            source_role: role.into(),
            start: 0,
            length,
            reverse: false,
        }],
        five_prime: if activated {
            p::FivePrime::Phosphate
        } else {
            p::FivePrime::Hydroxyl
        },
        three_prime: if activated {
            p::ThreePrime::Phosphotyrosyl {
                enzyme: p::Enzyme::VacciniaTypeIBTyr274,
                activation_id: role.into(),
            }
        } else {
            p::ThreePrime::Hydroxyl
        },
    }
}
#[tokio::test]
async fn real_old_topo_zero_product_failure_retains_full_bound_step_graph_and_plain_plans_stay_unchanged()
 {
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
    let text = ">vt\nAAGGGGCATCCCTT\n>vb\nGGTGAAGGGATGCCCCTT\n>it\nCACCATGACG\n>ib\nCGTCATGGTG\n>independent_target\nAAGGGGCATCCCTTCACCATGACG\n";
    let local = formats::import_text(text).unwrap().molecules;
    let imported = call(&c, "import_document", json!({"text":text})).await;
    let roles = ["vt", "vb", "it", "ib"];
    let sources = roles
        .iter()
        .map(|id| TopoSource {
            role: (*id).into(),
            source: Source::Input((*id).into()),
        })
        .collect::<Vec<_>>();
    let request=t::Request{recipe:t::GenerationRecipe{schema_version:1,materials:vec![t::MaterialRecipe::Prepared{declaration_id:"Explicit directional source strands".into(),assumption:"Declared source intervals, pairing and activation; no kit or repair is inferred".into(),strands:vec![strand("vt",14,true),strand("vb",18,true),strand("it",10,false),strand("ib",10,false)],pairing:vec![p::PairSpan{strand_a:0,start_a:0,strand_b:1,start_b:4,length:14},p::PairSpan{strand_a:2,start_a:0,strand_b:3,start_b:0,length:10}]}]},transfers:vec![p::Transfer{donor:0,acceptor:2,geometry:p::Geometry::DirectionalCacc},p::Transfer{donor:1,acceptor:3,geometry:p::Geometry::Flush}],engine_limits:p::Limits::default(),projection_options:projection::Options::default()};
    let models = local[..4].iter().collect::<Vec<_>>();
    let native = verification::topo_inputs(&sources, &models, &request).unwrap();
    let selection = t::Selection {
        batch_revision: t::batch_revision(&native, &request).unwrap(),
        component: 0,
        reference_seed: p::Call {
            strand: 0,
            index: 0,
        },
    };
    let expected_batch = selection.batch_revision.clone();
    let expected_product = t::product_revision(&native, &request, &selection).unwrap();
    let inputs = roles
        .iter()
        .map(|id| json!({"id":id,"document":reference(&imported,id)}))
        .collect::<Vec<_>>();
    let target = reference(&imported, "independent_target");
    let proposed=call(&c,"propose_plan",json!({"name":"Legacy directional flap remains graph-only","target":target,"inputs":inputs,"steps":[Operation::Topo{sources,request,selection}]})).await;
    let result = call(
        &c,
        "verify_plan",
        json!({"plan_revision":proposed["plan_revision"]}),
    )
    .await;
    assert_eq!(result["status"], "Incomplete");
    assert_eq!(result["products"], json!([]));
    assert_eq!(result["processed_step_count"], 0);
    assert_eq!(result["reference_step_count"], 0);
    assert!(result["artifact_id"].is_string());
    let report: VerificationReport =
        serde_json::from_value(artifact(&c, &result["artifact_id"]).await).unwrap();
    assert_eq!(
        report.plan_revision,
        proposed["plan_revision"].as_str().unwrap()
    );
    assert_eq!(report.target_revision, target["revision"].as_str().unwrap());
    assert_eq!(report.engine_version, verification::ENGINE_VERSION);
    assert_eq!(report.status, Status::Incomplete);
    assert!(report.products.is_empty());
    assert_eq!(report.step_reports.len(), 1);
    assert_eq!(report.step_reports[0].step_index, 0);
    let graph = report.step_reports[0].topo.as_ref().unwrap();
    assert!(graph.reference.is_none());
    assert_eq!(graph.batch_revision, expected_product);
    assert_eq!(graph.selection.batch_revision, expected_batch);
    assert_eq!(graph.sources.len(), 4);
    assert_eq!(graph.projection.physical.source_bases, 52);
    assert!(matches!(
        graph.projection.outcome,
        projection::Outcome::Withheld(projection::Withheld::UnpairedNucleotides { count: 4 })
    ));
    assert_eq!(
        graph
            .projection
            .physical
            .product
            .unpaired
            .iter()
            .map(|r| r.end - r.start)
            .sum::<usize>(),
        4
    );
    let plain=call(&c,"propose_plan",json!({"name":"Plain identity summary control","target":target,"inputs":[{"id":"vt","document":reference(&imported,"vt")}],"steps":[{"Overlap":{"fragments":[{"Input":"vt"}],"overlaps":[],"circular":false}}]})).await;
    let plain = call(
        &c,
        "verify_plan",
        json!({"plan_revision":plain["plan_revision"]}),
    )
    .await;
    assert!(plain["artifact_id"].is_null());
    assert_eq!(plain["products"].as_array().unwrap().len(), 1);
    assert_eq!(call(&c, "list_documents", json!({})).await["total"], 5);
    drop(c);
    server.abort();
}
