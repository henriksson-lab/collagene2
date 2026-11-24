//! Actual SDK TOPO declarations and independent target checks. The expected TA
//! contour is a literal, not copied from an imported alleged physical product.
use collagene_core::{Molecule, Topology, formats, topo as t};
use collagene_mcp::CollageneServer;
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
async fn import(c: &Client, name: &str, dna: &str, circular: bool) -> Value {
    let m = Molecule::dna(
        name,
        dna,
        if circular {
            Topology::Circular
        } else {
            Topology::Linear
        },
    )
    .unwrap();
    let r = call(
        c,
        "import_document",
        json!({"text":formats::export_genbank(&m).unwrap()}),
    )
    .await;
    json!({"id":r["documents"][0]["id"],"revision":r["documents"][0]["revision"]})
}
async fn evidence(c: &Client, id: &Value) -> Value {
    let mut text = String::new();
    let mut start = 0;
    loop {
        let page = call(
            c,
            "read_artifact",
            json!({"artifact_id":id,"start":start,"max_bytes":4096}),
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
fn strand(id: &str, n: usize, active: bool) -> t::PreparedStrand {
    use t::physical::{Enzyme, FivePrime, ThreePrime};
    t::PreparedStrand {
        id: id.into(),
        segments: vec![t::SourceInterval {
            source_role: id.into(),
            start: 0,
            length: n,
            reverse: false,
        }],
        five_prime: if active {
            FivePrime::Phosphate
        } else {
            FivePrime::Hydroxyl
        },
        three_prime: if active {
            ThreePrime::Phosphotyrosyl {
                enzyme: Enzyme::VacciniaTypeIBTyr274,
                activation_id: id.into(),
            }
        } else {
            ThreePrime::Hydroxyl
        },
    }
}
async fn fixture(c: &Client) -> (Value, t::Request) {
    use t::physical::{self as p, projection};
    let mut sources = vec![];
    for (role, dna) in [
        ("vt", "AGGGGCATCCCTT"),
        ("vb", "AGGGATGCCCCTT"),
        ("it", "CGAATGA"),
        ("ib", "CATTCGA"),
    ] {
        sources.push(json!({"role":role,"document":import(c,role,dna,false).await}));
    }
    let request=t::Request{recipe:t::GenerationRecipe{schema_version:1,materials:vec![t::MaterialRecipe::Prepared{
        declaration_id:"literal TA declaration".into(),assumption:"Prepared strands/activation are an explicit hypothesis, not inferred from record metadata".into(),
        strands:vec![strand("vt",13,true),strand("vb",13,true),strand("it",7,false),strand("ib",7,false)],
        pairing:vec![p::PairSpan{strand_a:0,start_a:0,strand_b:1,start_b:0,length:12},p::PairSpan{strand_a:2,start_a:0,strand_b:3,start_b:0,length:6}]}]},
        transfers:vec![p::Transfer{donor:0,acceptor:2,geometry:p::Geometry::ThreePrimeTA},p::Transfer{donor:1,acceptor:3,geometry:p::Geometry::ThreePrimeTA}],
        engine_limits:p::Limits::default(),projection_options:projection::Options::default()};
    (json!(sources), request)
}

#[tokio::test]
async fn topo_sdk_retains_nicks_all_orientations_and_explicit_independent_target_policy() {
    let (c, s) = setup().await;
    let tools = c.list_tools(Default::default()).await.unwrap().tools;
    assert_eq!(tools.len(), collagene_mcp::CURRENT_TOOL_COUNT);
    for name in [
        "simulate_topo",
        "inspect_topo_orientations",
        "project_topo_product",
        "check_topo_target",
    ] {
        assert_eq!(
            tools.iter().find(|t| t.name == name).unwrap().input_schema["additionalProperties"],
            false
        );
    }
    let (sources, request) = fixture(&c).await;
    let simulated = call(
        &c,
        "simulate_topo",
        json!({"sources":sources,"request":request}),
    )
    .await;
    assert_eq!(
        simulated["physical"]["all_requested_transfers_applied"],
        true
    );
    assert_eq!(simulated["physical"]["component_count"], 1);
    assert_eq!(simulated["physical"]["nick_count"], 2);
    assert!(
        simulated["physical"].get("product").is_none(),
        "large physical graph belongs in artifact"
    );
    let selection = json!({"batch_revision":simulated["batch_revision"],"component":0,"reference_seed":{"strand":0,"index":0}});
    let projected = call(
        &c,
        "project_topo_product",
        json!({"sources":sources,"request":request,"selection":selection}),
    )
    .await;
    assert_eq!(projected["projection"]["state"], "paired_contour");
    assert_eq!(projected["projection"]["nick_count"], 2);
    let artifact = evidence(&c, &projected["artifact_id"]).await;
    assert_eq!(artifact["source_bindings"], sources);
    assert_eq!(
        artifact["report"]["reference"]["sequence"]["bases"],
        "AGGGGCATCCCTTCGAATGA"
    );
    assert_eq!(
        artifact["report"]["reference"]["metadata"]["source:legacy:chemistry_unmapped"],
        "true"
    );
    let target = import(&c, "independent target", "AGGGGCATCCCTTCGAATGA", true).await;
    let mut policy = t::AcceptancePolicy::default();
    let args = json!({"sources":sources,"request":request,"selection":selection,"target":target,"policy":policy});
    let strict = call(&c, "check_topo_target", args.clone()).await;
    assert!(
        strict["outcome"].get("Rejected").is_some(),
        "DNA equality cannot seal two nicks"
    );
    policy.allow_nicks = true;
    let mut allowed = args.clone();
    allowed["policy"] = serde_json::to_value(&policy).unwrap();
    let accepted = call(&c, "check_topo_target", allowed.clone()).await;
    assert_eq!(accepted["outcome"], "AcceptedUnderDeclaredGraph");
    assert_ne!(accepted["batch_revision"], strict["batch_revision"]);
    let wrong = import(&c, "wrong independent target", "TTTTTTTTTTTTTTTTTTTT", true).await;
    allowed["target"] = wrong;
    let rejected = call(&c, "check_topo_target", allowed).await;
    assert!(rejected["outcome"].get("Rejected").is_some());
    let wrong_evidence = evidence(&c, &rejected["artifact_id"]).await;
    assert_eq!(
        wrong_evidence["report"]["product"]["reference"]["sequence"]["bases"],
        "AGGGGCATCCCTTCGAATGA",
        "target cannot supply generated DNA"
    );
    let mut orientations = request.clone();
    orientations.transfers.clear();
    let orientations=call(&c,"inspect_topo_orientations",json!({"sources":sources,"request":orientations,"vector":[0,1],"insert":[2,3],"directional":false})).await;
    assert_eq!(orientations["candidate_count"], 2);
    assert_eq!(orientations["orientation_search_complete"], true);
    assert_eq!(orientations["biochemical_mixture_exhaustive"], false);
    let full = evidence(&c, &orientations["artifact_id"]).await;
    assert_eq!(
        full["report"]["candidate_requests"]
            .as_array()
            .unwrap()
            .len(),
        2
    );
    assert_eq!(
        full["report"]["reference"]["candidates"]
            .as_array()
            .unwrap()
            .len(),
        2
    );
    c.cancel().await.unwrap();
    s.await.unwrap();
}

#[tokio::test]
async fn topo_sdk_rejects_foreign_claims_stale_sources_and_selection_without_repairing_state() {
    let (c, s) = setup().await;
    let (sources, request) = fixture(&c).await;
    let simulated = call(
        &c,
        "simulate_topo",
        json!({"sources":sources,"request":request}),
    )
    .await;
    let mut stale = sources.clone();
    stale[0]["document"]["revision"] = json!("0".repeat(64));
    reject(
        &c,
        "simulate_topo",
        json!({"sources":stale,"request":request}),
    )
    .await;
    let mut duplicate = sources.clone();
    duplicate[1]["role"] = duplicate[0]["role"].clone();
    reject(
        &c,
        "simulate_topo",
        json!({"sources":duplicate,"request":request}),
    )
    .await;
    reject(
        &c,
        "simulate_topo",
        json!({"sources":sources,"request":request,"foreign_physical_product":simulated}),
    )
    .await;
    let stale_selection = json!({"batch_revision":"0".repeat(64),"component":0,"reference_seed":{"strand":0,"index":0}});
    reject(
        &c,
        "project_topo_product",
        json!({"sources":sources,"request":request,"selection":stale_selection}),
    )
    .await;
    let mut unactivated = request.clone();
    if let t::MaterialRecipe::Prepared { strands, .. } = &mut unactivated.recipe.materials[0] {
        strands[0].three_prime = t::physical::ThreePrime::Hydroxyl;
    }
    let unactivated = call(
        &c,
        "simulate_topo",
        json!({"sources":sources,"request":unactivated}),
    )
    .await;
    assert_eq!(
        unactivated["physical"]["all_requested_transfers_applied"],
        false
    );
    reject(&c,"inspect_topo_orientations",json!({"sources":sources,"request":request,"vector":[0,1],"insert":[2,3],"directional":false})).await;
    c.cancel().await.unwrap();
    s.await.unwrap();
}
