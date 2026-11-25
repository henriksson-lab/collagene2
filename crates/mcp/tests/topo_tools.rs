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
    let sealing = json!({"ligase_assumption":"Explicit ideal ATP-dependent ligase nick sealing under compatible conditions; no host repair", "nicks":[
        {"three":{"strand":2,"index":6},"five":{"strand":0,"index":0}},
        {"three":{"strand":3,"index":6},"five":{"strand":1,"index":0}}
    ]});
    let processed = call(
        &c,
        "seal_topo_nicks",
        json!({"sources":sources,"request":{"preparation":request,"sealing":sealing}}),
    )
    .await;
    assert_eq!(processed["verified"], false);
    assert_eq!(processed["physical"]["nicks"], 0);
    assert_eq!(processed["physical"]["components"], 1);
    let proof = evidence(&c, &processed["artifact_id"]).await;
    assert_eq!(
        proof["report"]["sealing"]["parent"]["product"]["nicks"]
            .as_array()
            .unwrap()
            .len(),
        2
    );
    assert!(
        proof["report"]["sealing"]["product"]["chains"]
            .as_array()
            .unwrap()
            .iter()
            .all(|chain| chain["ends"] == "CovalentlyCircular")
    );
    let projected=call(&c,"project_processed_topo",json!({"sources":sources,"request":{
        "treatment":{"preparation":request,"sealing":sealing},
        "selection":{"sealing_context_revision":proof["report"]["sealing"]["context_revision"],"component":0,"reference_seed":{"strand":0,"index":0}},
        "options":{"require_all_requested_transfers_applied":true,"require_known_terminal_chemistry":true,"limits":t::physical::projection::ProjectionLimits::default()}
    }})).await;
    assert_eq!(projected["verified"], false);
    assert_eq!(projected["contour_proven"], true);
    assert!(!projected["document"].is_null());
    assert_eq!(projected["duplex"]["covalently_closed_in_model"], true);
    let product_evidence = evidence(&c, &projected["artifact_id"]).await;
    let molecule: Molecule =
        serde_json::from_value(product_evidence["report"]["product"].clone()).unwrap();
    assert_eq!(molecule.sequence.as_str(), "AGGGGCATCCCTTCGAATGA");
    assert!(collagene_core::chemistry::has_unmapped_physical_evidence(
        &molecule
    ));
    assert_eq!(
        product_evidence["report"]["projection"]["sealing"]["parent"]["product"]["nicks"]
            .as_array()
            .unwrap()
            .len(),
        2
    );
    let independent_target = import(
        &c,
        "Independent intended processed clone",
        "AGGGGCATCCCTTCGAATGA",
        true,
    )
    .await;
    let plan_inputs = sources
        .as_array()
        .unwrap()
        .iter()
        .map(|source| json!({"id":source["role"],"document":source["document"]}))
        .collect::<Vec<_>>();
    let role_sources = sources
        .as_array()
        .unwrap()
        .iter()
        .map(|source| json!({"role":source["role"],"source":{"Input":source["role"]}}))
        .collect::<Vec<_>>();
    let proposed=call(&c,"propose_plan",json!({"name":"Source-bound processed TA clone","target":independent_target,"inputs":plan_inputs,"steps":[{"ProcessedTopo":{"sources":role_sources,"request":product_evidence["report"]["request"]}}]})).await;
    assert_eq!(proposed["verified"], false);
    let verified = call(
        &c,
        "verify_plan",
        json!({"plan_revision":proposed["plan_revision"]}),
    )
    .await;
    assert_eq!(verified["status"], "ConsistentUnderAssumptions");
    let full_plan = evidence(&c, &verified["artifact_id"]).await;
    assert_eq!(
        full_plan["processed_step_reports"]
            .as_array()
            .unwrap()
            .len(),
        1
    );
    let mut wrong = sealing.clone();
    wrong["nicks"][0]["five"]["index"] = json!(1);
    let retained = call(
        &c,
        "seal_topo_nicks",
        json!({"sources":sources,"request":{"preparation":request,"sealing":wrong}}),
    )
    .await;
    assert!(
        retained["physical"]["nicks"].as_u64().unwrap() > 0,
        "wrong actual terminus cannot be sealed"
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

#[tokio::test]
async fn directional_preparation_compiles_actual_sources_without_inventing_repaired_contour() {
    use collagene_core::{analysis::Primer, topo_preparation as prep};
    use t::physical::{self as p, projection};
    let (c, s) = setup().await;
    let mut sources = vec![];
    for (role, dna) in [
        ("vt", "AAGGGGCATCCCTT"),
        ("vb", "GGTGAAGGGATGCCCCTT"),
        ("template", "CACCATGACG"),
    ] {
        sources.push(json!({"role":role,"document":import(&c,role,dna,false).await}));
    }
    let declaration = prep::PreparationRequest {
        mode: prep::PreparationMode::DirectionalCacc,
        vector: prep::PreparedVectorDeclaration {
            declaration_id: "literal directional fixture".into(),
            assumption: "Explicit synthetic activation hypothesis, not a stock observation".into(),
            strands: vec![strand("vt", 14, true), strand("vb", 18, true)],
            pairing: vec![p::PairSpan {
                strand_a: 0,
                start_a: 0,
                strand_b: 1,
                start_b: 4,
                length: 14,
            }],
        },
        insert: prep::ExactPcrDeclaration {
            material_id: "actual blunt PCR".into(),
            template_role: "template".into(),
            forward: Primer::new("Forward", "CACCA", "").unwrap(),
            reverse: Primer::new("Reverse", "CGTCA", "").unwrap(),
            tail: None,
        },
        engine_limits: p::Limits::default(),
        projection_options: projection::Options::default(),
    };
    let compiled = call(
        &c,
        "compile_topo_preparation",
        json!({"sources":sources,"declaration":declaration}),
    )
    .await;
    assert_eq!(compiled["verified"], false);
    assert_eq!(compiled["compiled"]["directional"], true);
    assert_eq!(compiled["compiled"]["request"]["transfers"], json!([]));
    let inspection=call(&c,"inspect_topo_orientations",json!({"sources":sources,"request":compiled["compiled"]["request"],"vector":compiled["compiled"]["vector_strands"],"insert":compiled["compiled"]["insert_strands"],"directional":true})).await;
    let full = evidence(&c, &inspection["artifact_id"]).await;
    assert!(
        full.to_string().contains("GGTGAAGGGATGCCCCTTCGTCATGGTG"),
        "physical graph must retain the actual displaced bases"
    );
    assert_eq!(inspection["biochemical_mixture_exhaustive"], false);
    let mut preparation = compiled["compiled"]["request"].clone();
    preparation["transfers"] = inspection["candidates"][0]["transfers"].clone();
    let processing = json!({"actions":[
        {"CleaveThreePrimeFlap":{"bond":{"upstream":{"strand":3,"index":5},"downstream":{"strand":3,"index":6}},"displacement_step":0,"profile":"IdealDuplexAnchoredThreePrimeHydrolysisV1","assumption":"Explicit ideal hydrolytic fixture only; no named nuclease or host guarantee"}},
        {"SealNick":{"selection":{"three":{"strand":2,"index":9},"five":{"strand":0,"index":0}},"ligase_assumption":"Explicit ideal templated ligase fixture"}},
        {"SealNick":{"selection":{"three":{"strand":3,"index":5},"five":{"strand":1,"index":0}},"ligase_assumption":"Explicit ideal templated ligase fixture"}}
    ],"limits":t::physical::processing::ProcessingLimits::default()});
    let inspected = call(
        &c,
        "inspect_topo_actions",
        json!({"sources":sources,"request":{"preparation":preparation,"processing":processing}}),
    )
    .await;
    assert_eq!(inspected["verified"], false);
    assert_eq!(inspected["source_calls"], 52);
    assert_eq!(inspected["component_count"], 2);
    let main = inspected["components"]
        .as_array()
        .unwrap()
        .iter()
        .find(|v| v["paired_calls"] == 48 && v["unpaired_calls"] == 0)
        .unwrap()["component"]
        .clone();
    let request = json!({"preparation":preparation,"processing":processing,
        "selection":{"processing_context_revision":inspected["processing_context_revision"],"component":main,"reference_seed":{"strand":0,"index":0}},
        "options":{"require_all_requested_transfers_applied":true,"require_known_terminal_chemistry":true,"limits":t::physical::projection::ProjectionLimits::default()}});
    let projected = call(
        &c,
        "project_topo_actions",
        json!({"sources":sources,"request":request}),
    )
    .await;
    assert_eq!(projected["contour_proven"], true);
    assert_eq!(projected["component_count"], 2);
    let data = evidence(&c, &projected["artifact_id"]).await;
    let product: Molecule = serde_json::from_value(data["report"]["product"].clone()).unwrap();
    assert_eq!(product.sequence.as_str(), "AAGGGGCATCCCTTCACCATGACG");
    assert!(collagene_core::chemistry::has_unmapped_physical_evidence(
        &product
    ));
    let independent_target = import(
        &c,
        "Independent intended directional clone",
        "AAGGGGCATCCCTTCACCATGACG",
        true,
    )
    .await;
    let bindings = sources
        .iter()
        .map(|s| json!({"id":s["role"],"document":s["document"]}))
        .collect::<Vec<_>>();
    let roles = sources
        .iter()
        .map(|s| json!({"role":s["role"],"source":{"Input":s["role"]}}))
        .collect::<Vec<_>>();
    let proposed=call(&c,"propose_plan",json!({"name":"Explicit accounted directional action fixture","target":independent_target,"inputs":bindings,"steps":[{"TopoActions":{"sources":roles,"request":request}}]})).await;
    let verified = call(
        &c,
        "verify_plan",
        json!({"plan_revision":proposed["plan_revision"]}),
    )
    .await;
    assert_eq!(verified["status"], "ConsistentUnderAssumptions");
    let full = evidence(&c, &verified["artifact_id"]).await;
    assert_eq!(
        full["physical_processing_reports"]
            .as_array()
            .unwrap()
            .len(),
        1
    );
    assert_eq!(
        full["physical_processing_reports"][0]["evidence"]["TopoActions"]["projection"]["processing"]
            ["source_bases"],
        52
    );
    let mut upstream_roles = roles.clone();
    upstream_roles
        .iter_mut()
        .find(|r| r["role"] == "template")
        .unwrap()["source"] = json!({"Step":0});
    let proposed = call(&c, "propose_plan", json!({
        "name":"Native modified template copied by PCR into TOPO preparation", "target":independent_target,
        "inputs":bindings, "steps":[
            {"DeclareSourceModifications":{"template":{"Input":"template"},"assumptions":{"dam":"Unknown","dcm":"Unknown","cpg":"FullyModified"}}},
            {"TopoActions":{"sources":upstream_roles,"request":request}}
        ]
    })).await;
    let verified = call(
        &c,
        "verify_plan",
        json!({"plan_revision":proposed["plan_revision"]}),
    )
    .await;
    assert_eq!(verified["status"], "ConsistentUnderAssumptions");
    let full = evidence(&c, &verified["artifact_id"]).await;
    assert_eq!(full["products"].as_array().unwrap().len(), 2);
    assert_eq!(full["physical_processing_reports"][0]["step_index"], 1);
    let before: collagene_core::modified_duplex::ModificationSummary = serde_json::from_str(
        full["products"][0]["metadata"]["tracked_modifications"]
            .as_str()
            .unwrap(),
    )
    .unwrap();
    assert!(before.known_modified > 0);
    let actual: collagene_core::modified_duplex::TrackedDuplex = serde_json::from_value(
        full["physical_processing_reports"][0]["evidence"]["TopoActions"]["tracked"].clone(),
    )
    .unwrap();
    assert_eq!(actual.top.as_str(), "AAGGGGCATCCCTTCACCATGACG");
    assert_eq!(
        collagene_core::modified_duplex::summary(&actual)
            .unwrap()
            .known_modified,
        0
    );
    assert!(actual.generation.is_empty());
    let mut stale = sources.clone();
    stale[0]["document"]["revision"] = json!("0".repeat(64));
    reject(
        &c,
        "compile_topo_preparation",
        json!({"sources":stale,"declaration":declaration}),
    )
    .await;
    c.cancel().await.unwrap();
    s.abort();
}
