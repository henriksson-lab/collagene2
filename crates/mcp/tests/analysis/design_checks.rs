//! Real in-memory SDK routing, immutable references and native evidence binding.
use collagene_core::{
    cloning::{assembly::standards, pcr::collection_binding},
    construct_qc,
};
use collagene_mcp::server::CollageneServer;
use collagene_model::{Feature, Location, Molecule, Sequence, Topology};
use rmcp::{RoleClient, ServiceExt, model::CallToolRequestParams, service::RunningService};
use serde_json::{Value, json};
type Client = RunningService<RoleClient, ()>;

async fn setup(molecule: &Molecule) -> (Client, tokio::task::JoinHandle<()>, Value) {
    let server = CollageneServer::default();
    let imported = server
        .preload_document_bytes(&collagene_core::protocol::save_document(molecule).unwrap())
        .unwrap()
        .structured_content
        .unwrap();
    let document = json!({"id":imported["documents"][0]["id"],"revision":imported["documents"][0]["revision"]});
    let (a, b) = tokio::io::duplex(65536);
    let service = tokio::spawn(async move {
        server.serve(a).await.unwrap().waiting().await.unwrap();
    });
    (().serve(b).await.unwrap(), service, document)
}
async fn call(client: &Client, name: &str, arguments: Value) -> Value {
    let result = client
        .call_tool(
            CallToolRequestParams::new(name.to_owned())
                .with_arguments(arguments.as_object().unwrap().clone()),
        )
        .await
        .unwrap();
    assert_ne!(result.is_error, Some(true), "{name}: {:?}", result.content);
    result.structured_content.unwrap()
}
async fn reject(client: &Client, name: &str, arguments: Value) {
    let result = client
        .call_tool(
            CallToolRequestParams::new(name.to_owned())
                .with_arguments(arguments.as_object().unwrap().clone()),
        )
        .await;
    assert!(result.is_err() || result.unwrap().is_error == Some(true));
}
async fn artifact(client: &Client, id: &Value) -> Value {
    let (mut offset, mut text) = (0, String::new());
    loop {
        let p = call(
            client,
            "read_artifact",
            json!({"artifact_id":id,"start":offset,"max_bytes":16384}),
        )
        .await;
        text.push_str(p["text"].as_str().unwrap());
        if p["done"] == true {
            break;
        }
        offset = p["next_start"].as_u64().unwrap();
    }
    serde_json::from_str(&text).unwrap()
}

#[tokio::test]
async fn construct_qc_routes_declared_roles_without_claiming_annotation_truth_or_global_success() {
    let mut molecule = Molecule::dna("explicit CDS", "ATGAAATAA", Topology::Linear).unwrap();
    molecule.features.push(Feature {
        kind: "CDS".into(),
        location: Location::range(0, 9),
        qualifiers: vec![],
        needs_review: false,
    });
    let mut options = construct_qc::Options::default();
    options
        .annotations
        .push(construct_qc::AnnotationAssumption {
            feature_index: 0,
            source_revision: molecule.revision(),
            annotation_revision: construct_qc::annotation_revision(&molecule.features[0]).unwrap(),
            rationale:
                "Literal source-role declaration for structural checking, not verified biology"
                    .into(),
        });
    let (client, service, document) = setup(&molecule).await;
    let before = call(&client, "list_documents", json!({})).await;
    let tools = client.list_tools(Default::default()).await.unwrap().tools;
    assert_eq!(tools.len(), collagene_mcp::CURRENT_TOOL_COUNT);
    for name in [
        "list_assembly_standards",
        "check_assembly_standard",
        "check_construct_qc",
        "search_oligo_collection",
    ] {
        assert!(tools.iter().any(|tool| tool.name.as_ref() == name));
    }
    let arguments = json!({"document":document,"options":options,"max_findings":1});
    let first = call(&client, "check_construct_qc", arguments.clone()).await;
    assert_eq!(first["source_revision"], molecule.revision());
    assert_eq!(first["findings"].as_array().unwrap().len(), 1);
    assert_eq!(first["done"], false);
    assert!(first["counts"]["pass"].as_u64().unwrap() > 0);
    let full = artifact(&client, &first["artifact_id"]).await;
    assert_eq!(full["statement"], construct_qc::STATEMENT);
    assert!(
        full["findings"]
            .as_array()
            .unwrap()
            .iter()
            .any(|f| f["rule"] == "CdsFrame"
                && f["verdict"] == "Pass"
                && f["under_declared_annotation_assumptions"] == true)
    );
    let report: construct_qc::Report = serde_json::from_value(full).unwrap();
    report.validate_for(&molecule, &options).unwrap();
    report
        .validate_binding(first["context_revision"].as_str().unwrap())
        .unwrap();
    let second = call(
        &client,
        "check_construct_qc",
        json!({"document":document,"options":options,"finding_start":1,"max_findings":1}),
    )
    .await;
    assert_eq!(second["artifact_id"], first["artifact_id"]);
    let mut stale_assumption = arguments.clone();
    stale_assumption["options"]["annotations"][0]["annotation_revision"] = json!("0".repeat(64));
    let uncertain = call(&client, "check_construct_qc", stale_assumption).await;
    assert!(uncertain["counts"]["indeterminate"].as_u64().unwrap() > 0);
    assert_ne!(uncertain["context_revision"], first["context_revision"]);
    let mut stale_source = arguments.clone();
    stale_source["document"]["revision"] = json!("0".repeat(64));
    reject(&client, "check_construct_qc", stale_source).await;
    let mut invalid = arguments;
    invalid["max_findings"] = json!(101);
    reject(&client, "check_construct_qc", invalid).await;
    assert_eq!(call(&client, "list_documents", json!({})).await, before);
    client.cancel().await.unwrap();
    service.await.unwrap();
}

#[tokio::test]
async fn standard_discovery_and_fragment_checks_preserve_profile_scope_and_original_coordinates() {
    let formatted = "GAATTCGCGGCCGCTTCTAGAGAAGCGGCCGCTTTACTAGTAGCGGCCGCTGCAG";
    let molecule = Molecule::dna(
        "vector context",
        &format!("TTTT{formatted}GGGG"),
        Topology::Linear,
    )
    .unwrap();
    let (client, service, document) = setup(&molecule).await;
    let before = call(&client, "list_documents", json!({})).await;
    let catalogue = call(&client, "list_assembly_standards", json!({})).await;
    assert_eq!(catalogue["profile_count"], 1);
    let profile = &catalogue["profiles"][0];
    assert_eq!(profile["profile_id"], "biobrick-rfc10-part-sequence-2007");
    let parts = json!([{"id":"selected-part","document":document,"start":4,"length":formatted.len(),"reverse":false,"slot":"noncoding","body_start":22,"body_end":34}]);
    let arguments = json!({"profile_id":profile["profile_id"],"profile_identity":profile["profile_identity"],"parts":parts,"limits":collagene_assembly_standards::assessment::Limits::default(),"max_findings":100});
    let result = call(&client, "check_assembly_standard", arguments.clone()).await;
    assert_eq!(result["outcome"], "Violations");
    assert_eq!(result["complete"], true);
    assert_eq!(result["sources"][0]["revision"], molecule.revision());
    assert!(
        result["findings"]
            .as_array()
            .unwrap()
            .iter()
            .any(|f| f["finding"]["kind"]["AvoidedSite"]["rule_id"] == "noti"
                && f["source_spans"] == json!([[28, 36]]))
    );
    let report: standards::Report =
        serde_json::from_value(artifact(&client, &result["artifact_id"]).await).unwrap();
    let request = standards::Request {
        profile_id: profile["profile_id"].as_str().unwrap().into(),
        profile_identity: profile["profile_identity"].as_str().unwrap().into(),
        parts: vec![standards::SourcePart {
            id: "selected-part".into(),
            molecule: molecule.clone(),
            start: 4,
            length: formatted.len(),
            reverse: false,
            slot: "noncoding".into(),
            body_start: 22,
            body_end: 34,
        }],
        limits: collagene_assembly_standards::assessment::Limits::default(),
    };
    report.validate_for(&request).unwrap();
    report
        .validate_binding(result["context_revision"].as_str().unwrap())
        .unwrap();
    let mut stale_profile = arguments.clone();
    stale_profile["profile_identity"] = json!("0".repeat(64));
    reject(&client, "check_assembly_standard", stale_profile).await;
    let mut unknown = arguments.clone();
    unknown["profile_id"] = json!("unverified-moclo");
    reject(&client, "check_assembly_standard", unknown).await;
    let mut forged_source = arguments.clone();
    forged_source["parts"][0]["document"]["revision"] = json!("0".repeat(64));
    reject(&client, "check_assembly_standard", forged_source).await;
    let mut inline_source = arguments;
    inline_source["parts"][0]["molecule"] = json!(molecule);
    reject(&client, "check_assembly_standard", inline_source).await;
    assert_eq!(call(&client, "list_documents", json!({})).await, before);
    client.cancel().await.unwrap();
    service.await.unwrap();
}

#[tokio::test]
async fn collection_search_pages_all_bound_sites_and_keeps_mismatch_thermodynamics_unavailable() {
    let molecule = Molecule::dna(
        "two literal sites",
        "TTTACGTCAGTTCGAGGGACGTCAGTTCGACCC",
        Topology::Linear,
    )
    .unwrap();
    let (client, service, document) = setup(&molecule).await;
    let before = call(&client, "list_documents", json!({})).await;
    let collection = collagene_oligo_binding::Collection {
        oligos: vec![collagene_oligo_binding::Oligo {
            id: "tail-four".into(),
            name: "declared tail".into(),
            ordered: Sequence::dna("GGGGACGTCAGTTCGA").unwrap(),
            declared_tail_length: 4,
        }],
    };
    let options = collagene_oligo_binding::Options {
        minimum_annealing_length: 12,
        maximum_extra_unpaired: 0,
        maximum_mismatches: 0,
        minimum_three_prime_pairing: 12,
        match_score: 2,
        mismatch_score: -3,
        minimum_score: 24,
    };
    let budget = collagene_oligo_binding::BudgetLimits {
        work: 1_000_000,
        hits: 100,
        output_charge: 1024 * 1024,
    };
    let conditions = collagene_core::cloning::pcr::primers::ThermoConditions::default();
    let arguments = json!({"document":document,"collection":collection,"options":options,"budget":budget,"conditions":conditions,"max_hits":1});
    let result = call(&client, "search_oligo_collection", arguments.clone()).await;
    assert_eq!(result["complete"], true);
    assert_eq!(result["hits"].as_array().unwrap().len(), 1);
    assert_eq!(result["done"], false);
    assert_eq!(result["hits"][0]["declared_tail_length"], 4);
    assert!(
        result["hits"][0]["thermodynamics"]["annealing"]
            .get("Available")
            .is_some()
    );
    assert!(
        result["hits"][0]["thermodynamics"]["whole_oligo"]
            .get("Available")
            .is_some()
    );
    let full = artifact(&client, &result["artifact_id"]).await;
    for start in [3, 18] {
        assert!(
            full["search"]["hits"]
                .as_array()
                .unwrap()
                .iter()
                .any(|h| h["orientation"] == "Forward"
                    && h["source_start"] == start
                    && h["length"] == 12
                    && h["annealing_start"] == 4)
        );
    }
    let report: collection_binding::Report = serde_json::from_value(full).unwrap();
    let request = collection_binding::Request {
        source: molecule.clone(),
        collection,
        options,
        budget,
        conditions,
    };
    report.validate_for(&request).unwrap();
    report
        .validate_binding(result["context_revision"].as_str().unwrap())
        .unwrap();
    let mut second = arguments.clone();
    second["hit_start"] = json!(1);
    let next = call(&client, "search_oligo_collection", second).await;
    assert_eq!(next["artifact_id"], result["artifact_id"]);
    let mut mismatched = arguments.clone();
    mismatched["collection"]["oligos"][0]["ordered"] =
        json!(Sequence::dna("GGGGACGTCAGTTCGG").unwrap());
    mismatched["options"]["maximum_mismatches"] = json!(1);
    mismatched["options"]["minimum_three_prime_pairing"] = json!(0);
    mismatched["options"]["minimum_score"] = json!(19);
    mismatched["max_hits"] = json!(100);
    let mismatch = call(&client, "search_oligo_collection", mismatched).await;
    assert!(mismatch["hits"].as_array().unwrap().iter().any(|h| {
        h["maximum_mismatches"] == 1
            && h["thermodynamics"]["annealing"]
                .get("Unavailable")
                .is_some()
    }));
    let mut limited = arguments.clone();
    limited["budget"]["work"] = json!(
        collagene_oligo_binding::Budget::admission_cost(&molecule.sequence, &request.collection)
            .unwrap()
            .work
    );
    assert_eq!(
        call(&client, "search_oligo_collection", limited).await["complete"],
        false
    );
    let mut below_admission = arguments.clone();
    below_admission["budget"]["work"] = json!(0);
    reject(&client, "search_oligo_collection", below_admission).await;
    let mut stale = arguments.clone();
    stale["document"]["revision"] = json!("0".repeat(64));
    reject(&client, "search_oligo_collection", stale).await;
    let mut unsupported = arguments;
    unsupported["options"]["gap_open"] = json!(-2);
    reject(&client, "search_oligo_collection", unsupported).await;
    assert_eq!(call(&client, "list_documents", json!({})).await, before);
    client.cancel().await.unwrap();
    service.await.unwrap();
}
