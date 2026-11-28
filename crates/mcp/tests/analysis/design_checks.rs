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
    let installed = collagene_assembly_standards::profile::BUILTIN_IDS;
    assert_eq!(catalogue["profile_count"], installed.len());
    let advertised = catalogue["profiles"]
        .as_array()
        .unwrap()
        .iter()
        .map(|profile| profile["profile_id"].as_str().unwrap())
        .collect::<std::collections::BTreeSet<_>>();
    assert_eq!(advertised, installed.iter().copied().collect());
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
async fn collection_search_pages_bound_sites_and_distinguishes_supported_mismatch_thermodynamics() {
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
    let mut internal = arguments.clone();
    internal["collection"]["oligos"][0]["ordered"] =
        json!(Sequence::dna("GGGGACGTAAGTTCGA").unwrap());
    internal["options"]["maximum_mismatches"] = json!(1);
    internal["options"]["minimum_three_prime_pairing"] = json!(4);
    internal["options"]["minimum_score"] = json!(19);
    internal["max_hits"] = json!(100);
    let internal_result = call(&client, "search_oligo_collection", internal.clone()).await;
    assert!(
        internal_result["hits"]
            .as_array()
            .unwrap()
            .iter()
            .any(|hit| {
                hit["thermodynamics"]["annealing"]
                    .get("MismatchAvailable")
                    .is_some()
            })
    );
    let internal_report: collection_binding::Report =
        serde_json::from_value(artifact(&client, &internal_result["artifact_id"]).await).unwrap();
    let mut internal_request = request.clone();
    internal_request.collection = serde_json::from_value(internal["collection"].clone()).unwrap();
    internal_request.options = serde_json::from_value(internal["options"].clone()).unwrap();
    internal_report.validate_for(&internal_request).unwrap();
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

fn declared_chain_options(molecule: &Molecule) -> construct_qc::Options {
    use construct_qc::{
        element_chain::{ChainElement, ElementChain, ElementRole},
        required_parts::{RequiredOrientation, RequiredPart},
    };
    let mut options = construct_qc::Options::default();
    options.annotations = molecule
        .features
        .iter()
        .enumerate()
        .map(|(i, feature)| construct_qc::AnnotationAssumption {
            feature_index: i,
            source_revision: molecule.revision(),
            annotation_revision: construct_qc::annotation_revision(feature).unwrap(),
            rationale: "Literal SDK fixture roles are explicit assumptions, not verified function."
                .into(),
        })
        .collect();
    for (id, dna) in [
        ("upstream-reference", "CC"),
        ("downstream-reference", "GG"),
        ("ambiguous-motif", "ATG"),
    ] {
        let sequence = Sequence::dna(dna).unwrap();
        options.references.push(construct_qc::ReferenceElement {
            id: id.into(),
            sequence_revision: construct_qc::reference_revision(&sequence).unwrap(),
            sequence,
            provenance: construct_qc::Provenance {
                source_url: "https://example.test/literal-sdk-fixtures".into(),
                accession: id.into(),
                retrieved_on: "2026-11-29".into(),
                license: "MIT".into(),
            },
        });
    }
    options.element_chains.push(ElementChain {
        id: "declared-forward-chain".into(),
        cds_feature: 0,
        upstream: Some(ChainElement {
            feature_index: 1,
            reference_id: Some("upstream-reference".into()),
            role: ElementRole::Promoter,
        }),
        downstream: Some(ChainElement {
            feature_index: 2,
            reference_id: Some("downstream-reference".into()),
            role: ElementRole::Terminator,
        }),
        max_upstream_distance: Some(0),
        max_downstream_distance: Some(0),
    });
    options.required_parts.push(RequiredPart {
        id: "declared-literal-marker".into(),
        role: "A literal interval label, not an origin/marker function claim".into(),
        reference_id: "upstream-reference".into(),
        minimum_copies: 1,
        maximum_copies: Some(1),
        orientation: RequiredOrientation::Forward,
    });
    options
}

#[tokio::test]
async fn construct_element_sdk_retains_literal_counts_and_never_promotes_incomplete_search_to_absence()
 {
    use construct_qc::required_parts::{RequiredOrientation, RequiredPart};
    let mut molecule = Molecule::dna(
        "forward chain plus unresolved motif",
        "CCATGAAATAAGGATN",
        Topology::Linear,
    )
    .unwrap();
    molecule.features = [
        ("CDS", 2, 11),
        ("misc_feature", 0, 2),
        ("misc_feature", 11, 13),
    ]
    .into_iter()
    .map(|(kind, start, end)| Feature {
        kind: kind.into(),
        location: Location::range(start, end),
        qualifiers: vec![],
        needs_review: false,
    })
    .collect();
    let options = declared_chain_options(&molecule);
    assert_eq!(options.schema_version, 3);
    let (client, service, document) = setup(&molecule).await;
    let before = call(&client, "list_documents", json!({})).await;
    let original = call(
        &client,
        "check_construct_qc",
        json!({"document":document,"options":options,"max_findings":100}),
    )
    .await;
    let full = artifact(&client, &original["artifact_id"]).await;
    for container in [&original, &full] {
        for forbidden in [
            "overall_pass",
            "overall_verdict",
            "verdict",
            "success",
            "pass",
        ] {
            assert!(
                container.get(forbidden).is_none(),
                "No global biological verdict: {forbidden}"
            );
        }
        assert_eq!(container["statement"], construct_qc::STATEMENT);
    }
    let report: construct_qc::Report = serde_json::from_value(full).unwrap();
    report.validate_for(&molecule, &options).unwrap();
    report
        .validate_binding(original["context_revision"].as_str().unwrap())
        .unwrap();
    assert_eq!(
        report.element_chains[0].verdict,
        construct_qc::Verdict::Pass
    );
    assert_eq!(report.element_chains[0].upstream.distance_bases, Some(0));
    assert_eq!(report.element_chains[0].downstream.distance_bases, Some(0));
    assert_eq!(
        report.required_parts[0].verdict,
        construct_qc::Verdict::Pass
    );
    assert_eq!(
        (
            report.required_parts[0].certain_copies,
            report.required_parts[0].possible_copies
        ),
        (1, 1)
    );
    let search = &report.required_part_searches[report.required_parts[0].search_index];
    assert!(search.complete);
    assert!(
        search
            .matches
            .iter()
            .any(|hit| hit.reverse == Some(false) && hit.coordinates == [Location::range(0, 2)])
    );

    let mut limited_options = options.clone();
    limited_options.required_part_limits.max_work = 1;
    let limited = call(
        &client,
        "check_construct_qc",
        json!({"document":document,"options":limited_options,"max_findings":100}),
    )
    .await;
    let limited_report: construct_qc::Report =
        serde_json::from_value(artifact(&client, &limited["artifact_id"]).await).unwrap();
    limited_report
        .validate_for(&molecule, &limited_options)
        .unwrap();
    assert!(!limited_report.required_part_searches[0].complete);
    assert_eq!(limited_report.required_parts[0].possible_copies, 0);
    assert_eq!(
        limited_report.required_parts[0].verdict,
        construct_qc::Verdict::Indeterminate
    );
    assert_eq!(
        limited_report.element_chains[0].verdict,
        construct_qc::Verdict::Pass
    );
    assert_ne!(limited["context_revision"], original["context_revision"]);

    let mut ambiguous_options = options.clone();
    ambiguous_options.required_parts.push(RequiredPart {
        id: "uncertain-forward-motif".into(),
        role: "Literal uncertainty fixture".into(),
        reference_id: "ambiguous-motif".into(),
        minimum_copies: 1,
        maximum_copies: Some(2),
        orientation: RequiredOrientation::Forward,
    });
    let ambiguous = call(
        &client,
        "check_construct_qc",
        json!({"document":document,"options":ambiguous_options,"max_findings":100}),
    )
    .await;
    let ambiguous_report: construct_qc::Report =
        serde_json::from_value(artifact(&client, &ambiguous["artifact_id"]).await).unwrap();
    ambiguous_report
        .validate_for(&molecule, &ambiguous_options)
        .unwrap();
    let uncertain = &ambiguous_report.required_parts[1];
    assert!(ambiguous_report.required_part_searches[uncertain.search_index].complete);
    assert_eq!(
        (uncertain.certain_copies, uncertain.possible_copies),
        (1, 2)
    );
    assert_eq!(uncertain.verdict, construct_qc::Verdict::Indeterminate);

    let mut changed_options = options.clone();
    changed_options.element_chains[0].upstream = None;
    let changed = call(
        &client,
        "check_construct_qc",
        json!({"document":document,"options":changed_options,"max_findings":100}),
    )
    .await;
    let changed_report: construct_qc::Report =
        serde_json::from_value(artifact(&client, &changed["artifact_id"]).await).unwrap();
    changed_report
        .validate_for(&molecule, &changed_options)
        .unwrap();
    assert_eq!(
        changed_report.element_chains[0].verdict,
        construct_qc::Verdict::Indeterminate
    );
    assert!(
        report
            .validate_binding(changed["context_revision"].as_str().unwrap())
            .is_err()
    );
    assert!(report.validate_for(&molecule, &changed_options).is_err());
    let mut changed_source = molecule.clone();
    changed_source.sequence = Sequence::dna("CCATGAAATAAGGATT").unwrap();
    assert!(report.validate_for(&changed_source, &options).is_err());
    let mut stale_document = document.clone();
    stale_document["revision"] = json!(changed_source.revision());
    reject(
        &client,
        "check_construct_qc",
        json!({"document":stale_document,"options":options}),
    )
    .await;
    let mut old_format = serde_json::to_value(&options).unwrap();
    old_format["schema_version"] = 1.into();
    reject(
        &client,
        "check_construct_qc",
        json!({"document":document,"options":old_format}),
    )
    .await;
    assert_eq!(call(&client, "list_documents", json!({})).await, before);
    client.cancel().await.unwrap();
    service.await.unwrap();
}

#[tokio::test]
async fn declared_context_and_orf_sdk_keep_individual_evidence_and_exact_source_scope() {
    use collagene_core::sequence::translation::{IncompletePolicy, StartPolicy};
    let mut molecule = Molecule::dna(
        "literal initiation context",
        "ACCATGGAATAA",
        Topology::Linear,
    )
    .unwrap();
    molecule.features.push(Feature {
        kind: "CDS".into(),
        location: Location::range(3, 12),
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
                "Explicit declared CDS role and contiguous transcript context for this fixture."
                    .into(),
        });
    options
        .initiation_contexts
        .push(construct_qc::initiation_context::ContextClaim {
            cds_feature: 0,
            profile: construct_qc::initiation_context::kozak_1986_pattern().unwrap(),
            rationale: "Declared mammalian sequence comparison; not an expression forecast.".into(),
        });
    options.orf_scan = Some(construct_qc::unannotated_orfs::ScanOptions {
        start_policy: StartPolicy::AtgOnly,
        incomplete: IncompletePolicy::CompleteOnly,
        minimum_amino_acids: 2,
        max_candidates: 8,
        both_strands: false,
    });
    let (client, service, document) = setup(&molecule).await;
    let before = call(&client, "list_documents", json!({})).await;
    let result = call(
        &client,
        "check_construct_qc",
        json!({"document":document,"options":options,"max_findings":100}),
    )
    .await;
    let full: construct_qc::Report =
        serde_json::from_value(artifact(&client, &result["artifact_id"]).await).unwrap();
    full.validate_for(&molecule, &options).unwrap();
    assert_eq!(
        full.initiation_contexts[0].observed.as_deref(),
        Some("ACCATGG")
    );
    assert_eq!(
        full.initiation_contexts[0].verdict,
        construct_qc::Verdict::Pass
    );
    let scan = full.orf_scan.as_ref().unwrap();
    assert!(scan.enumeration_complete);
    assert_eq!(scan.candidates.len(), 1);
    assert_eq!(
        scan.candidates[0].annotation_relation,
        construct_qc::unannotated_orfs::AnnotationRelation::CoveredByDeclaredCds
    );
    assert!(
        full.findings
            .iter()
            .filter(|f| f.rule == construct_qc::Rule::UnannotatedOrf)
            .all(|f| f.verdict != construct_qc::Verdict::Pass
                && f.verdict != construct_qc::Verdict::Fail)
    );
    assert_eq!(call(&client, "list_documents", json!({})).await, before);
    client.cancel().await.unwrap();
    service.await.unwrap();
}
