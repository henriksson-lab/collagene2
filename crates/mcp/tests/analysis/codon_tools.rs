//! Actual SDK routing and immutable source/candidate evidence for codon design.
use collagene_core::sequence::codon_design::*;
use collagene_mcp::server::CollageneServer;
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
    let r = c
        .call_tool(
            CallToolRequestParams::new(name.to_owned())
                .with_arguments(args.as_object().unwrap().clone()),
        )
        .await
        .unwrap();
    assert_ne!(r.is_error, Some(true), "{name}:{:?}", r.content);
    r.structured_content.unwrap()
}
async fn reject(c: &Client, name: &str, args: Value) {
    let r = c
        .call_tool(
            CallToolRequestParams::new(name.to_owned())
                .with_arguments(args.as_object().unwrap().clone()),
        )
        .await;
    assert!(r.is_err() || r.unwrap().is_error == Some(true));
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
#[tokio::test]
async fn codon_design_and_supplied_candidate_checks_use_bound_sources_and_full_artifacts() {
    let (c, server) = setup().await;
    let docs = call(
        &c,
        "import_document",
        json!({"text":">required [molecule=protein]\nMK*\n>wrong\nATGGAGTAA\n"}),
    )
    .await;
    let reference = |name: &str| {
        let d = docs["documents"]
            .as_array()
            .unwrap()
            .iter()
            .find(|d| {
                d["name"]
                    .as_str()
                    .is_some_and(|n| n.split_whitespace().next() == Some(name))
            })
            .unwrap();
        json!({"id":d["id"],"revision":d["revision"]})
    };
    let source = reference("required");
    let options = DesignOptions {
        preferences: CodonPreferences::Explicit {
            weights: vec![CodonWeight {
                codon: "AAG".into(),
                weight: 100,
            }],
        },
        ..Default::default()
    };
    let result = call(
        &c,
        "design_codons",
        json!({"document":source,"options":options}),
    )
    .await;
    assert_eq!(result["source_revision"], source["revision"]);
    assert_eq!(result["status"], "CandidateFound");
    assert_eq!(result["candidate_sequence_preview"], "ATGAAGTAA");
    assert_eq!(result["candidate_check"]["satisfies"], true);
    let full: DesignReport =
        serde_json::from_value(artifact(&c, &result["artifact_id"]).await).unwrap();
    assert_eq!(
        full.candidate.as_ref().unwrap().sequence.as_str(),
        "ATGAAGTAA"
    );
    let candidate =
        json!({"id":result["candidate"]["id"],"revision":result["candidate"]["revision"]});
    let check = call(
        &c,
        "check_codon_design",
        json!({"source":source,"candidate":candidate,"options":options}),
    )
    .await;
    assert_eq!(check["satisfies"], true);
    let full: CheckReport =
        serde_json::from_value(artifact(&c, &check["artifact_id"]).await).unwrap();
    assert_eq!(full.target_protein, "MK*");
    assert_eq!(full.translated_protein, "MK*");
    let failed = call(
        &c,
        "check_codon_design",
        json!({"source":source,"candidate":reference("wrong"),"options":options}),
    )
    .await;
    assert_eq!(failed["satisfies"], false);
    assert!(failed["violations"].as_array().unwrap().len() > 0);
    let contextual = DesignOptions {
        left_context: "A".into(),
        forbidden_motifs: vec!["AATG".into()],
        ..Default::default()
    };
    let blocked = call(
        &c,
        "check_codon_design",
        json!({"source":source,"candidate":candidate,"options":contextual}),
    )
    .await;
    assert_eq!(
        blocked["satisfies"], false,
        "Context-bound constraints must be checked again"
    );
    let mut stale = source.clone();
    stale["revision"] = json!("0".repeat(64));
    reject(
        &c,
        "design_codons",
        json!({"document":stale,"options":options}),
    )
    .await;
    let mut unsupported = serde_json::to_value(&options).unwrap();
    unsupported["organism"] = json!("unproven");
    reject(
        &c,
        "design_codons",
        json!({"document":source,"options":unsupported}),
    )
    .await;
    let options = DesignOptions {
        preferences: CodonPreferences::Explicit {
            weights: vec![
                CodonWeight {
                    codon: "AAA".into(),
                    weight: 0,
                },
                CodonWeight {
                    codon: "AAG".into(),
                    weight: 0,
                },
            ],
        },
        ..Default::default()
    };
    let no = call(
        &c,
        "design_codons",
        json!({"document":source,"options":options}),
    )
    .await;
    assert!(no["candidate"].is_null());
    assert_eq!(no["status"], "Unsatisfiable");
    c.cancel().await.unwrap();
    server.await.unwrap();
}

#[tokio::test]
async fn report_capacity_failure_never_installs_an_orphan_candidate_document() {
    let (c, server) = setup().await;
    let docs = call(
        &c,
        "import_document",
        json!({"text":">required [molecule=protein]\nMK*\n>candidate\nATGAAATAA\n"}),
    )
    .await;
    let reference = |name: &str| {
        let d = docs["documents"]
            .as_array()
            .unwrap()
            .iter()
            .find(|d| {
                d["name"]
                    .as_str()
                    .is_some_and(|n| n.split_whitespace().next() == Some(name))
            })
            .unwrap();
        json!({"id":d["id"],"revision":d["revision"]})
    };
    let source = reference("required");
    let candidate = reference("candidate");
    for beam_width in 1..=128 {
        let options = DesignOptions {
            beam_width,
            ..Default::default()
        };
        let result = call(
            &c,
            "check_codon_design",
            json!({"source":source,"candidate":candidate,"options":options}),
        )
        .await;
        assert_eq!(result["satisfies"], true);
    }
    let before = call(&c, "list_documents", json!({"start":0,"limit":128})).await;
    let options = DesignOptions {
        preferences: CodonPreferences::Explicit {
            weights: vec![CodonWeight {
                codon: "AAG".into(),
                weight: 100,
            }],
        },
        ..Default::default()
    };
    reject(
        &c,
        "design_codons",
        json!({"document":source,"options":options}),
    )
    .await;
    assert_eq!(
        call(&c, "list_documents", json!({"start":0,"limit":128})).await,
        before
    );
    c.cancel().await.unwrap();
    server.await.unwrap();
}

fn mg1655_reference() -> Value {
    json!({"id":"ecoli-k12-mg1655-reference","identity":"52d52265b4f997e5f95c2d4178bcfc575e5002bef278a7a327a9efe6f8d61f43"})
}
#[tokio::test]
async fn curated_host_tools_expose_pinned_observed_counts_and_original_source_identity() {
    let (c, server) = setup().await;
    let tools = c.list_tools(Default::default()).await.unwrap().tools;
    assert_eq!(tools.len(), collagene_mcp::CURRENT_TOOL_COUNT);
    for name in [
        "list_codon_hosts",
        "codon_host_profile",
        "design_codons",
        "check_codon_design",
    ] {
        assert_eq!(tools.iter().filter(|tool| tool.name == name).count(), 1);
    }
    let listed = call(&c, "list_codon_hosts", json!({})).await;
    assert_eq!(listed["expression_prediction"], false);
    assert_eq!(listed["calibrated"], false);
    let descriptor = listed["hosts"]
        .as_array()
        .unwrap()
        .iter()
        .find(|host| host["reference"] == mg1655_reference())
        .unwrap();
    assert_eq!(descriptor["assembly"], "GCF_000005845.2");
    assert_eq!(descriptor["genetic_code"], 11);
    assert_eq!(descriptor["admitted_loci"], 4294);
    assert_eq!(descriptor["excluded_loci"], 24);
    assert_eq!(descriptor["internal_codons"], 1_322_130);
    let selected = call(
        &c,
        "codon_host_profile",
        json!({"reference":mg1655_reference()}),
    )
    .await;
    assert_eq!(&selected["host"], descriptor);
    assert_eq!(
        selected["preferences"],
        json!({"ReferenceHost":{"reference":mg1655_reference()}})
    );
    let full = artifact(&c, &selected["artifact_id"]).await;
    assert_eq!(full["profile"]["identity"], mg1655_reference()["identity"]);
    assert_eq!(full["profile"]["observed_loci"], 4318);
    assert_eq!(
        full["profile"]["admitted_loci_blake3"],
        "fac662e86e6aa59b35869e05d31940dfdfaf417d0cfec0e67852c678f0340357"
    );
    assert_eq!(
        full["profile"]["sources"][0]["blake3"],
        "ddad8311d8431b59019fe664542829bb3ad841bb2a048b218cd3b9e6fa6bc1db"
    );
    assert_eq!(
        full["profile"]["sources"][0]["expanded_blake3"],
        "19821f380b04476aafc0c15ba20db7131b4e569742d7a8f16ce4e3a2a2eaf58d"
    );
    assert_eq!(
        full["profile"]["sources"][1]["blake3"],
        "741a2f1ee0efea79ed0076a0b51acd6a76ce42aa92e055cdce56bb6f24ed0598"
    );
    assert!(
        full["profile"]["sources"][0]["url"]
            .as_str()
            .unwrap()
            .ends_with("GCF_000005845.2_ASM584v2_cds_from_genomic.fna.gz")
    );
    let order = full["codon_order"].as_array().unwrap();
    assert_eq!(order.len(), 64);
    assert_eq!(full["resolved_preferences"].as_array().unwrap().len(), 64);
    let at = |codon: &str| {
        order
            .iter()
            .position(|value| value.as_str() == Some(codon))
            .unwrap()
    };
    for (codon, count) in [
        ("AAA", 44804),
        ("AAG", 13583),
        ("ATG", 33272),
        ("CTG", 70911),
        ("CGC", 29500),
        ("TGG", 20224),
    ] {
        assert_eq!(
            full["profile"]["counts"][at(codon)],
            count,
            "{codon} internal count"
        );
    }
    assert_eq!(full["profile"]["initiation_counts"][at("ATG")], 3870);
    assert_eq!(full["profile"]["initiation_counts"][at("CTG")], 2);
    for (codon, count) in [("TAA", 2749), ("TAG", 305), ("TGA", 1240)] {
        assert_eq!(full["profile"]["counts"][at(codon)], 0);
        assert_eq!(full["profile"]["terminal_stop_counts"][at(codon)], count);
        assert_eq!(full["resolved_preferences"][at(codon)]["weight"], 1);
    }
    assert_eq!(full["resolved_preferences"][at("AAA")]["weight"], 1_000_000);
    assert_eq!(
        full["resolved_preferences"][at("AAG")]["weight"],
        1 + 13583u64 * 999999 / 44804
    );
    let counts: u64 = full["profile"]["counts"]
        .as_array()
        .unwrap()
        .iter()
        .map(|value| value.as_u64().unwrap())
        .sum();
    assert_eq!(counts, 1_322_130);
    let mut stale = mg1655_reference();
    stale["identity"] = json!("0".repeat(64));
    reject(&c, "codon_host_profile", json!({"reference":stale})).await;
    let mut unknown = mg1655_reference();
    unknown["id"] = json!("unregistered-host");
    reject(&c, "codon_host_profile", json!({"reference":unknown})).await;
    c.cancel().await.unwrap();
    server.await.unwrap();
}

#[tokio::test]
async fn named_host_preferences_bind_design_and_independent_candidate_check_without_expression_claims()
 {
    let (c, server) = setup().await;
    let docs = call(
        &c,
        "import_document",
        json!({"text":">required [molecule=protein]\nMK*\n>wrong\nATGGAGTAA\n"}),
    )
    .await;
    let reference = |name: &str| {
        let doc = docs["documents"]
            .as_array()
            .unwrap()
            .iter()
            .find(|doc| doc["name"].as_str().unwrap().split_whitespace().next() == Some(name))
            .unwrap();
        json!({"id":doc["id"],"revision":doc["revision"]})
    };
    let source = reference("required");
    let before = call(&c, "list_documents", json!({})).await;
    let profile = call(
        &c,
        "codon_host_profile",
        json!({"reference":mg1655_reference()}),
    )
    .await;
    let mut options = serde_json::to_value(DesignOptions {
        genetic_code: 11,
        ..Default::default()
    })
    .unwrap();
    options["preferences"] = profile["preferences"].clone();
    let designed = call(
        &c,
        "design_codons",
        json!({"document":source,"options":options}),
    )
    .await;
    assert_eq!(designed["candidate_sequence_preview"], "ATGAAATAA");
    assert_eq!(designed["candidate_check"]["satisfies"], true);
    assert_eq!(designed["preferences"], profile["preferences"]);
    let full: DesignReport =
        serde_json::from_value(artifact(&c, &designed["artifact_id"]).await).unwrap();
    assert_eq!(
        serde_json::to_value(&full.options).unwrap()["preferences"],
        profile["preferences"]
    );
    assert_eq!(full.source_revision, source["revision"].as_str().unwrap());
    let candidate =
        json!({"id":designed["candidate"]["id"],"revision":designed["candidate"]["revision"]});
    let checked = call(
        &c,
        "check_codon_design",
        json!({"source":source,"candidate":candidate,"options":options}),
    )
    .await;
    assert_eq!(checked["satisfies"], true);
    assert_eq!(checked["options_revision"], designed["options_revision"]);
    assert_eq!(checked["preferences"], profile["preferences"]);
    let wrong = call(
        &c,
        "check_codon_design",
        json!({"source":source,"candidate":reference("wrong"),"options":options}),
    )
    .await;
    assert_eq!(wrong["satisfies"], false);
    let mut contextual = options.clone();
    contextual["left_context"] = json!("A");
    contextual["forbidden_motifs"] = json!(["AATG"]);
    let no = call(
        &c,
        "check_codon_design",
        json!({"source":source,"candidate":candidate,"options":contextual}),
    )
    .await;
    assert_eq!(no["satisfies"], false);
    let mut mismatch = options.clone();
    mismatch["genetic_code"] = json!(1);
    reject(
        &c,
        "design_codons",
        json!({"document":source,"options":mismatch}),
    )
    .await;
    reject(
        &c,
        "check_codon_design",
        json!({"source":source,"candidate":candidate,"options":mismatch}),
    )
    .await;
    let mut altered = options.clone();
    altered["preferences"]["ReferenceHost"]["reference"]["identity"] = json!("0".repeat(64));
    reject(
        &c,
        "design_codons",
        json!({"document":source,"options":altered}),
    )
    .await;
    let mut empty = options.clone();
    empty["preferences"] = json!({"Explicit":{"weights":[]}});
    reject(
        &c,
        "design_codons",
        json!({"document":source,"options":empty}),
    )
    .await;
    let mut explicit = options.clone();
    explicit["preferences"] = json!({"Explicit":{"weights":artifact(&c,&profile["artifact_id"]).await["resolved_preferences"]}});
    let explicit_check = call(
        &c,
        "check_codon_design",
        json!({"source":source,"candidate":candidate,"options":explicit}),
    )
    .await;
    assert_eq!(explicit_check["satisfies"], true);
    assert_ne!(
        explicit_check["options_revision"], checked["options_revision"],
        "equal arithmetic scores do not erase named-host provenance"
    );
    // Immutable source/candidate contents are still the same imported revisions.
    let after = call(&c, "list_documents", json!({})).await;
    for original in before["documents"].as_array().unwrap() {
        assert!(after["documents"].as_array().unwrap().contains(original));
    }
    c.cancel().await.unwrap();
    server.await.unwrap();
}
