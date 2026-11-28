//! Actual SDK routing and immutable source/candidate evidence for codon design.
use collagene_core::sequence::codon_design::*;
use collagene_mcp::server::CollageneServer;
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
        codon_weights: vec![CodonWeight {
            codon: "AAG".into(),
            weight: 100,
        }],
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
        codon_weights: vec![
            CodonWeight {
                codon: "AAA".into(),
                weight: 0,
            },
            CodonWeight {
                codon: "AAG".into(),
                weight: 0,
            },
        ],
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
        codon_weights: vec![CodonWeight {
            codon: "AAG".into(),
            weight: 100,
        }],
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
