//! Independent, hand-calculated scientific fixtures through real MCP JSON-RPC.
//! These tests exercise bindings and paginated evidence, not just core functions.
use base64::{Engine, engine::general_purpose::STANDARD};
use collagene_core::{
    Feature, Location, Molecule, Topology,
    alignment::{AlignmentMode, AlignmentOptions, ComparisonOptions},
    transfer::{PlanningWorkspace, save_workspace},
    translation::{IncompletePolicy, OrfOptions, TranslationOptions},
};
use collagene_mcp::CollageneServer;
use rmcp::{RoleClient, ServiceExt, model::CallToolRequestParams, service::RunningService};
use serde_json::{Value, json};

type Client = RunningService<RoleClient, ()>;
struct Session {
    client: Client,
    server: tokio::task::JoinHandle<()>,
}
impl Session {
    async fn new() -> Self {
        let (a, b) = tokio::io::duplex(65536);
        let server = tokio::spawn(async {
            CollageneServer::default()
                .serve(a)
                .await
                .unwrap()
                .waiting()
                .await
                .unwrap();
        });
        Self {
            client: ().serve(b).await.unwrap(),
            server,
        }
    }
    async fn call(&self, name: &'static str, args: Value) -> Value {
        let result = self
            .client
            .call_tool(
                CallToolRequestParams::new(name).with_arguments(args.as_object().unwrap().clone()),
            )
            .await
            .unwrap();
        assert_ne!(result.is_error, Some(true), "{name}: {:?}", result.content);
        result.structured_content.unwrap()
    }
    async fn reject(&self, name: &'static str, args: Value) {
        let result = self
            .client
            .call_tool(
                CallToolRequestParams::new(name).with_arguments(args.as_object().unwrap().clone()),
            )
            .await;
        assert!(
            result.is_err() || result.unwrap().is_error == Some(true),
            "accepted invalid {name}"
        );
    }
    async fn import(&self, documents: Vec<Molecule>) {
        let revisions: Vec<_> = documents.iter().map(Molecule::revision).collect();
        let bytes = save_workspace(&PlanningWorkspace {
            target_revision: revisions[0].clone(),
            documents,
        })
        .unwrap();
        let imported = self
            .call(
                "import_workspace",
                json!({"data_base64":STANDARD.encode(bytes)}),
            )
            .await;
        for revision in revisions {
            assert!(
                imported["documents"]
                    .as_array()
                    .unwrap()
                    .iter()
                    .any(|d| d["id"] == revision && d["revision"] == revision)
            );
        }
    }
    async fn artifact(&self, summary: &Value) -> Value {
        let id = summary["artifact_id"].as_str().unwrap();
        let mut text = String::new();
        let mut start = 0;
        loop {
            let chunk = self
                .call(
                    "read_artifact",
                    json!({"artifact_id":id,"start":start,"max_bytes":257}),
                )
                .await;
            assert_eq!(chunk["start"], start);
            assert_eq!(chunk["encoding"], "utf-8");
            assert_eq!(chunk["mime_type"], "application/json");
            text.push_str(chunk["text"].as_str().unwrap());
            let next = chunk["next_start"].as_u64().unwrap();
            assert!(next > start);
            start = next;
            if chunk["done"] == true {
                break;
            }
        }
        assert_eq!(blake3::hash(text.as_bytes()).to_hex().as_str(), id);
        let report: Value = serde_json::from_str(&text).unwrap();
        assert_eq!(
            report["engine_version"],
            collagene_core::verification::ENGINE_VERSION
        );
        report
    }
    async fn finish(self) {
        self.client.cancel().await.unwrap();
        self.server.await.unwrap();
    }
}
fn document(m: &Molecule) -> Value {
    json!({"id":m.revision(),"revision":m.revision()})
}
fn stale(mut reference: Value) -> Value {
    reference["revision"] = json!("stale-revision");
    reference
}
fn linear(name: &str, bases: &str) -> Molecule {
    Molecule::dna(name, bases, Topology::Linear).unwrap()
}
fn range(start: usize, end: usize) -> Location {
    Location::Range {
        start,
        end,
        before: false,
        after: false,
    }
}

#[tokio::test]
async fn translation_preserves_uncertainty_codes_and_full_revision_bound_evidence() {
    let s = Session::new().await;
    let ambiguous = linear("uncertain", "ATGGCNTAN");
    let codes = linear("code dependent", "ATATGAAGA");
    let long = linear("long preview", &"ATG".repeat(2105));
    s.import(vec![ambiguous.clone(), codes.clone(), long.clone()])
        .await;
    let opts = TranslationOptions::default();
    let args = json!({"document":document(&ambiguous),"options":opts});
    let summary = s.call("translate_sequence", args.clone()).await;
    assert_eq!(summary["protein_preview"], "MAX");
    assert_eq!(summary["amino_acids"], 3);
    assert_eq!(summary["unresolved_codons"], 1);
    assert_eq!(summary["possible_stop_codons"], 1);
    let full = s.artifact(&summary).await;
    assert_eq!(full["document"], document(&ambiguous));
    assert_eq!(full["translation"]["options"], json!(opts));
    assert_eq!(full["translation"]["ambiguous_codons"], 2);
    assert_eq!(
        s.call("translate_sequence", args).await["artifact_id"],
        summary["artifact_id"]
    );
    for (code, expected) in [(1, "I*R"), (2, "MW*")] {
        let options = TranslationOptions {
            genetic_code: code,
            ..opts
        };
        let translated = s
            .call(
                "translate_sequence",
                json!({"document":document(&codes),"options":options}),
            )
            .await;
        assert_eq!(translated["protein_preview"], expected);
    }
    let bounded = s
        .call(
            "translate_sequence",
            json!({"document":document(&long),"options":opts}),
        )
        .await;
    assert_eq!(bounded["preview_truncated"], true);
    assert_eq!(bounded["protein_preview"].as_str().unwrap().len(), 2000);
    assert_eq!(
        s.artifact(&bounded).await["translation"]["protein"],
        "M".repeat(2105)
    );
    s.reject(
        "translate_sequence",
        json!({"document":stale(document(&ambiguous)),"options":opts}),
    )
    .await;
    s.reject("translate_sequence", json!({"document":document(&ambiguous),"options":TranslationOptions{genetic_code:99,..opts}})).await;
    s.reject(
        "translate_sequence",
        json!({"document":document(&ambiguous),"options":TranslationOptions{frame:3,..opts}}),
    )
    .await;
    s.reject(
        "translate_sequence",
        json!({"document":document(&ambiguous),"options":opts,"verified":true}),
    )
    .await;
    s.finish().await;
}

#[tokio::test]
async fn cds_uses_compound_coordinates_and_refuses_unsupported_biology() {
    let s = Session::new().await;
    let mut m = Molecule::dna("origin CDS", "TAACCCATG", Topology::Circular).unwrap();
    let good = Feature {
        kind: "CDS".into(),
        location: Location::Join(vec![range(6, 9), range(0, 3)]),
        qualifiers: vec![
            ("transl_table".into(), Some("11".into())),
            ("translation".into(), Some("M".into())),
        ],
        needs_review: false,
    };
    let mut mismatch = good.clone();
    mismatch.qualifiers[1].1 = Some("W".into());
    let mut exception = good.clone();
    exception
        .qualifiers
        .push(("transl_except".into(), Some("(pos:7..9,aa:Sec)".into())));
    let mut noncoding = good.clone();
    noncoding.kind = "misc_feature".into();
    m.features = vec![good.clone(), mismatch, exception, noncoding];
    s.import(vec![m.clone()]).await;
    let summary = s
        .call(
            "translate_cds",
            json!({"document":document(&m),"feature_index":0}),
        )
        .await;
    assert_eq!(summary["protein_preview"], "M");
    assert_eq!(summary["code"], 11);
    assert_eq!(summary["complete_under_model"], true);
    assert_eq!(summary["stored_translation_matches"], true);
    assert_eq!(summary["location"], json!(good.location));
    let full = s.artifact(&summary).await;
    assert_eq!(full["document"], document(&m));
    assert_eq!(full["cds"]["source_revision"], m.revision());
    assert_eq!(full["cds"]["translation"]["protein"], "M*");
    assert_eq!(full["cds"]["has_terminal_stop"], true);
    let mismatch = s
        .call(
            "translate_cds",
            json!({"document":document(&m),"feature_index":1}),
        )
        .await;
    assert_eq!(mismatch["stored_translation_matches"], false);
    assert_eq!(mismatch["complete_under_model"], false);
    for index in [2, 3, 99] {
        s.reject(
            "translate_cds",
            json!({"document":document(&m),"feature_index":index}),
        )
        .await;
    }
    s.reject(
        "translate_cds",
        json!({"document":stale(document(&m)),"feature_index":0}),
    )
    .await;
    s.finish().await;
}

#[tokio::test]
async fn orf_origin_partial_and_search_limits_are_distinct_from_display_limits() {
    let s = Session::new().await;
    let circle = Molecule::dna("circular ORF", "TAACCCATG", Topology::Circular).unwrap();
    let partial = linear("uncertain stop", "ATGCCCTANAAA");
    let dense = linear("many ORFs", &"ATGAAATAA".repeat(120));
    let nested = linear("nested ORFs", "ATGATGTAA");
    s.import(vec![
        circle.clone(),
        partial.clone(),
        dense.clone(),
        nested.clone(),
    ])
    .await;
    let opts = OrfOptions {
        min_amino_acids: 1,
        both_strands: false,
        ..OrfOptions::default()
    };
    let summary = s
        .call(
            "find_orfs",
            json!({"document":document(&circle),"options":opts}),
        )
        .await;
    assert_eq!(summary["count"], 1);
    assert_eq!(summary["orfs"][0]["protein_preview"], "M");
    assert_eq!(summary["orfs"][0]["crosses_origin"], true);
    assert_eq!(summary["orfs"][0]["termination"], "StopCodon");
    assert_eq!(
        summary["orfs"][0]["location"],
        json!(Location::Join(vec![range(6, 9), range(0, 3)]))
    );
    let full = s.artifact(&summary).await;
    assert_eq!(full["document"], document(&circle));
    assert_eq!(full["orfs"]["source_revision"], circle.revision());
    let complete = s
        .call(
            "find_orfs",
            json!({"document":document(&partial),"options":opts}),
        )
        .await;
    assert_eq!(complete["count"], 0);
    let partial_opts = OrfOptions {
        incomplete: IncompletePolicy::IncludeAllPartials,
        ..opts.clone()
    };
    let uncertain = s
        .call(
            "find_orfs",
            json!({"document":document(&partial),"options":partial_opts}),
        )
        .await;
    assert!(
        uncertain["orfs"]
            .as_array()
            .unwrap()
            .iter()
            .any(|o| o["termination"] == "AmbiguousStop" && o["partial_3prime"] == true)
    );
    let bounded = s
        .call(
            "find_orfs",
            json!({"document":document(&nested),"options":OrfOptions{max_orfs:1,..opts.clone()}}),
        )
        .await;
    assert_eq!(bounded["count"], 1);
    assert_eq!(bounded["truncated_search"], true);
    assert_eq!(bounded["display_truncated"], false);
    let display = s
        .call(
            "find_orfs",
            json!({"document":document(&dense),"options":opts}),
        )
        .await;
    assert_eq!(display["count"], 120);
    assert_eq!(display["orfs"].as_array().unwrap().len(), 100);
    assert_eq!(display["display_truncated"], true);
    assert_eq!(display["truncated_search"], false);
    assert_eq!(
        s.artifact(&display).await["orfs"]["orfs"]
            .as_array()
            .unwrap()
            .len(),
        120
    );
    s.reject(
        "find_orfs",
        json!({"document":stale(document(&circle)),"options":opts}),
    )
    .await;
    s.finish().await;
}

#[tokio::test]
async fn pairwise_alignment_returns_exact_coordinates_gaps_and_both_source_bindings() {
    let s = Session::new().await;
    let r = linear("reference", "AAAACCCCGGGG");
    let q = linear("deletion", "AAAAGGGG");
    let flank = linear("flanked reference", "TTTACGTGGG");
    let insert = linear("insert", "ACGT");
    s.import(vec![r.clone(), q.clone(), flank.clone(), insert.clone()])
        .await;
    let opts = AlignmentOptions::default();
    let args = json!({"reference":document(&r),"query":document(&q),"options":opts});
    let summary = s.call("align_sequences", args.clone()).await;
    assert_eq!(summary["score"], 9); // 8 matches*2 - opening4 - extension3
    assert_eq!(summary["cigar_preview"], "4=4D4=");
    assert_eq!(summary["reference_range"], json!([0, 12]));
    assert_eq!(summary["query_range"], json!([0, 8]));
    let full = s.artifact(&summary).await;
    assert_eq!(full["reference"], document(&r));
    assert_eq!(full["query"], document(&q));
    assert_eq!(full["alignment"]["options"], json!(opts));
    assert_eq!(full["alignment"]["query_aligned"], "AAAA----GGGG");
    assert_eq!(
        full["alignment"]["differences"],
        json!([{"kind":"Deletion","reference_start":4,"reference_end":8,"query_start":4,"query_end":4,"reference_bases":"CCCC","query_bases":""}])
    );
    let embedded = s.call("align_sequences",json!({"reference":document(&flank),"query":document(&insert),"options":AlignmentOptions{mode:AlignmentMode::QueryInReference,..opts}})).await;
    assert_eq!(embedded["score"], 8);
    assert_eq!(embedded["reference_range"], json!([3, 7]));
    assert_eq!(embedded["query_range"], json!([0, 4]));
    assert_eq!(embedded["stats"]["reference_coverage"], 0.4);
    assert_eq!(embedded["stats"]["query_coverage"], 1.0);
    for key in ["reference", "query"] {
        let mut invalid = args.clone();
        invalid[key] = stale(invalid[key].clone());
        s.reject("align_sequences", invalid).await;
    }
    let mut limited = args;
    limited["options"]["max_cells"] = json!(10);
    s.reject("align_sequences", limited).await;
    s.finish().await;
}

#[tokio::test]
async fn circular_reverse_comparison_keeps_optimality_separate_from_enumeration() {
    let s = Session::new().await;
    let r = Molecule::dna("circle", "ACGATTCG", Topology::Circular).unwrap();
    let q = Molecule::dna("rotated reverse", "CGTCGAAT", Topology::Circular).unwrap();
    let edited = Molecule::dna("edited circle", "ATTCGACA", Topology::Circular).unwrap();
    s.import(vec![r.clone(), q.clone(), edited.clone()]).await;
    let opts = ComparisonOptions::default();
    let args = json!({"reference":document(&r),"query":document(&q),"options":opts});
    let summary = s.call("compare_molecules", args.clone()).await;
    assert_eq!(summary["score"], 16);
    assert_eq!(summary["reference_origin"], 3);
    assert_eq!(summary["query_reversed"], true);
    assert_eq!(summary["optimality_proven"], true);
    assert_eq!(summary["exhaustive_reference_origins"], false);
    let full = s.artifact(&summary).await;
    assert_eq!(full["reference"], document(&r));
    assert_eq!(full["query"], document(&q));
    assert_eq!(full["comparison"]["reference_revision"], r.revision());
    assert_eq!(full["comparison"]["query_revision"], q.revision());
    assert_eq!(full["comparison"]["alignment"]["cigar"], "8=");
    let bounded = s.call("compare_molecules",json!({"reference":document(&r),"query":document(&edited),"options":ComparisonOptions{max_origins:1,max_total_cells:81,..opts}})).await;
    assert_eq!(bounded["evaluated_candidates"], 1);
    assert_eq!(bounded["optimality_proven"], false);
    assert_eq!(bounded["exhaustive_reference_origins"], false);
    assert_eq!(bounded["exhaustive_query_orientations"], false);
    assert!(!bounded["warnings"].as_array().unwrap().is_empty());
    for key in ["reference", "query"] {
        let mut invalid = args.clone();
        invalid[key] = stale(invalid[key].clone());
        s.reject("compare_molecules", invalid).await;
    }
    s.finish().await;
}

#[tokio::test]
async fn scientific_tool_schemas_require_explicit_source_bindings_and_options() {
    let s = Session::new().await;
    let tools = s.client.list_all_tools().await.unwrap();
    for (name, required) in [
        ("translate_sequence", vec!["document", "options"]),
        ("translate_cds", vec!["document", "feature_index"]),
        ("find_orfs", vec!["document", "options"]),
        ("align_sequences", vec!["reference", "query", "options"]),
        ("compare_molecules", vec!["reference", "query", "options"]),
    ] {
        let tool = tools.iter().find(|t| t.name == name).unwrap();
        let schema = &tool.input_schema;
        assert_eq!(schema["additionalProperties"], false, "{name}");
        let fields = schema["required"].as_array().unwrap();
        for field in required {
            assert!(fields.contains(&json!(field)), "{name} requires {field}");
        }
        let serialized = serde_json::to_string(schema).unwrap();
        assert!(serialized.contains("revision") && serialized.contains("DocumentRef"));
        assert!(
            !schema["properties"]
                .as_object()
                .unwrap()
                .contains_key("verified")
        );
    }
    s.finish().await;
}
