//! Real SDK JSON-RPC exercises independent Type IIS, mismatch, inverse-PCR and gel fixtures.
use base64::{Engine, engine::general_purpose::STANDARD};
use collagene_core::documents::reference_edit::Operation as Edit;
use collagene_core::{
    cloning::{
        golden_gate::GoldenGateOptions,
        pcr::{
            mutagenesis::{InversePcrOptions, PcrOptions},
            primers::BindingOptions,
        },
    },
    documents::transfer::{PlanningWorkspace, save_workspace},
    output::gel::{Calibration, CalibrationPoint, DiagnosticOptions, Preparation},
    sequence::analysis::Primer,
};
use collagene_mcp::server::CollageneServer;
use collagene_model::{Molecule, Topology};
use rmcp::{RoleClient, ServiceExt, model::CallToolRequestParams, service::RunningService};
use serde_json::{Value, json};
struct Session {
    client: RunningService<RoleClient, ()>,
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
        let bytes = save_workspace(&PlanningWorkspace {
            target_revision: documents[0].revision(),
            documents,
        })
        .unwrap();
        self.call(
            "import_workspace",
            json!({"data_base64":STANDARD.encode(bytes)}),
        )
        .await;
    }
    async fn text(&self, id: &Value, mime: &str, encoding: &str) -> String {
        let mut text = String::new();
        let mut start = 0;
        loop {
            let chunk = self
                .call(
                    "read_artifact",
                    json!({"artifact_id":id,"start":start,"max_bytes":16381}),
                )
                .await;
            assert_eq!(chunk["start"], start);
            assert_eq!(chunk["mime_type"], mime);
            assert_eq!(chunk["encoding"], encoding);
            text.push_str(chunk["text"].as_str().unwrap());
            let next = chunk["next_start"].as_u64().unwrap();
            assert!(next > start);
            start = next;
            if chunk["done"] == true {
                break;
            }
        }
        text
    }
    async fn artifact(&self, summary: &Value) -> Value {
        let text = self
            .text(&summary["artifact_id"], "application/json", "utf-8")
            .await;
        assert_eq!(
            blake3::hash(text.as_bytes()).to_hex().as_str(),
            summary["artifact_id"].as_str().unwrap()
        );
        let value: Value = serde_json::from_str(&text).unwrap();
        assert_eq!(
            value["engine_version"],
            collagene_core::verification::ENGINE_VERSION
        );
        value
    }
    async fn pdf(&self, plan: &Value, status: &str, target: &Molecule) {
        let summary = self
            .call(
                "export_protocol_pdf",
                json!({"plan_revision":plan["plan_revision"]}),
            )
            .await;
        assert_eq!(summary["status"], status);
        assert_eq!(summary["target_revision"], target.revision());
        assert_eq!(summary["plan_revision"], plan["plan_revision"]);
        assert_eq!(
            summary["engine_version"],
            collagene_core::verification::ENGINE_VERSION
        );
        let text = self
            .text(&summary["artifact_id"], "application/pdf", "base64")
            .await;
        let bytes = STANDARD.decode(text).unwrap();
        assert!(bytes.starts_with(b"%PDF-"));
        assert_eq!(bytes.len() as u64, summary["bytes"].as_u64().unwrap());
        assert_eq!(
            blake3::hash(&bytes).to_hex().as_str(),
            summary["artifact_id"].as_str().unwrap()
        );
    }
    async fn finish(self) {
        self.client.cancel().await.unwrap();
        self.server.await.unwrap();
    }
}
fn reference(m: &Molecule) -> Value {
    json!({"id":m.revision(),"revision":m.revision()})
}
fn stale(m: &Molecule) -> Value {
    json!({"id":m.revision(),"revision":"stale"})
}
fn dna(name: &str, s: &str, t: Topology) -> Molecule {
    Molecule::dna(name, s, t).unwrap()
}
fn donor(name: &str, left: &str, payload: &str, right: &str) -> Molecule {
    dna(
        name,
        &format!("GGTCTCA{left}{payload}{right}AGAGACC"),
        Topology::Linear,
    )
}
fn calibration() -> Calibration {
    Calibration {
        label: "Independent synthetic test calibration".into(),
        points: vec![
            CalibrationPoint {
                size_bp: 1.,
                migration: 0.9,
            },
            CalibrationPoint {
                size_bp: 1000.,
                migration: 0.1,
            },
        ],
        resolution: 0.01,
    }
}

#[tokio::test]
async fn golden_gate_keeps_actual_strands_sources_completeness_and_separate_target_verification() {
    let s = Session::new().await;
    let a = donor("A", "AATG", "CCCC", "GCTT");
    let b = donor("B", "GCTT", "TTTT", "AATG");
    let target = dna("Golden Gate target", "AATGCCCCGCTTTTTT", Topology::Circular);
    let wrong = dna("Wrong target", "AATGCCACGCTTTTTT", Topology::Circular);
    s.import(vec![a.clone(), b.clone(), target.clone(), wrong.clone()])
        .await;
    let options = GoldenGateOptions::default();
    let unknown = s
        .call(
            "simulate_golden_gate",
            json!({"documents":[reference(&a),reference(&b)],"options":options}),
        )
        .await;
    assert_eq!(unknown["products_determined"], false);
    assert_eq!(unknown["candidate_count"], 0);
    assert_eq!(unknown["enumeration_complete"], false);
    let unknown_artifact = s.artifact(&unknown).await;
    assert_eq!(
        unknown_artifact["golden_gate"]["source_cleavage"]
            .as_array()
            .unwrap()
            .len(),
        2
    );
    assert!(
        unknown_artifact["golden_gate"]["fragments"]
            .as_array()
            .unwrap()
            .is_empty()
    );
    let generation = json!({"host":null,"dam":"Unmodified","dcm":"Unmodified","cpg":"Unmodified"});
    let inputs = json!([{"id":"A","document":reference(&a)},{"id":"B","document":reference(&b)}]);
    let request = json!({"name":"Declared unmodified reaction","rows":[
        {"role":"A","input_id":"A","generation":generation,"substrate":{"WholeSource":{"origin":null}}},
        {"role":"B","input_id":"B","generation":generation,"substrate":{"WholeSource":{"origin":null}}}
    ],"options":options});
    let response = s.call("inspect_golden_gate_reaction",json!({"inputs":inputs,"target":reference(&target),"origin_workspaces":[],"request":request})).await;
    assert_eq!(response["products_determined"], true);
    assert_eq!(response["candidate_count"], 2);
    assert_eq!(response["stable_physical_product_count"], 1);
    assert_eq!(response["enumeration_complete"], true);
    assert_eq!(response["exhaustive_mixture_search"], false);
    let artifact_text = s
        .text(&response["artifact_id"], "application/json", "utf-8")
        .await;
    assert_eq!(
        blake3::hash(artifact_text.as_bytes()).to_hex().as_str(),
        response["artifact_id"].as_str().unwrap()
    );
    let artifact: Value = serde_json::from_str(&artifact_text).unwrap();
    assert_eq!(
        artifact["report"]["verification"]["engine_version"],
        collagene_core::verification::ENGINE_VERSION
    );
    assert_eq!(artifact["inputs"], inputs);
    let pool = &artifact["report"]["pool"];
    let declared_a: Molecule =
        serde_json::from_value(artifact["report"]["verification"]["steps"][0]["product"].clone())
            .unwrap();
    assert_eq!(declared_a.sequence, a.sequence);
    assert_eq!(pool["sources"][0]["revision"], declared_a.revision());
    let fragments = pool["fragments"].as_array().unwrap();
    assert_eq!(fragments[1]["top_start"], 7);
    assert_eq!(fragments[1]["bottom_start"], 11);
    assert_eq!(fragments[1]["duplex"]["top"]["bases"], "AATGCCCC");
    assert_eq!(fragments[1]["duplex"]["bottom"]["bases"], "AAGCGGGG");
    let candidates = pool["candidates"].as_array().unwrap();
    let selected = candidates
        .iter()
        .position(|c| c["duplex"]["top"]["bases"] == "AATGCCCCGCTTTTTT")
        .unwrap();
    assert_eq!(candidates[selected]["restriction_site_free"], true);
    assert!(
        candidates[selected]["recut_sites"]
            .as_array()
            .unwrap()
            .is_empty()
    );
    for bad in [[stale(&a), reference(&b)], [reference(&a), stale(&b)]] {
        s.reject(
            "simulate_golden_gate",
            json!({"documents":bad,"options":options}),
        )
        .await;
    }
    let limited = GoldenGateOptions {
        max_states: 1,
        ..options.clone()
    };
    let mut capped_request = request.clone();
    capped_request["options"] = json!(limited);
    let truncated = s.call("inspect_golden_gate_reaction",json!({"inputs":inputs,"target":reference(&target),"origin_workspaces":[],"request":capped_request})).await;
    assert_eq!(truncated["enumeration_complete"], false);
    assert!(
        !truncated["truncation_reasons"]
            .as_array()
            .unwrap()
            .is_empty()
    );
    let incomplete=s.call("propose_plan",json!({"name":"Capped Type IIS search is not a verified design","target":reference(&target),"inputs":[{"id":"A","document":reference(&a)},{"id":"B","document":reference(&b)}],"steps":[{"DeclareSourceModifications":{"template":{"Input":"A"},"assumptions":generation}},{"DeclareSourceModifications":{"template":{"Input":"B"},"assumptions":generation}},{"GoldenGate":{"templates":[{"Step":0},{"Step":1}],"options":limited,"candidate_index":0}}]})).await;
    assert_eq!(
        s.call(
            "verify_plan",
            json!({"plan_revision":incomplete["plan_revision"]})
        )
        .await["status"],
        "Incomplete"
    );
    s.pdf(&incomplete, "Incomplete", &target).await;
    let operation = json!({"GoldenGate":{"templates":[{"Step":0},{"Step":1}],"options":options,"candidate_index":selected}});
    let proposal = |target: &Molecule| json!({"name":"Explicit Type IIS geometry","target":reference(target),"inputs":[{"id":"A","document":reference(&a)},{"id":"B","document":reference(&b)}],"steps":[{"DeclareSourceModifications":{"template":{"Input":"A"},"assumptions":generation}},{"DeclareSourceModifications":{"template":{"Input":"B"},"assumptions":generation}},operation]});
    let plan = s.call("propose_plan", proposal(&target)).await;
    assert_eq!(
        s.call(
            "verify_plan",
            json!({"plan_revision":plan["plan_revision"]})
        )
        .await["status"],
        "ConsistentUnderAssumptions"
    );
    s.pdf(&plan, "ConsistentUnderAssumptions", &target).await;
    let badplan = s.call("propose_plan", proposal(&wrong)).await;
    assert_eq!(
        s.call(
            "verify_plan",
            json!({"plan_revision":badplan["plan_revision"]})
        )
        .await["status"],
        "Invalid"
    );
    s.pdf(&badplan, "Invalid", &wrong).await;
    assert_eq!(
        s.call("list_documents", json!({})).await["documents"]
            .as_array()
            .unwrap()
            .len(),
        4
    );
    s.finish().await;
}

#[tokio::test]
async fn mismatched_oligos_reconstruct_both_introduced_changes_and_capped_search_stays_incomplete()
{
    let s = Session::new().await;
    let source = dna(
        "Mismatched PCR source",
        "AGTCCGATGCTACGTAGGACTTCAGACCTGTA",
        Topology::Linear,
    );
    let target = dna(
        "Mutant target",
        "AGACCGATGCTACGTAGGACTTCAGACCTGGA",
        Topology::Linear,
    );
    s.import(vec![source.clone(), target.clone()]).await;
    let forward = Primer::new("F", "AGACCGAT", "").unwrap();
    let reverse = Primer::new("R", "TCCAGGTC", "").unwrap();
    let options = PcrOptions {
        binding: BindingOptions {
            max_mismatches: 1,
            min_three_prime_exact: 2,
            ..Default::default()
        },
        ..Default::default()
    };
    let query = json!({"document":reference(&source),"forward":forward,"reverse":reverse,"options":options});
    let response = s.call("simulate_mismatch_pcr", query.clone()).await;
    assert_eq!(response["product_count"], 1);
    assert_eq!(response["search_complete"], true);
    assert_eq!(response["enumeration_complete"], true);
    assert_eq!(response["products"][0]["change_count"], 2);
    let artifact = s.artifact(&response).await;
    assert_eq!(artifact["document"], reference(&source));
    assert_eq!(
        artifact["pcr"]["products"][0]["sequence"]["bases"],
        target.sequence.as_str()
    );
    let changes = artifact["pcr"]["products"][0]["substitutions"]
        .as_array()
        .unwrap();
    assert_eq!(
        changes
            .iter()
            .map(|c| c["source_position"].as_u64().unwrap())
            .collect::<Vec<_>>(),
        vec![2, 30]
    );
    assert_eq!(changes[0]["template_base"], "T");
    assert_eq!(changes[0]["product_base"], "A");
    assert_eq!(changes[1]["template_base"], "T");
    assert_eq!(changes[1]["product_base"], "G");
    let mut bad = query.clone();
    bad["document"] = stale(&source);
    s.reject("simulate_mismatch_pcr", bad).await;
    let mut capped = query.clone();
    // The known32-base product cannot fit this materialization budget; absence must
    // remain incomplete rather than becoming a false "no amplification" assertion.
    capped["options"]["max_product_bases"] = json!(1);
    let capped = s.call("simulate_mismatch_pcr", capped).await;
    assert_eq!(capped["search_complete"], false);
    assert_eq!(capped["enumeration_complete"], false);
    let plan=s.call("propose_plan",json!({"name":"Two explicit primer-encoded substitutions","target":reference(&target),"inputs":[{"id":"template","document":reference(&source)}],"steps":[{"PcrWithMismatches":{"template":{"Input":"template"},"forward":forward,"reverse":reverse,"options":options}}]})).await;
    let verified = s
        .call(
            "verify_plan",
            json!({"plan_revision":plan["plan_revision"]}),
        )
        .await;
    assert_eq!(verified["status"], "ConsistentUnderAssumptions");
    s.pdf(&plan, "ConsistentUnderAssumptions", &target).await;
    s.finish().await;
}

#[tokio::test]
async fn inverse_pcr_creates_an_explicit_reviewable_phosphorylation_and_ligation_plan() {
    let s = Session::new().await;
    let source = dna(
        "Inverse source",
        "AGTCCGATGCTACGTAGGACTTCAGACCTGTA",
        Topology::Circular,
    );
    let target = dna(
        "Independent inverse target",
        "AGTGGAGATGCTACGTAGGACTTCAGACCTGTA",
        Topology::Circular,
    );
    s.import(vec![source.clone(), target.clone()]).await;
    let options = InversePcrOptions {
        edit_start: 3,
        edit_end: 5,
        inserted: "GGA".into(),
        forward_length: 8,
        reverse_length: 8,
        reverse_tail_bases: 0,
    };
    let response = s
        .call(
            "design_inverse_pcr",
            json!({"document":reference(&source),"target":reference(&target),"origin":null,"origin_workspace":null,"options":options}),
        )
        .await;
    assert_eq!(response["computed_status"], "ConsistentUnderAssumptions");
    assert_eq!(response["forward"]["tail"], "GGA");
    assert_eq!(response["forward"]["annealing"]["bases"], "GATGCTAC");
    assert_eq!(response["reverse"]["annealing"]["bases"], "ACTTACAG");
    let artifact = s.artifact(&response).await;
    assert_eq!(artifact["document"], reference(&source));
    assert_eq!(artifact["target"], reference(&target));
    let steps = artifact["inverse_pcr"]["plan"]["steps"].as_array().unwrap();
    assert_eq!(steps.len(), 3);
    assert!(steps[0].get("Pcr").is_some());
    assert_eq!(steps[1]["TreatEnds"]["treatment"], "Phosphorylate");
    assert_eq!(steps[2]["Ligate"]["circular"], true);
    assert_eq!(steps[2]["Ligate"]["allow_nicked"], false);
    let report = s
        .call(
            "verify_plan",
            json!({"plan_revision":response["plan_revision"]}),
        )
        .await;
    assert_eq!(report["status"], "ConsistentUnderAssumptions");
    assert_eq!(report["steps"].as_array().unwrap().len(), 3);
    s.pdf(&response, "ConsistentUnderAssumptions", &target)
        .await;
    for (document, target_ref) in [
        (stale(&source), reference(&target)),
        (reference(&source), stale(&target)),
    ] {
        s.reject(
            "design_inverse_pcr",
            json!({"document":document,"target":target_ref,"origin":null,"origin_workspace":null,"options":options}),
        )
        .await;
    }
    let wrong = s
        .call(
            "design_inverse_pcr",
            json!({"document":reference(&source),"target":reference(&source),"origin":null,"origin_workspace":null,"options":options}),
        )
        .await;
    assert_eq!(wrong["computed_status"], "Invalid");
    assert_eq!(
        s.call("list_documents", json!({})).await["documents"]
            .as_array()
            .unwrap()
            .len(),
        2
    );
    s.finish().await;
}

#[tokio::test]
async fn gel_json_svg_and_diagnostic_digests_preserve_unknowns_mass_and_source_bindings() {
    let s = Session::new().await;
    // EcoRV boundaries7/17 on24bp ->7,10,7 linear,10/14 circular.
    let linear = dna(
        "Linear reference",
        "AAAAGATATCCCCCGATATCGGGG",
        Topology::Linear,
    );
    let circle = dna(
        "Circular reference",
        "AAAAGATATCCCCCGATATCGGGG",
        Topology::Circular,
    );
    let alternative = dna(
        "Single-cut alternative",
        "AAAAGATATCCCCCGTATACGGGG",
        Topology::Circular,
    );
    s.import(vec![linear.clone(), circle.clone(), alternative.clone()])
        .await;
    let lane = |m: &Molecule, preparation: Preparation, mass: Option<f64>| json!({"label":m.name,"document":reference(m),"preparation":preparation,"mass_ng":mass});
    let query = json!({"lanes":[lane(&linear,Preparation::Digest{enzymes:vec!["EcoRV".into()]},Some(240.)),lane(&circle,Preparation::Uncut,None)],"calibration":calibration(),"image_format":"Svg"});
    let response = s.call("predict_gel", query.clone()).await;
    assert_eq!(response["lane_count"], 2);
    assert_eq!(response["lanes"][0]["fragment_count"], 3);
    assert_eq!(response["lanes"][0]["resolved_band_count"], 2);
    assert_eq!(response["lanes"][1]["resolved_band_count"], 0);
    assert_eq!(response["lanes"][1]["unpositioned_fragment_count"], 1);
    let artifact = s.artifact(&response).await;
    assert_eq!(
        artifact["documents"],
        json!([reference(&linear), reference(&circle)])
    );
    let fragments = artifact["gel"]["lanes"][0]["fragments"].as_array().unwrap();
    assert_eq!(
        fragments
            .iter()
            .map(|f| f["top_nt"].as_u64().unwrap())
            .collect::<Vec<_>>(),
        vec![7, 10, 7]
    );
    assert_eq!(
        fragments
            .iter()
            .map(|f| f["amount_ng"].as_f64().unwrap())
            .sum::<f64>(),
        240.
    );
    let svg = s
        .text(&response["image_artifact_id"], "image/svg+xml", "utf-8")
        .await;
    assert_eq!(
        blake3::hash(svg.as_bytes()).to_hex().as_str(),
        response["image_artifact_id"].as_str().unwrap()
    );
    assert!(svg.contains("collagene-gel-evidence"));
    assert!(svg.contains(&circle.revision()));
    assert!(svg.contains("1 unplaced"));
    let mut png_query = query.clone();
    png_query["image_format"] = json!({"Png":{"width":920}});
    let png_response = s.call("predict_gel", png_query.clone()).await;
    let png_bytes = STANDARD
        .decode(
            s.text(&png_response["image_artifact_id"], "image/png", "base64")
                .await,
        )
        .unwrap();
    assert_eq!(&png_bytes[..8], b"\x89PNG\r\n\x1a\n");
    assert_eq!(
        u32::from_be_bytes(png_bytes[16..20].try_into().unwrap()),
        920
    );
    assert_eq!(
        blake3::hash(&png_bytes).to_hex().as_str(),
        png_response["image_artifact_id"].as_str().unwrap()
    );
    let png_evidence = s.artifact(&png_response).await;
    assert_eq!(png_evidence["gel"], artifact["gel"]);
    assert_eq!(
        png_evidence["image_context_revision"],
        png_response["image_context_revision"]
    );
    png_query["image_format"] = json!({"Png":{"width":4097}});
    s.reject("predict_gel", png_query).await;

    for i in 0..2 {
        let mut bad = query.clone();
        bad["lanes"][i]["document"]["revision"] = json!("stale");
        s.reject("predict_gel", bad).await;
    }
    let options = DiagnosticOptions {
        enzymes: vec!["BamHI".into(), "EcoRV".into()],
        calibration: calibration(),
        ..Default::default()
    };
    let query = json!({"documents":[reference(&circle),reference(&alternative)],"options":options});
    let diagnostic = s.call("design_diagnostic_digest", query.clone()).await;
    assert_eq!(diagnostic["exhaustive_search"], true);
    assert_eq!(diagnostic["evaluated_digests"], 3);
    assert_eq!(diagnostic["candidates"][0]["enzymes"], json!(["EcoRV"]));
    assert_eq!(
        diagnostic["candidates"][0]["distinguished_pairs"],
        json!([[0, 1]])
    );
    let bam = diagnostic["candidates"]
        .as_array()
        .unwrap()
        .iter()
        .find(|c| c["enzymes"] == json!(["BamHI"]))
        .unwrap();
    assert_eq!(bam["unknown_pairs"], json!([[0, 1]]));
    let evidence = s.artifact(&diagnostic).await;
    assert_eq!(
        evidence["diagnostic_digests"]["source_revisions"],
        json!([circle.revision(), alternative.revision()])
    );
    let mut limited = query.clone();
    limited["options"]["max_evaluations"] = json!(1);
    let limited = s.call("design_diagnostic_digest", limited).await;
    assert_eq!(limited["exhaustive_search"], false);
    assert_eq!(limited["evaluated_digests"], 1);
    for i in 0..2 {
        let mut bad = query.clone();
        bad["documents"][i]["revision"] = json!("stale");
        s.reject("design_diagnostic_digest", bad).await;
    }
    s.finish().await;
}

#[tokio::test]
async fn pinned_gel_requires_exact_workspace_ownership_and_replays_historical_sources() {
    use collagene_core::{
        cloning::restriction::{
            catalog as catalog_digest,
            inspection::InspectionOptions,
            setups::{CleavageMode, Recipe},
        },
        documents::{
            history::DocumentHistory,
            workspace::{self, Workspace},
        },
        material::lineage::EventBinding,
    };
    let s = Session::new().await;
    let source = dna("historical AluI source", "TTTAGCTAAA", Topology::Linear);
    // Both owners deliberately share the same molecular revision. Sequence/hash
    // equality therefore cannot substitute for exact event ownership.
    let mut frame = Workspace::from_documents(
        "retained gel design",
        vec![source.clone(), source.clone()],
        0,
    )
    .unwrap();
    let old = &frame.documents[0];
    let binding = EventBinding {
        document_id: old.id.clone(),
        event_id: old.head_event.clone().unwrap(),
        revision: old.head.clone(),
    };
    let other = &frame.documents[1];
    let other_owner = other.id.clone();
    let other_event = other.head_event.clone().unwrap();
    let mut histories = DocumentHistory::restore_all(&frame).unwrap();
    let edit_end = histories[0].molecule.sequence.len();
    histories[0]
        .commit(
            vec![
                Edit::ReplaceSequence {
                    start: 0,
                    end: edit_end,
                    replacement: "TTTACCTAAA".into(),
                },
                Edit::Rename {
                    name: "later edited source".into(),
                },
            ],
            "Later source edit",
        )
        .unwrap();
    let later = histories[0].molecule.clone();
    frame.documents.clear();
    frame.history.clear();
    for history in histories {
        history.contribute(&mut frame).unwrap();
    }
    frame.validate().unwrap();
    assert_ne!(frame.documents[0].head, binding.revision);
    let bytes = workspace::save_workspace(&frame).unwrap();
    let imported = s
        .call(
            "import_workspace_archive",
            json!({"data_base64":STANDARD.encode(bytes)}),
        )
        .await;
    let workspace_ref =
        json!({"id":imported["workspace"]["id"],"revision":imported["workspace"]["revision"]});
    let mut alias = source.clone();
    alias.name = "same DNA with invented provenance".into();
    alias
        .metadata
        .insert("source_event".into(), binding.event_id.clone());
    alias
        .metadata
        .insert("source_document".into(), binding.document_id.clone());
    // Explicit immutable DocumentRefs remain separate from native owner/event IDs.
    s.import(vec![source.clone(), later.clone(), alias.clone()])
        .await;
    let recipe = Recipe {
        source: binding.clone(),
        origin: collagene_core::material::derivation::SourceOrigin::OriginalReference,
        selection: catalog_digest::select(
            &catalog_digest::builtin_catalog().unwrap(),
            &["alui".into()],
        )
        .unwrap(),
        options: InspectionOptions::default(),
        cleavage_mode: CleavageMode::AssumedCompleteMethylationUnassessed,
        conditions: None,
    };
    let lane = json!({"label":"Historical AluI digest","document":reference(&source),"workspace":workspace_ref,"preparation":Preparation::CatalogDigest {recipe:recipe.clone()},"mass_ng":100.});
    let query = json!({"lanes":[lane],"calibration":calibration(),"image_format":"Svg"});
    let response = s.call("predict_gel", query.clone()).await;
    assert_eq!(response["lane_count"], 1);
    assert_eq!(response["lanes"][0]["fragment_count"], 2);
    assert_eq!(response["lanes"][0]["complete_under_model"], true);
    let artifact = s.artifact(&response).await;
    assert_eq!(artifact["lanes"], query["lanes"]);
    assert_eq!(artifact["documents"], json!([reference(&source)]));
    let fragments = artifact["gel"]["lanes"][0]["fragments"].as_array().unwrap();
    assert_eq!(
        fragments
            .iter()
            .map(|f| (
                f["top_nt"].as_u64().unwrap(),
                f["bottom_nt"].as_u64().unwrap()
            ))
            .collect::<Vec<_>>(),
        vec![(5, 5), (5, 5)]
    );
    assert_eq!(
        fragments
            .iter()
            .map(|f| f["amount_ng"].as_f64().unwrap())
            .sum::<f64>(),
        100.
    );
    let mut missing = query.clone();
    missing["lanes"][0]
        .as_object_mut()
        .unwrap()
        .remove("workspace");
    s.reject("predict_gel", missing).await;
    let mut wrong_workspace_revision = query.clone();
    wrong_workspace_revision["lanes"][0]["workspace"]["revision"] = json!("0".repeat(64));
    s.reject("predict_gel", wrong_workspace_revision).await;
    let mut wrong_owner = query.clone();
    wrong_owner["lanes"][0]["preparation"]["CatalogDigest"]["recipe"]["source"]["document_id"] =
        json!(other_owner);
    s.reject("predict_gel", wrong_owner).await;
    let mut wrong_event = query.clone();
    wrong_event["lanes"][0]["preparation"]["CatalogDigest"]["recipe"]["source"]["event_id"] =
        json!(other_event);
    s.reject("predict_gel", wrong_event).await;
    let mut missing_event = query.clone();
    missing_event["lanes"][0]["preparation"]["CatalogDigest"]["recipe"]["source"]["event_id"] =
        json!(format!("event-{}", "a".repeat(64)));
    s.reject("predict_gel", missing_event).await;
    let mut wrong_revision = query.clone();
    wrong_revision["lanes"][0]["preparation"]["CatalogDigest"]["recipe"]["source"]["revision"] =
        json!(later.revision());
    s.reject("predict_gel", wrong_revision).await;
    for substituted in [&later, &alias] {
        let mut wrong_document = query.clone();
        wrong_document["lanes"][0]["document"] = reference(substituted);
        s.reject("predict_gel", wrong_document).await;
    }
    let unrelated =
        Workspace::from_documents("unrelated native workspace", vec![alias.clone()], 0).unwrap();
    let imported = s
        .call(
            "import_workspace_archive",
            json!({"data_base64":STANDARD.encode(workspace::save_workspace(&unrelated).unwrap())}),
        )
        .await;
    let mut wrong_workspace = query.clone();
    wrong_workspace["lanes"][0]["workspace"] =
        json!({"id":imported["workspace"]["id"],"revision":imported["workspace"]["revision"]});
    s.reject("predict_gel", wrong_workspace).await;
    let mut too_many = query.clone();
    too_many["lanes"] = json!(vec![query["lanes"][0].clone(); 33]);
    s.reject("predict_gel", too_many).await;
    // Valid ownership for the second document is acceptable even though the
    // molecular hash is shared: only a mismatched owner/event pair is rejected.
    let second = &frame.documents[1];
    let mut independently_owned = query.clone();
    independently_owned["lanes"][0]["preparation"]["CatalogDigest"]["recipe"]["source"] =
        json!({"document_id":second.id,"event_id":second.head_event,"revision":second.head});
    assert_eq!(
        s.call("predict_gel", independently_owned).await["lanes"][0]["fragment_count"],
        2
    );
    s.finish().await;
}

#[tokio::test]
async fn exact_pcr_long_concrete_late_mismatch_finishes_without_inventing_products() {
    let session = Session::new().await;
    let template = dna(
        "Long concrete template",
        &"A".repeat(200_000),
        Topology::Circular,
    );
    session.import(vec![template.clone()]).await;
    let forward = Primer::new("late mismatch", &format!("{}C", "A".repeat(39_999)), "").unwrap();
    let reverse = Primer::new("absent", &"C".repeat(40_000), "").unwrap();
    let report = session
        .call(
            "simulate_pcr",
            json!({
                "document": reference(&template), "forward": forward, "reverse": reverse
            }),
        )
        .await;
    assert_eq!(report["product_count"], 0);
    assert_eq!(report["products"], json!([]));
    assert_eq!(
        session.call("list_documents", json!({})).await["documents"]
            .as_array()
            .unwrap()
            .len(),
        1
    );
    session.finish().await;
}
