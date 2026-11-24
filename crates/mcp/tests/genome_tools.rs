//! Actual SDK file-backed reference lifecycle; no arbitrary paths or giant Molecule.
use base64::{Engine, engine::general_purpose::STANDARD};
use collagene_core::{
    Topology,
    crispr::{GuideOptions, GuideQuery},
    genome::*,
    genome_crispr::*,
};
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
    let client = ().serve(b).await.unwrap();
    (client, server)
}

const SPACER: &str = "ACGTACGTACGTACGTACGT";
fn fixture() -> (Vec<u8>, BuildInfo) {
    let mut dna = vec![b'A'; BLOCK_BASES + 100];
    let at = BLOCK_BASES - 10;
    dna[at..at + 20].copy_from_slice(SPACER.as_bytes());
    dna[at + 20..at + 23].copy_from_slice(b"AGG");
    let mut bytes = b">chr\n".to_vec();
    bytes.extend(&dna);
    bytes.push(b'\n');
    let mut build = BuildInfo::default();
    build.name = "Declared SDK fixture build".into();
    build.expected_components = vec![ExpectedComponent {
        id: "chr".into(),
        source_id: "chr".into(),
        accession: None,
        length: dna.len() as u64,
        sequence_hash: None,
        topology: Topology::Linear,
    }];
    (bytes, build)
}
async fn append(c: &Client, id: &Value, bytes: &[u8]) {
    let mut offset = 0;
    for chunk in bytes.chunks(7001) {
        let p = call(
            c,
            "append_genome_import",
            json!({"import_id":id,"offset":offset,"data_base64":STANDARD.encode(chunk)}),
        )
        .await;
        offset = p["received"].as_u64().unwrap();
    }
    assert_eq!(offset, bytes.len() as u64);
}
async fn finish(c: &Client, id: &Value) -> Value {
    for _ in 0..100 {
        let p = call(
            c,
            "finish_genome_import",
            json!({"import_id":id,"max_blocks":1}),
        )
        .await;
        if p["done"] == true {
            return p["dataset"].clone();
        }
        assert!(p["dataset"].is_null());
    }
    panic!("bounded test import did not finish")
}
#[tokio::test]
async fn native_sdk_chunked_reference_global_seam_guides_private_search_and_export_reimport() {
    let (c, s) = setup().await;
    assert_eq!(
        c.list_tools(Default::default()).await.unwrap().tools.len(),
        97
    );
    let (bytes, build) = fixture();
    let before = call(&c, "list_documents", json!({})).await;
    let begin=call(&c,"begin_genome_import",json!({"format":"Fasta","expected_bytes":bytes.len(),"expected_hash":blake3::hash(&bytes).to_hex().to_string(),"build":build})).await;
    let id = &begin["import_id"];
    reject(
        &c,
        "append_genome_import",
        json!({"import_id":id,"offset":1,"data_base64":STANDARD.encode(&bytes[..10])}),
    )
    .await;
    append(&c, id, &bytes).await;
    let dataset = finish(&c, id).await;
    assert_eq!(dataset["semantic_verified_at_import"], true);
    assert_eq!(dataset["expected_components_complete"], true);
    let genome = dataset["genome"].clone();
    let inspect = call(
        &c,
        "inspect_genome",
        json!({"genome":genome,"max_contigs":1}),
    )
    .await;
    assert_eq!(inspect["contigs"][0]["block_count"], 2);
    let start = BLOCK_BASES as u64 - 10;
    let read_options = GenomeIntervalReadOptions {
        source_role: "reference".into(),
        interval: GenomicInterval {
            contig_id: "chr".into(),
            start,
            end: start + 23,
        },
    };
    let read = call(
        &c,
        "read_genome_interval",
        json!({"genome":genome,"options":read_options}),
    )
    .await;
    assert_eq!(read["read"]["sequence"], format!("{SPACER}AGG"));
    let guide_options = GenomeGuideOptions {
        source_role: "reference".into(),
        interval: GenomicInterval {
            contig_id: "chr".into(),
            start,
            end: start + 1,
        },
        guide: GuideOptions::default(),
    };
    let guides = call(
        &c,
        "design_genome_interval_guides",
        json!({"genome":genome,"options":guide_options,"max_candidates":1}),
    )
    .await;
    let candidate = &guides["candidates"][0];
    assert!(
        candidate["guide_positions"]
            .as_array()
            .unwrap()
            .contains(&json!(BLOCK_BASES as u64))
    );
    let query = GuideQuery {
        id: "query".into(),
        sequence: SPACER.into(),
    };
    let mut options = GenomeSearchOptions::default();
    options.max_mismatches = 0;
    options.intended = Some(GenomeIntended {
        reference_role: "reference".into(),
        manifest_hash: genome["manifest_hash"].as_str().unwrap().into(),
        content_revision: genome["content_revision"].as_str().unwrap().into(),
        contig_id: "chr".into(),
        guide_positions: serde_json::from_value(candidate["guide_positions"].clone()).unwrap(),
        pam_positions: serde_json::from_value(candidate["pam_positions"].clone()).unwrap(),
        reverse: false,
    });
    let mut search = call(
        &c,
        "start_genome_off_target_search",
        json!({"genome":genome,"query":query,"options":options}),
    )
    .await;
    assert_eq!(search["progress"]["source_integrity_verified"], false);
    let job = search["job_id"].clone();
    while search["progress"]["finished"] != true {
        search = call(
            &c,
            "advance_genome_off_target_search",
            json!({"job_id":job,"max_blocks":1}),
        )
        .await;
    }
    assert_eq!(search["progress"]["source_integrity_verified"], true);
    assert_eq!(search["progress"]["complete"], true);
    assert_eq!(
        search["progress"]["bases_verified"],
        (BLOCK_BASES + 100) as u64
    );
    let hits = call(
        &c,
        "read_genome_off_target_hits",
        json!({"job_id":job,"limit":1}),
    )
    .await;
    assert_eq!(hits["hits"].as_array().unwrap().len(), 1);
    assert_eq!(hits["hits"][0]["is_intended"], true);
    reject(
        &c,
        "advance_genome_off_target_search",
        json!({"job_id":job,"checkpoint":{"finished":true,"source_integrity_verified":true}}),
    )
    .await;
    call(&c, "cancel_genome_job", json!({"id":job})).await;
    reject(&c, "read_genome_off_target_hits", json!({"job_id":job})).await;
    let mut native = vec![];
    let mut offset = 0;
    loop {
        let p = call(
            &c,
            "read_genome_export",
            json!({"genome":genome,"offset":offset,"max_bytes":4096}),
        )
        .await;
        native.extend(STANDARD.decode(p["data_base64"].as_str().unwrap()).unwrap());
        if p["done"] == true {
            break;
        }
        offset = p["next_offset"].as_u64().unwrap();
    }
    assert_eq!(native.len() as u64, dataset["file_length"]);
    assert!(native.starts_with(&HEADER));
    let begin=call(&c,"begin_genome_import",json!({"format":"Cgref","expected_bytes":native.len(),"expected_hash":blake3::hash(&native).to_hex().to_string()})).await;
    append(&c, &begin["import_id"], &native).await;
    let copy = finish(&c, &begin["import_id"]).await;
    assert_eq!(copy["genome"]["manifest_hash"], genome["manifest_hash"]);
    assert_eq!(
        copy["genome"]["content_revision"],
        genome["content_revision"]
    );
    let mut stale = genome.clone();
    stale["manifest_hash"] = json!("0".repeat(64));
    reject(&c, "inspect_genome", json!({"genome":stale})).await;
    reject(
        &c,
        "read_genome_export",
        json!({"genome":genome,"path":"/etc/passwd"}),
    )
    .await;
    reject(&c,"begin_genome_import",json!({"format":"Cgref","expected_bytes":native.len(),"filesystem_path":"/tmp/arbitrary.cgref"})).await;
    let mut corrupt = native;
    corrupt[8] ^= 1;
    let bad=call(&c,"begin_genome_import",json!({"format":"Cgref","expected_bytes":corrupt.len(),"expected_hash":blake3::hash(&corrupt).to_hex().to_string()})).await;
    append(&c, &bad["import_id"], &corrupt).await;
    reject(
        &c,
        "finish_genome_import",
        json!({"import_id":bad["import_id"],"max_blocks":64}),
    )
    .await;
    assert_eq!(call(&c, "list_documents", json!({})).await, before);
    call(
        &c,
        "cancel_genome_job",
        json!({"id":genome["id"],"genome":genome}),
    )
    .await;
    reject(&c, "inspect_genome", json!({"genome":dataset["genome"]})).await;
    call(
        &c,
        "cancel_genome_job",
        json!({"id":copy["genome"]["id"],"genome":copy["genome"]}),
    )
    .await;
    assert!(
        call(&c, "inspect_genome", json!({})).await["datasets"]
            .as_array()
            .unwrap()
            .is_empty()
    );
    c.cancel().await.unwrap();
    s.await.unwrap();
}
#[tokio::test]
async fn import_truncation_hash_failures_cancellation_and_build_labels_do_not_grant_coverage() {
    let (c, s) = setup().await;
    let bytes = format!(">label_only\n{SPACER}AGG\n").into_bytes();
    let mut build = BuildInfo::default();
    build.name = "A label claiming complete genome".into();
    let begin = call(
        &c,
        "begin_genome_import",
        json!({"format":"Fasta","expected_bytes":bytes.len(),"build":build}),
    )
    .await;
    let id = &begin["import_id"];
    let n = bytes.len() - 2;
    call(
        &c,
        "append_genome_import",
        json!({"import_id":id,"offset":0,"data_base64":STANDARD.encode(&bytes[..n])}),
    )
    .await;
    reject(&c, "finish_genome_import", json!({"import_id":id})).await;
    call(&c, "cancel_genome_job", json!({"id":id})).await;
    reject(
        &c,
        "append_genome_import",
        json!({"import_id":id,"offset":n,"data_base64":STANDARD.encode(&bytes[n..])}),
    )
    .await;
    let bad = call(
        &c,
        "begin_genome_import",
        json!({"format":"Fasta","expected_bytes":bytes.len(),"expected_hash":"0".repeat(64)}),
    )
    .await;
    append(&c, &bad["import_id"], &bytes).await;
    reject(
        &c,
        "finish_genome_import",
        json!({"import_id":bad["import_id"]}),
    )
    .await;
    assert!(
        call(&c, "inspect_genome", json!({})).await["datasets"]
            .as_array()
            .unwrap()
            .is_empty()
    );
    let begin = call(
        &c,
        "begin_genome_import",
        json!({"format":"Fasta","expected_bytes":bytes.len(),"build":build}),
    )
    .await;
    append(&c, &begin["import_id"], &bytes).await;
    let dataset = finish(&c, &begin["import_id"]).await;
    assert_eq!(dataset["expected_components_complete"], false);
    let mut options = GenomeSearchOptions::default();
    options.max_mismatches = 0;
    let q = GuideQuery {
        id: "label".into(),
        sequence: SPACER.into(),
    };
    let search = call(
        &c,
        "start_genome_off_target_search",
        json!({"genome":dataset["genome"],"query":q,"options":options}),
    )
    .await;
    let search = call(
        &c,
        "advance_genome_off_target_search",
        json!({"job_id":search["job_id"],"max_blocks":1}),
    )
    .await;
    assert_eq!(search["progress"]["exhaustive_for_manifest"], true);
    assert_eq!(search["progress"]["source_integrity_verified"], true);
    assert_eq!(search["progress"]["expected_components_complete"], false);
    assert_eq!(search["progress"]["complete"], false);
    c.cancel().await.unwrap();
    s.await.unwrap();
}
