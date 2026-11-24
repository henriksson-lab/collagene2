//! Actual SDK calls with independently authored variant evidence and DNA.
use base64::{Engine as _, engine::general_purpose::STANDARD};
use collagene_mcp::CollageneServer;
use rmcp::{RoleClient, ServiceExt, model::CallToolRequestParams, service::RunningService};
use serde_json::{Value, json};
type Client = RunningService<RoleClient, ()>;
const VCF:&[u8]=b"##fileformat=VCFv4.3\n##contig=<ID=chrExplicit,length=8>\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Declared genotype\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample_A\tsample_B\nchrExplicit\t2\tdeclared\tC\tTA,G\t01.50\t.\t.\tGT\t0|2\t./.\nchrExplicit\t1\tsymbolic\tA\t<DEL>\t.\t.\t.\tGT\t1/1\t.\n";
async fn call(c: &Client, name: &str, args: Value) -> Value {
    let response = c
        .call_tool(
            CallToolRequestParams::new(name.to_owned())
                .with_arguments(args.as_object().unwrap().clone()),
        )
        .await
        .unwrap();
    assert_ne!(
        response.is_error,
        Some(true),
        "{name}: {:?}",
        response.content
    );
    response.structured_content.unwrap()
}
async fn reject(c: &Client, name: &str, args: Value) {
    let response = c
        .call_tool(
            CallToolRequestParams::new(name.to_owned())
                .with_arguments(args.as_object().unwrap().clone()),
        )
        .await;
    assert!(
        response.is_err() || response.unwrap().is_error == Some(true),
        "{name} unexpectedly accepted invalid input"
    );
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
async fn artifact(c: &Client, id: &Value) -> (Vec<u8>, String) {
    let mut text = String::new();
    let mut start = 0;
    loop {
        let page = call(
            c,
            "read_artifact",
            json!({"artifact_id":id,"start":start,"max_bytes":1024}),
        )
        .await;
        text.push_str(page["text"].as_str().unwrap());
        if page["done"] == true {
            let encoding = page["encoding"].as_str().unwrap().to_owned();
            return (
                if encoding == "base64" {
                    STANDARD.decode(text).unwrap()
                } else {
                    text.into_bytes()
                },
                encoding,
            );
        }
        start = page["next_start"].as_u64().unwrap();
    }
}
async fn document(c: &Client, dna: &str) -> Value {
    let result = call(
        c,
        "import_document",
        json!({"text":format!(">explicit_actual_DNA\n{dna}\n")}),
    )
    .await;
    json!({"id":result["documents"][0]["id"],"revision":result["documents"][0]["revision"]})
}

#[tokio::test]
async fn private_variant_sources_fresh_reference_checks_and_atomic_reference_only_previews() {
    let (c, s) = setup().await;
    let tools = c.list_tools(Default::default()).await.unwrap().tools;
    assert_eq!(tools.len(), collagene_mcp::CURRENT_TOOL_COUNT);
    for name in [
        "import_variants",
        "inspect_variants",
        "check_variant_reference",
        "preview_variant_allele",
        "export_variants",
        "dispose_variants",
        "replay_variant_preview",
    ] {
        let t = tools.iter().find(|t| t.name == name).unwrap();
        assert_eq!(t.input_schema["additionalProperties"], false);
    }
    let reference = document(&c, "ACGTACGT").await;
    let wrong = document(&c, "AAAAAAAA").await;
    let imported = call(
        &c,
        "import_variants",
        json!({"data_base64":STANDARD.encode(VCF),"format":"vcf"}),
    )
    .await;
    assert_eq!(imported["record_count"], 2);
    assert_eq!(imported["sample_column_count"], 2);
    assert_eq!(
        imported["source_blake3"],
        blake3::hash(VCF).to_hex().to_string()
    );
    let collection = imported["collection"].clone();
    assert_eq!(collection["id"], collection["revision"]);
    let page = call(
        &c,
        "inspect_variants",
        json!({"collection":collection,"start":0,"limit":1}),
    )
    .await;
    assert_eq!(page["done"], false);
    assert_eq!(page["records"][0]["position_1based"], 2);
    let full: Value = serde_json::from_slice(&artifact(&c, &page["artifact_id"]).await.0).unwrap();
    assert!(
        full["records"][0]["vcf"]
            .as_str()
            .unwrap()
            .contains("0|2\t./.")
    );
    let args = json!({"collection":collection,"record_index":0,"reference":reference,"reference_name":"chrExplicit","alternate_index":1});
    let checked = call(&c, "check_variant_reference", args.clone()).await;
    assert_eq!(checked["outcome"], "ReferenceAlleleMatches");
    let full: Value =
        serde_json::from_slice(&artifact(&c, &checked["artifact_id"]).await.0).unwrap();
    assert_eq!(
        full["report"]["options"]["expected_reference_revision"],
        reference["revision"]
    );
    assert_eq!(full["report"]["source_blake3"], collection["id"]);
    let preview = call(&c, "preview_variant_allele", args.clone()).await;
    assert_eq!(preview["reference_only"], true);
    let product =
        json!({"id":preview["document"]["id"],"revision":preview["document"]["revision"]});
    let region = call(
        &c,
        "read_sequence_region",
        json!({"document":product,"start":0,"end":9}),
    )
    .await;
    assert_eq!(region["sequence"], "ATAGTACGT");
    let replay_args = json!({"source_id":collection["id"],"reference":reference,"design":product});
    let replay = call(&c, "replay_variant_preview", replay_args.clone()).await;
    assert_eq!(replay["literal_dna_and_topology_match"], true);
    assert_eq!(replay["claimed_revision"], product["revision"]);
    assert_eq!(replay["fresh_revision"], product["revision"]);
    let full: Value =
        serde_json::from_slice(&artifact(&c, &replay["artifact_id"]).await.0).unwrap();
    assert_eq!(full["source_id"], collection["id"]);
    assert_eq!(full["reference"], reference);
    assert_eq!(full["design"], product);
    let mut invalid = replay_args.clone();
    invalid["reference"] = wrong.clone();
    reject(&c, "replay_variant_preview", invalid).await;
    let mut invalid = replay_args.clone();
    invalid["design"] = reference.clone();
    reject(&c, "replay_variant_preview", invalid).await;
    let unchanged = call(
        &c,
        "read_sequence_region",
        json!({"document":reference,"start":0,"end":8}),
    )
    .await;
    assert_eq!(unchanged["sequence"], "ACGTACGT");
    let mut bad = args.clone();
    bad["reference"] = wrong;
    let mismatch = call(&c, "check_variant_reference", bad.clone()).await;
    assert_eq!(mismatch["outcome"], "ReferenceAlleleMismatch");
    reject(&c, "preview_variant_allele", bad).await;
    let mut bad = args.clone();
    bad["reference_name"] = json!("same_name_is_not_mapping");
    reject(&c, "check_variant_reference", bad).await;
    let mut bad = args.clone();
    bad["reference"]["revision"] = json!("0".repeat(64));
    reject(&c, "preview_variant_allele", bad).await;
    let mut bad = args.clone();
    bad["alternate_index"] = json!(0);
    reject(&c, "check_variant_reference", bad).await;
    let mut symbolic = args.clone();
    symbolic["record_index"] = json!(1);
    reject(&c, "preview_variant_allele", symbolic).await;
    reject(
        &c,
        "inspect_variants",
        json!({"collection":collection,"limit":101}),
    )
    .await;
    reject(
        &c,
        "import_variants",
        json!({"data_base64":STANDARD.encode(b"not a file"),"format":"vcf"}),
    )
    .await;
    reject(
        &c,
        "export_variants",
        json!({"collection":collection,"format":"bcf","max_bytes":10}),
    )
    .await;
    let original = call(
        &c,
        "export_variants",
        json!({"collection":collection,"format":"original"}),
    )
    .await;
    assert_eq!(artifact(&c, &original["artifact_id"]).await.0, VCF);
    let bcf = call(
        &c,
        "export_variants",
        json!({"collection":collection,"format":"bcf"}),
    )
    .await;
    assert!(!bcf["notices"].as_array().unwrap().is_empty());
    let binary = artifact(&c, &bcf["artifact_id"]).await.0;
    let other = call(
        &c,
        "import_variants",
        json!({"data_base64":STANDARD.encode(binary),"format":"bcf"}),
    )
    .await;
    assert_eq!(other["record_count"], 2);
    call(&c, "dispose_variants", json!({"collection":collection})).await;
    reject(&c, "inspect_variants", json!({"collection":collection})).await;
    let region = call(
        &c,
        "read_sequence_region",
        json!({"document":product,"start":0,"end":9}),
    )
    .await;
    assert_eq!(region["sequence"], "ATAGTACGT");
    drop(c);
    s.abort();
}
