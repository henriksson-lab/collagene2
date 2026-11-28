//! Actual SDK routes for source-bound CRISPR and independent donor checking.
use collagene_core::sequence::{
    crispr::*,
    template::{self as edit_template, EditTemplateCheck, EditTemplateOptions, EditTemplateReport},
};
use collagene_mcp::server::CollageneServer;
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
async fn artifact(c: &Client, id: &Value) -> Value {
    let mut text = String::new();
    let mut at = 0;
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
fn doc(imported: &Value, name: &str) -> Value {
    let d = imported["documents"]
        .as_array()
        .unwrap()
        .iter()
        .find(|d| {
            d["name"]
                .as_str()
                .is_some_and(|s| s.split_whitespace().next() == Some(name))
        })
        .unwrap();
    json!({"id":d["id"],"revision":d["revision"]})
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
fn guide_text() -> String {
    format!(">reference\n{SPACER}AGGTTT{SPACER}AGG\n")
}
const EDIT_TEXT: &str = ">source\nAAAACCCCGGGGTTTT\n>target\nAAAACCTCGGGGTTTT\n";
fn edit_options() -> EditTemplateOptions {
    EditTemplateOptions {
        edit_start: 6,
        delete_length: 1,
        replacement: "T".into(),
        left_arm: 3,
        right_arm: 3,
        reverse_donor: false,
        max_cut_distance: 100,
    }
}
#[tokio::test]
async fn guide_and_selected_reference_pages_retain_exact_intended_roles_and_all_same_sequence_hits()
{
    let (c, s) = setup().await;
    assert_eq!(
        c.list_tools(Default::default()).await.unwrap().tools.len(),
        collagene_mcp::CURRENT_TOOL_COUNT
    );
    let imported = call(&c, "import_document", json!({"text":guide_text()})).await;
    let document = doc(&imported, "reference");
    let options = GuideOptions::default();
    let before = call(&c, "list_documents", json!({})).await;
    let first = call(
        &c,
        "design_crispr_guides",
        json!({"document":document,"options":options,"max_candidates":1}),
    )
    .await;
    assert_eq!(first["candidates"].as_array().unwrap().len(), 1);
    assert_eq!(first["available_profiles"].as_array().unwrap().len(), 5);
    let guides: GuideReport =
        serde_json::from_value(artifact(&c, &first["artifact_id"]).await).unwrap();
    let molecule = collagene_core::formats::sequence::import_text(&guide_text())
        .unwrap()
        .molecules
        .remove(0);
    guides.validate_bindings_for(&molecule, &options).unwrap();
    let guide = guides
        .candidates
        .iter()
        .find(|g| !g.reverse && g.guide_positions[0] == 0 && g.guide_sequence == SPACER)
        .unwrap();
    let page = call(
        &c,
        "design_crispr_guides",
        json!({"document":document,"options":options,"candidate_start":1,"max_candidates":1}),
    )
    .await;
    assert_eq!(page["artifact_id"], first["artifact_id"]);
    let references = json!([{"id":"intended","document":document,"build":"declared-build","accession":"ref.1"},{"id":"same-sequence-other-role","document":document,"build":null,"accession":null}]);
    let query = GuideQuery {
        id: guide.id.clone(),
        sequence: guide.guide_sequence.clone(),
    };
    let mut off = OffTargetOptions::default();
    off.max_mismatches = 0;
    off.intended = Some(IntendedTarget {
        reference_id: "intended".into(),
        reference_revision: molecule.revision(),
        guide_positions: guide.guide_positions.clone(),
        pam_positions: guide.pam_positions.clone(),
        reverse: guide.reverse,
    });
    let result = call(
        &c,
        "search_crispr_off_targets",
        json!({"query":query,"references":references,"options":off,"max_hits":1}),
    )
    .await;
    assert_eq!(result["hit_count"], 4);
    assert_eq!(result["hits"].as_array().unwrap().len(), 1);
    let report: OffTargetReport =
        serde_json::from_value(artifact(&c, &result["artifact_id"]).await).unwrap();
    let inputs = vec![
        ReferenceInput {
            id: "intended".into(),
            molecule: molecule.clone(),
            build: Some("declared-build".into()),
            accession: Some("ref.1".into()),
        },
        ReferenceInput {
            id: "same-sequence-other-role".into(),
            molecule: molecule.clone(),
            build: None,
            accession: None,
        },
    ];
    report.validate_bindings_for(&query, &inputs, &off).unwrap();
    replay_off_targets(&query, &inputs, &off, &report).unwrap();
    assert_eq!(report.hits.iter().filter(|h| h.is_intended).count(), 1);
    assert!(
        report
            .hits
            .iter()
            .any(|h| h.reference_index == 1 && !h.is_intended)
    );
    let page = call(
        &c,
        "search_crispr_off_targets",
        json!({"query":query,"references":references,"options":off,"hit_start":1,"max_hits":1}),
    )
    .await;
    assert_eq!(page["artifact_id"], result["artifact_id"]);
    let mut changed = references.clone();
    changed[0]["build"] = json!("new-label");
    let changed = call(
        &c,
        "search_crispr_off_targets",
        json!({"query":query,"references":changed,"options":off}),
    )
    .await;
    assert_ne!(changed["batch_revision"], result["batch_revision"]);
    let mut partial = off.clone();
    partial.max_comparisons = 1;
    let partial = call(
        &c,
        "search_crispr_off_targets",
        json!({"query":query,"references":references,"options":partial}),
    )
    .await;
    assert_eq!(partial["complete"], false);
    assert_eq!(partial["exhaustive_for_selected_records"], false);
    let mut stale = document.clone();
    stale["revision"] = json!("0".repeat(64));
    reject(
        &c,
        "design_crispr_guides",
        json!({"document":stale,"options":options}),
    )
    .await;
    let mut duplicates = references.clone();
    duplicates[1]["id"] = duplicates[0]["id"].clone();
    reject(
        &c,
        "search_crispr_off_targets",
        json!({"query":query,"references":duplicates,"options":off}),
    )
    .await;
    let mut unknown = serde_json::to_value(&off).unwrap();
    unknown["ignore_every_identical_sequence"] = json!(true);
    reject(
        &c,
        "search_crispr_off_targets",
        json!({"query":query,"references":references,"options":unknown}),
    )
    .await;
    assert_eq!(call(&c, "list_documents", json!({})).await, before);
    c.cancel().await.unwrap();
    s.await.unwrap();
}
#[tokio::test]
async fn editing_donors_are_atomically_installed_only_after_independent_target_consistency() {
    let (c, s) = setup().await;
    let imported = call(&c, "import_document", json!({"text":EDIT_TEXT})).await;
    let source = doc(&imported, "source");
    let target = doc(&imported, "target");
    let options = edit_options();
    let design = call(
        &c,
        "design_edit_template",
        json!({"source":source,"target":target,"options":options}),
    )
    .await;
    assert_eq!(design["consistent"], true);
    assert_eq!(design["donor_preview"], "ACCTCGG");
    assert!(design["donor_document"].is_object());
    let full: EditTemplateReport =
        serde_json::from_value(artifact(&c, &design["artifact_id"]).await).unwrap();
    let molecules = collagene_core::formats::sequence::import_text(EDIT_TEXT)
        .unwrap()
        .molecules;
    full.validate_bindings_for(&molecules[0], &molecules[1], &options, None)
        .unwrap();
    assert!(
        edit_template::check(
            &molecules[0],
            &molecules[1],
            &options,
            None,
            full.donor.sequence.as_str()
        )
        .unwrap()
        .consistent
    );
    let donor = json!({"id":design["donor_document"]["id"],"revision":design["donor_document"]["revision"]});
    let checked = call(
        &c,
        "check_edit_template",
        json!({"source":source,"target":target,"donor":donor,"options":options}),
    )
    .await;
    assert_eq!(checked["check"]["consistent"], true);
    let check: EditTemplateCheck =
        serde_json::from_value(artifact(&c, &checked["artifact_id"]).await).unwrap();
    check
        .validate_bindings_for(&molecules[0], &molecules[1], &options, None, "ACCTCGG")
        .unwrap();
    let before = call(&c, "list_documents", json!({})).await;
    let draft = call(
        &c,
        "design_edit_template",
        json!({"source":source,"target":source,"options":options}),
    )
    .await;
    assert_eq!(draft["consistent"], false);
    assert!(draft["donor_document"].is_null());
    assert_eq!(call(&c, "list_documents", json!({})).await, before);
    let bad = call(
        &c,
        "import_document",
        json!({"text":">bad_donor\nACCTCGT\n"}),
    )
    .await;
    let checked = call(
        &c,
        "check_edit_template",
        json!({"source":source,"target":target,"donor":doc(&bad,"bad_donor"),"options":options}),
    )
    .await;
    assert_eq!(checked["check"]["consistent"], false);
    assert_eq!(checked["check"]["right_matches"], false);
    let mut reverse = options.clone();
    reverse.reverse_donor = true;
    let reversed = call(
        &c,
        "design_edit_template",
        json!({"source":source,"target":target,"options":reverse}),
    )
    .await;
    assert_eq!(reversed["consistent"], true);
    assert_eq!(reversed["donor_preview"], "CCGAGGT");
    let foreign_source = call(&c, "import_document", json!({"text":guide_text()})).await;
    let foreign_guides = call(
        &c,
        "design_crispr_guides",
        json!({"document":doc(&foreign_source,"reference"),"options":GuideOptions::default()}),
    )
    .await;
    let foreign =
        json!({"options":GuideOptions::default(),"candidate":foreign_guides["candidates"][0]});
    reject(
        &c,
        "design_edit_template",
        json!({"source":source,"target":target,"options":options,"guide":foreign}),
    )
    .await;
    let mut unknown = serde_json::to_value(&options).unwrap();
    unknown["guaranteed_hdr"] = json!(true);
    reject(
        &c,
        "design_edit_template",
        json!({"source":source,"target":target,"options":unknown}),
    )
    .await;
    c.cancel().await.unwrap();
    s.await.unwrap();
}
#[tokio::test]
async fn full_artifact_capacity_cannot_leave_an_orphan_donor_document() {
    let (c, s) = setup().await;
    let imported = call(&c, "import_document", json!({"text":EDIT_TEXT})).await;
    let guide_docs = call(&c, "import_document", json!({"text":guide_text()})).await;
    let document = doc(&guide_docs, "reference");
    for work in 1..=128 {
        let mut options = GuideOptions::default();
        options.max_comparisons = work;
        call(
            &c,
            "design_crispr_guides",
            json!({"document":document,"options":options}),
        )
        .await;
    }
    let before = call(&c, "list_documents", json!({})).await;
    reject(&c,"design_edit_template",json!({"source":doc(&imported,"source"),"target":doc(&imported,"target"),"options":edit_options()})).await;
    assert_eq!(call(&c, "list_documents", json!({})).await, before);
    c.cancel().await.unwrap();
    s.await.unwrap();
}

#[tokio::test]
async fn supplied_gap_and_screening_pam_scope_keeps_native_paths_coordinates_and_refusals() {
    use collagene_core::sequence::crispr::alignment::{
        AlignmentColumn, AlignmentLimitReached, MODEL_ID,
    };
    let (c, s) = setup().await;
    // Independent literal A-repeat fixture: a 13-base DNA interval has thirteen
    // distinct one-DNA-gap placements against the 12-base query, all retained.
    let text = format!(
        ">exact\n{}GGA\n>bulged\n{}GGA\n>alternate\n{}AGA\n",
        "A".repeat(12),
        "A".repeat(13),
        "A".repeat(12)
    );
    let imported = call(&c, "import_document", json!({"text":text})).await;
    let exact = doc(&imported, "exact");
    let bulged = doc(&imported, "bulged");
    let alternate = doc(&imported, "alternate");
    let before = call(&c, "list_documents", json!({})).await;
    let molecules = collagene_core::formats::sequence::import_text(&text)
        .unwrap()
        .molecules;
    let references = json!([
        {"id":"intended","document":exact,"build":null,"accession":null},
        {"id":"bulged","document":bulged,"build":null,"accession":null},
        {"id":"alternate","document":alternate,"build":null,"accession":null}
    ]);
    let inputs = ["intended", "bulged", "alternate"]
        .into_iter()
        .zip(&molecules)
        .map(|(id, molecule)| ReferenceInput {
            id: id.into(),
            molecule: molecule.clone(),
            build: None,
            accession: None,
        })
        .collect::<Vec<_>>();
    let query = GuideQuery {
        id: "literal12".into(),
        sequence: "A".repeat(12),
    };
    let options = OffTargetOptions {
        profile: NucleaseProfile {
            id: "sdk-gap-fixture".into(),
            version: "literal-1".into(),
            guide_length: 12,
            pam: "GGA".into(),
            pam_side: PamSide::ThreePrime,
            nominal_top_offset: 9,
            nominal_bottom_offset: 9,
            custom: true,
            sources: vec![],
        },
        bulges: BulgePolicy {
            max_dna_bulge_bases: 1,
            max_rna_bulge_bases: 1,
            max_total_bulge_bases: 1,
            ..Default::default()
        },
        pam_policy: PamPolicy::DeclaredScreening {
            pam: "RGA".into(),
            rationale: "Literal broader selected-record screening; design admission remains GGA"
                .into(),
            sources: vec![],
        },
        intended: Some(IntendedTarget {
            reference_id: "intended".into(),
            reference_revision: molecules[0].revision(),
            guide_positions: (0..12).collect(),
            pam_positions: vec![12, 13, 14],
            reverse: false,
        }),
        max_mismatches: 0,
        ..Default::default()
    };
    let result = call(
        &c,
        "search_crispr_off_targets",
        json!({"query":query,"references":references,"options":options,"max_hits":100}),
    )
    .await;
    assert_eq!(result["alignment_model"], MODEL_ID);
    assert_eq!(result["complete"], true);
    assert_eq!(result["exhaustive_for_selected_records"], true);
    let report: OffTargetReport =
        serde_json::from_value(artifact(&c, &result["artifact_id"]).await).unwrap();
    report
        .validate_bindings_for(&query, &inputs, &options)
        .unwrap();
    replay_off_targets(&query, &inputs, &options, &report).unwrap();
    assert_eq!(report.hits.len(), 7);
    assert_eq!(report.hits.iter().filter(|h| h.is_intended).count(), 1);
    let intended = report.hits.iter().find(|h| h.is_intended).unwrap();
    assert_eq!(intended.reference_index, 0);
    assert_eq!(intended.guide_positions, (0..12).collect::<Vec<_>>());
    assert_eq!(intended.pam_positions, vec![12, 13, 14]);
    assert_eq!(intended.design_pam, PamEvidence::Certain);
    assert!(
        intended
            .alignment
            .alignments
            .iter()
            .all(|a| a.dna_bulge_bases + a.rna_bulge_bases == 0)
    );
    let dna = report
        .hits
        .iter()
        .find(|h| h.reference_index == 1 && h.reference_guide.len() == 13)
        .unwrap();
    assert_eq!(dna.guide_positions, (0..13).collect::<Vec<_>>());
    assert_eq!(dna.pam_positions, vec![13, 14, 15]);
    assert!(!dna.reverse && !dna.is_intended && dna.alignment.complete);
    assert_eq!((dna.nominal_top_cut, dna.nominal_bottom_cut), (None, None));
    assert_eq!(dna.alignment.alignments.len(), 13);
    let gap_positions = dna
        .alignment
        .alignments
        .iter()
        .map(|path| {
            assert_eq!(
                (
                    path.min_mismatches,
                    path.max_mismatches,
                    path.dna_bulge_bases,
                    path.rna_bulge_bases
                ),
                (0, 0, 1, 0)
            );
            assert!(path.definite);
            assert_eq!(path.columns.len(), 13);
            path.columns
                .iter()
                .find_map(|column| match column {
                    AlignmentColumn::DnaBulge { target } => Some(*target),
                    _ => None,
                })
                .unwrap()
        })
        .collect::<std::collections::BTreeSet<_>>();
    assert_eq!(
        gap_positions,
        (0..13).collect::<std::collections::BTreeSet<_>>()
    );
    let compact = result["hits"]
        .as_array()
        .unwrap()
        .iter()
        .find(|row| row["id"] == dna.id)
        .unwrap();
    assert_eq!(compact["alignment"]["alignments_retained"], 13);
    assert_eq!(compact["alignment"]["preview_truncated"], true);
    assert_eq!(
        compact["alignment"]["alignments"].as_array().unwrap().len(),
        8
    );
    assert!(
        compact["alignment"]["alignments"]
            .as_array()
            .unwrap()
            .iter()
            .all(|path| path.get("columns").is_none())
    );
    assert!(
        report
            .hits
            .iter()
            .filter(|h| h.reference_index == 2)
            .all(|h| h.design_pam == PamEvidence::Excluded && !h.is_intended && h.pam_certain)
    );
    let mut capped = options.clone();
    capped.bulges.max_alignments_per_hit = 2;
    let partial = call(
        &c,
        "search_crispr_off_targets",
        json!({"query":query,"references":references,"options":capped}),
    )
    .await;
    assert_eq!(partial["complete"], false);
    assert_eq!(partial["exhaustive_for_selected_records"], false);
    let evidence: OffTargetReport =
        serde_json::from_value(artifact(&c, &partial["artifact_id"]).await).unwrap();
    evidence
        .validate_bindings_for(&query, &inputs, &capped)
        .unwrap();
    let hit = evidence
        .hits
        .iter()
        .find(|h| h.reference_index == 1 && h.reference_guide.len() == 13)
        .unwrap();
    assert_eq!(hit.alignment.alignments.len(), 2);
    assert_eq!(
        hit.alignment.limit_reached,
        Some(AlignmentLimitReached::Alignments)
    );
    let mut bounded = options.clone();
    bounded.max_comparisons = 1;
    let partial = call(
        &c,
        "search_crispr_off_targets",
        json!({"query":query,"references":references,"options":bounded}),
    )
    .await;
    assert_eq!(partial["complete"], false);
    assert_eq!(partial["exhaustive_for_selected_records"], false);
    assert_eq!(partial["hit_count"], 0);
    assert!(!partial["limits"].as_array().unwrap().is_empty());
    let partial_report: OffTargetReport =
        serde_json::from_value(artifact(&c, &partial["artifact_id"]).await).unwrap();
    partial_report
        .validate_bindings_for(&query, &inputs, &bounded)
        .unwrap();
    let mut stale = references.clone();
    stale[0]["document"]["revision"] = json!("0".repeat(64));
    reject(
        &c,
        "search_crispr_off_targets",
        json!({"query":query,"references":stale,"options":options}),
    )
    .await;
    let mut oversized = options.clone();
    oversized.bulges.max_total_bulge_bases = 5;
    reject(
        &c,
        "search_crispr_off_targets",
        json!({"query":query,"references":references,"options":oversized}),
    )
    .await;
    let mut oversized = options.clone();
    oversized.bulges.max_retained_alignment_columns = 1_000_001;
    reject(
        &c,
        "search_crispr_off_targets",
        json!({"query":query,"references":references,"options":oversized}),
    )
    .await;
    assert_eq!(call(&c, "list_documents", json!({})).await, before);
    c.cancel().await.unwrap();
    s.await.unwrap();
}
