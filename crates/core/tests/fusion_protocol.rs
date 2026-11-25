//! Real binary worker execution, selector/version gates and request binding.
use collagene_core::{
    Molecule, Topology, golden_gate_fusion_design as fusion,
    protocol::{self, AnalysisJob, AnalysisResult, Request},
    verification::PlanInput,
};
use serde_json::json;
fn fixture() -> (Vec<PlanInput>, fusion::Request) {
    let inputs = [("a", "ATGACCCCGGGGTTTCCAGA"), ("b", "ATCCGTTAGCGACCTTGGAA")]
        .into_iter()
        .map(|(id, dna)| {
            let molecule = Molecule::dna(id, dna, Topology::Linear).unwrap();
            PlanInput {
                id: id.into(),
                revision: molecule.revision(),
                molecule,
            }
        })
        .collect();
    let request = serde_json::from_value(json!({
        "base":{"base":{"name":"Fixed independent AB","fragments":[
            {"role":"left","input_id":"a","start":0,"end":20,"reverse":false},
            {"role":"right","input_id":"b","start":0,"end":20,"reverse":false}],"circular":true},
            "edits":[],"rows":[
                {"DesignPcr":{"role":"left","input_id":"a","reverse":false,"forward_length":8,"reverse_length":8}},
                {"DesignPcr":{"role":"right","input_id":"b","reverse":false,"forward_length":8,"reverse_length":8}}],
            "enzyme":"BsaI","clamp":"ATGC","spacer":"A","terminal_right_overhang":null,"allow_nicked":false},
        "windows":[{"junction":0,"min_shift":-1,"max_shift":1}],
        "constraints":{"reject_palindromes":false,"distinct_reverse_complements":false,"minimum_hamming":0},
        "limits":{"max_states":1000,"max_candidates":32}
    })).unwrap();
    (inputs, request)
}
fn selector(bytes: &[u8]) -> u32 {
    postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[42..])
        .unwrap()
        .0
        .2
}
fn exchange(job: AnalysisJob, revision: String, job_tag: u32, result_tag: u32) -> AnalysisResult {
    let encoded = protocol::encode_request(&Request {
        id: 601,
        revision,
        job,
    })
    .unwrap();
    assert_eq!(selector(&encoded), job_tag);
    let mut old = encoded.clone();
    old[4..6].copy_from_slice(&59u16.to_le_bytes());
    assert!(protocol::process_request(&old).is_err());
    let response = protocol::process_request(&encoded).unwrap();
    assert_eq!(selector(&response), result_tag);
    let mut old = response.clone();
    old[4..6].copy_from_slice(&59u16.to_le_bytes());
    assert!(protocol::decode_response(&old).is_err());
    protocol::decode_response(&response).unwrap().result
}
#[test]
fn signed_search_and_explicit_preparation_use_new_append_only_worker_selectors() {
    let (inputs, request) = fixture();
    let revision = fusion::inspect_revision(&inputs, &request).unwrap();
    let AnalysisResult::GoldenGateFusions(report) = exchange(
        AnalysisJob::SearchGoldenGateFusions {
            inputs: inputs.clone(),
            request: request.clone(),
        },
        revision,
        136,
        128,
    ) else {
        panic!("fusion search response required")
    };
    assert!(report.search_complete);
    let candidate = report
        .candidates
        .iter()
        .find(|c| c.shifts[0] == -1)
        .unwrap();
    let selection = fusion::Selection {
        context_revision: report.context_revision.clone(),
        candidate_id: candidate.id.clone(),
    };
    let revision = fusion::materialize_revision(&inputs, &request, &selection).unwrap();
    let AnalysisResult::GoldenGateFusionMaterialized(materialized) = exchange(
        AnalysisJob::MaterializeGoldenGateFusions {
            inputs: inputs.clone(),
            request: request.clone(),
            selection: selection.clone(),
        },
        revision,
        137,
        129,
    ) else {
        panic!("materialized original-source preparation required")
    };
    assert_eq!(
        materialized.original_intent.target.sequence.as_str(),
        "ATGACCCCGGGGTTTCCAGAATCCGTTAGCGACCTTGGAA"
    );
    assert_eq!(materialized.recipe.fragments.len(), 2);
    let stale = protocol::encode_request(&Request {
        id: 602,
        revision: "0".repeat(64),
        job: AnalysisJob::MaterializeGoldenGateFusions {
            inputs,
            request,
            selection,
        },
    })
    .unwrap();
    assert!(protocol::process_request(&stale).is_err());
}

#[test]
fn empirical_scores_are_condition_bound_and_keep_missing_counts_unknown_over_binary() {
    use collagene_core::ligation_fidelity as fidelity;
    let mut dataset: fidelity::Dataset = serde_json::from_value(json!({
        "schema_version":1,"overhang_length":3,"content_revision":"",
        "conditions":{"ligase":"synthetic","buffer":"synthetic","temperature_protocol":"not measured","incubation":"not measured","restriction_enzymes":[],"substrate_context":"independent artificial fixture"},
        "provenance":{"source_url":"https://example.org/synthetic","citation":"Independent fixture","license":"MIT","attribution":"Collagene 2 fixtures","source_blake3":blake3::hash(b"fixture").to_hex().to_string(),"normalization":"literal directed counts"},
        "cells":[{"end_5to3":"AAC","partner_5to3":"GTT","count":9}]
    })).unwrap();
    dataset.content_revision = fidelity::content_revision(&dataset).unwrap();
    let ends = vec!["AAC".into()];
    let revision = fidelity::revision(&dataset, &ends).unwrap();
    let AnalysisResult::LigationFidelity(report) = exchange(
        AnalysisJob::ScoreLigationFidelity {
            dataset: dataset.clone(),
            ends: ends.clone(),
        },
        revision,
        138,
        130,
    ) else {
        panic!("conditional estimate required")
    };
    assert_eq!(report.overall_fidelity, None);
    assert_eq!(report.pairs[0].orientations[0].correct_count, Some(9));
    let mut changed = dataset;
    changed.conditions.buffer = "different declared buffer".into();
    assert!(fidelity::revision(&changed, &ends).is_err());
}
