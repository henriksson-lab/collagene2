//! Development-only independent chain-identity assessment; no network/runtime tools.
use collagene_antibody_numbering::{
    domain_candidates as d, domain_tracks, profile, references, v_mapping as v,
};
use collagene_model::{Alphabet, Molecule, Sequence, Topology};
use serde::Deserialize;
use serde_json::{Value, json};
use std::collections::BTreeMap;
const DATA: &str = include_str!("../data/rcsb-assessment/panel.json");
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Input {
    url: String,
    sha256: String,
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Source {
    pdb: String,
    entity: String,
    description: String,
    sequence_source: String,
    entry_source: String,
}
#[derive(Clone, Copy, Deserialize)]
enum Truth {
    IgHeavy,
    IgLight,
    Other,
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Case {
    id: String,
    truth: Truth,
    sequence: String,
    sequence_sha256: String,
    sources: Vec<Source>,
    reason: Option<String>,
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Panel {
    version: u16,
    license: String,
    license_url: String,
    retrieved: String,
    truth_scope: String,
    inputs: BTreeMap<String, Input>,
    cases: Vec<Case>,
    excluded: Vec<Case>,
}
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let panel: Panel = serde_json::from_str(DATA)?;
    if panel.version != 1
        || panel.license != "CC0-1.0"
        || panel.truth_scope != "deposited-chain-identity-only-not-residue-numbering"
    {
        return Err("unsupported assessment panel".into());
    }
    for case in panel.cases.iter().chain(&panel.excluded) {
        if case.sources.is_empty() || case.sequence_sha256.len() != 64 {
            return Err("assessment provenance".into());
        }
        for source in &case.sources {
            if source.pdb.is_empty() || source.entity.is_empty() || source.description.is_empty() {
                return Err("source identity".into());
            }
            for key in [&source.sequence_source, &source.entry_source] {
                let input = panel.inputs.get(key).ok_or("missing pinned source")?;
                if !input.url.starts_with("https://data.rcsb.org/") || input.sha256.len() != 64 {
                    return Err("source URL/hash".into());
                }
            }
        }
    }
    if panel.cases.iter().any(|c| c.reason.is_some())
        || panel.excluded.iter().any(|c| c.reason.is_none())
    {
        return Err("assessment exclusion ledger".into());
    }
    let limits = d::Limits {
        v_mapping: v::Limits {
            max_cells: profile::MAX_CELLS,
            max_endpoints_per_locus: profile::MAX_ENDPOINTS,
            max_trace_states: profile::MAX_TRACE_STATES,
        },
        max_scan_windows: d::MAX_SCAN_WINDOWS,
        max_candidates: d::MAX_CANDIDATES,
    };
    let mut presets = vec![];
    for max_fr4_mismatches in [0, 1, 3] {
        let options = d::Options {
            v_scoring: v::Scoring {
                observed_reward: 5,
                frequency_bonus: 8,
                unseen_penalty: -6,
                gap_open: -12,
                gap_extend: -2,
                minimum_score: 300,
            },
            min_cdr3_residues: 2,
            max_cdr3_residues: 60,
            max_fr4_mismatches,
        };
        let mut rows = vec![];
        let mut missed_positives = vec![];
        let mut negative_calls = vec![];
        let mut positive_count = 0;
        let mut negative_count = 0;
        for case in &panel.cases {
            let molecule = Molecule {
                name: case.id.clone(),
                sequence: Sequence::new(Alphabet::Protein, &case.sequence)?,
                topology: Topology::Linear,
                features: vec![],
                metadata: Default::default(),
                original_genbank: None,
            };
            let track_report = domain_tracks::analyze(&molecule, &options, &limits, &mut || false)?;
            let report = &track_report.evidence;
            let positive = !matches!(case.truth, Truth::Other);
            let matched_class = report.candidates.iter().any(|c| match case.truth {
                Truth::IgHeavy => c.locus == v::Locus::Heavy,
                Truth::IgLight => matches!(c.locus, v::Locus::Kappa | v::Locus::Lambda),
                Truth::Other => false,
            });
            if positive {
                positive_count += 1;
                if !matched_class {
                    missed_positives.push(case.id.clone());
                }
            } else {
                negative_count += 1;
                if report.qualified_pairs_in_retained_v > 0 {
                    negative_calls.push(case.id.clone());
                }
            }
            let candidates=report.candidates.iter().map(|c| {
                let locus=report.v_evidence.mappings.iter().find(|m|m.locus==c.locus).unwrap();
                let path=locus.alignment.endpoints.iter().find(|e|e.end==c.v_endpoint).unwrap();
                json!({"locus":c.locus,"representative_v_start":c.representative_v_start,"v_projection_issue":c.v_projection_issue,"source_end":c.source_end,"v_end":c.v_endpoint,"v_score":path.upper.score,
                    "cdr1":c.cdr1,"cdr2":c.cdr2,"cdr3_start":c.cdr3.start,"cdr3_end":c.cdr3.end,"cdr3_sequence":&case.sequence[c.cdr3.start..c.cdr3.end],
                    "fr4_sequence":c.fr4.iter().map(|r|r.residue).collect::<String>(),"j_matches":c.j_matches,
                    "track":track_report.tracks.iter().find(|t|t.candidate_id==c.id)})
            }).collect::<Vec<Value>>();
            let loci=report.v_evidence.mappings.iter().map(|m|json!({"locus":m.locus,"reachable":m.alignment.reachable_endpoints,"qualified":m.alignment.qualified_endpoints,
                "retained":m.alignment.endpoints.len(),"best_retained_score":m.alignment.endpoints.iter().map(|e|e.upper.score).max(),
                "projected":m.projections.iter().filter(|p|p.projection.is_some()).count()})).collect::<Vec<_>>();
            rows.push(json!({"case":case.id,"positive_chain":positive,"matching_chain_class_candidate":matched_class,"length":case.sequence.len(),"source_revision":report.source_revision,
                "candidate_count_in_retained_v":report.qualified_pairs_in_retained_v,"v_truncated":report.v_retention_truncated,"pairs_truncated":report.pair_retention_truncated,
                "blocked_v_endpoints":report.blocked_v_endpoints,"unpaired_v_endpoints":report.v_endpoints_without_j,"loci":loci,"candidates":candidates}));
        }
        presets.push(json!({"options":options,"limits":limits,"positive_chains":positive_count,"negative_chains":negative_count,"missed_positive_chains":missed_positives,"negative_chains_with_candidates":negative_calls,"rows":rows}));
    }
    let output = json!({"panel_version":panel.version,"panel_revision":blake3::hash(DATA.as_bytes()).to_hex().to_string(),"license":panel.license,"license_url":panel.license_url,"retrieved":panel.retrieved,
        "truth_scope":panel.truth_scope,"excluded_cases":panel.excluded.len(),"engines":{"tracks":domain_tracks::TRACK_ENGINE,"candidates":d::CANDIDATE_ENGINE,"v_mapping":v::MAPPING_ENGINE,"endpoint":profile::ENDPOINT_ENGINE,"references":references::REFERENCE_ENGINE},
        "presets":presets});
    serde_json::to_writer_pretty(std::io::stdout().lock(), &output)?;
    Ok(())
}
