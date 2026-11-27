//! Conditional V/J pairing and IMGT CDR3 reconstruction. This composes existing
//! verified reference calculations; it is not a calibrated domain classifier.
use super::*;
use crate::{references, v_mapping};
use std::collections::BTreeMap;
pub const CANDIDATE_ENGINE: &str = "antibody-vj-candidates-1";
pub const MAX_CANDIDATES: usize = 64;
pub const MAX_SCAN_WINDOWS: usize = 2_000_000;
data! {pub struct Options {
    pub v_scoring:v_mapping::Scoring,
    pub min_cdr3_residues:usize,pub max_cdr3_residues:usize,
    /// Substitutions outside invariant J-anchor/motif positions only. No indels.
    pub max_fr4_mismatches:usize,
}}
data! {pub struct Limits {
    pub v_mapping:v_mapping::Limits,
    pub max_scan_windows:usize,pub max_candidates:usize,
}}
data! {pub struct Mismatch {pub source_index:usize,pub observed:char,pub reference:char}}
data! {pub struct JMatch {
    pub reference_protein:String,pub allele_refs:Vec<String>,pub mismatches:Vec<Mismatch>,
}}
data! {pub struct Candidate {
    pub id:String,pub locus:v_mapping::Locus,pub v_endpoint:usize,
    pub source_start:usize,pub source_end:usize,
    pub cdr3:ConditionalCdr3,pub fr4:Vec<NumberedResidue>,pub j_matches:Vec<JMatch>,
}}
data! {pub struct VEndpoint {pub locus:v_mapping::Locus,pub end:usize}}
data! {pub struct BlockedV {
    pub endpoint:VEndpoint,pub issue:v_mapping::ProjectionIssue,
}}
data! {pub struct Report {
    pub engine:String,pub source_revision:String,pub context_revision:String,
    pub options:Options,pub limits:Limits,pub v_evidence:v_mapping::Report,
    pub windows_examined:usize,
    /// Exact pair count within retained, uniquely projected V endpoints only.
    pub qualified_pairs_in_retained_v:usize,pub matching_j_groups:usize,
    pub v_retention_truncated:bool,pub pair_retention_truncated:bool,
    pub blocked_v_endpoints:Vec<BlockedV>,pub v_endpoints_without_j:Vec<VEndpoint>,
    pub candidates:Vec<Candidate>,
}}
fn validate_inputs(options: &Options, limits: &Limits) -> Result<()> {
    if options.min_cdr3_residues < 2
        || options.max_cdr3_residues < options.min_cdr3_residues
        || options.max_cdr3_residues > MAX_CDR3_RESIDUES
        || options.max_fr4_mismatches > 3
    {
        return Err(Error::Invalid("V/J loop or substitution scope"));
    }
    if limits.max_candidates > MAX_CANDIDATES
        || limits.max_scan_windows == 0
        || limits.max_scan_windows > MAX_SCAN_WINDOWS
    {
        return Err(Error::Limit("V/J retention or scan work"));
    }
    Ok(())
}
fn source_index(locus: &v_mapping::Locus) -> usize {
    match locus {
        v_mapping::Locus::Heavy => 0,
        v_mapping::Locus::Kappa => 1,
        v_mapping::Locus::Lambda => 2,
    }
}
/// Exhaustive complete-FR4 substitution windows inside the declared loop range
/// for retained, unique best-per-end V projections. Suboptimal/tied V paths and
/// omitted V endpoints are not silently treated as searched or absent domains.
pub fn find_candidates(
    source: &Molecule,
    options: &Options,
    limits: &Limits,
    cancel: &mut impl FnMut() -> bool,
) -> Result<Report> {
    if cancel() {
        return Err(Error::Cancelled);
    }
    validate_inputs(options, limits)?;
    let v_evidence = v_mapping::map_query(source, &options.v_scoring, &limits.v_mapping, cancel)?;
    let corpus = references::human_references()?;
    if corpus.revision()? != v_evidence.corpus_revision {
        return Err(Error::Binding);
    }
    let context_revision = hash(
        b"antibody-vj-candidate-context-1\0",
        &(
            CANDIDATE_ENGINE,
            v_mapping::MAPPING_ENGINE,
            profile::ENDPOINT_ENGINE,
            ENGINE,
            catalog_revision(),
            references::REFERENCE_ENGINE,
            &v_evidence.context_revision,
            options,
            limits,
            "retained unique V ends only; all complete J FR4 windows; literal locus anchor and G119/G121; bounded substitutions; no indels; source-order retention; no biological classification",
        ),
        65536,
    )?;
    let mut templates: [BTreeMap<String, Vec<String>>; 3] =
        std::array::from_fn(|_| BTreeMap::new());
    for entry in corpus.j_entries() {
        templates[entry.source]
            .entry(entry.fr4_protein.clone())
            .or_default()
            .push(entry.allele_ref.clone());
    }
    let bytes = source.sequence.as_str().as_bytes();
    let mut windows_examined = 0;
    let mut qualified_pairs_in_retained_v = 0;
    let mut matching_j_groups = 0;
    let mut blocked_v_endpoints = vec![];
    let mut v_endpoints_without_j = vec![];
    let mut candidates = vec![];
    let v_retention_truncated = v_evidence.mappings.iter().any(|m| m.alignment.truncated);
    for mapping in &v_evidence.mappings {
        let index = source_index(&mapping.locus);
        for (endpoint, projected) in mapping.alignment.endpoints.iter().zip(&mapping.projections) {
            if cancel() {
                return Err(Error::Cancelled);
            }
            if endpoint.end != projected.end {
                return Err(Error::Binding);
            }
            let Some(projection) = &projected.projection else {
                blocked_v_endpoints.push(BlockedV {
                    endpoint: VEndpoint {
                        locus: mapping.locus.clone(),
                        end: endpoint.end,
                    },
                    issue: projected.projection_issue.clone().ok_or(Error::Binding)?,
                });
                continue;
            };
            let cysteine = projection
                .matches
                .iter()
                .find(|r| r.reference_position == position(104))
                .ok_or(Error::Binding)?
                .source_index;
            if cysteine + 1 != endpoint.end || bytes[cysteine] != b'C' {
                return Err(Error::Binding);
            }
            let loop_start = cysteine + 1;
            let first = loop_start + options.min_cdr3_residues;
            let stop = bytes[loop_start..]
                .iter()
                .position(|&b| b == b'*')
                .map_or(bytes.len(), |n| loop_start + n);
            let last = (loop_start + options.max_cdr3_residues).min(stop);
            let mut found = false;
            for j_start in first..=last {
                if cancel() {
                    return Err(Error::Cancelled);
                }
                let mut j_matches = vec![];
                for (reference_protein, allele_refs) in &templates[index] {
                    let Some(observed) = bytes.get(j_start..j_start + reference_protein.len())
                    else {
                        continue;
                    };
                    windows_examined += 1;
                    if windows_examined > limits.max_scan_windows {
                        return Err(Error::Limit("V/J scan windows"));
                    }
                    if observed[0] != if index == 0 { b'W' } else { b'F' }
                        || observed[1] != b'G'
                        || observed[3] != b'G'
                        || observed.contains(&b'*')
                    {
                        continue;
                    }
                    let mismatches = observed
                        .iter()
                        .zip(reference_protein.bytes())
                        .enumerate()
                        .filter_map(|(i, (&a, b))| {
                            (a != b).then_some(Mismatch {
                                source_index: j_start + i,
                                observed: a as char,
                                reference: b as char,
                            })
                        })
                        .collect::<Vec<_>>();
                    if mismatches.len() <= options.max_fr4_mismatches {
                        j_matches.push(JMatch {
                            reference_protein: reference_protein.clone(),
                            allele_refs: allele_refs.clone(),
                            mismatches,
                        });
                    }
                }
                if j_matches.is_empty() {
                    continue;
                }
                found = true;
                qualified_pairs_in_retained_v += 1;
                matching_j_groups += j_matches.len();
                if candidates.len() >= limits.max_candidates {
                    continue;
                }
                let fr4_length = j_matches[0].reference_protein.len();
                // All current reference tails in a locus have one validated
                // length. Corpus validation refuses an unsupported extension.
                if j_matches
                    .iter()
                    .any(|j| j.reference_protein.len() != fr4_length)
                {
                    return Err(Error::Binding);
                }
                let end = j_start + fr4_length;
                let anchors = Anchors {
                    source_revision: v_evidence.source_revision.clone(),
                    cysteine_index: cysteine,
                    j_anchor_index: j_start,
                    provenance: format!(
                        "{}; conditional {:?} V end {}; V context {}; J corpus {}",
                        CANDIDATE_ENGINE,
                        mapping.locus,
                        endpoint.end,
                        v_evidence.context_revision,
                        v_evidence.corpus_revision
                    ),
                };
                let cdr3 =
                    number_validated_imgt_cdr3(source, &anchors, &v_evidence.source_revision)?;
                let fr4 = bytes[j_start..end]
                    .iter()
                    .enumerate()
                    .map(|(i, &b)| NumberedResidue {
                        source_index: j_start + i,
                        residue: b as char,
                        position: position(118 + i as u16),
                    })
                    .collect();
                let id = hash(
                    b"antibody-vj-candidate-id-1\0",
                    &(
                        &context_revision,
                        &mapping.locus,
                        endpoint.end,
                        j_start,
                        end,
                    ),
                    65536,
                )?;
                candidates.push(Candidate {
                    id,
                    locus: mapping.locus.clone(),
                    v_endpoint: endpoint.end,
                    source_start: endpoint.upper.start,
                    source_end: end,
                    cdr3,
                    fr4,
                    j_matches,
                });
            }
            if !found {
                v_endpoints_without_j.push(VEndpoint {
                    locus: mapping.locus.clone(),
                    end: endpoint.end,
                });
            }
        }
    }
    if cancel() {
        return Err(Error::Cancelled);
    }
    let pair_retention_truncated = qualified_pairs_in_retained_v > candidates.len();
    Ok(Report {
        engine: CANDIDATE_ENGINE.into(),
        source_revision: v_evidence.source_revision.clone(),
        context_revision,
        options: options.clone(),
        limits: limits.clone(),
        v_evidence,
        windows_examined,
        qualified_pairs_in_retained_v,
        matching_j_groups,
        v_retention_truncated,
        pair_retention_truncated,
        blocked_v_endpoints,
        v_endpoints_without_j,
        candidates,
    })
}
impl Report {
    pub fn validate_for(
        &self,
        source: &Molecule,
        options: &Options,
        limits: &Limits,
    ) -> Result<()> {
        if self.engine != CANDIDATE_ENGINE
            || self.candidates.len() > MAX_CANDIDATES
            || self.candidates.iter().any(|c| {
                c.fr4.len() > 11
                    || c.cdr3.residues.len() > MAX_CDR3_RESIDUES
                    || c.j_matches.len() > 23
            })
            || *self != find_candidates(source, options, limits, &mut || false)?
        {
            return Err(Error::Binding);
        }
        Ok(())
    }
}
