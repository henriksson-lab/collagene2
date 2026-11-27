//! Conditional source-coordinate tracks. Unresolved framework residues stay
//! unlabeled; a representative alignment never supplies missing assignments.
use super::*;
use crate::{domain_candidates as d, v_mapping};
use std::collections::BTreeSet;

pub const TRACK_ENGINE: &str = "antibody-domain-tracks-2";
const CHART: &str = "https://www.imgt.org/IMGTScientificChart/Numbering/IMGTIGVLsuperfamily.html";
data! {pub enum NumberingScheme {Imgt}}
data! {pub enum RegionDefinition {Imgt}}
data! {pub enum Region {Fr1,Cdr1,Fr2,Cdr2,Fr3,Cdr3,Fr4}}
data! {pub enum AssignmentEvidence {InvariantVReference,ConditionalLoopRule,JFr4Reference}}
data! {pub struct Assignment {
    pub position:Position,pub region:Region,pub evidence:AssignmentEvidence,
}}
data! {pub struct Residue {
    pub source_index:usize,pub residue:char,pub assignment:Option<Assignment>,
}}
data! {pub struct Track {
    pub candidate_id:String,pub locus:v_mapping::Locus,
    pub numbering_scheme:NumberingScheme,pub region_definition:RegionDefinition,
    /// Only present when reference position 1 has an invariant source match.
    /// Absence does not imply that the preceding source is outside the domain.
    pub invariant_position_one:Option<usize>,
    pub first_supported_source_index:usize,pub source_end:usize,
    pub unassigned_residues:usize,pub residues:Vec<Residue>,
}}
data! {pub struct Report {
    pub engine:String,pub source_revision:String,pub rule_revision:String,pub context_revision:String,
    pub evidence:d::Report,pub tracks:Vec<Track>,
}}
pub fn catalog_revision() -> String {
    hash(b"antibody-domain-track-rules-2\0", &(TRACK_ENGINE, CHART,
        "IMGT regions FR1 1..26 CDR1 27..38 FR2 39..55 CDR2 56..65 FR3 66..104 CDR3 105..117 FR4 118..128; decimal insertions at111/112 only; invariant reference framework assignments; conditional loop rules; complete J FR4; unresolved rows retained; no representative start or inferred insertion labels"),65536).expect("fixed track rules")
}
/// Definition selection is separate from numbering selection. No offset-based
/// conversion or unsupported scheme/definition pair is implicitly accepted.
pub fn region_for(
    scheme: &NumberingScheme,
    definition: &RegionDefinition,
    label: &Position,
) -> Result<Region> {
    label.validate_syntax()?;
    match (scheme, definition) {
        (NumberingScheme::Imgt, RegionDefinition::Imgt) => {
            if !matches!(label.insertion, Insertion::None)
                && !matches!(
                    (&label.insertion, label.number),
                    (Insertion::Decimal(_), 111 | 112)
                )
            {
                return Err(Error::Invalid("unsupported IMGT track insertion"));
            }
            Ok(match label.number {
                1..=26 => Region::Fr1,
                27..=38 => Region::Cdr1,
                39..=55 => Region::Fr2,
                56..=65 => Region::Cdr2,
                66..=104 => Region::Fr3,
                105..=117 => Region::Cdr3,
                118..=128 => Region::Fr4,
                _ => return Err(Error::Invalid("unsupported IMGT track position")),
            })
        }
    }
}
fn is_framework(label: &Position) -> bool {
    label.insertion == Insertion::None && matches!(label.number,1..=26|39..=55|66..=104)
}
fn assign(
    rows: &mut [Residue],
    start: usize,
    residue: &NumberedResidue,
    evidence: AssignmentEvidence,
    labels: &mut BTreeSet<String>,
) -> Result<()> {
    let index = residue
        .source_index
        .checked_sub(start)
        .ok_or(Error::Binding)?;
    let row = rows.get_mut(index).ok_or(Error::Binding)?;
    if row.assignment.is_some()
        || row.residue != residue.residue
        || !labels.insert(residue.position.to_string())
    {
        return Err(Error::Binding);
    }
    row.assignment = Some(Assignment {
        position: residue.position.clone(),
        region: region_for(
            &NumberingScheme::Imgt,
            &RegionDefinition::Imgt,
            &residue.position,
        )?,
        evidence,
    });
    Ok(())
}
fn build_track(
    bytes: &[u8],
    report: &d::Report,
    candidate: &d::Candidate,
    cancel: &mut impl FnMut() -> bool,
) -> Result<Track> {
    let invariants = report
        .v_evidence
        .mappings
        .iter()
        .find(|m| m.locus == candidate.locus)
        .and_then(|m| m.projections.iter().find(|p| p.end == candidate.v_endpoint))
        .and_then(|p| p.invariant_correspondences.as_ref())
        .ok_or(Error::Binding)?;
    let framework = invariants
        .iter()
        .filter(|m| is_framework(&m.reference_position))
        .collect::<Vec<_>>();
    let start = framework
        .iter()
        .map(|m| m.source_index)
        .min()
        .ok_or(Error::Binding)?;
    let segment = bytes
        .get(start..candidate.source_end)
        .ok_or(Error::Binding)?;
    let mut rows = Vec::with_capacity(segment.len());
    for (i, &residue) in segment.iter().enumerate() {
        if cancel() {
            return Err(Error::Cancelled);
        }
        rows.push(Residue {
            source_index: start + i,
            residue: residue as char,
            assignment: None,
        });
    }
    let mut labels = BTreeSet::new();
    for m in &framework {
        assign(
            &mut rows,
            start,
            &NumberedResidue {
                source_index: m.source_index,
                residue: m.residue,
                position: m.reference_position.clone(),
            },
            AssignmentEvidence::InvariantVReference,
            &mut labels,
        )?;
    }
    // Reference storage labels inside CDR1/2 are deliberately excluded above:
    // the conditional display rule is the only source of loop labels here.
    for outcome in [&candidate.cdr1, &candidate.cdr2] {
        if let d::LoopAssignment::Numbered(n) = outcome {
            for residue in &n.residues {
                assign(
                    &mut rows,
                    start,
                    residue,
                    AssignmentEvidence::ConditionalLoopRule,
                    &mut labels,
                )?;
            }
        }
    }
    for residue in &candidate.cdr3.residues {
        assign(
            &mut rows,
            start,
            residue,
            AssignmentEvidence::ConditionalLoopRule,
            &mut labels,
        )?;
    }
    for residue in &candidate.fr4 {
        assign(
            &mut rows,
            start,
            residue,
            AssignmentEvidence::JFr4Reference,
            &mut labels,
        )?;
    }
    Ok(Track {
        candidate_id: candidate.id.clone(),
        locus: candidate.locus.clone(),
        numbering_scheme: NumberingScheme::Imgt,
        region_definition: RegionDefinition::Imgt,
        invariant_position_one: framework
            .iter()
            .find(|m| m.reference_position.number == 1)
            .map(|m| m.source_index),
        first_supported_source_index: start,
        source_end: candidate.source_end,
        unassigned_residues: rows.iter().filter(|r| r.assignment.is_none()).count(),
        residues: rows,
    })
}
fn context_revision(
    source: &str,
    corpus: &str,
    options: &d::Options,
    limits: &d::Limits,
) -> Result<String> {
    hash(
        b"antibody-domain-track-context-2\0",
        &(
            TRACK_ENGINE,
            catalog_revision(),
            source,
            corpus,
            options,
            limits,
            d::CANDIDATE_ENGINE,
            v_mapping::MAPPING_ENGINE,
            profile::PROFILE_ENGINE,
            profile::ENDPOINT_ENGINE,
            references::REFERENCE_ENGINE,
        ),
        65536,
    )
}
/// Cheap source/options admission and request binding; no alignment is run.
pub fn request_revision(
    source: &Molecule,
    options: &d::Options,
    limits: &d::Limits,
) -> Result<String> {
    d::validate_inputs(options, limits)?;
    context_revision(
        &source_revision(source)?,
        &references::human_references()?.revision()?,
        options,
        limits,
    )
}
/// Compute the candidate evidence once, then derive all tracks from that exact
/// immutable result. These are conditional reference/rule assignments, not a
/// calibrated domain classifier or independently certified biological numbering.
pub fn analyze(
    source: &Molecule,
    options: &d::Options,
    limits: &d::Limits,
    cancel: &mut impl FnMut() -> bool,
) -> Result<Report> {
    let evidence = d::find_candidates(source, options, limits, cancel)?;
    let rule_revision = catalog_revision();
    let context_revision = context_revision(
        &evidence.source_revision,
        &evidence.v_evidence.corpus_revision,
        options,
        limits,
    )?;
    let mut tracks = Vec::with_capacity(evidence.candidates.len());
    for candidate in &evidence.candidates {
        if cancel() {
            return Err(Error::Cancelled);
        }
        tracks.push(build_track(
            source.sequence.as_str().as_bytes(),
            &evidence,
            candidate,
            cancel,
        )?);
    }
    if cancel() {
        return Err(Error::Cancelled);
    }
    Ok(Report {
        engine: TRACK_ENGINE.into(),
        source_revision: evidence.source_revision.clone(),
        rule_revision,
        context_revision,
        evidence,
        tracks,
    })
}
impl Report {
    /// Source-free current-format/structural checks for trusted local worker replies.
    /// This does not establish scores or source letters; imported evidence must
    /// still use validate_for against the independently resolved actual source.
    pub fn validate_binding(&self, expected: &str) -> Result<()> {
        let e = &self.evidence;
        d::validate_inputs(&e.options, &e.limits)?;
        let corpus = references::human_references()?.revision()?;
        if self.engine != TRACK_ENGINE
            || self.source_revision != e.source_revision
            || e.source_revision.len() != 64
            || !e.source_revision.bytes().all(|b| b.is_ascii_hexdigit())
            || self.rule_revision != catalog_revision()
            || self.context_revision != expected
            || context_revision(&e.source_revision, &corpus, &e.options, &e.limits)? != expected
            || e.engine != d::CANDIDATE_ENGINE
            || e.v_evidence.engine != v_mapping::MAPPING_ENGINE
            || e.v_evidence.source_revision != e.source_revision
            || e.v_evidence.corpus_revision != corpus
            || e.v_evidence.scoring != e.options.v_scoring
            || e.v_evidence.limits != e.limits.v_mapping
            || e.candidates.len() > e.limits.max_candidates
            || e.candidates.len() != self.tracks.len()
            || e.qualified_pairs_in_retained_v < e.candidates.len()
            || e.pair_retention_truncated != (e.qualified_pairs_in_retained_v > e.candidates.len())
            || e.windows_examined > e.limits.max_scan_windows
            || e.v_retention_truncated
                != e.v_evidence.mappings.iter().any(|m| m.alignment.truncated)
            || e.v_evidence.mappings.len() != 3
            || e.v_evidence.mappings.iter().any(|m| {
                m.projections.len() > profile::MAX_ENDPOINTS
                    || m.projections.len() != m.alignment.endpoints.len()
                    || m.alignment.endpoints.len() > e.limits.v_mapping.max_endpoints_per_locus
                    || m.alignment.qualified_endpoints < m.alignment.endpoints.len()
                    || m.alignment.truncated
                        != (m.alignment.qualified_endpoints > m.alignment.endpoints.len())
                    || m.projections.iter().any(|p| {
                        p.invariant_correspondences
                            .as_ref()
                            .is_some_and(|i| i.len() > 104)
                    })
            })
        {
            return Err(Error::Binding);
        }
        let mut observed = std::collections::BTreeMap::new();
        for (candidate, track) in e.candidates.iter().zip(&self.tracks) {
            if track.first_supported_source_index >= track.source_end
                || track.source_end > MAX_SOURCE_RESIDUES
                || track.residues.len() != track.source_end - track.first_supported_source_index
                || candidate.cdr3.residues.len() > MAX_SOURCE_RESIDUES
                || candidate.fr4.len() > 11
                || [&candidate.cdr1, &candidate.cdr2]
                    .iter()
                    .any(|o| matches!(o,d::LoopAssignment::Numbered(n) if n.residues.len()>12))
            {
                return Err(Error::Binding);
            }
            let mut bytes = vec![b'X'; track.source_end];
            for (i, row) in track.residues.iter().enumerate() {
                if row.source_index != track.first_supported_source_index + i
                    || !row.residue.is_ascii()
                    || row.residue == '*'
                    || collagene_model::protein_choices(row.residue as u8).is_none()
                {
                    return Err(Error::Binding);
                }
                if observed
                    .insert(row.source_index, row.residue)
                    .is_some_and(|old| old != row.residue)
                {
                    return Err(Error::Binding);
                }
                bytes[row.source_index] = row.residue as u8;
            }
            if *track != build_track(&bytes, e, candidate, &mut || false)? {
                return Err(Error::Binding);
            }
        }
        Ok(())
    }
    pub fn validate_for(
        &self,
        source: &Molecule,
        options: &d::Options,
        limits: &d::Limits,
    ) -> Result<()> {
        if self.engine != TRACK_ENGINE
            || self.tracks.len() > d::MAX_CANDIDATES
            || self
                .tracks
                .iter()
                .any(|t| t.residues.len() > MAX_SOURCE_RESIDUES)
            || *self != analyze(source, options, limits, &mut || false)?
        {
            return Err(Error::Binding);
        }
        Ok(())
    }
}
