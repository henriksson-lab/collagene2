//! Conditional source-coordinate tracks. Unresolved framework residues stay
//! unlabeled; a representative alignment never supplies missing assignments.
use super::*;
use crate::{domain_candidates as d, v_mapping};
use std::collections::BTreeSet;

pub const TRACK_ENGINE: &str = "antibody-domain-tracks-6";
const CHART: &str = "https://www.imgt.org/IMGTScientificChart/Numbering/IMGTIGVLsuperfamily.html";
data! {
    pub enum NumberingScheme {
        Imgt,
        Kabat1991,
        Chothia1997
    }
}
data! {
    pub enum RegionDefinition {
        Imgt,
        Kabat1991,
        ChothiaCanonical,
        North2011
    }
}
data! {
    pub struct Selection {
        pub numbering_scheme: NumberingScheme,
        pub region_definition: RegionDefinition
    }
}
impl Selection {
    pub fn validate_for_references(
        &self,
        references: &v_mapping::ReferenceSelection,
    ) -> Result<()> {
        references.validate()?;
        if references.species == crate::references::Species::MouseC57Bl6
            && (self.numbering_scheme != NumberingScheme::Imgt
                || self.region_definition != RegionDefinition::Imgt)
        {
            return Err(Error::Invalid(
                "mouse C57BL/6 admission supports IMGT numbering and IMGT regions only",
            ));
        }
        Ok(())
    }
}
data! {
    pub enum Region {
        Fr1,
        Cdr1,
        Fr2,
        Cdr2,
        Fr3,
        Cdr3,
        Fr4
    }
}
data! {
    pub enum AssignmentEvidence {
        InvariantVReference,
        ConditionalLoopRule,
        JFr4Reference,
        PublishedFrameworkCorrespondence,
        PublishedWindowRule
    }
}
data! {
    pub struct Assignment {
        pub position: Position,
        pub evidence: AssignmentEvidence,
    }
}
data! {
    pub struct Residue {
        pub source_index: usize,
        pub residue: char,
        pub assignment: Option<Assignment>,
        pub region: Option<Region>,
    }
}
data! {
    pub enum CdrStatus {
        Mapped { start: usize, end: usize },
        Unresolved
    }
}
data! {
    pub struct CdrRegion {
        pub kind: Region,
        pub status: CdrStatus
    }
}
data! {
    pub struct Track {
        pub candidate_id: String,
        pub locus: v_mapping::Locus,
        pub numbering_scheme: NumberingScheme,
        pub region_definition: RegionDefinition,
        /// Only present when reference position 1 has an invariant source match.
        /// Absence does not imply that the preceding source is outside the domain.
        pub invariant_position_one: Option<usize>,
        pub first_supported_source_index: usize,
        pub source_end: usize,
        pub unassigned_residues: usize,
        pub unassigned_regions: usize,
        pub issues: Vec<crate::alternate_tracks::Issue>,
        pub cdrs: Vec<CdrRegion>,
        pub residues: Vec<Residue>,
    }
}
data! {
    pub struct Report {
        pub engine: String,
        pub source_revision: String,
        pub rule_revision: String,
        pub context_revision: String,
        pub evidence: d::Report,
        pub selection: Selection,
        pub tracks: Vec<Track>,
    }
}
pub fn catalog_revision() -> String {
    hash(b"antibody-domain-track-rules-6\0", &(TRACK_ENGINE, CHART,crate::alternate_tracks::catalog_revision(),
        "IMGT regions FR1 1..26 CDR1 27..38 FR2 39..55 CDR2 56..65 FR3 66..104 CDR3 105..117 FR4 118..128; decimal insertions at111/112 only; invariant reference framework assignments; conditional loop rules; complete J FR4; unresolved rows retained; no representative start or inferred insertion labels"),65536).expect("fixed track rules")
}
/// Current IMGT region boundaries, applied to established IMGT labels only.
pub fn imgt_region(label: &Position) -> Result<Region> {
    label.validate_syntax()?;
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
/// Kabat regions are derived through separately established Kabat labels, not
/// numeric tests on a label expressed in another scheme.
pub fn kabat_region(locus: &v_mapping::Locus, label: &Position) -> Result<Region> {
    label.validate_syntax()?;
    if matches!(label.insertion, Insertion::Decimal(_)) {
        return Err(Error::Invalid("decimal Kabat insertion"));
    }
    let ranges = if *locus == v_mapping::Locus::Heavy {
        [
            (1, 30),
            (31, 35),
            (36, 49),
            (50, 65),
            (66, 94),
            (95, 102),
            (103, 113),
        ]
    } else {
        [
            (1, 23),
            (24, 34),
            (35, 49),
            (50, 56),
            (57, 88),
            (89, 97),
            (98, 107),
        ]
    };
    let regions = [
        Region::Fr1,
        Region::Cdr1,
        Region::Fr2,
        Region::Cdr2,
        Region::Fr3,
        Region::Cdr3,
        Region::Fr4,
    ];
    ranges
        .iter()
        .position(|&(a, b)| label.number >= a && label.number <= b)
        .map(|i| regions[i].clone())
        .ok_or(Error::Invalid("unsupported Kabat track position"))
}
/// Canonical-loop delimitations in the primary IMGT comparison charts, using
/// their later Chothia light-chain numbering. This is a specific definition,
/// not a synonym for every later convention called "Chothia".
pub fn chothia_canonical_region(locus: &v_mapping::Locus, label: &Position) -> Result<Region> {
    label.validate_syntax()?;
    if matches!(label.insertion, Insertion::Decimal(_)) {
        return Err(Error::Invalid("decimal Chothia insertion"));
    }
    let n = label.number;
    if *locus == v_mapping::Locus::Heavy {
        return Ok(match n {
            1..=25 => Region::Fr1,
            26..=32 => Region::Cdr1,
            33..=51 => Region::Fr2,
            52 if label.insertion == Insertion::None => Region::Fr2,
            52..=55 => Region::Cdr2,
            56..=95 => Region::Fr3,
            96..=101 => Region::Cdr3,
            102..=113 => Region::Fr4,
            _ => return Err(Error::Invalid("unsupported Chothia track position")),
        });
    }
    Ok(match n {
        1..=24 => Region::Fr1,
        25 if *locus == v_mapping::Locus::Lambda => Region::Cdr1,
        25 => Region::Fr1,
        26..=32 => Region::Cdr1,
        33..=49 => Region::Fr2,
        50..=52 => Region::Cdr2,
        53..=90 => Region::Fr3,
        91..=96 => Region::Cdr3,
        97..=107 => Region::Fr4,
        _ => return Err(Error::Invalid("unsupported Chothia track position")),
    })
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
    row.region = Some(imgt_region(&residue.position)?);
    row.assignment = Some(Assignment {
        position: residue.position.clone(),
        evidence,
    });
    Ok(())
}
fn build_track(
    bytes: &[u8],
    report: &d::Report,
    candidate: &d::Candidate,
    selection: &Selection,
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
            region: None,
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
    let (rows, issues, cdrs) =
        crate::alternate_tracks::assemble(&rows, candidate, selection, bytes, cancel)?;
    Ok(Track {
        candidate_id: candidate.id.clone(),
        locus: candidate.locus.clone(),
        numbering_scheme: selection.numbering_scheme.clone(),
        region_definition: selection.region_definition.clone(),
        invariant_position_one: framework
            .iter()
            .find(|m| m.reference_position.number == 1)
            .map(|m| m.source_index),
        first_supported_source_index: start,
        source_end: candidate.source_end,
        unassigned_residues: rows.iter().filter(|r| r.assignment.is_none()).count(),
        unassigned_regions: rows.iter().filter(|r| r.region.is_none()).count(),
        issues,
        cdrs,
        residues: rows,
    })
}
fn context_revision(
    source: &str,
    corpus: &str,
    options: &d::Options,
    limits: &d::Limits,
    selection: &Selection,
) -> Result<String> {
    hash(
        b"antibody-domain-track-context-5\0",
        &(
            TRACK_ENGINE,
            catalog_revision(),
            source,
            corpus,
            options,
            limits,
            selection,
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
    selection: &Selection,
) -> Result<String> {
    d::validate_inputs(options, limits)?;
    selection.validate_for_references(&options.references)?;
    context_revision(
        &source_revision(source)?,
        &references::load(&options.references.species)?.revision()?,
        options,
        limits,
        selection,
    )
}
/// Compute the candidate evidence once, then derive all tracks from that exact
/// immutable result. These are conditional reference/rule assignments, not a
/// calibrated domain classifier or independently certified biological numbering.
pub fn analyze(
    source: &Molecule,
    options: &d::Options,
    limits: &d::Limits,
    selection: &Selection,
    cancel: &mut impl FnMut() -> bool,
) -> Result<Report> {
    selection.validate_for_references(&options.references)?;
    let evidence = d::find_candidates(source, options, limits, cancel)?;
    let rule_revision = catalog_revision();
    let context_revision = context_revision(
        &evidence.source_revision,
        &evidence.v_evidence.corpus_revision,
        options,
        limits,
        selection,
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
            selection,
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
        selection: selection.clone(),
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
        for mapping in &e.v_evidence.mappings {
            mapping.validate_admission(&e.options.references.species)?;
        }
        self.selection
            .validate_for_references(&e.options.references)?;
        let corpus = references::load(&e.options.references.species)?.revision()?;
        if self.engine != TRACK_ENGINE
            || self.source_revision != e.source_revision
            || e.source_revision.len() != 64
            || !e.source_revision.bytes().all(|b| b.is_ascii_hexdigit())
            || self.rule_revision != catalog_revision()
            || self.context_revision != expected
            || context_revision(
                &e.source_revision,
                &corpus,
                &e.options,
                &e.limits,
                &self.selection,
            )? != expected
            || e.engine != d::CANDIDATE_ENGINE
            || e.v_evidence.engine != v_mapping::MAPPING_ENGINE
            || e.v_evidence.source_revision != e.source_revision
            || e.v_evidence.corpus_revision != corpus
            || e.v_evidence.selection != e.options.references
            || e.v_evidence.scoring != e.options.v_scoring
            || e.v_evidence.limits != e.limits.v_mapping
            || e.candidates.len() > e.limits.max_candidates
            || e.candidates.len() != self.tracks.len()
            || e.qualified_pairs_in_retained_v < e.candidates.len()
            || e.pair_retention_truncated != (e.qualified_pairs_in_retained_v > e.candidates.len())
            || e.windows_examined > e.limits.max_scan_windows
            || e.v_retention_truncated
                != e.v_evidence.mappings.iter().any(|m| m.alignment.truncated)
            || e.v_evidence.mappings.len() != e.options.references.loci.len()
            || e.v_evidence
                .mappings
                .iter()
                .zip(&e.options.references.loci)
                .any(|(mapping, locus)| &mapping.locus != locus)
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
            if !e.options.references.loci.contains(&candidate.locus)
                || track.cdrs.len() != 3
                || track.issues.len() > 6
                || track.first_supported_source_index >= track.source_end
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
            if *track != build_track(&bytes, e, candidate, &self.selection, &mut || false)? {
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
        selection: &Selection,
    ) -> Result<()> {
        if self.engine != TRACK_ENGINE
            || self.tracks.len() > d::MAX_CANDIDATES
            || self.tracks.iter().any(|t| {
                t.residues.len() > MAX_SOURCE_RESIDUES || t.cdrs.len() != 3 || t.issues.len() > 6
            })
            || *self != analyze(source, options, limits, selection, &mut || false)?
        {
            return Err(Error::Binding);
        }
        Ok(())
    }
}
