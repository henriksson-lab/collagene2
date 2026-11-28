//! Trace-to-reference evidence under explicit native-score thresholds. No basecaller,
//! calibrated Phred conversion, external aligner or silent reference-filled consensus.
use crate::{
    alignment::{self, Alignment, AlignmentMode, AlignmentOptions, Scoring},
    review::reads::traces::{self, QualityEvidence, Trace, TraceOrientation},
    sequence::analysis,
};
use collagene_model::{Alphabet, Error, Molecule, Sequence, Topology};
use serde::{Deserialize, Serialize};
use std::{
    collections::{BTreeMap, BTreeSet},
    io,
};

const MAX_READS: usize = 32;
const MAX_COLUMNS: usize = 200_000;
const MAX_OUTPUT: usize = 16 * 1024 * 1024;
const NEG: i32 = -1_000_000_000;
pub const MODEL: &str = "native-score/unique-affine-placement/unanimous-concrete-consensus-1";

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum TrimMode {
    AllCalls,
    QualityWindows,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum MissingQualityPolicy {
    Reject,
    KeepUnscored,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub struct TrimOptions {
    pub mode: TrimMode,
    pub raw_score_threshold: u8,
    pub window_size: usize,
    pub minimum_passing: usize,
    pub min_retained_bases: usize,
    pub missing_quality: MissingQualityPolicy,
}
impl Default for TrimOptions {
    fn default() -> Self {
        Self {
            mode: TrimMode::QualityWindows,
            raw_score_threshold: 20,
            window_size: 5,
            minimum_passing: 4,
            min_retained_bases: 20,
            missing_quality: MissingQualityPolicy::Reject,
        }
    }
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub struct TraceValidationOptions {
    pub trim: TrimOptions,
    pub scoring: Scoring,
    pub both_strands: bool,
    pub max_cells: usize,
    pub max_total_cells: usize,
    pub max_exact_comparisons: usize,
    pub max_consensus_sites: usize,
    pub min_support: usize,
    pub require_both_strands: bool,
}
impl Default for TraceValidationOptions {
    fn default() -> Self {
        Self {
            trim: TrimOptions::default(),
            scoring: Scoring::default(),
            both_strands: true,
            max_cells: 64_000_000,
            max_total_cells: 128_000_000,
            max_exact_comparisons: 50_000_000,
            max_consensus_sites: 100_000,
            min_support: 1,
            require_both_strands: false,
        }
    }
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum QualitySource {
    AbifRawBasecaller,
    ScfCalledBaseAccuracy,
    FastqPhred33,
    FastqPhred64,
    Missing,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ReadFormat {
    Trace,
    Fastq,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TrimReport {
    pub quality_source: QualitySource,
    pub start: usize,
    pub end: usize,
    pub original_calls: usize,
    pub passing_calls: usize,
    pub low_or_unknown_calls: Vec<usize>,
    pub actual_window_size: usize,
    pub warnings: Vec<String>,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ColumnKind {
    Pair,
    Insertion,
    Deletion,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ReadColumn {
    pub kind: ColumnKind,
    pub reference_position: Option<usize>,
    pub reference_boundary: usize,
    pub call_index: Option<usize>,
    pub stored_call_index: Option<usize>,
    pub sample_index: Option<u32>,
    pub stored_sample_index: Option<u32>,
    pub reference_base: Option<char>,
    pub observed_base: Option<char>,
    pub raw_confidence: Option<u8>,
    pub quality_passes: bool,
    pub flank_quality_passes: bool,
    pub eligible: bool,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ReadValidation {
    pub name: String,
    pub source_revision: String,
    pub source_hash: Option<String>,
    pub trace_orientation: Option<TraceOrientation>,
    pub read_format: Option<ReadFormat>,
    pub source_call_length: Option<usize>,
    pub source_retained_start: Option<usize>,
    pub source_retained_end: Option<usize>,
    pub trim: Option<TrimReport>,
    pub alignment: Option<Alignment>,
    pub reference_window_start: usize,
    pub reference_window_length: usize,
    pub query_reversed: bool,
    pub original_read_reverse_to_reference: bool,
    pub exact_placement_count: usize,
    pub optimal_paths_saturated: u8,
    pub placement_unique_proven: bool,
    pub optimality_proven: bool,
    pub search_complete: bool,
    pub paired_reference_fraction: f64,
    pub columns: Vec<ReadColumn>,
    pub warnings: Vec<String>,
    pub error: Option<String>,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum SiteState {
    SupportedMatch,
    SupportedDifference,
    ReferenceAmbiguous,
    ConflictingCalls,
    DeletionEvidence,
    InsufficientSupport,
    LowQualityOrAmbiguousOnly,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ConsensusSite {
    pub position: usize,
    pub reference_base: char,
    pub consensus_base: char,
    pub state: SiteState,
    pub counts: [usize; 4],
    pub eligible_read_indices: Vec<usize>,
    pub forward_reads: usize,
    pub reverse_reads: usize,
    pub low_quality_calls: usize,
    pub ambiguous_calls: usize,
    pub deletion_observations: usize,
    pub eligible_deletions: usize,
    pub deletion_forward_reads: usize,
    pub deletion_reverse_reads: usize,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct InsertionEvidence {
    pub reference_boundary: usize,
    pub sequence: String,
    pub read_index: usize,
    pub call_indices: Vec<usize>,
    pub original_read_reverse_to_reference: bool,
    pub eligible: bool,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ValidationStatus {
    Incomplete,
    Discrepancies,
    ConsistentUnderRecordedThresholds,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TraceValidationReport {
    pub model: String,
    pub reference_revision: String,
    pub read_revisions: Vec<String>,
    pub options: TraceValidationOptions,
    pub engine_version: String,
    pub status: ValidationStatus,
    pub reads: Vec<ReadValidation>,
    pub sites: Vec<ConsensusSite>,
    pub insertions: Vec<InsertionEvidence>,
    pub consensus: String,
    pub uncovered_ranges: Vec<(usize, usize)>,
    pub eligible_coverage: usize,
    pub reference_coverage_fraction: f64,
    pub cells_charged: usize,
    pub exact_comparisons: usize,
    pub warnings: Vec<String>,
}
fn bad(s: &str) -> Error {
    Error::Unsupported(s.into())
}
pub(crate) fn validate_options(o: &TraceValidationOptions) -> Result<(), Error> {
    if o.trim.window_size == 0
        || o.trim.window_size > 1000
        || o.trim.minimum_passing == 0
        || o.trim.minimum_passing > o.trim.window_size
        || o.trim.min_retained_bases == 0
        || o.trim.min_retained_bases > traces::MAX_TRACE_BASES
        || o.max_cells == 0
        || o.max_cells > alignment::HARD_MAX_CELLS
        || o.max_total_cells == 0
        || o.max_total_cells > 256_000_000
        || o.max_exact_comparisons == 0
        || o.max_exact_comparisons > 200_000_000
        || o.max_consensus_sites == 0
        || o.max_consensus_sites > 200_000
        || o.min_support == 0
        || o.min_support > MAX_READS
    {
        return Err(bad(
            "invalid trace validation, trimming, support or work bounds",
        ));
    }
    // Use the alignment boundary's authoritative scoring validation even for fast exact reads.
    alignment::align(
        &Sequence::dna("A")?,
        &Sequence::dna("A")?,
        &AlignmentOptions {
            global_storage: crate::alignment::grid::GlobalStorage::Full,
            mode: AlignmentMode::Global,
            scoring: o.scoring,
            max_cells: o.max_cells,
        },
    )?;
    Ok(())
}
fn native_scores(trace: &Trace) -> (QualitySource, Vec<Option<u8>>) {
    match &trace.qualities {
        Some(QualityEvidence::AbifBasecallerScores(q)) => (
            QualitySource::AbifRawBasecaller,
            q.iter().map(|&q| Some(q)).collect(),
        ),
        Some(QualityEvidence::ScfAccuracy(q)) => (
            QualitySource::ScfCalledBaseAccuracy,
            trace
                .sequence
                .as_str()
                .bytes()
                .enumerate()
                .map(|(i, b)| match b {
                    b'A' => Some(q.a[i]),
                    b'C' => Some(q.c[i]),
                    b'G' => Some(q.g[i]),
                    b'T' => Some(q.t[i]),
                    _ => None,
                })
                .collect(),
        ),
        None => (QualitySource::Missing, vec![None; trace.sequence.len()]),
    }
}
struct ReadData {
    sequence: Sequence,
    orientation: TraceOrientation,
    source_hash: String,
    format: ReadFormat,
    quality_source: QualitySource,
    scores: Vec<Option<u8>>,
    peaks: Option<Vec<u32>>,
    samples: usize,
    warnings: Vec<String>,
    source_call_length: usize,
    retained_start: usize,
    retained_end: usize,
}
impl ReadData {
    fn trace(trace: &Trace) -> Self {
        let (quality_source, scores) = native_scores(trace);
        Self {
            sequence: trace.sequence.clone(),
            orientation: trace.orientation,
            source_hash: trace.source_hash.clone(),
            format: ReadFormat::Trace,
            quality_source,
            scores,
            peaks: Some(trace.peak_positions.clone()),
            samples: trace.channels.len(),
            warnings: trace.warnings.clone(),
            source_call_length: trace.sequence.len(),
            retained_start: 0,
            retained_end: trace.sequence.len(),
        }
    }
    fn from_molecule(molecule: &Molecule) -> Result<Self, Error> {
        let is_trace = molecule.metadata.contains_key("trace_original_base64");
        let is_fastq = molecule
            .metadata
            .contains_key(crate::formats::oligos::FASTQ_ORIGINAL_KEY);
        if is_trace && is_fastq {
            return Err(bad(
                "multiple read-quality attachments require separate explicitly chosen read documents",
            ));
        }
        if is_trace {
            return Ok(Self::trace(&traces::trace_from_molecule(molecule)?));
        }
        if is_fastq {
            let evidence = crate::formats::oligos::quality_evidence(molecule)?;
            if !evidence.current_valid {
                return Err(bad(
                    "current FASTQ calls/topology no longer match immutable source quality evidence",
                ));
            }
            if evidence.called_sequence.len() > traces::MAX_TRACE_BASES {
                return Err(bad(
                    "read validation is bounded to 100,000 calls per source",
                ));
            }
            return Ok(Self {
                sequence: evidence.called_sequence,
                orientation: if evidence.reverse {
                    TraceOrientation::ReverseComplement
                } else {
                    TraceOrientation::AsStored
                },
                source_hash: evidence.original_record_blake3,
                format: ReadFormat::Fastq,
                quality_source: match evidence.encoding {
                    crate::formats::oligos::QualityEncoding::Phred33 => QualitySource::FastqPhred33,
                    crate::formats::oligos::QualityEncoding::Phred64 => QualitySource::FastqPhred64,
                },
                scores: evidence.qualities.into_iter().map(Some).collect(),
                peaks: None,
                samples: 0,
                warnings: evidence.warnings,
                source_call_length: evidence.source_call_length,
                retained_start: evidence.retained_start,
                retained_end: evidence.retained_end,
            });
        }
        Err(bad(
            "read has no supported immutable trace or explicit-encoding FASTQ quality attachment",
        ))
    }
}
/// The threshold addresses each source's native raw scale, never an inferred error probability.
pub fn trim_trace(trace: &Trace, options: &TrimOptions) -> Result<TrimReport, Error> {
    trace.validate()?;
    trim_data(&ReadData::trace(trace), options)
}
pub fn trim_read(molecule: &Molecule, options: &TrimOptions) -> Result<TrimReport, Error> {
    trim_data(&ReadData::from_molecule(molecule)?, options)
}
fn trim_data(trace: &ReadData, options: &TrimOptions) -> Result<TrimReport, Error> {
    if options.window_size == 0
        || options.window_size > 1000
        || options.minimum_passing == 0
        || options.minimum_passing > options.window_size
        || options.min_retained_bases == 0
        || options.min_retained_bases > traces::MAX_TRACE_BASES
    {
        return Err(bad("invalid trace trimming options"));
    }
    let (quality_source, scores) = (trace.quality_source, &trace.scores);
    let n = trace.sequence.len();
    if quality_source == QualitySource::Missing
        && options.missing_quality == MissingQualityPolicy::Reject
    {
        return Err(bad(
            "trace has no native confidence evidence; explicitly allow unscored alignment or use another read",
        ));
    }
    let pass: Vec<_> = scores
        .iter()
        .zip(trace.sequence.as_str().bytes())
        .map(|(q, b)| {
            q.is_some_and(|q| q >= options.raw_score_threshold)
                && analysis::base_mask(b).count_ones() == 1
        })
        .collect();
    let mut warnings = vec![if matches!(
        quality_source,
        QualitySource::FastqPhred33 | QualitySource::FastqPhred64
    ) {
        "FASTQ scores use the explicitly declared Phred ASCII offset; caller calibration is not independently verified and no chromatogram is inferred.".into()
    } else {
        "Raw confidence thresholds are source-scale policies, not calibrated Phred probabilities or validated experimental accuracy.".into()
    }];
    let window = options.window_size.min(n);
    let required = options.minimum_passing.min(window);
    let (mut start, mut end) = (0, n);
    if quality_source == QualitySource::Missing {
        warnings.push("Unscored calls retained for alignment only; they cannot contribute eligible consensus bases.".into());
    } else if options.mode == TrimMode::QualityWindows {
        let mut count = pass[..window].iter().filter(|&&p| p).count();
        let mut first = None;
        let mut last = 0;
        for at in 0..=n - window {
            if at > 0 {
                count -= usize::from(pass[at - 1]);
                count += usize::from(pass[at + window - 1]);
            }
            if count >= required {
                first.get_or_insert(at);
                last = at + window;
            }
        }
        (start, end) = first.map(|first| (first, last)).unwrap_or((0, 0));
    }
    if end - start < options.min_retained_bases {
        start = 0;
        end = 0;
        warnings.push(
            "No retained interval meets the requested minimum length and trimming policy.".into(),
        );
    }
    Ok(TrimReport {
        quality_source,
        start,
        end,
        original_calls: n,
        passing_calls: pass[start..end].iter().filter(|&&p| p).count(),
        low_or_unknown_calls: (start..end).filter(|&i| !pass[i]).collect(),
        actual_window_size: window,
        warnings,
    })
}
#[derive(Clone, Copy)]
struct Score {
    score: i32,
    count: u8,
}
const EMPTY: Score = Score {
    score: NEG,
    count: 0,
};
fn maximum(values: impl IntoIterator<Item = Score>) -> Score {
    let mut best = EMPTY;
    for value in values {
        if value.count == 0 {
            continue;
        }
        if value.score > best.score {
            best = value;
        } else if value.score == best.score {
            best.count = best.count.saturating_add(value.count).min(2);
        }
    }
    best
}
fn plus(v: Score, delta: i32) -> Score {
    if v.count == 0 {
        EMPTY
    } else {
        Score {
            score: v.score + delta,
            ..v
        }
    }
}
fn substitution(a: u8, b: u8, s: Scoring) -> i32 {
    let (a, b) = (analysis::base_mask(a), analysis::base_mask(b));
    if a & b == 0 {
        s.mismatch_score
    } else if a.count_ones() == 1 && b.count_ones() == 1 {
        s.match_score
    } else {
        s.ambiguous_score
    }
}
/// Independent rolling optimal-path counter. Saturation at two proves nonuniqueness;
/// circular free starts are restricted to canonical first-turn boundaries.
fn optimal_paths(
    reference: &Sequence,
    query: &Sequence,
    canonical_starts: usize,
    s: Scoring,
) -> Score {
    let m = query.len();
    let mut previous = vec![[EMPTY; 3]; m + 1];
    let mut current = previous.clone();
    previous[0][0] = Score { score: 0, count: 1 };
    for j in 1..=m {
        previous[j][1] = Score {
            score: s.gap_open + (j as i32 - 1) * s.gap_extend,
            count: 1,
        };
    }
    let mut endpoint = maximum(previous[m]);
    for (i, &r) in reference.as_str().as_bytes().iter().enumerate() {
        current.fill([EMPTY; 3]);
        if i + 1 < canonical_starts {
            current[0][0] = Score { score: 0, count: 1 };
        }
        for (j, &q) in query.as_str().as_bytes().iter().enumerate() {
            let j = j + 1;
            current[j][0] = plus(maximum(previous[j - 1]), substitution(r, q, s));
            current[j][1] = maximum([
                plus(current[j - 1][0], s.gap_open),
                plus(current[j - 1][1], s.gap_extend),
                plus(current[j - 1][2], s.gap_open),
            ]);
            current[j][2] = maximum([
                plus(previous[j][0], s.gap_open),
                plus(previous[j][1], s.gap_open),
                plus(previous[j][2], s.gap_extend),
            ]);
        }
        endpoint = maximum([endpoint, maximum(current[m])]);
        std::mem::swap(&mut previous, &mut current);
    }
    endpoint
}
struct Placement {
    alignment: Alignment,
    window_start: usize,
    window_length: usize,
    reversed: bool,
    count: u8,
    unique: bool,
    complete: bool,
    optimal: bool,
    exact_count: usize,
    warnings: Vec<String>,
}
fn alignment_for_read(
    reference: &Molecule,
    query: &Sequence,
    o: &TraceValidationOptions,
    cells: &mut usize,
    comparisons: &mut usize,
) -> Result<Placement, Error> {
    let n = reference.sequence.len();
    let m = query.len();
    if m > n && reference.topology == Topology::Circular {
        return Err(bad(
            "a trace cannot be mapped through more than one circular reference turn",
        ));
    }
    let directions = if o.both_strands { 2 } else { 1 };
    let mut orientations = vec![query.clone()];
    if directions == 2 {
        orientations.push(query.reverse_complement()?);
    }
    let mut exact = Vec::new();
    let mut scan_complete = true;
    // The alignment module's identity materializer is linear only when ambiguity
    // is nonnegative. Otherwise use the fully charged DP/proof path below.
    if o.scoring.ambiguous_score >= 0
        && o.scoring.ambiguous_score < o.scoring.match_score
        && query
            .as_str()
            .bytes()
            .all(|b| analysis::base_mask(b).count_ones() == 1)
    {
        'scan: for (reverse, q) in orientations.iter().enumerate() {
            let starts = if reference.topology == Topology::Circular {
                n
            } else {
                n.saturating_sub(m) + usize::from(m <= n)
            };
            for start in 0..starts {
                let mut matches = true;
                for (j, &b) in q.as_str().as_bytes().iter().enumerate() {
                    if *comparisons >= o.max_exact_comparisons {
                        scan_complete = false;
                        break 'scan;
                    }
                    *comparisons += 1;
                    if reference.sequence.as_str().as_bytes()[(start + j) % n] != b {
                        matches = false;
                        break;
                    }
                }
                if matches {
                    exact.push((start, reverse == 1));
                    if exact.len() >= 64 {
                        scan_complete = false;
                        break 'scan;
                    }
                }
            }
        }
        if let Some(&(start, reversed)) = exact.first() {
            let selected = &orientations[usize::from(reversed)];
            let a = alignment::align(
                selected,
                selected,
                &AlignmentOptions {
                    global_storage: crate::alignment::grid::GlobalStorage::Full,
                    mode: AlignmentMode::Global,
                    scoring: o.scoring,
                    max_cells: o.max_cells,
                },
            )?;
            return Ok(Placement{alignment:a,window_start:start,window_length:m,reversed,count:exact.len().min(2) as u8,unique:exact.len()==1&&scan_complete,complete:scan_complete,optimal:true,exact_count:exact.len(),warnings:vec!["Exact matching attains the scoring upper bound; canonical source placements, including requested read orientations, were scanned for placement uniqueness.".into()]});
        }
    }
    let extended = if reference.topology == Topology::Circular {
        Sequence::dna(&format!(
            "{}{}",
            reference.sequence.as_str(),
            reference.sequence.as_str()
        ))?
    } else {
        reference.sequence.clone()
    };
    let cost = (extended.len() + 1)
        .checked_mul(m + 1)
        .and_then(|v| v.checked_mul(2))
        .ok_or(Error::TooLarge)?;
    let mut best: Option<Placement> = None;
    let mut complete = true;
    let mut global = EMPTY;
    for (reverse, q) in orientations.iter().enumerate() {
        if cost / 2 > o.max_cells
            || cells
                .checked_add(cost)
                .is_none_or(|v| v > o.max_total_cells)
        {
            complete = false;
            continue;
        }
        let a = alignment::align(
            &extended,
            q,
            &AlignmentOptions {
                global_storage: crate::alignment::grid::GlobalStorage::Full,
                mode: AlignmentMode::QueryInReference,
                scoring: o.scoring,
                max_cells: o.max_cells,
            },
        )?;
        *cells += cost;
        let counted = optimal_paths(
            &extended,
            q,
            if reference.topology == Topology::Circular {
                n
            } else {
                n + 1
            },
            o.scoring,
        );
        let valid = reference.topology == Topology::Linear
            || (a.reference_start < n && a.reference_end - a.reference_start <= n);
        if !valid {
            return Err(bad(
                "best circular trace alignment exceeds one reference turn; a constrained alternative was not established",
            ));
        }
        if a.score != counted.score {
            return Err(bad(
                "alignment and canonical-placement proof disagree; no unique mapping is asserted",
            ));
        }
        global = maximum([global, counted]);
        if best.as_ref().is_none_or(|b| a.score > b.alignment.score) {
            best=Some(Placement{alignment:a,window_start:0,window_length:extended.len(),reversed:reverse==1,count:counted.count,unique:false,complete:false,optimal:false,exact_count:0,warnings:vec!["Affine end-free alignment and an independent optimal-path count establish placement only within the requested orientation and one-turn reference model; tied paths do not contribute consensus.".into()]});
        }
    }
    let mut best=best.ok_or_else(||bad("trace alignment/proof exceeds the remaining DP-cell budget; reduce reads or raise the recorded bounded limit"))?;
    best.count = global.count;
    best.unique = complete && global.count == 1;
    best.complete = complete;
    best.optimal = complete;
    if !complete {
        best.warnings.push("An orientation was skipped by the DP budget; alignment optimality and placement uniqueness remain unproven.".into());
    }
    Ok(best)
}
fn empty_read(m: &Molecule) -> ReadValidation {
    ReadValidation {
        name: m.name.clone(),
        source_revision: m.revision(),
        source_hash: None,
        trace_orientation: None,
        read_format: None,
        source_call_length: None,
        source_retained_start: None,
        source_retained_end: None,
        trim: None,
        alignment: None,
        reference_window_start: 0,
        reference_window_length: 0,
        query_reversed: false,
        original_read_reverse_to_reference: false,
        exact_placement_count: 0,
        optimal_paths_saturated: 0,
        placement_unique_proven: false,
        optimality_proven: false,
        search_complete: false,
        paired_reference_fraction: 0.,
        columns: vec![],
        warnings: vec![],
        error: None,
    }
}
fn call_index(trim: &TrimReport, q: usize, reverse: bool) -> usize {
    if reverse {
        trim.end - 1 - q
    } else {
        trim.start + q
    }
}
fn qpasses(
    trace: &ReadData,
    scores: &[Option<u8>],
    call: usize,
    o: &TraceValidationOptions,
) -> bool {
    scores[call].is_some_and(|v| v >= o.trim.raw_score_threshold)
        && analysis::base_mask(trace.sequence.as_str().as_bytes()[call]).count_ones() == 1
}
fn columns(
    reference: &Molecule,
    trace: &ReadData,
    trim: &TrimReport,
    placement: &Placement,
    o: &TraceValidationOptions,
) -> Result<Vec<ReadColumn>, Error> {
    let a = &placement.alignment;
    let mut r = a.reference_start;
    let mut q = a.query_start;
    let scores = &trace.scores;
    let mut result = Vec::new();
    let n = reference.sequence.len();
    if a.reference_aligned.len() > MAX_COLUMNS {
        return Err(Error::TooLarge);
    }
    for (&rb, &qb) in a
        .reference_aligned
        .as_bytes()
        .iter()
        .zip(a.query_aligned.as_bytes())
    {
        let reference_position =
            (rb != b'-').then_some(if reference.topology == Topology::Circular {
                (placement.window_start + r) % n
            } else {
                placement.window_start + r
            });
        let boundary = if reference.topology == Topology::Circular {
            (placement.window_start + r) % n
        } else {
            placement.window_start + r
        };
        let call = (qb != b'-').then(|| call_index(trim, q, placement.reversed));
        let stored_call = call.map(|i| {
            if trace.orientation == TraceOrientation::AsStored {
                trace.retained_start + i
            } else {
                trace.retained_end - 1 - i
            }
        });
        let sample = call.and_then(|i| trace.peaks.as_ref().map(|peaks| peaks[i]));
        let stored_sample = sample.map(|p| {
            if trace.orientation == TraceOrientation::AsStored {
                p
            } else {
                (trace.samples - 1) as u32 - p
            }
        });
        let quality = call.is_some_and(|i| qpasses(trace, &scores, i, o));
        let flank = if qb == b'-' && q > 0 && q < trim.end - trim.start {
            [q - 1, q]
                .into_iter()
                .all(|j| qpasses(trace, &scores, call_index(trim, j, placement.reversed), o))
        } else {
            false
        };
        let kind = if rb == b'-' {
            ColumnKind::Insertion
        } else if qb == b'-' {
            ColumnKind::Deletion
        } else {
            ColumnKind::Pair
        };
        result.push(ReadColumn {
            kind,
            reference_position,
            reference_boundary: boundary,
            call_index: call,
            stored_call_index: stored_call,
            sample_index: sample,
            stored_sample_index: stored_sample,
            reference_base: (rb != b'-').then_some(rb as char),
            observed_base: (qb != b'-').then_some(qb as char),
            raw_confidence: call.and_then(|i| scores[i]),
            quality_passes: quality,
            flank_quality_passes: flank,
            eligible: placement.unique
                && if kind == ColumnKind::Deletion {
                    flank
                } else {
                    quality
                },
        });
        r += usize::from(rb != b'-');
        q += usize::from(qb != b'-');
    }
    Ok(result)
}
#[derive(Default)]
struct SiteAccumulator {
    counts: [usize; 4],
    eligible: BTreeSet<usize>,
    forward: usize,
    reverse: usize,
    low: usize,
    ambiguous: usize,
    deletions: usize,
    eligible_deletions: usize,
    deletion_forward: usize,
    deletion_reverse: usize,
}
fn iupac(mask: u8) -> u8 {
    match mask {
        1 => b'A',
        2 => b'C',
        4 => b'G',
        8 => b'T',
        3 => b'M',
        5 => b'R',
        6 => b'S',
        7 => b'V',
        9 => b'W',
        10 => b'Y',
        11 => b'H',
        12 => b'K',
        13 => b'D',
        14 => b'B',
        _ => b'N',
    }
}
struct BoundedBytes(usize, usize);
impl io::Write for BoundedBytes {
    fn write(&mut self, b: &[u8]) -> io::Result<usize> {
        self.0 = self
            .0
            .checked_add(b.len())
            .filter(|&n| n <= self.1)
            .ok_or_else(|| io::Error::other("trace evidence output limit"))?;
        Ok(b.len())
    }
    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}
/// Revalidates immutable trace attachments, then reconstructs read evidence against a
/// separately bound reference. Missing and low-quality regions are never filled from it.
pub fn validate_reads(
    reference: &Molecule,
    reads: &[Molecule],
    options: &TraceValidationOptions,
) -> Result<TraceValidationReport, Error> {
    reference.validate()?;
    validate_options(options)?;
    if reference.sequence.alphabet() != Alphabet::Dna {
        return Err(bad("trace validation requires a DNA reference"));
    }
    if reads.is_empty() || reads.len() > MAX_READS {
        return Err(bad("trace validation needs 1–32 immutable read documents"));
    }
    let mut input = BoundedBytes(0, 32 * 1024 * 1024);
    serde_json::to_writer(&mut input, reference).map_err(|_| Error::TooLarge)?;
    for read in reads {
        read.validate()?;
        serde_json::to_writer(&mut input, read).map_err(|_| Error::TooLarge)?;
    }
    let mut result=TraceValidationReport{model:MODEL.into(),reference_revision:reference.revision(),read_revisions:reads.iter().map(Molecule::revision).collect(),options:options.clone(),engine_version:crate::verification::ENGINE_VERSION.into(),status:ValidationStatus::Incomplete,reads:vec![],sites:vec![],insertions:vec![],consensus:String::new(),uncovered_ranges:vec![],eligible_coverage:0,reference_coverage_fraction:0.,cells_charged:0,exact_comparisons:0,warnings:vec!["SCF accuracy and ABIF caller scores remain in native raw scales; FASTQ uses its explicitly selected Phred encoding. Caller calibration and sample identity are not independently confirmed.".into(),"Consensus uses only uniquely placed concrete calls that pass the recorded native-score threshold. Uncovered, ambiguous, conflicting and insufficiently supported sites remain N or explicit IUPAC conflicts; reference bases are never used to fill gaps.".into(),"Distinct original-file hashes identify evidence files, not independent molecules or biological samples. Paired alignment coverage and eligible consensus coverage are different quantities.".into()]};
    let mut seen = BTreeSet::new();
    let mut sites: BTreeMap<usize, SiteAccumulator> = BTreeMap::new();
    let mut total_columns = 0;
    for (read_index, molecule) in reads.iter().enumerate() {
        let mut read = empty_read(molecule);
        let process = (|| -> Result<(), Error> {
            let trace = ReadData::from_molecule(molecule)?;
            read.source_hash = Some(trace.source_hash.clone());
            read.trace_orientation =
                (trace.format == ReadFormat::Trace).then_some(trace.orientation);
            read.read_format = Some(trace.format);
            read.source_call_length = Some(trace.source_call_length);
            read.source_retained_start = Some(trace.retained_start);
            read.source_retained_end = Some(trace.retained_end);
            read.warnings.extend(trace.warnings.clone());
            if !seen.insert(trace.source_hash.clone()) {
                return Err(bad(
                    "duplicate original trace attachment is not independent read support",
                ));
            }
            let trim = trim_data(&trace, &options.trim)?;
            read.warnings.extend(trim.warnings.clone());
            read.trim = Some(trim.clone());
            if trim.end <= trim.start {
                return Err(bad(
                    "no trace calls retained by the requested trimming policy",
                ));
            }
            let query = Sequence::dna(&trace.sequence.as_str()[trim.start..trim.end])?;
            let placement = alignment_for_read(
                reference,
                &query,
                options,
                &mut result.cells_charged,
                &mut result.exact_comparisons,
            )?;
            read.reference_window_start = placement.window_start;
            read.reference_window_length = placement.window_length;
            read.query_reversed = placement.reversed;
            read.original_read_reverse_to_reference =
                placement.reversed ^ (trace.orientation == TraceOrientation::ReverseComplement);
            read.exact_placement_count = placement.exact_count;
            read.optimal_paths_saturated = placement.count;
            read.placement_unique_proven = placement.unique;
            read.optimality_proven = placement.optimal;
            read.search_complete = placement.complete;
            read.warnings.extend(placement.warnings.clone());
            read.columns = columns(reference, &trace, &trim, &placement, options)?;
            total_columns += read.columns.len();
            if total_columns > MAX_COLUMNS {
                return Err(Error::TooLarge);
            }
            read.paired_reference_fraction = read
                .columns
                .iter()
                .filter_map(|c| {
                    (c.kind == ColumnKind::Pair)
                        .then_some(c.reference_position)
                        .flatten()
                })
                .collect::<BTreeSet<_>>()
                .len() as f64
                / reference.sequence.len() as f64;
            if !placement.unique {
                read.warnings.push("Read placement or optimal alignment is tied/unproven; displayed coordinates are one preview and do not contribute consensus evidence.".into());
            }
            read.alignment = Some(placement.alignment);
            Ok(())
        })();
        if let Err(error) = process {
            read.error = Some(error.to_string());
            read.columns.clear();
        }
        let mut insertion: Option<InsertionEvidence> = None;
        for column in &read.columns {
            if column.kind == ColumnKind::Insertion {
                let call = column.call_index.expect("insertion call");
                if insertion
                    .as_ref()
                    .is_none_or(|i| i.reference_boundary != column.reference_boundary)
                {
                    if let Some(i) = insertion.take() {
                        result.insertions.push(i);
                    }
                    insertion = Some(InsertionEvidence {
                        reference_boundary: column.reference_boundary,
                        sequence: String::new(),
                        read_index,
                        call_indices: vec![],
                        original_read_reverse_to_reference: read.original_read_reverse_to_reference,
                        eligible: true,
                    });
                }
                let current = insertion.as_mut().unwrap();
                current.sequence.push(column.observed_base.unwrap());
                current.call_indices.push(call);
                current.eligible &= column.eligible;
                continue;
            }
            if let Some(i) = insertion.take() {
                result.insertions.push(i);
            }
            let Some(position) = column.reference_position else {
                continue;
            };
            if !sites.contains_key(&position) && sites.len() >= options.max_consensus_sites {
                return Err(Error::TooLarge);
            }
            let site = sites.entry(position).or_default();
            if column.kind == ColumnKind::Deletion {
                site.deletions += 1;
                site.eligible_deletions += usize::from(column.eligible);
                if column.eligible {
                    if read.original_read_reverse_to_reference {
                        site.deletion_reverse += 1
                    } else {
                        site.deletion_forward += 1
                    }
                }
                continue;
            }
            let base = column.observed_base.unwrap() as u8;
            let mask = analysis::base_mask(base);
            site.low += usize::from(!column.quality_passes);
            site.ambiguous += usize::from(mask.count_ones() != 1);
            if column.eligible && mask.count_ones() == 1 {
                site.counts[mask.trailing_zeros() as usize] += 1;
                if site.eligible.insert(read_index) {
                    if read.original_read_reverse_to_reference {
                        site.reverse += 1
                    } else {
                        site.forward += 1
                    }
                }
            }
        }
        if let Some(i) = insertion {
            result.insertions.push(i);
        }
        result.reads.push(read);
    }
    let mut insertion_groups: BTreeMap<(usize, String), BTreeSet<usize>> = BTreeMap::new();
    for insertion in result.insertions.iter().filter(|i| i.eligible) {
        insertion_groups
            .entry((insertion.reference_boundary, insertion.sequence.clone()))
            .or_default()
            .insert(insertion.read_index);
    }
    let mut discrepancy = insertion_groups.values().any(|reads| {
        reads.len() >= options.min_support
            && (!options.require_both_strands
                || (reads
                    .iter()
                    .any(|&i| !result.reads[i].original_read_reverse_to_reference)
                    && reads
                        .iter()
                        .any(|&i| result.reads[i].original_read_reverse_to_reference)))
    });
    let mut consensus = vec![b'N'; reference.sequence.len()];
    let mut covered = BTreeSet::new();
    for (position, site) in sites {
        let rb = reference.sequence.as_str().as_bytes()[position];
        let mask = site
            .counts
            .iter()
            .enumerate()
            .fold(0u8, |mask, (i, &n)| mask | if n > 0 { 1 << i } else { 0 });
        let support = site.eligible.len() >= options.min_support
            && (!options.require_both_strands || site.forward > 0 && site.reverse > 0);
        let state = if mask.count_ones() > 1 && support {
            consensus[position] = iupac(mask);
            discrepancy = true;
            SiteState::ConflictingCalls
        } else if site.eligible_deletions >= options.min_support
            && (!options.require_both_strands
                || site.deletion_forward > 0 && site.deletion_reverse > 0)
        {
            discrepancy = true;
            SiteState::DeletionEvidence
        } else if support && mask.count_ones() == 1 {
            covered.insert(position);
            consensus[position] = iupac(mask);
            if analysis::base_mask(rb) & mask == 0 {
                discrepancy = true;
                SiteState::SupportedDifference
            } else if analysis::base_mask(rb).count_ones() != 1 {
                SiteState::ReferenceAmbiguous
            } else {
                SiteState::SupportedMatch
            }
        } else if mask != 0 {
            SiteState::InsufficientSupport
        } else {
            SiteState::LowQualityOrAmbiguousOnly
        };
        result.sites.push(ConsensusSite {
            position,
            reference_base: rb as char,
            consensus_base: consensus[position] as char,
            state,
            counts: site.counts,
            eligible_read_indices: site.eligible.into_iter().collect(),
            forward_reads: site.forward,
            reverse_reads: site.reverse,
            low_quality_calls: site.low,
            ambiguous_calls: site.ambiguous,
            deletion_observations: site.deletions,
            eligible_deletions: site.eligible_deletions,
            deletion_forward_reads: site.deletion_forward,
            deletion_reverse_reads: site.deletion_reverse,
        });
    }
    let mut at = 0;
    for &position in &covered {
        if at < position {
            result.uncovered_ranges.push((at, position));
        }
        at = position + 1;
    }
    if at < reference.sequence.len() {
        result.uncovered_ranges.push((at, reference.sequence.len()));
    }
    result.eligible_coverage = covered.len();
    result.reference_coverage_fraction = covered.len() as f64 / reference.sequence.len() as f64;
    result.consensus = String::from_utf8(consensus).expect("IUPAC ASCII");
    result.status = if discrepancy {
        ValidationStatus::Discrepancies
    } else if covered.len() == reference.sequence.len()
        && result
            .sites
            .iter()
            .all(|s| s.state == SiteState::SupportedMatch)
        && result
            .reads
            .iter()
            .all(|r| r.error.is_none() && r.placement_unique_proven && r.search_complete)
    {
        ValidationStatus::ConsistentUnderRecordedThresholds
    } else {
        ValidationStatus::Incomplete
    };
    let mut size = BoundedBytes(0, MAX_OUTPUT);
    serde_json::to_writer(&mut size, &result).map_err(|_| Error::TooLarge)?;
    Ok(result)
}
