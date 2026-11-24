//! Independently implemented affine-gap nucleotide pairwise alignment.
//! All coordinates are zero-based and half-open on the input reference/query.
//! This module invokes no external alignment engine and does not infer read quality.
use crate::{Alphabet, Error, Molecule, Sequence, Topology, analysis};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

pub const DEFAULT_MAX_CELLS: usize = 16_000_000;
pub const HARD_MAX_CELLS: usize = 64_000_000;
pub const MAX_ALIGNMENT_COLUMNS: usize = 1_000_000;
const MAX_MATRIX_BYTES: usize = 96 * 1024 * 1024;
const NEG: i32 = -1_000_000_000;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum AlignmentMode {
    Global,
    Local,
    QueryInReference,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub struct Scoring {
    pub match_score: i32,
    /// Compatible IUPAC ambiguity is scored separately and never counted as a certain match.
    pub ambiguous_score: i32,
    pub mismatch_score: i32,
    /// A gap of length k costs gap_open + (k-1)*gap_extend.
    pub gap_open: i32,
    pub gap_extend: i32,
}
impl Default for Scoring {
    fn default() -> Self {
        Self {
            match_score: 2,
            ambiguous_score: 0,
            mismatch_score: -3,
            gap_open: -4,
            gap_extend: -1,
        }
    }
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub struct AlignmentOptions {
    pub mode: AlignmentMode,
    pub scoring: Scoring,
    pub max_cells: usize,
}
impl Default for AlignmentOptions {
    fn default() -> Self {
        Self {
            mode: AlignmentMode::Global,
            scoring: Scoring::default(),
            max_cells: DEFAULT_MAX_CELLS,
        }
    }
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum CigarKind {
    Match,
    Ambiguous,
    Mismatch,
    Insertion,
    Deletion,
}
impl CigarKind {
    fn symbol(self) -> char {
        match self {
            Self::Match => '=',
            Self::Ambiguous => 'M',
            Self::Mismatch => 'X',
            Self::Insertion => 'I',
            Self::Deletion => 'D',
        }
    }
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct CigarRun {
    pub kind: CigarKind,
    pub length: usize,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum DifferenceKind {
    AmbiguousCompatible,
    Substitution,
    Insertion,
    Deletion,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Difference {
    pub kind: DifferenceKind,
    pub reference_start: usize,
    pub reference_end: usize,
    pub query_start: usize,
    pub query_end: usize,
    pub reference_bases: String,
    pub query_bases: String,
}
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct AlignmentStats {
    pub columns: usize,
    pub paired_bases: usize,
    pub certain_matches: usize,
    pub ambiguous_compatible: usize,
    pub mismatches: usize,
    pub inserted_query_bases: usize,
    pub deleted_reference_bases: usize,
    pub reference_span_bases: usize,
    pub query_span_bases: usize,
    /// Paired columns divided by the full input length, not merely the aligned span.
    pub reference_coverage: f64,
    pub query_coverage: f64,
    /// Denominator is paired columns; None when there are no paired columns.
    pub certain_identity: Option<f64>,
    pub compatible_identity: Option<f64>,
}
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Alignment {
    pub algorithm: String,
    pub options: AlignmentOptions,
    pub score: i32,
    pub reference_start: usize,
    pub reference_end: usize,
    pub query_start: usize,
    pub query_end: usize,
    pub reference_aligned: String,
    pub query_aligned: String,
    /// Extended CIGAR: '=' certain, 'M' compatible ambiguity, 'X' disjoint, I/D gaps.
    pub cigar: String,
    pub operations: Vec<CigarRun>,
    pub differences: Vec<Difference>,
    pub stats: AlignmentStats,
    pub warnings: Vec<String>,
    pub reference_hash: String,
    pub query_hash: String,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(u8)]
enum State {
    Pair = 0,
    Insert = 1,
    Delete = 2,
    Stop = 3,
}
#[derive(Debug)]
struct Path {
    score: i32,
    rs: usize,
    re: usize,
    qs: usize,
    qe: usize,
    ops: Vec<State>,
}
fn invalid(message: &str) -> Error {
    Error::Format(message.into())
}
fn unsupported(message: &str) -> Error {
    Error::Unsupported(message.into())
}
fn validate(
    reference: &Sequence,
    query: &Sequence,
    options: &AlignmentOptions,
) -> Result<(), Error> {
    if reference.alphabet() == Alphabet::Protein || reference.alphabet() != query.alphabet() {
        return Err(unsupported(
            "pairwise alignment currently requires matching DNA or RNA alphabets",
        ));
    }
    let s = options.scoring;
    if !(1..=10_000).contains(&s.match_score)
        || !(-10_000..=0).contains(&s.mismatch_score)
        || s.ambiguous_score < s.mismatch_score
        || s.ambiguous_score > s.match_score
        || !(-10_000..=-1).contains(&s.gap_open)
        || !(-10_000..=-1).contains(&s.gap_extend)
    {
        return Err(invalid(
            "invalid pairwise scoring: positive bounded match, mismatch <= ambiguity <= match, and strictly negative bounded affine gap scores are required",
        ));
    }
    if options.max_cells == 0 || options.max_cells > HARD_MAX_CELLS {
        return Err(invalid(
            "max_cells must be between 1 and the 64 million cell hard limit",
        ));
    }
    let bound = (reference.len() as u64 + query.len() as u64)
        * [
            s.match_score,
            s.ambiguous_score,
            s.mismatch_score,
            s.gap_open,
            s.gap_extend,
        ]
        .into_iter()
        .map(i32::unsigned_abs)
        .max()
        .unwrap() as u64;
    if bound > 500_000_000 {
        return Err(unsupported(
            "input lengths and scoring exceed the safe pairwise score bound",
        ));
    }
    Ok(())
}
fn kind(a: u8, b: u8) -> CigarKind {
    let a = analysis::base_mask(a);
    let b = analysis::base_mask(b);
    if a & b == 0 {
        CigarKind::Mismatch
    } else if a.count_ones() == 1 && b.count_ones() == 1 {
        CigarKind::Match
    } else {
        CigarKind::Ambiguous
    }
}
fn substitution(a: u8, b: u8, s: Scoring) -> i32 {
    match kind(a, b) {
        CigarKind::Match => s.match_score,
        CigarKind::Ambiguous => s.ambiguous_score,
        _ => s.mismatch_score,
    }
}
fn add(score: i32, delta: i32) -> i32 {
    if score == NEG { NEG } else { score + delta }
}
fn best(values: &[(i32, State)]) -> (i32, State) {
    let mut best = values[0];
    for &candidate in &values[1..] {
        if candidate.0 > best.0 {
            best = candidate;
        }
    }
    best
}
fn set(trace: &mut [u8], index: usize, state: State, previous: State) {
    let shift = (state as u8) * 2;
    trace[index] = (trace[index] & !(3 << shift)) | ((previous as u8) << shift);
}
fn predecessor(trace: &[u8], index: usize, state: State) -> State {
    match (trace[index] >> ((state as u8) * 2)) & 3 {
        0 => State::Pair,
        1 => State::Insert,
        2 => State::Delete,
        _ => State::Stop,
    }
}
fn cells(n: usize, m: usize) -> Result<usize, Error> {
    n.checked_add(1)
        .and_then(|n| m.checked_add(1).and_then(|m| n.checked_mul(m)))
        .ok_or(Error::TooLarge)
}
fn checked_layout(n: usize, m: usize, options: &AlignmentOptions) -> Result<usize, Error> {
    let count = cells(n, m)?;
    if count > options.max_cells {
        return Err(unsupported(
            "pairwise alignment exceeds its requested matrix-cell budget; shorten the sequences or explicitly raise the bounded budget",
        ));
    }
    let bytes = (m + 1)
        .checked_mul(2 * 3 * std::mem::size_of::<i32>())
        .and_then(|rows| rows.checked_add(count))
        .ok_or(Error::TooLarge)?;
    if bytes > MAX_MATRIX_BYTES {
        return Err(unsupported(
            "pairwise matrix/traceback exceeds the 96 MiB allocation budget",
        ));
    }
    Ok(count)
}
fn compute(
    reference: &Sequence,
    query: &Sequence,
    options: &AlignmentOptions,
) -> Result<Path, Error> {
    let n = reference.len();
    let m = query.len();
    let count = checked_layout(n, m, options)?;
    let width = m + 1;
    let mut trace = Vec::new();
    trace
        .try_reserve_exact(count)
        .map_err(|_| Error::TooLarge)?;
    trace.resize(count, 0x3f);
    let mut previous = Vec::new();
    previous
        .try_reserve_exact(width)
        .map_err(|_| Error::TooLarge)?;
    previous.resize(width, [NEG; 3]);
    let mut current = previous.clone();
    let s = options.scoring;
    previous[0][0] = 0;
    for j in 1..=m {
        if options.mode == AlignmentMode::Local {
            previous[j][0] = 0;
        } else {
            previous[j][1] = s.gap_open + (j as i32 - 1) * s.gap_extend;
            set(
                &mut trace,
                j,
                State::Insert,
                if j == 1 { State::Pair } else { State::Insert },
            );
        }
    }
    let mut endpoint = if options.mode == AlignmentMode::QueryInReference {
        let value = best(&[
            (previous[m][0], State::Pair),
            (previous[m][1], State::Insert),
            (previous[m][2], State::Delete),
        ]);
        (value.0, 0, m, value.1)
    } else {
        (0, 0, 0, State::Stop)
    };
    let rb = reference.as_str().as_bytes();
    let qb = query.as_str().as_bytes();
    for i in 1..=n {
        current.fill([NEG; 3]);
        if options.mode == AlignmentMode::Global {
            current[0][2] = s.gap_open + (i as i32 - 1) * s.gap_extend;
            set(
                &mut trace,
                i * width,
                State::Delete,
                if i == 1 { State::Pair } else { State::Delete },
            );
        } else {
            current[0][0] = 0;
        }
        for j in 1..=m {
            let index = i * width + j;
            let mut pair = best(&[
                (previous[j - 1][0], State::Pair),
                (previous[j - 1][1], State::Insert),
                (previous[j - 1][2], State::Delete),
            ]);
            pair.0 = add(pair.0, substitution(rb[i - 1], qb[j - 1], s));
            if options.mode == AlignmentMode::Local && pair.0 <= 0 {
                pair = (0, State::Stop);
            }
            // Opposite-gap transitions are allowed, charging a new opening; this covers all
            // alignments even when a user supplies very costly mismatch scores.
            let insert = best(&[
                (add(current[j - 1][1], s.gap_extend), State::Insert),
                (add(current[j - 1][0], s.gap_open), State::Pair),
                (add(current[j - 1][2], s.gap_open), State::Delete),
            ]);
            let delete = best(&[
                (add(previous[j][2], s.gap_extend), State::Delete),
                (add(previous[j][0], s.gap_open), State::Pair),
                (add(previous[j][1], s.gap_open), State::Insert),
            ]);
            current[j] = [pair.0, insert.0, delete.0];
            set(&mut trace, index, State::Pair, pair.1);
            set(&mut trace, index, State::Insert, insert.1);
            set(&mut trace, index, State::Delete, delete.1);
            if options.mode == AlignmentMode::Local {
                for (state, score) in [
                    (State::Pair, pair.0),
                    (State::Insert, insert.0),
                    (State::Delete, delete.0),
                ] {
                    if score > endpoint.0 {
                        endpoint = (score, i, j, state);
                    }
                }
            }
        }
        if options.mode == AlignmentMode::QueryInReference {
            let value = best(&[
                (current[m][0], State::Pair),
                (current[m][1], State::Insert),
                (current[m][2], State::Delete),
            ]);
            if value.0 > endpoint.0 {
                endpoint = (value.0, i, m, value.1);
            }
        }
        std::mem::swap(&mut previous, &mut current);
    }
    if options.mode == AlignmentMode::Global {
        let value = best(&[
            (previous[m][0], State::Pair),
            (previous[m][1], State::Insert),
            (previous[m][2], State::Delete),
        ]);
        endpoint = (value.0, n, m, value.1);
    }
    let (score, re, qe, mut state) = endpoint;
    let (mut i, mut j) = (re, qe);
    let mut ops = Vec::new();
    while state != State::Stop {
        if i == 0 && j == 0 {
            break;
        }
        if options.mode == AlignmentMode::QueryInReference && j == 0 {
            break;
        }
        let before = predecessor(&trace, i * width + j, state);
        // A local zero cell is a restart boundary, not an emitted paired column.
        if options.mode == AlignmentMode::Local && state == State::Pair && before == State::Stop {
            break;
        }
        if ops.len() >= MAX_ALIGNMENT_COLUMNS {
            return Err(unsupported(
                "pairwise traceback exceeds the one million column output budget",
            ));
        }
        match state {
            State::Pair if i > 0 && j > 0 => {
                i -= 1;
                j -= 1;
            }
            State::Insert if j > 0 => j -= 1,
            State::Delete if i > 0 => i -= 1,
            _ => {
                return Err(invalid(
                    "pairwise traceback reached an inconsistent boundary",
                ));
            }
        }
        ops.push(state);
        state = before;
    }
    ops.reverse();
    Ok(Path {
        score,
        rs: i,
        re,
        qs: j,
        qe,
        ops,
    })
}
fn materialize(
    reference: &Sequence,
    query: &Sequence,
    options: AlignmentOptions,
    path: Path,
) -> Result<Alignment, Error> {
    let mut ra = String::with_capacity(path.ops.len());
    let mut qa = String::with_capacity(path.ops.len());
    let mut runs: Vec<CigarRun> = vec![];
    let (mut r, mut q) = (path.rs, path.qs);
    let rb = reference.as_str().as_bytes();
    let qb = query.as_str().as_bytes();
    for &op in &path.ops {
        let k = match op {
            State::Pair => {
                let k = kind(rb[r], qb[q]);
                ra.push(rb[r] as char);
                qa.push(qb[q] as char);
                r += 1;
                q += 1;
                k
            }
            State::Insert => {
                ra.push('-');
                qa.push(qb[q] as char);
                q += 1;
                CigarKind::Insertion
            }
            State::Delete => {
                ra.push(rb[r] as char);
                qa.push('-');
                r += 1;
                CigarKind::Deletion
            }
            State::Stop => unreachable!(),
        };
        if let Some(last) = runs.last_mut().filter(|last| last.kind == k) {
            last.length += 1;
        } else {
            runs.push(CigarRun { kind: k, length: 1 });
        }
    }
    if r != path.re || q != path.qe {
        return Err(invalid("pairwise materialization disagrees with traceback"));
    }
    let mut stats = AlignmentStats {
        columns: path.ops.len(),
        paired_bases: 0,
        certain_matches: 0,
        ambiguous_compatible: 0,
        mismatches: 0,
        inserted_query_bases: 0,
        deleted_reference_bases: 0,
        reference_span_bases: path.re - path.rs,
        query_span_bases: path.qe - path.qs,
        reference_coverage: 0.,
        query_coverage: 0.,
        certain_identity: None,
        compatible_identity: None,
    };
    let (mut r, mut q) = (path.rs, path.qs);
    let mut differences = vec![];
    let mut cigar = String::new();
    for run in &runs {
        use std::fmt::Write;
        write!(cigar, "{}{}", run.length, run.kind.symbol()).unwrap();
        let (rs, qs) = (r, q);
        let difference = match run.kind {
            CigarKind::Match => {
                stats.certain_matches += run.length;
                r += run.length;
                q += run.length;
                None
            }
            CigarKind::Ambiguous => {
                stats.ambiguous_compatible += run.length;
                r += run.length;
                q += run.length;
                Some(DifferenceKind::AmbiguousCompatible)
            }
            CigarKind::Mismatch => {
                stats.mismatches += run.length;
                r += run.length;
                q += run.length;
                Some(DifferenceKind::Substitution)
            }
            CigarKind::Insertion => {
                stats.inserted_query_bases += run.length;
                q += run.length;
                Some(DifferenceKind::Insertion)
            }
            CigarKind::Deletion => {
                stats.deleted_reference_bases += run.length;
                r += run.length;
                Some(DifferenceKind::Deletion)
            }
        };
        if let Some(kind) = difference {
            differences.push(Difference {
                kind,
                reference_start: rs,
                reference_end: r,
                query_start: qs,
                query_end: q,
                reference_bases: reference.as_str()[rs..r].into(),
                query_bases: query.as_str()[qs..q].into(),
            });
        }
    }
    stats.paired_bases = stats.certain_matches + stats.ambiguous_compatible + stats.mismatches;
    stats.reference_coverage = stats.paired_bases as f64 / reference.len() as f64;
    stats.query_coverage = stats.paired_bases as f64 / query.len() as f64;
    if stats.paired_bases > 0 {
        stats.certain_identity = Some(stats.certain_matches as f64 / stats.paired_bases as f64);
        stats.compatible_identity = Some(
            (stats.certain_matches + stats.ambiguous_compatible) as f64 / stats.paired_bases as f64,
        );
    }
    let mut warnings=vec!["One deterministic optimal alignment is returned; tied alternative alignments are not enumerated. No read-quality calibration or experimental validation is inferred.".into()];
    if stats.ambiguous_compatible > 0 {
        warnings.push(
            "Compatible IUPAC columns remain uncertain and are not counted as certain identity."
                .into(),
        );
    }
    if stats.columns > 0 && stats.paired_bases == 0 {
        warnings.push("No reference/query bases are paired; an all-gap path does not establish sequence identity or coverage.".into());
    }
    if stats.columns == 0 {
        warnings.push("No positive-scoring local alignment was found; zero differences does not establish identity or coverage.".into());
    }
    if path.rs > 0 || path.re < reference.len() || path.qs > 0 || path.qe < query.len() {
        warnings.push(
            "Unaligned input flanks are outside this alignment, not verified matching sequence."
                .into(),
        );
    }
    Ok(Alignment {
        algorithm: "collagene-affine-pairwise-1".into(),
        options,
        score: path.score,
        reference_start: path.rs,
        reference_end: path.re,
        query_start: path.qs,
        query_end: path.qe,
        reference_aligned: ra,
        query_aligned: qa,
        cigar,
        operations: runs,
        differences,
        stats,
        warnings,
        reference_hash: blake3::hash(reference.as_str().as_bytes())
            .to_hex()
            .to_string(),
        query_hash: blake3::hash(query.as_str().as_bytes()).to_hex().to_string(),
    })
}

pub fn align(
    reference: &Sequence,
    query: &Sequence,
    options: &AlignmentOptions,
) -> Result<Alignment, Error> {
    validate(reference, query, options)?;
    if options.mode == AlignmentMode::Global
        && options.scoring.ambiguous_score >= 0
        && reference == query
    {
        if reference.len() > MAX_ALIGNMENT_COLUMNS {
            return Err(unsupported(
                "pairwise identity materialization exceeds one million output columns",
            ));
        }
        let score = reference
            .as_str()
            .bytes()
            .zip(query.as_str().bytes())
            .map(|(a, b)| substitution(a, b, options.scoring))
            .sum();
        return materialize(
            reference,
            query,
            *options,
            Path {
                score,
                rs: 0,
                re: reference.len(),
                qs: 0,
                qe: query.len(),
                ops: vec![State::Pair; reference.len()],
            },
        );
    }
    // Global/local scoring is symmetric. Put the shorter sequence in row storage.
    let mut path =
        if query.len() > reference.len() && options.mode != AlignmentMode::QueryInReference {
            let mut path = compute(query, reference, options)?;
            std::mem::swap(&mut path.rs, &mut path.qs);
            std::mem::swap(&mut path.re, &mut path.qe);
            for op in &mut path.ops {
                *op = match *op {
                    State::Insert => State::Delete,
                    State::Delete => State::Insert,
                    other => other,
                };
            }
            path
        } else {
            compute(reference, query, options)?
        };
    // Empty local alignments have canonical zero coordinates on both inputs.
    if path.ops.is_empty() {
        path.rs = 0;
        path.re = 0;
        path.qs = 0;
        path.qe = 0;
    }
    materialize(reference, query, *options, path)
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub struct ComparisonOptions {
    pub alignment: AlignmentOptions,
    pub both_query_strands: bool,
    pub reference_rotations: bool,
    /// All origins are tried only for references this small; otherwise seed candidates are bounded.
    pub max_origins: usize,
    pub max_total_cells: usize,
}
impl Default for ComparisonOptions {
    fn default() -> Self {
        Self {
            alignment: AlignmentOptions::default(),
            both_query_strands: true,
            reference_rotations: true,
            max_origins: 16,
            max_total_cells: HARD_MAX_CELLS,
        }
    }
}
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct MoleculeComparison {
    pub alignment: Alignment,
    pub reference_revision: String,
    pub query_revision: String,
    /// Alignment coordinates address the reference rotated left by this many bases.
    /// Source reference boundary x is (reference_origin+x) mod reference length.
    pub reference_origin: usize,
    /// Query coordinates address its reverse complement when true; source interval is [len-end,len-start).
    pub query_reversed: bool,
    pub reference_topology: Topology,
    pub query_topology: Topology,
    pub evaluated_candidates: usize,
    pub exhaustive_reference_origins: bool,
    pub exhaustive_query_orientations: bool,
    /// Best score over requested orientations/origins is proven, not uniqueness of an alignment.
    pub optimality_proven: bool,
    pub warnings: Vec<String>,
}
fn rotate(sequence: &Sequence, origin: usize) -> Result<Sequence, Error> {
    if origin == 0 {
        return Ok(sequence.clone());
    }
    Sequence::new(
        sequence.alphabet(),
        &format!(
            "{}{}",
            &sequence.as_str()[origin..],
            &sequence.as_str()[..origin]
        ),
    )
}
fn candidate_origins(
    reference: &Molecule,
    query: &Sequence,
    maximum: usize,
) -> Result<(Vec<usize>, bool), Error> {
    let n = reference.sequence.len();
    if n <= maximum {
        return Ok(((0..n).collect(), true));
    }
    let k = 16.min(query.len());
    let mut votes: BTreeMap<usize, usize> = BTreeMap::new();
    // Spaced exact seeds are just origin candidates. No completeness is attributed to them.
    let mut offsets = vec![
        0,
        query.len() / 4,
        query.len() / 2,
        3 * query.len() / 4,
        query.len() - k,
    ];
    offsets.sort_unstable();
    offsets.dedup();
    for offset in offsets {
        let offset = offset.min(query.len() - k);
        let seed = &query.as_str()[offset..offset + k];
        if seed
            .bytes()
            .any(|b| analysis::base_mask(b).count_ones() != 1)
        {
            continue;
        }
        if k > n {
            continue;
        }
        let rb = reference.sequence.as_str().as_bytes();
        let sb = seed.as_bytes();
        let mut hits = 0usize;
        for start in 0..n {
            if (0..k).all(|j| rb[(start + j) % n] == sb[j]) {
                let origin = (start + n - offset % n) % n;
                *votes.entry(origin).or_default() += 1;
                hits += 1;
                if hits >= 4096 {
                    break;
                }
            }
        }
    }
    let mut ranked: Vec<_> = votes.into_iter().collect();
    ranked.sort_by_key(|&(origin, votes)| (std::cmp::Reverse(votes), origin));
    let mut origins = vec![0];
    for (origin, _) in ranked {
        if origin != 0 && origins.len() < maximum {
            origins.push(origin);
        }
    }
    Ok((origins, false))
}

/// Whole-sequence comparison. Circular origins are exact for identical rotated records,
/// exhaustive only for bounded small references, and otherwise explicitly seed-limited.
/// Global/local/semiglobal single-origin alignment is also available directly through align().
pub fn compare_molecules(
    reference: &Molecule,
    query: &Molecule,
    options: &ComparisonOptions,
) -> Result<MoleculeComparison, Error> {
    reference.validate()?;
    query.validate()?;
    validate(&reference.sequence, &query.sequence, &options.alignment)?;
    if options.alignment.mode != AlignmentMode::Global {
        return Err(unsupported(
            "whole-molecule orientation/origin comparison currently requires global mode; use align for local or query-in-reference alignment",
        ));
    }
    if options.max_origins == 0
        || options.max_origins > 64
        || options.max_total_cells == 0
        || options.max_total_cells > 2 * HARD_MAX_CELLS
    {
        return Err(invalid(
            "comparison requires 1–64 origins and a total budget of at most 128 million cells",
        ));
    }
    let circular = reference.topology == Topology::Circular && options.reference_rotations;
    let orientations = if options.both_query_strands { 2 } else { 1 };
    let mut warnings=vec!["Comparison verifies sequence-model differences, not sample identity or read coverage/quality. Alignment indexes refer to the explicitly rotated reference and oriented query.".into()];
    if reference.topology != query.topology {
        warnings.push("Input topologies differ; the sequence comparison does not establish equivalent physical topology.".into());
    }
    // Identical symbolic rotations attain an upper score bound when ambiguous columns are
    // nonnegative. This linear shortcut proves the best score but does not enumerate ties.
    if options.alignment.scoring.ambiguous_score >= 0
        && reference.sequence.len() == query.sequence.len()
    {
        for reversed in 0..orientations {
            let oriented = if reversed == 1 {
                query.sequence.reverse_complement()?
            } else {
                query.sequence.clone()
            };
            let origin = if circular {
                analysis::rotation_offset(reference.sequence.as_str(), oriented.as_str())
            } else {
                (reference.sequence == oriented).then_some(0)
            };
            if let Some(origin) = origin {
                if oriented.len() > MAX_ALIGNMENT_COLUMNS {
                    return Err(unsupported(
                        "pairwise identity materialization exceeds one million output columns",
                    ));
                }
                let rotated = rotate(&reference.sequence, origin)?;
                let score = rotated
                    .as_str()
                    .bytes()
                    .zip(oriented.as_str().bytes())
                    .map(|(a, b)| substitution(a, b, options.alignment.scoring))
                    .sum();
                let alignment = materialize(
                    &rotated,
                    &oriented,
                    options.alignment,
                    Path {
                        score,
                        rs: 0,
                        re: rotated.len(),
                        qs: 0,
                        qe: oriented.len(),
                        ops: vec![State::Pair; oriented.len()],
                    },
                )?;
                warnings.push("An identical symbolic orientation/rotation attains the scoring upper bound; other optimal placements were not enumerated.".into());
                return Ok(MoleculeComparison {
                    alignment,
                    reference_revision: reference.revision(),
                    query_revision: query.revision(),
                    reference_origin: origin,
                    query_reversed: reversed == 1,
                    reference_topology: reference.topology,
                    query_topology: query.topology,
                    evaluated_candidates: 1,
                    exhaustive_reference_origins: !circular || reference.sequence.len() == 1,
                    exhaustive_query_orientations: orientations == 1,
                    optimality_proven: true,
                    warnings,
                });
            }
        }
    }
    let one = cells(reference.sequence.len(), query.sequence.len())?;
    if one > options.alignment.max_cells || one > options.max_total_cells {
        return Err(unsupported(
            "comparison cannot fit even one affine alignment within the requested cell budgets",
        ));
    }
    let mut spent = 0usize;
    let mut evaluated = 0usize;
    let mut best_alignment: Option<(Alignment, usize, bool)> = None;
    let mut all_origins = true;
    let mut completed_orientations = 0usize;
    let mut budget_complete = true;
    'orientations: for reversed in 0..orientations {
        let oriented = if reversed == 1 {
            query.sequence.reverse_complement()?
        } else {
            query.sequence.clone()
        };
        let (origins, exhaustive) = if circular {
            candidate_origins(reference, &oriented, options.max_origins)?
        } else {
            (vec![0], true)
        };
        all_origins &= exhaustive;
        for origin in origins {
            if spent.saturating_add(one) > options.max_total_cells {
                budget_complete = false;
                break 'orientations;
            }
            let rotated = rotate(&reference.sequence, origin)?;
            let alignment = align(&rotated, &oriented, &options.alignment)?;
            spent += one;
            evaluated += 1;
            if best_alignment
                .as_ref()
                .is_none_or(|(best, _, _)| alignment.score > best.score)
            {
                best_alignment = Some((alignment, origin, reversed == 1));
            }
        }
        completed_orientations += 1;
    }
    let (alignment, origin, reversed) =
        best_alignment.ok_or_else(|| unsupported("no comparison candidate was evaluated"))?;
    let complete_orientation = completed_orientations == orientations && budget_complete;
    if circular && (!all_origins || !budget_complete) {
        warnings.push("Circular origin search is incomplete: exact seed candidates and the total compute budget bound the evaluated rotations. A better untested origin may exist.".into());
    }
    if !complete_orientation {
        warnings.push(
            "Not all requested query orientations/origin-candidate lists were fully evaluated within the total compute budget."
                .into(),
        );
    }
    Ok(MoleculeComparison {
        alignment,
        reference_revision: reference.revision(),
        query_revision: query.revision(),
        reference_origin: origin,
        query_reversed: reversed,
        reference_topology: reference.topology,
        query_topology: query.topology,
        evaluated_candidates: evaluated,
        exhaustive_reference_origins: (!circular || all_origins) && budget_complete,
        exhaustive_query_orientations: complete_orientation,
        optimality_proven: all_origins && complete_orientation && budget_complete,
        warnings,
    })
}
