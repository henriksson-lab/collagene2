//! Independently implemented affine-gap nucleotide pairwise alignment.
//! All coordinates are zero-based and half-open on the input reference/query.
//! This module invokes no external alignment engine and does not infer read quality.
use crate::sequence::analysis;
use collagene_model::{Alphabet, Error, Molecule, Sequence, Topology};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

pub const DEFAULT_MAX_CELLS: usize = 16_000_000;
pub const HARD_MAX_CELLS: usize = 64_000_000;
pub const MAX_ALIGNMENT_COLUMNS: usize = 1_000_000;
pub mod grid;
mod kernel;
use kernel::{Path, State};

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
    pub global_storage: grid::GlobalStorage,
}
impl Default for AlignmentOptions {
    fn default() -> Self {
        Self {
            mode: AlignmentMode::Global,
            scoring: Scoring::default(),
            max_cells: DEFAULT_MAX_CELLS,
            global_storage: grid::GlobalStorage::default(),
        }
    }
}
impl AlignmentOptions {
    /// Canonical score/index/storage admission before recurrence or boundary
    /// initialization. Also reusable by bounded documentary validators.
    pub fn validate_for_lengths(
        &self,
        reference_length: usize,
        query_length: usize,
    ) -> Result<(), Error> {
        self.global_storage.validate(self.mode)?;
        let s = self.scoring;
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
        if self.max_cells == 0 || self.max_cells > HARD_MAX_CELLS {
            return Err(invalid(
                "max_cells must be between 1 and the 64 million cell hard limit",
            ));
        }
        let magnitude = [
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
        let bound = u64::try_from(reference_length)
            .ok()
            .and_then(|n| {
                u64::try_from(query_length)
                    .ok()
                    .and_then(|m| n.checked_add(m))
            })
            .and_then(|length| length.checked_mul(magnitude))
            .ok_or(Error::TooLarge)?;
        if bound > 500_000_000 {
            return Err(unsupported(
                "input lengths and scoring exceed the safe pairwise score bound",
            ));
        }
        Ok(())
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
    pub search: grid::SearchEvidence,
    pub reference_length: usize,
    pub query_length: usize,
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
impl Alignment {
    pub fn validate_search(
        &self,
        reference_length: usize,
        query_length: usize,
    ) -> Result<(), Error> {
        self.options
            .validate_for_lengths(reference_length, query_length)?;
        if self.reference_length != reference_length
            || self.query_length != query_length
            || self.algorithm != "collagene-affine-pairwise-2"
        {
            return Err(invalid(
                "alignment source lengths or current algorithm changed",
            ));
        }
        if self.reference_aligned.len() > MAX_ALIGNMENT_COLUMNS
            || self.query_aligned.len() > MAX_ALIGNMENT_COLUMNS
            || self.reference_aligned.len() != self.query_aligned.len()
            || self.operations.len() > MAX_ALIGNMENT_COLUMNS
            || self.differences.len() > MAX_ALIGNMENT_COLUMNS
            || self.reference_start > self.reference_end
            || self.reference_end > reference_length
            || self.query_start > self.query_end
            || self.query_end > query_length
        {
            return Err(Error::TooLarge);
        }
        match self.options.mode {
            AlignmentMode::Global
                if self.reference_start != 0
                    || self.reference_end != reference_length
                    || self.query_start != 0
                    || self.query_end != query_length =>
            {
                return Err(invalid("global alignment must cover both complete inputs"));
            }
            AlignmentMode::QueryInReference
                if self.query_start != 0 || self.query_end != query_length =>
            {
                return Err(invalid(
                    "query-in-reference alignment must cover the complete query",
                ));
            }
            AlignmentMode::Local
                if (self.reference_aligned.is_empty()
                    && (self.reference_start != 0
                        || self.reference_end != 0
                        || self.query_start != 0
                        || self.query_end != 0))
                    || (!self.reference_aligned.is_empty() && self.score <= 0) =>
            {
                return Err(invalid(
                    "local alignment must have positive score or canonical empty coordinates",
                ));
            }
            _ => {}
        }
        review::validate_alignment(self, reference_length, query_length)?;
        self.search.validate(
            reference_length,
            query_length,
            self.options.mode,
            self.options.scoring,
            self.score,
            self.options.max_cells,
        )?;
        if self.search.stop_reason == grid::StopReason::ExactSymbolicIdentity {
            if self.reference_start != 0
                || self.query_start != 0
                || self.reference_end != reference_length
                || self.query_end != query_length
                || self.reference_aligned.len() != reference_length
                || self.reference_aligned != self.query_aligned
                || self.reference_aligned.contains('-')
            {
                return Err(invalid(
                    "exact-symbolic proof differs from selected complete rows",
                ));
            }
            let score: i64 = self
                .reference_aligned
                .bytes()
                .map(|call| i64::from(substitution(call, call, self.options.scoring)))
                .sum();
            if score != i64::from(self.score) {
                return Err(invalid("exact-symbolic score differs from native rows"));
            }
        }
        Ok(())
    }
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
    options.validate_for_lengths(reference.len(), query.len())
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
fn materialize(
    reference: &Sequence,
    query: &Sequence,
    options: AlignmentOptions,
    path: Path,
    search: grid::SearchEvidence,
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
    let mut warnings = vec![if search.unrestricted_optimum_proven {
        "One deterministic scoring-model optimum is returned; tied correspondences are not enumerated. No read-quality calibration or experimental validation is inferred.".into()
    } else {
        "One deterministic best path in the last completed row domain is returned; unrestricted optimality remains unproven. Tied correspondences are not enumerated, and no read-quality calibration or experimental validation is inferred.".into()
    }];
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
    let result = Alignment {
        algorithm: "collagene-affine-pairwise-2".into(),
        reference_length: reference.len(),
        query_length: query.len(),
        search,
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
    };
    result.validate_search(reference.len(), query.len())?;
    Ok(result)
}

pub fn align(
    reference: &Sequence,
    query: &Sequence,
    options: &AlignmentOptions,
) -> Result<Alignment, Error> {
    validate(reference, query, options)?;
    if let Some((path, evidence)) = identity_path(reference, query, options)? {
        return materialize(reference, query, *options, path, evidence);
    }
    let mut work = grid::CellWork::new(options.max_cells)?;
    let mut budget = grid::CellBudget::new(&mut work, options.max_cells)?;
    let (path, evidence) = kernel::align(reference, query, options, &mut budget)?;
    materialize(reference, query, *options, path, evidence)
}
fn identity_path(
    reference: &Sequence,
    query: &Sequence,
    options: &AlignmentOptions,
) -> Result<Option<(Path, grid::SearchEvidence)>, Error> {
    if options.mode != AlignmentMode::Global
        || options.scoring.ambiguous_score < 0
        || reference != query
    {
        return Ok(None);
    }
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
    Ok(Some((
        Path {
            score,
            rs: 0,
            re: reference.len(),
            qs: 0,
            qe: query.len(),
            ops: vec![State::Pair; reference.len()],
        },
        grid::SearchEvidence {
            storage: options.global_storage,
            attempts: vec![],
            evaluated_cells: 0,
            unrestricted_optimum_proven: true,
            stop_reason: grid::StopReason::ExactSymbolicIdentity,
        },
    )))
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
    pub evaluated_cells: usize,
    #[serde(
        deserialize_with = "collagene_model::decode::bounded_vec::<_, grid::CandidateEvidence, 128>"
    )]
    pub candidate_evidence: Vec<grid::CandidateEvidence>,
    pub exhaustive_reference_origins: bool,
    pub exhaustive_query_orientations: bool,
    /// Best score over requested orientations/origins is proven, not uniqueness of an alignment.
    pub optimality_proven: bool,
    pub warnings: Vec<String>,
}
impl MoleculeComparison {
    pub fn stored(
        reference: &Molecule,
        query: &Molecule,
        options: &AlignmentOptions,
    ) -> Result<Self, Error> {
        collagene_model::bounds::model_preflight(&[reference, query])?;
        reference.validate()?;
        query.validate()?;
        let alignment = align(&reference.sequence, &query.sequence, options)?;
        let search = alignment.search.clone();
        let score = alignment.score;
        Ok(Self {
            reference_revision: reference.revision(),
            query_revision: query.revision(),
            reference_origin: 0,
            query_reversed: false,
            reference_topology: reference.topology,
            query_topology: query.topology,
            evaluated_candidates: 1,
            evaluated_cells: search.evaluated_cells,
            candidate_evidence: vec![grid::CandidateEvidence {
                reference_origin: 0,
                query_reversed: false,
                score: Some(score),
                search: Some(search.clone()),
                disposition: if search.unrestricted_optimum_proven {
                    grid::CandidateDisposition::CertifiedOwnOptimum
                } else {
                    grid::CandidateDisposition::Limited
                },
            }],
            exhaustive_reference_origins: true,
            exhaustive_query_orientations: true,
            optimality_proven: search.unrestricted_optimum_proven,
            alignment,
            warnings: vec!["Stored origin/orientation policy: no rotation or reverse-complement search. Local/query-in-reference flanks are not assessed; unproven row bands do not establish unrestricted optimum.".into()],
        })
    }
    /// Bounded documentary arithmetic and correspondence of actual candidate
    /// receipts. This does not repeat DP, infer source ownership or prove truth.
    pub fn validate_search(&self, n: usize, m: usize, max_total_cells: usize) -> Result<(), Error> {
        use grid::CandidateDisposition as D;
        if self.candidate_evidence.is_empty()
            || self.candidate_evidence.len() > 128
            || max_total_cells == 0
            || max_total_cells > 2 * HARD_MAX_CELLS
        {
            return Err(invalid("invalid comparison candidate/work bounds"));
        }
        self.alignment.validate_search(n, m)?;
        let mut cells = 0usize;
        let mut completed = 0usize;
        let mut selected = 0usize;
        let mut keys = std::collections::BTreeSet::new();
        for candidate in &self.candidate_evidence {
            if candidate.reference_origin >= n.max(1)
                || (self.reference_topology != Topology::Circular
                    && candidate.reference_origin != 0)
                || !keys.insert((candidate.reference_origin, candidate.query_reversed))
            {
                return Err(invalid("comparison candidate coordinate/identity changed"));
            }
            match (&candidate.search, candidate.score) {
                (Some(search), Some(score)) => {
                    if search.storage != self.alignment.options.global_storage
                        || score > self.alignment.score
                    {
                        return Err(invalid("comparison candidate policy/score changed"));
                    }
                    search.validate(
                        n,
                        m,
                        self.alignment.options.mode,
                        self.alignment.options.scoring,
                        score,
                        self.alignment.options.max_cells,
                    )?;
                    cells = cells
                        .checked_add(search.evaluated_cells)
                        .filter(|cells| *cells <= max_total_cells)
                        .ok_or(Error::TooLarge)?;
                    completed += 1;
                    match candidate.disposition {
                        D::CertifiedOwnOptimum if search.unrestricted_optimum_proven => {}
                        D::DischargedByIncumbent
                            if !search.unrestricted_optimum_proven
                                && search.stop_reason == grid::StopReason::ComparisonIncumbent
                                && search.attempts.last().is_some_and(|attempt| {
                                    attempt.off_band_bound.certifies(self.alignment.score)
                                }) => {}
                        D::Limited if !search.unrestricted_optimum_proven => {}
                        _ => return Err(invalid("comparison candidate disposition/proof changed")),
                    }
                    if candidate.reference_origin == self.reference_origin
                        && candidate.query_reversed == self.query_reversed
                    {
                        if score != self.alignment.score || search != &self.alignment.search {
                            return Err(invalid("selected comparison attempt changed"));
                        }
                        selected += 1;
                    }
                }
                (None, None) if candidate.disposition == D::Unsearched => {}
                _ => return Err(invalid("comparison candidate completeness changed")),
            }
        }
        if selected != 1 || completed != self.evaluated_candidates || cells != self.evaluated_cells
        {
            return Err(invalid("comparison selected/cumulative work changed"));
        }
        let exact = self.alignment.search.stop_reason == grid::StopReason::ExactSymbolicIdentity;
        if self.optimality_proven
            && !exact
            && (!self.exhaustive_reference_origins
                || !self.exhaustive_query_orientations
                || !self.alignment.search.unrestricted_optimum_proven
                || self.candidate_evidence.iter().any(|candidate| {
                    !matches!(
                        candidate.disposition,
                        D::CertifiedOwnOptimum | D::DischargedByIncumbent
                    )
                }))
        {
            return Err(invalid(
                "unproven candidate/origin scope cannot claim comparison optimum",
            ));
        }
        Ok(())
    }
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
    collagene_model::bounds::model_preflight(&[reference, query])?;
    reference.validate()?;
    query.validate()?;
    validate(&reference.sequence, &query.sequence, &options.alignment)?;
    if options.alignment.mode != AlignmentMode::Global {
        return Err(unsupported(
            "whole-molecule comparison requires global mode",
        ));
    }
    if !(1..=64).contains(&options.max_origins)
        || options.max_total_cells == 0
        || options.max_total_cells > 2 * HARD_MAX_CELLS
    {
        return Err(invalid(
            "comparison requires1–64origins and at most128million totalcells",
        ));
    }
    use grid::{
        CandidateDisposition as Disposition, CandidateEvidence, Domain, GlobalStorage,
        SearchEvidence, StopReason,
    };
    let circular = reference.topology == Topology::Circular && options.reference_rotations;
    let orientations = if options.both_query_strands { 2 } else { 1 };
    let mut warnings=vec!["Comparison describes sequence-model differences, not sample identity or read depth/quality. Coordinates address the declared rotated reference and oriented query; best score does not enumerate tied correspondences.".into()];
    if reference.topology != query.topology {
        warnings.push("Input topology declarations differ; string comparison does not establish equivalent physical topology.".into());
    }
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
                let rotated = rotate(&reference.sequence, origin)?;
                let (path, search) = identity_path(&rotated, &oriented, &options.alignment)?
                    .expect("equal symbolic sequences");
                let score = path.score;
                let alignment =
                    materialize(&rotated, &oriented, options.alignment, path, search.clone())?;
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
                    evaluated_cells: 0,
                    candidate_evidence: vec![CandidateEvidence {
                        reference_origin: origin,
                        query_reversed: reversed == 1,
                        score: Some(score),
                        search: Some(search),
                        disposition: Disposition::CertifiedOwnOptimum,
                    }],
                    exhaustive_reference_origins: !circular || reference.sequence.len() == 1,
                    exhaustive_query_orientations: orientations == 1,
                    optimality_proven: true,
                    warnings,
                });
            }
        }
    }
    let mut candidates = Vec::new();
    let mut all_origins = true;
    for reversed in 0..orientations {
        let oriented = if reversed == 1 {
            query.sequence.reverse_complement()?
        } else {
            query.sequence.clone()
        };
        let (origins, complete) = if circular {
            candidate_origins(reference, &oriented, options.max_origins)?
        } else {
            (vec![0], true)
        };
        all_origins &= complete;
        for origin in origins {
            candidates.push(CandidateEvidence {
                reference_origin: origin,
                query_reversed: reversed == 1,
                score: None,
                search: None,
                disposition: Disposition::Unsearched,
            });
        }
    }
    let (n, m) = (reference.sequence.len(), query.sequence.len());
    let mut work = grid::CellWork::new(options.max_total_cells)?;
    let mut best: Option<(usize, Path)> = None;
    // Breadth-first band rounds expose cheap near-matches before any wrong
    // orientation consumes all work by widening. Only the winning Path is kept.
    let maximum_attempts = match options.alignment.global_storage {
        GlobalStorage::Full => 1,
        GlobalStorage::Adaptive { max_attempts, .. } => max_attempts,
    };
    for round in 0..maximum_attempts {
        let mut progressed = false;
        for (index, candidate) in candidates.iter_mut().enumerate() {
            if matches!(
                candidate.disposition,
                Disposition::CertifiedOwnOptimum | Disposition::DischargedByIncumbent
            ) {
                continue;
            }
            if round > 0 {
                let Some(search) = &mut candidate.search else {
                    continue;
                };
                if search.stop_reason != StopReason::ComparisonIncumbent {
                    continue;
                }
                let incumbent = best
                    .as_ref()
                    .expect("a completed search has an incumbent")
                    .1
                    .score;
                if search
                    .attempts
                    .last()
                    .unwrap()
                    .off_band_bound
                    .certifies(incumbent)
                {
                    candidate.disposition = Disposition::DischargedByIncumbent;
                    continue;
                }
            }
            let (domain, used) = if let Some(search) = &candidate.search {
                let GlobalStorage::Adaptive { maximum_radius, .. } =
                    options.alignment.global_storage
                else {
                    continue;
                };
                let Domain::Corridor { radius } = search.attempts.last().unwrap().domain else {
                    continue;
                };
                let next = kernel::next_radius(radius, maximum_radius, n, m)?;
                if next <= radius {
                    candidate.search.as_mut().unwrap().stop_reason = StopReason::RadiusLimit;
                    continue;
                }
                (Domain::Corridor { radius: next }, search.evaluated_cells)
            } else {
                (
                    match options.alignment.global_storage {
                        GlobalStorage::Full => Domain::Full,
                        GlobalStorage::Adaptive { initial_radius, .. } => Domain::Corridor {
                            radius: initial_radius.min(n.min(m)),
                        },
                    },
                    0,
                )
            };
            let mut budget =
                grid::CellBudget::resume(&mut work, options.alignment.max_cells, used)?;
            if domain.cell_count(n, m)? > budget.remaining() {
                if let Some(search) = &mut candidate.search {
                    search.stop_reason = StopReason::CellLimit;
                    candidate.disposition = Disposition::Limited;
                }
                continue;
            }
            let oriented = if candidate.query_reversed {
                query.sequence.reverse_complement()?
            } else {
                query.sequence.clone()
            };
            let rotated = rotate(&reference.sequence, candidate.reference_origin)?;
            let before = budget.used();
            let path =
                match kernel::attempt(&rotated, &oriented, &options.alignment, domain, &mut budget)
                {
                    Ok(path) => path,
                    Err(Error::TooLarge) => {
                        if let Some(search) = &mut candidate.search {
                            search.stop_reason = StopReason::MemoryLimit;
                            candidate.disposition = Disposition::Limited;
                        }
                        continue;
                    }
                    Err(error) => return Err(error),
                };
            progressed = true;
            let receipt = grid::Attempt::completed(
                domain,
                n,
                m,
                options.alignment.scoring,
                path.score,
                budget.used() - before,
            )?;
            let certified = receipt.unrestricted_optimum_proven;
            let search = candidate.search.get_or_insert_with(|| SearchEvidence {
                storage: options.alignment.global_storage,
                attempts: vec![],
                evaluated_cells: 0,
                unrestricted_optimum_proven: false,
                stop_reason: StopReason::ComparisonIncumbent,
            });
            search.attempts.push(receipt);
            search.evaluated_cells = budget.used();
            search.unrestricted_optimum_proven = certified;
            search.stop_reason = if certified {
                StopReason::Certified
            } else if round + 1 == maximum_attempts {
                StopReason::AttemptLimit
            } else {
                StopReason::ComparisonIncumbent
            };
            candidate.score = Some(path.score);
            candidate.disposition = if certified {
                Disposition::CertifiedOwnOptimum
            } else {
                Disposition::Limited
            };
            if best.as_ref().is_none_or(|(best_index, best)| {
                path.score > best.score || (path.score == best.score && index <= *best_index)
            }) {
                best = Some((index, path));
            }
        }
        if !progressed {
            break;
        }
    }
    let (best_index, path) =
        best.ok_or_else(|| unsupported("no comparison candidate fit its row-grid budgets"))?;
    let winning_score = path.score;
    // A common incumbent may discharge another candidate's omitted paths while
    // that candidate's OWN optimum stays explicitly unproven.
    for candidate in &mut candidates {
        if candidate.disposition != Disposition::CertifiedOwnOptimum {
            if let Some(search) = &mut candidate.search {
                if search
                    .attempts
                    .last()
                    .unwrap()
                    .off_band_bound
                    .certifies(winning_score)
                {
                    candidate.disposition = Disposition::DischargedByIncumbent;
                    search.stop_reason = StopReason::ComparisonIncumbent;
                }
            }
        }
    }
    let selected = &candidates[best_index];
    let origin = selected.reference_origin;
    let reversed = selected.query_reversed;
    let search = selected
        .search
        .clone()
        .expect("selected completed candidate");
    let rotated = rotate(&reference.sequence, origin)?;
    let oriented = if reversed {
        query.sequence.reverse_complement()?
    } else {
        query.sequence.clone()
    };
    let alignment = materialize(&rotated, &oriented, options.alignment, path, search)?;
    let completed = candidates.iter().all(|candidate| {
        matches!(
            candidate.disposition,
            Disposition::CertifiedOwnOptimum | Disposition::DischargedByIncumbent
        )
    });
    if !all_origins {
        warnings.push(
            "Circular origin search is incomplete and seed-limited: a better untested origin may exist.".into(),
        );
    }
    if !completed {
        warnings.push("Requested origin/orientation search is incomplete: some candidates remain unsearched or unproven under actual cumulative band/cell/memory limits; unrestricted optimum is not established.".into());
    }
    Ok(MoleculeComparison {
        alignment,
        reference_revision: reference.revision(),
        query_revision: query.revision(),
        reference_origin: origin,
        query_reversed: reversed,
        reference_topology: reference.topology,
        query_topology: query.topology,
        evaluated_candidates: candidates
            .iter()
            .filter(|candidate| candidate.score.is_some())
            .count(),
        evaluated_cells: work.used(),
        candidate_evidence: candidates,
        exhaustive_reference_origins: (!circular || all_origins) && completed,
        exhaustive_query_orientations: completed,
        optimality_proven: all_origins && completed,
        warnings,
    })
}

pub mod context;
pub mod dotplot;
pub mod export;
pub mod history;
pub mod review;
pub mod search;
