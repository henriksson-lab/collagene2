//! Bounded library similarity search using the native affine nucleotide aligner.
//! One best local score per subject, not BLAST statistics or sample verification.
use crate::alignment::{self, AlignmentMode, AlignmentOptions, AlignmentStats, Scoring};
use collagene_model::{Alphabet, Error, Molecule, Sequence, Topology};
use serde::{Deserialize, Serialize};
use std::collections::BTreeSet;
pub const ENGINE_VERSION: &str = "sequence-similarity-2";
pub const MAX_SUBJECTS: usize = 128;
pub const MAX_SOURCE_BASES: usize = 8 * 1024 * 1024;
pub const MAX_WORK: usize = 128_000_000;
#[derive(Debug, Clone, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Source {
    pub key: String,
    pub molecule: Molecule,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Options {
    pub scoring: Scoring,
    pub both_query_strands: bool,
    pub max_origins: usize,
    pub max_cells_per_alignment: usize,
    pub max_work: usize,
    pub max_hits: usize,
    pub minimum_score: i32,
}
impl Default for Options {
    fn default() -> Self {
        Self {
            scoring: Scoring::default(),
            both_query_strands: true,
            max_origins: 8,
            max_cells_per_alignment: 8_000_000,
            max_work: 64_000_000,
            max_hits: 100,
            minimum_score: 1,
        }
    }
}
impl Options {
    pub fn validate(&self) -> Result<(), Error> {
        crate::alignment::review::Options::Stored(AlignmentOptions {
            global_storage: crate::alignment::grid::GlobalStorage::Full,
            mode: AlignmentMode::Local,
            scoring: self.scoring,
            max_cells: self.max_cells_per_alignment,
        })
        .validate()?;
        if self.max_origins == 0
            || self.max_origins > 64
            || self.max_work == 0
            || self.max_work > MAX_WORK
            || self.max_hits == 0
            || self.max_hits > MAX_SUBJECTS
            || self.minimum_score < 1
        {
            return Err(bad("invalid origin, work, hit or score limits"));
        }
        Ok(())
    }
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Span {
    pub start: usize,
    pub end: usize,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum Method {
    ExactCanonicalSubstring,
    AffineLocal,
}
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Hit {
    pub subject_key: String,
    pub subject_revision: String,
    pub score: i32,
    pub method: Method,
    pub reference_origin: usize,
    pub reference_spans: Vec<Span>,
    pub query_span: Span,
    pub query_reversed: bool,
    #[cfg_attr(feature = "json-schema", schemars(with = "serde_json::Value"))]
    pub statistics: AlignmentStats,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum Limit {
    Cancelled,
    Work,
    AlignmentCells,
    CircularOrigins,
    AlignmentUnavailable(String),
    IncompatibleAlphabet,
}
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct SubjectReport {
    pub key: String,
    pub revision: String,
    pub length: usize,
    pub alphabet: Alphabet,
    pub topology: Topology,
    pub best: Option<Hit>,
    pub best_score_proven: bool,
    pub score_upper_bound_reached: bool,
    pub alignments_evaluated: usize,
    pub evaluated_origins: Vec<usize>,
    pub limits: Vec<Limit>,
}
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Report {
    pub engine_version: String,
    pub context_revision: String,
    pub query_key: String,
    pub query_revision: String,
    pub query_length: usize,
    pub query_alphabet: Alphabet,
    pub options: Options,
    pub subjects: Vec<SubjectReport>,
    pub hits: Vec<Hit>,
    pub complete: bool,
    pub hits_truncated: bool,
    pub work_used: usize,
    pub warnings: Vec<String>,
}
fn bad(s: &str) -> Error {
    Error::Format(format!("Sequence similarity: {s}"))
}
fn valid_key(s: &str) -> bool {
    !s.is_empty() && s.len() <= 128 && !s.chars().any(|c| c.is_control() || c.is_whitespace())
}
/// Bounds before source validation, hashing, rotating or allocating matrices.
pub fn preflight(query: &Source, subjects: &[Source], options: &Options) -> Result<(), Error> {
    if subjects.len() > MAX_SUBJECTS {
        return Err(Error::TooLarge);
    }
    preflight_borrowed(
        (&query.key, &query.molecule),
        &subjects
            .iter()
            .map(|s| (s.key.as_str(), &s.molecule))
            .collect::<Vec<_>>(),
        options,
    )
}
/// Native/UI callers check borrowed models before cloning a requested collection.
pub fn preflight_borrowed(
    query: (&str, &Molecule),
    subjects: &[(&str, &Molecule)],
    options: &Options,
) -> Result<(), Error> {
    options.validate()?;
    if subjects.is_empty() || subjects.len() > MAX_SUBJECTS || !valid_key(query.0) {
        return Err(bad("invalid query key or subject count"));
    }
    let mut budget = collagene_model::bounds::ModelBudget::default();
    budget.molecule(query.1)?;
    let mut bases = query.1.sequence.len();
    let mut keys = BTreeSet::new();
    for &(key, molecule) in subjects {
        if !valid_key(key) || !keys.insert(key) {
            return Err(bad("subject keys must be bounded and unique"));
        }
        budget.molecule(molecule)?;
        bases = bases
            .checked_add(molecule.sequence.len())
            .ok_or(Error::TooLarge)?;
    }
    if bases > MAX_SOURCE_BASES {
        return Err(Error::TooLarge);
    }
    Ok(())
}

pub fn context_revision(
    query: &Source,
    subjects: &[Source],
    options: &Options,
) -> Result<String, Error> {
    preflight(query, subjects, options)?;
    query.molecule.validate()?;
    if query.molecule.sequence.is_empty() || query.molecule.sequence.alphabet() == Alphabet::Protein
    {
        return Err(bad("a nonempty DNA or RNA query is required"));
    }
    for s in subjects {
        s.molecule.validate()?;
    }
    let query_revision = query.molecule.revision();
    let revisions = subjects
        .iter()
        .map(|s| s.molecule.revision())
        .collect::<Vec<_>>();
    bind_revisions(
        (&query.key, &query_revision),
        &subjects
            .iter()
            .zip(&revisions)
            .map(|(s, r)| (s.key.as_str(), r.as_str()))
            .collect::<Vec<_>>(),
        options,
    )
}
fn spans(start: usize, length: usize, n: usize) -> Vec<Span> {
    if length == 0 {
        return vec![];
    }
    let end = start + length;
    if end <= n {
        vec![Span { start, end }]
    } else {
        vec![
            Span { start, end: n },
            Span {
                start: 0,
                end: end - n,
            },
        ]
    }
}
fn better(a: &Hit, b: &Hit) -> bool {
    a.score > b.score
        || (a.score == b.score
            && (a.statistics.certain_matches > b.statistics.certain_matches
                || (a.statistics.certain_matches == b.statistics.certain_matches
                    && (a.statistics.paired_bases > b.statistics.paired_bases
                        || (a.statistics.paired_bases == b.statistics.paired_bases
                            && (
                                a.query_reversed,
                                a.reference_origin,
                                a.reference_spans[0].start,
                                a.query_span.start,
                            ) < (
                                b.query_reversed,
                                b.reference_origin,
                                b.reference_spans[0].start,
                                b.query_span.start,
                            ))))))
}
fn retain(best: &mut Option<Hit>, hit: Hit) {
    if best.as_ref().is_none_or(|b| better(&hit, b)) {
        *best = Some(hit);
    }
}
struct Work {
    used: usize,
    max: usize,
}
impl Work {
    fn take(&mut self, n: usize) -> bool {
        if n > self.max - self.used {
            false
        } else {
            self.used += n;
            true
        }
    }
}
fn exact(
    subject: &Source,
    query: &Sequence,
    reversed: bool,
    score: i32,
    work: &mut Work,
) -> Result<Option<Hit>, Limit> {
    let reference = &subject.molecule.sequence;
    let n = reference.len();
    let m = query.len();
    if m > n || !query.as_str().bytes().all(|b| b"ACGTU".contains(&b)) {
        return Ok(None);
    }
    let extra = if subject.molecule.topology == Topology::Circular {
        m - 1
    } else {
        0
    };
    if !work.take(n + extra + m) {
        return Err(Limit::Work);
    }
    let extended;
    let haystack = if extra > 0 {
        extended = format!("{}{}", reference.as_str(), &reference.as_str()[..extra]);
        extended.as_str()
    } else {
        reference.as_str()
    };
    let Some(start) = haystack.find(query.as_str()) else {
        return Ok(None);
    };
    let stats = AlignmentStats {
        columns: m,
        paired_bases: m,
        certain_matches: m,
        ambiguous_compatible: 0,
        mismatches: 0,
        inserted_query_bases: 0,
        deleted_reference_bases: 0,
        reference_span_bases: m,
        query_span_bases: m,
        reference_coverage: m as f64 / n as f64,
        query_coverage: 1.0,
        certain_identity: Some(1.0),
        compatible_identity: Some(1.0),
    };
    Ok(Some(Hit {
        subject_key: subject.key.clone(),
        subject_revision: subject.molecule.revision(),
        score,
        method: Method::ExactCanonicalSubstring,
        reference_origin: if subject.molecule.topology == Topology::Circular {
            start
        } else {
            0
        },
        reference_spans: spans(start, m, n),
        query_span: Span { start: 0, end: m },
        query_reversed: reversed,
        statistics: stats,
    }))
}
/// Cancellation is observed before each bounded substring scan and DP invocation.
/// Results describe best scores under the declared stored-query policy. They do
/// not enumerate every tied/repeated local occurrence or establish DNA identity.
pub fn search(
    query: &Source,
    subjects: &[Source],
    options: &Options,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Report, Error> {
    let context_revision = context_revision(query, subjects, options)?;
    let q = &query.molecule.sequence;
    let m = q.len();
    let upper = (m as u64)
        .checked_mul(options.scoring.match_score as u64)
        .filter(|n| *n <= 500_000_000)
        .ok_or_else(|| bad("query and scoring exceed safe score bounds"))? as i32;
    let reverse = if options.both_query_strands {
        Some(q.reverse_complement()?)
    } else {
        None
    };
    let orientations = std::iter::once((false, q))
        .chain(reverse.as_ref().map(|q| (true, q)))
        .collect::<Vec<_>>();
    let mut work = Work {
        used: 0,
        max: options.max_work,
    };
    let mut rows = Vec::new();
    let mut stopped = false;
    for subject in subjects {
        let n = subject.molecule.sequence.len();
        let mut row = SubjectReport {
            key: subject.key.clone(),
            revision: subject.molecule.revision(),
            length: n,
            alphabet: subject.molecule.sequence.alphabet(),
            topology: subject.molecule.topology,
            best: None,
            best_score_proven: false,
            score_upper_bound_reached: false,
            alignments_evaluated: 0,
            evaluated_origins: vec![],
            limits: vec![],
        };
        if stopped || cancelled() {
            row.limits.push(Limit::Cancelled);
            stopped = true;
            rows.push(row);
            continue;
        }
        if subject.molecule.sequence.alphabet() != q.alphabet() {
            row.limits.push(Limit::IncompatibleAlphabet);
            rows.push(row);
            continue;
        }
        if n == 0 {
            row.best_score_proven = true;
            rows.push(row);
            continue;
        }
        for &(reversed, oriented) in &orientations {
            if cancelled() {
                row.limits.push(Limit::Cancelled);
                stopped = true;
                break;
            }
            match exact(subject, oriented, reversed, upper, &mut work) {
                Ok(Some(hit)) => {
                    row.best = Some(hit);
                    row.score_upper_bound_reached = true;
                    row.best_score_proven = true;
                    break;
                }
                Ok(None) => (),
                Err(limit) => {
                    row.limits.push(limit);
                    break;
                }
            }
        }
        if row.best_score_proven || !row.limits.is_empty() {
            rows.push(row);
            continue;
        }
        let circular = subject.molecule.topology == Topology::Circular;
        let count = if circular {
            n.min(options.max_origins)
        } else {
            1
        };
        let origins = (0..count).map(|i| i * n / count).collect::<Vec<_>>();
        let cells = (n + 1).checked_mul(m + 1);
        if cells.is_none_or(|cells| cells > options.max_cells_per_alignment) {
            row.limits.push(Limit::AlignmentCells);
            rows.push(row);
            continue;
        }
        let cells = cells.expect("admitted bounded matrix cell count");
        'origins: for origin in origins {
            let rotated = if origin == 0 {
                subject.molecule.sequence.clone()
            } else {
                Sequence::new(
                    q.alphabet(),
                    &format!(
                        "{}{}",
                        &subject.molecule.sequence.as_str()[origin..],
                        &subject.molecule.sequence.as_str()[..origin]
                    ),
                )?
            };
            for &(reversed, oriented) in &orientations {
                if cancelled() {
                    row.limits.push(Limit::Cancelled);
                    stopped = true;
                    break 'origins;
                }
                if !work.take(cells) {
                    row.limits.push(Limit::Work);
                    break 'origins;
                }
                let alignment = match alignment::align(
                    &rotated,
                    oriented,
                    &AlignmentOptions {
                        global_storage: crate::alignment::grid::GlobalStorage::Full,
                        mode: AlignmentMode::Local,
                        scoring: options.scoring,
                        max_cells: options.max_cells_per_alignment,
                    },
                ) {
                    Ok(a) => a,
                    Err(e) => {
                        row.limits.push(Limit::AlignmentUnavailable(e.to_string()));
                        break 'origins;
                    }
                };
                row.alignments_evaluated += 1;
                if !row.evaluated_origins.contains(&origin) {
                    row.evaluated_origins.push(origin);
                }
                if alignment.score > 0 && alignment.stats.paired_bases > 0 {
                    let start = (origin + alignment.reference_start) % n;
                    let query_span = if reversed {
                        Span {
                            start: m - alignment.query_end,
                            end: m - alignment.query_start,
                        }
                    } else {
                        Span {
                            start: alignment.query_start,
                            end: alignment.query_end,
                        }
                    };
                    let hit = Hit {
                        subject_key: subject.key.clone(),
                        subject_revision: row.revision.clone(),
                        score: alignment.score,
                        method: Method::AffineLocal,
                        reference_origin: origin,
                        reference_spans: spans(
                            start,
                            alignment.reference_end - alignment.reference_start,
                            n,
                        ),
                        query_span,
                        query_reversed: reversed,
                        statistics: alignment.stats,
                    };
                    retain(&mut row.best, hit);
                }
                if row.best.as_ref().is_some_and(|h| h.score == upper) {
                    row.score_upper_bound_reached = true;
                    break 'origins;
                }
            }
        }
        if circular && count < n && !row.score_upper_bound_reached {
            row.limits.push(Limit::CircularOrigins);
        }
        row.best_score_proven = row.score_upper_bound_reached
            || (row.limits.is_empty() && row.alignments_evaluated == count * orientations.len());
        rows.push(row);
    }
    let complete = rows.iter().all(|r| r.best_score_proven);
    let (hits, hits_truncated) = ranked_hits(&rows, options);
    let warnings=vec!["One best-scoring local match per subject. Repeated/tied occurrences are not enumerated; scores are declared-model affine scores, not BLAST E-values or experimental identity.".into(),"The query is its stored linear sequence (and optional reverse complement). Subject circular origins are bounded; untested origins or resource/cancellation limits keep best-score completeness false unless the mathematical score upper bound was reached.".into(),"Canonical full-query exact matches attain the maximum possible score and use bounded substring scans instead of a DP matrix. Ambiguous bases never receive invented certain matches.".into()];
    let report = Report {
        engine_version: ENGINE_VERSION.into(),
        context_revision,
        query_key: query.key.clone(),
        query_revision: query.molecule.revision(),
        query_length: m,
        query_alphabet: q.alphabet(),
        options: options.clone(),
        subjects: rows,
        hits,
        complete,
        hits_truncated,
        work_used: work.used,
        warnings,
    };
    report.validate_binding(&report.context_revision)?;
    Ok(report)
}

fn digest(s: &str) -> bool {
    s.len() == 64
        && s.bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}
impl Report {
    /// Context/shape validation, not a new scientific verdict. Imported summaries
    /// cannot authorize material, sample identity or a cloning plan.
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        self.options.validate()?;
        if self.engine_version != ENGINE_VERSION
            || !digest(expected)
            || self.context_revision != expected
            || !digest(&self.query_revision)
            || !valid_key(&self.query_key)
            || self.query_length == 0
            || self.query_length > MAX_SOURCE_BASES
            || self.query_alphabet == Alphabet::Protein
            || self.subjects.is_empty()
            || self.subjects.len() > MAX_SUBJECTS
            || self.hits.len() > self.options.max_hits
            || self.work_used > self.options.max_work
            || self.warnings.len() > 16
            || self.warnings.iter().any(|w| w.len() > 4096)
        {
            return Err(bad("invalid bounded report or engine/context"));
        }
        let upper = (self.query_length as u64)
            .checked_mul(self.options.scoring.match_score as u64)
            .filter(|n| *n <= 500_000_000)
            .ok_or_else(|| bad("report score bound"))? as i32;
        let mut keys = BTreeSet::new();
        let mut bases = self.query_length;
        for row in &self.subjects {
            bases = bases
                .checked_add(row.length)
                .filter(|n| *n <= MAX_SOURCE_BASES)
                .ok_or(Error::TooLarge)?;
            if !valid_key(&row.key)
                || !keys.insert(row.key.as_str())
                || !digest(&row.revision)
                || row.evaluated_origins.len() > self.options.max_origins
                || row.alignments_evaluated > 2 * self.options.max_origins
                || row.limits.len() > 8
                || row
                    .limits
                    .iter()
                    .any(|l| matches!(l,Limit::AlignmentUnavailable(s) if s.len()>4096))
                || row.evaluated_origins.iter().any(|o| *o >= row.length)
                || row.evaluated_origins.windows(2).any(|w| w[0] >= w[1])
                || (row.score_upper_bound_reached
                    && (!row.best_score_proven
                        || row.best.as_ref().is_none_or(|h| h.score != upper)))
                || (row.best_score_proven
                    && !row.score_upper_bound_reached
                    && !row.limits.is_empty())
            {
                return Err(bad("invalid subject report"));
            }
            if row.alphabet != self.query_alphabet
                && (row.best.is_some() || row.best_score_proven || row.alignments_evaluated > 0)
            {
                return Err(bad("incompatible subject has search evidence"));
            }
            if row.best_score_proven && !row.score_upper_bound_reached {
                let required = if row.length == 0 {
                    0
                } else if row.topology == Topology::Circular {
                    row.length
                } else {
                    1
                };
                let orientations = if self.options.both_query_strands {
                    2
                } else {
                    1
                };
                if row.evaluated_origins.len() != required
                    || row.alignments_evaluated != required * orientations
                {
                    return Err(bad("unsearched origins/strands claimed complete"));
                }
            }
            if let Some(hit) = &row.best {
                self.validate_hit(row, hit, upper)?;
            }
        }
        if bind_revisions(
            (&self.query_key, &self.query_revision),
            &self
                .subjects
                .iter()
                .map(|s| (s.key.as_str(), s.revision.as_str()))
                .collect::<Vec<_>>(),
            &self.options,
        )? != self.context_revision
        {
            return Err(bad("report source/options binding changed"));
        }
        let (hits, truncated) = ranked_hits(&self.subjects, &self.options);
        if self.hits != hits
            || self.hits_truncated != truncated
            || self.complete != self.subjects.iter().all(|r| r.best_score_proven)
        {
            return Err(bad("report ranking/completeness changed"));
        }
        Ok(())
    }
    fn validate_hit(&self, row: &SubjectReport, h: &Hit, upper: i32) -> Result<(), Error> {
        if h.subject_key != row.key
            || h.subject_revision != row.revision
            || h.score < 1
            || h.score > upper
            || row.length == 0
            || h.reference_origin >= row.length
            || h.query_span.start >= h.query_span.end
            || h.query_span.end > self.query_length
            || h.reference_spans.is_empty()
            || h.reference_spans.len() > 2
            || h.reference_spans
                .iter()
                .any(|s| s.start >= s.end || s.end > row.length)
            || (h.reference_spans.len() == 2
                && (h.reference_spans[0].end != row.length
                    || h.reference_spans[1].start != 0
                    || h.reference_spans[0].start < h.reference_spans[1].end))
            || (!self.options.both_query_strands && h.query_reversed)
        {
            return Err(bad("invalid original-source hit coordinates"));
        }
        let s = &h.statistics;
        let counts = [
            s.columns,
            s.paired_bases,
            s.certain_matches,
            s.ambiguous_compatible,
            s.mismatches,
            s.inserted_query_bases,
            s.deleted_reference_bases,
            s.reference_span_bases,
            s.query_span_bases,
        ];
        if counts.iter().any(|n| *n > 2 * MAX_SOURCE_BASES) {
            return Err(Error::TooLarge);
        }
        let reference_span = h
            .reference_spans
            .iter()
            .map(|s| s.end - s.start)
            .sum::<usize>();
        if s.paired_bases == 0
            || s.paired_bases != s.certain_matches + s.ambiguous_compatible + s.mismatches
            || s.columns != s.paired_bases + s.inserted_query_bases + s.deleted_reference_bases
            || s.reference_span_bases != reference_span
            || reference_span != s.paired_bases + s.deleted_reference_bases
            || s.query_span_bases != h.query_span.end - h.query_span.start
            || s.query_span_bases != s.paired_bases + s.inserted_query_bases
            || s.reference_coverage != s.paired_bases as f64 / row.length as f64
            || s.query_coverage != s.paired_bases as f64 / self.query_length as f64
            || s.certain_identity != Some(s.certain_matches as f64 / s.paired_bases as f64)
            || s.compatible_identity
                != Some((s.certain_matches + s.ambiguous_compatible) as f64 / s.paired_bases as f64)
            || (h.method == Method::ExactCanonicalSubstring
                && (h.score != upper
                    || s.certain_matches != self.query_length
                    || s.columns != self.query_length))
        {
            return Err(bad("hit statistics changed"));
        }
        Ok(())
    }
}

fn ranked_hits(rows: &[SubjectReport], options: &Options) -> (Vec<Hit>, bool) {
    let mut hits = rows
        .iter()
        .filter_map(|r| r.best.clone())
        .filter(|h| h.score >= options.minimum_score)
        .collect::<Vec<_>>();
    hits.sort_by(|a, b| {
        b.score
            .cmp(&a.score)
            .then_with(|| {
                b.statistics
                    .certain_matches
                    .cmp(&a.statistics.certain_matches)
            })
            .then_with(|| b.statistics.paired_bases.cmp(&a.statistics.paired_bases))
            .then_with(|| a.subject_key.cmp(&b.subject_key))
    });
    let truncated = hits.len() > options.max_hits;
    hits.truncate(options.max_hits);
    (hits, truncated)
}

fn bind_revisions(
    query: (&str, &str),
    subjects: &[(&str, &str)],
    options: &Options,
) -> Result<String, Error> {
    crate::encoding::bounded_json::hash(&(ENGINE_VERSION, query.0, query.1, subjects, options))
}
