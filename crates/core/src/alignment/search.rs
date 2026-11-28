//! Bounded literal/IUPAC searches of retained aligned rows, skipping gap columns.
//! This neither realigns inputs nor searches their unreported flanks.
use crate::{alignment::review::Report, sequence::analysis::base_mask};
use collagene_model::Error;
use serde::{Deserialize, Serialize};
use std::collections::VecDeque;

pub const MAX_PATTERN: usize = 4096;
pub const MAX_WORK: usize = 160_000_000;
pub const MAX_HITS: usize = 10_000;
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum Row {
    Reference,
    Query,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum Rows {
    Reference,
    Query,
    Both,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum Mode {
    Literal,
    Compatible,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Options {
    pub rows: Rows,
    pub mode: Mode,
    pub max_work: usize,
    pub max_hits: usize,
}
impl Default for Options {
    fn default() -> Self {
        Self {
            rows: Rows::Both,
            mode: Mode::Literal,
            max_work: 32_000_000,
            max_hits: MAX_HITS,
        }
    }
}
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Hit {
    pub row: Row,
    /// Alignment columns include intervening gaps; zero-based half-open.
    pub columns: (usize, usize),
    /// Ungapped coordinates in the displayed rotated/oriented row.
    pub oriented: (usize, usize),
    /// Original native reference/query coordinates, including an origin crossing.
    pub source_spans: Vec<(usize, usize)>,
}
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Matches {
    pub pattern: String,
    pub options: Options,
    pub hits: Vec<Hit>,
    pub complete: bool,
    pub work: usize,
}
impl Matches {
    /// Standalone wire guard only; this does not repeat the search or establish
    /// source ownership. Hosts retain the exact report/request context.
    pub fn validate(&self) -> Result<(), Error> {
        validate_options(&self.pattern, &self.options)?;
        if self.pattern != self.pattern.to_ascii_uppercase()
            || self.work > self.options.max_work
            || self.hits.len() > self.options.max_hits
            || (!self.complete
                && self.work != self.options.max_work
                && self.hits.len() != self.options.max_hits)
        {
            return Err(Error::Format(
                "Invalid pairwise search result bounds.".into(),
            ));
        }
        for hit in &self.hits {
            if hit.columns.0 >= hit.columns.1
                || hit.columns.1 > crate::alignment::MAX_ALIGNMENT_COLUMNS
                || hit.columns.1 - hit.columns.0 < self.pattern.len()
                || hit.oriented.0 >= hit.oriented.1
                || hit.oriented.1 > collagene_model::MAX_SEQUENCE_LEN
                || hit.oriented.1 - hit.oriented.0 != self.pattern.len()
                || hit.source_spans.is_empty()
                || hit.source_spans.len() > if hit.row == Row::Reference { 2 } else { 1 }
                || hit
                    .source_spans
                    .iter()
                    .any(|(a, b)| a >= b || *b > collagene_model::MAX_SEQUENCE_LEN)
                || hit.source_spans.iter().map(|(a, b)| b - a).sum::<usize>() != self.pattern.len()
                || matches!(
                    (self.options.rows, hit.row),
                    (Rows::Reference, Row::Query) | (Rows::Query, Row::Reference)
                )
            {
                return Err(Error::Format(
                    "Invalid pairwise search hit geometry.".into(),
                ));
            }
        }
        Ok(())
    }
}

enum Matcher {
    Literal {
        prefix: Vec<usize>,
        matched: usize,
    },
    Compatible {
        masks: Vec<Vec<u64>>,
        active: Vec<u64>,
    },
}
impl Matcher {
    fn new(pattern: &[u8], mode: Mode) -> Self {
        match mode {
            Mode::Literal => {
                let mut prefix = vec![0; pattern.len()];
                for i in 1..pattern.len() {
                    let mut j = prefix[i - 1];
                    while j > 0 && pattern[i] != pattern[j] {
                        j = prefix[j - 1];
                    }
                    if pattern[i] == pattern[j] {
                        j += 1;
                    }
                    prefix[i] = j;
                }
                Self::Literal { prefix, matched: 0 }
            }
            Mode::Compatible => {
                let words = pattern.len().div_ceil(64);
                let mut masks = vec![vec![0u64; words]; 16];
                for (i, call) in pattern.iter().enumerate() {
                    for (mask, positions) in masks.iter_mut().enumerate().skip(1) {
                        if mask as u8 & base_mask(*call) != 0 {
                            positions[i / 64] |= 1u64 << (i % 64);
                        }
                    }
                }
                Self::Compatible {
                    masks,
                    active: vec![0; words],
                }
            }
        }
    }
    /// None means the explicit work bound interrupted this row.
    fn advance(
        &mut self,
        call: u8,
        pattern: &[u8],
        work: &mut usize,
        limit: usize,
    ) -> Option<bool> {
        match self {
            Self::Literal { prefix, matched } => {
                loop {
                    if *work == limit {
                        return None;
                    }
                    *work += 1;
                    if call == pattern[*matched] {
                        *matched += 1;
                        break;
                    }
                    if *matched == 0 {
                        break;
                    }
                    *matched = prefix[*matched - 1];
                }
                let hit = *matched == pattern.len();
                if hit {
                    *matched = prefix[*matched - 1];
                }
                Some(hit)
            }
            Self::Compatible { masks, active } => {
                let mut carry = 1;
                let mask = base_mask(call) as usize;
                for (word, allowed) in active.iter_mut().zip(&masks[mask]) {
                    if *work == limit {
                        return None;
                    }
                    *work += 1;
                    let next = *word >> 63;
                    *word = ((*word << 1) | carry) & allowed;
                    carry = next;
                }
                Some(active[(pattern.len() - 1) / 64] & (1u64 << ((pattern.len() - 1) % 64)) != 0)
            }
        }
    }
}

fn validate_options(pattern: &str, options: &Options) -> Result<(), Error> {
    if pattern.is_empty()
        || pattern.len() > MAX_PATTERN
        || !pattern.is_ascii()
        || pattern
            .bytes()
            .any(|call| base_mask(call.to_ascii_uppercase()) == 0)
        || options.max_work == 0
        || options.max_work > MAX_WORK
        || options.max_hits == 0
        || options.max_hits > MAX_HITS
    {
        return Err(Error::Format("Enter 1–4096 DNA/RNA IUPAC calls without gaps or whitespace, and use the bounded search limits.".into()));
    }
    Ok(())
}
pub fn preflight(report: &Report, pattern: &str, options: &Options) -> Result<(), Error> {
    validate_options(pattern, options)?;
    let mut count = CountBytes(0);
    serde_json::to_writer(&mut count, report).map_err(|_| Error::TooLarge)?;
    report.validate_binding(&report.context_revision)?;
    Ok(())
}
struct CountBytes(usize);
impl std::io::Write for CountBytes {
    fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
        self.0 = self
            .0
            .checked_add(bytes.len())
            .filter(|size| *size <= 32 * 1024 * 1024)
            .ok_or_else(|| std::io::Error::other("pairwise search input byte limit"))?;
        Ok(bytes.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}
pub fn revision(report: &Report, pattern: &str, options: &Options) -> Result<String, Error> {
    preflight(report, pattern, options)?;
    crate::encoding::bounded_json::hash(&("pairwise-row-search-1", report, pattern, options))
}
pub fn search(report: &Report, pattern: &str, options: Options) -> Result<Matches, Error> {
    preflight(report, pattern, &options)?;
    let pattern = pattern.to_ascii_uppercase();
    let mut result = Matches {
        pattern,
        options,
        hits: Vec::new(),
        complete: true,
        work: 0,
    };
    let alignment = &report.comparison.alignment;
    for row in [Row::Reference, Row::Query] {
        if matches!(
            (options.rows, row),
            (Rows::Reference, Row::Query) | (Rows::Query, Row::Reference)
        ) {
            continue;
        }
        let (text, mut coordinate) = match row {
            Row::Reference => (&alignment.reference_aligned, alignment.reference_start),
            Row::Query => (&alignment.query_aligned, alignment.query_start),
        };
        let mut matcher = Matcher::new(result.pattern.as_bytes(), options.mode);
        let mut recent = VecDeque::with_capacity(result.pattern.len());
        for (column, call) in text.bytes().enumerate() {
            if result.work == options.max_work {
                result.complete = false;
                return Ok(result);
            }
            result.work += 1;
            if call == b'-' {
                continue;
            }
            coordinate += 1;
            if recent.len() == result.pattern.len() {
                recent.pop_front();
            }
            recent.push_back(column);
            let Some(found) = matcher.advance(
                call,
                result.pattern.as_bytes(),
                &mut result.work,
                options.max_work,
            ) else {
                result.complete = false;
                return Ok(result);
            };
            if !found {
                continue;
            }
            if result.hits.len() == options.max_hits {
                result.complete = false;
                return Ok(result);
            }
            let start = coordinate - result.pattern.len();
            let source_spans = match row {
                Row::Reference => report.map_reference(start, coordinate)?,
                Row::Query => vec![report.map_query(start, coordinate)?],
            };
            result.hits.push(Hit {
                row,
                columns: (
                    *recent.front().ok_or_else(|| {
                        Error::Format("Search coordinate window unavailable".into())
                    })?,
                    column + 1,
                ),
                oriented: (start, coordinate),
                source_spans,
            });
        }
    }
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alignment::{
        AlignmentOptions, ComparisonOptions,
        review::{self as pairwise_review, Options as ReviewOptions},
    };
    use collagene_model::{Molecule, Topology};
    fn make_report(reference: &str, query: &str) -> Report {
        pairwise_review::compare(
            &Molecule::dna("reference", reference, Topology::Linear).unwrap(),
            &Molecule::dna("query", query, Topology::Linear).unwrap(),
            &ReviewOptions::Stored(AlignmentOptions::default()),
        )
        .unwrap()
    }
    #[test]
    fn overlapping_literal_hits_and_gapped_column_windows_keep_both_coordinate_systems() {
        let report = make_report("AAAACCCCGGGG", "AAAATCCCGGGG");
        let result = search(
            &report,
            "aaa",
            Options {
                rows: Rows::Reference,
                ..Options::default()
            },
        )
        .unwrap();
        assert!(result.complete);
        assert_eq!(
            result.hits.iter().map(|h| h.oriented).collect::<Vec<_>>(),
            vec![(0, 3), (1, 4)]
        );
        // Force an independently justified one-base insertion into a long matching region.
        let report = make_report("AACCGGTT", "AACCTGGTT");
        assert_eq!(report.comparison.alignment.reference_aligned, "AACC-GGTT");
        let result = search(
            &report,
            "CCGG",
            Options {
                rows: Rows::Reference,
                ..Options::default()
            },
        )
        .unwrap();
        assert_eq!(result.hits.len(), 1);
        assert_eq!(result.hits[0].columns, (2, 7));
        assert_eq!(result.hits[0].oriented, (2, 6));
        assert_eq!(result.hits[0].source_spans, vec![(2, 6)]);
    }
    #[test]
    fn compatibility_search_is_distinct_from_literal_and_supports_multiple_machine_words() {
        let report = make_report(
            &format!("{}C", "A".repeat(70)),
            &format!("{}C", "A".repeat(70)),
        );
        let pattern = format!("{}Y", "R".repeat(70));
        let literal = search(&report, &pattern, Options::default()).unwrap();
        assert!(literal.complete && literal.hits.is_empty());
        let compatible = search(
            &report,
            &pattern,
            Options {
                mode: Mode::Compatible,
                ..Options::default()
            },
        )
        .unwrap();
        assert!(compatible.complete);
        assert_eq!(compatible.hits.len(), 2);
        assert!(compatible.hits.iter().all(|h| h.columns == (0, 71)));
        let uncertain = make_report("ARNN", "AGCC");
        assert_eq!(
            search(
                &uncertain,
                "AG",
                Options {
                    mode: Mode::Compatible,
                    ..Options::default()
                }
            )
            .unwrap()
            .hits
            .iter()
            .filter(|h| h.row == Row::Reference)
            .count(),
            3
        );
    }
    #[test]
    fn reverse_query_and_rotated_origin_map_to_actual_original_bases() {
        // The non-symmetric circle has no exact forward rotation equal to q.
        // q is independently the reverse complement of GCCGTAAT (origin2).
        let r = Molecule::dna("circle", "ATGCCGTA", Topology::Circular).unwrap();
        let q = Molecule::dna("reverse", "ATTACGGC", Topology::Linear).unwrap();
        let report = pairwise_review::compare(
            &r,
            &q,
            &ReviewOptions::Whole(ComparisonOptions {
                both_query_strands: true,
                reference_rotations: true,
                ..ComparisonOptions::default()
            }),
        )
        .unwrap();
        assert!(report.comparison.query_reversed);
        let result = search(
            &report,
            "GCCG",
            Options {
                rows: Rows::Query,
                ..Options::default()
            },
        )
        .unwrap();
        assert_eq!(result.hits.len(), 1);
        assert_eq!(result.hits[0].source_spans, vec![(4, 8)]);
        // Directly test the report's declared circular display mapping on a valid rotated comparison.
        let q = Molecule::dna("rotated", "TAATGCCG", Topology::Linear).unwrap();
        let report = pairwise_review::compare(
            &r,
            &q,
            &ReviewOptions::Whole(ComparisonOptions {
                both_query_strands: false,
                reference_rotations: true,
                ..ComparisonOptions::default()
            }),
        )
        .unwrap();
        let hits = search(
            &report,
            "TAAT",
            Options {
                rows: Rows::Reference,
                ..Options::default()
            },
        )
        .unwrap();
        assert_eq!(hits.hits[0].source_spans, vec![(6, 8), (0, 2)]);
    }
    #[test]
    fn interrupted_search_never_claims_absence_and_hit_limit_is_not_silent() {
        let report = make_report("AAAAAA", "AAAAAA");
        let cut = search(
            &report,
            "T",
            Options {
                max_work: 3,
                ..Options::default()
            },
        )
        .unwrap();
        assert!(!cut.complete && cut.hits.is_empty());
        let cut = search(
            &report,
            "AA",
            Options {
                max_hits: 2,
                ..Options::default()
            },
        )
        .unwrap();
        assert!(!cut.complete);
        assert_eq!(cut.hits.len(), 2);
        let exact = search(
            &report,
            "AAAAAA",
            Options {
                rows: Rows::Reference,
                max_hits: 1,
                ..Options::default()
            },
        )
        .unwrap();
        assert!(exact.complete);
        assert_eq!(exact.hits.len(), 1);
        for invalid in ["", "A-C", "A C", "🙂"] {
            assert!(search(&report, invalid, Options::default()).is_err());
        }
        assert!(search(&report, &"A".repeat(MAX_PATTERN + 1), Options::default()).is_err());
    }
}
