//! Cached alignment-page coordinates. This is presentation of an already validated
//! alignment, not an alignment algorithm or proof of optimality.
use collagene_core::{alignment::Alignment, sequence::analysis};
use std::fmt::Write as _;

pub(crate) struct PageIndex {
    pages: Vec<(usize, usize)>,
    differences: Vec<Option<usize>>,
}
impl PageIndex {
    pub(crate) fn new(a: &Alignment) -> Self {
        let mut pages = Vec::new();
        let mut differences = vec![None; a.differences.len()];
        let (mut r, mut q) = (a.reference_start, a.query_start);
        let mut difference = 0usize;
        for (column, (rb, qb)) in a
            .reference_aligned
            .bytes()
            .zip(a.query_aligned.bytes())
            .enumerate()
        {
            if column % 600 == 0 {
                pages.push((r, q));
            }
            while let Some(d) = a.differences.get(difference) {
                if d.reference_start != r || d.query_start != q {
                    break;
                }
                differences[difference] = Some(column);
                difference += 1;
            }
            r += usize::from(rb != b'-');
            q += usize::from(qb != b'-');
        }
        if pages.is_empty() {
            pages.push((r, q));
        }
        Self { pages, differences }
    }
    pub(crate) fn difference_column(&self, difference: usize) -> Option<usize> {
        self.differences.get(difference).copied().flatten()
    }
    pub(crate) fn page_count(&self) -> usize {
        self.pages.len()
    }
}
pub(crate) fn render_page(a: &Alignment, index: &PageIndex, page: usize) -> String {
    let start = page.saturating_mul(600).min(a.reference_aligned.len());
    let end = (start + 600).min(a.reference_aligned.len());
    let (mut r, mut q) = index
        .pages
        .get(page)
        .copied()
        .unwrap_or((a.reference_end, a.query_end));
    let mut text = String::new();
    for column in (start..end).step_by(60) {
        let stop = (column + 60).min(end);
        let reference = &a.reference_aligned[column..stop];
        let query = &a.query_aligned[column..stop];
        let marks = reference
            .bytes()
            .zip(query.bytes())
            .map(|(rb, qb)| {
                if rb == b'-' || qb == b'-' {
                    ' '
                } else {
                    let r = analysis::base_mask(rb);
                    let q = analysis::base_mask(qb);
                    if r & q == 0 {
                        '×'
                    } else if r.count_ones() == 1 && q.count_ones() == 1 {
                        '|'
                    } else {
                        '?'
                    }
                }
            })
            .collect::<String>();
        let _ = writeln!(
            text,
            "R {r:>7} {reference}\n          {marks}\nQ {q:>7} {query}\n"
        );
        r += reference.bytes().filter(|b| *b != b'-').count();
        q += query.bytes().filter(|b| *b != b'-').count();
    }
    text
}

#[cfg(test)]
mod tests {
    use super::*;
    use collagene_core::alignment::{
        AlignmentMode, AlignmentOptions, AlignmentStats, CigarKind, CigarRun, Difference,
        DifferenceKind, Scoring,
    };

    #[test]
    fn gaps_across_a_page_boundary_keep_native_offsets_and_empty_difference_spans() {
        // Explicit valid alignment path: insert T at column599, match G at600,
        // delete C at601. No DP or helper under test constructs these coordinates.
        let reference_row = format!("{}-GC{}", "A".repeat(599), "T".repeat(64));
        let query_row = format!("{}TG-{}", "A".repeat(599), "T".repeat(64));
        // Source strings are only hash evidence here; the display path is explicit.
        let reference_source = format!(
            "{}{}GC{}{}",
            "N".repeat(10),
            "A".repeat(599),
            "T".repeat(64),
            "N".repeat(5)
        );
        let query_source = format!(
            "{}{}TG{}{}",
            "N".repeat(20),
            "A".repeat(599),
            "T".repeat(64),
            "N".repeat(7)
        );
        let alignment = Alignment {
            algorithm: "independently declared gap-containing presentation fixture".into(),
            options: AlignmentOptions {
                mode: AlignmentMode::Local,
                scoring: Scoring {
                    gap_open: -2,
                    ..Scoring::default()
                },
                ..AlignmentOptions::default()
            },
            score: 1324,
            reference_start: 10,
            reference_end: 675,
            query_start: 20,
            query_end: 685,
            reference_aligned: reference_row,
            query_aligned: query_row,
            cigar: "599=1I1=1D64=".into(),
            operations: vec![
                CigarRun {
                    kind: CigarKind::Match,
                    length: 599,
                },
                CigarRun {
                    kind: CigarKind::Insertion,
                    length: 1,
                },
                CigarRun {
                    kind: CigarKind::Match,
                    length: 1,
                },
                CigarRun {
                    kind: CigarKind::Deletion,
                    length: 1,
                },
                CigarRun {
                    kind: CigarKind::Match,
                    length: 64,
                },
            ],
            differences: vec![
                Difference {
                    kind: DifferenceKind::Insertion,
                    reference_start: 609,
                    reference_end: 609,
                    query_start: 619,
                    query_end: 620,
                    reference_bases: String::new(),
                    query_bases: "T".into(),
                },
                Difference {
                    kind: DifferenceKind::Deletion,
                    reference_start: 610,
                    reference_end: 611,
                    query_start: 621,
                    query_end: 621,
                    reference_bases: "C".into(),
                    query_bases: String::new(),
                },
            ],
            stats: AlignmentStats {
                columns: 666,
                paired_bases: 664,
                certain_matches: 664,
                ambiguous_compatible: 0,
                mismatches: 0,
                inserted_query_bases: 1,
                deleted_reference_bases: 1,
                reference_span_bases: 665,
                query_span_bases: 665,
                reference_coverage: 664.0 / 680.0,
                query_coverage: 664.0 / 692.0,
                certain_identity: Some(1.0),
                compatible_identity: Some(1.0),
            },
            warnings: vec![
                "Fixture records a declared path, without asserting DP optimality.".into(),
            ],
            reference_hash: blake3::hash(reference_source.as_bytes())
                .to_hex()
                .to_string(),
            query_hash: blake3::hash(query_source.as_bytes()).to_hex().to_string(),
        };
        let index = PageIndex::new(&alignment);
        assert_eq!(index.page_count(), 2);
        assert_eq!(index.difference_column(0), Some(599));
        assert_eq!(index.difference_column(1), Some(601));
        assert_eq!(index.difference_column(2), None);
        let next = render_page(&alignment, &index, 1);
        assert!(next.starts_with(&format!("R     609 GC{}\n", "T".repeat(58))));
        assert!(next.contains(&format!("Q     620 G-{}\n", "T".repeat(58))));
        assert!(next.contains("R     669 TTTTTT\n"));
        assert!(next.contains("Q     679 TTTTTT\n"));
        assert!(render_page(&alignment, &index, 0).starts_with("R      10 "));
        assert!(render_page(&alignment, &index, usize::MAX).is_empty());
    }

    #[test]
    fn columns_after_four_thousand_and_more_than_two_hundred_differences_remain_indexed() {
        // 250 isolated C→T substitutions at 19+20k; all 5,000 columns are
        // explicitly declared. This fixture needs no quadratic alignment run.
        let reference_row = "AC".repeat(2500);
        let query_row = "ACACACACACACACACACAT".repeat(250);
        let differences = (0..250)
            .map(|index| {
                let position = 19 + 20 * index;
                Difference {
                    kind: DifferenceKind::Substitution,
                    reference_start: position,
                    reference_end: position + 1,
                    query_start: position,
                    query_end: position + 1,
                    reference_bases: "C".into(),
                    query_bases: "T".into(),
                }
            })
            .collect();
        let operations = (0..250)
            .flat_map(|_| {
                [
                    CigarRun {
                        kind: CigarKind::Match,
                        length: 19,
                    },
                    CigarRun {
                        kind: CigarKind::Mismatch,
                        length: 1,
                    },
                ]
            })
            .collect();
        let alignment = Alignment {
            algorithm: "independently declared late-page presentation fixture".into(),
            options: AlignmentOptions::default(),
            score: 8750,
            reference_start: 0,
            reference_end: 5000,
            query_start: 0,
            query_end: 5000,
            reference_hash: blake3::hash(reference_row.as_bytes()).to_hex().to_string(),
            query_hash: blake3::hash(query_row.as_bytes()).to_hex().to_string(),
            reference_aligned: reference_row,
            query_aligned: query_row,
            cigar: "19=1X".repeat(250),
            operations,
            differences,
            stats: AlignmentStats {
                columns: 5000,
                paired_bases: 5000,
                certain_matches: 4750,
                ambiguous_compatible: 0,
                mismatches: 250,
                inserted_query_bases: 0,
                deleted_reference_bases: 0,
                reference_span_bases: 5000,
                query_span_bases: 5000,
                reference_coverage: 1.0,
                query_coverage: 1.0,
                certain_identity: Some(0.95),
                compatible_identity: Some(0.95),
            },
            warnings: vec![],
        };
        let index = PageIndex::new(&alignment);
        assert_eq!(index.page_count(), 9);
        assert_eq!(index.difference_column(200), Some(4019));
        assert_eq!(index.difference_column(249), Some(4999));
        assert_eq!(index.difference_column(250), None);
        let page = render_page(&alignment, &index, 6);
        assert!(page.starts_with("R    3600 "));
        assert!(page.contains(&format!("R    3960 {}\n", "AC".repeat(30))));
        assert!(page.contains(&format!("Q    3960 {}\n", "ACACACACACACACACACAT".repeat(3))));
        assert!(page.contains(&"|||||||||||||||||||×".repeat(3)));
        let last = render_page(&alignment, &index, 8);
        assert!(last.starts_with("R    4800 "));
        assert!(last.contains("Q    4980 ACACACACACACACACACAT\n"));
    }
}
