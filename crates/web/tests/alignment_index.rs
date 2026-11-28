#[path = "../src/alignment_index.rs"]
mod alignment_index;
use alignment_index::PreparedAlignment;
use collagene_core::alignment_formats::{AlignmentFormat, import_alignment};
use collagene_model::Alphabet;
#[test]
fn every_gap_pattern_has_literal_residue_coordinates() {
    for mask in 0..256 {
        let row = (0..8)
            .map(|i| if mask & (1 << i) == 0 { '-' } else { 'A' })
            .collect::<String>();
        let doc = import_alignment(
            format!(">row\n{row}\n").as_bytes(),
            AlignmentFormat::AlignedFasta,
            Alphabet::Dna,
            "gaps",
        )
        .unwrap();
        let p = PreparedAlignment::new(doc).unwrap();
        let mut ordinal = 0;
        for (i, c) in row.bytes().enumerate() {
            if c == b'-' {
                assert_eq!(p.position(0, i).unwrap(), None);
            } else {
                assert_eq!(p.position(0, i).unwrap(), Some(ordinal));
                ordinal += 1;
            }
        }
        assert_eq!(p.residue_count(0), Some(ordinal));
        assert_eq!(p.document().rows[0].aligned, row);
    }
}
#[test]
fn block_boundaries_preserve_counts_and_reject_invalid_coordinates() {
    let row = "--ACG".repeat(200);
    let doc = import_alignment(
        format!(">one\n{row}\n>two\n{}\n", "-".repeat(row.len())).as_bytes(),
        AlignmentFormat::AlignedFasta,
        Alphabet::Dna,
        "blocks",
    )
    .unwrap();
    let p = PreparedAlignment::new(doc).unwrap();
    for column in [0, 1, 2, 255, 256, 257, 511, 512, 768, 999] {
        let expected = if row.as_bytes()[column] == b'-' {
            None
        } else {
            Some(row[..column].bytes().filter(|c| *c != b'-').count())
        };
        assert_eq!(p.position(0, column).unwrap(), expected);
        assert_eq!(p.position(1, column).unwrap(), None);
    }
    assert_eq!(p.residue_count(0), Some(600));
    assert_eq!(p.residue_count(1), Some(0));
    assert!(p.position(2, 0).is_err());
    assert!(p.position(0, 1000).is_err());
}
#[test]
fn capped_scroll_still_reaches_the_last_of_eight_million_columns() {
    use alignment_index::{logical_left, native_left, spacer_width};
    let columns = 8_000_000;
    let viewport = 1000.;
    assert_eq!(spacer_width(columns), 8_000_000.);
    let at_end = logical_left(spacer_width(columns) - viewport, columns, viewport);
    assert_eq!(at_end, 79_999_180.);
    // Last base starts at79,999,990 and ends at80,000,000; visible sequence
    // width is820px because the immutable row-label gutter occupies180px.
    assert!(79_999_990. >= at_end && 80_000_000. <= at_end + 820.);
    assert_eq!(native_left(79_999_990., columns, viewport), 7_999_000.);
    assert_eq!(logical_left(125., 1000, viewport), 125.);
    assert_eq!(native_left(125., 1000, viewport), 125.);
    assert_eq!(logical_left(10., 5, viewport), 0.);
}
