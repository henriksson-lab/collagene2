//! Independent literal SAM coordinates and calls; no aligner output is an oracle.
use collagene_core::read_alignment_reference::*;
use collagene_model::{Alphabet, Molecule, Sequence, Topology};
use collagene_read_alignments::{Header, Limits, Record, SamReader};
use std::io::Cursor;

fn sam(header: &str, line: &str) -> (Header, Record) {
    let mut reader = SamReader::new(
        Cursor::new(format!("{header}{line}").into_bytes()),
        Limits::default(),
    )
    .unwrap();
    let header = reader.header().clone();
    let record = reader.next_record().unwrap().unwrap();
    (header, record)
}
fn reference(sequence: &str, topology: Topology) -> Molecule {
    Molecule::dna("explicit actual source", sequence, topology).unwrap()
}
fn options(reference: &Molecule) -> Options {
    Options {
        reference_name: "ref".into(),
        expected_reference_revision: reference.revision(),
        ..Options::default()
    }
}

#[test]
fn actual_calls_check_cigar_assertions_and_keep_insertions_deletions_skips_and_clips_distinct() {
    // Reference offsets1/2/3/4/8 contain C/G/T/A/A; stored aligned calls
    // C/T/A/A/A therefore have three agreements and two disagreements.
    let (header, record) = sam(
        "@SQ\tSN:ref\tLN:12\n",
        "r\t0\tref\t2\t60\t1S2M1I1=1X1D2N1P1M1S\t*\t0\t0\tTCTGAAAG\t*\n",
    );
    let source = reference("ACGTACGTACGT", Topology::Linear);
    let report = check_record(&header, &record, &source, &options(&source)).unwrap();
    assert!(report.complete && report.reference_bound);
    assert_eq!(report.binding.reference_revision, source.revision());
    assert_eq!(
        report
            .columns
            .iter()
            .map(|c| c.reference_position)
            .collect::<Vec<_>>(),
        vec![
            None,
            Some(1),
            Some(2),
            None,
            Some(3),
            Some(4),
            Some(8),
            None
        ]
    );
    assert_eq!(report.summary.checked_aligned_calls, 5);
    assert_eq!(report.summary.certain_agreements, 3);
    assert_eq!(report.summary.certain_disagreements, 2);
    assert_eq!(report.summary.cigar_contradictions, 2);
    assert_eq!(report.summary.inserted_calls, 1);
    assert_eq!(report.summary.soft_clipped_calls, 2);
    assert_eq!(report.summary.deleted_reference_bases, 1);
    assert_eq!(report.summary.skipped_reference_bases, 2);
    assert_eq!(report.summary.padding_units, 1);
    assert_eq!(
        report.runs.iter().map(|r| r.operation).collect::<Vec<_>>(),
        vec![
            Operation::SoftClip,
            Operation::Match,
            Operation::Insertion,
            Operation::Equal,
            Operation::Difference,
            Operation::Deletion,
            Operation::Skip,
            Operation::Pad,
            Operation::Match,
            Operation::SoftClip
        ]
    );
    for run in &report.runs {
        if matches!(
            run.operation,
            Operation::Deletion | Operation::Skip | Operation::Pad
        ) {
            assert_eq!(run.column_start, run.column_end);
        }
    }
}

#[test]
fn reverse_hard_clip_indexes_and_reference_equal_resolution_preserve_original_query_evidence() {
    let (header, record) = sam(
        "@SQ\tSN:ref\tLN:8\n",
        "r\t16\tref\t1\t30\t2H1S2M1I2M1H\t*\t0\t0\tG=R=GN\t*\n",
    );
    let source = reference("ACGTNACG", Topology::Linear);
    let original = record.original_bytes().to_vec();
    let report = check_record(&header, &record, &source, &options(&source)).unwrap();
    assert!(report.reference_bound && !report.complete && report.reverse);
    assert_eq!(report.implied_original_length, Some(9));
    assert_eq!(
        report
            .columns
            .iter()
            .map(|c| c.implied_original_query)
            .collect::<Vec<_>>(),
        vec![Some(6), Some(5), Some(4), Some(3), Some(2), Some(1)]
    );
    assert_eq!(report.columns[1].reported_symbol, '=');
    assert_eq!(report.columns[1].query_symbol, Some('A'));
    assert_eq!(
        report.columns[1].call_origin,
        CallOrigin::ReferenceEqualResolved
    );
    assert_eq!(report.columns[1].relation, Relation::ReferenceEqualResolved);
    assert_eq!(report.columns[3].query_symbol, None); // '=' inside insertion has no reference call.
    assert!(
        report
            .incomplete
            .contains(&Incomplete::ReferenceEqualityOutsideAlignment)
    );
    assert_eq!(report.summary.reference_equal_resolved, 1);
    assert_eq!(report.summary.certain_agreements, 1);
    assert_eq!(report.summary.certain_disagreements, 1);
    assert_eq!(report.summary.possible_agreements, 1);
    assert_eq!(report.summary.hard_clipped_absent_bases, 3);
    assert_eq!(record.original_bytes(), original);
    assert_eq!(record.fields().sequence.as_deref(), Some("G=R=GN"));
    assert!(record.dna_calls().is_err());
}

#[test]
fn circular_reference_wrap_and_iupac_uncertainty_do_not_become_exact_matches() {
    let (header, record) = sam(
        "@SQ\tSN:ref\tLN:4\tTP:circular\n",
        "r\t0\tref\t4\t30\t3M\t*\t0\t0\tNAT\t*\n",
    );
    let source = reference("ACGN", Topology::Circular);
    let report = check_record(&header, &record, &source, &options(&source)).unwrap();
    assert!(report.complete);
    assert_eq!(
        report
            .columns
            .iter()
            .map(|c| c.reference_position)
            .collect::<Vec<_>>(),
        vec![Some(3), Some(0), Some(1)]
    );
    assert_eq!(
        report
            .columns
            .iter()
            .map(|c| c.unwrapped_reference_position)
            .collect::<Vec<_>>(),
        vec![Some(3), Some(4), Some(5)]
    );
    assert_eq!(report.summary.possible_agreements, 1);
    assert_eq!(report.summary.certain_agreements, 1);
    assert_eq!(report.summary.certain_disagreements, 1);
    let wrong_topology = reference("ACGN", Topology::Linear);
    let withheld =
        check_record(&header, &record, &wrong_topology, &options(&wrong_topology)).unwrap();
    assert!(!withheld.complete && !withheld.reference_bound && withheld.columns.is_empty());
    assert!(
        withheld
            .incomplete
            .contains(&Incomplete::ReferenceTopologyMismatch)
    );
}

#[test]
fn explicit_revision_mapping_sealed_headers_missing_calls_unmapped_and_bounds_remain_required() {
    let text = "@SQ\tSN:ref\tLN:4\tM5:00000000000000000000000000000000\tUR:https://invalid.example/reference\n";
    let (header, record) = sam(text, "r\t0\tref\t1\t30\t2M\t*\t0\t0\tAC\t*\n");
    let source = reference("ACGT", Topology::Linear);
    let missing = check_record(&header, &record, &source, &Options::default()).unwrap();
    assert!(
        missing
            .incomplete
            .contains(&Incomplete::MissingExplicitReferenceMapping)
    );
    assert!(missing.columns.is_empty());
    let foreign = Header::from_sam(
        format!("{text}@CO\tdifferent immutable context\n").as_bytes(),
        &Limits::default(),
    )
    .unwrap();
    assert!(check_record(&foreign, &record, &source, &options(&source)).is_err());
    let actual_other = reference("TTTT", Topology::Linear);
    let wrong = check_record(&header, &record, &actual_other, &options(&source)).unwrap();
    assert!(
        wrong
            .incomplete
            .contains(&Incomplete::ReferenceRevisionMismatch)
    );
    let mut protein = source.clone();
    protein.sequence = Sequence::new(Alphabet::Protein, "ACDE").unwrap();
    let unsupported = check_record(&header, &record, &protein, &options(&protein)).unwrap();
    assert!(
        unsupported
            .incomplete
            .contains(&Incomplete::UnsupportedReferenceAlphabet)
    );
    let (header, absent) = sam(text, "r\t0\tref\t1\t30\t2M\t*\t0\t0\t*\t*\n");
    let missing = check_record(&header, &absent, &source, &options(&source)).unwrap();
    assert!(missing.reference_bound && !missing.complete && missing.columns.is_empty());
    assert!(missing.incomplete.contains(&Incomplete::MissingSequence));
    let (header, unmapped) = sam(text, "r\t4\tref\t1\t30\t2M\t*\t0\t0\tAC\t*\n");
    let missing = check_record(&header, &unmapped, &source, &options(&source)).unwrap();
    assert!(missing.incomplete.contains(&Incomplete::Unmapped));
    assert!(missing.columns.is_empty());
    let (conflict_header, conflict) = sam(text, "r\t0\tref\t1\t30\t1X\t*\t0\t0\t=\t*\n");
    let conflict = check_record(&conflict_header, &conflict, &source, &options(&source)).unwrap();
    assert!(!conflict.complete);
    assert!(
        conflict
            .incomplete
            .contains(&Incomplete::ConflictingImportedAssertions)
    );
    assert_eq!(conflict.summary.reference_equal_resolved, 1);
    assert_eq!(conflict.summary.certain_agreements, 0);
    let mut small = options(&source);
    small.max_columns = 1;
    assert!(check_record(&header, &unmapped, &source, &small).is_err());
}
