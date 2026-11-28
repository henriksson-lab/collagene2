//! Independently authored SCF2 fixtures and hand-specified reference/read evidence.
use collagene_core::review::reads::{trace_validation::*, traces};
use collagene_model::{Molecule, Topology};
fn scf(calls: &str, scores: &[u8], name: &str) -> Vec<u8> {
    assert_eq!(calls.len(), scores.len());
    let n = calls.len();
    let samples = n * 3 + 2;
    let bases_offset = 128 + samples * 8;
    let comment = format!("NAME={name}\n\0");
    let mut bytes = vec![0u8; 128];
    bytes[..4].copy_from_slice(b".scf");
    for (at, v) in [
        (4, samples),
        (8, 128),
        (12, n),
        (24, bases_offset),
        (28, comment.len()),
        (32, bases_offset + n * 12),
        (40, 2),
        (44, 2),
    ] {
        bytes[at..at + 4].copy_from_slice(&(v as u32).to_be_bytes());
    }
    bytes[36..40].copy_from_slice(b"2.02");
    bytes.resize(bases_offset, 0);
    for (i, b) in calls.bytes().enumerate() {
        bytes.extend_from_slice(&((i * 3 + 1) as u32).to_be_bytes());
        for base in b"ACGT" {
            bytes.push(if *base == b { scores[i] } else { 0 });
        }
        bytes.push(b);
        bytes.extend_from_slice(&[0; 3]);
    }
    bytes.extend_from_slice(comment.as_bytes());
    bytes
}
fn read(calls: &str, scores: &[u8], name: &str) -> Molecule {
    traces::import_trace(&scf(calls, scores, name))
        .unwrap()
        .attached_molecule()
        .unwrap()
}
fn refseq(s: &str, circular: bool) -> Molecule {
    Molecule::dna(
        "independent reference",
        s,
        if circular {
            Topology::Circular
        } else {
            Topology::Linear
        },
    )
    .unwrap()
}
fn options() -> TraceValidationOptions {
    TraceValidationOptions {
        trim: TrimOptions {
            mode: TrimMode::AllCalls,
            min_retained_bases: 1,
            ..Default::default()
        },
        ..Default::default()
    }
}

fn abif(calls: &str, quality: Option<&[u8]>) -> Vec<u8> {
    let samples = calls.len() * 3 + 2;
    let mut fields = vec![
        (*b"FWO_", 1u32, 2u16, 1u16, b"ACGT".to_vec()),
        (*b"PBAS", 2, 2, 1, calls.as_bytes().to_vec()),
        (
            *b"PLOC",
            2,
            4,
            2,
            (0..calls.len())
                .flat_map(|i| ((i * 3 + 1) as u16).to_be_bytes())
                .collect(),
        ),
    ];
    for tag in 9..=12 {
        fields.push((*b"DATA", tag, 4, 2, vec![0; samples * 2]));
    }
    if let Some(q) = quality {
        fields.push((*b"PCON", 2, 2, 1, q.to_vec()));
    }
    let offset = 128usize;
    let mut bytes = vec![0; offset + fields.len() * 28];
    bytes[..4].copy_from_slice(b"ABIF");
    bytes[4..6].copy_from_slice(&101u16.to_be_bytes());
    bytes[6..10].copy_from_slice(b"tdir");
    bytes[10..14].copy_from_slice(&1u32.to_be_bytes());
    bytes[14..16].copy_from_slice(&1023u16.to_be_bytes());
    bytes[16..18].copy_from_slice(&28u16.to_be_bytes());
    bytes[18..22].copy_from_slice(&(fields.len() as u32).to_be_bytes());
    bytes[22..26].copy_from_slice(&((fields.len() * 28) as u32).to_be_bytes());
    bytes[26..30].copy_from_slice(&(offset as u32).to_be_bytes());
    for (i, (name, number, kind, size, data)) in fields.into_iter().enumerate() {
        let at = offset + i * 28;
        bytes[at..at + 4].copy_from_slice(&name);
        bytes[at + 4..at + 8].copy_from_slice(&number.to_be_bytes());
        bytes[at + 8..at + 10].copy_from_slice(&kind.to_be_bytes());
        bytes[at + 10..at + 12].copy_from_slice(&size.to_be_bytes());
        bytes[at + 12..at + 16]
            .copy_from_slice(&((data.len() / size as usize) as u32).to_be_bytes());
        bytes[at + 16..at + 20].copy_from_slice(&(data.len() as u32).to_be_bytes());
        if data.len() <= 4 {
            bytes[at + 20..at + 20 + data.len()].copy_from_slice(&data);
        } else {
            let pos = bytes.len();
            bytes[at + 20..at + 24].copy_from_slice(&(pos as u32).to_be_bytes());
            bytes.extend(data);
        }
    }
    bytes
}

#[test]
fn abif_raw_scores_and_absent_scores_never_acquire_an_invented_phred_scale() {
    let reference = refseq("TTTCGATACGGG", false);
    let scored = traces::import_trace(&abif("CGATAC", Some(&[40; 6])))
        .unwrap()
        .attached_molecule()
        .unwrap();
    let report = validate_reads(&reference, &[scored], &options()).unwrap();
    assert_eq!(
        report.reads[0].trim.as_ref().unwrap().quality_source,
        QualitySource::AbifRawBasecaller
    );
    assert_eq!(report.eligible_coverage, 6);
    let unscored = traces::import_trace(&abif("CGATAC", None))
        .unwrap()
        .attached_molecule()
        .unwrap();
    let rejected = validate_reads(&reference, &[unscored.clone()], &options()).unwrap();
    assert!(
        rejected.reads[0]
            .error
            .as_ref()
            .unwrap()
            .contains("no native confidence")
    );
    let mut opts = options();
    opts.trim.missing_quality = MissingQualityPolicy::KeepUnscored;
    let kept = validate_reads(&reference, &[unscored], &opts).unwrap();
    assert!(kept.reads[0].placement_unique_proven);
    assert_eq!(kept.eligible_coverage, 0);
    assert!(
        kept.reads[0]
            .columns
            .iter()
            .all(|c| c.raw_confidence.is_none() && !c.eligible)
    );
    assert_eq!(kept.status, ValidationStatus::Incomplete);
}

fn all_paths(
    r: &[u8],
    q: &[u8],
    i: usize,
    j: usize,
    previous: u8,
    score: i32,
    best: &mut (i32, u8),
) {
    if i == r.len() && j == q.len() {
        if score > best.0 {
            *best = (score, 1);
        } else if score == best.0 {
            best.1 = best.1.saturating_add(1).min(2);
        }
        return;
    }
    if i < r.len() && j < q.len() {
        all_paths(
            r,
            q,
            i + 1,
            j + 1,
            0,
            score + if r[i] == q[j] { 2 } else { -3 },
            best,
        );
    }
    if i < r.len() {
        all_paths(
            r,
            q,
            i + 1,
            j,
            1,
            score + if previous == 1 { -1 } else { -4 },
            best,
        );
    }
    if j < q.len() {
        all_paths(
            r,
            q,
            i,
            j + 1,
            2,
            score + if previous == 2 { -1 } else { -4 },
            best,
        );
    }
}
fn binary(n: usize, bits: usize) -> String {
    (0..n)
        .map(|i| if bits & (1 << i) == 0 { 'A' } else { 'C' })
        .collect()
}
#[test]
fn optimal_alignment_and_uniqueness_match_independent_all_substring_path_enumeration() {
    let mut opts = options();
    opts.both_strands = false;
    for n in [3, 4] {
        for r in 0..1 << n {
            let reference = refseq(&binary(n, r), false);
            for m in [2, 3] {
                for q in 0..1 << m {
                    let query = binary(m, q);
                    let mut expected = (i32::MIN, 0u8);
                    for start in 0..=n {
                        for end in start..=n {
                            all_paths(
                                &reference.sequence.as_str().as_bytes()[start..end],
                                query.as_bytes(),
                                0,
                                0,
                                0,
                                0,
                                &mut expected,
                            );
                        }
                    }
                    let report =
                        validate_reads(&reference, &[read(&query, &vec![40; m], "oracle")], &opts)
                            .unwrap();
                    let actual = &report.reads[0];
                    assert!(actual.error.is_none(), "{query}: {:?}", actual.error);
                    assert_eq!(
                        actual.alignment.as_ref().unwrap().score,
                        expected.0,
                        "ref {} query {query}",
                        reference.sequence.as_str()
                    );
                    assert_eq!(
                        actual.optimal_paths_saturated,
                        expected.1,
                        "ref {} query {query}",
                        reference.sequence.as_str()
                    );
                    assert_eq!(actual.placement_unique_proven, expected.1 == 1);
                }
            }
        }
    }
}

#[test]
fn support_and_strand_requirements_apply_to_indels_as_well_as_base_calls() {
    let reference = refseq("TTTCGATACGGG", false);
    let insertion = read("CGAGTAC", &[40; 7], "insertion one file");
    let deletion = read("CGAAC", &[40; 5], "deletion one file");
    let mut opts = options();
    opts.min_support = 2;
    assert_eq!(
        validate_reads(&reference, &[insertion.clone()], &opts)
            .unwrap()
            .status,
        ValidationStatus::Incomplete
    );
    assert_eq!(
        validate_reads(&reference, &[deletion.clone()], &opts)
            .unwrap()
            .status,
        ValidationStatus::Incomplete
    );
    opts.min_support = 1;
    opts.require_both_strands = true;
    assert_eq!(
        validate_reads(&reference, &[insertion], &opts)
            .unwrap()
            .status,
        ValidationStatus::Incomplete
    );
    assert_eq!(
        validate_reads(&reference, &[deletion], &opts)
            .unwrap()
            .status,
        ValidationStatus::Incomplete
    );
}

#[test]
fn explicit_fastq_encodings_support_calls_without_inventing_chromatogram_coordinates() {
    use collagene_core::formats::oligos::{self as oligo_formats, QualityEncoding};
    let reference = refseq("TTTCGATACGGG", false);
    for (encoding, quality, kind) in [
        (
            QualityEncoding::Phred33,
            "IIIIII",
            QualitySource::FastqPhred33,
        ),
        (
            QualityEncoding::Phred64,
            "hhhhhh",
            QualitySource::FastqPhred64,
        ),
    ] {
        let bytes = format!("@caller evidence\nCGATAC\n+\n{quality}\n");
        let source = oligo_formats::import_fastq(bytes.as_bytes(), encoding)
            .unwrap()
            .molecules
            .remove(0);
        let report = validate_reads(&reference, &[source.clone()], &options()).unwrap();
        let read = &report.reads[0];
        assert_eq!(read.read_format, Some(ReadFormat::Fastq));
        assert_eq!(read.trace_orientation, None);
        assert_eq!(read.trim.as_ref().unwrap().quality_source, kind);
        assert_eq!(report.consensus, "NNNCGATACNNN");
        assert!(read.columns.iter().all(|c| c.sample_index.is_none()
            && c.stored_sample_index.is_none()
            && c.raw_confidence == Some(40)));
        let reverse = collagene_core::sequence::annotations::reverse_complement(&source).unwrap();
        let transformed = validate_reads(&reference, &[reverse], &options()).unwrap();
        assert_eq!(transformed.consensus, report.consensus);
        assert!(transformed.reads[0].query_reversed);
        assert!(!transformed.reads[0].original_read_reverse_to_reference);
        assert_eq!(transformed.reads[0].columns[0].stored_call_index, Some(0));
        let stale = source.replace(0, 1, "A").unwrap();
        let stale = validate_reads(&reference, &[stale], &options()).unwrap();
        assert!(stale.reads[0].error.as_ref().unwrap().contains("FASTQ"));
        assert_eq!(stale.eligible_coverage, 0);
    }
}

#[test]
fn trimmed_fastq_coordinates_remain_in_original_forward_record_after_reversal() {
    use collagene_core::formats::oligos::{self as oligo_formats, QualityEncoding};
    let source = oligo_formats::import_fastq(
        b"@asymmetric\nAACGATACG\n+\nIIIIIIIII\n",
        QualityEncoding::Phred33,
    )
    .unwrap()
    .molecules
    .remove(0);
    let trimmed = oligo_formats::trim_fastq(&source, 2, 8).unwrap();
    let reversed = collagene_core::sequence::annotations::reverse_complement(&trimmed).unwrap();
    for read in [trimmed, reversed] {
        let report = validate_reads(&refseq("TTTCGATACGGG", false), &[read], &options()).unwrap();
        let evidence = &report.reads[0];
        assert_eq!(report.consensus, "NNNCGATACNNN");
        assert_eq!(evidence.source_call_length, Some(9));
        assert_eq!(evidence.source_retained_start, Some(2));
        assert_eq!(evidence.source_retained_end, Some(8));
        assert_eq!(
            evidence
                .columns
                .iter()
                .map(|c| c.stored_call_index.unwrap())
                .collect::<Vec<_>>(),
            (2..8).collect::<Vec<_>>()
        );
        assert!(
            evidence
                .columns
                .iter()
                .all(|c| c.sample_index.is_none() && c.stored_sample_index.is_none())
        );
        assert!(!evidence.original_read_reverse_to_reference);
    }
}

#[test]
fn circular_unique_claims_are_sound_against_all_one_turn_substring_paths() {
    let mut opts = options();
    opts.both_strands = false;
    for n in [3, 4] {
        for r in 0..1 << n {
            let reference = refseq(&binary(n, r), true);
            for m in [2, 3] {
                for q in 0..1 << m {
                    let query = binary(m, q);
                    let mut expected = (i32::MIN, 0u8);
                    for start in 0..n {
                        for length in 0..=n {
                            let interval = (0..length)
                                .map(|i| reference.sequence.as_str().as_bytes()[(start + i) % n])
                                .collect::<Vec<_>>();
                            all_paths(&interval, query.as_bytes(), 0, 0, 0, 0, &mut expected);
                        }
                    }
                    let report = validate_reads(
                        &reference,
                        &[read(&query, &vec![40; m], "circle oracle")],
                        &opts,
                    )
                    .unwrap();
                    let actual = &report.reads[0];
                    if actual.error.is_none() {
                        assert_eq!(actual.alignment.as_ref().unwrap().score, expected.0);
                        if actual.placement_unique_proven {
                            assert_eq!(expected.1, 1, "{} / {query}", reference.sequence.as_str());
                        }
                    } else {
                        assert_eq!(report.status, ValidationStatus::Incomplete);
                        assert_eq!(report.eligible_coverage, 0);
                    }
                }
            }
        }
    }
}
#[test]
fn quality_window_trimming_retains_internal_weak_calls_as_evidence_without_phred() {
    let trace =
        traces::import_trace(&scf("ACGATTCG", &[0, 30, 40, 0, 40, 30, 0, 0], "trim")).unwrap();
    let options = TrimOptions {
        mode: TrimMode::QualityWindows,
        window_size: 3,
        minimum_passing: 2,
        min_retained_bases: 1,
        ..Default::default()
    };
    let report = trim_trace(&trace, &options).unwrap();
    assert_eq!((report.start, report.end), (0, 7));
    assert_eq!(report.low_or_unknown_calls, vec![0, 3, 6]);
    assert_eq!(report.passing_calls, 4);
    assert_eq!(report.quality_source, QualitySource::ScfCalledBaseAccuracy);
    assert!(report.warnings[0].contains("not calibrated Phred"));
    let uncertain = traces::import_trace(&scf("AN", &[40, 40], "uncertain")).unwrap();
    let report = trim_trace(
        &uncertain,
        &TrimOptions {
            mode: TrimMode::AllCalls,
            min_retained_bases: 1,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(report.low_or_unknown_calls, vec![1]);
}
#[test]
fn partial_exact_read_does_not_fill_uncovered_reference_or_claim_whole_plasmid_validation() {
    let reference = refseq("TTTCGATACGGG", false);
    let source = read("CGATAC", &[40; 6], "forward");
    let before = source.revision();
    let report = validate_reads(&reference, &[source.clone()], &options()).unwrap();
    assert_eq!(report.status, ValidationStatus::Incomplete);
    assert_eq!(report.consensus, "NNNCGATACNNN");
    assert_eq!(report.uncovered_ranges, vec![(0, 3), (9, 12)]);
    assert_eq!(report.eligible_coverage, 6);
    assert!(report.reads[0].placement_unique_proven);
    assert_eq!(
        report.reads[0]
            .columns
            .iter()
            .map(|c| c.reference_position.unwrap())
            .collect::<Vec<_>>(),
        vec![3, 4, 5, 6, 7, 8]
    );
    assert_eq!(source.revision(), before);
}
#[test]
fn raw_low_quality_substitution_is_visible_but_not_a_supported_discrepancy() {
    let reference = refseq("TTTCGATACGGG", false);
    let source = read("CGCTAC", &[40, 40, 0, 40, 40, 40], "low variant");
    let report = validate_reads(&reference, &[source], &options()).unwrap();
    assert!(report.reads[0].placement_unique_proven);
    assert_eq!(report.status, ValidationStatus::Incomplete);
    let column = &report.reads[0].columns[2];
    assert_eq!(column.observed_base, Some('C'));
    assert_eq!(column.reference_base, Some('A'));
    assert!(!column.eligible);
    assert_eq!(column.raw_confidence, Some(0));
    assert_eq!(report.consensus, "NNNCGNTACNNN");
}
#[test]
fn high_native_score_unique_mismatch_is_reconstructed_at_source_coordinates() {
    let reference = refseq("TTTCGATACGGG", false);
    let source = read("CGCTAC", &[40; 6], "variant");
    let report = validate_reads(&reference, &[source], &options()).unwrap();
    assert!(report.reads[0].placement_unique_proven);
    assert_eq!(report.status, ValidationStatus::Discrepancies);
    assert_eq!(report.consensus, "NNNCGCTACNNN");
    let site = report.sites.iter().find(|s| s.position == 5).unwrap();
    assert_eq!(site.state, SiteState::SupportedDifference);
    assert_eq!(site.counts, [0, 1, 0, 0]);
    assert_eq!(report.reads[0].alignment.as_ref().unwrap().cigar, "2=1X3=");
}
#[test]
fn reverse_alignment_and_reverse_display_keep_original_call_and_sample_coordinates() {
    let reference = refseq("TTTCGATACGGG", false);
    let source = read("GTATCG", &[40; 6], "reverse");
    let report = validate_reads(&reference, &[source.clone()], &options()).unwrap();
    let aligned = &report.reads[0];
    assert!(aligned.query_reversed);
    assert!(aligned.original_read_reverse_to_reference);
    assert_eq!(aligned.columns[0].call_index, Some(5));
    assert_eq!(aligned.columns[0].sample_index, Some(16));
    let view = traces::trace_from_molecule(&source)
        .unwrap()
        .reverse_complement()
        .unwrap()
        .attached_molecule()
        .unwrap();
    let transformed = validate_reads(&reference, &[view], &options()).unwrap();
    let aligned = &transformed.reads[0];
    assert!(!aligned.query_reversed);
    assert!(aligned.original_read_reverse_to_reference);
    assert_eq!(aligned.columns[0].call_index, Some(0));
    assert_eq!(aligned.columns[0].stored_call_index, Some(5));
    assert_eq!(aligned.columns[0].sample_index, Some(3));
    assert_eq!(aligned.columns[0].stored_sample_index, Some(16));
    assert_eq!(transformed.consensus, report.consensus);
}
#[test]
fn repeated_and_palindromic_placements_do_not_become_arbitrary_consensus_coverage() {
    let reference = refseq("CGATACAAACGATAC", false);
    let report = validate_reads(
        &reference,
        &[read("CGATAC", &[40; 6], "repeat")],
        &options(),
    )
    .unwrap();
    assert_eq!(report.reads[0].exact_placement_count, 2);
    assert!(!report.reads[0].placement_unique_proven);
    assert_eq!(report.eligible_coverage, 0);
    assert_eq!(report.consensus, "N".repeat(15));
    let reference = refseq("TTTACGTGGG", false);
    let source = read("ACGT", &[40; 4], "palindrome");
    let report = validate_reads(&reference, &[source.clone()], &options()).unwrap();
    assert!(!report.reads[0].placement_unique_proven);
    let mut opts = options();
    opts.both_strands = false;
    assert!(validate_reads(&reference, &[source], &opts).unwrap().reads[0].placement_unique_proven);
}
#[test]
fn circular_origin_maps_each_base_once_and_retains_unknown_gaps() {
    let reference = refseq("ACGATTCG", true);
    let report = validate_reads(
        &reference,
        &[read("CGACGA", &[40; 6], "origin")],
        &options(),
    )
    .unwrap();
    assert!(report.reads[0].placement_unique_proven);
    assert_eq!(
        report.reads[0]
            .columns
            .iter()
            .map(|c| c.reference_position.unwrap())
            .collect::<Vec<_>>(),
        vec![6, 7, 0, 1, 2, 3]
    );
    assert_eq!(report.consensus, "ACGANNCG");
    assert_eq!(report.uncovered_ranges, vec![(4, 6)]);
}
#[test]
fn overlapping_unique_reads_can_cover_reference_without_filling_it_from_target() {
    let reference = refseq("ACGATTCG", false);
    let report = validate_reads(
        &reference,
        &[
            read("ACGA", &[40; 4], "left"),
            read("TTCG", &[40; 4], "right"),
        ],
        &options(),
    )
    .unwrap();
    assert_eq!(
        report.status,
        ValidationStatus::ConsistentUnderRecordedThresholds
    );
    assert_eq!(report.consensus, "ACGATTCG");
    assert_eq!(report.reference_coverage_fraction, 1.);
    assert!(report.uncovered_ranges.is_empty());
}
#[test]
fn conflicting_eligible_calls_remain_iupac_disagreement_not_reference_or_majority_vote() {
    let reference = refseq("ACGATTCG", false);
    let report = validate_reads(
        &reference,
        &[
            read("ACGA", &[40; 4], "reference call"),
            read("ACCA", &[40; 4], "discordant call"),
            read("TTCG", &[40; 4], "right"),
        ],
        &options(),
    )
    .unwrap();
    assert_eq!(report.status, ValidationStatus::Discrepancies);
    assert_eq!(report.consensus, "ACSATTCG");
    assert_eq!(
        report.sites.iter().find(|s| s.position == 2).unwrap().state,
        SiteState::ConflictingCalls
    );
}
#[test]
fn duplicate_original_attachment_is_not_independent_support_after_rename() {
    let reference = refseq("CGATAC", false);
    let source = read("CGATAC", &[40; 6], "same attachment");
    let mut renamed = source.clone();
    renamed.name = "renamed copy".into();
    let mut opts = options();
    opts.min_support = 2;
    let report = validate_reads(&reference, &[source, renamed], &opts).unwrap();
    assert!(
        report.reads[1]
            .error
            .as_ref()
            .unwrap()
            .contains("duplicate original")
    );
    assert_eq!(report.status, ValidationStatus::Incomplete);
    assert_eq!(report.eligible_coverage, 0);
}
#[test]
fn explicit_strand_and_minimum_support_policies_are_respected() {
    let reference = refseq("CGATAC", false);
    let mut opts = options();
    opts.min_support = 2;
    opts.require_both_strands = true;
    let forward = read("CGATAC", &[40; 6], "F1");
    let another = read("CGATAC", &[40; 6], "F2");
    let reverse = read("GTATCG", &[40; 6], "R1");
    let same = validate_reads(&reference, &[forward.clone(), another], &opts).unwrap();
    assert_eq!(same.status, ValidationStatus::Incomplete);
    let both = validate_reads(&reference, &[forward, reverse], &opts).unwrap();
    assert_eq!(
        both.status,
        ValidationStatus::ConsistentUnderRecordedThresholds
    );
    assert!(
        both.sites
            .iter()
            .all(|s| s.forward_reads == 1 && s.reverse_reads == 1)
    );
}
#[test]
fn unique_indel_events_keep_read_and_reference_boundaries() {
    let reference = refseq("TTTCGATACGGG", false);
    let insertion = validate_reads(
        &reference,
        &[read("CGAGTAC", &[40; 7], "insertion")],
        &options(),
    )
    .unwrap();
    assert!(insertion.reads[0].placement_unique_proven);
    assert_eq!(insertion.status, ValidationStatus::Discrepancies);
    assert_eq!(insertion.insertions.len(), 1);
    assert_eq!(insertion.insertions[0].reference_boundary, 6);
    assert_eq!(insertion.insertions[0].sequence, "G");
    assert!(insertion.insertions[0].eligible);
    let deletion = validate_reads(
        &reference,
        &[read("CGAAC", &[40; 5], "deletion")],
        &options(),
    )
    .unwrap();
    assert!(deletion.reads[0].placement_unique_proven);
    assert_eq!(
        deletion
            .sites
            .iter()
            .find(|s| s.position == 6)
            .unwrap()
            .state,
        SiteState::DeletionEvidence
    );
    assert_eq!(deletion.status, ValidationStatus::Discrepancies);
}
#[test]
fn stale_calls_budget_limits_and_invalid_options_remain_unknown_not_negative_evidence() {
    let reference = refseq("TTTCGATACGGG", false);
    let source = read("CGATAC", &[40; 6], "source");
    let stale = source.replace(0, 1, "A").unwrap();
    let report = validate_reads(&reference, &[stale], &options()).unwrap();
    assert!(
        report.reads[0]
            .error
            .as_ref()
            .unwrap()
            .contains("differs from")
    );
    assert_eq!(report.eligible_coverage, 0);
    let mut opts = options();
    opts.max_exact_comparisons = 1;
    opts.max_total_cells = 10;
    let report = validate_reads(&reference, &[source], &opts).unwrap();
    assert!(report.reads[0].error.as_ref().unwrap().contains("budget"));
    assert_eq!(report.status, ValidationStatus::Incomplete);
    // A negative ambiguity score disables the alignment identity fast path;
    // exact reads must then charge their actual DP plus proof work as well.
    opts.max_exact_comparisons = 50_000_000;
    opts.scoring.ambiguous_score = -1;
    let report = validate_reads(
        &reference,
        &[read("CGATAC", &[40; 6], "negative ambiguity")],
        &opts,
    )
    .unwrap();
    assert!(report.reads[0].error.as_ref().unwrap().contains("budget"));
    assert_eq!(report.eligible_coverage, 0);
    opts.trim.window_size = 0;
    assert!(validate_reads(&reference, &[read("CGATAC", &[40; 6], "invalid")], &opts).is_err());
}
