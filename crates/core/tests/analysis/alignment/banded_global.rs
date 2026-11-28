//! Current-kernel row-band integration fixtures, registered with the analysis harness.
//! Literal edits define the expected differences, not another alignment engine.
use collagene_core::alignment::{
    self, AlignmentMode, AlignmentOptions, ComparisonOptions, DifferenceKind, grid::GlobalStorage,
};
use collagene_model::{Molecule, Sequence, Topology};

fn literal_backbone(length: usize) -> String {
    let mut state = 0x519b_d28a_671c_e043u64;
    (0..length)
        .map(|_| {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            b"ACGT"[((state >> 61) & 3) as usize] as char
        })
        .collect()
}
fn options() -> AlignmentOptions {
    AlignmentOptions {
        mode: AlignmentMode::Global,
        global_storage: GlobalStorage::Adaptive {
            initial_radius: 16,
            maximum_radius: 4096,
            max_attempts: 10,
        },
        ..AlignmentOptions::default()
    }
}
#[test]
fn ten_kilobase_literal_substitution_is_certified_without_a_full_rectangle() {
    let reference = literal_backbone(10_000);
    let position = 4321;
    let mut changed = reference.as_bytes().to_vec();
    changed[position] = match changed[position] {
        b'A' => b'C',
        _ => b'A',
    };
    let query = String::from_utf8(changed).unwrap();
    let result = alignment::align(
        &Sequence::dna(&reference).unwrap(),
        &Sequence::dna(&query).unwrap(),
        &options(),
    )
    .unwrap();
    assert_eq!(result.score, 19_995);
    assert!(result.search.unrestricted_optimum_proven);
    assert!(result.search.evaluated_cells < 2_000_000);
    assert_eq!(result.differences.len(), 1);
    let difference = &result.differences[0];
    assert_eq!(difference.kind, DifferenceKind::Substitution);
    assert_eq!(
        (difference.reference_start, difference.reference_end),
        (position, position + 1)
    );
    assert_eq!(
        difference.reference_bases,
        reference[position..position + 1]
    );
    assert_eq!(difference.query_bases, query[position..position + 1]);
    // A full rectangle cannot fit the same admitted cell budget.
    let full = AlignmentOptions {
        global_storage: GlobalStorage::Full,
        ..options()
    };
    assert!(
        alignment::align(
            &Sequence::dna(&reference).unwrap(),
            &Sequence::dna(&query).unwrap(),
            &full
        )
        .is_err()
    );
}
#[test]
fn twenty_kilobase_three_base_deletion_preserves_original_coordinates_and_affine_cost() {
    let reference = literal_backbone(20_000);
    // Avoid the immediate period-three boundary slippage that would make this
    // constructed deletion's exact placement nonunique under identical scoring.
    let position = (10_000..10_100)
        .find(|&p| {
            reference.as_bytes()[p - 1] != reference.as_bytes()[p + 2]
                && reference.as_bytes()[p] != reference.as_bytes()[p + 3]
        })
        .unwrap();
    let query = format!("{}{}", &reference[..position], &reference[position + 3..]);
    let result = alignment::align(
        &Sequence::dna(&reference).unwrap(),
        &Sequence::dna(&query).unwrap(),
        &options(),
    )
    .unwrap();
    assert_eq!(result.score, 39_988); // 19,997 matches, one opening and two extensions.
    assert!(result.search.unrestricted_optimum_proven);
    assert!(result.search.evaluated_cells < 2_000_000);
    assert_eq!(result.differences.len(), 1);
    let difference = &result.differences[0];
    assert_eq!(difference.kind, DifferenceKind::Deletion);
    assert_eq!(
        (difference.reference_start, difference.reference_end),
        (position, position + 3)
    );
    assert_eq!(
        (difference.query_start, difference.query_end),
        (position, position)
    );
    assert_eq!(
        difference.reference_bases,
        reference[position..position + 3]
    );
    assert!(difference.query_bases.is_empty());
}
#[test]
fn twenty_kilobase_reverse_rotated_variant_is_not_starved_by_wrong_initial_candidates() {
    let source = literal_backbone(20_000);
    let origin = 4093;
    let position = 12_345;
    let mut changed = source.as_bytes().to_vec();
    changed[position] = match changed[position] {
        b'T' => b'G',
        _ => b'T',
    };
    let changed = String::from_utf8(changed).unwrap();
    let rotated = format!("{}{}", &changed[origin..], &changed[..origin]);
    let query = Sequence::dna(&rotated)
        .unwrap()
        .reverse_complement()
        .unwrap();
    let design = Molecule::dna(
        "separately bound circular design",
        &source,
        Topology::Circular,
    )
    .unwrap();
    let consensus = Molecule::dna(
        "supplied reverse consensus",
        query.as_str(),
        Topology::Linear,
    )
    .unwrap();
    let compared = alignment::compare_molecules(
        &design,
        &consensus,
        &ComparisonOptions {
            alignment: options(),
            both_query_strands: true,
            reference_rotations: true,
            max_origins: 16,
            max_total_cells: 8_000_000,
        },
    )
    .unwrap();
    assert_eq!(compared.reference_origin, origin);
    assert!(compared.query_reversed);
    assert_eq!(compared.alignment.score, 39_995);
    assert!(compared.alignment.search.unrestricted_optimum_proven);
    assert_eq!(compared.alignment.differences.len(), 1);
    let difference = &compared.alignment.differences[0];
    assert_eq!(difference.kind, DifferenceKind::Substitution);
    assert_eq!(
        difference.reference_start,
        (position + source.len() - origin) % source.len()
    );
    // A certified tested path does not certify the seed-limited circular domain
    // or erase the separately retained physical-topology declaration mismatch.
    assert!(!compared.exhaustive_reference_origins);
    assert!(!compared.optimality_proven);
    assert_eq!(compared.reference_topology, Topology::Circular);
    assert_eq!(compared.query_topology, Topology::Linear);
}

#[test]
fn cumulative_retry_refusal_retains_actual_path_without_claiming_global_optimum() {
    use alignment::grid::StopReason;
    let scoring = alignment::Scoring {
        gap_open: -1,
        gap_extend: -1,
        ..Default::default()
    };
    let options = AlignmentOptions {
        mode: AlignmentMode::Global,
        scoring,
        max_cells: 6,
        global_storage: GlobalStorage::Adaptive {
            initial_radius: 0,
            maximum_radius: 4,
            max_attempts: 8,
        },
    };
    let reference = Sequence::dna("AAAA").unwrap();
    let query = Sequence::dna("CCCC").unwrap();
    let limited = alignment::align(&reference, &query, &options).unwrap();
    assert_eq!(limited.score, -12);
    assert_eq!(limited.search.evaluated_cells, 5);
    assert_eq!(limited.search.attempts.len(), 1);
    assert_eq!(limited.search.stop_reason, StopReason::CellLimit);
    assert!(!limited.search.unrestricted_optimum_proven);
    assert_eq!(limited.reference_aligned, "AAAA");
    assert_eq!(limited.query_aligned, "CCCC");
    let full = alignment::align(
        &reference,
        &query,
        &AlignmentOptions {
            max_cells: 100,
            global_storage: GlobalStorage::Full,
            ..options
        },
    )
    .unwrap();
    assert_eq!(full.score, -8);
    assert!(full.search.unrestricted_optimum_proven);
    let reference = Molecule::dna("original reference", "AAAA", Topology::Linear).unwrap();
    let query = Molecule::dna("supplied query", "CCCC", Topology::Linear).unwrap();
    let report = alignment::review::compare(
        &reference,
        &query,
        &alignment::review::Options::Stored(options),
    )
    .unwrap();
    assert!(!report.comparison.optimality_proven);
    report.validate_sources(&reference, &query).unwrap();
    let mut forged = report;
    forged.comparison.optimality_proven = true;
    assert!(forged.validate_sources(&reference, &query).is_err());
}
#[test]
fn whole_comparison_receipt_rejects_candidate_work_and_incumbent_discharge_forgery() {
    let reference = Molecule::dna("reference", "ACGATTCG", Topology::Circular).unwrap();
    let query = Molecule::dna("query", "ATTCGACA", Topology::Circular).unwrap();
    let options = ComparisonOptions::default();
    let report = alignment::compare_molecules(&reference, &query, &options).unwrap();
    report
        .validate_search(8, 8, options.max_total_cells)
        .unwrap();
    let mut changed = report.clone();
    changed.evaluated_cells += 1;
    assert!(
        changed
            .validate_search(8, 8, options.max_total_cells)
            .is_err()
    );
    let mut changed = report.clone();
    changed.candidate_evidence[0].disposition =
        alignment::grid::CandidateDisposition::DischargedByIncumbent;
    assert!(
        changed
            .validate_search(8, 8, options.max_total_cells)
            .is_err()
    );
    let mut changed = report;
    changed.candidate_evidence[0].reference_origin = 99;
    assert!(
        changed
            .validate_search(8, 8, options.max_total_cells)
            .is_err()
    );
}

#[test]
fn raw_worker_alignment_binds_both_sources_options_true_lengths_and_search_receipt() {
    use collagene_core::protocol::{
        execution::process_request,
        request::Request,
        response::Response,
        schema::{AnalysisJob, AnalysisResult},
    };
    let reference = Molecule::dna("separate design", "AAAA", Topology::Linear).unwrap();
    let query = Molecule::dna("supplied consensus", "CCCC", Topology::Linear).unwrap();
    let options = AlignmentOptions {
        mode: AlignmentMode::Global,
        scoring: alignment::Scoring {
            gap_open: -1,
            gap_extend: -1,
            ..Default::default()
        },
        max_cells: 6,
        global_storage: GlobalStorage::Adaptive {
            initial_radius: 0,
            maximum_radius: 4,
            max_attempts: 8,
        },
    };
    let revision = alignment::review::revision(
        &reference,
        &query,
        &alignment::review::Options::Stored(options),
    )
    .unwrap();
    let request = Request {
        id: 391,
        revision: revision.clone(),
        job: AnalysisJob::Align {
            molecule: reference.clone(),
            query: query.clone(),
            options,
        },
    };
    let response = Response::decode(&process_request(&request.encode().unwrap()).unwrap()).unwrap();
    assert_eq!(response.revision, revision);
    let AnalysisResult::Alignment(report) = &response.result else {
        panic!("native alignment response")
    };
    assert_eq!((report.reference_length, report.query_length), (4, 4));
    assert_eq!(report.options, options);
    assert!(!report.search.unrestricted_optimum_proven);
    let mut stale = request.clone();
    let AnalysisJob::Align { options, .. } = &mut stale.job else {
        unreachable!()
    };
    options.max_cells += 1;
    assert!(stale.encode().is_err());
    let mut changed = response.clone();
    let AnalysisResult::Alignment(report) = &mut changed.result else {
        unreachable!()
    };
    report.search.attempts[0].evaluated_cells -= 1;
    assert!(Response::decode(&changed.encode().unwrap()).is_err());
    let mut changed = response;
    let AnalysisResult::Alignment(report) = &mut changed.result else {
        unreachable!()
    };
    report.reference_length += 1;
    assert!(Response::decode(&changed.encode().unwrap()).is_err());
}
