use collagene_core::{
    alignment::{AlignmentMode, AlignmentOptions, ComparisonOptions, DifferenceKind},
    pairwise_review::{self, Options},
    protocol::{
        self,
        schema::{AnalysisJob, AnalysisResult},
    },
};
use collagene_model::{Molecule, Topology};
fn molecule(name: &str, dna: &str, topology: Topology) -> Molecule {
    Molecule::dna(name, dna, topology).unwrap()
}
#[test]
fn stored_differences_and_empty_spans_remain_navigable_in_native_coordinates() {
    let reference = molecule("reference", "ACGTT", Topology::Linear);
    let query = molecule("insertion", "ACGATT", Topology::Linear);
    let report = pairwise_review::compare(
        &reference,
        &query,
        &Options::Stored(AlignmentOptions::default()),
    )
    .unwrap();
    let insertion = report
        .comparison
        .alignment
        .differences
        .iter()
        .find(|d| d.kind == DifferenceKind::Insertion)
        .unwrap();
    assert_eq!((insertion.reference_start, insertion.reference_end), (3, 3));
    assert_eq!((insertion.query_start, insertion.query_end), (3, 4));
    assert_eq!(insertion.query_bases, "A");
    assert_eq!(report.map_reference(3, 3).unwrap(), vec![(3, 3)]);
    assert_eq!(report.map_query(3, 4).unwrap(), (3, 4));
    let deletion = pairwise_review::compare(
        &query,
        &reference,
        &Options::Stored(AlignmentOptions::default()),
    )
    .unwrap();
    let d = deletion
        .comparison
        .alignment
        .differences
        .iter()
        .find(|d| d.kind == DifferenceKind::Deletion)
        .unwrap();
    assert_eq!(
        (
            d.reference_start,
            d.reference_end,
            d.query_start,
            d.query_end
        ),
        (3, 4, 3, 3)
    );
    assert_eq!(deletion.map_query(3, 3).unwrap(), (3, 3));
    let changed = molecule("substitution", "ACATT", Topology::Linear);
    let difference = pairwise_review::compare(
        &reference,
        &changed,
        &Options::Stored(AlignmentOptions::default()),
    )
    .unwrap();
    assert_eq!(difference.comparison.alignment.score, 5);
    assert_eq!(
        difference.comparison.alignment.differences[0].reference_bases,
        "G"
    );
    assert_eq!(
        difference.comparison.alignment.differences[0].query_bases,
        "A"
    );
    let mut tampered = report.clone();
    tampered.comparison.alignment.differences[0].query_start = 0;
    assert!(tampered.validate_binding(&report.context_revision).is_err());
    let mut other = reference.clone();
    other.name = "other source event content".into();
    assert!(report.validate_sources(&other, &query).is_err());
    let mut tampered = report.clone();
    tampered.comparison.alignment.score += 1;
    assert!(tampered.validate_sources(&reference, &query).is_err());
}
#[test]
fn circular_origin_reverse_query_and_capped_search_keep_explicit_policy() {
    let reference = molecule("circle", "AGTCGATC", Topology::Circular);
    // Reverse complement of reference rotated left by three, specified literally.
    let query = molecule("reversed rotated", "ACTGATCG", Topology::Circular);
    let report = pairwise_review::compare(
        &reference,
        &query,
        &Options::Whole(ComparisonOptions::default()),
    )
    .unwrap();
    assert_eq!(report.comparison.reference_origin, 3);
    assert!(report.comparison.query_reversed);
    assert_eq!(report.comparison.alignment.reference_aligned, "CGATCAGT");
    assert_eq!(report.comparison.alignment.query_aligned, "CGATCAGT");
    assert_eq!(report.map_reference(4, 6).unwrap(), vec![(7, 8), (0, 1)]);
    assert_eq!(report.map_reference(5, 5).unwrap(), vec![(0, 0)]);
    assert_eq!(report.map_query(0, 2).unwrap(), (6, 8));
    assert_eq!(report.map_query(0, 0).unwrap(), (8, 8));
    assert_eq!(reference.sequence.as_str(), "AGTCGATC");
    assert_eq!(query.sequence.as_str(), "ACTGATCG");
    assert!(report.comparison.optimality_proven);
    assert!(!report.comparison.exhaustive_reference_origins);
    let r = molecule("bounded circle", "AGTCGATCGGACTTAGCCTA", Topology::Circular);
    let q = molecule(
        "one changed base",
        "AGTCGATCGGACTTAGCCTC",
        Topology::Circular,
    );
    let options = ComparisonOptions {
        max_origins: 1,
        max_total_cells: 441,
        ..Default::default()
    };
    let capped = pairwise_review::compare(&r, &q, &Options::Whole(options)).unwrap();
    assert_eq!(capped.comparison.evaluated_candidates, 1);
    assert!(!capped.comparison.optimality_proven);
    assert!(!capped.comparison.exhaustive_reference_origins);
    assert!(!capped.comparison.exhaustive_query_orientations);
    assert!(!capped.comparison.warnings.is_empty());
    capped.validate_sources(&r, &q).unwrap();
}
#[test]
fn local_flanks_and_zero_local_alignment_are_not_full_identity() {
    let reference = molecule("flanked reference", "GGGACGTTCCC", Topology::Linear);
    let query = molecule("local query", "ACGTT", Topology::Linear);
    let options = Options::Stored(AlignmentOptions {
        mode: AlignmentMode::Local,
        ..Default::default()
    });
    let report = pairwise_review::compare(&reference, &query, &options).unwrap();
    let a = &report.comparison.alignment;
    assert_eq!(
        (
            a.reference_start,
            a.reference_end,
            a.query_start,
            a.query_end
        ),
        (3, 8, 0, 5)
    );
    assert_eq!(a.stats.reference_coverage, 5. / 11.);
    assert_eq!(a.stats.certain_identity, Some(1.));
    assert_eq!(report.map_reference(3, 8).unwrap(), vec![(3, 8)]);
    let mut forged = report.clone();
    forged.options = Options::Stored(AlignmentOptions::default());
    forged.comparison.alignment.options = AlignmentOptions::default();
    assert!(forged.validate_binding(&report.context_revision).is_err());
    let r = molecule("no local pairs", "AAAA", Topology::Linear);
    let q = molecule("disjoint", "CCCC", Topology::Linear);
    let empty = pairwise_review::compare(&r, &q, &options).unwrap();
    assert_eq!(empty.comparison.alignment.stats.columns, 0);
    assert_eq!(empty.comparison.alignment.stats.certain_identity, None);
    assert_eq!(empty.comparison.alignment.stats.reference_coverage, 0.);
    empty.validate_sources(&r, &q).unwrap();
    let whole = Options::Whole(ComparisonOptions {
        alignment: options.alignment(),
        ..Default::default()
    });
    assert!(pairwise_review::revision(&r, &q, &whole).is_err());
}
#[test]
fn new_binary_selectors_bind_sources_and_reject_old_prefixes() {
    let reference = molecule("reference", "ACGTT", Topology::Linear);
    let query = molecule("query", "ACGATT", Topology::Linear);
    let options = Options::Stored(AlignmentOptions::default());
    let revision = pairwise_review::revision(&reference, &query, &options).unwrap();
    let request = protocol::request::Request {
        id: 52,
        revision: revision.clone(),
        job: AnalysisJob::CompareSequences {
            reference,
            query,
            options,
        },
    };
    let wire = protocol::request::Request::encode(&request).unwrap();
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&wire[42..]).unwrap();
    assert_eq!(tag, 117);
    let mut old = wire.clone();
    old[4..6].copy_from_slice(&51u16.to_le_bytes());
    assert!(matches!(
        protocol::execution::process_request(&old),
        Err(collagene_model::Error::Unsupported(_))
    ));
    let body = postcard::to_allocvec(&(request.id, request.revision.as_str(), 117u32)).unwrap();
    let mut prefix = wire[..42].to_vec();
    prefix[4..6].copy_from_slice(&51u16.to_le_bytes());
    prefix[6..10].copy_from_slice(&(body.len() as u32).to_le_bytes());
    prefix.extend_from_slice(&body);
    assert!(matches!(
        protocol::execution::process_request(&prefix),
        Err(collagene_model::Error::Unsupported(_))
    ));
    let response = protocol::execution::process_request(&wire).unwrap();
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&response[42..]).unwrap();
    assert_eq!(tag, 109);
    let decoded = protocol::response::Response::decode(&response).unwrap();
    let AnalysisResult::ReviewedComparison(report) = decoded.result else {
        panic!("wrong result")
    };
    report.validate_binding(&revision).unwrap();
    let mut old = response;
    old[4..6].copy_from_slice(&51u16.to_le_bytes());
    assert!(matches!(
        protocol::response::Response::decode(&old),
        Err(collagene_model::Error::Unsupported(_))
    ));
    let mut stale = request;
    stale.revision = "0".repeat(64);
    assert!(
        protocol::execution::process_request(&protocol::request::Request::encode(&stale).unwrap())
            .is_err()
    );
}
