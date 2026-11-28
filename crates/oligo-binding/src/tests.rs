//! Independently enumerated footprints: no calls to the extension recurrence.
use super::*;
use std::collections::BTreeSet;

fn collection(dna: &str, tail: usize) -> Collection {
    Collection {
        oligos: vec![Oligo {
            id: "o1".into(),
            name: "test".into(),
            ordered: Sequence::dna(dna).unwrap(),
            declared_tail_length: tail,
        }],
    }
}
fn options() -> Options {
    Options {
        minimum_annealing_length: 1,
        maximum_extra_unpaired: 200,
        maximum_mismatches: 200,
        minimum_three_prime_pairing: 0,
        match_score: 2,
        mismatch_score: -2,
        minimum_score: -400,
    }
}
fn run(source: &str, topology: Topology, c: &Collection, o: &Options) -> Report {
    search(
        &Sequence::dna(source).unwrap(),
        topology,
        c,
        o,
        &mut Budget::new(BudgetLimits::default()).unwrap(),
    )
    .unwrap()
}
fn complement_literal(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        _ => panic!("concrete oracle only"),
    }
}
type Key = (u8, usize, usize);
fn keys(report: &Report) -> BTreeSet<Key> {
    report
        .hits
        .iter()
        .map(|h| {
            (
                u8::from(h.orientation == Orientation::Reverse),
                h.three_prime_source_position,
                h.length,
            )
        })
        .collect()
}
fn brute(source: &str, topology: Topology, oligo: &Oligo, options: &Options) -> BTreeSet<Key> {
    let source = source.as_bytes();
    let primer = oligo.ordered.as_str().as_bytes();
    let minimum = options.minimum_annealing_length.max(
        (primer.len() - oligo.declared_tail_length).saturating_sub(options.maximum_extra_unpaired),
    );
    let mut result = BTreeSet::new();
    for reverse in [false, true] {
        for endpoint in 0..source.len() {
            let mut scored = Vec::new();
            for length in minimum..=primer.len() - oligo.declared_tail_length {
                if length > source.len() {
                    continue;
                }
                if topology == Topology::Linear
                    && ((!reverse && length > endpoint + 1)
                        || (reverse && endpoint + length > source.len()))
                {
                    continue;
                }
                // Literal alignment in oligo 5′→3′ order, unlike scanner's
                // incremental 3′→5′ recurrence. Recompute every suffix from zero.
                let expected: Vec<u8> = (0..length)
                    .map(|i| {
                        if reverse {
                            complement_literal(source[(endpoint + length - 1 - i) % source.len()])
                        } else {
                            source[(source.len() + endpoint + 1 - length + i) % source.len()]
                        }
                    })
                    .collect();
                let paired: Vec<bool> = primer[primer.len() - length..]
                    .iter()
                    .zip(expected)
                    .map(|(p, s)| *p == s)
                    .collect();
                let mismatches = paired.iter().filter(|p| !**p).count();
                let prime = paired.iter().rev().take_while(|p| **p).count();
                let score = (length - mismatches) as i32 * options.match_score
                    + mismatches as i32 * options.mismatch_score;
                if mismatches <= options.maximum_mismatches
                    && prime >= options.minimum_three_prime_pairing
                    && score >= options.minimum_score
                {
                    scored.push((length, score));
                }
            }
            if let Some(best) = scored.iter().map(|(_, s)| *s).max() {
                for (length, score) in scored {
                    if score == best {
                        result.insert((u8::from(reverse), endpoint, length));
                    }
                }
            }
        }
    }
    result
}
fn words(alphabet: &[u8], maximum: usize) -> Vec<String> {
    let mut result = Vec::new();
    let mut level = vec![String::new()];
    for _ in 0..maximum {
        level = level
            .iter()
            .flat_map(|prefix| {
                alphabet
                    .iter()
                    .map(move |base| format!("{prefix}{}", *base as char))
            })
            .collect();
        result.extend(level.iter().cloned());
    }
    result
}
#[test]
fn exhaustive_concrete_footprints_equal_independent_suffix_oracle() {
    let words = words(b"AT", 4);
    for source in &words {
        for primer in &words {
            for topology in [Topology::Linear, Topology::Circular] {
                let c = collection(primer, 0);
                let o = options();
                let report = run(source, topology, &c, &o);
                assert!(report.complete);
                report
                    .validate_binding(&Sequence::dna(source).unwrap(), topology, &c, &o)
                    .unwrap();
                assert_eq!(
                    keys(&report),
                    brute(source, topology, &c.oligos[0], &o),
                    "{source} {primer} {topology:?}"
                );
                let sequence = Sequence::dna(source).unwrap();
                for hit in &report.hits {
                    hit.validate_columns(&sequence, topology, &c.oligos[0], &o)
                        .unwrap();
                }
            }
        }
    }
}
#[test]
fn extension_recovers_longer_score_after_unfavorable_intermediate_prefix() {
    // 3′-anchored lengths 1..5 score 2,0,2,4,6: stopping at the
    // first decline would miss the actual best footprint.
    let c = collection("AAATA", 0);
    let report = run("AAAAA", Topology::Linear, &c, &options());
    let hits: Vec<_> = report
        .hits
        .iter()
        .filter(|h| h.orientation == Orientation::Forward && h.three_prime_source_position == 4)
        .collect();
    assert_eq!(hits.len(), 1);
    assert_eq!(hits[0].length, 5);
    assert_eq!(hits[0].maximum_mismatches, 1);
}
#[test]
fn terminal_mismatch_is_reported_without_inventing_extension_success() {
    let c = collection("AC", 0);
    let report = run("AA", Topology::Linear, &c, &options());
    assert!(
        report
            .hits
            .iter()
            .any(|h| h.orientation == Orientation::Forward
                && h.three_prime_certain == 0
                && h.minimum_mismatches > 0)
    );
    let mut o = options();
    o.minimum_three_prime_pairing = 1;
    assert!(
        run("AA", Topology::Linear, &c, &o)
            .hits
            .iter()
            .all(|h| h.three_prime_possible >= 1)
    );
}
#[test]
fn all_best_score_ties_retained_or_endpoint_rejected_atomically() {
    // At forward endpoint2, lengths1 and3 both score2.
    let c = collection("ATA", 0);
    let report = run("AAA", Topology::Linear, &c, &options());
    let lengths: Vec<_> = report
        .hits
        .iter()
        .filter(|h| h.orientation == Orientation::Forward && h.three_prime_source_position == 2)
        .map(|h| h.length)
        .collect();
    assert_eq!(lengths, vec![1, 3]);
    let mut budget = Budget::new(BudgetLimits {
        hits: 3,
        ..BudgetLimits::default()
    })
    .unwrap();
    let limited = search(
        &Sequence::dna("AAA").unwrap(),
        Topology::Linear,
        &c,
        &options(),
        &mut budget,
    )
    .unwrap();
    assert_eq!(
        limited.oligos[0].status,
        Status::Limited(Limit::HitRetention)
    );
    assert!(!limited.complete);
    assert!(
        !limited
            .hits
            .iter()
            .any(|h| h.orientation == Orientation::Forward && h.three_prime_source_position == 2)
    );
    assert_eq!(limited.oligos[0].alignments_found, 4);
}
#[test]
fn declared_tail_extra_prefix_and_circular_reverse_columns_remain_distinct() {
    let c = collection("GGTAC", 2);
    let mut o = options();
    o.maximum_mismatches = 0;
    let report = run("AGT", Topology::Circular, &c, &o);
    assert!(
        report
            .hits
            .iter()
            .any(|h| h.orientation == Orientation::Reverse && h.crosses_origin)
    );
    for h in &report.hits {
        assert_eq!(h.declared_tail_length, 2);
        assert_eq!(h.annealing_start, 2 + h.extra_unpaired_length);
        h.validate_columns(
            &Sequence::dna("AGT").unwrap(),
            Topology::Circular,
            &c.oligos[0],
            &o,
        )
        .unwrap();
    }
}
#[test]
fn ambiguity_envelope_contains_every_concrete_realization_optimum() {
    for source in words(b"AN", 3) {
        for primer in words(b"TN", 3) {
            let c = collection(&primer, 0);
            let mut o = options();
            o.maximum_mismatches = 1;
            let envelope = run(&source, Topology::Circular, &c, &o);
            let included = keys(&envelope);
            let expand = |text: &str| {
                let mut values = vec![String::new()];
                for b in text.bytes() {
                    let bases: &[u8] = if b == b'N' {
                        b"ACGT"
                    } else {
                        std::slice::from_ref(&b)
                    };
                    values = values
                        .iter()
                        .flat_map(|prefix| {
                            bases
                                .iter()
                                .map(move |base| format!("{prefix}{}", *base as char))
                        })
                        .collect();
                }
                values
            };
            for concrete_source in expand(&source) {
                for concrete_primer in expand(&primer) {
                    for key in brute(
                        &concrete_source,
                        Topology::Circular,
                        &collection(&concrete_primer, 0).oligos[0],
                        &o,
                    ) {
                        assert!(
                            included.contains(&key),
                            "missing {source}/{primer}: {key:?}"
                        );
                    }
                }
            }
            for hit in &envelope.hits {
                hit.validate_columns(
                    &Sequence::dna(&source).unwrap(),
                    Topology::Circular,
                    &c.oligos[0],
                    &o,
                )
                .unwrap();
            }
        }
    }
}
#[test]
fn exhausted_collection_does_not_claim_unsearched_oligos_are_nonbinding() {
    let mut c = collection("AAA", 0);
    let mut other = c.oligos[0].clone();
    other.id = "o2".into();
    c.oligos.push(other);
    let mut budget = Budget::new(BudgetLimits {
        hits: 0,
        ..BudgetLimits::default()
    })
    .unwrap();
    let r = search(
        &Sequence::dna("AAA").unwrap(),
        Topology::Linear,
        &c,
        &options(),
        &mut budget,
    )
    .unwrap();
    assert_eq!(r.oligos[0].status, Status::Limited(Limit::HitRetention));
    assert_eq!(r.oligos[1].status, Status::Unsearched);
    assert_eq!(r.oligos[1].endpoints_examined, 0);
    assert!(!r.complete);
    let used = budget.used();
    let second = search(
        &Sequence::dna("AAA").unwrap(),
        Topology::Linear,
        &c,
        &options(),
        &mut budget,
    )
    .unwrap();
    assert_eq!(second.budget_before, used);
    assert!(second.budget_after.work > used.work);
}
#[test]
fn impossible_geometry_is_complete_zero_sites_but_work_shortage_is_explicit() {
    let c = collection("AAA", 0);
    let mut o = options();
    o.minimum_annealing_length = 3;
    let r = run("AA", Topology::Circular, &c, &o);
    assert!(r.complete);
    assert!(r.hits.is_empty());
    assert_eq!(r.oligos[0].endpoints_total, 0);
    let input = 3 + c.oligos[0].ordered.len() + c.oligos[0].id.len() + c.oligos[0].name.len() + 256;
    let source = Sequence::dna("AAA").unwrap();
    let admission = Budget::admission_cost(&source, &c).unwrap();
    assert_eq!(admission.work, input as u64);
    for limits in [
        BudgetLimits {
            work: admission.work - 1,
            ..BudgetLimits::default()
        },
        BudgetLimits {
            output_charge: admission.output_charge - 1,
            ..BudgetLimits::default()
        },
    ] {
        let mut denied = Budget::new(limits).unwrap();
        assert!(matches!(
            search(&source, Topology::Linear, &c, &o, &mut denied),
            Err(Error::TooLarge)
        ));
        assert_eq!(denied.used(), Usage::default());
    }
    let mut budget = Budget::new(BudgetLimits {
        work: input as u64,
        ..BudgetLimits::default()
    })
    .unwrap();
    let r = search(
        &Sequence::dna("AAA").unwrap(),
        Topology::Linear,
        &c,
        &o,
        &mut budget,
    )
    .unwrap();
    assert_eq!(r.oligos[0].status, Status::Limited(Limit::Work));
    assert!(!r.complete);
    assert_eq!(r.oligos[0].endpoints_examined, 0);
}
#[test]
fn identities_and_aligned_facts_detect_context_edits() {
    let source = Sequence::dna("ACGT").unwrap();
    let c = collection("ACGT", 0);
    let o = options();
    assert_ne!(
        source_content_hash(&source, Topology::Linear).unwrap(),
        source_content_hash(&source, Topology::Circular).unwrap()
    );
    let mut changed = c.clone();
    changed.oligos[0].declared_tail_length = 1;
    assert_ne!(c.content_hash().unwrap(), changed.content_hash().unwrap());
    changed = c.clone();
    changed.oligos[0].name = "another".into();
    assert_ne!(c.content_hash().unwrap(), changed.content_hash().unwrap());
    let mut settings = o.clone();
    settings.maximum_mismatches = 0;
    assert_ne!(o.content_hash().unwrap(), settings.content_hash().unwrap());
    let mut h = run("ACGT", Topology::Circular, &c, &o).hits.remove(0);
    h.columns[0].source_position = 999;
    assert!(
        h.validate_columns(&source, Topology::Circular, &c.oligos[0], &o)
            .is_err()
    );
    let mut invalid = c;
    invalid.oligos[0].declared_tail_length = 4;
    assert!(invalid.validate().is_err());
}

#[test]
fn report_binding_rejects_wrong_model_context_and_completion_claim() {
    let source = Sequence::dna("ACGT").unwrap();
    let c = collection("ACGT", 0);
    let o = options();
    let original = run("ACGT", Topology::Circular, &c, &o);
    original
        .validate_binding(&source, Topology::Circular, &c, &o)
        .unwrap();
    let mut tampered = original.clone();
    tampered.model = "other-model".into();
    assert!(
        tampered
            .validate_binding(&source, Topology::Circular, &c, &o)
            .is_err()
    );
    let mut tampered = original.clone();
    tampered.oligos[0].endpoints_examined = 0;
    assert!(
        tampered
            .validate_binding(&source, Topology::Circular, &c, &o)
            .is_err()
    );
    let mut tampered = original.clone();
    tampered.complete = false;
    assert!(
        tampered
            .validate_binding(&source, Topology::Circular, &c, &o)
            .is_err()
    );
    let mut tampered = original;
    tampered.hits[0].oligo_index = usize::MAX;
    assert!(
        tampered
            .validate_binding(&source, Topology::Circular, &c, &o)
            .is_err()
    );
    assert!(matches!(
        collection(&"A".repeat(MAX_OLIGO_BASES + 1), 0).validate(),
        Err(Error::TooLarge)
    ));
}
