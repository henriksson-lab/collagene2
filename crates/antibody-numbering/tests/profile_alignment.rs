use collagene_antibody_numbering::{Error, profile::*};
use collagene_model::{Alphabet, Molecule, PROTEIN_RESIDUES, Sequence, Topology};
fn source(s: &str) -> Molecule {
    Molecule {
        name: "Authored profile input".into(),
        sequence: Sequence::new(Alphabet::Protein, s).unwrap(),
        topology: Topology::Linear,
        features: vec![],
        metadata: Default::default(),
        original_genbank: None,
    }
}
fn profile(s: &str) -> Profile {
    let bytes = s.as_bytes();
    Profile {
        version: PROFILE_VERSION,
        name: "Authored scoring profile, not a trained antibody model".into(),
        provenance: "Synthetic score/path fixture".into(),
        states: bytes
            .iter()
            .enumerate()
            .map(|(i, &b)| {
                let mut emissions = [Some(-4); 22];
                emissions[PROTEIN_RESIDUES.iter().position(|&a| a == b).unwrap()] = Some(5);
                MatchState {
                    emissions,
                    required: true,
                    allow_insert_after: i + 1 < bytes.len(),
                }
            })
            .collect(),
        gap_open: -3,
        gap_extend: -1,
        minimum_score: 5,
    }
}
fn limits() -> Limits {
    Limits {
        max_cells: MAX_CELLS,
    }
}
#[test]
fn embedded_affine_insertions_and_optional_deletions_preserve_original_source_paths() {
    let m = source("QQAQQCGTT");
    let p = profile("ACG");
    let r = align(&m, &p, &limits(), &mut || false).unwrap();
    let a = r.lower.as_ref().unwrap();
    assert_eq!(a.score, 11);
    assert_eq!((a.start, a.end), (2, 7));
    assert!(!a.alternative_optimal_paths);
    assert_eq!(r.lower, r.upper);
    assert_eq!(
        a.steps,
        vec![
            Step::Match {
                source_index: 2,
                model_index: 0
            },
            Step::Insert {
                source_index: 3,
                after_model: 0
            },
            Step::Insert {
                source_index: 4,
                after_model: 0
            },
            Step::Match {
                source_index: 5,
                model_index: 1
            },
            Step::Match {
                source_index: 6,
                model_index: 2
            }
        ]
    );
    r.validate_for(&m, &p, &limits()).unwrap();
    let mut p = profile("ACG");
    p.states[1].required = false;
    let r = align(&source("AG"), &p, &limits(), &mut || false).unwrap();
    assert_eq!(r.lower.as_ref().unwrap().score, 7);
    assert!(
        r.lower
            .unwrap()
            .steps
            .contains(&Step::Delete { model_index: 1 })
    );
    p.states[1].required = true;
    assert_eq!(
        align(&source("AG"), &p, &limits(), &mut || false)
            .unwrap()
            .threshold_status,
        ThresholdStatus::NoPath
    );
}
#[test]
fn ties_stops_and_thresholds_do_not_turn_score_hits_into_unique_domains() {
    let p = profile("A");
    let r = align(&source("AAA"), &p, &limits(), &mut || false).unwrap();
    let a = r.lower.unwrap();
    assert_eq!((a.start, a.end), (0, 1));
    assert!(a.alternative_optimal_paths);
    assert_eq!(r.threshold_status, ThresholdStatus::GuaranteedScore);
    assert_eq!(
        align(&source("C"), &p, &limits(), &mut || false)
            .unwrap()
            .threshold_status,
        ThresholdStatus::BelowThreshold
    );
    let p = profile("ACG");
    assert_eq!(
        align(&source("AC*G"), &p, &limits(), &mut || false)
            .unwrap()
            .threshold_status,
        ThresholdStatus::NoPath
    );
    let r = align(&source("M*ACG"), &p, &limits(), &mut || false).unwrap();
    assert_eq!(
        (
            r.lower.as_ref().unwrap().start,
            r.lower.as_ref().unwrap().end
        ),
        (2, 5)
    );
}
// Independent path-language enumeration: recurse over complete source/model
// consumption paths, without any production matrix, traceback or emission cache.
fn brute(s: &str, p: &Profile) -> Option<(i32, usize)> {
    brute_end(s, p, None).map(|r| (r.score, r.ways))
}
#[derive(Debug)]
struct Oracle {
    score: i32,
    ways: usize,
    matches: std::collections::BTreeSet<(usize, usize)>,
}
fn brute_end(s: &str, p: &Profile, wanted: Option<usize>) -> Option<Oracle> {
    fn visit(
        s: &[u8],
        p: &Profile,
        i: usize,
        j: usize,
        last: u8,
        score: i32,
        best: &mut Option<Oracle>,
        wanted: Option<usize>,
        path: &mut Vec<(usize, usize)>,
    ) {
        if j == p.states.len() && wanted.is_none_or(|end| end == i) {
            match best {
                Some(old) if old.score == score => {
                    old.ways = (old.ways + 1).min(2);
                    old.matches.retain(|pair| path.contains(pair));
                }
                Some(old) if old.score > score => (),
                _ => {
                    *best = Some(Oracle {
                        score,
                        ways: 1,
                        matches: path.iter().copied().collect(),
                    })
                }
            }
        }
        if i < s.len() && j < p.states.len() && s[i] != b'*' {
            let column = b"ACDEFGHIKLMNPQRSTVWYOU"
                .iter()
                .position(|&b| b == s[i])
                .unwrap();
            if let Some(e) = p.states[j].emissions[column] {
                path.push((j, i));
                visit(
                    s,
                    p,
                    i + 1,
                    j + 1,
                    0,
                    score + i32::from(e),
                    best,
                    wanted,
                    path,
                );
                path.pop();
            }
        }
        if j < p.states.len() && !p.states[j].required && last != 2 {
            let cost = if last == 1 { p.gap_extend } else { p.gap_open };
            visit(
                s,
                p,
                i,
                j + 1,
                1,
                score + i32::from(cost),
                best,
                wanted,
                path,
            );
        }
        if i < s.len() && j > 0 && p.states[j - 1].allow_insert_after && last != 1 && s[i] != b'*' {
            let cost = if last == 2 { p.gap_extend } else { p.gap_open };
            visit(
                s,
                p,
                i + 1,
                j,
                2,
                score + i32::from(cost),
                best,
                wanted,
                path,
            );
        }
    }
    let mut best = None;
    for i in 0..=s.len() {
        visit(s.as_bytes(), p, i, 0, 0, 0, &mut best, wanted, &mut vec![]);
    }
    best
}
#[test]
fn every_small_canonical_input_agrees_with_independent_complete_path_enumeration() {
    for template in ["A", "AC", "ACA"] {
        let mut p = profile(template);
        if p.states.len() == 3 {
            p.states[1].required = false;
        }
        for len in 1..=4 {
            for mut word in 0..4_usize.pow(len) {
                let mut s = String::new();
                for _ in 0..len {
                    s.push(b"ACG*"[word % 4] as char);
                    word /= 4;
                }
                let expected = brute(&s, &p);
                let r = align(&source(&s), &p, &limits(), &mut || false).unwrap();
                assert_eq!(
                    r.lower
                        .as_ref()
                        .map(|a| (a.score, if a.alternative_optimal_paths { 2 } else { 1 })),
                    expected,
                    "{template}/{s}"
                );
                assert_eq!(r.lower, r.upper);
            }
        }
    }
}
fn expand(s: &str) -> Vec<String> {
    let mut out = vec![String::new()];
    for b in s.chars() {
        let choices = match b {
            'B' => "DN".into(),
            'Z' => "EQ".into(),
            'J' => "IL".into(),
            'X' => "ACDEFGHIKLMNPQRSTVWYOU".into(),
            b => b.to_string(),
        };
        out = out
            .into_iter()
            .flat_map(|s| {
                choices
                    .chars()
                    .map(move |c| format!("{s}{c}"))
                    .collect::<Vec<_>>()
            })
            .collect();
    }
    out
}
#[test]
fn ambiguity_bounds_cover_every_concrete_best_score_and_impossible_emissions_remain_explicit() {
    for text in ["AX", "BX", "XCA", "JAC", "ZC", "XX"] {
        let p = profile("AC");
        let r = align(&source(text), &p, &limits(), &mut || false).unwrap();
        let concrete = expand(text)
            .into_iter()
            .map(|s| brute(&s, &p).map(|(score, _)| score))
            .collect::<Vec<_>>();
        let upper = concrete.iter().flatten().copied().max();
        assert_eq!(r.upper.as_ref().map(|a| a.score), upper, "{text}");
        if concrete.iter().any(Option::is_none) {
            assert!(r.lower.is_none());
        } else if let Some(lower) = &r.lower {
            assert!(concrete.iter().flatten().all(|&score| lower.score <= score));
        }
    }
    let mut p = profile("C");
    p.states[0].emissions = [None; 22];
    p.states[0].emissions[1] = Some(5);
    let r = align(&source("X"), &p, &limits(), &mut || false).unwrap();
    assert!(r.lower.is_none());
    assert_eq!(r.upper.unwrap().score, 5);
    assert_eq!(r.threshold_status, ThresholdStatus::PossibleScore);
    assert_eq!(
        align(&source("A"), &p, &limits(), &mut || false)
            .unwrap()
            .threshold_status,
        ThresholdStatus::NoPath
    );
}
#[test]
fn current_profiles_budgets_source_and_reconstructed_evidence_are_checked() {
    let m = source("AC");
    let p = profile("AC");
    let r = align(&m, &p, &limits(), &mut || false).unwrap();
    let round: Report = serde_json::from_slice(&serde_json::to_vec(&r).unwrap()).unwrap();
    round.validate_for(&m, &p, &limits()).unwrap();
    let mut changed = r.clone();
    changed.lower.as_mut().unwrap().score += 1;
    assert!(changed.validate_for(&m, &p, &limits()).is_err());
    let mut changed = p.clone();
    changed.provenance = "Different model origin".into();
    assert!(r.validate_for(&m, &changed, &limits()).is_err());
    let mut bad = p.clone();
    bad.version = 0;
    assert!(bad.revision().is_err());
    bad.version = 2;
    assert!(bad.revision().is_err());
    bad = p.clone();
    bad.states.iter_mut().for_each(|s| s.required = false);
    assert!(bad.revision().is_err());
    bad = p.clone();
    bad.gap_open = 1;
    assert!(bad.revision().is_err());
    bad = p.clone();
    bad.states[0].emissions[0] = Some(1001);
    assert!(bad.revision().is_err());
    assert!(align(&m, &p, &Limits { max_cells: 8 }, &mut || false).is_err());
    assert_eq!(
        align(&m, &p, &limits(), &mut || true).unwrap_err(),
        Error::Cancelled
    );
    let dna = Molecule::dna("DNA", "AC", Topology::Linear).unwrap();
    assert!(align(&dna, &p, &limits(), &mut || false).is_err());
}

#[test]
fn every_small_profile_endpoint_agrees_with_complete_path_enumeration() {
    for (template, free_gaps) in [
        ("A", false),
        ("AC", false),
        ("ACA", false),
        ("A", true),
        ("ACA", true),
    ] {
        let mut p = profile(template);
        if free_gaps {
            p.gap_open = 0;
            p.gap_extend = 0;
        }
        p.minimum_score = i32::MIN;
        if template == "A" {
            p.states[0].allow_insert_after = true;
        }
        if template == "ACA" {
            p.states[1].required = false;
        }
        for len in 1..=4 {
            for mut word in 0..4_usize.pow(len) {
                let mut s = String::new();
                for _ in 0..len {
                    s.push(b"ACG*"[word % 4] as char);
                    word /= 4;
                }
                let r = align_endpoints(
                    &source(&s),
                    &p,
                    &EndpointLimits {
                        max_cells: MAX_CELLS,
                        max_endpoints: MAX_ENDPOINTS,
                        max_trace_states: MAX_TRACE_STATES,
                    },
                    &mut || false,
                )
                .unwrap();
                let expected = (0..=s.len())
                    .filter_map(|end| brute_end(&s, &p, Some(end)).map(|r| (end, r.score, r.ways)))
                    .collect::<Vec<_>>();
                assert_eq!(r.reachable_endpoints, expected.len());
                assert_eq!(r.qualified_endpoints, expected.len());
                assert!(!r.truncated);
                assert_eq!(
                    r.endpoints
                        .iter()
                        .map(|e| (
                            e.end,
                            e.upper.score,
                            if e.upper.alternative_optimal_paths {
                                2
                            } else {
                                1
                            }
                        ))
                        .collect::<Vec<_>>(),
                    expected,
                    "{template}/{s}"
                );
                assert!(
                    r.endpoints
                        .iter()
                        .all(|e| e.lower.as_ref() == Some(&e.upper))
                );
                for endpoint in &r.endpoints {
                    let oracle = brute_end(&s, &p, Some(endpoint.end)).unwrap();
                    assert_eq!(
                        endpoint
                            .invariant_matches
                            .as_ref()
                            .unwrap()
                            .iter()
                            .map(|m| (m.model_index, m.source_index))
                            .collect::<std::collections::BTreeSet<_>>(),
                        oracle.matches,
                        "{template}/{s}/{}",
                        endpoint.end
                    );
                }
            }
        }
    }
}
#[test]
fn endpoint_retention_counts_ambiguity_and_evidence_are_source_bound() {
    let p = profile("A");
    let m = source("AAA");
    let limits = EndpointLimits {
        max_cells: MAX_CELLS,
        max_endpoints: 1,
        max_trace_states: MAX_TRACE_STATES,
    };
    let r = align_endpoints(&m, &p, &limits, &mut || false).unwrap();
    assert_eq!(
        (
            r.reachable_endpoints,
            r.qualified_endpoints,
            r.endpoints.len()
        ),
        (3, 3, 1)
    );
    assert!(r.truncated);
    assert_eq!((r.endpoints[0].upper.start, r.endpoints[0].end), (0, 1));
    r.validate_for(&m, &p, &limits).unwrap();
    let mut forged = r.clone();
    forged.endpoints[0].end = 2;
    assert!(forged.validate_for(&m, &p, &limits).is_err());
    let zero = align_endpoints(
        &m,
        &p,
        &EndpointLimits {
            max_cells: MAX_CELLS,
            max_endpoints: 0,
            max_trace_states: MAX_TRACE_STATES,
        },
        &mut || false,
    )
    .unwrap();
    assert_eq!(zero.qualified_endpoints, 3);
    assert!(zero.truncated && zero.endpoints.is_empty());
    let mut high = p.clone();
    high.minimum_score = 6;
    let filtered = align_endpoints(&m, &high, &limits, &mut || false).unwrap();
    assert_eq!(
        (filtered.reachable_endpoints, filtered.qualified_endpoints),
        (3, 0)
    );
    assert!(!filtered.truncated);
    assert!(
        align_endpoints(
            &m,
            &p,
            &EndpointLimits {
                max_cells: MAX_CELLS,
                max_endpoints: MAX_ENDPOINTS + 1,
                max_trace_states: MAX_TRACE_STATES,
            },
            &mut || false
        )
        .is_err()
    );
    let mut calls = 0;
    assert_eq!(
        align_endpoints(&m, &p, &limits, &mut || {
            calls += 1;
            calls == 3
        })
        .unwrap_err(),
        Error::Cancelled
    );
    for s in ["AX", "BX", "XCA", "JAC", "ZC", "XX"] {
        let mut p = profile("AC");
        p.minimum_score = i32::MIN;
        let r = align_endpoints(
            &source(s),
            &p,
            &EndpointLimits {
                max_cells: MAX_CELLS,
                max_endpoints: MAX_ENDPOINTS,
                max_trace_states: MAX_TRACE_STATES,
            },
            &mut || false,
        )
        .unwrap();
        for e in r.endpoints {
            assert!(e.invariant_matches.is_none());
            let values = expand(s)
                .iter()
                .map(|v| brute_end(v, &p, Some(e.end)).map(|r| r.score))
                .collect::<Vec<_>>();
            assert_eq!(Some(e.upper.score), values.iter().flatten().copied().max());
            if let Some(lower) = e.lower {
                assert!(values.iter().all(|v| v.is_some_and(|n| n >= lower.score)));
            }
        }
    }
}

#[test]
fn invariant_graph_work_is_bounded_and_counts_only_does_not_require_a_trace() {
    let m = source("AAA");
    let p = profile("A");
    let limits = EndpointLimits {
        max_cells: MAX_CELLS,
        max_endpoints: 1,
        max_trace_states: 1,
    };
    assert_eq!(
        align_endpoints(&m, &p, &limits, &mut || false).unwrap_err(),
        Error::Limit("profile invariant trace states")
    );
    let counted = align_endpoints(
        &m,
        &p,
        &EndpointLimits {
            max_cells: MAX_CELLS,
            max_endpoints: 0,
            max_trace_states: 0,
        },
        &mut || false,
    )
    .unwrap();
    assert_eq!(counted.qualified_endpoints, 3);
    assert_eq!(counted.invariant_trace_states, 0);
    let mut limits = limits;
    limits.max_trace_states = MAX_TRACE_STATES + 1;
    assert!(align_endpoints(&m, &p, &limits, &mut || false).is_err());
}
