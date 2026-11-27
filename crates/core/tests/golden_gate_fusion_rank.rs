//! Source/profile binding, explicit materialization and compact score oracles.
use collagene_core::{
    Molecule, Sequence, Topology, assembly_intent,
    golden_gate::GoldenGateEnzyme,
    golden_gate_fusion_design::{self as fusion, Constraints, Limits, Window},
    golden_gate_fusion_rank::{self as rank, PairEstimate, Score, ScoreClass},
    golden_gate_junction_edit::{self as gate, Row},
    ligation_fidelity::UnknownReason,
    ligation_profiles as profiles, ligation_recipe,
    overlap_primer_design::Fragment,
    verification::{PlanInput, Status},
};
const A: &str = "ATGACCCCGGGGTTTCCAGATCGC";
const B: &str = "ATGAGTTAGCGACCTTGGAAACGT";
fn profile(formulation: &str) -> profiles::Summary {
    profiles::summaries()
        .unwrap()
        .into_iter()
        .find(|p| p.enzyme_formulation == formulation)
        .expect("published formulation")
}
fn fixture(b: &str, circular: bool) -> (Vec<PlanInput>, rank::Request) {
    let profile = profile("BsaI-HFv2");
    let inputs = [("a", A), ("b", b)]
        .into_iter()
        .map(|(id, dna)| {
            let molecule = Molecule::dna(id, dna, Topology::Linear).unwrap();
            PlanInput {
                id: id.into(),
                revision: molecule.revision(),
                molecule,
            }
        })
        .collect();
    let fragments = vec![
        Fragment {
            role: "left".into(),
            input_id: "a".into(),
            start: 0,
            end: A.len(),
            reverse: false,
        },
        Fragment {
            role: "right".into(),
            input_id: "b".into(),
            start: 0,
            end: b.len(),
            reverse: false,
        },
    ];
    let rows = fragments
        .iter()
        .map(|f| Row::DesignPcr {
            role: f.role.clone(),
            input_id: f.input_id.clone(),
            reverse: f.reverse,
            forward_length: 8,
            reverse_length: 8,
        })
        .collect();
    let design = fusion::Request {
        base: gate::Request {
            base: assembly_intent::Request {
                name: "Independent source AB".into(),
                fragments,
                circular,
            },
            edits: vec![],
            rows,
            enzyme: GoldenGateEnzyme::BsaI,
            clamp: "ATGC".into(),
            spacer: "A".into(),
            terminal_right_overhang: if circular { None } else { Some("GGAT".into()) },
            allow_nicked: false,
        },
        windows: vec![Window {
            junction: 0,
            min_shift: -1,
            max_shift: 1,
        }],
        constraints: Constraints {
            reject_palindromes: false,
            distinct_reverse_complements: false,
            minimum_hamming: 0,
        },
        limits: Limits::default(),
    };
    (
        inputs,
        rank::Request {
            design,
            profile: profile.profile,
            enzyme_formulation: "BsaI-HFv2".into(),
        },
    )
}
#[test]
fn published_ranking_retains_original_candidates_and_never_defaults_a_selection() {
    let (inputs, request) = fixture(B, true);
    let before = serde_json::to_vec(&inputs).unwrap();
    let expected = rank::inspect_revision(&inputs, &request).unwrap();
    let report = rank::inspect(&inputs, &request).unwrap();
    report.validate_binding(&expected).unwrap();
    assert_eq!(report.metadata.enzyme_formulation, "BsaI-HFv2");
    assert_eq!(report.metadata.geometry_family, GoldenGateEnzyme::BsaI);
    assert_eq!(
        report.search.original_intent.target.sequence.as_str(),
        format!("{A}{B}")
    );
    assert_eq!(report.ranked.len(), 3);
    let known: Vec<_> = report
        .ranked
        .iter()
        .filter(|r| r.score.class != ScoreClass::Unknown)
        .collect();
    assert!(!known.is_empty());
    let unknown = report
        .ranked
        .iter()
        .position(|r| r.score.class == ScoreClass::Unknown)
        .unwrap();
    assert!(
        report.ranked[..unknown]
            .iter()
            .all(|r| r.score.class != ScoreClass::Unknown)
    );
    assert!(
        report.ranked[unknown..]
            .iter()
            .all(|r| r.score.class == ScoreClass::Unknown)
    );
    let original = report
        .search
        .candidates
        .iter()
        .find(|c| c.shifts == vec![1, 0])
        .unwrap();
    let selected = rank::select(&report, &expected, &original.id).unwrap();
    assert_eq!(selected.candidate_id, original.id);
    assert_eq!(selected.context_revision, report.search.context_revision);
    assert_ne!(selected.context_revision, report.context_revision);
    assert!(rank::select(&report, &expected, "").is_err());
    let materialized = fusion::materialize(&inputs, &request.design, &selected).unwrap();
    assert_eq!(
        serde_json::to_value(&materialized.request).unwrap(),
        serde_json::to_value(&request.design).unwrap(),
        "ranking must not rewrite generic recipe identity or formulation"
    );
    assert_eq!(
        materialized.original_intent.target.sequence.as_str(),
        format!("{A}{B}")
    );
    let actual = ligation_recipe::preview(
        &inputs,
        &materialized.original_intent.target,
        &materialized.recipe,
        collagene_core::assembly_recipe::CleavagePolicy::AssumedComplete,
    )
    .unwrap();
    assert_eq!(
        actual.verification.status,
        Status::ConsistentUnderAssumptions
    );
    assert_eq!(serde_json::to_vec(&inputs).unwrap(), before);
}
#[test]
fn fixed_linear_terminal_ends_participate_in_every_conditional_denominator() {
    let (inputs, request) = fixture(B, false);
    let report = rank::inspect(&inputs, &request).unwrap();
    assert_eq!(report.search.terminal_fusions, vec!["ATGA", "GGAT"]);
    let profile = profiles::load(&request.profile).unwrap();
    let prepared = profile.prepare().unwrap();
    for ranked in &report.ranked {
        let candidate = report
            .search
            .candidates
            .iter()
            .find(|c| c.id == ranked.candidate_id)
            .unwrap();
        let scope = fusion::end_scope(&report.search, &candidate.id).unwrap();
        assert_eq!(scope.intended_junctions, candidate.fusions);
        assert_eq!(scope.physical_competitors, vec!["ATGA", "ATCC"]);
        assert_eq!(
            rank::score_scope(&report, &report.context_revision, &candidate.id).unwrap(),
            scope
        );
        let detail = prepared
            .score_with_competitors(&scope.intended_junctions, &scope.physical_competitors)
            .unwrap();
        assert_eq!(ranked.score.overall_fidelity, detail.overall_fidelity);
        assert_eq!(
            ranked.score.log_overall_fidelity,
            detail.log_overall_fidelity
        );
        assert_eq!(
            ranked.score.pairs.len(),
            1,
            "only the actual intended junction is a multiplied event"
        );
        assert_eq!(ranked.score.pairs[0].estimate, detail.pairs[0].estimate);
        // Independent oracle: actual physical partners include both intended
        // orientations, the literal left end and literal right 5′ end. There
        // is no GGAT partner unless an intended endpoint independently supplies it.
        let end = &candidate.fusions[0];
        let complement = reverse(end);
        let physical = std::collections::BTreeSet::from([
            end.clone(),
            complement.clone(),
            "ATGA".into(),
            "ATCC".into(),
        ]);
        assert!(!physical.contains("GGAT"));
        assert_eq!(
            detail.physical_ends,
            physical.iter().cloned().collect::<Vec<_>>()
        );
        if end != &complement {
            let correct = u128::from(profile.weight(end, &complement).unwrap())
                + u128::from(profile.weight(&complement, end).unwrap());
            let mut total = 0u128;
            for orientation in [end, &complement] {
                for partner in &physical {
                    total += u128::from(profile.weight(orientation, partner).unwrap());
                }
            }
            let expected = correct as f64 / total as f64;
            assert_eq!(ranked.score.overall_fidelity, Some(expected));
            assert_eq!(
                ranked.score.pairs[0].estimate,
                profiles::Estimate::Known {
                    correct_weight: correct,
                    total_weight: total,
                    fidelity: expected
                }
            );
        }
    }
}
#[test]
fn outer_palindrome_is_only_a_competitor_and_no_join_has_no_probability() {
    let (inputs, mut request) = fixture(B, false);
    request.design.base.terminal_right_overhang = Some("GTAC".into());
    request.design.windows[0] = Window {
        junction: 0,
        min_shift: 0,
        max_shift: 0,
    };
    let report = rank::inspect(&inputs, &request).unwrap();
    assert_eq!(report.ranked.len(), 1);
    assert_ne!(
        report.ranked[0].score.class,
        ScoreClass::Unknown,
        "a palindromic singleton competitor is not an intended ambiguous pair"
    );
    request.design.base.base.fragments.truncate(1);
    request.design.base.rows.truncate(1);
    request.design.windows.clear();
    let single = rank::inspect(&inputs, &request).unwrap();
    assert_eq!(single.ranked.len(), 1);
    let scope = fusion::end_scope(&single.search, &single.ranked[0].candidate_id).unwrap();
    assert!(scope.intended_junctions.is_empty());
    assert_eq!(scope.physical_competitors, vec!["ATGA", "GTAC"]);
    assert_eq!(single.ranked[0].score.class, ScoreClass::Unknown);
    assert_eq!(single.ranked[0].score.overall_fidelity, None);
    assert_eq!(single.ranked[0].score.log_overall_fidelity, None);
    assert!(single.ranked[0].score.pairs.is_empty());
}

#[test]
fn explicit_formulation_geometry_profile_revision_and_sources_are_not_silently_substituted() {
    let (inputs, request) = fixture(B, true);
    let report = rank::inspect(&inputs, &request).unwrap();
    let mut changed = request.clone();
    changed.enzyme_formulation = "BsaI".into();
    assert!(rank::inspect_revision(&inputs, &changed).is_err());
    let mut changed = request.clone();
    changed.design.base.enzyme = GoldenGateEnzyme::SapI;
    assert!(rank::inspect_revision(&inputs, &changed).is_err());
    let mut missing = request.clone();
    missing.profile.profile_id = "missing-profile".into();
    assert!(rank::inspect_revision(&inputs, &missing).is_err());
    let mut stale = request.clone();
    stale.profile.profile_revision = "0".repeat(64);
    assert!(rank::inspect_revision(&inputs, &stale).is_err());
    let mut changed = inputs.clone();
    changed[0].molecule =
        Molecule::dna("changed original", &format!("{A}A"), Topology::Linear).unwrap();
    changed[0].revision = changed[0].molecule.revision();
    let current = rank::inspect_revision(&changed, &request).unwrap();
    assert!(rank::select(&report, &current, &report.ranked[0].candidate_id).is_err());
    let mut changed = report.clone();
    changed.metadata.conditions.buffer.push_str(" altered");
    assert!(changed.validate_binding(&report.context_revision).is_err());
    let esp3i = profile("Esp3I");
    let mut request = request;
    request.profile = esp3i.profile;
    request.enzyme_formulation = "Esp3I".into();
    request.design.base.enzyme = GoldenGateEnzyme::BsmBI;
    let explicit = rank::inspect(&inputs, &request).unwrap();
    assert_eq!(explicit.metadata.enzyme_formulation, "Esp3I");
    request.enzyme_formulation = "BsmBI-v2".into();
    assert!(rank::inspect_revision(&inputs, &request).is_err());
}
#[test]
fn equal_endsets_and_unknown_sets_tie_by_original_stable_candidate_id() {
    let (inputs, mut request) = fixture("AAAAAATGAGTTAGCGACCTTGGAA", true);
    request.design.windows[0] = Window {
        junction: 0,
        min_shift: 0,
        max_shift: 1,
    };
    let report = rank::inspect(&inputs, &request).unwrap();
    assert_eq!(report.ranked.len(), 2);
    assert_eq!(
        report.ranked[0].score.overall_fidelity,
        report.ranked[1].score.overall_fidelity
    );
    assert_eq!(
        report.ranked[0].score.log_overall_fidelity,
        report.ranked[1].score.log_overall_fidelity
    );
    assert!(report.ranked[0].candidate_id < report.ranked[1].candidate_id);
    let (inputs, mut request) = fixture("ATATATGAGTTAGCGACCTTGGAA", false);
    request.design.windows[0] = Window {
        junction: 0,
        min_shift: 0,
        max_shift: 2,
    };
    // These intended interfaces are ATAT/TATA/ATAT, all palindromic.
    // A palindromic outer competitor alone must never make a pair unknown.
    let report = rank::inspect(&inputs, &request).unwrap();
    assert!(
        report
            .ranked
            .iter()
            .all(|r| r.score.class == ScoreClass::Unknown)
    );
    let mut expected = report
        .search
        .candidates
        .iter()
        .map(|c| c.id.clone())
        .collect::<Vec<_>>();
    expected.sort();
    assert_eq!(
        report
            .ranked
            .iter()
            .map(|r| r.candidate_id.clone())
            .collect::<Vec<_>>(),
        expected
    );
}
fn reverse(s: &str) -> String {
    Sequence::dna(s)
        .unwrap()
        .reverse_complement()
        .unwrap()
        .as_str()
        .into()
}
fn single(
    estimate: profiles::Estimate,
    class: ScoreClass,
    overall: Option<f64>,
    log: Option<f64>,
) -> Score {
    Score {
        class,
        overall_fidelity: overall,
        log_overall_fidelity: log,
        pairs: vec![PairEstimate {
            canonical_end: "AAAC".into(),
            reverse_complement: "GTTT".into(),
            estimate,
        }],
    }
}
#[test]
fn compact_known_zero_missing_counts_and_palindromes_remain_distinct() {
    let empty = Score {
        class: ScoreClass::Unknown,
        overall_fidelity: None,
        log_overall_fidelity: None,
        pairs: Vec::new(),
    };
    empty.validate(&[], 4).unwrap();
    let false_identity = Score {
        class: ScoreClass::Positive,
        overall_fidelity: Some(1.0),
        log_overall_fidelity: Some(0.0),
        pairs: Vec::new(),
    };
    assert!(false_identity.validate(&[], 4).is_err());

    // Synthetic compact arithmetic fixtures; these do not claim to be observations
    // from an embedded dataset or native source-bound ranking receipts.
    let zero = single(
        profiles::Estimate::Known {
            correct_weight: 0,
            total_weight: 7,
            fidelity: 0.0,
        },
        ScoreClass::Zero,
        Some(0.0),
        None,
    );
    zero.validate(&["AAAC".into()], 4).unwrap();
    let missing = single(
        profiles::Estimate::Unknown {
            reason: UnknownReason::MissingCounts,
        },
        ScoreClass::Unknown,
        None,
        None,
    );
    missing.validate(&["AAAC".into()], 4).unwrap();
    let mut false_zero = missing;
    false_zero.class = ScoreClass::Zero;
    false_zero.overall_fidelity = Some(0.0);
    assert!(false_zero.validate(&["AAAC".into()], 4).is_err());
    let palindrome = Score {
        class: ScoreClass::Unknown,
        overall_fidelity: None,
        log_overall_fidelity: None,
        pairs: vec![PairEstimate {
            canonical_end: "GTAC".into(),
            reverse_complement: "GTAC".into(),
            estimate: profiles::Estimate::Unknown {
                reason: UnknownReason::PalindromicEnd,
            },
        }],
    };
    palindrome.validate(&["GTAC".into()], 4).unwrap();
}
#[test]
fn positive_products_that_underflow_are_still_known_positive_through_their_logs() {
    let mut ends = Vec::new();
    'bases: for a in b"ACGT" {
        for b in b"ACGT" {
            for c in b"ACGT" {
                for d in b"ACGT" {
                    let end = String::from_utf8(vec![*a, *b, *c, *d]).unwrap();
                    if end < reverse(&end) {
                        ends.push(end);
                        if ends.len() == 16 {
                            break 'bases;
                        }
                    }
                }
            }
        }
    }
    let total = u128::from(u64::MAX) * 64;
    let fidelity = 1.0 / total as f64;
    let score = Score {
        class: ScoreClass::Positive,
        overall_fidelity: Some(ends.iter().fold(1.0, |v, _| v * fidelity)),
        log_overall_fidelity: Some(ends.iter().fold(0.0, |v, _| v + fidelity.ln())),
        pairs: ends
            .iter()
            .map(|end| PairEstimate {
                canonical_end: end.clone(),
                reverse_complement: reverse(end),
                estimate: profiles::Estimate::Known {
                    correct_weight: 1,
                    total_weight: total,
                    fidelity,
                },
            })
            .collect(),
    };
    assert_eq!(score.overall_fidelity, Some(0.0));
    assert!(score.log_overall_fidelity.unwrap().is_finite());
    score.validate(&ends, 4).unwrap();
    let mut mislabeled = score;
    mislabeled.class = ScoreClass::Zero;
    assert!(mislabeled.validate(&ends, 4).is_err());
}
#[test]
fn incomplete_searches_and_tampered_rank_proofs_do_not_become_global_optima() {
    let (inputs, mut request) = fixture(B, true);
    request.design.limits.max_states = 1;
    let report = rank::inspect(&inputs, &request).unwrap();
    assert!(!report.search.search_complete);
    assert_eq!(report.ranked.len(), 1);
    assert!(!report.limits.is_empty());
    let mut altered = report.clone();
    altered.ranked[0].rank = 2;
    assert!(altered.validate_binding(&report.context_revision).is_err());
    let mut altered = report.clone();
    altered.ranked[0].candidate_id = "0".repeat(64);
    assert!(altered.validate_binding(&report.context_revision).is_err());
    let mut altered = report.clone();
    altered.metadata.data_revision = "not-a-hash".into();
    assert!(altered.validate_binding(&report.context_revision).is_err());
    request.design.limits.max_candidates = fusion::MAX_CANDIDATES + 1;
    assert!(rank::inspect_revision(&inputs, &request).is_err());
}
