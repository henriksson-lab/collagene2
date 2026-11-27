//! Independent signed partitions and actual PCR/digest/ligation oracles.
use collagene_core::{
    Molecule, Sequence, Topology, assembly_intent,
    assembly_junction_edit::Edit,
    golden_gate::GoldenGateEnzyme,
    golden_gate_fusion_design::{self as fusion, Constraints, Limits, Request, Selection, Window},
    golden_gate_junction_edit::{self as gate, Row},
    ligation_recipe,
    overlap_primer_design::Fragment,
    verification::{PlanInput, Status},
};
const A: &str = "ATGACCCCGGGGTTTCCAGA";
const B: &str = "ATGAGTTAGCGACCTTGGAA";
fn fixture(a: &str, b: &str, reverse_b: bool, circular: bool) -> (Vec<PlanInput>, Request) {
    let b_source = if reverse_b {
        Sequence::dna(b)
            .unwrap()
            .reverse_complement()
            .unwrap()
            .as_str()
            .to_owned()
    } else {
        b.into()
    };
    let inputs = [("a", a.to_owned()), ("b", b_source)]
        .into_iter()
        .map(|(id, dna)| {
            let molecule = Molecule::dna(id, &dna, Topology::Linear).unwrap();
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
            end: a.len(),
            reverse: false,
        },
        Fragment {
            role: "right".into(),
            input_id: "b".into(),
            start: 0,
            end: b.len(),
            reverse: reverse_b,
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
    (
        inputs,
        Request {
            base: gate::Request {
                base: assembly_intent::Request {
                    name: "independent literal AB".into(),
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
            windows: vec![],
            constraints: Constraints {
                reject_palindromes: false,
                distinct_reverse_complements: false,
                minimum_hamming: 0,
            },
            limits: Limits::default(),
        },
    )
}
fn choose(inputs: &[PlanInput], request: &Request, shifts: &[i32]) -> fusion::Materialized {
    let report = fusion::inspect(inputs, request).unwrap();
    let candidate = report
        .candidates
        .iter()
        .find(|c| c.shifts == shifts)
        .expect("declared signed placement");
    fusion::materialize(
        inputs,
        request,
        &Selection {
            context_revision: report.context_revision.clone(),
            candidate_id: candidate.id.clone(),
        },
    )
    .unwrap()
}
fn physical(inputs: &[PlanInput], m: &fusion::Materialized) {
    let report = ligation_recipe::preview(
        inputs,
        &m.original_intent.target,
        &m.recipe,
        collagene_core::assembly_recipe::CleavagePolicy::AssumedComplete,
    )
    .unwrap();
    assert_eq!(
        report.verification.status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        report.verification.findings
    );
    let actual = report.verification.products.last().unwrap();
    let declared = m
        .placements
        .iter()
        .map(|p| p.oriented_top.as_str())
        .collect::<String>();
    assert_eq!(actual.sequence.as_str(), declared);
}
#[test]
fn signed_linear_boundaries_transfer_literal_bases_without_changing_independent_target() {
    let (inputs, mut request) = fixture(A, B, false, false);
    request.windows = vec![Window {
        junction: 0,
        min_shift: -1,
        max_shift: 1,
    }];
    let snapshot = serde_json::to_vec(&inputs).unwrap();
    let left = choose(&inputs, &request, &[-1]);
    assert_eq!(
        left.original_intent.target.sequence.as_str(),
        format!("{A}{B}")
    );
    assert_eq!(left.placements[0].oriented_top, &A[..A.len() - 1]);
    assert_eq!(left.placements[1].prefix, &A[A.len() - 1..]);
    assert_eq!(
        left.placements[1].oriented_top,
        format!("{}{B}", &A[A.len() - 1..])
    );
    assert_eq!(
        left.candidate.fusions,
        vec![format!("{}{}", &A[A.len() - 1..], &B[..3])]
    );
    assert_eq!(
        left.terminal_fusions,
        vec![A[..4].to_owned(), "GGAT".into()]
    );
    physical(&inputs, &left);
    let right = choose(&inputs, &request, &[1]);
    assert_eq!(right.placements[0].suffix, &B[..1]);
    assert_eq!(right.placements[1].start, 1);
    assert_eq!(right.candidate.fusions, vec![B[1..5].to_owned()]);
    physical(&inputs, &right);
    assert_eq!(serde_json::to_vec(&inputs).unwrap(), snapshot);
}
#[test]
fn reverse_source_prefix_and_circle_closing_moves_use_actual_opposite_strand_geometry() {
    let a = format!("{A}TCGC");
    let b = format!("{B}ACGT");
    let (inputs, mut request) = fixture(&a, &b, true, true);
    request.windows = vec![
        Window {
            junction: 0,
            min_shift: -1,
            max_shift: 1,
        },
        Window {
            junction: 1,
            min_shift: -1,
            max_shift: 1,
        },
    ];
    for shifts in [[-1, 0], [1, 0], [0, -1], [0, 1]] {
        let materialized = choose(&inputs, &request, &shifts);
        assert_eq!(
            materialized.original_intent.target.sequence.as_str(),
            format!("{a}{b}")
        );
        if shifts == [-1, 0] {
            assert_eq!(materialized.placements[1].prefix, &a[a.len() - 1..]);
        }
        if shifts == [0, -1] {
            assert_eq!(materialized.placements[0].prefix, &b[b.len() - 1..]);
        }
        if shifts == [0, 1] {
            assert_eq!(materialized.placements[1].suffix, &a[..1]);
        }
        physical(&inputs, &materialized);
    }
}
#[test]
fn previous_trims_and_inserts_remain_literal_when_a_boundary_moves_left() {
    let a = format!("{A}TCGC");
    let b = format!("{B}ACGT");
    let (inputs, mut request) = fixture(&a, &b, true, false);
    request.base.edits = vec![Edit {
        junction: 0,
        remove_left: 1,
        remove_right: 2,
        inserted: "GT".into(),
    }];
    request.windows = vec![Window {
        junction: 0,
        min_shift: -3,
        max_shift: -3,
    }];
    let m = choose(&inputs, &request, &[-3]);
    assert_eq!(
        m.original_intent.target.sequence.as_str(),
        format!("{}GT{}", &a[..a.len() - 1], &b[2..])
    );
    assert_eq!(m.placements[0].oriented_top, &a[..a.len() - 2]);
    assert_eq!(
        m.placements[1].prefix,
        format!("{}GT", &a[a.len() - 2..a.len() - 1])
    );
    assert_eq!(m.original_intent.edits[0].inserted, "GT");
    assert_eq!(m.original_intent.edits[0].remove_right, 2);
    physical(&inputs, &m);
}
#[test]
fn manual_preparations_are_fixed_and_conflicting_windows_are_safely_pruned() {
    let (inputs, mut request) = fixture(A, B, false, true);
    let manual = gate::design(&inputs, &request.base)
        .unwrap()
        .recipe
        .fragments[0]
        .clone();
    request.base.rows[0] = Row::Retain(manual.clone());
    request.windows = vec![
        Window {
            junction: 0,
            min_shift: -128,
            max_shift: 128,
        },
        Window {
            junction: 1,
            min_shift: -128,
            max_shift: 128,
        },
    ];
    let report = fusion::inspect(&inputs, &request).unwrap();
    assert!(report.search_complete);
    assert_eq!(report.states_examined, 1);
    assert_eq!(report.candidates.len(), 1);
    let m = choose(&inputs, &request, &[0, 0]);
    assert!(m.designs[0].is_none());
    assert_eq!(
        serde_json::to_vec(&m.recipe.fragments[0]).unwrap(),
        serde_json::to_vec(&manual).unwrap()
    );
    physical(&inputs, &m);
    request.windows[0].min_shift = 1;
    let report = fusion::inspect(&inputs, &request).unwrap();
    assert!(report.search_complete);
    assert_eq!(report.states_examined, 0);
    assert!(report.candidates.is_empty());
}
#[test]
fn heuristic_constraints_include_fixed_terminal_ends_and_never_claim_empirical_success() {
    let (inputs, mut request) = fixture(A, B, false, false);
    let baseline = fusion::inspect(&inputs, &request).unwrap();
    assert_eq!(baseline.candidates.len(), 1);
    // Original left free end and the single intended junction both spell ATGA.
    request.constraints.distinct_reverse_complements = true;
    let constrained = fusion::inspect(&inputs, &request).unwrap();
    assert!(constrained.search_complete);
    assert!(constrained.candidates.is_empty());
    request.constraints.distinct_reverse_complements = false;
    request.constraints.minimum_hamming = 1;
    assert!(
        fusion::inspect(&inputs, &request)
            .unwrap()
            .candidates
            .is_empty()
    );
    request.constraints.minimum_hamming = 0;
    request.base.terminal_right_overhang = Some("ATAT".into());
    request.constraints.reject_palindromes = true;
    assert!(
        fusion::inspect(&inputs, &request)
            .unwrap()
            .candidates
            .is_empty()
    );
}
#[test]
fn caps_stale_context_and_tampered_geometry_withhold_authority() {
    let (inputs, mut request) = fixture(A, B, false, true);
    request.windows = vec![
        Window {
            junction: 0,
            min_shift: -128,
            max_shift: 128,
        },
        Window {
            junction: 1,
            min_shift: -128,
            max_shift: 128,
        },
    ];
    request.limits.max_states = 1;
    let report = fusion::inspect(&inputs, &request).unwrap();
    assert!(!report.search_complete);
    assert_eq!(report.states_examined, 1);
    assert_eq!(report.candidates[0].shifts, vec![0, 0]);
    assert!(!report.limits.is_empty());
    let selection = Selection {
        context_revision: report.context_revision.clone(),
        candidate_id: report.candidates[0].id.clone(),
    };
    let mut changed = request.clone();
    changed.limits.max_states = 2;
    assert!(fusion::materialize(&inputs, &changed, &selection).is_err());
    let mut forged = report.clone();
    forged.candidates[0].fusions[0] = "AAAA".into();
    assert!(forged.validate_binding(&report.context_revision).is_err());
    let mut m = fusion::materialize(&inputs, &request, &selection).unwrap();
    m.placements[0].suffix.push('A');
    let context = m.context_revision.clone();
    assert!(m.validate_binding(&context).is_err());
    let mut impossible = request;
    impossible.limits.max_states = 10;
    impossible.windows = vec![
        Window {
            junction: 0,
            min_shift: 8,
            max_shift: 8,
        },
        Window {
            junction: 1,
            min_shift: -8,
            max_shift: -8,
        },
    ];
    let report = fusion::inspect(&inputs, &impossible).unwrap();
    assert!(report.search_complete);
    assert!(
        report.candidates.is_empty(),
        "interacting shifts cannot erase an original annealing core"
    );
}
#[test]
fn redistribution_can_rescue_an_old_insert_but_actual_oligo_limits_still_apply() {
    let (inputs, mut request) = fixture(A, B, false, false);
    request.base.edits = vec![Edit {
        junction: 0,
        remove_left: 0,
        remove_right: 0,
        inserted: "G".repeat(200),
    }];
    assert!(
        gate::preflight(
            &inputs
                .iter()
                .map(|i| (i.id.as_str(), &i.molecule))
                .collect::<Vec<_>>(),
            &request.base
        )
        .is_err()
    );
    request.windows = vec![Window {
        junction: 0,
        min_shift: -23,
        max_shift: -23,
    }];
    let m = choose(&inputs, &request, &[-23]);
    assert_eq!(m.placements[1].prefix, "G".repeat(23));
    assert_eq!(
        m.designs[0]
            .as_ref()
            .unwrap()
            .reverse
            .ordered_sequence()
            .len(),
        200
    );
    physical(&inputs, &m);
    request.windows[0] = Window {
        junction: 0,
        min_shift: -22,
        max_shift: -22,
    };
    let report = fusion::inspect(&inputs, &request).unwrap();
    assert!(report.search_complete);
    assert!(
        report.candidates.is_empty(),
        "201-base actual oligo is not a valid placement"
    );
}
