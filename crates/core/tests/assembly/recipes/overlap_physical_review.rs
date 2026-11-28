//! Independent review: exact homology cannot remove a known blocked3′ group.
use collagene_core::{
    chemistry::{AnnealOptions, Terminus, ThreePrime},
    cloning::overlap::profiles::ProfileOptions,
    planning::plan::{DesignPlan, Operation, OrientedSource, PlanInput, Source},
    verification::{engine::verify, evidence::Status},
};
use collagene_model::{Molecule, Topology};
fn input(id: &str, bases: &str) -> PlanInput {
    let molecule = Molecule::dna(id, bases, Topology::Linear).unwrap();
    PlanInput {
        id: id.into(),
        revision: molecule.revision(),
        molecule,
    }
}
fn anneal(top: &str, bottom: &str, top_end: Terminus, bottom_end: Terminus) -> Operation {
    Operation::AnnealOligos {
        top: Source::Input(top.into()),
        bottom: Source::Input(bottom.into()),
        options: AnnealOptions {
            min_paired_bases: 6,
            top_chemistry: top_end,
            bottom_chemistry: bottom_end,
            ..Default::default()
        },
        candidate_index: 0,
    }
}
fn candidate(profiled: bool, top: Terminus, bottom: Terminus) -> (DesignPlan, Molecule) {
    let target = Molecule::dna("literal target", "AAAACCGGTT", Topology::Linear).unwrap();
    let operation = if profiled {
        Operation::ProfiledOverlap {
            fragments: vec![
                OrientedSource {
                    source: Source::Step(0),
                    reverse: false,
                },
                OrientedSource {
                    source: Source::Step(1),
                    reverse: false,
                },
            ],
            overlaps: vec![2],
            circular: false,
            options: ProfileOptions {
                scan_repeated_overlaps: false,
                ..Default::default()
            },
        }
    } else {
        Operation::Overlap {
            fragments: vec![Source::Step(0), Source::Step(1)],
            overlaps: vec![2],
            circular: false,
        }
    };
    (
        DesignPlan {
            name: "blocked end review".into(),
            target_revision: target.revision(),
            inputs: vec![
                input("a", "AAAACC"),
                input("b", "GGTTTT"),
                input("c", "CCGGTT"),
                input("d", "AACCGG"),
            ],
            steps: vec![
                anneal("a", "b", top, bottom),
                anneal("c", "d", Terminus::CLEAVED, Terminus::CLEAVED),
                operation,
            ],
        },
        target,
    )
}
#[test]
fn known_blocked_top_or_bottom_three_prime_end_never_disappears_in_homology_assembly() {
    let blocked = Terminus {
        three_prime: ThreePrime::Blocked,
        ..Terminus::CLEAVED
    };
    for profiled in [false, true] {
        for (top, bottom) in [(blocked, Terminus::CLEAVED), (Terminus::CLEAVED, blocked)] {
            let (plan, target) = candidate(profiled, top, bottom);
            let report = verify(&plan, &target);
            assert_eq!(
                report.status,
                Status::Incomplete,
                "profiled={profiled}: {:?}",
                report.findings
            );
            assert_eq!(
                report.products.len(),
                2,
                "the unsupported assembly must not be emitted as a reconstructed final product"
            );
            assert!(
                report
                    .findings
                    .iter()
                    .any(|f| f.contains("blocked") || f.contains("Blocked"))
            );
        }
    }
}
#[test]
fn accessible_hydroxyls_still_allow_the_literal_exact_candidate() {
    for profiled in [false, true] {
        let (plan, target) = candidate(profiled, Terminus::CLEAVED, Terminus::CLEAVED);
        let report = verify(&plan, &target);
        assert_eq!(
            report.status,
            Status::ConsistentUnderAssumptions,
            "{:?}",
            report.findings
        );
        assert_eq!(
            report.products.last().unwrap().sequence.as_str(),
            "AAAACCGGTT"
        );
    }
}
