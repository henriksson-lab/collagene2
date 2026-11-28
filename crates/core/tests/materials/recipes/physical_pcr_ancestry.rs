//! A sequence-only assembly is not a physical PCR template. Independent target
//! agreement must not repair missing strands, termini or modification ancestry.
use collagene_core::{
    analysis::Primer,
    chemistry::{DuplexSummary, EndTreatment},
    modified_duplex::ModificationSummary,
    mutagenesis::PcrOptions,
    overlap_profiles::ProfileOptions,
    tailing::TailingBase,
    verification::{self, DesignPlan, Operation, OrientedSource, PlanInput, Source, Status},
};
use collagene_model::{Molecule, Topology};

// Hand-composed end homology: AGTCCGATGCTA[CGTA] + [CGTA]TTCAGACCTGTA.
const EXPECTED: &str = "AGTCCGATGCTACGTATTCAGACCTGTA";
const MUTATED: &str = "TGTCCGATGCTACGTATTCAGACCTGTA";
fn input(id: &str, text: &str) -> PlanInput {
    let molecule = Molecule::dna(id, text, Topology::Linear).unwrap();
    PlanInput {
        id: id.into(),
        revision: molecule.revision(),
        molecule,
    }
}
fn primers(mutate: bool) -> (Primer, Primer) {
    (
        Primer::new(
            "Outer forward",
            if mutate { "TGTCCGAT" } else { "AGTCCGAT" },
            "",
        )
        .unwrap(),
        // Reverse complement of the independently specified final eight bases GACCTGTA.
        Primer::new("Outer reverse", "TACAGGTC", "").unwrap(),
    )
}
fn pcr(template: Source, mismatches: bool, mutate: bool) -> Operation {
    let (forward, reverse) = primers(mutate);
    if mismatches {
        let mut options = PcrOptions::default();
        options.binding.max_mismatches = usize::from(mutate);
        options.binding.min_three_prime_exact = 7;
        Operation::PcrWithMismatches {
            template,
            forward,
            reverse,
            options,
        }
    } else {
        assert!(!mutate);
        Operation::Pcr {
            template,
            forward,
            reverse,
            template_strand: collagene_core::chemistry::Strand::Top,
        }
    }
}
fn assembly(profiled: bool) -> Operation {
    if profiled {
        Operation::ProfiledOverlap {
            fragments: vec![
                OrientedSource {
                    source: Source::Input("left".into()),
                    reverse: false,
                },
                OrientedSource {
                    source: Source::Input("right".into()),
                    reverse: false,
                },
            ],
            overlaps: vec![4],
            circular: false,
            options: ProfileOptions::default(),
        }
    } else {
        Operation::Overlap {
            fragments: vec![Source::Input("left".into()), Source::Input("right".into())],
            overlaps: vec![4],
            circular: false,
        }
    }
}
fn plan(target: &Molecule, inputs: Vec<PlanInput>, steps: Vec<Operation>) -> DesignPlan {
    DesignPlan {
        name: "Explicit physical PCR ancestry".into(),
        target_revision: target.revision(),
        inputs,
        steps,
    }
}

#[test]
fn sequence_only_overlap_outputs_remain_diagnostic_but_cannot_supply_either_pcr_operation() {
    let target = Molecule::dna("Independent reference target", EXPECTED, Topology::Linear).unwrap();
    for profiled in [false, true] {
        let inputs = vec![
            input("left", "AGTCCGATGCTACGTA"),
            input("right", "CGTATTCAGACCTGTA"),
        ];
        let prefix = plan(&target, inputs, vec![assembly(profiled)]);
        let reference = verification::verify(&prefix, &target);
        assert_eq!(
            reference.status,
            Status::ConsistentUnderAssumptions,
            "{:?}",
            reference.findings
        );
        assert_eq!(reference.products.len(), 1);
        assert_eq!(reference.products[0].sequence.as_str(), EXPECTED);
        assert!(!reference.products[0].metadata.contains_key("duplex_model"));
        assert!(
            !reference.products[0]
                .metadata
                .contains_key("tracked_modifications")
        );
        for mismatches in [false, true] {
            let mut attempted = prefix.clone();
            attempted
                .steps
                .push(pcr(Source::Step(0), mismatches, false));
            let before = serde_json::to_value(&attempted).unwrap();
            let report = verification::verify(&attempted, &target);
            assert_eq!(report.status, Status::Incomplete, "{:?}", report.findings);
            assert_eq!(
                report.products.len(),
                1,
                "No PCR product may be invented from the reference prefix"
            );
            assert_eq!(report.products[0].sequence.as_str(), EXPECTED);
            assert!(!report.products[0].metadata.contains_key("duplex_model"));
            assert!(
                !report.products[0]
                    .metadata
                    .contains_key("tracked_modifications")
            );
            assert!(
                report.findings.iter().any(|f| f.contains(if mismatches {
                    "established physical template"
                } else {
                    "actual physical source state"
                })),
                "{:?}",
                report.findings
            );
            assert_eq!(serde_json::to_value(&attempted).unwrap(), before);
            assert_eq!(report.target_revision, target.revision());
        }
    }
}

#[test]
fn original_reference_input_retains_its_explicit_intact_template_hypothesis() {
    let target = Molecule::dna("Independent PCR target", EXPECTED, Topology::Linear).unwrap();
    for mismatches in [false, true] {
        let plan = plan(
            &target,
            vec![input("original", EXPECTED)],
            vec![pcr(Source::Input("original".into()), mismatches, false)],
        );
        let report = verification::verify(&plan, &target);
        assert_eq!(
            report.status,
            Status::ConsistentUnderAssumptions,
            "{:?}",
            report.findings
        );
        assert_eq!(report.products.len(), 1);
        assert_eq!(report.products[0].sequence.as_str(), EXPECTED);
        let duplex: DuplexSummary =
            serde_json::from_str(&report.products[0].metadata["duplex_model"]).unwrap();
        assert_eq!(duplex.topology, Topology::Linear);
        assert!(duplex.nicks.is_empty());
    }
}

#[test]
fn reconstructed_physical_pcr_feeds_exact_or_mutagenic_pcr_and_real_terminal_chemistry() {
    for mismatches in [false, true] {
        let expected = if mismatches { MUTATED } else { EXPECTED };
        let target =
            Molecule::dna("Independent sealed construct", expected, Topology::Circular).unwrap();
        let plan = plan(
            &target,
            vec![input("original", EXPECTED)],
            vec![
                pcr(Source::Input("original".into()), false, false),
                pcr(Source::Step(0), mismatches, mismatches),
                Operation::TreatEnds {
                    template: Source::Step(1),
                    treatment: EndTreatment::Phosphorylate,
                },
                Operation::Ligate {
                    fragments: vec![OrientedSource {
                        source: Source::Step(2),
                        reverse: false,
                    }],
                    circular: true,
                    allow_nicked: false,
                    options: Default::default(),
                },
            ],
        );
        let before = serde_json::to_value(&plan).unwrap();
        let report = verification::verify(&plan, &target);
        assert_eq!(
            report.status,
            Status::ConsistentUnderAssumptions,
            "{:?}",
            report.findings
        );
        assert_eq!(report.products.len(), 4);
        assert_eq!(report.products[0].sequence.as_str(), EXPECTED);
        assert_eq!(report.products[1].sequence.as_str(), expected);
        let final_product = report.products.last().unwrap();
        assert_eq!(final_product.sequence.as_str(), expected);
        assert_eq!(final_product.topology, Topology::Circular);
        let duplex: DuplexSummary =
            serde_json::from_str(&final_product.metadata["duplex_model"]).unwrap();
        assert!(duplex.covalently_closed_in_model);
        assert!(duplex.nicks.is_empty());
        let modifications: ModificationSummary =
            serde_json::from_str(&report.products[1].metadata["tracked_modifications"]).unwrap();
        assert_eq!(modifications.known_modified, 0);
        assert_eq!(modifications.unknown_eligible, 0);
        assert!(modifications.known_unmodified > 0);
        assert_eq!(serde_json::to_value(&plan).unwrap(), before);
    }
}

#[test]
fn imported_projection_metadata_is_not_reinterpreted_as_an_original_physical_template() {
    let target = Molecule::dna("Independent PCR target", EXPECTED, Topology::Linear).unwrap();
    for mismatches in [false, true] {
        let mut source = input("claimed-projection", EXPECTED);
        source.molecule.metadata.insert(
            "duplex_model".into(),
            "Reference-only imported projection; physical state unavailable".into(),
        );
        source.revision = source.molecule.revision();
        let plan = plan(
            &target,
            vec![source],
            vec![pcr(
                Source::Input("claimed-projection".into()),
                mismatches,
                false,
            )],
        );
        let report = verification::verify(&plan, &target);
        assert_eq!(report.status, Status::Incomplete, "{:?}", report.findings);
        assert!(report.products.is_empty());
        assert!(
            report.findings.iter().any(|f| f.contains(if mismatches {
                "imported physical"
            } else {
                "fresh native origin"
            })),
            "{:?}",
            report.findings
        );
    }
}

#[test]
fn exact_pcr_uses_covered_staggered_template_but_mismatch_profile_remains_unsupported() {
    let target = Molecule::dna("Independent PCR target", EXPECTED, Topology::Linear).unwrap();
    for mismatches in [false, true] {
        let plan = plan(
            &target,
            vec![input("original", EXPECTED)],
            vec![
                pcr(Source::Input("original".into()), false, false),
                Operation::Tail3Prime {
                    template: Source::Step(0),
                    base: TailingBase::A,
                },
                pcr(Source::Step(1), mismatches, false),
            ],
        );
        let report = verification::verify(&plan, &target);
        // Tailing changes the actual template; subsequent PCR must neither
        // overwrite that parent nor inherit its unpaired enzymatic tails.
        assert_eq!(report.products[1].sequence.as_str(), format!("{EXPECTED}A"));
        let template: DuplexSummary =
            serde_json::from_str(&report.products[1].metadata["duplex_model"]).unwrap();
        assert_eq!(template.top_length, EXPECTED.len() + 1);
        assert_eq!(template.bottom_length, EXPECTED.len() + 1);
        assert_eq!(template.bottom_start, -1);
        if mismatches {
            assert_eq!(report.status, Status::Incomplete, "{:?}", report.findings);
            assert_eq!(report.products.len(), 2);
            assert!(
                report.findings.iter().any(|f| f.contains("staggered")),
                "{:?}",
                report.findings
            );
        } else {
            assert_eq!(
                report.status,
                Status::ConsistentUnderAssumptions,
                "{:?}",
                report.findings
            );
            assert_eq!(report.products.len(), 3);
            assert_eq!(report.products[2].sequence.as_str(), EXPECTED);
            let product: DuplexSummary =
                serde_json::from_str(&report.products[2].metadata["duplex_model"]).unwrap();
            assert_eq!(product.top_length, EXPECTED.len());
            assert_eq!(product.bottom_length, EXPECTED.len());
            assert_eq!(product.bottom_start, 0);
            assert!(product.nicks.is_empty());
        }
    }
}
