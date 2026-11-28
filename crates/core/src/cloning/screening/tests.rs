use super::*;
use crate::{
    chemistry::{
        Strand,
        methylation::MethylationState,
        tracked_duplex::{GenerationAssumptions, checked_selection},
    },
    material::derivation::SourceOrigin,
    sequence::analysis::Primer,
};
use collagene_model::{Error, Molecule, Topology};
use std::collections::BTreeSet;

fn component(id: &str, sequence: &str, topology: Topology) -> Component {
    let molecule = Molecule::dna(id, sequence, topology).unwrap();
    Component {
        id: id.into(),
        label: id.into(),
        revision: molecule.revision(),
        molecule,
        origin: SourceOrigin::OriginalReference,
        original_generation: Some(GenerationAssumptions {
            host: None,
            dam: MethylationState::Unmodified,
            dcm: MethylationState::Unmodified,
            cpg: MethylationState::Unmodified,
        }),
    }
}
fn subject(id: &str, role: Role, components: &[&str]) -> Subject {
    Subject {
        id: id.into(),
        label: id.into(),
        role,
        components: components.iter().map(|s| (*s).into()).collect(),
    }
}
fn pcr_test(id: &str) -> Test {
    Test {
        id: id.into(),
        label: id.into(),
        method: Method::ExactPcr {
            forward: Primer::new("forward", "ACGTTGCA", "").unwrap(),
            reverse: Primer::new("reverse", "TGCATCGG", "").unwrap(),
            template_strand: Strand::Top,
        },
    }
}
fn digest_request() -> Request {
    Request {
        components: vec![
            component("intended", "AAAAGAATTCTTTT", Topology::Linear),
            component("parent", "CCCCCCCCCCCCCCCCCCCC", Topology::Linear),
        ],
        subjects: vec![
            subject("intended", Role::Intended, &["intended"]),
            subject(
                "parent carryover mixture",
                Role::Decoy(DecoyKind::ParentCarryover),
                &["intended", "parent"],
            ),
        ],
        tests: vec![Test {
            id: "EcoRI".into(),
            label: "EcoRI".into(),
            method: Method::Digest {
                selection: checked_selection("EcoRI").unwrap(),
            },
        }],
        readout: IdealSizeReadout {
            counterpart_tolerance_twice_bp: 0,
        },
        limits: Limits::default(),
    }
}

#[test]
fn literal_digest_sizes_combine_named_mixture_without_intensity_or_band_count_inference() {
    // Literal EcoRI site starts at4; top bond5/bottom9 on a14-base source.
    // Two strand lengths5+9 and9+5 each give one identical twice-mean size14.
    // The20-base uncut LINEAR parent gives size40, a separate mixture witness.
    let request = digest_request();
    let report = assess(&request).unwrap();
    report
        .validate_binding(&revision(&request).unwrap())
        .unwrap();
    assert_eq!(report.cells[0].native_product_count, 2);
    assert_eq!(report.cells[0].sizes_twice_bp, [14]);
    assert_eq!(report.cells[1].sizes_twice_bp, [40]);
    assert_eq!(report.tests[0].signatures[0].sizes_twice_bp, [14]);
    assert_eq!(report.tests[0].signatures[1].sizes_twice_bp, [14, 40]);
    assert_eq!(
        report.tests[0].comparisons[0].distinction,
        Distinction::Separated(Witness {
            side: WitnessSide::Decoy,
            size_twice_bp: 40
        })
    );
    assert_eq!(report.panel.test_indices, [0]);
    assert!(report.panel.complete_through_size_three);
    assert!(report.cleavage_work > 0);
}

#[test]
fn identical_dna_named_outcomes_are_retained_and_honestly_nondiscriminable() {
    let mut request = digest_request();
    request.components[1] = component("parent", "AAAAGAATTCTTTT", Topology::Linear);
    request.subjects[1] = subject(
        "independently named empty vector",
        Role::Decoy(DecoyKind::EmptyVector),
        &["parent"],
    );
    let report = assess(&request).unwrap();
    assert_eq!(report.components.len(), 2);
    assert_ne!(report.components[0].id, report.components[1].id);
    assert_eq!(report.subjects[1].label, "independently named empty vector");
    assert_eq!(
        report.tests[0].comparisons[0].distinction,
        Distinction::IdenticalSignature
    );
    assert!(report.panel.separated_decoy_indices.is_empty());
}

#[test]
fn complete_empty_pcr_signature_differs_from_unsupported_uncut_circular_digest() {
    let mut request = digest_request();
    request.components[0] = component("intended", "ACGTTGCACCGGAATTCCGATGCA", Topology::Linear);
    request.tests = vec![pcr_test("PCR")];
    let report = assess(&request).unwrap();
    assert_eq!(report.cells[0].state, CellState::Complete);
    assert!(report.cells[0].enumeration_complete);
    assert_eq!(report.cells[0].sizes_twice_bp, [48]);
    assert_eq!(report.cells[1].state, CellState::Complete);
    assert!(report.cells[1].enumeration_complete);
    assert_eq!(report.cells[1].native_product_count, 0);
    assert!(report.cells[1].sizes_twice_bp.is_empty());
    assert!(matches!(
        report.tests[0].comparisons[0].distinction,
        Distinction::IdenticalSignature
    ));
    // The decoy here is intended+parent: parent contributes no additional PCR size.
    let mut circular = digest_request();
    circular.components[1] = component("parent", "CCCCCCCCCCCCCCCCCCCC", Topology::Circular);
    let report = assess(&circular).unwrap();
    assert_eq!(report.cells[1].state, CellState::Unsupported);
    assert_eq!(
        report.tests[0].comparisons[0].distinction,
        Distinction::Unresolved
    );
}

#[test]
fn native_multiple_pcr_geometries_are_finite_not_unique_or_experimentally_guaranteed() {
    let mut request = digest_request();
    request.components[0] = component(
        "intended",
        "ACGTTGCACCGATGCAACGTTGCACCGATGCA",
        Topology::Linear,
    );
    request.tests = vec![pcr_test("PCR")];
    let report = assess(&request).unwrap();
    assert_eq!(report.cells[0].state, CellState::Complete);
    assert!(report.cells[0].enumeration_complete);
    assert!(report.cells[0].native_product_count > 1);
    assert_eq!(
        report.cells[0].native_product_count,
        report.cells[0].pcr_geometries.len()
    );
    assert!(
        report
            .statement
            .contains("not a claim that they amplify experimentally")
    );
    request.limits.max_products_per_cell = 1;
    let report = assess(&request).unwrap();
    assert_eq!(report.cells[0].state, CellState::Limited);
    assert!(report.cells[0].enumeration_complete);
    assert!(report.cells[0].sizes_twice_bp.is_empty());
    assert_eq!(
        report.tests[0].comparisons[0].distinction,
        Distinction::Unresolved
    );
}

#[test]
fn shared_budget_exhaustion_keeps_every_cell_and_withholds_all_missing_signal_claims() {
    let mut request = digest_request();
    request.tests = vec![pcr_test("first"), pcr_test("second")];
    request.limits.max_pcr_work = 1;
    let report = assess(&request).unwrap();
    assert_eq!(report.cells.len(), 4);
    assert_eq!(report.cells[0].state, CellState::Limited);
    assert_eq!(report.cells[1].state, CellState::Unsearched);
    assert_eq!(report.cells[2].state, CellState::Unsearched);
    assert_eq!(report.cells[3].state, CellState::Unsearched);
    assert_eq!(report.pcr_work, 1);
    assert!(
        report
            .tests
            .iter()
            .all(|test| test.comparisons[0].distinction == Distinction::Unresolved)
    );
}

#[test]
fn tolerance_requires_unmatched_size_and_never_infers_multiplicity_from_many_to_one_matches() {
    use super::assessment::distinguish;
    assert_eq!(
        distinguish(&[100, 104], &[102], 2),
        Distinction::WithinDeclaredTolerance
    );
    assert_eq!(distinguish(&[], &[], 0), Distinction::IdenticalSignature);
    assert_eq!(
        distinguish(&[100], &[], 10_000),
        Distinction::Separated(Witness {
            side: WitnessSide::Intended,
            size_twice_bp: 100
        })
    );
    // Exhaustive independent literal small-set oracle: distinction iff an element
    // of either set has no element in the other's closed tolerance interval.
    for amask in 0..16u32 {
        for bmask in 0..16u32 {
            let a: Vec<_> = (0..4)
                .filter(|i| amask & (1 << i) != 0)
                .map(|i| i * 2)
                .collect();
            let b: Vec<_> = (0..4)
                .filter(|i| bmask & (1 << i) != 0)
                .map(|i| i * 2)
                .collect();
            for tolerance in 0..=3 {
                let expected = a.iter().any(|x| {
                    b.iter()
                        .all(|y| (*x as isize - *y as isize).unsigned_abs() > tolerance)
                }) || b.iter().any(|x| {
                    a.iter()
                        .all(|y| (*x as isize - *y as isize).unsigned_abs() > tolerance)
                });
                assert_eq!(
                    matches!(distinguish(&a, &b, tolerance), Distinction::Separated(_)),
                    expected
                );
            }
        }
    }
}

#[test]
fn panel_search_matches_exhaustive_tiny_subset_oracle_and_discloses_limited_work() {
    use super::assessment::panel_search;
    // Four literal alternative tests covering three distinct decoys.
    let sets = [vec![1], vec![2], vec![1, 3], vec![2, 3]];
    let tests: Vec<_> = sets
        .iter()
        .enumerate()
        .map(|(test_index, set)| TestResult {
            test_index,
            definition: pcr_test(&format!("test{test_index}")),
            signatures: vec![],
            comparisons: (1..=3)
                .map(|decoy_subject_index| Comparison {
                    decoy_subject_index,
                    distinction: if set.contains(&decoy_subject_index) {
                        Distinction::Separated(Witness {
                            side: WitnessSide::Decoy,
                            size_twice_bp: 100,
                        })
                    } else {
                        Distinction::IdenticalSignature
                    },
                })
                .collect(),
            separated_decoys: set.len(),
            complete_comparisons: 3,
        })
        .collect();
    let panel = panel_search(&tests, 100_000);
    let mut oracle_best = 0;
    for mask in 1..16u32 {
        if mask.count_ones() > 3 {
            continue;
        }
        let covered: BTreeSet<_> = sets
            .iter()
            .enumerate()
            .filter(|(i, _)| mask & (1 << i) != 0)
            .flat_map(|(_, set)| set.iter().copied())
            .collect();
        oracle_best = oracle_best.max(covered.len());
    }
    assert_eq!(panel.separated_decoy_indices.len(), oracle_best);
    assert_eq!(panel.test_indices, [0, 3]);
    assert_eq!(panel.evaluated_combinations, 14); // C(4,1)+C(4,2)+C(4,3).
    assert!(panel.complete_through_size_three);
    let limited = panel_search(&tests, 1);
    assert!(!limited.complete_through_size_three);
    assert_eq!(limited.evaluated_combinations, 0);
    assert!(limited.test_indices.is_empty());
}

#[test]
fn exact_sources_request_policy_and_documentary_integrity_are_bound_without_mutation() {
    let request = digest_request();
    let original = serde_json::to_value(&request).unwrap();
    let report = assess(&request).unwrap();
    let context = revision(&request).unwrap();
    report.validate_binding(&context).unwrap();
    let mut stale = request.clone();
    stale.components[0].revision = "0".repeat(64);
    assert!(matches!(stale.validate(), Err(Error::Format(_))));
    let mut policy = request.clone();
    policy.readout.counterpart_tolerance_twice_bp = 2;
    assert!(
        report
            .validate_binding(&revision(&policy).unwrap())
            .is_err()
    );
    let mut altered = report.clone();
    altered.cells.swap(0, 1);
    assert!(altered.validate_binding(&context).is_err());
    let mut altered = report.clone();
    altered.tests[0].comparisons[0].distinction = Distinction::IdenticalSignature;
    assert!(altered.validate_binding(&context).is_err());
    // Current nullable fields must be present, even when no value is known.
    // Missing documentary fields must not silently become an older shape.
    for field in ["material_revision", "cleavage"] {
        let mut encoded = serde_json::to_value(&report).unwrap();
        encoded["cells"][0].as_object_mut().unwrap().remove(field);
        assert!(serde_json::from_value::<Report>(encoded).is_err());
    }
    let mut encoded = serde_json::to_value(&report).unwrap();
    encoded["components"][0]
        .as_object_mut()
        .unwrap()
        .remove("original_generation");
    assert!(serde_json::from_value::<Report>(encoded).is_err());
    assert_eq!(serde_json::to_value(&request).unwrap(), original);
    let mut oversized = request.clone();
    oversized.tests = (0..65).map(|i| pcr_test(&format!("test{i}"))).collect();
    assert!(matches!(oversized.validate(), Err(Error::TooLarge)));
}

#[test]
fn sorted_mixture_witnesses_preserve_closed_tolerance_at_integer_boundaries() {
    use super::assessment::distinguish;
    assert_eq!(
        distinguish(&[usize::MAX - 2], &[usize::MAX], 2),
        Distinction::WithinDeclaredTolerance
    );
    assert_eq!(
        distinguish(&[0], &[usize::MAX], usize::MAX),
        Distinction::WithinDeclaredTolerance
    );
    let intended: Vec<_> = (1..=1024).map(|n| n * 4).collect();
    let mut decoy = intended.clone();
    decoy.push(8192);
    assert_eq!(
        distinguish(&intended, &decoy, 1),
        Distinction::Separated(Witness {
            side: WitnessSide::Decoy,
            size_twice_bp: 8192,
        })
    );
}

// Recompute documentary identity after a deliberate forged body. These cases
// exercise shape consistency rather than the ordinary checksum mismatch guard.
fn reseal(report: &mut Report) {
    report.report_revision = crate::encoding::bounded_json::hash(&(
        &report.engine_version,
        &report.context_revision,
        &report.readout,
        &report.limits,
        &report.components,
        &report.subjects,
        &report.cells,
        &report.tests,
        &report.ranked_test_indices,
        &report.panel,
        report.cleavage_work,
        report.pcr_work,
        &report.statement,
    ))
    .unwrap();
}

#[test]
fn self_resealed_report_cannot_exceed_declared_work_or_contradict_pcr_lengths() {
    let mut request = digest_request();
    request.components[0] = component("intended", "ACGTTGCACCGGAATTCCGATGCA", Topology::Linear);
    request.tests = vec![pcr_test("PCR")];
    let report = assess(&request).unwrap();
    let context = revision(&request).unwrap();
    report.validate_binding(&context).unwrap();
    let mut forged = report.clone();
    forged.cells[0].pcr_geometries[0].product_length = 25;
    reseal(&mut forged);
    assert!(forged.validate_binding(&context).is_err());

    request.limits.max_pcr_work = 1;
    let mut report = assess(&request).unwrap();
    let context = revision(&request).unwrap();
    assert_eq!(report.pcr_work, 1);
    report.pcr_work = 2;
    reseal(&mut report);
    assert!(report.validate_binding(&context).is_err());

    let mut request = digest_request();
    request.limits.max_cleavage_work = 1;
    let mut report = assess(&request).unwrap();
    let context = revision(&request).unwrap();
    assert_eq!(report.cleavage_work, 1);
    report.cleavage_work = 2;
    reseal(&mut report);
    assert!(report.validate_binding(&context).is_err());
}
