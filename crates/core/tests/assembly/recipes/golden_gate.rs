//! Hand-derived Type IIS strand diagrams and a separate permutation/rotation oracle.
use collagene_core::{chemistry::Duplex, cloning::golden_gate::*};
use collagene_model::{Error, Molecule, Sequence, Topology};
use std::collections::BTreeSet;

fn donor(left: &str, payload: &str, right: &str) -> Molecule {
    Molecule::dna(
        format!("{left}-{payload}-{right}"),
        &format!("GGTCTCA{left}{payload}{right}AGAGACC"),
        Topology::Linear,
    )
    .unwrap()
}
fn generated_unmodified(
    inputs: &[Molecule],
    options: &GoldenGateOptions,
) -> Result<GoldenGateReport, Error> {
    use collagene_core::chemistry::{
        methylation::MethylationState,
        tracked_duplex::{CleavageBudget, GenerationAssumptions, TrackedDuplex},
    };
    let assumptions = GenerationAssumptions {
        dam: MethylationState::Unmodified,
        dcm: MethylationState::Unmodified,
        cpg: MethylationState::Unmodified,
    };
    let actual = inputs
        .iter()
        .map(|m| TrackedDuplex::from_input(m, &assumptions))
        .collect::<Result<Vec<_>, Error>>()?;
    let mut budget = CleavageBudget::new(100_000_000)?;
    Ok(simulate_tracked(inputs, &actual, options, None, &mut budget)?.report)
}
fn options() -> GoldenGateOptions {
    GoldenGateOptions {
        max_fragments: 32,
        max_products: 1000,
        ..Default::default()
    }
}
fn stable_sequences(report: &GoldenGateReport) -> BTreeSet<String> {
    report
        .candidates
        .iter()
        .filter(|c| c.restriction_site_free)
        .map(|c| canonical(c.duplex.top.as_str()))
        .collect()
}
// Enumerate every rotation and both reference orientations, intentionally not the implementation's
// linear-time rotation algorithm or hash. Small fixtures only.
fn canonical(s: &str) -> String {
    let rc = Sequence::dna(s).unwrap().reverse_complement().unwrap();
    [s, rc.as_str()]
        .into_iter()
        .flat_map(|x| (0..x.len()).map(move |i| format!("{}{}", &x[i..], &x[..i])))
        .min()
        .unwrap()
}

#[test]
fn hand_derived_nonpalindromic_overhangs_reconstruct_unique_physical_circle() {
    let inputs = [donor("AATG", "CCCC", "GCTT"), donor("GCTT", "TTTT", "AATG")];
    let report = generated_unmodified(&inputs, &options()).unwrap();
    assert!(report.enumeration_complete);
    assert!(!report.exhaustive_mixture_search);
    assert_eq!(report.fragments.len(), 6);
    assert_eq!(report.candidates.len(), 2); // Two directed references for one physical circle.
    assert_eq!(report.stable_physical_product_count, 1);
    let center = &report.fragments[1];
    assert_eq!(
        (
            center.source_index,
            center.fragment_index,
            center.top_start,
            center.bottom_start
        ),
        (0, 1, 7, 11)
    );
    assert_eq!(center.duplex.top.as_str(), "AATGCCCC");
    assert_eq!(center.duplex.bottom.as_str(), "AAGCGGGG");
    assert_eq!(center.duplex.bottom_start, 4);
    assert_eq!(center.left_end.as_ref().unwrap().overhang, "AATG");
    assert_eq!(center.right_end.as_ref().unwrap().overhang, "AAGC");
    let selected = report
        .candidates
        .iter()
        .find(|c| c.duplex.top.as_str() == "AATGCCCCGCTTTTTT")
        .unwrap();
    assert_eq!(
        selected.path,
        [
            OrientedFragment {
                fragment: 1,
                reverse: false
            },
            OrientedFragment {
                fragment: 4,
                reverse: false
            }
        ]
    );
    assert_eq!(
        selected.duplex.bottom,
        selected.duplex.top.reverse_complement().unwrap()
    );
    assert!(selected.duplex.nicks.is_empty());
    assert!(selected.recut_sites.is_empty());
    assert_ne!(
        report.candidates[0].directed_sequence_key,
        report.candidates[1].directed_sequence_key
    );
    assert_eq!(
        report.candidates[0].physical_sequence_key,
        report.candidates[1].physical_sequence_key
    );
    assert_eq!(
        stable_sequences(&report),
        BTreeSet::from([canonical("AATGCCCCGCTTTTTT")])
    );
    for (index, source) in report.sources.iter().enumerate() {
        assert_eq!(source.revision, inputs[index].revision());
        assert_eq!(source.recognition_sites.len(), 2);
    }
}

#[test]
fn reversing_and_rotating_sources_preserve_exact_assembly_geometry() {
    let left = donor("AATG", "CCCC", "GCTT");
    let right = donor("GCTT", "TTTT", "AATG");
    let expected = BTreeSet::from([canonical("AATGCCCCGCTTTTTT")]);
    let reversed = Molecule::dna(
        "reversed donor",
        left.sequence.reverse_complement().unwrap().as_str(),
        Topology::Linear,
    )
    .unwrap();
    let r = generated_unmodified(&[reversed, right.clone()], &options()).unwrap();
    assert_eq!(stable_sequences(&r), expected);
    assert!(r.candidates.iter().any(|c|c.restriction_site_free&&c.path.iter().any(|p|p.fragment==1&&p.reverse)));
    let mut right_circle = right;
    right_circle.topology = Topology::Circular;
    // Both donor circles also produce recuttable carrier/background candidates; these are retained.
    for origin in 0..left.sequence.len() {
        let s = left.sequence.as_str();
        let rotated = Molecule::dna(
            "rotated circular donor",
            &format!("{}{}", &s[origin..], &s[..origin]),
            Topology::Circular,
        )
        .unwrap();
        let report = generated_unmodified(&[rotated, right_circle.clone()], &options()).unwrap();
        assert!(
            report.enumeration_complete,
            "origin {origin}: {:?}",
            report.truncation_reasons
        );
        assert_eq!(stable_sequences(&report), expected, "origin {origin}");
        assert!(report.candidates.iter().any(|c| !c.restriction_site_free));
    }
}

#[test]
fn newly_created_recognition_crossing_junction_is_not_a_stable_gate_product() {
    // Neither retained top fragment contains GGTCTC: A ends G, B starts GTCT+C.
    let inputs = [
        donor("AATG", "CCCCG", "GTCT"),
        donor("GTCT", "CAAAA", "AATG"),
    ];
    let report = generated_unmodified(&inputs, &options()).unwrap();
    assert!(report.enumeration_complete);
    assert_eq!(report.candidates.len(), 2);
    assert_eq!(report.stable_physical_product_count, 0);
    assert!(
        report
            .sources
            .iter()
            .all(|s| s.recognition_sites.len() == 2)
    );
    let candidate = report
        .candidates
        .iter()
        .find(|c| c.duplex.top.as_str() == "AATGCCCCGGTCTCAAAA")
        .unwrap();
    assert!(!candidate.restriction_site_free);
    assert_eq!(candidate.recut_sites.len(), 1);
    let site = &candidate.recut_sites[0];
    assert_eq!(
        (
            site.recognition_start,
            site.reverse,
            site.top_cut,
            site.bottom_cut
        ),
        (8, false, 15, 1)
    );
    assert!(site.fully_double_stranded && site.cleavable_under_model);
}

#[test]
fn palindromic_overhangs_expose_self_ligation_and_distinct_histories() {
    let one = generated_unmodified(&[donor("AATT", "CCCC", "AATT")], &options()).unwrap();
    assert_eq!(one.candidates.len(), 2);
    assert!(one.candidates.iter().all(|c| c.path.len() == 1));
    assert_eq!(one.stable_physical_product_count, 1);
    let inputs = [
        donor("AATT", "AAAA", "AATT"),
        donor("AATT", "CCCC", "AATT"),
        donor("AATT", "GGGG", "AATT"),
    ];
    let report = generated_unmodified(&inputs, &options()).unwrap();
    assert!(report.enumeration_complete);
    // Every nonempty subset, every ordering, both orientations. Rotations identify cycles only.
    // n=1: 3*2, n=2: 3*4, n=3: 2*8 = 34 directed source histories.
    assert_eq!(report.candidates.len(), 34);
    let actual = report
        .candidates
        .iter()
        .map(|c| {
            c.path
                .iter()
                .map(|p| (p.fragment, p.reverse))
                .collect::<Vec<_>>()
        })
        .collect::<BTreeSet<_>>();
    fn all_paths(
        path: &mut Vec<(usize, bool)>,
        used: &mut Vec<usize>,
        out: &mut BTreeSet<Vec<(usize, bool)>>,
    ) {
        if !path.is_empty() {
            let cycle = (0..path.len())
                .map(|i| {
                    path[i..]
                        .iter()
                        .chain(&path[..i])
                        .copied()
                        .collect::<Vec<_>>()
                })
                .min()
                .unwrap();
            out.insert(cycle);
        }
        for source in 0..3 {
            if !used.contains(&source) {
                used.push(source);
                for reverse in [false, true] {
                    path.push((source * 3 + 1, reverse));
                    all_paths(path, used, out);
                    path.pop();
                }
                used.pop();
            }
        }
    }
    let mut expected = BTreeSet::new();
    all_paths(&mut vec![], &mut vec![], &mut expected);
    assert_eq!(actual, expected);
    // Same source-sequence bytes can have distinct original-source paths; these are not erased.
    let duplicate = generated_unmodified(
        &[donor("AATT", "CCCC", "AATT"), donor("AATT", "CCCC", "AATT")],
        &options(),
    )
    .unwrap();
    assert_eq!(
        duplicate
            .candidates
            .iter()
            .filter(|c| c.path.len() == 1)
            .count(),
        4
    );
}

#[test]
fn work_depth_and_output_caps_never_claim_complete_absence_of_alternatives() {
    let inputs = [
        donor("AATT", "AAAA", "AATT"),
        donor("AATT", "CCCC", "AATT"),
        donor("AATT", "GGGG", "AATT"),
    ];
    let complete = generated_unmodified(&inputs, &options()).unwrap();
    assert!(complete.enumeration_complete);
    for limited in [
        GoldenGateOptions {
            max_states: 1,
            ..options()
        },
        GoldenGateOptions {
            max_products: 1,
            ..options()
        },
        GoldenGateOptions {
            max_fragments: 1,
            ..options()
        },
    ] {
        let r = generated_unmodified(&inputs, &limited).unwrap();
        assert!(!r.enumeration_complete);
        assert!(!r.truncation_reasons.is_empty());
        assert!(!r.exhaustive_mixture_search);
    }
    let fixed = generated_unmodified(
        &inputs,
        &GoldenGateOptions {
            both_orientations: false,
            ..options()
        },
    )
    .unwrap();
    assert_eq!(fixed.candidates.len(), 8);
    assert!(fixed.enumeration_complete);
    assert!(
        fixed
            .warnings
            .iter()
            .any(|w| w.contains("Reverse fragment"))
    );
}

#[test]
fn wrong_overhangs_uncut_carryover_and_unknown_physical_states_are_explicit() {
    let inputs = [donor("AATG", "CCCC", "GCTT"), donor("CCAA", "TTTT", "AATG")];
    let report = generated_unmodified(&inputs, &options()).unwrap();
    assert!(report.enumeration_complete);
    assert!(report.candidates.is_empty());
    let intact = Molecule::dna(
        "uncut circular carryover",
        "AAAACCCCTTTT",
        Topology::Circular,
    )
    .unwrap();
    let report = generated_unmodified(std::slice::from_ref(&intact), &options()).unwrap();
    assert_eq!(report.carryover_sources, [0]);
    assert!(report.candidates.is_empty());
    assert!(report.enumeration_complete);
    let mut physical = intact;
    physical
        .metadata
        .insert("duplex_model".into(), "untrusted assertion".into());
    assert!(matches!(
        generated_unmodified(&[physical], &options()),
        Err(Error::Unsupported(_))
    ));
    let unknown = Molecule::dna("ambiguous", "GGTCTCNNNNNNNNGAGACC", Topology::Linear).unwrap();
    assert!(matches!(
        generated_unmodified(&[unknown], &options()),
        Err(Error::Unsupported(_))
    ));
    assert!(generated_unmodified(&[], &options()).is_err());
    assert!(
        generated_unmodified(
            &inputs,
            &GoldenGateOptions {
                max_fragments: 0,
                ..options()
            }
        )
        .is_err()
    );
}

#[test]
fn factual_enzyme_profiles_and_linear_edge_recognition_are_distinct_from_cleavability() {
    for (kind, motif, top, bottom) in [
        (GoldenGateEnzyme::BsaI, "GGTCTC", 7, 11),
        (GoldenGateEnzyme::BsmBI, "CGTCTC", 7, 11),
        (GoldenGateEnzyme::BbsI, "GAAGAC", 8, 12),
        (GoldenGateEnzyme::SapI, "GCTCTTC", 8, 11),
    ] {
        let p = kind.profile();
        assert_eq!(
            (p.recognition.as_str(), p.top_cut, p.bottom_cut),
            (motif, top, bottom)
        );
    }
    let linear = Duplex::from_molecule(
        &Molecule::dna("terminal motif", "GGTCTC", Topology::Linear).unwrap(),
    )
    .unwrap();
    let actual =
        collagene_core::chemistry::tracked_duplex::TrackedDuplex::annealed_unmodified(&linear)
            .unwrap();
    let mut budget =
        collagene_core::chemistry::tracked_duplex::CleavageBudget::new(100_000_000).unwrap();
    let (sites, _) = assess_product(&actual, GoldenGateEnzyme::BsaI, &mut budget).unwrap();
    assert_eq!(sites.len(), 1);
    assert!(sites[0].fully_double_stranded);
    assert!(!sites[0].cleavable_under_model);
    assert_eq!((sites[0].top_cut, sites[0].bottom_cut), (7, 11));
}

#[test]
fn graph_reconstructed_pcr_substrates_are_distinct_from_untrusted_imported_chemistry() {
    let originals = [donor("AATG", "CCCC", "GCTT"), donor("GCTT", "TTTT", "AATG")];
    let physical = originals
        .iter()
        .map(|m| {
            collagene_core::chemistry::tracked_duplex::TrackedDuplex::exact_pcr(&m.sequence)
                .unwrap()
        })
        .collect::<Vec<_>>();
    let snapshots = physical
        .iter()
        .enumerate()
        .map(|(i, d)| d.molecule(format!("PCR {i}")).unwrap())
        .collect::<Vec<_>>();
    assert!(simulate(&snapshots, &options()).is_err());
    let mut budget =
        collagene_core::chemistry::tracked_duplex::CleavageBudget::new(100_000_000).unwrap();
    let report = simulate_tracked(&snapshots, &physical, &options(), None, &mut budget)
        .unwrap()
        .report;
    assert_eq!(
        stable_sequences(&report),
        BTreeSet::from([canonical("AATGCCCCGCTTTTTT")])
    );
    let mut wrong = physical.clone();
    wrong.swap(0, 1);
    assert!(simulate_tracked(&snapshots, &wrong, &options(), None, &mut budget).is_err());
    assert!(simulate_tracked(&snapshots, &physical[..1], &options(), None, &mut budget).is_err());
}

#[test]
fn each_supported_type_iis_profile_reconstructs_its_hand_specified_overhang_geometry() {
    for (enzyme, motif, spacer, left, right, expected) in [
        (
            GoldenGateEnzyme::BsmBI,
            "CGTCTC",
            "A",
            "AATG",
            "GCTT",
            "AATGCCCCGCTTTTTT",
        ),
        (
            GoldenGateEnzyme::BbsI,
            "GAAGAC",
            "AA",
            "AATG",
            "GCTT",
            "AATGCCCCGCTTTTTT",
        ),
        (
            GoldenGateEnzyme::SapI,
            "GCTCTTC",
            "A",
            "AAT",
            "GCT",
            "AATCCCCGCTTTTT",
        ),
    ] {
        let reverse_motif = Sequence::dna(motif).unwrap().reverse_complement().unwrap();
        let make = |a: &str, payload: &str, b: &str| {
            Molecule::dna(
                "alternate profile",
                &format!(
                    "{motif}{spacer}{a}{payload}{b}{spacer}{}",
                    reverse_motif.as_str()
                ),
                Topology::Linear,
            )
            .unwrap()
        };
        let inputs = [make(left, "CCCC", right), make(right, "TTTT", left)];
        let report = generated_unmodified(
            &inputs,
            &GoldenGateOptions {
                enzyme,
                ..options()
            },
        )
        .unwrap();
        assert!(report.enumeration_complete);
        assert_eq!(report.stable_physical_product_count, 1);
        assert_eq!(
            stable_sequences(&report),
            BTreeSet::from([canonical(expected)])
        );
        assert_eq!(report.fragments[1].duplex.bottom_start, left.len() as i64);
        assert_eq!(report.fragments[1].top_start, motif.len() + spacer.len());
        assert_eq!(
            report.fragments[1].bottom_start,
            (motif.len() + spacer.len() + left.len()) as i64
        );
    }
}

#[test]
fn linear_time_sequence_keys_match_independent_full_rotation_equivalence_classes() {
    use std::collections::BTreeMap;
    let inputs = [
        donor("AATT", "AAAA", "AATT"),
        donor("AATT", "CCCC", "AATT"),
        donor("AATT", "AAAA", "AATT"),
        donor("AATT", "GGGG", "AATT"),
    ];
    let report = generated_unmodified(&inputs, &options()).unwrap();
    assert!(report.enumeration_complete);
    assert_eq!(report.candidates.len(), 192);
    let (mut by_sequence, mut by_key) = (
        BTreeMap::<String, String>::new(),
        BTreeMap::<String, String>::new(),
    );
    let mut directed = BTreeMap::<String, String>::new();
    for c in &report.candidates {
        let canonical = canonical(c.duplex.top.as_str());
        if let Some(old) = by_sequence.insert(canonical.clone(), c.physical_sequence_key.clone()) {
            assert_eq!(old, c.physical_sequence_key);
        }
        if let Some(old) = by_key.insert(c.physical_sequence_key.clone(), canonical.clone()) {
            assert_eq!(old, canonical);
        }
        let s = c.duplex.top.as_str();
        let top = (0..s.len())
            .map(|i| format!("{}{}", &s[i..], &s[..i]))
            .min()
            .unwrap();
        if let Some(old) = directed.insert(top, c.directed_sequence_key.clone()) {
            assert_eq!(old, c.directed_sequence_key);
        }
    }
    assert_eq!(report.stable_physical_product_count, by_sequence.len());
}

#[test]
fn thirty_two_linear_donors_are_not_accidentally_limited_to_sixty_four_digest_pieces() {
    let bases = [b'A', b'C', b'G', b'T'];
    let mut tags = vec![];
    for i in 0..256usize {
        let tag = String::from_utf8(vec![
            bases[(i >> 6) & 3],
            bases[(i >> 4) & 3],
            bases[(i >> 2) & 3],
            bases[i & 3],
        ])
        .unwrap();
        let rc = Sequence::dna(&tag).unwrap().reverse_complement().unwrap();
        if tag.as_str() < rc.as_str() {
            tags.push(tag);
            if tags.len() == 32 {
                break;
            }
        }
    }
    let inputs = (0..32)
        .map(|i| donor(&tags[i], "AAAA", &tags[(i + 1) % 32]))
        .collect::<Vec<_>>();
    let report = generated_unmodified(&inputs, &options()).unwrap();
    assert_eq!(report.fragments.len(), 96);
    assert!(
        report.enumeration_complete,
        "{:?}",
        report.truncation_reasons
    );
    assert_eq!(report.candidates.len(), 2);
    assert_eq!(report.stable_physical_product_count, 1);
    assert!(report.candidates.iter().all(|c| c.path.len() == 32));
    let expected = (0..32)
        .map(|i| format!("{}AAAA", tags[i]))
        .collect::<String>();
    assert_eq!(
        stable_sequences(&report),
        BTreeSet::from([canonical(&expected)])
    );
}

#[test]
fn full_u128_pool_boundary_is_addressable_and_next_fragment_rejects_before_shift() {
    let make = |count| {
        Molecule::dna(
            "pool boundary",
            &"GGTCTCAAAAACCCC".repeat(count),
            Topology::Circular,
        )
        .unwrap()
    };
    let opts = GoldenGateOptions {
        both_orientations: false,
        max_fragments: 1,
        ..options()
    };
    let report = generated_unmodified(&[make(128)], &opts).unwrap();
    assert_eq!(report.fragments.len(), 128);
    assert_eq!(report.candidates.len(), 128);
    assert!(report.candidates.iter().any(|c| c.path[0].fragment == 127));
    assert!(matches!(
        generated_unmodified(&[make(129)], &opts),
        Err(Error::Unsupported(_))
    ));
}
