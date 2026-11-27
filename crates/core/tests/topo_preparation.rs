//! Actual source calls and independently stated directional strand diagrams.
use collagene_core::{Molecule, Topology, analysis::Primer, topo, topo_preparation as prep};
use topo::physical::{self as graph, projection};
fn input(role: &str, dna: &str) -> topo::Input {
    topo::Input {
        role: role.into(),
        molecule: Molecule::dna(role, dna, Topology::Linear).unwrap(),
    }
}
fn strand(role: &str, length: usize) -> topo::PreparedStrand {
    topo::PreparedStrand {
        id: role.into(),
        segments: vec![topo::SourceInterval {
            source_role: role.into(),
            start: 0,
            length,
            reverse: false,
        }],
        five_prime: graph::FivePrime::Phosphate,
        three_prime: graph::ThreePrime::Phosphotyrosyl {
            enzyme: graph::Enzyme::VacciniaTypeIBTyr274,
            activation_id: format!("{role}-declared-activation"),
        },
    }
}
fn fixture() -> (Vec<topo::Input>, prep::PreparationRequest) {
    (
        vec![
            input("vt", "AAGGGGCATCCCTT"),
            input("vb", "GGTGAAGGGATGCCCCTT"),
            input("template", "CACCATGACG"),
        ],
        prep::PreparationRequest {
            mode: prep::PreparationMode::DirectionalCacc,
            source: prep::PreparationSource::Original {
                vector: prep::PreparedVectorDeclaration {
                    declaration_id: "declared-directional-vector".into(),
                    assumption:
                        "Explicit synthetic mechanism fixture, not stock activation evidence".into(),
                    strands: vec![strand("vt", 14), strand("vb", 18)],
                    pairing: vec![graph::PairSpan {
                        strand_a: 0,
                        start_a: 0,
                        strand_b: 1,
                        start_b: 4,
                        length: 14,
                    }],
                },
                insert: prep::ExactPcrDeclaration {
                    material_id: "actual-blunt-pcr".into(),
                    template_role: "template".into(),
                    forward: Primer::new("Forward", "CACCA", "").unwrap(),
                    reverse: Primer::new("Reverse", "CGTCA", "").unwrap(),
                    tail: None,
                },
            },
            engine_limits: graph::Limits::default(),
            projection_options: projection::Options::default(),
        },
    )
}
fn inspect(
    inputs: &[topo::Input],
    request: &prep::PreparationRequest,
) -> topo::OrientationInspection {
    let compiled = prep::compile_preparation(inputs, request).unwrap();
    topo::inspect_orientations(
        inputs,
        &compiled.request,
        compiled.vector_strands,
        compiled.insert_strands,
        compiled.directional,
    )
    .unwrap()
}
#[test]
fn literal_14_18_vector_and_actual_cacc_pcr_keep_the_flap_and_withhold_contour() {
    let (inputs, request) = fixture();
    let before = serde_json::to_value(&inputs).unwrap();
    let compiled = prep::compile_preparation(&inputs, &request).unwrap();
    assert!(compiled.directional);
    assert!(compiled.request.transfers.is_empty());
    let report = inspect(&inputs, &request);
    assert_eq!(report.reference.candidates.len(), 2);
    assert!(!report.reference.biochemical_mixture_exhaustive);
    assert_eq!(
        report.preparation.substrates.strands[2].sequence.as_str(),
        "CACCATGACG"
    );
    assert_eq!(
        report.preparation.substrates.strands[3].sequence.as_str(),
        "CGTCATGGTG"
    );
    let actual = &report.reference.candidates[0];
    assert!(actual.all_requested_transfers_applied);
    assert_eq!(actual.source_bases, 52);
    assert_eq!(
        actual
            .product
            .chains
            .iter()
            .map(|chain| chain.sequence.len())
            .sum::<usize>(),
        52
    );
    assert!(
        actual
            .product
            .chains
            .iter()
            .any(|chain| chain.sequence.as_str() == "GGTGAAGGGATGCCCCTTCGTCATGGTG")
    );
    assert!(
        actual
            .steps
            .iter()
            .any(|step| step
                .displacement
                .as_ref()
                .is_some_and(|displacement| matches!(
                    displacement.kind,
                    graph::DisplacementKind::DuplexAnchoredThreePrimeFlap
                )))
    );
    assert!(!actual.product.nicks.is_empty());
    let selected = &report.candidate_requests[0];
    let selection = topo::Selection {
        batch_revision: topo::batch_revision(&inputs, selected).unwrap(),
        component: 0,
        reference_seed: graph::Call {
            strand: 0,
            index: 0,
        },
    };
    let product = topo::product(&inputs, selected, &selection).unwrap();
    assert!(product.reference.is_none());
    assert!(matches!(
        product.projection.outcome,
        projection::Outcome::Withheld(projection::Withheld::UnpairedNucleotides { count: 4 })
    ));
    let target = Molecule::dna(
        "Separate intended reference",
        "AAGGGGCATCCCTTCACCATGACG",
        Topology::Circular,
    )
    .unwrap();
    let policy = topo::AcceptancePolicy {
        allow_nicks: true,
        allow_unknown_terminal_chemistry: true,
        allow_unapplied_transfers: true,
        allow_other_components: true,
        allow_circular_rotation: true,
        allow_reverse_complement: true,
    };
    assert!(matches!(
        topo::check_target(&inputs, selected, &selection, &target, &policy)
            .unwrap()
            .outcome,
        topo::TargetOutcome::Withheld { .. }
    ));
    assert_eq!(serde_json::to_value(&inputs).unwrap(), before);
}
#[test]
fn literal_cacc_primer_tail_is_distinct_from_cacc_already_in_the_template() {
    let (mut inputs, mut request) = fixture();
    inputs[2] = input("template", "CATGACGTAC");
    {
        let prep::PreparationSource::Original { insert, .. } = &mut request.source else {
            panic!("original declaration fixture")
        };
        insert.forward = Primer::new("User forward", "CATGA", "CACC").unwrap();
    }
    {
        let prep::PreparationSource::Original { insert, .. } = &mut request.source else {
            panic!("original declaration fixture")
        };
        insert.reverse = Primer::new("User reverse", "GTACG", "").unwrap();
    }
    let tailed = inspect(&inputs, &request);
    assert_eq!(
        tailed.preparation.substrates.strands[2].sequence.as_str(),
        "CACCCATGACGTAC"
    );
    assert!(tailed.preparation.maps.iter().any(|map| matches!(
        map.origin,
        topo::Origin::ExactPcrCopy {
            forward_tail_length: 4,
            ..
        }
    )));
    assert!(
        !tailed
            .preparation
            .maps
            .iter()
            .any(|map| matches!(map.origin, topo::Origin::SynthesizedTail { .. }))
    );
    inputs[2] = input("template", "CACCCATGACGTAC");
    {
        let prep::PreparationSource::Original { insert, .. } = &mut request.source else {
            panic!("original declaration fixture")
        };
        insert.forward = Primer::new("Already in template", "CACCC", "").unwrap();
    }
    let encoded = inspect(&inputs, &request);
    assert_eq!(
        encoded.preparation.substrates.strands[2].sequence,
        tailed.preparation.substrates.strands[2].sequence
    );
    assert!(encoded.preparation.maps.iter().any(|map| matches!(
        map.origin,
        topo::Origin::ExactPcrCopy {
            forward_tail_length: 0,
            ..
        }
    )));
    assert!(
        matches!(&request.source, prep::PreparationSource::Original { insert, .. } if insert.tail.is_none())
    );
}
#[test]
fn reversed_drawing_missing_activation_and_missing_cacc_never_gain_transfer_authority() {
    let (mut inputs, mut request) = fixture();
    inputs[1] = input("vb", "GTGGAAGGGATGCCCCTT");
    let wrong = inspect(&inputs, &request);
    assert!(!wrong.reference.candidates[0].all_requested_transfers_applied);
    inputs[1] = input("vb", "GGTGAAGGGATGCCCCTT");
    {
        let prep::PreparationSource::Original { vector, .. } = &mut request.source else {
            panic!("original declaration fixture")
        };
        vector.strands[0].three_prime = graph::ThreePrime::Hydroxyl;
    }
    let inactive = inspect(&inputs, &request);
    assert!(!inactive.reference.candidates[0].all_requested_transfers_applied);
    let (mut inputs, mut request) = fixture();
    inputs[2] = input("template", "AACCATGACG");
    {
        let prep::PreparationSource::Original { insert, .. } = &mut request.source else {
            panic!("original declaration fixture")
        };
        insert.forward = Primer::new("No CACC", "AACCA", "").unwrap();
    }
    let missing = inspect(&inputs, &request);
    assert!(!missing.reference.candidates[0].all_requested_transfers_applied);
    {
        let prep::PreparationSource::Original { insert, .. } = &mut request.source else {
            panic!("original declaration fixture")
        };
        insert.tail = Some(collagene_core::tailing::TailingBase::A);
    }
    assert!(prep::compile(&request).is_err());
    {
        let prep::PreparationSource::Original { insert, .. } = &mut request.source else {
            panic!("original declaration fixture")
        };
        insert.tail = None;
    }
    {
        let prep::PreparationSource::Original { vector, .. } = &mut request.source else {
            panic!("original declaration fixture")
        };
        vector.pairing[0].start_b = 0;
    }
    assert!(prep::compile(&request).is_err());
}
#[test]
fn native_declaration_bounds_and_ambiguous_actual_pcr_are_explicit() {
    let (inputs, mut request) = fixture();
    {
        let prep::PreparationSource::Original { insert, .. } = &mut request.source else {
            panic!("original declaration fixture")
        };
        insert.forward.tail = "A".repeat(196);
    }
    assert!(prep::compile(&request).is_err()); //201 ordered bases
    {
        let prep::PreparationSource::Original { insert, .. } = &mut request.source else {
            panic!("original declaration fixture")
        };
        insert.forward.tail.clear();
    }
    let mut absent = inputs.clone();
    absent.remove(1);
    assert!(prep::compile_preparation(&absent, &request).is_err());
    {
        let prep::PreparationSource::Original { vector, .. } = &mut request.source else {
            panic!("original declaration fixture")
        };
        vector.strands[1].three_prime = vector.strands[0].three_prime.clone();
    }
    assert!(prep::compile(&request).is_err());
    let (mut inputs, mut request) = fixture();
    inputs[2] = input("template", "CACCATGACGCACCATGACG");
    {
        let prep::PreparationSource::Original { insert, .. } = &mut request.source else {
            panic!("original declaration fixture")
        };
        insert.forward = Primer::new("Repeated forward", "CACCA", "").unwrap();
    }
    {
        let prep::PreparationSource::Original { insert, .. } = &mut request.source else {
            panic!("original declaration fixture")
        };
        insert.reverse = Primer::new("Repeated reverse", "CGTCA", "").unwrap();
    }
    let compiled = prep::compile_preparation(&inputs, &request).unwrap();
    assert!(
        topo::inspect_orientations(
            &inputs,
            &compiled.request,
            compiled.vector_strands,
            compiled.insert_strands,
            true
        )
        .is_err()
    );
}
