//! Literal source strands, actual transfer composition and explicit ligase scope.
use collagene_core::cloning::topo;
use collagene_core::cloning::topo::physical::{self as p, processing as process, projection};
use collagene_model::{Molecule, Topology};
fn fixture(directional: bool) -> (Vec<topo::Input>, topo::Request) {
    let dna = if directional {
        [
            "AAGGGGCATCCCTT",
            "GGTGAAGGGATGCCCCTT",
            "CACCATGACG",
            "CGTCATGGTG",
        ]
    } else {
        ["AGGGGCATCCCTT", "AGGGATGCCCCTT", "CGAATGA", "CATTCGA"]
    };
    let inputs = dna
        .iter()
        .enumerate()
        .map(|(i, dna)| topo::Input {
            role: format!("original-{i}"),
            molecule: Molecule::dna(format!("source{i}"), dna, Topology::Linear).unwrap(),
        })
        .collect();
    let strands = dna
        .iter()
        .enumerate()
        .map(|(i, dna)| topo::PreparedStrand {
            id: format!("strand{i}"),
            segments: vec![topo::SourceInterval {
                source_role: format!("original-{i}"),
                start: 0,
                length: dna.len(),
                reverse: false,
            }],
            five_prime: if i < 2 {
                p::FivePrime::Phosphate
            } else {
                p::FivePrime::Hydroxyl
            },
            three_prime: if i < 2 {
                p::ThreePrime::Phosphotyrosyl {
                    enzyme: p::Enzyme::VacciniaTypeIBTyr274,
                    activation_id: format!("enzyme{i}"),
                }
            } else {
                p::ThreePrime::Hydroxyl
            },
        })
        .collect();
    let pairing = if directional {
        vec![
            p::PairSpan {
                strand_a: 0,
                start_a: 0,
                strand_b: 1,
                start_b: 4,
                length: 14,
            },
            p::PairSpan {
                strand_a: 2,
                start_a: 0,
                strand_b: 3,
                start_b: 0,
                length: 10,
            },
        ]
    } else {
        vec![
            p::PairSpan {
                strand_a: 0,
                start_a: 0,
                strand_b: 1,
                start_b: 0,
                length: 12,
            },
            p::PairSpan {
                strand_a: 2,
                start_a: 0,
                strand_b: 3,
                start_b: 0,
                length: 6,
            },
        ]
    };
    let preparation=topo::Request{
        recipe:topo::GenerationRecipe{schema_version:1,materials:vec![topo::MaterialRecipe::Prepared{declaration_id:"declared original material".into(),assumption:"Activation and supplied strand geometry are explicit hypotheses, never inferred from reference DNA".into(),strands,pairing}]},
        steps:vec![p::Step::Transfer(p::Transfer{donor: p::Call {strand:0,index:if directional { 13 } else { 12 }},acceptor:p::Call{strand:2,index:0},geometry:if directional{p::Geometry::DirectionalCacc}else{p::Geometry::ThreePrimeTA}}),p::Step::Transfer(p::Transfer{donor: p::Call {strand:1,index:if directional { 17 } else { 12 }},acceptor:p::Call{strand:3,index:0},geometry:if directional{p::Geometry::Flush}else{p::Geometry::ThreePrimeTA}})],
        engine_limits:p::Limits::default(),projection_options:projection::Options::default()
    };
    let mut preparation = preparation;
    for strand in [2, 3] {
        preparation
            .steps
            .push(p::Step::Process(process::Action::SealNick {
                selection: p::nick_sealing::NickSelection {
                    three: p::Call { strand, index: 6 },
                    five: p::Call {
                        strand: strand - 2,
                        index: 0,
                    },
                },
                ligase_assumption:
                    "Explicit ligase hypothesis: only templated3OH/5P nick sealing, no repair"
                        .into(),
            }));
    }
    (inputs, preparation)
}

fn count_calls(graph: &p::ProductGraph) -> usize {
    graph.chains.iter().map(|c| c.sequence.len()).sum()
}
#[test]
fn source_backed_ta_two_nicks_to_zero_keeps_every_original_call_and_activation() {
    let (inputs, request) = fixture(false);
    let mut transfers = request.clone();
    transfers.steps.truncate(2);
    let before = topo::simulate(&inputs, &transfers).unwrap();
    let after = topo::simulate(&inputs, &request).unwrap();
    assert_eq!(before.physical.product.nicks.len(), 2);
    assert_eq!(after.physical.product.nicks.len(), 0);
    assert_eq!(count_calls(&before.physical.product), 40);
    assert_eq!(count_calls(&after.physical.product), 40);
    assert_eq!(
        before.physical.product.released_activations,
        after.physical.product.released_activations
    );
    assert_eq!(after.physical.product.released_activations.len(), 2);
    assert!(after.physical.all_requested_steps_applied);
    assert_eq!(after.physical.steps.len(), 4);
    for step in after.physical.steps.iter().skip(2) {
        assert_eq!(step.outcome(), p::Outcome::Applied);
    }
    assert!(
        after
            .physical
            .product
            .chains
            .iter()
            .all(|c| matches!(c.ends, p::ChainEnds::CovalentlyCircular))
    );
}
#[test]
fn declared_five_oh_blocks_sealing_without_phosphate_invention() {
    let (mut inputs, mut request) = fixture(false);
    let topo::MaterialRecipe::Prepared { strands, .. } = &mut request.recipe.materials[0] else {
        panic!()
    };
    strands[0].five_prime = p::FivePrime::Hydroxyl;
    strands[1].five_prime = p::FivePrime::Hydroxyl;
    let result = topo::simulate(&inputs, &request).unwrap();
    assert_eq!(result.physical.product.nicks.len(), 2);
    assert!(!result.physical.all_requested_steps_applied);
    assert!(
        result
            .physical
            .steps
            .iter()
            .skip(2)
            .all(|s| s.outcome() == p::Outcome::BlockedUnderDeclaredModel)
    );
    inputs[0].molecule.name.push_str(" changed");
    assert_ne!(
        result.batch_revision,
        topo::batch_revision(&inputs, &request).unwrap()
    );
}
#[test]
fn directional_sealing_never_deletes_or_repairs_displaced_flap() {
    let (inputs, mut request) = fixture(true);
    request.steps.truncate(2);
    let parent = topo::simulate(&inputs, &request).unwrap();
    assert_eq!(parent.physical.product.nicks.len(), 1);
    for nick in &parent.physical.product.nicks {
        request
            .steps
            .push(p::Step::Process(process::Action::SealNick {
                selection: p::nick_sealing::NickSelection {
                    three: nick.three,
                    five: nick.five,
                },
                ligase_assumption: "Explicit ideal ligase without flap hydrolysis".into(),
            }));
    }
    let result = topo::simulate(&inputs, &request).unwrap();
    assert_eq!(count_calls(&result.physical.product), 52);
    assert!(result.physical.product.unpaired.iter().any(|u| {
        &result.physical.product.chains[u.chain].sequence.as_str()[u.start..u.end] == "GGTG"
    }));
    assert_eq!(result.physical.product.released_activations.len(), 2);
}
#[test]
fn chronological_step_budget_and_exact_transcript_context_are_required() {
    let (inputs, mut request) = fixture(false);
    let context = topo::batch_revision(&inputs, &request).unwrap();
    let p::Step::Process(process::Action::SealNick {
        ligase_assumption, ..
    }) = &mut request.steps[2]
    else {
        panic!()
    };
    ligase_assumption.push_str(" explicitly changed");
    assert_ne!(topo::batch_revision(&inputs, &request).unwrap(), context);
    request.engine_limits.max_steps = 3;
    assert!(topo::source_preflight(&inputs, &request).is_err());
    request.engine_limits.max_steps = 4;
    topo::source_preflight(&inputs, &request).unwrap();
}
