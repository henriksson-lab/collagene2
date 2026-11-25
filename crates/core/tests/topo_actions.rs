//! Source-bound ideal processing is proven from literal original calls. Released
//! flap material remains in the native graph and never becomes silent deletion.
use collagene_core::{
    Feature, Location, Molecule, Topology, chemistry, modified_duplex, topo,
    topo_actions as adapter,
};
use topo::physical::{self as p, processing as actions, processing_projection as pp, projection};
fn call(strand: usize, index: usize) -> p::Call {
    p::Call { strand, index }
}
fn sources() -> (Vec<topo::Input>, topo::Request) {
    let dna = [
        "AAGGGGCATCCCTT",
        "GGTGAAGGGATGCCCCTT",
        "CACCATGACG",
        "CGTCATGGTG",
    ];
    let mut inputs = dna
        .iter()
        .enumerate()
        .map(|(i, dna)| topo::Input {
            role: format!("original{i}"),
            molecule: Molecule::dna(format!("source{i}"), dna, Topology::Linear).unwrap(),
        })
        .collect::<Vec<_>>();
    inputs[2].molecule.features.push(Feature {
        kind: "misc_feature".into(),
        location: Location::range(4, 10),
        qualifiers: vec![(
            "label".into(),
            Some("insert retained coding-side annotation".into()),
        )],
        needs_review: false,
    });
    let strands = dna
        .iter()
        .enumerate()
        .map(|(i, dna)| topo::PreparedStrand {
            id: format!("strand{i}"),
            segments: vec![topo::SourceInterval {
                source_role: format!("original{i}"),
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
                    activation_id: format!("activation{i}"),
                }
            } else {
                p::ThreePrime::Hydroxyl
            },
        })
        .collect();
    let request = topo::Request {
        recipe: topo::GenerationRecipe {
            schema_version: 1,
            materials: vec![topo::MaterialRecipe::Prepared {
                declaration_id: "explicit declared geometry".into(),
                assumption: "Original activation is an operator hypothesis".into(),
                strands,
                pairing: vec![
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
                ],
            }],
        },
        transfers: vec![
            p::Transfer {
                donor: 0,
                acceptor: 2,
                geometry: p::Geometry::DirectionalCacc,
            },
            p::Transfer {
                donor: 1,
                acceptor: 3,
                geometry: p::Geometry::Flush,
            },
        ],
        engine_limits: p::Limits::default(),
        projection_options: projection::Options::default(),
    };
    (inputs, request)
}
fn processing(cut: bool, seals: bool) -> actions::Request {
    let mut operations = vec![];
    if cut {
        operations.push(actions::Action::CleaveThreePrimeFlap {
            bond: actions::Bond {
                upstream: call(3, 5),
                downstream: call(3, 6),
            },
            displacement_step: 0,
            profile: actions::FlapCleavageProfile::IdealDuplexAnchoredThreePrimeHydrolysisV1,
            assumption: "Ideal explicitly requested junction hydrolysis, not named enzyme efficacy"
                .into(),
        });
    }
    if seals {
        for (three, five) in [(call(2, 9), call(0, 0)), (call(3, 5), call(1, 0))] {
            operations.push(actions::Action::SealNick {
                selection: p::nick_sealing::NickSelection { three, five },
                ligase_assumption:
                    "Explicit ideal ligase hypothesis for only exposed templated3OH5P".into(),
            });
        }
    }
    actions::Request {
        actions: operations,
        limits: actions::ProcessingLimits::default(),
    }
}
fn selected(
    inputs: &[topo::Input],
    preparation: topo::Request,
    processing: actions::Request,
) -> adapter::Request {
    let inspection = adapter::inspect_actions(
        inputs,
        &adapter::ActionsRequest {
            preparation: preparation.clone(),
            processing: processing.clone(),
        },
    )
    .unwrap();
    let component = inspection
        .processing
        .product
        .components
        .iter()
        .position(|c| {
            c.chains.iter().any(|&i| {
                inspection.processing.product.chains[i]
                    .sources
                    .iter()
                    .any(|s| s.source_strand == 0 && s.start == 0)
            })
        })
        .unwrap();
    adapter::Request {
        preparation,
        processing,
        selection: pp::Selection {
            processing_context_revision: inspection.processing.context_revision,
            component,
            reference_seed: call(0, 0),
        },
        options: projection::Options::default(),
    }
}
#[test]
fn literal_directional_cut_and_sealing_prove_48_main_calls_plus_free_four() {
    let (inputs, preparation) = sources();
    let request = selected(&inputs, preparation, processing(true, true));
    let out = adapter::inspect(&inputs, &request).unwrap();
    out.validate_binding(&adapter::revision(&inputs, &request).unwrap())
        .unwrap();
    let product = out.product.as_ref().unwrap();
    assert_eq!(product.sequence.as_str(), "AAGGGGCATCCCTTCACCATGACG");
    assert_eq!(product.topology, Topology::Circular);
    assert!(
        product
            .features
            .iter()
            .any(|f| f.label() == "insert retained coding-side annotation")
    );
    assert!(chemistry::has_unmapped_physical_evidence(product));
    assert!(chemistry::Duplex::from_molecule(product).is_err());
    let graph = &out.projection.processing.product;
    assert_eq!(
        graph.chains.iter().map(|c| c.sequence.len()).sum::<usize>(),
        52
    );
    assert_eq!(graph.components.len(), 2);
    assert!(
        graph.components.iter().any(|c| c.paired_nucleotides == 48
            && c.unpaired_nucleotides == 0
            && c.covalently_closed)
    );
    assert!(
        graph
            .components
            .iter()
            .any(|c| c.paired_nucleotides == 0 && c.unpaired_nucleotides == 4)
    );
    let free = graph
        .chains
        .iter()
        .find(|c| c.sequence.as_str() == "GGTG")
        .unwrap();
    assert_eq!(free.sources[0].source_strand, 3);
    assert_eq!((free.sources[0].start, free.sources[0].end), (6, 10));
    assert_eq!(
        free.ends,
        p::ChainEnds::Linear {
            five_prime: p::FivePrime::Phosphate,
            three_prime: p::ThreePrime::Hydroxyl
        }
    );
    let tracked = out.tracked.as_ref().unwrap();
    assert!(tracked.duplex.nicks.is_empty());
    assert!(tracked.generation.is_empty());
    assert_eq!(
        tracked.top_modifications.n6a_default,
        modified_duplex::CallState::Unknown
    );
    assert_eq!(
        tracked.bottom_modifications.c5c_default,
        modified_duplex::CallState::Unknown
    );
    adapter::replay(&inputs, &request, &out).unwrap();
}
#[test]
fn cleavage_without_ligation_retains_two_real_chemically_explicit_nicks() {
    let (inputs, preparation) = sources();
    let r = selected(&inputs, preparation, processing(true, false));
    let out = adapter::inspect(&inputs, &r).unwrap();
    assert_eq!(out.duplex.as_ref().unwrap().nicks.len(), 2);
    assert!(
        out.projection
            .processing
            .product
            .chains
            .iter()
            .any(|c| c.sequence.as_str() == "GGTG")
    );
    assert!(
        out.projection
            .processing
            .product
            .nicks
            .iter()
            .any(|n| n.three == call(3, 5) && n.three_prime == p::ThreePrime::Hydroxyl)
    );
}
#[test]
fn no_requested_cut_never_erases_the_flap_or_invents_an_ordinary_product() {
    let (inputs, preparation) = sources();
    let r = selected(&inputs, preparation, processing(false, false));
    let out = adapter::inspect(&inputs, &r).unwrap();
    assert!(matches!(
        out.projection.outcome,
        projection::Outcome::Withheld(projection::Withheld::UnpairedNucleotides { count: 4 })
    ));
    assert!(out.product.is_none());
    assert!(out.duplex.is_none());
    assert_eq!(out.projection.processing.source_bases, 52);
}
#[test]
fn actual_restriction_digest_consumes_the_bridged_duplex_and_keeps_modifications_unknown() {
    let (inputs, preparation) = sources();
    let r = selected(&inputs, preparation, processing(true, true));
    let out = adapter::inspect(&inputs, &r).unwrap();
    let enzyme = collagene_core::analysis::Enzyme {
        name: "literal-CATG".into(),
        recognition: "CATG".into(),
        top_cut: 1,
        bottom_cut: 3,
    };
    let fragments = modified_duplex::digest(out.tracked.as_ref().unwrap(), &[enzyme]).unwrap();
    assert!(!fragments.is_empty());
    assert_eq!(
        fragments
            .iter()
            .map(|f| f.duplex.duplex.top.len())
            .sum::<usize>(),
        24
    );
    assert!(
        fragments
            .iter()
            .all(|f| f.duplex.top_modifications.n6a_default == modified_duplex::CallState::Unknown)
    );
}
#[test]
fn source_context_action_order_and_imported_graph_tampering_never_authorize_replay() {
    let (inputs, preparation) = sources();
    let r = selected(&inputs, preparation, processing(true, true));
    let out = adapter::inspect(&inputs, &r).unwrap();
    let mut wrong = r.clone();
    wrong.processing.actions.swap(0, 1);
    assert!(adapter::inspect(&inputs, &wrong).is_err());
    let mut foreign = inputs.clone();
    foreign[0].molecule.name = "different exact source".into();
    assert!(adapter::replay(&foreign, &r, &out).is_err());
    let mut forged = out.clone();
    forged
        .projection
        .processing
        .product
        .chains
        .retain(|c| c.sequence.as_str() != "GGTG");
    assert!(forged.validate_binding(&out.context_revision).is_err());
    let raw = serde_json::to_vec(&out).unwrap();
    let restored: adapter::Report = serde_json::from_slice(&raw).unwrap();
    adapter::replay(&inputs, &r, &restored).unwrap();
}
