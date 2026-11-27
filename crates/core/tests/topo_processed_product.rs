//! Actual processed paired contours, conservative material bridge and digest oracle.
use collagene_core::{
    Feature, Location, Molecule, Topology,
    analysis::Enzyme,
    chemistry, modified_duplex as tracked,
    protocol::{self, AnalysisJob, AnalysisResult},
    topo, topo_nick_sealing as treatment, topo_processed_product as adapter,
};
use topo::physical::{self as p, nick_sealing as n, processed_projection as pp, projection};
fn fixture(directional: bool) -> (Vec<topo::Input>, treatment::Request) {
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
        transfers:vec![p::Transfer{donor:0,acceptor:2,geometry:if directional{p::Geometry::DirectionalCacc}else{p::Geometry::ThreePrimeTA}},p::Transfer{donor:1,acceptor:3,geometry:if directional{p::Geometry::Flush}else{p::Geometry::ThreePrimeTA}}],
        engine_limits:p::Limits::default(),projection_options:projection::Options::default()
    };
    let sealing = n::Request {
        ligase_assumption:
            "Explicit ligase hypothesis: only templated 3′OH/5′P nick sealing, no repair enzymes"
                .into(),
        nicks: vec![
            n::NickSelection {
                three: p::Call {
                    strand: 2,
                    index: 6,
                },
                five: p::Call {
                    strand: 0,
                    index: 0,
                },
            },
            n::NickSelection {
                three: p::Call {
                    strand: 3,
                    index: 6,
                },
                five: p::Call {
                    strand: 1,
                    index: 0,
                },
            },
        ],
    };
    (
        inputs,
        treatment::Request {
            preparation,
            sealing,
        },
    )
}
fn selected(inputs: &[topo::Input], treatment: treatment::Request) -> adapter::Request {
    let native = treatment::inspect(inputs, &treatment).unwrap();
    adapter::Request {
        treatment,
        selection: pp::Selection {
            sealing_context_revision: native.sealing.context_revision,
            component: 0,
            reference_seed: p::Call {
                strand: 0,
                index: 0,
            },
        },
        options: projection::Options {
            require_all_requested_transfers_applied: true,
            require_known_terminal_chemistry: true,
            ..Default::default()
        },
    }
}
fn feature() -> Feature {
    Feature {
        kind: "misc_feature".into(),
        location: Location::range(2, 5),
        qualifiers: vec![("label".into(), Some("Original vector region".into()))],
        needs_review: false,
    }
}
#[test]
fn fresh_ta_bridge_preserves_calls_annotations_and_actual_two_strand_taqi_digest() {
    let (mut inputs, treatment) = fixture(false);
    inputs[0].molecule.features.push(feature());
    inputs[0].molecule.metadata.insert(
        "methylation".into(),
        "supplier claims fully methylated".into(),
    );
    let request = selected(&inputs, treatment);
    let expected = adapter::revision(&inputs, &request).unwrap();
    let report = adapter::inspect(&inputs, &request).unwrap();
    report.validate_binding(&expected).unwrap();
    let d = report.duplex.as_ref().unwrap();
    assert_eq!(d.topology, Topology::Circular);
    assert!(d.nicks.is_empty());
    assert_eq!(d.top.as_str(), "AGGGGCATCCCTTCGAATGA");
    assert_eq!(d.bottom.as_str(), "TCATTCGAAGGGATGCCCCT");
    let product = report.product.as_ref().unwrap();
    assert_eq!(product.sequence, d.top);
    assert!(chemistry::has_unmapped_physical_evidence(product));
    assert!(chemistry::Duplex::from_molecule(product).is_err());
    assert!(
        product
            .features
            .iter()
            .any(|f| f.location == Location::range(2, 5))
    );
    assert!(report.annotations.is_some());
    let t = report.tracked.as_ref().unwrap();
    assert!(t.generation.is_empty());
    assert_eq!(t.top_modifications.n6a_default, tracked::CallState::Unknown);
    assert_eq!(
        t.bottom_modifications.c5c_default,
        tracked::CallState::Unknown
    );
    let enzyme = Enzyme {
        name: "TaqI geometry oracle".into(),
        recognition: "TCGA".into(),
        top_cut: 1,
        bottom_cut: 3,
    };
    let fragments = chemistry::digest(d, std::slice::from_ref(&enzyme)).unwrap();
    assert_eq!(fragments.len(), 1);
    let f = &fragments[0];
    assert_eq!(f.top_start, 13);
    assert_eq!(f.bottom_start, 15);
    assert_eq!(f.duplex.top.as_str(), "CGAATGAAGGGGCATCCCTT");
    assert_eq!(f.duplex.bottom.as_str(), "CGAAGGGATGCCCCTTCATT");
    assert_eq!(f.duplex.bottom_start, 2);
    let (left, right) = f.duplex.ends().unwrap();
    assert_eq!(left.overhang, "CG");
    assert_eq!(right.overhang, "CG");
    assert_eq!(left.chemistry, chemistry::Terminus::CLEAVED);
    let actual_tracked = tracked::digest(t, &[enzyme]).unwrap();
    assert_eq!(actual_tracked[0].duplex.duplex, f.duplex);
    assert_eq!(
        actual_tracked[0].duplex.top_modifications.n6a_default,
        tracked::CallState::Unknown
    );
    assert_eq!(
        report.projection.sealing.product.released_activations,
        vec!["enzyme0", "enzyme1"]
    );
    adapter::replay(&inputs, &request, &report).unwrap();
}
#[test]
fn remaining_ligatable_nick_is_mapped_but_five_oh_nick_is_not_collapsed() {
    let (inputs, mut treatment) = fixture(false);
    treatment.sealing.nicks.truncate(1);
    let request = selected(&inputs, treatment.clone());
    let report = adapter::inspect(&inputs, &request).unwrap();
    assert_eq!(report.duplex.as_ref().unwrap().nicks.len(), 1);
    assert_eq!(report.projection.sealing.product.nicks.len(), 1);
    let topo::MaterialRecipe::Prepared { strands, .. } =
        &mut treatment.preparation.recipe.materials[0]
    else {
        panic!()
    };
    strands[1].five_prime = p::FivePrime::Hydroxyl;
    let mut request = selected(&inputs, treatment);
    request.options.require_known_terminal_chemistry = false;
    let report = adapter::inspect(&inputs, &request).unwrap();
    assert!(matches!(
        report.projection.outcome,
        projection::Outcome::Proven(_)
    ));
    assert!(report.duplex.is_none());
    assert!(report.tracked.is_none());
    assert!(report.product.is_none());
    assert!(report.bridge_notes.iter().any(|s| s.contains("3′OH/5′P")));
}
#[test]
fn directional_flap_withholds_material_and_keeps_all_other_calls() {
    let (inputs, mut treatment) = fixture(true);
    let original = topo::simulate(&inputs, &treatment.preparation).unwrap();
    treatment.sealing.nicks = original
        .physical
        .product
        .nicks
        .iter()
        .map(|k| n::NickSelection {
            three: k.three,
            five: k.five,
        })
        .collect();
    let request = selected(&inputs, treatment);
    let report = adapter::inspect(&inputs, &request).unwrap();
    assert!(matches!(
        report.projection.outcome,
        projection::Outcome::Withheld(projection::Withheld::UnpairedNucleotides { count: 4 })
    ));
    assert!(report.duplex.is_none());
    assert!(report.tracked.is_none());
    assert!(report.product.is_none());
    assert!(report.annotations.is_none());
    assert_eq!(
        report
            .projection
            .sealing
            .product
            .chains
            .iter()
            .map(|c| c.sequence.len())
            .sum::<usize>(),
        52
    );
    assert!(report.projection.sealing.product.unpaired.iter().any(|u| {
        &report.projection.sealing.product.chains[u.chain]
            .sequence
            .as_str()[u.start..u.end]
            == "GGTG"
    }));
}
fn linear(five: p::FivePrime, three: p::ThreePrime) -> (Vec<topo::Input>, adapter::Request) {
    let inputs = vec![topo::Input {
        role: "original".into(),
        molecule: Molecule::dna("Source", "ACGTAC", Topology::Linear).unwrap(),
    }];
    let strands = (0..2)
        .map(|i| topo::PreparedStrand {
            id: format!("s{i}"),
            segments: vec![topo::SourceInterval {
                source_role: "original".into(),
                start: 0,
                length: 6,
                reverse: i == 1,
            }],
            five_prime: five,
            three_prime: three.clone(),
        })
        .collect();
    let treatment = treatment::Request {
        preparation: topo::Request {
            recipe: topo::GenerationRecipe {
                schema_version: 1,
                materials: vec![topo::MaterialRecipe::Prepared {
                    declaration_id: "Declared linear pair".into(),
                    assumption: "Explicit terminal states, no reference-derived chemistry".into(),
                    strands,
                    pairing: vec![p::PairSpan {
                        strand_a: 0,
                        start_a: 0,
                        strand_b: 1,
                        start_b: 0,
                        length: 6,
                    }],
                }],
            },
            transfers: vec![],
            engine_limits: p::Limits::default(),
            projection_options: projection::Options::default(),
        },
        sealing: n::Request {
            ligase_assumption: "No nick requested in this complete linear material".into(),
            nicks: vec![],
        },
    };
    let mut request = selected(&inputs, treatment);
    request.options = projection::Options::default();
    (inputs, request)
}
#[test]
fn unknown_linear_ends_remain_unknown_and_unsupported_terminal_species_are_withheld() {
    let (inputs, request) = linear(p::FivePrime::Unknown, p::ThreePrime::Unknown);
    let report = adapter::inspect(&inputs, &request).unwrap();
    let d = report.duplex.as_ref().unwrap();
    assert_eq!(d.left, chemistry::Terminus::UNKNOWN);
    assert_eq!(d.right, chemistry::Terminus::UNKNOWN);
    for (five, three, word) in [
        (p::FivePrime::Blocked, p::ThreePrime::Hydroxyl, "5′"),
        (
            p::FivePrime::Hydroxyl,
            p::ThreePrime::Phosphate,
            "phosphate",
        ),
    ] {
        let (inputs, request) = linear(five, three);
        let report = adapter::inspect(&inputs, &request).unwrap();
        assert!(matches!(
            report.projection.outcome,
            projection::Outcome::Proven(_)
        ));
        assert!(report.duplex.is_none());
        assert!(report.product.is_none());
        assert!(report.bridge_notes.iter().any(|s| s.contains(word)));
    }
    // A retained actual activation cannot be erased into a generic blocked end.
    let (inputs, mut request) = linear(p::FivePrime::Hydroxyl, p::ThreePrime::Hydroxyl);
    let topo::MaterialRecipe::Prepared { strands, .. } =
        &mut request.treatment.preparation.recipe.materials[0]
    else {
        panic!()
    };
    strands[0].three_prime = p::ThreePrime::Phosphotyrosyl {
        enzyme: p::Enzyme::VacciniaTypeIBTyr274,
        activation_id: "retained-activation".into(),
    };
    request = selected(&inputs, request.treatment);
    request.options = projection::Options::default();
    let report = adapter::inspect(&inputs, &request).unwrap();
    assert!(report.duplex.is_none());
    assert!(
        report
            .bridge_notes
            .iter()
            .any(|s| s.contains("phosphotyrosyl"))
    );
}
#[test]
fn source_selection_limits_and_imported_material_claims_need_fresh_reconstruction() {
    let (mut inputs, treatment) = fixture(false);
    let request = selected(&inputs, treatment);
    let report = adapter::inspect(&inputs, &request).unwrap();
    let mut tampered = report.clone();
    tampered
        .tracked
        .as_mut()
        .unwrap()
        .top_modifications
        .n6a_default = tracked::CallState::Unmodified;
    assert!(tampered.validate_binding(&report.context_revision).is_err());
    let mut tampered = report.clone();
    tampered
        .product
        .as_mut()
        .unwrap()
        .metadata
        .remove("duplex_model");
    assert!(tampered.validate_binding(&report.context_revision).is_err());
    let mut stale = request.clone();
    stale.selection.reference_seed.index = 1;
    assert_ne!(
        adapter::revision(&inputs, &stale).unwrap(),
        report.context_revision
    );
    assert!(adapter::replay(&inputs, &stale, &report).is_err());
    inputs[0]
        .molecule
        .metadata
        .insert("changed source".into(), "true".into());
    assert!(adapter::replay(&inputs, &request, &report).is_err());
    let fresh = adapter::inspect(&inputs, &request).unwrap();
    assert_ne!(
        fresh.context_revision, report.context_revision,
        "unchanged physical geometry may be reconstructed afresh, but metadata edits change source identity"
    );
    let mut changed_treatment = request.clone();
    changed_treatment
        .treatment
        .sealing
        .ligase_assumption
        .push_str(" changed condition declaration");
    assert!(
        adapter::inspect(&inputs, &changed_treatment).is_err(),
        "prior physical selection cannot bind a changed treatment context"
    );
    let mut invalid = request.clone();
    invalid.options.limits.max_pairs = 0;
    assert!(adapter::preflight(&inputs, &invalid).is_err());
    let (inputs, treatment) = fixture(false);
    let mut capped = selected(&inputs, treatment);
    capped.options.limits.max_pairs = 1;
    let bounded = adapter::inspect(&inputs, &capped).unwrap();
    assert!(bounded.duplex.is_none());
    assert!(matches!(
        bounded.projection.outcome,
        projection::Outcome::Withheld(_)
    ));
}
#[test]
fn binary_worker_148_140_rejects_pre64_and_stale_source_context() {
    let (inputs, treatment) = fixture(false);
    let request = selected(&inputs, treatment);
    let expected = adapter::revision(&inputs, &request).unwrap();
    let encoded = protocol::encode_request(&protocol::Request {
        id: 1801,
        revision: expected.clone(),
        job: AnalysisJob::ProjectProcessedTopo {
            inputs: inputs.clone(),
            request: request.clone(),
        },
    })
    .unwrap();
    assert_eq!(
        postcard::take_from_bytes::<(u32, &str, u32)>(&encoded[42..])
            .unwrap()
            .0
            .2,
        146
    );
    let mut old = encoded.clone();
    old[4..6].copy_from_slice(&63u16.to_le_bytes());
    assert!(protocol::process_request(&old).is_err());
    let response = protocol::process_request(&encoded).unwrap();
    assert_eq!(
        postcard::take_from_bytes::<(u32, &str, u32)>(&response[42..])
            .unwrap()
            .0
            .2,
        138
    );
    let mut old = response.clone();
    old[4..6].copy_from_slice(&63u16.to_le_bytes());
    assert!(protocol::decode_response(&old).is_err());
    let decoded = protocol::decode_response(&response).unwrap();
    assert_eq!(decoded.revision, expected);
    let AnalysisResult::ProcessedTopoProduct(report) = decoded.result else {
        panic!("wrong result")
    };
    assert!(report.tracked.is_some());
    let stale = protocol::encode_request(&protocol::Request {
        id: 1802,
        revision: "0".repeat(64),
        job: AnalysisJob::ProjectProcessedTopo { inputs, request },
    })
    .unwrap();
    assert!(protocol::process_request(&stale).is_err());
}
