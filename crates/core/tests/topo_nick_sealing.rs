//! Literal source strands, actual transfer composition and explicit ligase scope.
use collagene_core::{
    Molecule, Topology,
    protocol::{self, AnalysisJob, AnalysisResult},
    topo, topo_nick_sealing as adapter,
};
use topo::physical::{self as p, nick_sealing as n, projection};
fn fixture(directional: bool) -> (Vec<topo::Input>, adapter::Request) {
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
        adapter::Request {
            preparation,
            sealing,
        },
    )
}
fn count_calls(graph: &p::ProductGraph) -> usize {
    graph.chains.iter().map(|c| c.sequence.len()).sum()
}
#[test]
fn source_backed_ta_two_nicks_to_zero_keeps_every_original_call_and_activation() {
    let (inputs, r) = fixture(false);
    let raw = serde_json::to_vec(&inputs).unwrap();
    let expected = adapter::revision(&inputs, &r).unwrap();
    let report = adapter::inspect(&inputs, &r).unwrap();
    report.validate_binding(&expected).unwrap();
    assert_eq!(report.sealing.parent.product.nicks.len(), 2);
    assert_eq!(report.sealing.product.nicks.len(), 0);
    assert!(report.sealing.all_requested_nicks_sealed);
    assert_eq!(count_calls(&report.sealing.parent.product), 40);
    assert_eq!(count_calls(&report.sealing.product), 40);
    assert_eq!(
        report.sealing.product.released_activations,
        vec!["enzyme0", "enzyme1"]
    );
    assert_eq!(
        report.sealing.product.released_activations,
        report.sealing.parent.product.released_activations
    );
    assert!(
        report
            .sealing
            .product
            .chains
            .iter()
            .all(|c| matches!(c.ends, p::ChainEnds::CovalentlyCircular))
    );
    assert!(
        report
            .sealing
            .product
            .components
            .iter()
            .all(|c| c.covalently_closed)
    );
    assert_eq!(report.preparation.maps.len(), 4);
    let actual: std::collections::BTreeSet<_> = report
        .sealing
        .product
        .chains
        .iter()
        .map(|c| c.sequence.as_str())
        .collect();
    assert!(actual.contains("AGGGGCATCCCTTCGAATGA"));
    assert!(actual.contains("AGGGATGCCCCTTCATTCGA"));
    let mut retained = std::collections::BTreeSet::new();
    for c in &report.sealing.product.chains {
        for source in &c.sources {
            assert_eq!(source.start, 0);
            assert!(retained.insert(source.source_strand));
            assert_eq!(
                source.end,
                inputs[source.source_strand].molecule.sequence.len()
            );
        }
    }
    assert_eq!(retained, std::collections::BTreeSet::from([0, 1, 2, 3]));
    adapter::replay(&inputs, &r, &report).unwrap();
    assert_eq!(serde_json::to_vec(&inputs).unwrap(), raw);
}
#[test]
fn declared_five_oh_is_blocked_without_phosphate_invention() {
    let (inputs, mut r) = fixture(false);
    let topo::MaterialRecipe::Prepared { strands, .. } = &mut r.preparation.recipe.materials[0]
    else {
        panic!()
    };
    strands[0].five_prime = p::FivePrime::Hydroxyl;
    strands[1].five_prime = p::FivePrime::Hydroxyl;
    let report = adapter::inspect(&inputs, &r).unwrap();
    assert_eq!(report.sealing.parent.product.nicks.len(), 2);
    assert_eq!(report.sealing.product.nicks.len(), 2);
    assert!(!report.sealing.all_requested_nicks_sealed);
    assert!(
        report
            .sealing
            .steps
            .iter()
            .all(|s| s.outcome == p::Outcome::BlockedUnderDeclaredModel)
    );
    assert_eq!(report.sealing.product, report.sealing.parent.product);
}
#[test]
fn directional_cacc_sealing_never_deletes_or_repairs_displaced_flap() {
    let (inputs, mut r) = fixture(true);
    let parent = topo::simulate(&inputs, &r.preparation).unwrap();
    assert_eq!(parent.physical.product.unpaired.len(), 1);
    assert_eq!(parent.physical.product.nicks.len(), 1);
    r.sealing.nicks = parent
        .physical
        .product
        .nicks
        .iter()
        .map(|k| n::NickSelection {
            three: k.three,
            five: k.five,
        })
        .collect();
    let report = adapter::inspect(&inputs, &r).unwrap();
    let unpaired = &report.sealing.product.unpaired;
    assert_eq!(count_calls(&report.sealing.product), 52);
    assert!(!unpaired.is_empty());
    assert!(unpaired.iter().any(|u| u.terminal_three
        && &report.sealing.product.chains[u.chain].sequence.as_str()[u.start..u.end] == "GGTG"));
    assert!(
        report
            .sealing
            .product
            .components
            .iter()
            .any(|c| !c.covalently_closed)
    );
    assert_eq!(report.sealing.parent, parent.physical);
    assert_eq!(
        report.sealing.product.released_activations,
        vec!["enzyme0", "enzyme1"]
    );
}
#[test]
fn source_request_and_imported_evidence_never_restore_fresh_authority() {
    let (mut inputs, mut r) = fixture(false);
    let report = adapter::inspect(&inputs, &r).unwrap();
    inputs[0]
        .molecule
        .metadata
        .insert("unverified-supplier".into(), "activated".into());
    assert_ne!(
        adapter::revision(&inputs, &r).unwrap(),
        report.context_revision
    );
    assert!(adapter::replay(&inputs, &r, &report).is_err());
    let (inputs, _) = fixture(false);
    r.sealing.ligase_assumption.push_str(" changed");
    assert!(adapter::replay(&inputs, &r, &report).is_err());
    let mut changed = report.clone();
    changed
        .sealing
        .product
        .nicks
        .push(report.sealing.parent.product.nicks[0].clone());
    assert!(changed.validate_binding(&report.context_revision).is_err());
    let mut oversized = r.clone();
    oversized.sealing.nicks = vec![oversized.sealing.nicks[0].clone(); 65];
    assert!(adapter::preflight(&inputs, &oversized).is_err());
    let (_, mut capped) = fixture(false);
    capped.preparation.engine_limits.max_steps = 3;
    assert!(
        adapter::preflight(&inputs, &capped).is_err(),
        "two transfers plus two seals share the declared step budget"
    );
    capped.preparation.engine_limits.max_steps = 4;
    adapter::preflight(&inputs, &capped).unwrap();
}
#[test]
fn binary_worker_147_139_rejects_old_envelopes_and_stale_source_context() {
    let (inputs, request) = fixture(false);
    let expected = adapter::revision(&inputs, &request).unwrap();
    let encoded = protocol::encode_request(&protocol::Request {
        id: 1701,
        revision: expected.clone(),
        job: AnalysisJob::SealTopoNicks {
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
        145
    );
    let mut old = encoded.clone();
    old[4..6].copy_from_slice(&62u16.to_le_bytes());
    assert!(protocol::process_request(&old).is_err());
    let response = protocol::process_request(&encoded).unwrap();
    assert_eq!(
        postcard::take_from_bytes::<(u32, &str, u32)>(&response[42..])
            .unwrap()
            .0
            .2,
        137
    );
    let mut old = response.clone();
    old[4..6].copy_from_slice(&62u16.to_le_bytes());
    assert!(protocol::decode_response(&old).is_err());
    let decoded = protocol::decode_response(&response).unwrap();
    assert_eq!(decoded.revision, expected);
    let AnalysisResult::TopoNicksSealed(report) = decoded.result else {
        panic!("wrong result")
    };
    assert_eq!(report.sealing.product.nicks.len(), 0);
    let stale = protocol::encode_request(&protocol::Request {
        id: 1702,
        revision: "0".repeat(64),
        job: AnalysisJob::SealTopoNicks { inputs, request },
    })
    .unwrap();
    assert!(protocol::process_request(&stale).is_err());
}
