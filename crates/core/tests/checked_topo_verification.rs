use collagene_core::{Molecule, Topology, topo as t, verification::*};
use t::physical::{self as p, projection as v};
fn molecule(dna: &str) -> Molecule {
    Molecule::dna("same display name", dna, Topology::Linear).unwrap()
}
fn strand(role: &str, n: usize, activated: bool) -> t::PreparedStrand {
    t::PreparedStrand {
        id: role.into(),
        segments: vec![t::SourceInterval {
            source_role: role.into(),
            start: 0,
            length: n,
            reverse: false,
        }],
        five_prime: if activated {
            p::FivePrime::Phosphate
        } else {
            p::FivePrime::Hydroxyl
        },
        three_prime: if activated {
            p::ThreePrime::Phosphotyrosyl {
                enzyme: p::Enzyme::VacciniaTypeIBTyr274,
                activation_id: role.into(),
            }
        } else {
            p::ThreePrime::Hydroxyl
        },
    }
}
fn pair(a: usize, x: usize, b: usize, y: usize, n: usize) -> p::PairSpan {
    p::PairSpan {
        strand_a: a,
        start_a: x,
        strand_b: b,
        start_b: y,
        length: n,
    }
}
fn plan(directional: bool) -> (DesignPlan, Molecule) {
    let sequences = if directional {
        vec![
            ("vt", "AAGGGGCATCCCTT"),
            ("vb", "GGTGAAGGGATGCCCCTT"),
            ("it", "CACCATGACG"),
            ("ib", "CGTCATGGTG"),
        ]
    } else {
        vec![
            ("vt", "AGGGGCATCCCTT"),
            ("vb", "AGGGATGCCCCTT"),
            ("it", "CGAATGA"),
            ("ib", "CATTCGA"),
        ]
    };
    let inputs = sequences
        .iter()
        .map(|(id, dna)| {
            let molecule = molecule(dna);
            PlanInput {
                id: (*id).into(),
                revision: molecule.revision(),
                molecule,
            }
        })
        .collect::<Vec<_>>();
    let sources = sequences
        .iter()
        .map(|(role, _)| TopoSource {
            role: (*role).into(),
            source: Source::Input((*role).into()),
        })
        .collect::<Vec<_>>();
    let request=t::Request{recipe:t::GenerationRecipe{schema_version:1,materials:vec![t::MaterialRecipe::Prepared{
        declaration_id:"explicit TA strands".into(),assumption:"Declared source-corresponding strands and activation; no supplier state is inferred.".into(),
        strands:sequences.iter().enumerate().map(|(i,(role,dna))|strand(role,dna.len(),i<2)).collect(),
        pairing:if directional {vec![pair(0,0,1,4,14),pair(2,0,3,0,10)]}else{vec![pair(0,0,1,0,12),pair(2,0,3,0,6)]}}]},
        steps:vec![p::Step::Transfer(p::Transfer{donor: p::Call {strand:0,index:if directional { 13 } else { 12 }},acceptor:p::Call{strand:2,index:0},geometry:if directional {p::Geometry::DirectionalCacc}else{p::Geometry::ThreePrimeTA}}),p::Step::Transfer(p::Transfer{donor: p::Call {strand:1,index:if directional { 17 } else { 12 }},acceptor:p::Call{strand:3,index:0},geometry:if directional {p::Geometry::Flush}else{p::Geometry::ThreePrimeTA}})],engine_limits:p::Limits::default(),projection_options:v::Options::default()};
    let actual = inputs
        .iter()
        .map(|input| &input.molecule)
        .collect::<Vec<_>>();
    let bound = topo_inputs(&sources, &actual, &request).unwrap();
    let selection = t::Selection {
        batch_revision: t::batch_revision(&bound, &request).unwrap(),
        component: 0,
        reference_seed: p::Call {
            strand: 0,
            index: 0,
        },
    };
    let target = Molecule::dna(
        "independent target",
        if directional {
            "AAGGGGCATCCCTTCACCATGACG"
        } else {
            "AGGGGCATCCCTTCGAATGA"
        },
        Topology::Circular,
    )
    .unwrap();
    (
        DesignPlan {
            name: "TOPO plan".into(),
            target_revision: target.revision(),
            inputs,
            steps: vec![Operation::CheckedTopo {
                sources,
                request,
                selection,
                policy: t::AcceptancePolicy {
                    allow_nicks: true,
                    ..Default::default()
                },
            }],
        },
        target,
    )
}

#[test]
fn explicit_nick_policy_accepts_graph_without_creating_repaired_duplex_authority() {
    let (mut plan, target) = plan(false);
    let accepted = verify(&plan, &target);
    assert_eq!(
        accepted.status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        accepted.findings
    );
    assert!(
        accepted.step_reports[0]
            .topo
            .as_ref()
            .unwrap()
            .projection
            .physical
            .fully_specified_dna
    );
    assert!(
        accepted.step_reports[0]
            .topo
            .as_ref()
            .unwrap()
            .annotations
            .is_some()
    );
    assert!(collagene_core::chemistry::has_unmapped_physical_evidence(
        &accepted.products[0]
    ));
    assert_eq!(
        accepted.products[0].metadata[collagene_core::recombination::REFERENCE_ONLY_KEY],
        "true"
    );
    if let Operation::CheckedTopo { policy, .. } = &mut plan.steps[0] {
        policy.allow_nicks = false;
    }
    let rejected = verify(&plan, &target);
    assert_eq!(rejected.status, Status::Invalid, "{:?}", rejected.findings);
    assert!(rejected.products.is_empty());
    assert_eq!(rejected.step_reports.len(), 1);
    assert!(
        rejected.step_reports[0]
            .topo
            .as_ref()
            .unwrap()
            .annotations
            .is_some()
    );
    assert!(
        rejected
            .findings
            .iter()
            .any(|f| f.contains("retains nicks"))
    );
}

#[test]
fn explicit_rotation_and_reverse_policy_survives_final_generic_target_comparison() {
    let (base, original) = plan(false);
    for reverse in [false, true] {
        let dna = if reverse {
            original
                .sequence
                .reverse_complement()
                .unwrap()
                .as_str()
                .to_owned()
        } else {
            let dna = original.sequence.as_str();
            format!("{}{}", &dna[3..], &dna[..3])
        };
        let target = Molecule::dna(
            "independent alternative coordinate orientation",
            &dna,
            Topology::Circular,
        )
        .unwrap();
        let mut candidate = base.clone();
        candidate.target_revision = target.revision();
        assert!(!plan_target_matches(&candidate, &original, &target));
        let rejected = verify(&candidate, &target);
        assert_eq!(rejected.status, Status::Invalid, "{:?}", rejected.findings);
        assert_eq!(rejected.step_reports.len(), 1);
        if let Operation::CheckedTopo { policy, .. } = &mut candidate.steps[0] {
            if reverse {
                policy.allow_reverse_complement = true;
            } else {
                policy.allow_circular_rotation = true;
            }
        }
        assert!(plan_target_matches(&candidate, &original, &target));
        let accepted = verify(&candidate, &target);
        assert_eq!(
            accepted.status,
            Status::ConsistentUnderAssumptions,
            "{:?}",
            accepted.findings
        );
    }
}

#[test]
fn unresolved_contour_and_intermediate_checked_topo_are_incomplete() {
    let (mut candidate, target) = plan(true);
    let withheld = verify(&candidate, &target);
    assert_eq!(
        withheld.status,
        Status::Incomplete,
        "{:?}",
        withheld.findings
    );
    assert_eq!(withheld.step_reports.len(), 1);
    assert!(withheld.products.is_empty());
    candidate.steps.push(Operation::TreatEnds {
        template: Source::Step(0),
        treatment: collagene_core::chemistry::EndTreatment::Phosphorylate,
    });
    let intermediate = verify(&candidate, &target);
    assert_eq!(
        intermediate.status,
        Status::Incomplete,
        "{:?}",
        intermediate.findings
    );
    assert!(
        intermediate
            .findings
            .iter()
            .any(|f| f.contains("only as the final"))
    );
    assert!(intermediate.step_reports.is_empty());
}
