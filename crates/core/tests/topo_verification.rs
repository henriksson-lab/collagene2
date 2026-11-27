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
        transfers:vec![p::Transfer{donor:0,acceptor:2,geometry:if directional {p::Geometry::DirectionalCacc}else{p::Geometry::ThreePrimeTA}},p::Transfer{donor:1,acceptor:3,geometry:if directional {p::Geometry::Flush}else{p::Geometry::ThreePrimeTA}}],engine_limits:p::Limits::default(),projection_options:v::Options::default()};
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
            steps: vec![Operation::Topo {
                sources,
                request,
                selection,
            }],
        },
        target,
    )
}

#[test]
fn matching_contour_retains_full_graph_and_nicks_without_implicit_target_acceptance() {
    let (plan, target) = plan(false);
    let report = verify(&plan, &target);
    assert_eq!(report.status, Status::Incomplete, "{:?}", report.findings);
    assert_eq!(report.products.len(), 1);
    assert_eq!(report.products[0].sequence, target.sequence);
    assert_eq!(
        report.products[0].metadata[collagene_core::recombination::REFERENCE_ONLY_KEY],
        "true"
    );
    assert!(collagene_core::chemistry::has_unmapped_physical_evidence(
        &report.products[0]
    ));
    let evidence = report.step_reports[0].topo.as_ref().unwrap();
    assert_eq!(report.step_reports[0].step_index, 0);
    assert_eq!(evidence.projection.physical.source_bases, 40);
    assert_eq!(evidence.projection.physical.product.components.len(), 1);
    let v::Outcome::Proven(view) = &evidence.projection.outcome else {
        panic!("expected contour")
    };
    assert_eq!(view.nicks.len(), 2);
    assert!(evidence.reference.is_some());
    // The appended full-evidence field survives the actual binary worker seam.
    let request = collagene_core::protocol::Request {
        id: 9,
        revision: target.revision(),
        job: collagene_core::protocol::AnalysisJob::Verify {
            target: target.clone(),
            plan: plan.clone(),
        },
    };
    let mut bytes = collagene_core::protocol::encode_request(&request).unwrap();
    let response = collagene_core::protocol::decode_response(
        &collagene_core::protocol::process_request(&bytes).unwrap(),
    )
    .unwrap();
    let collagene_core::protocol::AnalysisResult::Verification(round) = response.result else {
        panic!("expected verifier report")
    };
    assert_eq!(round.step_reports.len(), 1);
    assert_eq!(
        round.step_reports[0]
            .topo
            .as_ref()
            .unwrap()
            .projection
            .physical
            .source_bases,
        40
    );
    bytes[4..6].copy_from_slice(&24u16.to_le_bytes());
    assert!(collagene_core::protocol::process_request(&bytes).is_err());
    let html = export_protocol_html(&plan, &target);
    assert!(html.contains("Declared source-corresponding strands and activation"));
    assert!(html.contains("TOPO physical step evidence") && html.contains("retained nicks"));
}

#[test]
fn directional_flap_withholds_dna_but_retains_complete_physical_step_report() {
    let (plan, target) = plan(true);
    let report = verify(&plan, &target);
    assert_eq!(report.status, Status::Incomplete, "{:?}", report.findings);
    assert!(report.products.is_empty());
    let evidence = report.step_reports[0].topo.as_ref().unwrap();
    assert!(evidence.reference.is_none());
    assert!(matches!(
        evidence.projection.outcome,
        v::Outcome::Withheld(v::Withheld::UnpairedNucleotides { count: 4 })
    ));
    assert_eq!(
        evidence
            .projection
            .physical
            .product
            .unpaired
            .iter()
            .map(|r| r.end - r.start)
            .sum::<usize>(),
        4
    );
}

#[test]
fn exact_roles_selection_and_source_dependencies_cannot_be_reinterpreted_or_promote_preview_chemistry()
 {
    let (mut changed, target) = plan(false);
    if let Operation::Topo { sources, .. } = &mut changed.steps[0] {
        sources[0].source = Source::Input("vb".into());
    }
    assert_eq!(verify(&changed, &target).status, Status::Invalid);
    let (mut duplicate, target) = plan(false);
    if let Operation::Topo { sources, .. } = &mut duplicate.steps[0] {
        sources[1].role = "vt".into();
    }
    assert_eq!(verify(&duplicate, &target).status, Status::Invalid);
    let (mut forward, target) = plan(false);
    if let Operation::Topo { sources, .. } = &mut forward.steps[0] {
        sources[0].source = Source::Step(0);
    }
    assert_eq!(verify(&forward, &target).status, Status::Invalid);
    let (mut downstream, target) = plan(false);
    downstream.steps.push(Operation::TreatEnds {
        template: Source::Step(0),
        treatment: collagene_core::chemistry::EndTreatment::Phosphorylate,
    });
    let report = verify(&downstream, &target);
    assert_eq!(report.status, Status::Incomplete);
    assert_eq!(report.products.len(), 1);
    assert_eq!(report.step_reports.len(), 1);
}

#[test]
fn html_pcr_preparation_includes_actual_ordered_oligos_and_explicit_physical_assumptions() {
    let template = molecule("GACTAGTCATGACCTGATCGTACGAGTCCGATGCA");
    let sources = vec![TopoSource {
        role: "template".into(),
        source: Source::Input("template".into()),
    }];
    let request = t::Request {
        recipe: t::GenerationRecipe {
            schema_version: 1,
            materials: vec![t::MaterialRecipe::ExactPcr {
                material_id: "new PCR".into(),
                template_role: "template".into(),
                forward: collagene_core::analysis::Primer::new("f", "GACTAGTC", "GG").unwrap(),
                reverse: collagene_core::analysis::Primer::new("r", "TGCATCGG", "CC").unwrap(),
                tail: None,
            }],
        },
        transfers: vec![],
        engine_limits: p::Limits::default(),
        projection_options: v::Options::default(),
    };
    let bound = topo_inputs(&sources, &[&template], &request).unwrap();
    let selection = t::Selection {
        batch_revision: t::batch_revision(&bound, &request).unwrap(),
        component: 0,
        reference_seed: p::Call {
            strand: 0,
            index: 0,
        },
    };
    let target = molecule("GGGACTAGTCATGACCTGATCGTACGAGTCCGATGCAGG");
    let plan = DesignPlan {
        name: "PCR preparation".into(),
        target_revision: target.revision(),
        inputs: vec![PlanInput {
            id: "template".into(),
            revision: template.revision(),
            molecule: template,
        }],
        steps: vec![Operation::Topo {
            sources,
            request,
            selection,
        }],
    };
    let report = verify(&plan, &target);
    assert_eq!(report.status, Status::Incomplete, "{:?}", report.findings);
    assert_eq!(report.products[0].sequence, target.sequence);
    let html = export_protocol_html(&plan, &target);
    assert!(html.contains("<code>GGGACTAGTC</code>"));
    assert!(html.contains("<code>CCTGCATCGG</code>"));
    assert!(html.contains("unique exact PCR geometry"));
}
