use collagene_core::topo as t;
use collagene_model::{Molecule, Topology};
use t::physical::{self as p, projection as v};
fn input(role: &str, dna: &str) -> t::Input {
    t::Input {
        role: role.into(),
        molecule: Molecule::dna(role, dna, Topology::Linear).unwrap(),
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
fn prepared(id: &str, role: &str, n: usize, activated: bool) -> t::PreparedStrand {
    t::PreparedStrand {
        id: id.into(),
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
                activation_id: id.into(),
            }
        } else {
            p::ThreePrime::Hydroxyl
        },
    }
}
fn ta() -> (Vec<t::Input>, t::Request) {
    let inputs = vec![
        input("vt", "AGGGGCATCCCTT"),
        input("vb", "AGGGATGCCCCTT"),
        input("it", "CGAATGA"),
        input("ib", "CATTCGA"),
    ];
    let request=t::Request{recipe:t::GenerationRecipe{schema_version:1,materials:vec![t::MaterialRecipe::Prepared{declaration_id:"explicit prepared material".into(),assumption:"Supplied strands and activation are declared, not inferred from these sequence references".into(),strands:vec![prepared("vt","vt",13,true),prepared("vb","vb",13,true),prepared("it","it",7,false),prepared("ib","ib",7,false)],pairing:vec![pair(0,0,1,0,12),pair(2,0,3,0,6)]}]},steps:vec![p::Step::Transfer(p::Transfer{donor: p::Call {strand:0,index:12},acceptor:p::Call{strand:2,index:0},geometry:p::Geometry::ThreePrimeTA}),p::Step::Transfer(p::Transfer{donor: p::Call {strand:1,index:12},acceptor:p::Call{strand:3,index:0},geometry:p::Geometry::ThreePrimeTA})],engine_limits:p::Limits::default(),projection_options:v::Options::default()};
    (inputs, request)
}
fn selected(inputs: &[t::Input], r: &t::Request) -> t::Selection {
    t::Selection {
        batch_revision: t::batch_revision(inputs, r).unwrap(),
        component: 0,
        reference_seed: p::Call {
            strand: 0,
            index: 0,
        },
    }
}
#[test]
fn literal_ta_contour_keeps_physical_graph_and_two_unsealed_nicks() {
    let (i, r) = ta();
    let s = selected(&i, &r);
    let product = t::product(&i, &r, &s).unwrap();
    let m = product.reference.as_ref().unwrap();
    assert_eq!(m.sequence.as_str(), "AGGGGCATCCCTTCGAATGA");
    assert_eq!(m.topology, Topology::Circular);
    assert!(collagene_core::chemistry::has_unmapped_physical_evidence(m));
    let v::Outcome::Proven(view) = &product.projection.outcome else {
        panic!()
    };
    assert_eq!(view.nicks.len(), 2);
    assert!(view.nicks.iter().all(
        |n| n.three_prime == p::ThreePrime::Hydroxyl && n.five_prime == p::FivePrime::Phosphate
    ));
    assert!(product.preparation.not_represented.is_empty());
    assert_eq!(product.preparation.maps.len(), 4);
    t::replay_product(&i, &r, &s, &product).unwrap();
}
#[test]
fn exact_target_match_cannot_silently_seal_nicks_or_replace_source_seed() {
    let (i, r) = ta();
    let s = selected(&i, &r);
    let target = Molecule::dna("independent", "AGGGGCATCCCTTCGAATGA", Topology::Circular).unwrap();
    let strict = t::check_target(&i, &r, &s, &target, &t::AcceptancePolicy::default()).unwrap();
    assert!(matches!(strict.outcome, t::TargetOutcome::Rejected { .. }));
    let permissive = t::AcceptancePolicy {
        allow_nicks: true,
        ..Default::default()
    };
    let checked = t::check_target(&i, &r, &s, &target, &permissive).unwrap();
    assert!(matches!(
        checked.outcome,
        t::TargetOutcome::AcceptedUnderDeclaredGraph
    ));
    assert_ne!(strict.batch_revision, checked.batch_revision);
    let rotated = Molecule::dna("rotated", "GGCATCCCTTCGAATGAAGG", Topology::Circular).unwrap();
    assert!(matches!(
        t::check_target(&i, &r, &s, &rotated, &permissive)
            .unwrap()
            .outcome,
        t::TargetOutcome::Rejected { .. }
    ));
    let rotate = t::AcceptancePolicy {
        allow_circular_rotation: true,
        ..permissive
    };
    assert!(matches!(
        t::check_target(&i, &r, &s, &rotated, &rotate)
            .unwrap()
            .outcome,
        t::TargetOutcome::AcceptedUnderDeclaredGraph
    ));
}
#[test]
fn recipe_source_roles_chemistry_and_selection_are_exact_pins() {
    let (mut i, mut r) = ta();
    let s = selected(&i, &r);
    i[0].molecule
        .metadata
        .insert("supplier:activated".into(), "true".into());
    assert!(t::product(&i, &r, &s).is_err());
    let new = selected(&i, &r);
    if let t::MaterialRecipe::Prepared { strands, .. } = &mut r.recipe.materials[0] {
        strands[0].three_prime = p::ThreePrime::Hydroxyl;
    }
    assert!(t::product(&i, &r, &new).is_err());
    let noactivation = t::simulate(&i, &r).unwrap();
    assert!(!noactivation.physical.all_requested_steps_applied);
    assert!(
        noactivation
            .physical
            .product
            .released_activations
            .iter()
            .all(|a| a != "vt")
    );
}
#[test]
fn circular_and_reverse_intervals_preserve_correspondence_without_inventing_waste() {
    let mut i = vec![input("reference", "ACGTTAGC")];
    i[0].molecule.topology = Topology::Circular;
    let mut a = prepared("a", "reference", 4, false);
    a.segments[0].start = 6;
    a.segments[0].reverse = true;
    let r = t::Request {
        recipe: t::GenerationRecipe {
            schema_version: 1,
            materials: vec![t::MaterialRecipe::Prepared {
                declaration_id: "declared".into(),
                assumption:
                    "only a reference correspondence, not a manufacturing conservation claim".into(),
                strands: vec![a],
                pairing: vec![],
            }],
        },
        steps: vec![],
        engine_limits: p::Limits::default(),
        projection_options: v::Options::default(),
    };
    let report = t::simulate(&i, &r).unwrap();
    assert_eq!(
        report.preparation.substrates.strands[0].sequence.as_str(),
        "GTGC"
    );
    assert_eq!(report.preparation.not_represented.len(), 1);
    assert_eq!(
        (
            report.preparation.not_represented[0].start,
            report.preparation.not_represented[0].end
        ),
        (2, 6)
    );
    assert_eq!(report.physical.source_bases, 4);
}
#[test]
fn unique_exact_pcr_and_explicit_tailing_feed_actual_strands_not_metadata() {
    let inputs = vec![input("template", "GACTAGTCATGACCTGATCGTACGAGTCCGATGCA")];
    let forward = collagene_core::analysis::Primer::new("f", "GACTAGTC", "").unwrap();
    let reverse = collagene_core::analysis::Primer::new("r", "TGCATCGG", "").unwrap();
    let request = t::Request {
        recipe: t::GenerationRecipe {
            schema_version: 1,
            materials: vec![t::MaterialRecipe::ExactPcr {
                material_id: "new PCR".into(),
                template_role: "template".into(),
                forward,
                reverse,
                tail: Some(collagene_core::tailing::TailingBase::A),
            }],
        },
        steps: vec![],
        engine_limits: p::Limits::default(),
        projection_options: v::Options::default(),
    };
    let report = t::simulate(&inputs, &request).unwrap();
    let s = &report.preparation.substrates;
    assert_eq!(
        s.strands[0].sequence.as_str(),
        "GACTAGTCATGACCTGATCGTACGAGTCCGATGCAA"
    );
    assert_eq!(
        s.strands[1].sequence.as_str(),
        "TGCATCGGACTCGTACGATCAGGTCATGACTAGTCA"
    );
    assert_eq!(s.pairing[0].length, 35);
    assert_eq!(
        report
            .physical
            .product
            .unpaired
            .iter()
            .map(|u| u.end - u.start)
            .sum::<usize>(),
        2
    );
    assert_eq!(report.preparation.maps.len(), 4);
    let mut forged = inputs.clone();
    forged[0]
        .molecule
        .metadata
        .insert("source:legacy:chemistry_unmapped".into(), "true".into());
    assert!(t::simulate(&forged, &request).is_err());
}

#[test]
fn pcr_search_work_is_reserved_for_all_passes_and_materials_before_scanning() {
    let inputs = vec![input("template", "GACTAGTCATGACCTGATCGTACGAGTCCGATGCA")];
    let material = t::MaterialRecipe::ExactPcr {
        material_id: "first PCR".into(),
        template_role: "template".into(),
        forward: collagene_core::analysis::Primer::new("f", "GACTAGTC", "").unwrap(),
        reverse: collagene_core::analysis::Primer::new("r", "TGCATCGG", "").unwrap(),
        tail: None,
    };
    let mut request = t::Request {
        recipe: t::GenerationRecipe {
            schema_version: 1,
            materials: vec![material.clone()],
        },
        steps: vec![],
        engine_limits: p::Limits::default(),
        projection_options: v::Options::default(),
    };
    // 35 template calls × (8 + 8) primer calls × two directions × two passes.
    let cost = 2240;
    request.engine_limits.max_work_cells = cost - 1;
    assert!(matches!(
        t::source_preflight(&inputs, &request),
        Err(collagene_model::Error::TooLarge)
    ));
    request.engine_limits.max_work_cells = cost;
    t::source_preflight(&inputs, &request).unwrap();
    request.recipe.materials.push(material);
    request.engine_limits.max_work_cells = 2 * cost - 1;
    assert!(matches!(
        t::source_preflight(&inputs, &request),
        Err(collagene_model::Error::TooLarge)
    ));
    request.engine_limits.max_work_cells = 2 * cost;
    t::source_preflight(&inputs, &request).unwrap();

    // Near-matching inputs would make naive searches expensive. The borrowed
    // preflight rejects their 800,000,000-cell upper bound without executing PCR.
    let long_inputs = vec![input("template", &"A".repeat(100_000))];
    request.recipe.materials = vec![t::MaterialRecipe::ExactPcr {
        material_id: "large PCR".into(),
        template_role: "template".into(),
        forward: collagene_core::analysis::Primer::new("f", &("A".repeat(999) + "C"), "").unwrap(),
        reverse: collagene_core::analysis::Primer::new("r", &("T".repeat(999) + "G"), "").unwrap(),
        tail: None,
    }];
    request.engine_limits.max_work_cells = 128_000_000;
    assert!(matches!(
        t::source_preflight(&long_inputs, &request),
        Err(collagene_model::Error::TooLarge)
    ));
}
#[test]
fn borrowed_native_trees_and_preparation_caps_fail_before_hashing() {
    let (mut i, mut r) = ta();
    let mut location = collagene_model::Location::range(0, 1);
    for _ in 0..128 {
        location = collagene_model::Location::Complement(Box::new(location));
    }
    i[0].molecule.features.push(collagene_model::Feature {
        kind: "misc_feature".into(),
        location,
        qualifiers: Default::default(),
        needs_review: false,
    });
    assert!(t::batch_revision(&i, &r).is_err());
    let (i, _) = ta();
    if let t::MaterialRecipe::Prepared { strands, .. } = &mut r.recipe.materials[0] {
        strands[0].segments = vec![strands[0].segments[0].clone(); 2049];
    }
    assert!(t::batch_revision(&i, &r).is_err());
}
#[test]
fn two_orientation_requests_retain_all_graphs_and_materialize_independently() {
    let (i, mut r) = ta();
    r.steps.clear();
    let all = t::inspect_orientations(&i, &r, [0, 1], [2, 3], false).unwrap();
    assert!(all.reference.orientation_search_complete);
    assert!(!all.reference.biochemical_mixture_exhaustive);
    assert_eq!(all.reference.candidates.len(), 2);
    assert_eq!(all.candidate_requests.len(), 2);
    for (request, candidate) in all.candidate_requests.iter().zip(&all.reference.candidates) {
        let s = selected(&i, request);
        let product = t::product(&i, request, &s).unwrap();
        assert_eq!(product.projection.physical, *candidate);
        assert_eq!(product.projection.physical.source_bases, 40);
        assert!(product.reference.is_some());
    }
    assert_ne!(
        all.candidate_requests[0].steps,
        all.candidate_requests[1].steps
    );
    assert!(
        t::inspect_orientations(&i, &all.candidate_requests[0], [0, 1], [2, 3], false).is_err()
    );
}
#[test]
fn directional_displacement_is_preserved_and_target_comparison_withheld() {
    let inputs = vec![
        input("vt", "AAGGGGCATCCCTT"),
        input("vb", "GGTGAAGGGATGCCCCTT"),
        input("it", "CACCATGACG"),
        input("ib", "CGTCATGGTG"),
    ];
    let mut r = ta().1;
    if let t::MaterialRecipe::Prepared {
        strands, pairing, ..
    } = &mut r.recipe.materials[0]
    {
        *strands = vec![
            prepared("vt", "vt", 14, true),
            prepared("vb", "vb", 18, true),
            prepared("it", "it", 10, false),
            prepared("ib", "ib", 10, false),
        ];
        *pairing = vec![pair(0, 0, 1, 4, 14), pair(2, 0, 3, 0, 10)];
    }
    let p::Step::Transfer(first) = &mut r.steps[0] else {
        panic!("transfer")
    };
    first.geometry = p::Geometry::DirectionalCacc;
    first.donor.index = 13;
    let p::Step::Transfer(second) = &mut r.steps[1] else {
        panic!("transfer")
    };
    second.geometry = p::Geometry::Flush;
    second.donor.index = 17;
    let s = selected(&inputs, &r);
    let product = t::product(&inputs, &r, &s).unwrap();
    assert!(product.reference.is_none());
    assert!(matches!(
        product.projection.outcome,
        v::Outcome::Withheld(v::Withheld::UnpairedNucleotides { count: 4 })
    ));
    assert_eq!(
        product.projection.physical.product.chains[1]
            .sequence
            .as_str(),
        "GGTGAAGGGATGCCCCTTCGTCATGGTG"
    );
    let target = Molecule::dna(
        "independent desired reference",
        "AAGGGGCATCCCTTCACCATGACG",
        Topology::Circular,
    )
    .unwrap();
    let policy = t::AcceptancePolicy {
        allow_nicks: true,
        allow_unknown_terminal_chemistry: true,
        allow_unapplied_steps: true,
        allow_other_components: true,
        allow_circular_rotation: true,
        allow_reverse_complement: true,
    };
    assert!(matches!(
        t::check_target(&inputs, &r, &s, &target, &policy)
            .unwrap()
            .outcome,
        t::TargetOutcome::Withheld { .. }
    ));
}
