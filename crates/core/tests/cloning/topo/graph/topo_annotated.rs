use collagene_core::cloning::topo::physical::{self as p, projection as v};
use collagene_core::{
    cloning::topo as t,
    planning::plan::{DesignPlan, Operation, PlanInput, Source, TopoSource},
    verification::{engine::verify, evidence::Status, output::export_protocol_html},
};
use collagene_model::{Feature, Location, Molecule, Topology};
fn source(sequence: &str) -> Molecule {
    Molecule::dna("actual source", sequence, Topology::Linear).unwrap()
}
fn feature(location: Location) -> Feature {
    Feature {
        kind: "CDS".into(),
        location,
        qualifiers: vec![
            ("label".into(), Some("source annotation".into())),
            ("translation".into(), Some("SOURCE_ASSERTION".into())),
            ("codon_start".into(), Some("1".into())),
        ],
        needs_review: false,
    }
}
fn selection(inputs: &[t::Input], request: &t::Request) -> t::Selection {
    t::Selection {
        batch_revision: t::batch_revision(inputs, request).unwrap(),
        component: 0,
        reference_seed: p::Call {
            strand: 0,
            index: 0,
        },
    }
}
fn prepared() -> (Vec<t::Input>, t::Request) {
    let mut molecule = source("AAAACCCC");
    molecule.features.push(feature(Location::range(0, 4)));
    molecule
        .features
        .push(feature(Location::Unresolved("REMOTE.1:1..4".into())));
    let inputs = vec![t::Input {
        role: "reference".into(),
        molecule,
    }];
    let strands = (0..2)
        .map(|i| t::PreparedStrand {
            id: format!("strand{i}"),
            segments: vec![t::SourceInterval {
                source_role: "reference".into(),
                start: 0,
                length: 8,
                reverse: i == 1,
            }],
            five_prime: p::FivePrime::Hydroxyl,
            three_prime: p::ThreePrime::Hydroxyl,
        })
        .collect();
    let request = t::Request {
        recipe: t::GenerationRecipe {
            schema_version: 1,
            materials: vec![t::MaterialRecipe::Prepared {
                declaration_id: "explicit strands".into(),
                assumption: "Source correspondence is not chemical ancestry".into(),
                strands,
                pairing: vec![p::PairSpan {
                    strand_a: 0,
                    start_a: 0,
                    strand_b: 1,
                    start_b: 0,
                    length: 8,
                }],
            }],
        },
        steps: vec![],
        engine_limits: p::Limits::default(),
        projection_options: v::Options::default(),
    };
    (inputs, request)
}

#[test]
fn canonical_product_projects_annotations_and_retains_unmapped_source_evidence() {
    let (inputs, request) = prepared();
    let selected = selection(&inputs, &request);
    let annotated = t::product(&inputs, &request, &selected).unwrap();
    assert_eq!(annotated.reference.as_ref().unwrap().features.len(), 1);
    assert!(collagene_core::chemistry::has_unmapped_physical_evidence(
        annotated.reference.as_ref().unwrap()
    ));
    let evidence = annotated.annotations.as_ref().unwrap();
    assert_eq!(evidence.sources[0].projected, 1);
    assert_eq!(evidence.sources[0].unmapped.len(), 1);
    assert_eq!(evidence.sources[0].unmapped[0].source_index, 1);
    assert!(evidence.sources[0].unmapped[0].feature.needs_review);
    assert_eq!(
        annotated.batch_revision,
        t::product_revision(&inputs, &request, &selected).unwrap()
    );
    let target = source("AAAACCCC");
    let policy = t::AcceptancePolicy::default();
    let checked = t::check_target(&inputs, &request, &selected, &target, &policy).unwrap();
    assert_eq!(
        checked.batch_revision,
        collagene_core::cloning::topo::worker::target_revision(
            &inputs, &request, &selected, &target, &policy
        )
        .unwrap()
    );
    assert!(matches!(
        checked.outcome,
        t::TargetOutcome::AcceptedUnderDeclaredGraph
    ));
    let plan = DesignPlan {
        name: "annotated TOPO".into(),
        target_revision: target.revision(),
        inputs: vec![PlanInput {
            id: "input".into(),
            revision: inputs[0].molecule.revision(),
            molecule: inputs[0].molecule.clone(),
        }],
        steps: vec![Operation::Topo {
            sources: vec![TopoSource {
                role: "reference".into(),
                source: Source::Input("input".into()),
            }],
            request,
            selection: selected,
        }],
    };
    let report = verify(&plan, &target);
    assert_eq!(report.status, Status::Incomplete);
    assert_eq!(report.products[0].features.len(), 1);
    assert_eq!(
        report.step_reports[0]
            .topo
            .as_ref()
            .unwrap()
            .annotations
            .as_ref()
            .unwrap()
            .sources[0]
            .projected,
        1
    );
    assert_eq!(
        report.step_reports[0]
            .topo
            .as_ref()
            .unwrap()
            .reference
            .as_ref()
            .unwrap()
            .features
            .len(),
        1
    );
    let html = export_protocol_html(&plan, &target);
    assert!(html.contains("1 projected (0 require review)"));
    assert!(html.contains("1 original annotations retained unmapped"));
    assert!(html.contains("neither inherited chemical ancestry"));
}

#[test]
fn pcr_template_features_require_review_and_unproven_tailed_contours_get_no_annotation_reference() {
    let mut template = source("GACTAGTCATGACCTGATCGTACGAGTCCGATGCA");
    template.features.push(feature(Location::range(0, 8)));
    let inputs = vec![t::Input {
        role: "template".into(),
        molecule: template,
    }];
    let mut request = t::Request {
        recipe: t::GenerationRecipe {
            schema_version: 1,
            materials: vec![t::MaterialRecipe::ExactPcr {
                material_id: "new PCR".into(),
                template_role: "template".into(),
                forward: collagene_core::sequence::analysis::Primer::new("F", "GACTAGTC", "GG")
                    .unwrap(),
                reverse: collagene_core::sequence::analysis::Primer::new("R", "TGCATCGG", "CC")
                    .unwrap(),
                tail: None,
            }],
        },
        steps: vec![],
        engine_limits: p::Limits::default(),
        projection_options: v::Options::default(),
    };
    let annotated = t::product(&inputs, &request, &selection(&inputs, &request)).unwrap();
    let feature = &annotated.reference.as_ref().unwrap().features[0];
    assert_eq!(feature.location.spans(), vec![(2, 10, false)]);
    assert!(feature.needs_review);
    assert!(
        !feature
            .qualifiers
            .iter()
            .any(|(key, _)| key == "translation" || key == "codon_start")
    );
    if let t::MaterialRecipe::ExactPcr { tail, .. } = &mut request.recipe.materials[0] {
        *tail = Some(collagene_core::cloning::ligation::tailing::TailingBase::A);
    }
    let selected = selection(&inputs, &request);
    let withheld = t::product(&inputs, &request, &selected).unwrap();
    assert!(withheld.reference.is_none() && withheld.annotations.is_none());
    assert_eq!(withheld.projection.physical.source_bases, 80);
    let checked = t::check_target(
        &inputs,
        &request,
        &selected,
        &source("GGGACTAGTCATGACCTGATCGTACGAGTCCGATGCAGG"),
        &t::AcceptancePolicy::default(),
    )
    .unwrap();
    assert!(matches!(checked.outcome, t::TargetOutcome::Withheld { .. }));
    assert!(checked.product.annotations.is_none());
}
