use collagene_core::{
    Feature, Location, Molecule, Topology,
    assembly_junctions::{self as junctions, Method, Options},
    golden_gate::{GoldenGateEnzyme, GoldenGateOptions},
    golden_gate_primer_design as primers, ligation_recipe, prepared_golden_gate,
    protocol::{self, AnalysisJob, AnalysisResult, Request},
    verification::{self, DesignPlan, Operation, PlanInput, Source, Status},
};
const A: &str = "ATGACCCCGGGGTTTCCAGA";
const B: &str = "ATGAGTTAGCGACCTTGGAA";
const AB: &str = "ATGACCCCGGGGTTTCCAGAATGAGTTAGCGACCTTGGAA";

fn raw(reverse: bool) -> (DesignPlan, Molecule) {
    // Two Type IIS sites: forward cuts11/15, reverse cuts31/35.
    // The right top fusion is ATGA and the opposite bottom carries its complement.
    let input = Molecule::dna(
        "original amplicon",
        &format!("ATGCGGTCTCA{A}ATGATGAGACCGCAT"),
        Topology::Linear,
    )
    .unwrap();
    let target = Molecule::dna(
        "independent target",
        if reverse { "TCTGGAAACCCCGGGGTCAT" } else { A },
        Topology::Circular,
    )
    .unwrap();
    let plan = DesignPlan {
        name: "raw cut circle".into(),
        inputs: vec![PlanInput {
            id: "source".into(),
            revision: input.revision(),
            molecule: input,
        }],
        target_revision: target.revision(),
        steps: vec![
            Operation::DeclareSourceModifications {
                template: Source::Input("source".into()),
                assumptions: collagene_core::modified_duplex::GenerationAssumptions {
                    dam: collagene_core::methylation::MethylationState::Unmodified,
                    dcm: collagene_core::methylation::MethylationState::Unmodified,
                    cpg: collagene_core::methylation::MethylationState::Unmodified,
                },
            },
            Operation::GoldenGate {
                templates: vec![Source::Step(0)],
                options: GoldenGateOptions {
                    max_fragments: 1,
                    ..Default::default()
                },
                candidate_index: usize::from(reverse),
            },
        ],
    };
    (plan, target)
}
fn options() -> Options {
    Options {
        flank_bases: 4,
        ..Default::default()
    }
}

#[test]
fn declared_raw_cuts_closure_and_window_codon_use_current_wire_and_reject_noncurrent_headers() {
    let (plan, target) = raw(false);
    let report = junctions::inspect_plan(&plan, &target, &options()).unwrap();
    assert_eq!(
        report.verification.status,
        Status::ConsistentUnderAssumptions
    );
    assert!(report.complete);
    let Some(Method::GoldenGate {
        path,
        candidate_index,
        ..
    }) = &report.method
    else {
        panic!("missing actual path")
    };
    assert_eq!(*candidate_index, 0);
    assert_eq!(path.len(), 1);
    assert_eq!((path[0].top_start, path[0].bottom_start), (11, 15));
    assert_eq!((path[0].output_start, path[0].output_end), (0, 20));
    assert!(matches!(&path[0].template, Source::Step(0)));
    assert_eq!(
        path[0].template_revision,
        report.verification.products[0].revision()
    );
    assert_eq!(
        report.verification.products[0].sequence,
        plan.inputs[0].molecule.sequence
    );
    assert!(!path[0].reverse);
    let j = &report.junctions[0];
    assert!(j.closing);
    assert_eq!(j.boundary, 0);
    assert_eq!(j.window.dna, "CAGAATGA");
    assert_eq!(j.frames[0].crossing_codon.as_ref().unwrap().dna, "AAT");
    assert_eq!(
        postcard::to_allocvec(report.method.as_ref().unwrap()).unwrap()[0],
        3
    );
    let revision = junctions::inspection_revision(&plan, &target, &options()).unwrap();
    let bytes = protocol::encode_request(&Request {
        id: 1,
        revision: revision.clone(),
        job: AnalysisJob::InspectAssemblyJunctions {
            plan: plan.clone(),
            target: target.clone(),
            options: options(),
        },
    })
    .unwrap();
    let bytes = protocol::process_request(&bytes).unwrap();
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[42..]).unwrap();
    assert_eq!(tag, 88);
    let AnalysisResult::AssemblyJunctionInspection(decoded) =
        protocol::decode_response(&bytes).unwrap().result
    else {
        panic!("wrong result")
    };
    decoded.validate_binding(&revision).unwrap();
    assert_eq!(
        u16::from_le_bytes(bytes[4..6].try_into().unwrap()),
        protocol::WORKER_SCHEMA_VERSION
    );
    let mut old = bytes;
    old[4..6].copy_from_slice(&(protocol::WORKER_SCHEMA_VERSION - 1).to_le_bytes());
    assert!(protocol::decode_response(&old).is_err());
    let mut changed = report.clone();
    let Some(Method::GoldenGate { path, .. }) = &mut changed.method else {
        unreachable!()
    };
    path[0].output_end = 19;
    assert!(changed.validate().is_err());
    let mut wrong = target.clone();
    wrong.sequence = collagene_core::Sequence::dna(AB).unwrap();
    wrong.features.push(Feature {
        kind: "CDS".into(),
        location: Location::range(0, 40),
        qualifiers: vec![],
        needs_review: false,
    });
    let mut plan = plan;
    plan.target_revision = wrong.revision();
    let failed = junctions::inspect_plan(&plan, &wrong, &options()).unwrap();
    assert_eq!(failed.verification.status, Status::Invalid);
    assert_eq!(failed.product_length, Some(20));
    assert_eq!(failed.junctions[0].window.dna, "CAGAATGA");
    assert!(
        failed.junctions[0].cds.is_empty(),
        "target-only CDS must not become actual product annotation"
    );
}

#[test]
fn reverse_path_uses_the_retained_bottom_strand_not_reverse_complement_of_the_whole_top() {
    let (plan, target) = raw(true);
    let report = junctions::inspect_plan(&plan, &target, &options()).unwrap();
    assert_eq!(
        report.verification.status,
        Status::ConsistentUnderAssumptions
    );
    // Actual reverse top is RC(A[4..] + ATGA), not RC(A). Target comparison permits circle rotation.
    assert_eq!(
        report
            .verification
            .products
            .last()
            .unwrap()
            .sequence
            .as_str(),
        "TCATTCTGGAAACCCCGGGG"
    );
    assert_ne!(
        report
            .verification
            .products
            .last()
            .unwrap()
            .sequence
            .as_str(),
        target.sequence.as_str()
    );
    let Some(Method::GoldenGate { path, .. }) = &report.method else {
        panic!("missing path")
    };
    assert!(path[0].reverse);
    assert_eq!((path[0].top_start, path[0].bottom_start), (11, 15));
    assert_eq!(report.junctions[0].window.dna, "GGGGTCAT");
}

fn prepared(circular: bool, capped: bool) -> (DesignPlan, Molecule) {
    let target = Molecule::dna(
        "independent A then B",
        AB,
        if circular {
            Topology::Circular
        } else {
            Topology::Linear
        },
    )
    .unwrap();
    let inputs = [("a", A), ("b", B)]
        .into_iter()
        .map(|(id, dna)| {
            let molecule = Molecule::dna(id, dna, Topology::Linear).unwrap();
            PlanInput {
                id: id.into(),
                revision: molecule.revision(),
                molecule,
            }
        })
        .collect::<Vec<_>>();
    let designed = primers::design(
        &inputs,
        &target,
        &primers::Request {
            name: "prepare originals".into(),
            fragments: inputs
                .iter()
                .map(|i| primers::Fragment {
                    role: i.id.clone(),
                    input_id: i.id.clone(),
                    start: 0,
                    end: 20,
                    reverse: false,
                    forward_length: 8,
                    reverse_length: 8,
                    right_overhang: "ATGA".into(),
                })
                .collect(),
            enzyme: GoldenGateEnzyme::BsaI,
            clamp: "ATGC".into(),
            spacer: "A".into(),
        },
    )
    .unwrap();
    let request = prepared_golden_gate::Request {
        recipe: ligation_recipe::Recipe {
            name: "prepared pool".into(),
            fragments: designed
                .fragments
                .into_iter()
                .map(|source| ligation_recipe::Fragment {
                    source,
                    treatment: None,
                })
                .collect(),
            circular,
            allow_nicked: false,
        },
        options: GoldenGateOptions {
            circular,
            min_fragments: if capped { 1 } else { 2 },
            max_fragments: 2,
            max_products: if capped { 2 } else { 128 },
            ..Default::default()
        },
    };
    let report = prepared_golden_gate::inspect(&inputs, &target, &request).unwrap();
    let a = report
        .assessments
        .iter()
        .find(|a| a.target_matches)
        .unwrap();
    let selected = prepared_golden_gate::select(
        &report,
        &prepared_golden_gate::Selection {
            context_revision: report.context_revision.clone(),
            candidate_index: a.candidate_index,
            directed_sequence_key: a.directed_sequence_key.clone(),
            physical_sequence_key: a.physical_sequence_key.clone(),
        },
    )
    .unwrap();
    (selected.plan, target)
}

#[test]
fn two_prepared_parts_map_to_their_real_steps_and_caps_remain_separate_from_junction_geometry() {
    for (circular, capped) in [(true, false), (false, false), (true, true)] {
        let (plan, target) = prepared(circular, capped);
        let report = junctions::inspect_plan(&plan, &target, &options()).unwrap();
        assert_eq!(report.product_length, Some(40));
        assert_eq!(
            report.verification.status,
            if capped {
                Status::Incomplete
            } else {
                Status::ConsistentUnderAssumptions
            }
        );
        let Some(Method::GoldenGate {
            path,
            enumeration_complete,
            ..
        }) = &report.method
        else {
            panic!("missing path")
        };
        assert_eq!(*enumeration_complete, !capped);
        assert_eq!(path.len(), 2);
        assert!(matches!(path[0].template, Source::Step(1)));
        assert!(matches!(path[1].template, Source::Step(3)));
        assert_eq!((path[0].output_end, path[1].output_end), (20, 40));
        assert_eq!(report.junctions[0].boundary, 20);
        assert_eq!(report.junctions[0].window.dna, "CAGAATGA");
        assert_eq!(report.junctions.len(), if circular { 2 } else { 1 });
        if circular {
            assert_eq!(report.junctions[1].boundary, 0);
            assert_eq!(report.junctions[1].window.dna, "GGAAATGA");
        }
        let clipped = junctions::inspect_plan(
            &plan,
            &target,
            &Options {
                max_junctions: 1,
                ..options()
            },
        )
        .unwrap();
        assert_eq!(clipped.omitted_junctions, usize::from(circular));
        assert_eq!(clipped.complete, !circular);
    }
}

#[test]
fn failed_final_selection_never_presents_a_preparation_as_a_final_product() {
    let (mut plan, target) = prepared(true, false);
    let Some(Operation::GoldenGate {
        candidate_index, ..
    }) = plan.steps.last_mut()
    else {
        unreachable!()
    };
    *candidate_index = 999;
    let report = junctions::inspect_plan(&plan, &target, &options()).unwrap();
    assert!(report.method.is_none());
    assert!(report.product_revision.is_none());
    assert!(report.junctions.is_empty());
    assert!(!report.complete);
    assert_eq!(report.verification.products.len() + 1, plan.steps.len());
}
