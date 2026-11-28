//! Literal reciprocal products and typed graph/native version boundaries.
use collagene_core::{
    analysis::Primer,
    chemistry::EndTreatment,
    protocol::{self, request::Request as WireRequest, schema::AnalysisJob},
    recombination::{self, ProfileSnapshot, Request, Selection, reference},
    run_lineage::{self, EventBinding, InputEventBinding},
    transfer::{self, Candidate},
    verification::{self, DesignPlan, Operation, PlanInput, Source, Status},
    workspace::{self, Workspace},
};
use collagene_model::{Molecule, Topology};
fn fixture(byproduct: bool) -> (Molecule, DesignPlan) {
    let d: String = include_str!("../../../../recombination/tests/fixtures/pdonr221_seq.txt")
        .split_whitespace()
        .collect();
    let b = "ACAAGTTTGTACAAAAAAGCAGGCTATGCCGTAAACCCAGCTTTCTTGTACAAAGTGGT";
    assert_eq!(b.len(), 59);
    let mut a = Molecule::dna("PCR", b, Topology::Linear).unwrap();
    a.features.push(collagene_model::Feature {
        kind: "misc_feature".into(),
        location: collagene_model::Location::range(25, 34),
        qualifiers: vec![("label".into(), Some("payload".into()))],
        needs_review: false,
    });
    a.original_genbank = Some(collagene_core::formats::export_genbank(&a).unwrap());
    let donor = Molecule::dna("donor", &d, Topology::Circular).unwrap();
    // Independent source-boundary oracle; neither sequence comes from a reported product.
    let (sequence, topology) = if byproduct {
        (
            format!("{}{}{}", &b[..8], &d[651..2902], &b[51..]),
            Topology::Linear,
        )
    } else {
        (
            format!("{}{}{}", &d[..651], &b[8..51], &d[2902..]),
            Topology::Circular,
        )
    };
    assert_eq!(sequence.len(), if byproduct { 2267 } else { 2553 });
    let target = Molecule::dna("independent requirement", &sequence, topology).unwrap();
    let sources = [Source::Input("a".into()), Source::Input("b".into())];
    let request = Request::new(
        ProfileSnapshot::manufacturer_context().unwrap(),
        reference::Reaction::BP,
    );
    let inputs = verification::recombination_inputs(&sources, [&a, &donor], &request).unwrap();
    let report = reference::recombine(
        &inputs,
        &request.profile.profile,
        request.reaction,
        &request.limits,
    )
    .unwrap();
    assert_eq!(report.candidates.len(), 1);
    let component_id = report.candidates[0]
        .products
        .iter()
        .position(|p| p.molecule.sequence.len() == sequence.len())
        .unwrap();
    let selection = Selection {
        batch_revision: recombination::batch_revision(&inputs, &request).unwrap(),
        candidate_id: 0,
        component_id,
    };
    let plan = DesignPlan {
        name: "pinned reference exchange".into(),
        target_revision: target.revision(),
        inputs: vec![
            PlanInput {
                id: "a".into(),
                revision: a.revision(),
                molecule: a,
            },
            PlanInput {
                id: "b".into(),
                revision: donor.revision(),
                molecule: donor,
            },
        ],
        steps: vec![Operation::RecombineReference {
            sources,
            request,
            selection,
        }],
    };
    (target, plan)
}
fn wire(magic: &[u8; 4], version: u16, v: &impl serde::Serialize) -> Vec<u8> {
    let p = postcard::to_allocvec(v).unwrap();
    let mut b = magic.to_vec();
    b.extend(version.to_le_bytes());
    b.extend((p.len() as u32).to_le_bytes());
    b.extend(blake3::hash(&p).as_bytes());
    b.extend(p);
    b
}
fn event(w: &Workspace, i: usize) -> EventBinding {
    let d = &w.documents[i];
    EventBinding {
        document_id: d.id.clone(),
        event_id: d.head_event.clone().unwrap(),
        revision: d.head.clone(),
    }
}
#[test]
fn literal_both_components_match_only_as_reference_models_and_protocols_retain_pins() {
    for byproduct in [false, true] {
        let (target, plan) = fixture(byproduct);
        let r = verification::verify(&plan, &target);
        assert_eq!(r.status, Status::Incomplete, "{:?}", r.findings);
        assert_eq!(r.products.len(), 1);
        assert!(!r.exhaustive_alternative_search);
        assert!(!r.products[0].metadata.contains_key("tracked_modifications"));
        assert!(!r.products[0].metadata.contains_key("duplex_model"));
        let html = verification::export_protocol_html(&plan, &target);
        assert!(
            html.contains("reference")
                && html.contains("candidate_id")
                && html.contains("max_comparisons")
                && html.contains("raw_sha256")
        );
        assert!(
            collagene_core::pdf::export_protocol_pdf(&plan, &target)
                .unwrap()
                .bytes
                .starts_with(b"%PDF")
        );
    }
}
#[test]
fn stale_source_profile_policy_slot_and_component_never_rebind_by_target() {
    let (target, plan) = fixture(false);
    for mode in 0..5 {
        let mut p = plan.clone();
        if mode == 0 {
            p.inputs[0].molecule.name.push('x');
            p.inputs[0].revision = p.inputs[0].molecule.revision();
        }
        if let Operation::RecombineReference {
            sources,
            request,
            selection,
        } = &mut p.steps[0]
        {
            match mode {
                1 => sources.swap(0, 1),
                2 => request.limits.max_comparisons -= 1,
                3 => {
                    request.profile.profile.edition.push('x');
                }
                4 => selection.component_id = 999,
                _ => {}
            }
        }
        assert_ne!(
            verification::verify(&p, &target).status,
            Status::ConsistentUnderAssumptions
        );
        assert!(verification::verify(&p, &target).products.is_empty());
    }
    let mut wrong = target.clone();
    wrong.sequence = collagene_model::Sequence::dna(&"A".repeat(target.sequence.len())).unwrap();
    let mut p = plan.clone();
    p.target_revision = wrong.revision();
    assert_eq!(verification::verify(&p, &wrong).status, Status::Invalid);
}
#[test]
fn reference_ancestry_cannot_gain_physical_authority_directly_or_after_overlap() {
    let (target, plan) = fixture(false);
    let primer = Primer::new("P", "ACGTACGT", "").unwrap();
    for through_overlap in [false, true] {
        for kind in 0..3 {
            let mut p = plan.clone();
            let source = if through_overlap {
                p.steps.push(Operation::Overlap {
                    fragments: vec![Source::Step(0)],
                    overlaps: vec![],
                    circular: false,
                });
                Source::Step(1)
            } else {
                Source::Step(0)
            };
            p.steps.push(match kind {
                0 => Operation::Pcr {
                    template: source,
                    forward: primer.clone(),
                    reverse: primer.clone(),
                    template_strand: collagene_core::chemistry::Strand::Top,
                },
                1 => Operation::PcrWithMismatches {
                    template: source,
                    forward: primer.clone(),
                    reverse: primer.clone(),
                    options: Default::default(),
                },
                _ => Operation::TreatEnds {
                    template: source,
                    treatment: EndTreatment::Phosphorylate,
                },
            });
            let r = verification::verify(&p, &target);
            assert_eq!(r.status, Status::Incomplete);
            assert!(r.products.len() < p.steps.len());
            assert!(
                r.products
                    .iter()
                    .all(|m| !m.metadata.contains_key("tracked_modifications"))
            );
        }
    }
}
#[test]
fn pinned_custom_profile_replays_without_installed_catalog_lookup() {
    let (target, mut p) = fixture(false);
    let molecules = [&p.inputs[0].molecule, &p.inputs[1].molecule];
    if let Operation::RecombineReference {
        sources,
        request,
        selection,
    } = &mut p.steps[0]
    {
        let mut custom = request.profile.profile.clone();
        custom.id = "explicit-custom-snapshot".into();
        custom.edition = "independent-edition".into();
        request.profile = ProfileSnapshot::new(custom).unwrap();
        let inputs = verification::recombination_inputs(sources, molecules, request).unwrap();
        selection.batch_revision = recombination::batch_revision(&inputs, request).unwrap();
    }
    let encoded = transfer::save_candidate(&Candidate::new(p)).unwrap();
    let decoded = transfer::load_candidate(&encoded).unwrap();
    assert_eq!(
        verification::verify(&decoded.plan, &target).status,
        Status::Incomplete
    );
}
#[test]
fn native_capture_is_explicitly_incomplete_has_exact_two_parents_and_replays() {
    let (target, p) = fixture(false);
    let w = Workspace::from_documents(
        "reference",
        vec![
            p.inputs[0].molecule.clone(),
            p.inputs[1].molecule.clone(),
            target,
            p.inputs[0].molecule.clone(),
        ],
        2,
    )
    .unwrap();
    let bindings = vec![
        InputEventBinding {
            input_id: "a".into(),
            event: event(&w, 0),
        },
        InputEventBinding {
            input_id: "b".into(),
            event: event(&w, 1),
        },
    ];
    assert!(run_lineage::record_run(&w, &p, event(&w, 2), bindings.clone(), false).is_err());
    let before = serde_json::to_value(&w).unwrap();
    let c = run_lineage::record_run(&w, &p, event(&w, 2), bindings, true).unwrap();
    c.validate().unwrap();
    assert_eq!(serde_json::to_value(&w).unwrap(), before);
    let run = run_lineage::load_run(&c.workspace, &c.receipt.run_id).unwrap();
    assert_eq!(
        run.plan_schema,
        collagene_core::workspace::PLAN_SCHEMA_VERSION
    );
    assert_eq!(run.assessment.status, Status::Incomplete);
    let pe = c
        .workspace
        .history
        .iter()
        .find(|e| e.id == run.products[0].event_id)
        .unwrap();
    let mut expected = vec![event(&w, 0).event_id, event(&w, 1).event_id];
    expected.sort();
    let mut actual = pe.parents.clone();
    actual.sort();
    assert_eq!(actual, expected);
    assert!(!actual.contains(&event(&w, 2).event_id));
    let loaded =
        workspace::load_workspace(&workspace::save_workspace(&c.workspace).unwrap()).unwrap();
    let replay = run_lineage::replay_run(&loaded, &run.run_id).unwrap();
    assert!(replay.products_match && replay.status_matches && replay.engine_matches);
    let other = run_lineage::record_run(
        &w,
        &p,
        event(&w, 2),
        vec![
            InputEventBinding {
                input_id: "a".into(),
                event: event(&w, 3),
            },
            InputEventBinding {
                input_id: "b".into(),
                event: event(&w, 1),
            },
        ],
        true,
    )
    .unwrap();
    assert_ne!(other.receipt.run_id, c.receipt.run_id);
    let job = WireRequest {
        id: 4,
        revision: "workspace".into(),
        job: AnalysisJob::SaveWorkspace {
            workspace: c.workspace,
        },
    };
    for v in [19, 20, 21] {
        assert!(protocol::execution::process_request(&wire(b"CGQ\0", v, &job)).is_err());
    }
}

#[test]
fn native_manifest_cannot_downcast_reference_graph_even_with_valid_step_blobs() {
    use std::io::{Cursor, Read, Write};
    let (target, plan) = fixture(false);
    let mut w = Workspace::from_documents(
        "pin",
        vec![
            plan.inputs[0].molecule.clone(),
            plan.inputs[1].molecule.clone(),
            target,
        ],
        2,
    )
    .unwrap();
    w.candidates.push(Candidate::new(plan));
    let bytes = workspace::save_workspace(&w).unwrap();
    for schema in 1..5 {
        let mut input = zip::ZipArchive::new(Cursor::new(&bytes)).unwrap();
        let mut output = zip::ZipWriter::new(Cursor::new(Vec::new()));
        for i in 0..input.len() {
            let mut entry = input.by_index(i).unwrap();
            let name = entry.name().to_owned();
            let mut data = Vec::new();
            entry.read_to_end(&mut data).unwrap();
            if name == "manifest.json" {
                let mut manifest: serde_json::Value = serde_json::from_slice(&data).unwrap();
                assert_eq!(
                    manifest["plan_schema"],
                    collagene_core::workspace::PLAN_SCHEMA_VERSION
                );
                manifest["plan_schema"] = schema.into();
                data = serde_json::to_vec(&manifest).unwrap();
            }
            output
                .start_file(name, zip::write::SimpleFileOptions::default())
                .unwrap();
            output.write_all(&data).unwrap();
        }
        assert!(workspace::load_workspace(&output.finish().unwrap().into_inner()).is_err());
    }
}

#[test]
fn forged_embedded_snapshot_is_rejected_before_plan_hashing() {
    let (target, mut plan) = fixture(false);
    if let Operation::RecombineReference { request, .. } = &mut plan.steps[0] {
        request.profile.profile.variants[0].sequence = "A".repeat(4097);
    }
    let r = verification::verify(&plan, &target);
    assert_eq!(r.status, Status::Invalid);
    assert!(r.products.is_empty());
    assert!(r.plan_revision.is_empty());
}

#[test]
fn text_roundtrips_keep_annotations_and_unverified_recipe_but_never_chemical_authority() {
    use collagene_core::{chemistry::Duplex, formats};
    let (target, plan) = fixture(false);
    let r = verification::verify(&plan, &target);
    let mut product = r.products[0].clone();
    product.name = "reference_product".into();
    assert!(product.original_genbank.is_none());
    assert!(plan.inputs[0].molecule.original_genbank.is_some());
    assert!(product.features.iter().any(|f| f.label() == "payload"));
    for format in [formats::TextFormat::GenBank, formats::TextFormat::Embl] {
        let out = formats::export_with_report(&product, format).unwrap();
        assert!(
            out.warnings
                .iter()
                .any(|w| w.contains("unverified") && w.contains("chemistry"))
        );
        let imported = formats::import_text(&out.text).unwrap().molecules.remove(0);
        assert!(imported.features.iter().any(|f| f.label() == "payload"));
        assert!(Duplex::from_molecule(&imported).is_err());
        assert_eq!(
            imported
                .metadata
                .get("source:recombination:reference_only")
                .map(String::as_str),
            Some("true")
        );
        assert_eq!(
            imported
                .metadata
                .get("source:recombination:unverified_recipe"),
            product.metadata.get(recombination::RECIPE_KEY)
        );
        assert!(!imported.metadata.contains_key(recombination::RECIPE_KEY));
        // The restrictive guard survives a second text export, with no native positive metadata.
        let twice =
            formats::import_text(&formats::export_with_report(&imported, format).unwrap().text)
                .unwrap()
                .molecules
                .remove(0);
        assert!(Duplex::from_molecule(&twice).is_err());
    }
    // A reference-only descendant can lack the original recipe; its restriction must survive.
    product.metadata.clear();
    product
        .metadata
        .insert(recombination::CHEMISTRY_UNMAPPED_KEY.into(), "true".into());
    let imported = formats::import_genbank(&formats::export_genbank(&product).unwrap())
        .unwrap()
        .molecules
        .remove(0);
    assert!(Duplex::from_molecule(&imported).is_err());
    let mut ordinary = Molecule::dna("ordinary", "ACGTACGT", Topology::Linear).unwrap();
    ordinary.features.push(collagene_model::Feature {
        kind: "misc_feature".into(),
        location: collagene_model::Location::range(0, 8),
        qualifiers: vec![(
            "note".into(),
            Some("Collagene reference-only v1: unverified sequence-only evidence".into()),
        )],
        needs_review: false,
    });
    let imported = formats::import_genbank(&formats::export_genbank(&ordinary).unwrap())
        .unwrap()
        .molecules
        .remove(0);
    assert!(Duplex::from_molecule(&imported).is_ok());
}

#[test]
fn native_source_and_target_trees_are_bounded_before_plan_or_target_hashing() {
    let (target, mut p) = fixture(false);
    let mut location = collagene_model::Location::range(0, 1);
    for _ in 0..128 {
        location = collagene_model::Location::Complement(Box::new(location));
    }
    let feature = collagene_model::Feature {
        kind: "misc_feature".into(),
        location,
        qualifiers: vec![],
        needs_review: false,
    };
    p.inputs[0].molecule.features.push(feature.clone());
    let r = verification::verify(&p, &target);
    assert_eq!(r.status, Status::Invalid);
    assert!(r.plan_revision.is_empty());
    let c = Candidate::new(p);
    assert!(c.plan_revision.is_empty());
    assert!(c.validate().is_err());
    let (_, p) = fixture(false);
    let mut t = target;
    t.features.push(feature);
    let r = verification::verify(&p, &t);
    assert!(r.plan_revision.is_empty());
    assert!(r.products.is_empty());
}

#[test]
fn imported_reference_guard_cannot_be_redeclared_as_materialized_annealing_input() {
    use collagene_core::{chemistry::AnnealOptions, formats};
    let (target, plan) = fixture(true);
    let mut reference = verification::verify(&plan, &target).products.remove(0);
    reference.name = "reference_source".into();
    for format in [formats::TextFormat::GenBank, formats::TextFormat::Embl] {
        let imported = formats::import_text(
            &formats::export_with_report(&reference, format)
                .unwrap()
                .text,
        )
        .unwrap()
        .molecules
        .remove(0);
        assert_eq!(
            imported
                .metadata
                .get("source:legacy:chemistry_unmapped")
                .map(String::as_str),
            Some("true")
        );
        let bottom = Molecule::dna(
            "independent complementary reference",
            imported.sequence.reverse_complement().unwrap().as_str(),
            Topology::Linear,
        )
        .unwrap();
        assert!(
            verification::simulate_oligo_annealing(&imported, &bottom, &AnnealOptions::default())
                .is_err()
        );
        let p = DesignPlan {
            name: "attempted import authority wash".into(),
            target_revision: target.revision(),
            inputs: vec![
                PlanInput {
                    id: "a".into(),
                    revision: imported.revision(),
                    molecule: imported,
                },
                PlanInput {
                    id: "b".into(),
                    revision: bottom.revision(),
                    molecule: bottom,
                },
            ],
            steps: vec![Operation::AnnealOligos {
                top: Source::Input("a".into()),
                bottom: Source::Input("b".into()),
                options: AnnealOptions::default(),
                candidate_index: 0,
            }],
        };
        let r = verification::verify(&p, &target);
        assert_eq!(r.status, Status::Incomplete);
        assert!(r.products.is_empty());
    }
    // Explicitly declared, ordinary ssDNA remains a supported independent input.
    let mut top = Molecule::dna("ordered top", "AACCGGTTAACCGGTT", Topology::Linear).unwrap();
    top.metadata
        .insert("genbank.molecule_type".into(), "ss-DNA".into());
    let bottom = Molecule::dna("ordered bottom", "AACCGGTTAACCGGTT", Topology::Linear).unwrap();
    assert!(
        verification::simulate_oligo_annealing(&top, &bottom, &AnnealOptions::default()).is_ok()
    );
}
