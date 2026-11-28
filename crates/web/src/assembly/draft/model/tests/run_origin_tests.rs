use super::*;
use collagene_core::{
    documents::workspace::Workspace,
    material::{
        lineage::{EventBinding, InputEventBinding},
        outputs as cloning_outputs, run as run_material,
    },
    planning::plan::{DesignPlan, Operation, PlanInput, Source},
    sequence::analysis::Primer,
};

fn saved_pcr() -> Draft {
    saved_pcr_workspace().0
}
fn saved_pcr_workspace() -> (Draft, Workspace) {
    let source = Molecule::dna("Original", "ACGACCCCGGGGTAAC", Topology::Linear).unwrap();
    let target = Molecule::dna("Independent", "ACGACCCCGGGGTAAC", Topology::Linear).unwrap();
    let workspace =
        Workspace::from_documents("origin", vec![source.clone(), target.clone()], 0).unwrap();
    let event = |index: usize| {
        let doc = &workspace.documents[index];
        EventBinding {
            document_id: doc.id.clone(),
            event_id: doc.head_event.clone().unwrap(),
            revision: doc.head.clone(),
        }
    };
    let plan = DesignPlan {
        standard: None,
        name: "Original PCR".into(),
        target_revision: target.revision(),
        inputs: vec![PlanInput {
            id: "original".into(),
            revision: source.revision(),
            molecule: source,
        }],
        steps: vec![Operation::Pcr {
            template_strand: collagene_core::chemistry::Strand::Top,
            template: Source::Input("original".into()),
            forward: Primer::new("Original forward", "ACGACCCC", "").unwrap(),
            reverse: Primer::new("Original reverse", "GTTACCCC", "").unwrap(),
        }],
    };
    let commit = cloning_outputs::finalize(
        &workspace,
        &plan,
        event(1),
        vec![InputEventBinding {
            input_id: "original".into(),
            event: event(0),
        }],
        &cloning_outputs::Options {
            name: "Actual saved PCR".into(),
            save_primers: false,
            show_prepared: false,
            allow_incomplete: false,
        },
    )
    .unwrap();
    let selection = run_material::Selection {
        run_id: commit.receipt.run.run_id.clone(),
        step_index: 0,
        event: commit.receipt.construct.clone(),
    };
    let molecule = commit
        .workspace
        .head(&selection.event.document_id)
        .unwrap()
        .clone();
    let captured = Arc::new(CapturedSource {
        source: SourceInfo {
            binding: selection.event.clone(),
            name: molecule.name.clone(),
            length: molecule.sequence.len(),
            alphabet: Alphabet::Dna,
            topology: molecule.topology,
        },
        molecule,
    });
    let mut draft = Draft::default();
    draft.methods.mode = Mode::Ligation;
    draft.circular = false;
    draft.slots.truncate(1);
    {
        let collagene_core::cloning::assembly::intent::Selector::ReferenceRange {
            end: value, ..
        } = &mut draft.slots[0].selector
        else {
            panic!("fixture requires reference selector")
        };
        *value = captured.molecule.sequence.len();
    }
    draft.slots[0].captured = Some(captured);
    draft.methods.rows.truncate(1);
    draft.methods.rows[0].preparation = OverlapPreparation::WholeLinear;
    draft.methods.rows[0].run_origin = Some(
        crate::history::run_origin_model::PinnedOrigin::from_material(
            collagene_core::material::run::MaterialOrigin::Run(
                collagene_core::material::capsule::capsule(&commit.workspace, &selection).unwrap(),
            ),
        )
        .unwrap(),
    );
    draft.slots[0].material_origin = draft.methods.rows[0].run_origin.clone();
    draft.sync_rows();
    draft.validate().unwrap();
    (draft, commit.workspace)
}

#[test]
fn saved_origin_preserves_inactive_original_generation_text_without_redeclaration() {
    let mut draft = saved_pcr();
    draft.methods.rows[0].source_generation = Some(
        crate::cloning::preparation::cloning_methods_model::SourceGenerationInput {
            host: None,
            dam: "unfinished".into(),
            dcm: "fully".into(),
            cpg: "unknown".into(),
        },
    );
    draft.validate().unwrap();
    let restored = decode(&encode_raw(&draft).unwrap()).unwrap();
    assert!(restored == draft);
    let (_, rows) =
        crate::cloning::preparation::cloning_preparations::prepared_rows(&restored.methods, 1)
            .unwrap();
    assert!(rows[0].generation.is_none());
    assert!(matches!(
        rows[0].preparation,
        collagene_core::cloning::assembly::recipe::Preparation::FromRun { .. }
    ));
    let mut original = restored;
    original.methods.rows[0].run_origin = None;
    original.slots[0].material_origin = original.methods.rows[0].run_origin.clone();
    assert!(
        crate::cloning::preparation::cloning_preparations::prepared_rows(&original.methods, 1)
            .is_err()
    );
}

#[test]
fn derived_digest_roundtrip_retains_actual_preparation_and_rejects_superseded_encoding() {
    let (original, mut workspace) = saved_pcr_workspace();
    let collagene_core::material::run::MaterialOrigin::Run(run) = original.methods.rows[0]
        .run_origin
        .as_ref()
        .unwrap()
        .as_ref()
    else {
        panic!("Fixture requires actual saved run material");
    };
    let run = run.clone();
    let source = run.selection.event.clone();
    workspace.active_document = source.document_id.clone();
    let catalog = collagene_core::cloning::restriction::catalog::builtin_catalog().unwrap();
    let enzyme = catalog
        .edition()
        .enzymes
        .iter()
        .find(|e| e.name == "SmaI")
        .unwrap();
    let selection =
        collagene_core::cloning::restriction::catalog::select(&catalog, &[enzyme.id.clone()])
            .unwrap();
    let installed = collagene_core::material::derivation::prepare_install(
        &workspace,
        &collagene_core::material::derivation::Request {
            sources: vec![collagene_core::material::derivation::BoundSource {
                role: collagene_core::material::derivation::SourceRole::Template,
                event: source.clone(),
                origin: collagene_core::material::derivation::SourceOrigin::Saved(
                    collagene_core::material::run::MaterialOrigin::Run(run),
                ),
            }],
            derivation: collagene_core::material::derivation::Derivation::CatalogDigest {
                conditions: None,
                selection,
                options: Default::default(),
                fragment_index: 0,
            },
            name: "Saved SmaI fragment".into(),
            save_primers: false,
        },
    )
    .unwrap();
    assert_eq!(
        installed.receipt.sources,
        vec![collagene_core::material::derivation::SourceBinding {
            role: collagene_core::material::derivation::SourceRole::Template,
            event: source.clone()
        }]
    );
    assert_eq!(
        (
            installed.receipt.duplex.top_length,
            installed.receipt.duplex.bottom_length,
            installed.receipt.duplex.bottom_start
        ),
        (8, 8, 0)
    );
    let (left, right) = installed.receipt.duplex.ends.as_ref().unwrap();
    assert!(left.polarity == collagene_core::chemistry::Polarity::Blunt);
    assert!(right.polarity == collagene_core::chemistry::Polarity::Blunt);
    assert!(left.chemistry.five_prime == collagene_core::chemistry::FivePrime::Hydroxyl);
    assert!(right.chemistry.five_prime == collagene_core::chemistry::FivePrime::Phosphate);
    assert!(installed.receipt.duplex.nicks.is_empty());
    let event = installed.receipt.product.clone();
    let capsule =
        collagene_core::material::derivation::capsule(&installed.workspace, &event).unwrap();
    let before = serde_json::to_vec(&installed.workspace).unwrap();
    let record = collagene_core::material::derivation::load_record(&capsule.record_bytes).unwrap();
    assert_eq!(record.request.sources[0].event, source);
    assert!(matches!(
        record.request.sources[0].origin,
        collagene_core::material::derivation::SourceOrigin::Saved(
            collagene_core::material::run::MaterialOrigin::Run(_)
        )
    ));
    assert!(matches!(
        record.request.derivation,
        collagene_core::material::derivation::Derivation::CatalogDigest {
            fragment_index: 0,
            ..
        }
    ));
    let material = collagene_core::material::run::MaterialOrigin::Derived(capsule);
    let pin = crate::history::run_origin_model::PinnedOrigin::from_material(material).unwrap();
    let molecule = installed.workspace.revisions[&event.revision].clone();
    assert_eq!(molecule.sequence.as_str(), "ACGACCCC");
    let mut draft = original;
    draft.slots[0].captured = Some(Arc::new(CapturedSource {
        source: SourceInfo {
            binding: event,
            name: molecule.name.clone(),
            length: molecule.sequence.len(),
            alphabet: Alphabet::Dna,
            topology: molecule.topology,
        },
        molecule,
    }));
    {
        let collagene_core::cloning::assembly::intent::Selector::ReferenceRange {
            end: value, ..
        } = &mut draft.slots[0].selector
        else {
            panic!("fixture requires reference selector")
        };
        *value = 8;
    }
    draft.methods.rows[0].run_origin = Some(pin);
    draft.slots[0].material_origin = draft.methods.rows[0].run_origin.clone();
    draft.methods.rows[0].preparation = OverlapPreparation::ExactPcr;
    draft.methods.rows[0].forward = "ACGACCCC".into();
    draft.methods.rows[0].reverse_primer = "GGGGTCGT".into();
    draft.sync_rows();
    assert_eq!(draft.version, CURRENT_VERSION);
    let bytes = encode_raw(&draft).unwrap();
    let restored = decode(&bytes).unwrap();
    assert!(restored == draft);
    let mut obsolete = serde_json::to_value(&draft).unwrap();
    let current = obsolete["methods"]["rows"][0]["run_origin"]
        .as_object_mut()
        .unwrap();
    let body = current.remove("Derived").unwrap();
    current.insert("Digest".into(), body);
    assert!(
        decode(&serde_json::to_vec(&obsolete).unwrap()).is_err(),
        "Removed origin tag must never be reinterpreted"
    );

    let (_, rows) =
        crate::cloning::preparation::cloning_preparations::prepared_rows(&restored.methods, 1)
            .unwrap();
    assert!(matches!(
        &rows[0].preparation,
        collagene_core::cloning::assembly::recipe::Preparation::FromDerived {
            preparation: collagene_core::cloning::assembly::recipe::RunPreparation::ExactPcr { .. },
            ..
        }
    ));
    let mut downgraded = draft.clone();
    downgraded.version = CURRENT_VERSION - 1;
    assert!(downgraded.validate().is_err());
    let mut future = draft;
    future.version = CURRENT_VERSION + 1;
    assert!(future.validate().is_err());
    assert_eq!(serde_json::to_vec(&installed.workspace).unwrap(), before);
}

#[test]
fn reusable_pcr_derived_origin_keeps_native_record_and_primers() {
    let (mut draft, mut workspace) = saved_pcr_workspace();
    let source_event = {
        let document = &workspace.documents[0];
        EventBinding {
            document_id: document.id.clone(),
            event_id: document.head_event.clone().unwrap(),
            revision: document.head.clone(),
        }
    };
    workspace.active_document = source_event.document_id.clone();
    let source = &workspace.revisions[&source_event.revision];
    let forward = Primer::new("F", "ACGACCCC", "").unwrap();
    let reverse = Primer::new("R", "GTTACCCC", "").unwrap();
    let preview = collagene_core::verification::pcr::simulate_pcr_documents(
        source,
        &forward,
        &reverse,
        collagene_core::chemistry::Strand::Top,
    )
    .unwrap()
    .remove(0);
    let request = collagene_core::material::derivation::Request {
        sources: vec![collagene_core::material::derivation::BoundSource {
            role: collagene_core::material::derivation::SourceRole::Template,
            event: source_event,
            origin: collagene_core::material::derivation::SourceOrigin::OriginalReference,
        }],
        derivation: collagene_core::material::derivation::Derivation::ExactPcr {
            template_strand: collagene_core::chemistry::Strand::Top,
            forward,
            reverse,
            selection: collagene_core::material::derivation::PcrSelection {
                selected_index: 0,
                expected_preview_revision: preview.molecule.revision(),
                geometry: preview.geometry,
            },
        },
        name: "Reusable exact PCR".into(),
        save_primers: true,
    };
    let installed =
        collagene_core::material::derivation::prepare_install(&workspace, &request).unwrap();
    assert!(installed.receipt.primer_attachment.is_some());
    let product = installed.receipt.product.clone();
    let molecule = installed.workspace.revisions[&product.revision].clone();
    let capsule =
        collagene_core::material::derivation::capsule(&installed.workspace, &product).unwrap();
    draft.slots[0].captured = Some(Arc::new(CapturedSource {
        source: SourceInfo {
            binding: product,
            name: molecule.name.clone(),
            length: molecule.sequence.len(),
            alphabet: Alphabet::Dna,
            topology: molecule.topology,
        },
        molecule,
    }));
    draft.methods.rows[0].run_origin = Some(
        crate::history::run_origin_model::PinnedOrigin::from_material(
            run_material::MaterialOrigin::Derived(capsule),
        )
        .unwrap(),
    );
    draft.slots[0].material_origin = draft.methods.rows[0].run_origin.clone();
    draft.sync_rows();
    assert_eq!(draft.version, CURRENT_VERSION);
    let bytes = encode_raw(&draft).unwrap();
    assert!(decode(&bytes).unwrap() == draft);

    let (_, rows) =
        crate::cloning::preparation::cloning_preparations::prepared_rows(&draft.methods, 1)
            .unwrap();
    assert!(matches!(
        &rows[0].preparation,
        collagene_core::cloning::assembly::recipe::Preparation::FromDerived { .. }
    ));
    draft.version = CURRENT_VERSION - 1;
    assert!(draft.validate().is_err());
    draft.version = CURRENT_VERSION + 1;
    assert!(draft.validate().is_err());
}

#[test]
fn saved_origin_mismatch_preparation_keeps_original_event_and_policy() {
    let mut draft = saved_pcr();
    let row = &mut draft.methods.rows[0];
    row.preparation = OverlapPreparation::MismatchPcr;
    row.forward = "ACGACCCC".into();
    row.reverse_primer = "GTTACCCC".into();
    row.mismatch_max_mismatches = "2".into();
    row.mismatch_three_prime_exact = "3".into();
    let restored = decode(&encode_raw(&draft).unwrap()).unwrap();
    let (_, rows) =
        crate::cloning::preparation::cloning_preparations::prepared_rows(&restored.methods, 1)
            .unwrap();
    let collagene_core::cloning::assembly::recipe::Preparation::FromRun {
        origin,
        preparation:
            collagene_core::cloning::assembly::recipe::RunPreparation::MismatchPcr { options, .. },
    } = &rows[0].preparation
    else {
        panic!("saved mismatch preparation")
    };
    assert_eq!(
        origin.selection.event,
        restored.slots[0].captured.as_ref().unwrap().source.binding
    );
    assert_eq!(options.binding.max_mismatches, 2);
    assert_eq!(options.binding.min_three_prime_exact, 3);
}

#[test]
fn exact_saved_source_survives_v5_roundtrip_and_compiles_explicit_reconstruction() {
    let mut draft = saved_pcr();
    assert_eq!(draft.version, CURRENT_VERSION);
    let pin = draft.methods.rows[0].run_origin.clone();
    draft.sync_rows();
    assert!(draft.methods.rows[0].run_origin == pin);
    let restored = decode(&encode_raw(&draft).unwrap()).unwrap();
    assert!(restored == draft);
    let (_, rows) =
        crate::cloning::preparation::cloning_preparations::prepared_rows(&restored.methods, 1)
            .unwrap();
    let collagene_core::cloning::assembly::recipe::Preparation::FromRun {
        origin,
        preparation,
    } = &rows[0].preparation
    else {
        panic!("saved material must remain explicit")
    };
    assert_eq!(
        origin.selection.event,
        restored.slots[0].captured.as_ref().unwrap().source.binding
    );
    assert!(matches!(
        preparation,
        collagene_core::cloning::assembly::recipe::RunPreparation::WholeLinear
    ));
    let mut old = draft.clone();
    old.version = CURRENT_VERSION - 1;
    assert!(old.validate().is_err());
    draft.methods.rows[0].preparation = OverlapPreparation::DeclaredDuplex;
    assert!(
        crate::cloning::preparation::cloning_preparations::prepared_rows(&draft.methods, 1)
            .is_err()
    );
}

#[test]
fn mismatched_persisted_origin_is_refused_while_explicit_source_edit_clears_it() {
    let mut draft = saved_pcr();
    {
        let collagene_core::cloning::assembly::intent::Selector::ReferenceRange {
            start: value,
            ..
        } = &mut draft.slots[0].selector
        else {
            panic!("fixture requires reference selector")
        };
        *value = 2;
    }
    {
        let collagene_core::cloning::assembly::intent::Selector::ReferenceRange {
            end: value, ..
        } = &mut draft.slots[0].selector
        else {
            panic!("fixture requires reference selector")
        };
        *value = 14;
    }
    draft.methods.rows[0].treatment =
        crate::cloning::preparation::cloning_methods_model::Treatment::Phosphorylate;
    draft.sync_rows();
    let mut changed = draft.clone();
    Arc::make_mut(changed.slots[0].captured.as_mut().unwrap())
        .source
        .binding
        .document_id = "unrelated-source".into();
    assert!(decode(&serde_json::to_vec(&changed).unwrap()).is_err());
    changed.sync_rows();
    assert!(changed.methods.rows[0].run_origin.is_none());
    assert_eq!(
        (
            changed.slots[0].selector.reference_interval().unwrap().0,
            changed.slots[0].selector.reference_interval().unwrap().1
        ),
        (2, 14)
    );
    assert!(changed.methods.rows[0].treatment == draft.methods.rows[0].treatment);
    assert_eq!(changed.name, draft.name);
    assert_eq!(changed.circular, draft.circular);
    assert!(changed.junction_edits == draft.junction_edits);
}
