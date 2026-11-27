//! Source-bound digest handoff data, shared by the view and native draft tests.
use crate::{run_origin_model::PinnedOrigin, source_identity::SourceInfo};
use collagene_core::{
    Molecule,
    catalog_digest::{CatalogSelection, InspectionOptions},
};
use std::sync::Arc;

/// An explicit fragment choice retains its original source and reconstruction
/// recipe. The displayed geometry is not serialized or installed as live stock.
#[derive(Clone)]
pub struct AssemblyFragment {
    pub source: SourceInfo,
    pub molecule: Arc<Molecule>,
    pub selection: CatalogSelection,
    pub options: InspectionOptions,
    pub fragment_index: usize,
    /// Ordinary catalog batch, or full origin-bound inspection context.
    pub batch_revision: String,
    pub top_start: usize,
    pub top_length: usize,
    pub bottom_start: i64,
    pub bottom_length: usize,
    pub epoch: u64,
    pub origin: Option<PinnedOrigin>,
    pub material_domain: Option<crate::assembly_draft_model::MaterialDomain>,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        assembly_draft_model::{self, Draft},
        cloning_methods_model::OverlapPreparation,
    };
    use collagene_core::{
        Topology, cloning_outputs,
        history::DocumentHistory,
        methylation::MethylationState,
        modified_duplex::GenerationAssumptions,
        run_lineage::{EventBinding, InputEventBinding},
        run_material,
        verification::{DesignPlan, Operation, PlanInput, Source},
        workspace::Workspace,
    };
    fn fixture() -> AssemblyFragment {
        let chosen = from_sources(
            vec![Molecule::dna("Original circle", "AAAGAATTCTTT", Topology::Circular).unwrap()],
            Molecule::dna("Independent target", "AAAGAATTCTTT", Topology::Circular).unwrap(),
            vec![Operation::DeclareSourceModifications {
                template: Source::Input("source0".into()),
                assumptions: GenerationAssumptions {
                    dam: MethylationState::Unmodified,
                    dcm: MethylationState::Unmodified,
                    cpg: MethylationState::Unmodified,
                },
            }],
            0,
        )
        .0;
        assert_eq!(
            (chosen.top_start, chosen.top_length, chosen.bottom_start),
            (4, 12, 8)
        );
        chosen
    }
    fn from_sources(
        originals: Vec<Molecule>,
        target: Molecule,
        steps: Vec<Operation>,
        fragment_index: usize,
    ) -> (AssemblyFragment, Workspace) {
        let target_index = originals.len();
        let plan = DesignPlan {
            name: "Explicit saved source".into(),
            target_revision: target.revision(),
            inputs: originals
                .iter()
                .enumerate()
                .map(|(i, m)| PlanInput {
                    id: format!("source{i}"),
                    revision: m.revision(),
                    molecule: m.clone(),
                })
                .collect(),
            steps,
        };
        let mut docs = originals;
        docs.push(target);
        let w = Workspace::from_documents("Saved material", docs, target_index).unwrap();
        let event = |index: usize| {
            let d = &w.documents[index];
            EventBinding {
                document_id: d.id.clone(),
                event_id: d.head_event.clone().unwrap(),
                revision: d.head.clone(),
            }
        };
        let output = cloning_outputs::finalize(
            &w,
            &plan,
            event(target_index),
            (0..target_index)
                .map(|i| InputEventBinding {
                    input_id: format!("source{i}"),
                    event: event(i),
                })
                .collect(),
            &cloning_outputs::Options {
                name: "Saved source".into(),
                ..Default::default()
            },
        )
        .unwrap();
        let binding = output.receipt.construct.clone();
        let molecule = output.workspace.head(&binding.document_id).unwrap().clone();
        let pin = PinnedOrigin::from_material(collagene_core::run_material::MaterialOrigin::Run(
            run_material::capsule(
                &output.workspace,
                &run_material::Selection {
                    run_id: output.receipt.run.run_id,
                    step_index: 0,
                    event: binding.clone(),
                },
            )
            .unwrap(),
        ))
        .unwrap();
        let info = SourceInfo {
            binding,
            name: molecule.name.clone(),
            length: molecule.sequence.len(),
            alphabet: molecule.sequence.alphabet(),
            topology: molecule.topology,
        };
        let catalog = collagene_core::catalog_digest::builtin_catalog().unwrap();
        let enzyme = catalog
            .edition()
            .enzymes
            .iter()
            .find(|e| e.name == "EcoRI")
            .unwrap();
        let selection =
            collagene_core::catalog_digest::select(&catalog, &[enzyme.id.clone()]).unwrap();
        let options = InspectionOptions::default();
        let role = collagene_core::restriction_history::source_role(&info.binding).unwrap();
        let inspected = collagene_core::material_derivation::inspect_digest(
            &role,
            &molecule,
            &collagene_core::material_derivation::SourceOrigin::Saved(pin.as_ref().clone()),
            &selection,
            &options,
        )
        .unwrap();
        let f = &inspected.inspection.fragments[fragment_index];
        (
            AssemblyFragment {
                source: info,
                molecule: Arc::new(molecule),
                selection,
                options,
                fragment_index,
                material_domain: Some(crate::assembly_draft_model::MaterialDomain {
                    top_length: inspected.source_duplex.top_length,
                    bottom_start: inspected.source_duplex.bottom_start,
                    bottom_length: inspected.source_duplex.bottom_length,
                    topology: inspected.source_duplex.topology,
                    state_revision: inspected.source_modifications.state_revision.clone(),
                }),
                batch_revision: inspected.context_revision,
                top_start: f.top_start,
                top_length: f.duplex.top_length,
                bottom_start: f.bottom_start,
                bottom_length: f.duplex.bottom_length,
                epoch: 4,
                origin: Some(pin),
            },
            output.workspace,
        )
    }
    #[test]
    fn digest_origin_handoff_keeps_saved_parent_and_reconstructable_row_after_roundtrip() {
        let chosen = fixture();
        let parent = chosen.source.clone();
        let original = chosen.molecule.clone();
        let pin = chosen.origin.clone();
        let mut draft = Draft::default();
        draft.name = "My intended construct".into();
        // This fixture exercises one saved preparation, not an unspecified
        // second source with the default empty exact-PCR primers.
        draft.slots.truncate(1);
        draft.sync_rows();
        assembly_draft_model::set_digest_fragment(&mut draft, chosen).unwrap();
        assert_eq!(draft.name, "My intended construct");
        assert_eq!(draft.version, crate::assembly_draft_model::CURRENT_VERSION);
        assert!(draft.methods.rows[0].run_origin == pin);
        assert_eq!(draft.slots[0].captured.as_ref().unwrap().source, parent);
        assert_eq!(
            draft.slots[0].captured.as_ref().unwrap().molecule,
            *original
        );
        assert_eq!(
            draft.slots[0]
                .captured
                .as_ref()
                .unwrap()
                .molecule
                .sequence
                .as_str(),
            "AAAGAATTCTTT"
        );
        assert!(draft.methods.rows[0].preparation == OverlapPreparation::CatalogDigest);
        let restored =
            assembly_draft_model::decode(&assembly_draft_model::encode_raw(&draft).unwrap())
                .unwrap();
        assert!(restored.methods.rows[0].run_origin == pin);
        let (_, rows) = crate::cloning_preparations::prepared_rows(&restored.methods, 1).unwrap();
        assert!(matches!(
            rows[0].preparation,
            collagene_core::assembly_recipe::Preparation::FromRun {
                preparation: collagene_core::assembly_recipe::RunPreparation::CatalogDigest { .. },
                ..
            }
        ));
    }
    #[test]
    fn digest_origin_handoff_rejects_changed_scope_or_equal_content_alias_atomically() {
        let chosen = fixture();
        let initial = Draft::default();
        let mut stale = chosen.clone();
        stale.batch_revision = "0".repeat(64);
        let mut draft = initial.clone();
        assert!(assembly_draft_model::set_digest_fragment(&mut draft, stale).is_err());
        assert!(draft == initial);
        let alias = DocumentHistory::new("equal-reference-alias", chosen.molecule.as_ref().clone())
            .unwrap();
        let mut alias_choice = chosen;
        alias_choice.source.binding = EventBinding {
            document_id: alias.id.clone(),
            event_id: alias.head_event().into(),
            revision: alias.head_revision().into(),
        };
        assert!(assembly_draft_model::set_digest_fragment(&mut draft, alias_choice).is_err());
        assert!(draft == initial);
    }
    #[test]
    fn saved_staggered_actual_bottom_intent_reverses_twice_and_ligates_after_archive_restart() {
        use collagene_core::chemistry::{AnnealOptions, Terminus};
        let (chosen, frame) = from_sources(
            vec![
                Molecule::dna("Top", "AATTGCTACGGAATTCTTACG", Topology::Linear).unwrap(),
                Molecule::dna("Bottom", "AATTCGTAAGAATTCCGTAGC", Topology::Linear).unwrap(),
            ],
            Molecule::dna("Independent top", "AATTGCTACGGAATTCTTACG", Topology::Linear).unwrap(),
            vec![Operation::AnnealOligos {
                top: Source::Input("source0".into()),
                bottom: Source::Input("source1".into()),
                options: AnnealOptions {
                    min_paired_bases: 17,
                    top_chemistry: Terminus::CLEAVED,
                    bottom_chemistry: Terminus::CLEAVED,
                    ..Default::default()
                },
                candidate_index: 0,
            }],
            1,
        );
        assert_eq!(
            (
                chosen.top_start,
                chosen.top_length,
                chosen.bottom_start,
                chosen.bottom_length
            ),
            (11, 10, 15, 10)
        );
        let derived = collagene_core::material_derivation::install(
            &frame,
            &collagene_core::material_derivation::Request {
                source: chosen.source.binding.clone(),
                origin: collagene_core::material_derivation::SourceOrigin::Saved(
                    chosen.origin.as_ref().unwrap().as_ref().clone(),
                ),
                derivation: collagene_core::material_derivation::Derivation::CatalogDigest {
                    selection: chosen.selection.clone(),
                    options: chosen.options.clone(),
                    fragment_index: chosen.fragment_index,
                },
                name: "Selected sticky fragment".into(),
                save_primers: false,
            },
        )
        .unwrap();
        let derived_frame = collagene_core::workspace::load_workspace(&derived.archive).unwrap();
        let derived_molecule = derived_frame.revisions[&derived.receipt.product.revision].clone();
        let derived_origin = collagene_core::run_material::MaterialOrigin::Derived(
            collagene_core::material_derivation::capsule(&derived_frame, &derived.receipt.product)
                .unwrap(),
        );
        let derived_role =
            collagene_core::restriction_history::source_role(&derived.receipt.product).unwrap();
        let derived_preview = collagene_core::material_derivation::inspect_digest(
            &derived_role,
            &derived_molecule,
            &collagene_core::material_derivation::SourceOrigin::Saved(derived_origin.clone()),
            &chosen.selection,
            &chosen.options,
        )
        .unwrap();
        let f = &derived_preview.inspection.fragments[0];
        let derived_choice = AssemblyFragment {
            source: SourceInfo {
                binding: derived.receipt.product.clone(),
                name: derived_molecule.name.clone(),
                length: derived_molecule.sequence.len(),
                alphabet: derived_molecule.sequence.alphabet(),
                topology: derived_molecule.topology,
            },
            molecule: Arc::new(derived_molecule),
            selection: chosen.selection.clone(),
            options: chosen.options.clone(),
            fragment_index: 0,
            batch_revision: derived_preview.context_revision.clone(),
            top_start: f.top_start,
            top_length: f.duplex.top_length,
            bottom_start: f.bottom_start,
            bottom_length: f.duplex.bottom_length,
            epoch: 4,
            origin: Some(PinnedOrigin::from_material(derived_origin).unwrap()),
            material_domain: Some(assembly_draft_model::MaterialDomain {
                top_length: derived_preview.source_duplex.top_length,
                bottom_start: derived_preview.source_duplex.bottom_start,
                bottom_length: derived_preview.source_duplex.bottom_length,
                topology: derived_preview.source_duplex.topology,
                state_revision: derived_preview.source_modifications.state_revision.clone(),
            }),
        };
        let mut derived_draft = Draft::default();
        derived_draft.slots.truncate(1);
        derived_draft.methods.mode = crate::cloning_methods_model::Mode::Ligation;
        derived_draft.sync_rows();
        assembly_draft_model::set_digest_fragment(&mut derived_draft, derived_choice).unwrap();
        let derived_draft = assembly_draft_model::orient_fragment(&derived_draft, 0, true).unwrap();
        let (derived_inputs, derived_request, _) =
            assembly_draft_model::inputs_request(&derived_draft, Some(&derived_frame)).unwrap();
        let derived_intent = collagene_core::assembly_intent::build(
            &derived_inputs,
            &derived_request,
            Some(&derived_frame),
        )
        .unwrap();
        assert_eq!(derived_intent.target.sequence.as_str(), "AATTCGTAAG");
        let (_, derived_rows) =
            crate::cloning_preparations::prepared_rows(&derived_draft.methods, 1).unwrap();
        assert!(matches!(
            &derived_rows[0].preparation,
            collagene_core::assembly_recipe::Preparation::FromDerived { .. }
        ));
        let derived_recipe = collagene_core::ligation_recipe::Recipe {
            name: "Derived actual strand circle".into(),
            fragments: derived_rows
                .into_iter()
                .map(|source| collagene_core::ligation_recipe::Fragment {
                    source,
                    treatment: None,
                })
                .collect(),
            circular: true,
            allow_nicked: false,
        };
        let derived_checked = collagene_core::ligation_recipe::preview(
            &derived_inputs,
            &derived_intent.target,
            &derived_recipe,
            collagene_core::assembly_recipe::CleavagePolicy::AssumedComplete,
        )
        .unwrap();
        assert!(
            derived_checked.verification.status
                == collagene_core::verification::Status::ConsistentUnderAssumptions
        );

        let source_binding = chosen.source.binding.clone();
        let original_documents = frame.documents.clone();
        let mut draft = Draft::default();
        draft.name = "Independent actual bottom contribution".into();
        draft.methods.mode = crate::cloning_methods_model::Mode::Ligation;
        draft.slots.truncate(1);
        draft.sync_rows();
        assembly_draft_model::set_digest_fragment(&mut draft, chosen).unwrap();
        let forward = draft.slots[0].selector.clone();
        assert!(matches!(
            &forward,
            collagene_core::assembly_intent::Selector::MaterialStrandRange {
                strand: collagene_core::chemistry::Strand::Top,
                start: 11,
                end: 21,
                ..
            }
        ));
        let reversed = assembly_draft_model::orient_fragment(&draft, 0, true).unwrap();
        assert!(matches!(
            &reversed.slots[0].selector,
            collagene_core::assembly_intent::Selector::MaterialStrandRange {
                strand: collagene_core::chemistry::Strand::Bottom,
                start: 0,
                end: 10,
                ..
            }
        ));
        assert!(
            assembly_draft_model::orient_fragment(&reversed, 0, false)
                .unwrap()
                .slots[0]
                .selector
                == forward
        );
        draft = reversed;
        let (inputs, request, _) =
            assembly_draft_model::inputs_request(&draft, Some(&frame)).unwrap();
        assert!(collagene_core::assembly_intent::build(&inputs, &request, None).is_err());
        let intent =
            collagene_core::assembly_intent::build(&inputs, &request, Some(&frame)).unwrap();
        assert_eq!(intent.target.sequence.as_str(), "AATTCGTAAG");
        assert_ne!(intent.target.sequence.as_str(), "CGTAAGAATT");
        let source_review =
            crate::assembly_edit_model::SourceReview::capture(&draft, &intent.fragments).unwrap();
        assert_eq!(
            crate::assembly_edit_model::source_text(&draft, 0, 0, 10, Some(&source_review))
                .unwrap(),
            "AATTCGTAAG"
        );
        assert!(
            crate::assembly_edit_model::top_regions(&draft, 0, 0, 4, Some(&source_review))
                .unwrap()
                .is_empty()
        );
        assert_eq!(
            crate::assembly_edit_model::top_regions(&draft, 0, 4, 10, Some(&source_review))
                .unwrap(),
            vec![(15, 21, true)]
        );
        let edit =
            crate::assembly_edit_model::declare(&draft, 0, "1", "2", "C", Some(&source_review))
                .unwrap();
        assert_eq!(edit.edits[0].expected, "A");
        assert_eq!(
            crate::assembly_edit_model::source_text(&edit, 0, 0, 4, Some(&source_review)).unwrap(),
            "AATT",
            "Original receipt must not read edited target DNA back"
        );
        let mut state_changed = draft.clone();
        if let collagene_core::assembly_intent::Selector::MaterialStrandRange {
            expected_state_revision,
            ..
        } = &mut state_changed.slots[0].selector
        {
            *expected_state_revision = "0".repeat(64);
        }
        assert!(source_review.get(&state_changed, 0).is_err());
        assert!(assembly_draft_model::inputs_request(&state_changed, Some(&frame)).is_err());
        let (_, prepared) = crate::cloning_preparations::prepared_rows(&draft.methods, 1).unwrap();
        let recipe = collagene_core::ligation_recipe::Recipe {
            name: "Reverse actual sticky circle".into(),
            fragments: prepared
                .into_iter()
                .map(|source| collagene_core::ligation_recipe::Fragment {
                    source,
                    treatment: None,
                })
                .collect(),
            circular: true,
            allow_nicked: false,
        };
        let checked = collagene_core::ligation_recipe::preview(
            &inputs,
            &intent.target,
            &recipe,
            collagene_core::assembly_recipe::CleavagePolicy::AssumedComplete,
        )
        .unwrap();
        assert!(
            checked.verification.status
                == collagene_core::verification::Status::ConsistentUnderAssumptions
        );
        let staged = crate::assembly_review::stage(
            &frame,
            &intent,
            &[source_binding.clone()],
            &collagene_core::transfer::Candidate::new(checked.compiled.plan.clone()),
        )
        .unwrap();
        let saved = cloning_outputs::finalize(
            &staged.workspace,
            &checked.compiled.plan,
            staged.target,
            staged.sources,
            &cloning_outputs::Options {
                name: "Actual strand ligation".into(),
                ..Default::default()
            },
        )
        .unwrap();
        let draft_bytes = assembly_draft_model::encode_raw(&draft).unwrap();
        let mut with_draft = saved.workspace;
        let asset = collagene_core::workspace::WorkspaceAttachment {
            kind: assembly_draft_model::ATTACHMENT_KIND.into(),
            media_type: "application/json".into(),
            name: "Actual strand draft".into(),
            bytes: draft_bytes,
        };
        let hash = blake3::hash(&asset.bytes).to_hex().to_string();
        with_draft.attachments.insert(hash.clone(), asset);
        let reopened = collagene_core::workspace::load_workspace(
            &collagene_core::workspace::save_workspace(&with_draft).unwrap(),
        )
        .unwrap();
        let reopened_draft =
            assembly_draft_model::decode(&reopened.attachments[&hash].bytes).unwrap();
        let (again, declared, _) =
            assembly_draft_model::inputs_request(&reopened_draft, Some(&reopened)).unwrap();
        let replayed_intent =
            collagene_core::assembly_intent::build(&again, &declared, Some(&reopened)).unwrap();
        assert_eq!(replayed_intent.target.sequence.as_str(), "AATTCGTAAG");
        for original in &original_documents {
            let retained = reopened
                .documents
                .iter()
                .find(|doc| doc.id == original.id)
                .unwrap();
            assert_eq!(retained.head, original.head);
            assert_eq!(retained.head_event, original.head_event);
        }
        let wrong =
            Molecule::dna("Wrong independent target", "AATTCGTAAA", Topology::Circular).unwrap();
        let rejected = collagene_core::ligation_recipe::preview(
            &inputs,
            &wrong,
            &recipe,
            collagene_core::assembly_recipe::CleavagePolicy::AssumedComplete,
        )
        .unwrap();
        assert!(
            rejected.verification.status
                != collagene_core::verification::Status::ConsistentUnderAssumptions
        );
    }
}
