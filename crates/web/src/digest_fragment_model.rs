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
        );
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
    ) -> AssemblyFragment {
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
        let pin = PinnedOrigin::new(
            run_material::capsule(
                &output.workspace,
                &run_material::Selection {
                    run_id: output.receipt.run.run_id,
                    step_index: 0,
                    event: binding.clone(),
                },
            )
            .unwrap(),
        )
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
        AssemblyFragment {
            source: info,
            molecule: Arc::new(molecule),
            selection,
            options,
            fragment_index,
            batch_revision: inspected.context_revision,
            top_start: f.top_start,
            top_length: f.duplex.top_length,
            bottom_start: f.bottom_start,
            bottom_length: f.duplex.bottom_length,
            epoch: 4,
            origin: Some(pin),
        }
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
        assert_eq!(draft.version, 5);
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
        let restored = assembly_draft_model::decode_version(
            &assembly_draft_model::encode_raw(&draft).unwrap(),
            5,
        )
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
    fn digest_origin_staggered_free_end_handoff_keeps_forward_reference_without_inventing_reverse_bases()
     {
        use collagene_core::chemistry::{AnnealOptions, Terminus};
        let chosen = from_sources(
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
        let mut draft = Draft::default();
        assembly_draft_model::set_digest_fragment(&mut draft, chosen).unwrap();
        assert_eq!((draft.slots[0].start, draft.slots[0].end), (11, 21));
        assert!(draft.methods.rows[0].run_origin.is_some());
        draft.validate().unwrap();
        draft.slots[0].reverse = true;
        draft.methods.rows[0].reverse = true;
        assert!(
            draft.validate().is_err(),
            "Out-of-reference reverse intent must not be filled with guessed bases"
        );
    }
}
