//! Full-provenance library fixtures use independent source documents, explicit
//! annealing choices and independently stated target DNA, never saved verdicts.
use super::*;
use crate::assembly::draft::model::CapturedSource;
use crate::assembly::library::{
    assembly_library_topology as topology, assembly_saved_designs_model as saved,
};
use crate::sources::identity::SourceInfo;
use collagene_core::{
    chemistry::{AnnealOptions, Terminus},
    documents::{history::DocumentHistory, workspace::Workspace},
    material::{lineage::InputEventBinding, outputs},
    planning::{combinatorial::release, plan::Operation},
    verification,
};
use collagene_model::{Molecule, Topology};
use std::sync::Arc;

fn captured(frame: &Workspace, index: usize) -> Arc<CapturedSource> {
    let document = &frame.documents[index];
    let molecule = frame.revisions[&document.head].clone();
    Arc::new(CapturedSource {
        source: SourceInfo {
            binding: EventBinding {
                document_id: document.id.clone(),
                event_id: document.head_event.clone().unwrap(),
                revision: molecule.revision(),
            },
            name: molecule.name.clone(),
            length: molecule.sequence.len(),
            alphabet: molecule.sequence.alphabet(),
            topology: molecule.topology,
        },
        molecule,
    })
}
fn fixture() -> (Workspace, Draft) {
    let frame = Workspace::from_documents(
        "Ordered oligo library",
        vec![
            Molecule::dna("First top", "ACGTTA", Topology::Linear).unwrap(),
            Molecule::dna("Second top", "AAGGCC", Topology::Linear).unwrap(),
            Molecule::dna("First bottom", "TAACGT", Topology::Linear).unwrap(),
            Molecule::dna("Second bottom", "GGCCTT", Topology::Linear).unwrap(),
        ],
        0,
    )
    .unwrap();
    let mut draft = Draft::default();
    draft.name = "Independently declared oligo concatenation".into();
    draft.circular = false;
    draft.methods.mode = Mode::Ligation;
    for index in 0..2 {
        draft.slots[index].captured = Some(captured(&frame, index));
        draft.slots[index].anneal_bottom = Some(captured(&frame, index + 2));
        draft.slots[index].selector =
            collagene_core::cloning::assembly::intent::Selector::ReferenceRange {
                start: 0,
                end: 6,
                reverse: false,
            };
        draft.methods.rows[index].preparation = OverlapPreparation::AnnealedOligos;
        draft.methods.rows[index].anneal_options = AnnealOptions {
            min_paired_bases: 6,
            top_chemistry: Terminus::CLEAVED,
            bottom_chemistry: Terminus::CLEAVED,
            ..Default::default()
        };
    }
    draft.sync_rows();
    for index in 0..2 {
        let top = &draft.slots[index]
            .captured
            .as_ref()
            .unwrap()
            .molecule
            .sequence;
        let bottom = &draft.slots[index]
            .anneal_bottom
            .as_ref()
            .unwrap()
            .molecule
            .sequence;
        let inspected = collagene_core::chemistry::anneal_oligos(
            top,
            bottom,
            &draft.methods.rows[index].anneal_options,
        )
        .unwrap();
        assert!(inspected.complete);
        assert_eq!(inspected.candidates.len(), 1);
        draft.methods.rows[index].anneal_candidate = Some(0);
    }
    draft.validate().unwrap();
    (frame, draft)
}
fn library() -> (Workspace, Draft) {
    let (frame, draft) = fixture();
    (
        frame,
        alternatives::enable(&draft, "paired-library").unwrap(),
    )
}
fn row(request: &batch::Request, index: usize) -> &SelectedConstruct {
    let RowRequest::Ready(row) = &request.rows[index] else {
        panic!("expected ready source-bound row")
    };
    row
}
fn expected_bindings(draft: &Draft) -> Vec<InputEventBinding> {
    let ids = draft.methods.all_input_ids().unwrap();
    draft
        .slots
        .iter()
        .filter_map(|slot| slot.captured.as_ref())
        .chain(
            draft
                .slots
                .iter()
                .filter_map(|slot| slot.anneal_bottom.as_ref()),
        )
        .zip(ids)
        .map(|(source, input_id)| InputEventBinding {
            input_id,
            event: source.source.binding.clone(),
        })
        .collect()
}

#[test]
fn paired_library_projection_compiles_fresh_annealing_with_both_ordered_sources() {
    let (frame, ordinary) = fixture();
    let contexts = ordinary
        .methods
        .rows
        .iter()
        .map(|raw| raw.anneal_source_context.clone())
        .collect::<Vec<_>>();
    let draft = alternatives::enable(&ordinary, "paired-library").unwrap();
    assert_eq!(
        draft
            .methods
            .rows
            .iter()
            .map(|raw| raw.anneal_candidate)
            .collect::<Vec<_>>(),
        vec![Some(0), Some(0)]
    );
    assert_eq!(
        draft
            .methods
            .rows
            .iter()
            .map(|raw| raw.anneal_source_context.clone())
            .collect::<Vec<_>>(),
        contexts
    );
    let projection = project(&draft).unwrap();
    let design = projection.request.design.as_ref().unwrap();
    let selected = row(&projection.request, 0);
    assert_eq!(
        design.bins[0].choices[0].anneal_bottom,
        Some(captured(&frame, 2).source.binding.clone())
    );
    assert_eq!(
        design.bins[1].choices[0].anneal_bottom,
        Some(captured(&frame, 3).source.binding.clone())
    );
    let expected = expected_bindings(&draft);
    assert_eq!(design.selected_input_bindings(selected).unwrap(), expected);
    let compiled = design::compile(design, selected, &frame).unwrap();
    assert_eq!(compiled.resolved.bindings, expected);
    assert_eq!(compiled.resolved.inputs.len(), 4);
    assert_eq!(compiled.resolved.target.sequence.as_str(), "ACGTTAAAGGCC");
    let anneals = compiled
        .plan
        .steps
        .iter()
        .filter_map(|step| match step {
            Operation::AnnealOligos {
                top,
                bottom,
                candidate_index,
                ..
            } => Some((top, bottom, *candidate_index)),
            _ => None,
        })
        .collect::<Vec<_>>();
    assert_eq!(anneals.len(), 2);
    assert!(anneals.iter().all(|(_, _, candidate)| *candidate == 0));
    let report = verification::engine::verify(&compiled.plan, &compiled.resolved.target);
    assert_eq!(
        report.final_product().unwrap().sequence.as_str(),
        "ACGTTAAAGGCC"
    );
    assert!(
        draft
            .alternatives
            .as_ref()
            .unwrap()
            .bins
            .iter()
            .flat_map(|bin| &bin.choices)
            .all(|choice| choice.row.anneal_bottom_document.is_none())
    );
}

#[test]
fn construct_namespace_and_bin_reorder_preserve_selected_pairs_and_reindex_all_bottoms() {
    let (frame, first) = library();
    let old = expected_bindings(&first);
    let choices = first.alternatives.as_ref().unwrap().rows[0].choices.clone();
    let second = alternatives::add_row(&first, choices).unwrap();
    let new = expected_bindings(&second);
    assert_ne!(old[0].input_id, new[0].input_id);
    assert_ne!(old[2].input_id, new[2].input_id);
    assert_eq!(
        old.iter().map(|pin| &pin.event).collect::<Vec<_>>(),
        new.iter().map(|pin| &pin.event).collect::<Vec<_>>()
    );
    assert!(
        second
            .methods
            .rows
            .iter()
            .all(|raw| raw.anneal_candidate == Some(0))
    );
    let bins = second
        .alternatives
        .as_ref()
        .unwrap()
        .bins
        .iter()
        .rev()
        .map(|bin| bin.id.clone())
        .collect();
    let proposal = topology::propose(&second, topology::Change::Reorder { bin_ids: bins }).unwrap();
    let reordered = proposal.after();
    reordered.validate().unwrap();
    assert_eq!(reordered.methods.rows[0].anneal_bottom_document, Some(2));
    assert_eq!(reordered.methods.rows[1].anneal_bottom_document, Some(3));
    assert!(reordered.slots[0].anneal_bottom == second.slots[1].anneal_bottom);
    assert_eq!(reordered.methods.rows[0].anneal_candidate, Some(0));
    let projected = project(reordered).unwrap();
    let design = projected.request.design.as_ref().unwrap();
    for index in 0..projected.request.rows.len() {
        let selected = row(&projected.request, index);
        let resolved = design::resolve(design, selected, &frame).unwrap();
        assert_eq!(
            resolved.bindings,
            design.selected_input_bindings(selected).unwrap()
        );
        assert_eq!(resolved.inputs.len(), 4);
    }
}

#[test]
fn equal_dna_from_a_distinct_bottom_event_invalidates_choice_freeze_and_requires_capture() {
    let (mut frame, original) = library();
    let same = captured(&frame, 3).molecule.clone();
    DocumentHistory::new("same-dna-other-bottom", same)
        .unwrap()
        .contribute(&mut frame)
        .unwrap();
    let alternate = captured(&frame, 4);
    assert_eq!(
        alternate.source.binding.revision,
        captured(&frame, 3).source.binding.revision
    );
    let before = project(&original).unwrap();
    let mut edited = original.clone();
    edited.active_slot = 1;
    edited.slots[1].anneal_bottom = Some(alternate.clone());
    edited.sync_rows();
    assert_eq!(edited.methods.rows[1].anneal_candidate, None);
    assert!(alternatives::require_active_choice_match(&edited).is_err());
    let changed = project(&edited).unwrap();
    assert_ne!(changed.revision, before.revision);
    assert!(matches!(
        changed.request.rows[0],
        RowRequest::DraftError { .. }
    ));
    let mut captured = alternatives::capture_active_choice(&edited).unwrap();
    assert_eq!(
        captured.alternatives.as_ref().unwrap().bins[1]
            .choices
            .len(),
        2
    );
    assert_eq!(captured.methods.rows[1].anneal_candidate, None);
    captured.methods.rows[1].anneal_candidate = Some(0);
    let projection = project(&captured).unwrap();
    let design = projection.request.design.as_ref().unwrap();
    let selected = row(&projection.request, 0);
    assert_eq!(
        design.selected_input_bindings(selected).unwrap()[3].event,
        alternate.source.binding
    );
    assert_eq!(
        design::resolve(design, selected, &frame)
            .unwrap()
            .inputs
            .len(),
        4
    );
}

#[test]
fn pair_or_options_changes_invalidate_intent_and_a_different_preparation_cannot_drop_bottoms() {
    let (_, first) = library();
    let mut changed = first.clone();
    changed.methods.rows[1].anneal_options.min_paired_bases = 5;
    changed.sync_rows();
    assert_eq!(changed.methods.rows[1].anneal_candidate, None);
    assert!(matches!(
        project(&changed).unwrap().request.rows[0],
        RowRequest::DraftError { .. }
    ));
    let swapped = first.swap_oligo_roles(1).unwrap();
    assert_eq!(swapped.methods.rows[1].anneal_candidate, None);
    assert!(alternatives::require_active_choice_match(&swapped).is_err());
    let mut ignored = first.clone();
    ignored.methods.rows[1].preparation = OverlapPreparation::WholeLinear;
    ignored.sync_rows();
    assert!(matches!(
        project(&ignored).unwrap().request.rows[0],
        RowRequest::DraftError { .. }
    ));
    let mut reference = first;
    reference.methods.mode = Mode::MultiSiteReference;
    assert!(reference.validate().is_err());
}

#[test]
fn saved_design_and_batch_release_keep_every_ordered_event_after_native_reload() {
    let (frame, draft) = library();
    let expected = expected_bindings(&draft);
    let (saved_frame, entry) = saved::save(&frame, "Paired oligo alternatives", &draft).unwrap();
    assert_eq!(entry.source_count, 4);
    let restored = saved::load(&saved_frame, &entry.id).unwrap();
    assert!(restored == draft);
    let projection = project(&restored).unwrap();
    let selected_id = row(&projection.request, 0).id.clone();
    let mut session = batch::Session::begin(saved_frame, projection.request).unwrap();
    let context = session.summary().context_revision;
    session.advance(&context, 0).unwrap();
    let batch = session.into_capture(&context).unwrap();
    let bytes = collagene_core::documents::workspace::save_workspace(&batch.workspace).unwrap();
    let reloaded = collagene_core::documents::workspace::load_workspace(&bytes).unwrap();
    let record = batch::load(&reloaded, &batch.attachment_hash).unwrap();
    assert_eq!(record.rows[0].input_bindings, expected);
    let request = release::Request {
        batch_hash: batch.attachment_hash,
        rows: vec![release::Selection {
            construct_id: selected_id,
            options: outputs::Options {
                name: "Released paired-oligo construct".into(),
                ..Default::default()
            },
        }],
    };
    let mut release = release::ReleaseSession::begin(reloaded, request).unwrap();
    let context = release.summary().context_revision;
    release.advance(&context, 0).unwrap();
    let inspected = release.inspect(0).unwrap();
    assert_eq!(inspected.row.inputs, expected);
    assert_eq!(inspected.product.sequence.as_str(), "ACGTTAAAGGCC");
}
