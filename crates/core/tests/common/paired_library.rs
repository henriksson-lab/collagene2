use collagene_core::{
    chemistry::{AnnealOptions, Terminus},
    cloning::assembly::{intent::Selector, recipe::Preparation},
    documents::workspace::Workspace,
    material::lineage::EventBinding,
    output::layout,
    planning::combinatorial::{self as library, batch},
};
use collagene_model::{Molecule, Topology};

pub fn event(frame: &Workspace, index: usize) -> EventBinding {
    let document = &frame.documents[index];
    EventBinding {
        document_id: document.id.clone(),
        event_id: document.head_event.clone().unwrap(),
        revision: document.head.clone(),
    }
}

/// Two bottom declarations have identical DNA/revision but distinct retained
/// events. The unselected declaration is still part of the saved design closure.
pub fn fixture() -> (Workspace, layout::Request, String) {
    let bottom = Molecule::dna("Independent bottom", "GGTT", Topology::Linear).unwrap();
    let frame = Workspace::from_documents(
        "Paired library sources",
        vec![
            Molecule::dna("Original top", "AACC", Topology::Linear).unwrap(),
            bottom.clone(),
            bottom,
            Molecule::dna("Private unrelated", "ACGTCGTA", Topology::Linear).unwrap(),
        ],
        0,
    )
    .unwrap();
    let options = AnnealOptions {
        min_paired_bases: 4,
        top_chemistry: Terminus::CLEAVED,
        bottom_chemistry: Terminus::CLEAVED,
        ..Default::default()
    };
    let choices = [1, 2]
        .into_iter()
        .map(|index| library::Choice {
            id: format!("choice-{index}"),
            label: format!("Pair {index}"),
            source: event(&frame, 0),
            anneal_bottom: Some(event(&frame, index)),
            selector: Selector::ReferenceRange {
                start: 0,
                end: 4,
                reverse: false,
            },
            preparation_reverse: false,
            preparation: Some(Preparation::AnnealedOligos {
                bottom_input_id: "explicit-template-bottom".into(),
                options: options.clone(),
                candidate_index: 0,
            }),
            treatment: None,
            generation: None,
        })
        .collect();
    let request = batch::Request {
        design: Some(library::Design {
            id: "paired-library".into(),
            name: "Explicit paired library".into(),
            circular: false,
            bins: vec![library::Bin {
                id: "insert".into(),
                label: "Oligo insert".into(),
                choices,
            }],
            method: library::Method::Ligation {
                allow_nicked: false,
            },
            material_bindings: vec![],
        }),
        rows: vec![batch::RowRequest::Ready(library::SelectedConstruct {
            standard: None,
            id: "pair-A".into(),
            name: "Independent AACC intent".into(),
            choices: vec![library::ChoiceSelection {
                bin_id: "insert".into(),
                choice_id: "choice-1".into(),
            }],
            target: library::TargetDeclaration::OriginalRegions {
                junction_edits: vec![],
                substitutions: vec![],
                edits: vec![],
            },
            overrides: None,
        })],
    };
    let unselected_bottom = event(&frame, 2).event_id;
    let mut session = batch::Session::begin(frame, request).unwrap();
    let context = session.summary().context_revision;
    let row = session.advance(&context, 0).unwrap().row.unwrap();
    assert_eq!(row.input_bindings.len(), 2);
    let captured = session.into_capture(&context).unwrap();
    let layout = layout::Request {
        name: "Paired oligo plate".into(),
        plate: layout::Plate::Wells96,
        assignments: vec![layout::Assignment {
            id: "pair-A1".into(),
            well: layout::Well { row: 0, column: 0 },
            label: "Paired AACC".into(),
            source: layout::Source::Planned {
                batch_hash: captured.attachment_hash,
                construct_id: "pair-A".into(),
                declaration_revision: row.declaration_revision.unwrap(),
            },
        }],
    };
    (captured.workspace, layout, unselected_bottom)
}
