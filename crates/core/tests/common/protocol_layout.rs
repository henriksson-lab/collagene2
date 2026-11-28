use collagene_core::{
    chemistry::Terminus,
    cloning::assembly::recipe::Preparation,
    documents::workspace::Workspace,
    material::lineage::EventBinding,
    output::layout,
    planning::combinatorial::{
        Bin, Choice, ChoiceSelection, Design, Method, SelectedConstruct, TargetDeclaration,
        batch::{self, RowRequest},
    },
};
use collagene_model::{Molecule, Topology};
pub fn fixture() -> (Workspace, layout::Request) {
    let w = Workspace::from_documents(
        "Explicit plate sources",
        vec![
            Molecule::dna("Stock", "AACCTTAA", Topology::Linear).unwrap(),
            Molecule::dna("Private unrelated", "GGGGCCCC", Topology::Linear).unwrap(),
        ],
        0,
    )
    .unwrap();
    let d = &w.documents[0];
    let source = EventBinding {
        document_id: d.id.clone(),
        event_id: d.head_event.clone().unwrap(),
        revision: d.head.clone(),
    };
    let request = batch::Request {
        design: Some(Design {
            material_bindings: vec![],
            id: "library".into(),
            name: "One literal clone".into(),
            circular: false,
            method: Method::Ligation {
                allow_nicked: false,
            },
            bins: vec![Bin {
                id: "insert".into(),
                label: "Insert".into(),
                choices: vec![Choice {
                    anneal_bottom: None,
                    id: "stock".into(),
                    label: "Stock".into(),
                    source,
                    selector: collagene_core::cloning::assembly::intent::Selector::ReferenceRange {
                        start: 0,
                        end: 8,
                        reverse: false,
                    },
                    preparation_reverse: false,
                    preparation: Some(Preparation::DeclaredDuplex {
                        left: Terminus::CLEAVED,
                        right: Terminus::CLEAVED,
                    }),
                    treatment: None,
                    generation: None,
                }],
            }],
        }),
        rows: vec![RowRequest::Ready(SelectedConstruct {
            id: "clone-A".into(),
            name: "Literal AACCTTAA".into(),
            choices: vec![ChoiceSelection {
                bin_id: "insert".into(),
                choice_id: "stock".into(),
            }],
            target: TargetDeclaration::OriginalRegions {
                substitutions: vec![],
                edits: vec![],
                junction_edits: vec![],
            },
            overrides: None,
        })],
    };
    let mut session = batch::Session::begin(w, request).unwrap();
    let context = session.summary().context_revision;
    let row = session.advance(&context, 0).unwrap().row.unwrap();
    let captured = session.into_capture(&context).unwrap();
    let request = layout::Request {
        name: "Explicit cloning plate".into(),
        plate: layout::Plate::Wells96,
        assignments: vec![layout::Assignment {
            id: "assignment-A1".into(),
            well: layout::Well { row: 0, column: 0 },
            label: "clone-A".into(),
            source: layout::Source::Planned {
                batch_hash: captured.attachment_hash,
                construct_id: "clone-A".into(),
                declaration_revision: row.declaration_revision.unwrap(),
            },
        }],
    };
    (captured.workspace, request)
}

pub fn released_fixture() -> (Workspace, layout::Request) {
    use collagene_core::{material::outputs as cloning_outputs, planning::combinatorial::release};
    let (w, planned) = fixture();
    let batch_hash = match &planned.assignments[0].source {
        layout::Source::Planned { batch_hash, .. } => batch_hash.clone(),
        _ => unreachable!(),
    };
    let request = release::Request {
        batch_hash,
        rows: vec![release::Selection {
            construct_id: "clone-A".into(),
            options: cloning_outputs::Options {
                name: "Released clone-A".into(),
                show_prepared: false,
                save_primers: false,
                allow_incomplete: false,
            },
        }],
    };
    let mut session = release::ReleaseSession::begin(w, request).unwrap();
    let context = session.summary().context_revision;
    let row = session.advance(&context, 0).unwrap().row.unwrap();
    let captured = session.into_commit(&context).unwrap();
    let mut request = planned;
    request.assignments[0].source = layout::Source::Released {
        manifest_hash: captured.manifest_hash,
        construct_id: row.construct_id,
        run_id: row.output.run.run_id,
        product: row.output.construct,
    };
    (captured.workspace, request)
}
