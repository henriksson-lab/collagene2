//! End-to-end inverse-PCR outputs against a separately declared literal target.
use collagene_core::{
    cloning::pcr::mutagenesis::{self, InversePcrOptions},
    documents::{history::DocumentHistory, workspace::Workspace},
    formats::oligos::{Delimiter, OligoTable},
    material::{
        lineage::{self as run_lineage, EventBinding, InputEventBinding},
        outputs::{self as cloning_outputs, Options},
    },
    planning::{
        oligos as plan_oligos,
        plan::{DesignPlan, Operation},
    },
    sequence::analysis,
    verification::evidence::Status,
};
use collagene_model::{Molecule, Topology};

const ORIGINAL: &str = "AGTCCGATGCTACGTAGGACTTCAGACCTGTA";
// Delete original positions 3..5, then insert GGA before original position 5.
// This expected reference is specified before the preparation simulator runs.
const INTENDED: &str = "GGAGATGCTACGTAGGACTTCAGACCTGTAAGT";
const PRODUCT_ORIGIN: &str = "GAGATGCTACGTAGGACTTCAGACCTGTAAGTG";

fn binding(frame: &Workspace, index: usize) -> EventBinding {
    let document = &frame.documents[index];
    EventBinding {
        document_id: document.id.clone(),
        event_id: document.head_event.clone().unwrap(),
        revision: document.head.clone(),
    }
}

fn fixture() -> (Workspace, DesignPlan, EventBinding, InputEventBinding) {
    let source = Molecule::dna("Original circle", ORIGINAL, Topology::Circular).unwrap();
    let target = Molecule::dna("Declared replacement", INTENDED, Topology::Circular).unwrap();
    let design = mutagenesis::design_inverse_pcr(
        &source,
        &target,
        None,
        &InversePcrOptions {
            edit_start: 3,
            edit_end: 5,
            inserted: "GGA".into(),
            forward_length: 8,
            reverse_length: 8,
            reverse_tail_bases: 1,
        },
    )
    .unwrap();
    assert_eq!(
        design.verification.status,
        Status::ConsistentUnderAssumptions
    );
    let mut plan = design.plan;
    let Operation::Pcr {
        forward, reverse, ..
    } = &mut plan.steps[0]
    else {
        panic!("inverse-PCR begins with the declared oligo pair")
    };
    forward.name = "Replacement forward".into();
    reverse.name = "Replacement reverse".into();
    assert_eq!(forward.tail, "GA");
    assert_eq!(reverse.tail, "C");
    // Same-content copies exist, but their history events remain different.
    let frame =
        Workspace::from_documents("inverse outputs", vec![source.clone(), target, source], 1)
            .unwrap();
    let target = binding(&frame, 1);
    let input = InputEventBinding {
        input_id: plan.inputs[0].id.clone(),
        event: binding(&frame, 0),
    };
    (frame, plan, target, input)
}

fn advance_original(frame: &mut Workspace, original: &EventBinding) -> EventBinding {
    let mut histories = DocumentHistory::restore_all(frame).unwrap();
    let source = histories
        .iter_mut()
        .find(|document| document.id == original.document_id)
        .unwrap();
    source
        .commit(
            Molecule::dna(
                "Later source edit",
                "TGTCCGATGCTACGTAGGACTTCAGACCTGTA",
                Topology::Circular,
            )
            .unwrap(),
            "Edit source after pinning inverse-PCR design",
        )
        .unwrap();
    let advanced = EventBinding {
        document_id: source.id.clone(),
        event_id: source.head_event().into(),
        revision: source.head_revision().into(),
    };
    frame.documents.clear();
    frame.history.clear();
    for history in &histories {
        history.contribute(frame).unwrap();
    }
    frame.validate().unwrap();
    advanced
}

#[test]
fn inverse_output_keeps_original_event_named_orders_and_hidden_preparations_after_source_edit() {
    let (mut frame, plan, target, input) = fixture();
    let advanced = advance_original(&mut frame, &input.event);
    let before = run_lineage::workspace_fingerprint(&frame).unwrap();
    let options = Options {
        name: "Reviewed replacement".into(),
        save_primers: true,
        ..Options::default()
    };
    let commit =
        cloning_outputs::finalize(&frame, &plan, target.clone(), vec![input.clone()], &options)
            .unwrap();
    assert_eq!(run_lineage::workspace_fingerprint(&frame).unwrap(), before);
    let construct = commit
        .workspace
        .head(&commit.receipt.construct.document_id)
        .unwrap();
    assert_eq!(construct.name, "Reviewed replacement");
    assert_eq!(construct.topology, Topology::Circular);
    assert_eq!(construct.sequence.as_str(), PRODUCT_ORIGIN);
    assert!(analysis::rotation_offset(PRODUCT_ORIGIN, INTENDED).is_some());
    assert_eq!(commit.receipt.run.products.len(), 3);
    assert_eq!(commit.receipt.hidden_products.len(), 2);
    assert_eq!(
        cloning_outputs::load_hidden(&commit.workspace)
            .unwrap()
            .len(),
        2
    );
    let saved = run_lineage::load_run(&commit.workspace, &commit.receipt.run.run_id).unwrap();
    assert_eq!(saved.target, target);
    assert_eq!(saved.inputs[0].event, input.event);
    assert_ne!(saved.inputs[0].event, advanced);
    assert_ne!(saved.inputs[0].event, binding(&frame, 2));
    assert_eq!(
        commit
            .workspace
            .head(&advanced.document_id)
            .unwrap()
            .revision(),
        advanced.revision
    );
    let first = &saved.products[0];
    let first_event = commit
        .workspace
        .history
        .iter()
        .find(|event| event.id == first.event_id)
        .unwrap();
    assert!(first_event.parents.contains(&input.event.event_id));
    let replay = run_lineage::replay_run(&commit.workspace, &saved.run_id).unwrap();
    assert!(replay.products_match && replay.status_matches);
    let table: OligoTable = serde_json::from_slice(
        &commit.workspace.attachments[commit.receipt.primer_attachment.as_ref().unwrap()].bytes,
    )
    .unwrap();
    assert_eq!(table.records.len(), 2);
    assert_eq!(table.records[0].name, "Replacement forward");
    assert_eq!(table.records[1].name, "Replacement reverse");
    assert_eq!(table.records[0].raw_tail.as_deref(), Some("GA"));
    assert_eq!(table.records[1].raw_tail.as_deref(), Some("C"));
    let orders = plan_oligos::export(
        &plan,
        frame.head(&target.document_id).unwrap(),
        &Default::default(),
        Delimiter::Comma,
    )
    .unwrap();
    assert!(orders.report.complete);
    assert!(orders.export.text.contains("Replacement forward"));
    assert!(orders.export.text.contains("Replacement reverse"));
    for row in &table.records {
        assert!(
            orders
                .export
                .text
                .contains(row.raw_ordered_sequence.as_ref().unwrap())
        );
    }
    let shown = cloning_outputs::finalize(
        &commit.workspace,
        &plan,
        target,
        vec![input],
        &Options {
            show_prepared: true,
            ..options
        },
    )
    .unwrap();
    assert!(shown.receipt.hidden_products.is_empty());
    assert_eq!(
        shown.workspace.history.len(),
        commit.workspace.history.len()
    );
}

#[test]
fn mismatched_exact_events_and_changed_input_heads_fail_without_output_mutation() {
    let (mut frame, plan, target, input) = fixture();
    let advanced = advance_original(&mut frame, &input.event);
    let before = run_lineage::workspace_fingerprint(&frame).unwrap();
    let mut bad = input.clone();
    // Mixing a same-content copy's identity with the original event is invalid.
    bad.event.document_id = frame.documents[2].id.clone();
    assert!(
        cloning_outputs::finalize(
            &frame,
            &plan,
            target.clone(),
            vec![bad],
            &Options::default()
        )
        .is_err()
    );
    let mut bad = input.clone();
    bad.event = advanced;
    assert!(
        cloning_outputs::finalize(
            &frame,
            &plan,
            target.clone(),
            vec![bad],
            &Options::default()
        )
        .is_err()
    );
    let mut bad_target = target;
    bad_target.event_id = input.event.event_id.clone();
    assert!(
        cloning_outputs::finalize(&frame, &plan, bad_target, vec![input], &Options::default())
            .is_err()
    );
    assert_eq!(run_lineage::workspace_fingerprint(&frame).unwrap(), before);
}
