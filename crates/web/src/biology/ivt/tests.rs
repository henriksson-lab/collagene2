use super::model::{Capture, Draft};
use collagene_core::{material::lineage::EventBinding, rna_chemistry::ivt};
use collagene_model::{Molecule, Topology};

fn fixture() -> (Molecule, EventBinding, Draft) {
    let m = Molecule::dna("IVT template", "TAATACGACTCACTATAAGAAAA", Topology::Linear).unwrap();
    let binding = EventBinding {
        document_id: "ivt-document".into(),
        event_id: format!("event-{}", "a".repeat(64)),
        revision: m.revision(),
    };
    let draft = Draft {
        first_base: "18".into(),
        length: "6".into(),
        cap: ivt::CapChoice::CleanCapAgN7113,
        ..Draft::default()
    };
    (m, binding, draft)
}

#[test]
fn identical_sequence_does_not_authorize_another_event_epoch_or_form() {
    let (m, binding, draft) = fixture();
    let capture = Capture::prepare(binding.clone(), 7, draft.clone(), &m).unwrap();
    capture.guard(&binding, 7, &draft).unwrap();
    let mut changed = binding.clone();
    changed.event_id = format!("event-{}", "b".repeat(64));
    assert!(capture.guard(&changed, 7, &draft).is_err());
    changed = binding.clone();
    changed.document_id = "another-document".into();
    assert!(capture.guard(&changed, 7, &draft).is_err());
    assert!(capture.guard(&binding, 8, &draft).is_err());
    let mut changed_draft = draft;
    changed_draft.cap = ivt::CapChoice::Uncapped;
    assert!(capture.guard(&binding, 7, &changed_draft).is_err());
}

#[test]
fn one_based_poly_a_entry_preserves_inclusive_end_and_optional_cut_pair() {
    let (m, _, mut draft) = fixture();
    draft.poly_a_first = "3".into();
    draft.poly_a_last = "6".into();
    draft.coding_cut = "6".into();
    assert!(
        draft.request(&m).is_err(),
        "one supplied strand cut cannot define an end"
    );
    draft.template_cut = "6".into();
    let request = draft.request(&m).unwrap();
    assert_eq!(request.transcript.start, 17);
    assert_eq!(request.poly_a.as_ref().unwrap().start_offset, 2);
    assert_eq!(request.poly_a.as_ref().unwrap().end_offset, 6);
    assert_eq!(request.linearization.as_ref().unwrap().coding_cut_offset, 6);
    draft.poly_a_first = "0".into();
    assert!(draft.request(&m).is_err());
    draft.poly_a_first = "7".into();
    assert!(draft.request(&m).is_err());
}

#[test]
fn captured_response_rejects_different_settings_and_mutated_findings() {
    let (m, binding, draft) = fixture();
    let capture = Capture::prepare(binding, 1, draft, &m).unwrap();
    let report = ivt::assess(&m, &capture.request).unwrap();
    capture.validate_report(&report).unwrap();
    let mut other_request = capture.request.clone();
    other_request.cap = ivt::CapChoice::Uncapped;
    let other_report = ivt::assess(&m, &other_request).unwrap();
    assert!(capture.validate_report(&other_report).is_err());
    let mut tampered = report;
    tampered.findings[0].message.push_str(" altered");
    assert!(capture.validate_report(&tampered).is_err());
}
