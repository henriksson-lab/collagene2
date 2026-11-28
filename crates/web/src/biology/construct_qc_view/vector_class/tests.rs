//! Synthetic UI declarations preserve exact parent premises and event identity.
use super::model::{Capture, Form};
use crate::{biology::construct_qc_model::Draft, sources::identity::SourceInfo};
use collagene_core::{
    construct_qc::{self as qc, vector_class as vector},
    material::lineage::EventBinding,
};
use collagene_model::{Feature, Location, Molecule, Sequence, Topology};

fn fixture() -> (Molecule, SourceInfo, Draft, Form) {
    let mut source = Molecule::dna(
        "Synthetic graph",
        "GGAACCGGCGTATTATGCCCGGTTAAAA",
        Topology::Linear,
    )
    .unwrap();
    source.features = [
        (2, 6, false),
        (8, 12, false),
        (14, 18, false),
        (20, 24, true),
    ]
    .into_iter()
    .map(|(start, end, reverse)| {
        let location = Location::range(start, end);
        Feature {
            kind: "misc_feature".into(),
            location: if reverse {
                Location::Complement(Box::new(location))
            } else {
                location
            },
            qualifiers: vec![],
            needs_review: false,
        }
    })
    .collect();
    let mut parent = Draft {
        assumed: vec![0, 1, 2, 3],
        ..Default::default()
    };
    parent.options.references = [
        ("repeat", "AACC"),
        ("regulatory", "CGTA"),
        ("payload", "ATGC"),
    ]
    .into_iter()
    .map(|(id, dna)| {
        let sequence = Sequence::dna(dna).unwrap();
        qc::ReferenceElement {
            id: id.into(),
            sequence_revision: qc::reference_revision(&sequence).unwrap(),
            sequence,
            provenance: qc::Provenance {
                source_url: "https://example.org/authored-fixture".into(),
                accession: format!("synthetic-{id}"),
                retrieved_on: "2026-11-28".into(),
                license: "CC0 authored synthetic fixture".into(),
            },
        }
    })
    .collect();
    let profile = vector::catalogue()
        .unwrap()
        .into_iter()
        .find(|profile| profile.reference.id == "raav-transfer-declared-graph-2026")
        .unwrap();
    let mut form = Form::for_profile(&profile).unwrap();
    for row in &mut form.bindings {
        let (index, reference) = match row.role.as_str() {
            "five-prime-itr" => (0, "repeat"),
            "promoter" => (1, "regulatory"),
            "payload" => (2, "payload"),
            "three-prime-itr" => (3, "repeat"),
            _ => unreachable!(),
        };
        row.feature = index.to_string();
        row.reference = reference.into();
    }
    form.first = form.bindings[0].id.clone();
    form.last = form
        .bindings
        .iter()
        .find(|row| row.role == "three-prime-itr")
        .unwrap()
        .id
        .clone();
    let info = SourceInfo {
        binding: EventBinding {
            document_id: "document-a".into(),
            event_id: format!("event-{}", "a".repeat(64)),
            revision: source.revision(),
        },
        name: source.name.clone(),
        length: source.sequence.len(),
        alphabet: source.sequence.alphabet(),
        topology: source.topology,
    };
    (source, info, parent, form)
}
#[test]
fn shared_premises_and_boundary_coordinates_survive_conversion() {
    let (source, info, parent, form) = fixture();
    let capture = Capture::prepare(info, 7, parent.clone(), form, &source).unwrap();
    let options = parent.options_for(&source).unwrap();
    assert_eq!(capture.request.annotations, options.annotations);
    assert_eq!(capture.request.references, options.references);
    assert_eq!(capture.request.literal_parts, options.required_parts);
    assert_eq!(
        capture.request.literal_part_limits,
        options.required_part_limits
    );
    let report = vector::assess(&source, &capture.request).unwrap();
    capture.validate_report(&report).unwrap();
    let boundary = report.boundaries.unwrap();
    assert_eq!((boundary.outer_bases, boundary.interior_bases), (22, 14));
}
#[test]
fn unfinished_rows_and_boundary_pairs_are_explicit_not_inferred() {
    let (source, _, parent, mut form) = fixture();
    form.bindings[0].reference.clear();
    assert!(form.request(&source, &parent).is_err());
    form.bindings[0].feature.clear();
    form.first.clear();
    form.last.clear();
    let request = form.request(&source, &parent).unwrap();
    assert_eq!(request.bindings.len(), 3);
    assert!(
        vector::assess(&source, &request)
            .unwrap()
            .boundaries
            .is_none()
    );
    form.first = "only-one".into();
    assert!(form.request(&source, &parent).is_err());
    assert!(Form::default().request(&source, &parent).is_err());
}
#[test]
fn unchanged_sequence_does_not_hide_event_epoch_or_settings_changes() {
    let (source, info, parent, form) = fixture();
    let capture = Capture::prepare(info.clone(), 7, parent.clone(), form.clone(), &source).unwrap();
    capture.guard(&info, 7, &parent, &form).unwrap();
    let mut changed = info.clone();
    changed.binding.event_id = format!("event-{}", "b".repeat(64));
    assert!(capture.guard(&changed, 7, &parent, &form).is_err());
    assert!(capture.guard(&info, 8, &parent, &form).is_err());
    let mut settings = parent.clone();
    settings.assumed.pop();
    assert!(capture.guard(&info, 7, &settings, &form).is_err());
    let mut edited = form.clone();
    edited.reverse = true;
    assert!(capture.guard(&info, 7, &parent, &edited).is_err());
}
#[test]
fn resource_and_spacing_constraints_refuse_before_dispatch() {
    let (source, _, parent, mut form) = fixture();
    form.add_spacing().unwrap();
    let row = &mut form.spacing[0];
    row.upstream = form.bindings[0].id.clone();
    row.downstream = form.bindings[1].id.clone();
    row.rationale = "Declared synthetic spacing".into();
    row.minimum = "-1".into();
    assert!(form.request(&source, &parent).is_err());
    form.spacing[0].minimum = "0".into();
    form.spacing[0].maximum = "1.5".into();
    assert!(form.request(&source, &parent).is_err());
    form.spacing.clear();
    while form.bindings.len() < 64 {
        form.add_binding("payload").unwrap();
    }
    assert!(form.add_binding("payload").is_err());
}
