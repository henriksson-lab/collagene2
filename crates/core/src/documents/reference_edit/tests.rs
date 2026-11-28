//! Literal native edit/event fixtures, including equal-DNA distinct sources.
use super::*;
use crate::{
    documents::history::DocumentHistory,
    material::lineage::EventBinding,
    verification::ancestry::{IntroductionKind, Limits, SourceHandle, TraceBudget},
};
use collagene_model::{Location, Sequence};
fn source() -> Molecule {
    let mut source =
        Molecule::dna("Literal native reference", "AACCGGTT", Topology::Linear).unwrap();
    source.features = vec![feature(0, 2), feature(4, 8)];
    source
}
fn feature(start: usize, end: usize) -> Feature {
    Feature {
        kind: "misc_feature".into(),
        location: Location::range(start, end),
        qualifiers: vec![("label".into(), Some("Repeated label".into()))],
        needs_review: false,
    }
}
fn binding(id: &str, source: &Molecule) -> EventBinding {
    let document = DocumentHistory::new(id, source.clone()).unwrap();
    EventBinding {
        document_id: document.id.clone(),
        event_id: document.head_event().into(),
        revision: source.revision(),
    }
}
fn captured(
    source: &Molecule,
    event: &EventBinding,
    actions: &[Operation],
) -> (Molecule, crate::verification::ancestry::ReconstructionTrace) {
    let (product, trace) = Operation::apply_all_with_ancestry(
        source,
        event,
        actions,
        TraceBudget::new(Limits::default()).unwrap(),
    )
    .unwrap();
    assert_eq!(product, Operation::apply_all(source, actions).unwrap());
    (product, trace.unwrap())
}
#[test]
fn normalized_literal_edit_and_reverse_follow_actual_intermediates_and_introductions() {
    let source = source();
    let event = binding("literal-source", &source);
    let actions = vec![
        Operation::ReplaceSequence {
            start: 2,
            end: 4,
            replacement: " t a\n".into(),
        },
        Operation::ReverseComplement,
    ];
    let (product, trace) = captured(&source, &event, &actions);
    assert_eq!(product.sequence.as_str(), "AACCTATT");
    let mut lookup = trace.lookup_session();
    let literal = lookup.at(0, 2, None).unwrap();
    assert!(literal.complete && literal.origins.is_empty());
    assert_eq!(literal.introductions.len(), 1);
    assert_eq!(
        *literal.introductions[0].kind,
        IntroductionKind::DeclaredLiteral
    );
    let reversed = lookup.at(1, 0, None).unwrap();
    assert!(reversed.complete);
    assert_eq!(reversed.origins.len(), 1);
    assert_eq!(reversed.origins[0].position, 7);
    assert!(reversed.origins[0].reverse);
    assert_eq!(
        reversed.origins[0].source.source,
        SourceHandle::ReferenceEvent(event)
    );
    let reversed_literal = lookup.at(1, 4, None).unwrap();
    assert!(reversed_literal.origins.is_empty());
    assert_eq!(reversed_literal.introductions.len(), 1);
    assert!(reversed_literal.introductions[0].reverse);
}
#[test]
fn native_deletion_retains_source_feature_ordinals_and_reassignment_severs_old_identity() {
    let source = source();
    let event = binding("annotated-source", &source);
    let (deleted, projection) = source.replace_with_projection(0, 4, "").unwrap();
    assert_eq!(deleted.sequence.as_str(), "GGTT");
    assert_eq!(projection.source_feature_indices, vec![1]);
    assert_eq!(projection.dropped_source_feature_indices, vec![0]);
    assert_eq!(
        (
            projection.map[0].source_start,
            projection.map[0].destination_start
        ),
        (4, 0)
    );
    let actions = vec![
        Operation::DeleteFeature { index: 0 },
        Operation::SetFeature {
            index: Some(0),
            feature: feature(0, 2),
        },
    ];
    let (_, trace) = captured(&source, &event, &actions);
    let mut lookup = trace.lookup_session();
    let inherited = lookup.feature_at(0, 0, None).unwrap();
    assert!(inherited.complete);
    assert_eq!(inherited.origins.len(), 1);
    assert_eq!(inherited.origins[0].feature_index, 1);
    let reassigned = lookup.feature_at(1, 0, None).unwrap();
    assert!(!reassigned.complete && reassigned.origins.is_empty());
    assert!(
        lookup.at(1, 0, None).unwrap().complete,
        "DNA remains copied even when annotation identity is reassigned"
    );
}
#[test]
fn palindrome_reverse_rotation_and_extraction_do_not_choose_identity_by_sequence() {
    let mut source = Molecule::dna("Palindrome", "ACGTACGT", Topology::Circular).unwrap();
    source.features.push(feature(0, 4));
    let event = binding("palindrome-source", &source);
    let actions = vec![
        Operation::ReverseComplement,
        Operation::RotateOrigin { new_origin: 2 },
        Operation::ExtractRegion { start: 6, end: 10 },
    ];
    let (product, trace) = captured(&source, &event, &actions);
    assert_eq!(product.sequence.as_str(), "ACGT");
    let mut lookup = trace.lookup_session();
    let first = lookup.at(0, 0, None).unwrap();
    assert_eq!(first.origins[0].position, 7);
    assert!(first.origins[0].reverse);
    let last = lookup.at(2, 0, None).unwrap();
    assert!(last.complete);
    assert_eq!(last.origins[0].position, 7);
    assert!(last.origins[0].reverse);
    assert_eq!(
        last.origins[0].source.source,
        SourceHandle::ReferenceEvent(event)
    );
}
#[test]
fn paste_binds_exact_donor_event_even_when_original_molecule_bytes_are_equal() {
    use crate::sequence::reference_transfer::{self as transfer, Part, RequestCopy, Selection};
    let mut source = Molecule::dna("Equal immutable reference", "ACGT", Topology::Linear).unwrap();
    source.features.push(feature(0, 2));
    let destination = binding("destination", &source);
    let donor = binding("donor", &source);
    assert_ne!(destination.event_id, donor.event_id);
    assert_eq!(destination.revision, donor.revision);
    let selection = Selection {
        parts: vec![Part {
            start: 0,
            end: 4,
            reverse: false,
        }],
        reverse: false,
    };
    let revision = transfer::copy_revision(&RequestCopy {
        source: source.clone(),
        selection: selection.clone(),
    })
    .unwrap();
    let actions = vec![Operation::PasteReference {
        source: source.clone(),
        source_event: donor.clone(),
        selection,
        fragment_revision: revision,
        start: 2,
        end: 2,
    }];
    let (product, trace) = captured(&source, &destination, &actions);
    assert_eq!(product.sequence.as_str(), "ACACGTGT");
    let mut lookup = trace.lookup_session();
    for (position, root, source_position) in
        [(0, &destination, 0), (2, &donor, 0), (6, &destination, 2)]
    {
        let call = lookup.at(0, position, None).unwrap();
        assert!(call.complete && call.introductions.is_empty());
        assert_eq!(call.origins.len(), 1);
        assert_eq!(call.origins[0].position, source_position);
        assert_eq!(
            call.origins[0].source.source,
            SourceHandle::ReferenceEvent(root.clone())
        );
    }
    let pasted_feature = lookup.feature_at(0, 1, None).unwrap();
    assert!(pasted_feature.complete);
    assert_eq!(
        pasted_feature.origins[0].source.source,
        SourceHandle::ReferenceEvent(donor)
    );
    assert_eq!(pasted_feature.origins[0].feature_index, 0);
}
#[test]
fn annotation_batch_and_import_keep_calls_but_never_recover_reassigned_identity() {
    let source = source();
    let event = binding("annotation-curation", &source);
    let request = annotations::batch::Request {
        source_revision: source.revision(),
        indices: vec![0],
        kind: None,
        flip_strand: true,
        color: annotations::batch::ColorEdit::Keep,
    };
    let (_, trace) = captured(&source, &event, &[Operation::EditAnnotations { request }]);
    let mut lookup = trace.lookup_session();
    assert!(!lookup.feature_at(0, 0, None).unwrap().complete);
    assert_eq!(
        lookup.feature_at(0, 1, None).unwrap().origins[0].feature_index,
        1
    );
    let (_, trace) = captured(
        &source,
        &event,
        &[Operation::FlipAnnotationStrands { indices: vec![1] }],
    );
    assert!(
        !trace
            .lookup_session()
            .feature_at(0, 1, None)
            .unwrap()
            .complete
    );
    let report = crate::formats::annotations::import_annotations(
        &source,
        "##gff-version 3\nref\tfixture\tmisc_feature\t3\t4\t.\t+\t.\tID=new\n",
        crate::formats::annotations::AnnotationFormat::Gff3,
        "ref",
    )
    .unwrap();
    let (_, trace) = captured(
        &source,
        &event,
        &[Operation::ImportAnnotations {
            report,
            selected: vec![0],
        }],
    );
    let mut lookup = trace.lookup_session();
    assert!(lookup.feature_at(0, 0, None).unwrap().complete);
    assert!(
        !lookup.feature_at(0, 2, None).unwrap().complete,
        "new annotation has no invented earlier feature identity"
    );
}
#[test]
fn exhausted_or_foreign_trace_context_preserves_native_reference_results() {
    let source = source();
    let event = binding("bounded-native", &source);
    let actions = vec![
        Operation::Rename {
            name: "Actual renamed reference".into(),
        },
        Operation::ReplaceSequence {
            start: 0,
            end: 1,
            replacement: "T".into(),
        },
    ];
    let expected = Operation::apply_all(&source, &actions).unwrap();
    let mut limits = Limits::default();
    limits.max_scanned_bytes = 1;
    let (actual, trace) = Operation::apply_all_with_ancestry(
        &source,
        &event,
        &actions,
        TraceBudget::new(limits).unwrap(),
    )
    .unwrap();
    assert_eq!(actual, expected);
    assert!(trace.is_err());
    let mut foreign = event;
    foreign.revision = "0".repeat(64);
    let (actual, trace) = Operation::apply_all_with_ancestry(
        &source,
        &foreign,
        &actions,
        TraceBudget::new(Limits::default()).unwrap(),
    )
    .unwrap();
    assert_eq!(actual, expected);
    assert!(trace.is_err());
    let (product, receipt) = source.replace_with_projection(0, 4, " a a c c\n").unwrap();
    assert_eq!(product, source);
    assert!(receipt.unchanged);
    assert_eq!(receipt.source_feature_indices, vec![0, 1]);
    assert_eq!(receipt.map.len(), 1);
    assert_eq!(receipt.map[0].source_end, 8);
    assert_eq!(Sequence::dna("AACCGGTT").unwrap(), product.sequence);
}

#[test]
fn library_proposals_manual_additions_and_reference_metadata_do_not_invent_prior_annotations() {
    let source = source();
    let event = binding("curated-reference", &source);
    let mut template = Molecule::dna("Library template", "CCGG", Topology::Linear).unwrap();
    template.features.push(feature(0, 4));
    let mut library = annotations::library::FeatureLibrary::new("Literal library", "1").unwrap();
    library.add_feature(&template, 0).unwrap();
    let report = annotations::library::scan(&source, &library, &Default::default()).unwrap();
    let selected = vec![
        report
            .hits
            .iter()
            .find(|hit| !hit.already_annotated)
            .unwrap()
            .id
            .clone(),
    ];
    let (_, trace) = captured(
        &source,
        &event,
        &[Operation::ApplyFeatureLibrary {
            library,
            report,
            selected,
        }],
    );
    let mut lookup = trace.lookup_session();
    assert!(lookup.feature_at(0, 0, None).unwrap().complete);
    assert!(!lookup.feature_at(0, 2, None).unwrap().complete);
    let actions = vec![
        Operation::Rename {
            name: "Explicit user name".into(),
        },
        Operation::SetTopology {
            topology: Topology::Circular,
        },
        Operation::SetFeature {
            index: None,
            feature: feature(2, 3),
        },
    ];
    let (product, trace) = captured(&source, &event, &actions);
    assert_eq!(product.name, "Explicit user name");
    assert_eq!(product.topology, Topology::Circular);
    assert_eq!(product.features.len(), 3);
    let mut lookup = trace.lookup_session();
    assert!(lookup.feature_at(1, 0, None).unwrap().complete);
    assert!(!lookup.feature_at(2, 2, None).unwrap().complete);
    let inherited = lookup.at(2, 7, None).unwrap();
    assert!(inherited.complete);
    assert_eq!(inherited.origins[0].position, 7);
    assert_eq!(
        inherited.origins[0].source.source,
        SourceHandle::ReferenceEvent(event)
    );
}

#[test]
fn compact_transform_receipt_retains_native_unmapped_review_disposition() {
    let mut source = source();
    source.features.push(Feature {
        kind: "misc_feature".into(),
        location: Location::Unresolved("remote:1..4".into()),
        qualifiers: vec![("label".into(), Some("Unresolved source evidence".into()))],
        needs_review: false,
    });
    let (product, receipt) = annotations::reverse_complement_with_projection(&source).unwrap();
    assert_eq!(receipt.unmapped_source_feature_indices, vec![2]);
    assert_eq!(receipt.unmapped_review, vec![true]);
    assert!(!source.features[2].needs_review);
    let evidence: Vec<Feature> =
        serde_json::from_str(&product.metadata["annotation_mapping:unmapped_features"]).unwrap();
    assert!(
        evidence
            .iter()
            .any(|feature| feature.label() == "Unresolved source evidence" && feature.needs_review)
    );
    let event = binding("unresolved-source", &source);
    let (captured, _) = captured(&source, &event, &[Operation::ReverseComplement]);
    assert_eq!(captured, product);
}
