//! Literal merge decisions, source identities and geometry; no function claims.
use collagene_core::{
    construct_qc::annotation_revision,
    documents::{history::DocumentHistory, reference_edit::Operation},
    material::lineage::EventBinding,
    sequence::annotations::merge::{
        self, Decision, Limits, RelationKind, Request, Selection, Source, Transform,
    },
};
use collagene_model::{Feature, Location, Molecule, Topology};

fn feature(location: Location, label: &str) -> Feature {
    Feature {
        kind: "misc_feature".into(),
        location,
        qualifiers: vec![("label".into(), Some(label.into()))],
        needs_review: false,
    }
}
fn event(id: &str, molecule: &Molecule) -> EventBinding {
    let history = DocumentHistory::new(id, molecule.clone()).unwrap();
    EventBinding {
        document_id: history.id.clone(),
        event_id: history.head_event().into(),
        revision: molecule.revision(),
    }
}
fn fixture() -> (Molecule, Request) {
    let mut target =
        Molecule::dna("Current annotation target", "AACCGGTTAA", Topology::Linear).unwrap();
    target
        .features
        .push(feature(Location::range(0, 2), "Repeated label"));
    let mut donor = target.clone();
    donor.name = "Independent annotation donor".into();
    donor
        .features
        .push(feature(Location::range(4, 8), "Repeated label"));
    donor
        .features
        .push(feature(Location::range(1, 5), "Overlapping annotation"));
    let request = Request {
        schema_version: 1,
        target_revision: target.revision(),
        source: Source::Donor {
            event: event("native-donor-A", &donor),
            molecule: donor,
            transform: Transform {
                origin: 0,
                reverse: false,
            },
            indices: vec![0, 1, 2],
        },
        selections: vec![
            Selection {
                proposal_id: "donor-0".into(),
                decision: Decision::Skip,
            },
            Selection {
                proposal_id: "donor-1".into(),
                decision: Decision::Append,
            },
            Selection {
                proposal_id: "donor-2".into(),
                decision: Decision::Append,
            },
        ],
        limits: Limits::default(),
    };
    (target, request)
}

#[test]
fn mixed_skip_append_preserves_existing_annotation_and_distinguishes_overlap_from_identity() {
    let (target, request) = fixture();
    let before = target.clone();
    let preview = merge::preview(&target, &request).unwrap();
    assert!(preview.complete);
    assert!(preview.source_complete);
    assert_eq!(preview.candidates.len(), 3);
    assert!(
        preview.candidates[0]
            .relations
            .iter()
            .any(|relation| matches!(
                relation.kind,
                RelationKind::ExactDuplicate | RelationKind::EquivalentAnnotation
            ))
    );
    assert!(preview.candidates[1].relations.is_empty());
    assert!(
        preview.candidates[2]
            .relations
            .iter()
            .any(|relation| relation.kind == RelationKind::Overlap)
    );
    assert!(
        preview
            .candidates
            .iter()
            .all(|candidate| candidate.relations_complete)
    );
    let applied = merge::apply(&target, &request).unwrap();
    assert_eq!(target, before);
    assert_eq!(applied.molecule.sequence, target.sequence);
    assert_eq!(applied.molecule.topology, target.topology);
    assert_eq!(applied.molecule.features[0], target.features[0]);
    assert_eq!(applied.molecule.features.len(), 3);
    assert_eq!(applied.molecule.features[1].location, Location::range(4, 8));
    assert_eq!(applied.molecule.features[2].location, Location::range(1, 5));
    assert!(
        applied.molecule.features[1..]
            .iter()
            .all(|feature| feature.needs_review)
    );
    assert_eq!(applied.review, preview);
    let native = Operation::MergeAnnotations {
        request: request.clone(),
    }
    .apply(&target)
    .unwrap();
    assert_eq!(native, applied.molecule);
    let reopened = collagene_core::protocol::load_document(
        &collagene_core::protocol::save_document(&native).unwrap(),
    )
    .unwrap();
    assert_eq!(reopened, native);
    // A late-invalid choice returns no partially merged result.
    let mut invalid = request.clone();
    invalid.selections.last_mut().unwrap().proposal_id = "donor-999".into();
    assert!(merge::apply(&target, &invalid).is_err());
    assert_eq!(target, before);
    let mut invalid = request;
    invalid.selections.last_mut().unwrap().decision = Decision::AddProvenance {
        target_index: 999,
        annotation_revision: "0".repeat(64),
    };
    assert!(merge::apply(&target, &invalid).is_err());
    assert_eq!(target, before);
}

#[test]
fn stale_target_donor_and_feature_pins_cannot_authorize_captured_decisions() {
    let (target, request) = fixture();
    let preview = merge::preview(&target, &request).unwrap();
    preview.validate_request(&target, &request).unwrap();
    let mut stale = request.clone();
    stale.target_revision = "0".repeat(64);
    assert!(merge::apply(&target, &stale).is_err());
    let mut stale = request.clone();
    if let Source::Donor { event, .. } = &mut stale.source {
        event.revision = "0".repeat(64);
    }
    assert!(merge::apply(&target, &stale).is_err());
    let mut stale = request.clone();
    if let Source::Donor { molecule, .. } = &mut stale.source {
        molecule.features[1]
            .qualifiers
            .push(("note".into(), Some("Changed independently".into())));
    }
    assert!(merge::apply(&target, &stale).is_err());
    let mut stale = request.clone();
    stale.selections[0].decision = Decision::AddProvenance {
        target_index: 0,
        annotation_revision: "0".repeat(64),
    };
    assert!(merge::apply(&target, &stale).is_err());
    let mut tampered = preview;
    tampered.candidates[1].feature.location = Location::range(8, 10);
    assert!(tampered.validate_request(&target, &request).is_err());
}

#[test]
fn explicit_reverse_then_rotation_maps_joined_origin_annotation_in_biological_order() {
    let mut donor = Molecule::dna("Circular source", "AACCGGTA", Topology::Circular).unwrap();
    donor.features.push(feature(
        Location::Complement(Box::new(Location::Join(vec![
            Location::range(6, 8),
            Location::range(0, 2),
        ]))),
        "Origin-spanning reverse feature",
    ));
    let target = Molecule::dna(
        "Declared reverse-rotated coordinate system",
        "CCGGTTTA",
        Topology::Circular,
    )
    .unwrap();
    let request = Request {
        schema_version: 1,
        target_revision: target.revision(),
        source: Source::Donor {
            event: event("circular-independent-donor", &donor),
            molecule: donor,
            transform: Transform {
                origin: 2,
                reverse: true,
            },
            indices: vec![0],
        },
        selections: vec![Selection {
            proposal_id: "donor-0".into(),
            decision: Decision::Append,
        }],
        limits: Limits::default(),
    };
    let output = merge::apply(&target, &request).unwrap();
    let actual = &output.molecule.features[0];
    assert_eq!(actual.location.extract(&target.sequence).unwrap(), "TTTA");
    let path = collagene_model::annotations::biological_path::BiologicalPath::from_location(
        &actual.location,
        8,
        Default::default(),
    )
    .unwrap();
    assert_eq!(
        path.calls()
            .map(|call| (call.position, call.reverse))
            .collect::<Vec<_>>(),
        vec![(4, false), (5, false), (6, false), (7, false)]
    );
    assert_eq!(output.molecule.sequence, target.sequence);
    let mut wrong_axis = request.clone();
    if let Source::Donor { transform, .. } = &mut wrong_axis.source {
        transform.reverse = false;
    }
    assert!(merge::apply(&target, &wrong_axis).is_err());
}

#[test]
fn provenance_addition_preserves_biological_claims_and_keeps_distinct_event_contributors() {
    let (target, mut request) = fixture();
    if let Source::Donor { indices, .. } = &mut request.source {
        *indices = vec![0];
    }
    request.selections = vec![Selection {
        proposal_id: "donor-0".into(),
        decision: Decision::AddProvenance {
            target_index: 0,
            annotation_revision: annotation_revision(&target.features[0]).unwrap(),
        },
    }];
    let first = merge::apply(&target, &request).unwrap().molecule;
    assert_eq!(first.features.len(), 1);
    assert_eq!(first.features[0].location, target.features[0].location);
    assert_eq!(
        first.features[0].needs_review,
        target.features[0].needs_review
    );
    assert!(
        first.features[0]
            .qualifiers
            .contains(&target.features[0].qualifiers[0])
    );
    let ledger = merge::Ledger::from_feature(&first.features[0]).unwrap();
    let encoded_first = serde_json::to_value(&ledger).unwrap();
    let mut independent = request.clone();
    independent.target_revision = first.revision();
    let independent_event = if let Source::Donor {
        event: binding,
        molecule,
        ..
    } = &mut independent.source
    {
        *binding = event("independently-bound-document-B", molecule);
        binding.clone()
    } else {
        unreachable!()
    };
    independent.selections[0].decision = Decision::AddProvenance {
        target_index: 0,
        annotation_revision: annotation_revision(&first.features[0]).unwrap(),
    };
    let second = merge::apply(&first, &independent).unwrap().molecule;
    let encoded_second =
        serde_json::to_value(merge::Ledger::from_feature(&second.features[0]).unwrap()).unwrap();
    assert_ne!(encoded_first, encoded_second);
    let evidence = serde_json::to_string(&encoded_second).unwrap();
    assert!(evidence.contains(&independent_event.document_id));
    if let Source::Donor { event, .. } = request.source {
        assert!(evidence.contains(&event.document_id));
    }
    assert_eq!(
        second.features[0].needs_review,
        target.features[0].needs_review
    );
    assert_eq!(second.sequence, target.sequence);
}

#[test]
fn duplicate_and_conflicting_claims_require_explicit_decisions_not_overlap_deduplication() {
    let (target, mut request) = fixture();
    request.selections[0].decision = Decision::Append;
    assert!(merge::apply(&target, &request).is_err());
    request.selections[0].decision = Decision::Skip;
    if let Source::Donor {
        molecule,
        event: binding,
        ..
    } = &mut request.source
    {
        molecule.features[0].qualifiers.push((
            "note".into(),
            Some("Conflicting source annotation claim".into()),
        ));
        *binding = event("conflicting-source", molecule);
    }
    request.selections[0].decision = Decision::Append;
    let result = merge::apply(&target, &request).unwrap();
    assert_eq!(result.molecule.features[0], target.features[0]);
    assert_eq!(result.molecule.features.len(), 4);
    assert!(result.molecule.features[1].needs_review);
    assert_eq!(
        result.molecule.features[1].location,
        target.features[0].location
    );
    assert_ne!(
        result.molecule.features[1].qualifiers,
        target.features[0].qualifiers
    );
}

#[test]
fn finite_relation_budget_reports_unknown_relationships_and_never_authorizes_absence() {
    let (mut target, mut request) = fixture();
    for index in 0..8 {
        target.features.push(feature(
            Location::range(index, index + 1),
            "Additional independent annotation",
        ));
    }
    request.target_revision = target.revision();
    let mut library = collagene_core::sequence::annotations::library::FeatureLibrary::new(
        "Literal relation budget source",
        "1",
    )
    .unwrap();
    if let Source::Donor { molecule, .. } = &request.source {
        library.add_feature(molecule, 1).unwrap();
    }
    let report = collagene_core::sequence::annotations::library::scan(
        &target,
        &library,
        &Default::default(),
    )
    .unwrap();
    request.selections = vec![Selection {
        proposal_id: report.hits[0].id.clone(),
        decision: Decision::Append,
    }];
    request.limits.max_comparisons = report.comparisons + 1;
    request.source = Source::Library { library, report };
    let preview = merge::preview(&target, &request).unwrap();
    assert!(!preview.complete);
    assert!(
        preview
            .candidates
            .iter()
            .any(|candidate| !candidate.relations_complete)
    );
    assert!(merge::apply(&target, &request).is_err());
    let before = target.clone();
    let mut invalid = request.clone();
    invalid.limits.max_output_bytes = 0;
    assert!(merge::apply(&target, &invalid).is_err());
    assert_eq!(target, before);
    // An early positive match does not make the rest of the selected review complete.
    let (_, mut attribution) = fixture();
    attribution.target_revision = target.revision();
    if let Source::Donor { indices, .. } = &mut attribution.source {
        *indices = vec![0];
    }
    attribution.selections = vec![Selection {
        proposal_id: "donor-0".into(),
        decision: Decision::AddProvenance {
            target_index: 0,
            annotation_revision: annotation_revision(&target.features[0]).unwrap(),
        },
    }];
    attribution.limits.max_comparisons = target.sequence.len() + 3;
    let partial = merge::preview(&target, &attribution).unwrap();
    assert!(!partial.complete);
    assert!(!partial.candidates[0].relations_complete);
    assert!(
        partial.candidates[0]
            .relations
            .iter()
            .any(|relation| matches!(
                relation.kind,
                RelationKind::ExactDuplicate | RelationKind::EquivalentAnnotation
            ))
    );
    assert!(merge::apply(&target, &attribution).is_err());
    let mut invalid = request;
    invalid.selections = vec![invalid.selections[0].clone(); 1001];
    assert!(merge::preview(&target, &invalid).is_err());
}

#[test]
fn explicit_all_skip_library_and_dataset_decisions_survive_history_undo_and_archive_reopen() {
    use collagene_core::documents::workspace::{self, Workspace};
    use collagene_core::sequence::annotations::library::{FeatureLibrary, scan};
    let (target, donor_request) = fixture();
    let mut library = FeatureLibrary::new("Explicit skipped source", "1").unwrap();
    if let Source::Donor { molecule, .. } = donor_request.source {
        library.add_feature(&molecule, 1).unwrap();
    }
    let report = scan(&target, &library, &Default::default()).unwrap();
    let library_id = report.hits[0].id.clone();
    let dataset = collagene_core::formats::annotations::import_annotations(
        &target,
        "ref\t2\t4\tIndependent skipped row\t0\t+\n",
        collagene_core::formats::annotations::AnnotationFormat::Bed,
        "ref",
    )
    .unwrap();
    let dataset_id = dataset.proposals[0].id.clone();
    for (source, proposal_id) in [
        (Source::Library { library, report }, library_id),
        (Source::Dataset { report: dataset }, dataset_id),
    ] {
        let mut history = DocumentHistory::new("explicit-skip-target", target.clone()).unwrap();
        let original_event = history.head_event().to_owned();
        let request = Request {
            schema_version: 1,
            target_revision: target.revision(),
            source,
            selections: vec![Selection {
                proposal_id,
                decision: Decision::Skip,
            }],
            limits: Default::default(),
        };
        history
            .commit(
                vec![Operation::MergeAnnotations { request }],
                "Explicit annotation Skip decision",
            )
            .unwrap();
        let decision_event = history.head_event().to_owned();
        assert_ne!(decision_event, original_event);
        assert_eq!(history.event_count(), 2);
        assert_eq!(history.molecule, target);
        assert_eq!(history.head_revision(), target.revision());
        let mut archive = Workspace {
            name: "Retained unchanged annotation decision".into(),
            active_document: history.id.clone(),
            documents: vec![],
            revisions: Default::default(),
            history: vec![],
            candidates: vec![],
            active_candidate: None,
            attachments: Default::default(),
        };
        history.contribute(&mut archive).unwrap();
        let reopened =
            workspace::load_workspace(&workspace::save_workspace(&archive).unwrap()).unwrap();
        let mut restored = DocumentHistory::restore_all(&reopened).unwrap().remove(0);
        assert_eq!(restored.head_event(), decision_event);
        assert!(restored.undo().unwrap());
        assert_eq!(restored.head_event(), original_event);
        assert_eq!(restored.molecule, target);
        assert!(restored.redo().unwrap());
        assert_eq!(restored.head_event(), decision_event);
        assert_eq!(restored.molecule, target);
    }
}

#[test]
fn fresh_library_rescan_adds_no_second_contribution_for_the_same_stable_source_evidence() {
    use collagene_core::sequence::annotations::library::{FeatureLibrary, scan};
    let (target, donor_request) = fixture();
    let mut library = FeatureLibrary::new("One actual stable source contribution", "1").unwrap();
    if let Source::Donor { molecule, .. } = donor_request.source {
        library.add_feature(&molecule, 0).unwrap();
    }
    let scan = |target: &Molecule| scan(target, &library, &Default::default()).unwrap();
    let first_scan = scan(&target);
    let first_id = first_scan
        .hits
        .iter()
        .find(|hit| hit.start == 0 && !hit.reverse)
        .unwrap()
        .id
        .clone();
    let first_request = Request {
        schema_version: 1,
        target_revision: target.revision(),
        source: Source::Library {
            library: library.clone(),
            report: first_scan,
        },
        selections: vec![Selection {
            proposal_id: first_id.clone(),
            decision: Decision::AddProvenance {
                target_index: 0,
                annotation_revision: annotation_revision(&target.features[0]).unwrap(),
            },
        }],
        limits: Default::default(),
    };
    let first = merge::apply(&target, &first_request).unwrap();
    assert!(first.receipt.actions[0].contribution_added);
    assert_ne!(first.molecule.revision(), target.revision());
    let original_ledger = merge::Ledger::from_feature(&first.molecule.features[0]).unwrap();
    assert_eq!(original_ledger.contributions.len(), 1);
    let next_scan = scan(&first.molecule);
    let next_id = next_scan
        .hits
        .iter()
        .find(|hit| hit.start == 0 && !hit.reverse)
        .unwrap()
        .id
        .clone();
    assert_ne!(
        first_id, next_id,
        "actual target-bound operation IDs change after attribution"
    );
    let next_request = Request {
        schema_version: 1,
        target_revision: first.molecule.revision(),
        source: Source::Library {
            library,
            report: next_scan,
        },
        selections: vec![Selection {
            proposal_id: next_id.clone(),
            decision: Decision::AddProvenance {
                target_index: 0,
                annotation_revision: annotation_revision(&first.molecule.features[0]).unwrap(),
            },
        }],
        limits: Default::default(),
    };
    let second = merge::apply(&first.molecule, &next_request).unwrap();
    assert!(!second.receipt.actions[0].contribution_added);
    assert_eq!(second.molecule, first.molecule);
    assert_eq!(
        second.molecule.features[0].needs_review,
        target.features[0].needs_review
    );
    assert_eq!(
        merge::Ledger::from_feature(&second.molecule.features[0]).unwrap(),
        original_ledger
    );
    let merge::SourcePin::Library { proposal_id, .. } = &original_ledger.contributions[0].source
    else {
        panic!("not library evidence")
    };
    assert_eq!(
        proposal_id, &first_id,
        "retain the original documentary operation pin"
    );
    let mut duplicate = original_ledger.clone();
    let mut repeated = duplicate.contributions[0].clone();
    if let merge::SourcePin::Library { proposal_id, .. } = &mut repeated.source {
        *proposal_id = next_id;
    }
    duplicate.contributions.push(repeated);
    assert!(
        duplicate.validate().is_err(),
        "transient proposal ID cannot conceal duplicate stable evidence"
    );
    let mut skipped = next_request;
    skipped.selections[0].decision = Decision::Skip;
    let skipped = merge::apply(&first.molecule, &skipped).unwrap();
    assert_eq!(skipped.molecule, first.molecule);
    assert!(!skipped.receipt.actions[0].contribution_added);
    assert_eq!(
        first.review.context_revision,
        merge::revision(&target, &first_request).unwrap()
    );
}

#[test]
fn fresh_same_file_reimport_preserves_one_dataset_contribution_and_its_first_operation_pin() {
    use collagene_core::formats::annotations::{AnnotationFormat, import_annotations};
    const FILE: &str = "ref\t0\t2\tIndependent supplied annotation\t0\t+\n";
    let mut target =
        Molecule::dna("Actual dataset target", "AACCGGTTAA", Topology::Linear).unwrap();
    let initial_import = import_annotations(&target, FILE, AnnotationFormat::Bed, "ref").unwrap();
    // An explicitly present annotation supplies the role; AddProvenance changes
    // attribution only, never its biological qualifiers or review state.
    target
        .features
        .push(initial_import.proposals[0].feature.clone());
    let first_import = import_annotations(&target, FILE, AnnotationFormat::Bed, "ref").unwrap();
    let first_id = first_import.proposals[0].id.clone();
    let first_request = Request {
        schema_version: 1,
        target_revision: target.revision(),
        source: Source::Dataset {
            report: first_import,
        },
        selections: vec![Selection {
            proposal_id: first_id.clone(),
            decision: Decision::AddProvenance {
                target_index: 0,
                annotation_revision: annotation_revision(&target.features[0]).unwrap(),
            },
        }],
        limits: Default::default(),
    };
    let first = merge::apply(&target, &first_request).unwrap();
    assert!(first.receipt.actions[0].contribution_added);
    let original_ledger = merge::Ledger::from_feature(&first.molecule.features[0]).unwrap();
    assert_eq!(original_ledger.contributions.len(), 1);
    let next_import =
        import_annotations(&first.molecule, FILE, AnnotationFormat::Bed, "ref").unwrap();
    let next_id = next_import.proposals[0].id.clone();
    assert_eq!(
        first_id, next_id,
        "dataset proposal IDs derive from source bytes/rows/content and are already target-independent"
    );
    let next_request = Request {
        schema_version: 1,
        target_revision: first.molecule.revision(),
        source: Source::Dataset {
            report: next_import,
        },
        selections: vec![Selection {
            proposal_id: next_id.clone(),
            decision: Decision::AddProvenance {
                target_index: 0,
                annotation_revision: annotation_revision(&first.molecule.features[0]).unwrap(),
            },
        }],
        limits: Default::default(),
    };
    let second = merge::apply(&first.molecule, &next_request).unwrap();
    assert!(!second.receipt.actions[0].contribution_added);
    assert_eq!(second.molecule, first.molecule);
    assert_eq!(
        second.molecule.features[0].needs_review,
        target.features[0].needs_review
    );
    for qualifier in &target.features[0].qualifiers {
        assert!(second.molecule.features[0].qualifiers.contains(qualifier));
    }
    assert_eq!(
        merge::Ledger::from_feature(&second.molecule.features[0]).unwrap(),
        original_ledger
    );
    let merge::SourcePin::Dataset { proposal_id, .. } = &original_ledger.contributions[0].source
    else {
        panic!("not dataset evidence")
    };
    assert_eq!(proposal_id, &first_id);
    let mut duplicate = original_ledger.clone();
    let mut repeated = duplicate.contributions[0].clone();
    if let merge::SourcePin::Dataset { proposal_id, .. } = &mut repeated.source {
        *proposal_id = "different-untrusted-operation-id".into();
    }
    duplicate.contributions.push(repeated);
    assert!(duplicate.validate().is_err());
    let mut skipped = next_request;
    skipped.selections[0].decision = Decision::Skip;
    let skipped = merge::apply(&first.molecule, &skipped).unwrap();
    assert_eq!(skipped.molecule, first.molecule);
    assert!(!skipped.receipt.actions[0].contribution_added);
    assert_eq!(
        first.review.context_revision,
        merge::revision(&target, &first_request).unwrap()
    );
}
