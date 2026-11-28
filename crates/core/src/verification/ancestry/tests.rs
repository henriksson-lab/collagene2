use super::*;
use collagene_model::{Feature, Location, Topology, annotations};

fn molecule(dna: &str) -> Molecule {
    Molecule::dna("literal native snapshot", dna, Topology::Linear).unwrap()
}
fn segment(start: usize, end: usize, destination: usize, reverse: bool) -> MappingSegment {
    MappingSegment {
        source_start: start,
        source_end: end,
        destination_start: destination,
        reverse,
    }
}
fn project(source: &Molecule, destination: &Molecule, map: &[MappingSegment]) -> Projection {
    annotations::project_features_with_topology(
        source,
        map,
        destination.sequence.len(),
        destination.topology,
    )
    .unwrap()
}
#[test]
fn identical_dna_and_labels_keep_exact_input_role_and_feature_ordinals() {
    let mut source = molecule("ACGT");
    source.features.push(Feature {
        kind: "CDS".into(),
        location: Location::range(0, 4),
        qualifiers: vec![("label".into(), Some("same label".into()))],
        needs_review: false,
    });
    let mut builder = TraceBuilder::new(1, Limits::default()).unwrap();
    let a = builder.register_input("a", &source).unwrap();
    let b = builder.register_input("b", &source).unwrap();
    assert_ne!(a, b);
    let mut step = builder.begin_step(0).unwrap();
    let mut product = molecule("ACGTACGT");
    let output = step.destination(&product).unwrap();
    let first = [segment(0, 4, 0, false)];
    let second = [segment(0, 4, 4, false)];
    let pa = project(&source, &product, &first);
    step.project(ProjectionCapture {
        source: a,
        source_molecule: &source,
        destination: output,
        destination_molecule: &product,
        occurrence: 0,
        kind: ProjectionKind::ReferenceCorrespondence,
        segments: &first,
        projection: &pa,
        destination_feature_indices: &[0],
        changes: &[],
    })
    .unwrap();
    product.features.extend(pa.features);
    let pb = project(&source, &product, &second);
    step.project(ProjectionCapture {
        source: b,
        source_molecule: &source,
        destination: output,
        destination_molecule: &product,
        occurrence: 1,
        kind: ProjectionKind::ReferenceCorrespondence,
        segments: &second,
        projection: &pb,
        destination_feature_indices: &[1],
        changes: &[],
    })
    .unwrap();
    product.features.extend(pb.features);
    step.seal(output, &product).unwrap();
    step.commit(output).unwrap();
    let trace = builder.finish();
    assert!(trace.step(0).unwrap().local_coverage_complete);
    let mut lookup = trace.lookup_session();
    let first = lookup.at(0, 0, None).unwrap();
    assert!(first.complete);
    assert_eq!(first.origins.len(), 1);
    assert_eq!(
        first.origins[0].source.source,
        SourceHandle::Plan(Source::Input("a".into()))
    );
    let second = lookup.at(0, 4, None).unwrap();
    assert_eq!(
        second.origins[0].source.source,
        SourceHandle::Plan(Source::Input("b".into()))
    );
    let feature = lookup.feature_at(0, 1, None).unwrap();
    assert!(feature.complete);
    assert_eq!(
        feature.origins[0].source.source,
        SourceHandle::Plan(Source::Input("b".into()))
    );
    assert_eq!(feature.origins[0].feature_index, 0);
}
#[test]
fn explicit_intermediate_tokens_compose_reverse_positions_without_content_aliases() {
    let source = molecule("ACGT"); // reverse complement is deliberately identical.
    let mut builder = TraceBuilder::new(1, Limits::default()).unwrap();
    let input = builder.register_input("original", &source).unwrap();
    let mut step = builder.begin_step(0).unwrap();
    let output_molecule = molecule("ACGT");
    // Destination creation can precede another actual source stage. Only sealing
    // and use order supplies the dependency; token equality never comes from DNA.
    let output = step.destination(&output_molecule).unwrap();
    let intermediate = step.destination(&output_molecule).unwrap();
    let reverse = [segment(0, 4, 0, true)];
    let projection = project(&source, &output_molecule, &reverse);
    step.project(ProjectionCapture {
        source: input,
        source_molecule: &source,
        destination: intermediate,
        destination_molecule: &output_molecule,
        occurrence: 0,
        kind: ProjectionKind::TemplateSynthesis,
        segments: &reverse,
        projection: &projection,
        destination_feature_indices: &[],
        changes: &[],
    })
    .unwrap();
    step.seal(intermediate, &output_molecule).unwrap();
    let identity = [segment(0, 4, 0, false)];
    let projection = project(&output_molecule, &output_molecule, &identity);
    step.project(ProjectionCapture {
        source: intermediate,
        source_molecule: &output_molecule,
        destination: output,
        destination_molecule: &output_molecule,
        occurrence: 1,
        kind: ProjectionKind::ReferenceCorrespondence,
        segments: &identity,
        projection: &projection,
        destination_feature_indices: &[],
        changes: &[],
    })
    .unwrap();
    step.seal(output, &output_molecule).unwrap();
    step.commit(output).unwrap();
    let trace = builder.finish();
    let mut lookup = trace.lookup_session();
    let result = lookup.at(0, 0, None).unwrap();
    assert!(result.complete);
    assert_eq!(result.origins[0].position, 3);
    assert!(result.origins[0].reverse);
    assert_eq!(result.origins[0].path.len(), 2);
}
#[test]
fn discarded_operation_tokens_never_rebind_when_the_same_slot_is_reused() {
    let source = molecule("ACGT");
    let mut builder = TraceBuilder::new(2, Limits::default()).unwrap();
    let input = builder.register_input("original", &source).unwrap();
    let stale = {
        let mut step = builder.begin_step(0).unwrap();
        step.destination(&source).unwrap()
    };
    let mut step = builder.begin_step(1).unwrap();
    let current = step.destination(&source).unwrap();
    assert_ne!(stale, current);
    let identity = [segment(0, 4, 0, false)];
    let projection = project(&source, &source, &identity);
    assert!(
        step.project(ProjectionCapture {
            source: stale,
            source_molecule: &source,
            destination: current,
            destination_molecule: &source,
            occurrence: 0,
            kind: ProjectionKind::ReferenceCorrespondence,
            segments: &identity,
            projection: &projection,
            destination_feature_indices: &[],
            changes: &[]
        })
        .is_err()
    );
    step.project(ProjectionCapture {
        source: input,
        source_molecule: &source,
        destination: current,
        destination_molecule: &source,
        occurrence: 0,
        kind: ProjectionKind::ReferenceCorrespondence,
        segments: &identity,
        projection: &projection,
        destination_feature_indices: &[],
        changes: &[],
    })
    .unwrap();
    step.seal(current, &source).unwrap();
    step.commit(current).unwrap();
    let trace = builder.finish();
    assert!(trace.step(0).is_none());
    assert_eq!(trace.step(1).unwrap().step_index, 1);
}
#[test]
fn missing_correspondence_is_incomplete_and_not_false_origin_absence() {
    let source = molecule("ACGT");
    let mut builder = TraceBuilder::new(1, Limits::default()).unwrap();
    let mut step = builder.begin_step(0).unwrap();
    let output = step.destination(&source).unwrap();
    step.seal(output, &source).unwrap();
    step.commit(output).unwrap();
    let trace = builder.finish();
    assert!(!trace.step(0).unwrap().local_coverage_complete);
    let mut lookup = trace.lookup_session();
    let result = lookup.at(0, 0, None).unwrap();
    assert!(!result.complete);
    assert!(result.origins.is_empty());
    assert!(result.reason.is_some());
}
#[test]
fn explicit_inserted_literal_is_not_fabricated_template_or_chemical_ancestry() {
    let product = molecule("AAA");
    let mut builder = TraceBuilder::new(1, Limits::default()).unwrap();
    let mut step = builder.begin_step(0).unwrap();
    let output = step.destination(&product).unwrap();
    step.introduce(output, 0, 3, IntroductionKind::DeclaredLiteral)
        .unwrap();
    step.seal(output, &product).unwrap();
    step.commit(output).unwrap();
    let trace = builder.finish();
    let mut lookup = trace.lookup_session();
    let result = lookup.at(0, 1, None).unwrap();
    assert!(result.complete);
    assert!(result.origins.is_empty());
    assert_eq!(result.introductions.len(), 1);
    assert_eq!(
        *result.introductions[0].kind,
        IntroductionKind::DeclaredLiteral
    );
}
#[test]
fn lookup_work_is_shared_and_failed_stage_capacity_is_not_reset() {
    let product = molecule("AAA");
    let mut limits = Limits::default();
    limits.max_lookup_work = 1;
    let mut builder = TraceBuilder::new(1, limits).unwrap();
    let mut step = builder.begin_step(0).unwrap();
    let output = step.destination(&product).unwrap();
    step.introduce(output, 0, 3, IntroductionKind::DeclaredLiteral)
        .unwrap();
    step.seal(output, &product).unwrap();
    step.commit(output).unwrap();
    let trace = builder.finish();
    let mut lookup = trace.lookup_session();
    assert!(!lookup.at(0, 0, None).unwrap().complete);
    assert_eq!(lookup.remaining(), 0);
    assert!(!lookup.at(0, 1, None).unwrap().complete);
    let mut limits = Limits::default();
    limits.max_snapshots = 1;
    let mut builder = TraceBuilder::new(2, limits).unwrap();
    {
        let mut first = builder.begin_step(0).unwrap();
        first.destination(&product).unwrap();
    }
    assert!(
        builder
            .begin_step(1)
            .unwrap()
            .destination(&product)
            .is_err()
    );
}

#[test]
fn an_unsealed_or_changed_native_source_cannot_supply_a_projection() {
    let source = molecule("ACGT");
    let mut builder = TraceBuilder::new(1, Limits::default()).unwrap();
    let input = builder.register_input("original", &source).unwrap();
    let mut step = builder.begin_step(0).unwrap();
    let draft = step.destination(&source).unwrap();
    let output = step.destination(&source).unwrap();
    let identity = [segment(0, 4, 0, false)];
    let projection = project(&source, &source, &identity);
    assert!(
        step.project(ProjectionCapture {
            source: draft,
            source_molecule: &source,
            destination: output,
            destination_molecule: &source,
            occurrence: 0,
            kind: ProjectionKind::ReferenceCorrespondence,
            segments: &identity,
            projection: &projection,
            destination_feature_indices: &[],
            changes: &[]
        })
        .is_err()
    );
    let changed = molecule("TCGT");
    assert!(
        step.project(ProjectionCapture {
            source: input,
            source_molecule: &changed,
            destination: output,
            destination_molecule: &source,
            occurrence: 0,
            kind: ProjectionKind::ReferenceCorrespondence,
            segments: &identity,
            projection: &projection,
            destination_feature_indices: &[],
            changes: &[]
        })
        .is_err()
    );
}
#[test]
fn mapped_substitutions_require_the_actual_native_change_coordinates_and_calls() {
    let source = molecule("ACGT");
    let product = molecule("ATGT");
    let mut builder = TraceBuilder::new(1, Limits::default()).unwrap();
    let input = builder.register_input("original", &source).unwrap();
    let mut step = builder.begin_step(0).unwrap();
    let output = step.destination(&product).unwrap();
    let identity = [segment(0, 4, 0, false)];
    let projection = project(&source, &product, &identity);
    assert!(
        step.project(ProjectionCapture {
            source: input,
            source_molecule: &source,
            destination: output,
            destination_molecule: &product,
            occurrence: 0,
            kind: ProjectionKind::ReferenceCorrespondence,
            segments: &identity,
            projection: &projection,
            destination_feature_indices: &[],
            changes: &[]
        })
        .is_err()
    );
    let changes = [MappedChange {
        source_position: 1,
        destination_position: 1,
        before: 'C',
        after: 'T',
    }];
    step.project(ProjectionCapture {
        source: input,
        source_molecule: &source,
        destination: output,
        destination_molecule: &product,
        occurrence: 0,
        kind: ProjectionKind::EditedTemplate,
        segments: &identity,
        projection: &projection,
        destination_feature_indices: &[],
        changes: &changes,
    })
    .unwrap();
    step.seal(output, &product).unwrap();
    step.commit(output).unwrap();
    let trace = builder.finish();
    let mut lookup = trace.lookup_session();
    let changed = lookup.at(0, 1, None).unwrap();
    assert!(changed.complete);
    assert_eq!(changed.origins[0].position, 1);
    assert_eq!(changed.origins[0].changes.len(), 1);
    assert_eq!(changed.origins[0].changes[0].change.after, 'T');
}

#[test]
fn a_foreign_trace_token_cannot_alias_an_equal_local_snapshot() {
    let source = molecule("ACGT");
    let mut foreign = TraceBuilder::new(1, Limits::default()).unwrap();
    let foreign_input = foreign.register_input("different role", &source).unwrap();
    let mut builder = TraceBuilder::new(1, Limits::default()).unwrap();
    let input = builder.register_input("local role", &source).unwrap();
    assert_ne!(input, foreign_input);
    let mut step = builder.begin_step(0).unwrap();
    let output = step.destination(&source).unwrap();
    let identity = [segment(0, 4, 0, false)];
    let projection = project(&source, &source, &identity);
    assert!(
        step.project(ProjectionCapture {
            source: foreign_input,
            source_molecule: &source,
            destination: output,
            destination_molecule: &source,
            occurrence: 0,
            kind: ProjectionKind::ReferenceCorrespondence,
            segments: &identity,
            projection: &projection,
            destination_feature_indices: &[],
            changes: &[],
        })
        .is_err()
    );
}

#[test]
fn a_homology_collapse_retains_both_actual_reference_origins() {
    let source = molecule("ACGT");
    let mut builder = TraceBuilder::new(1, Limits::default()).unwrap();
    let left = builder.register_input("left role", &source).unwrap();
    let right = builder.register_input("right role", &source).unwrap();
    let mut step = builder.begin_step(0).unwrap();
    let output = step.destination(&source).unwrap();
    let identity = [segment(0, 4, 0, false)];
    let projection = project(&source, &source, &identity);
    for (occurrence, token) in [left, right].into_iter().enumerate() {
        step.project(ProjectionCapture {
            source: token,
            source_molecule: &source,
            destination: output,
            destination_molecule: &source,
            occurrence: occurrence as u32,
            kind: ProjectionKind::HomologyCollapse,
            segments: &identity,
            projection: &projection,
            destination_feature_indices: &[],
            changes: &[],
        })
        .unwrap();
    }
    step.seal(output, &source).unwrap();
    step.commit(output).unwrap();
    let trace = builder.finish();
    let mut lookup = trace.lookup_session();
    let call = lookup.at(0, 2, None).unwrap();
    assert!(call.complete);
    assert_eq!(call.origins.len(), 2);
    let roles = call
        .origins
        .iter()
        .map(|origin| origin.source.source.clone())
        .collect::<Vec<_>>();
    assert!(roles.contains(&SourceHandle::Plan(Source::Input("left role".into()))));
    assert!(roles.contains(&SourceHandle::Plan(Source::Input("right role".into()))));
    assert!(call.origins.iter().all(|origin| origin.position == 2));
}

#[test]
fn an_appended_annotation_change_cannot_preserve_the_projection_pin() {
    let mut source = molecule("ACGT");
    source.features.push(Feature {
        kind: "CDS".into(),
        location: Location::range(0, 4),
        qualifiers: vec![],
        needs_review: false,
    });
    let mut builder = TraceBuilder::new(1, Limits::default()).unwrap();
    let input = builder.register_input("original", &source).unwrap();
    let mut step = builder.begin_step(0).unwrap();
    let mut product = molecule("ACGT");
    let output = step.destination(&product).unwrap();
    let identity = [segment(0, 4, 0, false)];
    let projection = project(&source, &product, &identity);
    step.project(ProjectionCapture {
        source: input,
        source_molecule: &source,
        destination: output,
        destination_molecule: &product,
        occurrence: 0,
        kind: ProjectionKind::ReferenceCorrespondence,
        segments: &identity,
        projection: &projection,
        destination_feature_indices: &[0],
        changes: &[],
    })
    .unwrap();
    product.features.extend(projection.features);
    product.features[0]
        .qualifiers
        .push(("note".into(), Some("edited after capture".into())));
    assert!(step.seal(output, &product).is_err());
    drop(step);
    assert!(builder.finish().step(0).is_none());
}

#[test]
fn a_forged_original_feature_ordinal_is_rejected_without_label_matching() {
    let mut source = molecule("ACGT");
    source.features.push(Feature {
        kind: "CDS".into(),
        location: Location::range(0, 4),
        qualifiers: vec![("label".into(), Some("plausible imported label".into()))],
        needs_review: true,
    });
    let mut builder = TraceBuilder::new(1, Limits::default()).unwrap();
    let input = builder.register_input("original", &source).unwrap();
    let mut step = builder.begin_step(0).unwrap();
    let output = step.destination(&source).unwrap();
    let identity = [segment(0, 4, 0, false)];
    let mut projection = project(&source, &source, &identity);
    projection.source_feature_indices[0] = 1;
    assert!(
        step.project(ProjectionCapture {
            source: input,
            source_molecule: &source,
            destination: output,
            destination_molecule: &source,
            occurrence: 0,
            kind: ProjectionKind::ReferenceCorrespondence,
            segments: &identity,
            projection: &projection,
            destination_feature_indices: &[1],
            changes: &[],
        })
        .is_err()
    );
    drop(step);
    assert!(builder.finish().step(0).is_none());
}

#[test]
fn same_executor_pcr_trace_preserves_report_and_literal_primer_tail_origins() {
    use crate::{
        chemistry::Strand,
        planning::plan::{DesignPlan, Operation, PlanInput},
        sequence::analysis::Primer,
    };
    let source = molecule("AGTCCGATGCTACGTAGGACTTCAGACCTGTA");
    let target = molecule("GGAGTCCGATGCTACGTAGGACTTCAGACCTGTACC");
    let plan = DesignPlan {
        standard: None,
        name: "native receipt oracle".into(),
        target_revision: target.revision(),
        inputs: vec![PlanInput {
            id: "template".into(),
            revision: source.revision(),
            molecule: source,
        }],
        steps: vec![Operation::Pcr {
            template: Source::Input("template".into()),
            forward: Primer::new("forward", "AGTCCGAT", "GG").unwrap(),
            reverse: Primer::new("reverse", "TACAGGTC", "GG").unwrap(),
            template_strand: Strand::Top,
        }],
    };
    let ordinary = crate::verification::engine::verify(&plan, &target);
    let (captured, trace) = crate::verification::engine::verify_with_ancestry(&plan, &target);
    assert_eq!(
        serde_json::to_value(&ordinary).unwrap(),
        serde_json::to_value(&captured).unwrap()
    );
    assert_eq!(
        captured.product(0).unwrap().sequence.as_str(),
        target.sequence.as_str()
    );
    let trace = trace.unwrap();
    let mut lookup = trace.lookup_session();
    for position in [0, 1, 34, 35] {
        let call = lookup.at(0, position, None).unwrap();
        assert!(call.complete);
        assert!(call.origins.is_empty());
        assert_eq!(call.introductions.len(), 1);
        assert_eq!(*call.introductions[0].kind, IntroductionKind::PrimerTail);
    }
    let call = lookup.at(0, 2, None).unwrap();
    assert!(call.complete);
    assert_eq!(call.origins[0].position, 0);
    assert_eq!(
        call.origins[0].source.source,
        SourceHandle::Plan(Source::Input("template".into()))
    );
}

#[test]
fn same_executor_overlap_retains_both_junction_sources_and_exact_feature_ordinals() {
    use crate::planning::plan::{DesignPlan, Operation, PlanInput};
    let mut left = molecule("AAAACCCC");
    let mut right = molecule("CCCCGGGG");
    for source in [&mut left, &mut right] {
        source.features.push(Feature {
            kind: "misc_feature".into(),
            location: Location::range(4, 8),
            qualifiers: vec![],
            needs_review: false,
        });
    }
    let target = molecule("AAAACCCCGGGG");
    let plan = DesignPlan {
        standard: None,
        name: "native overlap receipt oracle".into(),
        target_revision: target.revision(),
        inputs: vec![
            PlanInput {
                id: "left".into(),
                revision: left.revision(),
                molecule: left,
            },
            PlanInput {
                id: "right".into(),
                revision: right.revision(),
                molecule: right,
            },
        ],
        steps: vec![Operation::Overlap {
            fragments: vec![Source::Input("left".into()), Source::Input("right".into())],
            overlaps: vec![4],
            circular: false,
        }],
    };
    let (report, trace) = crate::verification::engine::verify_with_ancestry(&plan, &target);
    assert_eq!(report.product(0).unwrap().sequence.as_str(), "AAAACCCCGGGG");
    let trace = trace.unwrap();
    let mut lookup = trace.lookup_session();
    let call = lookup.at(0, 4, None).unwrap();
    assert!(call.complete);
    assert_eq!(call.origins.len(), 2);
    assert!(call.origins.iter().any(|origin| origin.source.source
        == SourceHandle::Plan(Source::Input("left".into()))
        && origin.position == 4));
    assert!(call.origins.iter().any(|origin| origin.source.source
        == SourceHandle::Plan(Source::Input("right".into()))
        && origin.position == 0));
    let feature = lookup.feature_at(0, 1, None).unwrap();
    assert!(feature.complete);
    assert_eq!(
        feature.origins[0].source.source,
        SourceHandle::Plan(Source::Input("right".into()))
    );
    assert_eq!(feature.origins[0].feature_index, 0);
}

#[test]
fn enclosing_execution_withholding_can_discard_the_latest_committed_trace_only() {
    let source = molecule("ACGT");
    let mut builder = TraceBuilder::new(2, Limits::default()).unwrap();
    let input = builder.register_input("original", &source).unwrap();
    let identity = [segment(0, 4, 0, false)];
    let projection = project(&source, &source, &identity);
    let mut step = builder.begin_step(0).unwrap();
    let withheld = step.destination(&source).unwrap();
    step.project(ProjectionCapture {
        source: input,
        source_molecule: &source,
        destination: withheld,
        destination_molecule: &source,
        occurrence: 0,
        kind: ProjectionKind::ReferenceCorrespondence,
        segments: &identity,
        projection: &projection,
        destination_feature_indices: &[],
        changes: &[],
    })
    .unwrap();
    step.seal(withheld, &source).unwrap();
    step.commit(withheld).unwrap();
    builder.discard_step(0).unwrap();
    let mut step = builder.begin_step(1).unwrap();
    let actual = step.destination(&source).unwrap();
    assert_ne!(withheld, actual);
    step.project(ProjectionCapture {
        source: input,
        source_molecule: &source,
        destination: actual,
        destination_molecule: &source,
        occurrence: 0,
        kind: ProjectionKind::ReferenceCorrespondence,
        segments: &identity,
        projection: &projection,
        destination_feature_indices: &[],
        changes: &[],
    })
    .unwrap();
    step.seal(actual, &source).unwrap();
    step.commit(actual).unwrap();
    assert!(builder.discard_step(0).is_err());
    let trace = builder.finish();
    assert!(trace.step(0).is_none());
    assert!(trace.step(1).is_some());
    assert!(trace.snapshot_source(withheld).is_none());
}

#[test]
fn actual_annotation_dedup_retains_both_declared_source_feature_origins() {
    use crate::planning::plan::{DesignPlan, Operation, PlanInput};
    let mut left = molecule("AAAACCCC");
    let mut right = molecule("CCCCGGGG");
    left.features.push(Feature {
        kind: "misc_feature".into(),
        location: Location::range(4, 8),
        qualifiers: vec![],
        needs_review: false,
    });
    right.features.push(Feature {
        kind: "misc_feature".into(),
        location: Location::range(0, 4),
        qualifiers: vec![],
        needs_review: false,
    });
    let target = molecule("AAAACCCCGGGG");
    let plan = DesignPlan {
        standard: None,
        name: "actual dedup receipt oracle".into(),
        target_revision: target.revision(),
        inputs: vec![
            PlanInput {
                id: "left".into(),
                revision: left.revision(),
                molecule: left,
            },
            PlanInput {
                id: "right".into(),
                revision: right.revision(),
                molecule: right,
            },
        ],
        steps: vec![Operation::Overlap {
            fragments: vec![Source::Input("left".into()), Source::Input("right".into())],
            overlaps: vec![4],
            circular: false,
        }],
    };
    let (report, trace) = crate::verification::engine::verify_with_ancestry(&plan, &target);
    assert_eq!(report.product(0).unwrap().features.len(), 1);
    let trace = trace.unwrap();
    let mut lookup = trace.lookup_session();
    let feature = lookup.feature_at(0, 0, None).unwrap();
    assert!(feature.complete);
    assert_eq!(feature.origins.len(), 2);
    assert!(feature.origins.iter().any(|origin| origin.source.source
        == SourceHandle::Plan(Source::Input("left".into()))
        && origin.feature_index == 0));
    assert!(feature.origins.iter().any(|origin| origin.source.source
        == SourceHandle::Plan(Source::Input("right".into()))
        && origin.feature_index == 0));
}

#[test]
fn retained_runs_with_identical_original_handles_keep_exact_native_namespaces() {
    use crate::material::{
        lineage::{EventBinding, InputEventBinding},
        run::Selection,
    };
    let mut source = molecule("ACGT");
    source.features.push(Feature {
        kind: "misc_feature".into(),
        location: Location::range(0, 4),
        qualifiers: vec![],
        needs_review: false,
    });
    let budget = TraceBudget::new(Limits::default()).unwrap();
    let child_pin = "a".repeat(64);
    let make_child = || {
        let mut builder = TraceBuilder::with_budget(1, budget.clone()).unwrap();
        builder.bind_plan(&child_pin).unwrap();
        let original = builder.register_input("same original id", &source).unwrap();
        let mut step = builder.begin_step(0).unwrap();
        let output = step.destination(&source).unwrap();
        let map = [segment(0, 4, 0, false)];
        let projection = project(&source, &source, &map);
        step.project(ProjectionCapture {
            source: original,
            source_molecule: &source,
            destination: output,
            destination_molecule: &source,
            occurrence: 0,
            kind: ProjectionKind::ReferenceCorrespondence,
            segments: &map,
            projection: &projection,
            destination_feature_indices: &[0],
            changes: &[],
        })
        .unwrap();
        step.seal(output, &source).unwrap();
        step.commit(output).unwrap();
        builder.finish()
    };
    let first = make_child();
    let second = make_child();
    let binding = |run_id: &str| {
        let event = EventBinding {
            document_id: format!("{run_id}-document"),
            event_id: format!("{run_id}-event"),
            revision: source.revision(),
        };
        OriginBinding::Run {
            selection: Selection {
                run_id: run_id.into(),
                step_index: 0,
                event: event.clone(),
            },
            plan_revision: child_pin.clone(),
            target_event: event.clone(),
            products: vec![crate::material::lineage::ProductBinding {
                step_index: 0,
                document_id: event.document_id.clone(),
                event_id: event.event_id.clone(),
                revision: event.revision.clone(),
            }],
            inputs: vec![InputEventBinding {
                input_id: "same original id".into(),
                event,
            }],
        }
    };
    let mut builder = TraceBuilder::with_budget(1, budget).unwrap();
    builder.bind_plan(&"b".repeat(64)).unwrap();
    let mut step = builder.begin_step(0).unwrap();
    let a = step.import(binding("retained A"), first, 0, 0).unwrap();
    let b = step.import(binding("retained B"), second, 0, 1).unwrap();
    let mut product = molecule("ACGTACGT");
    let output = step.destination(&product).unwrap();
    for (occurrence, token) in [a, b].into_iter().enumerate() {
        let map = [segment(0, 4, occurrence * 4, false)];
        let projection = project(&source, &product, &map);
        step.project(ProjectionCapture {
            source: token,
            source_molecule: &source,
            destination: output,
            destination_molecule: &product,
            occurrence: occurrence as u32,
            kind: ProjectionKind::ReferenceCorrespondence,
            segments: &map,
            projection: &projection,
            destination_feature_indices: &[occurrence],
            changes: &[],
        })
        .unwrap();
        product.features.extend(projection.features);
    }
    step.seal(output, &product).unwrap();
    step.commit(output).unwrap();
    let trace = builder.finish();
    let baseline = trace.sources().find(|(_, identity)| {
        identity.source == SourceHandle::Plan(Source::Step(0))
            && identity.scope.first().is_some_and(|scope| matches!(&*scope.binding, OriginBinding::Run { selection, .. } if selection.run_id == "retained A"))
    }).map(|(token, _)| token).unwrap();
    let mut lookup = trace.lookup_session();
    let stopped = lookup
        .at(0, 0, Some(LookupStop::Snapshot(baseline)))
        .unwrap();
    assert!(stopped.complete);
    assert_eq!(
        stopped.origins[0].source.source,
        SourceHandle::Plan(Source::Step(0))
    );
    assert_eq!(stopped.origins[0].path.len(), 1);
    assert!(
        lookup
            .at(
                0,
                0,
                Some(LookupStop::Snapshot(
                    trace.step(0).unwrap().final_snapshot()
                ))
            )
            .is_err()
    );
    for (position, run_id) in [(0, "retained A"), (4, "retained B")] {
        let evidence = lookup.at(0, position, None).unwrap();
        assert!(evidence.complete);
        assert_eq!(evidence.origins.len(), 1);
        let origin = &evidence.origins[0];
        assert_eq!(
            origin.source.source,
            SourceHandle::Plan(Source::Input("same original id".into()))
        );
        assert_eq!(origin.source.scope.len(), 1);
        let OriginBinding::Run { selection, .. } = &*origin.source.scope[0].binding else {
            panic!("run scope required")
        };
        assert_eq!(selection.run_id, run_id);
        assert!(origin.path.iter().any(|edge| !edge.scope.is_empty()));
    }
    let evidence = lookup.feature_at(0, 1, None).unwrap();
    assert!(evidence.complete);
    assert_eq!(evidence.origins[0].feature_index, 0);
    let OriginBinding::Run { selection, .. } = &*evidence.origins[0].source.scope[0].binding else {
        panic!("run feature scope required")
    };
    assert_eq!(selection.run_id, "retained B");
}

#[test]
fn material_receipts_name_actual_recipe_roles_and_share_nested_construction_limits() {
    use crate::material::derivation::SourceRole;
    let source = molecule("ACGT");
    let limits = Limits {
        max_snapshots: 1,
        ..Limits::default()
    };
    let budget = TraceBudget::new(limits).unwrap();
    let mut builder = TraceBuilder::for_material(budget.clone(), &"c".repeat(64)).unwrap();
    let template = builder
        .register_material_source(SourceRole::Template, &source)
        .unwrap();
    assert!(builder.register_input("fake plan input", &source).is_err());
    let mut sibling = TraceBuilder::with_budget(1, budget).unwrap();
    assert!(sibling.register_input("other replay", &source).is_err());
    let trace = builder.finish();
    assert_eq!(
        trace.snapshot_source(template).unwrap().source,
        SourceHandle::Material(SourceRole::Template)
    );
    assert!(
        matches!(trace.context(), TraceContext::Material { context_revision } if context_revision == &"c".repeat(64))
    );
}

#[test]
fn retained_graph_root_roles_keep_actual_scope_and_missing_event_provenance() {
    use crate::{
        cloning::topo::graph_history::Selection,
        material::lineage::EventBinding,
        verification::projection::{
            GraphProjectionSource, GraphProjections, NativeProjection, ReferenceProjection,
            SourceProjection, TopoScope,
        },
    };
    let source = molecule("ACGT");
    let product = molecule("ACGTACGT");
    let make_source = |record: char, event| {
        let binding = OriginBinding::TopoGraph {
            selection: Selection {
                record_hash: record.to_string().repeat(64),
                state_revision: "e".repeat(64),
            },
            context_revision: "f".repeat(64),
            sources: vec![TopoRootBinding {
                role: "same local role".into(),
                event,
            }],
        };
        GraphProjectionSource::RetainedGraph {
            role: "same local role".into(),
            scopes: vec![TopoScope {
                occurrence: 1,
                binding: Arc::new(binding),
            }],
            retained_index: None,
            molecule: source.clone(),
        }
    };
    let mut projections = Vec::new();
    for input in 0..2 {
        let map = vec![segment(0, 4, input * 4, false)];
        let projection = project(&source, &product, &map);
        projections.push(SourceProjection {
            input,
            reference: ReferenceProjection {
                map,
                native: NativeProjection {
                    projection,
                    destination_feature_indices: vec![],
                    previous_feature_indices: vec![],
                },
                introduced: vec![],
            },
        });
    }
    let native = GraphProjections {
        retained: vec![],
        sources: vec![
            make_source(
                'a',
                Some(EventBinding {
                    document_id: "original document".into(),
                    event_id: "exact event".into(),
                    revision: source.revision(),
                }),
            ),
            make_source('b', None),
        ],
        projections,
    };
    let mut builder = TraceBuilder::new(1, Limits::default()).unwrap();
    builder.bind_plan(&"c".repeat(64)).unwrap();
    let mut observer = native::NativeStep::new(Some(&mut builder), 0);
    let mut destination = None;
    observer.graph_sources(&[], &mut destination, &product, Some(native));
    observer.finish(destination, Some(&product)).unwrap();
    let trace = builder.finish();
    let mut lookup = trace.lookup_session();
    for (position, has_event, record) in [(0, true, 'a'), (4, false, 'b')] {
        let calls = lookup.at(0, position, None).unwrap();
        assert!(calls.complete);
        assert_eq!(calls.origins.len(), 1);
        let origin = &calls.origins[0];
        assert_eq!(
            origin.source.source,
            SourceHandle::TopoRole("same local role".into())
        );
        let OriginBinding::TopoGraph {
            selection, sources, ..
        } = &*origin.source.scope[0].binding
        else {
            panic!("actual graph scope required")
        };
        assert_eq!(selection.record_hash, record.to_string().repeat(64));
        assert_eq!(sources[0].event.is_some(), has_event);
    }
}

#[test]
fn native_reference_actions_bind_exact_events_without_synthetic_plan_sources() {
    use crate::material::lineage::EventBinding;
    let source = molecule("ACGT");
    let event = EventBinding {
        document_id: "native document".into(),
        event_id: "exact source event".into(),
        revision: source.revision(),
    };
    let budget = TraceBudget::new(Limits::default()).unwrap();
    let mut builder =
        TraceBuilder::for_reference_edits(budget, &event, &"a".repeat(64), 1).unwrap();
    let mut wrong = event.clone();
    wrong.event_id = "another event with equal DNA".into();
    assert!(builder.register_reference_source(&wrong, &source).is_err());
    let token = builder.register_reference_source(&event, &source).unwrap();
    let mut product = source.clone();
    product.name = "actual renamed presentation".into();
    let mut observer = native::NativeStep::new(Some(&mut builder), 0);
    let mut output = None;
    observer.unchanged_reference(Some(token), &source, &mut output, &product);
    observer.finish(output, Some(&product)).unwrap();
    let trace = builder.finish();
    assert!(
        matches!(trace.context(), TraceContext::ReferenceEdit { source, actions_revision } if source == &event && actions_revision == &"a".repeat(64))
    );
    assert_eq!(
        trace.snapshot_source(output.unwrap()).unwrap().source,
        SourceHandle::ReferenceProduct { action_index: 0 }
    );
    let mut lookup = trace.lookup_session();
    let call = lookup.at(0, 2, None).unwrap();
    assert!(call.complete);
    assert_eq!(
        call.origins[0].source.source,
        SourceHandle::ReferenceEvent(event)
    );
    assert_eq!(call.origins[0].position, 2);
    assert!(
        lookup
            .at(0, 0, Some(LookupStop::Snapshot(output.unwrap())))
            .is_err()
    );
}
