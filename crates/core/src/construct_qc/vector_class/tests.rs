//! Synthetic four-base elements exercise declared graph facts, never viral DNA.
use super::*;
use crate::construct_qc::{Provenance, annotation_revision, reference_revision};
use collagene_model::{Feature, Location, Sequence};

fn reference(id: &str, dna: &str) -> ReferenceElement {
    let sequence = Sequence::dna(dna).unwrap();
    ReferenceElement {
        id: id.into(),
        sequence_revision: reference_revision(&sequence).unwrap(),
        sequence,
        provenance: Provenance {
            source_url: "https://example.org/authored-synthetic-test-elements".into(),
            accession: format!("synthetic-{id}-1"),
            retrieved_on: "2026-11-28".into(),
            license: "CC0 authored non-biological fixture".into(),
        },
    }
}
fn rebind(source: &Molecule, request: &mut Request) {
    request.source_revision = source.revision();
    request.annotations = source
        .features
        .iter()
        .enumerate()
        .map(|(feature_index, feature)| AnnotationAssumption {
            feature_index,
            source_revision: source.revision(),
            annotation_revision: annotation_revision(feature).unwrap(),
            rationale:
                "Explicit synthetic fixture role and interval; not biological annotation truth"
                    .into(),
        })
        .collect();
}
fn fixture() -> (Molecule, Request) {
    // GG | AACC | GG | CGTA | TT | ATGC | CC | GGTT | AAAA
    let mut source = Molecule::dna(
        "Synthetic graph, not viral DNA",
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
        let range = Location::range(start, end);
        Feature {
            kind: "misc_feature".into(),
            location: if reverse {
                Location::Complement(Box::new(range))
            } else {
                range
            },
            qualifiers: vec![],
            needs_review: false,
        }
    })
    .collect();
    let mut request = Request {
        source_revision: source.revision(),
        profile: catalogue()
            .unwrap()
            .into_iter()
            .find(|profile| profile.reference.id == profiles::RAAV_ID)
            .unwrap()
            .reference,
        annotations: vec![],
        references: vec![
            reference("repeat", "AACC"),
            reference("regulatory", "CGTA"),
            reference("payload", "ATGC"),
        ],
        bindings: [
            ("left", "five-prime-itr", 0, "repeat"),
            ("promoter", "promoter", 1, "regulatory"),
            ("payload", "payload", 2, "payload"),
            ("right", "three-prime-itr", 3, "repeat"),
        ]
        .into_iter()
        .map(|(id, role_id, feature_index, reference_id)| RoleBinding {
            id: id.into(),
            role_id: role_id.into(),
            feature_index,
            reference_id: reference_id.into(),
        })
        .collect(),
        boundaries: Some(BoundarySelection {
            start_binding: "left".into(),
            end_binding: "right".into(),
            reverse: false,
        }),
        spacing: vec![],
        literal_parts: vec![],
        literal_part_limits: Default::default(),
    };
    rebind(&source, &mut request);
    (source, request)
}

#[test]
fn synthetic_literal_premises_yield_declared_geometry_not_biological_acceptance() {
    let (source, request) = fixture();
    assert_eq!(source.sequence.len(), 28);
    let report = assess(&source, &request).unwrap();
    assert!(
        report
            .bindings
            .iter()
            .all(|node| node.reference_verdict == Verdict::Pass && node.premise == Verdict::Pass)
    );
    let boundary = report.boundaries.as_ref().unwrap();
    assert_eq!(
        (
            boundary.start_bond,
            boundary.end_bond,
            boundary.outer_bases,
            boundary.interior_bases
        ),
        (2, 24, 22, 14)
    );
    assert_eq!(
        boundary.source_spans,
        vec![Span {
            start: 2,
            end: 24,
            reverse: false
        }]
    );
    assert!(
        report
            .findings
            .iter()
            .filter(|finding| matches!(
                finding.rule,
                Rule::Order | Rule::DeclaredMultiplicity | Rule::Orientation
            ))
            .all(|finding| finding.verdict == Verdict::Pass)
    );
    assert!(
        report
            .findings
            .iter()
            .filter(|finding| finding.rule == Rule::EmpiricalContext)
            .all(|finding| finding.verdict == Verdict::Informational)
    );
    report.validate_for(&source, &request).unwrap();
    assert!(
        serde_json::to_value(&report)
            .unwrap()
            .get("status")
            .is_none()
    );
}

#[test]
fn sin_labelled_graph_uses_the_same_engine_with_optional_roles_and_no_functional_certification() {
    let mut source = Molecule::dna(
        "Synthetic conventional graph; not lentiviral DNA",
        &"ACGTTT".repeat(6),
        Topology::Linear,
    )
    .unwrap();
    source.features = (0..6)
        .map(|index| Feature {
            kind: "misc_feature".into(),
            location: Location::range(index * 6, index * 6 + 4),
            qualifiers: vec![],
            needs_review: false,
        })
        .collect();
    let (_, mut request) = fixture();
    request.profile = catalogue()
        .unwrap()
        .into_iter()
        .find(|profile| profile.reference.id == profiles::SIN_LENTI_ID)
        .unwrap()
        .reference;
    request.references = vec![reference("synthetic-unit", "ACGT")];
    request.bindings = [
        "five-prime-ltr",
        "psi",
        "rre",
        "internal-promoter",
        "payload",
        "three-prime-sin-ltr",
    ]
    .into_iter()
    .enumerate()
    .map(|(feature_index, role)| RoleBinding {
        id: format!("node-{feature_index}"),
        role_id: role.into(),
        feature_index,
        reference_id: "synthetic-unit".into(),
    })
    .collect();
    request.boundaries = Some(BoundarySelection {
        start_binding: "node-0".into(),
        end_binding: "node-5".into(),
        reverse: false,
    });
    rebind(&source, &mut request);
    let report = assess(&source, &request).unwrap();
    assert_eq!(report.boundaries.as_ref().unwrap().outer_bases, 34);
    assert_eq!(report.source_length, 36);
    assert!(
        report
            .findings
            .iter()
            .filter(|finding| finding.rule == Rule::Order)
            .all(|finding| finding.verdict == Verdict::Pass)
    );
    assert_eq!(
        report
            .findings
            .iter()
            .filter(|finding| finding.rule == Rule::DeclaredMultiplicity
                && finding.verdict == Verdict::Indeterminate)
            .count(),
        2
    );
    assert!(report.profile.definition.numeric_context.is_empty());
    assert!(
        report
            .bindings
            .iter()
            .all(|node| node.reference_verdict == Verdict::Pass)
    );
    report.validate_for(&source, &request).unwrap();
}

#[test]
fn missing_or_unverified_role_premises_remain_unknown_and_identity_failure_remains_failed() {
    let (source, mut request) = fixture();
    request.bindings.clear();
    request.boundaries = None;
    request.annotations.clear();
    request.references.clear();
    let report = assess(&source, &request).unwrap();
    assert!(report.boundaries.is_none());
    assert!(
        report
            .findings
            .iter()
            .filter(|finding| finding.rule == Rule::DeclaredMultiplicity)
            .all(|finding| finding.verdict == Verdict::Indeterminate)
    );
    let (source, mut request) = fixture();
    request.references[0] = reference("repeat", "AAAA");
    let report = assess(&source, &request).unwrap();
    assert_eq!(report.bindings[0].reference_verdict, Verdict::Fail);
    assert_eq!(report.bindings[0].premise, Verdict::Indeterminate);
    assert!(report.boundaries.is_none());
    assert!(
        report
            .findings
            .iter()
            .any(|finding| finding.rule == Rule::ReferenceIntegrity
                && finding.verdict == Verdict::Fail
                && finding.binding_ids == ["left"])
    );
    request.references[0].sequence_revision = "0".repeat(64);
    let report = assess(&source, &request).unwrap();
    assert_eq!(report.bindings[0].reference_verdict, Verdict::Indeterminate);
}

#[test]
fn source_annotation_and_reference_pins_invalidate_geometry_without_role_inference() {
    let (mut source, mut request) = fixture();
    request.annotations[0].annotation_revision = "0".repeat(64);
    let report = assess(&source, &request).unwrap();
    assert_eq!(report.bindings[0].premise, Verdict::Indeterminate);
    assert!(report.boundaries.is_none());
    rebind(&source, &mut request);
    source.features[0].needs_review = true;
    request.source_revision = source.revision();
    let report = assess(&source, &request).unwrap();
    assert_eq!(report.bindings[0].premise, Verdict::Indeterminate);
    let (source, mut request) = fixture();
    request.source_revision = "0".repeat(64);
    assert!(assess(&source, &request).is_err());
}

#[test]
fn orientation_spacing_and_extra_qualified_nodes_have_explicit_independent_failures() {
    let (mut source, mut request) = fixture();
    source.features[3].location = Location::range(20, 24);
    request.references.push(reference("forward-right", "GGTT"));
    request.bindings[3].reference_id = "forward-right".into();
    request.spacing.push(DeclaredSpacing {
        id: "explicit-gap".into(),
        upstream_binding: "promoter".into(),
        downstream_binding: "payload".into(),
        minimum_gap: 3,
        maximum_gap: Some(9),
        rationale: "User's independent gap constraint, not an empirical threshold".into(),
    });
    rebind(&source, &mut request);
    let report = assess(&source, &request).unwrap();
    assert_eq!(report.boundaries.unwrap().outer_bases, 22);
    assert!(
        report
            .findings
            .iter()
            .any(|finding| finding.rule == Rule::Orientation
                && finding.verdict == Verdict::Fail
                && finding.binding_ids == ["right"])
    );
    assert!(
        report.findings.iter().any(
            |finding| finding.rule == Rule::DeclaredSpacing && finding.verdict == Verdict::Fail
        )
    );
    let (mut source, mut request) = fixture();
    source.sequence = Sequence::dna("GGAACCGGCGTATTATGCCCGGTTATGC").unwrap();
    source.features.push(Feature {
        kind: "misc_feature".into(),
        location: Location::range(24, 28),
        qualifiers: vec![],
        needs_review: false,
    });
    request.bindings.push(RoleBinding {
        id: "another-payload-claim".into(),
        feature_index: 4,
        ..request.bindings[2].clone()
    });
    rebind(&source, &mut request);
    let report = assess(&source, &request).unwrap();
    assert!(
        report
            .findings
            .iter()
            .any(|finding| finding.rule == Rule::DeclaredMultiplicity
                && finding.verdict == Verdict::Fail)
    );
    assert!(report.findings.iter().any(|finding| finding.rule == Rule::Order && finding.verdict == Verdict::Indeterminate));
}

#[test]
fn circular_and_reverse_axes_preserve_original_boundary_bonds_and_single_turn_ranges() {
    let (mut source, mut request) = fixture();
    // Rotate at a gap, so all four component annotations remain exact intervals.
    let dna = source.sequence.as_str();
    source.sequence = Sequence::dna(&format!("{}{}", &dna[7..], &dna[..7])).unwrap();
    source.topology = Topology::Circular;
    source.features[0].location = Location::range(23, 27);
    source.features[1].location = Location::range(1, 5);
    source.features[2].location = Location::range(7, 11);
    source.features[3].location = Location::Complement(Box::new(Location::range(13, 17)));
    rebind(&source, &mut request);
    let report = assess(&source, &request).unwrap();
    let boundary = report.boundaries.unwrap();
    assert_eq!(
        (
            boundary.start_bond,
            boundary.end_bond,
            boundary.outer_bases,
            boundary.interior_bases
        ),
        (23, 17, 22, 14)
    );
    assert_eq!(
        boundary.source_spans,
        vec![
            Span {
                start: 23,
                end: 28,
                reverse: false
            },
            Span {
                start: 0,
                end: 17,
                reverse: false
            }
        ]
    );
    let (mut source, mut request) = fixture();
    source.sequence = Sequence::dna("TTTTAACCGGGCATAATACGCCGGTTCC").unwrap();
    source.features[0].location = Location::Complement(Box::new(Location::range(22, 26)));
    source.features[1].location = Location::Complement(Box::new(Location::range(16, 20)));
    source.features[2].location = Location::Complement(Box::new(Location::range(10, 14)));
    source.features[3].location = Location::range(4, 8);
    request.boundaries.as_mut().unwrap().reverse = true;
    rebind(&source, &mut request);
    let report = assess(&source, &request).unwrap();
    let boundary = report.boundaries.unwrap();
    assert_eq!(
        (
            boundary.start_bond,
            boundary.end_bond,
            boundary.outer_bases,
            boundary.interior_bases
        ),
        (26, 4, 22, 14)
    );
    assert_eq!(
        boundary.source_spans,
        vec![Span {
            start: 4,
            end: 26,
            reverse: true
        }]
    );
}

#[test]
fn coincident_circular_outer_bonds_include_both_boundaries_in_one_full_turn() {
    let (mut source, mut request) = fixture();
    source.topology = Topology::Circular;
    source.features[3].location = Location::Complement(Box::new(Location::range(0, 2)));
    request.references.push(reference("short-right", "CC"));
    request.bindings[3].reference_id = "short-right".into();
    rebind(&source, &mut request);
    let report = assess(&source, &request).unwrap();
    let boundary = report.boundaries.unwrap();
    assert_eq!((boundary.start_bond, boundary.end_bond), (2, 2));
    assert_eq!((boundary.outer_bases, boundary.interior_bases), (28, 22));
    assert_eq!(
        boundary.source_spans,
        vec![
            Span {
                start: 2,
                end: 28,
                reverse: false
            },
            Span {
                start: 0,
                end: 2,
                reverse: false
            },
        ]
    );
}

#[test]
fn supplied_forbidden_literal_uses_current_exhaustive_search_and_partial_is_not_absence() {
    let (source, mut request) = fixture();
    request
        .literal_parts
        .push(super::super::required_parts::RequiredPart {
            id: "forbidden-synthetic-repeat".into(),
            role: "User excluded literal, not biological function".into(),
            reference_id: "repeat".into(),
            minimum_copies: 0,
            maximum_copies: Some(0),
            orientation: super::super::required_parts::RequiredOrientation::Either,
        });
    let report = assess(&source, &request).unwrap();
    assert_eq!(report.literal_parts[0].verdict, Verdict::Fail);
    assert_eq!(report.literal_parts[0].certain_copies, 2);
    assert!(report.literal_searches[0].complete);
    request.literal_part_limits.max_work = 1;
    let report = assess(&source, &request).unwrap();
    assert_eq!(report.literal_parts[0].verdict, Verdict::Indeterminate);
    assert!(!report.literal_searches[0].complete);
    assert!(report.literal_search_work <= 1);
}

#[test]
fn current_profile_report_and_required_nullable_boundaries_reject_tampering() {
    let (source, request) = fixture();
    let report = assess(&source, &request).unwrap();
    let reopened: Report = serde_json::from_slice(&serde_json::to_vec(&report).unwrap()).unwrap();
    reopened.validate_for(&source, &request).unwrap();
    let mut profile = report.profile.clone();
    profile.definition.roles[0].maximum = 5;
    assert!(profile.validate().is_err());
    let mut changed = report.clone();
    changed.boundaries.as_mut().unwrap().outer_bases += 1;
    assert!(changed.validate_binding(&report.context_revision).is_err());
    let mut changed_request = request.clone();
    changed_request.profile.identity = "0".repeat(64);
    assert!(assess(&source, &changed_request).is_err());
    let mut wire = serde_json::to_value(&request).unwrap();
    wire.as_object_mut().unwrap().remove("boundaries");
    assert!(serde_json::from_value::<Request>(wire).is_err());
}

#[test]
fn non_dna_deep_native_annotations_and_oversized_role_sets_fail_borrowed_admission() {
    let (mut source, request) = fixture();
    let mut location = Location::range(0, 1);
    for _ in 0..100 {
        location = Location::Complement(Box::new(location));
    }
    source.features[0].location = location;
    assert!(assess(&source, &request).is_err());
    let (source, mut request) = fixture();
    request.bindings = vec![request.bindings[0].clone(); 65];
    assert!(assess(&source, &request).is_err());
    let (mut protein, mut request) = fixture();
    protein.sequence = Sequence::new(Alphabet::Protein, "MPL").unwrap();
    protein.features.clear();
    protein.topology = Topology::Linear;
    protein.validate().unwrap();
    request.source_revision = protein.revision();
    assert!(
        assess(&protein, &request)
            .unwrap_err()
            .to_string()
            .contains("declared DNA vector required")
    );
}
