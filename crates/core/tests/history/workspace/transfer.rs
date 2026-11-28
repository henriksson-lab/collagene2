use collagene_core::{
    documents::transfer::*,
    planning::plan::{DesignPlan, Operation, PlanInput, Source},
    verification::{engine::verify, evidence::Status, output::export_protocol_html},
};
use collagene_model::{Molecule, Topology};

fn fixture() -> (PlanningWorkspace, Candidate) {
    let left = Molecule::dna("left", "AAAACCCC", Topology::Linear).unwrap();
    let right = Molecule::dna("right", "CCCCGGGG", Topology::Linear).unwrap();
    let target = Molecule::dna("target", "AAAACCCCGGGG", Topology::Linear).unwrap();
    let plan = DesignPlan {
        standard: None,
        name: "Independent target".into(),
        target_revision: target.revision(),
        inputs: vec![
            PlanInput {
                id: "left".into(),
                revision: left.revision(),
                molecule: left.clone(),
            },
            PlanInput {
                id: "right".into(),
                revision: right.revision(),
                molecule: right.clone(),
            },
        ],
        steps: vec![Operation::Overlap {
            fragments: vec![Source::Input("left".into()), Source::Input("right".into())],
            overlaps: vec![4],
            circular: false,
        }],
    };
    (
        PlanningWorkspace {
            target_revision: target.revision(),
            documents: vec![left, right, target],
        },
        Candidate::new(plan),
    )
}

#[test]
fn binary_exchange_preserves_exact_document_revisions_and_origin_annotations() {
    let molecule = collagene_core::formats::sequence::import_text(include_str!(
        "../../../../../fixtures/origin-spanning.gb"
    ))
    .unwrap()
    .molecules
    .remove(0);
    let workspace = PlanningWorkspace {
        target_revision: molecule.revision(),
        documents: vec![molecule.clone()],
    };
    let bytes = save_workspace(&workspace).unwrap();
    assert!(bytes.starts_with(b"CGD\0"));
    let restored = load_workspace(&bytes).unwrap();
    assert_eq!(restored.target_revision, molecule.revision());
    assert_eq!(restored.documents[0], molecule);
    assert_eq!(restored.documents[0].topology, Topology::Circular);
    let mut corrupt = bytes;
    let last = corrupt.len() - 1;
    corrupt[last] ^= 1;
    assert!(load_workspace(&corrupt).is_err());
}

#[test]
fn candidate_reconstructs_only_against_matching_current_documents() {
    let (workspace, candidate) = fixture();
    let bytes = save_candidate(&candidate).unwrap();
    assert!(bytes.starts_with(b"CGP\0"));
    let restored = load_candidate(&bytes).unwrap();
    let target = &workspace.documents[2];
    let plan = restored.bind(target, &workspace.documents).unwrap();
    assert_eq!(
        verify(&plan, target).status,
        Status::ConsistentUnderAssumptions
    );
    let mut edited = workspace.documents.clone();
    edited[0] = edited[0].replace(0, 1, "T").unwrap();
    assert!(
        restored
            .bind(target, &edited)
            .unwrap_err()
            .to_string()
            .contains("missing or has changed")
    );
    let mut renamed = workspace.documents.clone();
    renamed[0].name = "renamed".into();
    assert!(restored.bind(target, &renamed).is_err());
    let changed_target = target.replace(0, 1, "T").unwrap();
    assert!(
        restored
            .bind(&changed_target, &workspace.documents)
            .unwrap_err()
            .to_string()
            .contains("active target differs")
    );
    assert!(restored.bind(target, &workspace.documents[1..]).is_err());
}

#[test]
fn packages_reject_rebound_snapshots_and_do_not_promote_bad_science() {
    let (mut workspace, mut candidate) = fixture();
    workspace.target_revision = "missing".into();
    assert!(save_workspace(&workspace).is_err());
    workspace.target_revision = workspace.documents[2].revision();
    workspace.documents.push(workspace.documents[0].clone());
    assert!(save_workspace(&workspace).is_err());
    candidate.plan.name = "tampered".into();
    assert!(save_candidate(&candidate).is_err());
    candidate = Candidate::new(candidate.plan);
    candidate.plan.inputs[0].revision = "stale".into();
    candidate.plan_revision = collagene_core::planning::plan::DesignPlan::revision(&candidate.plan);
    assert!(save_candidate(&candidate).is_err());
    let (workspace, mut candidate) = fixture();
    candidate.plan.steps[0] = Operation::Overlap {
        fragments: vec![Source::Input("left".into()), Source::Input("right".into())],
        overlaps: vec![3],
        circular: false,
    };
    let bad = Candidate::new(candidate.plan);
    // A candidate package records a proposal, not a scientific verdict.
    let imported = load_candidate(&save_candidate(&bad).unwrap()).unwrap();
    let bound = imported
        .bind(&workspace.documents[2], &workspace.documents)
        .unwrap();
    assert_eq!(
        verify(&bound, &workspace.documents[2]).status,
        Status::Invalid
    );
}

#[test]
fn invalid_step_indices_render_diagnostics_without_integer_overflow() {
    let (workspace, mut candidate) = fixture();
    candidate.plan.steps = vec![Operation::Overlap {
        fragments: vec![Source::Step(usize::MAX)],
        overlaps: vec![],
        circular: false,
    }];
    let target = &workspace.documents[2];
    assert_eq!(verify(&candidate.plan, target).status, Status::Invalid);
    let html = export_protocol_html(&candidate.plan, target);
    assert!(html.contains("<strong>Invalid</strong>"));
    assert!(html.contains("invalid step index"));
}

fn claimed_count(payload: &mut Vec<u8>, count: usize) {
    payload.extend(postcard::to_allocvec(&count).unwrap());
}
// Postcard only exposes a sequence's known length when at least that many
// encoded bytes remain. Invalid varint bytes satisfy this admission-hint condition
// while ensuring that requesting the first element would fail during decoding.
fn oversized_collection(payload: &mut Vec<u8>, count: usize) {
    claimed_count(payload, count);
    payload.resize(payload.len() + count, 0xff);
}
fn current_frame_with_payload(template: &[u8], payload: &[u8]) -> Vec<u8> {
    let mut frame = template[..42].to_vec();
    frame[6..10].copy_from_slice(&u32::try_from(payload.len()).unwrap().to_le_bytes());
    frame[10..42].copy_from_slice(blake3::hash(payload).as_bytes());
    frame.extend_from_slice(payload);
    frame
}

#[test]
fn native_candidate_and_planning_counts_reject_before_decoding_any_claimed_element() {
    let (workspace, candidate) = fixture();
    let candidate_frame = save_candidate(&candidate).unwrap();
    let workspace_frame = save_workspace(&workspace).unwrap();

    // Supply enough invalid encoded bytes for Postcard to expose its count hint.
    // A custom admission error, rather than an invalid-varint error, proves the
    // count is rejected before an Operation, PlanInput or Molecule is requested.
    for count in [1001, 1_000_000] {
        // Candidate revision, required nullable declaration, then plan name/target.
        let prefix = postcard::to_allocvec(&(
            "a".repeat(64),
            None::<collagene_core::cloning::assembly::standards::declaration::Declaration>,
            "bounded",
            "b".repeat(64),
        ))
        .unwrap();
        for steps in [false, true] {
            let mut payload = prefix.clone();
            if steps {
                claimed_count(&mut payload, 0); // Empty input list precedes steps.
            }
            oversized_collection(&mut payload, count);
            assert_eq!(
                postcard::from_bytes::<Candidate>(&payload).unwrap_err(),
                postcard::Error::SerdeDeCustom
            );
            assert!(
                load_candidate(&current_frame_with_payload(&candidate_frame, &payload)).is_err()
            );
        }
        let mut payload = postcard::to_allocvec(&"b".repeat(64)).unwrap();
        oversized_collection(&mut payload, count);
        assert_eq!(
            postcard::from_bytes::<PlanningWorkspace>(&payload).unwrap_err(),
            postcard::Error::SerdeDeCustom
        );
        assert!(load_workspace(&current_frame_with_payload(&workspace_frame, &payload)).is_err());
    }

    let json = serde_json::to_vec(&candidate).unwrap();
    let restored: Candidate = serde_json::from_slice(&json).unwrap();
    assert_eq!(restored.plan_revision, candidate.plan_revision);
    assert_eq!(save_candidate(&restored).unwrap(), candidate_frame);
    assert_eq!(
        save_workspace(&load_workspace(&workspace_frame).unwrap()).unwrap(),
        workspace_frame
    );
}

#[test]
fn current_operation_collection_prefixes_admit_only_existing_semantic_capacities() {
    // Current Operation tags whose first field is a variable Source,
    // OrientedSource or TopoSource list. Invalid varint bytes after the count
    // expose Postcard's known count hint but cannot decode the first element.
    for tag in [1u32, 4, 6, 10, 14, 15, 18, 19, 20, 23, 25] {
        let mut payload = postcard::to_allocvec(&tag).unwrap();
        oversized_collection(&mut payload, 1001);
        assert_eq!(
            postcard::from_bytes::<Operation>(&payload).unwrap_err(),
            postcard::Error::SerdeDeCustom,
            "operation tag {tag} must bound its dependency list before allocation"
        );
    }
    for tag in [1u32, 10] {
        let mut payload = postcard::to_allocvec(&tag).unwrap();
        claimed_count(&mut payload, 0); // Empty fragment list.
        oversized_collection(&mut payload, 1001); // Overlap indices.
        assert_eq!(
            postcard::from_bytes::<Operation>(&payload).unwrap_err(),
            postcard::Error::SerdeDeCustom
        );
    }
    let mut digest = postcard::to_allocvec(&(2u32, Source::Step(0))).unwrap();
    oversized_collection(&mut digest, 129);
    assert_eq!(
        postcard::from_bytes::<Operation>(&digest).unwrap_err(),
        postcard::Error::SerdeDeCustom
    );

    // Catalog-pinned digest operations also decode their selected-enzyme list
    // before scientific validation; preserve its existing 64-enzyme limit.
    use collagene_core::cloning::restriction::catalog::{CatalogSelection, SCHEMA_VERSION};
    let mut selection =
        postcard::to_allocvec(&(SCHEMA_VERSION, "catalogue", "a".repeat(64))).unwrap();
    oversized_collection(&mut selection, 65);
    assert_eq!(
        postcard::from_bytes::<CatalogSelection>(&selection).unwrap_err(),
        postcard::Error::SerdeDeCustom
    );
    let mut nested = postcard::to_allocvec(&(17u32, Source::Step(0))).unwrap();
    nested.extend_from_slice(&selection);
    assert_eq!(
        postcard::from_bytes::<Operation>(&nested).unwrap_err(),
        postcard::Error::SerdeDeCustom
    );

    // Decoder guards retain the existing boundary; no stricter cloning limit is
    // introduced. This is shape admission, not a scientifically valid assembly.
    let operation = Operation::Overlap {
        fragments: vec![Source::Step(0); 1000],
        overlaps: vec![0; 1000],
        circular: true,
    };
    let bytes = postcard::to_allocvec(&operation).unwrap();
    let restored: Operation = postcard::from_bytes(&bytes).unwrap();
    assert_eq!(postcard::to_allocvec(&restored).unwrap(), bytes);
}

#[test]
fn unknown_json_counts_stop_at_the_bound_before_parsing_an_extra_element() {
    fn truncated_extra(prefix: &str, item: &str, count: usize) -> String {
        let mut json = prefix.to_owned();
        for index in 0..count {
            if index != 0 {
                json.push(',');
            }
            json.push_str(item);
        }
        // Without the collection guard the extra struct/enum would be requested
        // and reach EOF. The bound must instead produce a data-admission error.
        json.push_str(",{ ");
        json
    }
    let operation = r#"{"Overlap":{"fragments":[],"overlaps":[],"circular":false}}"#;
    let json = truncated_extra(
        r#"{"name":"bounded","target_revision":"","inputs":[],"steps":["#,
        operation,
        1000,
    );
    let error = serde_json::from_str::<DesignPlan>(&json).unwrap_err();
    assert!(
        error.is_data() && error.to_string().contains("resource budget"),
        "bound must precede the malformed extra step: {error}"
    );

    let molecule = Molecule::dna("small", "A", Topology::Linear).unwrap();
    let json = truncated_extra(
        r#"{"target_revision":"","documents":["#,
        &serde_json::to_string(&molecule).unwrap(),
        1000,
    );
    let error = serde_json::from_str::<PlanningWorkspace>(&json).unwrap_err();
    assert!(
        error.is_data() && error.to_string().contains("resource budget"),
        "bound must precede the malformed extra document: {error}"
    );

    for (variant, item) in [
        ("Overlap", r#"{"Step":0}"#),
        ("Ligate", r#"{"source":{"Step":0},"reverse":false}"#),
        ("Topo", r#"{"role":"r","source":{"Step":0}}"#),
    ] {
        let field = if variant == "Topo" {
            "sources"
        } else {
            "fragments"
        };
        let json = truncated_extra(&format!("{{\"{variant}\":{{\"{field}\":["), item, 1000);
        let error = serde_json::from_str::<Operation>(&json).unwrap_err();
        assert!(
            error.is_data() && error.to_string().contains("resource budget"),
            "{variant} must stop before its extra element: {error}"
        );
    }
}

#[test]
fn native_savers_reject_models_the_current_decoder_cannot_reopen() {
    use collagene_core::protocol;
    use collagene_model::{Feature, Location};
    for case in 0..4 {
        let mut molecule = Molecule::dna("bounded", "ACGT", Topology::Linear).unwrap();
        let mut feature = Feature {
            kind: "gene".into(),
            location: Location::range(0, 4),
            qualifiers: Vec::new(),
            needs_review: false,
        };
        match case {
            0 => feature.qualifiers = vec![(String::new(), None); 16_385],
            1 => feature.location = Location::Join(vec![Location::range(0, 4); 16_385]),
            2 => {
                molecule
                    .metadata
                    .insert("ordinary-text".into(), "x".repeat(1024 * 1024 + 1));
            }
            3 => molecule.name = "x".repeat(1025),
            _ => unreachable!(),
        }
        molecule.features.push(feature);
        // These remain molecularly valid: failure must be native resource
        // admission before saving, not an unrelated scientific validity error.
        molecule.validate().unwrap();
        assert!(
            protocol::save_document(&molecule).is_err(),
            "CG2 case {case}"
        );
        let workspace = PlanningWorkspace {
            target_revision: molecule.revision(),
            documents: vec![molecule.clone()],
        };
        assert!(save_workspace(&workspace).is_err(), "CGD case {case}");
        let (_, candidate) = fixture();
        let mut plan = candidate.plan;
        plan.inputs[0].revision = molecule.revision();
        plan.inputs[0].molecule = molecule;
        assert!(
            save_candidate(&Candidate::new(plan)).is_err(),
            "CGP case {case}"
        );
    }
}

#[test]
fn native_savers_preserve_import_scale_features_retained_sources_and_workspace_size() {
    use collagene_core::protocol;
    use collagene_model::{Feature, Location};
    let mut molecule = Molecule::dna("import-scale", "ACGT", Topology::Linear).unwrap();
    molecule.features = vec![
        Feature {
            kind: "gene".into(),
            location: Location::range(0, 4),
            qualifiers: Vec::new(),
            needs_review: false,
        };
        100_000
    ];
    molecule
        .metadata
        .insert("trace_original_base64".into(), "x".repeat(1024 * 1024 + 1));
    let bytes = protocol::save_document(&molecule).unwrap();
    assert_eq!(protocol::load_document(&bytes).unwrap(), molecule);

    let documents = (0..1000)
        .map(|index| {
            let mut molecule =
                Molecule::dna(format!("reference-{index}"), "ACGT", Topology::Linear).unwrap();
            molecule.metadata = (0..5)
                .map(|key| (format!("key-{key}"), "value".into()))
                .collect();
            molecule
        })
        .collect::<Vec<_>>();
    let workspace = PlanningWorkspace {
        target_revision: documents[0].revision(),
        documents,
    };
    let restored = load_workspace(&save_workspace(&workspace).unwrap()).unwrap();
    assert_eq!(restored.documents, workspace.documents);

    let (_, candidate) = fixture();
    let mut plan = candidate.plan;
    plan.inputs[0]
        .molecule
        .metadata
        .insert("fastq:original_base64".into(), "x".repeat(1024 * 1024 + 1));
    plan.inputs[0].revision = plan.inputs[0].molecule.revision();
    let candidate = Candidate::new(plan);
    let restored = load_candidate(&save_candidate(&candidate).unwrap()).unwrap();
    assert_eq!(restored.plan_revision, candidate.plan_revision);
    assert_eq!(
        restored.plan.inputs[0].molecule,
        candidate.plan.inputs[0].molecule
    );
}
