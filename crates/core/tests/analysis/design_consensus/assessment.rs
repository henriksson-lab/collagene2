//! Independent CDS/role/resource oracles against the one native assessment owner.
use super::literals::{CodingLiteral, coding_literals, expected_protein};
use collagene_core::{
    alignment::{
        AlignmentOptions, ComparisonOptions, design_consensus as native, grid::GlobalStorage,
    },
    construct_qc::{AnnotationAssumption, annotation_revision},
    documents::workspace::Workspace,
    material::lineage::EventBinding,
};
use collagene_model::Topology;

fn event(workspace: &Workspace, index: usize) -> EventBinding {
    let document = &workspace.documents[index];
    EventBinding {
        document_id: document.id.clone(),
        event_id: document.head_event.clone().unwrap(),
        revision: document.head.clone(),
    }
}
fn request(case: &CodingLiteral) -> (Workspace, native::Request) {
    let (design, consensus) = case.molecules();
    let annotation = AnnotationAssumption {
        feature_index: 0,
        source_revision: design.revision(),
        annotation_revision: annotation_revision(&design.features[0]).unwrap(),
        rationale:
            "Declared design CDS coordinates for conditional comparison, not annotation truth"
                .into(),
    };
    let workspace = Workspace::from_documents(
        "Independently supplied consensus",
        vec![design, consensus],
        1,
    )
    .unwrap();
    let request = native::Request {
        schema_version: 1,
        design: event(&workspace, 0),
        consensus: event(&workspace, 1),
        comparison: ComparisonOptions {
            alignment: AlignmentOptions {
                global_storage: GlobalStorage::Full,
                ..Default::default()
            },
            both_query_strands: false,
            reference_rotations: false,
            max_origins: 1,
            max_total_cells: 100_000,
        },
        consensus_axis: native::ConsensusAxis::StoredTopology,
        cds: vec![native::CdsSelection {
            id: "explicit CDS".into(),
            design_cds: annotation,
            genetic_code: 1,
        }],
        limits: Default::default(),
    };
    (workspace, request)
}

#[test]
fn literal_selected_path_coding_handles_substitution_indels_reverse_join_and_origin() {
    for case in coding_literals() {
        let (workspace, request) = request(&case);
        let report = native::assess(&workspace, &request).unwrap();
        assert!(report.evidence.selected_path_conditional, "{}", case.name);
        assert_eq!(report.evidence.cds.len(), 1);
        // The supplied observation had no CDS annotations. Conditional coding
        // evidence must not silently install a gene into that source model.
        assert!(
            workspace.revisions[&request.consensus.revision]
                .features
                .is_empty()
        );
        let coding = &report.evidence.cds[0];
        assert_eq!(
            coding.design_translation.as_ref().unwrap().protein,
            expected_protein(case.design_coding),
            "{} design",
            case.name
        );
        assert_eq!(
            coding.projected_translation.as_ref().unwrap().protein,
            expected_protein(case.consensus_coding),
            "{} supplied path",
            case.name
        );
        let stops: Vec<_> = coding
            .consequences
            .iter()
            .filter(|consequence| consequence.kind == native::ConsequenceKind::DefiniteInternalStop)
            .map(|consequence| consequence.projected_codon_index)
            .collect();
        assert_eq!(
            stops,
            case.new_stop.into_iter().collect::<Vec<_>>(),
            "{} conditional internal stops",
            case.name
        );
        report
            .pairwise
            .validate_sources(
                &workspace.revisions[&request.design.revision],
                &workspace.revisions[&request.consensus.revision],
            )
            .unwrap();
    }
}

#[test]
fn literal_synonymous_codon_keeps_explicit_parent_and_observed_calls() {
    let case = coding_literals().remove(0);
    let (workspace, request) = request(&case);
    let report = native::assess(&workspace, &request).unwrap();
    let consequence = report.evidence.cds[0]
        .consequences
        .iter()
        .find(|consequence| consequence.projected_codon_index == 1)
        .unwrap();
    assert_eq!(consequence.kind, native::ConsequenceKind::Synonymous);
    assert_eq!(consequence.design_codon.as_deref(), Some("GAA"));
    assert_eq!(consequence.projected_codon, "GAG");
    assert_eq!(consequence.projected_amino_acid, 'E');
}

#[test]
fn exact_event_ownership_and_annotation_source_changes_cannot_rebind_by_equal_dna() {
    let case = coding_literals().remove(1);
    let (workspace, request) = request(&case);
    let mut wrong_owner = request.clone();
    wrong_owner.design.document_id = wrong_owner.consensus.document_id.clone();
    assert!(native::assess(&workspace, &wrong_owner).is_err());
    let mut same_event = request.clone();
    same_event.consensus = same_event.design.clone();
    assert!(native::assess(&workspace, &same_event).is_err());
    let mut stale_source = request.clone();
    stale_source.design.revision = "0".repeat(64);
    assert!(native::assess(&workspace, &stale_source).is_err());
    let mut stale_annotation = request.clone();
    stale_annotation.cds[0].design_cds.annotation_revision = "0".repeat(64);
    let report = native::assess(&workspace, &stale_annotation).unwrap();
    assert!(!report.evidence.cds[0].mapping_complete);
    assert!(report.evidence.cds[0].projected_translation.is_none());
}

#[test]
fn certificate_limits_and_mapping_budget_withhold_unrestricted_or_complete_claims() {
    let case = CodingLiteral {
        name: "Explicit narrow band with a better omitted placement",
        design: "ATGCAACCCTAA",
        consensus: "TGCAACCCTAAA",
        location: collagene_model::Location::range(0, 12),
        topology: Topology::Linear,
        design_coding: "ATGCAACCCTAA",
        consensus_coding: "TGCAACCCTAAA",
        new_stop: None,
    };
    let (workspace, mut request) = request(&case);
    request.comparison.alignment.global_storage = GlobalStorage::Adaptive {
        initial_radius: 0,
        maximum_radius: 0,
        max_attempts: 1,
    };
    let report = native::assess(&workspace, &request).unwrap();
    assert!(!report.pairwise.comparison.optimality_proven);
    assert!(report.evidence.selected_path_conditional);
    request.comparison.alignment.global_storage = GlobalStorage::Full;
    request.limits.max_mapping_work = 1;
    let report = native::assess(&workspace, &request).unwrap();
    assert!(!report.evidence.coding_complete);
    assert!(
        !report.evidence.cds[0].mapping_complete || !report.evidence.cds[0].consequences_complete
    );
}

#[test]
fn register_segments_distinguish_coding_insertion_from_intron_and_later_restoration() {
    for (consensus, changed, restored) in [
        ("ATGAAAACCAGGTAACGAATAA", true, false),
        ("ATGCAAGCTCCCTAA", false, false),
        ("ATGCAAGAGCCCTAA", false, false),
        ("ATGAAAACCAGTAACGAATAA", true, true),
    ] {
        let case = coding_literals()
            .into_iter()
            .find(|case| case.consensus == consensus)
            .unwrap();
        let (workspace, request) = request(&case);
        let report = native::assess(&workspace, &request).unwrap();
        let coding = &report.evidence.cds[0];
        assert!(coding.register_complete, "{}", case.name);
        assert_eq!(
            coding
                .register_segments
                .iter()
                .any(|segment| segment.offset_mod3 != 0),
            changed,
            "{}",
            case.name
        );
        if restored {
            assert_eq!(coding.register_segments.last().unwrap().offset_mod3, 0);
        }
        if case.new_stop == Some(4) {
            let distant = coding
                .consequences
                .iter()
                .find(|consequence| consequence.projected_codon_index == 4)
                .unwrap();
            assert_eq!(distant.projected_codon, "TAA");
            assert!(distant.design_codon_index.is_none());
            assert!(distant.design_codon.is_none());
        }
    }
}

#[test]
fn joining_a_linear_consensus_origin_requires_an_explicit_single_turn_assumption() {
    let case = coding_literals()
        .into_iter()
        .find(|case| case.topology == Topology::Circular)
        .unwrap();
    let (design, mut consensus) = case.molecules();
    consensus.topology = Topology::Linear;
    let (unused, mut request) = request(&case);
    drop(unused);
    let workspace =
        Workspace::from_documents("Linearized supplied consensus", vec![design, consensus], 1)
            .unwrap();
    request.design = event(&workspace, 0);
    request.consensus = event(&workspace, 1);
    let stored = native::assess(&workspace, &request).unwrap();
    assert!(!stored.evidence.cds[0].mapping_complete);
    request.consensus_axis = native::ConsensusAxis::DeclaredSingleCircularTurn {
        rationale:"The supplied linear record is explicitly declared to contain one complete circular turn; this is an input assumption".into(),
    };
    let declared = native::assess(&workspace, &request).unwrap();
    assert_eq!(
        declared.pairwise.comparison.query_topology,
        Topology::Linear
    );
    assert_eq!(
        declared.evidence.cds[0]
            .projected_translation
            .as_ref()
            .unwrap()
            .protein,
        "M*P*"
    );
    assert!(declared.evidence.selected_path_conditional);
    assert_eq!(declared.evidence.provenance[1].topology, Topology::Linear);
}

#[test]
fn different_event_roles_accept_identical_content_without_inventing_provider_validation() {
    let case = coding_literals().remove(0);
    let (initial, mut request) = request(&case);
    let design = initial.revisions[&request.design.revision].clone();
    let workspace = Workspace::from_documents(
        "Two independent role declarations",
        vec![design.clone(), design],
        1,
    )
    .unwrap();
    request.design = event(&workspace, 0);
    request.consensus = event(&workspace, 1);
    assert_eq!(request.design.revision, request.consensus.revision);
    assert_ne!(request.design.event_id, request.consensus.event_id);
    let report = native::assess(&workspace, &request).unwrap();
    assert_eq!(report.evidence.provenance.len(), 2);
    assert_eq!(report.evidence.provenance[0].event, request.design);
    assert_eq!(report.evidence.provenance[1].event, request.consensus);
    assert!(
        report
            .evidence
            .provenance
            .iter()
            .all(|source| source.original.state == native::OriginalState::Unavailable)
    );
    assert!(report.evidence.selected_path_conditional);
    let encoded = serde_json::to_value(&report.evidence).unwrap();
    assert!(encoded.get("verified").is_none());
    assert!(encoded.get("sample_identity_confirmed").is_none());
}

#[test]
fn explicit_table6_changes_taa_and_tag_to_glutamine_without_a_definite_stop_claim() {
    let case = coding_literals()
        .into_iter()
        .find(|case| case.consensus == "ATGTARCCCTAA")
        .unwrap();
    let (workspace, mut request) = request(&case);
    request.cds[0].genetic_code = 6;
    let report = native::assess(&workspace, &request).unwrap();
    let coding = &report.evidence.cds[0];
    assert_eq!(
        coding.projected_translation.as_ref().unwrap().protein,
        "MQPQ"
    );
    assert_eq!(coding.code.descriptor.id, 6);
    assert_eq!(coding.code.ordinary_amino_acids.len(), 64);
    // TCAG enumeration places TAA and TAG at indices10 and11.
    assert_eq!(&coding.code.ordinary_amino_acids.as_bytes()[10..12], b"QQ");
    assert_eq!(
        &collagene_core::construct_qc::CodePin::for_id(1)
            .unwrap()
            .ordinary_amino_acids
            .as_bytes()[10..12],
        b"**"
    );
    assert!(
        coding
            .consequences
            .iter()
            .all(|consequence| consequence.kind != native::ConsequenceKind::DefiniteInternalStop)
    );
}

#[test]
fn terminal_stop_loss_is_a_conditional_codon_change_not_an_inferred_expression_effect() {
    let case = CodingLiteral {
        name: "Terminal TAA becomes CAA",
        design: "ATGCAACCCTAA",
        consensus: "ATGCAACCCCAA",
        location: collagene_model::Location::range(0, 12),
        topology: Topology::Linear,
        design_coding: "ATGCAACCCTAA",
        consensus_coding: "ATGCAACCCCAA",
        new_stop: None,
    };
    let (workspace, request) = request(&case);
    let report = native::assess(&workspace, &request).unwrap();
    let coding = &report.evidence.cds[0];
    assert_eq!(
        coding.projected_translation.as_ref().unwrap().protein,
        "MQPQ"
    );
    let lost = coding
        .consequences
        .iter()
        .find(|consequence| consequence.kind == native::ConsequenceKind::StopLost)
        .unwrap();
    assert_eq!(lost.projected_codon_index, 3);
    assert_eq!(lost.design_codon.as_deref(), Some("TAA"));
    assert_eq!(lost.projected_codon, "CAA");
    assert!(report.evidence.selected_path_conditional);
}

#[test]
fn original_file_integrity_is_separate_from_normalized_calls_and_provider_claims() {
    use base64::{Engine, engine::general_purpose::STANDARD};
    for (original, correct_hash, expected) in [
        (
            ">provider-label\nATGTAACCCTAA\n",
            true,
            native::OriginalState::IntegrityCheckedMatchesCalls,
        ),
        (
            ">provider-label\nATGCAACCCTAA\n",
            true,
            native::OriginalState::IntegrityCheckedDifferentCalls,
        ),
        (
            ">provider-label\nATGTAACCCTAA\n",
            false,
            native::OriginalState::InvalidOrUnsupported,
        ),
    ] {
        let case = coding_literals().remove(1);
        let (design, mut consensus) = case.molecules();
        consensus.metadata.insert(
            "source:original_fasta_base64".into(),
            STANDARD.encode(original),
        );
        consensus.metadata.insert(
            "source:original_fasta_blake3".into(),
            if correct_hash {
                blake3::hash(original.as_bytes()).to_hex().to_string()
            } else {
                "0".repeat(64)
            },
        );
        consensus.metadata.insert(
            "source:provider".into(),
            "An unverified provider label".into(),
        );
        let (_, mut request) = request(&case);
        let workspace = Workspace::from_documents(
            "Separate original integrity evidence",
            vec![design, consensus],
            1,
        )
        .unwrap();
        request.design = event(&workspace, 0);
        request.consensus = event(&workspace, 1);
        let report = native::assess(&workspace, &request).unwrap();
        assert_eq!(report.evidence.provenance[1].original.state, expected);
        assert_eq!(
            report.evidence.cds[0]
                .projected_translation
                .as_ref()
                .unwrap()
                .protein,
            "M*P*"
        );
        assert!(report.evidence.selected_path_conditional);
    }
}

#[test]
fn missing_or_unreviewed_cds_basis_withholds_that_row_without_erasing_independent_evidence() {
    let case = coding_literals().remove(1);
    let (workspace, mut request) = request(&case);
    let mut missing = request.cds[0].clone();
    missing.id = "Explicitly missing second feature".into();
    missing.design_cds.feature_index = 99;
    request.cds.push(missing);
    let report = native::assess(&workspace, &request).unwrap();
    assert_eq!(report.evidence.cds.len(), 2);
    assert_eq!(
        report.evidence.cds[0]
            .projected_translation
            .as_ref()
            .unwrap()
            .protein,
        "M*P*"
    );
    assert!(!report.evidence.cds[1].mapping_complete);
    assert!(report.evidence.cds[1].projected_translation.is_none());
    assert!(!report.evidence.coding_complete);

    let (mut design, consensus) = case.molecules();
    design.features[0].needs_review = true;
    let workspace =
        Workspace::from_documents("Unreviewed source CDS", vec![design.clone(), consensus], 1)
            .unwrap();
    request.cds.truncate(1);
    request.design = event(&workspace, 0);
    request.consensus = event(&workspace, 1);
    request.cds[0].design_cds.source_revision = design.revision();
    request.cds[0].design_cds.annotation_revision =
        annotation_revision(&design.features[0]).unwrap();
    let report = native::assess(&workspace, &request).unwrap();
    assert!(!report.evidence.cds[0].mapping_complete);
    assert!(report.evidence.cds[0].projected_translation.is_none());
}

#[test]
fn explicit_codon_start_uses_biological_prefix_before_comparing_observed_triplets() {
    let case = CodingLiteral {
        name: "One source call precedes declared codon_start2",
        design: "AATGCAACCCTAA",
        consensus: "AATGTAACCCTAA",
        location: collagene_model::Location::range(0, 13),
        topology: Topology::Linear,
        design_coding: "ATGCAACCCTAA",
        consensus_coding: "ATGTAACCCTAA",
        new_stop: Some(1),
    };
    let (mut design, consensus) = case.molecules();
    design.features[0].qualifiers = vec![("codon_start".into(), Some("2".into()))];
    let (_, mut request) = request(&case);
    let workspace = Workspace::from_documents(
        "Declared biological codon prefix",
        vec![design.clone(), consensus],
        1,
    )
    .unwrap();
    request.design = event(&workspace, 0);
    request.consensus = event(&workspace, 1);
    request.cds[0].design_cds.source_revision = design.revision();
    request.cds[0].design_cds.annotation_revision =
        annotation_revision(&design.features[0]).unwrap();
    let report = native::assess(&workspace, &request).unwrap();
    assert_eq!(
        report.evidence.cds[0]
            .design_translation
            .as_ref()
            .unwrap()
            .protein,
        "MQP*"
    );
    assert_eq!(
        report.evidence.cds[0]
            .projected_translation
            .as_ref()
            .unwrap()
            .protein,
        "M*P*"
    );
}

#[test]
fn a_definite_observed_stop_stays_absolute_when_the_design_parent_is_ambiguous() {
    let case = CodingLiteral {
        name: "Uncertain parent TAN, definite observed TAA",
        design: "ATGTANCCCTAA",
        consensus: "ATGTAACCCTAA",
        location: collagene_model::Location::range(0, 12),
        topology: Topology::Linear,
        design_coding: "ATGTANCCCTAA",
        consensus_coding: "ATGTAACCCTAA",
        new_stop: Some(1),
    };
    let (workspace, request) = request(&case);
    let report = native::assess(&workspace, &request).unwrap();
    let coding = &report.evidence.cds[0];
    assert_eq!(coding.design_translation.as_ref().unwrap().protein, "MXP*");
    assert_eq!(
        coding.projected_translation.as_ref().unwrap().protein,
        "M*P*"
    );
    let stop = coding
        .consequences
        .iter()
        .find(|consequence| consequence.projected_codon_index == 1)
        .unwrap();
    assert_eq!(stop.kind, native::ConsequenceKind::DefiniteInternalStop);
    assert_eq!(stop.design_amino_acid, Some('X'));
    assert!(report.evidence.selected_path_conditional);
}

#[test]
fn origin_turn_through_an_excluded_intron_still_requires_the_declared_query_axis() {
    let case = CodingLiteral {
        name: "Circular CDS crosses stored origin through skipped GGG",
        design: "CCCTAACGTATGCAAGGG",
        consensus: "CCCTAACGTATGTAAGGG",
        location: collagene_model::Location::Join(vec![
            collagene_model::Location::range(9, 15),
            collagene_model::Location::range(0, 6),
        ]),
        topology: Topology::Circular,
        design_coding: "ATGCAACCCTAA",
        consensus_coding: "ATGTAACCCTAA",
        new_stop: Some(1),
    };
    let (design, mut consensus) = case.molecules();
    consensus.topology = Topology::Linear;
    let (_, mut request) = request(&case);
    let workspace = Workspace::from_documents(
        "Origin-spanning excluded intron",
        vec![design, consensus],
        1,
    )
    .unwrap();
    request.design = event(&workspace, 0);
    request.consensus = event(&workspace, 1);
    let stored = native::assess(&workspace, &request).unwrap();
    assert!(!stored.evidence.cds[0].mapping_complete);
    request.consensus_axis = native::ConsensusAxis::DeclaredSingleCircularTurn {rationale:"Explicitly declared single circular query-coordinate turn; not an observed topology fact".into()};
    let declared = native::assess(&workspace, &request).unwrap();
    assert_eq!(
        declared.evidence.cds[0]
            .projected_translation
            .as_ref()
            .unwrap()
            .protein,
        "M*P*"
    );
    assert_eq!(
        declared.pairwise.comparison.query_topology,
        Topology::Linear
    );
}

#[test]
fn semantic_history_keeps_one_path_and_all_evidence_without_restoring_freshness() {
    use collagene_core::{
        alignment::{export, history},
        documents::workspace,
        material::lineage,
    };
    let case = coding_literals().remove(1);
    let (workspace, request) = request(&case);
    let report = native::assess(&workspace, &request).unwrap();
    let expected_evidence = report.evidence.evidence_revision.clone();
    let record = history::Record::from_design_consensus(report).unwrap();
    assert_eq!(record.schema_version, 3);
    let (saved, hash) = history::capture(&workspace, &record).unwrap();
    let reopened = workspace::load_workspace(&workspace::save_workspace(&saved).unwrap()).unwrap();
    let retained = history::load(&reopened, &hash).unwrap();
    let semantics = retained.design_consensus.as_ref().unwrap();
    assert_eq!(semantics.evidence_revision, expected_evidence);
    semantics.validate_for_pairwise(&retained.report).unwrap();
    assert_ne!(
        semantics.workspace_fingerprint,
        lineage::workspace_fingerprint(&reopened).unwrap()
    );
    assert_eq!(
        semantics.cds[0]
            .projected_translation
            .as_ref()
            .unwrap()
            .protein,
        "M*P*"
    );
    let html = export::export(
        &retained,
        "Design",
        "Supplied consensus",
        false,
        export::Format::Html,
    )
    .unwrap();
    let html = String::from_utf8(html.bytes).unwrap();
    assert!(html.contains(&expected_evidence));
    assert!(html.contains("M*P*"));
    let pdf = export::export(
        &retained,
        "Design",
        "Supplied consensus",
        false,
        export::Format::Pdf,
    )
    .unwrap();
    assert!(pdf.bytes.starts_with(b"%PDF"));
    let exchange = export::export(
        &retained,
        "Design",
        "Supplied consensus",
        false,
        export::Format::AlignedFasta,
    )
    .unwrap();
    assert!(!exchange.warnings.is_empty());
    let mut swapped = retained.clone();
    swapped.design_consensus.as_mut().unwrap().request.design = request.consensus;
    assert!(swapped.validate_documentary().is_err());
}

#[test]
fn independently_resealed_semantic_artifact_cannot_change_the_declared_code_data() {
    let case = coding_literals()
        .into_iter()
        .find(|case| case.consensus == "ATGTARCCCTAA")
        .unwrap();
    let (workspace, mut request) = request(&case);
    request.cds[0].genetic_code = 6;
    let mut report = native::assess(&workspace, &request).unwrap();
    report.evidence.cds[0]
        .code
        .ordinary_amino_acids
        .replace_range(10..11, "*");
    let evidence = &report.evidence;
    // Recompute the envelope's integrity independently so rejection exercises
    // canonical code-data association, not a stale outer checksum alone.
    let body = serde_json::to_vec(&(
        &evidence.engine_version,
        &evidence.context_revision,
        &evidence.workspace_fingerprint,
        &evidence.request_revision,
        &evidence.pairwise_context_revision,
        &evidence.pairwise_report_revision,
        &evidence.request,
        &evidence.provenance,
        &evidence.cds,
        &evidence.work,
        evidence.coding_complete,
        evidence.selected_path_conditional,
        &evidence.statement,
    ))
    .unwrap();
    report.evidence.evidence_revision = blake3::hash(&body).to_hex().to_string();
    assert!(
        report
            .validate_binding(&native::revision(&workspace, &request).unwrap())
            .is_err()
    );
}

#[test]
fn long_multibyte_or_sanitized_provider_fields_are_explicit_previews_not_whole_claims() {
    for provider in [
        "🧬".repeat(1024),
        "Lab\nAlias".into(),
        "Exact short label".into(),
    ] {
        let case = coding_literals().remove(1);
        let (design, mut consensus) = case.molecules();
        consensus
            .metadata
            .insert("source:provider".into(), provider.clone());
        let (_, mut request) = request(&case);
        let workspace = Workspace::from_documents(
            "Bounded Unicode provider display",
            vec![design, consensus],
            1,
        )
        .unwrap();
        request.design = event(&workspace, 0);
        request.consensus = event(&workspace, 1);
        let report = native::assess(&workspace, &request).unwrap();
        let field = report.evidence.provenance[1]
            .fields
            .iter()
            .find(|field| field.key == "source:provider")
            .unwrap();
        assert!(field.value.len() <= 2048);
        assert_eq!(field.value_complete, provider == "Exact short label");
        assert_eq!(
            report.evidence.cds[0]
                .projected_translation
                .as_ref()
                .unwrap()
                .protein,
            "M*P*"
        );
        assert_eq!(
            workspace.revisions[&request.consensus.revision].metadata["source:provider"],
            provider
        );
        report.validate_request(&workspace, &request).unwrap();
    }
}
