//! Literal four-parent reference products and independent annotation/identity checks.
use collagene_core::{
    Feature, Location, Molecule, Sequence, Topology, chemistry,
    multisite_recombination::{self as adapter, AnnotationLimits, Request, Selection, multisite},
    recombination::{
        self as old,
        reference::{Family, Input, Limits, Reaction, SiteVariant, SourceRecord},
    },
};
fn core(id: u8) -> &'static str {
    match id {
        1 => "GTACAAA",
        2 => "GTACAAG",
        3 => "GTACAAC",
        4 => "GTACAAT",
        _ => panic!("fixture specificity"),
    }
}
fn site(family: Family, id: u8) -> String {
    match family {
        Family::L => format!("ACGG{}TCC", core(id)),
        Family::R => format!("TGGA{}GCT", core(id)),
        _ => panic!("fixture family"),
    }
}
fn fixture() -> (Vec<Input>, Request) {
    let mut profile=multisite::Profile{id:"synthetic-four-parent-LR".into(),edition:"independent-literal-v1".into(),sources:vec![],variants:vec![],limitations:vec!["Synthetic reference topology fixture, not measured biological specificity or a manufacturer profile".into()]};
    for id in 1..=4 {
        for (family, prefix) in [(Family::L, "L"), (Family::R, "R")] {
            let key = format!("synthetic-{prefix}{id}");
            let sequence = site(family, id);
            profile.sources.push(SourceRecord {
                id: key.clone(),
                url: "https://example.org/synthetic-reference-fixture".into(),
                raw_sha256: "1".repeat(64),
                vector_name: "Explicit synthetic bookkeeping fixture".into(),
                vector_length: sequence.len(),
                start_1based: 1,
                end_1based: sequence.len(),
                reverse_extraction: false,
                map_url: "https://example.org/synthetic-reference-fixture-map".into(),
                map_sha256: "2".repeat(64),
            });
            profile.variants.push(SiteVariant {
                id: key.clone(),
                family,
                subtype: id,
                sequence,
                overlap_start: 4,
                source_id: key,
            });
        }
    }
    let dna = [
        format!("{}CCCCA{}AAAAC", site(Family::L, 4), site(Family::R, 1)),
        format!("{}GGGGT{}TTTTC", site(Family::L, 1), site(Family::L, 2)),
        format!("{}AAAAG{}CCCCT", site(Family::R, 2), site(Family::L, 3)),
        format!("{}TTTTA{}GGGGA", site(Family::R, 4), site(Family::R, 3)),
    ];
    (
        dna.iter()
            .enumerate()
            .map(|(i, dna)| Input {
                role: format!("original-{i}"),
                molecule: Molecule::dna(
                    format!("Original circular parent{i}"),
                    dna,
                    Topology::Circular,
                )
                .unwrap(),
            })
            .collect(),
        Request {
            reaction: multisite::Request {
                profile,
                reaction: Reaction::LR,
                specificities: vec![4, 1, 2, 3],
                limits: Limits::default(),
            },
            annotation_limits: AnnotationLimits::default(),
        },
    )
}
// Four directed source intervals, before any algorithm output or target is read.
const MAIN: &str = concat!(
    "GTACAAT", "TCC", "CCCCA", "TGGA", "GTACAAA", "TCC", "GGGGT", "ACGG", "GTACAAG", "GCT",
    "AAAAG", "ACGG", "GTACAAC", "GCT", "GGGGA", "TGGA"
);
const BACKBONES: &str = concat!(
    "GTACAAA", "GCT", "AAAAC", "ACGG", "GTACAAT", "GCT", "TTTTA", "TGGA", "GTACAAC", "TCC",
    "CCCCT", "TGGA", "GTACAAG", "TCC", "TTTTC", "ACGG"
);
fn rc(dna: &str) -> String {
    Sequence::dna(dna)
        .unwrap()
        .reverse_complement()
        .unwrap()
        .as_str()
        .into()
}
fn same_circle(a: &str, b: &str) -> bool {
    a.len() == b.len()
        && (format!("{a}{a}").contains(b) || format!("{}{}", rc(a), rc(a)).contains(b))
}
fn feature(label: &str, location: Location) -> Feature {
    Feature {
        kind: "misc_feature".into(),
        location,
        qualifiers: vec![("label".into(), Some(label.into()))],
        needs_review: false,
    }
}
fn label(feature: &Feature) -> &str {
    feature
        .qualifiers
        .iter()
        .find(|(key, _)| key == "label")
        .and_then(|(_, v)| v.as_deref())
        .unwrap()
}
#[test]
fn four_original_parents_produce_literal_main_and_reciprocal_circles_without_chemical_authority() {
    let (inputs, request) = fixture();
    let before = serde_json::to_vec(&inputs).unwrap();
    let batch = adapter::simulate(&inputs, &request).unwrap();
    batch
        .validate_binding(&adapter::revision(&inputs, &request).unwrap())
        .unwrap();
    assert!(batch.reference.complete);
    assert_eq!(batch.reference.candidates.len(), 1);
    assert_eq!(batch.sources.len(), 4);
    assert_eq!(batch.documents.len(), 2);
    assert_eq!(MAIN.len(), 76);
    assert_eq!(BACKBONES.len(), 76);
    assert!(
        batch
            .documents
            .iter()
            .any(|doc| same_circle(doc.molecule.sequence.as_str(), MAIN))
    );
    assert!(
        batch
            .documents
            .iter()
            .any(|doc| same_circle(doc.molecule.sequence.as_str(), BACKBONES))
    );
    assert_eq!(
        batch
            .documents
            .iter()
            .map(|doc| doc.molecule.sequence.len())
            .sum::<usize>(),
        152
    );
    for doc in &batch.documents {
        assert_eq!(doc.reciprocal_components.len(), 2);
        assert_eq!(doc.annotation_evidence.len(), 4);
        assert_eq!(doc.molecule.metadata[old::REFERENCE_ONLY_KEY], "true");
        assert_eq!(doc.molecule.metadata[old::CHEMISTRY_UNMAPPED_KEY], "true");
        assert!(doc.molecule.metadata.contains_key(adapter::RECIPE_KEY));
        assert!(!doc.molecule.metadata.contains_key("tracked_modifications"));
        assert!(chemistry::Duplex::from_molecule(&doc.molecule).is_err());
        doc.validate_binding(&batch.batch_revision).unwrap();
        let selected = batch.selected(&doc.selection).unwrap();
        assert_eq!(selected.molecule, doc.molecule);
    }
    for (i, input) in inputs.iter().enumerate() {
        let mut counts = vec![0; input.molecule.sequence.len()];
        for map in batch.reference.candidates[0]
            .products
            .iter()
            .flat_map(|p| &p.segments)
            .filter(|m| m.input == i)
        {
            for count in &mut counts[map.source_start..map.source_end] {
                *count += 1;
            }
        }
        assert!(counts.iter().all(|count| *count == 1));
    }
    assert_eq!(before, serde_json::to_vec(&inputs).unwrap());
}
#[test]
fn payload_origin_and_cosourced_annotations_keep_source_sequence_and_review_boundaries() {
    let (mut inputs, request) = fixture();
    inputs[1]
        .molecule
        .features
        .push(feature("payload", Location::range(14, 19)));
    inputs[0].molecule.features.push(feature(
        "origin",
        Location::Join(vec![Location::range(33, 38), Location::range(0, 4)]),
    ));
    let mut shared = feature("co-sourced-crossover", Location::range(23, 30));
    shared.kind = "CDS".into();
    shared.qualifiers.extend([
        ("translation".into(), Some("X".into())),
        ("codon_start".into(), Some("1".into())),
    ]);
    inputs[0].molecule.features.push(shared);
    let batch = adapter::simulate(&inputs, &request).unwrap();
    for (name, input) in [("payload", 1), ("origin", 0)] {
        let source = inputs[input]
            .molecule
            .features
            .iter()
            .find(|f| label(f) == name)
            .unwrap();
        let occurrences = batch
            .documents
            .iter()
            .flat_map(|doc| {
                doc.molecule
                    .features
                    .iter()
                    .filter(move |f| label(f) == name)
                    .map(move |f| (doc, f))
            })
            .collect::<Vec<_>>();
        assert!(!occurrences.is_empty());
        for (doc, mapped) in occurrences {
            assert_eq!(
                mapped.location.extract(&doc.molecule.sequence).unwrap(),
                source
                    .location
                    .extract(&inputs[input].molecule.sequence)
                    .unwrap()
            );
        }
    }
    let main = batch
        .documents
        .iter()
        .find(|doc| same_circle(doc.molecule.sequence.as_str(), MAIN))
        .unwrap();
    let shared = main
        .molecule
        .features
        .iter()
        .find(|f| label(f) == "co-sourced-crossover")
        .unwrap();
    assert!(shared.needs_review);
    assert!(
        !shared
            .qualifiers
            .iter()
            .any(|(key, _)| key == "translation" || key == "codon_start")
    );
    assert_eq!(
        shared
            .location
            .extract(&main.molecule.sequence)
            .unwrap()
            .as_str(),
        "GTACAAA"
    );
    assert!(
        main.annotation_evidence[0]
            .co_sourced_segments
            .iter()
            .any(|map| map.source_start <= 23 && map.source_end >= 30)
    );
}
#[test]
fn reversed_and_origin_shifted_parents_retain_annotation_correspondence() {
    let (mut inputs, request) = fixture();
    inputs[1]
        .molecule
        .features
        .push(feature("payload", Location::range(14, 19)));
    for input in &mut inputs {
        input.molecule = collagene_core::annotations::reverse_complement(&input.molecule).unwrap();
        input.molecule = collagene_core::annotations::rotate_origin(&input.molecule, 7).unwrap();
    }
    let batch = adapter::simulate(&inputs, &request).unwrap();
    assert_eq!(batch.documents.len(), 2);
    assert!(
        batch
            .documents
            .iter()
            .any(|doc| same_circle(doc.molecule.sequence.as_str(), MAIN))
    );
    let source = inputs[1]
        .molecule
        .features
        .iter()
        .find(|f| label(f) == "payload")
        .unwrap();
    let occurrences = batch
        .documents
        .iter()
        .flat_map(|doc| {
            doc.molecule
                .features
                .iter()
                .filter(|f| label(f) == "payload")
                .map(move |f| (doc, f))
        })
        .collect::<Vec<_>>();
    assert!(!occurrences.is_empty());
    for (doc, f) in occurrences {
        assert_eq!(
            f.location.extract(&doc.molecule.sequence).unwrap(),
            source
                .location
                .extract(&inputs[1].molecule.sequence)
                .unwrap()
        );
    }
}
#[test]
fn source_selection_and_independent_target_contexts_cannot_be_reused_or_rewritten() {
    let (inputs, request) = fixture();
    let target = Molecule::dna("Frozen independent target", MAIN, Topology::Circular).unwrap();
    let batch = adapter::simulate(&inputs, &request).unwrap();
    let main = batch
        .documents
        .iter()
        .find(|doc| same_circle(doc.molecule.sequence.as_str(), MAIN))
        .unwrap();
    let selected = adapter::product(&inputs, &request, &main.selection).unwrap();
    assert_eq!(selected.molecule, main.molecule);
    let checked = adapter::check_target(&inputs, &request, &main.selection, &target).unwrap();
    assert!(checked.assessment.target_matches);
    assert!(checked.assessment.complete);
    assert_eq!(
        checked.context_revision().unwrap(),
        adapter::check_revision(&inputs, &request, &main.selection, &target).unwrap()
    );
    let wrong = Molecule::dna(
        "Wrong independent target",
        &"A".repeat(76),
        Topology::Circular,
    )
    .unwrap();
    let rejected = adapter::check_target(&inputs, &request, &main.selection, &wrong).unwrap();
    assert!(!rejected.assessment.target_matches);
    assert_eq!(rejected.product.molecule, checked.product.molecule);
    assert_ne!(
        checked.context_revision().unwrap(),
        rejected.context_revision().unwrap()
    );
    let mut foreign = main.selection.clone();
    foreign.batch_revision = "0".repeat(64);
    assert!(batch.selected(&foreign).is_err());
    assert!(adapter::product(&inputs, &request, &foreign).is_err());
    let mut absent = main.selection.clone();
    absent.component_id = 31;
    assert!(batch.selected(&absent).is_err());
    let mut changed = inputs.clone();
    changed[0].role = "same-DNA-other-role".into();
    assert!(adapter::product(&changed, &request, &main.selection).is_err());
    let mut altered = batch.clone();
    altered.documents[0].molecule.sequence = Sequence::dna(&"C".repeat(76)).unwrap();
    assert!(altered.validate_binding(&batch.batch_revision).is_err());
    let mut altered = batch.clone();
    altered.documents[0]
        .molecule
        .metadata
        .remove(old::REFERENCE_ONLY_KEY);
    assert!(altered.validate_binding(&batch.batch_revision).is_err());
}
#[test]
fn aggregate_annotation_and_native_search_budgets_reject_before_copying_partial_documents() {
    let (mut inputs, mut request) = fixture();
    inputs[1]
        .molecule
        .features
        .push(feature("payload", Location::range(14, 19)));
    request.annotation_limits.max_features = 1;
    assert!(adapter::simulate(&inputs, &request).is_err());
    request.annotation_limits = AnnotationLimits::default();
    request.annotation_limits.max_copy_bytes = 1;
    assert!(adapter::simulate(&inputs, &request).is_err());
    request.annotation_limits = AnnotationLimits::default();
    request.annotation_limits.max_report_bytes = 1024;
    assert!(adapter::simulate(&inputs, &request).is_err());
    request.annotation_limits = AnnotationLimits::default();
    request.reaction.limits.max_comparisons = 1;
    let partial = adapter::simulate(&inputs, &request).unwrap();
    assert!(!partial.reference.complete);
    assert!(partial.documents.is_empty());
    let fake = Selection {
        batch_revision: partial.batch_revision.clone(),
        candidate_id: 0,
        component_id: 0,
    };
    assert!(partial.selected(&fake).is_err());
    request.reaction.limits.max_comparisons = 20_000_000;
    request.reaction.specificities.push(1);
    assert!(adapter::request_preflight(&request).is_err());
}

#[test]
fn typed_plan_preserves_whole_batch_and_withholds_physical_ancestry() {
    use collagene_core::verification::{
        self, DesignPlan, Operation, PlanInput, Source, Status, TopoSource,
    };
    let (inputs, request) = fixture();
    let target = Molecule::dna(
        "Independent literal intended target",
        MAIN,
        Topology::Circular,
    )
    .unwrap();
    let batch = adapter::simulate(&inputs, &request).unwrap();
    let selection = batch
        .documents
        .iter()
        .find(|d| same_circle(d.molecule.sequence.as_str(), MAIN))
        .unwrap()
        .selection
        .clone();
    let mut plan = DesignPlan {
        name: "Four original parents, not repeated two-parent exchanges".into(),
        target_revision: target.revision(),
        inputs: inputs
            .iter()
            .enumerate()
            .map(|(i, s)| PlanInput {
                id: format!("source-{i}"),
                revision: s.molecule.revision(),
                molecule: s.molecule.clone(),
            })
            .collect(),
        steps: vec![Operation::MultiSiteRecombineReference {
            sources: inputs
                .iter()
                .enumerate()
                .map(|(i, s)| TopoSource {
                    role: s.role.clone(),
                    source: Source::Input(format!("source-{i}")),
                })
                .collect(),
            request,
            selection,
        }],
    };
    let report = verification::verify(&plan, &target);
    assert_eq!(report.status, Status::Incomplete);
    assert_eq!(report.products.len(), 1);
    assert_eq!(report.reference_step_reports.len(), 1);
    assert!(report.physical_processing_reports.is_empty());
    let verification::ReferenceStepEvidence::MultiSiteRecombination { batch, .. } =
        &report.reference_step_reports[0].evidence;
    assert_eq!(batch.documents.len(), 2);
    assert_eq!(
        batch
            .documents
            .iter()
            .map(|d| d.molecule.sequence.len())
            .sum::<usize>(),
        152
    );
    let html = verification::export_protocol_html(&plan, &target);
    assert!(html.contains("multi-site"));
    let pdf = collagene_core::pdf::export_protocol_pdf(&plan, &target).unwrap();
    assert!(pdf.bytes.starts_with(b"%PDF-"));
    let wrong = Molecule::dna("Independent wrong target", "ACGT", Topology::Circular).unwrap();
    plan.target_revision = wrong.revision();
    let mismatch = verification::verify(&plan, &wrong);
    assert_eq!(mismatch.status, Status::Invalid);
    assert_eq!(mismatch.products.len(), 1);
    assert_eq!(mismatch.reference_step_reports.len(), 1);
    plan.target_revision = target.revision();
    plan.steps.push(Operation::TreatEnds {
        template: Source::Step(0),
        treatment: chemistry::EndTreatment::Phosphorylate,
    });
    let downstream = verification::verify(&plan, &target);
    assert_eq!(downstream.status, Status::Incomplete);
    assert_eq!(
        downstream.products.len(),
        1,
        "reference letters must not manufacture a physical duplex"
    );
}

#[test]
fn retained_multisite_run_requires_exact_immutable_reference_evidence() {
    use collagene_core::{
        run_lineage::{self, EventBinding, InputEventBinding},
        verification::{DesignPlan, Operation, PlanInput, Source, TopoSource},
        workspace::Workspace,
    };
    let (inputs, request) = fixture();
    let target = Molecule::dna(
        "Separate independent intended target",
        MAIN,
        Topology::Circular,
    )
    .unwrap();
    let batch = adapter::simulate(&inputs, &request).unwrap();
    let selection = batch
        .documents
        .iter()
        .find(|d| same_circle(d.molecule.sequence.as_str(), MAIN))
        .unwrap()
        .selection
        .clone();
    let plan = DesignPlan {
        name: "Retain whole reference network".into(),
        target_revision: target.revision(),
        inputs: inputs
            .iter()
            .enumerate()
            .map(|(i, s)| PlanInput {
                id: format!("p{i}"),
                revision: s.molecule.revision(),
                molecule: s.molecule.clone(),
            })
            .collect(),
        steps: vec![Operation::MultiSiteRecombineReference {
            sources: inputs
                .iter()
                .enumerate()
                .map(|(i, s)| TopoSource {
                    role: s.role.clone(),
                    source: Source::Input(format!("p{i}")),
                })
                .collect(),
            request,
            selection,
        }],
    };
    let mut molecules = inputs
        .iter()
        .map(|s| s.molecule.clone())
        .collect::<Vec<_>>();
    molecules.push(target);
    let workspace =
        Workspace::from_documents("Multi-parent immutable history", molecules, 4).unwrap();
    let binding = |index: usize| {
        let d = &workspace.documents[index];
        EventBinding {
            document_id: d.id.clone(),
            event_id: d.head_event.clone().unwrap(),
            revision: d.head.clone(),
        }
    };
    let commit = run_lineage::record_run(
        &workspace,
        &plan,
        binding(4),
        (0..4)
            .map(|i| InputEventBinding {
                input_id: format!("p{i}"),
                event: binding(i),
            })
            .collect(),
        true,
    )
    .unwrap();
    commit.validate().unwrap();
    let history =
        collagene_core::reference_run_evidence::load(&commit.workspace, &commit.receipt.run_id)
            .unwrap()
            .unwrap();
    assert_eq!(history.steps.len(), 1);
    let download =
        collagene_core::run_evidence::download(&commit.workspace, &commit.receipt.run_id).unwrap();
    let bundle: serde_json::Value =
        serde_json::from_slice(download.bytes.as_ref().unwrap()).unwrap();
    assert_eq!(bundle["format"], "scientific-run-bundle-v1");
    assert_eq!(download.step_count, 1);
    let bytes = collagene_core::workspace::save_workspace(&commit.workspace).unwrap();
    let roundtrip = collagene_core::workspace::load_workspace(&bytes).unwrap();
    roundtrip.validate().unwrap();
    let mut missing = roundtrip;
    missing
        .attachments
        .retain(|_, a| a.kind != collagene_core::reference_run_evidence::KIND);
    assert!(missing.validate().is_err());
}

#[test]
fn retained_reference_graphs_bind_exact_source_and_product_events_without_replaying_old_labels() {
    use collagene_core::{
        reference_run_evidence as history,
        run_lineage::{self, EventBinding, InputEventBinding},
        verification::{
            self, DesignPlan, Operation, PlanInput, ReferenceStepEvidence, Source, Status,
            TopoSource,
        },
        workspace::{self, Workspace},
    };
    fn event(w: &Workspace, index: usize) -> EventBinding {
        let d = &w.documents[index];
        EventBinding {
            document_id: d.id.clone(),
            event_id: d.head_event.clone().unwrap(),
            revision: d.head.clone(),
        }
    }
    fn replace(
        w: &mut Workspace,
        hash: &str,
        mut attachment: workspace::WorkspaceAttachment,
        evidence: &history::Evidence,
    ) {
        attachment.bytes = serde_json::to_vec(evidence).unwrap();
        w.attachments.remove(hash);
        w.attachments.insert(
            blake3::hash(&attachment.bytes).to_hex().to_string(),
            attachment,
        );
    }
    let (inputs, request) = fixture();
    let target =
        Molecule::dna("Independent literal main circle", MAIN, Topology::Circular).unwrap();
    let batch = adapter::simulate(&inputs, &request).unwrap();
    // Select by the intended original source interval, not by target comparison.
    let component = batch.reference.candidates[0]
        .products
        .iter()
        .position(|product| {
            product.segments.iter().any(|map| {
                map.input == 0 && map.source_start == 4 && map.source_end == 23 && !map.reverse
            })
        })
        .unwrap();
    let selection = Selection {
        batch_revision: batch.batch_revision,
        candidate_id: 0,
        component_id: component,
    };
    let plan = DesignPlan {
        name: "Four-parent reciprocal reference design".into(),
        target_revision: target.revision(),
        inputs: inputs
            .iter()
            .map(|input| PlanInput {
                id: input.role.clone(),
                revision: input.molecule.revision(),
                molecule: input.molecule.clone(),
            })
            .collect(),
        steps: vec![Operation::MultiSiteRecombineReference {
            sources: inputs
                .iter()
                .map(|input| TopoSource {
                    role: input.role.clone(),
                    source: Source::Input(input.role.clone()),
                })
                .collect(),
            request,
            selection,
        }],
    };
    let report = verification::verify(&plan, &target);
    assert_eq!(report.status, Status::Incomplete);
    assert_eq!(report.reference_step_reports.len(), 1);
    assert!(report.physical_processing_reports.is_empty());
    let mut documents = inputs
        .iter()
        .map(|input| input.molecule.clone())
        .collect::<Vec<_>>();
    documents.push(target.clone());
    let w = Workspace::from_documents("Original four-parent source events", documents, 4).unwrap();
    let bindings = inputs
        .iter()
        .enumerate()
        .map(|(i, input)| InputEventBinding {
            input_id: input.role.clone(),
            event: event(&w, i),
        })
        .collect();
    let commit = run_lineage::record_run(&w, &plan, event(&w, 4), bindings, true).unwrap();
    let archive = workspace::save_workspace(&commit.workspace).unwrap();
    let loaded = workspace::load_workspace(&archive).unwrap();
    let saved = history::load(&loaded, &commit.receipt.run_id)
        .unwrap()
        .unwrap();
    assert_eq!(saved.steps.len(), 1);
    assert_eq!(saved.inputs.len(), 4);
    assert_eq!(saved.products.len(), 1);
    let ReferenceStepEvidence::MultiSiteRecombination { batch, .. } = &saved.steps[0].evidence;
    assert_eq!(batch.documents.len(), 2);
    assert_eq!(
        batch
            .documents
            .iter()
            .map(|d| d.molecule.sequence.len())
            .sum::<usize>(),
        152
    );
    let (hash, attachment) = loaded
        .attachments
        .iter()
        .find(|(_, a)| a.kind == history::KIND)
        .map(|(h, a)| (h.clone(), a.clone()))
        .unwrap();
    let mut historical = loaded.clone();
    let mut older = saved.clone();
    let ReferenceStepEvidence::MultiSiteRecombination { batch, .. } = &mut older.steps[0].evidence;
    batch.engine_version = "historical-adapter-label".into();
    batch.reference.engine = "historical-native-label".into();
    for doc in &mut batch.documents {
        doc.engine_version = "historical-document-label".into();
    }
    replace(&mut historical, &hash, attachment.clone(), &older);
    historical.validate().unwrap();
    assert_eq!(
        history::load(&historical, &commit.receipt.run_id)
            .unwrap()
            .unwrap()
            .steps
            .len(),
        1
    );
    let replay = run_lineage::replay_run(&historical, &commit.receipt.run_id).unwrap();
    assert_eq!(replay.current_report.status, Status::Incomplete);
    assert!(replay.products_match);
    let mut missing = loaded.clone();
    missing.attachments.remove(&hash);
    assert!(missing.validate().is_err());
    let mut tampered = loaded.clone();
    let mut changed = saved;
    let ReferenceStepEvidence::MultiSiteRecombination { batch, .. } =
        &mut changed.steps[0].evidence;
    batch.sources[0].revision = target.revision();
    batch.reference.input_revisions[0] = target.revision();
    replace(&mut tampered, &hash, attachment, &changed);
    assert!(tampered.validate().is_err());
}
