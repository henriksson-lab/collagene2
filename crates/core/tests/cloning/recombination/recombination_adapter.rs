use collagene_core::{
    chemistry,
    cloning::recombination::{self as adapter, reference, *},
    sequence::annotations,
};
use collagene_model::{Feature, Location, Molecule, Sequence, Topology};
fn profile() -> ProfileSnapshot {
    ProfileSnapshot::manufacturer_context().unwrap()
}
fn rc(s: &str) -> String {
    Sequence::dna(s)
        .unwrap()
        .reverse_complement()
        .unwrap()
        .as_str()
        .into()
}
fn site(profile: &ProfileSnapshot, f: reference::Family, n: u8) -> String {
    profile
        .profile
        .variants
        .iter()
        .find(|v| v.family == f && v.subtype == n)
        .unwrap()
        .sequence
        .clone()
}
fn fixture() -> (Vec<reference::Input>, Request) {
    let profile = profile();
    let pcr = format!(
        "{}ATGAAATAA{}",
        site(&profile, reference::Family::B, 1),
        rc(&site(&profile, reference::Family::B, 2))
    );
    let donor = include_str!("../../../../recombination/tests/fixtures/pdonr221_seq.txt")
        .split_whitespace()
        .collect::<String>();
    (
        vec![
            reference::Input {
                role: "PCR".into(),
                molecule: Molecule::dna("PCR", &pcr, Topology::Linear).unwrap(),
            },
            reference::Input {
                role: "donor".into(),
                molecule: Molecule::dna("pDONR221", &donor, Topology::Circular).unwrap(),
            },
        ],
        Request::new(profile, reference::Reaction::BP),
    )
}
fn feature(label: &str, location: Location) -> Feature {
    Feature {
        kind: "CDS".into(),
        location,
        qualifiers: vec![
            ("label".into(), Some(label.into())),
            ("translation".into(), Some("MK".into())),
            ("codon_start".into(), Some("1".into())),
        ],
        needs_review: false,
    }
}
fn label(f: &Feature) -> &str {
    f.qualifiers
        .iter()
        .find(|(k, _)| k == "label")
        .and_then(|(_, v)| v.as_deref())
        .unwrap()
}
#[test]
fn fresh_adapter_preserves_all_components_and_blocks_physical_authority() {
    let (mut inputs, request) = fixture();
    inputs[0].molecule.metadata.insert(
        "tracked_modifications".into(),
        "untrusted-source-summary".into(),
    );
    let report = adapter::simulate(&inputs, &request).unwrap();
    assert_eq!(report.reference.candidates.len(), 1);
    assert_eq!(report.documents.len(), 2);
    assert_eq!(
        report.reference,
        reference::recombine(
            &inputs,
            &request.profile.profile,
            request.reaction,
            &request.limits
        )
        .unwrap()
    );
    report.validate_bindings_for(&inputs, &request).unwrap();
    for doc in &report.documents {
        doc.validate_bindings_for(&inputs, &request, &doc.selection)
            .unwrap();
        assert_eq!(doc.reciprocal_components.len(), 2);
        assert!(!doc.molecule.metadata.contains_key("tracked_modifications"));
        assert!(!doc.molecule.metadata.contains_key("duplex_model"));
        assert_eq!(doc.molecule.metadata[REFERENCE_ONLY_KEY], "true");
        assert!(chemistry::Duplex::from_molecule(&doc.molecule).is_err());
        let shared: collagene_model::Molecule = doc.molecule.clone();
        assert_eq!(shared.sequence, doc.molecule.sequence);
    }
}
#[test]
fn exact_payload_and_origin_annotation_projection_preserves_biological_sequence() {
    let (mut inputs, request) = fixture();
    inputs[0]
        .molecule
        .features
        .push(feature("payload", Location::range(25, 34)));
    inputs[1].molecule.features.push(feature(
        "origin",
        Location::Join(vec![Location::range(4700, 4761), Location::range(0, 30)]),
    ));
    let report = adapter::simulate(&inputs, &request).unwrap();
    let circular = report
        .documents
        .iter()
        .find(|d| d.molecule.topology == Topology::Circular)
        .unwrap();
    for (name, input) in [("payload", 0), ("origin", 1)] {
        let mapped = circular
            .molecule
            .features
            .iter()
            .find(|f| label(f) == name)
            .unwrap();
        let original = inputs[input]
            .molecule
            .features
            .iter()
            .find(|f| label(f) == name)
            .unwrap();
        assert_eq!(
            mapped
                .location
                .extract(&circular.molecule.sequence)
                .unwrap(),
            original
                .location
                .extract(&inputs[input].molecule.sequence)
                .unwrap()
        );
        assert!(!mapped.needs_review);
        assert!(mapped.qualifiers.iter().any(|(k, _)| k == "translation"));
    }
    assert!(
        report
            .documents
            .iter()
            .find(|d| d.molecule.topology == Topology::Linear)
            .unwrap()
            .molecule
            .features
            .is_empty()
    );
}
#[test]
fn reverse_and_origin_views_project_without_coordinate_or_strand_invention() {
    let (mut inputs, request) = fixture();
    inputs[0]
        .molecule
        .features
        .push(feature("payload", Location::range(25, 34)));
    inputs[1].molecule.features.push(feature(
        "origin",
        Location::Join(vec![Location::range(4700, 4761), Location::range(0, 30)]),
    ));
    let expected = inputs
        .iter()
        .flat_map(|i| {
            i.molecule.features.iter().map(|f| {
                (
                    label(f).to_owned(),
                    f.location.extract(&i.molecule.sequence).unwrap(),
                )
            })
        })
        .collect::<Vec<_>>();
    for i in &mut inputs {
        i.molecule = annotations::reverse_complement(&i.molecule).unwrap();
    }
    let report = adapter::simulate(&inputs, &request).unwrap();
    let circle = report
        .documents
        .iter()
        .find(|d| d.molecule.topology == Topology::Circular)
        .unwrap();
    for (name, s) in expected {
        let f = circle
            .molecule
            .features
            .iter()
            .find(|f| label(f) == name)
            .unwrap();
        assert_eq!(f.location.extract(&circle.molecule.sequence).unwrap(), s);
    }
}
#[test]
fn co_sourced_overlap_annotations_retain_both_parent_evidence_with_review() {
    let (mut inputs, request) = fixture();
    inputs[0]
        .molecule
        .features
        .push(feature("B1 overlap", Location::range(8, 15)));
    inputs[1]
        .molecule
        .features
        .push(feature("P1 overlap", Location::range(651, 658)));
    let report = adapter::simulate(&inputs, &request).unwrap();
    for d in &report.documents {
        let b = d
            .molecule
            .features
            .iter()
            .find(|f| label(f) == "B1 overlap")
            .unwrap();
        let p = d
            .molecule
            .features
            .iter()
            .find(|f| label(f) == "P1 overlap")
            .unwrap();
        assert_eq!(b.location.extract(&d.molecule.sequence).unwrap(), "GTACAAA");
        assert_eq!(p.location.extract(&d.molecule.sequence).unwrap(), "GTACAAA");
        assert_eq!(b.location.spans(), p.location.spans());
        assert!(b.needs_review || p.needs_review);
        let evidence = d
            .annotation_evidence
            .iter()
            .find(|e| !e.co_sourced_segments.is_empty())
            .unwrap();
        assert!(evidence.reviewed_features > 0);
        let reviewed = d
            .molecule
            .features
            .iter()
            .filter(|f| f.needs_review)
            .collect::<Vec<_>>();
        assert!(
            reviewed
                .iter()
                .all(|f| !f.qualifiers.iter().any(|(k, _)| k == "translation"))
        );
    }
}
#[test]
fn unsupported_projection_keeps_reviewed_source_features_and_other_parent_annotations() {
    let (mut inputs, request) = fixture();
    let mut unsupported = feature("coordinate qualifier", Location::range(25, 34));
    unsupported
        .qualifiers
        .push(("transl_except".into(), Some("(pos:25..27,aa:Sec)".into())));
    inputs[0].molecule.features.push(unsupported);
    inputs[1]
        .molecule
        .features
        .push(feature("donor backbone", Location::range(3500, 3512)));
    let report = adapter::simulate(&inputs, &request).unwrap();
    let circle = report
        .documents
        .iter()
        .find(|d| d.molecule.topology == Topology::Circular)
        .unwrap();
    assert!(
        circle
            .molecule
            .features
            .iter()
            .any(|f| label(f) == "donor backbone")
    );
    assert!(
        !circle
            .molecule
            .features
            .iter()
            .any(|f| label(f) == "coordinate qualifier")
    );
    assert_eq!(circle.annotation_evidence[0].unmapped_features.len(), 1);
    assert!(circle.annotation_evidence[0].unmapped_features[0].needs_review);
    assert!(
        circle.annotation_evidence[0]
            .warnings
            .iter()
            .any(|w| w.contains("unsupported"))
    );
}
#[test]
fn fresh_checked_selection_matches_exact_component_and_preserves_uncertainty() {
    let (inputs, request) = fixture();
    let report = adapter::simulate(&inputs, &request).unwrap();
    let selected = &report.documents[0];
    let target = &selected.molecule;
    let policy = reference::TargetPolicy {
        allow_reverse_complement: false,
    };
    let checked =
        adapter::checked_product(&inputs, &request, &selected.selection, target, &policy).unwrap();
    checked
        .validate_bindings_for(&inputs, &request, &selected.selection, target, &policy)
        .unwrap();
    let wrong = &report.documents[1].selection;
    assert!(adapter::checked_product(&inputs, &request, wrong, target, &policy).is_err());
    let mut foreign = selected.selection.clone();
    foreign.batch_revision = "0".repeat(64);
    assert!(adapter::product(&inputs, &request, &foreign).is_err());
    let mut forged = report.clone();
    forged.documents[0]
        .molecule
        .features
        .push(feature("forged claim", Location::range(0, 3)));
    assert!(forged.validate_bindings_for(&inputs, &request).is_err());
    let fresh = adapter::product(&inputs, &request, &selected.selection).unwrap();
    assert!(
        !fresh
            .molecule
            .features
            .iter()
            .any(|f| label(f) == "forged claim")
    );
}
#[test]
fn annotation_and_profile_limits_are_atomic_and_borrowed_inputs_are_guarded() {
    let (mut inputs, mut request) = fixture();
    inputs[0]
        .molecule
        .features
        .push(feature("payload", Location::range(25, 34)));
    request.annotation_limits.max_copy_bytes = 1;
    assert!(adapter::simulate(&inputs, &request).is_err());
    request.annotation_limits = AnnotationLimits::default();
    request.profile.profile.edition.push_str("-tampered");
    assert!(
        adapter::source_preflight([&inputs[0].molecule, &inputs[1].molecule], &request).is_err()
    );
    let (_, request) = fixture();
    let mut deep = Location::range(0, 1);
    for _ in 0..64 {
        deep = Location::Complement(Box::new(deep));
    }
    inputs[0]
        .molecule
        .features
        .push(feature("native deep", deep));
    assert!(
        adapter::borrowed_preflight(
            [("PCR", &inputs[0].molecule), ("donor", &inputs[1].molecule)],
            &request,
            None
        )
        .is_err()
    );
}

#[test]
fn native_report_duplicate_document_bases_rejected_before_per_document_validation() {
    let (inputs, request) = fixture();
    let report = adapter::simulate(&inputs, &request).unwrap();
    let mut forged = report.clone();
    let mut giant = report.documents[0].clone();
    giant.molecule.sequence = Sequence::dna(&"A".repeat(2_000_000)).unwrap();
    // A later per-document parse would fail Unsupported on the stale recipe/content
    // pin. The aggregate limit must reject TooLarge first, with no repeated copies,
    // hashes or serializations of64*2M bases.
    forged.documents = vec![giant; 64];
    assert!(matches!(
        forged.validate(),
        Err(collagene_model::Error::TooLarge)
    ));
    let mut forged = report.clone();
    forged.documents[0].annotation_evidence[0].warnings = vec!["x".repeat(65536)];
    assert!(matches!(
        forged.validate(),
        Err(collagene_model::Error::TooLarge)
    ));
}

#[test]
fn native_report_large_repeated_unmapped_text_is_capped_before_each_row_serializes() {
    let (inputs, request) = fixture();
    let original = adapter::simulate(&inputs, &request).unwrap();
    let mut forged = original.clone();
    forged.reference.candidates = (0..8)
        .map(|id| {
            let mut c = original.reference.candidates[0].clone();
            c.id = id;
            c
        })
        .collect();
    forged.reference.combinations = 8;
    forged.reference.graph_plans = 8;
    forged.reference.materializations = 8;
    forged.reference.materialized_bases = original.reference.materialized_bases * 8;
    forged.documents.clear();
    for candidate_id in 0..8 {
        for doc in &original.documents {
            let mut doc = doc.clone();
            doc.selection.candidate_id = candidate_id;
            let mut recipe: UnverifiedRecipe =
                serde_json::from_str(&doc.molecule.metadata[RECIPE_KEY]).unwrap();
            recipe.selection = doc.selection.clone();
            doc.molecule
                .metadata
                .insert(RECIPE_KEY.into(), serde_json::to_string(&recipe).unwrap());
            doc.annotation_evidence[0].unmapped_features = (0..6)
                .map(|_| Feature {
                    kind: "misc_feature".into(),
                    location: Location::range(0, 1),
                    qualifiers: vec![("note".into(), Some("x".repeat(1024 * 1024)))],
                    needs_review: true,
                })
                .collect();
            // Every individual row fits borrowed/tree/32MiB shape limits. Together
            // they exceed96MiB while raw reference DNA remains below40k bases.
            doc.validate().unwrap();
            forged.documents.push(doc);
        }
    }
    forged.reference.validate().unwrap();
    assert!(matches!(
        forged.validate(),
        Err(collagene_model::Error::TooLarge)
    ));
}

#[test]
fn report_current_binding_rejects_consistent_foreign_source_and_recipe_policy() {
    let (inputs, request) = fixture();
    let report = adapter::simulate(&inputs, &request).unwrap();
    let mut foreign = report.clone();
    foreign.documents[0].annotation_evidence[0].source.role =
        "self-consistent forged source".into();
    let mut recipe: UnverifiedRecipe =
        serde_json::from_str(&foreign.documents[0].molecule.metadata[RECIPE_KEY]).unwrap();
    recipe.source_identities[0].role = "self-consistent forged source".into();
    foreign.documents[0]
        .molecule
        .metadata
        .insert(RECIPE_KEY.into(), serde_json::to_string(&recipe).unwrap());
    foreign.validate().unwrap();
    assert!(foreign.validate_bindings_for(&inputs, &request).is_err());
    let mut foreign = report.clone();
    let mut recipe: UnverifiedRecipe =
        serde_json::from_str(&foreign.documents[0].molecule.metadata[RECIPE_KEY]).unwrap();
    recipe.request.annotation_limits.max_work -= 1;
    foreign.documents[0]
        .molecule
        .metadata
        .insert(RECIPE_KEY.into(), serde_json::to_string(&recipe).unwrap());
    foreign.validate().unwrap();
    assert!(foreign.validate_bindings_for(&inputs, &request).is_err());
    // Changing evidence and recipe revision together must also fail, even though
    // the unchanged outer/raw report still binds the actual scientific inputs.
    let mut foreign = report.clone();
    foreign.documents[0].annotation_evidence[1].source.revision = "0".repeat(64);
    let mut recipe: UnverifiedRecipe =
        serde_json::from_str(&foreign.documents[0].molecule.metadata[RECIPE_KEY]).unwrap();
    recipe.source_identities[1].revision = "0".repeat(64);
    foreign.documents[0]
        .molecule
        .metadata
        .insert(RECIPE_KEY.into(), serde_json::to_string(&recipe).unwrap());
    foreign.validate().unwrap();
    assert!(foreign.validate_bindings_for(&inputs, &request).is_err());
    report.validate_bindings_for(&inputs, &request).unwrap();
}

#[test]
fn standalone_dna_identity_target_and_annotation_counters_cannot_be_replaced() {
    let (inputs, request) = fixture();
    let report = adapter::simulate(&inputs, &request).unwrap();
    let document = &report.documents[0];
    let policy = reference::TargetPolicy {
        allow_reverse_complement: false,
    };
    let checked = adapter::checked_product(
        &inputs,
        &request,
        &document.selection,
        &document.molecule,
        &policy,
    )
    .unwrap();
    let mutate = |m: &mut Molecule| {
        let mut b = m.sequence.as_str().as_bytes().to_vec();
        b[20] = if b[20] == b'A' { b'C' } else { b'A' };
        m.sequence = Sequence::dna(std::str::from_utf8(&b).unwrap()).unwrap();
    };
    let mut wrong_alphabet = document.clone();
    wrong_alphabet.molecule.sequence = Sequence::new(
        collagene_model::Alphabet::Protein,
        document.molecule.sequence.as_str(),
    )
    .unwrap();
    assert!(wrong_alphabet.validate().is_err());
    let mut substituted = document.clone();
    mutate(&mut substituted.molecule);
    assert!(substituted.validate().is_err());
    assert!(
        substituted
            .validate_bindings_for(&inputs, &request, &document.selection)
            .is_err()
    );
    let mut substituted = checked.clone();
    mutate(&mut substituted.product.molecule);
    assert!(
        substituted
            .validate_bindings_for(
                &inputs,
                &request,
                &document.selection,
                &document.molecule,
                &policy
            )
            .is_err()
    );
    // Even a self-consistent replacement DNA identity/recipe is only an imported
    // claim. Bind the ACTUAL selected reference DNA to the independent target.
    let new_raw = Molecule::dna(
        format!(
            "reciprocal product {}",
            substituted.product.selection.component_id + 1
        ),
        substituted.product.molecule.sequence.as_str(),
        substituted.product.molecule.topology,
    )
    .unwrap()
    .revision();
    substituted.product.reference_product_revision = new_raw.clone();
    substituted.product.reciprocal_components[substituted.product.selection.component_id]
        .reference_revision = new_raw.clone();
    let mut recipe: UnverifiedRecipe =
        serde_json::from_str(&substituted.product.molecule.metadata[RECIPE_KEY]).unwrap();
    recipe.reference_product_revision = new_raw;
    recipe.reciprocal_components = substituted.product.reciprocal_components.clone();
    substituted
        .product
        .molecule
        .metadata
        .insert(RECIPE_KEY.into(), serde_json::to_string(&recipe).unwrap());
    substituted.validate().unwrap();
    assert!(
        substituted
            .validate_bindings_for(
                &inputs,
                &request,
                &document.selection,
                &document.molecule,
                &policy
            )
            .is_err()
    );
    let (mut annotated, request) = fixture();
    annotated[0]
        .molecule
        .features
        .push(feature("payload", Location::range(25, 34)));
    let report = adapter::simulate(&annotated, &request).unwrap();
    let mut false_budget = report.clone();
    false_budget.annotation_features = 0;
    assert!(
        false_budget
            .validate_bindings_for(&annotated, &request)
            .is_err()
    );
    let mut doc = report
        .documents
        .iter()
        .find(|d| !d.molecule.features.is_empty())
        .unwrap()
        .clone();
    doc.annotation_evidence[0].mapped_features = 0;
    assert!(doc.validate().is_err());
    let mut doc = report
        .documents
        .iter()
        .find(|d| d.molecule.features.is_empty())
        .unwrap()
        .clone();
    assert_eq!(doc.annotation_evidence[0].dropped_features, 1);
    doc.annotation_evidence[0].dropped_features = 0;
    doc.validate().unwrap();
    assert!(
        doc.validate_bindings_for(&annotated, &request, &doc.selection)
            .is_err()
    );
    let mut forged = report.clone();
    let index = forged
        .documents
        .iter()
        .position(|d| d.molecule.features.is_empty())
        .unwrap();
    forged.documents[index] = doc;
    assert!(forged.validate_bindings_for(&annotated, &request).is_err());
    annotated[0].molecule.features = vec![feature("co-source", Location::range(8, 15))];
    let report = adapter::simulate(&annotated, &request).unwrap();
    let mut doc = report
        .documents
        .iter()
        .find(|d| d.molecule.features.iter().any(|f| f.needs_review))
        .unwrap()
        .clone();
    for e in &mut doc.annotation_evidence {
        e.reviewed_features = 0;
    }
    assert!(doc.validate().is_err());
}
#[test]
fn checked_document_literal_comparison_respects_circle_origin_and_explicit_rc_policy() {
    let (inputs, request) = fixture();
    let report = adapter::simulate(&inputs, &request).unwrap();
    let doc = report
        .documents
        .iter()
        .find(|d| d.molecule.topology == Topology::Circular)
        .unwrap();
    let mut target = doc.molecule.clone();
    let s = target.sequence.as_str();
    target.sequence = Sequence::dna(&format!("{}{}", &s[17..], &s[..17])).unwrap();
    let strict = reference::TargetPolicy {
        allow_reverse_complement: false,
    };
    let checked =
        adapter::checked_product(&inputs, &request, &doc.selection, &target, &strict).unwrap();
    checked
        .validate_bindings_for(&inputs, &request, &doc.selection, &target, &strict)
        .unwrap();
    target.sequence = target.sequence.reverse_complement().unwrap();
    assert!(adapter::checked_product(&inputs, &request, &doc.selection, &target, &strict).is_err());
    let explicit = reference::TargetPolicy {
        allow_reverse_complement: true,
    };
    let checked =
        adapter::checked_product(&inputs, &request, &doc.selection, &target, &explicit).unwrap();
    checked
        .validate_bindings_for(&inputs, &request, &doc.selection, &target, &explicit)
        .unwrap();
}

#[test]
fn target_no_match_flags_cannot_contradict_results_or_actual_ambiguous_dna() {
    let (inputs, request) = fixture();
    let report = adapter::simulate(&inputs, &request).unwrap();
    let policy = reference::TargetPolicy {
        allow_reverse_complement: false,
    };
    let valid =
        adapter::check_target(&inputs, &request, &report.documents[0].molecule, &policy).unwrap();
    let mut contradictory = valid.clone();
    contradictory.reference.no_match_proven_under_model = true;
    assert!(contradictory.validate().is_err());
    let target = Molecule::dna(
        "independent uncertain requirement",
        "NNNNN",
        Topology::Linear,
    )
    .unwrap();
    let uncertain = adapter::check_target(&inputs, &request, &target, &policy).unwrap();
    assert!(uncertain.reference.ambiguous_target_dna);
    let mut contradictory = uncertain.clone();
    contradictory.reference.no_match_proven_under_model = true;
    assert!(contradictory.validate().is_err());
    let mut coherent_false_claim = uncertain.clone();
    coherent_false_claim.reference.ambiguous_target_dna = false;
    coherent_false_claim.reference.applicability =
        reference::TargetApplicability::DeclaredFullSiteProfileOnly;
    coherent_false_claim.reference.no_match_proven_under_model = true;
    coherent_false_claim.validate().unwrap();
    assert!(
        coherent_false_claim
            .validate_bindings_for(&inputs, &request, &target, &policy)
            .is_err()
    );
    let mut ambiguous_inputs = inputs.clone();
    let s = ambiguous_inputs[0]
        .molecule
        .sequence
        .as_str()
        .replacen("ATGAAATAA", "ATGNAATAA", 1);
    ambiguous_inputs[0].molecule.sequence = Sequence::dna(&s).unwrap();
    let exact_target = Molecule::dna(
        "independent exact requirement",
        "ACGTACGT",
        Topology::Linear,
    )
    .unwrap();
    let uncertain =
        adapter::check_target(&ambiguous_inputs, &request, &exact_target, &policy).unwrap();
    assert!(uncertain.reference.ambiguous_source_dna);
    let mut false_claim = uncertain;
    false_claim.reference.ambiguous_source_dna = false;
    false_claim.reference.applicability =
        reference::TargetApplicability::DeclaredFullSiteProfileOnly;
    false_claim.reference.no_match_proven_under_model = true;
    false_claim.validate().unwrap();
    assert!(
        false_claim
            .validate_bindings_for(&ambiguous_inputs, &request, &exact_target, &policy)
            .is_err()
    );
    let absent = adapter::check_target(&inputs, &request, &exact_target, &policy).unwrap();
    assert!(absent.reference.no_match_proven_under_model);
    let mut contradictory = absent.clone();
    contradictory.reference.search_complete = false;
    contradictory.reference.applicability = reference::TargetApplicability::Unresolved;
    assert!(contradictory.validate().is_err());
    let mut contradictory = absent;
    contradictory.reference.possible_sites = true;
    contradictory.reference.no_match_proven_under_model = false;
    contradictory.reference.applicability = reference::TargetApplicability::Unresolved;
    assert!(contradictory.reference.unique_under_declared_profile);
    assert!(contradictory.validate().is_err());
}
