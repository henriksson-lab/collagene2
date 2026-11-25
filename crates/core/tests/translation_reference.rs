use collagene_core::{
    Alphabet, Feature, Location, Molecule, Sequence, Topology,
    translation::{InitiationPolicy, TranslationOptions},
    translation_reference::{self as reference, Request, Selection},
};

fn request(bases: &str) -> Request {
    Request {
        source: Molecule::dna("actual source", bases, Topology::Linear).unwrap(),
        selection: Selection::Region {
            location: Location::range(0, bases.len()),
            options: TranslationOptions::default(),
        },
    }
}
fn cds(source: &mut Molecule, location: Location, qualifiers: Vec<(&str, &str)>) {
    source.features.push(Feature {
        kind: "CDS".into(),
        location,
        qualifiers: qualifiers
            .into_iter()
            .map(|(key, value)| (key.into(), Some(value.into())))
            .collect(),
        needs_review: false,
    });
}

#[test]
fn literal_region_keeps_stops_and_frame_skip_and_trailing_findings() {
    let mut r = request("AATGAAATAACC");
    let Selection::Region { options, .. } = &mut r.selection else {
        unreachable!()
    };
    options.frame = 1;
    let report = reference::inspect(&r).unwrap();
    assert_eq!(report.protein, "MK*");
    assert_eq!(report.translation.trailing_bases, 2);
    assert_eq!(report.translation.translated_codons, 3);
    assert_eq!(report.translation.nucleotide_length, 12);
    assert!(report.cds.is_none());
    assert!(!report.warnings.is_empty());
    report
        .validate_binding(&reference::inspect_revision(&r).unwrap())
        .unwrap();
    let outputs = reference::prepare(&r, "frame 2 protein").unwrap();
    assert_eq!(outputs.protein.sequence.as_str(), "MK*");
    let annotation = outputs.annotation.unwrap();
    assert!(annotation.needs_review);
    assert!(
        annotation
            .qualifiers
            .contains(&("codon_start".into(), Some("2".into())))
    );
    assert!(
        !annotation
            .qualifiers
            .iter()
            .any(|(key, _)| key == "translation" || key == "partial")
    );
}

#[test]
fn compound_extraction_and_global_reverse_follow_biological_order() {
    let mut r = request("CCCATGAAATAAGGG");
    let joined = Location::Join(vec![Location::range(3, 6), Location::range(6, 12)]);
    r.selection = Selection::Region {
        location: joined.clone(),
        options: TranslationOptions::default(),
    };
    assert_eq!(reference::inspect(&r).unwrap().protein, "MK*");
    let reversed = Location::Complement(Box::new(joined.clone()));
    r.selection = Selection::Region {
        location: reversed.clone(),
        options: TranslationOptions::default(),
    };
    assert_eq!(reference::inspect(&r).unwrap().protein, "LFH");
    let Selection::Region { options, .. } = &mut r.selection else {
        unreachable!()
    };
    options.reverse = true;
    let report = reference::inspect(&r).unwrap();
    assert_eq!(report.protein, "MK*");
    assert_eq!(report.effective_location, joined);
    assert!(report.translation.options.reverse);
    assert_eq!(
        reference::prepare(&r, "joined")
            .unwrap()
            .annotation
            .unwrap()
            .location,
        joined
    );
}

#[test]
fn circular_origin_selection_and_reverse_use_explicit_ordered_spans() {
    let mut r = request("AAATAACCCATG");
    r.source.topology = Topology::Circular;
    let location = Location::Join(vec![Location::range(9, 12), Location::range(0, 6)]);
    r.selection = Selection::Region {
        location: location.clone(),
        options: TranslationOptions::default(),
    };
    assert_eq!(reference::inspect(&r).unwrap().protein, "MK*");
    let Selection::Region { options, .. } = &mut r.selection else {
        unreachable!()
    };
    options.reverse = true;
    let report = reference::inspect(&r).unwrap();
    assert_eq!(report.protein, "LFH");
    assert_eq!(
        report.effective_location,
        Location::Complement(Box::new(location))
    );
    assert_eq!(
        reference::prepare(&r, "reverse origin")
            .unwrap()
            .protein
            .topology,
        Topology::Linear
    );
}

#[test]
fn exact_repeated_parts_are_bounded_but_not_silently_unioned() {
    let mut r = request("ATG");
    r.selection = Selection::Region {
        location: Location::Join(vec![Location::range(0, 3), Location::range(0, 3)]),
        options: TranslationOptions::default(),
    };
    let report = reference::inspect(&r).unwrap();
    assert_eq!(report.protein, "MM");
    assert_eq!(report.translation.nucleotide_length, 6);
}

#[test]
fn cds_uses_actual_source_qualifiers_and_omits_one_terminal_stop() {
    let mut r = request("GTGAAATAA");
    cds(
        &mut r.source,
        Location::range(0, 9),
        vec![
            ("label", "native CDS"),
            ("transl_table", "11"),
            ("translation", "QQ"),
        ],
    );
    r.selection = Selection::Cds { feature_index: 0 };
    let report = reference::inspect(&r).unwrap();
    assert_eq!(report.protein, "MK");
    assert_eq!(report.translation.protein, "MK*");
    let cds = report.cds.as_ref().unwrap();
    assert_eq!(cds.stored_translation_matches, Some(false));
    assert!(!cds.complete_under_model);
    assert!(
        report
            .warnings
            .iter()
            .any(|warning| warning.contains("Stored /translation differs"))
    );
    let outputs = reference::prepare(&r, "linked CDS protein").unwrap();
    assert_eq!(outputs.protein.sequence.as_str(), "MK");
    assert!(outputs.annotation.is_none());
    outputs.validate_sources(&r, "linked CDS protein").unwrap();
}

#[test]
fn region_annotation_keeps_literal_initiation_policy_without_inventing_partial_cds() {
    let mut r = request("GTGAAATAA");
    let Selection::Region { options, .. } = &mut r.selection else {
        unreachable!()
    };
    options.genetic_code = 11;
    options.initiation = InitiationPolicy::Ordinary;
    let outputs = reference::prepare(&r, "literal code 11").unwrap();
    assert_eq!(outputs.protein.sequence.as_str(), "VK*");
    let proposed = outputs.annotation.unwrap();
    assert!(proposed.needs_review);
    assert!(proposed.qualifiers.iter().any(|(key, value)| key == "note"
        && value.as_ref().is_some_and(
            |text| text.contains("Ordinary") && text.contains("semantics can differ")
        )));
    assert!(
        !proposed
            .qualifiers
            .iter()
            .any(|(key, _)| key == "partial" || key == "translation")
    );
    assert!(collagene_core::translation::translate_cds(&r.source, &proposed).is_err());
}

#[test]
fn rna_ambiguity_and_internal_stops_remain_explicit() {
    let mut r = request("ATGNNNTAAAC");
    r.source.sequence = Sequence::new(Alphabet::Rna, "AUGNNNUAAAC").unwrap();
    let report = reference::inspect(&r).unwrap();
    assert_eq!(report.protein, "MX*");
    assert_eq!(report.translation.ambiguous_codons, 1);
    assert_eq!(report.translation.unresolved_codons, 1);
    assert_eq!(report.translation.possible_stop_codons, 1);
    assert_eq!(report.translation.trailing_bases, 2);
    let mut internal = request("ATGTAATAG");
    cds(&mut internal.source, Location::range(0, 9), vec![]);
    internal.selection = Selection::Cds { feature_index: 0 };
    let report = reference::inspect(&internal).unwrap();
    assert_eq!(report.protein, "M*");
    assert_eq!(report.cds.as_ref().unwrap().internal_stops, 1);
    assert!(!report.cds.unwrap().complete_under_model);
}

#[test]
fn protein_output_has_only_scoped_provenance_and_clean_document_authority() {
    let mut r = request("ATGAAATAA");
    r.source.original_genbank = Some("original DNA header".into());
    r.source
        .metadata
        .insert("duplex_model".into(), "unrelated physical evidence".into());
    r.source
        .metadata
        .insert("trace:source".into(), "unrelated trace".into());
    let original = r.source.clone();
    let output = reference::prepare(&r, "new protein").unwrap();
    assert_eq!(r.source, original);
    assert_eq!(output.protein.sequence.alphabet(), Alphabet::Protein);
    assert!(output.protein.features.is_empty());
    assert!(output.protein.original_genbank.is_none());
    assert_eq!(output.protein.metadata.len(), 5);
    assert_eq!(
        output.protein.metadata.get("translation:source_revision"),
        Some(&r.source.revision())
    );
    assert_eq!(
        output.context_revision,
        reference::output_revision(&r, "new protein").unwrap()
    );
    output.validate_binding(&output.context_revision).unwrap();
}

#[test]
fn tampered_reports_outputs_and_source_bindings_are_rejected() {
    let r = request("ATGAAATAA");
    let report = reference::inspect(&r).unwrap();
    let mut tampered = report.clone();
    tampered.protein = "QQ*".into();
    assert!(tampered.validate_binding(&report.context_revision).is_err());
    tampered.translation.protein = "QQ*".into();
    // Matching lengths/context are not amino-acid authority. Checking actual
    // source calls detects coordinated replacement of both text fields.
    tampered.validate_binding(&report.context_revision).unwrap();
    assert!(tampered.validate_sources(&r).is_err());
    let mut changed = r.clone();
    changed.source = Molecule::dna("changed", "ATGCCCTAA", Topology::Linear).unwrap();
    assert!(report.validate_sources(&changed).is_err());
    let mut output = reference::prepare(&r, "protein").unwrap();
    output
        .protein
        .metadata
        .insert("duplex_model".into(), "must not inherit".into());
    assert!(output.validate_binding(&output.context_revision).is_err());
    let mut output = reference::prepare(&r, "protein").unwrap();
    output.annotation.as_mut().unwrap().needs_review = false;
    assert!(output.validate_binding(&output.context_revision).is_err());
}

#[test]
fn fuzzy_order_unknown_bad_frame_empty_and_protein_inputs_are_explicit_errors() {
    let original = request("ATGAAATAA");
    for location in [
        Location::Range {
            start: 0,
            end: 9,
            before: true,
            after: false,
        },
        Location::Order(vec![Location::range(0, 9)]),
        Location::Unresolved("remote:1..9".into()),
        Location::Join(vec![]),
        Location::range(3, 3),
        Location::range(0, 10),
    ] {
        let mut r = original.clone();
        r.selection = Selection::Region {
            location,
            options: TranslationOptions::default(),
        };
        assert!(reference::inspect_revision(&r).is_err());
        assert!(reference::inspect(&r).is_err());
    }
    let mut r = original.clone();
    let Selection::Region { options, .. } = &mut r.selection else {
        unreachable!()
    };
    options.frame = 3;
    assert!(reference::inspect(&r).is_err());
    let mut r = original.clone();
    r.source.sequence = Sequence::new(Alphabet::Protein, "MK").unwrap();
    assert!(reference::inspect(&r).is_err());
    let short = request("AT");
    assert!(reference::inspect(&short).unwrap().protein.is_empty());
    assert!(reference::prepare(&short, "no residues").is_err());
    let mut stop = request("TAA");
    cds(&mut stop.source, Location::range(0, 3), vec![]);
    stop.selection = Selection::Cds { feature_index: 0 };
    assert!(reference::inspect(&stop).unwrap().protein.is_empty());
    assert!(reference::prepare(&stop, "stop only CDS").is_err());
    assert!(reference::prepare(&original, "\n").is_err());
    assert!(reference::output_revision(&original, &"x".repeat(1025)).is_err());
}

#[test]
fn oversized_borrowed_models_locations_and_unknown_fields_are_rejected() {
    let mut r = request("ATG");
    r.source.features = vec![
        Feature {
            kind: "CDS".into(),
            location: Location::range(0, 3),
            qualifiers: vec![],
            needs_review: false
        };
        4097
    ];
    assert!(reference::inspect_revision(&r).is_err());
    let mut r = request("ATG");
    r.selection = Selection::Region {
        location: Location::Join(vec![Location::range(0, 3); 16385]),
        options: TranslationOptions::default(),
    };
    assert!(reference::inspect_revision(&r).is_err());
    let mut value = serde_json::to_value(request("ATG")).unwrap();
    value["verified"] = serde_json::json!(true);
    assert!(serde_json::from_value::<Request>(value).is_err());
}
