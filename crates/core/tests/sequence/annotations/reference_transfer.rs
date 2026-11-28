use collagene_core::sequence::reference_transfer::{
    self as transfer, Part, RequestCopy, RequestPaste, Selection,
};
use collagene_model::{Alphabet, Error, Feature, Location, Molecule, Sequence, Topology};
fn feature(kind: &str, location: Location, label: &str) -> Feature {
    Feature {
        kind: kind.into(),
        location,
        qualifiers: vec![("label".into(), Some(label.into()))],
        needs_review: false,
    }
}
fn selected(start: usize, end: usize, reverse: bool) -> Selection {
    Selection {
        parts: vec![Part {
            start,
            end,
            reverse: false,
        }],
        reverse,
    }
}
fn paste_request(
    source: Molecule,
    destination: Molecule,
    selection: Selection,
    start: usize,
    end: usize,
) -> RequestPaste {
    let fragment_revision = transfer::copy_revision(&RequestCopy {
        source: source.clone(),
        selection: selection.clone(),
    })
    .unwrap();
    RequestPaste {
        source,
        destination,
        selection,
        fragment_revision,
        start,
        end,
    }
}
#[test]
fn reverse_fragment_paste_projects_source_and_destination_features_without_source_mutation() {
    let mut source = Molecule::dna("original", "AACGTTAC", Topology::Linear).unwrap();
    source
        .features
        .push(feature("misc_feature", Location::range(3, 5), "copied"));
    source.features[0]
        .qualifiers
        .push(("CiTaTiOn".into(), Some("[2]".into())));
    source
        .metadata
        .insert("duplex_model".into(), "original physical descriptor".into());
    let mut destination = Molecule::dna("destination", "TTAAAACC", Topology::Linear).unwrap();
    destination
        .features
        .push(feature("CDS", Location::range(1, 5), "touched"));
    destination.features[0]
        .qualifiers
        .push(("translation".into(), Some("OLD".into())));
    destination.features[0].qualifiers.extend([
        ("TRANSLATION".into(), Some("STALE".into())),
        ("CODON_START".into(), Some("2".into())),
    ]);
    destination
        .features
        .push(feature("misc_feature", Location::range(4, 8), "shifted"));
    let source_before = source.clone();
    let destination_before = destination.clone();
    let request = paste_request(source, destination, selected(2, 6, true), 2, 4);
    let copied = transfer::copy(&RequestCopy {
        source: request.source.clone(),
        selection: request.selection.clone(),
    })
    .unwrap();
    assert_eq!(copied.molecule.sequence.as_str(), "AACG");
    assert_eq!(
        copied.molecule.features[0].location,
        Location::Complement(Box::new(Location::range(1, 3)))
    );
    let result = transfer::paste(&request).unwrap();
    assert_eq!(result.product.sequence.as_str(), "TTAACGAACC");
    assert_eq!((result.inserted_start, result.inserted_end), (2, 6));
    assert_eq!(result.product.features[0].location, Location::range(1, 7));
    assert!(result.product.features[0].needs_review);
    assert!(
        !result.product.features[0]
            .qualifiers
            .iter()
            .any(|(key, _)| key == "translation")
    );
    assert_eq!(result.product.features[1].location, Location::range(6, 10));
    assert_eq!(
        result.product.features[2].location,
        Location::Complement(Box::new(Location::range(3, 5)))
    );
    assert_eq!(
        result.product.features[2]
            .location
            .extract(&result.product.sequence)
            .unwrap(),
        "GT"
    );
    assert!(
        result.product.features[2]
            .qualifiers
            .contains(&("source_citation".into(), Some("[2]".into())))
    );
    assert!(
        !result.product.features[2]
            .qualifiers
            .iter()
            .any(|(key, _)| key == "citation")
    );
    assert!(result.product.metadata["duplex_model"].contains("Reference-only"));
    assert!(!result.product.metadata.contains_key("verification_status"));
    assert!(
        result
            .product
            .metadata
            .values()
            .any(|value| value.contains("original physical descriptor"))
    );
    assert_eq!(request.source, source_before);
    assert_eq!(request.destination, destination_before);
    result
        .validate_binding(&transfer::paste_revision(&request).unwrap())
        .unwrap();
    let reopened = collagene_core::protocol::load_document(
        &collagene_core::protocol::save_document(&result.product).unwrap(),
    )
    .unwrap();
    assert_eq!(reopened, result.product);
}
#[test]
fn circular_and_compound_biological_order_are_explicit_and_reversible() {
    let mut source = Molecule::dna("circle", "AACCGGTA", Topology::Circular).unwrap();
    let location = Location::Join(vec![Location::range(6, 8), Location::range(0, 2)]);
    source
        .features
        .push(feature("misc_feature", location.clone(), "origin feature"));
    let selection = Selection::from_location(&location, 8).unwrap();
    let request = RequestCopy {
        source: source.clone(),
        selection,
    };
    let fragment = transfer::copy(&request).unwrap();
    assert_eq!(fragment.molecule.sequence.as_str(), "TAAA");
    assert_eq!(
        fragment.molecule.features[0]
            .location
            .extract(&fragment.molecule.sequence)
            .unwrap(),
        "TAAA"
    );
    let mut reversed = request.clone();
    reversed.selection.reverse = true;
    let reverse = transfer::copy(&reversed).unwrap();
    assert_eq!(reverse.molecule.sequence.as_str(), "TTTA");
    assert_ne!(fragment.context_revision, reverse.context_revision);
    let compound = Location::Complement(Box::new(Location::Join(vec![
        Location::range(0, 2),
        Location::range(4, 6),
    ])));
    let selection = Selection::from_location(&compound, 8).unwrap();
    assert_eq!(
        selection.parts,
        vec![
            Part {
                start: 4,
                end: 6,
                reverse: true
            },
            Part {
                start: 0,
                end: 2,
                reverse: true
            }
        ]
    );
    let fragment = transfer::copy(&RequestCopy {
        source: source.clone(),
        selection,
    })
    .unwrap();
    assert_eq!(fragment.molecule.sequence.as_str(), "CCTT");
    assert_eq!(
        fragment.molecule.sequence.as_str(),
        compound.extract(&source.sequence).unwrap()
    );
}
#[test]
fn partial_unknown_annotations_and_source_evidence_remain_distinct() {
    let mut source = Molecule::dna("annotations", "AACCGGTT", Topology::Linear).unwrap();
    source
        .features
        .push(feature("CDS", Location::range(1, 7), "partial gene"));
    source.features[0].qualifiers.extend([
        ("translation".into(), Some("ASSERTED".into())),
        ("codon_start".into(), Some("1".into())),
    ]);
    source.features.push(feature(
        "misc_feature",
        Location::Unresolved("external:1..5".into()),
        "unknown",
    ));
    source
        .features
        .push(feature("misc_feature", Location::range(0, 1), "outside"));
    let destination = Molecule::dna("target", "ATAT", Topology::Linear).unwrap();
    let request = paste_request(source, destination, selected(3, 6, false), 2, 2);
    let copied = transfer::copy(&RequestCopy {
        source: request.source.clone(),
        selection: request.selection.clone(),
    })
    .unwrap();
    assert_eq!(copied.molecule.sequence.as_str(), "CGG");
    assert_eq!(copied.dropped_features, 1);
    assert_eq!(copied.unmapped_features.len(), 1);
    assert_eq!(copied.molecule.features.len(), 1);
    assert!(copied.molecule.features[0].needs_review);
    assert!(
        !copied.molecule.features[0]
            .qualifiers
            .iter()
            .any(|(key, _)| key == "translation" || key == "codon_start")
    );
    let output = transfer::paste(&request).unwrap();
    assert_eq!(output.product.sequence.as_str(), "ATCGGAT");
    assert_eq!(output.copied_features, 1);
    assert_eq!(output.unmapped_features.len(), 1);
    assert!(
        output
            .product
            .metadata
            .values()
            .any(|value| value.contains("ASSERTED"))
    );
    assert!(
        output
            .product
            .metadata
            .values()
            .any(|value| value.contains("external:1..5"))
    );
}
#[test]
fn unsupported_and_stale_requests_never_modify_either_document() {
    let source = Molecule::dna("source", "AACCGGTT", Topology::Linear).unwrap();
    let destination = Molecule::dna("destination", "ATAT", Topology::Linear).unwrap();
    let request = paste_request(
        source.clone(),
        destination.clone(),
        selected(2, 6, false),
        1,
        2,
    );
    let mut forged = transfer::copy(&RequestCopy {
        source: source.clone(),
        selection: request.selection.clone(),
    })
    .unwrap();
    forged.molecule.sequence = Sequence::dna("AAAA").unwrap();
    let output = transfer::paste(&request).unwrap();
    assert_eq!(output.product.sequence.as_str(), "ACCGGAT"); // copied CCGG, not forged AAAA
    let mut stale = request.clone();
    stale.fragment_revision = "0".repeat(64);
    assert!(transfer::paste(&stale).is_err());
    let mut unsupported = request.clone();
    unsupported
        .source
        .features
        .push(feature("CDS", Location::range(0, 8), "special"));
    unsupported.source.features[0]
        .qualifiers
        .push(("transl_except".into(), Some("(pos:2..4,aa:Sec)".into())));
    assert!(
        transfer::copy(&RequestCopy {
            source: unsupported.source,
            selection: unsupported.selection
        })
        .is_err()
    );
    let mut bad_destination = request.clone();
    bad_destination.destination.features.push(feature(
        "misc_feature",
        Location::Unresolved("unknown".into()),
        "unresolved",
    ));
    let before = bad_destination.destination.clone();
    assert!(transfer::paste(&bad_destination).is_err());
    assert_eq!(bad_destination.destination, before);
    let mut uppercase = request.clone();
    uppercase.destination.features.push(feature(
        "CDS",
        Location::range(0, 4),
        "coordinate qualifier",
    ));
    uppercase.destination.features[0]
        .qualifiers
        .push(("TRANSL_EXCEPT".into(), Some("(pos:1..3,aa:Sec)".into())));
    assert!(transfer::paste_revision(&uppercase).is_err());
    let overlap = Selection {
        parts: vec![
            Part {
                start: 0,
                end: 4,
                reverse: false,
            },
            Part {
                start: 2,
                end: 5,
                reverse: false,
            },
        ],
        reverse: false,
    };
    assert!(
        transfer::copy(&RequestCopy {
            source: source.clone(),
            selection: overlap
        })
        .is_err()
    );
    assert!(Selection::from_location(&Location::Order(vec![Location::range(0, 2)]), 8).is_err());
    assert!(
        Selection::from_location(
            &Location::Range {
                start: 0,
                end: 2,
                before: true,
                after: false
            },
            8
        )
        .is_err()
    );
    let many = Selection {
        parts: vec![
            Part {
                start: 0,
                end: 1,
                reverse: false
            };
            4097
        ],
        reverse: false,
    };
    assert!(matches!(
        transfer::copy_revision(&RequestCopy {
            source: source.clone(),
            selection: many
        }),
        Err(Error::Unsupported(_))
    ));
    assert_eq!(request.source, source);
    assert_eq!(request.destination, destination);
}
#[test]
fn supported_alphabets_and_reverse_restrictions_are_not_reinterpreted() {
    let mut rna = Molecule::dna("RNA", "AATTCC", Topology::Linear).unwrap();
    rna.sequence = Sequence::new(Alphabet::Rna, "AAUUCC").unwrap();
    let fragment = transfer::copy(&RequestCopy {
        source: rna,
        selection: selected(2, 4, true),
    })
    .unwrap();
    assert_eq!(fragment.molecule.sequence.alphabet(), Alphabet::Rna);
    assert_eq!(fragment.molecule.sequence.as_str(), "AA");
    let mut protein = Molecule::dna("protein", "ATG", Topology::Linear).unwrap();
    protein.sequence = Sequence::new(Alphabet::Protein, "MKL").unwrap();
    let copied = transfer::copy(&RequestCopy {
        source: protein.clone(),
        selection: selected(0, 2, false),
    })
    .unwrap();
    assert_eq!(copied.molecule.sequence.as_str(), "MK");
    assert!(
        transfer::copy(&RequestCopy {
            source: protein,
            selection: selected(0, 2, true)
        })
        .is_err()
    );
}

#[test]
fn final_circular_projection_uses_actual_continuity_not_linear_preview_flags() {
    let mut source = Molecule::dna("circle", "GGCCAAAT", Topology::Circular).unwrap();
    let original = Location::Join(vec![Location::range(6, 8), Location::range(0, 4)]);
    let mut cds = feature("CDS", original.clone(), "across origin");
    cds.qualifiers.extend([
        ("translation".into(), Some("MA".into())),
        ("codon_start".into(), Some("1".into())),
    ]);
    source.features.push(cds);
    let request = paste_request(
        source,
        Molecule::dna("circular destination", "TTTTTTTT", Topology::Circular).unwrap(),
        selected(0, 8, false),
        0,
        8,
    );
    let linear = transfer::copy(&RequestCopy {
        source: request.source.clone(),
        selection: request.selection.clone(),
    })
    .unwrap();
    assert!(linear.molecule.features[0].needs_review);
    assert!(
        !linear.molecule.features[0]
            .qualifiers
            .iter()
            .any(|(key, _)| key == "translation")
    );
    let result = transfer::paste(&request).unwrap();
    assert_eq!(result.product.sequence.as_str(), "GGCCAAAT");
    assert_eq!(result.product.features[0].location, original);
    assert!(!result.product.features[0].needs_review);
    assert_eq!(
        result.product.features[0]
            .location
            .extract(&result.product.sequence)
            .unwrap(),
        "ATGGCC"
    );
    assert!(
        result.product.features[0]
            .qualifiers
            .contains(&("translation".into(), Some("MA".into())))
    );
    assert_eq!(result.copied_features, 1);
    // Inserting the same complete linear clipboard letters into a larger circle
    // does not make the gene's former origin adjacency continuous.
    let larger = paste_request(
        request.source.clone(),
        Molecule::dna("larger circle", "TTTTTTTTTTTT", Topology::Circular).unwrap(),
        request.selection.clone(),
        2,
        2,
    );
    let interrupted = transfer::paste(&larger).unwrap();
    assert_eq!(
        interrupted.product.sequence.as_str(),
        "TTGGCCAAATTTTTTTTTTT"
    );
    assert!(interrupted.product.features[0].needs_review);
    assert!(
        !interrupted.product.features[0]
            .qualifiers
            .iter()
            .any(|(key, _)| key == "translation" || key == "codon_start")
    );
}
