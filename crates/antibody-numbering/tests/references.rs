use collagene_antibody_numbering::references::{Corpus, Species, load};
#[test]
fn pinned_sources_keep_explicit_numbering_and_complete_inclusion_exclusion_ledgers() {
    let corpus = load(&Species::Human).unwrap();
    let sources = corpus.sources();
    assert_eq!(
        sources
            .iter()
            .map(|s| (
                s.locus.as_str(),
                s.release_version,
                s.original_records,
                s.v_records,
                s.j_records
            ))
            .collect::<Vec<_>>(),
        vec![
            ("IGH", 10, 240, 194, 7),
            ("IGK", 5, 84, 74, 7),
            ("IGL", 4, 97, 86, 9)
        ]
    );
    assert_eq!(corpus.v_entries().len(), 354);
    let atypical = corpus
        .v_entries()
        .iter()
        .filter(|e| !e.non_cysteine_anchors.is_empty())
        .map(|e| (e.label.as_str(), e.non_cysteine_anchors.clone()))
        .collect::<Vec<_>>();
    assert_eq!(atypical.len(), 6);
    assert!(
        atypical
            .iter()
            .any(|(s, n)| *s == "IGHV1-38-4*01" && n.contains(&104))
    );
    assert!(
        atypical
            .iter()
            .any(|(s, n)| *s == "IGKV1D-37*01" && n.contains(&104))
    );
    assert!(
        atypical
            .iter()
            .any(|(s, n)| *s == "IGLV3-2*03" && n.contains(&104))
    );
    for e in corpus.v_entries() {
        assert_eq!(e.alignment_labels.first().unwrap(), "1");
        assert_eq!(e.alignment_labels.last().unwrap(), "104");
        assert_eq!(e.aligned_protein.len(), 104);
        assert_eq!(e.cdrs.iter().map(|c| c.kind).collect::<Vec<_>>(), [1, 2]);
    }
    let encoded = serde_json::to_vec(&corpus).unwrap();
    let copy: Corpus = serde_json::from_slice(&encoded).unwrap();
    assert_eq!(copy.revision().unwrap(), corpus.revision().unwrap());
}
#[test]
fn corpus_versions_column_ownership_and_annotation_ledgers_are_current_and_exact() {
    let c = load(&Species::Human).unwrap();
    let original = serde_json::to_value(&c).unwrap();
    for version in [0, 1, 2, 4] {
        let mut v = original.clone();
        v["version"] = version.into();
        assert!(
            serde_json::from_value::<Corpus>(v)
                .unwrap()
                .revision()
                .is_err()
        );
    }
    for (field, value) in [
        ("alignment_labels", serde_json::json!(["23"])),
        ("source", serde_json::json!(3)),
        ("non_cysteine_anchors", serde_json::json!([23, 104])),
        ("aligned_protein", serde_json::json!("C")),
    ] {
        let mut v = original.clone();
        v["v_entries"][0][field] = value;
        assert!(
            serde_json::from_value::<Corpus>(v)
                .unwrap()
                .revision()
                .is_err(),
            "{field}"
        );
    }
    let mut v = original.clone();
    v["sources"][0]["v_records"] = 193.into();
    assert!(
        serde_json::from_value::<Corpus>(v)
            .unwrap()
            .revision()
            .is_err()
    );
    let mut v = original;
    v["v_entries"][0]["cdrs"][0]["start"] = 0.into();
    assert!(
        serde_json::from_value::<Corpus>(v)
            .unwrap()
            .revision()
            .is_err()
    );
}

#[test]
fn j_reference_anchors_complete_codons_and_partial_splice_bases_are_explicit() {
    let corpus = load(&Species::Human).unwrap();
    assert_eq!(corpus.j_entries().len(), 23);
    let examples = [
        ("IGHJ4*02", "WGQGTLVTVSS", 128, "G"),
        ("IGKJ1*01", "FGQGTKVEIK", 127, "C"),
        ("IGLJ2*01", "FGGGTKLTVL", 127, "G"),
    ];
    for (label, protein, last, trailing) in examples {
        let entry = corpus
            .j_entries()
            .iter()
            .find(|e| e.label == label)
            .unwrap();
        assert_eq!(entry.fr4_protein, protein);
        assert_eq!(entry.fr4_positions.first().unwrap(), "118");
        assert_eq!(entry.fr4_positions.last().unwrap(), &last.to_string());
        assert_eq!(entry.trailing_nucleotides, trailing);
        assert_eq!(entry.fr4_nucleotides.len(), 3 * protein.len());
    }
    let original = serde_json::to_value(&corpus).unwrap();
    for (field, value) in [
        ("provider_codon_frame", serde_json::json!(0)),
        ("provider_j_cdr3_end", serde_json::json!(1)),
        ("anchor_offset", serde_json::json!(usize::MAX)),
        ("gene_start", serde_json::json!(0)),
        ("gene_end", serde_json::json!(usize::MAX)),
        ("fr4_positions", serde_json::json!(["118"])),
        ("trailing_nucleotides", serde_json::json!("")),
    ] {
        let mut value_json = original.clone();
        value_json["j_entries"][0][field] = value;
        let parsed: Corpus = serde_json::from_value(value_json).unwrap();
        assert!(parsed.revision().is_err(), "{field}");
    }
}
