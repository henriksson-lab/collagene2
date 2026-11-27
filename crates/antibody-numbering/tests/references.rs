use collagene_antibody_numbering::references::{Corpus, human_v_references};
#[test]
fn pinned_sources_keep_explicit_numbering_and_complete_inclusion_exclusion_ledgers() {
    let corpus = human_v_references().unwrap();
    let sources = corpus.sources();
    assert_eq!(
        sources
            .iter()
            .map(|s| (
                s.locus.as_str(),
                s.release_version,
                s.original_records,
                s.retained_records
            ))
            .collect::<Vec<_>>(),
        vec![
            ("IGH", 10, 240, 194),
            ("IGK", 5, 84, 74),
            ("IGL", 4, 97, 86)
        ]
    );
    assert_eq!(corpus.entries().len(), 354);
    let atypical = corpus
        .entries()
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
    for e in corpus.entries() {
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
    let c = human_v_references().unwrap();
    let original = serde_json::to_value(&c).unwrap();
    for version in [0, 2] {
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
        v["entries"][0][field] = value;
        assert!(
            serde_json::from_value::<Corpus>(v)
                .unwrap()
                .revision()
                .is_err(),
            "{field}"
        );
    }
    let mut v = original.clone();
    v["sources"][0]["retained_records"] = 193.into();
    assert!(
        serde_json::from_value::<Corpus>(v)
            .unwrap()
            .revision()
            .is_err()
    );
    let mut v = original;
    v["entries"][0]["cdrs"][0]["start"] = 0.into();
    assert!(
        serde_json::from_value::<Corpus>(v)
            .unwrap()
            .revision()
            .is_err()
    );
}
