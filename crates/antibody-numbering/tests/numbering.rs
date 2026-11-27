use collagene_antibody_numbering::*;
use collagene_model::{Alphabet, Molecule, Sequence, Topology};
fn labels(n: usize) -> Vec<String> {
    imgt_cdr3_positions(n)
        .unwrap()
        .iter()
        .map(ToString::to_string)
        .collect()
}
#[test]
fn literal_published_cdr3_lengths_preserve_gaps_and_reverse_right_insertions() {
    for (n, expected) in [
        (2, "105 117"),
        (6, "105 106 107 115 116 117"),
        (13, "105 106 107 108 109 110 111 112 113 114 115 116 117"),
        (
            14,
            "105 106 107 108 109 110 111 112.1 112 113 114 115 116 117",
        ),
        (
            15,
            "105 106 107 108 109 110 111 111.1 112.1 112 113 114 115 116 117",
        ),
        (
            17,
            "105 106 107 108 109 110 111 111.1 111.2 112.2 112.1 112 113 114 115 116 117",
        ),
        (
            21,
            "105 106 107 108 109 110 111 111.1 111.2 111.3 111.4 112.4 112.3 112.2 112.1 112 113 114 115 116 117",
        ),
    ] {
        assert_eq!(
            labels(n),
            expected.split_whitespace().collect::<Vec<_>>(),
            "length {n}"
        );
    }
    for n in 2..=MAX_CDR3_RESIDUES {
        let p = imgt_cdr3_positions(n).unwrap();
        assert_eq!(p.len(), n);
        let unique = p
            .iter()
            .map(ToString::to_string)
            .collect::<std::collections::BTreeSet<_>>();
        assert_eq!(unique.len(), n);
        assert_eq!(p.first().unwrap().to_string(), "105");
        assert_eq!(p.last().unwrap().to_string(), "117");
    }
}
#[test]
fn insertion_labels_are_exact_syntax_not_floats_or_automatic_scheme_conversions() {
    for text in ["1", "52A", "100AA", "111.1", "111.10", "112.2", "112.10"] {
        let p: Position = text.parse().unwrap();
        assert_eq!(p.to_string(), text);
    }
    assert_ne!(
        "111.1".parse::<Position>().unwrap(),
        "111.10".parse::<Position>().unwrap()
    );
    for text in [
        "0", "01", "52a", " 52A", "52.A", "111.0", "111.01", "+52", "52.1e2", "100ABCDE", "10000",
    ] {
        assert!(text.parse::<Position>().is_err(), "{text}");
    }
    for n in [0, 1, MAX_CDR3_RESIDUES + 1, usize::MAX] {
        assert!(imgt_cdr3_positions(n).is_err());
    }
}
fn source(s: &str) -> Molecule {
    Molecule {
        name: "Authored anchor context".into(),
        sequence: Sequence::new(Alphabet::Protein, s).unwrap(),
        topology: Topology::Linear,
        features: vec![],
        metadata: Default::default(),
        original_genbank: None,
    }
}
#[test]
fn supplied_anchor_allocation_binds_actual_source_without_claiming_domain_detection() {
    let m = source("ACNNNNNNWGQG*");
    let a = Anchors {
        source_revision: m.revision(),
        cysteine_index: 1,
        j_anchor_index: 8,
        provenance:
            "Manually provided positions for a conditional allocator fixture; not antibody identity"
                .into(),
    };
    let r = number_anchored_imgt_cdr3(&m, &a).unwrap();
    assert_eq!(r.source_revision, m.revision());
    assert_eq!((r.start, r.end), (2, 8));
    assert_eq!(r.residues[0].source_index, 2);
    assert_eq!(r.residues[0].residue, 'N');
    assert_eq!(r.residues[5].position.to_string(), "117");
    r.validate_for(&m, &a).unwrap();
    let decoded: ConditionalCdr3 =
        serde_json::from_slice(&serde_json::to_vec(&r).unwrap()).unwrap();
    decoded.validate_for(&m, &a).unwrap();
    let mut forged = r.clone();
    forged.residues[2].position = "111.10".parse().unwrap();
    assert!(forged.validate_for(&m, &a).is_err());
    let mut changed = m.clone();
    changed.name = "Different source metadata".into();
    assert!(r.validate_for(&changed, &a).is_err());
    let mut a = a.clone();
    a.provenance = "Changed assumption".into();
    assert!(r.validate_for(&m, &a).is_err());
}
#[test]
fn absent_ambiguous_or_out_of_bounds_anchors_and_interrupted_loops_refuse() {
    for text in [
        "AXNNNNNNWGQG",
        "ACNNNNNNXGQG",
        "ACNNNNNNWAQG",
        "ACNNNNNNWGQX",
        "ACNN*NNNWGQG",
    ] {
        let m = source(text);
        let a = Anchors {
            source_revision: m.revision(),
            cysteine_index: 1,
            j_anchor_index: 8,
            provenance: "Negative fixture".into(),
        };
        assert!(number_anchored_imgt_cdr3(&m, &a).is_err());
    }
    let m = source("ACNNNNNNFGQG");
    let a = Anchors {
        source_revision: m.revision(),
        cysteine_index: 1,
        j_anchor_index: 8,
        provenance: "Positive Phe anchor".into(),
    };
    number_anchored_imgt_cdr3(&m, &a).unwrap();
    for (c, j) in [(usize::MAX, 8), (1, usize::MAX), (8, 1), (1, 11)] {
        let mut b = a.clone();
        b.cysteine_index = c;
        b.j_anchor_index = j;
        assert!(number_anchored_imgt_cdr3(&m, &b).is_err());
    }
    let dna = Molecule::dna("DNA", "ACGT", Topology::Linear).unwrap();
    assert!(source_revision(&dna).is_err());
}
