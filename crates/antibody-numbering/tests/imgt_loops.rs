use collagene_antibody_numbering::{
    Insertion,
    imgt_loops::{self, Kind},
    references::human_references,
};
#[test]
fn published_literal_short_and_full_loop_placements_are_not_cdr3_gap_reuse() {
    for (length, expected) in [
        (5, vec![27, 28, 29, 37, 38]),
        (6, vec![27, 28, 29, 36, 37, 38]),
        (7, vec![27, 28, 29, 30, 36, 37, 38]),
        (8, vec![27, 28, 29, 30, 35, 36, 37, 38]),
        (9, vec![27, 28, 29, 30, 31, 35, 36, 37, 38]),
        (10, vec![27, 28, 29, 30, 31, 34, 35, 36, 37, 38]),
        (11, vec![27, 28, 29, 30, 31, 32, 34, 35, 36, 37, 38]),
        (12, (27..=38).collect()),
    ] {
        assert_eq!(
            imgt_loops::positions(Kind::Cdr1, length)
                .unwrap()
                .iter()
                .map(|p| p.number)
                .collect::<Vec<_>>(),
            expected
        );
    }
    for (length, expected) in [
        (0, vec![]),
        (1, vec![56]),
        (2, vec![56, 65]),
        (3, vec![56, 57, 65]),
        (4, vec![56, 57, 64, 65]),
        (5, vec![56, 57, 58, 64, 65]),
        (6, vec![56, 57, 58, 63, 64, 65]),
        (7, vec![56, 57, 58, 59, 63, 64, 65]),
        (8, vec![56, 57, 58, 59, 62, 63, 64, 65]),
        (9, vec![56, 57, 58, 59, 60, 62, 63, 64, 65]),
        (10, (56..=65).collect()),
    ] {
        assert_eq!(
            imgt_loops::positions(Kind::Cdr2, length)
                .unwrap()
                .iter()
                .map(|p| p.number)
                .collect::<Vec<_>>(),
            expected
        );
    }
}
#[test]
fn every_provider_v_reference_agrees_with_the_independent_loop_allocator() {
    let corpus = human_references().unwrap();
    for reference in corpus.v_entries() {
        for (kind, low, high) in [(Kind::Cdr1, 27, 38), (Kind::Cdr2, 56, 65)] {
            let observed = reference
                .alignment_labels
                .iter()
                .zip(reference.aligned_protein.bytes())
                .filter_map(|(label, aa)| {
                    let n = label.parse::<u16>().unwrap();
                    (aa != b'-' && (low..=high).contains(&n)).then_some(n)
                })
                .collect::<Vec<_>>();
            let allocated = imgt_loops::positions(kind, observed.len()).unwrap();
            assert!(allocated.iter().all(|p| p.insertion == Insertion::None));
            assert_eq!(
                allocated.iter().map(|p| p.number).collect::<Vec<_>>(),
                observed,
                "{}",
                reference.label
            );
        }
    }
}
#[test]
fn unsupported_lengths_are_refused_instead_of_inventing_extra_loop_positions() {
    for length in [0, 1, 4, 13, usize::MAX] {
        assert!(imgt_loops::positions(Kind::Cdr1, length).is_err());
    }
    for length in [11, usize::MAX] {
        assert!(imgt_loops::positions(Kind::Cdr2, length).is_err());
    }
    assert_eq!(
        imgt_loops::catalog_revision(),
        imgt_loops::catalog_revision()
    );
}

fn protein(s: &str) -> collagene_model::Molecule {
    use collagene_model::{Alphabet, Molecule, Sequence, Topology};
    Molecule {
        name: "Provided loop boundary fixture".into(),
        sequence: Sequence::new(Alphabet::Protein, s).unwrap(),
        topology: Topology::Linear,
        features: vec![],
        metadata: Default::default(),
        original_genbank: None,
    }
}
#[test]
fn provided_boundaries_reconstruct_exact_source_and_distinguish_empty_cdr2() {
    let source = protein("QACDEFGHIKLMNPQRSTVWY");
    let mut anchors = imgt_loops::Anchors {
        source_revision: collagene_antibody_numbering::source_revision(&source).unwrap(),
        left_index: 1,
        right_index: 10,
        provenance: "Provided framework26/39 hypothesis".into(),
    };
    let numbered = imgt_loops::number_anchored_loop(&source, &Kind::Cdr1, &anchors).unwrap();
    assert_eq!((numbered.start, numbered.end), (2, 10));
    assert_eq!(
        numbered
            .residues
            .iter()
            .map(|r| r.residue)
            .collect::<String>(),
        "CDEFGHIK"
    );
    assert_eq!(
        numbered
            .residues
            .iter()
            .map(|r| r.source_index)
            .collect::<Vec<_>>(),
        (2..10).collect::<Vec<_>>()
    );
    assert_eq!(
        numbered
            .residues
            .iter()
            .map(|r| r.position.number)
            .collect::<Vec<_>>(),
        [27, 28, 29, 30, 35, 36, 37, 38]
    );
    numbered
        .validate_for(&source, &Kind::Cdr1, &anchors)
        .unwrap();
    anchors.right_index = 2;
    anchors.provenance = "Provided adjacent framework55/66 hypothesis".into();
    let empty = imgt_loops::number_anchored_loop(&source, &Kind::Cdr2, &anchors).unwrap();
    assert_eq!((empty.start, empty.end), (2, 2));
    assert!(empty.residues.is_empty());
    assert!(imgt_loops::number_anchored_loop(&source, &Kind::Cdr1, &anchors).is_err());
}
#[test]
fn loop_reconstruction_refuses_stops_invalid_coordinates_stale_sources_and_forged_labels() {
    use collagene_antibody_numbering::{Error, source_revision};
    let source = protein("QACDEFGHIKLMNPQRSTVWY");
    let anchors = imgt_loops::Anchors {
        source_revision: source_revision(&source).unwrap(),
        left_index: 1,
        right_index: 10,
        provenance: "Provided framework boundary hypothesis".into(),
    };
    let numbered = imgt_loops::number_anchored_loop(&source, &Kind::Cdr1, &anchors).unwrap();
    let mut forged = numbered.clone();
    forged.residues[0].position.number += 1;
    assert_eq!(
        forged.validate_for(&source, &Kind::Cdr1, &anchors),
        Err(Error::Binding)
    );
    let mut renamed = source.clone();
    renamed.name = "Changed document".into();
    assert_eq!(
        numbered.validate_for(&renamed, &Kind::Cdr1, &anchors),
        Err(Error::Binding)
    );
    assert_eq!(
        numbered.validate_for(&source, &Kind::Cdr2, &anchors),
        Err(Error::Binding)
    );
    for (left, right) in [(10, 1), (usize::MAX, 10), (1, source.sequence.len())] {
        let mut a = anchors.clone();
        a.left_index = left;
        a.right_index = right;
        assert!(imgt_loops::number_anchored_loop(&source, &Kind::Cdr1, &a).is_err());
    }
    for sequence in ["QA*DEFGHIKLMNPQRSTVWY", "Q*CDEFGHIKLMNPQRSTVWY"] {
        let stopped = protein(sequence);
        let mut a = anchors.clone();
        a.source_revision = source_revision(&stopped).unwrap();
        assert_eq!(
            imgt_loops::number_anchored_loop(&stopped, &Kind::Cdr1, &a),
            Err(Error::Invalid("stop in conditional loop or boundary"))
        );
    }
    let mut a = anchors;
    a.provenance = "\n".into();
    assert!(imgt_loops::number_anchored_loop(&source, &Kind::Cdr1, &a).is_err());
}
