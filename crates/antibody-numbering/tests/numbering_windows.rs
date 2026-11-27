use collagene_antibody_numbering::{
    Error, numbering_windows::*, source_revision, v_mapping::Locus,
};
use collagene_model::{Alphabet, Molecule, Sequence, Topology};
use std::collections::BTreeSet;
fn text(layout: &Layout) -> String {
    layout
        .residues
        .iter()
        .map(ToString::to_string)
        .collect::<Vec<_>>()
        .join(" ")
}
fn source(sequence: &str) -> Molecule {
    Molecule {
        name: "Authored conditional window".into(),
        sequence: Sequence::new(Alphabet::Protein, sequence).unwrap(),
        topology: Topology::Linear,
        features: vec![],
        metadata: Default::default(),
        original_genbank: None,
    }
}
fn fixture(window: &Window, locus: &Locus, length: usize) -> (Molecule, Anchors) {
    let left = if *window == Window::MiddleLoop {
        'Q'
    } else {
        'C'
    };
    let tail = match window {
        Window::CysToTrp => "WAA",
        Window::MiddleLoop => "RAA",
        Window::CysToJ if *locus == Locus::Heavy => "WGQG",
        Window::CysToJ => "FGQG",
    };
    let m = source(&format!("AA{left}{}{tail}", "A".repeat(length)));
    let a = Anchors {
        source_revision: source_revision(&m).unwrap(),
        left_index: 2,
        right_index: length + 3,
        provenance: "Authored boundary interpretation; no biological numbering claim".into(),
    };
    (m, a)
}
#[test]
fn published_heavy_first_window_and_middle_loop_use_distinct_insertion_rules() {
    let k = layout(&Convention::Kabat1991, &Locus::Heavy, &Window::CysToTrp, 15).unwrap();
    assert_eq!(
        (k.left.to_string(), k.right.to_string()),
        ("22".into(), "36".into())
    );
    assert_eq!(text(&k), "23 24 25 26 27 28 29 30 31 32 33 34 35 35A 35B");
    let c = layout(
        &Convention::Chothia1997,
        &Locus::Heavy,
        &Window::CysToTrp,
        15,
    )
    .unwrap();
    assert_eq!(text(&c), "23 24 25 26 27 28 29 30 31 31A 31B 32 33 34 35");
    assert_eq!(
        text(
            &layout(
                &Convention::Kabat1991,
                &Locus::Heavy,
                &Window::MiddleLoop,
                6
            )
            .unwrap()
        ),
        "51 52 53 55 56 57"
    );
    assert_eq!(
        text(
            &layout(
                &Convention::Chothia1997,
                &Locus::Heavy,
                &Window::MiddleLoop,
                8
            )
            .unwrap()
        ),
        "51 52 52A 53 54 55 56 57"
    );
    assert!(
        layout(
            &Convention::Chothia1997,
            &Locus::Heavy,
            &Window::CysToTrp,
            16
        )
        .is_err()
    );
}
#[test]
fn kappa_and_lambda_chothia_patterns_are_not_a_shared_generic_insertion_counter() {
    let k = layout(
        &Convention::Chothia1997,
        &Locus::Kappa,
        &Window::CysToTrp,
        16,
    )
    .unwrap();
    assert_eq!(
        text(&k),
        "24 25 26 27 28 29 30 30A 30B 30C 30E 30F 31 32 33 34"
    );
    let l = layout(
        &Convention::Chothia1997,
        &Locus::Lambda,
        &Window::CysToTrp,
        12,
    )
    .unwrap();
    assert_eq!(text(&l), "24 25 26 27 30 30A 30B 30C 31 32 33 34");
    assert_eq!(
        text(
            &layout(
                &Convention::Chothia1997,
                &Locus::Kappa,
                &Window::CysToTrp,
                12
            )
            .unwrap()
        ),
        "24 25 26 27 28 29 30 30A 31 32 33 34"
    );
    assert_eq!(
        text(&layout(&Convention::Kabat1991, &Locus::Kappa, &Window::CysToTrp, 10).unwrap()),
        "24 25 26 27 29 30 31 32 33 34"
    );
    assert_eq!(
        text(
            &layout(
                &Convention::Chothia1997,
                &Locus::Kappa,
                &Window::CysToTrp,
                10
            )
            .unwrap()
        ),
        "24 25 26 27 28 29 30 32 33 34"
    );
    let hil = layout(
        &Convention::Chothia1997,
        &Locus::Lambda,
        &Window::CysToTrp,
        11,
    )
    .unwrap();
    assert_eq!(text(&hil), "24 25 26 27 30 30A 30B 31 32 33 34");
    assert!(
        layout(
            &Convention::Chothia1997,
            &Locus::Lambda,
            &Window::CysToTrp,
            10
        )
        .is_err()
    );
}
#[test]
fn junction_window_labels_are_not_cdr_boundaries_or_imgt_central_gap_reuse() {
    let h = layout(&Convention::Kabat1991, &Locus::Heavy, &Window::CysToJ, 13).unwrap();
    assert_eq!(
        (h.left.to_string(), h.right.to_string()),
        ("92".into(), "103".into())
    );
    assert_eq!(text(&h), "93 94 95 96 97 98 99 100 100A 100B 100C 101 102");
    assert_eq!(
        text(&layout(&Convention::Kabat1991, &Locus::Heavy, &Window::CysToJ, 9).unwrap()),
        "93 94 95 96 97 99 100 101 102"
    );
    let l = layout(
        &Convention::Chothia1997,
        &Locus::Lambda,
        &Window::CysToJ,
        11,
    )
    .unwrap();
    assert_eq!(text(&l), "89 90 91 92 93 94 95 95A 95B 96 97");
    assert_eq!(
        text(&layout(&Convention::Kabat1991, &Locus::Kappa, &Window::CysToJ, 7).unwrap()),
        "89 90 91 92 95 96 97"
    );
    assert!(layout(&Convention::Chothia1997, &Locus::Kappa, &Window::CysToJ, 11).is_err());
    assert!(layout(&Convention::Kabat1991, &Locus::Heavy, &Window::CysToJ, 8).is_err());
    assert!(serde_json::from_str::<Convention>("\"Chothia1987\"").is_err());
}
#[test]
fn supported_profile_sweep_has_exact_lengths_unique_labels_and_bound_source_order() {
    let mut supported = 0;
    for convention in [Convention::Kabat1991, Convention::Chothia1997] {
        for locus in [Locus::Heavy, Locus::Kappa, Locus::Lambda] {
            for window in [Window::CysToTrp, Window::MiddleLoop, Window::CysToJ] {
                for length in 0..=25 {
                    let Ok(labels) = layout(&convention, &locus, &window, length) else {
                        continue;
                    };
                    supported += 1;
                    assert_eq!(labels.residues.len(), length);
                    assert_eq!(
                        labels
                            .residues
                            .iter()
                            .map(ToString::to_string)
                            .collect::<BTreeSet<_>>()
                            .len(),
                        length
                    );
                    assert!(
                        !labels.residues.contains(&labels.left)
                            && !labels.residues.contains(&labels.right)
                    );
                    let (m, a) = fixture(&window, &locus, length);
                    let r = number_anchored_window(&m, &convention, &locus, &window, &a).unwrap();
                    assert_eq!((r.start, r.end), (3, 3 + length));
                    assert_eq!(
                        r.residues
                            .iter()
                            .map(|r| r.position.clone())
                            .collect::<Vec<_>>(),
                        labels.residues
                    );
                    assert!(
                        r.residues
                            .iter()
                            .enumerate()
                            .all(|(i, r)| r.source_index == 3 + i && r.residue == 'A')
                    );
                    r.validate_for(&m, &convention, &locus, &window, &a)
                        .unwrap();
                }
            }
        }
    }
    assert_eq!(supported, 80);
}
#[test]
fn source_binding_anchors_stops_provenance_and_supplied_context_are_reconstructed() {
    let convention = Convention::Kabat1991;
    let locus = Locus::Heavy;
    let window = Window::CysToJ;
    let (m, a) = fixture(&window, &locus, 13);
    let r = number_anchored_window(&m, &convention, &locus, &window, &a).unwrap();
    let mut forged = r.clone();
    forged.residues[0].position = "95".parse().unwrap();
    assert!(
        forged
            .validate_for(&m, &convention, &locus, &window, &a)
            .is_err()
    );
    let mut changed = m.clone();
    changed.name.push('x');
    assert!(
        r.validate_for(&changed, &convention, &locus, &window, &a)
            .is_err()
    );
    assert!(
        r.validate_for(&m, &Convention::Chothia1997, &locus, &window, &a)
            .is_err()
    );
    assert!(
        r.validate_for(&m, &convention, &Locus::Kappa, &window, &a)
            .is_err()
    );
    for (left, right) in [(usize::MAX, 3), (5, 4), (2, m.sequence.len()), (2, 2)] {
        let mut bad = a.clone();
        bad.left_index = left;
        bad.right_index = right;
        assert!(number_anchored_window(&m, &convention, &locus, &window, &bad).is_err());
    }
    for provenance in [
        "".to_string(),
        " \t".into(),
        "x\ny".into(),
        "x".repeat(1025),
    ] {
        let mut bad = a.clone();
        bad.provenance = provenance;
        assert!(number_anchored_window(&m, &convention, &locus, &window, &bad).is_err());
    }
    for index in [
        a.left_index,
        3,
        a.right_index,
        a.right_index + 1,
        a.right_index + 2,
        a.right_index + 3,
    ] {
        let mut s = m.sequence.as_str().as_bytes().to_vec();
        s[index] = b'*';
        let bad = source(std::str::from_utf8(&s).unwrap());
        let mut anchors = a.clone();
        anchors.source_revision = source_revision(&bad).unwrap();
        assert!(number_anchored_window(&bad, &convention, &locus, &window, &anchors).is_err());
    }
    let mut s = m.sequence.as_str().as_bytes().to_vec();
    s[a.right_index] = b'F';
    let bad = source(std::str::from_utf8(&s).unwrap());
    let mut anchors = a;
    anchors.source_revision = source_revision(&bad).unwrap();
    assert!(number_anchored_window(&bad, &convention, &locus, &window, &anchors).is_err());
    assert!(matches!(
        layout(&convention, &locus, &window, usize::MAX),
        Err(Error::Limit(_))
    ));
}
