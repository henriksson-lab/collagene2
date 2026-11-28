use collagene_core::{
    oligo_formats::{self as o, ColumnMapping, OligoImportOptions},
    primer_reuse as r,
};
use collagene_model::{Feature, Location, Molecule, Topology};
fn source() -> Molecule {
    Molecule::dna("actual DNA", "AAACCCGGGTTT", Topology::Linear).unwrap()
}
fn feature(location: Location, fields: &[(&str, &str)]) -> Feature {
    Feature {
        kind: "primer_bind".into(),
        location,
        qualifiers: fields
            .iter()
            .map(|(k, v)| (k.to_string(), Some(v.to_string())))
            .collect(),
        needs_review: false,
    }
}
#[test]
fn exact_forward_reverse_and_circular_origin_geometry_supply_only_missing_domains() {
    let mut m = source();
    m.features = vec![
        feature(Location::range(3, 6), &[("label", "forward")]),
        feature(
            Location::Complement(Box::new(Location::range(3, 6))),
            &[("label", "reverse")],
        ),
    ];
    assert_eq!(r::resolve(&m, 0).unwrap().primer.annealing.as_str(), "CCC");
    assert_eq!(r::resolve(&m, 1).unwrap().primer.annealing.as_str(), "GGG");
    let mut circle = Molecule::dna("circle", "ACGTAA", Topology::Circular).unwrap();
    circle.features.push(feature(
        Location::Join(vec![Location::range(4, 6), Location::range(0, 2)]),
        &[("label", "origin primer")],
    ));
    assert_eq!(
        r::resolve(&circle, 0).unwrap().primer.annealing.as_str(),
        "AAAC"
    );
    circle.features[0].location =
        Location::Join(vec![Location::range(0, 1), Location::range(3, 5)]);
    assert!(r::resolve(&circle, 0).is_err());
    circle.features[0].location = Location::range(0, 4);
    circle.features[0].needs_review = true;
    assert!(r::resolve(&circle, 0).is_err());
}
#[test]
fn stored_tails_and_intentional_mismatch_are_preserved_but_conflicts_duplicates_are_not() {
    let mut m = source();
    m.features.push(feature(
        Location::Complement(Box::new(Location::range(3, 6))),
        &[
            ("label", "named mismatch"),
            ("annealing_sequence", "GGA"),
            ("5prime_tail", "TT"),
            ("oligo_sequence", "TTGGA"),
        ],
    ));
    let resolved = r::resolve(&m, 0).unwrap();
    assert_eq!(resolved.primer.name, "named mismatch");
    assert_eq!(resolved.primer.ordered_sequence(), "TTGGA");
    assert!(resolved.warnings.iter().any(|w| w.contains("differs")));
    m.features[0]
        .qualifiers
        .push(("OLIGO_SEQUENCE".into(), Some("TTGGA".into())));
    assert!(r::resolve(&m, 0).is_err());
    m.features[0].qualifiers.pop();
    m.features[0]
        .qualifiers
        .iter_mut()
        .find(|(k, _)| k == "oligo_sequence")
        .unwrap()
        .1 = Some("CCGGA".into());
    assert!(r::resolve(&m, 0).is_err());
    m.features[0] = feature(
        Location::Unresolved("historical site".into()),
        &[
            ("label", "ambiguous stored oligo"),
            ("annealing_sequence", "ACGN"),
            ("5prime_tail", "RY"),
            ("oligo_sequence", "RYACGN"),
        ],
    );
    let resolved = r::resolve(&m, 0).unwrap();
    assert_eq!(resolved.primer.ordered_sequence(), "RYACGN");
    assert!(resolved.warnings.iter().any(|w| w.contains("IUPAC")));
}
#[test]
fn actual_oligo_import_feature_resolves_without_guessing_domains_from_bare_dna() {
    let options = OligoImportOptions {
        columns: ColumnMapping {
            name: 0,
            ordered_sequence: Some(1),
            annealing_sequence: Some(2),
            tail: Some(3),
            modification: None,
            amount: None,
            unit: None,
        },
        ..Default::default()
    };
    let table = o::import_oligos(
        "Name,Ordered,Annealing,Tail\nSaved oligo,TTACGN,ACGN,TT\n",
        &options,
    )
    .unwrap();
    let imported = table.selected_molecules(&[0]).unwrap();
    let m = &imported.molecules[0];
    let resolved = r::resolve_at_revision(m, 0, &m.revision()).unwrap();
    assert_eq!(resolved.primer.name, "Saved oligo");
    assert_eq!(resolved.primer.annealing.as_str(), "ACGN");
    assert_eq!(resolved.primer.tail, "TT");
    assert_eq!(resolved.primer.ordered_sequence(), "TTACGN");
    let bare = Molecule::dna("bare oligo", "TTACGN", Topology::Linear).unwrap();
    assert!(r::resolve(&bare, 0).is_err());
    assert_eq!(r::list(&bare, 100).unwrap().total, 0);
}
#[test]
fn bounded_lists_keep_actual_indices_and_stale_document_pins_cannot_resolve() {
    let mut m = source();
    m.features.push(Feature {
        kind: "gene".into(),
        location: Location::range(0, 2),
        qualifiers: vec![],
        needs_review: false,
    });
    m.features.push(feature(
        Location::range(0, 3),
        &[("label", "first actual primer")],
    ));
    m.features.push(feature(
        Location::Unresolved("unknown".into()),
        &[("label", "unusable")],
    ));
    m.features
        .push(feature(Location::range(6, 9), &[("label", "third")]));
    let list = r::list(&m, 2).unwrap();
    assert_eq!((list.total, list.omitted), (3, 1));
    assert_eq!(list.primers[0].feature_index, 1);
    assert!(list.primers[0].reusable);
    assert!(!list.primers[1].reusable);
    assert!(list.primers[1].diagnostic.is_some());
    let pin = m.revision();
    assert_eq!(list.source_revision, pin);
    m.name = "changed document revision".into();
    assert!(r::resolve_at_revision(&m, 1, &pin).is_err());
    m.features[1]
        .qualifiers
        .push(("source_revision".into(), Some("0".repeat(64))));
    assert!(
        r::resolve_at_revision(&m, 1, &m.revision())
            .unwrap()
            .warnings
            .iter()
            .any(|w| w.contains("historical"))
    );
    assert!(r::list(&m, 1001).is_err());
    m.features[1]
        .qualifiers
        .push(("label".into(), Some("conflicting name".into())));
    assert!(r::resolve(&m, 1).is_err());
}
