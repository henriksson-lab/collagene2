use collagene_model::{Error, Feature, Location, Molecule, Topology};
use collagene_recombination::{
    Limits, feature_preflight, manufacturer_context_profile, model_preflight, source_preflight,
};

fn molecule() -> Molecule {
    Molecule::dna("borrowed", "ACGT", Topology::Linear).unwrap()
}

#[test]
fn borrowed_guard_rejects_deep_native_trees_before_recursive_validation() {
    let mut source = molecule();
    let mut location = Location::Range {
        start: 0,
        end: 1,
        before: false,
        after: false,
    };
    for _ in 0..128 {
        location = Location::Complement(Box::new(location));
    }
    source.features.push(Feature {
        kind: "misc_feature".into(),
        location,
        qualifiers: vec![],
        needs_review: false,
    });
    assert!(matches!(model_preflight(&[&source]), Err(Error::TooLarge)));
    assert!(matches!(
        feature_preflight(&[&source.features]),
        Err(Error::TooLarge)
    ));
    assert!(matches!(
        source_preflight(
            [&source, &molecule()],
            &manufacturer_context_profile().unwrap(),
            &Limits::default(),
            None,
        ),
        Err(Error::TooLarge)
    ));
}

#[test]
fn source_and_target_share_the_borrowed_text_budget() {
    let mut source = molecule();
    let mut target = molecule();
    for i in 0..5 {
        source.metadata.insert(i.to_string(), "x".repeat(900_000));
        target.metadata.insert(i.to_string(), "y".repeat(900_000));
    }
    model_preflight(&[&source]).unwrap();
    model_preflight(&[&target]).unwrap();
    assert!(matches!(
        source_preflight(
            [&source, &molecule()],
            &manufacturer_context_profile().unwrap(),
            &Limits::default(),
            Some(&target),
        ),
        Err(Error::TooLarge)
    ));
}

#[test]
fn borrowed_guard_bounds_model_count_and_source_bases() {
    let source = molecule();
    assert!(matches!(
        model_preflight(&vec![&source; 1025]),
        Err(Error::TooLarge)
    ));
    let limits = Limits {
        max_input_bases: 7,
        ..Limits::default()
    };
    assert!(matches!(
        source_preflight(
            [&source, &source],
            &manufacturer_context_profile().unwrap(),
            &limits,
            None,
        ),
        Err(Error::TooLarge)
    ));
}
