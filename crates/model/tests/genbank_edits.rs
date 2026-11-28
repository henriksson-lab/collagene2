use collagene_model::*;

const SOURCE: &str = "LOCUS       reference                 12 bp    DNA     linear   SYN 01-JAN-2026\nREFERENCE   1  (bases 3 to 8; 10 to 12)\n  AUTHORS   Example,A.\n  TITLE     Literal citation\n  JOURNAL   Unpublished\nORIGIN\n        1 aaaaccccgg gg\n//\n";
fn source() -> Molecule {
    let mut m = Molecule::dna("reference", "AAAACCCCGGGG", Topology::Linear).unwrap();
    m.original_genbank = Some(SOURCE.into());
    m
}
fn recipe(m: &Molecule) -> GenbankEditRecipe {
    verified_genbank_edit_recipe(m, &source().sequence, 0)
        .unwrap()
        .unwrap()
}

#[test]
fn source_stays_exact_and_literal_joined_citation_spans_shift_and_clip() {
    let original = source();
    let inserted = original.replace(0, 0, "TT").unwrap();
    assert_eq!(inserted.original_genbank.as_deref(), Some(SOURCE));
    assert_eq!(
        remap_genbank_reference("1  (bases 3 to 8; 10 to 12)", &recipe(&inserted)).unwrap(),
        "1  (bases 5 to 10; 12 to 14)"
    );
    let deleted = original.replace(4, 6, "").unwrap();
    assert_eq!(
        remap_genbank_reference("1  (bases 3 to 8; 10 to 12)", &recipe(&deleted)).unwrap(),
        "1  (bases 3 to 6; 8 to 10)"
    );
    let substituted = original.replace(4, 6, "AT").unwrap();
    assert_eq!(
        remap_genbank_reference("1  (bases 3 to 8; 10 to 12)", &recipe(&substituted)).unwrap(),
        "1  (bases 3 to 8; 10 to 12)"
    );
    assert!(substituted.metadata.contains_key("genbank:citation_scope"));
    assert_eq!(original, source());
}

#[test]
fn consecutive_insert_typing_and_backspace_coalesce_without_losing_citation_boundaries() {
    let typed = source()
        .replace(2, 2, "A")
        .unwrap()
        .replace(3, 3, "T")
        .unwrap();
    let r = recipe(&typed);
    assert_eq!(r.edits.len(), 1);
    assert_eq!(r.edits[0].replacement, "AT");
    assert_eq!(
        remap_genbank_reference("1  (bases 3 to 8; 10 to 12)", &r).unwrap(),
        "1  (bases 5 to 10; 12 to 14)"
    );
    let erased = typed.replace(2, 4, "").unwrap();
    assert!(recipe(&erased).edits.is_empty());
    assert_eq!(erased.sequence, source().sequence);
    // A replacement followed by boundary insertion must stay separate: folding
    // that insertion into the replacement would silently extend citation scope.
    let replaced = source()
        .replace(2, 8, "TT")
        .unwrap()
        .replace(4, 4, "A")
        .unwrap();
    assert_eq!(recipe(&replaced).edits.len(), 2);
    assert_eq!(
        remap_genbank_reference("1  (bases 3 to 8; 10 to 12)", &recipe(&replaced)).unwrap(),
        "1  (bases 3 to 4; 7 to 9)"
    );
}

#[test]
fn unsupported_scope_deleted_coverage_and_coordinate_qualifiers_remain_explicit() {
    let mut narrow = source();
    narrow.original_genbank = Some(SOURCE.replace("3 to 8; 10 to 12", "3 to 4"));
    assert!(narrow.replace(2, 4, "").is_err());
    let mut unusual = source();
    unusual.original_genbank = Some(SOURCE.replace("1  (bases 3 to 8; 10 to 12)", "1  (sites)"));
    assert!(unusual.replace(0, 0, "A").is_err());
    let mut guarded = source();
    guarded.features.push(Feature {
        kind: "CDS".into(),
        location: Location::range(0, 3),
        qualifiers: vec![("transl_except".into(), Some("(pos:1..3,aa:Sec)".into()))],
        needs_review: false,
    });
    let before = guarded.clone();
    let outside = guarded.replace(5, 5, "A").unwrap();
    assert_eq!(outside.sequence.as_str(), "AAAACACCCGGGG");
    assert_eq!(outside.features, guarded.features);
    assert_eq!(outside.original_genbank.as_deref(), Some(SOURCE));
    assert_eq!(
        remap_genbank_reference("1  (bases 3 to 8; 10 to 12)", &recipe(&outside)).unwrap(),
        "1  (bases 3 to 9; 11 to 13)"
    );
    let shifted = guarded.replace(0, 0, "TT").unwrap();
    assert_eq!(shifted.features[0].location, Location::range(2, 5));
    assert!(!shifted.features[0].needs_review);
    assert_eq!(
        shifted.features[0].qualifiers,
        vec![("transl_except".into(), Some("(pos:3..5,aa:Sec)".into()))]
    );
    assert_eq!(shifted.original_genbank.as_deref(), Some(SOURCE));
    assert_eq!(
        remap_genbank_reference("1  (bases 3 to 8; 10 to 12)", &recipe(&shifted)).unwrap(),
        "1  (bases 5 to 10; 12 to 14)"
    );
    assert_eq!(guarded, before);
    let mut unsupported = guarded.clone();
    unsupported.features[0]
        .qualifiers
        .push(("rpt_unit_range".into(), Some("1..3".into())));
    let original = unsupported.clone();
    assert!(unsupported.replace(5, 5, "A").is_err());
    assert_eq!(unsupported, original);
}

#[test]
fn stale_raw_and_sequence_bindings_and_oversized_recipes_do_not_authorize_replay() {
    let edited = source().replace(0, 0, "TT").unwrap();
    let mut changed = edited.clone();
    changed.original_genbank = Some(SOURCE.replace("Literal citation", "Changed citation"));
    assert!(verified_genbank_edit_recipe(&changed, &source().sequence, 0).is_err());
    let mut forged = edited.clone();
    let mut r: GenbankEditRecipe =
        serde_json::from_str(&forged.metadata[GENBANK_EDIT_RECIPE_KEY]).unwrap();
    r.edits[0].replacement = "AA".into();
    forged.metadata.insert(
        GENBANK_EDIT_RECIPE_KEY.into(),
        serde_json::to_string(&r).unwrap(),
    );
    assert!(verified_genbank_edit_recipe(&forged, &source().sequence, 0).is_err());
    forged.metadata.insert(
        GENBANK_EDIT_RECIPE_KEY.into(),
        " ".repeat(MAX_GENBANK_EDIT_RECIPE_BYTES + 1),
    );
    assert!(verified_genbank_edit_recipe(&forged, &source().sequence, 0).is_err());
}
