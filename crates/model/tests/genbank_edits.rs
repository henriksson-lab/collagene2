use collagene_model::genbank_edit::{CitationProjection, EditLog, KEY, MAX_BYTES, VerifiedLog};
use collagene_model::*;

const SOURCE: &str = "LOCUS       reference                 12 bp    DNA     linear   SYN 01-JAN-2026\nREFERENCE   1  (bases 3 to 8; 10 to 12)\n  AUTHORS   Example,A.\n  TITLE     Literal citation\n  JOURNAL   Unpublished\nORIGIN\n        1 aaaaccccgg gg\n//\n";
fn source() -> Molecule {
    let mut m = Molecule::dna("reference", "AAAACCCCGGGG", Topology::Linear).unwrap();
    m.original_genbank = Some(SOURCE.into());
    m
}
fn recipe(m: &Molecule) -> VerifiedLog {
    VerifiedLog::verify(m, &source().sequence, 0)
        .unwrap()
        .unwrap()
}

fn mapped_description(description: &str, log: &VerifiedLog) -> String {
    match log.project_references(&[description]).unwrap().remove(0) {
        CitationProjection::Mapped { description, .. } => description,
        other => panic!("expected mapped citation, got {other:?}"),
    }
}

#[test]
fn source_stays_exact_and_literal_joined_citation_spans_shift_and_clip() {
    let original = source();
    let inserted = original.replace(0, 0, "TT").unwrap();
    assert_eq!(inserted.original_genbank.as_deref(), Some(SOURCE));
    assert_eq!(
        mapped_description("1  (bases 3 to 8; 10 to 12)", &recipe(&inserted)),
        "1  (bases 5 to 10; 12 to 14)"
    );
    let deleted = original.replace(4, 6, "").unwrap();
    assert_eq!(
        mapped_description("1  (bases 3 to 8; 10 to 12)", &recipe(&deleted)),
        "1  (bases 3 to 6; 8 to 10)"
    );
    let substituted = original.replace(4, 6, "AT").unwrap();
    assert_eq!(
        mapped_description("1  (bases 3 to 8; 10 to 12)", &recipe(&substituted)),
        "1  (bases 3 to 4; 7 to 8; 10 to 12)"
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
    assert_eq!(r.edits().len(), 1);
    assert_eq!(r.edits()[0].replacement, "AT");
    assert_eq!(
        mapped_description("1  (bases 3 to 8; 10 to 12)", &r),
        "1  (bases 5 to 10; 12 to 14)"
    );
    let erased = typed.replace(2, 4, "").unwrap();
    assert!(recipe(&erased).edits().is_empty());
    assert_eq!(erased.sequence, source().sequence);
    // A replacement and a later boundary insertion remain separate edit steps.
    // The replaced first cited interval contributes no surviving source bases.
    let replaced = source()
        .replace(2, 8, "TT")
        .unwrap()
        .replace(4, 4, "A")
        .unwrap();
    assert_eq!(recipe(&replaced).edits().len(), 2);
    assert_eq!(
        mapped_description("1  (bases 3 to 8; 10 to 12)", &recipe(&replaced)),
        "1  (bases 7 to 9)"
    );
}

#[test]
fn unsupported_scope_deleted_coverage_and_coordinate_qualifiers_remain_explicit() {
    let mut narrow = source();
    narrow.original_genbank = Some(SOURCE.replace("3 to 8; 10 to 12", "3 to 4"));
    assert!(narrow.replace(2, 4, "").is_ok());
    let mut unusual = source();
    unusual.original_genbank = Some(SOURCE.replace("1  (bases 3 to 8; 10 to 12)", "1  (sites)"));
    assert!(unusual.replace(0, 0, "A").is_ok());
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
        mapped_description("1  (bases 3 to 8; 10 to 12)", &recipe(&outside)),
        "1  (bases 3 to 5; 7 to 9; 11 to 13)"
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
        mapped_description("1  (bases 3 to 8; 10 to 12)", &recipe(&shifted)),
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
    assert!(VerifiedLog::verify(&changed, &source().sequence, 0).is_err());
    let mut forged = edited.clone();
    let mut r: EditLog = serde_json::from_str(&forged.metadata[KEY]).unwrap();
    r.edits[0].replacement = "AA".into();
    forged
        .metadata
        .insert(KEY.into(), serde_json::to_string(&r).unwrap());
    assert!(VerifiedLog::verify(&forged, &source().sequence, 0).is_err());
    forged
        .metadata
        .insert(KEY.into(), " ".repeat(MAX_BYTES + 1));
    assert!(VerifiedLog::verify(&forged, &source().sequence, 0).is_err());
}
