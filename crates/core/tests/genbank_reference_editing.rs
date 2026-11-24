use collagene_core::{formats::*, model::*};

const SOURCE: &str = "LOCUS       reference                 12 bp    DNA     linear   SYN 01-JAN-2026\nDEFINITION  Literal editing fixture.\nACCESSION   AB123456\nVERSION     AB123456.2\nREFERENCE   1  (bases 3 to 8; 10 to 12)\n  AUTHORS   Example,A.\n  TITLE     Literal citation\n  JOURNAL   Unpublished\nFEATURES             Location/Qualifiers\n     misc_feature    3..8\n                     /label=\"literal feature\"\nORIGIN\n        1 aaaaccccgg gg\n//\n";
fn source() -> Molecule {
    import_genbank(SOURCE).unwrap().molecules.remove(0)
}
fn record(text: &str) -> gb_io::seq::Seq {
    gb_io::reader::SeqReader::new(text.as_bytes())
        .next()
        .unwrap()
        .unwrap()
}

#[test]
fn export_replays_source_edits_then_maps_citations_and_clears_active_database_identity() {
    let original = source();
    let edited = original.replace(0, 0, "TT").unwrap();
    let exported = export_genbank(&edited).unwrap();
    let output = record(&exported);
    let expected = record(SOURCE);
    assert_eq!(
        output.references[0].description,
        "1  (bases 5 to 10; 12 to 14)"
    );
    assert_eq!(output.references[0].authors, expected.references[0].authors);
    assert_eq!(output.references[0].title, expected.references[0].title);
    assert_eq!(output.references[0].journal, expected.references[0].journal);
    assert!(
        output.references[0]
            .remark
            .as_deref()
            .unwrap()
            .contains("support for edited DNA is unverified")
    );
    assert!(output.accession.is_none() && output.version.is_none());
    assert!(output.comments.iter().any(|s| s.contains("AB123456.2")));
    assert_eq!(edited.original_genbank.as_deref(), Some(SOURCE));
    let round = import_genbank(&exported).unwrap().molecules.remove(0);
    assert_eq!(round.sequence, edited.sequence);
    assert_eq!(original, source());
}

#[test]
fn same_length_substitution_keeps_numbered_scope_but_marks_ancestral_citation_support() {
    let edited = source().replace(4, 6, "AT").unwrap();
    let output = record(&export_genbank(&edited).unwrap());
    assert_eq!(
        output.references[0].description,
        "1  (bases 3 to 8; 10 to 12)"
    );
    assert!(
        output.references[0]
            .remark
            .as_deref()
            .unwrap()
            .contains("projected coordinate scope")
    );
    assert_eq!(
        output
            .seq
            .iter()
            .map(|b| b.to_ascii_uppercase())
            .collect::<Vec<_>>(),
        b"AAAAATCCGGGG"
    );
}

#[test]
fn boolean_scope_flags_stale_logs_and_changed_original_citations_do_not_bypass_export_guard() {
    let mut bypass = source();
    bypass.sequence = Sequence::dna("TTAAAACCCCGGGG").unwrap();
    bypass
        .metadata
        .insert("genbank:citation_scope".into(), "trusted".into());
    assert!(export_genbank(&bypass).is_err());
    let mut edited = source().replace(0, 0, "TT").unwrap();
    let mut recipe: GenbankEditRecipe =
        serde_json::from_str(&edited.metadata[GENBANK_EDIT_RECIPE_KEY]).unwrap();
    recipe.edits[0].replacement = "AA".into();
    edited.metadata.insert(
        GENBANK_EDIT_RECIPE_KEY.into(),
        serde_json::to_string(&recipe).unwrap(),
    );
    assert!(export_genbank(&edited).is_err());
    let mut edited = source().replace(0, 0, "TT").unwrap();
    edited.original_genbank = Some(SOURCE.replace("Literal citation", "Altered citation"));
    assert!(export_genbank(&edited).is_err());
}

#[test]
fn transformed_edited_sources_archive_recipe_once_without_claiming_current_citation_coordinates() {
    let edited = source().replace(0, 0, "TT").unwrap();
    let recipe = edited.metadata[GENBANK_EDIT_RECIPE_KEY].clone();
    let derived = collagene_core::annotations::reverse_complement(&edited).unwrap();
    assert!(derived.original_genbank.is_none());
    assert!(!derived.metadata.contains_key(GENBANK_EDIT_RECIPE_KEY));
    assert!(!derived.metadata.contains_key("genbank:citation_scope"));
    assert_eq!(
        derived.metadata["annotation_mapping:unverified_genbank_edit_recipe"],
        recipe
    );
    assert_eq!(
        derived.metadata["annotation_mapping:original_genbank"],
        SOURCE
    );
    assert!(export_genbank(&derived).is_ok());
    let again = collagene_core::annotations::reverse_complement(&derived).unwrap();
    assert_eq!(
        again.metadata["annotation_mapping:unverified_genbank_edit_recipe"],
        recipe
    );
    assert_eq!(
        again
            .metadata
            .keys()
            .filter(|key| key.starts_with("annotation_mapping:unverified_genbank_edit_recipe"))
            .count(),
        1
    );
}
