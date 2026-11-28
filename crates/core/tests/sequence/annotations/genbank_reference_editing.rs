use collagene_core::formats::sequence::*;
use collagene_model::genbank_edit::{EditLog, KEY};
use collagene_model::*;

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
fn same_length_substitution_excludes_replacement_bases_and_marks_ancestral_citation_support() {
    let edited = source().replace(4, 6, "AT").unwrap();
    let output = record(&export_genbank(&edited).unwrap());
    assert_eq!(
        output.references[0].description,
        "1  (bases 3 to 4; 7 to 8; 10 to 12)"
    );
    // GenBank REMARK text is wrapped at column79; gb-io preserves those line
    // breaks when parsing. Test citation semantics independently of wrapping.
    let remark = output.references[0]
        .remark
        .as_deref()
        .unwrap()
        .split_whitespace()
        .collect::<Vec<_>>()
        .join(" ");
    assert!(remark.contains("projected coordinate scope"), "{remark}");
    assert!(
        remark.contains("citation support for edited DNA is unverified"),
        "{remark}"
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
    let mut recipe: EditLog = serde_json::from_str(&edited.metadata[KEY]).unwrap();
    recipe.edits[0].replacement = "AA".into();
    edited
        .metadata
        .insert(KEY.into(), serde_json::to_string(&recipe).unwrap());
    assert!(export_genbank(&edited).is_err());
    let mut edited = source().replace(0, 0, "TT").unwrap();
    edited.original_genbank = Some(SOURCE.replace("Literal citation", "Altered citation"));
    assert!(export_genbank(&edited).is_err());
}

#[test]
fn transformed_edited_sources_archive_recipe_once_without_claiming_current_citation_coordinates() {
    let edited = source().replace(0, 0, "TT").unwrap();
    let recipe = edited.metadata[KEY].clone();
    let derived = collagene_core::sequence::annotations::reverse_complement(&edited).unwrap();
    assert!(derived.original_genbank.is_none());
    assert!(!derived.metadata.contains_key(KEY));
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
    let again = collagene_core::sequence::annotations::reverse_complement(&derived).unwrap();
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

fn normalized(text: &str) -> String {
    text.split_whitespace().collect::<Vec<_>>().join(" ")
}

#[test]
fn legal_noncoordinate_sources_survive_native_reopen_without_claiming_destination_coverage() {
    for description in [
        "1  (sites)",
        "4294967296  (sites)",
        "1",
        "1  (annotated by the submitting laboratory)",
        "1  (bases 3 to 8) trailing material",
        "1  (bases 3 to 8; bases 10 to 12) trailing",
    ] {
        let raw = SOURCE.replace("1  (bases 3 to 8; 10 to 12)", description)
            .replace("  JOURNAL   Unpublished", "  JOURNAL   Unpublished\n   PUBMED   12345678\n  REMARK    Original laboratory context.");
        let original = import_genbank(&raw).unwrap().molecules.remove(0);
        let edited = original.replace(0, 0, "TT").unwrap();
        assert_eq!(edited.original_genbank.as_deref(), Some(raw.as_str()));
        let native = collagene_core::protocol::save_document(&edited).unwrap();
        let reopened = collagene_core::protocol::load_document(&native).unwrap();
        assert_eq!(reopened, edited);
        let exported = export_genbank(&reopened).unwrap();
        assert_eq!(export_genbank(&reopened).unwrap(), exported);
        let output = record(&exported);
        assert_eq!(
            output.references[0].description.trim(),
            description.split_whitespace().next().unwrap(),
            "{description}"
        );
        let expected = record(&raw);
        assert_eq!(output.references[0].authors, expected.references[0].authors);
        assert_eq!(output.references[0].title, expected.references[0].title);
        assert_eq!(output.references[0].journal, expected.references[0].journal);
        assert_eq!(output.references[0].pubmed, expected.references[0].pubmed);
        let remark = normalized(output.references[0].remark.as_deref().unwrap());
        assert!(remark.contains("Original laboratory context."), "{remark}");
        assert!(
            remark.contains(&normalized(description)),
            "missing original scope: {remark}"
        );
        assert!(
            remark.contains("source") && remark.contains("coverage"),
            "{remark}"
        );
        // Import/export of the emitted bare reference must not append a second
        // ancestral-source remark when there has been no further edit.
        let imported = import_genbank(&exported).unwrap().molecules.remove(0);
        let repeated = record(&export_genbank(&imported).unwrap());
        assert_eq!(
            normalized(repeated.references[0].remark.as_deref().unwrap()),
            remark
        );
        assert_eq!(original.original_genbank.as_deref(), Some(raw.as_str()));
    }
}

#[test]
fn deleted_and_unscoped_sources_keep_reference_numbers_and_feature_citations() {
    let raw = SOURCE.replace("REFERENCE   1  (bases 3 to 8; 10 to 12)\n", concat!(
        "REFERENCE   1  (bases 3 to 4)\n",
        "  AUTHORS   Deleted,A.\n  TITLE     Deleted interval\n  JOURNAL   Unpublished\n",
        "REFERENCE   2  (sites)\n",
        "  AUTHORS   Annotator,A.\n  TITLE     Annotation source\n  JOURNAL   Unpublished\n",
        "REFERENCE   3  (bases 7 to 12)\n"))
        .replace("/label=\"literal feature\"", "/label=\"literal feature\"\n                     /citation=[1]\n                     /citation=[2]\n                     /citation=[3]");
    let original = import_genbank(&raw).unwrap().molecules.remove(0);
    let deleted = original.replace(2, 4, "").unwrap();
    for edited in [deleted.clone(), deleted.replace(0, 0, "TT").unwrap()] {
        let output = record(&export_genbank(&edited).unwrap());
        assert_eq!(output.references.len(), 3);
        assert_eq!(output.references[0].description.trim(), "1");
        assert_eq!(output.references[1].description.trim(), "2");
        let expected = if edited.sequence.len() == 10 {
            "3  (bases 5 to 10)"
        } else {
            "3  (bases 7 to 12)"
        };
        assert_eq!(output.references[2].description, expected);
        // Qualifier quoting is a writer choice; the parsed citation identities
        // must survive unchanged and continue to name the bibliography entries.
        assert_eq!(output.features.len(), 1);
        assert_eq!(
            output.features[0]
                .qualifier_values("citation")
                .collect::<Vec<_>>(),
            ["[1]", "[2]", "[3]"]
        );
        let remark = normalized(output.references[0].remark.as_deref().unwrap());
        assert!(remark.contains("3 to 4"), "{remark}");
        assert_eq!(edited.original_genbank.as_deref(), Some(raw.as_str()));
    }
    assert_eq!(original.sequence.as_str(), "AAAACCCCGGGG");
}

#[test]
fn repeated_bases_keywords_are_fully_consumed_and_projected() {
    let raw = SOURCE.replace("3 to 8; 10 to 12", "3 to 8; bases 10 to 12");
    let original = import_genbank(&raw).unwrap().molecules.remove(0);
    let output = record(&export_genbank(&original.replace(0, 0, "TT").unwrap()).unwrap());
    assert_eq!(
        output.references[0].description,
        "1  (bases 5 to 10; 12 to 14)"
    );
}

#[test]
fn source_only_citations_do_not_disable_export_proof_checks() {
    let raw = SOURCE.replace("1  (bases 3 to 8; 10 to 12)", "1  (sites)");
    let original = import_genbank(&raw).unwrap().molecules.remove(0);
    let edited = original.replace(0, 0, "TT").unwrap();
    let mut changed_raw = edited.clone();
    changed_raw.original_genbank = Some(raw.replace("Literal citation", "Changed citation"));
    assert!(export_genbank(&changed_raw).is_err());
    let mut changed_dna = edited.clone();
    changed_dna.sequence = Sequence::dna("AAAAAACCCCGGGG").unwrap();
    assert!(export_genbank(&changed_dna).is_err());
    let mut bad_index = edited.clone();
    bad_index
        .metadata
        .insert("source_record_index".into(), "1".into());
    assert!(export_genbank(&bad_index).is_err());
}

#[test]
fn rna_source_only_citation_export_uses_rna_baseline_and_clears_ancestral_accession() {
    let raw = SOURCE
        .replace("DNA     linear", "RNA     linear")
        .replace("aaaaccccgg gg", "uaaaccccgg gg")
        .replace("1  (bases 3 to 8; 10 to 12)", "1  (sites)");
    let original = import_genbank(&raw).unwrap().molecules.remove(0);
    assert_eq!(original.sequence.alphabet(), Alphabet::Rna);
    assert_eq!(original.sequence.as_str(), "UAAACCCCGGGG");
    let edited = original.replace(0, 0, "UU").unwrap();
    assert_eq!(edited.sequence.as_str(), "UUUAAACCCCGGGG");
    let exported = export_genbank(&edited).unwrap();
    let output = record(&exported);
    assert_eq!(output.references[0].description.trim(), "1");
    assert!(output.accession.is_none() && output.version.is_none());
    assert!(
        output
            .comments
            .iter()
            .any(|comment| comment.contains("AB123456.2"))
    );
    assert_eq!(edited.original_genbank.as_deref(), Some(raw.as_str()));
    let restored = import_genbank(&exported).unwrap().molecules.remove(0);
    assert_eq!(restored.sequence, edited.sequence);
}
