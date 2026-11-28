use collagene_core::embl::*;
use collagene_model::{Alphabet, Error, Feature, Location, Molecule, Sequence, Topology};

const EMBL: &str = "ID   synthetic; SV 3; circular; genomic DNA; STD; SYN; 16 BP.\nXX\nAC   SYNTHETIC0001;\nDE   Synthetic coordinate fixture,\nDE   not a functional construct.\nXX\nRN   [1]\nRP   1-16\nRA   Independent fixture author;\nRT   \"Synthetic evidence\";\nRL   Unpublished.\nZZ   uninterpreted header evidence\nXX\nFH   Key             Location/Qualifiers\nFH\nFT   CDS             complement(join(13..16,1..4))\nFT                   /label=\"Across origin\"\nFT                   /note=\"first\"\nFT                   /note=\"second line\nFT                   continued\"\nFT                   /pseudo\nFT   misc_feature    <5..>8\nFT                   /note=\"fuzzy ends\"\nXX\nSQ   Sequence 16 BP; 4 A; 4 C; 4 G; 4 T; 0 other;\n     aaaaccccgg ggtttt       16\n//\n";

#[test]
fn modern_embl_preserves_compound_features_qualifiers_and_original_evidence() {
    let imported = import_embl(EMBL).unwrap();
    assert_eq!(imported.molecules.len(), 1);
    assert!(imported.warnings.iter().any(|w| w.contains("ZZ")));
    let mol = &imported.molecules[0];
    assert_eq!(mol.name, "synthetic");
    assert_eq!(mol.topology, Topology::Circular);
    assert_eq!(mol.sequence.as_str(), "AAAACCCCGGGGTTTT");
    assert_eq!(mol.metadata["embl.original_record"], EMBL);
    assert_eq!(
        mol.metadata["definition"],
        "Synthetic coordinate fixture, not a functional construct."
    );
    assert_eq!(mol.metadata["embl.sequence_version"], "3");
    assert!(mol.original_genbank.is_none());
    assert_eq!(
        mol.features[0].location.extract(&mol.sequence).unwrap(),
        "TTTTAAAA"
    );
    assert_eq!(
        mol.features[0].qualifiers,
        vec![
            ("label".into(), Some("Across origin".into())),
            ("note".into(), Some("first".into())),
            ("note".into(), Some("second line\ncontinued".into())),
            ("pseudo".into(), None)
        ]
    );
    assert!(mol.features[1].location.extract(&mol.sequence).is_err());
    assert_eq!(
        mol.features[1].location,
        Location::Range {
            start: 4,
            end: 8,
            before: true,
            after: true
        }
    );
    assert_eq!(export_embl(mol).unwrap(), EMBL);
}

#[test]
fn regenerated_export_uses_current_annotations_and_retains_supported_headers() {
    let mut mol = import_embl(EMBL).unwrap().molecules.remove(0);
    mol.name = "renamed fixture".into();
    mol.metadata
        .insert("definition".into(), "Current definition".into());
    mol.features[0]
        .qualifiers
        .push(("experiment".into(), Some("independent annotation".into())));
    let exported = export_embl(&mol).unwrap();
    assert!(
        exported.starts_with("ID   renamed_fixture; SV 3; circular; genomic DNA; STD; SYN; 16 BP.")
    );
    assert!(
        exported.contains("AC   SYNTHETIC0001;")
            && exported.contains("RP   1-16")
            && exported.contains("ZZ   uninterpreted header evidence")
    );
    let restored = import_embl(&exported).unwrap().molecules.remove(0);
    assert_eq!(restored.sequence, mol.sequence);
    assert_eq!(restored.features, mol.features);
    assert_eq!(restored.metadata["definition"], "Current definition");
}

#[test]
fn multi_record_evidence_is_not_multiplied_and_lengths_are_checked() {
    let report = import_embl(&format!("{EMBL}\n{EMBL}")).unwrap();
    assert_eq!(report.molecules.len(), 2);
    for mol in report.molecules {
        assert_eq!(mol.metadata["embl.original_record"], EMBL);
    }
    for malformed in [
        EMBL.replace("16 BP.", "15 BP."),
        EMBL.replace("Sequence 16 BP;", "Sequence 17 BP;"),
        EMBL.replace("       16", "       15"),
        EMBL.replace("//\n", ""),
        EMBL.replace("ggtttt", "ggttt!"),
    ] {
        assert!(import_embl(&malformed).is_err(), "{malformed}");
    }
}

#[test]
fn edits_with_unmapped_reference_or_unknown_headers_are_not_exported_as_remapped() {
    let original = import_embl(EMBL).unwrap().molecules.remove(0);
    let edited = original.replace(4, 5, "A").unwrap();
    assert!(matches!(export_embl(&edited), Err(Error::Unsupported(_))));
    let simple = "ID   simple; SV 1; linear; genomic DNA; STD; SYN; 8 BP.\nDE   Simple editable record\nFT   misc_feature    3..6\nSQ   Sequence 8 BP; 4 A; 4 C; 0 G; 0 T; 0 other;\n     aaaacccc        8\n//\n";
    let mol = import_embl(simple).unwrap().molecules.remove(0);
    let edited = mol.replace(0, 1, "T").unwrap();
    let exported = export_embl(&edited).unwrap();
    assert!(exported.contains("sequence differs from original EMBL evidence"));
    assert_eq!(
        import_embl(&exported).unwrap().molecules[0]
            .sequence
            .as_str(),
        "TAAACCCC"
    );
}

#[test]
fn rna_encoding_is_explicit_and_proteins_are_not_misidentified() {
    let rna = "ID   rna; SV 1; linear; mRNA; STD; SYN; 8 BP.\nDE   RNA fixture\nSQ   Sequence 8 BP; 2 A; 2 C; 2 G; 2 T; 0 other;\n     acgttgca        8\n//\n";
    let imported = import_embl(rna).unwrap();
    assert_eq!(imported.molecules[0].sequence.alphabet(), Alphabet::Rna);
    assert_eq!(imported.molecules[0].sequence.as_str(), "ACGUUGCA");
    assert!(imported.warnings.iter().any(|w| w.contains("T to U")));
    let mut molecule = imported.molecules[0].clone();
    molecule.name = "renamed".into();
    let exported = export_embl(&molecule).unwrap();
    assert!(exported.contains("acgttgca"));
    assert_eq!(
        import_embl(&exported).unwrap().molecules[0].sequence,
        molecule.sequence
    );
    assert!(import_embl(&rna.replace("mRNA", "PRT")).is_err());
    assert!(
        export_embl(&Molecule {
            sequence: Sequence::new(Alphabet::Protein, "ACDE").unwrap(),
            ..Molecule::dna("protein", "ACGT", Topology::Linear).unwrap()
        })
        .is_err()
    );
}

#[test]
fn native_export_and_legacy_id_round_trip_without_external_programs() {
    let mut mol = Molecule::dna("new molecule", "ACGTN", Topology::Circular).unwrap();
    mol.features.push(Feature {
        kind: "misc_feature".into(),
        location: Location::Order(vec![Location::range(0, 2), Location::range(3, 5)]),
        qualifiers: vec![("note".into(), Some("ordered location".into()))],
        needs_review: false,
    });
    let text = export_embl(&mol).unwrap();
    let restored = import_embl(&text).unwrap().molecules.remove(0);
    assert_eq!(restored.sequence, mol.sequence);
    assert_eq!(restored.features, mol.features);
    assert_eq!(restored.topology, mol.topology);
    let legacy = "ID   legacy standard; DNA; SYN; 4 BP.\nSQ   Sequence 4 BP; 1 A; 1 C; 1 G; 1 T; 0 other;\n     acgt        4\n//\n";
    let imported = import_embl(legacy).unwrap();
    assert_eq!(imported.molecules[0].name, "legacy");
    assert_eq!(imported.molecules[0].topology, Topology::Linear);
    assert!(imported.warnings.iter().any(|w| w.contains("legacy")));
    let legacy_version = legacy
        .replace("DNA; SYN", "circular DNA; SYN")
        .replace("SQ   Sequence", "SV   legacy.7\nSQ   Sequence");
    let mut imported = import_embl(&legacy_version).unwrap().molecules.remove(0);
    assert_eq!(imported.topology, Topology::Circular);
    assert_eq!(imported.metadata["embl.sequence_version"], "7");
    imported.name = "renamed".into();
    let modern = export_embl(&imported).unwrap();
    assert!(modern.starts_with("ID   renamed; SV 7; circular; DNA;"));
    assert!(!modern.contains("SV   legacy.7"));
    assert_eq!(
        import_embl(&modern).unwrap().molecules[0].topology,
        Topology::Circular
    );
}

#[test]
fn location_depth_and_reference_only_records_are_rejected_before_deep_parsing() {
    let nested = format!("{}1..4{}", "complement(".repeat(40), ")".repeat(40));
    let text = format!(
        "ID   deep; SV 1; linear; genomic DNA; STD; SYN; 4 BP.\nFT   misc_feature    {nested}\nSQ   Sequence 4 BP;\n     acgt        4\n//\n"
    );
    assert!(import_embl(&text).is_err());
    let contig =
        "ID   contig; SV 1; linear; genomic DNA; CON; SYN; 4 BP.\nCO   join(OTHER:1..4)\n//\n";
    assert!(matches!(import_embl(contig), Err(Error::Unsupported(_))));
}
