use collagene_core::alignment_formats::*;
use collagene_model::{Alphabet, Topology};

fn import(text: &str, format: AlignmentFormat) -> AlignmentDocument {
    import_alignment(text.as_bytes(), format, Alphabet::Dna, "fixture").unwrap()
}
const FASTA: &str = ">alpha first row\nAC-G\nT-\n>beta second row\nA-TGT-\n";
const CLUSTAL: &str = "CLUSTAL W (2.1) multiple sequence alignment\n\nalpha  AC-G 3\nbeta   A-TG 3\n       *  *\n\nalpha  T- 4\nbeta   T- 4\n       * \n";
const STOCKHOLM: &str = "# STOCKHOLM 1.0\n#=GF ID fixture\n#=GF CC first comment\n#=GF CC second comment\n#=GS alpha DE source one\n#=GS beta   DE source two\nalpha AC.G\n#=GR alpha PP 99.8\nbeta A-TG\n#=GR beta PP 9.88\n#=GC SS_cons HH.H\n\nalpha T-\n#=GR alpha PP 7.\nbeta T.\n#=GR beta PP 7.\n#=GC SS_cons H.\n//\n";

#[test]
fn aligned_fasta_and_ungapped_mapping_are_explicit() {
    let doc = import(FASTA, AlignmentFormat::AlignedFasta);
    assert_eq!(doc.column_count(), 6);
    assert_eq!(doc.rows[0].aligned, "AC-GT-");
    assert_eq!(doc.rows[0].description, "first row");
    assert!(doc.source_is_current().unwrap());
    assert_eq!(doc.original_bytes().unwrap(), FASTA.as_bytes());
    let id = &doc.rows[0].id;
    for (column, position) in [Some(0), Some(1), None, Some(2), Some(3), None]
        .into_iter()
        .enumerate()
    {
        assert_eq!(
            doc.column_to_ungapped(id, column)
                .unwrap()
                .map(|p| p.position),
            position
        );
    }
    assert!(doc.column_to_ungapped(id, 6).is_err());
    let extracted = doc.extract_rows(&[id.clone()]).unwrap();
    assert_eq!(extracted[0].sequence.as_str(), "ACGT");
    assert_eq!(extracted[0].topology, Topology::Linear);
    assert!(extracted[0].features.is_empty());
    assert!(extracted[0].original_genbank.is_none());
}
#[test]
fn clustal_numbered_blocks_consensus_and_exact_evidence() {
    let doc = import(CLUSTAL, AlignmentFormat::Clustal);
    assert_eq!(doc.rows[0].aligned, "AC-GT-");
    assert_eq!(doc.metadata.clustal_consensus.as_deref(), Some("*  ** "));
    let out = export_alignment(&doc, AlignmentFormat::Clustal).unwrap();
    let round = import(&out.text, AlignmentFormat::Clustal);
    assert_eq!(
        round.metadata.clustal_consensus,
        doc.metadata.clustal_consensus
    );
    assert_eq!(round.rows[1].aligned, doc.rows[1].aligned);
    assert!(out.renamings.is_empty());
    assert!(
        import_alignment(
            CLUSTAL.replace("AC-G 3", "AC-G 4").as_bytes(),
            AlignmentFormat::Clustal,
            Alphabet::Dna,
            "bad"
        )
        .is_err()
    );
}
#[test]
fn stockholm_annotations_and_dot_gaps_are_retained() {
    let doc = import(STOCKHOLM, AlignmentFormat::Stockholm);
    assert_eq!(doc.rows[0].aligned, "AC-GT-");
    assert_eq!(doc.rows[1].aligned, "A-TGT-");
    assert_eq!(
        doc.metadata
            .annotations
            .iter()
            .filter(|a| a.scope == "GF" && a.tag == "CC")
            .count(),
        2
    );
    assert_eq!(
        doc.metadata
            .annotations
            .iter()
            .find(|a| a.scope == "GR" && a.row_identifier.as_deref() == Some("alpha"))
            .unwrap()
            .value,
        "99.87."
    );
    assert_eq!(
        doc.metadata
            .annotations
            .iter()
            .find(|a| a.scope == "GC")
            .unwrap()
            .value,
        "HH.HH."
    );
    let out = export_alignment(&doc, AlignmentFormat::Stockholm).unwrap();
    let round = import(&out.text, AlignmentFormat::Stockholm);
    assert_eq!(doc.metadata.annotations, round.metadata.annotations);
    assert_eq!(doc.original_bytes().unwrap(), STOCKHOLM.as_bytes());
    let fasta = export_alignment(&doc, AlignmentFormat::AlignedFasta).unwrap();
    assert!(
        fasta
            .losses
            .iter()
            .any(|s| s.contains("Stockholm annotations"))
    );
}
#[test]
fn explicit_alphabet_and_gap_profiles() {
    assert!(
        import_alignment(
            b">a\nACGT\n",
            AlignmentFormat::AlignedFasta,
            Alphabet::Protein,
            "protein"
        )
        .is_ok()
    );
    assert!(
        import_alignment(
            b">a\nACGU\n",
            AlignmentFormat::AlignedFasta,
            Alphabet::Rna,
            "rna"
        )
        .is_ok()
    );
    assert!(
        import_alignment(
            b">a\nACGU\n",
            AlignmentFormat::AlignedFasta,
            Alphabet::Dna,
            "dna"
        )
        .is_err()
    );
    assert!(
        import_alignment(
            b">a\nAC.G\n",
            AlignmentFormat::AlignedFasta,
            Alphabet::Dna,
            "dots"
        )
        .is_err()
    );
    assert!(
        import_alignment(
            b"CLUSTAL\n\na  AC.G\n",
            AlignmentFormat::Clustal,
            Alphabet::Dna,
            "dots"
        )
        .is_err()
    );
    assert!(
        import_alignment(
            b">a\nA?G\n",
            AlignmentFormat::AlignedFasta,
            Alphabet::Dna,
            "bad"
        )
        .is_err()
    );
}
#[test]
fn malformed_framing_order_lengths_and_annotation_references_rejected() {
    for input in [
        ">a\nAA\n>b\nA\n",
        ">a\nAA\n>a\nAA\n",
        "AA\n>a\nAA\n",
        ">a\n\n",
    ] {
        assert!(
            import_alignment(
                input.as_bytes(),
                AlignmentFormat::AlignedFasta,
                Alphabet::Dna,
                "bad"
            )
            .is_err(),
            "{input}"
        );
    }
    for input in [
        "# STOCKHOLM 1.0\na AA\n",
        "# STOCKHOLM 1.0\na AA\n//\na AA\n",
        "# STOCKHOLM 1.0\na AA\n#=GC SS H\n//\n",
        "# STOCKHOLM 1.0\na AA\n#=GR absent SS HH\n//\n",
        "# STOCKHOLM 1.0\na AA\n#=GS absent DE text\n//\n",
        "# STOCKHOLM 1.0\na AA\nb AA\n\nb TT\na TT\n//\n",
        "# STOCKHOLM 1.0\na AA\n#=GC SS HH\n\na TT\n//\n",
        "# STOCKHOLM 1.0\na AA\n#=GC SS HH\n#=GC SS HH\n//\n",
        "# STOCKHOLM 1.0\na AA\n#=XY bad\n//\n",
    ] {
        assert!(
            import_alignment(
                input.as_bytes(),
                AlignmentFormat::Stockholm,
                Alphabet::Dna,
                "bad"
            )
            .is_err(),
            "{input}"
        );
    }
    assert!(
        import_alignment(
            b"CLUSTAL\n\na  AA\nb  AA\n\nb  TT\na  TT\n",
            AlignmentFormat::Clustal,
            Alphabet::Dna,
            "bad"
        )
        .is_err()
    );
}
#[test]
fn workspace_evidence_checks_and_edits() {
    let mut doc = import(FASTA, AlignmentFormat::AlignedFasta);
    assert_eq!(load_alignment(&save_alignment(&doc).unwrap()).unwrap(), doc);
    doc.rows[0].aligned = "AT-GT-".into();
    assert!(!doc.source_is_current().unwrap());
    let loaded = load_alignment(&save_alignment(&doc).unwrap()).unwrap();
    assert!(!loaded.source_is_current().unwrap());
    assert_eq!(loaded.original_bytes().unwrap(), FASTA.as_bytes());
    let mut broken = doc.clone();
    broken.rows[0].id = "forged".into();
    assert!(save_alignment(&broken).is_err());
    let mut broken = doc.clone();
    broken.source.blake3 = "0".repeat(64);
    assert!(save_alignment(&broken).is_err());
    let mut raw = serde_json::to_value(doc).unwrap();
    raw["source_is_current"] = true.into();
    assert!(load_alignment(&serde_json::to_vec(&raw).unwrap()).is_err());
    let mut annotated = import(STOCKHOLM, AlignmentFormat::Stockholm);
    annotated.rows[0].aligned = "AT-GT-".into();
    assert!(save_alignment(&annotated).is_err());
    annotated
        .metadata
        .annotations
        .retain(|a| a.scope != "GR" && a.scope != "GC");
    assert!(!annotated.source_is_current().unwrap());
}
#[test]
fn all_gap_rows_remain_alignment_rows_but_cannot_be_molecules() {
    let doc = import(">a\n--\n>b\nAC\n", AlignmentFormat::AlignedFasta);
    assert!(doc.extract_rows(&[doc.rows[0].id.clone()]).is_err());
    assert!(
        doc.column_to_ungapped(&doc.rows[0].id, 1)
            .unwrap()
            .is_none()
    );
}
#[test]
fn limits_are_bounded_before_parsing() {
    let oversized = vec![b'A'; MAX_ALIGNMENT_BYTES + 1];
    assert!(
        import_alignment(
            &oversized,
            AlignmentFormat::AlignedFasta,
            Alphabet::Dna,
            "large"
        )
        .is_err()
    );
    let long = format!(">a\n{}\n", "A".repeat(MAX_ALIGNMENT_LINE + 1));
    assert!(
        import_alignment(
            long.as_bytes(),
            AlignmentFormat::AlignedFasta,
            Alphabet::Dna,
            "long"
        )
        .is_err()
    );
    let rows = (0..=MAX_ALIGNMENT_ROWS)
        .map(|i| format!(">r{i}\nA\n"))
        .collect::<String>();
    assert!(
        import_alignment(
            rows.as_bytes(),
            AlignmentFormat::AlignedFasta,
            Alphabet::Dna,
            "rows"
        )
        .is_err()
    );
}
#[test]
fn independent_gap_position_oracle_and_all_exports() {
    // Exhaust all gap patterns of width eight; oracle is a separately maintained residue counter.
    for mask in 0u16..256 {
        let row: String = (0..8)
            .map(|c| if mask & (1 << c) == 0 { 'A' } else { '-' })
            .collect();
        let doc = import(&format!(">r\n{row}\n"), AlignmentFormat::AlignedFasta);
        let mut n = 0;
        for (column, c) in row.chars().enumerate() {
            let actual = doc.column_to_ungapped(&doc.rows[0].id, column).unwrap();
            if c == '-' {
                assert!(actual.is_none());
            } else {
                assert_eq!(actual.unwrap().position, n);
                n += 1;
            }
        }
        for format in [
            AlignmentFormat::AlignedFasta,
            AlignmentFormat::Clustal,
            AlignmentFormat::Stockholm,
        ] {
            let out = export_alignment(&doc, format).unwrap();
            let round = import(&out.text, format);
            assert_eq!(round.rows[0].aligned, row);
        }
    }
}
#[test]
fn identifiers_are_not_silently_truncated_or_renamed() {
    let id = "name".repeat(40);
    let doc = import(&format!(">{id}\nACGT\n"), AlignmentFormat::AlignedFasta);
    for format in [
        AlignmentFormat::AlignedFasta,
        AlignmentFormat::Clustal,
        AlignmentFormat::Stockholm,
    ] {
        let out = export_alignment(&doc, format).unwrap();
        assert!(out.renamings.is_empty());
        assert_eq!(import(&out.text, format).rows[0].identifier, id);
    }
    let mut doc = doc;
    doc.rows[0].identifier = "#reserved".into();
    assert!(export_alignment(&doc, AlignmentFormat::Stockholm).is_err());
}

#[test]
fn wrapped_export_preserves_column_metadata_across_80_column_boundary() {
    let source = format!(
        "# STOCKHOLM 1.0\na {}\n#=GR a PP {}\nb {}\n#=GR b PP {}\n#=GC RF {}\n//\n",
        "AC-".repeat(30),
        "9.8".repeat(30),
        "-GT".repeat(30),
        ".87".repeat(30),
        "H_H".repeat(30)
    );
    let doc = import(&source, AlignmentFormat::Stockholm);
    let map = doc.row_position_map(&doc.rows[0].id).unwrap();
    assert_eq!(map.len(), 90);
    assert_eq!(map[88].unwrap().position, 59);
    assert!(map[89].is_none());
    let exported = export_alignment(&doc, AlignmentFormat::Stockholm).unwrap();
    let round = import(&exported.text, AlignmentFormat::Stockholm);
    assert_eq!(round.metadata.annotations, doc.metadata.annotations);
    assert_eq!(round.rows[0].aligned, doc.rows[0].aligned);
    let clustal = export_alignment(&doc, AlignmentFormat::Clustal).unwrap();
    assert_eq!(
        import(&clustal.text, AlignmentFormat::Clustal).rows[0].aligned,
        doc.rows[0].aligned
    );
}

#[test]
fn annotation_delimiters_can_be_spaces_or_tabs() {
    let doc = import(
        "# STOCKHOLM 1.0\n#=GF\t  ID name\n#=GS\t a\t DE text\na AC\n#=GR\ta\tPP\t99\n#=GC\tRF\tHH\n//\n",
        AlignmentFormat::Stockholm,
    );
    assert_eq!(doc.metadata.annotations.len(), 4);
    assert_eq!(doc.metadata.annotations[0].value, "name");
}
