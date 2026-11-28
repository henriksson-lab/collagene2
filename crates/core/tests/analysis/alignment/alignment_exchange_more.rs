use collagene_core::formats::alignment::{
    self as base, AlignmentDocument, AlignmentFormat as Format,
};
use collagene_core::formats::alignment::{pir_gde_selex as exchange, profiles::ExportOptions};
use collagene_model::{Alphabet, Error, Topology};
fn read(text: &str, format: Format, alphabet: Alphabet) -> AlignmentDocument {
    exchange::import_alignment(
        text.as_bytes(),
        format,
        alphabet,
        "independent literal alignment",
    )
    .unwrap()
}
fn rows(doc: &AlignmentDocument) -> Vec<(&str, &str)> {
    doc.rows
        .iter()
        .map(|r| (r.identifier.as_str(), r.aligned.as_str()))
        .collect()
}
const PIR: &str = ">P1;first\ncomplete protein source title\nACD-EF GH*\nC;Species: documentary source text\nF;2/Binding site: annotation only\n>F1;second\nfragment title\nACDE-FGH\n*\n";
#[test]
fn pir_type_description_text_and_terminal_star_have_literal_meanings() {
    let doc = read(PIR, Format::Pir, Alphabet::Protein);
    assert_eq!(
        rows(&doc),
        vec![("first", "ACD-EFGH"), ("second", "ACDE-FGH")]
    );
    assert_eq!(doc.rows[0].description, "complete protein source title");
    assert!(doc.metadata.annotations.iter().any(|a| a.tag == "PIR_TYPE"
        && a.row_identifier.as_deref() == Some("second")
        && a.value == "F1"));
    assert!(
        doc.metadata
            .annotations
            .iter()
            .any(|a| a.tag == "PIR_TEXT" && a.value == "F;2/Binding site: annotation only")
    );
    assert_eq!(doc.original_bytes().unwrap(), PIR.as_bytes());
    assert!(doc.source_is_current().unwrap());
    let out = exchange::export_alignment(&doc, Format::Pir).unwrap();
    let again = read(&out.text, Format::Pir, Alphabet::Protein);
    assert_eq!(rows(&again), rows(&doc));
    assert_eq!(again.rows[1].description, "fragment title");
    assert!(out.text.contains(">P1;first") && out.text.contains(">F1;second"));
    let dna = read(">DL;HLA*01\nDNA title\nAC-G*\n", Format::Pir, Alphabet::Dna);
    assert_eq!(rows(&dna), vec![("HLA*01", "AC-G")]);
    assert_eq!(
        dna.extract_rows(&[dna.rows[0].id.clone()]).unwrap()[0]
            .sequence
            .as_str(),
        "ACG"
    );
    let rna = read(">RL;r\nRNA title\nAU-G*\n", Format::Pir, Alphabet::Rna);
    assert_eq!(rna.alphabet, Alphabet::Rna);
    assert!(
        exchange::import_alignment(PIR.as_bytes(), Format::Pir, Alphabet::Dna, "wrong alphabet")
            .is_err()
    );
}
#[test]
fn pir_rejects_circular_state_and_stop_confusion_without_discarding_biology() {
    for text in [">DC;x\ncircularDNA\nACGT*\n", ">RC;x\ncircularRNA\nACGU*\n"] {
        assert!(matches!(
            exchange::import_alignment(text.as_bytes(), Format::Pir, Alphabet::Dna, "circular"),
            Err(Error::Unsupported(_))
        ));
    }
    for text in [
        ">P1;x\ntitle\nACD",
        ">P1;x\ntitle\nAC**\n",
        ">P1;x\ntitle\nAC. D*\n",
        ">P1;x\ntitle\nAC2D*\n",
    ] {
        assert!(
            exchange::import_alignment(
                text.as_bytes(),
                Format::Pir,
                Alphabet::Protein,
                "malformed"
            )
            .is_err()
        );
    }
    let stopped = read("stop    AC*D\n", Format::Selex, Alphabet::Protein);
    assert_eq!(stopped.rows[0].aligned, "AC*D");
    assert!(exchange::export_alignment(&stopped, Format::Pir).is_err());
}
#[test]
fn gde_offsets_are_explicit_alignment_gaps_and_roundtrip_without_external_coordinates() {
    let source = "#early(0)\nAC--G\n#late(2)\nA~G\n";
    let doc = read(source, Format::Gde, Alphabet::Dna);
    assert_eq!(rows(&doc), vec![("early", "AC--G"), ("late", "--A-G")]);
    assert!(
        doc.metadata
            .annotations
            .iter()
            .any(|a| a.tag == "GDE_OFFSET"
                && a.row_identifier.as_deref() == Some("late")
                && a.value == "2")
    );
    assert!(
        doc.column_to_ungapped(&doc.rows[1].id, 0)
            .unwrap()
            .is_none()
    );
    assert_eq!(
        doc.column_to_ungapped(&doc.rows[1].id, 2)
            .unwrap()
            .unwrap()
            .position,
        0
    );
    assert_eq!(doc.original_bytes().unwrap(), source.as_bytes());
    let out = exchange::export_alignment(&doc, Format::Gde).unwrap();
    assert!(out.text.contains("#late(2)"));
    assert_eq!(
        rows(&read(&out.text, Format::Gde, Alphabet::Dna)),
        rows(&doc)
    );
    let m = &doc.extract_rows(&[doc.rows[1].id.clone()]).unwrap()[0];
    assert_eq!(m.sequence.as_str(), "AG");
    assert_eq!(m.topology, Topology::Linear);
    assert!(m.features.is_empty());
    assert!(m.metadata["alignment_mapping"].contains("no external reference coordinates"));
    let mut changed = doc.clone();
    changed.rows[1].aligned = "C-A-G".into();
    assert!(exchange::export_alignment(&changed, Format::Gde).is_err());
}
#[test]
fn gde_author_manual_rna_example_and_unsupported_profiles_are_not_guessed() {
    // GDE2.2 author's manual Appendix A: two E.coli tRNA flat records. The
    // second source has75residues; native alignment adds one explicit trailinggap.
    let source = "#ECOTRNT4\nGGGUCGUUAGCUCAGUUGGUAGAGCAGUUGACUUUUAAUCAAUUGGNCGCAGGUUCGAAU\nCCUGCACGACCCACCA\n#ECOTRQ1\nUGGGGUAUCGCCAAGCGGUAAGGCACCGGUUUUUGAUACCGGCAUUCCCUGGUUCGAAUC\nCAGGUACCCCAGCCA\n";
    let doc = read(source, Format::Gde, Alphabet::Rna);
    assert_eq!(doc.rows[0].aligned.len(), 76);
    assert_eq!(doc.rows[1].aligned.len(), 76);
    assert!(doc.rows[1].aligned.ends_with("CAGGUACCCCAGCCA-"));
    assert_eq!(doc.original_bytes().unwrap(), source.as_bytes());
    assert!(
        exchange::import_alignment(source.as_bytes(), Format::Gde, Alphabet::Dna, "wrong").is_err()
    );
    for source in [
        "#x(-1)\nAC\n",
        "@mask\n0011\n",
        "\"text\ncomment\n",
        "{name \"tagged\"}\n",
    ] {
        assert!(matches!(
            exchange::import_alignment(
                source.as_bytes(),
                Format::Gde,
                Alphabet::Dna,
                "unsupported"
            ),
            Err(Error::Unsupported(_))
        ));
    }
    assert!(
        exchange::import_alignment(b"%protein\nACGT\n", Format::Gde, Alphabet::Dna, "wrongtype")
            .is_err()
    );
    assert!(
        exchange::import_alignment(
            b"#x(9000000)\nA\n",
            Format::Gde,
            Alphabet::Dna,
            "offsetbound"
        )
        .is_err()
    );
}
const SELEX: &str = "#=ID independent\n#=DE Stored annotation, not a computed confidence\n#=SQ one 1.5 source ACC1 5..7:100 source row one\n#=SQ two - - - 0 source row two\n#=RF   x.xx\n#=CS   H..L\none    Ac g\n#=SS   H..L\n#=SA   1234\ntwo    A--G\n#=SS   H..L\n#=SA   5678\n\n#=RF   xxxx\n#=CS   LLLL\none    U-Ac\n#=SS   LLLL\n#=SA   5678\ntwo    U--C\n#=SS   LLLL\n#=SA   1234\n";
#[test]
fn selex_interleaved_spaces_case_and_opaque_column_row_annotations_roundtrip() {
    let doc = read(SELEX, Format::Selex, Alphabet::Rna);
    assert_eq!(rows(&doc), vec![("one", "AC-GU-AC"), ("two", "A--GU--C")]);
    assert_eq!(doc.rows[0].description, "source row one");
    let find = |scope: &str, row: Option<&str>, tag: &str| {
        doc.metadata
            .annotations
            .iter()
            .find(|a| a.scope == scope && a.row_identifier.as_deref() == row && a.tag == tag)
            .unwrap()
            .value
            .as_str()
    };
    assert_eq!(find("GC", None, "RF"), "x.xxxxxx");
    assert_eq!(find("GR", Some("one"), "SS"), "H..LLLLL");
    assert_eq!(find("GR", Some("one"), "SA"), "12345678");
    assert_eq!(find("GR", Some("one"), "SELEX_CASE"), "ul.lu.ul");
    assert_eq!(
        find("GS", Some("one"), "SELEX_SQ"),
        "1.5 source ACC1 5..7:100 source row one"
    );
    assert_eq!(doc.original_bytes().unwrap(), SELEX.as_bytes());
    let out = exchange::export_alignment_with_options(
        &doc,
        Format::Selex,
        &ExportOptions {
            block_width: 3,
            ..Default::default()
        },
    )
    .unwrap();
    let again = read(&out.text, Format::Selex, Alphabet::Rna);
    assert_eq!(rows(&again), rows(&doc));
    assert_eq!(again.rows[0].description, doc.rows[0].description);
    for annotation in doc
        .metadata
        .annotations
        .iter()
        .filter(|a| a.scope == "GR" || a.scope == "GC")
    {
        assert!(again.metadata.annotations.contains(annotation));
    }
    let molecule = &doc.extract_rows(&[doc.rows[0].id.clone()]).unwrap()[0];
    assert_eq!(molecule.sequence.as_str(), "ACGUAC");
    assert!(molecule.features.is_empty());
    let fasta = exchange::export_alignment(&doc, Format::Pir).unwrap();
    assert!(fasta.losses.iter().any(|loss| loss.contains("Annotations")));
}
#[test]
fn selex_author_guide_variable_width_blocks_are_padded_not_reinterpreted_as_coordinates() {
    // Sean Eddy's HMMER2.2 user guide simple SELEX example, complete data lines.
    let input = "# Example selex file\nseq1 ACGACGACGACG.\nseq2 ..GGGAAAGG.GA\nseq3 UUU..AAAUUU.A\n\nseq1 ..ACG\nseq2 AAGGG\nseq3 AA...UUU\n";
    let doc = read(input, Format::Selex, Alphabet::Rna);
    assert_eq!(doc.rows[0].aligned, "ACGACGACGACG---ACG---");
    assert_eq!(doc.rows[1].aligned, "--GGGAAAGG-GAAAGGG---");
    assert_eq!(doc.rows[2].aligned, "UUU--AAAUUU-AAA---UUU");
    assert_eq!(doc.rows[0].aligned.len(), 21);
    assert!(
        doc.metadata
            .annotations
            .iter()
            .any(|a| a.tag == "SELEX_LAYOUT" && a.value.contains("padded"))
    );
}
#[test]
fn malformed_selex_headers_coverage_tabs_and_alignment_bounds_fail_closed() {
    for source in [
        "#=SS HHH\none ACG\n",
        "one ACG\ntwo ACG\n\ntwo UUU\none UUU\n",
        "#=SQ one 1.0 - - 0 desc\none ACG\ntwo ACG\n",
        "one\tACG\n",
        "#=XX unsupported\none ACG\n",
        "#=ID a\n#=ID b\none ACG\n",
    ] {
        assert!(
            exchange::import_alignment(
                source.as_bytes(),
                Format::Selex,
                Alphabet::Dna,
                "malformed"
            )
            .is_err(),
            "accepted{source}"
        );
    }
    let mismatch = SELEX.replacen("#=CS   LLLL\n", "", 1);
    assert!(
        exchange::import_alignment(
            mismatch.as_bytes(),
            Format::Selex,
            Alphabet::Rna,
            "missingcoverage"
        )
        .is_err()
    );
    assert!(
        exchange::import_alignment(
            &vec![b' '; base::MAX_ALIGNMENT_BYTES + 1],
            Format::Selex,
            Alphabet::Dna,
            "size"
        )
        .is_err()
    );
    let doc = read(PIR, Format::Pir, Alphabet::Protein);
    assert!(
        exchange::export_alignment_with_options(
            &doc,
            Format::Pir,
            &ExportOptions {
                block_width: 0,
                ..Default::default()
            }
        )
        .is_err()
    );
}
