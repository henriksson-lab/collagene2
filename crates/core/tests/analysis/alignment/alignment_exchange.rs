use collagene_core::formats::alignment::profiles::{self as exchange, ExportOptions};
use collagene_core::formats::alignment::{self as native, AlignmentDocument, AlignmentFormat as F};
use collagene_model::{Alphabet, Error};

fn import(s: &str, f: F, a: Alphabet) -> AlignmentDocument {
    exchange::import_alignment(s.as_bytes(), f, a, "Independent exchange fixture").unwrap()
}
const STRICT_I: &str = "2 6\nH. SapiensAC-G\nMouse     A-TG\n\nT-\nT-\n";
const STRICT_S: &str = "2 6\nH. SapiensAC-G\nT-\nMouse     A-TG\nT-\n";
const RELAXED_I: &str = "2 6\nlong_first_identifier AC-G\nsecond A-TG\n\nT-\nT-\n";
const MSF: &str = "!!NA_MULTIPLE_ALIGNMENT 1.0\n\nfixture MSF: 6 Type: N Check: 2706 ..\n Name: alpha Len: 6 Check: 1317 Weight: 0.5\n Name: beta Len: 6 Check: 1389 Weight: 1.0\n//\n\n 1 4\nalpha AC.G\nbeta A.TG\n\n 5 6\nalpha T.\nbeta T.\n";
const NEXUS: &str = "#NEXUS\n[outer [nested] comment]\nBEGIN TAXA;\n DIMENSIONS NTAX=2;\n TAXLABELS Homo_sapiens 'mouse''s';\nEND;\nBEGIN CHARACTERS;\n DIMENSIONS NCHAR=6;\n FORMAT DATATYPE=DNA GAP=- MISSING=? MATCHCHAR=. INTERLEAVE=YES;\n MATRIX\n Homo_sapiens AC-G\n 'mouse''s' .-..\n\n Homo_sapiens T-\n 'mouse''s' ..\n ;\nEND;\nBEGIN TREES;\n TREE 'name;with;semicolons'=(Homo_sapiens,'mouse''s');\nEND;\nBEGIN PRIVATE;\n arbitrary opaque=(one,two);\nEND;\n";

#[test]
fn phylip_profiles_preserve_real_columns_names_and_native_originals() {
    let formats = [
        (STRICT_I, F::PhylipInterleaved),
        (STRICT_S, F::PhylipSequential),
    ];
    for (s, f) in formats {
        let d = import(s, f, Alphabet::Dna);
        assert_eq!(d.rows[0].identifier, "H. Sapiens");
        assert_eq!(d.rows[0].aligned, "AC-GT-");
        assert_eq!(d.rows[1].aligned, "A-TGT-");
        assert_eq!(d.column_to_ungapped(&d.rows[0].id, 2).unwrap(), None);
        assert_eq!(
            d.column_to_ungapped(&d.rows[0].id, 4)
                .unwrap()
                .unwrap()
                .position,
            3
        );
        let bytes = native::save_alignment(&d).unwrap();
        let loaded = native::load_alignment(&bytes).unwrap();
        assert_eq!(loaded, d);
        assert_eq!(loaded.original_bytes().unwrap(), s.as_bytes());
        assert!(loaded.source_is_current().unwrap());
        let export = exchange::export_alignment(&d, f).unwrap();
        let round = import(&export.text, f, Alphabet::Dna);
        assert_eq!(round.rows[0].identifier, "H. Sapiens");
        assert_eq!(round.rows[1].aligned, "A-TGT-");
        assert!(export.renamings.is_empty());
    }
    let relaxed = import(RELAXED_I, F::PhylipRelaxedInterleaved, Alphabet::Dna);
    assert_eq!(relaxed.rows[0].identifier, "long_first_identifier");
    let sequential = exchange::export_alignment(&relaxed, F::PhylipRelaxedSequential).unwrap();
    let round = import(&sequential.text, F::PhylipRelaxedSequential, Alphabet::Dna);
    assert_eq!(round.rows[0].aligned, "AC-GT-");
    assert_eq!(round.rows[1].identifier, "second");
}

#[test]
fn strict_phylip_truncation_is_opt_in_reported_and_never_merges_taxa() {
    let d = import(RELAXED_I, F::PhylipRelaxedInterleaved, Alphabet::Dna);
    assert!(exchange::export_alignment(&d, F::PhylipSequential).is_err());
    let options = ExportOptions {
        truncate_phylip_identifiers: true,
        block_width: 3,
    };
    let out = exchange::export_alignment_with_options(&d, F::PhylipInterleaved, &options).unwrap();
    assert_eq!(out.renamings["long_first_identifier"], "long_first");
    assert!(out.losses.iter().any(|s| s.contains("truncation")));
    let round = import(&out.text, F::PhylipInterleaved, Alphabet::Dna);
    assert_eq!(round.rows[0].aligned, "AC-GT-");
    let mut collision = d.clone();
    collision.rows[1].identifier = "long_first_other".into();
    assert!(
        exchange::export_alignment_with_options(&collision, F::PhylipInterleaved, &options)
            .is_err()
    );
    // Names containing spaces fit the ten-column profile but cannot be silently
    // rewritten for whitespace-delimited relaxed PHYLIP or MSF.
    let spaced = import(STRICT_I, F::PhylipInterleaved, Alphabet::Dna);
    assert!(exchange::export_alignment(&spaced, F::PhylipRelaxedInterleaved).is_err());
    assert!(exchange::export_alignment(&spaced, F::Msf).is_err());
    for legacy in [F::AlignedFasta, F::Clustal, F::Stockholm] {
        assert!(native::export_alignment(&spaced, legacy).is_err());
    }
}

#[test]
fn phylip_rejects_bad_dimensions_block_geometry_and_ambiguous_extensions() {
    for malformed in [
        "2 6\na         ACGT\nb         AC-G\nT-\n",
        "2 4\na         ACGT\na         ACGT\n",
        "2 4\na         ACGT\nb         ACGT\nextra\n",
        "1 4\na         AC.G\n",
        "1 4\na         AC?G\n",
        "1 4\na ACGT\n",
    ] {
        assert!(
            exchange::import_alignment(
                malformed.as_bytes(),
                F::PhylipInterleaved,
                Alphabet::Dna,
                "bad"
            )
            .is_err(),
            "{malformed}"
        );
    }
    let huge = "1000 18446744073709551615\n";
    assert!(
        exchange::import_alignment(huge.as_bytes(), F::PhylipSequential, Alphabet::Dna, "bad")
            .is_err()
    );
    // Explicit alphabets, never a heuristic based on A/C/G/T-looking protein data.
    assert_eq!(
        import("1 4\nrna       ACGU\n", F::PhylipSequential, Alphabet::Rna).alphabet,
        Alphabet::Rna
    );
    assert_eq!(
        import(
            "1 4\nprotein   ACGT\n",
            F::PhylipSequential,
            Alphabet::Protein
        )
        .alphabet,
        Alphabet::Protein
    );
    assert!(
        exchange::import_alignment(
            b"1 4\nrna       ACGU\n",
            F::PhylipSequential,
            Alphabet::Dna,
            "bad"
        )
        .is_err()
    );
}

#[test]
fn msf_checksums_weights_and_gap_mapping_are_independently_pinned() {
    assert_eq!(exchange::msf_checksum("AC.GT."), 1317);
    assert_eq!(exchange::msf_checksum("A.TGT."), 1389);
    let d = import(MSF, F::Msf, Alphabet::Dna);
    assert_eq!(d.rows[0].aligned, "AC-GT-");
    assert_eq!(d.rows[1].aligned, "A-TGT-");
    assert!(
        d.metadata
            .annotations
            .iter()
            .any(|a| a.scope == "GS" && a.tag == "WT" && a.value == "0.5")
    );
    assert_eq!(d.original_bytes().unwrap(), MSF.as_bytes());
    let out = exchange::export_alignment(&d, F::Msf).unwrap();
    assert!(out.text.contains("Check: 2706"));
    let round = import(&out.text, F::Msf, Alphabet::Dna);
    assert_eq!(round.rows[0].aligned, d.rows[0].aligned);
    assert!(out.text.contains("Weight: 0.5"));
    assert!(
        exchange::import_alignment(
            MSF.replace("alpha AC.G", "alpha TC.G").as_bytes(),
            F::Msf,
            Alphabet::Dna,
            "tampered"
        )
        .is_err()
    );
    assert!(
        exchange::import_alignment(MSF.as_bytes(), F::Msf, Alphabet::Protein, "wrong type")
            .is_err()
    );
    let unchecked = MSF
        .replace("Check: 2706", "Check: 0")
        .replace("Check: 1317", "Check: 0")
        .replace("Check: 1389", "Check: 0");
    let zero = import(&unchecked, F::Msf, Alphabet::Dna);
    assert!(
        zero.metadata
            .annotations
            .iter()
            .any(|a| a.value.contains("unspecified"))
    );
    assert!(
        exchange::import_alignment(
            MSF.replace("Weight: 0.5", "Weight: NaN").as_bytes(),
            F::Msf,
            Alphabet::Dna,
            "bad"
        )
        .is_err()
    );
}

#[test]
fn nexus_quotes_matchchar_and_opaque_blocks_retain_exact_evidence() {
    let d = import(NEXUS, F::Nexus, Alphabet::Dna);
    assert_eq!(d.rows[0].identifier, "Homo sapiens");
    assert_eq!(d.rows[1].identifier, "mouse's");
    assert_eq!(d.rows[0].aligned, "AC-GT-");
    assert_eq!(d.rows[1].aligned, "A--GT-");
    assert_eq!(d.column_to_ungapped(&d.rows[1].id, 2).unwrap(), None);
    assert_eq!(d.original_bytes().unwrap(), NEXUS.as_bytes());
    let saved = native::save_alignment(&d).unwrap();
    let reopened = native::load_alignment(&saved).unwrap();
    assert_eq!(reopened, d);
    assert!(
        d.metadata
            .annotations
            .iter()
            .any(|a| a.value.contains("TREES") && a.value.contains("not interpreted"))
    );
    assert!(
        d.metadata
            .annotations
            .iter()
            .any(|a| a.value.contains("PRIVATE") && a.value.contains("not interpreted"))
    );
    let out = exchange::export_alignment(&d, F::Nexus).unwrap();
    assert!(!out.text.contains("BEGIN TREES"));
    assert!(out.text.contains("'mouse''s'"));
    assert!(out.losses.iter().any(|s| s.contains("TREES")));
    let round = import(&out.text, F::Nexus, Alphabet::Dna);
    assert_eq!(round.rows[0].identifier, "Homo sapiens");
    assert_eq!(round.rows[1].aligned, "A--GT-");
    let underscore = import(
        "#NEXUS\nBEGIN DATA; DIMENSIONS NTAX=1 NCHAR=4; FORMAT DATATYPE=DNA; MATRIX\n'alpha_beta' ACGT\n; END;",
        F::Nexus,
        Alphabet::Dna,
    );
    assert_eq!(underscore.rows[0].identifier, "alpha_beta");
}

#[test]
fn nexus_sequential_rna_protein_and_unnamed_continuations_are_explicit() {
    let rna = "#NEXUS\nBEGIN DATA; DIMENSIONS NTAX=2 NCHAR=4; FORMAT DATATYPE=RNA MATCHCHAR=.; MATRIX\na ACGU\nb .C..\n; END;";
    assert_eq!(import(rna, F::Nexus, Alphabet::Rna).rows[1].aligned, "ACGU");
    assert!(exchange::import_alignment(rna.as_bytes(), F::Nexus, Alphabet::Dna, "bad").is_err());
    let protein = rna
        .replace("RNA", "PROTEIN")
        .replace("ACGU", "ACDE")
        .replace(".C..", ".A..");
    assert_eq!(
        import(&protein, F::Nexus, Alphabet::Protein).rows[1].aligned,
        "AADE"
    );
    let unnamed = "#NEXUS\nBEGIN DATA; DIMENSIONS NTAX=2 NCHAR=6; FORMAT DATATYPE=DNA INTERLEAVE; MATRIX\na AC-G\nb A-TG\n\nT-\nT-\n; END;";
    assert_eq!(
        import(unnamed, F::Nexus, Alphabet::Dna).rows[1].aligned,
        "A-TGT-"
    );
}

#[test]
fn nexus_unsupported_matrix_semantics_and_malformed_framing_fail_closed() {
    let valid = "#NEXUS\nBEGIN DATA; DIMENSIONS NTAX=1 NCHAR=4; FORMAT DATATYPE=DNA; MATRIX\na ACGT\n; END;";
    for malformed in [
        valid.replace("ACGT", "A?GT"),
        valid.replace("DATATYPE=DNA", "DATATYPE=DNA TRANSPOSE"),
        valid.replace("DATATYPE=DNA", "DATATYPE=MIXED(DNA:1-4)"),
        valid.replace("ACGT", "{AG}CGT"),
        valid.replace("NCHAR=4", "NCHAR=3"),
        valid.replace("; END;", ";"),
        valid.replace("a ACGT", "'unterminated ACGT"),
        format!("{valid}\nBEGIN DATA; END;"),
    ] {
        assert!(
            exchange::import_alignment(malformed.as_bytes(), F::Nexus, Alphabet::Dna, "bad")
                .is_err(),
            "{malformed}"
        );
    }
    let mut d = import(valid, F::Nexus, Alphabet::Dna);
    d.rows[0].aligned = "ACGA".into();
    assert!(!d.source_is_current().unwrap());
    assert_eq!(d.original_bytes().unwrap(), valid.as_bytes());
    let out = exchange::export_alignment(&d, F::Nexus).unwrap();
    assert!(out.text.contains("ACGA"));
}

#[test]
fn input_output_and_profile_budgets_reject_without_partial_documents() {
    let oversized = vec![b' '; native::MAX_ALIGNMENT_BYTES + 1];
    assert!(matches!(
        exchange::import_alignment(&oversized, F::Msf, Alphabet::Dna, "bad"),
        Err(Error::TooLarge)
    ));
    let d = import(RELAXED_I, F::PhylipRelaxedInterleaved, Alphabet::Dna);
    for width in [0, 1001] {
        assert!(
            exchange::export_alignment_with_options(
                &d,
                F::Nexus,
                &ExportOptions {
                    truncate_phylip_identifiers: false,
                    block_width: width
                }
            )
            .is_err()
        );
    }
    let non_ascii = b"1 4\na         AC\xffT\n";
    assert!(
        exchange::import_alignment(non_ascii, F::PhylipSequential, Alphabet::Dna, "bad").is_err()
    );
    assert!(
        exchange::import_alignment(b"#NEXUS [unclosed", F::Nexus, Alphabet::Dna, "bad").is_err()
    );
}
