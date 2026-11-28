use base64::{Engine, engine::general_purpose::STANDARD};
use collagene_core::{formats::oligos::*, sequence::analysis::Primer};
use collagene_model::Location;

fn options() -> OligoImportOptions {
    OligoImportOptions {
        delimiter: Delimiter::Comma,
        has_header: true,
        columns: ColumnMapping {
            name: 0,
            ordered_sequence: Some(1),
            annealing_sequence: Some(2),
            tail: Some(3),
            modification: Some(4),
            amount: Some(5),
            unit: Some(6),
        },
    }
}
#[test]
fn mapped_roles_duplicate_names_quantities_and_opaque_modifications_are_distinct() {
    let text = "Name,Ordered,Domain,Tail,Modification,Amount,Unit,Notes\r\nF,GGACGT,ACGT,GG,,25,nmol,plain\r\nF,/5Phos/ACGT,ACGT,,5-prime phosphate,10,uM,modified\r\nR,TTACGA,ACGA,TT,,5,unknown,opaque\r\n";
    let table = import_oligos(text, &options()).unwrap();
    table.validate().unwrap();
    assert_eq!(table.records.len(), 3);
    assert_eq!(
        table.records[0].primer.as_ref().unwrap().ordered_sequence(),
        "GGACGT"
    );
    let q = table.records[0].quantity.as_ref().unwrap();
    assert_eq!(q.dimension.as_deref(), Some("amount"));
    assert!((q.si_value.unwrap() - 25e-9).abs() < 1e-20);
    assert_eq!(table.records[1].name, "F");
    assert_ne!(table.records[0].row_id, table.records[1].row_id);
    assert!(table.records[1].primer.is_none());
    assert_eq!(
        table.records[1].modification.as_deref(),
        Some("5-prime phosphate")
    );
    assert_eq!(
        table.records[1]
            .quantity
            .as_ref()
            .unwrap()
            .dimension
            .as_deref(),
        Some("concentration")
    );
    assert!(
        table.records[2]
            .quantity
            .as_ref()
            .unwrap()
            .si_value
            .is_none()
    );
    assert!(table.records[2].primer.is_some());
    assert_eq!(
        STANDARD.decode(&table.original_base64).unwrap(),
        text.as_bytes()
    );
    assert!(table.selected_molecules(&[1]).is_err());
    let report = table.selected_molecules(&[0, 2]).unwrap();
    assert_eq!(report.molecules.len(), 2);
    let first = &report.molecules[0];
    assert_eq!(first.sequence.as_str(), "GGACGT");
    assert_eq!(first.features[0].location, Location::range(2, 6));
    assert_eq!(first.metadata["genbank.molecule_type"], "ss-DNA");
    assert!(collagene_core::chemistry::Duplex::from_molecule(first).is_err());
    assert!(!first.metadata.contains_key("oligo:original_table_base64"));
    assert_eq!(
        STANDARD
            .decode(&first.metadata["oligo:original_record_base64"])
            .unwrap(),
        b"F,GGACGT,ACGT,GG,,25,nmol,plain\r\n"
    );
    let native = collagene_core::protocol::save_document(first).unwrap();
    assert_eq!(
        collagene_core::protocol::load_document(&native)
            .unwrap()
            .revision(),
        first.revision()
    );
}
#[test]
fn csv_and_tsv_quotes_unicode_bom_and_multiline_extra_fields_retain_source() {
    let text = "\u{feff}\"Name\",\"Sequence\",\"Extra\"\r\n\"α, primer\",\"acgtn\",\"line one\r\nline \"\"two\"\"\"\r\n";
    let table = import_oligos(text, &OligoImportOptions::default()).unwrap();
    assert_eq!(table.records[0].name, "α, primer");
    assert_eq!(
        table.records[0].original_fields[2],
        "line one\r\nline \"two\""
    );
    assert_eq!(
        table.records[0].primer.as_ref().unwrap().annealing.as_str(),
        "ACGTN"
    );
    assert!(
        table.records[0]
            .diagnostics
            .iter()
            .any(|w| w.contains("full ordered"))
    );
    table.validate().unwrap();
    let opts = OligoImportOptions {
        delimiter: Delimiter::Tab,
        has_header: false,
        columns: Default::default(),
    };
    let table = import_oligos("x\tACGA\n\"y\tname\"\tTTGC\n", &opts).unwrap();
    assert_eq!(table.records.len(), 2);
    assert_eq!(table.records[1].name, "y\tname");
    assert!(
        table.records[1].primer.is_none(),
        "control in a mapped name must not become a DNA document"
    );
    assert_eq!(table.headers, vec!["Column 1", "Column 2"]);
}
#[test]
fn mapped_inconsistency_missing_cells_and_unknown_sequence_notation_never_guess_dna() {
    let table=import_oligos("Name,Order,Domain,Tail,Mod,Amount,Unit\na,GGACGT,ACGA,GG,,,\nb,GGACGT,,GG,,,\nc,GGACGT,ACGT\nd,[FAM]ACGT,ACGT,,,,\ne,ACGU,ACGU,,,,\n",&options()).unwrap();
    assert!(table.records.iter().all(|r| r.primer.is_none()));
    assert_eq!(
        table.records[3].raw_ordered_sequence.as_deref(),
        Some("[FAM]ACGT")
    );
    let mut bad = options();
    bad.columns.tail = Some(0);
    assert!(import_oligos("a,b\n", &bad).is_err());
    for text in [
        "name,sequence\na,\"ACGT\n",
        "name,sequence\na,AC\"GT\n",
        "name,sequence\na,\"ACGT\"oops\n",
        "name,sequence\ra,ACGT\r",
    ] {
        assert!(
            inspect_table(text, Delimiter::Comma, true).is_err(),
            "{text:?}"
        );
    }
}
#[test]
fn collection_validation_rejects_fabricated_computable_rows_or_source_evidence() {
    let text = "name,sequence\na,/5Phos/ACGT\n";
    let mut table = import_oligos(text, &Default::default()).unwrap();
    table.records[0].primer = Some(Primer::new("a", "ACGT", "").unwrap());
    assert!(table.validate().is_err());
    assert!(table.selected_molecules(&[0]).is_err());
    let mut table = import_oligos("name,sequence\na,ACGT\n", &Default::default()).unwrap();
    table.source_blake3 = "0".repeat(64);
    assert!(table.validate().is_err());
    let table = import_oligos("name,sequence\na,ACGT\n", &Default::default()).unwrap();
    assert!(table.selected_molecules(&[0, 0]).is_err());
    assert!(table.selected_molecules(&[1]).is_err());
}
#[test]
fn order_sheet_quoting_and_formula_policy_do_not_erase_original_values() {
    let table = import_oligos(
        "Name,Sequence\n\"=HYPERLINK(\"\"x\"\")\",ACGT\n\"  +SUM(1)\",TGCA\n",
        &Default::default(),
    )
    .unwrap();
    let safe = export_oligos(&table, Delimiter::Comma, TableExportMode::SpreadsheetSafe).unwrap();
    let rows = csv::Reader::from_reader(safe.text.as_bytes())
        .records()
        .map(Result::unwrap)
        .collect::<Vec<_>>();
    assert_eq!(&rows[0][0], "'=HYPERLINK(\"x\")");
    assert_eq!(&rows[1][0], "'  +SUM(1)");
    assert!(
        safe.warnings
            .iter()
            .any(|s| s.contains("intentionally differ"))
    );
    let machine = export_oligos(&table, Delimiter::Tab, TableExportMode::MachineData).unwrap();
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_reader(machine.text.as_bytes());
    assert_eq!(
        &reader.records().next().unwrap().unwrap()[0],
        "=HYPERLINK(\"x\")"
    );
    assert!(machine.warnings.iter().any(|s| s.contains("quoting alone")));
    assert_eq!(table.records[0].name, "=HYPERLINK(\"x\")");
}
#[test]
fn table_limits_and_positive_case_sensitive_units_are_explicit() {
    assert!(import_oligos(&"x".repeat(MAX_TABLE_BYTES + 1), &Default::default()).is_err());
    let many = format!("name,sequence\n{}", "x,ACGT\n".repeat(MAX_TABLE_ROWS + 1));
    assert!(import_oligos(&many, &Default::default()).is_err());
    let table = import_oligos(
        "N,O,D,T,M,A,U\na,ACGT,ACGT,,,25,nM\nb,ACGT,ACGT,,,25,NM\nc,ACGT,ACGT,,,-1,nmol\n",
        &options(),
    )
    .unwrap();
    assert_eq!(
        table.records[0]
            .quantity
            .as_ref()
            .unwrap()
            .dimension
            .as_deref(),
        Some("concentration")
    );
    assert!(
        table.records[1..]
            .iter()
            .all(|r| r.quantity.as_ref().unwrap().si_value.is_none())
    );
}
#[test]
fn fastq_wrapped_records_descriptions_quality_markers_and_duplicates_roundtrip() {
    let first = b"@read1 description\r\nacg\r\ntn\r\n+read1\r\n!@+\r\nI~\r\n";
    let second = b"@read1 other\nNN\n+\n@@\n";
    let bytes = [first.as_slice(), second.as_slice()].concat();
    let report = import_fastq(&bytes, QualityEncoding::Phred33).unwrap();
    assert_eq!(report.molecules.len(), 2);
    let m = &report.molecules[0];
    assert_eq!(m.name, "read1");
    assert_eq!(m.sequence.as_str(), "ACGTN");
    let evidence = quality_evidence(m).unwrap();
    assert_eq!(evidence.qualities, vec![0, 31, 10, 40, 93]);
    assert_eq!(evidence.header, "read1 description");
    assert!(evidence.current_valid);
    assert_eq!(original_fastq_bytes(m).unwrap(), first);
    assert_eq!(original_fastq_bytes(&report.molecules[1]).unwrap(), second);
    assert_eq!(
        export_fastq(m, QualityEncoding::Phred33).unwrap().text,
        "@read1 description\nACGTN\n+\n!@+I~\n"
    );
    assert_ne!(m.revision(), report.molecules[1].revision());
    let native = collagene_core::protocol::save_document(m).unwrap();
    let m = collagene_core::protocol::load_document(&native).unwrap();
    assert_eq!(
        quality_evidence(&m).unwrap().qualities,
        vec![0, 31, 10, 40, 93]
    );
}
#[test]
fn fastq_encoding_is_selected_never_inferred_and_conversion_is_lossless_or_rejected() {
    let source = b"@r\nAC\n+\n@h\n";
    let a = import_fastq(source, QualityEncoding::Phred33).unwrap();
    let b = import_fastq(source, QualityEncoding::Phred64).unwrap();
    assert_eq!(
        quality_evidence(&a.molecules[0]).unwrap().qualities,
        vec![31, 71]
    );
    assert_eq!(
        quality_evidence(&b.molecules[0]).unwrap().qualities,
        vec![0, 40]
    );
    assert_eq!(
        export_fastq(&b.molecules[0], QualityEncoding::Phred33)
            .unwrap()
            .text,
        "@r\nAC\n+\n!I\n"
    );
    assert!(export_fastq(&a.molecules[0], QualityEncoding::Phred64).is_err());
    assert!(import_fastq(b"@r\nA\n+\n!\n", QualityEncoding::Phred64).is_err());
}
#[test]
fn fastq_malformed_framing_scores_bases_and_header_mismatches_are_rejected() {
    for bytes in [
        b"@x\nAC\n+\n!\n".as_slice(),
        b"@x\nAC\n+\n!!!\n",
        b"@x\nAC\n+y\n!!\n",
        b"@x\nA C\n+\n!!!\n",
        b"@x\nAU\n+\n!!\n",
        b"@x\nAC\n+\n !\n",
        b"@x\nAC\n+\n!!\n\n",
        b"@x\nAC\n+\n!!\r",
    ] {
        assert!(
            import_fastq(bytes, QualityEncoding::Phred33).is_err(),
            "{bytes:?}"
        );
    }
}
#[test]
fn edited_calls_keep_original_quality_evidence_without_inventing_current_scores() {
    let source = b"@r\nACGT\n+\n!5I~\n";
    let m = import_fastq(source, QualityEncoding::Phred33)
        .unwrap()
        .molecules
        .remove(0);
    let edited = m.replace(1, 2, "A").unwrap();
    assert!(!quality_evidence(&edited).unwrap().current_valid);
    assert!(export_fastq(&edited, QualityEncoding::Phred33).is_err());
    assert_eq!(original_fastq_bytes(&edited).unwrap(), source);
    let mut tampered = m.clone();
    tampered
        .metadata
        .insert("fastq:quality_base64".into(), STANDARD.encode([1, 2, 3, 4]));
    assert!(quality_evidence(&tampered).is_err());
    let mut tampered = m.clone();
    tampered
        .metadata
        .insert("fastq:original_blake3".into(), "0".repeat(64));
    assert!(original_fastq_bytes(&tampered).is_err());
}
#[test]
fn palindromic_reverse_complement_reverses_scores_and_other_identical_coordinate_transforms_stale()
{
    let source = b"@palindrome\nACGT\n+\n!5I~\n";
    let m = import_fastq(source, QualityEncoding::Phred33)
        .unwrap()
        .molecules
        .remove(0);
    let mut rc = collagene_core::sequence::annotations::reverse_complement(&m).unwrap();
    project_fastq_reverse_complement(&m, &mut rc).unwrap();
    assert_eq!(rc.sequence, m.sequence);
    let e = quality_evidence(&rc).unwrap();
    assert!(e.current_valid);
    assert_eq!(e.qualities, vec![93, 40, 20, 0]);
    assert_eq!(
        export_fastq(&rc, QualityEncoding::Phred33).unwrap().text,
        "@palindrome\nACGT\n+\n~I5!\n"
    );
    let mut back = collagene_core::sequence::annotations::reverse_complement(&rc).unwrap();
    project_fastq_reverse_complement(&rc, &mut back).unwrap();
    assert_eq!(
        quality_evidence(&back).unwrap().qualities,
        vec![0, 20, 40, 93]
    );
    assert!(quality_evidence(&back).unwrap().current_valid);
    let extracted = collagene_core::sequence::annotations::extract_region(&m, 0, 4).unwrap();
    assert_eq!(extracted.sequence, m.sequence);
    assert!(!quality_evidence(&extracted).unwrap().current_valid);
    assert!(export_fastq(&extracted, QualityEncoding::Phred33).is_err());
}

#[test]
fn explicit_quality_trim_and_reversal_keep_actual_original_score_coordinates() {
    let raw = b"@palindromic_source desc\nAACGTT\n+\n!+5?I~\n";
    let source = import_fastq(raw, QualityEncoding::Phred33)
        .unwrap()
        .molecules
        .remove(0);
    let trimmed = trim_fastq(&source, 1, 4).unwrap();
    let e = quality_evidence(&trimmed).unwrap();
    assert!(e.current_valid);
    assert_eq!(e.called_sequence.as_str(), "ACG");
    assert_eq!(e.qualities, vec![10, 20, 30]);
    assert_eq!((e.retained_start, e.retained_end), (1, 4));
    let rc = collagene_core::sequence::annotations::reverse_complement(&trimmed).unwrap();
    let e = quality_evidence(&rc).unwrap();
    assert!(e.current_valid);
    assert_eq!(e.called_sequence.as_str(), "CGT");
    assert_eq!(e.qualities, vec![30, 20, 10]);
    assert_eq!((e.retained_start, e.retained_end), (1, 4));
    let shorter = trim_fastq(&rc, 1, 2).unwrap();
    let e = quality_evidence(&shorter).unwrap();
    assert_eq!(e.called_sequence.as_str(), "G");
    assert_eq!(e.qualities, vec![20]);
    let back = collagene_core::sequence::annotations::reverse_complement(&shorter).unwrap();
    let e = quality_evidence(&back).unwrap();
    assert!(e.current_valid);
    assert_eq!(e.called_sequence.as_str(), "C");
    assert_eq!(e.qualities, vec![20]);
    assert_eq!(original_fastq_bytes(&back).unwrap(), raw);
    let current = export_fastq(&back, QualityEncoding::Phred33).unwrap();
    assert!(current.text.contains("\nC\n+\n5\n"));
    assert!(trim_fastq(&source, 0, 0).is_err());
    assert!(trim_fastq(&source, 0, 7).is_err());
}

#[test]
fn independent_fastq_framing_matrix_keeps_at_and_plus_as_quality_characters() {
    let bases = b"ACGT";
    let scores = b"!@+~";
    for a in 0..64usize {
        for q in 0..64usize {
            let sequence = [bases[a / 16], bases[a / 4 % 4], bases[a % 4]];
            let quality = [scores[q / 16], scores[q / 4 % 4], scores[q % 4]];
            let raw = format!(
                "@matrix desc\n{}\n{}\n+matrix desc\n{}\n{}\n",
                sequence[0] as char,
                String::from_utf8(sequence[1..].to_vec()).unwrap(),
                quality[0] as char,
                String::from_utf8(quality[1..].to_vec()).unwrap()
            );
            let report = import_fastq(raw.as_bytes(), QualityEncoding::Phred33).unwrap();
            let e = quality_evidence(&report.molecules[0]).unwrap();
            assert_eq!(e.called_sequence.as_str().as_bytes(), sequence);
            assert_eq!(e.qualities, quality.map(|c| c - 33));
            let output = export_fastq(&report.molecules[0], QualityEncoding::Phred33).unwrap();
            assert_eq!(
                output.text,
                format!(
                    "@matrix desc\n{}\n+\n{}\n",
                    String::from_utf8(sequence.to_vec()).unwrap(),
                    String::from_utf8(quality.to_vec()).unwrap()
                )
            );
        }
    }
}

#[test]
fn independent_quoted_cell_matrix_preserves_both_delimiters_and_record_ranges() {
    let alphabet = ['a', ',', '"', '\n'];
    let mut cells = vec![String::new()];
    for length in 1..=3 {
        for number in 0..4usize.pow(length) {
            let mut n = number;
            let mut cell = String::new();
            for _ in 0..length {
                cell.push(alphabet[n % 4]);
                n /= 4;
            }
            cells.push(cell);
        }
    }
    for delimiter in [Delimiter::Comma, Delimiter::Tab] {
        let sep = if delimiter == Delimiter::Comma {
            ','
        } else {
            '\t'
        };
        for cell in &cells {
            let quoted = format!("\"{}\"", cell.replace('"', "\"\""));
            let raw = format!("Name{sep}Sequence{sep}Extra\r\n{quoted}{sep}ACGT{sep}{quoted}\r\n");
            let options = OligoImportOptions {
                delimiter,
                ..Default::default()
            };
            let table = import_oligos(&raw, &options).unwrap();
            assert_eq!(
                table.records[0].original_fields,
                vec![cell.clone(), "ACGT".into(), cell.clone()]
            );
            assert_eq!(
                STANDARD
                    .decode(&table.records[0].original_record_base64)
                    .unwrap(),
                format!("{quoted}{sep}ACGT{sep}{quoted}\r\n").as_bytes()
            );
        }
    }
}

#[test]
fn oligo_collection_json_roundtrip_preserves_mapped_quantities_and_validation() {
    // Independently specified SI expectations catch both conversion mistakes and
    // JSON parsers that round a serialized f64 to its adjacent representable value.
    let cases: &[(&str, &str, &str, f64)] = &[
        ("25", "nmol", "amount", 0.000000025),
        ("50", "nmol", "amount", 0.000000050),
        ("100", "pmol", "amount", 0.000000000100),
        ("1.25", "umol", "amount", 0.00000125),
        ("0.125", "mmol", "amount", 0.000125),
        ("2", "mol", "amount", 2.0),
        ("2.5", "nM", "concentration", 0.0000000025),
        ("7.25", "µM", "concentration", 0.00000725),
        ("0.2", "M", "concentration", 0.2),
        ("3", "mM", "concentration", 0.003),
        ("5", "pM", "concentration", 0.000000000005),
    ];
    let mut source = String::from("Name,Ordered,Domain,Tail,Modification,Amount,Unit\r\n");
    for (index, (value, unit, _, _)) in cases.iter().enumerate() {
        source.push_str(&format!("primer-{index},ACGT,ACGT,,,{value},{unit}\r\n"));
    }
    let table = import_oligos(&source, &options()).unwrap();
    table.validate().unwrap();
    let bytes = serde_json::to_vec(&table).unwrap();
    let decoded: OligoTable = serde_json::from_slice(&bytes).unwrap();
    for ((original, restored), (raw_value, raw_unit, dimension, expected)) in
        table.records.iter().zip(&decoded.records).zip(cases)
    {
        let quantity = original.quantity.as_ref().unwrap();
        let restored_quantity = restored.quantity.as_ref().unwrap();
        let actual = quantity.si_value.unwrap();
        assert_eq!(quantity.raw_value, *raw_value);
        assert_eq!(quantity.raw_unit, *raw_unit);
        assert_eq!(quantity.dimension.as_deref(), Some(*dimension));
        assert!(
            (actual - expected).abs() <= expected.abs() * f64::EPSILON * 2.0,
            "{raw_value} {raw_unit}: {actual} != independent SI expectation {expected}"
        );
        assert_eq!(
            actual.to_bits(),
            restored_quantity.si_value.unwrap().to_bits(),
            "{raw_value} {raw_unit}: JSON changed retained SI value"
        );
    }
    decoded.validate().unwrap();
    assert_eq!(decoded.original_bytes().unwrap(), source.as_bytes());
    assert_eq!(serde_json::to_vec(&decoded).unwrap(), bytes);
    assert_eq!(
        decoded.selected_molecules(&[0, 1]).unwrap().molecules.len(),
        2
    );
}
