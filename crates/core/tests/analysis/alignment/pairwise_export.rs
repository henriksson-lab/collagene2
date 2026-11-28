use collagene_core::{
    alignment::{
        AlignmentMode, AlignmentOptions, ComparisonOptions,
        export::{self as pairwise_export, Format},
        history::Record,
        review::{self as pairwise_review, Options},
    },
    documents::workspace::Workspace,
    formats::alignment::{self as alignment_formats, AlignmentFormat},
    material::lineage::EventBinding,
};
use collagene_model::{Molecule, Topology};
use std::collections::BTreeMap;

fn record(r: &Molecule, q: &Molecule, options: Options) -> Record {
    let frame = Workspace::from_documents("export", vec![r.clone(), q.clone()], 0).unwrap();
    let binding = |index: usize| {
        let d = &frame.documents[index];
        EventBinding {
            document_id: d.id.clone(),
            event_id: d.head_event.clone().unwrap(),
            revision: d.head.clone(),
        }
    };
    Record {
        design_consensus: None,
        schema_version: collagene_core::alignment::history::VERSION,
        reference: binding(0),
        query: binding(1),
        report: pairwise_review::compare(r, q, &options).unwrap(),
    }
}
fn dna(name: &str, calls: &str) -> Molecule {
    Molecule::dna(name, calls, Topology::Linear).unwrap()
}

#[test]
fn literal_insertion_rows_survive_exchange_and_rich_reports_keep_native_carets() {
    let r = dna("R <script>alert('x')</script>", "ACGTT");
    let q = dna("query", "ACGATT");
    let record = record(&r, &q, Options::Stored(AlignmentOptions::default()));
    assert_eq!(
        record.report.comparison.alignment.reference_aligned,
        "ACG-TT"
    );
    assert_eq!(record.report.comparison.alignment.query_aligned, "ACGATT");
    for (format, exchange) in [
        (Format::AlignedFasta, AlignmentFormat::AlignedFasta),
        (Format::Clustal, AlignmentFormat::Clustal),
    ] {
        let out = pairwise_export::export(&record, &r.name, &q.name, false, format).unwrap();
        let parsed = alignment_formats::import_alignment(
            &out.bytes,
            exchange,
            collagene_model::Alphabet::Dna,
            "computed output",
        )
        .unwrap();
        assert_eq!(parsed.rows[0].identifier, "reference");
        assert_eq!(parsed.rows[0].aligned, "ACG-TT");
        assert_eq!(parsed.rows[1].aligned, "ACGATT");
        assert!(
            out.warnings
                .iter()
                .any(|w| w.contains("omit native coordinate"))
        );
    }
    let out = pairwise_export::export(&record, &r.name, &q.name, false, Format::Html).unwrap();
    let text = String::from_utf8(out.bytes).unwrap();
    assert!(text.contains("Historical comparison - not recomputed on opening"));
    assert!(text.contains("R &lt;script&gt;alert(&#39;x&#39;)&lt;/script&gt;"));
    assert!(!text.contains("<script>"));
    assert!(text.contains("Native R [(3, 3)]; native Q [3, 4)"));
    assert!(text.contains(&record.reference.event_id));
    assert!(text.contains("ACG-TT"));
    assert!(text.contains("3=1I2="));
    assert!(text.contains("max_cells"));
}

#[test]
fn rotated_reverse_rows_are_not_silently_reoriented_to_native_sources() {
    let r = Molecule::dna("circle", "AGTCGATC", Topology::Circular).unwrap();
    let q = Molecule::dna("reverse circle", "ACTGATCG", Topology::Circular).unwrap();
    let record = record(&r, &q, Options::Whole(ComparisonOptions::default()));
    assert_eq!(record.report.comparison.reference_origin, 3);
    assert!(record.report.comparison.query_reversed);
    let out =
        pairwise_export::export(&record, &r.name, &q.name, true, Format::AlignedFasta).unwrap();
    let text = String::from_utf8(out.bytes).unwrap();
    assert_eq!(text.matches("CGATCAGT").count(), 2);
    assert!(!text.contains("\nAGTCGATC\n"));
    assert!(text.contains("reference_origin=3 query_reverse=true"));
    let rich = pairwise_export::export(&record, &r.name, &q.name, true, Format::Html).unwrap();
    let rich = String::from_utf8(rich.bytes).unwrap();
    assert!(rich.contains("Reference origin offset 3; query reverse complement true"));
    assert!(rich.contains("Exhaustive requested reference origins false"));
    assert!(rich.contains("best score proven true"));
}

#[test]
fn empty_local_result_exports_uncertainty_instead_of_making_up_exchange_rows() {
    let r = dna("all A", "AAAA");
    let q = dna("all C", "CCCC");
    let record = record(
        &r,
        &q,
        Options::Stored(AlignmentOptions {
            mode: AlignmentMode::Local,
            global_storage: collagene_core::alignment::grid::GlobalStorage::Full,
            ..Default::default()
        }),
    );
    assert!(
        record
            .report
            .comparison
            .alignment
            .reference_aligned
            .is_empty()
    );
    assert!(
        pairwise_export::export(&record, &r.name, &q.name, false, Format::AlignedFasta).is_err()
    );
    assert!(pairwise_export::export(&record, &r.name, &q.name, false, Format::Clustal).is_err());
    let out = pairwise_export::export(&record, &r.name, &q.name, false, Format::Html).unwrap();
    assert!(
        String::from_utf8(out.bytes)
            .unwrap()
            .contains("No aligned columns were returned")
    );
    let pdf = pairwise_export::export(&record, &r.name, &q.name, false, Format::Pdf).unwrap();
    assert!(pdf.bytes.starts_with(b"%PDF-"));
    let text = pdf_cells(&pdf.bytes)
        .into_iter()
        .map(|(_, _, text)| text)
        .collect::<String>();
    assert!(text.contains("SEQUENCE COMPARISON"));
    assert!(!text.contains("CLONING PROTOCOL"));
}

#[test]
fn embedded_pdf_keeps_source_identity_and_fixed_base_cell_positions() {
    let r = dna("reference", "ACGTT");
    let q = dna("query", "ACGATT");
    let record = record(&r, &q, Options::Stored(AlignmentOptions::default()));
    let out = pairwise_export::export(&record, &r.name, &q.name, false, Format::Pdf).unwrap();
    out.validate().unwrap();
    let bytes = String::from_utf8_lossy(&out.bytes);
    let cells = pdf_cells(&out.bytes);
    let text = cells
        .iter()
        .map(|(_, _, text)| text.as_str())
        .collect::<String>();
    assert!(text.contains("Historical comparison - not recomputed on opening"));
    assert!(text.contains(&record.reference.event_id));
    assert!(!text.contains("CLONING PROTOCOL"));
    // Inspect actual text operators and ToUnicode decoding, independently of
    // the font/renderer: each base and gap advances exactly six points.
    for (role, expected) in [("R", "R         0 ACG-TT"), ("Q", "Q         0 ACGATT")] {
        let y = cells
            .iter()
            .find(|(x, _, text)| *x == 42.0 && text == role)
            .unwrap()
            .1;
        let row = cells
            .iter()
            .filter(|(_, cy, _)| *cy == y)
            .collect::<Vec<_>>();
        assert_eq!(
            row.iter()
                .map(|(_, _, text)| text.as_str())
                .collect::<String>(),
            expected
        );
        for (index, (x, _, _)) in row.iter().enumerate() {
            assert_eq!(*x, 42.0 + 6.0 * index as f32);
        }
    }
    assert!(bytes.contains("/ToUnicode"));
    assert!(bytes.contains("/FontFile2"));
}

#[test]
fn source_and_output_shape_guards_reject_tampering_and_large_prints_fail_explicitly() {
    let r = dna("reference", "ACGTT");
    let q = dna("query", "ACGATT");
    let record = record(&r, &q, Options::Stored(AlignmentOptions::default()));
    let revision =
        pairwise_export::revision(&record, &r.name, &q.name, false, Format::Html).unwrap();
    assert_ne!(
        revision,
        pairwise_export::revision(&record, &r.name, &q.name, true, Format::Html).unwrap()
    );
    assert_ne!(
        revision,
        pairwise_export::revision(&record, &r.name, &q.name, false, Format::Pdf).unwrap()
    );
    let mut bad = record.clone();
    bad.query.revision = "0".repeat(64);
    assert!(pairwise_export::preflight(&bad, &r.name, &q.name).is_err());
    let mut bad = record.clone();
    bad.report.comparison.alignment.score += 1;
    assert!(pairwise_export::preflight(&bad, &r.name, &q.name).is_err());
    assert!(pairwise_export::preflight(&record, "x\n<script>", &q.name).is_err());
    let mut out = pairwise_export::export(&record, &r.name, &q.name, false, Format::Html).unwrap();
    out.extension = "pdf".into();
    assert!(out.validate().is_err());
    // Identity uses the native exact-match shortcut; no quadratic alignment is
    // needed to test print resource refusal rather than silent truncation.
    let long = dna("long exact reference", &"ACGT".repeat(40_000));
    let big = record_for_long(&long);
    assert!(pairwise_export::export(&big, &long.name, &long.name, false, Format::Pdf).is_err());
}
fn record_for_long(m: &Molecule) -> Record {
    record(m, m, Options::Whole(ComparisonOptions::default()))
}

fn literals(bytes: &[u8]) -> Vec<(Vec<u8>, usize)> {
    let mut result = vec![];
    let mut i = 0;
    while i < bytes.len() {
        if bytes[i] != b'(' {
            i += 1;
            continue;
        }
        let mut out = vec![];
        let mut depth = 1;
        i += 1;
        while i < bytes.len() && depth > 0 {
            let b = bytes[i];
            i += 1;
            match b {
                b'\\' => {
                    if i >= bytes.len() {
                        break;
                    }
                    let escaped = bytes[i];
                    i += 1;
                    if escaped.is_ascii_digit() && escaped <= b'7' {
                        let mut n = escaped - b'0';
                        for _ in 0..2 {
                            if i < bytes.len() && (b'0'..=b'7').contains(&bytes[i]) {
                                n = n.wrapping_mul(8).wrapping_add(bytes[i] - b'0');
                                i += 1;
                            } else {
                                break;
                            }
                        }
                        out.push(n);
                    } else {
                        match escaped {
                            b'n' => out.push(b'\n'),
                            b'r' => out.push(b'\r'),
                            b't' => out.push(b'\t'),
                            b'b' => out.push(8),
                            b'f' => out.push(12),
                            b'\n' => {}
                            b'\r' => {
                                if bytes.get(i) == Some(&b'\n') {
                                    i += 1;
                                }
                            }
                            x => out.push(x),
                        }
                    }
                }
                b'(' => {
                    depth += 1;
                    out.push(b);
                }
                b')' => {
                    depth -= 1;
                    if depth > 0 {
                        out.push(b);
                    }
                }
                _ => out.push(b),
            }
        }
        if depth == 0 {
            result.push((out, i));
        }
    }
    result
}

// Narrow independent reader for uncompressed content streams and their CMap.
fn pdf_cells(bytes: &[u8]) -> Vec<(f32, f32, String)> {
    let syntax = String::from_utf8_lossy(bytes);
    let mut map = BTreeMap::<u16, String>::new();
    for line in syntax.lines() {
        let fields = line.split_whitespace().collect::<Vec<_>>();
        if fields.len() != 2
            || fields[0].len() != 6
            || !fields[0].starts_with('<')
            || !fields[1].starts_with('<')
        {
            continue;
        }
        let Ok(cid) = u16::from_str_radix(&fields[0][1..5], 16) else {
            continue;
        };
        let unicode = &fields[1][1..fields[1].len() - 1];
        if unicode.len() % 4 != 0 {
            continue;
        }
        let units = (0..unicode.len())
            .step_by(4)
            .map(|i| u16::from_str_radix(&unicode[i..i + 4], 16).ok())
            .collect::<Option<Vec<_>>>();
        if let Some(Ok(text)) = units.map(|v| String::from_utf16(&v)) {
            map.insert(cid, text);
        }
    }
    let (mut start, mut result) = (0, Vec::new());
    while let Some(relative) = bytes[start..].windows(8).position(|s| s == b"\nstream\n") {
        let content = start + relative + 8;
        let dictionary = String::from_utf8_lossy(&bytes[start..content]);
        let Some(length) = dictionary
            .rsplit("/Length ")
            .next()
            .and_then(|s| s.split_whitespace().next())
            .and_then(|s| s.parse::<usize>().ok())
        else {
            break;
        };
        if content + length > bytes.len() {
            break;
        }
        let stream = &bytes[content..content + length];
        if stream.windows(3).any(|s| s == b"BT\n") {
            for block in stream
                .split(|b| *b == b'\n')
                .collect::<Vec<_>>()
                .split(|line| *line == b"ET")
            {
                let mut position = None;
                for line in block {
                    let text = String::from_utf8_lossy(line);
                    let fields = text.split_whitespace().collect::<Vec<_>>();
                    if fields.len() == 3 && fields[2] == "Td" {
                        position = fields[0]
                            .parse::<f32>()
                            .ok()
                            .zip(fields[1].parse::<f32>().ok());
                    }
                    if let Some((x, y)) = position {
                        for (literal, end) in literals(line) {
                            if !line[end..]
                                .iter()
                                .skip_while(|b| b.is_ascii_whitespace())
                                .copied()
                                .take(2)
                                .eq(b"Tj".iter().copied())
                            {
                                continue;
                            }
                            let text = literal
                                .chunks_exact(2)
                                .filter_map(|v| map.get(&u16::from_be_bytes([v[0], v[1]])))
                                .cloned()
                                .collect::<String>();
                            result.push((x, y, text));
                        }
                    }
                }
            }
        }
        start = content + length;
    }
    result
}
