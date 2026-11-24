use collagene_core::{
    Feature, Location, Molecule, Topology,
    analysis::Primer,
    pdf,
    verification::{self, DesignPlan, Operation, PlanInput, Source, Status},
};
use std::collections::BTreeMap;

fn molecule(name: &str, bases: &str, topology: Topology) -> Molecule {
    Molecule::dna(name, bases, topology).unwrap()
}
fn fixture() -> (DesignPlan, Molecule) {
    let left = molecule("left", "AAAACCCC", Topology::Linear);
    let right = molecule("right", "CCCCGGGG", Topology::Linear);
    let mut target = molecule("pβ_target", "AAAACCCCGGGG", Topology::Linear);
    target.features.push(Feature {
        kind: "CDS".into(),
        location: Location::Complement(Box::new(Location::Join(vec![
            Location::range(0, 4),
            Location::range(8, 12),
        ]))),
        qualifiers: vec![("label".into(), Some("β reverse compound".into()))],
        needs_review: false,
    });
    let plan = DesignPlan {
        name: "Assembly α → β / café".into(),
        target_revision: target.revision(),
        inputs: vec![
            PlanInput {
                id: "L".into(),
                revision: left.revision(),
                molecule: left,
            },
            PlanInput {
                id: "R".into(),
                revision: right.revision(),
                molecule: right,
            },
        ],
        steps: vec![Operation::Overlap {
            fragments: vec![Source::Input("L".into()), Source::Input("R".into())],
            overlaps: vec![4],
            circular: false,
        }],
    };
    (plan, target)
}
/// Independent narrow PDF reader for these uncompressed development outputs.
/// It decodes literal strings and the PDF ToUnicode mapping, not exporter internals.
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
fn extracted(bytes: &[u8]) -> String {
    let text = String::from_utf8_lossy(bytes);
    let mut map = BTreeMap::<u16, String>::new();
    for line in text.lines() {
        let tokens: Vec<_> = line.split_whitespace().collect();
        if tokens.len() != 2
            || !tokens[0].starts_with('<')
            || tokens[0].len() != 6
            || !tokens[1].starts_with('<')
        {
            continue;
        }
        let Ok(cid) = u16::from_str_radix(&tokens[0][1..5], 16) else {
            continue;
        };
        let unicode = &tokens[1][1..tokens[1].len() - 1];
        if unicode.len() % 4 != 0 {
            continue;
        }
        let units: Option<Vec<_>> = (0..unicode.len())
            .step_by(4)
            .map(|i| u16::from_str_radix(&unicode[i..i + 4], 16).ok())
            .collect();
        if let Some(units) = units {
            if let Ok(s) = String::from_utf16(&units) {
                map.insert(cid, s);
            }
        }
    }
    let mut out = String::new();
    // Delimit real streams by their declared /Length, independent of embedded bytes.
    let mut start = 0;
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
            for (s, end) in literals(stream) {
                if !stream[end..]
                    .iter()
                    .skip_while(|b| b.is_ascii_whitespace())
                    .copied()
                    .take(2)
                    .eq(b"Tj".iter().copied())
                {
                    continue;
                }
                for pair in s.chunks_exact(2) {
                    if let Some(t) = map.get(&u16::from_be_bytes([pair[0], pair[1]])) {
                        out.push_str(t);
                    }
                }
                out.push('\n');
            }
        }
        start = content + length;
    }
    out
}
fn check_xref(bytes: &[u8]) {
    let text = String::from_utf8_lossy(bytes);
    let start = text
        .rsplit("startxref\n")
        .next()
        .unwrap()
        .lines()
        .next()
        .unwrap()
        .parse::<usize>()
        .unwrap();
    assert_eq!(&bytes[start..start + 5], b"xref\n");
    let xref = std::str::from_utf8(&bytes[start..]).unwrap();
    let mut lines = xref.lines();
    assert_eq!(lines.next(), Some("xref"));
    let header = lines.next().unwrap().split_whitespace().collect::<Vec<_>>();
    assert_eq!(header[0], "0");
    let count = header[1].parse::<usize>().unwrap();
    for id in 0..count {
        let row = lines.next().unwrap();
        if id == 0 {
            assert_eq!(row.split_whitespace().last(), Some("f"));
            continue;
        }
        let offset = row[..10].parse::<usize>().unwrap();
        assert!(
            bytes[offset..].starts_with(format!("{id} 0 obj\n").as_bytes()),
            "xref {id} -> {offset}"
        );
    }
    assert!(bytes.ends_with(b"%%EOF") || bytes.ends_with(b"%%EOF\n"));
}
#[test]
fn verified_embedded_pdf_has_unicode_font_vector_map_and_independent_cross_reference() {
    let (plan, target) = fixture();
    let output = pdf::export_protocol_pdf(&plan, &target).unwrap();
    if std::env::var_os("COLLAGENE_PDF_FIXTURE").is_some() {
        std::fs::write("/tmp/collagene2-embedded-protocol.pdf", &output.bytes).unwrap();
    }
    assert_eq!(output.status, Status::ConsistentUnderAssumptions);
    assert!(output.pages >= 2);
    assert_eq!(output.plan_revision, verification::plan_revision(&plan));
    assert_eq!(output.target_revision, target.revision());
    assert_eq!(output.engine_version, verification::ENGINE_VERSION);
    assert!(output.warnings.is_empty(), "{:?}", output.warnings);
    assert!(output.bytes.starts_with(b"%PDF-1.7"));
    check_xref(&output.bytes);
    let syntax = String::from_utf8_lossy(&output.bytes);
    assert!(syntax.contains("/FontFile2"));
    assert!(syntax.contains("/ToUnicode"));
    assert!(syntax.contains("/CIDToGIDMap"));
    assert!(syntax.contains("/CMapType 2 def"));
    assert!(!syntax.contains("/Subtype /Image"));
    assert!(!syntax.contains("/JavaScript"));
    assert!(!syntax.contains("/URI"));
    let text = extracted(&output.bytes);
    assert!(text.contains("Assembly α → β / café"), "{text}");
    assert!(text.contains("β reverse compound"));
    assert!(text.contains("Junction overlaps [4] bp"));
    assert!(text.contains("Reconstructed product matches the separately bound target"));
    assert!(
        text.split_whitespace()
            .collect::<Vec<_>>()
            .join(" ")
            .contains("Reagent conditions and quantities are unspecified")
    );
    assert!(text.contains("Exhaustive alternative-product search: false"));
    assert_eq!(
        text.matches("COLLAGENE 2 / CLONING PROTOCOL").count(),
        output.pages
    );
    assert_eq!(
        text.matches("laboratory review required").count(),
        output.pages
    );
    assert!(text.contains(&format!("Page {} of {}", output.pages, output.pages)));
    if std::env::var_os("COLLAGENE_PDF_FIXTURE").is_some() {
        std::fs::write("/tmp/collagene2-embedded-protocol.pdf", &output.bytes).unwrap();
    }
}
#[test]
fn pdf_reverifies_current_binding_and_marks_every_invalid_page() {
    let (mut plan, target) = fixture();
    plan.steps[0] = Operation::Overlap {
        fragments: vec![Source::Input("L".into()), Source::Input("R".into())],
        overlaps: vec![3],
        circular: false,
    };
    let output = pdf::export_protocol_pdf(&plan, &target).unwrap();
    assert_eq!(output.status, Status::Invalid);
    let text = extracted(&output.bytes);
    assert_eq!(text.matches("DRAFT - INVALID").count(), output.pages);
    assert!(text.contains("Reconstructed product does not match"));
    assert!(!text.contains("Reconstructed product matches"));
    let (plan, mut target) = fixture();
    target.name.push_str(" changed");
    let output = pdf::export_protocol_pdf(&plan, &target).unwrap();
    assert_eq!(output.status, Status::Invalid);
    assert!(
        extracted(&output.bytes).contains("target revision changed or was not bound separately")
    );
}
#[test]
fn incomplete_drafts_retain_prior_products_and_scientific_limits() {
    let (mut plan, target) = fixture();
    plan.steps.push(Operation::Digest {
        template: Source::Step(0),
        enzymes: vec!["UnknownEnzyme".into()],
        fragment_index: 0,
    });
    let output = pdf::export_protocol_pdf(&plan, &target).unwrap();
    assert_eq!(output.status, Status::Incomplete);
    let text = extracted(&output.bytes);
    assert_eq!(text.matches("DRAFT - INCOMPLETE").count(), output.pages);
    assert!(text.contains("Reconstructed:"));
    assert!(text.contains("UnknownEnzyme"));
    assert!(text.contains("No reconstructed product for this step"));
}
#[test]
fn text_shaping_missing_glyphs_and_untrusted_pdf_syntax_are_explicit_data() {
    let (mut plan, target) = fixture();
    plan.name = "x ) Tj /JavaScript ( café β 中 🧬 ا a\u{301} \u{202e}".into();
    let output = pdf::export_protocol_pdf(&plan, &target).unwrap();
    assert_eq!(output.warnings.len(), 1);
    let text = extracted(&output.bytes);
    assert!(
        text.split_whitespace()
            .collect::<Vec<_>>()
            .join(" ")
            .contains("x ) Tj /JavaScript ( café β [U+4E2D] [U+1F9EC] [U+0627] a[U+0301] [U+202E]"),
        "{text}"
    );
    assert!(text.contains("PDF text rendering limits"));
    let syntax = String::from_utf8_lossy(&output.bytes);
    assert!(!syntax.contains("/S /JavaScript"));
    assert!(!syntax.contains("/OpenAction"));
    check_xref(&output.bytes);
}
#[test]
fn long_protocol_paginates_without_losing_final_steps_and_is_deterministic() {
    let (mut plan, target) = fixture();
    for i in 0..28 {
        plan.inputs[0]
            .molecule
            .metadata
            .insert(format!("irrelevant{i}"), "not printed".into());
    }
    // Reconstruct unsupported later steps as visibly unresolved; all requested steps still print.
    for i in 0..45 {
        plan.steps.push(Operation::TreatEnds {
            template: Source::Step(i),
            treatment: collagene_core::chemistry::EndTreatment::Phosphorylate,
        });
    }
    let output = pdf::export_protocol_pdf(&plan, &target).unwrap();
    assert!(output.pages > 4);
    let text = extracted(&output.bytes);
    assert!(text.contains("Step 46."));
    assert!(text.contains(&format!("Page {} of {}", output.pages, output.pages)));
    assert_eq!(
        output.bytes,
        pdf::export_protocol_pdf(&plan, &target).unwrap().bytes
    );
    check_xref(&output.bytes);
}
#[test]
fn physical_pcr_and_ligation_steps_print_oligos_tails_chemistry_nicks_and_source_provenance() {
    let template = molecule("template", "ACGACCCCGGGGTAAC", Topology::Linear);
    let target = molecule("product", "TTACGACCCCGGGGTAACGG", Topology::Linear);
    let mut input = template.clone();
    input
        .metadata
        .insert("source:provider".into(), "NCBI nuccore".into());
    input
        .metadata
        .insert("source:accession_version".into(), "SYN000001.1".into());
    let plan = DesignPlan {
        name: "PCR evidence".into(),
        target_revision: target.revision(),
        inputs: vec![PlanInput {
            id: "template".into(),
            revision: input.revision(),
            molecule: input,
        }],
        steps: vec![Operation::Pcr {
            template: Source::Input("template".into()),
            forward: Primer::new("F", "ACGA", "TT").unwrap(),
            reverse: Primer::new("R", "GTTA", "CC").unwrap(),
        }],
    };
    let output = pdf::export_protocol_pdf(&plan, &target).unwrap();
    let text = extracted(&output.bytes);
    assert!(text.contains("F: TTACGA (annealing ACGA; 5′ tail TT)"));
    assert!(text.contains("R: CCGTTA (annealing GTTA; 5′ tail CC)"));
    assert!(text.contains("source:accession_version: SYN000001.1"));
    assert!(text.contains("Top strand:"));
    assert!(text.contains("5′ terminus Hydroxyl"));
    assert!(text.contains("Retained nicks"));
}
#[test]
fn large_protocol_text_and_input_count_fail_instead_of_silently_truncating() {
    let (mut plan, target) = fixture();
    plan.name = "x".repeat(pdf::MAX_PDF_TEXT_BYTES + 1);
    assert!(matches!(
        pdf::export_protocol_pdf(&plan, &target),
        Err(collagene_core::Error::TooLarge)
    ));
    let (mut plan, target) = fixture();
    plan.name = "x\n".repeat(7000);
    assert!(matches!(
        pdf::export_protocol_pdf(&plan, &target),
        Err(collagene_core::Error::TooLarge)
    ));
    let (mut plan, target) = fixture();
    plan.inputs = vec![plan.inputs[0].clone(); 1001];
    assert!(matches!(
        pdf::export_protocol_pdf(&plan, &target),
        Err(collagene_core::Error::TooLarge)
    ));
}

#[test]
fn circular_map_keeps_origin_crossing_and_uncertain_review_geometry() {
    let source = molecule("closure", "AAAACCCCGGGGAAAA", Topology::Linear);
    let mut target = molecule("circular target", "AAAACCCCGGGG", Topology::Circular);
    target.features.push(Feature {
        kind: "CDS".into(),
        location: Location::Join(vec![Location::range(10, 12), Location::range(0, 4)]),
        qualifiers: vec![("label".into(), Some("across origin".into()))],
        needs_review: false,
    });
    target.features.push(Feature {
        kind: "misc_feature".into(),
        location: Location::Range {
            start: 4,
            end: 8,
            before: true,
            after: false,
        },
        qualifiers: vec![("label".into(), Some("uncertain feature".into()))],
        needs_review: true,
    });
    let plan = DesignPlan {
        name: "Circular closure".into(),
        target_revision: target.revision(),
        inputs: vec![PlanInput {
            id: "source".into(),
            revision: source.revision(),
            molecule: source,
        }],
        steps: vec![Operation::Overlap {
            fragments: vec![Source::Input("source".into())],
            overlaps: vec![4],
            circular: true,
        }],
    };
    let output = pdf::export_protocol_pdf(&plan, &target).unwrap();
    assert_eq!(output.status, Status::ConsistentUnderAssumptions);
    let text = extracted(&output.bytes);
    assert!(text.contains("Circular Dna"));
    assert!(text.contains("across origin"));
    assert!(text.contains("uncertain feature"));
    assert!(
        output
            .warnings
            .iter()
            .any(|w| w.contains("Approximate boundaries"))
    );
    let syntax = String::from_utf8_lossy(&output.bytes);
    assert!(syntax.contains("[3 2] 0 d"));
    assert!(syntax.contains(" c\n"));
    if std::env::var_os("COLLAGENE_PDF_FIXTURE").is_some() {
        std::fs::write("/tmp/collagene2-embedded-circle.pdf", &output.bytes).unwrap();
    }
}
