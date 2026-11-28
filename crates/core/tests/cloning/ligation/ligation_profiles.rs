//! Published-source facts and an independent paper-level arithmetic oracle.
use collagene_core::{
    golden_gate::GoldenGateEnzyme, ligation_fidelity::UnknownReason, ligation_profiles as profiles,
};
use std::io::{Cursor, Read, Write};
fn reference(id: &str) -> profiles::ProfileRef {
    profiles::list()
        .unwrap()
        .into_iter()
        .find(|p| p.profile.profile_id == id)
        .unwrap()
        .profile
}
fn ends(values: &[&str]) -> Vec<String> {
    values.iter().map(|s| (*s).into()).collect()
}
fn alphabet(k: usize) -> Vec<String> {
    (0..4usize.pow(k as u32))
        .map(|mut n| {
            let mut text = vec![b'A'; k];
            for b in text.iter_mut().rev() {
                *b = b"ACGT"[n % 4];
                n /= 4;
            }
            String::from_utf8(text).unwrap()
        })
        .collect()
}
#[test]
fn five_profiles_pin_original_bytes_formulation_conditions_license_and_complete_matrices() {
    let expected = [
        (
            "pryor2020-bsai-hfv2",
            "BsaI-HFv2",
            202549,
            4,
            203364u64,
            789,
            61531,
        ),
        (
            "pryor2020-bsmbi-v2",
            "BsmBI-v2",
            203006,
            4,
            219378,
            920,
            61406,
        ),
        ("pryor2020-esp3i", "Esp3I", 202302, 4, 224222, 988, 61485),
        (
            "pryor2020-bbsi-hf",
            "BbsI-HF",
            203060,
            4,
            189282,
            712,
            61306,
        ),
        ("pryor2020-sapi", "SapI", 24232, 3, 172332, 3622, 3357),
    ];
    let all = profiles::list().unwrap();
    assert_eq!(all.len(), 5);
    assert_eq!(profiles::summaries().unwrap(), all);
    let catalog = profiles::catalog_revision().unwrap();
    assert_eq!(catalog.len(), 64);
    assert_eq!(profiles::catalog_revision().unwrap(), catalog);
    for (id, formulation, bytes, k, sum, maximum, zeros) in expected {
        let r = reference(id);
        let raw = profiles::source_bytes(&r).unwrap();
        assert_eq!(raw.len(), bytes);
        let p = profiles::load(&r).unwrap();
        let m = p.metadata();
        assert_eq!(m.enzyme_formulation, formulation);
        assert_eq!(m.overhang_length, k);
        assert_eq!(m.cells, 4usize.pow((2 * k) as u32));
        assert_eq!(
            m.provenance.source_blake3,
            blake3::hash(raw).to_hex().as_str()
        );
        assert!(m.provenance.license.starts_with("CC-BY-4.0"));
        assert!(
            m.provenance
                .citation
                .contains("10.1371/journal.pone.0238592")
        );
        assert!(
            m.provenance
                .normalization
                .contains("published ligation-frequency")
        );
        assert!(m.conditions.ligase.contains("T4"));
        assert!(m.conditions.temperature_protocol.contains(
            if k == 4 && formulation == "BsmBI-v2" {
                "42°C"
            } else {
                "37°C"
            }
        ));
        assert!(m.conditions.temperature_protocol.contains("30 cycles"));
        assert!(m.conditions.substrate_context.contains("100 nM"));
        if formulation == "BbsI-HF" {
            assert!(m.conditions.buffer.contains("CutSmart"));
            assert!(m.conditions.buffer.contains("DTT"));
        } else {
            assert!(m.conditions.buffer.contains("T4 DNA ligase"));
        }
        assert_eq!(m.data_revision.len(), 64);
        assert_eq!(m.model, profiles::MODEL);
        let names = alphabet(k);
        let mut total = 0;
        let mut largest = 0;
        let mut zero = 0;
        for a in &names {
            for b in &names {
                let value = p.weight(a, b).unwrap();
                assert_eq!(Some(value), p.weight(b, a));
                total += value;
                largest = largest.max(value);
                zero += usize::from(value == 0);
            }
        }
        assert_eq!(total, sum);
        assert_eq!(largest, maximum);
        assert_eq!(zero, zeros);
        assert!(
            std::sync::Arc::ptr_eq(&p, &profiles::load(&r).unwrap()),
            "parsed immutable profile is reused"
        );
    }
    let esp = profiles::load(&reference("pryor2020-esp3i")).unwrap();
    assert_eq!(esp.metadata().geometry_family, GoldenGateEnzyme::BsmBI);
    assert_eq!(esp.metadata().enzyme_formulation, "Esp3I");
}
#[test]
fn actual_axis_labels_place_watson_crick_counts_not_self_counts_on_the_visual_diagonal() {
    let r = reference("pryor2020-bsai-hfv2");
    let matrix = profiles::parse_published_matrix(profiles::source_bytes(&r).unwrap(), 4).unwrap();
    // Original B2 is635, but its row label isTTTT and its column isAAAA.
    assert_eq!(matrix.weight("TTTT", "AAAA"), Some(635));
    assert_eq!(matrix.weight("AAAA", "TTTT"), Some(635));
    assert_eq!(matrix.weight("TTTT", "TTTT"), Some(0));
    assert_eq!(matrix.weight("AAAA", "AAAA"), Some(0));
    assert_eq!(matrix.weight("TTTT", "AAAC"), Some(8));
    assert_eq!(matrix.weight("AAAC", "TTTT"), Some(8));
    assert_eq!(matrix.weight("ATAT", "ATAT"), Some(638));
    assert_eq!(matrix.cell_count(), 65536);
    let r = reference("pryor2020-sapi");
    let p = profiles::load(&r).unwrap();
    assert_eq!(p.weight("TTT", "AAA"), Some(1421));
    assert_eq!(p.weight("TTT", "TTT"), Some(0));
    assert_eq!(p.weight("TTTT", "AAAA"), None);
}
#[test]
fn independent_paper_fig4_set_reproduces_81_percent_and_92_percent_without_normalization() {
    let r = reference("pryor2020-bsmbi-v2");
    let all = ends(&[
        "GGAG", "TGAC", "TCCC", "TACT", "CCAT", "AATG", "AGCC", "TTCG", "GCTT", "GGTA", "CGCT",
    ]);
    let p = profiles::load(&r).unwrap();
    let prepared = p.prepare().unwrap();
    let report = prepared.score(&all).unwrap();
    assert!((report.overall_fidelity.unwrap() - 0.8092502583687903).abs() < 1e-14);
    assert!((report.log_overall_fidelity.unwrap() - 0.8092502583687903f64.ln()).abs() < 1e-14);
    report
        .validate_binding(&profiles::score_revision(&r, &all).unwrap())
        .unwrap();
    let reduced = all
        .iter()
        .filter(|s| s.as_str() != "GGTA")
        .cloned()
        .collect::<Vec<_>>();
    let report = prepared.score(&reduced).unwrap();
    assert!((report.overall_fidelity.unwrap() - 0.9211385107679805).abs() < 1e-14);
    assert_eq!(report.metadata.model, profiles::MODEL);
    assert_ne!(
        report.metadata.model,
        collagene_core::ligation_fidelity::MODEL
    );
}
#[test]
fn duplicate_reverse_complements_palindromes_serde_and_stale_receipts_remain_explicit() {
    let r = reference("pryor2020-bsai-hfv2");
    let e = ends(&["AAAA", "TTTT", "AAAA"]);
    let report = profiles::score(&r, &e).unwrap();
    assert_eq!(report.pairs.len(), 1);
    assert_eq!(report.pairs[0].occurrences, vec![0, 1, 2]);
    assert_eq!(report.overall_fidelity, Some(1.0));
    assert_eq!(
        report.pairs[0].estimate,
        profiles::Estimate::Known {
            correct_weight: 1270,
            total_weight: 1270,
            fidelity: 1.0
        }
    );
    let value = serde_json::to_value(&report).unwrap();
    assert_eq!(
        value["pairs"][0]["estimate"]["Known"]["correct_weight"],
        "1270"
    );
    let decoded: profiles::Report = serde_json::from_value(value).unwrap();
    decoded.validate_binding(&report.context_revision).unwrap();
    let decoded: profiles::Report =
        postcard::from_bytes(&postcard::to_allocvec(&report).unwrap()).unwrap();
    decoded.validate_binding(&report.context_revision).unwrap();
    let palindrome = profiles::score(&r, &ends(&["ATAT"])).unwrap();
    assert_eq!(palindrome.overall_fidelity, None);
    assert_eq!(
        palindrome.pairs[0].orientations[0].correct_weight,
        Some(638)
    );
    assert_eq!(
        palindrome.pairs[0].estimate,
        profiles::Estimate::Unknown {
            reason: UnknownReason::PalindromicEnd
        }
    );
    let mut stale = r.clone();
    stale.profile_revision = "0".repeat(64);
    assert!(profiles::load(&stale).is_err());
    assert!(profiles::revision(&stale, &e).is_err());
    let mut changed = report.clone();
    changed.metadata.enzyme_formulation = "BsaI".into();
    assert!(changed.validate_binding(&report.context_revision).is_err());
    let mut changed = report.clone();
    changed.metadata.provenance.license = "MIT".into();
    assert!(changed.validate_binding(&report.context_revision).is_err());
    let mut changed = report.clone();
    changed.pairs[0].orientations[0].total_weight = Some(1);
    assert!(changed.validate_binding(&report.context_revision).is_err());
    let mut changed = report.clone();
    changed.pairs[0].orientations[0].missing_partners = vec!["AAAA".into()];
    assert!(changed.validate_binding(&report.context_revision).is_err());
    assert!(profiles::score(&r, &ends(&["AAA"])).is_err());
    assert!(profiles::score(&r, &vec!["AAAA".into(); 257]).is_err());
}
fn rewrite(part: &str, edit: impl Fn(&str) -> String) -> Vec<u8> {
    let r = reference("pryor2020-sapi");
    let mut source =
        zip::ZipArchive::new(Cursor::new(profiles::source_bytes(&r).unwrap())).unwrap();
    let mut writer = zip::ZipWriter::new(Cursor::new(Vec::new()));
    for i in 0..source.len() {
        let mut f = source.by_index(i).unwrap();
        let name = f.name().to_owned();
        let mut raw = vec![];
        f.read_to_end(&mut raw).unwrap();
        if name == part {
            raw = edit(std::str::from_utf8(&raw).unwrap()).into_bytes();
        }
        writer
            .start_file(
                name,
                zip::write::SimpleFileOptions::default()
                    .compression_method(zip::CompressionMethod::Deflated),
            )
            .unwrap();
        writer.write_all(&raw).unwrap();
    }
    writer.finish().unwrap().into_inner()
}
fn modify_cell(xml: &str, coordinate: &str, action: impl Fn(&str) -> String) -> String {
    let start = xml.find(&format!("<c r=\"{coordinate}\"")).unwrap();
    let end = start + xml[start..].find("</c>").unwrap() + 4;
    format!(
        "{}{}{}",
        &xml[..start],
        action(&xml[start..end]),
        &xml[end..]
    )
}
#[test]
fn bounded_importer_rejects_formulas_missing_or_duplicate_cells_and_noninteger_weights() {
    const SHEET: &str = "xl/worksheets/sheet1.xml";
    for bytes in [
        rewrite(SHEET, |xml| {
            modify_cell(xml, "B2", |c| c.replace("<v>", "<f>1+1</f><v>"))
        }),
        rewrite(SHEET, |xml| modify_cell(xml, "B2", |_| String::new())),
        rewrite(SHEET, |xml| modify_cell(xml, "B2", |c| format!("{c}{c}"))),
        rewrite(SHEET, |xml| {
            modify_cell(xml, "B2", |c| c.replace("1421", "1421.5"))
        }),
        rewrite(SHEET, |xml| {
            modify_cell(xml, "B2", |c| c.replace("1421", "NaN"))
        }),
        rewrite(SHEET, |xml| {
            modify_cell(xml, "B2", |c| c.replace("1421", "1422"))
        }),
        rewrite(SHEET, |xml| {
            modify_cell(xml, "C1", |c| c.replace("<v>2</v>", "<v>1</v>"))
        }),
        rewrite(SHEET, |xml| {
            modify_cell(xml, "B2", |c| c.replace("B2", "XFD1048576"))
        }),
    ] {
        assert!(profiles::parse_published_matrix(&bytes, 3).is_err());
    }
}
#[test]
fn external_relationships_dtd_and_inflation_bombs_fail_without_external_runtime() {
    let external = rewrite("xl/_rels/workbook.xml.rels", |x| {
        x.replace(
            "Target=\"worksheets/sheet1.xml\"",
            "TargetMode=\"External\" Target=\"https://invalid.example/matrix\"",
        )
    });
    assert!(profiles::parse_published_matrix(&external, 3).is_err());
    let dtd = rewrite("xl/worksheets/sheet1.xml", |x| {
        x.replace(
            "<worksheet",
            "<!DOCTYPE worksheet [<!ENTITY bad 'external'>]><worksheet",
        )
    });
    assert!(profiles::parse_published_matrix(&dtd, 3).is_err());
    let bomb = rewrite("xl/sharedStrings.xml", |_| {
        format!("<sst>{}</sst>", "A".repeat((4 << 20) + 1))
    });
    assert!(bomb.len() < profiles::MAX_SOURCE_BYTES);
    assert!(profiles::parse_published_matrix(&bomb, 3).is_err());
    assert!(profiles::parse_published_matrix(&vec![0; profiles::MAX_SOURCE_BYTES + 1], 3).is_err());
    assert!(profiles::parse_published_matrix(b"not XLSX", 4).is_err());
}

#[test]
fn published_free_ends_are_actual_singleton_competitors_not_extra_successful_junctions() {
    let r = reference("pryor2020-bsai-hfv2");
    let intended = ends(&["AAAA"]);
    let extras = ends(&["AAAC"]);
    let report = profiles::score_with_competitors(&r, &intended, &extras).unwrap();
    assert_eq!(report.requested_ends, intended);
    assert_eq!(report.physical_competitors, extras);
    assert_eq!(report.physical_ends, ends(&["AAAA", "AAAC", "TTTT"]));
    assert!(!report.physical_ends.contains(&"GTTT".into()));
    // Literal S1 source cells: C(AAAA,TTTT)=C(TTTT,AAAA)=635;
    // C(TTTT,AAAC)=8 and C(AAAA,AAAC)=0. Invented GTTT adds3.
    assert_eq!(report.pairs.len(), 1);
    assert_eq!(
        report.pairs[0].estimate,
        profiles::Estimate::Known {
            correct_weight: 1270,
            total_weight: 1278,
            fidelity: 635.0 / 639.0
        }
    );
    assert_eq!(report.overall_fidelity, Some(635.0 / 639.0));
    report
        .validate_binding(&profiles::revision_with_competitors(&r, &intended, &extras).unwrap())
        .unwrap();
    let p = profiles::load(&r).unwrap();
    let palindrome = p
        .prepare()
        .unwrap()
        .score_with_competitors(&intended, &ends(&["AAAC", "ATAT"]))
        .unwrap();
    assert_eq!(palindrome.pairs.len(), 1);
    assert_eq!(
        palindrome.overall_fidelity, report.overall_fidelity,
        "unpaired ATAT never receives an intended-join multiplier"
    );
    let empty = profiles::score_with_competitors(&r, &[], &ends(&["AAAC"])).unwrap();
    assert!(empty.pairs.is_empty());
    assert_eq!(empty.physical_ends, ends(&["AAAC"]));
    assert_eq!(empty.overall_fidelity, None);
    assert_eq!(empty.log_overall_fidelity, None);
    empty
        .validate_binding(&profiles::revision_with_competitors(&r, &[], &ends(&["AAAC"])).unwrap())
        .unwrap();
    assert_eq!(
        serde_json::from_value::<profiles::Report>(serde_json::to_value(&report).unwrap()).unwrap(),
        report
    );
    let mut forged = report.clone();
    forged.physical_competitors = ends(&["GTTT"]);
    assert!(forged.validate_binding(&report.context_revision).is_err());
}
