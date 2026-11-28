use collagene_core::formats::annotations::{self as a, AnnotationFormat as F};
use collagene_model::{Feature, Location, Molecule, Topology};
fn reference() -> Molecule {
    Molecule::dna(
        "explicit actual DNA",
        "AAACCCGGGTTTAAACCCGGGTTT",
        Topology::Linear,
    )
    .unwrap()
}
fn q<'a>(f: &'a Feature, key: &str) -> Option<&'a str> {
    f.qualifiers
        .iter()
        .find(|(k, _)| k == key)
        .and_then(|(_, v)| v.as_deref())
}
#[test]
fn reverse_cds_and_repeated_attributes_remain_independent_rows_with_exact_original() {
    let r = reference();
    let text = "# GENCODE-shaped literal fixture\r\next\tHAVANA\tCDS\t4\t9\t0.25\t-\t2\tgene_id \"g\"; transcript_id \"t\"; gene_name \"semi; colon\"; tag \"basic\"; tag \"basic\"; tag \"other\"; exon_number 1;\r\next\tHAVANA\tCDS\t16\t21\t.\t-\t0\tgene_id \"g\"; transcript_id \"t\";\r\nother\t.\texon\t900\t999\t.\t+\t.\tgene_id \"off\"; transcript_id \"off\";\r\n";
    let report = a::import_annotations(&r, text, F::Gtf, "ext").unwrap();
    assert_eq!(report.proposals.len(), 2);
    assert_eq!(report.proposals[0].source_rows, vec![2]);
    assert_eq!(report.skipped_seqids, vec!["other"]);
    assert_eq!(report.original_bytes().unwrap(), text.as_bytes());
    let f = &report.proposals[0].feature;
    assert_eq!(f.location.extract(&r.sequence).unwrap(), "CCCGGG");
    assert_eq!(q(f, "codon_start"), Some("3"));
    assert!(f.needs_review);
    let tags = f
        .qualifiers
        .iter()
        .filter(|(k, _)| k == "gtf:attribute:tag")
        .map(|(_, v)| v.as_deref().unwrap())
        .collect::<Vec<_>>();
    assert_eq!(tags, vec!["basic", "basic", "other"]);
    let applied = a::apply_annotations(&r, &report, &[0, 1]).unwrap();
    let exported = a::export_annotations(&applied, F::Gtf, "renamed").unwrap();
    assert!(exported.text.contains("\t4\t9\t0.25\t-\t2\t"));
    assert!(
        exported
            .text
            .contains("tag \"basic\"; tag \"basic\"; tag \"other\";")
    );
    assert!(exported.text.contains("exon_number \"1\";"));
    let again = a::import_annotations(&r, &exported.text, F::Gtf, "renamed").unwrap();
    for (i, p) in again.proposals.iter().enumerate() {
        assert_eq!(p.feature.location, report.proposals[i].feature.location);
        assert_eq!(
            q(&p.feature, "codon_start"),
            q(&report.proposals[i].feature, "codon_start")
        );
    }
    assert_eq!(
        a::original_annotation_source(&applied, &report.source_blake3).unwrap(),
        text.as_bytes()
    );
}
#[test]
fn single_base_and_unknown_strand_or_missing_phase_are_not_guessed() {
    let r = reference();
    let text = "ext\t.\tCDS\t2\t2\t.\t+\t.\tgene_id \"g\"; transcript_id \"t\";\next\t.\texon\t4\t6\t.\t.\t.\tgene_id \"g\"; transcript_id \"t\";\n";
    let report = a::import_annotations(&r, text, F::Gtf, "ext").unwrap();
    assert_eq!(report.proposals[0].feature.location, Location::range(1, 2));
    assert!(q(&report.proposals[0].feature, "codon_start").is_none());
    assert!(
        report.proposals[1]
            .feature
            .location
            .extract(&r.sequence)
            .is_err()
    );
    let applied = a::apply_annotations(&r, &report, &[0, 1]).unwrap();
    assert!(
        a::export_annotations(&applied, F::Gtf, "ext")
            .unwrap()
            .text
            .contains("\t4\t6\t.\t.\t.\t")
    );
}
#[test]
fn invalid_syntax_unmapped_rows_and_stale_source_evidence_fail_atomically() {
    let r = reference();
    for attributes in [
        "gene_id=bad;",
        "gene_id \"unterminated;",
        "gene_id \"g\"",
        "gene_id \"g\";transcript_id \"t\";",
        "gene_id \"bad\\nvalue\";",
        "",
    ] {
        let text = format!("ext\t.\texon\t1\t2\t.\t+\t.\t{attributes}\n");
        assert!(
            a::import_annotations(&r, &text, F::Gtf, "ext").is_err(),
            "{attributes}"
        );
    }
    for (start, end, strand, phase) in [
        (0, 2, "+", "."),
        (4, 3, "+", "."),
        (24, 25, "+", "."),
        (1, 2, "?", "."),
        (1, 2, "+", "3"),
    ] {
        let text = format!(
            "ext\t.\tCDS\t{start}\t{end}\t.\t{strand}\t{phase}\tgene_id \"g\"; transcript_id \"t\";\n"
        );
        assert!(a::import_annotations(&r, &text, F::Gtf, "ext").is_err());
    }
    let valid = "ext\t.\texon\t1\t2\t.\t+\t.\tgene_id \"g\"; transcript_id \"t\";\n";
    let malformed = format!("{valid}other\t.\texon\t0\t2\t.\t+\t.\tgene_id \"off\";\n");
    assert!(a::import_annotations(&r, &malformed, F::Gtf, "ext").is_err());
    let report = a::import_annotations(&r, valid, F::Gtf, "ext").unwrap();
    let mut applied = a::apply_annotations(&r, &report, &[0]).unwrap();
    applied.features[0].location = Location::range(2, 4);
    assert!(a::export_annotations(&applied, F::Gtf, "ext").is_err());
}
#[test]
fn generic_joined_cds_export_uses_explicit_ids_and_reverse_biological_phase() {
    let mut r = reference();
    let qualifiers = vec![
        ("label".into(), Some("coding".into())),
        ("gtf:attribute:gene_id".into(), Some("g".into())),
        ("gtf:attribute:transcript_id".into(), Some("t".into())),
        ("codon_start".into(), Some("2".into())),
    ];
    r.features.push(Feature {
        kind: "CDS".into(),
        location: Location::Complement(Box::new(Location::Join(vec![
            Location::range(1, 5),
            Location::range(9, 14),
        ]))),
        qualifiers,
        needs_review: true,
    });
    let out = a::export_annotations(&r, F::Gtf, "ext").unwrap();
    let rows = out.text.lines().collect::<Vec<_>>();
    assert!(rows[0].contains("\t10\t14\t.\t-\t1\t"));
    assert!(rows[1].contains("\t2\t5\t.\t-\t2\t"));
    r.features[0]
        .qualifiers
        .retain(|(k, _)| k != "gtf:attribute:gene_id");
    assert!(a::export_annotations(&r, F::Gtf, "ext").is_err());
}
