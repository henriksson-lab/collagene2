use collagene_core::{
    Alphabet, Feature, Location, Molecule, Sequence, Topology,
    annotation_formats::{self as a, AnnotationFormat as F},
};
fn reference() -> Molecule {
    Molecule::dna(
        "actual_reference_name",
        "AAACCCGGGTTTAAACCCGGGTTT",
        Topology::Linear,
    )
    .unwrap()
}
fn gff(body: &str) -> String {
    format!("##gff-version 3\n{body}")
}
fn q<'a>(feature: &'a Feature, key: &str) -> Option<&'a str> {
    feature
        .qualifiers
        .iter()
        .find(|(k, _)| k == key)
        .and_then(|(_, v)| v.as_deref())
}
fn apply(reference: &Molecule, text: &str, format: F) -> Molecule {
    let report = a::import_annotations(reference, text, format, "external_chr").unwrap();
    let selection = (0..report.proposals.len()).collect::<Vec<_>>();
    a::apply_annotations(reference, &report, &selection).unwrap()
}

#[test]
fn explicit_mapping_percent_attributes_parent_grouping_and_original_crlf_are_retained() {
    let r = reference();
    let text = "##gff-version 3\r\n##sequence-region external_chr 1 24\r\nexternal_chr\ttool%20A\tgene\t1\t21\t.\t+\t.\tID=g;Name=gene%3B%3D%2C%25+literal;Note=one%2Cvalue,two\r\nexternal_chr\ttool%20A\tCDS\t2\t5\t0.001\t+\t1\tID=c;Parent=g\r\nexternal_chr\ttool%20A\tCDS\t10\t14\t.\t+\t0\tID=c;Parent=g\r\nother\t.\tgene\t500\t600\t.\t+\t.\tID=other_gene\r\n";
    assert!(a::import_annotations(&r, text, F::Gff3, &r.name).is_err());
    let report = a::import_annotations(&r, text, F::Gff3, "external_chr").unwrap();
    assert_eq!(report.proposals.len(), 2);
    assert_eq!(report.skipped_seqids, vec!["other"]);
    assert_eq!(report.original_bytes().unwrap(), text.as_bytes());
    assert_eq!(report.proposals[0].feature.label(), "gene;=,%+literal");
    let notes = report.proposals[0]
        .feature
        .qualifiers
        .iter()
        .filter(|(k, _)| k == "gff3:attribute:Note")
        .map(|(_, v)| v.as_deref().unwrap())
        .collect::<Vec<_>>();
    assert_eq!(notes, vec!["one,value", "two"]);
    let cds = &report.proposals[1];
    assert_eq!(cds.source_rows, vec![4, 5]);
    assert_eq!(
        cds.feature.location,
        Location::Join(vec![Location::range(1, 5), Location::range(9, 14)])
    );
    assert_eq!(
        cds.feature.location.extract(&r.sequence).unwrap(),
        "AACCTTTAA"
    );
    assert_eq!(q(&cds.feature, "codon_start"), Some("2"));
    assert!(cds.feature.needs_review);
    report.validate(&r).unwrap();
    let restored: a::AnnotationImportReport =
        serde_json::from_slice(&serde_json::to_vec(&report).unwrap()).unwrap();
    restored.validate(&r).unwrap();
    let applied = a::apply_annotations(&r, &restored, &[0, 1]).unwrap();
    assert_eq!(applied.sequence, r.sequence);
    assert_eq!(applied.features.len(), 2);
    let native = collagene_core::protocol::save_document(&applied).unwrap();
    assert_eq!(
        collagene_core::protocol::load_document(&native).unwrap(),
        applied
    );
    assert_eq!(
        applied
            .metadata
            .iter()
            .filter(|(k, _)| k.ends_with(":original_base64"))
            .count(),
        1
    );
    let exported = a::export_annotations(&applied, F::Gff3, "new_chr").unwrap();
    assert!(exported.text.contains("Name=gene%3B%3D%2C%25+literal"));
    let again = a::import_annotations(&r, &exported.text, F::Gff3, "new_chr").unwrap();
    assert_eq!(again.proposals[1].feature.location, cds.feature.location);
    assert_eq!(q(&again.proposals[1].feature, "codon_start"), Some("2"));
}
#[test]
fn reverse_compound_cds_phase_is_counted_in_biological_order() {
    let r = reference();
    let text = gff(
        "external_chr\t.\tCDS\t2\t5\t.\t-\t0\tID=rev\nexternal_chr\t.\tCDS\t10\t14\t.\t-\t2\tID=rev\n",
    );
    let report = a::import_annotations(&r, &text, F::Gff3, "external_chr").unwrap();
    let f = &report.proposals[0].feature;
    assert_eq!(
        f.location,
        Location::Complement(Box::new(Location::Join(vec![
            Location::range(1, 5),
            Location::range(9, 14)
        ])))
    );
    assert_eq!(f.location.extract(&r.sequence).unwrap(), "TTAAAGGTT");
    assert_eq!(q(f, "codon_start"), Some("3"));
    let inconsistent = text.replace("\t-\t0\tID=rev", "\t-\t1\tID=rev");
    let report = a::import_annotations(&r, &inconsistent, F::Gff3, "external_chr").unwrap();
    assert_eq!(q(&report.proposals[0].feature, "codon_start"), None);
    assert_eq!(
        q(&report.proposals[0].feature, "gff3:inconsistent_phase"),
        Some("true")
    );
}
#[test]
fn single_circular_gff_row_crosses_origin_without_swapping_strand() {
    let mut r = reference();
    r.topology = Topology::Circular;
    for (strand, expected) in [("+", "TTAAAC"), ("-", "GTTTAA")] {
        let text = gff(&format!(
            "external_chr\t.\tpromoter\t23\t28\t.\t{strand}\t.\tID=origin\n"
        ));
        let report = a::import_annotations(&r, &text, F::Gff3, "external_chr").unwrap();
        let f = &report.proposals[0].feature;
        assert_eq!(f.location.extract(&r.sequence).unwrap(), expected);
        let applied = a::apply_annotations(&r, &report, &[0]).unwrap();
        let export = a::export_annotations(&applied, F::Gff3, "external_chr").unwrap();
        assert!(export.text.contains("\t23\t28\t"));
    }
    let text = gff("external_chr\t.\tgene\t23\t28\t.\t+\t.\tID=o\n");
    r.topology = Topology::Linear;
    assert!(a::import_annotations(&r, &text, F::Gff3, "external_chr").is_err());
}
#[test]
fn unknown_strands_and_zero_length_sites_do_not_invent_extractable_bases() {
    let r = reference();
    for strand in [".", "?"] {
        let text = gff(&format!(
            "external_chr\t.\tgene\t2\t5\t.\t{strand}\t.\tID=u\n"
        ));
        let report = a::import_annotations(&r, &text, F::Gff3, "external_chr").unwrap();
        let f = &report.proposals[0].feature;
        assert_eq!(f.location.spans(), vec![(1, 5, false)]);
        assert!(f.location.extract(&r.sequence).is_err());
    }
    for (text, format) in [
        (
            gff("external_chr\t.\tinsertion_site\t3\t3\t.\t+\t.\tID=p\n"),
            F::Gff3,
        ),
        ("external_chr\t3\t3\n".into(), F::Bed),
    ] {
        let report = a::import_annotations(&r, &text, format, "external_chr").unwrap();
        assert!(matches!(
            report.proposals[0].feature.location,
            Location::Unresolved(_)
        ));
        assert!(
            report.proposals[0]
                .feature
                .location
                .extract(&r.sequence)
                .is_err()
        );
    }
}
#[test]
fn multi_row_circular_anchor_and_remote_parts_are_explicit_not_flattened() {
    let mut r = reference();
    r.topology = Topology::Circular;
    let text = gff(
        "external_chr\t.\tCDS\t1\t3\t.\t+\t0\tID=c\nexternal_chr\t.\tCDS\t22\t24\t.\t+\t0\tID=c\n",
    );
    let report = a::import_annotations(&r, &text, F::Gff3, "external_chr").unwrap();
    assert!(matches!(
        report.proposals[0].feature.location,
        Location::Unresolved(_)
    ));
    assert_eq!(q(&report.proposals[0].feature, "codon_start"), None);
    let remote = text.replace(
        "external_chr\t.\tCDS\t22",
        "different_reference\t.\tCDS\t22",
    );
    assert!(a::import_annotations(&r, &remote, F::Gff3, "external_chr").is_err());
}
#[test]
fn bed12_blocks_strand_visualization_and_exact_original_are_preserved() {
    let r = reference();
    let text = "browser position external_chr:1-18\ntrack name=example itemRgb=On\nexternal_chr\t0\t18\tblocked reverse\t900\t-\t2\t17\t255,0,128\t2\t3,3,\t0,15,\n";
    let report = a::import_annotations(&r, text, F::Bed, "external_chr").unwrap();
    assert_eq!(report.directives.len(), 2);
    let f = &report.proposals[0].feature;
    assert_eq!(
        f.location,
        Location::Complement(Box::new(Location::Join(vec![
            Location::range(0, 3),
            Location::range(15, 18)
        ])))
    );
    assert_eq!(f.location.extract(&r.sequence).unwrap(), "GGGTTT");
    assert_eq!(f.kind, "misc_feature");
    assert_eq!(q(f, "codon_start"), None);
    let applied = a::apply_annotations(&r, &report, &[0]).unwrap();
    let export = a::export_annotations(&applied, F::Bed, "explicit_new_name").unwrap();
    assert_eq!(
        export.text,
        "explicit_new_name\t0\t18\tblocked reverse\t900\t-\t2\t17\t255,0,128\t2\t3,3,\t0,15,\n"
    );
    assert!(export.warnings.iter().any(|w| w.contains("attachments")));
    assert_eq!(report.original_bytes().unwrap(), text.as_bytes());
}
#[test]
fn bed_tracks_require_consistent_width_and_validate_all_rows_including_unmapped() {
    let r = reference();
    assert!(
        a::import_annotations(
            &r,
            "external_chr\t0\t3\nother\t0\t3\tlabel\t0\t+\n",
            F::Bed,
            "external_chr"
        )
        .is_err()
    );
    let text = "external_chr\t0\t3\ntrack name=next\nother\t0\t3\tlabel\t0\t+\n";
    let report = a::import_annotations(&r, text, F::Bed, "external_chr").unwrap();
    assert_eq!(report.skipped_seqids, vec!["other"]);
    assert!(
        a::import_annotations(
            &r,
            "external_chr\t0\t3\nother\t9\t3\n",
            F::Bed,
            "external_chr"
        )
        .is_err()
    );
    assert!(
        a::import_annotations(
            &r,
            "external_chr\t0\t3\tbad\t1001\t+\n",
            F::Bed,
            "external_chr"
        )
        .is_err()
    );
    let applied = apply(&r, "external_chr\t0\t3\n", F::Bed);
    assert!(applied.features[0].location.extract(&r.sequence).is_err());
    let export = a::export_annotations(&applied, F::Bed, "external_chr").unwrap();
    assert!(export.text.ends_with("\t0\t0\t0\t1\t3\t0\n"));
    assert!(export.warnings.iter().any(|w| w.contains("placeholders")));
}
#[test]
fn malformed_bed_block_geometry_and_optional_fields_fail_atomically() {
    let r = reference();
    let rows = [
        "external_chr\t0\t18\tx\t0\t-\t0\t18\t0\t2\t3,3\t1,15", // first block not0
        "external_chr\t0\t18\tx\t0\t-\t0\t18\t0\t2\t3,3\t0,14", // final notend
        "external_chr\t0\t18\tx\t0\t-\t0\t18\t0\t2\t10,9\t0,9", // overlap
        "external_chr\t0\t18\tx\t0\t-\t0\t18\t0\t3\t3,3\t0,15", // count
        "external_chr\t0\t18\tx\t0\t-\t0\t18\t0\t2\t3,0\t0,18", // zero block
        "external_chr\t0\t18\tx\t0\t-\t0\t18\t0\t2\t3,3\t0,9999999999999999999999999999999",
        "external_chr\t0\t18\tx\t0\t-\t0\t19",
        "external_chr\t0\t18\tx\t0\t-\t0\t18\t256,0,0",
        "external_chr\t0\t18\tx\t0\t?",
        "external_chr\t0\t25\tx\t0\t+",
    ];
    for text in rows {
        assert!(
            a::import_annotations(&r, text, F::Bed, "external_chr").is_err(),
            "{text}"
        );
    }
    assert!(r.features.is_empty());
}
#[test]
fn malformed_gff_attributes_phase_relationships_and_bounds_are_rejected() {
    let r = reference();
    let rows = [
        "external_chr\t.\tgene\t2\t5\t.\t+\t.\tID=x;Name=bad%ZZ\n",
        "external_chr\t.\tgene\t2\t5\t.\t+\t.\tID=x;Name=%FF\n",
        "external_chr\t.\tgene\t2\t5\t.\t+\t.\tID=x;ID=y\n",
        "external_chr\t.\tgene\t2\t5\t.\t+\t.\tID=x;Parent=missing\n",
        "external_chr\t.\tgene\t2\t5\t.\t+\t.\tID=x;Parent=y\nexternal_chr\t.\tgene\t7\t9\t.\t+\t.\tID=y;Parent=x\n",
        "external_chr\t.\tCDS\t2\t5\t.\t+\t.\tID=x\n",
        "external_chr\t.\tgene\t2\t5\t.\t+\t0\tID=x\n",
        "external_chr\t.\tgene\t0\t5\t.\t+\t.\tID=x\n",
        "external_chr\t.\tgene\t9\t5\t.\t+\t.\tID=x\n",
        "external_chr\t.\tgene\t2\t25\t.\t+\t.\tID=x\n",
        "external_chr\t.\tgene\t2\t5\tNaN\t+\t.\tID=x\n",
        "external_chr\t.\tgene\t2\t5\t.\t+\t.\tID=x\nexternal_chr\t.\tgene\t4\t8\t.\t+\t.\tID=x\n",
    ];
    for body in rows {
        assert!(
            a::import_annotations(&r, &gff(body), F::Gff3, "external_chr").is_err(),
            "{body}"
        );
    }
    assert!(
        a::import_annotations(
            &r,
            "external_chr\t.\tgene\t2\t5\t.\t+\t.\tID=x\n",
            F::Gff3,
            "external_chr"
        )
        .is_err()
    );
    assert!(a::import_annotations(&r,&gff("external_chr\t.\tgene\t2\t5\t.\t+\t.\tID=x\n###\nexternal_chr\t.\tgene\t7\t9\t.\t+\t.\tID=x\n"),F::Gff3,"external_chr").is_err());
}
#[test]
fn embedded_fasta_is_checked_against_explicit_reference_not_silently_ignored() {
    let r = reference();
    let text = gff(&format!(
        "external_chr\t.\tgene\t2\t5\t.\t+\t.\tID=x\n##FASTA\n>external_chr\n{}\n",
        r.sequence.as_str()
    ));
    a::import_annotations(&r, &text, F::Gff3, "external_chr").unwrap();
    let wrong = text.replace(r.sequence.as_str(), "TTTCCCGGGTTTAAACCCGGGTTT");
    assert!(a::import_annotations(&r, &wrong, F::Gff3, "external_chr").is_err());
}
#[test]
fn proposals_are_source_revision_bound_and_selection_and_tampering_are_atomic() {
    let r = reference();
    let text = gff("external_chr\t.\tpromoter\t2\t5\t.\t+\t.\tID=x\n");
    let report = a::import_annotations(&r, &text, F::Gff3, "external_chr").unwrap();
    let changed = r.replace(0, 1, "T").unwrap();
    assert!(a::apply_annotations(&changed, &report, &[0]).is_err());
    assert!(a::apply_annotations(&r, &report, &[0, 0]).is_err());
    assert!(a::apply_annotations(&r, &report, &[1]).is_err());
    assert!(a::apply_annotations(&r, &report, &[]).is_err());
    let mut tampered = report.clone();
    tampered.proposals[0].feature.location = Location::range(3, 5);
    assert!(tampered.validate(&r).is_err());
    let mut tampered = report.clone();
    tampered.warnings.clear();
    assert!(tampered.validate(&r).is_err());
    let mut tampered = report.clone();
    tampered.source_blake3 = "0".repeat(64);
    assert!(tampered.original_bytes().is_err());
    assert!(r.features.is_empty());
    let applied = a::apply_annotations(&r, &report, &[0]).unwrap();
    let rebound = a::import_annotations(&applied, &text, F::Gff3, "external_chr").unwrap();
    assert!(a::apply_annotations(&applied, &rebound, &[0]).is_err());
}
#[test]
fn source_record_export_guard_blocks_stale_edit_reverse_origin_and_forged_geometry() {
    let mut r = reference();
    r.topology = Topology::Circular;
    let annotated = apply(
        &r,
        &gff("external_chr\t.\tpromoter\t5\t9\t.\t+\t.\tID=x\n"),
        F::Gff3,
    );
    let shifted = annotated.replace(0, 0, "AAA").unwrap();
    assert!(a::export_annotations(&shifted, F::Gff3, "external_chr").is_err());
    let reverse = collagene_core::annotations::reverse_complement(&annotated).unwrap();
    assert!(a::export_annotations(&reverse, F::Gff3, "external_chr").is_err());
    let rotated = collagene_core::annotations::rotate_origin(&annotated, 4).unwrap();
    assert!(a::export_annotations(&rotated, F::Gff3, "external_chr").is_err());
    let mut forged = annotated.clone();
    let raw = forged.features[0]
        .qualifiers
        .iter_mut()
        .find(|(k, _)| k == "annotation_exchange:source_record")
        .unwrap()
        .1
        .as_mut()
        .unwrap();
    let mut value: serde_json::Value = serde_json::from_str(raw).unwrap();
    value["gff"][0]["start"] = serde_json::json!(6);
    *raw = serde_json::to_string(&value).unwrap();
    assert!(a::export_annotations(&forged, F::Gff3, "external_chr").is_err());
    let mut renamed = annotated.clone();
    renamed.features[0].qualifiers[0].1 = Some("changed label".into());
    assert!(a::export_annotations(&renamed, F::Gff3, "external_chr").is_err());
}
#[test]
fn selected_children_cannot_export_fabricated_absent_parent_hierarchies() {
    let r = reference();
    let text = gff(
        "external_chr\t.\tgene\t1\t20\t.\t+\t.\tID=g\nexternal_chr\t.\tCDS\t2\t5\t.\t+\t0\tID=c;Parent=g\n",
    );
    let report = a::import_annotations(&r, &text, F::Gff3, "external_chr").unwrap();
    let child = a::apply_annotations(&r, &report, &[1]).unwrap();
    assert!(a::export_annotations(&child, F::Gff3, "external_chr").is_err());
    assert_eq!(report.original_bytes().unwrap(), text.as_bytes());
}
#[test]
fn generic_export_requires_explicit_cds_phase_and_does_not_flatten_origin_to_bed() {
    let mut r = reference();
    r.features.push(Feature {
        kind: "CDS".into(),
        location: Location::Complement(Box::new(Location::Join(vec![
            Location::range(1, 5),
            Location::range(9, 14),
        ]))),
        qualifiers: vec![("label".into(), Some("independent CDS".into()))],
        needs_review: true,
    });
    assert!(a::export_annotations(&r, F::Gff3, "explicit").is_err());
    r.features[0]
        .qualifiers
        .push(("codon_start".into(), Some("3".into())));
    let report = a::export_annotations(&r, F::Gff3, "explicit").unwrap();
    let phases = report
        .text
        .lines()
        .filter(|l| !l.starts_with('#'))
        .map(|l| l.split('\t').nth(7).unwrap())
        .collect::<Vec<_>>();
    assert_eq!(phases, vec!["2", "0"]);
    let bed = a::export_annotations(&r, F::Bed, "explicit").unwrap();
    assert!(bed.text.contains("\t-\t1\t1\t0\t2\t4,5\t0,8\n"));
    r.topology = Topology::Circular;
    r.features[0].kind = "promoter".into();
    r.features[0].location = Location::Join(vec![Location::range(21, 24), Location::range(0, 3)]);
    assert!(a::export_annotations(&r, F::Bed, "explicit").is_err());
    let gff = a::export_annotations(&r, F::Gff3, "explicit").unwrap();
    assert!(gff.text.contains("\t22\t27\t"));
}
#[test]
fn literal_small_coordinate_matrix_covers_both_format_boundaries_and_strands() {
    // Independent substring/letter complement oracle; no implementation coordinate helpers.
    let bases = "AAACCCGGGTTTAAACCCGGGTTT";
    let r = reference();
    let mut cases = 0;
    for start in 0..bases.len() {
        for end in start + 2..=bases.len() {
            for reverse in [false, true] {
                let expected = if reverse {
                    bases.as_bytes()[start..end]
                        .iter()
                        .rev()
                        .map(|b| match b {
                            b'A' => 'T',
                            b'C' => 'G',
                            b'G' => 'C',
                            b'T' => 'A',
                            _ => unreachable!(),
                        })
                        .collect::<String>()
                } else {
                    bases[start..end].into()
                };
                let strand = if reverse { "-" } else { "+" };
                let bed = format!("external_chr\t{start}\t{end}\tx\t0\t{strand}\n");
                let gff = gff(&format!(
                    "external_chr\t.\tgene\t{}\t{end}\t.\t{strand}\t.\tID=x\n",
                    start + 1
                ));
                for (text, format) in [(&bed, F::Bed), (&gff, F::Gff3)] {
                    let report = a::import_annotations(&r, text, format, "external_chr").unwrap();
                    assert_eq!(
                        report.proposals[0]
                            .feature
                            .location
                            .extract(&r.sequence)
                            .unwrap(),
                        expected
                    );
                    cases += 1;
                }
            }
        }
    }
    assert_eq!(cases, 1104);
}
#[test]
fn missing_mapping_protein_controls_bounds_and_large_sources_are_explicit_errors() {
    let r = reference();
    for seqid in ["", ".", "external\nchr"] {
        assert!(a::import_annotations(&r, "external_chr\t0\t3\n", F::Bed, seqid).is_err());
    }
    let mut protein = r.clone();
    protein.sequence = Sequence::new(Alphabet::Protein, "MKKW").unwrap();
    assert!(
        a::import_annotations(&protein, "external_chr\t0\t3\n", F::Bed, "external_chr").is_err()
    );
    assert!(a::import_annotations(&r, "external_chr\t0\t3\0\n", F::Bed, "external_chr").is_err());
    assert!(
        a::import_annotations(
            &r,
            &"x".repeat(a::MAX_ANNOTATION_BYTES + 1),
            F::Bed,
            "external_chr"
        )
        .is_err()
    );
    assert!(
        a::import_annotations(
            &r,
            "external_chr\t0\t999999999999999999999999999999999999\n",
            F::Bed,
            "external_chr"
        )
        .is_err()
    );
}

#[test]
fn generic_types_do_not_invent_ontology_equivalence_and_directives_are_bounded() {
    let mut r = reference();
    r.features.push(Feature {
        kind: "laboratory custom label".into(),
        location: Location::range(1, 5),
        qualifiers: vec![("label".into(), Some("display".into()))],
        needs_review: true,
    });
    let output = a::export_annotations(&r, F::Gff3, "external_chr").unwrap();
    assert!(output.text.contains("\tsequence_feature\t2\t5\t"));
    assert!(
        output
            .text
            .contains("cg2_original_type=laboratory custom label")
    );
    assert!(
        output
            .warnings
            .iter()
            .any(|w| w.contains("ontology equivalence"))
    );
    r.features[0].location = Location::range(1, 2);
    let single = a::export_annotations(&r, F::Gff3, "external_chr").unwrap();
    assert!(single.text.contains("\tsequence_feature\t2\t2\t"));
    let roundtrip =
        a::import_annotations(&reference(), &single.text, F::Gff3, "external_chr").unwrap();
    assert_eq!(
        roundtrip.proposals[0].feature.location.spans(),
        vec![(1, 2, false)]
    );
    assert!(a::export_annotations(&r, F::Bed, "external_chr").is_ok());
    let source = format!("{}external_chr\t0\t3\n", "# comment\n".repeat(20_001));
    assert!(a::import_annotations(&reference(), &source, F::Bed, "external_chr").is_err());
    // Repeated completed-boundary directives do not rescan every prior row.
    let source = format!(
        "##gff-version 3\nexternal_chr\t.\tgene\t2\t5\t.\t+\t.\tID=x\n{}",
        "###\n".repeat(10_000)
    );
    assert!(a::import_annotations(&reference(), &source, F::Gff3, "external_chr").is_ok());
}

#[test]
fn gff_single_base_intervals_and_insertion_boundaries_are_distinct() {
    use collagene_core::annotation_formats::{
        AnnotationFormat, apply_annotations, export_annotations, import_annotations,
    };
    let reference = Molecule::dna("onebase", "ACGT", Topology::Linear).unwrap();
    for kind in ["exon", "SNP", "SNV", "gene", "CDS"] {
        let phase = if kind == "CDS" { "0" } else { "." };
        let source =
            format!("##gff-version 3\nonebase\ttest\t{kind}\t2\t2\t.\t-\t{phase}\tID=one\n");
        let report =
            import_annotations(&reference, &source, AnnotationFormat::Gff3, "onebase").unwrap();
        assert_eq!(
            report.proposals[0]
                .feature
                .location
                .extract(&reference.sequence)
                .unwrap(),
            "G"
        );
        assert_eq!(
            report.proposals[0].feature.location.spans(),
            vec![(1, 2, true)]
        );
        let annotated = apply_annotations(&reference, &report, &[0]).unwrap();
        let output = export_annotations(&annotated, AnnotationFormat::Gff3, "onebase").unwrap();
        let reopened =
            import_annotations(&reference, &output.text, AnnotationFormat::Gff3, "onebase")
                .unwrap();
        assert_eq!(
            reopened.proposals[0].feature.location,
            report.proposals[0].feature.location
        );
    }
    let boundary = import_annotations(
        &reference,
        "##gff-version 3\nonebase\ttest\tinsertion_site\t2\t2\t.\t+\t.\tID=site\n",
        AnnotationFormat::Gff3,
        "onebase",
    )
    .unwrap();
    assert!(matches!(
        boundary.proposals[0].feature.location,
        Location::Unresolved(_)
    ));
    let mut generic = reference.clone();
    generic.features.push(Feature {
        kind: "exon".into(),
        location: Location::range(1, 2),
        qualifiers: vec![("label".into(), Some("one base".into()))],
        needs_review: false,
    });
    let text = export_annotations(&generic, AnnotationFormat::Gff3, "onebase").unwrap();
    let restored =
        import_annotations(&reference, &text.text, AnnotationFormat::Gff3, "onebase").unwrap();
    assert_eq!(
        restored.proposals[0]
            .feature
            .location
            .extract(&reference.sequence)
            .unwrap(),
        "C"
    );
}
