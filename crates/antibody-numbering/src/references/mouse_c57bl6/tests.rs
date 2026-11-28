//! Maintainer raw-pair curation checks, separate from compact runtime loading.
use super::*;

fn sources() -> (input::RawSet, Vec<input::GappedRecord>) {
    input::decode(JSON, GAPPED).unwrap()
}

#[test]
fn pinned_capture_and_independent_inventory_are_not_numbering_model_claims() {
    assert_eq!(JSON.len(), 254_707);
    assert_eq!(GAPPED.len(), 35_235);
    let prepared = Prepared::new(JSON, GAPPED).unwrap();
    let d = prepared.dataset();
    assert_eq!(d.source.set_id, "OGRDB:G00004");
    assert_eq!(d.source.release_version, 5);
    assert_eq!(d.source.taxonomy_id, "NCBITAXON:10090");
    assert_eq!(d.species_subgroup, "C57BL/6");
    assert_eq!(
        (d.joins.len(), d.v_entries.len(), d.j_segments.len()),
        (114, 95, 4)
    );
    assert_eq!(d.source.excluded.len(), 15);
    assert_eq!(
        d.source
            .excluded
            .iter()
            .filter(|a| a.sequence_type == "D")
            .count(),
        8
    );
    assert_eq!(
        d.source
            .excluded
            .iter()
            .filter(|a| a.reason == "provider-not-functional")
            .count(),
        7
    );
    assert_eq!(d.license, "CC0-1.0");
    assert!(d.v_label_evidence.contains("JSON labels are absent"));
    assert_eq!(prepared.identity(), d.identity().unwrap());
    for c in &d.captures {
        assert_eq!(c.blake3.len(), 64);
        assert!(c.url.ends_with("/5") || c.url.ends_with("/5/gapped"));
    }
}

#[test]
fn literal_provider_gapped_columns_match_source_cdrs_without_repairing_atypical_anchors() {
    let prepared = Prepared::new(JSON, GAPPED).unwrap();
    let d = prepared.dataset();
    let v = d.v_entries.iter().find(|v| v.label == "IGHV-24BS").unwrap();
    let join = d
        .joins
        .iter()
        .find(|j| j.allele_ref == v.allele_ref)
        .unwrap();
    assert_eq!(join.gapped_label, "IGHV0-24BS*00");
    assert_eq!(
        join.allele_ref,
        "OGRDB:Mus musculus/C57BL/6_IGH:IGHV-24BS.1"
    );
    for (column, codon, residue) in [(23, "TGC", b'C'), (41, "TGG", b'W'), (104, "TGT", b'C')] {
        assert_eq!(&v.aligned_nucleotides[(column - 1) * 3..column * 3], codon);
        assert_eq!(v.aligned_protein.as_bytes()[column - 1], residue);
    }
    assert_eq!(d.atypical_anchors.len(), 4);
    let observed = d
        .atypical_anchors
        .iter()
        .map(|a| {
            let label = &d
                .v_entries
                .iter()
                .find(|v| v.allele_ref == a.allele_ref)
                .unwrap()
                .label;
            (label.as_str(), a.imgt_column, a.codon.as_str(), a.residue)
        })
        .collect::<Vec<_>>();
    assert!(observed.contains(&("IGHV-ZPQB", 23, "AGC", 'S')));
    assert!(observed.contains(&("IGHV-IQE5", 41, "TGT", 'C')));
    assert!(observed.contains(&("IGHV-OSZZ", 41, "AGG", 'R')));
    assert!(observed.contains(&("IGHV-OSZZ", 104, "TCT", 'S')));
    assert_eq!(
        d.atypical_anchors
            .iter()
            .map(|a| &a.allele_ref)
            .collect::<BTreeSet<_>>()
            .len(),
        3
    );
    assert_eq!(
        d.v_entries
            .iter()
            .filter(|v| !d
                .atypical_anchors
                .iter()
                .any(|a| a.allele_ref == v.allele_ref))
            .count(),
        92
    );
}

#[test]
fn j_facts_preserve_real_provider_frame_disagreement_and_do_not_supply_numbering_labels() {
    let prepared = Prepared::new(JSON, GAPPED).unwrap();
    for j in &prepared.dataset().j_segments {
        assert_eq!(j.fr4_protein.len(), 11);
        assert!(j.fr4_protein.starts_with("WG"));
        assert_eq!(j.fr4_protein.as_bytes()[3], b'G');
        assert_eq!(j.trailing_nucleotides, "G");
        assert_eq!(j.anchor_offset + 1, j.provider_j_cdr3_end);
        assert_eq!(j.anchor_phase, j.provider_codon_frame - 1);
        assert_eq!(j.literal_airr_frame_phase, (4 - j.provider_codon_frame) % 3);
        assert_ne!(j.anchor_phase, j.literal_airr_frame_phase);
        let json = serde_json::to_value(j).unwrap();
        assert!(json.get("fr4_positions").is_none());
    }
    let j = prepared
        .dataset()
        .j_segments
        .iter()
        .find(|j| j.label == "IGHJ-32C2")
        .unwrap();
    assert_eq!(
        (
            j.provider_codon_frame,
            j.provider_j_cdr3_end,
            j.anchor_offset
        ),
        (2, 20, 19)
    );
    assert_eq!(&j.fr4_nucleotides[..3], "TGG");
}

#[test]
fn raw_integrity_bounds_and_line_endings_are_checked_before_curation() {
    let mut altered = JSON.to_vec();
    altered[0] = b'[';
    assert!(matches!(
        Prepared::new(&altered, GAPPED),
        Err(Error::Binding)
    ));
    let crlf = String::from_utf8(GAPPED.to_vec())
        .unwrap()
        .replace('\n', "\r\n");
    assert!(matches!(
        Prepared::new(JSON, crlf.as_bytes()),
        Err(Error::Binding)
    ));
    assert!(matches!(
        Prepared::new(&vec![b' '; input::MAX_JSON + 1], GAPPED),
        Err(Error::Limit(_))
    ));
    assert!(matches!(
        Prepared::new(JSON, &vec![b'A'; input::MAX_FASTA + 1]),
        Err(Error::Limit(_))
    ));
}

#[test]
fn joins_use_exact_unique_coding_dna_never_opaque_names_or_first_match() {
    let (mut raw, records) = sources();
    raw.allele_descriptions[1].coding_sequence = raw.allele_descriptions[0].coding_sequence.clone();
    assert!(joins(&raw, &records).is_err());
    let (raw, mut records) = sources();
    records[0] = records[1].clone();
    assert!(joins(&raw, &records).is_err());
    let (mut raw, records) = sources();
    let label = raw.allele_descriptions[0].label.clone();
    raw.allele_descriptions[1].label = label.clone();
    let d = derive(&raw, &records).unwrap();
    assert_eq!(
        d.joins.iter().filter(|j| j.provider_label == label).count(),
        2
    );
    let refs = d
        .joins
        .iter()
        .filter(|j| j.provider_label == label)
        .map(|j| &j.allele_ref)
        .collect::<BTreeSet<_>>();
    assert_eq!(refs.len(), 2);
}

#[test]
fn unchanged_ungapped_dna_does_not_hide_partial_codon_gaps_or_coordinate_disagreement() {
    let (raw, mut records) = sources();
    let record = records
        .iter_mut()
        .find(|r| r.label.starts_with("IGHV"))
        .unwrap();
    let mut bytes = record.sequence.as_bytes().to_vec();
    let dot = bytes.iter().position(|&b| b == b'.').unwrap();
    assert!(dot > 0 && bytes[dot - 1] != b'.');
    bytes.swap(dot - 1, dot);
    record.sequence = String::from_utf8(bytes).unwrap();
    assert!(joins(&raw, &records).is_ok());
    assert!(derive(&raw, &records).is_err());
    let (mut raw, records) = sources();
    let a = raw
        .allele_descriptions
        .iter_mut()
        .find(|a| a.functional && a.sequence_type == "V")
        .unwrap();
    a.v_gene_delineations.as_mut().unwrap()[0].cdr1_start += 1;
    assert!(derive(&raw, &records).is_err());
    let (mut raw, records) = sources();
    raw.allele_descriptions
        .iter_mut()
        .find(|a| a.sequence_type == "J")
        .unwrap()
        .j_codon_frame = Some(1);
    assert!(derive(&raw, &records).is_err());
}

#[test]
fn identical_reference_proteins_retain_all_source_alleles_and_ties() {
    let prepared = Prepared::new(JSON, GAPPED).unwrap();
    let d = prepared.dataset();
    for group in &d.duplicate_groups {
        assert!(group.allele_refs.len() > 1);
        let members = group.allele_refs.iter().collect::<BTreeSet<_>>();
        assert_eq!(members.len(), group.allele_refs.len());
        for member in &group.allele_refs {
            assert!(
                d.v_entries.iter().any(|v| &v.allele_ref == member)
                    || d.j_segments.iter().any(|j| &j.allele_ref == member)
            );
        }
    }
    // A downstream scorer may encounter tied reference proteins; curation does
    // not discard an allele or turn the group into a unique origin claim.
    let first = d.v_entries[0].clone();
    let mut duplicate = first.clone();
    duplicate.allele_ref = "synthetic-distinct-allele-same-protein".into();
    let expected = BTreeSet::from([first.allele_ref.clone(), duplicate.allele_ref.clone()]);
    let groups = duplicate_groups(&[first, duplicate], &[]);
    assert_eq!(groups.len(), 1);
    assert_eq!(groups[0].kind, DuplicateKind::VAlignedProtein);
    assert_eq!(
        groups[0]
            .allele_refs
            .iter()
            .cloned()
            .collect::<BTreeSet<_>>(),
        expected
    );
}

#[test]
fn every_retained_translation_agrees_with_independently_decoded_primary_code_one() {
    let oracle: serde_json::Value = serde_json::from_str(include_str!(
        "../../../data/ogrdb-mouse-c57bl6-igh-v5/independent-oracle.json"
    ))
    .unwrap();
    assert_eq!(
        oracle["ncbi_code1_gc_prt_sha256"],
        "2aecbdd09ecc5e0d233c105090a593737298ef93d88f6bc0c1bef26df6051a04"
    );
    let prepared = Prepared::new(JSON, GAPPED).unwrap();
    let d = prepared.dataset();
    assert_eq!(oracle["v_entries"].as_array().unwrap().len(), 95);
    assert_eq!(oracle["j_entries"].as_array().unwrap().len(), 4);
    assert!(d.duplicate_groups.is_empty());
    for row in oracle["v_entries"].as_array().unwrap() {
        let v = d
            .v_entries
            .iter()
            .find(|v| v.label == row["label"].as_str().unwrap())
            .unwrap();
        assert_eq!(
            v.aligned_protein,
            row["aligned104_protein"].as_str().unwrap()
        );
        let join = d
            .joins
            .iter()
            .find(|j| j.allele_ref == v.allele_ref)
            .unwrap();
        assert_eq!(join.gapped_label, row["fasta_label"].as_str().unwrap());
    }
    for row in oracle["j_entries"].as_array().unwrap() {
        let j = d
            .j_segments
            .iter()
            .find(|j| j.label == row["label"].as_str().unwrap())
            .unwrap();
        assert_eq!(j.fr4_protein, row["fr4_protein"].as_str().unwrap());
        assert_eq!(
            u64::from(j.provider_codon_frame),
            row["provider_frame"].as_u64().unwrap()
        );
        assert_eq!(
            j.provider_j_cdr3_end as u64,
            row["provider_anchor_1based"].as_u64().unwrap()
        );
    }
}

#[test]
fn compact_runtime_is_exact_prepared_projection_not_a_dropdown_raw_parser() {
    let prepared = Prepared::new(JSON, GAPPED).unwrap();
    let expected = prepared.runtime_corpus().unwrap();
    let compact = crate::references::load(&crate::references::Species::MouseC57Bl6).unwrap();
    assert_eq!(
        serde_json::to_value(&expected).unwrap(),
        serde_json::to_value(compact).unwrap()
    );
    assert_eq!(expected.revision().unwrap(), compact.revision().unwrap());
}
