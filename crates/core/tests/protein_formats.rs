use collagene_core::protein_formats::{self as protein, ExportOptions, InsdcDate, Profile};
use collagene_core::{Alphabet, Location, Molecule, Sequence, Topology};
use std::collections::BTreeMap;
const GP: &str = "LOCUS       GP_FIXTURE                    12 aa            linear   PRI 24-NOV-2026\nDEFINITION  Independent protein fixture\n            with continued description.\nACCESSION   NP_000001 ALT_000001\nVERSION     NP_000001.1\nDBSOURCE    unrelated nucleotide reference retained as documentary text\nFEATURES             Location/Qualifiers\n     source          1..12\n                     /organism=\"Synthetic construct\"\n     Protein         1..12\n                     /product=\"Independent protein\"\n     Region          join(3..4,8..9)\n                     /note=\"two joined\n                     annotated pieces with \"\"quoted\"\" text\"\n     Site            <2..>4\n                     /site_type=\"active\"\n     CDS             1..12\n                     /coded_by=\"complement(NC_000001.1:10..45)\"\nORIGIN\n        1 mkwvtfisll fl\n//\n";
// Freely readable UniProt P69905 excerpt: selected header/FT lines and complete
// sequence/SQ checksum, not a claim to include its full annotation/citation set.
// https://rest.uniprot.org/uniprotkb/P69905.txt (retrieved2026-11-26).
const HBA: &str = "ID   HBA_HUMAN               Reviewed;         142 AA.\nAC   P69905; P01922;\nDE   RecName: Full=Hemoglobin subunit alpha;\nOS   Homo sapiens (Human).\nFT   INIT_MET        1\nFT                   /note=\"Removed\"\nFT   CHAIN           2..142\nFT                   /note=\"Hemoglobin subunit alpha\"\nFT                   /id=\"PRO_0000052653\"\nFT   PEPTIDE         96..104\nFT                   /note=\"Hemopressin\"\nSQ   SEQUENCE   142 AA;  15258 MW;  15E13666573BBBAE CRC64;\n     MVLSPADKTN VKAAWGKVGA HAGEYGAEAL ERMFLSFPTT KTYFPHFDLS HGSAQVKGHG\n     KKVADALTNA VAHVDDMPNA LSALSDLHAH KLRVDPVNFK LLSHCLLVTL AAHLPAEFTP\n     AVHASLDKFL ASVSTVLTSK YR\n//\n";
fn range(start: usize, end: usize) -> Location {
    Location::Range {
        start,
        end,
        before: false,
        after: false,
    }
}
#[test]
fn genpept_is_an_explicit_protein_profile_with_joined_fuzzy_features_and_coded_by_documentary_only()
{
    assert_eq!(protein::detect(GP.as_bytes()), Some(Profile::GenPept));
    let report = protein::import_bytes(GP.as_bytes(), Profile::GenPept).unwrap();
    let m = &report.molecules[0];
    assert_eq!(m.sequence.alphabet(), Alphabet::Protein);
    assert_eq!(m.sequence.as_str(), "MKWVTFISLLFL");
    assert_eq!(m.topology, Topology::Linear);
    assert_eq!(m.metadata["protein:accessions"], "NP_000001;ALT_000001");
    assert_eq!(
        m.metadata["protein:description"],
        "Independent protein fixture with continued description."
    );
    assert_eq!(
        m.features[2].location,
        Location::Join(vec![range(2, 4), range(7, 9)])
    );
    assert!(m.features[2].qualifiers.contains(&(
        "note".into(),
        Some("two joined annotated pieces with \"quoted\" text".into())
    )));
    assert_eq!(
        m.features[3].location,
        Location::Range {
            start: 1,
            end: 4,
            before: true,
            after: true
        }
    );
    assert!(m.features[3].needs_review);
    assert_eq!(m.features[4].location, range(0, 12));
    assert!(m.features[4].qualifiers.contains(&(
        "coded_by".into(),
        Some("complement(NC_000001.1:10..45)".into())
    )));
    assert!(m.original_genbank.is_none());
    assert!(
        m.metadata
            .keys()
            .all(|k| !k.starts_with("source:physical") && !k.contains("verified_material"))
    );
    assert_eq!(protein::original_bytes(m).unwrap(), GP.as_bytes());
    let out = protein::export(&report.molecules, Profile::GenPept).unwrap();
    assert!(out.text.contains("12 aa") && !out.text.contains("12 bp"));
    let again = protein::import_bytes(out.text.as_bytes(), Profile::GenPept).unwrap();
    assert_eq!(again.molecules[0].sequence, m.sequence);
    assert_eq!(
        again.molecules[0].features[2].location,
        m.features[2].location
    );
    assert_eq!(
        again.molecules[0].features[3].location,
        m.features[3].location
    );
}
#[test]
fn actual_uniprot_sequence_crc64_and_literal_residue_coordinates_are_independent_oracles() {
    let report = protein::import_bytes(HBA.as_bytes(), Profile::UniProt).unwrap();
    let m = &report.molecules[0];
    assert_eq!(m.name, "HBA_HUMAN");
    assert_eq!(m.sequence.len(), 142);
    assert_eq!(&m.sequence.as_str()[95..104], "PVNFKLLSH");
    assert_eq!(m.features[0].location, range(0, 1));
    assert_eq!(m.features[1].location, range(1, 142));
    assert_eq!(m.features[2].location, range(95, 104));
    assert_eq!(
        protein::crc64(m.sequence.as_str().as_bytes()),
        "15E13666573BBBAE"
    );
    assert_eq!(protein::original_bytes(m).unwrap(), HBA.as_bytes());
    let exported = protein::export(&report.molecules, Profile::UniProt).unwrap();
    assert!(exported.text.contains("Unreviewed;"));
    assert!(!exported.text.contains("Reviewed;"));
    assert!(exported.text.contains("15E13666573BBBAE CRC64"));
    let again = protein::import_bytes(exported.text.as_bytes(), Profile::UniProt).unwrap();
    assert_eq!(again.molecules[0].sequence, m.sequence);
    assert_eq!(again.molecules[0].features[1].location, range(1, 142));
    assert!(
        protein::import_bytes(
            HBA.replace("15E13666573BBBAE", "0000000000000000")
                .as_bytes(),
            Profile::UniProt
        )
        .is_err()
    );
    assert!(
        protein::import_bytes(
            HBA.replace("142 AA.", "143 AA.").as_bytes(),
            Profile::UniProt
        )
        .is_err()
    );
}
#[test]
fn legacy_ft_columns_modern_quoted_continuations_and_unknown_isoform_locations_remain_distinct() {
    let text = "ID   SYNTH                Unreviewed;       12 AA.\nAC   P12345;\nFT   CHAIN         1     12       Mature protein.\nFT                                Continued legacy prose.\nFT                                /FTId=PRO_0000000001\nFT   REGION          ?3..>10\nFT                   /note=\"modern wrapped\nFT                   qualifier text\"\nFT   DOMAIN          P12345-2:2..5\nSQ   SEQUENCE   12 AA;\n     MKWVTFISLL FL\n//\n";
    let report = protein::import_bytes(text.as_bytes(), Profile::UniProt).unwrap();
    let m = &report.molecules[0];
    assert_eq!(m.features[0].location, range(0, 12));
    assert!(m.features[0].qualifiers.contains(&(
        "note".into(),
        Some("Mature protein. Continued legacy prose.".into())
    )));
    assert_eq!(
        m.features[1].location,
        Location::Unresolved("?3..>10".into())
    );
    assert!(
        m.features[1]
            .qualifiers
            .contains(&("note".into(), Some("modern wrapped qualifier text".into())))
    );
    assert_eq!(
        m.features[2].location,
        Location::Unresolved("P12345-2:2..5".into())
    );
    assert!(m.features.iter().all(|f| f.needs_review));
    let out = protein::export(&report.molecules, Profile::UniProt).unwrap();
    assert!(out.text.contains("/id=\"PRO_0000000001\""));
    assert!(out.text.contains("?3..>10") && out.text.contains("P12345-2:2..5"));
}
#[test]
fn crosslinks_are_endpoint_pairs_and_never_contiguous_protein_feature_intervals() {
    let text = "ID   BOND_LITERAL         Unreviewed;        12 AA.\nAC   P12345;\nFT   DISULFID        2..11\nSQ   SEQUENCE   12 AA;\n     ACDEFGHIKL CN\n//\n";
    let mut report = protein::import_bytes(text.as_bytes(), Profile::UniProt).unwrap();
    let m = &report.molecules[0];
    assert_eq!(
        m.features[0].location,
        Location::Join(vec![range(1, 2), range(10, 11)])
    );
    assert_eq!(m.features[0].location.extract(&m.sequence).unwrap(), "CC");
    let actual_material =
        collagene_core::protein_properties::analyze_protein(m, &Default::default()).unwrap();
    assert!(actual_material.mass.is_none());
    let out = protein::export(&report.molecules, Profile::UniProt).unwrap();
    let round_trip = protein::import_bytes(out.text.as_bytes(), Profile::UniProt).unwrap();
    assert_eq!(
        round_trip.molecules[0].features[0].location,
        Location::Join(vec![range(1, 2), range(10, 11)])
    );
    assert_eq!(
        round_trip.molecules[0].features[0]
            .location
            .extract(&round_trip.molecules[0].sequence)
            .unwrap(),
        "CC"
    );
    report.molecules[0].features[0].location = range(1, 11);
    assert!(protein::export(&report.molecules, Profile::UniProt).is_err());
}
#[test]
fn multirecord_crlf_original_chunks_are_exact_and_edited_export_does_not_replay_stale_source() {
    let padding = "COMMENT     original-only field\r\n".repeat(2500);
    let first = GP
        .replacen("FEATURES", &format!("{padding}FEATURES"), 1)
        .replace('\n', "\r\n");
    let input = format!("{first}\r\n{}\r\n", GP.replace("GP_FIXTURE", "SECOND"));
    let report = protein::import_bytes(input.as_bytes(), Profile::GenPept).unwrap();
    assert_eq!(report.molecules.len(), 2);
    let original = report
        .molecules
        .iter()
        .map(|m| protein::original_bytes(m).unwrap())
        .collect::<Vec<_>>()
        .concat();
    assert_eq!(original, input.as_bytes());
    assert!(
        report.molecules[0].metadata["protein:original:chunks"]
            .parse::<usize>()
            .unwrap()
            > 1
    );
    let mut edited = report.molecules.clone();
    edited[0].sequence = Sequence::new(Alphabet::Protein, "ACDEFGHIKLMN").unwrap();
    let out = protein::export(&edited, Profile::GenPept).unwrap();
    assert!(out.text.contains("acdefghikl mn"));
    assert!(!out.text.contains("original-only field"));
    assert_eq!(
        protein::original_bytes(&edited[0]).unwrap(),
        first.as_bytes()
    );
}
#[test]
fn explicit_synthetic_export_identifiers_and_date_are_validated_without_inventing_database_claims()
{
    let m = Molecule {
        name: "SyntheticProtein".into(),
        sequence: Sequence::new(Alphabet::Protein, "ACDEFGHIKLMN").unwrap(),
        topology: Topology::Linear,
        features: vec![],
        metadata: BTreeMap::new(),
        original_genbank: None,
    };
    assert!(protein::export(&[m.clone()], Profile::GenPept).is_err());
    assert!(protein::export(&[m.clone()], Profile::UniProt).is_err());
    let options = ExportOptions {
        locus_date: Some(InsdcDate::new("26-NOV-2026").unwrap()),
        uniprot_accession: Some("P12345".into()),
        entry_name: Some("LOCAL_PROTEIN".into()),
    };
    let gp = protein::export_with_options(&[m.clone()], Profile::GenPept, &options).unwrap();
    assert!(gp.text.contains("26-NOV-2026"));
    let sp = protein::export_with_options(&[m.clone()], Profile::UniProt, &options).unwrap();
    assert!(sp.text.contains("AC   P12345;") && sp.text.contains("Unreviewed;"));
    for date in [
        "29-FEB-2025",
        "31-APR-2026",
        "00-JAN-2026",
        "01-jan-2026",
        "1-JAN-2026",
    ] {
        assert!(InsdcDate::new(date).is_err());
    }
    assert!(InsdcDate::new("29-FEB-2024").is_ok());
    let mut ambiguous = m.clone();
    ambiguous.sequence = Sequence::new(Alphabet::Protein, "ABZX").unwrap();
    assert!(protein::export_with_options(&[ambiguous], Profile::UniProt, &options).is_err());
    assert!(protein::export_with_options(&[m.clone(), m], Profile::UniProt, &options).is_err());
}
#[test]
fn malformed_profiles_counters_records_nesting_and_bounds_fail_without_guessing_locations() {
    assert!(
        protein::import_bytes(
            GP.replace("12 aa", "12 bp DNA").as_bytes(),
            Profile::GenPept
        )
        .is_err()
    );
    assert!(
        protein::import_bytes(GP.replace("1 mkwv", "2 mkwv").as_bytes(), Profile::GenPept).is_err()
    );
    assert!(
        protein::import_bytes(
            GP.replace("mkwvtfisll", "mkwv2fisll").as_bytes(),
            Profile::GenPept
        )
        .is_err()
    );
    assert!(
        protein::import_bytes(GP.trim_end_matches("//\n").as_bytes(), Profile::GenPept).is_err()
    );
    let unknown = protein::import_bytes(
        GP.replace("join(3..4,8..9)", "50..60").as_bytes(),
        Profile::GenPept,
    )
    .unwrap();
    assert_eq!(
        unknown.molecules[0].features[2].location,
        Location::Unresolved("50..60".into())
    );
    let deep = format!("{}1{}", "join(".repeat(20), ")".repeat(20));
    assert!(
        protein::import_bytes(
            GP.replace("join(3..4,8..9)", &deep).as_bytes(),
            Profile::GenPept
        )
        .is_err()
    );
    assert!(protein::import_bytes(&vec![b'x'; protein::MAX_BYTES + 1], Profile::UniProt).is_err());
    let mut report = protein::import_bytes(GP.as_bytes(), Profile::GenPept).unwrap();
    report.molecules[0]
        .metadata
        .insert("protein:original:blake3".into(), "0".repeat(64));
    assert!(protein::original_bytes(&report.molecules[0]).is_err());
}

#[test]
fn generic_dispatch_preserves_original_protein_bytes_and_oversize_aa_never_falls_back_to_dna() {
    let imported = collagene_core::formats::import_bytes(GP.as_bytes()).unwrap();
    assert_eq!(imported.molecules[0].sequence.alphabet(), Alphabet::Protein);
    assert_eq!(
        protein::original_bytes(&imported.molecules[0]).unwrap(),
        GP.as_bytes()
    );
    assert!(
        collagene_core::formats::export_with_report(
            &imported.molecules[0],
            collagene_core::formats::TextFormat::GenPept
        )
        .unwrap()
        .text
        .contains("12 aa")
    );
    let imported = collagene_core::formats::import_text(HBA).unwrap();
    assert_eq!(imported.molecules[0].sequence.alphabet(), Alphabet::Protein);
    let mut oversized = GP.as_bytes().to_vec();
    oversized.resize(protein::MAX_BYTES + 1, b' ');
    assert_eq!(protein::detect(&oversized), Some(Profile::GenPept));
    assert!(matches!(
        collagene_core::formats::import_bytes(&oversized),
        Err(collagene_core::Error::TooLarge)
    ));
}
