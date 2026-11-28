//! Literal local-file bytes, independent of FASTA export/reconstruction.
use base64::{Engine, engine::general_purpose::STANDARD};
use collagene_core::{
    formats::sequence::{self as formats, TextFormat},
    protocol,
};
use collagene_model::{Alphabet, Topology};

fn original(molecule: &collagene_model::Molecule) -> Vec<u8> {
    let bytes = STANDARD
        .decode(&molecule.metadata["source:original_fasta_base64"])
        .unwrap();
    assert_eq!(
        molecule.metadata["source:original_fasta_blake3"],
        blake3::hash(&bytes).to_hex().to_string()
    );
    bytes
}

#[test]
fn literal_record_slices_keep_bom_crlf_indentation_and_all_surrounding_whitespace() {
    const FIRST: &[u8] = b"\xef\xbb\xbf \r\n\t >alpha literal header \r\n ac gt \r\n\t\r\n";
    const SECOND: &[u8] = b"  >beta [molecule=RNA]  \r\n au gc \r\n  \t";
    let raw = [FIRST, SECOND].concat();
    let report = formats::import_bytes(&raw).unwrap();
    assert_eq!(report.molecules.len(), 2);
    assert_eq!(report.molecules[0].name, "alpha literal header");
    assert_eq!(report.molecules[0].sequence.as_str(), "ACGT");
    assert_eq!(report.molecules[0].sequence.alphabet(), Alphabet::Dna);
    assert_eq!(report.molecules[1].name, "beta [molecule=RNA]");
    assert_eq!(report.molecules[1].sequence.as_str(), "AUGC");
    assert_eq!(report.molecules[1].sequence.alphabet(), Alphabet::Rna);
    for (molecule, expected) in report.molecules.iter().zip([FIRST, SECOND]) {
        assert_eq!(molecule.topology, Topology::Linear);
        assert!(molecule.features.is_empty());
        assert_eq!(original(molecule), expected);
        assert!(molecule.original_genbank.is_none());
        assert!(!molecule.metadata.contains_key("source:provider"));
        let reopened =
            protocol::load_document(&protocol::save_document(molecule).unwrap()).unwrap();
        assert_eq!(reopened, *molecule);
        assert_eq!(original(&reopened), expected);
    }
    assert_eq!(
        report
            .molecules
            .iter()
            .map(|molecule| original(molecule).len())
            .sum::<usize>(),
        raw.len()
    );
    let direct = formats::import_fasta(std::str::from_utf8(&raw).unwrap()).unwrap();
    assert_eq!(direct.molecules, report.molecules);
    let text = formats::import_text(std::str::from_utf8(&raw).unwrap()).unwrap();
    assert_eq!(text.molecules, report.molecules);
}

#[test]
fn raw_evidence_does_not_reinterpret_alphabet_markers_or_literal_header_names() {
    const RAW: &str = "\u{feff} \r\n>protein label [molecule=protein] \r\n ac gt \r\n \t";
    let explicit = formats::import_fasta_with_alphabet(RAW, Alphabet::Protein).unwrap();
    assert_eq!(explicit.molecules[0].sequence.alphabet(), Alphabet::Protein);
    assert_eq!(
        explicit.molecules[0].name,
        "protein label [molecule=protein]"
    );
    assert_eq!(original(&explicit.molecules[0]), RAW.as_bytes());
    assert!(formats::import_fasta_with_alphabet(RAW, Alphabet::Dna).is_err());
    assert!(formats::import_fasta(">DNA [molecule=DNA] [molecule=protein]\r\nACGT\r\n").is_err());
    assert!(formats::import_fasta(">valid\r\nACGT\r\n>empty\r\n \t\r\n").is_err());
}

#[test]
fn many_records_retain_a_disjoint_partition_not_a_copy_of_the_entire_file_per_record() {
    let mut raw = Vec::new();
    let mut expected = Vec::new();
    for index in 0..257 {
        let record = format!("  >independent-{index}\r\n{}\r\n \t\r\n", "ACGT".repeat(31));
        expected.push(record.as_bytes().to_vec());
        raw.extend_from_slice(record.as_bytes());
    }
    let report = formats::import_bytes(&raw).unwrap();
    assert_eq!(report.molecules.len(), 257);
    let mut encoded_bytes = 0usize;
    let mut decoded_bytes = 0usize;
    for (molecule, expected) in report.molecules.iter().zip(expected) {
        let actual = original(molecule);
        assert_eq!(actual, expected);
        decoded_bytes += actual.len();
        encoded_bytes += molecule.metadata["source:original_fasta_base64"].len();
    }
    assert_eq!(decoded_bytes, raw.len());
    assert!(encoded_bytes <= raw.len().div_ceil(3) * 4 + report.molecules.len() * 4);
}

#[test]
fn ordinary_large_retained_record_survives_native_metadata_bounds_without_being_plain_text() {
    let raw = format!(">large literal record\r\n{}\r\n", "A".repeat(1_000_001));
    let report = formats::import_bytes(raw.as_bytes()).unwrap();
    let molecule = &report.molecules[0];
    assert_eq!(molecule.sequence.len(), 1_000_001);
    assert_eq!(original(molecule), raw.as_bytes());
    assert!(molecule.metadata["source:original_fasta_base64"].len() > 1_048_576);
    let reopened = protocol::load_document(&protocol::save_document(molecule).unwrap()).unwrap();
    assert_eq!(original(&reopened), raw.as_bytes());
    // A new FASTA export is a normalized presentation, never the retained input.
    let exported = formats::export_with_report(molecule, TextFormat::Fasta).unwrap();
    assert_ne!(exported.text.as_bytes(), raw.as_bytes());
    assert_eq!(original(molecule), raw.as_bytes());
}

#[test]
fn atomic_input_and_aggregate_evidence_bounds_refuse_valid_looking_amplification() {
    let too_many = (0..1001)
        .map(|index| format!(">record-{index}\nA\n"))
        .collect::<String>();
    assert!(matches!(
        formats::import_bytes(too_many.as_bytes()),
        Err(collagene_model::Error::TooLarge)
    ));
    let sequence = "A".repeat(1_000_000);
    let mut raw = String::new();
    for index in 0..15 {
        raw.push_str(&format!(">valid-bounded-record-{index}\n"));
        raw.push_str(&sequence);
        raw.push('\n');
    }
    // Fifteen individually ordinary records fit the32MiB raw file cap, while
    // normalized calls plus their base64 originals would exceed the32MiB report.
    assert!(raw.len() < formats::MAX_IMPORT_BYTES);
    assert!(matches!(
        formats::import_bytes(raw.as_bytes()),
        Err(collagene_model::Error::TooLarge)
    ));
    let report = formats::ImportReport {
        molecules: vec![
            collagene_model::Molecule::dna("native row", "A", Topology::Linear).unwrap();
            1001
        ],
        warnings: vec![],
    };
    assert!(matches!(
        report.validate_bounds(),
        Err(collagene_model::Error::TooLarge)
    ));
}
