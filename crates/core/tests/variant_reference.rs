use collagene_core::{
    Molecule, Topology,
    variant_reference::{self as v, exchange},
};
const SOURCE:&[u8]=b"##fileformat=VCFv4.3\n##contig=<ID=chr1,length=7>\n##INFO=<ID=END,Number=1,Type=Integer,Description=\"Declared end\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t3\t.\tCC\tCTC\t.\tPASS\t.\nchr1\t3\t.\tC\t<DEL>\t.\tPASS\tEND=7\n";
fn options(c: &exchange::Collection, m: &Molecule) -> v::Options {
    v::Options {
        source_blake3: c.source_blake3().into(),
        record_index: 0,
        reference_name: "chr1".into(),
        expected_reference_revision: m.revision(),
        alternate_index: 1,
    }
}
#[test]
fn literal_allele_requires_actual_ref_and_keeps_independent_inputs() {
    let c =
        exchange::Collection::import(SOURCE, exchange::Format::Vcf, &exchange::Limits::default())
            .unwrap();
    let original = Molecule::dna("arbitrary display name", "AACCGTT", Topology::Linear).unwrap();
    let selected = options(&c, &original);
    let report = v::check_record(&c, &original, &selected).unwrap();
    assert_eq!(report.outcome, v::Outcome::ReferenceAlleleMatches);
    assert_eq!((report.start, report.end), (Some(2), Some(4)));
    let preview = v::preview(&c, &original, &selected).unwrap();
    assert_eq!(preview.molecule.sequence.as_str(), "AACTCGTT");
    assert_eq!(
        preview.molecule.metadata["source:legacy:chemistry_unmapped"],
        "true"
    );
    assert!(collagene_core::chemistry::Duplex::from_molecule(&preview.molecule).is_err());
    assert_eq!(original.sequence.as_str(), "AACCGTT");
    assert_eq!(c.original_bytes(), SOURCE);
    let mismatch = Molecule::dna("chr1", "AATCGTT", Topology::Linear).unwrap();
    assert_eq!(
        v::check_record(&c, &mismatch, &options(&c, &mismatch))
            .unwrap()
            .outcome,
        v::Outcome::ReferenceAlleleMismatch
    );
    assert!(v::preview(&c, &mismatch, &options(&c, &mismatch)).is_err());
}
#[test]
fn names_and_symbolic_claims_cannot_authorize_an_edit() {
    let c =
        exchange::Collection::import(SOURCE, exchange::Format::Vcf, &exchange::Limits::default())
            .unwrap();
    let original = Molecule::dna("chr1", "AACCGTT", Topology::Linear).unwrap();
    let mut selected = options(&c, &original);
    selected.record_index = 1;
    assert!(matches!(
        v::check_record(&c, &original, &selected).unwrap().outcome,
        v::Outcome::Unsupported { .. }
    ));
    assert!(v::preview(&c, &original, &selected).is_err());
    selected = options(&c, &original);
    selected.reference_name = "other".into();
    assert!(v::check_record(&c, &original, &selected).is_err());
    selected = options(&c, &original);
    selected.source_blake3 = "0".repeat(64);
    assert!(v::check_record(&c, &original, &selected).is_err());
    let ambiguous = Molecule::dna("chr1", "AANCGTT", Topology::Linear).unwrap();
    assert_eq!(
        v::check_record(&c, &ambiguous, &options(&c, &ambiguous))
            .unwrap()
            .outcome,
        v::Outcome::AmbiguousReference
    );
    assert!(v::preview(&c, &ambiguous, &options(&c, &ambiguous)).is_err());
}
