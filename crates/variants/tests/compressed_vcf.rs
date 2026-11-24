use collagene_variants::{Collection, ExportFormat, Format, Limits};
use std::io::Write;
const VCF:&[u8]=b"##fileformat=VCFv4.3\n##contig=<ID=c,length=2>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nc\t1\t.\tA\tC\t.\tPASS\t.\n";
fn gzip(bytes: &[u8]) -> Vec<u8> {
    let mut w = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
    w.write_all(bytes).unwrap();
    w.finish().unwrap()
}
#[test]
fn gzip_members_preserve_exact_source_and_reject_crc_trailing_truncation_and_expansion() {
    let mut bytes = gzip(&VCF[..40]);
    bytes.extend(gzip(&VCF[40..]));
    let c = Collection::import(&bytes, Format::Vcf, &Limits::default()).unwrap();
    assert_eq!(c.len(), 1);
    assert_eq!(
        c.export(ExportFormat::Original, bytes.len()).unwrap().bytes,
        bytes
    );
    assert!(c.retained_charge() >= bytes.len());
    let mut corrupt = bytes.clone();
    let n = corrupt.len();
    corrupt[n - 8] ^= 1;
    assert!(Collection::import(&corrupt, Format::Vcf, &Limits::default()).is_err());
    let mut trailing = bytes.clone();
    trailing.push(0);
    assert!(Collection::import(&trailing, Format::Vcf, &Limits::default()).is_err());
    assert!(
        Collection::import(&bytes[..bytes.len() - 1], Format::Vcf, &Limits::default()).is_err()
    );
    let limits = Limits {
        decoded_bytes: VCF.len() - 1,
        ..Default::default()
    };
    assert!(Collection::import(&bytes, Format::Vcf, &limits).is_err());
    assert!(Collection::import(&gzip(VCF), Format::Bcf, &Limits::default()).is_err());
}
