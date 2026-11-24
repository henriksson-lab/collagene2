#![doc = include_str!("../README.md")]
//! Bounded SAM/BAM exchange. Alignments are imported assertions, never validated variants or coverage.
mod bam;
mod bgzf;
mod mapping;
mod sam;
pub use bam::{BamReader, BamWriter};
pub use bgzf::{BgzfReader, BgzfWriter, VirtualOffset};
pub use mapping::{CoordinateMap, MapRun, QueryProjection};
pub use sam::{SamReader, SamWriter};
use std::{fmt, io};
pub type Result<T> = std::result::Result<T, Error>;
#[derive(Debug)]
pub enum Error {
    Io(io::Error),
    Invalid(&'static str),
    Unsupported(&'static str),
    Limit(&'static str),
    Poisoned,
}
impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Io(e) => write!(f, "{e}"),
            Self::Invalid(s) => write!(f, "invalid alignment: {s}"),
            Self::Unsupported(s) => write!(f, "unsupported alignment: {s}"),
            Self::Limit(s) => write!(f, "alignment limit: {s}"),
            Self::Poisoned => write!(f, "reader failed; no completion receipt is available"),
        }
    }
}
impl std::error::Error for Error {}
impl From<std::fmt::Error> for Error {
    fn from(_: std::fmt::Error) -> Self {
        Self::Limit("encoded SAM record bytes")
    }
}
impl From<io::Error> for Error {
    fn from(e: io::Error) -> Self {
        Self::Io(e)
    }
}
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Limits {
    pub max_header_bytes: usize,
    pub max_record_bytes: usize,
    pub max_references: usize,
    pub max_sequence_bases: usize,
    pub max_cigar_ops: usize,
    pub max_aux_tags: usize,
    pub max_array_values: usize,
    pub max_records: u64,
    pub max_input_bytes: u64,
    pub max_uncompressed_bytes: u64,
    pub max_bgzf_blocks: u64,
}
impl Default for Limits {
    fn default() -> Self {
        Self {
            max_header_bytes: 8 << 20,
            max_record_bytes: 4 << 20,
            max_references: 100_000,
            max_sequence_bases: 1_000_000,
            max_cigar_ops: 100_000,
            max_aux_tags: 1024,
            max_array_values: 1_000_000,
            max_records: 10_000_000,
            max_input_bytes: 64 << 30,
            max_uncompressed_bytes: 256 << 30,
            max_bgzf_blocks: 10_000_000,
        }
    }
}
impl Limits {
    pub fn validate(&self) -> Result<()> {
        if self.max_header_bytes == 0
            || self.max_header_bytes > 16 << 20
            || self.max_record_bytes == 0
            || self.max_record_bytes > 32 << 20
            || self.max_references == 0
            || self.max_references > 1_000_000
            || self.max_sequence_bases == 0
            || self.max_sequence_bases > 10_000_000
            || self.max_cigar_ops == 0
            || self.max_cigar_ops > 1_000_000
            || self.max_aux_tags == 0
            || self.max_aux_tags > 4096
            || self.max_array_values == 0
            || self.max_array_values > 4_000_000
            || self.max_records == 0
            || self.max_input_bytes == 0
            || self.max_uncompressed_bytes == 0
            || self.max_bgzf_blocks == 0
        {
            return Err(Error::Limit("invalid resource policy"));
        }
        Ok(())
    }
}
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Format {
    Sam,
    Bam,
}
/// Quick triage only: gzip is a candidate, not proof of BAM/BGZF validity.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Probe {
    SamCandidate,
    BgzfCandidate,
    Unknown,
}
pub fn detect(prefix: &[u8]) -> Result<Probe> {
    if prefix.starts_with(b"CRAM") && prefix.get(4).is_some_and(|v| *v < 32 && *v != b'\t') {
        return Err(Error::Unsupported(
            "CRAM reference resolution is not implemented",
        ));
    }
    if prefix.starts_with(b"BAI\x01") || prefix.starts_with(b"CSI\x01") {
        return Err(Error::Unsupported(
            "alignment indexes are not alignment records",
        ));
    }
    if prefix.starts_with(b"BCF") && prefix.get(3).is_some_and(|v| *v < 32 && *v != b'\t')
        || prefix.starts_with(b"##fileformat=VCF") && !prefix.contains(&b'\t')
    {
        return Err(Error::Unsupported("variant files are not read alignments"));
    }
    if prefix.starts_with(&[31, 139]) {
        Ok(Probe::BgzfCandidate)
    } else if prefix.starts_with(b"BAM\x01") {
        Err(Error::Unsupported(
            "uncompressed BAM transport; BGZF is required",
        ))
    } else if prefix.starts_with(b"@") || prefix.contains(&b'\t') {
        Ok(Probe::SamCandidate)
    } else {
        Ok(Probe::Unknown)
    }
}
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Reference {
    pub name: String,
    pub length: u32,
    pub circular: bool,
}
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Header {
    pub(crate) raw: Vec<u8>,
    pub(crate) references: std::sync::Arc<[Reference]>,
    pub(crate) bam_original: Option<Vec<u8>>,
    pub(crate) dictionary_identity: [u8; 32],
    pub(crate) reference_index: std::sync::Arc<std::collections::BTreeMap<String, usize>>,
    pub(crate) sam_dictionary_matches: bool,
}
impl Header {
    pub(crate) fn new(raw: Vec<u8>, references: Vec<Reference>) -> Self {
        let mut header = Self {
            raw,
            references: [].into(),
            bam_original: None,
            dictionary_identity: [0; 32],
            reference_index: Default::default(),
            sam_dictionary_matches: true,
        };
        header.set_references(references);
        header
    }
    pub(crate) fn set_references(&mut self, references: Vec<Reference>) {
        let mut hash = blake3::Hasher::new();
        hash.update(b"alignment-reference-dictionary-v1");
        hash.update(&(self.raw.len() as u64).to_le_bytes());
        hash.update(&self.raw);
        let mut reference_index = std::collections::BTreeMap::new();
        for (i, r) in references.iter().enumerate() {
            hash.update(&(r.name.len() as u64).to_le_bytes());
            hash.update(r.name.as_bytes());
            hash.update(&r.length.to_le_bytes());
            hash.update(&[u8::from(r.circular)]);
            reference_index.insert(r.name.clone(), i);
        }
        self.reference_index = std::sync::Arc::new(reference_index);
        self.dictionary_identity = *hash.finalize().as_bytes();
        self.references = references.into();
    }
    pub fn reference(&self, name: &str) -> Option<&Reference> {
        self.reference_index
            .get(name)
            .and_then(|i| self.references.get(*i))
    }
    pub(crate) fn validate_policy(&self, limits: &Limits) -> Result<()> {
        limits.validate()?;
        if self.raw.len() > limits.max_header_bytes
            || self.references.len() > limits.max_references
            || self
                .raw
                .split_inclusive(|b| *b == b'\n')
                .any(|line| line.len() > limits.max_record_bytes)
        {
            return Err(Error::Limit("header bytes/reference count/line bytes"));
        }
        Ok(())
    }
    pub fn original_bam_header(&self) -> Option<&[u8]> {
        self.bam_original.as_deref()
    }
    pub fn sam_bytes(&self) -> &[u8] {
        &self.raw
    }
    pub fn references(&self) -> &[Reference] {
        &self.references
    }
    pub fn from_sam(bytes: &[u8], limits: &Limits) -> Result<Self> {
        sam::parse_header(bytes, limits)
    }
}
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CigarKind {
    Match,
    Insertion,
    Deletion,
    Skip,
    SoftClip,
    HardClip,
    Pad,
    Equal,
    Difference,
}
impl CigarKind {
    pub fn consumes_query(self) -> bool {
        matches!(
            self,
            Self::Match | Self::Insertion | Self::SoftClip | Self::Equal | Self::Difference
        )
    }
    pub fn consumes_reference(self) -> bool {
        matches!(
            self,
            Self::Match | Self::Deletion | Self::Skip | Self::Equal | Self::Difference
        )
    }
    pub fn symbol(self) -> char {
        b"MIDNSHP=X"[self.code() as usize] as char
    }
    pub(crate) fn code(self) -> u32 {
        match self {
            Self::Match => 0,
            Self::Insertion => 1,
            Self::Deletion => 2,
            Self::Skip => 3,
            Self::SoftClip => 4,
            Self::HardClip => 5,
            Self::Pad => 6,
            Self::Equal => 7,
            Self::Difference => 8,
        }
    }
    pub(crate) fn from_code(c: u32) -> Result<Self> {
        Ok(match c {
            0 => Self::Match,
            1 => Self::Insertion,
            2 => Self::Deletion,
            3 => Self::Skip,
            4 => Self::SoftClip,
            5 => Self::HardClip,
            6 => Self::Pad,
            7 => Self::Equal,
            8 => Self::Difference,
            _ => return Err(Error::Unsupported("CIGAR operation code")),
        })
    }
}
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CigarOp {
    pub len: u32,
    pub kind: CigarKind,
}
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum AuxValue {
    Char(u8),
    Integer(i64),
    FloatBits(u32),
    String(Vec<u8>),
    Hex(Vec<u8>),
    Array { element_type: u8, values: Vec<u32> },
}
/// `wire_type` retains BAM signedness/width; array values retain original integer/float bits.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Aux {
    pub tag: [u8; 2],
    pub wire_type: u8,
    pub value: AuxValue,
}
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Fields {
    pub name: Option<String>,
    pub flags: u16,
    pub reference: Option<String>,
    pub start: Option<u32>,
    pub mapping_quality: Option<u8>,
    pub cigar: Vec<CigarOp>,
    pub mate_reference: Option<String>,
    pub mate_start: Option<u32>,
    pub template_length: i32,
    pub sequence: Option<String>,
    pub qualities: Option<Vec<u8>>,
    pub aux: Vec<Aux>,
}
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Record {
    pub(crate) fields: Fields,
    pub(crate) format: Format,
    pub(crate) original: Vec<u8>,
    pub(crate) dictionary: [u8; 32],
    pub(crate) virtual_offset: Option<VirtualOffset>,
    pub(crate) long_cigar: bool,
}
/// A reported symbol is not automatically a called DNA base.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SequenceSymbol {
    Iupac(u8),
    ReferenceEqual,
    Unresolved(u8),
}
/// Borrowed, loss-reviewed handoff to a host's shared DNA Sequence parser.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct DnaCallView<'a> {
    sequence: &'a str,
    pub case_normalization_required: bool,
}
impl DnaCallView<'_> {
    pub fn as_str(&self) -> &str {
        self.sequence
    }
}
impl Record {
    pub fn reported_symbol(&self, index: usize) -> Option<SequenceSymbol> {
        let b = *self.fields.sequence.as_ref()?.as_bytes().get(index)?;
        Some(if b == b'=' {
            SequenceSymbol::ReferenceEqual
        } else if b"ACMGRSVTWYHKDBN".contains(&b.to_ascii_uppercase()) {
            SequenceSymbol::Iupac(b)
        } else {
            SequenceSymbol::Unresolved(b)
        })
    }
    /// Never fills reference-equal, absent or non-DNA calls. Host Sequence conversion must remain explicit.
    pub fn dna_calls(&self) -> Result<DnaCallView<'_>> {
        let sequence = self
            .fields
            .sequence
            .as_deref()
            .ok_or(Error::Unsupported("SEQ is absent; no called DNA sequence"))?;
        if !sequence
            .bytes()
            .all(|b| b"ACMGRSVTWYHKDBN".contains(&b.to_ascii_uppercase()))
        {
            return Err(Error::Unsupported(
                "reference-equal or non-IUPAC symbols require explicit reference evidence or loss-reviewed conversion",
            ));
        }
        Ok(DnaCallView {
            sequence,
            case_normalization_required: sequence.bytes().any(|b| b.is_ascii_lowercase()),
        })
    }
    pub fn fields(&self) -> &Fields {
        &self.fields
    }
    pub fn format(&self) -> Format {
        self.format
    }
    /// Exact SAM line (terminator included) or BAM block_size+record payload, not BGZF framing.
    pub fn original_bytes(&self) -> &[u8] {
        &self.original
    }
    pub fn original_blake3(&self) -> String {
        blake3::hash(&self.original).to_hex().to_string()
    }
    pub fn virtual_offset(&self) -> Option<VirtualOffset> {
        self.virtual_offset
    }
    pub fn coordinate_map(&self, header: &Header) -> Result<CoordinateMap> {
        mapping::project(self, header)
    }
}
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StreamState {
    Reading,
    Complete,
    Failed,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct StreamProgress {
    pub state: StreamState,
    pub records_emitted: u64,
    pub source_bytes_accounted: u64,
}
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Completion {
    pub format: Format,
    pub records: u64,
    pub source_bytes: u64,
    pub source_blake3: String,
    pub warnings: Vec<String>,
}
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ExportNotice {
    TextFormattingNotPreserved,
    BgzfFramingNotPreserved,
    SequenceCaseNormalized,
    AuxIntegerWidthNotPreserved,
    AuxIntegerWidthSelected,
    LongCigarTransportConverted,
}
pub(crate) fn checked_add(a: u64, b: u64, cap: u64, label: &'static str) -> Result<u64> {
    a.checked_add(b)
        .filter(|v| *v <= cap)
        .ok_or(Error::Limit(label))
}
pub(crate) fn validate_fields(f: &Fields, limits: &Limits) -> Result<()> {
    if f.template_length == i32::MIN {
        return Err(Error::Invalid("template length cannot be i32::MIN"));
    }
    if f.sequence
        .as_ref()
        .is_some_and(|s| s.len() > limits.max_sequence_bases)
        || f.cigar.len() > limits.max_cigar_ops
        || f.aux.len() > limits.max_aux_tags
    {
        return Err(Error::Limit("record elements"));
    }
    if let Some(seq) = &f.sequence {
        if seq.is_empty()
            || !seq
                .bytes()
                .all(|b| b.is_ascii_alphabetic() || b == b'=' || b == b'.')
        {
            return Err(Error::Invalid("sequence alphabet"));
        }
        if let Some(q) = &f.qualities
            && (q.len() != seq.len() || q.iter().any(|v| *v > 93))
        {
            return Err(Error::Invalid("quality length or score"));
        }
    } else if f.qualities.is_some() {
        return Err(Error::Invalid("quality without sequence"));
    }
    let mut query = 0u64;
    let mut tags = std::collections::BTreeSet::new();
    for (i, op) in f.cigar.iter().enumerate() {
        if op.len == 0 {
            return Err(Error::Unsupported("zero-length CIGAR operations"));
        }
        if op.kind.consumes_query() {
            query = checked_add(
                query,
                u64::from(op.len),
                limits.max_sequence_bases as u64,
                "query CIGAR span",
            )?;
        }
        if op.kind == CigarKind::HardClip && i != 0 && i + 1 != f.cigar.len() {
            return Err(Error::Invalid("internal hard clip"));
        }
        if op.kind == CigarKind::SoftClip
            && !(i == 0
                || i + 1 == f.cigar.len()
                || i == 1 && f.cigar[0].kind == CigarKind::HardClip
                || i + 2 == f.cigar.len() && f.cigar[i + 1].kind == CigarKind::HardClip)
        {
            return Err(Error::Invalid("internal soft clip"));
        }
    }
    if !f.cigar.is_empty() && f.sequence.as_ref().is_some_and(|s| s.len() as u64 != query) {
        return Err(Error::Invalid("CIGAR/query length mismatch"));
    }
    for a in &f.aux {
        if !a.tag[0].is_ascii_alphabetic()
            || !a.tag[1].is_ascii_alphanumeric()
            || !tags.insert(a.tag)
        {
            return Err(Error::Invalid("auxiliary tag name or duplicate"));
        }
        if matches!(&a.value,AuxValue::Array{values,..} if values.len()>limits.max_array_values) {
            return Err(Error::Limit("auxiliary array"));
        }
    }
    Ok(())
}
