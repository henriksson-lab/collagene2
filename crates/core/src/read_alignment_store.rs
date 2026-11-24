//! Bounded, immutable SAM/BAM collections. File completion is not biological truth.
//! Kept separate from plasmid documents and alignment/assembly inference.
use crate::Error;
use collagene_read_alignments::{
    self as exchange, AuxValue, BamReader, Completion, Header, Record, SamReader,
};
pub use collagene_read_alignments::{Format, Limits as ReaderLimits};
use std::{io::Cursor, mem::size_of, sync::Arc};

#[derive(Clone, Debug)]
pub struct ImportLimits {
    pub max_source_bytes: usize,
    /// Conservative managed-buffer charge, not a measurement of process RSS.
    pub max_retained_bytes: usize,
    pub max_records: usize,
    pub reader: exchange::Limits,
}
impl Default for ImportLimits {
    fn default() -> Self {
        Self {
            max_source_bytes: 64 << 20,
            max_retained_bytes: 128 << 20,
            max_records: 100_000,
            reader: exchange::Limits {
                max_header_bytes: 1 << 20,
                max_record_bytes: 1 << 20,
                max_references: 4096,
                max_sequence_bases: 250_000,
                max_cigar_ops: 16_384,
                max_aux_tags: 128,
                max_array_values: 250_000,
                max_records: 100_000,
                max_input_bytes: 64 << 20,
                max_uncompressed_bytes: 128 << 20,
                max_bgzf_blocks: 65_536,
            },
        }
    }
}
impl ImportLimits {
    pub fn validate(&self) -> Result<(), Error> {
        if self.max_source_bytes == 0
            || self.max_source_bytes > 256 << 20
            || self.max_retained_bytes == 0
            || self.max_retained_bytes > 512 << 20
            || self.max_records == 0
            || self.max_records > 1_000_000
        {
            return Err(Error::TooLarge);
        }
        self.reader.validate().map_err(error)?;
        if self.reader.max_input_bytes > self.max_source_bytes as u64
            || self.reader.max_records > self.max_records as u64
            || self.reader.max_header_bytes > self.max_retained_bytes
            || self.reader.max_record_bytes > self.max_retained_bytes
        {
            return Err(Error::TooLarge);
        }
        Ok(())
    }
}
#[derive(Clone, Debug, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ReferenceSummary {
    pub name: String,
    pub length: u32,
    pub circular: bool,
}
#[derive(Clone, Debug, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Summary {
    pub source_blake3: String,
    pub source_bytes: u64,
    pub format: String,
    pub records: usize,
    pub flag_unmapped: usize,
    pub flag_secondary: usize,
    pub flag_supplementary: usize,
    pub references: Vec<ReferenceSummary>,
    pub warnings: Vec<String>,
    pub interpretation: String,
}
/// No Deserialize or mutable record/header access: a JSON receipt cannot recreate
/// the private completed stream or replace its source/record dictionary binding.
pub struct Collection {
    source: Arc<[u8]>,
    header: Header,
    records: Vec<Record>,
    completion: Completion,
    retained_charge: usize,
    flag_counts: [usize; 3],
}
fn error(e: impl std::fmt::Display) -> Error {
    Error::Unsupported(format!("read alignment collection: {e}"))
}
fn add(total: &mut usize, amount: usize, cap: usize) -> Result<(), Error> {
    *total = total
        .checked_add(amount)
        .filter(|v| *v <= cap)
        .ok_or(Error::TooLarge)?;
    Ok(())
}
fn string_capacity(value: &Option<String>) -> usize {
    value.as_ref().map_or(0, String::capacity)
}
fn header_charge(header: &Header) -> Result<usize, Error> {
    let mut total = size_of::<Header>();
    add(
        &mut total,
        header.sam_bytes().len().saturating_mul(2),
        usize::MAX,
    )?;
    add(
        &mut total,
        header
            .original_bam_header()
            .map_or(0, |b| b.len().saturating_mul(2)),
        usize::MAX,
    )?;
    // Includes shared reference storage and a conservative dictionary/node charge.
    for r in header.references() {
        add(
            &mut total,
            r.name
                .capacity()
                .saturating_mul(4)
                .saturating_add(size_of::<exchange::Reference>() * 2 + 256),
            usize::MAX,
        )?;
    }
    Ok(total)
}
fn record_charge(record: &Record) -> Result<usize, Error> {
    let f = record.fields();
    // Account for Vec<Record> geometric growth and the private original-byte
    // buffer; allocator bookkeeping and the caller's source slice are external.
    let mut total = size_of::<Record>() * 2 + 64;
    for n in [
        record.original_bytes().len().saturating_mul(2),
        string_capacity(&f.name),
        string_capacity(&f.reference),
        string_capacity(&f.mate_reference),
        string_capacity(&f.sequence),
        f.qualities.as_ref().map_or(0, Vec::capacity),
        f.cigar
            .capacity()
            .saturating_mul(size_of::<exchange::CigarOp>()),
        f.aux.capacity().saturating_mul(size_of::<exchange::Aux>()),
    ] {
        add(&mut total, n, usize::MAX)?;
    }
    for aux in &f.aux {
        match &aux.value {
            AuxValue::String(v) | AuxValue::Hex(v) => add(&mut total, v.capacity(), usize::MAX)?,
            AuxValue::Array { values, .. } => add(
                &mut total,
                values.capacity().saturating_mul(size_of::<u32>()),
                usize::MAX,
            )?,
            _ => {}
        }
    }
    Ok(total)
}
enum Reader<'a> {
    Sam(SamReader<Cursor<&'a [u8]>>),
    Bam(BamReader<Cursor<&'a [u8]>>),
}
impl Reader<'_> {
    fn header(&self) -> &Header {
        match self {
            Self::Sam(r) => r.header(),
            Self::Bam(r) => r.header(),
        }
    }
    fn next(&mut self) -> Result<Option<Record>, Error> {
        match self {
            Self::Sam(r) => r.next_record(),
            Self::Bam(r) => r.next_record(),
        }
        .map_err(error)
    }
    fn completion(&self) -> Option<Completion> {
        match self {
            Self::Sam(r) => r.completion(),
            Self::Bam(r) => r.completion(),
        }
    }
}
impl Collection {
    /// Publish only after strict EOF/CRC/limits and an independent exact original
    /// source hash agree. Failure discards all locally built records.
    pub fn import(source: &[u8], format: Format, limits: &ImportLimits) -> Result<Self, Error> {
        limits.validate()?;
        if source.is_empty() || source.len() > limits.max_source_bytes {
            return Err(Error::TooLarge);
        }
        let mut reader = match format {
            Format::Sam => Reader::Sam(
                SamReader::new(Cursor::new(source), limits.reader.clone()).map_err(error)?,
            ),
            Format::Bam => Reader::Bam(
                BamReader::new(Cursor::new(source), limits.reader.clone()).map_err(error)?,
            ),
        };
        let mut charge = source.len();
        add(
            &mut charge,
            header_charge(reader.header())?,
            limits.max_retained_bytes,
        )?;
        add(
            &mut charge,
            size_of::<Record>() * 4,
            limits.max_retained_bytes,
        )?;
        let mut records = Vec::new();
        let mut flags = [0usize; 3];
        while let Some(record) = reader.next()? {
            if records.len() == limits.max_records {
                return Err(Error::TooLarge);
            }
            add(
                &mut charge,
                record_charge(&record)?,
                limits.max_retained_bytes,
            )?;
            let bits = record.fields().flags;
            flags[0] += usize::from(bits & 4 != 0);
            flags[1] += usize::from(bits & 256 != 0);
            flags[2] += usize::from(bits & 2048 != 0);
            records.push(record);
        }
        let completion = reader
            .completion()
            .ok_or_else(|| error("stream has no verified completion"))?;
        if completion.records != records.len() as u64
            || completion.source_bytes != source.len() as u64
            || completion.source_blake3 != blake3::hash(source).to_hex().as_str()
        {
            return Err(error("completion/source mismatch"));
        }
        // Cloning happens after the aggregate retention preflight. The reader's
        // one bounded header is the only temporary header copy during publish.
        let header = reader.header().clone();
        Ok(Self {
            source: Arc::from(source),
            header,
            records,
            completion,
            retained_charge: charge,
            flag_counts: flags,
        })
    }
    pub fn format(&self) -> Format {
        self.completion.format
    }
    pub fn header(&self) -> &Header {
        &self.header
    }
    pub fn record(&self, index: usize) -> Option<&Record> {
        self.records.get(index)
    }
    pub fn len(&self) -> usize {
        self.records.len()
    }
    pub fn is_empty(&self) -> bool {
        self.records.is_empty()
    }
    pub fn retained_charge(&self) -> usize {
        self.retained_charge
    }
    pub fn original_bytes(&self) -> &[u8] {
        &self.source
    }
    pub fn source_blake3(&self) -> &str {
        &self.completion.source_blake3
    }
    pub fn page(&self, start: usize, count: usize) -> Result<&[Record], Error> {
        if count == 0 || count > 100 || start > self.records.len() {
            return Err(Error::Coordinates("read alignment page".into()));
        }
        Ok(&self.records[start..start.saturating_add(count).min(self.records.len())])
    }
    pub fn summary(&self) -> Summary {
        Summary{source_blake3:self.completion.source_blake3.clone(),source_bytes:self.completion.source_bytes,format:match self.completion.format{Format::Sam=>"SAM",Format::Bam=>"BAM"}.into(),records:self.records.len(),flag_unmapped:self.flag_counts[0],flag_secondary:self.flag_counts[1],flag_supplementary:self.flag_counts[2],references:self.header.references().iter().map(|r|ReferenceSummary{name:r.name.clone(),length:r.length,circular:r.circular}).collect(),warnings:self.completion.warnings.clone(),interpretation:"Complete file parsing only. Flags, reference names/build labels, MAPQ and imported matches remain source assertions; no actual-reference checking, unique alignment, coverage, consensus or variant inference is established.".into()}
    }
}

/// Source receipt wraps the actual-reference check; it cannot claim unique
/// mapping or turn an imported report into a privately completed collection.
#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
#[serde(deny_unknown_fields)]
pub struct BoundReferenceCheck {
    pub source_blake3: String,
    pub record_index: usize,
    pub record_blake3: String,
    pub report: crate::read_alignment_reference::Report,
}
impl Collection {
    pub fn check_reference(
        &self,
        index: usize,
        reference: &crate::Molecule,
        options: &crate::read_alignment_reference::Options,
    ) -> Result<BoundReferenceCheck, Error> {
        let record = self
            .record(index)
            .ok_or_else(|| Error::Coordinates("read alignment record index".into()))?;
        let report = crate::read_alignment_reference::check_record(
            &self.header,
            record,
            reference,
            options,
        )?;
        Ok(BoundReferenceCheck {
            source_blake3: self.source_blake3().into(),
            record_index: index,
            record_blake3: record.original_blake3(),
            report,
        })
    }
}

pub struct Export {
    pub format: Format,
    pub bytes: Vec<u8>,
    pub notices: Vec<exchange::ExportNotice>,
}
struct BoundedOutput {
    bytes: Vec<u8>,
    cap: usize,
}
impl std::io::Write for BoundedOutput {
    fn write(&mut self, source: &[u8]) -> std::io::Result<usize> {
        let wanted = self
            .bytes
            .len()
            .checked_add(source.len())
            .filter(|n| *n <= self.cap)
            .ok_or_else(|| std::io::Error::other("alignment export byte limit"))?;
        if wanted > self.bytes.capacity() {
            let capacity = wanted.max(
                self.bytes
                    .capacity()
                    .max(4096)
                    .saturating_mul(2)
                    .min(self.cap),
            );
            self.bytes
                .try_reserve_exact(capacity - self.bytes.len())
                .map_err(|_| std::io::Error::other("alignment export allocation"))?;
        }
        self.bytes.extend_from_slice(source);
        Ok(source.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}
impl Collection {
    /// Convert through the sealed readers/writers; missing sequence/quality and
    /// unsupported cross-format values are never patched to make an export fit.
    /// BAM output is re-framed BGZF; original_bytes() is the exact input download.
    pub fn export(
        &self,
        format: Format,
        reader_limits: exchange::Limits,
        max_bytes: usize,
    ) -> Result<Export, Error> {
        if max_bytes == 0 || max_bytes > 256 << 20 {
            return Err(Error::TooLarge);
        }
        reader_limits.validate().map_err(error)?;
        let sink = BoundedOutput {
            bytes: Vec::new(),
            cap: max_bytes,
        };
        let (sink, notices) = match format {
            Format::Sam => {
                let mut writer =
                    exchange::SamWriter::new(sink, &self.header, reader_limits).map_err(error)?;
                for record in &self.records {
                    writer.write_record(record).map_err(error)?;
                }
                let notices = writer.notices().to_vec();
                (writer.finish().map_err(error)?, notices)
            }
            Format::Bam => {
                let mut writer =
                    exchange::BamWriter::new(sink, &self.header, reader_limits).map_err(error)?;
                for record in &self.records {
                    writer.write_record(record).map_err(error)?;
                }
                let notices = writer.notices().to_vec();
                (writer.finish().map_err(error)?, notices)
            }
        };
        Ok(Export {
            format,
            bytes: sink.bytes,
            notices,
        })
    }
}
