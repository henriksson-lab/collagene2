//! Persistent bounded browser-worker alignment collections. Decoded DTOs are
//! historical descriptions; only the private Collection grants record access.
use crate::{
    Error, Molecule,
    protocol::{AnalysisJob, AnalysisResult},
    read_alignment_reference as reference, read_alignment_store as store,
};
use collagene_read_alignments::{AuxValue, Record};
use serde::{Deserialize, Serialize};
use std::{collections::BTreeSet, io::Write};

pub const MAX_SOURCE_BYTES: usize = 16 * 1024 * 1024;
pub const MAX_EXPORT_BYTES: usize = 16 * 1024 * 1024;
const MAX_RETAINED_BYTES: usize = 64 * 1024 * 1024;
const MAX_RECORDS: usize = 100_000;
const MAX_PAGE_BYTES: usize = 2 * 1024 * 1024;
const DOMAIN: &str = "read-alignment-worker-1";
macro_rules! dto {($($t:tt)*)=>{#[derive(Clone,Debug,Serialize,Deserialize)]#[serde(deny_unknown_fields)]$($t)*};}
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum FileFormat {
    Sam,
    Bam,
}
impl FileFormat {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Sam => "SAM",
            Self::Bam => "BAM",
        }
    }
    fn exchange(self) -> store::Format {
        match self {
            Self::Sam => store::Format::Sam,
            Self::Bam => store::Format::Bam,
        }
    }
    fn label(value: &str) -> Result<Self, Error> {
        match value {
            "SAM" => Ok(Self::Sam),
            "BAM" => Ok(Self::Bam),
            _ => Err(bad("invalid alignment format")),
        }
    }
}
dto! {pub struct Prepared {pub revision:String,pub summary:store::Summary}}
dto! {pub struct CigarPreview {pub length:u32,pub operation:char}}
dto! {pub struct TagPreview {pub tag:[u8;2],pub wire_type:u8,pub value_summary:String,pub truncated:bool}}
dto! {#[derive(Default)] pub struct Truncated {pub name:bool,pub reference:bool,pub cigar:bool,pub sequence:bool,pub qualities:bool,pub tags:bool}}
dto! {pub struct RecordPreview {
    pub index:usize,pub original_blake3:String,pub name:Option<String>,pub flags:u16,
    pub reference:Option<String>,pub start:Option<u32>,pub mapping_quality:Option<u8>,
    pub cigar:Vec<CigarPreview>,pub cigar_count:usize,
    pub sequence:Option<String>,pub sequence_length:Option<usize>,
    pub qualities:Option<Vec<u8>>,pub qualities_length:Option<usize>,
    pub tags:Vec<TagPreview>,pub tag_count:usize,pub truncated:Truncated,
}}
dto! {pub struct Page {pub revision:String,pub source_hash:String,pub total:usize,pub start:usize,pub count:usize,pub records:Vec<RecordPreview>}}
dto! {pub struct Checked {pub revision:String,pub options:reference::Options,pub check:store::BoundReferenceCheck}}
dto! {pub struct Exported {pub revision:String,pub source_hash:String,pub format:FileFormat,pub max_output_bytes:usize,pub bytes:Vec<u8>,pub notices:Vec<String>}}
dto! {pub struct Disposed {pub revision:String,pub source_hash:String}}
#[derive(Default)]
pub struct Session {
    collection: Option<store::Collection>,
}
fn bad(message: &str) -> Error {
    Error::Unsupported(format!("read alignment worker: {message}"))
}
fn hex(value: &str) -> bool {
    value.len() == 64
        && value
            .bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}
fn source_hash(value: &str) -> Result<(), Error> {
    if !hex(value) {
        Err(bad("invalid source hash"))
    } else {
        Ok(())
    }
}
struct Hash {
    hash: blake3::Hasher,
    used: usize,
    cap: usize,
}
impl Write for Hash {
    fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
        self.used = self
            .used
            .checked_add(bytes.len())
            .filter(|&n| n <= self.cap)
            .ok_or_else(|| std::io::Error::other("bounded alignment DTO"))?;
        self.hash.update(bytes);
        Ok(bytes.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}
fn hash(value: &impl Serialize, cap: usize) -> Result<String, Error> {
    let mut out = Hash {
        hash: blake3::Hasher::new(),
        used: 0,
        cap,
    };
    serde_json::to_writer(&mut out, value).map_err(|_| Error::TooLarge)?;
    Ok(out.hash.finalize().to_hex().to_string())
}
fn context(value: &impl Serialize) -> Result<String, Error> {
    hash(&(DOMAIN, value), 64 * 1024)
}
fn prepare_context(source: &str, format: FileFormat) -> Result<String, Error> {
    source_hash(source)?;
    context(&(
        "prepare",
        source,
        format,
        "source16MiB-retained64MiB-profile1",
    ))
}
pub fn prepare_revision(bytes: &[u8], format: FileFormat) -> Result<String, Error> {
    if bytes.is_empty() || bytes.len() > MAX_SOURCE_BYTES {
        return Err(Error::TooLarge);
    }
    prepare_context(blake3::hash(bytes).to_hex().as_str(), format)
}
pub fn page_revision(source: &str, start: usize, count: usize) -> Result<String, Error> {
    source_hash(source)?;
    if start > MAX_RECORDS || count == 0 || count > 100 {
        return Err(bad("invalid page bounds"));
    }
    context(&("page", source, start, count))
}
/// Cheap borrowed native-model guard before the UI clones a reference. Exact
/// canonical JSON/hash and wire bounds are independently enforced in the worker.
pub fn reference_preflight(molecule: &Molecule) -> Result<(), Error> {
    collagene_model::bounds::model_preflight(&[molecule])?;
    if molecule.sequence.len() > crate::model::MAX_SEQUENCE_LEN {
        return Err(Error::TooLarge);
    }
    molecule.validate()
}
pub fn check_revision_for_revision(
    source: &str,
    index: usize,
    reference_revision: &str,
    options: &reference::Options,
) -> Result<String, Error> {
    source_hash(source)?;
    options.validate()?;
    if index >= MAX_RECORDS || !hex(reference_revision) {
        return Err(bad("invalid reference context"));
    }
    context(&("check", source, index, reference_revision, options))
}
pub fn check_revision(
    source: &str,
    index: usize,
    molecule: &Molecule,
    options: &reference::Options,
) -> Result<String, Error> {
    reference_preflight(molecule)?;
    check_revision_for_revision(source, index, &hash(molecule, 32 * 1024 * 1024)?, options)
}
pub fn export_revision(
    source: &str,
    format: FileFormat,
    max_output_bytes: usize,
) -> Result<String, Error> {
    source_hash(source)?;
    if max_output_bytes == 0 || max_output_bytes > MAX_EXPORT_BYTES {
        return Err(Error::TooLarge);
    }
    context(&("export", source, format, max_output_bytes))
}
pub fn dispose_revision(source: &str) -> Result<String, Error> {
    source_hash(source)?;
    context(&("dispose", source))
}
fn binding(actual: &str, expected: &str) -> Result<(), Error> {
    if actual != expected {
        Err(bad("request/response context changed"))
    } else {
        Ok(())
    }
}
fn limits() -> store::ImportLimits {
    let mut limits = store::ImportLimits::default();
    limits.max_source_bytes = MAX_SOURCE_BYTES;
    limits.max_retained_bytes = MAX_RETAINED_BYTES;
    limits.max_records = MAX_RECORDS;
    limits.reader.max_input_bytes = MAX_SOURCE_BYTES as u64;
    limits.reader.max_uncompressed_bytes = MAX_RETAINED_BYTES as u64;
    limits
}
fn prefix(value: &str, cap: usize) -> (String, bool) {
    let mut end = value.len().min(cap);
    while !value.is_char_boundary(end) {
        end -= 1;
    }
    (value[..end].to_owned(), end < value.len())
}
fn optional_prefix(value: &Option<String>, cap: usize) -> (Option<String>, bool) {
    match value {
        Some(v) => {
            let (v, t) = prefix(v, cap);
            (Some(v), t)
        }
        None => (None, false),
    }
}
fn tag_preview(aux: &collagene_read_alignments::Aux) -> TagPreview {
    let (value_summary, truncated) = match &aux.value {
        AuxValue::Char(v) => (format!("character byte {v:02x}"), false),
        AuxValue::Integer(v) => (format!("integer {v}"), false),
        AuxValue::FloatBits(v) => (format!("float bits {v:08x}"), false),
        AuxValue::String(v) | AuxValue::Hex(v) => {
            let hex = v
                .iter()
                .take(64)
                .map(|b| format!("{b:02x}"))
                .collect::<String>();
            (
                format!("{} source bytes; hexadecimal prefix {hex}", v.len()),
                v.len() > 64,
            )
        }
        AuxValue::Array {
            element_type,
            values,
        } => {
            let bits = values
                .iter()
                .take(8)
                .map(|v| format!("{v:08x}"))
                .collect::<Vec<_>>()
                .join(" ");
            (
                format!(
                    "array type {}: {} values; raw bits {bits}",
                    *element_type as char,
                    values.len()
                ),
                values.len() > 8,
            )
        }
    };
    TagPreview {
        tag: aux.tag,
        wire_type: aux.wire_type,
        value_summary,
        truncated,
    }
}
fn preview(index: usize, record: &Record) -> RecordPreview {
    let f = record.fields();
    let (name, name_cut) = optional_prefix(&f.name, 256);
    let (reference, ref_cut) = optional_prefix(&f.reference, 256);
    let (sequence, seq_cut) = optional_prefix(&f.sequence, 200);
    RecordPreview {
        index,
        original_blake3: record.original_blake3(),
        name,
        flags: f.flags,
        reference,
        start: f.start,
        mapping_quality: f.mapping_quality,
        cigar: f
            .cigar
            .iter()
            .take(64)
            .map(|op| CigarPreview {
                length: op.len,
                operation: op.kind.symbol(),
            })
            .collect(),
        cigar_count: f.cigar.len(),
        sequence,
        sequence_length: f.sequence.as_ref().map(String::len),
        qualities: f
            .qualities
            .as_ref()
            .map(|v| v.iter().take(200).copied().collect()),
        qualities_length: f.qualities.as_ref().map(Vec::len),
        tags: f.aux.iter().take(32).map(tag_preview).collect(),
        tag_count: f.aux.len(),
        truncated: Truncated {
            name: name_cut,
            reference: ref_cut,
            cigar: f.cigar.len() > 64,
            sequence: seq_cut,
            qualities: f.qualities.as_ref().is_some_and(|v| v.len() > 200),
            tags: f.aux.len() > 32,
        },
    }
}
fn validate_summary(summary: &store::Summary) -> Result<(), Error> {
    source_hash(&summary.source_blake3)?;
    FileFormat::label(&summary.format)?;
    if summary.source_bytes == 0
        || summary.source_bytes > MAX_SOURCE_BYTES as u64
        || summary.records > MAX_RECORDS
        || summary.flag_unmapped > summary.records
        || summary.flag_secondary > summary.records
        || summary.flag_supplementary > summary.records
        || summary.references.len() > 4096
        || summary.warnings.len() > 16
        || summary.warnings.iter().any(|s| s.len() > 4096)
        || summary.interpretation.len() > 4096
    {
        return Err(Error::TooLarge);
    }
    let mut names = BTreeSet::new();
    let mut bytes = 0usize;
    for reference in &summary.references {
        bytes = bytes
            .checked_add(reference.name.len())
            .ok_or(Error::TooLarge)?;
        if reference.name.is_empty()
            || reference.length == 0
            || bytes > 2 * 1024 * 1024
            || !names.insert(reference.name.as_str())
        {
            return Err(bad("invalid reference summary"));
        }
    }
    hash(summary, 4 * 1024 * 1024)?;
    Ok(())
}
impl Prepared {
    pub fn validate(&self) -> Result<(), Error> {
        validate_summary(&self.summary)?;
        binding(
            &self.revision,
            &prepare_context(
                &self.summary.source_blake3,
                FileFormat::label(&self.summary.format)?,
            )?,
        )
    }
}
impl RecordPreview {
    fn validate(&self) -> Result<(), Error> {
        source_hash(&self.original_blake3)?;
        if self.index >= MAX_RECORDS
            || self.name.as_ref().is_some_and(|s| s.len() > 256)
            || self.reference.as_ref().is_some_and(|s| s.len() > 256)
            || self.cigar_count > 16_384
            || self.cigar.len() != self.cigar_count.min(64)
            || self.truncated.cigar != (self.cigar_count > 64)
            || self
                .cigar
                .iter()
                .any(|op| op.length == 0 || !"MIDNSHP=X".contains(op.operation))
            || self.sequence.as_ref().is_some_and(|s| {
                s.len() > 200
                    || !s
                        .bytes()
                        .all(|b| b.is_ascii_alphabetic() || b == b'=' || b == b'.')
            })
            || self.sequence_length.is_some_and(|n| n == 0 || n > 250_000)
            || self
                .qualities
                .as_ref()
                .is_some_and(|v| v.len() > 200 || v.iter().any(|q| *q > 93))
            || self.tag_count > 128
            || self.tags.len() != self.tag_count.min(32)
            || self.truncated.tags != (self.tag_count > 32)
            || self.tags.iter().any(|t| t.value_summary.len() > 256)
            || self.qualities_length.is_some() && self.qualities_length != self.sequence_length
        {
            return Err(bad("invalid record preview"));
        }
        for (preview, length, truncated) in [
            (
                self.sequence.as_ref().map(String::len),
                self.sequence_length,
                self.truncated.sequence,
            ),
            (
                self.qualities.as_ref().map(Vec::len),
                self.qualities_length,
                self.truncated.qualities,
            ),
        ] {
            if preview != length.map(|n| n.min(200)) || truncated != length.is_some_and(|n| n > 200)
            {
                return Err(bad("preview presence/truncation mismatch"));
            }
        }
        Ok(())
    }
}
impl Page {
    pub fn validate(&self) -> Result<(), Error> {
        binding(
            &self.revision,
            &page_revision(&self.source_hash, self.start, self.count)?,
        )?;
        if self.total > MAX_RECORDS
            || self.start > self.total
            || self.records.len() != self.count.min(self.total - self.start)
        {
            return Err(bad("invalid page shape"));
        }
        for (offset, record) in self.records.iter().enumerate() {
            record.validate()?;
            if record.index != self.start + offset {
                return Err(bad("page record index mismatch"));
            }
        }
        hash(self, MAX_PAGE_BYTES)?;
        Ok(())
    }
}
impl Checked {
    pub fn validate(&self) -> Result<(), Error> {
        let check = &self.check;
        let report = &check.report;
        source_hash(&check.source_blake3)?;
        source_hash(&check.record_blake3)?;
        for digest in [
            &report.revision,
            &report.binding.reference_revision,
            &report.binding.header_revision,
            &report.binding.record_original_blake3,
            &report.binding.options_revision,
        ] {
            source_hash(digest)?;
        }
        self.options.validate()?;
        binding(
            &self.revision,
            &check_revision_for_revision(
                &check.source_blake3,
                check.record_index,
                &report.binding.reference_revision,
                &self.options,
            )?,
        )?;
        if report.engine_version != reference::ENGINE_VERSION
            || check.record_blake3 != report.binding.record_original_blake3
            || report.binding.reference_name != self.options.reference_name
            || report.binding.options_revision != hash(&self.options, 16384)?
            || report.revision != hash(&(reference::ENGINE_VERSION, &report.binding), 16384)?
            || !matches!(report.binding.format.as_str(), "SAM" | "BAM")
            || report.runs.len() > self.options.max_runs
            || report.columns.len() > self.options.max_columns
            || report.incomplete.len() > 32
            || report.complete != report.incomplete.is_empty()
            || report.complete && !report.reference_bound
            || !report.reference_bound && (!report.columns.is_empty() || !report.runs.is_empty())
            || report.reference_length == 0
            || report.reference_length > crate::model::MAX_SEQUENCE_LEN
            || report.reverse != (report.reported_flags & 16 != 0)
            || report.warnings.len() > 16
            || report.warnings.iter().any(|s| s.len() > 4096)
        {
            return Err(bad("invalid bound reference report"));
        }
        let mut column_end = 0;
        for run in &report.runs {
            if run.query_start > run.query_end
                || run.reference_start > run.reference_end
                || run.column_start != column_end
                || run.column_end < run.column_start
                || run.column_end > report.columns.len()
            {
                return Err(bad("invalid checked CIGAR run"));
            }
            column_end = run.column_end;
        }
        if column_end != report.columns.len() {
            return Err(bad("unbound checked columns"));
        }
        for (index, column) in report.columns.iter().enumerate() {
            if column.query_index != index as u64
                || column
                    .reference_position
                    .is_some_and(|p| p >= report.reference_length)
                || column.raw_phred.is_some_and(|q| q > 93)
            {
                return Err(bad("invalid checked call coordinates"));
            }
        }
        hash(report, self.options.max_report_bytes)?;
        Ok(())
    }
    /// Cheap current-context checks only; never replay calls or import authority.
    pub fn validate_bindings(
        &self,
        expected_revision: &str,
        source: &str,
        index: usize,
        record_hash: &str,
        reference_revision: &str,
        options: &reference::Options,
    ) -> Result<(), Error> {
        // Full readonly shape validation already runs in protocol::decode_response.
        binding(&self.revision, expected_revision)?;
        binding(
            &self.revision,
            &check_revision_for_revision(source, index, reference_revision, options)?,
        )?;
        if self.check.source_blake3 != source
            || self.check.record_index != index
            || self.check.record_blake3 != record_hash
            || self.check.report.binding.record_original_blake3 != record_hash
            || self.check.report.binding.reference_revision != reference_revision
            || self.options != *options
        {
            return Err(bad("checked source/reference context mismatch"));
        }
        Ok(())
    }
}
impl Exported {
    pub fn validate(&self) -> Result<(), Error> {
        binding(
            &self.revision,
            &export_revision(&self.source_hash, self.format, self.max_output_bytes)?,
        )?;
        if self.bytes.is_empty()
            || self.bytes.len() > self.max_output_bytes
            || self.notices.len() > 32
            || self.notices.iter().any(|s| s.len() > 512)
        {
            return Err(Error::TooLarge);
        }
        Ok(())
    }
}
impl Disposed {
    pub fn validate(&self) -> Result<(), Error> {
        binding(&self.revision, &dispose_revision(&self.source_hash)?)
    }
}
impl Session {
    fn collection(&self, source: &str) -> Result<&store::Collection, Error> {
        source_hash(source)?;
        self.collection
            .as_ref()
            .filter(|c| c.source_blake3() == source)
            .ok_or_else(|| bad("complete source is not retained in this worker"))
    }
    pub fn process(&mut self, job: AnalysisJob, revision: &str) -> Result<AnalysisResult, Error> {
        match job {
            AnalysisJob::PrepareReadAlignments { bytes, format } => {
                binding(revision, &prepare_revision(&bytes, format)?)?;
                let collection = store::Collection::import(&bytes, format.exchange(), &limits())?;
                let result = Prepared {
                    revision: revision.into(),
                    summary: collection.summary(),
                };
                result.validate()?;
                // Publish atomically only after complete parsing and bounded summary.
                self.collection = Some(collection);
                Ok(AnalysisResult::ReadAlignmentsPrepared(result))
            }
            AnalysisJob::ReadAlignmentPage {
                source_hash,
                start,
                count,
            } => {
                binding(revision, &page_revision(&source_hash, start, count)?)?;
                let collection = self.collection(&source_hash)?;
                let records = collection
                    .page(start, count)?
                    .iter()
                    .enumerate()
                    .map(|(offset, record)| preview(start + offset, record))
                    .collect();
                let result = Page {
                    revision: revision.into(),
                    source_hash,
                    total: collection.len(),
                    start,
                    count,
                    records,
                };
                result.validate()?;
                Ok(AnalysisResult::ReadAlignmentPage(result))
            }
            AnalysisJob::CheckReadAlignment {
                source_hash,
                index,
                reference,
                options,
            } => {
                binding(
                    revision,
                    &check_revision(&source_hash, index, &reference, &options)?,
                )?;
                let check = self
                    .collection(&source_hash)?
                    .check_reference(index, &reference, &options)?;
                let result = Checked {
                    revision: revision.into(),
                    options,
                    check,
                };
                result.validate()?;
                Ok(AnalysisResult::ReadAlignmentChecked(result))
            }
            AnalysisJob::ExportReadAlignments {
                source_hash,
                format,
                max_output_bytes,
            } => {
                binding(
                    revision,
                    &export_revision(&source_hash, format, max_output_bytes)?,
                )?;
                let exported = self.collection(&source_hash)?.export(
                    format.exchange(),
                    limits().reader,
                    max_output_bytes,
                )?;
                let result = Exported {
                    revision: revision.into(),
                    source_hash,
                    format,
                    max_output_bytes,
                    bytes: exported.bytes,
                    notices: exported.notices.iter().map(|n| format!("{n:?}")).collect(),
                };
                result.validate()?;
                Ok(AnalysisResult::ReadAlignmentsExported(result))
            }
            AnalysisJob::DisposeReadAlignments { source_hash } => {
                binding(revision, &dispose_revision(&source_hash)?)?;
                self.collection(&source_hash)?;
                let result = Disposed {
                    revision: revision.into(),
                    source_hash,
                };
                result.validate()?;
                self.collection = None;
                Ok(AnalysisResult::ReadAlignmentsDisposed(result))
            }
            _ => Err(bad("job is not an alignment session operation")),
        }
    }
}
