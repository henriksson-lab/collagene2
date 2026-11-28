//! Bounded persistent variant exchange. DTOs are readonly descriptions: decoding
//! one never restores the private completed Collection or authorizes an edit.
use crate::protocol::schema::{AnalysisJob, AnalysisResult};
use crate::variant_reference as reference;
use collagene_model::{Alphabet, Error, Molecule};
use reference::exchange;
use serde::{Deserialize, Serialize};
use std::io::Write;
pub const MAX_SOURCE_BYTES: usize = 16 << 20;
pub const MAX_EXPORT_BYTES: usize = 16 << 20;
const DOMAIN: &str = "variant-worker-1";
macro_rules! dto {($($t:tt)*)=>{#[derive(Clone,Debug,Serialize,Deserialize)]#[serde(deny_unknown_fields)]$($t)*};}
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(rename_all = "lowercase")]
pub enum FileFormat {
    Vcf,
    Bcf,
}
impl FileFormat {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Vcf => "VCF",
            Self::Bcf => "BCF",
        }
    }
    fn exchange(self) -> exchange::Format {
        match self {
            Self::Vcf => exchange::Format::Vcf,
            Self::Bcf => exchange::Format::Bcf,
        }
    }
}
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(rename_all = "lowercase")]
pub enum ExportFormat {
    Original,
    Vcf,
    Bcf,
}
impl ExportFormat {
    fn exchange(self) -> exchange::ExportFormat {
        match self {
            Self::Original => exchange::ExportFormat::Original,
            Self::Vcf => exchange::ExportFormat::Vcf,
            Self::Bcf => exchange::ExportFormat::Bcf,
        }
    }
}
dto! {pub struct ContigSummary {pub name:String,pub length:Option<usize>}}
dto! {pub struct Summary {pub source_hash:String,pub format:FileFormat,pub source_bytes:usize,pub records:usize,pub retained_charge:usize,pub samples:Vec<String>,pub contigs:Vec<ContigSummary>}}
dto! {pub struct Prepared {pub revision:String,pub summary:Summary}}
dto! {pub struct Row {pub index:usize,pub reference_name:String,pub position_1based:Option<usize>,pub vcf:String}}
dto! {pub struct Page {pub revision:String,pub source_hash:String,pub total:usize,pub start:usize,pub count:usize,pub rows:Vec<Row>}}
dto! {pub struct Checked {pub revision:String,pub source_hash:String,pub options:reference::Options,pub report:reference::Report}}
dto! {pub struct Previewed {pub revision:String,pub source_hash:String,pub options:reference::Options,pub preview:reference::Preview}}
dto! {pub struct Exported {pub revision:String,pub source_hash:String,pub format:ExportFormat,pub max_output_bytes:usize,pub bytes:Vec<u8>,pub notices:Vec<String>}}
dto! {pub struct Disposed {pub revision:String,pub source_hash:String}}
#[derive(Default)]
pub struct Session {
    collection: Option<exchange::Collection>,
}
fn bad(s: &str) -> Error {
    Error::Unsupported(format!("variant worker: {s}"))
}
fn hex(s: &str) -> bool {
    s.len() == 64
        && s.bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}
fn source(s: &str) -> Result<(), Error> {
    if hex(s) {
        Ok(())
    } else {
        Err(bad("invalid source hash"))
    }
}
fn binding(a: &str, b: &str) -> Result<(), Error> {
    if a == b {
        Ok(())
    } else {
        Err(bad("request/response context changed"))
    }
}
struct Hash {
    hasher: blake3::Hasher,
    used: usize,
    cap: usize,
}
impl Write for Hash {
    fn write(&mut self, b: &[u8]) -> std::io::Result<usize> {
        self.used = self
            .used
            .checked_add(b.len())
            .filter(|n| *n <= self.cap)
            .ok_or_else(|| std::io::Error::other("bounded variant DTO"))?;
        self.hasher.update(b);
        Ok(b.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}
fn hash(v: &impl Serialize, cap: usize) -> Result<String, Error> {
    let mut w = Hash {
        hasher: blake3::Hasher::new(),
        used: 0,
        cap,
    };
    serde_json::to_writer(&mut w, v).map_err(|_| Error::TooLarge)?;
    Ok(w.hasher.finalize().to_hex().to_string())
}
fn context(v: &impl Serialize) -> Result<String, Error> {
    hash(&(DOMAIN, v), 64 << 10)
}
pub fn browser_limits() -> exchange::Limits {
    exchange::Limits {
        source_bytes: 16 << 20,
        decoded_bytes: 16 << 20,
        records: 50_000,
        retained_bytes: 96 << 20,
        ..Default::default()
    }
}
pub fn mcp_limits() -> exchange::Limits {
    exchange::Limits {
        source_bytes: 4 << 20,
        decoded_bytes: 8 << 20,
        records: 25_000,
        retained_bytes: 48 << 20,
        ..Default::default()
    }
}
fn prepare_context(s: &str, f: FileFormat) -> Result<String, Error> {
    source(s)?;
    context(&(
        "prepare",
        s,
        f,
        "source16-decoded16-records50000-retained96-v1",
    ))
}
pub fn prepare_revision(bytes: &[u8], format: FileFormat) -> Result<String, Error> {
    if bytes.is_empty() || bytes.len() > MAX_SOURCE_BYTES {
        return Err(Error::TooLarge);
    }
    prepare_context(blake3::hash(bytes).to_hex().as_str(), format)
}
pub fn page_revision(s: &str, start: usize, count: usize) -> Result<String, Error> {
    source(s)?;
    if start > 50_000 || count == 0 || count > 100 {
        return Err(bad("invalid page bounds"));
    }
    context(&("page", s, start, count))
}
fn options_preflight(o: &reference::Options) -> Result<(), Error> {
    source(&o.source_blake3)?;
    source(&o.expected_reference_revision)?;
    if o.record_index >= 50_000
        || o.alternate_index == 0
        || o.alternate_index > 64
        || o.reference_name.is_empty()
        || o.reference_name.len() > 1024
        || o.reference_name.chars().any(char::is_control)
    {
        return Err(bad("invalid explicit allele/reference mapping"));
    }
    Ok(())
}
/// Borrowed resource guard before a host clones an actual mapped DNA document.
pub fn reference_preflight(m: &Molecule, o: &reference::Options) -> Result<(), Error> {
    options_preflight(o)?;
    collagene_model::bounds::model_preflight(&[m])?;
    if m.sequence.alphabet() != Alphabet::Dna {
        return Err(bad("mapped reference must be DNA"));
    }
    m.validate()
}
fn reference_context(
    kind: &str,
    s: &str,
    r: &str,
    o: &reference::Options,
) -> Result<String, Error> {
    source(s)?;
    source(r)?;
    options_preflight(o)?;
    binding(s, &o.source_blake3)?;
    binding(r, &o.expected_reference_revision)?;
    context(&(kind, s, r, o))
}
pub fn check_revision_for_revision(
    s: &str,
    r: &str,
    o: &reference::Options,
) -> Result<String, Error> {
    reference_context("check", s, r, o)
}
pub fn preview_revision_for_revision(
    s: &str,
    r: &str,
    o: &reference::Options,
) -> Result<String, Error> {
    reference_context("preview", s, r, o)
}
pub fn check_revision(s: &str, m: &Molecule, o: &reference::Options) -> Result<String, Error> {
    reference_preflight(m, o)?;
    check_revision_for_revision(s, &hash(m, 32 << 20)?, o)
}
pub fn preview_revision(s: &str, m: &Molecule, o: &reference::Options) -> Result<String, Error> {
    reference_preflight(m, o)?;
    preview_revision_for_revision(s, &hash(m, 32 << 20)?, o)
}
pub fn export_revision(s: &str, f: ExportFormat, cap: usize) -> Result<String, Error> {
    source(s)?;
    if cap == 0 || cap > MAX_EXPORT_BYTES {
        return Err(Error::TooLarge);
    }
    context(&("export", s, f, cap))
}
pub fn dispose_revision(s: &str) -> Result<String, Error> {
    source(s)?;
    context(&("dispose", s))
}
pub fn summary(c: &exchange::Collection) -> Result<Summary, Error> {
    if c.header().contigs().len() > 4096 || c.header().sample_names().len() > 256 {
        return Err(Error::TooLarge);
    }
    let s = Summary {
        source_hash: c.source_blake3().into(),
        format: match c.format() {
            exchange::Format::Vcf => FileFormat::Vcf,
            exchange::Format::Bcf => FileFormat::Bcf,
        },
        source_bytes: c.original_bytes().len(),
        records: c.len(),
        retained_charge: c.retained_charge(),
        samples: c.header().sample_names().iter().cloned().collect(),
        contigs: c
            .header()
            .contigs()
            .iter()
            .map(|(n, v)| ContigSummary {
                name: n.clone(),
                length: v.length(),
            })
            .collect(),
    };
    s.validate()?;
    Ok(s)
}
impl Summary {
    pub fn validate(&self) -> Result<(), Error> {
        source(&self.source_hash)?;
        if self.source_bytes == 0
            || self.source_bytes > MAX_SOURCE_BYTES
            || self.records > 50_000
            || self.retained_charge > 96 << 20
            || self.samples.len() > 256
            || self.contigs.len() > 4096
        {
            return Err(Error::TooLarge);
        }
        if self
            .samples
            .iter()
            .chain(self.contigs.iter().map(|c| &c.name))
            .any(|s| s.is_empty() || s.len() > 1024 || s.chars().any(char::is_control))
        {
            return Err(bad("invalid bounded header names"));
        }
        hash(self, 1 << 20)?;
        Ok(())
    }
}
impl Prepared {
    pub fn validate(&self) -> Result<(), Error> {
        self.summary.validate()?;
        binding(
            &self.revision,
            &prepare_context(&self.summary.source_hash, self.summary.format)?,
        )
    }
}
impl Page {
    pub fn validate(&self) -> Result<(), Error> {
        binding(
            &self.revision,
            &page_revision(&self.source_hash, self.start, self.count)?,
        )?;
        if self.total > 50_000
            || self.start > self.total
            || self.rows.len() != self.count.min(self.total - self.start)
        {
            return Err(bad("invalid page shape"));
        }
        let mut bytes = 0usize;
        for (i, row) in self.rows.iter().enumerate() {
            bytes = bytes
                .checked_add(row.vcf.len())
                .and_then(|n| n.checked_add(row.reference_name.len()))
                .ok_or(Error::TooLarge)?;
            if row.index != self.start + i
                || row.reference_name.len() > 1024
                || row.vcf.len() > 256 << 10
                || bytes > 4 << 20
            {
                return Err(Error::TooLarge);
            }
        }
        hash(self, 30 << 20)?;
        Ok(())
    }
}
fn report_validate(r: &reference::Report, s: &str, o: &reference::Options) -> Result<(), Error> {
    options_preflight(o)?;
    binding(s, &o.source_blake3)?;
    binding(&r.source_blake3, s)?;
    binding(&r.actual_reference_revision, &o.expected_reference_revision)?;
    binding(&r.reference_name, &o.reference_name)?;
    binding(&hash(&r.options, 64 << 10)?, &hash(o, 64 << 10)?)?;
    if r.engine_version != reference::ENGINE_VERSION
        || !hex(&r.revision)
        || r.declared_reference.len() > 65_536
        || r.declared_alternative.len() > 65_536
        || r.actual_reference
            .as_ref()
            .is_some_and(|v| v.len() > 65_536)
        || r.interpretation.len() > 64 << 10
        || matches!(&r.outcome,reference::Outcome::Unsupported{reason} if reason.len()>64<<10)
    {
        return Err(Error::TooLarge);
    }
    let revision = hash(
        &(reference::ENGINE_VERSION, o, &r.actual_reference_revision),
        64 << 10,
    )?;
    binding(&r.revision, &revision)?;
    hash(r, 2 << 20)?;
    Ok(())
}
impl Checked {
    pub fn validate(&self) -> Result<(), Error> {
        report_validate(&self.report, &self.source_hash, &self.options)?;
        binding(
            &self.revision,
            &check_revision_for_revision(
                &self.source_hash,
                &self.options.expected_reference_revision,
                &self.options,
            )?,
        )
    }
    pub fn validate_bindings(&self, s: &str, o: &reference::Options) -> Result<(), Error> {
        self.validate()?;
        binding(&self.source_hash, s)?;
        binding(
            &self.revision,
            &check_revision_for_revision(s, &o.expected_reference_revision, o)?,
        )
    }
}
impl Previewed {
    pub fn validate(&self) -> Result<(), Error> {
        report_validate(&self.preview.report, &self.source_hash, &self.options)?;
        binding(
            &self.revision,
            &preview_revision_for_revision(
                &self.source_hash,
                &self.options.expected_reference_revision,
                &self.options,
            )?,
        )?;
        collagene_model::bounds::model_preflight(&[&self.preview.molecule])?;
        self.preview.molecule.validate()?;
        if self.preview.report.outcome != reference::Outcome::ReferenceAlleleMatches
            || self
                .preview
                .molecule
                .metadata
                .get("variant:reference_check")
                != Some(&self.preview.report.revision)
            || self
                .preview
                .molecule
                .metadata
                .get("variant:reference_only")
                .map(String::as_str)
                != Some("true")
        {
            return Err(bad("invalid readonly preview bindings"));
        }
        hash(&self.preview, 30 << 20)?;
        Ok(())
    }
    pub fn validate_bindings(&self, s: &str, o: &reference::Options) -> Result<(), Error> {
        self.validate()?;
        binding(&self.source_hash, s)?;
        binding(
            &self.revision,
            &preview_revision_for_revision(s, &o.expected_reference_revision, o)?,
        )
    }
}
impl Exported {
    pub fn validate(&self) -> Result<(), Error> {
        binding(
            &self.revision,
            &export_revision(&self.source_hash, self.format, self.max_output_bytes)?,
        )?;
        if self.bytes.len() > self.max_output_bytes
            || self.notices.len() > 16
            || self.notices.iter().any(|s| s.len() > 8192)
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
    fn get(&self, s: &str) -> Result<&exchange::Collection, Error> {
        source(s)?;
        let c = self
            .collection
            .as_ref()
            .ok_or_else(|| bad("no completed variant collection prepared"))?;
        binding(c.source_blake3(), s)?;
        Ok(c)
    }
    pub fn process(&mut self, job: AnalysisJob, revision: &str) -> Result<AnalysisResult, Error> {
        match job {
            AnalysisJob::PrepareVariants { bytes, format } => {
                binding(revision, &prepare_revision(&bytes, format)?)?;
                let c = exchange::Collection::import(&bytes, format.exchange(), &browser_limits())
                    .map_err(|e| bad(&e.to_string()))?;
                let result = Prepared {
                    revision: revision.into(),
                    summary: summary(&c)?,
                };
                result.validate()?;
                self.collection = Some(c);
                Ok(AnalysisResult::VariantsPrepared(result))
            }
            AnalysisJob::VariantPage {
                source_hash,
                start,
                count,
            } => {
                binding(revision, &page_revision(&source_hash, start, count)?)?;
                let c = self.get(&source_hash)?;
                let rows = c
                    .page(start, count)
                    .map_err(|e| bad(&e.to_string()))?
                    .into_iter()
                    .map(|r| Row {
                        index: r.index,
                        reference_name: r.reference_name,
                        position_1based: r.position_1based,
                        vcf: r.vcf,
                    })
                    .collect();
                let result = Page {
                    revision: revision.into(),
                    source_hash,
                    total: c.len(),
                    start,
                    count,
                    rows,
                };
                result.validate()?;
                Ok(AnalysisResult::VariantPage(result))
            }
            AnalysisJob::CheckVariantReference {
                source_hash,
                reference,
                options,
            } => {
                binding(
                    revision,
                    &check_revision(&source_hash, &reference, &options)?,
                )?;
                let c = self.get(&source_hash)?;
                let result = Checked {
                    revision: revision.into(),
                    source_hash,
                    options: options.clone(),
                    report: reference::check_record(c, &reference, &options)?,
                };
                result.validate()?;
                Ok(AnalysisResult::VariantReferenceChecked(result))
            }
            AnalysisJob::PreviewVariantAllele {
                source_hash,
                reference,
                options,
            } => {
                binding(
                    revision,
                    &preview_revision(&source_hash, &reference, &options)?,
                )?;
                let c = self.get(&source_hash)?;
                let result = Previewed {
                    revision: revision.into(),
                    source_hash,
                    options: options.clone(),
                    preview: reference::preview(c, &reference, &options)?,
                };
                result.validate()?;
                Ok(AnalysisResult::VariantAllelePreview(result))
            }
            AnalysisJob::ExportVariants {
                source_hash,
                format,
                max_output_bytes,
            } => {
                binding(
                    revision,
                    &export_revision(&source_hash, format, max_output_bytes)?,
                )?;
                let e = self
                    .get(&source_hash)?
                    .export(format.exchange(), max_output_bytes)
                    .map_err(|e| bad(&e.to_string()))?;
                let result = Exported {
                    revision: revision.into(),
                    source_hash,
                    format,
                    max_output_bytes,
                    bytes: e.bytes,
                    notices: e.notices,
                };
                result.validate()?;
                Ok(AnalysisResult::VariantsExported(result))
            }
            AnalysisJob::DisposeVariants { source_hash } => {
                binding(revision, &dispose_revision(&source_hash)?)?;
                self.get(&source_hash)?;
                self.collection = None;
                Ok(AnalysisResult::VariantsDisposed(Disposed {
                    revision: revision.into(),
                    source_hash,
                }))
            }
            _ => Err(bad("job is not a variant collection operation")),
        }
    }
}
