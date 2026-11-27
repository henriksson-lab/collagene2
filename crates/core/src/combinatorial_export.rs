//! Portable archives of exact released native events. Historical manifests never
//! supply current scientific authority. Protocols use one fresh private native
//! reconstruction shared by the existing HTML/PDF renderers. The retained typed
//! finishing route (reference, 5′ homology, or unsealed 3′ annealing) is rebuilt
//! from independent declarations; export never replaces it with a reference merge.
use crate::{
    Error, combinatorial_batch as batch, combinatorial_design as design,
    combinatorial_release::{self as release, Manifest, ReleasedRow},
    formats::{self, TextFormat},
    oligo_formats::Delimiter,
    plan_oligos, run_lineage,
    verification::{self, Status},
    workspace::Workspace,
};
use serde::{Deserialize, Serialize};
use std::{
    collections::BTreeSet,
    io::{self, Cursor, Seek, SeekFrom, Write},
};

pub const MAX_ARCHIVE_BYTES: usize = 120 * 1024 * 1024;
pub const MAX_FILE_BYTES: usize = 16 * 1024 * 1024;
pub const MAX_UNCOMPRESSED_BYTES: usize = 256 * 1024 * 1024;
pub const MAX_FILES: usize = 256 * 6 + 1;
pub const MAX_METADATA_BYTES: usize = 4 * 1024 * 1024;
const CONTRACT: &str = "combinatorial-output-archive-1";

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Options {
    pub genbank: bool,
    pub fasta: bool,
    pub primer_csv: bool,
    pub primer_tsv: bool,
    pub protocol_html: bool,
    pub protocol_pdf: bool,
    /// Membership subset; entry order retains original manifest ordinals.
    #[serde(default)]
    pub construct_ids: Option<Vec<String>>,
}
impl Default for Options {
    fn default() -> Self {
        Self {
            genbank: true,
            fasta: false,
            primer_csv: true,
            primer_tsv: false,
            protocol_html: true,
            protocol_pdf: false,
            construct_ids: None,
        }
    }
}
impl Options {
    pub fn validate(&self) -> Result<(), Error> {
        if !(self.genbank
            || self.fasta
            || self.primer_csv
            || self.primer_tsv
            || self.protocol_html
            || self.protocol_pdf)
        {
            return Err(bad("select at least one output format"));
        }
        if let Some(ids) = &self.construct_ids {
            if ids.is_empty() || ids.len() > 256 {
                return Err(Error::TooLarge);
            }
            let mut seen = BTreeSet::new();
            for id in ids {
                if id.is_empty()
                    || id.len() > 48
                    || !id
                        .bytes()
                        .all(|b| b.is_ascii_alphanumeric() || matches!(b, b'_' | b'-' | b'.'))
                    || !seen.insert(id)
                {
                    return Err(bad("select distinct native construct IDs"));
                }
            }
        }
        Ok(())
    }
    fn protocols(&self) -> bool {
        self.protocol_html || self.protocol_pdf
    }
}
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct FileSummary {
    pub path: String,
    pub bytes: usize,
    pub blake3: String,
    pub mime: String,
    pub construct_id: Option<String>,
    pub warnings: Vec<String>,
}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Archive {
    pub bytes: Vec<u8>,
    pub blake3: String,
    pub mime: String,
    pub extension: String,
    pub context_revision: String,
    pub manifest_hash: String,
    pub files: Vec<FileSummary>,
    pub warnings: Vec<String>,
}
fn bad(message: &str) -> Error {
    Error::Format(format!("Combinatorial export: {message}"))
}
fn digest(value: &str) -> bool {
    value.len() == 64
        && value
            .bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}
fn safe_path(path: &str) -> bool {
    path.is_ascii()
        && !path.contains('\\')
        && !path.contains(':')
        && !path.bytes().any(|b| b.is_ascii_control())
        && (path == "manifest.json" || path.starts_with("constructs/"))
        && path
            .split('/')
            .all(|part| !part.is_empty() && part != "." && part != "..")
}
struct Counter {
    size: usize,
    limit: usize,
}
impl Write for Counter {
    fn write(&mut self, b: &[u8]) -> io::Result<usize> {
        self.size = self
            .size
            .checked_add(b.len())
            .filter(|n| *n <= self.limit)
            .ok_or_else(|| io::Error::other("archive metadata byte limit"))?;
        Ok(b.len())
    }
    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}
fn count(value: &impl Serialize, limit: usize) -> Result<(), Error> {
    serde_json::to_writer(&mut Counter { size: 0, limit }, value).map_err(|_| Error::TooLarge)
}
impl Archive {
    /// Bounded wire/content shape, not a reconstructed scientific verdict.
    pub fn validate(&self) -> Result<(), Error> {
        if self.bytes.len() > MAX_ARCHIVE_BYTES
            || !self.bytes.starts_with(b"PK\x03\x04")
            || self.mime != "application/zip"
            || self.extension != "zip"
            || !digest(&self.blake3)
            || blake3::hash(&self.bytes).to_hex().as_str() != self.blake3
            || !digest(&self.context_revision)
            || !digest(&self.manifest_hash)
            || self.files.is_empty()
            || self.files.len() > MAX_FILES
        {
            return Err(bad("invalid archive shape, digest or limits"));
        }
        let mut paths = BTreeSet::new();
        let mut total = 0usize;
        for file in &self.files {
            if !safe_path(&file.path)
                || !paths.insert(&file.path)
                || file.bytes > MAX_FILE_BYTES
                || !digest(&file.blake3)
                || file.mime.len() > 128
            {
                return Err(bad("invalid archive file summary"));
            }
            total = total
                .checked_add(file.bytes)
                .filter(|n| *n <= MAX_UNCOMPRESSED_BYTES)
                .ok_or(Error::TooLarge)?;
        }
        if !paths.iter().any(|path| path.as_str() == "manifest.json") {
            return Err(bad("archive manifest is missing"));
        }
        count(&(&self.files, &self.warnings), MAX_METADATA_BYTES)?;
        let mut zip = zip::ZipArchive::new(Cursor::new(self.bytes.as_slice()))
            .map_err(|e| bad(&e.to_string()))?;
        if zip.len() != self.files.len() {
            return Err(bad("ZIP entries differ from file summaries"));
        }
        for (index, summary) in self.files.iter().enumerate() {
            let file = zip.by_index(index).map_err(|e| bad(&e.to_string()))?;
            if file.name() != summary.path
                || file.size() != summary.bytes as u64
                || !safe_path(file.name())
            {
                return Err(bad(
                    "ZIP path or size differs from its bounded file summary",
                ));
            }
        }
        Ok(())
    }
}
fn selected(manifest: &Manifest, options: &Options) -> Result<Vec<usize>, Error> {
    match &options.construct_ids {
        None => Ok((0..manifest.rows.len()).collect()),
        Some(ids) => {
            for id in ids {
                if !manifest.rows.iter().any(|row| row.construct_id == *id) {
                    return Err(bad("a requested construct is not in this release manifest"));
                }
            }
            Ok(manifest
                .rows
                .iter()
                .enumerate()
                .filter(|(_, row)| ids.contains(&row.construct_id))
                .map(|(i, _)| i)
                .collect())
        }
    }
}
fn prepare(
    workspace: &Workspace,
    manifest_hash: &str,
    options: &Options,
) -> Result<(Manifest, Vec<usize>), Error> {
    options.validate()?;
    if !digest(manifest_hash) {
        return Err(bad("an exact release manifest hash is required"));
    }
    let manifest = release::load(workspace, manifest_hash)?;
    let indices = selected(&manifest, options)?;
    if options.protocols() {
        // A native run can be internally consistent yet forged away from the
        // original independent declaration. Rebuild declarations before science.
        // The same fresh independent-declaration context supplies actual strand
        // costs. Do not replay origins a second time merely to estimate work.
        let costs = manifest.validate_current_declarations(workspace)?;
        let mut work = 0usize;
        for &index in &indices {
            let id = &manifest.rows[index].construct_id;
            let charge = *costs
                .get(id)
                .ok_or_else(|| bad("released declaration cost unavailable"))?;
            work=work.checked_add(charge).filter(|n|*n<=batch::MAX_BATCH_WORK).ok_or_else(||bad("aggregate protocol reconstruction work limit reached; select fewer constructs"))?;
        }
    }
    Ok((manifest, indices))
}
pub fn preflight(
    workspace: &Workspace,
    manifest_hash: &str,
    options: &Options,
) -> Result<(), Error> {
    prepare(workspace, manifest_hash, options).map(|_| ())
}
fn context(workspace: &Workspace, manifest_hash: &str, options: &Options) -> Result<String, Error> {
    crate::assembly_recipe::hash(&(
        CONTRACT,
        verification::ENGINE_VERSION,
        run_lineage::workspace_fingerprint(workspace)?,
        manifest_hash,
        options,
    ))
}
pub fn revision(
    workspace: &Workspace,
    manifest_hash: &str,
    options: &Options,
) -> Result<String, Error> {
    preflight(workspace, manifest_hash, options)?;
    context(workspace, manifest_hash, options)
}

struct BoundedCursor(Cursor<Vec<u8>>);
impl Write for BoundedCursor {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        if self
            .0
            .position()
            .checked_add(bytes.len() as u64)
            .is_none_or(|n| n > MAX_ARCHIVE_BYTES as u64)
        {
            return Err(io::Error::other(
                "ZIP exceeds120MiB; choose fewer rows or output formats",
            ));
        }
        self.0.write(bytes)
    }
    fn flush(&mut self) -> io::Result<()> {
        self.0.flush()
    }
}
impl Seek for BoundedCursor {
    fn seek(&mut self, position: SeekFrom) -> io::Result<u64> {
        let n = self.0.seek(position)?;
        if n > MAX_ARCHIVE_BYTES as u64 {
            return Err(io::Error::other("ZIP seek exceeds archive limit"));
        }
        Ok(n)
    }
}
struct Builder {
    zip: zip::ZipWriter<BoundedCursor>,
    files: Vec<FileSummary>,
    uncompressed: usize,
}
impl Builder {
    fn add(
        &mut self,
        path: String,
        bytes: &[u8],
        mime: &str,
        construct_id: Option<&str>,
        warnings: Vec<String>,
    ) -> Result<(), Error> {
        if !safe_path(&path) || self.files.iter().any(|f| f.path == path) {
            return Err(bad("unsafe or duplicate archive path"));
        }
        if bytes.len() > MAX_FILE_BYTES || self.files.len() >= MAX_FILES {
            return Err(bad(
                "file exceeds16MiB or archive file-count limit; choose fewer formats/constructs",
            ));
        }
        self.uncompressed = self
            .uncompressed
            .checked_add(bytes.len())
            .filter(|n| *n <= MAX_UNCOMPRESSED_BYTES)
            .ok_or_else(|| bad("archive exceeds256MiB uncompressed; select fewer constructs"))?;
        let summary = FileSummary {
            path: path.clone(),
            bytes: bytes.len(),
            blake3: blake3::hash(bytes).to_hex().to_string(),
            mime: mime.into(),
            construct_id: construct_id.map(str::to_owned),
            warnings,
        };
        count(&(&self.files, &summary), MAX_METADATA_BYTES)?;
        self.zip
            .start_file(
                path,
                zip::write::SimpleFileOptions::default()
                    .compression_method(zip::CompressionMethod::Deflated),
            )
            .map_err(|e| bad(&e.to_string()))?;
        self.zip.write_all(bytes).map_err(|e| bad(&e.to_string()))?;
        self.files.push(summary);
        Ok(())
    }
}
#[derive(Serialize)]
struct ExportedRow {
    released: ReleasedRow,
    fresh_status: Option<Status>,
    fresh_engine: Option<String>,
    data_state: &'static str,
}
#[derive(Serialize)]
struct ExportManifest<'a> {
    schema_version: u16,
    release_manifest_hash: &'a str,
    release_context_revision: &'a str,
    release_recorded_engine: &'a str,
    source_workspace_fingerprint: &'a str,
    export_context_revision: &'a str,
    options: &'a Options,
    rows: &'a [ExportedRow],
    files: &'a [FileSummary],
    manifest_self_entry: &'static str,
    warnings: &'a [String],
}
struct LimitedJson(Vec<u8>);
impl Write for LimitedJson {
    fn write(&mut self, b: &[u8]) -> io::Result<usize> {
        if self
            .0
            .len()
            .checked_add(b.len())
            .is_none_or(|n| n > MAX_FILE_BYTES)
        {
            return Err(io::Error::other("export manifest exceeds16MiB"));
        }
        self.0.extend_from_slice(b);
        Ok(b.len())
    }
    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}

pub fn export(
    workspace: &Workspace,
    manifest_hash: &str,
    options: &Options,
) -> Result<Archive, Error> {
    let (manifest, indices) = prepare(workspace, manifest_hash, options)?;
    let context_revision = context(workspace, manifest_hash, options)?;
    let mut out = Builder {
        zip: zip::ZipWriter::new(BoundedCursor(Cursor::new(Vec::new()))),
        files: vec![],
        uncompressed: 0,
    };
    let mut rows = vec![];
    let warnings=vec![if options.protocols(){"Protocols use fresh native reconstruction of each selected independent declaration, retained plan and target. Named sequence events remain historical native outputs; current computational consistency is not experimental success.".into()}else{"Historical sequence/declared-primer exchange only. No current assembly computation was performed, and imported release records do not restore scientific authority.".into()},"Every requested file is included or archive generation fails. Resource ceilings:120MiB ZIP,256MiB uncompressed,16MiB per file, plus inherited protocol PDF text/page bounds. Choose a construct subset or fewer formats for large releases.".into()];
    for index in indices {
        let row = &manifest.rows[index];
        let product = design::resolve_event(workspace, &row.output.construct)?;
        let target = design::resolve_event(workspace, &row.target)?;
        let asset = workspace
            .attachments
            .get(&row.output.run.attachment_hash)
            .ok_or_else(|| bad("released run attachment unavailable"))?;
        let run = run_lineage::load_run_record(&asset.bytes)?;
        let checked = options
            .protocols()
            .then(|| verification::checked_protocol(&run.plan, target));
        if let Some(checked) = &checked {
            let report = checked.report();
            if report.engine_version != verification::ENGINE_VERSION
                || report.plan_revision != run.plan_revision
                || report.target_revision != row.target.revision
                || report.status == Status::Invalid
                || (report.status == Status::Incomplete && !row.options.allow_incomplete)
                || report.products.len() != run.products.len()
                || report.products.is_empty()
                || report
                    .products
                    .iter()
                    .zip(&run.products)
                    .any(|(m, b)| m.revision() != b.revision)
                || !report
                    .products
                    .last()
                    .is_some_and(|m| verification::plan_target_matches(&run.plan, m, target))
            {
                return Err(bad(
                    "fresh reconstruction differs from this retained release or independent target; no archive was published",
                ));
            }
        }
        for (enabled, format, extension, mime) in [
            (
                options.genbank,
                TextFormat::GenBank,
                "gb",
                "text/plain;charset=utf-8",
            ),
            (
                options.fasta,
                TextFormat::Fasta,
                "fasta",
                "text/plain;charset=utf-8",
            ),
        ] {
            if enabled {
                let exported = formats::export_with_report(product, format)?;
                out.add(
                    release::archive_entry(row.ordinal, &row.construct_id, extension)?,
                    exported.text.as_bytes(),
                    mime,
                    Some(&row.construct_id),
                    exported.warnings,
                )?;
            }
        }
        for (enabled, delimiter, extension) in [
            (options.primer_csv, Delimiter::Comma, "csv"),
            (options.primer_tsv, Delimiter::Tab, "tsv"),
        ] {
            if enabled {
                let exported = plan_oligos::export(
                    &run.plan,
                    target,
                    &plan_oligos::Options::default(),
                    delimiter,
                )?;
                if !exported.report.complete {
                    return Err(bad(
                        "complete declared-primer export is unavailable for a selected row; no partial sheet was substituted",
                    ));
                }
                let mut losses = exported.export.warnings;
                losses.extend(exported.report.warnings);
                if exported.report.oligos.records.is_empty() {
                    losses.push("No oligos are declared by this plan; this sheet contains column headers only.".into());
                }
                out.add(
                    release::archive_entry(row.ordinal, &row.construct_id, extension)?,
                    exported.export.text.as_bytes(),
                    "text/plain;charset=utf-8",
                    Some(&row.construct_id),
                    losses,
                )?;
            }
        }
        if options.protocol_html {
            let text = verification::export_protocol_html_checked(
                checked.as_ref().expect("protocol context prepared"),
            );
            out.add(
                release::archive_entry(row.ordinal, &row.construct_id, "html")?,
                text.as_bytes(),
                "text/html;charset=utf-8",
                Some(&row.construct_id),
                vec![],
            )?;
        }
        if options.protocol_pdf {
            let pdf = crate::pdf::export_protocol_pdf_checked(
                checked.as_ref().expect("protocol context prepared"),
            )?;
            out.add(
                release::archive_entry(row.ordinal, &row.construct_id, "pdf")?,
                &pdf.bytes,
                "application/pdf",
                Some(&row.construct_id),
                pdf.warnings,
            )?;
        }
        rows.push(ExportedRow {
            released: row.clone(),
            fresh_status: checked.as_ref().map(|c| c.report().status),
            fresh_engine: checked.as_ref().map(|c| c.report().engine_version.clone()),
            data_state: "exact historical native output events",
        });
    }
    let metadata = ExportManifest {
        schema_version: 1,
        release_manifest_hash: manifest_hash,
        release_context_revision: &manifest.context_revision,
        release_recorded_engine: &manifest.engine_version,
        source_workspace_fingerprint: &manifest.source_workspace_fingerprint,
        export_context_revision: &context_revision,
        options,
        rows: &rows,
        files: &out.files,
        manifest_self_entry: "manifest.json",
        warnings: &warnings,
    };
    let mut encoded = LimitedJson(vec![]);
    serde_json::to_writer_pretty(&mut encoded, &metadata).map_err(|e| bad(&e.to_string()))?;
    out.add(
        "manifest.json".into(),
        &encoded.0,
        "application/json",
        None,
        vec![],
    )?;
    let bytes = out
        .zip
        .finish()
        .map_err(|e| bad(&e.to_string()))?
        .0
        .into_inner();
    let archive = Archive {
        blake3: blake3::hash(&bytes).to_hex().to_string(),
        bytes,
        mime: "application/zip".into(),
        extension: "zip".into(),
        context_revision,
        manifest_hash: manifest_hash.into(),
        files: out.files,
        warnings,
    };
    archive.validate()?;
    Ok(archive)
}
