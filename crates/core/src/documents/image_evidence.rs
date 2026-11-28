//! Immutable documentary image attachments. Images never alter sequence, features,
//! material origin, history events, experimental status or verification authority.
//! Codec work belongs in the private compute worker, never UI/archive validation.
use crate::documents::workspace::{MAX_ATTACHMENTS, Workspace, WorkspaceAttachment};
use collagene_image_evidence::{Evidence, Limits, PREVIEW_ENGINE_VERSION};
use collagene_model::Error;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

pub const MANIFEST_KIND: &str = "image-evidence-record-v1";
/// One binary kind permits exact-byte deduplication when source PNG equals preview.
/// Original/preview roles are distinct, explicit fields of the immutable manifest.
pub const IMAGE_KIND: &str = "image-evidence-image-v1";
pub const VERSION: u16 = 1;
pub const MAX_RECORD_BYTES: usize = 16 * 1024;
pub const MAX_IMAGES: usize = 100;
pub const MAX_SOURCE_BYTES: usize = 16 * 1024 * 1024;
pub const MAX_PREVIEW_BYTES: usize = 5 * 1024 * 1024;

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Binding {
    pub document_id: String,
    pub revision: String,
}
impl Binding {
    /// Identity of one explicit source/name/association preparation, without
    /// transferring unrelated workspace documents to the codec worker.
    pub fn preparation_revision(&self, name: &str, original_hash: &str) -> Result<String, Error> {
        self.validate_fields()?;
        if !name_valid(name) || !digest(original_hash) {
            return Err(bad("invalid preparation name/source hash"));
        }
        identity(&("documentary-image-prepare-v1", name, self, original_hash))
    }
    fn validate_fields(&self) -> Result<(), Error> {
        if self.document_id.is_empty()
            || self.document_id.len() > 128
            || !self
                .document_id
                .bytes()
                .all(|b| b.is_ascii_alphanumeric() || matches!(b, b'_' | b'-' | b'.' | b':'))
            || !digest(&self.revision)
        {
            return Err(bad("invalid associated document/revision"));
        }
        Ok(())
    }
    pub fn validate_retained(&self, workspace: &Workspace) -> Result<(), Error> {
        self.validate_fields()?;
        let document = workspace
            .documents
            .iter()
            .find(|d| d.id == self.document_id)
            .ok_or_else(|| bad("associated document is unavailable"))?;
        let belongs_to_document = document.head == self.revision
            || document
                .undo
                .iter()
                .chain(&document.redo)
                .any(|revision| revision == &self.revision)
            || workspace.history.iter().any(|event| {
                event.document_id == self.document_id && event.revision == self.revision
            });
        if !workspace.revisions.contains_key(&self.revision) || !belongs_to_document {
            return Err(bad("associated revision is not retained for this document"));
        }
        Ok(())
    }
    pub fn validate_selected(&self, workspace: &Workspace) -> Result<(), Error> {
        self.validate_retained(workspace)?;
        if !workspace
            .documents
            .iter()
            .any(|d| d.id == self.document_id && d.head == self.revision)
        {
            return Err(bad(
                "selected document revision changed before image attachment",
            ));
        }
        Ok(())
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum Format {
    Jpeg,
    Tiff,
    Png,
}
impl Format {
    pub fn media_type(self) -> &'static str {
        match self {
            Self::Jpeg => "image/jpeg",
            Self::Tiff => "image/tiff",
            Self::Png => "image/png",
        }
    }
    fn extension(self) -> &'static str {
        match self {
            Self::Jpeg => "jpg",
            Self::Tiff => "tif",
            Self::Png => "png",
        }
    }
    fn from_image(format: collagene_image_evidence::Format) -> Self {
        match format {
            collagene_image_evidence::Format::Jpeg => Self::Jpeg,
            collagene_image_evidence::Format::Tiff => Self::Tiff,
            collagene_image_evidence::Format::Png => Self::Png,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum Warning {
    DocumentaryPreview,
    MetadataNotTransferred,
    UncalibratedDisplay,
    AdditionalImagesNotPreviewed,
    PrecisionReducedForPreview,
    DownsampledForPreview,
}
impl Warning {
    fn from_image(warning: collagene_image_evidence::Warning) -> Self {
        match warning {
            collagene_image_evidence::Warning::DocumentaryPreview => Self::DocumentaryPreview,
            collagene_image_evidence::Warning::MetadataNotTransferred => {
                Self::MetadataNotTransferred
            }
            collagene_image_evidence::Warning::UncalibratedDisplay => Self::UncalibratedDisplay,
            collagene_image_evidence::Warning::AdditionalImagesNotPreviewed => {
                Self::AdditionalImagesNotPreviewed
            }
            collagene_image_evidence::Warning::PrecisionReducedForPreview => {
                Self::PrecisionReducedForPreview
            }
            collagene_image_evidence::Warning::DownsampledForPreview => Self::DownsampledForPreview,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Asset {
    pub hash: String,
    pub media_type: String,
    pub bytes: u64,
    pub width: u32,
    pub height: u32,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Record {
    pub version: u16,
    pub engine_version: String,
    pub name: String,
    pub association: Binding,
    pub format: Format,
    pub source_orientation: u8,
    pub original: Asset,
    pub preview: Asset,
    #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, Warning, 6>")]
    pub warnings: Vec<Warning>,
}
impl Record {
    fn from_evidence(name: &str, association: Binding, evidence: &Evidence) -> Self {
        let format = Format::from_image(evidence.format());
        let (width, height) = evidence.source_dimensions();
        let (preview_width, preview_height) = evidence.preview_dimensions();
        Self {
            version: VERSION,
            engine_version: PREVIEW_ENGINE_VERSION.into(),
            name: name.into(),
            association,
            format,
            source_orientation: evidence.source_orientation(),
            original: Asset {
                hash: blake3::hash(evidence.original_bytes()).to_hex().to_string(),
                media_type: format.media_type().into(),
                bytes: evidence.original_bytes().len() as u64,
                width,
                height,
            },
            preview: Asset {
                hash: blake3::hash(evidence.preview_png()).to_hex().to_string(),
                media_type: "image/png".into(),
                bytes: evidence.preview_png().len() as u64,
                width: preview_width,
                height: preview_height,
            },
            warnings: evidence
                .warnings()
                .iter()
                .copied()
                .map(Warning::from_image)
                .collect(),
        }
    }
    fn validate_fields(&self) -> Result<(), Error> {
        if self.version != VERSION || self.engine_version != PREVIEW_ENGINE_VERSION {
            return Err(Error::Unsupported(
                "image evidence record/preview engine version".into(),
            ));
        }
        if !name_valid(&self.name) || !(1..=8).contains(&self.source_orientation) {
            return Err(bad("invalid name/orientation"));
        }
        self.association.validate_fields()?;
        let limits = Limits::default();
        if self.original.media_type != self.format.media_type()
            || self.preview.media_type != "image/png"
            || !digest(&self.original.hash)
            || !digest(&self.preview.hash)
            || self.original.bytes == 0
            || self.original.bytes > limits.input_bytes as u64
            || self.preview.bytes == 0
            || self.preview.bytes > limits.preview_bytes as u64
            || self.original.width == 0
            || self.original.height == 0
            || self.original.width > limits.dimension
            || self.original.height > limits.dimension
            || u64::from(self.original.width) * u64::from(self.original.height)
                > limits.source_pixels
            || self.preview.width == 0
            || self.preview.height == 0
            || self.preview.width > limits.preview_dimension
            || self.preview.height > limits.preview_dimension
            || u64::from(self.preview.width) * u64::from(self.preview.height)
                > limits.preview_pixels
        {
            return Err(bad("invalid source/preview descriptor"));
        }
        if self.warnings.len() < 3
            || self.warnings.len() > 6
            || self.warnings[..3]
                != [
                    Warning::DocumentaryPreview,
                    Warning::MetadataNotTransferred,
                    Warning::UncalibratedDisplay,
                ]
        {
            return Err(bad("documentary warnings are missing"));
        }
        // Canonical warning order also rejects duplicates and forged extra enum values.
        let tail = [
            Warning::AdditionalImagesNotPreviewed,
            Warning::PrecisionReducedForPreview,
            Warning::DownsampledForPreview,
        ];
        let mut last = None;
        for warning in &self.warnings[3..] {
            let index = tail
                .iter()
                .position(|w| w == warning)
                .ok_or_else(|| bad("invalid warning order"))?;
            if last.is_some_and(|previous| previous >= index) {
                return Err(bad("duplicate/unordered warning"));
            }
            last = Some(index);
        }
        Ok(())
    }
    fn validate(&self, workspace: &Workspace) -> Result<(), Error> {
        self.validate_fields()?;
        self.association.validate_retained(workspace)?;
        self.asset(workspace, &self.original)?;
        self.asset(workspace, &self.preview)?;
        Ok(())
    }
    // Private accessor: ImageEvidence construction / workspace validation hashes
    // each immutable image once. Do not rehash a shared 16 MiB source per record.
    fn asset<'w>(
        &self,
        workspace: &'w Workspace,
        descriptor: &Asset,
    ) -> Result<&'w WorkspaceAttachment, Error> {
        let asset = workspace
            .attachments
            .get(&descriptor.hash)
            .ok_or_else(|| bad("image asset is missing"))?;
        if asset.kind != IMAGE_KIND
            || asset.media_type != descriptor.media_type
            || !name_valid(&asset.name)
            || asset.bytes.len() as u64 != descriptor.bytes
        {
            return Err(bad("image asset role/media/name/size/hash differs"));
        }
        Ok(asset)
    }
    fn load(workspace: &Workspace, hash: &str) -> Result<Self, Error> {
        if !digest(hash) {
            return Err(bad("invalid manifest hash"));
        }
        let asset = workspace
            .attachments
            .get(hash)
            .ok_or_else(|| bad("image manifest is missing"))?;
        if asset.kind != MANIFEST_KIND
            || asset.media_type != "application/json"
            || asset.bytes.len() > MAX_RECORD_BYTES
            || !name_valid(&asset.name)
            || blake3::hash(&asset.bytes).to_hex().as_str() != hash
        {
            return Err(bad("manifest role/media/name/size/hash differs"));
        }
        let record = Self::decode(&asset.bytes, hash)?;
        record.validate(workspace)?;
        Ok(record)
    }
    fn decode(bytes: &[u8], hash: &str) -> Result<Self, Error> {
        if bytes.len() > MAX_RECORD_BYTES
            || !digest(hash)
            || blake3::hash(bytes).to_hex().as_str() != hash
        {
            return Err(bad("manifest size/hash differs"));
        }
        let record: Self = serde_json::from_slice(bytes).map_err(bad)?;
        record.validate_fields()?;
        Ok(record)
    }
    fn summary(&self, manifest_hash: String) -> Summary {
        Summary {
            manifest_hash,
            record: self.clone(),
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Summary {
    pub manifest_hash: String,
    pub record: Record,
}
impl Summary {
    pub fn reopen_revision(&self) -> Result<String, Error> {
        self.record.validate_fields()?;
        if !digest(&self.manifest_hash) {
            return Err(bad("invalid manifest hash"));
        }
        identity(&("documentary-image-reopen-v1", self))
    }
}
/// A bounded documentary display result. Archive reopening constructs this
/// only after fresh codec recomputation; it grants no experimental authority.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Preview {
    pub summary: Summary,
    #[serde(
        serialize_with = "collagene_model::decode::serialize_bytes",
        deserialize_with = "collagene_model::decode::bounded_bytes::<_, MAX_PREVIEW_BYTES>"
    )]
    pub png: Vec<u8>,
}
impl Preview {
    /// Hash/metadata admission on the UI thread; full codec work stays private.
    pub fn validate_binding(&self, expected: &Summary) -> Result<(), Error> {
        self.summary.reopen_revision()?;
        if self.summary.manifest_hash != expected.manifest_hash
            || self.summary.record != expected.record
            || self.png.len() > MAX_PREVIEW_BYTES
            || self.png.len() as u64 != expected.record.preview.bytes
            || blake3::hash(&self.png).to_hex().as_str() != expected.record.preview.hash
            || self.png.len() < 33
            || &self.png[..8] != b"\x89PNG\r\n\x1a\n"
            || &self.png[12..16] != b"IHDR"
            || u32::from_be_bytes(self.png[16..20].try_into().unwrap())
                != expected.record.preview.width
            || u32::from_be_bytes(self.png[20..24].try_into().unwrap())
                != expected.record.preview.height
        {
            return Err(bad("preview bytes or selected manifest differ"));
        }
        Ok(())
    }
    pub fn validate_transport_binding(&self, revision: &str) -> Result<(), Error> {
        if self.summary.reopen_revision()? != revision {
            return Err(bad("preview request context differs"));
        }
        self.validate_binding(&self.summary)
    }
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Capture {
    pub source_workspace_fingerprint: String,
    pub summary: Summary,
    #[cfg_attr(feature = "json-schema", schemars(with = "serde_json::Value"))]
    pub workspace: Workspace,
}
/// The decoded owner is freshly constructed locally, never restored from JSON.
#[derive(Debug)]
pub struct Opened {
    pub summary: Summary,
    pub evidence: Evidence,
}
impl Opened {
    pub fn into_preview(self) -> Preview {
        Preview {
            summary: self.summary,
            png: self.evidence.preview_png().to_vec(),
        }
    }
}
#[derive(Debug)]
pub struct Original<'w> {
    pub name: String,
    pub media_type: &'static str,
    pub content_hash: String,
    pub bytes: &'w [u8],
}

/// Bounded worker packet. Preparation and reopening need no Workspace. A received
/// packet remains untrusted; every operation validates hashes and current manifest.
/// It carries no verification receipt or physical-material authority.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Prepared {
    manifest_hash: String,
    #[serde(
        serialize_with = "collagene_model::decode::serialize_bytes",
        deserialize_with = "collagene_model::decode::bounded_bytes::<_, MAX_RECORD_BYTES>"
    )]
    manifest: Vec<u8>,
    #[serde(
        serialize_with = "collagene_model::decode::serialize_bytes",
        deserialize_with = "collagene_model::decode::bounded_bytes::<_, MAX_SOURCE_BYTES>"
    )]
    original: Vec<u8>,
    #[serde(
        serialize_with = "collagene_model::decode::serialize_bytes",
        deserialize_with = "collagene_model::decode::bounded_bytes::<_, MAX_PREVIEW_BYTES>"
    )]
    preview: Vec<u8>,
}
impl Prepared {
    pub fn validate_preparation_binding(&self, revision: &str) -> Result<(), Error> {
        let record = self.record()?;
        if record
            .association
            .preparation_revision(&record.name, &record.original.hash)?
            != revision
        {
            return Err(bad("prepared source/name/association differs"));
        }
        Ok(())
    }
    /// Use after validating the fresh worker response against its captured
    /// source. Existing archive previews require `reopen` in the worker first.
    pub fn preview(&self) -> Result<Preview, Error> {
        let preview = Preview {
            summary: self.summary()?,
            png: self.preview.clone(),
        };
        preview.validate_binding(&preview.summary)?;
        Ok(preview)
    }
    /// CPU-heavy construction, for private workers. Send only original+binding.
    pub fn new(name: &str, association: Binding, original: Vec<u8>) -> Result<Self, Error> {
        if !name_valid(name) {
            return Err(bad("a 1–256 byte image name is required"));
        }
        association.validate_fields()?;
        let evidence = Evidence::prepare(original, Limits::default()).map_err(bad)?;
        let record = Record::from_evidence(name, association, &evidence);
        let manifest = serde_json::to_vec(&record).map_err(bad)?;
        if manifest.len() > MAX_RECORD_BYTES {
            return Err(Error::TooLarge);
        }
        let manifest_hash = blake3::hash(&manifest).to_hex().to_string();
        let preview = evidence.preview_png().to_vec();
        Ok(Self {
            manifest_hash,
            manifest,
            original: evidence.into_original(),
            preview,
        })
    }
    fn record(&self) -> Result<Record, Error> {
        let record = Record::decode(&self.manifest, &self.manifest_hash)?;
        if self.original.len() > MAX_SOURCE_BYTES
            || self.preview.len() > MAX_PREVIEW_BYTES
            || self.original.len() as u64 != record.original.bytes
            || self.preview.len() as u64 != record.preview.bytes
            || blake3::hash(&self.original).to_hex().as_str() != record.original.hash
            || blake3::hash(&self.preview).to_hex().as_str() != record.preview.hash
        {
            return Err(bad("prepared image asset size/hash differs"));
        }
        if !self.preview.starts_with(b"\x89PNG\r\n\x1a\n") {
            return Err(bad("prepared preview is not PNG"));
        }
        Ok(record)
    }
    pub fn summary(&self) -> Result<Summary, Error> {
        Ok(self.record()?.summary(self.manifest_hash.clone()))
    }
    /// Cheap metadata/hash admission plus atomic local workspace replacement.
    /// No codec is invoked. The host must bind the worker response to its request
    /// and current workspace CAS; this method also rejects a changed document head.
    pub fn apply_to(&self, workspace: &Workspace) -> Result<Capture, Error> {
        workspace.validate()?;
        let record = self.record()?;
        record.association.validate_selected(workspace)?;
        let source_workspace_fingerprint =
            crate::material::lineage::workspace_fingerprint(workspace)?;
        let mut additions = BTreeMap::new();
        stage(
            &mut additions,
            record.original.hash.clone(),
            WorkspaceAttachment {
                kind: IMAGE_KIND.into(),
                name: format!(
                    "image-{}.{}",
                    record.original.hash,
                    record.format.extension()
                ),
                media_type: record.original.media_type.clone(),
                bytes: self.original.clone(),
            },
        )?;
        stage(
            &mut additions,
            record.preview.hash.clone(),
            WorkspaceAttachment {
                kind: IMAGE_KIND.into(),
                name: format!("image-{}.png", record.preview.hash),
                media_type: "image/png".into(),
                bytes: self.preview.clone(),
            },
        )?;
        stage(
            &mut additions,
            self.manifest_hash.clone(),
            WorkspaceAttachment {
                kind: MANIFEST_KIND.into(),
                name: format!("image-evidence-{}.json", self.manifest_hash),
                media_type: "application/json".into(),
                bytes: self.manifest.clone(),
            },
        )?;
        let mut new_count = 0;
        for (key, asset) in &additions {
            if let Some(existing) = workspace.attachments.get(key) {
                same_asset(existing, asset)?;
            } else {
                new_count += 1;
            }
        }
        let existing_records = workspace
            .attachments
            .values()
            .filter(|a| a.kind == MANIFEST_KIND)
            .count();
        if workspace.attachments.len() + new_count > MAX_ATTACHMENTS
            || existing_records
                + usize::from(!workspace.attachments.contains_key(&self.manifest_hash))
                > MAX_IMAGES
        {
            return Err(Error::TooLarge);
        }
        let mut next = workspace.clone();
        for (key, asset) in additions {
            next.attachments.entry(key).or_insert(asset);
        }
        next.validate()?;
        record.validate(&next)?;
        Ok(Capture {
            source_workspace_fingerprint,
            summary: record.summary(self.manifest_hash.clone()),
            workspace: next,
        })
    }
    /// Fresh decode of only this bounded packet, for a private worker. Exact
    /// recomputation rejects tampered metadata and previews without unrelated docs.
    pub fn reopen(self) -> Result<Opened, Error> {
        let record = self.record()?;
        let evidence = Evidence::prepare(self.original, Limits::default()).map_err(bad)?;
        let reconstructed =
            Record::from_evidence(&record.name, record.association.clone(), &evidence);
        if reconstructed != record || evidence.preview_png() != self.preview.as_slice() {
            return Err(bad(
                "current decode/preview differs from archived descriptor",
            ));
        }
        Ok(Opened {
            summary: record.summary(self.manifest_hash),
            evidence,
        })
    }
}

/// Workspace attachment owner: metadata/hash admission, no codecs. Historical
/// bindings remain at the explicitly recorded revision when document heads change.
pub struct ImageEvidence<'w> {
    workspace: &'w Workspace,
}
impl<'w> ImageEvidence<'w> {
    pub fn new(workspace: &'w Workspace) -> Result<Self, Error> {
        workspace.validate()?;
        Ok(Self { workspace })
    }
    pub fn list(&self) -> Result<Vec<Summary>, Error> {
        self.workspace
            .attachments
            .iter()
            .filter(|(_, a)| a.kind == MANIFEST_KIND)
            .map(|(hash, _)| Ok(Record::load(self.workspace, hash)?.summary(hash.clone())))
            .collect()
    }
    /// Extract only the bounded image packet to send to the private reopen worker.
    pub fn prepared(&self, manifest_hash: &str) -> Result<Prepared, Error> {
        let record = Record::load(self.workspace, manifest_hash)?;
        Ok(Prepared {
            manifest_hash: manifest_hash.into(),
            manifest: self.workspace.attachments[manifest_hash].bytes.clone(),
            original: record
                .asset(self.workspace, &record.original)?
                .bytes
                .clone(),
            preview: record.asset(self.workspace, &record.preview)?.bytes.clone(),
        })
    }
    /// Exact source export after metadata/hash checks, without promoting validity.
    pub fn original(&self, manifest_hash: &str) -> Result<Original<'w>, Error> {
        let record = Record::load(self.workspace, manifest_hash)?;
        let asset = record.asset(self.workspace, &record.original)?;
        Ok(Original {
            name: format!(
                "image-source-{}.{}",
                record.original.hash,
                record.format.extension()
            ),
            media_type: record.format.media_type(),
            content_hash: record.original.hash,
            bytes: &asset.bytes,
        })
    }
}

/// Hook for archive validation: no codecs, no Workspace::validate recursion, no
/// reconstructed physical/history authority. Keep associated revisions when pruning.
pub fn validate_workspace_records(workspace: &Workspace) -> Result<(), Error> {
    let mut count = 0;
    for (hash, asset) in &workspace.attachments {
        if asset.kind.starts_with("image-evidence-record-") && asset.kind != MANIFEST_KIND
            || asset.kind.starts_with("image-evidence-image-") && asset.kind != IMAGE_KIND
        {
            return Err(Error::Unsupported(
                "internal image evidence attachment version".into(),
            ));
        }
        if asset.kind == MANIFEST_KIND {
            count += 1;
            if count > MAX_IMAGES {
                return Err(Error::TooLarge);
            }
            Record::load(workspace, hash)?;
        } else if asset.kind == IMAGE_KIND {
            let limits = Limits::default();
            if !digest(hash)
                || asset.bytes.is_empty()
                || asset.bytes.len() > limits.input_bytes
                || !matches!(
                    asset.media_type.as_str(),
                    "image/jpeg" | "image/tiff" | "image/png"
                )
                || !name_valid(&asset.name)
                || blake3::hash(&asset.bytes).to_hex().as_str() != hash
            {
                return Err(bad("invalid documentary image asset"));
            }
        }
    }
    Ok(())
}
fn name_valid(value: &str) -> bool {
    !value.trim().is_empty() && value.len() <= 256 && !value.chars().any(char::is_control)
}
fn digest(value: &str) -> bool {
    value.len() == 64
        && value
            .bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}
fn bad(value: impl std::fmt::Display) -> Error {
    Error::Format(format!("image evidence: {value}"))
}
fn identity(value: &impl Serialize) -> Result<String, Error> {
    let bytes = serde_json::to_vec(value).map_err(bad)?;
    if bytes.len() > MAX_RECORD_BYTES {
        return Err(Error::TooLarge);
    }
    Ok(blake3::hash(&bytes).to_hex().to_string())
}
fn same_asset(a: &WorkspaceAttachment, b: &WorkspaceAttachment) -> Result<(), Error> {
    if a.kind != b.kind || a.media_type != b.media_type || a.name != b.name || a.bytes != b.bytes {
        return Err(bad(
            "content hash already has a conflicting attachment role/name/media",
        ));
    }
    Ok(())
}
fn stage(
    additions: &mut BTreeMap<String, WorkspaceAttachment>,
    hash: String,
    asset: WorkspaceAttachment,
) -> Result<(), Error> {
    if let Some(existing) = additions.get(&hash) {
        same_asset(existing, &asset)
    } else {
        additions.insert(hash, asset);
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use collagene_model::{Molecule, Topology};

    fn workspace() -> Workspace {
        Workspace::from_documents(
            "image fixture",
            vec![Molecule::dna("p", "ACGT", Topology::Circular).unwrap()],
            0,
        )
        .unwrap()
    }
    fn binding(w: &Workspace) -> Binding {
        Binding {
            document_id: w.documents[0].id.clone(),
            revision: w.documents[0].head.clone(),
        }
    }
    // Hand-built uncompressed TIFF, one 8-bit grayscale pixel, nine entries.
    fn source() -> Vec<u8> {
        let offset = 8 + 2 + 9 * 12 + 4;
        let mut bytes = b"II\x2a\0\x08\0\0\0".to_vec();
        bytes.extend(9u16.to_le_bytes());
        for (tag, kind, value) in [
            (256u16, 4u16, 1u32),
            (257, 4, 1),
            (258, 3, 8),
            (259, 3, 1),
            (262, 3, 1),
            (273, 4, offset as u32),
            (277, 3, 1),
            (278, 4, 1),
            (279, 4, 1),
        ] {
            bytes.extend(tag.to_le_bytes());
            bytes.extend(kind.to_le_bytes());
            bytes.extend(1u32.to_le_bytes());
            bytes.extend(value.to_le_bytes());
        }
        bytes.extend(0u32.to_le_bytes());
        assert_eq!(bytes.len(), offset);
        bytes.push(127);
        bytes
    }
    fn prepare(w: &Workspace) -> Prepared {
        Prepared::new("documentary image", binding(w), source()).unwrap()
    }

    #[test]
    fn prepare_apply_list_reopen_and_exact_export_never_change_molecules_or_history() {
        let w = workspace();
        let before = serde_json::to_vec(&(
            &w.documents,
            &w.revisions,
            &w.history,
            &w.candidates,
            &w.active_candidate,
        ))
        .unwrap();
        let prepared = prepare(&w);
        let capture = prepared.apply_to(&w).unwrap();
        assert!(w.attachments.is_empty());
        assert_eq!(capture.workspace.attachments.len(), 3);
        let after = serde_json::to_vec(&(
            &capture.workspace.documents,
            &capture.workspace.revisions,
            &capture.workspace.history,
            &capture.workspace.candidates,
            &capture.workspace.active_candidate,
        ))
        .unwrap();
        assert_eq!(before, after);
        let owner = ImageEvidence::new(&capture.workspace).unwrap();
        let list = owner.list().unwrap();
        assert_eq!(list.len(), 1);
        assert_eq!(list[0].record.association, binding(&w));
        let exported = owner.original(&capture.summary.manifest_hash).unwrap();
        assert_eq!(exported.bytes, source());
        assert_eq!(exported.media_type, "image/tiff");
        let opened = owner
            .prepared(&capture.summary.manifest_hash)
            .unwrap()
            .reopen()
            .unwrap();
        assert_eq!(opened.evidence.original_bytes(), source());
    }
    #[test]
    fn duplicate_preparation_is_idempotent_and_original_equal_preview_deduplicates() {
        let w = workspace();
        let first = prepare(&w);
        let capture = first.apply_to(&w).unwrap();
        let second = first.apply_to(&capture.workspace).unwrap();
        assert_eq!(
            capture.workspace.attachments.len(),
            second.workspace.attachments.len()
        );
        let png = first.preview.clone();
        let png_prepared = Prepared::new("same PNG", binding(&w), png).unwrap();
        assert_eq!(png_prepared.original, png_prepared.preview);
        let capture = png_prepared.apply_to(&w).unwrap();
        assert_eq!(capture.workspace.attachments.len(), 2);
    }
    #[test]
    fn apply_rejects_conflicts_and_stale_selection_without_partial_insertion() {
        let mut w = workspace();
        let prepared = prepare(&w);
        let record = prepared.record().unwrap();
        w.attachments.insert(
            record.original.hash.clone(),
            WorkspaceAttachment {
                kind: "other-image-role".into(),
                name: "other".into(),
                media_type: "image/tiff".into(),
                bytes: source(),
            },
        );
        let before = serde_json::to_vec(&w).unwrap();
        assert!(prepared.apply_to(&w).is_err());
        assert_eq!(serde_json::to_vec(&w).unwrap(), before);
        let mut w = workspace();
        let newer = Molecule::dna("changed", "TGCA", Topology::Circular).unwrap();
        let revision = newer.revision();
        w.revisions.insert(revision.clone(), newer);
        w.documents[0].head = revision.clone();
        let event = crate::documents::workspace::RevisionEvent::new(
            &w.documents[0].id,
            &revision,
            vec![w.documents[0].head_event.clone().unwrap()],
            "Edited document",
        )
        .unwrap();
        w.documents[0].head_event = Some(event.id.clone());
        w.history.push(event);
        assert!(prepared.apply_to(&w).is_err());
        assert!(w.attachments.is_empty());
    }
    #[test]
    fn historical_image_binding_survives_changed_document_head() {
        let w = workspace();
        let capture = prepare(&w).apply_to(&w).unwrap();
        let original_binding = capture.summary.record.association.clone();
        let mut w = capture.workspace;
        let newer = Molecule::dna("changed", "TGCA", Topology::Circular).unwrap();
        let revision = newer.revision();
        w.revisions.insert(revision.clone(), newer);
        w.documents[0].head = revision.clone();
        let event = crate::documents::workspace::RevisionEvent::new(
            &w.documents[0].id,
            &revision,
            vec![w.documents[0].head_event.clone().unwrap()],
            "Edited document",
        )
        .unwrap();
        w.documents[0].head_event = Some(event.id.clone());
        w.history.push(event);
        let owner = ImageEvidence::new(&w).unwrap();
        assert_eq!(
            owner.list().unwrap()[0].record.association,
            original_binding
        );
        assert!(
            owner
                .prepared(&capture.summary.manifest_hash)
                .unwrap()
                .reopen()
                .is_ok()
        );
    }
    #[test]
    fn manifest_versions_hashes_and_attachment_roles_are_strict() {
        let w = workspace();
        let prepared = prepare(&w);
        let mut unsupported = prepared.clone();
        let mut record = unsupported.record().unwrap();
        record.version += 1;
        unsupported.manifest = serde_json::to_vec(&record).unwrap();
        unsupported.manifest_hash = blake3::hash(&unsupported.manifest).to_hex().to_string();
        assert!(matches!(
            unsupported.apply_to(&w),
            Err(Error::Unsupported(_))
        ));
        let mut tampered = prepared.clone();
        tampered.original[0] ^= 1;
        assert!(tampered.apply_to(&w).is_err());
        let capture = prepared.apply_to(&w).unwrap();
        let mut invalid = capture.workspace.clone();
        invalid
            .attachments
            .get_mut(&capture.summary.record.original.hash)
            .unwrap()
            .media_type = "image/png".into();
        assert!(ImageEvidence::new(&invalid).is_err());
        let mut invalid = capture.workspace;
        invalid
            .attachments
            .get_mut(&capture.summary.manifest_hash)
            .unwrap()
            .kind = "image-evidence-record-v0".into();
        assert!(validate_workspace_records(&invalid).is_err());
    }
    #[test]
    fn reopen_recomputes_and_rejects_plausible_but_forged_preview_metadata() {
        let w = workspace();
        let mut prepared = prepare(&w);
        let mut record = prepared.record().unwrap();
        record.original.width = 2;
        prepared.manifest = serde_json::to_vec(&record).unwrap();
        prepared.manifest_hash = blake3::hash(&prepared.manifest).to_hex().to_string();
        // Metadata admission is deliberately not a restored decoder receipt.
        assert!(prepared.summary().is_ok());
        assert!(prepared.reopen().is_err());
    }
    #[test]
    fn forged_association_cannot_bind_another_documents_retained_revision() {
        let w = Workspace::from_documents(
            "two documents",
            vec![
                Molecule::dna("first", "ACGT", Topology::Circular).unwrap(),
                Molecule::dna("second", "TGCA", Topology::Circular).unwrap(),
            ],
            0,
        )
        .unwrap();
        let capture = prepare(&w).apply_to(&w).unwrap();
        let mut forged = capture.workspace;
        let mut record = capture.summary.record;
        record.association.revision = forged.documents[1].head.clone();
        // Both identities independently exist, but this revision is retained only
        // by document two. Rehashing the forged manifest must not admit the link.
        assert!(forged.revisions.contains_key(&record.association.revision));
        assert!(record.association.validate_retained(&forged).is_err());
        let bytes = serde_json::to_vec(&record).unwrap();
        let hash = blake3::hash(&bytes).to_hex().to_string();
        forged.attachments.remove(&capture.summary.manifest_hash);
        forged.attachments.insert(
            hash.clone(),
            WorkspaceAttachment {
                kind: MANIFEST_KIND.into(),
                name: format!("image-evidence-{hash}.json"),
                media_type: "application/json".into(),
                bytes,
            },
        );
        assert!(validate_workspace_records(&forged).is_err());
        assert!(ImageEvidence::new(&forged).is_err());
    }
    #[test]
    fn current_worker_transport_pins_source_and_reopens_exact_archived_pixels() {
        use crate::protocol::{
            execution,
            request::Request,
            response::Response,
            schema::{AnalysisJob, AnalysisResult, WORKER_SCHEMA_VERSION},
            session::AnalysisSession,
        };
        let w = workspace();
        let association = binding(&w);
        let original = source();
        let revision = association
            .preparation_revision("image", &blake3::hash(&original).to_hex())
            .unwrap();
        let request = Request {
            id: 7,
            revision: revision.clone(),
            job: AnalysisJob::PrepareImageEvidence {
                name: "image".into(),
                association,
                original,
            },
        };
        let wire = request.encode().unwrap();
        let ((id, context, tag), _) =
            postcard::take_from_bytes::<(u32, &str, u32)>(&wire[42..]).unwrap();
        assert_eq!((id, context, tag), (7, revision.as_str(), 200));
        let reply = execution::process_request(&wire).unwrap();
        assert_eq!(
            reply,
            AnalysisSession::default().process_request(&wire).unwrap()
        );
        let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&reply[42..]).unwrap();
        assert_eq!(tag, 190);
        let response = Response::decode(&reply).unwrap();
        assert_eq!(response.id, 7);
        let AnalysisResult::ImageEvidencePrepared(prepared) = response.result else {
            panic!("wrong result")
        };
        prepared.validate_preparation_binding(&revision).unwrap();
        let expected = prepared.summary().unwrap();
        let fresh = prepared.preview().unwrap();
        fresh.validate_binding(&expected).unwrap();
        let reopened_request = Request {
            id: 8,
            revision: expected.reopen_revision().unwrap(),
            job: AnalysisJob::ReopenImageEvidence { prepared },
        };
        let wire = reopened_request.encode().unwrap();
        let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&wire[42..]).unwrap();
        assert_eq!(tag, 201);
        let reply = execution::process_request(&wire).unwrap();
        let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&reply[42..]).unwrap();
        assert_eq!(tag, 191);
        let response = Response::decode(&reply).unwrap();
        let AnalysisResult::ImageEvidencePreview(preview) = response.result else {
            panic!("wrong result")
        };
        preview.validate_binding(&expected).unwrap();
        assert_eq!(preview.png, fresh.png);
        let mut old = wire;
        old[4..6].copy_from_slice(&(WORKER_SCHEMA_VERSION - 1).to_le_bytes());
        assert!(execution::process_request(&old).is_err());
    }
    #[test]
    fn transport_refuses_rebound_source_name_manifest_and_preview_bytes() {
        use crate::protocol::{
            request::Request,
            response::Response,
            schema::{AnalysisJob, AnalysisResult},
        };
        let w = workspace();
        let prepared = prepare(&w);
        let summary = prepared.summary().unwrap();
        let context = summary
            .record
            .association
            .preparation_revision(&summary.record.name, &summary.record.original.hash)
            .unwrap();
        let mut wrong = context.clone();
        wrong.replace_range(..1, if context.starts_with('a') { "b" } else { "a" });
        assert!(prepared.validate_preparation_binding(&wrong).is_err());
        let mut request = Request {
            id: 1,
            revision: context,
            job: AnalysisJob::PrepareImageEvidence {
                name: summary.record.name.clone(),
                association: summary.record.association.clone(),
                original: source(),
            },
        };
        request.encode().unwrap();
        if let AnalysisJob::PrepareImageEvidence { original, .. } = &mut request.job {
            *original.last_mut().unwrap() ^= 1;
        }
        assert!(request.encode().is_err());
        let mut preview = prepared.preview().unwrap();
        preview.png[0] ^= 1;
        assert!(preview.validate_binding(&summary).is_err());
        let response = Response {
            id: 1,
            revision: summary.reopen_revision().unwrap(),
            result: AnalysisResult::ImageEvidencePreview(preview),
        };
        assert!(Response::decode(&response.encode().unwrap()).is_err());
        let mut wrong_manifest = summary.clone();
        wrong_manifest.manifest_hash = "a".repeat(64);
        assert!(
            prepared
                .preview()
                .unwrap()
                .validate_binding(&wrong_manifest)
                .is_err()
        );
    }
    #[test]
    fn native_workspace_roundtrip_retains_original_preview_and_historical_association() {
        use crate::documents::workspace::{load_workspace, save_workspace};
        let w = workspace();
        let capture = prepare(&w).apply_to(&w).unwrap();
        let bytes = save_workspace(&capture.workspace).unwrap();
        let restored = load_workspace(&bytes).unwrap();
        assert_eq!(save_workspace(&restored).unwrap(), bytes);
        let owner = ImageEvidence::new(&restored).unwrap();
        let summary = &owner.list().unwrap()[0];
        assert_eq!(summary.record.association, binding(&w));
        assert_eq!(
            owner.original(&summary.manifest_hash).unwrap().bytes,
            source()
        );
        let preview = owner
            .prepared(&summary.manifest_hash)
            .unwrap()
            .reopen()
            .unwrap()
            .into_preview();
        preview.validate_binding(summary).unwrap();
        assert_eq!(
            restored.head(&w.documents[0].id).unwrap(),
            w.head(&w.documents[0].id).unwrap()
        );
    }
}
