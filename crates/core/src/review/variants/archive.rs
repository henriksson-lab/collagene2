//! Raw VCF/BCF source attachments, without restored parser/reference authority.
//! Kind/media describe the caller's declared format; neither proves file validity.
use crate::{
    documents::workspace::WorkspaceAttachment, review::variants::session::FileFormat as Format,
};
use collagene_model::Error;

pub const VCF_SOURCE_KIND: &str = "vcf-source-v1";
pub const BCF_SOURCE_KIND: &str = "bcf-source-v1";
pub const VCF_MEDIA_TYPE: &str = "application/octet-stream";
pub const BCF_MEDIA_TYPE: &str = "application/octet-stream";
pub const MAX_SOURCE_BYTES: usize = 16 * 1024 * 1024;
pub const MAX_NAME_BYTES: usize = 1024;

/// Source metadata only. No Deserialize, parser-completion receipt, reference
/// check or biological validity flag can be restored from this summary.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SourceInfo {
    pub name: String,
    pub format: Format,
    pub source_blake3: String,
    pub source_bytes: u64,
}

fn name_valid(name: &str) -> bool {
    name.len() <= MAX_NAME_BYTES && !name.trim().is_empty() && !name.chars().any(char::is_control)
}

fn metadata(format: Format) -> (&'static str, &'static str) {
    match format {
        Format::Vcf => (VCF_SOURCE_KIND, VCF_MEDIA_TYPE),
        Format::Bcf => (BCF_SOURCE_KIND, BCF_MEDIA_TYPE),
    }
}

/// Retain exact bytes, including malformed/empty sources. This does not parse,
/// normalize line endings, reframe BGZF, infer format from a name, or fetch a URI.
/// On reopening, run Collection::import against these bytes and explicit format.
/// All source/name bounds are checked before the source vector is copied.
pub fn make_attachment(
    bytes: &[u8],
    format: Format,
    name: &str,
) -> Result<(String, WorkspaceAttachment), Error> {
    if bytes.len() > MAX_SOURCE_BYTES {
        return Err(Error::TooLarge);
    }
    if !name_valid(name) {
        return Err(Error::Format("invalid variant source name".into()));
    }
    let (kind, media_type) = metadata(format);
    let source_blake3 = blake3::hash(bytes).to_hex().to_string();
    Ok((
        source_blake3,
        WorkspaceAttachment {
            kind: kind.into(),
            name: name.into(),
            media_type: media_type.into(),
            bytes: bytes.to_vec(),
        },
    ))
}

/// Inspect bounded source metadata and compute its actual content hash. This
/// accepts bytes without file validation, including an empty file; declared VCF
/// can therefore still fail VCF parsing. A workspace owner must separately match
/// the returned hash to its attachment key. No parsed/verified state is restored.
pub fn inspect(attachment: &WorkspaceAttachment) -> Result<SourceInfo, Error> {
    if attachment.bytes.len() > MAX_SOURCE_BYTES {
        return Err(Error::TooLarge);
    }
    if !name_valid(&attachment.name) {
        return Err(Error::Format("invalid variant source name".into()));
    }
    let format = match attachment.kind.as_str() {
        VCF_SOURCE_KIND if attachment.media_type == VCF_MEDIA_TYPE => Format::Vcf,
        BCF_SOURCE_KIND if attachment.media_type == BCF_MEDIA_TYPE => Format::Bcf,
        _ => {
            return Err(Error::Format(
                "unsupported variant source kind/media type".into(),
            ));
        }
    };
    Ok(SourceInfo {
        name: attachment.name.clone(),
        format,
        source_blake3: blake3::hash(&attachment.bytes).to_hex().to_string(),
        source_bytes: attachment.bytes.len() as u64,
    })
}
