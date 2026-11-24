//! Raw SAM/BAM source attachments, without restored parser/reference authority.
//! Kind/media describe the caller's declared format; neither proves file validity.
use crate::{Error, read_alignment_store::Format, workspace::WorkspaceAttachment};

pub const SAM_SOURCE_KIND: &str = "sam-source-v1";
pub const BAM_SOURCE_KIND: &str = "bam-source-v1";
pub const SAM_MEDIA_TYPE: &str = "text/x-sam";
pub const BAM_MEDIA_TYPE: &str = "application/x-bam";
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
        Format::Sam => (SAM_SOURCE_KIND, SAM_MEDIA_TYPE),
        Format::Bam => (BAM_SOURCE_KIND, BAM_MEDIA_TYPE),
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
        return Err(Error::Format("invalid read alignment source name".into()));
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
/// accepts bytes without file validation, including an empty file; declared SAM
/// can therefore still fail SAM parsing. A workspace owner must separately match
/// the returned hash to its attachment key. No parsed/verified state is restored.
pub fn inspect(attachment: &WorkspaceAttachment) -> Result<SourceInfo, Error> {
    if attachment.bytes.len() > MAX_SOURCE_BYTES {
        return Err(Error::TooLarge);
    }
    if !name_valid(&attachment.name) {
        return Err(Error::Format("invalid read alignment source name".into()));
    }
    let format = match attachment.kind.as_str() {
        SAM_SOURCE_KIND if attachment.media_type == SAM_MEDIA_TYPE => Format::Sam,
        BAM_SOURCE_KIND if attachment.media_type == BAM_MEDIA_TYPE => Format::Bam,
        _ => {
            return Err(Error::Format(
                "unsupported read alignment source kind/media type".into(),
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn literal_sources_retain_lexical_bytes_names_and_content_keys() {
        let sam = b"@SQ\tSN:ref\tLN:4\r\nr\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*";
        for (bytes, format, name) in [
            (sam.as_slice(), Format::Sam, "literal.sam"),
            (b"\x1f\x8b\0\xff".as_slice(), Format::Bam, "literal.bam"),
        ] {
            let (hash, attachment) = make_attachment(bytes, format, name).unwrap();
            assert_eq!(attachment.bytes, bytes);
            assert_eq!(hash, blake3::hash(bytes).to_hex().to_string());
            let info = inspect(&attachment).unwrap();
            assert_eq!(info.source_blake3, hash);
            assert_eq!(info.source_bytes, bytes.len() as u64);
            assert_eq!(info.name, name);
            assert_eq!(info.format, format);
        }
    }

    #[test]
    fn metadata_inspection_does_not_masquerade_as_file_validation() {
        let (_, attachment) = make_attachment(b"not valid SAM", Format::Sam, "source.bam").unwrap();
        assert_eq!(
            inspect(&attachment).unwrap().format,
            Format::Sam,
            "the display name never infers format"
        );
        assert!(
            crate::read_alignment_store::Collection::import(
                &attachment.bytes,
                Format::Sam,
                &Default::default()
            )
            .is_err()
        );
        let (_, empty) = make_attachment(&[], Format::Bam, "empty.bam").unwrap();
        assert_eq!(inspect(&empty).unwrap().source_bytes, 0);
        assert!(
            crate::read_alignment_store::Collection::import(&[], Format::Bam, &Default::default())
                .is_err()
        );
    }

    #[test]
    fn strict_metadata_and_size_bounds_fail_before_cloning_source() {
        for name in ["", " \t ", "bad\nname", "bad\0name"] {
            assert!(make_attachment(b"source", Format::Sam, name).is_err());
        }
        assert!(make_attachment(b"source", Format::Sam, &"n".repeat(MAX_NAME_BYTES + 1)).is_err());
        let (_, mut attachment) = make_attachment(b"source", Format::Sam, "reads.sam").unwrap();
        attachment.media_type = BAM_MEDIA_TYPE.into();
        assert!(inspect(&attachment).is_err());
        attachment.media_type = SAM_MEDIA_TYPE.into();
        attachment.kind = "sam-source-v999".into();
        assert!(inspect(&attachment).is_err());
        let at_cap = vec![0u8; MAX_SOURCE_BYTES];
        assert!(make_attachment(&at_cap, Format::Bam, "at-cap.bam").is_ok());
        let over_cap = vec![0u8; MAX_SOURCE_BYTES + 1];
        assert!(matches!(
            make_attachment(&over_cap, Format::Bam, "over.bam"),
            Err(Error::TooLarge)
        ));
    }
}
