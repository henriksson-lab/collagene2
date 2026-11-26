//! Regenerated SCF 3.00 from actual retained trace evidence. No basecalling,
//! synthetic signal, normalization, clipping or quality-scale conversion.
use crate::{
    Error,
    traces::{self, QualityEvidence, Trace, TraceFormat, TraceOrientation},
};
use serde::{Deserialize, Serialize};
pub const ENGINE_VERSION: &str = "scf-export-1";
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum SamplePrecision {
    Smallest,
    EightBit,
    SixteenBit,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum LossPolicy {
    Reject,
    AllowDocumentaryOmissions,
}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Options {
    pub precision: SamplePrecision,
    pub loss_policy: LossPolicy,
}
impl Default for Options {
    fn default() -> Self {
        Self {
            precision: SamplePrecision::Smallest,
            loss_policy: LossPolicy::Reject,
        }
    }
}
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Loss {
    pub field: String,
    pub bytes: usize,
    pub detail: String,
}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Exported {
    pub context_revision: String,
    pub bytes: Vec<u8>,
    pub source_hash: String,
    pub orientation: TraceOrientation,
    pub sample_size: u32,
    pub losses: Vec<Loss>,
    pub warnings: Vec<String>,
}
fn preflight(trace: &Trace) -> Result<(), Error> {
    trace.validate()?;
    if trace.metadata.len() > 256 {
        return Err(Error::TooLarge);
    }
    let mut text = 0usize;
    for (key, value) in &trace.metadata {
        if key.len() > 4096 || value.len() > 1024 * 1024 {
            return Err(Error::TooLarge);
        }
        text = text
            .checked_add(key.len())
            .and_then(|n| n.checked_add(value.len()))
            .filter(|n| *n <= 8 * 1024 * 1024)
            .ok_or(Error::TooLarge)?;
    }
    Ok(())
}
/// All supported mutable output parameters are pinned to the actual attachment.
/// Signal/call edits are not parameters: export independently replays the source
/// and refuses any such edits even when its declared source hash is unchanged.
pub fn revision(trace: &Trace, options: &Options) -> Result<String, Error> {
    preflight(trace)?;
    let bytes = serde_json::to_vec(&(
        ENGINE_VERSION,
        &trace.source_hash,
        trace.orientation,
        &trace.name,
        &trace.metadata,
        options,
    ))
    .map_err(|e| Error::Format(e.to_string()))?;
    Ok(blake3::hash(&bytes).to_hex().to_string())
}
impl Exported {
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        let hash = |s: &str| {
            s.len() == 64
                && s.bytes()
                    .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
        };
        if !hash(expected)
            || self.context_revision != expected
            || !hash(&self.source_hash)
            || self.bytes.len() < 128
            || self.bytes.len() > crate::formats::MAX_IMPORT_BYTES
            || self.bytes.get(..4) != Some(b".scf")
            || self.bytes.get(36..40) != Some(b"3.00")
            || !matches!(self.sample_size, 1 | 2)
            || word(&self.bytes, 40)? != self.sample_size as usize
            || self.losses.len() > 32
            || self.warnings.len() > 32
            || self.losses.iter().any(|l| {
                l.field.len() > 4096
                    || l.detail.len() > 8192
                    || l.bytes > crate::formats::MAX_IMPORT_BYTES
            })
            || self.warnings.iter().any(|w| w.len() > 8192)
        {
            return Err(Error::Format(
                "SCF export response/source context changed".into(),
            ));
        }
        Ok(())
    }
}
fn bad(s: &str) -> Error {
    Error::Unsupported(format!("SCF export: {s}"))
}
fn word(bytes: &[u8], at: usize) -> Result<usize, Error> {
    let b = bytes
        .get(at..at.checked_add(4).ok_or(Error::TooLarge)?)
        .ok_or_else(|| bad("truncated source header"))?;
    Ok(u32::from_be_bytes(b.try_into().unwrap()) as usize)
}
fn range(bytes: &[u8], at: usize, len: usize) -> Result<&[u8], Error> {
    if len == 0 {
        return Ok(&bytes[..0]);
    }
    bytes
        .get(at..at.checked_add(len).ok_or(Error::TooLarge)?)
        .ok_or_else(|| bad("source section outside original attachment"))
}
fn set_word(bytes: &mut [u8], at: usize, n: usize) -> Result<(), Error> {
    bytes[at..at + 4]
        .copy_from_slice(&u32::try_from(n).map_err(|_| Error::TooLarge)?.to_be_bytes());
    Ok(())
}
fn loss(losses: &mut Vec<Loss>, field: &str, bytes: usize, detail: &str) {
    if bytes > 0 {
        losses.push(Loss {
            field: field.into(),
            bytes,
            detail: detail.into(),
        });
    }
}
fn documentary_losses(trace: &Trace, comments: &[u8]) -> Result<Vec<Loss>, Error> {
    // The source was fully parsed before this bounded container inspection.
    let raw = &trace.original_bytes;
    let samples = word(raw, 4)?;
    let bases = word(raw, 12)?;
    let sample_size = if &raw[36..40] == b"2.00" {
        1
    } else {
        word(raw, 40)?
    };
    let sections = [
        (
            word(raw, 8)?,
            samples
                .checked_mul(4)
                .and_then(|n| n.checked_mul(sample_size))
                .ok_or(Error::TooLarge)?,
        ),
        (
            word(raw, 24)?,
            bases.checked_mul(12).ok_or(Error::TooLarge)?,
        ),
        (word(raw, 32)?, word(raw, 28)?),
        (word(raw, 52)?, word(raw, 48)?),
    ];
    let mut losses = vec![];
    if range(raw, sections[2].0, sections[2].1)? != comments {
        loss(
            &mut losses,
            "comments",
            sections[2].1,
            "Original comment bytes are omitted; only the current sample NAME is regenerated. The original attachment remains unchanged.",
        );
    }
    loss(
        &mut losses,
        "private_data",
        sections[3].1,
        "Opaque private data has no retained portable interpretation and is omitted.",
    );
    let clips = usize::from(word(raw, 16)? != 0) + usize::from(word(raw, 20)? != 0);
    loss(
        &mut losses,
        "obsolete_clips",
        clips * 4,
        "Obsolete clipping declarations are not applied or transformed; regenerated fields are zero.",
    );
    let header_reserved = raw[56..128].iter().filter(|&&b| b != 0).count();
    loss(
        &mut losses,
        "header_reserved",
        header_reserved,
        "Nonzero unknown header-reserved bytes are omitted.",
    );
    let base = range(raw, sections[1].0, sections[1].1)?;
    let reserved = if &raw[36..40] == b"3.00" {
        base[bases * 9..].iter().filter(|&&b| b != 0).count()
    } else {
        base.chunks_exact(12)
            .map(|r| r[9..12].iter().filter(|&&b| b != 0).count())
            .sum()
    };
    loss(
        &mut losses,
        "base_reserved",
        reserved,
        "Nonzero unknown base-reserved bytes are omitted; regenerated reserved arrays are zero.",
    );
    let mut spans = vec![(0, 128)];
    for (offset, len) in sections {
        if len > 0 {
            spans.push((offset, offset.checked_add(len).ok_or(Error::TooLarge)?));
        }
    }
    spans.sort_unstable();
    let mut stop = 0;
    let mut unclassified = 0;
    for (start, end) in spans {
        unclassified += range(raw, stop, start.saturating_sub(stop))?.len();
        stop = end;
    }
    unclassified += raw
        .len()
        .checked_sub(stop)
        .ok_or_else(|| bad("invalid source section extent"))?;
    loss(
        &mut losses,
        "unclassified_container_bytes",
        unclassified,
        "Padding/trailing bytes outside declared source sections are omitted, without guessing their semantics.",
    );
    let expected_metadata = [
        (
            "scf_reported_left_clip_obsolete",
            word(raw, 16)?.to_string(),
        ),
        (
            "scf_reported_right_clip_obsolete",
            word(raw, 20)?.to_string(),
        ),
        ("scf_private_bytes", sections[3].1.to_string()),
        ("scf_comment_bytes", sections[2].1.to_string()),
    ];
    let extra = trace
        .metadata
        .iter()
        .filter(|(key, value)| {
            !expected_metadata
                .iter()
                .any(|(k, v)| key.as_str() == *k && *value == v)
        })
        .map(|(key, value)| key.len() + value.len())
        .sum();
    loss(
        &mut losses,
        "additional_metadata",
        extra,
        "Additional native trace metadata is not serialized into SCF comments.",
    );
    Ok(losses)
}
/// Export the displayed, explicitly oriented actual evidence. Scalar ABIF scores
/// cannot be converted to SCF's four accuracy vectors without inventing evidence.
pub fn export(trace: &Trace, options: &Options) -> Result<Exported, Error> {
    let context_revision = revision(trace, options)?;
    let QualityEvidence::ScfAccuracy(accuracy)=trace.qualities.as_ref().ok_or_else(||bad("four retained SCF accuracy vectors are required; missing qualities are not replaced with zeroes"))? else{return Err(bad("ABIF scalar basecaller scores cannot be mapped to four SCF accuracy estimates without calibration assumptions"));};
    let TraceFormat::Scf { code_set, .. } = trace.format else {
        return Err(bad(
            "this profile requires retained SCF signal and four-channel accuracy evidence",
        ));
    };
    if trace.name.chars().any(char::is_control) {
        return Err(bad(
            "sample name cannot contain comment delimiters/control characters",
        ));
    }
    // A source hash alone does not bind the editable Trace fields to the signal.
    // Re-read the actual attachment and apply only the declared view transform.
    let original = traces::import_trace(&trace.original_bytes)?;
    let actual = match trace.orientation {
        TraceOrientation::AsStored => original,
        TraceOrientation::ReverseComplement => original.reverse_complement()?,
    };
    if trace.channels != actual.channels
        || trace.sequence != actual.sequence
        || trace.called_symbols != actual.called_symbols
        || trace.peak_positions != actual.peak_positions
        || trace.qualities != actual.qualities
        || trace.format != actual.format
    {
        return Err(bad(
            "trace fields differ from the actual retained attachment/orientation; edited or fabricated samples/calls/qualities are unsupported",
        ));
    }
    let channels = [
        &trace.channels.a,
        &trace.channels.c,
        &trace.channels.g,
        &trace.channels.t,
    ];
    let mut maximum = 0;
    for value in channels.iter().flat_map(|c| c.iter()) {
        if !(0..=65535).contains(value) {
            return Err(bad(
                "sample signal cannot be represented as unsigned 8/16-bit SCF without rescaling or clipping",
            ));
        }
        maximum = maximum.max(*value);
    }
    let precision = match options.precision {
        SamplePrecision::Smallest => {
            if maximum <= 255 {
                1
            } else {
                2
            }
        }
        SamplePrecision::EightBit if maximum > 255 => {
            return Err(bad("8-bit sample precision would lose retained signal"));
        }
        SamplePrecision::EightBit => 1,
        SamplePrecision::SixteenBit => 2,
    };
    let comments = format!("NAME={}\n\0", trace.name).into_bytes();
    let losses = documentary_losses(trace, &comments)?;
    if options.loss_policy == LossPolicy::Reject && !losses.is_empty() {
        return Err(bad(
            "regeneration would omit documentary container fields; review losses with AllowDocumentaryOmissions or export the unchanged original attachment",
        ));
    }
    let n = trace.channels.len();
    let b = trace.sequence.len();
    let sample_bytes = n
        .checked_mul(4)
        .and_then(|n| n.checked_mul(precision))
        .ok_or(Error::TooLarge)?;
    let base_offset = 128usize.checked_add(sample_bytes).ok_or(Error::TooLarge)?;
    let comment_offset = base_offset
        .checked_add(b.checked_mul(12).ok_or(Error::TooLarge)?)
        .ok_or(Error::TooLarge)?;
    let size = comment_offset
        .checked_add(comments.len())
        .filter(|n| *n <= crate::formats::MAX_IMPORT_BYTES)
        .ok_or(Error::TooLarge)?;
    let mut bytes = Vec::new();
    bytes.try_reserve_exact(size).map_err(|_| Error::TooLarge)?;
    bytes.resize(128, 0);
    bytes[..4].copy_from_slice(b".scf");
    set_word(&mut bytes, 4, n)?;
    set_word(&mut bytes, 8, 128)?;
    set_word(&mut bytes, 12, b)?;
    set_word(&mut bytes, 24, base_offset)?;
    set_word(&mut bytes, 28, comments.len())?;
    set_word(&mut bytes, 32, comment_offset)?;
    bytes[36..40].copy_from_slice(b"3.00");
    set_word(&mut bytes, 40, precision)?;
    set_word(&mut bytes, 44, code_set as usize)?;
    let mask = if precision == 1 { 255u32 } else { 65535u32 };
    for channel in channels {
        let mut previous_sample = 0u32;
        let mut previous_difference = 0u32;
        for &value in channel {
            let value = value as u32;
            let difference = value.wrapping_sub(previous_sample) & mask;
            let second = difference.wrapping_sub(previous_difference) & mask;
            if precision == 1 {
                bytes.push(second as u8);
            } else {
                bytes.extend_from_slice(&(second as u16).to_be_bytes());
            }
            previous_sample = value;
            previous_difference = difference;
        }
    }
    for &peak in &trace.peak_positions {
        bytes.extend_from_slice(&peak.to_be_bytes());
    }
    for channel in [&accuracy.a, &accuracy.c, &accuracy.g, &accuracy.t] {
        bytes.extend_from_slice(channel);
    }
    bytes.extend_from_slice(&trace.called_symbols);
    bytes.resize(comment_offset, 0);
    bytes.extend_from_slice(&comments);
    let report = Exported{context_revision:context_revision.clone(),bytes,source_hash:trace.source_hash.clone(),orientation:trace.orientation,sample_size:precision as u32,losses,warnings:vec!["SCF 3.00 regenerated from actual retained signal, calls, peaks and four accuracy vectors. No basecalling, quality recalibration, normalization, implicit clipping or biological verification is asserted. This is not a byte-identical original export.".into()]};
    report.validate_binding(&context_revision)?;
    Ok(report)
}
