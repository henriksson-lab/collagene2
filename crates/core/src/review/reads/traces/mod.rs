//! Independently implemented ABIF/SCF readers from the vendor and Staden specifications.
//! No basecalling, quality recalibration, implicit clipping or external programs.
pub mod display;
use crate::formats::sequence::MAX_IMPORT_BYTES;
use base64::{Engine as _, engine::general_purpose::STANDARD};
use collagene_model::{Alphabet, Error, Molecule, Sequence, Topology};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

pub const MAX_TRACE_SAMPLES: usize = 2_000_000;
pub const MAX_TRACE_BASES: usize = 100_000;
const MAX_ABIF_ENTRIES: usize = 4096;

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub enum TraceFormat {
    Abif { version: u16 },
    Scf { version: String, code_set: u32 },
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum TraceOrientation {
    AsStored,
    ReverseComplement,
}
/// Canonical channel labels. Intensities are not rescaled: signed ABIF shorts stay signed.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct TraceChannels {
    pub a: Vec<i32>,
    pub c: Vec<i32>,
    pub g: Vec<i32>,
    pub t: Vec<i32>,
}
impl TraceChannels {
    pub fn len(&self) -> usize {
        self.a.len()
    }
    pub fn is_empty(&self) -> bool {
        self.a.is_empty()
    }
    /// Borrow channels in nucleotide order: A, C, G, T.
    pub fn all(&self) -> [&[i32]; 4] {
        [&self.a, &self.c, &self.g, &self.t]
    }
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct TraceAccuracy {
    pub a: Vec<u8>,
    pub c: Vec<u8>,
    pub g: Vec<u8>,
    pub t: Vec<u8>,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub enum QualityEvidence {
    /// PCON bytes, in the original caller's scale. They are not asserted to be Phred.
    AbifBasecallerScores(Vec<u8>),
    /// Four SCF accuracy estimates per call, retained without calibration/conversion.
    ScfAccuracy(TraceAccuracy),
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Trace {
    pub name: String,
    pub format: TraceFormat,
    pub orientation: TraceOrientation,
    pub sequence: Sequence,
    /// Exact stored call symbols, transformed alongside explicit reverse-complement views.
    /// SCF's documented '-' no-call symbol is projected to N in sequence, never deleted.
    pub called_symbols: Vec<u8>,
    /// Zero-based sample indexes, in the displayed orientation. Equal positions are retained.
    pub peak_positions: Vec<u32>,
    pub channels: TraceChannels,
    pub qualities: Option<QualityEvidence>,
    pub metadata: BTreeMap<String, String>,
    pub warnings: Vec<String>,
    /// Unchanged original attachment. It is never executed or treated as instructions.
    pub original_bytes: Vec<u8>,
    pub source_hash: String,
}
fn format_error(message: impl Into<String>) -> Error {
    Error::Format(message.into())
}
fn unsupported(message: impl Into<String>) -> Error {
    Error::Unsupported(message.into())
}
fn span(bytes: &[u8], offset: usize, length: usize) -> Result<&[u8], Error> {
    let end = offset.checked_add(length).ok_or(Error::TooLarge)?;
    bytes
        .get(offset..end)
        .ok_or_else(|| format_error("trace section exceeds file bounds"))
}
fn u16be(bytes: &[u8], offset: usize) -> Result<u16, Error> {
    Ok(u16::from_be_bytes(
        span(bytes, offset, 2)?.try_into().unwrap(),
    ))
}
fn u32be(bytes: &[u8], offset: usize) -> Result<u32, Error> {
    Ok(u32::from_be_bytes(
        span(bytes, offset, 4)?.try_into().unwrap(),
    ))
}
fn nonnegative(bytes: &[u8], offset: usize) -> Result<usize, Error> {
    let value = i32::from_be_bytes(span(bytes, offset, 4)?.try_into().unwrap());
    usize::try_from(value)
        .map_err(|_| format_error("negative ABIF count, size or directory offset"))
}
fn product(a: usize, b: usize) -> Result<usize, Error> {
    a.checked_mul(b).ok_or(Error::TooLarge)
}
fn source_limit(bytes: &[u8]) -> Result<(), Error> {
    if bytes.len() > MAX_IMPORT_BYTES {
        Err(Error::TooLarge)
    } else {
        Ok(())
    }
}
pub(super) fn calls(raw: &[u8], dash_is_no_call: bool) -> Result<Sequence, Error> {
    if raw.is_empty() {
        return Err(format_error("trace has no base calls"));
    }
    if raw.len() > MAX_TRACE_BASES {
        return Err(Error::TooLarge);
    }
    let mut normalized = String::with_capacity(raw.len());
    for &b in raw {
        let upper = b.to_ascii_uppercase();
        if dash_is_no_call && b == b'-' {
            normalized.push('N');
        } else if b"ACGTRYSWKMBDHVN".contains(&upper) {
            normalized.push(upper as char);
        } else {
            return Err(unsupported(format!(
                "unsupported trace call symbol 0x{b:02x}; no bases were invented or removed"
            )));
        }
    }
    Sequence::dna(&normalized)
}
impl Trace {
    pub fn validate(&self) -> Result<(), Error> {
        if self.name.trim().is_empty() || self.name.len() > 1024 {
            return Err(format_error("invalid trace name"));
        }
        let samples = self.channels.len();
        let n = self.sequence.len();
        if samples == 0
            || samples > MAX_TRACE_SAMPLES
            || n > MAX_TRACE_BASES
            || self.original_bytes.len() > MAX_IMPORT_BYTES
        {
            return Err(Error::TooLarge);
        }
        if self.sequence.alphabet() != Alphabet::Dna
            || self.channels.all().iter().any(|v| v.len() != samples)
            || self.peak_positions.len() != n
            || self.called_symbols.len() != n
        {
            return Err(format_error("trace channel/call/peak dimensions disagree"));
        }
        if self.peak_positions.iter().any(|&p| p as usize >= samples)
            || self.peak_positions.windows(2).any(|p| p[0] > p[1])
        {
            return Err(format_error(
                "trace peak indexes are out of bounds or descending",
            ));
        }
        match (&self.format, &self.qualities) {
            (TraceFormat::Abif { version }, _) if version / 100 != 1 => {
                return Err(unsupported("unsupported ABIF trace-model version"));
            }
            (TraceFormat::Abif { .. }, Some(QualityEvidence::ScfAccuracy(_)))
            | (TraceFormat::Scf { .. }, Some(QualityEvidence::AbifBasecallerScores(_))) => {
                return Err(format_error(
                    "quality evidence has the wrong source-format semantics",
                ));
            }
            _ => {}
        }
        let expected_magic = match self.format {
            TraceFormat::Abif { .. } => b"ABIF",
            TraceFormat::Scf { .. } => b".scf",
        };
        if self.original_bytes.len() < 128 || !self.original_bytes.starts_with(expected_magic) {
            return Err(format_error(
                "trace attachment does not match its declared source format",
            ));
        }
        let dash = matches!(self.format, TraceFormat::Scf { code_set: 0, .. });
        if calls(&self.called_symbols, dash)? != self.sequence {
            return Err(format_error(
                "trace call symbols disagree with the normalized sequence",
            ));
        }
        match &self.qualities {
            Some(QualityEvidence::AbifBasecallerScores(v)) if v.len() != n => {
                return Err(format_error("ABIF quality count differs from base count"));
            }
            Some(QualityEvidence::ScfAccuracy(q))
                if [&q.a, &q.c, &q.g, &q.t].iter().any(|v| v.len() != n) =>
            {
                return Err(format_error("SCF accuracy count differs from base count"));
            }
            _ => {}
        }
        if blake3::hash(&self.original_bytes).to_hex().as_str() != self.source_hash {
            return Err(format_error("trace original-evidence hash mismatch"));
        }
        Ok(())
    }
    /// A called-sequence projection only: this molecule does not contain the chromatogram.
    pub fn molecule(&self) -> Result<Molecule, Error> {
        self.validate()?;
        let mut m = Molecule::dna(&self.name, self.sequence.as_str(), Topology::Linear)?;
        m.metadata
            .insert("trace_source_hash".into(), self.source_hash.clone());
        m.metadata.insert(
            "trace_format".into(),
            match self.format {
                TraceFormat::Abif { .. } => "ABIF",
                TraceFormat::Scf { .. } => "SCF",
            }
            .into(),
        );
        m.metadata.insert("trace_projection".into(),"Called sequence only. Preserve the original trace attachment for signals, peaks, scores and metadata.".into());
        Ok(m)
    }
    /// Lossless original-attachment path for the current single-document native format.
    /// A full native serialization, including base64 overhead, must fit the 32 MiB budget.
    pub fn attached_molecule(&self) -> Result<Molecule, Error> {
        self.validate()?;
        let source = import_trace(&self.original_bytes)?;
        let expected = match self.orientation {
            TraceOrientation::AsStored => source,
            TraceOrientation::ReverseComplement => source.reverse_complement()?,
        };
        if self.format != expected.format
            || self.sequence != expected.sequence
            || self.called_symbols != expected.called_symbols
            || self.peak_positions != expected.peak_positions
            || self.channels != expected.channels
            || self.qualities != expected.qualities
        {
            return Err(format_error(
                "trace view differs from its original evidence; edited calls/signals require a separate derived record",
            ));
        }
        let encoded_len = self
            .original_bytes
            .len()
            .div_ceil(3)
            .checked_mul(4)
            .ok_or(Error::TooLarge)?;
        if encoded_len > MAX_IMPORT_BYTES {
            return Err(Error::TooLarge);
        }
        let mut molecule = self.molecule()?;
        molecule.metadata.remove("trace_projection");
        molecule.metadata.insert(
            "trace_original_base64".into(),
            STANDARD.encode(&self.original_bytes),
        );
        molecule.metadata.insert(
            "trace_orientation".into(),
            match self.orientation {
                TraceOrientation::AsStored => "as_stored",
                TraceOrientation::ReverseComplement => "reverse_complement",
            }
            .into(),
        );
        molecule.metadata.insert("trace_attachment_note".into(),"Original chromatogram retained unchanged. Signal/peak/score views must be reconstructed from the attachment and checked against current calls.".into());
        if crate::protocol::save_document(&molecule)?.len() > MAX_IMPORT_BYTES {
            return Err(Error::TooLarge);
        }
        Ok(molecule)
    }
    /// Explicit view transformation; original bytes, source hash and raw metadata stay unchanged.
    pub fn reverse_complement(&self) -> Result<Self, Error> {
        self.validate()?;
        let mut view = self.clone();
        view.sequence = self.sequence.reverse_complement()?;
        view.called_symbols = self
            .called_symbols
            .iter()
            .rev()
            .map(|&b| {
                if b == b'-' {
                    b'-'
                } else {
                    let c = collagene_model::complement(b.to_ascii_uppercase(), Alphabet::Dna);
                    if b.is_ascii_lowercase() {
                        c.to_ascii_lowercase()
                    } else {
                        c
                    }
                }
            })
            .collect();
        view.peak_positions = self
            .peak_positions
            .iter()
            .rev()
            .map(|&p| (self.channels.len() - 1) as u32 - p)
            .collect();
        let reversed = |v: &Vec<i32>| v.iter().rev().copied().collect();
        view.channels = TraceChannels {
            a: reversed(&self.channels.t),
            c: reversed(&self.channels.g),
            g: reversed(&self.channels.c),
            t: reversed(&self.channels.a),
        };
        view.qualities = match &self.qualities {
            Some(QualityEvidence::AbifBasecallerScores(v)) => Some(
                QualityEvidence::AbifBasecallerScores(v.iter().rev().copied().collect()),
            ),
            Some(QualityEvidence::ScfAccuracy(q)) => {
                let reverse = |v: &Vec<u8>| v.iter().rev().copied().collect();
                Some(QualityEvidence::ScfAccuracy(TraceAccuracy {
                    a: reverse(&q.t),
                    c: reverse(&q.g),
                    g: reverse(&q.c),
                    t: reverse(&q.a),
                }))
            }
            None => None,
        };
        view.orientation = match self.orientation {
            TraceOrientation::AsStored => TraceOrientation::ReverseComplement,
            TraceOrientation::ReverseComplement => TraceOrientation::AsStored,
        };
        view.validate()?;
        Ok(view)
    }
}

pub fn import_trace(bytes: &[u8]) -> Result<Trace, Error> {
    source_limit(bytes)?;
    let trace = if bytes.starts_with(b"ABIF") {
        abif(bytes)?
    } else if bytes.starts_with(b".scf") {
        scf(bytes)?
    } else {
        return Err(unsupported("not a supported ABIF or SCF trace"));
    };
    trace.validate()?;
    Ok(trace)
}

/// Download path for immutable original evidence, independent of current called-sequence edits.
/// Checks encoding, size, trace signature and stored source hash; it does not reinterpret calls.
pub fn original_trace_bytes(molecule: &Molecule) -> Result<Vec<u8>, Error> {
    let encoded = molecule
        .metadata
        .get("trace_original_base64")
        .ok_or_else(|| unsupported("document has no original trace attachment"))?;
    if encoded.len() > MAX_IMPORT_BYTES {
        return Err(Error::TooLarge);
    }
    let bytes = STANDARD
        .decode(encoded)
        .map_err(|_| format_error("trace attachment is not valid base64"))?;
    let hash = molecule
        .metadata
        .get("trace_source_hash")
        .ok_or_else(|| format_error("trace attachment source hash is missing"))?;
    if blake3::hash(&bytes).to_hex().as_str() != hash {
        return Err(format_error("trace attachment source hash mismatch"));
    }
    if bytes.len() < 128 || !(bytes.starts_with(b"ABIF") || bytes.starts_with(b".scf")) {
        return Err(format_error(
            "original trace attachment has no supported trace signature",
        ));
    }
    source_limit(&bytes)?;
    Ok(bytes)
}

/// Reconstruct attached observed evidence, checking exact source hash, view orientation and calls.
/// Sequence edits do not rewrite an immutable trace; stale projections are rejected explicitly.
pub fn trace_from_molecule(molecule: &Molecule) -> Result<Trace, Error> {
    molecule.validate()?;
    if molecule.topology != Topology::Linear || molecule.sequence.alphabet() != Alphabet::Dna {
        return Err(unsupported(
            "trace views require a linear DNA called-sequence document",
        ));
    }
    let bytes = original_trace_bytes(molecule)?;
    let source = import_trace(&bytes)?;
    let mut trace = match molecule
        .metadata
        .get("trace_orientation")
        .map(String::as_str)
    {
        Some("as_stored") => source,
        Some("reverse_complement") => source.reverse_complement()?,
        _ => {
            return Err(unsupported(
                "trace attachment orientation is missing or unsupported",
            ));
        }
    };
    if trace.sequence != molecule.sequence {
        return Err(format_error(
            "current sequence differs from the immutable trace attachment; retain it as original evidence and use a separate derived sequence record",
        ));
    }
    trace.name = molecule.name.clone();
    trace.validate()?;
    Ok(trace)
}

#[derive(Clone)]
struct Entry<'a> {
    kind: u16,
    size: usize,
    count: usize,
    data: &'a [u8],
}
impl Entry<'_> {
    fn bytes(&self, kinds: &[u16]) -> Result<&[u8], Error> {
        if !kinds.contains(&self.kind) || self.size != 1 || self.count != self.data.len() {
            return Err(unsupported(
                "unsupported ABIF byte/character entry type or dimensions",
            ));
        }
        Ok(self.data)
    }
    fn numbers(&self) -> Result<Vec<i32>, Error> {
        let expected = match self.kind {
            3 | 4 => 2,
            5 => 4,
            _ => return Err(unsupported("unsupported ABIF integer entry type")),
        };
        if self.size != expected || product(self.count, self.size)? != self.data.len() {
            return Err(format_error("ABIF integer entry dimensions disagree"));
        }
        if self.count > MAX_TRACE_SAMPLES {
            return Err(Error::TooLarge);
        }
        self.data
            .chunks_exact(expected)
            .map(|c| {
                Ok(match self.kind {
                    3 => u16::from_be_bytes(c.try_into().unwrap()) as i32,
                    4 => i16::from_be_bytes(c.try_into().unwrap()) as i32,
                    5 => i32::from_be_bytes(c.try_into().unwrap()),
                    _ => unreachable!(),
                })
            })
            .collect()
    }
    fn string(&self) -> Result<String, Error> {
        if self.size != 1 || self.count != self.data.len() {
            return Err(format_error("ABIF text entry dimensions disagree"));
        }
        let data = match self.kind {
            18 => {
                let (&len, rest) = self
                    .data
                    .split_first()
                    .ok_or_else(|| format_error("empty ABIF Pascal string"))?;
                if len as usize != rest.len() {
                    return Err(format_error("ABIF Pascal string length mismatch"));
                }
                rest
            }
            19 => self
                .data
                .strip_suffix(&[0])
                .ok_or_else(|| format_error("unterminated ABIF C string"))?,
            2 => self.data,
            _ => return Err(unsupported("unsupported ABIF text type")),
        };
        if !data.is_ascii() || data.contains(&0) {
            return Err(unsupported("ABIF metadata string is not supported ASCII"));
        }
        String::from_utf8(data.to_vec()).map_err(|_| format_error("invalid ABIF string"))
    }
}
fn abif(bytes: &[u8]) -> Result<Trace, Error> {
    span(bytes, 0, 128)?;
    let version = u16be(bytes, 4)?;
    if version / 100 != 1 {
        return Err(unsupported(format!(
            "ABIF major version {} is not implemented",
            version / 100
        )));
    }
    // The vendor explicitly says to obtain the root count/offset and ignore its other fields.
    let count = nonnegative(bytes, 18)?;
    let offset = nonnegative(bytes, 26)?;
    if count == 0 || count > MAX_ABIF_ENTRIES {
        return Err(Error::TooLarge);
    }
    let directory_length = product(count, 28)?;
    if offset < 128 {
        return Err(format_error("ABIF directory overlaps its header"));
    }
    span(bytes, offset, directory_length)?;
    let directory_end = offset + directory_length;
    let mut entries = BTreeMap::new();
    for i in 0..count {
        let p = offset + i * 28;
        let tag: [u8; 4] = span(bytes, p, 4)?.try_into().unwrap();
        if !tag.iter().all(|&b| (0x21..=0x7e).contains(&b)) {
            return Err(format_error("invalid ABIF tag name"));
        }
        let number = i32::from_be_bytes(span(bytes, p + 4, 4)?.try_into().unwrap());
        let size = u16be(bytes, p + 10)? as usize;
        let count = nonnegative(bytes, p + 12)?;
        let length = nonnegative(bytes, p + 16)?;
        let data = if length <= 4 {
            span(bytes, p + 20, length)?
        } else {
            let data_offset = nonnegative(bytes, p + 20)?;
            let data_end = data_offset.checked_add(length).ok_or(Error::TooLarge)?;
            if data_offset < 128 || (data_offset < directory_end && data_end > offset) {
                return Err(format_error("ABIF payload overlaps header or directory"));
            }
            span(bytes, data_offset, length)?
        };
        if entries
            .insert(
                (tag, number),
                Entry {
                    kind: u16be(bytes, p + 8)?,
                    size,
                    count,
                    data,
                },
            )
            .is_some()
        {
            return Err(format_error("duplicate ABIF tag and number"));
        }
    }
    let get = |name: [u8; 4], n: i32| {
        entries.get(&(name, n)).ok_or_else(|| {
            unsupported(format!(
                "ABIF requires {} {n}",
                String::from_utf8_lossy(&name)
            ))
        })
    };
    let mut warnings=vec!["ABIF intensities and caller scores are preserved without rescaling, quality calibration or automatic trimming.".into()];
    if let Some(flag) = entries.get(&(*b"RevC", 1)) {
        let values = flag.numbers()?;
        if values.len() != 1 {
            return Err(format_error("ABIF RevC flag must be a scalar"));
        }
        if values[0] != 0 {
            return Err(unsupported(
                "ABIF RevC is set; stored call/trace orientation semantics require an audited adapter before import",
            ));
        }
    }
    let edition = if entries.contains_key(&(*b"PBAS", 2)) {
        2
    } else {
        1
    };
    let raw = get(*b"PBAS", edition)?.bytes(&[2])?.to_vec();
    let sequence = calls(&raw, false)?;
    let positions = get(*b"PLOC", edition)?.numbers()?;
    let peak_positions = positions
        .into_iter()
        .map(|p| u32::try_from(p).map_err(|_| format_error("negative ABIF peak position")))
        .collect::<Result<Vec<_>, _>>()?;
    let order = get(*b"FWO_", 1)?.bytes(&[2])?;
    if order.len() != 4
        || b"ACGT"
            .iter()
            .any(|base| order.iter().filter(|&&b| b == *base).count() != 1)
    {
        return Err(unsupported(
            "ABIF FWO_ is not an unambiguous four-base channel permutation",
        ));
    }
    let mut channels: [Vec<i32>; 4] = std::array::from_fn(|_| Vec::new());
    for (i, &base) in order.iter().enumerate() {
        let entry = get(*b"DATA", 9 + i as i32)?;
        if !matches!(entry.kind, 3 | 4) {
            return Err(unsupported(
                "ABIF analyzed DATA channels require 16-bit word or short elements",
            ));
        }
        let values = entry.numbers()?;
        let index = match base {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => unreachable!(),
        };
        channels[index] = values;
    }
    let [a, c, g, t] = channels;
    let qualities = match entries.get(&(*b"PCON", edition)) {
        Some(q) => Some(QualityEvidence::AbifBasecallerScores(
            q.bytes(&[1, 2])?.to_vec(),
        )),
        None => {
            warnings.push("No PCON scores matching the selected call edition were present; no qualities were inferred.".into());
            None
        }
    };
    if edition == 2 && entries.contains_key(&(*b"PBAS", 1)) {
        warnings.push("PBAS/PLOC edition 2 was selected. Edited edition 1 data remains in the original attachment and was not mixed with edition 2.".into());
    }
    if peak_positions.windows(2).any(|p| p[0] == p[1]) {
        warnings.push("Coincident called-peak positions were retained.".into());
    }
    let mut metadata = BTreeMap::new();
    metadata.insert("abif_call_edition".into(), edition.to_string());
    metadata.insert(
        "abif_channel_order".into(),
        String::from_utf8(order.to_vec()).unwrap(),
    );
    metadata.insert("abif_directory_entries".into(), count.to_string());
    let name = if let Some(entry) = entries.get(&(*b"SMPL", 1)) {
        match entry.string() {
            Ok(name) if !name.trim().is_empty() && name.len() <= 1024 => name,
            _ => {
                warnings.push(
                    "Unsupported sample-name metadata remains in the original attachment.".into(),
                );
                "ABIF trace".into()
            }
        }
    } else {
        "ABIF trace".into()
    };
    Ok(Trace {
        name,
        format: TraceFormat::Abif { version },
        orientation: TraceOrientation::AsStored,
        sequence,
        called_symbols: raw,
        peak_positions,
        channels: TraceChannels { a, c, g, t },
        qualities,
        metadata,
        warnings,
        original_bytes: bytes.to_vec(),
        source_hash: blake3::hash(bytes).to_hex().to_string(),
    })
}

fn scf(bytes: &[u8]) -> Result<Trace, Error> {
    span(bytes, 0, 128)?;
    let version = span(bytes, 36, 4)?;
    if !matches!(version, b"2.00" | b"2.01" | b"2.02" | b"3.00") {
        return Err(unsupported(
            "SCF version is not one of 2.00, 2.01, 2.02 or 3.00",
        ));
    }
    let major3 = version[0] == b'3';
    let samples = u32be(bytes, 4)? as usize;
    let bases = u32be(bytes, 12)? as usize;
    if samples == 0 || bases == 0 || samples > MAX_TRACE_SAMPLES || bases > MAX_TRACE_BASES {
        return Err(Error::TooLarge);
    }
    let declared_precision = u32be(bytes, 40)?;
    let declared_code = u32be(bytes, 44)?;
    let (precision, code_set) = if version == b"2.00" {
        if !matches!(declared_precision, 0 | 1) || declared_code != 0 {
            return Err(unsupported(
                "SCF 2.00 requires the documented default 8-bit samples and default call code set",
            ));
        }
        (1, 0)
    } else {
        (declared_precision, declared_code)
    };
    if !matches!(precision, 1 | 2) {
        return Err(unsupported("SCF sample precision must be one or two bytes"));
    }
    if !matches!(code_set, 0 | 2 | 4) {
        return Err(unsupported(format!(
            "SCF uncertainty code set {code_set} is not implemented"
        )));
    }
    let sample_offset = u32be(bytes, 8)? as usize;
    let base_offset = u32be(bytes, 24)? as usize;
    let comment_len = u32be(bytes, 28)? as usize;
    let comment_offset = u32be(bytes, 32)? as usize;
    let private_len = u32be(bytes, 48)? as usize;
    let private_offset = u32be(bytes, 52)? as usize;
    let sample_len = product(product(samples, 4)?, precision as usize)?;
    let base_len = product(bases, 12)?;
    let mut occupied = Vec::new();
    for (offset, length) in [
        (sample_offset, sample_len),
        (base_offset, base_len),
        (comment_offset, comment_len),
        (private_offset, private_len),
    ] {
        if length == 0 {
            continue;
        }
        if offset < 128 {
            return Err(format_error("SCF section overlaps its header"));
        }
        span(bytes, offset, length)?;
        let end = offset.checked_add(length).ok_or(Error::TooLarge)?;
        if occupied
            .iter()
            .any(|&(start, stop)| offset < stop && end > start)
        {
            return Err(format_error("SCF sections overlap"));
        }
        occupied.push((offset, end));
    }
    let raw_samples = span(bytes, sample_offset, sample_len)?;
    let mut channels: [Vec<i32>; 4] = std::array::from_fn(|_| Vec::with_capacity(samples));
    for (channel, values) in channels.iter_mut().enumerate() {
        for sample in 0..samples {
            let index = if major3 {
                channel * samples + sample
            } else {
                sample * 4 + channel
            };
            let offset = index * precision as usize;
            values.push(if precision == 1 {
                raw_samples[offset] as i32
            } else {
                u16be(raw_samples, offset)? as i32
            });
        }
        if major3 {
            let mask = if precision == 1 { 255 } else { 65535 };
            // Two prefix sums invert the format's modulo-width second differences.
            for _ in 0..2 {
                let mut previous = 0;
                for value in values.iter_mut() {
                    previous = (previous + *value) & mask;
                    *value = previous;
                }
            }
        }
    }
    let mut peak_positions = Vec::with_capacity(bases);
    let mut called_symbols = Vec::with_capacity(bases);
    let mut scores: [Vec<u8>; 4] = std::array::from_fn(|_| Vec::with_capacity(bases));
    let base_data = span(bytes, base_offset, base_len)?;
    for i in 0..bases {
        let peak = u32be(base_data, if major3 { i * 4 } else { i * 12 })?;
        peak_positions.push(peak);
        for (channel, values) in scores.iter_mut().enumerate() {
            values.push(
                base_data[if major3 {
                    bases * (4 + channel) + i
                } else {
                    i * 12 + 4 + channel
                }],
            );
        }
        called_symbols.push(base_data[if major3 { bases * 8 + i } else { i * 12 + 8 }]);
    }
    // Do not reinterpret an undeclared legacy encoding as IUPAC.
    for &b in &called_symbols {
        let b = b.to_ascii_uppercase();
        let permitted = match code_set {
            0 => b"ACGT-".contains(&b),
            2 => b"ACGTRYSWKMBDHVN".contains(&b),
            4 => b"ACGTN".contains(&b),
            _ => false,
        };
        if !permitted {
            return Err(unsupported(
                "SCF call symbol is incompatible with its declared uncertainty code set",
            ));
        }
    }
    let sequence = calls(&called_symbols, code_set == 0)?;
    let [a, c, g, t] = channels;
    let [qa, qc, qg, qt] = scores;
    let mut warnings=vec!["SCF stores four accuracy estimates per call; no Phred scale, automatic trimming or recalibration is inferred.".into()];
    if called_symbols.contains(&b'-') {
        warnings.push("SCF '-' no-call symbols were preserved and projected to N in the DNA sequence, without deleting positions.".into());
    }
    if peak_positions.windows(2).any(|p| p[0] == p[1]) {
        warnings.push("Coincident called-peak positions were retained.".into());
    }
    let mut metadata = BTreeMap::new();
    metadata.insert(
        "scf_reported_left_clip_obsolete".into(),
        u32be(bytes, 16)?.to_string(),
    );
    metadata.insert(
        "scf_reported_right_clip_obsolete".into(),
        u32be(bytes, 20)?.to_string(),
    );
    metadata.insert("scf_private_bytes".into(), private_len.to_string());
    metadata.insert("scf_comment_bytes".into(), comment_len.to_string());
    let mut name = "SCF trace".to_owned();
    if comment_len > 0 {
        let comments = span(bytes, comment_offset, comment_len)?;
        let content = comments.split(|&b| b == 0).next().unwrap_or_default();
        if let Ok(text) = std::str::from_utf8(content) {
            if let Some(value) = text.lines().find_map(|line| line.strip_prefix("NAME=")) {
                if !value.trim().is_empty() && value.len() <= 1024 {
                    name = value.into();
                }
            }
        } else {
            warnings.push(
                "Non-UTF-8 comment bytes were retained only in the original attachment.".into(),
            );
        }
    }
    Ok(Trace {
        name,
        format: TraceFormat::Scf {
            version: String::from_utf8(version.to_vec()).unwrap(),
            code_set,
        },
        orientation: TraceOrientation::AsStored,
        sequence,
        called_symbols,
        peak_positions,
        channels: TraceChannels { a, c, g, t },
        qualities: Some(QualityEvidence::ScfAccuracy(TraceAccuracy {
            a: qa,
            c: qc,
            g: qg,
            t: qt,
        })),
        metadata,
        warnings,
        original_bytes: bytes.to_vec(),
        source_hash: blake3::hash(bytes).to_hex().to_string(),
    })
}
