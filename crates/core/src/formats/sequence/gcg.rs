//! Independent bounded Wisconsin/GCG single-sequence exchange profile.
//! Text evidence never creates native cloning history or chemical authority.
use crate::formats::sequence::{ExportReport, ImportReport};
use base64::{Engine as _, engine::general_purpose::STANDARD};
use collagene_model::{Alphabet, Error, Molecule, Sequence, Topology};
use serde::{Deserialize, Serialize};
use std::{collections::BTreeMap, fmt::Write};
pub const MAX_BYTES: usize = 4 * 1024 * 1024;
pub const MAX_SEQUENCE: usize = 1024 * 1024;
const MAX_HEADER: usize = 65536;
const CHUNK: usize = 65536;
const ORIGINAL: &str = "gcg:original:";
#[derive(Debug, Clone, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct ImportOptions {
    pub alphabet: Alphabet,
    pub fallback_name: String,
}
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct ExportOptions {
    pub record_id: Option<String>,
    pub date_label: Option<String>,
}
fn bad(message: &str) -> Error {
    Error::Format(format!("GCG single sequence: {message}"))
}
/// The documented Wisconsin checksum: uppercase sequence byte times a cycling
/// weight 1..57, summed modulo 10000. No whitespace/gap stripping occurs here.
pub fn checksum(sequence: &str) -> Result<u16, Error> {
    let mut sum = 0u32;
    for (i, b) in sequence.bytes().enumerate() {
        if !(b.is_ascii_alphabetic() || b == b'*') {
            return Err(bad("checksum requires ungapped ASCII sequence calls"));
        }
        sum = (sum + ((i % 57 + 1) as u32) * u32::from(b.to_ascii_uppercase())) % 10000;
    }
    Ok(sum as u16)
}
fn numeric(s: &str) -> Result<usize, Error> {
    if s.is_empty() || !s.bytes().all(|b| b.is_ascii_digit()) {
        return Err(bad("invalid unsigned header/position number"));
    }
    s.parse()
        .map_err(|_| bad("header/position number overflow"))
}
fn date_label(s: &str) -> Result<(), Error> {
    if s.trim().is_empty()
        || s.len() > 128
        || s.chars().any(char::is_control)
        || ["..", "Type:", "Length:", "Check:", "Checksum:"]
            .iter()
            .any(|p| s.contains(p))
    {
        return Err(bad("a bounded explicit documentary date label is required"));
    }
    Ok(())
}
fn record_id(s: &str) -> Result<(), Error> {
    if s.is_empty()
        || s.len() > 256
        || !s.bytes().all(|b| b.is_ascii_graphic())
        || ["..", "Length:", "Type:", "Check:"]
            .iter()
            .any(|p| s.contains(p))
    {
        return Err(bad(
            "record ID must be a bounded non-whitespace ASCII token",
        ));
    }
    Ok(())
}
struct Header<'a> {
    id: &'a str,
    date: &'a str,
    length: usize,
    kind: &'a str,
    check: u16,
}
fn header(line: &str) -> Result<Header<'_>, Error> {
    let line = line
        .trim()
        .strip_suffix("..")
        .ok_or_else(|| bad("missing dividing-line terminator"))?
        .trim_end();
    if line.matches("Length:").count() != 1 || line.matches("Type:").count() != 1 {
        return Err(bad("exactly one Length and Type are required"));
    }
    let (id, after) = line
        .split_once("Length:")
        .ok_or_else(|| bad("missing Length"))?;
    let (length_date, typed) = after
        .split_once("Type:")
        .ok_or_else(|| bad("missing explicit Type: N or P; no alphabet inference"))?;
    let length_date = length_date.trim();
    let boundary = length_date
        .find(char::is_whitespace)
        .unwrap_or(length_date.len());
    let length = numeric(&length_date[..boundary])?;
    if length == 0 || length > MAX_SEQUENCE {
        return Err(Error::TooLarge);
    }
    let date = length_date[boundary..].trim();
    if !date.is_empty() {
        date_label(date)?;
    }
    let (kind, check) = if typed.matches("Check:").count() == 1 && !typed.contains("Checksum:") {
        typed.split_once("Check:").unwrap()
    } else if typed.matches("Checksum:").count() == 1 && !typed.contains("Check:") {
        typed.split_once("Checksum:").unwrap()
    } else {
        return Err(bad("exactly one Check or Checksum field is required"));
    };
    let kind = kind.trim();
    if !matches!(kind, "N" | "P") {
        return Err(bad("Type must be N or P"));
    }
    let check = numeric(check.trim())?;
    if check >= 10000 {
        return Err(bad("checksum must be between 0 and 9999"));
    }
    let id = id.trim();
    if !id.is_empty() {
        record_id(id)?;
    }
    Ok(Header {
        id,
        date,
        length,
        kind,
        check: check as u16,
    })
}
fn retain(m: &mut Molecule, bytes: &[u8]) {
    let encoded = STANDARD.encode(bytes);
    m.metadata.insert(
        format!("{ORIGINAL}chunks"),
        encoded.len().div_ceil(CHUNK).to_string(),
    );
    m.metadata.insert(
        format!("{ORIGINAL}hash"),
        blake3::hash(bytes).to_hex().to_string(),
    );
    for (i, chunk) in encoded.as_bytes().chunks(CHUNK).enumerate() {
        m.metadata.insert(
            format!("{ORIGINAL}{i:04}"),
            std::str::from_utf8(chunk).expect("base64").into(),
        );
    }
}
pub fn original_bytes(m: &Molecule) -> Result<Vec<u8>, Error> {
    let count: usize = m
        .metadata
        .get(&format!("{ORIGINAL}chunks"))
        .ok_or_else(|| bad("no exact original retained"))?
        .parse()
        .map_err(|_| bad("original chunk count"))?;
    if count == 0 || count > MAX_BYTES.div_ceil(3) * 4 / CHUNK + 1 {
        return Err(Error::TooLarge);
    }
    let mut encoded = String::new();
    for i in 0..count {
        let part = m
            .metadata
            .get(&format!("{ORIGINAL}{i:04}"))
            .ok_or_else(|| bad("missing original chunk"))?;
        if part.len() > CHUNK
            || encoded
                .len()
                .checked_add(part.len())
                .is_none_or(|n| n > MAX_BYTES.div_ceil(3) * 4)
        {
            return Err(Error::TooLarge);
        }
        encoded.push_str(part);
    }
    let bytes = STANDARD
        .decode(encoded)
        .map_err(|_| bad("original base64"))?;
    if bytes.len() > MAX_BYTES
        || m.metadata.get(&format!("{ORIGINAL}hash"))
            != Some(&blake3::hash(&bytes).to_hex().to_string())
    {
        return Err(bad("original hash mismatch"));
    }
    Ok(bytes)
}
/// Narrow header admission only. Never chooses DNA versus RNA from Type: N,
/// and never treats GCG multiple-alignment magic as a single sequence.
pub fn is_gcg(bytes: &[u8]) -> bool {
    if bytes.is_empty() || bytes.len() > MAX_BYTES {
        return false;
    }
    let prefix = &bytes[..bytes.len().min(MAX_HEADER)];
    let Ok(text) = std::str::from_utf8(prefix) else {
        return false;
    };
    let mut lines = text.lines();
    if let Some(first) = lines.find(|line| !line.trim().is_empty()) {
        if matches!(first.trim_end(), "!!NA_SEQUENCE 1.0" | "!!AA_SEQUENCE 1.0") {
            return true;
        }
        if first.starts_with("!!") {
            return false;
        }
        if header(first).is_ok() {
            return true;
        }
    }
    lines.any(|line| line.trim_end().ends_with("..") && header(line).is_ok())
}
pub fn import_gcg(bytes: &[u8], options: &ImportOptions) -> Result<ImportReport, Error> {
    if bytes.is_empty() || bytes.len() > MAX_BYTES {
        return Err(Error::TooLarge);
    }
    if options.fallback_name.len() > 1024 || options.fallback_name.chars().any(char::is_control) {
        return Err(bad("fallback name bounds"));
    }
    let text = std::str::from_utf8(bytes).map_err(|_| bad("UTF-8 text required"))?;
    let mut lines = text.lines();
    let mut magic = None;
    let mut preamble = 0usize;
    let mut found = None;
    for line in lines.by_ref() {
        preamble = preamble
            .checked_add(line.len() + 1)
            .filter(|n| *n <= MAX_HEADER)
            .ok_or(Error::TooLarge)?;
        if line.starts_with("!!") {
            let k = match line.trim_end() {
                "!!NA_SEQUENCE 1.0" => "N",
                "!!AA_SEQUENCE 1.0" => "P",
                _ => {
                    return Err(bad(
                        "unsupported GCG magic/version; alignments need the alignment reader",
                    ));
                }
            };
            if magic.replace(k).is_some() {
                return Err(bad("multiple sequence headers are unsupported"));
            }
        }
        if line.trim_end().ends_with("..") && line.contains("Length:") {
            found = Some(header(line)?);
            break;
        }
    }
    let h = found.ok_or_else(|| bad("no bounded typed sequence dividing line"))?;
    if magic.is_some_and(|k| k != h.kind) {
        return Err(bad("magic and Type disagree"));
    }
    if (h.kind == "P") != (options.alphabet == Alphabet::Protein) {
        return Err(bad("declared type and explicit alphabet disagree"));
    }
    let mut sequence = String::with_capacity(h.length);
    let mut rows = 0usize;
    for line in lines {
        rows += 1;
        if rows > 100000 || line.len() > 65536 {
            return Err(Error::TooLarge);
        }
        let mut tokens = line.split_ascii_whitespace().peekable();
        if let Some(first) = tokens.peek() {
            if first.bytes().all(|b| b.is_ascii_digit()) {
                let start = numeric(tokens.next().unwrap())?;
                if start != sequence.len() + 1 {
                    return Err(bad(
                        "numbered row does not start at the exact next one-based call",
                    ));
                }
                if tokens.peek().is_none() {
                    return Err(bad("numbered row has no sequence calls"));
                }
            }
        }
        for group in tokens {
            if !group.bytes().all(|b| b.is_ascii_alphabetic() || b == b'*') {
                return Err(bad(
                    "gaps, inline annotations, trailing positions and additional records are outside this profile",
                ));
            }
            if sequence
                .len()
                .checked_add(group.len())
                .is_none_or(|n| n > h.length)
            {
                return Err(bad(
                    "sequence exceeds declared length or contains multiple records",
                ));
            }
            sequence.push_str(group);
        }
    }
    if sequence.len() != h.length || checksum(&sequence)? != h.check {
        return Err(bad(
            "declared length or checksum disagrees with actual calls",
        ));
    }
    let name = if h.id.is_empty() {
        options.fallback_name.trim()
    } else {
        h.id
    };
    if name.is_empty() {
        return Err(bad(
            "record has no ID; an explicit fallback name is required",
        ));
    }
    let mut m = Molecule {
        name: name.into(),
        sequence: Sequence::new(options.alphabet, &sequence)?,
        topology: Topology::Linear,
        features: vec![],
        metadata: BTreeMap::new(),
        original_genbank: None,
    };
    m.metadata.insert("gcg:record_id".into(), h.id.into());
    m.metadata.insert("gcg:date_label".into(), h.date.into());
    m.metadata.insert("gcg:source_type".into(), h.kind.into());
    retain(&mut m, bytes);
    crate::formats::sequence::retain_digest_producer_comments(&mut m, text);
    m.validate()?;
    collagene_model::bounds::model_preflight(&[&m])?;
    let mut warnings=vec!["GCG imports sequence calls only. Linear display does not establish original topology; documentary headers remain exact native evidence, not mapped features or cloning history.".into()];
    if magic.is_none() {
        warnings.push("Typed legacy GCG dividing line imported without modern file magic.".into());
    }
    if h.date.is_empty() {
        warnings.push(
            "Source has no documentary date; regenerated GCG requires an explicit date label."
                .into(),
        );
    }
    Ok(ImportReport {
        molecules: vec![m],
        warnings,
    })
}
pub fn export_gcg(m: &Molecule, options: &ExportOptions) -> Result<ExportReport, Error> {
    collagene_model::bounds::model_preflight(&[m])?;
    m.validate()?;
    if m.sequence.is_empty() || m.sequence.len() > MAX_SEQUENCE {
        return Err(Error::TooLarge);
    }
    // Share the same restrictive-evidence gate as the selected SnapGene writer;
    // unlike declaration preflight this does not restrict ordinary RNA/protein,
    // circular reference topology, or legitimate ambiguity calls.
    crate::formats::fidelity::reference_export_guard(m)?;
    if m.sequence.alphabet() == Alphabet::Protein
        && crate::protein::properties::has_modification_evidence(m)
    {
        return Err(Error::Unsupported("GCG cannot retain protein modification/crosslink evidence when dropping annotations; keep the native protein archive".into()));
    }
    if m.metadata
        .keys()
        .any(|key| key.starts_with("recombination:") || key.starts_with("source:recombination:"))
    {
        return Err(Error::Unsupported(
            "GCG cannot preserve restrictive reference-only ancestry; retain the native archive"
                .into(),
        ));
    }
    let id = options
        .record_id
        .as_deref()
        .or_else(|| m.metadata.get("gcg:record_id").map(String::as_str))
        .filter(|s| !s.is_empty())
        .ok_or_else(|| bad("declare a record ID for synthetic export"))?;
    record_id(id)?;
    let date = options
        .date_label
        .as_deref()
        .or_else(|| m.metadata.get("gcg:date_label").map(String::as_str))
        .ok_or_else(|| {
            bad("declare an actual export date label; no core clock or placeholder is invented")
        })?;
    date_label(date)?;
    let kind = if m.sequence.alphabet() == Alphabet::Protein {
        "P"
    } else {
        "N"
    };
    let magic = if kind == "P" {
        "!!AA_SEQUENCE 1.0"
    } else {
        "!!NA_SEQUENCE 1.0"
    };
    let mut text = format!(
        "{magic}\n\n{id}  Length: {}  {date}  Type: {kind}  Check: {}  ..\n\n",
        m.sequence.len(),
        checksum(m.sequence.as_str())?
    );
    for (start, row) in m.sequence.as_str().as_bytes().chunks(50).enumerate() {
        write!(&mut text, "{:>9}  ", start * 50 + 1).map_err(|_| bad("export text"))?;
        for (i, group) in row.chunks(10).enumerate() {
            if i > 0 {
                text.push(' ');
            }
            text.push_str(std::str::from_utf8(group).expect("ASCII native calls"));
        }
        text.push('\n');
    }
    if text.len() > MAX_BYTES {
        return Err(Error::TooLarge);
    }
    let mut warnings=vec!["GCG is a single-sequence projection: native features, metadata, originals, primers, chemistry and cloning history are omitted. Keep the native archive for complete evidence.".into()];
    if m.topology == Topology::Circular {
        warnings.push("Circular topology is not encoded by this GCG profile; reimport requires a separate explicit topology decision.".into());
    }
    if m.sequence.alphabet() == Alphabet::Rna {
        warnings.push(
            "GCG Type: N does not distinguish RNA from DNA; select RNA explicitly when importing."
                .into(),
        );
    }
    Ok(ExportReport { text, warnings })
}
