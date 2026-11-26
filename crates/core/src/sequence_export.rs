//! Source-bound reviewed reference-file exports. No UI, paths, clocks, scientific
//! reconstruction or physical-authority restoration belongs in this facade.
use crate::{Alphabet, Error, Molecule, Topology, gcg_sequence, snapgene_export, xdna_export};
use serde::{Deserialize, Serialize};
use std::io::Write;
pub const ENGINE_VERSION: &str = "sequence-export-1";
pub const GCG_PROFILE: &str = "collagene-gcg-single-sequence-1";
pub const MAX_BYTES: usize = 32 * 1024 * 1024;
const MAX_LOSSES: usize = 1024;
const MAX_WARNINGS: usize = 256;
#[derive(Debug, Clone, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub enum Format {
    SnapGene,
    Xdna,
    Gcg {
        options: gcg_sequence::ExportOptions,
    },
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Loss {
    pub path: String,
    pub message: String,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Assessment {
    pub engine_version: String,
    pub context_revision: String,
    pub source_revision: String,
    pub profile: String,
    pub filename: String,
    pub mime: String,
    pub losses: Vec<Loss>,
    pub warnings: Vec<String>,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Prepared {
    pub assessment: Assessment,
    pub bytes: Vec<u8>,
}
fn bad(s: &str) -> Error {
    Error::Unsupported(format!("Sequence export: {s}"))
}
fn digest(s: &str) -> bool {
    s.len() == 64
        && s.bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}
fn source_preflight(m: &Molecule, format: &Format) -> Result<(), Error> {
    crate::model::bounds::model_preflight(&[m])?;
    m.validate()?;
    crate::format_fidelity::reference_export_guard(m)?;
    match format {
        Format::SnapGene => {
            if m.sequence.alphabet() == Alphabet::Rna
                || m.sequence.alphabet() == Alphabet::Protein && m.topology != Topology::Linear
            {
                return Err(bad(
                    "Selected SnapGene profile supports DNA or linear protein only.",
                ));
            }
        }
        Format::Xdna => {
            if m.sequence.alphabet() == Alphabet::Rna
                || m.sequence.alphabet() == Alphabet::Protein && m.topology != Topology::Linear
            {
                return Err(bad(
                    "Selected XDNA profile supports DNA or linear protein only.",
                ));
            }
            if m.features.len() > xdna_export::MAX_FEATURES {
                return Err(Error::TooLarge);
            }
            if ["legacy:xdna_left_length", "legacy:xdna_right_length"]
                .iter()
                .any(|key| {
                    m.metadata
                        .get(*key)
                        .is_some_and(|value| value.parse::<i64>().ok() != Some(0))
                })
            {
                return Err(bad(
                    "XDNA overhang evidence cannot be discarded into an intact reference.",
                ));
            }
        }
        Format::Gcg { options } => {
            if m.sequence.len() > gcg_sequence::MAX_SEQUENCE {
                return Err(Error::TooLarge);
            }
            if m.sequence.alphabet() == Alphabet::Protein
                && crate::protein_properties::has_modification_evidence(m)
            {
                return Err(bad(
                    "GCG cannot omit protein modification/crosslink evidence.",
                ));
            }
            if m.metadata.keys().any(|key| {
                key.starts_with("recombination:") || key.starts_with("source:recombination:")
            }) {
                return Err(bad("GCG cannot omit restrictive recombination evidence."));
            }
            // Mirror the existing writer's structural options, before source
            // hashing. The writer remains authoritative for its actual content.
            let id = options
                .record_id
                .as_deref()
                .or_else(|| m.metadata.get("gcg:record_id").map(String::as_str))
                .ok_or_else(|| bad("Declare a GCG record ID."))?;
            if id.is_empty()
                || id.len() > 256
                || !id.bytes().all(|b| b.is_ascii_graphic())
                || ["..", "Length:", "Type:", "Check:"]
                    .iter()
                    .any(|p| id.contains(p))
            {
                return Err(bad(
                    "GCG record ID must be a bounded ASCII token without header delimiters.",
                ));
            }
            let date = options
                .date_label
                .as_deref()
                .or_else(|| m.metadata.get("gcg:date_label").map(String::as_str))
                .ok_or_else(|| bad("Declare an explicit GCG documentary date label."))?;
            if date.trim().is_empty()
                || date.len() > 128
                || date.chars().any(char::is_control)
                || ["..", "Type:", "Length:", "Check:", "Checksum:"]
                    .iter()
                    .any(|p| date.contains(p))
            {
                return Err(bad(
                    "GCG date label is invalid or contains a header delimiter.",
                ));
            }
        }
    }
    if m.sequence.is_empty() {
        return Err(Error::TooLarge);
    }
    Ok(())
}
fn source_revision(m: &Molecule) -> Result<String, Error> {
    struct Hash {
        hash: blake3::Hasher,
        bytes: usize,
    }
    impl Write for Hash {
        fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
            self.bytes = self
                .bytes
                .checked_add(bytes.len())
                .filter(|n| *n <= 128 * 1024 * 1024)
                .ok_or_else(|| std::io::Error::other("source export binding exceeds128MiB"))?;
            self.hash.update(bytes);
            Ok(bytes.len())
        }
        fn flush(&mut self) -> std::io::Result<()> {
            Ok(())
        }
    }
    let mut writer = Hash {
        hash: blake3::Hasher::new(),
        bytes: 0,
    };
    serde_json::to_writer(&mut writer, m).map_err(|_| Error::TooLarge)?;
    Ok(writer.hash.finalize().to_hex().to_string())
}
fn context(source: &str, format: &Format) -> Result<String, Error> {
    crate::assembly_recipe::hash(&(
        ENGINE_VERSION,
        crate::verification::ENGINE_VERSION,
        source,
        format,
    ))
}
pub fn revision(m: &Molecule, format: &Format) -> Result<String, Error> {
    source_preflight(m, format)?;
    context(&source_revision(m)?, format)
}
fn safe_filename(name: &str, extension: &str, source: &str, warnings: &mut Vec<String>) -> String {
    let mut base = name
        .chars()
        .map(|c| {
            if c.is_control() || "/\\:*?\"<>|".contains(c) {
                '_'
            } else {
                c
            }
        })
        .collect::<String>();
    base = base.trim().trim_matches('.').to_string();
    if base.is_empty() || base.len() > 200 {
        base = format!("sequence-{}", &source[..16]);
    }
    let stem = base.split('.').next().unwrap_or("").to_ascii_uppercase();
    if matches!(stem.as_str(), "CON" | "PRN" | "AUX" | "NUL")
        || (stem.len() == 4
            && (stem.starts_with("COM") || stem.starts_with("LPT"))
            && matches!(stem.as_bytes()[3], b'1'..=b'9'))
    {
        base = format!("sequence-{base}");
    }
    if base != name {
        warnings.push("Portable output basename differs from the native display name: unsafe characters were replaced, reserved names prefixed, or an oversized/empty name replaced by a stable source-hash basename. No native document name was changed.".into());
    }
    format!("{base}{extension}")
}
fn assessment(
    m: &Molecule,
    format: &Format,
    source: String,
    profile: &str,
    extension: &str,
    mime: &str,
    losses: Vec<Loss>,
    mut warnings: Vec<String>,
) -> Result<Assessment, Error> {
    let filename = safe_filename(&m.name, extension, &source, &mut warnings);
    let result = Assessment {
        engine_version: ENGINE_VERSION.into(),
        context_revision: context(&source, format)?,
        source_revision: source,
        profile: profile.into(),
        filename,
        mime: mime.into(),
        losses,
        warnings,
    };
    result.validate_binding(&result.context_revision)?;
    Ok(result)
}
fn gcg_losses(
    m: &Molecule,
    options: &gcg_sequence::ExportOptions,
    warnings: &[String],
) -> Vec<Loss> {
    // Existing UI required explicit acknowledgement of every GCG writer warning.
    // Retain that contract; a migration must not reclassify omitted information
    // as a purely advisory status. No source values are printed into messages.
    let mut losses = warnings
        .iter()
        .enumerate()
        .map(|(i, message)| Loss {
            path: format!("gcg.warnings[{i}]"),
            message: message.clone(),
        })
        .collect::<Vec<_>>();
    let id = options
        .record_id
        .as_deref()
        .or_else(|| m.metadata.get("gcg:record_id").map(String::as_str));
    if id != Some(m.name.as_str()) {
        losses.push(Loss{path:"name".into(),message:"Explicit GCG record ID differs from the native display name; only that ID is embedded in the regenerated header.".into()});
    }
    losses
}
pub fn review(m: &Molecule, format: &Format) -> Result<Assessment, Error> {
    source_preflight(m, format)?;
    let source = source_revision(m)?;
    match format {
        Format::SnapGene => {
            let r = snapgene_export::review(m)?;
            if r.source_revision != source {
                return Err(bad("Writer source binding changed."));
            }
            assessment(
                m,
                format,
                source,
                &r.profile,
                &r.extension,
                "application/octet-stream",
                r.losses
                    .into_iter()
                    .map(|l| Loss {
                        path: l.path,
                        message: l.message,
                    })
                    .collect(),
                r.warnings,
            )
        }
        Format::Xdna => {
            let r = xdna_export::review(m)?;
            if r.source_revision != source {
                return Err(bad("Writer source binding changed."));
            }
            assessment(
                m,
                format,
                source,
                &r.profile,
                &r.extension,
                "application/octet-stream",
                r.losses
                    .into_iter()
                    .map(|l| Loss {
                        path: l.path,
                        message: l.message,
                    })
                    .collect(),
                r.warnings,
            )
        }
        Format::Gcg { options } => {
            let r = gcg_sequence::export_gcg(m, options)?;
            let losses = gcg_losses(m, options, &r.warnings);
            assessment(
                m,
                format,
                source,
                GCG_PROFILE,
                ".gcg",
                "text/plain;charset=utf-8",
                losses,
                r.warnings,
            )
        }
    }
}
pub fn prepare(
    m: &Molecule,
    format: &Format,
    expected_context_revision: &str,
    allow_omissions: bool,
) -> Result<Prepared, Error> {
    if !digest(expected_context_revision) {
        return Err(bad("Invalid expected export context."));
    }
    source_preflight(m, format)?;
    let source = source_revision(m)?;
    if context(&source, format)? != expected_context_revision {
        return Err(bad(
            "Source or export options changed; review the current export again.",
        ));
    }
    let (assessment, bytes) = match format {
        Format::SnapGene => {
            let r = snapgene_export::export_snapgene(m, snapgene_export::Policy::AllowOmissions)?;
            if r.source_revision != source {
                return Err(bad("Writer source binding changed."));
            }
            let a = assessment(
                m,
                format,
                source,
                &r.profile,
                &r.extension,
                &r.mime,
                r.losses
                    .into_iter()
                    .map(|l| Loss {
                        path: l.path,
                        message: l.message,
                    })
                    .collect(),
                r.warnings,
            )?;
            (a, r.bytes)
        }
        Format::Xdna => {
            let r = xdna_export::export_xdna(m, xdna_export::Policy::AllowOmissions)?;
            if r.source_revision != source {
                return Err(bad("Writer source binding changed."));
            }
            let a = assessment(
                m,
                format,
                source,
                &r.profile,
                &r.extension,
                &r.mime,
                r.losses
                    .into_iter()
                    .map(|l| Loss {
                        path: l.path,
                        message: l.message,
                    })
                    .collect(),
                r.warnings,
            )?;
            (a, r.bytes)
        }
        Format::Gcg { options } => {
            let r = gcg_sequence::export_gcg(m, options)?;
            let losses = gcg_losses(m, options, &r.warnings);
            let a = assessment(
                m,
                format,
                source,
                GCG_PROFILE,
                ".gcg",
                "text/plain;charset=utf-8",
                losses,
                r.warnings,
            )?;
            (a, r.text.into_bytes())
        }
    };
    if !allow_omissions && !assessment.losses.is_empty() {
        return Err(bad(
            "This export omits information; explicitly acknowledge its current loss review before preparing bytes.",
        ));
    }
    let result = Prepared { assessment, bytes };
    result.validate_binding(expected_context_revision)?;
    Ok(result)
}
impl Assessment {
    /// Documentary shape/context checks only; no source or chemistry restoration.
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        if self.engine_version != ENGINE_VERSION
            || !digest(expected)
            || self.context_revision != expected
            || !digest(&self.source_revision)
            || self.filename.is_empty()
            || self.filename.len() > 240
            || self
                .filename
                .chars()
                .any(|c| c.is_control() || "/\\:*?\"<>|".contains(c))
            || self.filename.trim() != self.filename
            || self.filename.starts_with('.')
            || self.filename.ends_with('.')
            || self.losses.len() > MAX_LOSSES
            || self.warnings.len() > MAX_WARNINGS
            || self.losses.iter().any(|l| {
                l.path.is_empty()
                    || l.path.len() > 8192
                    || l.message.is_empty()
                    || l.message.len() > 32768
            })
            || self.warnings.iter().any(|w| w.len() > 32768)
        {
            return Err(bad("Export assessment context or bounds are invalid."));
        }
        let stem = self
            .filename
            .split('.')
            .next()
            .unwrap_or("")
            .to_ascii_uppercase();
        if matches!(stem.as_str(), "CON" | "PRN" | "AUX" | "NUL")
            || (stem.len() == 4
                && (stem.starts_with("COM") || stem.starts_with("LPT"))
                && matches!(stem.as_bytes()[3], b'1'..=b'9'))
        {
            return Err(bad("Output basename is reserved on a supported platform."));
        }
        let known = match self.profile.as_str() {
            snapgene_export::PROFILE => {
                self.mime == "application/octet-stream"
                    && (self.filename.ends_with(".dna") || self.filename.ends_with(".prot"))
            }
            xdna_export::PROFILE => {
                self.mime == "application/octet-stream" && self.filename.ends_with(".xdna")
            }
            GCG_PROFILE => {
                self.mime == "text/plain;charset=utf-8" && self.filename.ends_with(".gcg")
            }
            _ => false,
        };
        if !known {
            return Err(bad("Unknown export profile or inconsistent filename/MIME."));
        }
        Ok(())
    }
}
impl Prepared {
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        self.assessment.validate_binding(expected)?;
        if self.bytes.is_empty() || self.bytes.len() > MAX_BYTES {
            return Err(Error::TooLarge);
        }
        let known = match self.assessment.profile.as_str() {
            snapgene_export::PROFILE => {
                self.bytes.first() == Some(&9)
                    && self.bytes.get(5..13) == Some(b"SnapGene".as_slice())
                    && self.bytes.get(13..15)
                        == Some(if self.assessment.filename.ends_with(".prot") {
                            &[0, 2][..]
                        } else {
                            &[0, 1][..]
                        })
            }
            xdna_export::PROFILE => {
                self.bytes.len() >= 112
                    && self.bytes[0] == 0
                    && matches!(self.bytes[1], 1 | 4)
                    && self.bytes[2] <= 1
                    && (self.bytes[1] != 4 || self.bytes[2] == 0)
                    && self.bytes[111] == 255
            }
            GCG_PROFILE => {
                self.bytes.starts_with(b"!!NA_SEQUENCE 1.0\n")
                    || self.bytes.starts_with(b"!!AA_SEQUENCE 1.0\n")
            }
            _ => false,
        };
        if !known {
            return Err(bad(
                "Prepared bytes do not match the declared export profile.",
            ));
        }
        Ok(())
    }
}
