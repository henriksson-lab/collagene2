//! Regenerated selected XDNA v0 reference profile, independently specified against
//! our bounded v0 reader. No opaque original fields or physical state are replayed.
use collagene_model::{Alphabet, Error, Feature, Location, Molecule, Topology};
use serde::{Deserialize, Serialize};
use std::io::Write;
pub const PROFILE: &str = "collagene-xdna-reference-v0-1";
pub const MAX_BYTES: usize = crate::formats::sequence::MAX_IMPORT_BYTES;
pub const MAX_FEATURES: usize = 255;
const MAX_LOSSES: usize = 1024;
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum Policy {
    Strict,
    AllowOmissions,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum LossKind {
    Name,
    Metadata,
    OriginalEvidence,
    FeatureLocation,
    FeatureReview,
    FeatureText,
    FeatureColor,
    Qualifier,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Loss {
    pub path: String,
    pub kind: LossKind,
    pub message: String,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Assessment {
    pub profile: String,
    pub source_revision: String,
    pub extension: String,
    pub suggested_name: String,
    pub losses: Vec<Loss>,
    pub warnings: Vec<String>,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Report {
    pub profile: String,
    pub source_revision: String,
    pub extension: String,
    pub suggested_name: String,
    pub mime: String,
    pub bytes: Vec<u8>,
    pub losses: Vec<Loss>,
    pub warnings: Vec<String>,
}
fn bad(s: impl Into<String>) -> Error {
    Error::Unsupported(format!("XDNA export: {}", s.into()))
}
fn loss(
    out: &mut Vec<Loss>,
    path: impl Into<String>,
    kind: LossKind,
    message: &str,
) -> Result<(), Error> {
    if out.len() >= MAX_LOSSES {
        return Err(Error::TooLarge);
    }
    out.push(Loss {
        path: path.into(),
        kind,
        message: message.into(),
    });
    Ok(())
}
fn ascii(s: &str) -> bool {
    s.is_ascii()
        && s.bytes()
            .all(|b| !b.is_ascii_control() || b"\t\r\n".contains(&b))
}
// Flatten only the exact biological traversal, with a maximum of two actual
// ranges. An origin wrap is not a license to collapse a discontinuous feature.
fn parts(location: &Location, reverse: bool, out: &mut Vec<(usize, usize, bool)>) -> bool {
    match location {
        Location::Range {
            start,
            end,
            before: false,
            after: false,
        } => {
            out.push((*start, *end, reverse));
            out.len() <= 2
        }
        Location::Complement(child) => parts(child, !reverse, out),
        Location::Join(children) => {
            if reverse {
                children
                    .iter()
                    .rev()
                    .all(|child| parts(child, reverse, out))
            } else {
                children.iter().all(|child| parts(child, reverse, out))
            }
        }
        _ => false,
    }
}
fn coordinates(m: &Molecule, location: &Location) -> Option<(usize, usize, bool)> {
    let mut spans = Vec::new();
    if !parts(location, false, &mut spans) || spans.is_empty() {
        return None;
    }
    let reverse = spans[0].2;
    if spans
        .iter()
        .any(|(a, b, r)| a >= b || *b > m.sequence.len() || *r != reverse)
        || reverse && m.sequence.alphabet() == Alphabet::Protein
    {
        return None;
    }
    if reverse {
        spans.reverse()
    }
    let (start, end) = match spans.as_slice() {
        [(a, b, _)] => (*a, *b),
        [(a, end, _), (zero, b, _)]
            if m.topology == Topology::Circular
                && *end == m.sequence.len()
                && *zero == 0
                && *b <= *a =>
        {
            (*a, *b)
        }
        _ => return None,
    };
    let first = start.checked_add(1)?;
    // Inclusive v0 strand0 stores the high/right endpoint first; wrapping
    // complements preserve the complete complement-of-join biological order.
    Some(if reverse {
        (end, first, true)
    } else {
        (first, end, false)
    })
}
fn chemical_protein_feature(f: &Feature) -> bool {
    [
        "MOD_RES",
        "modified_residue",
        "cross-link",
        "cross_link",
        "bond",
        "LIPID",
        "CARBOHYD",
        "DISULFID",
        "CROSSLNK",
        "NON_STD",
    ]
    .iter()
    .any(|kind| f.kind.eq_ignore_ascii_case(kind))
        || f.qualifiers.iter().any(|(key, _)| {
            ["modification", "modifications", "cross_link", "disulfide"]
                .iter()
                .any(|known| key.eq_ignore_ascii_case(known))
        })
}
struct EncodedFeature {
    label: String,
    note: String,
    kind: String,
    a: usize,
    b: usize,
    reverse: bool,
    extra: String,
}
fn field(f: &Feature, index: usize, key: &str, out: &mut Vec<Loss>) -> Result<String, Error> {
    let mut selected = None;
    for (n, (name, value)) in f.qualifiers.iter().enumerate().filter(|(_, q)| q.0 == key) {
        if selected.is_none() {
            if let Some(value) = value.as_deref().filter(|v| v.len() <= 255 && ascii(v)) {
                selected = Some(value.to_string());
                continue;
            }
        }
        loss(
            out,
            format!("features[{index}].qualifiers[{n}]"),
            LossKind::FeatureText,
            "Qualifier omitted: v0 has one ASCII field of at most255 bytes here. Duplicate, absent, non-ASCII or oversized values are never merged or truncated.",
        )?;
        let _ = name;
    }
    Ok(selected.unwrap_or_default())
}
fn feature(
    m: &Molecule,
    f: &Feature,
    index: usize,
    out: &mut Vec<Loss>,
) -> Result<Option<EncodedFeature>, Error> {
    if m.sequence.alphabet() == Alphabet::Protein && chemical_protein_feature(f) {
        return Err(bad(
            "Modified/crosslinked protein annotations cannot be mapped to an unmodified XDNA protein reference; retain the native archive.",
        ));
    }
    let Some((a, b, reverse)) = coordinates(m, &f.location) else {
        loss(
            out,
            format!("features[{index}].location"),
            LossKind::FeatureLocation,
            "Whole feature omitted: fuzzy, ordered, unresolved, mixed-strand, reverse-protein or discontinuous location is outside the exact range/origin-wrap profile.",
        )?;
        return Ok(None);
    };
    if f.kind.is_empty()
        || f.kind.len() > 255
        || !f
            .kind
            .bytes()
            .all(|b| b.is_ascii_alphanumeric() || b"_:-".contains(&b))
    {
        loss(
            out,
            format!("features[{index}].kind"),
            LossKind::FeatureText,
            "Whole feature omitted: its literal type cannot fit the ASCII v0 key field. No replacement kind was invented.",
        )?;
        return Ok(None);
    }
    if !f.needs_review {
        loss(
            out,
            format!("features[{index}].needs_review"),
            LossKind::FeatureReview,
            "XDNA has no mapped native review-state field; our v0 importer marks every imported annotation for review. The current reviewed state cannot be retained.",
        )?;
    }
    let label = field(f, index, "label", out)?;
    let note = field(f, index, "note", out)?;
    let colors = f
        .qualifiers
        .iter()
        .enumerate()
        .filter(|(_, q)| q.0 == "color")
        .collect::<Vec<_>>();
    let mut extra = String::new();
    if colors.len() == 1 {
        let (n, (_, value)) = colors[0];
        if let Some(color) = value.as_deref().filter(|s| {
            s.len() == 7 && s.starts_with('#') && s[1..].bytes().all(|b| b.is_ascii_hexdigit())
        }) {
            let rgb = [
                u8::from_str_radix(&color[1..3], 16).unwrap(),
                u8::from_str_radix(&color[3..5], 16).unwrap(),
                u8::from_str_radix(&color[5..7], 16).unwrap(),
            ];
            extra = format!("{},{},{}", rgb[0], rgb[1], rgb[2]);
            if color != format!("#{:02x}{:02x}{:02x}", rgb[0], rgb[1], rgb[2]) {
                loss(
                    out,
                    format!("features[{index}].qualifiers[{n}]"),
                    LossKind::FeatureColor,
                    "RGB values are retained numerically, but v0 import reconstructs lowercase hexadecimal text; the original qualifier spelling is not retained.",
                )?;
            }
        } else {
            loss(
                out,
                format!("features[{index}].qualifiers[{n}]"),
                LossKind::FeatureColor,
                "Color omitted: v0 selected profile requires exactly one literal #RRGGBB value; no color was inferred.",
            )?
        }
    } else if !colors.is_empty() {
        loss(
            out,
            format!("features[{index}].color"),
            LossKind::FeatureColor,
            "All duplicate color fields omitted; a single v0 RGB extra cannot preserve their multiplicity or conflicting values.",
        )?
    }
    for (n, (key, _)) in f.qualifiers.iter().enumerate() {
        if !matches!(key.as_str(), "label" | "note" | "color") {
            loss(
                out,
                format!("features[{index}].qualifiers[{n}]"),
                LossKind::Qualifier,
                "Qualifier omitted: arbitrary qualifiers, original coordinates, display flags and opaque extras have no mapped field in this regenerated profile. Current geometry determines strand; display/arrow flags use documented defaults.",
            )?;
        }
    }
    Ok(Some(EncodedFeature {
        label,
        note,
        kind: f.kind.clone(),
        a,
        b,
        reverse,
        extra,
    }))
}
fn source_revision(m: &Molecule) -> Result<String, Error> {
    struct Hash {
        h: blake3::Hasher,
        n: usize,
    }
    impl Write for Hash {
        fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
            self.n = self
                .n
                .checked_add(bytes.len())
                .filter(|n| *n <= 128 * 1024 * 1024)
                .ok_or_else(|| std::io::Error::other("bounded XDNA source"))?;
            self.h.update(bytes);
            Ok(bytes.len())
        }
        fn flush(&mut self) -> std::io::Result<()> {
            Ok(())
        }
    }
    let mut writer = Hash {
        h: blake3::Hasher::new(),
        n: 0,
    };
    serde_json::to_writer(&mut writer, m).map_err(|_| Error::TooLarge)?;
    Ok(writer.h.finalize().to_hex().to_string())
}
struct Prepared {
    assessment: Assessment,
    comment: String,
    features: Vec<EncodedFeature>,
}
fn prepare(m: &Molecule) -> Result<Prepared, Error> {
    collagene_model::bounds::model_preflight(&[m])?;
    m.validate()?;
    crate::formats::fidelity::reference_export_guard(m)?;
    // Known retained XDNA overhang evidence is restrictive even if a producer
    // failed to add the conventional ms-DNA summary marker.
    if ["legacy:xdna_left_length", "legacy:xdna_right_length"]
        .iter()
        .any(|key| {
            m.metadata
                .get(*key)
                .is_some_and(|value| value.parse::<i64>().ok() != Some(0))
        })
    {
        return Err(bad(
            "Nonzero or unresolved original XDNA overhang descriptors cannot be discarded into a clean reference. Retain the native archive.",
        ));
    }
    if m.sequence.alphabet() == Alphabet::Rna {
        return Err(bad(
            "RNA type3 is outside this selected DNA/protein export profile.",
        ));
    }
    if m.sequence.alphabet() == Alphabet::Protein && m.topology != Topology::Linear {
        return Err(bad(
            "Circular protein is not a supported biological XDNA topology.",
        ));
    }
    if m.sequence.is_empty() || m.features.len() > MAX_FEATURES {
        return Err(Error::TooLarge);
    }
    let mut losses = Vec::new();
    loss(
        &mut losses,
        "name",
        LossKind::Name,
        "The selected v0 reader has no mapped embedded molecule-name field. The complete native name is retained only in suggested_name; the caller must choose the filename/context explicitly. It is not hidden in the comment or truncated.",
    )?;
    let mut comment = String::new();
    for (key, value) in &m.metadata {
        if key == "description" && ascii(value) {
            comment = value.clone()
        } else {
            loss(
                &mut losses,
                format!("metadata.{key}"),
                if key.starts_with("legacy:") || key.starts_with("source:") {
                    LossKind::OriginalEvidence
                } else {
                    LossKind::Metadata
                },
                "Metadata omitted: only a literal ASCII description maps to the v0 comment. Original source bytes, nonzero origin, overhang descriptors, verification, history, display state and unrecognized fields are never replayed.",
            )?
        }
    }
    if m.original_genbank.is_some() {
        loss(
            &mut losses,
            "original_genbank",
            LossKind::OriginalEvidence,
            "Original INSDC text is not emitted as current XDNA evidence; retain the native archive for source text and provenance.",
        )?
    }
    let mut features = Vec::new();
    for (index, f) in m.features.iter().enumerate() {
        if let Some(f) = feature(m, f, index, &mut losses)? {
            features.push(f)
        }
    }
    Ok(Prepared{assessment:Assessment{profile:PROFILE.into(),source_revision:source_revision(m)?,extension:".xdna".into(),suggested_name:m.name.clone(),losses,warnings:vec!["Selected independently written v0 reference profile only; compatibility with vendor-produced application files is not established by our own reader roundtrip.".into(),"Origin is explicitly zero; sequence order/topology are current native values. Zero overhang descriptors assert no mapped termini, phosphorylation, material ancestry or experimental verification.".into(),"Feature count is at most255; string fields at most255 ASCII bytes. Strand follows exact current coordinates, visibility=1, reserved=0 and arrow=1. All imported annotations require review.".into()]},comment,features})
}
pub fn review(m: &Molecule) -> Result<Assessment, Error> {
    Ok(prepare(m)?.assessment)
}
fn append(bytes: &mut Vec<u8>, data: &[u8]) -> Result<(), Error> {
    if bytes
        .len()
        .checked_add(data.len())
        .is_none_or(|n| n > MAX_BYTES)
    {
        return Err(Error::TooLarge);
    }
    bytes.extend_from_slice(data);
    Ok(())
}
fn string(bytes: &mut Vec<u8>, value: &str) -> Result<(), Error> {
    let len: u8 = value.len().try_into().map_err(|_| Error::TooLarge)?;
    append(bytes, &[len])?;
    append(bytes, value.as_bytes())
}
pub fn export_xdna(m: &Molecule, policy: Policy) -> Result<Report, Error> {
    let prepared = prepare(m)?;
    if policy == Policy::Strict && !prepared.assessment.losses.is_empty() {
        return Err(bad(format!(
            "{} explicit loss(es) require AllowOmissions after review(); retain a native archive for complete fidelity.",
            prepared.assessment.losses.len()
        )));
    }
    let mut bytes = vec![0u8; 112];
    bytes[1] = if m.sequence.alphabet() == Alphabet::Protein {
        4
    } else {
        1
    };
    bytes[2] = u8::from(m.topology == Topology::Circular);
    bytes[28..32].copy_from_slice(
        &u32::try_from(m.sequence.len())
            .map_err(|_| Error::TooLarge)?
            .to_be_bytes(),
    );
    bytes[96..100].copy_from_slice(
        &u32::try_from(prepared.comment.len())
            .map_err(|_| Error::TooLarge)?
            .to_be_bytes(),
    );
    bytes[111] = 255;
    append(&mut bytes, m.sequence.as_str().as_bytes())?;
    append(&mut bytes, prepared.comment.as_bytes())?;
    append(&mut bytes, &[0])?;
    string(&mut bytes, "0")?;
    string(&mut bytes, "0")?;
    append(
        &mut bytes,
        &[u8::try_from(prepared.features.len()).map_err(|_| Error::TooLarge)?],
    )?;
    for f in prepared.features {
        string(&mut bytes, &f.label)?;
        string(&mut bytes, &f.note)?;
        string(&mut bytes, &f.kind)?;
        string(&mut bytes, &f.a.to_string())?;
        string(&mut bytes, &f.b.to_string())?;
        append(&mut bytes, &[u8::from(!f.reverse), 1, 0, 1])?;
        string(&mut bytes, &f.extra)?
    }
    let a = prepared.assessment;
    Ok(Report {
        profile: a.profile,
        source_revision: a.source_revision,
        extension: a.extension,
        suggested_name: a.suggested_name,
        mime: "application/octet-stream".into(),
        bytes,
        losses: a.losses,
        warnings: a.warnings,
    })
}
