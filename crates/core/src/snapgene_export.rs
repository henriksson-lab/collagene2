//! Independently written, bounded SnapGene reference exchange profile. No opaque
//! source packet, physical material, experimental status or history is re-emitted.
use crate::{Alphabet, Error, Feature, Location, Molecule, Topology};
use serde::{Deserialize, Serialize};
use std::{collections::BTreeSet, io::Write};
pub const PROFILE: &str = "collagene-snapgene-reference-15-v1";
pub const MAX_BYTES: usize = 32 * 1024 * 1024;
const MAX_XML: usize = 8 * 1024 * 1024;
const MAX_NODES: usize = 100_000;
const MAX_FEATURES: usize = 4096;
const MAX_LOSSES: usize = 256;
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum Policy {
    Strict,
    AllowOmissions,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum LossKind {
    Metadata,
    OriginalEvidence,
    FeatureLocation,
    FeatureReview,
    FeatureColor,
    FeatureDisplay,
    QualifierEncoding,
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
    pub mime: String,
    pub bytes: Vec<u8>,
    pub losses: Vec<Loss>,
    pub warnings: Vec<String>,
}
fn unsupported(s: impl Into<String>) -> Error {
    Error::Unsupported(format!("SnapGene export: {}", s.into()))
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
fn xml_char(c: char) -> bool {
    matches!(c, '\t' | '\n' | '\r')
        || matches!(c as u32,0x20..=0xd7ff|0xe000..=0xfffd|0x10000..=0x10ffff)
}
fn xml_ok(s: &str) -> bool {
    s.chars().all(xml_char)
}
struct Xml {
    data: String,
    nodes: usize,
}
impl Xml {
    fn new() -> Self {
        Self {
            data: String::new(),
            nodes: 0,
        }
    }
    fn push(&mut self, s: &str) -> Result<(), Error> {
        if self
            .data
            .len()
            .checked_add(s.len())
            .is_none_or(|n| n > MAX_XML)
        {
            return Err(Error::TooLarge);
        }
        self.data.push_str(s);
        Ok(())
    }
    fn node(&mut self) -> Result<(), Error> {
        self.nodes += 1;
        if self.nodes > MAX_NODES {
            return Err(Error::TooLarge);
        }
        Ok(())
    }
    fn escaped(&mut self, s: &str) -> Result<(), Error> {
        for c in s.chars() {
            match c {
                '&' => self.push("&amp;")?,
                '<' => self.push("&lt;")?,
                '>' => self.push("&gt;")?,
                '"' => self.push("&quot;")?,
                '\'' => self.push("&apos;")?,
                '\n' => self.push("&#10;")?,
                '\r' => self.push("&#13;")?,
                '\t' => self.push("&#9;")?,
                _ => {
                    if !xml_char(c) {
                        return Err(unsupported("XML-invalid character"));
                    }
                    let mut b = [0; 4];
                    self.push(c.encode_utf8(&mut b))?
                }
            }
        }
        Ok(())
    }
}

fn preflight(m: &Molecule) -> Result<(), Error> {
    if m.name.len() > 1024
        || m.sequence.is_empty()
        || m.sequence.len() > MAX_BYTES
        || m.features.len() > MAX_FEATURES
        || m.metadata.len() > 4096
    {
        return Err(Error::TooLarge);
    }
    let mut text = m.name.len();
    for (k, v) in &m.metadata {
        if k.len() > 4096 || v.len() > MAX_BYTES {
            return Err(Error::TooLarge);
        }
        text = text
            .checked_add(k.len())
            .and_then(|n| n.checked_add(v.len()))
            .filter(|n| *n <= 64 * 1024 * 1024)
            .ok_or(Error::TooLarge)?;
    }
    if let Some(raw) = &m.original_genbank {
        if raw.len() > MAX_BYTES
            || text
                .checked_add(raw.len())
                .is_none_or(|n| n > 64 * 1024 * 1024)
        {
            return Err(Error::TooLarge);
        }
    }
    crate::model::bounds::ModelBudget::default().features(&m.features)?;
    m.validate()?;
    if !xml_ok(&m.name) {
        return Err(unsupported(
            "document name contains XML-invalid control characters",
        ));
    }
    if m.sequence.alphabet() == Alphabet::Rna {
        return Err(unsupported(
            "RNA packet32 is outside this selected DNA/protein writer profile",
        ));
    }
    if m.sequence.alphabet() == Alphabet::Protein && m.topology != Topology::Linear {
        return Err(unsupported("circular proteins are not representable"));
    }
    crate::format_fidelity::reference_export_guard(m)
}
fn source_revision(m: &Molecule) -> Result<String, Error> {
    struct Hash {
        h: blake3::Hasher,
        n: usize,
    }
    impl Write for Hash {
        fn write(&mut self, b: &[u8]) -> std::io::Result<usize> {
            self.n = self
                .n
                .checked_add(b.len())
                .filter(|n| *n <= 128 * 1024 * 1024)
                .ok_or_else(|| std::io::Error::other("bounded SnapGene source binding"))?;
            self.h.update(b);
            Ok(b.len())
        }
        fn flush(&mut self) -> std::io::Result<()> {
            Ok(())
        }
    }
    let mut out = Hash {
        h: blake3::Hasher::new(),
        n: 0,
    };
    serde_json::to_writer(&mut out, m).map_err(|_| Error::TooLarge)?;
    Ok(out.h.finalize().to_hex().to_string())
}
// Collect in literal biological extraction order, including complement-of-join.
// A single global SnapGene direction can represent only one uniform orientation.
fn locations(l: &Location, reverse: bool, out: &mut Vec<(usize, usize, bool)>) -> bool {
    match l {
        Location::Range {
            start,
            end,
            before: false,
            after: false,
        } => {
            out.push((*start, *end, reverse));
            out.len() <= 10_000
        }
        Location::Complement(child) => locations(child, !reverse, out),
        Location::Join(parts) => {
            if reverse {
                parts.iter().rev().all(|p| locations(p, reverse, out))
            } else {
                parts.iter().all(|p| locations(p, reverse, out))
            }
        }
        _ => false,
    }
}
fn color(s: &str) -> bool {
    s.len() == 7 && s.starts_with('#') && s[1..].bytes().all(|b| b.is_ascii_hexdigit())
}
fn feature(
    xml: &mut Xml,
    f: &Feature,
    index: usize,
    protein: bool,
    losses: &mut Vec<Loss>,
) -> Result<(), Error> {
    let path = format!("features[{index}]");
    let mut spans = vec![];
    if !locations(&f.location, false, &mut spans)
        || spans.is_empty()
        || spans.iter().any(|s| s.2 != spans[0].2)
        || protein && spans[0].2
    {
        loss(
            losses,
            format!("{path}.location"),
            LossKind::FeatureLocation,
            "Whole feature omitted: fuzzy, ordered, unresolved, mixed-strand or reverse-protein location is outside the exact uniform-strand profile.",
        )?;
        return Ok(());
    }
    if !xml_ok(&f.kind) || f.kind.trim().is_empty() || !xml_ok(f.label()) {
        loss(
            losses,
            &path,
            LossKind::QualifierEncoding,
            "Whole feature omitted because its kind/name cannot be represented as valid XML.",
        )?;
        return Ok(());
    }
    if f.needs_review {
        loss(
            losses,
            format!("{path}.needs_review"),
            LossKind::FeatureReview,
            "The native annotation-review flag is not a mapped SnapGene field. Retain the native document for its review state.",
        )?;
    }
    let mut direction = if spans[0].2 { "2" } else { "1" };
    let hints = f
        .qualifiers
        .iter()
        .filter(|(k, _)| k == "snapgene_directionality")
        .map(|(_, v)| v.as_deref())
        .collect::<Vec<_>>();
    if !hints.is_empty() {
        if hints.iter().all(|h| *h == hints[0]) {
            match hints[0] {
                Some("0" | "3") if !spans[0].2 => direction = hints[0].unwrap_or("1"),
                Some(value) if value == direction => (),
                _ => loss(
                    losses,
                    format!("{path}.directionality"),
                    LossKind::FeatureDisplay,
                    "Original directionality hint conflicts with the current location or is unsupported. Current location determines the arrow; hint qualifiers remain literal.",
                )?,
            }
        } else {
            loss(
                losses,
                format!("{path}.directionality"),
                LossKind::FeatureDisplay,
                "Conflicting directionality qualifiers remain literal; current location determines the arrow.",
            )?
        }
    }
    let colors = f
        .qualifiers
        .iter()
        .filter(|(k, _)| matches!(k.as_str(), "ApEinfo_fwdcolor" | "ApEinfo_revcolor"))
        .map(|(_, v)| v.as_deref())
        .collect::<Vec<_>>();
    let uniform = colors
        .first()
        .copied()
        .flatten()
        .filter(|c| color(c) && colors.iter().all(|v| *v == Some(*c)));
    if !colors.is_empty() && uniform.is_none() {
        loss(
            losses,
            format!("{path}.color"),
            LossKind::FeatureColor,
            "Distinct, absent or invalid strand colors cannot become one uniform segment color. Original color qualifiers remain literal; no color was normalized or chosen silently.",
        )?;
    }
    xml.node()?;
    xml.push("<Feature name=\"")?;
    xml.escaped(f.label())?;
    xml.push("\" type=\"")?;
    xml.escaped(&f.kind)?;
    xml.push("\" directionality=\"")?;
    xml.push(direction)?;
    xml.push("\">")?;
    if spans[0].2 {
        spans.reverse();
    }
    for (start, end, _) in spans {
        xml.node()?;
        xml.push(&format!(
            "<Segment range=\"{}-{end}\" type=\"standard\"",
            start + 1
        ))?;
        if let Some(c) = uniform {
            xml.push(" color=\"")?;
            xml.escaped(c)?;
            xml.push("\"")?;
        }
        xml.push("/>")?;
    }
    for (q, (key, value)) in f.qualifiers.iter().enumerate() {
        if key.is_empty()
            || !xml_ok(key)
            || value
                .as_ref()
                .is_some_and(|s| !xml_ok(s) || s.contains('<') || s.contains('>'))
        {
            loss(
                losses,
                format!("{path}.qualifiers[{q}]"),
                LossKind::QualifierEncoding,
                "Qualifier omitted: empty/XML-invalid name/value or markup-like rich text is outside the literal-text exchange profile.",
            )?;
            continue;
        }
        xml.node()?;
        xml.push("<Q name=\"")?;
        xml.escaped(key)?;
        if let Some(value) = value {
            xml.push("\"><V text=\"")?;
            xml.node()?;
            xml.escaped(value)?;
            xml.push("\"/></Q>")?;
        } else {
            xml.push("\"/>")?;
        }
    }
    xml.push("</Feature>")?;
    Ok(())
}
struct Prepared {
    assessment: Assessment,
    features: Vec<u8>,
    notes: Vec<u8>,
}
fn prepare(m: &Molecule) -> Result<Prepared, Error> {
    preflight(m)?;
    let mut losses = vec![];
    let mut features = Xml::new();
    features.node()?;
    features.push("<Features>")?;
    for (i, f) in m.features.iter().enumerate() {
        let before = losses.len();
        feature(
            &mut features,
            f,
            i,
            m.sequence.alphabet() == Alphabet::Protein,
            &mut losses,
        )?;
        let modification = m.sequence.alphabet() == Alphabet::Protein
            && [
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
            .any(|kind| f.kind.eq_ignore_ascii_case(kind));
        if modification && losses.len() != before {
            return Err(unsupported(
                "A modified/crosslinked protein annotation would lose required evidence. Retain native protein state; AllowOmissions cannot turn it into an unmodified free-chain reference.",
            ));
        }
    }
    features.push("</Features>")?;
    let mut notes = Xml::new();
    notes.node()?;
    notes.push("<Notes><CustomMapLabel>")?;
    notes.node()?;
    notes.escaped(&m.name)?;
    notes.push("</CustomMapLabel><UseCustomMapLabel>1</UseCustomMapLabel><ConfirmedExperimentally>0</ConfirmedExperimentally>")?;
    notes.node()?;
    notes.node()?;
    let mut emitted = BTreeSet::new();
    for (k, v) in &m.metadata {
        let field = match k.as_str() {
            "definition" | "description" => Some("Description"),
            "comment" => Some("Comments"),
            "accession" => Some("AccessionNumber"),
            _ => None,
        };
        if let Some(field) = field {
            if xml_ok(v) && !v.contains('<') && !v.contains('>') && emitted.insert(field) {
                notes.node()?;
                notes.push(&format!("<{field}>"))?;
                notes.escaped(v)?;
                notes.push(&format!("</{field}>"))?;
            } else {
                loss(
                    &mut losses,
                    format!("metadata.{k}"),
                    LossKind::Metadata,
                    "Metadata omitted because its XML/literal-text value is unsupported (including markup-like text), or another native field already supplied this single SnapGene note; no value was merged silently.",
                )?;
            }
        } else {
            loss(
                &mut losses,
                format!("metadata.{k}"),
                if k.starts_with("source:snapgene:") || k.starts_with("protein:original:") {
                    LossKind::OriginalEvidence
                } else {
                    LossKind::Metadata
                },
                "Native/source metadata is outside the mapped current-document notes. Original packets, source claims, verification, history and caches are not copied into this new file.",
            )?;
        }
    }
    if m.original_genbank.is_some() {
        loss(
            &mut losses,
            "original_genbank",
            LossKind::OriginalEvidence,
            "Original INSDC source text is not embedded or replayed as current SnapGene data; retain the native archive for exact evidence.",
        )?;
    }
    notes.push("</Notes>")?;
    let warnings=vec!["Selected SnapGene reference profile15; vendor-produced file corpus and application acceptance remain unverified. This is not a complete workspace/primer/history conversion.".into(),"Feature joins are flattened only when exact biological order and one uniform orientation are preserved. Display names may be derived from existing labels; literal qualifiers are retained separately.".into(),"Double-stranded DNA is a sequence-reference representation only; no end state, methylation, material ancestry, experimental confirmation or current verification is asserted. Original opaque packets are never emitted.".into()];
    Ok(Prepared {
        assessment: Assessment {
            profile: PROFILE.into(),
            source_revision: source_revision(m)?,
            extension: if m.sequence.alphabet() == Alphabet::Protein {
                ".prot"
            } else {
                ".dna"
            }
            .into(),
            losses,
            warnings,
        },
        features: features.data.into_bytes(),
        notes: notes.data.into_bytes(),
    })
}
/// Bounded loss review, without serializing current sequence into a new packet.
pub fn review(m: &Molecule) -> Result<Assessment, Error> {
    Ok(prepare(m)?.assessment)
}
fn packet(bytes: &mut Vec<u8>, kind: u8, data: &[u8]) -> Result<(), Error> {
    if bytes
        .len()
        .checked_add(5)
        .and_then(|n| n.checked_add(data.len()))
        .is_none_or(|n| n > MAX_BYTES)
    {
        return Err(Error::TooLarge);
    }
    let len: u32 = data.len().try_into().map_err(|_| Error::TooLarge)?;
    bytes.push(kind);
    bytes.extend_from_slice(&len.to_be_bytes());
    bytes.extend_from_slice(data);
    Ok(())
}
pub fn export_snapgene(m: &Molecule, policy: Policy) -> Result<Report, Error> {
    let prepared = prepare(m)?;
    if policy == Policy::Strict && !prepared.assessment.losses.is_empty() {
        return Err(unsupported(format!(
            "{} explicit loss(es) require AllowOmissions after review; use review() for typed details",
            prepared.assessment.losses.len()
        )));
    }
    let protein = m.sequence.alphabet() == Alphabet::Protein;
    let mut cookie = b"SnapGene".to_vec();
    cookie.extend_from_slice(&(if protein { 2u16 } else { 1u16 }).to_be_bytes());
    cookie.extend_from_slice(&15u16.to_be_bytes());
    cookie.extend_from_slice(&15u16.to_be_bytes());
    let mut bytes = vec![];
    packet(&mut bytes, 9, &cookie)?;
    let flags = if protein {
        0
    } else {
        2 | if m.topology == Topology::Circular {
            1
        } else {
            0
        }
    };
    let mut sequence = Vec::with_capacity(m.sequence.len() + 1);
    sequence.push(flags);
    sequence.extend_from_slice(m.sequence.as_str().as_bytes());
    packet(&mut bytes, if protein { 21 } else { 0 }, &sequence)?;
    packet(&mut bytes, 6, &prepared.notes)?;
    if !m.features.is_empty() {
        packet(&mut bytes, 10, &prepared.features)?;
    }
    Ok(Report {
        profile: prepared.assessment.profile,
        source_revision: prepared.assessment.source_revision,
        extension: prepared.assessment.extension,
        mime: "application/octet-stream".into(),
        bytes,
        losses: prepared.assessment.losses,
        warnings: prepared.assessment.warnings,
    })
}
