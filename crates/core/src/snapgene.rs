//! Independently implemented, deliberately limited SnapGene DNA/protein packet/XML import.
//! See docs/SNAPGENE_FORMAT.md for layout provenance, unsupported semantics and limits.
use crate::formats::{ImportReport, MAX_IMPORT_BYTES};
use base64::{Engine as _, engine::general_purpose::STANDARD};
use collagene_model::{Alphabet, Error, Feature, Location, Molecule, Sequence, Topology};
use roxmltree::{Document, Node, ParsingOptions};
use std::collections::{BTreeMap, BTreeSet};

pub const ORIGINAL_BASE64_KEY: &str = "source:snapgene:original_base64";
const ORIGINAL_HASH_KEY: &str = "source:snapgene:original_blake3";
const MAX_PACKETS: usize = 4096;
const MAX_XML_BYTES: usize = 8 * 1024 * 1024;
const MAX_XML_NODES: u32 = 100_000;
const MAX_FEATURES: usize = 10_000;
const MAX_WARNINGS: usize = 100;
fn malformed(message: impl Into<String>) -> Error {
    Error::Format(message.into())
}
fn unsupported(message: impl Into<String>) -> Error {
    Error::Unsupported(message.into())
}
fn warn(warnings: &mut BTreeSet<String>, message: impl Into<String>) {
    let message = message.into();
    if warnings.contains(&message) {
        return;
    }
    if warnings.len() < MAX_WARNINGS - 1 {
        warnings.insert(message);
    } else {
        warnings.insert("Additional SnapGene fidelity warnings were truncated; packet_counts and the retained original source describe the remaining data.".into());
    }
}

fn xml<'a>(bytes: &'a [u8], expected: &str) -> Result<Document<'a>, Error> {
    if bytes.len() > MAX_XML_BYTES {
        return Err(Error::TooLarge);
    }
    let text = std::str::from_utf8(bytes).map_err(|_| malformed("SnapGene XML is not UTF-8"))?;
    if text.contains("<!DOCTYPE") || text.contains("<!ENTITY") {
        return Err(unsupported(
            "SnapGene XML DTD/entity declarations are not supported",
        ));
    }
    let document = Document::parse_with_options(
        text,
        ParsingOptions {
            allow_dtd: false,
            nodes_limit: MAX_XML_NODES,
            entity_resolver: None,
        },
    )
    .map_err(|e| malformed(format!("invalid SnapGene XML: {e}")))?;
    if document.root_element().tag_name().name() != expected {
        return Err(malformed(format!("expected SnapGene XML root {expected}")));
    }
    for node in document.descendants() {
        if node.ancestors().take(34).count() > 33 {
            return Err(malformed("SnapGene XML nesting exceeds 32 element levels"));
        }
        if node.is_element() && node.tag_name().namespace().is_some() {
            return Err(unsupported(
                "namespaced SnapGene XML is not mapped by this importer",
            ));
        }
    }
    Ok(document)
}
fn text_value(node: Node<'_, '_>) -> String {
    node.descendants()
        .filter(|n| n.is_text())
        .filter_map(|n| n.text())
        .collect()
}
fn range(raw: &str, len: usize, topology: Topology) -> Result<Location, Error> {
    let (start, end) = raw
        .split_once('-')
        .ok_or_else(|| malformed("invalid SnapGene feature range"))?;
    if start.is_empty()
        || end.is_empty()
        || !start.bytes().chain(end.bytes()).all(|b| b.is_ascii_digit())
    {
        return Err(malformed(
            "SnapGene feature range needs two positive integer coordinates",
        ));
    }
    let start = start.parse::<usize>().map_err(|_| Error::TooLarge)?;
    let end = end.parse::<usize>().map_err(|_| Error::TooLarge)?;
    if start == 0 || end == 0 || start > len || end > len {
        return Err(Error::Coordinates(format!(
            "SnapGene feature range {raw} outside {len} bases"
        )));
    }
    if start <= end {
        Ok(Location::range(start - 1, end))
    } else if topology == Topology::Circular {
        Ok(Location::Join(vec![
            Location::range(start - 1, len),
            Location::range(0, end),
        ]))
    } else {
        Err(Error::Coordinates(
            "origin-spanning SnapGene feature on a linear sequence".into(),
        ))
    }
}
fn color_valid(value: &str) -> bool {
    value.len() == 7 && value.starts_with('#') && value[1..].bytes().all(|b| b.is_ascii_hexdigit())
}
fn parse_features(
    bytes: &[u8],
    molecule: &mut Molecule,
    warnings: &mut BTreeSet<String>,
) -> Result<(), Error> {
    let doc = xml(bytes, "Features")?;
    let root = doc.root_element();
    for child in root.children().filter(|n| n.is_element()) {
        if child.tag_name().name() != "Feature" {
            warn(
                warnings,
                "Unknown SnapGene feature-list elements were retained only in the original source.",
            );
            continue;
        }
        if molecule.features.len() >= MAX_FEATURES {
            return Err(Error::TooLarge);
        }
        let kind = child
            .attribute("type")
            .filter(|s| !s.trim().is_empty())
            .unwrap_or("misc_feature");
        let name = child.attribute("name");
        let direction = child.attribute("directionality").unwrap_or("0");
        if !matches!(direction, "0" | "1" | "2" | "3") {
            return Err(unsupported(format!(
                "unknown SnapGene feature directionality {direction}"
            )));
        }
        if molecule.sequence.alphabet() == Alphabet::Protein && direction == "2" {
            return Err(unsupported(
                "reverse protein feature direction is outside this linear amino-acid import profile",
            ));
        }
        let mut review = direction == "3";
        let mut qualifiers = vec![("snapgene_directionality".into(), Some(direction.into()))];
        if let Some(name) = name {
            qualifiers.push(("snapgene_name".into(), Some(name.into())));
        }
        if direction == "0" || direction == "3" {
            warn(
                warnings,
                "SnapGene non-directional/bidirectional display semantics are retained as qualifiers; the current location model cannot render both arrow styles exactly.",
            );
        }
        let mut locations = vec![];
        let mut segment_colors = BTreeSet::new();
        let mut gap = false;
        for attribute in child.attributes() {
            if !matches!(attribute.name(), "name" | "type" | "directionality") {
                review = true;
                warn(
                    warnings,
                    "Some SnapGene feature attributes have no mapped biological/display semantics; review flagged annotations and retain the original source.",
                );
            }
        }
        for (segment_count, segment) in child
            .children()
            .filter(|n| n.has_tag_name("Segment"))
            .enumerate()
        {
            if segment_count >= 10_000 {
                return Err(Error::TooLarge);
            }
            let location = range(
                segment
                    .attribute("range")
                    .ok_or_else(|| malformed("SnapGene feature segment lacks range"))?,
                molecule.sequence.len(),
                molecule.topology,
            )?;
            match segment.attribute("type").unwrap_or("standard") {
                "standard" => locations.push(location),
                "gap" => {
                    gap = true;
                }
                _ => {
                    gap = true;
                    warn(
                        warnings,
                        "Unknown SnapGene segment types were retained as unresolved locations; no joined sequence is asserted.",
                    );
                }
            }
            if let Some(color) = segment.attribute("color") {
                if color_valid(color) {
                    segment_colors.insert(color.to_string());
                } else {
                    review = true;
                    warn(
                        warnings,
                        "A SnapGene segment has an unsupported color value; source bytes retain it.",
                    );
                }
            }
            if segment
                .attributes()
                .any(|a| !matches!(a.name(), "range" | "type" | "color"))
            {
                review = true;
                warn(
                    warnings,
                    "Named/styled SnapGene subsegments are retained in original source only; affected features require review.",
                );
            }
        }
        for qualifier in child.children().filter(|n| n.has_tag_name("Q")) {
            let qname = qualifier
                .attribute("name")
                .filter(|s| !s.is_empty())
                .ok_or_else(|| malformed("SnapGene qualifier lacks a name"))?;
            let values: Vec<_> = qualifier.children().filter(|n| n.is_element()).collect();
            if values.is_empty() {
                qualifiers.push((qname.into(), None));
            }
            for value in values {
                if !value.has_tag_name("V") {
                    return Err(unsupported("unmapped SnapGene qualifier element"));
                }
                let mut found = None;
                for attribute in ["text", "int", "predef"] {
                    if let Some(text) = value.attribute(attribute) {
                        if found.is_some() {
                            return Err(malformed(
                                "SnapGene qualifier has conflicting value attributes",
                            ));
                        }
                        found = Some(text.to_string());
                    }
                }
                let value = found
                    .ok_or_else(|| unsupported("unmapped SnapGene qualifier value encoding"))?;
                if value.contains('<') || value.contains('>') {
                    warn(
                        warnings,
                        "SnapGene rich-text qualifiers are retained as literal text/markup; rendered HTML formatting is not reproduced.",
                    );
                }
                qualifiers.push((qname.into(), Some(value)));
            }
        }
        if child
            .children()
            .any(|n| n.is_element() && !matches!(n.tag_name().name(), "Segment" | "Q"))
        {
            review = true;
            warn(
                warnings,
                "Unmapped SnapGene feature children are retained in the original source; affected annotations require review.",
            );
        }
        if !qualifiers.iter().any(|(k, _)| k == "label") {
            if let Some(name) = name {
                qualifiers.insert(0, ("label".into(), Some(name.into())));
            }
        }
        if segment_colors.len() == 1 {
            let color = segment_colors.into_iter().next().unwrap();
            qualifiers.push(("ApEinfo_fwdcolor".into(), Some(color.clone())));
            qualifiers.push(("ApEinfo_revcolor".into(), Some(color)));
        } else if segment_colors.len() > 1 {
            review = true;
            warn(
                warnings,
                "Per-segment SnapGene colors are not representable by one feature color; source bytes retain them.",
            );
        }
        let mut location = if gap {
            review = true;
            warn(
                warnings,
                "SnapGene gap/nonstandard segment features are retained as unresolved locations; no joined sequence is asserted.",
            );
            Location::Unresolved(format!("snapgene-original-feature:{}", child.range().start))
        } else {
            if locations.is_empty() {
                return Err(malformed("SnapGene feature lacks a standard segment"));
            }
            if locations.len() == 1 {
                locations.remove(0)
            } else {
                Location::Join(locations)
            }
        };
        if direction == "2" {
            location = Location::Complement(Box::new(location));
        }
        if molecule.sequence.alphabet() == Alphabet::Protein
            && ["DISULFID", "cross-link", "cross_link", "bond"]
                .iter()
                .any(|k| kind.eq_ignore_ascii_case(k))
        {
            review = true;
            location = Location::Unresolved(format!(
                "snapgene-original-protein-crosslink:{}",
                child.range().start
            ));
            warn(
                warnings,
                "Protein crosslink/bond segment semantics are not established by this profile. The feature stays unresolved/reviewed with its original source; no contiguous peptide or physical bond is inferred.",
            );
        }
        molecule.features.push(Feature {
            kind: kind.into(),
            location,
            qualifiers,
            needs_review: review,
        });
    }
    Ok(())
}
fn parse_notes(
    bytes: &[u8],
    molecule: &mut Molecule,
    warnings: &mut BTreeSet<String>,
) -> Result<(), Error> {
    let doc = xml(bytes, "Notes")?;
    let mut known = BTreeMap::new();
    for node in doc.root_element().children().filter(|n| n.is_element()) {
        let key = node.tag_name().name();
        let value = text_value(node);
        if known.insert(key.to_string(), value.clone()).is_some() {
            warn(
                warnings,
                "Duplicate SnapGene note fields are retained only in original source; the mapped summary uses the last value.",
            );
        }
        if matches!(
            key,
            "Type"
                | "Created"
                | "LastModified"
                | "AccessionNumber"
                | "Description"
                | "Comments"
                | "UUID"
                | "CustomMapLabel"
                | "UseCustomMapLabel"
                | "SequenceClass"
                | "ConfirmedExperimentally"
        ) {
            // These are source assertions only; none grants verified status or tool authority.
            molecule
                .metadata
                .insert(format!("source:snapgene:note:{key}"), value);
        } else {
            warn(
                warnings,
                "Some SnapGene notes are retained in original bytes only.",
            );
        }
    }
    if known.get("UseCustomMapLabel").is_some_and(|s| s == "1") {
        if let Some(label) = known.get("CustomMapLabel").filter(|s| !s.trim().is_empty()) {
            molecule.name = label.clone();
        }
    } else if let Some(accession) = known
        .get("AccessionNumber")
        .filter(|s| !s.trim().is_empty())
    {
        molecule.name = accession.clone();
    }
    Ok(())
}

pub fn import_snapgene(bytes: &[u8]) -> Result<ImportReport, Error> {
    if bytes.len() > MAX_IMPORT_BYTES {
        return Err(Error::TooLarge);
    }
    if bytes.len() < 19 || bytes[..5] != [9, 0, 0, 0, 14] || &bytes[5..13] != b"SnapGene" {
        return Err(malformed("not a supported SnapGene cookie"));
    }
    let sequence_type = u16::from_be_bytes([bytes[13], bytes[14]]);
    if !matches!(sequence_type, 1 | 2) {
        return Err(unsupported(format!(
            "SnapGene sequence type {sequence_type}; this profile supports uncompressed DNA and linear protein only"
        )));
    }
    let protein = sequence_type == 2;
    let export_version = u16::from_be_bytes([bytes[15], bytes[16]]);
    let import_version = u16::from_be_bytes([bytes[17], bytes[18]]);
    if export_version > 15 || import_version > 15 {
        return Err(unsupported(
            "SnapGene producer/import revisions above 15 have not been audited by this import profile",
        ));
    }
    let mut packets = vec![];
    let mut offset = 19;
    let mut seen = BTreeSet::new();
    let mut packet_counts = BTreeMap::<u8, usize>::new();
    while offset < bytes.len() {
        if packets.len() >= MAX_PACKETS {
            return Err(Error::TooLarge);
        }
        let header = bytes
            .get(offset..offset + 5)
            .ok_or_else(|| malformed("truncated SnapGene packet header"))?;
        let kind = header[0];
        let len = u32::from_be_bytes(header[1..5].try_into().unwrap()) as usize;
        let start = offset + 5;
        let end = start.checked_add(len).ok_or(Error::TooLarge)?;
        let payload = bytes
            .get(start..end)
            .ok_or_else(|| malformed("SnapGene packet length exceeds remaining file"))?;
        if matches!(kind, 0 | 1 | 6 | 9 | 10 | 21 | 32) && !seen.insert(kind) {
            return Err(malformed(format!(
                "duplicate singleton SnapGene packet {kind}"
            )));
        }
        if kind == 9 {
            return Err(malformed(
                "SnapGene cookie must occur only at the beginning",
            ));
        }
        packets.push((kind, payload));
        *packet_counts.entry(kind).or_default() += 1;
        offset = end;
    }
    if seen.contains(&1) || seen.contains(&32) {
        return Err(unsupported(
            "compressed sequence or RNA packets are not supported by this SnapGene import profile",
        ));
    }
    if protein && seen.contains(&0) || !protein && seen.contains(&21) {
        return Err(malformed(
            "SnapGene cookie alphabet conflicts with its sequence packet, or mixed DNA/protein packets are present",
        ));
    }
    let sequence_tag = if protein { 21 } else { 0 };
    let data = packets
        .iter()
        .find(|(kind, _)| *kind == sequence_tag)
        .map(|(_, data)| *data)
        .ok_or_else(|| {
            malformed(if protein {
                "SnapGene protein file has no uncompressed protein packet21"
            } else {
                "SnapGene file has no DNA packet"
            })
        })?;
    let (&flags, bases) = data
        .split_first()
        .ok_or_else(|| malformed("empty SnapGene sequence packet"))?;
    if protein {
        if flags & 1 != 0 {
            return Err(unsupported(
                "circular SnapGene proteins are not representable",
            ));
        }
        if flags != 0 {
            return Err(unsupported(
                "nonzero/unknown protein property flags are outside the selected import profile",
            ));
        }
    } else {
        if flags & 2 == 0 {
            return Err(unsupported(
                "single-stranded SnapGene DNA is not representable by this molecule import profile",
            ));
        }
        if flags & !0x1f != 0 {
            return Err(unsupported("unknown SnapGene DNA property bits"));
        }
    }
    if bases.iter().any(|b| b.is_ascii_whitespace()) {
        return Err(malformed("SnapGene sequence packet contains whitespace"));
    }
    let bases = std::str::from_utf8(bases)
        .map_err(|_| malformed("SnapGene sequence is not ASCII/UTF-8"))?;
    let mut molecule = Molecule {
        name: if protein {
            "SnapGene protein"
        } else {
            "SnapGene sequence"
        }
        .into(),
        sequence: Sequence::new(
            if protein {
                Alphabet::Protein
            } else {
                Alphabet::Dna
            },
            bases,
        )?,
        topology: if !protein && flags & 1 != 0 {
            Topology::Circular
        } else {
            Topology::Linear
        },
        features: vec![],
        metadata: BTreeMap::new(),
        original_genbank: None,
    };
    let mut warnings = BTreeSet::new();
    warn(
        &mut warnings,
        "Experimental SnapGene import profile: sequence, topology and supported annotations are mapped; this is not complete SnapGene workspace/history fidelity. Keep the original source.",
    );
    if bases.bytes().any(|b| b.is_ascii_lowercase()) {
        warn(
            &mut warnings,
            "SnapGene sequence case was normalized; exact original case remains in retained source bytes.",
        );
    }
    if protein {
        molecule.metadata.insert("protein:authority".into(),"Imported SnapGene amino acids and annotations are documentary reference data, not modeled PTM/disulfide/terminal chemistry, physical material or experimental verification.".into());
        warn(
            &mut warnings,
            "Protein sequence and amino-acid coordinates are reference data only. No PTM, disulfide, terminal processing, folding or material state is inferred from source annotations.",
        );
    }
    if flags & 0x1c != 0 || packet_counts.contains_key(&8) {
        molecule
            .metadata
            .insert("source:snapgene:chemistry_unmapped".into(), "true".into());
        if protein {
            molecule.metadata.insert("protein:modifications".into(),"Unmapped original SnapGene additional sequence properties; retain original bytes. No unmodified/free-chain chemical state is asserted.".into());
        }
        warn(
            &mut warnings,
            if protein {
                "SnapGene additional protein properties remain unmapped source evidence, not PTM/terminal or unmodified free-chain chemistry."
            } else {
                "SnapGene methylation/end/overhang properties are not converted to verified chemistry. Cloning checks must not assume those source assertions are modeled."
            },
        );
    }
    for (kind, payload) in &packets {
        match kind {
            0 | 21 => {}
            6 => parse_notes(payload, &mut molecule, &mut warnings)?,
            10 => parse_features(payload, &mut molecule, &mut warnings)?,
            _ => {}
        }
    }
    for (kind, count) in &packet_counts {
        if !matches!(kind, 0 | 6 | 10 | 21) {
            let meaning = match kind {
                5 => "primers and binding predictions",
                7 | 11 | 29 | 30 => "cloning history",
                8 => "end/overhang properties",
                16 | 18 | 27 => "traces/trace alignments",
                23 => "attachments",
                _ => "unmapped data",
            };
            warn(
                &mut warnings,
                format!(
                    "SnapGene packet {kind} ({count} packet(s), {meaning}) retained in original bytes only; not converted or scientifically verified."
                ),
            );
        }
    }
    molecule.metadata.insert(
        "source_format".into(),
        if protein {
            "SnapGene protein (experimental import)"
        } else {
            "SnapGene DNA (experimental import)"
        }
        .into(),
    );
    molecule.metadata.insert(
        "source:snapgene:export_version".into(),
        export_version.to_string(),
    );
    molecule.metadata.insert(
        "source:snapgene:import_version".into(),
        import_version.to_string(),
    );
    molecule.metadata.insert(
        if protein {
            "source:snapgene:protein_flags"
        } else {
            "source:snapgene:dna_flags"
        }
        .into(),
        flags.to_string(),
    );
    molecule.metadata.insert(
        "source:snapgene:packet_counts".into(),
        serde_json::to_string(&packet_counts).map_err(|e| malformed(e.to_string()))?,
    );
    molecule
        .metadata
        .insert(ORIGINAL_BASE64_KEY.into(), STANDARD.encode(bytes));
    molecule.metadata.insert(
        ORIGINAL_HASH_KEY.into(),
        blake3::hash(bytes).to_hex().to_string(),
    );
    molecule.validate()?;
    // A successful import must remain saveable in the bounded .cg2 single-document format.
    if serde_json::to_vec(&molecule)
        .map_err(|e| malformed(e.to_string()))?
        .len()
        > MAX_IMPORT_BYTES
    {
        return Err(unsupported(
            "SnapGene source retention and mapped annotations exceed the 32 MiB document budget; retain the source and use a smaller/intermediate export",
        ));
    }
    Ok(ImportReport {
        molecules: vec![molecule],
        warnings: warnings.into_iter().collect(),
    })
}

/// Recover the unmodified attachment, NOT a SnapGene serialization of current edits.
pub fn original_snapgene_bytes(molecule: &Molecule) -> Result<Vec<u8>, Error> {
    let encoded = molecule
        .metadata
        .get(ORIGINAL_BASE64_KEY)
        .ok_or_else(|| malformed("no original SnapGene source attachment"))?;
    if encoded.len() > MAX_IMPORT_BYTES.div_ceil(3).saturating_mul(4) {
        return Err(Error::TooLarge);
    }
    let bytes = STANDARD
        .decode(encoded)
        .map_err(|_| malformed("invalid original SnapGene base64 attachment"))?;
    if bytes.len() > MAX_IMPORT_BYTES {
        return Err(Error::TooLarge);
    }
    let expected = molecule
        .metadata
        .get(ORIGINAL_HASH_KEY)
        .ok_or_else(|| malformed("original SnapGene attachment lacks integrity hash"))?;
    if blake3::hash(&bytes).to_hex().as_str() != expected {
        return Err(malformed("original SnapGene attachment checksum mismatch"));
    }
    Ok(bytes)
}
