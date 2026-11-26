//! Explicit lower-strand-only reference projection for original Collagene XML.
//! Ordinary migration is delegated unchanged to legacy_formats. No complementary
//! strand, termini, history or material authority is synthesized by this adapter.
use crate::{
    Error, Location, Molecule, Topology,
    formats::{ImportReport, MAX_IMPORT_BYTES},
    legacy_formats as legacy,
};
use base64::{Engine as _, engine::general_purpose::STANDARD};
use roxmltree::{Document, Node, ParsingOptions};
use serde::{Deserialize, Serialize};
use std::io::Write;
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum ReferencePolicy {
    UpperOnly,
    AllowLowerOnly,
}
fn bad(message: &str) -> Error {
    Error::Format(message.into())
}
fn unsupported(message: &str) -> Error {
    Error::Unsupported(message.into())
}
fn document(bytes: &[u8]) -> Result<Document<'_>, Error> {
    if bytes.len() > MAX_IMPORT_BYTES {
        return Err(Error::TooLarge);
    }
    let text = std::str::from_utf8(bytes).map_err(|_| bad("Collagene XML must be UTF-8"))?;
    let doc = Document::parse_with_options(
        text,
        ParsingOptions {
            allow_dtd: false,
            nodes_limit: 300_000,
            ..Default::default()
        },
    )
    .map_err(|e| Error::Format(format!("invalid Collagene XML: {e}")))?;
    let mut pending = vec![(doc.root_element(), 1usize)];
    while let Some((node, depth)) = pending.pop() {
        if depth > 64 {
            return Err(unsupported("Collagene XML nesting exceeds64 levels"));
        }
        if node.tag_name().namespace().is_some() {
            return Err(unsupported(
                "namespaced Collagene elements require another schema",
            ));
        }
        if node.tag_name().name().len() > 128 {
            return Err(Error::TooLarge);
        }
        pending.extend(
            node.children()
                .filter(Node::is_element)
                .map(|child| (child, depth + 1)),
        );
    }
    Ok(doc)
}
fn records<'a>(doc: &'a Document<'a>) -> Result<Vec<Node<'a, 'a>>, Error> {
    let root = doc.root_element();
    let result = match root.tag_name().name() {
        "madgene" => {
            let children = root.children().filter(Node::is_element).collect::<Vec<_>>();
            if children.iter().any(|n| !n.has_tag_name("oneseq")) {
                return Err(unsupported(
                    "unknown Collagene project record; no partial import",
                ));
            }
            children
        }
        "f" | "oneseq" => vec![root],
        "sequence" if root.children().any(|n| n.has_tag_name("sequence")) => vec![root],
        _ => return Err(unsupported("unrecognized Collagene XML container")),
    };
    if result.is_empty() || result.len() > legacy::MAX_RECORDS {
        return Err(bad("Collagene project requires1–1000 records"));
    }
    Ok(result)
}
fn rows<'a>(record: Node<'a, 'a>) -> Result<(Node<'a, 'a>, &'a str, &'a str), Error> {
    let mut children = record.children().filter(|n| n.has_tag_name("sequence"));
    let sequence = children
        .next()
        .ok_or_else(|| bad("Collagene record lacks sequence child"))?;
    if children.next().is_some() {
        return Err(bad("Collagene record has duplicate sequence children"));
    }
    let upper = sequence
        .attribute("su")
        .ok_or_else(|| bad("Collagene sequence lacks upper row"))?;
    let lower = sequence
        .attribute("sl")
        .ok_or_else(|| bad("Collagene sequence lacks lower row"))?;
    if upper.is_empty()
        || upper.len() != lower.len()
        || upper.len() > crate::model::MAX_SEQUENCE_LEN
    {
        return Err(bad(
            "Collagene aligned row lengths differ, are empty or exceed bounds",
        ));
    }
    if !upper.is_ascii()
        || !lower.is_ascii()
        || upper
            .bytes()
            .chain(lower.bytes())
            .any(|b| b != b' ' && !b"ACGTRYSWKMBDHVN".contains(&b.to_ascii_uppercase()))
    {
        return Err(bad(
            "Collagene strand rows require literal IUPAC DNA or space cells",
        ));
    }
    Ok((sequence, upper, lower))
}
/// A review hint only. Malformed containers fail explicitly; no source is
/// imported and no ordinary record is dropped while inspecting this condition.
pub fn needs_lower_only_review(bytes: &[u8]) -> Result<bool, Error> {
    let doc = document(bytes)?;
    let records = records(&doc)?;
    let mut needs = false;
    for record in records {
        let (_, upper, lower) = rows(record)?;
        needs |= upper.bytes().all(|b| b == b' ') && lower.bytes().any(|b| b != b' ');
    }
    Ok(needs)
}
fn value_range(
    bytes: &[u8],
    attribute: roxmltree::Attribute<'_, '_>,
) -> Result<std::ops::Range<usize>, Error> {
    // roxmltree's convenience range_value compresses whitespace-around-equals
    // length to u8. Locate the actual quote so unusually spaced XML cannot
    // redirect an edit into another attribute.
    let range = attribute.range();
    let raw = bytes
        .get(range.clone())
        .ok_or_else(|| bad("invalid XML attribute range"))?;
    let equal = raw
        .iter()
        .position(|b| *b == b'=')
        .ok_or_else(|| bad("XML attribute lacks equals"))?;
    let quote = raw[equal + 1..]
        .iter()
        .position(|b| !matches!(*b, b' ' | b'\t' | b'\r' | b'\n'))
        .map(|n| equal + 1 + n)
        .ok_or_else(|| bad("XML attribute lacks quoted value"))?;
    if !matches!(raw[quote], b'\'' | b'"')
        || raw.last() != Some(&raw[quote])
        || quote + 1 >= raw.len()
    {
        return Err(bad("invalid XML attribute quoting"));
    }
    Ok(range.start + quote + 1..range.end - 1)
}
fn temporary_lower_record(
    bytes: &[u8],
    record: Node<'_, '_>,
    sequence: Node<'_, '_>,
    lower: &str,
) -> Result<Vec<u8>, Error> {
    let range = record.range();
    let mut edits = vec![];
    for attribute in sequence.attributes() {
        match attribute.name() {
            "su" => edits.push((
                value_range(bytes, attribute)?,
                lower
                    .bytes()
                    .rev()
                    .map(|b| b.to_ascii_uppercase() as char)
                    .collect::<String>(),
            )),
            "sl" => edits.push((value_range(bytes, attribute)?, " ".repeat(lower.len()))),
            _ => (),
        }
    }
    edits.sort_by_key(|(range, _)| range.start);
    let mut output = Vec::new();
    let mut offset = range.start;
    for (span, replacement) in edits {
        output.extend_from_slice(
            bytes
                .get(offset..span.start)
                .ok_or_else(|| bad("invalid XML attribute span"))?,
        );
        output.extend_from_slice(replacement.as_bytes());
        offset = span.end;
        if output.len() > MAX_IMPORT_BYTES {
            return Err(Error::TooLarge);
        }
    }
    output.extend_from_slice(
        bytes
            .get(offset..range.end)
            .ok_or_else(|| bad("invalid XML record span"))?,
    );
    if output.len() > MAX_IMPORT_BYTES {
        return Err(Error::TooLarge);
    }
    Ok(output)
}
fn source_interval(feature: &crate::Feature) -> Result<(i64, i64, &str), Error> {
    let range = feature
        .qualifiers
        .iter()
        .find_map(|(k, v)| {
            (k == "legacy:aligned_range")
                .then_some(v.as_deref())
                .flatten()
        })
        .ok_or_else(|| bad("legacy parser omitted source interval"))?;
    let body = range
        .strip_prefix('[')
        .and_then(|s| s.strip_suffix(')'))
        .ok_or_else(|| bad("invalid validated source interval"))?;
    let (a, b) = body
        .split_once(',')
        .ok_or_else(|| bad("invalid source interval pair"))?;
    let direction = feature
        .qualifiers
        .iter()
        .find_map(|(k, v)| {
            (k == "legacy:orientation")
                .then_some(v.as_deref())
                .flatten()
        })
        .ok_or_else(|| bad("legacy parser omitted source orientation"))?;
    Ok((
        a.parse().map_err(|_| bad("invalid source start"))?,
        b.parse().map_err(|_| bad("invalid source end"))?,
        direction,
    ))
}
fn lower_location(
    from: i64,
    to: i64,
    direction: &str,
    prefix: &[usize],
    topology: Topology,
) -> Option<Location> {
    let n = prefix.len() - 1;
    if from < 0 || to < 0 || from as u64 > n as u64 || to as u64 > n as u64 || from == to {
        return None;
    }
    let (a, b) = (from as usize, to as usize);
    let spans = if a < b {
        vec![(a, b)]
    } else if topology == Topology::Circular {
        vec![(a, n), (0, b)]
    } else {
        return None;
    };
    let length = *prefix.last()?;
    let mut mapped = Vec::new();
    for (a, b) in spans.into_iter().rev() {
        if a == b {
            continue;
        }
        // Projection can omit source-gap columns outside a feature, but never
        // fills or joins missing cells inside its claimed readable footprint.
        if prefix[b] - prefix[a] != b - a {
            return None;
        }
        let (start, end) = (length - prefix[b], length - prefix[a]);
        if start >= end {
            return None;
        }
        mapped.push(Location::range(start, end));
    }
    let location = match mapped.len() {
        0 => return None,
        1 => mapped.pop()?,
        _ => Location::Join(mapped),
    };
    Some(match direction {
        "fwd" => Location::Complement(Box::new(location)),
        "none" => Location::Order(vec![location]),
        _ => location,
    })
}
fn restore_source(molecule: &mut Molecule, bytes: &[u8], index: usize) -> Result<(), Error> {
    molecule
        .metadata
        .insert(legacy::ORIGINAL_BASE64_KEY.into(), STANDARD.encode(bytes));
    molecule.metadata.insert(
        legacy::SOURCE_HASH_KEY.into(),
        blake3::hash(bytes).to_hex().to_string(),
    );
    molecule
        .metadata
        .insert(legacy::KIND_KEY.into(), "collagene-xml".into());
    molecule
        .metadata
        .insert("legacy:record_index".into(), index.to_string());
    molecule.validate()?;
    crate::protocol::save_document(molecule)?;
    Ok(())
}
struct Count(usize);
impl Write for Count {
    fn write(&mut self, b: &[u8]) -> std::io::Result<usize> {
        self.0 = self
            .0
            .checked_add(b.len())
            .filter(|n| *n <= MAX_IMPORT_BYTES)
            .ok_or_else(|| std::io::Error::other("Collagene aggregate output bound"))?;
        Ok(b.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}
/// Atomic opt-in view of original lower-only material. Only original lower calls
/// are reversed into5′→3′; the missing upper strand is never complemented/filled.
pub fn import_with_policy(bytes: &[u8], policy: ReferencePolicy) -> Result<ImportReport, Error> {
    if policy == ReferencePolicy::UpperOnly {
        return legacy::import_collagene_xml(bytes);
    }
    let doc = document(bytes)?;
    let records = records(&doc)?;
    let mut lower_only = false;
    for record in &records {
        let (_, upper, lower) = rows(*record)?;
        lower_only |= upper.bytes().all(|b| b == b' ') && lower.bytes().any(|b| b != b' ');
    }
    if !lower_only {
        return legacy::import_collagene_xml(bytes);
    }
    // Every output retains the complete original project. Reject impossible
    // aggregate attachment sizes before making a per-record base64 copy.
    if bytes
        .len()
        .div_ceil(3)
        .checked_mul(4)
        .and_then(|n| n.checked_mul(records.len()))
        .is_none_or(|n| n > MAX_IMPORT_BYTES)
    {
        return Err(Error::TooLarge);
    }
    // The producer emits no namespaces. Unknown namespace-dependent fragments
    // cannot be safely detached from ancestor declarations for this opt-in view.
    if doc
        .descendants()
        .filter(Node::is_element)
        .flat_map(|n| n.attributes())
        .any(|a| a.namespace().is_some())
    {
        return Err(unsupported(
            "namespaced attributes are not interpreted by the lower-reference profile",
        ));
    }
    let root_unknown = doc.root_element().has_tag_name("madgene")
        && doc.root_element().attributes().next().is_some();
    let mut output=ImportReport{molecules:Vec::with_capacity(records.len()),warnings:vec!["Explicit lower-only reference policy selected. Actual lower-row calls may be reversed into5′→3′; missing strands, end chemistry, nicks and workflow are not reconstructed.".into()]};
    let mut total = Count(0);
    for (index, record) in records.into_iter().enumerate() {
        let (sequence, upper, lower) = rows(record)?;
        let projected = upper.bytes().all(|b| b == b' ') && lower.bytes().any(|b| b != b' ');
        let raw = if projected {
            temporary_lower_record(bytes, record, sequence, lower)?
        } else {
            bytes[record.range()].to_vec()
        };
        let mut parsed = legacy::import_collagene_xml(&raw)?;
        if parsed.molecules.len() != 1 {
            return Err(bad("record parser returned another scope"));
        }
        let mut molecule = parsed
            .molecules
            .pop()
            .ok_or_else(|| bad("record parser returned no sequence"))?;
        if projected {
            let mut prefix = Vec::with_capacity(lower.len() + 1);
            prefix.push(0);
            for b in lower.bytes() {
                prefix.push(prefix.last().copied().unwrap_or(0) + usize::from(b != b' '));
            }
            for feature in &mut molecule.features {
                let (from, to, direction) = source_interval(feature)?;
                let direction = direction.to_owned();
                feature.location=lower_location(from,to,&direction,&prefix,molecule.topology).unwrap_or_else(||Location::Unresolved(format!("Collagene lower-row aligned interval [{from},{to}); gapped, out-of-range or ambiguous footprint")));
                feature.qualifiers.push((
                    "legacy:projected_reference".into(),
                    Some("lower strand; reversed5′→3′".into()),
                ));
                feature.qualifiers.push((
                    "legacy:projected_orientation".into(),
                    Some(
                        match direction.as_str() {
                            "fwd" => "rev",
                            "rev" => "fwd",
                            _ => "none",
                        }
                        .into(),
                    ),
                ));
                feature.needs_review = true;
            }
            molecule
                .metadata
                .insert("legacy:upper_aligned".into(), upper.into());
            molecule
                .metadata
                .insert("legacy:lower_aligned".into(), lower.into());
            molecule
                .metadata
                .insert("legacy:reference_strand".into(), "lower-5prime".into());
            molecule
                .metadata
                .insert("genbank.molecule_type".into(), "ss-DNA".into());
            molecule
                .metadata
                .insert("source:legacy:chemistry_unmapped".into(), "true".into());
            output.warnings.push(format!("Record {} uses only actual reversed lower calls as its reference; original upper/lower rows and coordinates remain attached. Feature directions are reversed only as documentary coordinates; uncertain footprints remain unresolved. No binding, complement, duplex or reusable material is established.",index+1));
        }
        if root_unknown {
            molecule
                .metadata
                .insert("source:legacy:chemistry_unmapped".into(), "true".into());
            output.warnings.push("Unknown project-root attributes remain original-source evidence; physical assumptions are blocked.".into());
        }
        restore_source(&mut molecule, bytes, index)?;
        serde_json::to_writer(&mut total, &molecule).map_err(|_| Error::TooLarge)?;
        output.warnings.extend(parsed.warnings);
        output.molecules.push(molecule);
    }
    output.warnings.sort();
    output.warnings.dedup();
    Ok(output)
}
