//! Independently implemented Collagene XML and version-0 XDNA migration.
//! Stored source geometry/chemistry is evidence, never a trusted reconstructed duplex.
use crate::formats::sequence::{ImportReport, MAX_IMPORT_BYTES};
use base64::{Engine, engine::general_purpose::STANDARD};
use collagene_model::{Alphabet, Error, Feature, Location, Molecule, Sequence, Topology};
use roxmltree::{Document, Node, ParsingOptions};
pub const ORIGINAL_BASE64_KEY: &str = "legacy:original_base64";
pub const KIND_KEY: &str = "legacy:format";
pub const SOURCE_HASH_KEY: &str = "legacy:source_blake3";
pub const MAX_RECORDS: usize = 1000;
pub const MAX_FEATURES: usize = 100_000;
fn bad(s: impl Into<String>) -> Error {
    Error::Format(s.into())
}
fn unsupported(s: impl Into<String>) -> Error {
    Error::Unsupported(s.into())
}
fn bounded(bytes: &[u8]) -> Result<(), Error> {
    if bytes.len() > MAX_IMPORT_BYTES {
        return Err(Error::TooLarge);
    }
    Ok(())
}
fn xml(bytes: &[u8]) -> Result<Document<'_>, Error> {
    bounded(bytes)?;
    let text = std::str::from_utf8(bytes).map_err(|_| bad("legacy XML must be UTF-8"))?;
    Document::parse_with_options(
        text,
        ParsingOptions {
            allow_dtd: false,
            nodes_limit: 300_000,
            ..ParsingOptions::default()
        },
    )
    .map_err(|e| bad(format!("invalid legacy XML: {e}")))
}
pub fn is_collagene_xml(bytes: &[u8]) -> bool {
    xml(bytes).is_ok_and(|d| {
        matches!(
            d.root_element().tag_name().name(),
            "madgene" | "f" | "oneseq"
        ) || d.root_element().tag_name().name() == "sequence"
            && d.root_element()
                .children()
                .any(|n| n.has_tag_name("sequence"))
    })
}
pub fn is_xdna(bytes: &[u8]) -> bool {
    bytes.len() >= 112
        && bytes[0] == 0
        && (1..=4).contains(&bytes[1])
        && bytes[2] <= 1
        && bytes[111] == 255
}
pub fn original_legacy_bytes(molecule: &Molecule) -> Result<Vec<u8>, Error> {
    molecule.validate()?;
    let raw = molecule
        .metadata
        .get(ORIGINAL_BASE64_KEY)
        .ok_or_else(|| bad("no original legacy file retained"))?;
    if raw.len() > MAX_IMPORT_BYTES.div_ceil(3) * 4 {
        return Err(Error::TooLarge);
    }
    let bytes = STANDARD
        .decode(raw)
        .map_err(|_| bad("invalid retained legacy encoding"))?;
    bounded(&bytes)?;
    if molecule.metadata.get(SOURCE_HASH_KEY).map(String::as_str)
        != Some(blake3::hash(&bytes).to_hex().as_str())
    {
        return Err(bad("retained legacy source hash mismatch"));
    }
    Ok(bytes)
}
fn attach(molecule: &mut Molecule, bytes: &[u8], kind: &str, index: usize) -> Result<(), Error> {
    molecule
        .metadata
        .insert(ORIGINAL_BASE64_KEY.into(), STANDARD.encode(bytes));
    molecule.metadata.insert(KIND_KEY.into(), kind.into());
    molecule.metadata.insert(
        SOURCE_HASH_KEY.into(),
        blake3::hash(bytes).to_hex().to_string(),
    );
    molecule
        .metadata
        .insert("legacy:record_index".into(), index.to_string());
    molecule.validate()?;
    crate::protocol::save_document(molecule)?;
    Ok(())
}
fn attr<'a>(node: Node<'a, 'a>, name: &str) -> Result<&'a str, Error> {
    node.attribute(name)
        .ok_or_else(|| bad(format!("legacy {} lacks {name}", node.tag_name().name())))
}
fn integer(s: &str) -> Result<i64, Error> {
    if s.len() > 20 || s.trim() != s || s.is_empty() {
        return Err(bad("invalid legacy signed coordinate"));
    }
    s.parse().map_err(|_| bad("invalid legacy signed integer"))
}
fn orient(s: &str) -> Result<&str, Error> {
    if matches!(s, "fwd" | "rev" | "none") {
        Ok(s)
    } else {
        Err(unsupported("unknown legacy orientation"))
    }
}
fn children<'a, 'input: 'a>(node: Node<'a, 'input>) -> impl Iterator<Item = Node<'a, 'input>> {
    node.children().filter(Node::is_element)
}
fn text(s: &str, limit: usize, label: &str) -> Result<String, Error> {
    if s.len() > limit
        || s.chars()
            .any(|c| c.is_control() && !matches!(c, '\t' | '\r' | '\n'))
    {
        return Err(bad(format!("invalid {label}")));
    }
    Ok(s.to_string())
}
fn old_location(
    from: i64,
    to: i64,
    columns: usize,
    prefix: &[usize],
    topology: Topology,
    direction: &str,
    warnings: &mut Vec<String>,
) -> Result<Location, Error> {
    let n = columns as i64;
    let (from, to) = if from < 0 || to < 0 || from > n || to > n {
        if topology != Topology::Circular
            || from.unsigned_abs() > (n * 2) as u64
            || to.unsigned_abs() > (n * 2) as u64
        {
            warnings.push(
                "Legacy range lies outside supported reference columns; retained unresolved."
                    .into(),
            );
            return Ok(Location::Unresolved(format!(
                "Collagene aligned columns [{from},{to})"
            )));
        }
        let length = if to >= from { to - from } else { to - from + n };
        if length < 0 || length > n {
            warnings.push("Multiple-turn legacy range retained unresolved.".into());
            return Ok(Location::Unresolved(format!(
                "Collagene multiple-turn [{from},{to})"
            )));
        }
        let start = from.rem_euclid(n);
        (start, start + length)
    } else {
        (from, to)
    };
    let spans = if from == to {
        warnings.push("Legacy equal-endpoint/empty range retained unresolved.".into());
        return Ok(Location::Unresolved(format!(
            "Collagene empty [{from},{to})"
        )));
    } else if to > n {
        vec![(from as usize, columns), (0, (to - n) as usize)]
    } else if from < to {
        vec![(from as usize, to as usize)]
    } else if topology == Topology::Circular {
        vec![(from as usize, columns), (0, to as usize)]
    } else {
        warnings.push("Wrapped annotation on linear legacy sequence retained unresolved.".into());
        return Ok(Location::Unresolved(format!(
            "Collagene linear wrapped [{from},{to})"
        )));
    };
    let mut parts = vec![];
    for (a, b) in spans {
        let (a, b) = (prefix[a], prefix[b]);
        if a < b {
            parts.push(Location::range(a, b));
        }
    }
    if parts.is_empty() {
        warnings.push(
            "Legacy annotation covers only absent upper-strand columns; retained unresolved."
                .into(),
        );
        return Ok(Location::Unresolved(
            "Collagene lower-strand-only interval".into(),
        ));
    }
    let location = if parts.len() == 1 {
        parts.remove(0)
    } else {
        Location::Join(parts)
    };
    Ok(match direction {
        "rev" => Location::Complement(Box::new(location)),
        "none" => {
            warnings.push("Unoriented legacy feature retains spans but not a biological extraction direction.".into());
            Location::Order(vec![location])
        }
        _ => location,
    })
}
fn xml_one(
    node: Node<'_, '_>,
    bytes: &[u8],
    index: usize,
    warnings: &mut Vec<String>,
) -> Result<Molecule, Error> {
    let sequences = children(node)
        .filter(|child| child.has_tag_name("sequence"))
        .collect::<Vec<_>>();
    if sequences.len() != 1 {
        return Err(bad("legacy record requires exactly one sequence child"));
    }
    let sequence = sequences[0];
    let name = text(attr(sequence, "name")?, 1024, "legacy sequence name")?;
    let upper = attr(sequence, "su")?;
    let lower = attr(sequence, "sl")?;
    if upper.len() != lower.len()
        || upper.is_empty()
        || upper.len() > collagene_model::MAX_SEQUENCE_LEN
    {
        return Err(bad(
            "legacy aligned strands differ in length or exceed limits",
        ));
    }
    if !upper.is_ascii() || !lower.is_ascii() {
        return Err(unsupported(
            "legacy strand cells must be ASCII DNA or space gaps",
        ));
    }
    let mut bases = String::new();
    let mut prefix = Vec::with_capacity(upper.len() + 1);
    prefix.push(0);
    for (position, base) in upper.bytes().enumerate() {
        if base != b' ' {
            let normalized = base.to_ascii_uppercase();
            if !b"ACGTRYSWKMBDHVN".contains(&normalized) {
                return Err(Error::InvalidBase {
                    alphabet: Alphabet::Dna,
                    character: base as char,
                    position: position + 1,
                });
            }
            bases.push(normalized as char);
        }
        prefix.push(bases.len());
    }
    for (position, base) in lower.bytes().enumerate() {
        if base != b' ' && !b"ACGTRYSWKMBDHVN".contains(&base.to_ascii_uppercase()) {
            return Err(Error::InvalidBase {
                alphabet: Alphabet::Dna,
                character: base as char,
                position: position + 1,
            });
        }
    }
    let topology = match attr(sequence, "circular")? {
        "true" => Topology::Circular,
        "false" => Topology::Linear,
        _ => return Err(bad("legacy circular flag is not true/false")),
    };
    let mut molecule = Molecule::dna(name, &bases, topology)?;
    molecule
        .metadata
        .insert("legacy:upper_aligned".into(), upper.into());
    molecule
        .metadata
        .insert("legacy:lower_aligned".into(), lower.into());
    if let Some(notes) = sequence.attribute("notes") {
        molecule.metadata.insert(
            "description".into(),
            text(notes, 1_048_576, "legacy notes")?,
        );
    }
    let actual = upper
        .bytes()
        .map(|b| collagene_model::complement(b.to_ascii_uppercase(), Alphabet::Dna))
        .collect::<Vec<_>>();
    if upper.contains(' ') || lower.contains(' ') || actual != lower.to_ascii_uppercase().as_bytes()
    {
        molecule
            .metadata
            .insert("genbank.molecule_type".into(), "ms-DNA".into());
        warnings.push("Legacy aligned strands contain gaps/mismatches. Upper-reference coordinates are projected, both source strands remain evidence, and intact-duplex scientific models are blocked; no terminal chemistry/nick state is inferred.".into());
    }
    let mut count = 0;
    for child in children(node) {
        match child.tag_name().name(){"sequence"=>{},"annotation"|"primer"=>{
        count+=1;if count>MAX_FEATURES{return Err(Error::TooLarge);}let is_primer=child.has_tag_name("primer");let name=text(attr(child,"name")?,1024,"legacy feature name")?;let direction=orient(attr(child,"orientation")?)?;
        let mut qualifiers=vec![("label".into(),Some(name)),("legacy:orientation".into(),Some(direction.into()))];
        let(from,to)=if is_primer{let ordered=Sequence::dna(attr(child,"sequence")?)?;let target=integer(attr(child,"target")?)?;let len=ordered.len()as i64;qualifiers.push(("oligo_sequence".into(),Some(ordered.as_str().into())));qualifiers.push(("legacy:primer_target".into(),Some(target.to_string())));warnings.push("Legacy primer records have a target and full oligo but no separate tail/annealing or binding proof. The oligo remains evidence and inferred primer features require review.".into());if direction=="fwd"{(target.checked_sub(len).ok_or(Error::TooLarge)?,target)}else if direction=="rev"{(target,target.checked_add(len).ok_or(Error::TooLarge)?)}else{(0,0)}}else{let from=integer(attr(child,"from")?)?;let to=integer(attr(child,"to")?)?;let mut rgb=vec![];for key in ["colr","colg","colb"]{let value=integer(attr(child,key)?)?;if !(0..=255).contains(&value){return Err(bad("legacy color component is outside 0..255"));}rgb.push(value as u8);}qualifiers.push(("color".into(),Some(format!("#{:02x}{:02x}{:02x}",rgb[0],rgb[1],rgb[2]))));(from,to)};
        qualifiers.push(("legacy:aligned_range".into(),Some(format!("[{from},{to})"))));let location=old_location(from,to,upper.len(),&prefix,topology,direction,warnings)?;
        molecule.features.push(Feature{kind:if is_primer{"primer_bind"}else{"misc_feature"}.into(),location,qualifiers,needs_review:true});
    },unknown=>warnings.push(format!("Legacy child {unknown:?} has no interpreted migration semantics; its complete XML remains attached."))}
    }
    // Preserve all unknown attributes/elements in the complete original file; expose omissions.
    let known = ["name", "su", "sl", "circular", "notes"];
    if sequence
        .attributes()
        .any(|attribute| !known.contains(&attribute.name()))
    {
        warnings.push("Unknown legacy sequence attributes remain original-source evidence.".into());
    }
    let unknown_material = sequence
        .attributes()
        .any(|attribute| !known.contains(&attribute.name()))
        || children(sequence).next().is_some()
        || node.attributes().next().is_some()
        || children(node).any(|child| {
            !matches!(
                child.tag_name().name(),
                "sequence" | "annotation" | "primer"
            )
        });
    if unknown_material {
        molecule
            .metadata
            .insert("source:legacy:chemistry_unmapped".into(), "true".into());
        warnings.push("Uninterpreted legacy material fields may describe physical states. Source data are retained, but intact-duplex assumptions are blocked until those fields are modeled.".into());
    }
    attach(&mut molecule, bytes, "collagene-xml", index)?;
    Ok(molecule)
}
pub fn import_collagene_xml(bytes: &[u8]) -> Result<ImportReport, Error> {
    let doc = xml(bytes)?;
    let root = doc.root_element();
    let mut traversal = vec![(root, 1usize)];
    while let Some((node, depth)) = traversal.pop() {
        if depth > 64 {
            return Err(unsupported("legacy XML nesting exceeds 64 levels"));
        }
        traversal.extend(children(node).map(|child| (child, depth + 1)));
    }
    if root.tag_name().namespace().is_some() {
        return Err(unsupported(
            "namespaced legacy XML is not the original serializer schema",
        ));
    }
    let records = match root.tag_name().name() {
        "madgene" => {
            let records = children(root)
                .filter(|child| child.has_tag_name("oneseq"))
                .collect::<Vec<_>>();
            if children(root).any(|child| !child.has_tag_name("oneseq")) {
                return Err(unsupported(
                    "unknown legacy project record type; no project was partially imported",
                ));
            }
            records
        }
        "f" | "oneseq" => vec![root],
        "sequence" if children(root).any(|n| n.has_tag_name("sequence")) => vec![root],
        _ => return Err(unsupported("unrecognized Collagene XML root")),
    };
    if records.is_empty() || records.len() > MAX_RECORDS {
        return Err(bad("legacy XML record count is empty or exceeds 1000"));
    }
    let mut warnings=vec!["Collagene XML migrates serialized sequence pairs, notes, annotations/colors/orientations and primer evidence. It did not serialize typed cloning history, terminal chemistry or nicks; these are not invented. Complete input bytes remain attached.".into()];
    for node in doc.descendants().filter(Node::is_element) {
        if node.tag_name().namespace().is_some() {
            return Err(unsupported(
                "namespaced legacy elements require a separate schema adapter",
            ));
        }
        if node.tag_name().name().len() > 128 {
            return Err(Error::TooLarge);
        }
        let known: &[&str] = match node.tag_name().name() {
            "sequence" => &["name", "su", "sl", "circular", "notes"],
            "annotation" => &["name", "from", "to", "colr", "colg", "colb", "orientation"],
            "primer" => &["name", "sequence", "target", "orientation"],
            _ => &[],
        };
        if node
            .attributes()
            .any(|attribute| attribute.namespace().is_some() || !known.contains(&attribute.name()))
        {
            warnings.push(format!("Unknown attributes on legacy {} remain in original-source XML; they are not interpreted.",node.tag_name().name()));
        }
        if matches!(node.tag_name().name(), "sequence" | "annotation" | "primer")
            && children(node).next().is_some()
        {
            warnings.push(format!(
                "Nested metadata on legacy {} remains original-source evidence only.",
                node.tag_name().name()
            ));
        }
    }
    let mut molecules = vec![];
    let mut total = 0usize;
    for (index, record) in records.into_iter().enumerate() {
        let molecule = xml_one(record, bytes, index, &mut warnings)?;
        total += serde_json::to_vec(&molecule)
            .map_err(|e| bad(e.to_string()))?
            .len();
        if total > MAX_IMPORT_BYTES {
            return Err(Error::TooLarge);
        }
        molecules.push(molecule);
    }
    warnings.sort();
    warnings.dedup();
    Ok(ImportReport {
        molecules,
        warnings,
    })
}

struct Cursor<'a> {
    bytes: &'a [u8],
    position: usize,
}
impl<'a> Cursor<'a> {
    fn take(&mut self, len: usize) -> Result<&'a [u8], Error> {
        let end = self.position.checked_add(len).ok_or(Error::TooLarge)?;
        let bytes = self
            .bytes
            .get(self.position..end)
            .ok_or_else(|| bad("truncated XDNA field"))?;
        self.position = end;
        Ok(bytes)
    }
    fn byte(&mut self) -> Result<u8, Error> {
        Ok(self.take(1)?[0])
    }
    fn string(&mut self) -> Result<String, Error> {
        let len = self.byte()? as usize;
        ascii(self.take(len)?)
    }
    fn overhang(&mut self) -> Result<(i64, String), Error> {
        let value = integer(&self.string()?)?;
        if value.unsigned_abs() > collagene_model::MAX_SEQUENCE_LEN as u64 {
            return Err(Error::TooLarge);
        }
        let sequence = ascii(self.take(value.unsigned_abs() as usize)?)?;
        if !sequence.is_empty() {
            Sequence::dna(&sequence)?;
        }
        Ok((value, sequence))
    }
}
fn ascii(bytes: &[u8]) -> Result<String, Error> {
    if !bytes.is_ascii() {
        return Err(unsupported(
            "XDNA text encoding is not ASCII; original bytes are required for an explicit encoding adapter",
        ));
    }
    String::from_utf8(bytes.to_vec()).map_err(|_| bad("invalid XDNA text"))
}
pub fn import_xdna(bytes: &[u8], name: &str) -> Result<ImportReport, Error> {
    bounded(bytes)?;
    if !is_xdna(bytes) {
        return Err(unsupported("unsupported XDNA header/version/type/topology"));
    }
    let length = u32::from_be_bytes(bytes[28..32].try_into().unwrap()) as usize;
    let origin = u32::from_be_bytes(bytes[32..36].try_into().unwrap()) as usize;
    let comment_length = u32::from_be_bytes(bytes[96..100].try_into().unwrap()) as usize;
    if length == 0 || length > collagene_model::MAX_SEQUENCE_LEN || origin > length {
        return Err(bad("XDNA length/origin is invalid"));
    }
    let mut cursor = Cursor {
        bytes,
        position: 112,
    };
    let raw_sequence = ascii(cursor.take(length)?)?;
    let alphabet = match bytes[1] {
        1 | 2 => Alphabet::Dna,
        3 => Alphabet::Rna,
        4 => Alphabet::Protein,
        _ => unreachable!(),
    };
    let sequence = Sequence::new(alphabet, &raw_sequence)?;
    if sequence.len() != length {
        return Err(bad(
            "XDNA sequence contains whitespace rather than declared residues",
        ));
    }
    let comment = ascii(cursor.take(comment_length)?)?;
    let topology = if bytes[2] == 1 {
        Topology::Circular
    } else {
        Topology::Linear
    };
    if alphabet == Alphabet::Protein && topology == Topology::Circular {
        return Err(unsupported(
            "circular protein XDNA is not a supported biological topology",
        ));
    }
    let mut molecule = Molecule {
        name: text(name, 1024, "XDNA name")?,
        sequence,
        topology,
        features: vec![],
        metadata: Default::default(),
        original_genbank: None,
    };
    molecule.metadata.insert("description".into(), comment);
    molecule
        .metadata
        .insert("legacy:xdna_origin_offset".into(), origin.to_string());
    let mut warnings=vec!["XDNA has no authenticated source, cloning history or free-end chemistry. Header padding and all original bytes are retained; nonzero negative-length origin is recorded without rotating bases or annotation coordinates.".into()];
    if cursor.position < bytes.len() {
        let marker = cursor.byte()?;
        molecule
            .metadata
            .insert("legacy:xdna_extension_marker".into(), marker.to_string());
        let (right, right_sequence) = cursor.overhang()?;
        let (left, left_sequence) = cursor.overhang()?;
        for (key, value) in [
            ("right_length", right.to_string()),
            ("right_sequence", right_sequence),
            ("left_length", left.to_string()),
            ("left_sequence", left_sequence),
        ] {
            molecule
                .metadata
                .insert(format!("legacy:xdna_{key}"), value);
        }
        if right != 0 || left != 0 {
            if alphabet != Alphabet::Dna {
                return Err(unsupported("XDNA overhangs require DNA"));
            }
            molecule
                .metadata
                .insert("genbank.molecule_type".into(), "ms-DNA".into());
            warnings.push("XDNA overhang descriptors remain unmapped physical source evidence. Intact-duplex modeling is blocked; no ligatable ends or phosphorylation are asserted.".into());
        }
        let count = cursor.byte()? as usize;
        for _ in 0..count {
            let label = cursor.string()?;
            let description = cursor.string()?;
            let kind = cursor.string()?;
            if !kind.is_empty()
                && !kind
                    .bytes()
                    .all(|byte| byte.is_ascii_alphanumeric() || b"_:-".contains(&byte))
            {
                return Err(unsupported(
                    "XDNA feature type is not a supported GenBank-style key",
                ));
            }
            let a = integer(&cursor.string()?)?;
            let b = integer(&cursor.string()?)?;
            let strand = cursor.byte()?;
            let display = cursor.byte()?;
            let unknown = cursor.byte()?;
            let arrow = cursor.byte()?;
            let extra = cursor.string()?;
            if strand > 1 || display > 1 || arrow > 1 {
                return Err(unsupported("unsupported XDNA feature flags"));
            }
            let (start, end) = if strand == 1 { (a, b) } else { (b, a) };
            if start < 1 || end < 1 || start > length as i64 || end > length as i64 {
                return Err(bad("XDNA feature coordinates are outside sequence"));
            }
            let mut parts = if start <= end {
                vec![Location::range((start - 1) as usize, end as usize)]
            } else if topology == Topology::Circular {
                vec![
                    Location::range((start - 1) as usize, length),
                    Location::range(0, end as usize),
                ]
            } else {
                return Err(bad("XDNA wrapped feature requires circular sequence"));
            };
            let location = if parts.len() == 1 {
                parts.remove(0)
            } else {
                Location::Join(parts)
            };
            let location = if origin != 0 {
                warnings.push("XDNA feature coordinates with a nonzero display-origin offset require a producer fixture; source values are retained unresolved instead of guessing their anchor.".into());
                Location::Unresolved(format!(
                    "XDNA origin-offset {origin}; source feature [{a},{b}], strand {strand}"
                ))
            } else if strand == 0 && alphabet == Alphabet::Protein {
                warnings.push("Reverse-strand flags on protein XDNA features have no mapped biological meaning; source coordinates remain unresolved.".into());
                Location::Unresolved(format!(
                    "XDNA protein reverse flag; source feature [{a},{b}]"
                ))
            } else if strand == 0 {
                Location::Complement(Box::new(location))
            } else {
                location
            };
            let mut qualifiers = vec![
                ("label".into(), Some(label)),
                ("note".into(), Some(description)),
                (
                    "legacy:xdna_flags".into(),
                    Some(format!("{strand},{display},{unknown},{arrow}")),
                ),
                ("legacy:xdna_extra".into(), Some(extra.clone())),
                ("legacy:xdna_coordinates".into(), Some(format!("{a},{b}"))),
            ];
            let color = extra
                .trim_end_matches(',')
                .split(',')
                .map(str::parse::<u8>)
                .collect::<Result<Vec<_>, _>>();
            if let Ok(rgb) = color {
                if rgb.len() == 3 {
                    qualifiers.push((
                        "color".into(),
                        Some(format!("#{:02x}{:02x}{:02x}", rgb[0], rgb[1], rgb[2])),
                    ));
                }
            }
            molecule.features.push(Feature {
                kind: if kind.is_empty() {
                    "misc_feature".into()
                } else {
                    kind
                },
                location,
                qualifiers,
                needs_review: true,
            });
        }
        if cursor.position != bytes.len() {
            return Err(unsupported(
                "unrecognized trailing XDNA extension data; no partial file was imported",
            ));
        }
        warnings.push("XDNA feature descriptions/flags/colors are retained. GenBank-like text within descriptions remains literal source evidence rather than unvalidated qualifier/phase claims.".into());
    }
    attach(&mut molecule, bytes, "xdna", 0)?;
    Ok(ImportReport {
        molecules: vec![molecule],
        warnings,
    })
}

pub mod collagene_lower_strand;
pub mod xdna_export;
