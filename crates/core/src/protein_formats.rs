//! Independent annotated protein flat-file adapters. Source text is documentary;
//! these formats never establish native physical-material or verification authority.
use crate::formats::{ExportReport, ImportReport};
use crate::{Alphabet, Error, Feature, Location, Molecule, Sequence, Topology};
use base64::{Engine as _, engine::general_purpose::STANDARD};
use std::{collections::BTreeMap, fmt::Write};

pub const MAX_BYTES: usize = 4 * 1024 * 1024;
pub const MAX_RECORD_BYTES: usize = 1024 * 1024;
pub const MAX_RECORDS: usize = 128;
const MAX_FEATURES: usize = 4096;
const MAX_QUALIFIERS: usize = 16384;
const MAX_LINE: usize = 65536;
const CHUNK: usize = 65536;
const ORIGINAL: &str = "protein:original:";
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Profile {
    GenPept,
    UniProt,
}
impl Profile {
    pub fn name(self) -> &'static str {
        match self {
            Self::GenPept => "genpept",
            Self::UniProt => "uniprot",
        }
    }
}
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct InsdcDate(String);
impl InsdcDate {
    pub fn new(value: &str) -> Result<Self, Error> {
        let parts: Vec<_> = value.split('-').collect();
        if parts.len() != 3 || parts[0].len() != 2 || parts[2].len() != 4 {
            return Err(bad("INSDC date must be DD-MMM-YYYY"));
        }
        let day = parts[0]
            .parse::<u32>()
            .map_err(|_| bad("invalid INSDC date"))?;
        let year = parts[2]
            .parse::<u32>()
            .map_err(|_| bad("invalid INSDC date"))?;
        let month = [
            "JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC",
        ]
        .iter()
        .position(|m| *m == parts[1])
        .ok_or_else(|| bad("invalid INSDC month"))?;
        let leap = year % 4 == 0 && (year % 100 != 0 || year % 400 == 0);
        let days = [
            31,
            if leap { 29 } else { 28 },
            31,
            30,
            31,
            30,
            31,
            31,
            30,
            31,
            30,
            31,
        ][month];
        if year == 0 || day == 0 || day > days {
            return Err(bad("invalid calendar date"));
        }
        Ok(Self(value.into()))
    }
    pub fn as_str(&self) -> &str {
        &self.0
    }
}
#[derive(Debug, Clone, Default)]
pub struct ExportOptions {
    /// Caller-selected actual export date, not a biological event timestamp.
    pub locus_date: Option<InsdcDate>,
    /// Explicit UniProt identifier. No database allocation/review is inferred.
    pub uniprot_accession: Option<String>,
    pub entry_name: Option<String>,
}
const UNIPROT_KEYS: &[&str] = &[
    "INIT_MET", "SIGNAL", "PROPEP", "TRANSIT", "CHAIN", "PEPTIDE", "TOPO_DOM", "TRANSMEM",
    "INTRAMEM", "DOMAIN", "REPEAT", "CA_BIND", "ZN_FING", "DNA_BIND", "NP_BIND", "REGION",
    "COILED", "MOTIF", "COMPBIAS", "ACT_SITE", "BINDING", "SITE", "NON_STD", "MOD_RES", "LIPID",
    "CARBOHYD", "DISULFID", "CROSSLNK", "VAR_SEQ", "VARIANT", "MUTAGEN", "UNSURE", "CONFLICT",
    "NON_CONS", "NON_TER", "HELIX", "STRAND", "TURN",
];
fn known_key(profile: Profile, key: &str) -> bool {
    match profile {
        Profile::UniProt => UNIPROT_KEYS.contains(&key),
        Profile::GenPept => [
            "source",
            "Protein",
            "Region",
            "Site",
            "CDS",
            "Bond",
            "bond",
            "mat_peptide",
            "sig_peptide",
            "propeptide",
        ]
        .contains(&key),
    }
}
fn endpoint_feature(kind: &str) -> bool {
    matches!(kind, "DISULFID" | "CROSSLNK")
}
fn bad(message: &str) -> Error {
    Error::Format(format!("protein flat file: {message}"))
}
fn unsupported(message: &str) -> Error {
    Error::Unsupported(format!(
        "{message}; retain the native/original record or export protein FASTA"
    ))
}
/// Content hint only. Actual import validates every record under the selected profile.
pub fn detect(bytes: &[u8]) -> Option<Profile> {
    // Recognize up to the generic import bound, even when the stricter protein
    // parser will reject the record, rather than falling through to DNA parsing.
    if bytes.len() > crate::formats::MAX_IMPORT_BYTES {
        return None;
    }
    let text = std::str::from_utf8(bytes).ok()?;
    let line = text
        .trim_start_matches('\u{feff}')
        .lines()
        .find(|l| !l.trim().is_empty())?;
    if line.starts_with("LOCUS ") && line.split_whitespace().any(|v| v == "aa") {
        Some(Profile::GenPept)
    } else if line.starts_with("ID   ") && line.trim_end().ends_with("AA.") {
        Some(Profile::UniProt)
    } else {
        None
    }
}
fn retain(m: &mut Molecule, bytes: &[u8], profile: Profile) {
    let encoded = STANDARD.encode(bytes);
    m.metadata.insert(
        format!("{ORIGINAL}chunks"),
        encoded.len().div_ceil(CHUNK).to_string(),
    );
    m.metadata.insert(
        format!("{ORIGINAL}blake3"),
        blake3::hash(bytes).to_hex().to_string(),
    );
    m.metadata
        .insert("protein:profile".into(), profile.name().into());
    for (index, chunk) in encoded.as_bytes().chunks(CHUNK).enumerate() {
        m.metadata.insert(
            format!("{ORIGINAL}{index:04}"),
            std::str::from_utf8(chunk).expect("base64").into(),
        );
    }
}
/// Immutable source attachment, not a normalized export of the current molecule.
pub fn original_bytes(molecule: &Molecule) -> Result<Vec<u8>, Error> {
    let count: usize = molecule
        .metadata
        .get(&format!("{ORIGINAL}chunks"))
        .ok_or_else(|| bad("no original protein record"))?
        .parse()
        .map_err(|_| bad("original chunk count"))?;
    if count == 0
        || count
            > MAX_RECORD_BYTES
                .div_ceil(3)
                .saturating_mul(4)
                .div_ceil(CHUNK)
    {
        return Err(Error::TooLarge);
    }
    let mut encoded = String::new();
    for index in 0..count {
        let part = molecule
            .metadata
            .get(&format!("{ORIGINAL}{index:04}"))
            .ok_or_else(|| bad("missing original chunk"))?;
        if part.len() > CHUNK || (index + 1 < count && part.len() != CHUNK) {
            return Err(bad("original chunk size"));
        }
        encoded.push_str(part);
    }
    let bytes = STANDARD
        .decode(&encoded)
        .map_err(|_| bad("original encoding"))?;
    if bytes.len() > MAX_RECORD_BYTES || STANDARD.encode(&bytes) != encoded {
        return Err(Error::TooLarge);
    }
    if molecule
        .metadata
        .get(&format!("{ORIGINAL}blake3"))
        .map(String::as_str)
        != Some(blake3::hash(&bytes).to_hex().as_str())
    {
        return Err(bad("original digest mismatch"));
    }
    Ok(bytes)
}
#[derive(Default)]
struct Builder {
    kind: String,
    location: String,
    qualifiers: Vec<String>,
    legacy_note: String,
    review: bool,
    legacy: bool,
}
fn quoted_complete(value: &str) -> bool {
    let Some((_, value)) = value.split_once('=') else {
        return true;
    };
    let value = value.trim();
    if !value.starts_with('"') {
        return true;
    }
    let bytes = value.as_bytes();
    let mut index = 1;
    while index < bytes.len() {
        if bytes[index] == b'\\' && bytes.get(index + 1) == Some(&b'"') {
            index += 2;
            continue;
        }
        if bytes[index] == b'"' {
            if bytes.get(index + 1) == Some(&b'"') {
                index += 2;
                continue;
            }
            return value[index + 1..].trim().is_empty();
        }
        index += 1;
    }
    false
}
impl Builder {
    fn continuation(&mut self, value: &str, total_qualifiers: &mut usize) -> Result<(), Error> {
        if self.qualifiers.last().is_some_and(|v| !quoted_complete(v)) {
            let last = self.qualifiers.last_mut().expect("checked");
            last.push(' ');
            last.push_str(value);
        } else if value.starts_with('/') {
            *total_qualifiers += 1;
            if *total_qualifiers > MAX_QUALIFIERS || self.qualifiers.len() >= MAX_QUALIFIERS {
                return Err(Error::TooLarge);
            }
            self.qualifiers.push(value.into());
        } else if self.legacy || !self.legacy_note.is_empty() || !self.qualifiers.is_empty() {
            self.review = true;
            if !self.legacy_note.is_empty() {
                self.legacy_note.push(' ');
            }
            self.legacy_note.push_str(value);
        } else {
            self.location.push_str(value);
        }
        Ok(())
    }
    fn finish(self, length: usize, profile: Profile) -> Result<Feature, Error> {
        let mut review = self.review;
        let mut location = parse_location(&self.location, length, profile, 0)?;
        review |= !known_key(profile, &self.kind);
        if profile == Profile::UniProt && endpoint_feature(&self.kind) {
            location = match location {
                Location::Range {
                    start,
                    end,
                    before: false,
                    after: false,
                } if end > start + 1 => Location::Join(vec![
                    Location::Range {
                        start,
                        end: start + 1,
                        before: false,
                        after: false,
                    },
                    Location::Range {
                        start: end - 1,
                        end,
                        before: false,
                        after: false,
                    },
                ]),
                _ => Location::Unresolved(self.location.clone()),
            };
        }
        review |= !exact_location(&location);
        let mut qualifiers = vec![("protein:source_location".into(), Some(self.location))];
        if !self.legacy_note.is_empty() {
            qualifiers.push(("note".into(), Some(self.legacy_note)));
        }
        for raw in self.qualifiers {
            if !quoted_complete(&raw) {
                review = true;
                qualifiers.push(("protein:unparsed_qualifier".into(), Some(raw)));
                continue;
            }
            let field = raw.trim_start_matches('/');
            let (key, value) = field
                .split_once('=')
                .map_or((field, None), |(k, v)| (k, Some(v.trim())));
            if !token(key) {
                review = true;
                qualifiers.push(("protein:unparsed_qualifier".into(), Some(raw)));
                continue;
            }
            let value = value.map(|v| {
                v.strip_prefix('"')
                    .and_then(|v| v.strip_suffix('"'))
                    .unwrap_or(v)
                    .replace("\"\"", "\"")
                    .replace("\\\"", "\"")
            });
            qualifiers.push((key.into(), value));
        }
        Ok(Feature {
            kind: self.kind,
            location,
            qualifiers,
            needs_review: review,
        })
    }
}
fn exact_location(location: &Location) -> bool {
    match location {
        Location::Range { before, after, .. } => !before && !after,
        Location::Join(v) | Location::Order(v) => v.iter().all(exact_location),
        _ => false,
    }
}
fn parse_location(
    raw: &str,
    length: usize,
    profile: Profile,
    depth: usize,
) -> Result<Location, Error> {
    if raw.len() > 65536 || depth > 16 {
        return Err(Error::TooLarge);
    }
    if profile == Profile::GenPept {
        for (prefix, order) in [("join(", false), ("order(", true)] {
            if let Some(inner) = raw.strip_prefix(prefix).and_then(|v| v.strip_suffix(')')) {
                let mut level = 0usize;
                let mut start = 0;
                let mut parts = Vec::new();
                for (i, byte) in inner.bytes().enumerate() {
                    match byte {
                        b'(' => {
                            level += 1;
                            if level + depth > 16 {
                                return Err(Error::TooLarge);
                            }
                        }
                        b')' => {
                            let Some(v) = level.checked_sub(1) else {
                                return Ok(Location::Unresolved(raw.into()));
                            };
                            level = v;
                        }
                        b',' if level == 0 => {
                            parts.push(&inner[start..i]);
                            start = i + 1;
                        }
                        _ => (),
                    }
                }
                parts.push(&inner[start..]);
                if parts.len() > 1024 {
                    return Err(Error::TooLarge);
                }
                if level != 0 || parts.iter().any(|v| v.is_empty()) {
                    return Ok(Location::Unresolved(raw.into()));
                }
                let locations = parts
                    .into_iter()
                    .map(|v| parse_location(v, length, profile, depth + 1))
                    .collect::<Result<Vec<_>, _>>()?;
                return Ok(if order {
                    Location::Order(locations)
                } else {
                    Location::Join(locations)
                });
            }
        }
    }
    let (left, right) = raw.split_once("..").unwrap_or((raw, raw));
    let before = left.starts_with('<');
    let after = right.starts_with('>');
    let a = left.strip_prefix('<').unwrap_or(left).parse::<usize>();
    let b = right.strip_prefix('>').unwrap_or(right).parse::<usize>();
    match (a, b) {
        (Ok(a), Ok(b)) if a > 0 && a <= b && b <= length => Ok(Location::Range {
            start: a - 1,
            end: b,
            before,
            after,
        }),
        _ => Ok(Location::Unresolved(raw.into())),
    }
}
fn token(value: &str) -> bool {
    !value.is_empty()
        && value.len() <= 128
        && value
            .bytes()
            .all(|b| b.is_ascii_alphanumeric() || matches!(b, b'_' | b'-' | b'.'))
}
fn sequence_line(sequence: &mut String, line: &str, numbered: bool) -> Result<(), Error> {
    let mut words = line.split_whitespace();
    if numbered {
        let position = words
            .next()
            .ok_or_else(|| bad("missing sequence counter"))?
            .parse::<usize>()
            .map_err(|_| bad("invalid sequence counter"))?;
        if position != sequence.len() + 1 {
            return Err(bad("sequence counter disagrees with residue position"));
        }
    }
    for word in words {
        for b in word.bytes() {
            if !b.is_ascii_alphabetic() {
                return Err(bad("non-residue sequence character"));
            }
            sequence.push(b.to_ascii_uppercase() as char);
        }
    }
    Ok(())
}
fn declared_length(words: &[&str], unit: &str) -> Result<usize, Error> {
    let unit_index = words
        .iter()
        .position(|v| v.trim_end_matches([';', '.']) == unit)
        .ok_or_else(|| bad("missing amino-acid length/unit"))?;
    let length = unit_index
        .checked_sub(1)
        .and_then(|i| words[i].parse::<usize>().ok())
        .filter(|v| *v > 0)
        .ok_or_else(|| bad("invalid amino-acid length"))?;
    if length > MAX_RECORD_BYTES {
        return Err(Error::TooLarge);
    }
    Ok(length)
}
fn parse_record(bytes: &[u8], profile: Profile) -> Result<Molecule, Error> {
    let text = std::str::from_utf8(bytes).map_err(|_| bad("UTF-8 text required"))?;
    let mut lines = text
        .trim_start_matches('\u{feff}')
        .lines()
        .filter(|l| !l.trim().is_empty());
    let header = lines.next().ok_or_else(|| bad("empty record"))?;
    let words: Vec<_> = header.split_whitespace().collect();
    let length = declared_length(
        &words,
        if profile == Profile::GenPept {
            "aa"
        } else {
            "AA"
        },
    )?;
    let entry = words
        .get(1)
        .ok_or_else(|| bad("missing entry identifier"))?
        .to_string();
    if !token(&entry)
        || (profile == Profile::GenPept && !header.starts_with("LOCUS "))
        || (profile == Profile::UniProt && !header.starts_with("ID   "))
    {
        return Err(bad("wrong protein profile/header"));
    }
    if words.iter().any(|v| v.eq_ignore_ascii_case("circular")) {
        return Err(unsupported("circular protein records are unsupported"));
    }
    let mut sequence = String::new();
    let mut builders: Vec<Builder> = Vec::new();
    let mut total_qualifiers = 0;
    let mut in_features = false;
    let mut in_sequence = false;
    let mut description = String::new();
    let mut accessions = Vec::new();
    let mut metadata = BTreeMap::new();
    if profile == Profile::GenPept {
        if let Some(date) = words.last().and_then(|v| InsdcDate::new(v).ok()) {
            metadata.insert("protein:locus_date".into(), date.0);
        }
    }

    let mut crc = None;
    let mut sq = false;
    let mut header_field = String::new();
    for raw in lines {
        let line = raw.trim_end_matches('\r');
        if line.trim() == "//" {
            break;
        }
        if in_sequence {
            sequence_line(&mut sequence, line, profile == Profile::GenPept)?;
            continue;
        }
        if profile == Profile::GenPept {
            if line.starts_with("FEATURES ") {
                in_features = true;
                continue;
            }
            if line.trim_end() == "ORIGIN" {
                in_sequence = true;
                continue;
            }
            if in_features && line.starts_with("     ") {
                let kind = line
                    .get(5..21)
                    .ok_or_else(|| bad("short GenPept feature line"))?
                    .trim();
                let data = line.get(21..).unwrap_or("").trim();
                if kind.is_empty() {
                    builders
                        .last_mut()
                        .ok_or_else(|| bad("orphan feature continuation"))?
                        .continuation(data, &mut total_qualifiers)?;
                } else {
                    if !token(kind) {
                        return Err(bad("invalid protein feature key"));
                    }
                    builders.push(Builder {
                        kind: kind.into(),
                        location: data.into(),
                        ..Default::default()
                    });
                }
            } else {
                let key = line.get(..12).unwrap_or(line).trim();
                let data = line.get(12..).unwrap_or("").trim();
                if !key.is_empty() {
                    header_field = key.into();
                }
                match key {
                    "ACCESSION" => accessions.extend(data.split_whitespace().map(str::to_string)),
                    "VERSION" => {
                        metadata.insert(
                            "protein:version".into(),
                            data.split_whitespace().next().unwrap_or("").into(),
                        );
                    }
                    "DEFINITION" => {
                        description.push_str(data);
                    }
                    "" if header_field == "ACCESSION" => {
                        accessions.extend(data.split_whitespace().map(str::to_string));
                    }
                    "" if header_field == "DEFINITION" => {
                        description.push(' ');
                        description.push_str(data);
                    }
                    "SOURCE" => {
                        metadata.insert("protein:organism".into(), data.into());
                    }
                    _ => (),
                }
            }
        } else if let Some(ft) = line.strip_prefix("FT   ") {
            if !ft.starts_with(' ') {
                let (kind, rest) = ft
                    .split_once(char::is_whitespace)
                    .ok_or_else(|| bad("FT location missing"))?;
                if !token(kind) {
                    return Err(bad("invalid FT key"));
                }
                let mut words = rest.split_whitespace();
                let a = words.next().ok_or_else(|| bad("FT location missing"))?;
                let b = words.next();
                let legacy = b.is_some_and(|v| {
                    v.trim_start_matches(['<', '>', '?'])
                        .parse::<usize>()
                        .is_ok()
                        || v == "?"
                });
                let (location, note) = if legacy {
                    (
                        format!("{a}..{}", b.expect("checked")),
                        words.collect::<Vec<_>>().join(" "),
                    )
                } else {
                    (
                        a.into(),
                        b.into_iter().chain(words).collect::<Vec<_>>().join(" "),
                    )
                };
                builders.push(Builder {
                    kind: kind.into(),
                    location,
                    legacy_note: note,
                    review: legacy,
                    legacy,
                    ..Default::default()
                });
            } else {
                builders
                    .last_mut()
                    .ok_or_else(|| bad("orphan FT continuation"))?
                    .continuation(ft.trim(), &mut total_qualifiers)?;
            }
        } else if line.starts_with("SQ   ") {
            if sq {
                return Err(bad("duplicate SQ header"));
            }
            sq = true;
            let tokens: Vec<_> = line.split_whitespace().collect();
            if declared_length(&tokens, "AA")? != length {
                return Err(bad("ID/SQ lengths disagree"));
            }
            if let Some(i) = tokens
                .iter()
                .position(|v| v.trim_end_matches(';') == "CRC64")
            {
                let value = i
                    .checked_sub(1)
                    .and_then(|i| tokens.get(i))
                    .ok_or_else(|| bad("missing CRC64"))?;
                if value.len() != 16 || !value.bytes().all(|b| b.is_ascii_hexdigit()) {
                    return Err(bad("invalid CRC64"));
                }
                crc = Some(value.to_ascii_uppercase());
            }
            metadata.insert("protein:source_sq".into(), line.into());
            in_sequence = true;
        } else {
            let key = line.get(..2).unwrap_or("");
            let data = line.get(5..).unwrap_or("").trim();
            match key {
                "AC" => accessions.extend(
                    data.split(';')
                        .map(str::trim)
                        .filter(|v| !v.is_empty())
                        .map(str::to_string),
                ),
                "DE" => {
                    if !description.is_empty() {
                        description.push(' ');
                    }
                    description.push_str(data);
                }
                "OS" => {
                    let value = metadata
                        .entry("protein:organism".into())
                        .or_insert_with(String::new);
                    if !value.is_empty() {
                        value.push(' ');
                    }
                    value.push_str(data);
                }
                _ => (),
            }
        }
        if builders.len() > MAX_FEATURES {
            return Err(Error::TooLarge);
        }
    }
    if !in_sequence || sequence.len() != length {
        return Err(bad("sequence/header length mismatch or missing sequence"));
    }
    let sequence = Sequence::new(Alphabet::Protein, &sequence)?;
    if let Some(crc) = crc {
        if crc != crc64(sequence.as_str().as_bytes()) {
            return Err(bad("UniProt CRC64 mismatch"));
        }
        metadata.insert("protein:source_crc64_verified".into(), "true".into());
    }
    if !accessions.iter().all(|v| token(v)) {
        return Err(bad("invalid accession token"));
    }
    metadata.insert("protein:entry_name".into(), entry.clone());
    metadata.insert("protein:accessions".into(), accessions.join(";"));
    metadata.insert("protein:description".into(), description);
    metadata.insert("protein:authority".into(),"Imported protein annotations are documentary; no physical material or verified cloning history is inferred.".into());
    let features = builders
        .into_iter()
        .map(|b| b.finish(length, profile))
        .collect::<Result<Vec<_>, _>>()?;
    let mut molecule = Molecule {
        name: entry,
        sequence,
        topology: Topology::Linear,
        features,
        metadata,
        original_genbank: None,
    };
    retain(&mut molecule, bytes, profile);
    molecule.validate()?;
    Ok(molecule)
}
pub fn import_bytes(bytes: &[u8], profile: Profile) -> Result<ImportReport, Error> {
    if bytes.is_empty() || bytes.len() > MAX_BYTES {
        return Err(Error::TooLarge);
    }
    let text = std::str::from_utf8(bytes).map_err(|_| bad("UTF-8 text required"))?;
    let mut ranges = Vec::new();
    let mut start = 0;
    let mut offset = 0;
    for line in text.split_inclusive('\n') {
        if line.len() > MAX_LINE {
            return Err(Error::TooLarge);
        }
        offset += line.len();
        if line.trim() == "//" {
            ranges.push(start..offset);
            start = offset;
            if ranges.len() > MAX_RECORDS {
                return Err(Error::TooLarge);
            }
        }
    }
    if !text[start..].trim().is_empty() || ranges.is_empty() {
        return Err(bad("unterminated record"));
    }
    if let Some(last) = ranges.last_mut() {
        last.end = bytes.len();
    }
    let mut report = ImportReport { molecules:Vec::new(),warnings:vec!["Protein annotations imported as documentary data. Exact source records retained; unmodeled headers remain in the original attachment.".into()] };
    for range in ranges {
        if range.len() > MAX_RECORD_BYTES {
            return Err(Error::TooLarge);
        }
        let molecule = parse_record(&bytes[range], profile)?;
        if profile == Profile::UniProt
            && !molecule
                .metadata
                .contains_key("protein:source_crc64_verified")
        {
            report.warnings.push(format!(
                "{}: source SQ has no CRC64; sequence integrity is not checksum-verified.",
                molecule.name
            ));
        }
        if profile == Profile::UniProt
            && molecule
                .metadata
                .get("protein:accessions")
                .is_none_or(|v| v.is_empty())
        {
            report.warnings.push(format!("{}: source record has no accession; normalized UniProt export needs an explicit identifier.",molecule.name));
        }
        if molecule.features.iter().any(|f| f.needs_review) {
            report.warnings.push(format!("{}: legacy, fuzzy, unresolved, or unparsed feature annotations require review; coordinates were not guessed.",molecule.name));
        }
        report.molecules.push(molecule);
        let feature_count: usize = report.molecules.iter().map(|m| m.features.len()).sum();
        let qualifier_count: usize = report
            .molecules
            .iter()
            .flat_map(|m| &m.features)
            .map(|f| f.qualifiers.len())
            .sum();
        if feature_count > MAX_FEATURES || qualifier_count > MAX_QUALIFIERS {
            return Err(Error::TooLarge);
        }
    }
    crate::model::bounds::model_preflight(&report.molecules.iter().collect::<Vec<_>>())?;
    Ok(report)
}
/// UniProt/Swiss-Prot CRC64: reflected ISO polynomial, initial0, final XOR0;
/// distinct from CRC-64/ECMA and CRC-64/GO-ISO's all-ones initialization/XOR.
pub fn crc64(bytes: &[u8]) -> String {
    let mut crc = 0u64;
    for byte in bytes {
        crc ^= u64::from(*byte);
        for _ in 0..8 {
            crc = if crc & 1 != 0 {
                (crc >> 1) ^ 0xd800_0000_0000_0000
            } else {
                crc >> 1
            };
        }
    }
    format!("{crc:016X}")
}

fn text_field(value: &str) -> Result<(), Error> {
    if value.len() > MAX_LINE || value.chars().any(|c| c.is_control()) {
        return Err(unsupported(
            "text field cannot be represented safely in a protein flat file",
        ));
    }
    Ok(())
}
fn location_text(location: &Location, profile: Profile, depth: usize) -> Result<String, Error> {
    if depth > 16 {
        return Err(Error::TooLarge);
    }
    match location {
        Location::Range {
            start,
            end,
            before,
            after,
        } => {
            let a = format!("{}{}", if *before { "<" } else { "" }, start + 1);
            let b = format!("{}{}", if *after { ">" } else { "" }, end);
            Ok(if *end == start + 1 && !before && !after {
                a
            } else {
                format!("{a}..{b}")
            })
        }
        Location::Join(parts) | Location::Order(parts) if profile == Profile::GenPept => {
            let children = parts
                .iter()
                .map(|p| location_text(p, profile, depth + 1))
                .collect::<Result<Vec<_>, _>>()?;
            Ok(format!(
                "{}({})",
                if matches!(location, Location::Join(_)) {
                    "join"
                } else {
                    "order"
                },
                children.join(",")
            ))
        }
        Location::Unresolved(raw) => {
            text_field(raw)?;
            // Exact unknown/isoform source syntax is allowed only under its
            // original profile. The outer export checks that source binding.
            if raw.is_empty() || raw.bytes().any(|b| b.is_ascii_whitespace()) {
                return Err(unsupported("unresolved location cannot be rendered safely"));
            }
            Ok(raw.clone())
        }
        _ => Err(unsupported(
            "this location tree is not representable in the chosen protein profile",
        )),
    }
}
fn uniprot_accession(value: &str) -> bool {
    let b = value.as_bytes();
    if b.len() == 6 && matches!(b[0], b'O' | b'P' | b'Q') {
        return b[1].is_ascii_digit()
            && b[2..5]
                .iter()
                .all(|v| v.is_ascii_uppercase() || v.is_ascii_digit())
            && b[5].is_ascii_digit();
    }
    if !matches!(b.len(), 6 | 10)
        || !b[0].is_ascii_uppercase()
        || matches!(b[0], b'O' | b'P' | b'Q')
        || !b[1].is_ascii_digit()
    {
        return false;
    }
    b[2..].chunks_exact(4).all(|p| {
        p[0].is_ascii_uppercase()
            && p[1..3]
                .iter()
                .all(|v| v.is_ascii_uppercase() || v.is_ascii_digit())
            && p[3].is_ascii_digit()
    })
}
fn lines_wrapped(out: &mut String, prefix: &str, value: &str, width: usize) -> Result<(), Error> {
    text_field(value)?;
    // Break only at existing whitespace: never insert a space into evidence IDs
    // or protein strings. An unbreakable token wider than the profile is refused.
    let mut column = prefix.len();
    out.push_str(prefix);
    for (i, word) in value.split_whitespace().enumerate() {
        if prefix.len() + word.len() > width {
            return Err(unsupported(
                "unbreakable text token exceeds flat-file line width",
            ));
        }
        if i > 0 && column + 1 + word.len() > width {
            out.push('\n');
            out.push_str(if matches!(prefix, "DEFINITION  " | "ACCESSION   ") {
                "            "
            } else {
                prefix
            });
            column = prefix.len();
        } else if column > prefix.len() {
            out.push(' ');
            column += 1;
        }
        out.push_str(word);
        column += word.len();
    }
    out.push('\n');
    Ok(())
}
/// Normalized exchange export, not a byte-identical reconstruction or database
/// submission. Original evidence is available separately via original_bytes().
pub fn export(molecules: &[Molecule], profile: Profile) -> Result<ExportReport, Error> {
    export_with_options(molecules, profile, &ExportOptions::default())
}
pub fn export_with_options(
    molecules: &[Molecule],
    profile: Profile,
    options: &ExportOptions,
) -> Result<ExportReport, Error> {
    if molecules.len() != 1 && (options.entry_name.is_some() || options.uniprot_accession.is_some())
    {
        return Err(unsupported(
            "explicit record identifier overrides require a single protein record",
        ));
    }

    if molecules.is_empty() || molecules.len() > MAX_RECORDS {
        return Err(Error::TooLarge);
    }
    crate::model::bounds::model_preflight(&molecules.iter().collect::<Vec<_>>())?;
    let mut output = String::new();
    let mut warnings = vec!["Normalized protein export does not reproduce all source headers/citations and does not establish database review or physical-material authority; retain the native workspace and original attachment.".into()];
    for molecule in molecules {
        molecule.validate()?;
        crate::format_fidelity::reference_export_guard(molecule)?;
        if molecule.sequence.alphabet() != Alphabet::Protein
            || molecule.topology != Topology::Linear
        {
            return Err(unsupported(
                "annotated protein export requires an explicit linear Protein document",
            ));
        }
        let sequence = molecule.sequence.as_str();
        if sequence.is_empty() || sequence.contains('*') {
            return Err(unsupported(
                "empty/stopped protein sequence is not representable in annotated protein exchange",
            ));
        }
        let entry = options
            .entry_name
            .as_ref()
            .or_else(|| molecule.metadata.get("protein:entry_name"))
            .map(String::as_str)
            .unwrap_or(&molecule.name);
        if !token(entry) || entry.len() > 20 {
            return Err(unsupported("a stable entry identifier is required"));
        }
        let accession_text = options
            .uniprot_accession
            .as_ref()
            .filter(|_| profile == Profile::UniProt)
            .or_else(|| molecule.metadata.get("protein:accessions"))
            .map(String::as_str)
            .unwrap_or("");
        let accessions: Vec<_> = accession_text
            .split(';')
            .filter(|v| !v.is_empty())
            .collect();
        if !accessions.iter().all(|a| token(a)) {
            return Err(unsupported("invalid source accession identifier"));
        }
        let description = molecule
            .metadata
            .get("protein:description")
            .map(String::as_str)
            .filter(|v| !v.is_empty())
            .unwrap_or(&molecule.name);
        text_field(description)?;
        match profile {
            Profile::GenPept => {
                let date = options.locus_date.clone().or_else(||molecule.metadata.get("protein:locus_date").and_then(|v|InsdcDate::new(v).ok())).ok_or_else(||unsupported("GenPept export requires an explicit valid export date or retained source LOCUS date"))?;
                writeln!(
                    output,
                    "LOCUS       {entry:<20} {:>7} aa            linear   UNK {}",
                    sequence.len(),
                    date.as_str()
                )
                .expect("String");
                lines_wrapped(&mut output, "DEFINITION  ", description, 79)?;
                if !accessions.is_empty() {
                    lines_wrapped(&mut output, "ACCESSION   ", &accessions.join(" "), 79)?;
                }
                if let Some(version) = molecule.metadata.get("protein:version") {
                    text_field(version)?;
                    writeln!(output, "VERSION     {version}").expect("String");
                }
                output.push_str("COMMENT     Collagene normalized protein export; source identifiers are documentary.\nFEATURES             Location/Qualifiers\n");
            }
            Profile::UniProt => {
                if accessions.is_empty() || !accessions.iter().all(|a| uniprot_accession(a)) {
                    return Err(unsupported(
                        "UniProt export requires explicit valid UniProt accession metadata; no accession is invented",
                    ));
                }
                writeln!(
                    output,
                    "ID   {entry:<20} Unreviewed; {:>9} AA.",
                    sequence.len()
                )
                .expect("String");
                lines_wrapped(
                    &mut output,
                    "AC   ",
                    &format!("{};", accessions.join("; ")),
                    79,
                )?;
                if description.starts_with("RecName:") || description.starts_with("SubName:") {
                    lines_wrapped(&mut output, "DE   ", description, 79)?;
                } else {
                    lines_wrapped(
                        &mut output,
                        "DE   ",
                        &format!("SubName: Full={};", description.trim_end_matches(';')),
                        79,
                    )?;
                }
                output.push_str("CC   -!- CAUTION: Collagene generated exchange record; not a database release.\n");
            }
        }
        for feature in &molecule.features {
            if !token(&feature.kind) || feature.kind.len() > 15 {
                return Err(unsupported(
                    "protein feature key exceeds the selected layout",
                ));
            }
            if matches!(feature.location, Location::Unresolved(_))
                && molecule.metadata.get("protein:profile").map(String::as_str)
                    != Some(profile.name())
            {
                return Err(unsupported(
                    "unresolved source coordinates cannot be translated between profiles",
                ));
            }
            if profile == Profile::UniProt && !known_key(profile, &feature.kind) {
                return Err(unsupported("feature key has no declared UniProt mapping"));
            }
            if profile == Profile::UniProt {
                if let Location::Unresolved(raw) = &feature.location {
                    let (isoform, location) = raw
                        .rsplit_once(':')
                        .map_or((None, raw.as_str()), |(id, v)| (Some(id), v));
                    if let Some(id) = isoform {
                        let valid = id.rsplit_once('-').is_some_and(|(a, n)| {
                            uniprot_accession(a) && n.parse::<usize>().is_ok_and(|v| v > 0)
                        });
                        if !valid {
                            return Err(unsupported(
                                "unresolved isoform identifier is not representable by current UniProt FT grammar",
                            ));
                        }
                    }
                    let parts: Vec<_> = location.split("..").collect();
                    let valid = parts.len() <= 2
                        && parts.iter().enumerate().all(|(index, value)| {
                            if *value == "?" {
                                return true;
                            }
                            let value = value
                                .strip_prefix(if index == 0 { '<' } else { '>' })
                                .or_else(|| value.strip_prefix('?'))
                                .unwrap_or(value);
                            value
                                .parse::<usize>()
                                .is_ok_and(|n| n > 0 && (isoform.is_some() || n <= sequence.len()))
                        });
                    if !valid {
                        return Err(unsupported(
                            "unparsed location is not representable by current UniProt FT grammar",
                        ));
                    }
                }
            }
            let location = if profile == Profile::UniProt && endpoint_feature(&feature.kind) {
                match &feature.location {
                    Location::Join(parts) if parts.len() == 2 => match (&parts[0], &parts[1]) {
                        (
                            Location::Range {
                                start: a,
                                end: b,
                                before: false,
                                after: false,
                            },
                            Location::Range {
                                start: c,
                                end: d,
                                before: false,
                                after: false,
                            },
                        ) if *b == a + 1 && *d == c + 1 && a < c => format!("{}..{}", a + 1, c + 1),
                        _ => {
                            return Err(unsupported(
                                "crosslink requires two exact endpoint residues",
                            ));
                        }
                    },
                    Location::Unresolved(raw) => {
                        location_text(&Location::Unresolved(raw.clone()), profile, 0)?
                    }
                    _ => {
                        return Err(unsupported(
                            "crosslink interval cannot stand in for two participating residues",
                        ));
                    }
                }
            } else {
                location_text(&feature.location, profile, 0)?
            };
            if location.len() > 58 {
                return Err(unsupported(
                    "location expression exceeds flat-file line width",
                ));
            }
            match profile {
                Profile::GenPept => {
                    writeln!(output, "     {:<16}{location}", feature.kind).expect("String")
                }
                Profile::UniProt => {
                    writeln!(output, "FT   {:<15}{location}", feature.kind).expect("String")
                }
            }
            for (key, value) in &feature.qualifiers {
                if key.starts_with("protein:") {
                    if key == "protein:unparsed_qualifier" {
                        return Err(unsupported(
                            "unparsed source qualifier needs original/native export",
                        ));
                    }
                    continue;
                }
                if !token(key) {
                    return Err(unsupported("qualifier name cannot be represented"));
                }
                let key = if profile == Profile::UniProt && key == "FTId" {
                    "id"
                } else {
                    key
                };
                let field = match value {
                    None => format!("/{key}"),
                    Some(value) => {
                        text_field(value)?;
                        format!("/{key}=\"{}\"", value.replace('"', "\"\""))
                    }
                };
                lines_wrapped(
                    &mut output,
                    if profile == Profile::GenPept {
                        "                     "
                    } else {
                        "FT                   "
                    },
                    &field,
                    79,
                )?;
            }
            if feature.needs_review {
                warnings.push(format!("{}: feature {} still requires review; exported coordinates remain documentary.",molecule.name,feature.kind));
            }
        }
        if profile == Profile::GenPept {
            output.push_str("ORIGIN\n");
        } else {
            let options = crate::protein_properties::ProteinOptions {
                stop_policy: crate::protein_properties::StopPolicy::Reject,
                ..Default::default()
            };
            // SQ defines the unmodified primary-chain MW, not the mature/PTM
            // material's mass. Analyze an explicit sequence-only projection;
            // retain the original annotations and their native chemistry guard.
            let primary_chain = Molecule {
                name: "UniProt unmodified primary-chain mass projection".into(),
                sequence: molecule.sequence.clone(),
                topology: Topology::Linear,
                features: vec![],
                metadata: BTreeMap::new(),
                original_genbank: None,
            };
            let properties = crate::protein_properties::analyze_protein(&primary_chain, &options)?;
            let mass = properties
                .mass
                .filter(|m| m.exact_identity && m.average_da.minimum == m.average_da.maximum)
                .ok_or_else(|| {
                    unsupported(
                        "UniProt molecular mass is unavailable for unresolved residue identities",
                    )
                })?;
            writeln!(
                output,
                "SQ   SEQUENCE   {} AA;  {:.0} MW;  {} CRC64;",
                sequence.len(),
                mass.average_da.minimum,
                crc64(sequence.as_bytes())
            )
            .expect("String");
            warnings.push("UniProt SQ MW is computed only from an explicit unmodified reduced primary-sequence projection under the native protein mass model. Retained modification/crosslink annotations are not quantitative chemical-state authority; processed/PTM material masses remain unavailable without a supported model.".into());
        }
        for (i, line) in sequence.as_bytes().chunks(60).enumerate() {
            if profile == Profile::GenPept {
                write!(output, "{:>9} ", i * 60 + 1).expect("String");
            } else {
                output.push_str("     ");
            }
            for (j, group) in line.chunks(10).enumerate() {
                if j > 0 {
                    output.push(' ');
                }
                let group = std::str::from_utf8(group).expect("protein ASCII");
                if profile == Profile::GenPept {
                    output.push_str(&group.to_ascii_lowercase());
                } else {
                    output.push_str(group);
                }
            }
            output.push('\n');
        }
        output.push_str("//\n");
        if output.len() > MAX_BYTES {
            return Err(Error::TooLarge);
        }
    }
    Ok(ExportReport {
        text: output,
        warnings,
    })
}
