//! EMBL/ENA nucleotide flat-file framing with a shared INSDC feature grammar.
//! This adapter is independently implemented; gb-io parses the common location/qualifier syntax.
use crate::{
    Alphabet, Error, Molecule, Sequence, Topology,
    formats::{self, ImportReport, MAX_IMPORT_BYTES},
    model::MAX_SEQUENCE_LEN,
};
use std::collections::BTreeSet;

const ORIGINAL: &str = "embl.original_record";
const MAX_RECORDS: usize = 1000;
const MAX_LINE_BYTES: usize = 1024 * 1024;

fn bad(message: impl Into<String>) -> Error {
    Error::Format(format!("EMBL: {}", message.into()))
}
fn unsupported(message: impl Into<String>) -> Error {
    Error::Unsupported(format!("EMBL: {}", message.into()))
}

struct Id {
    name: String,
    length: usize,
    alphabet: Alphabet,
    topology: Topology,
    molecule_type: String,
    version: String,
    data_class: String,
    division: String,
    legacy: bool,
}
fn length_field(value: &str) -> Result<usize, Error> {
    let value = value.trim().trim_end_matches('.');
    let fields: Vec<_> = value.split_ascii_whitespace().collect();
    if fields.len() != 2 || fields[1] != "BP" {
        return Err(unsupported(
            "ID/SQ length must be a nucleotide count followed by BP",
        ));
    }
    let length = fields[0]
        .parse::<usize>()
        .map_err(|_| bad("invalid declared sequence length"))?;
    if length == 0 || length > MAX_SEQUENCE_LEN {
        return Err(Error::TooLarge);
    }
    Ok(length)
}
fn molecule_alphabet(value: &str) -> Option<Alphabet> {
    let value = value.trim().to_ascii_lowercase();
    let value = value
        .strip_prefix("circular ")
        .or_else(|| value.strip_prefix("linear "))
        .unwrap_or(&value);
    match value {
        "dna" | "genomic dna" | "other dna" | "unassigned dna" | "ss-dna" | "ds-dna" => {
            Some(Alphabet::Dna)
        }
        "rna" | "genomic rna" | "mrna" | "rrna" | "trna" | "snrna" | "snorna" | "scrna"
        | "pre-rna" | "transcribed rna" | "viral crna" | "other rna" | "unassigned rna"
        | "crna" | "ss-rna" | "ds-rna" => Some(Alphabet::Rna),
        _ => None,
    }
}
fn parse_id(value: &str) -> Result<Id, Error> {
    let fields: Vec<_> = value.split(';').map(str::trim).collect();
    let modern = fields.len() == 7;
    if !modern && fields.len() != 4 {
        return Err(unsupported(
            "expected a seven-field modern ID or four-field legacy nucleotide ID",
        ));
    }
    let name = if modern {
        fields[0]
    } else {
        fields[0].split_ascii_whitespace().next().unwrap_or("")
    };
    if name.is_empty()
        || name.len() > 1024
        || !name.is_ascii()
        || name
            .bytes()
            .any(|b| b.is_ascii_whitespace() || b.is_ascii_control() || b == b';')
    {
        return Err(bad("invalid entry identifier"));
    }
    let molecule_type = if modern { fields[3] } else { fields[1] };
    let upper = molecule_type.to_ascii_uppercase();
    let alphabet=molecule_alphabet(molecule_type).ok_or_else(||unsupported("protein/Swiss-Prot and unspecified molecule types are not supported nucleotide EMBL records"))?;
    let topology = if modern {
        match fields[2] {
            "linear" => Topology::Linear,
            "circular" => Topology::Circular,
            _ => return Err(unsupported("unknown or missing topology")),
        }
    } else if upper
        .split_ascii_whitespace()
        .any(|word| word == "CIRCULAR")
    {
        Topology::Circular
    } else {
        Topology::Linear
    };
    let version = if modern {
        let version = fields[1].strip_prefix("SV ").ok_or_else(|| {
            unsupported("modern ID sequence version must be SV followed by a positive integer")
        })?;
        if version.parse::<u64>().ok().filter(|n| *n > 0).is_none() {
            return Err(bad("invalid ID sequence version"));
        }
        version.to_owned()
    } else {
        "1".into()
    };
    let data_class = if modern { fields[4] } else { "STD" };
    let division = if modern { fields[5] } else { fields[2] };
    for token in [data_class, division] {
        if token.is_empty()
            || token.len() > 16
            || !token
                .bytes()
                .all(|b| b.is_ascii_alphanumeric() || b == b'_')
        {
            return Err(bad("invalid ID class or division"));
        }
    }
    Ok(Id {
        name: name.into(),
        length: length_field(fields.last().unwrap())?,
        alphabet,
        topology,
        molecule_type: molecule_type.into(),
        version,
        data_class: data_class.into(),
        division: division.into(),
        legacy: !modern,
    })
}
fn tag(line: &str) -> Result<(&str, &str), Error> {
    let code = line.get(..2).ok_or_else(|| bad("short tagged line"))?;
    if !code.bytes().all(|b| b.is_ascii_uppercase()) {
        return Err(bad(
            "header line must start with a two-letter uppercase field code",
        ));
    }
    if line.len() == 2 {
        return Ok((code, ""));
    }
    let separator = line
        .get(2..5)
        .ok_or_else(|| bad("header field lacks three separating spaces"))?;
    if separator != "   " {
        return Err(bad(
            "header fields require three spaces after the field code",
        ));
    }
    Ok((
        code,
        line.get(5..).ok_or_else(|| bad("invalid field encoding"))?,
    ))
}
fn sequence_line(line: &str, bases: &mut String) -> Result<(), Error> {
    if !line.starts_with("     ") {
        return Err(bad("sequence data must start with five spaces"));
    }
    let mut words = line.split_ascii_whitespace().collect::<Vec<_>>();
    let count = words
        .pop()
        .ok_or_else(|| bad("empty sequence line"))?
        .parse::<usize>()
        .map_err(|_| bad("sequence line needs a trailing cumulative base count"))?;
    if words.is_empty() {
        return Err(bad("sequence line has no bases"));
    }
    for word in words {
        if !word.bytes().all(|b| b.is_ascii_alphabetic()) {
            return Err(bad("invalid character in sequence block"));
        }
        if bases.len().saturating_add(word.len()) > MAX_SEQUENCE_LEN {
            return Err(Error::TooLarge);
        }
        bases.push_str(word);
    }
    if count != bases.len() {
        return Err(bad(
            "sequence-line cumulative count differs from parsed bases",
        ));
    }
    Ok(())
}
fn known_header(code: &str) -> bool {
    matches!(
        code,
        "AC" | "SV"
            | "PR"
            | "DT"
            | "DE"
            | "KW"
            | "OS"
            | "OC"
            | "OG"
            | "OX"
            | "RN"
            | "RP"
            | "RC"
            | "RA"
            | "RT"
            | "RL"
            | "RX"
            | "DR"
            | "CC"
            | "AH"
            | "AS"
            | "FH"
            | "FT"
            | "XX"
    )
}
fn safe_after_sequence_change(code: &str) -> bool {
    // Retained references, comments, assembly maps and unknown headers can contain coordinates.
    matches!(
        code,
        "ID" | "AC"
            | "SV"
            | "PR"
            | "DT"
            | "DE"
            | "KW"
            | "OS"
            | "OC"
            | "OG"
            | "OX"
            | "DR"
            | "FH"
            | "FT"
            | "XX"
            | "SQ"
    )
}
fn sq_length(value: &str) -> Result<Option<usize>, Error> {
    if value.is_empty() {
        return Ok(None);
    }
    let length = value
        .strip_prefix("Sequence ")
        .ok_or_else(|| bad("unsupported SQ header"))?
        .split(';')
        .next()
        .unwrap();
    Ok(Some(length_field(length)?))
}

fn validate_composition(value: &str, bases: &str) -> Result<(), Error> {
    let mut expected = std::collections::BTreeMap::new();
    for item in value
        .split(';')
        .skip(1)
        .map(str::trim)
        .filter(|s| !s.is_empty())
    {
        let words: Vec<_> = item.split_ascii_whitespace().collect();
        if words.len() != 2 || !matches!(words[1], "A" | "C" | "G" | "T" | "U" | "other") {
            return Err(unsupported("unrecognized SQ composition field"));
        }
        let count = words[0]
            .parse::<usize>()
            .map_err(|_| bad("invalid SQ composition count"))?;
        if expected.insert(words[1], count).is_some() {
            return Err(bad("duplicate SQ composition count"));
        }
    }
    if !expected.is_empty() {
        let upper = bases.to_ascii_uppercase();
        let mut total = 0usize;
        for (label, count) in expected {
            let actual = if label == "other" {
                upper.bytes().filter(|b| !b"ACGTU".contains(b)).count()
            } else {
                upper.bytes().filter(|b| *b == label.as_bytes()[0]).count()
            };
            if count != actual {
                return Err(bad(format!(
                    "SQ {label} composition count differs from sequence"
                )));
            }
            total = total.checked_add(count).ok_or(Error::TooLarge)?;
        }
        if total != bases.len() {
            return Err(bad("SQ composition counts do not cover the sequence"));
        }
    }
    Ok(())
}

fn parse_record(record: &str) -> Result<(Molecule, Vec<String>), Error> {
    let mut lines = record.lines();
    let first = lines
        .next()
        .ok_or_else(|| bad("empty record"))?
        .trim_end_matches('\r');
    let (code, value) = tag(first)?;
    if code != "ID" {
        return Err(bad("record must begin with ID"));
    }
    let mut id = parse_id(value)?;
    let mut version_header_seen = false;
    let mut features = String::new();
    let mut definition = Vec::new();
    let mut bases = String::new();
    let mut in_sequence = false;
    let mut ended = false;
    let mut declared_sq = None;
    let mut composition = String::new();
    let mut warnings = Vec::new();
    let mut unknown = BTreeSet::new();
    for line in lines {
        let line = line.trim_end_matches('\r');
        if line == "//" {
            ended = true;
            break;
        }
        if line.is_empty() {
            return Err(bad(
                "blank lines inside records are not EMBL field separators; use XX",
            ));
        }
        if in_sequence {
            sequence_line(line, &mut bases)?;
            continue;
        }
        let (code, value) = tag(line)?;
        match code {
            "ID" => return Err(bad("new ID before record terminator")),
            "SV" => {
                if version_header_seen {
                    return Err(bad("duplicate SV header"));
                }
                version_header_seen = true;
                let version = value
                    .rsplit_once('.')
                    .map(|(_, v)| v)
                    .ok_or_else(|| bad("SV must contain an accession.version identifier"))?;
                if version.parse::<u64>().ok().filter(|n| *n > 0).is_none() {
                    return Err(bad("invalid SV sequence version"));
                }
                if !id.legacy && version != id.version {
                    return Err(bad("ID and SV sequence versions disagree"));
                }
                id.version = version.into();
            }
            "SQ" => {
                declared_sq = sq_length(value)?;
                composition = value.to_owned();
                in_sequence = true;
            }
            "CO" => {
                return Err(unsupported(
                    "CO assembly/reference-resolution records are not implemented",
                ));
            }
            "FT" => {
                if line.len() < 21 {
                    return Err(bad(
                        "feature line is shorter than the INSDC feature columns",
                    ));
                }
                // EMBL FT and GenBank feature keys/qualifiers share the same absolute columns.
                features.push_str("  ");
                features.push_str(&line[2..]);
                features.push('\n');
            }
            "DE" => definition.push(value.trim().to_owned()),
            _ => {
                if !known_header(code) {
                    unknown.insert(code.to_owned());
                }
            }
        }
    }
    if !ended {
        return Err(bad("record lacks // terminator"));
    }
    if !in_sequence {
        return Err(unsupported(
            "record has no SQ sequence; external reference resolution is not implemented",
        ));
    }
    if bases.len() != id.length || declared_sq.is_some_and(|n| n != bases.len()) {
        return Err(bad("ID/SQ declared length differs from sequence length"));
    }
    validate_composition(&composition, &bases)?;
    // The application stores RNA as U. This explicit normalization is reversible on export;
    // original evidence retains the actual T/U encoding supplied by the file.
    if id.alphabet == Alphabet::Rna && bases.bytes().any(|b| b == b't' || b == b'T') {
        bases = bases
            .chars()
            .map(|c| {
                if c == 't' {
                    'u'
                } else if c == 'T' {
                    'U'
                } else {
                    c
                }
            })
            .collect();
        warnings.push(format!("{}: RNA-declared sequence normalized from T to U; original EMBL text retains the supplied encoding.",id.name));
    }
    let sequence = Sequence::new(id.alphabet, &bases)?;
    let topology = if id.topology == Topology::Circular {
        "circular"
    } else {
        "linear"
    };
    let alphabet = if id.alphabet == Alphabet::Rna {
        "RNA"
    } else {
        "DNA"
    };
    let genbank = format!(
        "LOCUS       embl_adapter              {} bp    {}     {} SYN 01-JAN-2026\nFEATURES             Location/Qualifiers\n{}ORIGIN\n        1 {}\n//\n",
        sequence.len(),
        alphabet,
        topology,
        features,
        sequence.as_str()
    );
    let mut parsed =
        formats::import_genbank(&genbank).map_err(|e| bad(format!("feature table: {e}")))?;
    // Discard the adapter's generated syntax, never mislabel it as original evidence.
    let mut molecule = parsed.molecules.remove(0);
    molecule.name = id.name;
    molecule.sequence = sequence;
    molecule.original_genbank = None;
    molecule.metadata.remove("source_record_index");
    molecule.metadata.insert(ORIGINAL.into(), record.to_owned());
    if formats::retain_digest_producer_comments(&mut molecule, record) {
        warnings.push("Sequence-only digest producer comments were retained as unverified source data; strand/end chemistry was not restored.".into());
    }
    let molecule_type = if id.legacy {
        id.molecule_type
            .split_ascii_whitespace()
            .filter(|word| {
                !word.eq_ignore_ascii_case("circular") && !word.eq_ignore_ascii_case("linear")
            })
            .collect::<Vec<_>>()
            .join(" ")
    } else {
        id.molecule_type
    };
    molecule
        .metadata
        .insert("embl.molecule_type".into(), molecule_type);
    molecule
        .metadata
        .insert("embl.sequence_version".into(), id.version);
    molecule
        .metadata
        .insert("embl.data_class".into(), id.data_class);
    molecule
        .metadata
        .insert("embl.division".into(), id.division);
    if !definition.is_empty() {
        molecule
            .metadata
            .insert("definition".into(), definition.join(" "));
    }
    for warning in parsed.warnings {
        if !warning.starts_with("Original GenBank text") {
            warnings.push(warning.replace("embl_adapter", &molecule.name));
        }
    }
    if id.legacy {
        warnings.push(format!(
            "{}: legacy four-field ID imported; {}. Regenerated IDs use modern syntax.",
            molecule.name,
            if id.topology == Topology::Circular {
                "explicit circular topology retained"
            } else {
                "omitted topology defaults to linear"
            }
        ));
    }
    if !unknown.is_empty() {
        warnings.push(format!(
            "{}: unrecognized EMBL headers {} retained as original text but not interpreted.",
            molecule.name,
            unknown.into_iter().collect::<Vec<_>>().join(", ")
        ));
    }
    molecule.validate()?;
    Ok((molecule, warnings))
}

/// Read bounded, terminated multi-record EMBL/ENA nucleotide text. No remote references are fetched.
pub fn import_embl(text: &str) -> Result<ImportReport, Error> {
    if text.len() > MAX_IMPORT_BYTES {
        return Err(Error::TooLarge);
    }
    let text = text.trim_start_matches('\u{feff}');
    let mut molecules = Vec::new();
    let mut warnings = Vec::new();
    let mut start = None;
    let mut offset = 0usize;
    for line in text.split_inclusive('\n') {
        if line.len() > MAX_LINE_BYTES {
            return Err(Error::TooLarge);
        }
        let body = line.trim_end_matches(['\r', '\n']);
        if start.is_none() {
            if body.is_empty() {
                offset += line.len();
                continue;
            }
            if !body.starts_with("ID   ") {
                return Err(bad("expected ID at the start of each record"));
            }
            start = Some(offset);
        }
        offset += line.len();
        if body == "//" {
            if molecules.len() >= MAX_RECORDS {
                return Err(Error::TooLarge);
            }
            let record = &text[start.take().unwrap()..offset];
            let (molecule, mut record_warnings) = parse_record(record)?;
            molecules.push(molecule);
            warnings.append(&mut record_warnings);
        }
    }
    if start.is_some() {
        return Err(bad("record lacks // terminator"));
    }
    if molecules.is_empty() {
        return Err(bad("no records"));
    }
    warnings.push("Original EMBL records are retained individually in native metadata. Common INSDC features/qualifiers are interpreted; retained reference/header text is not an editable submission model.".into());
    Ok(ImportReport {
        molecules,
        warnings,
    })
}

fn one_line(value: &str) -> String {
    value.replace(['\r', '\n'], " ")
}
fn identifier(value: &str) -> Result<String, Error> {
    let value = value.replace(char::is_whitespace, "_");
    if value.is_empty()
        || value.len() > 1024
        || !value.is_ascii()
        || value.bytes().any(|b| b == b';' || b.is_ascii_control())
    {
        return Err(unsupported(
            "entry name cannot be represented as an EMBL identifier",
        ));
    }
    Ok(value)
}
fn metadata_token(molecule: &Molecule, key: &str, default: &str) -> Result<String, Error> {
    let value = molecule
        .metadata
        .get(key)
        .map(String::as_str)
        .unwrap_or(default);
    if value.is_empty()
        || value.len() > 16
        || !value
            .bytes()
            .all(|b| b.is_ascii_alphanumeric() || b == b'_')
    {
        return Err(bad(format!("invalid export metadata token {key}")));
    }
    Ok(value.into())
}

/// Export the nucleotide model plus retained header evidence. This is not an ENA submission validator.
/// A sequence edit with reference/comment/assembly/unknown headers is conservatively rejected.
pub fn export_embl(molecule: &Molecule) -> Result<String, Error> {
    molecule.validate()?;
    if molecule.sequence.alphabet() == Alphabet::Protein {
        return Err(unsupported("protein/Swiss-Prot export is not implemented"));
    }
    let mut headers = Vec::new();
    let mut sequence_changed = false;
    if let Some(original) = molecule.metadata.get(ORIGINAL) {
        if original.len() > MAX_IMPORT_BYTES {
            return Err(Error::TooLarge);
        }
        let parsed = import_embl(original)?;
        if parsed.molecules.len() != 1 {
            return Err(bad(
                "retained original evidence must contain exactly one record",
            ));
        }
        let old = &parsed.molecules[0];
        if old == molecule {
            return Ok(original.clone());
        }
        sequence_changed = old.sequence != molecule.sequence;
        for line in original.lines() {
            let line = line.trim_end_matches('\r');
            if line.starts_with("SQ") {
                break;
            }
            let (code, _) = tag(line)?;
            if sequence_changed && !safe_after_sequence_change(code) {
                return Err(unsupported(format!(
                    "retained {code} metadata cannot be remapped after a sequence edit; use the native format to keep the evidence"
                )));
            }
            if !matches!(code, "ID" | "SV" | "DE" | "FH" | "FT" | "XX") {
                headers.push(line.to_owned());
            }
        }
    }
    let id = identifier(&molecule.name)?;
    let topology = if molecule.topology == Topology::Circular {
        "circular"
    } else {
        "linear"
    };
    let mol_type = match molecule.sequence.alphabet() {
        Alphabet::Dna => "genomic DNA",
        Alphabet::Rna => "RNA",
        Alphabet::Protein => unreachable!(),
    };
    let mol_type = if let Some(retained) = molecule.metadata.get("embl.molecule_type") {
        if molecule_alphabet(retained) != Some(molecule.sequence.alphabet())
            || !retained.is_ascii()
            || retained.contains(['\r', '\n', ';'])
            || retained.len() > 64
        {
            return Err(bad(
                "retained molecule type conflicts with the current alphabet",
            ));
        }
        retained.as_str()
    } else {
        mol_type
    };
    let class = metadata_token(molecule, "embl.data_class", "STD")?;
    let division = metadata_token(molecule, "embl.division", "SYN")?;
    let version = metadata_token(molecule, "embl.sequence_version", "1")?;
    if version.parse::<u64>().ok().filter(|v| *v > 0).is_none() {
        return Err(bad("export sequence version must be a positive integer"));
    }
    let mut output = format!(
        "ID   {id}; SV {version}; {topology}; {mol_type}; {class}; {division}; {} BP.\nXX\n",
        molecule.sequence.len()
    );
    for line in headers {
        output.push_str(&line);
        output.push('\n');
    }
    if let Some(definition) = molecule.metadata.get("definition") {
        output.push_str("DE   ");
        output.push_str(&one_line(definition));
        output.push('\n');
    }
    if sequence_changed {
        output.push_str("CC   Collagene: sequence differs from original EMBL evidence; retained\nCC   identifiers are provenance, not a new database submission.\n");
    }
    for comment in crate::formats::provenance_comments(molecule) {
        output.push_str("CC   ");
        output.push_str(&comment);
        output.push('\n');
    }
    output.push_str("XX\nFH   Key             Location/Qualifiers\nFH\n");
    let mut feature_model = molecule.clone();
    feature_model.original_genbank = None;
    // Exporting features uses the same standalone INSDC writer; generated headers are ignored.
    let genbank = formats::export_genbank(&feature_model)?;
    let mut in_features = false;
    for line in genbank.lines() {
        if line.starts_with("FEATURES") {
            in_features = true;
            continue;
        }
        if line.starts_with("ORIGIN") {
            break;
        }
        if in_features {
            output.push_str("FT");
            output.push_str(
                line.get(2..)
                    .ok_or_else(|| bad("feature writer produced a short line"))?,
            );
            output.push('\n');
        }
    }
    let bases = molecule
        .sequence
        .as_str()
        .replace('U', "T")
        .to_ascii_lowercase();
    let counts = |letter: u8| bases.bytes().filter(|b| *b == letter).count();
    let a = counts(b'a');
    let c = counts(b'c');
    let g = counts(b'g');
    let t = counts(b't');
    let other = bases.len() - a - c - g - t;
    output.push_str(&format!(
        "XX\nSQ   Sequence {} BP; {a} A; {c} C; {g} G; {t} T; {other} other;\n",
        bases.len()
    ));
    let mut position = 0;
    for chunk in bases.as_bytes().chunks(60) {
        output.push_str("     ");
        for group in chunk.chunks(10) {
            output.push_str(std::str::from_utf8(group).expect("ASCII sequence"));
            output.push(' ');
        }
        position += chunk.len();
        output.push_str(&format!("{position:>9}\n"));
    }
    output.push_str("//\n");
    if output.len() > MAX_IMPORT_BYTES {
        return Err(Error::TooLarge);
    }
    Ok(output)
}
