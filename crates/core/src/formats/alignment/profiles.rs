//! Independently implemented bounded molecular-alignment exchange profiles.
//! No alignment engine, trees, filesystem, network, or external program is used.
use crate::formats::alignment::{
    self as base, AlignmentAnnotation, AlignmentDocument, AlignmentExport, AlignmentFormat,
    AlignmentMetadata, AlignmentRow, MAX_ALIGNMENT_BYTES, MAX_ALIGNMENT_CELLS, MAX_ALIGNMENT_LINE,
    MAX_ALIGNMENT_ROWS,
};
use collagene_model::{Alphabet, Error, Sequence};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};

type Rows = Vec<(String, String, String)>;
const MAX_TOKENS: usize = 300_000;
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct ExportOptions {
    /// Explicit opt-in. Colliding truncated names remain an error.
    pub truncate_phylip_identifiers: bool,
    pub block_width: usize,
}
impl Default for ExportOptions {
    fn default() -> Self {
        Self {
            truncate_phylip_identifiers: false,
            block_width: 80,
        }
    }
}
fn bad(s: impl Into<String>) -> Error {
    Error::Format(s.into())
}
fn note(meta: &mut AlignmentMetadata, s: impl Into<String>) -> Result<(), Error> {
    if meta.annotations.len() >= base::MAX_ALIGNMENT_ANNOTATIONS {
        return Err(Error::TooLarge);
    }
    meta.annotations.push(AlignmentAnnotation {
        scope: "GF".into(),
        row_identifier: None,
        tag: "CG2_EXCHANGE_NOTE".into(),
        value: s.into().replace('\t', " "),
    });
    Ok(())
}
fn text(bytes: &[u8]) -> Result<&str, Error> {
    if bytes.len() > MAX_ALIGNMENT_BYTES {
        return Err(Error::TooLarge);
    }
    let s =
        std::str::from_utf8(bytes).map_err(|_| bad("alignment exchange requires ASCII UTF-8"))?;
    if !s.is_ascii()
        || s.lines().any(|l| l.len() > MAX_ALIGNMENT_LINE)
        || s.bytes()
            .any(|b| b.is_ascii_control() && !b"\r\n\t".contains(&b))
    {
        return Err(bad(
            "invalid control, non-ASCII input or oversized alignment line",
        ));
    }
    Ok(s)
}
fn dimensions(n: usize, w: usize) -> Result<(), Error> {
    if n == 0
        || n > MAX_ALIGNMENT_ROWS
        || w == 0
        || n.checked_mul(w)
            .is_none_or(|cells| cells > MAX_ALIGNMENT_CELLS)
    {
        return Err(Error::TooLarge);
    }
    Ok(())
}
fn number(s: &str) -> Result<usize, Error> {
    if s.is_empty() || !s.bytes().all(|b| b.is_ascii_digit()) {
        return Err(bad("expected an unsigned alignment dimension"));
    }
    s.parse().map_err(|_| Error::TooLarge)
}
fn label(s: &str) -> Result<(), Error> {
    if s.trim().is_empty() || s.len() > 4096 || !s.is_ascii() || s.chars().any(char::is_control) {
        return Err(bad("invalid alignment identifier"));
    }
    Ok(())
}
fn molecular(s: &str, alphabet: Alphabet) -> Result<String, Error> {
    let s = s.to_ascii_uppercase();
    let residues: String = s.chars().filter(|c| *c != '-').collect();
    if !residues.is_empty() {
        Sequence::new(alphabet, &residues)?;
    }
    Ok(s)
}
fn phylip(format: AlignmentFormat) -> Option<(bool, bool)> {
    match format {
        AlignmentFormat::PhylipSequential => Some((false, false)),
        AlignmentFormat::PhylipInterleaved => Some((false, true)),
        AlignmentFormat::PhylipRelaxedSequential => Some((true, false)),
        AlignmentFormat::PhylipRelaxedInterleaved => Some((true, true)),
        _ => None,
    }
}
/// Same provenance/identity/coordinate document as existing alignment formats.
/// Shared import_parsed validation is deliberately used rather than a second
/// imported-record validation implementation.
pub fn import_alignment(
    bytes: &[u8],
    format: AlignmentFormat,
    alphabet: Alphabet,
    name: &str,
) -> Result<AlignmentDocument, Error> {
    if matches!(
        format,
        AlignmentFormat::Pir | AlignmentFormat::Gde | AlignmentFormat::Selex
    ) {
        return crate::formats::alignment::pir_gde_selex::import_alignment(
            bytes, format, alphabet, name,
        );
    }

    let s = text(bytes)?;
    let (rows, meta) = if let Some((relaxed, interleaved)) = phylip(format) {
        parse_phylip(s, alphabet, relaxed, interleaved)?
    } else {
        match format {
            AlignmentFormat::Msf => parse_msf(s, alphabet)?,
            AlignmentFormat::Nexus => parse_nexus(s, alphabet)?,
            _ => return Err(bad("not an alignment-exchange profile")),
        }
    };
    base::import_parsed(bytes, format, alphabet, name, rows, meta)
}
fn phylip_named(line: &str, relaxed: bool) -> Result<(String, &str), Error> {
    if relaxed {
        let line = line.trim_start();
        let at = line
            .find(char::is_whitespace)
            .ok_or_else(|| bad("relaxed PHYLIP needs a whitespace-delimited identifier"))?;
        let id = &line[..at];
        label(id)?;
        Ok((id.into(), &line[at..]))
    } else {
        if line.len() < 10 {
            return Err(bad(
                "strict PHYLIP requires the complete ten-byte name field",
            ));
        }
        let id = line[..10].trim();
        label(id)?;
        Ok((id.into(), &line[10..]))
    }
}
fn phylip_segment(s: &str, alphabet: Alphabet) -> Result<String, Error> {
    let joined: String = s
        .chars()
        .filter(|c| !c.is_ascii_whitespace() && !c.is_ascii_digit())
        .collect();
    if joined.is_empty() {
        return Err(bad("empty PHYLIP sequence segment"));
    }
    if alphabet == Alphabet::Protein && joined.bytes().any(|b| b"JOUjou".contains(&b)) {
        return Err(bad(
            "PHYLIP protein J/O/U extensions are outside this profile",
        ));
    }
    if joined.contains('.') || joined.contains('?') || joined.contains('O') || joined.contains('o')
    {
        return Err(bad(
            "PHYLIP dot/missing/deletion extensions outside the explicit molecular profile",
        ));
    }
    molecular(&joined, alphabet)
}
fn append(row: &mut String, part: &str, width: usize) -> Result<(), Error> {
    if part.len() > width.saturating_sub(row.len()) {
        return Err(bad("alignment row exceeds its declared width"));
    }
    row.push_str(part);
    Ok(())
}
fn parse_phylip(
    s: &str,
    alphabet: Alphabet,
    relaxed: bool,
    interleaved: bool,
) -> Result<(Rows, AlignmentMetadata), Error> {
    let mut lines = s.lines().filter(|l| !l.trim().is_empty());
    let header = lines
        .next()
        .ok_or_else(|| bad("missing PHYLIP dimensions"))?;
    let fields: Vec<_> = header.split_whitespace().collect();
    if fields.len() != 2 {
        return Err(bad(
            "PHYLIP header must contain exactly NTAX NCHAR; concatenated data sets/options are unsupported",
        ));
    }
    let n = number(fields[0])?;
    let width = number(fields[1])?;
    dimensions(n, width)?;
    let mut rows = Rows::with_capacity(n);
    let mut ids = BTreeSet::new();
    if !interleaved {
        for _ in 0..n {
            let line = lines.next().ok_or_else(|| bad("truncated PHYLIP rows"))?;
            let (id, rest) = phylip_named(line, relaxed)?;
            if !ids.insert(id.clone()) {
                return Err(bad("duplicate PHYLIP identifier"));
            }
            let mut aligned = String::new();
            if !rest.trim().is_empty() {
                append(&mut aligned, &phylip_segment(rest, alphabet)?, width)?;
            }
            while aligned.len() < width {
                let line = lines
                    .next()
                    .ok_or_else(|| bad("truncated sequential PHYLIP row"))?;
                append(&mut aligned, &phylip_segment(line, alphabet)?, width)?;
            }
            rows.push((id, String::new(), aligned));
        }
    } else {
        let mut first_width = None;
        for _ in 0..n {
            let (id, rest) = phylip_named(
                lines
                    .next()
                    .ok_or_else(|| bad("truncated PHYLIP first block"))?,
                relaxed,
            )?;
            if !ids.insert(id.clone()) {
                return Err(bad("duplicate PHYLIP identifier"));
            }
            let part = phylip_segment(rest, alphabet)?;
            if first_width.is_some_and(|w| w != part.len()) {
                return Err(bad("unequal PHYLIP interleaved block widths"));
            }
            first_width = Some(part.len());
            if part.len() > width {
                return Err(bad("PHYLIP block exceeds NCHAR"));
            }
            rows.push((id, String::new(), part));
        }
        while rows[0].2.len() < width {
            let mut block_width = None;
            for row in &mut rows {
                let line = lines
                    .next()
                    .ok_or_else(|| bad("truncated PHYLIP continuation block"))?;
                // Official PHYLIP continuation rows contain sequence only, not
                // repeated names. This avoids ambiguous DNA-like name guessing.
                let part = phylip_segment(line, alphabet)?;
                if block_width.is_some_and(|w| w != part.len()) {
                    return Err(bad("unequal PHYLIP continuation widths"));
                }
                block_width = Some(part.len());
                append(&mut row.2, &part, width)?;
            }
        }
    }
    if lines.next().is_some() {
        return Err(bad("extra PHYLIP rows or concatenated alignments"));
    }
    let mut meta = AlignmentMetadata {
        source_header: Some(header.trim().into()),
        ..Default::default()
    };
    note(
        &mut meta,
        "PHYLIP profile: explicit alphabet, '-' gaps; numeric site counters/spacing ignored; no inferred tree or aligner confidence.",
    )?;
    Ok((rows, meta))
}

/// GCG's weighted, case-insensitive 57-position checksum, modulo 10000.
/// Gap glyphs are checked as actually written, before canonical '-' normalization.
pub fn msf_checksum(sequence: &str) -> u32 {
    sequence.bytes().enumerate().fold(0u32, |sum, (i, b)| {
        (sum + ((i % 57 + 1) as u32) * u32::from(b.to_ascii_uppercase())) % 10000
    })
}
fn field<'a>(fields: &[&'a str], key: &str) -> Result<&'a str, Error> {
    let mut hits = fields.windows(2).filter(|p| p[0] == key);
    let value = hits
        .next()
        .map(|p| p[1])
        .ok_or_else(|| bad(format!("MSF field {key} is missing")))?;
    if hits.next().is_some() {
        return Err(bad("duplicate MSF header field"));
    }
    Ok(value)
}
fn parse_msf(s: &str, alphabet: Alphabet) -> Result<(Rows, AlignmentMetadata), Error> {
    let mut lines = s.lines();
    let first = lines.next().ok_or_else(|| bad("empty MSF"))?.trim();
    let declared = match first {
        "!!AA_MULTIPLE_ALIGNMENT 1.0" => true,
        "!!NA_MULTIPLE_ALIGNMENT 1.0" => false,
        _ => return Err(bad("MSF requires its explicit AA/NA version-1 header")),
    };
    if declared != (alphabet == Alphabet::Protein) {
        return Err(bad("MSF AA/NA type conflicts with caller alphabet"));
    }
    let mut width = None;
    let mut total_check = None;
    let mut rows = Rows::new();
    let mut checks = Vec::new();
    let mut meta = AlignmentMetadata {
        source_header: Some(first.into()),
        ..Default::default()
    };
    let mut seen = BTreeSet::new();
    let mut separator = false;
    for line in lines.by_ref() {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        if line == "//" {
            separator = true;
            break;
        }
        if line.split_whitespace().any(|v| v == "MSF:") {
            if width.is_some() || !line.ends_with("..") {
                return Err(bad("duplicate or unterminated MSF dimensions header"));
            }
            let fields: Vec<_> = line.split_whitespace().collect();
            let w = number(field(&fields, "MSF:")?)?;
            let kind = field(&fields, "Type:")?;
            if kind != if declared { "P" } else { "N" } {
                return Err(bad("MSF Type conflicts with header/alphabet"));
            }
            let key = if fields.contains(&"Check:") {
                "Check:"
            } else {
                "CompCheck:"
            };
            let check = number(field(&fields, key)?)?;
            if check >= 10000 {
                return Err(bad("MSF checksum out of range"));
            }
            width = Some(w);
            total_check = Some(check as u32);
        } else if line.starts_with("Name:") {
            let fields: Vec<_> = line.split_whitespace().collect();
            let id = field(&fields, "Name:")?;
            label(id)?;
            if !seen.insert(id.to_owned()) || rows.len() >= MAX_ALIGNMENT_ROWS {
                return Err(bad("duplicate/excess MSF row"));
            }
            let len = number(field(&fields, "Len:")?)?;
            let check = number(field(&fields, "Check:")?)?;
            let weight = field(&fields, "Weight:")?;
            let parsed: f64 = weight.parse().map_err(|_| bad("invalid MSF weight"))?;
            if !parsed.is_finite() || parsed < 0.0 || check >= 10000 {
                return Err(bad("invalid MSF weight/checksum"));
            }
            checks.push((len, check as u32));
            rows.push((id.into(), String::new(), String::new()));
            meta.annotations.push(AlignmentAnnotation {
                scope: "GS".into(),
                row_identifier: Some(id.into()),
                tag: "WT".into(),
                value: weight.into(),
            });
        } else {
            note(
                &mut meta,
                format!("MSF preamble retained in original bytes: {line}"),
            )?;
        }
    }
    let width = width.ok_or_else(|| bad("missing MSF dimensions"))?;
    dimensions(rows.len(), width)?;
    if !separator || checks.iter().any(|(len, _)| *len != width) {
        return Err(bad(
            "MSF separator or declared per-row length disagrees with alignment width",
        ));
    }
    let mut raw = vec![String::new(); rows.len()];
    let mut at = 0;
    let mut block_width = None;
    for line in lines {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        if line
            .bytes()
            .all(|b| b.is_ascii_whitespace() || b.is_ascii_digit())
        {
            continue;
        }
        let mut fields = line.split_whitespace();
        let id = fields.next().ok_or_else(|| bad("missing MSF row"))?;
        if id != rows[at].0 {
            return Err(bad("MSF block row order/name disagrees with declared rows"));
        }
        let part: String = fields.collect();
        if part.is_empty() || part.bytes().any(|b| b.is_ascii_digit()) {
            return Err(bad("invalid MSF sequence block"));
        }
        if block_width.is_some_and(|w| w != part.len()) {
            return Err(bad("unequal MSF block widths"));
        }
        block_width = Some(part.len());
        append(&mut raw[at], &part.to_ascii_uppercase(), width)?;
        let normalized = part.replace(['.', '~'], "-");
        append(&mut rows[at].2, &molecular(&normalized, alphabet)?, width)?;
        at += 1;
        if at == rows.len() {
            at = 0;
            block_width = None;
        }
    }
    if at != 0 || rows.iter().any(|r| r.2.len() != width) {
        return Err(bad("incomplete MSF block/row"));
    }
    let mut sum = 0;
    let mut unspecified = false;
    for ((_, expected), raw) in checks.iter().zip(&raw) {
        let actual = msf_checksum(raw);
        sum = (sum + actual) % 10000;
        if *expected == 0 {
            unspecified = true
        } else if *expected != actual {
            return Err(bad("MSF sequence checksum mismatch"));
        }
    }
    if let Some(total) = total_check {
        if total != 0 && total != sum {
            return Err(bad("MSF aggregate checksum mismatch"));
        }
        if total == 0 {
            unspecified = true
        }
    }
    if unspecified {
        note(
            &mut meta,
            "MSF zero checksum fields were unspecified, not verified checksum claims.",
        )?;
    }
    note(
        &mut meta,
        "MSF '.' and '~' gap glyphs normalize to '-'; declared nonzero checksums verified against original gap glyphs.",
    )?;
    Ok((rows, meta))
}

#[derive(Debug)]
struct Token {
    value: String,
    line: usize,
    quoted: bool,
}
fn lex_nexus(s: &str) -> Result<Vec<Token>, Error> {
    let b = s.as_bytes();
    let mut at = 0;
    let mut line = 0;
    let mut out = Vec::new();
    while at < b.len() {
        if b[at].is_ascii_whitespace() {
            if b[at] == b'\n' {
                line += 1;
            }
            at += 1;
            continue;
        }
        if b[at] == b'[' {
            let mut depth = 1usize;
            at += 1;
            while at < b.len() && depth > 0 {
                match b[at] {
                    b'[' => {
                        depth += 1;
                        if depth > 128 {
                            return Err(bad("NEXUS comment nesting exceeds profile"));
                        }
                    }
                    b']' => depth -= 1,
                    b'\n' => line += 1,
                    _ => {}
                }
                at += 1;
            }
            if depth != 0 {
                return Err(bad("unterminated NEXUS comment"));
            }
            continue;
        }
        if out.len() >= MAX_TOKENS {
            return Err(Error::TooLarge);
        }
        let first_line = line;
        if b[at] == b'\'' {
            at += 1;
            let mut value = String::new();
            let mut closed = false;
            while at < b.len() {
                if b[at] == b'\'' {
                    if b.get(at + 1) == Some(&b'\'') {
                        value.push('\'');
                        at += 2;
                    } else {
                        at += 1;
                        closed = true;
                        break;
                    }
                } else {
                    if b[at].is_ascii_control() {
                        return Err(bad("control/newline inside a quoted NEXUS token"));
                    }
                    value.push(b[at] as char);
                    at += 1;
                }
            }
            if !closed {
                return Err(bad("unterminated quoted NEXUS token"));
            }
            out.push(Token {
                value,
                line: first_line,
                quoted: true,
            });
        } else if b";=(),{}:".contains(&b[at]) {
            out.push(Token {
                value: (b[at] as char).to_string(),
                line: first_line,
                quoted: false,
            });
            at += 1;
        } else {
            let start = at;
            while at < b.len() && !b[at].is_ascii_whitespace() && !b"[]';=(),{}:".contains(&b[at]) {
                at += 1;
            }
            if at == start {
                return Err(bad("unexpected NEXUS punctuation"));
            }
            out.push(Token {
                value: s[start..at].into(),
                line: first_line,
                quoted: false,
            });
        }
    }
    Ok(out)
}
fn is(t: &Token, s: &str) -> bool {
    !t.quoted && t.value.eq_ignore_ascii_case(s)
}
fn nexus_label(t: &Token) -> Result<String, Error> {
    let s = if t.quoted {
        t.value.clone()
    } else {
        t.value.replace('_', " ")
    };
    label(&s)?;
    Ok(s)
}
fn assignment<'a>(tokens: &'a [Token], key: &str) -> Result<Option<&'a Token>, Error> {
    let mut found = None;
    for (i, t) in tokens.iter().enumerate() {
        if is(t, key) {
            if found.is_some() || tokens.get(i + 1).is_none_or(|t| !is(t, "=")) {
                return Err(bad(format!("invalid or duplicate NEXUS {key} assignment")));
            }
            found = Some(
                tokens
                    .get(i + 2)
                    .ok_or_else(|| bad("missing NEXUS assignment value"))?,
            );
        }
    }
    Ok(found)
}
fn optional_dimension(
    tokens: &[Token],
    key: &str,
    current: &mut Option<usize>,
) -> Result<(), Error> {
    if let Some(value) = assignment(tokens, key)? {
        let n = number(&value.value)?;
        if current.is_some_and(|v| v != n) {
            return Err(bad("conflicting NEXUS dimensions"));
        }
        *current = Some(n);
    }
    Ok(())
}
fn dimensions_command(tokens: &[Token], taxa: bool) -> Result<(), Error> {
    let mut at = 0;
    let mut seen = BTreeSet::new();
    while at < tokens.len() {
        let key = tokens[at].value.to_ascii_uppercase();
        if tokens[at].quoted || !seen.insert(key.clone()) {
            return Err(bad("duplicate/quoted NEXUS dimension key"));
        }
        if key == "NEWTAXA" && !taxa {
            at += 1;
            continue;
        }
        if key != "NTAX" && (key != "NCHAR" || taxa) {
            return Err(bad("unsupported NEXUS DIMENSIONS option"));
        }
        if tokens.get(at + 1).is_none_or(|t| !is(t, "=")) {
            return Err(bad("NEXUS dimension requires an explicit value"));
        }
        number(
            &tokens
                .get(at + 2)
                .ok_or_else(|| bad("missing NEXUS dimension value"))?
                .value,
        )?;
        at += 3;
    }
    Ok(())
}
#[derive(Clone)]
struct NexusFormat {
    gap: u8,
    missing: u8,
    matchchar: Option<u8>,
    interleave: bool,
}
fn nexus_format(tokens: &[Token], alphabet: Alphabet) -> Result<NexusFormat, Error> {
    let mut f = NexusFormat {
        gap: b'-',
        missing: b'?',
        matchchar: None,
        interleave: false,
    };
    let mut datatype = false;
    let mut seen = BTreeSet::new();
    let mut at = 0;
    while at < tokens.len() {
        let key = tokens[at].value.to_ascii_uppercase();
        if tokens[at].quoted || !seen.insert(key.clone()) {
            return Err(bad("duplicate or quoted NEXUS FORMAT key"));
        }
        at += 1;
        let has_value = tokens.get(at).is_some_and(|t| is(t, "="));
        let value = if has_value {
            at += 1;
            let value = tokens
                .get(at)
                .ok_or_else(|| bad("missing NEXUS FORMAT value"))?;
            at += 1;
            Some(value)
        } else {
            None
        };
        match key.as_str() {
            "DATATYPE" => {
                let value = value.ok_or_else(|| bad("NEXUS DATATYPE needs a value"))?;
                let matches = match alphabet {
                    Alphabet::Dna => value.value.eq_ignore_ascii_case("DNA"),
                    Alphabet::Rna => value.value.eq_ignore_ascii_case("RNA"),
                    Alphabet::Protein => value.value.eq_ignore_ascii_case("PROTEIN"),
                };
                if !matches {
                    return Err(bad(
                        "NEXUS datatype conflicts with explicit DNA/RNA/protein alphabet",
                    ));
                }
                datatype = true;
            }
            "GAP" | "MISSING" | "MATCHCHAR" => {
                let value = value.ok_or_else(|| bad("NEXUS symbol option needs a value"))?;
                if value.value.len() != 1 {
                    return Err(bad(
                        "NEXUS gap/missing/match symbols must be one ASCII byte",
                    ));
                }
                let c = value.value.as_bytes()[0];
                if c.is_ascii_alphanumeric() || c.is_ascii_whitespace() || b"[]();=,'".contains(&c)
                {
                    return Err(bad("NEXUS symbol conflicts with molecular/token syntax"));
                }
                match key.as_str() {
                    "GAP" => f.gap = c,
                    "MISSING" => f.missing = c,
                    _ => f.matchchar = Some(c),
                }
            }
            "INTERLEAVE" => {
                f.interleave = match value {
                    None => true,
                    Some(t) if t.value.eq_ignore_ascii_case("YES") => true,
                    Some(t) if t.value.eq_ignore_ascii_case("NO") => false,
                    _ => return Err(bad("invalid NEXUS INTERLEAVE value")),
                }
            }
            "NOINTERLEAVE" if value.is_none() => f.interleave = false,
            "LABELS" if value.is_none() => {}
            _ => {
                return Err(bad(format!(
                    "unsupported NEXUS FORMAT option {key}; original semantics cannot be guessed"
                )));
            }
        }
    }
    if !datatype {
        return Err(bad("explicit NEXUS molecular DATATYPE is required"));
    }
    if f.gap == f.missing || f.matchchar.is_some_and(|c| c == f.gap || c == f.missing) {
        return Err(bad("overlapping NEXUS gap/missing/match symbols"));
    }
    Ok(f)
}
fn nexus_segment(
    tokens: &[Token],
    alphabet: Alphabet,
    f: &NexusFormat,
    first: Option<&str>,
    offset: usize,
) -> Result<String, Error> {
    let mut out = String::new();
    for t in tokens {
        if t.quoted {
            return Err(bad(
                "quoted NEXUS sequence symbols are outside molecular profile",
            ));
        }
        for c in t.value.bytes() {
            if c == f.missing {
                return Err(bad(
                    "NEXUS missing-data cells are not equivalent to bases or gaps; unsupported by this alignment profile",
                ));
            }
            let c = if c == f.gap {
                b'-'
            } else if Some(c) == f.matchchar {
                *first
                    .and_then(|s| s.as_bytes().get(offset + out.len()))
                    .ok_or_else(|| bad("NEXUS MATCHCHAR has no prior first-row cell"))?
            } else {
                c
            };
            out.push(c as char);
        }
    }
    if out.is_empty() {
        return Err(bad("empty NEXUS matrix row segment"));
    }
    molecular(&out, alphabet)
}
fn parse_nexus_matrix(
    tokens: &[Token],
    alphabet: Alphabet,
    n: usize,
    w: usize,
    f: &NexusFormat,
    taxa: &[String],
) -> Result<Rows, Error> {
    dimensions(n, w)?;
    let mut lines: Vec<&[Token]> = Vec::new();
    let mut start = 0;
    while start < tokens.len() {
        let mut end = start + 1;
        while end < tokens.len() && tokens[end].line == tokens[start].line {
            end += 1;
        }
        lines.push(&tokens[start..end]);
        start = end;
    }
    let mut at = 0;
    let mut rows = Rows::with_capacity(n);
    let mut names = BTreeSet::new();
    if !f.interleave {
        for row_index in 0..n {
            let line = lines
                .get(at)
                .ok_or_else(|| bad("truncated sequential NEXUS matrix"))?;
            at += 1;
            let id = nexus_label(&line[0])?;
            if !names.insert(id.clone()) || (!taxa.is_empty() && taxa.get(row_index) != Some(&id)) {
                return Err(bad("duplicate NEXUS row or TAXLABELS/order mismatch"));
            }
            let mut aligned = String::new();
            if line.len() > 1 {
                let p = nexus_segment(
                    &line[1..],
                    alphabet,
                    f,
                    rows.first().map(|r| r.2.as_str()),
                    0,
                )?;
                append(&mut aligned, &p, w)?;
            }
            while aligned.len() < w {
                let line = lines
                    .get(at)
                    .ok_or_else(|| bad("truncated NEXUS sequential row"))?;
                at += 1;
                let p = nexus_segment(
                    line,
                    alphabet,
                    f,
                    rows.first().map(|r| r.2.as_str()),
                    aligned.len(),
                )?;
                append(&mut aligned, &p, w)?;
            }
            rows.push((id, String::new(), aligned));
        }
    } else {
        let mut width = None;
        for row_index in 0..n {
            let line = lines
                .get(at)
                .ok_or_else(|| bad("truncated first NEXUS interleaved block"))?;
            at += 1;
            let id = nexus_label(&line[0])?;
            if !names.insert(id.clone()) || (!taxa.is_empty() && taxa.get(row_index) != Some(&id)) {
                return Err(bad("duplicate NEXUS row or TAXLABELS/order mismatch"));
            }
            let p = nexus_segment(
                &line[1..],
                alphabet,
                f,
                rows.first().map(|r| r.2.as_str()),
                0,
            )?;
            if width.is_some_and(|v| v != p.len()) {
                return Err(bad("unequal NEXUS block widths"));
            }
            width = Some(p.len());
            if p.len() > w {
                return Err(bad("NEXUS block exceeds NCHAR"));
            }
            rows.push((id, String::new(), p));
        }
        while rows[0].2.len() < w {
            let mut width = None;
            let mut named = None;
            for i in 0..n {
                let line = lines
                    .get(at)
                    .ok_or_else(|| bad("truncated NEXUS continuation block"))?;
                at += 1;
                let name = if line[0].quoted {
                    line[0].value.clone()
                } else {
                    line[0].value.replace('_', " ")
                };
                let has_name = name == rows[i].0;
                if named.is_some_and(|v| v != has_name) {
                    return Err(bad("mixed labelled/unlabelled NEXUS continuation block"));
                }
                named = Some(has_name);
                if !has_name && names.contains(&name) {
                    return Err(bad("NEXUS continuation row order mismatch"));
                }
                let body = if has_name { &line[1..] } else { line };
                let p = nexus_segment(
                    body,
                    alphabet,
                    f,
                    if i == 0 {
                        None
                    } else {
                        Some(rows[0].2.as_str())
                    },
                    rows[i].2.len(),
                )?;
                if width.is_some_and(|v| v != p.len()) {
                    return Err(bad("unequal NEXUS continuation widths"));
                }
                width = Some(p.len());
                append(&mut rows[i].2, &p, w)?;
            }
        }
    }
    if at != lines.len() || rows.iter().any(|r| r.2.len() != w) {
        return Err(bad("extra/incomplete NEXUS matrix rows"));
    }
    Ok(rows)
}
fn parse_nexus(s: &str, alphabet: Alphabet) -> Result<(Rows, AlignmentMetadata), Error> {
    let tokens = lex_nexus(s)?;
    if tokens.first().is_none_or(|t| !is(t, "#NEXUS")) {
        return Err(bad("missing #NEXUS header"));
    }
    let mut at = 1;
    let mut meta = AlignmentMetadata {
        source_header: Some("#NEXUS".into()),
        ..Default::default()
    };
    let mut ntax = None;
    let mut nchar = None;
    let mut taxa = Vec::new();
    let mut matrix = None;
    let mut format = None;
    let mut saw_characters = false;
    let mut saw_taxa = false;
    while at < tokens.len() {
        if !is(&tokens[at], "BEGIN") || tokens.get(at + 2).is_none_or(|t| !is(t, ";")) {
            return Err(bad("NEXUS expects BEGIN block;"));
        }
        let block = tokens
            .get(at + 1)
            .ok_or_else(|| bad("missing NEXUS block name"))?
            .value
            .to_ascii_uppercase();
        at += 3;
        let supported = matches!(block.as_str(), "DATA" | "CHARACTERS" | "TAXA");
        if block == "TAXA" {
            if saw_taxa {
                return Err(bad("duplicate NEXUS TAXA block"));
            }
            saw_taxa = true;
        }
        if matches!(block.as_str(), "DATA" | "CHARACTERS") {
            if saw_characters {
                return Err(bad("multiple molecular matrices are unsupported"));
            }
            saw_characters = true;
        }
        if !supported {
            note(
                &mut meta,
                format!(
                    "NEXUS {block} block retained in original bytes, not interpreted or exported canonically."
                ),
            )?;
        }
        let mut ended = false;
        while at < tokens.len() {
            let start = at;
            while at < tokens.len() && !is(&tokens[at], ";") {
                at += 1;
            }
            if at == tokens.len() {
                return Err(bad("unterminated NEXUS command/block"));
            }
            let command = &tokens[start..at];
            at += 1;
            if command.is_empty() {
                continue;
            }
            if is(&command[0], "END") || is(&command[0], "ENDBLOCK") {
                if command.len() != 1 {
                    return Err(bad("invalid NEXUS block terminator"));
                }
                ended = true;
                break;
            }
            if !supported {
                continue;
            }
            let head = command[0].value.to_ascii_uppercase();
            match (block.as_str(), head.as_str()) {
                (_, "DIMENSIONS") => {
                    dimensions_command(&command[1..], block == "TAXA")?;
                    optional_dimension(&command[1..], "NTAX", &mut ntax)?;
                    if block != "TAXA" {
                        optional_dimension(&command[1..], "NCHAR", &mut nchar)?;
                    }
                }
                ("TAXA", "TAXLABELS") => {
                    if !taxa.is_empty() {
                        return Err(bad("duplicate NEXUS TAXLABELS"));
                    }
                    if command.len() - 1 > MAX_ALIGNMENT_ROWS {
                        return Err(Error::TooLarge);
                    }
                    for t in &command[1..] {
                        taxa.push(nexus_label(t)?);
                    }
                }
                ("DATA" | "CHARACTERS", "FORMAT") => {
                    if format.is_some() {
                        return Err(bad("duplicate NEXUS FORMAT"));
                    }
                    format = Some(nexus_format(&command[1..], alphabet)?);
                }
                ("DATA" | "CHARACTERS", "MATRIX") => {
                    if matrix.is_some() {
                        return Err(bad("duplicate NEXUS MATRIX"));
                    }
                    matrix = Some(&command[1..]);
                }
                (_, "TITLE") => note(
                    &mut meta,
                    format!(
                        "NEXUS {block} {head} command retained in original bytes, not applied."
                    ),
                )?,
                _ => {
                    return Err(bad(format!(
                        "unsupported {head} command in NEXUS {block}; cannot assume it has no matrix semantics"
                    )));
                }
            }
        }
        if !ended {
            return Err(bad("NEXUS block has no END"));
        }
    }
    let n = ntax.ok_or_else(|| bad("NEXUS NTAX required"))?;
    let w = nchar.ok_or_else(|| bad("NEXUS NCHAR required"))?;
    dimensions(n, w)?;
    if !taxa.is_empty()
        && (taxa.len() != n || taxa.iter().collect::<BTreeSet<_>>().len() != taxa.len())
    {
        return Err(bad("NEXUS TAXLABELS count/identity mismatch"));
    }
    let format = format.ok_or_else(|| bad("NEXUS FORMAT required"))?;
    let rows = parse_nexus_matrix(
        matrix.ok_or_else(|| bad("NEXUS molecular MATRIX required"))?,
        alphabet,
        n,
        w,
        &format,
        &taxa,
    )?;
    note(
        &mut meta,
        "NEXUS molecular matrix only; trees, sets and private block semantics are not inferred. Exact original bytes retain all blocks/comments.",
    )?;
    Ok((rows, meta))
}

struct Output {
    text: String,
}
impl Output {
    fn new() -> Self {
        Self {
            text: String::new(),
        }
    }
    fn line(&mut self, line: &str) -> Result<(), Error> {
        if line.len() > MAX_ALIGNMENT_LINE
            || line.len() + 1 > MAX_ALIGNMENT_BYTES.saturating_sub(self.text.len())
        {
            return Err(Error::TooLarge);
        }
        self.text.push_str(line);
        self.text.push('\n');
        Ok(())
    }
}
pub fn export_alignment(
    doc: &AlignmentDocument,
    format: AlignmentFormat,
) -> Result<AlignmentExport, Error> {
    export_alignment_with_options(doc, format, &ExportOptions::default())
}
pub fn export_alignment_with_options(
    doc: &AlignmentDocument,
    format: AlignmentFormat,
    options: &ExportOptions,
) -> Result<AlignmentExport, Error> {
    doc.validate()?;
    export_rows(&doc.rows, &doc.metadata, doc.alphabet, format, options)
}
/// Shared dispatcher hook also supports scientifically computed rows without
/// fabricating an imported source record. Structural validation remains shared.
pub(crate) fn export_rows(
    rows: &[AlignmentRow],
    metadata: &AlignmentMetadata,
    alphabet: Alphabet,
    format: AlignmentFormat,
    options: &ExportOptions,
) -> Result<AlignmentExport, Error> {
    if matches!(
        format,
        AlignmentFormat::Pir | AlignmentFormat::Gde | AlignmentFormat::Selex
    ) {
        return crate::formats::alignment::pir_gde_selex::export_rows(
            rows, metadata, alphabet, format, options,
        );
    }

    base::validate_export_content(alphabet, rows, metadata)?;
    if options.block_width == 0 || options.block_width > 1000 {
        return Err(bad("exchange block width must be 1..1000"));
    }
    let mut out = Output::new();
    let mut renamings = BTreeMap::new();
    let mut losses=vec!["Canonical export changes original byte layout, case, comments and header text; use original_bytes for exact evidence.".into()];
    if rows.iter().any(|r| !r.description.is_empty()) {
        losses.push("Row descriptions are not represented by this exchange profile.".into());
    }
    if metadata.clustal_consensus.is_some() {
        losses.push("Imported CLUSTAL consensus symbols are not represented.".into());
    }
    for a in &metadata.annotations {
        if a.tag == "CG2_EXCHANGE_NOTE" && a.value.contains("not interpreted") {
            losses.push(a.value.clone());
        }
    }
    let only_weights = format == AlignmentFormat::Msf
        && metadata
            .annotations
            .iter()
            .all(|a| a.scope == "GS" && a.tag == "WT" || a.tag == "CG2_EXCHANGE_NOTE");
    if !metadata.annotations.is_empty() && !only_weights {
        losses.push("General alignment annotations and unsupported original block semantics are not emitted; they remain in the native/original document.".into());
    }
    if let Some((relaxed, interleaved)) = phylip(format) {
        if alphabet == Alphabet::Protein
            && rows
                .iter()
                .any(|r| r.aligned.bytes().any(|b| b"JOU".contains(&b)))
        {
            return Err(bad(
                "PHYLIP protein J/O/U extensions are outside this profile",
            ));
        }
        let mut ids = Vec::with_capacity(rows.len());
        let mut seen = BTreeSet::new();
        for row in rows {
            let mut id = row.identifier.clone();
            if relaxed {
                if id.bytes().any(|b| b.is_ascii_whitespace()) {
                    return Err(bad("relaxed PHYLIP identifiers cannot contain whitespace"));
                }
            } else {
                if id.trim() != id {
                    return Err(bad(
                        "strict PHYLIP cannot preserve leading/trailing name whitespace",
                    ));
                }
                if id.len() > 10 {
                    if !options.truncate_phylip_identifiers {
                        return Err(bad(
                            "strict PHYLIP names exceed ten bytes; explicitly opt in to reported truncation or choose relaxed PHYLIP",
                        ));
                    }
                    id.truncate(10);
                    if id.trim() != id || id.is_empty() {
                        return Err(bad(
                            "truncated PHYLIP name ends in padding whitespace; rename explicitly",
                        ));
                    }
                    renamings.insert(row.identifier.clone(), id.clone());
                }
            }
            if !seen.insert(id.clone()) {
                return Err(bad(
                    "PHYLIP output identifiers collide, including after ten-byte truncation",
                ));
            }
            ids.push(id);
        }
        if !renamings.is_empty() {
            losses.push("Explicitly requested strict-PHYLIP ten-byte identifier truncation; consult renamings.".into());
        }
        losses.push("PHYLIP does not declare the molecular alphabet; preserve the explicit alphabet separately.".into());
        out.line(&format!("{} {}", rows.len(), rows[0].aligned.len()))?;
        if interleaved {
            for start in (0..rows[0].aligned.len()).step_by(options.block_width) {
                if start > 0 {
                    out.line("")?;
                }
                let end = (start + options.block_width).min(rows[0].aligned.len());
                for (row, id) in rows.iter().zip(&ids) {
                    let segment = &row.aligned[start..end];
                    if start == 0 {
                        out.line(&if relaxed {
                            format!("{id} {segment}")
                        } else {
                            format!("{id:<10}{segment}")
                        })?;
                    } else {
                        out.line(segment)?;
                    }
                }
            }
        } else {
            for (row, id) in rows.iter().zip(&ids) {
                for start in (0..row.aligned.len()).step_by(options.block_width) {
                    let end = (start + options.block_width).min(row.aligned.len());
                    let segment = &row.aligned[start..end];
                    if start == 0 {
                        out.line(&if relaxed {
                            format!("{id} {segment}")
                        } else {
                            format!("{id:<10}{segment}")
                        })?;
                    } else {
                        out.line(segment)?;
                    }
                }
            }
        }
    } else {
        match format {
            AlignmentFormat::Msf => {
                if rows
                    .iter()
                    .any(|r| r.identifier.bytes().any(|b| b.is_ascii_whitespace()))
                {
                    return Err(bad("MSF identifiers cannot contain whitespace"));
                }
                let sequences: Vec<_> = rows.iter().map(|r| r.aligned.replace('-', ".")).collect();
                let checks: Vec<_> = sequences.iter().map(|s| msf_checksum(s)).collect();
                let sum = checks.iter().fold(0u32, |n, v| (n + v) % 10000);
                out.line(if alphabet == Alphabet::Protein {
                    "!!AA_MULTIPLE_ALIGNMENT 1.0"
                } else {
                    "!!NA_MULTIPLE_ALIGNMENT 1.0"
                })?;
                out.line("")?;
                out.line(&format!(
                    " collagene MSF: {} Type: {} Check: {sum} ..",
                    rows[0].aligned.len(),
                    if alphabet == Alphabet::Protein {
                        'P'
                    } else {
                        'N'
                    }
                ))?;
                let row_indices: std::collections::HashMap<&str, usize> = rows
                    .iter()
                    .enumerate()
                    .map(|(i, row)| (row.identifier.as_str(), i))
                    .collect();
                let mut weights = vec![None; rows.len()];
                for annotation in &metadata.annotations {
                    if annotation.scope == "GS" && annotation.tag == "WT" {
                        let id = annotation
                            .row_identifier
                            .as_deref()
                            .ok_or_else(|| bad("MSF weight has no row identifier"))?;
                        let i = *row_indices
                            .get(id)
                            .ok_or_else(|| bad("MSF weight refers to absent row"))?;
                        if weights[i].replace(annotation.value.as_str()).is_some() {
                            return Err(bad("ambiguous repeated MSF weight metadata"));
                        }
                    }
                }
                for (row_index, (row, check)) in rows.iter().zip(&checks).enumerate() {
                    let weight = weights[row_index].unwrap_or("1.00");
                    let parsed: f64 = weight
                        .parse()
                        .map_err(|_| bad("MSF weight metadata is not numeric"))?;
                    if !parsed.is_finite() || parsed < 0.0 {
                        return Err(bad("MSF weight metadata is not finite/nonnegative"));
                    }
                    out.line(&format!(
                        " Name: {} Len: {} Check: {check} Weight: {weight}",
                        row.identifier,
                        row.aligned.len()
                    ))?;
                }
                out.line("")?;
                out.line("//")?;
                for start in (0..rows[0].aligned.len()).step_by(options.block_width) {
                    out.line("")?;
                    let end = (start + options.block_width).min(rows[0].aligned.len());
                    for (row, sequence) in rows.iter().zip(&sequences) {
                        out.line(&format!("{} {}", row.identifier, &sequence[start..end]))?;
                    }
                }
                losses.push("MSF checksums are regenerated over '.' gaps; N type does not distinguish explicit DNA from RNA.".into());
            }
            AlignmentFormat::Nexus => {
                out.line("#NEXUS")?;
                out.line("BEGIN DATA;")?;
                out.line(&format!(
                    "  DIMENSIONS NTAX={} NCHAR={};",
                    rows.len(),
                    rows[0].aligned.len()
                ))?;
                out.line(&format!(
                    "  FORMAT DATATYPE={} GAP=- MISSING=? INTERLEAVE=YES;",
                    match alphabet {
                        Alphabet::Dna => "DNA",
                        Alphabet::Rna => "RNA",
                        Alphabet::Protein => "PROTEIN",
                    }
                ))?;
                out.line("  MATRIX")?;
                for start in (0..rows[0].aligned.len()).step_by(options.block_width) {
                    if start > 0 {
                        out.line("")?;
                    }
                    let end = (start + options.block_width).min(rows[0].aligned.len());
                    for row in rows {
                        let quoted = row.identifier.replace('\'', "''");
                        out.line(&format!("    '{quoted}' {}", &row.aligned[start..end]))?;
                    }
                }
                out.line("  ;")?;
                out.line("END;")?;
                losses.push("NEXUS canonical output contains one molecular DATA matrix only; original trees, sets, commands, comments and private blocks are not reconstructed.".into());
            }
            _ => return Err(bad("not an alignment-exchange export profile")),
        }
    }
    Ok(AlignmentExport {
        text: out.text,
        losses,
        renamings,
    })
}
