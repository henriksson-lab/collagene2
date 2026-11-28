//! Independent bounded PIR/NBRF, GDE flat, and SELEX alignment exchange.
//! This module parses stored alignments; it never computes an alignment.
use crate::formats::alignment::profiles::ExportOptions;
use crate::formats::alignment::{
    self as base, AlignmentAnnotation, AlignmentDocument, AlignmentExport, AlignmentFormat,
    AlignmentMetadata, AlignmentRow, MAX_ALIGNMENT_ANNOTATIONS, MAX_ALIGNMENT_BYTES,
    MAX_ALIGNMENT_CELLS, MAX_ALIGNMENT_LINE, MAX_ALIGNMENT_ROWS,
};
use collagene_model::{Alphabet, Error, Sequence};
use std::collections::{BTreeMap, BTreeSet};
type Rows = Vec<(String, String, String)>;
fn bad(message: &str) -> Error {
    Error::Format(format!("alignment exchange: {message}"))
}
fn unsupported(message: &str) -> Error {
    Error::Unsupported(format!(
        "{message}; retain the native alignment and exact original source"
    ))
}
fn token(value: &str) -> bool {
    !value.is_empty()
        && value.is_ascii()
        && !value
            .bytes()
            .any(|b| b.is_ascii_whitespace() || b.is_ascii_control())
}
fn text(bytes: &[u8]) -> Result<&str, Error> {
    if bytes.len() > MAX_ALIGNMENT_BYTES {
        return Err(Error::TooLarge);
    }
    let text = std::str::from_utf8(bytes).map_err(|_| bad("ASCII text required"))?;
    if !text.is_ascii()
        || text.lines().any(|l| l.len() > MAX_ALIGNMENT_LINE)
        || text
            .bytes()
            .any(|b| b.is_ascii_control() && !b"\r\n\t".contains(&b))
    {
        return Err(bad("unsupported control, encoding, or line bound"));
    }
    Ok(text)
}
fn dimensions(rows: usize, width: usize) -> Result<(), Error> {
    if rows == 0
        || rows > MAX_ALIGNMENT_ROWS
        || width == 0
        || rows
            .checked_mul(width)
            .is_none_or(|n| n > MAX_ALIGNMENT_CELLS)
    {
        Err(Error::TooLarge)
    } else {
        Ok(())
    }
}
fn annotation(
    meta: &mut AlignmentMetadata,
    scope: &str,
    row: Option<&str>,
    tag: &str,
    value: impl Into<String>,
) -> Result<(), Error> {
    if meta.annotations.len() >= MAX_ALIGNMENT_ANNOTATIONS {
        return Err(Error::TooLarge);
    }
    meta.annotations.push(AlignmentAnnotation {
        scope: scope.into(),
        row_identifier: row.map(str::to_string),
        tag: tag.into(),
        value: value.into(),
    });
    Ok(())
}
fn segment(
    value: &str,
    alphabet: Alphabet,
    gaps: &[u8],
    strip_space: bool,
) -> Result<String, Error> {
    let mut sequence = String::new();
    for b in value.bytes() {
        if strip_space && b.is_ascii_whitespace() {
            continue;
        }
        sequence.push(if gaps.contains(&b) {
            '-'
        } else {
            b.to_ascii_uppercase() as char
        });
    }
    let residues: String = sequence.chars().filter(|c| *c != '-').collect();
    if !residues.is_empty() {
        Sequence::new(alphabet, &residues)?;
    }
    Ok(sequence)
}
fn type_code(code: &str, alphabet: Alphabet) -> Result<(), Error> {
    let expected = match code {
        "P1" | "F1" => Alphabet::Protein,
        "DL" => Alphabet::Dna,
        "RL" => Alphabet::Rna,
        "DC" | "RC" => {
            return Err(unsupported(
                "circular PIR rows cannot be represented by this linear extraction profile",
            ));
        }
        _ => return Err(unsupported("unsupported PIR/NBRF sequence type")),
    };
    if expected != alphabet {
        return Err(bad(
            "PIR type code disagrees with explicitly selected alphabet",
        ));
    }
    Ok(())
}
fn parse_pir(text: &str, alphabet: Alphabet) -> Result<(Rows, AlignmentMetadata), Error> {
    let lines: Vec<_> = text.lines().collect();
    let mut i = 0;
    let mut rows = Vec::new();
    let mut meta = AlignmentMetadata::default();
    let mut cells = 0usize;
    while i < lines.len() {
        if lines[i].trim().is_empty() {
            i += 1;
            continue;
        }
        let header = lines[i]
            .strip_prefix('>')
            .ok_or_else(|| bad("PIR record must start >TYPE;identifier"))?;
        let (code, id) = header
            .split_once(';')
            .ok_or_else(|| bad("PIR type separator missing"))?;
        type_code(code, alphabet)?;
        if !token(id) {
            return Err(bad("PIR identifier must be one nonempty token"));
        }
        if rows.len() >= MAX_ALIGNMENT_ROWS {
            return Err(Error::TooLarge);
        }
        i += 1;
        let description = lines
            .get(i)
            .ok_or_else(|| bad("PIR description line missing"))?
            .to_string();
        i += 1;
        let mut raw = String::new();
        let mut terminated = false;
        while i < lines.len() {
            let line = lines[i];
            i += 1;
            if line.starts_with('>') {
                return Err(bad("PIR record missing terminal asterisk"));
            }
            if let Some((before, after)) = line.split_once('*') {
                if !after.trim().is_empty() {
                    return Err(bad("PIR data after terminal asterisk"));
                }
                raw.push_str(before);
                terminated = true;
                break;
            }
            raw.push_str(line);
            if raw.len() > MAX_ALIGNMENT_CELLS {
                return Err(Error::TooLarge);
            }
        }
        if !terminated {
            return Err(bad("PIR terminal asterisk missing"));
        }
        let aligned = segment(&raw, alphabet, b"-", true)?;
        if aligned.is_empty() {
            return Err(bad("empty PIR row"));
        }
        cells = cells
            .checked_add(aligned.len())
            .filter(|n| *n <= MAX_ALIGNMENT_CELLS)
            .ok_or(Error::TooLarge)?;
        annotation(&mut meta, "GS", Some(id), "PIR_TYPE", code)?;
        while i < lines.len() && !lines[i].starts_with('>') {
            if !lines[i].trim().is_empty() {
                annotation(
                    &mut meta,
                    "GS",
                    Some(id),
                    "PIR_TEXT",
                    lines[i].replace('\t', " "),
                )?;
            }
            i += 1;
        }
        rows.push((id.into(), description, aligned));
    }
    if let Some(row) = rows.first() {
        dimensions(rows.len(), row.2.len())?;
    } else {
        return Err(bad("no PIR rows"));
    }
    Ok((rows, meta))
}
fn parse_gde(text: &str, alphabet: Alphabet) -> Result<(Rows, AlignmentMetadata), Error> {
    if text.trim_start().starts_with('{') {
        return Err(unsupported(
            "tagged GDE records require a different explicit profile",
        ));
    }
    let mut rows: Rows = Vec::new();
    let mut offsets = Vec::new();
    let mut current: Option<(String, usize, String)> = None;
    let mut meta = AlignmentMetadata::default();
    let mut charged = 0usize;
    fn finish(
        current: &mut Option<(String, usize, String)>,
        rows: &mut Rows,
        offsets: &mut Vec<usize>,
        meta: &mut AlignmentMetadata,
    ) -> Result<(), Error> {
        if let Some((id, offset, aligned)) = current.take() {
            if aligned.is_empty() {
                return Err(bad("empty GDE sequence"));
            }
            annotation(meta, "GS", Some(&id), "GDE_OFFSET", offset.to_string())?;
            annotation(
                meta,
                "GS",
                Some(&id),
                "GDE_SOURCE_WIDTH",
                aligned.len().to_string(),
            )?;
            rows.push((id, String::new(), aligned));
            offsets.push(offset);
        }
        Ok(())
    }
    for line in text.lines() {
        if line.trim().is_empty() {
            continue;
        }
        if matches!(line.as_bytes()[0], b'#' | b'%' | b'@' | b'"') {
            finish(&mut current, &mut rows, &mut offsets, &mut meta)?;
            if rows.len() >= MAX_ALIGNMENT_ROWS {
                return Err(Error::TooLarge);
            }
            let code = line.as_bytes()[0];
            if matches!(code, b'@' | b'"') {
                return Err(unsupported(
                    "GDE text/mask rows are not biological sequence rows",
                ));
            }
            if (code == b'%') != (alphabet == Alphabet::Protein) {
                return Err(bad("GDE type marker disagrees with explicit alphabet"));
            }
            let header = &line[1..];
            let (id, offset) = if let Some((id, tail)) = header.rsplit_once('(') {
                let value = tail
                    .strip_suffix(')')
                    .ok_or_else(|| bad("invalid GDE offset suffix"))?;
                if value.starts_with('-') {
                    return Err(unsupported(
                        "negative GDE offsets are outside the nonnegative alignment-column profile",
                    ));
                }
                (
                    id.trim_end(),
                    value
                        .parse::<usize>()
                        .map_err(|_| bad("invalid GDE offset"))?,
                )
            } else {
                (header, 0)
            };
            if !token(id) || id.contains(['(', ')', '#', '%', '@', '"']) {
                return Err(bad("GDE short name is ambiguous"));
            }
            if offset > MAX_ALIGNMENT_CELLS {
                return Err(Error::TooLarge);
            }
            current = Some((id.into(), offset, String::new()));
        } else {
            let value = segment(line, alphabet, b"-~", true)?;
            let row = current
                .as_mut()
                .ok_or_else(|| bad("GDE sequence before type/name"))?;
            charged = charged
                .checked_add(value.len())
                .filter(|n| *n <= MAX_ALIGNMENT_CELLS)
                .ok_or(Error::TooLarge)?;
            row.2.push_str(&value);
        }
    }
    finish(&mut current, &mut rows, &mut offsets, &mut meta)?;
    let width = rows
        .iter()
        .zip(&offsets)
        .map(|(r, o)| r.2.len().checked_add(*o).ok_or(Error::TooLarge))
        .collect::<Result<Vec<_>, _>>()?
        .into_iter()
        .max()
        .ok_or_else(|| bad("no GDE rows"))?;
    dimensions(rows.len(), width)?;
    for (row, offset) in rows.iter_mut().zip(offsets) {
        let mut aligned = "-".repeat(offset);
        aligned.push_str(&row.2);
        aligned.extend(std::iter::repeat_n('-', width - aligned.len()));
        row.2 = aligned;
    }
    annotation(
        &mut meta,
        "GF",
        None,
        "GDE_LAYOUT",
        "Nonnegative offsets are leading alignment gaps; shorter source extents are padded with trailing gaps. Offsets are not external biological coordinates.",
    )?;
    Ok((rows, meta))
}
// A SELEX block is interpreted at common text columns. Exterior whitespace is
// layout; internal spaces are gaps. Explicit glyphs anchor outer all-gap columns.
struct SelexLine {
    scope: &'static str,
    row: Option<String>,
    tag: String,
    raw: String,
    prefix_end: usize,
    data_start: Option<usize>,
}
fn datum(
    line: &str,
    scope: &'static str,
    row: Option<String>,
    tag: &str,
    prefix_end: usize,
) -> Result<SelexLine, Error> {
    if line.len() >= 4096 || line.contains('\t') {
        return Err(unsupported(
            "SELEX data lines must be <4096 columns without ambiguous tab stops",
        ));
    }
    let start = line[prefix_end..]
        .bytes()
        .position(|b| !b.is_ascii_whitespace())
        .map(|i| i + prefix_end);
    Ok(SelexLine {
        scope,
        row,
        tag: tag.into(),
        raw: line.trim_end().into(),
        prefix_end,
        data_start: start,
    })
}
fn word<'a>(value: &mut &'a str) -> Option<&'a str> {
    let s = value.trim_start();
    if s.is_empty() {
        return None;
    }
    let end = s.find(char::is_whitespace).unwrap_or(s.len());
    let result = &s[..end];
    *value = &s[end..];
    Some(result)
}
#[derive(Default)]
struct SelexState {
    rows: Rows,
    columns: BTreeMap<(String, Option<String>, String), String>,
    keys: Option<BTreeSet<(String, Option<String>, String)>>,
    case: Vec<String>,
    blocks: usize,
    cells: usize,
    column_bytes: usize,
}
impl SelexState {
    fn flush(&mut self, block: &mut Vec<SelexLine>, alphabet: Alphabet) -> Result<(), Error> {
        if block.is_empty() {
            return Ok(());
        }
        let left = block
            .iter()
            .filter_map(|l| l.data_start)
            .min()
            .ok_or_else(|| bad("SELEX block contains no explicit alignment symbols"))?;
        let right = block.iter().map(|l| l.raw.len()).max().unwrap_or(left);
        if right <= left || block.iter().any(|l| l.prefix_end >= left) {
            return Err(unsupported(
                "SELEX labels overlap alignment columns; align names or use explicit gap glyphs",
            ));
        }
        let width = right - left;
        let row_count = block.iter().filter(|l| l.scope == "ROW").count();
        dimensions(row_count, width)?;
        let prospective = self
            .cells
            .checked_add(row_count.checked_mul(width).ok_or(Error::TooLarge)?)
            .filter(|n| *n <= MAX_ALIGNMENT_CELLS)
            .ok_or(Error::TooLarge)?;
        let mut names = Vec::new();
        let mut segments = Vec::new();
        let mut keys = BTreeSet::new();
        let mut masks = Vec::new();
        for line in block.iter() {
            let slice = if line.raw.len() > left {
                &line.raw[left..]
            } else {
                ""
            };
            let mut value = slice.to_string();
            value.extend(std::iter::repeat_n(' ', width - value.len()));
            if line.scope == "ROW" {
                let id = line.row.as_ref().expect("row");
                if names.contains(id) {
                    return Err(bad("duplicate SELEX row in block"));
                }
                names.push(id.clone());
                let aligned = segment(&value, alphabet, b"-_. ", false)?;
                let mask: String = value
                    .bytes()
                    .map(|b| {
                        if b"-_. ".contains(&b) {
                            '.'
                        } else if b.is_ascii_lowercase() {
                            'l'
                        } else {
                            'u'
                        }
                    })
                    .collect();
                segments.push(aligned);
                masks.push(mask);
            } else {
                let key = (line.scope.into(), line.row.clone(), line.tag.clone());
                if !keys.insert(key.clone()) {
                    return Err(bad("duplicate SELEX column annotation"));
                }
                self.column_bytes = self
                    .column_bytes
                    .checked_add(width)
                    .filter(|n| *n <= MAX_ALIGNMENT_BYTES)
                    .ok_or(Error::TooLarge)?;
                self.columns
                    .entry(key)
                    .or_default()
                    .extend(value.chars().map(|c| if c == ' ' { '.' } else { c }));
            }
        }
        dimensions(names.len(), width)?;
        if self.blocks == 0 {
            self.rows = names
                .iter()
                .cloned()
                .zip(segments)
                .map(|(id, seq)| (id, String::new(), seq))
                .collect();
            self.case = masks;
            self.keys = Some(keys);
        } else {
            if self
                .rows
                .iter()
                .map(|r| r.0.as_str())
                .ne(names.iter().map(String::as_str))
                || self.keys.as_ref() != Some(&keys)
            {
                return Err(bad("SELEX block row order or annotation coverage changed"));
            }
            for ((row, seq), (mask, part)) in self
                .rows
                .iter_mut()
                .zip(segments)
                .zip(self.case.iter_mut().zip(masks))
            {
                row.2.push_str(&seq);
                mask.push_str(&part);
            }
        }
        self.cells = prospective;
        self.blocks += 1;
        block.clear();
        Ok(())
    }
}
fn selex_header_value(tag: &str, value: &str) -> Result<String, Error> {
    if matches!(tag, "ID" | "AC") && !token(value) {
        return Err(bad("SELEX ID/AC must contain one token"));
    }
    if matches!(tag, "GA" | "TC" | "NC") {
        let words: Vec<_> = value.trim_end_matches(';').split_whitespace().collect();
        if words.len() != 2
            || !words
                .iter()
                .all(|v| v.parse::<f64>().is_ok_and(f64::is_finite))
        {
            return Err(bad("SELEX cutoffs need two finite source numbers"));
        }
        return Ok(words.join(" "));
    }
    Ok(value.into())
}
fn parse_selex(text: &str, alphabet: Alphabet) -> Result<(Rows, AlignmentMetadata), Error> {
    let mut state = SelexState::default();
    let mut block = Vec::new();
    let mut meta = AlignmentMetadata::default();
    let mut headers = Vec::new();
    let mut last_row = None;
    let mut seen_data = false;
    let mut global_headers = BTreeSet::new();
    for line in text.lines() {
        if block.len() > MAX_ALIGNMENT_ROWS * 3 + 2 {
            return Err(Error::TooLarge);
        }
        if line.trim().is_empty() {
            state.flush(&mut block, alphabet)?;
            last_row = None;
            continue;
        }
        if let Some(machine) = line.strip_prefix("#=") {
            let end = machine
                .find(char::is_whitespace)
                .ok_or_else(|| bad("empty SELEX machine comment"))?;
            let tag = &machine[..end];
            let prefix_end = 2 + end;
            match tag {
                "RF" | "CS" => {
                    if block.iter().any(|l: &SelexLine| l.scope == "ROW") {
                        return Err(bad("SELEX consensus/reference line must precede rows"));
                    }
                    block.push(datum(
                        line,
                        "GC",
                        None,
                        if tag == "CS" { "SS_cons" } else { "RF" },
                        prefix_end,
                    )?);
                    seen_data = true;
                    last_row = None;
                }
                "SS" | "SA" => {
                    let id = last_row.clone().ok_or_else(|| {
                        bad("SELEX row annotation must immediately follow a sequence")
                    })?;
                    block.push(datum(line, "GR", Some(id), tag, prefix_end)?);
                }
                "SQ" => {
                    if seen_data {
                        return Err(bad("SELEX SQ header must precede alignment blocks"));
                    }
                    let mut fields = &machine[end..];
                    let id = word(&mut fields).ok_or_else(|| bad("missing SQ row name"))?;
                    if !token(id) || id.len() >= 32 {
                        return Err(unsupported(
                            "SELEX row names must be <32 non-whitespace bytes",
                        ));
                    }
                    let remaining = fields.trim_start();
                    let mut check = remaining;
                    let weight = word(&mut check).ok_or_else(|| bad("missing SQ weight"))?;
                    if weight != "-"
                        && !weight
                            .parse::<f64>()
                            .is_ok_and(|v| v.is_finite() && v > 0.0)
                    {
                        return Err(bad("invalid SQ source weight"));
                    }
                    for _ in 0..3 {
                        word(&mut check).ok_or_else(|| bad("missing SQ source fields"))?;
                    }
                    if check.trim().is_empty() {
                        return Err(bad("SQ description field missing (use - when unknown)"));
                    }
                    headers.push((
                        id.to_string(),
                        remaining.to_string(),
                        check.trim().to_string(),
                    ));
                    if headers.len() > MAX_ALIGNMENT_ROWS {
                        return Err(Error::TooLarge);
                    }
                }
                "ID" | "AC" | "DE" | "AU" | "GA" | "TC" | "NC" => {
                    if seen_data {
                        return Err(bad("SELEX header after data"));
                    }
                    if !global_headers.insert(tag.to_string()) {
                        return Err(bad("repeated unique SELEX header field"));
                    }
                    let value = machine[end..].trim();
                    if value.is_empty() {
                        return Err(bad("empty SELEX header field"));
                    }
                    let value = selex_header_value(tag, &value.replace('\t', " "))?;
                    annotation(&mut meta, "GF", None, tag, value)?;
                }
                _ => return Err(unsupported("unsupported SELEX machine-comment class")),
            }
            continue;
        }
        if line.starts_with(['#', '%']) {
            annotation(
                &mut meta,
                "COMMENT",
                None,
                "",
                if line.starts_with('#') {
                    line.replace('\t', " ")
                } else {
                    format!("#SELEX {line}").replace('\t', " ")
                },
            )?;
            continue;
        }
        if line.starts_with(char::is_whitespace) {
            return Err(bad("SELEX row names must start in column1"));
        }
        let end = line
            .find(char::is_whitespace)
            .ok_or_else(|| bad("SELEX row needs an explicit data separator"))?;
        let id = &line[..end];
        if !token(id) || id.len() >= 32 {
            return Err(unsupported(
                "SELEX row names must be <32 non-whitespace bytes",
            ));
        }
        if block
            .iter()
            .any(|l| l.scope == "ROW" && l.row.as_deref() == Some(id))
        {
            state.flush(&mut block, alphabet)?;
        }
        block.push(datum(line, "ROW", Some(id.into()), "", end)?);
        last_row = Some(id.to_string());
        seen_data = true;
        if block.len() > MAX_ALIGNMENT_ROWS + MAX_ALIGNMENT_ANNOTATIONS {
            return Err(Error::TooLarge);
        }
    }
    state.flush(&mut block, alphabet)?;
    if state.rows.is_empty() {
        return Err(bad("no SELEX rows"));
    }
    if !headers.is_empty() {
        if headers
            .iter()
            .map(|h| h.0.as_str())
            .ne(state.rows.iter().map(|r| r.0.as_str()))
        {
            return Err(bad(
                "SELEX SQ headers must cover every row once in exact order",
            ));
        }
        for ((id, value, description), row) in headers.into_iter().zip(state.rows.iter_mut()) {
            annotation(&mut meta, "GS", Some(&id), "SELEX_SQ", value)?;
            row.1 = if description == "-" {
                String::new()
            } else {
                description
            };
        }
    }
    for ((scope, row, tag), value) in state.columns {
        annotation(&mut meta, &scope, row.as_deref(), &tag, value)?;
    }
    for (row, mask) in state.rows.iter().zip(state.case) {
        if mask.contains('l') {
            annotation(&mut meta, "GR", Some(&row.0), "SELEX_CASE", mask)?;
        }
    }
    annotation(
        &mut meta,
        "GF",
        None,
        "SELEX_LAYOUT",
        "Internal spaces are alignment gaps; missing trailing cells are padded. Exterior all-blank columns require explicit gap glyphs. Source letter case is retained as SELEX_CASE, not inferred confidence or biological annotation.",
    )?;
    Ok((state.rows, meta))
}
pub fn import_alignment(
    bytes: &[u8],
    format: AlignmentFormat,
    alphabet: Alphabet,
    name: &str,
) -> Result<AlignmentDocument, Error> {
    let text = text(bytes)?;
    let (rows, metadata) = match format {
        AlignmentFormat::Pir => parse_pir(text, alphabet)?,
        AlignmentFormat::Gde => parse_gde(text, alphabet)?,
        AlignmentFormat::Selex => parse_selex(text, alphabet)?,
        _ => return Err(unsupported("not a PIR/GDE/SELEX exchange profile")),
    };
    base::import_parsed(bytes, format, alphabet, name, rows, metadata)
}

fn row_value<'a>(
    meta: &'a AlignmentMetadata,
    id: &str,
    tag: &str,
) -> Result<Option<&'a str>, Error> {
    let mut values = meta
        .annotations
        .iter()
        .filter(|a| a.scope == "GS" && a.row_identifier.as_deref() == Some(id) && a.tag == tag);
    let first = values.next().map(|a| a.value.as_str());
    if values.next().is_some() {
        return Err(bad("ambiguous repeated structured row metadata"));
    }
    Ok(first)
}
struct Output(String);
impl Output {
    fn line(&mut self, value: &str) -> Result<(), Error> {
        if value.len() > MAX_ALIGNMENT_LINE
            || self
                .0
                .len()
                .checked_add(value.len() + 1)
                .is_none_or(|n| n > MAX_ALIGNMENT_BYTES)
        {
            return Err(Error::TooLarge);
        }
        self.0.push_str(value);
        self.0.push('\n');
        Ok(())
    }
}
fn handled(annotation: &AlignmentAnnotation, format: AlignmentFormat) -> bool {
    match format {
        AlignmentFormat::Pir => {
            annotation.scope == "GS" && matches!(annotation.tag.as_str(), "PIR_TYPE" | "PIR_TEXT")
        }
        AlignmentFormat::Gde => {
            annotation.scope == "GS"
                && matches!(annotation.tag.as_str(), "GDE_OFFSET" | "GDE_SOURCE_WIDTH")
                || annotation.scope == "GF" && annotation.tag == "GDE_LAYOUT"
        }
        AlignmentFormat::Selex => {
            annotation.scope == "COMMENT"
                || annotation.scope == "GS" && annotation.tag == "SELEX_SQ"
                || annotation.scope == "GF"
                    && matches!(
                        annotation.tag.as_str(),
                        "ID" | "AC" | "DE" | "AU" | "GA" | "TC" | "NC" | "SELEX_LAYOUT"
                    )
                || annotation.scope == "GC" && matches!(annotation.tag.as_str(), "RF" | "SS_cons")
                || annotation.scope == "GR"
                    && matches!(annotation.tag.as_str(), "SS" | "SA" | "SELEX_CASE")
        }
        _ => false,
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
/// Dispatcher hook for an already scientifically computed alignment; no imported
/// original is fabricated when the caller supplies validated current rows.
pub(crate) fn export_rows(
    rows: &[AlignmentRow],
    meta: &AlignmentMetadata,
    alphabet: Alphabet,
    format: AlignmentFormat,
    options: &ExportOptions,
) -> Result<AlignmentExport, Error> {
    base::validate_export_content(alphabet, rows, meta)?;
    if options.block_width == 0 || options.block_width > 1000 {
        return Err(bad("exchange block width must be1..1000"));
    }
    if !matches!(
        format,
        AlignmentFormat::Pir | AlignmentFormat::Gde | AlignmentFormat::Selex
    ) {
        return Err(unsupported("not a PIR/GDE/SELEX exchange profile"));
    }
    let mut out = Output(String::new());
    let mut losses=vec!["Canonical exchange changes original layout/case/gap glyphs and does not encode the native alignment document name. Exact source bytes remain available separately.".into()];
    if meta.clustal_consensus.is_some() {
        losses.push("Imported CLUSTAL consensus is not emitted by this exchange profile.".into());
    }
    if meta.annotations.iter().any(|a| !handled(a, format)) {
        losses.push("Annotations with no declared representation in this profile remain only in the native/original alignment; no biological interpretation is substituted.".into());
    }
    if meta.source_header.is_some() {
        losses.push("Original format-specific header text is not reconstructed.".into());
    }
    match format {
        AlignmentFormat::Pir => {
            for row in rows {
                if !token(&row.identifier)
                    || row.identifier.contains(';')
                    || row.identifier.len() + 4 > 500
                {
                    return Err(unsupported(
                        "PIR identifier is not representable without renaming",
                    ));
                }
                if row.description.len() > 500 {
                    return Err(unsupported(
                        "PIR description exceeds its single physical title line",
                    ));
                }
                if row.aligned.contains('*') {
                    return Err(unsupported(
                        "protein stop symbols cannot be confused with PIR record terminators",
                    ));
                }
                let code =
                    row_value(meta, &row.identifier, "PIR_TYPE")?.unwrap_or(match alphabet {
                        Alphabet::Protein => "F1",
                        Alphabet::Dna => "DL",
                        Alphabet::Rna => "RL",
                    });
                type_code(code, alphabet)?;
                if alphabet == Alphabet::Protein
                    && row_value(meta, &row.identifier, "PIR_TYPE")?.is_none()
                {
                    losses.push(format!("{}: no completeness declaration; emitted F1 fragment profile rather than claiming a complete protein.",row.identifier));
                }
                out.line(&format!(">{code};{}", row.identifier))?;
                out.line(&row.description)?;
                let width = options.block_width.min(500);
                for line in row.aligned.as_bytes().chunks(width) {
                    out.line(std::str::from_utf8(line).expect("validatedASCII"))?;
                }
                out.line("*")?;
                for annotation in meta.annotations.iter().filter(|a| {
                    a.scope == "GS"
                        && a.row_identifier.as_deref() == Some(&row.identifier)
                        && a.tag == "PIR_TEXT"
                }) {
                    if annotation.value.starts_with('>') || annotation.value.len() > 500 {
                        return Err(unsupported(
                            "PIR source text is not representable without changing its record structure",
                        ));
                    }
                    out.line(&annotation.value)?;
                }
            }
        }
        AlignmentFormat::Gde => {
            if rows.iter().any(|r| !r.description.is_empty()) {
                losses.push("GDE flat has no row-description field; descriptions remain in the native/original alignment.".into());
            }
            losses.push("GDE exports retain the complete common-width trailing gap columns. Original shorter row extents remain source metadata, not inferred biological endpoints.".into());
            for row in rows {
                if !token(&row.identifier)
                    || row.identifier.contains(['(', ')', '#', '%', '@', '"'])
                {
                    return Err(unsupported(
                        "GDE short name cannot be represented without renaming",
                    ));
                }
                let offset = row_value(meta, &row.identifier, "GDE_OFFSET")?
                    .map(str::parse::<usize>)
                    .transpose()
                    .map_err(|_| bad("invalid retained GDE offset"))?
                    .unwrap_or(0);
                if offset >= row.aligned.len()
                    || !row.aligned.as_bytes()[..offset].iter().all(|b| *b == b'-')
                {
                    return Err(unsupported(
                        "retained GDE offset no longer agrees with the current alignment prefix",
                    ));
                }
                let marker = if alphabet == Alphabet::Protein {
                    '%'
                } else {
                    '#'
                };
                out.line(&format!("{marker}{}({offset})", row.identifier))?;
                for line in row.aligned.as_bytes()[offset..].chunks(options.block_width) {
                    out.line(std::str::from_utf8(line).expect("validatedASCII"))?;
                }
            }
        }
        AlignmentFormat::Selex => {
            for row in rows {
                if !token(&row.identifier)
                    || row.identifier.len() >= 32
                    || row.identifier.starts_with(['#', '%'])
                {
                    return Err(unsupported(
                        "SELEX identifiers must be <32 bytes and cannot resemble comment directives",
                    ));
                }
            }
            let mut global = BTreeSet::new();
            for annotation in &meta.annotations {
                if annotation.scope == "GF"
                    && matches!(
                        annotation.tag.as_str(),
                        "ID" | "AC" | "DE" | "AU" | "GA" | "TC" | "NC"
                    )
                {
                    if !global.insert(annotation.tag.clone()) {
                        return Err(unsupported(
                            "repeated global fields cannot be emitted as unique SELEX headers",
                        ));
                    }
                    let value = selex_header_value(&annotation.tag, &annotation.value)?;
                    out.line(&format!("#={} {}", annotation.tag, value))?;
                } else if annotation.scope == "COMMENT" {
                    out.line(&annotation.value)?;
                }
            }
            let headers = rows.iter().any(|r| {
                !r.description.is_empty()
                    || meta.annotations.iter().any(|a| {
                        a.scope == "GS"
                            && a.row_identifier.as_deref() == Some(&r.identifier)
                            && a.tag == "SELEX_SQ"
                    })
            });
            if headers {
                for row in rows {
                    let header = if let Some(value) = row_value(meta, &row.identifier, "SELEX_SQ")?
                    {
                        let mut rest = value;
                        let mut values = Vec::new();
                        for _ in 0..4 {
                            values.push(
                                word(&mut rest)
                                    .ok_or_else(|| bad("incomplete retained SQ fields"))?,
                            );
                        }
                        if values[0] != "-"
                            && !values[0]
                                .parse::<f64>()
                                .is_ok_and(|v| v.is_finite() && v > 0.0)
                        {
                            return Err(bad("invalid retained SQ weight"));
                        }
                        values.join(" ")
                    } else {
                        "- - - 0".into()
                    };
                    out.line(&format!(
                        "#=SQ {} {} {}",
                        row.identifier,
                        header,
                        if row.description.is_empty() {
                            "-"
                        } else {
                            &row.description
                        }
                    ))?;
                }
            }
            let width = rows[0].aligned.len();
            for start in (0..width).step_by(options.block_width) {
                let end = (start + options.block_width).min(width);
                out.line("")?;
                for (tag, source) in [("RF", "RF"), ("CS", "SS_cons")] {
                    if let Some(annotation) = meta
                        .annotations
                        .iter()
                        .find(|a| a.scope == "GC" && a.tag == source)
                    {
                        out.line(&format!(
                            "{:<32}{}",
                            format!("#={tag}"),
                            &annotation.value[start..end]
                        ))?;
                    }
                }
                for row in rows {
                    let mut sequence = row.aligned.as_bytes()[start..end].to_vec();
                    if let Some(mask) = meta.annotations.iter().find(|a| {
                        a.scope == "GR"
                            && a.row_identifier.as_deref() == Some(&row.identifier)
                            && a.tag == "SELEX_CASE"
                    }) {
                        for (residue, case) in sequence
                            .iter_mut()
                            .zip(mask.value.as_bytes()[start..end].iter())
                        {
                            match (*residue, *case) {
                                (b'-', b'.') => (),
                                (b'-', _) | (_, b'.') => {
                                    return Err(bad(
                                        "SELEX case mask disagrees with gap positions",
                                    ));
                                }
                                (_, b'l') => *residue = residue.to_ascii_lowercase(),
                                (_, b'u') => (),
                                _ => return Err(bad("invalid SELEX case mask")),
                            }
                        }
                    }
                    out.line(&format!(
                        "{:<32}{}",
                        row.identifier,
                        std::str::from_utf8(&sequence).expect("validatedASCII")
                    ))?;
                    for tag in ["SS", "SA"] {
                        if let Some(annotation) = meta.annotations.iter().find(|a| {
                            a.scope == "GR"
                                && a.row_identifier.as_deref() == Some(&row.identifier)
                                && a.tag == tag
                        }) {
                            out.line(&format!(
                                "{:<32}{}",
                                format!("#={tag}"),
                                &annotation.value[start..end]
                            ))?;
                        }
                    }
                }
            }
        }
        _ => unreachable!("formatchecked"),
    }
    Ok(AlignmentExport {
        text: out.0,
        losses,
        renamings: BTreeMap::new(),
    })
}
