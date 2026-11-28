//! Bounded alignment exchange. This module does not compute alignments.
use base64::{Engine, engine::general_purpose::STANDARD};
use collagene_model::{Alphabet, Error, Molecule, Sequence, Topology};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};

pub const WORKSPACE_ATTACHMENT_KIND: &str = "alignment-document-v1";
pub const MAX_ALIGNMENT_BYTES: usize = 8 * 1024 * 1024;
pub const MAX_ALIGNMENT_ROWS: usize = 1000;
pub const MAX_ALIGNMENT_CELLS: usize = 8_000_000;
pub const MAX_ALIGNMENT_LINE: usize = 65_536;
pub const MAX_ALIGNMENT_ANNOTATIONS: usize = 20_000;
const MAX_ATTACHMENT_BYTES: usize = 32 * 1024 * 1024;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum AlignmentFormat {
    AlignedFasta,
    Clustal,
    Stockholm,
    PhylipSequential,
    PhylipInterleaved,
    PhylipRelaxedSequential,
    PhylipRelaxedInterleaved,
    Msf,
    Nexus,
    Pir,
    Gde,
    Selex,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct AlignmentRow {
    /// Stable source-derived identity, independent of editable display identifiers.
    pub id: String,
    pub source_row_index: usize,
    pub identifier: String,
    pub description: String,
    pub aligned: String,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct AlignmentAnnotation {
    /// GF, GS, GR, GC, or COMMENT. Repeated free-text GF/GS entries remain ordered.
    pub scope: String,
    pub row_identifier: Option<String>,
    pub tag: String,
    pub value: String,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize, Default)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct AlignmentMetadata {
    pub annotations: Vec<AlignmentAnnotation>,
    /// Imported CLUSTAL symbols, never a computed confidence score.
    pub clustal_consensus: Option<String>,
    pub source_header: Option<String>,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct AlignmentSource {
    pub format: AlignmentFormat,
    /// Explicit import provenance; no alphabet guessing.
    pub alphabet: Alphabet,
    pub import_name: String,
    pub blake3: String,
    pub original_base64: String,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct AlignmentDocument {
    pub schema_version: u32,
    pub name: String,
    pub alphabet: Alphabet,
    pub rows: Vec<AlignmentRow>,
    pub metadata: AlignmentMetadata,
    pub source: AlignmentSource,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct AlignmentExport {
    pub text: String,
    pub losses: Vec<String>,
    /// Source identifier to output identifier; currently always empty (no silent renaming).
    pub renamings: BTreeMap<String, String>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum RowOrientation {
    AsStored,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct UngappedPosition {
    /// Zero-based index in the extracted row, not an external reference coordinate.
    pub position: usize,
    pub orientation: RowOrientation,
}

fn bad(message: impl Into<String>) -> Error {
    Error::Format(message.into())
}
fn token(s: &str) -> bool {
    !s.is_empty()
        && s.is_ascii()
        && !s
            .bytes()
            .any(|c| c.is_ascii_whitespace() || c.is_ascii_control())
}
fn row_id(hash: &str, i: usize) -> String {
    format!("aln-{}-{i}", &hash[..24])
}
fn normalized(s: &str, alphabet: Alphabet, dots: bool) -> Result<String, Error> {
    if s.is_empty() || !s.is_ascii() || s.bytes().any(|b| b.is_ascii_whitespace()) {
        return Err(bad("empty or non-ASCII alignment segment"));
    }
    let value = if dots {
        s.replace('.', "-").to_ascii_uppercase()
    } else {
        s.to_ascii_uppercase()
    };
    let residues: String = value.chars().filter(|c| *c != '-').collect();
    if !residues.is_empty() {
        Sequence::new(alphabet, &residues)?;
    }
    Ok(value)
}

/// Caller must explicitly supply the alphabet, including for ambiguous A/C/G/T proteins.
pub fn import_alignment(
    bytes: &[u8],
    format: AlignmentFormat,
    alphabet: Alphabet,
    name: &str,
) -> Result<AlignmentDocument, Error> {
    if bytes.len() > MAX_ALIGNMENT_BYTES {
        return Err(Error::TooLarge);
    }
    let text = std::str::from_utf8(bytes).map_err(|_| bad("alignment must be UTF-8 ASCII text"))?;
    if !text.is_ascii() || text.lines().any(|s| s.len() > MAX_ALIGNMENT_LINE) {
        return Err(bad("non-ASCII alignment or line exceeds 65536 bytes"));
    }
    if text
        .bytes()
        .any(|b| b.is_ascii_control() && !b"\r\n\t".contains(&b))
    {
        return Err(bad("unsupported control character"));
    }
    let (raw_rows, metadata) = match format {
        AlignmentFormat::AlignedFasta => parse_fasta(text, alphabet)?,
        AlignmentFormat::Clustal => parse_blocks(text, alphabet, false)?,
        AlignmentFormat::Stockholm => parse_blocks(text, alphabet, true)?,
        AlignmentFormat::Pir | AlignmentFormat::Gde | AlignmentFormat::Selex => {
            return crate::alignment_exchange_more::import_alignment(bytes, format, alphabet, name);
        }
        _ => return crate::alignment_exchange::import_alignment(bytes, format, alphabet, name),
    };
    import_parsed(bytes, format, alphabet, name, raw_rows, metadata)
}

/// Shared native provenance constructor for independently parsed exchange rows.
pub(crate) fn import_parsed(
    bytes: &[u8],
    format: AlignmentFormat,
    alphabet: Alphabet,
    name: &str,
    raw_rows: Vec<(String, String, String)>,
    metadata: AlignmentMetadata,
) -> Result<AlignmentDocument, Error> {
    if bytes.len() > MAX_ALIGNMENT_BYTES {
        return Err(Error::TooLarge);
    }
    let hash = blake3::hash(bytes).to_hex().to_string();
    let rows = raw_rows
        .into_iter()
        .enumerate()
        .map(|(i, (identifier, description, aligned))| AlignmentRow {
            id: row_id(&hash, i),
            source_row_index: i,
            identifier,
            description,
            aligned,
        })
        .collect();
    let doc = AlignmentDocument {
        schema_version: 1,
        name: name.into(),
        alphabet,
        rows,
        metadata,
        source: AlignmentSource {
            format,
            alphabet,
            import_name: name.into(),
            blake3: hash,
            original_base64: STANDARD.encode(bytes),
        },
    };
    doc.validate_content()?;
    Ok(doc)
}
type RawRows = Vec<(String, String, String)>;
fn parse_fasta(text: &str, alphabet: Alphabet) -> Result<(RawRows, AlignmentMetadata), Error> {
    let mut rows: RawRows = vec![];
    let mut ids = BTreeSet::new();
    for line in text.lines() {
        if line.trim().is_empty() {
            continue;
        }
        if let Some(header) = line.strip_prefix('>') {
            let mut p = header.splitn(2, char::is_whitespace);
            let id = p.next().unwrap_or("");
            if !token(id) || !ids.insert(id.to_string()) {
                return Err(bad("empty or duplicate FASTA identifier"));
            }
            if rows.len() >= MAX_ALIGNMENT_ROWS {
                return Err(Error::TooLarge);
            }
            rows.push((
                id.into(),
                p.next().unwrap_or("").trim().into(),
                String::new(),
            ));
        } else {
            let row = rows
                .last_mut()
                .ok_or_else(|| bad("FASTA residues before first header"))?;
            row.2.push_str(&normalized(line.trim(), alphabet, false)?);
            if row.2.len() > MAX_ALIGNMENT_CELLS {
                return Err(Error::TooLarge);
            }
        }
    }
    Ok((rows, AlignmentMetadata::default()))
}

#[derive(Default)]
struct Block {
    rows: Vec<(String, String)>,
    gr: Vec<(String, String, String)>,
    gc: Vec<(String, String)>,
    consensus: Option<String>,
    offset: Option<usize>,
}
fn split_fields(line: &str, n: usize) -> Result<Vec<&str>, Error> {
    let fields: Vec<_> = line.split_whitespace().collect();
    if fields.len() != n {
        return Err(bad(format!("expected {n} fields: {line}")));
    }
    Ok(fields)
}
fn annotation_rest<'a>(line: &'a str, class: &str) -> Option<&'a str> {
    let rest = line.strip_prefix(class)?;
    if rest.is_empty() || rest.starts_with(char::is_whitespace) {
        Some(rest.trim_start())
    } else {
        None
    }
}
fn annotation(scope: &str, row: Option<String>, tag: &str, value: &str) -> AlignmentAnnotation {
    AlignmentAnnotation {
        scope: scope.into(),
        row_identifier: row,
        tag: tag.into(),
        value: value.into(),
    }
}
fn finish_block(
    block: &mut Block,
    rows: &mut RawRows,
    meta: &mut AlignmentMetadata,
    first: &mut bool,
    stockholm: bool,
) -> Result<(), Error> {
    if block.rows.is_empty() {
        if !block.gr.is_empty() || !block.gc.is_empty() || block.consensus.is_some() {
            return Err(bad("annotation block has no rows"));
        }
        return Ok(());
    }
    let width = block.rows[0].1.len();
    if width == 0 || block.rows.iter().any(|(_, s)| s.len() != width) {
        return Err(bad("unequal alignment block widths"));
    }
    let mut unique = BTreeSet::new();
    if block.rows.iter().any(|(id, _)| !unique.insert(id)) {
        return Err(bad("duplicate row in block"));
    }
    if *first {
        if block.rows.len() > MAX_ALIGNMENT_ROWS {
            return Err(Error::TooLarge);
        }
        for (id, s) in &block.rows {
            rows.push((id.clone(), String::new(), s.clone()));
        }
    } else {
        if block.rows.len() != rows.len()
            || block
                .rows
                .iter()
                .zip(rows.iter())
                .any(|((id, _), (old, _, _))| id != old)
        {
            return Err(bad(
                "alignment blocks must have identical row identifiers and order",
            ));
        }
        for ((_, s), (_, _, full)) in block.rows.iter().zip(rows.iter_mut()) {
            full.push_str(s);
        }
    }
    if rows.len().saturating_mul(rows[0].2.len()) > MAX_ALIGNMENT_CELLS {
        return Err(Error::TooLarge);
    }
    if stockholm {
        let mut keys = BTreeSet::new();
        for (row, tag, value) in &block.gr {
            if !unique.contains(row)
                || value.len() != width
                || !token(tag)
                || !token(value)
                || !keys.insert(("GR", row.clone(), tag.clone()))
            {
                return Err(bad("invalid or duplicate GR annotation"));
            }
        }
        for (tag, value) in &block.gc {
            if value.len() != width
                || !token(tag)
                || !token(value)
                || !keys.insert(("GC", String::new(), tag.clone()))
            {
                return Err(bad("invalid or duplicate GC annotation"));
            }
        }
        let oldkeys: BTreeSet<_> = meta
            .annotations
            .iter()
            .filter(|a| a.scope == "GR" || a.scope == "GC")
            .map(|a| {
                (
                    a.scope.as_str(),
                    a.row_identifier.clone().unwrap_or_default(),
                    a.tag.clone(),
                )
            })
            .collect();
        if !*first && oldkeys != keys {
            return Err(bad(
                "column annotations must occur exactly once in every block",
            ));
        }
        for (row, tag, value) in &block.gr {
            append_annotation(meta, *first, "GR", Some(row.clone()), tag, value)?;
        }
        for (tag, value) in &block.gc {
            append_annotation(meta, *first, "GC", None, tag, value)?;
        }
    } else {
        let c = block.consensus.clone().unwrap_or_else(|| " ".repeat(width));
        if c.len() != width || !c.bytes().all(|b| b" *:.".contains(&b)) {
            return Err(bad("invalid CLUSTAL consensus"));
        }
        meta.clustal_consensus
            .get_or_insert_with(String::new)
            .push_str(&c);
    }
    *first = false;
    *block = Block::default();
    Ok(())
}
fn append_annotation(
    meta: &mut AlignmentMetadata,
    first: bool,
    scope: &str,
    row: Option<String>,
    tag: &str,
    value: &str,
) -> Result<(), Error> {
    if first {
        meta.annotations.push(annotation(scope, row, tag, value));
    } else {
        meta.annotations
            .iter_mut()
            .find(|a| a.scope == scope && a.row_identifier == row && a.tag == tag)
            .ok_or_else(|| bad("annotation absent in first block"))?
            .value
            .push_str(value);
    }
    Ok(())
}
fn parse_blocks(
    text: &str,
    alphabet: Alphabet,
    stockholm: bool,
) -> Result<(RawRows, AlignmentMetadata), Error> {
    let mut lines = text.lines();
    let header = lines
        .next()
        .ok_or_else(|| bad("missing alignment header"))?;
    if (stockholm && header != "# STOCKHOLM 1.0")
        || (!stockholm
            && !(header == "CLUSTAL"
                || header.starts_with("CLUSTAL ")
                || header.starts_with("CLUSTALW ")))
    {
        return Err(bad("unsupported alignment header/version"));
    }
    let mut meta = AlignmentMetadata {
        source_header: Some(header.into()),
        ..Default::default()
    };
    let mut rows = vec![];
    let mut block = Block::default();
    let mut first = true;
    let mut ended = false;
    let mut annotation_count = 0;
    for line in lines {
        if ended {
            if !line.trim().is_empty() {
                return Err(Error::Unsupported(
                    "multiple Stockholm records or trailing data".into(),
                ));
            }
            continue;
        }
        if line.trim().is_empty() {
            finish_block(&mut block, &mut rows, &mut meta, &mut first, stockholm)?;
            continue;
        }
        if stockholm && line == "//" {
            finish_block(&mut block, &mut rows, &mut meta, &mut first, true)?;
            ended = true;
            continue;
        }
        if stockholm && line.starts_with('#') {
            annotation_count += 1;
            if annotation_count > MAX_ALIGNMENT_ANNOTATIONS {
                return Err(Error::TooLarge);
            }
            if let Some(rest) = annotation_rest(line, "#=GF") {
                let (tag, value) = rest
                    .split_once(char::is_whitespace)
                    .ok_or_else(|| bad("GF requires tag and value"))?;
                if !token(tag) || value.trim().is_empty() {
                    return Err(bad("invalid GF annotation"));
                }
                meta.annotations
                    .push(annotation("GF", None, tag, value.trim_start()));
            } else if let Some(rest) = annotation_rest(line, "#=GS") {
                let rest = rest.trim_start();
                let (row, rest) = rest
                    .split_once(char::is_whitespace)
                    .ok_or_else(|| bad("GS requires row tag value"))?;
                let rest = rest.trim_start();
                let (tag, value) = rest
                    .split_once(char::is_whitespace)
                    .ok_or_else(|| bad("GS requires tag value"))?;
                let value = value.trim_start();
                if !token(row) || !token(tag) || value.is_empty() {
                    return Err(bad("invalid GS annotation"));
                }
                meta.annotations
                    .push(annotation("GS", Some(row.into()), tag, value));
            } else if let Some(rest) = annotation_rest(line, "#=GR") {
                let p = split_fields(rest, 3)?;
                block.gr.push((p[0].into(), p[1].into(), p[2].into()));
            } else if let Some(rest) = annotation_rest(line, "#=GC") {
                let p = split_fields(rest, 2)?;
                block.gc.push((p[0].into(), p[1].into()));
            } else if line.starts_with("#=") {
                return Err(Error::Unsupported(
                    "unknown Stockholm annotation class".into(),
                ));
            } else {
                meta.annotations.push(annotation("COMMENT", None, "", line));
            }
            continue;
        }
        if !stockholm && line.starts_with(char::is_whitespace) {
            if block.rows.is_empty() || block.consensus.is_some() || line.contains('\t') {
                return Err(bad("misplaced/ambiguous CLUSTAL consensus"));
            }
            let offset = block.offset.unwrap();
            let width = block.rows[0].1.len();
            if line.len() < offset
                || !line[..offset].bytes().all(|b| b == b' ')
                || line.len() > offset + width && !line[offset + width..].trim().is_empty()
            {
                return Err(bad("CLUSTAL consensus is not aligned with residues"));
            }
            let mut c = line[offset..].to_string();
            if c.len() > width {
                c.truncate(width);
            }
            while c.len() < width {
                c.push(' ');
            }
            block.consensus = Some(c);
            continue;
        }
        let p: Vec<_> = line.split_whitespace().collect();
        if p.len() < 2 || p.len() > if stockholm { 2 } else { 3 } || !token(p[0]) {
            return Err(bad("invalid alignment row fields"));
        }
        if block.rows.first().is_some_and(|(id, _)| id == p[0]) {
            finish_block(&mut block, &mut rows, &mut meta, &mut first, stockholm)?;
        }
        if block.consensus.is_some() {
            return Err(bad("row follows consensus without block boundary"));
        }
        let aligned = normalized(p[1], alphabet, stockholm)?;
        if !stockholm {
            let after_id = line
                .find(char::is_whitespace)
                .ok_or_else(|| bad("missing residues"))?;
            let offset = after_id
                + line[after_id..]
                    .find(|c: char| !c.is_whitespace())
                    .ok_or_else(|| bad("missing residues"))?;
            if block.offset.is_some_and(|old| old != offset) {
                return Err(bad("CLUSTAL residue columns disagree"));
            }
            block.offset = Some(offset);
            if p.len() == 3 {
                let count = p[2]
                    .parse::<usize>()
                    .map_err(|_| bad("invalid CLUSTAL residue count"))?;
                let previous = rows
                    .iter()
                    .find(|(id, _, _)| id == p[0])
                    .map(|(_, _, s)| s.bytes().filter(|b| *b != b'-').count())
                    .unwrap_or(0);
                if count != previous + aligned.bytes().filter(|b| *b != b'-').count() {
                    return Err(bad("CLUSTAL residue count mismatch"));
                }
            }
        }
        if block.rows.len() >= MAX_ALIGNMENT_ROWS {
            return Err(Error::TooLarge);
        }
        block.rows.push((p[0].into(), aligned));
    }
    if stockholm && !ended {
        return Err(bad("Stockholm record missing // terminator"));
    }
    if !stockholm {
        finish_block(&mut block, &mut rows, &mut meta, &mut first, false)?;
    }
    Ok((rows, meta))
}

pub(crate) fn validate_export_content(
    alphabet: Alphabet,
    rows: &[AlignmentRow],
    metadata: &AlignmentMetadata,
) -> Result<(), Error> {
    if rows.is_empty() || rows.len() > MAX_ALIGNMENT_ROWS {
        return Err(bad("alignment needs 1..1000 rows"));
    }
    let width = rows[0].aligned.len();
    if width == 0 || width.saturating_mul(rows.len()) > MAX_ALIGNMENT_CELLS {
        return Err(Error::TooLarge);
    }
    let mut ids = BTreeSet::new();
    let mut labels = BTreeSet::new();
    let mut origins = BTreeSet::new();
    for row in rows {
        if row.identifier.trim().is_empty()
            || !row.identifier.is_ascii()
            || row.identifier.chars().any(char::is_control)
            || row.identifier.len() > 4096
            || !ids.insert(&row.id)
            || !labels.insert(&row.identifier)
            || !origins.insert(row.source_row_index)
        {
            return Err(bad("invalid/duplicate row identity"));
        }
        if row.description.len() > MAX_ALIGNMENT_LINE
            || row.description.chars().any(char::is_control)
            || !row.description.is_ascii()
        {
            return Err(bad("invalid row description"));
        }
        if row.aligned.len() != width || normalized(&row.aligned, alphabet, false)? != row.aligned {
            return Err(bad(
                "rows must have equal width and canonical uppercase residues/gaps",
            ));
        }
    }
    if metadata.annotations.len() > MAX_ALIGNMENT_ANNOTATIONS {
        return Err(Error::TooLarge);
    }
    if metadata
        .annotations
        .iter()
        .map(|a| a.value.len() + a.tag.len() + a.row_identifier.as_ref().map_or(0, String::len))
        .sum::<usize>()
        > MAX_ALIGNMENT_BYTES
    {
        return Err(Error::TooLarge);
    }
    let mut column_keys = BTreeSet::new();
    for a in &metadata.annotations {
        if a.value.len() > MAX_ALIGNMENT_CELLS
            || !a.value.is_ascii()
            || a.value.bytes().any(|b| b.is_ascii_control())
        {
            return Err(bad("invalid annotation text"));
        }
        if a.row_identifier
            .as_ref()
            .is_some_and(|s| !labels.contains(s))
        {
            return Err(bad("annotation refers to absent row"));
        }
        if (a.scope == "GF" || a.scope == "GS") && a.value.trim().is_empty() {
            return Err(bad("free-text annotations require a nonempty value"));
        }
        match a.scope.as_str() {
            "GF" | "GC" if a.row_identifier.is_none() && token(&a.tag) => {}
            "GS" | "GR" if a.row_identifier.is_some() && token(&a.tag) => {}
            "COMMENT"
                if a.row_identifier.is_none()
                    && a.tag.is_empty()
                    && a.value.starts_with('#')
                    && !a.value.starts_with("#=") => {}
            _ => return Err(bad("invalid annotation scope/tag")),
        }
        if a.scope == "GR" || a.scope == "GC" {
            if a.value.len() != width
                || !token(&a.value)
                || !column_keys.insert((&a.scope, &a.row_identifier, &a.tag))
            {
                return Err(bad("column annotation width or identity mismatch"));
            }
        }
    }
    if let Some(c) = &metadata.clustal_consensus {
        if c.len() != width || !c.bytes().all(|b| b" *:.".contains(&b)) {
            return Err(bad("invalid consensus width/symbol"));
        }
    }
    if metadata.source_header.as_ref().is_some_and(|s| {
        s.len() > MAX_ALIGNMENT_LINE || s.chars().any(char::is_control) || !s.is_ascii()
    }) {
        return Err(bad("invalid source header"));
    }
    Ok(())
}

impl AlignmentDocument {
    fn validate_content(&self) -> Result<(), Error> {
        if self.schema_version != 1 {
            return Err(Error::Unsupported("alignment schema version".into()));
        }
        if self.name.trim().is_empty()
            || self.name.len() > 4096
            || self.name.chars().any(char::is_control)
        {
            return Err(bad("invalid alignment name"));
        }
        validate_export_content(self.alphabet, &self.rows, &self.metadata)
    }
    /// Validate structure, immutable source hash and stable row derivation. Edited content
    /// remains usable but cannot claim that the original source describes current content.
    pub fn validate(&self) -> Result<(), Error> {
        self.validate_content()?;
        let original = self.source_document()?;
        if self.alphabet != self.source.alphabet {
            return Err(bad(
                "current alphabet disagrees with explicit source provenance",
            ));
        }
        for row in &self.rows {
            let old = original
                .rows
                .get(row.source_row_index)
                .ok_or_else(|| bad("unknown source row index"))?;
            if row.id != old.id {
                return Err(bad("stable row identity does not derive from source"));
            }
        }
        let sequence_changed = self.rows.len() != original.rows.len()
            || self.rows.iter().any(|row| {
                original.rows.get(row.source_row_index).is_none_or(|old| {
                    old.aligned != row.aligned || old.identifier != row.identifier
                })
            });
        if sequence_changed
            && (self.metadata.clustal_consensus.is_some()
                || self
                    .metadata
                    .annotations
                    .iter()
                    .any(|a| a.scope == "GR" || a.scope == "GC"))
        {
            return Err(bad(
                "edited alignment must clear original per-column annotations/consensus before reuse",
            ));
        }
        Ok(())
    }
    fn source_document(&self) -> Result<Self, Error> {
        let bytes = self.original_bytes()?;
        import_alignment(
            &bytes,
            self.source.format,
            self.source.alphabet,
            &self.source.import_name,
        )
    }
    /// Exact imported bytes, even after current content has been edited.
    pub fn original_bytes(&self) -> Result<Vec<u8>, Error> {
        if self.source.original_base64.len() > MAX_ALIGNMENT_BYTES.div_ceil(3) * 4 {
            return Err(Error::TooLarge);
        }
        let bytes = STANDARD
            .decode(&self.source.original_base64)
            .map_err(|_| bad("invalid source base64"))?;
        if bytes.len() > MAX_ALIGNMENT_BYTES {
            return Err(Error::TooLarge);
        }
        if STANDARD.encode(&bytes) != self.source.original_base64
            || blake3::hash(&bytes).to_hex().as_str() != self.source.blake3
        {
            return Err(bad("source hash/base64 mismatch"));
        }
        Ok(bytes)
    }
    pub fn source_is_current(&self) -> Result<bool, Error> {
        self.validate()?;
        Ok(*self == self.source_document()?)
    }
    pub fn column_count(&self) -> usize {
        self.rows.first().map_or(0, |r| r.aligned.len())
    }
    pub fn column_to_ungapped(
        &self,
        row_id: &str,
        column: usize,
    ) -> Result<Option<UngappedPosition>, Error> {
        self.validate()?;
        let row = self
            .rows
            .iter()
            .find(|r| r.id == row_id)
            .ok_or_else(|| bad("unknown alignment row"))?;
        if column >= row.aligned.len() {
            return Err(Error::Coordinates("alignment column outside row".into()));
        }
        Ok(
            (row.aligned.as_bytes()[column] != b'-').then(|| UngappedPosition {
                position: row.aligned.as_bytes()[..column]
                    .iter()
                    .filter(|b| **b != b'-')
                    .count(),
                orientation: RowOrientation::AsStored,
            }),
        )
    }
    /// Build a complete row map once for viewers; linear in aligned row length.
    /// Unlike repeated single-column queries this validates source provenance only once.
    pub fn row_position_map(&self, row_id: &str) -> Result<Vec<Option<UngappedPosition>>, Error> {
        self.validate()?;
        let row = self
            .rows
            .iter()
            .find(|r| r.id == row_id)
            .ok_or_else(|| bad("unknown alignment row"))?;
        let mut position = 0;
        Ok(row
            .aligned
            .bytes()
            .map(|b| {
                if b == b'-' {
                    None
                } else {
                    let mapped = UngappedPosition {
                        position,
                        orientation: RowOrientation::AsStored,
                    };
                    position += 1;
                    Some(mapped)
                }
            })
            .collect())
    }
    /// Linear ordinary molecules. No external coordinates, confidence or annotations
    /// are copied from an alignment; all-gap rows cannot become empty molecules.
    pub fn extract_rows(&self, row_ids: &[String]) -> Result<Vec<Molecule>, Error> {
        self.validate()?;
        if row_ids.len() > self.rows.len() {
            return Err(bad("too many selected rows"));
        }
        let mut seen = BTreeSet::new();
        let mut result = vec![];
        let current_hash = blake3::hash(&serde_json::to_vec(self).map_err(|e| bad(e.to_string()))?)
            .to_hex()
            .to_string();
        for id in row_ids {
            if !seen.insert(id) {
                return Err(bad("duplicate selected row"));
            }
            let row = self
                .rows
                .iter()
                .find(|r| &r.id == id)
                .ok_or_else(|| bad("unknown selected row"))?;
            let residues = row.aligned.replace('-', "");
            let molecule = Molecule {
                name: row.identifier.clone(),
                sequence: Sequence::new(self.alphabet, &residues)?,
                topology: Topology::Linear,
                features: vec![],
                metadata: BTreeMap::from([
                    ("alignment_row_id".into(), row.id.clone()),
                    ("alignment_source_blake3".into(), self.source.blake3.clone()),
                    ("alignment_current_blake3".into(), current_hash.clone()),
                    ("alignment_mapping".into(), "ungapped row as stored; no external reference coordinates or alignment confidence".into()),
                ]),
                original_genbank: None,
            };
            molecule.validate()?;
            result.push(molecule);
        }
        Ok(result)
    }
}

pub fn save_alignment(doc: &AlignmentDocument) -> Result<Vec<u8>, Error> {
    doc.validate()?;
    let bytes = serde_json::to_vec(doc).map_err(|e| bad(e.to_string()))?;
    if bytes.len() > MAX_ATTACHMENT_BYTES {
        return Err(Error::TooLarge);
    }
    Ok(bytes)
}
pub fn load_alignment(bytes: &[u8]) -> Result<AlignmentDocument, Error> {
    if bytes.len() > MAX_ATTACHMENT_BYTES {
        return Err(Error::TooLarge);
    }
    let raw: serde_json::Value = serde_json::from_slice(bytes).map_err(|e| bad(e.to_string()))?;
    let doc: AlignmentDocument = serde_json::from_slice(bytes).map_err(|e| bad(e.to_string()))?;
    if serde_json::to_value(&doc).map_err(|e| bad(e.to_string()))? != raw {
        return Err(bad("alignment JSON is not canonical typed content"));
    }
    doc.validate()?;
    Ok(doc)
}

/// Canonical sequence-preserving export. Losses describe metadata/formatting changes;
/// exact evidence export is `original_bytes`, not this function.
pub fn export_alignment(
    doc: &AlignmentDocument,
    format: AlignmentFormat,
) -> Result<AlignmentExport, Error> {
    doc.validate()?;
    let mut result = export_rows(&doc.rows, &doc.metadata, doc.alphabet, format)?;
    if matches!(
        format,
        AlignmentFormat::AlignedFasta | AlignmentFormat::Clustal | AlignmentFormat::Stockholm
    ) {
        result.losses.insert(0, "Canonical export does not preserve original byte layout, letter case, Stockholm dot gaps, header text, or document name; use original_bytes for exact evidence.".into());
    }
    Ok(result)
}

/// Export computed rows without inventing an imported original-file record.
/// Row identities, canonical residues/gaps and metadata are checked independently.
pub(crate) fn export_rows(
    rows: &[AlignmentRow],
    metadata: &AlignmentMetadata,
    alphabet: Alphabet,
    format: AlignmentFormat,
) -> Result<AlignmentExport, Error> {
    if matches!(
        format,
        AlignmentFormat::Pir | AlignmentFormat::Gde | AlignmentFormat::Selex
    ) {
        return crate::alignment_exchange_more::export_rows(
            rows,
            metadata,
            alphabet,
            format,
            &Default::default(),
        );
    }
    if !matches!(
        format,
        AlignmentFormat::AlignedFasta | AlignmentFormat::Clustal | AlignmentFormat::Stockholm
    ) {
        return crate::alignment_exchange::export_rows(
            rows,
            metadata,
            alphabet,
            format,
            &Default::default(),
        );
    }
    validate_export_content(alphabet, rows, metadata)?;
    if rows.iter().any(|row| !token(&row.identifier)) {
        return Err(bad(
            "this unquoted alignment format cannot preserve whitespace in row identifiers; rename explicitly or choose NEXUS/strict PHYLIP",
        ));
    }
    let mut text = String::new();
    let mut losses = Vec::new();
    match format {
        AlignmentFormat::AlignedFasta => {
            for row in rows {
                text.push('>');
                text.push_str(&row.identifier);
                if !row.description.is_empty() {
                    text.push(' ');
                    text.push_str(&row.description);
                }
                text.push('\n');
                for chunk in row.aligned.as_bytes().chunks(80) {
                    text.push_str(std::str::from_utf8(chunk).unwrap());
                    text.push('\n');
                }
            }
            if !metadata.annotations.is_empty() {
                losses.push("FASTA omits all Stockholm annotations/comments.".into());
            }
            if metadata.clustal_consensus.is_some() {
                losses.push("FASTA omits imported CLUSTAL consensus.".into());
            }
        }
        AlignmentFormat::Clustal => {
            text.push_str("CLUSTAL Collagene alignment exchange\n\n");
            let pad = rows.iter().map(|r| r.identifier.len()).max().unwrap() + 2;
            if pad + 80 > MAX_ALIGNMENT_LINE {
                return Err(Error::TooLarge);
            }
            for start in (0..rows[0].aligned.len()).step_by(80) {
                let end = (start + 80).min(rows[0].aligned.len());
                for row in rows {
                    text.push_str(&format!(
                        "{:<pad$}{}\n",
                        row.identifier,
                        &row.aligned[start..end]
                    ));
                }
                if let Some(c) = &metadata.clustal_consensus {
                    text.push_str(&" ".repeat(pad));
                    text.push_str(&c[start..end]);
                    text.push('\n');
                }
                text.push('\n');
            }
            if !metadata.annotations.is_empty() {
                losses.push("CLUSTAL omits all Stockholm annotations/comments.".into());
            }
            if rows.iter().any(|r| !r.description.is_empty()) {
                losses.push("CLUSTAL omits FASTA row descriptions.".into());
            }
        }
        AlignmentFormat::Stockholm => {
            if rows
                .iter()
                .any(|r| r.identifier.starts_with('#') || r.identifier == "//")
            {
                return Err(bad(
                    "Stockholm identifier conflicts with reserved record/annotation syntax; explicitly rename before export",
                ));
            }
            text.push_str("# STOCKHOLM 1.0\n");
            for a in &metadata.annotations {
                match a.scope.as_str() {
                    "GF" => text.push_str(&format!("#=GF {} {}\n", a.tag, a.value)),
                    "GS" => text.push_str(&format!(
                        "#=GS {} {} {}\n",
                        a.row_identifier.as_ref().unwrap(),
                        a.tag,
                        a.value
                    )),
                    "COMMENT" => {
                        text.push_str(&a.value);
                        text.push('\n');
                    }
                    _ => {}
                }
            }
            for start in (0..rows[0].aligned.len()).step_by(80) {
                let end = (start + 80).min(rows[0].aligned.len());
                for row in rows {
                    text.push_str(&format!(
                        "{} {}\n",
                        row.identifier,
                        &row.aligned[start..end]
                    ));
                    for a in metadata.annotations.iter().filter(|a| {
                        a.scope == "GR" && a.row_identifier.as_deref() == Some(&row.identifier)
                    }) {
                        text.push_str(&format!(
                            "#=GR {} {} {}\n",
                            row.identifier,
                            a.tag,
                            &a.value[start..end]
                        ));
                    }
                }
                for a in metadata.annotations.iter().filter(|a| a.scope == "GC") {
                    text.push_str(&format!("#=GC {} {}\n", a.tag, &a.value[start..end]));
                }
                text.push('\n');
            }
            text.push_str("//\n");
            if rows.iter().any(|r| !r.description.is_empty()) {
                losses.push("Stockholm omits FASTA row descriptions.".into());
            }
            if metadata.clustal_consensus.is_some() {
                losses.push("Stockholm omits imported CLUSTAL consensus; no biological annotation is inferred.".into());
            }
        }
        _ => return Err(bad("unsupported alignment output profile")),
    }
    if text.len() > MAX_ALIGNMENT_BYTES || text.lines().any(|line| line.len() > MAX_ALIGNMENT_LINE)
    {
        return Err(Error::TooLarge);
    }
    Ok(AlignmentExport {
        text,
        losses,
        renamings: BTreeMap::new(),
    })
}
