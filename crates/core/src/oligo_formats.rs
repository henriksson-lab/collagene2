//! Independently implemented typed oligo-table and explicit-encoding FASTQ adapters.
//! Source records are evidence, never authority for physical chemistry or verification.
use crate::{
    analysis::Primer,
    formats::{ExportReport, ImportReport, MAX_IMPORT_BYTES},
};
use base64::{Engine, engine::general_purpose::STANDARD};
use collagene_model::{Error, Feature, Location, Molecule, Sequence, Topology};
use serde::{Deserialize, Serialize};
use std::{
    collections::{BTreeMap, BTreeSet},
    io,
};

pub const MAX_TABLE_BYTES: usize = 8 * 1024 * 1024;
pub const MAX_TABLE_ROWS: usize = 10_000;
pub const MAX_TABLE_COLUMNS: usize = 64;
const MAX_CELL: usize = 8192;
const MAX_FASTQ_RECORDS: usize = 10_000;
pub const FASTQ_ORIGINAL_KEY: &str = "fastq:original_base64";

fn invalid(s: impl Into<String>) -> Error {
    Error::Format(s.into())
}
fn unsupported(s: impl Into<String>) -> Error {
    Error::Unsupported(s.into())
}
struct SizeLimit(usize);
impl io::Write for SizeLimit {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        self.0 = self
            .0
            .checked_add(bytes.len())
            .filter(|n| *n <= MAX_IMPORT_BYTES)
            .ok_or_else(|| io::Error::other("format output limit"))?;
        Ok(bytes.len())
    }
    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}
fn bounded_json(value: &impl Serialize) -> Result<(), Error> {
    serde_json::to_writer(&mut SizeLimit(0), value).map_err(|_| Error::TooLarge)
}
fn decode_attachment(raw: &str) -> Result<Vec<u8>, Error> {
    if raw.len() > MAX_IMPORT_BYTES.div_ceil(3) * 4 {
        return Err(Error::TooLarge);
    }
    let bytes = STANDARD
        .decode(raw)
        .map_err(|_| invalid("invalid source attachment base64"))?;
    if bytes.len() > MAX_IMPORT_BYTES {
        return Err(Error::TooLarge);
    }
    Ok(bytes)
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum Delimiter {
    Comma,
    Tab,
}
impl Delimiter {
    fn byte(self) -> u8 {
        match self {
            Self::Comma => b',',
            Self::Tab => b'\t',
        }
    }
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub struct ColumnMapping {
    pub name: usize,
    pub ordered_sequence: Option<usize>,
    pub annealing_sequence: Option<usize>,
    pub tail: Option<usize>,
    pub modification: Option<usize>,
    pub amount: Option<usize>,
    pub unit: Option<usize>,
}
impl Default for ColumnMapping {
    fn default() -> Self {
        Self {
            name: 0,
            ordered_sequence: Some(1),
            annealing_sequence: None,
            tail: None,
            modification: None,
            amount: None,
            unit: None,
        }
    }
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub struct OligoImportOptions {
    pub delimiter: Delimiter,
    pub has_header: bool,
    pub columns: ColumnMapping,
}
impl Default for OligoImportOptions {
    fn default() -> Self {
        Self {
            delimiter: Delimiter::Comma,
            has_header: true,
            columns: Default::default(),
        }
    }
}
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TablePreview {
    pub headers: Vec<String>,
    pub rows: Vec<Vec<String>>,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OligoQuantity {
    pub raw_value: String,
    pub raw_unit: String,
    /// Supported amount units normalize to mol; supported concentration units to mol/litre.
    pub dimension: Option<String>,
    pub si_value: Option<f64>,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OligoRecord {
    /// Zero-based data-record index, not a physical line number (quoted rows may span lines).
    pub row_index: usize,
    pub row_id: String,
    pub name: String,
    pub original_fields: Vec<String>,
    pub original_record_base64: String,
    pub original_record_blake3: String,
    pub raw_ordered_sequence: Option<String>,
    pub raw_annealing_sequence: Option<String>,
    pub raw_tail: Option<String>,
    pub modification: Option<String>,
    pub quantity: Option<OligoQuantity>,
    /// None means unknown modification/notation or inconsistent fields are retained as evidence only.
    pub primer: Option<Primer>,
    pub diagnostics: Vec<String>,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OligoTable {
    pub schema_version: u16,
    pub source_blake3: String,
    pub original_base64: String,
    pub options: OligoImportOptions,
    pub headers: Vec<String>,
    pub records: Vec<OligoRecord>,
    pub warnings: Vec<String>,
}

// csv intentionally accepts many dialects. This profile rejects malformed quotes rather
// than silently repairing a changed sequence cell; LF/CRLF and quoted newlines are supported.
fn strict_quotes(text: &str, delimiter: u8) -> Result<(), Error> {
    let bytes = text.strip_prefix('\u{feff}').unwrap_or(text).as_bytes();
    let mut state = 0u8;
    let mut i = 0;
    while i < bytes.len() {
        let b = bytes[i];
        if b < 32 && !matches!(b, b'\t' | b'\r' | b'\n') {
            return Err(invalid("control character in oligo table"));
        }
        if state != 2 && b == b'\r' {
            if bytes.get(i + 1) != Some(&b'\n') {
                return Err(invalid("bare CR outside a quoted table cell"));
            }
            state = 0;
            i += 2;
            continue;
        }
        state = match (state, b) {
            (0, b'"') => 2,
            (0, b'\n') => 0,
            (0, b) if b == delimiter => 0,
            (0, _) => 1,
            (1, b'"') => return Err(invalid("quote inside unquoted table cell")),
            (1, b'\n') => 0,
            (1, b) if b == delimiter => 0,
            (1, _) => 1,
            (2, b'"') => 3,
            (2, _) => 2,
            (3, b'"') => 2,
            (3, b'\n') => 0,
            (3, b) if b == delimiter => 0,
            (3, _) => return Err(invalid("text after closing table quote")),
            _ => unreachable!(),
        };
        i += 1;
    }
    if state == 2 {
        return Err(invalid("unterminated quoted table cell"));
    }
    Ok(())
}
struct ParsedTable {
    headers: Vec<String>,
    rows: Vec<Vec<String>>,
    ranges: Vec<(usize, usize)>,
}
fn parse_table(text: &str, delimiter: Delimiter, has_header: bool) -> Result<ParsedTable, Error> {
    if text.len() > MAX_TABLE_BYTES {
        return Err(Error::TooLarge);
    }
    strict_quotes(text, delimiter.byte())?;
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(delimiter.byte())
        .has_headers(false)
        .flexible(true)
        .from_reader(text.as_bytes());
    let mut rows = Vec::new();
    let mut starts = Vec::new();
    for row in reader.records() {
        let row = row.map_err(|e| invalid(format!("oligo table syntax: {e}")))?;
        if rows.len() >= MAX_TABLE_ROWS + usize::from(has_header)
            || row.len() > MAX_TABLE_COLUMNS
            || row.iter().any(|s| s.len() > MAX_CELL)
        {
            return Err(Error::TooLarge);
        }
        // csv's CRLF reader can report the separator LF as the next record's
        // position. Keep that LF with the preceding raw record, not its successor.
        let mut start = row.position().map_or(0, |p| p.byte() as usize);
        while text
            .as_bytes()
            .get(start)
            .is_some_and(|b| matches!(b, b'\r' | b'\n'))
        {
            start += 1;
        }
        starts.push(start);
        rows.push(row.iter().map(String::from).collect::<Vec<_>>());
    }
    if rows.is_empty() {
        return Err(invalid("oligo table contains no records"));
    }
    let count = rows.iter().map(Vec::len).max().unwrap_or(0);
    let mut ranges = starts
        .iter()
        .enumerate()
        .map(|(i, &a)| (a, *starts.get(i + 1).unwrap_or(&text.len())))
        .collect::<Vec<_>>();
    let headers = if has_header {
        ranges.remove(0);
        rows.remove(0)
    } else {
        (1..=count).map(|i| format!("Column {i}")).collect()
    };
    Ok(ParsedTable {
        headers,
        rows,
        ranges,
    })
}
pub fn inspect_table(
    text: &str,
    delimiter: Delimiter,
    has_header: bool,
) -> Result<TablePreview, Error> {
    let parsed = parse_table(text, delimiter, has_header)?;
    Ok(TablePreview {
        headers: parsed.headers,
        rows: parsed.rows,
    })
}
fn mapping(columns: &ColumnMapping, width: usize) -> Result<(), Error> {
    if columns.ordered_sequence.is_none() && columns.annealing_sequence.is_none() {
        return Err(invalid(
            "map an ordered sequence or annealing-domain column",
        ));
    }
    let mut seen = BTreeSet::new();
    for index in [
        Some(columns.name),
        columns.ordered_sequence,
        columns.annealing_sequence,
        columns.tail,
        columns.modification,
        columns.amount,
        columns.unit,
    ]
    .into_iter()
    .flatten()
    {
        if index >= width || !seen.insert(index) {
            return Err(invalid("mapped columns must exist and have distinct roles"));
        }
    }
    Ok(())
}
fn quantity(
    value: Option<&str>,
    unit: Option<&str>,
    warnings: &mut Vec<String>,
) -> Option<OligoQuantity> {
    let (value, unit) = (value.unwrap_or(""), unit.unwrap_or(""));
    if value.is_empty() && unit.is_empty() {
        return None;
    }
    let scale = match unit.trim() {
        "mol" => Some(("amount", 1.)),
        "mmol" => Some(("amount", 1e-3)),
        "umol" | "µmol" | "μmol" => Some(("amount", 1e-6)),
        "nmol" => Some(("amount", 1e-9)),
        "pmol" => Some(("amount", 1e-12)),
        "M" => Some(("concentration", 1.)),
        "mM" => Some(("concentration", 1e-3)),
        "uM" | "µM" | "μM" => Some(("concentration", 1e-6)),
        "nM" => Some(("concentration", 1e-9)),
        "pM" => Some(("concentration", 1e-12)),
        _ => None,
    };
    let parsed = value
        .trim()
        .parse::<f64>()
        .ok()
        .filter(|v| v.is_finite() && *v > 0.);
    let converted = parsed.zip(scale).and_then(|(v, (kind, s))| {
        let v = v * s;
        (v.is_finite() && v > 0.).then_some((kind, v))
    });
    if converted.is_none() {
        warnings.push("Amount/unit retained literally; unsupported, absent or non-positive units are not inferred.".into());
    }
    Some(OligoQuantity {
        raw_value: value.into(),
        raw_unit: unit.into(),
        dimension: converted.map(|v| v.0.into()),
        si_value: converted.map(|v| v.1),
    })
}
pub fn import_oligos(text: &str, options: &OligoImportOptions) -> Result<OligoTable, Error> {
    let parsed = parse_table(text, options.delimiter, options.has_header)?;
    mapping(
        &options.columns,
        parsed
            .headers
            .len()
            .max(parsed.rows.iter().map(Vec::len).max().unwrap_or(0)),
    )?;
    let source_blake3 = blake3::hash(text.as_bytes()).to_hex().to_string();
    let mut records = Vec::new();
    let mut names = BTreeSet::new();
    let mut warnings=vec!["Mapped fields are explicit import choices. Modified/unsupported oligos remain collection evidence and cannot become unmodified DNA primers.".into()];
    for (row_index, fields) in parsed.rows.into_iter().enumerate() {
        let mut diagnostics = Vec::new();
        let get = |column: Option<usize>| column.and_then(|i| fields.get(i).map(String::as_str));
        let name = fields
            .get(options.columns.name)
            .cloned()
            .unwrap_or_default();
        if !names.insert(name.clone()) {
            diagnostics
                .push("Duplicate name retained with a distinct source-row identifier.".into());
        }
        let raw_ordered_sequence = get(options.columns.ordered_sequence).map(String::from);
        let raw_annealing_sequence = get(options.columns.annealing_sequence).map(String::from);
        let raw_tail = get(options.columns.tail).map(String::from);
        let modification = get(options.columns.modification)
            .map(String::from)
            .filter(|s| !s.trim().is_empty());
        let quantity = quantity(
            get(options.columns.amount),
            get(options.columns.unit),
            &mut diagnostics,
        );
        let primer=(||->Result<Primer,Error>{
            if name.trim().is_empty()||name.len()>1024||name.chars().any(char::is_control){return Err(invalid("oligo name is empty, too long or contains controls"));}
            for column in [Some(options.columns.name),options.columns.ordered_sequence,options.columns.annealing_sequence,options.columns.tail,options.columns.modification,options.columns.amount,options.columns.unit].into_iter().flatten(){if column>=fields.len(){return Err(invalid("data row lacks a mapped column"));}}
            if modification.is_some(){return Err(unsupported("explicit modification retained; chemical semantics are not an unmodified DNA oligo"));}
            let tail=raw_tail.as_deref().unwrap_or("");let tail=if tail.trim().is_empty(){String::new()}else{Sequence::dna(tail).map_err(|_|unsupported("tail contains non-DNA/modification notation; retained literally"))?.as_str().to_string()};
            let ordered=raw_ordered_sequence.as_deref().map(Sequence::dna).transpose().map_err(|_|unsupported("ordered sequence contains non-DNA/modification notation; retained literally"))?;
            let annealing=if let Some(raw)=raw_annealing_sequence.as_deref(){Sequence::dna(raw).map_err(|_|unsupported("annealing domain contains non-DNA/modification notation; retained literally"))?.as_str().to_string()}
                else{let ordered=ordered.as_ref().ok_or_else(||invalid("missing ordered/domain sequence"))?.as_str();ordered.strip_prefix(&tail).filter(|s|!s.is_empty()).ok_or_else(||invalid("mapped tail is not a proper prefix of the ordered oligo"))?.to_string()};
            let p=Primer::new(name.clone(),&annealing,&tail)?;
            if p.ordered_sequence().len()>200{return Err(unsupported("oligo exceeds the initial 200 ordered-base primer profile"));}
            if let Some(ordered)=ordered{if ordered.as_str()!=p.ordered_sequence(){return Err(invalid("ordered sequence disagrees with separate tail/domain fields"));}}
            Ok(p)
        })().map_err(|e|diagnostics.push(e.to_string())).ok();
        if primer.is_some()
            && options.columns.annealing_sequence.is_none()
            && options.columns.tail.is_none()
        {
            diagnostics.push("No separate annealing domain/tail mapped: the full ordered sequence is used as its annealing domain.".into());
        }
        if primer
            .as_ref()
            .is_some_and(|p| p.ordered_sequence().bytes().any(|b| !b"ACGT".contains(&b)))
        {
            diagnostics.push("IUPAC ambiguity is retained as degenerate DNA, not a concrete single oligo. Modification shorthand is not inferred from valid alphabet letters.".into());
        }
        let (a, b) = parsed.ranges[row_index];
        let raw = &text.as_bytes()[a..b];
        records.push(OligoRecord {
            row_index,
            row_id: format!("{source_blake3}:{row_index}"),
            name,
            original_fields: fields,
            original_record_base64: STANDARD.encode(raw),
            original_record_blake3: blake3::hash(raw).to_hex().to_string(),
            raw_ordered_sequence,
            raw_annealing_sequence,
            raw_tail,
            modification,
            quantity,
            primer,
            diagnostics,
        });
    }
    if records.iter().any(|r| r.primer.is_none()) {
        warnings.push("Some source rows are evidence-only. They remain in the collection/original table and are not converted to DNA documents.".into());
    }
    let table = OligoTable {
        schema_version: 1,
        source_blake3,
        original_base64: STANDARD.encode(text),
        options: options.clone(),
        headers: parsed.headers,
        records,
        warnings,
    };
    bounded_json(&table)?;
    Ok(table)
}
fn serialized_hash(value: &impl Serialize) -> Result<blake3::Hash, Error> {
    struct Hash(blake3::Hasher);
    impl io::Write for Hash {
        fn write(&mut self, b: &[u8]) -> io::Result<usize> {
            self.0.update(b);
            Ok(b.len())
        }
        fn flush(&mut self) -> io::Result<()> {
            Ok(())
        }
    }
    let mut hash = Hash(blake3::Hasher::new());
    serde_json::to_writer(&mut hash, value).map_err(|e| invalid(e.to_string()))?;
    Ok(hash.0.finalize())
}
impl OligoTable {
    pub fn original_bytes(&self) -> Result<Vec<u8>, Error> {
        self.validate()?;
        decode_attachment(&self.original_base64)
    }
    pub fn validate(&self) -> Result<(), Error> {
        if self.schema_version != 1 {
            return Err(unsupported("unknown oligo collection schema"));
        }
        bounded_json(self)?;
        let original = decode_attachment(&self.original_base64)?;
        if blake3::hash(&original).to_hex().as_str() != self.source_blake3 {
            return Err(invalid("oligo table source hash mismatch"));
        }
        let source =
            std::str::from_utf8(&original).map_err(|_| invalid("oligo source is not UTF-8"))?;
        let rebuilt = import_oligos(source, &self.options)?;
        if self.headers != rebuilt.headers
            || serialized_hash(&self.records)? != serialized_hash(&rebuilt.records)?
        {
            return Err(invalid(
                "oligo records disagree with retained source/mapping",
            ));
        }
        Ok(())
    }
    pub fn selected_molecules(&self, selection: &[usize]) -> Result<ImportReport, Error> {
        self.validate()?;
        let mut seen = BTreeSet::new();
        let mut molecules = Vec::new();
        let mut total = 0usize;
        for &index in selection {
            if !seen.insert(index) {
                return Err(invalid("duplicate selected oligo row"));
            }
            let r = self
                .records
                .get(index)
                .ok_or_else(|| invalid("selected oligo row does not exist"))?;
            let p=r.primer.as_ref().ok_or_else(||unsupported(format!("source row {} has unsupported modification/sequence fields; keep it as collection evidence",index+1)))?;
            let mut m = Molecule::dna(&r.name, &p.ordered_sequence(), Topology::Linear)?;
            m.features.push(Feature {
                kind: "primer_bind".into(),
                location: Location::range(p.tail.len(), p.ordered_sequence().len()),
                qualifiers: vec![
                    ("label".into(), Some(p.name.clone())),
                    ("oligo_sequence".into(), Some(p.ordered_sequence())),
                    (
                        "annealing_sequence".into(),
                        Some(p.annealing.as_str().into()),
                    ),
                    ("5prime_tail".into(), Some(p.tail.clone())),
                ],
                needs_review: false,
            });
            m.metadata = BTreeMap::from([
                ("genbank.molecule_type".into(), "ss-DNA".into()),
                ("reference_kind".into(), "oligo".into()),
                ("oligo:table_blake3".into(), self.source_blake3.clone()),
                ("oligo:row_id".into(), r.row_id.clone()),
                (
                    "oligo:original_record_base64".into(),
                    r.original_record_base64.clone(),
                ),
                (
                    "oligo:original_record_blake3".into(),
                    r.original_record_blake3.clone(),
                ),
                (
                    "oligo:source_fields".into(),
                    serde_json::to_string(&r.original_fields)
                        .map_err(|e| invalid(e.to_string()))?,
                ),
                (
                    "oligo:source_mapping".into(),
                    serde_json::to_string(&self.options).map_err(|e| invalid(e.to_string()))?,
                ),
            ]);
            m.validate()?;
            total = total.saturating_add(
                serde_json::to_vec(&m)
                    .map_err(|e| invalid(e.to_string()))?
                    .len(),
            );
            if total > MAX_IMPORT_BYTES {
                return Err(Error::TooLarge);
            }
            molecules.push(m);
        }
        Ok(ImportReport{molecules,warnings:vec!["Imported oligos are single-stranded reference documents, not intact duplex cloning templates. Use their primer fields for PCR.".into(),"The complete original table is retained once in the oligo collection. A single-row document retains its exact row and table hash, not the entire table; save the workspace/collection to retain every source row.".into()]})
    }
}
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum TableExportMode {
    SpreadsheetSafe,
    MachineData,
}
fn formula_cell(cell: &str) -> bool {
    let significant = cell.trim_start_matches(|c: char| c.is_whitespace() || c == '\u{feff}');
    significant.starts_with(['=', '+', '-', '@']) || cell.starts_with(['\t', '\r', '\n'])
}
pub fn export_oligos(
    table: &OligoTable,
    delimiter: Delimiter,
    mode: TableExportMode,
) -> Result<ExportReport, Error> {
    table.validate()?;
    let mut writer = csv::WriterBuilder::new()
        .delimiter(delimiter.byte())
        .terminator(csv::Terminator::CRLF)
        .quote_style(csv::QuoteStyle::Always)
        .from_writer(Vec::new());
    writer
        .write_record([
            "Name",
            "Ordered sequence",
            "Annealing domain",
            "5prime tail",
            "Modifications",
            "Amount",
            "Unit",
            "Original fields JSON",
            "Source row",
        ])
        .map_err(|e| invalid(e.to_string()))?;
    let mut escaped = 0;
    for r in &table.records {
        let mut row = vec![
            r.name.clone(),
            r.raw_ordered_sequence.clone().unwrap_or_else(|| {
                r.primer
                    .as_ref()
                    .map(Primer::ordered_sequence)
                    .unwrap_or_default()
            }),
            r.primer
                .as_ref()
                .map(|p| p.annealing.as_str().into())
                .or_else(|| r.raw_annealing_sequence.clone())
                .unwrap_or_default(),
            r.primer
                .as_ref()
                .map(|p| p.tail.clone())
                .or_else(|| r.raw_tail.clone())
                .unwrap_or_default(),
            r.modification.clone().unwrap_or_default(),
            r.quantity
                .as_ref()
                .map(|q| q.raw_value.clone())
                .unwrap_or_default(),
            r.quantity
                .as_ref()
                .map(|q| q.raw_unit.clone())
                .unwrap_or_default(),
            serde_json::to_string(&r.original_fields).map_err(|e| invalid(e.to_string()))?,
            (r.row_index + 1).to_string(),
        ];
        for cell in &mut row {
            if formula_cell(cell) {
                escaped += 1;
                if matches!(mode, TableExportMode::SpreadsheetSafe) {
                    cell.insert(0, '\'');
                }
            }
        }
        writer
            .write_record(row)
            .map_err(|e| invalid(e.to_string()))?;
        if writer.get_ref().len() > MAX_IMPORT_BYTES {
            return Err(Error::TooLarge);
        }
    }
    let bytes = writer.into_inner().map_err(|e| invalid(e.to_string()))?;
    if bytes.len() > MAX_IMPORT_BYTES {
        return Err(Error::TooLarge);
    }
    let mut warnings=vec!["Order-sheet columns use explicit mapped roles. Unmapped source columns remain in Original fields JSON and the retained original table.".into()];
    if escaped > 0 {
        warnings.push(match mode{TableExportMode::SpreadsheetSafe=>format!("Prefixed {escaped} formula-like cells with an apostrophe for spreadsheet safety. Those display-cell values intentionally differ; original fields/source remain retained."),TableExportMode::MachineData=>format!("Machine-data export retains {escaped} formula-like cells literally. Treat this file as data; quoting alone does not disable spreadsheet formulas.")});
    }
    Ok(ExportReport {
        text: String::from_utf8(bytes).map_err(|e| invalid(e.to_string()))?,
        warnings,
    })
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum QualityEncoding {
    Phred33,
    Phred64,
}
impl QualityEncoding {
    fn offset(self) -> u8 {
        match self {
            Self::Phred33 => 33,
            Self::Phred64 => 64,
        }
    }
    fn text(self) -> &'static str {
        match self {
            Self::Phred33 => "phred33",
            Self::Phred64 => "phred64",
        }
    }
}
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FastqEvidence {
    pub header: String,
    pub encoding: QualityEncoding,
    pub called_sequence: Sequence,
    pub qualities: Vec<u8>,
    pub source_call_length: usize,
    /// Half-open coordinates in the original full forward FASTQ record, regardless of orientation.
    pub retained_start: usize,
    pub retained_end: usize,
    pub reverse: bool,
    pub original_record_blake3: String,
    pub current_valid: bool,
    pub warnings: Vec<String>,
}
struct Read {
    header: String,
    id: String,
    sequence: Sequence,
    qualities: Vec<u8>,
    start: usize,
    end: usize,
}
fn next_line<'a>(lines: &mut impl Iterator<Item = &'a str>, offset: &mut usize) -> Option<&'a str> {
    let raw = lines.next()?;
    *offset += raw.len();
    let line = if let Some(line) = raw.strip_suffix('\n') {
        line.strip_suffix('\r').unwrap_or(line)
    } else {
        raw
    };
    Some(line)
}
fn parse_fastq(bytes: &[u8], encoding: QualityEncoding) -> Result<Vec<Read>, Error> {
    if bytes.is_empty() || bytes.len() > MAX_IMPORT_BYTES {
        return Err(if bytes.is_empty() {
            invalid("FASTQ has no records")
        } else {
            Error::TooLarge
        });
    }
    let text = std::str::from_utf8(bytes)
        .map_err(|_| invalid("FASTQ must be UTF-8 header text with ASCII DNA/quality"))?;
    let mut lines = text.split_inclusive('\n');
    let mut offset = 0;
    let mut reads = Vec::new();
    let mut total = 0usize;
    while offset < bytes.len() {
        let start = offset;
        let header =
            next_line(&mut lines, &mut offset).ok_or_else(|| invalid("missing FASTQ header"))?;
        let header = header.strip_prefix('@').ok_or_else(|| {
            invalid("FASTQ record must begin with @; blank/inter-record lines are not a record")
        })?;
        if header.len() > 8192
            || header.is_empty()
            || header.chars().any(|c| c.is_control() && c != '\t')
        {
            return Err(invalid("invalid FASTQ identifier/description"));
        }
        let id = header
            .split_ascii_whitespace()
            .next()
            .filter(|id| !id.is_empty())
            .ok_or_else(|| invalid("FASTQ identifier is empty"))?
            .to_string();
        if header.starts_with(char::is_whitespace) {
            return Err(invalid("FASTQ identifier cannot begin with whitespace"));
        }
        let mut bases = String::new();
        loop {
            let line = next_line(&mut lines, &mut offset)
                .ok_or_else(|| invalid("FASTQ lacks + delimiter"))?;
            if let Some(plus) = line.strip_prefix('+') {
                if !plus.is_empty() && plus != header && plus != id {
                    return Err(invalid("FASTQ + identifier differs from header"));
                }
                break;
            }
            if line.is_empty()
                || line
                    .bytes()
                    .any(|b| !b"ACGTRYSWKMBDHVNacgtryswkmbdhvn".contains(&b))
            {
                return Err(invalid(
                    "FASTQ sequence contains whitespace, RNA, gaps or unsupported bases",
                ));
            }
            if bases.len().saturating_add(line.len()) > collagene_model::MAX_SEQUENCE_LEN {
                return Err(Error::TooLarge);
            }
            bases.push_str(line);
        }
        let sequence = Sequence::dna(&bases)?;
        let mut qualities = Vec::with_capacity(sequence.len());
        while qualities.len() < sequence.len() {
            let line = next_line(&mut lines, &mut offset)
                .ok_or_else(|| invalid("FASTQ quality string is truncated"))?;
            if line.is_empty() || line.len() > sequence.len() - qualities.len() {
                return Err(invalid("FASTQ quality length differs from sequence"));
            }
            for byte in line.bytes() {
                if !(encoding.offset()..=126).contains(&byte) {
                    return Err(invalid(
                        "FASTQ score is outside the explicitly selected ASCII Phred encoding",
                    ));
                }
                qualities.push(byte - encoding.offset());
            }
        }
        if reads.len() >= MAX_FASTQ_RECORDS {
            return Err(Error::TooLarge);
        }
        total = total.checked_add(sequence.len()).ok_or(Error::TooLarge)?;
        if total > collagene_model::MAX_SEQUENCE_LEN {
            return Err(Error::TooLarge);
        }
        reads.push(Read {
            header: header.into(),
            id,
            sequence,
            qualities,
            start,
            end: offset,
        });
    }
    Ok(reads)
}
pub fn import_fastq(bytes: &[u8], encoding: QualityEncoding) -> Result<ImportReport, Error> {
    let reads = parse_fastq(bytes, encoding)?;
    let source_hash = blake3::hash(bytes).to_hex().to_string();
    let mut molecules = Vec::new();
    let mut size = 0usize;
    for (index, r) in reads.into_iter().enumerate() {
        let original = &bytes[r.start..r.end];
        let estimate = r
            .sequence
            .len()
            .saturating_add(r.qualities.len().div_ceil(3) * 4)
            .saturating_add(original.len().div_ceil(3) * 4);
        if estimate > MAX_IMPORT_BYTES {
            return Err(Error::TooLarge);
        }
        let mut molecule = Molecule::dna(r.id, r.sequence.as_str(), Topology::Linear)?;
        molecule.metadata=BTreeMap::from([(FASTQ_ORIGINAL_KEY.into(),STANDARD.encode(original)),("fastq:original_blake3".into(),blake3::hash(original).to_hex().to_string()),("fastq:encoding".into(),encoding.text().into()),("fastq:header".into(),r.header),("fastq:quality_base64".into(),STANDARD.encode(&r.qualities)),("fastq:called_sequence_blake3".into(),blake3::hash(r.sequence.as_str().as_bytes()).to_hex().to_string()),("fastq:source_blake3".into(),source_hash.clone()),("fastq:record_index".into(),index.to_string()),("fastq:orientation".into(),"forward".into()),("fastq:quality_model".into(),"Phred values under an explicitly selected ASCII offset; caller calibration is not independently confirmed.".into())]);
        molecule.validate()?;
        size = size.saturating_add(
            serde_json::to_vec(&molecule)
                .map_err(|e| invalid(e.to_string()))?
                .len(),
        );
        if size > MAX_IMPORT_BYTES {
            return Err(Error::TooLarge);
        }
        molecules.push(molecule);
    }
    Ok(ImportReport{molecules,warnings:vec![format!("FASTQ imported with explicit {} encoding; no offset was inferred. Original record bytes and supplied scores are retained; paired-read relationships/platform calibration are not inferred.",encoding.text()),"Edits do not invent or silently realign quality values. Current FASTQ export requires the attached called sequence/orientation to remain valid.".into()]})
}
pub fn original_fastq_bytes(molecule: &Molecule) -> Result<Vec<u8>, Error> {
    let raw = molecule
        .metadata
        .get(FASTQ_ORIGINAL_KEY)
        .ok_or_else(|| unsupported("document has no FASTQ source attachment"))?;
    let bytes = decode_attachment(raw)?;
    if molecule
        .metadata
        .get("fastq:original_blake3")
        .map(String::as_str)
        != Some(blake3::hash(&bytes).to_hex().as_str())
    {
        return Err(invalid("FASTQ original attachment hash mismatch"));
    }
    Ok(bytes)
}
pub fn quality_evidence(molecule: &Molecule) -> Result<FastqEvidence, Error> {
    molecule.validate()?;
    let original = original_fastq_bytes(molecule)?;
    let encoding = match molecule.metadata.get("fastq:encoding").map(String::as_str) {
        Some("phred33") => QualityEncoding::Phred33,
        Some("phred64") => QualityEncoding::Phred64,
        _ => return Err(invalid("unknown FASTQ source encoding")),
    };
    let mut records = parse_fastq(&original, encoding)?;
    if records.len() != 1 {
        return Err(invalid(
            "FASTQ attachment must contain exactly its own source record",
        ));
    }
    let r = records.remove(0);
    if molecule.metadata.get("fastq:header") != Some(&r.header)
        || molecule
            .metadata
            .get("fastq:called_sequence_blake3")
            .map(String::as_str)
            != Some(
                blake3::hash(r.sequence.as_str().as_bytes())
                    .to_hex()
                    .as_str(),
            )
        || decode_attachment(
            molecule
                .metadata
                .get("fastq:quality_base64")
                .ok_or_else(|| invalid("FASTQ qualities missing"))?,
        )? != r.qualities
    {
        return Err(invalid(
            "FASTQ quality/call/header metadata disagrees with retained source",
        ));
    }
    let source_call_length = r.sequence.len();
    let reverse = match molecule
        .metadata
        .get("fastq:orientation")
        .map(String::as_str)
    {
        Some("forward") => false,
        Some("reverse") => true,
        _ => return Err(invalid("invalid FASTQ source orientation")),
    };
    let retained_start = molecule
        .metadata
        .get("fastq:retained_start")
        .map(|s| s.parse::<usize>())
        .transpose()
        .map_err(|_| invalid("invalid FASTQ retained start"))?
        .unwrap_or(0);
    let retained_end = molecule
        .metadata
        .get("fastq:retained_end")
        .map(|s| s.parse::<usize>())
        .transpose()
        .map_err(|_| invalid("invalid FASTQ retained end"))?
        .unwrap_or(source_call_length);
    if retained_start >= retained_end || retained_end > source_call_length {
        return Err(invalid("FASTQ retained calls are out of bounds"));
    }
    let mut called_sequence = Sequence::dna(&r.sequence.as_str()[retained_start..retained_end])?;
    let mut qualities = r.qualities[retained_start..retained_end].to_vec();
    if reverse {
        called_sequence = called_sequence.reverse_complement()?;
        qualities.reverse();
    }
    let mapped = match molecule.metadata.get("annotation_mapping:operation") {
        None => true,
        Some(op) => {
            (op == "reverse complement" || op.starts_with("source interval "))
                && molecule.metadata.get("fastq:mapped_operation") == Some(op)
                && molecule.metadata.get("fastq:mapped_operation_source")
                    == molecule.metadata.get("annotation_mapping:source_revision")
        }
    };
    let current_valid =
        molecule.sequence == called_sequence && molecule.topology == Topology::Linear && mapped;
    Ok(FastqEvidence {
        header: r.header,
        encoding,
        called_sequence,
        qualities,
        source_call_length,
        retained_start,
        retained_end,
        reverse,
        original_record_blake3: blake3::hash(&original).to_hex().to_string(),
        current_valid,
        warnings: if current_valid {
            vec!["Quality values are source caller evidence under the selected encoding, not independently calibrated probabilities.".into()]
        } else {
            vec!["Current bases/topology/coordinate transform no longer establish a valid one-to-one FASTQ quality mapping. Original source remains downloadable.".into()]
        },
    })
}
/// Call after the generic reverse-complement transformation has recorded its operation/source.
pub fn project_fastq_reverse_complement(
    source: &Molecule,
    result: &mut Molecule,
) -> Result<(), Error> {
    if !source.metadata.contains_key(FASTQ_ORIGINAL_KEY) {
        return Ok(());
    }
    let evidence = quality_evidence(source)?;
    if result.sequence != source.sequence.reverse_complement()? {
        return Err(invalid(
            "FASTQ reverse projection does not match reference transformation",
        ));
    }
    if evidence.current_valid {
        let orientation =
            if source.metadata.get("fastq:orientation").map(String::as_str) == Some("reverse") {
                "forward"
            } else {
                "reverse"
            };
        result
            .metadata
            .insert("fastq:orientation".into(), orientation.into());
        result.metadata.insert(
            "fastq:retained_start".into(),
            evidence.retained_start.to_string(),
        );
        result.metadata.insert(
            "fastq:retained_end".into(),
            evidence.retained_end.to_string(),
        );
        result
            .metadata
            .insert("fastq:mapped_operation_source".into(), source.revision());
        result
            .metadata
            .insert("fastq:mapped_operation".into(), "reverse complement".into());
    } else {
        result.metadata.remove("fastq:mapped_operation_source");
        result.metadata.remove("fastq:mapped_operation");
    }
    Ok(())
}
/// Trim a valid called read by explicit half-open current call coordinates, keeping
/// the original attachment and exactly the corresponding caller quality slice.
pub fn trim_fastq(source: &Molecule, start: usize, end: usize) -> Result<Molecule, Error> {
    let evidence = quality_evidence(source)?;
    if !evidence.current_valid {
        return Err(unsupported("cannot trim stale FASTQ quality coordinates"));
    }
    if start >= end || end > evidence.called_sequence.len() {
        return Err(Error::Coordinates(
            "FASTQ trim is empty or outside current calls".into(),
        ));
    }
    let mut result = crate::annotations::extract_region(source, start, end)?;
    let (retained_start, retained_end) = if evidence.reverse {
        (evidence.retained_end - end, evidence.retained_end - start)
    } else {
        (
            evidence.retained_start + start,
            evidence.retained_start + end,
        )
    };
    result
        .metadata
        .insert("fastq:retained_start".into(), retained_start.to_string());
    result
        .metadata
        .insert("fastq:retained_end".into(), retained_end.to_string());
    result
        .metadata
        .insert("fastq:mapped_operation_source".into(), source.revision());
    let operation = result
        .metadata
        .get("annotation_mapping:operation")
        .cloned()
        .ok_or_else(|| invalid("trim lost its coordinate provenance"))?;
    result
        .metadata
        .insert("fastq:mapped_operation".into(), operation);
    result.validate()?;
    crate::protocol::save_document(&result)?;
    Ok(result)
}
pub fn export_fastq(molecule: &Molecule, encoding: QualityEncoding) -> Result<ExportReport, Error> {
    let evidence = quality_evidence(molecule)?;
    if !evidence.current_valid {
        return Err(unsupported(
            "current FASTQ quality mapping is stale; export original evidence or sequence-only FASTA",
        ));
    }
    if molecule.name.trim().is_empty()
        || molecule.name.starts_with(char::is_whitespace)
        || molecule.name.len() > 8192
        || molecule.name.chars().any(char::is_control)
    {
        return Err(invalid(
            "FASTQ output header is empty, too long or contains controls",
        ));
    }
    let original_id = evidence
        .header
        .split_ascii_whitespace()
        .next()
        .unwrap_or("");
    let header = if molecule.name == original_id {
        evidence.header.clone()
    } else {
        molecule.name.clone()
    };
    let mut quality = Vec::with_capacity(evidence.qualities.len());
    for q in evidence.qualities {
        quality.push(q.checked_add(encoding.offset()).filter(|q|*q<=126).ok_or_else(||unsupported("source quality cannot be represented in requested FASTQ encoding without loss"))?);
    }
    let text = format!(
        "@{header}\n{}\n+\n{}\n",
        molecule.sequence.as_str(),
        String::from_utf8(quality).map_err(|e| invalid(e.to_string()))?
    );
    if text.len() > MAX_IMPORT_BYTES {
        return Err(Error::TooLarge);
    }
    Ok(ExportReport{text,warnings:vec!["Current FASTQ export uses normalized uppercase DNA, canonical four-line framing and supplied per-base qualities. Original wrapping/repeated + identifiers remain in original-source download; rich annotations are not represented.".into()]})
}
