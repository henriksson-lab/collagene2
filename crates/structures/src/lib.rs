//! Imported deposited coordinates, never a fold prediction or cloning verdict.
//! Independent implementation of the documented wwPDB PDB/mmCIF profiles.
mod cif;
mod pdb;
pub mod scene;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};

pub const ENGINE: &str = "structure-import-1";
pub const MAX_SOURCE_BYTES: usize = 16 << 20;
pub const MAX_ATOMS: usize = 100_000;
pub const MAX_CONNECTIONS: usize = 400_000;
pub const MAX_MODELS: usize = 1024;
pub(crate) const MAX_FIELD: usize = 256;

#[derive(Debug, thiserror::Error, PartialEq, Eq)]
pub enum Error {
    #[error("structure resource bound exceeded: {0}")]
    Limit(&'static str),
    #[error("invalid structure at line {line}: {message}")]
    Invalid { line: usize, message: &'static str },
    #[error("unsupported structure profile: {0}")]
    Unsupported(&'static str),
    #[error("structure does not match its original source or current parser")]
    Binding,
}
pub type Result<T> = std::result::Result<T, Error>;
pub(crate) fn bad(line: usize, message: &'static str) -> Error {
    Error::Invalid { line, message }
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum Format {
    Pdb,
    Mmcif,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum Record {
    Atom,
    Hetatm,
}

/// A missing author identifier is not replaced with a label identifier.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Identity {
    pub atom: Option<String>,
    pub component: Option<String>,
    pub chain: Option<String>,
    pub sequence: Option<String>,
}
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Atom {
    pub id: String,
    pub model: u32,
    pub record: Record,
    pub label: Identity,
    pub author: Identity,
    pub entity: Option<String>,
    pub alternate: Option<String>,
    pub insertion: Option<String>,
    /// Angstrom, original f64 values. Rendering transforms do not change these.
    pub position: [f64; 3],
    pub element: Option<String>,
    pub occupancy: Option<f64>,
    pub b_factor: Option<f64>,
    pub formal_charge: Option<i8>,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub enum CifValue {
    Text(String),
    Unknown,
    NotApplicable,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub enum Connection {
    /// PDB CONECT names serials, not vector offsets. None applies to all models.
    Pdb {
        model: Option<u32>,
        from: String,
        to: String,
    },
    /// Entire deposited _struct_conn row. Symmetry/alternate selectors remain
    /// explicit; this import never guesses a bond from atom distance or names.
    Mmcif { fields: BTreeMap<String, CifValue> },
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ImportNotes {
    pub ignored_pdb_records: BTreeMap<String, usize>,
    pub unprojected_cif_categories: BTreeSet<String>,
    pub atoms_without_element: usize,
    pub atoms_with_alternates: usize,
    pub connection_declarations: usize,
}
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Structure {
    engine: String,
    format: Format,
    source_hash: String,
    original: Vec<u8>,
    block: Option<String>,
    atoms: Vec<Atom>,
    connections: Vec<Connection>,
    notes: ImportNotes,
}
impl Structure {
    pub fn atoms(&self) -> &[Atom] {
        &self.atoms
    }
    pub fn connections(&self) -> &[Connection] {
        &self.connections
    }
    pub fn notes(&self) -> &ImportNotes {
        &self.notes
    }
    pub fn original(&self) -> &[u8] {
        &self.original
    }
    pub fn source_hash(&self) -> &str {
        &self.source_hash
    }
    pub fn format(&self) -> Format {
        self.format
    }
    pub fn block(&self) -> Option<&str> {
        self.block.as_deref()
    }
    pub fn models(&self) -> BTreeSet<u32> {
        self.atoms.iter().map(|a| a.model).collect()
    }
    /// Imported/transported data grants no scientific authority. Reparse retained
    /// bytes before treating a deserialized projection as current import evidence.
    pub fn validate(&self) -> Result<()> {
        if self.engine != ENGINE
            || self.atoms.len() > MAX_ATOMS
            || self.connections.len() > MAX_CONNECTIONS
        {
            return Err(Error::Binding);
        }
        if import(&self.original, self.format)? != *self {
            return Err(Error::Binding);
        }
        Ok(())
    }
}

pub fn import(bytes: &[u8], format: Format) -> Result<Structure> {
    if bytes.is_empty() || bytes.len() > MAX_SOURCE_BYTES {
        return Err(Error::Limit("source bytes"));
    }
    if bytes
        .iter()
        .any(|&b| !matches!(b, b'\t' | b'\n' | b'\r' | 32..=126))
    {
        return Err(Error::Unsupported(
            "ASCII coordinate file required; CIF2 is not supported",
        ));
    }
    let text =
        std::str::from_utf8(bytes).map_err(|_| Error::Unsupported("ASCII source required"))?;
    let mut out = Structure {
        engine: ENGINE.into(),
        format,
        source_hash: blake3::hash(bytes).to_hex().to_string(),
        original: bytes.to_vec(),
        block: None,
        atoms: vec![],
        connections: vec![],
        notes: ImportNotes {
            ignored_pdb_records: BTreeMap::new(),
            unprojected_cif_categories: BTreeSet::new(),
            atoms_without_element: 0,
            atoms_with_alternates: 0,
            connection_declarations: 0,
        },
    };
    match format {
        Format::Pdb => pdb::parse(text, &mut out)?,
        Format::Mmcif => cif::parse(text, &mut out)?,
    }
    if out.atoms.is_empty() {
        return Err(bad(1, "no coordinate atoms"));
    }
    let mut ids = BTreeSet::new();
    let mut models = BTreeSet::new();
    for a in &out.atoms {
        if !ids.insert((a.model, &a.id)) {
            return Err(bad(1, "duplicate atom ID within a model"));
        }
        models.insert(a.model);
        out.notes.atoms_without_element += usize::from(a.element.is_none());
        out.notes.atoms_with_alternates += usize::from(a.alternate.is_some());
    }
    if models.len() > MAX_MODELS {
        return Err(Error::Limit("models"));
    }
    out.notes.connection_declarations = out.connections.len();
    Ok(out)
}
pub(crate) fn push_atom(out: &mut Structure, atom: Atom) -> Result<()> {
    if out.atoms.len() >= MAX_ATOMS {
        return Err(Error::Limit("atoms"));
    }
    out.atoms.push(atom);
    Ok(())
}
pub(crate) fn push_connection(out: &mut Structure, value: Connection) -> Result<()> {
    if out.connections.len() >= MAX_CONNECTIONS {
        return Err(Error::Limit("connections"));
    }
    out.connections.push(value);
    Ok(())
}
pub(crate) fn field(s: &str, line: usize) -> Result<String> {
    if s.len() > MAX_FIELD || s.bytes().any(|b| b.is_ascii_control()) {
        return Err(bad(line, "invalid identifier field"));
    }
    Ok(s.into())
}
pub(crate) fn number(s: &str, line: usize) -> Result<f64> {
    let n: f64 = s.parse().map_err(|_| bad(line, "invalid numeric field"))?;
    if !n.is_finite() {
        return Err(bad(line, "nonfinite numeric field"));
    }
    Ok(n)
}
pub(crate) fn positive(s: &str, line: usize) -> Result<u32> {
    let n = s
        .parse::<u32>()
        .map_err(|_| bad(line, "invalid positive integer"))?;
    if n == 0 {
        return Err(bad(line, "zero model or serial number"));
    }
    Ok(n)
}
pub(crate) fn element(s: &str, line: usize) -> Result<Option<String>> {
    if s.is_empty() {
        return Ok(None);
    }
    // Element comes only from its declared field, never the atom label (CA).
    const SYMBOLS: &str = "H HE LI BE B C N O F NE NA MG AL SI P S CL AR K CA SC TI V CR MN FE CO NI CU ZN GA GE AS SE BR KR RB SR Y ZR NB MO TC RU RH PD AG CD IN SN SB TE I XE CS BA LA CE PR ND PM SM EU GD TB DY HO ER TM YB LU HF TA W RE OS IR PT AU HG TL PB BI PO AT RN FR RA AC TH PA U NP PU AM CM BK CF ES FM MD NO LR RF DB SG BH HS MT DS RG CN NH FL MC LV TS OG D";
    let name = s.to_ascii_uppercase();
    if !SYMBOLS.split(' ').any(|e| e == name) {
        return Err(bad(line, "unknown element symbol"));
    }
    Ok(Some(name))
}
