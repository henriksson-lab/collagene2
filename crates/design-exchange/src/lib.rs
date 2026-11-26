//! Explicit structural SBOL3.1 JSON-LD profile; no remote context fetches.
mod jsonld;
mod project;
mod write;
pub use collagene_model::{Alphabet, Feature, Location, Molecule, Sequence, Topology};
pub use project::{Projection, ProjectionPolicy, Selection};
use serde::{Deserialize, Serialize};
use serde_json::Value;
use std::collections::BTreeMap;
pub use write::{Export, ExportIds};
pub const PROFILE: &str = "collagene-sbol-3.1-structural-jsonld-1";
pub const SBOL: &str = "http://sbols.org/v3#";
pub const DNA: &str = "https://identifiers.org/SBO:0000251";
pub const RNA: &str = "https://identifiers.org/SBO:0000250";
pub const PROTEIN: &str = "https://identifiers.org/SBO:0000252";
pub const LINEAR: &str = "http://identifiers.org/SO:0000987";
pub const CIRCULAR: &str = "http://identifiers.org/SO:0000988";
pub const INLINE: &str = "https://identifiers.org/SO:0001030";
pub const REVERSE: &str = "https://identifiers.org/SO:0001031";
// Four explicit identifier aliases only. Never rewrite retained graph IRIs,
// follow redirects, or generalize HTTP/HTTPS equivalence to arbitrary terms.
pub(crate) fn semantic_so_type(iri: &str) -> Option<&'static str> {
    match iri {
        "http://identifiers.org/SO:0000987" | "https://identifiers.org/SO:0000987" => Some(LINEAR),
        "http://identifiers.org/SO:0000988" | "https://identifiers.org/SO:0000988" => {
            Some(CIRCULAR)
        }
        "http://identifiers.org/SO:0000984" | "https://identifiers.org/SO:0000984" => {
            Some("http://identifiers.org/SO:0000984")
        }
        "http://identifiers.org/SO:0000985" | "https://identifiers.org/SO:0000985" => {
            Some("http://identifiers.org/SO:0000985")
        }
        _ => None,
    }
}

pub const NUCLEIC_ENCODING: &str = "https://identifiers.org/edam:format_1207";
pub const PROTEIN_ENCODING: &str = "https://identifiers.org/edam:format_1208";
pub const MAX_BYTES: usize = 16 << 20;
pub const MAX_OBJECTS: usize = 10_000;
pub const MAX_VALUES: usize = 200_000;
pub const MAX_SEQUENCE_BASES: usize = 8 << 20;
#[derive(Debug, thiserror::Error)]
pub enum Error {
    #[error("SBOL profile limit exceeded")]
    Limit,
    #[error("invalid JSON-LD profile: {0}")]
    Json(&'static str),
    #[error("invalid SBOL structural profile: {0}")]
    Invalid(String),
    #[error("unsupported or lossy selected conversion: {0}")]
    Projection(String),
    #[error("native molecule: {0}")]
    Native(#[from] collagene_model::Error),
}
pub type Result<T> = std::result::Result<T, Error>;
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Issue {
    pub subject: String,
    pub code: String,
    pub message: String,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Identified {
    pub iri: String,
    pub name: Option<String>,
    pub display_id: Option<String>,
    pub description: Option<String>,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Component {
    pub types: Vec<String>,
    pub roles: Vec<String>,
    pub sequences: Vec<String>,
    pub features: Vec<String>,
    pub constraints: Vec<String>,
    pub namespace: String,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct SequenceRecord {
    pub elements: Option<String>,
    pub encoding: Option<String>,
    pub namespace: String,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub enum FeatureKind {
    SequenceFeature,
    SubComponent {
        instance_of: String,
        source_locations: Vec<String>,
        role_integration: Option<String>,
    },
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct FeatureRecord {
    pub kind: FeatureKind,
    pub roles: Vec<String>,
    pub orientation: Option<String>,
    pub locations: Vec<String>,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub enum LocationKind {
    Range { start: usize, end: usize },
    Cut { at: usize },
    EntireSequence,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct LocationRecord {
    pub kind: LocationKind,
    pub sequence: String,
    pub orientation: Option<String>,
    pub order: Option<i64>,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Constraint {
    pub subject: String,
    pub object: String,
    pub restriction: String,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub enum Kind {
    Component(Component),
    Sequence(SequenceRecord),
    Feature(FeatureRecord),
    Location(LocationRecord),
    Constraint(Constraint),
    Collection {
        namespace: String,
        members: Vec<String>,
    },
    Unsupported,
}
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Object {
    pub identity: Identified,
    pub classes: Vec<String>,
    pub kind: Kind,
    pub retained_properties: BTreeMap<String, Value>,
}
/// Immutable imported graph plus original evidence. No field is a physical
/// material receipt, executable URL, current verification or native event ID.
#[derive(Debug, Clone)]
pub struct Document {
    original: Vec<u8>,
    objects: BTreeMap<String, Object>,
    expanded: Vec<Value>,
    issues: Vec<Issue>,
}
impl Document {
    pub fn objects(&self) -> &BTreeMap<String, Object> {
        &self.objects
    }
    pub fn issues(&self) -> &[Issue] {
        &self.issues
    }
    pub fn original_bytes(&self) -> &[u8] {
        &self.original
    }
    pub fn revision(&self) -> String {
        blake3::hash(&self.original).to_hex().to_string()
    }
    pub fn expanded_bytes(&self) -> Result<Vec<u8>> {
        write::encode(&serde_json::json!({"@graph":self.expanded}))
    }
    pub fn project(&self, selection: &Selection, policy: ProjectionPolicy) -> Result<Projection> {
        project::project(self, selection, policy)
    }
}
pub fn import(bytes: &[u8]) -> Result<Document> {
    jsonld::import(bytes)
}
pub fn export_molecule(
    molecule: &Molecule,
    ids: &ExportIds,
    policy: ProjectionPolicy,
) -> Result<Export> {
    write::export(molecule, ids, policy)
}
fn invalid(message: impl Into<String>) -> Error {
    Error::Invalid(message.into())
}
fn issue(subject: &str, code: &str, message: &str) -> Issue {
    Issue {
        subject: subject.into(),
        code: code.into(),
        message: message.into(),
    }
}
