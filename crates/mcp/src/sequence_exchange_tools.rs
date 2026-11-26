//! Typed compatibility profiles, without host paths or restored chemistry authority.
use crate::DocumentRef;
use collagene_core::{collagene_xml_compat, gcg_sequence};
use schemars::JsonSchema;
use serde::Deserialize;
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub enum ImportProfile {
    Gcg(gcg_sequence::ImportOptions),
    Fasta {
        /// None retains per-record declarations and defaults unmarked records to DNA.
        /// Some explicitly sets all records, refusing contradictory declarations.
        alphabet: Option<collagene_core::Alphabet>,
    },
    CollageneXml {
        policy: collagene_xml_compat::ReferencePolicy,
    },
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ImportQuery {
    pub text: String,
    pub profile: ImportProfile,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct SnapGeneExportQuery {
    pub document: DocumentRef,
    /// Default false returns a review instead of dropping unrepresentable fields.
    #[serde(default)]
    pub accept_omissions: bool,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct GcgExportQuery {
    pub document: DocumentRef,
    pub options: gcg_sequence::ExportOptions,
    /// Must explicitly acknowledge annotations/topology/history omitted by GCG.
    #[serde(default)]
    pub accept_omissions: bool,
}
