//! Read-only identity and display fields for an exact retained source event.
use collagene_core::{Alphabet, Topology, run_lineage::EventBinding};
use serde::{Deserialize, Serialize};
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct SourceInfo {
    pub binding: EventBinding,
    pub name: String,
    pub length: usize,
    pub alphabet: Alphabet,
    pub topology: Topology,
}
