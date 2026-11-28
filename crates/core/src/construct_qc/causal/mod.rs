//! Causal structural CDS comparisons over freshly replayed native history.
//! Declared feature pairing never substitutes for DNA/annotation correspondence.
use crate::documents::workspace::Workspace;
use collagene_model::Error;
mod coding;
mod comparison;
mod contract;
mod history;
mod register;
mod stops;
#[cfg(test)]
mod tests;
mod types;
pub use types::*;
pub const ENGINE_VERSION: &str = "construct-causal-qc-1";
pub const MAX_REQUEST_BYTES: usize = 128 * 1024;
pub const MAX_REPORT_BYTES: usize = 2 * 1024 * 1024;
pub const STATEMENT: &str = "Native-history construct QC compares explicitly declared CDS roles and genetic code against actual retained reference ancestry. Annotation declarations do not verify biological identity. Structural findings do not establish expression, function or experimental success; missing recipes or incomplete correspondence remain indeterminate.";

pub fn revision(workspace: &Workspace, request: &Request) -> Result<String, Error> {
    contract::revision(workspace, request)
}
pub fn assess(workspace: &Workspace, request: &Request) -> Result<Report, Error> {
    comparison::assess(workspace, request)
}

/// Actual receipt traversal/comparison work for one complete assessment.
/// Native scientific engines and source-copy admission retain their own explicit bounds.
pub(super) struct Work {
    maximum: usize,
    used: usize,
}
impl Work {
    pub fn new(maximum: usize) -> Self {
        Self { maximum, used: 0 }
    }
    pub fn charge(&mut self, amount: usize) -> bool {
        match self
            .used
            .checked_add(amount)
            .filter(|used| *used <= self.maximum)
        {
            Some(used) => {
                self.used = used;
                true
            }
            None => {
                self.used = self.maximum;
                false
            }
        }
    }
    pub fn remaining(&self) -> usize {
        self.maximum - self.used
    }
    pub fn used(&self) -> usize {
        self.used
    }
}
