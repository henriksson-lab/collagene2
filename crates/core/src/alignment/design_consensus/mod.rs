//! Supplied-consensus sequence/CDS observations conditional on one actual path.
//! This owner does not replay cloning history or infer read quality/material state.
mod assessment;
mod coding;
mod contract;
mod model;
mod projection;
mod provenance;
use crate::documents::workspace::Workspace;
use collagene_model::Error;
pub use model::*;
pub const ENGINE_VERSION: &str = "design-consensus-1";
pub const MAX_REQUEST_BYTES: usize = 128 * 1024;
pub const MAX_REPORT_BYTES: usize = 32 * 1024 * 1024;
pub const MAX_CDS: usize = 32;
pub const STATEMENT: &str = "Supplied consensus and independently selected design are immutable sequence observations. Differences and coding consequences are conditional on the actual selected affine alignment path and declared annotation/code/coordinate assumptions; a certified score does not establish unique correspondence. Provenance availability is separate from sequence agreement. No read depth, Phred confidence, native edit cause, physical material, expression/function or experimental confirmation is inferred.";
pub fn revision(workspace: &Workspace, request: &Request) -> Result<String, Error> {
    contract::revision(workspace, request)
}
pub fn assess(workspace: &Workspace, request: &Request) -> Result<Report, Error> {
    assessment::assess(workspace, request)
}

pub(super) struct Work {
    pub limits: Limits,
    pub summary: WorkSummary,
}
impl Work {
    pub fn new(limits: Limits) -> Self {
        Self {
            limits,
            summary: WorkSummary {
                mapping_work: 0,
                coding_calls: 0,
                segments: 0,
                consequences: 0,
            },
        }
    }
    pub fn mapping(&mut self, count: usize) -> bool {
        match self
            .summary
            .mapping_work
            .checked_add(count)
            .filter(|value| *value <= self.limits.max_mapping_work)
        {
            Some(value) => {
                self.summary.mapping_work = value;
                true
            }
            None => {
                self.summary.mapping_work = self.limits.max_mapping_work;
                false
            }
        }
    }
    pub fn coding(&mut self, count: usize) -> bool {
        match self
            .summary
            .coding_calls
            .checked_add(count)
            .filter(|value| *value <= self.limits.max_coding_calls)
        {
            Some(value) => {
                self.summary.coding_calls = value;
                true
            }
            None => {
                self.summary.coding_calls = self.limits.max_coding_calls;
                false
            }
        }
    }
}
