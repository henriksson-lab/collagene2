pub mod catalog_tools;

use crate::server::CollageneServer;

use crate::server::response::output;
use rmcp::ErrorData;
use rmcp::handler::server::wrapper::Parameters;
use rmcp::model::CallToolResult;
use rmcp::tool;
use rmcp::tool_router;

#[tool_router(router = catalog_router, vis = "pub(crate)")]
impl CollageneServer {
    #[tool(
        description = "Search or exactly look up the embedded versioned restriction catalog. Returns edition/hash, stable enzyme geometry identities, separate formulation records, aliases and documentary source citations. Search/exact ambiguity pages are bounded to 100 records. Optional edition/hash pin rejects stale editions. Enzyme identity and formulation assay facts are distinct; lookup does not predict cleavage or transfer supplier conditions between related enzymes."
    )]
    pub(crate) fn lookup_enzyme_catalog(
        &self,
        Parameters(args): Parameters<catalog_tools::LookupQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        output(catalog_tools::lookup(args)?)
    }
    #[tool(
        description = "Resolve one explicit formulation/condition/assay context from the pinned embedded documentary catalog. Keeps Known, Unknown, Conflict and ContextRequired distinct, with every original qualifier, citation, source URL and context preserved. No evidence or conflicting evidence never becomes compatibility or activity. Buffer activity is a contextual vendor assay observation, not digestion efficiency; methylation categories are not site-specific permission. Optional edition/hash pin rejects stale catalog claims."
    )]
    pub(crate) fn resolve_enzyme_evidence(
        &self,
        Parameters(args): Parameters<catalog_tools::EvidenceQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        output(catalog_tools::evidence(args)?)
    }
}
