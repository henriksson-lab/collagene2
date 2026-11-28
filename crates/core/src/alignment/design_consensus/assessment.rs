//! One actual pairwise comparison followed by bounded selected-path semantics.
use super::*;
use crate::{
    alignment::review, documents::workspace::Workspace, encoding::bounded_json, material::lineage,
};
use collagene_model::Error;
pub(super) fn assess(workspace: &Workspace, request: &Request) -> Result<Report, Error> {
    let context_revision = super::revision(workspace, request)?;
    let (design, consensus) = request.sources(workspace)?;
    let pairwise = review::compare(
        design,
        consensus,
        &review::Options::Whole(request.comparison),
    )?;
    let mut work = Work::new(request.limits.clone());
    let index = if request.cds.is_empty() {
        None
    } else {
        Some(projection::Index::new(&pairwise, &mut work))
    };
    let mut cds = Vec::with_capacity(request.cds.len());
    for selection in &request.cds {
        let code = crate::construct_qc::CodePin::for_id(selection.genetic_code)?;
        let row = match &index {
            Some(Ok(index)) => coding::assess(
                selection, code, request, design, consensus, &pairwise, index, &mut work,
            ),
            _ => CdsReport::unavailable(
                selection,
                code,
                "Selected-path index could not be admitted within the declared work/memory budget; no coding absence claim is made",
            ),
        };
        cds.push(row);
    }
    let coding_complete = cds
        .iter()
        .all(|row| row.mapping_complete && row.register_complete && row.consequences_complete);
    let mut evidence = Evidence {
        engine_version: ENGINE_VERSION.into(),
        context_revision,
        workspace_fingerprint: lineage::workspace_fingerprint(workspace)?,
        request_revision: bounded_json::hash(request)?,
        pairwise_context_revision: String::new(),
        pairwise_report_revision: String::new(),
        request: request.clone(),
        provenance: vec![
            provenance::capture(SourceRole::Design, &request.design, design),
            provenance::capture(SourceRole::Consensus, &request.consensus, consensus),
        ],
        cds,
        work: work.summary,
        coding_complete,
        selected_path_conditional: true,
        statement: STATEMENT.into(),
        evidence_revision: String::new(),
    };
    evidence.seal(&pairwise)?;
    let report = Report { pairwise, evidence };
    report.validate_binding(&report.evidence.context_revision)?;
    Ok(report)
}
