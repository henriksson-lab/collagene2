//! Readable rendering of the same fresh source-format advice in HTML and PDF.
use crate::{
    cloning::assembly::standards::declaration::DeclaredReviewState, planning::plan::DesignPlan,
    verification::evidence::VerificationReport,
};
use collagene_model::Error;

pub(crate) fn lines(
    plan: &DesignPlan,
    verification: &VerificationReport,
) -> Result<Vec<String>, Error> {
    verification.validate_standard_for_plan(plan)?;
    let Some(review) = &verification.standard_review else {
        return Ok(vec![]);
    };
    let mut lines = vec![
        "Declared source-part standard — advisory format review".into(),
        review.scope.clone(),
        format!(
            "Plan {}; declaration {}; review {}.",
            review.plan_revision, review.declaration_revision, review.review_revision
        ),
        format!(
            "Comparison allowance: requested {}, effective {:?}, remaining before {}; charge {:?}. Reservations after refusal are not measured comparisons.",
            review.requested_limits.max_comparisons,
            review.effective_limits.as_ref().map(|l| l.max_comparisons),
            review.shared_remaining_before,
            review.charge
        ),
    ];
    match &review.state {
        DeclaredReviewState::Unavailable { reasons } => {
            lines.push(
                "Source-format advice unavailable; no scoped format pass is established.".into(),
            );
            for reason in reasons {
                lines.push(format!(
                    "Part {:?}: {:?}. {}",
                    reason.part_index.map(|i| i + 1),
                    reason.kind,
                    reason.message
                ));
            }
        }
        DeclaredReviewState::Assessed(report) => {
            lines.push(format!(
                "Profile {} ({}) identity {}; data revision {}; catalogue {} / {}.",
                report.profile.name,
                report.profile.profile_id,
                report.profile.profile_identity,
                report.profile.definition.revision,
                report.profile.catalog_pin.edition,
                report.profile.catalog_pin.content_hash
            ));
            lines.push(format!("Scoped format outcome {:?}; complete {}; observed violations {}, unknowns {}, advisories {}; budget hits {:?}. This outcome is not full standard or physical acceptance.", report.assessment.outcome, report.assessment.complete, report.assessment.violations_observed, report.assessment.indeterminate_observed, report.assessment.advisories_observed, report.assessment.budgets_hit));
            lines.push(format!(
                "Native request {}; assessment {}; native report {}.",
                report.context_revision, report.assessment.request_revision, report.report_revision
            ));
            lines.push(report.profile.scope.clone());
            for (index, (source, part)) in report
                .sources
                .iter()
                .zip(&report.assessment.parts)
                .enumerate()
            {
                lines.push(format!("Part {} role {}: source revision {}, source start {} length {} reverse {}; slot {}; explicit oriented body [{}..{}).", index+1, source.id, source.revision, source.start, source.length, source.reverse, part.slot, part.body_start, part.body_end));
            }
            for finding in &report.assessment.findings {
                lines.push(format!(
                    "Part {:?}: {:?}; {:?}.",
                    finding.part_index.map(|i| i + 1),
                    finding.severity,
                    finding.kind
                ));
            }
            lines.push(format!(
                "Legal next declared roles {:?}; may finish {}; {}",
                report.successors.allowed_slots,
                report.successors.may_finish,
                report.successors.scope
            ));
            lines.extend(report.profile.definition.limitations.iter().cloned());
            for source in &report.profile.definition.sources {
                lines.push(format!("Source {}: {}; revision {}; locator {}; reviewed {}; rights {}; capture SHA256 {:?}; {}.", source.id, source.url, source.revision, source.locator, source.reviewed_on, source.source_license, source.upstream_sha256, source.snapshot_note));
            }
        }
    }
    Ok(lines)
}
