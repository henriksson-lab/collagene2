//! Validation shared by inbound reports and execution preflight.
use collagene_model::Error;
pub(super) fn validate_verification_report(
    report: &crate::verification::evidence::VerificationReport,
) -> Result<(), Error> {
    if report.step_reports.len() > 1000 {
        return Err(Error::TooLarge);
    }
    let mut seen_steps = std::collections::BTreeSet::new();
    for step in &report.step_reports {
        step.validate_documentary()?;
        if !seen_steps.insert(step.step_index) {
            return Err(Error::Format("duplicate retained step evidence".into()));
        }
    }
    if !report.reference_step_reports.is_empty() {
        if report.reference_step_reports.len() > 64 {
            return Err(Error::TooLarge);
        }
        let mut seen = std::collections::BTreeSet::new();
        for step in &report.reference_step_reports {
            if step.step_index >= 1000 || !seen.insert(step.step_index) {
                return Err(Error::Format("Invalid reference step index".into()));
            }
            match &step.evidence {
                crate::verification::evidence::ReferenceStepEvidence::MultiSiteRecombination {
                    selection,
                    batch,
                } => {
                    batch.validate_binding(&selection.batch_revision)?;
                    if !batch.documents.iter().any(|p| p.selection == *selection) {
                        return Err(Error::Format(
                            "Reference selection is missing from its batch".into(),
                        ));
                    }
                }
            }
        }
    }
    if !report.physical_processing_reports.is_empty() {
        if report.physical_processing_reports.len() > 64 {
            return Err(Error::TooLarge);
        }
        let mut seen = std::collections::BTreeSet::new();
        for step in &report.physical_processing_reports {
            if step.step_index >= 1000 || !seen.insert(step.step_index) {
                return Err(Error::Format(
                    "Invalid physical processing step index".into(),
                ));
            }
            match &step.evidence {
                crate::verification::evidence::PhysicalProcessingEvidence::LigationOutcomes(r) => {
                    r.validate_documentary()?;
                }
                crate::verification::evidence::PhysicalProcessingEvidence::DerivedOrigin(r) => {
                    r.validate_documentary()?;
                }
                crate::verification::evidence::PhysicalProcessingEvidence::RunProduct(r) => {
                    r.validate_documentary()?;
                }
                crate::verification::evidence::PhysicalProcessingEvidence::TopoActions(r) => {
                    r.validate_binding(&r.context_revision)?
                }
                crate::verification::evidence::PhysicalProcessingEvidence::OverlapExtension(r) => {
                    r.validate_documentary()?;
                }
                crate::verification::evidence::PhysicalProcessingEvidence::ThreePrimeResection(
                    report,
                ) => {
                    report.validate_documentary()?;
                }
                crate::verification::evidence::PhysicalProcessingEvidence::HomologyAnnealing(
                    report,
                ) => {
                    report.validate_documentary()?;
                }
                crate::verification::evidence::PhysicalProcessingEvidence::HomologyAssembly(
                    report,
                ) => {
                    report.validate_documentary()?;
                }
            }
        }
    }

    Ok(())
}
