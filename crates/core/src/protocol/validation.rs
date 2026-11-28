//! Validation shared by inbound reports and execution preflight.
use collagene_model::Error;
pub(super) fn validate_verification_report(
    report: &crate::verification::VerificationReport,
) -> Result<(), Error> {
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
                crate::verification::ReferenceStepEvidence::MultiSiteRecombination {
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
                crate::verification::PhysicalProcessingEvidence::LigationOutcomes(r) => {
                    r.validate_documentary()?;
                }
                crate::verification::PhysicalProcessingEvidence::DerivedOrigin(r) => {
                    r.validate_documentary()?;
                }
                crate::verification::PhysicalProcessingEvidence::RunProduct(r) => {
                    r.validate_documentary()?;
                }
                crate::verification::PhysicalProcessingEvidence::TopoActions(r) => {
                    r.validate_binding(&r.context_revision)?
                }
                crate::verification::PhysicalProcessingEvidence::OverlapExtension(r) => {
                    r.validate_documentary()?;
                }
                crate::verification::PhysicalProcessingEvidence::ThreePrimeResection(report) => {
                    report.validate_documentary()?;
                }
                crate::verification::PhysicalProcessingEvidence::HomologyAnnealing(report) => {
                    report.validate_documentary()?;
                }
                crate::verification::PhysicalProcessingEvidence::HomologyAssembly(report) => {
                    report.validate_documentary()?;
                }
            }
        }
    }

    Ok(())
}
