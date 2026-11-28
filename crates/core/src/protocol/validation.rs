//! Documentary validation delegates to the report's single current owner.
use collagene_model::Error;
pub(super) fn validate_verification_report(
    report: &crate::verification::evidence::VerificationReport,
) -> Result<(), Error> {
    report.validate_documentary()
}
