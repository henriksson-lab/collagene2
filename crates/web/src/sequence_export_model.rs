//! Pure context and approval checks for worker-prepared reference exports.
use collagene_core::{Molecule, gcg_sequence, sequence_export};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Mode {
    SnapGene,
    Xdna,
    Gcg,
}
impl Mode {
    pub fn format(self, id: &str, date: &str) -> sequence_export::Format {
        match self {
            Self::SnapGene => sequence_export::Format::SnapGene,
            Self::Xdna => sequence_export::Format::Xdna,
            Self::Gcg => sequence_export::Format::Gcg {
                options: gcg_sequence::ExportOptions {
                    record_id: Some(id.into()),
                    date_label: Some(date.into()),
                },
            },
        }
    }
    pub fn label(self) -> &'static str {
        match self {
            Self::SnapGene => "SnapGene",
            Self::Xdna => "XDNA",
            Self::Gcg => "GCG",
        }
    }
}
#[derive(Debug, Clone)]
pub struct Capture {
    pub source_revision: String,
    pub context_revision: String,
    pub format: sequence_export::Format,
}
impl Capture {
    pub fn new(source: &Molecule, format: sequence_export::Format) -> Result<Self, String> {
        Ok(Self {
            source_revision: source.revision(),
            context_revision: sequence_export::revision(source, &format)
                .map_err(|e| e.to_string())?,
            format,
        })
    }
    /// Request identity includes all options, rather than filename or sequence alone.
    pub fn matches(&self, source: &Molecule, format: &sequence_export::Format) -> bool {
        source.revision() == self.source_revision
            && sequence_export::revision(source, format).ok().as_ref()
                == Some(&self.context_revision)
    }
}
#[derive(Debug, Clone)]
pub struct Review {
    pub capture: Capture,
    pub assessment: sequence_export::Assessment,
}
impl Review {
    pub fn new(capture: Capture, assessment: sequence_export::Assessment) -> Result<Self, String> {
        assessment
            .validate_binding(&capture.context_revision)
            .map_err(|e| e.to_string())?;
        if assessment.source_revision != capture.source_revision {
            return Err("Export assessment belongs to another source revision.".into());
        }
        Ok(Self {
            capture,
            assessment,
        })
    }
    pub fn allow_prepare(&self, acknowledged: bool) -> Result<(), String> {
        if !self.assessment.losses.is_empty() && !acknowledged {
            return Err(
                "Review and explicitly accept the listed omissions before preparing this export."
                    .into(),
            );
        }
        Ok(())
    }
    pub fn validate_prepared(
        &self,
        prepared: &sequence_export::Prepared,
        acknowledged: bool,
    ) -> Result<(), String> {
        self.allow_prepare(acknowledged)?;
        prepared
            .validate_binding(&self.capture.context_revision)
            .map_err(|e| e.to_string())?;
        if serde_json::to_vec(&prepared.assessment).map_err(|e| e.to_string())?
            != serde_json::to_vec(&self.assessment).map_err(|e| e.to_string())?
        {
            return Err("Prepared export differs from the exact reviewed assessment.".into());
        }
        Ok(())
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    use collagene_core::{Topology, sequence_export};
    fn source() -> Molecule {
        Molecule::dna("literal", "ACGTACGT", Topology::Circular).unwrap()
    }
    #[test]
    fn consent_and_prepared_bytes_remain_bound_to_format_options_and_source() {
        let source = source();
        let format = Mode::Gcg.format("literal", "2026-11-26");
        let capture = Capture::new(&source, format.clone()).unwrap();
        let assessment = sequence_export::review(&source, &format).unwrap();
        let review = Review::new(capture, assessment).unwrap();
        assert!(!review.assessment.losses.is_empty());
        assert!(review.allow_prepare(false).is_err());
        let prepared =
            sequence_export::prepare(&source, &format, &review.capture.context_revision, true)
                .unwrap();
        review.validate_prepared(&prepared, true).unwrap();
        assert!(review.validate_prepared(&prepared, false).is_err());
        assert!(
            !review
                .capture
                .matches(&source, &Mode::Gcg.format("changed", "2026-11-26"))
        );
        assert!(
            !review
                .capture
                .matches(&source, &Mode::Gcg.format("literal", "2026-11-27"))
        );
        assert!(
            !review
                .capture
                .matches(&source, &Mode::Xdna.format("literal", "2026-11-26"))
        );
        let mut edited = source.clone();
        edited.name = "renamed".into();
        assert!(!review.capture.matches(&edited, &format));
        let mut changed = prepared.clone();
        changed.assessment.filename = "other.gcg".into();
        assert!(review.validate_prepared(&changed, true).is_err());
        let wrong = Capture::new(&source, Mode::Xdna.format("", "")).unwrap();
        assert!(Review::new(wrong, review.assessment).is_err());
    }
    #[test]
    fn compatible_snapgene_and_xdna_are_distinct_reviewed_profiles() {
        let source = source();
        let before = source.clone();
        for mode in [Mode::SnapGene, Mode::Xdna] {
            let format = mode.format("", "");
            let capture = Capture::new(&source, format.clone()).unwrap();
            let assessment = sequence_export::review(&source, &format).unwrap();
            let review = Review::new(capture, assessment).unwrap();
            let prepared =
                sequence_export::prepare(&source, &format, &review.capture.context_revision, true)
                    .unwrap();
            review.validate_prepared(&prepared, true).unwrap();
            assert!(!prepared.bytes.is_empty());
        }
        assert_eq!(source, before);
    }
}
