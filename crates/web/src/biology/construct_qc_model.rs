//! Committed UI declarations. Annotation choices are assumptions, not receipts.
use collagene_core::construct_qc as native;
use collagene_model::Molecule;

#[derive(Clone, PartialEq, Eq, Default)]
pub(crate) struct Draft {
    pub options: native::Options,
    pub assumed: Vec<usize>,
}
impl Draft {
    pub fn options_for(&self, m: &Molecule) -> Result<native::Options, String> {
        let mut options = self.options.clone();
        options.annotations.clear();
        let source_revision = m.revision();
        for &index in &self.assumed {
            let feature = m
                .features
                .get(index)
                .ok_or("An annotation selection is stale.")?;
            if feature.needs_review {
                return Err("An annotation marked for review cannot be accepted as a current role assumption.".into());
            }
            options.annotations.push(native::AnnotationAssumption {feature_index:index,source_revision:source_revision.clone(),
                annotation_revision:native::annotation_revision(feature).map_err(|e|e.to_string())?,
                rationale:"User explicitly declared this current annotation role and its coordinates for structural checks; biological identity is not verified.".into()});
        }
        options.validate(m).map_err(|e| e.to_string())?;
        Ok(options)
    }
}
