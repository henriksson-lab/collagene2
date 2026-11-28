//! Tracked-material digest preview and its complete documentary validation.
use super::ENGINE_VERSION;
use super::MAX_RECORD_BYTES;
use super::SourceOrigin;
use super::bad;
use super::encoded;
use super::fresh_source;
use super::hex;
use super::limitations;
use super::summary_shape;
use crate::chemistry;
use crate::chemistry::tracked_duplex as modified_duplex;
use crate::chemistry::tracked_duplex::ModificationSummary;
use crate::cloning::restriction::catalog as catalog_digest;
use crate::cloning::restriction::catalog::CatalogSelection;
use crate::cloning::restriction::inspection as restriction_inspection;
use crate::cloning::restriction::inspection::InspectionOptions;
use crate::cloning::restriction::inspection::RestrictionInspection;
use crate::verification;
use collagene_model::Error;
use collagene_model::Molecule;
use serde::Deserialize;
use serde::Serialize;

dto! {
/// Summaries from actual private native source construction/replay, never input authority.
pub struct DigestSourceMaterial{pub duplex:chemistry::DuplexSummary,pub modifications:ModificationSummary}}
dto! {pub struct DigestPreview{pub engine_version:String,pub context_revision:String,pub selection:CatalogSelection,pub inspection:RestrictionInspection,
#[serde(deserialize_with="crate::encoding::required_option")]
#[cfg_attr(feature="json-schema",schemars(required))]
pub source_material:Option<DigestSourceMaterial>,pub fragment_modifications:Vec<ModificationSummary>,pub incomplete:bool,pub findings:Vec<String>,pub findings_total:usize,pub findings_truncated:bool}}
pub fn digest_inspection_revision(
    role: &str,
    source: &Molecule,
    origin: &SourceOrigin,
    selection: &CatalogSelection,
    options: &InspectionOptions,
) -> Result<String, Error> {
    collagene_model::bounds::model_preflight(&[source])?;
    source.validate()?;
    options.validate()?;
    let catalog = catalog_digest::builtin_catalog()?;
    catalog_digest::source_preflight(&catalog, selection, role, source, options)?;
    match origin {
        // Recognition applies to a bounded DNA reference even when no intact
        // physical substrate can be established. Actual admission remains below.
        SourceOrigin::OriginalReference => {}
        SourceOrigin::Saved(origin) => {
            origin.validate_documentary()?;
            if origin.event().revision != source.revision() {
                return Err(bad("inspection origin/source revision differs"));
            }
        }
    }
    digest_inspection_context(role, source, origin, selection, options)
}
/// Current digest declaration identity shared by live admission and historical evidence.
/// Callers establish their own current-catalog or documentary guards first.
pub(crate) fn digest_inspection_context(
    role: &str,
    source: &Molecule,
    origin: &SourceOrigin,
    selection: &CatalogSelection,
    options: &InspectionOptions,
) -> Result<String, Error> {
    crate::encoding::bounded_json::hash(&(
        ENGINE_VERSION,
        "digest-inspection",
        verification::ENGINE_VERSION,
        role,
        source,
        origin,
        selection,
        options,
    ))
}

pub fn inspect_digest(
    role: &str,
    source: &Molecule,
    origin: &SourceOrigin,
    selection: &CatalogSelection,
    options: &InspectionOptions,
) -> Result<DigestPreview, Error> {
    let context_revision = digest_inspection_revision(role, source, origin, selection, options)?;
    let catalog = catalog_digest::builtin_catalog()?;
    let enzymes = selection.resolved_enzymes(&catalog)?;
    let fresh = match fresh_source(source, origin) {
        Ok(fresh) => fresh,
        Err(Error::Unsupported(reason)) if matches!(origin, SourceOrigin::OriginalReference) => {
            let inspection = restriction_inspection::inspect_with_digest(
                role,
                source,
                &enzymes,
                options,
                || Err(Error::Unsupported(reason.clone())),
            )?;
            let (findings, findings_total, findings_truncated) = limitations(
                vec![format!(
                    "Recognition-only reference: actual source material is unavailable ({reason}). No physical products or source chemistry are asserted."
                )],
                1,
                false,
            );
            let report = DigestPreview {
                engine_version: ENGINE_VERSION.into(),
                context_revision,
                selection: selection.clone(),
                inspection,
                source_material: None,
                fragment_modifications: Vec::new(),
                incomplete: true,
                findings,
                findings_total,
                findings_truncated,
            };
            report.validate_binding(&report.context_revision)?;
            return Ok(report);
        }
        Err(error) => return Err(error),
    };
    let mut fragment_modifications = Vec::new();
    let inspection = restriction_inspection::inspect_with_digest_for_source(
        role,
        source,
        &enzymes,
        options,
        &fresh.tracked.duplex.summary()?,
        || {
            let mut result = Vec::new();
            for fragment in modified_duplex::digest(&fresh.tracked, &enzymes)? {
                fragment_modifications.push(modified_duplex::summary(&fragment.duplex)?);
                result.push(chemistry::DigestFragment {
                    duplex: fragment.duplex.duplex,
                    top_start: fragment.top_start,
                    bottom_start: fragment.bottom_start,
                })
            }
            Ok(result)
        },
    )?;
    if !inspection.products_determined {
        fragment_modifications.clear()
    }
    let report = DigestPreview {
        engine_version: ENGINE_VERSION.into(),
        context_revision,
        selection: selection.clone(),
        source_material: Some(DigestSourceMaterial {
            duplex: fresh.tracked.duplex.summary()?,
            modifications: modified_duplex::summary(&fresh.tracked)?,
        }),
        fragment_modifications,
        incomplete: fresh.incomplete || !inspection.scan_complete || !inspection.limits.is_empty(),
        inspection,
        findings: fresh.findings,
        findings_total: fresh.findings_total,
        findings_truncated: fresh.findings_truncated,
    };
    report.validate_binding(&report.context_revision)?;
    Ok(report)
}
impl DigestPreview {
    pub fn selected_fragment(
        &self,
        choice: &crate::cloning::restriction::cut_selection::Choice,
    ) -> Result<&crate::cloning::restriction::inspection::FragmentSummary, Error> {
        self.validate_binding(&self.context_revision)?;
        let material = self.source_material.as_ref().ok_or_else(|| {
            bad("recognition-only reference has no established source material or fragment")
        })?;
        crate::cloning::restriction::cut_selection::selected_fragment_raw_with_source(
            &self.inspection,
            choice,
            &material.duplex,
        )
    }
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        if self.engine_version != ENGINE_VERSION
            || !hex(expected)
            || self.context_revision != expected
            || self.findings.len() != self.findings_total.min(128)
            || self.findings.iter().any(|s| s.chars().count() > 2048)
            || (self.findings_total > self.findings.len() && !self.findings_truncated)
            || self.findings_total > crate::material::lineage::MAX_RUN_BYTES
        {
            return Err(bad("digest preview binding/bounds differ"));
        }
        self.selection.validate()?;
        if self.inspection.enzymes
            != self
                .selection
                .enzymes
                .iter()
                .map(|pin| pin.enzyme.clone())
                .collect::<Vec<_>>()
        {
            return Err(bad(
                "Inspection enzyme geometry differs from its pinned catalog selection.",
            ));
        }
        let Some(material) = &self.source_material else {
            self.inspection.validate()?;
            if !self.incomplete
                || self.inspection.products_determined
                || !self.inspection.fragments.is_empty()
                || !self.fragment_modifications.is_empty()
            {
                return Err(bad(
                    "recognition-only source cannot assert actual material or fragments",
                ));
            }
            encoded(self, MAX_RECORD_BYTES)?;
            return Ok(());
        };
        summary_shape(&material.modifications)?;
        self.inspection.validate_with_source(&material.duplex)?;
        crate::cloning::restriction::cut_selection::source_geometry_preflight(
            &self.inspection,
            &material.duplex,
        )?;
        if material.modifications.top_length != self.inspection.source_length
            || material.modifications.bottom_length != material.duplex.bottom_length
        {
            return Err(bad(
                "Source chemistry summaries differ from the actual source geometry.",
            ));
        }
        if self.fragment_modifications.len()
            != if self.inspection.products_determined {
                self.inspection.fragments.len()
            } else {
                0
            }
        {
            return Err(bad(
                "fragment chemistry summaries differ from actual inspection",
            ));
        }
        for (s, f) in self
            .fragment_modifications
            .iter()
            .zip(&self.inspection.fragments)
        {
            summary_shape(s)?;
            if s.top_length != f.duplex.top_length
                || s.bottom_length != f.duplex.bottom_length
                || s.generation != material.modifications.generation
            {
                return Err(bad("fragment summary coordinate lengths differ"));
            }
        }
        encoded(self, MAX_RECORD_BYTES)?;
        Ok(())
    }
}
