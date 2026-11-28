//! Reciprocal polymerase-extension preview and exact product selection.
use super::ENGINE_VERSION;
use super::OverlapExtensionSelection;
use super::SourceOrigin;
use super::bad;
use super::encoded;
use super::fresh_source;
use super::hex;
use super::limitations;
use super::workspace_source;
use crate::chemistry;
use crate::documents::workspace::Workspace;
use crate::material::lineage::EventBinding;
use crate::material::run as run_material;
use crate::verification;
use collagene_model::Error;
use collagene_model::Molecule;
use serde::Deserialize;
use serde::Serialize;

// A bounded source-backed OE reference preview. Its documents and geometry are
// documentary; only native recipe installation/replay establishes material.
dto! {pub struct OverlapExtensionPreview{
pub engine_version:String,pub context_revision:String,pub source_revision:String,
pub source:EventBinding,pub origin_revision:Option<String>,pub origin_incomplete:bool,
pub preparation:crate::cloning::overlap::extension::pcr::Request,
pub pcr_stages:Vec<crate::cloning::overlap::extension::digest::PcrStage>,
#[cfg_attr(feature="json-schema",schemars(with="serde_json::Value"))]
pub extension:Option<crate::cloning::overlap::extension::Geometry>,pub product:Option<Molecule>,
pub incomplete:bool,pub findings:Vec<String>,pub findings_total:usize,pub findings_truncated:bool}}
pub(super) fn overlap_extension_context(
    source: &EventBinding,
    origin: &SourceOrigin,
    preparation: &crate::cloning::overlap::extension::pcr::Request,
) -> Result<String, Error> {
    overlap_extension_context_fields(
        source,
        &match origin {
            SourceOrigin::OriginalReference => None,
            SourceOrigin::Saved(origin) => Some(crate::encoding::bounded_json::hash(origin)?),
        },
        preparation,
    )
}
fn overlap_extension_context_fields(
    source: &EventBinding,
    origin_revision: &Option<String>,
    preparation: &crate::cloning::overlap::extension::pcr::Request,
) -> Result<String, Error> {
    crate::encoding::bounded_json::hash(&(
        ENGINE_VERSION,
        "overlap-extension-preview",
        verification::ENGINE_VERSION,
        source,
        origin_revision,
        preparation,
    ))
}
/// Admission is workspace-bound even for original references: content equality
/// cannot substitute for the source event selected by the user.
pub fn overlap_extension_preview_revision(
    w: &Workspace,
    source: &EventBinding,
    origin: &SourceOrigin,
    preparation: &crate::cloning::overlap::extension::pcr::Request,
) -> Result<String, Error> {
    preparation.validate_shape()?;
    let molecule = workspace_source(w, source, origin)?;
    if matches!(origin, SourceOrigin::OriginalReference)
        && [&preparation.left, &preparation.right]
            .iter()
            .any(|pair| pair.template_strand != chemistry::Strand::Top)
    {
        return Err(bad(
            "Original-reference OE PCR requires Top source strands; Bottom requires native source replay",
        ));
    }
    crate::cloning::overlap::extension::execution::preflight(molecule, preparation)?;
    overlap_extension_context(source, origin, preparation)
}
pub fn preview_overlap_extension(
    w: &Workspace,
    source: &EventBinding,
    origin: &SourceOrigin,
    preparation: &crate::cloning::overlap::extension::pcr::Request,
) -> Result<OverlapExtensionPreview, Error> {
    let context_revision = overlap_extension_preview_revision(w, source, origin, preparation)?;
    let molecule = run_material::resolve(w, source)?;
    let fresh = fresh_source(molecule, origin)?;
    let executed = crate::cloning::overlap::extension::execution::execute(
        molecule,
        Some(&fresh.tracked),
        preparation,
    )?;
    let total = fresh
        .findings_total
        .checked_add(executed.findings.len())
        .ok_or(Error::TooLarge)?;
    let (findings, findings_total, findings_truncated) = limitations(
        executed
            .findings
            .into_iter()
            .chain(fresh.findings)
            .collect(),
        total,
        fresh.findings_truncated,
    );
    let result = OverlapExtensionPreview {
        engine_version: ENGINE_VERSION.into(),
        context_revision,
        source_revision: molecule.revision(),
        source: source.clone(),
        origin_revision: match origin {
            SourceOrigin::OriginalReference => None,
            SourceOrigin::Saved(origin) => Some(crate::encoding::bounded_json::hash(origin)?),
        },
        origin_incomplete: fresh.incomplete,
        preparation: preparation.clone(),
        pcr_stages: executed.stages,
        extension: executed.extension,
        product: executed.product,
        incomplete: fresh.incomplete || executed.incomplete,
        findings,
        findings_total,
        findings_truncated,
    };
    result.validate_binding(&result.context_revision)?;
    Ok(result)
}
impl OverlapExtensionPreview {
    /// No source replay, PCR search or target evaluation occurs here.
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        self.preparation.validate_shape()?;
        run_material::event_shape(&self.source)?;
        if self.engine_version != ENGINE_VERSION
            || !hex(expected)
            || self.context_revision != expected
            || self.source_revision != self.source.revision
            || self.origin_revision.as_ref().is_some_and(|s| !hex(s))
            || self.origin_revision.is_none() && self.origin_incomplete
            || self.context_revision
                != overlap_extension_context_fields(
                    &self.source,
                    &self.origin_revision,
                    &self.preparation,
                )?
            || self.findings.len() != self.findings_total.min(128)
            || self.findings_total > crate::material::lineage::MAX_RUN_BYTES
            || self.findings.iter().any(|s| s.chars().count() > 2048)
            || self.findings_total > self.findings.len() && !self.findings_truncated
        {
            return Err(bad(
                "OE preview source, origin, preparation or bounded context changed",
            ));
        }
        crate::cloning::overlap::extension::execution::validate_evidence(
            &self.source_revision,
            &self.preparation,
            &self.pcr_stages,
            &self.extension,
            &self.product,
        )?;
        let expected_incomplete = self.origin_incomplete
            || self.product.is_none()
            || self
                .extension
                .as_ref()
                .is_none_or(|g| !g.terminal_pairings.complete || !g.terminal_pairings.unique);
        if self.incomplete != expected_incomplete {
            return Err(bad(
                "OE preview lost inherited or terminal-pair uncertainty",
            ));
        }
        encoded(self, crate::material::lineage::MAX_RUN_BYTES)?;
        Ok(())
    }
    /// Conditional source/terminal uncertainty is retained, never upgraded by selection.
    pub fn reusable(&self) -> bool {
        self.product.is_some()
            && self.extension.is_some()
            && self.pcr_stages.len() == 3
            && self
                .pcr_stages
                .iter()
                .all(|p| p.geometries.len() == 1 && p.product_revision.is_some())
    }
    pub fn selection(&self) -> Result<OverlapExtensionSelection, Error> {
        self.validate_binding(&self.context_revision)?;
        if !self.reusable() {
            return Err(bad(
                "No reusable OE outer-PCR product: each source PCR and outer PCR must have one actual geometry",
            ));
        }
        Ok(OverlapExtensionSelection {
            expected_preview_revision: self.context_revision.clone(),
            expected_product_revision: self
                .product
                .as_ref()
                .ok_or_else(|| bad("outer-PCR product absent"))?
                .revision(),
        })
    }
}
