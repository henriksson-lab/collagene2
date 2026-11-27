//! Target-free digest inspection of freshly reconstructed overlap-extension PCR.
//! Imported summaries never restore physical state; saved origins are replayed.
use crate::{
    Error, Molecule, assembly_recipe,
    catalog_digest::{self, CatalogSelection},
    chemistry,
    modified_duplex::TrackedDuplex,
    overlap_extension::Geometry,
    overlap_extension_pcr,
    pcr_digest_inspection::GeometrySummary,
    restriction_inspection::{self, InspectionOptions, RestrictionInspection},
    run_lineage::EventBinding,
    run_material::{MaterialOrigin, MaterialSelection},
    verification,
};
use serde::{Deserialize, Serialize};
use std::io::Write;
pub const ENGINE_VERSION: &str = "overlap-extension-digest-inspection-2";
const MAX_BYTES: usize = 32 << 20;
macro_rules! dto{($($t:tt)*)=>{#[derive(Clone,Debug,Serialize,Deserialize)]#[cfg_attr(feature="json-schema",derive(schemars::JsonSchema))]#[serde(deny_unknown_fields)]$($t)*};}
dto! {pub struct Request {pub source_role:String,pub preparation:overlap_extension_pcr::Request,pub selection:CatalogSelection,pub options:InspectionOptions}}
dto! {pub struct PcrStage {pub stage:String,pub template_revision:String,pub geometries:Vec<GeometrySummary>,pub product_revision:Option<String>,pub product_length:Option<usize>}}
dto! {pub struct OriginSummary {pub selection:MaterialSelection,pub source:EventBinding,pub material_revision:String,pub incomplete:bool,pub findings:Vec<String>,pub findings_total:usize,pub findings_truncated:bool}}
dto! {pub struct Preview {pub engine_version:String,pub context_revision:String,pub source_revision:String,pub origin_revision:Option<String>,pub origin_incomplete:Option<bool>,pub request:Request,pub pcr_stages:Vec<PcrStage>,pub extension:Option<Geometry>,pub origin_evidence:Option<OriginSummary>,pub amplicon:Option<Molecule>,pub inspection:Option<RestrictionInspection>,pub incomplete:bool,pub findings:Vec<String>}}
fn bad(message: &str) -> Error {
    Error::Unsupported(format!("overlap-extension digest inspection: {message}"))
}
fn hex(s: &str) -> bool {
    s.len() == 64
        && s.bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}
struct Count(usize);
impl Write for Count {
    fn write(&mut self, b: &[u8]) -> std::io::Result<usize> {
        self.0 = self
            .0
            .checked_add(b.len())
            .filter(|n| *n <= MAX_BYTES)
            .ok_or_else(|| std::io::Error::other("bounded extension digest evidence"))?;
        Ok(b.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}
fn size(value: &impl Serialize) -> Result<(), Error> {
    serde_json::to_writer(Count(0), value).map_err(|_| Error::TooLarge)
}
fn shape(request: &Request) -> Result<(), Error> {
    if request.source_role.trim().is_empty()
        || request.source_role.len() > 128
        || request.source_role.chars().any(char::is_control)
    {
        return Err(bad("invalid bounded source role"));
    }
    request.preparation.validate_shape()?;
    request.selection.validate()?;
    request.options.validate()
}
pub fn source_preflight(
    source: &Molecule,
    origin: Option<&MaterialOrigin>,
    request: &Request,
    workspace: Option<&crate::workspace::Workspace>,
) -> Result<(), Error> {
    shape(request)?;
    crate::model::bounds::model_preflight(&[source])?;
    crate::run_material::material_origins_preflight(origin)?;
    size(&(source, origin, request))?;
    if let Some(origin) = origin {
        let workspace = workspace
            .ok_or_else(|| bad("saved template requires its exact owning native workspace"))?;
        origin.validate_for_workspace(workspace, origin.event())?;
        if source.revision() != origin.event().revision
            || crate::run_material::resolve(workspace, origin.event())? != source
        {
            return Err(bad(
                "saved source differs from the exact pinned origin event",
            ));
        }
    } else if chemistry::has_unmapped_physical_evidence(source)
        || source
            .metadata
            .get(crate::recombination::REFERENCE_ONLY_KEY)
            .is_some()
    {
        return Err(bad(
            "imported physical/reference-only material requires reconstruction of its originating graph",
        ));
    }
    if origin.is_none()
        && [&request.preparation.left, &request.preparation.right]
            .iter()
            .any(|p| p.template_strand != chemistry::Strand::Top)
    {
        return Err(bad(
            "Bottom PCR frames require exact freshly reconstructed saved material",
        ));
    }
    // The borrowed guards ran before these bounded request-field copies. There
    // is no intended target and no copied/sliced original declared as new stock.
    let row = assembly_recipe::FragmentSpec {
        generation: None,
        role: request.source_role.clone(),
        input_id: request.source_role.clone(),
        preparation: assembly_recipe::Preparation::OverlapExtensionPcrCatalogDigest {
            request: request.preparation.clone(),
            selection: request.selection.clone(),
            fragment_index: 0,
            repair: None,
        },
        reverse: false,
    };
    assembly_recipe::preparation_rows_preflight(
        &[(&request.source_role, source)],
        &request.source_role,
        &[row],
    )?;
    let catalog = catalog_digest::builtin_catalog()?;
    catalog_digest::source_preflight(
        &catalog,
        &request.selection,
        &request.source_role,
        source,
        &request.options,
    )
}
fn context(
    source_revision: &str,
    origin_revision: &Option<String>,
    request: &Request,
) -> Result<String, Error> {
    assembly_recipe::hash(&(
        ENGINE_VERSION,
        verification::ENGINE_VERSION,
        restriction_inspection::ENGINE_VERSION,
        source_revision,
        origin_revision,
        request,
    ))
}
pub fn inspection_revision(
    source: &Molecule,
    origin: Option<&MaterialOrigin>,
    request: &Request,
    workspace: Option<&crate::workspace::Workspace>,
) -> Result<String, Error> {
    source_preflight(source, origin, request, workspace)?;
    context(
        &source.revision(),
        &origin.map(assembly_recipe::hash).transpose()?,
        request,
    )
}
pub fn inspect(
    source: &Molecule,
    origin: Option<&MaterialOrigin>,
    request: &Request,
    workspace: Option<&crate::workspace::Workspace>,
) -> Result<Preview, Error> {
    let context_revision = inspection_revision(source, origin, request, workspace)?;
    let mut findings=vec!["Fresh ideal exact-PCR branches, selected reciprocal extension and outer PCR, followed by the pinned assumed-complete digest. There is no intended target input. Geometry uniqueness is conditional on this model and is not experimental isolation or reaction success.".into()];
    let material = origin.map(|origin| origin.replay(source)).transpose()?;
    let origin_evidence = origin
        .zip(material.as_ref())
        .map(|(origin, material)| {
            Ok::<_, Error>(OriginSummary {
                selection: origin.selection(),
                source: origin.event().clone(),
                material_revision: crate::modified_duplex::summary(material.tracked())?
                    .state_revision,
                incomplete: material.incomplete(),
                findings: material.findings().to_vec(),
                findings_total: material.findings_total(),
                findings_truncated: material.findings_truncated(),
            })
        })
        .transpose()?;
    let origin_incomplete = material.as_ref().map(|m| m.incomplete());
    let mut incomplete = origin_incomplete.unwrap_or(false);
    if let Some(material) = &material {
        findings.push("Source strands and chemistry were freshly reconstructed from the exact retained native Run/Derived event; inherited limitations remain applicable. This does not establish experimental isolation.".into());
        findings.extend(material.findings().iter().map(|s| format!("Origin: {s}")));
    }
    let executed = crate::overlap_extension_execution::execute(
        source,
        material.as_ref().map(|m| m.tracked()),
        &request.preparation,
    )?;
    incomplete |= executed.incomplete;
    findings.extend(executed.findings);
    let pcr_stages = executed.stages;
    let extension = executed.extension;
    let amplicon = executed.product;
    let mut inspection = None;
    if let Some(product) = &amplicon {
        if incomplete {
            findings.push("The selected actual digest is conditional: originating/terminal-pair uncertainty remains incomplete. Selecting a fragment never resolves alternative pairings or upgrades whole-plan status.".into());
        }
        let physical = TrackedDuplex::exact_pcr(&product.sequence)?;
        let catalog = catalog_digest::builtin_catalog()?;
        let enzymes = request.selection.resolved_enzymes(&catalog)?;
        let result = restriction_inspection::inspect_fresh_duplex(
            &request.source_role,
            product,
            &enzymes,
            &request.options,
            &physical.duplex,
        )?;
        incomplete |= !result.scan_complete || !result.products_determined;
        inspection = Some(result);
    }
    let preview = Preview {
        engine_version: ENGINE_VERSION.into(),
        context_revision,
        source_revision: source.revision(),
        origin_revision: origin.map(assembly_recipe::hash).transpose()?,
        origin_incomplete,
        request: request.clone(),
        pcr_stages,
        extension,
        origin_evidence,
        amplicon,
        inspection,
        incomplete,
        findings,
    };
    preview.validate()?;
    Ok(preview)
}
impl Preview {
    /// Bounded documentary correspondence only. No PCR, source-origin replay,
    /// polymerase extension or digest is performed by this validator.
    pub fn validate(&self) -> Result<(), Error> {
        shape(&self.request)?;
        if self.engine_version != ENGINE_VERSION
            || !hex(&self.context_revision)
            || !hex(&self.source_revision)
            || self.context_revision
                != context(&self.source_revision, &self.origin_revision, &self.request)?
            || self.pcr_stages.len() < 2
            || self.pcr_stages.len() > 3
            || self.findings.len() > 1000
            || self.findings.iter().any(|s| s.len() > 8192)
        {
            return Err(bad(
                "invalid preview context, source or bounded stage evidence",
            ));
        }
        crate::overlap_extension_execution::validate_evidence(
            &self.source_revision,
            &self.request.preparation,
            &self.pcr_stages,
            &self.extension,
            &self.amplicon,
        )?;
        size(self)?;
        match (
            &self.origin_revision,
            self.origin_incomplete,
            &self.origin_evidence,
        ) {
            (None, None, None) => {}
            (Some(revision), Some(incomplete), Some(e))
                if hex(revision)
                    && e.incomplete == incomplete
                    && hex(&e.material_revision)
                    && e.source.revision == self.source_revision
                    && e.selection.event() == &e.source
                    && e.findings.len() <= 128
                    && e.findings.iter().all(|s| s.chars().count() <= 2048)
                    && e.findings_total >= e.findings.len()
                    && (e.findings_truncated || e.findings_total == e.findings.len()) =>
            {
                e.selection.validate()?;
            }
            _ => return Err(bad("saved-origin documentary binding changed")),
        }
        if self.amplicon.is_none() && self.inspection.is_some() {
            return Err(bad("digest fragments require a native singleton outer PCR"));
        }
        let pipeline_available = self.extension.is_some() && self.amplicon.is_some();
        let inherited_incomplete = self.origin_incomplete == Some(true)
            || self
                .extension
                .as_ref()
                .is_none_or(|g| !g.terminal_pairings.complete || !g.terminal_pairings.unique);
        if let Some(report) = &self.inspection {
            if !pipeline_available {
                return Err(bad(
                    "digest evidence requires actual singleton PCRs and the explicit selected extension",
                ));
            }
            let catalog = catalog_digest::builtin_catalog()?;
            let enzymes = self.request.selection.resolved_enzymes(&catalog)?;
            report.validate_bindings_for(
                &self.request.source_role,
                self.amplicon
                    .as_ref()
                    .ok_or_else(|| bad("missing amplicon"))?,
                &enzymes,
                &self.request.options,
            )?;
        }
        let incomplete = inherited_incomplete
            || !pipeline_available
            || self
                .inspection
                .as_ref()
                .is_none_or(|r| !r.scan_complete || !r.products_determined);
        if self.incomplete != incomplete {
            return Err(bad("inherited or capped incompleteness was changed"));
        }
        size(self)
    }
    pub fn validate_binding(&self, revision: &str) -> Result<(), Error> {
        if self.context_revision != revision {
            return Err(bad("stale source/preparation/catalog/origin context"));
        }
        self.validate()
    }
    pub fn selectable(&self) -> bool {
        self.amplicon.is_some()
            && self.inspection.as_ref().is_some_and(|r| {
                r.scan_complete && r.products_determined && !r.fragments.is_empty()
            })
    }
}
