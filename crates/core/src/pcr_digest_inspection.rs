//! Fresh exact-PCR followed by pinned complete-digest inspection. Displayed
//! fragments are documentary choices, never serialized physical stock.
use crate::{
    Error, Molecule, Topology,
    analysis::Primer,
    assembly_recipe,
    catalog_digest::{self, CatalogSelection},
    chemistry,
    restriction_inspection::{self, InspectionOptions, RestrictionInspection},
    verification,
};
use serde::{Deserialize, Serialize};
pub const ENGINE_VERSION: &str = "pcr-digest-inspection-1";
const MAX_GEOMETRIES: usize = 10_000;
macro_rules! dto {($($x:tt)*)=>{#[derive(Clone,Debug,Serialize,Deserialize)]#[cfg_attr(feature="json-schema",derive(schemars::JsonSchema))]#[serde(deny_unknown_fields)]$($x)*};}
dto! {pub struct Request {pub template_strand:chemistry::Strand,pub source_role:String,pub forward:Primer,pub reverse:Primer,pub selection:CatalogSelection,pub options:InspectionOptions}}
dto! {pub struct GeometrySummary {pub template_strand:chemistry::Strand,pub index:usize,pub forward_start:usize,pub reverse_start:usize,pub product_length:usize}}
dto! {pub struct Preview {pub engine_version:String,pub context_revision:String,pub source_revision:String,pub request:Request,pub pcr_geometries:Vec<GeometrySummary>,pub amplicon:Option<Molecule>,pub inspection:Option<RestrictionInspection>,pub findings:Vec<String>}}
fn bad(message: &str) -> Error {
    Error::Unsupported(format!("PCR digest inspection: {message}"))
}
fn hex(value: &str) -> bool {
    value.len() == 64
        && value
            .bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}
pub(crate) fn shape(request: &Request) -> Result<(), Error> {
    if request.source_role.is_empty()
        || request.source_role.len() > 128
        || request.source_role.chars().any(char::is_control)
    {
        return Err(bad("invalid source role"));
    }
    for p in [&request.forward, &request.reverse] {
        if p.name.is_empty()
            || p.name.len() > 256
            || p.name.chars().any(char::is_control)
            || p.annealing.is_empty()
            || p.annealing.len() > 1000
            || p.tail.len() > 1000
            || p.annealing.alphabet() != crate::Alphabet::Dna
            || !p
                .annealing
                .as_str()
                .bytes()
                .chain(p.tail.bytes())
                .all(|b| b"ACGT".contains(&b))
        {
            return Err(bad(
                "exact-PCR primers require bounded canonical concrete DNA",
            ));
        }
    }
    request.options.validate()?;
    request.selection.validate()?;
    Ok(())
}
pub fn source_preflight(source: &Molecule, request: &Request) -> Result<(), Error> {
    preflight(source, request, true)
}
fn preflight(source: &Molecule, request: &Request, reject_imported: bool) -> Result<(), Error> {
    shape(request)?;
    if reject_imported && request.template_strand != chemistry::Strand::Top {
        return Err(bad(
            "Bottom template PCR requires freshly reconstructed native material",
        ));
    }
    crate::model::bounds::model_preflight(&[source])?;
    if reject_imported && chemistry::has_unmapped_physical_evidence(source) {
        return Err(bad(
            "imported physical intermediates require reconstruction of their original graph",
        ));
    }
    let fragment = assembly_recipe::FragmentSpec {
        generation: None,
        role: request.source_role.clone(),
        input_id: request.source_role.clone(),
        preparation: assembly_recipe::Preparation::PcrCatalogDigest {
            forward: request.forward.clone(),
            reverse: request.reverse.clone(),
            selection: request.selection.clone(),
            fragment_index: 0,
            repair: None,
            template_strand: request.template_strand,
        },
        reverse: false,
    };
    assembly_recipe::preparation_rows_preflight(
        &[(&request.source_role, source)],
        &request.source_role,
        &[fragment],
    )?;
    let catalog = catalog_digest::builtin_catalog()?;
    catalog_digest::source_preflight(
        &catalog,
        &request.selection,
        &request.source_role,
        source,
        &request.options,
    )?;
    Ok(())
}
fn context(source_revision: &str, request: &Request) -> Result<String, Error> {
    assembly_recipe::hash(&(
        ENGINE_VERSION,
        verification::ENGINE_VERSION,
        restriction_inspection::ENGINE_VERSION,
        source_revision,
        request,
    ))
}
pub fn inspection_revision(source: &Molecule, request: &Request) -> Result<String, Error> {
    source_preflight(source, request)?;
    context(&source.revision(), request)
}
pub fn preview_revision(source: &Molecule, request: &Request) -> Result<String, Error> {
    inspection_revision(source, request)
}
pub fn inspect(source: &Molecule, request: &Request) -> Result<Preview, Error> {
    let context_revision = inspection_revision(source, request)?;
    let documents = verification::simulate_pcr_documents(
        source,
        &request.forward,
        &request.reverse,
        request.template_strand,
    )?;
    finish(source, request, context_revision, documents)
}
fn finish(
    source: &Molecule,
    request: &Request,
    context_revision: String,
    documents: Vec<verification::PcrDocument>,
) -> Result<Preview, Error> {
    if documents.len() > MAX_GEOMETRIES {
        return Err(Error::TooLarge);
    }
    let geometries = documents
        .iter()
        .enumerate()
        .map(|(index, p)| GeometrySummary {
            template_strand: p.geometry.template_strand,
            index,
            forward_start: p.geometry.forward_start,
            reverse_start: p.geometry.reverse_start,
            product_length: p.molecule.sequence.len(),
        })
        .collect();
    let mut findings=vec!["Fresh exact-PCR geometry under ideal synthesis followed by assumed-complete cleavage. PCR efficiency, reagent compatibility, methylation susceptibility and experimental success are not established. Newly synthesized outer termini are 5′OH/3′OH; internal restriction cuts retain their actual cleaved chemistry. Displayed fragments are not transferable physical stock.".into()];
    let (amplicon, inspection) = if documents.len() == 1 {
        let document = documents.into_iter().next().unwrap();
        findings.extend(document.findings);
        // Construct only here, from this pass's unique native PCR result.
        let duplex = chemistry::Duplex::exact_pcr(&document.molecule.sequence)?;
        let catalog = catalog_digest::builtin_catalog()?;
        let enzymes = request.selection.resolved_enzymes(&catalog)?;
        let inspection = restriction_inspection::inspect_fresh_duplex(
            &request.source_role,
            &document.molecule,
            &enzymes,
            &request.options,
            &duplex,
        )?;
        (Some(document.molecule), Some(inspection))
    } else {
        findings.push(format!("{} exact PCR geometries: digest products withheld. A unique geometry is required by the existing preparation operation; no geometry was selected implicitly.",documents.len()));
        (None, None)
    };
    let preview = Preview {
        engine_version: ENGINE_VERSION.into(),
        context_revision,
        source_revision: source.revision(),
        request: request.clone(),
        pcr_geometries: geometries,
        amplicon,
        inspection,
        findings,
    };
    preview.validate()?;
    Ok(preview)
}
impl Preview {
    /// Readonly binding/shape checks; no imported receipt becomes fresh science.
    pub fn validate(&self) -> Result<(), Error> {
        shape(&self.request)?;
        if self.engine_version != ENGINE_VERSION
            || !hex(&self.source_revision)
            || self.context_revision != context(&self.source_revision, &self.request)?
            || self.pcr_geometries.len() > MAX_GEOMETRIES
            || self.findings.len() > 10_000
            || self.findings.iter().any(|s| s.len() > 65_536)
        {
            return Err(bad("invalid preview context or bounds"));
        }
        for (i, g) in self.pcr_geometries.iter().enumerate() {
            if g.index != i
                || g.template_strand != self.request.template_strand
                || g.product_length == 0
                || g.product_length > crate::model::MAX_SEQUENCE_LEN
                || g.forward_start > crate::model::MAX_SEQUENCE_LEN
                || g.reverse_start > crate::model::MAX_SEQUENCE_LEN
            {
                return Err(bad("invalid geometry summary"));
            }
        }
        match (&self.amplicon, &self.inspection) {
            (Some(m), Some(report)) if self.pcr_geometries.len() == 1 => {
                crate::model::bounds::model_preflight(&[m])?;
                m.validate()?;
                if m.topology != Topology::Linear
                    || m.sequence.len() != self.pcr_geometries[0].product_length
                    || m.metadata.get("pcr:template_revision") != Some(&self.source_revision)
                    || m.metadata.get("pcr:engine_version").map(String::as_str)
                        != Some(verification::ENGINE_VERSION)
                    || m.metadata.get("pcr:primers")
                        != Some(
                            &serde_json::to_string(&(&self.request.forward, &self.request.reverse))
                                .map_err(|_| Error::TooLarge)?,
                        )
                {
                    return Err(bad("amplicon source/primer binding differs"));
                }
                let catalog = catalog_digest::builtin_catalog()?;
                let enzymes = self.request.selection.resolved_enzymes(&catalog)?;
                report.validate_bindings_for(
                    &self.request.source_role,
                    m,
                    &enzymes,
                    &self.request.options,
                )?;
            }
            (None, None) if self.pcr_geometries.len() != 1 => {}
            _ => return Err(bad("only unique PCR may expose digest fragments")),
        }
        assembly_recipe::hash(self)?;
        Ok(())
    }
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        self.validate()?;
        if self.context_revision != expected {
            return Err(bad("stale source/primer/catalog/options context"));
        }
        Ok(())
    }
}

// Native material is attested against its owning archive, then replayed once.
// The inner preview remains documentary and never supplies physical stock.
dto! {pub struct MaterialPreview {pub context_revision:String,pub origin_revision:String,pub origin_incomplete:bool,pub origin_findings_total:usize,pub origin_findings_truncated:bool,pub preview:Preview}}
fn material_context(origin_revision: &str, preview_revision: &str) -> Result<String, Error> {
    assembly_recipe::hash(&(
        "saved-material-pcr-digest-inspection-1",
        ENGINE_VERSION,
        verification::ENGINE_VERSION,
        origin_revision,
        preview_revision,
    ))
}
pub fn preflight_from_material(
    source: &Molecule,
    origin: &crate::run_material::MaterialOrigin,
    request: &Request,
    workspace: &crate::workspace::Workspace,
) -> Result<(), Error> {
    crate::run_material::material_origins_preflight([origin])?;
    preflight(source, request, false)?;
    origin.validate_for_workspace(workspace, origin.event())?;
    if crate::run_material::resolve(workspace, origin.event())? != source {
        return Err(bad("source snapshot differs from its exact retained event"));
    }
    if source.revision() != origin.event().revision {
        return Err(bad("saved source differs from the pinned origin event"));
    }
    Ok(())
}
pub fn inspection_revision_from_material(
    source: &Molecule,
    origin: &crate::run_material::MaterialOrigin,
    request: &Request,
    workspace: &crate::workspace::Workspace,
) -> Result<String, Error> {
    preflight_from_material(source, origin, request, workspace)?;
    material_context(
        &assembly_recipe::hash(origin)?,
        &context(&source.revision(), request)?,
    )
}
pub fn inspect_from_material(
    source: &Molecule,
    origin: &crate::run_material::MaterialOrigin,
    request: &Request,
    workspace: &crate::workspace::Workspace,
) -> Result<MaterialPreview, Error> {
    let context_revision = inspection_revision_from_material(source, origin, request, workspace)?;
    let material = origin.replay(source)?;
    let documents = verification::simulate_pcr_documents_from_fresh(
        source,
        &request.forward,
        &request.reverse,
        material.tracked(),
        request.template_strand,
    )?;
    let mut preview = finish(
        source,
        request,
        context(&source.revision(), request)?,
        documents,
    )?;
    preview.findings.push("Template material was freshly reconstructed from its exact retained native origin; historical assumptions remain applicable.".into());
    if material.incomplete() {
        preview.findings.push("The native origin remains incomplete; fresh template reconstruction does not establish experimental success or resolve its retained limitations.".into());
    }
    preview.findings.push(format!("Native origin findings: showing {} of {}; excerpt truncated: {}. These retained limitations do not establish reaction success.", material.findings().len(),material.findings_total(),material.findings_truncated()));
    preview.findings.extend(material.findings().iter().cloned());
    let report = MaterialPreview {
        context_revision,
        origin_revision: assembly_recipe::hash(origin)?,
        origin_incomplete: material.incomplete(),
        origin_findings_total: material.findings_total(),
        origin_findings_truncated: material.findings_truncated(),
        preview,
    };
    report.validate_binding(&report.context_revision)?;
    Ok(report)
}
impl MaterialPreview {
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        self.preview.validate()?;
        if self.origin_findings_total > 32 * 1024 * 1024
            || (self.origin_findings_total > 128 && !self.origin_findings_truncated)
            || !hex(&self.origin_revision)
            || self.context_revision != expected
            || self.context_revision
                != material_context(&self.origin_revision, &self.preview.context_revision)?
        {
            return Err(bad("stale saved-origin/source/primer/catalog context"));
        }
        Ok(())
    }
}
