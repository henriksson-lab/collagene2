//! Conditional intentional-mismatch PCR followed by complete-digest geometry.
//! An incomplete singleton remains conditional; no target or unique design is inferred.
use crate::{
    assembly_recipe, catalog_digest, chemistry, mutagenesis, pcr_digest_inspection,
    restriction_inspection::{self, RestrictionInspection},
    run_material::MaterialOrigin,
    verification,
};
use collagene_model::{Error, Molecule, Topology};
use serde::{Deserialize, Serialize};
pub const ENGINE_VERSION: &str = "mismatch-pcr-digest-inspection-2";
macro_rules! dto{($($t:tt)*)=>{#[derive(Clone,Debug,Serialize,Deserialize)]#[cfg_attr(feature="json-schema",derive(schemars::JsonSchema))]#[serde(deny_unknown_fields)]$($t)*};}
dto! {pub struct Request{pub base:pcr_digest_inspection::Request,pub pcr_options:mutagenesis::PcrOptions}}
dto! {pub struct Preview{pub engine_version:String,pub context_revision:String,pub source_revision:String,pub origin_revision:Option<String>,pub origin_incomplete:Option<bool>,pub origin_findings:Vec<String>,pub origin_findings_total:usize,pub origin_findings_truncated:bool,pub request:Request,#[cfg_attr(feature="json-schema",schemars(with="serde_json::Value"))] pub pcr:verification::MutagenesisPreview,pub amplicon:Option<Molecule>,pub inspection:Option<RestrictionInspection>,pub incomplete:bool,pub findings:Vec<String>}}
fn bad(s: &str) -> Error {
    Error::Unsupported(format!("Mismatch PCR digest inspection: {s}"))
}
fn hex(s: &str) -> bool {
    s.len() == 64
        && s.bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}
fn same(a: &impl Serialize, b: &impl Serialize) -> Result<bool, Error> {
    Ok(assembly_recipe::hash(a)? == assembly_recipe::hash(b)?)
}
fn shape(r: &Request) -> Result<(), Error> {
    pcr_digest_inspection::shape(&r.base)?;
    if r.base.template_strand != chemistry::Strand::Top {
        return Err(bad(
            "Mismatch PCR supports the explicit Top frame only; Bottom mismatch seeding is not modeled",
        ));
    }
    assembly_recipe::mismatch_pcr_shape(&r.base.forward, &r.base.reverse, &r.pcr_options)?;
    Ok(())
}
fn context(source: &str, origin: &Option<String>, request: &Request) -> Result<String, Error> {
    assembly_recipe::hash(&(
        ENGINE_VERSION,
        verification::ENGINE_VERSION,
        restriction_inspection::ENGINE_VERSION,
        source,
        origin,
        request,
    ))
}
pub fn source_preflight(
    source: &Molecule,
    origin: Option<&MaterialOrigin>,
    request: &Request,
    workspace: Option<&crate::workspace::Workspace>,
) -> Result<(), Error> {
    shape(request)?;
    collagene_model::bounds::model_preflight(&[source])?;
    source.validate()?;
    if let Some(origin) = origin {
        crate::run_material::material_origins_preflight([origin])?;
        let workspace = workspace
            .ok_or_else(|| bad("the exact owning workspace is required for saved material"))?;
        origin.validate_for_workspace(workspace, origin.event())?;
        if crate::run_material::resolve(workspace, origin.event())? != source {
            return Err(bad("source snapshot differs from its exact retained event"));
        }
        if source.revision() != origin.event().revision {
            return Err(bad("source differs from pinned originating event"));
        }
    } else if chemistry::has_unmapped_physical_evidence(source) {
        return Err(bad(
            "imported physical state requires fresh origin reconstruction",
        ));
    }
    let fragment = assembly_recipe::FragmentSpec {
        generation: None,
        role: request.base.source_role.clone(),
        input_id: request.base.source_role.clone(),
        reverse: false,
        preparation: assembly_recipe::Preparation::MismatchPcrCatalogDigest {
            forward: request.base.forward.clone(),
            reverse: request.base.reverse.clone(),
            options: request.pcr_options.clone(),
            selection: request.base.selection.clone(),
            fragment_index: 0,
            repair: None,
        },
    };
    assembly_recipe::preparation_rows_preflight(
        &[(&request.base.source_role, source)],
        &request.base.source_role,
        &[fragment],
    )?;
    let catalog = catalog_digest::builtin_catalog()?;
    catalog_digest::source_preflight(
        &catalog,
        &request.base.selection,
        &request.base.source_role,
        source,
        &request.base.options,
    )?;
    Ok(())
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
fn selectable(p: &verification::MutagenesisPreview) -> bool {
    p.report.products.len() == 1
        && p.report.products[0].certain
        && !p.report.products[0].overlapping_domains
}
pub fn inspect(
    source: &Molecule,
    origin: Option<&MaterialOrigin>,
    request: &Request,
    workspace: Option<&crate::workspace::Workspace>,
) -> Result<Preview, Error> {
    let context_revision = inspection_revision(source, origin, request, workspace)?;
    let mut origin_findings = Vec::new();
    let mut origin_findings_total = 0;
    let mut origin_findings_truncated = false;
    let (native, origin_incomplete) = if let Some(origin) = origin {
        let material = origin.replay(source)?;
        origin_findings = material.findings().to_vec();
        origin_findings_total = material.findings_total();
        origin_findings_truncated = material.findings_truncated();
        (
            mutagenesis::simulate_pcr_reconstructed(
                source,
                &material.tracked().duplex,
                &request.base.forward,
                &request.base.reverse,
                &request.pcr_options,
            )?,
            Some(material.incomplete()),
        )
    } else {
        (
            mutagenesis::simulate_pcr(
                source,
                &request.base.forward,
                &request.base.reverse,
                &request.pcr_options,
            )?,
            None,
        )
    };
    let pcr = verification::project_mutagenesis_preview(
        source,
        &request.base.forward,
        &request.base.reverse,
        &request.pcr_options,
        native,
    )?;
    let mut findings=vec!["Intentional-mismatch PCR under an explicit binding/search policy, followed by assumed-complete cleavage. Primer-encoded substitutions are conditional sequence predictions; amplification, modification susceptibility, reagent compatibility and laboratory success are not established. New PCR outer termini are 5′OH/3′OH and restriction cuts retain cleaved phosphate chemistry.".into()];
    if origin.is_some() {
        findings.push(format!("Native origin findings: showing {} of {}; excerpt truncated: {}. Retained source limitations are not resolved by this PCR prediction.",origin_findings.len(),origin_findings_total,origin_findings_truncated));
        findings.extend(origin_findings.iter().cloned());
    }
    let (amplicon, inspection) = if selectable(&pcr) {
        let molecule = pcr.documents[0].clone();
        let duplex = chemistry::Duplex::exact_pcr(&molecule.sequence)?;
        let catalog = catalog_digest::builtin_catalog()?;
        let enzymes = request.base.selection.resolved_enzymes(&catalog)?;
        let inspection = restriction_inspection::inspect_fresh_duplex(
            &request.base.source_role,
            &molecule,
            &enzymes,
            &request.base.options,
            &duplex,
        )?;
        (Some(molecule), Some(inspection))
    } else {
        findings.push("No sole certain nonoverlapping PCR geometry is available; digest products are withheld without an implicit candidate choice.".into());
        (None, None)
    };
    let incomplete = !pcr.report.enumeration_complete
        || origin_incomplete == Some(true)
        || inspection
            .as_ref()
            .is_none_or(|r| !r.scan_complete || !r.products_determined);
    if incomplete {
        findings.push("This result is incomplete. A displayed singleton is conditional and does not establish uniqueness; retained native-origin assumptions remain applicable.".into());
    }
    let result = Preview {
        engine_version: ENGINE_VERSION.into(),
        context_revision,
        source_revision: source.revision(),
        origin_revision: origin.map(assembly_recipe::hash).transpose()?,
        origin_incomplete,
        origin_findings,
        origin_findings_total,
        origin_findings_truncated,
        request: request.clone(),
        pcr,
        amplicon,
        inspection,
        incomplete,
        findings,
    };
    result.validate_binding(&result.context_revision)?;
    Ok(result)
}
impl Preview {
    /// Readonly bounded documentary validation. Fresh inspection always reruns science.
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        shape(&self.request)?;
        if self.engine_version != ENGINE_VERSION
            || !hex(&self.source_revision)
            || self.origin_revision.as_ref().is_some_and(|s| !hex(s))
            || self.origin_revision.is_some() != self.origin_incomplete.is_some()
            || self.context_revision != expected
            || self.context_revision
                != context(&self.source_revision, &self.origin_revision, &self.request)?
        {
            return Err(bad("stale source/origin/policy/catalog context"));
        }
        if self.origin_findings.len() > 128
            || self
                .origin_findings
                .iter()
                .any(|f| f.chars().count() > 2048)
            || self.origin_findings_total > 32 * 1024 * 1024
            || self.origin_findings.len() > self.origin_findings_total
            || (self.origin_findings_total > self.origin_findings.len()
                && !self.origin_findings_truncated)
            || (self.origin_revision.is_none()
                && (!self.origin_findings.is_empty()
                    || self.origin_findings_total != 0
                    || self.origin_findings_truncated))
        {
            return Err(bad("invalid inherited native finding excerpt"));
        }
        let r = &self.pcr.report;
        let limit = self.request.pcr_options.max_products;
        if r.products.len() > limit
            || self.pcr.documents.len() != r.products.len()
            || r.bindings.iter().any(|b| {
                b.hits.len() > self.request.pcr_options.binding.max_hits || b.warnings.len() > 10000
            })
            || r.warnings.len() > 10000
            || self.pcr.findings.len() > 10000
            || self.findings.len() > 10000
        {
            return Err(Error::TooLarge);
        }
        if r.source_revision != self.source_revision
            || !same(&r.options, &self.request.pcr_options)?
        {
            return Err(bad("PCR source or policy differs"));
        }
        let complete = r.search_complete
            && r.unsupported_geometries == 0
            && r.conflicting_overlaps == 0
            && r.ambiguous_bindings == 0
            && r.ambiguous_products == 0;
        if r.enumeration_complete != complete
            || (r.search_complete && r.bindings.iter().any(|binding| binding.truncated))
            || r.pair_checks > self.request.pcr_options.max_pair_checks
        {
            return Err(bad("PCR completeness evidence is inconsistent"));
        }
        let primers = [&self.request.base.forward, &self.request.base.reverse];
        for (binding, primer) in r.bindings.iter().zip(primers) {
            if binding.source_revision != self.source_revision
                || binding.primer_name != primer.name
                || binding.annealing != primer.annealing.as_str()
                || !same(&binding.options, &self.request.pcr_options.binding)?
            {
                return Err(bad("primer binding evidence differs"));
            }
            for hit in &binding.hits {
                if hit.length > 1000
                    || hit.mismatches.len() > 1000
                    || hit.ambiguous_positions.len() > 1000
                {
                    return Err(Error::TooLarge);
                }
            }
        }
        let models = self
            .pcr
            .documents
            .iter()
            .chain(self.amplicon.iter())
            .collect::<Vec<_>>();
        collagene_model::bounds::model_preflight(&models)?;
        let mut bases = 0usize;
        for (geometry, document) in r.products.iter().zip(&self.pcr.documents) {
            if geometry.forward_primer_index > 1 || geometry.reverse_primer_index > 1 {
                return Err(bad("invalid native oligo-role index"));
            }
            let certain = geometry
                .sequence
                .as_str()
                .bytes()
                .all(|base| b"ACGT".contains(&base))
                && geometry.forward_binding.ambiguous_positions.is_empty()
                && geometry.reverse_binding.ambiguous_positions.is_empty();
            if geometry.certain != certain
                || geometry.overlapping_domains
                    != (geometry.template_length
                        < geometry
                            .forward_binding
                            .length
                            .saturating_add(geometry.reverse_binding.length))
                || geometry.forward_tail_length != primers[geometry.forward_primer_index].tail.len()
                || geometry.reverse_tail_length != primers[geometry.reverse_primer_index].tail.len()
            {
                return Err(bad("native certainty/overlap/domain evidence differs"));
            }
            bases = bases
                .checked_add(geometry.sequence.len())
                .filter(|n| *n <= self.request.pcr_options.max_product_bases)
                .ok_or(Error::TooLarge)?;
            document.validate()?;
            if geometry.substitutions.len() > geometry.sequence.len()
                || geometry.forward_binding.mismatches.len() > 1000
                || geometry.reverse_binding.mismatches.len() > 1000
                || geometry.forward_binding.ambiguous_positions.len() > 1000
                || geometry.reverse_binding.ambiguous_positions.len() > 1000
            {
                return Err(Error::TooLarge);
            }
            if document.topology != Topology::Linear
                || document.sequence != geometry.sequence
                || document.metadata.get("pcr:template_revision") != Some(&self.source_revision)
                || document
                    .metadata
                    .get("pcr:engine_version")
                    .map(String::as_str)
                    != Some(verification::ENGINE_VERSION)
                || document.metadata.get("pcr:mismatch_options")
                    != Some(
                        &serde_json::to_string(&self.request.pcr_options)
                            .map_err(|_| Error::TooLarge)?,
                    )
                || document.metadata.get("pcr:primers")
                    != Some(
                        &serde_json::to_string(&(
                            &self.request.base.forward,
                            &self.request.base.reverse,
                        ))
                        .map_err(|_| Error::TooLarge)?,
                    )
            {
                return Err(bad(
                    "projected PCR document differs from native geometry/request",
                ));
            }
        }
        match (&self.amplicon, &self.inspection) {
            (Some(m), Some(report)) if selectable(&self.pcr) => {
                if !same(m, &self.pcr.documents[0])? {
                    return Err(bad("selected amplicon differs"));
                }
                let catalog = catalog_digest::builtin_catalog()?;
                let enzymes = self.request.base.selection.resolved_enzymes(&catalog)?;
                report.validate_bindings_for(
                    &self.request.base.source_role,
                    m,
                    &enzymes,
                    &self.request.base.options,
                )?;
            }
            (None, None) if !selectable(&self.pcr) => {}
            _ => {
                return Err(bad(
                    "only a sole certain nonoverlapping geometry can expose a digest",
                ));
            }
        }
        let incomplete = !r.enumeration_complete
            || self.origin_incomplete == Some(true)
            || self
                .inspection
                .as_ref()
                .is_none_or(|r| !r.scan_complete || !r.products_determined);
        if self.incomplete != incomplete {
            return Err(bad("incomplete evidence cannot be cleared"));
        }
        if self
            .findings
            .iter()
            .chain(&self.pcr.findings)
            .chain(&r.warnings)
            .chain(r.bindings.iter().flat_map(|b| &b.warnings))
            .any(|s| s.len() > 65536)
        {
            return Err(Error::TooLarge);
        }
        assembly_recipe::hash(self)?;
        Ok(())
    }
}
