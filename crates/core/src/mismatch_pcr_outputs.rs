//! Explicit selected outputs from freshly reconstructed conditional mismatch PCR.
//! The host separately binds immutable source events and its current editable draft.
use crate::{
    Error, Molecule, Topology,
    analysis::Primer,
    mutagenesis::{MutantPcrProduct, PcrOptions, PcrReport},
    oligo_formats::OligoTable,
    pcr_outputs::Options,
    verification,
};
use serde::{Deserialize, Serialize};
pub const ENGINE_VERSION: &str = "mismatch-pcr-outputs-1";
macro_rules! dto {($($t:tt)*)=>{#[derive(Clone,Debug,Serialize,Deserialize)]#[cfg_attr(feature="json-schema",derive(schemars::JsonSchema))]#[serde(deny_unknown_fields)]$($t)*};}
dto! {pub struct Request {
    pub source:Molecule,pub forward:Primer,pub reverse:Primer,pub options:PcrOptions,
    pub selected_index:usize,pub expected_product_revision:String,pub outputs:Options,
}}
dto! {pub struct SearchSummary {
    pub source_revision:String,pub source_length:usize,pub options:PcrOptions,
    pub products_reported:usize,pub pair_checks:usize,pub search_complete:bool,
    pub enumeration_complete:bool,pub unsupported_geometries:usize,
    pub conflicting_overlaps:usize,pub ambiguous_bindings:usize,pub ambiguous_products:usize,
    pub rejected_geometries:usize,pub warnings:Vec<String>,
}}
dto! {pub struct Prepared {
    pub context_revision:String,pub selected_index:usize,
    #[cfg_attr(feature="json-schema",schemars(with="serde_json::Value"))]
    pub geometry:MutantPcrProduct,
    pub product:Option<Molecule>,
    #[cfg_attr(feature="json-schema",schemars(with="serde_json::Value"))]
    pub primers:Option<OligoTable>,
    pub search:SearchSummary,pub findings:Vec<String>,
}}
fn bad(message: &str) -> Error {
    Error::Format(format!("Mismatch-PCR outputs: {message}"))
}
fn digest(value: &str) -> bool {
    value.len() == 64
        && value
            .bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}
fn options(value: &PcrOptions) -> Result<(), Error> {
    if value.max_pair_checks == 0
        || value.max_pair_checks > 1_000_000
        || value.max_products == 0
        || value.max_products > 10_000
        || value.max_product_bases == 0
        || value.max_product_bases > 32_000_000
        || value.binding.max_hits == 0
        || value.binding.max_hits > 10_000
        || value.binding.max_work == 0
        || value.binding.max_work > 200_000_000
    {
        return Err(bad("invalid native search/materialization bounds"));
    }
    Ok(())
}
/// Bounded content identity, not authority for a historical event or physical PCR.
pub fn revision(request: &Request) -> Result<String, Error> {
    crate::model::bounds::model_preflight(&[&request.source])?;
    request.source.validate()?;
    request.outputs.validate()?;
    options(&request.options)?;
    if request.selected_index >= request.options.max_products
        || !digest(&request.expected_product_revision)
    {
        return Err(bad("invalid selected geometry or preview revision"));
    }
    for primer in [&request.forward, &request.reverse] {
        if primer.name.trim().is_empty()
            || primer.name.len() > 256
            || primer.name.chars().any(char::is_control)
        {
            return Err(bad(
                "oligo names require 1–256 bytes without control characters",
            ));
        }
        // Scientific domain/tail and mismatch validation remains in the native
        // simulator; this preflight avoids unbounded serialized oligo fields.
        if primer.annealing.len() > 200 || primer.tail.len() > 200 {
            return Err(bad(
                "oligos exceed the existing 200-base mismatch-PCR profile",
            ));
        }
    }
    crate::assembly_recipe::hash(&(ENGINE_VERSION, request))
}
fn summary(report: &PcrReport, source_length: usize) -> SearchSummary {
    SearchSummary {
        source_revision: report.source_revision.clone(),
        source_length,
        options: report.options.clone(),
        products_reported: report.products.len(),
        pair_checks: report.pair_checks,
        search_complete: report.search_complete,
        enumeration_complete: report.enumeration_complete,
        unsupported_geometries: report.unsupported_geometries,
        conflicting_overlaps: report.conflicting_overlaps,
        ambiguous_bindings: report.ambiguous_bindings,
        ambiguous_products: report.ambiguous_products,
        rejected_geometries: report.rejected_geometries,
        warnings: report.warnings.clone(),
    }
}
/// Run in the worker. A prior preview/index alone never authorizes an output.
pub fn prepare(request: &Request) -> Result<Prepared, Error> {
    let context_revision = revision(request)?;
    let preview = verification::simulate_mutagenesis_preview(
        &request.source,
        &request.forward,
        &request.reverse,
        &request.options,
    )?;
    let mut product = preview
        .documents
        .get(request.selected_index)
        .ok_or_else(|| bad("selected geometry is no longer available"))?
        .clone();
    if product.revision() != request.expected_product_revision {
        return Err(bad(
            "selected preview changed; compute and choose its geometry again",
        ));
    }
    let geometry = preview
        .report
        .products
        .get(request.selected_index)
        .ok_or_else(|| bad("selected product has no matching native geometry"))?
        .clone();
    product.name = request.outputs.name.clone();
    let primers = if request.outputs.save_primers {
        Some(crate::oligo_orders::from_primers(&[
            request.forward.clone(),
            request.reverse.clone(),
        ])?)
    } else {
        None
    };
    let mut findings = preview.findings;
    if primers
        .as_ref()
        .is_some_and(|table| table.records.iter().any(|row| row.primer.is_none()))
    {
        findings.push("Some oligos remain order-table evidence only; inspect collection row diagnostics before reuse.".into());
    }
    let prepared = Prepared {
        context_revision,
        selected_index: request.selected_index,
        geometry,
        product: request.outputs.save_product.then_some(product),
        primers,
        search: summary(&preview.report, request.source.sequence.len()),
        findings,
    };
    prepared.validate_binding(&prepared.context_revision)?;
    Ok(prepared)
}
impl Prepared {
    /// Readonly bounds/consistency guard. Fresh native reconstruction and host
    /// source-event/draft checks, rather than serialized receipts, supply authority.
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        options(&self.search.options)?;
        if !digest(expected)
            || self.context_revision != expected
            || !digest(&self.search.source_revision)
            || self.search.source_length == 0
            || self.selected_index >= self.search.products_reported
            || self.search.products_reported > self.search.options.max_products
            || self.search.pair_checks > self.search.options.max_pair_checks
            || (self.product.is_none() && self.primers.is_none())
        {
            return Err(bad(
                "output receipt differs from its selection/search bounds",
            ));
        }
        let complete = self.search.search_complete
            && self.search.unsupported_geometries == 0
            && self.search.conflicting_overlaps == 0
            && self.search.ambiguous_bindings == 0
            && self.search.ambiguous_products == 0;
        if self.search.enumeration_complete != complete {
            return Err(bad("inconsistent search completeness"));
        }
        let g = &self.geometry;
        if g.forward_primer_index > 1
            || g.reverse_primer_index > 1
            || g.forward_start >= self.search.source_length
            || g.reverse_start >= self.search.source_length
            || g.template_length > self.search.source_length
            || g.sequence.len() > self.search.options.max_product_bases
            || g.forward_tail_length
                .saturating_add(g.template_length)
                .saturating_add(g.reverse_tail_length)
                != g.sequence.len()
            || g.substitutions.iter().any(|c| {
                c.source_position >= self.search.source_length
                    || c.product_position >= g.sequence.len()
            })
        {
            return Err(bad("invalid selected native geometry bounds"));
        }
        if let Some(product) = &self.product {
            crate::model::bounds::model_preflight(&[product])?;
            product.validate()?;
            if product.sequence != g.sequence
                || product.topology != Topology::Linear
                || product.metadata.get("pcr:template_revision")
                    != Some(&self.search.source_revision)
            {
                return Err(bad(
                    "reference output does not match selected geometry/source",
                ));
            }
        }
        if let Some(table) = &self.primers {
            table.validate()?;
            if table.records.len() != 2
                || table.records.iter().any(|r| {
                    r.raw_ordered_sequence.is_none()
                        || r.raw_annealing_sequence.is_none()
                        || r.raw_tail.is_none()
                })
            {
                return Err(bad("collection must retain both declared oligos"));
            }
        }
        crate::assembly_recipe::hash(self)?;
        Ok(())
    }
}
