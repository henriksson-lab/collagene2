//! Paired original-oligo preview and explicit annealing candidate selection.
use super::AnnealSelection;
use super::ENGINE_VERSION;
use super::SourceBinding;
use super::SourceRole;
use super::bad;
use super::declared_oligo_source;
use super::encoded;
use super::hex;
use super::validate_source_bindings;
use crate::chemistry;
use crate::documents::workspace::Workspace;
use crate::material::lineage::EventBinding;
use crate::material::run as run_material;
use crate::verification;
use collagene_model::Error;
use collagene_model::Molecule;
use serde::Deserialize;
use serde::Serialize;
use std::collections::BTreeSet;

pub(super) fn anneal_options_shape(options: &chemistry::AnnealOptions) -> Result<(), Error> {
    crate::cloning::assembly::recipe::anneal_shape("bottom", options, 0)
}
pub(super) fn candidate_id_shape(id: &str) -> bool {
    id.strip_prefix("anneal-").is_some_and(hex)
}
pub(super) fn anneal_context(
    sources: &[SourceBinding],
    options: &chemistry::AnnealOptions,
) -> Result<String, Error> {
    crate::encoding::bounded_json::hash(&(
        ENGINE_VERSION,
        "annealed-oligo-preview",
        verification::ENGINE_VERSION,
        sources,
        options,
    ))
}
pub(super) fn protect_anneal_preview(product: &mut Molecule) -> Result<(), Error> {
    product.metadata.insert("source:material_derivation:scope".into(),"Source-backed standard-oligo annealing reference preview. Both exact native source events and the selected geometry recipe are required for reusable material.".into());
    product.validate()
}
dto! {pub struct AnnealedOligosPreview{
pub engine_version:String,pub context_revision:String,pub sources:Vec<SourceBinding>,pub options:chemistry::AnnealOptions,
#[cfg_attr(feature="json-schema",schemars(with="serde_json::Value"))]pub preview:verification::terminal::AnnealPreview,
pub products:Vec<Molecule>,pub incomplete:bool,pub findings:Vec<String>}}
pub fn annealed_oligos_preview_revision(
    w: &Workspace,
    top: &EventBinding,
    bottom: &EventBinding,
    options: &chemistry::AnnealOptions,
) -> Result<String, Error> {
    anneal_options_shape(options)?;
    run_material::event_shape(top)?;
    run_material::event_shape(bottom)?;
    w.validate()?;
    for event in [top, bottom] {
        declared_oligo_source(run_material::resolve(w, event)?)?;
    }
    anneal_context(
        &[
            SourceBinding {
                role: SourceRole::TopOligo,
                event: top.clone(),
            },
            SourceBinding {
                role: SourceRole::BottomOligo,
                event: bottom.clone(),
            },
        ],
        options,
    )
}
pub fn preview_annealed_oligos(
    w: &Workspace,
    top: &EventBinding,
    bottom: &EventBinding,
    options: &chemistry::AnnealOptions,
) -> Result<AnnealedOligosPreview, Error> {
    let context_revision = annealed_oligos_preview_revision(w, top, bottom, options)?;
    let a = run_material::resolve(w, top)?;
    let b = run_material::resolve(w, bottom)?;
    let preview = verification::terminal::simulate_oligo_annealing(a, b, options)?;
    let mut findings = preview.annealing.warnings.clone();
    let mut products = Vec::new();
    for index in 0..preview.annealing.candidates.len() {
        let mut product = verification::terminal::annealed_document_inner(
            a,
            b,
            &preview.annealing,
            index,
            &mut findings,
        )?
        .molecule;
        protect_anneal_preview(&mut product)?;
        products.push(product);
    }
    let result = AnnealedOligosPreview {
        engine_version: ENGINE_VERSION.into(),
        context_revision,
        sources: vec![
            SourceBinding {
                role: SourceRole::TopOligo,
                event: top.clone(),
            },
            SourceBinding {
                role: SourceRole::BottomOligo,
                event: bottom.clone(),
            },
        ],
        options: options.clone(),
        incomplete: !preview.annealing.complete || !preview.annealing.unique,
        preview,
        products,
        findings,
    };
    result.validate_binding(&result.context_revision)?;
    Ok(result)
}
impl AnnealedOligosPreview {
    /// Documentary binding only: does not rerun annealing or establish private state.
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        validate_source_bindings(&self.sources)?;
        anneal_options_shape(&self.options)?;
        if self.sources.len() != 2
            || self.engine_version != ENGINE_VERSION
            || !hex(expected)
            || self.context_revision != expected
            || self.context_revision != anneal_context(&self.sources, &self.options)?
            || self.preview.top_revision != self.sources[0].event.revision
            || self.preview.bottom_revision != self.sources[1].event.revision
            || self.preview.annealing.options != self.options
            || self.preview.annealing.engine_version != chemistry::ANNEAL_ENGINE_VERSION
            || self.products.len() != self.preview.annealing.candidates.len()
            || self.products.len() > self.options.max_candidates
            || self.preview.annealing.comparisons > self.options.max_comparisons
            || self.preview.annealing.unique
                != (self.preview.annealing.complete && self.products.len() == 1)
            || self.incomplete
                != (!self.preview.annealing.complete || !self.preview.annealing.unique)
            || self.findings.len() > 1000
            || self.findings.iter().any(|s| s.len() > 8192)
        {
            return Err(bad(
                "Annealed preview changed source roles, candidate evidence or context",
            ));
        }
        let mut ids = BTreeSet::new();
        for (candidate, product) in self.preview.annealing.candidates.iter().zip(&self.products) {
            candidate.duplex.validate()?;
            product.validate()?;
            if !candidate_id_shape(&candidate.id)
                || !ids.insert(&candidate.id)
                || product.sequence != candidate.duplex.top
                || product.topology != collagene_model::Topology::Linear
                || product.metadata.get("anneal:top_revision")
                    != Some(&self.sources[0].event.revision)
                || product.metadata.get("anneal:bottom_revision")
                    != Some(&self.sources[1].event.revision)
                || product.metadata.get("anneal:candidate_id") != Some(&candidate.id)
                || product.metadata.get("duplex_model")
                    != Some(
                        &serde_json::to_string(&candidate.duplex.summary()?)
                            .map_err(|_| Error::TooLarge)?,
                    )
                || !product
                    .metadata
                    .contains_key("source:material_derivation:scope")
            {
                return Err(bad(
                    "Annealed product differs from its retained candidate/source geometry",
                ));
            }
        }
        encoded(self, crate::material::lineage::MAX_RUN_BYTES)?;
        Ok(())
    }
    pub fn selection(&self, candidate_id: &str) -> Result<AnnealSelection, Error> {
        self.validate_binding(&self.context_revision)?;
        let index = self
            .preview
            .annealing
            .candidates
            .iter()
            .position(|c| c.id == candidate_id)
            .ok_or_else(|| bad("choose a retained native annealing candidate ID"))?;
        Ok(AnnealSelection {
            expected_preview_revision: self.context_revision.clone(),
            candidate_id: candidate_id.into(),
            expected_product_revision: self.products[index].revision(),
        })
    }
}
