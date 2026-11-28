//! Exact and intentional-mismatch PCR previews, with one native source authority.
use super::ENGINE_VERSION;
use super::FreshMaterial;
use super::SourceOrigin;
use super::bad;
use super::encoded;
use super::fresh_source;
use super::limitations;
use super::original_source;
use super::preview_receipt_shape;
use super::preview_source_shape;
use super::source_origin_preflight;
use crate::chemistry;
use crate::chemistry::tracked_duplex::ModificationSummary;
use crate::sequence::analysis::Primer;
use crate::verification;
use collagene_model::Error;
use collagene_model::Molecule;
use serde::Deserialize;
use serde::Serialize;

/// Source-bound reference previews only; this helper creates no native product.
pub fn pcr_preview_revision(
    source: &Molecule,
    origin: &SourceOrigin,
    forward: &Primer,
    reverse: &Primer,
    template_strand: chemistry::Strand,
) -> Result<String, Error> {
    collagene_model::bounds::model_preflight(&[source])?;
    source.validate()?;
    for primer in [forward, reverse] {
        if primer.annealing.len() > 1000 || primer.tail.len() > 1000 {
            return Err(Error::TooLarge);
        }
        crate::cloning::pcr::outputs::primer(primer)?;
    }
    match origin {
        SourceOrigin::OriginalReference => {
            if template_strand != chemistry::Strand::Top {
                return Err(bad(
                    "Original reference PCR uses explicit Top; Bottom requires native source replay",
                ));
            }
            original_source(source)?;
        }
        SourceOrigin::Saved(origin) => {
            origin.validate_documentary()?;
            if origin.event().revision != source.revision() {
                return Err(bad("PCR preview source differs from saved origin."));
            }
        }
    }
    crate::encoding::bounded_json::hash(&(
        ENGINE_VERSION,
        "pcr-preview",
        verification::ENGINE_VERSION,
        source,
        origin,
        forward,
        reverse,
        template_strand,
    ))
}
pub fn preview_pcr(
    source: &Molecule,
    origin: &SourceOrigin,
    forward: &Primer,
    reverse: &Primer,
    template_strand: chemistry::Strand,
) -> Result<Vec<verification::pcr::PcrDocument>, Error> {
    pcr_preview_revision(source, origin, forward, reverse, template_strand)?;
    match origin {
        SourceOrigin::OriginalReference => {
            verification::pcr::simulate_pcr_documents(source, forward, reverse, template_strand)
        }
        SourceOrigin::Saved(origin) => crate::cloning::pcr::material::simulate(
            source,
            origin,
            forward,
            reverse,
            template_strand,
        ),
    }
}

// Native source-bound previews share the same origin replay and mutation kernels
// as installation. They are documentary objects, never deserialized material.
dto! {pub struct MismatchPreview{pub engine_version:String,pub context_revision:String,pub source_revision:String,
#[cfg_attr(feature="json-schema",schemars(with="serde_json::Value"))]pub preview:verification::pcr::MutagenesisPreview,
pub source_duplex:chemistry::DuplexSummary,pub source_modifications:ModificationSummary,
pub incomplete:bool,pub findings:Vec<String>,pub findings_total:usize,pub findings_truncated:bool}}
pub fn mismatch_pcr_preview_revision(
    source: &Molecule,
    origin: &SourceOrigin,
    forward: &Primer,
    reverse: &Primer,
    options: &crate::cloning::pcr::mutagenesis::PcrOptions,
) -> Result<String, Error> {
    crate::cloning::assembly::recipe::mismatch_pcr_shape(forward, reverse, options)?;
    source_origin_preflight(source, origin)?;
    if crate::cloning::assembly::recipe::mismatch_pcr_work(
        source.sequence.len(),
        forward,
        reverse,
        options,
    )? > 200_000_000
    {
        return Err(Error::TooLarge);
    }
    crate::encoding::bounded_json::hash(&(
        ENGINE_VERSION,
        "mismatch-pcr-preview",
        verification::ENGINE_VERSION,
        source,
        origin,
        forward,
        reverse,
        options,
    ))
}
pub(super) fn mismatch_from_fresh(
    source: &Molecule,
    origin: &SourceOrigin,
    fresh: &FreshMaterial,
    forward: &Primer,
    reverse: &Primer,
    options: &crate::cloning::pcr::mutagenesis::PcrOptions,
) -> Result<MismatchPreview, Error> {
    let context_revision =
        mismatch_pcr_preview_revision(source, origin, forward, reverse, options)?;
    let report = crate::cloning::pcr::mutagenesis::simulate_pcr_reconstructed(
        source,
        &fresh.tracked.duplex,
        forward,
        reverse,
        options,
    )?;
    let mut preview =
        verification::pcr::project_mutagenesis_preview(source, forward, reverse, options, report)?;
    for document in &mut preview.documents {
        document.metadata.insert("duplex_model".into(),"Intentional-mismatch PCR reference preview only; replay the exact native derivation to establish reusable material.".into());
        document.metadata.insert(
            "source:material_derivation:preview_context".into(),
            context_revision.clone(),
        );
    }
    let local = preview
        .findings
        .iter()
        .chain(&preview.report.warnings)
        .cloned()
        .collect::<Vec<_>>();
    let total = fresh
        .findings_total
        .checked_add(local.len())
        .ok_or(Error::TooLarge)?;
    let (findings, findings_total, findings_truncated) = limitations(
        local
            .into_iter()
            .chain(fresh.findings.iter().cloned())
            .collect(),
        total,
        fresh.findings_truncated,
    );
    let result = MismatchPreview {
        engine_version: ENGINE_VERSION.into(),
        context_revision,
        source_revision: source.revision(),
        source_duplex: fresh.tracked.duplex.summary()?,
        source_modifications: fresh.tracked.summary()?,
        incomplete: fresh.incomplete || !preview.report.enumeration_complete,
        preview,
        findings,
        findings_total,
        findings_truncated,
    };
    result.validate_binding(&result.context_revision)?;
    Ok(result)
}
pub fn preview_mismatch_pcr(
    source: &Molecule,
    origin: &SourceOrigin,
    forward: &Primer,
    reverse: &Primer,
    options: &crate::cloning::pcr::mutagenesis::PcrOptions,
) -> Result<MismatchPreview, Error> {
    mismatch_pcr_preview_revision(source, origin, forward, reverse, options)?;
    let fresh = fresh_source(source, origin)?;
    mismatch_from_fresh(source, origin, &fresh, forward, reverse, options)
}
impl MismatchPreview {
    /// Exact same eligibility as installation. Incomplete enumeration stays conditional;
    /// two reported geometries never establish a selected isolated product.
    pub fn reusable(&self) -> bool {
        self.validate_binding(&self.context_revision).is_ok()
            && self.preview.report.products.len() == 1
            && self.preview.documents.len() == 1
            && self.preview.report.products[0].certain
            && !self.preview.report.products[0].overlapping_domains
    }
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        preview_receipt_shape(
            &self.engine_version,
            &self.context_revision,
            expected,
            &self.source_revision,
            &self.findings,
            self.findings_total,
            self.findings_truncated,
        )?;
        preview_source_shape(&self.source_duplex, &self.source_modifications)?;
        let report = &self.preview.report;
        if self.source_duplex.bottom_start != 0
            || self.source_duplex.top_length != self.source_duplex.bottom_length
            || !self.source_duplex.nicks.is_empty()
        {
            return Err(bad(
                "mismatch preview requires actual fully paired intact source geometry",
            ));
        }
        for binding in &report.bindings {
            if binding.source_revision != self.source_revision
                || crate::encoding::bounded_json::hash(&binding.options)?
                    != crate::encoding::bounded_json::hash(&report.options.binding)?
            {
                return Err(bad(
                    "mismatch preview binding report differs from declared policy/source",
                ));
            }
        }
        if report.source_revision != self.source_revision
            || report.products.len() != self.preview.documents.len()
            || report.products.len() > 10000
            || report.pair_checks > report.options.max_pair_checks
            || (!report.enumeration_complete && !self.incomplete)
        {
            return Err(bad("mismatch preview report correspondence differs"));
        }
        for (document, geometry) in self.preview.documents.iter().zip(&report.products) {
            collagene_model::bounds::model_preflight(&[document])?;
            document.validate()?;
            if document.sequence != geometry.sequence
                || document.topology != collagene_model::Topology::Linear
                || document
                    .metadata
                    .get("source:material_derivation:preview_context")
                    .map(String::as_str)
                    != Some(expected)
            {
                return Err(bad(
                    "mismatch preview geometry/document correspondence differs",
                ));
            }
        }
        encoded(self, 32 * 1024 * 1024)?;
        Ok(())
    }
}
