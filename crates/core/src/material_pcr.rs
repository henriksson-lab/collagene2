//! Exact-PCR reference previews from freshly replayed, explicitly selected native
//! material. No target is fabricated and no serialized state grants chemistry.
use crate::{
    Error, Molecule,
    analysis::Primer,
    pcr_outputs,
    run_material::{MaterialOrigin, MaterialSelection},
    verification,
};
use serde::{Deserialize, Serialize};
pub const ENGINE_VERSION: &str = "material-pcr-1";
pub const ORIGIN_METADATA: &str = "pcr:material_origin";
const MAX_FINDINGS: usize = 32;
const MAX_FINDING_CHARS: usize = 512;
#[derive(Clone, Debug, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct OriginEvidence {
    pub engine_version: String,
    pub context_revision: String,
    pub source_revision: String,
    /// Hash of the complete typed native origin, not a chemical receipt.
    pub origin_revision: String,
    pub selection: MaterialSelection,
    pub origin_incomplete: bool,
    pub origin_findings: Vec<String>,
    pub origin_findings_total: usize,
    pub origin_findings_truncated: bool,
}
fn bad(s: &str) -> Error {
    Error::Unsupported(format!("Saved-material PCR: {s}"))
}
fn preflight(
    source: &Molecule,
    origin: &MaterialOrigin,
    forward: &Primer,
    reverse: &Primer,
) -> Result<(), Error> {
    crate::model::bounds::model_preflight(&[source])?;
    crate::run_material_capsule::material_origins_preflight([origin])?;
    source.validate()?;
    if source.revision() != origin.event().revision {
        return Err(bad("source differs from the exact saved material event"));
    }
    for primer in [forward, reverse] {
        if primer.name.len() > 256 || primer.annealing.len() > 1000 || primer.tail.len() > 1000 {
            return Err(Error::TooLarge);
        }
        pcr_outputs::primer(primer)?;
    }
    Ok(())
}
/// Full native source/origin/oligo binding, without evaluating PCR or a target.
pub fn revision(
    source: &Molecule,
    origin: &MaterialOrigin,
    forward: &Primer,
    reverse: &Primer,
) -> Result<String, Error> {
    preflight(source, origin, forward, reverse)?;
    crate::assembly_recipe::hash(&(
        ENGINE_VERSION,
        verification::ENGINE_VERSION,
        source,
        origin,
        forward,
        reverse,
    ))
}
fn evidence(
    source: &Molecule,
    origin: &MaterialOrigin,
    context: String,
    incomplete: bool,
    findings: Vec<String>,
    total: usize,
    truncated: bool,
) -> Result<OriginEvidence, Error> {
    let mut cut = truncated || findings.len() > MAX_FINDINGS;
    let findings = findings
        .into_iter()
        .take(MAX_FINDINGS)
        .map(|f| {
            if f.chars().count() > MAX_FINDING_CHARS {
                cut = true;
                f.chars().take(MAX_FINDING_CHARS).collect()
            } else {
                f
            }
        })
        .collect();
    Ok(OriginEvidence {
        engine_version: ENGINE_VERSION.into(),
        context_revision: context,
        source_revision: source.revision(),
        origin_revision: crate::assembly_recipe::hash(origin)?,
        selection: origin.selection(),
        origin_incomplete: incomplete,
        origin_findings: findings,
        origin_findings_total: total,
        origin_findings_truncated: cut,
    })
}
/// Every call freshly replays the closed original native graph and then asks the
/// exact PCR engine for all supported geometries on its actual private strands.
pub fn simulate(
    source: &Molecule,
    origin: &MaterialOrigin,
    forward: &Primer,
    reverse: &Primer,
) -> Result<Vec<verification::PcrDocument>, Error> {
    let context = revision(source, origin, forward, reverse)?;
    let (mut previews, receipt) = match origin {
        MaterialOrigin::Run(origin_run) => {
            let fresh = crate::run_material_capsule::replay(origin_run, source)?;
            let report = fresh.evidence.decode_report()?;
            let total = report.findings.len();
            let previews = verification::simulate_pcr_documents_from_fresh(
                source,
                forward,
                reverse,
                &fresh.tracked,
            )?;
            (
                previews,
                evidence(
                    source,
                    origin,
                    context,
                    fresh.incomplete,
                    report.findings,
                    total,
                    false,
                )?,
            )
        }
        MaterialOrigin::Digest(origin_digest) => {
            let fresh = crate::run_digest_history::replay(origin_digest, source)?;
            let previews = verification::simulate_pcr_documents_from_fresh(
                source,
                forward,
                reverse,
                fresh.tracked(),
            )?;
            let (_, _, incomplete, findings, total, truncated) = fresh.into_parts();
            (
                previews,
                evidence(
                    source, origin, context, incomplete, findings, total, truncated,
                )?,
            )
        }
    };
    let encoded = serde_json::to_string(&receipt).map_err(|e| bad(&e.to_string()))?;
    let mut bytes = 0usize;
    for preview in &mut previews {
        // This document remains a reference preview, not a replacement stock.
        // The original marker is never removed from the template for replay.
        preview
            .molecule
            .metadata
            .insert(ORIGIN_METADATA.into(), encoded.clone());
        preview.molecule.metadata.insert("duplex_model".into(),"Saved-material PCR reference preview; retain original source event and native preparation recipe. No target verification or independently reusable physical state is asserted.".into());
        if receipt.origin_incomplete {
            preview.findings.push("The native material origin remains incomplete; this PCR reference preview does not resolve its assumptions or search limits.".into());
        }
        preview.findings.extend(
            receipt
                .origin_findings
                .iter()
                .map(|f| format!("Origin finding: {f}")),
        );
        if receipt.origin_findings_truncated {
            preview.findings.push(format!("Origin findings are an excerpt: {} of {} entries; displayed text may also be truncated.",receipt.origin_findings.len(),receipt.origin_findings_total));
        }
        crate::model::bounds::model_preflight(&[&preview.molecule])?;
        bytes = bytes
            .checked_add(
                serde_json::to_vec(preview)
                    .map_err(|e| bad(&e.to_string()))?
                    .len(),
            )
            .filter(|n| *n <= crate::formats::MAX_IMPORT_BYTES)
            .ok_or(Error::TooLarge)?;
    }
    Ok(previews)
}
/// Pins the complete origin in addition to the ordinary selected-output request.
pub fn outputs_revision(
    origin: &MaterialOrigin,
    request: &pcr_outputs::Request,
) -> Result<String, Error> {
    preflight(&request.source, origin, &request.forward, &request.reverse)?;
    pcr_outputs::revision(request)?;
    crate::assembly_recipe::hash(&(
        ENGINE_VERSION,
        "outputs",
        verification::ENGINE_VERSION,
        origin,
        request,
    ))
}
/// Replays again before selection; preview hashes/indices alone never restore
/// private state. The unchanged output DTO is documentary and host-event-bound.
pub fn prepare_outputs(
    origin: &MaterialOrigin,
    request: &pcr_outputs::Request,
) -> Result<pcr_outputs::Prepared, Error> {
    let context = outputs_revision(origin, request)?;
    let previews = simulate(&request.source, origin, &request.forward, &request.reverse)?;
    pcr_outputs::prepare_from_previews(request, context, previews)
}
