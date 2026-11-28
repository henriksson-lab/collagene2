//! One native target-free OE execution path shared by reusable material and digest inspection.
//! Private execution takes only freshly reconstructed source strands; reports never supply state.
use crate::{
    annotations::{self, MappingSegment},
    assembly_recipe, chemistry,
    modified_duplex::TrackedDuplex,
    overlap_extension::{self, Geometry},
    overlap_extension_digest_inspection::PcrStage,
    overlap_extension_pcr,
    pcr_digest_inspection::GeometrySummary,
    verification::{self, PcrDocument},
};
use collagene_model::{Error, Molecule};
use std::collections::BTreeSet;
const MAX_GEOMETRIES: usize = 10_000;
fn bad(s: &str) -> Error {
    Error::Unsupported(format!("overlap-extension execution: {s}"))
}
fn hex(s: &str) -> bool {
    s.len() == 64
        && s.bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}
pub(crate) struct Execution {
    pub stages: Vec<PcrStage>,
    pub extension: Option<Geometry>,
    pub product: Option<Molecule>,
    pub incomplete: bool,
    pub findings: Vec<String>,
}
pub(crate) fn preflight(
    source: &Molecule,
    preparation: &overlap_extension_pcr::Request,
) -> Result<(), Error> {
    collagene_model::bounds::model_preflight(&[source])?;
    source.validate()?;
    if preparation.work(source.sequence.len())?.0 > 100_000_000 {
        return Err(Error::TooLarge);
    }
    Ok(())
}
fn pcr_stage(
    stage: &str,
    template: &Molecule,
    result: Result<Vec<PcrDocument>, Error>,
    findings: &mut Vec<String>,
) -> Result<(PcrStage, Option<PcrDocument>), Error> {
    let mut documents = match result {
        Ok(documents) => documents,
        Err(error) => {
            findings.push(format!("{stage} PCR could not complete: {error}"));
            return Ok((
                PcrStage {
                    stage: stage.into(),
                    template_revision: template.revision(),
                    geometries: vec![],
                    product_revision: None,
                    product_length: None,
                },
                None,
            ));
        }
    };
    if stage == "outer" {
        for document in &mut documents {
            document.molecule.metadata.insert("duplex_model".into(),"OE outer-PCR reference preview; retain the exact source event and native recipe. No standalone stock chemistry or intended-target verification is asserted.".into());
        }
    }
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
    let product = if documents.len() == 1 {
        documents.into_iter().next()
    } else {
        findings.push(format!("{stage} PCR retained{} geometries; downstream material is withheld rather than selecting a sequence matching an intended product.",documents.len()));
        None
    };
    let summary = PcrStage {
        stage: stage.into(),
        template_revision: template.revision(),
        geometries,
        product_revision: product.as_ref().map(|p| p.molecule.revision()),
        product_length: product.as_ref().map(|p| p.molecule.sequence.len()),
    };
    if let Some(p) = &product {
        findings.extend(p.findings.iter().take(100).cloned());
    }
    Ok((summary, product))
}
fn project(
    source: &Molecule,
    offset: usize,
    product: &mut Molecule,
    findings: &mut Vec<String>,
) -> Result<(), Error> {
    if source.features.is_empty() {
        return Ok(());
    }
    let projected = annotations::project_features_with_topology(
        source,
        &[MappingSegment {
            source_start: 0,
            source_end: source.sequence.len(),
            destination_start: offset,
            reverse: false,
        }],
        product.sequence.len(),
        product.topology,
    )?;
    product.features.extend(projected.features);
    if !projected.unmapped_features.is_empty() {
        product.metadata.insert(
            format!("annotation_mapping:unmapped:{}", source.revision()),
            serde_json::to_string(&projected.unmapped_features).map_err(|_| Error::TooLarge)?,
        );
    }
    findings.extend(
        projected
            .warnings
            .iter()
            .take(100)
            .map(|s| s.chars().take(2048).collect()),
    );
    let mut seen = BTreeSet::new();
    let mut features = Vec::with_capacity(product.features.len());
    for feature in product.features.drain(..) {
        if seen.insert(assembly_recipe::hash(&feature)?) {
            features.push(feature);
        }
    }
    product.features = features;
    Ok(())
}

pub(crate) fn execute(
    source: &Molecule,
    tracked: Option<&TrackedDuplex>,
    preparation: &overlap_extension_pcr::Request,
) -> Result<Execution, Error> {
    preflight(source, preparation)?;
    if let Some(tracked) = tracked.filter(|m| m.duplex.nicks.is_empty()) {
        let declaration = assembly_recipe::Preparation::OverlapExtensionPcr {
            request: preparation.clone(),
        };
        for pair in [&preparation.left, &preparation.right] {
            let template =
                crate::pcr_source::Template::fresh(source, tracked, pair.template_strand)?;
            if assembly_recipe::preparation_work_from_frame(template, &declaration)? > 100_000_000 {
                return Err(Error::TooLarge);
            }
        }
    }
    let mut findings=vec!["Fresh ideal exact-PCR branches, selected reciprocal extension and outer PCR. There is no intended target input. Geometry uniqueness is conditional on this model and is not experimental isolation or reaction success.".into()];
    let mut incomplete = false;
    let initial = |pair: &overlap_extension_pcr::PrimerPair| {
        if let Some(tracked) = tracked {
            verification::simulate_pcr_documents_from_fresh(
                source,
                &pair.forward,
                &pair.reverse,
                tracked,
                pair.template_strand,
            )
        } else {
            verification::simulate_pcr_documents(
                source,
                &pair.forward,
                &pair.reverse,
                pair.template_strand,
            )
        }
    };
    let (left_stage, left) = pcr_stage("left", source, initial(&preparation.left), &mut findings)?;
    let (right_stage, right) =
        pcr_stage("right", source, initial(&preparation.right), &mut findings)?;
    let mut pcr_stages = vec![left_stage, right_stage];
    let mut extension = None;
    let mut product_document = None;
    if let (Some(left), Some(right)) = (left, right) {
        let a = TrackedDuplex::exact_pcr(&left.molecule.sequence)?;
        let b = TrackedDuplex::exact_pcr(&right.molecule.sequence)?;
        match a.overlap_extend(&b, &preparation.extension) {
            Ok((extended, geometry)) => {
                incomplete |=
                    !geometry.terminal_pairings.complete || !geometry.terminal_pairings.unique;
                findings.extend(geometry.warnings.iter().take(100).cloned());
                let mut reference = extended
                    .duplex
                    .molecule("Fresh reciprocal overlap-extension substrate")?;
                project(&left.molecule, 0, &mut reference, &mut findings)?;
                project(
                    &right.molecule,
                    a.duplex.top.len() - preparation.extension.overlap,
                    &mut reference,
                    &mut findings,
                )?;
                let outer = &preparation.outer;
                let (outer_stage, product) = pcr_stage(
                    "outer",
                    &reference,
                    verification::simulate_pcr_documents_from_fresh(
                        &reference,
                        &outer.forward,
                        &outer.reverse,
                        &extended,
                        outer.template_strand,
                    ),
                    &mut findings,
                )?;
                pcr_stages.push(outer_stage);
                if let Some(product) = product {
                    product_document = Some(product.molecule);
                } else {
                    incomplete = true;
                }
                extension = Some(geometry);
            }
            Err(error) => {
                incomplete = true;
                findings.push(format!("Reciprocal extension unavailable: {error}"));
            }
        }
    } else {
        incomplete = true;
    }

    Ok(Execution {
        stages: pcr_stages,
        extension,
        product: product_document,
        incomplete,
        findings,
    })
}

/// Documentary shape/correspondence only: no geometry search, source replay or extension.
pub(crate) fn validate_evidence(
    source_revision: &str,
    preparation: &overlap_extension_pcr::Request,
    stages: &[PcrStage],
    extension: &Option<Geometry>,
    product: &Option<Molecule>,
) -> Result<(), Error> {
    preparation.validate_shape()?;
    if !hex(source_revision) || stages.len() < 2 || stages.len() > 3 {
        return Err(bad("invalid stage count/source revision"));
    }
    for (i, p) in stages.iter().enumerate() {
        if p.stage != ["left", "right", "outer"][i]
            || !hex(&p.template_revision)
            || i < 2 && p.template_revision != source_revision
            || p.geometries.len() > MAX_GEOMETRIES
        {
            return Err(bad("invalid native PCR stage binding"));
        }
        for (index, g) in p.geometries.iter().enumerate() {
            let expected_strand = [
                preparation.left.template_strand,
                preparation.right.template_strand,
                preparation.outer.template_strand,
            ][i];
            if g.template_strand != expected_strand
                || g.index != index
                || g.product_length == 0
                || g.product_length > collagene_model::MAX_SEQUENCE_LEN
                || g.forward_start > collagene_model::MAX_SEQUENCE_LEN
                || g.reverse_start > collagene_model::MAX_SEQUENCE_LEN
            {
                return Err(bad("invalid native PCR geometry summary"));
            }
        }
        match (&p.product_revision, p.product_length) {
            (Some(revision), Some(length))
                if p.geometries.len() == 1
                    && hex(revision)
                    && length == p.geometries[0].product_length => {}
            (None, None) if p.geometries.len() != 1 => {}
            _ => {
                return Err(bad(
                    "only a native singleton PCR may expose a selected stage product",
                ));
            }
        }
    }
    if let Some(g) = &extension {
        g.validate_documentary()?;
        if g.engine_version != overlap_extension::ENGINE_VERSION
            || g.options != preparation.extension
            || g.homology.len() != g.options.overlap
            || g.mappings.len() > 8
            || g.fill.mappings.len() > 8
            || g.source_geometry_revisions.iter().any(|s| !hex(s))
            || g.terminal_pairings.overlaps.len() > g.options.max_candidates
            || g.terminal_pairings.comparisons > g.options.max_comparisons
            || g.terminal_pairings.unique
                && (!g.terminal_pairings.complete
                    || g.terminal_pairings.first_alternative.is_some())
            || g.fill.source != g.annealed
            || g.fill.method != chemistry::EndRepairMethod::FillFivePrimeOverhangs
            || stages[0].product_length != Some(g.annealed.top.len())
            || stages[1].product_length != Some(g.annealed.bottom.len())
            || stages.len() != 3
        {
            return Err(bad(
                "extension evidence differs from actual singleton branch lengths/options",
            ));
        }
    } else if stages.len() != 2 || product.is_some() {
        return Err(bad(
            "downstream products require retained native extension evidence",
        ));
    }
    if let Some(m) = &product {
        collagene_model::bounds::model_preflight(&[m])?;
        m.validate()?;
        let outer = &stages[2];
        let amplicon_revision = m.revision();
        if m.topology != collagene_model::Topology::Linear
            || outer.product_length != Some(m.sequence.len())
            || outer.product_revision.as_deref() != Some(amplicon_revision.as_str())
            || m.metadata.get("pcr:template_revision") != Some(&outer.template_revision)
            || m.metadata.get("pcr:engine_version").map(String::as_str)
                != Some(verification::ENGINE_VERSION)
            || m.metadata.get("pcr:primers")
                != Some(
                    &serde_json::to_string(&(
                        &preparation.outer.forward,
                        &preparation.outer.reverse,
                    ))
                    .map_err(|_| Error::TooLarge)?,
                )
        {
            return Err(bad("outer amplicon differs from its native PCR stage"));
        }
    } else if stages.get(2).is_some_and(|p| p.product_revision.is_some()) {
        return Err(bad("digest fragments require a native singleton outer PCR"));
    }
    Ok(())
}
