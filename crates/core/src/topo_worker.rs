//! Bounded worker bindings for fresh TOPO science. Revision helpers hash borrowed
//! inputs only; readonly decoded reports cannot restore graph authority.
use crate::{Alphabet, Error, Molecule, Topology, topo};
use serde::Serialize;
use std::{collections::BTreeSet, io::Write};
use topo::physical::{self, projection};

const MAX_BYTES: usize = 32 * 1024 * 1024;
fn bad(message: &str) -> Error {
    Error::Unsupported(format!("TOPO worker: {message}"))
}
fn hex(value: &str) -> bool {
    value.len() == 64
        && value
            .bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}
struct JsonHash {
    hash: blake3::Hasher,
    bytes: usize,
    max: usize,
}
impl Write for JsonHash {
    fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
        self.bytes = self
            .bytes
            .checked_add(bytes.len())
            .filter(|&n| n <= self.max)
            .ok_or_else(|| std::io::Error::other("bounded TOPO worker JSON"))?;
        self.hash.update(bytes);
        Ok(bytes.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}
fn hash(value: &impl Serialize, max: usize) -> Result<String, Error> {
    let mut writer = JsonHash {
        hash: blake3::Hasher::new(),
        bytes: 0,
        max,
    };
    serde_json::to_writer(&mut writer, value).map_err(|_| Error::TooLarge)?;
    Ok(writer.hash.finalize().to_hex().to_string())
}
fn identities(inputs: &[topo::Input]) -> Result<Vec<topo::SourceIdentity>, Error> {
    inputs
        .iter()
        .map(|input| {
            Ok(topo::SourceIdentity {
                role: input.role.clone(),
                revision: hash(&input.molecule, MAX_BYTES)?,
            })
        })
        .collect()
}
/// Source adapter's borrowed count/tree/interval/policy guard; no preparation,
/// primer-site scanning, physical graph reconstruction or target comparison.
pub fn source_preflight(inputs: &[topo::Input], request: &topo::Request) -> Result<(), Error> {
    topo::source_preflight(inputs, request)
}
pub fn target_preflight(target: &Molecule) -> Result<(), Error> {
    collagene_model::bounds::model_preflight(&[target])?;
    if target.sequence.alphabet() != Alphabet::Dna || target.sequence.len() > 100_000 {
        return Err(Error::TooLarge);
    }
    target.validate()
}
pub fn simulation_revision(
    inputs: &[topo::Input],
    request: &topo::Request,
) -> Result<String, Error> {
    source_preflight(inputs, request)?;
    hash(
        &(
            topo::ENGINE_VERSION,
            physical::ENGINE,
            physical::PRIMARY_MODEL,
            projection::ENGINE,
            identities(inputs)?,
            request,
        ),
        MAX_BYTES,
    )
}
pub fn orientations_revision(
    inputs: &[topo::Input],
    request: &topo::Request,
    vector: [usize; 2],
    insert: [usize; 2],
    directional: bool,
) -> Result<String, Error> {
    source_preflight(inputs, request)?;
    if !request.transfers.is_empty() {
        return Err(bad("orientation inspection requires no initial transfers"));
    }
    let strands = request
        .recipe
        .materials
        .iter()
        .map(|m| match m {
            topo::MaterialRecipe::Prepared { strands, .. } => strands.len(),
            topo::MaterialRecipe::ExactPcr { .. } => 2,
        })
        .sum::<usize>();
    if vector.iter().chain(&insert).any(|&i| i >= strands) {
        return Err(bad("orientation strand outside preparation"));
    }
    hash(
        &(
            topo::ENGINE_VERSION,
            "two-declared-orientations",
            identities(inputs)?,
            request,
            vector,
            insert,
            directional,
        ),
        MAX_BYTES,
    )
}
pub fn projection_revision(
    inputs: &[topo::Input],
    request: &topo::Request,
    selection: &topo::Selection,
) -> Result<String, Error> {
    let batch = simulation_revision(inputs, request)?;
    if !hex(&selection.batch_revision)
        || selection.batch_revision != batch
        || selection.component >= 64
    {
        return Err(bad(
            "selection does not bind these source inputs and request",
        ));
    }
    hash(
        &(
            topo::ENGINE_VERSION,
            "selected-physical-component",
            batch,
            selection,
        ),
        MAX_BYTES,
    )
}
pub fn target_revision(
    inputs: &[topo::Input],
    request: &topo::Request,
    selection: &topo::Selection,
    target: &Molecule,
    policy: &topo::AcceptancePolicy,
) -> Result<String, Error> {
    target_preflight(target)?;
    let product = projection_revision(inputs, request, selection)?;
    hash(
        &(
            topo::ENGINE_VERSION,
            "independent-target",
            product,
            hash(target, MAX_BYTES)?,
            policy,
        ),
        MAX_BYTES,
    )
}
fn binding(actual: &str, expected: &str) -> Result<(), Error> {
    if !hex(actual) || actual != expected {
        Err(bad("source/options/target response binding changed"))
    } else {
        Ok(())
    }
}
fn source_identities(sources: &[topo::SourceIdentity]) -> Result<(), Error> {
    if sources.is_empty() || sources.len() > 32 {
        return Err(Error::TooLarge);
    }
    let mut roles = BTreeSet::new();
    for source in sources {
        if source.role.trim().is_empty()
            || source.role.len() > 256
            || !hex(&source.revision)
            || !roles.insert(source.role.as_str())
        {
            return Err(bad("invalid source identity description"));
        }
    }
    Ok(())
}
fn preparation(preparation: &topo::PreparationEvidence) -> Result<(), Error> {
    let strands = &preparation.substrates.strands;
    if strands.is_empty()
        || strands.len() > 32
        || preparation.substrates.pairing.len() > 1024
        || preparation.maps.len() > 4096
        || preparation.not_represented.len() > 4096
        || preparation.assumptions.len() > 32
        || preparation.assumptions.iter().any(|s| s.len() > 4096)
    {
        return Err(Error::TooLarge);
    }
    let mut total = 0usize;
    for strand in strands {
        total = total
            .checked_add(strand.sequence.len())
            .filter(|&n| n <= 100_000)
            .ok_or(Error::TooLarge)?;
        if strand.sequence.alphabet() != Alphabet::Dna
            || strand.role.len() > 256
            || strand.id.len() > 256
        {
            return Err(bad("invalid prepared strand description"));
        }
    }
    for pair in &preparation.substrates.pairing {
        if pair.strand_a >= strands.len()
            || pair.strand_b >= strands.len()
            || pair.length == 0
            || pair
                .start_a
                .checked_add(pair.length)
                .is_none_or(|end| end > strands[pair.strand_a].sequence.len())
            || pair
                .start_b
                .checked_add(pair.length)
                .is_none_or(|end| end > strands[pair.strand_b].sequence.len())
        {
            return Err(bad("pairing description outside prepared strands"));
        }
    }
    for map in &preparation.maps {
        if map.strand >= strands.len()
            || map.length == 0
            || map
                .strand_start
                .checked_add(map.length)
                .is_none_or(|end| end > strands[map.strand].sequence.len())
        {
            return Err(bad("source map outside prepared strand"));
        }
    }
    Ok(())
}
fn physical_shape(
    report: &physical::Report,
    preparation: &topo::PreparationEvidence,
) -> Result<(), Error> {
    let graph = &report.product;
    if report.engine != physical::ENGINE
        || report.model != physical::PRIMARY_MODEL
        || !hex(&report.batch_revision)
        || report.source_revisions.len() != preparation.substrates.strands.len()
        || report.source_revisions.iter().any(|s| !hex(s))
        || report.steps.len() > 16
        || report.source_bases > 100_000
        || graph.chains.len() > 64
        || graph.components.len() > 64
        || graph.nicks.len() > 100_000
        || graph.unpaired.len() > 100_000
        || graph.source_pairing.len() > 100_000
        || graph.released_activations.len() > 32
        || report.assumptions.len() > 32
    {
        return Err(bad("invalid physical report bounds/engine"));
    }
    let source_bases = preparation
        .substrates
        .strands
        .iter()
        .map(|s| s.sequence.len())
        .sum::<usize>();
    let mut product_bases = 0usize;
    let mut maps = 0usize;
    for chain in &graph.chains {
        product_bases = product_bases
            .checked_add(chain.sequence.len())
            .filter(|&n| n <= 100_000)
            .ok_or(Error::TooLarge)?;
        maps = maps
            .checked_add(chain.sources.len())
            .filter(|&n| n <= 100_000)
            .ok_or(Error::TooLarge)?;
    }
    if report.source_bases != source_bases || product_bases != source_bases {
        return Err(bad("physical report nucleotide counts disagree"));
    }
    for component in &graph.components {
        if component.chains.len() > 64 || component.chains.iter().any(|&i| i >= graph.chains.len())
        {
            return Err(bad("component chain index outside graph"));
        }
    }
    for run in &graph.unpaired {
        if run.chain >= graph.chains.len()
            || run.start >= run.end
            || run.end > graph.chains[run.chain].sequence.len()
            || run.displaced_by.len() > 16
        {
            return Err(bad("unpaired description outside graph"));
        }
    }
    Ok(())
}
/// Cheap readonly bounds/binding validation only: never replays enzyme activity
/// or accepts the report as input to a fresh scientific operation.
pub fn validate_simulation(
    report: &topo::TopoReport,
    expected_revision: &str,
) -> Result<(), Error> {
    binding(&report.batch_revision, expected_revision)?;
    source_identities(&report.sources)?;
    preparation(&report.preparation)?;
    if report.engine_version != topo::ENGINE_VERSION {
        return Err(bad("adapter engine mismatch"));
    }
    physical_shape(&report.physical, &report.preparation)?;
    hash(report, MAX_BYTES)?;
    Ok(())
}
pub fn validate_orientations(
    report: &topo::OrientationInspection,
    expected_revision: &str,
) -> Result<(), Error> {
    binding(&report.batch_revision, expected_revision)?;
    source_identities(&report.sources)?;
    preparation(&report.preparation)?;
    if report.engine_version != topo::ENGINE_VERSION
        || report.reference.engine != physical::ENGINE
        || report.reference.candidates.len() > 2
        || report.candidate_requests.len() != report.reference.candidates.len()
        || report.reference.biochemical_mixture_exhaustive
    {
        return Err(bad("invalid orientation report description"));
    }
    for candidate in &report.reference.candidates {
        physical_shape(candidate, &report.preparation)?;
    }
    for request in &report.candidate_requests {
        if request.recipe.materials.len() > 16 || request.transfers.len() > 16 {
            return Err(Error::TooLarge);
        }
        hash(request, 1024 * 1024)?;
    }
    hash(report, MAX_BYTES)?;
    Ok(())
}
pub fn validate_projection(product: &topo::Product, expected_revision: &str) -> Result<(), Error> {
    binding(&product.batch_revision, expected_revision)?;
    source_identities(&product.sources)?;
    preparation(&product.preparation)?;
    let proof = &product.projection;
    if product.engine_version != topo::ENGINE_VERSION
        || proof.engine != projection::ENGINE
        || !hex(&proof.batch_revision)
        || !hex(&product.selection.batch_revision)
        || product.selection.component >= proof.physical.product.components.len()
        || proof.selection.component != product.selection.component
        || proof.selection.reference_seed != product.selection.reference_seed
        || proof.selection.source_batch_revision != proof.physical.batch_revision
    {
        return Err(bad("invalid selected contour binding"));
    }
    binding(
        &product.batch_revision,
        &hash(
            &(
                topo::ENGINE_VERSION,
                "selected-physical-component",
                &product.selection.batch_revision,
                &product.selection,
            ),
            MAX_BYTES,
        )?,
    )?;
    physical_shape(&proof.physical, &product.preparation)?;
    match (&proof.outcome, &product.reference) {
        (projection::Outcome::Proven(view), Some(reference)) => {
            collagene_model::bounds::model_preflight(&[reference])?;
            if reference.sequence.alphabet() != Alphabet::Dna
                || reference.sequence.len() > 100_000
                || reference.sequence != view.top
                || view.top.len() != view.bottom.len()
                || view.source_maps.len() > 16384
                || view.nicks.len() > 100_000
                || reference.topology
                    != match view.topology {
                        projection::ContourTopology::Linear => Topology::Linear,
                        projection::ContourTopology::Circular => Topology::Circular,
                    }
                || reference
                    .metadata
                    .get("source:legacy:chemistry_unmapped")
                    .map(String::as_str)
                    != Some("true")
                || reference
                    .metadata
                    .get("topo:physical_graph_required")
                    .map(String::as_str)
                    != Some("true")
            {
                return Err(bad(
                    "reference view disagrees with restrictive contour description",
                ));
            }
            reference.validate()?;
        }
        (projection::Outcome::Withheld(_), None) => {}
        _ => {
            return Err(bad(
                "withheld/proven contour and reference presence disagree",
            ));
        }
    }
    hash(product, MAX_BYTES)?;
    Ok(())
}
pub fn validate_target(
    report: &topo::TargetAssessment,
    expected_revision: &str,
) -> Result<(), Error> {
    binding(&report.batch_revision, expected_revision)?;
    if report.engine_version != topo::ENGINE_VERSION || !hex(&report.target_revision) {
        return Err(bad("invalid target description"));
    }
    validate_projection(&report.product, &report.product.batch_revision)?;
    binding(
        &report.batch_revision,
        &hash(
            &(
                topo::ENGINE_VERSION,
                "independent-target",
                &report.product.batch_revision,
                &report.target_revision,
                &report.policy,
            ),
            MAX_BYTES,
        )?,
    )?;
    match &report.outcome {
        topo::TargetOutcome::AcceptedUnderDeclaredGraph => {
            let projection::Outcome::Proven(view) = &report.product.projection.outcome else {
                return Err(bad("accepted target has no proven contour"));
            };
            if !report.policy.allow_nicks && !view.nicks.is_empty()
                || !report.policy.allow_unknown_terminal_chemistry
                    && !view.fully_known_terminal_chemistry
                || !report.policy.allow_unapplied_transfers && !view.all_requested_transfers_applied
                || !report.policy.allow_other_components
                    && report.product.projection.physical.product.components.len() != 1
            {
                return Err(bad("accepted target conflicts with its explicit policy"));
            }
        }
        topo::TargetOutcome::Rejected { reasons } => {
            if reasons.is_empty() || reasons.len() > 16 || reasons.iter().any(|s| s.len() > 4096) {
                return Err(Error::TooLarge);
            }
        }
        topo::TargetOutcome::Withheld { reason } => {
            if reason.len() > 4096
                || !matches!(
                    report.product.projection.outcome,
                    projection::Outcome::Withheld(_)
                )
            {
                return Err(bad("withheld target description disagrees with contour"));
            }
        }
    }
    hash(report, MAX_BYTES)?;
    Ok(())
}

/// Annotated revisions bind the existing physical request plus the independent
/// coordinate-projection engine. These helpers do not reconstruct either layer.
pub fn annotated_projection_revision(
    inputs: &[topo::Input],
    request: &topo::Request,
    selection: &topo::Selection,
) -> Result<String, Error> {
    topo::annotated_product_revision(inputs, request, selection)
}
pub fn annotated_target_revision(
    inputs: &[topo::Input],
    request: &topo::Request,
    selection: &topo::Selection,
    target: &Molecule,
    policy: &topo::AcceptancePolicy,
) -> Result<String, Error> {
    topo::annotated_target_revision(inputs, request, selection, target, policy)
}
fn annotation_shape(
    evidence: &crate::topo_annotations::AnnotationEvidence,
    physical: &topo::Product,
    length: usize,
) -> Result<(), Error> {
    if evidence.sources.len() != physical.sources.len()
        || evidence.sources.len() > 32
        || evidence.features.len() > 4096
        || evidence.warnings.len() > 65
        || evidence.warnings.iter().any(|s| s.len() > 4096)
    {
        return Err(Error::TooLarge);
    }
    let mut feature_budget = collagene_model::bounds::ModelBudget::default();
    feature_budget.features(&evidence.features)?;
    let (mut maps, mut projected, mut retained) = (0usize, 0usize, 0usize);
    for (source, identity) in evidence.sources.iter().zip(&physical.sources) {
        if source.role != identity.role
            || source.source_revision != identity.revision
            || source.maps.len() > 4096
            || source.unmapped.len() > 4096
            || source.projected > 4096
            || source.dropped > 4096
            || source.warnings.len() > 65
            || source.warnings.iter().any(|s| s.len() > 4096)
            || !source.projection_supported && source.projected != 0
        {
            return Err(bad("annotation role/status/count binding changed"));
        }
        maps = maps
            .checked_add(source.maps.len())
            .filter(|n| *n <= 4096)
            .ok_or(Error::TooLarge)?;
        projected = projected
            .checked_add(source.projected)
            .filter(|n| *n <= 4096)
            .ok_or(Error::TooLarge)?;
        retained = retained
            .checked_add(source.unmapped.len())
            .filter(|n| *n <= 4096)
            .ok_or(Error::TooLarge)?;
        for map in &source.maps {
            if map.source_start >= map.source_end
                || map.source_end > 100_000
                || map
                    .destination_start
                    .checked_add(map.source_end - map.source_start)
                    .is_none_or(|end| end > length)
            {
                return Err(bad("invalid annotation correspondence map"));
            }
        }
        for feature in &source.unmapped {
            if feature.source_index >= 4096
                || feature.reason.len() > 4096
                || feature.reason.trim().is_empty()
                || !feature.feature.needs_review
            {
                return Err(bad("unmapped source annotation lost review evidence"));
            }
            feature_budget.features(std::slice::from_ref(&feature.feature))?;
        }
    }
    if projected != evidence.features.len() {
        return Err(bad("projected feature count differs from evidence"));
    }
    hash(&evidence.features, 4 * 1024 * 1024)?;
    hash(evidence, 16 * 1024 * 1024)?;
    Ok(())
}
/// Readonly shape/source/recipe consistency. A decoded wrapper never authorizes
/// annotation projection or restores terminal/duplex/modification chemistry.
pub fn validate_annotated_projection(
    product: &topo::AnnotatedProduct,
    expected_revision: &str,
) -> Result<(), Error> {
    binding(&product.batch_revision, expected_revision)?;
    if product.engine_version != topo::ANNOTATION_ENGINE_VERSION {
        return Err(bad("annotation engine mismatch"));
    }
    // Guard every native feature tree before recursive validation/hashing.
    if let Some(reference) = &product.reference {
        collagene_model::bounds::model_preflight(&[reference])?;
        if reference.sequence.len() > 100_000 {
            return Err(Error::TooLarge);
        }
    }
    if let Some(evidence) = &product.annotations {
        annotation_shape(
            evidence,
            &product.physical,
            product.reference.as_ref().map_or(0, |r| r.sequence.len()),
        )?;
    }
    validate_projection(&product.physical, &product.physical.batch_revision)?;
    binding(
        &product.batch_revision,
        &hash(
            &(
                topo::ANNOTATION_ENGINE_VERSION,
                "annotated-product",
                &product.physical.batch_revision,
            ),
            MAX_BYTES,
        )?,
    )?;
    match (
        &product.physical.reference,
        &product.annotations,
        &product.reference,
    ) {
        (Some(bare), Some(evidence), Some(reference)) => {
            if reference.name != bare.name
                || reference.sequence != bare.sequence
                || reference.topology != bare.topology
                || reference.features != evidence.features
                || reference.original_genbank.is_some()
                || reference
                    .metadata
                    .get("topo:annotation_engine")
                    .map(String::as_str)
                    != Some(topo::ANNOTATION_ENGINE_VERSION)
                || reference.metadata.get("topo:annotation_revision")
                    != Some(&product.batch_revision)
            {
                return Err(bad(
                    "annotated reference disagrees with its physical view/evidence",
                ));
            }
            for (key, value) in &bare.metadata {
                if reference.metadata.get(key) != Some(value) {
                    return Err(bad(
                        "annotated reference removed physical restriction metadata",
                    ));
                }
            }
            for (key, value) in &reference.metadata {
                if !matches!(
                    key.as_str(),
                    "topo:annotation_engine" | "topo:annotation_revision" | "topo:annotation_scope"
                ) && bare.metadata.get(key) != Some(value)
                {
                    return Err(bad(
                        "annotation metadata cannot add new inherited authority",
                    ));
                }
            }
            reference.validate()?;
        }
        (None, None, None) => {}
        _ => {
            return Err(bad(
                "annotated evidence/reference presence disagrees with contour",
            ));
        }
    }
    hash(product, MAX_BYTES)?;
    Ok(())
}
pub fn validate_annotated_target(
    report: &topo::AnnotatedAssessment,
    expected_revision: &str,
) -> Result<(), Error> {
    binding(&report.batch_revision, expected_revision)?;
    if report.engine_version != topo::ANNOTATION_ENGINE_VERSION || !hex(&report.target_revision) {
        return Err(bad("invalid annotated target engine/revision"));
    }
    validate_annotated_projection(&report.product, &report.product.batch_revision)?;
    binding(
        &report.batch_revision,
        &hash(
            &(
                topo::ANNOTATION_ENGINE_VERSION,
                "annotated-target",
                &report.product.batch_revision,
                &report.target_revision,
                &report.policy,
            ),
            MAX_BYTES,
        )?,
    )?;
    match &report.outcome {
        topo::TargetOutcome::AcceptedUnderDeclaredGraph => {
            let projection::Outcome::Proven(view) = &report.product.physical.projection.outcome
            else {
                return Err(bad("accepted annotated target has no paired contour"));
            };
            if !report.policy.allow_nicks && !view.nicks.is_empty()
                || !report.policy.allow_unknown_terminal_chemistry
                    && !view.fully_known_terminal_chemistry
                || !report.policy.allow_unapplied_transfers && !view.all_requested_transfers_applied
                || !report.policy.allow_other_components
                    && report
                        .product
                        .physical
                        .projection
                        .physical
                        .product
                        .components
                        .len()
                        != 1
            {
                return Err(bad(
                    "accepted annotated target conflicts with explicit physical policy",
                ));
            }
        }
        topo::TargetOutcome::Rejected { reasons } => {
            if reasons.is_empty() || reasons.len() > 16 || reasons.iter().any(|s| s.len() > 4096) {
                return Err(Error::TooLarge);
            }
        }
        topo::TargetOutcome::Withheld { reason } => {
            if reason.len() > 4096
                || !matches!(
                    &report.product.physical.projection.outcome,
                    projection::Outcome::Withheld(_)
                )
            {
                return Err(bad(
                    "withheld annotated target disagrees with physical contour",
                ));
            }
        }
    }
    hash(report, MAX_BYTES)?;
    Ok(())
}
