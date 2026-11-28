//! Bounded worker bindings for fresh TOPO science. Revision helpers hash borrowed
//! inputs only; readonly decoded reports cannot restore graph authority.
use crate::cloning::topo;
use crate::cloning::topo::physical::{self, projection};
use collagene_model::{Alphabet, Error, Molecule, Topology};
use serde::Serialize;
use std::{collections::BTreeSet, io::Write};

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
    topo::batch_revision(inputs, request)
}
pub fn orientations_revision(
    inputs: &[topo::Input],
    request: &topo::Request,
    vector: [usize; 2],
    insert: [usize; 2],
    directional: bool,
) -> Result<String, Error> {
    topo::orientations_revision(inputs, request, vector, insert, directional)
}
pub fn projection_revision(
    inputs: &[topo::Input],
    request: &topo::Request,
    selection: &topo::Selection,
) -> Result<String, Error> {
    topo::product_revision(inputs, request, selection)
}
pub fn target_revision(
    inputs: &[topo::Input],
    request: &topo::Request,
    selection: &topo::Selection,
    target: &Molecule,
    policy: &topo::AcceptancePolicy,
) -> Result<String, Error> {
    topo::target_revision(inputs, request, selection, target, policy)
}
fn binding(actual: &str, expected: &str) -> Result<(), Error> {
    if !hex(actual) || actual != expected {
        Err(bad("source/options/target response binding changed"))
    } else {
        Ok(())
    }
}
fn source_identities(sources: &[topo::SourceIdentity]) -> Result<(), Error> {
    if sources.is_empty() || sources.len() > 64 {
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
pub(crate) fn preparation(preparation: &topo::PreparationEvidence) -> Result<(), Error> {
    let strands = &preparation.substrates.strands;
    if strands.is_empty()
        || strands.len() > 32
        || preparation.substrates.pairing.len() > 1024
        || preparation.maps.len() > 4096
        || preparation.not_represented.len() > 4096
        || preparation.assumptions.len() > 128
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
    if preparation.retained_steps.len() > 64
        || preparation.layout.len() > 256
        || preparation.graph_origins.len() > 128
    {
        return Err(Error::TooLarge);
    }
    let mut limits = physical::Limits::default();
    limits.max_steps = 64;
    limits.max_source_bases = 500_000;
    limits.max_strands = 64;
    limits.max_work_cells = 128_000_000;
    limits.max_work_bytes = 128 * 1024 * 1024;
    limits.max_report_bytes = 32 * 1024 * 1024;
    physical::preflight(
        &preparation.substrates,
        &preparation.retained_steps,
        &limits,
    )
    .map_err(|e| bad(&e.to_string()))?;
    for layout in &preparation.layout {
        if layout.material >= 64
            || layout.owner.as_ref().is_some_and(|owner| !hex(owner))
            || layout.source_strands.is_empty()
            || layout.source_strands.len() > 64
            || layout.source_strands.iter().any(|&i| i >= strands.len())
            || layout.source_roles.len() > 64
            || layout.source_roles.iter().any(|(a, b)| {
                a.trim().is_empty() || b.trim().is_empty() || a.len() > 256 || b.len() > 256
            })
        {
            return Err(bad("invalid native material/source layout"));
        }
        if layout.source_strands.iter().collect::<BTreeSet<_>>().len()
            != layout.source_strands.len()
        {
            return Err(bad("native material layout repeats a source strand"));
        }
    }
    for origin in &preparation.graph_origins {
        if !hex(&origin.selection.record_hash)
            || !hex(&origin.selection.state_revision)
            || origin.strand_start >= origin.strand_end
            || origin.strand_end > strands.len()
            || origin.after_steps > preparation.retained_steps.len()
            || origin.step_map.len() > 64
            || origin.step_map.iter().any(|&i| i >= origin.after_steps)
            || origin.activation_unmap.len() > 64
            || origin.activation_unmap.iter().any(|(a, b)| {
                a.trim().is_empty() || b.trim().is_empty() || a.len() > 256 || b.len() > 256
            })
        {
            return Err(bad("invalid retained graph source checkpoint description"));
        }
        if origin.step_map.windows(2).any(|pair| pair[0] >= pair[1])
            || origin
                .activation_unmap
                .iter()
                .map(|(a, _)| a)
                .collect::<BTreeSet<_>>()
                .len()
                != origin.activation_unmap.len()
        {
            return Err(bad("retained graph source mapping is not one-to-one"));
        }
    }
    Ok(())
}
fn physical_shape(
    report: &physical::Report,
    preparation: &topo::PreparationEvidence,
) -> Result<(), Error> {
    physical::validate_report(
        &preparation.substrates,
        &report.request,
        &report.limits,
        report,
    )
    .map_err(|e| bad(&e.to_string()))?;
    if report.request.len() < preparation.retained_steps.len()
        || report.request[..preparation.retained_steps.len()] != preparation.retained_steps
    {
        return Err(bad("native retained source transcript differs"));
    }
    Ok(())
}
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
    if (!report.physical.all_requested_steps_applied
        || report.physical.unsupported_or_unknown
        || report
            .preparation
            .graph_origins
            .iter()
            .any(|o| o.incomplete))
        && !report.incomplete
    {
        return Err(bad("native source or transfer incompleteness was lost"));
    }
    if report.findings.len() > 128 || report.findings.iter().any(|s| s.len() > 4096) {
        return Err(Error::TooLarge);
    }
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
        || report.reference.pending_steps.len() != report.reference.candidates.len()
        || report.reference.satisfied_prefix_steps.len() != report.reference.candidates.len()
        || report.reference.checkpoints.len() != report.preparation.graph_origins.len()
        || report.reference.biochemical_mixture_exhaustive
    {
        return Err(bad("invalid orientation report description"));
    }
    for candidate in &report.reference.candidates {
        physical_shape(candidate, &report.preparation)?;
    }
    for (index, request) in report.candidate_requests.iter().enumerate() {
        if request.recipe.materials.len() > 64 || request.steps.len() > 16 {
            return Err(Error::TooLarge);
        }
        crate::material::capsule::native_origins_preflight(
            request.recipe.materials.iter().filter_map(|m| {
                if let topo::MaterialRecipe::SavedGraph { origin, .. } = m {
                    Some(crate::material::capsule::NativeOriginRef::TopoGraph(origin))
                } else {
                    None
                }
            }),
        )?;
        let candidate = &report.reference.candidates[index];
        let pending = &report.reference.pending_steps[index];
        let prefix = report.preparation.retained_steps.len();
        let satisfied = &report.reference.satisfied_prefix_steps[index];
        if &request.steps != pending
            || pending.len() > 2
            || candidate.steps.len() != prefix + pending.len()
            || candidate.request[prefix..] != pending[..]
            || satisfied.len() > 2
            || satisfied
                .iter()
                .any(|&i| i >= prefix || candidate.steps[i].outcome() != physical::Outcome::Applied)
            || satisfied.iter().collect::<BTreeSet<_>>().len() != satisfied.len()
        {
            return Err(bad(
                "orientation suffix/prefix satisfaction description changed",
            ));
        }
        hash(request, MAX_BYTES)?;
    }
    for (origin, graph) in report
        .preparation
        .graph_origins
        .iter()
        .zip(&report.reference.checkpoints)
    {
        if graph.chains.len() > 128
            || graph.components.len() > 128
            || graph.nicks.len() > 100_000
            || graph.unpaired.len() > 100_000
            || graph.source_pairing.len() > 100_000
            || graph.released_activations.len() > 64
            || graph
                .chains
                .iter()
                .try_fold(0usize, |n, c| n.checked_add(c.sequence.len()))
                .is_none_or(|n| n > 100_000)
        {
            return Err(Error::TooLarge);
        }
        if hash(
            &(crate::cloning::topo::graph_history::ENGINE_VERSION, graph),
            MAX_BYTES,
        )? != origin.selection.state_revision
        {
            return Err(bad("retained native graph checkpoint state changed"));
        }
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
    projection::validate_report(
        &product.preparation.substrates,
        &proof.physical.request,
        &proof.physical.limits,
        &proof.selection,
        &proof.options,
        proof,
    )
    .map_err(|e| bad(&e.to_string()))?;
    match (&proof.outcome, &product.reference, &product.annotations) {
        (projection::Outcome::Proven(view), Some(reference), Some(annotations)) => {
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
            annotation_sources_shape(annotations, &product.sources, reference.sequence.len())?;
            if reference.features != annotations.features || reference.original_genbank.is_some() {
                return Err(bad("canonical reference annotation evidence changed"));
            }
            let metadata = std::collections::BTreeMap::from([
                ("source:legacy:chemistry_unmapped".to_owned(), "true".to_owned()),
                ("topo:physical_graph_required".to_owned(), "true".to_owned()),
                ("topo:adapter_batch_revision".to_owned(), product.selection.batch_revision.clone()),
                (topo::RECIPE_KEY.to_owned(), serde_json::to_string(&(&product.selection.batch_revision, &product.selection, &product.sources)).map_err(|_| Error::TooLarge)?),
                ("topo:annotation_engine".to_owned(), topo::ANNOTATION_ENGINE_VERSION.to_owned()),
                ("topo:annotation_revision".to_owned(), product.batch_revision.clone()),
                ("topo:annotation_scope".to_owned(), "Source-coordinate correspondence only; no inherited chemistry or biological validity claim.".to_owned()),
            ]);
            if reference.metadata != metadata {
                return Err(bad(
                    "canonical reference restriction/correspondence metadata changed",
                ));
            }
            reference.validate()?;
        }
        (projection::Outcome::Withheld(_), None, None) => {}
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
                || !report.policy.allow_unapplied_steps && !view.all_requested_steps_applied
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

pub(crate) fn annotation_sources_shape(
    evidence: &crate::cloning::topo::annotations::AnnotationEvidence,
    sources: &[topo::SourceIdentity],
    length: usize,
) -> Result<(), Error> {
    if evidence.sources.len() != sources.len()
        || evidence.sources.len() > 64
        || evidence.features.len() > 4096
        || evidence.warnings.len() > 65
        || evidence.warnings.iter().any(|s| s.len() > 4096)
    {
        return Err(Error::TooLarge);
    }
    let mut feature_budget = collagene_model::bounds::ModelBudget::default();
    feature_budget.features(&evidence.features)?;
    let (mut maps, mut projected, mut retained) = (0usize, 0usize, 0usize);
    for (source, identity) in evidence.sources.iter().zip(sources) {
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
