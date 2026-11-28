//! Source-bound plan preflight and role binding.
use super::checks::RuleId;
use super::evidence::evidence_bytes;
use crate::{
    chemistry,
    planning::plan::{DesignPlan, Operation, Source, TopoSource},
};
use collagene_model::{Alphabet, Error, Molecule, Topology};
use std::cell::Cell;
use std::collections::BTreeSet;
/// Cheap borrowed snapshot/resource guard before plan hashing or source clones.
/// Actual sequence intervals, strand pairing and scalar engine policies are
/// independently checked by the source adapter against freshly bound molecules.
pub(super) fn retain_topo_source_incompleteness(
    preparation: &crate::cloning::topo::PreparationEvidence,
    index: usize,
    incomplete: &mut Vec<String>,
) {
    for origin in &preparation.graph_origins {
        if origin.incomplete {
            incomplete.push(format!("Step {} retains incomplete saved native graph {} source history; fresh scoped replay attests this flag, and independent target/policy acceptance does not establish completed source history or isolation.",index+1,origin.selection.record_hash));
        }
    }
}
pub(super) fn validate_topo_shape(
    sources: &[TopoSource],
    request: &crate::cloning::topo::Request,
) -> Result<(), Error> {
    crate::cloning::topo::physical::steps_preflight(&request.steps, &request.engine_limits)
        .map_err(|e| Error::Unsupported(format!("TOPO transcript: {e}")))?;
    if sources.len() > 32
        || request.recipe.schema_version != 1
        || request.recipe.materials.is_empty()
        || request.recipe.materials.len() > 16
        || request.steps.len() > 16
    {
        return Err(Error::TooLarge);
    }
    let mut roles = std::collections::BTreeSet::new();
    for source in sources {
        if source.role.len() > 256
            || source.role.trim().is_empty()
            || !roles.insert(source.role.as_str())
            || matches!(&source.source,Source::Input(id) if id.len()>128 || id.trim().is_empty())
        {
            return Err(Error::Format(
                "invalid or duplicate TOPO source role".into(),
            ));
        }
    }
    let mut strands = 0usize;
    let mut pairing = 0usize;
    let mut segments = 0usize;
    for material in &request.recipe.materials {
        match material {
            crate::cloning::topo::MaterialRecipe::SavedGraph {
                material_id,
                origin,
            } => {
                if material_id.trim().is_empty() || material_id.len() > 256 {
                    return Err(Error::TooLarge);
                }
                crate::material::capsule::native_origins_preflight([
                    crate::material::capsule::NativeOriginRef::TopoGraph(origin),
                ])?;
                strands = strands
                    .checked_add(origin.strand_count()?)
                    .filter(|n| *n <= 64)
                    .ok_or(Error::TooLarge)?;
            }
            crate::cloning::topo::MaterialRecipe::Prepared {
                declaration_id,
                assumption,
                strands: items,
                pairing: spans,
            } => {
                strands = strands.checked_add(items.len()).ok_or(Error::TooLarge)?;
                pairing = pairing.checked_add(spans.len()).ok_or(Error::TooLarge)?;
                if strands > 64
                    || pairing > 1024
                    || declaration_id.len() > 256
                    || assumption.len() > 4096
                {
                    return Err(Error::TooLarge);
                }
                for strand in items {
                    if strand.id.len() > 256
                        || strand.segments.is_empty()
                        || strand.segments.len() > 1024
                    {
                        return Err(Error::TooLarge);
                    }
                    if let crate::cloning::topo::physical::ThreePrime::Phosphotyrosyl {
                        activation_id,
                        ..
                    } = &strand.three_prime
                    {
                        if activation_id.len() > 256 {
                            return Err(Error::TooLarge);
                        }
                    }
                    segments = segments
                        .checked_add(strand.segments.len())
                        .ok_or(Error::TooLarge)?;
                    if segments > 2048 {
                        return Err(Error::TooLarge);
                    }
                    for segment in &strand.segments {
                        if segment.source_role.len() > 256
                            || !roles.contains(segment.source_role.as_str())
                        {
                            return Err(Error::Format(
                                "TOPO interval uses an unbound source role".into(),
                            ));
                        }
                    }
                }
            }
            crate::cloning::topo::MaterialRecipe::ExactPcr {
                material_id,
                template_role,
                forward,
                reverse,
                ..
            } => {
                strands = strands.checked_add(2).ok_or(Error::TooLarge)?;
                if strands > 64
                    || material_id.len() > 256
                    || template_role.len() > 256
                    || !roles.contains(template_role.as_str())
                {
                    return Err(Error::Format(
                        "TOPO PCR uses an unbound role or exceeds preparation bounds".into(),
                    ));
                }
                for primer in [forward, reverse] {
                    if primer.name.len() > 256
                        || primer.annealing.len() > 1000
                        || primer.tail.len() > 1000
                    {
                        return Err(Error::TooLarge);
                    }
                }
            }
        }
    }
    evidence_bytes(request, 1024 * 1024)?;
    Ok(())
}

/// Borrowed literal request bounds; actual source geometry is independently
/// reconstructed only after source models and aggregate resource caps pass.

/// Bind the exact ordered plan sources to explicit recipe roles. Display names
/// and imported producer metadata cannot resolve a role or establish chemistry.
pub fn topo_inputs(
    sources: &[TopoSource],
    molecules: &[&Molecule],
    request: &crate::cloning::topo::Request,
) -> Result<Vec<crate::cloning::topo::Input>, Error> {
    validate_topo_shape(sources, request)?;
    if sources.len() != molecules.len() {
        return Err(Error::Format("TOPO source count mismatch".into()));
    }
    collagene_model::bounds::model_preflight(molecules)?;
    let mut bases = 0usize;
    for molecule in molecules {
        bases = bases
            .checked_add(molecule.sequence.len())
            .filter(|&n| n <= 100_000)
            .ok_or(Error::TooLarge)?;
        if molecule.sequence.alphabet() != Alphabet::Dna {
            return Err(Error::Format("TOPO sources must be DNA".into()));
        }
    }
    let inputs = sources
        .iter()
        .zip(molecules)
        .map(|(source, molecule)| crate::cloning::topo::Input {
            role: source.role.clone(),
            molecule: (*molecule).clone(),
        })
        .collect::<Vec<_>>();
    crate::cloning::topo::source_preflight(&inputs, request)
        .map_err(|error| Error::Format(format!("TOPO source preparation rejected: {error}")))?;
    Ok(inputs)
}

/// A declaration may constrain an original intact reference hypothesis only.
/// This guard never interprets producer summaries as chemical authority.
pub fn declared_source_termini_preflight(source: &Molecule) -> Result<(), Error> {
    collagene_model::bounds::model_preflight(&[source])?;
    source.validate()?;
    if source.topology != Topology::Linear
        || source.sequence.alphabet() != Alphabet::Dna
        || source
            .sequence
            .as_str()
            .bytes()
            .any(|b| !b"ACGT".contains(&b))
    {
        return Err(Error::Unsupported("Declared source termini require an original intact linear unambiguous DNA reference; topology or missing calls are not repaired.".into()));
    }
    let raw_producer_marker = source.original_genbank.as_ref().is_some_and(|raw| {
        raw.lines()
            .filter_map(|line| {
                line.strip_prefix("COMMENT     ")
                    .or_else(|| line.strip_prefix("CC   "))
            })
            .any(|value| {
                matches!(
                    value.trim_end(),
                    "Collagene physical product v1: unverified sequence-only evidence"
                        | "Collagene reference-only v1: unverified sequence-only evidence"
                )
            })
    });
    if raw_producer_marker
        || chemistry::has_unmapped_physical_evidence(source)
        || source.metadata.keys().any(|key| {
            matches!(
                key.as_str(),
                "tracked_modifications"
                    | "source:declared_termini"
                    | "source:legacy:chemistry_unmapped"
                    | "source:snapgene:chemistry_unmapped"
                    | "topo:physical_graph_required"
                    | "source:recombination:reference_only"
                    | "source:recombination:unverified_recipe"
                    | "source:restriction_inspection:recipe"
            ) || key == crate::cloning::recombination::REFERENCE_ONLY_KEY
                || key == crate::cloning::recombination::RECIPE_KEY
                || key == crate::cloning::restriction::inspection::RECIPE_KEY
                || key.starts_with("variant:")
                || key.starts_with("topo:")
                || key.starts_with("source:restriction_inspection:")
                || key.starts_with("restriction_catalog:")
                || key.starts_with("source:restriction_catalog:")
        })
    {
        return Err(Error::Unsupported("Source terminal declarations cannot promote imported physical, mixed/single-strand, variant or reference-only evidence. Reconstruct the actual originating graph; this declaration is not materialization, nick repair or chemical recovery.".into()));
    }
    Ok(())
}

/// Current per-operation semantic validation. Global resource admission is
/// separate so one rejected profile does not suppress independent branches.
pub(super) struct ReferenceChecks {
    pinned_catalog: Option<crate::cloning::restriction::catalog::catalog::Catalog>,
}
impl ReferenceChecks {
    pub(super) fn new(plan: &DesignPlan) -> Result<Self, Error> {
        let count = plan
            .steps
            .iter()
            .filter(|step| {
                matches!(
                    step,
                    Operation::ResectThreePrime { .. }
                        | Operation::HomologyAnnealing { .. }
                        | Operation::HomologyAssembly { .. }
                        | Operation::Topo { .. }
                        | Operation::CheckedTopo { .. }
                        | Operation::MultiSiteRecombineReference { .. }
                        | Operation::TopoActions { .. }
                        | Operation::RecombineReference { .. }
                )
            })
            .count();
        if count > 64 {
            return Err(Error::TooLarge);
        }
        let pinned_catalog = if plan.steps.iter().any(|step| {
            matches!(
                step,
                Operation::CatalogDigest { .. } | Operation::CheckedDigest { .. }
            )
        }) {
            Some(crate::cloning::restriction::catalog::builtin_catalog()?)
        } else {
            None
        };
        Ok(Self { pinned_catalog })
    }
    pub(super) fn validate_step(&self, step: &Operation, rule: &Cell<RuleId>) -> Result<(), Error> {
        if let Operation::CatalogDigest {
            selection,
            conditions: Some(conditions),
            ..
        }
        | Operation::CheckedDigest {
            selection,
            conditions: Some(conditions),
            ..
        } = step
        {
            rule.set(RuleId::ReactionConditions);
            conditions.validate_for(
                self.pinned_catalog.as_ref().ok_or(Error::TooLarge)?,
                selection,
            )?;
        }
        if let Operation::CheckedDigest {
            selection,
            options,
            fragment_index,
            ..
        } = step
        {
            rule.set(RuleId::RestrictionProfile);
            selection.validate_for(
                self.pinned_catalog
                    .as_ref()
                    .ok_or_else(|| Error::Format("missing checked catalog".into()))?,
            )?;
            options.validate()?;
            if *fragment_index >= options.max_fragments {
                return Err(Error::TooLarge);
            }
        }
        if let Operation::ResectThreePrime { cuts, .. } = step {
            rule.set(RuleId::EndRepair);
            if cuts.left > crate::cloning::overlap::resection::MAX_SOURCE_BASES / 2
                || cuts.right > crate::cloning::overlap::resection::MAX_SOURCE_BASES / 2
            {
                return Err(Error::TooLarge);
            }
        }
        if let Operation::HomologyAnnealing {
            fragments,
            options,
            policy,
        } = step
        {
            rule.set(RuleId::HomologyPolicy);
            policy.validate_shape(fragments.len())?;
            policy.validate_profile(options.method)?;
            options.validate()?;
            if fragments.len() == 1 && policy.topology != Topology::Circular {
                return Err(Error::Unsupported("Single linear3′ processing uses ResectThreePrime, not a joining assembly operation".into()));
            }
            if fragments
                .iter()
                .any(|f| matches!(&f.source,Source::Input(id)if id.trim().is_empty()||id.len()>128))
            {
                return Err(Error::Format("Invalid annealing source identity".into()));
            }
        }
        if let Operation::HomologyAssembly {
            fragments,
            options,
            policy,
        } = step
        {
            rule.set(RuleId::HomologyPolicy);
            policy.validate_shape(fragments.len())?;

            options.validate()?;
            if options.method == crate::cloning::overlap::profiles::OverlapMethod::InFusion {
                return Err(Error::Unsupported("The declared 5′ resection/fill/seal model cannot represent In-Fusion 3′ processing and host rescue; retain its reference-only interpretation.".into()));
            }
            if fragments.iter().any(
                |f| matches!(&f.source, Source::Input(id) if id.trim().is_empty() || id.len()>128),
            ) {
                return Err(Error::Format(
                    "Invalid physical homology source identity".into(),
                ));
            }
        }
        if let Operation::SelectedGoldenGate {
            templates,
            options,
            candidate_index,
            directed_sequence_key,
            physical_sequence_key,
        } = step
        {
            rule.set(RuleId::GoldenGatePolicy);
            if templates.is_empty()
                || templates.len() > 32
                || options.min_fragments == 0
                || options.min_fragments > options.max_fragments
                || options.max_fragments > 32
                || options.max_products == 0
                || options.max_products > 1000
                || options.max_states == 0
                || options.max_states > 1_000_000
                || *candidate_index >= options.max_products
                || [directed_sequence_key, physical_sequence_key]
                    .iter()
                    .any(|key| {
                        key.len() != 64
                            || !key
                                .bytes()
                                .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
                    })
                || templates.iter().any(
                    |source| matches!(source, Source::Input(id) if id.is_empty() || id.len() > 128),
                )
            {
                return Err(Error::Format(
                    "Invalid bounded key-bound Golden Gate selection".into(),
                ));
            }
        }
        if let Operation::CatalogDigest {
            template,
            selection,
            fragment_index,
            ..
        } = step
        {
            rule.set(RuleId::RestrictionProfile);
            selection
                .validate_for(self.pinned_catalog.as_ref().ok_or(Error::TooLarge)?)
                .map_err(|error| Error::Format(format!("Catalog digest pin rejected: {error}")))?;
            if *fragment_index > 100_000
                || matches!(template,Source::Input(id) if id.is_empty()||id.len()>128)
            {
                return Err(Error::TooLarge);
            }
        }
        if let Operation::DeclareSourceTermini {
            template: Source::Input(id),
            ..
        } = step
        {
            rule.set(RuleId::SourceChemistry);
            if id.is_empty() || id.len() > 128 {
                return Err(Error::TooLarge);
            }
        }
        if let Operation::Topo {
            sources,
            request,
            selection,
        }
        | Operation::CheckedTopo {
            sources,
            request,
            selection,
            ..
        } = step
        {
            rule.set(RuleId::TopoDeclaration);
            validate_topo_shape(sources, request)?;
            if selection.batch_revision.len() != 64
                || !selection
                    .batch_revision
                    .bytes()
                    .all(|b| b.is_ascii_hexdigit())
                || selection.component >= 64
            {
                return Err(Error::Format("invalid TOPO source selection".into()));
            }
        }
        if let Operation::MultiSiteRecombineReference {
            sources,
            request,
            selection,
        } = step
        {
            rule.set(RuleId::RecombinationReference);
            if !(2..=5).contains(&sources.len()) {
                return Err(Error::TooLarge);
            }
            let mut roles = BTreeSet::new();
            for source in sources {
                if source.role.trim().is_empty()
                    || source.role.len() > 256
                    || !roles.insert(source.role.as_str())
                {
                    return Err(Error::Format(
                        "Multi-site roles must be unique and bounded".into(),
                    ));
                }
            }
            crate::cloning::recombination::multisite::request_preflight(request)?;
            if selection.batch_revision.len() != 64
                || !selection
                    .batch_revision
                    .bytes()
                    .all(|b| b.is_ascii_hexdigit())
                || selection.candidate_id >= request.reaction.limits.max_candidates
                || selection.component_id >= request.reaction.limits.max_components_per_candidate
            {
                return Err(Error::Format("Invalid multi-site selection".into()));
            }
        }
        if let Operation::TopoActions { sources, request } = step {
            rule.set(RuleId::TopoDeclaration);
            validate_topo_shape(sources, &request.preparation)?;
            crate::cloning::topo::actions::request_preflight(request)?;
        }

        if let Operation::RecombineReference {
            sources,
            request,
            selection,
        } = step
        {
            rule.set(RuleId::RecombinationReference);
            if selection.batch_revision.len() != 64
                || !selection
                    .batch_revision
                    .bytes()
                    .all(|b| b.is_ascii_hexdigit())
            {
                return Err(Error::TooLarge);
            }
            request.validate()?;
            if selection.candidate_id >= request.limits.max_candidates
                || selection.component_id >= request.limits.max_components_per_candidate
                || sources.iter().any(
                    |s| matches!(s, Source::Input(id) if id.len() > 128 || id.trim().is_empty()),
                )
            {
                return Err(Error::Format(
                    "invalid reference-exchange selection or source shape".into(),
                ));
            }
        }
        Ok(())
    }
}

pub(crate) fn validate_reference_steps(plan: &DesignPlan) -> Result<(), Error> {
    crate::material::capsule::plan_preflight(plan)?;
    if plan.steps.len() > 1000
        || plan.inputs.len() > 1000
        || plan.name.len() > 1024
        || plan.target_revision.len() != 64
    {
        return Err(Error::TooLarge);
    }
    let checks = ReferenceChecks::new(plan)?;
    let rule = Cell::new(RuleId::Admission);
    for step in &plan.steps {
        checks.validate_step(step, &rule)?;
    }
    // Portable candidates preserve import-scale annotations and original records.
    // Scientific replay applies its tighter shared computation budget in admit().
    collagene_model::bounds::native_frame_preflight(
        plan.inputs.iter().map(|input| &input.molecule),
        crate::formats::sequence::MAX_IMPORT_BYTES,
    )?;
    Ok(())
}

/// Roles bind ordered slots and typed ancestry, independently of molecule display names.
pub fn recombination_inputs(
    sources: &[Source; 2],
    molecules: [&Molecule; 2],
    request: &crate::cloning::recombination::Request,
) -> Result<Vec<crate::cloning::recombination::reference::Input>, Error> {
    crate::cloning::recombination::source_preflight(molecules, request)?;
    if sources
        .iter()
        .any(|s| matches!(s, Source::Input(id) if id.len() > 128 || id.trim().is_empty()))
    {
        return Err(Error::TooLarge);
    }
    sources
        .iter()
        .zip(molecules)
        .enumerate()
        .map(|(slot, (source, molecule))| {
            let bytes = serde_json::to_vec(&("plan-reference-exchange-source-v1", slot, source))
                .map_err(|e| Error::Format(e.to_string()))?;
            Ok(crate::cloning::recombination::reference::Input {
                role: blake3::hash(&bytes).to_hex().to_string(),
                molecule: molecule.clone(),
            })
        })
        .collect()
}
