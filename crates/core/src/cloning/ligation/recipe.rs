//! Ordered ligation from original-source preparations. Shared compilation emits
//! existing operations; only a fresh verifier pass supplies chemical evidence.
use crate::{
    chemistry::{Duplex, DuplexSummary, End, EndTreatment},
    cloning::assembly::recipe::{self as preparation, FragmentSpec, PreparedRow},
    planning::plan::{DesignPlan, Operation, OrientedSource, PlanInput, Source},
    verification::{self, evidence::VerificationReport},
};
use collagene_model::{Error, Molecule};
use serde::{Deserialize, Serialize};
pub const ENGINE_VERSION: &str = "ligation-recipe-1";
macro_rules! dto{($($t:tt)*)=>{#[derive(Clone,Debug,Serialize,Deserialize)]#[cfg_attr(feature="json-schema",derive(schemars::JsonSchema))]#[serde(deny_unknown_fields)]$($t)*};}
dto! {pub struct Fragment {pub source:FragmentSpec,pub treatment:Option<EndTreatment>}}
dto! {pub struct Recipe {pub name:String,pub fragments:Vec<Fragment>,pub circular:bool,pub allow_nicked:bool}}
dto! {pub struct Row {pub preparation:PreparedRow,pub treatment_step:Option<usize>,pub prepared_source:Source}}
dto! {pub struct Compiled {pub plan:DesignPlan,pub rows:Vec<Row>,pub ligation_step:usize}}
dto! {pub struct Junction {pub index:usize,pub upstream:usize,pub downstream:usize,pub closing:bool,pub upstream_end:Option<End>,pub downstream_end:Option<End>}}
dto! {pub struct Preview {pub context_revision:String,pub compiled:Compiled,pub verification:VerificationReport,
/// Original-row-indexed prepared references, before assembly orientation.
pub prepared:Vec<Option<Molecule>>,
/// Fresh physical summaries in assembly orientation; absent state is explicit.
pub physical:Vec<Option<DuplexSummary>>,pub product_physical:Option<DuplexSummary>,pub junctions:Vec<Junction>,pub preparation_complete:bool}}
fn bad(s: &str) -> Error {
    Error::Format(format!("ligation recipe: {s}"))
}
fn hex(s: &str) -> bool {
    s.len() == 64
        && s.bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}
fn sources(recipe: &Recipe) -> Result<Vec<FragmentSpec>, Error> {
    if recipe.fragments.is_empty() || recipe.fragments.len() > preparation::MAX_ROWS {
        return Err(Error::TooLarge);
    }
    crate::material::capsule::material_capsules_preflight(recipe.fragments.iter().filter_map(
        |row| match &row.source.preparation {
            preparation::Preparation::FromRun { origin, .. } => {
                Some(crate::material::capsule::MaterialOriginRef::Run(origin))
            }
            preparation::Preparation::FromDerived { origin, .. } => {
                Some(crate::material::capsule::MaterialOriginRef::Derived(origin))
            }
            _ => None,
        },
    ))?;
    for row in &recipe.fragments {
        if let preparation::Preparation::FromDerived {
            origin,
            preparation: child,
        } = &row.source.preparation
        {
            preparation::from_derived_shape(origin, child)?;
        }
    }
    Ok(recipe.fragments.iter().map(|f| f.source.clone()).collect())
}
/// Borrowed host guard before document copies, hashing or scientific work.
pub fn source_preflight(
    inputs: &[(&str, &Molecule)],
    target: &Molecule,
    recipe: &Recipe,
) -> Result<(), Error> {
    // Guard before cloning the small serialized recipe fields: large untrusted
    // primer/tail/name strings must not be copied merely to call the shared helper.
    if recipe.fragments.is_empty()
        || recipe.fragments.len() > preparation::MAX_ROWS
        || recipe.name.len() > 1024
    {
        return Err(Error::TooLarge);
    }
    let operations = recipe
        .fragments
        .iter()
        .try_fold(1usize, |total, fragment| {
            total
                .checked_add(fragment.source.preparation.operation_count())
                .and_then(|n| {
                    n.checked_add(
                        usize::from(fragment.source.generation.is_some())
                            + usize::from(fragment.treatment.is_some()),
                    )
                })
                .ok_or(Error::TooLarge)
        })?;
    if operations > 1000 {
        return Err(Error::TooLarge);
    }
    for f in &recipe.fragments {
        if f.source.role.len() > 128 || f.source.input_id.len() > 128 {
            return Err(Error::TooLarge);
        }
        if let preparation::Preparation::PcrCatalogDigest { selection, .. } = &f.source.preparation
        {
            selection.validate()?;
        }
        match &f.source.preparation {
            preparation::Preparation::OverlapExtensionPcrCatalogDigest {
                request,
                selection,
                fragment_index,
                ..
            } => {
                preparation::overlap_extension_catalog_shape(request, selection, *fragment_index)?;
            }
            preparation::Preparation::OverlapExtensionPcr { request } => {
                request.validate_shape()?
            }
            preparation::Preparation::MismatchPcrCatalogDigest {
                forward,
                reverse,
                options,
                selection,
                fragment_index,
                ..
            } => preparation::mismatch_pcr_catalog_shape(
                forward,
                reverse,
                options,
                selection,
                *fragment_index,
            )?,
            preparation::Preparation::MismatchPcr {
                forward,
                reverse,
                options,
            } => preparation::mismatch_pcr_shape(forward, reverse, options)?,
            preparation::Preparation::ExactPcr {
                forward, reverse, ..
            }
            | preparation::Preparation::PcrCatalogDigest {
                forward, reverse, ..
            } => {
                for p in [forward, reverse] {
                    if p.name.len() > 256 || p.annealing.len() > 1000 || p.tail.len() > 1000 {
                        return Err(Error::TooLarge);
                    }
                }
            }
            preparation::Preparation::FromRun {
                origin,
                preparation: child,
            } => preparation::from_run_shape(origin, child)?,

            preparation::Preparation::FromDerived {
                origin,
                preparation: child,
            } => preparation::from_derived_shape(origin, child)?,
            preparation::Preparation::AnnealedOligos {
                bottom_input_id,
                options,
                candidate_index,
            } => {
                preparation::anneal_shape(bottom_input_id, options, *candidate_index)?;
            }
            preparation::Preparation::CatalogDigest { selection, .. } => {
                selection.validate()?;
            }
            preparation::Preparation::Digest { enzymes, .. } => {
                if enzymes.len() > 32 || enzymes.iter().any(|s| s.len() > 128) {
                    return Err(Error::TooLarge);
                }
            }
            _ => {}
        }
    }
    crate::encoding::bounded_json::hash(recipe)?;
    preparation::preparation_source_preflight(inputs, target, &recipe.name, &sources(recipe)?)?;
    Ok(())
}
pub fn preflight(inputs: &[PlanInput], target: &Molecule, recipe: &Recipe) -> Result<(), Error> {
    if inputs.is_empty() || inputs.len() > preparation::MAX_ROWS {
        return Err(Error::TooLarge);
    }
    source_preflight(
        &inputs
            .iter()
            .map(|i| (i.id.as_str(), &i.molecule))
            .collect::<Vec<_>>(),
        target,
        recipe,
    )?;
    if inputs.iter().any(|i| !hex(&i.revision)) {
        return Err(bad("explicit canonical original revisions required"));
    }
    Ok(())
}
pub fn preview_revision(
    inputs: &[PlanInput],
    target: &Molecule,
    recipe: &Recipe,
    policy: preparation::CleavagePolicy,
) -> Result<String, Error> {
    preflight(inputs, target, recipe)?;
    preparation::bound_revisions(inputs, target)?;
    crate::encoding::bounded_json::hash(&(
        ENGINE_VERSION,
        "preview",
        inputs,
        target,
        recipe,
        policy,
    ))
}
pub fn compile(
    inputs: &[PlanInput],
    target: &Molecule,
    recipe: &Recipe,
    policy: preparation::CleavagePolicy,
) -> Result<Compiled, Error> {
    preflight(inputs, target, recipe)?;
    let prepared =
        preparation::compile_preparations(inputs, target, &recipe.name, &sources(recipe)?, policy)?;
    let mut steps = prepared.steps;
    let mut rows = Vec::with_capacity(prepared.rows.len());
    let mut fragments = Vec::with_capacity(prepared.rows.len());
    for (row, request) in prepared.rows.into_iter().zip(&recipe.fragments) {
        let (mut source, mut treatment_step) = (row.prepared_source.clone(), None);
        if let Some(treatment) = request.treatment {
            let index = steps.len();
            steps.push(Operation::TreatEnds {
                template: source,
                treatment,
            });
            source = Source::Step(index);
            treatment_step = Some(index);
        }
        fragments.push(OrientedSource {
            source: source.clone(),
            reverse: row.reverse,
        });
        rows.push(Row {
            preparation: row,
            treatment_step,
            prepared_source: source,
        });
    }
    let ligation_step = steps.len();
    steps.push(crate::cloning::ligation::outcomes::operation(
        fragments,
        recipe.circular,
        recipe.allow_nicked,
    ));
    let compiled = Compiled {
        plan: DesignPlan {
            name: recipe.name.clone(),
            target_revision: prepared.target_revision,
            inputs: prepared.inputs,
            steps,
        },
        rows,
        ligation_step,
    };
    crate::encoding::bounded_json::hash(&compiled)?;
    Ok(compiled)
}
fn summary(duplex: &Duplex, reverse: bool) -> Result<DuplexSummary, Error> {
    if reverse {
        duplex.reverse()?.summary()
    } else {
        duplex.summary()
    }
}
fn junctions(physical: &[Option<DuplexSummary>], rows: usize, circular: bool) -> Vec<Junction> {
    let count = rows.saturating_sub(1) + usize::from(circular);
    (0..count)
        .map(|index| {
            let upstream = index;
            let downstream = (index + 1) % rows;
            Junction {
                index,
                upstream,
                downstream,
                closing: downstream == 0,
                upstream_end: physical
                    .get(upstream)
                    .and_then(Option::as_ref)
                    .and_then(|s| s.ends.as_ref())
                    .map(|(_, right)| right.clone()),
                downstream_end: physical
                    .get(downstream)
                    .and_then(Option::as_ref)
                    .and_then(|s| s.ends.as_ref())
                    .map(|(left, _)| left.clone()),
            }
        })
        .collect()
}
pub fn preview(
    inputs: &[PlanInput],
    target: &Molecule,
    recipe: &Recipe,
    policy: preparation::CleavagePolicy,
) -> Result<Preview, Error> {
    let mut budget = crate::chemistry::tracked_duplex::CleavageBudget::new(100_000_000)?;
    Ok(preview_with_tracked(inputs, target, recipe, policy, &mut budget)?.0)
}
/// Same fresh verifier pass as public preview; private raw-orientation duplexes
/// cannot be reconstructed from the serialized summaries or molecule metadata.
pub(crate) fn preview_with_tracked(
    inputs: &[PlanInput],
    target: &Molecule,
    recipe: &Recipe,
    policy: preparation::CleavagePolicy,
    budget: &mut crate::chemistry::tracked_duplex::CleavageBudget,
) -> Result<
    (
        Preview,
        Vec<Option<crate::chemistry::tracked_duplex::TrackedDuplex>>,
    ),
    Error,
> {
    let context_revision = preview_revision(inputs, target, recipe, policy)?;
    let compiled = compile(inputs, target, recipe, policy)?;
    let mut scientific = verification::replay_budget::ReplayBudget::plan_local(budget);
    let (verification, mut trace) =
        verification::engine::verify_with_tracked_scoped(&compiled.plan, target, &mut scientific);
    let available = compiled
        .rows
        .iter()
        .map(|row| match &row.prepared_source {
            Source::Input(id) => compiled
                .plan
                .inputs
                .iter()
                .find(|input| &input.id == id)
                .map(|input| &input.molecule),
            Source::Step(index) => verification.product(*index),
        })
        .collect::<Vec<_>>();
    let models = available
        .iter()
        .filter_map(|model| *model)
        .collect::<Vec<_>>();
    collagene_model::bounds::model_preflight(&models)?;
    let mut bases = 0usize;
    for model in &models {
        bases = bases
            .checked_add(model.sequence.len())
            .filter(|n| *n <= preparation::MAX_ROLE_BASES)
            .ok_or(Error::TooLarge)?;
    }
    crate::encoding::bounded_json::hash(&models)?;
    let preparation_complete = available.iter().all(Option::is_some);
    let mut physical = Vec::with_capacity(available.len());
    let mut fresh: Vec<Option<crate::chemistry::tracked_duplex::TrackedDuplex>> =
        vec![None; compiled.rows.len()];
    for (index, (row, m)) in compiled.rows.iter().zip(&available).enumerate() {
        let Some(m) = m else {
            physical.push(None);
            continue;
        };
        let raw = match &row.prepared_source {
            Source::Step(i) => trace.get_mut(*i).and_then(Option::take),
            // Existing original intact-reference assumption only. Markers may
            // restrict it; serialized end/graph summaries never restore state.
            Source::Input(_) => Duplex::from_molecule(m)
                .and_then(|duplex| {
                    crate::chemistry::tracked_duplex::TrackedDuplex::from_unknown_duplex(&duplex)
                })
                .ok(),
        };
        physical.push(
            raw.as_ref()
                .map(|d| summary(&d.duplex, row.preparation.reverse))
                .transpose()?,
        );
        fresh[index] = raw;
    }
    let product_physical = trace
        .get(compiled.ligation_step)
        .and_then(Option::as_ref)
        .map(|d| d.duplex.summary())
        .transpose()?;
    let junctions = junctions(&physical, compiled.rows.len(), recipe.circular);
    let prepared = available.into_iter().map(|model| model.cloned()).collect();
    let result = Preview {
        context_revision,
        compiled,
        verification,
        prepared,
        physical,
        product_physical,
        junctions,
        preparation_complete,
    };
    result.validate()?;
    Ok((result, fresh))
}
fn same(a: &Source, b: &Source) -> bool {
    match (a, b) {
        (Source::Input(a), Source::Input(b)) => a == b,
        (Source::Step(a), Source::Step(b)) => a == b,
        _ => false,
    }
}
fn summary_bounds(s: &DuplexSummary) -> Result<(), Error> {
    if s.model.len() > 256
        || s.top_length > preparation::MAX_ROLE_BASES
        || s.bottom_length > preparation::MAX_ROLE_BASES
        || s.nicks.len() > 10_000
        || s.bottom_start.unsigned_abs() > preparation::MAX_ROLE_BASES as u64
    {
        return Err(Error::TooLarge);
    }
    if s.ends.as_ref().is_some_and(|(a, b)| {
        a.overhang.len() > preparation::MAX_ROLE_BASES
            || b.overhang.len() > preparation::MAX_ROLE_BASES
    }) {
        return Err(Error::TooLarge);
    }
    Ok(())
}
impl Preview {
    /// Readonly bounds/correspondence checks, never chemical replay or authority.
    pub fn validate(&self) -> Result<(), Error> {
        self.verification.validate_documentary()?;
        let rows = self.compiled.rows.len();
        let steps = self.compiled.plan.steps.len();
        if !hex(&self.context_revision)
            || rows == 0
            || rows > preparation::MAX_ROWS
            || steps > (preparation::MAX_STEPS_PER_ROW + 1) * preparation::MAX_ROWS + 1
            || self.compiled.ligation_step.checked_add(1) != Some(steps)
            || self.compiled.plan.inputs.len() > preparation::MAX_ROWS
            || self.verification.planned_steps != steps
            || self.prepared.len() != rows
            || self.physical.len() != self.prepared.len()
            || self.preparation_complete != self.prepared.iter().all(Option::is_some)
            || self.verification.step_reports.len() > steps
        {
            return Err(bad("invalid preview shape"));
        }
        for m in self
            .compiled
            .plan
            .inputs
            .iter()
            .map(|i| &i.molecule)
            .chain(self.prepared.iter().filter_map(Option::as_ref))
            .chain(
                self.verification
                    .steps
                    .iter()
                    .filter_map(|step| step.product.as_ref()),
            )
        {
            collagene_model::bounds::model_preflight(&[m])?;
        }
        for report in &self.verification.step_reports {
            report.validate_documentary()?;
            if report.step_index >= steps {
                return Err(Error::TooLarge);
            }
        }
        preparation::validate_prepared_inputs(
            &self.compiled.plan.inputs,
            &self.compiled.plan.steps,
        )?;
        if self.verification.plan_revision
            != crate::encoding::bounded_json::hash(&self.compiled.plan)?
            || self.verification.target_revision != self.compiled.plan.target_revision
        {
            return Err(bad("verification plan/target binding mismatch"));
        }
        let (fragments, circular, _allow_nicked) = match self.compiled.plan.steps.last() {
            Some(Operation::Ligate {
                fragments,
                circular,
                allow_nicked,
                options,
            }) => {
                if crate::encoding::bounded_json::hash(options)?
                    != crate::encoding::bounded_json::hash(
                        &crate::cloning::ligation::outcomes::Options::default(),
                    )?
                {
                    return Err(bad("compiled outcome bounds changed"));
                }
                (fragments, circular, allow_nicked)
            }
            _ => return Err(bad("preview final ligation/outcomes profile changed")),
        };
        if fragments.len() != rows
            || self.compiled.plan.steps[..self.compiled.ligation_step]
                .iter()
                .any(|s| {
                    !matches!(
                        s,
                        Operation::Pcr { .. }
                            | Operation::PcrWithMismatches { .. }
                            | Operation::OverlapExtension { .. }
                            | Operation::ReconstructRunProduct { .. }
                            | Operation::ReconstructMaterial { .. }
                            | Operation::AnnealOligos { .. }
                            | Operation::DeclareSourceTermini { .. }
                            | Operation::Digest { .. }
                            | Operation::CatalogDigest { .. }
                            | Operation::CheckedDigest { .. }
                            | Operation::DeclareSourceModifications { .. }
                            | Operation::RepairEnds { .. }
                            | Operation::TreatEnds { .. }
                    )
                })
        {
            return Err(bad("preview operation shape changed"));
        }
        let preparation_boundary = self
            .compiled
            .rows
            .iter()
            .filter_map(|row| row.treatment_step)
            .min()
            .unwrap_or(self.compiled.ligation_step);
        let preparations = self
            .compiled
            .rows
            .iter()
            .map(|row| row.preparation.clone())
            .collect::<Vec<_>>();
        preparation::validate_preparation_rows(
            &self.compiled.plan,
            &preparations,
            preparation_boundary,
        )?;
        for (index, (row, fragment)) in self.compiled.rows.iter().zip(fragments).enumerate() {
            let p = &row.preparation;
            preparation::validate_annealed_row(p, &self.compiled.plan.steps)?;
            if p.row_index != index
                || p.role.len() > 128
                || p.input_id.len() > 128
                || !hex(&p.input_revision)
                || !same(&fragment.source, &row.prepared_source)
                || fragment.reverse != p.reverse
            {
                return Err(bad("row source/orientation binding changed"));
            }
            let input = self
                .compiled
                .plan
                .inputs
                .iter()
                .find(|i| i.id == p.input_id)
                .ok_or_else(|| bad("row lost original source"))?;
            if input.revision != p.input_revision {
                return Err(bad("row original revision changed"));
            }
            if let Some(t) = row.treatment_step {
                let Some(Operation::TreatEnds { template, .. }) = self.compiled.plan.steps.get(t)
                else {
                    return Err(bad("row lost explicit treatment"));
                };
                if t >= self.compiled.ligation_step
                    || !same(template, &p.prepared_source)
                    || !same(&row.prepared_source, &Source::Step(t))
                {
                    return Err(bad("treatment source mapping changed"));
                }
            } else if !same(&row.prepared_source, &p.prepared_source) {
                return Err(bad("unexpected row preparation mapping"));
            }
        }
        for s in self
            .physical
            .iter()
            .filter_map(Option::as_ref)
            .chain(self.product_physical.iter())
        {
            summary_bounds(s)?;
        }
        let expected = junctions(&self.physical, rows, *circular);
        if crate::encoding::bounded_json::hash(&expected)?
            != crate::encoding::bounded_json::hash(&self.junctions)?
        {
            return Err(bad(
                "junction endpoints disagree with actual prepared summaries",
            ));
        }
        if self.product_physical.is_some()
            && self
                .verification
                .product(self.compiled.ligation_step)
                .is_none()
        {
            return Err(bad(
                "product physical evidence has no reconstructed product",
            ));
        }
        crate::encoding::bounded_json::hash(self)?;
        Ok(())
    }
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        self.validate()?;
        if self.context_revision != expected {
            return Err(bad("preview source/target/recipe context changed"));
        }
        Ok(())
    }
}
