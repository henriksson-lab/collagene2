//! Ordered single-base TA/GC tailing from original-source preparations. Compilation emits
//! existing operations; only a fresh verifier pass supplies chemical evidence.
use crate::{
    Error, Molecule,
    assembly_recipe::{self as preparation, FragmentSpec, PreparedRow},
    chemistry::{Duplex, DuplexSummary, EndTreatment},
    tailing::TailingBase,
    verification::{
        self, DesignPlan, Operation, OrientedSource, PlanInput, Source, VerificationReport,
    },
};
use serde::{Deserialize, Serialize};
pub const ENGINE_VERSION: &str = "tailing-recipe-1";
macro_rules! dto{($($t:tt)*)=>{#[derive(Clone,Debug,Serialize,Deserialize)]#[cfg_attr(feature="json-schema",derive(schemars::JsonSchema))]#[serde(deny_unknown_fields)]$($t)*};}
dto! {pub struct Fragment {pub source:FragmentSpec,pub base:TailingBase,pub treatment:Option<EndTreatment>}}
dto! {pub struct Recipe {pub name:String,pub fragments:Vec<Fragment>,pub circular:bool,pub allow_nicked:bool}}
dto! {pub struct Row {pub preparation:PreparedRow,pub tail_step:usize,pub base:TailingBase,pub treatment_step:Option<usize>,pub prepared_source:Source}}
dto! {pub struct Compiled {pub plan:DesignPlan,pub rows:Vec<Row>,pub ligation_step:usize}}
pub use crate::ligation_recipe::Junction;
dto! {pub struct Preview {pub context_revision:String,pub recipe:Recipe,pub compiled:Compiled,pub verification:VerificationReport,
/// Leading prepared reference-document prefix, before assembly orientation.
pub prepared:Vec<Molecule>,
/// Fresh physical summaries in assembly orientation; absent state is explicit.
pub physical:Vec<Option<DuplexSummary>>,pub product_physical:Option<DuplexSummary>,pub junctions:Vec<Junction>,pub preparation_complete:bool}}
fn bad(s: &str) -> Error {
    Error::Format(format!("tailing recipe: {s}"))
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
    Ok(recipe.fragments.iter().map(|f| f.source.clone()).collect())
}
fn recipe_shape(recipe: &Recipe) -> Result<(), Error> {
    // Guard before cloning the small serialized recipe fields: large untrusted
    // primer/tail/name strings must not be copied merely to call the shared helper.
    if recipe.fragments.is_empty()
        || recipe.fragments.len() > preparation::MAX_ROWS
        || recipe.name.len() > 1024
    {
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
            preparation::Preparation::ExactPcr { forward, reverse }
            | preparation::Preparation::PcrCatalogDigest {
                forward, reverse, ..
            } => {
                for p in [forward, reverse] {
                    if p.name.len() > 256 || p.annealing.len() > 1000 || p.tail.len() > 1000 {
                        return Err(Error::TooLarge);
                    }
                }
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
    preparation::hash(recipe)?;
    Ok(())
}
/// Borrowed host guard before document copies, hashing or scientific work.
pub fn source_preflight(
    inputs: &[(&str, &Molecule)],
    target: &Molecule,
    recipe: &Recipe,
) -> Result<(), Error> {
    recipe_shape(recipe)?;
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
fn context_revision(
    inputs: &[PlanInput],
    target_revision: &str,
    recipe: &Recipe,
) -> Result<String, Error> {
    preparation::hash(&(
        ENGINE_VERSION,
        verification::ENGINE_VERSION,
        "preview",
        inputs,
        target_revision,
        recipe,
    ))
}
pub fn preview_revision(
    inputs: &[PlanInput],
    target: &Molecule,
    recipe: &Recipe,
) -> Result<String, Error> {
    preflight(inputs, target, recipe)?;
    let target_revision = preparation::bound_revisions(inputs, target)?;
    context_revision(inputs, &target_revision, recipe)
}
pub fn compile(
    inputs: &[PlanInput],
    target: &Molecule,
    recipe: &Recipe,
) -> Result<Compiled, Error> {
    preflight(inputs, target, recipe)?;
    let prepared =
        preparation::compile_preparations(inputs, target, &recipe.name, &sources(recipe)?)?;
    let mut steps = prepared.steps;
    let mut rows = Vec::with_capacity(prepared.rows.len());
    let mut fragments = Vec::with_capacity(prepared.rows.len());
    for (row, request) in prepared.rows.into_iter().zip(&recipe.fragments) {
        let tail_step = steps.len();
        steps.push(Operation::Tail3Prime {
            template: row.prepared_source.clone(),
            base: request.base,
        });
        let (mut source, mut treatment_step) = (Source::Step(tail_step), None);
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
            tail_step,
            base: request.base,
            treatment_step,
            prepared_source: source,
        });
    }
    let ligation_step = steps.len();
    steps.push(Operation::Ligate {
        fragments,
        circular: recipe.circular,
        allow_nicked: recipe.allow_nicked,
    });
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
    preparation::hash(&compiled)?;
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
pub fn preview(inputs: &[PlanInput], target: &Molecule, recipe: &Recipe) -> Result<Preview, Error> {
    Ok(preview_with_physical(inputs, target, recipe)?.0)
}
/// Same fresh verifier pass as public preview; private raw-orientation duplexes
/// cannot be reconstructed from the serialized summaries or molecule metadata.
pub(crate) fn preview_with_physical(
    inputs: &[PlanInput],
    target: &Molecule,
    recipe: &Recipe,
) -> Result<(Preview, Vec<Option<Duplex>>), Error> {
    let context_revision = preview_revision(inputs, target, recipe)?;
    let compiled = compile(inputs, target, recipe)?;
    let (verification, mut trace) = verification::verify_with_physical(&compiled.plan, target);
    let mut available = Vec::with_capacity(compiled.rows.len());
    for row in &compiled.rows {
        let m = match &row.prepared_source {
            Source::Input(id) => compiled
                .plan
                .inputs
                .iter()
                .find(|i| &i.id == id)
                .map(|i| &i.molecule),
            Source::Step(i) => verification.products.get(*i),
        };
        if let Some(m) = m {
            available.push(m)
        } else {
            break;
        }
    }
    crate::model::bounds::model_preflight(&available)?;
    let mut bases = 0usize;
    for m in &available {
        bases = bases
            .checked_add(m.sequence.len())
            .filter(|n| *n <= preparation::MAX_ROLE_BASES)
            .ok_or(Error::TooLarge)?;
    }
    preparation::hash(&available)?;
    let preparation_complete = available.len() == compiled.rows.len();
    let mut physical = Vec::with_capacity(available.len());
    let mut fresh: Vec<Option<Duplex>> = vec![None; compiled.rows.len()];
    for (index, (row, m)) in compiled.rows.iter().zip(&available).enumerate() {
        let raw = match &row.prepared_source {
            Source::Step(i) => trace.get_mut(*i).and_then(Option::take),
            // Existing original intact-reference assumption only. Markers may
            // restrict it; serialized end/graph summaries never restore state.
            Source::Input(_) => Duplex::from_molecule(m).ok(),
        };
        physical.push(
            raw.as_ref()
                .map(|d| summary(d, row.preparation.reverse))
                .transpose()?,
        );
        fresh[index] = raw;
    }
    let product_physical = trace
        .get(compiled.ligation_step)
        .and_then(Option::as_ref)
        .map(Duplex::summary)
        .transpose()?;
    let junctions = junctions(&physical, compiled.rows.len(), recipe.circular);
    let prepared = available.into_iter().cloned().collect();
    let result = Preview {
        context_revision,
        recipe: recipe.clone(),
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
/// Structural correspondence only: compare each declared preparation operation
/// with the common compiler's retained prefix, without simulating any DNA.
fn validate_preparation_prefix(
    compiled: &Compiled,
    recipe: &Recipe,
    tail_start: usize,
) -> Result<(), Error> {
    use preparation::Preparation;
    let mut next = 0usize;
    let mut check = |expected: Operation| -> Result<usize, Error> {
        let index = next;
        if index >= tail_start
            || preparation::hash(
                compiled
                    .plan
                    .steps
                    .get(index)
                    .ok_or_else(|| bad("missing preparation operation"))?,
            )? != preparation::hash(&expected)?
        {
            return Err(bad(
                "preparation operation differs from declared original-source recipe",
            ));
        }
        next += 1;
        Ok(index)
    };
    for (i, (row, request)) in compiled.rows.iter().zip(&recipe.fragments).enumerate() {
        let spec = &request.source;
        let original = Source::Input(spec.input_id.clone());
        let operation = match &spec.preparation {
            Preparation::WholeLinear => None,
            Preparation::DeclaredDuplex { left, right } => Some(Operation::DeclareSourceTermini {
                template: original.clone(),
                left: *left,
                right: *right,
            }),
            Preparation::Digest {
                enzymes,
                fragment_index,
                ..
            } => Some(Operation::Digest {
                template: original.clone(),
                enzymes: enzymes.clone(),
                fragment_index: *fragment_index,
            }),
            Preparation::CatalogDigest {
                selection,
                fragment_index,
                ..
            } => Some(Operation::CatalogDigest {
                template: original.clone(),
                selection: selection.clone(),
                fragment_index: *fragment_index,
            }),
            Preparation::ExactPcr { forward, reverse }
            | Preparation::PcrCatalogDigest {
                forward, reverse, ..
            } => Some(Operation::Pcr {
                template: original.clone(),
                forward: forward.clone(),
                reverse: reverse.clone(),
            }),
        };
        let preparation_step = operation.map(&mut check).transpose()?;
        let mut source = preparation_step.map(Source::Step).unwrap_or(original);
        if let Preparation::PcrCatalogDigest {
            selection,
            fragment_index,
            ..
        } = &spec.preparation
        {
            source = Source::Step(check(Operation::CatalogDigest {
                template: source,
                selection: selection.clone(),
                fragment_index: *fragment_index,
            })?);
        }
        let repair = match &spec.preparation {
            Preparation::Digest { repair, .. }
            | Preparation::CatalogDigest { repair, .. }
            | Preparation::PcrCatalogDigest { repair, .. } => *repair,
            _ => None,
        };
        let repair_step = repair
            .map(|method| {
                check(Operation::RepairEnds {
                    template: source.clone(),
                    method,
                })
            })
            .transpose()?;
        if let Some(index) = repair_step {
            source = Source::Step(index);
        }
        let input = compiled
            .plan
            .inputs
            .iter()
            .find(|p| p.id == spec.input_id)
            .ok_or_else(|| bad("preparation input missing"))?;
        let expected = PreparedRow {
            row_index: i,
            role: spec.role.clone(),
            input_id: spec.input_id.clone(),
            input_revision: input.revision.clone(),
            preparation_step,
            repair_step,
            prepared_source: source,
            reverse: spec.reverse,
        };
        if preparation::hash(&expected)? != preparation::hash(&row.preparation)? {
            return Err(bad(
                "prepared row metadata differs from its declared operations",
            ));
        }
    }
    if next != tail_start {
        return Err(bad("unclaimed preparation operations"));
    }
    Ok(())
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
        recipe_shape(&self.recipe)?;
        let rows = self.compiled.rows.len();
        let steps = self.compiled.plan.steps.len();
        if !hex(&self.context_revision)
            || rows == 0
            || rows > preparation::MAX_ROWS
            || steps > 5 * preparation::MAX_ROWS + 1
            || self.compiled.ligation_step.checked_add(1) != Some(steps)
            || self.compiled.plan.inputs.len() > preparation::MAX_ROWS
            || self.verification.products.len() > steps
            || self.prepared.len() > rows
            || self.physical.len() != self.prepared.len()
            || self.preparation_complete != (self.prepared.len() == rows)
            || self.verification.step_reports.len() > steps
            || self.verification.annotation_reports.len() > steps
        {
            return Err(bad("invalid preview shape"));
        }
        for m in self
            .compiled
            .plan
            .inputs
            .iter()
            .map(|i| &i.molecule)
            .chain(self.prepared.iter())
            .chain(self.verification.products.iter())
        {
            crate::model::bounds::model_preflight(&[m])?;
        }
        for report in &self.verification.annotation_reports {
            if report.step_index >= steps || report.evidence.sources.len() > preparation::MAX_ROWS {
                return Err(Error::TooLarge);
            }
            let mut budget = crate::model::bounds::ModelBudget::default();
            budget.features(&report.evidence.features)?;
            for role in &report.evidence.sources {
                if role.unmapped.len() > 4096 {
                    return Err(Error::TooLarge);
                }
                for f in &role.unmapped {
                    budget.features(std::slice::from_ref(&f.feature))?;
                }
            }
        }
        for input in &self.compiled.plan.inputs {
            if preparation::hash(&input.molecule)? != input.revision {
                return Err(bad("original source content/revision differs"));
            }
        }
        if self.verification.engine_version != verification::ENGINE_VERSION
            || self.recipe.fragments.len() != rows
            || self.compiled.plan.name != self.recipe.name
            || !hex(&self.compiled.plan.target_revision)
            || self.context_revision
                != context_revision(
                    &self.compiled.plan.inputs,
                    &self.compiled.plan.target_revision,
                    &self.recipe,
                )?
            || self.verification.plan_revision != preparation::hash(&self.compiled.plan)?
            || self.verification.target_revision != self.compiled.plan.target_revision
        {
            return Err(bad("verification plan/target binding mismatch"));
        }
        let Some(Operation::Ligate {
            fragments,
            circular,
            allow_nicked,
        }) = self.compiled.plan.steps.last()
        else {
            return Err(bad("preview has no final ordered ligation operation"));
        };
        if *circular != self.recipe.circular
            || *allow_nicked != self.recipe.allow_nicked
            || fragments.len() != rows
            || self.compiled.plan.steps[..self.compiled.ligation_step]
                .iter()
                .any(|s| {
                    !matches!(
                        s,
                        Operation::Pcr { .. }
                            | Operation::DeclareSourceTermini { .. }
                            | Operation::Digest { .. }
                            | Operation::CatalogDigest { .. }
                            | Operation::RepairEnds { .. }
                            | Operation::TreatEnds { .. }
                            | Operation::Tail3Prime { .. }
                    )
                })
        {
            return Err(bad("preview operation shape changed"));
        }
        let tail_start = self
            .compiled
            .rows
            .first()
            .map(|r| r.tail_step)
            .ok_or_else(|| bad("no tail rows"))?;
        if tail_start > self.compiled.ligation_step {
            return Err(bad("tail prefix outside plan"));
        }
        validate_preparation_prefix(&self.compiled, &self.recipe, tail_start)?;
        let mut next_tail = tail_start;
        if self.compiled.plan.steps[..tail_start].iter().any(|op| {
            matches!(
                op,
                Operation::Tail3Prime { .. } | Operation::TreatEnds { .. }
            )
        }) {
            return Err(bad("unexpected tail/treatment in preparation prefix"));
        }

        for (index, (row, fragment)) in self.compiled.rows.iter().zip(fragments).enumerate() {
            let request = &self.recipe.fragments[index];
            let p = &row.preparation;
            if p.row_index != index
                || p.role.len() > 128
                || p.input_id.len() > 128
                || p.input_id != request.source.input_id
                || p.role != request.source.role
                || p.reverse != request.source.reverse
                || row.base != request.base
                || row.tail_step != next_tail
                || row.tail_step >= self.compiled.ligation_step
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
            let Some(Operation::Tail3Prime { template, base }) =
                self.compiled.plan.steps.get(row.tail_step)
            else {
                return Err(bad("row lost its explicit tail operation"));
            };
            if *base != request.base
                || !same(template, &p.prepared_source)
                || matches!(&p.prepared_source,Source::Step(i) if *i>=tail_start)
            {
                return Err(bad("tail base or prepared-source binding changed"));
            }
            next_tail = next_tail.checked_add(1).ok_or(Error::TooLarge)?;
            if let Some(t) = row.treatment_step {
                let Some(Operation::TreatEnds {
                    template,
                    treatment,
                }) = self.compiled.plan.steps.get(t)
                else {
                    return Err(bad("row lost explicit treatment"));
                };
                if t >= self.compiled.ligation_step
                    || t != next_tail
                    || Some(*treatment) != request.treatment
                    || !same(template, &Source::Step(row.tail_step))
                    || !same(&row.prepared_source, &Source::Step(t))
                {
                    return Err(bad("treatment source mapping changed"));
                }
                next_tail = next_tail.checked_add(1).ok_or(Error::TooLarge)?;
            } else if request.treatment.is_some()
                || !same(&row.prepared_source, &Source::Step(row.tail_step))
            {
                return Err(bad("unexpected tail/treatment source mapping"));
            }
        }
        if next_tail != self.compiled.ligation_step {
            return Err(bad("extra or missing tail/treatment steps"));
        }
        let expected_prepared = self
            .compiled
            .rows
            .iter()
            .map(|r| match &r.prepared_source {
                Source::Step(i) => self.verification.products.get(*i),
                Source::Input(id) => self
                    .compiled
                    .plan
                    .inputs
                    .iter()
                    .find(|i| &i.id == id)
                    .map(|i| &i.molecule),
            })
            .take_while(Option::is_some)
            .collect::<Vec<_>>();
        if expected_prepared.len() != self.prepared.len()
            || expected_prepared
                .iter()
                .zip(&self.prepared)
                .any(|(a, b)| *a != Some(b))
        {
            return Err(bad(
                "prepared documents differ from retained reconstructed row products",
            ));
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
        if preparation::hash(&expected)? != preparation::hash(&self.junctions)? {
            return Err(bad(
                "junction endpoints disagree with actual prepared summaries",
            ));
        }
        if self.product_physical.is_some()
            && self.verification.products.len() <= self.compiled.ligation_step
        {
            return Err(bad(
                "product physical evidence has no reconstructed product",
            ));
        }
        preparation::hash(self)?;
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
