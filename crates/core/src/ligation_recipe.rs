//! Ordered ligation from original-source preparations. Shared compilation emits
//! existing operations; only a fresh verifier pass supplies chemical evidence.
use crate::{
    Error, Molecule,
    assembly_recipe::{self as preparation, FragmentSpec, PreparedRow},
    chemistry::{Duplex, DuplexSummary, End, EndTreatment},
    verification::{
        self, DesignPlan, Operation, OrientedSource, PlanInput, Source, VerificationReport,
    },
};
use serde::{Deserialize, Serialize};
pub const ENGINE_VERSION: &str = "ligation-recipe-1";
macro_rules! dto{($($t:tt)*)=>{#[derive(Clone,Debug,Serialize,Deserialize)]#[cfg_attr(feature="json-schema",derive(schemars::JsonSchema))]#[serde(deny_unknown_fields)]$($t)*};}
dto! {pub struct Fragment {pub source:FragmentSpec,pub treatment:Option<EndTreatment>}}
dto! {pub struct Recipe {pub name:String,pub fragments:Vec<Fragment>,pub circular:bool,pub allow_nicked:bool}}
dto! {pub struct Row {pub preparation:PreparedRow,pub treatment_step:Option<usize>,pub prepared_source:Source}}
dto! {pub struct Compiled {pub plan:DesignPlan,pub rows:Vec<Row>,pub ligation_step:usize}}
dto! {pub struct Junction {pub index:usize,pub upstream:usize,pub downstream:usize,pub closing:bool,pub upstream_end:Option<End>,pub downstream_end:Option<End>}}
dto! {pub struct Preview {pub context_revision:String,pub compiled:Compiled,pub verification:VerificationReport,
/// Leading prepared reference-document prefix, before assembly orientation.
pub prepared:Vec<Molecule>,
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
    preparation::preparation_source_preflight(inputs, target, &recipe.name, &sources(recipe)?)?;
    preparation::hash(recipe)?;
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
) -> Result<String, Error> {
    preflight(inputs, target, recipe)?;
    preparation::bound_revisions(inputs, target)?;
    preparation::hash(&(ENGINE_VERSION, "preview", inputs, target, recipe))
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
        let rows = self.compiled.rows.len();
        let steps = self.compiled.plan.steps.len();
        if !hex(&self.context_revision)
            || rows == 0
            || rows > preparation::MAX_ROWS
            || steps > 4 * preparation::MAX_ROWS + 1
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
        if self.verification.plan_revision != preparation::hash(&self.compiled.plan)?
            || self.verification.target_revision != self.compiled.plan.target_revision
        {
            return Err(bad("verification plan/target binding mismatch"));
        }
        let Some(Operation::Ligate {
            fragments,
            circular,
            ..
        }) = self.compiled.plan.steps.last()
        else {
            return Err(bad("preview has no final ordered ligation operation"));
        };
        if fragments.len() != rows
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
                    )
                })
        {
            return Err(bad("preview operation shape changed"));
        }
        for (index, (row, fragment)) in self.compiled.rows.iter().zip(fragments).enumerate() {
            let p = &row.preparation;
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
