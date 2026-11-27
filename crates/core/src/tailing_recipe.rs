//! Explicit per-row hold/add-tail TA/GC recipes; no implied commercial chemistry.
use crate::{
    Error, Molecule,
    assembly_recipe::{self as preparation, FragmentSpec, PreparedRow},
    chemistry::{Duplex, DuplexSummary, EndTreatment, Polarity},
    tailing::TailingBase,
    verification::{
        self, DesignPlan, Operation, OrientedSource, PlanInput, Source, VerificationReport,
    },
};
use serde::{Deserialize, Serialize};
pub const ENGINE_VERSION: &str = "tailing-recipe-1";
macro_rules! dto{($($t:tt)*)=>{#[derive(Clone,Debug,Serialize,Deserialize)]#[cfg_attr(feature="json-schema",derive(schemars::JsonSchema))]#[serde(deny_unknown_fields)]$($t)*};}
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum Profile {
    Ta,
    Gc,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub enum Action {
    KeepExisting,
    AddTail { base: TailingBase },
}
/// Geometry-family assessment from the fresh native prefix; independent of ligation chemistry/status.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum ProfileStatus {
    Matched,
    Incompatible,
    Undetermined,
}
dto! {pub struct Fragment{pub source:FragmentSpec,pub action:Action,pub treatment:Option<EndTreatment>}}
dto! {pub struct Recipe{pub name:String,pub profile:Profile,pub fragments:Vec<Fragment>,pub circular:bool,pub allow_nicked:bool}}
dto! {pub struct Row{pub preparation:PreparedRow,pub action:Action,pub tail_step:Option<usize>,pub treatment_step:Option<usize>,pub prepared_source:Source}}
dto! {pub struct Compiled{pub plan:DesignPlan,pub rows:Vec<Row>,pub ligation_step:usize}}
pub use crate::ligation_recipe::Junction;
dto! {pub struct Preview{pub context_revision:String,pub recipe:Recipe,pub compiled:Compiled,pub verification:VerificationReport,
pub prepared:Vec<Molecule>,pub physical:Vec<Option<DuplexSummary>>,pub product_physical:Option<DuplexSummary>,
pub junctions:Vec<Junction>,pub preparation_complete:bool,pub profile_status:ProfileStatus,pub profile_findings:Vec<String>,pub warnings:Vec<String>}}
fn bad(s: &str) -> Error {
    Error::Format(format!("tailing recipe: {s}"))
}
fn shape(recipe: &Recipe) -> Result<(), Error> {
    if recipe.fragments.is_empty()
        || recipe.fragments.len() > preparation::MAX_ROWS
        || recipe.name.len() > 1024
    {
        return Err(Error::TooLarge);
    }
    let mut operations = 1usize;
    for fragment in &recipe.fragments {
        operations = operations
            .checked_add(preparation::preparation_operation_count(
                &fragment.source.preparation,
            ))
            .and_then(|n| {
                n.checked_add(
                    usize::from(matches!(fragment.action, Action::AddTail { .. }))
                        + usize::from(fragment.treatment.is_some()),
                )
            })
            .filter(|n| *n <= 1000)
            .ok_or(Error::TooLarge)?;
        if let Action::AddTail { base } = fragment.action {
            let allowed = match recipe.profile {
                Profile::Ta => matches!(base, TailingBase::A | TailingBase::T),
                Profile::Gc => matches!(base, TailingBase::G | TailingBase::C),
            };
            if !allowed {
                return Err(bad(
                    "added base does not belong to the explicitly selected TA/GC profile",
                ));
            }
        }
    }
    // Streamed recipe bound precedes spec/primer/origin cloning. Shared preparation
    // preflight below owns every preparation's more specific shape/model limits.
    preparation::hash(recipe)?;
    crate::run_material::material_capsules_preflight(recipe.fragments.iter().filter_map(|row| {
        match &row.source.preparation {
            preparation::Preparation::FromRun { origin, .. } => {
                Some(crate::run_material::MaterialOriginRef::Run(origin))
            }
            preparation::Preparation::FromDigest { origin, .. } => {
                Some(crate::run_material::MaterialOriginRef::Digest(origin))
            }
            preparation::Preparation::FromDerived { origin, .. } => {
                Some(crate::run_material::MaterialOriginRef::Derived(origin))
            }
            _ => None,
        }
    }))?;
    Ok(())
}
fn sources(recipe: &Recipe) -> Vec<FragmentSpec> {
    recipe.fragments.iter().map(|f| f.source.clone()).collect()
}
pub fn source_preflight(
    inputs: &[(&str, &Molecule)],
    target: &Molecule,
    recipe: &Recipe,
) -> Result<(), Error> {
    shape(recipe)?;
    preparation::preparation_source_preflight(inputs, target, &recipe.name, &sources(recipe))
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
    preparation::bound_revisions(inputs, target)?;
    Ok(())
}
pub fn preview_revision(
    inputs: &[PlanInput],
    target: &Molecule,
    recipe: &Recipe,
) -> Result<String, Error> {
    preflight(inputs, target, recipe)?;
    preparation::hash(&(
        ENGINE_VERSION,
        verification::ENGINE_VERSION,
        inputs,
        target.revision(),
        recipe,
    ))
}
pub fn compile(
    inputs: &[PlanInput],
    target: &Molecule,
    recipe: &Recipe,
) -> Result<Compiled, Error> {
    preflight(inputs, target, recipe)?;
    compile_bound(inputs, &target.revision(), recipe)
}
/// Compile declared operations against an independently pinned target revision.
/// No target DNA or physical replay is needed to check a documentary receipt.
fn compile_bound(
    inputs: &[PlanInput],
    target_revision: &str,
    recipe: &Recipe,
) -> Result<Compiled, Error> {
    shape(recipe)?;
    if target_revision.len() != 64 || !target_revision.bytes().all(|b| b.is_ascii_hexdigit()) {
        return Err(bad("invalid target revision"));
    }
    let specs = sources(recipe);
    preparation::preparation_rows_preflight(
        &inputs
            .iter()
            .map(|i| (i.id.as_str(), &i.molecule))
            .collect::<Vec<_>>(),
        &recipe.name,
        &specs,
    )?;
    for input in inputs {
        if preparation::hash(&input.molecule)? != input.revision {
            return Err(bad("original source revision changed"));
        }
    }
    let (mut steps, prepared_rows) = preparation::preparation_steps(inputs, &specs)?;
    let mut rows = Vec::with_capacity(prepared_rows.len());
    let mut fragments = Vec::with_capacity(prepared_rows.len());
    for (row, request) in prepared_rows.into_iter().zip(&recipe.fragments) {
        let mut source = row.prepared_source.clone();
        let mut tail_step = None;
        let mut treatment_step = None;
        if let Action::AddTail { base } = request.action {
            let index = steps.len();
            steps.push(Operation::Tail3Prime {
                template: source,
                base,
            });
            source = Source::Step(index);
            tail_step = Some(index);
        }
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
            action: request.action,
            tail_step,
            treatment_step,
            prepared_source: source,
        });
    }
    let ligation_step = steps.len();
    steps.push(crate::ligation_outcomes::operation(
        fragments,
        recipe.circular,
        recipe.allow_nicked,
    ));
    let result = Compiled {
        plan: DesignPlan {
            name: recipe.name.clone(),
            target_revision: target_revision.into(),
            inputs: inputs.to_vec(),
            steps,
        },
        rows,
        ligation_step,
    };
    preparation::hash(&result)?;
    Ok(result)
}
fn source_document<'a>(
    plan: &'a DesignPlan,
    report: &'a VerificationReport,
    source: &Source,
) -> Option<&'a Molecule> {
    match source {
        Source::Input(id) => plan
            .inputs
            .iter()
            .find(|i| &i.id == id)
            .map(|i| &i.molecule),
        Source::Step(i) => report.products.get(*i),
    }
}
fn junctions(physical: &[Option<DuplexSummary>], rows: usize, circular: bool) -> Vec<Junction> {
    (0..rows.saturating_sub(1) + usize::from(circular))
        .map(|index| {
            let next = (index + 1) % rows;
            Junction {
                index,
                upstream: index,
                downstream: next,
                closing: next == 0,
                upstream_end: physical
                    .get(index)
                    .and_then(Option::as_ref)
                    .and_then(|s| s.ends.as_ref())
                    .map(|(_, r)| r.clone()),
                downstream_end: physical
                    .get(next)
                    .and_then(Option::as_ref)
                    .and_then(|s| s.ends.as_ref())
                    .map(|(l, _)| l.clone()),
            }
        })
        .collect()
}
// Preview invokes this only on summaries just produced from native duplex values.
// Readonly DTO validation can check documentary correspondence, never recreate authority.
fn profile_assessment(
    physical: &[Option<DuplexSummary>],
    rows: usize,
    profile: Profile,
    circular: bool,
) -> (ProfileStatus, Vec<String>) {
    let mut findings = Vec::new();
    let mut incompatible = false;
    let mut unknown = physical.len() != rows;
    for row in 0..rows {
        let Some(s) = physical.get(row).and_then(Option::as_ref) else {
            unknown = true;
            findings.push(format!(
                "Row {} has no fresh native end geometry; the TA/GC profile is not established.",
                row + 1
            ));
            continue;
        };
        let Some((left, right)) = &s.ends else {
            incompatible = true;
            findings.push(format!(
                "Row {} has no accessible linear ends; it does not match this tailing profile.",
                row + 1
            ));
            continue;
        };
        for (side, end) in [("left", left), ("right", right)] {
            let family = match profile {
                Profile::Ta => matches!(end.overhang.as_str(), "A" | "T"),
                Profile::Gc => matches!(end.overhang.as_str(), "G" | "C"),
            };
            if end.polarity != Polarity::ThreePrime || !family {
                incompatible = true;
                findings.push(format!("Row {} {side} end is not an actual single-base 3′ overhang in the declared TA/GC family.",row+1));
            }
        }
    }
    for junction in junctions(physical, rows, circular) {
        if let (Some(a), Some(b)) = (junction.upstream_end, junction.downstream_end) {
            let complementary = matches!(
                (a.overhang.as_str(), b.overhang.as_str()),
                ("A", "T") | ("T", "A") | ("G", "C") | ("C", "G")
            );
            if !complementary {
                incompatible = true;
                findings.push(format!("Declared junction {} lacks complementary one-base tails in the fresh prepared geometry.",junction.index+1));
            }
        }
    }
    if incompatible {
        (ProfileStatus::Incompatible, findings)
    } else if unknown {
        (ProfileStatus::Undetermined, findings)
    } else {
        findings.push("Fresh prepared end geometry matches the declared TA/GC family and intended junction complements. Phosphate/OH eligibility, nicking and product compatibility remain separate native verifier conclusions.".into());
        (ProfileStatus::Matched, findings)
    }
}
pub fn preview(inputs: &[PlanInput], target: &Molecule, recipe: &Recipe) -> Result<Preview, Error> {
    let context_revision = preview_revision(inputs, target, recipe)?;
    let compiled = compile(inputs, target, recipe)?;
    let (verification, trace) = verification::verify_with_physical(&compiled.plan, target);
    let mut prepared = Vec::new();
    let mut physical = Vec::new();
    let mut bases = 0usize;
    for row in &compiled.rows {
        let Some(molecule) = source_document(&compiled.plan, &verification, &row.prepared_source)
        else {
            break;
        };
        bases = bases
            .checked_add(molecule.sequence.len())
            .filter(|n| *n <= preparation::MAX_ROLE_BASES)
            .ok_or(Error::TooLarge)?;
        crate::model::bounds::model_preflight(&[molecule])?;
        let raw = match &row.prepared_source {
            Source::Step(i) => trace.get(*i).and_then(Option::as_ref).cloned(),
            Source::Input(_) => Duplex::from_molecule(molecule).ok(),
        };
        let summary = raw
            .map(|d| {
                if row.preparation.reverse {
                    d.reverse()?.summary()
                } else {
                    d.summary()
                }
            })
            .transpose()?;
        prepared.push(molecule.clone());
        physical.push(summary);
    }
    let preparation_complete = prepared.len() == compiled.rows.len();
    let product_physical = trace
        .get(compiled.ligation_step)
        .and_then(Option::as_ref)
        .map(Duplex::summary)
        .transpose()?;
    let junctions = junctions(&physical, compiled.rows.len(), recipe.circular);
    let (profile_status, profile_findings) = profile_assessment(
        &physical,
        compiled.rows.len(),
        recipe.profile,
        recipe.circular,
    );
    let result=Preview{context_revision,recipe:recipe.clone(),compiled,verification,prepared,physical,product_physical,junctions,preparation_complete,profile_status,profile_findings,
        warnings:vec!["KeepExisting adds no tail and makes no assertion that the existing ends are complementary, phosphorylated, intact or ligatable. Only fresh native source reconstruction and ligation assess those states.".into(),
        "TA/GC profile geometry requires fresh one-base 3′ tails of the declared family on every prepared row and complementary intended junctions; KeepExisting alone proves none of these. Profiles do not identify an enzyme, commercial kit, reaction yield, isolated clone or exclusive product; retained alternatives remain within their declared one-copy search scope.".into()]};
    preparation::hash(&result)?;
    Ok(result)
}
impl Preview {
    /// Readonly declaration/source correspondence, never chemical replay or authority.
    pub fn validate(&self) -> Result<(), Error> {
        self.validate_binding(&self.context_revision)
    }
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        let recipe = &self.recipe;
        let inputs = &self.compiled.plan.inputs;
        let target_revision = &self.compiled.plan.target_revision;
        let context = preparation::hash(&(
            ENGINE_VERSION,
            verification::ENGINE_VERSION,
            inputs,
            target_revision,
            recipe,
        ))?;
        if self.context_revision != expected
            || context != expected
            || preparation::hash(&self.compiled)?
                != preparation::hash(&compile_bound(inputs, target_revision, recipe)?)?
            || self.verification.engine_version != verification::ENGINE_VERSION
            || self.verification.plan_revision != verification::plan_revision(&self.compiled.plan)
            || self.verification.target_revision != *target_revision
            || self.verification.products.len() > self.compiled.plan.steps.len()
            || self.prepared.len() > self.compiled.rows.len()
            || self.physical.len() != self.prepared.len()
            || self.preparation_complete != (self.prepared.len() == self.compiled.rows.len())
        {
            return Err(bad(
                "preview source, recipe, compiler or verifier binding differs",
            ));
        }
        let prefix = self
            .compiled
            .rows
            .iter()
            .take_while(|row| {
                source_document(
                    &self.compiled.plan,
                    &self.verification,
                    &row.prepared_source,
                )
                .is_some()
            })
            .count();
        if self.prepared.len() != prefix {
            return Err(bad("prepared prefix was truncated or extended"));
        }
        for (row, molecule) in self.compiled.rows.iter().zip(&self.prepared) {
            let expected = source_document(
                &self.compiled.plan,
                &self.verification,
                &row.prepared_source,
            )
            .ok_or_else(|| bad("prepared source missing"))?;
            if preparation::hash(molecule)? != preparation::hash(expected)? {
                return Err(bad("prepared source snapshot differs"));
            }
        }
        for s in self
            .physical
            .iter()
            .filter_map(Option::as_ref)
            .chain(self.product_physical.iter())
        {
            if s.model.len() > 256
                || s.top_length > preparation::MAX_ROLE_BASES
                || s.bottom_length > preparation::MAX_ROLE_BASES
                || s.bottom_start.unsigned_abs() > preparation::MAX_ROLE_BASES as u64
                || s.nicks.len() > 10_000
                || s.ends.as_ref().is_some_and(|(a, b)| {
                    a.overhang.len() > preparation::MAX_ROLE_BASES
                        || b.overhang.len() > preparation::MAX_ROLE_BASES
                })
            {
                return Err(Error::TooLarge);
            }
        }
        if preparation::hash(&self.junctions)?
            != preparation::hash(&junctions(
                &self.physical,
                self.compiled.rows.len(),
                recipe.circular,
            ))?
        {
            return Err(bad("junction correspondence differs"));
        }
        let (profile_status, profile_findings) = profile_assessment(
            &self.physical,
            self.compiled.rows.len(),
            recipe.profile,
            recipe.circular,
        );
        if self.profile_status != profile_status || self.profile_findings != profile_findings {
            return Err(bad("profile geometry correspondence differs"));
        }
        preparation::hash(self)?;
        Ok(())
    }
    /// Additionally binds a receipt to the caller's actual independently declared request.
    pub fn validate_for(
        &self,
        inputs: &[PlanInput],
        target: &Molecule,
        recipe: &Recipe,
        expected: &str,
    ) -> Result<(), Error> {
        let context = preview_revision(inputs, target, recipe)?;
        if context != expected
            || preparation::hash(&inputs)? != preparation::hash(&self.compiled.plan.inputs)?
            || target.revision() != self.compiled.plan.target_revision
            || preparation::hash(recipe)? != preparation::hash(&self.recipe)?
        {
            return Err(bad("preview belongs to another source/target request"));
        }
        self.validate_binding(expected)
    }
}
