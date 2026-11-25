//! Explicit whole-template, whole-amplicon and isolated-substrate reaction inputs.
//! The final native enumeration is captured once, including failed probe selection.
use crate::{
    Alphabet, Error, Molecule,
    assembly_recipe::{self, FragmentSpec, Preparation},
    chemistry::EndTreatment,
    golden_gate::{self, GoldenGateOptions, GoldenGateReport},
    prepared_golden_gate::{Assessment, Selection, SelectionReport},
    verification::{self, DesignPlan, Operation, PlanInput, Source, VerificationReport},
};
use serde::{Deserialize, Serialize};
use std::collections::BTreeSet;
pub const ENGINE_VERSION: &str = "golden-gate-reaction-1";
macro_rules! dto{($($t:tt)*)=>{#[derive(Clone,Debug,Serialize,Deserialize)]#[cfg_attr(feature="json-schema",derive(schemars::JsonSchema))]#[serde(deny_unknown_fields)]$($t)*};}
dto! {pub enum Substrate{WholeSource,Prepared{preparation:Preparation,treatment:Option<EndTreatment>}}}
dto! {pub struct Row{pub role:String,pub input_id:String,pub substrate:Substrate}}
dto! {pub struct Request{pub name:String,pub rows:Vec<Row>,pub options:GoldenGateOptions}}
dto! {pub struct Material{pub index:usize,pub role:String,pub input_id:String,pub source:Source}}
dto! {pub struct Report{pub context_revision:String,pub request:Request,pub target:Molecule,pub probe_plan:DesignPlan,pub materials:Vec<Material>,pub verification:VerificationReport,#[cfg_attr(feature="json-schema",schemars(with="serde_json::Value"))]pub pool:Option<GoldenGateReport>,pub assessments:Vec<Assessment>,pub warnings:Vec<String>}}
fn bad(s: &str) -> Error {
    Error::Unsupported(format!("Golden Gate reaction: {s}"))
}
fn token(s: &str, n: usize) -> bool {
    !s.is_empty() && s.len() <= n && !s.chars().any(char::is_control)
}
fn options(o: &GoldenGateOptions) -> Result<(), Error> {
    if o.min_fragments == 0
        || o.min_fragments > o.max_fragments
        || o.max_fragments > 32
        || o.max_products == 0
        || o.max_products > 1000
        || o.max_states == 0
        || o.max_states > 1_000_000
    {
        return Err(bad("invalid finite fragment/product/state limits"));
    }
    Ok(())
}
fn specs(r: &Request) -> Vec<FragmentSpec> {
    r.rows
        .iter()
        .filter_map(|row| match &row.substrate {
            Substrate::WholeSource => None,
            Substrate::Prepared { preparation, .. } => Some(FragmentSpec {
                role: row.role.clone(),
                input_id: row.input_id.clone(),
                preparation: preparation.clone(),
                reverse: false,
            }),
        })
        .collect()
}
pub fn preflight(
    inputs: &[(&str, &Molecule)],
    target: &Molecule,
    r: &Request,
) -> Result<(), Error> {
    options(&r.options)?;
    if inputs.is_empty()
        || inputs.len() > 32
        || r.rows.is_empty()
        || r.rows.len() > 32
        || !token(&r.name, 1024)
    {
        return Err(bad("choose1..32 bound inputs and explicit reaction roles"));
    }
    let mut models = inputs.iter().map(|(_, m)| *m).collect::<Vec<_>>();
    models.push(target);
    crate::model::bounds::model_preflight(&models)?;
    let mut bases = 0usize;
    let mut ids = BTreeSet::new();
    for (id, m) in inputs {
        if !token(id, 128) || !ids.insert(*id) {
            return Err(bad("invalid/duplicate original input ID"));
        }
        m.validate()?;
        if m.sequence.alphabet() != Alphabet::Dna {
            return Err(bad("reaction inputs must be DNA"));
        }
        bases = bases.checked_add(m.sequence.len()).ok_or(Error::TooLarge)?;
    }
    target.validate()?;
    if target.sequence.alphabet() != Alphabet::Dna
        || bases
            .checked_add(target.sequence.len())
            .is_none_or(|n| n > 8 << 20)
    {
        return Err(Error::TooLarge);
    }
    let mut roles = BTreeSet::new();
    for row in &r.rows {
        if !token(&row.role, 128)
            || !roles.insert(&row.role)
            || !ids.contains(row.input_id.as_str())
        {
            return Err(bad("invalid/duplicate role or missing original input"));
        }
    }
    assembly_recipe::hash(r)?;
    let prepared = specs(r);
    if !prepared.is_empty() {
        assembly_recipe::preparation_source_preflight(inputs, target, &r.name, &prepared)?;
    }
    assembly_recipe::hash(&(inputs, target, r))?;
    Ok(())
}
pub fn inspect_revision(
    inputs: &[PlanInput],
    target: &Molecule,
    r: &Request,
) -> Result<String, Error> {
    preflight(
        &inputs
            .iter()
            .map(|i| (i.id.as_str(), &i.molecule))
            .collect::<Vec<_>>(),
        target,
        r,
    )?;
    for i in inputs {
        if i.revision != i.molecule.revision() {
            return Err(bad("original source revision changed"));
        }
    }
    assembly_recipe::hash(&(
        ENGINE_VERSION,
        inputs
            .iter()
            .map(|i| (&i.id, &i.revision))
            .collect::<Vec<_>>(),
        target.revision(),
        r,
    ))
}
fn compile(
    inputs: &[PlanInput],
    target: &Molecule,
    r: &Request,
) -> Result<(DesignPlan, Vec<Material>), Error> {
    let prepared = specs(r);
    let (mut steps, prepared_rows) = if prepared.is_empty() {
        (vec![], vec![])
    } else {
        let c = assembly_recipe::compile_preparations(inputs, target, &r.name, &prepared)?;
        (c.steps, c.rows)
    };
    let mut materials = Vec::with_capacity(r.rows.len());
    for (index, row) in r.rows.iter().enumerate() {
        let mut source = match &row.substrate {
            Substrate::WholeSource => Source::Input(row.input_id.clone()),
            Substrate::Prepared { .. } => prepared_rows
                .iter()
                .find(|p| p.role == row.role)
                .ok_or_else(|| bad("compiled preparation role missing"))?
                .prepared_source
                .clone(),
        };
        if let Substrate::Prepared {
            treatment: Some(treatment),
            ..
        } = &row.substrate
        {
            steps.push(Operation::TreatEnds {
                template: source,
                treatment: *treatment,
            });
            source = Source::Step(steps.len() - 1);
        }
        materials.push(Material {
            index,
            role: row.role.clone(),
            input_id: row.input_id.clone(),
            source,
        });
    }
    steps.push(Operation::GoldenGate {
        templates: materials.iter().map(|m| m.source.clone()).collect(),
        options: r.options.clone(),
        candidate_index: 0,
    });
    Ok((
        DesignPlan {
            name: r.name.clone(),
            inputs: inputs.to_vec(),
            target_revision: target.revision(),
            steps,
        },
        materials,
    ))
}
fn matches(
    plan: &DesignPlan,
    duplex: &crate::chemistry::Duplex,
    target: &Molecule,
) -> Result<bool, Error> {
    let m = Molecule::dna("candidate comparison", duplex.top.as_str(), duplex.topology)?;
    Ok(verification::plan_target_matches(plan, &m, target))
}
fn assessments(
    plan: &DesignPlan,
    target: &Molecule,
    pool: &GoldenGateReport,
) -> Result<Vec<Assessment>, Error> {
    pool.candidates
        .iter()
        .enumerate()
        .map(|(candidate_index, c)| {
            let mut roles = BTreeSet::new();
            for p in &c.path {
                roles.insert(
                    pool.fragments
                        .get(p.fragment)
                        .ok_or_else(|| bad("candidate path fragment missing"))?
                        .source_index,
                );
            }
            Ok(Assessment {
                candidate_index,
                directed_sequence_key: c.directed_sequence_key.clone(),
                physical_sequence_key: c.physical_sequence_key.clone(),
                target_matches: matches(plan, &c.duplex, target)?,
                source_roles: roles.into_iter().collect(),
            })
        })
        .collect()
}
pub fn inspect(inputs: &[PlanInput], target: &Molecule, r: &Request) -> Result<Report, Error> {
    let context_revision = inspect_revision(inputs, target, r)?;
    let (probe_plan, materials) = compile(inputs, target, r)?;
    let (verification, _, pool) = verification::verify_with_golden_gate_trace(&probe_plan, target);
    let assessments = pool
        .as_ref()
        .map(|p| assessments(&probe_plan, target, p))
        .transpose()?
        .unwrap_or_default();
    let report=Report{context_revision:context_revision.clone(),request:r.clone(),target:target.clone(),probe_plan,materials,verification,pool,assessments,warnings:vec![
        "Each declared role contributes one actual reaction substrate. WholeSource includes the complete original circular/linear DNA; Prepared ExactPcr includes the complete amplicon, without discarding flanking digest fragments. Other explicit preparations may intentionally select isolated material.".into(),
        "The probe plan uses index0 only to collect preparation and enumeration diagnostics. Its target assessment is not an automatic choice or saved candidate. Select an exact path and freshly verify the selected plan.".into(),
        "Complete cleavage and the finite one-copy/orientation limits are assumptions, not kinetic cycling, repeated-copy stoichiometry, methylation clearance, real-mixture exhaustiveness or experimental yield. Carryover is reported separately from newly assembled products.".into()]};
    report.validate_binding(&context_revision)?;
    Ok(report)
}
impl Report {
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        let actual = inspect_revision(&self.probe_plan.inputs, &self.target, &self.request)?;
        if expected != actual || self.context_revision != actual {
            return Err(bad("reaction context changed"));
        }
        let (plan, materials) = compile(&self.probe_plan.inputs, &self.target, &self.request)?;
        if assembly_recipe::hash(&plan)? != assembly_recipe::hash(&self.probe_plan)?
            || assembly_recipe::hash(&materials)? != assembly_recipe::hash(&self.materials)?
            || self.verification.plan_revision != verification::plan_revision(&plan)
            || self.verification.target_revision != self.target.revision()
            || self.verification.engine_version != verification::ENGINE_VERSION
            || self.verification.products.len() > plan.steps.len()
        {
            return Err(bad("reaction plan/preparation receipt changed"));
        }
        for product in &self.verification.products {
            crate::model::bounds::model_preflight(&[product])?;
            product.validate()?;
            if product.metadata.get("plan_revision") != Some(&self.verification.plan_revision) {
                return Err(bad("preparation product belongs to another plan"));
            }
        }
        if let Some(pool) = &self.pool {
            if pool.model != golden_gate::MODEL
                || assembly_recipe::hash(&pool.options)?
                    != assembly_recipe::hash(&self.request.options)?
                || pool.sources.len() != materials.len()
                || pool.fragments.len() > 128
                || pool.candidates.len() > self.request.options.max_products
                || pool.states_examined > self.request.options.max_states + 1
                || pool.exhaustive_mixture_search
                || self.verification.products.len() + 1 < plan.steps.len()
            {
                return Err(bad("reaction enumeration scope changed"));
            }
            for (index, (material, s)) in materials.iter().zip(&pool.sources).enumerate() {
                let actual = match &material.source {
                    Source::Input(id) => {
                        &plan
                            .inputs
                            .iter()
                            .find(|i| &i.id == id)
                            .ok_or_else(|| bad("reaction original missing"))?
                            .molecule
                    }
                    Source::Step(i) => self
                        .verification
                        .products
                        .get(*i)
                        .ok_or_else(|| bad("actual preparation missing"))?,
                };
                if s.source_index != index || s.revision != actual.revision() {
                    return Err(bad(
                        "reaction substrate revision differs from its actual source",
                    ));
                }
            }
            for f in &pool.fragments {
                f.duplex.validate()?;
                if f.source_index >= materials.len()
                    || f.source_revision != pool.sources[f.source_index].revision
                {
                    return Err(bad("reaction fragment source changed"));
                }
            }
            let mut stable = BTreeSet::new();
            for c in &pool.candidates {
                c.duplex.validate()?;
                if c.restriction_site_free != c.recut_sites.is_empty()
                    || (c.duplex.topology == crate::Topology::Circular)
                        != self.request.options.circular
                {
                    return Err(bad("reaction candidate topology or recut scope changed"));
                }
                if c.restriction_site_free {
                    stable.insert(&c.physical_sequence_key);
                }
                let (directed, physical) = golden_gate::sequence_keys(&c.duplex)?;
                if c.directed_sequence_key != directed
                    || c.physical_sequence_key != physical
                    || c.path.len() < self.request.options.min_fragments
                    || c.path.len() > self.request.options.max_fragments
                {
                    return Err(bad("reaction candidate identity or path scope changed"));
                }
                let mut used = BTreeSet::new();
                for p in &c.path {
                    if p.fragment >= pool.fragments.len()
                        || !used.insert(p.fragment)
                        || p.reverse && !self.request.options.both_orientations
                    {
                        return Err(bad("reaction path violates one-copy/orientation scope"));
                    }
                }
            }
            if stable.len() != pool.stable_physical_product_count
                || pool.carryover_sources.iter().any(|i| *i >= materials.len())
                || pool.enumeration_complete && !pool.truncation_reasons.is_empty()
            {
                return Err(bad("reaction completeness or carryover scope changed"));
            }
            if assembly_recipe::hash(&self.assessments)?
                != assembly_recipe::hash(&assessments(&plan, &self.target, pool)?)?
            {
                return Err(bad("reaction candidate assessments changed"));
            }
        } else if !self.assessments.is_empty() {
            return Err(bad("withheld reaction cannot retain candidate assessments"));
        }
        if self.warnings.len() > 64 || self.warnings.iter().any(|w| w.len() > 8192) {
            return Err(Error::TooLarge);
        }
        assembly_recipe::hash(self)?;
        Ok(())
    }
}
fn selected_plan(report: &Report, s: &Selection) -> Result<DesignPlan, Error> {
    report.validate_binding(&s.context_revision)?;
    let c = report
        .pool
        .as_ref()
        .and_then(|p| p.candidates.get(s.candidate_index))
        .ok_or_else(|| bad("reaction candidate unavailable"))?;
    if c.directed_sequence_key != s.directed_sequence_key
        || c.physical_sequence_key != s.physical_sequence_key
        || !c.restriction_site_free
    {
        return Err(bad("candidate keys changed or recognition sites remain"));
    }
    let mut plan = report.probe_plan.clone();
    let Some(Operation::GoldenGate {
        candidate_index, ..
    }) = plan.steps.last_mut()
    else {
        return Err(bad("final reaction operation missing"));
    };
    *candidate_index = s.candidate_index;
    Ok(plan)
}
pub fn selection_revision(report: &Report, s: &Selection) -> Result<String, Error> {
    assembly_recipe::hash(&(
        "prepared-golden-gate-selection-1",
        &report.context_revision,
        s,
        selected_plan(report, s)?,
    ))
}
pub fn select(report: &Report, s: &Selection) -> Result<SelectionReport, Error> {
    let plan = selected_plan(report, s)?;
    let context_revision = assembly_recipe::hash(&(
        "prepared-golden-gate-selection-1",
        &report.context_revision,
        s,
        &plan,
    ))?;
    let selected=SelectionReport{context_revision:context_revision.clone(),parent_context:report.context_revision.clone(),selection:s.clone(),plan,warnings:vec!["Unverified whole-reaction candidate. All original source and explicit preparation steps remain bound. Fresh native verification must reconstruct this exact candidate; target correspondence is not a laboratory or uniqueness verdict.".into()]};
    selected.validate_binding(&context_revision)?;
    Ok(selected)
}
