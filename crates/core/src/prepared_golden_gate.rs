//! Finite Golden Gate alternatives from the same fresh preparation verifier pass.
//! Serialized summaries and pool duplexes never restore source physical authority.
use crate::{
    Error, Molecule, Topology, assembly_recipe,
    golden_gate::{self, GoldenGateOptions, GoldenGateReport},
    ligation_recipe,
    verification::{self, DesignPlan, Operation, PlanInput, Source},
};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};
pub const ENGINE_VERSION: &str = "prepared-golden-gate-1";
macro_rules! dto{($($t:tt)*)=>{#[derive(Clone,Debug,Serialize,Deserialize)]#[cfg_attr(feature="json-schema",derive(schemars::JsonSchema))]#[serde(deny_unknown_fields)]$($t)*};}
dto! {pub struct Request {pub recipe:ligation_recipe::Recipe,pub options:GoldenGateOptions}}
dto! {pub struct Assessment {pub candidate_index:usize,pub directed_sequence_key:String,pub physical_sequence_key:String,pub target_matches:bool,pub source_roles:Vec<usize>}}
dto! {pub struct Report {pub context_revision:String,pub request:Request,pub target:Molecule,pub preview:ligation_recipe::Preview,pub missing_physical:Vec<usize>,#[cfg_attr(feature="json-schema",schemars(with="serde_json::Value"))] pub pool:Option<GoldenGateReport>,pub assessments:Vec<Assessment>,pub warnings:Vec<String>}}
dto! {pub struct Selection {pub context_revision:String,pub candidate_index:usize,pub directed_sequence_key:String,pub physical_sequence_key:String}}
dto! {pub struct SelectionReport {pub context_revision:String,pub parent_context:String,pub selection:Selection,pub plan:DesignPlan,pub warnings:Vec<String>}}
fn bad(s: &str) -> Error {
    Error::Unsupported(format!("prepared Golden Gate: {s}"))
}
fn hex(s: &str) -> bool {
    s.len() == 64
        && s.bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
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
        return Err(bad(
            "finite pool requires fragments1..32, products1..1000, states1..1000000",
        ));
    }
    Ok(())
}
pub fn material_origins(
    r: &Request,
) -> Vec<(&str, &str, crate::run_material::MaterialOriginRef<'_>)> {
    r.recipe
        .fragments
        .iter()
        .filter_map(|f| {
            let origin = match &f.source.preparation {
                assembly_recipe::Preparation::FromRun { origin, .. } => {
                    crate::run_material::MaterialOriginRef::Run(origin)
                }
                assembly_recipe::Preparation::FromDerived { origin, .. } => {
                    crate::run_material::MaterialOriginRef::Derived(origin)
                }
                _ => return None,
            };
            Some((f.source.role.as_str(), f.source.input_id.as_str(), origin))
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
        || r.recipe.fragments.is_empty()
        || r.recipe.fragments.len() > 32
    {
        return Err(Error::TooLarge);
    }
    ligation_recipe::source_preflight(inputs, target, &r.recipe)?;
    assembly_recipe::hash(&(inputs, target, r))?;
    Ok(())
}
pub fn source_preflight(
    inputs: &[(&str, &Molecule)],
    target: &Molecule,
    r: &Request,
) -> Result<(), Error> {
    preflight(inputs, target, r)
}
pub fn inspect_revision(
    inputs: &[PlanInput],
    target: &Molecule,
    r: &Request,
) -> Result<String, Error> {
    if inputs.len() > 32 {
        return Err(Error::TooLarge);
    }
    preflight(
        &inputs
            .iter()
            .map(|i| (i.id.as_str(), &i.molecule))
            .collect::<Vec<_>>(),
        target,
        r,
    )?;
    for i in inputs {
        if i.revision != assembly_recipe::hash(&i.molecule)? {
            return Err(bad("original input revision changed"));
        }
    }
    assembly_recipe::hash(&(
        ENGINE_VERSION,
        inputs
            .iter()
            .map(|i| (&i.id, &i.revision))
            .collect::<Vec<_>>(),
        assembly_recipe::hash(target)?,
        r,
    ))
}
fn roles(pool: &GoldenGateReport, candidate: usize) -> Result<Vec<usize>, Error> {
    let mut result = Vec::new();
    let mut seen = BTreeSet::new();
    for part in &pool
        .candidates
        .get(candidate)
        .ok_or_else(|| bad("candidate is unavailable"))?
        .path
    {
        let source = pool
            .fragments
            .get(part.fragment)
            .ok_or_else(|| bad("pool fragment index is unavailable"))?
            .source_index;
        if seen.insert(source) {
            result.push(source);
        }
    }
    Ok(result)
}
fn target_matches(plan: &DesignPlan, duplex: &crate::chemistry::Duplex, target: &Molecule) -> bool {
    // Readonly reference view for the existing coordinate comparison only. This is
    // never fed back as a physical substrate, preparation, target or simulated input.
    let reference = Molecule {
        name: "Candidate reference comparison".into(),
        sequence: duplex.top.clone(),
        topology: duplex.topology,
        features: vec![],
        metadata: BTreeMap::new(),
        original_genbank: None,
    };
    verification::plan_target_matches(plan, &reference, target)
}
pub fn inspect(inputs: &[PlanInput], target: &Molecule, r: &Request) -> Result<Report, Error> {
    let context_revision = inspect_revision(inputs, target, r)?;
    let mut budget = crate::modified_duplex::CleavageBudget::new(100_000_000)?;
    let (preview, physical) = ligation_recipe::preview_with_tracked(
        inputs,
        target,
        &r.recipe,
        assembly_recipe::CleavagePolicy::PositionalChecked,
        &mut budget,
    )?;
    let missing_physical = (0..r.recipe.fragments.len())
        .filter(|i| *i >= preview.prepared.len() || physical.get(*i).is_none_or(Option::is_none))
        .collect::<Vec<_>>();
    let pool = if preview.preparation_complete && missing_physical.is_empty() {
        let substrates = physical
            .into_iter()
            .map(|d| d.ok_or_else(|| bad("fresh prepared physical state unavailable")))
            .collect::<Result<Vec<_>, _>>()?;
        Some(
            golden_gate::simulate_tracked(
                &preview.prepared,
                &substrates,
                &r.options,
                None,
                &mut budget,
            )?
            .report,
        )
    } else {
        None
    };
    let mut assessments = Vec::new();
    if let Some(pool) = &pool {
        for (i, c) in pool.candidates.iter().enumerate() {
            assessments.push(Assessment {
                candidate_index: i,
                directed_sequence_key: c.directed_sequence_key.clone(),
                physical_sequence_key: c.physical_sequence_key.clone(),
                target_matches: target_matches(&preview.compiled.plan, &c.duplex, target),
                source_roles: roles(pool, i)?,
            });
        }
    }
    let mut warnings=vec!["Finite prepared-pool alternatives use actual private duplexes from the same fresh preparation pass. Ordered ligation and independent target correspondence do not establish a unique mixture outcome or experimental success.".into(),"both_orientations=false refers to the raw stored prepared substrates, not the ordered recipe's reverse flags. All source roles and all reported candidate paths, including subsets, self-products and recuttable products, remain visible.".into(),"Complete cleavage, tracked modification uncertainty and every originating verifier limitation still apply. No serialized endpoint summary restores physical substrate state; selection produces an unverified native operation requiring fresh verification.".into()];
    warnings.push("Pool scope is the explicitly selected prepared roles after any requested digest isolation, repair and treatment. Discarded digest fragments and uncut original-source/amplicon material are not automatically included. This is not a complete one-pot original reaction/carryover model; reported carryover refers only to these prepared substrates.".into());
    if pool.is_none() {
        warnings.push("Pool withheld: one or more requested preparations or actual physical duplex states are unavailable. No available subset was silently simulated as the complete requested pool.".into());
    }
    let result = Report {
        context_revision: context_revision.clone(),
        request: r.clone(),
        target: target.clone(),
        preview,
        missing_physical,
        pool,
        assessments,
        warnings,
    };
    result.validate_binding(&context_revision)?;
    Ok(result)
}
fn texts(v: &[String]) -> Result<(), Error> {
    if v.len() > 10_000 || v.iter().any(|s| s.len() > 8192) {
        return Err(Error::TooLarge);
    }
    Ok(())
}
fn pool_guard(report: &Report, pool: &GoldenGateReport) -> Result<(), Error> {
    if pool.profile_revision != crate::modified_duplex::profile_revision()
        || pool.comparisons > 100_000_000
    {
        return Err(bad("checked cleavage profile/budget differs"));
    }
    for source in &pool.source_cleavage {
        source.validate()?;
    }
    for candidate in &pool.candidates {
        candidate.cleavage.validate()?;
    }
    let n = report.preview.prepared.len();
    if pool.options != report.request.options
        || pool.model != golden_gate::MODEL
        || pool.sources.len() != n
        || pool.fragments.len() > 128
        || pool.candidates.len() > pool.options.max_products
        || pool.states_examined > pool.options.max_states
        || pool.exhaustive_mixture_search
        || pool.enumeration_complete != pool.truncation_reasons.is_empty()
    {
        return Err(bad("invalid finite pool report shape or options"));
    }
    texts(&pool.warnings)?;
    texts(&pool.truncation_reasons)?;
    let motif = pool.options.enzyme.profile().recognition.len();
    for (i, s) in pool.sources.iter().enumerate() {
        let actual = &report.preview.prepared[i];
        if s.source_index != i
            || s.revision != assembly_recipe::hash(actual)?
            || s.name != actual.name.chars().take(1024).collect::<String>()
            || s.recognition_sites.len() > 1_000_000
            || s.recognition_sites
                .iter()
                .any(|m| m.start >= actual.sequence.len() || m.length != motif)
        {
            return Err(bad(
                "pool sources differ from actual ordered prepared references",
            ));
        }
    }
    let mut material = 0usize;
    for f in &pool.fragments {
        if f.source_index >= n
            || f.fragment_index >= 128
            || f.source_revision != pool.sources[f.source_index].revision
            || f.top_start > report.preview.prepared[f.source_index].sequence.len()
        {
            return Err(bad("invalid pool source/fragment binding"));
        }
        f.duplex.validate()?;
        material = material
            .checked_add(f.duplex.top.len())
            .and_then(|n| n.checked_add(f.duplex.bottom.len()))
            .filter(|n| *n <= 16 << 20)
            .ok_or(Error::TooLarge)?;
        match f.duplex.topology {
            Topology::Linear => {
                let (left, right) = f.duplex.ends()?;
                if f.left_end.as_ref() != Some(&left) || f.right_end.as_ref() != Some(&right) {
                    return Err(bad(
                        "fragment endpoint summaries disagree with retained duplex geometry",
                    ));
                }
            }
            Topology::Circular => {
                if f.left_end.is_some() || f.right_end.is_some() {
                    return Err(bad(
                        "uncut circular carryover has no linear endpoint summaries",
                    ));
                }
            }
        }
    }
    let mut carry = BTreeSet::new();
    for i in &pool.carryover_sources {
        if *i >= n || !carry.insert(*i) {
            return Err(bad("invalid carryover source index"));
        }
    }
    let mut stable = BTreeSet::new();
    for c in &pool.candidates {
        if c.path.len() < pool.options.min_fragments
            || c.path.len() > pool.options.max_fragments
            || c.duplex.topology
                != if pool.options.circular {
                    Topology::Circular
                } else {
                    Topology::Linear
                }
            || c.recut_sites.len() > 1_000_000
            || c.restriction_site_free != (c.cleavage.complete && c.cleavage.sites.is_empty())
            || !c.duplex.nicks.is_empty()
        {
            return Err(bad("invalid candidate path/topology/recut shape"));
        }
        let mut used = BTreeSet::new();
        for part in &c.path {
            if part.fragment >= pool.fragments.len()
                || !used.insert(part.fragment)
                || part.reverse && !pool.options.both_orientations
            {
                return Err(bad("candidate violates one-copy/raw-orientation scope"));
            }
        }
        c.duplex.validate()?;
        material = material
            .checked_add(c.duplex.top.len())
            .and_then(|n| n.checked_add(c.duplex.bottom.len()))
            .filter(|n| *n <= 16 << 20)
            .ok_or(Error::TooLarge)?;
        let (directed, physical) = golden_gate::sequence_keys(&c.duplex)?;
        if c.directed_sequence_key != directed || c.physical_sequence_key != physical {
            return Err(bad(
                "candidate sequence keys disagree with its actual retained reference",
            ));
        }
        for s in &c.recut_sites {
            if s.recognition_start >= c.duplex.top.len()
                || s.top_cut.unsigned_abs() > c.duplex.top.len() as u64 + 16
                || s.bottom_cut.unsigned_abs() > c.duplex.top.len() as u64 + 16
            {
                return Err(bad("invalid recut coordinates"));
            }
        }
        if c.restriction_site_free {
            stable.insert(&c.physical_sequence_key);
        }
    }
    if pool.stable_physical_product_count != stable.len() {
        return Err(bad("physical product group count changed"));
    }
    Ok(())
}
impl Report {
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        if !hex(expected) || self.context_revision != expected {
            return Err(bad("inspection context changed"));
        }
        crate::model::bounds::model_preflight(&[&self.target])?;
        self.target.validate()?;
        let inputs = &self.preview.compiled.plan.inputs;
        if self.preview.verification.target_revision != assembly_recipe::hash(&self.target)?
            || inspect_revision(inputs, &self.target, &self.request)? != self.context_revision
        {
            return Err(bad("actual source/target/request binding changed"));
        }
        self.preview
            .validate_binding(&ligation_recipe::preview_revision(
                inputs,
                &self.target,
                &self.request.recipe,
                assembly_recipe::CleavagePolicy::PositionalChecked,
            )?)?;
        let expected_compiled = ligation_recipe::compile(
            inputs,
            &self.target,
            &self.request.recipe,
            assembly_recipe::CleavagePolicy::PositionalChecked,
        )?;
        if assembly_recipe::hash(&expected_compiled)?
            != assembly_recipe::hash(&self.preview.compiled)?
        {
            return Err(bad(
                "ordered preview no longer compiles from its declared original-source recipe",
            ));
        }
        let n = self.request.recipe.fragments.len();
        let mut missing = BTreeSet::new();
        for i in &self.missing_physical {
            if *i >= n || !missing.insert(*i) {
                return Err(bad("invalid missing-physical role index"));
            }
        }
        for i in 0..n {
            if (i >= self.preview.prepared.len()
                || self.preview.physical.get(i).is_none_or(Option::is_none))
                != missing.contains(&i)
            {
                return Err(bad(
                    "missing physical roles disagree with the fresh prepared prefix",
                ));
            }
        }
        match &self.pool {
            None => {
                if missing.is_empty() && self.preview.preparation_complete
                    || !self.assessments.is_empty()
                {
                    return Err(bad(
                        "withheld pool has no missing preparation or contains invented assessments",
                    ));
                }
            }
            Some(pool) => {
                if !missing.is_empty()
                    || !self.preview.preparation_complete
                    || self.assessments.len() != pool.candidates.len()
                {
                    return Err(bad(
                        "partial preparations cannot be substituted for the complete requested pool",
                    ));
                }
                pool_guard(self, pool)?;
                for (i, a) in self.assessments.iter().enumerate() {
                    let c = &pool.candidates[i];
                    if a.candidate_index != i
                        || a.directed_sequence_key != c.directed_sequence_key
                        || a.physical_sequence_key != c.physical_sequence_key
                        || a.source_roles != roles(pool, i)?
                        || a.target_matches
                            != target_matches(&self.preview.compiled.plan, &c.duplex, &self.target)
                    {
                        return Err(bad(
                            "candidate assessment changed its reference, source roles or target correspondence",
                        ));
                    }
                }
            }
        }
        texts(&self.warnings)?;
        assembly_recipe::hash(self)?;
        Ok(())
    }
}
fn selected_plan(report: &Report, selection: &Selection) -> Result<DesignPlan, Error> {
    report.validate_binding(&selection.context_revision)?;
    if !hex(&selection.directed_sequence_key) || !hex(&selection.physical_sequence_key) {
        return Err(bad("exact candidate sequence keys are required"));
    }
    let pool = report
        .pool
        .as_ref()
        .ok_or_else(|| bad("prepared pool is withheld"))?;
    let candidate = pool
        .candidates
        .get(selection.candidate_index)
        .ok_or_else(|| bad("candidate index is outside this retained pool"))?;
    if candidate.directed_sequence_key != selection.directed_sequence_key
        || candidate.physical_sequence_key != selection.physical_sequence_key
        || !candidate.restriction_site_free
    {
        return Err(bad(
            "selection key changed or candidate retains recognition sites",
        ));
    }
    let mut plan = report.preview.compiled.plan.clone();
    if !matches!(plan.steps.pop(), Some(Operation::Ligate { .. })) {
        return Err(bad(
            "selection requires the old final ordered ligation operation",
        ));
    }
    plan.steps.push(Operation::GoldenGate {
        templates: report
            .preview
            .compiled
            .rows
            .iter()
            .map(|r| r.prepared_source.clone())
            .collect(),
        options: report.request.options.clone(),
        candidate_index: selection.candidate_index,
    });
    Ok(plan)
}
pub fn selection_revision(report: &Report, selection: &Selection) -> Result<String, Error> {
    let plan = selected_plan(report, selection)?;
    assembly_recipe::hash(&(
        "prepared-golden-gate-selection-1",
        &report.context_revision,
        selection,
        &plan,
    ))
}
pub fn select(report: &Report, selection: &Selection) -> Result<SelectionReport, Error> {
    let plan = selected_plan(report, selection)?;
    let context_revision = assembly_recipe::hash(&(
        "prepared-golden-gate-selection-1",
        &report.context_revision,
        selection,
        &plan,
    ))?;
    let mut warnings=vec!["Unverified explicit candidate selection only. All original inputs and preparation/treatment steps are retained; the sole change replaces ordered Ligate with the existing native GoldenGate operation over raw prepared-source references.".into(),"Target correspondence is not a selection rule, unique-product claim or experimental result. Freshly verify this entire native plan against its separately bound target before review or output capture.".into()];
    if report
        .pool
        .as_ref()
        .is_some_and(|p| !p.enumeration_complete)
    {
        warnings.push("This candidate comes from a truncated finite search. Fresh verification retains that incompleteness; omitted products were not ruled out.".into());
    }
    let result = SelectionReport {
        context_revision: context_revision.clone(),
        parent_context: report.context_revision.clone(),
        selection: selection.clone(),
        plan,
        warnings,
    };
    result.validate_binding(&context_revision)?;
    Ok(result)
}
impl SelectionReport {
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        if !hex(expected)
            || self.context_revision != expected
            || !hex(&self.parent_context)
            || self.selection.context_revision != self.parent_context
            || !hex(&self.selection.directed_sequence_key)
            || !hex(&self.selection.physical_sequence_key)
            || self.plan.inputs.is_empty()
            || self.plan.inputs.len() > 32
            || self.plan.steps.is_empty()
            || self.plan.steps.len() > 1000
            || !hex(&self.plan.target_revision)
        {
            return Err(bad("invalid selection context or native plan shape"));
        }
        crate::model::bounds::model_preflight(
            &self
                .plan
                .inputs
                .iter()
                .map(|i| &i.molecule)
                .collect::<Vec<_>>(),
        )?;
        verification::validate_reference_steps(&self.plan)?;
        assembly_recipe::validate_prepared_inputs(&self.plan.inputs, &self.plan.steps)?;
        let mut ids = BTreeSet::new();
        for i in &self.plan.inputs {
            if i.id.is_empty()
                || i.id.len() > 128
                || !ids.insert(i.id.as_str())
                || i.revision != assembly_recipe::hash(&i.molecule)?
            {
                return Err(bad("selected plan original source changed"));
            }
            i.molecule.validate()?;
        }
        let last = self.plan.steps.len() - 1;
        for (index, op) in self.plan.steps.iter().enumerate() {
            if index < last
                && !matches!(
                    op,
                    Operation::Pcr { .. }
                        | Operation::PcrWithMismatches { .. }
                        | Operation::OverlapExtension { .. }
                        | Operation::ReconstructRunProduct { .. }
                        | Operation::ReconstructMaterial { .. }
                        | Operation::AnnealOligos { .. }
                        | Operation::DeclareSourceModifications { .. }
                        | Operation::DeclareSourceTermini { .. }
                        | Operation::CheckedDigest { .. }
                        | Operation::RepairEnds { .. }
                        | Operation::TreatEnds { .. }
                )
            {
                return Err(bad(
                    "selection altered the preparation/reference prefix scope",
                ));
            }
            if let Operation::DeclareSourceModifications { template, .. }
            | Operation::DeclareSourceTermini { template, .. } = op
            {
                let Source::Input(id) = template else {
                    return Err(bad(
                        "original chemistry declaration must retain its exact original input",
                    ));
                };
                let input = self
                    .plan
                    .inputs
                    .iter()
                    .find(|i| &i.id == id)
                    .ok_or_else(|| bad("original chemistry declaration lost its input"))?;
                match op {
                    Operation::DeclareSourceModifications { .. } => {
                        crate::methylation::source_preflight(&input.molecule, &Default::default())?
                    }
                    Operation::DeclareSourceTermini { .. } => {
                        verification::declared_source_termini_preflight(&input.molecule)?
                    }
                    _ => unreachable!(),
                }
            }
            for source in verification::operation_sources(op) {
                match source {
                    Source::Input(id) if ids.contains(id.as_str()) => {}
                    Source::Step(i) if *i < index => {}
                    _ => return Err(bad("selected plan has unknown or forward source")),
                }
            }
        }
        let Operation::GoldenGate {
            templates,
            options: pool_options,
            candidate_index,
        } = self
            .plan
            .steps
            .last()
            .ok_or_else(|| bad("missing selected operation"))?
        else {
            return Err(bad("selected plan must end with native GoldenGate"));
        };
        options(pool_options)?;
        if templates.is_empty()
            || templates.len() > 32
            || *candidate_index != self.selection.candidate_index
            || *candidate_index >= pool_options.max_products
        {
            return Err(bad("selected candidate/template mapping changed"));
        }
        if assembly_recipe::hash(&(
            "prepared-golden-gate-selection-1",
            &self.parent_context,
            &self.selection,
            &self.plan,
        ))? != self.context_revision
        {
            return Err(bad("selected plan revision changed"));
        }
        texts(&self.warnings)?;
        assembly_recipe::hash(self)?;
        Ok(())
    }
}
