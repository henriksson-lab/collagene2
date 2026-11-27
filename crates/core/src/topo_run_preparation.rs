//! Fresh native TOPO sources: original declarations, saved ordinary PCR templates,
//! and exact retained elementary graph origins. No serialized graph supplies state.
use crate::{
    Error, Molecule,
    run_material::{self, MaterialOrigin},
    topo, topo_actions, topo_worker,
    verification::{DesignPlan, Operation, PlanInput, Source, TopoSource},
};
use serde::{Deserialize, Serialize};
pub const ENGINE_VERSION: &str = "topo-native-sources-1";
macro_rules! dto {($($t:tt)*)=>{#[derive(Debug,Clone,Serialize,Deserialize)]#[cfg_attr(feature="json-schema",derive(schemars::JsonSchema))]#[serde(deny_unknown_fields)]$($t)*};}
dto! {pub struct BoundTemplate{pub role:String,pub origin:MaterialOrigin}}
dto! {pub enum Task {
 Orientations{request:topo::Request,vector:[usize;2],insert:[usize;2],directional:bool},
 Project{request:topo::Request,selection:topo::Selection},
 Check{request:topo::Request,selection:topo::Selection,target:Molecule,policy:topo::AcceptancePolicy},
 InspectActions{request:topo_actions::ActionsRequest},
 ProcessActions{request:topo_actions::Request},
}}
dto! {pub enum Reply {Orientations(topo::OrientationInspection),Projected(topo::Product),Checked(topo::TargetAssessment),ActionsInspected(topo_actions::Inspection),ActionsProcessed(topo_actions::Report)}}
fn bad(s: &str) -> Error {
    Error::Unsupported(format!("TOPO native sources: {s}"))
}
fn hash(v: &impl Serialize) -> Result<String, Error> {
    crate::assembly_recipe::hash(v)
}
impl Task {
    pub fn preparation(&self) -> &topo::Request {
        match self {
            Self::Orientations { request, .. }
            | Self::Project { request, .. }
            | Self::Check { request, .. } => request,
            Self::InspectActions { request } => &request.preparation,
            Self::ProcessActions { request } => &request.preparation,
        }
    }
}
pub fn source_inputs(
    sources: &[TopoSource],
    molecules: &[&Molecule],
    request: &topo::Request,
    template: Option<&BoundTemplate>,
) -> Result<Vec<topo::Input>, Error> {
    if sources.len() > 32 || sources.len() != molecules.len() {
        return Err(Error::TooLarge);
    }
    crate::model::bounds::model_preflight(molecules)?;
    for s in sources {
        if s.role.is_empty()
            || s.role.len() > 256
            || !matches!(&s.source,Source::Input(id)if !id.is_empty()&&id.len()<=128)
        {
            return Err(bad("bounded original input role required"));
        }
    }
    let roots = template
        .into_iter()
        .map(|t| match &t.origin {
            MaterialOrigin::Run(c) => run_material::NativeOriginRef::Run(c),
            MaterialOrigin::Derived(c) => run_material::NativeOriginRef::Derived(c),
        })
        .chain(request.recipe.materials.iter().filter_map(|m| {
            if let topo::MaterialRecipe::SavedGraph { origin, .. } = m {
                Some(run_material::NativeOriginRef::TopoGraph(origin))
            } else {
                None
            }
        }));
    run_material::native_origins_preflight(roots)?;
    let inputs = sources
        .iter()
        .zip(molecules)
        .map(|(s, m)| topo::Input {
            role: s.role.clone(),
            molecule: (*m).clone(),
        })
        .collect::<Vec<_>>();
    if let Some(t) = template {
        topo::declared_template_preflight(&inputs, request, &t.role)?;
        if inputs
            .iter()
            .find(|i| i.role == t.role)
            .map(|i| i.molecule.revision())
            != Some(t.origin.event().revision.clone())
        {
            return Err(bad("saved template exact revision differs"));
        }
    } else {
        topo::source_preflight(&inputs, request)?
    }
    Ok(inputs)
}
pub fn preflight(
    inputs: &[topo::Input],
    template: Option<&BoundTemplate>,
    task: &Task,
) -> Result<(), Error> {
    let roots = template
        .into_iter()
        .map(|t| match &t.origin {
            MaterialOrigin::Run(c) => run_material::NativeOriginRef::Run(c),
            MaterialOrigin::Derived(c) => run_material::NativeOriginRef::Derived(c),
        })
        .chain(task.preparation().recipe.materials.iter().filter_map(|m| {
            if let topo::MaterialRecipe::SavedGraph { origin, .. } = m {
                Some(run_material::NativeOriginRef::TopoGraph(origin))
            } else {
                None
            }
        }));
    run_material::native_origins_preflight(roots)?;
    if let Some(t) = template {
        topo::declared_template_preflight(inputs, task.preparation(), &t.role)?;
        if inputs
            .iter()
            .find(|i| i.role == t.role)
            .map(|i| i.molecule.revision())
            != Some(t.origin.event().revision.clone())
        {
            return Err(bad("saved template source changed"));
        }
    } else {
        topo::source_preflight(inputs, task.preparation())?
    }
    if let Task::Check { target, .. } = task {
        topo_worker::target_preflight(target)?;
    }
    hash(&(inputs, template, task))?;
    Ok(())
}
pub fn revision(
    inputs: &[topo::Input],
    template: Option<&BoundTemplate>,
    task: &Task,
) -> Result<String, Error> {
    preflight(inputs, template, task)?;
    hash(&(ENGINE_VERSION, topo::identities(inputs), template, task))
}
pub fn preparation_revision(
    inputs: &[topo::Input],
    template: Option<&BoundTemplate>,
    request: &topo::Request,
) -> Result<String, Error> {
    if let Some(t) = template {
        topo::declared_template_preflight(inputs, request, &t.role)?;
    } else {
        topo::source_preflight(inputs, request)?;
    }
    topo::documentary_batch_revision(inputs, request)
}
pub fn documentary_revision(
    inputs: &[topo::Input],
    template: Option<&BoundTemplate>,
    task: &Task,
) -> Result<String, Error> {
    preflight(inputs, template, task)?;
    let prep = task.preparation();
    let batch = topo::documentary_batch_revision(inputs, prep)?;
    let sources = topo::documentary_sources(inputs, prep)?;
    let product = |s: &topo::Selection| -> Result<String, Error> {
        if s.batch_revision != batch || s.component >= 64 {
            return Err(bad("selected native preparation differs"));
        }
        hash(&(
            topo::ENGINE_VERSION,
            "selected-physical-component",
            &batch,
            s,
        ))
    };
    match task {
        Task::Orientations {
            vector,
            insert,
            directional,
            ..
        } => hash(&(
            topo::ENGINE_VERSION,
            "two-declared-orientations",
            topo::identities(inputs),
            prep,
            vector,
            insert,
            directional,
        )),
        Task::Project { selection, .. } => product(selection),
        Task::Check {
            selection,
            target,
            policy,
            ..
        } => hash(&(
            topo::ENGINE_VERSION,
            "independent-target",
            product(selection)?,
            target.revision(),
            policy,
        )),
        Task::InspectActions { request } => hash(&(
            topo_actions::ENGINE_VERSION,
            "actions-discovery",
            topo::physical::processing::ENGINE,
            &batch,
            &sources,
            request,
        )),
        Task::ProcessActions { request } => hash(&(
            topo_actions::ENGINE_VERSION,
            topo::physical::processing_projection::ENGINE,
            &batch,
            &sources,
            request,
        )),
    }
}
pub fn inspect(
    inputs: &[topo::Input],
    template: Option<&BoundTemplate>,
    task: &Task,
) -> Result<Inspection, Error> {
    preflight(inputs, template, task)?;
    let fresh = template
        .map(|t| {
            let m = &inputs
                .iter()
                .find(|i| i.role == t.role)
                .ok_or_else(|| bad("saved template missing"))?
                .molecule;
            t.origin.replay(m)
        })
        .transpose()?;
    let tracked = template
        .zip(fresh.as_ref())
        .map(|(t, f)| vec![(t.role.clone(), f.tracked())])
        .unwrap_or_default();
    let reply = topo::inspect_reconstructed_task(inputs, &tracked, task)?;
    let mut findings = vec![];
    let mut incomplete = false;
    let mut omitted_source_findings = 0usize;
    let mut source_findings_truncated = false;
    if let Some(f) = &fresh {
        incomplete |= f.incomplete();
        findings.extend(f.findings().iter().cloned());
        omitted_source_findings = f
            .findings_total()
            .checked_sub(f.findings().len())
            .ok_or_else(|| bad("native source findings count is inconsistent"))?;
        source_findings_truncated = f.findings_truncated() || omitted_source_findings > 0;
    }
    for m in &task.preparation().recipe.materials {
        if let topo::MaterialRecipe::SavedGraph { origin, .. } = m {
            let r = origin.validate_shallow()?;
            incomplete |= r.incomplete;
            findings.extend(r.state.findings);
        }
    }
    let evidence = match &reply {
        Reply::Orientations(r) => &r.preparation,
        Reply::Projected(r) => &r.preparation,
        Reply::Checked(r) => &r.product.preparation,
        Reply::ActionsInspected(r) => &r.preparation,
        Reply::ActionsProcessed(r) => &r.preparation,
    };
    for origin in &evidence.graph_origins {
        if origin.incomplete {
            incomplete = true;
            findings.push(format!("Saved native graph {} retains incomplete source history, attested against the same fresh scoped execution. Target/policy agreement does not upgrade that source.",origin.selection.record_hash));
        }
    }
    let total = findings
        .len()
        .checked_add(omitted_source_findings)
        .ok_or(Error::TooLarge)?;
    let mut truncated = source_findings_truncated || total > 128;
    let findings = findings
        .into_iter()
        .take(128)
        .map(|s| {
            if s.chars().count() > 2048 {
                truncated = true;
            }
            s.chars().take(2048).collect()
        })
        .collect();
    let out = Inspection {
        context_revision: revision(inputs, template, task)?,
        origin_incomplete: incomplete,
        origin_findings: findings,
        origin_findings_total: total,
        origin_findings_truncated: truncated,
        reply,
    };
    out.validate_for(inputs, template, task, &out.context_revision)?;
    Ok(out)
}
impl Reply {
    pub fn validate_binding(
        &self,
        inputs: &[topo::Input],
        template: Option<&BoundTemplate>,
        task: &Task,
        expected: &str,
    ) -> Result<(), Error> {
        if revision(inputs, template, task)? != expected {
            return Err(bad("native source/task receipt changed"));
        }
        let inner = documentary_revision(inputs, template, task)?;
        match (self, task) {
            (Self::Orientations(r), Task::Orientations { .. }) => {
                topo_worker::validate_orientations(r, &inner)
            }
            (Self::Projected(r), Task::Project { .. }) => {
                topo_worker::validate_projection(r, &inner)
            }
            (Self::Checked(r), Task::Check { .. }) => topo_worker::validate_target(r, &inner),
            (Self::ActionsInspected(r), Task::InspectActions { .. }) => r.validate_binding(&inner),
            (Self::ActionsProcessed(r), Task::ProcessActions { .. }) => r.validate_binding(&inner),
            _ => Err(bad("reply kind differs")),
        }
    }
}
pub fn compile(
    inputs: &[topo::Input],
    target: &Molecule,
    template: Option<&BoundTemplate>,
    task: &Task,
) -> Result<DesignPlan, Error> {
    preflight(inputs, template, task)?;
    topo_worker::target_preflight(target)?;
    let mut steps = vec![];
    if let Some(t) = template {
        steps.push(match &t.origin {
            MaterialOrigin::Run(origin) => Operation::ReconstructRunProduct {
                template: Source::Input(t.role.clone()),
                origin: origin.clone(),
            },
            MaterialOrigin::Derived(origin) => Operation::ReconstructMaterial {
                template: Source::Input(t.role.clone()),
                origin: origin.clone(),
            },
        });
    }
    let sources = inputs
        .iter()
        .map(|i| TopoSource {
            role: i.role.clone(),
            source: if template.is_some_and(|t| t.role == i.role) {
                Source::Step(0)
            } else {
                Source::Input(i.role.clone())
            },
        })
        .collect();
    steps.push(match task {
        Task::Check {
            request,
            selection,
            target: checked,
            policy,
        } => {
            if checked != target {
                return Err(bad(
                    "checked independent target differs from requested plan target",
                ));
            }
            Operation::CheckedTopo {
                sources,
                request: request.clone(),
                selection: selection.clone(),
                policy: policy.clone(),
            }
        }
        Task::ProcessActions { request } => Operation::TopoActions {
            sources,
            request: request.clone(),
        },
        _ => {
            return Err(bad(
                "only checked target or native processed product can compile a candidate",
            ));
        }
    });
    let plan = DesignPlan {
        name: "TOPO from declared native sources".into(),
        target_revision: target.revision(),
        inputs: inputs
            .iter()
            .map(|i| PlanInput {
                id: i.role.clone(),
                revision: i.molecule.revision(),
                molecule: i.molecule.clone(),
            })
            .collect(),
        steps,
    };
    run_material::plan_preflight(&plan)?;
    Ok(plan)
}
dto! {pub struct Inspection{pub context_revision:String,pub origin_incomplete:bool,pub origin_findings:Vec<String>,pub origin_findings_total:usize,pub origin_findings_truncated:bool,pub reply:Reply}}
impl Inspection {
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        if self.origin_findings_total > 32 * 1024 * 1024
            || self.origin_findings.len() != self.origin_findings_total.min(128)
            || self.origin_findings_total > self.origin_findings.len()
                && !self.origin_findings_truncated
        {
            return Err(Error::TooLarge);
        }
        if self.origin_findings.len() > 128
            || self
                .origin_findings
                .iter()
                .any(|s| s.chars().count() > 2048 || s.len() > 8192)
        {
            return Err(Error::TooLarge);
        }
        if expected.len() != 64
            || !expected
                .bytes()
                .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
            || self.context_revision != expected
        {
            return Err(bad("outer saved-template inspection binding changed"));
        }
        match &self.reply {
            Reply::Orientations(r) => topo_worker::validate_orientations(r, &r.batch_revision),
            Reply::Projected(r) => topo_worker::validate_projection(r, &r.batch_revision),
            Reply::Checked(r) => topo_worker::validate_target(r, &r.batch_revision),
            Reply::ActionsInspected(r) => r.validate_binding(&r.context_revision),
            Reply::ActionsProcessed(r) => r.validate_binding(&r.context_revision),
        }
    }
    pub fn validate_for(
        &self,
        inputs: &[topo::Input],
        template: Option<&BoundTemplate>,
        task: &Task,
        expected: &str,
    ) -> Result<(), Error> {
        self.validate_binding(expected)?;
        self.reply
            .validate_binding(inputs, template, task, expected)
    }
}
