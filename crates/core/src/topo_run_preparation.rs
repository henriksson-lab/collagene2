//! Exact saved-material replay as a template for new TOPO insert PCR.
//! Documentary bounds/bindings never reconstruct or return physical authority.
use crate::{
    Error, Molecule,
    run_material::MaterialOrigin,
    topo, topo_actions, topo_worker,
    verification::{DesignPlan, Operation, PlanInput, Source, TopoSource},
};
use serde::{Deserialize, Serialize};
pub const ENGINE_VERSION: &str = "topo-saved-template-1";
macro_rules! dto {($($t:tt)*)=>{#[derive(Debug,Clone,Serialize,Deserialize)]#[cfg_attr(feature="json-schema",derive(schemars::JsonSchema))]#[serde(deny_unknown_fields)]$($t)*};}
dto! {pub enum Task {
 Orientations{request:topo::Request,vector:[usize;2],insert:[usize;2],directional:bool},
 Project{request:topo::Request,selection:topo::Selection},
 Check{request:topo::Request,selection:topo::Selection,target:Molecule,policy:topo::AcceptancePolicy},
 InspectActions{request:topo_actions::ActionsRequest},
 ProcessActions{request:topo_actions::Request},
}}
dto! {pub enum Reply {Orientations(topo::OrientationInspection),Projected(topo::AnnotatedProduct),Checked(topo::AnnotatedAssessment),ActionsInspected(topo_actions::Inspection),ActionsProcessed(topo_actions::Report)}}
fn bad(s: &str) -> Error {
    Error::Unsupported(format!("TOPO saved template: {s}"))
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
/// Borrowed bounds before cloning rich source models. This is declaration
/// validation, not a grant of source chemistry; inspect always replays the capsule.
pub fn source_inputs(
    sources: &[TopoSource],
    molecules: &[&Molecule],
    request: &topo::Request,
    template_role: &str,
    origin: Option<&MaterialOrigin>,
) -> Result<Vec<topo::Input>, Error> {
    if sources.is_empty() || sources.len() > 32 || sources.len() != molecules.len() {
        return Err(Error::TooLarge);
    }
    crate::model::bounds::model_preflight(molecules)?;
    let mut sum = 0usize;
    for (s, m) in sources.iter().zip(molecules) {
        if s.role.is_empty()
            || s.role.len() > 256
            || !matches!(&s.source,Source::Input(id) if !id.is_empty()&&id.len()<=128)
        {
            return Err(bad("capture requires bounded original input roles"));
        }
        sum = sum
            .checked_add(m.sequence.len())
            .filter(|n| *n <= 100_000)
            .ok_or(Error::TooLarge)?;
    }
    if let Some(origin) = origin {
        origin.validate_documentary()?;
        let index = sources
            .iter()
            .position(|s| s.role == template_role)
            .ok_or_else(|| bad("missing selected template role"))?;
        if molecules[index].revision() != origin.event().revision {
            return Err(bad("selected saved event/revision changed"));
        }
    }
    let inputs = sources
        .iter()
        .zip(molecules)
        .map(|(s, m)| topo::Input {
            role: s.role.clone(),
            molecule: (*m).clone(),
        })
        .collect::<Vec<_>>();
    if origin.is_some() {
        topo::declared_template_preflight(&inputs, request, template_role)?
    } else {
        topo::source_preflight(&inputs, request)?
    }
    Ok(inputs)
}
pub fn preflight(
    inputs: &[topo::Input],
    template_role: &str,
    origin: &MaterialOrigin,
    task: &Task,
) -> Result<(), Error> {
    topo::declared_template_preflight(inputs, task.preparation(), template_role)?;
    origin.validate_documentary()?;
    let source = inputs
        .iter()
        .find(|i| i.role == template_role)
        .ok_or_else(|| bad("missing selected template"))?;
    if source.molecule.revision() != origin.event().revision {
        return Err(bad("selected saved template revision changed"));
    }
    if let Task::Check { target, .. } = task {
        topo_worker::target_preflight(target)?;
    }
    // Bounded serializer protects literal option trees before any hashing/replay.
    crate::assembly_recipe::hash(&(origin, task))?;
    Ok(())
}
fn hash(v: &impl Serialize) -> Result<String, Error> {
    crate::assembly_recipe::hash(v)
}
pub fn revision(
    inputs: &[topo::Input],
    role: &str,
    origin: &MaterialOrigin,
    task: &Task,
) -> Result<String, Error> {
    preflight(inputs, role, origin, task)?;
    hash(&(ENGINE_VERSION, topo::identities(inputs), role, origin, task))
}
/// Original report binding, for readonly validation only. The outer worker
/// revision additionally binds the entire retained native origin and task.
pub fn documentary_revision(
    inputs: &[topo::Input],
    role: &str,
    origin: &MaterialOrigin,
    task: &Task,
) -> Result<String, Error> {
    preflight(inputs, role, origin, task)?;
    let sources = topo::identities(inputs);
    let prep = task.preparation();
    let batch = hash(&(
        topo::ENGINE_VERSION,
        topo::physical::ENGINE,
        topo::physical::PRIMARY_MODEL,
        topo::physical::projection::ENGINE,
        &sources,
        prep,
    ))?;
    let product = |selection: &topo::Selection| -> Result<String, Error> {
        if selection.batch_revision != batch || selection.component >= 64 {
            return Err(bad("selection belongs to another preparation"));
        }
        hash(&(
            topo::ENGINE_VERSION,
            "selected-physical-component",
            &batch,
            selection,
        ))
    };
    let annotated = |selection: &topo::Selection| -> Result<String, Error> {
        hash(&(
            topo::ANNOTATION_ENGINE_VERSION,
            "annotated-product",
            product(selection)?,
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
            &sources,
            prep,
            vector,
            insert,
            directional,
        )),
        Task::Project { selection, .. } => annotated(selection),
        Task::Check {
            selection,
            target,
            policy,
            ..
        } => hash(&(
            topo::ANNOTATION_ENGINE_VERSION,
            "annotated-target",
            annotated(selection)?,
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
    role: &str,
    origin: &MaterialOrigin,
    task: &Task,
) -> Result<Inspection, Error> {
    preflight(inputs, role, origin, task)?;
    let source = &inputs
        .iter()
        .find(|i| i.role == role)
        .ok_or_else(|| bad("missing selected template"))?
        .molecule;
    let fresh = origin.replay(source)?;
    let reply =
        topo::inspect_reconstructed_task(inputs, &[(role.to_string(), fresh.tracked())], task)?;
    reply.validate_binding(
        inputs,
        role,
        origin,
        task,
        &revision(inputs, role, origin, task)?,
    )?;
    Ok(Inspection {
        context_revision: revision(inputs, role, origin, task)?,
        origin_incomplete: fresh.incomplete(),
        origin_findings_total: fresh.findings_total(),
        origin_findings_truncated: fresh.findings_truncated(),
        origin_findings: fresh.findings().to_vec(),
        reply,
    })
}
impl Reply {
    pub fn validate_binding(
        &self,
        inputs: &[topo::Input],
        role: &str,
        origin: &MaterialOrigin,
        task: &Task,
        expected: &str,
    ) -> Result<(), Error> {
        if revision(inputs, role, origin, task)? != expected {
            return Err(bad("saved source/task binding changed"));
        }
        let inner = documentary_revision(inputs, role, origin, task)?;
        match (self, task) {
            (Self::Orientations(r), Task::Orientations { .. }) => {
                topo_worker::validate_orientations(r, &inner)
            }
            (Self::Projected(r), Task::Project { .. }) => {
                topo_worker::validate_annotated_projection(r, &inner)
            }
            (Self::Checked(r), Task::Check { .. }) => {
                topo_worker::validate_annotated_target(r, &inner)
            }
            (Self::ActionsInspected(r), Task::InspectActions { .. }) => r.validate_binding(&inner),
            (Self::ActionsProcessed(r), Task::ProcessActions { .. }) => r.validate_binding(&inner),
            _ => Err(bad("reply kind differs from requested inspection")),
        }
    }
}
/// Compile only declarations; ordinary verification reconstructs the entire origin
/// and new PCR again. No graph, report, metadata or saved sequence grants material.
pub fn compile(
    inputs: &[topo::Input],
    target: &Molecule,
    role: &str,
    origin: &MaterialOrigin,
    request: &topo_actions::Request,
) -> Result<DesignPlan, Error> {
    preflight(
        inputs,
        role,
        origin,
        &Task::ProcessActions {
            request: request.clone(),
        },
    )?;
    topo_worker::target_preflight(target)?;
    let plan = DesignPlan {
        name: "TOPO PCR from saved native material".into(),
        target_revision: target.revision(),
        inputs: inputs
            .iter()
            .map(|i| PlanInput {
                id: i.role.clone(),
                revision: i.molecule.revision(),
                molecule: i.molecule.clone(),
            })
            .collect(),
        steps: vec![
            match origin {
                MaterialOrigin::Run(origin) => Operation::ReconstructRunProduct {
                    template: Source::Input(role.into()),
                    origin: origin.clone(),
                },
                MaterialOrigin::Derived(origin) => Operation::ReconstructMaterial {
                    template: Source::Input(role.into()),
                    origin: origin.clone(),
                },
            },
            Operation::TopoActions {
                sources: inputs
                    .iter()
                    .map(|i| TopoSource {
                        role: i.role.clone(),
                        source: if i.role == role {
                            Source::Step(0)
                        } else {
                            Source::Input(i.role.clone())
                        },
                    })
                    .collect(),
                request: request.clone(),
            },
        ],
    };
    crate::run_material::plan_preflight(&plan)?;
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
            Reply::Projected(r) => topo_worker::validate_annotated_projection(r, &r.batch_revision),
            Reply::Checked(r) => topo_worker::validate_annotated_target(r, &r.batch_revision),
            Reply::ActionsInspected(r) => r.validate_binding(&r.context_revision),
            Reply::ActionsProcessed(r) => r.validate_binding(&r.context_revision),
        }
    }
    pub fn validate_for(
        &self,
        inputs: &[topo::Input],
        role: &str,
        origin: &MaterialOrigin,
        task: &Task,
        expected: &str,
    ) -> Result<(), Error> {
        self.validate_binding(expected)?;
        self.reply
            .validate_binding(inputs, role, origin, task, expected)
    }
}
/// Documentary preparation identity. Does not execute preparation or replay.
pub fn preparation_revision(
    inputs: &[topo::Input],
    role: &str,
    origin: &MaterialOrigin,
    request: &topo::Request,
) -> Result<String, Error> {
    preflight(
        inputs,
        role,
        origin,
        &Task::Orientations {
            request: request.clone(),
            vector: [0, 1],
            insert: [2, 3],
            directional: false,
        },
    )?;
    hash(&(
        topo::ENGINE_VERSION,
        topo::physical::ENGINE,
        topo::physical::PRIMARY_MODEL,
        topo::physical::projection::ENGINE,
        topo::identities(inputs),
        request,
    ))
}
