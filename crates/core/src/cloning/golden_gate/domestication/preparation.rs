//! Exact-source physical preparation for explicitly checked domestication.
//! A reference-only recoded target never becomes an invented PCR template.
use crate::{
    cloning::{
        golden_gate::domestication,
        pcr::mutagenesis::{self, InversePcrOptions},
    },
    documents::workspace::Workspace,
    material::run::MaterialOrigin,
    planning::plan::{Operation, PlanInput, Source},
    verification,
};
use collagene_model::{Error, Molecule, Topology};
use serde::{Deserialize, Serialize};
use std::{collections::BTreeSet, io};

pub const ENGINE_VERSION: &str = "domestication-preparation-3";
/// Current documentary review receipt, not executable material authority.
pub const REVIEW_EVIDENCE_KIND: &str = "domestication-intent-v3";
pub const REVIEW_EVIDENCE_VERSION: u16 = 3;
pub const INPUT_ID: &str = "template";
const MAX_SOURCE_BASES: usize = 8 << 20;
const MAX_WIRE_BYTES: usize = 32 << 20;
macro_rules! dto {($($tokens:tt)*)=>{#[derive(Clone,Debug,Serialize,Deserialize)]#[cfg_attr(feature="json-schema",derive(schemars::JsonSchema))]#[serde(deny_unknown_fields)]$($tokens)*};}
dto! {pub enum Route {
    InversePcr {edit_start:usize,edit_end:usize,forward_length:usize,reverse_length:usize,reverse_tail_bases:usize},
    OverlapExtension {options:crate::cloning::overlap::extension::primers::Options},
    InversePcrRounds {rounds:Vec<InverseRound>},
}}
dto! {pub struct InverseRound {
    pub edit_start:usize,pub edit_end:usize,pub forward_length:usize,pub reverse_length:usize,pub reverse_tail_bases:usize,
}}
dto! {pub struct RoundEvidence {
    pub round:InverseRound,pub template:Source,pub source_origin:usize,pub options:InversePcrOptions,
    pub forward:crate::sequence::analysis::Primer,pub reverse:crate::sequence::analysis::Primer,
    pub forward_location:collagene_model::Location,pub reverse_location:collagene_model::Location,pub retained_template_length:usize,
}}
dto! {pub struct Request {pub design:domestication::Request,pub selection:domestication::Selection,pub route:Route,#[serde(deserialize_with="crate::encoding::required_option")]pub origin:Option<MaterialOrigin>}}
dto! {pub enum MethodEvidence {
    InversePcr {options:InversePcrOptions,forward:crate::sequence::analysis::Primer,reverse:crate::sequence::analysis::Primer,forward_location:collagene_model::Location,reverse_location:collagene_model::Location,retained_template_length:usize},
    OverlapExtension {context_revision:String,source_frame:crate::cloning::overlap::extension::primers::SourceFrame,preparation:crate::cloning::overlap::extension::pcr::Request,#[cfg_attr(feature="json-schema",schemars(with="serde_json::Value"))]designs:Vec<crate::cloning::pcr::primers::DesignReport>,search:Vec<crate::cloning::overlap::extension::primers::PairSearch>},
    InversePcrRounds {rounds:Vec<RoundEvidence>},
}}
dto! {pub struct Prepared {
    pub context_revision:String,pub request:Request,pub source_revision:String,pub target:Molecule,
    pub check:domestication::CheckReport,pub plan:crate::planning::plan::DesignPlan,pub method:MethodEvidence,
    pub verification:verification::evidence::VerificationReport,pub proposal_ready:bool,pub proposal_findings:Vec<String>,pub warnings:Vec<String>,
}}
fn bad(message: &str) -> Error {
    Error::Unsupported(format!("domestication preparation: {message}"))
}
fn digest(value: &str) -> bool {
    value.len() == 64
        && value
            .bytes()
            .all(|byte| byte.is_ascii_digit() || (b'a'..=b'f').contains(&byte))
}
struct Bytes(usize);
impl io::Write for Bytes {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        self.0 = self
            .0
            .checked_add(bytes.len())
            .filter(|size| *size <= MAX_WIRE_BYTES)
            .ok_or_else(|| io::Error::other("domestication preparation wire budget"))?;
        Ok(bytes.len())
    }
    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}
fn dimensions(
    input: &PlanInput,
    request: &Request,
) -> Result<(usize, usize, usize, usize, usize), Error> {
    let Route::InversePcr {
        edit_start,
        edit_end,
        forward_length,
        reverse_length,
        reverse_tail_bases,
    } = &request.route
    else {
        return Err(bad(
            "inverse-PCR dimensions require the declared inverse-PCR route",
        ));
    };
    let (edit_start, edit_end, forward_length, reverse_length, reverse_tail_bases) = (
        *edit_start,
        *edit_end,
        *forward_length,
        *reverse_length,
        *reverse_tail_bases,
    );
    let n = input.molecule.sequence.len();
    if n == 0 || n > MAX_SOURCE_BASES || input.molecule.topology != Topology::Circular {
        return Err(bad(
            "this explicit physical route requires a bounded circular DNA original; other routes remain unsupported",
        ));
    }
    let span = edit_end
        .checked_sub(edit_start)
        .filter(|span| *span > 0 && *span < n)
        .ok_or_else(|| bad("replace a nonempty interval shorter than one original turn"))?;
    if edit_start >= n
        || edit_start.checked_add(n).is_none_or(|end| edit_end > end)
        || !(8..=80).contains(&forward_length)
        || !(8..=80).contains(&reverse_length)
        || forward_length + reverse_length > n - span
    {
        return Err(bad(
            "choose original-source coordinates and nonoverlapping8–80 base annealing domains",
        ));
    }
    if reverse_tail_bases > span
        || reverse_tail_bases > 200 - reverse_length
        || span - reverse_tail_bases > 200 - forward_length
    {
        return Err(bad(
            "every substituted base must fit the explicitly split actual primer tails, with each complete ordered primer at most200 bases; distant edits require another physical route or reviewed rounds",
        ));
    }
    if request.selection.patches.is_empty() {
        return Err(bad(
            "no substitutions were selected for a physical preparation",
        ));
    }
    let mut seen = BTreeSet::new();
    for patch in &request.selection.patches {
        if patch.position >= n {
            return Err(bad(
                "a selected substitution lies outside the original source",
            ));
        }
        let offset = if patch.position >= edit_start {
            patch.position - edit_start
        } else {
            patch.position + n - edit_start
        };
        if patch.position >= n || offset >= span || !seen.insert(patch.position) {
            return Err(bad(
                "each selected substitution must occur exactly once inside the explicit replacement interval",
            ));
        }
        if patch.expected.len() != 1
            || patch.replacement.len() != 1
            || !patch.expected.bytes().all(|base| b"ACGT".contains(&base))
            || !patch
                .replacement
                .bytes()
                .all(|base| b"ACGT".contains(&base))
            || patch.expected == patch.replacement
            || input.molecule.sequence.as_str().as_bytes()[patch.position]
                != patch.expected.as_bytes()[0]
        {
            return Err(bad(
                "selected substitutions must name the actual concrete original base and a distinct concrete replacement",
            ));
        }
    }
    Ok((
        edit_start,
        edit_end,
        forward_length,
        reverse_length,
        reverse_tail_bases,
    ))
}
/// Cheap bounded identity/geometry guards; no genetic checking or PCR is replayed.
fn documentary_preflight(input: &PlanInput, request: &Request) -> Result<(), Error> {
    crate::material::capsule::material_origins_preflight(request.origin.iter())?;
    if input.id != INPUT_ID
        || input.revision != request.design.source_revision
        || !digest(&input.revision)
    {
        return Err(bad(
            "the explicit template input and design must bind the same original source revision",
        ));
    }
    collagene_model::bounds::model_preflight(&[&input.molecule])?;
    if input.molecule.sequence.len() > MAX_SOURCE_BASES {
        return Err(Error::TooLarge);
    }
    // Bound the borrowed request before hashing it or copying route options.
    serde_json::to_writer(&mut Bytes(0), &(input, request)).map_err(|_| Error::TooLarge)?;
    input.molecule.validate()?;
    if input.molecule.revision() != input.revision {
        return Err(bad("original template content changed"));
    }
    // Native domestication preflight owns catalog, feature, patch and search bounds.
    domestication::check_revision(&input.molecule, &request.design, &request.selection)?;
    if let Some(origin) = &request.origin {
        origin.validate_documentary()?;
        if origin.event().revision != input.revision {
            return Err(bad("saved origin differs from exact input revision"));
        }
    } else if crate::chemistry::has_unmapped_physical_evidence(&input.molecule) {
        return Err(bad(
            "these routes require an original DNA reference; saved/imported physical material must retain a separately supported native origin, never lose its markers",
        ));
    }
    match &request.route {
        Route::InversePcr { .. } => {
            dimensions(input, request)?;
        }
        Route::OverlapExtension { .. } => {
            let native = oe_request(input, request)?;
            crate::cloning::overlap::extension::primers::documentary_preflight(
                &input.molecule,
                &native,
            )?;
            if request.origin.is_none() {
                let (work, material) = crate::cloning::overlap::extension::primers::work_preflight(
                    &input.molecule,
                    &native,
                )?;
                let extra = if input.molecule.topology == Topology::Circular {
                    input
                        .molecule
                        .sequence
                        .len()
                        .checked_mul(64)
                        .ok_or(Error::TooLarge)?
                } else {
                    0
                };
                if work.checked_add(extra).is_none_or(|n| n > 100_000_000)
                    || material
                        .checked_add(extra)
                        .is_none_or(|n| n > 64 * 1024 * 1024)
                {
                    return Err(Error::TooLarge);
                }
            }
        }
        Route::InversePcrRounds { rounds } => {
            crate::cloning::golden_gate::domestication::rounds::preflight(
                input,
                &request.selection.patches,
                rounds,
            )?
        }
    }
    Ok(())
}
pub fn preflight(
    input: &PlanInput,
    request: &Request,
    workspace: Option<&Workspace>,
) -> Result<(), Error> {
    collagene_model::bounds::model_preflight(&[&input.molecule])?;
    crate::material::capsule::material_origins_preflight(request.origin.iter())?;
    if let Some(origin) = &request.origin {
        let workspace =
            workspace.ok_or_else(|| bad("saved material requires its exact owning workspace"))?;
        origin.validate_for_workspace(workspace, origin.event())?;
        if crate::material::run::resolve(workspace, origin.event())? != &input.molecule {
            return Err(bad(
                "saved source snapshot differs from its exact native event",
            ));
        }
    }
    documentary_preflight(input, request)
}
fn documentary_revision(input: &PlanInput, request: &Request) -> Result<String, Error> {
    documentary_preflight(input, request)?;
    crate::encoding::bounded_json::hash(&(ENGINE_VERSION, input, request))
}
pub fn revision(
    input: &PlanInput,
    request: &Request,
    workspace: Option<&Workspace>,
) -> Result<String, Error> {
    preflight(input, request, workspace)?;
    crate::encoding::bounded_json::hash(&(ENGINE_VERSION, input, request))
}
pub(crate) fn append_origin(
    input: &str,
    origin: Option<&MaterialOrigin>,
    steps: &mut Vec<Operation>,
) -> Source {
    let template = Source::Input(input.into());
    if let Some(origin) = origin {
        steps.push(match origin {
            MaterialOrigin::Run(origin) => Operation::ReconstructRunProduct {
                template,
                origin: origin.clone(),
            },
            MaterialOrigin::Derived(origin) => Operation::ReconstructMaterial {
                template,
                origin: origin.clone(),
            },
        });
        Source::Step(steps.len() - 1)
    } else {
        template
    }
}
fn options(input: &PlanInput, request: &Request) -> Result<InversePcrOptions, Error> {
    let (edit_start, edit_end, forward_length, reverse_length, reverse_tail_bases) =
        dimensions(input, request)?;
    let n = input.molecule.sequence.len();
    let original = input.molecule.sequence.as_str().as_bytes();
    let mut inserted = (edit_start..edit_end)
        .map(|position| original[position % n])
        .collect::<Vec<_>>();
    for patch in &request.selection.patches {
        let offset = if patch.position >= edit_start {
            patch.position - edit_start
        } else {
            patch.position + n - edit_start
        };
        inserted[offset] = patch.replacement.as_bytes()[0];
    }
    Ok(InversePcrOptions {
        edit_start,
        edit_end,
        inserted: String::from_utf8(inserted).map_err(|_| bad("non-ASCII original DNA"))?,
        forward_length,
        reverse_length,
        reverse_tail_bases,
    })
}
fn oe_request(
    input: &PlanInput,
    request: &Request,
) -> Result<crate::cloning::overlap::extension::primers::Request, Error> {
    let Route::OverlapExtension { options } = &request.route else {
        return Err(bad("expected overlap-extension route"));
    };
    if !options.outer_forward_tail.is_empty() || !options.outer_reverse_tail.is_empty() {
        return Err(bad(
            "whole-source substitution-only domestication does not declare external primer tails",
        ));
    }
    let n = input.molecule.sequence.len();
    if request.selection.patches.is_empty() {
        return Err(bad(
            "select actual literal substitutions before primer planning",
        ));
    }
    let mut seen = BTreeSet::new();
    let mut edits = Vec::with_capacity(request.selection.patches.len());
    for patch in &request.selection.patches {
        if patch.position >= n
            || patch.expected.len() != 1
            || patch.replacement.len() != 1
            || !patch.expected.bytes().all(|b| b"ACGT".contains(&b))
            || !patch.replacement.bytes().all(|b| b"ACGT".contains(&b))
            || patch.expected == patch.replacement
            || input.molecule.sequence.as_str().as_bytes()[patch.position]
                != patch.expected.as_bytes()[0]
            || !seen.insert(patch.position)
        {
            return Err(bad(
                "each selected substitution must bind one distinct concrete original base",
            ));
        }
        edits.push(crate::cloning::assembly::edit_intent::Edit {
            role: INPUT_ID.into(),
            start: patch.position,
            end: patch.position + 1,
            expected: patch.expected.clone(),
            replacement: patch.replacement.clone(),
        });
    }
    Ok(crate::cloning::overlap::extension::primers::Request {
        template_strand: crate::chemistry::Strand::Top,
        material_binding: request.origin.as_ref().map(|origin| {
            crate::cloning::overlap::extension::primers::BoundOrigin {
                event: origin.event().clone(),
                origin: origin.clone(),
            }
        }),
        fragment: crate::cloning::overlap::primer_design::Fragment {
            role: INPUT_ID.into(),
            input_id: INPUT_ID.into(),
            start: 0,
            end: n,
            reverse: false,
        },
        edits,
        options: options.clone(),
    })
}
fn inverse_evidence(
    p: mutagenesis::InversePcrPreparation,
) -> (
    crate::planning::plan::DesignPlan,
    MethodEvidence,
    Vec<String>,
) {
    (
        p.plan,
        MethodEvidence::InversePcr {
            options: p.options,
            forward: p.forward,
            reverse: p.reverse,
            forward_location: p.forward_location,
            reverse_location: p.reverse_location,
            retained_template_length: p.retained_template_length,
        },
        p.warnings,
    )
}
fn readiness(method: &MethodEvidence) -> (bool, Vec<String>) {
    match method {
        MethodEvidence::InversePcr{..}=>(true,vec!["Explicit original-source inverse-PCR oligos and domains are declared; inspect actual physical verification separately.".into()]),
        MethodEvidence::InversePcrRounds { rounds } => (true,vec![format!("{} explicitly ordered inverse-PCR rounds use actual preceding closed products. Oligos are declarations from independently intended intermediate genotypes; inspect native amplification/ligation alternatives and final verification separately.",rounds.len())]),
        MethodEvidence::OverlapExtension{search,..}=>{
            let ready=search.len()==3 && search.iter().all(|s|s.unique_geometry);
            let mut notes=vec![if ready{"All three selected native pairs passed the declared filters and established their declared source/extended-substrate geometry. Physical verification, terminal-search limits and isolation remain separate.".into()}else{"Diagnostic proposal only: at least one ranked native pair search did not establish its unique declared PCR geometry. A matching final sequence cannot promote this fallback to a ready primer proposal.".into()}];
            for item in search.iter().filter(|s|s.capped||!s.unique_geometry){
                notes.push(format!("{}: examined {}/{} ranked pairs; rejected {}; capped {}; selected geometry established {}. {}",item.call,item.examined,item.candidates,item.rejected,item.capped,item.unique_geometry,item.notes.join(" ")));
            }
            (ready,notes)
        },
    }
}
/// Freeze the independently checked literal target BEFORE source-bound primer
/// planning, then verify the one composed canonical operation graph once.
pub fn materialize(
    input: &PlanInput,
    request: &Request,
    workspace: Option<&Workspace>,
) -> Result<Prepared, Error> {
    let context_revision = revision(input, request, workspace)?;
    let checked = domestication::check(&input.molecule, &request.design, &request.selection)?;
    let target = checked.target().clone();
    let check = checked.report().clone();
    let (plan, method, mut warnings) = match &request.route {
        Route::InversePcr { .. } => inverse_evidence(mutagenesis::compile_inverse_pcr(
            &input.molecule,
            &target,
            request.origin.as_ref(),
            &options(input, request)?,
        )?),
        Route::InversePcrRounds { rounds } => {
            let (steps, evidence) = crate::cloning::golden_gate::domestication::rounds::compile(
                input,
                &request.selection.patches,
                rounds,
                request.origin.as_ref(),
            )?;
            (crate::planning::plan::DesignPlan {
                standard: None,
                name: format!("Sequential inverse-PCR domestication of {}", input.molecule.name),
                target_revision: target.revision(), inputs: vec![input.clone()], steps,
            }, MethodEvidence::InversePcrRounds { rounds: evidence }, vec!["Every round consumes the preceding real step product. Intermediate intended sequence frames declare primer bases only; they are never invented stocks, saved origins or additional plan inputs. Template removal, isolation, transformation and empirical success are not modeled.".into()])
        }
        Route::OverlapExtension { .. } => {
            let native = oe_request(input, request)?;
            let extra = if input.molecule.topology == Topology::Circular {
                input
                    .molecule
                    .sequence
                    .len()
                    .checked_mul(64)
                    .ok_or(Error::TooLarge)?
            } else {
                0
            };
            let mut budget = crate::cloning::overlap::extension::primers::Budget {
                used: extra,
                limit: 100_000_000,
                materialized: extra,
            };
            let fresh = request
                .origin
                .as_ref()
                .map(|o| o.replay(&input.molecule))
                .transpose()?;
            let template = match &fresh {
                Some(f) => crate::cloning::pcr::source::Template::fresh(
                    &input.molecule,
                    f.tracked(),
                    crate::chemistry::Strand::Top,
                )?,
                None => crate::cloning::pcr::source::Template::reference(&input.molecule)?,
            };
            let (work, material) =
                crate::cloning::overlap::extension::primers::source_work_preflight(
                    template, &native,
                )?;
            if work.checked_add(extra).is_none_or(|n| n > 100_000_000)
                || material
                    .checked_add(extra)
                    .is_none_or(|n| n > 64 * 1024 * 1024)
            {
                return Err(Error::TooLarge);
            }
            let source_frame = crate::cloning::overlap::extension::primers::SourceFrame {
                template_strand: template.strand,
                sequence: template.sequence.clone(),
                topology: template.topology,
                context_revision: template.context_revision()?,
                state_revision: match &fresh {
                    Some(f) => f.tracked().summary()?.state_revision,
                    None => input.revision.clone(),
                },
            };
            let compiled = crate::cloning::overlap::extension::primers::compile_source_with_budget(
                template,
                &native,
                &mut budget,
            )?;
            // This independently source/literal-built reference is a consistency
            // check only. Neither it nor the checked target selects a PCR pair.
            if compiled.expected.sequence != target.sequence {
                return Err(bad(
                    "literal OE declaration differs from independently checked substitutions",
                ));
            }
            let mut steps = vec![];
            let template = append_origin(INPUT_ID, request.origin.as_ref(), &mut steps);
            let outer = compiled.preparation.append(template, &mut steps)?;
            let mut plan = crate::planning::plan::DesignPlan {
                standard: None,
                name: format!("Overlap-extension domestication of {}", input.molecule.name),
                target_revision: target.revision(),
                inputs: vec![input.clone()],
                steps,
            };
            if input.molecule.topology == Topology::Circular {
                plan.steps.push(Operation::TreatEnds {
                    template: outer.clone(),
                    treatment: crate::chemistry::EndTreatment::Phosphorylate,
                });
                plan.steps.push(Operation::Ligate {
                    fragments: vec![crate::planning::plan::OrientedSource {
                        source: Source::Step(plan.steps.len() - 1),
                        reverse: false,
                    }],
                    circular: true,
                    allow_nicked: false,
                    options: Default::default(),
                });
            }
            (
                plan,
                MethodEvidence::OverlapExtension {
                    context_revision: compiled.context_revision,
                    source_frame,
                    preparation: compiled.preparation,
                    designs: compiled.designs,
                    search: compiled.search,
                },
                compiled.warnings,
            )
        }
    };
    if plan.inputs.len() != 1
        || plan.inputs[0].id != INPUT_ID
        || plan.inputs[0].revision != input.revision
    {
        return Err(bad("physical compiler changed the explicit original input"));
    }
    let verification = verification::engine::verify(&plan, &target);
    warnings.push("The target was independently checked and frozen before primer planning. Inspect actual material, search limits, all oligos and intended-path diagnostics. No recoded stock, template clearance, host rescue or experimental success is inferred.".into());
    let (proposal_ready, proposal_findings) = readiness(&method);
    let prepared = Prepared {
        context_revision: context_revision.clone(),
        request: request.clone(),
        source_revision: input.revision.clone(),
        target,
        check,
        plan,
        method,
        verification,
        proposal_ready,
        proposal_findings,
        warnings,
    };
    serde_json::to_writer(&mut Bytes(0), &prepared).map_err(|_| Error::TooLarge)?;
    prepared.validate_binding(&context_revision)?;
    Ok(prepared)
}
impl Prepared {
    /// Bounded readonly source/target/method/graph checks, without genetic search,
    /// primer ranking, PCR simulation, native replay, or status adoption.
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        self.verification.validate_documentary()?;
        self.verification.validate_standard_for_plan(&self.plan)?;
        serde_json::to_writer(&mut Bytes(0), self).map_err(|_| Error::TooLarge)?;
        if self.context_revision != expected || !digest(expected) || self.plan.inputs.len() != 1 {
            return Err(bad("prepared context/input malformed"));
        }
        let input = &self.plan.inputs[0];
        if documentary_revision(input, &self.request)? != expected
            || input.revision != self.source_revision
            || self.plan.target_revision != self.target.revision()
        {
            return Err(bad("prepared original, request or target binding changed"));
        }
        self.target.validate()?;
        if self.target.topology != input.molecule.topology
            || self.target.sequence.len() != input.molecule.sequence.len()
            || self.target.sequence.alphabet() != collagene_model::Alphabet::Dna
        {
            return Err(bad(
                "substitution-only target must retain original topology/length",
            ));
        }
        let patches = self
            .request
            .selection
            .patches
            .iter()
            .map(|p| (p.position, p))
            .collect::<std::collections::BTreeMap<_, _>>();
        for (position, (original, target)) in input
            .molecule
            .sequence
            .as_str()
            .bytes()
            .zip(self.target.sequence.as_str().bytes())
            .enumerate()
        {
            let literal = patches
                .get(&position)
                .map_or(original, |p| p.replacement.as_bytes()[0]);
            if target != literal {
                return Err(bad(
                    "target contains DNA outside exact selected original substitutions",
                ));
            }
        }
        let check_context = domestication::check_revision(
            &input.molecule,
            &self.request.design,
            &self.request.selection,
        )?;
        self.check.validate_binding(&check_context)?;
        if self.check.source_revision != self.source_revision
            || self.check.target_revision != self.target.revision()
            || self.check.selection.candidate_id != self.request.selection.candidate_id
        {
            return Err(bad("checked substitutions changed source/target"));
        }
        match (&self.request.route, &self.method) {
            (
                Route::InversePcr { .. },
                MethodEvidence::InversePcr {
                    options: actual,
                    forward,
                    reverse,
                    forward_location,
                    reverse_location,
                    retained_template_length,
                },
            ) => {
                let desired = options(input, &self.request)?;
                let compiled = mutagenesis::compile_inverse_pcr(
                    &input.molecule,
                    &self.target,
                    self.request.origin.as_ref(),
                    &desired,
                )?;
                if crate::encoding::bounded_json::hash(&self.plan)?
                    != crate::encoding::bounded_json::hash(&compiled.plan)?
                    || crate::encoding::bounded_json::hash(&(
                        actual,
                        forward,
                        reverse,
                        forward_location,
                        reverse_location,
                        retained_template_length,
                    ))? != crate::encoding::bounded_json::hash(&(
                        &compiled.options,
                        &compiled.forward,
                        &compiled.reverse,
                        &compiled.forward_location,
                        &compiled.reverse_location,
                        compiled.retained_template_length,
                    ))?
                {
                    return Err(bad(
                        "inverse-PCR route, primer evidence or canonical graph changed",
                    ));
                }
            }
            (
                Route::OverlapExtension { .. },
                MethodEvidence::OverlapExtension {
                    context_revision,
                    source_frame,
                    preparation,
                    designs,
                    search,
                },
            ) => {
                let native = oe_request(input, &self.request)?;
                let independently_edited =
                    crate::cloning::overlap::extension::primers::validate_documentary_preparation(
                        &input.molecule,
                        &native,
                        source_frame,
                        context_revision,
                        preparation,
                        designs,
                        search,
                    )?;
                if independently_edited.sequence != self.target.sequence {
                    return Err(bad("OE literals differ from checked intended reference"));
                }
                let mut steps = vec![];
                let template = append_origin(INPUT_ID, self.request.origin.as_ref(), &mut steps);
                let outer = preparation.append(template, &mut steps)?;
                if input.molecule.topology == Topology::Circular {
                    steps.push(Operation::TreatEnds {
                        template: outer.clone(),
                        treatment: crate::chemistry::EndTreatment::Phosphorylate,
                    });
                    steps.push(Operation::Ligate {
                        fragments: vec![crate::planning::plan::OrientedSource {
                            source: Source::Step(steps.len() - 1),
                            reverse: false,
                        }],
                        circular: true,
                        allow_nicked: false,
                        options: Default::default(),
                    });
                }
                if self.plan.name
                    != format!("Overlap-extension domestication of {}", input.molecule.name)
                    || crate::encoding::bounded_json::hash(&steps)?
                        != crate::encoding::bounded_json::hash(&self.plan.steps)?
                {
                    return Err(bad(
                        "OE physical branch/extension/outer/closure graph changed",
                    ));
                }
            }
            (
                Route::InversePcrRounds { rounds },
                MethodEvidence::InversePcrRounds { rounds: evidence },
            ) => {
                let (steps, expected_rounds) =
                    crate::cloning::golden_gate::domestication::rounds::compile(
                        input,
                        &self.request.selection.patches,
                        rounds,
                        self.request.origin.as_ref(),
                    )?;
                if self.plan.name
                    != format!(
                        "Sequential inverse-PCR domestication of {}",
                        input.molecule.name
                    )
                    || crate::encoding::bounded_json::hash(&self.plan.steps)?
                        != crate::encoding::bounded_json::hash(&steps)?
                    || crate::encoding::bounded_json::hash(evidence)?
                        != crate::encoding::bounded_json::hash(&expected_rounds)?
                {
                    return Err(bad(
                        "sequential original-coordinate rounds, oligos or preceding-product graph changed",
                    ));
                }
            }
            _ => {
                return Err(bad(
                    "method evidence differs from requested preparation route",
                ));
            }
        }
        if self.verification.plan_revision
            != crate::planning::plan::DesignPlan::revision(&self.plan)
            || self.verification.target_revision != self.target.revision()
            || self.verification.engine_version != verification::ENGINE_VERSION
            || self.verification.planned_steps != self.plan.steps.len()
        {
            return Err(bad(
                "verification belongs to another canonical plan/independent target",
            ));
        }
        let (ready, findings) = readiness(&self.method);
        if self.proposal_ready != ready || self.proposal_findings != findings {
            return Err(bad(
                "primer-proposal readiness differs from actual ranked search evidence",
            ));
        }
        if self.warnings.len() > 64 || self.warnings.iter().any(|w| w.len() > 65536) {
            return Err(Error::TooLarge);
        }
        Ok(())
    }
}
