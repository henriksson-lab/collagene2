//! Original-template physical preparation for explicitly checked domestication.
//! A reference-only recoded target never becomes an invented PCR template.
use crate::{
    Error, Molecule, Topology, assembly_recipe, domestication,
    mutagenesis::{self, InversePcrOptions},
    verification::{self, Operation, PlanInput, Source},
};
use serde::{Deserialize, Serialize};
use std::{collections::BTreeSet, io};

pub const ENGINE_VERSION: &str = "domestication-preparation-2";
/// Current documentary review receipt, not executable material authority.
pub const REVIEW_EVIDENCE_KIND: &str = "domestication-intent-v2";
pub const REVIEW_EVIDENCE_VERSION: u16 = 2;
pub const INPUT_ID: &str = "template";
const MAX_SOURCE_BASES: usize = 8 << 20;
const MAX_WIRE_BYTES: usize = 32 << 20;
macro_rules! dto {($($tokens:tt)*)=>{#[derive(Clone,Debug,Serialize,Deserialize)]#[cfg_attr(feature="json-schema",derive(schemars::JsonSchema))]#[serde(deny_unknown_fields)]$($tokens)*};}
dto! {pub enum Route {
    InversePcr {edit_start:usize,edit_end:usize,forward_length:usize,reverse_length:usize,reverse_tail_bases:usize},
    OverlapExtension {options:crate::overlap_extension_primer_design::Options},
}}
dto! {pub struct Request {pub design:domestication::Request,pub selection:domestication::Selection,pub route:Route}}
dto! {pub enum MethodEvidence {
    InversePcr {options:InversePcrOptions,forward:crate::analysis::Primer,reverse:crate::analysis::Primer,forward_location:crate::Location,reverse_location:crate::Location,retained_template_length:usize},
    OverlapExtension {context_revision:String,request:crate::overlap_extension_primer_design::Request,preparation:crate::overlap_extension_pcr::Request,#[cfg_attr(feature="json-schema",schemars(with="serde_json::Value"))]designs:Vec<crate::primers::DesignReport>,search:Vec<crate::overlap_extension_primer_design::PairSearch>},
}}
dto! {pub struct Prepared {
    pub context_revision:String,pub request:Request,pub source_revision:String,pub target:Molecule,
    pub check:domestication::CheckReport,pub plan:verification::DesignPlan,pub method:MethodEvidence,
    pub verification:verification::VerificationReport,pub proposal_ready:bool,pub proposal_findings:Vec<String>,pub warnings:Vec<String>,
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
pub fn preflight(input: &PlanInput, request: &Request) -> Result<(), Error> {
    if input.id != INPUT_ID
        || input.revision != request.design.source_revision
        || !digest(&input.revision)
    {
        return Err(bad(
            "the explicit template input and design must bind the same original source revision",
        ));
    }
    crate::model::bounds::model_preflight(&[&input.molecule])?;
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
    if crate::chemistry::has_unmapped_physical_evidence(&input.molecule) {
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
            crate::overlap_extension_primer_design::preflight(&input.molecule, &native)?;
            let (work, material) =
                crate::overlap_extension_primer_design::work_preflight(&input.molecule, &native)?;
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
    Ok(())
}
pub fn revision(input: &PlanInput, request: &Request) -> Result<String, Error> {
    preflight(input, request)?;
    assembly_recipe::hash(&(ENGINE_VERSION, input, request))
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
) -> Result<crate::overlap_extension_primer_design::Request, Error> {
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
        edits.push(crate::assembly_edit_intent::Edit {
            role: INPUT_ID.into(),
            start: patch.position,
            end: patch.position + 1,
            expected: patch.expected.clone(),
            replacement: patch.replacement.clone(),
        });
    }
    Ok(crate::overlap_extension_primer_design::Request {
        fragment: crate::overlap_primer_design::Fragment {
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
) -> (verification::DesignPlan, MethodEvidence, Vec<String>) {
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
pub fn materialize(input: &PlanInput, request: &Request) -> Result<Prepared, Error> {
    let context_revision = revision(input, request)?;
    let checked = domestication::check(&input.molecule, &request.design, &request.selection)?;
    let target = checked.target().clone();
    let check = checked.report().clone();
    let (plan, method, mut warnings) = match &request.route {
        Route::InversePcr { .. } => inverse_evidence(mutagenesis::compile_inverse_pcr(
            &input.molecule,
            &target,
            None,
            &options(input, request)?,
        )?),
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
            let mut budget = crate::overlap_extension_primer_design::Budget {
                used: extra,
                limit: 100_000_000,
                materialized: extra,
            };
            let compiled = crate::overlap_extension_primer_design::compile_with_budget(
                &input.molecule,
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
            let mut plan = compiled.plan;
            plan.name = format!("Overlap-extension domestication of {}", input.molecule.name);
            plan.target_revision = target.revision();
            if input.molecule.topology == Topology::Circular {
                plan.steps.push(Operation::TreatEnds {
                    template: Source::Step(3),
                    treatment: crate::chemistry::EndTreatment::Phosphorylate,
                });
                plan.steps.push(Operation::Ligate {
                    fragments: vec![verification::OrientedSource {
                        source: Source::Step(4),
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
                    request: compiled.request,
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
    let verification = verification::verify(&plan, &target);
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
        serde_json::to_writer(&mut Bytes(0), self).map_err(|_| Error::TooLarge)?;
        if self.context_revision != expected || !digest(expected) || self.plan.inputs.len() != 1 {
            return Err(bad("prepared context/input malformed"));
        }
        let input = &self.plan.inputs[0];
        if revision(input, &self.request)? != expected
            || input.revision != self.source_revision
            || self.plan.target_revision != self.target.revision()
        {
            return Err(bad("prepared original, request or target binding changed"));
        }
        self.target.validate()?;
        if self.target.topology != input.molecule.topology
            || self.target.sequence.len() != input.molecule.sequence.len()
            || self.target.sequence.alphabet() != crate::Alphabet::Dna
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
                    None,
                    &desired,
                )?;
                if assembly_recipe::hash(&self.plan)? != assembly_recipe::hash(&compiled.plan)?
                    || assembly_recipe::hash(&(
                        actual,
                        forward,
                        reverse,
                        forward_location,
                        reverse_location,
                        retained_template_length,
                    ))? != assembly_recipe::hash(&(
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
                    request,
                    preparation,
                    designs,
                    search,
                },
            ) => {
                let native = oe_request(input, &self.request)?;
                if assembly_recipe::hash(request)? != assembly_recipe::hash(&native)? {
                    return Err(bad("OE request differs from literal original patches"));
                }
                let independently_edited =
                    crate::overlap_extension_primer_design::validate_preparation(
                        &input.molecule,
                        request,
                        context_revision,
                        preparation,
                        designs,
                        search,
                    )?;
                if independently_edited.sequence != self.target.sequence {
                    return Err(bad("OE literals differ from checked intended reference"));
                }
                let mut steps = vec![];
                preparation.append(Source::Input(INPUT_ID.into()), &mut steps)?;
                if input.molecule.topology == Topology::Circular {
                    steps.push(Operation::TreatEnds {
                        template: Source::Step(3),
                        treatment: crate::chemistry::EndTreatment::Phosphorylate,
                    });
                    steps.push(Operation::Ligate {
                        fragments: vec![verification::OrientedSource {
                            source: Source::Step(4),
                            reverse: false,
                        }],
                        circular: true,
                        allow_nicked: false,
                        options: Default::default(),
                    });
                }
                if self.plan.name
                    != format!("Overlap-extension domestication of {}", input.molecule.name)
                    || assembly_recipe::hash(&steps)? != assembly_recipe::hash(&self.plan.steps)?
                {
                    return Err(bad(
                        "OE physical branch/extension/outer/closure graph changed",
                    ));
                }
            }
            _ => {
                return Err(bad(
                    "method evidence differs from requested preparation route",
                ));
            }
        }
        if self.verification.plan_revision != verification::plan_revision(&self.plan)
            || self.verification.target_revision != self.target.revision()
            || self.verification.engine_version != verification::ENGINE_VERSION
            || self.verification.products.len() > self.plan.steps.len()
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
