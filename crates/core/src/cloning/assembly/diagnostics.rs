//! Readonly structural locations from bound current-engine recipe receipts.
//! No molecular operation is rerun; no human finding text determines a location.
use crate::{
    cloning::{
        assembly::recipe::{self as assembly_recipe, PreparedRow},
        ligation::recipe as ligation_recipe,
    },
    planning::plan::{DesignPlan, Operation, Source},
    verification::{
        self,
        evidence::{Status, StepState, VerificationReport},
    },
};
use collagene_model::{Error, Molecule};
use serde::{Deserialize, Serialize};
use std::collections::BTreeSet;
macro_rules! dto{($($t:tt)*)=>{#[derive(Clone,Debug,Serialize,Deserialize)]#[cfg_attr(feature="json-schema",derive(schemars::JsonSchema))]#[serde(deny_unknown_fields)]$($t)*};}
dto! {#[derive(Copy,PartialEq,Eq)] pub enum RowState {InputReference,Reconstructed,Rejected,Withheld,Blocked,NotEvaluated}}
dto! {#[derive(Copy,PartialEq,Eq)] pub enum Stage {Pcr,Digest,Repair,EndTreatment,SourceDeclaration,OverlapAssembly,LigationAssembly,TargetComparison,Annealing,MaterialReconstruction,OverlapExtension}}
dto! {#[derive(Copy,PartialEq,Eq)] pub enum JunctionState {ExactHomology,HomologyMismatch,HomologyUnavailable,EndpointsAvailable,EndpointsUnavailable}}
dto! {pub struct Row {pub index:usize,pub role:String,pub input_id:String,pub reverse:bool,pub state:RowState,pub preparation_step:Option<usize>,pub repair_step:Option<usize>,pub completed_steps:usize,pub total_steps:usize}}
dto! {pub struct Failure {pub step_state:StepState,pub blocked_by:Vec<Source>,pub step:usize,pub stage:Stage,pub row:Option<usize>,pub upstream:Option<usize>,pub downstream:Option<usize>,pub message:String,pub recovery:String}}
dto! {pub struct Junction {pub index:usize,pub upstream:usize,pub downstream:usize,pub closing:bool,pub state:JunctionState,pub message:String,pub recovery:String}}
dto! {pub struct Report {pub context_revision:String,pub plan_revision:String,pub target_revision:String,pub status:Status,pub rows:Vec<Row>,pub failures:Vec<Failure>,pub junctions:Vec<Junction>,pub findings:Vec<crate::verification::checks::CheckFinding>,pub warnings:Vec<String>}}
fn bad(s: &str) -> Error {
    Error::Format(format!("assembly diagnostics: {s}"))
}
fn hex(s: &str) -> bool {
    s.len() == 64
        && s.bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}
fn same(a: &Source, b: &Source) -> bool {
    match (a, b) {
        (Source::Input(a), Source::Input(b)) => a == b,
        (Source::Step(a), Source::Step(b)) => a == b,
        _ => false,
    }
}

fn stage(op: &Operation) -> Result<Stage, Error> {
    Ok(match op {
        Operation::ReconstructRunProduct { .. } | Operation::ReconstructMaterial { .. } => {
            Stage::MaterialReconstruction
        }
        Operation::Pcr { .. } | Operation::PcrWithMismatches { .. } => Stage::Pcr,
        Operation::OverlapExtension { .. } => Stage::OverlapExtension,
        Operation::AnnealOligos { .. } => Stage::Annealing,
        Operation::Digest { .. }
        | Operation::CatalogDigest { .. }
        | Operation::CheckedDigest { .. } => Stage::Digest,
        Operation::RepairEnds { .. } | Operation::ResectThreePrime { .. } => Stage::Repair,
        Operation::TreatEnds { .. } => Stage::EndTreatment,
        Operation::DeclareSourceTermini { .. } | Operation::DeclareSourceModifications { .. } => {
            Stage::SourceDeclaration
        }
        Operation::ProfiledOverlap { .. }
        | Operation::HomologyAssembly { .. }
        | Operation::HomologyAnnealing { .. } => Stage::OverlapAssembly,
        Operation::Ligate { .. } => Stage::LigationAssembly,
        _ => return Err(bad("operation is outside shared recipe compiler scope")),
    })
}
fn guard(
    plan: &DesignPlan,
    report: &VerificationReport,
    context: &str,
    assembly: usize,
) -> Result<(), Error> {
    report.validate_documentary()?;
    if !hex(context)
        || !hex(&plan.target_revision)
        || report.engine_version != verification::ENGINE_VERSION
        || plan.steps.len() != assembly + 1
        || plan.steps.is_empty()
        || plan.steps.len() > 1000
        || plan.inputs.len() > assembly_recipe::MAX_ROWS
        || report.planned_steps != plan.steps.len()
        || report.target_revision != plan.target_revision
        || report.plan_revision != crate::encoding::bounded_json::hash(plan)?
        || report.status == Status::ConsistentUnderAssumptions
            && report.reconstructed_count() != plan.steps.len()
    {
        return Err(bad(
            "receipt context, engine or original-indexed outcome shape is invalid",
        ));
    }
    verification::sources::validate_reference_steps(plan)?;
    assembly_recipe::validate_prepared_inputs(&plan.inputs, &plan.steps)?;
    let mut ids = BTreeSet::new();
    for input in &plan.inputs {
        if input.id.trim().is_empty() || input.id.len() > 128 || !ids.insert(&input.id) {
            return Err(bad("invalid or duplicate original input role"));
        }
        input.molecule.validate()?;
        if input.revision != crate::encoding::bounded_json::hash(&input.molecule)? {
            return Err(bad("stale original input revision"));
        }
    }
    for product in report.steps.iter().filter_map(|step| step.product.as_ref()) {
        product.validate()?;
        if product.metadata.get("plan_revision") != Some(&report.plan_revision)
            || product.metadata.get("engine_version") != Some(&report.engine_version)
        {
            return Err(bad(
                "retained product is not bound to this current original step",
            ));
        }
    }
    if report.findings.len() > 10_000
        || report
            .findings
            .iter()
            .any(|finding| finding.validate().is_err())
    {
        return Err(Error::TooLarge);
    }
    Ok(())
}
// A row's native preparation is Input, one initial Input operation, or annealing
// of its own top original plus a separately bound bottom original; optionally
// PCR->CatalogDigest, then explicit repair. No steps may be skipped or shared.
fn preparation(
    plan: &DesignPlan,
    row: &PreparedRow,
    index: usize,
    cursor: &mut usize,
    assembly: usize,
    owners: &mut [Option<usize>],
) -> Result<Vec<usize>, Error> {
    assembly_recipe::prepared_chain(plan, row, index, cursor, assembly, owners)
}

fn resolve<'a>(
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
        Source::Step(i) => report.product(*i),
    }
}
fn prepared_binding(
    plan: &DesignPlan,
    report: &VerificationReport,
    sources: &[&Source],
    prepared: &[Option<Molecule>],
    complete: bool,
) -> Result<(), Error> {
    if prepared.len() != sources.len() || complete != prepared.iter().all(Option::is_some) {
        return Err(bad("prepared original-row slots or completeness changed"));
    }
    for (source, shown) in sources.iter().zip(prepared) {
        match (resolve(plan, report, source), shown) {
            (Some(actual), Some(shown))
                if crate::encoding::bounded_json::hash(actual)?
                    == crate::encoding::bounded_json::hash(shown)? => {}
            (None, None) => {}
            _ => {
                return Err(bad(
                    "prepared row differs from its original indexed reconstruction outcome",
                ));
            }
        }
    }
    Ok(())
}
fn failure(step: usize, stage: Stage, row: Option<usize>) -> Failure {
    let (message, recovery) = match stage {
        Stage::MaterialReconstruction => (
            "No reusable material was retained from the saved origin reconstruction.",
            "Review the exact saved source event and unchanged originating plan, target and inputs; for a saved digest fragment also review the pinned catalog, fragment index and original digest source event. Reconstruct the complete origin; imported physical summaries do not supply material authority.",
        ),
        Stage::OverlapExtension => (
            "No product was retained for this reciprocal extension step.",
            "Inspect both actual PCR products, explicit orientations, paired terminal homology, 3′ chemistry and search bounds; no intended-target sequence supplies missing material.",
        ),
        Stage::Annealing => (
            "No product was retained for this oligo-annealing step.",
            "Review both independently bound original oligos, explicit terminal assumptions, selected placement and search limits; no stock duplex is substituted.",
        ),
        Stage::Pcr => (
            "No product was retained for this PCR step.",
            "Review this row's actual domains, tails and selected source; rerun preparation to inspect all possible placements.",
        ),
        Stage::Digest => (
            "No product was retained for this digest step.",
            "Review catalog/enzyme pins, actual cut sites, fragment index and source chemical evidence; rerun the declared digest.",
        ),
        Stage::Repair => (
            "No product was retained for this end-repair step.",
            "Review the actual prepared ends and explicit repair method; an unreached repair is not performed chemistry.",
        ),
        Stage::EndTreatment => (
            "No product was retained for this end-treatment step.",
            "Review the bound preparation and explicit end treatment; keep performed and unreached steps separate.",
        ),
        Stage::SourceDeclaration => (
            "No product was retained for this source declaration.",
            "Review the intact original reference and explicit terminal assumptions; imported physical summaries do not supply authority.",
        ),
        Stage::OverlapAssembly => (
            "The overlap assembly step retained no final product.",
            "Inspect typed homology locations and the actual prepared substrates. Keep the declared inputs and rerun full reconstruction after edits.",
        ),
        Stage::LigationAssembly => (
            "The ligation step retained no final product.",
            "Inspect actual oriented endpoint summaries and the explicit nick policy. This receipt does not identify one incompatible junction.",
        ),
        Stage::TargetComparison => (
            "The designated final reference does not satisfy the independent-target comparison.",
            "Compare the actual final reference with the separately declared target and topology/origin policy; do not replace the target with the simulator's answer.",
        ),
    };
    Failure {
        step_state: StepState::Rejected,
        blocked_by: vec![],
        step,
        stage,
        row,
        upstream: None,
        downstream: None,
        message: message.into(),
        recovery: recovery.into(),
    }
}
fn finish(
    context: &str,
    plan: &DesignPlan,
    verification: &VerificationReport,
    prepared_rows: &[&PreparedRow],
    chains: Vec<Vec<usize>>,
    owners: &[Option<usize>],
    assembly: usize,
    junctions: Vec<Junction>,
) -> Result<Report, Error> {
    if owners.iter().any(Option::is_none) {
        return Err(bad(
            "compiled receipt contains an unassigned preparation/treatment step",
        ));
    }
    let mut rows = Vec::with_capacity(chains.len());
    for (index, (row, chain)) in prepared_rows.iter().zip(chains).enumerate() {
        let completed_steps = chain
            .iter()
            .filter(|index| verification.product(**index).is_some())
            .count();
        let state = if chain.is_empty() {
            RowState::InputReference
        } else if completed_steps == chain.len() {
            RowState::Reconstructed
        } else if chain
            .iter()
            .any(|index| verification.steps[*index].state == StepState::Rejected)
        {
            RowState::Rejected
        } else if chain
            .iter()
            .any(|index| verification.steps[*index].state == StepState::Withheld)
        {
            RowState::Withheld
        } else if chain
            .iter()
            .any(|index| verification.steps[*index].state == StepState::Blocked)
        {
            RowState::Blocked
        } else {
            RowState::NotEvaluated
        };
        rows.push(Row {
            index,
            role: row.role.clone(),
            input_id: row.input_id.clone(),
            reverse: row.reverse,
            state,
            preparation_step: row.preparation_step,
            repair_step: row.repair_step,
            completed_steps,
            total_steps: chain.len(),
        });
    }
    let mut failures = Vec::new();
    for outcome in &verification.steps {
        if outcome.state != StepState::Reconstructed {
            let mut entry = failure(
                outcome.step_index,
                stage(&plan.steps[outcome.step_index])?,
                owners.get(outcome.step_index).copied().flatten(),
            );
            entry.step_state = outcome.state;
            entry.blocked_by = outcome.blocked_by.clone();
            failures.push(entry);
        }
    }
    if verification.findings.iter().any(|finding| {
        finding.rule == crate::verification::checks::RuleId::TargetAgreement
            && finding.verdict == crate::verification::checks::CheckVerdict::Fail
    }) {
        let mut entry = failure(assembly, Stage::TargetComparison, None);
        entry.step_state = verification.steps[assembly].state;
        failures.push(entry);
    }
    let mut warnings=vec!["Readonly locations from a bound current-engine receipt, not a new scientific reconstruction or experimental verdict. Reconstructed means that original step retained a reference product, not proof of valid chemistry or successful preparation.".into()];
    if verification.status == Status::Incomplete
        && verification.reconstructed_count() == plan.steps.len()
    {
        warnings.push("All operation products are retained while the overall result remains incomplete under explicit named checks. No missing operation or local cause is inferred from finding text.".into());
    }
    let result = Report {
        context_revision: context.into(),
        plan_revision: verification.plan_revision.clone(),
        target_revision: verification.target_revision.clone(),
        status: verification.status,
        rows,
        failures,
        junctions,
        findings: verification.findings.clone(),
        warnings,
    };
    crate::encoding::bounded_json::hash(&result)?;
    Ok(result)
}
fn span_matches(
    m: &Molecule,
    span: &crate::cloning::overlap::profiles::ReferenceSpan,
    text: &str,
) -> bool {
    if span.start >= span.end || span.end > m.sequence.len() || span.end - span.start != text.len()
    {
        return false;
    }
    text.bytes().enumerate().all(|(i, b)| {
        let actual = if span.reverse {
            collagene_model::complement(
                m.sequence.as_str().as_bytes()[span.end - 1 - i],
                collagene_model::Alphabet::Dna,
            )
        } else {
            m.sequence.as_str().as_bytes()[span.start + i]
        };
        actual == b
    })
}
pub fn overlap(p: &assembly_recipe::Preview) -> Result<Report, Error> {
    p.validate_binding(&p.context_revision)?;
    let plan = &p.compiled.plan;
    let assembly = p.compiled.assembly_step;
    guard(plan, &p.verification, &p.context_revision, assembly)?;
    let mut cursor = 0usize;
    let mut owners = vec![None; assembly];
    let mut chains = Vec::new();
    let mut roles = BTreeSet::new();
    for (i, row) in p.compiled.rows.iter().enumerate() {
        if !roles.insert(&row.role) {
            return Err(bad("duplicate row role"));
        }
        chains.push(preparation(
            plan,
            row,
            i,
            &mut cursor,
            assembly,
            &mut owners,
        )?);
    }
    if cursor != assembly {
        return Err(bad(
            "unexpected operations between preparation and assembly",
        ));
    }
    prepared_binding(
        plan,
        &p.verification,
        &p.compiled
            .rows
            .iter()
            .map(|r| &r.prepared_source)
            .collect::<Vec<_>>(),
        &p.prepared,
        p.preparation_complete,
    )?;
    if let Some(profile) = &p.profile {
        if profile.fragment_lengths
            != p.prepared
                .iter()
                .map(|m| m.as_ref().map(|m| m.sequence.len()))
                .collect::<Option<Vec<_>>>()
                .ok_or_else(|| bad("profile requires every original prepared row"))?
        {
            return Err(bad("profile fragment lengths do not match prepared rows"));
        }
        for j in &profile.junctions {
            if j.upstream_span.as_ref().is_some_and(|span| {
                !p.prepared[j.upstream]
                    .as_ref()
                    .is_some_and(|model| span_matches(model, span, &j.homology))
            }) || j.downstream_span.as_ref().is_some_and(|span| {
                !p.prepared[j.downstream]
                    .as_ref()
                    .is_some_and(|model| span_matches(model, span, &j.downstream_homology))
            }) {
                return Err(bad(
                    "typed homology strings do not match their actual prepared source spans",
                ));
            }
        }
        if p.verification.final_product().is_some() && !profile.geometry_valid {
            return Err(bad(
                "retained final overlap product contradicts typed mismatch evidence",
            ));
        }
    }
    let (fragments, circular, overlaps) = match &plan.steps[assembly] {
        Operation::ProfiledOverlap {
            fragments,
            circular,
            overlaps,
            ..
        } => (fragments, *circular, overlaps),
        Operation::HomologyAssembly {
            fragments, policy, ..
        } => (
            fragments,
            policy.topology == collagene_model::Topology::Circular,
            &policy.overlaps,
        ),
        Operation::HomologyAnnealing {
            fragments, policy, ..
        } => (
            fragments,
            policy.topology == collagene_model::Topology::Circular,
            &policy.overlaps,
        ),
        _ => return Err(bad("expected final overlap operation")),
    };
    let mut junctions = Vec::with_capacity(overlaps.len());
    for index in 0..overlaps.len() {
        let upstream = index;
        let downstream = (index + 1) % fragments.len();
        let typed = p.profile.as_ref().and_then(|r| r.junctions.get(index));
        let (state, message, recovery) = match typed {
            Some(j)
                if j.upstream_span.is_some() && j.downstream_span.is_some() && j.exact_match =>
            {
                (
                    JunctionState::ExactHomology,
                    "The bound profile reports exact oriented end homology.",
                    "Retain its occurrence/completeness and physical-state limitations; sequence agreement does not establish successful assembly.",
                )
            }
            Some(j) if j.upstream_span.is_some() && j.downstream_span.is_some() => (
                JunctionState::HomologyMismatch,
                "The bound profile reports unequal homology at this row pair.",
                "Review this upstream suffix, downstream prefix, orientation and declared overlap length; change explicit choices before rerunning.",
            ),
            _ => (
                JunctionState::HomologyUnavailable,
                "No bound homology location is available for this junction.",
                "Prepare the required rows and inspect a fresh profile; no coordinates or matching bases are invented.",
            ),
        };
        junctions.push(Junction {
            index,
            upstream,
            downstream,
            closing: circular && downstream == 0,
            state,
            message: message.into(),
            recovery: recovery.into(),
        });
    }
    finish(
        &p.context_revision,
        plan,
        &p.verification,
        &p.compiled.rows.iter().collect::<Vec<_>>(),
        chains,
        &owners,
        assembly,
        junctions,
    )
}
pub fn ligation(p: &ligation_recipe::Preview) -> Result<Report, Error> {
    p.validate_binding(&p.context_revision)?;
    let plan = &p.compiled.plan;
    let assembly = p.compiled.ligation_step;
    guard(plan, &p.verification, &p.context_revision, assembly)?;
    let mut cursor = 0usize;
    let mut owners = vec![None; assembly];
    let mut chains = Vec::new();
    let mut roles = BTreeSet::new();
    for (i, row) in p.compiled.rows.iter().enumerate() {
        if !roles.insert(&row.preparation.role) {
            return Err(bad("duplicate row role"));
        }
        chains.push(preparation(
            plan,
            &row.preparation,
            i,
            &mut cursor,
            assembly,
            &mut owners,
        )?);
    }
    // The shared compiler appends treatments only after ALL preparation chains.
    for (i, row) in p.compiled.rows.iter().enumerate() {
        if let Some(step) = row.treatment_step {
            if step != cursor
                || step >= assembly
                || !matches!(plan.steps.get(step),Some(Operation::TreatEnds{template,..})if same(template,&row.preparation.prepared_source))
                || !same(&row.prepared_source, &Source::Step(step))
            {
                return Err(bad("end-treatment chain is not in compiler order"));
            }
            if owners[step].replace(i).is_some() {
                return Err(bad("treatment operation assigned twice"));
            }
            chains[i].push(step);
            cursor += 1;
        } else if !same(&row.prepared_source, &row.preparation.prepared_source) {
            return Err(bad("row final source changed without treatment"));
        }
    }
    if cursor != assembly {
        return Err(bad("unassigned operations before ligation"));
    }
    prepared_binding(
        plan,
        &p.verification,
        &p.compiled
            .rows
            .iter()
            .map(|r| &r.prepared_source)
            .collect::<Vec<_>>(),
        &p.prepared,
        p.preparation_complete,
    )?;
    let junctions=p.junctions.iter().map(|j|{let ready=j.upstream_end.is_some()&&j.downstream_end.is_some();Junction{index:j.index,upstream:j.upstream,downstream:j.downstream,closing:j.closing,state:if ready{JunctionState::EndpointsAvailable}else{JunctionState::EndpointsUnavailable},message:if ready{"Both actual oriented endpoint summaries are retained; compatibility has not been re-evaluated here."}else{"One or both oriented endpoint summaries are unavailable."}.into(),recovery:if ready{"Review these ends with the fresh ligation findings and explicit nick policy; do not infer a failed pair from generic text."}else{"Resolve the missing preparation/treatment or unsupported physical state, then rerun the complete recipe."}.into()}}).collect();
    finish(
        &p.context_revision,
        plan,
        &p.verification,
        &p.compiled
            .rows
            .iter()
            .map(|r| &r.preparation)
            .collect::<Vec<_>>(),
        chains,
        &owners,
        assembly,
        junctions,
    )
}
