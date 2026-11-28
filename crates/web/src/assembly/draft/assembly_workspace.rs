//! Persistent declared assembly slots. Frozen intent/scientific evidence is
//! component-private and is never restored from a serialized draft.
use crate::assembly::draft::model::CapturedSource;
use crate::assembly::draft::model::Draft;
use crate::assembly::draft::model::IntentIdentity;
use crate::assembly::review::assembly_review::IntendedReference;
use crate::cloning::preparation::cloning_methods::CloningMethodsPanel;
use crate::cloning::preparation::cloning_methods_model::MethodsState;
use crate::cloning::preparation::cloning_methods_model::Mode;
use crate::cloning::preparation::cloning_methods_model::OverlapInput;
use crate::cloning::preparation::cloning_methods_model::OverlapPreparation;
use crate::cloning::preparation::review::TargetContext;
use crate::platform::worker::CachedWorker;
use crate::primers::primer_picker::ResolveRequest;
use crate::sources::identity::SourceInfo;
use collagene_core::{
    cloning::{
        assembly::{
            edit_intent as assembly_edit_intent, intent as assembly_intent,
            junction_edit as assembly_junction_edit, selective_edit as selective_junction_edit,
            substitution_intent as assembly_substitution_intent,
        },
        golden_gate::junction_edit as golden_gate_junction_edit,
        pcr::reuse as primer_reuse,
    },
    documents::transfer::Candidate,
    material::lineage::EventBinding,
    planning::plan::PlanInput,
    protocol::{
        request::Request,
        schema::{AnalysisJob, AnalysisResult},
    },
};
use collagene_model::{Alphabet, Molecule, Topology};
use leptos::prelude::*;
use std::{
    rc::Rc,
    sync::{
        Arc,
        atomic::{AtomicBool, Ordering},
    },
};
use wasm_bindgen_futures::spawn_local;
#[derive(Clone)]
pub struct SelectionNavigation {
    pub source: EventBinding,
    pub regions: Vec<(usize, usize, bool)>,
    pub description: String,
}
#[derive(Clone)]
pub struct ReviewRequest {
    pub candidate: Candidate,
    pub sources: Vec<EventBinding>,
    pub intent: IntendedReference,
}
struct Frozen {
    identity: IntentIdentity,
    epoch: u64,
    intent: IntendedReference,
    sources: Vec<Arc<CapturedSource>>,
    input_ids: Vec<String>,
    target_revision: String,
}

fn method_overview(methods: &MethodsState) -> String {
    match methods.mode {
        Mode::Overlap => format!(
            "{} overlap assembly",
            match methods.method {
                collagene_core::cloning::overlap::profiles::OverlapMethod::Gibson => "Gibson",
                collagene_core::cloning::overlap::profiles::OverlapMethod::HiFi => "HiFi",
                collagene_core::cloning::overlap::profiles::OverlapMethod::InFusion => "In-Fusion",
            }
        ),
        Mode::Ligation => "Ordered ligation".into(),
        Mode::GoldenGate => format!(
            "{} · {:?}",
            if methods.golden_gate_pool.is_some() {
                "Finite Golden Gate reaction"
            } else {
                "Ordered Type IIS restriction/ligation"
            },
            methods.golden_gate.enzyme
        ),
        Mode::MultiSiteReference => "Multi-parent BP/LR reference exchange".into(),
        Mode::Tailing => "TA / GC tailing and ligation".into(),
    }
}
fn preparation_overview(methods: &MethodsState, index: usize) -> String {
    if methods.mode == Mode::MultiSiteReference {
        return "Whole original reference parent; contribution region defines intent only".into();
    }
    if methods.mode == Mode::GoldenGate {
        if let Some(pool) = methods
            .golden_gate_pool
            .as_ref()
            .filter(|pool| pool.reaction)
        {
            match pool.choices.get(index).copied() {
                Some(1) => return "Whole original source in the declared reaction pool".into(),
                Some(2) => {
                    return "Complete PCR amplicon before digestion in the declared reaction pool"
                        .into();
                }
                _ => {}
            }
        }
    }
    methods
        .rows
        .get(index)
        .map(OverlapInput::summary)
        .unwrap_or_else(|| "Preparation has not been set".into())
}

#[component]
pub fn AssemblyWorkspace(
    web_mcp: RwSignal<Option<crate::ai::web_mcp_assembly::bridge::Endpoint>>,
    draft: RwSignal<Draft>,
    previous_draft: RwSignal<Option<Arc<Draft>>>,
    epoch: RwSignal<u64>,
    snapshot: Callback<(), Result<collagene_core::documents::workspace::Workspace, String>>,
    workspace_busy: RwSignal<bool>,
    apply_saved_attachments: Callback<
        (
            collagene_core::documents::workspace::Workspace,
            collagene_core::documents::workspace::Workspace,
        ),
        Result<(), String>,
    >,
    sources: Memo<Vec<SourceInfo>>,
    active_source: Memo<SourceInfo>,
    selected_region: Callback<(), Option<(usize, usize)>>,
    get_source: Callback<EventBinding, Result<Molecule, String>>,
    on_navigate: Callback<SelectionNavigation>,
    on_review: Callback<ReviewRequest, Result<(), String>>,
    on_finalize: Callback<ReviewRequest, Result<(), String>>,
    on_products: Callback<Vec<Molecule>, Result<(), String>>,
    primer_sources: Memo<Vec<SourceInfo>>,
    list_primers: Callback<EventBinding, Result<primer_reuse::List, String>>,
    resolve_primer: Callback<ResolveRequest, Result<primer_reuse::Resolved, String>>,
) -> impl IntoView {
    let strategy_chosen =
        RwSignal::new(draft.with_untracked(|d| d.slots.iter().any(|s| s.captured.is_some())));
    let identity = Memo::new(move |_| draft.with(Draft::identity));
    let methods = RwSignal::new(draft.with_untracked(|d| d.methods.clone()));
    Effect::new(move |_| {
        let value = draft.with(|d| d.methods.clone());
        if methods.get_untracked() != value {
            methods.set(value);
        }
    });
    Effect::new(move |_| {
        let value = methods.get();
        if draft.with_untracked(|d| d.methods != value) {
            draft.update(|d| {
                d.methods = value;
                d.sync_rows();
            });
        }
    });
    let frozen = RwSignal::new(None::<Arc<Frozen>>);
    let review_pin = RwSignal::new(None::<crate::ai::web_mcp_assembly::Pin>);
    let source_receipt = RwSignal::new(
        None::<(
            u64,
            Arc<crate::assembly::edit::assembly_edit_model::SourceReview>,
        )>,
    );
    let source_review = Memo::new(move |_| {
        source_receipt
            .get()
            .filter(|(held, _)| *held == epoch.get())
            .map(|(_, review)| review)
    });
    let undo_junction = RwSignal::new(None::<(Draft, Draft)>);
    let busy = RwSignal::new(false);
    let status = RwSignal::new(String::new());
    let generation = RwSignal::new(0u32);
    let worker = StoredValue::new_local(None::<Rc<CachedWorker>>);
    let alive = Arc::new(AtomicBool::new(true));
    let cleanup = alive.clone();
    on_cleanup(move || {
        cleanup.store(false, Ordering::Relaxed);
        if let Some(w) = worker.get_value() {
            w.close();
        }
        worker.set_value(None);
    });
    Effect::new(move |_| {
        let (current, current_epoch) = (identity.get(), epoch.get());
        if frozen.with_untracked(|value| {
            value
                .as_ref()
                .is_some_and(|f| f.identity == current && f.epoch == current_epoch)
        }) {
            return;
        }
        generation.update(|g| *g = g.wrapping_add(1));
        if let Some(w) = worker.get_value() {
            w.close();
        }
        worker.set_value(None);
        frozen.set(None);
        busy.set(false);
        status.set(String::new());
    });
    let pin = Callback::new(move |kind: u8| {
        let source = active_source.get_untracked();
        let n = source.length;
        if source.alphabet != Alphabet::Dna || n == 0 {
            status.set("Choose a nonempty DNA source in the sequence viewer.".into());
            return;
        }
        let (start, end) = if kind == 0 {
            (0, n)
        } else {
            let Some((a, b)) = selected_region.run(()) else {
                status
                    .set("Select a contiguous source region in the sequence viewer first.".into());
                return;
            };
            if a >= b || b > n {
                status.set("Selection is outside the active source. Select it again.".into());
                return;
            }
            if kind == 2 {
                if source.topology != Topology::Circular {
                    status.set("Keeping outside a selection requires a circular source; choose an explicit region for linear DNA.".into());
                    return;
                }
                let length = n - (b - a);
                if length == 0 {
                    status.set("The outside of a whole-source selection is empty.".into());
                    return;
                }
                let start = b % n;
                let Some(end) = start.checked_add(length) else {
                    status.set("Selection range overflow.".into());
                    return;
                };
                (start, end)
            } else {
                (a, b)
            }
        };
        let molecule = match get_source.run(source.binding.clone()) {
            Ok(v) => v,
            Err(e) => {
                status.set(e);
                return;
            }
        };
        if active_source.get_untracked().binding != source.binding {
            status.set("The active source changed before capture. Select the region again.".into());
            return;
        }
        let mut next = draft.get_untracked();
        let index = next.active_slot;
        if index >= next.slots.len() {
            status.set("Choose an assembly slot first.".into());
            return;
        }
        if let Err(error) = next.clear_slot_substitutions(index) {
            status.set(error);
            return;
        }
        next.slots[index].material_origin = None;
        next.slots[index].captured = Some(Arc::new(CapturedSource { source, molecule }));
        next.slots[index].anneal_bottom = None;
        let orientation = next.methods.rows[index].reverse;
        next.slots[index].selector = assembly_intent::Selector::ReferenceRange {
            start,
            end,
            reverse: orientation,
        };
        next.slots[index].digest_geometry = None;
        next.junction_edits.clear();
        let orientation = next.methods.rows[index].reverse;
        next.methods.rows[index] = OverlapInput {
            preparation: OverlapPreparation::ExactPcr,
            reverse: orientation,
            ..Default::default()
        };
        next.sync_rows();
        match next.validate() {
            Ok(()) => {
                draft.set(next);
                status.set("Source region pinned to this slot. Review the intended construct, then choose or design its PCR primers/preparation.".into());
            }
            Err(e) => status.set(e),
        }
    });
    let pin_bottom = Callback::new(move |()| {
        let source = active_source.get_untracked();
        if source.alphabet != Alphabet::Dna
            || source.topology != Topology::Linear
            || source.length == 0
            || source.length > collagene_core::chemistry::MAX_ANNEAL_OLIGO_BASES
        {
            status.set(
                "Choose the actual bottom DNA oligo as a linear 5′→3′ source (at most10000 bases)."
                    .into(),
            );
            return;
        }
        let molecule = match get_source.run(source.binding.clone()) {
            Ok(m) => m,
            Err(error) => {
                status.set(error);
                return;
            }
        };
        if active_source.get_untracked().binding != source.binding {
            status.set("The viewer source changed; capture the bottom oligo again.".into());
            return;
        }
        let mut next = draft.get_untracked();
        let index = next.active_slot;
        if next.methods.mode != Mode::Ligation
            || next.alternatives.is_some()
            || index >= next.slots.len()
            || next.slots[index].captured.is_none()
        {
            status.set("Use ordered ligation and pin this fragment’s top oligo before capturing its bottom strand.".into());
            return;
        }
        next.slots[index].anneal_bottom = Some(Arc::new(CapturedSource { source, molecule }));
        next.methods.rows[index].preparation = OverlapPreparation::AnnealedOligos;
        next.methods.rows[index].anneal_candidate = None;
        next.methods.rows[index].repair = None;
        next.sync_rows();
        match next.validate() {
            Ok(()) => {
                draft.set(next);
                status.set("Bottom oligo pinned. Review the independently selected top-strand contribution, then inspect and choose an annealing geometry in preparation.".into());
            }
            Err(error) => status.set(error),
        }
    });
    let review_intent = Callback::new(move |()| {
        let declared = draft.get_untracked();
        let state = identity.get_untracked();
        let bound_epoch = epoch.get_untracked();
        review_pin.set(crate::ai::web_mcp_assembly::pin(&declared, bound_epoch).ok());
        let frame = if declared.slots.iter().any(|slot| {
            matches!(
                &slot.selector,
                assembly_intent::Selector::MaterialStrandRange { .. }
            )
        }) {
            match snapshot.run(()) {
                Ok(frame) => Some(frame),
                Err(error) => {
                    status.set(error);
                    return;
                }
            }
        } else {
            None
        };
        let (inputs, request, captures) = match declared.inputs_request(frame.as_ref()) {
            Ok(v) => v,
            Err(e) => {
                status.set(e);
                return;
            }
        };
        let revision = match if !declared.edits.is_empty() {
            assembly_edit_intent::intent_revision(
                &inputs,
                &declared.edit_request(request.clone()),
                frame.as_ref(),
            )
        } else if !declared.substitutions.is_empty() {
            assembly_substitution_intent::intent_revision(
                &inputs,
                &assembly_substitution_intent::Request {
                    base: request.clone(),
                    junction_edits: declared.junction_edits.clone(),
                    substitutions: declared.substitutions.clone(),
                },
                frame.as_ref(),
            )
        } else if declared.junction_edits.is_empty() {
            assembly_intent::intent_revision(&inputs, &request, frame.as_ref())
        } else {
            assembly_junction_edit::intent_revision(
                &inputs,
                &request,
                &declared.junction_edits,
                frame.as_ref(),
            )
        } {
            Ok(value) => value,
            Err(error) => {
                status.set(error.to_string());
                return;
            }
        };
        let input_ids = inputs
            .iter()
            .map(|input| input.id.clone())
            .collect::<Vec<_>>();
        let job = if !declared.edits.is_empty() {
            AnalysisJob::BuildContributionEditIntent {
                request: declared.edit_request(request),
                inputs,
                workspace: frame.clone(),
            }
        } else if !declared.substitutions.is_empty() {
            AnalysisJob::BuildSubstitutedAssemblyIntent {
                inputs,
                workspace: frame.clone(),
                request: assembly_substitution_intent::Request {
                    base: request,
                    junction_edits: declared.junction_edits.clone(),
                    substitutions: declared.substitutions.clone(),
                },
            }
        } else if declared.junction_edits.is_empty() {
            AnalysisJob::BuildAssemblyIntent {
                inputs,
                request,
                workspace: frame.clone(),
            }
        } else {
            AnalysisJob::BuildEditedAssemblyIntent {
                inputs,
                base: request,
                edits: declared.junction_edits.clone(),
                workspace: frame,
            }
        };
        let current = match CachedWorker::new() {
            Ok(v) => Rc::new(v),
            Err(e) => {
                status.set(e);
                return;
            }
        };
        if let Some(old) = worker.get_value() {
            old.close();
        }
        worker.set_value(Some(current.clone()));
        let token = generation.get_untracked().wrapping_add(1);
        generation.set(token);
        busy.set(true);
        frozen.set(None);
        status.set("Building the intended reference from your pinned source regions…".into());
        let alive = alive.clone();
        spawn_local(async move {
            let response = current
                .ask(Request {
                    id: token,
                    revision: revision.clone(),
                    job,
                })
                .await;
            current.close();
            if !alive.load(Ordering::Relaxed)
                || generation.get_untracked() != token
                || identity.get_untracked() != state
                || epoch.get_untracked() != bound_epoch
            {
                return;
            }
            worker.set_value(None);
            busy.set(false);
            match response {
                Ok(response) if response.id == token && response.revision == revision => {
                    let result = match response.result {
                        AnalysisResult::AssemblyIntent(intent) => intent
                            .validate_binding(&revision)
                            .map(|()| IntendedReference::Regions(intent)),
                        AnalysisResult::EditedAssemblyIntent(intent) => intent
                            .validate_binding(&revision)
                            .map(|()| IntendedReference::Junctions(intent)),
                        AnalysisResult::SubstitutedAssemblyIntent(intent) => intent
                            .validate_binding(&revision)
                            .map(|()| IntendedReference::Substitutions(intent)),
                        AnalysisResult::ContributionEditIntent(intent) => intent
                            .validate_binding(&revision)
                            .map(|()| IntendedReference::Edits(intent)),
                        AnalysisResult::Error(error) => {
                            status.set(error);
                            return;
                        }
                        _ => {
                            status.set("Unexpected intended-construct response.".into());
                            return;
                        }
                    };
                    match result {
                        Ok(intent) => {
                            let original = match crate::assembly::edit::assembly_edit_model::SourceReview::capture(
                                &declared,
                                intent.fragments(),
                            ) {
                                Ok(value) => value,
                                Err(error) => {
                                    status.set(error);
                                    return;
                                }
                            };
                            source_receipt.set(Some((bound_epoch, Arc::new(original))));
                            let target_revision = intent.target().revision();
                            frozen.set(Some(Arc::new(Frozen {
                                identity: state,
                                epoch: bound_epoch,
                                intent,
                                sources: captures,
                                input_ids,
                                target_revision,
                            })));
                            status.set("Intended construct ready for review. No cloning method has been verified.".into());
                        }
                        Err(error) => status.set(error.to_string()),
                    }
                }
                Ok(_) => status.set("Intended-construct response identity changed.".into()),
                Err(error) => status.set(error),
            }
        });
    });
    let apply_junction = Callback::new(
        move |proposal: Arc<JunctionProposal>| -> Result<(), String> {
            if epoch.get_untracked() != proposal.epoch || draft.get_untracked() != proposal.before {
                return Err(
                    "Source slots or settings changed. Preview the junction edit again.".into(),
                );
            }
            let mut next = proposal.before.clone();
            if !next.substitutions.is_empty() || !next.edits.is_empty() {
                return Err("This primer proposal does not include internal substitutions, insertions or deletions. Keep the literal intended edits and enter primers explicitly.".into());
            }
            if proposal.report.len() != next.slots.len() {
                return Err(
                    "Proposed PCR preparation does not cover every pinned source slot.".into(),
                );
            }
            for i in 0..proposal.report.len() {
                let fragment = proposal
                    .report
                    .fragment(i)
                    .ok_or("Proposed fragment is missing.")?;
                let roles = next.methods.role_ids(next.slots.len())?;
                if fragment.input_id != roles[i].input_id
                    || fragment.role != roles[i].role
                    || fragment.reverse != next.methods.rows[i].reverse
                {
                    return Err("Proposed source order or orientation changed.".into());
                }
                match proposal.report.retained(i) {
                    Some(true) => continue,
                    Some(false) => {}
                    None => return Err("Proposed row preservation choices are incomplete.".into()),
                }
                let (forward, reverse, selection, fragment_index, template_strand, conditions) = match &fragment
                    .preparation
                {
                    collagene_core::cloning::assembly::recipe::Preparation::ExactPcr {
                        template_strand,
                        forward,
                        reverse,
                    } if matches!(&proposal.report, JunctionReport::Overlap(_)) => {
                        (forward, reverse, None, None, *template_strand, None)
                    }
                    collagene_core::cloning::assembly::recipe::Preparation::PcrCatalogDigest {
                        conditions,
                        template_strand,
                        forward,
                        reverse,
                        selection,
                        fragment_index,
                        repair: None,
                    } if matches!(
                        &proposal.report,
                        JunctionReport::GoldenGate(_) | JunctionReport::Fusion { .. }
                    ) =>
                    {
                        (
                            forward,
                            reverse,
                            Some(selection.clone()),
                            Some(*fragment_index),
                            *template_strand,
                            conditions.clone(),
                        )
                    }
                    _ => {
                        return Err(
                            "Selected row does not contain its proposed preparation chain.".into(),
                        );
                    }
                };
                let row = &mut next.methods.rows[i];
                row.preparation = if selection.is_some() {
                    OverlapPreparation::PcrCatalogDigest
                } else {
                    OverlapPreparation::ExactPcr
                };
                row.template_strand = template_strand;
                row.catalog_selection = selection;
                if let Some(conditions) = conditions {
                    row.reaction_conditions =
                        crate::restriction::reaction_conditions_model::Draft::from_saved(Some(
                            conditions,
                        ));
                }
                if let Some(index) = fragment_index {
                    row.fragment = (index + 1).to_string();
                }
                row.repair = None;
                row.treatment = crate::cloning::preparation::cloning_methods_model::Treatment::Keep;
                row.forward_name = forward.name.clone();
                row.forward = forward.annealing.as_str().into();
                row.forward_tail = forward.tail.clone();
                row.reverse_name = reverse.name.clone();
                row.reverse_primer = reverse.annealing.as_str().into();
                row.reverse_tail = reverse.tail.clone();
            }
            next.junction_edits = proposal.report.edits().to_vec();
            if let JunctionReport::Overlap(report) = &proposal.report {
                next.methods.overlaps = report
                    .recipe
                    .overlaps
                    .iter()
                    .map(usize::to_string)
                    .collect::<Vec<_>>()
                    .join(",");
            }
            next.sync_rows();
            next.validate()?;
            let identity = next.identity();
            let intent = IntendedReference::Junctions(proposal.report.intent().clone());
            let target_revision = intent.target().revision();
            undo_junction.set(Some((proposal.before.clone(), next.clone())));
            methods.set(next.methods.clone());
            let input_ids = next.methods.all_input_ids()?;
            draft.set(next);
            frozen.set(Some(Arc::new(Frozen {
                identity,
                epoch: proposal.epoch,
                intent,
                sources: proposal.sources.clone(),
                input_ids,
                target_revision,
            })));
            status.set("Junction edits and selected PCR primers applied together; other preparations were preserved. Compute fresh method evidence against the revised intended reference.".into());
            Ok(())
        },
    );
    let cancel = Callback::new(move |()| {
        generation.update(|g| *g = g.wrapping_add(1));
        if let Some(w) = worker.get_value() {
            w.close();
        }
        worker.set_value(None);
        busy.set(false);
        frozen.set(None);
        status.set("Intended-construct review cancelled.".into());
    });
    let replace_draft = Callback::new(move |next: Draft| {
        generation.update(|g| *g = g.wrapping_add(1));
        if let Some(w) = worker.get_value() {
            w.close();
        }
        worker.set_value(None);
        batch(|| {
            busy.set(false);
            frozen.set(None);
            undo_junction.set(None);
            strategy_chosen.set(next.slots.iter().any(|s| s.captured.is_some()));
            methods.set(next.methods.clone());
            draft.set(next);
            status.set(String::new());
        });
    });
    let load_saved = Callback::new(
        move |(before, next): (Draft, Draft)| -> Result<(), String> {
            if workspace_busy.get_untracked() || busy.get_untracked() {
                return Err(
                    "Wait for current workspace or intent work to finish before loading a design."
                        .into(),
                );
            }
            if draft.with_untracked(|current| current != &before) {
                return Err(
                    "The working draft changed before replacement. Choose the saved design again."
                        .into(),
                );
            }
            next.validate()?;
            // The source-bound loader already re-attested this declaration in the
            // current workspace. Replacement clears private intent/method receipts.
            batch(|| {
                previous_draft.set(Some(Arc::new(before)));
                replace_draft.run(next);
            });
            Ok(())
        },
    );
    let saved_busy = Memo::new(move |_| workspace_busy.get() || busy.get());
    // Bind browser commands to this mounted workspace, using the same worker and
    // cancel commands as the visible controls. No Frozen receipt is supplied by JS.
    use crate::ai::web_mcp_assembly as browser_assembly;
    use crate::ai::web_mcp_assembly::bridge::Command;
    use crate::ai::web_mcp_assembly::bridge::Endpoint;
    static NEXT_BROWSER_OWNER: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(1);
    let browser_owner = NEXT_BROWSER_OWNER.fetch_add(1, Ordering::Relaxed);
    let browser_alive = Arc::new(AtomicBool::new(true));
    let design_generation = RwSignal::new(0u32);
    let design_proposal = RwSignal::new(
        None::<
            Arc<(
                crate::ai::web_mcp_design::Proposal,
                crate::ai::web_mcp_design::Token,
            )>,
        >,
    );
    let status_command = Callback::new(
        move |expected: browser_assembly::Pin| -> Result<serde_json::Value, String> {
            draft.with_untracked(|draft| {
                browser_assembly::require_pin(draft, epoch.get_untracked(), &expected)
            })?;
            let same_review = review_pin.with_untracked(|held| held.as_ref() == Some(&expected));
            let target = if same_review {
                frozen.with_untracked(|held|held.as_ref().filter(|f|f.epoch==epoch.get_untracked()&&f.identity==identity.get_untracked()).map(|f|serde_json::json!({
            "name":f.intent.target().name,"length":f.intent.target().sequence.len(),"revision":f.target_revision,"intent_context_revision":f.intent.context()
        })))
            } else {
                None
            };
            let phase = if !same_review {
                "not_reviewed"
            } else if busy.get_untracked() {
                "running"
            } else if target.is_some() {
                "intended_reference_ready"
            } else {
                "stopped"
            };
            let proposal = design_proposal.with_untracked(|held| held.as_ref().map(|held| serde_json::json!({
                "token":held.1,"phase":if draft.with_untracked(|draft|held.0.is_current(draft,epoch.get_untracked())) {"awaiting_user_review"} else {"stale"},
                "changes":held.0.changes(),"application":"Only the visible user review controls can apply these declarations. No primers or physical method have been accepted."
            })));
            Ok(
                serde_json::json!({"phase":phase,"design_proposal":proposal,"review":if same_review {Some(browser_assembly::ReviewToken{owner:browser_owner,generation:generation.get_untracked()})}else{None},
            "message":if same_review {status.get_untracked()}else{"Review this exact draft to build its intended reference.".into()},
            "target":target,"scope":"intended_reference_review","method_verification":"not_performed_by_this_action","primer_apply":"Use the existing reviewed proposal controls; no automatic application is provided."}),
            )
        },
    );
    let command_alive = browser_alive.clone();
    let command = Callback::new(
        move |command: Command| -> Result<serde_json::Value, String> {
            if !command_alive.load(Ordering::Relaxed) {
                return Err("The assembly workspace was closed.".into());
            }
            if workspace_busy.get_untracked() {
                return Err("Wait for the workspace operation to finish.".into());
            }
            match command {
                Command::Status(expected) => status_command.run(expected),
                Command::Review(expected) => {
                    draft.with_untracked(|draft| {
                        browser_assembly::require_pin(draft, epoch.get_untracked(), &expected)
                    })?;
                    if busy.get_untracked() {
                        return Err("An intended-reference review is already running.".into());
                    }
                    review_intent.run(());
                    if !busy.get_untracked() {
                        return Err(status.get_untracked());
                    }
                    status_command.run(expected)
                }
                Command::Propose(request) => {
                    if busy.get_untracked() {
                        return Err("Wait for intended-reference review to finish before proposing changes.".into());
                    }
                    if design_proposal.with_untracked(|held| held.is_some()) {
                        return Err("A design proposal is already awaiting review. Discard that exact proposal before replacing it.".into());
                    }
                    let expected = request.expected.clone();
                    let review = source_review.get_untracked();
                    let proposal = draft.with_untracked(|draft| {
                        crate::ai::web_mcp_design::Proposal::prepare(
                            draft,
                            epoch.get_untracked(),
                            request,
                            review.as_deref(),
                        )
                    })?;
                    let generation = design_generation
                        .get_untracked()
                        .checked_add(1)
                        .ok_or("Design proposal generation exhausted; reopen the workspace.")?;
                    let token = crate::ai::web_mcp_design::Token {
                        owner: browser_owner,
                        generation,
                        revision: proposal.revision().into(),
                    };
                    design_generation.set(generation);
                    design_proposal.set(Some(Arc::new((proposal, token))));
                    status_command.run(expected)
                }
                Command::Discard(request) => {
                    draft.with_untracked(|draft| {
                        browser_assembly::require_pin(
                            draft,
                            epoch.get_untracked(),
                            &request.expected,
                        )
                    })?;
                    let held = design_proposal
                        .get_untracked()
                        .ok_or("No design proposal is awaiting review.")?;
                    if held.1 != request.proposal {
                        return Err(
                            "The pending proposal changed; read the current assembly context."
                                .into(),
                        );
                    }
                    design_proposal.set(None);
                    status_command.run(request.expected)
                }
                Command::Cancel(expected, token) => {
                    draft.with_untracked(|draft| {
                        browser_assembly::require_pin(draft, epoch.get_untracked(), &expected)
                    })?;
                    if token.owner != browser_owner
                        || token.generation != generation.get_untracked()
                        || !busy.get_untracked()
                        || review_pin.with_untracked(|held| held.as_ref() != Some(&expected))
                    {
                        return Err("This intended-reference review is no longer running; read its current context.".into());
                    }
                    cancel.run(());
                    status_command.run(expected)
                }
            }
        },
    );
    web_mcp.set(Some(Endpoint {
        alive: browser_alive.clone(),
        command,
    }));
    on_cleanup(move || {
        browser_alive.store(false, Ordering::Relaxed);
        let _ = web_mcp.try_update(|held| {
            if held
                .as_ref()
                .is_some_and(|endpoint| Arc::ptr_eq(&endpoint.alive, &browser_alive))
            {
                *held = None;
            }
        });
    });

    view! {<section class="assembly-workspace"><header class="assembly-workspace-header"><h2>"Design assembly"</h2><p class="helper">"Choose a fragment, then select its source in the sequence viewer above."</p><crate::assembly::draft::assembly_draft_actions::AssemblyDraftActions draft previous=previous_draft epoch locked=Memo::new(move |_| workspace_busy.get()||busy.get()) on_replace=replace_draft/></header>
        <crate::assembly::library::assembly_saved_designs_view::SavedDesigns draft epoch busy=saved_busy snapshot
            apply_attachments=apply_saved_attachments replace_draft=load_saved/>
        <Show when=move||!strategy_chosen.get()>
            <crate::assembly::draft::assembly_strategy_view::AssemblyStrategyPicker initial=draft.with_untracked(|d|d.methods.mode) on_start=Callback::new(move|mode|{
                draft.update(|d|{
                    d.methods.mode=mode;
                    d.sync_rows();
                    if mode==Mode::MultiSiteReference && d.methods.multisite.is_none(){d.methods.multisite=Some(crate::cloning::methods::multisite_form::Form::default());}
                });
                strategy_chosen.set(true);
            })/>
        </Show>
        <Show when=move||strategy_chosen.get()>
        <details class="assembly-design-settings"><summary>{move||draft.with(|d|format!("{} · {} · {}",d.name,method_overview(&d.methods),if d.circular{"circular"}else{"linear"}))}</summary>
        <div class="assembly-draft-settings"><label>"Assembly name"<input maxlength="256" aria-label="Assembly draft name" prop:value=move||draft.with(|d|d.name.clone()) on:input=move|e|draft.update(|d|d.name=event_target_value(&e))/></label><label class="checkbox-label"><input type="checkbox" aria-label="Intended assembly circular" prop:checked=move||draft.with(|d|d.circular) on:change=move|e|draft.update(|d|{d.circular=event_target_checked(&e);d.junction_edits.clear();d.sync_rows();})/>"Circular intended construct"</label></div>
        <label>"Cloning method"<select aria-label="Ordinary assembly cloning method" disabled=move||workspace_busy.get() prop:value=move||draft.with(|d|match d.methods.mode{Mode::Overlap=>"overlap",Mode::Ligation=>"ligation",Mode::GoldenGate=>"golden_gate",Mode::MultiSiteReference=>"multisite",Mode::Tailing=>"tailing"}) on:change=move|e|draft.update(|d|{
            d.methods.mode=match event_target_value(&e).as_str(){"ligation"=>Mode::Ligation,"golden_gate"=>Mode::GoldenGate,"multisite"=>Mode::MultiSiteReference,"tailing"=>Mode::Tailing,_=>Mode::Overlap};
            if d.methods.mode==Mode::MultiSiteReference && d.methods.multisite.is_none(){d.methods.multisite=Some(crate::cloning::methods::multisite_form::Form::default());}
            d.sync_rows();
        })><option value="overlap">"Gibson / HiFi / In-Fusion overlap"</option><option value="ligation">"Ordered ligation"</option><option value="golden_gate">"Type IIS / Golden Gate"</option><option value="multisite">"BP/LR multi-parent reference exchange"</option><option value="tailing">"TA / GC tailing and ligation"</option></select></label>
        <Show when=move||draft.with(|d|d.methods.mode==Mode::MultiSiteReference)><crate::assembly::draft::assembly_multisite_view::MultisiteSettings draft epoch busy=workspace_busy/></Show>
        <p class="helper" aria-label="Current assembly method"><strong>"Current method: "</strong>{move||draft.with(|d|method_overview(&d.methods))}". Review the intended construct to edit preparation settings in step 3."</p>
        </details>
        <nav class="assembly-slot-ribbon" aria-label="Assembly fragment order">
        <For each=move||draft.with(|d|(0..d.slots.len()).collect::<Vec<_>>()) key=|i|*i children=move|i|view!{
            <button class="assembly-slot-tab" class:active=move||draft.with(|d|d.active_slot==i) aria-pressed=move||draft.with(|d|d.active_slot==i) on:click=move|_|draft.update(|d|d.active_slot=i)>
                <strong>{move||draft.with(|d|format!("{}. {}",i+1,d.slots.get(i).map(|s|s.label.as_str()).unwrap_or("Fragment")))}</strong>
                <span>{move||draft.with(|d|d.slots.get(i).map(|s|match s.captured.as_deref(){Some(c)=>format!("{} · {}",c.source.name,crate::assembly::draft::model::selector_label(&s.selector)),None=>"Choose source…".into()}).unwrap_or_default())}</span>
            </button>
        }/>
        </nav>
        <details class="assembly-source-edit" prop:open=move||frozen.with(|value|value.is_none())><summary>"Source region and orientation"</summary><div class="assembly-source-stage"><div class="assembly-slot-list" role="group" aria-label="Ordered assembly fragment slots">
        <For each=move||draft.with(|d|(0..d.slots.len()).collect::<Vec<_>>()) key=|i|*i children=move|i|view!{<Show when=move||draft.with(|d|d.active_slot==i)><article class="assembly-slot"><h3>{move||draft.with(|d|d.slots.get(i).map(|s|s.label.clone()).unwrap_or_default())}</h3><details class="assembly-fragment-options"><summary>"Fragment options"</summary><label>"Slot label"<input maxlength="80" aria-label=format!("Assembly slot {} label",i+1) prop:value=move||draft.with(|d|d.slots.get(i).map(|s|s.label.clone()).unwrap_or_default()) on:input=move|e|draft.update(|d|if let Some(s)=d.slots.get_mut(i){s.label=event_target_value(&e)})/></label><div class="method-row-actions"><button class="text-tool" disabled=move||i==0||draft.with(|d|d.alternatives.is_some()) on:click=move|_|draft.update(|d|{if d.alternatives.is_none()&&i>0&&i<d.slots.len(){d.slots.swap(i,i-1);d.methods.swap_rows(i,i-1);d.junction_edits.clear();d.active_slot=i-1;d.sync_rows();}})>"Move earlier"</button><button class="text-tool" disabled=move||draft.with(|d|d.alternatives.is_some()||i+1>=d.slots.len()) on:click=move|_|draft.update(|d|{if d.alternatives.is_none()&&i+1<d.slots.len(){d.slots.swap(i,i+1);d.methods.swap_rows(i,i+1);d.junction_edits.clear();d.active_slot=i+1;d.sync_rows();}})>"Move later"</button><button class="text-tool" disabled=move||draft.with(|d|d.alternatives.is_some()||d.slots.len()<=if d.methods.mode==Mode::Tailing{2}else{1}) on:click=move|_|draft.update(|d|{if d.alternatives.is_none()&&d.slots.len()>if d.methods.mode==Mode::Tailing{2}else{1}&&i<d.slots.len(){if let Err(error)=d.clear_slot_substitutions(i){status.set(error);return;}d.slots.remove(i);d.methods.remove_row(i);d.junction_edits.clear();d.active_slot=d.active_slot.min(d.slots.len()-1);d.sync_rows();}})>"Remove slot"</button></div></details>
            <p>{move||draft.with(|d|d.slots.get(i).map(|s|match s.captured.as_deref(){Some(c)=>format!("{} · {}",c.source.name,crate::assembly::draft::model::selector_label(&s.selector)),None=>"No source fragment set".into()}).unwrap_or_default())}</p>
            <p class="helper"><strong>"Declared preparation: "</strong>{move||draft.with(|d|preparation_overview(&d.methods,i))}</p>
            <crate::history::run_origin_view::SourceOriginPicker source=Memo::new(move|_|draft.with(|d|d.slots.get(i).and_then(|s|s.captured.as_ref().map(|c|c.source.clone())))) epoch get_workspace=snapshot selected=Memo::new(move|_|draft.with(|d|d.slots.get(i).and_then(|s|s.material_origin.clone()))) on_select=Callback::new(move|pin|{let next=draft.get_untracked().set_material_origin(i,pin)?;draft.set(next);Ok(())})/>
            <details class="assembly-source-details"><summary>"Source details"</summary><p class="helper">{move||draft.with(|d|d.slots.get(i).and_then(|s|s.captured.as_ref()).map(|c|format!("Pinned source event {} · revision {}",c.source.binding.event_id,c.source.binding.revision)).unwrap_or_default())}</p></details>
            <p class="helper">{move||{let captured=draft.with(|d|d.slots.get(i).and_then(|s|s.captured.clone()));captured.map(|c|sources.with(|current|if current.iter().any(|s|s.binding==c.source.binding){"Current source version is pinned"}else{"Using the captured source version; current tabs do not replace it"})).unwrap_or("")}}</p>
            <label class="checkbox-label"><input type="checkbox" aria-label=format!("Reverse assembly slot {}",i+1) disabled=move||draft.with(|d|d.methods.rows.get(i).is_some_and(|row|row.preparation==OverlapPreparation::AnnealedOligos)) prop:checked=move||draft.with(|d|d.slots.get(i).is_some_and(|slot|match &slot.selector{assembly_intent::Selector::ReferenceRange{reverse,..}=>*reverse,assembly_intent::Selector::MaterialStrandRange{strand,..}=>*strand==collagene_core::chemistry::Strand::Bottom})) on:change=move|e|{match draft.get_untracked().orient_fragment(i,event_target_checked(&e)){Ok(next)=>draft.set(next),Err(error)=>status.set(error)}}/>"Reverse source contribution"</label>
            <Show when=move||draft.with(|d|d.methods.rows.get(i).is_some_and(|row|row.preparation==OverlapPreparation::AnnealedOligos))><p class="helper">"Swap the ordered oligo roles to reverse a paired insert. This uses the full previous bottom as the new top contribution and resets its neighboring junction declarations. Inspect and select the geometry again."</p><button class="button secondary" on:click=move|_|{
                let next=draft.with_untracked(|d|d.swap_oligo_roles(i));
                match next{Ok(next)=>{draft.set(next);status.set("Oligo roles and chemistry swapped; neighboring junction declarations reset. Review the full new top contribution and reconstruct annealing and ligation again.".into());},Err(error)=>status.set(error)}
            }>"Swap ordered oligo roles"</button></Show>
            <Show when=move||draft.with(|d|d.slots.get(i).is_some_and(|slot|matches!(&slot.selector,assembly_intent::Selector::MaterialStrandRange{..})))><label class="checkbox-label"><input type="checkbox" aria-label=format!("Reverse prepared fragment {}",i+1) prop:checked=move||draft.with(|d|d.methods.rows.get(i).is_some_and(|row|row.reverse)) on:change=move|event|draft.update(|d|{if let Some(row)=d.methods.rows.get_mut(i){row.reverse=event_target_checked(&event);}})/>"Reverse prepared fragment"</label><p class="helper">"Preparation orientation is separate from the selected source strand. Primer proposals declare their own source frame."</p></Show>
            <div class="method-row-actions"><button class="text-tool" on:click=move|_|{let current=draft.get_untracked();let result=(||->Result<_,String>{let slot=current.slots.get(i).ok_or("Choose a source first")?;let capture=slot.captured.as_ref().ok_or("Pin a source first")?;let review=source_review.get_untracked();let length=slot.selector.length().map_err(|e|e.to_string())?;let regions=crate::assembly::edit::assembly_edit_model::top_regions(&current,i,0,length,review.as_deref())?;Ok(SelectionNavigation{source:capture.source.binding.clone(),regions,description:crate::assembly::draft::model::selector_label(&slot.selector)})})();match result{Ok(nav)=>on_navigate.run(nav),Err(error)=>status.set(error)}}>"View source region"</button></div>
        </article></Show>}/></div>
        <button class="button secondary" disabled=move||draft.with(|d|d.alternatives.is_some()||d.slots.len()>=if d.methods.mode==Mode::Tailing{2}else{12}) on:click=move|_|{
            let mut next=draft.get_untracked();
            if next.alternatives.is_some()||next.slots.len()>=if next.methods.mode==Mode::Tailing{2}else{12}{return;}
            match crate::assembly::edit::assembly_substitution_model::append_slot(&next) {Ok(next)=>draft.set(next),Err(error)=>status.set(error)}
        }>"Add fragment slot"</button>
        <details class="assembly-library-options"><summary>"Source alternatives / library"</summary><crate::assembly::draft::assembly_alternatives_view::AssemblyAlternativesView draft epoch/></details>
        <Show when=move||draft.with(|d|d.methods.mode==Mode::Ligation&&d.alternatives.is_none())><details class="assembly-fragment-options"><summary>"Use an annealed oligo insert"</summary><p class="helper">"Pin the top oligo for this slot using Use whole source, then open the bottom oligo in the viewer and capture it below. Both ordered sequences are 5′→3′. The top source region declares your intended reference contribution; add any intended extra junction bases separately, before reaction simulation."</p><button class="button secondary" on:click=move|_|pin_bottom.run(())>"Capture viewer source as bottom oligo"</button><p class="helper">{move||draft.with(|d|d.slots.get(d.active_slot).and_then(|slot|slot.anneal_bottom.as_ref()).map(|capture|format!("Pinned bottom: {} · {} bases",capture.source.name,capture.molecule.sequence.len())).unwrap_or_else(||"No bottom oligo captured for this slot.".into()))}</p></details></Show>
        <section class="assembly-set-fragment"><p>{move||format!("Active viewer source: {} · set {}",active_source.get().name,draft.with(|d|d.slots.get(d.active_slot).map(|s|s.label.as_str()).unwrap_or("selected slot").to_owned()))}</p><button class="button primary" on:click=move|_|pin.run(1)>"Set from sequence selection"</button><button class="button secondary" on:click=move|_|pin.run(0)>"Use whole source"</button><button class="button secondary" on:click=move|_|pin.run(2)>"Keep outside selection"</button><p class="helper">"Keep outside selection retains the complementary path on circular DNA, including an origin crossing. This sets intended source coordinates; it does not digest, amplify or prove a physical fragment."</p></section></div></details>
        {move||design_proposal.get().map(|held| {
            let apply=held.clone();
            view!{<section class="assembly-intended-review"><h3>"AI proposed design changes"</h3>
                <p class="helper">"Review these original contribution coordinates and literals. Accepting changes updates only this draft's intended declarations, then runs the normal intended-reference review. Primers and cloning-method verification remain separate."</p>
                {held.0.changes().iter().map(|change|view!{<pre>{change.description()}</pre>}).collect_view()}
                <button class="button secondary" disabled=move||busy.get()||workspace_busy.get() on:click=move|_|{
                    let result=draft.with_untracked(|draft|apply.0.apply(draft,epoch.get_untracked()));
                    match result {Ok(next)=>{draft.set(next);design_proposal.set(None);review_intent.run(());},Err(error)=>status.set(error)}
                }>"Apply declarations and review intended construct"</button>
                <button class="text-tool" on:click=move|_|design_proposal.set(None)>"Discard proposal"</button>
            </section>}
        })}
        <button class="button primary" disabled=move||busy.get() on:click=move|_|review_intent.run(())>"Review intended construct"</button><Show when=move||busy.get()><button class="text-tool" on:click=move|_|cancel.run(())>"Cancel intended construct"</button></Show><p class="assembly-workspace-status" role="status">{move||status.get()}</p>
        <Show when=move||draft.with(|d|matches!(d.methods.mode,Mode::Tailing|Mode::Ligation)||!d.substitutions.is_empty()||!d.edits.is_empty())><crate::assembly::edit::literal_intent_view::LiteralIntentEditor draft/></Show>
        <ActualStrandSource draft source_review/>
        <crate::assembly::edit::assembly_substitution_view::SubstitutionEditor draft active_source selected_region on_navigate source_review/>
        <crate::assembly::edit::assembly_edit_view::ContributionEditEditor draft active_source selected_region on_navigate source_review/>
        <Show when=move||draft.with(|d|!matches!(d.methods.mode,Mode::Tailing|Mode::Ligation)&&d.substitutions.is_empty()&&d.edits.is_empty())><details class="assembly-junction-options"><summary>"Edit intended junctions"</summary><JunctionEditor draft epoch on_apply=apply_junction/></details></Show>
        <Show when=move||undo_junction.with(|undo|undo.is_some())><button class="text-tool" disabled=move||undo_junction.with(|undo|undo.as_ref().is_none_or(|(_,applied)|!draft.get().same_declaration(applied))) on:click=move|_|{
            if let Some((mut previous,applied))=undo_junction.get_untracked(){
                let current=draft.get_untracked();
                if !current.same_declaration(&applied){status.set("Undo is available only before further draft changes. Your later changes were preserved.".into());return;}
                if current.active_slot<previous.slots.len(){previous.active_slot=current.active_slot;}
                methods.set(previous.methods.clone());draft.set(previous);frozen.set(None);undo_junction.set(None);status.set("Previous assembly draft restored. Review its intent and recompute evidence; no positive result was restored.".into());
            }
        }>"Undo applied junction edit"</button></Show>
        <Show when=move||undo_junction.with(|undo|undo.as_ref().is_some_and(|(_,applied)|!draft.get().same_declaration(applied)))><p class="helper">"Undo available before further draft changes. Later source, primer or design edits are preserved."</p></Show>
        {move||frozen.get().map(|bound|view!{<FrozenWorkspace bound draft methods epoch snapshot workspace_busy on_navigate on_review on_finalize on_products primer_sources list_primers resolve_primer/>})}
    </Show></section>}
}

fn intended_window(molecule: &Molecule, boundary: usize) -> String {
    let dna = molecule.sequence.as_str().as_bytes();
    let n = dna.len();
    if n == 0 {
        return "Intended sequence unavailable".into();
    }
    let (left, right) = if molecule.topology == Topology::Circular {
        let b = boundary % n;
        let count = 20.min(n);
        (
            (0..count)
                .map(|i| dna[(b + n - count + i) % n] as char)
                .collect::<String>(),
            (0..count)
                .map(|i| dna[(b + i) % n] as char)
                .collect::<String>(),
        )
    } else {
        let b = boundary.min(n);
        (
            String::from_utf8_lossy(&dna[b.saturating_sub(20)..b]).into_owned(),
            String::from_utf8_lossy(&dna[b..b.saturating_add(20).min(n)]).into_owned(),
        )
    };
    format!("{left}|{right}")
}
#[component]
fn FrozenWorkspace(
    bound: Arc<Frozen>,
    draft: RwSignal<Draft>,
    methods: RwSignal<MethodsState>,
    epoch: RwSignal<u64>,
    snapshot: Callback<(), Result<collagene_core::documents::workspace::Workspace, String>>,
    workspace_busy: RwSignal<bool>,
    on_navigate: Callback<SelectionNavigation>,
    on_review: Callback<ReviewRequest, Result<(), String>>,
    on_finalize: Callback<ReviewRequest, Result<(), String>>,
    on_products: Callback<Vec<Molecule>, Result<(), String>>,
    primer_sources: Memo<Vec<SourceInfo>>,
    list_primers: Callback<EventBinding, Result<primer_reuse::List, String>>,
    resolve_primer: Callback<ResolveRequest, Result<primer_reuse::Resolved, String>>,
) -> impl IntoView {
    let bound = StoredValue::new(bound);
    let junction = RwSignal::new(0usize);
    let active_preparation = Memo::new(move |_| draft.with(|d| d.active_slot));
    let edit_preparation = Callback::new(move |row: usize| {
        draft.update(|d| {
            if row < d.slots.len() {
                d.active_slot = row;
            }
        });
    });
    let status = RwSignal::new(String::new());
    let target = Memo::new(move |_| bound.with_value(|f| f.intent.target().clone()));
    let target_info = Memo::new(move |_| {
        bound.with_value(|f| TargetContext::Intent {
            context_revision: f.intent.context().to_string(),
            revision: f.target_revision.clone(),
        })
    });
    let query_options = Memo::new(move |_| {
        bound.with_value(|f| {
            f.sources
                .iter()
                .enumerate()
                .map(|(i, c)| {
                    (
                        i,
                        format!(
                            "{} · {}",
                            f.identity.slot_label(i).unwrap_or("Bottom oligo"),
                            c.source.name
                        ),
                    )
                })
                .collect::<Vec<_>>()
        })
    });
    let get_query = Callback::new(move |i: usize| {
        bound.with_value(|f| {
            f.sources.get(i).and_then(|c| {
                collagene_model::bounds::model_preflight(&[&c.molecule])
                    .ok()
                    .map(|_| c.molecule.clone())
            })
        })
    });
    let source_info = Callback::new(move |i: usize| {
        bound.with_value(|f| f.sources.get(i).map(|c| c.source.clone()))
    });
    let prepare_review = Callback::new(move |candidate: Candidate| {
        let outcome = bound.with_value(|f| {
            if epoch.get_untracked() != f.epoch
                || draft.with_untracked(Draft::identity) != f.identity
            {
                return Err("Source slots changed; review the intended construct again.".into());
            }
            if candidate.plan.target_revision != f.target_revision
                || candidate.plan.inputs.len() != f.sources.len()
                || candidate.plan.inputs.iter().enumerate().any(|(i, p)| {
                    p.id != f.input_ids[i] || p.revision != f.sources[i].source.binding.revision
                })
            {
                return Err(
                    "Candidate inputs differ from the pinned source order; no candidate opened."
                        .into(),
                );
            }
            Ok(ReviewRequest {
                candidate,
                sources: f.sources.iter().map(|c| c.source.binding.clone()).collect(),
                intent: f.intent.clone(),
            })
        });
        outcome
    });
    let review = Callback::new(move |candidate: Candidate| {
        let result = prepare_review
            .run(candidate)
            .and_then(|request| on_review.run(request));
        match result { Ok(()) => status.set("Candidate opened for fresh verification against the intended reference. No experimental success is established.".into()), Err(error) => status.set(error) }
    });
    let finalize = Callback::new(move |candidate: Candidate| {
        prepare_review
            .run(candidate)
            .and_then(|request| on_finalize.run(request))
    });
    let navigate = Callback::new(move |i: usize| {
        if let Some(nav) = bound.with_value(|f| {
            f.intent.fragments().get(i).and_then(|region| {
                f.sources.get(i).map(|source| SelectionNavigation {
                    source: source.source.binding.clone(),
                    regions: region
                        .top_mappings
                        .iter()
                        .map(|m| (m.source_start, m.source_end, m.reverse))
                        .collect(),
                    description: crate::assembly::draft::model::selector_label(&region.selector),
                })
            })
        }) {
            draft.update(|d| d.active_slot = i);
            on_navigate.run(nav);
        }
    });
    let diagnostic_source = Callback::new(move |(document, row): (usize, usize)| {
        let nav = bound.with_value(|f| {
            if epoch.get_untracked() != f.epoch
                || draft.with_untracked(Draft::identity) != f.identity
                || document != row
            {
                return Err(
                    "The pinned source slots changed; review the intended construct again."
                        .to_string(),
                );
            }
            let region = f
                .intent
                .fragments()
                .get(row)
                .ok_or("The preparation row is unavailable.")?;
            let source = f
                .sources
                .get(row)
                .ok_or("The original source is unavailable.")?;
            Ok(SelectionNavigation {
                source: source.source.binding.clone(),
                regions: region
                    .top_mappings
                    .iter()
                    .map(|m| (m.source_start, m.source_end, m.reverse))
                    .collect(),
                description: crate::assembly::draft::model::selector_label(&region.selector),
            })
        })?;
        draft.update(|d| d.active_slot = row);
        on_navigate.run(nav);
        Ok(())
    });
    let reference_guard = Callback::new(move |()| -> Result<(), String> {
        bound.with_value(|f| {
            if workspace_busy.get_untracked()
                || epoch.get_untracked() != f.epoch
                || draft.with_untracked(Draft::identity) != f.identity
            {
                return Err(
                    "The intended reference or pinned source changed; review intent again.".into(),
                );
            }
            Ok(())
        })
    });
    let use_assembly_primers = Callback::new(
        move |(before, next): (Draft, Draft)| -> Result<(), String> {
            reference_guard.run(())?;
            let count = draft.with_untracked(|current| {
                crate::assembly::primer_design::assembly_edited_primer_commit::validate(
                    &before, &next, current,
                )
            })?;
            if before.identity() != next.identity() {
                return Err(
                    "The proposed PCR oligos changed the independently reviewed construct.".into(),
                );
            }
            batch(|| {
                methods.set(next.methods.clone());
                draft.set(next);
                status.set(format!("Suggested PCR preparations applied to {count} contributions. Reconstruct the complete assembly against the unchanged intended reference."));
            });
            Ok(())
        },
    );
    let use_extension_primers = Callback::new(
        move |(before, next): (Draft, Draft)| -> Result<(), String> {
            reference_guard.run(())?;
            let index = draft.with_untracked(|current| {
                crate::assembly::primer_design::assembly_extension_primer_commit::validate(
                    &before, &next, current,
                )
            })?;
            // The proposal component re-attests its exact historical source at
            // use time. This host guard independently preserves the frozen
            // declaration and permits only the native transition above.
            if before.identity() != next.identity() {
                return Err(
                    "Primer proposal changed the independently reviewed source declaration.".into(),
                );
            }
            batch(|| {
                methods.set(next.methods.clone());
                draft.set(next);
                status.set(format!("Overlap-extension primers applied to fragment {}. Reconstruct the complete assembly against the unchanged intended reference.", index + 1));
            });
            Ok(())
        },
    );
    let reference_inputs = Callback::new(move |()| -> Result<Vec<PlanInput>, String> {
        reference_guard.run(())?;
        bound.with_value(|f| {
            collagene_model::bounds::model_preflight(
                &f.sources.iter().map(|s| &s.molecule).collect::<Vec<_>>(),
            )
            .map_err(|e| e.to_string())?;
            Ok(f.input_ids
                .iter()
                .zip(&f.sources)
                .map(|(r, c)| PlanInput {
                    id: r.clone(),
                    revision: c.source.binding.revision.clone(),
                    molecule: c.molecule.clone(),
                })
                .collect())
        })
    });
    let reference_current = Memo::new(move |_| {
        bound.with_value(|f| {
            !workspace_busy.get()
                && epoch.get() == f.epoch
                && draft.with(Draft::identity) == f.identity
        })
    });
    view! {<section class="assembly-intended-review"><details class="assembly-intent-details"><summary>{bound.with_value(|f|format!("Intended construct: {} · {} bp",f.intent.target().name,f.intent.target().sequence.len()))}</summary><h3>"2. Review the intended construct"</h3><p>{bound.with_value(|f|format!("{} · {} bp · {}",f.intent.target().name,f.intent.target().sequence.len(),if f.intent.target().topology==Topology::Circular{"circular"}else{"linear"}))}</p><p class="helper">"This sequence is your declared design requirement, assembled from chosen source regions and explicit junction or internal sequence edits. It is separate from any simulated cloning outcome. Editing a slot clears this review."</p>
        {bound.with_value(|f|match &f.intent{IntendedReference::Junctions(intent)=>Some(view!{<p class="helper">{format!("{} explicit junction edits are included. Ribbon output spans include each adjusted source core and its outgoing synthetic insert; source navigation selects only the adjusted original-source core.",intent.edits.len())}</p>}),_=>None})}
        <div class="assembly-fragment-ribbon" role="group" aria-label="Intended construct fragments">{bound.with_value(|f|f.intent.fragments().iter().enumerate().map(|(i,fragment)|view!{<button class="assembly-ribbon-fragment" on:click=move|_|navigate.run(i)>{format!("{}{} · {} bp · {}–{}",f.identity.slot_label(i).expect("Frozen intent source slot"),match &fragment.selector{assembly_intent::Selector::ReferenceRange{reverse:true,..}|assembly_intent::Selector::MaterialStrandRange{strand:collagene_core::chemistry::Strand::Bottom,..}=>" ←",_=>" →"},fragment.output_end-fragment.output_start,fragment.output_start+1,fragment.output_end)}</button>}).collect_view())}</div>
        {bound.with_value(|f|{let count=f.intent.fragments().len().saturating_sub(1)+usize::from(f.intent.target().topology==Topology::Circular);if count==0{return view!{<p>"One linear fragment has no internal junction."</p>}.into_any();}view!{<div class="assembly-intended-junction">{move||bound.with_value(|f|{let i=junction.get().min(count-1);let next=(i+1)%f.intent.fragments().len();let boundary=f.intent.fragments()[i].output_end;view!{<h4>{format!("Intended junction {} · {} → {}{}",i+1,f.identity.slot_label(i).expect("Frozen intent source slot"),f.identity.slot_label(next).expect("Frozen intent downstream slot"),if i+1==f.intent.fragments().len(){" · circle closure"}else{""})}</h4><pre>{intended_window(f.intent.target(),boundary)}</pre><p class="helper">{format!("| marks intended product boundary {} (zero-based). These letters are not a digest/ligation/homology validation.",if boundary==f.intent.target().sequence.len()&&f.intent.target().topology==Topology::Circular{0}else{boundary})}</p><button class="text-tool" on:click=move|_|navigate.run(i)>"View upstream source"</button><button class="text-tool" on:click=move|_|navigate.run(next)>"View downstream source"</button>}})}<div class="method-row-actions"><button class="text-tool" disabled=move||junction.get()==0 on:click=move|_|junction.update(|i|*i=i.saturating_sub(1))>"Previous intended junction"</button><button class="text-tool" disabled=move||{junction.get()+1>=count} on:click=move|_|junction.update(|i|*i=(*i+1).min(count-1))>"Next intended junction"</button></div></div>}.into_any()})}
        </details>
        {bound.with_value(|f|f.intent.warnings().iter().map(|warning|view!{<p class="helper">{warning.clone()}</p>}).collect_view())}
        <h3>{move||draft.with(|d|if d.methods.mode==Mode::MultiSiteReference{"3. Reconstruct and choose a reference component"}else{"3. Prepare fragments and verify the method"})}</h3><p class="helper">{move||draft.with(|d|if d.methods.mode==Mode::MultiSiteReference{"Reconstruct the whole original parents under the reviewed site profile, then choose a reciprocal component explicitly. The intended reference remains a separate requirement; sequence agreement does not establish physical recombination."}else{"PCR rows begin with no primer sequences. Suggest and explicitly apply primers, or choose another declared preparation. The pinned original sources remain the inputs; the intended reference is never used as physical source material."})}</p>
        <Show when=move||draft.with(|d|d.methods.mode==Mode::MultiSiteReference)><crate::assembly::draft::assembly_multisite_view::MultisiteReview draft epoch busy=workspace_busy snapshot capture=reference_inputs guard=reference_guard target_revision=bound.with_value(|f|f.target_revision.clone()) on_review=review on_finalize=finalize/></Show>
        {bound.with_value(|f| {
            let edited=matches!(&f.intent,IntendedReference::Edits(_)|IntendedReference::Substitutions(_));
            let material=f.identity.has_material_sources();
            if !edited&&!material{return None;}
            let selected_intent=StoredValue::new(Arc::new(f.intent.clone()));
            Some(view!{
                <Show when=move||edited><crate::assembly::primer_design::assembly_inverse_pcr_view::InversePcrPreparationView draft intent=selected_intent.get_value() epoch snapshot current=reference_current guard=reference_guard on_review=review on_finalize=finalize/></Show>
                <crate::assembly::primer_design::assembly_edited_primer_view::EditedAssemblyPrimerView draft intent=selected_intent.get_value() epoch snapshot current=reference_current guard=reference_guard on_use=use_assembly_primers/>
                <Show when=move||edited><crate::assembly::primer_design::assembly_extension_primer_view::ExtensionPrimerProposalView draft intent=selected_intent.get_value() epoch snapshot current=reference_current guard=reference_guard on_use=use_extension_primers/></Show>
            })
        })}
        <Show when=move||draft.with(|d|d.methods.mode!=Mode::MultiSiteReference)>        <CloningMethodsPanel molecule=target query_options get_query source_info target_info selected_region=Callback::new(|()|None) config=methods epoch on_candidate=review on_finalize=Some(finalize) on_products primer_sources list_primers resolve_primer active_row=Some(active_preparation) on_edit_row=Some(edit_preparation) fixed_sources=true get_workspace=Some(snapshot) disable_plain_primer_design=bound.with_value(|f|f.identity.has_sequence_edits()||f.identity.has_material_sources()) on_source=Some(diagnostic_source)/></Show>
        <p role="status">{move||status.get()}</p>
    </section>}
}

fn default_regeneration(draft: &Draft, proposed: &[assembly_junction_edit::Edit]) -> Vec<bool> {
    let n = draft.slots.len();
    if n == 0 || n > 12 {
        return Vec::new();
    }
    let mut selected = (0..n)
        .map(|i| {
            draft.methods.rows.get(i).is_some_and(|row| {
                matches!(
                    row.preparation,
                    OverlapPreparation::ExactPcr | OverlapPreparation::PcrCatalogDigest
                ) && (row.forward.trim().is_empty() || row.reverse_primer.trim().is_empty())
            })
        })
        .collect::<Vec<_>>();
    let junctions = draft
        .junction_edits
        .iter()
        .chain(proposed)
        .map(|edit| edit.junction)
        .collect::<std::collections::BTreeSet<_>>();
    for junction in junctions {
        if junction >= n {
            continue;
        }
        let old = draft
            .junction_edits
            .iter()
            .find(|edit| edit.junction == junction);
        let new = proposed.iter().find(|edit| edit.junction == junction);
        if old != new {
            selected[junction] = true;
            if old.is_some_and(|edit| edit.remove_right > 0)
                || new.is_some_and(|edit| edit.remove_right > 0)
            {
                selected[(junction + 1) % n] = true;
            }
        }
    }
    selected
}
#[derive(Clone, PartialEq, Eq)]
struct JunctionSettings {
    overlap: String,
    minimum: String,
    maximum: String,
    tm: String,
}
impl Default for JunctionSettings {
    fn default() -> Self {
        Self {
            overlap: "20".into(),
            minimum: "18".into(),
            maximum: "30".into(),
            tm: "60".into(),
        }
    }
}
enum JunctionReport {
    Overlap(selective_junction_edit::Report),
    GoldenGate(golden_gate_junction_edit::Report),
    Fusion {
        report: collagene_core::cloning::golden_gate::fusion_design::Materialized,
        retained_rows: Vec<usize>,
    },
}
impl JunctionReport {
    fn intent(&self) -> &assembly_junction_edit::Intent {
        match self {
            Self::Overlap(r) => &r.intent,
            Self::GoldenGate(r) => &r.intent,
            Self::Fusion { report, .. } => &report.original_intent,
        }
    }
    fn edits(&self) -> &[assembly_junction_edit::Edit] {
        match self {
            Self::Overlap(r) => &r.request.edits,
            Self::GoldenGate(r) => &r.request.edits,
            Self::Fusion { report, .. } => &report.request.base.edits,
        }
    }
    fn len(&self) -> usize {
        match self {
            Self::Overlap(r) => r.recipe.fragments.len(),
            Self::GoldenGate(r) => r.recipe.fragments.len(),
            Self::Fusion { report, .. } => report.recipe.fragments.len(),
        }
    }
    fn fragment(
        &self,
        index: usize,
    ) -> Option<&collagene_core::cloning::assembly::recipe::FragmentSpec> {
        match self {
            Self::Overlap(r) => r.recipe.fragments.get(index),
            Self::GoldenGate(r) => r.recipe.fragments.get(index).map(|f| &f.source),
            Self::Fusion { report, .. } => report
                .recipe
                .fragments
                .get(index)
                .map(|fragment| &fragment.source),
        }
    }
    fn retained(&self, index: usize) -> Option<bool> {
        match self {
            Self::Fusion { report, .. } => report
                .request
                .base
                .rows
                .get(index)
                .map(|row| matches!(row, golden_gate_junction_edit::Row::Retain(_))),
            Self::Overlap(r) => r
                .request
                .rows
                .get(index)
                .map(|row| matches!(row, selective_junction_edit::Row::Retain(_))),
            Self::GoldenGate(r) => r
                .request
                .rows
                .get(index)
                .map(|row| matches!(row, golden_gate_junction_edit::Row::Retain(_))),
        }
    }
    fn designed_count(&self) -> usize {
        match self {
            Self::Overlap(r) => r.designs.iter().filter(|d| d.is_some()).count(),
            Self::GoldenGate(r) => r.designs.iter().filter(|d| d.is_some()).count(),
            Self::Fusion { report, .. } => {
                report.designs.iter().filter(|row| row.is_some()).count()
            }
        }
    }
    fn retained_rows(&self) -> &[usize] {
        match self {
            Self::Overlap(r) => &r.retained_rows,
            Self::GoldenGate(r) => &r.retained_rows,
            Self::Fusion { retained_rows, .. } => retained_rows,
        }
    }
    fn warnings(&self) -> &[String] {
        match self {
            Self::Overlap(r) => &r.warnings,
            Self::GoldenGate(r) => &r.warnings,
            Self::Fusion { report, .. } => &report.warnings,
        }
    }
    fn validate_binding(&self, revision: &str) -> Result<(), collagene_model::Error> {
        match self {
            Self::Overlap(r) => r.validate_binding(revision),
            Self::GoldenGate(r) => r.validate_binding(revision),
            Self::Fusion { report, .. } => report.validate_binding(revision),
        }
    }
}
struct JunctionProposal {
    before: Draft,
    report: JunctionReport,
    sources: Vec<Arc<CapturedSource>>,
    settings: JunctionSettings,
    regenerate: Vec<bool>,
    epoch: u64,
}
#[component]
fn JunctionEditor(
    draft: RwSignal<Draft>,
    epoch: RwSignal<u64>,
    on_apply: Callback<Arc<JunctionProposal>, Result<(), String>>,
) -> impl IntoView {
    let proposed = RwSignal::new(draft.with_untracked(|d| d.junction_edits.clone()));
    let junction = RwSignal::new(0usize);
    let remove_left = RwSignal::new("0".to_string());
    let remove_right = RwSignal::new("0".to_string());
    let inserted = RwSignal::new(String::new());
    let settings = RwSignal::new(JunctionSettings {
        overlap: draft.with_untracked(|d| d.methods.overlaps.clone()),
        ..Default::default()
    });
    let regenerate = RwSignal::new(default_regeneration(
        &draft.get_untracked(),
        &proposed.get_untracked(),
    ));
    let manual_regeneration = RwSignal::new(false);
    let result = RwSignal::new(None::<Arc<JunctionProposal>>);
    let status = RwSignal::new(String::new());
    let pending = RwSignal::new(false);
    let generation = RwSignal::new(0u32);
    let worker = StoredValue::new_local(None::<Rc<CachedWorker>>);
    let alive = Arc::new(AtomicBool::new(true));
    let cleanup = alive.clone();
    on_cleanup(move || {
        cleanup.store(false, Ordering::Relaxed);
        if let Some(w) = worker.get_value() {
            w.close();
        }
    });
    let source_context = Memo::new(move |_| draft.with(Draft::identity));
    Effect::new(move |_| {
        let _ = (source_context.get(), epoch.get());
        let edits = draft.with_untracked(|d| d.junction_edits.clone());
        let first = edits.iter().find(|edit| edit.junction == 0);
        remove_left.set(
            first
                .map(|e| e.remove_left.to_string())
                .unwrap_or_else(|| "0".into()),
        );
        remove_right.set(
            first
                .map(|e| e.remove_right.to_string())
                .unwrap_or_else(|| "0".into()),
        );
        inserted.set(first.map(|e| e.inserted.clone()).unwrap_or_default());
        regenerate.set(default_regeneration(&draft.get_untracked(), &edits));
        manual_regeneration.set(false);
        settings.update(|settings| {
            settings.overlap = draft.with_untracked(|d| d.methods.overlaps.clone())
        });
        proposed.set(edits);
        junction.set(0);
    });
    Effect::new(move |_| {
        let before = draft.get();
        let edits = proposed.get();
        if !manual_regeneration.get() {
            let defaults = default_regeneration(&before, &edits);
            if regenerate.get_untracked() != defaults {
                regenerate.set(defaults);
            }
        }
    });
    Effect::new(move |_| {
        let _ = (
            draft.get(),
            proposed.get(),
            settings.get(),
            regenerate.get(),
            epoch.get(),
            remove_left.get(),
            remove_right.get(),
            inserted.get(),
            junction.get(),
        );
        generation.update(|g| *g = g.wrapping_add(1));
        if let Some(w) = worker.get_value() {
            w.close();
        }
        worker.set_value(None);
        pending.set(false);
        result.set(None);
    });
    let load_edit = Callback::new(move |index: usize| {
        junction.set(index);
        let edit =
            proposed.with_untracked(|edits| edits.iter().find(|e| e.junction == index).cloned());
        remove_left.set(
            edit.as_ref()
                .map(|e| e.remove_left.to_string())
                .unwrap_or_else(|| "0".into()),
        );
        remove_right.set(
            edit.as_ref()
                .map(|e| e.remove_right.to_string())
                .unwrap_or_else(|| "0".into()),
        );
        inserted.set(edit.map(|e| e.inserted).unwrap_or_default());
    });
    let unstaged = Memo::new(move |_| {
        let current = junction.get();
        proposed.with(|edits| {
            let edit = edits.iter().find(|edit| edit.junction == current);
            remove_left.get().parse::<usize>().ok()
                != Some(edit.map(|e| e.remove_left).unwrap_or(0))
                || remove_right.get().parse::<usize>().ok()
                    != Some(edit.map(|e| e.remove_right).unwrap_or(0))
                || inserted.get() != edit.map(|e| e.inserted.as_str()).unwrap_or("")
        })
    });
    let set_edit = move |_| {
        let outcome = (|| -> Result<_, String> {
            let junction = junction.get_untracked();
            let count =
                draft.with_untracked(|d| d.slots.len().saturating_sub(1) + usize::from(d.circular));
            if junction >= count {
                return Err("Choose an available junction.".into());
            }
            let left = remove_left.get_untracked();
            let right = remove_right.get_untracked();
            let inserted = inserted.get_untracked();
            if left.len() > 20 || right.len() > 20 || inserted.len() > 200 {
                return Err(
                    "Trim fields or inserted text exceed their bounds (200 raw inserted bytes)."
                        .into(),
                );
            }
            Ok(assembly_junction_edit::Edit {
                junction,
                remove_left: left.parse().map_err(|_| "Enter an upstream trim length.")?,
                remove_right: right
                    .parse()
                    .map_err(|_| "Enter a downstream trim length.")?,
                inserted,
            })
        })();
        match outcome {
            Ok(edit) => {
                proposed.update(|items| {
                    if let Some(old) = items.iter_mut().find(|e| e.junction == edit.junction) {
                        *old = edit;
                    } else {
                        items.push(edit);
                        items.sort_by_key(|e| e.junction);
                    }
                });
                status.set("Proposed edit updated. Preview all junction edits and PCR primers before applying.".into());
            }
            Err(error) => status.set(error),
        }
    };
    let design_alive = alive.clone();
    let preview = Callback::new(move |()| {
        let before = draft.get_untracked();
        let edits = proposed.get_untracked();
        let form = settings.get_untracked();
        let selected_rows = regenerate.get_untracked();
        let captured_epoch = epoch.get_untracked();
        let bound_edits = edits.clone();
        let prepared = (|| -> Result<_, String> {
            if !matches!(before.methods.mode, Mode::Overlap | Mode::GoldenGate) {
                return Err(
                    "Choose overlap or Golden Gate mode before preparing junction edits.".into(),
                );
            }
            if edits.is_empty() && before.junction_edits.is_empty() {
                return Err("Set at least one proposed junction edit first.".into());
            }
            let mut candidate = before.clone();
            candidate.junction_edits = edits.clone();
            crate::assembly::draft::model::require_reference_selectors(&candidate)?;
            let (inputs, base, sources) = candidate.inputs_request(None)?;
            if selected_rows.len() != before.slots.len() {
                return Err("Select a preparation choice for every slot.".into());
            }
            let mut retained = before.methods.clone();
            for (i, selected) in selected_rows.iter().enumerate() {
                if *selected {
                    retained.rows[i] = OverlapInput {
                        source_generation: before.methods.rows[i].source_generation.clone(),
                        document: i,
                        reverse: before.methods.rows[i].reverse,
                        ..Default::default()
                    };
                }
            }
            let (revision, job) = if before.methods.mode == Mode::GoldenGate {
                let request = crate::assembly::edit::assembly_fusion_design::base_request(
                    &before,
                    base,
                    edits,
                    &selected_rows,
                )?;
                let revision = golden_gate_junction_edit::design_revision(&inputs, &request)
                    .map_err(|e| e.to_string())?;
                (
                    revision,
                    AnalysisJob::DesignGoldenGateJunctionPrimers { inputs, request },
                )
            } else {
                if form.overlap.len() > 256
                    || [&form.minimum, &form.maximum, &form.tm]
                        .iter()
                        .any(|text| text.len() > 32)
                {
                    return Err("Primer settings exceed their field bounds.".into());
                }
                let (_, fragments) =
                    crate::cloning::preparation::cloning_preparations::prepared_rows(
                        &retained,
                        inputs.len(),
                    )?;
                let rows = fragments
                    .into_iter()
                    .enumerate()
                    .map(|(i, fragment)| {
                        if selected_rows[i] {
                            selective_junction_edit::Row::DesignPcr {
                                role: fragment.role,
                                input_id: fragment.input_id,
                                reverse: fragment.reverse,
                            }
                        } else {
                            selective_junction_edit::Row::Retain(fragment)
                        }
                    })
                    .collect();
                let request = selective_junction_edit::Request {
                    base,
                    edits,
                    rows,
                    overlaps: form
                        .overlap
                        .split(',')
                        .map(|part| {
                            part.trim().parse::<usize>().map_err(
                                |_| "Enter one comma-separated homology length per junction.",
                            )
                        })
                        .collect::<Result<Vec<_>, _>>()?,
                    primer_options: collagene_core::cloning::pcr::primers::DesignOptions {
                        min_length: form
                            .minimum
                            .parse()
                            .map_err(|_| "Enter a minimum annealing length.")?,
                        max_length: form
                            .maximum
                            .parse()
                            .map_err(|_| "Enter a maximum annealing length.")?,
                        target_tm_c: form
                            .tm
                            .parse()
                            .map_err(|_| "Enter an annealing target Tm.")?,
                        max_pairs: 1,
                        ..Default::default()
                    },
                    profile_options: before.methods.profile_options()?,
                };
                let revision = selective_junction_edit::design_revision(&inputs, &request)
                    .map_err(|e| e.to_string())?;
                (
                    revision,
                    AnalysisJob::DesignSelectiveAssemblyPrimers { inputs, request },
                )
            };
            let w = Rc::new(CachedWorker::new()?);
            Ok((job, sources, revision, w))
        })();
        let (job, sources, revision, w) = match prepared {
            Ok(value) => value,
            Err(error) => {
                status.set(error);
                return;
            }
        };
        if let Some(old) = worker.get_value() {
            old.close();
        }
        worker.set_value(Some(w.clone()));
        let token = generation.get_untracked().wrapping_add(1);
        generation.set(token);
        pending.set(true);
        result.set(None);
        status.set(
            "Building the edited intended reference and proposing original-source PCR primers…"
                .into(),
        );
        let alive = design_alive.clone();
        spawn_local(async move {
            let response = w
                .ask(Request {
                    id: token,
                    revision: revision.clone(),
                    job,
                })
                .await;
            w.close();
            if !alive.load(Ordering::Relaxed)
                || generation.get_untracked() != token
                || epoch.get_untracked() != captured_epoch
                || draft.get_untracked() != before
                || settings.get_untracked() != form
                || proposed.get_untracked() != bound_edits
                || regenerate.get_untracked() != selected_rows
            {
                return;
            }
            pending.set(false);
            worker.set_value(None);
            let outcome = (|| -> Result<_, String> {
                let response = response?;
                if response.id != token || response.revision != revision {
                    return Err("Junction proposal response identity changed.".into());
                }
                let report = match response.result {
                    AnalysisResult::SelectiveAssemblyPrimers(report)
                        if before.methods.mode == Mode::Overlap =>
                    {
                        JunctionReport::Overlap(report)
                    }
                    AnalysisResult::GoldenGateJunctionPrimers(report)
                        if before.methods.mode == Mode::GoldenGate =>
                    {
                        JunctionReport::GoldenGate(report)
                    }
                    AnalysisResult::Error(error) => return Err(error),
                    _ => return Err("Unexpected junction primer proposal result.".into()),
                };
                report
                    .validate_binding(&revision)
                    .map_err(|e| e.to_string())?;
                if report.edits() != bound_edits.as_slice() {
                    return Err("Junction proposal edits differ from the requested changes.".into());
                }
                Ok(Arc::new(JunctionProposal {
                    before,
                    report,
                    sources,
                    settings: form,
                    regenerate: selected_rows,
                    epoch: captured_epoch,
                }))
            })();
            match outcome {
                Ok(proposal) => {
                    result.set(Some(proposal));
                    status.set("Proposal ready. Only checked rows will receive new original-source PCR primers; retained preparations and oligos stay unchanged.".into());
                }
                Err(error) => status.set(error),
            }
        });
    });
    let apply = move |_| {
        let outcome = (|| -> Result<(), String> {
            let proof = result
                .get_untracked()
                .ok_or("Preview the proposed junction edits first.")?;
            if proof.before != draft.get_untracked()
                || proof.epoch != epoch.get_untracked()
                || proof.settings != settings.get_untracked()
                || proof.report.edits() != proposed.get_untracked().as_slice()
                || proof.regenerate != regenerate.get_untracked()
            {
                return Err(
                    "Sources, proposed edits or primer settings changed. Preview again.".into(),
                );
            }
            on_apply.run(proof)
        })();
        if let Err(error) = outcome {
            status.set(error);
        }
    };
    let fusion_capture = Callback::new(
        move |()| -> Result<crate::assembly::edit::assembly_fusion_view::Capture, String> {
            if unstaged.get_untracked() {
                return Err("Set the literal junction edit before fusion-position search.".into());
            }
            let before = draft.get_untracked();
            let edits = proposed.get_untracked();
            let selected_rows = regenerate.get_untracked();
            let mut declared = before.clone();
            declared.junction_edits = edits.clone();
            crate::assembly::draft::model::require_reference_selectors(&declared)?;
            let (inputs, base, sources) = declared.inputs_request(None)?;
            let base = crate::assembly::edit::assembly_fusion_design::base_request(
                &before,
                base,
                edits,
                &selected_rows,
            )?;
            Ok(crate::assembly::edit::assembly_fusion_view::Capture {
                before,
                epoch: epoch.get_untracked(),
                inputs,
                base,
                sources,
                regenerate: selected_rows,
            })
        },
    );
    let fusion_guard = Callback::new(
        move |saved: Arc<crate::assembly::edit::assembly_fusion_view::Capture>| -> Result<(), String> {
            if saved.epoch != epoch.get_untracked()
                || saved.before != draft.get_untracked()
                || saved.base.edits != proposed.get_untracked()
                || saved.regenerate != regenerate.get_untracked()
                || unstaged.get_untracked()
            {
                return Err("The original sources, literal edits or regeneration choices changed; search fusion placements again.".into());
            }
            Ok(())
        },
    );
    let fusion_apply = Callback::new(
        move |saved: Arc<crate::assembly::edit::assembly_fusion_view::MaterializedEvidence>| -> Result<(), String> {
            fusion_guard.run(saved.capture.clone())?;
            let retained_rows = saved
                .capture
                .regenerate
                .iter()
                .enumerate()
                .filter_map(|(index, regenerate)| (!regenerate).then_some(index))
                .collect();
            let proposal = Arc::new(JunctionProposal {
                before: saved.capture.before.clone(),
                report: JunctionReport::Fusion {
                    report: saved.report.clone(),
                    retained_rows,
                },
                sources: saved.capture.sources.clone(),
                settings: settings.get_untracked(),
                regenerate: saved.capture.regenerate.clone(),
                epoch: saved.capture.epoch,
            });
            on_apply.run(proposal)
        },
    );
    let fusion_labels = Memo::new(move |_| {
        draft.with(|draft| {
            let count = draft.slots.len().saturating_sub(1) + usize::from(draft.circular);
            (0..count)
                .map(|index| {
                    format!(
                        "{} → {}{}",
                        draft.slots[index].label,
                        draft.slots[(index + 1) % draft.slots.len()].label,
                        if index + 1 == draft.slots.len() {
                            " (circle closure)"
                        } else {
                            ""
                        }
                    )
                })
                .collect::<Vec<_>>()
        })
    });
    view! {<details class="assembly-junction-editor"><summary>"Edit a junction and redesign its PCR preparation"</summary>
        <p class="helper">"Trim bases from the upstream core's assembly-facing end and/or the downstream core's beginning, then insert literal DNA between them. Original source regions stay pinned. Choose which rows receive newly designed PCR primers. Apply replaces only those checked rows; retained manual primers and digest preparations stay unchanged. This does not establish cloning success."</p>
        <label>"Junction"<select aria-label="Assembly junction to edit" prop:value=move||junction.get().to_string() on:change=move|e|if let Ok(index)=event_target_value(&e).parse(){load_edit.run(index)}>{move||draft.with(|d|{let count=d.slots.len().saturating_sub(1)+usize::from(d.circular);(0..count).map(|i|view!{<option value=i.to_string()>{format!("{} → {}{}",d.slots[i].label,d.slots[(i+1)%d.slots.len()].label,if i+1==d.slots.len(){" (circle closure)"}else{""})}</option>}).collect_view()})}</select></label>
        <div class="primer-grid"><label>"Remove upstream bases"<input type="number" min="0" aria-label="Junction upstream trim" prop:value=move||remove_left.get() on:input=move|e|remove_left.set(event_target_value(&e))/></label><label>"Remove downstream bases"<input type="number" min="0" aria-label="Junction downstream trim" prop:value=move||remove_right.get() on:input=move|e|remove_right.set(event_target_value(&e))/></label></div>
        <label>"Insert literal DNA at this junction"<textarea maxlength="200" aria-label="Junction inserted DNA" prop:value=move||inserted.get() on:input=move|e|inserted.set(event_target_value(&e))></textarea></label>
        <button class="button secondary" on:click=set_edit>"Set proposed junction edit"</button>
        {move||proposed.get().into_iter().map(|edit|{let index=edit.junction;view!{<p>{format!("Junction {}: remove {} upstream / {} downstream; insert {}",index+1,edit.remove_left,edit.remove_right,if edit.inserted.is_empty(){"(none)"}else{&edit.inserted})}<button class="text-tool" on:click=move|_|{proposed.update(|items|items.retain(|e|e.junction!=index));if junction.get_untracked()==index{remove_left.set("0".into());remove_right.set("0".into());inserted.set(String::new());}}>"Remove proposed edit"</button></p>}}).collect_view()}
        <h4>"Preparation choices"</h4><p class="helper">"Recommended boxes start with changed upstream fragments, downstream trims and unconfigured PCR rows. After you change a box, later edit changes do not override that choice. Retaining an incompatible preparation can prevent reconstruction; it is never silently replaced."</p>
        {move||draft.with(|d|d.slots.iter().enumerate().map(|(i,slot)|view!{<label class="checkbox-label"><input type="checkbox" aria-label=format!("Regenerate PCR primers for slot {}",i+1) prop:checked=move||regenerate.with(|rows|rows.get(i).copied().unwrap_or(false)) on:change=move|e|{manual_regeneration.set(true);regenerate.update(|rows|if let Some(row)=rows.get_mut(i){*row=event_target_checked(&e)});}/>{format!("{} — regenerate PCR primers",slot.label)}</label><p class="helper">{move||draft.with(|current|current.methods.rows.get(i).map(|row|format!("{}: {}",if regenerate.with(|rows|rows.get(i).copied().unwrap_or(false)){"Replace after explicit Apply"}else{"Retain exact current preparation"},row.summary())).unwrap_or_default())}</p>}).collect_view())}
        <Show when=move||draft.with(|d|d.methods.mode==Mode::Overlap)>
        <p class="helper">{move||draft.with(|d|format!("Enter {} explicit homology lengths in junction order, including circle closure when present. Existing values are preserved; no uniform length is inferred.",d.slots.len().saturating_sub(1)+usize::from(d.circular)))}</p>
        <div class="primer-grid"><label>"Per-junction homology lengths (comma-separated)"<input aria-label="Edited junction homology" prop:value=move||settings.with(|s|s.overlap.clone()) on:input=move|e|settings.update(|s|s.overlap=event_target_value(&e))/></label><label>"Minimum annealing length"<input type="number" min="1" max="200" aria-label="Edited junction primer minimum" prop:value=move||settings.with(|s|s.minimum.clone()) on:input=move|e|settings.update(|s|s.minimum=event_target_value(&e))/></label><label>"Maximum annealing length"<input type="number" min="1" max="200" aria-label="Edited junction primer maximum" prop:value=move||settings.with(|s|s.maximum.clone()) on:input=move|e|settings.update(|s|s.maximum=event_target_value(&e))/></label><label>"Target annealing Tm (°C)"<input type="number" aria-label="Edited junction primer target Tm" prop:value=move||settings.with(|s|s.tm.clone()) on:input=move|e|settings.update(|s|s.tm=event_target_value(&e))/></label></div>
        <p class="helper">"The selected Gibson/HiFi/In-Fusion profile is retained. Primer metrics are model calculations, not reagent conditions or demonstrated binding. Concrete source and inserted DNA are required for these primer proposals; ambiguous calls are not silently resolved."</p>
        </Show>
        <Show when=move||draft.with(|d|d.methods.mode==Mode::GoldenGate)>
            <h4>"Golden Gate junction preparation"</h4><p class="helper">"Checked rows receive PCR-added Type IIS sites and the literal junction insert. Internal outgoing fusion bases come from the next adjusted source core; no overlap Tm or homology-length setting is used. Retained rows keep their exact digest/primer/end-treatment choices."</p>
            <div class="primer-grid"><label>"Type IIS enzyme"<select aria-label="Junction Type IIS enzyme" prop:value=move||format!("{:?}",draft.with(|d|d.methods.golden_gate.enzyme)) on:change=move|e|draft.update(|d|d.methods.golden_gate.enzyme=match event_target_value(&e).as_str(){"BsmBI"=>collagene_core::cloning::golden_gate::GoldenGateEnzyme::BsmBI,"BbsI"=>collagene_core::cloning::golden_gate::GoldenGateEnzyme::BbsI,"SapI"=>collagene_core::cloning::golden_gate::GoldenGateEnzyme::SapI,_=>collagene_core::cloning::golden_gate::GoldenGateEnzyme::BsaI})><option value="BsaI">"BsaI · 4-base fusion"</option><option value="BsmBI">"BsmBI · 4-base fusion"</option><option value="BbsI">"BbsI · 4-base fusion"</option><option value="SapI">"SapI · 3-base fusion"</option></select></label>
            <label>"Pre-recognition clamp"<input maxlength="12" aria-label="Junction Golden Gate clamp" prop:value=move||draft.with(|d|d.methods.golden_gate.clamp.clone()) on:input=move|e|draft.update(|d|d.methods.golden_gate.clamp=event_target_value(&e))/></label><label>"Cleavage spacer"<input maxlength="2" aria-label="Junction Golden Gate spacer" prop:value=move||draft.with(|d|d.methods.golden_gate.spacer.clone()) on:input=move|e|draft.update(|d|d.methods.golden_gate.spacer=event_target_value(&e))/></label>
            <label>"Forward annealing length"<input type="number" min="1" max="200" aria-label="Junction Golden Gate forward length" prop:value=move||draft.with(|d|d.methods.golden_gate.forward_length.clone()) on:input=move|e|draft.update(|d|d.methods.golden_gate.forward_length=event_target_value(&e))/></label><label>"Reverse annealing length"<input type="number" min="1" max="200" aria-label="Junction Golden Gate reverse length" prop:value=move||draft.with(|d|d.methods.golden_gate.reverse_length.clone()) on:input=move|e|draft.update(|d|d.methods.golden_gate.reverse_length=event_target_value(&e))/></label></div>
            <Show when=move||draft.with(|d|!d.circular)&&regenerate.with(|selected|selected.last()==Some(&true))><label>"Terminal right fusion (assembly direction, 5′→3′)"<input maxlength="4" aria-label="Junction Golden Gate terminal fusion" prop:value=move||draft.with(|d|d.methods.rows.last().map(|row|row.right_fusion.clone()).unwrap_or_default()) on:input=move|e|draft.update(|d|if let Some(row)=d.methods.rows.last_mut(){row.right_fusion=event_target_value(&e)})/></label><p class="helper">"A selected terminal linear fragment requires this explicit fusion DNA (3 bases for SapI; 4 for the other listed enzymes). It is not filled from the intended target."</p></Show>
            <label class="checkbox-label"><input type="checkbox" aria-label="Junction Golden Gate allow retained nicks" prop:checked=move||draft.with(|d|d.methods.allow_nicked) on:change=move|e|draft.update(|d|d.methods.allow_nicked=event_target_checked(&e))/>"Allow explicitly retained nicks"</label><p class="helper">"Clamp/spacer and terminal fusion are explicit design choices, not supplier conditions. Preview checks the actual PCR/digest ends and the full ordered ligation; it does not establish reaction efficiency or mixture uniqueness."</p>
        </Show>
        <Show when=move||unstaged.get()><p class="helper">"Set the proposed junction edit before previewing, so the changed trim and inserted bases are included."</p></Show>
        <button class="button primary" disabled=move||pending.get()||unstaged.get()||draft.with(|d|!matches!(d.methods.mode,Mode::Overlap|Mode::GoldenGate)) on:click=move|_|preview.run(())>"Preview edited construct and PCR primers"</button>
        <button class="text-tool" disabled=move||!pending.get() on:click=move|_|{generation.update(|g|*g=g.wrapping_add(1));if let Some(w)=worker.get_value(){w.close();}worker.set_value(None);pending.set(false);result.set(None);status.set("Junction proposal cancelled.".into());}>"Cancel junction proposal"</button>
        <p role="status">{move||status.get()}</p>
        <Show when=move||draft.with(|draft|draft.methods.mode==Mode::GoldenGate)>
            <crate::assembly::edit::assembly_fusion_view::AssemblyFusionView epoch invalidation=generation junction_labels=fusion_labels capture=fusion_capture guard=fusion_guard on_apply=fusion_apply/>
        </Show>
        {move||result.get().map(|proof|view!{<section class="assembly-junction-proposal"><h3>"Edited construct and proposed PCR preparation"</h3><p>{format!("{} · {} bases",proof.report.intent().target.name,proof.report.intent().target.sequence.len())}</p><pre>{let dna=proof.report.intent().target.sequence.as_str();format!("{}{}",dna.chars().take(240).collect::<String>(),if dna.len()>240{"… (bounded preview)"}else{""})}</pre>
            {proof.report.intent().fragments.iter().enumerate().map(|(i,fragment)|view!{<p>{format!("{}: {} → {}; output span {}–{} includes any outgoing synthetic insert",proof.before.slots[i].label,crate::assembly::draft::model::selector_label(&proof.before.slots[i].selector),crate::assembly::draft::model::selector_label(&fragment.selector),fragment.output_start+1,fragment.output_end)}</p>}).collect_view()}
            <p>{format!("Retained preparations: {} · regenerated PCR rows: {}",proof.report.retained_rows().len(),proof.report.designed_count())}</p>
            {proof.report.retained_rows().iter().filter_map(|&i|proof.before.methods.rows.get(i).map(|row|(i,row))).map(|(i,row)|view!{<p class="helper">{format!("Retained {}: {}",proof.before.slots[i].label,row.summary())}</p>}).collect_view()}
            {junction_primers(&proof.report)}
            {proof.report.warnings().iter().chain(proof.report.intent().warnings.iter()).map(|warning|view!{<p class="helper">{warning.clone()}</p>}).collect_view()}
            <button class="button primary" on:click=apply>"Apply junction edits and selected PCR primers"</button>
        </section>})}
    </details>}
}

fn junction_primers(report: &JunctionReport) -> AnyView {
    match report {
        JunctionReport::Fusion{..}=>().into_any(),
        JunctionReport::Overlap(report)=>report.designs.iter().enumerate().filter_map(|(i,design)|design.as_ref().map(|design|(i,design))).map(|(i,design)|view!{<article><h4>{format!("Fragment {} proposed overlap PCR primers",i+1)}</h4>{design.pairs.first().map(|pair|view!{<pre>{format!("{}\nOrdered 5′→3′: {}{}\nAnnealing: {}\nAdded tail: {}\n\n{}\nOrdered 5′→3′: {}{}\nAnnealing: {}\nAdded tail: {}",pair.forward.primer.name,pair.forward.primer.tail,pair.forward.primer.annealing.as_str(),pair.forward.primer.annealing.as_str(),pair.forward.primer.tail,pair.reverse.primer.name,pair.reverse.primer.tail,pair.reverse.primer.annealing.as_str(),pair.reverse.primer.annealing.as_str(),pair.reverse.primer.tail)}</pre><p>{format!("Annealing Tm {} / {} · pair ΔTm {:.1}°C",pair.forward.analysis.tm.as_ref().map(|t|format!("{:.1}°C",t.tm_c)).unwrap_or_else(||"unavailable".into()),pair.reverse.analysis.tm.as_ref().map(|t|format!("{:.1}°C",t.tm_c)).unwrap_or_else(||"unavailable".into()),pair.tm_difference_c)}</p>})}{design.warnings.iter().map(|warning|view!{<p class="helper">{warning.clone()}</p>}).collect_view()}</article>}).collect_view().into_any(),
        JunctionReport::GoldenGate(report)=>report.designs.iter().enumerate().filter_map(|(i,design)|design.as_ref().map(|design|(i,design))).map(|(i,design)|view!{<article><h4>{format!("Fragment {} proposed Type IIS PCR/digest chain",i+1)}</h4><pre>{format!("{}\nOrdered 5′→3′: {}{}\nAnnealing: {}\nTail: clamp [{}] · motif [{}] · spacer [{}] · fusion [{}] · inserted [{}]\n\n{}\nOrdered 5′→3′: {}{}\nAnnealing: {}\nTail: clamp [{}] · motif [{}] · spacer [{}] · fusion [{}] · inserted [{}]",design.forward.name,design.forward.tail,design.forward.annealing.as_str(),design.forward.annealing.as_str(),design.forward_tail.clamp,design.forward_tail.motif,design.forward_tail.spacer,design.forward_tail.fusion,design.forward_tail.inserted,design.reverse.name,design.reverse.tail,design.reverse.annealing.as_str(),design.reverse.annealing.as_str(),design.reverse_tail.clamp,design.reverse_tail.motif,design.reverse_tail.spacer,design.reverse_tail.fusion,design.reverse_tail.inserted)}</pre><p>{format!("Actual digested fragment {} · assembly-oriented top DNA: {}{}",design.fragment_index+1,design.oriented_top.chars().take(240).collect::<String>(),if design.oriented_top.len()>240{"… (bounded preview)"}else{""})}</p><p>{format!("Left end: {}",crate::cloning::preparation::diagnostics::method_end(Some(&design.left)))}</p><p>{format!("Right end: {}",crate::cloning::preparation::diagnostics::method_end(Some(&design.right)))}</p><p class="helper">"Tail segments are shown in each actual oligo's own 5′→3′ direction. Literal insert and fusion segments may therefore be reverse-complemented; the displayed digested top DNA is in assembly direction."</p></article>}).collect_view().into_any(),
    }
}

#[component]
fn ActualStrandSource(
    draft: RwSignal<Draft>,
    source_review: Memo<Option<Arc<crate::assembly::edit::assembly_edit_model::SourceReview>>>,
) -> impl IntoView {
    let page = RwSignal::new(0usize);
    let context = Memo::new(move |_| {
        draft.with(|d| {
            d.slots.get(d.active_slot).map(|slot| {
                (
                    d.active_slot,
                    slot.captured.as_ref().map(|c| c.source.binding.clone()),
                    slot.selector.clone(),
                )
            })
        })
    });
    Effect::new(move || {
        let _ = context.get();
        page.set(0);
    });
    let material = move || {
        draft.with(|d| {
            d.slots.get(d.active_slot).is_some_and(|slot| {
                matches!(
                    &slot.selector,
                    assembly_intent::Selector::MaterialStrandRange { .. }
                )
            })
        })
    };
    let text = move || {
        draft.with(|d| {
            source_review.with(|review| {
                review
                    .as_ref()
                    .ok_or("Review intended DNA to inspect this actual strand".to_string())
                    .and_then(|review| review.get(d, d.active_slot))
                    .map(|fragment| {
                        let length = fragment.original_sequence.len();
                        let start = (page.get().saturating_mul(100)).min(length);
                        let end = (start + 100).min(length);
                        format!(
                            "{}–{} · {}",
                            start + 1,
                            end,
                            &fragment.original_sequence.as_str()[start..end]
                        )
                    })
                    .unwrap_or_else(|error| error)
            })
        })
    };
    let pages = move || {
        draft.with(|d| {
            source_review.with(|review| {
                review
                    .as_ref()
                    .and_then(|review| review.get(d, d.active_slot).ok())
                    .map(|f| f.original_sequence.len().div_ceil(100).max(1))
                    .unwrap_or(1)
            })
        })
    };
    view! {<Show when=material><details><summary>"Actual source strand"</summary><p class="helper">"Sequence is shown in the selected strand’s own 5′→3′ direction. Extra strand bases keep their actual coordinates; the original top-reference viewer is unchanged."</p><p class="helper">"Use an explicit preparation for this actual-strand contribution. Primer suggestions currently use original reference ranges."</p><pre>{text}</pre><button class="text-tool" disabled=move||page.get()==0 on:click=move|_|page.update(|p|*p=p.saturating_sub(1))>"Previous bases"</button><button class="text-tool" disabled={move||page.get()+1>=pages()} on:click=move|_|page.update(|p|*p+=1)>"Next bases"</button>{move||draft.with(|d|source_review.with(|review|review.as_ref().and_then(|review|review.get(d,d.active_slot).ok()).map(|f|view!{<p>{format!("Original contribution {} bases · inherited source limits {}",f.original_sequence.len(),if f.incomplete{"present"}else{"none reported"})}</p>{f.source_domain.as_ref().map(|domain|view!{<p>{format!("Source strand domain [{}, {}) · {} stored bases",domain.start,domain.end,domain.length)}</p>})}<details><summary>"Original top-reference mappings"</summary>{f.top_mappings.iter().map(|m|view!{<p>{format!("Reference {}–{} → contribution {}–{}{}",m.source_start+1,m.source_end,m.destination_start+1,m.destination_start+m.source_end-m.source_start,if m.reverse{" · reverse"}else{""})}</p>}).collect_view()}<p class="helper">"Bases outside these mappings have no original top-reference position."</p></details>})))}</details></Show>}
}
