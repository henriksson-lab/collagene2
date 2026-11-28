use collagene_core::documents::transfer;
use collagene_core::material::lineage::EventBinding;
use collagene_core::planning::plan::DesignPlan;
use collagene_core::planning::plan::Operation;
use collagene_core::planning::plan::PlanInput;
use collagene_core::planning::plan::Source;
use collagene_core::protocol;
use collagene_core::verification;
use collagene_core::verification::evidence::VerificationReport;
use collagene_core::{
    chemistry::{self, tracked_duplex as modified_duplex},
    cloning::{
        ligation::outcomes as ligation_outcomes, pcr::mutagenesis, recombination,
        restriction::catalog as catalog_digest,
    },
    documents::workspace,
    material::lineage as run_lineage,
};
use collagene_model::Alphabet;
use collagene_model::Molecule;
use collagene_model::Topology;
use leptos::prelude::*;
use wasm_bindgen::JsCast;

use crate::workspace::session::OpenDocument;
use crate::workspace::session::RestoredWorkspace;
use crate::workspace::session::WorkspaceSession;

use crate::workspace::session::current_event_binding;

use crate::viewer::sequence::selection::ViewerSelection;
use crate::workspace::shell_state::ShellState;

#[derive(Clone)]
pub(crate) struct PendingFinalize {
    pub(crate) workspace: std::sync::Arc<collagene_core::documents::workspace::Workspace>,
    pub(crate) plan: DesignPlan,
    pub(crate) target: collagene_core::material::lineage::EventBinding,
    pub(crate) sources: Vec<collagene_core::material::lineage::InputEventBinding>,
    pub(crate) source_fingerprint: String,
}
#[derive(Clone, Copy)]
pub(crate) struct AssemblyController {
    pub(crate) session: WorkspaceSession,
    pub(crate) viewer: ViewerSelection,
    pub(crate) shell: ShellState,
}
impl AssemblyController {
    pub(crate) fn verify_candidate(
        self,
        candidate: RwSignal<Option<transfer::Candidate>>,
        candidate_events: RwSignal<
            Option<(
                String,
                u64,
                collagene_core::material::lineage::EventBinding,
                Vec<collagene_core::material::lineage::InputEventBinding>,
            )>,
        >,
        last_plan: RwSignal<Option<(DesignPlan, Molecule)>>,
        last_plan_target: RwSignal<
            Option<(String, u64, collagene_core::material::lineage::EventBinding)>,
        >,
        plan_report: RwSignal<Option<VerificationReport>>,
        run_plan: impl Fn(
            DesignPlan,
            Molecule,
            Option<Vec<collagene_core::material::lineage::InputEventBinding>>,
            Option<collagene_core::material::lineage::EventBinding>,
        ) + Copy
        + 'static,
    ) -> impl Fn(web_sys::MouseEvent) + Copy + 'static {
        let session = self.session;
        let WorkspaceSession {
            documents,
            active,
            workspace_epoch,
            ..
        } = self.session;
        let ViewerSelection { molecule, .. } = self.viewer;
        let ShellState { message, .. } = self.shell;
        move |_: web_sys::MouseEvent| {
            let Some(candidate) = candidate.get_untracked() else {
                return;
            };
            let preferences = candidate_events
                .get_untracked()
                .filter(|(revision, _, _, _)| *revision == candidate.plan_revision);
            let is_reference = candidate.plan.steps.iter().any(|step| {
                matches!(
                    step,
                    Operation::RecombineReference { .. }
                        | Operation::MultiSiteRecombineReference { .. }
                        | Operation::Topo { .. }
                        | Operation::CheckedTopo { .. }
                        | Operation::TopoActions { .. }
                )
            });
            // Exact native source preferences may name retained events even when a
            // different version is currently open. The native verifier reconstructs
            // their declarations; reading history grants no physical authority.
            let exact_preferences = preferences.is_some();
            if let Some((_, epoch, target, sources)) = &preferences {
                let valid = *epoch == workspace_epoch.get_untracked()
                    && documents.with_untracked(|docs| {
                        crate::protocol::plan_editing_sources::resolve_target(
                            docs,
                            target,
                            &candidate.plan.target_revision,
                        )
                        .is_ok()
                            && sources.iter().all(|s| {
                                docs.iter().any(|d| {
                                    if exact_preferences {
                                        d.id == s.event.document_id
                                            && d.revision_molecule(&s.event.event_id)
                                                .is_ok_and(|m| m.revision() == s.event.revision)
                                    } else {
                                        current_event_binding(d) == s.event
                                    }
                                })
                            })
                    });
                if !valid {
                    message.set("Pinned source or target events are unavailable or changed. Inspect the sources again before verification.".into());
                    return;
                }
            }
            let captured = if exact_preferences {
                documents.with_untracked(|docs| -> Result<(Molecule, Vec<Molecule>), String> {
                    let target_binding = &preferences
                        .as_ref()
                        .ok_or("Missing exact intended target")?
                        .2;
                    let target = crate::protocol::plan_editing_sources::resolve_target(
                        docs,
                        target_binding,
                        &candidate.plan.target_revision,
                    )?;
                    let sources = &preferences
                        .as_ref()
                        .ok_or("Missing pinned source events")?
                        .3;
                    let originals = crate::assembly::review::assembly_review::retained_inputs(
                        docs, sources, &candidate,
                    )?;
                    Ok((target.clone(), originals))
                })
            } else if is_reference {
                documents.with_untracked(|open| -> Result<(Molecule, Vec<Molecule>), String> {
                let target = open.get(active.get_untracked()).ok_or("No active target")?;
                if target.head_revision() != candidate.plan.target_revision {
                    return Err(
                        "The active target differs from this candidate's bound revision.".into(),
                    );
                }
                let bound = candidate
                    .plan
                    .inputs
                    .iter()
                    .map(|input| {
                        if let Some((_, _, _, events)) = &preferences {
                            let event = events
                                .iter()
                                .find(|e| e.input_id == input.id)
                                .ok_or("Missing exact source event preference")?;
                            if event.event.revision != input.revision {
                                return Err("Candidate/source-event revision mismatch".to_string());
                            }
                            open.iter()
                                .find(|d| current_event_binding(d) == event.event)
                                .ok_or_else(|| "The exact source event changed".to_string())
                        } else {
                            open.iter()
                                .find(|d| d.head_revision() == input.revision)
                                .ok_or_else(|| {
                                    format!("Input '{}' is missing or changed", input.id)
                                })
                        }
                    })
                    .collect::<Result<Vec<_>, String>>()?;
                let mut models = bound.iter().map(|d| &d.molecule).collect::<Vec<_>>();
                models.push(&target.molecule);
                models.extend(candidate.plan.inputs.iter().map(|i| &i.molecule));
                collagene_core::cloning::recombination::reference::model_preflight(&models)
                    .map_err(|e| e.to_string())?;
                for step in &candidate.plan.steps {
                    if let Operation::Topo{sources,..}|Operation::CheckedTopo{sources,..}|Operation::TopoActions{sources,..}=step {
                        let mut source_bases=0usize;
                        for source in sources {
                            if let Source::Input(id)=&source.source {
                                let molecule=candidate.plan.inputs.iter().position(|input|input.id==*id).and_then(|index|bound.get(index)).map(|doc|&doc.molecule).ok_or("Missing TOPO source binding")?;
                                if molecule.sequence.alphabet()!=Alphabet::Dna{return Err("TOPO sources must be DNA.".into());}
                                source_bases=source_bases.checked_add(molecule.sequence.len()).filter(|&n|n<=100_000).ok_or("TOPO source DNA exceeds its 100,000-base profile; no reference is copied.")?;
                            }
                        }
                    }
                    if let Operation::RecombineReference {
                        sources, request, ..
                    } = step
                    {
                        let resolve = |source: &Source| match source {
                            Source::Input(id) => candidate
                                .plan
                                .inputs
                                .iter()
                                .position(|i| i.id == *id)
                                .and_then(|i| bound.get(i))
                                .map(|d| &d.molecule),
                            Source::Step(_) => None,
                        };
                        if let (Some(a), Some(b)) = (resolve(&sources[0]), resolve(&sources[1])) {
                            collagene_core::cloning::recombination::source_target_preflight(
                                [a, b],
                                request,
                                &target.molecule,
                            )
                            .map_err(|e| e.to_string())?;
                        }
                    }
                }
                Ok((
                    target.molecule.clone(),
                    bound.iter().map(|d| d.molecule.clone()).collect(),
                ))
            })
            } else {
                Ok((
                    molecule.get_untracked(),
                    documents
                        .with_untracked(|docs| docs.iter().map(|d| d.molecule.clone()).collect()),
                ))
            };
            let (target, docs) = match captured {
                Ok(value) => value,
                Err(error) => {
                    message.set(error);
                    return;
                }
            };
            match candidate.bind(&target, &docs) {
                Ok(plan) => {
                    let target_event = preferences.as_ref().map(|p| p.2.clone());
                    run_plan(
                        plan,
                        target,
                        preferences.map(|(_, _, _, sources)| sources),
                        target_event,
                    )
                }
                Err(e) => {
                    plan_report.set(None);
                    last_plan.set(None);
                    last_plan_target.set(None);
                    message.set(e.to_string());
                }
            }
        }
    }
    pub(crate) fn do_plan(
        self,
        allow_nicked: RwSignal<bool>,
        closure_overlap: RwSignal<String>,
        first_enzymes: RwSignal<String>,
        first_fragment: RwSignal<String>,
        first_input: RwSignal<usize>,
        first_modifications: crate::cloning::methods::checked_cloning::AssumptionSignals,
        first_repair: RwSignal<String>,
        first_reverse: RwSignal<bool>,
        first_treatment: RwSignal<String>,
        overlap: RwSignal<String>,
        plan_method: RwSignal<&'static str>,
        restriction_final: RwSignal<&'static str>,
        restriction_policy: RwSignal<&'static str>,
        run_plan: impl Fn(
            DesignPlan,
            Molecule,
            Option<Vec<collagene_core::material::lineage::InputEventBinding>>,
            Option<collagene_core::material::lineage::EventBinding>,
        ) + Copy
        + 'static,
        second_enzymes: RwSignal<String>,
        second_fragment: RwSignal<String>,
        second_input: RwSignal<usize>,
        second_modifications: crate::cloning::methods::checked_cloning::AssumptionSignals,
        second_repair: RwSignal<String>,
        second_restriction: RwSignal<bool>,
        second_reverse: RwSignal<bool>,
        second_treatment: RwSignal<String>,
    ) -> impl Fn(web_sys::MouseEvent) + Copy + 'static {
        let session = self.session;
        let WorkspaceSession { documents, .. } = self.session;
        let ViewerSelection {
            selection,
            molecule,
            ..
        } = self.viewer;
        let ShellState { message, .. } = self.shell;
        move |_: web_sys::MouseEvent| {
            let docs = documents.get_untracked();
            let a = first_input.get_untracked();
            let b = second_input.get_untracked();
            let Some(first) = docs.get(a) else {
                message.set("Select an open first source document.".into());
                return;
            };
            let second = if plan_method.get_untracked() == "restriction"
                && (!second_restriction.get_untracked()
                    || restriction_final.get_untracked() != "ligate")
            {
                first
            } else {
                let Some(second) = docs.get(b) else {
                    message.set("Select an open second source document.".into());
                    return;
                };
                second
            };
            if plan_method.get_untracked() == "restriction" {
                let parse_index =
                    |s: String| s.parse::<usize>().ok().and_then(|i| i.checked_sub(1));
                let Some(a_index) = parse_index(first_fragment.get_untracked()) else {
                    message.set("Choose fragment numbers starting at 1.".into());
                    return;
                };
                let ligate = restriction_final.get_untracked() == "ligate";
                let use_second = second_restriction.get_untracked() && ligate;
                let b_index = if use_second {
                    let Some(i) = parse_index(second_fragment.get_untracked()) else {
                        message.set("Choose fragment numbers starting at 1.".into());
                        return;
                    };
                    i
                } else {
                    0
                };
                let checked = restriction_policy.get_untracked() == "checked";
                let names = |s: String| {
                    s.split(',')
                        .map(str::trim)
                        .filter(|s| !s.is_empty())
                        .map(str::to_owned)
                        .collect::<Vec<_>>()
                };
                let mut steps = Vec::new();
                let mut references = Vec::new();
                let mut inputs = Vec::new();
                let mut bindings = Vec::new();
                for (
                    ordinal,
                    id,
                    doc,
                    enzymes,
                    fragment_index,
                    treatment,
                    repair,
                    reverse,
                    assumptions,
                ) in [
                    (
                        0,
                        "first",
                        first,
                        first_enzymes.get_untracked(),
                        a_index,
                        first_treatment.get_untracked(),
                        first_repair.get_untracked(),
                        first_reverse.get_untracked(),
                        first_modifications.current(),
                    ),
                    (
                        1,
                        "second",
                        second,
                        second_enzymes.get_untracked(),
                        b_index,
                        second_treatment.get_untracked(),
                        second_repair.get_untracked(),
                        second_reverse.get_untracked(),
                        second_modifications.current(),
                    ),
                ] {
                    if ordinal == 1 && !use_second {
                        continue;
                    }
                    inputs.push(PlanInput {
                        id: id.into(),
                        revision: doc.molecule.revision(),
                        molecule: doc.molecule.clone(),
                    });
                    bindings.push(collagene_core::material::lineage::InputEventBinding {
                        input_id: id.into(),
                        event: collagene_core::material::lineage::EventBinding {
                            document_id: doc.id.clone(),
                            event_id: doc.head_event().into(),
                            revision: doc.molecule.revision(),
                        },
                    });
                    let source = if checked {
                        let declaration = steps.len();
                        steps.push(Operation::DeclareSourceModifications {
                            template: Source::Input(id.into()),
                            assumptions,
                        });
                        Source::Step(declaration)
                    } else {
                        Source::Input(id.into())
                    };
                    let digest_index = steps.len();
                    steps.push(if checked {
                    let selection = match (|| -> Result<_, String> {
                            let catalog=collagene_core::cloning::restriction::catalog::builtin_catalog().map_err(|e|e.to_string())?;
                            let ids=names(enzymes).into_iter().map(|name|{
                                let mut found=catalog.edition().enzymes.iter().filter(|e|e.name==name);
                                let record=found.next().ok_or_else(||format!("Choose an exact catalog enzyme name; {name} is unavailable."))?;
                                if found.next().is_some(){return Err(format!("{name} is ambiguous; choose a unique canonical formulation."));}Ok(record.id.clone())
                            }).collect::<Result<Vec<_>,String>>()?;
                            collagene_core::cloning::restriction::catalog::select(&catalog,&ids).map_err(|e|e.to_string())
                    })() {
                        Ok(selection) => selection,
                        Err(error) => { message.set(error); return; }
                    };
                    Operation::CheckedDigest {
                        template: source,
                        selection,
                        options:collagene_core::chemistry::tracked_duplex::checked_inspection_options(),
                        fragment_index,
                    }
                } else {
                    Operation::Digest {
                        template: source,
                        enzymes: names(enzymes),
                        fragment_index,
                    }
                });
                    let mut product_index = digest_index;
                    if repair != "none" {
                        let i = steps.len();
                        steps.push(Operation::RepairEnds {
                            template: Source::Step(product_index),
                            method: match repair.as_str() {
                                "fill" => chemistry::EndRepairMethod::FillFivePrimeOverhangs,
                                "chew" => chemistry::EndRepairMethod::RemoveThreePrimeOverhangs,
                                _ => chemistry::EndRepairMethod::Blunt,
                            },
                        });
                        product_index = i;
                    }
                    if treatment != "none" {
                        let i = steps.len();
                        steps.push(Operation::TreatEnds {
                            template: Source::Step(product_index),
                            treatment: if treatment == "phosphorylate" {
                                chemistry::EndTreatment::Phosphorylate
                            } else {
                                chemistry::EndTreatment::Dephosphorylate
                            },
                        });
                        product_index = i;
                    }
                    references.push(collagene_core::planning::plan::OrientedSource {
                        source: Source::Step(product_index),
                        reverse,
                    });
                }
                let target = molecule.get_untracked();
                if ligate {
                    steps.push(collagene_core::cloning::ligation::outcomes::operation(
                        references,
                        target.topology == Topology::Circular,
                        allow_nicked.get_untracked(),
                    ));
                }
                let plan = DesignPlan {
                    name: format!(
                        "{} restriction {} to {}",
                        if checked {
                            "Checked"
                        } else {
                            "Assumption-based"
                        },
                        if ligate { "ligation" } else { "fragment" },
                        target.name
                    ),
                    target_revision: target.revision(),
                    inputs,
                    steps,
                };
                run_plan(plan, target, Some(bindings), None);
                return;
            }
            let Ok(overlap) = overlap.get_untracked().parse::<usize>() else {
                message.set("Overlap must be a positive integer.".into());
                return;
            };
            let target = molecule.get_untracked();
            let circular = target.topology == Topology::Circular;
            let mut overlaps = vec![overlap];
            if circular {
                let Ok(close) = closure_overlap.get_untracked().parse::<usize>() else {
                    message.set("Closure overlap must be a positive integer.".into());
                    return;
                };
                overlaps.push(close);
            }
            let plan = DesignPlan {
                name: format!("Assembly to {}", target.name),
                target_revision: target.revision(),
                inputs: vec![
                    PlanInput {
                        id: "first".into(),
                        revision: first.molecule.revision(),
                        molecule: first.molecule.clone(),
                    },
                    PlanInput {
                        id: "second".into(),
                        revision: second.molecule.revision(),
                        molecule: second.molecule.clone(),
                    },
                ],
                steps: vec![Operation::Overlap {
                    fragments: vec![
                        Source::Input("first".into()),
                        Source::Input("second".into()),
                    ],
                    overlaps,
                    circular,
                }],
            };
            let events = [("first", first), ("second", second)]
                .into_iter()
                .map(
                    |(id, doc)| collagene_core::material::lineage::InputEventBinding {
                        input_id: id.into(),
                        event: collagene_core::material::lineage::EventBinding {
                            document_id: doc.id.clone(),
                            event_id: doc.head_event().into(),
                            revision: doc.molecule.revision(),
                        },
                    },
                )
                .collect();
            run_plan(plan, target, Some(events), None);
        }
    }
    pub(crate) fn finalize_assembly(
        &self,
        candidate: RwSignal<Option<transfer::Candidate>>,
        pending_finalize: RwSignal<Option<std::sync::Arc<PendingFinalize>>>,
        stored_candidates: RwSignal<Vec<transfer::Candidate>>,
    ) -> Callback<crate::assembly::draft::assembly_workspace::ReviewRequest, Result<(), String>>
    {
        let session = self.session;
        let WorkspaceSession {
            documents,
            active,
            workspace_restoring,
            workspace_name,
            attachments,
            ..
        } = self.session;
        let ShellState { modal, .. } = self.shell;
        Callback::new(
        move |request: crate::assembly::draft::assembly_workspace::ReviewRequest| -> Result<(), String> {
            if workspace_restoring.get_untracked() {
                return Err("Wait for the workspace operation to finish.".into());
            }
            let frame = documents.with_untracked(|docs| {
                WorkspaceSession::snapshot(
                    docs,
                    active.get_untracked(),
                    candidate.get_untracked().as_ref(),
                    &stored_candidates.get_untracked(),
                    &attachments.get_untracked(),
                    &workspace_name.get_untracked(),
                )
            })?;
            let fingerprint = collagene_core::material::lineage::workspace_fingerprint(&frame)
                .map_err(|e| e.to_string())?;
            let staged = crate::assembly::review::assembly_review::stage_any(
                &frame,
                &request.intent,
                &request.sources,
                &request.candidate,
            )?;
            RestoredWorkspace::from_frame(&staged.workspace)?;
            pending_finalize.set(Some(std::sync::Arc::new(PendingFinalize {
                workspace: std::sync::Arc::new(staged.workspace),
                plan: request.candidate.plan,
                target: staged.target,
                sources: staged.sources,
                source_fingerprint: fingerprint,
            })));
            modal.set("assembly_finalize");
            Ok(())
        },
    )
    }
    pub(crate) fn finalize_reviewed_candidate(
        &self,
        candidate: RwSignal<Option<transfer::Candidate>>,
        capture_junction_plan: Callback<(), Result<(DesignPlan, Molecule), String>>,
        last_plan_events: RwSignal<
            Option<(
                String,
                Vec<collagene_core::material::lineage::InputEventBinding>,
            )>,
        >,
        last_plan_target: RwSignal<
            Option<(String, u64, collagene_core::material::lineage::EventBinding)>,
        >,
        pending_finalize: RwSignal<Option<std::sync::Arc<PendingFinalize>>>,
        stored_candidates: RwSignal<Vec<transfer::Candidate>>,
    ) -> Callback<transfer::Candidate, Result<(), String>> {
        let session = self.session;
        let WorkspaceSession {
            documents,
            active,
            workspace_epoch,
            workspace_name,
            attachments,
            ..
        } = self.session;
        let ShellState { modal, .. } = self.shell;
        Callback::new(move |value: transfer::Candidate| -> Result<(), String> {
            let (plan, _) = capture_junction_plan.run(())?;
            if value.plan_revision != collagene_core::planning::plan::DesignPlan::revision(&plan) {
                return Err("The reviewed reconstruction changed.".into());
            }
            let (_, epoch, target) = last_plan_target
                .get_untracked()
                .ok_or("Missing exact target event")?;
            if epoch != workspace_epoch.get_untracked() {
                return Err("The workspace changed; reconstruct again.".into());
            }
            let sources=last_plan_events.get_untracked().filter(|(revision,_)|*revision==value.plan_revision)
            .map(|(_,sources)|sources).ok_or("Exact source events are unavailable. Review a source-bound design before creating outputs.")?;
            let frame = documents.with_untracked(|docs| {
                WorkspaceSession::snapshot(
                    docs,
                    active.get_untracked(),
                    candidate.get_untracked().as_ref(),
                    &stored_candidates.get_untracked(),
                    &attachments.get_untracked(),
                    &workspace_name.get_untracked(),
                )
            })?;
            let fingerprint = collagene_core::material::lineage::workspace_fingerprint(&frame)
                .map_err(|e| e.to_string())?;
            let staged = crate::assembly::review::assembly_review::stage_existing(
                &frame, &value, &target, &sources,
            )?;
            RestoredWorkspace::from_frame(&staged.workspace)?;
            pending_finalize.set(Some(std::sync::Arc::new(PendingFinalize {
                workspace: std::sync::Arc::new(staged.workspace),
                plan: value.plan,
                target: staged.target,
                sources: staged.sources,
                source_fingerprint: fingerprint,
            })));
            modal.set("assembly_finalize");
            Ok(())
        })
    }
    pub(crate) fn review_assembly(
        &self,
        assembly_visible: RwSignal<bool>,
        candidate: RwSignal<Option<transfer::Candidate>>,
        candidate_events: RwSignal<
            Option<(
                String,
                u64,
                collagene_core::material::lineage::EventBinding,
                Vec<collagene_core::material::lineage::InputEventBinding>,
            )>,
        >,
        last_plan: RwSignal<Option<(DesignPlan, Molecule)>>,
        last_plan_events: RwSignal<
            Option<(
                String,
                Vec<collagene_core::material::lineage::InputEventBinding>,
            )>,
        >,
        last_plan_target: RwSignal<
            Option<(String, u64, collagene_core::material::lineage::EventBinding)>,
        >,
        plan_generation: RwSignal<u32>,
        plan_report: RwSignal<Option<VerificationReport>>,
        stored_candidates: RwSignal<Vec<transfer::Candidate>>,
    ) -> Callback<crate::assembly::draft::assembly_workspace::ReviewRequest, Result<(), String>>
    {
        let session = self.session;
        let WorkspaceSession {
            documents,
            active,
            workspace_restoring,
            workspace_epoch,
            workspace_name,
            attachments,
            ..
        } = self.session;
        let ShellState { panel, message, .. } = self.shell;
        Callback::new(
        move |request: crate::assembly::draft::assembly_workspace::ReviewRequest| -> Result<(), String> {
            if workspace_restoring.get_untracked() {
                return Err("Wait for the workspace operation to finish.".into());
            }
            let frame = documents.with_untracked(|docs| {
                WorkspaceSession::snapshot(
                    docs,
                    active.get_untracked(),
                    candidate.get_untracked().as_ref(),
                    &stored_candidates.get_untracked(),
                    &attachments.get_untracked(),
                    &workspace_name.get_untracked(),
                )
            })?;
            let staged = crate::assembly::review::assembly_review::stage_any(
                &frame,
                &request.intent,
                &request.sources,
                &request.candidate,
            )?;
            let state = RestoredWorkspace::from_frame(&staged.workspace)?;
            let preferred = (
                request.candidate.plan_revision.clone(),
                workspace_epoch.get_untracked(),
                staged.target,
                staged.sources,
            );
            batch(|| {
                documents.set(state.documents);
                active.set(state.active);
                stored_candidates.set(state.candidates);
                candidate_events.set(None);
                candidate.set(Some(request.candidate));
                candidate_events.set(Some(preferred));
                plan_report.set(None);
                last_plan.set(None);
                last_plan_target.set(None);
                last_plan_events.set(None);
                plan_generation.update(|g| *g = g.wrapping_add(1));
                assembly_visible.set(false);
                panel.set("plan");
                message.set("Opened the intended reference and proposed preparation. Verify the assembly to reconstruct it from the pinned source events.".into());
            });
            Ok(())
        },
    )
    }
    pub(crate) fn navigate_assembly(
        &self,
        candidate: RwSignal<Option<transfer::Candidate>>,
        stored_candidates: RwSignal<Vec<transfer::Candidate>>,
    ) -> Callback<crate::assembly::draft::assembly_workspace::SelectionNavigation, ()> {
        let session = self.session;
        let WorkspaceSession {
            documents,
            active,
            workspace_restoring,
            workspace_name,
            attachments,
            ..
        } = self.session;
        let ViewerSelection {
            selection,
            selected_feature,
            search_regions,
            molecule,
            ..
        } = self.viewer;
        let ShellState { message, .. } = self.shell;
        Callback::new(
            move |request: crate::assembly::draft::assembly_workspace::SelectionNavigation| {
                if workspace_restoring.get_untracked() {
                    return;
                }
                let opened = documents.with_untracked(
                    |docs| -> Result<(Vec<OpenDocument>, usize), String> {
                        let index = docs
                            .iter()
                            .position(|d| d.id == request.source.document_id)
                            .ok_or("The pinned source document is unavailable.")?;
                        let source = &docs[index];
                        let exact = source
                            .revision_molecule(&request.source.event_id)
                            .map_err(|e| e.to_string())?;
                        if exact.revision() != request.source.revision {
                            return Err(
                                "The source event no longer matches the selected fragment.".into(),
                            );
                        }
                        if current_event_binding(source) == request.source {
                            return Ok((docs.clone(), index));
                        }
                        if docs.len() >= 1000 {
                            return Err("The workspace document limit has been reached.".into());
                        }
                        let prefix = format!("assembly-source-{}", &request.source.revision[..16]);
                        let id = (1..=1001)
                            .map(|n| format!("{prefix}-{n}"))
                            .find(|id| !docs.iter().any(|d| d.id == *id))
                            .ok_or("No source document identity is available.")?;
                        let fork = source
                            .fork_event(&request.source.event_id, id)
                            .map_err(|e| e.to_string())?;
                        let mut staged = docs.clone();
                        let index = staged.len();
                        staged.push(fork);
                        WorkspaceSession::snapshot(
                            &staged,
                            index,
                            candidate.get_untracked().as_ref(),
                            &stored_candidates.get_untracked(),
                            &attachments.get_untracked(),
                            &workspace_name.get_untracked(),
                        )?
                        .validate()
                        .map_err(|e| e.to_string())?;
                        Ok((staged, index))
                    },
                );
                let (staged, index) = match opened {
                    Ok(v) => v,
                    Err(error) => {
                        message.set(error);
                        return;
                    }
                };
                let expected_document = staged[index].id.clone();
                let expected_event = staged[index].head_event().to_owned();
                let expected_revision = staged[index].head_revision().to_owned();
                batch(|| {
                    documents.set(staged);
                    active.set(index);
                });
                let apply = wasm_bindgen::closure::Closure::once_into_js(move || {
                    if documents.with_untracked(|docs| {
                        docs.get(active.get_untracked()).is_none_or(|d| {
                            d.id != expected_document
                                || d.head_event() != expected_event
                                || d.head_revision() != expected_revision
                        })
                    }) {
                        return;
                    }
                    let n = molecule.with_untracked(|m| m.sequence.len());
                    let topology = molecule.with_untracked(|m| m.topology);
                    let mut parts = Vec::new();
                    for (start, end, _) in &request.regions {
                        if n == 0
                            || *start >= n
                            || end <= start
                            || end - start > n
                            || (topology == Topology::Linear && *end > n)
                        {
                            message.set("The source projection lies outside the displayed reference. Review the strand contribution again.".into());
                            return;
                        }
                        if *end <= n {
                            parts.push((*start, *end));
                        } else {
                            parts.push((*start, n));
                            parts.push((0, *end - n));
                        }
                    }
                    selected_feature.set(None);
                    selection.set(if parts.len() == 1 {
                        Some(parts[0])
                    } else {
                        None
                    });
                    search_regions.set(if parts.is_empty() {
                        None
                    } else {
                        Some((expected_revision, parts))
                    });
                    message.set(if request.regions.is_empty(){format!("{} · These actual strand bases have no displayed top-reference positions. Inspect the strand sequence in the assembly source details.",request.description)}else{format!("{} · Highlighted only mapped original top-reference bases; actual strand-only bases remain in assembly source details.",request.description)});
                });
                if let Some(window) = web_sys::window() {
                    let _ = window.request_animation_frame(apply.unchecked_ref());
                }
            },
        )
    }
    pub(crate) fn open_mutation_intent(
        &self,
        contig_snapshot: Callback<
            (),
            Result<collagene_core::documents::workspace::Workspace, String>,
        >,
        last_plan: RwSignal<Option<(DesignPlan, Molecule)>>,
        last_plan_target: RwSignal<
            Option<(String, u64, collagene_core::material::lineage::EventBinding)>,
        >,
        plan_report: RwSignal<Option<VerificationReport>>,
    ) -> Callback<crate::cloning::methods::mutagenesis::IntentOpen, Result<(), String>> {
        let session = self.session;
        let WorkspaceSession {
            documents,
            active,
            workspace_restoring,
            attachments,
            ..
        } = self.session;
        let ViewerSelection {
            selection,
            selected_feature,
            ..
        } = self.viewer;
        let ShellState { message, .. } = self.shell;
        Callback::new(
            move |intent: crate::cloning::methods::mutagenesis::IntentOpen| -> Result<(), String> {
                if workspace_restoring.get_untracked() {
                    return Err("Wait for the workspace operation to finish.".into());
                }
                let frame = contig_snapshot.run(())?;
                if !documents.with_untracked(|docs| {
                    docs.iter()
                        .any(|doc| current_event_binding(doc) == intent.source)
                }) {
                    return Err(
                        "The mutation source event changed; preview the intended edit again."
                            .into(),
                    );
                }
                let staged = crate::cloning::methods::mutation_review::stage(
                    &frame,
                    &intent.source,
                    &intent.report,
                )?;
                RestoredWorkspace::from_frame(&staged)?;
                let restored = OpenDocument::restore_all(&staged).map_err(|e| e.to_string())?;
                let index = restored
                    .iter()
                    .position(|doc| doc.id == staged.active_document)
                    .ok_or("The intended reference is unavailable.")?;
                batch(|| {
                    documents.set(restored);
                    active.set(index);
                    attachments.set(staged.attachments);
                    selection.set(None);
                    selected_feature.set(None);
                    plan_report.set(None);
                    last_plan.set(None);
                    last_plan_target.set(None);
                    message.set("Opened the declared intended reference with its exact source history. Design and verify a preparation against this target next.".into());
                });
                Ok(())
            },
        )
    }
    pub(crate) fn open_domestication(
        &self,
        candidate: RwSignal<Option<transfer::Candidate>>,
        candidate_events: RwSignal<
            Option<(
                String,
                u64,
                collagene_core::material::lineage::EventBinding,
                Vec<collagene_core::material::lineage::InputEventBinding>,
            )>,
        >,
        contig_snapshot: Callback<
            (),
            Result<collagene_core::documents::workspace::Workspace, String>,
        >,
        last_plan: RwSignal<Option<(DesignPlan, Molecule)>>,
        last_plan_events: RwSignal<
            Option<(
                String,
                Vec<collagene_core::material::lineage::InputEventBinding>,
            )>,
        >,
        last_plan_target: RwSignal<
            Option<(String, u64, collagene_core::material::lineage::EventBinding)>,
        >,
        plan_generation: RwSignal<u32>,
        plan_report: RwSignal<Option<VerificationReport>>,
        stored_candidates: RwSignal<Vec<transfer::Candidate>>,
    ) -> Callback<crate::cloning::methods::domestication_review::Open, Result<(), String>> {
        let session = self.session;
        let WorkspaceSession {
            documents,
            active,
            workspace_restoring,
            workspace_epoch,
            attachments,
            ..
        } = self.session;
        let ViewerSelection {
            selection,
            selected_feature,
            search_regions,
            ..
        } = self.viewer;
        let ShellState { panel, message, .. } = self.shell;
        Callback::new(
        move |opened: crate::cloning::methods::domestication_review::Open| -> Result<(), String> {
            if workspace_restoring.get_untracked()
                || opened.epoch != workspace_epoch.get_untracked()
            {
                return Err(
                    "The domestication workspace changed; inspect the current source again.".into(),
                );
            }
            if !documents.with_untracked(|docs| {
                docs.get(active.get_untracked())
                    .is_some_and(|doc| current_event_binding(doc) == opened.source)
            }) {
                return Err("The exact original source event changed.".into());
            }
            let staged = crate::cloning::methods::domestication_review::stage(&contig_snapshot.run(())?, &opened)?;
            RestoredWorkspace::from_frame(&staged.workspace)?;
            let restored =
                OpenDocument::restore_all(&staged.workspace).map_err(|e| e.to_string())?;
            let index = restored
                .iter()
                .position(|doc| doc.id == staged.target.document_id)
                .ok_or("The independently checked target was not retained.")?;
            let preference = staged.candidate.as_ref().map(|value| {
                (
                    value.plan_revision.clone(),
                    opened.epoch,
                    staged.target,
                    staged.inputs,
                )
            });
            batch(|| {
                documents.set(restored);
                active.set(index);
                attachments.set(staged.workspace.attachments);
                stored_candidates.set(staged.workspace.candidates);
                candidate_events.set(None);
                candidate.set(staged.candidate);
                candidate_events.set(preference);
                selection.set(None);
                selected_feature.set(None);
                search_regions.set(None);
                plan_report.set(None);
                last_plan.set(None);
                last_plan_target.set(None);
                last_plan_events.set(None);
                plan_generation.update(|g| *g = g.wrapping_add(1));
                if opened.prepared.is_some() {
                    panel.set("plan");
                    message.set("Opened the original-source mutagenic preparation and independently checked intended reference. Verify the plan, inspect actual products, then choose outputs.".into());
                } else {
                    message.set("Opened checked substitution intent with unchanged annotation coordinates and exact source ancestry. This reference is not prepared DNA.".into());
                }
            });
            Ok(())
        },
    )
    }
}
