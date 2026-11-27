//! Atomic installation of an intended reference and its proposed preparation.
//! Exact source events are checked; no scientific success is restored or inferred.
use collagene_core::{
    assembly_intent::{self, Report},
    history::DocumentHistory,
    run_lineage::{EventBinding, InputEventBinding},
    transfer::Candidate,
    workspace::Workspace,
};
/// Independently declared reference intent, never a target copied from a simulator.
#[derive(Clone)]
pub enum IntendedReference {
    Regions(assembly_intent::Report),
    Junctions(collagene_core::assembly_junction_edit::Intent),
    Substitutions(collagene_core::assembly_substitution_intent::Report),
    Edits(collagene_core::assembly_edit_intent::Report),
}
impl IntendedReference {
    pub fn target(&self) -> &collagene_core::Molecule {
        match self {
            Self::Regions(intent) => &intent.target,
            Self::Junctions(intent) => &intent.target,
            Self::Substitutions(intent) => &intent.target,
            Self::Edits(intent) => &intent.target,
        }
    }
    pub fn fragments(&self) -> &[assembly_intent::PlacedFragment] {
        match self {
            Self::Regions(intent) => &intent.fragments,
            Self::Junctions(intent) => &intent.fragments,
            Self::Substitutions(intent) => &intent.fragments,
            Self::Edits(intent) => &intent.fragments,
        }
    }
}
pub fn stage_any(
    frame: &Workspace,
    intent: &IntendedReference,
    sources: &[EventBinding],
    candidate: &Candidate,
) -> Result<Staged, String> {
    match intent {
        IntendedReference::Regions(intent) => stage(frame, intent, sources, candidate),
        IntendedReference::Junctions(intent) => stage_edited(frame, intent, sources, candidate),
        IntendedReference::Substitutions(intent) => {
            stage_substituted(frame, intent, sources, candidate)
        }
        IntendedReference::Edits(intent) => {
            stage_contribution_edits(frame, intent, sources, candidate)
        }
    }
}
pub struct Staged {
    pub workspace: Workspace,
    pub target: EventBinding,
    pub sources: Vec<InputEventBinding>,
}
/// Stage an existing independently specified target and explicitly pinned inputs.
/// No target is constructed, source head rewound or scientific check restored.
/// Historical source events are valid only when supplied with their exact IDs;
/// equal content in another document is never used as a replacement.
pub fn stage_existing(
    frame: &Workspace,
    candidate: &Candidate,
    target: &EventBinding,
    sources: &[InputEventBinding],
) -> Result<Staged, String> {
    frame.validate().map_err(|e| e.to_string())?;
    candidate.validate().map_err(|e| e.to_string())?;
    if sources.len() != candidate.plan.inputs.len() || sources.len() > 128 {
        return Err("Pinned input count differs from the candidate output profile.".into());
    }
    let docs = DocumentHistory::restore_all(frame).map_err(|e| e.to_string())?;
    crate::plan_editing_sources::resolve_target(&docs, target, &candidate.plan.target_revision)?;
    for (input, pinned) in candidate.plan.inputs.iter().zip(sources) {
        if input.id != pinned.input_id {
            return Err("Ordered pinned input identities differ from the candidate.".into());
        }
        let source = docs
            .iter()
            .find(|doc| doc.id == pinned.event.document_id)
            .ok_or("A pinned source document is unavailable.")?;
        let exact = source
            .revision_molecule(&pinned.event.event_id)
            .map_err(|e| e.to_string())?;
        if exact.revision() != pinned.event.revision
            || pinned.event.revision != input.revision
            || exact != input.molecule
        {
            return Err("An exact source event differs from its candidate input.".into());
        }
    }
    let mut staged = frame.clone();
    if !staged
        .candidates
        .iter()
        .any(|existing| existing.plan_revision == candidate.plan_revision)
    {
        staged.candidates.push(candidate.clone());
    }
    staged.active_candidate = Some(candidate.plan_revision.clone());
    staged.validate().map_err(|e| e.to_string())?;
    Ok(Staged {
        workspace: staged,
        target: target.clone(),
        sources: sources.to_vec(),
    })
}
pub fn retained_inputs(
    docs: &[DocumentHistory],
    sources: &[InputEventBinding],
    candidate: &Candidate,
) -> Result<Vec<collagene_core::Molecule>, String> {
    if sources.len() != candidate.plan.inputs.len() {
        return Err("Pinned input count differs from the candidate.".into());
    }
    candidate
        .plan
        .inputs
        .iter()
        .map(|input| {
            let pinned = sources
                .iter()
                .find(|source| source.input_id == input.id)
                .ok_or("A pinned input event is missing.")?;
            let doc = docs
                .iter()
                .find(|d| d.id == pinned.event.document_id)
                .ok_or("The pinned source document is unavailable.")?;
            let molecule = doc
                .revision_molecule(&pinned.event.event_id)
                .map_err(|e| e.to_string())?;
            if molecule.revision() != pinned.event.revision
                || pinned.event.revision != input.revision
            {
                return Err("Pinned source content differs from the preparation input.".into());
            }
            Ok(molecule)
        })
        .collect()
}
pub fn stage(
    frame: &Workspace,
    intent: &Report,
    sources: &[EventBinding],
    candidate: &Candidate,
) -> Result<Staged, String> {
    frame.validate().map_err(|e| e.to_string())?;
    candidate.validate().map_err(|e| e.to_string())?;
    if sources.len() != candidate.plan.inputs.len() || sources.is_empty() || sources.len() > 32 {
        return Err("Assembly source bindings do not match the candidate.".into());
    }
    let docs = DocumentHistory::restore_all(frame).map_err(|e| e.to_string())?;
    let mut parents = Vec::with_capacity(sources.len());
    for (binding, input) in sources.iter().zip(&candidate.plan.inputs) {
        let doc = docs
            .iter()
            .find(|d| d.id == binding.document_id)
            .ok_or("An assembly source document is unavailable.")?;
        let exact = doc
            .revision_molecule(&binding.event_id)
            .map_err(|e| e.to_string())?;
        if exact.revision() != binding.revision
            || binding.revision != input.revision
            || exact.revision() != input.molecule.revision()
        {
            return Err("A pinned assembly source does not match its preparation input.".into());
        }
        parents.push((doc, binding.event_id.as_str()));
    }
    let fresh = assembly_intent::build(&candidate.plan.inputs, &intent.request, Some(frame))
        .map_err(|e| e.to_string())?;
    intent
        .validate_binding(&fresh.context_revision)
        .map_err(|e| e.to_string())?;
    if intent.target.revision() != fresh.target.revision()
        || candidate.plan.target_revision != fresh.target.revision()
    {
        return Err(
            "The candidate does not target the independently declared source regions.".into(),
        );
    }
    let expected_target_revision = fresh.target.revision();
    let expected_parents = sources
        .iter()
        .map(|s| s.event_id.as_str())
        .collect::<std::collections::BTreeSet<_>>();
    let existing = docs.iter().find(|doc| {
        doc.head_revision() == expected_target_revision
            && doc
                .summaries()
                .iter()
                .find(|event| event.current)
                .is_some_and(|event| {
                    event
                        .parents
                        .iter()
                        .map(String::as_str)
                        .collect::<std::collections::BTreeSet<_>>()
                        == expected_parents
                })
    });
    let mut staged = frame.clone();
    let target = if let Some(existing) = existing {
        EventBinding {
            document_id: existing.id.clone(),
            event_id: existing.head_event().into(),
            revision: existing.head_revision().into(),
        }
    } else {
        if docs.len() >= 1000 {
            return Err("The workspace document limit has been reached.".into());
        }
        let prefix = format!("assembly-intent-{}", &fresh.context_revision[..16]);
        let id = (1..=1001)
            .map(|n| format!("{prefix}-{n}"))
            .find(|id| !docs.iter().any(|d| d.id == *id))
            .ok_or("No assembly document identity is available.")?;
        let derived=DocumentHistory::derive_from_events(&parents,id,fresh.target,"Declared assembly design intent from pinned source regions; expected reference only, preparation and experimental outcome unverified").map_err(|e|e.to_string())?;
        let target = EventBinding {
            document_id: derived.id.clone(),
            event_id: derived.head_event().into(),
            revision: derived.head_revision().into(),
        };
        derived.contribute(&mut staged).map_err(|e| e.to_string())?;
        target
    };
    staged.active_document = target.document_id.clone();
    if !staged
        .candidates
        .iter()
        .any(|c| c.plan_revision == candidate.plan_revision)
    {
        staged.candidates.push(candidate.clone());
    }
    staged.active_candidate = Some(candidate.plan_revision.clone());
    staged.validate().map_err(|e| e.to_string())?;
    Ok(Staged {
        workspace: staged,
        target,
        sources: sources
            .iter()
            .zip(&candidate.plan.inputs)
            .map(|(event, input)| InputEventBinding {
                input_id: input.id.clone(),
                event: event.clone(),
            })
            .collect(),
    })
}

fn stage_edited(
    frame: &Workspace,
    intent: &collagene_core::assembly_junction_edit::Intent,
    sources: &[EventBinding],
    candidate: &Candidate,
) -> Result<Staged, String> {
    if candidate.plan.inputs.is_empty() || candidate.plan.inputs.len() > 32 {
        return Err("Edited assembly inputs exceed the bounded staging profile.".into());
    }
    collagene_core::assembly_junction_edit::intent_preflight(
        &candidate
            .plan
            .inputs
            .iter()
            .map(|input| (input.id.as_str(), &input.molecule))
            .collect::<Vec<_>>(),
        &intent.base_request,
        &intent.edits,
        Some(frame),
    )
    .map_err(|e| e.to_string())?;
    frame.validate().map_err(|e| e.to_string())?;
    candidate.validate().map_err(|e| e.to_string())?;
    if sources.len() != candidate.plan.inputs.len() || sources.is_empty() || sources.len() > 32 {
        return Err("Assembly source bindings do not match the candidate.".into());
    }
    let docs = DocumentHistory::restore_all(frame).map_err(|e| e.to_string())?;
    let mut parents = Vec::with_capacity(sources.len());
    for (binding, input) in sources.iter().zip(&candidate.plan.inputs) {
        let doc = docs
            .iter()
            .find(|d| d.id == binding.document_id)
            .ok_or("An assembly source document is unavailable.")?;
        let exact = doc
            .revision_molecule(&binding.event_id)
            .map_err(|e| e.to_string())?;
        if exact.revision() != binding.revision
            || binding.revision != input.revision
            || exact.revision() != input.molecule.revision()
        {
            return Err("A pinned assembly source does not match its preparation input.".into());
        }
        parents.push((doc, binding.event_id.as_str()));
    }
    let fresh = collagene_core::assembly_junction_edit::build_intent(
        &candidate.plan.inputs,
        &intent.base_request,
        &intent.edits,
        Some(frame),
    )
    .map_err(|e| e.to_string())?;
    intent
        .validate_binding(&fresh.context_revision)
        .map_err(|e| e.to_string())?;
    if intent.target.revision() != fresh.target.revision()
        || candidate.plan.target_revision != fresh.target.revision()
    {
        return Err(
            "The candidate does not target the independently declared source regions and junction edits.".into(),
        );
    }
    let expected_target_revision = fresh.target.revision();
    let expected_parents = sources
        .iter()
        .map(|s| s.event_id.as_str())
        .collect::<std::collections::BTreeSet<_>>();
    let existing = docs.iter().find(|doc| {
        doc.head_revision() == expected_target_revision
            && doc
                .summaries()
                .iter()
                .find(|event| event.current)
                .is_some_and(|event| {
                    event
                        .parents
                        .iter()
                        .map(String::as_str)
                        .collect::<std::collections::BTreeSet<_>>()
                        == expected_parents
                })
    });
    let mut staged = frame.clone();
    let target = if let Some(existing) = existing {
        EventBinding {
            document_id: existing.id.clone(),
            event_id: existing.head_event().into(),
            revision: existing.head_revision().into(),
        }
    } else {
        if docs.len() >= 1000 {
            return Err("The workspace document limit has been reached.".into());
        }
        let prefix = format!("assembly-junction-intent-{}", &fresh.context_revision[..16]);
        let id = (1..=1001)
            .map(|n| format!("{prefix}-{n}"))
            .find(|id| !docs.iter().any(|d| d.id == *id))
            .ok_or("No assembly document identity is available.")?;
        let derived=DocumentHistory::derive_from_events(&parents,id,fresh.target,"Declared assembly design intent from pinned source regions and explicit junction edits; synthetic junction bases are design intent, not observed DNA or verified physical preparation").map_err(|e|e.to_string())?;
        let target = EventBinding {
            document_id: derived.id.clone(),
            event_id: derived.head_event().into(),
            revision: derived.head_revision().into(),
        };
        derived.contribute(&mut staged).map_err(|e| e.to_string())?;
        target
    };
    staged.active_document = target.document_id.clone();
    if !staged
        .candidates
        .iter()
        .any(|c| c.plan_revision == candidate.plan_revision)
    {
        staged.candidates.push(candidate.clone());
    }
    staged.active_candidate = Some(candidate.plan_revision.clone());
    staged.validate().map_err(|e| e.to_string())?;
    Ok(Staged {
        workspace: staged,
        target,
        sources: sources
            .iter()
            .zip(&candidate.plan.inputs)
            .map(|(event, input)| InputEventBinding {
                input_id: input.id.clone(),
                event: event.clone(),
            })
            .collect(),
    })
}

pub fn stage_substituted(
    frame: &Workspace,
    intent: &collagene_core::assembly_substitution_intent::Report,
    sources: &[EventBinding],
    candidate: &Candidate,
) -> Result<Staged, String> {
    stage_internal_edits(frame, sources, candidate, || {
        let fresh = collagene_core::assembly_substitution_intent::build(
            &candidate.plan.inputs,
            &intent.request,
            Some(frame),
        )
        .map_err(|e| e.to_string())?;
        intent
            .validate_binding(&fresh.context_revision)
            .map_err(|e| e.to_string())?;
        if intent.target.revision() != fresh.target.revision() {
            return Err(
                "Declared substitution target differs from fresh source-bound intent.".into(),
            );
        }
        Ok((fresh.context_revision, fresh.target))
    })
}

pub fn stage_contribution_edits(
    frame: &Workspace,
    intent: &collagene_core::assembly_edit_intent::Report,
    sources: &[EventBinding],
    candidate: &Candidate,
) -> Result<Staged, String> {
    stage_internal_edits(frame, sources, candidate, || {
        let expected = collagene_core::assembly_edit_intent::intent_revision(
            &candidate.plan.inputs,
            &intent.request,
            Some(frame),
        )
        .map_err(|e| e.to_string())?;
        intent
            .validate_sources(&candidate.plan.inputs, &expected, Some(frame))
            .map_err(|e| e.to_string())?;
        Ok((expected, intent.target.clone()))
    })
}

/// Both declarations use the same exact-event and atomic workspace installation.
/// The supplied builder must reconstruct intent from the candidate's original inputs.
fn stage_internal_edits(
    frame: &Workspace,
    sources: &[EventBinding],
    candidate: &Candidate,
    build: impl FnOnce() -> Result<(String, collagene_core::Molecule), String>,
) -> Result<Staged, String> {
    if candidate.plan.inputs.is_empty() || candidate.plan.inputs.len() > 32 {
        return Err("Edited assembly inputs exceed the bounded staging profile.".into());
    }
    frame.validate().map_err(|e| e.to_string())?;
    candidate.validate().map_err(|e| e.to_string())?;
    if sources.len() != candidate.plan.inputs.len() || sources.is_empty() || sources.len() > 32 {
        return Err("Assembly source bindings do not match the candidate.".into());
    }
    let docs = DocumentHistory::restore_all(frame).map_err(|e| e.to_string())?;
    let mut parents = Vec::with_capacity(sources.len());
    for (binding, input) in sources.iter().zip(&candidate.plan.inputs) {
        let doc = docs
            .iter()
            .find(|d| d.id == binding.document_id)
            .ok_or("An assembly source document is unavailable.")?;
        let exact = doc
            .revision_molecule(&binding.event_id)
            .map_err(|e| e.to_string())?;
        if exact.revision() != binding.revision
            || binding.revision != input.revision
            || exact.revision() != input.molecule.revision()
        {
            return Err("A pinned assembly source does not match its preparation input.".into());
        }
        parents.push((doc, binding.event_id.as_str()));
    }
    let (fresh_context, fresh_target) = build()?;
    if candidate.plan.target_revision != fresh_target.revision() {
        return Err("The candidate does not target the independently declared source regions and literal edits.".into());
    }
    let expected_target_revision = fresh_target.revision();
    let expected_parents = sources
        .iter()
        .map(|s| s.event_id.as_str())
        .collect::<std::collections::BTreeSet<_>>();
    let existing = docs.iter().find(|doc| {
        doc.head_revision() == expected_target_revision
            && doc
                .summaries()
                .iter()
                .find(|event| event.current)
                .is_some_and(|event| {
                    event
                        .parents
                        .iter()
                        .map(String::as_str)
                        .collect::<std::collections::BTreeSet<_>>()
                        == expected_parents
                })
    });
    let mut staged = frame.clone();
    let target = if let Some(existing) = existing {
        EventBinding {
            document_id: existing.id.clone(),
            event_id: existing.head_event().into(),
            revision: existing.head_revision().into(),
        }
    } else {
        if docs.len() >= 1000 {
            return Err("The workspace document limit has been reached.".into());
        }
        let prefix = format!("assembly-edited-intent-{}", &fresh_context[..16]);
        let id = (1..=1001)
            .map(|n| format!("{prefix}-{n}"))
            .find(|id| !docs.iter().any(|d| d.id == *id))
            .ok_or("No assembly document identity is available.")?;
        let derived=DocumentHistory::derive_from_events(&parents,id,fresh_target,"Declared assembly design intent from pinned source regions and explicit junction edits and literal internal edits; declared bases are design intent, not observed DNA or verified physical preparation").map_err(|e|e.to_string())?;
        let target = EventBinding {
            document_id: derived.id.clone(),
            event_id: derived.head_event().into(),
            revision: derived.head_revision().into(),
        };
        derived.contribute(&mut staged).map_err(|e| e.to_string())?;
        target
    };
    staged.active_document = target.document_id.clone();
    if !staged
        .candidates
        .iter()
        .any(|c| c.plan_revision == candidate.plan_revision)
    {
        staged.candidates.push(candidate.clone());
    }
    staged.active_candidate = Some(candidate.plan_revision.clone());
    staged.validate().map_err(|e| e.to_string())?;
    Ok(Staged {
        workspace: staged,
        target,
        sources: sources
            .iter()
            .zip(&candidate.plan.inputs)
            .map(|(event, input)| InputEventBinding {
                input_id: input.id.clone(),
                event: event.clone(),
            })
            .collect(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use collagene_core::{
        Molecule, Topology,
        verification::{DesignPlan, Operation, PlanInput, Source},
    };
    fn fixture() -> (Workspace, Report, Vec<EventBinding>, Candidate) {
        let mut frame = Workspace {
            name: "assembly workspace".into(),
            active_document: String::new(),
            documents: vec![],
            revisions: Default::default(),
            history: vec![],
            candidates: vec![],
            active_candidate: None,
            attachments: Default::default(),
        };
        let doc = DocumentHistory::new(
            "source",
            Molecule::dna("template", "ACGTACGT", Topology::Linear).unwrap(),
        )
        .unwrap();
        let binding = EventBinding {
            document_id: doc.id.clone(),
            event_id: doc.head_event().into(),
            revision: doc.head_revision().into(),
        };
        let inputs = vec![PlanInput {
            id: "source-1".into(),
            revision: binding.revision.clone(),
            molecule: doc.molecule.clone(),
        }];
        doc.contribute(&mut frame).unwrap();
        frame.active_document = doc.id.clone();
        let intent = assembly_intent::build(
            &inputs,
            &assembly_intent::Request {
                material_bindings: Vec::new(),
                name: "intended prefix".into(),
                fragments: vec![collagene_core::assembly_intent::Contribution {
                    role: "backbone".into(),
                    input_id: "source-1".into(),
                    selector: collagene_core::assembly_intent::Selector::ReferenceRange {
                        start: 0,
                        end: 4,
                        reverse: false,
                    },
                }],
                circular: false,
            },
            None,
        )
        .unwrap();
        let candidate = Candidate::new(DesignPlan {
            name: "proposal".into(),
            target_revision: intent.target.revision(),
            inputs,
            steps: vec![Operation::Overlap {
                fragments: vec![Source::Input("source-1".into())],
                overlaps: vec![],
                circular: false,
            }],
        });
        (frame, intent, vec![binding], candidate)
    }
    fn existing_fixture() -> (Workspace, EventBinding, Vec<InputEventBinding>, Candidate) {
        let (frame, intent, bindings, candidate) = fixture();
        let staged = stage(&frame, &intent, &bindings, &candidate).unwrap();
        let mut frame = staged.workspace;
        frame.candidates.clear();
        frame.active_candidate = None;
        (frame, staged.target, staged.sources, candidate)
    }
    #[test]
    fn paired_insert_stages_two_contributions_and_finalizes_all_three_exact_parents() {
        use crate::cloning_methods_model::{MethodsState, Mode, OverlapPreparation, Treatment};
        use collagene_core::{
            catalog_digest, cloning_outputs, ligation_recipe, run_lineage, verification,
        };
        let originals = vec![
            Molecule::dna("Backbone original", "GGAATTCTTACG", Topology::Circular).unwrap(),
            Molecule::dna("Top original", "AATTGCTACG", Topology::Linear).unwrap(),
            Molecule::dna("Bottom original", "AATTCGTAGC", Topology::Linear).unwrap(),
        ];
        let frame = Workspace::from_documents("Three original events", originals, 0).unwrap();
        let mut state = MethodsState::default();
        state.mode = Mode::Ligation;
        let catalog = catalog_digest::builtin_catalog().unwrap();
        let ecori = catalog
            .edition()
            .enzymes
            .iter()
            .find(|enzyme| enzyme.name == "EcoRI")
            .unwrap()
            .id
            .clone();
        state.rows[0].preparation = OverlapPreparation::CatalogDigest;
        state.rows[0].catalog_selection = Some(catalog_digest::select(&catalog, &[ecori]).unwrap());
        state.rows[1].preparation = OverlapPreparation::AnnealedOligos;
        state.rows[1].anneal_bottom_document = Some(2);
        state.rows[1].anneal_options.min_paired_bases = 6;
        state.rows[1].anneal_candidate = Some(0);
        state.rows[1].treatment = Treatment::Phosphorylate;
        let ids = state.all_input_ids().unwrap();
        let bindings = frame
            .documents
            .iter()
            .map(|document| EventBinding {
                document_id: document.id.clone(),
                event_id: document.head_event.clone().unwrap(),
                revision: document.head.clone(),
            })
            .collect::<Vec<_>>();
        let inputs = ids
            .iter()
            .zip(&bindings)
            .map(|(id, binding)| PlanInput {
                id: id.clone(),
                revision: binding.revision.clone(),
                molecule: frame.revisions[&binding.revision].clone(),
            })
            .collect::<Vec<_>>();
        let roles = state.role_ids(2).unwrap();
        let declaration = assembly_intent::Request {
            material_bindings: Vec::new(),
            name: "Independent two-contribution intent".into(),
            circular: true,
            fragments: vec![
                collagene_core::assembly_intent::Contribution {
                    role: roles[0].role.clone(),
                    input_id: ids[0].clone(),
                    selector: collagene_core::assembly_intent::Selector::ReferenceRange {
                        start: 2,
                        end: 14,
                        reverse: false,
                    },
                },
                collagene_core::assembly_intent::Contribution {
                    role: roles[1].role.clone(),
                    input_id: ids[1].clone(),
                    selector: collagene_core::assembly_intent::Selector::ReferenceRange {
                        start: 0,
                        end: 10,
                        reverse: false,
                    },
                },
            ],
        };
        let intent = IntendedReference::Regions(
            assembly_intent::build(&inputs, &declaration, None).unwrap(),
        );
        // Literal independent expectation precedes actual digestion/annealing.
        assert_eq!(intent.target().sequence.as_str(), "AATTCTTACGGGAATTGCTACG");
        assert_eq!(intent.fragments().len(), 2);
        let (_, rows) = crate::cloning_preparations::prepared_rows(&state, 3).unwrap();
        let recipe = ligation_recipe::Recipe {
            name: "Backbone and annealed insert".into(),
            circular: true,
            allow_nicked: false,
            fragments: rows
                .into_iter()
                .zip(&state.rows)
                .map(|(source, row)| ligation_recipe::Fragment {
                    source,
                    treatment: row.treatment.operation(),
                })
                .collect(),
        };
        let compiled = ligation_recipe::compile(
            &inputs,
            intent.target(),
            &recipe,
            collagene_core::assembly_recipe::CleavagePolicy::AssumedComplete,
        )
        .unwrap();
        assert_eq!(compiled.plan.inputs.len(), 3);
        assert_eq!(
            verification::verify(&compiled.plan, intent.target()).status,
            verification::Status::ConsistentUnderAssumptions
        );
        assert!(compiled.plan.steps.iter().any(|step|matches!(step,Operation::AnnealOligos{top:Source::Input(top),bottom:Source::Input(bottom),..} if top==&ids[1]&&bottom==&ids[2])));
        let candidate = Candidate::new(compiled.plan);
        let before = serde_json::to_value(&frame).unwrap();
        let staged = stage_any(&frame, &intent, &bindings, &candidate).unwrap();
        assert_eq!(serde_json::to_value(&frame).unwrap(), before);
        assert_eq!(staged.sources.len(), 3);
        for (index, pinned) in staged.sources.iter().enumerate() {
            assert_eq!(pinned.input_id, ids[index]);
            assert_eq!(pinned.event, bindings[index]);
        }
        let intended_event = staged
            .workspace
            .history
            .iter()
            .find(|event| event.id == staged.target.event_id)
            .unwrap();
        assert_eq!(
            intended_event
                .parents
                .iter()
                .collect::<std::collections::BTreeSet<_>>(),
            bindings.iter().map(|source| &source.event_id).collect()
        );
        let committed = cloning_outputs::finalize(
            &staged.workspace,
            &candidate.plan,
            staged.target.clone(),
            staged.sources.clone(),
            &cloning_outputs::Options {
                name: "Named adaptor clone".into(),
                show_prepared: false,
                save_primers: false,
                allow_incomplete: false,
            },
        )
        .unwrap();
        let recorded =
            run_lineage::load_run(&committed.workspace, &committed.receipt.run.run_id).unwrap();
        assert_eq!(recorded.inputs.len(), 3);
        for original in &bindings {
            assert!(recorded.inputs.iter().any(|input| &input.event == original));
            assert!(
                committed
                    .workspace
                    .history
                    .iter()
                    .any(|event| event.id == original.event_id)
            );
            assert!(
                committed
                    .workspace
                    .revisions
                    .contains_key(&original.revision)
            );
        }
        assert_eq!(
            committed.workspace.revisions[&committed.receipt.construct.revision]
                .sequence
                .as_str(),
            "AATTCTTACGGGAATTGCTACG"
        );
        assert_eq!(
            committed.workspace.revisions[&committed.receipt.construct.revision].name,
            "Named adaptor clone"
        );
    }
    fn paired_edit_fixture() -> (
        Workspace,
        crate::cloning_methods_model::MethodsState,
        Vec<PlanInput>,
        Vec<EventBinding>,
        assembly_intent::Request,
    ) {
        use crate::cloning_methods_model::{MethodsState, Mode, OverlapPreparation, Treatment};
        use collagene_core::catalog_digest;
        let originals = vec![
            Molecule::dna("Backbone original", "GGAATTCTTACG", Topology::Circular).unwrap(),
            Molecule::dna("Top original", "AATTGCTACG", Topology::Linear).unwrap(),
            Molecule::dna("Bottom original", "AATTCGTAGC", Topology::Linear).unwrap(),
        ];
        let frame = Workspace::from_documents("Three original events", originals, 0).unwrap();
        let mut state = MethodsState::default();
        state.mode = Mode::Ligation;
        let catalog = catalog_digest::builtin_catalog().unwrap();
        let ecori = catalog
            .edition()
            .enzymes
            .iter()
            .find(|enzyme| enzyme.name == "EcoRI")
            .unwrap()
            .id
            .clone();
        state.rows[0].preparation = OverlapPreparation::CatalogDigest;
        state.rows[0].catalog_selection = Some(catalog_digest::select(&catalog, &[ecori]).unwrap());
        state.rows[1].preparation = OverlapPreparation::AnnealedOligos;
        state.rows[1].anneal_bottom_document = Some(2);
        state.rows[1].anneal_options.min_paired_bases = 6;
        state.rows[1].anneal_candidate = Some(0);
        state.rows[1].treatment = Treatment::Phosphorylate;
        let ids = state.all_input_ids().unwrap();
        let bindings = frame
            .documents
            .iter()
            .map(|document| EventBinding {
                document_id: document.id.clone(),
                event_id: document.head_event.clone().unwrap(),
                revision: document.head.clone(),
            })
            .collect::<Vec<_>>();
        let inputs = ids
            .iter()
            .zip(&bindings)
            .map(|(id, binding)| PlanInput {
                id: id.clone(),
                revision: binding.revision.clone(),
                molecule: frame.revisions[&binding.revision].clone(),
            })
            .collect::<Vec<_>>();
        let roles = state.role_ids(2).unwrap();
        let declaration = assembly_intent::Request {
            material_bindings: Vec::new(),
            name: "Independent two-contribution intent".into(),
            circular: true,
            fragments: vec![
                collagene_core::assembly_intent::Contribution {
                    role: roles[0].role.clone(),
                    input_id: ids[0].clone(),
                    selector: collagene_core::assembly_intent::Selector::ReferenceRange {
                        start: 2,
                        end: 14,
                        reverse: false,
                    },
                },
                collagene_core::assembly_intent::Contribution {
                    role: roles[1].role.clone(),
                    input_id: ids[1].clone(),
                    selector: collagene_core::assembly_intent::Selector::ReferenceRange {
                        start: 0,
                        end: 10,
                        reverse: false,
                    },
                },
            ],
        };
        (frame, state, inputs, bindings, declaration)
    }
    #[test]
    fn literal_ligation_intent_is_declared_before_annealing_and_wrong_literal_cannot_finalize() {
        use collagene_core::{
            assembly_junction_edit::{self, Edit},
            cloning_outputs, ligation_recipe, run_lineage, verification,
        };
        let (frame, state, inputs, bindings, mut declaration) = paired_edit_fixture();
        // Selected source core omits the four top overhang bases. They are an
        // explicit independent requirement, physically supplied by the actual
        // original oligo rather than copied from any reaction output.
        if let collagene_core::assembly_intent::Selector::ReferenceRange { start, .. } =
            &mut declaration.fragments[1].selector
        {
            *start = 4;
        } else {
            panic!("fixture requires reference selection")
        }
        let expected = "AATTCTTACGGGAATTGCTACG";
        let build = |literal: &str| {
            IntendedReference::Junctions(
                assembly_junction_edit::build_intent(
                    &inputs,
                    &declaration,
                    &[Edit {
                        junction: 0,
                        remove_left: 0,
                        remove_right: 0,
                        inserted: literal.into(),
                    }],
                    None,
                )
                .unwrap(),
            )
        };
        let intent = build("AATT");
        assert_eq!(intent.target().sequence.as_str(), expected);
        assert_eq!(intent.fragments().len(), 2);
        let (_, rows) = crate::cloning_preparations::prepared_rows(&state, 3).unwrap();
        let recipe = ligation_recipe::Recipe {
            name: "Literal intent, actual paired oligos".into(),
            circular: true,
            allow_nicked: false,
            fragments: rows
                .into_iter()
                .zip(&state.rows)
                .map(|(source, row)| ligation_recipe::Fragment {
                    source,
                    treatment: row.treatment.operation(),
                })
                .collect(),
        };
        let candidate = Candidate::new(
            ligation_recipe::compile(
                &inputs,
                intent.target(),
                &recipe,
                collagene_core::assembly_recipe::CleavagePolicy::AssumedComplete,
            )
            .unwrap()
            .plan,
        );
        assert_eq!(
            verification::verify(&candidate.plan, intent.target()).status,
            verification::Status::ConsistentUnderAssumptions
        );
        let before = serde_json::to_value(&frame).unwrap();
        let staged = stage_any(&frame, &intent, &bindings, &candidate).unwrap();
        assert_eq!(staged.sources.len(), 3);
        let parents = &staged
            .workspace
            .history
            .iter()
            .find(|e| e.id == staged.target.event_id)
            .unwrap()
            .parents;
        assert_eq!(
            parents.iter().collect::<std::collections::BTreeSet<_>>(),
            bindings.iter().map(|b| &b.event_id).collect()
        );
        let options = cloning_outputs::Options {
            name: "Literal-intent clone".into(),
            show_prepared: false,
            save_primers: false,
            allow_incomplete: false,
        };
        let commit = cloning_outputs::finalize(
            &staged.workspace,
            &candidate.plan,
            staged.target,
            staged.sources,
            &options,
        )
        .unwrap();
        let run = run_lineage::load_run(&commit.workspace, &commit.receipt.run.run_id).unwrap();
        assert_eq!(run.inputs.len(), 3);
        for binding in &bindings {
            assert!(run.inputs.iter().any(|i| &i.event == binding));
        }
        assert_eq!(
            commit.workspace.revisions[&commit.receipt.construct.revision]
                .sequence
                .as_str(),
            expected
        );
        let wrong = build("AATC");
        assert_eq!(wrong.target().sequence.as_str(), "AATTCTTACGGGAATCGCTACG");
        let candidate = Candidate::new(
            ligation_recipe::compile(
                &inputs,
                wrong.target(),
                &recipe,
                collagene_core::assembly_recipe::CleavagePolicy::AssumedComplete,
            )
            .unwrap()
            .plan,
        );
        assert_eq!(
            verification::verify(&candidate.plan, wrong.target()).status,
            verification::Status::Invalid
        );
        let staged = stage_any(&frame, &wrong, &bindings, &candidate).unwrap();
        let unchanged = serde_json::to_value(&staged.workspace).unwrap();
        assert!(
            cloning_outputs::finalize(
                &staged.workspace,
                &candidate.plan,
                staged.target,
                staged.sources,
                &options
            )
            .is_err()
        );
        assert_eq!(serde_json::to_value(&staged.workspace).unwrap(), unchanged);
        assert_eq!(serde_json::to_value(&frame).unwrap(), before);
    }
    #[test]
    fn swapping_actual_ordered_oligos_changes_independent_intent_and_requires_fresh_selection() {
        use crate::{
            assembly_draft_model::{CapturedSource, Draft},
            source_identity::SourceInfo,
        };
        use collagene_core::{cloning_outputs, ligation_recipe, run_lineage, verification};
        use std::sync::Arc;
        let (frame, state, inputs, bindings, _) = paired_edit_fixture();
        let captures = inputs
            .iter()
            .zip(&bindings)
            .map(|(input, binding)| {
                Arc::new(CapturedSource {
                    source: SourceInfo {
                        binding: binding.clone(),
                        name: input.molecule.name.clone(),
                        length: input.molecule.sequence.len(),
                        alphabet: input.molecule.sequence.alphabet(),
                        topology: input.molecule.topology,
                    },
                    molecule: input.molecule.clone(),
                })
            })
            .collect::<Vec<_>>();
        let mut draft = Draft::default();
        draft.methods = state;
        draft.slots[0].captured = Some(captures[0].clone());
        {
            let collagene_core::assembly_intent::Selector::ReferenceRange { start: value, .. } =
                &mut draft.slots[0].selector
            else {
                panic!("fixture requires reference selector")
            };
            *value = 2;
        }
        {
            let collagene_core::assembly_intent::Selector::ReferenceRange { end: value, .. } =
                &mut draft.slots[0].selector
            else {
                panic!("fixture requires reference selector")
            };
            *value = 14;
        }
        draft.slots[1].captured = Some(captures[1].clone());
        {
            let collagene_core::assembly_intent::Selector::ReferenceRange { start: value, .. } =
                &mut draft.slots[1].selector
            else {
                panic!("fixture requires reference selector")
            };
            *value = 0;
        }
        {
            let collagene_core::assembly_intent::Selector::ReferenceRange { end: value, .. } =
                &mut draft.slots[1].selector
            else {
                panic!("fixture requires reference selector")
            };
            *value = 10;
        }
        draft.slots[1].anneal_bottom = Some(captures[2].clone());
        draft.sync_rows();
        draft.methods.rows[1].anneal_candidate = Some(0);
        draft.validate().unwrap();
        let original = crate::assembly_draft_model::encode_raw(&draft).unwrap();
        let mut swapped = draft.swap_oligo_roles(1).unwrap();
        assert_eq!(
            crate::assembly_draft_model::encode_raw(&draft).unwrap(),
            original
        );
        assert_eq!(
            swapped.slots[1].captured.as_ref().unwrap().source.binding,
            bindings[2]
        );
        assert_eq!(
            swapped.slots[1]
                .anneal_bottom
                .as_ref()
                .unwrap()
                .source
                .binding,
            bindings[1]
        );
        assert_eq!(swapped.methods.rows[1].anneal_candidate, None);
        assert!(crate::cloning_preparations::prepared_rows(&swapped.methods, 3).is_err());
        let options = swapped.methods.rows[1].anneal_options.clone();
        let top = &swapped.slots[1].captured.as_ref().unwrap().molecule;
        let bottom = &swapped.slots[1].anneal_bottom.as_ref().unwrap().molecule;
        let fresh = verification::simulate_oligo_annealing(top, bottom, &options).unwrap();
        assert!(fresh.annealing.complete && fresh.annealing.unique);
        assert_eq!(fresh.annealing.candidates.len(), 1);
        assert_eq!(
            fresh.annealing.candidates[0].duplex.top.as_str(),
            "AATTCGTAGC"
        );
        swapped.methods.rows[1].anneal_candidate = Some(0);
        swapped.validate().unwrap();
        let ids = swapped.methods.all_input_ids().unwrap();
        let ordered_bindings = vec![
            bindings[0].clone(),
            bindings[2].clone(),
            bindings[1].clone(),
        ];
        let ordered = ids
            .iter()
            .zip(&ordered_bindings)
            .map(|(id, b)| PlanInput {
                id: id.clone(),
                revision: b.revision.clone(),
                molecule: frame.revisions[&b.revision].clone(),
            })
            .collect::<Vec<_>>();
        let roles = swapped.methods.role_ids(2).unwrap();
        let declaration = assembly_intent::Request {
            material_bindings: Vec::new(),
            name: "Explicit opposite oligo contribution".into(),
            circular: true,
            fragments: vec![
                collagene_core::assembly_intent::Contribution {
                    role: roles[0].role.clone(),
                    input_id: ids[0].clone(),
                    selector: collagene_core::assembly_intent::Selector::ReferenceRange {
                        start: 2,
                        end: 14,
                        reverse: false,
                    },
                },
                collagene_core::assembly_intent::Contribution {
                    role: roles[1].role.clone(),
                    input_id: ids[1].clone(),
                    selector: collagene_core::assembly_intent::Selector::ReferenceRange {
                        start: 0,
                        end: 10,
                        reverse: false,
                    },
                },
            ],
        };
        let intent = IntendedReference::Regions(
            assembly_intent::build(&ordered, &declaration, None).unwrap(),
        );
        assert_eq!(intent.target().sequence.as_str(), "AATTCTTACGGGAATTCGTAGC");
        let (_, rows) = crate::cloning_preparations::prepared_rows(&swapped.methods, 3).unwrap();
        let recipe = ligation_recipe::Recipe {
            name: "Opposite ordered oligos".into(),
            circular: true,
            allow_nicked: false,
            fragments: rows
                .into_iter()
                .zip(&swapped.methods.rows)
                .map(|(source, row)| ligation_recipe::Fragment {
                    source,
                    treatment: row.treatment.operation(),
                })
                .collect(),
        };
        let candidate = Candidate::new(
            ligation_recipe::compile(
                &ordered,
                intent.target(),
                &recipe,
                collagene_core::assembly_recipe::CleavagePolicy::AssumedComplete,
            )
            .unwrap()
            .plan,
        );
        assert_eq!(
            verification::verify(&candidate.plan, intent.target()).status,
            verification::Status::ConsistentUnderAssumptions
        );
        let staged = stage_any(&frame, &intent, &ordered_bindings, &candidate).unwrap();
        let commit = cloning_outputs::finalize(
            &staged.workspace,
            &candidate.plan,
            staged.target,
            staged.sources,
            &cloning_outputs::Options {
                name: "Opposite insert clone".into(),
                show_prepared: false,
                save_primers: false,
                allow_incomplete: false,
            },
        )
        .unwrap();
        assert_eq!(
            commit.workspace.revisions[&commit.receipt.construct.revision]
                .sequence
                .as_str(),
            "AATTCTTACGGGAATTCGTAGC"
        );
        let run = run_lineage::load_run(&commit.workspace, &commit.receipt.run.run_id).unwrap();
        assert_eq!(run.inputs.len(), 3);
        for binding in &bindings {
            assert!(run.inputs.iter().any(|input| &input.event == binding));
        }
    }
    #[test]
    fn existing_target_outputs_stage_only_the_candidate_without_inventing_intent() {
        let (frame, target, sources, candidate) = existing_fixture();
        let before = serde_json::to_value(&frame).unwrap();
        let staged = stage_existing(&frame, &candidate, &target, &sources).unwrap();
        assert_eq!(staged.target, target);
        assert_eq!(
            serde_json::to_value(&staged.sources).unwrap(),
            serde_json::to_value(&sources).unwrap()
        );
        assert_eq!(staged.workspace.documents.len(), 2);
        assert_eq!(staged.workspace.active_document, frame.active_document);
        assert_eq!(
            staged.workspace.active_candidate.as_ref(),
            Some(&candidate.plan_revision)
        );
        assert_eq!(staged.workspace.candidates.len(), 1);
        assert!(staged.workspace.attachments.is_empty());
        let again = stage_existing(&staged.workspace, &candidate, &target, &sources).unwrap();
        assert_eq!(
            serde_json::to_value(&again.workspace).unwrap(),
            serde_json::to_value(&staged.workspace).unwrap()
        );
        let mut without_candidate = staged.workspace;
        without_candidate.candidates.clear();
        without_candidate.active_candidate = None;
        assert_eq!(serde_json::to_value(&without_candidate).unwrap(), before);
        assert_eq!(serde_json::to_value(&frame).unwrap(), before);
    }
    #[test]
    fn existing_target_outputs_use_original_source_events_without_rewinding_heads() {
        let (mut frame, target, sources, candidate) = existing_fixture();
        let mut docs = DocumentHistory::restore_all(&frame).unwrap();
        let source = docs
            .iter_mut()
            .find(|doc| doc.id == sources[0].event.document_id)
            .unwrap();
        source
            .commit(
                source.molecule.replace(0, 1, "T").unwrap(),
                "Changed original after inverse-PCR design",
            )
            .unwrap();
        let advanced = source.head_event().to_owned();
        frame.documents.clear();
        frame.history.clear();
        frame.revisions.clear();
        for doc in &docs {
            doc.contribute(&mut frame).unwrap();
        }
        let staged = stage_existing(&frame, &candidate, &target, &sources).unwrap();
        assert_eq!(staged.sources[0].event, sources[0].event);
        assert_eq!(
            staged
                .workspace
                .documents
                .iter()
                .find(|doc| doc.id == sources[0].event.document_id)
                .unwrap()
                .head_event
                .as_deref(),
            Some(advanced.as_str())
        );
        assert_eq!(
            staged.workspace.revisions[&sources[0].event.revision]
                .sequence
                .as_str(),
            "ACGTACGT"
        );
    }
    #[test]
    fn existing_target_outputs_reject_alias_substitution_and_keep_the_exact_retained_target() {
        let (mut frame, target, sources, candidate) = existing_fixture();
        let alias = DocumentHistory::new(
            "equal-source-alias",
            candidate.plan.inputs[0].molecule.clone(),
        )
        .unwrap();
        alias.contribute(&mut frame).unwrap();
        let before = serde_json::to_value(&frame).unwrap();
        let mut swapped = sources.clone();
        // Equal content in this other document cannot supply the original event.
        swapped[0].event.document_id = alias.id;
        assert!(stage_existing(&frame, &candidate, &target, &swapped).is_err());
        assert_eq!(serde_json::to_value(&frame).unwrap(), before);
        let mut docs = DocumentHistory::restore_all(&frame).unwrap();
        let intended = docs
            .iter_mut()
            .find(|doc| doc.id == target.document_id)
            .unwrap();
        let original = intended.molecule.clone();
        intended
            .commit(
                original.replace(0, 1, "T").unwrap(),
                "Temporary target edit",
            )
            .unwrap();
        intended
            .commit(original, "Reverted DNA in a distinct target event")
            .unwrap();
        assert_eq!(intended.head_revision(), target.revision);
        let next = EventBinding {
            document_id: intended.id.clone(),
            event_id: intended.head_event().into(),
            revision: intended.head_revision().into(),
        };
        assert_ne!(next.event_id, target.event_id);
        frame.documents.clear();
        frame.history.clear();
        frame.revisions.clear();
        for doc in &docs {
            doc.contribute(&mut frame).unwrap();
        }
        let before = serde_json::to_value(&frame).unwrap();
        let retained = stage_existing(&frame, &candidate, &target, &sources).unwrap();
        assert_eq!(retained.target, target);
        assert_eq!(serde_json::to_value(&frame).unwrap(), before);
        assert!(stage_existing(&frame, &candidate, &next, &sources).is_ok());
        let mut foreign = target.clone();
        foreign.document_id = sources[0].event.document_id.clone();
        assert!(stage_existing(&frame, &candidate, &foreign, &sources).is_err());
    }
    #[test]
    fn existing_target_outputs_require_ordered_unique_input_ids_and_a_retained_target() {
        let (mut frame, target, mut sources, candidate) = existing_fixture();
        let mut plan = candidate.plan;
        let mut second = plan.inputs[0].clone();
        second.id = "second-source-role".into();
        plan.inputs.push(second);
        sources.push(InputEventBinding {
            input_id: "second-source-role".into(),
            event: sources[0].event.clone(),
        });
        let candidate = Candidate::new(plan);
        assert!(stage_existing(&frame, &candidate, &target, &sources).is_ok());
        let mut incorrect = sources.clone();
        incorrect.swap(0, 1);
        assert!(stage_existing(&frame, &candidate, &target, &incorrect).is_err());
        incorrect = sources.clone();
        incorrect[1].input_id = incorrect[0].input_id.clone();
        assert!(stage_existing(&frame, &candidate, &target, &incorrect).is_err());
        assert!(stage_existing(&frame, &candidate, &target, &sources[..1]).is_err());
        frame.active_document = sources[0].event.document_id.clone();
        let staged = stage_existing(&frame, &candidate, &target, &sources).unwrap();
        assert_eq!(staged.target, target);
        assert_eq!(staged.workspace.active_document, frame.active_document);
    }
    #[test]
    fn target_and_candidate_install_atomically_without_a_positive_check() {
        let (frame, intent, bindings, candidate) = fixture();
        let staged = stage(&frame, &intent, &bindings, &candidate).unwrap();
        assert_eq!(frame.documents.len(), 1);
        assert!(frame.candidates.is_empty());
        assert_eq!(staged.workspace.documents.len(), 2);
        assert_eq!(
            staged.workspace.revisions[&staged.target.revision]
                .sequence
                .as_str(),
            "ACGT"
        );
        assert_eq!(staged.sources[0].event, bindings[0]);
        assert!(staged.workspace.attachments.is_empty());
        assert_eq!(
            staged
                .workspace
                .history
                .iter()
                .find(|e| e.id == staged.target.event_id)
                .unwrap()
                .parents,
            vec![bindings[0].event_id.clone()]
        );
    }
    #[test]
    fn a_pinned_historical_source_is_resolved_without_rewinding_its_head() {
        let (mut frame, intent, bindings, candidate) = fixture();
        let mut docs = DocumentHistory::restore_all(&frame).unwrap();
        let edited = docs[0].molecule.replace(0, 1, "T").unwrap();
        docs[0]
            .commit(edited, "Edited after selecting assembly region")
            .unwrap();
        let advanced = docs[0].head_revision().to_owned();
        frame.documents.clear();
        frame.revisions.clear();
        frame.history.clear();
        docs[0].contribute(&mut frame).unwrap();
        let inputs = vec![InputEventBinding {
            input_id: "source-1".into(),
            event: bindings[0].clone(),
        }];
        assert_eq!(
            retained_inputs(&docs, &inputs, &candidate).unwrap()[0]
                .sequence
                .as_str(),
            "ACGTACGT"
        );
        let staged = stage(&frame, &intent, &bindings, &candidate).unwrap();
        assert_eq!(staged.workspace.documents[0].head, advanced);
        assert_eq!(staged.sources[0].event.event_id, bindings[0].event_id);
    }
    #[test]
    fn repeated_review_reuses_target_only_for_the_same_source_events() {
        let (frame, intent, bindings, candidate) = fixture();
        let first = stage(&frame, &intent, &bindings, &candidate).unwrap();
        let second = stage(&first.workspace, &intent, &bindings, &candidate).unwrap();
        assert_eq!(first.target, second.target);
        assert_eq!(second.workspace.documents.len(), 2);
        let mut aliased = second.workspace;
        let original = DocumentHistory::restore_all(&aliased).unwrap()[0]
            .molecule
            .clone();
        let alias = DocumentHistory::new("distinct-source", original).unwrap();
        let source = EventBinding {
            document_id: alias.id.clone(),
            event_id: alias.head_event().into(),
            revision: alias.head_revision().into(),
        };
        alias.contribute(&mut aliased).unwrap();
        let distinct = stage(&aliased, &intent, &[source.clone()], &candidate).unwrap();
        assert_ne!(first.target, distinct.target);
        assert_eq!(
            distinct
                .workspace
                .history
                .iter()
                .find(|event| event.id == distinct.target.event_id)
                .unwrap()
                .parents,
            vec![source.event_id]
        );
    }
    #[test]
    fn stale_binding_or_forged_intended_dna_never_changes_workspace() {
        let (frame, mut intent, mut bindings, candidate) = fixture();
        bindings[0].revision = "0".repeat(64);
        assert!(stage(&frame, &intent, &bindings, &candidate).is_err());
        let (_, _, bindings, _) = fixture();
        intent.target.sequence = collagene_core::Sequence::dna("TTTT").unwrap();
        assert!(stage(&frame, &intent, &bindings, &candidate).is_err());
        assert_eq!(frame.documents.len(), 1);
    }
    fn edited_fixture() -> (Workspace, IntendedReference, Vec<EventBinding>, Candidate) {
        let (frame, plain, bindings, candidate) = fixture();
        let mut base = plain.request;
        base.name = "Explicit inserted and trimmed junction".into();
        base.fragments = vec![
            collagene_core::assembly_intent::Contribution {
                role: "left".into(),
                input_id: "source-1".into(),
                selector: collagene_core::assembly_intent::Selector::ReferenceRange {
                    start: 0,
                    end: 4,
                    reverse: false,
                },
            },
            collagene_core::assembly_intent::Contribution {
                role: "right".into(),
                input_id: "source-1".into(),
                selector: collagene_core::assembly_intent::Selector::ReferenceRange {
                    start: 4,
                    end: 8,
                    reverse: false,
                },
            },
        ];
        let edits = vec![collagene_core::assembly_junction_edit::Edit {
            junction: 0,
            remove_left: 1,
            remove_right: 1,
            inserted: "GG".into(),
        }];
        let edited = collagene_core::assembly_junction_edit::build_intent(
            &candidate.plan.inputs,
            &base,
            &edits,
            None,
        )
        .unwrap();
        assert_eq!(edited.target.sequence.as_str(), "ACGGGCGT");
        let mut plan = candidate.plan;
        plan.target_revision = edited.target.revision();
        (
            frame,
            IntendedReference::Junctions(edited),
            bindings,
            Candidate::new(plan),
        )
    }
    #[test]
    fn edited_junction_target_and_candidate_stage_with_exact_ancestry_and_retry_reuse() {
        let (frame, intent, bindings, candidate) = edited_fixture();
        let before = serde_json::to_value(&frame).unwrap();
        assert_eq!(intent.target().sequence.as_str(), "ACGGGCGT");
        assert_eq!(intent.fragments().len(), 2);
        let staged = stage_any(&frame, &intent, &bindings, &candidate).unwrap();
        assert_eq!(serde_json::to_value(&frame).unwrap(), before);
        assert_eq!(
            staged.workspace.revisions[&staged.target.revision]
                .sequence
                .as_str(),
            "ACGGGCGT"
        );
        assert_eq!(staged.workspace.documents[0].head, frame.documents[0].head);
        assert_eq!(
            staged
                .workspace
                .history
                .iter()
                .find(|e| e.id == staged.target.event_id)
                .unwrap()
                .parents,
            vec![bindings[0].event_id.clone()]
        );
        assert_eq!(
            staged.workspace.candidates[0].plan_revision,
            candidate.plan_revision
        );
        let again = stage_any(&staged.workspace, &intent, &bindings, &candidate).unwrap();
        assert_eq!(again.target, staged.target);
        assert_eq!(again.workspace.documents.len(), 2);
        // The original region path still dispatches to its unchanged staging contract.
        let (old, plain, sources, c) = fixture();
        let through = stage_any(
            &old,
            &IntendedReference::Regions(plain.clone()),
            &sources,
            &c,
        )
        .unwrap();
        let direct = stage(&old, &plain, &sources, &c).unwrap();
        assert_eq!(through.target, direct.target);
    }
    #[test]
    fn edited_equal_dna_aliases_do_not_reuse_a_target_from_other_source_events() {
        let (frame, intent, bindings, candidate) = edited_fixture();
        let first = stage_any(&frame, &intent, &bindings, &candidate).unwrap();
        let mut next = first.workspace;
        let original = next.revisions[&bindings[0].revision].clone();
        let alias = DocumentHistory::new("edited-intent-distinct-source", original).unwrap();
        let binding = EventBinding {
            document_id: alias.id.clone(),
            event_id: alias.head_event().into(),
            revision: alias.head_revision().into(),
        };
        alias.contribute(&mut next).unwrap();
        let other = stage_any(&next, &intent, &[binding.clone()], &candidate).unwrap();
        assert_ne!(first.target, other.target);
        assert_eq!(
            other
                .workspace
                .history
                .iter()
                .find(|e| e.id == other.target.event_id)
                .unwrap()
                .parents,
            vec![binding.event_id]
        );
    }
    #[test]
    fn edited_staging_rejects_forged_target_or_context_and_stale_source_events() {
        let (frame, intent, bindings, candidate) = edited_fixture();
        let before = serde_json::to_value(&frame).unwrap();
        let mut stale = bindings.clone();
        stale[0].event_id = "0".repeat(64);
        assert!(stage_any(&frame, &intent, &stale, &candidate).is_err());
        let IntendedReference::Junctions(original) = intent else {
            unreachable!()
        };
        let mut forged = original.clone();
        forged.target.sequence = collagene_core::Sequence::dna("TTTTTTTT").unwrap();
        assert!(
            stage_any(
                &frame,
                &IntendedReference::Junctions(forged),
                &bindings,
                &candidate
            )
            .is_err()
        );
        let mut forged = original;
        forged.context_revision = "0".repeat(64);
        assert!(
            stage_any(
                &frame,
                &IntendedReference::Junctions(forged),
                &bindings,
                &candidate
            )
            .is_err()
        );
        assert_eq!(serde_json::to_value(&frame).unwrap(), before);
    }
    #[test]
    fn substituted_intent_stages_actual_mismatch_pcr_and_rejects_stale_sources_or_forged_target() {
        use collagene_core::{
            analysis::Primer, assembly_substitution_intent as substitutions, verification,
        };
        let source =
            Molecule::dna("Original parent", "ACGACCCCGGGGTAAC", Topology::Linear).unwrap();
        let frame =
            Workspace::from_documents("Explicit mutation", vec![source.clone()], 0).unwrap();
        let doc = &frame.documents[0];
        let event = EventBinding {
            document_id: doc.id.clone(),
            event_id: doc.head_event.clone().unwrap(),
            revision: doc.head.clone(),
        };
        let inputs = vec![PlanInput {
            id: "original".into(),
            revision: source.revision(),
            molecule: source,
        }];
        let request = substitutions::Request {
            base: assembly_intent::Request {
                material_bindings: Vec::new(),
                name: "Independent internal substitution".into(),
                circular: false,
                fragments: vec![collagene_core::assembly_intent::Contribution {
                    role: "backbone".into(),
                    input_id: "original".into(),
                    selector: collagene_core::assembly_intent::Selector::ReferenceRange {
                        start: 0,
                        end: 16,
                        reverse: false,
                    },
                }],
            },
            junction_edits: vec![],
            substitutions: vec![substitutions::Substitution {
                role: "backbone".into(),
                offset: 2,
                expected: "G".into(),
                replacement: "T".into(),
            }],
        };
        let intent = substitutions::build(&inputs, &request, None).unwrap();
        assert_eq!(intent.target.sequence.as_str(), "ACTACCCCGGGGTAAC");
        let candidate = Candidate::new(DesignPlan {
            name: "Explicit mismatch PCR".into(),
            target_revision: intent.target.revision(),
            inputs,
            steps: vec![Operation::PcrWithMismatches {
                template: Source::Input("original".into()),
                forward: Primer::new("Mutation F", "ACTACCCC", "").unwrap(),
                reverse: Primer::new("Mutation R", "GTTACCCC", "").unwrap(),
                options: Default::default(),
            }],
        });
        let checked = verification::verify(&candidate.plan, &intent.target);
        assert_eq!(
            checked.status,
            verification::Status::ConsistentUnderAssumptions,
            "{:?}",
            checked.findings
        );
        assert_eq!(
            checked.products.last().unwrap().sequence.as_str(),
            "ACTACCCCGGGGTAAC"
        );
        let before = serde_json::to_value(&frame).unwrap();
        let staged = stage_any(
            &frame,
            &IntendedReference::Substitutions(intent.clone()),
            std::slice::from_ref(&event),
            &candidate,
        )
        .unwrap();
        assert_eq!(staged.sources[0].event, event);
        let target = DocumentHistory::restore_all(&staged.workspace)
            .unwrap()
            .into_iter()
            .find(|doc| doc.id == staged.target.document_id)
            .unwrap();
        assert_eq!(
            target
                .summaries()
                .iter()
                .find(|event| event.current)
                .unwrap()
                .parents,
            vec![event.event_id.clone()]
        );
        let mut changed = intent.clone();
        changed
            .target
            .metadata
            .insert("tampered".into(), "yes".into());
        assert!(
            stage_substituted(&frame, &changed, std::slice::from_ref(&event), &candidate).is_err()
        );
        let mut stale = candidate.plan.clone();
        stale.inputs[0].molecule.sequence =
            collagene_core::Sequence::dna("TTTTTTTTTTTTTTTT").unwrap();
        stale.inputs[0].revision = stale.inputs[0].molecule.revision();
        assert!(
            stage_substituted(
                &frame,
                &intent,
                std::slice::from_ref(&event),
                &Candidate::new(stale)
            )
            .is_err()
        );
        let mut changed = intent;
        changed.request.substitutions[0].replacement = "C".into();
        assert!(
            stage_substituted(&frame, &changed, std::slice::from_ref(&event), &candidate).is_err()
        );
        assert_eq!(serde_json::to_value(&frame).unwrap(), before);
    }
    #[test]
    fn contribution_indels_stage_original_events_without_upgrading_wrong_pcr_output() {
        use collagene_core::{
            analysis::Primer, assembly_edit_intent as edits,
            assembly_substitution_intent as substitutions, verification,
        };
        let source =
            Molecule::dna("Original parent", "ACGACCCCGGGGTAAC", Topology::Linear).unwrap();
        let frame =
            Workspace::from_documents("Independent literal indels", vec![source.clone()], 0)
                .unwrap();
        let doc = &frame.documents[0];
        let original_event = EventBinding {
            document_id: doc.id.clone(),
            event_id: doc.head_event.clone().unwrap(),
            revision: doc.head.clone(),
        };
        let inputs = vec![PlanInput {
            id: "original".into(),
            revision: source.revision(),
            molecule: source.clone(),
        }];
        let request = edits::Request {
            base: substitutions::Request {
                base: assembly_intent::Request {
                    material_bindings: Vec::new(),
                    name: "Independent insertion and deletion".into(),
                    circular: false,
                    fragments: vec![collagene_core::assembly_intent::Contribution {
                        role: "backbone".into(),
                        input_id: "original".into(),
                        selector: collagene_core::assembly_intent::Selector::ReferenceRange {
                            start: 0,
                            end: 16,
                            reverse: false,
                        },
                    }],
                },
                junction_edits: vec![],
                substitutions: vec![],
            },
            edits: vec![
                edits::Edit {
                    role: "backbone".into(),
                    start: 4,
                    end: 4,
                    expected: "".into(),
                    replacement: "TT".into(),
                },
                edits::Edit {
                    role: "backbone".into(),
                    start: 8,
                    end: 10,
                    expected: "GG".into(),
                    replacement: "".into(),
                },
            ],
        };
        let intent = edits::build(&inputs, &request, None).unwrap();
        assert_eq!(intent.target.sequence.as_str(), "ACGATTCCCCGGTAAC");
        let candidate = Candidate::new(DesignPlan {
            name: "Unedited PCR cannot satisfy literal indels".into(),
            target_revision: intent.target.revision(),
            inputs,
            steps: vec![Operation::Pcr {
                template_strand: collagene_core::chemistry::Strand::Top,
                template: Source::Input("original".into()),
                forward: Primer::new("Original F", "ACGACCCC", "").unwrap(),
                reverse: Primer::new("Original R", "GTTACCCC", "").unwrap(),
            }],
        });
        let checked = verification::verify(&candidate.plan, &intent.target);
        assert_eq!(
            checked.status,
            verification::Status::Invalid,
            "{:?}",
            checked.findings
        );
        assert_eq!(
            checked.products.last().unwrap().sequence.as_str(),
            source.sequence.as_str()
        );
        let before = serde_json::to_value(&frame).unwrap();
        let original_inputs = serde_json::to_value(&candidate.plan.inputs).unwrap();
        let staged = stage_any(
            &frame,
            &IntendedReference::Edits(intent.clone()),
            std::slice::from_ref(&original_event),
            &candidate,
        )
        .unwrap();
        assert_eq!(staged.sources[0].event, original_event);
        assert_eq!(staged.sources[0].input_id, "original");
        assert_eq!(staged.workspace.documents[0].head, frame.documents[0].head);
        assert_eq!(
            serde_json::to_value(&staged.workspace.candidates[0].plan.inputs).unwrap(),
            original_inputs
        );
        let target_document = DocumentHistory::restore_all(&staged.workspace)
            .unwrap()
            .into_iter()
            .find(|doc| doc.id == staged.target.document_id)
            .unwrap();
        assert_eq!(target_document.molecule, intent.target);
        assert_eq!(
            target_document
                .summaries()
                .iter()
                .find(|e| e.current)
                .unwrap()
                .parents,
            vec![original_event.event_id.clone()]
        );
        let rechecked = verification::verify(
            &staged.workspace.candidates[0].plan,
            &target_document.molecule,
        );
        assert_eq!(rechecked.status, verification::Status::Invalid);
        assert_eq!(
            rechecked.products.last().unwrap().sequence.as_str(),
            source.sequence.as_str()
        );
        assert_eq!(serde_json::to_value(&frame).unwrap(), before);
        let again = stage_contribution_edits(
            &staged.workspace,
            &intent,
            std::slice::from_ref(&original_event),
            &candidate,
        )
        .unwrap();
        assert_eq!(again.target, staged.target);
        assert_eq!(again.workspace.documents.len(), 2);

        let mut forged = intent.clone();
        forged
            .target
            .metadata
            .insert("tampered".into(), "yes".into());
        assert!(
            stage_contribution_edits(
                &frame,
                &forged,
                std::slice::from_ref(&original_event),
                &candidate
            )
            .is_err()
        );
        forged = intent.clone();
        forged.context_revision = "0".repeat(64);
        assert!(
            stage_contribution_edits(
                &frame,
                &forged,
                std::slice::from_ref(&original_event),
                &candidate
            )
            .is_err()
        );
        forged = intent.clone();
        forged.applied[0].output_start += 1;
        assert!(
            stage_contribution_edits(
                &frame,
                &forged,
                std::slice::from_ref(&original_event),
                &candidate
            )
            .is_err()
        );
        let mut stale = original_event.clone();
        stale.event_id = "0".repeat(64);
        assert!(stage_contribution_edits(&frame, &intent, &[stale], &candidate).is_err());
        let mut changed = candidate.plan.clone();
        changed.inputs[0].molecule.sequence =
            collagene_core::Sequence::dna("TTTTTTTTTTTTTTTT").unwrap();
        changed.inputs[0].revision = changed.inputs[0].molecule.revision();
        assert!(
            stage_contribution_edits(
                &frame,
                &intent,
                std::slice::from_ref(&original_event),
                &Candidate::new(changed)
            )
            .is_err()
        );
        assert_eq!(serde_json::to_value(&frame).unwrap(), before);
        assert_eq!(
            serde_json::to_value(&candidate.plan.inputs).unwrap(),
            original_inputs
        );
    }
    #[test]
    fn declared_inverse_pcr_indel_stages_and_finalizes_original_source_ancestry() {
        use collagene_core::{
            analysis, assembly_edit_intent as edits, assembly_substitution_intent as substitutions,
            chemistry::EndTreatment, cloning_outputs, mutagenesis, run_lineage,
        };
        let original = Molecule::dna(
            "Original circular parent",
            "AGTCCGATGCTACGTAGGACTTCAGACCTGTA",
            Topology::Circular,
        )
        .unwrap();
        let frame =
            Workspace::from_documents("Inverse-PCR dock", vec![original.clone()], 0).unwrap();
        let doc = &frame.documents[0];
        let original_event = EventBinding {
            document_id: doc.id.clone(),
            event_id: doc.head_event.clone().unwrap(),
            revision: doc.head.clone(),
        };
        let inputs = vec![PlanInput {
            id: "pinned-original".into(),
            revision: original.revision(),
            molecule: original.clone(),
        }];
        let intent = edits::build(
            &inputs,
            &edits::Request {
                base: substitutions::Request {
                    base: assembly_intent::Request {
                        material_bindings: Vec::new(),
                        name: "Independent four-base insertion".into(),
                        circular: true,
                        fragments: vec![collagene_core::assembly_intent::Contribution {
                            role: "backbone".into(),
                            input_id: "pinned-original".into(),
                            selector: collagene_core::assembly_intent::Selector::ReferenceRange {
                                start: 0,
                                end: 32,
                                reverse: false,
                            },
                        }],
                    },
                    junction_edits: vec![],
                    substitutions: vec![],
                },
                edits: vec![edits::Edit {
                    role: "backbone".into(),
                    start: 10,
                    end: 10,
                    expected: "".into(),
                    replacement: "AACG".into(),
                }],
            },
            None,
        )
        .unwrap();
        assert_eq!(
            intent.target.sequence.as_str(),
            "AGTCCGATGCAACGTACGTAGGACTTCAGACCTGTA"
        );
        let target_before = intent.target.clone();
        let mut preparation = mutagenesis::compile_inverse_pcr(
            &original,
            &intent.target,
            None,
            &mutagenesis::InversePcrOptions {
                edit_start: 10,
                edit_end: 10,
                inserted: "AACG".into(),
                forward_length: 8,
                reverse_length: 8,
                reverse_tail_bases: 2,
            },
        )
        .unwrap();
        assert_eq!(preparation.plan.inputs[0].molecule, original);
        preparation.plan.inputs[0].id = "pinned-original".into();
        let Operation::Pcr {
            template: Source::Input(id),
            ..
        } = &mut preparation.plan.steps[0]
        else {
            panic!("inverse-PCR must start from the actual original source");
        };
        *id = "pinned-original".into();
        let candidate = Candidate::new(preparation.plan);
        assert_eq!(candidate.plan.steps.len(), 3);
        assert!(matches!(
            &candidate.plan.steps[1],
            Operation::TreatEnds {
                template: Source::Step(0),
                treatment: EndTreatment::Phosphorylate,
            }
        ));
        assert!(matches!(&candidate.plan.steps[2], Operation::Ligate {
            fragments, circular: true, allow_nicked: false, ..
        } if fragments.len()==1 && matches!(&fragments[0].source, Source::Step(1)) && !fragments[0].reverse));
        let assessed = collagene_core::verification::verify(&candidate.plan, &intent.target);
        assert_eq!(
            assessed.status,
            collagene_core::verification::Status::ConsistentUnderAssumptions,
            "{:?}",
            assessed.findings
        );
        assert_eq!(assessed.products.len(), 3);
        let actual = assessed.products.last().unwrap();
        assert_eq!(actual.topology, Topology::Circular);
        assert!(
            analysis::rotation_offset(
                actual.sequence.as_str(),
                "AGTCCGATGCAACGTACGTAGGACTTCAGACCTGTA"
            )
            .is_some()
        );
        let original_workspace = serde_json::to_value(&frame).unwrap();
        let original_inputs = serde_json::to_value(&candidate.plan.inputs).unwrap();
        let staged = stage_any(
            &frame,
            &IntendedReference::Edits(intent.clone()),
            std::slice::from_ref(&original_event),
            &candidate,
        )
        .unwrap();
        assert_eq!(staged.sources[0].input_id, "pinned-original");
        assert_eq!(staged.sources[0].event, original_event);
        assert_eq!(
            staged.workspace.revisions[&original_event.revision],
            original
        );
        let target_event = staged
            .workspace
            .history
            .iter()
            .find(|e| e.id == staged.target.event_id)
            .unwrap();
        assert_eq!(target_event.parents, vec![original_event.event_id.clone()]);
        assert_eq!(
            staged.workspace.revisions[&staged.target.revision],
            target_before
        );
        let committed = cloning_outputs::finalize(
            &staged.workspace,
            &candidate.plan,
            staged.target.clone(),
            staged.sources.clone(),
            &cloning_outputs::Options {
                name: "Named inverse-PCR insertion clone".into(),
                show_prepared: false,
                save_primers: true,
                allow_incomplete: false,
            },
        )
        .unwrap();
        assert!(committed.receipt.primer_attachment.is_some());
        let recorded =
            run_lineage::load_run(&committed.workspace, &committed.receipt.run.run_id).unwrap();
        assert_eq!(recorded.inputs.len(), 1);
        assert_eq!(recorded.inputs[0].event, original_event);
        assert_eq!(
            serde_json::to_value(&recorded.plan.inputs).unwrap(),
            original_inputs
        );
        assert_eq!(recorded.target, staged.target);
        let clone = &committed.workspace.revisions[&committed.receipt.construct.revision];
        assert_eq!(clone.name, "Named inverse-PCR insertion clone");
        assert_eq!(clone.topology, Topology::Circular);
        assert!(
            analysis::rotation_offset(clone.sequence.as_str(), target_before.sequence.as_str())
                .is_some()
        );
        let replay =
            run_lineage::replay_run(&committed.workspace, &committed.receipt.run.run_id).unwrap();
        assert!(replay.products_match && replay.status_matches && replay.engine_matches);
        assert_eq!(
            replay.current_report.status,
            collagene_core::verification::Status::ConsistentUnderAssumptions
        );
        assert_eq!(
            committed.workspace.revisions[&original_event.revision],
            original
        );
        assert_eq!(intent.target, target_before);
        assert_eq!(serde_json::to_value(&frame).unwrap(), original_workspace);
        assert_eq!(
            serde_json::to_value(&candidate.plan.inputs).unwrap(),
            original_inputs
        );
    }
}
