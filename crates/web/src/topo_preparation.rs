//! Explicit TA/directional preparation form over the shared native builder.
use collagene_core::{Sequence, analysis::Primer, tailing::TailingBase, topo};
use topo::physical::{self as physical, projection};

#[derive(Clone, Debug, PartialEq, Eq, Default)]
pub struct Draft {
    pub mode: collagene_core::topo_preparation::PreparationMode,
    pub vector_top: String,
    pub vector_bottom: String,
    pub insert_template: String,
    pub vector_top_length: usize,
    pub vector_bottom_length: usize,
    pub forward_sequence: String,
    pub forward_annealing: usize,
    pub reverse_sequence: String,
    pub reverse_annealing: usize,
    pub preparation_assumption: String,
    pub activation_confirmed: bool,
    pub add_a_tail: bool,
}
fn primer(name: &str, ordered: &str, annealing: usize) -> Result<Primer, String> {
    if ordered.len() > 1024 {
        return Err("An ordered primer exceeds the input text bound.".into());
    }
    let ordered = Sequence::dna(ordered).map_err(|e| e.to_string())?;
    if ordered.is_empty()
        || ordered.len() > collagene_core::primers::MAX_OLIGO
        || annealing == 0
        || annealing > ordered.len()
    {
        return Err(format!(
            "Choose a nonempty 3′ annealing suffix within the ordered {name} primer (at most {} bases).",
            collagene_core::primers::MAX_OLIGO
        ));
    }
    let boundary = ordered.len() - annealing;
    Primer::new(
        name,
        &ordered.as_str()[boundary..],
        &ordered.as_str()[..boundary],
    )
    .map_err(|e| e.to_string())
}

/// Explicit source roles remain separate from model identities and stock claims.
pub fn preparation_request(
    draft: &Draft,
) -> Result<
    (
        Vec<(String, String)>,
        collagene_core::topo_preparation::PreparationRequest,
    ),
    String,
> {
    use collagene_core::topo_preparation::{
        ExactPcrDeclaration, PreparationMode, PreparationRequest, PreparedVectorDeclaration,
    };
    for id in [
        &draft.vector_top,
        &draft.vector_bottom,
        &draft.insert_template,
    ] {
        if id.trim().is_empty() || id.len() > 1024 || id.chars().any(char::is_control) {
            return Err("Choose all three source documents.".into());
        }
    }
    if !draft.activation_confirmed {
        return Err(
            "Explicitly declare the prepared vector's5′ phosphate and3′ phosphotyrosyl chemistry."
                .into(),
        );
    }
    if draft.preparation_assumption.trim().is_empty()
        || draft.preparation_assumption.len() > 4096
        || draft
            .preparation_assumption
            .chars()
            .any(|c| c.is_control() && !matches!(c, '\n' | '\r' | '\t'))
    {
        return Err("Describe the actual vector preparation assumption (up to4,096 bytes).".into());
    }
    let top = draft.vector_top_length;
    let bottom = draft.vector_bottom_length;
    let (paired, bottom_start, declaration_id, material_id, tail) = match draft.mode {
        PreparationMode::Ta => {
            if top < 2 || top != bottom || top > 100_000 {
                return Err("TA requires equal vector strand lengths of2–100,000 bases, with one explicit3′ overhang.".into());
            }
            if !draft.add_a_tail {
                return Err(
                    "TA requires an explicit single A addition to both PCR strand3′ ends.".into(),
                );
            }
            (
                top - 1,
                0,
                "user-declared-ta-vector-1",
                "explicit-pcr-a-tail-insert",
                Some(TailingBase::A),
            )
        }
        PreparationMode::DirectionalCacc => {
            if top == 0 || top > 99_996 || top.checked_add(4) != Some(bottom) {
                return Err("Directional preparation requires actual top/bottom strand sources differing by four bases: bottom5′ GGTG, paired from bottom offset4. No strand is synthesized from the other.".into());
            }
            if draft.add_a_tail {
                return Err("Directional CACC uses blunt exact PCR. Explicitly disable A addition; no tail is silently removed.".into());
            }
            (
                top,
                4,
                "user-declared-directional-vector-1",
                "explicit-blunt-pcr-insert",
                None,
            )
        }
    };
    let strands = [
        ("vector_top", "declared-vector-top-3prime", top),
        ("vector_bottom", "declared-vector-bottom-3prime", bottom),
    ]
    .into_iter()
    .map(|(role, activation, length)| topo::PreparedStrand {
        id: role.into(),
        segments: vec![topo::SourceInterval {
            source_role: role.into(),
            start: 0,
            length,
            reverse: false,
        }],
        five_prime: physical::FivePrime::Phosphate,
        three_prime: physical::ThreePrime::Phosphotyrosyl {
            enzyme: physical::Enzyme::VacciniaTypeIBTyr274,
            activation_id: activation.into(),
        },
    })
    .collect();
    let request = PreparationRequest {
        mode: draft.mode,
        vector: PreparedVectorDeclaration {
            declaration_id: declaration_id.into(),
            assumption: draft.preparation_assumption.clone(),
            strands,
            pairing: vec![physical::PairSpan {
                strand_a: 0,
                start_a: 0,
                strand_b: 1,
                start_b: bottom_start,
                length: paired,
            }],
        },
        insert: ExactPcrDeclaration {
            material_id: material_id.into(),
            template_role: "insert_template".into(),
            forward: primer("Forward", &draft.forward_sequence, draft.forward_annealing)?,
            reverse: primer("Reverse", &draft.reverse_sequence, draft.reverse_annealing)?,
            tail,
        },
        engine_limits: physical::Limits::default(),
        projection_options: projection::Options::default(),
    };
    Ok((
        vec![
            ("vector_top".into(), draft.vector_top.clone()),
            ("vector_bottom".into(), draft.vector_bottom.clone()),
            ("insert_template".into(), draft.insert_template.clone()),
        ],
        request,
    ))
}
pub fn compile_request(
    draft: &Draft,
) -> Result<
    (
        Vec<(String, String)>,
        collagene_core::topo_preparation::CompiledPreparation,
    ),
    String,
> {
    let (roles, request) = preparation_request(draft)?;
    let compiled =
        collagene_core::topo_preparation::compile(&request).map_err(|error| error.to_string())?;
    Ok((roles, compiled))
}

/// Propose a literal five-prime tail edit; caller must separately guard the
/// exact original source events and explicitly apply this returned draft.
pub fn cacc_prefix_proposal(draft: &Draft) -> Result<Draft, String> {
    if draft.mode != collagene_core::topo_preparation::PreparationMode::DirectionalCacc {
        return Err("CACC prefix proposal requires directional mode.".into());
    }
    let original = Sequence::dna(&draft.forward_sequence).map_err(|error| error.to_string())?;
    let mut next = draft.clone();
    next.forward_sequence = format!("CACC{}", original.as_str());
    primer("Forward", &next.forward_sequence, next.forward_annealing)?;
    compile_request(&next)?;
    Ok(next)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn declaration_preserves_ordered_oligos_and_never_selects_a_transfer() {
        let draft = Draft {
            mode: collagene_core::topo_preparation::PreparationMode::Ta,
            vector_top: "vt".into(),
            vector_bottom: "vb".into(),
            insert_template: "template".into(),
            vector_top_length: 13,
            vector_bottom_length: 13,
            forward_sequence: "GGACGT".into(),
            forward_annealing: 4,
            reverse_sequence: "TTTGCA".into(),
            reverse_annealing: 4,
            preparation_assumption:
                "Independently supplied prepared vector; explicit chemistry assumption".into(),
            activation_confirmed: true,
            add_a_tail: true,
        };
        let (roles, compiled) = compile_request(&draft).unwrap();
        let request = compiled.request;
        assert_eq!(roles[0], ("vector_top".into(), "vt".into()));
        assert!(request.transfers.is_empty());
        let topo::MaterialRecipe::ExactPcr {
            forward,
            reverse,
            tail,
            ..
        } = &request.recipe.materials[1]
        else {
            panic!()
        };
        assert_eq!(forward.annealing.as_str(), "ACGT");
        assert_eq!(forward.tail, "GG");
        assert_eq!(reverse.ordered_sequence(), "TTTGCA");
        assert_eq!(*tail, Some(TailingBase::A));
        let mut missing = draft.clone();
        missing.activation_confirmed = false;
        assert!(compile_request(&missing).is_err());
        missing = draft;
        missing.forward_annealing = 7;
        assert!(compile_request(&missing).is_err());
    }

    fn directional() -> Draft {
        Draft {
            mode: collagene_core::topo_preparation::PreparationMode::DirectionalCacc,
            vector_top: "vt".into(),
            vector_bottom: "vb".into(),
            insert_template: "original".into(),
            vector_top_length: 14,
            vector_bottom_length: 18,
            forward_sequence: "CATGA".into(),
            forward_annealing: 5,
            reverse_sequence: "GTACG".into(),
            reverse_annealing: 5,
            preparation_assumption:
                "Explicit actual strand preparation; no kit activation is inferred".into(),
            activation_confirmed: true,
            add_a_tail: false,
        }
    }
    #[test]
    fn directional_form_declares_actual_geometry_and_never_reuses_ta_or_adds_a_tail() {
        let draft = directional();
        let original = draft.clone();
        let (roles, compiled) = compile_request(&draft).unwrap();
        assert!(compiled.directional);
        assert!(compiled.request.transfers.is_empty());
        assert_eq!(roles[2], ("insert_template".into(), "original".into()));
        let topo::MaterialRecipe::Prepared {
            pairing, strands, ..
        } = &compiled.request.recipe.materials[0]
        else {
            panic!("missing declared vector");
        };
        assert_eq!(pairing[0].start_b, 4);
        assert_eq!(pairing[0].length, 14);
        assert_eq!(strands[0].segments[0].length, 14);
        assert_eq!(strands[1].segments[0].length, 18);
        let topo::MaterialRecipe::ExactPcr { tail, forward, .. } =
            &compiled.request.recipe.materials[1]
        else {
            panic!("missing actual PCR");
        };
        assert!(tail.is_none());
        assert_eq!(forward.ordered_sequence(), "CATGA");
        assert!(forward.tail.is_empty());
        assert_eq!(draft, original);
        let mut invalid = draft.clone();
        invalid.add_a_tail = true;
        assert!(compile_request(&invalid).is_err());
        invalid = draft.clone();
        invalid.vector_bottom_length = 14;
        assert!(compile_request(&invalid).is_err());
        invalid = draft;
        invalid.activation_confirmed = false;
        assert!(compile_request(&invalid).is_err());
    }
    #[test]
    fn cacc_prefix_is_an_explicit_unchanged_annealing_proposal_not_a_source_edit() {
        let draft = directional();
        let original = draft.clone();
        let proposed = cacc_prefix_proposal(&draft).unwrap();
        assert_eq!(draft, original);
        assert_eq!(proposed.forward_sequence, "CACCCATGA");
        assert_eq!(proposed.forward_annealing, 5);
        assert_eq!(proposed.reverse_sequence, original.reverse_sequence);
        assert_eq!(proposed.insert_template, original.insert_template);
        assert!(!proposed.add_a_tail);
        let (_, request) = preparation_request(&proposed).unwrap();
        assert_eq!(request.insert.forward.tail, "CACC");
        assert_eq!(request.insert.forward.annealing.as_str(), "CATGA");
        let mut too_long = draft;
        too_long.forward_sequence = "A".repeat(collagene_core::primers::MAX_OLIGO);
        too_long.forward_annealing = 20;
        assert!(cacc_prefix_proposal(&too_long).is_err());
    }
}

/// Shared native/standalone dispatch; retained origins always use the native
/// capsule replay endpoint. Original report identities remain unchanged.
pub(crate) fn task_revision(
    inputs: &[topo::Input],
    origin: Option<&collagene_core::run_material::MaterialOrigin>,
    task: &collagene_core::topo_run_preparation::Task,
) -> Result<String, collagene_core::Error> {
    use collagene_core::{
        topo_actions,
        topo_run_preparation::{self, Task},
        topo_worker,
    };
    if let Some(origin) = origin {
        return topo_run_preparation::documentary_revision(inputs, "insert_template", origin, task);
    }
    match task {
        Task::Orientations {
            request,
            vector,
            insert,
            directional,
        } => topo_worker::orientations_revision(inputs, request, *vector, *insert, *directional),
        Task::Project { request, selection } => {
            topo_worker::annotated_projection_revision(inputs, request, selection)
        }
        Task::Check {
            request,
            selection,
            target,
            policy,
        } => topo_worker::annotated_target_revision(inputs, request, selection, target, policy),
        Task::InspectActions { request } => topo_actions::revision_actions(inputs, request),
        Task::ProcessActions { request } => topo_actions::revision(inputs, request),
    }
}
pub(crate) fn batch_revision(
    inputs: &[topo::Input],
    origin: Option<&collagene_core::run_material::MaterialOrigin>,
    request: &topo::Request,
) -> Result<String, collagene_core::Error> {
    match origin {
        Some(origin) => collagene_core::topo_run_preparation::preparation_revision(
            inputs,
            "insert_template",
            origin,
            request,
        ),
        None => collagene_core::topo_worker::simulation_revision(inputs, request),
    }
}
#[derive(Clone)]
pub(crate) struct OriginReview {
    pub incomplete: bool,
    pub findings: Vec<String>,
    pub total: usize,
    pub truncated: bool,
}
#[cfg(target_arch = "wasm32")]
pub(crate) async fn execute(
    worker: std::rc::Rc<crate::worker::CachedWorker>,
    id: u32,
    inputs: Vec<topo::Input>,
    origin: Option<collagene_core::run_material::MaterialOrigin>,
    task: collagene_core::topo_run_preparation::Task,
) -> Result<
    (
        collagene_core::protocol::AnalysisResult,
        Option<OriginReview>,
    ),
    String,
> {
    use collagene_core::{
        protocol::{AnalysisJob, AnalysisResult, Request},
        topo_run_preparation::{self, Reply, Task},
    };
    let (revision, job) = if let Some(origin) = &origin {
        (
            topo_run_preparation::revision(&inputs, "insert_template", origin, &task)
                .map_err(|e| e.to_string())?,
            AnalysisJob::InspectTopoRun {
                inputs: inputs.clone(),
                template_role: "insert_template".into(),
                origin: origin.clone(),
                task: task.clone(),
            },
        )
    } else {
        let revision = task_revision(&inputs, None, &task).map_err(|e| e.to_string())?;
        let job = match task.clone() {
            Task::Orientations {
                request,
                vector,
                insert,
                directional,
            } => AnalysisJob::InspectTopoOrientations {
                inputs: inputs.clone(),
                request,
                vector,
                insert,
                directional,
            },
            Task::Project { request, selection } => AnalysisJob::ProjectAnnotatedTopo {
                inputs: inputs.clone(),
                request,
                selection,
            },
            Task::Check {
                request,
                selection,
                target,
                policy,
            } => AnalysisJob::CheckAnnotatedTopoTarget {
                inputs: inputs.clone(),
                request,
                selection,
                target,
                policy,
            },
            Task::InspectActions { request } => AnalysisJob::InspectTopoActions {
                inputs: inputs.clone(),
                request,
            },
            Task::ProcessActions { request } => AnalysisJob::ProcessTopoActions {
                inputs: inputs.clone(),
                request,
            },
        };
        (revision, job)
    };
    let reply = worker
        .ask(Request {
            id,
            revision: revision.clone(),
            job,
        })
        .await;
    worker.close();
    let reply = reply?;
    if reply.id != id || reply.revision != revision {
        return Err("TOPO response belongs to a different source, origin or declaration.".into());
    }
    match (reply.result, origin) {
        (AnalysisResult::Error(e), _) => Err(e),
        (AnalysisResult::TopoRunInspection(report), Some(origin)) => {
            report
                .validate_for(&inputs, "insert_template", &origin, &task, &revision)
                .map_err(|e| e.to_string())?;
            let review = OriginReview {
                incomplete: report.origin_incomplete,
                findings: report.origin_findings,
                total: report.origin_findings_total,
                truncated: report.origin_findings_truncated,
            };
            Ok((
                match report.reply {
                    Reply::Orientations(r) => AnalysisResult::TopoOrientations(r),
                    Reply::Projected(r) => AnalysisResult::AnnotatedTopoProjection(r),
                    Reply::Checked(r) => AnalysisResult::AnnotatedTopoTarget(r),
                    Reply::ActionsInspected(r) => AnalysisResult::TopoActionsInspected(r),
                    Reply::ActionsProcessed(r) => AnalysisResult::TopoActionsProduct(r),
                },
                Some(review),
            ))
        }
        (_, Some(_)) => Err("Unexpected saved-template inspection response.".into()),
        (value, None) => Ok((value, None)),
    }
}

#[cfg(test)]
mod generic_material_tests {
    use super::*;
    use collagene_core::{
        Molecule, Topology, material_derivation as derivation,
        run_lineage::EventBinding,
        run_material::{self, MaterialOrigin, MaterialSelection},
        verification::{self, Source, TopoSource},
        workspace::Workspace,
    };
    #[test]
    fn native_derived_insert_template_keeps_exact_origin_and_original_vector_declarations() {
        let parent =
            Molecule::dna("Original insert parent", "CACCATGACG", Topology::Linear).unwrap();
        let workspace =
            Workspace::from_documents("Retained insert", vec![parent.clone()], 0).unwrap();
        let doc = &workspace.documents[0];
        let event = EventBinding {
            document_id: doc.id.clone(),
            event_id: doc.head_event.clone().unwrap(),
            revision: doc.head.clone(),
        };
        let forward = Primer::new("F", "CACCA", "").unwrap();
        let reverse = Primer::new("R", "CGTCA", "").unwrap();
        let preview = verification::simulate_pcr_documents(&parent, &forward, &reverse)
            .unwrap()
            .remove(0);
        let installed = derivation::prepare_install(
            &workspace,
            &derivation::Request {
                source: event,
                origin: derivation::SourceOrigin::OriginalReference,
                derivation: derivation::Derivation::ExactPcr {
                    forward,
                    reverse,
                    selection: derivation::PcrSelection {
                        selected_index: 0,
                        expected_preview_revision: preview.molecule.revision(),
                        geometry: preview.geometry,
                    },
                },
                name: "Saved insert PCR".into(),
                save_primers: false,
            },
        )
        .unwrap();
        let product = installed.receipt.product.clone();
        let selected = MaterialSelection::Derived(product.clone());
        let origin =
            run_material::inspect_material(&installed.workspace, &product, &Some(selected.clone()))
                .unwrap()
                .selected
                .unwrap();
        assert!(matches!(&origin, MaterialOrigin::Derived(_)));
        assert_eq!(origin.selection(), selected);
        let template = installed.workspace.revisions[&product.revision].clone();
        assert_eq!(template.sequence.as_str(), "CACCATGACG");
        let original_workspace = serde_json::to_vec(&installed.workspace).unwrap();
        let draft = Draft {
            mode: collagene_core::topo_preparation::PreparationMode::DirectionalCacc,
            vector_top: "vector-top-doc".into(),
            vector_bottom: "vector-bottom-doc".into(),
            insert_template: product.document_id.clone(),
            vector_top_length: 14,
            vector_bottom_length: 18,
            forward_sequence: "CACCA".into(),
            forward_annealing: 5,
            reverse_sequence: "CGTCA".into(),
            reverse_annealing: 5,
            preparation_assumption:
                "Original explicit vector activation; no saved-vector substitution".into(),
            activation_confirmed: true,
            add_a_tail: false,
        };
        let (_, compiled) = compile_request(&draft).unwrap();
        let top = Molecule::dna("Original vector top", "AAGGGGCATCCCTT", Topology::Linear).unwrap();
        let bottom = Molecule::dna(
            "Original vector bottom",
            "GGTGAAGGGATGCCCCTT",
            Topology::Linear,
        )
        .unwrap();
        let roles = ["vector_top", "vector_bottom", "insert_template"].map(|role| TopoSource {
            role: role.into(),
            source: Source::Input(role.into()),
        });
        let inputs = collagene_core::topo_run_preparation::source_inputs(
            &roles,
            &[&top, &bottom, &template],
            &compiled.request,
            "insert_template",
            Some(&origin),
        )
        .unwrap();
        assert_eq!(inputs[0].molecule, top);
        assert_eq!(inputs[1].molecule, bottom);
        assert_eq!(inputs[2].molecule, template);
        let task = collagene_core::topo_run_preparation::Task::Orientations {
            request: compiled.request.clone(),
            vector: compiled.vector_strands,
            insert: compiled.insert_strands,
            directional: compiled.directional,
        };
        assert_eq!(
            task_revision(&inputs, Some(&origin), &task).unwrap(),
            collagene_core::topo_run_preparation::documentary_revision(
                &inputs,
                "insert_template",
                &origin,
                &task
            )
            .unwrap()
        );
        let wrong = Molecule::dna(
            "Equal DNA, different native revision",
            "CACCATGACG",
            Topology::Linear,
        )
        .unwrap();
        assert!(
            collagene_core::topo_run_preparation::source_inputs(
                &roles,
                &[&top, &bottom, &wrong],
                &compiled.request,
                "insert_template",
                Some(&origin)
            )
            .is_err()
        );
        let mut wrong_event = product;
        wrong_event.event_id = format!("event-{}", "f".repeat(64));
        assert!(
            run_material::inspect_material(
                &installed.workspace,
                &wrong_event,
                &Some(MaterialSelection::Derived(wrong_event.clone()))
            )
            .is_err()
        );
        assert_eq!(
            serde_json::to_vec(&installed.workspace).unwrap(),
            original_workspace
        );
    }
}
