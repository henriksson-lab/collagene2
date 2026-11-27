//! Explicit TA/directional preparation form over the shared native builder.
use collagene_core::{Sequence, analysis::Primer, tailing::TailingBase, topo};
use topo::physical::{self as physical, projection};

#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub enum SourceMode {
    #[default]
    Original,
    SavedVector,
    ResumeGraph,
}
#[derive(Clone, Debug, PartialEq, Eq, Default)]
pub struct Draft {
    pub source_mode: SourceMode,
    pub graph: Option<collagene_core::topo_graph_history::Selection>,
    pub vector_strands: [String; 2],
    pub insert_strands: [String; 2],
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

pub(crate) fn active_graph(draft: &Draft) -> Option<collagene_core::topo_graph_history::Selection> {
    if draft.source_mode == SourceMode::Original {
        None
    } else {
        draft.graph.clone()
    }
}
pub fn external_sources(draft: &Draft) -> Vec<(String, String)> {
    match draft.source_mode {
        SourceMode::Original => vec![
            ("vector_top".into(), draft.vector_top.clone()),
            ("vector_bottom".into(), draft.vector_bottom.clone()),
            ("insert_template".into(), draft.insert_template.clone()),
        ],
        SourceMode::SavedVector => vec![("insert_template".into(), draft.insert_template.clone())],
        SourceMode::ResumeGraph => vec![],
    }
}
pub fn original_vector(
    draft: &Draft,
) -> Result<collagene_core::topo_preparation::PreparedVectorDeclaration, String> {
    use collagene_core::topo_preparation::{PreparationMode, PreparedVectorDeclaration};
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
    let (paired, bottom_start, declaration_id) = match draft.mode {
        PreparationMode::Ta => {
            if top < 2 || top != bottom || top > 100_000 {
                return Err("TA requires equal vector strand lengths of2–100,000 bases, with one explicit3′ overhang.".into());
            }
            (top - 1, 0, "user-declared-ta-vector-1")
        }
        PreparationMode::DirectionalCacc => {
            if top == 0 || top > 99_996 || top.checked_add(4) != Some(bottom) {
                return Err("Directional preparation requires actual top/bottom strand sources differing by four bases: bottom5′ GGTG, paired from bottom offset4. No strand is synthesized from the other.".into());
            }
            (top, 4, "user-declared-directional-vector-1")
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
    Ok(PreparedVectorDeclaration {
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
    })
}
pub fn initial_vector_request(draft: &Draft) -> Result<topo::Request, String> {
    let vector = original_vector(draft)?;
    Ok(topo::Request {
        recipe: topo::GenerationRecipe {
            schema_version: 1,
            materials: vec![topo::MaterialRecipe::Prepared {
                declaration_id: vector.declaration_id,
                assumption: vector.assumption,
                strands: vector.strands,
                pairing: vector.pairing,
            }],
        },
        transfers: vec![],
        engine_limits: physical::Limits::default(),
        projection_options: projection::Options::default(),
    })
}
fn exact_insert(
    draft: &Draft,
) -> Result<collagene_core::topo_preparation::ExactPcrDeclaration, String> {
    use collagene_core::topo_preparation::{ExactPcrDeclaration, PreparationMode};
    let (material_id, tail) = match draft.mode {
        PreparationMode::Ta => {
            if !draft.add_a_tail {
                return Err("TA requires an explicit A addition on both PCR strands".into());
            }
            ("explicit-pcr-a-tail-insert", Some(TailingBase::A))
        }
        PreparationMode::DirectionalCacc => {
            if draft.add_a_tail {
                return Err(
                    "Directional CACC uses blunt exact PCR; disable A addition explicitly".into(),
                );
            }
            ("explicit-blunt-pcr-insert", None)
        }
    };
    Ok(ExactPcrDeclaration {
        material_id: material_id.into(),
        template_role: "insert_template".into(),
        forward: primer("Forward", &draft.forward_sequence, draft.forward_annealing)?,
        reverse: primer("Reverse", &draft.reverse_sequence, draft.reverse_annealing)?,
        tail,
    })
}
fn handles(raw: &[String; 2], label: &str) -> Result<[usize; 2], String> {
    let parse = |text: &str| {
        if text.is_empty() || text.len() > 8 || !text.bytes().all(|b| b.is_ascii_digit()) {
            Err(format!("Choose an existing {label} strand"))
        } else {
            text.parse::<usize>()
                .map_err(|_| format!("Invalid {label} strand"))
        }
    };
    let value = [parse(&raw[0])?, parse(&raw[1])?];
    if value[0] == value[1] {
        return Err(format!("Choose two different {label} strands"));
    }
    Ok(value)
}
pub fn preparation_request(
    draft: &Draft,
    graph: Option<&collagene_core::topo_graph_history::Capsule>,
) -> Result<
    (
        Vec<(String, String)>,
        collagene_core::topo_preparation::PreparationRequest,
    ),
    String,
> {
    use collagene_core::topo_preparation::{PreparationRequest, PreparationSource};
    let roles = external_sources(draft);
    if roles
        .iter()
        .any(|(_, id)| id.trim().is_empty() || id.len() > 1024 || id.chars().any(char::is_control))
    {
        return Err("Choose every requested source document".into());
    }
    let source = match draft.source_mode {
        SourceMode::Original => {
            if graph.is_some() {
                return Err("Original declaration cannot consume a saved graph".into());
            }
            PreparationSource::Original {
                vector: original_vector(draft)?,
                insert: exact_insert(draft)?,
            }
        }
        SourceMode::SavedVector | SourceMode::ResumeGraph => {
            let origin = graph.ok_or("Select and inspect a saved graph first")?;
            if draft.graph.as_ref() != Some(&origin.selection) {
                return Err("The exact saved graph selection changed".into());
            }
            let vector_strands = handles(&draft.vector_strands, "vector")?;
            match draft.source_mode {
                SourceMode::SavedVector => PreparationSource::SavedVector {
                    material_id: "saved-vector-graph".into(),
                    origin: origin.clone(),
                    vector_strands,
                    insert: exact_insert(draft)?,
                },
                SourceMode::ResumeGraph => PreparationSource::ResumeGraph {
                    material_id: "saved-reaction-graph".into(),
                    origin: origin.clone(),
                    vector_strands,
                    insert_strands: handles(&draft.insert_strands, "existing insert")?,
                },
                SourceMode::Original => {
                    return Err("Choose the active graph source mode explicitly".into());
                }
            }
        }
    };
    Ok((
        roles,
        PreparationRequest {
            mode: draft.mode,
            source,
            engine_limits: physical::Limits::default(),
            projection_options: projection::Options::default(),
        },
    ))
}
pub fn compile_request(
    draft: &Draft,
    graph: Option<&collagene_core::topo_graph_history::Capsule>,
) -> Result<
    (
        Vec<(String, String)>,
        collagene_core::topo_preparation::CompiledPreparation,
    ),
    String,
> {
    let (roles, request) = preparation_request(draft, graph)?;
    Ok((
        roles,
        collagene_core::topo_preparation::compile(&request).map_err(|e| e.to_string())?,
    ))
}

/// Exact transient request; capsules are hydrated only during explicit capture.
pub struct CaptureRequest {
    pub draft: Draft,
    pub sources: Vec<(String, collagene_core::run_lineage::EventBinding)>,
    pub target: collagene_core::run_lineage::EventBinding,
    pub request: Option<topo::Request>,
    pub template_origin: Option<collagene_core::run_material::MaterialSelection>,
}
pub struct Captured {
    pub inputs: Vec<topo::Input>,
    pub target: collagene_core::Molecule,
    pub template: Option<collagene_core::topo_run_preparation::BoundTemplate>,
    pub preparation: collagene_core::topo_preparation::CompiledPreparation,
}
fn declaration_hash(value: &impl serde::Serialize) -> Result<String, String> {
    struct Hash {
        hash: blake3::Hasher,
        count: usize,
    }
    impl std::io::Write for Hash {
        fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
            self.count = self
                .count
                .checked_add(bytes.len())
                .filter(|n| *n <= 32 * 1024 * 1024)
                .ok_or_else(|| std::io::Error::other("TOPO declaration exceeds32MiB"))?;
            self.hash.update(bytes);
            Ok(bytes.len())
        }
        fn flush(&mut self) -> std::io::Result<()> {
            Ok(())
        }
    }
    let mut out = Hash {
        hash: blake3::Hasher::new(),
        count: 0,
    };
    serde_json::to_writer(&mut out, value).map_err(|e| e.to_string())?;
    Ok(out.hash.finalize().to_hex().to_string())
}
pub fn capture(
    workspace: &collagene_core::workspace::Workspace,
    chosen: &CaptureRequest,
) -> Result<Captured, String> {
    use collagene_core::{
        Alphabet, run_material,
        topo_run_preparation::{self, BoundTemplate},
        verification::{Source, TopoSource},
    };
    workspace.validate().map_err(|e| e.to_string())?;
    let resolve=|binding:&collagene_core::run_lineage::EventBinding|->Result<&collagene_core::Molecule,String>{
        let doc=workspace.documents.iter().find(|doc|doc.id==binding.document_id).ok_or("A selected source document is unavailable")?;
        if doc.head!=binding.revision||doc.head_event.as_ref()!=Some(&binding.event_id){return Err("A selected source event changed; inspect again".into());}
        workspace.head(&doc.id).map_err(|e|e.to_string())
    };
    if workspace.active_document != chosen.target.document_id {
        return Err("The independent target changed".into());
    }
    let target = resolve(&chosen.target)?;
    if target.sequence.alphabet() != Alphabet::Dna || target.sequence.len() > 100_000 {
        return Err("Choose an independent DNA target within100000 bases".into());
    }
    let expected = external_sources(&chosen.draft);
    if expected.len() != chosen.sources.len()
        || chosen.sources.len() > 32
        || expected
            .iter()
            .zip(&chosen.sources)
            .any(|((role, id), (given, event))| role != given || id != &event.document_id)
    {
        return Err("The exact external source roles changed".into());
    }
    let mut molecules = Vec::new();
    let mut roles = Vec::new();
    for (role, binding) in &chosen.sources {
        molecules.push(resolve(binding)?);
        roles.push(TopoSource {
            role: role.clone(),
            source: Source::Input(role.clone()),
        });
    }
    let mut all = molecules.clone();
    all.push(target);
    collagene_core::model::bounds::model_preflight(&all).map_err(|e| e.to_string())?;
    if let Some(selection) = &chosen.template_origin {
        if chosen.draft.source_mode == SourceMode::ResumeGraph {
            return Err("Continuing a saved graph does not consume a new PCR template".into());
        }
        let (_, event) = chosen
            .sources
            .iter()
            .find(|(role, _)| role == "insert_template")
            .ok_or("Select the exact insert template")?;
        if selection.event() != event {
            return Err("Saved insert origin belongs to another exact source event".into());
        }
    }
    let graph_selections = active_graph(&chosen.draft).into_iter().collect::<Vec<_>>();
    let material_selections = chosen
        .template_origin
        .clone()
        .into_iter()
        .collect::<Vec<_>>();
    let closed =
        run_material::documentary_origins(workspace, &material_selections, &graph_selections)
            .map_err(|e| e.to_string())?;
    if closed.materials.len() != material_selections.len()
        || closed.graphs.len() != graph_selections.len()
    {
        return Err("Native source selections changed during capture".into());
    }
    let graph = closed.graphs.into_iter().next();
    let template = closed
        .materials
        .into_iter()
        .next()
        .map(|origin| BoundTemplate {
            role: "insert_template".into(),
            origin,
        });
    let mut draft = chosen.draft.clone();
    if draft.source_mode == SourceMode::Original {
        draft.vector_top_length = molecules[0].sequence.len();
        draft.vector_bottom_length = molecules[1].sequence.len();
    }
    let (_, mut preparation) = compile_request(&draft, graph.as_ref())?;
    if let Some(request) = &chosen.request {
        if declaration_hash(&(
            &request.recipe,
            &request.engine_limits,
            &request.projection_options,
        ))? != declaration_hash(&(
            &preparation.request.recipe,
            &preparation.request.engine_limits,
            &preparation.request.projection_options,
        ))? {
            return Err("Inspected request differs from the exact current preparation".into());
        }
        preparation.request = request.clone();
    }
    let inputs = topo_run_preparation::source_inputs(
        &roles,
        &molecules,
        &preparation.request,
        template.as_ref(),
    )
    .map_err(|e| e.to_string())?;
    Ok(Captured {
        inputs,
        target: target.clone(),
        template,
        preparation,
    })
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
    if next.source_mode == SourceMode::ResumeGraph {
        return Err("Continuing a saved reaction does not synthesize new PCR primers".into());
    }
    exact_insert(&next)?;
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
            ..Default::default()
        };
        let (roles, compiled) = compile_request(&draft, None).unwrap();
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
        assert!(compile_request(&missing, None).is_err());
        missing = draft;
        missing.forward_annealing = 7;
        assert!(compile_request(&missing, None).is_err());
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
            ..Default::default()
        }
    }
    #[test]
    fn directional_form_declares_actual_geometry_and_never_reuses_ta_or_adds_a_tail() {
        let draft = directional();
        let original = draft.clone();
        let (roles, compiled) = compile_request(&draft, None).unwrap();
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
        assert!(compile_request(&invalid, None).is_err());
        invalid = draft.clone();
        invalid.vector_bottom_length = 14;
        assert!(compile_request(&invalid, None).is_err());
        invalid = draft;
        invalid.activation_confirmed = false;
        assert!(compile_request(&invalid, None).is_err());
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
        let (_, request) = preparation_request(&proposed, None).unwrap();
        let collagene_core::topo_preparation::PreparationSource::Original { insert, .. } =
            request.source
        else {
            panic!("explicit original source expected")
        };
        assert_eq!(insert.forward.tail, "CACC");
        assert_eq!(insert.forward.annealing.as_str(), "CATGA");
        let mut too_long = draft;
        too_long.forward_sequence = "A".repeat(collagene_core::primers::MAX_OLIGO);
        too_long.forward_annealing = 20;
        assert!(cacc_prefix_proposal(&too_long).is_err());
    }
}

/// Documentary inner context for the single native TOPO task endpoint.
pub(crate) fn task_revision(
    inputs: &[topo::Input],
    template: Option<&collagene_core::topo_run_preparation::BoundTemplate>,
    task: &collagene_core::topo_run_preparation::Task,
) -> Result<String, collagene_core::Error> {
    collagene_core::topo_run_preparation::documentary_revision(inputs, template, task)
}
pub(crate) fn batch_revision(
    inputs: &[topo::Input],
    template: Option<&collagene_core::topo_run_preparation::BoundTemplate>,
    request: &topo::Request,
) -> Result<String, collagene_core::Error> {
    collagene_core::topo_run_preparation::preparation_revision(inputs, template, request)
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
    template: Option<collagene_core::topo_run_preparation::BoundTemplate>,
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
        topo_run_preparation::{self, Reply},
    };
    let revision = topo_run_preparation::revision(&inputs, template.as_ref(), &task)
        .map_err(|e| e.to_string())?;
    let has_origin = template.is_some()
        || task
            .preparation()
            .recipe
            .materials
            .iter()
            .any(|material| matches!(material, topo::MaterialRecipe::SavedGraph { .. }));
    let response = worker
        .ask(Request {
            id,
            revision: revision.clone(),
            job: AnalysisJob::InspectTopoRun {
                inputs: inputs.clone(),
                template: template.clone(),
                task: task.clone(),
            },
        })
        .await?;
    if response.id != id || response.revision != revision {
        return Err("TOPO response belongs to another exact source context".into());
    }
    let report = match response.result {
        AnalysisResult::TopoRunInspection(report) => report,
        AnalysisResult::Error(e) => return Err(e),
        _ => return Err("Unexpected native TOPO task response".into()),
    };
    report
        .validate_for(&inputs, template.as_ref(), &task, &revision)
        .map_err(|e| e.to_string())?;
    let review = has_origin.then(|| OriginReview {
        incomplete: report.origin_incomplete,
        findings: report.origin_findings,
        total: report.origin_findings_total,
        truncated: report.origin_findings_truncated,
    });
    Ok((
        match report.reply {
            Reply::Orientations(r) => AnalysisResult::TopoOrientations(r),
            Reply::Projected(r) => AnalysisResult::TopoProjection(r),
            Reply::Checked(r) => AnalysisResult::TopoTarget(r),
            Reply::ActionsInspected(r) => AnalysisResult::TopoActionsInspected(r),
            Reply::ActionsProcessed(r) => AnalysisResult::TopoActionsProduct(r),
        },
        review,
    ))
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
        let preview = verification::simulate_pcr_documents(
            &parent,
            &forward,
            &reverse,
            collagene_core::chemistry::Strand::Top,
        )
        .unwrap()
        .remove(0);
        let installed = derivation::prepare_install(
            &workspace,
            &derivation::Request {
                source: event,
                origin: derivation::SourceOrigin::OriginalReference,
                derivation: derivation::Derivation::ExactPcr {
                    template_strand: collagene_core::chemistry::Strand::Top,
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
        let saved_template = collagene_core::topo_run_preparation::BoundTemplate {
            role: "insert_template".into(),
            origin: origin.clone(),
        };
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
            ..Default::default()
        };
        let (_, compiled) = compile_request(&draft, None).unwrap();
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
            Some(&saved_template),
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
            task_revision(&inputs, Some(&saved_template), &task).unwrap(),
            collagene_core::topo_run_preparation::documentary_revision(
                &inputs,
                Some(&saved_template),
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
                Some(&saved_template)
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
