//! Explicit preparation graphs for TA/GC ligation and oriented overlap profiles.
use crate::{app::download, worker::CachedWorker};
use collagene_core::{
    Molecule,
    analysis::Primer,
    chemistry::{EndRepairMethod, EndTreatment, FivePrime, Terminus, ThreePrime},
    overlap_profiles::{OverlapMethod, ProfileOptions, ProfileReport},
    protocol::{AnalysisJob, AnalysisResult, Request},
    tailing::TailingBase,
    transfer::Candidate,
    verification::{
        DesignPlan, Operation, OrientedSource, PlanInput, Source, TailPreview, TailSubstrate,
    },
};
use leptos::{html, prelude::*};
use serde::{Deserialize, Serialize};
use std::{
    rc::Rc,
    sync::{
        Arc,
        atomic::{AtomicBool, Ordering},
    },
};
use wasm_bindgen_futures::spawn_local;

#[derive(Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum Mode {
    Tailing,
    Overlap,
    Ligation,
    GoldenGate,
}
#[derive(Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum Treatment {
    Keep,
    Phosphorylate,
    Dephosphorylate,
}
#[derive(Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct TailInput {
    pub document: usize,
    pub pcr: bool,
    pub enzymes: String,
    pub fragment: String,
    pub blunt: bool,
    pub forward_name: String,
    pub reverse_name: String,
    pub forward: String,
    pub reverse: String,
    pub forward_tail: String,
    pub reverse_tail: String,
    pub base: TailingBase,
    pub treatment: Treatment,
    pub reversed: bool,
}
impl TailInput {
    fn new(pcr: bool, base: TailingBase) -> Self {
        Self {
            document: 0,
            pcr,
            enzymes: "EcoRV".into(),
            fragment: "1".into(),
            blunt: false,
            forward_name: "Forward".into(),
            reverse_name: "Reverse".into(),
            forward: String::new(),
            reverse: String::new(),
            forward_tail: String::new(),
            reverse_tail: String::new(),
            base,
            treatment: Treatment::Keep,
            reversed: false,
        }
    }
    fn recipe(&self) -> Result<TailSubstrate, String> {
        if self.pcr {
            Ok(TailSubstrate::ExactPcr {
                forward: Primer::new(&self.forward_name, &self.forward, &self.forward_tail)
                    .map_err(|e| e.to_string())?,
                reverse: Primer::new(&self.reverse_name, &self.reverse, &self.reverse_tail)
                    .map_err(|e| e.to_string())?,
            })
        } else {
            let enzymes = self
                .enzymes
                .split(',')
                .map(str::trim)
                .filter(|s| !s.is_empty())
                .map(str::to_owned)
                .collect::<Vec<_>>();
            if enzymes.is_empty() {
                return Err("Choose explicit restriction enzymes".into());
            }
            let fragment_index = self
                .fragment
                .parse::<usize>()
                .ok()
                .and_then(|n| n.checked_sub(1))
                .ok_or("Fragment numbers start at 1")?;
            Ok(TailSubstrate::Digest {
                enzymes,
                fragment_index,
                blunt_before: self.blunt,
            })
        }
    }
}
#[derive(Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum OverlapPreparation {
    WholeLinear,
    DeclaredDuplex,
    Digest,
    CatalogDigest,
    ExactPcr,
    PcrCatalogDigest,
}
#[derive(Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct OverlapInput {
    pub document: usize,
    pub interval_start: String,
    pub interval_end: String,
    pub reverse: bool,
    pub treatment: Treatment,
    pub preparation: OverlapPreparation,
    #[serde(default)]
    pub catalog_selection: Option<collagene_core::catalog_digest::CatalogSelection>,
    pub left_terminus: Terminus,
    pub right_terminus: Terminus,
    pub enzymes: String,
    pub fragment: String,
    pub repair: Option<EndRepairMethod>,
    pub forward_name: String,
    pub reverse_name: String,
    pub forward: String,
    pub reverse_primer: String,
    pub forward_tail: String,
    pub reverse_tail: String,
    #[serde(default)]
    pub right_fusion: String,
}
impl Default for OverlapInput {
    fn default() -> Self {
        Self {
            document: 0,
            interval_start: "1".into(),
            interval_end: String::new(),
            reverse: false,
            treatment: Treatment::Keep,
            preparation: OverlapPreparation::WholeLinear,
            catalog_selection: None,
            left_terminus: Terminus::UNKNOWN,
            right_terminus: Terminus::UNKNOWN,
            enzymes: "EcoRV".into(),
            fragment: "1".into(),
            repair: None,
            forward_name: "Forward".into(),
            reverse_name: "Reverse".into(),
            forward: String::new(),
            reverse_primer: String::new(),
            forward_tail: String::new(),
            reverse_tail: String::new(),
            right_fusion: String::new(),
        }
    }
}
impl OverlapInput {
    pub(crate) fn summary(&self) -> String {
        let preparation = match self.preparation {
            OverlapPreparation::WholeLinear => {
                "Whole linear source (circular DNA requires preparation)".into()
            }
            OverlapPreparation::DeclaredDuplex => format!(
                "Declared intact linear duplex assumption; left {}; right {}",
                declared_terminus(self.left_terminus),
                declared_terminus(self.right_terminus)
            ),
            OverlapPreparation::CatalogDigest | OverlapPreparation::PcrCatalogDigest => format!(
                "{}: {}; fragment {}",
                if self.preparation == OverlapPreparation::PcrCatalogDigest {
                    "Exact PCR followed by pinned catalog digest"
                } else {
                    "Pinned catalog digest"
                },
                self.catalog_selection
                    .as_ref()
                    .map(|s| s
                        .enzymes
                        .iter()
                        .map(|e| e.enzyme.name.as_str())
                        .collect::<Vec<_>>()
                        .join(", "))
                    .unwrap_or_else(|| "select a virtual digest fragment".into()),
                self.fragment
            ),
            OverlapPreparation::Digest => format!(
                "Restriction digest: {}; fragment {} (1-based); {}",
                self.enzymes,
                self.fragment,
                match self.repair {
                    None => "cohesive ends retained",
                    Some(EndRepairMethod::FillFivePrimeOverhangs) => "explicit 5′ overhang fill-in",
                    Some(EndRepairMethod::RemoveThreePrimeOverhangs) =>
                        "explicit 3′ overhang removal",
                    Some(EndRepairMethod::Blunt) => "explicit blunting",
                }
            ),
            OverlapPreparation::ExactPcr => format!(
                "Exact PCR: annealing {} / {} bp; 5′ extensions {} / {} bp",
                self.forward.len(),
                self.reverse_primer.len(),
                self.forward_tail.len(),
                self.reverse_tail.len()
            ),
        };
        format!(
            "{} · {}",
            preparation,
            if self.reverse {
                "reverse complement after preparation"
            } else {
                "forward after preparation"
            }
        )
    }
}
#[derive(Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct MethodsState {
    pub mode: Mode,
    pub vector: TailInput,
    pub insert: TailInput,
    pub circular: bool,
    pub allow_nicked: bool,
    #[serde(deserialize_with = "deserialize_method_rows")]
    pub rows: Vec<OverlapInput>,
    pub overlaps: String,
    pub method: OverlapMethod,
    pub scan_repeats: bool,
    pub comparisons: String,
    pub occurrences: String,
    #[serde(default)]
    pub golden_gate: GoldenGateForm,
}
#[derive(Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct GoldenGateForm {
    pub enzyme: collagene_core::golden_gate::GoldenGateEnzyme,
    pub clamp: String,
    pub spacer: String,
    pub forward_length: String,
    pub reverse_length: String,
}
impl Default for GoldenGateForm {
    fn default() -> Self {
        Self {
            enzyme: collagene_core::golden_gate::GoldenGateEnzyme::BsaI,
            clamp: "ATGC".into(),
            spacer: "A".into(),
            forward_length: "20".into(),
            reverse_length: "20".into(),
        }
    }
}
fn deserialize_method_rows<'de, D: serde::Deserializer<'de>>(
    deserializer: D,
) -> Result<Vec<OverlapInput>, D::Error> {
    struct Visitor;
    impl<'de> serde::de::Visitor<'de> for Visitor {
        type Value = Vec<OverlapInput>;
        fn expecting(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
            f.write_str("at most12 assembly preparation rows")
        }
        fn visit_seq<A: serde::de::SeqAccess<'de>>(
            self,
            mut seq: A,
        ) -> Result<Self::Value, A::Error> {
            if seq.size_hint().is_some_and(|n| n > 12) {
                return Err(serde::de::Error::custom(
                    "at most12 assembly preparation rows",
                ));
            }
            let mut rows = Vec::with_capacity(2);
            while rows.len() < 12 {
                match seq.next_element()? {
                    Some(row) => rows.push(row),
                    None => return Ok(rows),
                }
            }
            if seq.next_element::<serde::de::IgnoredAny>()?.is_some() {
                return Err(serde::de::Error::custom(
                    "at most12 assembly preparation rows",
                ));
            }
            Ok(rows)
        }
    }
    deserializer.deserialize_seq(Visitor)
}
impl Default for MethodsState {
    fn default() -> Self {
        Self {
            mode: Mode::Overlap,
            vector: TailInput::new(false, TailingBase::T),
            insert: TailInput::new(true, TailingBase::A),
            circular: true,
            allow_nicked: false,
            rows: vec![OverlapInput::default(), OverlapInput::default()],
            overlaps: "20".into(),
            method: OverlapMethod::HiFi,
            scan_repeats: true,
            comparisons: "5000000".into(),
            occurrences: "1000".into(),
            golden_gate: GoldenGateForm::default(),
        }
    }
}
impl MethodsState {
    fn input(&self, vector: bool) -> &TailInput {
        if vector { &self.vector } else { &self.insert }
    }
    fn input_mut(&mut self, vector: bool) -> &mut TailInput {
        if vector {
            &mut self.vector
        } else {
            &mut self.insert
        }
    }
    pub(crate) fn profile(&self) -> Result<(Vec<usize>, ProfileOptions), String> {
        let overlaps = self
            .overlaps
            .split(',')
            .map(|s| {
                s.trim()
                    .parse::<usize>()
                    .map_err(|_| "Supply comma-separated overlap lengths".to_string())
            })
            .collect::<Result<Vec<_>, _>>()?;
        Ok((overlaps, self.profile_options()?))
    }
    pub(crate) fn profile_options(&self) -> Result<ProfileOptions, String> {
        let options = ProfileOptions {
            method: self.method,
            scan_repeated_overlaps: self.scan_repeats,
            max_comparisons: self
                .comparisons
                .parse()
                .map_err(|_| "Invalid comparison budget")?,
            max_occurrences: self
                .occurrences
                .parse()
                .map_err(|_| "Invalid occurrence limit")?,
        };
        options.validate().map_err(|e| e.to_string())?;
        Ok(options)
    }
}
#[derive(Clone, PartialEq, Eq, Serialize, Deserialize)]
struct PlannerForm {
    overlap: String,
    minimum: String,
    maximum: String,
    tm: String,
}
impl Default for PlannerForm {
    fn default() -> Self {
        Self {
            overlap: "20".into(),
            minimum: "18".into(),
            maximum: "30".into(),
            tm: "60".into(),
        }
    }
}
/// Actual document identity or an explicitly reviewed reference intention.
/// An intent has no invented document/event ancestry.
#[derive(Clone, PartialEq, Eq, Serialize, Deserialize)]
pub enum TargetContext {
    Document(crate::digest::SourceInfo),
    Intent {
        context_revision: String,
        revision: String,
    },
}
impl TargetContext {
    pub fn revision(&self) -> &str {
        match self {
            Self::Document(info) => &info.binding.revision,
            Self::Intent { revision, .. } => revision,
        }
    }
}
struct Suggestions {
    report: collagene_core::overlap_primer_design::Report,
    config: MethodsState,
    form: PlannerForm,
    context: Vec<crate::digest::SourceInfo>,
    target: String,
    target_context: TargetContext,
    epoch: u64,
}
struct GoldenGateSuggestions {
    report: collagene_core::golden_gate_primer_design::Report,
    config: MethodsState,
    context: Vec<crate::digest::SourceInfo>,
    target_context: TargetContext,
    epoch: u64,
}
#[derive(Clone)]
struct Evidence {
    config: MethodsState,
    revisions: Vec<String>,
    context: Vec<crate::digest::SourceInfo>,
    target_context: TargetContext,
    target: String,
    epoch: u64,
    tails: Vec<TailPreview>,
    profile: Option<ProfileReport>,
    overlap: Option<collagene_core::assembly_recipe::Preview>,
    ligation: Option<collagene_core::ligation_recipe::Preview>,
}
fn guidance(minimum: Option<usize>, maximum: Option<usize>) -> String {
    match (minimum, maximum) {
        (Some(a), Some(b)) => format!("{a}–{b} bp"),
        (Some(a), None) => format!("at least {a} bp; no upper bound specified"),
        _ => "not established".into(),
    }
}
fn short(sequence: &str) -> String {
    let s = sequence.chars().take(300).collect::<String>();
    if s.len() < sequence.len() {
        format!("{s}…")
    } else {
        s
    }
}
fn append_preparation(
    steps: &mut Vec<Operation>,
    input_id: &str,
    input: &TailInput,
) -> Result<Source, String> {
    let source = Source::Input(input_id.into());
    match input.recipe()? {
        TailSubstrate::Digest {
            enzymes,
            fragment_index,
            blunt_before,
        } => {
            steps.push(Operation::Digest {
                template: source,
                enzymes,
                fragment_index,
            });
            if blunt_before {
                let template = Source::Step(steps.len() - 1);
                steps.push(Operation::RepairEnds {
                    template,
                    method: EndRepairMethod::Blunt,
                });
            }
        }
        TailSubstrate::ExactPcr { forward, reverse } => steps.push(Operation::Pcr {
            template: source,
            forward,
            reverse,
        }),
    };
    let template = Source::Step(steps.len() - 1);
    steps.push(Operation::Tail3Prime {
        template,
        base: input.base,
    });
    let mut source = Source::Step(steps.len() - 1);
    let treatment = match input.treatment {
        Treatment::Keep => None,
        Treatment::Phosphorylate => Some(EndTreatment::Phosphorylate),
        Treatment::Dephosphorylate => Some(EndTreatment::Dephosphorylate),
    };
    if let Some(treatment) = treatment {
        steps.push(Operation::TreatEnds {
            template: source,
            treatment,
        });
        source = Source::Step(steps.len() - 1);
    }
    Ok(source)
}
async fn ask(
    worker: &CachedWorker,
    job: AnalysisJob,
    revision: String,
    id: u32,
) -> Result<AnalysisResult, String> {
    let response = worker.ask(Request { id, revision, job }).await?;
    match response.result {
        AnalysisResult::Error(e) => Err(e),
        result => Ok(result),
    }
}
#[component]
pub fn CloningMethodsPanel(
    molecule: Memo<Molecule>,
    query_options: Memo<Vec<(usize, String)>>,
    get_query: Callback<usize, Option<Molecule>>,
    source_info: Callback<usize, Option<crate::digest::SourceInfo>>,
    target_info: Memo<TargetContext>,
    selected_region: Callback<(), Option<(usize, usize, usize)>>,
    config: RwSignal<MethodsState>,
    epoch: RwSignal<u64>,
    on_candidate: Callback<Candidate>,
    on_products: Callback<Vec<Molecule>, Result<(), String>>,
    primer_sources: Memo<Vec<crate::digest::SourceInfo>>,
    list_primers: Callback<
        collagene_core::run_lineage::EventBinding,
        Result<collagene_core::primer_reuse::List, String>,
    >,
    resolve_primer: Callback<
        crate::primer_picker::ResolveRequest,
        Result<collagene_core::primer_reuse::Resolved, String>,
    >,
    #[prop(optional_no_strip)] on_finalize: Option<Callback<Candidate, Result<(), String>>>,
    #[prop(optional_no_strip)] on_source: Option<Callback<(usize, usize), Result<(), String>>>,
    #[prop(default = false)] fixed_sources: bool,
    #[prop(default = false)] disable_plain_primer_design: bool,
) -> impl IntoView {
    let row_nodes = StoredValue::new_local(
        (0..12)
            .map(|_| NodeRef::<html::Div>::new())
            .collect::<Vec<_>>(),
    );
    let source_indices = Memo::new(move |_| {
        config.with(|c| {
            if c.mode == Mode::Tailing {
                vec![c.vector.document, c.insert.document]
            } else {
                c.rows.iter().map(|r| r.document).collect()
            }
        })
    });
    let source_context = Memo::new(move |_| {
        source_indices
            .get()
            .into_iter()
            .map(|i| source_info.run(i))
            .collect::<Option<Vec<_>>>()
    });
    let sources = Memo::new(move |_| {
        source_indices
            .get()
            .into_iter()
            .map(|i| get_query.run(i))
            .collect::<Option<Vec<_>>>()
    });
    let source_revisions = Memo::new(move |_| {
        sources.with(|s| {
            s.as_ref()
                .map(|m| m.iter().map(Molecule::revision).collect::<Vec<_>>())
        })
    });
    let target = Memo::new(move |_| molecule.with(Molecule::revision));
    let busy = RwSignal::new(false);
    let status = RwSignal::new(String::new());
    let generation = RwSignal::new(0u32);
    let evidence = RwSignal::new(None::<Arc<Evidence>>);
    let planner_form = RwSignal::new(PlannerForm::default());
    let suggestions = RwSignal::new(None::<Arc<Suggestions>>);
    let gate_suggestions = RwSignal::new(None::<Arc<GoldenGateSuggestions>>);

    let active = StoredValue::new_local(None::<Rc<CachedWorker>>);
    let alive = Arc::new(AtomicBool::new(true));
    let cleanup = alive.clone();
    on_cleanup(move || {
        cleanup.store(false, Ordering::Relaxed);
        if let Some(w) = active.get_value() {
            w.close();
        }
        active.set_value(None);
    });
    Effect::new(move |_| {
        let _ = (
            config.get(),
            source_revisions.get(),
            source_context.get(),
            target_info.get(),
            planner_form.get(),
            target.get(),
            epoch.get(),
        );
        generation.update(|n| *n = n.wrapping_add(1));
        evidence.set(None);
        suggestions.set(None);
        gate_suggestions.set(None);
        if let Some(w) = active.get_value() {
            w.close();
        }
        active.set_value(None);
        busy.set(false);
        status.set(String::new());
    });
    let gate_alive = alive.clone();
    let suggest_gate = Callback::new(move |()| {
        if busy.get_untracked() {
            return;
        }
        let state = config.get_untracked();
        if state.mode != Mode::GoldenGate {
            return;
        }
        let (Some(originals), Some(context)) =
            (sources.get_untracked(), source_context.get_untracked())
        else {
            status.set("Choose the original source documents.".into());
            return;
        };
        let target_molecule = molecule.get_untracked();
        let target_context = target_info.get_untracked();
        let captured_epoch = epoch.get_untracked();
        let (inputs, request) = match golden_gate_request(&state, originals, &target_molecule) {
            Ok(value) => value,
            Err(error) => {
                status.set(error);
                return;
            }
        };
        let revisions = inputs
            .iter()
            .map(|p| p.revision.clone())
            .collect::<Vec<_>>();
        if revisions
            != context
                .iter()
                .map(|s| s.binding.revision.clone())
                .collect::<Vec<_>>()
            || target_context.revision() != target_molecule.revision()
        {
            status.set("Source or target identity changed before primer planning.".into());
            return;
        }
        let revision = match collagene_core::golden_gate_primer_design::preview_revision(
            &inputs,
            &target_molecule,
            &request,
        ) {
            Ok(value) => value,
            Err(error) => {
                status.set(error.to_string());
                return;
            }
        };
        let worker = match CachedWorker::new() {
            Ok(worker) => Rc::new(worker),
            Err(error) => {
                status.set(error);
                return;
            }
        };
        if let Some(old) = active.get_value() {
            old.close();
        }
        active.set_value(Some(worker.clone()));
        let token = generation.get_untracked().wrapping_add(1);
        generation.set(token);
        busy.set(true);
        evidence.set(None);
        suggestions.set(None);
        gate_suggestions.set(None);
        status.set(
            "Proposing explicit Type IIS PCR tails and freshly checking their digested ends…"
                .into(),
        );
        let alive = gate_alive.clone();
        spawn_local(async move {
            let result = ask(
                &worker,
                AnalysisJob::DesignGoldenGatePrimers {
                    inputs,
                    target: target_molecule,
                    request,
                },
                revision.clone(),
                token,
            )
            .await;
            worker.close();
            if !alive.load(Ordering::Relaxed)
                || generation.get_untracked() != token
                || config.get_untracked() != state
                || source_context.get_untracked() != Some(context.clone())
                || source_revisions.get_untracked() != Some(revisions.clone())
                || target_info.get_untracked() != target_context
                || target.get_untracked() != target_context.revision()
                || epoch.get_untracked() != captured_epoch
            {
                return;
            }
            active.set_value(None);
            busy.set(false);
            match result {
                Ok(AnalysisResult::GoldenGatePrimersDesigned(report)) => {
                    match report.validate_binding(&revision) {
                        Ok(()) if report.source_revisions == revisions => {
                            gate_suggestions.set(Some(Arc::new(GoldenGateSuggestions {
                                report,
                                config: state,
                                context,
                                target_context,
                                epoch: captured_epoch,
                            })));
                            status.set("Type IIS primer proposals ready. Review tails and actual digested ends, then explicitly apply.".into());
                        }
                        Ok(()) => status
                            .set("Primer proposals belong to different original sources.".into()),
                        Err(error) => status.set(error.to_string()),
                    }
                }
                Ok(AnalysisResult::Error(error)) | Err(error) => status.set(error),
                _ => status.set("Unexpected Golden Gate primer response.".into()),
            }
        });
    });
    let apply_gate = Callback::new(move |()| {
        let Some(saved) = gate_suggestions.get_untracked() else {
            return;
        };
        if config.get_untracked() != saved.config
            || source_context.get_untracked() != Some(saved.context.clone())
            || target_info.get_untracked() != saved.target_context
            || epoch.get_untracked() != saved.epoch
        {
            status.set("Sources or settings changed; request fresh Type IIS primers.".into());
            return;
        }
        let mut next = saved.config.clone();
        for fragment in &saved.report.fragments {
            let Some(index) = (0..next.rows.len()).find(|index| {
                fragment.role == format!("fragment-{}", index + 1)
                    && fragment.input_id == format!("source-{}", index + 1)
            }) else {
                status.set("Proposal has an unknown source role; nothing applied.".into());
                return;
            };
            let collagene_core::assembly_recipe::Preparation::PcrCatalogDigest {
                forward,
                reverse,
                selection,
                fragment_index,
                repair,
            } = &fragment.preparation
            else {
                status.set("Proposal is not an explicit PCR/digest chain; nothing applied.".into());
                return;
            };
            let row = &mut next.rows[index];
            if row.preparation != OverlapPreparation::PcrCatalogDigest
                || row.reverse != fragment.reverse
            {
                status.set("Preparation or orientation changed; nothing applied.".into());
                return;
            }
            row.forward_name = forward.name.clone();
            row.forward = forward.annealing.as_str().into();
            row.forward_tail = forward.tail.clone();
            row.reverse_name = reverse.name.clone();
            row.reverse_primer = reverse.annealing.as_str().into();
            row.reverse_tail = reverse.tail.clone();
            row.catalog_selection = Some(selection.clone());
            row.fragment = (fragment_index + 1).to_string();
            row.repair = *repair;
        }
        config.set(next);
        gate_suggestions.set(None);
        evidence.set(None);
        status.set("Proposed PCR/digest preparations applied explicitly. Existing-site rows were preserved. Compute method evidence to check actual junctions and the intended target.".into());
    });
    let suggestion_alive = alive.clone();
    let suggest = Callback::new(move |()| {
        let state = config.get_untracked();
        let selected = planner_form.get_untracked();
        if state.mode != Mode::Overlap || disable_plain_primer_design {
            status.set("Choose overlap assembly before suggesting primers.".into());
            return;
        }
        let (Some(inputs), Some(context)) =
            (sources.get_untracked(), source_context.get_untracked())
        else {
            status.set("Choose available original source documents.".into());
            return;
        };
        let bound_target = molecule.get_untracked();
        let target_context = target_info.get_untracked();
        let bound_epoch = epoch.get_untracked();
        let (inputs, request) =
            match primer_design_request(&state, inputs, &bound_target, &selected) {
                Ok(v) => v,
                Err(e) => {
                    status.set(e);
                    return;
                }
            };
        let revisions = inputs
            .iter()
            .map(|p| p.revision.clone())
            .collect::<Vec<_>>();
        if revisions
            != context
                .iter()
                .map(|s| s.binding.revision.clone())
                .collect::<Vec<_>>()
        {
            status.set("Original sources changed before primer planning.".into());
            return;
        }
        let revision = match collagene_core::overlap_primer_design::design_revision(
            &inputs,
            &bound_target,
            &request,
        ) {
            Ok(v) => v,
            Err(e) => {
                status.set(e.to_string());
                return;
            }
        };
        let target_revision = bound_target.revision();
        if target_revision != target_context.revision() {
            status.set("Target changed before primer planning.".into());
            return;
        }
        if source_context.get_untracked() != Some(context.clone())
            || target_info.get_untracked() != target_context
            || config.get_untracked() != state
            || planner_form.get_untracked() != selected
            || epoch.get_untracked() != bound_epoch
        {
            status.set("Sources or settings changed before primer planning.".into());
            return;
        }
        let worker = match CachedWorker::new() {
            Ok(v) => Rc::new(v),
            Err(e) => {
                status.set(e);
                return;
            }
        };
        if let Some(old) = active.get_value() {
            old.close();
        }
        active.set_value(Some(worker.clone()));
        let token = generation.get_untracked().wrapping_add(1);
        generation.set(token);
        busy.set(true);
        evidence.set(None);
        suggestions.set(None);
        gate_suggestions.set(None);
        status.set("Suggesting deterministic annealing domains and homology tails from original source regions…".into());
        let alive = suggestion_alive.clone();
        spawn_local(async move {
            let response = ask(
                &worker,
                AnalysisJob::DesignOverlapPrimers {
                    inputs,
                    target: bound_target,
                    request,
                },
                revision.clone(),
                token,
            )
            .await;
            worker.close();
            if !alive.load(Ordering::Relaxed)
                || generation.get_untracked() != token
                || config.get_untracked() != state
                || planner_form.get_untracked() != selected
                || source_context.get_untracked() != Some(context.clone())
                || target_info.get_untracked() != target_context
                || target.get_untracked() != target_revision
                || epoch.get_untracked() != bound_epoch
            {
                return;
            }
            active.set_value(None);
            busy.set(false);
            match response {
                Ok(AnalysisResult::OverlapPrimersDesigned(report)) => match report
                    .validate_binding(&revision)
                {
                    Ok(()) if report.source_revisions == revisions => {
                        suggestions.set(Some(Arc::new(Suggestions {
                            report,
                            config: state,
                            form: selected,
                            context,
                            target: target_revision,
                            target_context,
                            epoch: bound_epoch,
                        })));
                        status.set("Primer proposals ready. Review ordered sequences and explicitly apply before computing assembly evidence.".into());
                    }
                    Ok(()) => status.set(
                        "Primer suggestions did not match the original source revisions.".into(),
                    ),
                    Err(error) => status.set(error.to_string()),
                },
                Ok(_) => status.set("Unexpected overlap primer design response.".into()),
                Err(error) => status.set(error),
            }
        });
    });
    let apply_suggestions = Callback::new(move |()| {
        let Some(proposal) = suggestions.get_untracked() else {
            return;
        };
        if busy.get_untracked()
            || config.get_untracked() != proposal.config
            || planner_form.get_untracked() != proposal.form
            || source_context.get_untracked() != Some(proposal.context.clone())
            || target_info.get_untracked() != proposal.target_context
            || target.get_untracked() != proposal.target
            || epoch.get_untracked() != proposal.epoch
        {
            status.set("Sources, target or settings changed; request fresh primer suggestions before applying.".into());
            return;
        }
        if proposal.report.recipe.fragments.len() != proposal.config.rows.len() {
            status.set("Suggested primer rows do not match the original fragment order.".into());
            return;
        }
        let mut updated = proposal.config.clone();
        for (i, fragment) in proposal.report.recipe.fragments.iter().enumerate() {
            if fragment.input_id != format!("source-{}", i + 1)
                || fragment.reverse != updated.rows[i].reverse
            {
                status
                    .set("Suggested source order or orientation changed; nothing applied.".into());
                return;
            }
            let collagene_core::assembly_recipe::Preparation::ExactPcr { forward, reverse } =
                &fragment.preparation
            else {
                status
                    .set("Suggestion is not an explicit PCR preparation; nothing applied.".into());
                return;
            };
            let row = &mut updated.rows[i];
            row.preparation = OverlapPreparation::ExactPcr;
            row.forward_name = forward.name.clone();
            row.forward = forward.annealing.as_str().into();
            row.forward_tail = forward.tail.clone();
            row.reverse_name = reverse.name.clone();
            row.reverse_primer = reverse.annealing.as_str().into();
            row.reverse_tail = reverse.tail.clone();
        }
        updated.overlaps = proposal
            .report
            .recipe
            .overlaps
            .iter()
            .map(usize::to_string)
            .collect::<Vec<_>>()
            .join(",");
        config.set(updated);
        suggestions.set(None);
        gate_suggestions.set(None);
        evidence.set(None);
        status.set("Suggested primers copied explicitly. Compute method evidence and verify the prepared assembly against the independent target.".into());
    });
    let compute = Callback::new(move |()| {
        let state = config.get_untracked();
        let Some(inputs) = sources.get_untracked() else {
            status.set("Choose available source documents".into());
            return;
        };
        let Some(bound_context) = source_context.get_untracked() else {
            status.set("Choose original source events before computing method evidence.".into());
            return;
        };
        let bound_target_context = target_info.get_untracked();
        if source_revisions.get_untracked()
            != Some(
                bound_context
                    .iter()
                    .map(|info| info.binding.revision.clone())
                    .collect(),
            )
            || target.get_untracked() != bound_target_context.revision()
        {
            status.set("Source or target event revisions differ from their current models. Select the current documents again.".into());
            return;
        }

        let recipes = if state.mode == Mode::Tailing {
            match state
                .vector
                .recipe()
                .and_then(|v| state.insert.recipe().map(|i| vec![v, i]))
            {
                Ok(r) => r,
                Err(e) => {
                    status.set(e);
                    return;
                }
            }
        } else {
            vec![]
        };
        let profile = if state.mode == Mode::Overlap {
            match state.profile() {
                Ok(p) => Some(p),
                Err(e) => {
                    status.set(e);
                    return;
                }
            }
        } else {
            None
        };
        let worker = match CachedWorker::new() {
            Ok(w) => Rc::new(w),
            Err(e) => {
                status.set(e);
                return;
            }
        };
        if let Some(old) = active.get_value() {
            old.close();
        }
        active.set_value(Some(worker.clone()));
        let token = generation.get_untracked().wrapping_add(1);
        generation.set(token);
        busy.set(true);
        evidence.set(None);
        status.set("Reconstructing method evidence…".into());
        let bound_target = target.get_untracked();
        let bound_target_molecule = molecule.get_untracked();
        let bound_epoch = epoch.get_untracked();
        let alive = alive.clone();
        spawn_local(async move {
            let revs = inputs.iter().map(Molecule::revision).collect::<Vec<_>>();
            let result = async {
                if matches!(state.mode, Mode::Ligation | Mode::GoldenGate) {
                    let (bound_inputs, recipe) =
                        ligation_request(&state, inputs, &bound_target_molecule)?;
                    let revision = collagene_core::ligation_recipe::preview_revision(
                        &bound_inputs,
                        &bound_target_molecule,
                        &recipe,
                    )
                    .map_err(|e| e.to_string())?;
                    let response = ask(
                        &worker,
                        AnalysisJob::PreviewLigationRecipe {
                            inputs: bound_inputs,
                            target: bound_target_molecule,
                            recipe,
                        },
                        revision.clone(),
                        token,
                    )
                    .await?;
                    if let AnalysisResult::LigationRecipePreview(preview) = response {
                        preview
                            .validate_binding(&revision)
                            .map_err(|e| e.to_string())?;
                        Ok((vec![], None, None, Some(preview)))
                    } else {
                        Err("Unexpected ordered ligation response".into())
                    }
                } else if let Some((overlaps, options)) = profile {
                    let _ = (overlaps, options); // shared recipe reuses the exact parsed profile.
                    let (bound_inputs, recipe) =
                        overlap_request(&state, inputs, &bound_target_molecule)?;
                    let revision = collagene_core::assembly_recipe::preview_revision(
                        &bound_inputs,
                        &bound_target_molecule,
                        &recipe,
                    )
                    .map_err(|e| e.to_string())?;
                    let response = ask(
                        &worker,
                        AnalysisJob::PreviewOverlapRecipe {
                            inputs: bound_inputs,
                            target: bound_target_molecule,
                            recipe,
                        },
                        revision.clone(),
                        token,
                    )
                    .await?;
                    if let AnalysisResult::OverlapRecipePreview(preview) = response {
                        preview
                            .validate_binding(&revision)
                            .map_err(|e| e.to_string())?;
                        Ok((vec![], preview.profile.clone(), Some(preview), None))
                    } else {
                        Err("Unexpected overlap preparation response".into())
                    }
                } else {
                    let mut previews = vec![];
                    for (i, (input, recipe)) in inputs.into_iter().zip(recipes).enumerate() {
                        let base = if i == 0 {
                            state.vector.base
                        } else {
                            state.insert.base
                        };
                        let result = ask(
                            &worker,
                            AnalysisJob::TailPreview {
                                molecule: input,
                                recipe,
                                base,
                            },
                            revs[i].clone(),
                            token.wrapping_add(i as u32),
                        )
                        .await
                        .map_err(|e| {
                            format!("{}: {e}", if i == 0 { "Vector" } else { "Insert" })
                        })?;
                        if let AnalysisResult::TailPreview(preview) = result {
                            previews.push(preview);
                        } else {
                            return Err("Unexpected tailing preview response".into());
                        }
                    }
                    Ok((previews, None, None, None))
                }
            }
            .await;
            worker.close();
            if !alive.load(Ordering::Relaxed)
                || generation.get_untracked() != token
                || config.get_untracked() != state
                || source_revisions.get_untracked() != Some(revs.clone())
                || source_context.get_untracked() != Some(bound_context.clone())
                || target_info.get_untracked() != bound_target_context
                || target.get_untracked() != bound_target
                || epoch.get_untracked() != bound_epoch
            {
                return;
            }
            active.set_value(None);
            busy.set(false);
            match result {
                Ok((tails, profile, overlap, ligation)) => {
                    evidence.set(Some(Arc::new(Evidence {
                        config: state,
                        revisions: revs,
                        context: bound_context,
                        target_context: bound_target_context,
                        target: bound_target,
                        epoch: bound_epoch,
                        tails,
                        profile,
                        overlap,
                        ligation,
                    })));
                    status.set("Evidence ready. Review assumptions and build a candidate against the independently selected target.".into());
                }
                Err(e) => status.set(e),
            }
        });
    });
    let diagnostic_guard = Callback::new(move |saved: Arc<Evidence>| -> Result<(), String> {
        if !evidence.with_untracked(|current| {
            current
                .as_ref()
                .is_some_and(|value| Arc::ptr_eq(value, &saved))
        }) || config.get_untracked() != saved.config
            || source_context.get_untracked() != Some(saved.context.clone())
            || target_info.get_untracked() != saved.target_context
            || source_revisions.get_untracked() != Some(saved.revisions.clone())
            || target.get_untracked() != saved.target
            || epoch.get_untracked() != saved.epoch
        {
            return Err("Sources or preparation settings changed. Compute fresh evidence before navigating its diagnostics.".into());
        }
        Ok(())
    });
    let review = Callback::new(move |()| {
        let (Some(proof), Some(inputs)) = (evidence.get_untracked(), sources.get_untracked())
        else {
            return;
        };
        if let Err(error) = diagnostic_guard.run(proof.clone()) {
            status.set(error);
            return;
        }
        let state = config.get_untracked();
        let current_target = molecule.get_untracked();
        if proof.config != state
            || proof.revisions != inputs.iter().map(Molecule::revision).collect::<Vec<_>>()
            || proof.target != current_target.revision()
            || proof.epoch != epoch.get_untracked()
        {
            status.set(
                "Inputs changed; recompute method evidence before reviewing a candidate".into(),
            );
            return;
        }
        let mut steps = vec![];
        let mut bound_inputs = vec![];
        if state.mode == Mode::Tailing {
            for (id, mol) in ["vector", "insert"].into_iter().zip(inputs) {
                bound_inputs.push(PlanInput {
                    id: id.into(),
                    revision: mol.revision(),
                    molecule: mol,
                });
            }
            let vector = match append_preparation(&mut steps, "vector", &state.vector) {
                Ok(s) => s,
                Err(e) => {
                    status.set(e);
                    return;
                }
            };
            let insert = match append_preparation(&mut steps, "insert", &state.insert) {
                Ok(s) => s,
                Err(e) => {
                    status.set(e);
                    return;
                }
            };
            steps.push(Operation::Ligate {
                fragments: vec![
                    OrientedSource {
                        source: vector,
                        reverse: state.vector.reversed,
                    },
                    OrientedSource {
                        source: insert,
                        reverse: state.insert.reversed,
                    },
                ],
                circular: state.circular,
                allow_nicked: state.allow_nicked,
            });
        } else if matches!(state.mode, Mode::Ligation | Mode::GoldenGate) {
            let (bound, recipe) = match ligation_request(&state, inputs, &current_target) {
                Ok(v) => v,
                Err(e) => {
                    status.set(e);
                    return;
                }
            };
            match collagene_core::ligation_recipe::compile(&bound, &current_target, &recipe) {
                Ok(compiled) => {
                    on_candidate.run(Candidate::new(compiled.plan));
                    return;
                }
                Err(error) => {
                    status.set(error.to_string());
                    return;
                }
            }
        } else {
            let (bound, recipe) = match overlap_request(&state, inputs, &current_target) {
                Ok(v) => v,
                Err(e) => {
                    status.set(e);
                    return;
                }
            };
            match collagene_core::assembly_recipe::compile(&bound, &current_target, &recipe) {
                Ok(compiled) => {
                    on_candidate.run(Candidate::new(compiled.plan));
                    return;
                }
                Err(error) => {
                    status.set(error.to_string());
                    return;
                }
            }
        }
        on_candidate.run(Candidate::new(DesignPlan {
            name: format!(
                "{} to {}",
                if state.mode == Mode::Tailing {
                    "Declared tailing and ligation"
                } else {
                    "Profiled overlap assembly"
                },
                current_target.name
            ),
            target_revision: current_target.revision(),
            inputs: bound_inputs,
            steps,
        }));
    });
    view! {<div class="inspector-content cloning-methods-panel"><h2>"Cloning methods"</h2><p class="helper">"The active document is the independently chosen target. Method evidence never selects a target or establishes experimental success."</p>
        <label>"Method workflow"<select aria-label="Cloning workflow" prop:value=move||match config.get().mode{Mode::Tailing=>"tailing",Mode::Overlap=>"overlap",Mode::Ligation=>"ligation",Mode::GoldenGate=>"golden_gate"} on:change=move|e|config.update(|c|c.mode=match event_target_value(&e).as_str(){"overlap"=>Mode::Overlap,"ligation"=>Mode::Ligation,"golden_gate"=>Mode::GoldenGate,_=>Mode::Tailing})><option value="tailing" disabled=fixed_sources>"TA / GC tailing and ligation"</option><option value="overlap">"Gibson / HiFi / In-Fusion overlap"</option><option value="ligation" disabled=disable_plain_primer_design>"Ordered multi-fragment ligation"</option><option value="golden_gate" disabled=disable_plain_primer_design>"Golden Gate — prepared ordered assembly"</option></select></label>
        <Show when=move||config.get().mode==Mode::Tailing><p class="helper">"Append one declared nucleotide to both actual 3′ ends of an established blunt duplex. This is an ideal endpoint operation; polymerase activity, mixtures and proprietary kit chemistry are not predicted. TA/GC products include two added junction bases."</p><div class="method-shortcuts"><button class="text-tool" on:click=move|_|config.update(|c|{c.vector.base=TailingBase::T;c.insert.base=TailingBase::A;})>"Use T/A pair"</button><button class="text-tool" on:click=move|_|config.update(|c|{c.vector.base=TailingBase::C;c.insert.base=TailingBase::G;})>"Use C/G pair"</button></div><TailFields vector=true config=config query_options=query_options primer_sources list_primers resolve_primer epoch/><TailFields vector=false config=config query_options=query_options primer_sources list_primers resolve_primer epoch/><label class="checkbox-label"><input type="checkbox" aria-label="Allow retained tailing ligation nicks" prop:checked=move||config.get().allow_nicked on:change=move|e|config.update(|c|c.allow_nicked=event_target_checked(&e))/>"Allow explicitly retained nicks"</label><p class="helper">"Exact PCR has 5′ OH termini. A sealed ligation may need an explicit phosphorylation step. Factory-vector chemistry and host nick repair are not inferred."</p></Show>
        <Show when=move||config.get().mode!=Mode::Tailing><Show when=move||config.get().mode==Mode::Overlap><label>"Overlap profile"<select aria-label="Overlap profile method" prop:value=move||format!("{:?}",config.get().method) on:change=move|e|config.update(|c|c.method=match event_target_value(&e).as_str(){"Gibson"=>OverlapMethod::Gibson,"InFusion"=>OverlapMethod::InFusion,_=>OverlapMethod::HiFi})><option value="Gibson">"Gibson"</option><option value="HiFi">"HiFi"</option><option value="InFusion">"In-Fusion"</option></select></label><p class="helper">"Fragments are listed in product order with explicit orientation. Homology geometry is checked exactly; published length guidance is advisory. No reaction completion or unique mixture product is established."</p></Show>
            <Show when=move||matches!(config.get().mode,Mode::Ligation|Mode::GoldenGate)><p class="helper">"List fragments in intended ligation order. Preparation, end repair and terminal treatment are explicit source-bound operations. Whole linear DNA retains unknown terminal chemistry; its sequence alone does not imply cleaved or PCR ends. Unknown ends can withhold ligation even after phosphorylation."</p><label class="checkbox-label"><input type="checkbox" aria-label="Allow retained ordered ligation nicks" prop:checked=move||config.get().allow_nicked on:change=move|e|config.update(|c|c.allow_nicked=event_target_checked(&e))/>"Allow explicitly retained nicks"</label></Show>
            <h3>"1. Choose source fragments and regions"</h3>
            <Show when=move||config.with(|c|c.rows.len()>12)><p role="alert">"This preparation exceeds the supported 12 rows. No partial row list is shown."</p></Show><For each=move||config.with(|c|if c.rows.len()<=12{(0..c.rows.len()).collect::<Vec<_>>()}else{Vec::new()}) key=|i|*i children=move|i|view!{<div node_ref=row_nodes.with_value(|nodes|nodes[i]) tabindex="-1" class="method-preparation-focus"><OverlapFields config index=i query_options primer_sources list_primers resolve_primer epoch selected_region fixed_sources/></div>}/>
            <button class="button secondary full" disabled=move||{fixed_sources||config.with(|c|c.rows.len()>=12)} on:click=move|_|config.update(|c|c.rows.push(OverlapInput::default()))>{move||if matches!(config.get().mode,Mode::Ligation|Mode::GoldenGate){"Add ligation fragment"}else{"Add overlap fragment"}}</button><p class="helper">"This panel supports up to 12 ordered fragments; the typed core has its separate documented work bounds."</p><Show when=move||config.get().mode==Mode::Overlap><label>"Overlap lengths in junction order (bp)"<input aria-label="Profile overlap lengths" prop:value=move||config.get().overlaps on:input=move|e|config.update(|c|c.overlaps=event_target_value(&e))/></label><label class="checkbox-label"><input type="checkbox" aria-label="Scan repeated overlap words" prop:checked=move||config.get().scan_repeats on:change=move|e|config.update(|c|c.scan_repeats=event_target_checked(&e))/>"Search supplied fragments for repeated homology words"</label><label>"Repeat comparison budget"<input type="number" min="1" aria-label="Overlap repeat budget" prop:value=move||config.get().comparisons on:input=move|e|config.update(|c|c.comparisons=event_target_value(&e))/></label><label>"Maximum reported homology occurrences"<input type="number" min="1" aria-label="Overlap occurrence limit" prop:value=move||config.get().occurrences on:input=move|e|config.update(|c|c.occurrences=event_target_value(&e))/></label></Show>
        </Show>
        <Show when=move||config.get().mode==Mode::GoldenGate><section class="golden-gate-primer-planner"><h3>"2. Prepare Type IIS fragments"</h3>
            <p class="helper">"Use existing pinned digest fragments, or choose Exact PCR → pinned catalog digest for rows needing new sites. Proposals use each chosen original source region and explicit right fusion. They do not replace existing-site rows or choose a target. This ordered reconstruction does not enumerate every mixture outcome."</p>
            <div class="primer-grid">
                <label>"Type IIS enzyme"<select aria-label="Golden Gate preparation enzyme" prop:value=move||format!("{:?}",config.with(|c|c.golden_gate.enzyme)) on:change=move|e|config.update(|c|c.golden_gate.enzyme=match event_target_value(&e).as_str(){"BsmBI"=>collagene_core::golden_gate::GoldenGateEnzyme::BsmBI,"BbsI"=>collagene_core::golden_gate::GoldenGateEnzyme::BbsI,"SapI"=>collagene_core::golden_gate::GoldenGateEnzyme::SapI,_=>collagene_core::golden_gate::GoldenGateEnzyme::BsaI})><option value="BsaI">"BsaI — 4-base fusion"</option><option value="BsmBI">"BsmBI — 4-base fusion"</option><option value="BbsI">"BbsI — 4-base fusion"</option><option value="SapI">"SapI — 3-base fusion"</option></select></label>
                <label>"Pre-recognition clamp"<input maxlength="12" aria-label="Golden Gate clamp" prop:value=move||config.with(|c|c.golden_gate.clamp.clone()) on:input=move|e|config.update(|c|c.golden_gate.clamp=event_target_value(&e))/></label>
                <label>"Cleavage spacer"<input maxlength="2" aria-label="Golden Gate spacer" prop:value=move||config.with(|c|c.golden_gate.spacer.clone()) on:input=move|e|config.update(|c|c.golden_gate.spacer=event_target_value(&e))/></label>
                <label>"Forward annealing length"<input type="number" min="1" max="200" aria-label="Golden Gate forward annealing length" prop:value=move||config.with(|c|c.golden_gate.forward_length.clone()) on:input=move|e|config.update(|c|c.golden_gate.forward_length=event_target_value(&e))/></label>
                <label>"Reverse annealing length"<input type="number" min="1" max="200" aria-label="Golden Gate reverse annealing length" prop:value=move||config.with(|c|c.golden_gate.reverse_length.clone()) on:input=move|e|config.update(|c|c.golden_gate.reverse_length=event_target_value(&e))/></label>
            </div>
            <p class="helper">"The clamp and spacer are explicit synthesized DNA. BbsI needs a 2-base spacer; the other listed enzymes need 1. These settings do not establish reagent compatibility, primer thermodynamics or amplification success."</p>
            <button class="button secondary full" disabled=move||busy.get() on:click=move|_|suggest_gate.run(())>"Suggest Type IIS PCR primers"</button>
            {move||gate_suggestions.get().map(|saved|view!{<section class="golden-gate-primer-proposals"><h4>"Proposed PCR/digest chains — not applied"</h4>
                {saved.report.rows.iter().map(|row|view!{<article><h4>{row.role.clone()}</h4>
                    <pre>{format!("{}\nOrdered 5′→3′: {}{}\nAnnealing: {}\nTail: clamp [{}] · motif [{}] · spacer [{}] · fusion [{}]\n\n{}\nOrdered 5′→3′: {}{}\nAnnealing: {}\nTail: clamp [{}] · motif [{}] · spacer [{}] · fusion [{}]",row.forward.name,row.forward.tail,row.forward.annealing.as_str(),row.forward.annealing.as_str(),row.forward_tail.clamp,row.forward_tail.motif,row.forward_tail.spacer,row.forward_tail.fusion,row.reverse.name,row.reverse.tail,row.reverse.annealing.as_str(),row.reverse.annealing.as_str(),row.reverse_tail.clamp,row.reverse_tail.motif,row.reverse_tail.spacer,row.reverse_tail.fusion)}</pre>
                    <p>{format!("Actual digested fragment {} · assembly-oriented top DNA: {}",row.fragment_index+1,short(&row.oriented_top))}</p><p>{format!("Left: {}",method_end(Some(&row.left)))}</p><p>{format!("Right: {}",method_end(Some(&row.right)))}</p>
                </article>}).collect_view()}
                {saved.report.warnings.iter().map(|warning|view!{<p class="helper">{warning.clone()}</p>}).collect_view()}
                <p class="helper">"Review each own-5′→3′ tail segment and the actual cleaved strand ends. Applying copies these named primers, the exact catalog pin and PCR-product fragment choice; preparation evidence must then be recomputed."</p>
                <button class="button primary full" disabled=move||busy.get() on:click=move|_|apply_gate.run(())>"Apply Type IIS primer proposals"</button>
            </section>})}
        </section></Show>
        <Show when=move||config.get().mode==Mode::Overlap&&!disable_plain_primer_design><section class="overlap-primer-planner"><h3>"2. Suggest and review PCR primers"</h3><p class="helper">"Use the selected source regions and order to suggest one primer pair per fragment, including added homology tails. Your current primers stay unchanged until you explicitly apply the proposal. The independent target is never generated or selected by this step."</p><div class="primer-grid"><label>"Junction homology (bp)"<input type="number" min="1" aria-label="Suggested overlap homology bases" prop:value=move||planner_form.with(|f|f.overlap.clone()) on:input=move|e|planner_form.update(|f|f.overlap=event_target_value(&e))/></label><label>"Minimum annealing length"<input type="number" min="1" max="200" aria-label="Suggested primer minimum length" prop:value=move||planner_form.with(|f|f.minimum.clone()) on:input=move|e|planner_form.update(|f|f.minimum=event_target_value(&e))/></label><label>"Maximum annealing length"<input type="number" min="1" max="200" aria-label="Suggested primer maximum length" prop:value=move||planner_form.with(|f|f.maximum.clone()) on:input=move|e|planner_form.update(|f|f.maximum=event_target_value(&e))/></label><label>"Target annealing Tm (°C)"<input type="number" step="0.1" aria-label="Suggested primer target Tm" prop:value=move||planner_form.with(|f|f.tm.clone()) on:input=move|e|planner_form.update(|f|f.tm=event_target_value(&e))/></label></div><p class="helper">"Default filters: GC 20–80%, pair ΔTm ≤5°C, homopolymer run ≤5 and 3′ complementary run ≤3. Tm assumes 50 mM sodium and 500 nM total strands; these are calculations, not reagent or experimental performance predictions."</p><button class="button secondary full" disabled=move||busy.get() on:click=move|_|suggest.run(())>"Suggest PCR primers"</button>
        {move||suggestions.get().map(|proposal|{let saved=proposal.clone();view!{<div class="overlap-primer-suggestions"><h4>"Proposed primer pairs — not yet applied"</h4>{proposal.report.designs.iter().enumerate().map(|(i,design)|{let interval=proposal.report.request.fragments.get(i);view!{<article><h4>{format!("Fragment {} · {}",i+1,proposal.context.get(i).map(|s|s.name.as_str()).unwrap_or("source unavailable"))}</h4>{interval.map(|region|view!{<p>{format!("Source interval {}–{} (1-based inclusive, unrolled) · {}",region.start+1,region.end,if region.reverse{"reverse in assembly"}else{"forward in assembly"})}</p>})}{design.pairs.first().map(|pair|view!{<pre>{format!("{}\nOrdered 5′→3′: {}{}\nAnnealing 5′→3′: {}\n5′ tail: {}\n\n{}\nOrdered 5′→3′: {}{}\nAnnealing 5′→3′: {}\n5′ tail: {}",pair.forward.primer.name,pair.forward.primer.tail,pair.forward.primer.annealing.as_str(),pair.forward.primer.annealing.as_str(),pair.forward.primer.tail,pair.reverse.primer.name,pair.reverse.primer.tail,pair.reverse.primer.annealing.as_str(),pair.reverse.primer.annealing.as_str(),pair.reverse.primer.tail)}</pre><p>{format!("Annealing Tm: {} / {} · pair ΔTm {:.1}°C · proposed PCR length {} bp",pair.forward.analysis.tm.as_ref().map(|t|format!("{:.1}°C",t.tm_c)).unwrap_or_else(||"unavailable".into()),pair.reverse.analysis.tm.as_ref().map(|t|format!("{:.1}°C",t.tm_c)).unwrap_or_else(||"unavailable".into()),pair.tm_difference_c,pair.product_length)}</p>})}{design.warnings.iter().map(|w|view!{<p>{w.clone()}</p>}).collect_view()}</article>}}).collect_view()}{proposal.report.warnings.iter().map(|w|view!{<p>{w.clone()}</p>}).collect_view()}<p class="helper">"A deterministic proposal under these filters is not a unique or experimentally validated primer design. Applying changes the preparation to exact PCR; verify the resulting assembly next."</p><button class="button primary full" disabled=move||busy.get() on:click=move|_|apply_suggestions.run(())>"Apply suggested primers"</button><button class="button secondary full" on:click=move|_|{if config.get_untracked()!=saved.config||planner_form.get_untracked()!=saved.form||source_context.get_untracked()!=Some(saved.context.clone())||target_info.get_untracked()!=saved.target_context||target.get_untracked()!=saved.target||epoch.get_untracked()!=saved.epoch{status.set("Sources or settings changed; request fresh suggestions.".into());return;}match serde_json::to_vec_pretty(&serde_json::json!({"sources":saved.context,"target":saved.target_context,"report":saved.report})){Ok(bytes)=>if let Err(e)=download("overlap-primer-suggestions.json","application/json",&bytes){status.set(e)},Err(e)=>status.set(e.to_string())}}>"Download suggested primer evidence"</button></div>}})}
        </section></Show><Show when=move||config.get().mode!=Mode::Tailing><h3>"3. Verify the prepared assembly"</h3><p class="helper">"Review the applied preparation and orientation, then compute method evidence against your independent target. Suggested primers alone do not verify a cloning design."</p></Show>
        <label class="checkbox-label"><input disabled=fixed_sources type="checkbox" aria-label="Method product circular" prop:checked=move||config.get().circular on:change=move|e|config.update(|c|c.circular=event_target_checked(&e))/>"Circular candidate product"</label><button class="button primary full" disabled=move||busy.get() on:click=move|_|compute.run(())>"Compute method evidence"</button><button class="text-tool" disabled=move||!busy.get() on:click=move|_|{generation.update(|n|*n=n.wrapping_add(1));if let Some(w)=active.get_value(){w.close();}active.set_value(None);busy.set(false);suggestions.set(None);
        gate_suggestions.set(None);status.set("Method analysis cancelled".into());}>"Cancel method analysis"</button><p class="method-status" role="status">{move||status.get()}</p>
        {move||evidence.get().map(|proof|{
            let report=proof.clone();
            let current=StoredValue::new(proof.clone());
            let edit=Callback::new(move|index:usize|{
                let outcome=current.with_value(|saved|{diagnostic_guard.run(saved.clone())?;
                    if index>=saved.config.rows.len(){return Err("The preparation row is unavailable.".into());}
                    let node=row_nodes.with_value(|nodes|nodes.get(index).copied()).and_then(|node|node.get()).ok_or("The preparation editor is not currently visible.")?;
                    node.scroll_into_view();node.focus().map_err(|_|"The preparation editor could not receive focus.".to_string())?;
                    Ok::<_,String>(())
                });if let Err(error)=outcome{status.set(error);}
            });
            let view_source=on_source.map(|navigate|Callback::new(move|row:usize|{
                let outcome=current.with_value(|saved|{diagnostic_guard.run(saved.clone())?;
                    let input=saved.config.rows.get(row).ok_or("The original source row is unavailable.")?;
                    navigate.run((input.document,row))
                });if let Err(error)=outcome{status.set(error);}
            }));
            let diagnostics=proof.overlap.as_ref().map(collagene_core::assembly_diagnostics::overlap).or_else(||proof.ligation.as_ref().map(collagene_core::assembly_diagnostics::ligation));
            let names=proof.context.iter().map(|source|source.name.clone()).collect::<Vec<_>>();
            let preparations=proof.config.rows.iter().map(OverlapInput::summary).collect::<Vec<_>>();
            view!{<div class="method-evidence">
            {diagnostics.map(|result|match result{Ok(diagnostics)=>view!{<MethodDiagnostics diagnostics names preparations on_edit=edit on_source=view_source/>}.into_any(),Err(error)=>view!{<p class="helper">{format!("Structured preparation diagnostics unavailable: {error}. The full fresh findings remain below.")}</p>}.into_any()})}
            {proof.ligation.as_ref().map(|p|(&p.compiled.plan,&p.verification)).or_else(||proof.overlap.as_ref().map(|p|(&p.compiled.plan,&p.verification))).map(|(plan,report)|{
                let captured = StoredValue::new(proof.clone());
                let captured_sources = StoredValue::new(Some(proof.context.clone()));
                let captured_target = StoredValue::new(proof.target_context.clone());
                let guarded_finalize = on_finalize.map(|finalize|Callback::new(move|candidate:Candidate| {
                    let fresh = captured.with_value(|saved| {
                        evidence.with_untracked(|current|current.as_ref().is_some_and(|current|Arc::ptr_eq(current,saved)))
                            && config.get_untracked()==saved.config
                            && source_revisions.get_untracked()==Some(saved.revisions.clone())
                            && target.get_untracked()==saved.target
                            && epoch.get_untracked()==saved.epoch
                    });
                    if !fresh || source_context.get_untracked()!=captured_sources.get_value() || target_info.get_untracked()!=captured_target.get_value() {
                        return Err("Sources or preparation settings changed. Compute fresh method evidence before creating outputs.".into());
                    }
                    finalize.run(candidate)
                }));
                view!{<crate::assembly_product_view::AssemblyProductView plan=plan.clone() target=molecule.get_untracked() report=report.clone() epoch on_products on_finalize=guarded_finalize/>}
            })}

            {proof.tails.iter().enumerate().map(|(i,preview)|{let opened=preview.molecule.clone();view!{<article class="tail-evidence"><h3>{format!("{} · {}-tailed duplex",if i==0{"Vector"}else{"Insert"},preview.base.as_str())}</h3><pre>{format!("Top 5′→3′: {}\nBottom 5′→3′: {}\n{}",short(preview.tailing.duplex.top.as_str()),short(preview.tailing.duplex.bottom.as_str()),preview.tailing.duplex.summary().map(|s|s.describe()).unwrap_or_else(|e|e.to_string()))}</pre><details><summary>"Source, additions and assumptions"</summary><pre>{format!("Source top 5′→3′: {}\nSource bottom 5′→3′: {}",short(preview.tailing.source.top.as_str()),short(preview.tailing.source.bottom.as_str()))}</pre>{preview.tailing.additions.iter().map(|a|view!{<p>{format!("New {:?} residue {}–{}: {} (no template ancestry)",a.strand,a.product_start+1,a.product_end,a.sequence)}</p>}).collect_view()}{preview.findings.iter().chain(&preview.tailing.warnings).map(|w|view!{<p>{w.clone()}</p>}).collect_view()}</details><button class="button secondary full" on:click=move|_|match on_products.run(vec![opened.clone()]){Ok(())=>status.set("Reference preview opened; reconstruct its original graph for later physical operations.".into()),Err(e)=>status.set(e)}>{format!("Open {} reference preview",if i==0{"vector"}else{"insert"})}</button></article>}}).collect_view()}
            {proof.ligation.as_ref().map(|preview|view!{<article class="ordered-ligation-evidence"><h3>"Reconstructed ligation evidence"</h3><p>{format!("Preparation reconstructed: {} · independent-target assessment: {}",preview.preparation_complete,method_status(&preview.verification.status))}</p><p class="helper">"These ends and junctions were freshly reconstructed from the original source operations. An end-compatible ordered candidate is not a unique mixture product or demonstrated reaction outcome."</p>{preview.prepared.iter().enumerate().map(|(i,fragment)|view!{<details><summary>{format!("Fragment {} · {} bp · {}",i+1,fragment.sequence.len(),proof.config.rows.get(i).map(OverlapInput::summary).unwrap_or_default())}</summary><pre>{format!("Prepared source before orientation: {}",short(fragment.sequence.as_str()))}</pre><p>{preview.physical.get(i).and_then(Option::as_ref).map(method_duplex).unwrap_or_else(||"Physical end geometry unavailable; no chemistry inferred".into())}</p><p class="helper">"End summary is in assembly orientation, after any explicit terminal treatment."</p></details>}).collect_view()}{preview.junctions.iter().map(|junction|view!{<details><summary>{format!("Junction {} · fragment {} → {}{}",junction.index+1,junction.upstream+1,junction.downstream+1,if junction.closing{" · circle closure"}else{""})}</summary><p>{format!("Upstream right: {}",method_end(junction.upstream_end.as_ref()))}</p><p>{format!("Downstream left: {}",method_end(junction.downstream_end.as_ref()))}</p></details>}).collect_view()}{preview.product_physical.as_ref().map(|physical|view!{<section><h4>"Reconstructed product"</h4><p>{method_duplex(physical)}</p>{physical.nicks.iter().map(|nick|view!{<p>{format!("Retained {} strand nick at top-reference boundary {}",match nick.strand{collagene_core::chemistry::Strand::Top=>"top",collagene_core::chemistry::Strand::Bottom=>"bottom"},nick.position)}</p>}).collect_view()}</section>})}{preview.verification.findings.iter().map(|finding|view!{<p>{finding.clone()}</p>}).collect_view()}</article>})}
            {proof.overlap.as_ref().map(|preview|view!{<article class="overlap-preparation-evidence"><h3>"Prepared fragments from original sources"</h3><p>{format!("Preparation reconstructed: {} · independent-target assessment: {}",preview.preparation_complete,method_status(&preview.verification.status))}</p><p class="helper">"Prepared reference sequences below are reconstructed from the declared originals. They are not substituted as plan inputs, confirmed amplicons or physical reaction outcomes."</p>{preview.prepared.iter().enumerate().map(|(i,fragment)|view!{<details><summary>{format!("Fragment {} · {} bp · {}",i+1,fragment.sequence.len(),proof.config.rows.get(i).map(OverlapInput::summary).unwrap_or_default())}</summary><pre>{format!("Prepared source 5′→3′ (before orientation): {}",short(fragment.sequence.as_str()))}</pre></details>}).collect_view()}{preview.verification.findings.iter().map(|finding|view!{<p>{finding.clone()}</p>}).collect_view()}</article>})}
            {proof.profile.as_ref().map(|p|view!{<article class="overlap-profile-evidence"><h3>"Overlap evidence"</h3><p>{format!("Geometry consistent: {} · analysis complete: {} · published length guidance: {}",p.geometry_valid,p.complete,match p.recommendation_compliant{Some(true)=>"met",Some(false)=>"outside published guidance",None=>"not established"})}</p><p class="helper">{format!("{} comparisons; no exhaustive alternative-product search",p.comparisons)}</p>{p.junctions.iter().map(|j|view!{<details class="profile-junction"><summary>{format!("Junction {} · {}bp · {}",j.index+1,j.overlap_bases,if j.exact_match{"exact end homology"}else{"incompatible homology"})}</summary><pre>{format!("Upstream {}\nDownstream {}",j.homology,j.downstream_homology)}</pre><p>{format!("GC: {} · advisory length: {} · repeat search complete: {}",j.gc_fraction.map(|f|format!("{:.1}%",f*100.)).unwrap_or_else(||"unavailable".into()),guidance(j.recommended_min,j.recommended_max),j.repeat_search_complete)}</p>{j.occurrences.iter().take(40).map(|o|view!{<p>{format!("Source {} · {}–{}{}{}",o.fragment+1,o.start+1,o.end,if o.reverse{" reverse"}else{" forward"},if o.assigned_junction{" · assigned end"}else{" · additional occurrence"})}</p>}).collect_view()}<p class="helper">"Up to 40 occurrence positions shown; full evidence download retains every reported position."</p></details>}).collect_view()}<details><summary>"Profile findings and assumptions"</summary>{p.findings.iter().chain(&p.warnings).map(|w|view!{<p>{w.clone()}</p>}).collect_view()}</details></article>})}
            <button class="button primary full" on:click=move|_|review.run(())>"Review method candidate"</button><button class="button secondary full" on:click=move|_|{let result=if let Some(preview)=&report.ligation{serde_json::to_vec_pretty(preview)}else if let Some(preview)=&report.overlap{serde_json::to_vec_pretty(preview)}else if let Some(profile)=&report.profile{serde_json::to_vec_pretty(profile)}else{serde_json::to_vec_pretty(&report.tails)};match result{Ok(bytes)=>if let Err(e)=download("cloning-method-evidence.json","application/json",&bytes){status.set(e)},Err(e)=>status.set(e.to_string())}}>"Download method evidence"</button>
        </div>}})}
    </div>}
}

#[component]
fn TailFields(
    vector: bool,
    config: RwSignal<MethodsState>,
    query_options: Memo<Vec<(usize, String)>>,
    primer_sources: Memo<Vec<crate::digest::SourceInfo>>,
    list_primers: Callback<
        collagene_core::run_lineage::EventBinding,
        Result<collagene_core::primer_reuse::List, String>,
    >,
    resolve_primer: Callback<
        crate::primer_picker::ResolveRequest,
        Result<collagene_core::primer_reuse::Resolved, String>,
    >,
    epoch: RwSignal<u64>,
) -> impl IntoView {
    let label = if vector { "Vector" } else { "Insert" };
    view! {<fieldset class="tail-source-fields"><legend>{label}</legend><label>"Source document"<select aria-label=format!("{label} tailing source") prop:value=move||config.with(|c|c.input(vector).document.to_string()) on:change=move|e|config.update(|c|c.input_mut(vector).document=event_target_value(&e).parse().unwrap_or(0))>{move||query_options.get().into_iter().map(|(i,name)|view!{<option value=i.to_string()>{name}</option>}).collect_view()}</select></label><label>"Preparation"<select aria-label=format!("{label} tailing preparation") prop:value=move||if config.with(|c|c.input(vector).pcr){"pcr"}else{"digest"} on:change=move|e|config.update(|c|c.input_mut(vector).pcr=event_target_value(&e)=="pcr")><option value="digest">"Restriction fragment"</option><option value="pcr">"Unique exact PCR product"</option></select></label>
        <Show when=move||!config.with(|c|c.input(vector).pcr)><label>"Restriction enzymes"<input aria-label=format!("{label} tailing enzymes") prop:value=move||config.with(|c|c.input(vector).enzymes.clone()) on:input=move|e|config.update(|c|c.input_mut(vector).enzymes=event_target_value(&e))/></label><label>"Digest fragment (cut order, 1-based)"<input type="number" min="1" aria-label=format!("{label} tailing fragment") prop:value=move||config.with(|c|c.input(vector).fragment.clone()) on:input=move|e|config.update(|c|c.input_mut(vector).fragment=event_target_value(&e))/></label><label class="checkbox-label"><input type="checkbox" aria-label=format!("{label} blunt before tailing") prop:checked=move||config.with(|c|c.input(vector).blunt) on:change=move|e|config.update(|c|c.input_mut(vector).blunt=event_target_checked(&e))/>"Blunt staggered ends before tailing"</label></Show>
        <Show when=move||config.with(|c|c.input(vector).pcr)>
        <crate::primer_picker::PrimerPicker label=format!("{label} tailing forward") sources=primer_sources list_primers resolve_primer epoch on_copy=Callback::new(move|primer:Primer|config.update(|c|{let r=c.input_mut(vector);r.forward_name=primer.name;r.forward=primer.annealing.as_str().into();r.forward_tail=primer.tail;}))/>
        <label>"Forward primer name"<input maxlength="256" aria-label=format!("{label} tailing forward primer name") prop:value=move||config.with(|c|c.input(vector).forward_name.clone()) on:input=move|e|config.update(|c|c.input_mut(vector).forward_name=event_target_value(&e))/></label>
        <crate::primer_picker::PrimerPicker label=format!("{label} tailing reverse") sources=primer_sources list_primers resolve_primer epoch on_copy=Callback::new(move|primer:Primer|config.update(|c|{let r=c.input_mut(vector);r.reverse_name=primer.name;r.reverse=primer.annealing.as_str().into();r.reverse_tail=primer.tail;}))/>
        <label>"Reverse primer name"<input maxlength="256" aria-label=format!("{label} tailing reverse primer name") prop:value=move||config.with(|c|c.input(vector).reverse_name.clone()) on:input=move|e|config.update(|c|c.input_mut(vector).reverse_name=event_target_value(&e))/></label>
    <label>"Forward annealing sequence 5′→3′"<input aria-label=format!("{label} tailing forward primer") prop:value=move||config.with(|c|c.input(vector).forward.clone()) on:input=move|e|config.update(|c|c.input_mut(vector).forward=event_target_value(&e))/></label><label>"Reverse annealing sequence 5′→3′"<input aria-label=format!("{label} tailing reverse primer") prop:value=move||config.with(|c|c.input(vector).reverse.clone()) on:input=move|e|config.update(|c|c.input_mut(vector).reverse=event_target_value(&e))/></label><label>"Forward primer 5′ extension"<input aria-label=format!("{label} tailing forward extension") prop:value=move||config.with(|c|c.input(vector).forward_tail.clone()) on:input=move|e|config.update(|c|c.input_mut(vector).forward_tail=event_target_value(&e))/></label><label>"Reverse primer 5′ extension"<input aria-label=format!("{label} tailing reverse extension") prop:value=move||config.with(|c|c.input(vector).reverse_tail.clone()) on:input=move|e|config.update(|c|c.input_mut(vector).reverse_tail=event_target_value(&e))/></label><details><summary>"Review current primer domains"</summary><pre>{move||config.with(|c|{let r=c.input(vector);format!("{}\nOrdered 5′→3′: {}{}\nAnnealing: {}\n5′ tail: {}\n\n{}\nOrdered 5′→3′: {}{}\nAnnealing: {}\n5′ tail: {}",r.forward_name,r.forward_tail,r.forward,r.forward,r.forward_tail,r.reverse_name,r.reverse_tail,r.reverse,r.reverse,r.reverse_tail)})}</pre></details></Show>
        <label>"One nucleotide added to BOTH 3′ ends"<select aria-label=format!("{label} tailing base") prop:value=move||config.with(|c|c.input(vector).base.as_str()) on:change=move|e|config.update(|c|c.input_mut(vector).base=match event_target_value(&e).as_str(){"C"=>TailingBase::C,"G"=>TailingBase::G,"T"=>TailingBase::T,_=>TailingBase::A})><option value="A">"A"</option><option value="C">"C"</option><option value="G">"G"</option><option value="T">"T"</option></select></label><label>"Explicit terminal treatment after tailing"<select aria-label=format!("{label} tailing terminal treatment") prop:value=move||config.with(|c|match c.input(vector).treatment{Treatment::Keep=>"keep",Treatment::Phosphorylate=>"phosphorylate",Treatment::Dephosphorylate=>"dephosphorylate"}) on:change=move|e|config.update(|c|c.input_mut(vector).treatment=match event_target_value(&e).as_str(){"phosphorylate"=>Treatment::Phosphorylate,"dephosphorylate"=>Treatment::Dephosphorylate,_=>Treatment::Keep})><option value="keep">"Keep reconstructed 5′ chemistry"</option><option value="phosphorylate">"Add 5′ phosphorylation"</option><option value="dephosphorylate">"Add 5′ dephosphorylation"</option></select></label><label class="checkbox-label"><input type="checkbox" aria-label=format!("Reverse {label} tailing product") prop:checked=move||config.with(|c|c.input(vector).reversed) on:change=move|e|config.update(|c|c.input_mut(vector).reversed=event_target_checked(&e))/>"Reverse fragment for ligation"</label></fieldset>}
}

#[component]
fn OverlapFields(
    config: RwSignal<MethodsState>,
    index: usize,
    query_options: Memo<Vec<(usize, String)>>,
    primer_sources: Memo<Vec<crate::digest::SourceInfo>>,
    list_primers: Callback<
        collagene_core::run_lineage::EventBinding,
        Result<collagene_core::primer_reuse::List, String>,
    >,
    resolve_primer: Callback<
        crate::primer_picker::ResolveRequest,
        Result<collagene_core::primer_reuse::Resolved, String>,
    >,
    epoch: RwSignal<u64>,
    selected_region: Callback<(), Option<(usize, usize, usize)>>,
    fixed_sources: bool,
) -> impl IntoView {
    let row = move || config.with(|c| c.rows.get(index).cloned().unwrap_or_default());
    let region_status = RwSignal::new(String::new());
    let catalog = StoredValue::new(collagene_core::catalog_digest::builtin_catalog().ok());

    view! {<fieldset class="overlap-source-row"><legend>{format!("Fragment {}",index+1)}</legend>
        <label>"Original source"<select disabled=fixed_sources aria-label=move||format!("{} source {}",if matches!(config.get().mode,Mode::Ligation|Mode::GoldenGate){"Ligation"}else{"Overlap"},index+1) prop:value=move||row().document.to_string() on:change=move|e|config.update(|c|if let Some(r)=c.rows.get_mut(index){r.document=event_target_value(&e).parse().unwrap_or(0)})>{move||query_options.get().into_iter().map(|(n,label)|view!{<option value=n.to_string()>{label}</option>}).collect_view()}</select></label>
        <Show when=move||matches!(config.get().mode,Mode::Overlap|Mode::GoldenGate)><button class="button secondary" disabled=fixed_sources aria-label=format!("Set fragment {} from sequence selection",index+1) on:click=move|_|{match selected_region.run(()){Some((document,start,end)) if start<end&&start.checked_add(1).is_some()=>{config.update(|c|if let Some(row)=c.rows.get_mut(index){row.document=document;row.interval_start=(start+1).to_string();row.interval_end=end.to_string();});region_status.set("Selected source region copied. Review preparation and orientation; primers have not been changed.".into());},_=>region_status.set("Select one contiguous region in the sequence view first. Compound or origin-spanning feature selections need an explicit interval.".into())}} >"Set from sequence selection"</button><p class="helper" role="status">{move||region_status.get()}</p><label>"Primer design interval start (1-based inclusive)"<input disabled=fixed_sources type="number" min="1" aria-label=format!("Overlap interval start {}",index+1) prop:value=move||row().interval_start on:input=move|e|config.update(|c|if let Some(r)=c.rows.get_mut(index){r.interval_start=event_target_value(&e)})/></label><label>"Primer design interval end (1-based inclusive; blank = source end)"<input disabled=fixed_sources type="number" min="1" aria-label=format!("Overlap interval end {}",index+1) prop:value=move||row().interval_end on:input=move|e|config.update(|c|if let Some(r)=c.rows.get_mut(index){r.interval_end=event_target_value(&e)})/></label><p class="helper">"Intervals declare intended regions for primer proposals. For a circular source, an end beyond its length describes one unrolled origin crossing. Existing whole-source, digest and manual PCR preparations are not silently sliced."</p></Show>
        <Show when=move||config.get().mode==Mode::GoldenGate><label>"Right fusion end (assembly direction, 5′→3′)"<input maxlength="4" aria-label=format!("Golden Gate right fusion {}",index+1) prop:value=move||row().right_fusion on:input=move|e|config.update(|c|if let Some(row)=c.rows.get_mut(index){row.right_fusion=event_target_value(&e)})/></label><p class="helper">{move||format!("Enter {} concrete bases explicitly. The left fusion comes from this oriented source region; enter the next fragment’s first fusion bases explicitly for a matching junction. The right fusion is added design DNA. Existing-site compatibility is checked from its actual ends.",if config.with(|c|c.golden_gate.enzyme)==collagene_core::golden_gate::GoldenGateEnzyme::SapI{3}else{4})}</p></Show>
        <label>"Prepare fragment"<select aria-label=move||format!("{} preparation {}",if matches!(config.get().mode,Mode::Ligation|Mode::GoldenGate){"Ligation"}else{"Overlap"},index+1) prop:value=move||match row().preparation{OverlapPreparation::WholeLinear=>"whole",OverlapPreparation::DeclaredDuplex=>"declared",OverlapPreparation::Digest=>"digest",OverlapPreparation::CatalogDigest=>"catalog",OverlapPreparation::ExactPcr=>"pcr",OverlapPreparation::PcrCatalogDigest=>"pcr_catalog"} on:change=move|e|config.update(|c|if let Some(r)=c.rows.get_mut(index){r.preparation=match event_target_value(&e).as_str(){"declared"=>OverlapPreparation::DeclaredDuplex,"digest"=>OverlapPreparation::Digest,"catalog"=>OverlapPreparation::CatalogDigest,"pcr"=>OverlapPreparation::ExactPcr,"pcr_catalog"=>OverlapPreparation::PcrCatalogDigest,_=>OverlapPreparation::WholeLinear}})><option value="whole">"Whole linear DNA"</option><option value="declared">"Declare intact linear duplex"</option><option value="digest">"Manual starter-enzyme digest"</option><option value="catalog">"Pinned catalog digest fragment"</option><option value="pcr">"Exact PCR product"</option><option value="pcr_catalog">"Exact PCR → pinned catalog digest"</option></select></label>
        <Show when=move||row().preparation==OverlapPreparation::WholeLinear><p class="helper">"Use an already linear source. Circular DNA must first be prepared by restriction digest or exact PCR; its origin is not an assumed cut."</p></Show>
        <Show when=move||row().preparation==OverlapPreparation::DeclaredDuplex><p class="helper">"Operator assumption for an original, intact linear duplex: declare its two terminal groups explicitly. This does not measure chemistry, cut or repair DNA, change its reference sequence, or reuse producer metadata. Circular, single-stranded and imported physical/reference-only intermediate sources are rejected by the core."</p><DeclaredEndFields config index left=true/><DeclaredEndFields config index left=false/><p class="helper">"All groups initially remain unknown. Orientation is applied after this declaration; any ligation terminal treatment remains a separate explicit step."</p></Show>
        <Show when=move||matches!(row().preparation,OverlapPreparation::Digest|OverlapPreparation::CatalogDigest|OverlapPreparation::PcrCatalogDigest)><Show when=move||row().preparation==OverlapPreparation::Digest>
            <label>"Restriction enzymes (comma-separated)"<input aria-label=move||format!("{} digest enzymes {}",if matches!(config.get().mode,Mode::Ligation|Mode::GoldenGate){"Ligation"}else{"Overlap"},index+1) prop:value=move||row().enzymes on:input=move|e|config.update(|c|if let Some(r)=c.rows.get_mut(index){r.enzymes=event_target_value(&e)})/></label>
            <label>"Digest fragment number (1-based)"<input type="number" min="1" aria-label=move||format!("{} digest fragment {}",if matches!(config.get().mode,Mode::Ligation|Mode::GoldenGate){"Ligation"}else{"Overlap"},index+1) prop:value=move||row().fragment on:input=move|e|config.update(|c|if let Some(r)=c.rows.get_mut(index){r.fragment=event_target_value(&e)})/></label>
            </Show><Show when=move||row().preparation==OverlapPreparation::PcrCatalogDigest><label>"Enzyme for PCR product digestion"<select aria-label=format!("PCR catalog enzyme {}",index+1) prop:value=move||row().catalog_selection.as_ref().and_then(|s|s.enzymes.first().map(|e|e.id.clone())).unwrap_or_default() on:change=move|e|{let id=event_target_value(&e);if let Some(result)=catalog.with_value(|catalog|catalog.as_ref().map(|c|collagene_core::catalog_digest::select(c,&[id]))){match result{Ok(selection)=>config.update(|c|if let Some(row)=c.rows.get_mut(index){row.catalog_selection=Some(selection);}),Err(error)=>region_status.set(error.to_string())}}}><option value="">"Choose a pinned catalog enzyme…"</option>{catalog.with_value(|catalog|catalog.as_ref().map(|c|c.edition().enzymes.iter().map(|e|view!{<option value=e.id.clone()>{e.name.clone()}</option>}).collect_view()))}</select></label><label>"Digest fragment of PCR product (1-based)"<input type="number" min="1" aria-label=format!("PCR catalog fragment {}",index+1) prop:value=move||row().fragment on:input=move|e|config.update(|c|if let Some(row)=c.rows.get_mut(index){row.fragment=event_target_value(&e)})/></label><p class="helper">"Digest coordinates refer to the newly reconstructed PCR product, not the original source. Applied Type IIS proposals pin this selection and fragment automatically; manual edits require fresh reconstruction."</p></Show><Show when=move||matches!(row().preparation,OverlapPreparation::CatalogDigest|OverlapPreparation::PcrCatalogDigest)>{move||row().catalog_selection.map(|selection|view!{<p>{format!("Selected fragment {} · {}",row().fragment,selection.enzymes.iter().map(|e|e.enzyme.name.as_str()).collect::<Vec<_>>().join(", "))}</p><details><summary>"Pinned catalog identity"</summary><p>{format!("Edition {} · content {}",selection.catalog_edition,selection.catalog_content_hash)}</p>{selection.enzymes.iter().map(|entry|view!{<p>{format!("{} · {} · {} · cuts {}/{}",entry.id,entry.enzyme.name,entry.enzyme.recognition,entry.enzyme.top_cut,entry.enzyme.bottom_cut)}</p>}).collect_view()}</details>})}<p class="helper">"The exact catalog selection and fragment index are retained without starter-enzyme name substitution. For an original-source fragment, choose it in Virtual digest. PCR-product fragments use the explicit PCR/digest settings here. Saved pins describe a request; they do not restore current digest status or prove cleavage."</p></Show>
            <label>"Optional explicit end repair"<select aria-label=move||format!("{} end repair {}",if matches!(config.get().mode,Mode::Ligation|Mode::GoldenGate){"Ligation"}else{"Overlap"},index+1) prop:value=move||match row().repair{None=>"none",Some(EndRepairMethod::FillFivePrimeOverhangs)=>"fill",Some(EndRepairMethod::RemoveThreePrimeOverhangs)=>"remove",Some(EndRepairMethod::Blunt)=>"blunt"} on:change=move|e|config.update(|c|if let Some(r)=c.rows.get_mut(index){r.repair=match event_target_value(&e).as_str(){"fill"=>Some(EndRepairMethod::FillFivePrimeOverhangs),"remove"=>Some(EndRepairMethod::RemoveThreePrimeOverhangs),"blunt"=>Some(EndRepairMethod::Blunt),_=>None}})><option value="none">"None — retain cohesive ends"</option><option value="fill">"Fill 5′ overhangs"</option><option value="remove">"Remove 3′ overhangs"</option><option value="blunt">"Blunt ends"</option></select></label>
            <p class="helper">"Fragments are reconstructed in core digest order using the stated complete-cleavage assumption. This is not evidence of enzyme activity, methylation permission or reagent compatibility."</p>
        </Show>
        <Show when=move||matches!(row().preparation,OverlapPreparation::ExactPcr|OverlapPreparation::PcrCatalogDigest)>
            <crate::primer_picker::PrimerPicker label=format!("Fragment {} forward",index+1) sources=primer_sources list_primers resolve_primer epoch on_copy=Callback::new(move|primer:Primer|config.update(|c|if let Some(r)=c.rows.get_mut(index){r.forward_name=primer.name;r.forward=primer.annealing.as_str().into();r.forward_tail=primer.tail;}))/>
            <label>"Forward primer name"<input maxlength="256" aria-label=move||format!("{} forward primer name {}",if matches!(config.get().mode,Mode::Ligation|Mode::GoldenGate){"Ligation"}else{"Overlap"},index+1) prop:value=move||row().forward_name on:input=move|e|config.update(|c|if let Some(r)=c.rows.get_mut(index){r.forward_name=event_target_value(&e)})/></label>
            <crate::primer_picker::PrimerPicker label=format!("Fragment {} reverse",index+1) sources=primer_sources list_primers resolve_primer epoch on_copy=Callback::new(move|primer:Primer|config.update(|c|if let Some(r)=c.rows.get_mut(index){r.reverse_name=primer.name;r.reverse_primer=primer.annealing.as_str().into();r.reverse_tail=primer.tail;}))/>
            <label>"Reverse primer name"<input maxlength="256" aria-label=move||format!("{} reverse primer name {}",if matches!(config.get().mode,Mode::Ligation|Mode::GoldenGate){"Ligation"}else{"Overlap"},index+1) prop:value=move||row().reverse_name on:input=move|e|config.update(|c|if let Some(r)=c.rows.get_mut(index){r.reverse_name=event_target_value(&e)})/></label>

            <label>"Forward annealing sequence 5′→3′"<input aria-label=move||format!("{} forward primer {}",if matches!(config.get().mode,Mode::Ligation|Mode::GoldenGate){"Ligation"}else{"Overlap"},index+1) prop:value=move||row().forward on:input=move|e|config.update(|c|if let Some(r)=c.rows.get_mut(index){r.forward=event_target_value(&e)})/></label>
            <label>"Forward primer 5′ extension"<input aria-label=move||format!("{} forward extension {}",if matches!(config.get().mode,Mode::Ligation|Mode::GoldenGate){"Ligation"}else{"Overlap"},index+1) prop:value=move||row().forward_tail on:input=move|e|config.update(|c|if let Some(r)=c.rows.get_mut(index){r.forward_tail=event_target_value(&e)})/></label>
            <label>"Reverse annealing sequence 5′→3′"<input aria-label=move||format!("{} reverse primer {}",if matches!(config.get().mode,Mode::Ligation|Mode::GoldenGate){"Ligation"}else{"Overlap"},index+1) prop:value=move||row().reverse_primer on:input=move|e|config.update(|c|if let Some(r)=c.rows.get_mut(index){r.reverse_primer=event_target_value(&e)})/></label>
            <label>"Reverse primer 5′ extension"<input aria-label=move||format!("{} reverse extension {}",if matches!(config.get().mode,Mode::Ligation|Mode::GoldenGate){"Ligation"}else{"Overlap"},index+1) prop:value=move||row().reverse_tail on:input=move|e|config.update(|c|if let Some(r)=c.rows.get_mut(index){r.reverse_tail=event_target_value(&e)})/></label>
            <details><summary>"Review current primer domains"</summary><pre>{move||{let r=row();format!("{}\nOrdered 5′→3′: {}{}\nAnnealing: {}\n5′ tail: {}\n\n{}\nOrdered 5′→3′: {}{}\nAnnealing: {}\n5′ tail: {}",r.forward_name,r.forward_tail,r.forward,r.forward,r.forward_tail,r.reverse_name,r.reverse_tail,r.reverse_primer,r.reverse_primer,r.reverse_tail)}}</pre></details>
            <p class="helper">"Enter annealing domains and added tails separately. Both oligos are written 5′→3′; no primer design or unique amplicon is assumed."</p>
        </Show>
        <label class="checkbox-label"><input disabled=fixed_sources type="checkbox" aria-label=format!("Reverse overlap source {}",index+1) prop:checked=move||row().reverse on:change=move|e|config.update(|c|if let Some(r)=c.rows.get_mut(index){r.reverse=event_target_checked(&e)})/>"Reverse complement prepared fragment"</label>
        <Show when=move||matches!(config.get().mode,Mode::Ligation|Mode::GoldenGate)><label>"Explicit 5′ terminal treatment"<select aria-label=format!("Ligation terminal treatment {}",index+1) prop:value=move||match row().treatment{Treatment::Keep=>"none",Treatment::Phosphorylate=>"phosphorylate",Treatment::Dephosphorylate=>"dephosphorylate"} on:change=move|e|config.update(|c|if let Some(r)=c.rows.get_mut(index){r.treatment=match event_target_value(&e).as_str(){"phosphorylate"=>Treatment::Phosphorylate,"dephosphorylate"=>Treatment::Dephosphorylate,_=>Treatment::Keep}})><option value="none">"None — keep reconstructed chemistry"</option><option value="phosphorylate">"Add 5′ phosphorylation"</option><option value="dephosphorylate">"Add 5′ dephosphorylation"</option></select></label><p class="helper">"Treatment does not resolve unknown 3′ chemistry, repair ends or imply host nick repair."</p></Show>
        <p class="helper">{move||row().summary()}</p>
        <div class="method-row-actions"><button class="text-tool" disabled=move||fixed_sources||index==0 on:click=move|_|config.update(|c|{if index>0&&index<c.rows.len(){c.rows.swap(index,index-1);}})>"Move up"</button><button class="text-tool" disabled=move||fixed_sources||config.with(|c|index+1>=c.rows.len()) on:click=move|_|config.update(|c|{if index+1<c.rows.len(){c.rows.swap(index,index+1);}})>"Move down"</button><button class="text-tool" disabled=move||fixed_sources||config.with(|c|c.rows.len()<=1) on:click=move|_|config.update(|c|{if c.rows.len()>1&&index<c.rows.len(){c.rows.remove(index);}})>"Remove fragment"</button></div>
    </fieldset>}
}

pub(crate) fn prepared_rows(
    state: &MethodsState,
    input_count: usize,
) -> Result<
    (
        Vec<String>,
        Vec<collagene_core::assembly_recipe::FragmentSpec>,
    ),
    String,
> {
    use collagene_core::assembly_recipe::{FragmentSpec, Preparation};
    if input_count != state.rows.len() || input_count == 0 || input_count > 12 {
        return Err("Choose 1–12 ordered original sources".into());
    }
    let ids = (0..input_count)
        .map(|i| format!("source-{}", i + 1))
        .collect::<Vec<_>>();
    let mut fragments = Vec::with_capacity(input_count);
    for (i, row) in state.rows.iter().enumerate() {
        if row.fragment.len() > 20 || row.forward_name.len() > 256 || row.reverse_name.len() > 256 {
            return Err("Preparation names or fragment number exceed their field bounds".into());
        }
        if [
            &row.enzymes,
            &row.forward,
            &row.reverse_primer,
            &row.forward_tail,
            &row.reverse_tail,
        ]
        .iter()
        .any(|s| s.len() > 65536)
        {
            return Err("Preparation fields exceed the 64 KiB field limit".into());
        }
        let preparation = match row.preparation {
            OverlapPreparation::WholeLinear => Preparation::WholeLinear,
            OverlapPreparation::DeclaredDuplex => Preparation::DeclaredDuplex {
                left: row.left_terminus,
                right: row.right_terminus,
            },
            OverlapPreparation::CatalogDigest | OverlapPreparation::PcrCatalogDigest => {
                let selection = row.catalog_selection.as_ref().ok_or(
                    "Pin a catalog selection using Virtual digest or the PCR-product enzyme selector; then choose its fragment.",
                )?;
                selection.validate().map_err(|e| e.to_string())?;
                let fragment_index = row
                    .fragment
                    .trim()
                    .parse::<usize>()
                    .ok()
                    .and_then(|n| n.checked_sub(1))
                    .ok_or("Digest fragment numbers start at 1")?;
                if row.preparation == OverlapPreparation::PcrCatalogDigest {
                    Preparation::PcrCatalogDigest {
                        forward: Primer::new(&row.forward_name, &row.forward, &row.forward_tail)
                            .map_err(|e| e.to_string())?,
                        reverse: Primer::new(
                            &row.reverse_name,
                            &row.reverse_primer,
                            &row.reverse_tail,
                        )
                        .map_err(|e| e.to_string())?,
                        selection: selection.clone(),
                        fragment_index,
                        repair: row.repair,
                    }
                } else {
                    Preparation::CatalogDigest {
                        selection: selection.clone(),
                        fragment_index,
                        repair: row.repair,
                    }
                }
            }
            OverlapPreparation::Digest => {
                let enzymes = row
                    .enzymes
                    .split(',')
                    .map(str::trim)
                    .filter(|s| !s.is_empty())
                    .map(str::to_owned)
                    .collect::<Vec<_>>();
                if enzymes.is_empty() || enzymes.len() > 64 || enzymes.iter().any(|e| e.len() > 256)
                {
                    return Err("Choose 1–64 explicit restriction enzyme names".into());
                }
                let fragment_index = row
                    .fragment
                    .trim()
                    .parse::<usize>()
                    .ok()
                    .and_then(|n| n.checked_sub(1))
                    .ok_or("Digest fragment numbers start at 1")?;
                Preparation::Digest {
                    enzymes,
                    fragment_index,
                    repair: row.repair,
                }
            }
            OverlapPreparation::ExactPcr => Preparation::ExactPcr {
                forward: Primer::new(&row.forward_name, &row.forward, &row.forward_tail)
                    .map_err(|e| e.to_string())?,
                reverse: Primer::new(&row.reverse_name, &row.reverse_primer, &row.reverse_tail)
                    .map_err(|e| e.to_string())?,
            },
        };
        fragments.push(FragmentSpec {
            role: format!("fragment-{}", i + 1),
            input_id: ids[i].clone(),
            preparation,
            reverse: row.reverse,
        });
    }
    Ok((ids, fragments))
}

fn overlap_request(
    state: &MethodsState,
    inputs: Vec<Molecule>,
    target: &Molecule,
) -> Result<(Vec<PlanInput>, collagene_core::assembly_recipe::Recipe), String> {
    use collagene_core::assembly_recipe::Recipe;
    let (overlaps, options) = state.profile()?;
    let (ids, fragments) = prepared_rows(state, inputs.len())?;
    let recipe = Recipe {
        name: format!("Profiled overlap assembly to {}", target.name),
        fragments,
        overlaps,
        circular: state.circular,
        options,
    };
    let borrowed = ids
        .iter()
        .zip(&inputs)
        .map(|(id, m)| (id.as_str(), m))
        .collect::<Vec<_>>();
    collagene_core::assembly_recipe::source_preflight(&borrowed, target, &recipe)
        .map_err(|e| e.to_string())?;
    let bound = ids
        .into_iter()
        .zip(inputs)
        .map(|(id, molecule)| PlanInput {
            revision: molecule.revision(),
            id,
            molecule,
        })
        .collect();
    Ok((bound, recipe))
}

fn ligation_request(
    state: &MethodsState,
    inputs: Vec<Molecule>,
    target: &Molecule,
) -> Result<(Vec<PlanInput>, collagene_core::ligation_recipe::Recipe), String> {
    use collagene_core::ligation_recipe::{Fragment, Recipe};
    let (ids, fragments) = prepared_rows(state, inputs.len())?;
    if state.mode == Mode::GoldenGate {
        let geometry = state.golden_gate.enzyme.profile();
        for fragment in &fragments {
            let selection=match &fragment.preparation {
                collagene_core::assembly_recipe::Preparation::CatalogDigest{selection,..} | collagene_core::assembly_recipe::Preparation::PcrCatalogDigest{selection,..}=>selection,
                _=>return Err("Golden Gate rows require an existing pinned Type IIS digest fragment or an explicit PCR → catalog digest preparation. Use ordered ligation for other preparation chains.".into()),
            };
            if selection.enzymes.len() != 1
                || !selection.enzymes.iter().all(|entry| {
                    entry.enzyme.recognition == geometry.recognition
                        && entry.enzyme.top_cut == geometry.top_cut
                        && entry.enzyme.bottom_cut == geometry.bottom_cut
                })
            {
                return Err("A row's pinned enzyme geometry differs from the chosen Type IIS enzyme. Choose the matching existing digest fragment or apply fresh PCR-site proposals.".into());
            }
        }
    }
    let recipe = Recipe {
        name: format!(
            "{} to {}",
            if state.mode == Mode::GoldenGate {
                "Ordered Golden Gate preparation/ligation"
            } else {
                "Ordered ligation"
            },
            target.name
        ),
        fragments: fragments
            .into_iter()
            .zip(&state.rows)
            .map(|(source, row)| Fragment {
                source,
                treatment: match row.treatment {
                    Treatment::Keep => None,
                    Treatment::Phosphorylate => Some(EndTreatment::Phosphorylate),
                    Treatment::Dephosphorylate => Some(EndTreatment::Dephosphorylate),
                },
            })
            .collect(),
        circular: state.circular,
        allow_nicked: state.allow_nicked,
    };
    let borrowed = ids
        .iter()
        .zip(&inputs)
        .map(|(id, m)| (id.as_str(), m))
        .collect::<Vec<_>>();
    collagene_core::ligation_recipe::source_preflight(&borrowed, target, &recipe)
        .map_err(|e| e.to_string())?;
    let bound = ids
        .into_iter()
        .zip(inputs)
        .map(|(id, molecule)| PlanInput {
            revision: molecule.revision(),
            id,
            molecule,
        })
        .collect();
    Ok((bound, recipe))
}
fn method_status(status: &collagene_core::verification::Status) -> &'static str {
    match status {
        collagene_core::verification::Status::ConsistentUnderAssumptions => {
            "consistent under stated assumptions"
        }
        collagene_core::verification::Status::Incomplete => "incomplete",
        collagene_core::verification::Status::Invalid => "inconsistent",
    }
}
fn method_end(end: Option<&collagene_core::chemistry::End>) -> String {
    use collagene_core::chemistry::{FivePrime, Polarity, ThreePrime};
    let Some(end) = end else {
        return "End unavailable; no terminal chemistry inferred".into();
    };
    let shape = match end.polarity {
        Polarity::Blunt => "blunt".into(),
        Polarity::FivePrime => format!("5′ overhang {} (5′→3′)", end.overhang),
        Polarity::ThreePrime => format!("3′ overhang {} (5′→3′)", end.overhang),
    };
    format!(
        "{shape}; 5′ {}; 3′ {}",
        match end.chemistry.five_prime {
            FivePrime::Phosphate => "phosphate",
            FivePrime::Hydroxyl => "hydroxyl",
            FivePrime::Unknown => "chemistry unknown",
        },
        match end.chemistry.three_prime {
            ThreePrime::Hydroxyl => "hydroxyl",
            ThreePrime::Blocked => "blocked",
            ThreePrime::Unknown => "chemistry unknown",
        }
    )
}
fn method_duplex(summary: &collagene_core::chemistry::DuplexSummary) -> String {
    format!(
        "Top {} bases; bottom {} bases; bottom reference offset {}. {}. Retained nicks: {}. Covalently closed in model: {}.",
        summary.top_length,
        summary.bottom_length,
        summary.bottom_start,
        match summary.ends.as_ref() {
            Some((left, right)) => format!(
                "Left: {}. Right: {}",
                method_end(Some(left)),
                method_end(Some(right))
            ),
            None => "No free ends on the reconstructed circular duplex".into(),
        },
        summary.nicks.len(),
        summary.covalently_closed_in_model
    )
}

fn declared_terminus(end: Terminus) -> String {
    format!(
        "5′ {}; 3′ {}",
        match end.five_prime {
            FivePrime::Phosphate => "phosphate",
            FivePrime::Hydroxyl => "hydroxyl",
            FivePrime::Unknown => "unknown",
        },
        match end.three_prime {
            ThreePrime::Hydroxyl => "hydroxyl",
            ThreePrime::Blocked => "blocked",
            ThreePrime::Unknown => "unknown",
        }
    )
}
#[component]
fn DeclaredEndFields(config: RwSignal<MethodsState>, index: usize, left: bool) -> impl IntoView {
    let end = move || {
        config.with(|c| {
            c.rows
                .get(index)
                .map(|row| {
                    if left {
                        row.left_terminus
                    } else {
                        row.right_terminus
                    }
                })
                .unwrap_or(Terminus::UNKNOWN)
        })
    };
    let side = if left { "Left" } else { "Right" };
    view! {<fieldset><legend>{format!("{side} source end")}</legend><p class="helper">{if left{"Top-strand 5′ terminus and bottom-strand 3′ terminus, before orientation."}else{"Bottom-strand 5′ terminus and top-strand 3′ terminus, before orientation."}}</p>
        <label>"Declared 5′ group"<select aria-label=format!("Fragment {} {side} declared 5 prime",index+1) prop:value=move||match end().five_prime{FivePrime::Phosphate=>"phosphate",FivePrime::Hydroxyl=>"hydroxyl",FivePrime::Unknown=>"unknown"} on:change=move|e|{let value=match event_target_value(&e).as_str(){"phosphate"=>FivePrime::Phosphate,"hydroxyl"=>FivePrime::Hydroxyl,_=>FivePrime::Unknown};config.update(|c|if let Some(row)=c.rows.get_mut(index){if left{row.left_terminus.five_prime=value}else{row.right_terminus.five_prime=value}})}><option value="unknown">"Unknown"</option><option value="phosphate">"Phosphate (P)"</option><option value="hydroxyl">"Hydroxyl (OH)"</option></select></label>
        <label>"Declared 3′ group"<select aria-label=format!("Fragment {} {side} declared 3 prime",index+1) prop:value=move||match end().three_prime{ThreePrime::Hydroxyl=>"hydroxyl",ThreePrime::Blocked=>"blocked",ThreePrime::Unknown=>"unknown"} on:change=move|e|{let value=match event_target_value(&e).as_str(){"hydroxyl"=>ThreePrime::Hydroxyl,"blocked"=>ThreePrime::Blocked,_=>ThreePrime::Unknown};config.update(|c|if let Some(row)=c.rows.get_mut(index){if left{row.left_terminus.three_prime=value}else{row.right_terminus.three_prime=value}})}><option value="unknown">"Unknown"</option><option value="hydroxyl">"Hydroxyl (OH)"</option><option value="blocked">"Blocked"</option></select></label>
    </fieldset>}
}

fn primer_design_request(
    state: &MethodsState,
    inputs: Vec<Molecule>,
    target: &Molecule,
    form: &PlannerForm,
) -> Result<
    (
        Vec<PlanInput>,
        collagene_core::overlap_primer_design::Request,
    ),
    String,
> {
    use collagene_core::overlap_primer_design::{Fragment, Request};
    if state.mode != Mode::Overlap
        || inputs.len() != state.rows.len()
        || inputs.is_empty()
        || inputs.len() > 12
    {
        return Err("Choose 1–12 ordered source regions for overlap assembly.".into());
    }
    if [&form.overlap, &form.minimum, &form.maximum, &form.tm]
        .iter()
        .any(|s| s.len() > 32)
    {
        return Err("Primer settings exceed their field limits.".into());
    }
    let ids = (0..inputs.len())
        .map(|i| format!("source-{}", i + 1))
        .collect::<Vec<_>>();
    let mut fragments = Vec::with_capacity(inputs.len());
    for (i, (row, molecule)) in state.rows.iter().zip(&inputs).enumerate() {
        if row.interval_start.len() > 20 || row.interval_end.len() > 20 {
            return Err("Invalid primer design interval.".into());
        }
        let start = row
            .interval_start
            .trim()
            .parse::<usize>()
            .ok()
            .and_then(|n| n.checked_sub(1))
            .ok_or_else(|| {
                format!(
                    "Fragment {} start must be a 1-based source position.",
                    i + 1
                )
            })?;
        let end = if row.interval_end.trim().is_empty() {
            molecule.sequence.len()
        } else {
            row.interval_end
                .trim()
                .parse::<usize>()
                .map_err(|_| format!("Invalid end position for fragment {}.", i + 1))?
        };
        fragments.push(Fragment {
            role: format!("fragment-{}", i + 1),
            input_id: ids[i].clone(),
            start,
            end,
            reverse: row.reverse,
        });
    }
    let profile_options = ProfileOptions {
        method: state.method,
        scan_repeated_overlaps: state.scan_repeats,
        max_comparisons: state
            .comparisons
            .parse()
            .map_err(|_| "Invalid repeat comparison budget")?,
        max_occurrences: state
            .occurrences
            .parse()
            .map_err(|_| "Invalid repeat occurrence limit")?,
    };
    profile_options.validate().map_err(|e| e.to_string())?;
    let request = Request {
        name: format!("Suggested overlap PCR primers for {}", target.name),
        fragments,
        overlap: form
            .overlap
            .parse()
            .map_err(|_| "Invalid homology length")?,
        circular: state.circular,
        primer_options: collagene_core::primers::DesignOptions {
            min_length: form
                .minimum
                .parse()
                .map_err(|_| "Invalid minimum annealing length")?,
            max_length: form
                .maximum
                .parse()
                .map_err(|_| "Invalid maximum annealing length")?,
            target_tm_c: form.tm.parse().map_err(|_| "Invalid target annealing Tm")?,
            max_pairs: 1,
            ..Default::default()
        },
        profile_options,
    };
    let borrowed = ids
        .iter()
        .zip(&inputs)
        .map(|(id, m)| (id.as_str(), m))
        .collect::<Vec<_>>();
    collagene_core::overlap_primer_design::preflight(&borrowed, target, &request)
        .map_err(|e| e.to_string())?;
    let bound = ids
        .into_iter()
        .zip(inputs)
        .map(|(id, molecule)| PlanInput {
            revision: molecule.revision(),
            id,
            molecule,
        })
        .collect();
    Ok((bound, request))
}

fn golden_gate_request(
    state: &MethodsState,
    inputs: Vec<Molecule>,
    target: &Molecule,
) -> Result<
    (
        Vec<PlanInput>,
        collagene_core::golden_gate_primer_design::Request,
    ),
    String,
> {
    use collagene_core::golden_gate_primer_design::{Fragment, Request};
    if inputs.len() != state.rows.len() || inputs.is_empty() || inputs.len() > 12 {
        return Err("Choose 1–12 ordered original sources.".into());
    }
    let form = &state.golden_gate;
    if form.clamp.len() > 12
        || form.spacer.len() > 2
        || form.forward_length.len() > 4
        || form.reverse_length.len() > 4
    {
        return Err("Type IIS primer settings exceed their field bounds.".into());
    }
    let forward_length = form
        .forward_length
        .parse()
        .map_err(|_| "Enter a forward annealing length.")?;
    let reverse_length = form
        .reverse_length
        .parse()
        .map_err(|_| "Enter a reverse annealing length.")?;
    let ids = (0..inputs.len())
        .map(|i| format!("source-{}", i + 1))
        .collect::<Vec<_>>();
    let mut fragments = Vec::new();
    for (i, row) in state.rows.iter().enumerate() {
        if row.preparation != OverlapPreparation::PcrCatalogDigest {
            continue;
        }
        if row.interval_start.len() > 20
            || row.interval_end.len() > 20
            || row.right_fusion.len() > 16
        {
            return Err("Primer source interval or fusion field exceeds its bound.".into());
        }
        let start = row
            .interval_start
            .trim()
            .parse::<usize>()
            .ok()
            .and_then(|n| n.checked_sub(1))
            .ok_or("Source interval starts at a 1-based position.")?;
        let end = if row.interval_end.trim().is_empty() {
            inputs[i].sequence.len()
        } else {
            row.interval_end
                .trim()
                .parse()
                .map_err(|_| "Invalid source interval end.")?
        };
        fragments.push(Fragment {
            role: format!("fragment-{}", i + 1),
            input_id: ids[i].clone(),
            start,
            end,
            reverse: row.reverse,
            forward_length,
            reverse_length,
            right_overhang: row.right_fusion.clone(),
        });
    }
    if fragments.is_empty() {
        return Err("Choose Exact PCR → pinned catalog digest for at least one row before proposing primers. Existing-site rows remain unchanged.".into());
    }
    let request = Request {
        name: format!("Type IIS primer preparation for {}", target.name),
        fragments,
        enzyme: form.enzyme,
        clamp: form.clamp.clone(),
        spacer: form.spacer.clone(),
    };
    let borrowed = ids
        .iter()
        .zip(&inputs)
        .map(|(id, m)| (id.as_str(), m))
        .collect::<Vec<_>>();
    collagene_core::golden_gate_primer_design::source_preflight(&borrowed, target, &request)
        .map_err(|e| e.to_string())?;
    let inputs = ids
        .into_iter()
        .zip(inputs)
        .map(|(id, molecule)| PlanInput {
            revision: molecule.revision(),
            id,
            molecule,
        })
        .collect();
    Ok((inputs, request))
}

#[component]
fn MethodDiagnostics(
    diagnostics: collagene_core::assembly_diagnostics::Report,
    names: Vec<String>,
    preparations: Vec<String>,
    on_edit: Callback<usize>,
    #[prop(optional_no_strip)] on_source: Option<Callback<usize>>,
) -> impl IntoView {
    use collagene_core::assembly_diagnostics::{JunctionState, RowState, Stage};
    let status = method_status(&diagnostics.status);
    let stage_name = |stage: Stage| match stage {
        Stage::Pcr => "PCR",
        Stage::Digest => "restriction digest",
        Stage::Repair => "end repair",
        Stage::EndTreatment => "terminal treatment",
        Stage::SourceDeclaration => "source assumption",
        Stage::OverlapAssembly => "overlap assembly",
        Stage::LigationAssembly => "ligation",
        Stage::TargetComparison => "independent-target comparison",
    };
    view! {<section class="method-diagnostics" aria-label="Preparation and junction diagnostics">

        <h3>"Preparation and junction review"</h3><p>{status}</p>
        <p class="helper">"These locations come from the fresh typed step graph and retained product prefix. Later blocked rows have not been reconstructed; they are not additional failed experiments. Editing a field clears this evidence; viewing its editor does not change the draft."</p>
        {diagnostics.failure.map(|failure|view!{<article class="diagnostic-failure"><h4>{format!("First blocking result: step {} · {}",failure.step+1,stage_name(failure.stage))}</h4><p>{failure.message}</p><p class="helper">{failure.recovery}</p>
            <div class="method-row-actions">{failure.row.map(|row|view!{<button class="button secondary" on:click=move|_|on_edit.run(row)>"Edit preparation"</button>})}{failure.upstream.map(|row|view!{<button class="text-tool" on:click=move|_|on_edit.run(row)>"Edit upstream preparation"</button>})}{failure.downstream.map(|row|view!{<button class="text-tool" on:click=move|_|on_edit.run(row)>"Edit downstream preparation"</button>})}</div>
        </article>})}
        <div class="diagnostic-table"><table><thead><tr><th>"Original fragment / direction"</th><th>"Declared preparation"</th><th>"Reconstruction"</th><th>"Review"</th></tr></thead><tbody>
            {diagnostics.rows.into_iter().map(|row|{let index=row.index;let state=match row.state{RowState::InputReference=>"Original reference; no preparation step",RowState::Reconstructed=>"Reconstructed",RowState::Failed=>"Failed at its preparation",RowState::Blocked=>"Blocked by an earlier result"};view!{<tr><td>{format!("{} · {} · {}",index+1,names.get(index).map(String::as_str).unwrap_or(&row.input_id),if row.reverse{"reverse after preparation"}else{"forward"})}</td><td>{preparations.get(index).cloned().unwrap_or_else(||"Preparation unavailable".into())}</td><td>{state}<p>{format!("{}/{} declared steps reconstructed",row.completed_steps,row.total_steps)}</p>{row.preparation_step.map(|step|view!{<p>{format!("Preparation step {}",step+1)}</p>})}{row.repair_step.map(|step|view!{<p>{format!("Repair step {}",step+1)}</p>})}</td><td><div class="method-row-actions"><button class="text-tool" aria-label=format!("Edit fragment {} preparation",index+1) on:click=move|_|on_edit.run(index)>"Edit preparation"</button>{on_source.map(|source|view!{<button class="text-tool" on:click=move|_|source.run(index)>"View original source"</button>})}</div></td></tr>}}).collect_view()}
        </tbody></table></div>
        {(!diagnostics.junctions.is_empty()).then(||view!{<h4>"Fresh junction evidence"</h4>})}
        {diagnostics.junctions.into_iter().map(|junction|{let upstream=junction.upstream;let downstream=junction.downstream;let state=match junction.state{JunctionState::ExactHomology=>"Exact homology",JunctionState::HomologyMismatch=>"Incompatible homology",JunctionState::HomologyUnavailable=>"Homology unavailable",JunctionState::EndpointsAvailable=>"Ends available; compatibility is not inferred here",JunctionState::EndpointsUnavailable=>"Physical ends unavailable"};view!{<details><summary>{format!("Junction {} · fragment {} → {}{} · {}",junction.index+1,upstream+1,downstream+1,if junction.closing{" (circle closure)"}else{""},state)}</summary><p>{junction.message}</p><p class="helper">{junction.recovery}</p><div class="method-row-actions"><button class="text-tool" on:click=move|_|on_edit.run(upstream)>"Edit upstream preparation"</button><button class="text-tool" on:click=move|_|on_edit.run(downstream)>"Edit downstream preparation"</button>{on_source.map(|source|view!{<button class="text-tool" on:click=move|_|source.run(upstream)>"View upstream source"</button><button class="text-tool" on:click=move|_|source.run(downstream)>"View downstream source"</button>})}</div></details>}}).collect_view()}
        <ul>{diagnostics.warnings.into_iter().map(|warning|view!{<li class="helper">{warning}</li>}).collect_view()}</ul>
        <details><summary>"Full fresh findings"</summary><ul>{diagnostics.findings.into_iter().map(|finding|view!{<li>{finding}</li>}).collect_view()}</ul></details>
    </section>}
}
