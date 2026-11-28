//! Exact-source virtual restriction inspection. Scientific work stays in the worker.
pub use crate::source_identity::SourceInfo;
use crate::{
    app::prepare_download as download, run_origin_model::PinnedOrigin, worker::CachedWorker,
};
use collagene_core::{
    analysis,
    catalog_digest::{self, CatalogSelection, catalog},
    chemistry::{End, Polarity},
    digest_setups::{Recipe as SavedRecipe, Record as SavedDigest, Summary as SavedSummary},
    formats,
    protocol::{
        request::Request,
        schema::{AnalysisJob, AnalysisResult},
    },
    restriction_history,
    restriction_inspection::{
        self, InspectionLimit, InspectionOptions, InspectionSite, RestrictionInspection,
        SiteDiagnostic,
    },
    run_lineage::EventBinding,
};
use collagene_model::{Alphabet, Molecule, Topology};
use leptos::prelude::*;
use serde::Serialize;
use std::{
    rc::Rc,
    sync::{
        Arc,
        atomic::{AtomicBool, Ordering},
    },
};
use wasm_bindgen_futures::spawn_local;
#[derive(Clone, PartialEq, Eq)]
pub struct BoundaryPick {
    pub source: EventBinding,
    pub enzyme: restriction_inspection::EnzymeSpec,
    pub recognition_start: usize,
    pub reverse: bool,
}
pub use crate::digest_fragment_model::AssemblyFragment;
#[derive(Clone)]
pub struct Navigation {
    pub source: EventBinding,
    pub epoch: u64,
    pub spans: Vec<(usize, usize)>,
}
#[derive(Clone, PartialEq, Eq)]
struct Limits {
    comparisons: u64,
    sites: usize,
    fragments: usize,
    bytes: usize,
}
impl Limits {
    fn options(&self) -> InspectionOptions {
        InspectionOptions {
            max_comparisons: self.comparisons,
            max_sites: self.sites,
            max_fragments: self.fragments,
            max_report_bytes: self.bytes,
        }
    }
}
#[derive(Clone, PartialEq, Eq)]
struct Inputs {
    names: Vec<String>,
    limits: Limits,
    selection: CatalogSelection,
}
#[derive(Clone)]
struct Bound {
    origin: Option<PinnedOrigin>,
    origin_inspection: Option<Arc<collagene_core::material_derivation::DigestPreview>>,
    source: SourceInfo,
    molecule: Arc<Molecule>,
    inputs: Inputs,
    report: Arc<RestrictionInspection>,
    catalog_revision: String,
    catalog_engine_version: String,
    epoch: u64,
}
fn number<T: std::str::FromStr>(text: &str, name: &str) -> Result<T, String> {
    text.trim()
        .parse()
        .map_err(|_| format!("Enter a whole number for {name}."))
}
fn end_label(end: &End) -> String {
    match end.polarity {
        Polarity::Blunt => "Blunt".into(),
        Polarity::FivePrime => format!("5′ overhang {}", end.overhang),
        Polarity::ThreePrime => format!("3′ overhang {}", end.overhang),
    }
}
fn span_label(site: &InspectionSite) -> String {
    site.spans
        .iter()
        .map(|s| format!("{}–{}", s.start + 1, s.end))
        .collect::<Vec<_>>()
        .join(", ")
}
fn cut_label(value: Option<usize>, raw: i64) -> String {
    value
        .map(|n| format!("{n}"))
        .unwrap_or_else(|| format!("outside substrate ({raw})"))
}
fn limit_label(limit: InspectionLimit) -> &'static str {
    match limit {
        InspectionLimit::Comparisons => "Comparison limit",
        InspectionLimit::Sites => "Recognition-site limit",
        InspectionLimit::Fragments => "Fragment limit",
        InspectionLimit::ReportBytes => "Report size limit",
    }
}
fn site_diagnostic(d: SiteDiagnostic) -> &'static str {
    match d {
        SiteDiagnostic::PossibleRecognition => "Possible recognition: ambiguous bases",
        SiteDiagnostic::AmbiguousStrandGeometry => {
            "Palindromic recognition cannot resolve the declared asymmetric strand cuts"
        }
        SiteDiagnostic::CutsOutsideSubstrate => {
            "Recognition found; modeled cuts extend outside this linear substrate"
        }
    }
}
struct Json(Vec<u8>);
impl std::io::Write for Json {
    fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
        if bytes.len() > 32 * 1024 * 1024 - self.0.len() {
            return Err(std::io::Error::other(
                "Digest evidence exceeds the 32 MiB export bound",
            ));
        }
        self.0.extend_from_slice(bytes);
        Ok(bytes.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}
#[component]
pub fn DigestPanel(
    get_workspace: Callback<(), Result<collagene_core::workspace::Workspace, String>>,
    origin: RwSignal<Option<PinnedOrigin>>,
    source: Memo<SourceInfo>,
    get_molecule: Callback<(), Molecule>,
    epoch: RwSignal<u64>,
    on_navigate: Callback<Navigation>,
    on_open: Callback<crate::digest_product::DigestOpenRequest, Result<(), String>>,
    on_assembly: Callback<AssemblyFragment, Result<(), String>>,
    boundary_pick: RwSignal<Option<BoundaryPick>>,
    initial_enzymes: Vec<String>,
    on_methylation: Callback<()>,
    saved: Memo<Vec<SavedSummary>>,
    on_save: Callback<
        (
            String,
            catalog_digest::CatalogRestrictionInspection,
            SourceInfo,
        ),
        Result<(), String>,
    >,
    on_load: Callback<String, Result<SavedDigest, String>>,
    on_delete: Callback<String, Result<(), String>>,
    on_gel: Callback<(String, SavedRecipe), Result<(), String>>,
    on_fork: Callback<EventBinding, Result<SourceInfo, String>>,
) -> impl IntoView {
    let catalog = StoredValue::new(
        catalog_digest::builtin_catalog()
            .map(Arc::new)
            .map_err(|e| e.to_string()),
    );
    let enzyme_choices = StoredValue::new(catalog.with_value(|c| {
        c.as_ref()
            .map(|c| {
                let starter = analysis::starter_enzymes();
                let mut choices = c.edition().enzymes.clone();
                choices.sort_by_key(|enzyme| {
                    (
                        starter
                            .iter()
                            .position(|e| e.name == enzyme.name)
                            .unwrap_or(usize::MAX),
                        enzyme.name.clone(),
                    )
                });
                choices
            })
            .unwrap_or_default()
    }));
    let names = RwSignal::new(if initial_enzymes.is_empty() {
        vec!["EcoRI".to_string()]
    } else {
        initial_enzymes
    });
    let filter = RwSignal::new(String::new());
    let pinned_selection = RwSignal::new(None::<CatalogSelection>);
    let setup_name = RwSignal::new(String::new());
    let historical = RwSignal::new(None::<Arc<SavedDigest>>);
    let historical_hash = RwSignal::new(None::<String>);
    let saved_page = RwSignal::new(0usize);
    let comparisons = RwSignal::new("50000000".to_string());
    let max_sites = RwSignal::new("10000".to_string());
    let max_fragments = RwSignal::new("10000".to_string());
    let report_bytes = RwSignal::new("4194304".to_string());
    let inputs = Memo::new(move |_| -> Result<Inputs, String> {
        let catalog = catalog.get_value()?;
        let names = names.get();
        let ids = enzyme_choices.with_value(|choices| {
            choices
                .iter()
                .filter(|e| names.contains(&e.name))
                .map(|e| e.id.clone())
                .collect::<Vec<_>>()
        });
        let selection = if let Some(pin) = pinned_selection.get() {
            pin.validate_for(&catalog).map_err(|e| e.to_string())?;
            let pinned_names = pin
                .enzymes
                .iter()
                .map(|e| e.enzyme.name.clone())
                .collect::<Vec<_>>();
            if pinned_names != names {
                return Err("Loaded enzyme choices changed; explicitly reselect enzymes.".into());
            }
            pin
        } else {
            catalog_digest::select(&catalog, &ids).map_err(|e| e.to_string())?
        };
        Ok(Inputs {
            names,
            selection,
            limits: Limits {
                comparisons: number(&comparisons.get(), "comparison limit")?,
                sites: number(&max_sites.get(), "site limit")?,
                fragments: number(&max_fragments.get(), "fragment limit")?,
                bytes: number(&report_bytes.get(), "report byte limit")?,
            },
        })
    });
    let bound = RwSignal::new(None::<Arc<Bound>>);
    let status = RwSignal::new("Choose enzymes and run a virtual digest.".to_string());
    let busy = RwSignal::new(false);
    let product_busy = RwSignal::new(None::<usize>);
    let generation = RwSignal::new(0u64);
    let first_boundary = RwSignal::new(None::<collagene_core::digest_cut_selection::Boundary>);
    let second_boundary = RwSignal::new(None::<collagene_core::digest_cut_selection::Boundary>);
    let complement = RwSignal::new(false);

    let worker = StoredValue::new_local(None::<Rc<CachedWorker>>);
    let site_page = RwSignal::new(0usize);
    let fragment_page = RwSignal::new(0usize);
    let origin_finding_page = RwSignal::new(0usize);
    let alive = Arc::new(AtomicBool::new(true));
    let stop = Callback::new(move |()| {
        generation.update(|g| *g = g.wrapping_add(1));
        if let Some(w) = worker.get_value() {
            w.close();
        }
        worker.set_value(None);
        busy.set(false);
        product_busy.set(None);
    });
    {
        let alive = alive.clone();
        on_cleanup(move || {
            alive.store(false, Ordering::Relaxed);
            if let Some(w) = worker.get_value() {
                w.close();
            }
        });
    }
    Effect::new(move |_| {
        let _ = (source.get(), inputs.get(), epoch.get(), origin.get());
        stop.run(());
        bound.set(None);
        first_boundary.set(None);
        second_boundary.set(None);
        complement.set(false);
        site_page.set(0);
        fragment_page.set(0);
        origin_finding_page.set(0);
        status.set(
            "Choose enzymes and run a virtual digest. Results refresh only after an explicit run."
                .into(),
        );
    });
    let current = Callback::new(move |data: Arc<Bound>| {
        data.origin == origin.get_untracked()
            && data.epoch == epoch.get_untracked()
            && source.get_untracked() == data.source
            && inputs.get_untracked() == Ok(data.inputs.clone())
    });
    let choose_boundary = Callback::new(
        move |(first, boundary): (bool, collagene_core::digest_cut_selection::Boundary)| {
            use collagene_core::digest_cut_selection::Boundary;
            let Some(data) = bound.get_untracked() else {
                status.set(
                    "Select enzymes and run a complete digest before choosing cut boundaries."
                        .into(),
                );
                return;
            };
            if !current.run(data.clone())
                || !data.report.scan_complete
                || !data.report.products_determined
            {
                status.set(
                    "A fresh complete digest must determine the selected fragment first.".into(),
                );
                return;
            }
            let available = match &boundary {
                Boundary::SourceStart | Boundary::SourceEnd => {
                    data.source.topology == Topology::Linear
                }
                Boundary::Cut { site_index } => {
                    data.report.sites.get(*site_index).is_some_and(|s| {
                        s.index == *site_index
                            && s.certain
                            && s.top_cut.is_some()
                            && s.bottom_cut.is_some()
                    })
                }
            };
            if !available {
                status.set("This boundary is uncertain, unavailable or outside the linear source. Choose a certain cut or an actual linear DNA end.".into());
                return;
            }
            stop.run(());
            if first {
                first_boundary.set(Some(boundary));
            } else {
                second_boundary.set(Some(boundary));
            }
            status.set("Boundary selected. Review both strand cuts and the retained circular path before assigning the fragment.".into());
        },
    );
    Effect::new(move |_| {
        let Some(pick) = boundary_pick.get() else {
            return;
        };
        let active = source.get();
        let data = bound.get();
        let settings = inputs.get();
        let current_epoch = epoch.get();
        if pick.source != active.binding {
            status.set("The picked restriction site belongs to another document version. Pick a site on the current source or clear the selection.".into());
            return;
        }
        let Some(data) = data else {
            status.set(format!(
                "Select {} in Virtual digest and run a digest to choose this cut.",
                pick.enzyme.name
            ));
            return;
        };
        if data.epoch != current_epoch
            || data.source != active
            || settings != Ok(data.inputs.clone())
        {
            status
                .set("Digest settings changed; run the digest before using the picked cut.".into());
            return;
        }
        stop.run(());
        if !data.report.scan_complete || !data.report.products_determined {
            status.set("Complete the digest before selecting cut boundaries; recognition alone does not determine a physical fragment.".into());
            return;
        }
        let site = data.report.sites.iter().find(|site| {
            site.recognition_start == pick.recognition_start
                && site.reverse == pick.reverse
                && data.report.enzymes.get(site.enzyme_index) == Some(&pick.enzyme)
        });
        let Some(site) = site else {
            status.set(format!("The picked {} site is not in this digest. Select the enzyme and explicitly run the digest.",pick.enzyme.name));
            return;
        };
        stop.run(());
        if first_boundary.get_untracked().is_none() {
            choose_boundary.run((
                true,
                collagene_core::digest_cut_selection::Boundary::Cut {
                    site_index: site.index,
                },
            ));
        } else if second_boundary.get_untracked().is_none() {
            choose_boundary.run((
                false,
                collagene_core::digest_cut_selection::Boundary::Cut {
                    site_index: site.index,
                },
            ));
        } else {
            status.set("Both cut boundaries are set. Clear boundaries to start a new pair, or replace one using its site button.".into());
        }
        boundary_pick.set(None);
    });
    let run = {
        let alive = alive.clone();
        move |_| {
            let prepared = (|| -> Result<_, String> {
                let source = source.get_untracked();
                let inputs = inputs.get_untracked()?;
                let selection = inputs.selection.clone();
                let options = inputs.limits.options();
                let molecule = get_molecule.run(());
                let role =
                    restriction_history::source_role(&source.binding).map_err(|e| e.to_string())?;
                let bound_catalog = catalog.get_value()?;
                let origin_snapshot = origin.get_untracked();
                if molecule.revision() != source.binding.revision
                    || origin_snapshot
                        .as_ref()
                        .is_some_and(|pin| *pin.event() != source.binding)
                {
                    return Err("The source document/event changed. Run the digest again.".into());
                }
                let revision = if let Some(pin) = &origin_snapshot {
                    collagene_core::material_derivation::digest_inspection_revision(
                        &role,
                        &molecule,
                        &collagene_core::material_derivation::SourceOrigin::Saved(
                            pin.as_ref().clone(),
                        ),
                        &selection,
                        &options,
                    )
                } else {
                    catalog_digest::source_preflight(
                        &bound_catalog,
                        &selection,
                        &role,
                        &molecule,
                        &options,
                    )
                    .and_then(|()| {
                        catalog_digest::batch_revision(
                            &bound_catalog,
                            &selection,
                            &role,
                            &molecule,
                            &options,
                        )
                    })
                }
                .map_err(|e| e.to_string())?;
                stop.run(());
                let w = Rc::new(CachedWorker::new()?);
                worker.set_value(Some(w.clone()));
                busy.set(true);
                bound.set(None);
                site_page.set(0);
                fragment_page.set(0);
                Ok((
                    origin_snapshot,
                    source,
                    inputs,
                    selection,
                    options,
                    molecule,
                    role,
                    revision,
                    w,
                    generation.get_untracked(),
                    epoch.get_untracked(),
                ))
            })();
            let (
                origin_snapshot,
                source_snapshot,
                settings,
                selection,
                options,
                molecule,
                role,
                revision,
                w,
                expected_generation,
                expected_epoch,
            ) = match prepared {
                Ok(v) => v,
                Err(e) => {
                    status.set(e);
                    return;
                }
            };
            status.set("Scanning recognition sites and reconstructing virtual fragments in the background…".into());
            let alive = alive.clone();
            spawn_local(async move {
                let response = w
                    .ask(Request {
                        id: 1,
                        revision: revision.clone(),
                        job: if let Some(pin) = &origin_snapshot {
                            AnalysisJob::InspectMaterialDigest {
                                source_role: role.clone(),
                                source: molecule.clone(),
                                origin: collagene_core::material_derivation::SourceOrigin::Saved(
                                    pin.as_ref().clone(),
                                ),
                                selection: selection.clone(),
                                options: options.clone(),
                            }
                        } else {
                            AnalysisJob::InspectCatalogRestriction {
                                source_role: role.clone(),
                                molecule: molecule.clone(),
                                selection: selection.clone(),
                                options: options.clone(),
                            }
                        },
                    })
                    .await;
                if !alive.load(Ordering::Relaxed)
                    || generation.try_get_untracked() != Some(expected_generation)
                    || epoch.get_untracked() != expected_epoch
                    || source.get_untracked() != source_snapshot
                    || origin.get_untracked() != origin_snapshot
                    || inputs.get_untracked() != Ok(settings.clone())
                {
                    return;
                }
                w.close();
                worker.set_value(None);
                busy.set(false);
                let bound_catalog = match catalog.get_value() {
                    Ok(c) => c,
                    Err(e) => {
                        status.set(e);
                        return;
                    }
                };
                let (report, origin_inspection, catalog_revision, catalog_engine_version) =
                    match response {
                        Ok(response) if response.id == 1 && response.revision == revision => {
                            match response.result {
                                AnalysisResult::CatalogRestrictionInspection(report)
                                    if origin_snapshot.is_none() =>
                                {
                                    if let Err(e) = report.validate_bindings_for(
                                        &bound_catalog,
                                        &selection,
                                        &role,
                                        &molecule,
                                        &options,
                                    ) {
                                        status.set(format!("Digest report binding failed: {e}"));
                                        return;
                                    }
                                    (
                                        report.inspection,
                                        None,
                                        report.batch_revision,
                                        report.engine_version,
                                    )
                                }
                                AnalysisResult::MaterialDigestInspection(report)
                                    if origin_snapshot.is_some() =>
                                {
                                    if let Err(e) = report.validate_binding(&revision) {
                                        status.set(format!("Saved digest binding failed: {e}"));
                                        return;
                                    }
                                    if serde_json::to_vec(&report.selection).ok()
                                        != serde_json::to_vec(&selection).ok()
                                    {
                                        status.set("Saved digest catalog changed.".into());
                                        return;
                                    }
                                    (
                                        report.inspection.clone(),
                                        Some(Arc::new(report.clone())),
                                        report.context_revision,
                                        report.engine_version,
                                    )
                                }
                                AnalysisResult::Error(e) => {
                                    status.set(e);
                                    return;
                                }
                                _ => {
                                    status
                                        .set("Unexpected restriction inspection response.".into());
                                    return;
                                }
                            }
                        }
                        Ok(_) => {
                            status.set(
                                "Restriction response envelope does not match this request.".into(),
                            );
                            return;
                        }
                        Err(e) => {
                            status.set(e);
                            return;
                        }
                    };
                status.set(if report.products_determined{format!("{} recognition sites · {} virtual fragment(s), assuming complete cleavage.",report.sites.len(),report.fragments.len())}else if report.scan_complete{"Recognition scan finished; a unique complete product set is not determined. Review the evidence below.".into()}else{"Recognition scan is incomplete. Unsearched sequence cannot establish absence or a complete product set.".into()});
                bound.set(Some(Arc::new(Bound {
                    origin: origin_snapshot,
                    origin_inspection,
                    source: source_snapshot,
                    molecule: Arc::new(molecule),
                    inputs: settings,
                    report: Arc::new(report),
                    catalog_revision,
                    catalog_engine_version,
                    epoch: expected_epoch,
                })));
            });
        }
    };
    let navigate = Callback::new(move |index: usize| {
        let Some(data) = bound.get_untracked() else {
            return;
        };
        if !current.run(data.clone()) {
            status.set("Source or digest options changed. Recompute before navigation.".into());
            return;
        }
        if let Some(site) = data.report.sites.get(index) {
            on_navigate.run(Navigation {
                source: data.source.binding.clone(),
                epoch: data.epoch,
                spans: site.spans.iter().map(|s| (s.start, s.end)).collect(),
            });
        }
    });
    let open = Callback::new(move |index: usize| {
        let Some(data) = bound.get_untracked() else {
            return;
        };
        if !current.run(data.clone())
            || !data.report.products_determined
            || data.report.fragments.get(index).is_none()
        {
            status
                .set("A fresh current digest must determine this fragment before opening.".into());
            return;
        }
        let request = crate::digest_product::DigestOpenRequest {
            source: data.source.binding.clone(),
            origin: data.origin.clone(),
            selection: data.inputs.selection.clone(),
            options: data.inputs.limits.options(),
            fragment_index: index,
            epoch: data.epoch,
            guard: Callback::new(move |()| {
                if data.origin != origin.get()
                    || data.epoch != epoch.get()
                    || source.get() != data.source
                    || inputs.get() != Ok(data.inputs.clone())
                {
                    return Err("Digest source, origin or options changed; inspect again.".into());
                }
                Ok(())
            }),
        };
        if let Err(e) = on_open.run(request) {
            status.set(e);
        }
    });
    let assembly_alive = alive.clone();
    let use_in_assembly = Callback::new(move |index: usize| {
        if !assembly_alive.load(Ordering::Relaxed) {
            return;
        }
        let Some(data) = bound.get_untracked() else {
            return;
        };
        if busy.get_untracked()
            || !current.run(data.clone())
            || !data.report.scan_complete
            || !data.report.products_determined
        {
            status.set("Run a fresh complete digest before choosing an assembly fragment.".into());
            return;
        }
        let Some(fragment) = data.report.fragments.get(index) else {
            status.set("This digest fragment is unavailable.".into());
            return;
        };
        let n = data.molecule.sequence.len();
        let end = fragment.top_start.checked_add(fragment.duplex.top_length);
        if fragment.index != index
            || fragment.duplex.topology != Topology::Linear
            || fragment.duplex.top_length == 0
            || fragment.duplex.top_length > n
            || fragment.top_start >= n
            || end.is_none_or(|end| data.molecule.topology == Topology::Linear && end > n)
        {
            status.set(
                "Choose a nonempty linear digest product; an uncut circle has no assembly ends."
                    .into(),
            );
            return;
        }
        if fragment.duplex.bottom_length == 0
            || fragment.duplex.bottom_length
                > data
                    .origin_inspection
                    .as_ref()
                    .map(|saved| saved.source_duplex.bottom_length)
                    .unwrap_or(n)
            || i64::try_from(fragment.top_start)
                .ok()
                .and_then(|top| fragment.bottom_start.checked_sub(top))
                != Some(fragment.duplex.bottom_start)
        {
            status.set("Digest strand coordinates are unavailable or inconsistent.".into());
            return;
        }
        let payload = AssemblyFragment {
            material_domain: data.origin_inspection.as_ref().map(|preview| {
                crate::assembly_draft_model::MaterialDomain {
                    top_length: preview.source_duplex.top_length,
                    bottom_start: preview.source_duplex.bottom_start,
                    bottom_length: preview.source_duplex.bottom_length,
                    topology: preview.source_duplex.topology,
                    state_revision: preview.source_modifications.state_revision.clone(),
                }
            }),
            origin: data.origin.clone(),
            source: data.source.clone(),
            molecule: data.molecule.clone(),
            selection: data.inputs.selection.clone(),
            options: data.inputs.limits.options(),
            fragment_index: index,
            batch_revision: data.catalog_revision.clone(),
            top_start: fragment.top_start,
            top_length: fragment.duplex.top_length,
            bottom_start: fragment.bottom_start,
            bottom_length: fragment.duplex.bottom_length,
            epoch: data.epoch,
        };
        if !current.run(data.clone()) {
            status.set("Source or digest settings changed; run the digest again.".into());
            return;
        }
        match on_assembly.run(payload){Ok(())=>status.set("Digest fragment assigned to the active assembly slot with its original source and exact catalog selection. Review the intended construct and preparation.".into()),Err(error)=>status.set(error)}
    });
    let navigate_fragment = Callback::new(move |index: usize| {
        let Some(data) = bound.get_untracked() else {
            return;
        };
        if !current.run(data.clone()) {
            status.set("Source or digest settings changed; run the digest again.".into());
            return;
        }
        let Some(fragment) = data.report.fragments.get(index) else {
            return;
        };
        let n = data.molecule.sequence.len();
        let Some(end) = fragment.top_start.checked_add(fragment.duplex.top_length) else {
            return;
        };
        if n == 0
            || fragment.top_start >= n
            || fragment.duplex.top_length == 0
            || fragment.duplex.top_length > n
        {
            return;
        }
        let spans = if end <= n {
            vec![(fragment.top_start, end)]
        } else if data.molecule.topology == Topology::Circular {
            vec![(fragment.top_start, n), (0, end - n)]
        } else {
            return;
        };
        on_navigate.run(Navigation {
            source: data.source.binding.clone(),
            epoch: data.epoch,
            spans,
        });
    });
    let selection_alive = alive.clone();
    let use_boundaries = Callback::new(move |()| {
        if !selection_alive.load(Ordering::Relaxed) {
            return;
        }
        use collagene_core::digest_cut_selection::Choice;
        let Some(data) = bound.get_untracked() else {
            return;
        };
        if !current.run(data.clone())
            || !data.report.scan_complete
            || !data.report.products_determined
        {
            status.set("Run a complete current digest before using cut boundaries.".into());
            return;
        }
        let (Some(first), Some(second)) = (
            first_boundary.get_untracked(),
            second_boundary.get_untracked(),
        ) else {
            status.set("Choose both cut boundaries first.".into());
            return;
        };
        let choice = Choice {
            first,
            second,
            complement: complement.get_untracked(),
        };
        if let Some(saved) = &data.origin_inspection {
            let summary = match saved.selected_fragment(&choice) {
                Ok(fragment) => fragment,
                Err(error) => {
                    status.set(error.to_string());
                    return;
                }
            };
            let payload = AssemblyFragment {
                material_domain: data.origin_inspection.as_ref().map(|preview| {
                    crate::assembly_draft_model::MaterialDomain {
                        top_length: preview.source_duplex.top_length,
                        bottom_start: preview.source_duplex.bottom_start,
                        bottom_length: preview.source_duplex.bottom_length,
                        topology: preview.source_duplex.topology,
                        state_revision: preview.source_modifications.state_revision.clone(),
                    }
                }),
                origin: data.origin.clone(),
                source: data.source.clone(),
                molecule: data.molecule.clone(),
                selection: data.inputs.selection.clone(),
                options: data.inputs.limits.options(),
                fragment_index: summary.index,
                batch_revision: data.catalog_revision.clone(),
                top_start: summary.top_start,
                top_length: summary.duplex.top_length,
                bottom_start: summary.bottom_start,
                bottom_length: summary.duplex.bottom_length,
                epoch: data.epoch,
            };
            if !current.run(data.clone()) {
                return;
            }
            match on_assembly.run(payload){Ok(())=>status.set("Saved-material cut fragment assigned with its exact origin; assembly will reconstruct it again.".into()),Err(error)=>status.set(error)};
            return;
        }
        let selected_catalog = match catalog.get_value() {
            Ok(c) => c,
            Err(e) => {
                status.set(e);
                return;
            }
        };
        let role = match restriction_history::source_role(&data.source.binding) {
            Ok(v) => v,
            Err(e) => {
                status.set(e.to_string());
                return;
            }
        };
        let selection = data.inputs.selection.clone();
        let options = data.inputs.limits.options();
        let revision = match collagene_core::digest_cut_selection::selection_revision(
            &selected_catalog,
            &selection,
            &role,
            &data.molecule,
            &options,
            &data.catalog_revision,
            &choice,
        ) {
            Ok(v) => v,
            Err(e) => {
                status.set(e.to_string());
                return;
            }
        };
        stop.run(());
        let w = match CachedWorker::new() {
            Ok(v) => Rc::new(v),
            Err(e) => {
                status.set(e);
                return;
            }
        };
        worker.set_value(Some(w.clone()));
        busy.set(true);
        let token = generation.get_untracked();
        status.set(
            "Reconstructing the selected cut-defined fragment before assigning it to assembly…"
                .into(),
        );
        let alive = selection_alive.clone();
        spawn_local(async move {
            let result = w
                .ask(Request {
                    id: token as u32,
                    revision: revision.clone(),
                    job: AnalysisJob::SelectDigestFragment {
                        selection: selection.clone(),
                        source_role: role.clone(),
                        molecule: data.molecule.as_ref().clone(),
                        options: options.clone(),
                        expected_batch_revision: data.catalog_revision.clone(),
                        choice: choice.clone(),
                    },
                })
                .await;
            w.close();
            if !alive.load(Ordering::Relaxed)
                || generation.get_untracked() != token
                || !current.run(data.clone())
                || first_boundary.get_untracked() != Some(choice.first.clone())
                || second_boundary.get_untracked() != Some(choice.second.clone())
                || complement.get_untracked() != choice.complement
            {
                return;
            }
            worker.set_value(None);
            busy.set(false);
            let chosen = match result {
                Ok(response) => match response.result {
                    AnalysisResult::DigestFragmentSelected(value) => value,
                    AnalysisResult::Error(error) => {
                        status.set(error);
                        return;
                    }
                    _ => {
                        status.set("Unexpected cut-boundary selection response.".into());
                        return;
                    }
                },
                Err(error) => {
                    status.set(error);
                    return;
                }
            };
            if let Err(error) = chosen.validate_binding(&revision) {
                status.set(error.to_string());
                return;
            }
            if chosen.choice != choice || chosen.inspection.batch_revision != data.catalog_revision
            {
                status.set(
                    "Selected fragment belongs to a different digest or boundary pair.".into(),
                );
                return;
            }
            let index = chosen.product.product.fragment_index;
            if let Err(error) = chosen.inspection.validate_bindings_for(
                &selected_catalog,
                &selection,
                &role,
                &data.molecule,
                &options,
            ) {
                status.set(error.to_string());
                return;
            }
            if let Err(error) = chosen.product.validate_bindings_for(
                &selected_catalog,
                &selection,
                &role,
                &data.molecule,
                &options,
                index,
            ) {
                status.set(error.to_string());
                return;
            }
            let summary = &chosen.product.product.summary;
            let n = data.molecule.sequence.len();
            let end = summary.top_start.checked_add(summary.duplex.top_length);
            if summary.index != index
                || summary.duplex.topology != Topology::Linear
                || summary.top_start >= n
                || summary.duplex.top_length == 0
                || summary.duplex.top_length > n
                || summary.duplex.bottom_length == 0
                || summary.duplex.bottom_length > n
                || end.is_none_or(|end| data.source.topology == Topology::Linear && end > n)
                || i64::try_from(summary.top_start)
                    .ok()
                    .and_then(|top| summary.bottom_start.checked_sub(top))
                    != Some(summary.duplex.bottom_start)
            {
                status.set(
                    "The selected path does not describe one supported linear digest fragment."
                        .into(),
                );
                return;
            }
            let payload = AssemblyFragment {
                material_domain: data.origin_inspection.as_ref().map(|preview| {
                    crate::assembly_draft_model::MaterialDomain {
                        top_length: preview.source_duplex.top_length,
                        bottom_start: preview.source_duplex.bottom_start,
                        bottom_length: preview.source_duplex.bottom_length,
                        topology: preview.source_duplex.topology,
                        state_revision: preview.source_modifications.state_revision.clone(),
                    }
                }),
                origin: data.origin.clone(),
                source: data.source.clone(),
                molecule: data.molecule.clone(),
                selection,
                options,
                fragment_index: index,
                batch_revision: chosen.inspection.batch_revision,
                top_start: summary.top_start,
                top_length: summary.duplex.top_length,
                bottom_start: summary.bottom_start,
                bottom_length: summary.duplex.bottom_length,
                epoch: data.epoch,
            };
            if !current.run(data.clone()) || generation.get_untracked() != token {
                return;
            }
            match on_assembly.run(payload){Ok(())=>status.set("Cut-defined fragment assigned to the active assembly slot. Its original source, exact enzyme selection and both strand boundaries are retained.".into()),Err(error)=>status.set(error)}
        });
    });
    let clear_boundaries = move |_| {
        stop.run(());
        first_boundary.set(None);
        second_boundary.set(None);
        complement.set(false);
        boundary_pick.set(None);
        status.set(
            "Cut boundaries cleared. Choose the first and second boundary for a new pair.".into(),
        );
    };
    let export_product = {
        let alive = alive.clone();
        Callback::new(move |index: usize| {
            if origin.get_untracked().is_some() {
                status.set("Saved-material GenBank export is unavailable; use origin-preserving assembly output.".into());
                return;
            }
            let Some(data) = bound.get_untracked() else {
                return;
            };
            if busy.get_untracked()
                || !current.run(data.clone())
                || !data.report.products_determined
                || data.report.fragments.get(index).is_none()
            {
                status.set(
                    "Run a current digest that determines this fragment before exporting.".into(),
                );
                return;
            }
            let role = match restriction_history::source_role(&data.source.binding) {
                Ok(v) => v,
                Err(e) => {
                    status.set(e.to_string());
                    return;
                }
            };
            let selection = data.inputs.selection.clone();
            let bound_catalog = match catalog.get_value() {
                Ok(c) => c,
                Err(e) => {
                    status.set(e);
                    return;
                }
            };
            let options = data.inputs.limits.options();
            let revision = match catalog_digest::product_revision(
                &bound_catalog,
                &selection,
                &role,
                &data.molecule,
                &options,
                index,
            ) {
                Ok(v) => v,
                Err(e) => {
                    status.set(e.to_string());
                    return;
                }
            };
            stop.run(());
            let w = match CachedWorker::new() {
                Ok(w) => Rc::new(w),
                Err(e) => {
                    status.set(e);
                    return;
                }
            };
            worker.set_value(Some(w.clone()));
            busy.set(true);
            product_busy.set(Some(index));
            let expected_generation = generation.get_untracked();
            let alive = alive.clone();
            status.set("Reconstructing this fragment afresh before GenBank export…".into());
            spawn_local(async move {
                let response = w
                    .ask(Request {
                        id: 2,
                        revision: revision.clone(),
                        job: AnalysisJob::MaterializeCatalogDigestFragment {
                            source_role: role.clone(),
                            molecule: (*data.molecule).clone(),
                            selection: selection.clone(),
                            options: options.clone(),
                            fragment_index: index,
                        },
                    })
                    .await;
                if !alive.load(Ordering::Relaxed)
                    || generation.try_get_untracked() != Some(expected_generation)
                    || !current.run(data.clone())
                {
                    return;
                }
                w.close();
                worker.set_value(None);
                busy.set(false);
                product_busy.set(None);
                let product = match response {
                    Ok(response) if response.id == 2 && response.revision == revision => {
                        match response.result {
                            AnalysisResult::CatalogDigestProduct(product) => product,
                            AnalysisResult::Error(e) => {
                                status.set(e);
                                return;
                            }
                            _ => {
                                status.set("Unexpected digest-product response.".into());
                                return;
                            }
                        }
                    }
                    Ok(_) => {
                        status.set("Digest product response does not match this request.".into());
                        return;
                    }
                    Err(e) => {
                        status.set(e);
                        return;
                    }
                };
                if let Err(e) = product.validate_bindings_for(
                    &bound_catalog,
                    &selection,
                    &role,
                    &data.molecule,
                    &options,
                    index,
                ) {
                    status.set(format!("Fresh digest product binding failed: {e}"));
                    return;
                }
                let product = product.product;
                if product.inspection_revision != data.report.batch_revision {
                    status.set("Fresh product belongs to a different inspection.".into());
                    return;
                }
                match formats::export_with_report(&product.molecule, formats::TextFormat::GenBank) {
                    Ok(export) => {
                        let dispatch_data = data.clone();
                        let dispatch_alive = alive.clone();
                        crate::app::download_with_completion_guarded(
                            &format!("digest-fragment-{}.gb", index + 1),
                            "text/plain",
                            export.text.as_bytes(),
                            move || {
                                if !dispatch_alive.load(Ordering::Relaxed)
                                    || generation.try_get_untracked() != Some(expected_generation)
                                    || !current.run(dispatch_data.clone())
                                {
                                    return Err(
                                        "The digest source changed before export acceptance."
                                            .into(),
                                    );
                                }
                                Ok(())
                            },
                            move |result| {
                                if status.try_get_untracked().is_none()
                                    || !current.run(data.clone())
                                {
                                    return;
                                }
                                match result {
                                    Ok(()) => {
                                        let mut findings = product.projection_findings;
                                        findings.extend(export.warnings);
                                        status.set(format!("Exported freshly reconstructed fragment {}. GenBank producer/end comments are unverified provenance; keep the native workspace for exact source history. {}",index+1,findings.join(" ")));
                                    }
                                    Err(e) => status.set(e),
                                }
                            },
                        );
                    }
                    Err(e) => status.set(e.to_string()),
                }
            });
        })
    };
    let save_setup = Callback::new(move |()| {
        if origin.get_untracked().is_some() {
            status.set("Saved-material mode supports origin-preserving assembly handoff. Standalone fragment output, named setups and gel recipes are not yet supported.".into());
            return;
        }
        let Some(data) = bound.get_untracked() else {
            status.set("Run a fresh digest before saving its setup and evidence.".into());
            return;
        };
        if busy.get_untracked() || !current.run(data.clone()) {
            status.set("Digest source/options changed. Recompute before saving.".into());
            return;
        }
        let report = catalog_digest::CatalogRestrictionInspection {
            engine_version: data.catalog_engine_version.clone(),
            selection: data.inputs.selection.clone(),
            batch_revision: data.catalog_revision.clone(),
            inspection: data.report.as_ref().clone(),
        };
        status.set(match on_save.run((setup_name.get_untracked(), report, data.source.clone())) {
            Ok(()) => "Saved the named digest recipe and historical evidence. Reopening does not restore a current result.".into(), Err(error) => error,
        });
    });
    let add_current_to_gel = Callback::new(move |()| {
        if origin.get_untracked().is_some() {
            status.set("Saved-material mode supports origin-preserving assembly handoff. Standalone fragment output, named setups and gel recipes are not yet supported.".into());
            return;
        }
        let result = (|| -> Result<(), String> {
            let data = bound
                .get_untracked()
                .ok_or("Run a fresh digest before sending its exact recipe to gel.")?;
            if busy.get_untracked()
                || !current.run(data.clone())
                || epoch.get_untracked() != data.epoch
            {
                return Err(
                    "Digest source/options changed. Recompute before sending the recipe to gel."
                        .into(),
                );
            }
            let entered = setup_name.get_untracked();
            let name = if entered.is_empty() {
                "Current digest".to_string()
            } else {
                if entered.trim().is_empty()
                    || entered.trim() != entered
                    || entered.len() > 256
                    || entered.chars().any(char::is_control)
                {
                    return Err("Use a trimmed digest name of 1–256 bytes without control characters, or leave it blank for Current digest.".into());
                }
                entered
            };
            let recipe = SavedRecipe {
                source: data.source.binding.clone(),
                selection: data.inputs.selection.clone(),
                options: data.report.options.clone(),
                cleavage_mode: collagene_core::digest_setups::CleavageMode::AssumedCompleteMethylationUnassessed,
            };
            on_gel.run((name, recipe))
        })();
        status.set(match result { Ok(()) => "Added the exact current digest recipe to predicted gel. Prediction reconstructs it freshly and withholds unresolved products; cached fragment lengths are not reused.".into(), Err(error) => error });
    });
    let apply_saved_settings = Callback::new(
        move |(record, expected_source): (Arc<SavedDigest>, EventBinding)| -> Result<(), String> {
            if busy.get_untracked() {
                return Err("Wait for or cancel the current digest first.".into());
            }
            if source.get_untracked().binding != expected_source {
                return Err(
                    "The active source event changed; saved settings were not applied.".into(),
                );
            }
            let catalog = catalog.get_value()?;
            record
                .recipe
                .selection
                .validate_for(&catalog)
                .map_err(|e| e.to_string())?;
            let options = &record.recipe.options;
            batch(|| {
                origin.set(None);
                names.set(
                    record
                        .recipe
                        .selection
                        .enzymes
                        .iter()
                        .map(|e| e.enzyme.name.clone())
                        .collect(),
                );
                pinned_selection.set(Some(record.recipe.selection.clone()));
                comparisons.set(options.max_comparisons.to_string());
                max_sites.set(options.max_sites.to_string());
                max_fragments.set(options.max_fragments.to_string());
                report_bytes.set(options.max_report_bytes.to_string());
                setup_name.set(record.name.clone());
            });
            Ok(())
        },
    );
    let use_saved_settings = Callback::new(move |record: Arc<SavedDigest>| {
        let result = if source.get_untracked().binding != record.recipe.source {
            Err("These settings belong to a different source event. Use the explicit fork action to continue from that saved version.".into())
        } else {
            apply_saved_settings.run((record.clone(), record.recipe.source.clone()))
        };
        status.set(match result { Ok(()) => "Loaded the exact saved catalog pin and limits. Run virtual digest explicitly to create fresh evidence.".into(), Err(error) => error });
    });
    let fork_saved_source = Callback::new(move |record: Arc<SavedDigest>| -> Result<(), String> {
        if busy.get_untracked() {
            return Err("Cancel the current digest before forking a saved source.".into());
        }
        let fork = on_fork.run(record.recipe.source.clone())?;
        let result = (|| -> Result<(), String> {
            if fork.binding == record.recipe.source {
                return Err("The fork callback did not create a new source event.".into());
            }
            if fork.binding.revision != record.recipe.source.revision {
                return Err(
                    "The new fork does not contain the exact saved source revision.".into(),
                );
            }
            if source.get_untracked().binding != fork.binding {
                return Err("The new fork is not the exact active source event.".into());
            }
            apply_saved_settings.run((record, fork.binding))
        })();
        result.map_err(|error| {
            format!("The source fork exists, but saved settings were not applied: {error}")
        })
    });
    let export_report = move |_| {
        let Some(data) = bound.get_untracked() else {
            return;
        };
        if !current.run(data.clone()) {
            status
                .set("Digest source/options changed. Recompute before exporting evidence.".into());
            return;
        }
        #[derive(Serialize)]
        struct Evidence<'a> {
            kind: &'static str,
            authority: &'static str,
            source: &'a SourceInfo,
            saved_origin: Option<&'a collagene_core::run_material::MaterialOrigin>,
            saved_inspection: Option<&'a collagene_core::material_derivation::DigestPreview>,
            report: &'a RestrictionInspection,
            catalog_pin: &'a CatalogSelection,
            catalog_engine_version: &'a str,
            catalog_batch_revision: &'a str,
        }
        let artifact = Evidence {
            kind: "virtual-digest-review-v2",
            authority: "Sequence/geometry prediction under stated cleavage assumptions; saved evidence never restores current verification",
            source: &data.source,
            saved_origin: data.origin.as_ref().map(PinnedOrigin::as_ref),
            saved_inspection: data.origin_inspection.as_deref(),
            report: &data.report,
            catalog_pin: &data.inputs.selection,
            catalog_engine_version: &data.catalog_engine_version,
            catalog_batch_revision: &data.catalog_revision,
        };
        let mut bytes = Json(vec![]);
        match serde_json::to_writer_pretty(&mut bytes, &artifact) {
            Ok(()) => crate::app::download_with_completion(
                "virtual-digest.json",
                "application/json",
                &bytes.0,
                move |result| {
                    if status.try_get_untracked().is_none() {
                        return;
                    }
                    if let Err(e) = result {
                        status.set(e)
                    }
                },
            ),
            Err(e) => status.set(e.to_string()),
        }
    };
    view! {
            <section class="inspector-content digest-panel"><h2>"Virtual digest"</h2><crate::run_origin_view::SourceOriginPicker source=Memo::new(move |_|Some(source.get())) epoch get_workspace selected=Memo::new(move |_|origin.get()) on_select=Callback::new(move |pin:Option<PinnedOrigin>|{if pin.as_ref().is_some_and(|p|*p.event()!=source.get_untracked().binding){return Err("Saved origin source changed.".into());}origin.set(pin);Ok(())})/><Show when=move||origin.with(Option::is_some)><p class="helper">"Saved-material mode: assembly retains the exact originating run. Fragment opening and assembly preserve the saved origin. GenBank export, named setups and predicted-gel recipes are unavailable for saved material."</p></Show><p class="helper">"Choose enzymes to inspect recognition sites and predicted fragments. Products assume complete cleavage; methylation, reaction conditions and laboratory success are not assessed."</p>
            <details class="saved-digest-setups"><summary>"Named digest setups"</summary>
                <p class="helper">"Save the exact source event, ordered catalog geometry, cleavage assumption, limits and historical evidence together. Reopening a saved setup does not move document heads or restore a current result."</p>
                <label>"Digest name"<input aria-label="Saved digest setup name" maxlength="256" placeholder="Name this digest" prop:value=move||setup_name.get() on:input=move|e|{setup_name.set(event_target_value(&e));}/></label>
                <button class="button secondary full" disabled=move||{busy.get()||origin.with(Option::is_some)||bound.with(Option::is_none)} on:click=move|_|{save_setup.run(());}>"Save fresh digest setup and evidence"</button>
                <button class="button secondary full" disabled=move||{busy.get()||origin.with(Option::is_some)||bound.with(Option::is_none)} on:click=move|_|{add_current_to_gel.run(());}>"Add current digest to predicted gel"</button>
                <p class="helper">"A saved setup is optional for gel prediction. The exact current recipe is reconstructed again; incomplete or undetermined products cannot supply trusted band lengths."</p>
                {move||{
                    let entries=saved.get();let count=entries.len();let pages=count.div_ceil(100).max(1);let page=saved_page.get().min(pages-1);
                    let rows=entries.into_iter().skip(page*100).take(100).map(|entry|{
                        let load_hash=entry.attachment_hash.clone();let delete_hash=entry.attachment_hash.clone();
                        view!{<article class="saved-digest-entry"><strong>{entry.name}</strong><p class="helper">{format!("{} · {} · saved scan {} · products {}",entry.source_name,entry.catalog_edition,if entry.scan_complete{"complete"}else{"incomplete"},if entry.products_determined{"determined"}else{"undetermined"})}</p><div class="digest-actions">
                            <button class="text-tool" on:click=move|_|{match on_load.run(load_hash.clone()){Ok(record)=>{historical.set(Some(Arc::new(record)));historical_hash.set(Some(load_hash.clone()));},Err(error)=>status.set(error)}}>"Open saved digest"</button>
                            <button class="text-tool" on:click=move|_|{match on_delete.run(delete_hash.clone()){Ok(())=>{if historical_hash.get_untracked().as_ref()==Some(&delete_hash){historical.set(None);historical_hash.set(None);}status.set("Deleted saved digest setup.".into());},Err(error)=>status.set(error)}}>"Delete"</button>
                        </div></article>}
                    }).collect_view();
                    view!{<div>{rows}<p class="helper">{format!("{} saved digests · page {} of {}",count,page+1,pages)}</p><div class="digest-actions"><button class="text-tool" disabled={page==0} on:click=move|_|{saved_page.set(page.saturating_sub(1));}>"Previous saved digests"</button><button class="text-tool" disabled={page+1>=pages} on:click=move|_|{saved_page.set((page+1).min(pages-1));}>"Next saved digests"</button></div></div>}
                }}
                {move||historical.get().map(|record|view!{<SavedDigestView record on_settings=use_saved_settings on_gel on_fork=fork_saved_source/>})}
            </details>
            <label>"Find an enzyme"<input aria-label="Find digest enzyme" prop:value=move||filter.get() on:input=move|e|filter.set(event_target_value(&e)) placeholder="Name or recognition sequence"/></label>
            <div class="digest-enzyme-picker">{move||{let needle=filter.get().trim().to_ascii_lowercase();enzyme_choices.get_value().into_iter().filter(|e|e.name.to_ascii_lowercase().contains(&needle)||e.geometry.recognition.to_ascii_lowercase().contains(&needle)).map(|enzyme|{let name=enzyme.name.clone();let change=enzyme.name.clone();view!{<label class="checkbox-label digest-enzyme-choice"><input type="checkbox" data-enzyme=enzyme.name.clone() data-catalog-id=enzyme.id.clone() aria-label=format!("Digest enzyme {}",enzyme.name) prop:checked=move||names.with(|names|names.contains(&name)) on:change=move|e|{let checked=event_target_checked(&e);pinned_selection.set(None);names.update(|names|{names.retain(|n|n!=&change);if checked{names.push(change.clone());}});}/><span><strong>{enzyme.name}</strong><small>{enzyme.geometry.recognition}</small></span></label>}}).collect_view()}}</div>
            <p class="helper">{move||format!("Selected: {}",names.get().join(", "))}</p>
            <details class="digest-catalog" data-catalog-edition=move||catalog.with_value(|c|c.as_ref().map(|c|c.edition().edition.clone()).unwrap_or_default()) data-catalog-hash=move||catalog.with_value(|c|c.as_ref().map(|c|c.content_hash().to_string()).unwrap_or_default())><summary>"Catalog and formulation evidence"</summary><p class="helper">{move||catalog.with_value(|c|match c {Ok(c)=>format!("{} enzyme identities · edition {} · content hash {}",c.edition().enzymes.len(),c.edition().edition,c.content_hash()),Err(e)=>format!("Catalog unavailable: {e}")})}</p><p class="helper">"Recognition and strand-cut geometry determine the virtual digest. Same-cut enzymes remain distinct identities; vendor formulations and their documentary conditions are reviewed separately and never change or validate predicted cleavage."</p>{move||catalog.with_value(|c|c.as_ref().ok().map(|c|catalog_evidence(c,&names.get())))}</details><div class="digest-actions"><button class="button primary" disabled=move||busy.get()||names.with(Vec::is_empty)||source.get().alphabet!=Alphabet::Dna on:click=run>"Run virtual digest"</button><button class="button secondary" disabled=move||!busy.get() on:click=move|_|{stop.run(());status.set("Digest request cancelled; no pending result or product was accepted.".into());}>"Cancel digest"</button></div>
            <p class="digest-status helper" role="status">{move||status.get()}</p>
            <details class="digest-options"><summary>"Assumptions and analysis limits"</summary><p class="helper">"The pinned catalog selection describes recognition and strand-cut geometry. Documentary vendor formulation evidence is separate and is not a reaction-compatibility assessment. Recognition can be certain while cuts lack sufficient substrate, and ambiguous sites can prevent a complete product set."</p><button class="text-tool" on:click=move|_|on_methylation.run(())>"Review methylation assumptions"</button><label>"Maximum comparisons"<input aria-label="Digest comparison limit" prop:value=move||comparisons.get() on:input=move|e|comparisons.set(event_target_value(&e))/></label><label>"Maximum recognition sites"<input aria-label="Digest site limit" prop:value=move||max_sites.get() on:input=move|e|max_sites.set(event_target_value(&e))/></label><label>"Maximum fragments"<input aria-label="Digest fragment limit" prop:value=move||max_fragments.get() on:input=move|e|max_fragments.set(event_target_value(&e))/></label><label>"Maximum report bytes"<input aria-label="Digest report byte limit" prop:value=move||report_bytes.get() on:input=move|e|report_bytes.set(event_target_value(&e))/></label></details>
            {move||bound.get().map(|data|{let sites=data.report.sites.len();let fragments=data.report.fragments.len();let first_report=data.report.clone();let second_report=data.report.clone();let linear_source=data.source.topology==Topology::Linear;view!{<div class="digest-results" data-scan-complete=data.report.scan_complete.to_string() data-products-determined=data.report.products_determined.to_string()>{data.origin_inspection.as_ref().map(|saved|{let saved=saved.clone();let total=saved.findings.len();view!{<section class="saved-digest-material"><h3>"Fresh saved-material reconstruction"</h3>{crate::end_view::summary(&saved.source_duplex)}<p>{format!("Inherited incomplete assessment: {}. Original source modifications: {} known modified · {} eligible unknown · {} known unmodified.",saved.incomplete,saved.source_modifications.known_modified,saved.source_modifications.unknown_eligible,saved.source_modifications.known_unmodified)}</p><p class="helper">"Cleavage geometry assumes complete digestion; modification susceptibility and reaction performance remain unassessed. Origin limitations remain applicable."</p>{move||saved.findings.iter().skip(origin_finding_page.get()*20).take(20).map(|finding|view!{<p class="digest-warning">{finding.clone()}</p>}).collect_view()}<div class="digest-actions"><button class="text-tool" disabled=move||origin_finding_page.get()==0 on:click=move|_|origin_finding_page.update(|page|*page=page.saturating_sub(1))>"Previous origin findings"</button><button class="text-tool" disabled={move||(origin_finding_page.get()+1)*20>=total} on:click=move|_|origin_finding_page.update(|page|*page+=1)>"Next origin findings"</button></div></section>}})}<h3>"Recognition and products"</h3><p>{format!("Recognition scan: {}. Product set: {}.",if data.report.scan_complete{"finished"}else{"incomplete"},if data.report.products_determined{"determined under complete-cleavage assumptions"}else{"not determined"})}</p>{data.report.limits.iter().map(|limit|view!{<p class="digest-warning">{format!("{} reached. Incomplete evidence cannot establish absence or complete products.",limit_label(*limit))}</p>}).collect_view()}{data.report.diagnostics.iter().map(|d|view!{<p class="digest-warning">{d.clone()}</p>}).collect_view()}
            <section class="digest-cut-selection"><h3>"Choose a cut-defined assembly fragment"</h3><p class="helper">"Pick the first cut, then the second in the source sequence/map, or use the buttons below. Circular paths follow first → second; keep the opposite path for a retained backbone. A same-cut pair can describe a full circle opened at that cut. Internal selected-enzyme cuts must not split the retained path."</p><p>{move||cut_boundary_label(first_boundary.get().as_ref(),&first_report,"First")}</p><p>{move||cut_boundary_label(second_boundary.get().as_ref(),&second_report,"Second")}</p><Show when=move||linear_source><div class="digest-actions"><button class="text-tool" disabled=move||busy.get() on:click=move|_|choose_boundary.run((true,collagene_core::digest_cut_selection::Boundary::SourceStart))>"First: source start"</button><button class="text-tool" disabled=move||busy.get() on:click=move|_|choose_boundary.run((true,collagene_core::digest_cut_selection::Boundary::SourceEnd))>"First: source end"</button><button class="text-tool" disabled=move||busy.get() on:click=move|_|choose_boundary.run((false,collagene_core::digest_cut_selection::Boundary::SourceStart))>"Second: source start"</button><button class="text-tool" disabled=move||busy.get() on:click=move|_|choose_boundary.run((false,collagene_core::digest_cut_selection::Boundary::SourceEnd))>"Second: source end"</button></div></Show><Show when=move||!linear_source><label class="checkbox-label"><input type="checkbox" aria-label="Keep complementary digest cut path" prop:checked=move||complement.get() on:change=move|e|{stop.run(());complement.set(event_target_checked(&e));}/>"Keep opposite path (outside selected cut interval)"</label></Show><div class="digest-actions"><button class="text-tool" on:click=clear_boundaries>"Clear boundaries / start new pair"</button><button class="button primary" disabled=move||busy.get()||first_boundary.with(Option::is_none)||second_boundary.with(Option::is_none) on:click=move|_|use_boundaries.run(())>"Use selected cut boundaries in assembly"</button></div></section>
    <h3>{format!("{sites} recognition sites")}</h3><p class="helper">"Recognition ranges are original source bases, 1-based inclusive. Strand cuts are boundaries measured from the source origin; 0 is before the first base."</p>{data.report.sites.iter().skip(site_page.get()*25).take(25).map(|site|{let index=site.index;view!{<article class="digest-site"><button class="result-entry" on:click=move|_|navigate.run(index)>{format!("{} · {} · {}{}",site.enzyme,span_label(site),if site.reverse{"reverse"}else{"forward"},if site.certain{""}else{" · possible"})}</button><p>{format!("Top / bottom cut boundaries: {} / {}",cut_label(site.top_cut,site.top_boundary),cut_label(site.bottom_cut,site.bottom_boundary))}</p><div class="digest-actions"><button class="text-tool" disabled=move||busy.get() on:click=move|_|{choose_boundary.run((true,collagene_core::digest_cut_selection::Boundary::Cut{site_index:index}));navigate.run(index);}>"Use as first boundary"</button><button class="text-tool" disabled=move||busy.get() on:click=move|_|{choose_boundary.run((false,collagene_core::digest_cut_selection::Boundary::Cut{site_index:index}));navigate.run(index);}>"Use as second boundary"</button></div>{site.diagnostics.iter().map(|d|view!{<p class="digest-warning">{site_diagnostic(*d)}</p>}).collect_view()}</article>}}).collect_view()}<div class="digest-actions"><button class="button secondary" disabled=move||site_page.get()==0 on:click=move|_|site_page.update(|p|*p=p.saturating_sub(1))>"Previous digest sites"</button><button class="button secondary" disabled={move||(site_page.get()+1)*25>=sites} on:click=move|_|site_page.update(|p|*p+=1)>"Next digest sites"</button></div>
            <h3>"Reference-strand fragment lengths"</h3><p class="fragment-sizes">{data.report.ref_top_lengths.as_ref().map(|lengths|lengths.iter().take(100).map(|n|format!("{n} nt")).collect::<Vec<_>>().join(" · ")).unwrap_or_else(||"Not determined from this recognition evidence.".into())}</p><p class="helper">"First 100 lengths shown, following the reference/top strand and sorted by size. The JSON evidence retains every determined length. Physical summaries below use ascending reference cut order; staggered ends can give different top and bottom lengths."</p>
            <h3>{format!("{fragments} virtual fragment(s)")}</h3>{data.report.fragments.iter().skip(fragment_page.get()*25).take(25).map(|fragment|{let index=fragment.index;let circular=fragment.duplex.topology==Topology::Circular;let ends=fragment.duplex.ends.as_ref().map(|(a,b)|format!("Left: {} · right: {}",end_label(a),end_label(b))).unwrap_or_else(||"Uncut circular substrate: no free ends".into());view!{<article class="digest-fragment"><h4>{format!("Fragment {} · {}",index+1,if fragment.duplex.topology==Topology::Circular{"circular"}else{"linear"})}</h4><p>{format!("{} nt top / {} nt bottom · original reference start {}",fragment.duplex.top_length,fragment.duplex.bottom_length,fragment.top_start+1)}</p><p>{ends}</p>{data.origin_inspection.as_ref().and_then(|saved|saved.fragment_modifications.get(index)).map(|mods|view!{<p class="helper">{format!("Retained positional calls: {} modified · {} eligible unknown · {} unmodified",mods.known_modified,mods.unknown_eligible,mods.known_unmodified)}</p>})}<div class="digest-actions"><button class="button secondary" disabled=move||busy.get() on:click=move|_|navigate_fragment.run(index)>"View fragment in source"</button><button class="button primary" disabled=move||busy.get()||circular on:click=move|_|use_in_assembly.run(index)>"Use in active assembly slot"</button><button class="button secondary" disabled=move||busy.get() on:click=move|_|open.run(index)>"Open fragment"</button><button class="button secondary" disabled=move||busy.get()||origin.with(Option::is_some) on:click=move|_|export_product.run(index)>{move||if product_busy.get()==Some(index){"Reconstructing…"}else{"Export fragment GenBank"}}</button></div></article>}}).collect_view()}<div class="digest-actions"><button class="button secondary" disabled=move||fragment_page.get()==0 on:click=move|_|fragment_page.update(|p|*p=p.saturating_sub(1))>"Previous digest fragments"</button><button class="button secondary" disabled={move||(fragment_page.get()+1)*25>=fragments} on:click=move|_|fragment_page.update(|p|*p+=1)>"Next digest fragments"</button></div><button class="button secondary full" on:click=export_report>"Download digest evidence"</button><p class="helper">"Opening reconstructs the fragment again and retains its exact source event in the workspace. GenBank comments preserve descriptive, unverified producer/end information; native workspace files retain full source history."</p></div>}})}
            </section>
        }.into_any()
}

#[component]
fn SavedDigestView(
    record: Arc<SavedDigest>,
    on_settings: Callback<Arc<SavedDigest>>,
    on_gel: Callback<(String, SavedRecipe), Result<(), String>>,
    on_fork: Callback<Arc<SavedDigest>, Result<(), String>>,
) -> impl IntoView {
    let bound = StoredValue::new(record);
    let site_page = RwSignal::new(0usize);
    let fragment_page = RwSignal::new(0usize);
    let status = RwSignal::new(String::new());
    let catalog_status = bound.with_value(|record| {
        catalog_digest::builtin_catalog()
            .map_err(|e| e.to_string())
            .and_then(|catalog| {
                record
                    .recipe
                    .selection
                    .validate_for(&catalog)
                    .map_err(|e| e.to_string())
            })
            .err()
    });
    view! {
        <section class="saved-digest-review">
            <h3>{bound.with_value(|r|format!("{} · saved historical digest",r.name))}</h3>
            <p class="helper">"Historical inspection evidence only. This view does not establish current freshness, measured cleavage or a verified cloning product. Source heads and current digest results are unchanged."</p>
            <details><summary>"Exact source, catalog and options"</summary>
                <p>{bound.with_value(|r|format!("Document {} · event {} · revision {}",r.recipe.source.document_id,r.recipe.source.event_id,r.recipe.source.revision))}</p>
                <p>{bound.with_value(|r|format!("Catalog {} · content {} · saved engine {} · batch {}",r.recipe.selection.catalog_edition,r.recipe.selection.catalog_content_hash,r.report.engine_version,r.report.batch_revision))}</p>
                <p>{bound.with_value(|r|format!("Assumed complete cleavage; methylation and reaction conditions unassessed. Limits: {} comparisons · {} sites · {} fragments · {} report bytes.",r.recipe.options.max_comparisons,r.recipe.options.max_sites,r.recipe.options.max_fragments,r.recipe.options.max_report_bytes))}</p>
                <ul>{bound.with_value(|r|r.recipe.selection.enzymes.iter().map(|enzyme|view!{<li>{format!("{} [{}] · {} · top cut {} · bottom cut {}",enzyme.enzyme.name,enzyme.id,enzyme.enzyme.recognition,enzyme.enzyme.top_cut,enzyme.enzyme.bottom_cut)}</li>}).collect_view())}</ul>
            </details>
            {catalog_status.map(|error|view!{<p class="digest-warning">{format!("The saved catalog pin cannot currently be replayed: {error}. Historical evidence remains inspectable; no enzyme substitution is made.")}</p>})}
            <p class="helper">{bound.with_value(|r|format!("Saved recognition scan: {} · saved product set: {} · {} sites · {} fragments.",if r.report.inspection.scan_complete{"complete"}else{"incomplete"},if r.report.inspection.products_determined{"determined under its stated model"}else{"undetermined"},r.report.inspection.sites.len(),r.report.inspection.fragments.len()))}</p>
            {bound.with_value(|r|r.report.inspection.limits.iter().map(|limit|view!{<p class="digest-warning">{format!("Saved limit reached: {}",limit_label(*limit))}</p>}).collect_view())}
            {bound.with_value(|r|r.report.inspection.diagnostics.iter().map(|diagnostic|view!{<p class="digest-warning">{diagnostic.clone()}</p>}).collect_view())}
            <div class="digest-actions">
                <button class="button secondary" on:click=move|_|{on_settings.run(bound.get_value());}>"Use settings on exact current source"</button>
                <button class="button secondary" on:click=move|_|{let record=bound.get_value();status.set(match on_fork.run(record){Ok(())=>"Opened the exact source fork and loaded the saved catalog pin and limits. Run virtual digest explicitly; the saved historical result is unchanged.".into(),Err(error)=>error});}>"Fork saved source and load exact settings"</button>
                <button class="button secondary" on:click=move|_|{let record=bound.get_value();status.set(match on_gel.run((record.name.clone(),record.recipe.clone())){Ok(())=>"Added the exact recipe to predicted gel. Prediction must freshly reconstruct this source and catalog geometry; cached fragments are not lane authority.".into(),Err(error)=>error});}>"Add recipe to predicted gel"</button>
            </div>
            <p class="helper">"Use settings requires the same exact active source event. The explicit fork action can copy that saved version and apply its exact catalog pin and limits to the new fork; a fresh digest run is still required. Historical gel lanes retain their saved source version."</p>
            <details><summary>"Saved recognition-site evidence"</summary>
                {move||{
                    let count=bound.with_value(|r|r.report.inspection.sites.len());let pages=count.div_ceil(100).max(1);let page=site_page.get().min(pages-1);
                    let sites=bound.with_value(|r|r.report.inspection.sites.iter().skip(page*100).take(100).map(|site|view!{<article><strong>{site.enzyme.clone()}</strong><p>{format!("{} · {} · top boundary {} · bottom boundary {}",span_label(site),if site.certain{"certain recognition"}else{"possible recognition"},cut_label(site.top_cut,site.top_boundary),cut_label(site.bottom_cut,site.bottom_boundary))}</p>{site.diagnostics.iter().map(|diagnostic|view!{<p class="digest-warning">{site_diagnostic(*diagnostic)}</p>}).collect_view()}</article>}).collect_view());
                    view!{<div>{sites}<p class="helper">{format!("{} saved sites · page {} of {}",count,page+1,pages)}</p><div class="digest-actions"><button class="text-tool" disabled={page==0} on:click=move|_|{site_page.set(page.saturating_sub(1));}>"Previous saved sites"</button><button class="text-tool" disabled={page+1>=pages} on:click=move|_|{site_page.set((page+1).min(pages-1));}>"Next saved sites"</button></div></div>}
                }}
            </details>
            <details><summary>"Saved fragment evidence"</summary>
                {move||{
                    let count=bound.with_value(|r|r.report.inspection.fragments.len());let pages=count.div_ceil(25).max(1);let page=fragment_page.get().min(pages-1);
                    let fragments=bound.with_value(|r|r.report.inspection.fragments.iter().skip(page*25).take(25).map(|fragment|{
                        let ends=fragment.duplex.ends.as_ref().map(|(a,b)|format!("Left: {} · right: {}",end_label(a),end_label(b))).unwrap_or_else(||"No free ends on circular substrate".into());
                        view!{<article><strong>{format!("Saved fragment {} · {:?}",fragment.index+1,fragment.duplex.topology)}</strong><p>{format!("{} nt top / {} nt bottom · original reference start {}",fragment.duplex.top_length,fragment.duplex.bottom_length,fragment.top_start+1)}</p><p>{ends}</p></article>}
                    }).collect_view());
                    view!{<div>{fragments}<p class="helper">{format!("{} saved fragments · page {} of {}",count,page+1,pages)}</p><div class="digest-actions"><button class="text-tool" disabled={page==0} on:click=move|_|{fragment_page.set(page.saturating_sub(1));}>"Previous saved fragments"</button><button class="text-tool" disabled={page+1>=pages} on:click=move|_|{fragment_page.set((page+1).min(pages-1));}>"Next saved fragments"</button></div></div>}
                }}
            </details>
            <button class="text-tool" on:click=move|_|{
                let mut bytes=Json(Vec::new());
                let result=bound.with_value(|record|serde_json::to_writer_pretty(&mut bytes,record.as_ref()).map_err(|e|e.to_string())).and_then(|()|download("saved-digest.json","application/json",&bytes.0));
                crate::app::complete_download(result,move|result|{if status.try_get_untracked().is_some(){if let Err(error)=result{status.set(error);}}});
            }>"Download historical digest recipe and evidence"</button>
            <p class="helper" role="status">{move||status.get()}</p>
        </section>
    }.into_any()
}

/// A retained recipe is descriptive evidence; it cannot restore a current digest.
#[component]
pub fn DigestProvenance(molecule: Memo<Molecule>) -> impl IntoView {
    let retained = Memo::new(move |_| {
        molecule.with(|m| {
            m.metadata.contains_key(catalog_digest::RECIPE_KEY)
                || m.metadata.contains_key(catalog_digest::IMPORTED_RECIPE_KEY)
                || m.metadata.contains_key(restriction_inspection::RECIPE_KEY)
                || m.metadata
                    .contains_key("source:restriction_inspection:recipe")
        })
    });
    let recipe = Memo::new(move |_| {
        molecule.with(|m| {
            restriction_inspection::unverified_recipe(m)
                .map(|r| {
                    r.and_then(|_| m.metadata.get(restriction_inspection::RECIPE_KEY).cloned())
                })
                .map_err(|e| e.to_string())
        })
    });
    let catalog_recipe = Memo::new(move |_| -> Result<Option<String>, String> {
        molecule.with(|m| {
            let Some(raw) = m
                .metadata
                .get(catalog_digest::RECIPE_KEY)
                .or_else(|| m.metadata.get(catalog_digest::IMPORTED_RECIPE_KEY))
            else {
                return Ok(None);
            };
            catalog_digest::unverified_selection(m).map_err(|e| e.to_string())?;
            Ok(Some(raw.clone()))
        })
    });
    let status = RwSignal::new(String::new());
    view! {<Show when=move||retained.get()><details class="digest-provenance"><summary>"Saved virtual-digest provenance · unverified"</summary><p>"Saved recipes and end summaries do not restore a current digest or prove experimental cleavage. Open the original source from workspace history and run a new digest before reusing its predictions. Text exports may retain only bounded producer comments; keep the native workspace for exact source history."</p><Show when=move||recipe.get().is_ok_and(|r|r.is_some())><button class="button secondary" on:click=move|_|match recipe.get_untracked(){Ok(Some(raw))=>crate::app::download_with_completion("digest-retained-recipe.json","application/json",raw.as_bytes(),move|result|{if status.try_get_untracked().is_none(){return;}if let Err(e)=result {status.set(e)}}),Err(e)=>status.set(e),_=>status.set("A complete native producer recipe is not present.".into())}>"Download retained digest recipe"</button></Show><Show when=move||catalog_recipe.get().is_ok_and(|r|r.is_some())><p class="helper">"The retained catalog edition, content hash, ordered enzyme identities and geometry are unverified producer claims. They do not select current enzymes or restore a current result."</p><button class="button secondary" on:click=move|_|match catalog_recipe.get_untracked(){Ok(Some(raw))=>crate::app::download_with_completion("digest-retained-catalog.json","application/json",raw.as_bytes(),move|result|{if status.try_get_untracked().is_none(){return;}if let Err(e)=result {status.set(e)}}),Err(e)=>status.set(e),_=>{}}>"Download retained catalog selection"</button></Show><p class="digest-warning">{move||catalog_recipe.get().err().unwrap_or_default()}</p><p class="helper" role="status">{move||status.get()}</p></details></Show>}
}

fn catalog_temperature_label(t: i32) -> String {
    if t % 1000 == 0 {
        format!("{} °C", t / 1000)
    } else {
        format!("{}.{:03} °C", t / 1000, t % 1000)
    }
}
fn catalog_condition_label(condition: &catalog::Condition) -> String {
    match condition {
        catalog::Condition::IncubationTemperature => "Incubation temperature".into(),
        catalog::Condition::BufferActivity {
            buffer,
            temperature_milli_c,
        } => {
            format!(
                "{buffer} activity at {}",
                catalog_temperature_label(*temperature_milli_c)
            )
        }
        catalog::Condition::Methylation { system } => format!(
            "{} methylation sensitivity",
            match system {
                catalog::MethylationSystem::Dam => "Dam",
                catalog::MethylationSystem::Dcm => "Dcm",
                catalog::MethylationSystem::CpG => "CpG",
            }
        ),
    }
}
fn catalog_value_label(value: &catalog::FactValue) -> String {
    match value {
        catalog::FactValue::TemperatureMilliC(t) => catalog_temperature_label(*t),
        catalog::FactValue::Percent(p) => format!("{p}% relative activity"),
        catalog::FactValue::LessThanPercent(p) => format!("<{p}% relative activity"),
        catalog::FactValue::Unknown(reason) => format!("Unknown: {reason}"),
        catalog::FactValue::Methylation(sensitivity) => match sensitivity {
            catalog::Sensitivity::NotSensitive => "Not sensitive under the documentary profile",
            catalog::Sensitivity::Blocked => "Blocked",
            catalog::Sensitivity::Impaired => "Impaired",
            catalog::Sensitivity::BlockedOverlapping => "Blocked by overlapping methylation",
            catalog::Sensitivity::ImpairedOverlapping => "Impaired by overlapping methylation",
            catalog::Sensitivity::BlockedSomeOverlappingCombinations => {
                "Blocked in some overlapping combinations"
            }
            catalog::Sensitivity::ImpairedSomeOverlappingCombinations => {
                "Impaired in some overlapping combinations"
            }
        }
        .into(),
    }
}
fn catalog_citations(
    catalog: &catalog::Catalog,
    citations: &[catalog::Citation],
) -> impl IntoView + use<> {
    citations.iter().filter_map(|citation| {
        let source = catalog.edition().sources.iter().find(|s| s.id == citation.source)?;
        Some(view!{<li><a href=source.url.clone() target="_blank" rel="noopener noreferrer">{format!("{} · {}",source.publisher,source.edition)}</a><span>{format!(" — {} · reviewed {}",citation.locator,source.retrieved_on)}</span></li>})
    }).collect_view()
}
fn catalog_evidence(catalog: &catalog::Catalog, names: &[String]) -> impl IntoView + use<> {
    catalog.edition().enzymes.iter().filter(|enzyme| names.contains(&enzyme.name)).map(|enzyme| {
        let forms = catalog.edition().formulations.iter().filter(|f|f.enzyme==enzyme.id).collect::<Vec<_>>();
        let geometry = &enzyme.geometry;
        let form_views = forms.iter().map(|formulation| {
            let mut conditions = formulation.evidence.iter().map(|v|v.condition.clone()).collect::<Vec<_>>();
            conditions.sort();conditions.dedup();
            let observations = conditions.into_iter().map(|condition| {
                let mut contexts = formulation.evidence.iter().filter(|v|v.condition==condition).map(|v|v.context.clone()).collect::<Vec<_>>();
                contexts.sort();contexts.dedup();
                let context_notice = if contexts.len()>1 { "Multiple documentary contexts: review each separately." } else { "" };
                let context_views = contexts.into_iter().map(|context| {
                    let (state, value) = match catalog.resolve(&formulation.id, &condition, Some(&context)) {
                        Ok(catalog::Resolution::Known { value, .. }) => ("Known documentary observation",catalog_value_label(value)),
                        Ok(catalog::Resolution::Unknown(_)) => ("Unknown", "No known value in this context.".into()),
                        Ok(catalog::Resolution::Conflict(_)) => ("Conflicting evidence", "The observations disagree; no single value is selected.".into()),
                        Ok(catalog::Resolution::ContextRequired(_)) => ("Context required", "Choose an exact documentary context.".into()),
                        Err(error) => ("Unavailable",error.to_string()),
                    };
                    let evidence = formulation.evidence.iter().filter(|v|v.condition==condition&&v.context==context).map(|fact| {
                        view!{<div class="digest-catalog-observation"><p>{catalog_value_label(&fact.value)}</p>{fact.qualifiers.iter().map(|q|view!{<p class="digest-warning">{q.clone()}</p>}).collect_view()}<ul class="digest-catalog-citations">{catalog_citations(catalog,&fact.citations)}</ul></div>}
                    }).collect_view();
                    view!{<div class="digest-catalog-context" data-evidence-state=state><p><strong>{state}</strong>{format!(" · {value}")}</p><p class="helper">{format!("Context: {context}")}</p>{evidence}</div>}
                }).collect_view();
                view!{<section class="digest-catalog-fact"><h5>{catalog_condition_label(&condition)}</h5><p class="helper">{context_notice}</p>{context_views}</section>}
            }).collect_view();
            view!{<details class="digest-formulation" data-formulation-id=formulation.id.clone()><summary>{format!("{} · {} {}",formulation.name,formulation.vendor,formulation.product_code)}</summary><p class="helper">{formulation.formulation_scope.clone()}</p><p class="helper">"Documentary values apply to this formulation and context. They do not establish cleavage, buffer compatibility or reaction success for this digest."</p>{observations}</details>}
        }).collect_view();
        view!{<section class="digest-catalog-identity" data-catalog-enzyme=enzyme.id.clone()><h4>{enzyme.name.clone()}</h4><p>{format!("{} · strand-cut boundaries {} / {} relative to recognition start",geometry.recognition,geometry.top_cut,geometry.bottom_cut)}</p><ul class="digest-catalog-citations">{catalog_citations(catalog,&enzyme.citations)}</ul>{if forms.is_empty(){Some(view!{<p class="digest-catalog-unknown">"Unknown: no vendor formulation or condition evidence has been curated for this identity. Geometry alone does not supply reaction conditions."</p>})}else{None}}{form_views}</section>}
    }).collect_view()
}

fn cut_boundary_label(
    boundary: Option<&collagene_core::digest_cut_selection::Boundary>,
    report: &RestrictionInspection,
    label: &str,
) -> String {
    use collagene_core::digest_cut_selection::Boundary;
    match boundary {
        None => format!("{label} boundary: not selected"),
        Some(Boundary::SourceStart) => format!("{label}: source start · top/bottom boundary 0"),
        Some(Boundary::SourceEnd) => format!(
            "{label}: source end · top/bottom boundary {}",
            report.source_length
        ),
        Some(Boundary::Cut { site_index }) => report
            .sites
            .get(*site_index)
            .map(|s| {
                format!(
                    "{label}: {} at {} · top/bottom cuts {} / {}{}",
                    s.enzyme,
                    span_label(s),
                    cut_label(s.top_cut, s.top_boundary),
                    cut_label(s.bottom_cut, s.bottom_boundary),
                    if s.certain {
                        ""
                    } else {
                        " · possible recognition"
                    }
                )
            })
            .unwrap_or_else(|| format!("{label}: site no longer available; choose again")),
    }
}
