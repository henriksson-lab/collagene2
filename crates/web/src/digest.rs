//! Exact-source virtual restriction inspection. Scientific work stays in the worker.
use crate::{app::download, worker::CachedWorker};
use collagene_core::{
    Alphabet, Molecule, Topology, analysis,
    catalog_digest::{self, CatalogSelection, catalog},
    chemistry::{End, Polarity},
    formats,
    protocol::{AnalysisJob, AnalysisResult, Request},
    restriction_history,
    restriction_inspection::{
        self, InspectionLimit, InspectionOptions, InspectionSite, RestrictionInspection,
        SiteDiagnostic,
    },
    run_lineage::EventBinding,
};
use leptos::prelude::*;
use serde::{Deserialize, Serialize};
use std::{
    rc::Rc,
    sync::{
        Arc,
        atomic::{AtomicBool, Ordering},
    },
};
use wasm_bindgen_futures::spawn_local;
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct SourceInfo {
    pub binding: EventBinding,
    pub name: String,
    pub length: usize,
    pub alphabet: Alphabet,
    pub topology: Topology,
}
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
    source: Memo<SourceInfo>,
    get_molecule: Callback<(), Molecule>,
    epoch: RwSignal<u64>,
    on_navigate: Callback<Navigation>,
    on_open: Callback<crate::digest_product::DigestOpenRequest, Result<(), String>>,
    on_methylation: Callback<()>,
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
    let names = RwSignal::new(vec!["EcoRI".to_string()]);
    let filter = RwSignal::new(String::new());
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
        let selection = catalog_digest::select(&catalog, &ids).map_err(|e| e.to_string())?;
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
    let worker = StoredValue::new_local(None::<Rc<CachedWorker>>);
    let site_page = RwSignal::new(0usize);
    let fragment_page = RwSignal::new(0usize);
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
        let _ = (source.get(), inputs.get(), epoch.get());
        stop.run(());
        bound.set(None);
        site_page.set(0);
        fragment_page.set(0);
        status.set(
            "Choose enzymes and run a virtual digest. Results refresh only after an explicit run."
                .into(),
        );
    });
    let current = Callback::new(move |data: Arc<Bound>| {
        data.epoch == epoch.get_untracked()
            && source.get_untracked() == data.source
            && inputs.get_untracked() == Ok(data.inputs.clone())
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
                catalog_digest::source_preflight(
                    &bound_catalog,
                    &selection,
                    &role,
                    &molecule,
                    &options,
                )
                .map_err(|e| e.to_string())?;
                if molecule.revision() != source.binding.revision {
                    return Err("The source document/event changed. Run the digest again.".into());
                }
                let revision = catalog_digest::batch_revision(
                    &bound_catalog,
                    &selection,
                    &role,
                    &molecule,
                    &options,
                )
                .map_err(|e| e.to_string())?;
                stop.run(());
                let w = Rc::new(CachedWorker::new()?);
                worker.set_value(Some(w.clone()));
                busy.set(true);
                bound.set(None);
                site_page.set(0);
                fragment_page.set(0);
                Ok((
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
                        job: AnalysisJob::InspectCatalogRestriction {
                            source_role: role.clone(),
                            molecule: molecule.clone(),
                            selection: selection.clone(),
                            options: options.clone(),
                        },
                    })
                    .await;
                if !alive.load(Ordering::Relaxed)
                    || generation.try_get_untracked() != Some(expected_generation)
                    || epoch.get_untracked() != expected_epoch
                    || source.get_untracked() != source_snapshot
                    || inputs.get_untracked() != Ok(settings.clone())
                {
                    return;
                }
                w.close();
                worker.set_value(None);
                busy.set(false);
                let report = match response {
                    Ok(response) if response.id == 1 && response.revision == revision => {
                        match response.result {
                            AnalysisResult::CatalogRestrictionInspection(report) => report,
                            AnalysisResult::Error(e) => {
                                status.set(e);
                                return;
                            }
                            _ => {
                                status.set("Unexpected restriction inspection response.".into());
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
                let bound_catalog = match catalog.get_value() {
                    Ok(c) => c,
                    Err(e) => {
                        status.set(e);
                        return;
                    }
                };
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
                let catalog_revision = report.batch_revision.clone();
                let catalog_engine_version = report.engine_version.clone();
                let report = report.inspection;
                status.set(if report.products_determined{format!("{} recognition sites · {} virtual fragment(s), assuming complete cleavage.",report.sites.len(),report.fragments.len())}else if report.scan_complete{"Recognition scan finished; a unique complete product set is not determined. Review the evidence below.".into()}else{"Recognition scan is incomplete. Unsearched sequence cannot establish absence or a complete product set.".into()});
                bound.set(Some(Arc::new(Bound {
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
            selection: data.inputs.selection.clone(),
            options: data.inputs.limits.options(),
            fragment_index: index,
            epoch: data.epoch,
        };
        if let Err(e) = on_open.run(request) {
            status.set(e);
        }
    });
    let export_product = {
        let alive = alive.clone();
        Callback::new(move |index: usize| {
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
                    Ok(export) => match download(
                        &format!("digest-fragment-{}.gb", index + 1),
                        "text/plain",
                        export.text.as_bytes(),
                    ) {
                        Ok(()) => {
                            let mut findings = product.projection_findings;
                            findings.extend(export.warnings);
                            status.set(format!("Exported freshly reconstructed fragment {}. GenBank producer/end comments are unverified provenance; keep the native workspace for exact source history. {}",index+1,findings.join(" ")));
                        }
                        Err(e) => status.set(e),
                    },
                    Err(e) => status.set(e.to_string()),
                }
            });
        })
    };
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
            report: &'a RestrictionInspection,
            catalog_pin: &'a CatalogSelection,
            catalog_engine_version: &'a str,
            catalog_batch_revision: &'a str,
        }
        let artifact = Evidence {
            kind: "virtual-digest-review-v2",
            authority: "Sequence/geometry prediction under stated cleavage assumptions; saved evidence never restores current verification",
            source: &data.source,
            report: &data.report,
            catalog_pin: &data.inputs.selection,
            catalog_engine_version: &data.catalog_engine_version,
            catalog_batch_revision: &data.catalog_revision,
        };
        let mut bytes = Json(vec![]);
        match serde_json::to_writer_pretty(&mut bytes, &artifact) {
            Ok(()) => {
                if let Err(e) = download("virtual-digest.json", "application/json", &bytes.0) {
                    status.set(e)
                }
            }
            Err(e) => status.set(e.to_string()),
        }
    };
    view! {
        <section class="inspector-content digest-panel"><h2>"Virtual digest"</h2><p class="helper">"Choose enzymes to inspect recognition sites and predicted fragments. Products assume complete cleavage; methylation, reaction conditions and laboratory success are not assessed."</p>
        <label>"Find an enzyme"<input aria-label="Find digest enzyme" prop:value=move||filter.get() on:input=move|e|filter.set(event_target_value(&e)) placeholder="Name or recognition sequence"/></label>
        <div class="digest-enzyme-picker">{move||{let needle=filter.get().trim().to_ascii_lowercase();enzyme_choices.get_value().into_iter().filter(|e|e.name.to_ascii_lowercase().contains(&needle)||e.geometry.recognition.to_ascii_lowercase().contains(&needle)).map(|enzyme|{let name=enzyme.name.clone();let change=enzyme.name.clone();view!{<label class="checkbox-label digest-enzyme-choice"><input type="checkbox" data-enzyme=enzyme.name.clone() data-catalog-id=enzyme.id.clone() aria-label=format!("Digest enzyme {}",enzyme.name) prop:checked=move||names.with(|names|names.contains(&name)) on:change=move|e|{let checked=event_target_checked(&e);names.update(|names|{names.retain(|n|n!=&change);if checked{names.push(change.clone());}});}/><span><strong>{enzyme.name}</strong><small>{enzyme.geometry.recognition}</small></span></label>}}).collect_view()}}</div>
        <p class="helper">{move||format!("Selected: {}",names.get().join(", "))}</p>
        <details class="digest-catalog" data-catalog-edition=move||catalog.with_value(|c|c.as_ref().map(|c|c.edition().edition.clone()).unwrap_or_default()) data-catalog-hash=move||catalog.with_value(|c|c.as_ref().map(|c|c.content_hash().to_string()).unwrap_or_default())><summary>"Catalog and formulation evidence"</summary><p class="helper">{move||catalog.with_value(|c|match c {Ok(c)=>format!("{} enzyme identities · edition {} · content hash {}",c.edition().enzymes.len(),c.edition().edition,c.content_hash()),Err(e)=>format!("Catalog unavailable: {e}")})}</p><p class="helper">"Recognition and strand-cut geometry determine the virtual digest. Same-cut enzymes remain distinct identities; vendor formulations and their documentary conditions are reviewed separately and never change or validate predicted cleavage."</p>{move||catalog.with_value(|c|c.as_ref().ok().map(|c|catalog_evidence(c,&names.get())))}</details><div class="digest-actions"><button class="button primary" disabled=move||busy.get()||names.with(Vec::is_empty)||source.get().alphabet!=Alphabet::Dna on:click=run>"Run virtual digest"</button><button class="button secondary" disabled=move||!busy.get() on:click=move|_|{stop.run(());status.set("Digest request cancelled; no pending result or product was accepted.".into());}>"Cancel digest"</button></div>
        <p class="digest-status helper" role="status">{move||status.get()}</p>
        <details class="digest-options"><summary>"Assumptions and analysis limits"</summary><p class="helper">"The pinned catalog selection describes recognition and strand-cut geometry. Documentary vendor formulation evidence is separate and is not a reaction-compatibility assessment. Recognition can be certain while cuts lack sufficient substrate, and ambiguous sites can prevent a complete product set."</p><button class="text-tool" on:click=move|_|on_methylation.run(())>"Review methylation assumptions"</button><label>"Maximum comparisons"<input aria-label="Digest comparison limit" prop:value=move||comparisons.get() on:input=move|e|comparisons.set(event_target_value(&e))/></label><label>"Maximum recognition sites"<input aria-label="Digest site limit" prop:value=move||max_sites.get() on:input=move|e|max_sites.set(event_target_value(&e))/></label><label>"Maximum fragments"<input aria-label="Digest fragment limit" prop:value=move||max_fragments.get() on:input=move|e|max_fragments.set(event_target_value(&e))/></label><label>"Maximum report bytes"<input aria-label="Digest report byte limit" prop:value=move||report_bytes.get() on:input=move|e|report_bytes.set(event_target_value(&e))/></label></details>
        {move||bound.get().map(|data|{let sites=data.report.sites.len();let fragments=data.report.fragments.len();view!{<div class="digest-results" data-scan-complete=data.report.scan_complete.to_string() data-products-determined=data.report.products_determined.to_string()><h3>"Recognition and products"</h3><p>{format!("Recognition scan: {}. Product set: {}.",if data.report.scan_complete{"finished"}else{"incomplete"},if data.report.products_determined{"determined under complete-cleavage assumptions"}else{"not determined"})}</p>{data.report.limits.iter().map(|limit|view!{<p class="digest-warning">{format!("{} reached. Incomplete evidence cannot establish absence or complete products.",limit_label(*limit))}</p>}).collect_view()}{data.report.diagnostics.iter().map(|d|view!{<p class="digest-warning">{d.clone()}</p>}).collect_view()}
        <h3>{format!("{sites} recognition sites")}</h3><p class="helper">"Recognition ranges are original source bases, 1-based inclusive. Strand cuts are boundaries measured from the source origin; 0 is before the first base."</p>{data.report.sites.iter().skip(site_page.get()*25).take(25).map(|site|{let index=site.index;view!{<article class="digest-site"><button class="result-entry" on:click=move|_|navigate.run(index)>{format!("{} · {} · {}{}",site.enzyme,span_label(site),if site.reverse{"reverse"}else{"forward"},if site.certain{""}else{" · possible"})}</button><p>{format!("Top / bottom cut boundaries: {} / {}",cut_label(site.top_cut,site.top_boundary),cut_label(site.bottom_cut,site.bottom_boundary))}</p>{site.diagnostics.iter().map(|d|view!{<p class="digest-warning">{site_diagnostic(*d)}</p>}).collect_view()}</article>}}).collect_view()}<div class="digest-actions"><button class="button secondary" disabled=move||site_page.get()==0 on:click=move|_|site_page.update(|p|*p=p.saturating_sub(1))>"Previous digest sites"</button><button class="button secondary" disabled={move||(site_page.get()+1)*25>=sites} on:click=move|_|site_page.update(|p|*p+=1)>"Next digest sites"</button></div>
        <h3>"Reference-strand fragment lengths"</h3><p class="fragment-sizes">{data.report.ref_top_lengths.as_ref().map(|lengths|lengths.iter().take(100).map(|n|format!("{n} nt")).collect::<Vec<_>>().join(" · ")).unwrap_or_else(||"Not determined from this recognition evidence.".into())}</p><p class="helper">"First 100 lengths shown, following the reference/top strand and sorted by size. The JSON evidence retains every determined length. Physical summaries below use ascending reference cut order; staggered ends can give different top and bottom lengths."</p>
        <h3>{format!("{fragments} virtual fragment(s)")}</h3>{data.report.fragments.iter().skip(fragment_page.get()*25).take(25).map(|fragment|{let index=fragment.index;let ends=fragment.duplex.ends.as_ref().map(|(a,b)|format!("Left: {} · right: {}",end_label(a),end_label(b))).unwrap_or_else(||"Uncut circular substrate: no free ends".into());view!{<article class="digest-fragment"><h4>{format!("Fragment {} · {}",index+1,if fragment.duplex.topology==Topology::Circular{"circular"}else{"linear"})}</h4><p>{format!("{} nt top / {} nt bottom · original reference start {}",fragment.duplex.top_length,fragment.duplex.bottom_length,fragment.top_start+1)}</p><p>{ends}</p><div class="digest-actions"><button class="button secondary" disabled=move||busy.get() on:click=move|_|open.run(index)>"Open fragment"</button><button class="button secondary" disabled=move||busy.get() on:click=move|_|export_product.run(index)>{move||if product_busy.get()==Some(index){"Reconstructing…"}else{"Export fragment GenBank"}}</button></div></article>}}).collect_view()}<div class="digest-actions"><button class="button secondary" disabled=move||fragment_page.get()==0 on:click=move|_|fragment_page.update(|p|*p=p.saturating_sub(1))>"Previous digest fragments"</button><button class="button secondary" disabled={move||(fragment_page.get()+1)*25>=fragments} on:click=move|_|fragment_page.update(|p|*p+=1)>"Next digest fragments"</button></div><button class="button secondary full" on:click=export_report>"Download digest evidence"</button><p class="helper">"Opening reconstructs the fragment again and retains its exact source event in the workspace. GenBank comments preserve descriptive, unverified producer/end information; native workspace files retain full source history."</p></div>}})}
        </section>
    }
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
    view! {<Show when=move||retained.get()><details class="digest-provenance"><summary>"Saved virtual-digest provenance · unverified"</summary><p>"Saved recipes and end summaries do not restore a current digest or prove experimental cleavage. Open the original source from workspace history and run a new digest before reusing its predictions. Text exports may retain only bounded producer comments; keep the native workspace for exact source history."</p><Show when=move||recipe.get().is_ok_and(|r|r.is_some())><button class="button secondary" on:click=move|_|match recipe.get_untracked(){Ok(Some(raw))=>if let Err(e)=download("digest-retained-recipe.json","application/json",raw.as_bytes()){status.set(e)},Err(e)=>status.set(e),_=>status.set("A complete native producer recipe is not present.".into())}>"Download retained digest recipe"</button></Show><Show when=move||catalog_recipe.get().is_ok_and(|r|r.is_some())><p class="helper">"The retained catalog edition, content hash, ordered enzyme identities and geometry are unverified producer claims. They do not select current enzymes or restore a current result."</p><button class="button secondary" on:click=move|_|match catalog_recipe.get_untracked(){Ok(Some(raw))=>if let Err(e)=download("digest-retained-catalog.json","application/json",raw.as_bytes()){status.set(e)},Err(e)=>status.set(e),_=>{}}>"Download retained catalog selection"</button></Show><p class="digest-warning">{move||catalog_recipe.get().err().unwrap_or_default()}</p><p class="helper" role="status">{move||status.get()}</p></details></Show>}
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
