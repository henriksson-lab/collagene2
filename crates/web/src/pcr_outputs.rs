//! Explicit selection and output review for exact-PCR reference previews.
//! Reaction reconstruction remains in the worker; this view never invents a target.
use crate::pcr_overlay::{self, Overlay, Publication, Role};
use collagene_core::{
    Molecule,
    analysis::Primer,
    oligo_formats::{self, Delimiter, TableExportMode},
    run_lineage::EventBinding,
    verification::PcrDocument,
};
use leptos::prelude::*;
use std::sync::{
    Arc,
    atomic::{AtomicU32, Ordering},
};
static NEXT_OVERLAY_OWNER: AtomicU32 = AtomicU32::new(1);

#[derive(Clone)]
pub struct Review {
    pub source: EventBinding,
    pub origin: Option<crate::run_origin_model::PinnedOrigin>,
    pub epoch: u64,
    pub generation: u32,
    pub inputs: (String, String, String, String, String, String),
    pub forward: Primer,
    pub reverse: Primer,
    pub products: Arc<Vec<PcrDocument>>,
}

#[derive(Clone)]
pub struct Choice {
    pub review: Review,
    pub product_index: usize,
    pub name: String,
    pub save_product: bool,
    pub save_primers: bool,
}

fn valid_name(name: &str) -> bool {
    !name.trim().is_empty()
        && name.trim() == name
        && name.len() <= 256
        && !name.chars().any(char::is_control)
}

#[component]
pub fn PcrOutputsPanel(
    review: Review,
    on_create: Callback<Choice>,
    guard: Callback<Review, bool>,
    busy: RwSignal<bool>,
    source_molecule: Memo<Molecule>,
    on_overlay: Callback<Publication>,
    on_navigate: Callback<(EventBinding, u64, Vec<(usize, usize)>), Result<(), String>>,
) -> impl IntoView {
    let saved_material = review.origin.is_some();
    let count = review.products.len();
    let selected = RwSignal::new((count == 1).then_some(0usize));
    let initial_name = if count == 1 {
        review.products[0].molecule.name.clone()
    } else {
        String::new()
    };
    let name = RwSignal::new(initial_name);
    let save_product = RwSignal::new(true);
    let save_primers = RwSignal::new(false);
    let page = RwSignal::new(0usize);
    let sequence_page = RwSignal::new(0usize);
    let status = RwSignal::new(String::new());
    let bound = StoredValue::new(review);
    let owner = StoredValue::new(format!(
        "exact-{}",
        NEXT_OVERLAY_OWNER.fetch_add(1, Ordering::Relaxed)
    ));
    let projected = RwSignal::new(None::<Arc<Overlay>>);
    let project = Callback::new(move |()| -> Result<Arc<Overlay>, String> {
        let review = bound.get_value();
        if !guard.run(review.clone()) {
            return Err("The template or primer pair changed; compute PCR products again.".into());
        }
        let index = selected
            .get_untracked()
            .ok_or("Choose a reported binding geometry first.")?;
        let product = review
            .products
            .get(index)
            .ok_or("The selected PCR geometry is unavailable.")?;
        pcr_overlay::exact(
            &source_molecule.get_untracked(),
            &review.source,
            review.epoch,
            &review.forward,
            &review.reverse,
            product,
        )
        .map(Arc::new)
    });
    Effect::new(move |_| {
        let index = selected.get();
        let _ = source_molecule.get();
        let overlay = if index.is_some() {
            match project.run(()) {
                Ok(overlay) => Some(overlay),
                Err(error) => {
                    status.set(error);
                    None
                }
            }
        } else {
            None
        };
        projected.set(overlay.clone());
        on_overlay.run(Publication {
            owner: owner.get_value(),
            overlay,
        });
    });
    on_cleanup(move || {
        on_overlay.run(Publication {
            owner: owner.get_value(),
            overlay: None,
        })
    });
    let navigate = Callback::new(move |role: Role| {
        let result = project.run(()).and_then(|overlay| {
            let spans = overlay
                .spans
                .iter()
                .filter(|span| span.role == role)
                .map(|span| (span.start, span.end))
                .collect::<Vec<_>>();
            if spans.is_empty() {
                return Err("The selected geometry has no mapped bases for this role.".into());
            }
            on_overlay.run(Publication {
                owner: owner.get_value(),
                overlay: Some(overlay.clone()),
            });
            on_navigate.run((overlay.source.clone(), overlay.epoch, spans))
        });
        if let Err(error) = result {
            status.set(error);
        }
    });
    let export = Callback::new(move |delimiter: Delimiter| {
        let review = bound.get_value();
        if busy.get_untracked() || !guard.run(review.clone()) {
            status.set("The template or primer pair changed; compute PCR products again.".into());
            return;
        }
        let result = collagene_core::pcr_outputs::primer_table(&review.forward, &review.reverse)
            .and_then(|table| {
                oligo_formats::export_oligos(&table, delimiter, TableExportMode::SpreadsheetSafe)
            });
        match result {
            Ok(report) => {
                let (filename, mime) = match delimiter {
                    Delimiter::Comma => ("pcr-primers.csv", "text/csv;charset=utf-8"),
                    Delimiter::Tab => {
                        ("pcr-primers.tsv", "text/tab-separated-values;charset=utf-8")
                    }
                };
                if !guard.run(review) {
                    status.set("The PCR inputs changed; no primer sheet downloaded.".into());
                    return;
                }
                crate::app::download_with_completion(
                    filename,
                    mime,
                    report.text.as_bytes(),
                    move |result| {
                        if status.try_get_untracked().is_none() {
                            return;
                        }
                        match result {
                            Ok(()) => status.set(report.warnings.join(" ")),
                            Err(error) => status.set(error),
                        }
                    },
                )
            }
            Err(error) => status.set(error.to_string()),
        }
    });
    let create = move |_| {
        let review = bound.get_value();
        if busy.get_untracked() || !guard.run(review.clone()) {
            status.set("The template or primer pair changed; compute PCR products again.".into());
            return;
        }
        let Some(product_index) = selected.get_untracked() else {
            status.set("Choose a reported binding geometry before creating outputs.".into());
            return;
        };
        let output_name = name.get_untracked();
        if !valid_name(&output_name) {
            status.set(
                "Use a trimmed nonempty output name, at most 256 bytes, with no control characters.".into(),
            );
            return;
        }
        if !save_product.get_untracked() && !save_primers.get_untracked() {
            status.set("Select at least one output to create.".into());
            return;
        }
        status.set(String::new());
        on_create.run(Choice {
            review,
            product_index,
            name: output_name,
            save_product: save_product.get_untracked(),
            save_primers: save_primers.get_untracked(),
        });
    };
    view! {
        <section class="pcr-output-review">
            <h3>"Review PCR products"</h3>
            {saved_material.then(||view!{<p class="helper">"The saved material was reconstructed before PCR geometry was calculated. Opened outputs remain documentary reference previews; use a saved-origin PCR preparation in assembly and record its verified run to create reusable material."</p>})}
            <p class="helper">{format!("{count} exact-annealing binding geometries. Choose the reference preview to retain; alternative geometries remain separate. Physical end chemistry, target verification and experimental success are not established.")}</p>
            <details><summary>"Original template event"</summary><p>{bound.with_value(|r|format!("Document {} · event {} · revision {}",r.source.document_id,r.source.event_id,r.source.revision))}</p></details>
            <details open><summary>"Declared oligos · ordered 5′→3′"</summary>
                {bound.with_value(|r|[r.forward.clone(),r.reverse.clone()]).into_iter().enumerate().map(|(index,primer)|view!{
                    <article><strong>{format!("{}: {}",if index==0 {"Forward"} else {"Reverse"},primer.name)}</strong>
                        <p class="sequence-preview">{primer.ordered_sequence()}</p>
                        <p class="helper">{format!("5′ tail: {} · annealing domain: {}",if primer.tail.is_empty(){"(none)"}else{&primer.tail},primer.annealing.as_str())}</p>
                    </article>
                }).collect_view()}
                <p class="helper">"The search can use either oligo in either binding role, including single-oligo amplification. Product rows show the roles actually found. Order sheets retain both declared oligos; they are reference data, with no quantities or modifications inferred."</p>
                <div class="method-row-actions"><button class="text-tool" disabled=move||busy.get() on:click=move|_|export.run(Delimiter::Comma)>"Primer CSV"</button><button class="text-tool" disabled=move||busy.get() on:click=move|_|export.run(Delimiter::Tab)>"Primer TSV"</button></div>
            </details>
            {move || {
                let current = page.get().min(count.saturating_sub(1)/10);
                let pages = count.div_ceil(10).max(1);
                let rows=bound.with_value(|r|r.products.iter().enumerate().skip(current*10).take(10).map(|(index,product)|{
                    let geometry=&product.geometry;
                    let default_name=product.molecule.name.clone();
                    view!{
                        <article class="mutation-product">
                            <label><input type="radio" name="pcr-output-geometry" prop:checked=move||selected.get()==Some(index) disabled=move||busy.get() on:change=move|_|{selected.set(Some(index));name.set(default_name.clone());sequence_page.set(0);status.set(String::new());}/>{format!("Product {} · {} bp",index+1,product.molecule.sequence.len())}</label>
                            <p class="helper">{format!("Forward-role {} at base {} · reverse-role {} at base {}{} · {} projected annotations",geometry.forward_primer,geometry.forward_start+1,geometry.reverse_primer,geometry.reverse_start+1,if geometry.crosses_origin{" · crosses template origin"}else{""},product.molecule.features.len())}</p>
                        </article>
                    }
                }).collect_view());
                view!{<div>{rows}</div><div class="method-row-actions"><button class="text-tool" disabled=move||busy.get()||current==0 on:click=move|_|page.set(current.saturating_sub(1))>"Previous products"</button><span>{format!("Page {} of {pages}",current+1)}</span><button class="text-tool" disabled=move||{busy.get()||current+1>=pages} on:click=move|_|page.set((current+1).min(pages-1))>"Next products"</button></div>}
            }}
            {move||selected.get().and_then(|index|bound.with_value(|r|r.products.get(index).map(|product|{
                let length=product.molecule.sequence.len();
                let pages=length.div_ceil(500).max(1);
                let current=sequence_page.get().min(pages-1);
                let start=(current*500).min(length);
                let end=(start+500).min(length);
                let sequence=product.molecule.sequence.as_str()[start..end].to_string();
                view!{<details open><summary>"Selected reference preview"</summary><p class="helper">{format!("Product bases {}–{} of {length}; page {} of {pages}",start+1,end,current+1)}</p><pre class="sequence-preview">{sequence}</pre><div class="method-row-actions"><button class="text-tool" disabled=move||busy.get()||current==0 on:click=move|_|sequence_page.set(current.saturating_sub(1))>"Previous bases"</button><button class="text-tool" disabled=move||{busy.get()||current+1>=pages} on:click=move|_|sequence_page.set((current+1).min(pages-1))>"Next bases"</button></div>{product.findings.iter().cloned().map(|finding|view!{<p class="helper">{finding}</p>}).collect_view()}</details>}
            })))}
            {move||projected.get().map(|overlay|view!{
                <section class="pcr-placement-review"><h4>"Selected template placement"</h4>
                    <p class="helper">{overlay.summary.clone()}</p>
                    <div class="method-row-actions"><button class="text-tool" disabled=move||{busy.get()} on:click=move|_|{navigate.run(Role::Amplified);}>"View template interval"</button><button class="text-tool" disabled=move||{busy.get()} on:click=move|_|{navigate.run(Role::Forward);}>"View forward-role primer"</button><button class="text-tool" disabled=move||{busy.get()} on:click=move|_|{navigate.run(Role::Reverse);}>"View reverse-role primer"</button></div>
                    <p class="helper">"Only the selected reported geometry is overlaid. Tail bases are synthetic additions and have no template coordinates."</p>
                    {overlay.tails.iter().cloned().map(|tail|view!{<p class="sequence-preview">{tail}</p>}).collect_view()}
                    {overlay.warnings.iter().cloned().map(|warning|view!{<p class="helper">{warning}</p>}).collect_view()}
                </section>
            })}
            <label>"Output name"<input maxlength="256" disabled=move||busy.get() prop:value=move||name.get() on:input=move|event|name.set(event_target_value(&event))/></label>
            <label><input type="checkbox" prop:checked=move||save_product.get() disabled=move||busy.get() on:change=move|event|save_product.set(event_target_checked(&event))/>"Create selected PCR reference preview"</label>
            <label><input type="checkbox" prop:checked=move||save_primers.get() disabled=move||busy.get() on:change=move|event|save_primers.set(event_target_checked(&event))/>"Save declared primer pair as an oligo collection"</label>
            <button class="button primary full" disabled=move||busy.get()||selected.get().is_none()||!valid_name(&name.get())||(!save_product.get()&&!save_primers.get()) on:click=create>"Create selected outputs"</button>
            <p class="helper" role="status">{move||status.get()}</p>
        </section>
    }.into_any()
}
