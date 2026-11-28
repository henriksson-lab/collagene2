//! Read-only original coordinates for historical or origin-spanning comparisons.
use collagene_core::material::lineage::EventBinding;
use collagene_model::Molecule;
use leptos::prelude::*;
use std::sync::Arc;

#[derive(Clone)]
pub struct Snapshot {
    pub source: EventBinding,
    pub epoch: u64,
    pub molecule: Arc<Molecule>,
    pub spans: Vec<(usize, usize)>,
}

#[component]
pub fn PairwiseSourceView(
    snapshot: Snapshot,
    epoch: RwSignal<u64>,
    on_back: Callback<()>,
) -> impl IntoView {
    let current = Memo::new(move |_| epoch.get() == snapshot.epoch);
    let identity = format!(
        "comparison-version:{}:{}:{}:{:?}",
        snapshot.source.document_id,
        snapshot.source.event_id,
        snapshot.source.revision,
        snapshot.spans
    );
    let status = RwSignal::new(String::new());
    view! {
        <section class="comparison-source-view">
            <h3>{snapshot.molecule.name.clone()}</h3>
            <p class="helper">"Exact comparison source version · read-only. Coordinates refer to the original stored document. An empty interval marks an insertion boundary. No source head is changed."</p>
            <p style="overflow-wrap:anywhere">{format!("Original boundary intervals: {:?}",snapshot.spans)}</p>
            <details style="overflow-wrap:anywhere"><summary>"Source version"</summary>
                <p>{format!("Document {} · event {} · revision {}",snapshot.source.document_id,snapshot.source.event_id,snapshot.source.revision)}</p>
            </details>
            <button class="button secondary" on:click=move|_|on_back.run(())>"Back to comparison"</button>
            <p class="helper" role="status">{move||if current.get(){status.get()}else{"Workspace context changed; reopen the comparison.".into()}}</p>
            <crate::assembly::review::assembly_product_view::ReferenceCanvas product=snapshot.molecule.as_ref().clone()
                identity current initial_spans=snapshot.spans on_message=Callback::new(move|message|status.set(message))/>
        </section>
    }.into_any()
}
