//! Inspect a PCR template without retargeting the independent design target.
use crate::pcr_overlay::Overlay;
use collagene_model::Molecule;
use leptos::prelude::*;
use std::sync::Arc;

#[derive(Clone)]
pub struct Snapshot {
    pub owner: String,
    pub overlay: Arc<Overlay>,
    pub molecule: Arc<Molecule>,
    pub spans: Vec<(usize, usize)>,
}

#[component]
pub fn PcrSourceView(
    snapshot: Snapshot,
    current: Memo<bool>,
    on_back: Callback<()>,
) -> impl IntoView {
    let identity = format!(
        "{}:{}:{}",
        snapshot.overlay.epoch,
        snapshot.overlay.source.document_id,
        snapshot.overlay.source.event_id
    );
    let overlay = StoredValue::new(snapshot.overlay.clone());
    let visible = Memo::new(move |_| current.get().then(|| overlay.get_value()));
    let status = RwSignal::new(String::new());
    view! {
        <section>
            <h3>{snapshot.molecule.name.clone()}</h3>
            <p class="helper">"PCR template version · read-only. The independent target and source document heads are unchanged."</p>
            <p class="helper">{snapshot.overlay.summary.clone()}</p>
            <button class="button secondary" on:click=move|_|on_back.run(())>"Back to PCR review"</button>
            <p class="helper" role="status">{move||if current.get(){status.get()}else{"PCR settings or source context changed; select the product again.".into()}}</p>
            <crate::assembly_product_view::ReferenceCanvas product=snapshot.molecule.as_ref().clone()
                identity current initial_spans=snapshot.spans pcr_overlay=Some(visible)
                on_message=Callback::new(move|text|status.set(text))/>
        </section>
    }.into_any()
}
