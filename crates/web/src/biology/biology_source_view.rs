//! Historical translation origin is a reference view, not current protein validation.
use collagene_core::material::lineage::EventBinding;
use collagene_model::Molecule;
use leptos::prelude::*;

#[component]
pub fn TranslationSourceLink(
    molecule: Memo<Molecule>,
    on_open: Callback<EventBinding, Result<(), String>>,
) -> impl IntoView {
    let status = RwSignal::new(String::new());
    let source = Memo::new(move |_| {
        molecule.with(|m| {
            Some(EventBinding {
                document_id: m.metadata.get("translation:source_document")?.clone(),
                event_id: m.metadata.get("translation:source_event")?.clone(),
                revision: m.metadata.get("translation:source_revision")?.clone(),
            })
        })
    });
    let policy = Memo::new(move |_| {
        molecule.with(|m| {
        let policy=m.metadata.get("translation:policy").map(String::as_str).unwrap_or("Declared translation policy unavailable");
        let options=m.metadata.get("translation:options").filter(|text|text.len()<=4096)
            .and_then(|text|serde_json::from_str::<collagene_core::sequence::translation::TranslationOptions>(text).ok());
        match options {
            Some(options)=>format!("{policy}. Genetic code {} · skip {} bases at the oriented 5′ end · {} · {}.",
                options.genetic_code,options.frame,if options.reverse{"reverse complement of selected location"}else{"selected biological orientation"},
                match options.initiation {
                    collagene_core::sequence::translation::InitiationPolicy::Ordinary=>"ordinary codon translation",
                    collagene_core::sequence::translation::InitiationPolicy::RecognizedStartAsMethionine=>"recognized first initiator recoded as methionine",
                    collagene_core::sequence::translation::InitiationPolicy::RequireRecognizedStart=>"definite recognized first initiator required",
                }),
            None=>format!("{policy}. Exact options unavailable."),
        }
    })
    });
    let warnings = Memo::new(move |_| {
        molecule.with(|m| {
            m.metadata
                .get("translation:warnings")
                .filter(|text| text.len() <= 256 * 1024)
                .and_then(|text| serde_json::from_str::<Vec<String>>(text).ok())
                .unwrap_or_default()
                .into_iter()
                .take(64)
                .map(|text| text.chars().take(4096).collect::<String>())
                .collect::<Vec<_>>()
        })
    });
    view! {
        <Show when=move||source.with(Option::is_some)>
            <section class="inspector-content translation-source-link">
                <h3>"Translation origin"</h3>
                <p class="helper">{move||policy.get()}</p>
                <ul>{move||warnings.get().into_iter().map(|warning|view!{<li>{warning}</li>}).collect_view()}</ul>
                <p class="helper">"Open the exact retained nucleotide source version. This origin link does not verify later protein edits or establish expression."</p>
                <button class="button secondary full" on:click=move|_|{
                    if let Some(source)=source.get_untracked(){
                        status.set(match on_open.run(source){Ok(())=>String::new(),Err(error)=>error});
                    }
                }>"View original translation source"</button>
                <p class="helper" role="status">{move||status.get()}</p>
            </section>
        </Show>
    }.into_any()
}

#[component]
pub fn TranslationSourceView(
    snapshot: crate::viewer::alignment::pairwise_source_view::Snapshot,
    epoch: RwSignal<u64>,
    on_back: Callback<()>,
) -> impl IntoView {
    let current = Memo::new(move |_| epoch.get() == snapshot.epoch);
    let identity = format!(
        "translation-origin:{}:{}:{}",
        snapshot.source.document_id, snapshot.source.event_id, snapshot.source.revision
    );
    let status = RwSignal::new(String::new());
    view! {
        <section>
            <h3>{snapshot.molecule.name.clone()}</h3>
            <p class="helper">"Original translation source version · read-only. Opening this view does not move the source document's current head."</p>
            <details style="overflow-wrap:anywhere"><summary>"Source identity"</summary>
                <p>{format!("Document {} · event {} · revision {}",snapshot.source.document_id,snapshot.source.event_id,snapshot.source.revision)}</p>
            </details>
            <button class="button secondary" on:click=move|_|on_back.run(())>"Back to protein"</button>
            <p class="helper" role="status">{move||if current.get(){status.get()}else{"Workspace changed; reopen the source from the protein.".into()}}</p>
            <crate::assembly::review::assembly_product_view::ReferenceCanvas product=snapshot.molecule.as_ref().clone()
                identity current initial_spans=snapshot.spans on_message=Callback::new(move|text|status.set(text))/>
        </section>
    }.into_any()
}
