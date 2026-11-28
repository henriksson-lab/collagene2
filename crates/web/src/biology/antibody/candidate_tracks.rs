//! Simultaneous supported-interval overview; overlapping candidates are alternatives.
use collagene_core::antibody_numbering::domain_tracks::Report;
use leptos::prelude::*;
use std::sync::Arc;

#[component]
pub(super) fn CandidateTracks(
    report: Arc<Report>,
    selected: RwSignal<usize>,
    source_length: usize,
) -> impl IntoView {
    let source_length = source_length.max(1);
    view! {
        <details>
            <summary>"Retained candidate intervals"</summary>
            <p class="helper">
                "All retained supported intervals are shown together. Overlapping candidates can be alternative hypotheses, not separate domains. Bars do not establish unresolved N-terminal boundaries or biological classification. Select one for detailed labels."
            </p>
            {report.tracks.iter().enumerate().map(|(index, track)| {
                let start = track.first_supported_source_index;
                let end = track.source_end;
                let left = start as f64 * 100. / source_length as f64;
                let width = end.saturating_sub(start) as f64 * 100. / source_length as f64;
                view! {
                    <button
                        class="text-tool"
                        style="display:block;width:100%;text-align:left;margin-block:6px"
                        aria-pressed=move || (selected.get() == index).to_string()
                        on:click=move |_| selected.set(index)
                    >
                        {format!("{} · {:?} · supported {}–{}", index + 1, track.locus, start + 1, end)}
                        <span style="display:block;position:relative;height:6px;background:#e4ebe8;margin-top:4px">
                            <span style=format!("position:absolute;left:{left}%;width:{width}%;height:6px;background:#6b8c80")></span>
                        </span>
                    </button>
                }
            }).collect_view()}
        </details>
    }
}
