//! Conditional coding evidence alongside the same native pairwise result view.
use crate::viewer::alignment::{compare::Navigation, pairwise_result_view::Review};
use collagene_core::alignment::design_consensus as native;
use leptos::prelude::*;
use std::sync::Arc;

#[component]
pub(crate) fn DesignConsensusEvidenceView(
    review: Arc<Review>,
    guard: Callback<Arc<Review>, Result<(), String>>,
    on_navigate: Callback<Navigation, Result<(), String>>,
) -> impl IntoView {
    let bound = StoredValue::new(review);
    let status = RwSignal::new(String::new());
    let cds_page = RwSignal::new(0usize);
    let navigate = Callback::new(move |(row, consequence, consensus): (usize, usize, bool)| {
        let outcome = (|| -> Result<(), String> {
            let shown = bound.get_value();
            guard.run(shown.clone())?;
            let evidence = shown
                .record
                .design_consensus
                .as_ref()
                .ok_or("The retained comparison has no design/consensus evidence.")?;
            let consequence = evidence
                .cds
                .get(row)
                .and_then(|row| row.consequences.get(consequence))
                .ok_or("The selected codon consequence is unavailable.")?;
            let (event, coordinates) = if consensus {
                (
                    &evidence.request.consensus,
                    &consequence.consensus_coordinates,
                )
            } else {
                (&evidence.request.design, &consequence.design_coordinates)
            };
            if coordinates.is_empty() {
                return Err(
                    "No exact original-source coordinates are available for this side.".into(),
                );
            }
            // The native owner supplies original-source spans. Do not infer
            // correspondence from codon text, feature labels or in-range offsets.
            on_navigate.run(Navigation {
                event: event.clone(),
                epoch: shown.epoch,
                spans: coordinates
                    .iter()
                    .map(|span| (span.start, span.end))
                    .collect(),
            })
        })();
        if let Err(error) = outcome {
            status.set(error);
        }
    });
    view! {
        <section class="design-consensus-review">
            <h3>"Design / supplied-consensus interpretation"</h3>
            <p class="helper">
                "These coding consequences apply to one selected alignment path. Even a proven best score
                does not establish a unique correspondence, sample identity or biological function."
            </p>
            {bound.with_value(|shown| {
                shown.record.design_consensus.as_ref().map(|evidence| {
                    view! {
                        <p>{evidence.statement.clone()}</p>
                        <p>{format!(
                            "Coding summaries complete: {} · conditional on selected path: {}",
                            evidence.coding_complete, evidence.selected_path_conditional,
                        )}</p>
                        <details>
                            <summary>"Exact roles, assumptions and source provenance"</summary>
                            <p>{format!("Design: {} · {} · {}",
                                evidence.request.design.document_id,
                                evidence.request.design.event_id,
                                evidence.request.design.revision,
                            )}</p>
                            <p>{format!("Supplied consensus: {} · {} · {}",
                                evidence.request.consensus.document_id,
                                evidence.request.consensus.event_id,
                                evidence.request.consensus.revision,
                            )}</p>
                            <p>{match &evidence.request.consensus_axis {
                                native::ConsensusAxis::StoredTopology =>
                                    "Consensus axis: stored source topology".to_owned(),
                                native::ConsensusAxis::DeclaredSingleCircularTurn { rationale } =>
                                    format!("Declared single circular turn: {rationale}"),
                            }}</p>
                            {evidence.provenance.iter().map(|source| view! {
                                <article>
                                    <h4>{match &source.role {
                                        native::SourceRole::Design => "Design source",
                                        native::SourceRole::Consensus => "Supplied consensus source",
                                    }}</h4>
                                    <p>{source.original.message.clone()}</p>
                                    <p>{format!("Original bytes: {} · original evidence {:?}",
                                        source.original.byte_count, source.original.state,
                                    )}</p>
                                    {source.fields.iter().map(|field| view! {
                                        <p>{format!("{}: {}{}", field.key, field.value,
                                            if field.value_complete { "" } else { " (shortened or cleaned for display; exact source remains revision-bound)" },
                                        )}</p>
                                    }).collect_view()}
                                </article>
                            }).collect_view()}
                            <p>{format!("Context {} · evidence {} · workspace {}",
                                evidence.context_revision,
                                evidence.evidence_revision,
                                evidence.workspace_fingerprint,
                            )}</p>
                            <p>{format!("Mapping work {} · coding calls {} · segments {} · consequences {}",
                                evidence.work.mapping_work, evidence.work.coding_calls,
                                evidence.work.segments, evidence.work.consequences,
                            )}</p>
                        </details>
                    }
                })
            })}
            {move || {
                let shown = bound.get_value();
                let Some(evidence) = shown.record.design_consensus.as_ref() else {
                    return view! { <p>"Design/consensus evidence is unavailable."</p> }.into_any();
                };
                let pages = evidence.cds.len().div_ceil(4).max(1);
                let current = cds_page.get().min(pages - 1);
                let rows = evidence.cds.iter().enumerate().skip(current * 4).take(4)
                    .map(|(row_index, row)| view! {
                        <CodingRow row=row.clone() row_index navigate/>
                    }).collect_view();
                view! {
                    <section>
                        <h4>"Conditional design CDS summaries"</h4>
                        <p class="helper">"An empty list means no CDS was declared. Unsupported or incomplete projection is not evidence that a change or stop is absent."</p>
                        {rows}
                        <div class="method-row-actions">
                            <button class="text-tool" disabled=current == 0
                                on:click=move |_| cds_page.set(current.saturating_sub(1))>
                                "Previous CDS summaries"
                            </button>
                            <span>{format!("Page {} of {pages}", current + 1)}</span>
                            <button class="text-tool" disabled=current + 1 >= pages
                                on:click=move |_| cds_page.set((current + 1).min(pages - 1))>
                                "Next CDS summaries"
                            </button>
                        </div>
                    </section>
                }.into_any()
            }}
            <p class="helper" role="status">{move || status.get()}</p>

        </section>
    }
}

#[component]
fn CodingRow(
    row: native::CdsReport,
    row_index: usize,
    navigate: Callback<(usize, usize, bool)>,
) -> impl IntoView {
    let row = StoredValue::new(row);
    let page = RwSignal::new(0usize);
    let segment_page = RwSignal::new(0usize);
    view! {
        <article>
            <h4>{row.with_value(|row| format!("{} · genetic code {}", row.id, row.code.descriptor.id))}</h4>
            <details>
                <summary>"Exact genetic-code source"</summary>
                <p>{row.with_value(|row| format!("NCBI table {}: {} · pin {}",
                    row.code.descriptor.id, row.code.descriptor.name, row.code.revision,
                ))}</p>
                <p>{row.with_value(|row| format!("Start codons: {} · context-dependent stop codons: {}",
                    row.code.descriptor.start_codons.join(", "),
                    row.code.descriptor.context_dependent_stop_codons.join(", "),
                ))}</p>
                <p class="helper">"Ordinary amino-acid assignments in TCAG codon enumeration order"</p>
                <code>{row.with_value(|row| row.code.ordinary_amino_acids.clone())}</code>
            </details>
            <p>{row.with_value(|row| format!(
                "Mapping complete: {} · register complete: {} · consequences complete: {}",
                row.mapping_complete, row.register_complete, row.consequences_complete,
            ))}</p>
            <p>{row.with_value(|row| format!(
                "Design translation: {} · projected translation: {}",
                row.design_translation.as_ref().map(|translation| format!(
                    "{} codons; {} ambiguous, {} unresolved; {} trailing bases",
                    translation.translated_codons, translation.ambiguous_codons,
                    translation.unresolved_codons, translation.trailing_bases,
                )).unwrap_or_else(|| "unavailable".into()),
                row.projected_translation.as_ref().map(|translation| format!(
                    "{} codons; {} ambiguous, {} unresolved; {} trailing bases",
                    translation.translated_codons, translation.ambiguous_codons,
                    translation.unresolved_codons, translation.trailing_bases,
                )).unwrap_or_else(|| "unavailable".into()),
            ))}</p>
            {row.with_value(|row| row.findings.iter().map(|finding| view! {
                <p>{format!("{:?}: {}", finding.verdict, finding.message)}</p>
            }).collect_view())}
            <details>
                <summary>"Coding register segments"</summary>
                <p class="helper">"Intervals are coding ordinals after the declared codon-start prefix, not genomic source positions."</p>
                {move || row.with_value(|row| {
                    let pages = row.register_segments.len().div_ceil(20).max(1);
                    let current = segment_page.get().min(pages - 1);
                    let segments = row.register_segments.iter().skip(current * 20).take(20)
                        .map(|segment| view! {
                            <p>{format!("Design [{}, {}) → projected [{}, {}); offset modulo 3 = {}",
                                segment.design_start, segment.design_end,
                                segment.projected_start, segment.projected_end, segment.offset_mod3,
                            )}</p>
                        }).collect_view();
                    view! {
                        {segments}
                        <button class="text-tool" disabled=current == 0
                            on:click=move |_| segment_page.set(current.saturating_sub(1))>
                            "Previous register segments"
                        </button>
                        <span>{format!("Page {} of {pages}", current + 1)}</span>
                        <button class="text-tool" disabled=current + 1 >= pages
                            on:click=move |_| segment_page.set((current + 1).min(pages - 1))>
                            "Next register segments"
                        </button>
                    }
                })}
            </details>
            <details>
                <summary>"Codon consequences on this path"</summary>
                {move || row.with_value(|row| {
                    let pages = row.consequences.len().div_ceil(20).max(1);
                    let current = page.get().min(pages - 1);
                    let consequences = row.consequences.iter().enumerate()
                        .skip(current * 20).take(20).map(|(index, consequence)| {
                            let design_missing = consequence.design_coordinates.is_empty();
                            let consensus_missing = consequence.consensus_coordinates.is_empty();
                            view! {
                                <article>
                                    <p>{format!("Projected codon {}: {} ({}) · {:?}",
                                        consequence.projected_codon_index + 1,
                                        consequence.projected_codon,
                                        consequence.projected_amino_acid,
                                        consequence.kind,
                                    )}</p>
                                    <p>{consequence.message.clone()}</p>
                                    <button class="text-tool" disabled=design_missing
                                        on:click=move |_| navigate.run((row_index, index, false))>
                                        "Show original design bases"
                                    </button>
                                    <button class="text-tool" disabled=consensus_missing
                                        on:click=move |_| navigate.run((row_index, index, true))>
                                        "Show supplied-consensus bases"
                                    </button>
                                </article>
                            }
                        }).collect_view();
                    view! {
                        {consequences}
                        <button class="text-tool" disabled=current == 0
                            on:click=move |_| page.set(current.saturating_sub(1))>
                            "Previous consequences"
                        </button>
                        <span>{format!("Page {} of {pages}", current + 1)}</span>
                        <button class="text-tool" disabled=current + 1 >= pages
                            on:click=move |_| page.set((current + 1).min(pages - 1))>
                            "Next consequences"
                        </button>
                    }
                })}
            </details>
        </article>
    }
}
