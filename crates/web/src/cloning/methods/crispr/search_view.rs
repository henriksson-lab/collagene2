//! Paged selected-reference evidence with explicit per-path sequence columns.
use super::{BoundSearch, boundary, positions};
use collagene_core::sequence::crispr::{
    OffTargetHit,
    alignment::{AlignmentColumn, GuideAlignment},
};
use leptos::prelude::*;
use std::sync::Arc;

#[component]
pub(super) fn WholeField(label: &'static str, value: RwSignal<String>) -> impl IntoView {
    view! {
        <label>{label}<input type="number" min="0" step="1" aria-label=label
            prop:value=move || value.get()
            on:input=move |event| value.set(event_target_value(&event))/></label>
    }
}

fn base(sequence: &str, index: usize) -> char {
    sequence
        .as_bytes()
        .get(index)
        .copied()
        .map(char::from)
        .unwrap_or('?')
}

fn columns(query: &str, hit: &OffTargetHit, alignment: &GuideAlignment) -> String {
    alignment
        .columns
        .iter()
        .map(|column| match *column {
            AlignmentColumn::Paired { guide, target } => format!(
                "guide {} {} ↔ oriented DNA {} {} · source {}",
                guide + 1,
                base(query, guide),
                target + 1,
                base(&hit.reference_guide, target),
                hit.guide_positions
                    .get(target)
                    .map(|position| (position + 1).to_string())
                    .unwrap_or_else(|| "unavailable".into()),
            ),
            AlignmentColumn::DnaBulge { target } => format!(
                "guide gap ↔ oriented DNA {} {} · source {} · extra target DNA",
                target + 1,
                base(&hit.reference_guide, target),
                hit.guide_positions
                    .get(target)
                    .map(|position| (position + 1).to_string())
                    .unwrap_or_else(|| "unavailable".into()),
            ),
            AlignmentColumn::RnaBulge { guide } => format!(
                "guide {} {} ↔ DNA gap · extra guide RNA; no source base assigned",
                guide + 1,
                base(query, guide),
            ),
        })
        .collect::<Vec<_>>()
        .join("\n")
}

#[component]
pub(super) fn SearchResults(
    search: RwSignal<Option<Arc<BoundSearch>>>,
    page: RwSignal<usize>,
    selected: RwSignal<Option<String>>,
    navigate: Callback<String>,
    export: Callback<()>,
) -> impl IntoView {
    view! {
        {move || search.get().map(|bound| view! {
            <SearchPage bound page selected navigate export/>
        })}
    }
}

#[component]
fn SearchPage(
    bound: Arc<BoundSearch>,
    page: RwSignal<usize>,
    selected: RwSignal<Option<String>>,
    navigate: Callback<String>,
    export: Callback<()>,
) -> impl IntoView {
    let report = &bound.report;
    let summary = format!(
        "{} selected records · {} retained loci · {} loci include the declared exact match · {} definite · {} charged work units",
        report.references.len(),
        report.hits.len(),
        report.hits.iter().filter(|h| h.is_intended).count(),
        report.hits.iter().filter(|h| h.definite).count(),
        report.comparisons,
    );
    let completion = if report.exhaustive_for_selected_records {
        "Enumeration is complete for these records under the exact declared alignment and PAM policy. No absence claim applies to unselected records or an unrepresented genome."
    } else {
        "Enumeration is incomplete. Retained loci and paths cannot establish absence of additional matches."
    };
    let policy = format!(
        "Alignment model: {}. Screening PAM: {:?}. Design PAM: {}.",
        report.alignment_model, report.options.pam_policy, report.options.profile.pam
    );
    let acceleration = format!(
        "Index {} · {} entries. Bulged searches use bounded variable-span scanning.",
        if report.index_used {
            "used"
        } else {
            "not used"
        },
        report.index_entries
    );
    let warnings = report
        .limits
        .iter()
        .chain(&report.warnings)
        .map(|warning| {
            view! {<p class="helper">{warning.clone()}</p>}
        })
        .collect_view();
    let count = report.hits.len();
    let bound = StoredValue::new(bound);
    view! {
        <div class="crispr-search-result">
            <p>{summary}</p><p class="helper">{completion}</p>
            <details><summary>"Declared search model"</summary><p>{policy}</p><p>{acceleration}</p>
                <p class="helper">"Columns describe full-input sequence alignments, including terminal and mixed gap kinds. They are not an assertion about physical heteroduplex structure, cellular activity or cleavage."</p>
            </details>
            <button class="button secondary" on:click=move |_| export.run(())>"Download off-target report"</button>
            {move || bound.with_value(|bound| bound.report.hits.iter().skip(page.get() * 25).take(25).map(|hit| {
                let id = hit.id.clone();
                let attr_id = id.clone();
                let name = bound.report.references.get(hit.reference_index).map(|r| r.id.as_str()).unwrap_or("missing reference");
                let text = format!("{} · {} · mismatches {}–{} · {} · PAM {} · {} paths{}", name,
                    if hit.is_intended { "includes exact intended match" } else { "additional locus" },
                    hit.min_mismatches, hit.max_mismatches, if hit.reverse { "reverse" } else { "forward" },
                    hit.pam_sequence, hit.alignment.alignments.len(),
                    if hit.alignment.complete { "" } else { " · incomplete path enumeration" });
                view! {
                    <button class="button secondary crispr-hit-entry" data-hit-id=attr_id
                        on:click=move |_| selected.set(Some(id.clone()))>{text}</button>
                }
            }).collect_view())}
            <div class="crispr-actions">
                <button class="button secondary" disabled=move || page.get() == 0
                    on:click=move |_| page.update(|p| *p = p.saturating_sub(1))>"Previous loci"</button>
                <button class="button secondary" disabled=move || (page.get() + 1) * 25 >= count
                    on:click=move |_| page.update(|p| *p += 1)>"Next loci"</button>
            </div>
            {move || selected.get().and_then(|id| bound.with_value(|bound| {
                bound.report.hits.iter().find(|hit| hit.id == id).map(|hit| (bound.clone(), hit.clone()))
            })).map(|(bound, hit)| view! {<HitDetail bound hit navigate/>})}
            {warnings}
        </div>
    }
}

#[component]
fn HitDetail(
    bound: Arc<BoundSearch>,
    hit: OffTargetHit,
    navigate: Callback<String>,
) -> impl IntoView {
    let name = bound
        .report
        .references
        .get(hit.reference_index)
        .map(|r| format!("{} · {}", r.id, r.name))
        .unwrap_or_default();
    let summary = format!(
        "Oriented reference spacer {} · PAM {} · mismatches {}–{} · design-PAM evidence {:?}",
        hit.reference_guide,
        hit.pam_sequence,
        hit.min_mismatches,
        hit.max_mismatches,
        hit.design_pam
    );
    let source_positions = format!(
        "Target DNA bases {} · PAM bases {} (1-based, ordered)",
        positions(&hit.guide_positions),
        positions(&hit.pam_positions)
    );
    let gapped = hit
        .alignment
        .alignments
        .iter()
        .any(|path| path.dna_bulge_bases > 0 || path.rna_bulge_bases > 0);
    let cut = if gapped {
        "Nominal cut boundaries are unassigned: no cleavage-shift model is admitted for gapped paths.".to_string()
    } else {
        format!(
            "Nominal top boundary: {}. Nominal bottom boundary: {}. These are declared geometry, not cleavage predictions.",
            boundary(hit.nominal_top_cut),
            boundary(hit.nominal_bottom_cut)
        )
    };
    let complete = hit.alignment.complete;
    let limit = format!(
        "Path enumeration complete: {}; limit: {:?}; charged alignment work: {}.",
        complete, hit.alignment.limit_reached, hit.alignment.work
    );
    let paths = hit
        .alignment
        .alignments
        .iter()
        .enumerate()
        .map(|(index, alignment)| {
            let heading = format!(
                "Path {} · mismatches {}–{} · extra DNA {} / RNA {} · {}",
                index + 1,
                alignment.min_mismatches,
                alignment.max_mismatches,
                alignment.dna_bulge_bases,
                alignment.rna_bulge_bases,
                if alignment.definite {
                    "certain within declared bounds"
                } else {
                    "possible under ambiguity"
                }
            );
            let text = columns(&bound.report.query.sequence, &hit, alignment);
            view! {<details><summary>{heading}</summary><pre>{text}</pre></details>}
        })
        .collect_view();
    let id = hit.id;
    view! {
        <article class="crispr-hit-detail">
            <h4>{name}</h4><p>{summary}</p><p>{source_positions}</p><p>{cut}</p><p>{limit}</p>
            <p class="helper">"An intended-locus tag requires the exact unbulged design match. Other retained alignments at that locus remain explicit; the tag does not remove them."</p>
            <button class="button secondary" on:click=move |_| navigate.run(id.clone())>"Open reference spacer position"</button>
            {paths}
        </article>
    }
}
