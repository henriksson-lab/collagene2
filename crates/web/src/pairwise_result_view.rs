//! Complete paged pairwise evidence and native-coordinate discrepancy navigation.
use crate::{
    app::download,
    compare::Navigation,
    pairwise_pages::{PageIndex, render_page},
};
use collagene_core::pairwise_history::Record;
use leptos::prelude::*;
use std::{io, sync::Arc};

pub struct Review {
    pub record: Record,
    pub reference_name: String,
    pub query_name: String,
    pub fresh: bool,
    pub epoch: u64,
}
fn clipped(value: &str) -> String {
    if value.is_empty() {
        "∅".into()
    } else if value.len() <= 100 {
        value.into()
    } else {
        format!(
            "{}… ({} bases; complete values retained in JSON)",
            &value[..100],
            value.len()
        )
    }
}
struct Limited(Vec<u8>);
impl io::Write for Limited {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        if self.0.len().saturating_add(bytes.len())
            > collagene_core::pairwise_history::MAX_SAVED_BYTES
        {
            return Err(io::Error::other("comparison JSON limit"));
        }
        self.0.extend_from_slice(bytes);
        Ok(bytes.len())
    }
    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}

#[component]
pub fn PairwiseResultView(
    review: Arc<Review>,
    guard: Callback<Arc<Review>, Result<(), String>>,
    on_save: Callback<Record, Result<(), String>>,
    on_navigate: Callback<Navigation, Result<(), String>>,
) -> impl IntoView {
    let index = StoredValue::new(PageIndex::new(&review.record.report.comparison.alignment));
    let count = review.record.report.comparison.alignment.differences.len();
    let columns = review
        .record
        .report
        .comparison
        .alignment
        .reference_aligned
        .len();
    let alignment_pages = index.with_value(PageIndex::page_count);
    let difference_pages = count.div_ceil(20).max(1);
    let cigar_pages = review
        .record
        .report
        .comparison
        .alignment
        .cigar
        .len()
        .div_ceil(2000)
        .max(1);
    let bound = StoredValue::new(review);
    let page = RwSignal::new(0usize);
    let differences_page = RwSignal::new(0usize);
    let cigar_page = RwSignal::new(0usize);
    let focused = RwSignal::new((count > 0).then_some(0usize));
    let jump = RwSignal::new(String::new());
    let column_jump = RwSignal::new(String::new());
    let status = RwSignal::new(String::new());
    let focus = Callback::new(move |difference: usize| {
        if difference >= count {
            return;
        }
        focused.set(Some(difference));
        differences_page.set(difference / 20);
        jump.set((difference + 1).to_string());
        if let Some(column) = index.with_value(|index| index.difference_column(difference)) {
            page.set(column / 600);
        }
    });
    let navigate = Callback::new(move |(difference, query): (usize, bool)| {
        let review = bound.get_value();
        let result = (|| -> Result<(), String> {
            guard.run(review.clone())?;
            let d = review
                .record
                .report
                .comparison
                .alignment
                .differences
                .get(difference)
                .ok_or("Difference is unavailable")?;
            let (event, spans) = if query {
                (
                    review.record.query.clone(),
                    vec![
                        review
                            .record
                            .report
                            .map_query(d.query_start, d.query_end)
                            .map_err(|e| e.to_string())?,
                    ],
                )
            } else {
                (
                    review.record.reference.clone(),
                    review
                        .record
                        .report
                        .map_reference(d.reference_start, d.reference_end)
                        .map_err(|e| e.to_string())?,
                )
            };
            on_navigate.run(Navigation {
                event,
                epoch: review.epoch,
                spans,
            })
        })();
        if let Err(error) = result {
            status.set(error);
        }
    });
    let save = move |_| {
        let review = bound.get_value();
        let result = (|| -> Result<(), String> {
            guard.run(review.clone())?;
            if !review.fresh {
                return Err(
                    "Recompute this historical comparison before saving a current result.".into(),
                );
            }
            on_save.run(review.record.clone())
        })();
        status.set(match result{Ok(())=>"Comparison retained with exact input versions. Reopening remains historical until explicit recomputation.".into(),Err(error)=>error});
    };
    let export = move |_| {
        let review = bound.get_value();
        let result = (|| -> Result<(), String> {
            guard.run(review.clone())?;
            let mut out = Limited(Vec::new());
            serde_json::to_writer_pretty(&mut out, &review.record).map_err(|e| e.to_string())?;
            guard.run(review)?;
            download("pairwise-comparison.json", "application/json", &out.0)
        })();
        status.set(match result{Ok(())=>"Complete comparison evidence downloaded; recorded data is not experimental verification.".into(),Err(error)=>error});
    };
    view!{
        <section class="compare-result">
            <h3>{bound.with_value(|review|format!("{} → {}",review.reference_name,review.query_name))}</h3>
            <p class="helper">{bound.with_value(|review|if review.fresh{"Fresh native computation on the pinned input versions"}else{"Historical saved comparison · not recomputed on opening"})}</p>
            <p class="helper">"Differences describe sequence models, not sample identity, read confidence or experimental success. Navigation uses original document coordinates; alignment text uses the declared comparison orientation/origin."</p>
            {bound.with_value(|review|{
                let report=&review.record.report;let comparison=&report.comparison;let a=&comparison.alignment;let stats=&a.stats;
                view!{
                    <p>{format!("Score {} · reference [{}, {}) · query [{}, {})",a.score,a.reference_start,a.reference_end,a.query_start,a.query_end)}</p>
                    <p class="helper">{match report.options{collagene_core::pairwise_review::Options::Stored(_)=>"Stored origin/orientation: no reference-origin or query-strand search requested.",collagene_core::pairwise_review::Options::Whole(_)=>"Global whole-construct comparison under the explicitly requested origin/strand policy and compute budgets."}}</p>
                    <p>{format!("Reference origin {} · query reverse complement {} · evaluated candidates {}",comparison.reference_origin,comparison.query_reversed,comparison.evaluated_candidates)}</p>
                    <p>{format!("Requested reference origins exhaustive: {} · query orientations exhaustive: {} · best score proven: {}",comparison.exhaustive_reference_origins,comparison.exhaustive_query_orientations,comparison.optimality_proven)}</p>
                    <ul><li>{format!("Paired coverage: reference {:.1}% · query {:.1}%",stats.reference_coverage*100.,stats.query_coverage*100.)}</li><li>{format!("Certain matches {} · compatible ambiguity {} · mismatches {}",stats.certain_matches,stats.ambiguous_compatible,stats.mismatches)}</li><li>{format!("Inserted query bases {} · deleted reference bases {}",stats.inserted_query_bases,stats.deleted_reference_bases)}</li><li>{format!("Unaligned comparison flanks: reference {} left / {} right; query {} left / {} right",a.reference_start,report.reference_length.saturating_sub(a.reference_end),a.query_start,report.query_length.saturating_sub(a.query_end))}</li></ul>
                    {comparison.warnings.iter().chain(&a.warnings).cloned().map(|warning|view!{<p class="compare-warning">{warning}</p>}).collect_view()}
                    <details><summary>"Exact source and calculation details"</summary><p>{format!("Reference {} · {} · {}",review.record.reference.document_id,review.record.reference.event_id,review.record.reference.revision)}</p><p>{format!("Query {} · {} · {}",review.record.query.document_id,review.record.query.event_id,review.record.query.revision)}</p><p>{format!("{} · {} · context {}",report.engine_version,a.algorithm,report.context_revision)}</p></details>
                }
            })}
            <details open><summary>"Aligned sequence · all columns"</summary><p class="helper">"R/Q coordinates are zero-based boundaries in the rotated/oriented alignment. | certain; ? compatible ambiguity; × mismatch; gaps remain explicit."</p>
                {move||{let current=page.get().min(alignment_pages-1);let text=bound.with_value(|review|index.with_value(|index|render_page(&review.record.report.comparison.alignment,index,current)));view!{<pre>{text}</pre><p>{format!("Columns {}–{} of {columns}; page {} of {alignment_pages}",if columns==0{0}else{current*600+1},((current+1)*600).min(columns),current+1)}</p><div class="method-row-actions"><button class="text-tool" disabled=current==0 on:click=move|_|page.set(current.saturating_sub(1))>"Previous columns"</button><button class="text-tool" disabled={current+1>=alignment_pages} on:click=move|_|page.set((current+1).min(alignment_pages-1))>"Next columns"</button></div>}}}
                <label>"Jump to alignment column (1-based)"<input type="number" min="1" prop:value=move||column_jump.get() on:input=move|event|column_jump.set(event_target_value(&event))/></label><button class="text-tool" on:click=move|_|{match column_jump.get_untracked().parse::<usize>(){Ok(column)if column>0&&column<=columns=>page.set((column-1)/600),_=>status.set("Enter a reported alignment column.".into())}}>"Show column"</button>
            </details>
            <h4>{format!("Differences and uncertain regions · {count}")}</h4>
            <div class="method-row-actions"><button class="text-tool" disabled=move||focused.get().is_none_or(|index|index==0) on:click=move|_|if let Some(index)=focused.get_untracked(){focus.run(index.saturating_sub(1));}>"Previous difference"</button><button class="text-tool" disabled=move||focused.get().is_none_or(|index|index+1>=count) on:click=move|_|if let Some(index)=focused.get_untracked(){focus.run((index+1).min(count.saturating_sub(1)));}>"Next difference"</button></div>
            <label>"Difference number (1-based)"<input type="number" min="1" prop:value=move||jump.get() on:input=move|event|jump.set(event_target_value(&event))/></label><button class="text-tool" on:click=move|_|{match jump.get_untracked().parse::<usize>(){Ok(index)if index>0&&index<=count=>focus.run(index-1),_=>status.set("Enter a reported difference number.".into())}}>"Show difference"</button>
            {move||{
                let current=differences_page.get().min(difference_pages-1);
                let rows=bound.with_value(|review|review.record.report.comparison.alignment.differences.iter().enumerate().skip(current*20).take(20).map(|(index,d)|{
                    let before=clipped(&d.reference_bases);let after=clipped(&d.query_bases);let label=format!("{} · {:?}: R [{}, {}), Q [{}, {})",index+1,d.kind,d.reference_start,d.reference_end,d.query_start,d.query_end);
                    view!{<article class="compare-difference"><button class="text-tool" aria-pressed=move||(focused.get()==Some(index)).to_string() on:click=move|_|focus.run(index)>{label}</button><code>{format!("R: {before}\nQ: {after}")}</code><div class="method-row-actions"><button class="text-tool" on:click=move|_|navigate.run((index,false))>"Show reference bases / caret"</button><button class="text-tool" on:click=move|_|navigate.run((index,true))>"Show query bases / caret"</button></div></article>}
                }).collect_view());
                view!{<div>{rows}</div><p>{format!("Difference page {} of {difference_pages}",current+1)}</p><div class="method-row-actions"><button class="text-tool" disabled=current==0 on:click=move|_|differences_page.set(current.saturating_sub(1))>"Previous differences"</button><button class="text-tool" disabled={current+1>=difference_pages} on:click=move|_|differences_page.set((current+1).min(difference_pages-1))>"Next differences"</button></div>}
            }}
            <details><summary>"Complete CIGAR"</summary>{move||{let current=cigar_page.get().min(cigar_pages-1);let value=bound.with_value(|review|{let cigar=&review.record.report.comparison.alignment.cigar;let start=(current*2000).min(cigar.len());cigar[start..(start+2000).min(cigar.len())].to_string()});view!{<pre>{value}</pre><p>{format!("CIGAR text page {} of {cigar_pages}",current+1)}</p><button class="text-tool" disabled=current==0 on:click=move|_|cigar_page.set(current.saturating_sub(1))>"Previous CIGAR"</button><button class="text-tool" disabled={current+1>=cigar_pages} on:click=move|_|cigar_page.set((current+1).min(cigar_pages-1))>"Next CIGAR"</button>}}}</details>
            <button class="button secondary full" disabled=bound.with_value(|review|!review.fresh) on:click=save>"Save event-bound comparison"</button>
            <button class="button secondary full" on:click=export>"Download complete comparison JSON"</button>
            <p class="helper" role="status">{move||status.get()}</p>
        </section>
    }.into_any()
}
