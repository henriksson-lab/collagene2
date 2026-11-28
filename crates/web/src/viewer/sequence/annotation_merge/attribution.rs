//! Inspection of the native per-feature attribution ledger; no qualifier parser.
use collagene_core::sequence::annotations::merge::{Contribution, Ledger, SourcePin};
use collagene_model::Molecule;
use leptos::prelude::*;

pub(super) fn contribution_text(contribution: &Contribution) -> String {
    let source = match &contribution.source {
        SourcePin::Library {
            library_revision,
            library_name,
            library_version,
            entry_id,
            entry_revision,
            source_revision,
            source_feature_index,
            proposal_id,
            concrete_matches,
            compatible_ambiguous,
            mismatches,
            template_length,
        } => format!(
            "Library {library_name} v{library_version} ({library_revision});
            entry {entry_id} ({entry_revision}); source {source_revision}, annotation {};
            proposal {proposal_id}. Template {template_length} bases: {concrete_matches} concrete,
            {compatible_ambiguous} ambiguity-compatible, {mismatches} mismatches.",
            source_feature_index + 1,
        ),
        SourcePin::Dataset {
            source_blake3,
            format,
            seqid,
            proposal_id,
            source_rows,
        } => format!(
            "Imported {format:?}; source {source_blake3}; seqid {seqid}; proposal {proposal_id}; {} source rows (first ten: {:?}).",
            source_rows.len(),
            source_rows.iter().take(10).collect::<Vec<_>>(),
        ),
        SourcePin::Donor {
            event,
            feature_index,
            annotation_revision,
            transform,
        } => format!(
            "Retained donor {} / {} / {}; annotation {} ({annotation_revision}); reverse {}, origin {} after optional reverse.",
            event.document_id,
            event.event_id,
            event.revision,
            feature_index + 1,
            transform.reverse,
            transform.origin,
        ),
    };
    format!(
        "{source} Contribution feature revision: {}.",
        contribution.feature_revision
    )
}

#[component]
pub(crate) fn AttributionInspection(
    molecule: Memo<Molecule>,
    selected_feature: RwSignal<Option<usize>>,
) -> impl IntoView {
    let ledger = Memo::new(move |_| {
        molecule.with(|molecule| {
            selected_feature.get().map(|index| {
                let feature = molecule
                    .features
                    .get(index)
                    .ok_or("Select a current annotation.")?;
                Ledger::from_feature(feature).map_err(|error| error.to_string())
            })
        })
    });
    view! {
        <details>
            <summary>"Source attribution for selected annotation"</summary>
            <p class="helper">"Attribution records documentary contributions. It does not authenticate a source or verify biological function."</p>
            {move || match ledger.get() {
                None => view! { <p>"Select an annotation in the feature list or map."</p> }.into_any(),
                Some(Err(error)) => view! { <p>{error}</p> }.into_any(),
                Some(Ok(ledger)) => {
                    let count = ledger.contributions.len();
                    let rows = ledger.contributions.iter().map(|contribution| view! {
                        <p>{contribution_text(contribution)}</p>
                    }).collect_view();
                    view! { <p>{format!("{count} retained contributions")}</p>{rows} }.into_any()
                }
            }}
        </details>
    }
}
