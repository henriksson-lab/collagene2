use collagene_core::cloning::screening as native;
use leptos::prelude::*;
use std::sync::Arc;

pub(super) fn state(state: &native::CellState) -> &'static str {
    match state {
        native::CellState::Complete => "Complete modeled signature",
        native::CellState::Ambiguous => "Unresolved / ambiguous",
        native::CellState::Limited => "Limited search or retained output",
        native::CellState::Unsearched => "Not searched",
        native::CellState::Unsupported => "Outside the supported model",
    }
}
fn sizes(values: &[usize]) -> String {
    if values.is_empty() {
        return "No sizes retained (see completeness)".into();
    }
    values
        .iter()
        .map(|value| {
            if value % 2 == 0 {
                format!("{} bp", value / 2)
            } else {
                format!("{}.5 bp", value / 2)
            }
        })
        .collect::<Vec<_>>()
        .join(", ")
}
fn distinction(value: &native::Distinction) -> String {
    match value {
        native::Distinction::Separated(witness) => format!(
            "Distinguished by {}-only modeled size {}",
            match witness.side {
                native::WitnessSide::Intended => "intended",
                native::WitnessSide::Decoy => "alternative",
            },
            sizes(&[witness.size_twice_bp])
        ),
        native::Distinction::IdenticalSignature => "Identical ideal size signature".into(),
        native::Distinction::WithinDeclaredTolerance => {
            "No unmatched size under the declared tolerance".into()
        }
        native::Distinction::Unresolved => {
            "Not established: one or more component searches unresolved".into()
        }
    }
}
pub(super) fn results(
    report: Arc<native::Report>,
    selected: usize,
    navigate: Callback<(usize, usize)>,
) -> impl IntoView {
    let selected = selected.min(report.tests.len().saturating_sub(1));
    let panel_labels = report
        .panel
        .test_indices
        .iter()
        .filter_map(|index| report.tests.get(*index).map(|t| t.definition.label.clone()))
        .collect::<Vec<_>>()
        .join(" + ");
    let panel_summary = format!(
        "{} of {} named alternatives distinguished; {} combinations evaluated / {} comparison-work units. {}",
        report.panel.separated_decoy_indices.len(),
        report.subjects.len().saturating_sub(1),
        report.panel.evaluated_combinations,
        report.panel.comparison_work,
        if report.panel.complete_through_size_three {
            "All candidate panels of one to three tests were compared."
        } else {
            "Search limited: best among the panels actually compared; unsearched panels may do better."
        }
    );
    let selected_report = report.tests.get(selected);
    view! {
        <p class="helper">{report.statement.clone()}</p>
        <table class="analysis-table">
            <thead><tr><th>"Ranked candidate"</th><th>"Distinguished alternatives"</th><th>"Complete comparisons"</th></tr></thead>
            <tbody>{report.ranked_test_indices.iter().filter_map(|index| report.tests.get(*index)).map(|test| view! {
                <tr><td>{test.definition.label.clone()}</td><td>{test.separated_decoys}</td><td>{test.complete_comparisons}</td></tr>
            }).collect_view()}</tbody>
        </table>
        <details class="biology-section">
            <summary>"Modeled panel of up to three tests"</summary>
            <p>{if panel_labels.is_empty() { "No evaluated panel establishes a distinction.".into() } else { panel_labels }}</p>
            <p>{panel_summary}</p>
            <p class="helper">"Coverage is a deterministic size-model comparison; it is not a recommendation about yield or experimental success."</p>
        </details>
        {selected_report.map(|test| view! {
            <h4>{test.definition.label.clone()}</h4>
            <table class="analysis-table">
                <thead><tr><th>"Named alternative"</th><th>"Evidence"</th></tr></thead>
                <tbody>{test.comparisons.iter().map(|comparison| {
                    let label = report.subjects.get(comparison.decoy_subject_index).map(|s| s.label.clone()).unwrap_or_default();
                    view! { <tr><td>{label}</td><td>{distinction(&comparison.distinction)}</td></tr> }
                }).collect_view()}</tbody>
            </table>
            <details class="biology-section">
                <summary>"Named outcome signatures"</summary>
                {test.signatures.iter().map(|signature| {
                    let label = report.subjects.get(signature.subject_index).map(|s| s.label.clone()).unwrap_or_default();
                    let size_label = if signature.state == native::CellState::Complete && signature.sizes_twice_bp.is_empty() {
                        "Complete empty modeled signal".into()
                    } else {
                        sizes(&signature.sizes_twice_bp)
                    };
                    view! { <p><strong>{label}</strong>{format!(" · {} · {}", state(&signature.state), size_label)}</p> }
                }).collect_view()}
            </details>
        })}
        <details class="biology-section">
            <summary>"Component searches and stored strand coordinates"</summary>
            {report.cells.iter().enumerate().filter(|(_, cell)| cell.test_index == selected).map(|(cell_index, cell)| {
                let component = report.components.get(cell.component_index);
                let label = component.map(|c| c.label.clone()).unwrap_or_default();
                let original = component.is_some_and(|c| c.origin == native::ComponentOrigin::OriginalReference);
                let geometry = format!(
                    "{} modeled products · finite enumeration {} · {}",
                    cell.native_product_count,
                    if cell.enumeration_complete { "complete" } else { "not established" },
                    sizes(&cell.sizes_twice_bp),
                );
                view! {
                    <details>
                        <summary>{format!("{} · {}", label, state(&cell.state))}</summary>
                        <p>{geometry}</p><p>{cell.reason.clone()}</p>
                        <p class="helper">{format!("Exact source revision: {}", cell.source_revision)}</p>
                        {cell.cleavage.as_ref().map(|evidence| view! {
                            {evidence.sites.iter().enumerate().map(|(site_index, site)| {
                                let text = format!(
                                    "{} · recognition starts at base {} · top bond {:?} / bottom bond {:?} · {:?}",
                                    site.enzyme, site.recognition_start + 1,
                                    site.top_cut, site.bottom_cut, site.outcome,
                                );
                                view! {
                                    <p>{text}
                                        {original.then(|| view! { <button on:click=move |_| navigate.run((cell_index, site_index))>"Select recognition site"</button> })}
                                    </p>
                                }
                            }).collect_view()}
                            <p class="helper">"Bond coordinates are zero-based stored strand bonds. Saved material may have a different frame; only original reference sites offer map selection."</p>
                            {evidence.warnings.iter().map(|warning| view! { <p class="helper">{warning.clone()}</p> }).collect_view()}
                        })}
                        {cell.pcr_geometries.iter().map(|geometry| view! {
                            <p>{format!(
                                "Primer starts {} / {} (zero-based selected template strand); {} bp product{}",
                                geometry.forward_start, geometry.reverse_start, geometry.product_length,
                                if geometry.crosses_origin { "; crosses circular origin" } else { "" },
                            )}</p>
                        }).collect_view()}
                    </details>
                }
            }).collect_view()}
        </details>
    }
}
