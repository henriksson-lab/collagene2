//! Presentation of fresh preparation diagnostics and physical ends.
use leptos::prelude::*;
pub(super) fn guidance(minimum: Option<usize>, maximum: Option<usize>) -> String {
    match (minimum, maximum) {
        (Some(a), Some(b)) => format!("{a}–{b} bp"),
        (Some(a), None) => format!("at least {a} bp; no upper bound specified"),
        _ => "not established".into(),
    }
}
pub(super) fn short(sequence: &str) -> String {
    let s = sequence.chars().take(300).collect::<String>();
    if s.len() < sequence.len() {
        format!("{s}…")
    } else {
        s
    }
}
pub(super) fn method_status(
    status: &collagene_core::verification::evidence::Status,
) -> &'static str {
    match status {
        collagene_core::verification::evidence::Status::ConsistentUnderAssumptions => {
            "consistent under stated assumptions"
        }
        collagene_core::verification::evidence::Status::Incomplete => "incomplete",
        collagene_core::verification::evidence::Status::Invalid => "inconsistent",
    }
}
pub(crate) fn method_end(end: Option<&collagene_core::chemistry::End>) -> String {
    use collagene_core::chemistry::{FivePrime, Polarity, ThreePrime};
    let Some(end) = end else {
        return "End unavailable; no terminal chemistry inferred".into();
    };
    let shape = match end.polarity {
        Polarity::Blunt => "blunt".into(),
        Polarity::FivePrime => format!("5′ overhang {} (5′→3′)", end.overhang),
        Polarity::ThreePrime => format!("3′ overhang {} (5′→3′)", end.overhang),
    };
    format!(
        "{shape}; 5′ {}; 3′ {}",
        match end.chemistry.five_prime {
            FivePrime::Phosphate => "phosphate",
            FivePrime::Hydroxyl => "hydroxyl",
            FivePrime::Unknown => "chemistry unknown",
        },
        match end.chemistry.three_prime {
            ThreePrime::Hydroxyl => "hydroxyl",
            ThreePrime::Blocked => "blocked",
            ThreePrime::Unknown => "chemistry unknown",
        }
    )
}
pub(super) fn method_duplex(summary: &collagene_core::chemistry::DuplexSummary) -> String {
    format!(
        "Top {} bases; bottom {} bases; bottom reference offset {}. {}. Retained nicks: {}. Covalently closed in model: {}.",
        summary.top_length,
        summary.bottom_length,
        summary.bottom_start,
        match summary.ends.as_ref() {
            Some((left, right)) => format!(
                "Left: {}. Right: {}",
                method_end(Some(left)),
                method_end(Some(right))
            ),
            None => "No free ends on the reconstructed circular duplex".into(),
        },
        summary.nicks.len(),
        summary.covalently_closed_in_model
    )
}

#[component]
pub(super) fn MethodDiagnostics(
    diagnostics: collagene_core::cloning::assembly::diagnostics::Report,
    names: Vec<String>,
    preparations: Vec<String>,
    on_edit: Callback<usize>,
    #[prop(optional_no_strip)] on_source: Option<Callback<usize>>,
) -> impl IntoView {
    use collagene_core::cloning::assembly::diagnostics::{JunctionState, RowState, Stage};
    let status = method_status(&diagnostics.status);
    let stage_name = |stage: Stage| match stage {
        Stage::Pcr => "PCR",
        Stage::OverlapExtension => "reciprocal overlap extension",
        Stage::Digest => "restriction digest",
        Stage::Repair => "end repair",
        Stage::EndTreatment => "terminal treatment",
        Stage::SourceDeclaration => "source assumption",
        Stage::Annealing => "oligo annealing",
        Stage::MaterialReconstruction => "saved material reconstruction",
        Stage::OverlapAssembly => "overlap assembly",
        Stage::LigationAssembly => "ligation",
        Stage::TargetComparison => "independent-target comparison",
    };
    view! {<section class="method-diagnostics" aria-label="Preparation and junction diagnostics">

        <h3>"Preparation and junction review"</h3><p>{status}</p>
        <p class="helper">"These locations come from the fresh typed step graph and retained product prefix. Later blocked rows have not been reconstructed; they are not additional failed experiments. Editing a field clears this evidence; viewing its editor does not change the draft."</p>
        {diagnostics.failure.map(|failure|view!{<article class="diagnostic-failure"><h4>{format!("First blocking result: step {} · {}",failure.step+1,stage_name(failure.stage))}</h4><p>{failure.message}</p><p class="helper">{failure.recovery}</p>
            <div class="method-row-actions">{failure.row.map(|row|view!{<button class="button secondary" on:click=move|_|on_edit.run(row)>"Edit preparation"</button>})}{failure.upstream.map(|row|view!{<button class="text-tool" on:click=move|_|on_edit.run(row)>"Edit upstream preparation"</button>})}{failure.downstream.map(|row|view!{<button class="text-tool" on:click=move|_|on_edit.run(row)>"Edit downstream preparation"</button>})}</div>
        </article>})}
        <div class="diagnostic-table"><table><thead><tr><th>"Original fragment / direction"</th><th>"Declared preparation"</th><th>"Reconstruction"</th><th>"Review"</th></tr></thead><tbody>
            {diagnostics.rows.into_iter().map(|row|{let index=row.index;let state=match row.state{RowState::InputReference=>"Original reference; no preparation step",RowState::Reconstructed=>"Reconstructed",RowState::Failed=>"Failed at its preparation",RowState::Blocked=>"Blocked by an earlier result"};view!{<tr><td>{format!("{} · {} · {}",index+1,names.get(index).map(String::as_str).unwrap_or(&row.input_id),if row.reverse{"reverse after preparation"}else{"forward"})}</td><td>{preparations.get(index).cloned().unwrap_or_else(||"Preparation unavailable".into())}</td><td>{state}<p>{format!("{}/{} declared steps reconstructed",row.completed_steps,row.total_steps)}</p>{row.preparation_step.map(|step|view!{<p>{format!("Preparation step {}",step+1)}</p>})}{row.repair_step.map(|step|view!{<p>{format!("Repair step {}",step+1)}</p>})}</td><td><div class="method-row-actions"><button class="text-tool" aria-label=format!("Edit fragment {} preparation",index+1) on:click=move|_|on_edit.run(index)>"Edit preparation"</button>{on_source.map(|source|view!{<button class="text-tool" on:click=move|_|source.run(index)>"View original source"</button>})}</div></td></tr>}}).collect_view()}
        </tbody></table></div>
        {(!diagnostics.junctions.is_empty()).then(||view!{<h4>"Fresh junction evidence"</h4>})}
        {diagnostics.junctions.into_iter().map(|junction|{let upstream=junction.upstream;let downstream=junction.downstream;let state=match junction.state{JunctionState::ExactHomology=>"Exact homology",JunctionState::HomologyMismatch=>"Incompatible homology",JunctionState::HomologyUnavailable=>"Homology unavailable",JunctionState::EndpointsAvailable=>"Ends available; compatibility is not inferred here",JunctionState::EndpointsUnavailable=>"Physical ends unavailable"};view!{<details><summary>{format!("Junction {} · fragment {} → {}{} · {}",junction.index+1,upstream+1,downstream+1,if junction.closing{" (circle closure)"}else{""},state)}</summary><p>{junction.message}</p><p class="helper">{junction.recovery}</p><div class="method-row-actions"><button class="text-tool" on:click=move|_|on_edit.run(upstream)>"Edit upstream preparation"</button><button class="text-tool" on:click=move|_|on_edit.run(downstream)>"Edit downstream preparation"</button>{on_source.map(|source|view!{<button class="text-tool" on:click=move|_|source.run(upstream)>"View upstream source"</button><button class="text-tool" on:click=move|_|source.run(downstream)>"View downstream source"</button>})}</div></details>}}).collect_view()}
        <ul>{diagnostics.warnings.into_iter().map(|warning|view!{<li class="helper">{warning}</li>}).collect_view()}</ul>
        <details><summary>"Full fresh findings"</summary><ul>{diagnostics.findings.into_iter().map(|finding|view!{<li>{finding}</li>}).collect_view()}</ul></details>
    </section>}
}
