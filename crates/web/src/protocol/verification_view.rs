//! Typed verification presentation; source/step identities remain explicit.
use collagene_core::verification::checks::{CheckFinding, CheckScope};
use leptos::prelude::*;

pub(crate) fn finding(value: &CheckFinding) -> impl IntoView + use<> {
    let scope = match &value.scope {
        CheckScope::Plan => "Plan".into(),
        CheckScope::Input { input_id } => format!("Input {input_id}"),
        CheckScope::Step { step_index } => format!("Step {}", step_index + 1),
        CheckScope::Target { step_index } => {
            format!("Target comparison at step {}", step_index + 1)
        }
    };
    view! {
        <details class="verification-finding">
            <summary>{format!("{:?} · {} · {}", value.verdict, scope, value.rule.id())}</summary>
            <p>{value.message.clone()}</p>
            {(!value.subjects.is_empty()).then(|| view! {
                <p class="helper">{format!("Subjects: {}", value.subjects.iter().map(|source| match source {
                    collagene_core::planning::plan::Source::Input(id) => format!("input {id}"),
                    collagene_core::planning::plan::Source::Step(index) => index.checked_add(1).map(|number| format!("step {number}")).unwrap_or_else(|| "invalid step index (overflow)".into()),
                }).collect::<Vec<_>>().join(", "))}</p>
            })}
            {value.citations.iter().map(|citation| view! {
                <p class="helper">{format!("{} · {}", citation.source, citation.locator)}</p>
            }).collect_view()}
        </details>
    }
}

/// Keep the full typed evidence while rendering one bounded page of checks.
pub(crate) fn findings(values: &[CheckFinding]) -> impl IntoView + use<> {
    let values = StoredValue::new(values.to_vec());
    let page = RwSignal::new(0usize);
    let count = values.with_value(Vec::len);
    view! {
        <section class="verification-findings">
            <p class="helper">{format!("{} named checks · full evidence remains in the report", count)}</p>
            {move || values.with_value(|values| values.iter().skip(page.get()*25).take(25).map(finding).collect_view())}
            <button class="text-tool" disabled=move || page.get()==0
                on:click=move |_| page.update(|index| *index=index.saturating_sub(1))>"Previous checks"</button>
            <button class="text-tool" disabled=move || (page.get()+1)*25>=count
                on:click=move |_| page.update(|index| *index+=1)>"Next checks"</button>
        </section>
    }
}
