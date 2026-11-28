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

/// Format advice is deliberately separate from status-changing cloning checks.
pub(crate) fn standard_review(
    review: Option<&collagene_core::cloning::assembly::standards::declaration::DeclaredReview>,
) -> impl IntoView + use<> {
    use collagene_core::cloning::assembly::standards::declaration::DeclaredReviewState;
    let content = review.map(|review| {
        let body = match &review.state {
            DeclaredReviewState::Unavailable { reasons } => view! {
                <p>"Source-format advice unavailable; no scoped pass is established."</p>
                {reasons.iter().map(|reason| view! {
                    <p>{format!("Part {:?}: {:?} · {}", reason.part_index.map(|i|i+1), reason.kind, reason.message)}</p>
                }).collect_view()}
            }.into_any(),
            DeclaredReviewState::Assessed(report) => {
                use collagene_assembly_standards::assessment::{FindingKind, Outcome};
                let outcome = match report.assessment.outcome {
                    Outcome::Conforms => "Scoped source format: pass",
                    Outcome::Violations => "Scoped source format: fail",
                    Outcome::Indeterminate => "Scoped source format: unresolved",
                };
                let findings = StoredValue::new(report.assessment.findings.clone());
                let sources = StoredValue::new(report.sources.clone());
                let page = RwSignal::new(0usize);
                let count = findings.with_value(Vec::len);
                view! {
                    <p><strong>{outcome}</strong>{format!(
                        " · complete {} · observed violations {}, unknowns {}, advisories {}",
                        report.assessment.complete, report.assessment.violations_observed,
                        report.assessment.indeterminate_observed, report.assessment.advisories_observed,
                    )}</p>
                    <p>{report.profile.scope.clone()}</p>
                    <p class="helper">{format!(
                        "{} · profile {} · data {} · catalogue {} / {}",
                        report.profile.name, report.profile.profile_identity, report.profile.definition.revision,
                        report.profile.catalog_pin.edition, report.profile.catalog_pin.content_hash,
                    )}</p>
                    <p class="helper">{format!("Request {} · report {} · budget limits {:?}", report.context_revision, report.report_revision, report.assessment.budgets_hit)}</p>
                    <details><summary>"Exact ordered source and body declarations"</summary>
                        {report.sources.iter().zip(&report.assessment.parts).enumerate().map(|(index, (source, part))| view! {
                            <p>{format!(
                                "{}: {} · slot {} · source revision {} · source interval start {}, length {}, reverse {} · oriented body [{}..{})",
                                index+1, source.id, part.slot, source.revision, source.start, source.length,
                                source.reverse, part.body_start, part.body_end,
                            )}</p>
                        }).collect_view()}
                    </details>
                    <p>{format!("{} retained findings; full evidence is retained in the report/protocol.", count)}</p>
                    {move || findings.with_value(|items| items.iter().skip(page.get()*25).take(25).map(|finding| {
                        let coordinates = match (&finding.kind, finding.part_index) {
                            (FindingKind::AvoidedSite { start, end, .. }, Some(index)) => sources.with_value(|sources| {
                                sources.get(index)
                                    .and_then(|source| source.map_interval(*start, *end).ok())
                                    .map(|spans| format!("Original source spans (one-based inclusive): {:?}",
                                        spans.into_iter().map(|(start,end)| (start+1,end)).collect::<Vec<_>>()))
                                    .unwrap_or_else(|| "Original source coordinate mapping unavailable.".into())
                            }),
                            _ => String::new(),
                        };
                        view! { <p>{format!("Part {:?}: {:?} · {:?}. {}", finding.part_index.map(|i|i+1), finding.severity, finding.kind, coordinates)}</p> }
                    }).collect_view())}
                    <button class="text-tool" disabled=move || page.get()==0 on:click=move |_| page.update(|p| *p=p.saturating_sub(1))>"Previous format findings"</button>
                    <button class="text-tool" disabled=move || (page.get()+1)*25>=count on:click=move |_| page.update(|p| *p+=1)>"Next format findings"</button>
                    <p class="helper">{format!("Allowed next declared roles {:?}; may finish {}. {}", report.successors.allowed_slots, report.successors.may_finish, report.successors.scope)}</p>
                    <details><summary>"Scope, source provenance and limitations"</summary>
                        {report.profile.definition.limitations.iter().map(|text| view! { <p>{text.clone()}</p> }).collect_view()}
                        {report.profile.definition.sources.iter().map(|source| view! {
                            <p>{format!("{} · {} · revision {} · {} · rights {} · capture SHA256 {:?}", source.id, source.url, source.revision, source.locator, source.source_license, source.upstream_sha256)}</p>
                        }).collect_view()}
                    </details>
                }.into_any()
            }
        };
        view! {
            <details class="verification-source-format">
                <summary>"Declared source-part format — advisory"</summary>
                <p class="helper">{review.scope.clone()}</p>
                <p class="helper">{format!("Declaration {} · review {} · comparison charge {:?}", review.declaration_revision, review.review_revision, review.charge)}</p>
                {body}
            </details>
        }
    });
    view! { {content} }
}

pub(crate) fn report(
    report: &collagene_core::verification::evidence::VerificationReport,
) -> impl IntoView + use<> {
    view! {
        {standard_review(report.standard_review.as_ref())}
        {findings(&report.findings)}
    }
}
