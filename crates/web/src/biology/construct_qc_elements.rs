//! Explicit element-chain and literal reference-copy declarations for construct QC.
use super::{construct_qc_model::Draft, construct_qc_view::feature_choice};
use crate::sources::identity::SourceInfo;
use collagene_core::construct_qc::{
    element_chain::{ChainElement, ElementChain, ElementRole},
    required_parts::{RequiredOrientation, RequiredPart},
};
use collagene_model::Molecule;
use leptos::prelude::*;

#[component]
fn TextField(label: &'static str, value: RwSignal<String>) -> impl IntoView {
    view! {
        <label>{label}<input maxlength="128" prop:value=move || value.get()
            on:input=move |e| value.set(event_target_value(&e))/></label>
    }
}

#[component]
fn ReferenceChoice(
    label: &'static str,
    value: RwSignal<String>,
    draft: RwSignal<Draft>,
) -> impl IntoView {
    view! {
        <label>{label}
            <select prop:value=move || value.get()
                on:change=move |e| value.set(event_target_value(&e))>
                <option value="">"No reference selected"</option>
                {move || {
                    let selected = value.get();
                    let missing = !selected.is_empty()
                        && draft.with(|d| !d.options.references.iter().any(|r| r.id == selected));
                    missing.then(|| view! {
                        <option value=selected.clone()>{format!("{selected} · reference unavailable")}</option>
                    })
                }}
                {move || draft.with(|d| d.options.references.iter().map(|reference| {
                    let id = reference.id.clone();
                    let label = id.clone();
                    view! {<option value=id>{label}</option>}
                }).collect_view())}
            </select>
        </label>
    }
}

fn optional_number(value: &str) -> Result<Option<usize>, String> {
    if value.trim().is_empty() {
        Ok(None)
    } else {
        value.trim().parse().map(Some).map_err(|_| {
            "Enter a nonnegative whole count or leave the optional field empty.".into()
        })
    }
}
fn annotation(value: RwSignal<String>) -> Result<usize, String> {
    value
        .get_untracked()
        .parse()
        .map_err(|_| "Choose a CDS annotation.".into())
}

#[component]
pub(super) fn ElementDeclarations(
    draft: RwSignal<Draft>,
    rows: Memo<Vec<(usize, String, bool)>>,
    get_molecule: Callback<(), Molecule>,
    source: Memo<SourceInfo>,
    epoch: RwSignal<u64>,
    status: RwSignal<String>,
) -> impl IntoView {
    let chain_id = RwSignal::new(String::new());
    let cds = RwSignal::new(String::new());
    let upstream = RwSignal::new(String::new());
    let downstream = RwSignal::new(String::new());
    let upstream_reference = RwSignal::new(String::new());
    let downstream_reference = RwSignal::new(String::new());
    let downstream_role = RwSignal::new("terminator".to_string());
    let max_upstream = RwSignal::new(String::new());
    let max_downstream = RwSignal::new(String::new());
    let part_id = RwSignal::new(String::new());
    let part_role = RwSignal::new(String::new());
    let part_reference = RwSignal::new(String::new());
    let part_minimum = RwSignal::new("1".to_string());
    let part_maximum = RwSignal::new(String::new());
    let part_orientation = RwSignal::new("either".to_string());
    Effect::new(move |_| {
        source.get();
        epoch.get();
        for field in [
            chain_id,
            cds,
            upstream,
            downstream,
            upstream_reference,
            downstream_reference,
            max_upstream,
            max_downstream,
            part_id,
            part_role,
            part_reference,
            part_maximum,
        ] {
            field.set(String::new());
        }
        downstream_role.set("terminator".into());
        part_minimum.set("1".into());
        part_orientation.set("either".into());
    });
    let add_chain = Callback::new(move |(): ()| {
        let result = (|| -> Result<(), String> {
            let element = |feature: RwSignal<String>, reference: RwSignal<String>, role| {
                let feature = feature.get_untracked();
                if feature.is_empty() {
                    if !reference.get_untracked().is_empty() {
                        return Err("Choose the annotation for the selected chain reference, or clear that reference.".into());
                    }
                    return Ok(None);
                }
                let feature_index = feature
                    .parse()
                    .map_err(|_| "Choose an element annotation.")?;
                let reference = reference.get_untracked();
                Ok::<_, String>(Some(ChainElement {
                    feature_index,
                    reference_id: (!reference.is_empty()).then_some(reference),
                    role,
                }))
            };
            let chain = ElementChain {
                id: chain_id.get_untracked(),
                cds_feature: annotation(cds)?,
                upstream: element(upstream, upstream_reference, ElementRole::Promoter)?,
                downstream: element(
                    downstream,
                    downstream_reference,
                    if downstream_role.get_untracked() == "polyadenylation" {
                        ElementRole::PolyadenylationSignal
                    } else {
                        ElementRole::Terminator
                    },
                )?,
                max_upstream_distance: optional_number(&max_upstream.get_untracked())?,
                max_downstream_distance: optional_number(&max_downstream.get_untracked())?,
            };
            let mut next = draft.get_untracked();
            next.options.element_chains.push(chain);
            next.options_for(&get_molecule.run(()))?;
            draft.set(next);
            chain_id.set(String::new());
            Ok(())
        })();
        if let Err(error) = result {
            status.set(error);
        }
    });
    let add_part = Callback::new(move |(): ()| {
        let result = (|| -> Result<(), String> {
            let part = RequiredPart {
                id: part_id.get_untracked(),
                role: part_role.get_untracked(),
                reference_id: part_reference.get_untracked(),
                minimum_copies: optional_number(&part_minimum.get_untracked())?
                    .ok_or("Enter a minimum literal copy count, including zero if appropriate.")?,
                maximum_copies: optional_number(&part_maximum.get_untracked())?,
                orientation: match part_orientation.get_untracked().as_str() {
                    "forward" => RequiredOrientation::Forward,
                    "reverse" => RequiredOrientation::Reverse,
                    _ => RequiredOrientation::Either,
                },
            };
            let mut next = draft.get_untracked();
            next.options.required_parts.push(part);
            next.options_for(&get_molecule.run(()))?;
            draft.set(next);
            part_id.set(String::new());
            Ok(())
        })();
        if let Err(error) = result {
            status.set(error);
        }
    });
    view! {
        <details>
            <summary>"Declared element chains"</summary>
            <p class="helper">"Check orientation, ordering and declared distance limits around a CDS. Annotation roles remain assumptions; literal reference identity does not establish promoter or termination function."</p>
            <TextField label="Chain ID" value=chain_id/>
            {feature_choice(rows, cds, "Chain CDS")}
            {feature_choice(rows, upstream, "Upstream promoter annotation (optional)")}
            <ReferenceChoice label="Upstream reference (optional)" value=upstream_reference draft/>
            <TextField label="Maximum upstream distance in bases (optional)" value=max_upstream/>
            {feature_choice(rows, downstream, "Downstream annotation (optional)")}
            <label>"Declared downstream role"
                <select prop:value=move || downstream_role.get()
                    on:change=move |e| downstream_role.set(event_target_value(&e))>
                    <option value="terminator">"Terminator"</option>
                    <option value="polyadenylation">"Polyadenylation signal"</option>
                </select>
            </label>
            <ReferenceChoice label="Downstream reference (optional)" value=downstream_reference draft/>
            <TextField label="Maximum downstream distance in bases (optional)" value=max_downstream/>
            <button class="text-tool" on:click=move |_| add_chain.run(())>"Add chain check"</button>
            {move || draft.with(|d| d.options.element_chains.iter().map(|chain| {
                let id = chain.id.clone();
                let remove_id = id.clone();
                let label = format!("{} · CDS {}", id, chain.cds_feature + 1);
                view! {<p>{label}<button class="text-tool" on:click=move |_| draft.update(|d| {
                    d.options.element_chains.retain(|c| c.id != remove_id);
                })>"Remove chain"</button></p>}
            }).collect_view())}
        </details>
        <details>
            <summary>"Literal reference copy requirements"</summary>
            <p class="helper">"Search the complete construct for a pinned reference sequence. A zero count is conclusive only within a completed literal search; uncertainty and work limits remain visible. The declared role is a label, not a functional classification."</p>
            <TextField label="Requirement ID" value=part_id/>
            <TextField label="Declared role label" value=part_role/>
            <ReferenceChoice label="Pinned reference" value=part_reference draft/>
            <TextField label="Minimum literal copies" value=part_minimum/>
            <TextField label="Maximum literal copies (optional)" value=part_maximum/>
            <label>"Allowed reference orientation"
                <select prop:value=move || part_orientation.get()
                    on:change=move |e| part_orientation.set(event_target_value(&e))>
                    <option value="either">"Either orientation"</option>
                    <option value="forward">"Forward"</option>
                    <option value="reverse">"Reverse"</option>
                </select>
            </label>
            <button class="text-tool" on:click=move |_| add_part.run(())>"Add copy-count check"</button>
            {move || draft.with(|d| d.options.required_parts.iter().map(|part| {
                let id = part.id.clone();
                let label = format!("{} · {} · reference {} · minimum {} · maximum {} · {:?}",
                    part.id, part.role, part.reference_id, part.minimum_copies,
                    part.maximum_copies.map(|v| v.to_string()).unwrap_or_else(|| "unspecified".into()),
                    part.orientation);
                view! {<p>{label}<button class="text-tool" on:click=move |_| draft.update(|d| {
                    d.options.required_parts.retain(|p| p.id != id);
                })>"Remove requirement"</button></p>}
            }).collect_view())}
            <details><summary>"Search limits"</summary>
                <p class="helper">"These limits are shared across required reference searches."</p>
                <label>"Work units"<input type="number" min="1" max="50000000"
                    prop:value=move || draft.with(|d| d.options.required_part_limits.max_work.to_string())
                    on:change=move |e| match event_target_value(&e).parse() {
                        Ok(v) => draft.update(|d| d.options.required_part_limits.max_work = v),
                        Err(_) => status.set("Enter a whole search work limit.".into()),
                    }/></label>
                <label>"Total retained matches"<input type="number" min="1" max="1024"
                    prop:value=move || draft.with(|d| d.options.required_part_limits.max_matches.to_string())
                    on:change=move |e| match event_target_value(&e).parse() {
                        Ok(v) => draft.update(|d| d.options.required_part_limits.max_matches = v),
                        Err(_) => status.set("Enter a whole retained-match limit.".into()),
                    }/></label>
                <label>"Matches per reference"<input type="number" min="1" max="256"
                    prop:value=move || draft.with(|d| d.options.required_part_limits.max_matches_per_reference.to_string())
                    on:change=move |e| match event_target_value(&e).parse() {
                        Ok(v) => draft.update(|d| d.options.required_part_limits.max_matches_per_reference = v),
                        Err(_) => status.set("Enter a whole per-reference match limit.".into()),
                    }/></label>
            </details>
        </details>
    }
}
