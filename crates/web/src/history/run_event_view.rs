//! Frozen, read-only navigation into a version retained by a cloning run.
//! A separate host-owned action creates an editable fork; viewing never replays
//! the design, changes a document head or restores a scientific assessment.
use crate::history::run_event_model::Request;
use crate::history::run_event_model::Role;
use crate::history::run_event_model::Snapshot;
use leptos::prelude::*;
use std::sync::Arc;

fn role_label(role: &Role) -> String {
    match role {
        Role::Target => "Intended target".into(),
        Role::Input { input_id } => format!("Source: {input_id}"),
        Role::Product { step_index } => format!("Product from step {}", step_index + 1),
    }
}

#[component]
pub fn RunEventPanel(
    snapshot: Arc<Snapshot>,
    epoch: RwSignal<u64>,
    on_fork: Callback<Request, Result<(), String>>,
    on_back: Callback<String>,
) -> impl IntoView {
    let captured_epoch = snapshot.request.epoch;
    let current = Memo::new(move |_| epoch.get() == captured_epoch);
    let title = snapshot.molecule.name.clone();
    let role = role_label(&snapshot.request.role);
    let identity = format!(
        "run-version:{}:{}:{}:{}:{}",
        snapshot.request.run_id,
        role,
        snapshot.request.binding.document_id,
        snapshot.request.binding.event_id,
        snapshot.request.binding.revision,
    );
    let product = snapshot.molecule.clone();
    let snapshot = StoredValue::new(snapshot);
    let status = RwSignal::new(String::new());
    let fork = move |_| {
        if epoch.get_untracked() != captured_epoch {
            status.set(
                "The workspace changed. Reopen the cloning run before forking this version.".into(),
            );
            return;
        }
        let request = snapshot.with_value(|snapshot| snapshot.request.clone());
        if let Err(error) = on_fork.run(request) {
            status.set(error);
        }
    };
    let back = move |_| {
        if epoch.get_untracked() != captured_epoch {
            status.set(
                "The workspace changed. Reopen the cloning run from its retained history.".into(),
            );
            return;
        }
        on_back.run(snapshot.with_value(|snapshot| snapshot.request.run_id.clone()));
    };
    view! {
        <section class="run-event-view" aria-label="Recorded cloning version">
            <header style="overflow-wrap:anywhere">
                <p class="eyebrow">"Recorded version · read-only"</p>
                <h2>{title}</h2>
                <p>{role}</p>
            </header>
            <p class="helper">"This is the exact reference version retained by the cloning run. Viewing it does not change the original document or establish a current verification result."</p>
            <p class="helper">{snapshot.with_value(|snapshot|if snapshot.current_head {
                "This version was the document head when opened. The viewer remains a frozen snapshot."
            } else {
                "The original document has a different head. Its current version remains unchanged."
            })}</p>
            <details style="overflow-wrap:anywhere">
                <summary>"Recorded source details"</summary>
                <dl>{snapshot.with_value(|snapshot|view!{
                    <dt>"Cloning run"</dt><dd>{snapshot.request.run_id.clone()}</dd>
                    <dt>"Document"</dt><dd>{snapshot.request.binding.document_id.clone()}</dd>
                    <dt>"Version event"</dt><dd>{snapshot.request.binding.event_id.clone()}</dd>
                    <dt>"Content revision"</dt><dd>{snapshot.request.binding.revision.clone()}</dd>
                })}</dl>
            </details>
            <div class="method-row-actions">
                <button class="button secondary" disabled=move||!current.get() on:click=back>"Back to cloning run"</button>
                <button class="button primary" disabled=move||!current.get() on:click=fork>"Fork editable copy"</button>
            </div>
            <p class="helper">"Forking creates a separate document from this recorded version and preserves its ancestry. It does not replace the original head."</p>
            <p class="helper" role="status">{move||if current.get(){status.get()}else{"Workspace context changed. Reopen this version from its cloning run.".into()}}</p>
            <crate::assembly::review::assembly_product_view::ReferenceCanvas product identity current on_message=Callback::new(move|message|status.set(message))/>
        </section>
    }.into_any()
}
