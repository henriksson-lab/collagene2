//! Reversible working-draft replacement. This stores declarations, never receipts.
use crate::assembly_draft_model::Draft;
use leptos::prelude::*;
use std::sync::Arc;

#[derive(Clone)]
struct Pending {
    current: Arc<Draft>,
    previous: Option<Arc<Draft>>,
    epoch: u64,
    restore: bool,
}

#[component]
pub fn AssemblyDraftActions(
    draft: RwSignal<Draft>,
    previous: RwSignal<Option<Arc<Draft>>>,
    epoch: RwSignal<u64>,
    locked: Memo<bool>,
    on_replace: Callback<Draft>,
) -> impl IntoView {
    let pending = RwSignal::new(None::<Pending>);
    let message = RwSignal::new(String::new());
    let propose = Callback::new(move |restore: bool| {
        if locked.get_untracked() {
            return;
        }
        pending.set(Some(Pending {
            current: Arc::new(draft.get_untracked()),
            previous: previous.get_untracked(),
            epoch: epoch.get_untracked(),
            restore,
        }));
        message.set(String::new());
    });
    let apply = Callback::new(move |saved: Pending| {
        if locked.get_untracked() {
            return;
        }
        if epoch.get_untracked() != saved.epoch
            || !draft.with_untracked(|current| current == saved.current.as_ref())
            || previous.get_untracked() != saved.previous
        {
            pending.set(None);
            message.set("The working draft changed. Choose the action again to review its current contents.".into());
            return;
        }
        let next = if saved.restore {
            let Some(value) = saved.previous else {
                return;
            };
            value.as_ref().clone()
        } else {
            Draft::default()
        };
        batch(|| {
            previous.set(Some(saved.current));
            on_replace.run(next);
            pending.set(None);
            message.set(
                if saved.restore {
                    "Previous draft restored; review intent and recompute its evidence."
                } else {
                    "New assembly draft started."
                }
                .into(),
            );
        });
    });
    view! {
        <details class="assembly-draft-menu" data-transient-menu="">
            <summary>"Draft ▾"</summary>
            <div class="assembly-draft-menu-content">
                <button class="text-tool" data-menu-keep-open="" disabled=move||locked.get() on:click=move|_|propose.run(false)>"New assembly…"</button>
                <button class="text-tool" data-menu-keep-open="" disabled=move||locked.get()||previous.with(Option::is_none) on:click=move|_|propose.run(true)>"Switch to previous draft…"</button>
                <p class="helper">"The working draft is saved with the workspace. One previous draft is kept for this session, until another draft replaces it, a workspace update clears it, or this session ends."</p>
                {move||pending.get().map(|saved|{
                    let title=if saved.restore { "Switch working drafts?" } else { "Start a new assembly?" };
                    let action=if saved.restore { "Switch drafts" } else { "Start new assembly" };
                    let name=saved.current.name.clone();
                    let count=saved.current.slots.iter().filter(|s|s.captured.is_some()).count();
                    let restore_name=saved.previous.as_ref().map(|d|d.name.clone());
                    view!{<section class="assembly-draft-confirm" aria-label=title><h3>{title}</h3>
                        <p>{format!("Current draft: {name} · {count} captured fragments.")}</p>
                        {saved.restore.then(||view!{<p>{format!("Restore: {}",restore_name.unwrap_or_default())}</p>})}
                        <p class="helper">"Your current draft becomes the previous draft. Saved constructs and recorded plans stay in the workspace; preparation and verification will be recomputed."</p>
                        <div class="method-row-actions"><button class="button primary" disabled=move||locked.get() on:click=move|_|apply.run(saved.clone())>{action}</button><button class="text-tool" on:click=move|_|pending.set(None)>"Cancel"</button></div>
                    </section>}
                })}
                <p class="assembly-draft-action-status" role="status">{move||message.get()}</p>
            </div>
        </details>
    }
}
