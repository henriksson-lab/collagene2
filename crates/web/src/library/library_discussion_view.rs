use crate::library::library_discussion_client::Client;
use crate::library::library_discussion_model::Review;
use collagene_library_protocol::{self as wire, discussion as api};
use leptos::prelude::*;
use std::{cell::Cell, rc::Rc};
#[derive(Clone)]
enum Command {
    Load { more: bool },
    Append(api::Action),
}
#[component]
pub fn LibraryDiscussion(target: api::Target, client: Client) -> impl IntoView {
    let target = StoredValue::new(target);
    let session = StoredValue::new_local(client);
    let epoch = session.with_value(|s| s.publication_epoch());
    let alive = StoredValue::new_local(Rc::new(Cell::new(true)));
    on_cleanup(move || alive.with_value(|a| a.set(false)));
    let review = RwSignal::new(None::<Review>);
    let body = RwSignal::new(String::new());
    let reply_to = RwSignal::new(None::<u64>);
    let busy = RwSignal::new(false);
    let status = RwSignal::new(String::new());
    let send = Callback::new(move |command: Command| {
        if busy.get_untracked() {
            return;
        }
        let client = session.get_value();
        if client.publication_epoch() != epoch {
            status.set("The member session changed. Reopen revision history.".into());
            return;
        }
        let Some(lab) = client.lab_id() else {
            status.set("Connect to the lab library first.".into());
            return;
        };
        let previous = review.get_untracked();
        let operation = match &command {
            Command::Load { more } => {
                let (offset, expected_generation) = if *more {
                    let Some(p) = &previous else {
                        return;
                    };
                    let Some(offset) = p.next_offset else {
                        return;
                    };
                    (offset, Some(p.generation))
                } else {
                    (0, None)
                };
                api::Request::List {
                    target: target.get_value(),
                    offset,
                    limit: 200,
                    expected_generation,
                }
            }
            Command::Append(action) => {
                let Some(p) = &previous else {
                    return;
                };
                match p.append_request(action.clone()) {
                    Ok(r) => r,
                    Err(e) => {
                        status.set(e);
                        return;
                    }
                }
            }
        };
        busy.set(true);
        status.set(String::new());
        leptos::task::spawn_local(async move {
            let mut acknowledged = false;
            let result: Result<Review, String> = async {
                let mut operation = operation;
                let mut response = client.request(operation.clone()).await?;
                if !alive.with_value(|a| a.get()) || client.publication_epoch() != epoch {
                    return Err("Discussion session changed.".into());
                }
                if let api::Outcome::Failure(f) = &response.outcome {
                    return Err(failure(f));
                }
                let previous = if matches!(command, Command::Append(_)) {
                    // Clear only after native acknowledgment; failed writes keep the draft.
                    acknowledged = true;
                    body.set(String::new());
                    reply_to.set(None);
                    operation = api::Request::List {
                        target: target.get_value(),
                        offset: 0,
                        limit: 200,
                        expected_generation: None,
                    };
                    response = client.request(operation.clone()).await?;
                    if let api::Outcome::Failure(f) = &response.outcome {
                        return Err(failure(f));
                    }
                    None
                } else {
                    previous.as_ref()
                };
                Review::page(
                    previous,
                    &api::RequestEnvelope {
                        request_id: response.request_id,
                        request: operation,
                    },
                    &response,
                    &lab,
                )
            }
            .await;
            if !alive.with_value(|a| a.get()) || client.publication_epoch() != epoch {
                return;
            }
            match result {
                Ok(p) => {
                    review.set(Some(p));
                    status.set("Discussion loaded for this exact revision.".into());
                }
                Err(e) => {
                    review.update(|p| {
                        if let Some(p) = p {
                            p.fresh = false;
                        }
                    });
                    status.set(if acknowledged {
                        format!("Saved successfully, but refreshing the discussion failed: {e}")
                    } else {
                        e
                    });
                }
            }
            busy.set(false);
        });
    });
    let post = Callback::new(move |decision: Option<api::ReviewDecision>| {
        let text = match wire::Text::new(body.get_untracked()) {
            Ok(t) => t,
            Err(_) => {
                status.set("Comment exceeds 8,192 UTF-8 bytes.".into());
                return;
            }
        };
        let action = if let Some(decision) = decision {
            api::Action::Review {
                decision,
                body: text,
            }
        } else if let Some(thread) = reply_to.get_untracked() {
            if !review.with_untracked(|r| r.as_ref().is_some_and(|r| r.can_reply(thread))) {
                status.set("Reload or reopen the thread before replying.".into());
                return;
            }
            api::Action::Reply { thread, body: text }
        } else {
            api::Action::OpenThread { body: text }
        };
        send.run(Command::Append(action));
    });
    view! {
        <details class="library-discussion">
            <summary>"Comments and human review"</summary>
            <p class="helper">"Applies only to this saved revision. Human approval does not replace molecular verification."</p>
            <button class="text-tool" disabled=move||busy.get() on:click=move|_|send.run(Command::Load{more:false})>"Load / refresh discussion"</button>
            {move||review.get().map(|r|{
                let label=format!("{} · {:?} · discussion generation {}",r.member.author.display_name.as_str(),r.member.role,r.generation);
                let reviews=if r.fresh&&r.next_offset.is_none(){r.latest_reviews()}else{Vec::new()};
                view!{
                    <p class="helper">{label}</p>
                    {(!r.fresh).then(||view!{<p class="helper">"Discussion changed or the request failed. Refresh before posting; your draft is retained."</p>})}
                    <ul>{reviews.into_iter().map(|e|{
                        let api::Action::Review{decision,..}=e.action else{unreachable!()};
                        view!{<li>{format!("{}: {:?}",e.author.display_name.as_str(),decision)}</li>}
                    }).collect_view()}</ul>
                    <ol class="library-discussion-events">{r.events.iter().cloned().map(|event|{
                        let sequence=event.sequence;
                        let is_root=matches!(event.action,api::Action::OpenThread{..});
                        let label=format!("#{} · {} · {}",sequence,event.author.display_name.as_str(),action_label(&event.action));
                        let text=match event.action{api::Action::OpenThread{body}|api::Action::Reply{body,..}|api::Action::Review{body,..}=>body.into_string(),_=>String::new()};
                        view!{<li><strong>{label}</strong><p class="discussion-body">{text}</p>
                            {is_root.then(||view!{
                                <button class="text-tool" disabled=move||busy.get()||!review.with(|r|r.as_ref().is_some_and(|r|r.can_reply(sequence))) on:click=move|_|reply_to.set(Some(sequence))>"Reply"</button>
                                <button class="text-tool" disabled=move||busy.get()||!review.with(|r|r.as_ref().is_some_and(|r|r.can_resolve(sequence))) on:click=move|_|{
                                    if let Some(resolved)=review.with_untracked(|r|r.as_ref().and_then(|r|r.resolved(sequence))){send.run(Command::Append(api::Action::SetResolved{thread:sequence,resolved:!resolved}));}
                                }>{move||if review.with(|r|r.as_ref().and_then(|r|r.resolved(sequence)))==Some(true){"Reopen"}else{"Resolve"}}</button>
                            })}
                        </li>}
                    }).collect_view()}</ol>
                    {r.next_offset.is_some().then(||view!{<button class="text-tool" disabled=move||busy.get() on:click=move|_|send.run(Command::Load{more:true})>"Load more events"</button><p class="helper">"Load all current events before posting."</p>})}
                }
            })}
            <form on:submit=move|e|{e.prevent_default();post.run(None);}>
                <label>{move||reply_to.get().map(|n|format!("Reply to thread #{n}")).unwrap_or_else(||"New comment / review note".into())}
                    <textarea maxlength="8192" prop:value=move||body.get() on:input=move|e|body.set(event_target_value(&e)) disabled=move||busy.get()/>
                </label>
                <Show when=move||reply_to.get().is_some()><button type="button" class="text-tool" on:click=move|_|reply_to.set(None)>"Cancel reply"</button></Show>
                <button type="submit" class="text-tool" disabled=move||busy.get()||body.with(|b|b.trim().is_empty())||!review.with(|r|r.as_ref().is_some_and(Review::can_write))>"Post comment"</button>
                <button type="button" class="text-tool" disabled=move||busy.get()||!review.with(|r|r.as_ref().is_some_and(Review::can_write)) on:click=move|_|post.run(Some(api::ReviewDecision::Approved))>"Approve revision"</button>
                <button type="button" class="text-tool" disabled=move||busy.get()||!review.with(|r|r.as_ref().is_some_and(Review::can_write)) on:click=move|_|post.run(Some(api::ReviewDecision::ChangesRequested))>"Request changes"</button>
                <button type="button" class="text-tool" disabled=move||busy.get()||!review.with(|r|r.as_ref().is_some_and(Review::can_write)) on:click=move|_|post.run(Some(api::ReviewDecision::Withdrawn))>"Withdraw my review"</button>
            </form>
            <p role="status" aria-live="polite">{move||status.get()}</p>
        </details>
    }
}
fn failure(f: &api::Failure) -> String {
    match f.code {
        wire::ErrorCode::Conflict => format!(
            "Discussion changed (generation {}). Refresh and review before resubmitting.",
            f.actual_generation.unwrap_or(0)
        ),
        wire::ErrorCode::Forbidden => "Your member identity or role does not permit this action. Refresh the review before posting.".into(),
        wire::ErrorCode::NotFound => {
            "This exact revision or thread is unavailable in this lab.".into()
        }
        wire::ErrorCode::TooLarge => "The discussion reached a storage bound.".into(),
        _ => "Discussion request failed. Refresh before retrying; your draft is retained.".into(),
    }
}
fn action_label(action: &api::Action) -> String {
    match action {
        api::Action::OpenThread { .. } => "opened thread".into(),
        api::Action::Reply { thread, .. } => format!("replied to #{thread}"),
        api::Action::SetResolved { thread, resolved } => format!(
            "{} #{thread}",
            if *resolved { "resolved" } else { "reopened" }
        ),
        api::Action::Review { decision, .. } => format!("review: {decision:?}"),
    }
}
