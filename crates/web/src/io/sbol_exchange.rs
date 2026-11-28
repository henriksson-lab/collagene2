//! Explicit SBOL design-graph review. Native export must not erase restrictive
//! physical/reference-only source evidence; imported SBOL remains a reference.
use crate::io::download::ExportArtifact;
use crate::io::download::prepare_download;
use collagene_core::documents::workspace::WorkspaceAttachment;
use collagene_design_exchange::{
    self as exchange, Document, Export, ExportIds, FeatureKind, Issue, Kind, Projection,
    ProjectionPolicy, Selection,
};
use collagene_model::{Molecule, Topology};
use leptos::prelude::*;
use std::sync::Arc;
pub(crate) const ATTACHMENT_KIND: &str = "sbol-design-exchange-v1";
pub(crate) const MEDIA_TYPE: &str = "application/ld+json";
/// Generic JSON must be a successfully decoded supported SBOL graph, not merely
/// mention an SBOL URL in arbitrary text. Explicit .sbol inputs can show errors.
pub(crate) fn looks_like_sbol(bytes: &[u8]) -> bool {
    exchange::import(bytes).is_ok_and(|d| {
        d.objects()
            .values()
            .any(|o| o.classes.iter().any(|c| c.starts_with(exchange::SBOL)))
    })
}
fn policy(allow: bool) -> ProjectionPolicy {
    if allow {
        ProjectionPolicy::AllowOmissions
    } else {
        ProjectionPolicy::Strict
    }
}
fn selection(component: String, sequence: String, topology: String) -> Result<Selection, String> {
    if component.is_empty() || sequence.is_empty() {
        return Err("Choose an exact Component and one of its linked Sequences.".into());
    }
    Ok(Selection {
        component,
        sequence,
        topology: match topology.as_str() {
            "linear" => Some(Topology::Linear),
            "circular" => Some(Topology::Circular),
            "" => None,
            _ => return Err("Choose a supported topology policy.".into()),
        },
    })
}
fn complete(
    result: Result<ExportArtifact, String>,
    check: impl Fn() -> Result<(), String> + Clone + 'static,
    status: RwSignal<String>,
    message: &'static str,
) {
    match result {
        Ok(artifact) => {
            let finished = check.clone();
            artifact.complete_guarded(check, move |result| {
                if status.try_get_untracked().is_none() || finished().is_err() {
                    return;
                }
                status.set(match result {
                    Ok(()) => message.into(),
                    Err(error) => error,
                });
            });
        }
        Err(error) => status.set(error),
    }
}
fn object_label(doc: &Document, id: &str) -> String {
    let name = doc
        .objects()
        .get(id)
        .and_then(|o| o.identity.name.as_ref().or(o.identity.display_id.as_ref()));
    match name {
        Some(name) => format!("{} · {}", name.chars().take(96).collect::<String>(), id),
        None => id.into(),
    }
}
#[component]
fn IriPicker(
    label: &'static str,
    choices: Memo<Vec<(String, String)>>,
    selected: RwSignal<String>,
) -> impl IntoView {
    let filter = RwSignal::new(String::new());
    let page = RwSignal::new(0usize);
    let filtered = Memo::new(move |_| {
        let query = filter.get().to_lowercase();
        choices
            .get()
            .into_iter()
            .filter(|(id, label)| {
                query.is_empty()
                    || id.to_lowercase().contains(&query)
                    || label.to_lowercase().contains(&query)
            })
            .collect::<Vec<_>>()
    });
    Effect::new(move |_| {
        let _ = filter.get();
        let _ = choices.get();
        page.set(0);
    });
    view! {<div class="sbol-iri-picker"><label>{format!("Find {label}")}<input type="search" prop:value=move||filter.get() on:input=move|e|filter.set(event_target_value(&e))/></label><label>{label}<select prop:value=move||selected.get() on:change=move|e|selected.set(event_target_value(&e))><option value="">"Choose explicitly…"</option>{move||{let all=filtered.get();let count=all.len();let start=page.get().min(count.div_ceil(100).max(1)-1)*100;let shown=all.into_iter().skip(start).take(100).collect::<Vec<_>>();let current=selected.get();let mut rows=Vec::new();if !current.is_empty()&&!shown.iter().any(|(id,_)|id==&current){rows.push(view!{<option value=current.clone()>{current.clone()}</option>});}rows.extend(shown.into_iter().map(|(id,label)|view!{<option value=id>{label}</option>}));rows.collect_view()}}</select></label><p class="helper">{move||format!("{} matching choices · page {} of {}",filtered.with(Vec::len),page.get()+1,filtered.with(|v|v.len().div_ceil(100).max(1)))}</p><div class="method-row-actions"><button class="text-tool" disabled={move||page.get()==0} on:click=move|_|page.update(|v|*v=v.saturating_sub(1))>"Previous choices"</button><button class="text-tool" disabled={move||(page.get()+1)*100>=filtered.with(Vec::len)} on:click=move|_|page.update(|v|*v+=1)>"Next choices"</button></div></div>}
}
#[component]
fn Issues(issues: Arc<Vec<Issue>>) -> impl IntoView {
    let page = RwSignal::new(0usize);
    let total = issues.len();
    let issues = StoredValue::new(issues);
    view! {<div class="sbol-loss-review"><p>{format!("{total} retained-data limitations / representation losses")}</p>{move||issues.with_value(|all|all.iter().skip(page.get()*40).take(40).map(|i|view!{<article><strong>{i.code.clone()}</strong><p style="overflow-wrap:anywhere">{i.subject.clone()}</p><p>{i.message.clone()}</p></article>}).collect_view())}<div class="method-row-actions"><button class="text-tool" disabled={move||page.get()==0} on:click=move|_|page.update(|v|*v=v.saturating_sub(1))>"Previous notices"</button><button class="text-tool" disabled={move||(page.get()+1)*40>=total} on:click=move|_|page.update(|v|*v+=1)>"Next notices"</button></div></div>}
}
#[component]
fn GraphReview(doc: Arc<Document>) -> impl IntoView {
    let page = RwSignal::new(0usize);
    let counts = doc.objects().values().fold([0usize; 6], |mut c, o| {
        c[match o.kind {
            Kind::Component(_) => 0,
            Kind::Sequence(_) => 1,
            Kind::Feature(_) => 2,
            Kind::Location(_) => 3,
            Kind::Constraint(_) => 4,
            _ => 5,
        }] += 1;
        c
    });
    let total = doc.objects().len();
    let issues = Arc::new(doc.issues().to_vec());
    let doc = StoredValue::new(doc);
    view! {<details><summary>"Retained design hierarchy and relationships"</summary><p>{format!("{} Components · {} Sequences · {} Features · {} Locations · {} Constraints · {} Collections/other objects",counts[0],counts[1],counts[2],counts[3],counts[4],counts[5])}</p><p class="helper">"Relationships are design declarations. Constraints, provenance and imported reports are not current cloning verification or physical material."</p>{move||doc.with_value(|d|d.objects().iter().skip(page.get()*40).take(40).map(|(id,o)|{let summary=match &o.kind{Kind::Component(c)=>format!("Component: {} sequence links, {} owned features, {} constraints",c.sequences.len(),c.features.len(),c.constraints.len()),Kind::Sequence(s)=>format!("Sequence: {} explicit characters",s.elements.as_ref().map_or(0,String::len)),Kind::Feature(f)=>match &f.kind{FeatureKind::SequenceFeature=>format!("SequenceFeature: {} locations",f.locations.len()),FeatureKind::SubComponent{instance_of,source_locations,..}=>format!("SubComponent → {instance_of}; {} placements / {} source locations",f.locations.len(),source_locations.len())},Kind::Location(l)=>format!("Location {:?} on {} · order {:?}",l.kind,l.sequence,l.order),Kind::Constraint(c)=>format!("Constraint {} → {} via {}",c.subject,c.object,c.restriction),Kind::Collection{members,..}=>format!("Collection: {} members",members.len()),Kind::Unsupported=>"Unsupported object retained without interpretation".into()};view!{<article><p style="overflow-wrap:anywhere">{id.clone()}</p><p style="overflow-wrap:anywhere">{summary}</p></article>}}).collect_view())}<div class="method-row-actions"><button class="text-tool" disabled={move||page.get()==0} on:click=move|_|page.update(|v|*v=v.saturating_sub(1))>"Previous graph objects"</button><button class="text-tool" disabled={move||(page.get()+1)*40>=total} on:click=move|_|page.update(|v|*v+=1)>"Next graph objects"</button></div><Issues issues/></details>}
}
#[derive(Clone)]
struct ProjectionReview {
    selection: Selection,
    policy: ProjectionPolicy,
    value: Arc<Projection>,
}
#[derive(Clone)]
struct ExportReview {
    ids: ExportIds,
    policy: ProjectionPolicy,
    value: Arc<Export>,
}
#[component]
pub fn SbolExchangeDialog(
    bytes: Option<Arc<Vec<u8>>>,
    current: Arc<Molecule>,
    document_id: String,
    epoch: RwSignal<u64>,
    guard: Callback<(), Result<(), String>>,
    on_import: Callback<(Molecule, Option<WorkspaceAttachment>), Result<(), String>>,
    on_close: Callback<()>,
) -> impl IntoView {
    let status = RwSignal::new(String::new());
    let captured_epoch = epoch.get_untracked();
    let checked = Callback::new(move |()| {
        if status.try_get_untracked().is_none() {
            return Err("SBOL review closed.".into());
        }
        if epoch.get_untracked() != captured_epoch {
            return Err(
                "The workspace changed. Reopen SBOL exchange for the current document.".into(),
            );
        }
        guard.run(())
    });
    let parsed = bytes.as_ref().map(|b| exchange::import(b));
    let graph = StoredValue::new(match parsed {
        Some(Ok(d)) => Some(Arc::new(d)),
        Some(Err(e)) => {
            status.set(e.to_string());
            None
        }
        None => None,
    });
    let original = StoredValue::new(bytes);
    let source_present = original.with_value(Option::is_some);
    let current = StoredValue::new(current);
    let component = RwSignal::new(String::new());
    let sequence = RwSignal::new(String::new());
    let topology = RwSignal::new(String::new());
    let allow_projection_loss = RwSignal::new(false);
    let preview = RwSignal::new(None::<Arc<ProjectionReview>>);
    let components = Memo::new(move |_| {
        graph.with_value(|doc| {
            doc.as_ref()
                .map(|d| {
                    d.objects()
                        .iter()
                        .filter(|(_, o)| matches!(o.kind, Kind::Component(_)))
                        .map(|(id, _)| (id.clone(), object_label(d, id)))
                        .collect::<Vec<_>>()
                })
                .unwrap_or_default()
        })
    });
    let sequences = Memo::new(move |_| {
        let id = component.get();
        graph.with_value(|doc| {
            doc.as_ref()
                .and_then(|d| d.objects().get(&id).map(|o| (d, o)))
                .and_then(|(d, o)| {
                    if let Kind::Component(c) = &o.kind {
                        Some(
                            c.sequences
                                .iter()
                                .map(|id| (id.clone(), object_label(d, id)))
                                .collect::<Vec<_>>(),
                        )
                    } else {
                        None
                    }
                })
                .unwrap_or_default()
        })
    });
    Effect::new(move |_| {
        let _ = component.get();
        sequence.set(String::new());
        topology.set(String::new());
        preview.set(None);
    });
    Effect::new(move |_| {
        let _ = sequence.get();
        let _ = topology.get();
        let _ = allow_projection_loss.get();
        preview.set(None);
    });
    let get_selection = move || {
        selection(
            component.get_untracked(),
            sequence.get_untracked(),
            topology.get_untracked(),
        )
    };
    let review = move |_| {
        let result = (|| {
            checked.run(())?;
            let chosen = get_selection()?;
            let p = policy(allow_projection_loss.get_untracked());
            let doc = graph
                .get_value()
                .ok_or("Load a supported SBOL graph first.")?;
            let value = doc.project(&chosen, p).map_err(|e| e.to_string())?;
            checked.run(())?;
            Ok::<_, String>(Arc::new(ProjectionReview {
                selection: chosen,
                policy: p,
                value: Arc::new(value),
            }))
        })();
        match result {
            Ok(value) => {
                preview.set(Some(value));
                status.set("Review this explicit sequence projection and every listed loss before opening it.".into());
            }
            Err(e) => status.set(e),
        }
    };
    let use_projection = move |_| {
        let result = (|| {
            checked.run(())?;
            let saved = preview
                .get_untracked()
                .ok_or("Review a projection first.")?;
            if get_selection()? != saved.selection
                || policy(allow_projection_loss.get_untracked()) != saved.policy
            {
                return Err(
                    "Selection or loss policy changed. Review the projection again.".into(),
                );
            }
            let doc = graph
                .get_value()
                .ok_or("The retained source graph is unavailable.")?;
            let attachment = WorkspaceAttachment {
                kind: ATTACHMENT_KIND.into(),
                name: "original-design.sbol.jsonld".into(),
                media_type: MEDIA_TYPE.into(),
                bytes: doc.original_bytes().to_vec(),
            };
            let molecule = saved.value.molecule.clone();
            checked.run(())?;
            on_import.run((molecule, Some(attachment)))
        })();
        if status.try_get_untracked().is_some() {
            status.set(match result{Ok(())=>"Opened the selected sequence reference and retained its original SBOL graph. No physical or current verification authority was imported.".into(),Err(e)=>e});
        }
    };
    let graph_download = Callback::new(move |expanded: bool| {
        let fingerprint = (
            component.get_untracked(),
            sequence.get_untracked(),
            topology.get_untracked(),
        );
        let result = (|| {
            checked.run(())?;
            if expanded {
                let doc = graph
                    .get_value()
                    .ok_or("No validated retained graph is available.")?;
                let bytes = doc.expanded_bytes().map_err(|e| e.to_string())?;
                prepare_download("retained-design.expanded.jsonld", MEDIA_TYPE, &bytes)
            } else {
                let bytes = original.get_value().ok_or("No source file is loaded.")?;
                prepare_download("original-design.sbol.jsonld", MEDIA_TYPE, &bytes)
            }
        })();
        complete(
            result,
            move || {
                checked.run(())?;
                if fingerprint
                    != (
                        component.get_untracked(),
                        sequence.get_untracked(),
                        topology.get_untracked(),
                    )
                {
                    return Err("The selected SBOL context changed.".into());
                }
                Ok(())
            },
            status,
            if expanded {
                "Saved the complete retained expanded design graph; no unsupported relationships were flattened."
            } else {
                "Saved the exact original SBOL file."
            },
        );
    });
    let default = current.with_value(|m| {
        let stable = blake3::hash(document_id.as_bytes()).to_hex().to_string();
        let ns = "https://example.org/collagene".to_string();
        let root = format!("{ns}/document/{stable}/{}", m.revision());
        ExportIds {
            namespace: ns,
            component: root.clone(),
            sequence: format!("{root}/sequence"),
        }
    });
    let namespace = RwSignal::new(default.namespace);
    let export_component = RwSignal::new(default.component);
    let export_sequence = RwSignal::new(default.sequence);
    let allow_export_loss = RwSignal::new(false);
    let exported = RwSignal::new(None::<Arc<ExportReview>>);
    let get_ids = move || ExportIds {
        namespace: namespace.get_untracked(),
        component: export_component.get_untracked(),
        sequence: export_sequence.get_untracked(),
    };
    Effect::new(move |_| {
        let _ = namespace.get();
        let _ = export_component.get();
        let _ = export_sequence.get();
        let _ = allow_export_loss.get();
        exported.set(None);
    });
    let review_export = move |_| {
        let result = (|| {
            checked.run(())?;
            let ids = get_ids();
            let p = policy(allow_export_loss.get_untracked());
            let value = current.with_value(|m| {
                collagene_core::formats::fidelity::reference_export_guard(m)
                    .map_err(|e| e.to_string())?;
                exchange::export_molecule(m, &ids, p).map_err(|e| e.to_string())
            })?;
            checked.run(())?;
            Ok::<_, String>(Arc::new(ExportReview {
                ids,
                policy: p,
                value: Arc::new(value),
            }))
        })();
        match result {
            Ok(value) => {
                exported.set(Some(value));
                status.set("Export prepared. Review the declared identities and any omitted native fields before saving.".into());
            }
            Err(e) => status.set(e),
        }
    };
    let save_export = move |_| {
        let Some(saved) = exported.get_untracked() else {
            return;
        };
        let saved_check = saved.clone();
        let check = move || {
            checked.run(())?;
            current.with_value(|m| {
                collagene_core::formats::fidelity::reference_export_guard(m)
                    .map_err(|e| e.to_string())
            })?;
            if get_ids() != saved_check.ids
                || policy(allow_export_loss.get_untracked()) != saved_check.policy
            {
                return Err(
                    "SBOL export identities or loss policy changed. Prepare it again.".into(),
                );
            }
            Ok(())
        };
        let result = check().and_then(|()| {
            prepare_download("current-design.sbol.jsonld", MEDIA_TYPE, &saved.value.bytes)
        });
        complete(
            result,
            check,
            status,
            "Saved the reviewed SBOL component/sequence/feature graph. This is design exchange, not transferred physical or verification authority.",
        );
    };
    view! {<section class="sbol-exchange"><h2>"SBOL design exchange"</h2><p class="helper">"SBOL 3.1 structural JSON-LD profile. Design relationships and original evidence remain separate from native sequence references, cloning chemistry and current verification."</p>
     <p role="status" aria-live="polite" style="overflow-wrap:anywhere">{move||status.get()}</p>
     <Show when=move||source_present><section><h3>"Loaded design graph"</h3><button class="text-tool" on:click=move|_|graph_download.run(false)>"Save exact original file"</button>{graph.get_value().map(|doc|view!{<><button class="text-tool" on:click=move|_|graph_download.run(true)>"Save retained expanded graph"</button><GraphReview doc/><IriPicker label="Component" choices=components selected=component/><IriPicker label="Sequence" choices=sequences selected=sequence/><label>"Native topology choice"<select prop:value=move||topology.get() on:change=move|e|topology.set(event_target_value(&e))><option value="">"Use declared topology / applicable default"</option><option value="linear">"Explicit linear reference"</option><option value="circular">"Explicit circular reference"</option></select></label><p class="helper">"Unspecified DNA topology needs an explicit choice; a choice cannot override conflicting source topology. Child sequences are never assembled automatically."</p><label class="checkbox-label"><input type="checkbox" prop:checked=move||allow_projection_loss.get() on:change=move|e|allow_projection_loss.set(event_target_checked(&e))/>"Allow explicit representation losses for review"</label><button class="button secondary" on:click=review>"Review selected sequence projection"</button>{move||preview.get().map(|saved|{let value=&saved.value;let issues=Arc::new(value.issues.clone());view!{<section><h4>"Selected native reference"</h4><p style="overflow-wrap:anywhere">{format!("Component {} · Sequence {}",value.component,value.sequence)}</p><p>{format!("{} · {:?} · {} symbols · {:?} · {} mapped features",value.molecule.name,value.molecule.sequence.alphabet(),value.molecule.sequence.len(),value.molecule.topology,value.molecule.features.len())}</p>{value.warnings.iter().map(|w|view!{<p class="helper">{w.clone()}</p>}).collect_view()}<Issues issues/><p class="helper">"The full original graph will be retained locally with this sequence reference. Publishing a selected database item requires a separate scoped graph-closure check; unrelated source objects must not be silently shared."</p><button class="button primary" on:click=use_projection>"Open reviewed sequence reference"</button></section>}})}</>})}</section></Show>
     <details open=!source_present><summary>"Export current native sequence as SBOL"</summary><p>{current.with_value(|m|format!("Current reference: {} · {:?} · {} symbols",m.name,m.sequence.alphabet(),m.sequence.len()))}</p><p class="helper">"Choose stable identities under a namespace you control. Defaults are deterministic placeholders based on the native document identity and exact revision; they do not register a public identifier. Physical/reference-only products that this profile cannot preserve require a native archive instead."</p><label>"Namespace URL"<input type="text" prop:value=move||namespace.get() on:input=move|e|namespace.set(event_target_value(&e))/></label><label>"Component IRI"<input type="text" prop:value=move||export_component.get() on:input=move|e|export_component.set(event_target_value(&e))/></label><label>"Sequence IRI"<input type="text" prop:value=move||export_sequence.get() on:input=move|e|export_sequence.set(event_target_value(&e))/></label><label class="checkbox-label"><input type="checkbox" prop:checked=move||allow_export_loss.get() on:change=move|e|allow_export_loss.set(event_target_checked(&e))/>"Allow explicit native-field omissions for review"</label><button class="button secondary" on:click=review_export>"Prepare SBOL export"</button>{move||exported.get().map(|saved|{let issues=Arc::new(saved.value.issues.clone());view!{<section><p>{format!("{} reviewed JSON-LD bytes",saved.value.bytes.len())}</p><Issues issues/><button class="button primary" on:click=save_export>"Save reviewed SBOL graph"</button></section>}})}</details><button class="button secondary" on:click=move|_|on_close.run(())>"Close"</button>
    </section>}
}
