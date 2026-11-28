//! Document and feature editors commit against an explicit base revision.
use collagene_core::documents::reference_edit::Operation;
use collagene_model::{Alphabet, Error, Feature, Location, Molecule, Topology};
use leptos::prelude::*;
type Commit = Callback<(String, Vec<Operation>), Result<(), String>>;

#[derive(Clone)]
struct QualifierRow {
    id: usize,
    key: ArcRwSignal<String>,
    value: ArcRwSignal<String>,
    flag: ArcRwSignal<bool>,
}
impl QualifierRow {
    fn new(id: usize, key: String, value: Option<String>) -> Self {
        Self {
            id,
            key: ArcRwSignal::new(key),
            flag: ArcRwSignal::new(value.is_none()),
            value: ArcRwSignal::new(value.unwrap_or_default()),
        }
    }
}
#[component]
fn QualifierFields(row: QualifierRow, on_remove: Callback<usize>) -> impl IntoView {
    let key_value = row.key.clone();
    let key_input = row.key.clone();
    let value = row.value.clone();
    let input = row.value.clone();
    let flag = row.flag.clone();
    let flag_input = row.flag.clone();
    let show_value = row.flag.clone();
    let id = row.id;
    view! {<div class="qualifier-row"><label>"Qualifier"<input aria-label="Qualifier name" prop:value=move||key_value.get() on:input=move|e|key_input.set(event_target_value(&e))/></label><label class="checkbox-label"><input type="checkbox" aria-label="Qualifier has no value" prop:checked=move||flag.get() on:change=move|e|flag_input.set(event_target::<web_sys::HtmlInputElement>(&e).checked())/>"Flag (no value)"</label><Show when=move||!show_value.get()><label>"Value"<textarea aria-label="Qualifier value" rows="2" prop:value={let value=value.clone();move||value.get()} on:input={let input=input.clone();move|e|input.set(event_target_value(&e))}/></label></Show><button class="text-tool" on:click=move|_|on_remove.run(id)>"Remove qualifier"</button></div>}
}
#[component]
pub fn FeatureEditor(
    molecule: Molecule,
    index: Option<usize>,
    selection: Option<(usize, usize)>,
    #[prop(optional_no_strip)] initial_location: Option<Location>,
    on_commit: Commit,
) -> impl IntoView {
    let original = index.and_then(|i| molecule.features.get(i).cloned());
    let old_label = original
        .as_ref()
        .map(|f| f.label().to_string())
        .unwrap_or_default();
    let label = RwSignal::new(old_label.clone());
    let kind = RwSignal::new(
        original
            .as_ref()
            .map(|f| f.kind.clone())
            .unwrap_or_else(|| "misc_feature".into()),
    );
    let initial_location = original
        .as_ref()
        .map(|f| f.location.clone())
        .or(initial_location)
        .unwrap_or_else(|| {
            let (a, b) = selection.unwrap_or((0, molecule.sequence.len()));
            Location::range(a, b)
        });
    let location = RwSignal::new(
        collagene_model::annotations::insdc::location_text(&initial_location).unwrap_or_default(),
    );
    let review = RwSignal::new(original.as_ref().is_some_and(|f| f.needs_review));
    let error = RwSignal::new(String::new());
    let base_revision = molecule.revision();
    let source = StoredValue::new(molecule);
    let label_index = original.as_ref().and_then(|f| {
        f.qualifiers
            .iter()
            .position(|(key, value)| key == "label" && value.is_some())
    });
    let initial: Vec<_> = original
        .as_ref()
        .map(|f| f.qualifiers.clone())
        .unwrap_or_default()
        .into_iter()
        .enumerate()
        .filter(|(i, _)| Some(*i) != label_index)
        .map(|(id, (key, value))| QualifierRow::new(id, key, value))
        .collect();
    let next = RwSignal::new(original.as_ref().map_or(0, |f| f.qualifiers.len()));
    let rows = RwSignal::new(initial);
    let save = move |_| {
        let source = source.get_value();
        let text = location.get_untracked();
        let parsed = if collagene_model::annotations::insdc::location_text(&initial_location)
            .ok()
            .as_deref()
            == Some(text.trim())
        {
            Ok(initial_location.clone())
        } else {
            collagene_model::annotations::insdc::parse_location(&text, source.sequence.len())
        };
        let parsed = match parsed {
            Ok(location) => location,
            Err(e) => {
                error.set(e.to_string());
                return;
            }
        };
        let mut qualifiers: Vec<_> = rows
            .get_untracked()
            .into_iter()
            .map(|row| {
                (
                    row.key.get_untracked(),
                    if row.flag.get_untracked() {
                        None
                    } else {
                        Some(row.value.get_untracked())
                    },
                )
            })
            .collect();
        let new_label = label.get_untracked();
        if let Some(i) = label_index {
            qualifiers.insert(i.min(qualifiers.len()), ("label".into(), Some(new_label)));
        } else if !new_label.trim().is_empty() && (index.is_none() || new_label != old_label) {
            qualifiers.insert(0, ("label".into(), Some(new_label)));
        }
        let mut feature = Feature {
            kind: kind.get_untracked(),
            location: parsed,
            qualifiers,
            needs_review: review.get_untracked(),
        };
        if original
            .as_ref()
            .is_some_and(|f| f.location != feature.location)
        {
            feature.qualifiers.retain(|(key, _)| key != "translation");
            feature.needs_review = true;
        }
        if let Err(e) = on_commit.run((
            base_revision.clone(),
            vec![Operation::SetFeature { index, feature }],
        )) {
            error.set(e);
        }
    };
    view! {<div class="feature-editor"><p>"Locations use one-based inclusive coordinates. Compound examples: join(1..4,13..16), complement(join(13..16,1..4)). Fuzzy and unresolved locations stay distinct from exact sequence."</p><Show when=move||!error.get().is_empty()><p class="dialog-error" role="alert">{move||error.get()}</p></Show><label>"Display name"<input aria-label="Annotation name" prop:value=move||label.get() on:input=move|e|label.set(event_target_value(&e))/></label><label>"Feature type"<input aria-label="Annotation type" prop:value=move||kind.get() on:input=move|e|kind.set(event_target_value(&e))/></label><label>"Location"<input aria-label="Annotation location" prop:value=move||location.get() on:input=move|e|location.set(event_target_value(&e))/></label><label class="checkbox-label"><input type="checkbox" prop:checked=move||review.get() on:change=move|e|review.set(event_target::<web_sys::HtmlInputElement>(&e).checked())/>"Review required"</label><h3>"Qualifiers"</h3><For each=move||rows.get() key=|row|row.id children=move|row|view!{<QualifierFields row=row on_remove=Callback::new(move|id|rows.update(|rows|rows.retain(|row|row.id!=id)))/>}/><button class="button secondary" on:click=move|_|{let id=next.get_untracked();next.set(id+1);rows.update(|rows|rows.push(QualifierRow::new(id,"note".into(),Some(String::new()))));}>"Add qualifier"</button><button class="button primary" on:click=save>"Save annotation"</button></div>}
}

#[component]
pub fn DocumentEditor(
    molecule: Molecule,
    selection: Option<(usize, usize)>,
    on_commit: Commit,
) -> impl IntoView {
    let name = RwSignal::new(molecule.name.clone());
    let topology = RwSignal::new(if molecule.topology == Topology::Circular {
        "circular"
    } else {
        "linear"
    });
    let origin = RwSignal::new(selection.map_or(1, |(a, _)| a + 1).to_string());
    let error = RwSignal::new(String::new());
    let base = molecule.revision();
    let source = StoredValue::new(molecule);
    let commit = move |result: Result<Vec<Operation>, Error>| match result {
        Ok(value) => {
            if let Err(e) = on_commit.run((base.clone(), value)) {
                error.set(e);
            }
        }
        Err(e) => error.set(e.to_string()),
    };
    let commit_settings = commit.clone();
    let commit_reverse = commit.clone();
    let commit_origin = commit.clone();
    view! {<div class="document-editor"><p>"These controls edit the reference document. Physical cleavage, end chemistry and ligation belong in a cloning plan. Changes can be undone."</p><Show when=move||!error.get().is_empty()><p class="dialog-error" role="alert">{move||error.get()}</p></Show><label>"Document name"<input aria-label="Document name" prop:value=move||name.get() on:input=move|e|name.set(event_target_value(&e))/></label><label>"Reference topology"<select aria-label="Reference topology" prop:value=move||topology.get() on:change=move|e|topology.set(if event_target_value(&e)=="circular"{"circular"}else{"linear"})><option value="linear">"Linear"</option><option value="circular">"Circular"</option></select></label><button class="button primary" on:click=move|_|{let result=Ok(vec![Operation::Rename {name:name.get_untracked()},Operation::SetTopology {topology:if topology.get_untracked()=="circular"{Topology::Circular}else{Topology::Linear}}]);commit_settings(result);}>"Apply document settings"</button><hr/><button class="button secondary" disabled=move||source.with_value(|m|m.sequence.alphabet()==Alphabet::Protein) on:click=move|_|commit_reverse(Ok(vec![Operation::ReverseComplement]))>"Reverse complement document"</button><label>"New origin (one-based base number)"<input aria-label="New origin" type="number" min="1" prop:value=move||origin.get() on:input=move|e|origin.set(event_target_value(&e))/></label><button class="button secondary" disabled=move||source.with_value(|m|m.topology!=Topology::Circular) on:click=move|_|{let result=origin.get_untracked().parse::<usize>().ok().and_then(|n|n.checked_sub(1)).ok_or_else(||Error::Coordinates("origin must be a one-based base number".into())).map(|n|vec![Operation::RotateOrigin {new_origin:n}]);commit_origin(result);}>"Rotate origin"</button></div>}
}
