//! The shared reaction editor writes directly into the selected preparation row.
use super::cloning_methods_model::MethodsState;
use crate::{
    restriction::{reaction_conditions::ReactionConditionsPanel, reaction_conditions_model::Draft},
    sources::identity::SourceInfo,
};
use collagene_core::cloning::restriction::catalog::{self, CatalogSelection};
use leptos::prelude::*;
use std::sync::Arc;

#[component]
pub(crate) fn PreparationReactionConditions(
    config: RwSignal<MethodsState>,
    index: usize,
    source_info: Callback<usize, Option<SourceInfo>>,
    epoch: RwSignal<u64>,
) -> impl IntoView {
    let catalog = StoredValue::new(
        catalog::builtin_catalog()
            .map(Arc::new)
            .map_err(|e| e.to_string()),
    );
    let selection = Memo::new(move |_| -> Result<CatalogSelection, String> {
        config.with(|state| {
            state
                .rows
                .get(index)
                .and_then(|row| row.catalog_selection.clone())
                .ok_or_else(|| {
                    "Pin the exact preparation enzymes before declaring reaction conditions.".into()
                })
        })
    });
    let source = Memo::new(move |_| {
        config.with(|state| {
            state
                .rows
                .get(index)
                .and_then(|row| source_info.run(row.document))
        })
    });
    let origin = Memo::new(move |_| {
        config.with(|state| state.rows.get(index).and_then(|row| row.run_origin.clone()))
    });
    let draft = Memo::new(move |_| {
        config.with(|state| {
            state
                .rows
                .get(index)
                .map(|row| row.reaction_conditions.clone())
                .unwrap_or_default()
        })
    });
    let on_change = Callback::new(move |value: Draft| {
        config.update(|state| {
            if let Some(row) = state.rows.get_mut(index) {
                row.reaction_conditions = value;
            }
        });
    });
    let conditions = RwSignal::new(None);
    let error = RwSignal::new(None);
    view! {
        <ReactionConditionsPanel
            catalog selection source epoch origin draft=Signal::from(draft) on_change conditions error
        />
    }
}
