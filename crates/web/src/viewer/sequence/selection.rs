use collagene_model::Molecule;
use leptos::prelude::*;

#[derive(Clone, Copy)]
pub(crate) struct ViewerSelection {
    pub(crate) selection: RwSignal<Option<(usize, usize)>>,
    pub(crate) selected_feature: RwSignal<Option<usize>>,
    pub(crate) search_regions: RwSignal<Option<(String, Vec<(usize, usize)>)>>,
    pub(crate) molecule: Memo<Molecule>,
    pub(crate) has_open_document: Memo<bool>,
}
