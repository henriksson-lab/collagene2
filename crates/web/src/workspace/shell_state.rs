//! Shared presentation notices and panel/modal routing state.
use leptos::prelude::*;
#[derive(Clone, Copy)]
pub(crate) struct ShellState {
    pub(crate) panel: RwSignal<&'static str>,
    pub(crate) message: RwSignal<String>,
    pub(crate) modal: RwSignal<&'static str>,
    pub(crate) modal_error: RwSignal<String>,
}
impl ShellState {
    pub(crate) fn workspace_replaced(&self) {
        self.modal.set("");
    }
    pub(crate) fn new() -> Self {
        Self {
            panel: RwSignal::new(""),
            message: RwSignal::new(String::new()),
            modal: RwSignal::new(""),
            modal_error: RwSignal::new(String::new()),
        }
    }
}
