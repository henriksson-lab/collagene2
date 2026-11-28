//! Working tab identities and their explorer projections; retained histories stay in WorkspaceSession.
use crate::workspace::session::WorkspaceSession;
use leptos::prelude::*;
#[derive(Clone, Copy)]
pub(crate) struct DocumentTabs {
    session: WorkspaceSession,
    pub(crate) open_tabs: RwSignal<crate::workspace::document_tabs_model::Tabs>,
    pub(crate) has_open_document: Memo<bool>,
}
pub(crate) struct DocumentItems {
    pub(crate) document_items: Memo<Vec<crate::workspace::document_explorer::DocumentItem>>,
    pub(crate) explorer_active: Memo<String>,
    pub(crate) explorer_documents: Memo<Vec<crate::workspace::document_explorer::DocumentItem>>,
    pub(crate) tab_items: Memo<Vec<crate::workspace::document_explorer::DocumentItem>>,
    pub(crate) explorer_locked: Memo<bool>,
    pub(crate) explorer_activate: Callback<String>,
}
impl DocumentTabs {
    pub(crate) fn new(session: WorkspaceSession) -> Self {
        let WorkspaceSession {
            documents,
            active,
            workspace_epoch,
            ..
        } = session;
        let open_tabs = RwSignal::new(crate::workspace::document_tabs_model::Tabs::default());
        let tab_epoch = StoredValue::new(None::<u64>);
        let tab_document = StoredValue::new(None::<String>);
        let retained_document_ids = Memo::new(move |_| {
            documents.with(|docs| docs.iter().map(|doc| doc.id.clone()).collect::<Vec<_>>())
        });
        let has_open_document = Memo::new(move |_| open_tabs.with(|tabs| tabs.active().is_some()));
        Effect::new(move |_| {
            let epoch = workspace_epoch.get();
            let index = active.get();
            let id = documents.with_untracked(|docs| docs.get(index).map(|doc| doc.id.clone()));
            let current = open_tabs.get_untracked();
            let mut next = current.clone();
            if tab_epoch.get_value() != Some(epoch) {
                // A finalized assembly also replaces the workspace frame. Keep
                // working source tabs whose IDs remain in that frame; a different
                // workspace naturally removes IDs it does not retain.
                let retained = documents.with_untracked(|docs| {
                    docs.iter().map(|doc| doc.id.clone()).collect::<Vec<_>>()
                });
                next.reconcile(&retained);
                if next.ids().is_empty() {
                    next.reset(id.as_deref());
                } else if let Some(id) = &id {
                    next.open(id);
                }
            } else if let Some(id) = &id {
                next.open(id);
            }
            tab_epoch.set_value(Some(epoch));
            tab_document.set_value(id);
            if next != current {
                open_tabs.set(next);
            }
        });
        Effect::new(move |_| {
            let retained = retained_document_ids.get();
            let id = documents
                .with_untracked(|docs| docs.get(active.get_untracked()).map(|doc| doc.id.clone()));
            let current = open_tabs.get_untracked();
            let mut next = current.clone();
            next.reconcile(&retained);
            // A replacement can change the active document ID without changing its
            // array index. Ordinary revision edits do not reopen closed tabs.
            if tab_document.get_value() != id {
                if let Some(id) = &id {
                    next.open(id);
                }
                tab_document.set_value(id);
            }
            if next != current {
                open_tabs.set(next);
            }
        });
        Self {
            session,
            open_tabs,
            has_open_document,
        }
    }
    pub(crate) fn items(
        &self,
        show_prepared_history: RwSignal<bool>,
        hidden_products: RwSignal<std::collections::BTreeSet<String>>,
    ) -> DocumentItems {
        let WorkspaceSession {
            documents,
            active,
            workspace_restoring,
            ..
        } = self.session;
        let open_tabs = self.open_tabs;
        let document_items = Memo::new(move |_| {
            documents.with(|docs| {
                docs.iter()
                    .map(|doc| {
                        let provider = doc.molecule.metadata.get("source:provider");
                        let accession = doc
                            .molecule
                            .metadata
                            .get("source:accession_version")
                            .or_else(|| doc.molecule.metadata.get("source:identifier"));
                        crate::workspace::document_explorer::DocumentItem {
                            id: doc.id.clone(),
                            name: doc.molecule.name.clone(),
                            length: doc.molecule.sequence.len(),
                            alphabet: doc.molecule.sequence.alphabet(),
                            topology: doc.molecule.topology,
                            source: if doc
                                .molecule
                                .metadata
                                .get("assembly:intent:reference_only")
                                .or_else(|| {
                                    doc.molecule.metadata.get("mutation:intent:reference_only")
                                })
                                .or_else(|| doc.molecule.metadata.get("read_review:reference_only"))
                                .is_some_and(|value| value == "true")
                            {
                                Some("Intended reference".into())
                            } else {
                                provider.map(|provider| {
                                    accession
                                        .map(|id| format!("{provider} · {id}"))
                                        .unwrap_or_else(|| provider.clone())
                                })
                            },
                        }
                    })
                    .collect::<Vec<_>>()
            })
        });
        let explorer_active =
            Memo::new(move |_| open_tabs.with(|tabs| tabs.active().unwrap_or_default().to_owned()));
        let explorer_documents = Memo::new(move |_| {
            let show_all = show_prepared_history.get();
            let hidden = hidden_products.get();
            let current = explorer_active.get();
            document_items.with(|items| {
                items
                    .iter()
                    .filter(|doc| show_all || !hidden.contains(&doc.id) || doc.id == current)
                    .cloned()
                    .collect::<Vec<_>>()
            })
        });
        let tab_items = Memo::new(move |_| {
            let index = document_items.with(|items| {
                items
                    .iter()
                    .map(|item| (item.id.clone(), item.clone()))
                    .collect::<std::collections::BTreeMap<_, _>>()
            });
            open_tabs.with(|tabs| {
                tabs.ids()
                    .iter()
                    .filter_map(|id| index.get(id).cloned())
                    .collect::<Vec<_>>()
            })
        });
        let explorer_locked = Memo::new(move |_| workspace_restoring.get());
        let explorer_activate = Callback::new(move |id: String| {
            if workspace_restoring.get_untracked() {
                return;
            }
            if let Some(index) =
                documents.with_untracked(|docs| docs.iter().position(|doc| doc.id == id))
            {
                batch(|| {
                    open_tabs.update(|tabs| tabs.open(&id));
                    active.set(index);
                });
            }
        });
        DocumentItems {
            document_items,
            explorer_active,
            explorer_documents,
            tab_items,
            explorer_locked,
            explorer_activate,
        }
    }
}
