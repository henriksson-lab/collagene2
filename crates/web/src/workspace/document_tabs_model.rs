//! Presentation-only open document tabs, separate from retained documents/history.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Tabs {
    ids: Vec<String>,
    active: Option<String>,
}
impl Tabs {
    pub fn ids(&self) -> &[String] {
        &self.ids
    }
    pub fn active(&self) -> Option<&str> {
        self.active.as_deref()
    }
    pub fn open(&mut self, id: &str) {
        if id.is_empty() {
            return;
        }
        if !self.ids.iter().any(|open| open == id) {
            self.ids.push(id.into());
        }
        self.active = Some(id.into());
    }
    pub fn close(&mut self, id: &str) {
        let Some(index) = self.ids.iter().position(|open| open == id) else {
            return;
        };
        let selected = self.active.as_deref() == Some(id);
        self.ids.remove(index);
        if selected {
            self.active = self.ids.get(index).or_else(|| self.ids.last()).cloned();
        }
    }
    pub fn reconcile(&mut self, retained: &[String]) {
        self.ids.retain(|id| retained.contains(id));
        if self
            .active
            .as_ref()
            .is_some_and(|id| !self.ids.contains(id))
        {
            self.active = self.ids.first().cloned();
        }
    }
    pub fn reset(&mut self, id: Option<&str>) {
        self.ids.clear();
        self.active = None;
        if let Some(id) = id {
            self.open(id);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn opening_is_idempotent_preserves_order_and_selects_existing_tabs() {
        let mut tabs = Tabs::default();
        tabs.open("");
        assert!(tabs.ids().is_empty());
        assert_eq!(tabs.active(), None);
        for id in ["backbone", "insert", "target", "backbone"] {
            tabs.open(id);
        }
        assert_eq!(tabs.ids(), &["backbone", "insert", "target"]);
        assert_eq!(tabs.active(), Some("backbone"));
    }
    #[test]
    fn closing_background_keeps_selection_and_active_prefers_right_then_left() {
        let mut tabs = Tabs::default();
        for id in ["a", "b", "c", "d"] {
            tabs.open(id);
        }
        tabs.open("b");
        tabs.close("a");
        assert_eq!(tabs.active(), Some("b"));
        tabs.close("absent");
        assert_eq!(tabs.ids(), &["b", "c", "d"]);
        tabs.close("b");
        assert_eq!(tabs.active(), Some("c"));
        tabs.open("d");
        tabs.close("d");
        assert_eq!(tabs.active(), Some("c"));
        tabs.close("c");
        assert!(tabs.ids().is_empty());
        assert_eq!(tabs.active(), None);
        tabs.open("b");
        assert_eq!(tabs.ids(), &["b"]);
        assert_eq!(tabs.active(), Some("b"));
    }
    #[test]
    fn reconciliation_keeps_tab_order_and_closed_state_without_mutating_retained_documents() {
        let retained = vec!["c".into(), "a".into(), "new".into()];
        let before = retained.clone();
        let mut tabs = Tabs::default();
        for id in ["a", "b", "c"] {
            tabs.open(id);
        }
        tabs.open("b");
        tabs.reconcile(&retained);
        assert_eq!(tabs.ids(), &["a", "c"]);
        assert_eq!(tabs.active(), Some("a"));
        assert_eq!(retained, before);
        tabs.open("c");
        tabs.reconcile(&retained);
        assert_eq!(tabs.active(), Some("c"));
        tabs.close("a");
        tabs.close("c");
        tabs.reconcile(&retained);
        assert!(tabs.ids().is_empty());
        assert_eq!(tabs.active(), None);
        assert_eq!(retained, before);
    }
    #[test]
    fn reset_discards_only_old_presentation_and_opens_the_requested_workspace_tab() {
        let mut tabs = Tabs::default();
        tabs.open("old-backbone");
        tabs.open("old-target");
        tabs.reset(Some("new-workspace-first"));
        assert_eq!(tabs.ids(), &["new-workspace-first"]);
        assert_eq!(tabs.active(), Some("new-workspace-first"));
        tabs.reset(None);
        assert_eq!(tabs, Tabs::default());
        tabs.reset(Some(""));
        assert_eq!(tabs, Tabs::default());
    }
}
