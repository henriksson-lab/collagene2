//! A named digest belongs to one exact retained source event. Saved material
//! dependencies follow that event's existing native run/derivation owner.
use super::{Error, Graph};
use collagene_core::{cloning::restriction::setups, documents::workspace::Workspace};
use std::collections::BTreeMap;

pub(super) struct Index<'a> {
    source: &'a Workspace,
    by_event: BTreeMap<String, Vec<String>>,
}

impl<'a> Index<'a> {
    pub(super) fn new(source: &'a Workspace) -> Result<Self, Error> {
        let mut by_event = BTreeMap::<String, Vec<String>>::new();
        // Validate/list once, rather than validating the whole workspace anew
        // for each setup. Listing attests exact Run/Derived ancestry without
        // replaying material or making historical evidence current.
        for record in setups::list(source)? {
            by_event
                .entry(record.source.event_id)
                .or_default()
                .push(record.attachment_hash);
        }
        Ok(Self { source, by_event })
    }

    pub(super) fn follow_event(&self, graph: &mut Graph<'a>, event: &str) -> Result<(), Error> {
        if !std::ptr::eq(self.source, graph.source) {
            return Err(Error::Unresolved(
                "digest setup index belongs to another workspace",
            ));
        }
        if let Some(hashes) = self.by_event.get(event) {
            for hash in hashes {
                graph.edge()?;
                graph.kept_assets.insert(hash.clone());
            }
        }
        Ok(())
    }
}
