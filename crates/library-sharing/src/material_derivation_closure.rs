//! Exact target-free material-record dependencies. This child module is wired
//! into Graph only when the matching native schema/validator is registered.
//! Indexing an unrelated record never retains its bytes or source documents.
use super::*;
use collagene_core::material_derivation::{self as native, Capsule};

struct Owner {
    attachment: String,
    product: EventBinding,
}
/// Immutable ownership index; Graph remains the single bounded retention queue.
/// The source borrow prevents an index outliving its validated native archive.
pub(super) struct Index<'a> {
    source: &'a Workspace,
    owners: BTreeMap<String, Owner>,
}
impl<'a> Index<'a> {
    /// Native documentary validation is historical-engine tolerant. This does
    /// not replay chemistry, grant current material, or admit arbitrary schemas.
    pub(super) fn new(source: &'a Workspace) -> Result<Self, Error> {
        native::validate_workspace_records(source)?;
        let mut owners = BTreeMap::new();
        let events: BTreeMap<_, _> = source
            .history
            .iter()
            .map(|event| (event.id.as_str(), event))
            .collect();
        for (hash, attachment) in &source.attachments {
            if attachment.kind != native::WORKSPACE_ATTACHMENT_KIND {
                continue;
            }
            if owners.len() >= native::MAX_RECORDS
                || attachment.bytes.len() > native::MAX_RECORD_BYTES
                || attachment.media_type != native::MEDIA_TYPE
                || blake3::hash(&attachment.bytes).to_hex().as_str() != hash
            {
                return Err(Error::Unresolved("material derivation attachment binding"));
            }
            let record = native::load_record(&attachment.bytes)?;
            let event = events
                .get(record.product.event_id.as_str())
                .ok_or(Error::Unresolved("material derivation product event"))?;
            // No producer-label substring or sequence matching establishes an
            // owner. Both the exact native label and complete pins are required.
            if event.operation != native::EVENT_OPERATION
                || event.document_id != record.product.document_id
                || event.revision != record.product.revision
                || event.parents.len() != 1
                || event.parents.first() != Some(&record.request.source.event_id)
                || source.revisions.get(&record.product.revision).is_none()
                || owners
                    .insert(
                        record.product.event_id.clone(),
                        Owner {
                            attachment: hash.clone(),
                            product: record.product,
                        },
                    )
                    .is_some()
            {
                return Err(Error::Unresolved("material derivation exact product owner"));
            }
        }
        Ok(Self { source, owners })
    }
    /// Invoke after ordinary run and saved-digest indexes have been built. A
    /// product cannot acquire a second native owner depending on attachment order.
    pub(super) fn validate_existing_owners(
        &self,
        run_owners: &BTreeMap<String, String>,
        digest_owners: &BTreeMap<String, String>,
    ) -> Result<(), Error> {
        if self
            .owners
            .keys()
            .any(|id| run_owners.contains_key(id) || digest_owners.contains_key(id))
        {
            return Err(Error::Unresolved(
                "multiple native owners for material product event",
            ));
        }
        Ok(())
    }
    pub(super) fn contains_product(&self, event: &str) -> bool {
        self.owners.contains_key(event)
    }
    /// Follow only a retained exact derived product. All dependencies enter the
    /// ordinary Graph queues, sharing their global edge limits/deduplication.
    /// Derived dependencies are not recursively decoded by this sharing module.
    pub(super) fn follow_product(&self, graph: &mut Graph<'a>, event: &str) -> Result<bool, Error> {
        let Some(owner) = self.owners.get(event) else {
            return Ok(false);
        };
        self.same_source(graph)?;
        let attachment = self
            .source
            .attachments
            .get(&owner.attachment)
            .ok_or(Error::Unresolved("material derivation record bytes"))?;
        let record = native::load_record(&attachment.bytes)?;
        if record.product != owner.product {
            return Err(Error::Unresolved(
                "material derivation indexed product changed",
            ));
        }
        graph.edge()?;
        graph.kept_assets.insert(owner.attachment.clone());
        graph.binding(&record.product)?;
        for binding in native::source_bindings(&record)? {
            graph.binding(&binding)?;
        }
        for run in native::required_run_ids(&record)? {
            graph.run(&run)?;
        }
        for binding in native::required_derivation_products(&record)? {
            let expected = self
                .owners
                .get(&binding.event_id)
                .ok_or(Error::Unresolved("nested material derivation owner"))?;
            if expected.product != binding {
                return Err(Error::Unresolved("nested material derivation exact event"));
            }
            graph.binding(&binding)?;
        }
        Ok(true)
    }
    /// A reconstruction operation's documentary capsule must equal the retained
    /// canonical record. It cannot introduce a new private dependency contract.
    pub(super) fn follow_capsule(
        &self,
        graph: &mut Graph<'a>,
        capsule: &Capsule,
    ) -> Result<(), Error> {
        self.same_source(graph)?;
        native::validate_capsule_for_workspace(self.source, capsule, &capsule.product)?;
        let owner = self
            .owners
            .get(&capsule.product.event_id)
            .ok_or(Error::Unresolved("material capsule has no retained owner"))?;
        let asset = self
            .source
            .attachments
            .get(&owner.attachment)
            .ok_or(Error::Unresolved("material capsule record attachment"))?;
        if owner.product != capsule.product || asset.bytes != capsule.record_bytes {
            return Err(Error::Unresolved(
                "material capsule differs from retained native record",
            ));
        }
        graph.binding(&capsule.product)?;
        self.follow_product(graph, &capsule.product.event_id)?;
        Ok(())
    }
    fn same_source(&self, graph: &Graph<'a>) -> Result<(), Error> {
        if !std::ptr::eq(self.source, graph.source) {
            return Err(Error::Unresolved(
                "material ownership index belongs to another workspace",
            ));
        }
        Ok(())
    }
    /// Call for every retained event after the known-owner branches. This exact
    /// reserved producer label without its typed record is not a plain reference.
    pub(super) fn reject_orphan(&self, event: &workspace::RevisionEvent) -> Result<(), Error> {
        if event.operation == native::EVENT_OPERATION && !self.contains_product(&event.id) {
            return Err(Error::Unresolved(
                "orphan material derivation product event",
            ));
        }
        Ok(())
    }
}
