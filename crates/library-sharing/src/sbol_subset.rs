//! Explicitly reviewed documentary SBOL projection. Native identities never change.
use crate::{Error, Export, Inventory, Projection, Selection, Workspace, sbol_source};
use collagene_core::documents::workspace::{self, WorkspaceAttachment};
use collagene_design_exchange::{self as sbol, subset};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};
pub const SUBSET_KIND: &str = "sbol-selected-subset-v1";
pub const RECEIPT_KIND: &str = "sbol-selected-source-receipt-v1";
pub const MAX_RECEIPT_BYTES: usize = 1024 * 1024;
pub type PublicationChoices = BTreeMap<String, subset::SequenceChoices>;
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Receipt {
    pub schema_version: u32,
    pub source_profile: String,
    pub subset_policy: String,
    pub original_hash: String,
    pub subset_hash: String,
    pub request: subset::Request,
    pub omissions: subset::Omissions,
    pub original_bytes_included: bool,
    pub derivation: String,
}
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ReviewEntry {
    pub original_hash: String,
    pub subset_hash: String,
    pub request: subset::Request,
    pub omissions: subset::Omissions,
    pub retained_objects: usize,
}
/// Local review only. No omitted names, identities, or values are exposed here.
#[derive(Debug, Clone)]
pub struct PublicationReview {
    entries: Vec<ReviewEntry>,
    inventory: Inventory,
    binding: String,
    choices: PublicationChoices,
}
/// An opaque native/UI acknowledgement, never deserializable wire authority.
#[derive(Debug, Clone)]
pub struct ReviewedPublication {
    binding: String,
    choices: PublicationChoices,
}
impl PublicationReview {
    pub fn entries(&self) -> &[ReviewEntry] {
        &self.entries
    }
    pub fn inventory(&self) -> &Inventory {
        &self.inventory
    }
    pub fn approve(&self) -> ReviewedPublication {
        ReviewedPublication {
            binding: self.binding.clone(),
            choices: self.choices.clone(),
        }
    }
}
fn denied(reason: &'static str) -> Error {
    Error::SbolSource(reason)
}
fn digest(bytes: &[u8]) -> String {
    blake3::hash(bytes).to_hex().to_string()
}
fn hash(s: &str) -> bool {
    s.len() == 64
        && s.bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}
fn encode(receipt: &Receipt) -> Result<Vec<u8>, Error> {
    let bytes = serde_json::to_vec(receipt).map_err(|_| denied("receipt encoding failed"))?;
    if bytes.len() > MAX_RECEIPT_BYTES {
        return Err(Error::TooLarge);
    }
    Ok(bytes)
}
fn read_receipt(a: &WorkspaceAttachment) -> Result<Receipt, Error> {
    if a.kind != RECEIPT_KIND
        || a.media_type != "application/json"
        || a.bytes.len() > MAX_RECEIPT_BYTES
    {
        return Err(denied("invalid documentary receipt attachment"));
    }
    let r: Receipt =
        serde_json::from_slice(&a.bytes).map_err(|_| denied("invalid documentary receipt"))?;
    if r.schema_version != 1
        || r.source_profile != sbol::PROFILE
        || r.subset_policy != subset::POLICY
        || !hash(&r.original_hash)
        || !hash(&r.subset_hash)
        || r.original_hash == r.subset_hash
        || r.original_bytes_included
        || r.derivation != "unverified-documentary"
        || r.omissions.unrelated_objects > sbol::MAX_OBJECTS
        || r.omissions.alternative_sequence_edges > sbol::MAX_VALUES
        || encode(&r)? != a.bytes
    {
        return Err(denied("unsupported or noncanonical documentary receipt"));
    }
    Ok(r)
}
fn pairs(roots: &subset::SequenceChoices) -> BTreeSet<(String, String)> {
    roots
        .iter()
        .flat_map(|(c, ss)| ss.iter().map(move |s| (c.clone(), s.clone())))
        .collect()
}
/// An original-absent receipt checks scope and exact canonical bytes only. Its
/// claimed derivation and omission counts cannot be proved without the original.
pub(crate) fn received_assets(
    source: &Workspace,
    original: &str,
    claims: &sbol_source::Selections,
    limit: usize,
) -> Result<(BTreeSet<String>, usize), Error> {
    let required = pairs(claims);
    let mut covered = BTreeSet::new();
    let mut assets = BTreeSet::new();
    let mut used = 0usize;
    for (id, a) in &source.attachments {
        if a.kind != RECEIPT_KIND {
            continue;
        }
        let Ok(r) = read_receipt(a) else {
            continue;
        };
        if r.original_hash != original {
            continue;
        }
        let roots = pairs(&r.request.roots);
        if roots.is_empty() || !roots.is_subset(&required) {
            continue;
        }
        if digest(&a.bytes) != *id {
            return Err(denied("documentary receipt hash mismatch"));
        }
        let part = source
            .attachments
            .get(&r.subset_hash)
            .ok_or_else(|| denied("documentary subset is missing"))?;
        if part.kind != SUBSET_KIND
            || part.media_type != "application/ld+json"
            || digest(&part.bytes) != r.subset_hash
        {
            return Err(denied("documentary subset kind/hash mismatch"));
        }
        let doc = subset::validate_received(&part.bytes, &r.request)
            .map_err(|_| denied("documentary subset is not its exact supported minimal closure"))?;
        used = used
            .checked_add(doc.objects().len())
            .and_then(|n| n.checked_add(roots.len()))
            .filter(|n| *n <= limit)
            .ok_or(Error::TooLarge)?;
        if roots.iter().any(|p| covered.contains(p)) {
            return Err(denied(
                "overlapping documentary receipts ambiguously cover native claims",
            ));
        }
        covered.extend(roots);
        assets.insert(id.clone());
        assets.insert(r.subset_hash);
    }
    if covered != required {
        return Err(denied(
            "original SBOL source absent without complete documentary subset receipts",
        ));
    }
    Ok((assets, used))
}
fn make(
    source: &Workspace,
    selection: &Selection,
    choices: &PublicationChoices,
) -> Result<(PublicationReview, Projection), Error> {
    let native = crate::inspect_native_closure(source, selection)?;
    let (claims, _) = sbol_source::claims(source, &native.revisions, crate::MAX_REFERENCE_EDGES)?;
    if claims.is_empty() || choices.keys().any(|h| !claims.contains_key(h)) {
        return Err(denied(
            "subset review requires selected native source claims and no unrelated choices",
        ));
    }
    let mut projection = crate::project_inventory(source, selection, native)?;
    let mut entries = Vec::new();
    for (original, roots) in claims {
        let a=source.attachments.get(&original).ok_or_else(||denied("review requires exact original bytes; received documentary subsets can be republished through ordinary closure"))?;
        if a.kind != sbol_source::KIND
            || a.media_type != "application/ld+json"
            || digest(&a.bytes) != original
        {
            return Err(denied("original SBOL source kind/hash mismatch"));
        }
        let doc = sbol::import(&a.bytes).map_err(|_| denied("unsupported original SBOL source"))?;
        let request = subset::Request {
            roots,
            alternatives: choices.get(&original).cloned().unwrap_or_default(),
        };
        let built=subset::build(&doc,&request).map_err(|_|denied("unsupported selected structure or missing explicit nested alternative sequence choices"))?;
        if built.subset_hash == original {
            return Err(denied(
                "subset equals exact original; use ordinary exact-original sharing",
            ));
        }
        let receipt = Receipt {
            schema_version: 1,
            source_profile: sbol::PROFILE.into(),
            subset_policy: subset::POLICY.into(),
            original_hash: original.clone(),
            subset_hash: built.subset_hash.clone(),
            request: request.clone(),
            omissions: built.omissions.clone(),
            original_bytes_included: false,
            derivation: "unverified-documentary".into(),
        };
        let bytes = encode(&receipt)?;
        let receipt_hash = digest(&bytes);
        for (key, attachment) in [
            (
                built.subset_hash.clone(),
                WorkspaceAttachment {
                    kind: SUBSET_KIND.into(),
                    name: "selected-sbol-subset.jsonld".into(),
                    media_type: "application/ld+json".into(),
                    bytes: built.bytes,
                },
            ),
            (
                receipt_hash,
                WorkspaceAttachment {
                    kind: RECEIPT_KIND.into(),
                    name: "selected-sbol-source-receipt.json".into(),
                    media_type: "application/json".into(),
                    bytes,
                },
            ),
        ] {
            if projection
                .workspace
                .attachments
                .insert(key, attachment)
                .is_some()
            {
                return Err(denied(
                    "generated documentary attachment collides with native dependency",
                ));
            }
        }
        entries.push(ReviewEntry {
            original_hash: original,
            subset_hash: built.subset_hash,
            request,
            omissions: built.omissions,
            retained_objects: built.retained_objects.len(),
        });
    }
    projection.workspace.validate()?;
    let checked = crate::inspect_document_closure(&projection.workspace, selection)?;
    if checked.attachments != projection.workspace.attachments.keys().cloned().collect()
        || checked.revisions != projection.inventory.revisions
        || checked.events != projection.inventory.events
        || checked.document_heads != projection.inventory.document_heads
    {
        return Err(denied("reviewed closure scope changed native identities"));
    }
    projection.inventory.attachments = checked.attachments;
    projection.inventory.omitted.attachments = source
        .attachments
        .keys()
        .filter(|id| !projection.inventory.attachments.contains(*id))
        .count();
    // Bind exact selected native closure, selected original bytes via hashes,
    // declared choices, and omission counts. No unrelated private bytes are copied.
    let binding = digest(
        &serde_json::to_vec(&(
            selection,
            &projection.inventory,
            &choices,
            &entries
                .iter()
                .map(|e| (&e.original_hash, &e.subset_hash, &e.request, &e.omissions))
                .collect::<Vec<_>>(),
        ))
        .map_err(|_| denied("review binding encoding failed"))?,
    );
    let review = PublicationReview {
        entries,
        inventory: projection.inventory.clone(),
        binding,
        choices: choices.clone(),
    };
    Ok((review, projection))
}
pub fn prepare_documentary_subset(
    source: &Workspace,
    selection: &Selection,
    choices: &PublicationChoices,
) -> Result<PublicationReview, Error> {
    Ok(make(source, selection, choices)?.0)
}
pub fn project_document_closure_reviewed(
    source: &Workspace,
    selection: &Selection,
    approved: &ReviewedPublication,
) -> Result<Projection, Error> {
    let (review, projection) = make(source, selection, &approved.choices)?;
    if review.binding != approved.binding {
        return Err(denied(
            "publication review is stale; review the current exact closure again",
        ));
    }
    Ok(projection)
}
pub fn document_closure_reviewed(
    source: &Workspace,
    selection: &Selection,
    approved: &ReviewedPublication,
) -> Result<Export, Error> {
    let p = project_document_closure_reviewed(source, selection, approved)?;
    let archive = workspace::save_workspace(&p.workspace)?;
    Ok(Export {
        workspace: p.workspace,
        archive,
        inventory: p.inventory,
    })
}

/// Local-only choice UI data. Never copied into a documentary receipt unless a
/// sequence is explicitly chosen. Reachability follows structure, not all graphs.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ChoiceRequirement {
    pub original_hash: String,
    pub component_iri: String,
    pub available_sequences: BTreeSet<String>,
    pub required_sequences: BTreeSet<String>,
}
pub fn choice_requirements(
    source: &Workspace,
    selection: &Selection,
) -> Result<Vec<ChoiceRequirement>, Error> {
    let native = crate::inspect_native_closure(source, selection)?;
    let (claims, _) = sbol_source::claims(source, &native.revisions, crate::MAX_REFERENCE_EDGES)?;
    let mut result = Vec::new();
    for (hash, roots) in claims {
        let a = source
            .attachments
            .get(&hash)
            .ok_or_else(|| denied("choice review requires exact original source"))?;
        if a.kind != sbol_source::KIND
            || a.media_type != "application/ld+json"
            || digest(&a.bytes) != hash
        {
            return Err(denied("original SBOL source kind/hash mismatch"));
        }
        let doc = sbol::import(&a.bytes).map_err(|_| denied("unsupported original source"))?;
        let mut pending: Vec<_> = roots.keys().cloned().collect();
        let mut seen = BTreeSet::new();
        let mut edges = 0usize;
        while let Some(id) = pending.pop() {
            edges = edges
                .checked_add(1)
                .filter(|v| *v <= sbol::MAX_VALUES)
                .ok_or(Error::TooLarge)?;
            if !seen.insert(id.clone()) {
                continue;
            }
            let o = doc
                .objects()
                .get(&id)
                .ok_or_else(|| denied("unresolved selected component"))?;
            let sbol::Kind::Component(c) = &o.kind else {
                return Err(denied("invalid selected component"));
            };
            if c.sequences.len() > 1 {
                result.push(ChoiceRequirement {
                    original_hash: hash.clone(),
                    component_iri: id.clone(),
                    available_sequences: c.sequences.iter().cloned().collect(),
                    required_sequences: roots.get(&id).cloned().unwrap_or_default(),
                });
            }
            for f in &c.features {
                let o = doc
                    .objects()
                    .get(f)
                    .ok_or_else(|| denied("unresolved selected feature"))?;
                if let sbol::Kind::Feature(sbol::FeatureRecord {
                    kind: sbol::FeatureKind::SubComponent { instance_of, .. },
                    ..
                }) = &o.kind
                {
                    pending.push(instance_of.clone());
                }
            }
        }
    }
    Ok(result)
}
