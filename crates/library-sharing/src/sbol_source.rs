//! Exact-original SBOL documentary scope. Never rewrites native revisions or
//! publishes an opaque graph merely because one projected plasmid refers to it.
use crate::{Error, Workspace};
use collagene_design_exchange::{self as sbol, FeatureKind, Kind};
use std::collections::{BTreeMap, BTreeSet, VecDeque};

pub const KIND: &str = "sbol-design-exchange-v1";
const SOURCE_KEYS: [&str; 4] = [
    "sbol:component_iri",
    "sbol:sequence_iri",
    "sbol:source_hash",
    "sbol:exchange_profile",
];
pub(crate) type Selections = BTreeMap<String, BTreeSet<String>>;
pub(crate) type Claims = BTreeMap<String, Selections>;
#[derive(Clone, Copy)]
enum Expected {
    Component,
    Sequence,
    Feature,
    Location,
    Constraint,
}
fn denied(reason: &'static str) -> Error {
    Error::SbolSource(reason)
}
struct Closure<'a> {
    doc: &'a sbol::Document,
    selections: &'a Selections,
    pending: VecDeque<(String, Expected)>,
    kept: BTreeSet<String>,
    edges: usize,
    limit: usize,
}
impl<'a> Closure<'a> {
    fn edge(&mut self, id: &str, kind: Expected) -> Result<(), Error> {
        self.edges = self
            .edges
            .checked_add(1)
            .filter(|n| *n <= self.limit)
            .ok_or(Error::TooLarge)?;
        self.pending.push_back((id.into(), kind));
        Ok(())
    }
    fn resolve(mut self) -> Result<(BTreeSet<String>, usize), Error> {
        while let Some((id, expected)) = self.pending.pop_front() {
            let doc = self.doc;
            let object=doc.objects().get(&id).ok_or_else(||denied("selected structural reference is unresolved; no external SBOL document is fetched"))?;
            if !matches!(
                (&object.kind, expected),
                (Kind::Component(_), Expected::Component)
                    | (Kind::Sequence(_), Expected::Sequence)
                    | (Kind::Feature(_), Expected::Feature)
                    | (Kind::Location(_), Expected::Location)
                    | (Kind::Constraint(_), Expected::Constraint)
            ) {
                return Err(denied(
                    "selected structural reference has an unsupported or conflicting object class",
                ));
            }
            if object.classes.len() != 1 || !object.retained_properties.is_empty() {
                return Err(denied(
                    "selected SBOL objects have additional classes or properties without a publication reference contract",
                ));
            }
            if !self.kept.insert(id.clone()) {
                continue;
            }
            match &object.kind {
                Kind::Component(c) => {
                    if c.sequences.len() > 1
                        && self.selections.get(&id).is_none_or(|selected| {
                            c.sequences
                                .iter()
                                .any(|sequence| !selected.contains(sequence))
                        })
                    {
                        return Err(denied(
                            "unselected alternative SBOL sequences would be published; exact-original sharing requires explicit retained molecule selections for each alternative",
                        ));
                    }
                    for sequence in &c.sequences {
                        self.edge(sequence, Expected::Sequence)?;
                    }
                    for feature in &c.features {
                        self.edge(feature, Expected::Feature)?;
                    }
                    for constraint in &c.constraints {
                        self.edge(constraint, Expected::Constraint)?;
                    }
                }
                Kind::Feature(f) => {
                    for location in &f.locations {
                        self.edge(location, Expected::Location)?;
                    }
                    if let FeatureKind::SubComponent {
                        instance_of,
                        source_locations,
                        ..
                    } = &f.kind
                    {
                        self.edge(instance_of, Expected::Component)?;
                        for location in source_locations {
                            self.edge(location, Expected::Location)?;
                        }
                    }
                }
                Kind::Location(l) => self.edge(&l.sequence, Expected::Sequence)?,
                Kind::Constraint(c) => {
                    self.edge(&c.subject, Expected::Feature)?;
                    self.edge(&c.object, Expected::Feature)?;
                }
                Kind::Sequence(_) => {}
                _ => {
                    return Err(denied(
                        "unsupported SBOL object in the selected structural closure",
                    ));
                }
            }
        }
        Ok((self.kept, self.edges))
    }
}
/// Read claims from exactly the already-derived native revision closure. A
/// source hash is documentary provenance, not proof that current edited DNA is
/// still the original SBOL sequence, and never an executable/material receipt.
pub(crate) fn claims(
    source: &Workspace,
    revisions: &BTreeSet<String>,
    edge_limit: usize,
) -> Result<(Claims, usize), Error> {
    let mut claims: BTreeMap<String, Selections> = BTreeMap::new();
    let mut edges = 0usize;
    for revision in revisions {
        let molecule = source
            .revisions
            .get(revision)
            .ok_or(Error::Unresolved("SBOL native revision"))?;
        if !SOURCE_KEYS
            .iter()
            .any(|key| molecule.metadata.contains_key(*key))
        {
            continue;
        }
        let fields: Vec<_> = SOURCE_KEYS
            .iter()
            .map(|key| {
                molecule
                    .metadata
                    .get(*key)
                    .ok_or_else(|| denied("selected molecule has an incomplete SBOL source claim"))
            })
            .collect::<Result<_, _>>()?;
        let component = fields[0];
        let sequence = fields[1];
        let hash = fields[2];
        let profile = fields[3];
        if profile != sbol::PROFILE
            || hash.len() != 64
            || !hash
                .bytes()
                .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
        {
            return Err(denied(
                "selected molecule has an unknown SBOL source profile or invalid source hash",
            ));
        }
        edges = edges
            .checked_add(1)
            .filter(|n| *n <= edge_limit)
            .ok_or(Error::TooLarge)?;
        claims
            .entry(hash.clone())
            .or_default()
            .entry(component.clone())
            .or_default()
            .insert(sequence.clone());
    }
    Ok((claims, edges))
}
pub(crate) fn required_assets(
    source: &Workspace,
    revisions: &BTreeSet<String>,
    edge_limit: usize,
) -> Result<(BTreeSet<String>, usize), Error> {
    let (claims, mut edges) = claims(source, revisions, edge_limit)?;
    let mut assets = BTreeSet::new();
    for (hash, selections) in claims {
        let Some(attachment) = source.attachments.get(&hash) else {
            let (required, used) = crate::sbol_subset::received_assets(
                source,
                &hash,
                &selections,
                edge_limit.saturating_sub(edges),
            )?;
            assets.extend(required);
            edges = edges
                .checked_add(used)
                .filter(|v| *v <= edge_limit)
                .ok_or(Error::TooLarge)?;
            continue;
        };
        if attachment.kind != KIND || attachment.media_type != "application/ld+json" {
            return Err(denied(
                "SBOL source hash refers to a different attachment contract",
            ));
        }
        let doc = sbol::import(&attachment.bytes).map_err(|_| {
            denied("the selected original SBOL graph is not a valid bounded structural profile")
        })?;
        if doc.revision() != hash {
            return Err(denied(
                "SBOL original graph hash disagrees with the selected source claim",
            ));
        }
        let mut closure = Closure {
            doc: &doc,
            selections: &selections,
            pending: VecDeque::new(),
            kept: BTreeSet::new(),
            edges: 0,
            limit: edge_limit.saturating_sub(edges),
        };
        for (component, sequences) in &selections {
            let object = doc.objects().get(component).ok_or_else(|| {
                denied("selected component is absent from its original SBOL graph")
            })?;
            let Kind::Component(c) = &object.kind else {
                return Err(denied("selected component IRI does not name a Component"));
            };
            for sequence in sequences {
                if !c.sequences.contains(sequence) {
                    return Err(denied(
                        "selected Sequence is not declared by its selected Component",
                    ));
                }
                closure.edge(sequence, Expected::Sequence)?;
            }
            closure.edge(component, Expected::Component)?;
        }
        let (kept, used) = closure.resolve()?;
        edges = edges
            .checked_add(used)
            .filter(|n| *n <= edge_limit)
            .ok_or(Error::TooLarge)?;
        if kept.len() != doc.objects().len() {
            return Err(denied(
                "the original SBOL file contains unrelated graph objects; selected subgraph export or explicit omission/loss review is required and not yet implemented",
            ));
        }
        assets.insert(hash);
    }
    Ok((assets, edges))
}

/// Local-only documentary scope: preserve exact original bytes for retained
/// native claims even when that original contains other graph objects. This is
/// not a permitted-publication inventory and must never authorize an upload.
pub(crate) fn required_private_assets(
    source: &Workspace,
    revisions: &BTreeSet<String>,
    edge_limit: usize,
) -> Result<BTreeSet<String>, Error> {
    let (claims, mut used) = claims(source, revisions, edge_limit)?;
    let mut assets = BTreeSet::new();
    for (hash, selections) in claims {
        if let Some(a) = source.attachments.get(&hash) {
            if a.kind != KIND
                || a.media_type != "application/ld+json"
                || blake3::hash(&a.bytes).to_hex().as_str() != hash
            {
                return Err(denied(
                    "private original SBOL attachment kind/hash mismatch",
                ));
            }
            // Bounded syntax/profile validation retains opaque unsupported graph
            // objects as documentary bytes; it grants no native material status.
            let doc = sbol::import(&a.bytes).map_err(|_| {
                denied("private original SBOL source is not a bounded supported document")
            })?;
            used = used
                .checked_add(doc.objects().len())
                .filter(|n| *n <= edge_limit)
                .ok_or(Error::TooLarge)?;
            assets.insert(hash);
        } else {
            let (required, count) = crate::sbol_subset::received_assets(
                source,
                &hash,
                &selections,
                edge_limit.saturating_sub(used),
            )?;
            used = used
                .checked_add(count)
                .filter(|n| *n <= edge_limit)
                .ok_or(Error::TooLarge)?;
            assets.extend(required);
        }
    }
    Ok(assets)
}
