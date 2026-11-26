//! Deterministic typed structural subsets, not claims to reproduce original RDF.
use crate::{Document, Error, FeatureKind, Kind, MAX_OBJECTS, MAX_VALUES, Result, SBOL};
use serde::{Deserialize, Serialize};
use serde_json::Value;
use std::collections::{BTreeMap, BTreeSet, VecDeque};
pub const POLICY: &str = "collagene-sbol-selected-structural-subset-1";
pub type SequenceChoices = BTreeMap<String, BTreeSet<String>>;
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize, Default)]
#[serde(deny_unknown_fields)]
pub struct Request {
    /// Actual selected roots and required source sequences; library publication
    /// derives these from retained native claims, never from arbitrary graph text.
    pub roots: SequenceChoices,
    /// Explicit sequence choices for reachable components. Required for a child
    /// with alternatives; root choices must include every required root sequence.
    pub alternatives: SequenceChoices,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Omissions {
    pub unrelated_objects: usize,
    pub alternative_sequence_edges: usize,
}
#[derive(Debug, Clone)]
pub struct Export {
    pub bytes: Vec<u8>,
    pub original_hash: String,
    pub subset_hash: String,
    pub request: Request,
    pub omissions: Omissions,
    /// Local review data only; contains only retained identities, never omitted IDs.
    pub retained_objects: BTreeSet<String>,
}
fn denied(s: &str) -> Error {
    Error::Projection(format!("selected documentary subset: {s}"))
}
#[derive(Clone, Copy)]
enum Expected {
    Component,
    Sequence,
    Feature,
    Location,
    Constraint,
}
struct State<'a> {
    doc: &'a Document,
    request: &'a Request,
    pending: VecDeque<(String, Expected)>,
    kept: BTreeSet<String>,
    components: SequenceChoices,
    edges: usize,
    omitted_edges: usize,
}
impl State<'_> {
    fn edge(&mut self, id: &str, kind: Expected) -> Result<()> {
        self.edges = self
            .edges
            .checked_add(1)
            .filter(|v| *v <= MAX_VALUES)
            .ok_or(Error::Limit)?;
        self.pending.push_back((id.into(), kind));
        Ok(())
    }
    fn resolve(&mut self) -> Result<()> {
        while let Some((id, expected)) = self.pending.pop_front() {
            let doc = self.doc;
            let object = doc.objects().get(&id).ok_or_else(|| {
                denied("unresolved structural reference; external documents are never fetched")
            })?;
            if !matches!(
                (&object.kind, expected),
                (Kind::Component(_), Expected::Component)
                    | (Kind::Sequence(_), Expected::Sequence)
                    | (Kind::Feature(_), Expected::Feature)
                    | (Kind::Location(_), Expected::Location)
                    | (Kind::Constraint(_), Expected::Constraint)
            ) {
                return Err(denied("conflicting/unsupported structural object class"));
            }
            if object.classes.len() != 1 || !object.retained_properties.is_empty() {
                return Err(denied(
                    "selected objects contain extra classes/properties without an omission/reference contract",
                ));
            }
            if !self.kept.insert(id.clone()) {
                continue;
            }
            match &object.kind {
                Kind::Component(component) => {
                    let required = self.request.roots.get(&id);
                    let chosen = if let Some(choice) = self.request.alternatives.get(&id) {
                        choice.clone()
                    } else if let Some(required) = required {
                        required.clone()
                    } else if component.sequences.len() <= 1 {
                        component.sequences.iter().cloned().collect()
                    } else {
                        return Err(denied(
                            "a nested component has alternative sequences; explicit choices are required",
                        ));
                    };
                    if required.is_some_and(|v| !v.is_subset(&chosen))
                        || chosen.iter().any(|v| !component.sequences.contains(v))
                        || (!component.sequences.is_empty() && chosen.is_empty())
                    {
                        return Err(denied(
                            "sequence choice does not cover required native claims or is absent/empty",
                        ));
                    }
                    self.omitted_edges = self
                        .omitted_edges
                        .checked_add(component.sequences.len() - chosen.len())
                        .ok_or(Error::Limit)?;
                    for sequence in &chosen {
                        self.edge(sequence, Expected::Sequence)?;
                    }
                    for feature in &component.features {
                        self.edge(feature, Expected::Feature)?;
                    }
                    for constraint in &component.constraints {
                        self.edge(constraint, Expected::Constraint)?;
                    }
                    self.components.insert(id, chosen);
                }
                Kind::Feature(feature) => {
                    for location in &feature.locations {
                        self.edge(location, Expected::Location)?;
                    }
                    if let FeatureKind::SubComponent {
                        instance_of,
                        source_locations,
                        ..
                    } = &feature.kind
                    {
                        self.edge(instance_of, Expected::Component)?;
                        for location in source_locations {
                            self.edge(location, Expected::Location)?;
                        }
                    }
                }
                Kind::Location(location) => self.edge(&location.sequence, Expected::Sequence)?,
                Kind::Constraint(constraint) => {
                    self.edge(&constraint.subject, Expected::Feature)?;
                    self.edge(&constraint.object, Expected::Feature)?;
                }
                Kind::Sequence(_) => (),
                _ => return Err(denied("unsupported structural closure object")),
            }
        }
        if self
            .request
            .alternatives
            .keys()
            .any(|id| !self.components.contains_key(id))
        {
            return Err(denied(
                "unused/unrelated component choice is not a publication root",
            ));
        }
        let allowed: BTreeSet<_> = self
            .components
            .values()
            .flat_map(|s| s.iter().cloned())
            .collect();
        if self.kept.iter().any(|id| {
            matches!(self.doc.objects()[id].kind, Kind::Sequence(_)) && !allowed.contains(id)
        }) {
            return Err(denied(
                "a structural sequence is not declared by a selected reachable component",
            ));
        }
        Ok(())
    }
}
fn check_internal_refs(value: &Value, document: &Document, kept: &BTreeSet<String>) -> Result<()> {
    match value {
        Value::Object(map) => {
            for (key, value) in map {
                if key == "@id" || key == "@type" {
                    let refs: Vec<&str> = match value {
                        Value::String(s) => vec![s.as_str()],
                        Value::Array(items) => items.iter().filter_map(Value::as_str).collect(),
                        _ => vec![],
                    };
                    if refs
                        .iter()
                        .any(|id| document.objects().contains_key(*id) && !kept.contains(*id))
                    {
                        return Err(denied(
                            "retained vocabulary reference names an omitted original object",
                        ));
                    }
                }
                check_internal_refs(value, document, kept)?;
            }
        }
        Value::Array(items) => {
            for item in items {
                check_internal_refs(item, document, kept)?;
            }
        }
        _ => (),
    }
    Ok(())
}
fn canonical(value: &Value) -> Value {
    match value {
        Value::Object(object) => {
            let sorted: BTreeMap<_, _> = object
                .iter()
                .map(|(k, v)| (k.clone(), canonical(v)))
                .collect();
            Value::Object(sorted.into_iter().collect())
        }
        Value::Array(values) => Value::Array(values.iter().map(canonical).collect()),
        other => other.clone(),
    }
}
/// Canonical only for this typed structural profile: object keys and graph nodes
/// are sorted; declared array ordering is preserved. This is not RDF URDNA2015.
pub fn build(document: &Document, request: &Request) -> Result<Export> {
    if request.roots.is_empty()
        || request.roots.len() > MAX_OBJECTS
        || request.alternatives.len() > MAX_OBJECTS
    {
        return Err(Error::Limit);
    }
    let choices = request.roots.iter().chain(&request.alternatives);
    let mut entries = 0usize;
    let mut choice_bytes = 0usize;
    for (component, sequences) in choices {
        crate::jsonld::iri(component)?;
        choice_bytes = choice_bytes
            .checked_add(component.len())
            .filter(|v| *v <= crate::MAX_BYTES)
            .ok_or(Error::Limit)?;
        if sequences.is_empty() {
            return Err(denied("root/alternative sequence choices must be nonempty"));
        }
        for sequence in sequences {
            crate::jsonld::iri(sequence)?;
            choice_bytes = choice_bytes
                .checked_add(sequence.len())
                .filter(|v| *v <= crate::MAX_BYTES)
                .ok_or(Error::Limit)?;
            entries = entries
                .checked_add(1)
                .filter(|v| *v <= MAX_VALUES)
                .ok_or(Error::Limit)?;
        }
    }
    let mut state = State {
        doc: document,
        request,
        pending: VecDeque::new(),
        kept: BTreeSet::new(),
        components: BTreeMap::new(),
        edges: 0,
        omitted_edges: 0,
    };
    for component in request.roots.keys() {
        state.edge(component, Expected::Component)?;
    }
    state.resolve()?;
    let mut nodes = BTreeMap::new();
    for raw in &document.expanded {
        let id = raw
            .get("@id")
            .and_then(Value::as_str)
            .ok_or_else(|| denied("expanded node identity missing"))?;
        if !state.kept.contains(id) {
            continue;
        }
        let mut node = raw.clone();
        if let Some(sequences) = state.components.get(id) {
            let component = match &document.objects()[id].kind {
                Kind::Component(c) => c,
                _ => unreachable!(),
            };
            let values = component
                .sequences
                .iter()
                .filter(|id| sequences.contains(*id))
                .map(|id| serde_json::json!({"@id":id}))
                .collect();
            if !component.sequences.is_empty() {
                node.as_object_mut()
                    .ok_or_else(|| denied("expanded object expected"))?
                    .insert(format!("{SBOL}hasSequence"), Value::Array(values));
            }
        }
        // Open vocabulary terms are preserved, but must not disclose an omitted
        // object identity from this very source through a non-structural edge.
        check_internal_refs(&node, document, &state.kept)?;
        nodes.insert(id.to_string(), canonical(&node));
    }
    let bytes = crate::write::encode(
        &serde_json::json!({"@graph":nodes.into_values().collect::<Vec<_>>()}),
    )?;
    let imported = crate::import(&bytes)?;
    if imported.objects().keys().cloned().collect::<BTreeSet<_>>() != state.kept {
        return Err(denied("canonical subset changed structural identities"));
    }
    Ok(Export {
        original_hash: document.revision(),
        subset_hash: blake3::hash(&bytes).to_hex().to_string(),
        bytes,
        request: request.clone(),
        omissions: Omissions {
            unrelated_objects: document.objects().len() - state.kept.len(),
            alternative_sequence_edges: state.omitted_edges,
        },
        retained_objects: state.kept,
    })
}
/// Independently check a received canonical subset without claiming access to,
/// or derivation from, the missing original graph.
pub fn validate_received(bytes: &[u8], request: &Request) -> Result<Document> {
    let document = crate::import(bytes)?;
    let subset = build(&document, request)?;
    if subset.omissions.unrelated_objects != 0
        || subset.omissions.alternative_sequence_edges != 0
        || subset.bytes != bytes
    {
        return Err(denied(
            "received subset is not the exact minimal canonical selected closure",
        ));
    }
    Ok(document)
}
