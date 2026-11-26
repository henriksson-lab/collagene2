use crate::*;
use serde::de::{self, DeserializeSeed, MapAccess, SeqAccess, Visitor};
use serde_json::{Map, Number};
use std::{
    collections::{BTreeMap, BTreeSet, VecDeque},
    fmt,
};
struct Checked<'a> {
    count: &'a mut usize,
    depth: usize,
}
impl<'de> DeserializeSeed<'de> for Checked<'_> {
    type Value = Value;
    fn deserialize<D: de::Deserializer<'de>>(self, d: D) -> std::result::Result<Value, D::Error> {
        *self.count += 1;
        if *self.count > MAX_VALUES || self.depth > 32 {
            return Err(de::Error::custom("JSON value/depth bound"));
        }
        d.deserialize_any(self)
    }
}
impl<'de> Visitor<'de> for Checked<'_> {
    type Value = Value;
    fn expecting(&self, f: &mut fmt::Formatter) -> fmt::Result {
        f.write_str("bounded duplicate-free JSON")
    }
    fn visit_bool<E: de::Error>(self, v: bool) -> std::result::Result<Value, E> {
        Ok(Value::Bool(v))
    }
    fn visit_i64<E: de::Error>(self, v: i64) -> std::result::Result<Value, E> {
        Ok(Value::Number(v.into()))
    }
    fn visit_u64<E: de::Error>(self, v: u64) -> std::result::Result<Value, E> {
        Ok(Value::Number(v.into()))
    }
    fn visit_f64<E: de::Error>(self, v: f64) -> std::result::Result<Value, E> {
        Number::from_f64(v)
            .map(Value::Number)
            .ok_or_else(|| E::custom("nonfinite number"))
    }
    fn visit_str<E: de::Error>(self, v: &str) -> std::result::Result<Value, E> {
        Ok(Value::String(v.into()))
    }
    fn visit_string<E: de::Error>(self, v: String) -> std::result::Result<Value, E> {
        Ok(Value::String(v))
    }
    fn visit_unit<E: de::Error>(self) -> std::result::Result<Value, E> {
        Ok(Value::Null)
    }
    fn visit_none<E: de::Error>(self) -> std::result::Result<Value, E> {
        Ok(Value::Null)
    }
    fn visit_seq<A: SeqAccess<'de>>(self, mut a: A) -> std::result::Result<Value, A::Error> {
        let mut out = vec![];
        while let Some(v) = a.next_element_seed(Checked {
            count: self.count,
            depth: self.depth + 1,
        })? {
            out.push(v);
        }
        Ok(Value::Array(out))
    }
    fn visit_map<A: MapAccess<'de>>(self, mut a: A) -> std::result::Result<Value, A::Error> {
        let mut out = Map::new();
        while let Some(k) = a.next_key::<String>()? {
            if out.contains_key(&k) {
                return Err(de::Error::custom("duplicate JSON property"));
            }
            let v = a.next_value_seed(Checked {
                count: self.count,
                depth: self.depth + 1,
            })?;
            out.insert(k, v);
        }
        Ok(Value::Object(out))
    }
}
pub(crate) fn iri(value: &str) -> Result<()> {
    if value.len() > 2048
        || value
            .chars()
            .any(|c| c.is_whitespace() || c.is_control() || "<>\"{}|\\^`".contains(c))
    {
        return Err(invalid("invalid or overlong absolute IRI"));
    }
    let (scheme, rest) = value
        .split_once(':')
        .ok_or_else(|| invalid("relative IRI is outside the profile"))?;
    if scheme.is_empty()
        || !scheme.as_bytes()[0].is_ascii_alphabetic()
        || !scheme
            .bytes()
            .all(|b| b.is_ascii_alphanumeric() || b"+.-".contains(&b))
        || rest.is_empty()
    {
        return Err(invalid("invalid absolute IRI"));
    }
    if let Some(authority) = rest.strip_prefix("//") {
        if authority
            .split(['/', '?', '#'])
            .next()
            .is_none_or(|s| s.is_empty())
        {
            return Err(invalid("IRI authority is empty"));
        }
    }
    let b = value.as_bytes();
    let mut i = 0;
    while i < b.len() {
        if b[i] == b'%' {
            if i + 2 >= b.len() || !b[i + 1].is_ascii_hexdigit() || !b[i + 2].is_ascii_hexdigit() {
                return Err(invalid("invalid IRI percent escape"));
            }
            i += 2;
        }
        i += 1;
    }
    Ok(())
}
fn context(v: Option<&Value>) -> Result<BTreeMap<String, String>> {
    let mut out = BTreeMap::new();
    let Some(v) = v else { return Ok(out) };
    let map = v
        .as_object()
        .ok_or(Error::Json("remote/array contexts are not resolved"))?;
    for (k, v) in map {
        let expected = match k.as_str() {
            "sbol" => SBOL,
            "prov" => "http://www.w3.org/ns/prov#",
            "xsd" => "http://www.w3.org/2001/XMLSchema#",
            "rdf" => "http://www.w3.org/1999/02/22-rdf-syntax-ns#",
            _ => {
                return Err(Error::Json(
                    "context must contain only pinned prefix definitions",
                ));
            }
        };
        if v.as_str() != Some(expected) {
            return Err(Error::Json("context prefix redefinition is forbidden"));
        }
        out.insert(k.clone(), expected.into());
    }
    Ok(out)
}
fn expand(s: &str, prefixes: &BTreeMap<String, String>) -> Result<String> {
    let expanded = if let Some((p, t)) = s.split_once(':') {
        prefixes
            .get(p)
            .map(|base| format!("{base}{t}"))
            .unwrap_or_else(|| s.into())
    } else {
        s.into()
    };
    iri(&expanded)?;
    Ok(expanded)
}
fn item(v: &Value, prefixes: &BTreeMap<String, String>) -> Result<Value> {
    if let Some(o) = v.as_object() {
        if o.len() == 1 && o.contains_key("@id") {
            return Ok(
                serde_json::json!({"@id":expand(o["@id"].as_str().ok_or(Error::Json("reference IRI must be a string"))?,prefixes)?}),
            );
        }
        if o.contains_key("@value")
            && o.keys()
                .all(|k| matches!(k.as_str(), "@value" | "@type" | "@language"))
            && !(o.contains_key("@type") && o.contains_key("@language"))
        {
            if !matches!(
                &o["@value"],
                Value::String(_) | Value::Number(_) | Value::Bool(_)
            ) {
                return Err(Error::Json("only scalar RDF literals are supported"));
            }
            let mut copy = o.clone();
            if let Some(t) = copy.get_mut("@type") {
                *t = Value::String(expand(
                    t.as_str().ok_or(Error::Json("literal datatype"))?,
                    prefixes,
                )?);
            }
            if copy
                .get("@language")
                .is_some_and(|v| v.as_str().is_none_or(|s| s.is_empty()))
            {
                return Err(Error::Json("literal language"));
            }
            return Ok(Value::Object(copy));
        }
        return Err(Error::Json(
            "nested definitions, lists, sets, reverse and scoped contexts are outside the profile",
        ));
    }
    if matches!(v, Value::String(_) | Value::Number(_) | Value::Bool(_)) {
        return Ok(serde_json::json!({"@value":v}));
    }
    Err(Error::Json(
        "null/array property items are outside the profile",
    ))
}
fn node(v: &Value, prefixes: &BTreeMap<String, String>) -> Result<Value> {
    let map = v
        .as_object()
        .ok_or(Error::Json("graph nodes must be objects"))?;
    let mut out = Map::new();
    let id = expand(
        map.get("@id")
            .and_then(Value::as_str)
            .ok_or(Error::Json("every graph node needs a stable absolute @id"))?,
        prefixes,
    )?;
    out.insert("@id".into(), Value::String(id));
    let class = map
        .get("@type")
        .ok_or(Error::Json("every graph node needs @type"))?;
    let types = class
        .as_array()
        .map(|a| a.iter().collect())
        .unwrap_or_else(|| vec![class]);
    let mut classes = vec![];
    for t in types {
        classes.push(Value::String(expand(
            t.as_str().ok_or(Error::Json("node type IRI"))?,
            prefixes,
        )?));
    }
    if classes.is_empty() {
        return Err(Error::Json("empty node types"));
    }
    out.insert("@type".into(), Value::Array(classes));
    for (k, v) in map {
        if matches!(k.as_str(), "@id" | "@type") {
            continue;
        }
        if k.starts_with('@') {
            return Err(Error::Json("unsupported JSON-LD keyword"));
        }
        let key = expand(k, prefixes)?;
        if out.contains_key(&key) {
            return Err(Error::Json("aliased duplicate predicate"));
        }
        let values = v
            .as_array()
            .map(|a| a.iter().collect())
            .unwrap_or_else(|| vec![v]);
        let values = values
            .into_iter()
            .map(|v| item(v, prefixes))
            .collect::<Result<Vec<_>>>()?;
        out.insert(key, Value::Array(values));
    }
    Ok(Value::Object(out))
}
fn key(name: &str) -> String {
    format!("{SBOL}{name}")
}
fn vals<'a>(m: &'a Map<String, Value>, name: &str) -> &'a [Value] {
    m.get(&key(name))
        .and_then(Value::as_array)
        .map(Vec::as_slice)
        .unwrap_or(&[])
}
fn one<'a>(m: &'a Map<String, Value>, name: &str, required: bool) -> Result<Option<&'a Value>> {
    let v = vals(m, name);
    if v.len() > 1 || (required && v.is_empty()) {
        return Err(invalid(format!("{name} cardinality")));
    }
    Ok(v.first())
}
fn reference(v: &Value) -> Result<String> {
    v.as_object()
        .filter(|o| o.len() == 1)
        .and_then(|o| o.get("@id"))
        .and_then(Value::as_str)
        .map(str::to_owned)
        .ok_or_else(|| invalid("expected an IRI reference, not a literal"))
}
fn refs(m: &Map<String, Value>, name: &str) -> Result<Vec<String>> {
    let out = vals(m, name)
        .iter()
        .map(reference)
        .collect::<Result<Vec<_>>>()?;
    let unique = out.iter().collect::<BTreeSet<_>>();
    if unique.len() != out.len() {
        return Err(invalid(format!("duplicate {name} values")));
    }
    Ok(out)
}
fn opt_ref(m: &Map<String, Value>, name: &str, required: bool) -> Result<Option<String>> {
    one(m, name, required)?.map(reference).transpose()
}
fn text(m: &Map<String, Value>, name: &str) -> Result<Option<String>> {
    one(m, name, false)?
        .map(|v| {
            let o = v.as_object().ok_or_else(|| invalid("text literal"))?;
            if o.contains_key("@language")
                || o.get("@type")
                    .is_some_and(|t| t.as_str() != Some("http://www.w3.org/2001/XMLSchema#string"))
            {
                return Err(invalid(format!(
                    "localized/typed {name} text is retained only outside this structural profile"
                )));
            }
            o.get("@value")
                .and_then(Value::as_str)
                .map(str::to_owned)
                .ok_or_else(|| invalid(format!("{name} is not string text")))
        })
        .transpose()
}
fn integer(m: &Map<String, Value>, name: &str, required: bool) -> Result<Option<i64>> {
    one(m, name, required)?
        .map(|v| {
            let o = v.as_object().ok_or_else(|| invalid("integer literal"))?;
            if o.contains_key("@language")
                || o.get("@type")
                    .is_some_and(|t| t.as_str() != Some("http://www.w3.org/2001/XMLSchema#integer"))
            {
                return Err(invalid("integer datatype"));
            }
            let v = o.get("@value").ok_or_else(|| invalid("integer value"))?;
            v.as_i64()
                .or_else(|| {
                    if o.get("@type").is_some() {
                        v.as_str().and_then(|s| s.parse().ok())
                    } else {
                        None
                    }
                })
                .ok_or_else(|| invalid(format!("{name} must be a bounded integer")))
        })
        .transpose()
}
fn namespace(m: &Map<String, Value>, id: &str) -> Result<String> {
    let ns = opt_ref(m, "hasNamespace", true)?.unwrap();
    if !(ns.starts_with("https://") || ns.starts_with("http://")) {
        return Err(invalid(
            "hasNamespace must be an absolute HTTP(S) URL in this profile",
        ));
    }
    if (id.starts_with("https://") || id.starts_with("http://")) && !id.starts_with(&ns) {
        return Err(invalid(
            "TopLevel identity must have declared namespace prefix",
        ));
    }
    Ok(ns)
}
fn typed(v: &Value) -> Result<Object> {
    let m = v.as_object().ok_or(Error::Json("graph object"))?;
    let id = m["@id"].as_str().unwrap();
    let classes = m["@type"]
        .as_array()
        .unwrap()
        .iter()
        .map(|v| v.as_str().unwrap().to_owned())
        .collect::<Vec<_>>();
    let identity = Identified {
        iri: id.into(),
        name: text(m, "name")?,
        display_id: text(m, "displayId")?,
        description: text(m, "description")?,
    };
    let known = classes
        .iter()
        .filter_map(|c| c.strip_prefix(SBOL))
        .filter(|t| {
            matches!(
                *t,
                "Component"
                    | "Sequence"
                    | "SequenceFeature"
                    | "SubComponent"
                    | "Range"
                    | "Cut"
                    | "EntireSequence"
                    | "Constraint"
                    | "Collection"
            )
        })
        .collect::<Vec<_>>();
    if known.len() > 1 {
        return Err(invalid("conflicting or repeated supported SBOL classes"));
    }
    let mut consumed = vec!["name", "displayId", "description"];
    let kind = match known.first().copied() {
        Some("Component") => {
            consumed.extend([
                "type",
                "role",
                "hasSequence",
                "hasFeature",
                "hasConstraint",
                "hasNamespace",
            ]);
            let types = refs(m, "type")?;
            if types.is_empty() {
                return Err(invalid("Component.type is required"));
            }
            Kind::Component(Component {
                types,
                roles: refs(m, "role")?,
                sequences: refs(m, "hasSequence")?,
                features: refs(m, "hasFeature")?,
                constraints: refs(m, "hasConstraint")?,
                namespace: namespace(m, id)?,
            })
        }
        Some("Sequence") => {
            consumed.extend(["elements", "encoding", "hasNamespace"]);
            let elements = text(m, "elements")?;
            let encoding = opt_ref(m, "encoding", elements.is_some())?;
            Kind::Sequence(SequenceRecord {
                elements,
                encoding,
                namespace: namespace(m, id)?,
            })
        }
        Some(t @ ("SequenceFeature" | "SubComponent")) => {
            consumed.extend(["role", "orientation", "hasLocation"]);
            let locations = refs(m, "hasLocation")?;
            let kind = if t == "SequenceFeature" {
                if locations.is_empty() {
                    return Err(invalid("SequenceFeature needs locations"));
                }
                FeatureKind::SequenceFeature
            } else {
                consumed.extend(["instanceOf", "sourceLocation", "roleIntegration"]);
                FeatureKind::SubComponent {
                    instance_of: opt_ref(m, "instanceOf", true)?.unwrap(),
                    source_locations: refs(m, "sourceLocation")?,
                    role_integration: opt_ref(m, "roleIntegration", false)?,
                }
            };
            let roles = refs(m, "role")?;
            if let FeatureKind::SubComponent {
                role_integration, ..
            } = &kind
            {
                if (!roles.is_empty() && role_integration.is_none())
                    || role_integration.as_ref().is_some_and(|v| {
                        !matches!(v.strip_prefix(SBOL), Some("mergeRoles" | "overrideRoles"))
                    })
                {
                    return Err(invalid(
                        "SubComponent roleIntegration is required for explicit roles and must be mergeRoles/overrideRoles",
                    ));
                }
            }
            Kind::Feature(FeatureRecord {
                kind,
                roles,
                orientation: opt_ref(m, "orientation", false)?,
                locations,
            })
        }
        Some(t @ ("Range" | "Cut" | "EntireSequence")) => {
            consumed.extend(["hasSequence", "orientation", "order"]);
            let kind = match t {
                "Range" => {
                    consumed.extend(["start", "end"]);
                    let a = integer(m, "start", true)?.unwrap();
                    let b = integer(m, "end", true)?.unwrap();
                    if a <= 0 || b < a {
                        return Err(invalid(
                            "Range requires1-based positive inclusive start<=end",
                        ));
                    }
                    LocationKind::Range {
                        start: usize::try_from(a).map_err(|_| Error::Limit)?,
                        end: usize::try_from(b).map_err(|_| Error::Limit)?,
                    }
                }
                "Cut" => {
                    consumed.push("at");
                    let at = integer(m, "at", true)?.unwrap();
                    LocationKind::Cut {
                        at: usize::try_from(at)
                            .map_err(|_| invalid("Cut boundary must be nonnegative"))?,
                    }
                }
                _ => LocationKind::EntireSequence,
            };
            Kind::Location(LocationRecord {
                kind,
                sequence: opt_ref(m, "hasSequence", true)?.unwrap(),
                orientation: opt_ref(m, "orientation", false)?,
                order: integer(m, "order", false)?,
            })
        }
        Some("Constraint") => {
            consumed.extend(["subject", "object", "restriction"]);
            let subject = opt_ref(m, "subject", true)?.unwrap();
            let object = opt_ref(m, "object", true)?.unwrap();
            if subject == object {
                return Err(invalid("Constraint endpoints must differ"));
            }
            Kind::Constraint(Constraint {
                subject,
                object,
                restriction: opt_ref(m, "restriction", true)?.unwrap(),
            })
        }
        Some("Collection") => {
            consumed.extend(["member", "hasNamespace"]);
            Kind::Collection {
                namespace: namespace(m, id)?,
                members: refs(m, "member")?,
            }
        }
        _ => Kind::Unsupported,
    };
    let consumed = consumed.into_iter().map(key).collect::<BTreeSet<_>>();
    let retained_properties = m
        .iter()
        .filter(|(k, _)| !k.starts_with('@') && !consumed.contains(*k))
        .map(|(k, v)| (k.clone(), v.clone()))
        .collect();
    Ok(Object {
        identity,
        classes,
        kind,
        retained_properties,
    })
}
fn note(issues: &mut Vec<Issue>, subject: &str, code: &str, message: &str) -> Result<()> {
    if issues.len() >= 1000 {
        return Err(Error::Limit);
    }
    issues.push(issue(subject, code, message));
    Ok(())
}
fn graph_validate(objects: &BTreeMap<String, Object>, issues: &mut Vec<Issue>) -> Result<()> {
    let mut owners = BTreeMap::<String, String>::new();
    let mut component_edges = BTreeMap::<String, Vec<String>>::new();
    let mut total_bases = 0usize;
    for (id, object) in objects {
        if !object.retained_properties.is_empty() {
            note(
                issues,
                id,
                "retained-properties",
                "Additional RDF properties remain in the exchange graph; native sequence projection does not implement them.",
            )?;
        }
        if matches!(object.kind, Kind::Unsupported) {
            note(
                issues,
                id,
                "unsupported-class",
                "Unsupported typed object is retained verbatim in the graph; no behavior or scientific authority is imported.",
            )?;
        }
        if object.classes.len() > 1 && !matches!(object.kind, Kind::Unsupported) {
            note(
                issues,
                id,
                "retained-classes",
                "Additional RDF classes are retained; complete class semantics are outside this structural profile.",
            )?;
        }
        if let Kind::Collection { members, .. } = &object.kind {
            for member in members {
                match objects.get(member) {
                    None => note(
                        issues,
                        id,
                        "unresolved-member",
                        "Collection member is absent; no external document is fetched.",
                    )?,
                    Some(o)
                        if matches!(
                            o.kind,
                            Kind::Feature(_) | Kind::Location(_) | Kind::Constraint(_)
                        ) =>
                    {
                        return Err(invalid(
                            "Collection member must be a TopLevel design object",
                        ));
                    }
                    _ => {}
                }
            }
        }
        if let Kind::Sequence(s) = &object.kind {
            if let Some(text) = &s.elements {
                total_bases = total_bases
                    .checked_add(text.len())
                    .filter(|n| *n <= MAX_SEQUENCE_BASES)
                    .ok_or(Error::Limit)?;
            }
        }
        if let Kind::Component(c) = &object.kind {
            let basic = c
                .types
                .iter()
                .filter(|t| {
                    matches!(
                        t.as_str(),
                        DNA | RNA
                            | PROTEIN
                            | "https://identifiers.org/SBO:0000247"
                            | "https://identifiers.org/SBO:0000253"
                            | "https://identifiers.org/SBO:0000241"
                    )
                })
                .count();
            if basic > 1
                || (c.types.iter().any(|s| semantic_so_type(s) == Some(LINEAR))
                    && c.types
                        .iter()
                        .any(|s| semantic_so_type(s) == Some(CIRCULAR)))
                || (c
                    .types
                    .iter()
                    .any(|s| semantic_so_type(s) == Some("http://identifiers.org/SO:0000984"))
                    && c.types
                        .iter()
                        .any(|s| semantic_so_type(s) == Some("http://identifiers.org/SO:0000985")))
            {
                return Err(invalid(
                    "Component contains conflicting basic types/topologies/strandedness",
                ));
            }
            component_edges.insert(id.clone(), vec![]);
            for child in c.features.iter().chain(&c.constraints) {
                if owners.insert(child.clone(), id.clone()).is_some() {
                    return Err(invalid(
                        "a feature/constraint has multiple Component owners",
                    ));
                }
            }
            for seq in &c.sequences {
                match objects.get(seq) {
                    None => note(
                        issues,
                        id,
                        "unresolved-sequence",
                        "Referenced sequence is not present; it will not be downloaded.",
                    )?,
                    Some(o) if !matches!(o.kind, Kind::Sequence(_)) => {
                        return Err(invalid("hasSequence does not reference a Sequence"));
                    }
                    _ => {}
                }
            }
        }
    }
    for (id, object) in objects {
        if let Kind::Component(c) = &object.kind {
            for fid in &c.features {
                let Some(fobj) = objects.get(fid) else {
                    note(
                        issues,
                        id,
                        "unresolved-feature",
                        "Referenced feature is absent; selected projection must refuse or explicitly omit it.",
                    )?;
                    continue;
                };
                let Kind::Feature(f) = &fobj.kind else {
                    if matches!(fobj.kind, Kind::Unsupported) {
                        continue;
                    }
                    return Err(invalid("hasFeature references the wrong supported class"));
                };
                let mut all_locations = f.locations.clone();
                let child_sequences = if let FeatureKind::SubComponent {
                    instance_of,
                    source_locations,
                    ..
                } = &f.kind
                {
                    all_locations.extend(source_locations.iter().cloned());
                    match objects.get(instance_of) {
                        Some(o) => match &o.kind {
                            Kind::Component(child) => {
                                component_edges
                                    .get_mut(id)
                                    .unwrap()
                                    .push(instance_of.clone());
                                Some(&child.sequences)
                            }
                            _ => return Err(invalid("instanceOf must reference a Component")),
                        },
                        None => {
                            note(
                                issues,
                                fid,
                                "unresolved-instance",
                                "Referenced child Component is absent; hierarchy remains documentary only.",
                            )?;
                            None
                        }
                    }
                } else {
                    None
                };
                for lid in &all_locations {
                    if owners.insert(lid.clone(), fid.clone()).is_some() {
                        return Err(invalid("a Location has multiple feature owners/uses"));
                    }
                    let Some(lobj) = objects.get(lid) else {
                        note(
                            issues,
                            fid,
                            "unresolved-location",
                            "Referenced location is absent.",
                        )?;
                        continue;
                    };
                    let Kind::Location(l) = &lobj.kind else {
                        if matches!(lobj.kind, Kind::Unsupported) {
                            continue;
                        }
                        return Err(invalid(
                            "hasLocation/sourceLocation references the wrong supported class",
                        ));
                    };
                    let is_source = matches!(&f.kind,FeatureKind::SubComponent{source_locations,..}if source_locations.contains(lid));
                    let allowed = if is_source {
                        child_sequences
                    } else {
                        Some(&c.sequences)
                    };
                    if allowed.is_some_and(|sequences| !sequences.contains(&l.sequence)) {
                        return Err(invalid(
                            "Location sequence is not owned by its parent Component (or source child)",
                        ));
                    }
                    if let Some(s) = objects.get(&l.sequence) {
                        let Kind::Sequence(s) = &s.kind else {
                            return Err(invalid("Location.hasSequence references wrong class"));
                        };
                        if let Some(text) = &s.elements {
                            match l.kind {
                                LocationKind::Range { end, .. } if end > text.len() => {
                                    return Err(invalid("Range exceeds explicit Sequence length"));
                                }
                                LocationKind::Cut { at } if at > text.len() => {
                                    return Err(invalid("Cut exceeds explicit Sequence length"));
                                }
                                _ => {}
                            }
                        }
                    } else {
                        note(
                            issues,
                            lid,
                            "unresolved-location-sequence",
                            "Location sequence is absent; no remote resolution is performed.",
                        )?;
                    }
                }
            }
            for cid in &c.constraints {
                match objects.get(cid) {
                    None => note(
                        issues,
                        id,
                        "unresolved-constraint",
                        "Referenced Constraint is absent.",
                    )?,
                    Some(o) => match &o.kind {
                        Kind::Constraint(k) => {
                            if !c.features.contains(&k.subject) || !c.features.contains(&k.object) {
                                return Err(invalid(
                                    "Constraint endpoints must be distinct features of the same owning Component",
                                ));
                            }
                        }
                        Kind::Unsupported => {}
                        _ => return Err(invalid("hasConstraint references wrong class")),
                    },
                }
            }
        }
    }
    // Detect recursive instanceOf chains in O(V+E), without recursive traversal.
    let mut incoming = component_edges
        .keys()
        .map(|k| (k.clone(), 0usize))
        .collect::<BTreeMap<_, _>>();
    for edges in component_edges.values() {
        for child in edges {
            *incoming.get_mut(child).unwrap() += 1;
        }
    }
    let mut queue = incoming
        .iter()
        .filter(|(_, n)| **n == 0)
        .map(|(id, _)| id.clone())
        .collect::<VecDeque<_>>();
    let mut visited = 0;
    while let Some(id) = queue.pop_front() {
        visited += 1;
        for child in &component_edges[&id] {
            let n = incoming.get_mut(child).unwrap();
            *n -= 1;
            if *n == 0 {
                queue.push_back(child.clone());
            }
        }
    }
    if visited != component_edges.len() {
        return Err(invalid(
            "cyclic Component/SubComponent instanceOf hierarchy",
        ));
    }
    Ok(())
}
pub(crate) fn import(bytes: &[u8]) -> Result<Document> {
    if bytes.is_empty() || bytes.len() > MAX_BYTES {
        return Err(Error::Limit);
    }
    let mut deserializer = serde_json::Deserializer::from_slice(bytes);
    let mut count = 0;
    let root = Checked {
        count: &mut count,
        depth: 0,
    }
    .deserialize(&mut deserializer)
    .map_err(|_| Error::Json("invalid, duplicate, oversized or too deeply nested JSON"))?;
    deserializer
        .end()
        .map_err(|_| Error::Json("trailing JSON input"))?;
    let (graph, prefixes) = match &root {
        Value::Array(a) => (a, context(None)?),
        Value::Object(o) => {
            if !o
                .keys()
                .all(|k| matches!(k.as_str(), "@context" | "@graph"))
            {
                return Err(Error::Json(
                    "document wrapper may contain only @context and @graph",
                ));
            }
            (
                o.get("@graph")
                    .and_then(Value::as_array)
                    .ok_or(Error::Json("expected flat @graph array"))?,
                context(o.get("@context"))?,
            )
        }
        _ => return Err(Error::Json("expected expanded node array or graph wrapper")),
    };
    if graph.is_empty() || graph.len() > MAX_OBJECTS {
        return Err(Error::Limit);
    }
    let mut expanded = Vec::with_capacity(graph.len());
    let mut objects = BTreeMap::new();
    for value in graph {
        let raw = node(value, &prefixes)?;
        let object = typed(&raw)?;
        if objects
            .insert(object.identity.iri.clone(), object)
            .is_some()
        {
            return Err(invalid(
                "duplicate graph object identity; RDF node merging is outside this profile",
            ));
        }
        expanded.push(raw);
    }
    let mut issues = vec![];
    graph_validate(&objects, &mut issues)?;
    Ok(Document {
        original: bytes.to_vec(),
        objects,
        expanded,
        issues,
    })
}
