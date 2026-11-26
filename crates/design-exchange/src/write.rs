use crate::*;
use serde::Serialize;
use std::io::Write;
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ExportIds {
    pub namespace: String,
    pub component: String,
    pub sequence: String,
}
#[derive(Debug, Clone)]
pub struct Export {
    pub bytes: Vec<u8>,
    pub issues: Vec<Issue>,
}
struct Bounded(Vec<u8>);
impl Write for Bounded {
    fn write(&mut self, b: &[u8]) -> std::io::Result<usize> {
        if b.len() > MAX_BYTES.saturating_sub(self.0.len()) {
            return Err(std::io::Error::other("SBOL output bound"));
        }
        self.0.extend_from_slice(b);
        Ok(b.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}
pub(crate) fn encode(value: &impl Serialize) -> Result<Vec<u8>> {
    let mut out = Bounded(vec![]);
    serde_json::to_writer_pretty(&mut out, value).map_err(|_| Error::Limit)?;
    Ok(out.0)
}
fn node(id: &str, class: &str) -> serde_json::Map<String, Value> {
    let mut n = serde_json::Map::new();
    n.insert("@id".into(), Value::String(id.into()));
    n.insert(
        "@type".into(),
        serde_json::json!([format!("{SBOL}{class}")]),
    );
    n
}
fn refs(
    n: &mut serde_json::Map<String, Value>,
    property: &str,
    values: impl IntoIterator<Item = String>,
) {
    n.insert(
        format!("{SBOL}{property}"),
        Value::Array(
            values
                .into_iter()
                .map(|v| serde_json::json!({"@id":v}))
                .collect(),
        ),
    );
}
fn literal(n: &mut serde_json::Map<String, Value>, property: &str, value: impl Serialize) {
    n.insert(
        format!("{SBOL}{property}"),
        serde_json::json!([{"@value":value}]),
    );
}
fn parts(
    location: &Location,
    reverse: bool,
    depth: usize,
    out: &mut Vec<(usize, usize, bool)>,
) -> Result<()> {
    if depth > 32 || out.len() > MAX_OBJECTS {
        return Err(Error::Limit);
    }
    match location {
        Location::Range {
            start,
            end,
            before: false,
            after: false,
        } => out.push((*start, *end, reverse)),
        Location::Complement(inner) => parts(inner, !reverse, depth + 1, out)?,
        Location::Join(children) => {
            if reverse {
                for p in children.iter().rev() {
                    parts(p, reverse, depth + 1, out)?;
                }
            } else {
                for p in children {
                    parts(p, reverse, depth + 1, out)?;
                }
            }
        }
        _ => {
            return Err(Error::Projection(
                "fuzzy, unresolved or unordered native feature location is not an exact SBOL Range"
                    .into(),
            ));
        }
    }
    Ok(())
}
fn loss(issues: &mut Vec<Issue>, subject: &str, message: &str) -> Result<()> {
    if issues.len() >= 1000 {
        return Err(Error::Limit);
    }
    issues.push(issue(subject, "sbol-export-loss", message));
    Ok(())
}
pub(crate) fn export(
    molecule: &Molecule,
    ids: &ExportIds,
    policy: ProjectionPolicy,
) -> Result<Export> {
    collagene_model::bounds::model_preflight(&[molecule])?;
    molecule.validate()?;
    if molecule.sequence.len() > MAX_SEQUENCE_BASES || molecule.features.len() > MAX_OBJECTS / 2 {
        return Err(Error::Limit);
    }
    for iri in [&ids.namespace, &ids.component, &ids.sequence] {
        jsonld::iri(iri)?;
    }
    if !(ids.namespace.starts_with("https://") || ids.namespace.starts_with("http://"))
        || ids.component == ids.sequence
        || [&ids.component, &ids.sequence]
            .iter()
            .any(|id| !id.starts_with(&ids.namespace))
    {
        return Err(invalid(
            "export requires distinct stable Component/Sequence identities under explicit HTTP(S) namespace",
        ));
    }
    let mut component = node(&ids.component, "Component");
    refs(&mut component, "hasNamespace", [ids.namespace.clone()]);
    literal(&mut component, "name", &molecule.name);
    refs(&mut component, "hasSequence", [ids.sequence.clone()]);
    let alphabet = molecule.sequence.alphabet();
    let mut types = vec![
        match alphabet {
            Alphabet::Dna => DNA,
            Alphabet::Rna => RNA,
            Alphabet::Protein => PROTEIN,
        }
        .into(),
    ];
    if alphabet != Alphabet::Protein {
        types.push(
            if molecule.topology == Topology::Circular {
                CIRCULAR
            } else {
                LINEAR
            }
            .into(),
        );
    }
    refs(&mut component, "type", types);
    let mut sequence = node(&ids.sequence, "Sequence");
    refs(&mut sequence, "hasNamespace", [ids.namespace.clone()]);
    refs(
        &mut sequence,
        "encoding",
        [if alphabet == Alphabet::Protein {
            PROTEIN_ENCODING
        } else {
            NUCLEIC_ENCODING
        }
        .into()],
    );
    literal(&mut sequence, "elements", molecule.sequence.as_str());
    let mut nodes = vec![];
    let mut features = vec![];
    let mut issues = vec![];
    if !molecule.metadata.is_empty() || molecule.original_genbank.is_some() {
        loss(
            &mut issues,
            &ids.component,
            "Native metadata/original file/evidence are not promoted into SBOL design semantics or physical verification. Keep the original native archive separately.",
        )?;
    }
    for (index, feature) in molecule.features.iter().enumerate() {
        let fid = format!("{}/feature/{}", ids.component, index + 1);
        let mut ranges = vec![];
        if let Err(error) = parts(&feature.location, false, 0, &mut ranges) {
            loss(&mut issues, &fid, &format!("Feature omitted: {error}"))?;
            continue;
        }
        if nodes
            .len()
            .checked_add(ranges.len() + 1)
            .is_none_or(|n| n + 2 > MAX_OBJECTS)
        {
            return Err(Error::Limit);
        }
        if alphabet == Alphabet::Protein && ranges.iter().any(|(_, _, r)| *r) {
            loss(
                &mut issues,
                &fid,
                "Reverse-complement protein location is unsupported; feature omitted.",
            )?;
            continue;
        }
        let mut f = node(&fid, "SequenceFeature");
        literal(&mut f, "name", feature.label());
        let mut locations = vec![];
        for (i, (start, end, reverse)) in ranges.into_iter().enumerate() {
            let lid = format!("{fid}/location/{}", i + 1);
            let mut l = node(&lid, "Range");
            literal(&mut l, "start", start + 1);
            literal(&mut l, "end", end);
            literal(&mut l, "order", i as i64 + 1);
            refs(&mut l, "hasSequence", [ids.sequence.clone()]);
            refs(
                &mut l,
                "orientation",
                [if reverse { REVERSE } else { INLINE }.into()],
            );
            locations.push(lid);
            nodes.push(Value::Object(l));
        }
        refs(&mut f, "hasLocation", locations);
        let roles = feature
            .qualifiers
            .iter()
            .filter_map(|(k, v)| if k == "sbol:role" { v.clone() } else { None })
            .collect::<Vec<_>>();
        for role in &roles {
            jsonld::iri(role)?;
        }
        if !roles.is_empty() {
            refs(&mut f, "role", roles);
        }
        if feature.kind != "misc_feature"
            || feature
                .qualifiers
                .iter()
                .any(|(k, _)| !matches!(k.as_str(), "label" | "sbol:role"))
            || feature.needs_review
        {
            loss(
                &mut issues,
                &fid,
                "Native feature kind/qualifiers/review flags have no complete structural-profile equivalent; explicit location/name/role are exported, other fields remain native-only.",
            )?;
        }
        features.push(fid);
        nodes.push(Value::Object(f));
    }
    if policy == ProjectionPolicy::Strict && !issues.is_empty() {
        return Err(Error::Projection(format!(
            "{} explicit SBOL export losses require AllowOmissions review",
            issues.len()
        )));
    }
    if !features.is_empty() {
        refs(&mut component, "hasFeature", features);
    }
    nodes.push(Value::Object(component));
    nodes.push(Value::Object(sequence));
    let bytes = encode(&serde_json::json!({"@graph":nodes}))?;
    // Validate generated graph with the independent importer profile before return.
    crate::import(&bytes)?;
    Ok(Export { bytes, issues })
}
