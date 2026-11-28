use crate::*;
use std::collections::BTreeSet;
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ProjectionPolicy {
    Strict,
    AllowOmissions,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Selection {
    pub component: String,
    pub sequence: String,
    pub topology: Option<Topology>,
}
#[derive(Debug, Clone)]
pub struct Projection {
    pub molecule: Molecule,
    pub component: String,
    pub sequence: String,
    pub issues: Vec<Issue>,
    pub warnings: Vec<String>,
}
fn loss(issues: &mut Vec<Issue>, subject: &str, message: &str) -> Result<()> {
    if issues.len() >= 1000 {
        return Err(Error::Limit);
    }
    issues.push(issue(subject, "native-projection-loss", message));
    Ok(())
}
fn orientation(iri: &str) -> Result<bool> {
    match iri {
        INLINE | "http://sbols.org/v3#inline" => Ok(false),
        REVERSE | "http://sbols.org/v3#reverseComplement" => Ok(true),
        _ => Err(Error::Projection("unsupported orientation IRI".into())),
    }
}
fn retained(object: &Object, issues: &mut Vec<Issue>) -> Result<()> {
    if !object.retained_properties.is_empty() {
        loss(
            issues,
            &object.identity.iri,
            "Additional properties remain in the SBOL graph and are not native molecule semantics.",
        )?;
    }
    if object.classes.len() != 1 {
        loss(
            issues,
            &object.identity.iri,
            "Additional RDF classes remain in the SBOL graph and are not interpreted by the selected conversion.",
        )?;
    }
    Ok(())
}
fn feature(
    doc: &Document,
    id: &str,
    selected: &str,
    alphabet: Alphabet,
    issues: &mut Vec<Issue>,
) -> Result<Feature> {
    let object = doc
        .objects
        .get(id)
        .ok_or_else(|| Error::Projection("feature reference is unresolved".into()))?;
    let Kind::Feature(f) = &object.kind else {
        return Err(Error::Projection("unsupported feature class".into()));
    };
    retained(object, issues)?;
    let mut locations = f
        .locations
        .iter()
        .map(|lid| {
            let o = doc
                .objects
                .get(lid)
                .ok_or_else(|| Error::Projection("location reference is unresolved".into()))?;
            let Kind::Location(l) = &o.kind else {
                return Err(Error::Projection("unsupported location class".into()));
            };
            if l.sequence != selected {
                return Err(Error::Projection(
                    "feature locations belong to a different selected Sequence".into(),
                ));
            }
            Ok((o, l))
        })
        .collect::<Result<Vec<_>>>()?;
    if locations.is_empty() {
        return Err(Error::Projection(
            "feature has no explicit location on selected Sequence".into(),
        ));
    }
    if locations.len() > 1 {
        let orders = locations
            .iter()
            .filter_map(|(_, l)| l.order)
            .collect::<BTreeSet<_>>();
        if orders.len() != locations.len() {
            return Err(Error::Projection(
                "multipart location ordering is absent or ambiguous".into(),
            ));
        }
        locations.sort_by_key(|(_, l)| l.order.unwrap());
    }
    let mut parts = vec![];
    for (o, l) in locations {
        retained(o, issues)?;
        let feature_orientation = f.orientation.as_deref().map(orientation).transpose()?;
        let location_orientation = l.orientation.as_deref().map(orientation).transpose()?;
        let effective = match (feature_orientation, location_orientation) {
            (Some(a), Some(b)) if a != b => {
                return Err(Error::Projection(
                    "feature/location orientations conflict".into(),
                ));
            }
            (_, Some(value)) | (Some(value), None) => value,
            _ => {
                return Err(Error::Projection(
                    "orientation is unspecified; no implicit strand is invented".into(),
                ));
            }
        };
        if alphabet == Alphabet::Protein && effective {
            return Err(Error::Projection(
                "reverse-complement protein location is unsupported".into(),
            ));
        }
        let location=match l.kind{
   LocationKind::Range{start,end}=>Location::Range{start:start-1,end,before:false,after:false},
   LocationKind::Cut{..}=>return Err(Error::Projection("Cut boundary is retained in graph; native feature range conversion is not supported".into())),
   LocationKind::EntireSequence=>{
    let Some(Object{kind:Kind::Sequence(s),..})=doc.objects.get(selected)else{return Err(Error::Projection("EntireSequence reference is unresolved".into()));};
    Location::Range{start:0,end:s.elements.as_ref().ok_or_else(||Error::Projection("EntireSequence has no elements".into()))?.len(),before:false,after:false}
   }
  };
        parts.push(if effective {
            Location::Complement(Box::new(location))
        } else {
            location
        });
    }
    let location = if parts.len() == 1 {
        parts.pop().unwrap()
    } else {
        Location::Join(parts)
    };
    let mut qualifiers = vec![("sbol:identity".into(), Some(id.into()))];
    if let Some(name) = object
        .identity
        .name
        .as_ref()
        .or(object.identity.display_id.as_ref())
    {
        qualifiers.push(("label".into(), Some(name.clone())));
    }
    for role in &f.roles {
        qualifiers.push(("sbol:role".into(), Some(role.clone())));
    }
    if let Some(text) = &object.identity.description {
        qualifiers.push(("note".into(), Some(text.clone())));
    }
    if let Some(text) = &object.identity.display_id {
        qualifiers.push(("sbol:display_id".into(), Some(text.clone())));
    }
    if let FeatureKind::SubComponent { instance_of, .. } = &f.kind {
        qualifiers.push(("sbol:instance_of".into(), Some(instance_of.clone())));
        loss(
            issues,
            id,
            "SubComponent is only an explicitly located annotation here. Child hierarchy, role integration and source locations remain in the exchange graph; no child DNA is synthesized or flattened.",
        )?;
    }
    Ok(Feature {
        kind: "misc_feature".into(),
        location,
        qualifiers,
        needs_review: true,
    })
}
pub(crate) fn project(
    doc: &Document,
    selection: &Selection,
    policy: ProjectionPolicy,
) -> Result<Projection> {
    let object = doc
        .objects
        .get(&selection.component)
        .ok_or_else(|| Error::Projection("selected Component is absent".into()))?;
    let Kind::Component(component) = &object.kind else {
        return Err(Error::Projection(
            "selection is not a supported Component".into(),
        ));
    };
    if !component.sequences.contains(&selection.sequence) {
        return Err(Error::Projection(
            "selected Sequence is not linked by selected Component.hasSequence".into(),
        ));
    }
    let seq_object = doc
        .objects
        .get(&selection.sequence)
        .ok_or_else(|| Error::Projection("selected Sequence is unresolved".into()))?;
    let Kind::Sequence(sequence) = &seq_object.kind else {
        return Err(Error::Projection("selection is not a Sequence".into()));
    };
    let alphabets = component
        .types
        .iter()
        .filter_map(|t| match t.as_str() {
            DNA => Some(Alphabet::Dna),
            RNA => Some(Alphabet::Rna),
            PROTEIN => Some(Alphabet::Protein),
            _ => None,
        })
        .collect::<Vec<_>>();
    if alphabets.len() != 1 {
        return Err(Error::Projection(
            "selected Component needs one supported explicit DNA/RNA/protein category".into(),
        ));
    }
    let alphabet = alphabets[0];
    let expected = if alphabet == Alphabet::Protein {
        PROTEIN_ENCODING
    } else {
        NUCLEIC_ENCODING
    };
    if sequence.encoding.as_deref() != Some(expected) {
        return Err(Error::Projection(
            "selected Sequence encoding is unsupported or inconsistent with Component type".into(),
        ));
    }
    let text = sequence.elements.as_deref().ok_or_else(|| {
        Error::Projection(
            "selected Sequence has no explicit elements; hierarchy assembly is not inferred".into(),
        )
    })?;
    if text.bytes().any(|b| b.is_ascii_whitespace()) || !text.is_ascii() {
        return Err(Error::Projection(
            "selected IUPAC elements must contain only literal ASCII symbols, without whitespace"
                .into(),
        ));
    }
    let declared = if component
        .types
        .iter()
        .any(|t| semantic_so_type(t) == Some(CIRCULAR))
    {
        Some(Topology::Circular)
    } else if component
        .types
        .iter()
        .any(|t| semantic_so_type(t) == Some(LINEAR))
        || alphabet != Alphabet::Dna
    {
        Some(Topology::Linear)
    } else {
        None
    };
    let topology = match (declared, selection.topology) {
        (Some(a), Some(b)) if a != b => {
            return Err(Error::Projection(
                "requested topology conflicts with declared Component topology".into(),
            ));
        }
        (Some(a), _) => a,
        (None, Some(b)) => b,
        _ => return Err(Error::Projection(
            "DNA topology is unspecified; choose explicitly without modifying the retained graph"
                .into(),
        )),
    };
    let mut issues = vec![];
    let mut warnings = vec![];
    retained(object, &mut issues)?;
    retained(seq_object, &mut issues)?;
    if component.sequences.len() > 1 {
        loss(
            &mut issues,
            &selection.component,
            "Only the explicitly selected Sequence is projected; other sequence representations remain in the graph.",
        )?;
    }
    if !component.constraints.is_empty() {
        loss(
            &mut issues,
            &selection.component,
            "Design constraints are retained in the SBOL graph, not checked or promoted into native cloning verification.",
        )?;
    }
    if !component.roles.is_empty() {
        loss(
            &mut issues,
            &selection.component,
            "Component ontology roles remain in the graph; they are not executable biological or chemical claims.",
        )?;
    }
    let supported_types: &[&str] = match alphabet {
        Alphabet::Dna => &[DNA, LINEAR, CIRCULAR, "http://identifiers.org/SO:0000985"],
        Alphabet::Rna => &[RNA, LINEAR, CIRCULAR, "http://identifiers.org/SO:0000984"],
        Alphabet::Protein => &[PROTEIN],
    };
    if component
        .types
        .iter()
        .any(|t| !supported_types.contains(&semantic_so_type(t).unwrap_or(t.as_str())))
    {
        return Err(Error::Projection("additional/unsupported physical Component type cannot be discarded safely, even with AllowOmissions (including single-stranded DNA)".into()));
    }
    if declared.is_none() {
        warnings.push("Native topology follows the explicit caller choice; the source graph's unspecified topology remains unchanged.".into());
    }
    let sequence_model = Sequence::new(alphabet, text)?;
    if sequence_model.as_str() != text {
        warnings.push("IUPAC letter case was normalized to uppercase in the native reference; original SBOL bytes remain exact.".into());
    }
    let name = object
        .identity
        .name
        .as_ref()
        .or(object.identity.display_id.as_ref())
        .cloned()
        .unwrap_or_else(|| selection.component.clone());
    let mut molecule = Molecule {
        name,
        sequence: sequence_model,
        topology,
        features: vec![],
        metadata: BTreeMap::new(),
        original_genbank: None,
    };
    for fid in &component.features {
        match feature(doc, fid, &selection.sequence, alphabet, &mut issues) {
            Ok(f) => molecule.features.push(f),
            Err(e) => loss(&mut issues, fid, &format!("Feature omitted: {e}"))?,
        }
    }
    if policy == ProjectionPolicy::Strict && !issues.is_empty() {
        return Err(Error::Projection(format!(
            "{} explicit representation losses require AllowOmissions review; inspect retained graph",
            issues.len()
        )));
    }
    if let Some(text) = &object.identity.description {
        molecule
            .metadata
            .insert("sbol:description".into(), text.clone());
    }
    if let Some(text) = &object.identity.display_id {
        molecule
            .metadata
            .insert("sbol:display_id".into(), text.clone());
    }
    molecule
        .metadata
        .insert("sbol:component_iri".into(), selection.component.clone());
    molecule
        .metadata
        .insert("sbol:sequence_iri".into(), selection.sequence.clone());
    molecule
        .metadata
        .insert("sbol:source_hash".into(), doc.revision());
    molecule
        .metadata
        .insert("sbol:exchange_profile".into(), PROFILE.into());
    molecule.validate()?;
    collagene_model::bounds::model_preflight(&[&molecule])?;
    Ok(Projection {
        molecule,
        component: selection.component.clone(),
        sequence: selection.sequence.clone(),
        issues,
        warnings,
    })
}
