//! Validated annotation editing and explicit reference-document metadata changes.
use collagene_model::{Error, Feature, Molecule, Topology};

pub fn set_feature(
    molecule: &Molecule,
    index: Option<usize>,
    mut feature: Feature,
) -> Result<Molecule, Error> {
    // Bound both inputs before recursive validation, evidence copies or cloning.
    collagene_model::bounds::model_preflight(&[molecule])?;
    collagene_model::bounds::feature_preflight(&[std::slice::from_ref(&feature)])?;
    molecule.validate()?;
    if feature.kind.is_empty()
        || feature.kind.len() > 64
        || !feature
            .kind
            .bytes()
            .all(|b| b.is_ascii_alphanumeric() || b"_'-".contains(&b))
    {
        return Err(Error::Format(
            "feature type needs 1–64 letters, digits, underscores, apostrophes or hyphens".into(),
        ));
    }
    feature.location.validate(molecule.sequence.len())?;
    if feature.qualifiers.len() > 10_000
        || feature.qualifiers.iter().any(|(key, value)| {
            key.is_empty()
                || key.len() > 128
                || !key
                    .bytes()
                    .all(|b| b.is_ascii_alphanumeric() || b"_:-".contains(&b))
                || value.as_ref().is_some_and(|v| v.len() > 1_048_576)
        })
    {
        return Err(Error::Format(
            "invalid qualifier key or annotation size limit exceeded".into(),
        ));
    }
    if let Some(old) = index.and_then(|i| molecule.features.get(i)) {
        if old.location != feature.location {
            if old
                .qualifiers
                .iter()
                .chain(&feature.qualifiers)
                .any(|(key, _)| crate::sequence::annotations::coordinate_qualifier(key))
            {
                return Err(Error::Unsupported(
                    "Changing an annotation location cannot reinterpret coordinate-bearing qualifiers. Remove or review those qualifiers at the current location before relocating the annotation, then explicitly curate the new annotation.".into(),
                ));
            }
            let mut evidence = Vec::new();
            feature.qualifiers.retain(|(key, value)| {
                if key.eq_ignore_ascii_case("translation") || key.eq_ignore_ascii_case("codon_start") {
                    evidence.push(("note".into(), Some(format!(
                        "Collagene source qualifier /{key}={} retained only as source evidence; recalculate after annotation location review.",
                        value.as_deref().unwrap_or("(flag)")))));
                    false
                } else { true }
            });
            feature.qualifiers.extend(evidence);
            feature.needs_review = true;
        }
    }
    collagene_model::bounds::feature_preflight(&[std::slice::from_ref(&feature)])?;
    let mut result = molecule.clone();
    match index {
        Some(index) => {
            *result
                .features
                .get_mut(index)
                .ok_or_else(|| Error::Coordinates("annotation no longer exists".into()))? = feature
        }
        None => result.features.push(feature),
    }
    collagene_model::bounds::model_preflight(&[&result])?;
    result.validate()?;
    Ok(result)
}
pub fn delete_feature(molecule: &Molecule, index: usize) -> Result<Molecule, Error> {
    molecule.validate()?;
    if index >= molecule.features.len() {
        return Err(Error::Coordinates("annotation no longer exists".into()));
    }
    let mut result = molecule.clone();
    result.features.remove(index);
    Ok(result)
}
pub fn rename(molecule: &Molecule, name: &str) -> Result<Molecule, Error> {
    molecule.validate()?;
    if name.trim().is_empty() || name.len() > 1024 || name.chars().any(char::is_control) {
        return Err(Error::Format(
            "name must be 1–1024 bytes without control characters".into(),
        ));
    }
    let mut result = molecule.clone();
    result.name = name.trim().into();
    Ok(result)
}
/// Reference topology change; does not assert a physical cleavage or covalent ligation.
pub fn set_reference_topology(molecule: &Molecule, topology: Topology) -> Result<Molecule, Error> {
    crate::sequence::annotations::change_reference_topology(molecule, topology)
}
