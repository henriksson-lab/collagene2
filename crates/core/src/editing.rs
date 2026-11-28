//! Validated annotation editing and explicit reference-document metadata changes.
use collagene_model::{Error, Feature, Location, Molecule, Topology};

pub fn set_feature(
    molecule: &Molecule,
    index: Option<usize>,
    mut feature: Feature,
) -> Result<Molecule, Error> {
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
    if index
        .and_then(|i| molecule.features.get(i))
        .is_some_and(|old| old.location != feature.location)
    {
        feature.qualifiers.retain(|(key, _)| key != "translation");
        feature.needs_review = true;
    }
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
    crate::annotations::change_reference_topology(molecule, topology)
}

pub fn location_text(location: &Location) -> Result<String, Error> {
    crate::formats::location_text(location)
}
pub fn parse_location(text: &str, length: usize) -> Result<Location, Error> {
    crate::formats::parse_location(text, length)
}
