//! Atomic, source-bound edits of existing annotations. Sequence calls never change.
use collagene_model::{Alphabet, Error, Feature, Location, Molecule};
use serde::{Deserialize, Serialize};
use std::collections::BTreeSet;

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Request {
    pub source_revision: String,
    pub indices: Vec<usize>,
    pub kind: Option<String>,
    /// Reverse the biological location using one outer complement toggle.
    pub flip_strand: bool,
    pub color: ColorEdit,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub enum ColorEdit {
    Keep,
    Set(String),
    Clear,
}

fn bad(message: &str) -> Error {
    Error::Unsupported(format!("annotation batch: {message}"))
}

/// Only a literal six-digit hexadecimal RGB value is accepted. No CSS parsing,
/// named colors, alpha, whitespace, or executable/style text is interpreted.
pub fn parse_rgb(value: &str) -> Option<[u8; 3]> {
    let bytes = value.as_bytes();
    if bytes.len() != 7 || bytes[0] != b'#' {
        return None;
    }
    fn digit(byte: u8) -> Option<u8> {
        match byte {
            b'0'..=b'9' => Some(byte - b'0'),
            b'a'..=b'f' => Some(byte - b'a' + 10),
            b'A'..=b'F' => Some(byte - b'A' + 10),
            _ => None,
        }
    }
    Some([
        digit(bytes[1])? * 16 + digit(bytes[2])?,
        digit(bytes[3])? * 16 + digit(bytes[4])?,
        digit(bytes[5])? * 16 + digit(bytes[6])?,
    ])
}

fn color_key(key: &str) -> bool {
    ["color", "ApEinfo_fwdcolor", "ApEinfo_revcolor"]
        .iter()
        .any(|name| key.eq_ignore_ascii_case(name))
}

/// Bounds and validates the source and entire request before cloning once.
/// Unknown qualifier multiplicity, unselected features, fuzzy/compound locations,
/// original source evidence, and reference calls are preserved.
///
/// Strand reversal changes the annotation's biological meaning without changing
/// source calls. Coordinate, translation and phase assertions therefore become
/// bounded source-evidence notes and require review; their positions are never
/// reinterpreted as preserved biological calls. Type and color edits leave other
/// qualifiers intact.
pub fn apply(source: &Molecule, request: &Request) -> Result<Molecule, Error> {
    collagene_model::bounds::model_preflight(&[source])?;
    source.validate()?;
    if request.source_revision.len() != 64
        || !request
            .source_revision
            .bytes()
            .all(|b| b.is_ascii_hexdigit())
        || request.indices.is_empty()
        || request.indices.len() > source.features.len()
    {
        return Err(bad(
            "choose existing annotations on an exact source revision",
        ));
    }
    if request.kind.as_ref().is_some_and(|kind| {
        kind.is_empty()
            || kind.len() > 64
            || !kind
                .bytes()
                .all(|b| b.is_ascii_alphanumeric() || b"_'-".contains(&b))
    }) {
        return Err(bad(
            "feature type needs 1–64 letters, digits, underscores, apostrophes or hyphens",
        ));
    }
    let color = match &request.color {
        ColorEdit::Set(value) => {
            let [r, g, b] = parse_rgb(value).ok_or_else(|| bad("color must be #RRGGBB"))?;
            Some(format!("#{r:02x}{g:02x}{b:02x}"))
        }
        _ => None,
    };
    if request.flip_strand && source.sequence.alphabet() == Alphabet::Protein {
        return Err(bad("protein annotations do not have a nucleotide strand"));
    }
    let mut seen = BTreeSet::new();
    for &index in &request.indices {
        source
            .features
            .get(index)
            .ok_or_else(|| bad("annotation index is out of bounds"))?;
        if !seen.insert(index) {
            return Err(bad("duplicate annotation index"));
        }
    }
    if request.source_revision != source.revision() {
        return Err(bad("source revision changed"));
    }
    let mut result = source.clone();
    for &index in &request.indices {
        let feature = &mut result.features[index];
        if let Some(kind) = &request.kind {
            feature.kind.clone_from(kind);
        }
        if request.flip_strand {
            let old = std::mem::replace(&mut feature.location, Location::range(0, 0));
            feature.location = match old {
                Location::Complement(inner) => *inner,
                other => Location::Complement(Box::new(other)),
            };
            collagene_model::annotations::deactivate_assertions(feature)?;
        }
        if !matches!(request.color, ColorEdit::Keep) {
            feature.qualifiers.retain(|(key, _)| !color_key(key));
            if let Some(value) = &color {
                feature
                    .qualifiers
                    .push(("color".into(), Some(value.clone())));
            }
        }
    }
    collagene_model::bounds::model_preflight(&[&result])?;
    result.validate()?;
    Ok(result)
}

fn qualifier_rgb(feature: &Feature, name: &str) -> Result<Option<[u8; 3]>, ()> {
    if feature.qualifiers.len() > 16_384 {
        return Err(());
    }
    let mut found = None;
    for (key, value) in &feature.qualifiers {
        if !key.eq_ignore_ascii_case(name) {
            continue;
        }
        let Some(rgb) = value.as_deref().and_then(parse_rgb) else {
            continue;
        };
        if found.is_some_and(|previous| previous != rgb) {
            return Err(());
        }
        found = Some(rgb);
    }
    Ok(found)
}

/// Safe color for a rendered strand/segment. A valid shared `/color` takes
/// precedence over ApE directional colors. A missing directional color falls
/// back to the other direction; conflicting duplicate values return None.
pub fn feature_rgb_for_strand(feature: &Feature, reverse: bool) -> Option<[u8; 3]> {
    if let Some(shared) = qualifier_rgb(feature, "color").ok()? {
        return Some(shared);
    }
    let forward = qualifier_rgb(feature, "ApEinfo_fwdcolor").ok()?;
    let backward = qualifier_rgb(feature, "ApEinfo_revcolor").ok()?;
    if reverse {
        backward.or(forward)
    } else {
        forward.or(backward)
    }
}

/// Safe color for a whole feature. Mixed/unknown strand features with different
/// ApE directional colors have no single color; use `feature_rgb_for_strand`
/// per segment instead. Unsupported color text is ignored, never emitted as CSS.
pub fn feature_rgb(feature: &Feature) -> Option<[u8; 3]> {
    if let Some(shared) = qualifier_rgb(feature, "color").ok()? {
        return Some(shared);
    }
    let forward = feature_rgb_for_strand(feature, false);
    let backward = feature_rgb_for_strand(feature, true);
    if forward == backward {
        return forward;
    }
    let mut stack = vec![(&feature.location, false, 1usize)];
    let mut nodes = 0usize;
    let mut orientation = None;
    while let Some((location, reverse, depth)) = stack.pop() {
        nodes += 1;
        if depth > 32 || nodes > 16_384 {
            return None;
        }
        match location {
            Location::Range { .. } => {
                if orientation.is_some_and(|old| old != reverse) {
                    return None;
                }
                orientation = Some(reverse);
            }
            Location::Complement(inner) => stack.push((inner, !reverse, depth + 1)),
            Location::Join(parts) | Location::Order(parts) => {
                if parts.len() > 16_384usize.saturating_sub(nodes + stack.len()) {
                    return None;
                }
                stack.extend(parts.iter().map(|part| (part, reverse, depth + 1)));
            }
            Location::Unresolved(_) => return None,
        }
    }
    match orientation {
        Some(false) => forward,
        Some(true) => backward,
        None => None,
    }
}
