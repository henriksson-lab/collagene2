//! Independently declared internal substitutions in original source contributions.
//! No primer, preparation or simulated product supplies this reference target.
use crate::{
    Alphabet, Error, Feature, Location, Molecule, Sequence, Topology,
    assembly_intent::{self, PlacedFragment},
    assembly_junction_edit, assembly_recipe,
    verification::PlanInput,
};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};
pub const ENGINE_VERSION: &str = "assembly-substitution-intent-1";
pub const MAX_SUBSTITUTIONS: usize = 256;
const MAX_BASES: usize = 8 << 20;
const PREFIX: &str = "assembly_substitution_intent:";
macro_rules! dto{($($t:tt)*)=>{#[derive(Clone,Debug,Serialize,Deserialize)]#[cfg_attr(feature="json-schema",derive(schemars::JsonSchema))]#[serde(deny_unknown_fields)]$($t)*};}
dto! {#[derive(PartialEq,Eq)] pub struct Substitution{pub role:String,pub offset:usize,pub expected:String,pub replacement:String}}
dto! {pub struct Request{pub base:assembly_intent::Request,pub junction_edits:Vec<assembly_junction_edit::Edit>,pub substitutions:Vec<Substitution>}}
dto! {#[derive(PartialEq,Eq)] pub struct AppliedSubstitution{pub role:String,pub input_id:String,pub input_revision:String,pub offset:usize,pub source_position:usize,pub output_position:usize,pub expected:String,pub replacement:String,pub reverse:bool}}
dto! {pub struct Report{pub context_revision:String,pub request:Request,pub target:Molecule,pub fragments:Vec<PlacedFragment>,pub applied:Vec<AppliedSubstitution>,pub warnings:Vec<String>}}
dto! {struct Pin{input_id:String,revision:String,length:usize,topology:Topology}}
fn bad(s: &str) -> Error {
    Error::Unsupported(format!("assembly substitution intent: {s}"))
}
fn hex(s: &str) -> bool {
    s.len() == 64
        && s.bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}
fn token(s: &str) -> bool {
    !s.trim().is_empty() && s.len() <= 128 && !s.chars().any(char::is_control)
}
fn base(s: &str) -> bool {
    s.len() == 1 && b"ACGT".contains(&s.as_bytes()[0])
}
fn shape(r: &Request) -> Result<(), Error> {
    if r.substitutions.len() > MAX_SUBSTITUTIONS
        || r.base.fragments.is_empty()
        || r.base.fragments.len() > 32
        || r.junction_edits.len() > 32
    {
        return Err(Error::TooLarge);
    }
    if r.base.name.trim().is_empty()
        || r.base.name.len() > 1024
        || r.base.name.chars().any(char::is_control)
    {
        return Err(bad("invalid intended construct name"));
    }
    let mut roles = BTreeSet::new();
    for f in &r.base.fragments {
        if !token(&f.role) || !token(&f.input_id) || !roles.insert(&f.role) {
            return Err(bad("invalid or repeated contribution identity"));
        }
    }
    let mut seen = BTreeSet::new();
    for edit in &r.substitutions {
        if !token(&edit.role)
            || !base(&edit.expected)
            || !base(&edit.replacement)
            || edit.expected == edit.replacement
            || !seen.insert((edit.role.as_str(), edit.offset))
        {
            return Err(bad(
                "declare unique role/offset substitutions with distinct single uppercase concrete DNA bases",
            ));
        }
    }
    Ok(())
}
fn pins(inputs: &[PlanInput]) -> Vec<Pin> {
    inputs
        .iter()
        .map(|i| Pin {
            input_id: i.id.clone(),
            revision: i.revision.clone(),
            length: i.molecule.sequence.len(),
            topology: i.molecule.topology,
        })
        .collect()
}
fn pin_shape(pins: &[Pin]) -> Result<(), Error> {
    if pins.is_empty() || pins.len() > 32 {
        return Err(Error::TooLarge);
    }
    let mut seen = BTreeSet::new();
    let mut total = 0usize;
    for p in pins {
        if !token(&p.input_id)
            || !hex(&p.revision)
            || p.length == 0
            || p.length > MAX_BASES
            || !seen.insert(&p.input_id)
        {
            return Err(bad("invalid original source pins"));
        }
        total = total
            .checked_add(p.length)
            .filter(|n| *n <= MAX_BASES)
            .ok_or(Error::TooLarge)?;
    }
    Ok(())
}
fn context(pins: &[Pin], r: &Request) -> Result<String, Error> {
    assembly_recipe::hash(&(
        ENGINE_VERSION,
        assembly_junction_edit::ENGINE_VERSION,
        pins,
        r,
    ))
}
fn json(value: &impl Serialize) -> Result<String, Error> {
    let raw = serde_json::to_string(value).map_err(|_| Error::TooLarge)?;
    if raw.len() > 256 * 1024 {
        return Err(Error::TooLarge);
    }
    Ok(raw)
}
/// Derive source/output correspondence only. Original offsets are in construct
/// orientation BEFORE trims; outgoing inserted junction bases are not editable.
fn geometry(
    pins: &[Pin],
    r: &Request,
) -> Result<(Vec<PlacedFragment>, Vec<AppliedSubstitution>), Error> {
    shape(r)?;
    pin_shape(pins)?;
    let lookup = pins
        .iter()
        .map(|p| (p.input_id.as_str(), p))
        .collect::<BTreeMap<_, _>>();
    let lengths = r
        .base
        .fragments
        .iter()
        .map(|f| {
            lookup
                .get(f.input_id.as_str())
                .map(|p| p.length)
                .ok_or_else(|| bad("contribution has no original input"))
        })
        .collect::<Result<Vec<_>, _>>()?;
    let adjusted = assembly_junction_edit::intent_adjustment(&r.base, &r.junction_edits, &lengths)?;
    let mut inserts = vec![0usize; r.base.fragments.len()];
    for edit in &r.junction_edits {
        inserts[edit.junction] = if edit.inserted.bytes().all(|b| b.is_ascii_whitespace()) {
            0
        } else {
            Sequence::dna(&edit.inserted)?.len()
        };
    }
    let mut fragments = Vec::with_capacity(adjusted.fragments.len());
    let mut offset = 0usize;
    for (f, extra) in adjusted.fragments.iter().zip(&inserts) {
        let p = lookup[f.input_id.as_str()];
        if p.topology == Topology::Linear && f.end > p.length {
            return Err(bad("linear contribution exceeds its actual source"));
        }
        let end = offset
            .checked_add(f.end - f.start)
            .and_then(|n| n.checked_add(*extra))
            .filter(|n| *n <= MAX_BASES)
            .ok_or(Error::TooLarge)?;
        fragments.push(PlacedFragment {
            role: f.role.clone(),
            input_id: f.input_id.clone(),
            input_revision: p.revision.clone(),
            start: f.start,
            end: f.end,
            reverse: f.reverse,
            output_start: offset,
            output_end: end,
        });
        offset = end;
    }
    let roles = r
        .base
        .fragments
        .iter()
        .enumerate()
        .map(|(i, f)| (f.role.as_str(), i))
        .collect::<BTreeMap<_, _>>();
    let mut applied = Vec::with_capacity(r.substitutions.len());
    let mut outputs = BTreeSet::new();
    for edit in &r.substitutions {
        let i = *roles
            .get(edit.role.as_str())
            .ok_or_else(|| bad("substitution refers to an unknown contribution role"))?;
        let original = &r.base.fragments[i];
        let f = &fragments[i];
        let p = lookup[original.input_id.as_str()];
        let length = original.end - original.start;
        if edit.offset >= length {
            return Err(bad(
                "substitution lies outside its original selected contribution",
            ));
        }
        // The shared junction builder normalizes a retained circular core back
        // into the first source turn. Restore the ORIGINAL unrolled coordinate
        // frame before measuring trims; contribution offsets never rotate.
        let turn = if p.topology == Topology::Circular && f.start < original.start {
            p.length
        } else {
            0
        };
        let retained_start = f.start.checked_add(turn).ok_or(Error::TooLarge)?;
        let retained_end = f.end.checked_add(turn).ok_or(Error::TooLarge)?;
        let retained_begin = if original.reverse {
            original.end.checked_sub(retained_end)
        } else {
            retained_start.checked_sub(original.start)
        }
        .ok_or_else(|| bad("adjusted source geometry changed direction"))?;
        let local = edit
            .offset
            .checked_sub(retained_begin)
            .filter(|n| *n < f.end - f.start)
            .ok_or_else(|| bad("substitution is removed by a declared junction trim"))?;
        let output_position = f.output_start.checked_add(local).ok_or(Error::TooLarge)?;
        if !outputs.insert(output_position) {
            return Err(bad("substitutions overlap in the intended reference"));
        }
        let source_offset = if original.reverse {
            length - 1 - edit.offset
        } else {
            edit.offset
        };
        let source_position = original
            .start
            .checked_add(source_offset)
            .ok_or(Error::TooLarge)?
            % p.length;
        applied.push(AppliedSubstitution {
            role: edit.role.clone(),
            input_id: original.input_id.clone(),
            input_revision: p.revision.clone(),
            offset: edit.offset,
            source_position,
            output_position,
            expected: edit.expected.clone(),
            replacement: edit.replacement.clone(),
            reverse: original.reverse,
        });
    }
    Ok((fragments, applied))
}
/// Borrowed original inputs only. Source reference/physical markers do not turn
/// declared substitutions into new source material or physical ancestry.
pub fn source_preflight(inputs: &[(&str, &Molecule)], r: &Request) -> Result<(), Error> {
    shape(r)?;
    if inputs.is_empty() || inputs.len() > 32 {
        return Err(Error::TooLarge);
    }
    assembly_junction_edit::intent_preflight(inputs, &r.base, &r.junction_edits)?;
    let pins = inputs
        .iter()
        .map(|(id, m)| Pin {
            input_id: (*id).into(),
            revision: m.revision(),
            length: m.sequence.len(),
            topology: m.topology,
        })
        .collect::<Vec<_>>();
    let (_, applied) = geometry(&pins, r)?;
    let lookup = inputs.iter().copied().collect::<BTreeMap<_, _>>();
    for a in applied {
        let m = lookup[a.input_id.as_str()];
        let mut expected = m.sequence.as_str().as_bytes()[a.source_position];
        if a.reverse {
            expected = crate::model::complement(expected, Alphabet::Dna);
        }
        if expected != a.expected.as_bytes()[0] {
            return Err(bad(
                "expected source base is stale or ambiguous in construct orientation",
            ));
        }
    }
    assembly_recipe::hash(&(inputs, r))?;
    Ok(())
}
pub fn intent_revision(inputs: &[PlanInput], r: &Request) -> Result<String, Error> {
    if inputs.is_empty() || inputs.len() > 32 {
        return Err(Error::TooLarge);
    }
    source_preflight(
        &inputs
            .iter()
            .map(|i| (i.id.as_str(), &i.molecule))
            .collect::<Vec<_>>(),
        r,
    )?;
    for i in inputs {
        if !hex(&i.revision) || i.revision != i.molecule.revision() {
            return Err(bad("original source revision changed"));
        }
    }
    context(&pins(inputs), r)
}
fn affected(location: &Location, positions: &BTreeSet<usize>) -> Result<bool, Error> {
    match location {
        Location::Range {
            start,
            end,
            before,
            after,
        } => Ok(*before || *after || positions.range(*start..*end).next().is_some()),
        Location::Complement(inner) => affected(inner, positions),
        Location::Join(parts) | Location::Order(parts) => {
            let mut touched = false;
            for p in parts {
                touched |= affected(p, positions)?;
            }
            Ok(touched)
        }
        Location::Unresolved(_) => Err(bad(
            "unresolved projected annotation locations prevent safe substitution review; original annotations remain retained",
        )),
    }
}
pub fn build(inputs: &[PlanInput], r: &Request) -> Result<Report, Error> {
    let context_revision = intent_revision(inputs, r)?;
    let source_pins = pins(inputs);
    let (fragments, applied) = geometry(&source_pins, r)?;
    let initial = assembly_junction_edit::build_intent(inputs, &r.base, &r.junction_edits)?;
    if assembly_recipe::hash(&fragments)? != assembly_recipe::hash(&initial.fragments)? {
        return Err(bad(
            "source/junction correspondence differs from native intent",
        ));
    }
    let mut target = initial.target;
    let mut warnings = initial.warnings;
    let positions = applied
        .iter()
        .map(|a| a.output_position)
        .collect::<BTreeSet<_>>();
    if !positions.is_empty() {
        for feature in &mut target.features {
            if affected(&feature.location, &positions)? {
                if feature.qualifiers.iter().any(|(k, _)| {
                    ["transl_except", "anticodon", "rpt_unit_range", "coded_by"]
                        .iter()
                        .any(|key| k.eq_ignore_ascii_case(key))
                }) {
                    return Err(bad(
                        "affected coordinate-bearing annotation qualifiers need explicit review before substitution; original sources remain intact",
                    ));
                }
                feature.needs_review = true;
                feature
                    .qualifiers
                    .retain(|(key, _)| !key.eq_ignore_ascii_case("translation"));
                feature.qualifiers.push(("note".into(),Some("Operator-declared reference substitution intersects this annotation; recalculate sequence/function assertions. Complete original qualifiers remain in the original source input.".into())));
            }
        }
    }
    let mut bases = target.sequence.as_str().as_bytes().to_vec();
    for edit in &applied {
        if bases[edit.output_position] != edit.expected.as_bytes()[0] {
            return Err(bad(
                "expected base differs after native source/junction construction",
            ));
        }
        bases[edit.output_position] = edit.replacement.as_bytes()[0];
    }
    target.sequence = Sequence::dna(std::str::from_utf8(&bases).map_err(|_| Error::TooLarge)?)?;
    if target
        .features
        .len()
        .checked_add(applied.len())
        .is_none_or(|n| n > 4096)
    {
        return Err(Error::TooLarge);
    }
    for edit in &applied {
        target.features.push(Feature{kind:"misc_feature".into(),location:Location::range(edit.output_position,edit.output_position+1),qualifiers:vec![("label".into(),Some(format!("Declared {}{}{}",edit.expected,edit.offset+1,edit.replacement))),("note".into(),Some(format!("Independent operator-declared substitution in role {}, original construct-oriented offset {}. This is reference intent, not an inherited base or physical synthesis claim.",edit.role,edit.offset)))],needs_review:true});
    }
    target.metadata.retain(|key, _| {
        !key.starts_with("assembly_intent:") && !key.starts_with("assembly_junction_edit:")
    });
    for (key, value) in [
        ("engine", ENGINE_VERSION.into()),
        ("context_revision", context_revision.clone()),
        ("request", json(r)?),
        ("inputs", json(&source_pins)?),
        ("applied", json(&applied)?),
        ("contributions", json(&fragments)?),
        (
            "sequence_blake3",
            blake3::hash(target.sequence.as_str().as_bytes())
                .to_hex()
                .to_string(),
        ),
        (
            "features_revision",
            assembly_recipe::hash(&target.features)?,
        ),
    ] {
        target.metadata.insert(format!("{PREFIX}{key}"), value);
    }
    target.metadata.insert(
        "duplex_model".into(),
        "Reference-only assembly substitution design intent; physical material state unavailable"
            .into(),
    );
    target
        .metadata
        .insert("assembly:intent:reference_only".into(), "true".into());
    warnings.push("Internal substitutions are explicitly declared in original contribution orientation before junction trims. Original sources/preparations remain unchanged; no PCR result or simulated product determines this target.".into());
    warnings.push("Affected annotations require review and stale translations are removed from the intended reference. Original features and qualifiers remain retained in exact original input documents; this reference asserts no inherited chemical state, biological function or experimental success.".into());
    let report = Report {
        context_revision: context_revision.clone(),
        request: r.clone(),
        target,
        fragments,
        applied,
        warnings,
    };
    report.validate_binding(&context_revision)?;
    Ok(report)
}
impl Report {
    /// Bounded readonly consistency, never an independent rebuild or scientific proof.
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        shape(&self.request)?;
        crate::model::bounds::model_preflight(&[&self.target])?;
        self.target.validate()?;
        if !hex(expected)
            || self.context_revision != expected
            || self.target.sequence.alphabet() != Alphabet::Dna
            || self.target.sequence.len() > MAX_BASES
            || self.target.name != self.request.base.name
            || self.target.topology
                != if self.request.base.circular {
                    Topology::Circular
                } else {
                    Topology::Linear
                }
            || self.warnings.len() > 2048
            || self.warnings.iter().any(|s| s.len() > 16384)
        {
            return Err(bad("invalid intent target/context/bounds"));
        }
        let get = |key: &str| {
            self.target
                .metadata
                .get(&format!("{PREFIX}{key}"))
                .filter(|s| s.len() <= 256 * 1024)
                .ok_or_else(|| bad("missing bounded intent provenance"))
        };
        let source_pins: Vec<Pin> = serde_json::from_str(get("inputs")?)
            .map_err(|_| bad("invalid original source pins"))?;
        let (fragments, applied) = geometry(&source_pins, &self.request)?;
        if context(&source_pins, &self.request)? != self.context_revision
            || assembly_recipe::hash(&fragments)? != assembly_recipe::hash(&self.fragments)?
            || applied != self.applied
            || get("engine")? != ENGINE_VERSION
            || get("context_revision")? != &self.context_revision
            || get("request")? != &json(&self.request)?
            || get("applied")? != &json(&self.applied)?
            || get("contributions")? != &json(&self.fragments)?
            || get("sequence_blake3")?
                != blake3::hash(self.target.sequence.as_str().as_bytes())
                    .to_hex()
                    .as_str()
            || get("features_revision")? != &assembly_recipe::hash(&self.target.features)?
            || self
                .target
                .metadata
                .get("assembly:intent:reference_only")
                .map(String::as_str)
                != Some("true")
            || self.target.metadata.get("duplex_model").map(String::as_str)
                != Some(
                    "Reference-only assembly substitution design intent; physical material state unavailable",
                )
        {
            return Err(bad(
                "intent provenance/source/request correspondence changed",
            ));
        }
        let length = self.fragments.last().map_or(0, |f| f.output_end);
        if length != self.target.sequence.len() {
            return Err(bad("intended reference length changed"));
        }
        for edit in &self.applied {
            if self
                .target
                .sequence
                .as_str()
                .as_bytes()
                .get(edit.output_position)
                != Some(&edit.replacement.as_bytes()[0])
            {
                return Err(bad("applied substitution differs from target"));
            }
        }
        assembly_recipe::hash(self)?;
        Ok(())
    }
}
