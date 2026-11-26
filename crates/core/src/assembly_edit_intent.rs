//! Literal internal contribution edits declared independently of preparation or products.
//! Offsets refer to original construct-oriented contributions, before junction trims.
use crate::{
    Alphabet, Error, Molecule, Topology, assembly_intent::PlacedFragment, assembly_junction_edit,
    assembly_recipe, assembly_substitution_intent, verification::PlanInput,
};
use serde::{Deserialize, Serialize};
use std::{collections::BTreeMap, io::Write};
pub const ENGINE_VERSION: &str = "assembly-edit-intent-1";
pub const MAX_EDITS: usize = 256;
const MAX_EDIT_BASES: usize = 100_000;
const MAX_BASES: usize = 8 << 20;
const MAX_COPY_WORK: usize = 128 << 20;
const PREFIX: &str = "assembly_edit_intent:";
macro_rules! dto {($($t:tt)*)=>{#[derive(Clone,Debug,Serialize,Deserialize)]#[cfg_attr(feature="json-schema",derive(schemars::JsonSchema))]#[serde(deny_unknown_fields)]$($t)*};}
dto! {#[derive(PartialEq,Eq)] pub struct Edit {pub role:String,pub start:usize,pub end:usize,pub expected:String,pub replacement:String}}
dto! {pub struct Request {pub base:assembly_substitution_intent::Request,pub edits:Vec<Edit>}}
dto! {#[derive(PartialEq,Eq)] pub struct AppliedEdit {pub role:String,pub input_id:String,pub input_revision:String,pub start:usize,pub end:usize,pub reverse:bool,pub original_output_start:usize,pub original_output_end:usize,pub output_start:usize,pub output_end:usize,pub expected:String,pub replacement:String}}
dto! {pub struct Report {pub context_revision:String,pub request:Request,pub target:Molecule,pub fragments:Vec<PlacedFragment>,pub applied:Vec<AppliedEdit>,pub warnings:Vec<String>}}
dto! {struct Pin {input_id:String,revision:String,length:usize,topology:Topology}}
fn pins(inputs: &[(&str, &Molecule)]) -> Vec<Pin> {
    inputs
        .iter()
        .map(|(id, m)| Pin {
            input_id: (*id).into(),
            revision: m.revision(),
            length: m.sequence.len(),
            topology: m.topology,
        })
        .collect()
}
fn bad(message: &str) -> Error {
    Error::Unsupported(format!("assembly edit intent: {message}"))
}
fn literal(s: &str) -> bool {
    s.bytes().all(|b| b"ACGTRYSWKMBDHVN".contains(&b))
}
fn conflicts(a: usize, b: usize, c: usize, d: usize) -> bool {
    if a == b {
        c <= a && a <= d
    } else if c == d {
        a <= c && c <= b
    } else {
        a < d && c < b
    }
}
fn shape(r: &Request) -> Result<(), Error> {
    if r.edits.len() > MAX_EDITS
        || r.base.substitutions.len() > 256
        || r.base.base.fragments.len() > 32
        || r.base.junction_edits.len() > 32
    {
        return Err(Error::TooLarge);
    }
    let name = &r.base.base.name;
    if name.trim().is_empty() || name.len() > 1024 || name.chars().any(char::is_control) {
        return Err(bad("invalid bounded intended name"));
    }
    let mut roles = std::collections::BTreeSet::new();
    for f in &r.base.base.fragments {
        if f.role.trim().is_empty()
            || f.role.len() > 128
            || f.role.chars().any(char::is_control)
            || f.input_id.trim().is_empty()
            || f.input_id.len() > 128
            || f.input_id.chars().any(char::is_control)
            || !roles.insert(f.role.as_str())
        {
            return Err(bad("invalid or duplicate source role"));
        }
    }
    let mut positions = std::collections::BTreeSet::new();
    for e in &r.base.substitutions {
        if !roles.contains(e.role.as_str())
            || e.expected.len() != 1
            || e.replacement.len() != 1
            || !e.expected.bytes().all(|b| b"ACGT".contains(&b))
            || !e.replacement.bytes().all(|b| b"ACGT".contains(&b))
            || e.expected == e.replacement
            || !positions.insert((e.role.as_str(), e.offset))
        {
            return Err(bad("invalid or duplicate substitution declaration"));
        }
        let f = r
            .base
            .base
            .fragments
            .iter()
            .find(|f| f.role == e.role)
            .ok_or_else(|| bad("unknown role"))?;
        if f.end.checked_sub(f.start).is_none_or(|n| e.offset >= n) {
            return Err(bad(
                "substitution lies outside original source contribution",
            ));
        }
    }
    let mut total = 0usize;
    for (i, e) in r.edits.iter().enumerate() {
        if e.role.is_empty()
            || e.role.len() > 128
            || e.role.chars().any(char::is_control)
            || e.end.checked_sub(e.start) != Some(e.expected.len())
            || !literal(&e.expected)
            || !literal(&e.replacement)
            || e.expected == e.replacement
        {
            return Err(bad(
                "declare a bounded role and distinct uppercase DNA literals matching the original half-open interval",
            ));
        }
        total = total
            .checked_add(e.expected.len())
            .and_then(|n| n.checked_add(e.replacement.len()))
            .filter(|n| *n <= MAX_EDIT_BASES)
            .ok_or(Error::TooLarge)?;
        if r.edits[..i]
            .iter()
            .any(|p| p.role == e.role && conflicts(e.start, e.end, p.start, p.end))
        {
            return Err(bad(
                "internal edits overlap or share an ambiguous insertion boundary",
            ));
        }
        if r.base.substitutions.iter().any(|p| {
            p.role == e.role && conflicts(e.start, e.end, p.offset, p.offset.saturating_add(1))
        }) {
            return Err(bad("internal edit conflicts with a declared substitution"));
        }
    }
    Ok(())
}
fn geometry(inputs: &[Pin], r: &Request) -> Result<Vec<AppliedEdit>, Error> {
    let sources = inputs
        .iter()
        .map(|p| (p.input_id.as_str(), p))
        .collect::<BTreeMap<_, _>>();
    let originals = &r.base.base.fragments;
    let lengths = originals
        .iter()
        .map(|f| {
            sources
                .get(f.input_id.as_str())
                .map(|m| m.length)
                .ok_or_else(|| bad("unknown original input"))
        })
        .collect::<Result<Vec<_>, _>>()?;
    let adjusted =
        assembly_junction_edit::intent_adjustment(&r.base.base, &r.base.junction_edits, &lengths)?;
    let mut offsets = Vec::with_capacity(originals.len());
    let mut output = 0usize;
    for (i, f) in adjusted.fragments.iter().enumerate() {
        offsets.push(output);
        output = output.checked_add(f.end - f.start).ok_or(Error::TooLarge)?;
        if let Some(j) = r.base.junction_edits.iter().find(|e| e.junction == i) {
            output = output
                .checked_add(
                    j.inserted
                        .bytes()
                        .filter(|b| !b.is_ascii_whitespace())
                        .count(),
                )
                .ok_or(Error::TooLarge)?;
        }
    }
    let mut out = Vec::with_capacity(r.edits.len());
    for e in &r.edits {
        let i = originals
            .iter()
            .position(|f| f.role == e.role)
            .ok_or_else(|| bad("edit refers to an unknown contribution role"))?;
        let original = &originals[i];
        let retained = &adjusted.fragments[i];
        let source = sources[original.input_id.as_str()];
        let turn = if source.topology == Topology::Circular && retained.start < original.start {
            source.length
        } else {
            0
        };
        let begin = if original.reverse {
            original
                .end
                .checked_sub(retained.end.checked_add(turn).ok_or(Error::TooLarge)?)
        } else {
            retained
                .start
                .checked_add(turn)
                .and_then(|n| n.checked_sub(original.start))
        }
        .ok_or_else(|| bad("trimmed contribution changed coordinate direction"))?;
        let a = e
            .start
            .checked_sub(begin)
            .ok_or_else(|| bad("edit is removed by a junction trim"))?;
        let b = e
            .end
            .checked_sub(begin)
            .filter(|n| *n <= retained.end - retained.start)
            .ok_or_else(|| bad("edit extends beyond the retained contribution"))?;
        let start = offsets[i].checked_add(a).ok_or(Error::TooLarge)?;
        let end = offsets[i].checked_add(b).ok_or(Error::TooLarge)?;
        out.push(AppliedEdit {
            role: e.role.clone(),
            input_id: original.input_id.clone(),
            input_revision: source.revision.clone(),
            start: e.start,
            end: e.end,
            reverse: original.reverse,
            original_output_start: start,
            original_output_end: end,
            output_start: 0,
            output_end: 0,
            expected: e.expected.clone(),
            replacement: e.replacement.clone(),
        });
    }
    for (i, e) in out.iter().enumerate() {
        if out[..i].iter().any(|p| {
            conflicts(
                e.original_output_start,
                e.original_output_end,
                p.original_output_start,
                p.original_output_end,
            )
        }) {
            return Err(bad("edits share an ambiguous output boundary"));
        }
    }
    // Preserve request order in evidence while computing final coordinates in output order.
    let mut order = (0..out.len()).collect::<Vec<_>>();
    order.sort_by_key(|i| out[*i].original_output_start);
    let mut removed = 0usize;
    let mut added = 0usize;
    for i in order {
        let e = &mut out[i];
        e.output_start = e
            .original_output_start
            .checked_sub(removed)
            .and_then(|n| n.checked_add(added))
            .ok_or(Error::TooLarge)?;
        e.output_end = e
            .output_start
            .checked_add(e.replacement.len())
            .ok_or(Error::TooLarge)?;
        removed = removed
            .checked_add(e.end - e.start)
            .ok_or(Error::TooLarge)?;
        added = added
            .checked_add(e.replacement.len())
            .ok_or(Error::TooLarge)?;
    }
    if output.checked_add(added).is_none_or(|n| n > MAX_BASES)
        || output
            .checked_sub(removed)
            .and_then(|n| n.checked_add(added))
            .is_none_or(|n| n == 0)
    {
        return Err(Error::TooLarge);
    }
    Ok(out)
}
/// Borrowed bounded validation before constructing the intended reference.
pub fn source_preflight(inputs: &[(&str, &Molecule)], r: &Request) -> Result<(), Error> {
    shape(r)?;
    assembly_substitution_intent::source_preflight(inputs, &r.base)?;
    geometry(&pins(inputs), r)?;
    let sources = inputs.iter().copied().collect::<BTreeMap<_, _>>();
    for e in &r.edits {
        let original = r
            .base
            .base
            .fragments
            .iter()
            .find(|f| f.role == e.role)
            .ok_or_else(|| bad("unknown role"))?;
        let source = sources[original.input_id.as_str()];
        // Check literal against immutable original letters, never primers or a mutated target.
        for (offset, expected) in (e.start..e.end).zip(e.expected.bytes()) {
            let local = if original.reverse {
                original.end - original.start - 1 - offset
            } else {
                offset
            };
            let pos =
                original.start.checked_add(local).ok_or(Error::TooLarge)? % source.sequence.len();
            let raw = source.sequence.as_str().as_bytes()[pos];
            let actual = if original.reverse {
                crate::model::complement(raw, Alphabet::Dna)
            } else {
                raw
            };
            if actual != expected {
                return Err(bad(
                    "expected literal differs from the original construct-oriented source",
                ));
            }
        }
    }

    Ok(())
}
pub fn intent_revision(inputs: &[PlanInput], r: &Request) -> Result<String, Error> {
    source_preflight(
        &inputs
            .iter()
            .map(|i| (i.id.as_str(), &i.molecule))
            .collect::<Vec<_>>(),
        r,
    )?;
    let base = assembly_substitution_intent::intent_revision(inputs, &r.base)?;
    assembly_recipe::hash(&(ENGINE_VERSION, base, r))
}
struct Count(usize);
impl Write for Count {
    fn write(&mut self, b: &[u8]) -> std::io::Result<usize> {
        self.0 = self
            .0
            .checked_add(b.len())
            .filter(|n| *n <= MAX_COPY_WORK)
            .ok_or_else(|| std::io::Error::other("intended-edit copy budget"))?;
        Ok(b.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}
pub fn build(inputs: &[PlanInput], r: &Request) -> Result<Report, Error> {
    let context_revision = intent_revision(inputs, r)?;
    let source_pins = pins(
        &inputs
            .iter()
            .map(|i| (i.id.as_str(), &i.molecule))
            .collect::<Vec<_>>(),
    );
    let mut applied = geometry(&source_pins, r)?;
    let initial = assembly_substitution_intent::build(inputs, &r.base)?;
    let mut target = initial.target;
    let mut count = Count(0);
    serde_json::to_writer(&mut count, &target).map_err(|_| Error::TooLarge)?;
    if count
        .0
        .checked_mul(applied.len())
        .is_none_or(|n| n > MAX_COPY_WORK)
    {
        return Err(bad(
            "bounded annotation/sequence copy budget exceeded; split this design into smaller explicit edit sets",
        ));
    }
    // Remove old target-coordinate receipts; the full original declaration is embedded below.
    target
        .metadata
        .retain(|k, _| !k.starts_with("assembly_substitution_intent:"));
    let mut order = (0..applied.len()).collect::<Vec<_>>();
    order.sort_by_key(|i| std::cmp::Reverse(applied[*i].original_output_start));
    for i in order {
        let e = &applied[i];
        if target
            .sequence
            .as_str()
            .get(e.original_output_start..e.original_output_end)
            != Some(e.expected.as_str())
        {
            return Err(bad(
                "native source/junction reference differs from expected edit literal",
            ));
        }
        target = target.replace(
            e.original_output_start,
            e.original_output_end,
            &e.replacement,
        )?;
    }
    let mut fragments = initial.fragments;
    let mut offset = 0usize;
    for f in &mut fragments {
        let mut len = f.output_end - f.output_start;
        for e in applied.iter().filter(|e| e.role == f.role) {
            len = len
                .checked_sub(e.end - e.start)
                .and_then(|n| n.checked_add(e.replacement.len()))
                .ok_or(Error::TooLarge)?;
        }
        f.output_start = offset;
        offset = offset.checked_add(len).ok_or(Error::TooLarge)?;
        f.output_end = offset;
    }
    if offset != target.sequence.len() {
        return Err(bad(
            "edited contribution accounting differs from intended DNA",
        ));
    }
    for (key, value) in [
        ("engine", ENGINE_VERSION.to_owned()),
        (
            "inputs",
            serde_json::to_string(&source_pins).map_err(|_| Error::TooLarge)?,
        ),
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
        ("context_revision", context_revision.clone()),
        (
            "request",
            serde_json::to_string(r).map_err(|_| Error::TooLarge)?,
        ),
        (
            "contributions",
            serde_json::to_string(&fragments).map_err(|_| Error::TooLarge)?,
        ),
        (
            "applied",
            serde_json::to_string(&applied).map_err(|_| Error::TooLarge)?,
        ),
    ] {
        if value.len() > 1024 * 1024 {
            return Err(Error::TooLarge);
        }
        target.metadata.insert(format!("{PREFIX}{key}"), value);
    }
    target.metadata.insert(
        "duplex_model".into(),
        "Reference-only assembly edit design intent; physical material state unavailable".into(),
    );
    target
        .metadata
        .insert("assembly:intent:reference_only".into(), "true".into());
    crate::model::bounds::model_preflight(&[&target])?;
    target.validate()?;
    let mut warnings = initial.warnings;
    warnings.push("Internal insertions, deletions and replacements are literal operator declarations in original contribution orientation before junction trims. Original inputs remain unchanged; no primer or reconstructed product determines this target.".into());
    warnings.push("Native reference editing remaps retained annotations, flags intersected features and removes stale translations. Deleted features remain available in original source documents. This target is not physical material or experimental evidence.".into());
    Ok(Report {
        context_revision,
        request: r.clone(),
        target,
        fragments,
        applied: std::mem::take(&mut applied),
        warnings,
    })
}
impl Report {
    /// Replay against exact original inputs, rejecting changed target, coordinates or receipts.
    pub fn validate_sources(&self, inputs: &[PlanInput], expected: &str) -> Result<(), Error> {
        crate::model::bounds::model_preflight(&[&self.target])?;
        if self.context_revision != expected
            || self.applied.len() > MAX_EDITS
            || self.fragments.len() > 32
            || self.warnings.len() > 2048
            || self.warnings.iter().any(|s| s.len() > 16384)
        {
            return Err(bad("stale or oversized edit report"));
        }
        self.validate_binding(expected)?;
        let rebuilt = build(inputs, &self.request)?;
        if rebuilt.context_revision != expected
            || assembly_recipe::hash(self)? != assembly_recipe::hash(&rebuilt)?
        {
            return Err(bad(
                "edit report differs from fresh original-source reconstruction",
            ));
        }
        Ok(())
    }
    /// Bounded self-consistency only; use validate_sources for fresh source replay.
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        shape(&self.request)?;
        crate::model::bounds::model_preflight(&[&self.target])?;
        self.target.validate()?;
        if self.context_revision != expected
            || self.fragments.len() > 32
            || self.applied.len() > MAX_EDITS
            || self.warnings.len() > 2048
            || self.warnings.iter().any(|s| s.len() > 16384)
        {
            return Err(bad("stale or oversized report"));
        }
        let get = |key: &str| {
            self.target
                .metadata
                .get(&format!("{PREFIX}{key}"))
                .filter(|s| s.len() <= 1024 * 1024)
                .ok_or_else(|| bad("missing bounded provenance"))
        };
        let source_pins: Vec<Pin> =
            serde_json::from_str(get("inputs")?).map_err(|_| bad("invalid source pins"))?;
        if source_pins.is_empty()
            || source_pins.len() > 32
            || source_pins.iter().any(|p| {
                p.input_id.trim().is_empty()
                    || p.input_id.len() > 128
                    || p.input_id.chars().any(char::is_control)
                    || p.length == 0
                    || p.length > MAX_BASES
                    || p.revision.len() != 64
                    || !p
                        .revision
                        .bytes()
                        .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
            })
        {
            return Err(bad("invalid original source pins"));
        }
        let base = assembly_recipe::hash(&(
            assembly_substitution_intent::ENGINE_VERSION,
            assembly_junction_edit::ENGINE_VERSION,
            &source_pins,
            &self.request.base,
        ))?;
        let context = assembly_recipe::hash(&(ENGINE_VERSION, base, &self.request))?;
        let sources = source_pins
            .iter()
            .map(|p| (p.input_id.as_str(), p))
            .collect::<BTreeMap<_, _>>();
        if sources.len() != source_pins.len()
            || source_pins
                .iter()
                .try_fold(0usize, |n, p| n.checked_add(p.length))
                .is_none_or(|n| n > MAX_BASES)
        {
            return Err(bad("duplicate or oversized original sources"));
        }
        let lengths = self
            .request
            .base
            .base
            .fragments
            .iter()
            .map(|f| {
                sources
                    .get(f.input_id.as_str())
                    .map(|p| p.length)
                    .ok_or_else(|| bad("unknown original input"))
            })
            .collect::<Result<Vec<_>, _>>()?;
        let adjusted = assembly_junction_edit::intent_adjustment(
            &self.request.base.base,
            &self.request.base.junction_edits,
            &lengths,
        )?;
        let mut expected_fragments = Vec::with_capacity(adjusted.fragments.len());
        let mut offset = 0usize;
        for (i, f) in adjusted.fragments.iter().enumerate() {
            let source = sources[f.input_id.as_str()];
            if source.topology == Topology::Linear && f.end > source.length {
                return Err(bad("linear source exceeds its actual length"));
            }
            let mut length = f.end - f.start;
            if let Some(j) = self
                .request
                .base
                .junction_edits
                .iter()
                .find(|e| e.junction == i)
            {
                length = length
                    .checked_add(
                        j.inserted
                            .bytes()
                            .filter(|b| !b.is_ascii_whitespace())
                            .count(),
                    )
                    .ok_or(Error::TooLarge)?;
            }
            for e in self.request.edits.iter().filter(|e| e.role == f.role) {
                length = length
                    .checked_sub(e.end - e.start)
                    .and_then(|n| n.checked_add(e.replacement.len()))
                    .ok_or(Error::TooLarge)?;
            }
            let end = offset.checked_add(length).ok_or(Error::TooLarge)?;
            expected_fragments.push(PlacedFragment {
                role: f.role.clone(),
                input_id: f.input_id.clone(),
                input_revision: source.revision.clone(),
                start: f.start,
                end: f.end,
                reverse: f.reverse,
                output_start: offset,
                output_end: end,
            });
            offset = end;
        }
        if assembly_recipe::hash(&expected_fragments)? != assembly_recipe::hash(&self.fragments)? {
            return Err(bad("edited contribution geometry changed"));
        }
        let applied = geometry(&source_pins, &self.request)?;
        if context != expected
            || applied != self.applied
            || get("engine")? != ENGINE_VERSION
            || get("context_revision")? != expected
            || get("request")?
                != &serde_json::to_string(&self.request).map_err(|_| Error::TooLarge)?
            || get("applied")?
                != &serde_json::to_string(&self.applied).map_err(|_| Error::TooLarge)?
            || get("contributions")?
                != &serde_json::to_string(&self.fragments).map_err(|_| Error::TooLarge)?
            || get("sequence_blake3")?
                != blake3::hash(self.target.sequence.as_str().as_bytes())
                    .to_hex()
                    .as_str()
            || get("features_revision")? != &assembly_recipe::hash(&self.target.features)?
            || self.fragments.last().map_or(0, |f| f.output_end) != self.target.sequence.len()
            || self.target.sequence.len() > MAX_BASES
            || self.target.sequence.alphabet() != Alphabet::Dna
            || self.target.name != self.request.base.base.name
            || self.target.topology
                != if self.request.base.base.circular {
                    Topology::Circular
                } else {
                    Topology::Linear
                }
            || self.target.metadata.get("duplex_model").map(String::as_str)
                != Some(
                    "Reference-only assembly edit design intent; physical material state unavailable",
                )
            || self
                .target
                .metadata
                .get("assembly:intent:reference_only")
                .map(String::as_str)
                != Some("true")
        {
            return Err(bad("report source/request/target correspondence changed"));
        }
        for e in &self.applied {
            if self
                .target
                .sequence
                .as_str()
                .get(e.output_start..e.output_end)
                != Some(e.replacement.as_str())
            {
                return Err(bad("applied edit differs from target"));
            }
        }
        assembly_recipe::hash(self)?;
        Ok(())
    }
}
