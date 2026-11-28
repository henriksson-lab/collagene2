//! Original oriented-contribution indels, independent of preparation/assembly products.
use crate::{
    cloning::assembly::{
        intent::{self as assembly_intent, PlacedFragment, Selector},
        junction_edit as assembly_junction_edit,
        substitution_intent as assembly_substitution_intent,
    },
    documents::workspace::Workspace,
    planning::plan::PlanInput,
};
use collagene_model::{Alphabet, Error, Molecule, Topology};
use serde::{Deserialize, Serialize};
pub const ENGINE_VERSION: &str = "assembly-edit-intent-2";
pub const MAX_EDITS: usize = 256;
const PREFIX: &str = "assembly_edit_intent:";
macro_rules! dto{($($t:tt)*)=>{#[derive(Clone,Debug,Serialize,Deserialize)]#[cfg_attr(feature="json-schema",derive(schemars::JsonSchema))]#[serde(deny_unknown_fields)]$($t)*};}
dto! {#[derive(PartialEq,Eq)]pub struct Edit{pub role:String,pub start:usize,pub end:usize,pub expected:String,pub replacement:String}}
dto! {pub struct Request{pub base:assembly_substitution_intent::Request,pub edits:Vec<Edit>}}
dto! {#[derive(PartialEq,Eq)]pub struct AppliedEdit{pub role:String,pub input_id:String,pub input_revision:String,pub start:usize,pub end:usize,pub selector:Selector,pub original_output_start:usize,pub original_output_end:usize,pub output_start:usize,pub output_end:usize,pub expected:String,pub replacement:String}}
dto! {pub struct Report{pub context_revision:String,pub source_context_revision:String,pub request:Request,pub target:Molecule,pub fragments:Vec<PlacedFragment>,pub applied:Vec<AppliedEdit>,pub warnings:Vec<String>}}
fn bad(s: &str) -> Error {
    Error::Unsupported(format!("assembly edit intent: {s}"))
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
    assembly_substitution_intent::shape(&r.base)?;
    if r.edits.len() > MAX_EDITS {
        return Err(Error::TooLarge);
    }
    let mut total = 0usize;
    for (i, e) in r.edits.iter().enumerate() {
        let f = r
            .base
            .base
            .fragments
            .iter()
            .find(|f| f.role == e.role)
            .ok_or_else(|| bad("unknown edit role"))?;
        if e.end.checked_sub(e.start) != Some(e.expected.len())
            || e.end > f.selector.length()?
            || !literal(&e.expected)
            || !literal(&e.replacement)
            || e.expected == e.replacement
        {
            return Err(bad(
                "distinct bounded uppercase literals must match original contribution interval",
            ));
        }
        total = total
            .checked_add(e.expected.len())
            .and_then(|n| n.checked_add(e.replacement.len()))
            .filter(|n| *n <= 100_000)
            .ok_or(Error::TooLarge)?;
        if r.edits[..i]
            .iter()
            .any(|p| p.role == e.role && conflicts(e.start, e.end, p.start, p.end))
            || r.base
                .substitutions
                .iter()
                .any(|s| s.role == e.role && conflicts(e.start, e.end, s.offset, s.offset + 1))
        {
            return Err(bad(
                "overlapping edits or ambiguous insertion/substitution boundary",
            ));
        }
    }
    Ok(())
}
pub fn source_preflight(
    inputs: &[(&str, &Molecule)],
    r: &Request,
    w: Option<&Workspace>,
) -> Result<(), Error> {
    shape(r)?;
    assembly_substitution_intent::source_preflight(inputs, &r.base, w)?;
    let layout = assembly_junction_edit::retained_layout(&r.base.base, &r.base.junction_edits)?;
    for e in &r.edits {
        let i = r
            .base
            .base
            .fragments
            .iter()
            .position(|f| f.role == e.role)
            .ok_or(Error::TooLarge)?;
        if e.start < layout[i].0 || e.end > layout[i].1 {
            return Err(bad("edit is outside retained original source contribution"));
        }
        if let Selector::ReferenceRange { .. } = &r.base.base.fragments[i].selector {
            let f = &r.base.base.fragments[i];
            let m = inputs
                .iter()
                .find(|(id, _)| *id == f.input_id)
                .ok_or(Error::TooLarge)?
                .1;
            let original = assembly_intent::reference_sequence(m, &f.selector)?;
            if original.as_str().get(e.start..e.end) != Some(e.expected.as_str()) {
                return Err(bad(
                    "expected edit literal differs from original oriented reference",
                ));
            }
        }
    }
    crate::encoding::bounded_json::hash(&(inputs, r))?;
    Ok(())
}
fn context(source: &str, r: &Request) -> Result<String, Error> {
    crate::encoding::bounded_json::hash(&(ENGINE_VERSION, source, r))
}
pub fn intent_revision(
    inputs: &[PlanInput],
    r: &Request,
    w: Option<&Workspace>,
) -> Result<String, Error> {
    source_preflight(
        &inputs
            .iter()
            .map(|i| (i.id.as_str(), &i.molecule))
            .collect::<Vec<_>>(),
        r,
        w,
    )?;
    context(
        &assembly_intent::intent_revision(inputs, &r.base.base, w)?,
        r,
    )
}
fn applied(rows: &[PlacedFragment], r: &Request) -> Result<Vec<AppliedEdit>, Error> {
    let layout = assembly_junction_edit::retained_layout(&r.base.base, &r.base.junction_edits)?;
    let insert = assembly_junction_edit::intent_inserts(&r.base.base, &r.base.junction_edits)?;
    let mut positions = vec![];
    let mut pos = 0usize;
    for ((a, b), s) in layout.iter().zip(insert) {
        positions.push(pos);
        pos = pos
            .checked_add(b - a)
            .and_then(|n| n.checked_add(s.len()))
            .ok_or(Error::TooLarge)?;
    }
    let mut result = vec![];
    for e in &r.edits {
        let i = rows
            .iter()
            .position(|f| f.role == e.role)
            .ok_or(Error::TooLarge)?;
        let row = &rows[i];
        let (a, b) = layout[i];
        if e.start < a
            || e.end > b
            || row.original_sequence.as_str().get(e.start..e.end) != Some(e.expected.as_str())
        {
            return Err(bad(
                "edit literal differs from actual original oriented strand or was removed by trimming",
            ));
        }
        result.push(AppliedEdit {
            role: row.role.clone(),
            input_id: row.input_id.clone(),
            input_revision: row.input_revision.clone(),
            start: e.start,
            end: e.end,
            selector: row.selector.clone(),
            original_output_start: positions[i] + e.start - a,
            original_output_end: positions[i] + e.end - a,
            output_start: 0,
            output_end: 0,
            expected: e.expected.clone(),
            replacement: e.replacement.clone(),
        });
    }
    for (i, e) in result.iter().enumerate() {
        if result[..i].iter().any(|p| {
            conflicts(
                e.original_output_start,
                e.original_output_end,
                p.original_output_start,
                p.original_output_end,
            )
        }) {
            return Err(bad(
                "edits conflict at a shared construct boundary; declare one explicit junction edit",
            ));
        }
    }
    let mut order = (0..result.len()).collect::<Vec<_>>();
    order.sort_by_key(|i| result[*i].original_output_start);
    let mut shift = 0i64;
    for i in order {
        let e = &mut result[i];
        e.output_start = usize::try_from(
            (e.original_output_start as i64)
                .checked_add(shift)
                .ok_or(Error::TooLarge)?,
        )
        .map_err(|_| Error::TooLarge)?;
        e.output_end = e
            .output_start
            .checked_add(e.replacement.len())
            .ok_or(Error::TooLarge)?;
        shift = shift
            .checked_add(e.replacement.len() as i64 - e.expected.len() as i64)
            .ok_or(Error::TooLarge)?;
    }
    Ok(result)
}
fn request_json(r: &Request) -> Result<String, Error> {
    #[derive(Serialize)]
    struct D<'a> {
        base: serde_json::Value,
        edits: &'a [Edit],
    }
    serde_json::to_string(&D {
        base: serde_json::from_str(&assembly_substitution_intent::request_json(&r.base)?)
            .map_err(|_| Error::TooLarge)?,
        edits: &r.edits,
    })
    .map_err(|_| Error::TooLarge)
}
pub(crate) fn build_resolved(
    resolved: &assembly_intent::Resolved<'_>,
    r: &Request,
) -> Result<Report, Error> {
    shape(r)?;
    let context_revision = context(&resolved.context, r)?;
    let applied = applied(&resolved.rows, r)?;
    let initial = assembly_substitution_intent::build_resolved(resolved, &r.base)?;
    let mut target = initial.target;
    let warnings = initial.warnings;
    if target
        .sequence
        .len()
        .checked_mul(applied.len())
        .is_none_or(|n| n > 128 << 20)
    {
        return Err(bad("bounded intended-reference copy work exceeded"));
    }
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
                "independent source/substitution reference differs from edit literal",
            ));
        }
        target = target.replace(
            e.original_output_start,
            e.original_output_end,
            &e.replacement,
        )?;
    }
    let mut fragments = initial.fragments;
    let mut cursor = 0usize;
    for f in &mut fragments {
        let mut length = f.output_end - f.output_start;
        for e in applied.iter().filter(|e| e.role == f.role) {
            length = length
                .checked_sub(e.expected.len())
                .and_then(|n| n.checked_add(e.replacement.len()))
                .ok_or(Error::TooLarge)?;
        }
        if length == 0 {
            return Err(bad("edited row must retain nonempty intended DNA"));
        }
        f.output_start = cursor;
        cursor = cursor.checked_add(length).ok_or(Error::TooLarge)?;
        f.output_end = cursor;
    }
    if cursor != target.sequence.len() {
        return Err(bad("edited contribution accounting differs from target"));
    }
    target
        .metadata
        .retain(|k, _| !k.starts_with("assembly_substitution_intent:"));
    for (key, value) in [
        ("engine", ENGINE_VERSION.into()),
        ("context_revision", context_revision.clone()),
        ("source_context_revision", resolved.context.clone()),
        ("request", request_json(r)?),
        (
            "applied",
            serde_json::to_string(&applied).map_err(|_| Error::TooLarge)?,
        ),
        (
            "contributions",
            assembly_intent::placements_json(&fragments)?,
        ),
        (
            "sequence_blake3",
            blake3::hash(target.sequence.as_str().as_bytes())
                .to_hex()
                .to_string(),
        ),
        (
            "features_revision",
            crate::encoding::bounded_json::hash(&target.features)?,
        ),
    ] {
        target.metadata.insert(format!("{PREFIX}{key}"), value);
    }
    let report = Report {
        context_revision,
        source_context_revision: resolved.context.clone(),
        request: r.clone(),
        target,
        fragments,
        applied,
        warnings,
    };
    report.validate_binding(&report.context_revision)?;
    Ok(report)
}
pub fn build(inputs: &[PlanInput], r: &Request, w: Option<&Workspace>) -> Result<Report, Error> {
    source_preflight(
        &inputs
            .iter()
            .map(|i| (i.id.as_str(), &i.molecule))
            .collect::<Vec<_>>(),
        r,
        w,
    )?;
    let resolved = assembly_intent::resolve(inputs, &r.base.base, w)?;
    build_resolved(&resolved, r)
}
impl Report {
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        shape(&self.request)?;
        if self.context_revision != expected
            || !assembly_intent::hex(&self.source_context_revision)
            || context(&self.source_context_revision, &self.request)? != expected
        {
            return Err(bad("source/edit context changed"));
        }
        collagene_model::bounds::model_preflight(&[&self.target])?;
        self.target.validate()?;
        assembly_intent::validate_rows(&self.fragments, &self.request.base.base, &self.target)?;
        let desired = applied(&self.fragments, &self.request)?;
        if desired != self.applied {
            return Err(bad("edit original-strand correspondence changed"));
        }
        let layout = assembly_junction_edit::retained_layout(
            &self.request.base.base,
            &self.request.base.junction_edits,
        )?;
        let insert = assembly_junction_edit::intent_inserts(
            &self.request.base.base,
            &self.request.base.junction_edits,
        )?;
        for (i, row) in self.fragments.iter().enumerate() {
            let (a, b) = layout[i];
            if row.retained_start != a || row.retained_end != b {
                return Err(bad("retained original strand window changed"));
            }
            let mut expected =
                format!("{}{}", &row.original_sequence.as_str()[a..b], insert[i]).into_bytes();
            for s in self
                .request
                .base
                .substitutions
                .iter()
                .filter(|s| s.role == row.role)
            {
                if row.original_sequence.as_str().get(s.offset..s.offset + 1)
                    != Some(s.expected.as_str())
                {
                    return Err(bad("expected substitution differs from original strand"));
                }
                expected[s.offset - a] = s.replacement.as_bytes()[0];
            }
            let mut edits = self
                .applied
                .iter()
                .filter(|e| e.role == row.role)
                .collect::<Vec<_>>();
            edits.sort_by_key(|e| std::cmp::Reverse(e.start));
            for e in edits {
                let start = e.start - a;
                let end = e.end - a;
                if expected.get(start..end) != Some(e.expected.as_bytes()) {
                    return Err(bad("conflicting edit/substitution literals"));
                }
                expected.splice(start..end, e.replacement.bytes());
            }
            if self
                .target
                .sequence
                .as_str()
                .as_bytes()
                .get(row.output_start..row.output_end)
                != Some(expected.as_slice())
            {
                return Err(bad(
                    "edited target differs from independently declared original-strand literals",
                ));
            }
        }
        for (key, value) in [
            ("engine", ENGINE_VERSION.into()),
            ("context_revision", self.context_revision.clone()),
            (
                "source_context_revision",
                self.source_context_revision.clone(),
            ),
            ("request", request_json(&self.request)?),
            (
                "applied",
                serde_json::to_string(&self.applied).map_err(|_| Error::TooLarge)?,
            ),
            (
                "contributions",
                assembly_intent::placements_json(&self.fragments)?,
            ),
            (
                "sequence_blake3",
                blake3::hash(self.target.sequence.as_str().as_bytes())
                    .to_hex()
                    .to_string(),
            ),
            (
                "features_revision",
                crate::encoding::bounded_json::hash(&self.target.features)?,
            ),
        ] {
            if self.target.metadata.get(&format!("{PREFIX}{key}")) != Some(&value) {
                return Err(bad("retained edit evidence changed"));
            }
        }
        if self.target.name != self.request.base.base.name
            || self.target.topology
                != if self.request.base.base.circular {
                    Topology::Circular
                } else {
                    Topology::Linear
                }
            || self.target.sequence.alphabet() != Alphabet::Dna
        {
            return Err(bad("edited reference identity changed"));
        }
        assembly_intent::validate_material_descriptors(
            &self.target,
            &self.request.base.base,
            &self.fragments,
        )?;
        assembly_intent::report_bound(self)
    }
    pub fn validate_sources(
        &self,
        inputs: &[PlanInput],
        expected: &str,
        w: Option<&Workspace>,
    ) -> Result<(), Error> {
        if intent_revision(inputs, &self.request, w)? != expected {
            return Err(bad("original source declaration changed"));
        }
        self.validate_binding(expected)?;
        assembly_intent::validate_reference_rows(inputs, &self.fragments)
    }
}
