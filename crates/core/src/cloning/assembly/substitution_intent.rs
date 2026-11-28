//! Literal substitutions in original oriented contributions, never edited-target readback.
use crate::{
    cloning::assembly::{
        intent::{self as assembly_intent, ContributionPosition, PlacedFragment},
        junction_edit as assembly_junction_edit,
    },
    documents::workspace::Workspace,
    planning::plan::PlanInput,
};
use collagene_model::{Alphabet, Error, Feature, Location, Molecule, Sequence, Topology};
use serde::{Deserialize, Serialize};
use std::collections::BTreeSet;
pub const ENGINE_VERSION: &str = "assembly-substitution-intent-2";
pub const MAX_SUBSTITUTIONS: usize = 256;
const PREFIX: &str = "assembly_substitution_intent:";
macro_rules! dto{($($t:tt)*)=>{#[derive(Clone,Debug,Serialize,Deserialize)]#[cfg_attr(feature="json-schema",derive(schemars::JsonSchema))]#[serde(deny_unknown_fields)]$($t)*};}
dto! {#[derive(PartialEq,Eq)]pub struct Substitution{pub role:String,pub offset:usize,pub expected:String,pub replacement:String}}
dto! {pub struct Request{pub base:assembly_intent::Request,pub junction_edits:Vec<assembly_junction_edit::Edit>,pub substitutions:Vec<Substitution>}}
dto! {#[derive(PartialEq,Eq)]pub struct AppliedSubstitution{pub role:String,pub input_id:String,pub input_revision:String,pub offset:usize,pub source_position:ContributionPosition,pub output_position:usize,pub expected:String,pub replacement:String}}
dto! {pub struct Report{pub context_revision:String,pub source_context_revision:String,pub request:Request,pub target:Molecule,pub fragments:Vec<PlacedFragment>,pub applied:Vec<AppliedSubstitution>,pub warnings:Vec<String>}}
fn bad(s: &str) -> Error {
    Error::Unsupported(format!("assembly substitution intent: {s}"))
}
pub(crate) fn shape(r: &Request) -> Result<(), Error> {
    assembly_intent::shape(&r.base)?;
    assembly_junction_edit::intent_inserts(&r.base, &r.junction_edits)?;
    if r.substitutions.len() > MAX_SUBSTITUTIONS {
        return Err(Error::TooLarge);
    }
    let mut positions = BTreeSet::new();
    for s in &r.substitutions {
        let f = r
            .base
            .fragments
            .iter()
            .find(|f| f.role == s.role)
            .ok_or_else(|| bad("unknown substitution role"))?;
        if s.offset >= f.selector.length()?
            || s.expected.len() != 1
            || s.replacement.len() != 1
            || !s.expected.bytes().all(|b| b"ACGT".contains(&b))
            || !s.replacement.bytes().all(|b| b"ACGT".contains(&b))
            || s.expected == s.replacement
            || !positions.insert((&s.role, s.offset))
        {
            return Err(bad("invalid/duplicate concrete literal substitution"));
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
    assembly_junction_edit::intent_preflight(inputs, &r.base, &r.junction_edits, w)?;
    let layout = assembly_junction_edit::retained_layout(&r.base, &r.junction_edits)?;
    for s in &r.substitutions {
        let i = r
            .base
            .fragments
            .iter()
            .position(|f| f.role == s.role)
            .ok_or(Error::TooLarge)?;
        if s.offset < layout[i].0 || s.offset >= layout[i].1 {
            return Err(bad(
                "substitution was removed by declared junction trimming",
            ));
        }
    }
    for edit in &r.substitutions {
        let f = r
            .base
            .fragments
            .iter()
            .find(|f| f.role == edit.role)
            .ok_or(Error::TooLarge)?;
        if matches!(f.selector, assembly_intent::Selector::ReferenceRange { .. }) {
            let m = inputs
                .iter()
                .find(|(id, _)| *id == f.input_id)
                .ok_or(Error::TooLarge)?
                .1;
            let original = assembly_intent::reference_sequence(m, &f.selector)?;
            if original.as_str().get(edit.offset..edit.offset + 1) != Some(edit.expected.as_str()) {
                return Err(bad(
                    "expected source base differs from original oriented reference",
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
    context(&assembly_intent::intent_revision(inputs, &r.base, w)?, r)
}
pub(crate) fn applied(
    rows: &[PlacedFragment],
    r: &Request,
) -> Result<Vec<AppliedSubstitution>, Error> {
    let layout = assembly_junction_edit::retained_layout(&r.base, &r.junction_edits)?;
    let inserts = assembly_junction_edit::intent_inserts(&r.base, &r.junction_edits)?;
    let mut starts = vec![];
    let mut cursor = 0usize;
    for ((a, b), insert) in layout.iter().zip(&inserts) {
        starts.push(cursor);
        cursor = cursor
            .checked_add(b - a)
            .and_then(|n| n.checked_add(insert.len()))
            .ok_or(Error::TooLarge)?;
    }
    let mut result = vec![];
    for edit in &r.substitutions {
        let i = rows
            .iter()
            .position(|f| f.role == edit.role)
            .ok_or_else(|| bad("missing original contribution"))?;
        let row = &rows[i];
        if row
            .original_sequence
            .as_str()
            .get(edit.offset..edit.offset + 1)
            != Some(edit.expected.as_str())
        {
            return Err(bad(
                "expected substitution base differs from actual original oriented source strand",
            ));
        }
        let (a, b) = layout[i];
        if edit.offset < a || edit.offset >= b {
            return Err(bad("substitution outside retained original strand"));
        }
        result.push(AppliedSubstitution {
            role: row.role.clone(),
            input_id: row.input_id.clone(),
            input_revision: row.input_revision.clone(),
            offset: edit.offset,
            source_position: row.source_position(edit.offset)?,
            output_position: starts[i] + edit.offset - a,
            expected: edit.expected.clone(),
            replacement: edit.replacement.clone(),
        });
    }
    Ok(result)
}
pub(crate) fn request_json(r: &Request) -> Result<String, Error> {
    #[derive(Serialize)]
    struct D<'a> {
        base: serde_json::Value,
        junction_edits: &'a [assembly_junction_edit::Edit],
        substitutions: &'a [Substitution],
    }
    serde_json::to_string(&D {
        base: serde_json::from_str(&assembly_intent::request_json(&r.base)?)
            .map_err(|_| Error::TooLarge)?,
        junction_edits: &r.junction_edits,
        substitutions: &r.substitutions,
    })
    .map_err(|_| Error::TooLarge)
}
fn affected(l: &Location, positions: &BTreeSet<usize>) -> Result<bool, Error> {
    assembly_junction_edit::location_guard(l)?;
    Ok(l.spans()
        .iter()
        .any(|(a, b, _)| positions.range(*a..*b).next().is_some()))
}
pub(crate) fn build_resolved(
    resolved: &assembly_intent::Resolved<'_>,
    r: &Request,
) -> Result<assembly_intent::Built, Error> {
    shape(r)?;
    let context_revision = context(&resolved.context, r)?;
    let applied = applied(&resolved.rows, r)?;
    let initial = assembly_junction_edit::build_resolved(resolved, &r.base, &r.junction_edits)?;
    let mut target = initial.target;
    let mut warnings = initial.warnings;
    let positions = applied
        .iter()
        .map(|a| a.output_position)
        .collect::<BTreeSet<_>>();
    for feature in &mut target.features {
        if affected(&feature.location, &positions)? {
            if feature.qualifiers.iter().any(|(k, _)| {
                ["transl_except", "anticodon", "rpt_unit_range", "coded_by"]
                    .iter()
                    .any(|v| k.eq_ignore_ascii_case(v))
            }) {
                return Err(bad(
                    "affected coordinate-bearing annotation qualifier requires explicit review",
                ));
            }
            feature.needs_review = true;
            feature
                .qualifiers
                .retain(|(k, _)| !k.eq_ignore_ascii_case("translation"));
            feature.qualifiers.push(("note".into(),Some("Operator-declared substitution changes this reference feature; original source assertions remain in the unchanged source document.".into())));
        }
    }
    let mut dna = target.sequence.as_str().as_bytes().to_vec();
    for a in &applied {
        if dna[a.output_position] != a.expected.as_bytes()[0] {
            return Err(bad(
                "independent source/junction reference differs from expected literal",
            ));
        }
        dna[a.output_position] = a.replacement.as_bytes()[0];
    }
    target.sequence = Sequence::dna(std::str::from_utf8(&dna).map_err(|_| Error::TooLarge)?)?;
    if target.features.len() + applied.len() > 4096 {
        return Err(Error::TooLarge);
    }
    for a in &applied {
        target.features.push(Feature {
            kind: "misc_feature".into(),
            location: Location::range(a.output_position, a.output_position + 1),
            qualifiers: vec![
                (
                    "label".into(),
                    Some(format!(
                        "Declared {}{}{}",
                        a.expected,
                        a.offset + 1,
                        a.replacement
                    )),
                ),
                (
                    "note".into(),
                    Some(
                        "Reference intention, not inherited source DNA or synthesis evidence."
                            .into(),
                    ),
                ),
            ],
            needs_review: true,
        });
    }
    target.metadata.retain(|k, _| {
        !k.starts_with("assembly_intent:") && !k.starts_with("assembly_junction_edit:")
    });
    let fragments = initial.fragments;
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
    warnings.push("Substitution offsets refer to original resolved oriented contribution calls, before junction trims. Actual strand/state/event provenance remains separate from edited intended DNA.".into());
    Ok(assembly_intent::Built {
        target,
        fragments,
        warnings,
    })
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
    let resolved = assembly_intent::resolve(inputs, &r.base, w)?;
    report_resolved(&resolved, r)
}
pub(crate) fn report_resolved(
    resolved: &assembly_intent::Resolved<'_>,
    r: &Request,
) -> Result<Report, Error> {
    let built = build_resolved(resolved, r)?;
    let report = Report {
        context_revision: context(&resolved.context, r)?,
        source_context_revision: resolved.context.clone(),
        request: r.clone(),
        target: built.target,
        fragments: built.fragments,
        applied: applied(&resolved.rows, r)?,
        warnings: built.warnings,
    };
    report.validate_binding(&report.context_revision)?;
    Ok(report)
}
impl Report {
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        shape(&self.request)?;
        if self.context_revision != expected
            || !assembly_intent::hex(&self.source_context_revision)
            || context(&self.source_context_revision, &self.request)? != expected
        {
            return Err(bad("source/substitution context changed"));
        }
        collagene_model::bounds::model_preflight(&[&self.target])?;
        self.target.validate()?;
        assembly_intent::validate_rows(&self.fragments, &self.request.base, &self.target)?;
        let desired = applied(&self.fragments, &self.request)?;
        if desired != self.applied {
            return Err(bad("source-position/substitution evidence changed"));
        }
        let inserts = assembly_junction_edit::intent_inserts(
            &self.request.base,
            &self.request.junction_edits,
        )?;
        let layout = assembly_junction_edit::retained_layout(
            &self.request.base,
            &self.request.junction_edits,
        )?;
        for (i, row) in self.fragments.iter().enumerate() {
            let (a, b) = layout[i];
            if row.retained_start != a || row.retained_end != b {
                return Err(bad("retained original interval changed"));
            }
            let mut expected =
                format!("{}{}", &row.original_sequence.as_str()[a..b], inserts[i]).into_bytes();
            for s in self.applied.iter().filter(|s| s.role == row.role) {
                expected[s.offset - a] = s.replacement.as_bytes()[0];
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
                    "substituted reference differs from original literal declarations",
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
                return Err(bad("retained substitution evidence changed"));
            }
        }
        if self.target.name != self.request.base.name
            || self.target.topology
                != if self.request.base.circular {
                    Topology::Circular
                } else {
                    Topology::Linear
                }
            || self.target.sequence.alphabet() != Alphabet::Dna
        {
            return Err(bad("intended reference identity changed"));
        }
        assembly_intent::validate_material_descriptors(
            &self.target,
            &self.request.base,
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
