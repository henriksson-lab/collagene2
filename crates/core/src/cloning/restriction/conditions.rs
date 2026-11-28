//! Exact formulation/unit-assay documentary assessment. This does not validate
//! digestion, predict kinetics/products or assess methylation. `source_revision`
//! binds catalogue source evidence, never a molecule or experimental sample.
use super::catalog::{CatalogSelection, catalog};
use collagene_model::Error;
use serde::{Deserialize, Serialize};
use std::{collections::BTreeSet, io::Write};

pub const ENGINE_VERSION: &str = "restriction-conditions-1";
pub const MAX_OPTIONS_BYTES: usize = 32 * 1024;
pub const MAX_REPORT_BYTES: usize = 512 * 1024;
const MAX_WITNESS_BYTES: usize = 64 * 1024;
const MAX_SUGGESTIONS: usize = 256;
macro_rules! dto { ($($item:tt)*) => {
    #[derive(Debug,Clone,PartialEq,Eq,Serialize,Deserialize)]
    #[cfg_attr(feature="json-schema",derive(schemars::JsonSchema))]
    #[serde(deny_unknown_fields)] $($item)*
}; }
dto! {pub struct ChosenFormulation {pub enzyme_id:String,pub formulation_id:Option<String>}}
dto! {pub struct ReactionConditions {
    pub schema_version:u16,
    #[serde(deserialize_with="collagene_model::decode::bounded_vec::<_, ChosenFormulation, 64>")]
    pub formulations:Vec<ChosenFormulation>,
    pub buffer:String,
    pub temperature_milli_c:i32,
    pub incubation_minutes:u32,
    pub minimum_activity_percent:u8,
    pub activity_context:Option<String>,
    pub temperature_context:Option<String>,
}}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum EvidenceStatus {
    Known,
    BelowThreshold,
    Unknown,
    Conflict,
    ContextRequired,
    TemperatureMismatch,
}
dto! {pub struct EvidenceAssessment {
    pub status:EvidenceStatus,
    /// Exact indices into `FormulationAssessment.formulation.evidence`.
    #[serde(deserialize_with="collagene_model::decode::bounded_vec::<_, usize, 64>")]
    pub evidence_indices:Vec<usize>,
    #[serde(deserialize_with="collagene_model::decode::bounded_vec::<_, String, 64>")]
    pub available_contexts:Vec<String>,
}}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum Unassessed {
    Unknown,
}
dto! {pub struct FormulationAssessment {
    /// Full documentary record preserves qualifiers, citations, scope and
    /// context, including star warnings. Unassessed evidence is not a prediction.
    pub enzyme_id:String,
    pub formulation:Option<catalog::Formulation>,
    pub activity:EvidenceAssessment,
    pub temperature:EvidenceAssessment,
    pub star_risk:Unassessed,
    pub incubation_duration:Unassessed,
}}
dto! {pub struct CommonCondition {
    pub buffer:String,pub temperature_milli_c:i32,pub activity_context:String,
    pub temperature_context:String,
}}
dto! {pub struct ConditionsReport {
    pub engine_version:String,
    pub catalog_edition:String,pub catalog_content_hash:String,
    pub selection:CatalogSelection,pub options:ReactionConditions,
    pub selection_revision:String,pub options_revision:String,
    pub request_revision:String,
    /// Hash of the included catalogue provenance records, not molecular input.
    pub source_revision:String,
    #[serde(deserialize_with="collagene_model::decode::bounded_vec::<_, catalog::Source, 128>")]
    pub sources:Vec<catalog::Source>,
    #[serde(deserialize_with="collagene_model::decode::bounded_vec::<_, FormulationAssessment, 64>")]
    pub formulations:Vec<FormulationAssessment>,
    /// Combined unit-assay/temperature evidence only. Even Known leaves star
    /// risk and duration unknown and does not establish reaction compatibility.
    pub documentary_status:EvidenceStatus,
    #[serde(deserialize_with="collagene_model::decode::bounded_vec::<_, CommonCondition, 256>")]
    pub common_conditions:Vec<CommonCondition>,
    /// Completeness of the recorded candidate enumeration only. An empty list
    /// does not establish that no real reaction conditions could work.
    pub suggestions_complete:bool,
    pub batch_revision:String,
}}
fn bad(s: &str) -> Error {
    Error::Unsupported(format!("restriction conditions: {s}"))
}
fn text(s: &str, n: usize) -> bool {
    !s.trim().is_empty() && s.len() <= n && !s.chars().any(char::is_control)
}
struct Budget(usize);
impl Write for Budget {
    fn write(&mut self, b: &[u8]) -> std::io::Result<usize> {
        self.0 = self
            .0
            .checked_sub(b.len())
            .ok_or_else(|| std::io::Error::other("bounded reaction evidence"))?;
        Ok(b.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}
fn bounded(v: &impl Serialize, max: usize) -> Result<(), Error> {
    serde_json::to_writer(Budget(max), v).map_err(|_| Error::TooLarge)
}
fn hash(v: &impl Serialize) -> Result<String, Error> {
    bounded(v, MAX_REPORT_BYTES)?;
    Ok(
        blake3::hash(&serde_json::to_vec(v).map_err(|_| Error::TooLarge)?)
            .to_hex()
            .to_string(),
    )
}
impl ReactionConditions {
    pub fn validate(&self, selection: &CatalogSelection) -> Result<(), Error> {
        bounded(selection, super::catalog::MAX_SELECTION_BYTES)?;
        selection.validate()?;
        bounded(self, MAX_OPTIONS_BYTES)?;
        if self.schema_version != 1
            || self.formulations.is_empty()
            || self.formulations.len() > 64
            || self.formulations.len() != selection.enzymes.len()
            || !text(&self.buffer, 128)
            || !(0..=100_000).contains(&self.temperature_milli_c)
            || !(1..=100_800).contains(&self.incubation_minutes)
            || !(1..=100).contains(&self.minimum_activity_percent)
            || self
                .activity_context
                .as_ref()
                .is_some_and(|s| !text(s, 256))
            || self
                .temperature_context
                .as_ref()
                .is_some_and(|s| !text(s, 256))
        {
            return Err(bad("invalid bounded options"));
        }
        let mut ids = BTreeSet::new();
        for (f, e) in self.formulations.iter().zip(&selection.enzymes) {
            if f.enzyme_id != e.id
                || !text(&f.enzyme_id, 96)
                || f.formulation_id.as_ref().is_some_and(|id| !text(id, 96))
                || !ids.insert(&f.enzyme_id)
            {
                return Err(bad("invalid/duplicate exact formulation choice"));
            }
        }
        Ok(())
    }
    pub fn validate_for(&self, c: &catalog::Catalog, s: &CatalogSelection) -> Result<(), Error> {
        self.validate(s)?;
        // Measure before selection validation's normalization allocations.
        bounded(s, super::catalog::MAX_SELECTION_BYTES)?;
        s.validate_for(c)?;
        if self.formulations.len() != s.enzymes.len() {
            return Err(bad("one exact formulation per selected enzyme required"));
        }
        for (f, e) in self.formulations.iter().zip(&s.enzymes) {
            if f.enzyme_id != e.id
                || f.formulation_id.as_ref().is_some_and(|id| {
                    !c.edition()
                        .formulations
                        .iter()
                        .any(|x| &x.id == id && x.enzyme == e.id)
                })
            {
                return Err(bad("ordered exact enzyme/formulation membership mismatch"));
            }
        }
        Ok(())
    }
}
fn assessment(
    f: &catalog::Formulation,
    condition: &catalog::Condition,
    context: Option<&str>,
    o: &ReactionConditions,
) -> EvidenceAssessment {
    let indices = f
        .evidence
        .iter()
        .enumerate()
        .filter(|(_, e)| &e.condition == condition)
        .map(|(i, _)| i)
        .collect::<Vec<_>>();
    let contexts = indices
        .iter()
        .map(|&i| f.evidence[i].context.as_str())
        .collect::<BTreeSet<_>>();
    let available_contexts = contexts.iter().map(|s| (*s).to_string()).collect();
    if context.is_none() && contexts.len() > 1 {
        return EvidenceAssessment {
            status: EvidenceStatus::ContextRequired,
            evidence_indices: indices,
            available_contexts,
        };
    }
    let evidence_indices = indices
        .into_iter()
        .filter(|&i| context.is_none_or(|c| f.evidence[i].context == c))
        .collect::<Vec<_>>();
    let known = evidence_indices
        .iter()
        .map(|&i| &f.evidence[i].value)
        .filter(|v| !matches!(v, catalog::FactValue::Unknown(_)))
        .collect::<Vec<_>>();
    let status = if known.is_empty() {
        EvidenceStatus::Unknown
    } else if known.iter().any(|v| *v != known[0]) {
        EvidenceStatus::Conflict
    } else {
        match known[0] {
            catalog::FactValue::Percent(p) if *p >= o.minimum_activity_percent => {
                EvidenceStatus::Known
            }
            catalog::FactValue::Percent(_) => EvidenceStatus::BelowThreshold,
            catalog::FactValue::LessThanPercent(p) if *p <= o.minimum_activity_percent => {
                EvidenceStatus::BelowThreshold
            }
            catalog::FactValue::TemperatureMilliC(t) if *t == o.temperature_milli_c => {
                EvidenceStatus::Known
            }
            catalog::FactValue::TemperatureMilliC(_) => EvidenceStatus::TemperatureMismatch,
            _ => EvidenceStatus::Unknown,
        }
    };
    EvidenceAssessment {
        status,
        evidence_indices,
        available_contexts,
    }
}
fn combined(statuses: impl Iterator<Item = EvidenceStatus>) -> EvidenceStatus {
    let all = statuses.collect::<Vec<_>>();
    // Preserve hard documentary failures even when other facts are missing.
    [
        EvidenceStatus::Conflict,
        EvidenceStatus::BelowThreshold,
        EvidenceStatus::TemperatureMismatch,
        EvidenceStatus::ContextRequired,
        EvidenceStatus::Unknown,
    ]
    .into_iter()
    .find(|x| all.contains(x))
    .unwrap_or(EvidenceStatus::Known)
}
impl ConditionsReport {
    fn revision(&self) -> Result<String, Error> {
        hash(&(
            ENGINE_VERSION,
            &self.catalog_edition,
            &self.catalog_content_hash,
            &self.selection_revision,
            &self.options_revision,
            &self.source_revision,
            &self.formulations,
            &self.documentary_status,
            &self.common_conditions,
            self.suggestions_complete,
        ))
    }
    /// Historical/structural integrity only. Fresh catalogue authority requires
    /// `validate_for`; serializable checksums cannot establish documentary truth.
    pub fn validate(&self) -> Result<(), Error> {
        bounded(self, MAX_REPORT_BYTES)?;
        bounded(&self.selection, super::catalog::MAX_SELECTION_BYTES)?;
        self.selection.validate()?;
        self.options.validate(&self.selection)?;
        if self.engine_version != ENGINE_VERSION
            || self.formulations.len() != self.selection.enzymes.len()
            || self.options.formulations.len() != self.formulations.len()
            || self.catalog_edition != self.selection.catalog_edition
            || self.catalog_content_hash != self.selection.catalog_content_hash
            || self.sources.len() > 128
            || self.common_conditions.len() > MAX_SUGGESTIONS
            || self.request_revision != documentary_revision(&self.selection, &self.options)?
            || self.selection_revision != hash(&self.selection)?
            || self.options_revision != hash(&self.options)?
            || self.source_revision != hash(&self.sources)?
            || self.batch_revision != self.revision()?
        {
            return Err(bad("report shape/version/binding mismatch"));
        }
        for ((a, e), choice) in self
            .formulations
            .iter()
            .zip(&self.selection.enzymes)
            .zip(&self.options.formulations)
        {
            let evidence_len = a.formulation.as_ref().map_or(0, |f| f.evidence.len());
            if a.enzyme_id != e.id
                || choice.enzyme_id != e.id
                || a.formulation.as_ref().map(|f| &f.id) != choice.formulation_id.as_ref()
                || a.formulation.as_ref().is_some_and(|f| f.enzyme != e.id)
                || evidence_len > 64
                || [&a.activity, &a.temperature].iter().any(|v| {
                    v.evidence_indices.len() > 64
                        || v.available_contexts.len() > 64
                        || v.evidence_indices.iter().any(|&i| i >= evidence_len)
                })
                || (a.formulation.is_none()
                    && (a.activity != unknown() || a.temperature != unknown()))
            {
                return Err(bad("invalid formulation/evidence references"));
            }
        }
        Ok(())
    }
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        self.validate()?;
        if self.request_revision != expected {
            return Err(bad("assessment request binding differs"));
        }
        Ok(())
    }
    /// Reconstruct every documentary decision against the currently pinned
    /// catalogue and exact caller options, including qualifiers and suggestions.
    pub fn validate_for(
        &self,
        c: &catalog::Catalog,
        s: &CatalogSelection,
        o: &ReactionConditions,
    ) -> Result<(), Error> {
        self.validate()?;
        if *self != assess(c, s, o)? {
            return Err(bad("report is stale or differs from fresh assessment"));
        }
        Ok(())
    }
}
fn unknown() -> EvidenceAssessment {
    EvidenceAssessment {
        status: EvidenceStatus::Unknown,
        evidence_indices: Vec::new(),
        available_contexts: Vec::new(),
    }
}
/// Deterministic historical request identity. This does not validate membership
/// against an installed catalogue and grants no live assessment authority.
pub fn documentary_revision(s: &CatalogSelection, o: &ReactionConditions) -> Result<String, Error> {
    o.validate(s)?;
    hash(&(
        ENGINE_VERSION,
        &s.catalog_edition,
        &s.catalog_content_hash,
        s,
        o,
    ))
}
pub fn revision(
    c: &catalog::Catalog,
    s: &CatalogSelection,
    o: &ReactionConditions,
) -> Result<String, Error> {
    o.validate_for(c, s)?;
    documentary_revision(s, o)
}
pub fn assess(
    c: &catalog::Catalog,
    s: &CatalogSelection,
    o: &ReactionConditions,
) -> Result<ConditionsReport, Error> {
    o.validate_for(c, s)?;
    let forms = o
        .formulations
        .iter()
        .map(|choice| {
            choice
                .formulation_id
                .as_ref()
                .and_then(|id| c.edition().formulations.iter().find(|f| &f.id == id))
        })
        .collect::<Vec<_>>();
    let mut source_ids = BTreeSet::new();
    for selected in &s.enzymes {
        if let Some(enzyme) = c.edition().enzymes.iter().find(|e| e.id == selected.id) {
            for citation in &enzyme.citations {
                source_ids.insert(citation.source.as_str());
            }
        }
    }
    for f in forms.iter().flatten() {
        for citation in f
            .citations
            .iter()
            .chain(f.aliases.iter().flat_map(|a| &a.citations))
            .chain(f.evidence.iter().flat_map(|e| &e.citations))
        {
            source_ids.insert(citation.source.as_str());
        }
    }
    let sources = c
        .edition()
        .sources
        .iter()
        .filter(|v| source_ids.contains(v.id.as_str()))
        .collect::<Vec<_>>();
    // Conservative shared witness cap BEFORE cloning full records, contexts or
    // qualifiers. Catalogue validation already bounds each record to 64 rows.
    // Subsequent suggestion strings are separately capped by count and length.
    bounded(&(s, o, &forms, &sources), MAX_WITNESS_BYTES)?;
    let condition = catalog::Condition::BufferActivity {
        buffer: o.buffer.clone(),
        temperature_milli_c: o.temperature_milli_c,
    };
    let formulations = forms
        .iter()
        .zip(&s.enzymes)
        .map(|(f, e)| FormulationAssessment {
            enzyme_id: e.id.clone(),
            activity: f.map_or_else(unknown, |f| {
                assessment(f, &condition, o.activity_context.as_deref(), o)
            }),
            temperature: f.map_or_else(unknown, |f| {
                assessment(
                    f,
                    &catalog::Condition::IncubationTemperature,
                    o.temperature_context.as_deref(),
                    o,
                )
            }),
            formulation: f.map(Clone::clone),
            star_risk: Unassessed::Unknown,
            incubation_duration: Unassessed::Unknown,
        })
        .collect::<Vec<_>>();
    let documentary_status = combined(
        formulations
            .iter()
            .flat_map(|f| [f.activity.status, f.temperature.status]),
    );
    // Every common candidate is a literal buffer/temp/context triple from the
    // first formulation and is then independently checked in every formulation.
    let candidates = forms[0]
        .into_iter()
        .flat_map(|f| f.evidence.iter())
        .filter_map(|e| match &e.condition {
            catalog::Condition::BufferActivity {
                buffer,
                temperature_milli_c,
            } if o
                .activity_context
                .as_ref()
                .is_none_or(|ctx| ctx == &e.context) =>
            {
                Some((buffer.as_str(), *temperature_milli_c, e.context.as_str()))
            }
            _ => None,
        })
        .collect::<BTreeSet<_>>();
    let mut common_conditions = Vec::new();
    let mut suggestions_complete = true;
    let temperature_context = o.temperature_context.as_deref().or_else(|| {
        let contexts = &formulations[0].temperature.available_contexts;
        (contexts.len() == 1).then(|| contexts[0].as_str())
    });
    for (buffer, temperature, context) in candidates {
        let Some(temperature_context) = temperature_context else {
            continue;
        };
        let candidate = ReactionConditions {
            buffer: buffer.into(),
            temperature_milli_c: temperature,
            activity_context: Some(context.into()),
            temperature_context: Some(temperature_context.into()),
            ..o.clone()
        };
        let condition = catalog::Condition::BufferActivity {
            buffer: buffer.into(),
            temperature_milli_c: temperature,
        };
        if forms.iter().all(|f| {
            f.is_some_and(|f| {
                assessment(f, &condition, Some(context), &candidate).status == EvidenceStatus::Known
                    && assessment(
                        f,
                        &catalog::Condition::IncubationTemperature,
                        Some(temperature_context),
                        &candidate,
                    )
                    .status
                        == EvidenceStatus::Known
            })
        }) {
            if common_conditions.len() == MAX_SUGGESTIONS {
                suggestions_complete = false;
                break;
            }
            common_conditions.push(CommonCondition {
                buffer: buffer.into(),
                temperature_milli_c: temperature,
                activity_context: context.into(),
                temperature_context: temperature_context.into(),
            });
        }
    }
    let sources = sources.into_iter().cloned().collect::<Vec<_>>();
    let mut report = ConditionsReport {
        engine_version: ENGINE_VERSION.into(),
        catalog_edition: c.edition().edition.clone(),
        catalog_content_hash: c.content_hash().into(),
        selection: s.clone(),
        options: o.clone(),
        request_revision: revision(c, s, o)?,
        selection_revision: hash(s)?,
        options_revision: hash(o)?,
        source_revision: hash(&sources)?,
        sources,
        formulations,
        documentary_status,
        common_conditions,
        suggestions_complete,
        batch_revision: String::new(),
    };
    bounded(&report, MAX_REPORT_BYTES)?;
    report.batch_revision = report.revision()?;
    report.validate()?;
    Ok(report)
}

#[cfg(test)]
mod tests;
