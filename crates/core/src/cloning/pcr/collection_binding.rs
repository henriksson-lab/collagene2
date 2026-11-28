//! Source/event adapter and existing perfect-duplex thermodynamics for collection hits.
use super::primers::{ThermoConditions, TmReport};
use crate::encoding::bounded_json;
use collagene_model::{Alphabet, Error, Molecule};
use collagene_oligo_binding::{Budget, BudgetLimits, Collection, Options, Report as SearchReport};
use serde::{Deserialize, Serialize};

pub const ENGINE_VERSION: &str = "oligo-collection-review-1";
pub const MAX_REPORT_BYTES: usize = 24 * 1024 * 1024;
#[derive(Debug, Clone, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Request {
    pub source: Molecule,
    pub collection: Collection,
    pub options: Options,
    pub budget: BudgetLimits,
    pub conditions: ThermoConditions,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum TmValue {
    Available(TmReport),
    Unavailable(String),
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct HitThermodynamics {
    pub hit_index: usize,
    pub annealing: TmValue,
    pub whole_oligo: TmValue,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Report {
    pub engine_version: String,
    pub context_revision: String,
    pub source_revision: String,
    pub search: SearchReport,
    #[serde(
        deserialize_with = "collagene_model::decode::bounded_vec::<_, HitThermodynamics, {collagene_oligo_binding::MAX_HITS}>"
    )]
    pub thermodynamics: Vec<HitThermodynamics>,
    pub conditions: ThermoConditions,
    pub report_revision: String,
    #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, String, 16>")]
    pub warnings: Vec<String>,
}
fn bad(s: &str) -> Error {
    Error::Unsupported(format!("oligo collection: {s}"))
}
impl Request {
    pub fn validate(&self) -> Result<(), Error> {
        collagene_model::bounds::model_preflight(&[&self.source])?;
        self.source.validate()?;
        if self.source.sequence.alphabet() != Alphabet::Dna
            || self.source.sequence.is_empty()
            || self.source.sequence.len() > 8 * 1024 * 1024
        {
            return Err(bad("bounded DNA source required"));
        }
        self.collection.validate()?;
        self.options.validate()?;
        Budget::new(self.budget)?;
        let admission = Budget::admission_cost(&self.source.sequence, &self.collection)?;
        if self.budget.work < admission.work || self.budget.output_charge < admission.output_charge
        {
            return Err(bad(&format!(
                "budget must allow at least {} work units and {} output bytes to bind these inputs; no endpoints have been searched",
                admission.work, admission.output_charge
            )));
        }
        self.conditions.validate()?;
        bounded_json::preflight(self, 16 * 1024 * 1024)
    }
}
pub fn revision(request: &Request) -> Result<String, Error> {
    request.validate()?;
    bounded_json::hash(&(ENGINE_VERSION, request))
}
fn tm(sequence: &str, conditions: &ThermoConditions) -> TmValue {
    match conditions.perfect_duplex_tm(sequence) {
        Ok(v) => TmValue::Available(v),
        Err(e) => TmValue::Unavailable(e.to_string()),
    }
}
pub fn assess(request: &Request) -> Result<Report, Error> {
    let context_revision = revision(request)?;
    let search = collagene_oligo_binding::search(
        &request.source.sequence,
        request.source.topology,
        &request.collection,
        &request.options,
        &mut Budget::new(request.budget)?,
    )?;
    let mut whole = vec![None::<TmValue>; request.collection.oligos.len()];
    let mut thermodynamics = Vec::with_capacity(search.hits.len());
    for (index, hit) in search.hits.iter().enumerate() {
        let oligo = request
            .collection
            .oligos
            .get(hit.oligo_index)
            .ok_or_else(|| bad("hit outside the selected collection"))?;
        hit.validate_columns(
            &request.source.sequence,
            request.source.topology,
            oligo,
            &request.options,
        )?;
        let annealing = if hit.maximum_mismatches != 0 || !hit.certainly_admissible {
            TmValue::Unavailable("No mismatch- or mixture-duplex thermodynamic model is admitted; a perfect-duplex estimate is not assigned to this binding footprint.".into())
        } else {
            tm(
                &oligo.ordered.as_str()[hit.annealing_start..],
                &request.conditions,
            )
        };
        thermodynamics.push(HitThermodynamics {
            hit_index: index,
            annealing,
            whole_oligo: whole[hit.oligo_index]
                .get_or_insert_with(|| tm(oligo.ordered.as_str(), &request.conditions))
                .clone(),
        });
    }
    let mut report = Report {
        engine_version: ENGINE_VERSION.into(),
        context_revision,
        source_revision: request.source.revision(),
        search,
        thermodynamics,
        conditions: request.conditions.clone(),
        report_revision: String::new(),
        warnings: vec![
            "Binding evidence uses a declared ungapped suffix alignment model, not PCR extension, amplification, experimental annealing or absence beyond the actually completed search.".into(),
            "Annealing Tm applies only to an unambiguous perfectly paired footprint under the declared equal-concentration sodium-only DNA-duplex model. Whole-oligo Tm is a hypothetical perfect complement of the entire ordered oligo, including its tail, not the current template-binding Tm or a recommended PCR temperature. Magnesium, additives, modifications and mismatched duplexes are not modeled.".into(),
        ],
    };
    report.report_revision = report.content_revision()?;
    report.validate_binding(&report.context_revision)?;
    Ok(report)
}
impl Report {
    fn content_revision(&self) -> Result<String, Error> {
        bounded_json::hash(&(
            &self.engine_version,
            &self.context_revision,
            &self.source_revision,
            &self.search,
            &self.thermodynamics,
            &self.conditions,
            &self.warnings,
        ))
    }
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        bounded_json::preflight(self, MAX_REPORT_BYTES)?;
        self.conditions.validate()?;
        if self.engine_version != ENGINE_VERSION
            || self.context_revision != expected
            || expected.len() != 64
            || self.search.engine_version != collagene_oligo_binding::ENGINE_VERSION
            || self.search.model != collagene_oligo_binding::MODEL
            || self.source_revision.len() != 64
            || self.thermodynamics.len() != self.search.hits.len()
            || self.report_revision != self.content_revision()?
            || self.warnings.len() > 16
            || self
                .thermodynamics
                .iter()
                .enumerate()
                .any(|(i, t)| i != t.hit_index)
        {
            return Err(bad("report source/options/thermodynamic identity differs"));
        }
        for (hit, values) in self.search.hits.iter().zip(&self.thermodynamics) {
            if (hit.maximum_mismatches != 0 || !hit.certainly_admissible)
                && matches!(&values.annealing, TmValue::Available(_))
            {
                return Err(bad(
                    "mismatch/uncertain footprint cannot claim perfect-duplex binding Tm",
                ));
            }
            for value in [&values.annealing, &values.whole_oligo] {
                match value {
                    TmValue::Available(value) => value.validate_for_conditions(&self.conditions)?,
                    TmValue::Unavailable(reason) if reason.is_empty() || reason.len() > 2048 => {
                        return Err(bad("invalid Tm unavailability reason"));
                    }
                    _ => (),
                }
            }
        }
        Ok(())
    }
    pub fn validate_for(&self, request: &Request) -> Result<(), Error> {
        self.validate_binding(&revision(request)?)?;
        self.search.validate_binding(
            &request.source.sequence,
            request.source.topology,
            &request.collection,
            &request.options,
        )?;
        if bounded_json::hash(self)? != bounded_json::hash(&assess(request)?)? {
            return Err(bad("fresh assessment differs"));
        }
        Ok(())
    }
}
