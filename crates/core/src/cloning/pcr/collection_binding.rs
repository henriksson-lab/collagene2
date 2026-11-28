//! Source/event adapter and existing perfect-duplex thermodynamics for collection hits.
use super::primers::{ThermoConditions, TmReport, canonical_dna_stack};
use crate::encoding::bounded_json;
use collagene_model::{Alphabet, Error, Molecule};
use collagene_oligo_binding::{
    Budget, BudgetLimits, Collection, Options, Report as SearchReport, thermodynamics as mismatch,
};
use serde::{Deserialize, Serialize};

pub const ENGINE_VERSION: &str = "oligo-collection-review-2";
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
    MismatchAvailable(MismatchTm),
    Unavailable(String),
}
pub const MISMATCH_TM_MODEL: &str =
    "single internal DNA mismatch; Allawi/SantaLucia/Peyret; pH7; sodium entropy correction";
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct MismatchTm {
    pub parameters: mismatch::Assessment,
    pub estimate: TmReport,
}
impl MismatchTm {
    fn from_parameters(
        parameters: mismatch::Assessment,
        conditions: &ThermoConditions,
    ) -> Result<Self, Error> {
        parameters.validate()?;
        let mismatch::Outcome::Replacements(plan) = &parameters.outcome else {
            return Err(bad("mismatch parameter replacements required"));
        };
        let baseline = conditions.perfect_duplex_tm(&parameters.top_5_to_3)?;
        let mut h = baseline.delta_h_kcal_per_mol;
        let mut s = baseline.delta_s_cal_per_mol_kelvin_at_1m;
        // A palindromic top does not make this admitted, distinct-partner
        // mismatched duplex self-complementary. Undo that baseline term.
        if baseline.self_complementary {
            s += 1.4;
        }
        for replacement in &plan.replacements {
            let pair = replacement.canonical.top.map(|base| base as u8);
            let (canonical_h, canonical_s) = canonical_dna_stack(&pair)?;
            h += f64::from(replacement.delta_h_tenths_kcal) / 10.0 - canonical_h;
            s += f64::from(replacement.delta_s_tenths_cal_per_kelvin) / 10.0 - canonical_s;
        }
        let estimate = TmReport::from_totals(
            MISMATCH_TM_MODEL,
            conditions,
            parameters.top_5_to_3.len(),
            h,
            s,
            false,
        )?;
        Ok(Self {
            parameters,
            estimate,
        })
    }
    fn validate(&self, conditions: &ThermoConditions) -> Result<(), Error> {
        self.parameters.validate()?;
        if self.parameters.engine != mismatch::ENGINE
            || self.parameters.profile_revision != mismatch::profile_revision()
            || self.parameters.assumptions != mismatch_assumptions(conditions)
            || !matches!(self.parameters.outcome, mismatch::Outcome::Replacements(_))
        {
            return Err(bad("mismatch thermodynamics profile/assumptions differ"));
        }
        let expected = Self::from_parameters(self.parameters.clone(), conditions)?;
        if bounded_json::hash(&expected.estimate)? != bounded_json::hash(&self.estimate)? {
            return Err(bad("mismatch thermodynamic totals differ"));
        }
        Ok(())
    }
}
fn mismatch_assumptions(conditions: &ThermoConditions) -> mismatch::Assumptions {
    mismatch::Assumptions {
        sodium_molar: conditions.sodium_molar,
        total_strand_molar: conditions.total_strand_molar,
        concentration: mismatch::ConcentrationConvention::EqualConcentrationTotalStrands,
        salt: mismatch::SaltConvention::SodiumEntropy1998,
        ph: mismatch::PhConvention::Ph7,
    }
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
    let mismatch_context =
        mismatch::DuplexContext::new(&request.source.sequence, request.source.topology)?;
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
        let parameters = mismatch_context.assess(
            oligo,
            &request.options,
            hit,
            mismatch_assumptions(&request.conditions),
        )?;
        let annealing = match parameters.outcome {
            mismatch::Outcome::PerfectDuplex => tm(&parameters.top_5_to_3, &request.conditions),
            mismatch::Outcome::Replacements(_) => {
                match MismatchTm::from_parameters(parameters, &request.conditions) {
                    Ok(value) => TmValue::MismatchAvailable(value),
                    Err(error) => TmValue::Unavailable(error.to_string()),
                }
            }
            mismatch::Outcome::Unsupported(reason) => TmValue::Unavailable(format!(
                "Mismatch thermodynamics unsupported for this footprint: {reason:?}."
            )),
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
            "Annealing Tm uses either the concrete perfect-duplex model or the admitted single internal DNA-mismatch parameter profile (pH7, two distinct equal-concentration strands and sodium-only entropy correction). Unsupported footprints remain unavailable. Whole-oligo Tm is a hypothetical perfect complement of the entire ordered oligo, including its tail, not the current template-binding Tm or a recommended PCR temperature. Magnesium, additives, modifications, multiple or near-end mismatches and competing self-association are not modeled.".into(),
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
            if matches!(values.whole_oligo, TmValue::MismatchAvailable(_)) {
                return Err(bad(
                    "hypothetical whole-oligo Tm must use a perfect complement",
                ));
            }
            if matches!(values.annealing, TmValue::MismatchAvailable(_))
                && (hit.minimum_mismatches != 1
                    || hit.maximum_mismatches != 1
                    || !hit.certainly_admissible)
            {
                return Err(bad(
                    "mismatch thermodynamics requires one concrete admitted mismatch",
                ));
            }
            if let TmValue::MismatchAvailable(value) = &values.annealing {
                let facts = &value.parameters;
                if facts.top_5_to_3.len() != hit.columns.len()
                    || facts.source_positions.len() != hit.columns.len()
                    || facts.oligo_start != hit.annealing_start
                    || !facts
                        .top_5_to_3
                        .chars()
                        .zip(&hit.columns)
                        .all(|(base, c)| base == c.oligo_base)
                    || !facts
                        .bottom_3_to_5
                        .chars()
                        .zip(&hit.columns)
                        .all(|(base, c)| base == c.template_base)
                    || !facts
                        .source_positions
                        .iter()
                        .zip(&hit.columns)
                        .all(|(position, c)| *position == c.source_position)
                {
                    return Err(bad(
                        "mismatch thermodynamic footprint differs from retained hit columns",
                    ));
                }
            }
            for value in [&values.annealing, &values.whole_oligo] {
                match value {
                    TmValue::Available(value) => value.validate_for_conditions(&self.conditions)?,
                    TmValue::MismatchAvailable(value) => value.validate(&self.conditions)?,
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
