//! Sequence alignment facts, not PCR extension, kinetics, thermodynamics or success.
//! Explicit ungapped, 3′-anchored annealing suffixes only. No gap-aware claim.
mod decode;
mod identity;
mod scan;

use collagene_model::{Alphabet, Error, Sequence, Topology};
pub use scan::search;
use serde::{Deserialize, Serialize};

pub const ENGINE_VERSION: &str = "oligo-binding-1";
pub const MODEL: &str = "ungapped-three-prime-suffix-score-envelope-1";
pub const MAX_OLIGOS: usize = 1024;
pub const MAX_OLIGO_BASES: usize = 200;
pub const MAX_COLLECTION_BASES: usize = 200_000;
pub const MAX_WORK: u64 = 200_000_000;
pub const MAX_HITS: usize = 10_000;
pub const MAX_OUTPUT_CHARGE: usize = 16 * 1024 * 1024;
pub const MAX_RETAINED_COLUMNS: usize = MAX_OUTPUT_CHARGE / 256;

#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Oligo {
    #[serde(deserialize_with = "decode::identifier")]
    pub id: String,
    #[serde(deserialize_with = "decode::name")]
    pub name: String,
    pub ordered: Sequence,
    /// Declared nonannealing 5′ tail, never silently aligned as template.
    pub declared_tail_length: usize,
}
impl Oligo {
    pub fn validate(&self) -> Result<(), Error> {
        if self.id.len() > 128 || self.name.len() > 256 || self.ordered.len() > MAX_OLIGO_BASES {
            return Err(Error::TooLarge);
        }
        if self.id.is_empty()
            || !self
                .id
                .bytes()
                .all(|b| b.is_ascii_alphanumeric() || b"_.:-".contains(&b))
            || self.name.trim().is_empty()
            || self.name.chars().any(char::is_control)
            || self.ordered.alphabet() != Alphabet::Dna
            || self.declared_tail_length >= self.ordered.len()
        {
            return Err(invalid("invalid bounded DNA oligo identity/sequence/tail"));
        }
        Ok(())
    }
}
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Collection {
    #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, Oligo, MAX_OLIGOS>")]
    pub oligos: Vec<Oligo>,
}
impl Collection {
    pub fn validate(&self) -> Result<(), Error> {
        if self.oligos.is_empty() || self.oligos.len() > MAX_OLIGOS {
            return Err(Error::TooLarge);
        }
        let mut bases = 0usize;
        let mut ids = std::collections::BTreeSet::new();
        for oligo in &self.oligos {
            oligo.validate()?;
            bases += oligo.ordered.len();
            if bases > MAX_COLLECTION_BASES || !ids.insert(&oligo.id) {
                return Err(invalid("collection size/duplicate oligo id"));
            }
        }
        Ok(())
    }
    pub fn content_hash(&self) -> Result<String, Error> {
        self.validate()?;
        Ok(identity::collection(self))
    }
}
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Options {
    pub minimum_annealing_length: usize,
    /// Additional unpaired 5′ bases after the oligo's separately declared tail.
    pub maximum_extra_unpaired: usize,
    pub maximum_mismatches: usize,
    /// Required possible contiguous pairing at the physical 3′ terminus.
    /// Zero explicitly permits terminal mismatches. Certainty is reported separately.
    pub minimum_three_prime_pairing: usize,
    pub match_score: i32,
    pub mismatch_score: i32,
    pub minimum_score: i32,
}
impl Default for Options {
    fn default() -> Self {
        Self {
            minimum_annealing_length: 12,
            maximum_extra_unpaired: 200,
            maximum_mismatches: 2,
            minimum_three_prime_pairing: 0,
            match_score: 2,
            mismatch_score: -3,
            minimum_score: 12,
        }
    }
}
impl Options {
    pub fn validate(&self) -> Result<(), Error> {
        if !(1..=MAX_OLIGO_BASES).contains(&self.minimum_annealing_length)
            || self.maximum_extra_unpaired > MAX_OLIGO_BASES
            || self.maximum_mismatches > MAX_OLIGO_BASES
            || self.minimum_three_prime_pairing > MAX_OLIGO_BASES
            || !(1..=1000).contains(&self.match_score)
            || !(-1000..=0).contains(&self.mismatch_score)
            || !(-200_000..=200_000).contains(&self.minimum_score)
        {
            return Err(invalid("invalid ungapped alignment options"));
        }
        Ok(())
    }
    pub fn content_hash(&self) -> Result<String, Error> {
        self.validate()?;
        Ok(identity::options(self))
    }
}
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum Orientation {
    /// Oligo anneals to the reference complement; extension direction is forward.
    Forward,
    /// Oligo anneals to the reference strand; extension direction is reverse.
    Reverse,
}
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum Pairing {
    Certain,
    Possible,
    Mismatch,
}
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Column {
    pub oligo_index: usize,
    /// Original plus-reference coordinate, including descending reverse footprints.
    pub source_position: usize,
    pub oligo_base: char,
    pub source_base: char,
    /// Actual antiparallel partner: complement of source for Forward, source for Reverse.
    pub template_base: char,
    pub pairing: Pairing,
}
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum ScoreScope {
    /// All admissible best-score suffix ties at this endpoint/orientation.
    ExactConcreteOptimum,
    /// All possibly admissible suffixes not ruled out by a guaranteed competitor.
    /// A conservative score-bound superset, not exact ambiguity resolution.
    ConservativeAmbiguousEnvelope,
}
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Hit {
    pub oligo_index: usize,
    #[serde(deserialize_with = "decode::identifier")]
    pub oligo_id: String,
    pub orientation: Orientation,
    pub three_prime_source_position: usize,
    /// Forward-coordinate footprint start, traversed once at most; reverse columns descend.
    pub source_start: usize,
    pub length: usize,
    pub crosses_origin: bool,
    pub declared_tail_length: usize,
    pub extra_unpaired_length: usize,
    pub annealing_start: usize,
    pub minimum_mismatches: usize,
    pub maximum_mismatches: usize,
    pub minimum_score: i32,
    pub maximum_score: i32,
    pub three_prime_certain: usize,
    pub three_prime_possible: usize,
    pub certainly_admissible: bool,
    pub score_scope: ScoreScope,
    /// Oligo 5′→3′ order; contains the entire admitted annealing footprint.
    #[serde(deserialize_with = "decode::columns")]
    pub columns: Vec<Column>,
}
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum Limit {
    Work,
    HitRetention,
    OutputCharge,
}
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum Status {
    Unsearched,
    Complete,
    Limited(Limit),
}
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct OligoProgress {
    #[serde(deserialize_with = "decode::identifier")]
    pub oligo_id: String,
    pub status: Status,
    pub endpoints_total: u64,
    pub endpoints_examined: u64,
    /// Includes a fully scored endpoint rejected atomically by the output budget.
    pub alignments_found: u64,
    pub alignments_retained: usize,
}
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct BudgetLimits {
    pub work: u64,
    pub hits: usize,
    pub output_charge: usize,
}
impl Default for BudgetLimits {
    fn default() -> Self {
        Self {
            work: 20_000_000,
            hits: 1000,
            output_charge: 8 * 1024 * 1024,
        }
    }
}
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Usage {
    pub work: u64,
    pub hits: usize,
    pub output_charge: usize,
}
/// Shared across calls/collections. Limits are never reset implicitly.
#[derive(Clone, Debug)]
pub struct Budget {
    limits: BudgetLimits,
    used: Usage,
}
impl Budget {
    /// Charge needed to bind validated inputs and retain progress metadata,
    /// before any endpoint can be examined. Smaller allowances cannot produce a
    /// source-bound report; meeting this charge does not promise a complete scan.
    pub fn admission_cost(source: &Sequence, collection: &Collection) -> Result<Usage, Error> {
        check_source(source)?;
        collection.validate()?;
        Ok(Usage {
            work: source.len() as u64
                + collection
                    .oligos
                    .iter()
                    .map(|o| (o.ordered.len() + o.id.len() + o.name.len()) as u64)
                    .sum::<u64>()
                + 256,
            hits: 0,
            output_charge: 1024
                + collection
                    .oligos
                    .iter()
                    .map(|o| 256 + o.id.len())
                    .sum::<usize>(),
        })
    }
    pub fn new(limits: BudgetLimits) -> Result<Self, Error> {
        if limits.work > MAX_WORK
            || limits.hits > MAX_HITS
            || limits.output_charge > MAX_OUTPUT_CHARGE
        {
            return Err(Error::TooLarge);
        }
        Ok(Self {
            limits,
            used: Usage::default(),
        })
    }
    pub fn limits(&self) -> BudgetLimits {
        self.limits
    }
    pub fn used(&self) -> Usage {
        self.used
    }
    pub(crate) fn work(&mut self, amount: u64) -> bool {
        if amount > self.limits.work - self.used.work {
            false
        } else {
            self.used.work += amount;
            true
        }
    }
    pub(crate) fn output(&mut self, hits: usize, charge: usize) -> Result<(), Limit> {
        if hits > self.limits.hits - self.used.hits {
            return Err(Limit::HitRetention);
        }
        if charge > self.limits.output_charge - self.used.output_charge {
            return Err(Limit::OutputCharge);
        }
        self.used.hits += hits;
        self.used.output_charge += charge;
        Ok(())
    }
}
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Report {
    #[serde(deserialize_with = "decode::identifier")]
    pub engine_version: String,
    #[serde(deserialize_with = "decode::identifier")]
    pub model: String,
    /// Sequence+topology identity only; never a molecule/event/material authority.
    #[serde(deserialize_with = "decode::identifier")]
    pub source_content_hash: String,
    #[serde(deserialize_with = "decode::identifier")]
    pub collection_hash: String,
    #[serde(deserialize_with = "decode::identifier")]
    pub options_hash: String,
    pub options: Options,
    pub budget_limits: BudgetLimits,
    pub budget_before: Usage,
    pub budget_after: Usage,
    pub complete: bool,
    #[serde(
        deserialize_with = "collagene_model::decode::bounded_vec::<_, OligoProgress, MAX_OLIGOS>"
    )]
    pub oligos: Vec<OligoProgress>,
    #[serde(deserialize_with = "decode::hits")]
    pub hits: Vec<Hit>,
}
pub fn source_content_hash(source: &Sequence, topology: Topology) -> Result<String, Error> {
    check_source(source)?;
    Ok(identity::source(source, topology))
}
fn check_source(source: &Sequence) -> Result<(), Error> {
    if source.alphabet() != Alphabet::Dna {
        return Err(invalid("binding source requires DNA"));
    }
    if source.is_empty() || source.len() > collagene_model::MAX_SEQUENCE_LEN {
        return Err(Error::TooLarge);
    }
    Ok(())
}
fn invalid(message: &str) -> Error {
    Error::Unsupported(format!("oligo binding: {message}"))
}
fn mask(base: u8) -> u8 {
    match base {
        b'A' => 1,
        b'C' => 2,
        b'G' => 4,
        b'T' => 8,
        b'R' => 5,
        b'Y' => 10,
        b'S' => 6,
        b'W' => 9,
        b'K' => 12,
        b'M' => 3,
        b'B' => 14,
        b'D' => 13,
        b'H' => 11,
        b'V' => 7,
        b'N' => 15,
        _ => 0,
    }
}
fn complement(base: u8) -> u8 {
    match base {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        b'R' => b'Y',
        b'Y' => b'R',
        b'S' => b'S',
        b'W' => b'W',
        b'K' => b'M',
        b'M' => b'K',
        b'B' => b'V',
        b'V' => b'B',
        b'D' => b'H',
        b'H' => b'D',
        b'N' => b'N',
        _ => 0,
    }
}
fn pairing(oligo: u8, template: u8) -> Pairing {
    let a = mask(oligo);
    let b = mask(complement(template));
    if a & b == 0 {
        Pairing::Mismatch
    } else if a == b && a.count_ones() == 1 {
        Pairing::Certain
    } else {
        Pairing::Possible
    }
}

#[cfg(test)]
mod tests;

impl Report {
    /// Validate current context, quotas, counters and individual aligned facts.
    /// Does not authenticate an untrusted claim of completeness or optimality:
    /// current scientific claims require trusted recomputation with `search`.
    pub fn validate_binding(
        &self,
        source: &Sequence,
        topology: Topology,
        collection: &Collection,
        options: &Options,
    ) -> Result<(), Error> {
        check_source(source)?;
        collection.validate()?;
        options.validate()?;
        if self.engine_version != ENGINE_VERSION || self.model != MODEL {
            return Err(invalid("unsupported report engine/model version"));
        }
        if self.source_content_hash != identity::source(source, topology)
            || self.collection_hash != identity::collection(collection)
            || self.options_hash != identity::options(options)
            || self.options != *options
            || self.oligos.len() != collection.oligos.len()
            || self.hits.len() > MAX_HITS
        {
            return Err(invalid(
                "report does not bind current source/collection/options",
            ));
        }
        Budget::new(self.budget_limits)?;
        let before = self.budget_before;
        let after = self.budget_after;
        if before.work > after.work
            || before.hits > after.hits
            || before.output_charge > after.output_charge
            || after.work > self.budget_limits.work
            || after.hits > self.budget_limits.hits
            || after.output_charge > self.budget_limits.output_charge
            || after.hits - before.hits != self.hits.len()
        {
            return Err(invalid("invalid report budget accounting"));
        }
        let mut retained = vec![0usize; collection.oligos.len()];
        let mut columns = 0usize;
        for hit in &self.hits {
            let oligo = collection
                .oligos
                .get(hit.oligo_index)
                .ok_or_else(|| invalid("unknown hit oligo index"))?;
            columns = columns
                .checked_add(hit.columns.len())
                .ok_or(Error::TooLarge)?;
            if columns > MAX_RETAINED_COLUMNS {
                return Err(Error::TooLarge);
            }
            hit.validate_columns(source, topology, oligo, options)?;
            retained[hit.oligo_index] += 1;
        }
        let mut incomplete_seen = false;
        for (index, (progress, oligo)) in self.oligos.iter().zip(&collection.oligos).enumerate() {
            if progress.oligo_id != oligo.id
                || progress.alignments_retained != retained[index]
                || progress.alignments_found < progress.alignments_retained as u64
                || progress.endpoints_total
                    != scan::total_endpoints(source, topology, oligo, options)
                || progress.endpoints_examined > progress.endpoints_total
                || (incomplete_seen && progress.status != Status::Unsearched)
            {
                return Err(invalid("invalid per-oligo progress"));
            }
            match progress.status {
                Status::Complete
                    if progress.endpoints_examined != progress.endpoints_total
                        || progress.alignments_found != progress.alignments_retained as u64 =>
                {
                    return Err(invalid("invalid completed oligo counters"));
                }
                Status::Unsearched
                    if progress.endpoints_examined != 0 || progress.alignments_found != 0 =>
                {
                    return Err(invalid("unsearched oligo contains search claims"));
                }
                Status::Complete => {}
                _ => incomplete_seen = true,
            }
        }
        if self.complete == incomplete_seen {
            return Err(invalid("invalid collection completion claim"));
        }
        Ok(())
    }
}
