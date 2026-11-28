//! Independently implemented, condition-specific empirical count scoring.
//! No dataset, optimizer code, chemical verification or experimental authority.
//!
//! MODEL pools two directed observation rows for each unique reverse-complement
//! pair. It is our declared model, not a claim of SnapGene/NEB equivalence. For
//! nonpalindromic x with complement r and present physical ends U:
//! p(x)=(C[x,r]+C[r,x])/(sum(y in U) C[x,y]+sum(y in U) C[r,y]).
//! Missing cells (including diagonals) are unknown, never implicit zero. Multiply
//! once per unique pair only when every pair is estimable. Counts are interpreted
//! literally as directed observations; no transpose, pseudocount or abundance
//! weighting is invented. Correct complementarity is not an efficiency estimate.
use crate::cloning::assembly::recipe as assembly_recipe;
use collagene_model::Error;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};
pub const MODEL: &str = "directed-count-pooled-physical-pairs-1";
pub const MAX_CELLS: usize = 65_536;
pub const MAX_ENDS: usize = 256;
pub const MAX_BYTES: usize = 16 << 20;
macro_rules! dto {($($t:tt)*)=>{#[derive(Clone,Debug,PartialEq,Serialize,Deserialize)]#[cfg_attr(feature="json-schema",derive(schemars::JsonSchema))]#[serde(deny_unknown_fields)]$($t)*};}
dto! {pub struct Conditions {
    pub ligase:String,pub buffer:String,pub temperature_protocol:String,
    pub incubation:String,pub restriction_enzymes:Vec<String>,pub substrate_context:String,
}}
dto! {pub struct Provenance {
    pub source_url:String,pub citation:String,pub license:String,pub attribution:String,
    /// Documentary hash of the supplied source, not a claim that it was fetched.
    pub source_blake3:String,
    /// Explain how source counts were converted to directed cells, including
    /// any symmetry expansion already performed by the data provider.
    pub normalization:String,
}}
dto! {pub struct Cell {pub end_5to3:String,pub partner_5to3:String,pub count:u64}}
dto! {pub struct Dataset {
    pub schema_version:u16,pub overhang_length:usize,pub content_revision:String,
    pub conditions:Conditions,pub provenance:Provenance,
    #[serde(deserialize_with="decode_cells")]pub cells:Vec<Cell>,
}}
dto! {pub struct Orientation {
    pub end_5to3:String,pub correct_partner_5to3:String,
    pub correct_count:Option<u64>,
    /// None if any required partner cell is absent; known zero is Some(0).
    #[serde(with="optional_decimal_u128")]
    #[cfg_attr(feature="json-schema",schemars(with="Option<String>"))]
    pub total_count:Option<u128>,pub missing_partners:Vec<String>,
}}
dto! {#[derive(Copy,Eq)] pub enum UnknownReason {MissingCounts,ZeroTotal,PalindromicEnd}}
dto! {pub enum Estimate {
    Known{
        #[serde(with="decimal_u128")]
        #[cfg_attr(feature="json-schema",schemars(with="String"))]
        correct_count:u128,
        #[serde(with="decimal_u128")]
        #[cfg_attr(feature="json-schema",schemars(with="String"))]
        total_count:u128,fidelity:f64
    },
    Unknown{reason:UnknownReason},
}}
dto! {pub struct PairScore {
    pub canonical_end:String,pub reverse_complement:String,
    /// Original request positions, including repeated/complementary declarations.
    pub occurrences:Vec<usize>,pub orientations:Vec<Orientation>,pub estimate:Estimate,
}}
dto! {pub struct Report {
    pub model:String,pub context_revision:String,pub dataset_revision:String,
    pub overhang_length:usize,pub conditions:Conditions,pub provenance:Provenance,
    pub requested_ends:Vec<String>,pub physical_ends:Vec<String>,pub pairs:Vec<PairScore>,
    pub overall_fidelity:Option<f64>,
    /// Natural log of the product when all estimated pair ratios are positive.
    /// Preserves comparison when multiplication underflows; None for unknown or exact zero.
    pub log_overall_fidelity:Option<f64>,pub warnings:Vec<String>,
}}
// JSON numbers cannot portably carry these exact aggregate counts. Human
// formats use bounded canonical decimal strings; binary formats retain rawu128.
pub(crate) mod decimal_u128 {
    use serde::{Deserialize, Serialize};
    pub fn serialize<S: serde::Serializer>(value: &u128, s: S) -> Result<S::Ok, S::Error> {
        if s.is_human_readable() {
            s.serialize_str(&value.to_string())
        } else {
            value.serialize(s)
        }
    }
    pub fn deserialize<'de, D: serde::Deserializer<'de>>(d: D) -> Result<u128, D::Error> {
        if !d.is_human_readable() {
            return u128::deserialize(d);
        }
        struct Decimal;
        impl<'de> serde::de::Visitor<'de> for Decimal {
            type Value = u128;
            fn expecting(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
                f.write_str("a canonical unsigned decimal string of at most39 digits")
            }
            fn visit_str<E: serde::de::Error>(self, v: &str) -> Result<u128, E> {
                if v.is_empty()
                    || v.len() > 39
                    || !v.bytes().all(|b| b.is_ascii_digit())
                    || (v.len() > 1 && v.starts_with('0'))
                {
                    return Err(E::custom("invalid canonical decimal count"));
                }
                v.parse::<u128>()
                    .map_err(|_| E::custom("aggregate count overflow"))
            }
        }
        d.deserialize_str(Decimal)
    }
    pub(super) struct Value(pub u128);
    impl Serialize for Value {
        fn serialize<S: serde::Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
            serialize(&self.0, s)
        }
    }
    impl<'de> Deserialize<'de> for Value {
        fn deserialize<D: serde::Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
            deserialize(d).map(Self)
        }
    }
}
pub(crate) mod optional_decimal_u128 {
    use super::decimal_u128::Value;
    use serde::{Deserialize, Serialize};
    pub fn serialize<S: serde::Serializer>(value: &Option<u128>, s: S) -> Result<S::Ok, S::Error> {
        value.map(Value).serialize(s)
    }
    pub fn deserialize<'de, D: serde::Deserializer<'de>>(d: D) -> Result<Option<u128>, D::Error> {
        Option::<Value>::deserialize(d).map(|v| v.map(|n| n.0))
    }
}
fn bad(s: &str) -> Error {
    Error::Format(format!("ligation fidelity: {s}"))
}
fn text(s: &str, n: usize) -> bool {
    !s.trim().is_empty() && s.len() <= n && !s.chars().any(char::is_control)
}
fn hex(s: &str) -> bool {
    s.len() == 64
        && s.bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}
fn dna(s: &str, k: usize) -> bool {
    s.len() == k && s.bytes().all(|b| matches!(b, b'A' | b'C' | b'G' | b'T'))
}
fn rc(s: &str) -> String {
    s.bytes()
        .rev()
        .map(|b| match b {
            b'A' => 'T',
            b'C' => 'G',
            b'G' => 'C',
            _ => 'A',
        })
        .collect()
}
fn metadata(c: &Conditions, p: &Provenance) -> Result<(), Error> {
    if !text(&c.ligase, 256)
        || !text(&c.buffer, 1024)
        || !text(&c.temperature_protocol, 2048)
        || !text(&c.incubation, 1024)
        || !text(&c.substrate_context, 2048)
        || c.restriction_enzymes.len() > 16
        || c.restriction_enzymes.iter().any(|s| !text(s, 128))
        || !text(&p.source_url, 2048)
        || !text(&p.citation, 2048)
        || !text(&p.license, 512)
        || !text(&p.attribution, 2048)
        || !text(&p.normalization, 2048)
        || !hex(&p.source_blake3)
    {
        return Err(bad(
            "bounded nonempty conditions/provenance and source hash required",
        ));
    }
    Ok(())
}
fn decode_cells<'de, D: serde::Deserializer<'de>>(d: D) -> Result<Vec<Cell>, D::Error> {
    struct V;
    impl<'de> serde::de::Visitor<'de> for V {
        type Value = Vec<Cell>;
        fn expecting(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
            write!(f, "at most {MAX_CELLS} directed count cells")
        }
        fn visit_seq<A: serde::de::SeqAccess<'de>>(
            self,
            mut a: A,
        ) -> Result<Self::Value, A::Error> {
            if a.size_hint().is_some_and(|n| n > MAX_CELLS) {
                return Err(serde::de::Error::custom("matrix cell bound"));
            }
            let mut out = Vec::with_capacity(a.size_hint().unwrap_or(0).min(MAX_CELLS));
            while out.len() < MAX_CELLS {
                match a.next_element()? {
                    Some(v) => out.push(v),
                    None => return Ok(out),
                }
            }
            if a.next_element::<serde::de::IgnoredAny>()?.is_some() {
                return Err(serde::de::Error::custom("matrix cell bound"));
            }
            Ok(out)
        }
    }
    d.deserialize_seq(V)
}
fn cells(data: &Dataset) -> Result<BTreeMap<(&str, &str), u64>, Error> {
    metadata(&data.conditions, &data.provenance)?;
    if data.schema_version != 1
        || !matches!(data.overhang_length, 3 | 4)
        || data.cells.is_empty()
        || data.cells.len() > MAX_CELLS
    {
        return Err(bad(
            "schema1 requires a nonempty bounded3/4-base directed count matrix",
        ));
    }
    let mut map = BTreeMap::new();
    for cell in &data.cells {
        if !dna(&cell.end_5to3, data.overhang_length)
            || !dna(&cell.partner_5to3, data.overhang_length)
            || map
                .insert(
                    (cell.end_5to3.as_str(), cell.partner_5to3.as_str()),
                    cell.count,
                )
                .is_some()
        {
            return Err(bad(
                "invalid concrete uppercase end or duplicate directed cell",
            ));
        }
    }
    Ok(map)
}
/// Canonical data identity ignores supplied cell order and excludes the claimed
/// content_revision itself. Conditions and attribution are part of its binding.
pub fn content_revision(data: &Dataset) -> Result<String, Error> {
    let map = cells(data)?;
    assembly_recipe::hash(&(
        MODEL,
        data.schema_version,
        data.overhang_length,
        &data.conditions,
        &data.provenance,
        map.into_iter().collect::<Vec<_>>(),
    ))
}
impl Dataset {
    pub fn validate(&self) -> Result<(), Error> {
        if !hex(&self.content_revision) || content_revision(self)? != self.content_revision {
            return Err(bad("dataset content revision changed"));
        }
        Ok(())
    }
}
/// Bounded JSON convenience for caller-supplied data. No bundled data or fetch.
pub fn decode(bytes: &[u8]) -> Result<Dataset, Error> {
    if bytes.len() > MAX_BYTES {
        return Err(Error::TooLarge);
    }
    let data: Dataset = serde_json::from_slice(bytes).map_err(|e| bad(&e.to_string()))?;
    data.validate()?;
    Ok(data)
}
fn request(
    ends: &[String],
    k: usize,
) -> Result<(BTreeSet<String>, BTreeMap<String, Vec<usize>>), Error> {
    if ends.is_empty() || ends.len() > MAX_ENDS || ends.iter().any(|s| !dna(s, k)) {
        return Err(bad("choose1..256 concrete uppercase profile-length ends"));
    }
    let mut physical = BTreeSet::new();
    let mut pairs = BTreeMap::<String, Vec<usize>>::new();
    for (i, end) in ends.iter().enumerate() {
        let complement = rc(end);
        physical.insert(end.clone());
        physical.insert(complement.clone());
        pairs
            .entry(end.min(&complement).clone())
            .or_default()
            .push(i);
    }
    Ok((physical, pairs))
}
fn context(
    data_revision: &str,
    k: usize,
    c: &Conditions,
    p: &Provenance,
    ends: &[String],
) -> Result<String, Error> {
    assembly_recipe::hash(&(MODEL, data_revision, k, c, p, ends))
}
pub fn revision(data: &Dataset, ends: &[String]) -> Result<String, Error> {
    data.validate()?;
    request(ends, data.overhang_length)?;
    context(
        &data.content_revision,
        data.overhang_length,
        &data.conditions,
        &data.provenance,
        ends,
    )
}
/// Validated borrowed matrix reusable across a bounded caller-owned search.
/// It cannot be constructed from a historical score report.
pub struct Prepared<'a> {
    data: &'a Dataset,
    map: BTreeMap<(&'a str, &'a str), u64>,
}
impl<'a> Prepared<'a> {
    pub fn new(data: &'a Dataset) -> Result<Self, Error> {
        data.validate()?;
        Ok(Self {
            data,
            map: cells(data)?,
        })
    }
    pub fn score(&self, ends: &[String]) -> Result<Report, Error> {
        let numeric = calculate(self.data.overhang_length, ends, |a, b| {
            self.map.get(&(a, b)).copied()
        })?;
        let scores = numeric.pairs;
        let physical = numeric.physical_ends;
        let overall = numeric.overall_fidelity;
        let log_overall_fidelity = numeric.log_overall_fidelity;
        let report=Report{model:MODEL.into(),context_revision:context(&self.data.content_revision,self.data.overhang_length,&self.data.conditions,&self.data.provenance,ends)?,
            dataset_revision:self.data.content_revision.clone(),overhang_length:self.data.overhang_length,conditions:self.data.conditions.clone(),provenance:self.data.provenance.clone(),
            requested_ends:ends.to_vec(),physical_ends:physical,pairs:scores,overall_fidelity:overall,log_overall_fidelity,
            warnings:vec!["Condition-specific count-based ranking only, not ligation efficiency, reaction yield, chemical verification or experimental success.".into(),
                "Each unique reverse-complement pair contributes once. Both directed rows and every physically present partner, including self-mispairs, enter its conditional denominator. Equal partner abundance and independence of unique pair outcomes are modeling assumptions; duplicate requested ends do not add stoichiometric weighting.".into(),
                "Caller-supplied provenance, source hash and licensing are documentary declarations; this module neither fetches nor authenticates the original experiment. Floating ratios can round; positive products can underflow to0, so log_overall_fidelity retains their logarithm. No uncertainty intervals or pseudocounts are inferred. A numerical1 means no mismatches observed in the supplied competing counts, not guaranteed specificity.".into()]};
        report.validate_binding(&report.context_revision)?;
        Ok(report)
    }
}
pub(crate) struct NumericScore {
    pub physical_ends: Vec<String>,
    pub pairs: Vec<PairScore>,
    pub overall_fidelity: Option<f64>,
    pub log_overall_fidelity: Option<f64>,
}
pub(crate) fn calculate(
    k: usize,
    ends: &[String],
    lookup: impl Fn(&str, &str) -> Option<u64>,
) -> Result<NumericScore, Error> {
    request(ends, k)?; // Preserve the original API's nonempty intended-set guard.
    calculate_with_competitors(k, ends, &[], lookup)
}
pub(crate) fn physical_scope(
    k: usize,
    ends: &[String],
    extras: &[String],
) -> Result<(BTreeSet<String>, BTreeMap<String, Vec<usize>>), Error> {
    if !matches!(k, 3 | 4)
        || ends.len() > MAX_ENDS
        || extras.len() > MAX_ENDS
        || ends.iter().chain(extras).any(|s| !dna(s, k))
    {
        return Err(bad(
            "invalid bounded intended joins or literal physical competitors",
        ));
    }
    let (mut physical, pairs) = if ends.is_empty() {
        (BTreeSet::new(), BTreeMap::new())
    } else {
        request(ends, k)?
    };
    physical.extend(extras.iter().cloned()); // Actual singleton ends: never invent RC.
    Ok((physical, pairs))
}
pub(crate) fn calculate_with_competitors(
    k: usize,
    ends: &[String],
    extras: &[String],
    lookup: impl Fn(&str, &str) -> Option<u64>,
) -> Result<NumericScore, Error> {
    let (physical, pairs) = physical_scope(k, ends, extras)?;
    let scores = pairs
        .into_iter()
        .map(|(end, occurrences)| {
            let complement = rc(&end);
            let oriented = if end == complement {
                vec![end.clone()]
            } else {
                vec![end.clone(), complement.clone()]
            };
            let mut orientations = vec![];
            for source in oriented {
                let partner = rc(&source);
                let mut total = 0u128;
                let mut missing = vec![];
                for destination in &physical {
                    match lookup(source.as_str(), destination.as_str()) {
                        Some(count) => {
                            total = total
                                .checked_add(u128::from(count))
                                .ok_or(Error::TooLarge)?
                        }
                        None => missing.push(destination.clone()),
                    }
                }
                orientations.push(Orientation {
                    correct_count: lookup(source.as_str(), partner.as_str()),
                    end_5to3: source,
                    correct_partner_5to3: partner,
                    total_count: missing.is_empty().then_some(total),
                    missing_partners: missing,
                });
            }
            let estimate = estimate(&orientations, end == complement)?;
            Ok(PairScore {
                canonical_end: end,
                reverse_complement: complement,
                occurrences,
                orientations,
                estimate,
            })
        })
        .collect::<Result<Vec<_>, Error>>()?;
    let overall = if ends.is_empty() {
        None
    } else {
        overall(&scores)
    };
    let log_overall_fidelity = if ends.is_empty() {
        None
    } else {
        log_overall(&scores)
    };
    Ok(NumericScore {
        physical_ends: physical.into_iter().collect(),
        pairs: scores,
        overall_fidelity: overall,
        log_overall_fidelity,
    })
}
fn estimate(rows: &[Orientation], palindrome: bool) -> Result<Estimate, Error> {
    if palindrome {
        return Ok(Estimate::Unknown {
            reason: UnknownReason::PalindromicEnd,
        });
    }
    if rows.iter().any(|r| {
        !r.missing_partners.is_empty() || r.correct_count.is_none() || r.total_count.is_none()
    }) {
        return Ok(Estimate::Unknown {
            reason: UnknownReason::MissingCounts,
        });
    }
    let mut correct = 0u128;
    let mut total = 0u128;
    for row in rows {
        correct = correct
            .checked_add(u128::from(row.correct_count.ok_or(Error::TooLarge)?))
            .ok_or(Error::TooLarge)?;
        total = total
            .checked_add(row.total_count.ok_or(Error::TooLarge)?)
            .ok_or(Error::TooLarge)?;
    }
    Ok(if total == 0 {
        Estimate::Unknown {
            reason: UnknownReason::ZeroTotal,
        }
    } else {
        Estimate::Known {
            correct_count: correct,
            total_count: total,
            fidelity: correct as f64 / total as f64,
        }
    })
}
fn overall(pairs: &[PairScore]) -> Option<f64> {
    pairs.iter().try_fold(1.0, |product, p| match p.estimate {
        Estimate::Known { fidelity, .. } => Some(product * fidelity),
        Estimate::Unknown { .. } => None,
    })
}
fn log_overall(pairs: &[PairScore]) -> Option<f64> {
    pairs.iter().try_fold(0.0, |sum, p| match &p.estimate {
        Estimate::Known { fidelity, .. } if *fidelity > 0.0 => Some(sum + fidelity.ln()),
        _ => None,
    })
}
pub fn score(data: &Dataset, ends: &[String]) -> Result<Report, Error> {
    Prepared::new(data)?.score(ends)
}
impl Report {
    /// Read-only consistency/binding check, never a scientific recomputation or
    /// authentication of imported observation counts.
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        metadata(&self.conditions, &self.provenance)?;
        if !matches!(self.overhang_length, 3 | 4) {
            return Err(bad("unsupported end length"));
        }
        request(&self.requested_ends, self.overhang_length)?;
        if self.model != MODEL
            || !hex(&self.dataset_revision)
            || !matches!(self.overhang_length, 3 | 4)
            || self.context_revision != expected
            || self.context_revision
                != context(
                    &self.dataset_revision,
                    self.overhang_length,
                    &self.conditions,
                    &self.provenance,
                    &self.requested_ends,
                )?
            || self.warnings.len() > 16
            || self.warnings.iter().any(|w| w.len() > 4096)
        {
            return Err(bad("score binding or metadata changed"));
        }
        validate_numeric(
            self.overhang_length,
            &self.requested_ends,
            &self.physical_ends,
            &self.pairs,
            self.overall_fidelity,
            self.log_overall_fidelity,
        )?;
        Ok(())
    }
}

pub(crate) fn validate_numeric(
    k: usize,
    requested: &[String],
    physical_ends: &[String],
    scores: &[PairScore],
    overall_fidelity: Option<f64>,
    log_overall_fidelity: Option<f64>,
) -> Result<(), Error> {
    request(requested, k)?; // Original DTO keeps its original nonempty contract.
    validate_numeric_with_competitors(
        k,
        requested,
        &[],
        physical_ends,
        scores,
        overall_fidelity,
        log_overall_fidelity,
    )
}
pub(crate) fn validate_numeric_with_competitors(
    k: usize,
    requested: &[String],
    extras: &[String],
    physical_ends: &[String],
    scores: &[PairScore],
    overall_fidelity: Option<f64>,
    log_overall_fidelity: Option<f64>,
) -> Result<(), Error> {
    let (physical, pairs) = physical_scope(k, requested, extras)?;
    if physical_ends != physical.iter().cloned().collect::<Vec<_>>().as_slice()
        || scores.len() != pairs.len()
    {
        return Err(bad("physical end scope changed"));
    }
    for (pair, (canonical, occurrences)) in scores.iter().zip(pairs) {
        let complement = rc(&canonical);
        let sources = if canonical == complement {
            vec![canonical.clone()]
        } else {
            vec![canonical.clone(), complement.clone()]
        };
        if pair.canonical_end != canonical
            || pair.reverse_complement != complement
            || pair.occurrences != occurrences
            || pair.orientations.len() != sources.len()
        {
            return Err(bad("pair identity changed"));
        }
        for (row, source) in pair.orientations.iter().zip(sources) {
            if row.end_5to3 != source
                || row.correct_partner_5to3 != rc(&source)
                || row.missing_partners.len() > physical.len()
                || row.missing_partners.iter().any(|v| !physical.contains(v))
                || row.missing_partners.windows(2).any(|w| w[0] >= w[1])
                || row.total_count.is_some() != row.missing_partners.is_empty()
                || row.correct_count.is_none()
                    != row.missing_partners.contains(&row.correct_partner_5to3)
                || row.total_count.is_some_and(|n| {
                    n > u128::from(u64::MAX) * (physical.len() as u128)
                        || row.correct_count.is_some_and(|v| u128::from(v) > n)
                })
            {
                return Err(bad("orientation count shape changed"));
            }
        }
        if pair.estimate != estimate(&pair.orientations, canonical == complement)? {
            return Err(bad("pair conditional estimate changed"));
        }
    }
    let expected_overall = if requested.is_empty() {
        None
    } else {
        overall(scores)
    };
    let expected_log = if requested.is_empty() {
        None
    } else {
        log_overall(scores)
    };
    if overall_fidelity != expected_overall || log_overall_fidelity != expected_log {
        return Err(bad("overall estimate changed"));
    }
    Ok(())
}

pub const COMPETITION_MODEL: &str = "directed-count-literal-physical-competition-1";
// Separate versioned DTO: original Report/postcard layout and meaning unchanged.
dto! {pub struct CompetitionReport{
    pub model:String,pub context_revision:String,pub dataset_revision:String,
    pub overhang_length:usize,pub conditions:Conditions,pub provenance:Provenance,
    pub requested_ends:Vec<String>,pub physical_competitors:Vec<String>,pub physical_ends:Vec<String>,pub pairs:Vec<PairScore>,
    pub overall_fidelity:Option<f64>,pub log_overall_fidelity:Option<f64>,pub warnings:Vec<String>,
}}
fn competition_context(
    revision: &str,
    k: usize,
    c: &Conditions,
    p: &Provenance,
    joins: &[String],
    extras: &[String],
) -> Result<String, Error> {
    assembly_recipe::hash(&(COMPETITION_MODEL, revision, k, c, p, joins, extras))
}
pub fn competition_revision(
    data: &Dataset,
    joins: &[String],
    extras: &[String],
) -> Result<String, Error> {
    data.validate()?;
    physical_scope(data.overhang_length, joins, extras)?;
    competition_context(
        &data.content_revision,
        data.overhang_length,
        &data.conditions,
        &data.provenance,
        joins,
        extras,
    )
}
impl Prepared<'_> {
    pub fn score_competition(
        &self,
        joins: &[String],
        extras: &[String],
    ) -> Result<CompetitionReport, Error> {
        let score =
            calculate_with_competitors(self.data.overhang_length, joins, extras, |a, b| {
                self.map.get(&(a, b)).copied()
            })?;
        Ok(CompetitionReport{model:COMPETITION_MODEL.into(),context_revision:competition_context(&self.data.content_revision,self.data.overhang_length,&self.data.conditions,&self.data.provenance,joins,extras)?,dataset_revision:self.data.content_revision.clone(),overhang_length:self.data.overhang_length,conditions:self.data.conditions.clone(),provenance:self.data.provenance.clone(),requested_ends:joins.to_vec(),physical_competitors:extras.to_vec(),physical_ends:score.physical_ends,pairs:score.pairs,overall_fidelity:score.overall_fidelity,log_overall_fidelity:score.log_overall_fidelity,
            warnings:vec!["Condition-specific directed-observation estimate only; not chemical verification, reaction yield or experimental success.".into(),"Intended junctions contribute unique reverse-complement pairs. Extra physical5′→3′ ends enter denominators literally as singletons; their complements are not invented and they are not multiplied as intended successful joins. Empty intended scope has no overall estimate.".into(),"Equal physical-partner weighting and independent intended-pair outcomes are declared assumptions; duplicate declarations do not supply stoichiometry. Missing counts remain unknown and no pseudocount is invented. Documentary provenance is not authenticated experimental evidence.".into()]})
    }
}
pub fn score_competition(
    data: &Dataset,
    joins: &[String],
    extras: &[String],
) -> Result<CompetitionReport, Error> {
    Prepared::new(data)?.score_competition(joins, extras)
}
impl CompetitionReport {
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        metadata(&self.conditions, &self.provenance)?;
        physical_scope(
            self.overhang_length,
            &self.requested_ends,
            &self.physical_competitors,
        )?;
        if self.model != COMPETITION_MODEL
            || !hex(&self.dataset_revision)
            || self.context_revision != expected
            || expected
                != competition_context(
                    &self.dataset_revision,
                    self.overhang_length,
                    &self.conditions,
                    &self.provenance,
                    &self.requested_ends,
                    &self.physical_competitors,
                )?
            || self.warnings.len() > 16
            || self.warnings.iter().any(|s| s.len() > 4096)
        {
            return Err(bad("competition scope or metadata binding changed"));
        }
        validate_numeric_with_competitors(
            self.overhang_length,
            &self.requested_ends,
            &self.physical_competitors,
            &self.physical_ends,
            &self.pairs,
            self.overall_fidelity,
            self.log_overall_fidelity,
        )
    }
}
