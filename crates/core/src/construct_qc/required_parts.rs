//! Required literal-reference occurrences, using the canonical bounded search.
//! Sequence identity/cardinality never establishes an origin or marker's function.
use super::*;
use crate::sequence::search::{self, SearchMode, SearchOptions};

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum RequiredOrientation {
    Either,
    Forward,
    Reverse,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct RequiredPart {
    pub id: String,
    /// User-declared role only; no role inference is made from this text.
    pub role: String,
    pub reference_id: String,
    pub minimum_copies: usize,
    pub maximum_copies: Option<usize>,
    pub orientation: RequiredOrientation,
}

impl RequiredPart {
    pub(super) fn validate(&self) -> Result<(), Error> {
        if !text(&self.id, 128)
            || !text(&self.role, 256)
            || !text(&self.reference_id, 128)
            || self.minimum_copies > 1024
            || self
                .maximum_copies
                .is_some_and(|n| n > 1024 || n < self.minimum_copies)
        {
            return Err(bad("invalid required literal-part declaration"));
        }
        Ok(())
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct RequiredPartLimits {
    pub max_work: usize,
    pub max_matches: usize,
    pub max_matches_per_reference: usize,
}

impl Default for RequiredPartLimits {
    fn default() -> Self {
        Self {
            max_work: 20_000_000,
            max_matches: 512,
            max_matches_per_reference: 128,
        }
    }
}

impl RequiredPartLimits {
    pub(super) fn validate(&self) -> Result<(), Error> {
        if !(1..=50_000_000).contains(&self.max_work)
            || !(1..=1024).contains(&self.max_matches)
            || !(1..=256).contains(&self.max_matches_per_reference)
        {
            return Err(bad("invalid shared required-part search limits"));
        }
        Ok(())
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct PartMatch {
    #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, Location, 2>")]
    pub coordinates: Vec<Location>,
    /// None for a reverse-complement palindrome: sequence alone cannot orient it.
    pub reverse: Option<bool>,
    pub certain: bool,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct ReferenceSearch {
    pub reference_id: String,
    pub reference_revision: Option<String>,
    pub complete: bool,
    pub work: usize,
    /// One shared result per reference, rather than copying a search per role.
    #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, PartMatch, 256>")]
    pub matches: Vec<PartMatch>,
    pub reason: String,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct RequiredPartReport {
    pub id: String,
    pub role: String,
    pub reference_id: String,
    pub search_index: usize,
    pub verdict: Verdict,
    /// Bounds among retained intervals; incomplete searches cannot prove absence.
    pub certain_copies: usize,
    pub possible_copies: usize,
    pub reason: String,
}

fn scan(
    checks: &Checks<'_>,
    id: &str,
    remaining_work: usize,
    remaining_matches: usize,
) -> Result<ReferenceSearch, Error> {
    let mut result = ReferenceSearch {
        reference_id: id.into(),
        reference_revision: None,
        complete: false,
        work: 0,
        matches: vec![],
        reason: String::new(),
    };
    let Some(reference) = checks.o.references.iter().find(|r| r.id == id) else {
        result.reason = "No pinned reference was supplied; missing reference data cannot establish part absence.".into();
        return Ok(result);
    };
    result.reference_revision = Some(reference.sequence_revision.clone());
    if !checks
        .reference_facts
        .get(id)
        .is_some_and(|r| r.pin_matches)
    {
        result.reason = "Captured reference does not match its declared sequence identity.".into();
        return Ok(result);
    }
    if reference.sequence.alphabet() != checks.m.sequence.alphabet()
        || reference.sequence.len() > 4096
        || !checks.reference_facts.get(id).is_some_and(|r| r.concrete)
    {
        result.reason = "Exhaustive cardinality requires a concrete same-alphabet reference of at most 4,096 bases, within the current canonical search profile.".into();
        return Ok(result);
    }
    if remaining_work == 0 || remaining_matches == 0 {
        result.reason = "The shared job work or retained-match budget is exhausted; no absence/cardinality conclusion is made.".into();
        return Ok(result);
    }
    let search_options = SearchOptions {
        query: reference.sequence.as_str().into(),
        mode: SearchMode::Sequence,
        both_strands: true,
        max_work: remaining_work,
        max_hits: remaining_matches.min(checks.o.required_part_limits.max_matches_per_reference),
        ..Default::default()
    };
    let report = search::search_document(checks.m, &search_options)?;
    result.complete = report.complete;
    result.work = report.work;
    let palindrome = reference.sequence == reference.sequence.reverse_complement()?;
    // The current canonical search already avoids a second palindromic scan.
    // Physical-interval keys also prevent orientation aliases counting twice.
    let mut intervals = BTreeMap::<Vec<(usize, usize)>, PartMatch>::new();
    for hit in report.hits {
        let mut key = hit
            .parts
            .iter()
            .map(|p| (p.start, p.end))
            .collect::<Vec<_>>();
        key.sort_unstable();
        if let Some(old) = intervals.get_mut(&key) {
            old.certain |= hit.certain;
            if old.reverse != hit.parts.first().map(|p| p.reverse) {
                old.reverse = None;
            }
        } else {
            intervals.insert(
                key,
                PartMatch {
                    coordinates: hit
                        .parts
                        .iter()
                        .map(|p| Location::range(p.start, p.end))
                        .collect(),
                    reverse: if palindrome {
                        None
                    } else {
                        hit.parts.first().map(|p| p.reverse)
                    },
                    certain: hit.certain,
                },
            );
        }
    }
    result.matches = intervals.into_values().collect();
    result.reason = if result.complete {
        "Every permitted start in the exact current source was examined on both strands, including one-turn circular matches. Counts concern literal reference intervals, not functional copy number."
    } else {
        "Canonical search stopped at a work or match bound. Retained matches are observations only; zero retained matches is not evidence of absence."
    }.into();
    Ok(result)
}

pub(super) struct Results {
    pub reports: Vec<RequiredPartReport>,
    pub searches: Vec<ReferenceSearch>,
    pub work: usize,
}

pub(super) fn inspect(checks: &mut Checks<'_>) -> Result<Results, Error> {
    let options = checks.o;
    let mut result = Results {
        reports: vec![],
        searches: vec![],
        work: 0,
    };
    let mut cache = BTreeMap::<String, usize>::new();
    let mut retained = 0usize;
    if options.required_parts.is_empty() {
        checks.push(Rule::RequiredPart, Verdict::Indeterminate, &[],
            "No explicit role/reference/cardinality requirements were declared; this does not establish that replication origins or selection markers are present or absent.")?;
    }
    for requirement in &options.required_parts {
        let index = if let Some(&i) = cache.get(&requirement.reference_id) {
            i
        } else {
            let search = scan(
                checks,
                &requirement.reference_id,
                options
                    .required_part_limits
                    .max_work
                    .checked_sub(result.work)
                    .ok_or(Error::TooLarge)?,
                options
                    .required_part_limits
                    .max_matches
                    .checked_sub(retained)
                    .ok_or(Error::TooLarge)?,
            )?;
            result.work = result
                .work
                .checked_add(search.work)
                .filter(|n| *n <= options.required_part_limits.max_work)
                .ok_or(Error::TooLarge)?;
            retained = retained
                .checked_add(search.matches.len())
                .filter(|n| *n <= options.required_part_limits.max_matches)
                .ok_or(Error::TooLarge)?;
            let i = result.searches.len();
            cache.insert(requirement.reference_id.clone(), i);
            result.searches.push(search);
            i
        };
        let search = &result.searches[index];
        let mut certain = 0usize;
        let mut possible = 0usize;
        let mut orientation_unresolved = false;
        let mut coordinates = Vec::new();
        for hit in &search.matches {
            let allowed = match (&requirement.orientation, hit.reverse) {
                (RequiredOrientation::Either, _) => true,
                (RequiredOrientation::Forward, Some(false))
                | (RequiredOrientation::Reverse, Some(true)) => true,
                (_, None) => {
                    orientation_unresolved = true;
                    true
                }
                _ => false,
            };
            if allowed {
                possible += 1;
                certain += usize::from(
                    hit.certain
                        && (requirement.orientation == RequiredOrientation::Either
                            || hit.reverse.is_some()),
                );
                // Only a bounded selectable sample is copied to a finding. The
                // complete retained interval list remains in the shared search.
                if coordinates.len() + hit.coordinates.len() <= 16 {
                    coordinates.extend(hit.coordinates.iter().cloned());
                }
            }
        }
        let verdict = if !search.complete || orientation_unresolved || certain != possible {
            Verdict::Indeterminate
        } else if certain < requirement.minimum_copies
            || requirement.maximum_copies.is_some_and(|n| certain > n)
        {
            Verdict::Fail
        } else {
            Verdict::Pass
        };
        let reason = format!(
            "Declared role '{}': {} certain and {} possible retained literal-reference interval(s); required minimum {}, maximum {:?}, orientation {:?}. {}{} No role/function is inferred from the reference or its absence.",
            requirement.role,
            certain,
            possible,
            requirement.minimum_copies,
            requirement.maximum_copies,
            requirement.orientation,
            search.reason,
            if orientation_unresolved {
                " Palindromic/aliased sequence does not establish the declared orientation."
            } else {
                ""
            },
        );
        checks.push_literal_finding(Rule::RequiredPart, verdict.clone(), &coordinates, &reason)?;
        result.reports.push(RequiredPartReport {
            id: requirement.id.clone(),
            role: requirement.role.clone(),
            reference_id: requirement.reference_id.clone(),
            search_index: index,
            verdict,
            certain_copies: certain,
            possible_copies: possible,
            reason,
        });
    }
    Ok(result)
}

pub(super) fn validate_reports(report: &Report) -> Result<(), Error> {
    let mut work = 0usize;
    let mut matches = 0usize;
    let mut references = BTreeSet::new();
    for search in &report.required_part_searches {
        if !text(&search.reference_id, 128)
            || !references.insert(&search.reference_id)
            || search.reference_revision.as_ref().is_some_and(|s| !hex(s))
            || search.matches.len() > 256
            || !text(&search.reason, 4096)
        {
            return Err(bad("invalid retained reference-search shape"));
        }
        match &search.reference_revision {
            Some(revision) => {
                if !report
                    .references
                    .iter()
                    .any(|r| r.id == search.reference_id && &r.sequence_revision == revision)
                {
                    return Err(bad("reference search differs from captured reference pins"));
                }
            }
            None if search.complete || search.work != 0 || !search.matches.is_empty() => {
                return Err(bad("missing reference cannot yield a completed search"));
            }
            None => (),
        }
        work = work
            .checked_add(search.work)
            .filter(|n| *n <= 50_000_000)
            .ok_or(Error::TooLarge)?;
        matches = matches
            .checked_add(search.matches.len())
            .filter(|n| *n <= 1024)
            .ok_or(Error::TooLarge)?;
        for occurrence in &search.matches {
            if !(1..=2).contains(&occurrence.coordinates.len()) {
                return Err(bad("invalid retained literal-match geometry"));
            }
            for l in &occurrence.coordinates {
                l.validate(report.source_length)?;
                if !matches!(
                    l,
                    Location::Range {
                        before: false,
                        after: false,
                        ..
                    }
                ) {
                    return Err(bad("literal matches require exact ascending source spans"));
                }
            }
        }
    }
    if work != report.required_part_search_work {
        return Err(bad(
            "shared required-part work differs from its reference searches",
        ));
    }
    let mut ids = BTreeSet::new();
    let mut used_searches = BTreeSet::new();
    for required in &report.required_parts {
        let search = report
            .required_part_searches
            .get(required.search_index)
            .ok_or_else(|| bad("required-part search index is outside retained searches"))?;
        used_searches.insert(required.search_index);
        if !text(&required.id, 128)
            || !ids.insert(&required.id)
            || !text(&required.role, 256)
            || required.reference_id != search.reference_id
            || required.certain_copies > required.possible_copies
            || required.possible_copies > search.matches.len()
            || !text(&required.reason, 4096)
            || (!search.complete && required.verdict != Verdict::Indeterminate)
        {
            return Err(bad("invalid retained required-part bounds/binding"));
        }
    }
    if used_searches.len() != report.required_part_searches.len() {
        return Err(bad(
            "retained reference search has no corresponding declaration",
        ));
    }
    Ok(())
}
