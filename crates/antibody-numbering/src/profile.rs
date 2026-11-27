//! Bounded semiglobal affine profile alignment, independently implemented.
//! Scores/paths under an explicit profile are not calibrated antibody identity,
//! numbering, species or CDR assignments. No reference models are bundled here.
use super::*;
use collagene_model::{PROTEIN_RESIDUES, protein_choices};
pub const PROFILE_ENGINE: &str = "antibody-profile-alignment-1";
pub const PROFILE_VERSION: u16 = 1;
pub const MAX_STATES: usize = 512;
pub const MAX_CELLS: usize = 2_000_000;
const MAX_WEIGHT: i16 = 1000;
const M: usize = 0;
const D: usize = 1;
const I: usize = 2;
data! {pub struct MatchState {
    /// Same column order as collagene_model::PROTEIN_RESIDUES. None forbids that
    /// emission; finite weights are integer scores, not probabilities.
    pub emissions:[Option<i16>;22],pub required:bool,pub allow_insert_after:bool,
}}
data! {pub struct Profile {
    pub version:u16,pub name:String,pub provenance:String,pub states:Vec<MatchState>,
    pub gap_open:i16,pub gap_extend:i16,pub minimum_score:i32,
}}
data! {pub struct Limits {pub max_cells:usize}}
data! {pub enum Step {
    Match{source_index:usize,model_index:usize},
    Insert{source_index:usize,after_model:usize},
    Delete{model_index:usize},
}}
data! {pub struct Alignment {
    pub score:i32,pub start:usize,pub end:usize,pub steps:Vec<Step>,
    /// More than one optimal state path for this scoring bound. The stored
    /// representative is deterministic, never a unique-numbering certificate.
    pub alternative_optimal_paths:bool,
}}
data! {pub enum ThresholdStatus {NoPath,BelowThreshold,PossibleScore,GuaranteedScore}}
data! {pub struct Report {
    pub engine:String,pub source_revision:String,pub profile_revision:String,pub context_revision:String,
    pub source_length:usize,pub limits:Limits,pub cells:usize,
    /// Conservative lower bound of the best score across residue assignments.
    pub lower:Option<Alignment>,
    /// Upper bound is the maximum attainable best score over assignments.
    pub upper:Option<Alignment>,
    pub ambiguous_residues:usize,pub threshold_status:ThresholdStatus,
}}
impl Profile {
    pub fn revision(&self) -> Result<String> {
        if self.version != PROFILE_VERSION {
            return Err(Error::Invalid("unsupported profile version"));
        }
        for s in [&self.name, &self.provenance] {
            if s.trim().is_empty() || s.len() > 1024 || s.chars().any(char::is_control) {
                return Err(Error::Invalid("profile label/provenance"));
            }
        }
        if self.states.is_empty()
            || self.states.len() > MAX_STATES
            || !self.states.iter().any(|r| r.required)
        {
            return Err(Error::Invalid(
                "bounded profile with at least one required match state needed",
            ));
        }
        if !(-MAX_WEIGHT..=0).contains(&self.gap_open)
            || !(-MAX_WEIGHT..=0).contains(&self.gap_extend)
        {
            return Err(Error::Invalid("gap scores must be nonpositive and bounded"));
        }
        for r in &self.states {
            if !r.emissions.iter().any(Option::is_some)
                || r.emissions
                    .iter()
                    .flatten()
                    .any(|s| !(-MAX_WEIGHT..=MAX_WEIGHT).contains(s))
            {
                return Err(Error::Invalid("emission scores"));
            }
        }
        hash(
            b"antibody-profile-1\0",
            &(PROFILE_ENGINE, PROTEIN_RESIDUES, self),
            1 << 20,
        )
    }
}
fn prepare(
    source: &Molecule,
    profile: &Profile,
    limits: &Limits,
) -> Result<(String, String, usize)> {
    let source = source_revision(source)?;
    let revision = profile.revision()?;
    if limits.max_cells == 0 || limits.max_cells > MAX_CELLS {
        return Err(Error::Limit("profile cell budget"));
    }
    // Source length is supplied separately below after validating the full model.
    Ok((source, revision, limits.max_cells))
}
fn request(source_hash: &str, profile_hash: &str, limits: &Limits) -> Result<String> {
    hash(
        b"antibody-profile-context-1\0",
        &(
            PROFILE_ENGINE,
            source_hash,
            profile_hash,
            limits,
            "semiglobal source flanks free; whole profile; M/D/I affine; no I-D edges; required states cannot delete; M,D,I deterministic tie order; stop barrier; lower=min-emission upper=max-emission",
        ),
        65536,
    )
}
#[derive(Clone, Copy, Default)]
struct Node {
    score: i32,
    ways: u8,
    parents: u8,
}
impl Node {
    fn add(&mut self, parent: Node, cost: i32, parent_state: usize) {
        if parent.ways == 0 {
            return;
        }
        let score = parent.score + cost;
        if self.ways == 0 || score > self.score {
            *self = Self {
                score,
                ways: parent.ways,
                parents: 1 << parent_state,
            };
        } else if score == self.score {
            self.ways = (self.ways + parent.ways).min(2);
            self.parents |= 1 << parent_state;
        }
    }
}
fn emissions(profile: &Profile) -> Vec<[[Option<i16>; 2]; 128]> {
    profile
        .states
        .iter()
        .map(|r| {
            let mut cache = [[None; 2]; 128];
            for symbol in PROTEIN_RESIDUES
                .iter()
                .copied()
                .chain(b"BZJX".iter().copied())
            {
                let choices = protein_choices(symbol).expect("protein alphabet");
                let mut low = Some(i16::MAX);
                let mut high = None::<i16>;
                for &b in choices {
                    let e = r.emissions[PROTEIN_RESIDUES.iter().position(|&a| a == b).unwrap()];
                    low = low.zip(e).map(|(a, b)| a.min(b));
                    if let Some(e) = e {
                        high = Some(high.map_or(e, |v| v.max(e)));
                    }
                }
                cache[symbol as usize] = [low, high];
            }
            cache
        })
        .collect()
}
fn fill(
    source: &[u8],
    profile: &Profile,
    cache: &[[[Option<i16>; 2]; 128]],
    mode: usize,
    matrix: &mut [[Node; 3]],
    cancel: &mut impl FnMut() -> bool,
) -> Result<()> {
    matrix.fill([Node::default(); 3]);
    let width = profile.states.len() + 1;
    for i in 0..=source.len() {
        matrix[i * width][M] = Node {
            score: 0,
            ways: 1,
            parents: 0,
        };
    }
    for i in 0..=source.len() {
        if i % 16 == 0 && cancel() {
            return Err(Error::Cancelled);
        }
        for j in 1..width {
            let index = i * width + j;
            let state = &profile.states[j - 1];
            let mut cell = [Node::default(); 3];
            if i > 0 {
                if let Some(e) = cache[j - 1][source[i - 1] as usize][mode] {
                    for k in [M, D, I] {
                        cell[M].add(matrix[(i - 1) * width + j - 1][k], i32::from(e), k);
                    }
                }
                if state.allow_insert_after && source[i - 1] != b'*' {
                    cell[I].add(
                        matrix[(i - 1) * width + j][M],
                        i32::from(profile.gap_open),
                        M,
                    );
                    cell[I].add(
                        matrix[(i - 1) * width + j][I],
                        i32::from(profile.gap_extend),
                        I,
                    );
                }
            }
            if !state.required {
                cell[D].add(matrix[index - 1][M], i32::from(profile.gap_open), M);
                cell[D].add(matrix[index - 1][D], i32::from(profile.gap_extend), D);
            }
            matrix[index] = cell;
        }
    }
    Ok(())
}
fn traceback(
    matrix: &[[Node; 3]],
    width: usize,
    end: usize,
    state: usize,
    best: Node,
) -> Option<Alignment> {
    if best.ways == 0 {
        return None;
    }
    let mut i = end;
    let mut j = width - 1;
    let mut k = state;
    let mut steps = vec![];
    while j > 0 {
        let n = matrix[i * width + j][k];
        match k {
            M => {
                steps.push(Step::Match {
                    source_index: i - 1,
                    model_index: j - 1,
                });
                i -= 1;
                j -= 1;
            }
            D => {
                steps.push(Step::Delete { model_index: j - 1 });
                j -= 1;
            }
            I => {
                steps.push(Step::Insert {
                    source_index: i - 1,
                    after_model: j - 1,
                });
                i -= 1;
            }
            _ => unreachable!(),
        }
        k = n.parents.trailing_zeros() as usize;
    }
    steps.reverse();
    Some(Alignment {
        score: best.score,
        start: i,
        end,
        steps,
        alternative_optimal_paths: best.ways > 1,
    })
}
fn endpoint(matrix: &[[Node; 3]], width: usize, end: usize) -> (Node, usize) {
    let mut best = Node::default();
    let mut state = M;
    for k in [M, D, I] {
        let n = matrix[end * width + width - 1][k];
        if n.ways == 0 {
            continue;
        }
        if best.ways == 0 || n.score > best.score {
            best = n;
            state = k;
        } else if n.score == best.score {
            best.ways = (best.ways + n.ways).min(2);
        }
    }
    (best, state)
}
fn solve(
    source: &[u8],
    profile: &Profile,
    cache: &[[[Option<i16>; 2]; 128]],
    mode: usize,
    matrix: &mut [[Node; 3]],
    cancel: &mut impl FnMut() -> bool,
) -> Result<Option<Alignment>> {
    fill(source, profile, cache, mode, matrix, cancel)?;
    let width = profile.states.len() + 1;
    let mut best = Node::default();
    let mut selected = (0, M);
    for end in 0..=source.len() {
        let (node, state) = endpoint(matrix, width, end);
        if node.ways == 0 {
            continue;
        }
        if best.ways == 0 || node.score > best.score {
            best = node;
            selected = (end, state);
        } else if node.score == best.score {
            best.ways = (best.ways + node.ways).min(2);
        }
    }
    Ok(traceback(matrix, width, selected.0, selected.1, best))
}
pub fn align(
    source: &Molecule,
    profile: &Profile,
    limits: &Limits,
    cancel: &mut impl FnMut() -> bool,
) -> Result<Report> {
    if cancel() {
        return Err(Error::Cancelled);
    }
    let (source_revision, profile_revision, budget) = prepare(source, profile, limits)?;
    let cells = (source.sequence.len() + 1)
        .checked_mul(profile.states.len() + 1)
        .filter(|&n| n <= budget)
        .ok_or(Error::Limit("profile cells"))?;
    let context_revision = request(&source_revision, &profile_revision, limits)?;
    let cache = emissions(profile);
    let mut matrix = vec![[Node::default(); 3]; cells];
    let bytes = source.sequence.as_str().as_bytes();
    let lower = solve(bytes, profile, &cache, 0, &mut matrix, cancel)?;
    let upper = solve(bytes, profile, &cache, 1, &mut matrix, cancel)?;
    if cancel() {
        return Err(Error::Cancelled);
    }
    let threshold_status = if upper.is_none() {
        ThresholdStatus::NoPath
    } else if upper.as_ref().unwrap().score < profile.minimum_score {
        ThresholdStatus::BelowThreshold
    } else if lower
        .as_ref()
        .is_some_and(|l| l.score >= profile.minimum_score)
    {
        ThresholdStatus::GuaranteedScore
    } else {
        ThresholdStatus::PossibleScore
    };
    Ok(Report {
        engine: PROFILE_ENGINE.into(),
        source_revision,
        profile_revision,
        context_revision,
        source_length: bytes.len(),
        limits: limits.clone(),
        cells,
        lower,
        upper,
        ambiguous_residues: bytes
            .iter()
            .filter(|&&b| matches!(b, b'B' | b'Z' | b'J' | b'X'))
            .count(),
        threshold_status,
    })
}
impl Report {
    /// Reconstruct score bounds, representative paths and ties. This verifies
    /// only the declared score model, never its biological calibration.
    pub fn validate_for(
        &self,
        source: &Molecule,
        profile: &Profile,
        limits: &Limits,
    ) -> Result<()> {
        if self.engine != PROFILE_ENGINE
            || self
                .lower
                .as_ref()
                .is_some_and(|a| a.steps.len() > MAX_SOURCE_RESIDUES + MAX_STATES)
            || self
                .upper
                .as_ref()
                .is_some_and(|a| a.steps.len() > MAX_SOURCE_RESIDUES + MAX_STATES)
            || align(source, profile, limits, &mut || false)? != *self
        {
            return Err(Error::Binding);
        }
        Ok(())
    }
}

// Traverse the optimal predecessor DAG, not an exponential list of paths.
// Intersection of model/source matches is equivalent to one source value and
// no deletion among every reachable optimal consumption of that model state.
fn invariant_matches(
    matrix: &[[Node; 3]],
    width: usize,
    end: usize,
    score: i32,
    visited: &mut [u8],
    remaining: &mut usize,
    cancel: &mut impl FnMut() -> bool,
) -> Result<Vec<InvariantMatch>> {
    fn push(
        cell: usize,
        state: usize,
        visited: &mut [u8],
        remaining: &mut usize,
        stack: &mut Vec<(usize, usize)>,
    ) -> Result<()> {
        let bit = 1_u8 << state;
        if visited[cell] & bit == 0 {
            if *remaining == 0 {
                return Err(Error::Limit("profile invariant trace states"));
            }
            *remaining -= 1;
            visited[cell] |= bit;
            stack.push((cell, state));
        }
        Ok(())
    }
    visited.fill(0);
    let mut stack = vec![];
    for state in [M, D, I] {
        let cell = end * width + width - 1;
        let node = matrix[cell][state];
        if node.ways > 0 && node.score == score {
            push(cell, state, visited, remaining, &mut stack)?;
        }
    }
    let mut sources = vec![None; width - 1];
    let mut conflicts = vec![false; width - 1];
    let mut deleted = vec![false; width - 1];
    let mut processed = 0;
    while let Some((cell, state)) = stack.pop() {
        if processed % 256 == 0 && cancel() {
            return Err(Error::Cancelled);
        }
        processed += 1;
        let i = cell / width;
        let j = cell % width;
        if j == 0 {
            continue;
        }
        let predecessor = match state {
            M => {
                let index = i - 1;
                match sources[j - 1] {
                    None => sources[j - 1] = Some(index),
                    Some(other) if other != index => conflicts[j - 1] = true,
                    _ => (),
                };
                cell - width - 1
            }
            D => {
                deleted[j - 1] = true;
                cell - 1
            }
            I => cell - width,
            _ => unreachable!(),
        };
        let parents = matrix[cell][state].parents;
        for parent in [M, D, I] {
            if parents & (1 << parent) != 0 {
                push(predecessor, parent, visited, remaining, &mut stack)?;
            }
        }
    }
    Ok(sources
        .into_iter()
        .enumerate()
        .filter_map(|(model_index, index)| {
            index
                .filter(|_| !conflicts[model_index] && !deleted[model_index])
                .map(|source_index| InvariantMatch {
                    model_index,
                    source_index,
                })
        })
        .collect())
}

pub const ENDPOINT_ENGINE: &str = "antibody-profile-endpoints-2";
pub const MAX_ENDPOINTS: usize = 128;
pub const MAX_TRACE_STATES: usize = 2_000_000;
data! {pub struct EndpointLimits {pub max_cells:usize,pub max_endpoints:usize,pub max_trace_states:usize}}
data! {pub struct InvariantMatch {pub model_index:usize,pub source_index:usize}}
data! {pub struct Endpoint {
    /// Zero-based exclusive source endpoint; each endpoint has its own optimum.
    pub end:usize,pub lower:Option<Alignment>,pub upper:Alignment,pub threshold_status:ThresholdStatus,
    /// Matches shared by every optimal path. None for residue-ambiguous sources.
    pub invariant_matches:Option<Vec<InvariantMatch>>,
}}
data! {pub struct EndpointReport {
    pub engine:String,pub source_revision:String,pub profile_revision:String,pub context_revision:String,
    pub source_length:usize,pub limits:EndpointLimits,pub cells:usize,pub invariant_trace_states:usize,
    pub reachable_endpoints:usize,pub qualified_endpoints:usize,pub truncated:bool,
    pub ambiguous_residues:usize,pub endpoints:Vec<Endpoint>,
}}
/// Every full-profile endpoint with an attainable score meeting the declared
/// threshold. Retention is bounded and ordered by source endpoint, not score.
/// Multiple starts/state paths at one end remain tied; they are not enumerated.
pub fn align_endpoints(
    source: &Molecule,
    profile: &Profile,
    limits: &EndpointLimits,
    cancel: &mut impl FnMut() -> bool,
) -> Result<EndpointReport> {
    if cancel() {
        return Err(Error::Cancelled);
    }
    if limits.max_endpoints > MAX_ENDPOINTS || limits.max_trace_states > MAX_TRACE_STATES {
        return Err(Error::Limit("profile endpoint retention"));
    }
    let (source_revision, profile_revision, budget) = prepare(
        source,
        profile,
        &Limits {
            max_cells: limits.max_cells,
        },
    )?;
    let cells = (source.sequence.len() + 1)
        .checked_mul(profile.states.len() + 1)
        .filter(|&n| n <= budget)
        .ok_or(Error::Limit("profile cells"))?;
    let context_revision = hash(
        b"antibody-profile-endpoints-context-2\0",
        &(
            ENDPOINT_ENGINE,
            PROFILE_ENGINE,
            &source_revision,
            &profile_revision,
            limits,
            "all-optimal-path invariant matches for resolved source only; bounded DAG state union; per-source-end optimum; upper score threshold; source-end order; exact counts beyond retention; ties retained as flags",
        ),
        65536,
    )?;
    let bytes = source.sequence.as_str().as_bytes();
    let width = profile.states.len() + 1;
    let cache = emissions(profile);
    let mut matrix = vec![[Node::default(); 3]; cells];
    fill(bytes, profile, &cache, 1, &mut matrix, cancel)?;
    let mut reachable_endpoints = 0;
    let mut qualified_endpoints = 0;
    let ambiguous_residues = bytes
        .iter()
        .filter(|&&b| matches!(b, b'B' | b'Z' | b'J' | b'X'))
        .count();
    let mut remaining_trace = limits.max_trace_states;
    let mut visited = vec![0_u8; cells];
    let mut endpoints = vec![];
    for end in 0..=bytes.len() {
        if end % 16 == 0 && cancel() {
            return Err(Error::Cancelled);
        }
        let (best, state) = endpoint(&matrix, width, end);
        if best.ways == 0 {
            continue;
        }
        reachable_endpoints += 1;
        if best.score < profile.minimum_score {
            continue;
        }
        qualified_endpoints += 1;
        if endpoints.len() < limits.max_endpoints {
            endpoints.push(Endpoint {
                end,
                lower: None,
                upper: traceback(&matrix, width, end, state, best).unwrap(),
                threshold_status: ThresholdStatus::PossibleScore,
                invariant_matches: if ambiguous_residues == 0 {
                    Some(invariant_matches(
                        &matrix,
                        width,
                        end,
                        best.score,
                        &mut visited,
                        &mut remaining_trace,
                        cancel,
                    )?)
                } else {
                    None
                },
            });
        }
    }
    if !endpoints.is_empty() {
        fill(bytes, profile, &cache, 0, &mut matrix, cancel)?;
        for value in &mut endpoints {
            if cancel() {
                return Err(Error::Cancelled);
            }
            let (best, state) = endpoint(&matrix, width, value.end);
            value.lower = traceback(&matrix, width, value.end, state, best);
            if value
                .lower
                .as_ref()
                .is_some_and(|a| a.score >= profile.minimum_score)
            {
                value.threshold_status = ThresholdStatus::GuaranteedScore;
            }
        }
    }
    if cancel() {
        return Err(Error::Cancelled);
    }
    Ok(EndpointReport {
        engine: ENDPOINT_ENGINE.into(),
        source_revision,
        profile_revision,
        context_revision,
        source_length: bytes.len(),
        limits: limits.clone(),
        cells,
        invariant_trace_states: limits.max_trace_states - remaining_trace,
        reachable_endpoints,
        qualified_endpoints,
        truncated: qualified_endpoints > endpoints.len(),
        ambiguous_residues: bytes
            .iter()
            .filter(|&&b| matches!(b, b'B' | b'Z' | b'J' | b'X'))
            .count(),
        endpoints,
    })
}
impl EndpointReport {
    pub fn validate_for(
        &self,
        source: &Molecule,
        profile: &Profile,
        limits: &EndpointLimits,
    ) -> Result<()> {
        if self.engine != ENDPOINT_ENGINE
            || self.endpoints.len() > MAX_ENDPOINTS
            || self.endpoints.iter().any(|e| {
                e.invariant_matches
                    .as_ref()
                    .is_some_and(|m| m.len() > MAX_STATES)
                    || e.upper.steps.len() > MAX_SOURCE_RESIDUES + MAX_STATES
                    || e.lower
                        .as_ref()
                        .is_some_and(|a| a.steps.len() > MAX_SOURCE_RESIDUES + MAX_STATES)
            })
            || *self != align_endpoints(source, profile, limits, &mut || false)?
        {
            return Err(Error::Binding);
        }
        Ok(())
    }
}
