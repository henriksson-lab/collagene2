//! Incremental suffix-extension recurrence, O(maximum oligo length) per endpoint.
use super::*;

#[derive(Clone, Copy)]
struct Summary {
    length: usize,
    minimum_mismatches: usize,
    maximum_mismatches: usize,
    minimum_score: i32,
    maximum_score: i32,
    prime_certain: usize,
    prime_possible: usize,
    certain: bool,
}
fn minimum(oligo: &Oligo, options: &Options) -> usize {
    options.minimum_annealing_length.max(
        (oligo.ordered.len() - oligo.declared_tail_length)
            .saturating_sub(options.maximum_extra_unpaired),
    )
}
pub(super) fn total_endpoints(
    source: &Sequence,
    topology: Topology,
    oligo: &Oligo,
    options: &Options,
) -> u64 {
    let minimum = minimum(oligo, options);
    if minimum > oligo.ordered.len() - oligo.declared_tail_length || minimum > source.len() {
        return 0;
    }
    match topology {
        Topology::Circular => source.len() as u64 * 2,
        Topology::Linear => (source.len() - minimum + 1) as u64 * 2,
    }
}
fn source_position(
    endpoint: usize,
    offset_from_three_prime: usize,
    n: usize,
    orientation: Orientation,
) -> usize {
    match orientation {
        Orientation::Forward => (endpoint + n - offset_from_three_prime) % n,
        Orientation::Reverse => (endpoint + offset_from_three_prime) % n,
    }
}
fn template_base(source: u8, orientation: Orientation) -> u8 {
    match orientation {
        Orientation::Forward => complement(source),
        Orientation::Reverse => source,
    }
}
fn scores(
    length: usize,
    minimum_mismatches: usize,
    maximum_mismatches: usize,
    options: &Options,
) -> (i32, i32) {
    let score =
        |m: usize| (length - m) as i32 * options.match_score + m as i32 * options.mismatch_score;
    (score(maximum_mismatches), score(minimum_mismatches))
}
fn summaries(
    source: &Sequence,
    oligo: &Oligo,
    options: &Options,
    orientation: Orientation,
    endpoint: usize,
    maximum: usize,
) -> (Vec<Summary>, bool) {
    let mut candidates = Vec::new();
    let (mut low, mut high, mut certain_run, mut possible_run) = (0, 0, 0, 0);
    let (mut certain_open, mut possible_open, mut concrete) = (true, true, true);
    let bases = source.as_str().as_bytes();
    let ordered = oligo.ordered.as_str().as_bytes();
    for offset in 0..maximum {
        let position = source_position(endpoint, offset, bases.len(), orientation);
        let primer = ordered[ordered.len() - 1 - offset];
        let physical = template_base(bases[position], orientation);
        concrete &= mask(primer).count_ones() == 1 && mask(physical).count_ones() == 1;
        let relation = pairing(primer, physical);
        low += usize::from(relation == Pairing::Mismatch);
        high += usize::from(relation != Pairing::Certain);
        certain_open &= relation == Pairing::Certain;
        possible_open &= relation != Pairing::Mismatch;
        certain_run += usize::from(certain_open);
        possible_run += usize::from(possible_open);
        let length = offset + 1;
        let (score_min, score_max) = scores(length, low, high, options);
        if length < minimum(oligo, options)
            || low > options.maximum_mismatches
            || possible_run < options.minimum_three_prime_pairing
            || score_max < options.minimum_score
        {
            continue;
        }
        candidates.push(Summary {
            length,
            minimum_mismatches: low,
            maximum_mismatches: high,
            minimum_score: score_min,
            maximum_score: score_max,
            prime_certain: certain_run,
            prime_possible: possible_run,
            certain: high <= options.maximum_mismatches
                && certain_run >= options.minimum_three_prime_pairing
                && score_min >= options.minimum_score,
        });
    }
    // A lower bound from a merely possible competitor cannot safely eliminate
    // another alignment: that competitor might fail the mismatch/anchor policy.
    let guaranteed = candidates
        .iter()
        .filter(|s| s.certain)
        .map(|s| s.minimum_score)
        .max();
    candidates.retain(|candidate| guaranteed.is_none_or(|floor| candidate.maximum_score >= floor));
    (candidates, concrete)
}
fn build_hit(
    source: &Sequence,
    oligo: &Oligo,
    oligo_index: usize,
    orientation: Orientation,
    endpoint: usize,
    summary: Summary,
    concrete: bool,
) -> Hit {
    let n = source.len();
    let m = oligo.ordered.len();
    let length = summary.length;
    let start = match orientation {
        Orientation::Forward => (endpoint + n + 1 - length) % n,
        Orientation::Reverse => endpoint,
    };
    let bases = source.as_str().as_bytes();
    let ordered = oligo.ordered.as_str().as_bytes();
    let mut columns = Vec::with_capacity(length);
    for i in 0..length {
        let position = source_position(endpoint, length - 1 - i, n, orientation);
        let primer = ordered[m - length + i];
        let physical = template_base(bases[position], orientation);
        columns.push(Column {
            oligo_index: m - length + i,
            source_position: position,
            oligo_base: primer as char,
            source_base: bases[position] as char,
            template_base: physical as char,
            pairing: pairing(primer, physical),
        });
    }
    Hit {
        oligo_index,
        oligo_id: oligo.id.clone(),
        orientation,
        three_prime_source_position: endpoint,
        source_start: start,
        length,
        crosses_origin: start + length > n,
        declared_tail_length: oligo.declared_tail_length,
        extra_unpaired_length: m - length - oligo.declared_tail_length,
        annealing_start: m - length,
        minimum_mismatches: summary.minimum_mismatches,
        maximum_mismatches: summary.maximum_mismatches,
        minimum_score: summary.minimum_score,
        maximum_score: summary.maximum_score,
        three_prime_certain: summary.prime_certain,
        three_prime_possible: summary.prime_possible,
        certainly_admissible: summary.certain,
        score_scope: if concrete {
            ScoreScope::ExactConcreteOptimum
        } else {
            ScoreScope::ConservativeAmbiguousEnvelope
        },
        columns,
    }
}
/// Retain every best-score concrete tie, or the conservative ambiguity envelope,
/// independently at each physical 3′ endpoint/orientation. No gaps are supported.
/// Work/hit/charge exhaustion is explicit per oligo; later oligos stay Unsearched.
pub fn search(
    source: &Sequence,
    topology: Topology,
    collection: &Collection,
    options: &Options,
    budget: &mut Budget,
) -> Result<Report, Error> {
    let admission = Budget::admission_cost(source, collection)?;
    options.validate()?;
    let before = budget.used;
    if admission.work > budget.limits.work - budget.used.work
        || admission.output_charge > budget.limits.output_charge - budget.used.output_charge
    {
        // No report can be constructed honestly inside the caller's input/metadata
        // allowance. Error is explicit; it is never an empty completed scan.
        return Err(Error::TooLarge);
    }
    budget.work(admission.work);
    budget
        .output(0, admission.output_charge)
        .map_err(|_| Error::TooLarge)?;
    let mut report = Report {
        engine_version: ENGINE_VERSION.into(),
        model: MODEL.into(),
        source_content_hash: identity::source(source, topology),
        collection_hash: identity::collection(collection),
        options_hash: identity::options(options),
        options: options.clone(),
        budget_limits: budget.limits,
        budget_before: before,
        budget_after: before,
        complete: true,
        oligos: collection
            .oligos
            .iter()
            .map(|o| OligoProgress {
                oligo_id: o.id.clone(),
                status: Status::Unsearched,
                endpoints_total: total_endpoints(source, topology, o, options),
                endpoints_examined: 0,
                alignments_found: 0,
                alignments_retained: 0,
            })
            .collect(),
        hits: vec![],
    };
    'collection: for (oligo_index, oligo) in collection.oligos.iter().enumerate() {
        let minimum = minimum(oligo, options);
        let domain = oligo.ordered.len() - oligo.declared_tail_length;
        if report.oligos[oligo_index].endpoints_total == 0 {
            report.oligos[oligo_index].status = Status::Complete;
            continue;
        }
        for orientation in [Orientation::Forward, Orientation::Reverse] {
            let (start, end) = match (topology, orientation) {
                (Topology::Circular, _) => (0, source.len()),
                (Topology::Linear, Orientation::Forward) => (minimum - 1, source.len()),
                (Topology::Linear, Orientation::Reverse) => (0, source.len() - minimum + 1),
            };
            for endpoint in start..end {
                let available = match (topology, orientation) {
                    (Topology::Circular, _) => source.len(),
                    (Topology::Linear, Orientation::Forward) => endpoint + 1,
                    (Topology::Linear, Orientation::Reverse) => source.len() - endpoint,
                };
                let maximum = domain.min(available);
                // Reserve the entire endpoint before scoring. Never publish a short
                // suffix as optimal merely because extension work ran out halfway.
                if !budget.work(1 + maximum as u64) {
                    report.oligos[oligo_index].status = Status::Limited(Limit::Work);
                    report.complete = false;
                    break 'collection;
                }
                let (candidates, concrete) =
                    summaries(source, oligo, options, orientation, endpoint, maximum);
                let progress = &mut report.oligos[oligo_index];
                progress.endpoints_examined += 1;
                progress.alignments_found += candidates.len() as u64;
                // Conservative fixed allocation charge includes column and string/
                // container overhead. Preflight ALL ties before any result clone.
                let charge = candidates
                    .iter()
                    .map(|s| 1024 + oligo.id.len() + s.length * 256)
                    .sum();
                if let Err(limit) = budget.output(candidates.len(), charge) {
                    progress.status = Status::Limited(limit);
                    report.complete = false;
                    break 'collection;
                }
                progress.alignments_retained += candidates.len();
                for candidate in candidates {
                    report.hits.push(build_hit(
                        source,
                        oligo,
                        oligo_index,
                        orientation,
                        endpoint,
                        candidate,
                        concrete,
                    ));
                }
            }
        }
        report.oligos[oligo_index].status = Status::Complete;
    }
    report.budget_after = budget.used;
    Ok(report)
}
impl Hit {
    /// Check exact recorded footprint/columns and admissibility independently of
    /// the extension recurrence. Does NOT certify optimality, ties or completeness.
    pub fn validate_columns(
        &self,
        source: &Sequence,
        topology: Topology,
        oligo: &Oligo,
        options: &Options,
    ) -> Result<(), Error> {
        check_source(source)?;
        oligo.validate()?;
        options.validate()?;
        let m = oligo.ordered.len();
        let n = source.len();
        if self.length < minimum(oligo, options)
            || self.length > m - oligo.declared_tail_length
            || self.length > n
            || self.columns.len() != self.length
            || self.three_prime_source_position >= n
            || self.oligo_id != oligo.id
            || self.annealing_start != m - self.length
            || self.declared_tail_length != oligo.declared_tail_length
            || self.extra_unpaired_length != self.annealing_start - self.declared_tail_length
        {
            return Err(invalid("invalid hit footprint/oligo bounds"));
        }
        let expected_start = match self.orientation {
            Orientation::Forward => (self.three_prime_source_position + n + 1 - self.length) % n,
            Orientation::Reverse => self.three_prime_source_position,
        };
        if self.source_start != expected_start
            || self.crosses_origin != (self.source_start + self.length > n)
            || topology == Topology::Linear && self.crosses_origin
        {
            return Err(invalid("invalid strand/circular footprint"));
        }
        let bases = source.as_str().as_bytes();
        let ordered = oligo.ordered.as_str().as_bytes();
        let (mut low, mut high) = (0, 0);
        for (i, column) in self.columns.iter().enumerate() {
            let position = match self.orientation {
                Orientation::Forward => (self.source_start + i) % n,
                Orientation::Reverse => (self.source_start + self.length - 1 - i) % n,
            };
            let physical = template_base(bases[position], self.orientation);
            let relation = pairing(ordered[self.annealing_start + i], physical);
            if column.source_position != position
                || column.oligo_index != self.annealing_start + i
                || column.oligo_base != ordered[self.annealing_start + i] as char
                || column.source_base != bases[position] as char
                || column.template_base != physical as char
                || column.pairing != relation
            {
                return Err(invalid("aligned column differs from retained source/oligo"));
            }
            low += usize::from(relation == Pairing::Mismatch);
            high += usize::from(relation != Pairing::Certain);
        }
        let certain = self
            .columns
            .iter()
            .rev()
            .take_while(|c| c.pairing == Pairing::Certain)
            .count();
        let possible = self
            .columns
            .iter()
            .rev()
            .take_while(|c| c.pairing != Pairing::Mismatch)
            .count();
        let (score_min, score_max) = scores(self.length, low, high, options);
        let admissible = high <= options.maximum_mismatches
            && certain >= options.minimum_three_prime_pairing
            && score_min >= options.minimum_score;
        if self.minimum_mismatches != low
            || self.maximum_mismatches != high
            || self.minimum_score != score_min
            || self.maximum_score != score_max
            || self.three_prime_certain != certain
            || self.three_prime_possible != possible
            || self.certainly_admissible != admissible
            || self.score_scope == ScoreScope::ExactConcreteOptimum
                && self.columns.iter().any(|c| {
                    mask(c.oligo_base as u8).count_ones() != 1
                        || mask(c.template_base as u8).count_ones() != 1
                })
            || low > options.maximum_mismatches
            || possible < options.minimum_three_prime_pairing
            || score_max < options.minimum_score
        {
            return Err(invalid("score/mismatch/3′ pairing evidence differs"));
        }
        Ok(())
    }
}
