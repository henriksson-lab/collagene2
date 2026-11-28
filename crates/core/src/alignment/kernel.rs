//! One affine recurrence over checked full or endpoint-corridor row storage.
use super::grid::{
    Attempt, CellBudget, Domain, GlobalStorage, Grid, NEG, RowLayout, SearchEvidence, StopReason,
};
use super::{AlignmentMode, AlignmentOptions, MAX_ALIGNMENT_COLUMNS, substitution};
use collagene_model::{Error, Sequence};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(u8)]
pub(super) enum State {
    Pair = 0,
    Insert = 1,
    Delete = 2,
    Stop = 3,
}
#[derive(Debug)]
pub(super) struct Path {
    pub score: i32,
    pub rs: usize,
    pub re: usize,
    pub qs: usize,
    pub qe: usize,
    pub ops: Vec<State>,
}
fn add(score: i32, delta: i32) -> i32 {
    if score == NEG { NEG } else { score + delta }
}
fn best(values: &[(i32, State)]) -> (i32, State) {
    let mut best = values[0];
    for &candidate in &values[1..] {
        if candidate.0 > best.0 {
            best = candidate;
        }
    }
    best
}
fn state(tag: u8) -> State {
    match tag {
        0 => State::Pair,
        1 => State::Insert,
        2 => State::Delete,
        _ => State::Stop,
    }
}
fn invalid(message: &str) -> Error {
    Error::Format(message.into())
}
pub(super) fn next_radius(
    radius: usize,
    maximum: usize,
    n: usize,
    m: usize,
) -> Result<usize, Error> {
    radius
        .checked_add(1)
        .and_then(|one| radius.checked_mul(2).map(|two| one.max(two)))
        .map(|next| next.min(maximum).min(n.min(m)))
        .ok_or(Error::TooLarge)
}
/// One admitted attempt. Every evaluated row charges the same borrowed ledger.
pub(super) fn attempt(
    reference: &Sequence,
    query: &Sequence,
    options: &AlignmentOptions,
    domain: Domain,
    budget: &mut CellBudget<'_>,
) -> Result<Path, Error> {
    let mut path =
        if query.len() > reference.len() && options.mode != AlignmentMode::QueryInReference {
            let mut path = compute(query, reference, options, domain, budget)?;
            std::mem::swap(&mut path.rs, &mut path.qs);
            std::mem::swap(&mut path.re, &mut path.qe);
            for op in &mut path.ops {
                *op = match *op {
                    State::Insert => State::Delete,
                    State::Delete => State::Insert,
                    other => other,
                };
            }
            path
        } else {
            compute(reference, query, options, domain, budget)?
        };
    if path.ops.is_empty() {
        path.rs = 0;
        path.re = 0;
        path.qs = 0;
        path.qe = 0;
    }
    Ok(path)
}
fn compute(
    reference: &Sequence,
    query: &Sequence,
    options: &AlignmentOptions,
    domain: Domain,
    budget: &mut CellBudget<'_>,
) -> Result<Path, Error> {
    let (n, m) = (reference.len(), query.len());
    let layout = RowLayout::new(n, m, domain, budget.remaining())?;
    budget.admit(layout.cells())?;
    let mut grid = Grid::new(layout)?;
    let s = options.scoring;
    let (start, end) = grid.layout().row(0)?;
    budget.charge_row(grid.layout().row_width(0)?)?;
    grid.begin(0)?;
    for j in start..=end {
        let value = if j == 0 {
            [0, NEG, NEG]
        } else if options.mode == AlignmentMode::Local {
            [0, NEG, NEG]
        } else {
            [NEG, s.gap_open + (j as i32 - 1) * s.gap_extend, NEG]
        };
        grid.rows().1.set(j, value)?;
        if j > 0 && options.mode != AlignmentMode::Local {
            grid.set_predecessor(
                0,
                j,
                State::Insert as u8,
                if j == 1 { State::Pair } else { State::Insert } as u8,
            )?;
        }
    }
    grid.finish_row();
    let mut endpoint = if options.mode == AlignmentMode::QueryInReference {
        let row = grid.final_row().get(m);
        let value = best(&[
            (row[0], State::Pair),
            (row[1], State::Insert),
            (row[2], State::Delete),
        ]);
        (value.0, 0, m, value.1)
    } else {
        (0, 0, 0, State::Stop)
    };
    let (rb, qb) = (reference.as_str().as_bytes(), query.as_str().as_bytes());
    for i in 1..=n {
        let (start, end) = grid.layout().row(i)?;
        budget.charge_row(grid.layout().row_width(i)?)?;
        grid.begin(i)?;
        for j in start..=end {
            if j == 0 {
                if options.mode == AlignmentMode::Global {
                    grid.rows()
                        .1
                        .set(0, [NEG, NEG, s.gap_open + (i as i32 - 1) * s.gap_extend])?;
                    grid.set_predecessor(
                        i,
                        0,
                        State::Delete as u8,
                        if i == 1 { State::Pair } else { State::Delete } as u8,
                    )?;
                } else {
                    grid.rows().1.set(0, [0, NEG, NEG])?;
                }
                continue;
            }
            let (previous, current) = grid.rows();
            let diagonal = previous.get(j - 1);
            let left = current.get(j - 1);
            let up = previous.get(j);
            let mut pair = best(&[
                (diagonal[0], State::Pair),
                (diagonal[1], State::Insert),
                (diagonal[2], State::Delete),
            ]);
            pair.0 = add(pair.0, substitution(rb[i - 1], qb[j - 1], s));
            if options.mode == AlignmentMode::Local && pair.0 <= 0 {
                pair = (0, State::Stop);
            }
            // Same native affine model and tie order, including opposite-gap opens.
            let insert = best(&[
                (add(left[1], s.gap_extend), State::Insert),
                (add(left[0], s.gap_open), State::Pair),
                (add(left[2], s.gap_open), State::Delete),
            ]);
            let delete = best(&[
                (add(up[2], s.gap_extend), State::Delete),
                (add(up[0], s.gap_open), State::Pair),
                (add(up[1], s.gap_open), State::Insert),
            ]);
            current.set(j, [pair.0, insert.0, delete.0])?;
            grid.set_predecessor(i, j, State::Pair as u8, pair.1 as u8)?;
            grid.set_predecessor(i, j, State::Insert as u8, insert.1 as u8)?;
            grid.set_predecessor(i, j, State::Delete as u8, delete.1 as u8)?;
            if options.mode == AlignmentMode::Local {
                for (state, score) in [
                    (State::Pair, pair.0),
                    (State::Insert, insert.0),
                    (State::Delete, delete.0),
                ] {
                    if score > endpoint.0 {
                        endpoint = (score, i, j, state);
                    }
                }
            }
        }
        grid.finish_row();
        if options.mode == AlignmentMode::QueryInReference {
            let row = grid.final_row().get(m);
            let value = best(&[
                (row[0], State::Pair),
                (row[1], State::Insert),
                (row[2], State::Delete),
            ]);
            if value.0 > endpoint.0 {
                endpoint = (value.0, i, m, value.1);
            }
        }
    }
    if options.mode == AlignmentMode::Global {
        let row = grid.final_row().get(m);
        let value = best(&[
            (row[0], State::Pair),
            (row[1], State::Insert),
            (row[2], State::Delete),
        ]);
        endpoint = (value.0, n, m, value.1);
    }
    let (score, re, qe, mut current_state) = endpoint;
    if score == NEG {
        return Err(invalid(
            "Admitted global corridor has no reachable endpoint",
        ));
    }
    let (mut i, mut j) = (re, qe);
    let mut ops = Vec::new();
    while current_state != State::Stop {
        if i == 0 && j == 0 {
            break;
        }
        if options.mode == AlignmentMode::QueryInReference && j == 0 {
            break;
        }
        let before = state(grid.predecessor(i, j, current_state as u8)?);
        if options.mode == AlignmentMode::Local
            && current_state == State::Pair
            && before == State::Stop
        {
            break;
        }
        if ops.len() >= MAX_ALIGNMENT_COLUMNS {
            return Err(Error::Unsupported(
                "pairwise traceback exceeds the one million column output budget".into(),
            ));
        }
        match current_state {
            State::Pair if i > 0 && j > 0 => {
                i -= 1;
                j -= 1;
            }
            State::Insert if j > 0 => j -= 1,
            State::Delete if i > 0 => i -= 1,
            _ => {
                return Err(invalid(
                    "pairwise traceback reached an inconsistent boundary",
                ));
            }
        }
        ops.push(current_state);
        current_state = before;
    }
    ops.reverse();
    Ok(Path {
        score,
        rs: i,
        re,
        qs: j,
        qe,
        ops,
    })
}
pub(super) fn align(
    reference: &Sequence,
    query: &Sequence,
    options: &AlignmentOptions,
    budget: &mut CellBudget<'_>,
) -> Result<(Path, SearchEvidence), Error> {
    let mut attempts = Vec::new();
    let mut retained = None;
    let (mut domain, maximum, maximum_attempts) = match options.global_storage {
        GlobalStorage::Full => (Domain::Full, 0, 1),
        GlobalStorage::Adaptive {
            initial_radius,
            maximum_radius,
            max_attempts,
        } => (
            Domain::Corridor {
                radius: initial_radius.min(reference.len().min(query.len())),
            },
            maximum_radius,
            max_attempts,
        ),
    };
    let stop = loop {
        if domain.cell_count(reference.len(), query.len())? > budget.remaining() {
            if retained.is_none() {
                return Err(Error::TooLarge);
            }
            break StopReason::CellLimit;
        }
        let before = budget.used();
        let path = match attempt(reference, query, options, domain, budget) {
            Ok(path) => path,
            Err(Error::TooLarge) if retained.is_some() => break StopReason::MemoryLimit,
            Err(error) => return Err(error),
        };
        let receipt = Attempt::completed(
            domain,
            reference.len(),
            query.len(),
            options.scoring,
            path.score,
            budget.used() - before,
        )?;
        let certified = receipt.unrestricted_optimum_proven;
        attempts.push(receipt);
        retained = Some(path);
        if certified {
            break StopReason::Certified;
        }
        if attempts.len() == maximum_attempts {
            break StopReason::AttemptLimit;
        }
        let Domain::Corridor { radius } = domain else {
            unreachable!()
        };
        let next = next_radius(radius, maximum, reference.len(), query.len())?;
        if next <= radius {
            break StopReason::RadiusLimit;
        }
        domain = Domain::Corridor { radius: next };
    };
    let path = retained.expect("a completed attempt or an error");
    let evidence = SearchEvidence {
        storage: options.global_storage,
        attempts,
        evaluated_cells: budget.used(),
        unrestricted_optimum_proven: stop == StopReason::Certified,
        stop_reason: stop,
    };
    evidence.validate(
        reference.len(),
        query.len(),
        options.mode,
        options.scoring,
        path.score,
        options.max_cells,
    )?;
    Ok((path, evidence))
}
