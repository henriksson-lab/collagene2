//! A bounded index of the actual retained CIGAR, not another alignment format.
use super::{Work, contract::bad};
use crate::alignment::{CigarKind, review};
use collagene_model::Error;
struct ReferenceRun {
    start: usize,
    end: usize,
    query: usize,
    paired: bool,
}
struct Insertion {
    boundary: usize,
    query: usize,
    length: usize,
}
pub(super) struct Index {
    runs: Vec<ReferenceRun>,
    insertions: Vec<Insertion>,
}
impl Index {
    pub fn new(report: &review::Report, work: &mut Work) -> Result<Self, Error> {
        let operations = &report.comparison.alignment.operations;
        if operations.len() > work.limits.max_mapping_work - work.summary.mapping_work
            || operations
                .len()
                .checked_mul(std::mem::size_of::<ReferenceRun>() + std::mem::size_of::<Insertion>())
                .is_none_or(|bytes| bytes > 64 * 1024 * 1024)
        {
            return Err(Error::TooLarge);
        }
        let mut runs = Vec::new();
        let mut insertions = Vec::new();
        runs.try_reserve_exact(operations.len())
            .map_err(|_| Error::TooLarge)?;
        insertions
            .try_reserve_exact(operations.len())
            .map_err(|_| Error::TooLarge)?;
        let (mut r, mut q) = (0usize, 0usize);
        for operation in operations {
            if !work.mapping(1) {
                return Err(Error::TooLarge);
            }
            match operation.kind {
                CigarKind::Insertion => {
                    insertions.push(Insertion {
                        boundary: r,
                        query: q,
                        length: operation.length,
                    });
                    q = q.checked_add(operation.length).ok_or(Error::TooLarge)?;
                }
                CigarKind::Deletion => {
                    let end = r.checked_add(operation.length).ok_or(Error::TooLarge)?;
                    runs.push(ReferenceRun {
                        start: r,
                        end,
                        query: q,
                        paired: false,
                    });
                    r = end;
                }
                _ => {
                    let end = r.checked_add(operation.length).ok_or(Error::TooLarge)?;
                    runs.push(ReferenceRun {
                        start: r,
                        end,
                        query: q,
                        paired: true,
                    });
                    r = end;
                    q = q.checked_add(operation.length).ok_or(Error::TooLarge)?;
                }
            }
        }
        if r != report.reference_length || q != report.query_length {
            return Err(bad("Whole selected path does not cover its exact inputs"));
        }
        Ok(Self { runs, insertions })
    }
    pub fn at(&self, position: usize) -> Option<usize> {
        let index = self
            .runs
            .partition_point(|run| run.start <= position)
            .checked_sub(1)?;
        let run = &self.runs[index];
        (position < run.end && run.paired).then(|| run.query + position - run.start)
    }
    pub fn insertion(&self, boundary: usize) -> Option<(usize, usize)> {
        let index = self
            .insertions
            .partition_point(|run| run.boundary < boundary);
        self.insertions
            .get(index)
            .filter(|run| run.boundary == boundary)
            .map(|run| (run.query, run.length))
    }
}
