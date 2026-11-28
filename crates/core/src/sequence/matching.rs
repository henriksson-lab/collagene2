//! One bounded sequence-matching kernel for exact PCR and document searches.
use super::analysis::base_mask;
use collagene_model::{Alphabet, Error, Sequence};

/// Work includes source/pattern classification, prefix initialization and symbol comparisons.
pub(super) struct Work {
    used: usize,
    limit: usize,
}
impl Work {
    pub(super) fn new(limit: usize) -> Self {
        Self { used: 0, limit }
    }
    pub(super) fn resume(used: usize, limit: usize) -> Self {
        Self { used, limit }
    }
    pub(super) fn remaining(&self) -> usize {
        self.limit.saturating_sub(self.used)
    }
    pub(super) fn used(&self) -> usize {
        self.used
    }
    pub(super) fn take(&mut self, amount: usize) -> bool {
        if amount > self.limit.saturating_sub(self.used) {
            self.used = self.limit;
            return false;
        }
        self.used += amount;
        true
    }
}
#[derive(Clone, Copy)]
pub(super) struct Placement {
    pub(super) start: usize,
    pub(super) certain: bool,
}
fn symbol(base: u8, protein: bool) -> u8 {
    if protein { base } else { base_mask(base) }
}

/// Returns false when the shared work bound or the caller's output bound stops the scan.
/// A true result visited every permitted one-turn start. No partial result is called complete.
pub(super) fn scan(
    source: &Sequence,
    pattern: &Sequence,
    circular: bool,
    work: &mut Work,
    mut emit: impl FnMut(Placement) -> bool,
) -> Result<bool, Error> {
    if source.alphabet() != pattern.alphabet() {
        return Err(Error::Unsupported("search alphabets differ".into()));
    }
    let text = source.as_str().as_bytes();
    let query = pattern.as_str().as_bytes();
    let (n, m) = (text.len(), query.len());
    if m > n {
        return Ok(true);
    }
    // Sequence construction guarantees nonempty, bounded validated input.
    let protein = source.alphabet() == Alphabet::Protein;
    let mut concrete = true;
    if !protein {
        for &base in query.iter().chain(text) {
            if !work.take(1) {
                return Ok(false);
            }
            if base_mask(base).count_ones() != 1 {
                concrete = false;
                break;
            }
        }
    }
    if concrete {
        // A u32 prefix entry covers every permitted 10-million-base pattern on native/WASM.
        // Charge initialization before allocating; neither an adversarial prefix nor a
        // circular late mismatch can restore the former O(n*m) comparison behavior.
        if !work.take(m) {
            return Ok(false);
        }
        let mut prefix = Vec::<u32>::new();
        prefix.try_reserve_exact(m).map_err(|_| Error::TooLarge)?;
        prefix.resize(m, 0);
        for i in 1..m {
            let mut j = prefix[i - 1] as usize;
            loop {
                if !work.take(1) {
                    return Ok(false);
                }
                if symbol(query[i], protein) == symbol(query[j], protein) {
                    j += 1;
                    break;
                }
                if j == 0 {
                    break;
                }
                j = prefix[j - 1] as usize;
            }
            prefix[i] = u32::try_from(j).map_err(|_| Error::TooLarge)?;
        }
        let length = if circular {
            n.checked_add(m - 1).ok_or(Error::TooLarge)?
        } else {
            n
        };
        let mut matched = 0usize;
        for i in 0..length {
            let base = symbol(text[i % n], protein);
            loop {
                if !work.take(1) {
                    return Ok(false);
                }
                if base == symbol(query[matched], protein) {
                    matched += 1;
                    break;
                }
                if matched == 0 {
                    break;
                }
                matched = prefix[matched - 1] as usize;
            }
            if matched == m {
                if !emit(Placement {
                    start: i + 1 - m,
                    certain: true,
                }) {
                    return Ok(false);
                }
                matched = prefix[matched - 1] as usize;
            }
        }
    } else {
        // Compatible IUPAC matching is not transitive, so ordinary KMP failure
        // links would omit real hits. Keep the exact compatibility/subset test,
        // with every comparison charged against the same finite job budget.
        let starts = if circular { n } else { n - m + 1 };
        for start in 0..starts {
            let mut certain = true;
            let mut matched = true;
            for (j, &base) in query.iter().enumerate() {
                if !work.take(1) {
                    return Ok(false);
                }
                let a = base_mask(text[(start + j) % n]);
                let b = base_mask(base);
                if a & b == 0 {
                    matched = false;
                    break;
                }
                certain &= a & !b == 0;
            }
            if matched && !emit(Placement { start, certain }) {
                return Ok(false);
            }
        }
    }
    Ok(true)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn concrete_late_mismatch_proves_absence_with_linear_work() {
        let source = Sequence::dna(&"A".repeat(200_000)).unwrap();
        let pattern = Sequence::dna(&format!("{}C", "A".repeat(39_999))).unwrap();
        let mut work = Work::new(8 * (source.len() + pattern.len()));
        let mut hits = Vec::new();
        assert!(
            scan(&source, &pattern, true, &mut work, |p| {
                hits.push(p.start);
                true
            })
            .unwrap()
        );
        assert!(hits.is_empty());
        assert!(work.used() <= 8 * (source.len() + pattern.len()));
    }

    #[test]
    fn compatible_fallback_consumes_shared_budget_without_claiming_absence() {
        let source = Sequence::dna(&"N".repeat(200)).unwrap();
        let pattern = Sequence::dna(&"A".repeat(80)).unwrap();
        let mut work = Work::new(512);
        let mut hits = Vec::new();
        assert!(
            !scan(&source, &pattern, true, &mut work, |p| {
                hits.push((p.start, p.certain));
                true
            })
            .unwrap()
        );
        assert_eq!(work.used(), 512);
        assert_eq!(
            hits,
            vec![(0, false), (1, false), (2, false), (3, false), (4, false)]
        );
        // Another scan cannot reset the same job's exhausted work allowance.
        assert!(
            !scan(&source, &pattern, true, &mut work, |_| panic!(
                "exhausted scan emitted a hit"
            ))
            .unwrap()
        );
        assert_eq!(work.used(), 512);
    }
}
