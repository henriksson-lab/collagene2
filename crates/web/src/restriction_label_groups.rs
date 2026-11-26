//! Bounded screen-density grouping of actual restriction-site anchors.
//! A group never exceeds3 screen pixels of total contour span: neighboring
//! distances cannot chain an entire dense plasmid into one invented site.
const MAX_SITES: usize = 10_000;
const SPAN: f64 = 3.;
#[derive(Clone, Copy, Debug)]
pub(crate) struct Anchor {
    pub row: usize,
    pub boundary: u64,
}
pub(crate) fn group(
    anchors: &[Anchor],
    length: u64,
    backbone_radius: f64,
) -> Result<Vec<Vec<usize>>, ()> {
    if anchors.len() > MAX_SITES
        || length == 0
        || !backbone_radius.is_finite()
        || backbone_radius <= 0.
        || backbone_radius > 1_000_000_000.
        || anchors.iter().any(|a| a.boundary >= length)
        || anchors.windows(2).any(|p| p[0].boundary > p[1].boundary)
    {
        return Err(());
    }
    let circumference = std::f64::consts::TAU * backbone_radius;
    let close = |bases: u64| bases as f64 / length as f64 * circumference <= SPAN;
    let mut groups = Vec::<Vec<usize>>::new();
    let mut first = 0u64;
    for anchor in anchors {
        if groups.is_empty() || !close(anchor.boundary - first) {
            first = anchor.boundary;
            groups.push(Vec::new());
        }
        groups.last_mut().ok_or(())?.push(anchor.row);
    }
    // Only the two origin-adjacent groups can merge. The whole merged contour
    // span is tested, including their extremities, rather than one nearest pair.
    if groups.len() > 1 {
        let last_count = groups.last().ok_or(())?.len();
        let first_count = groups[0].len();
        let last_start = anchors[anchors.len() - last_count].boundary;
        let first_end = anchors[first_count - 1].boundary;
        if close(length - last_start + first_end) {
            let last = groups.pop().ok_or(())?;
            groups[0].extend(last);
        }
    }
    Ok(groups)
}
pub(crate) fn budget(height: f64) -> usize {
    if !height.is_finite() {
        return 16;
    }
    (2. * (height - 16.).max(0.) / 16.).floor().clamp(16., 128.) as usize
}
/// Deterministic contour-stratified sampling with an explicit selected group.
/// It chooses only presentation candidates and never removes actual site marks.
pub(crate) fn candidates(total: usize, selected: Option<usize>, budget: usize) -> Vec<usize> {
    if total > MAX_SITES {
        return Vec::new();
    }
    let budget = budget.clamp(1, 128);
    if total <= budget {
        return (0..total).collect();
    }
    let selected = selected.filter(|index| *index < total);
    let available = (0..total)
        .filter(|index| Some(*index) != selected)
        .collect::<Vec<_>>();
    let slots = budget - usize::from(selected.is_some());
    let mut out = Vec::with_capacity(budget);
    out.extend(selected);
    for i in 0..slots {
        let start = i * available.len() / slots;
        let end = (i + 1) * available.len() / slots;
        out.push(available[start + (end - start) / 2]);
    }
    out.sort_unstable();
    out
}
/// Keep already-visible enzyme names fixed under the pointer. Only an otherwise
/// hidden selected member replaces the last named entry of a compact variant.
/// The caller supplies canonical unique-name representatives, with any selected
/// same-name member substituted in place so its exact identity is retained.
pub(crate) fn visible_members(
    canonical: &[usize],
    selected: Option<usize>,
    count: usize,
) -> Vec<usize> {
    if canonical.len() > MAX_SITES {
        return Vec::new();
    }
    let mut visible = canonical
        .iter()
        .copied()
        .take(count.min(3))
        .collect::<Vec<_>>();
    if let Some(selected) = selected {
        if canonical.contains(&selected) && !visible.contains(&selected) {
            if let Some(last) = visible.last_mut() {
                *last = selected;
            }
        }
    }
    visible
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn nearby_anchors_use_total_span_not_transitive_neighbor_chaining() {
        // Circumference1000px: boundaries are also pixel offsets.
        let radius = 1000. / std::f64::consts::TAU;
        let anchors = (0..10)
            .map(|i| Anchor {
                row: i,
                boundary: (100 + 2 * i) as u64,
            })
            .collect::<Vec<_>>();
        let groups = group(&anchors, 1000, radius).unwrap();
        assert_eq!(
            groups,
            vec![vec![0, 1], vec![2, 3], vec![4, 5], vec![6, 7], vec![8, 9]]
        );
        assert_eq!(groups.iter().flatten().count(), anchors.len());
    }
    #[test]
    fn origin_neighbors_merge_without_moving_or_losing_member_anchors() {
        let radius = 1000. / std::f64::consts::TAU;
        let anchors = [
            Anchor {
                row: 9,
                boundary: 0,
            },
            Anchor {
                row: 4,
                boundary: 1,
            },
            Anchor {
                row: 6,
                boundary: 500,
            },
            Anchor {
                row: 2,
                boundary: 999,
            },
        ];
        assert_eq!(
            group(&anchors, 1000, radius).unwrap(),
            vec![vec![9, 4, 2], vec![6]]
        );
        let dispersed = [
            Anchor {
                row: 9,
                boundary: 0,
            },
            Anchor {
                row: 4,
                boundary: 2,
            },
            Anchor {
                row: 6,
                boundary: 500,
            },
            Anchor {
                row: 2,
                boundary: 998,
            },
        ];
        assert_eq!(
            group(&dispersed, 1000, radius).unwrap(),
            vec![vec![9, 4], vec![6], vec![2]]
        );
        assert!(group(&anchors, 1000, f64::NAN).is_err());
    }
    #[test]
    fn visible_hover_names_stay_in_place_and_hidden_names_promote_only_when_needed() {
        let canonical = [9, 4, 2, 8];
        assert_eq!(visible_members(&canonical, None, 3), vec![9, 4, 2]);
        assert_eq!(visible_members(&canonical, Some(4), 3), vec![9, 4, 2]);
        assert_eq!(visible_members(&canonical, Some(2), 3), vec![9, 4, 2]);
        assert_eq!(visible_members(&canonical, Some(8), 3), vec![9, 4, 8]);
        assert_eq!(visible_members(&canonical, Some(4), 2), vec![9, 4]);
        assert_eq!(visible_members(&canonical, Some(2), 2), vec![9, 2]);
        assert_eq!(visible_members(&canonical, Some(8), 1), vec![8]);
        assert_eq!(visible_members(&canonical, Some(99), 3), vec![9, 4, 2]);
    }
    #[test]
    fn responsive_budgets_keep_selected_actual_group_and_distinct_strata() {
        assert_eq!(budget(80.), 16);
        assert_eq!(budget(528.), 64);
        assert_eq!(budget(2000.), 128);
        let chosen = candidates(10_000, Some(9876), 128);
        assert_eq!(chosen.len(), 128);
        assert!(chosen.contains(&9876));
        assert!(chosen.windows(2).all(|pair| pair[0] < pair[1]));
        assert!(chosen[0] < 100);
        assert!(chosen[127] > 9900);
        assert_eq!(candidates(4, Some(2), 16), vec![0, 1, 2, 3]);
    }
}
