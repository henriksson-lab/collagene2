//! Pure CSS-pixel layout. These helpers never infer enzyme cleavage or annotation biology.
use std::collections::BTreeMap;

pub const CELL: f64 = 8.5;
pub const LEFT: f64 = 64.;
pub const LABEL_STEP: f64 = 17.;
pub const MAX_LABEL_LANES: usize = 5;

#[derive(Clone, Debug, PartialEq)]
pub struct Label {
    pub index: usize,
    pub row: usize,
    pub anchor: f64,
    pub left: f64,
    pub width: f64,
    pub lane: usize,
}
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Labels {
    pub rows: BTreeMap<usize, Vec<Label>>,
    pub overflow: BTreeMap<usize, usize>,
    pub lanes: usize,
}
impl Labels {
    pub fn height(&self) -> f64 {
        self.lanes as f64 * LABEL_STEP
    }
}

/// Stable whole-source packing; scrolling never changes row height. Dense rows
/// get an explicit count/inspector affordance instead of overlapping labels.
pub fn labels(sites: impl IntoIterator<Item = (usize, usize, usize)>, columns: usize) -> Labels {
    let mut result = Labels::default();
    if columns == 0 {
        return result;
    }
    let available = columns as f64 * CELL;
    let mut sorted: Vec<_> = sites.into_iter().collect();
    sorted.sort_unstable_by_key(|&(index, position, _)| (position, index));
    let mut ends = BTreeMap::<usize, [f64; MAX_LABEL_LANES]>::new();
    for (index, position, characters) in sorted {
        let row = position / columns;
        let anchor = LEFT + (position % columns) as f64 * CELL;
        // Deliberately conservative for 12px system-ui; paint clips at this width.
        let width = (characters as f64 * 8. + 10.).min(available);
        let left = (anchor - width * 0.35).clamp(LEFT, LEFT + available - width);
        let row_ends = ends
            .entry(row)
            .or_insert([f64::NEG_INFINITY; MAX_LABEL_LANES]);
        if let Some(lane) = row_ends.iter().position(|end| *end + 10. <= left) {
            row_ends[lane] = left + width;
            result.lanes = result.lanes.max(lane + 1);
            result.rows.entry(row).or_default().push(Label {
                index,
                row,
                anchor,
                left,
                width,
                lane,
            });
        } else {
            *result.overflow.entry(row).or_default() += 1;
            result.lanes = MAX_LABEL_LANES + 1;
        }
    }
    result
}

#[derive(Debug, PartialEq)]
pub struct FeaturePiece {
    pub left: usize,
    pub right: usize,
    pub head_left: bool,
    pub head_right: bool,
    pub continues_left: bool,
    pub continues_right: bool,
}
pub fn feature_piece(
    a: usize,
    b: usize,
    start: usize,
    end: usize,
    reverse: Option<bool>,
) -> Option<FeaturePiece> {
    let left = a.max(start);
    let right = b.min(end);
    (left < right).then_some(FeaturePiece {
        left,
        right,
        head_left: reverse == Some(true) && left == a,
        head_right: reverse == Some(false) && right == b,
        continues_left: a < start,
        continues_right: b > end,
    })
}

/// Biological 3′ boundary in source coordinates. Complement reverses Join
/// order; unordered/fuzzy endpoints do not acquire a definite arrowhead.
pub fn terminal(location: &collagene_model::Location) -> Option<(usize, bool)> {
    use collagene_model::Location;
    fn walk(location: &Location, reverse: bool, depth: usize) -> Option<(usize, bool)> {
        if depth >= 32 {
            return None;
        }
        match location {
            Location::Range {
                start,
                end,
                before,
                after,
            } => {
                if (reverse && *before) || (!reverse && *after) {
                    None
                } else {
                    Some((if reverse { *start } else { *end }, reverse))
                }
            }
            Location::Complement(inner) => walk(inner, !reverse, depth + 1),
            Location::Join(parts) => walk(
                if reverse {
                    parts.first()?
                } else {
                    parts.last()?
                },
                reverse,
                depth + 1,
            ),
            _ => None,
        }
    }
    walk(location, false, 0)
}

/// Footprint between the actual raw top/bottom interbase cuts. Raw coordinates
/// retain which side of the circular origin the cut belongs to; choosing the
/// shorter arc from normalized endpoints would invent enzyme geometry.
pub fn cut_footprint(top: i64, bottom: i64, length: usize, circular: bool) -> Vec<(usize, usize)> {
    if length == 0 {
        return vec![];
    }
    let (a, b) = (top.min(bottom), top.max(bottom));
    let Some(size) = b.checked_sub(a) else {
        return vec![];
    };
    if size == 0 || size as u64 > length as u64 {
        return vec![];
    }
    if !circular {
        return if a >= 0 && b <= length as i64 {
            vec![(a as usize, b as usize)]
        } else {
            vec![]
        };
    }
    let start = a.rem_euclid(length as i64) as usize;
    let size = size as usize;
    if size <= length - start {
        vec![(start, start + size)]
    } else {
        vec![(start, length), (0, size - (length - start))]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn labels_stagger_without_collisions_or_crossing_row_margins() {
        let packed = labels((0..12).map(|i| (i, 25 + i, 8)), 60);
        assert!(!packed.overflow.is_empty());
        for row in packed.rows.values() {
            for label in row {
                assert!(label.left >= LEFT && label.left + label.width <= LEFT + 60. * CELL);
                for other in row
                    .iter()
                    .filter(|other| other.index != label.index && other.lane == label.lane)
                {
                    assert!(
                        label.left + label.width <= other.left
                            || other.left + other.width <= label.left
                    );
                }
            }
        }
    }
    #[test]
    fn arrow_heads_only_at_directional_range_ends_not_row_clips() {
        let first = feature_piece(10, 140, 0, 60, Some(false)).unwrap();
        assert!(!first.head_right && first.continues_right);
        let last = feature_piece(10, 140, 120, 180, Some(false)).unwrap();
        assert!(last.head_right && last.continues_left);
        let reverse = feature_piece(10, 140, 0, 60, Some(true)).unwrap();
        assert!(reverse.head_left && !reverse.head_right);
        assert!(!feature_piece(10, 140, 0, 60, None).unwrap().head_left);
    }
    #[test]
    fn cut_footprint_preserves_raw_stagger_and_origin() {
        assert_eq!(cut_footprint(2, 6, 100, false), vec![(2, 6)]);
        assert_eq!(cut_footprint(98, 102, 100, true), vec![(98, 100), (0, 2)]);
        assert_eq!(cut_footprint(-2, 2, 100, true), vec![(98, 100), (0, 2)]);
        assert_eq!(cut_footprint(2, 98, 100, true), vec![(2, 98)]);
        assert!(cut_footprint(-2, 2, 100, false).is_empty());
        assert!(cut_footprint(7, 7, 100, false).is_empty());
    }
    #[test]
    fn joined_origin_and_complement_heads_follow_biological_order() {
        use collagene_model::Location;
        let joined = Location::Join(vec![Location::range(90, 100), Location::range(0, 20)]);
        assert_eq!(terminal(&joined), Some((20, false)));
        assert_eq!(
            terminal(&Location::Complement(Box::new(joined))),
            Some((90, true))
        );
        assert_eq!(
            terminal(&Location::Range {
                start: 0,
                end: 20,
                before: false,
                after: true
            }),
            None
        );
    }
}
