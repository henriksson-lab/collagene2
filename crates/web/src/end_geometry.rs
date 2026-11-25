//! Presentation of one declared duplex end, not ligation compatibility or physical inference.
//! Both input strands use their own5′→3′ convention. A bottom protrusion is
//! displayed in reverse order on a left-to-right diagram, never complemented.
use collagene_core::chemistry::{End, FivePrime, Polarity, ThreePrime};

#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum Side {
    Left,
    Right,
}
#[derive(Clone, Debug, PartialEq)]
pub struct Glyph {
    pub column: f64,
    pub symbol: char,
    pub top: bool,
}
#[derive(Clone, Debug, PartialEq)]
pub struct Layout {
    pub width_columns: usize,
    pub top_tip: f64,
    pub bottom_tip: f64,
    pub body_edge: f64,
    pub glyphs: Vec<Glyph>,
    pub top_terminus: &'static str,
    pub bottom_terminus: &'static str,
    pub overhang_length: usize,
    pub abbreviated: bool,
}
fn five(value: FivePrime) -> &'static str {
    match value {
        FivePrime::Phosphate => "5′ P",
        FivePrime::Hydroxyl => "5′ OH",
        FivePrime::Unknown => "5′ ?",
    }
}
fn three(value: ThreePrime) -> &'static str {
    match value {
        ThreePrime::Hydroxyl => "3′ OH",
        ThreePrime::Blocked => "3′ blocked",
        ThreePrime::Unknown => "3′ ?",
    }
}
/// Columns describe only the reported protrusion and a generic paired-body stub.
/// No complementary letters, phosphate treatment or compatibility are synthesized.
/// ASCII letter case is retained; all letters must represent IUPAC DNA bases.
pub fn layout(end: &End, side: Side) -> Result<Layout, String> {
    let n = end.overhang.len();
    // Check length before inspecting or copying any source letters. The display
    // allocates at most24 glyphs even for a maximum-sized declared protrusion.
    if n > collagene_core::model::MAX_SEQUENCE_LEN {
        return Err("End overhang exceeds the core sequence-length bound.".into());
    }
    if (end.polarity == Polarity::Blunt) != (n == 0) {
        return Err(
            "Blunt ends require no overhang; staggered ends require a nonempty overhang.".into(),
        );
    }
    if !end
        .overhang
        .bytes()
        .all(|b| b"ACGTRYSWKMBDHVN".contains(&b.to_ascii_uppercase()))
    {
        return Err("End overhang must contain only ASCII IUPAC DNA letters.".into());
    }
    let top = matches!(
        (side, end.polarity),
        (Side::Left, Polarity::FivePrime) | (Side::Right, Polarity::ThreePrime)
    );
    let abbreviated = n > 24;
    let shown = if abbreviated { 21 } else { n };
    let (protruding, recessed, body_edge, glyph_start) = match side {
        Side::Left => (1., 1. + shown as f64, 7. + shown as f64, 1.),
        Side::Right => (7. + shown as f64, 7., 1., 7.),
    };
    let (top_tip, bottom_tip) = if top {
        (protruding, recessed)
    } else {
        (recessed, protruding)
    };
    let mut glyphs = Vec::with_capacity(shown);
    for i in 0..shown {
        let symbol = if abbreviated && i == 10 {
            '…'
        } else {
            let display_index = if abbreviated && i > 10 {
                n - 10 + i - 11
            } else {
                i
            };
            let source_index = if top {
                display_index
            } else {
                n - 1 - display_index
            };
            end.overhang.as_bytes()[source_index] as char
        };
        glyphs.push(Glyph {
            column: glyph_start + i as f64 + 0.5,
            symbol,
            top,
        });
    }
    let (top_terminus, bottom_terminus) = match side {
        Side::Left => (
            five(end.chemistry.five_prime),
            three(end.chemistry.three_prime),
        ),
        Side::Right => (
            three(end.chemistry.three_prime),
            five(end.chemistry.five_prime),
        ),
    };
    Ok(Layout {
        width_columns: top_tip.max(bottom_tip).max(body_edge).ceil() as usize + 1,
        top_tip,
        bottom_tip,
        body_edge,
        glyphs,
        top_terminus,
        bottom_terminus,
        overhang_length: n,
        abbreviated,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use collagene_core::chemistry::Terminus;
    fn end(polarity: Polarity, letters: &str) -> End {
        End {
            polarity,
            overhang: letters.into(),
            chemistry: Terminus::CLEAVED,
        }
    }
    fn text(layout: &Layout) -> String {
        layout.glyphs.iter().map(|g| g.symbol).collect()
    }
    #[test]
    fn four_nonpalindromic_end_orientations_have_literal_own_strand_letters_and_true_tips() {
        // AAGC's reverse is CGAA; its reverse complement is GCTT. These must
        // not be interchangeable when displaying an actual bottom protrusion.
        for (side, polarity, top, letters, top_tip, bottom_tip, body) in [
            (Side::Left, Polarity::FivePrime, true, "AAGC", 1., 5., 11.),
            (Side::Left, Polarity::ThreePrime, false, "CGAA", 5., 1., 11.),
            (Side::Right, Polarity::FivePrime, false, "CGAA", 7., 11., 1.),
            (Side::Right, Polarity::ThreePrime, true, "AAGC", 11., 7., 1.),
        ] {
            let original = end(polarity, "AAGC");
            let before = original.clone();
            let actual = layout(&original, side).unwrap();
            assert_eq!(original, before);
            assert_eq!(text(&actual), letters);
            assert_eq!(
                (actual.top_tip, actual.bottom_tip, actual.body_edge),
                (top_tip, bottom_tip, body)
            );
            assert_eq!(actual.width_columns, 12);
            assert_eq!(actual.overhang_length, 4);
            assert!(!actual.abbreviated);
            assert!(actual.glyphs.iter().all(|g| g.top == top));
            let first = if side == Side::Left { 1.5 } else { 7.5 };
            assert_eq!(
                actual.glyphs.iter().map(|g| g.column).collect::<Vec<_>>(),
                vec![first, first + 1., first + 2., first + 3.]
            );
            assert_eq!(
                (actual.top_terminus, actual.bottom_terminus),
                if side == Side::Left {
                    ("5′ P", "3′ OH")
                } else {
                    ("3′ OH", "5′ P")
                }
            );
        }
    }
    #[test]
    fn blunt_unknown_and_blocked_chemistry_remain_labels_without_invented_bases() {
        for side in [Side::Left, Side::Right] {
            let actual = layout(&end(Polarity::Blunt, ""), side).unwrap();
            assert_eq!(actual.width_columns, 8);
            assert!(actual.glyphs.is_empty());
            assert_eq!(actual.overhang_length, 0);
            assert_eq!(
                (actual.top_tip, actual.bottom_tip, actual.body_edge),
                if side == Side::Left {
                    (1., 1., 7.)
                } else {
                    (7., 7., 1.)
                }
            );
            let unknown = End {
                chemistry: Terminus::UNKNOWN,
                ..end(Polarity::Blunt, "")
            };
            let a = layout(&unknown, side).unwrap();
            assert_eq!(
                (a.top_terminus, a.bottom_terminus),
                if side == Side::Left {
                    ("5′ ?", "3′ ?")
                } else {
                    ("3′ ?", "5′ ?")
                }
            );
            let blocked = End {
                chemistry: Terminus {
                    five_prime: FivePrime::Hydroxyl,
                    three_prime: ThreePrime::Blocked,
                },
                ..end(Polarity::FivePrime, "ATG")
            };
            let a = layout(&blocked, side).unwrap();
            assert_eq!(
                (a.top_terminus, a.bottom_terminus),
                if side == Side::Left {
                    ("5′ OH", "3′ blocked")
                } else {
                    ("3′ blocked", "5′ OH")
                }
            );
        }
    }
    #[test]
    fn malformed_shapes_and_non_dna_calls_are_rejected_before_display() {
        assert!(layout(&end(Polarity::Blunt, "A"), Side::Left).is_err());
        for p in [Polarity::FivePrime, Polarity::ThreePrime] {
            assert!(layout(&end(p, ""), Side::Right).is_err());
        }
        for invalid in ["AC U", "ACU", "AC-G", "AC?G", "AΩG", "A\nG"] {
            assert!(layout(&end(Polarity::FivePrime, invalid), Side::Left).is_err());
        }
        let a = layout(&end(Polarity::ThreePrime, "aRykmbvdhn"), Side::Left).unwrap();
        assert_eq!(text(&a), "nhdvbmkyRa");
        let huge = end(
            Polarity::FivePrime,
            &"A".repeat(collagene_core::model::MAX_SEQUENCE_LEN + 1),
        );
        assert!(layout(&huge, Side::Left).is_err());
    }
    #[test]
    fn long_elision_is_bounded_and_applied_after_display_orientation() {
        let letters = "ACGTRYSWKMNNNNNNNNNNNNNNNNNNNNBDHVNACGTA";
        let first = "ACGTRYSWKM";
        let last = "BDHVNACGTA";
        assert!(letters.len() > 24);
        assert!(letters.starts_with(first) && letters.ends_with(last));
        let top = layout(&end(Polarity::FivePrime, letters), Side::Left).unwrap();
        assert_eq!(text(&top), "ACGTRYSWKM…BDHVNACGTA");
        let bottom = layout(&end(Polarity::ThreePrime, letters), Side::Left).unwrap();
        assert_eq!(text(&bottom), "ATGCANVHDB…MKWSYRTGCA");
        for a in [top, bottom] {
            assert_eq!(a.glyphs.len(), 21);
            assert_eq!(a.overhang_length, letters.len());
            assert!(a.abbreviated);
            assert_eq!(a.width_columns, 29);
            assert_eq!(a.body_edge, 28.);
        }
        let exact = layout(
            &end(Polarity::ThreePrime, "ACGTACGTACGTACGTACGTACGT"),
            Side::Right,
        )
        .unwrap();
        assert_eq!(exact.glyphs.len(), 24);
        assert!(!exact.abbreviated);
        assert_eq!(exact.width_columns, 32);
        let maximum = layout(
            &end(
                Polarity::FivePrime,
                &"A".repeat(collagene_core::model::MAX_SEQUENCE_LEN),
            ),
            Side::Right,
        )
        .unwrap();
        assert_eq!(maximum.glyphs.len(), 21);
        assert_eq!(
            maximum.overhang_length,
            collagene_core::model::MAX_SEQUENCE_LEN
        );
        assert_eq!(maximum.width_columns, 29);
    }
}
