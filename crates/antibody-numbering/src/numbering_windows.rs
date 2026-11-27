//! Conditional published Kabat/Chothia window rules. Window boundaries are
//! supplied assumptions, not inferred antibody identity or a CDR definition.
use super::*;
use crate::v_mapping::Locus;
pub const WINDOW_ENGINE: &str = "antibody-numbering-windows-1";
pub const MAX_WINDOW_RESIDUES: usize = 21;
data! {pub enum Convention {Kabat1991,Chothia1997}}
data! {pub enum Window {CysToTrp,MiddleLoop,CysToJ}}
data! {pub struct Layout {pub left:Position,pub residues:Vec<Position>,pub right:Position}}
data! {pub struct Anchors {
    pub source_revision:String,pub left_index:usize,pub right_index:usize,pub provenance:String,
}}
data! {pub struct Report {
    pub engine:String,pub source_revision:String,pub rule_revision:String,pub context_revision:String,
    pub convention:Convention,pub locus:Locus,pub window:Window,pub anchors:Anchors,
    pub start:usize,pub end:usize,pub left_position:Position,pub right_position:Position,
    pub residues:Vec<NumberedResidue>,
}}
pub fn reference(locus: &Locus) -> &'static str {
    match locus {
        Locus::Heavy => {
            "https://www.imgt.org/IMGTScientificChart/Numbering/IMGTnumberingCDR_VH.html"
        }
        Locus::Kappa => {
            "https://www.imgt.org/IMGTScientificChart/Numbering/IMGTnumberingCDR_VK.html"
        }
        Locus::Lambda => {
            "https://www.imgt.org/IMGTScientificChart/Numbering/IMGTnumberingCDR_VL.html"
        }
    }
}
pub fn catalog_revision() -> String {
    hash(b"antibody-numbering-window-rules-1\0", &(
        WINDOW_ENGINE,[reference(&Locus::Heavy),reference(&Locus::Kappa),reference(&Locus::Lambda),"https://www.imgt.org/IMGTScientificChart/Numbering/CDR1-IMGTgaps.html"],
        "2008 correspondence charts and CDR1-IMGTgaps corrected Hil1997; Kabat1991 and Chothia1997 only; uppercase letter labels; CysToTrp has IMGT loop length plus five supplied flank residues; heavy CDR1 Kabat8..12 Chothia8..10; kappa Kabat5..12 Chothia5,6,7,10,11,12; lambda Kabat6..9 Chothia6,7,8,9; MiddleLoop heavy6..10 light3; CysToJ Kabat heavy9..21 kappa7..15 lambda9..15; Chothia heavy11,13,19,21 kappa9 lambda9,11; exact lambda6,7 and kappa11 Chothia patterns; no unsupported extrapolation or CDR definitions",
    ),65536).expect("fixed published window catalog")
}
fn sequence(first: u16, last: u16) -> Vec<Position> {
    (first..=last).map(position).collect()
}
fn insert(number: u16, letters: &str) -> Vec<Position> {
    letters
        .chars()
        .map(|c| Position {
            number,
            insertion: Insertion::Letters(c.to_string()),
        })
        .collect()
}
fn alphabet(n: usize) -> &'static str {
    &"ABCDEFGHIJKLMNOPQRSTUVWXYZ"[..n]
}
fn unsupported<T>() -> Result<T> {
    Err(Error::Invalid(
        "unsupported published numbering window profile",
    ))
}
/// Source-order labels between provided boundaries. CysToTrp includes flanking
/// framework residues, and CysToJ includes positions outside Kabat heavy CDR3.
/// Do not interpret these whole windows as CDR intervals.
pub fn layout(
    convention: &Convention,
    locus: &Locus,
    window: &Window,
    length: usize,
) -> Result<Layout> {
    if length > MAX_WINDOW_RESIDUES {
        return Err(Error::Limit("published numbering window residues"));
    }
    let heavy = *locus == Locus::Heavy;
    let (left, right, residues) = match window {
        Window::CysToTrp => {
            let Some(n) = length.checked_sub(5) else {
                return unsupported();
            };
            let residues = match (convention, locus, n) {
                (Convention::Kabat1991, Locus::Heavy, 8..=12) => {
                    let mut p = sequence(23, 35);
                    p.extend(insert(35, alphabet(n - 8)));
                    p
                }
                (Convention::Chothia1997, Locus::Heavy, 8..=10) => {
                    let mut p = sequence(23, 31);
                    p.extend(insert(31, alphabet(n - 8)));
                    p.extend(sequence(32, 35));
                    p
                }
                (Convention::Kabat1991, Locus::Kappa, 5) => sequence(24, 34)
                    .into_iter()
                    .filter(|p| p.number != 28)
                    .collect(),
                (Convention::Kabat1991, Locus::Kappa, 6..=12)
                | (Convention::Kabat1991, Locus::Lambda, 6..=9) => {
                    let mut p = sequence(24, 27);
                    p.extend(insert(27, alphabet(n - 6)));
                    p.extend(sequence(28, 34));
                    p
                }
                (Convention::Chothia1997, Locus::Kappa, 5) => sequence(24, 34)
                    .into_iter()
                    .filter(|p| p.number != 31)
                    .collect(),
                (Convention::Chothia1997, Locus::Kappa, 6) => sequence(24, 34),
                (Convention::Chothia1997, Locus::Kappa, 7 | 10 | 11 | 12) => {
                    let letters = match n {
                        7 => "A",
                        10 => "ABCD",
                        11 => "ABCEF",
                        12 => "ABCDEF",
                        _ => unreachable!(),
                    };
                    let mut p = sequence(24, 30);
                    p.extend(insert(30, letters));
                    p.extend(sequence(31, 34));
                    p
                }
                (Convention::Chothia1997, Locus::Lambda, 6 | 7) => {
                    let mut p = sequence(24, 27);
                    p.push(position(30));
                    p.extend(insert(30, alphabet(n - 4)));
                    p.extend(sequence(31, 34));
                    p
                }
                (Convention::Chothia1997, Locus::Lambda, 8 | 9) => {
                    let mut p = sequence(24, 30);
                    p.extend(insert(30, alphabet(n - 6)));
                    p.extend(sequence(31, 34));
                    p
                }
                _ => return unsupported(),
            };
            (
                if heavy { 22 } else { 23 },
                if heavy { 36 } else { 35 },
                residues,
            )
        }
        Window::MiddleLoop => {
            let p = if heavy {
                match length {
                    6 => [51, 52, 53, 55, 56, 57].into_iter().map(position).collect(),
                    7..=10 => {
                        let mut p = sequence(51, 52);
                        p.extend(insert(52, alphabet(length - 7)));
                        p.extend(sequence(53, 57));
                        p
                    }
                    _ => return unsupported(),
                }
            } else if length == 3 {
                sequence(50, 52)
            } else {
                return unsupported();
            };
            (if heavy { 50 } else { 49 }, if heavy { 58 } else { 53 }, p)
        }
        Window::CysToJ => {
            let supported = match (convention, locus) {
                (Convention::Kabat1991, Locus::Heavy) => (9..=21).contains(&length),
                (Convention::Kabat1991, Locus::Kappa) => (7..=15).contains(&length),
                (Convention::Kabat1991, Locus::Lambda) => (9..=15).contains(&length),
                (Convention::Chothia1997, Locus::Heavy) => matches!(length, 11 | 13 | 19 | 21),
                (Convention::Chothia1997, Locus::Kappa) => length == 9,
                (Convention::Chothia1997, Locus::Lambda) => matches!(length, 9 | 11),
            };
            if !supported {
                return unsupported();
            }
            let p = if heavy {
                if length == 9 {
                    let mut p = sequence(93, 97);
                    p.extend(sequence(99, 102));
                    p
                } else {
                    let mut p = sequence(93, 100);
                    p.extend(insert(100, alphabet(length - 10)));
                    p.extend(sequence(101, 102));
                    p
                }
            } else if length < 9 {
                let mut p = sequence(89, 92);
                p.extend(sequence(if length == 7 { 95 } else { 94 }, 97));
                p
            } else {
                let mut p = sequence(89, 95);
                p.extend(insert(95, alphabet(length - 9)));
                p.extend(sequence(96, 97));
                p
            };
            (if heavy { 92 } else { 88 }, if heavy { 103 } else { 98 }, p)
        }
    };
    if residues.len() != length {
        return Err(Error::Binding);
    }
    Ok(Layout {
        left: position(left),
        residues,
        right: position(right),
    })
}
/// Conditional under provided counterpart IMGT boundaries 23/41, 55/66 or
/// 104/118 and the catalog's flank assumptions. Literal anchor identities cannot
/// prove the numbering, locus, domain classification or a CDR definition.
pub fn number_anchored_window(
    source: &Molecule,
    convention: &Convention,
    locus: &Locus,
    window: &Window,
    anchors: &Anchors,
) -> Result<Report> {
    let revision = source_revision(source)?;
    if revision != anchors.source_revision {
        return Err(Error::Binding);
    }
    if anchors.provenance.trim().is_empty()
        || anchors.provenance.len() > 1024
        || anchors.provenance.chars().any(char::is_control)
    {
        return Err(Error::Invalid("numbering window provenance"));
    }
    let start = anchors
        .left_index
        .checked_add(1)
        .ok_or(Error::Invalid("numbering window boundary"))?;
    let end = anchors.right_index;
    let bytes = source.sequence.as_str().as_bytes();
    if start > end || end >= bytes.len() {
        return Err(Error::Invalid("numbering window boundary order/bounds"));
    }
    if bytes[anchors.left_index..=end].contains(&b'*') {
        return Err(Error::Invalid("stop in numbering window/boundary"));
    }
    match window {
        Window::CysToTrp if bytes[anchors.left_index] != b'C' || bytes[end] != b'W' => {
            return Err(Error::Invalid("provided Cys/Trp window anchors"));
        }
        Window::CysToJ => {
            let motif = bytes
                .get(
                    end..end
                        .checked_add(4)
                        .ok_or(Error::Invalid("J window boundary"))?,
                )
                .ok_or(Error::Invalid("incomplete J window context"))?;
            let expected = if *locus == Locus::Heavy { b'W' } else { b'F' };
            if bytes[anchors.left_index] != b'C'
                || motif[0] != expected
                || motif[1] != b'G'
                || motif[2] == b'*'
                || motif[3] != b'G'
            {
                return Err(Error::Invalid("provided Cys/J window anchors"));
            }
        }
        _ => (),
    }
    let labels = layout(convention, locus, window, end - start)?;
    let rule_revision = catalog_revision();
    let context_revision = hash(
        b"antibody-numbering-window-context-1\0",
        &(
            WINDOW_ENGINE,
            &revision,
            &rule_revision,
            convention,
            locus,
            window,
            anchors,
        ),
        65536,
    )?;
    Ok(Report {
        engine: WINDOW_ENGINE.into(),
        source_revision: revision,
        rule_revision,
        context_revision,
        convention: convention.clone(),
        locus: locus.clone(),
        window: window.clone(),
        anchors: anchors.clone(),
        start,
        end,
        left_position: labels.left,
        right_position: labels.right,
        residues: labels
            .residues
            .into_iter()
            .enumerate()
            .map(|(i, position)| NumberedResidue {
                source_index: start + i,
                residue: bytes[start + i] as char,
                position,
            })
            .collect(),
    })
}
impl Report {
    pub fn validate_for(
        &self,
        source: &Molecule,
        convention: &Convention,
        locus: &Locus,
        window: &Window,
        anchors: &Anchors,
    ) -> Result<()> {
        if self.engine != WINDOW_ENGINE
            || self.residues.len() > MAX_WINDOW_RESIDUES
            || *self != number_anchored_window(source, convention, locus, window, anchors)?
        {
            return Err(Error::Binding);
        }
        Ok(())
    }
}
