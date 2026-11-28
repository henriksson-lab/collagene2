//! Exact two-strand restriction and ordered-ligation geometry.
//! Strands are stored 5′→3′; bottom_start is the left boundary of the antiparallel
//! bottom strand in top-reference coordinates. No kinetics or methylation inference.
use crate::analysis::{self, Enzyme};
use collagene_model::{Alphabet, Error, Molecule, Sequence, Topology};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum FivePrime {
    Phosphate,
    Hydroxyl,
    Unknown,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum ThreePrime {
    Hydroxyl,
    Blocked,
    Unknown,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum EndTreatment {
    Phosphorylate,
    Dephosphorylate,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Terminus {
    pub five_prime: FivePrime,
    pub three_prime: ThreePrime,
}
impl Terminus {
    pub const UNKNOWN: Self = Self {
        five_prime: FivePrime::Unknown,
        three_prime: ThreePrime::Unknown,
    };
    pub const CLEAVED: Self = Self {
        five_prime: FivePrime::Phosphate,
        three_prime: ThreePrime::Hydroxyl,
    };
    pub const PCR: Self = Self {
        five_prime: FivePrime::Hydroxyl,
        three_prime: ThreePrime::Hydroxyl,
    };
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum Polarity {
    Blunt,
    FivePrime,
    ThreePrime,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub struct End {
    pub polarity: Polarity,
    /// The protruding strand, 5′→3′; empty for blunt ends.
    pub overhang: String,
    pub chemistry: Terminus,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum Strand {
    Top,
    Bottom,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub struct Nick {
    pub strand: Strand,
    /// Boundary in top-reference coordinates; circular boundaries are modulo length.
    pub position: i64,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub struct Duplex {
    pub top: Sequence,
    pub bottom: Sequence,
    pub bottom_start: i64,
    pub topology: Topology,
    pub left: Terminus,
    pub right: Terminus,
    pub nicks: Vec<Nick>,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub struct DuplexSummary {
    pub model: String,
    pub top_length: usize,
    pub bottom_length: usize,
    pub bottom_start: i64,
    pub topology: Topology,
    pub ends: Option<(End, End)>,
    pub nicks: Vec<Nick>,
    pub covalently_closed_in_model: bool,
}
impl DuplexSummary {
    pub fn describe(&self) -> String {
        let end = |e: &End| {
            format!(
                "{}; 5′ terminus {:?}, 3′ terminus {:?}",
                match e.polarity {
                    Polarity::Blunt => "blunt".into(),
                    Polarity::FivePrime => format!("5′ overhang {} (5′→3′)", e.overhang),
                    Polarity::ThreePrime => format!("3′ overhang {} (5′→3′)", e.overhang),
                },
                e.chemistry.five_prime,
                e.chemistry.three_prime
            )
        };
        let ends = match &self.ends {
            Some((left, right)) => format!("Left end: {}\nRight end: {}", end(left), end(right)),
            None => "No free ends.".into(),
        };
        let nicks = if self.nicks.is_empty() {
            "none".into()
        } else {
            self.nicks
                .iter()
                .map(|n| format!("{:?} strand at {}", n.strand, n.position))
                .collect::<Vec<_>>()
                .join(", ")
        };
        format!(
            "Top strand: {} nt; bottom strand: {} nt.\n{}\nRetained nicks (zero-based top-reference boundaries): {}\nCovalently closed in this model: {}",
            self.top_length, self.bottom_length, ends, nicks, self.covalently_closed_in_model
        )
    }
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub struct DigestFragment {
    pub duplex: Duplex,
    /// Top-strand reference boundary in the source (zero-based).
    pub top_start: usize,
    /// Unwrapped bottom-strand reference boundary in the source.
    pub bottom_start: i64,
}

fn unsupported(message: &str) -> Error {
    Error::Unsupported(message.into())
}
fn invalid(message: &str) -> Error {
    Error::Format(message.into())
}
fn exact_dna(sequence: &Sequence) -> Result<(), Error> {
    if sequence.alphabet() != Alphabet::Dna
        || sequence.as_str().bytes().any(|b| !b"ACGT".contains(&b))
    {
        return Err(unsupported(
            "duplex chemistry currently requires unambiguous DNA",
        ));
    }
    Ok(())
}
/// Source metadata may restrict a model, but cannot establish phosphorylation or erase nicks.
pub fn has_unmapped_physical_evidence(molecule: &Molecule) -> bool {
    let non_duplex = |kind: &str| {
        let kind = kind.to_ascii_uppercase();
        kind.starts_with("SS-") || kind.starts_with("MS-")
    };
    molecule.metadata.contains_key("duplex_model")
        || molecule
            .metadata
            .contains_key(crate::material_pcr::ORIGIN_METADATA)
        || molecule.metadata.keys().any(|key| {
            key.starts_with("source:material_derivation:")
                || key.starts_with("source:ligation_outcomes:")
        })
        || molecule
            .metadata
            .get("source:legacy:chemistry_unmapped")
            .is_some_and(|value| value == "true")
        || molecule
            .metadata
            .get("source:snapgene:chemistry_unmapped")
            .is_some_and(|s| s == "true")
        || molecule
            .metadata
            .get("genbank.molecule_type")
            .is_some_and(|s| non_duplex(s))
        || molecule
            .original_genbank
            .as_ref()
            .and_then(|s| s.lines().find(|line| line.starts_with("LOCUS")))
            .is_some_and(|line| line.split_whitespace().any(non_duplex))
}
impl Duplex {
    /// A starting sequence is modeled as intact duplex DNA, with unknown free-end chemistry.
    /// Imported metadata cannot assert phosphorylation or restore a physical intermediate.
    pub fn from_molecule(molecule: &Molecule) -> Result<Self, Error> {
        molecule.validate()?;
        if has_unmapped_physical_evidence(molecule) {
            return Err(unsupported(
                "single/mixed-strand or imported physical-intermediate evidence is not an established intact duplex; reconstruct a supported originating graph",
            ));
        }
        exact_dna(&molecule.sequence)?;
        Ok(Self {
            top: molecule.sequence.clone(),
            bottom: molecule.sequence.reverse_complement()?,
            bottom_start: 0,
            topology: molecule.topology,
            left: Terminus::UNKNOWN,
            right: Terminus::UNKNOWN,
            nicks: vec![],
        })
    }
    /// Ideal blunt product of the exact-PCR model, with unmodified oligos.
    pub fn exact_pcr(sequence: &Sequence) -> Result<Self, Error> {
        exact_dna(sequence)?;
        Ok(Self {
            top: sequence.clone(),
            bottom: sequence.reverse_complement()?,
            bottom_start: 0,
            topology: Topology::Linear,
            left: Terminus::PCR,
            right: Terminus::PCR,
            nicks: vec![],
        })
    }
    fn bottom_end(&self) -> i64 {
        self.bottom_start + self.bottom.len() as i64
    }
    pub fn validate(&self) -> Result<(), Error> {
        exact_dna(&self.top)?;
        exact_dna(&self.bottom)?;
        let n = self.top.len() as i64;
        if self.bottom_start.unsigned_abs() > collagene_model::MAX_SEQUENCE_LEN as u64 {
            return Err(Error::TooLarge);
        }
        let start = self.bottom_start.max(0);
        let end = self.bottom_end().min(n);
        if start >= end {
            return Err(unsupported(
                "single-stranded or non-overlapping strand fragments require a more detailed model",
            ));
        }
        if self.topology == Topology::Circular
            && (self.bottom_start != 0 || self.bottom.len() != self.top.len())
        {
            return Err(invalid(
                "circular duplex must have aligned, equal-length strands",
            ));
        }
        let bottom = self.bottom.as_str().as_bytes();
        for i in start..end {
            let b = bottom[self.bottom.len() - 1 - (i - self.bottom_start) as usize];
            if self.top.as_str().as_bytes()[i as usize]
                != collagene_model::complement(b, Alphabet::Dna)
            {
                return Err(invalid("duplex strands disagree within the paired region"));
            }
        }
        if self.nicks.len() > 10_000 {
            return Err(Error::TooLarge);
        }
        let mut seen_nicks = BTreeSet::new();
        for nick in &self.nicks {
            if !seen_nicks.insert((nick.strand == Strand::Bottom, nick.position)) {
                return Err(invalid("duplicate nick boundary on the same strand"));
            }
            let (low, high) = match nick.strand {
                Strand::Top => (0, n),
                Strand::Bottom => (self.bottom_start, self.bottom_end()),
            };
            let valid = if self.topology == Topology::Circular {
                nick.position >= 0 && nick.position < n
            } else {
                nick.position > low && nick.position < high
            };
            if !valid {
                return Err(invalid("nick boundary lies outside its strand"));
            }
        }
        Ok(())
    }
    pub fn ends(&self) -> Result<(End, End), Error> {
        self.validate()?;
        if self.topology == Topology::Circular {
            return Err(unsupported("circular duplex has no free ends"));
        }
        let s = self.bottom_start;
        let e = self.bottom_end();
        let n = self.top.len() as i64;
        let (polarity, overhang) = if s > 0 {
            (Polarity::FivePrime, self.top.as_str()[..s as usize].into())
        } else if s < 0 {
            (
                Polarity::ThreePrime,
                self.bottom.as_str()[self.bottom.len() - (-s) as usize..].into(),
            )
        } else {
            (Polarity::Blunt, String::new())
        };
        let left = End {
            polarity,
            overhang,
            chemistry: self.left,
        };
        let (polarity, overhang) = if e > n {
            (
                Polarity::FivePrime,
                self.bottom.as_str()[..(e - n) as usize].into(),
            )
        } else if e < n {
            (Polarity::ThreePrime, self.top.as_str()[e as usize..].into())
        } else {
            (Polarity::Blunt, String::new())
        };
        Ok((
            left,
            End {
                polarity,
                overhang,
                chemistry: self.right,
            },
        ))
    }
    pub fn reverse(&self) -> Result<Self, Error> {
        self.validate()?;
        let end = self.bottom_end();
        let n = self.top.len() as i64;
        let mut reversed = Self {
            top: self.bottom.clone(),
            bottom: self.top.clone(),
            bottom_start: end - n,
            topology: self.topology,
            left: self.right,
            right: self.left,
            nicks: self
                .nicks
                .iter()
                .map(|nick| Nick {
                    strand: match nick.strand {
                        Strand::Top => Strand::Bottom,
                        Strand::Bottom => Strand::Top,
                    },
                    position: if self.topology == Topology::Circular {
                        (end - nick.position).rem_euclid(n)
                    } else {
                        end - nick.position
                    },
                })
                .collect(),
        };
        if self.topology == Topology::Circular {
            reversed.bottom_start = 0;
        }
        reversed.validate()?;
        Ok(reversed)
    }
    pub fn treat(&self, treatment: EndTreatment) -> Result<Self, Error> {
        self.validate()?;
        if self.topology != Topology::Linear || !self.nicks.is_empty() {
            return Err(unsupported(
                "terminal end treatment requires a linear duplex without internal nicks",
            ));
        }
        let mut result = self.clone();
        let state = match treatment {
            EndTreatment::Phosphorylate => FivePrime::Phosphate,
            EndTreatment::Dephosphorylate => FivePrime::Hydroxyl,
        };
        result.left.five_prime = state;
        result.right.five_prime = state;
        Ok(result)
    }
    pub fn molecule(&self, name: impl Into<String>) -> Result<Molecule, Error> {
        self.validate()?;
        let mut m = Molecule::dna(name, self.top.as_str(), self.topology)?;
        // Inspectable evidence only. Reading a molecule never trusts this metadata as chemistry.
        m.metadata.insert(
            "duplex_model".into(),
            serde_json::to_string(&self.summary()?).map_err(|e| Error::Format(e.to_string()))?,
        );
        Ok(m)
    }
    pub fn summary(&self) -> Result<DuplexSummary, Error> {
        self.validate()?;
        let ends = if self.topology == Topology::Linear {
            self.ends().ok()
        } else {
            None
        };
        Ok(DuplexSummary {
            model: "exact duplex geometry".into(),
            top_length: self.top.len(),
            bottom_length: self.bottom.len(),
            bottom_start: self.bottom_start,
            topology: self.topology,
            ends,
            nicks: self.nicks.clone(),
            covalently_closed_in_model: self.topology == Topology::Circular
                && self.nicks.is_empty(),
        })
    }
}

#[derive(Clone, Copy)]
struct Boundary {
    top: i64,
    bottom: i64,
    chemistry: Terminus,
}
fn extract(
    sequence: &Sequence,
    start: i64,
    length: usize,
    circular: bool,
    reference_start: i64,
) -> Result<Sequence, Error> {
    if length == 0 || length > collagene_model::MAX_SEQUENCE_LEN {
        return Err(Error::TooLarge);
    }
    if circular {
        if length > sequence.len() {
            return Err(unsupported(
                "strand interval spans more than one molecule turn",
            ));
        }
        let index = (start - reference_start).rem_euclid(sequence.len() as i64) as usize;
        let first = length.min(sequence.len() - index);
        let mut bases = sequence.as_str()[index..index + first].to_owned();
        bases.push_str(&sequence.as_str()[..length - first]);
        Sequence::dna(&bases)
    } else {
        let from = start
            .checked_sub(reference_start)
            .and_then(|n| usize::try_from(n).ok())
            .ok_or_else(|| invalid("strand interval outside source"))?;
        let bases = sequence
            .as_str()
            .get(from..from.saturating_add(length))
            .ok_or_else(|| invalid("strand interval outside source"))?;
        Sequence::dna(bases)
    }
}

/// Complete digestion, sorted by ascending top-reference cut. The last circular fragment wraps.
/// Cut-free circular input remains circular; it is never presented as a linear fragment.
pub fn digest(duplex: &Duplex, enzymes: &[Enzyme]) -> Result<Vec<DigestFragment>, Error> {
    duplex.validate()?;
    if enzymes.is_empty() || enzymes.len() > 64 {
        return Err(invalid("digest needs 1–64 enzymes"));
    }
    if !duplex.nicks.is_empty() {
        return Err(unsupported(
            "digest propagation of existing nicks is not yet modeled",
        ));
    }
    let molecule = Molecule::dna("digest source", duplex.top.as_str(), duplex.topology)?;
    let mut cuts = BTreeMap::new();
    for enzyme in enzymes {
        let recognition = Sequence::dna(&enzyme.recognition)?;
        let delta = i64::from(enzyme.bottom_cut) - i64::from(enzyme.top_cut);
        for site in analysis::restriction_sites(&molecule, enzyme)? {
            add_digest_cut(duplex, &site, &recognition.len(), &delta, &mut cuts)?;
        }
    }
    materialize_digest(duplex, cuts)
}

/// Shared geometry for already assessed actual cuts; does not assess methylation.
pub(crate) fn digest_sites(
    duplex: &Duplex,
    sites: &[(analysis::CutSite, usize, i64)],
) -> Result<Vec<DigestFragment>, Error> {
    duplex.validate()?;
    if !duplex.nicks.is_empty() {
        return Err(unsupported(
            "digest propagation of existing nicks is not yet modeled",
        ));
    }
    let mut cuts = BTreeMap::new();
    for (site, recognition_length, delta) in sites {
        add_digest_cut(duplex, site, recognition_length, delta, &mut cuts)?;
    }
    materialize_digest(duplex, cuts)
}
fn add_digest_cut(
    duplex: &Duplex,
    site: &analysis::CutSite,
    recognition_length: &usize,
    delta: &i64,
    cuts: &mut BTreeMap<i64, i64>,
) -> Result<(), Error> {
    let n = duplex.top.len() as i64;
    let circular = duplex.topology == Topology::Circular;
    if delta.unsigned_abs() >= n as u64 {
        return Err(unsupported(
            "cleavage stagger is at least one source-molecule length",
        ));
    }
    if !site.certain {
        return Err(unsupported(
            "ambiguous restriction sites cannot establish a complete duplex digest",
        ));
    }
    let top = site.top_cut as i64;
    let bottom = top + delta;
    if !circular {
        let recognition_start = site.recognition_start as i64;
        if recognition_start < duplex.bottom_start.max(0)
            || recognition_start + *recognition_length as i64 > duplex.bottom_end().min(n)
        {
            return Err(unsupported(
                "restriction recognition site is not fully double-stranded",
            ));
        }
        if top <= 0 || top >= n || bottom <= duplex.bottom_start || bottom >= duplex.bottom_end() {
            return Err(unsupported(
                "cleavage coincides with a free end or lacks a complete strand substrate",
            ));
        }
    }
    if let Some(previous) = cuts.insert(top, bottom) {
        if previous != bottom {
            return Err(unsupported(
                "enzymes cut the same top boundary with conflicting bottom boundaries",
            ));
        }
    }
    if cuts.len() > 10_000 {
        return Err(unsupported(
            "duplex digest exceeds the 10,000-fragment budget",
        ));
    }
    Ok(())
}
fn materialize_digest(
    duplex: &Duplex,
    cuts: BTreeMap<i64, i64>,
) -> Result<Vec<DigestFragment>, Error> {
    let n = duplex.top.len() as i64;
    let circular = duplex.topology == Topology::Circular;
    if cuts.is_empty() {
        return Ok(vec![DigestFragment {
            duplex: duplex.clone(),
            top_start: 0,
            bottom_start: duplex.bottom_start,
        }]);
    }
    let mut boundaries: Vec<_> = cuts
        .into_iter()
        .map(|(top, bottom)| Boundary {
            top,
            bottom,
            chemistry: Terminus::CLEAVED,
        })
        .collect();
    if circular {
        let first = boundaries[0];
        boundaries.push(Boundary {
            top: first.top + n,
            bottom: first.bottom + n,
            ..first
        });
    } else {
        boundaries.insert(
            0,
            Boundary {
                top: 0,
                bottom: duplex.bottom_start,
                chemistry: duplex.left,
            },
        );
        boundaries.push(Boundary {
            top: n,
            bottom: duplex.bottom_end(),
            chemistry: duplex.right,
        });
    }
    let bottom_forward = duplex.bottom.reverse_complement()?;
    let mut fragments = Vec::new();
    for pair in boundaries.windows(2) {
        let (a, b) = (pair[0], pair[1]);
        if b.top <= a.top || b.bottom <= a.bottom || a.top.max(a.bottom) >= b.top.min(b.bottom) {
            return Err(unsupported(
                "crossing or overlapping strand cuts produce a single-stranded intermediate outside this model",
            ));
        }
        let top = extract(&duplex.top, a.top, (b.top - a.top) as usize, circular, 0)?;
        let bottom = extract(
            &bottom_forward,
            a.bottom,
            (b.bottom - a.bottom) as usize,
            circular,
            duplex.bottom_start,
        )?
        .reverse_complement()?;
        let fragment = Duplex {
            top,
            bottom,
            bottom_start: a.bottom - a.top,
            topology: Topology::Linear,
            left: a.chemistry,
            right: b.chemistry,
            nicks: vec![],
        };
        fragment.validate()?;
        fragments.push(DigestFragment {
            duplex: fragment,
            top_start: a.top.rem_euclid(n) as usize,
            bottom_start: a.bottom,
        });
    }
    Ok(fragments)
}

fn compatible(right: &End, left: &End) -> Result<(), Error> {
    if right.polarity != left.polarity || right.overhang.len() != left.overhang.len() {
        return Err(invalid(
            "ligation ends have incompatible overhang polarity or length",
        ));
    }
    if right.polarity != Polarity::Blunt
        && Sequence::dna(&right.overhang)?
            .reverse_complement()?
            .as_str()
            != left.overhang
    {
        return Err(invalid("ligation overhang sequences are not complementary"));
    }
    Ok(())
}
fn seal(five: FivePrime, three: ThreePrime) -> Result<bool, Error> {
    if five == FivePrime::Unknown || three == ThreePrime::Unknown {
        return Err(unsupported(
            "junction phosphorylation or 3′ terminal chemistry is unknown",
        ));
    }
    Ok(five == FivePrime::Phosphate && three == ThreePrime::Hydroxyl)
}
fn junction(
    right: &End,
    left: &End,
    top_position: i64,
    bottom_position: i64,
    allow_nicked: bool,
) -> Result<Vec<Nick>, Error> {
    compatible(right, left)?;
    let top = seal(left.chemistry.five_prime, right.chemistry.three_prime)?;
    let bottom = seal(right.chemistry.five_prime, left.chemistry.three_prime)?;
    if !top && !bottom {
        return Err(invalid(
            "neither strand can form a phosphodiester bond at the ligation junction",
        ));
    }
    if (!top || !bottom) && !allow_nicked {
        return Err(invalid(
            "ligation would retain nicks but this operation requires sealed junctions",
        ));
    }
    let mut nicks = vec![];
    if !top {
        nicks.push(Nick {
            strand: Strand::Top,
            position: top_position,
        });
    }
    if !bottom {
        nicks.push(Nick {
            strand: Strand::Bottom,
            position: bottom_position,
        });
    }
    Ok(nicks)
}

/// Check one explicitly ordered, oriented ligation candidate, not all mixture products.
/// A junction must form at least one bond. Remaining nicks require explicit permission.
pub fn ligate(fragments: &[Duplex], circular: bool, allow_nicked: bool) -> Result<Duplex, Error> {
    if fragments.is_empty() || fragments.len() > 1000 {
        return Err(invalid("ligation needs 1–1000 fragments"));
    }
    let nick_count = fragments
        .iter()
        .try_fold(0usize, |n, f| n.checked_add(f.nicks.len()))
        .ok_or(Error::TooLarge)?;
    if nick_count > 10_000 {
        return Err(Error::TooLarge);
    }
    let total_top = fragments
        .iter()
        .try_fold(0usize, |n, f| n.checked_add(f.top.len()))
        .ok_or(Error::TooLarge)?;
    let total_bottom = fragments
        .iter()
        .try_fold(0usize, |n, f| n.checked_add(f.bottom.len()))
        .ok_or(Error::TooLarge)?;
    if total_top > collagene_model::MAX_SEQUENCE_LEN
        || total_bottom > collagene_model::MAX_SEQUENCE_LEN
    {
        return Err(Error::TooLarge);
    }
    let mut top = String::new();
    let mut bottom_parts = Vec::new();
    let mut nicks = vec![];
    for (i, fragment) in fragments.iter().enumerate() {
        fragment.validate()?;
        if fragment.topology != Topology::Linear {
            return Err(invalid(
                "ligation substrates must be linear duplex fragments",
            ));
        }
        if top.len().saturating_add(fragment.top.len()) > collagene_model::MAX_SEQUENCE_LEN {
            return Err(Error::TooLarge);
        }
        let offset = top.len() as i64;
        if i > 0 {
            let previous = &fragments[i - 1];
            let right = previous.ends()?.1;
            let left = fragment.ends()?.0;
            nicks.extend(junction(
                &right,
                &left,
                offset,
                offset + fragment.bottom_start,
                allow_nicked,
            )?);
        }
        nicks.extend(fragment.nicks.iter().map(|nick| Nick {
            strand: nick.strand,
            position: nick.position + offset,
        }));
        if nicks.len() > 10_000 {
            return Err(Error::TooLarge);
        }
        top.push_str(fragment.top.as_str());
        bottom_parts.push(fragment.bottom.as_str());
    }
    let first = &fragments[0];
    let last = fragments.last().unwrap();
    let mut result = Duplex {
        top: Sequence::dna(&top)?,
        bottom: Sequence::dna(&bottom_parts.into_iter().rev().collect::<String>())?,
        bottom_start: first.bottom_start,
        topology: Topology::Linear,
        left: first.left,
        right: last.right,
        nicks,
    };
    // End compatibility guarantees both strand seams align; validate also checks all paired bases.
    result.validate()?;
    if circular {
        let (left, right) = result.ends()?;
        let n = result.top.len() as i64;
        if result.bottom.len() != result.top.len() {
            return Err(invalid(
                "circular ligation does not close both strand lengths",
            ));
        }
        result.nicks.extend(junction(
            &right,
            &left,
            n,
            result.bottom_end(),
            allow_nicked,
        )?);
        let forward = result.bottom.reverse_complement()?;
        let aligned = extract(&forward, 0, result.top.len(), true, result.bottom_start)?;
        if aligned != result.top {
            return Err(invalid(
                "circular ligation leaves discordant strand sequences",
            ));
        }
        result.bottom = aligned.reverse_complement()?;
        result.bottom_start = 0;
        result.topology = Topology::Circular;
        for nick in &mut result.nicks {
            nick.position = nick.position.rem_euclid(n);
        }
        result.left = Terminus::UNKNOWN;
        result.right = Terminus::UNKNOWN;
    }
    if !allow_nicked && !result.nicks.is_empty() {
        return Err(invalid("ligation product retains pre-existing nicks"));
    }
    result.validate()?;
    Ok(result)
}

/// Ideal end geometry, not enzyme kinetics or a guarantee of laboratory completion.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum EndRepairMethod {
    FillFivePrimeOverhangs,
    RemoveThreePrimeOverhangs,
    Blunt,
}

/// Source/product indices always address the actual strand stored 5′→3′.
/// Cross-strand synthesized mappings reverse-complement their source bases.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub struct StrandMapping {
    pub source_strand: Strand,
    pub source_start: usize,
    pub source_end: usize,
    pub product_strand: Strand,
    pub product_start: usize,
    pub product_end: usize,
    pub reverse: bool,
    pub synthesized: bool,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub struct EndRepairReport {
    pub method: EndRepairMethod,
    pub source: Duplex,
    pub duplex: Duplex,
    pub mappings: Vec<StrandMapping>,
    pub top_added_3prime: String,
    pub top_removed_3prime: String,
    pub bottom_added_3prime: String,
    pub bottom_removed_3prime: String,
    pub warnings: Vec<String>,
}

pub const MAX_CHEMISTRY_REPORT_BYTES: usize = 32 * 1024 * 1024;
pub const MAX_ANNEAL_OLIGO_BASES: usize = 10_000;
pub const MAX_ANNEAL_COMPARISONS: usize = 200_000_000;
pub const MAX_ANNEAL_CANDIDATES: usize = 1000;
pub const ANNEAL_ENGINE_VERSION: &str = concat!(env!("CARGO_PKG_VERSION"), "/annealing-1");

fn bounded_chemistry_json(value: &impl Serialize) -> Result<(), Error> {
    struct Count(usize);
    impl std::io::Write for Count {
        fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
            self.0 = self
                .0
                .checked_add(bytes.len())
                .filter(|n| *n <= MAX_CHEMISTRY_REPORT_BYTES)
                .ok_or_else(|| std::io::Error::other("chemistry report limit"))?;
            Ok(bytes.len())
        }
        fn flush(&mut self) -> std::io::Result<()> {
            Ok(())
        }
    }
    serde_json::to_writer(&mut Count(0), value).map_err(|_| Error::TooLarge)
}
fn identity_mapping(strand: Strand, length: usize) -> StrandMapping {
    StrandMapping {
        source_strand: strand,
        source_start: 0,
        source_end: length,
        product_strand: strand,
        product_start: 0,
        product_end: length,
        reverse: false,
        synthesized: false,
    }
}
fn reverse_complement_text(bases: &str) -> Result<String, Error> {
    Ok(Sequence::dna(bases)?.reverse_complement()?.as_str().into())
}

/// Fill recessed 3′ ends from their existing opposite-strand template and/or
/// remove protruding 3′ ends to the initial paired boundary. This never modifies
/// 5′ chemistry, seals a nick, or repairs an unknown/blocked 3′ terminus.
pub fn repair_ends(source: &Duplex, method: EndRepairMethod) -> Result<EndRepairReport, Error> {
    source.validate()?;
    if source.topology != Topology::Linear {
        return Err(unsupported(
            "end repair requires free ends on a linear duplex",
        ));
    }
    let n = source.top.len() as i64;
    let s = source.bottom_start;
    let e = source.bottom_end();
    let paired_start = s.max(0);
    let paired_end = e.min(n);
    if source
        .nicks
        .iter()
        .any(|nick| nick.position <= paired_start || nick.position >= paired_end)
    {
        return Err(unsupported(
            "end repair cannot establish an intact overhang/template when a nick is at or outside its paired boundary",
        ));
    }
    let fill = matches!(
        method,
        EndRepairMethod::FillFivePrimeOverhangs | EndRepairMethod::Blunt
    );
    let trim = matches!(
        method,
        EndRepairMethod::RemoveThreePrimeOverhangs | EndRepairMethod::Blunt
    );
    let change_left = (s > 0 && fill) || (s < 0 && trim);
    let change_right = (e > n && fill) || (e < n && trim);
    if (change_left && source.left.three_prime != ThreePrime::Hydroxyl)
        || (change_right && source.right.three_prime != ThreePrime::Hydroxyl)
    {
        return Err(unsupported(
            "affected end repair requires a known accessible 3′ hydroxyl; blocked/unknown ends are not repaired by assumption",
        ));
    }
    let new_top_length = if change_right {
        e as usize
    } else {
        source.top.len()
    };
    let new_bottom_start = if change_left { 0 } else { s };
    let new_bottom_length = (e - new_bottom_start) as usize;
    if new_top_length > collagene_model::MAX_SEQUENCE_LEN
        || new_bottom_length > collagene_model::MAX_SEQUENCE_LEN
    {
        return Err(Error::TooLarge);
    }
    let total = source
        .top
        .len()
        .checked_add(source.bottom.len())
        .and_then(|n| n.checked_add(new_top_length))
        .and_then(|n| n.checked_add(new_bottom_length))
        .ok_or(Error::TooLarge)?;
    if total > MAX_CHEMISTRY_REPORT_BYTES {
        return Err(Error::TooLarge);
    }
    let top_retained = source.top.len().min(new_top_length);
    let bottom_retained = source.bottom.len().min(new_bottom_length);
    let top_added = if new_top_length > source.top.len() {
        reverse_complement_text(&source.bottom.as_str()[..new_top_length - source.top.len()])?
    } else {
        String::new()
    };
    let bottom_added = if new_bottom_length > source.bottom.len() {
        reverse_complement_text(&source.top.as_str()[..new_bottom_length - source.bottom.len()])?
    } else {
        String::new()
    };
    let top_removed = source.top.as_str()[top_retained..].to_string();
    let bottom_removed = source.bottom.as_str()[bottom_retained..].to_string();
    let mut top = source.top.as_str()[..top_retained].to_string();
    top.push_str(&top_added);
    let mut bottom = source.bottom.as_str()[..bottom_retained].to_string();
    bottom.push_str(&bottom_added);
    let mut result = source.clone();
    result.top = Sequence::dna(&top)?;
    result.bottom = Sequence::dna(&bottom)?;
    result.bottom_start = new_bottom_start;
    if change_left {
        result.left.three_prime = ThreePrime::Hydroxyl;
    }
    if change_right {
        result.right.three_prime = ThreePrime::Hydroxyl;
    }
    result.validate()?;
    let mut mappings = vec![
        identity_mapping(Strand::Top, top_retained),
        identity_mapping(Strand::Bottom, bottom_retained),
    ];
    if !top_added.is_empty() {
        mappings.push(StrandMapping {
            source_strand: Strand::Bottom,
            source_start: 0,
            source_end: top_added.len(),
            product_strand: Strand::Top,
            product_start: source.top.len(),
            product_end: new_top_length,
            reverse: true,
            synthesized: true,
        });
    }
    if !bottom_added.is_empty() {
        mappings.push(StrandMapping {
            source_strand: Strand::Top,
            source_start: 0,
            source_end: bottom_added.len(),
            product_strand: Strand::Bottom,
            product_start: source.bottom.len(),
            product_end: new_bottom_length,
            reverse: true,
            synthesized: true,
        });
    }
    let mut warnings=vec!["Ideal specified end geometry only. Polymerase/exonuclease kinetics, excessive resection, reagent availability and laboratory completion are not predicted. 5′ chemistry and internal nicks are preserved; no nick repair, phosphorylation or blocking-group removal is inferred.".into()];
    if !change_left && !change_right {
        warnings.push(
            "No end has the selected repair polarity; sequence and chemistry are unchanged.".into(),
        );
    }
    if !source.nicks.is_empty() {
        warnings.push("Existing internal paired-region nicks remain unsealed.".into());
    }
    let report = EndRepairReport {
        method,
        source: source.clone(),
        duplex: result,
        mappings,
        top_added_3prime: top_added,
        top_removed_3prime: top_removed,
        bottom_added_3prime: bottom_added,
        bottom_removed_3prime: bottom_removed,
        warnings,
    };
    bounded_chemistry_json(&report)?;
    Ok(report)
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct AnnealOptions {
    pub min_paired_bases: usize,
    pub max_comparisons: usize,
    pub max_candidates: usize,
    pub top_chemistry: Terminus,
    pub bottom_chemistry: Terminus,
}
impl Default for AnnealOptions {
    fn default() -> Self {
        Self {
            min_paired_bases: 12,
            max_comparisons: 20_000_000,
            max_candidates: 100,
            top_chemistry: Terminus::PCR,
            bottom_chemistry: Terminus::PCR,
        }
    }
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub struct AnnealCandidate {
    pub id: String,
    pub duplex: Duplex,
    pub paired_start: usize,
    pub paired_end: usize,
    pub paired_bases: usize,
    pub mappings: Vec<StrandMapping>,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub struct AnnealReport {
    pub top_sequence_blake3: String,
    pub bottom_sequence_blake3: String,
    pub options: AnnealOptions,
    pub engine_version: String,
    pub complete: bool,
    pub comparisons: usize,
    pub unique: bool,
    pub candidates: Vec<AnnealCandidate>,
    pub warnings: Vec<String>,
}

/// Enumerate every antiparallel offset whose *entire* overlapping region is
/// concretely complementary. Nonoverlapping tails are retained as free overhangs.
/// No thermodynamic ranking, melting prediction or arbitrary best alignment.
pub fn anneal_oligos(
    top: &Sequence,
    bottom: &Sequence,
    options: &AnnealOptions,
) -> Result<AnnealReport, Error> {
    if top.len() > MAX_ANNEAL_OLIGO_BASES || bottom.len() > MAX_ANNEAL_OLIGO_BASES {
        return Err(Error::TooLarge);
    }
    exact_dna(top)?;
    exact_dna(bottom)?;
    if options.min_paired_bases == 0
        || options.min_paired_bases > MAX_ANNEAL_OLIGO_BASES
        || options.max_comparisons == 0
        || options.max_comparisons > MAX_ANNEAL_COMPARISONS
        || options.max_candidates == 0
        || options.max_candidates > MAX_ANNEAL_CANDIDATES
    {
        return Err(invalid("invalid oligo-annealing limits"));
    }
    let top_hash = blake3::hash(top.as_str().as_bytes()).to_hex().to_string();
    let bottom_hash = blake3::hash(bottom.as_str().as_bytes())
        .to_hex()
        .to_string();
    let mut report=AnnealReport{top_sequence_blake3:top_hash,bottom_sequence_blake3:bottom_hash,options:options.clone(),engine_version:ANNEAL_ENGINE_VERSION.into(),complete:true,comparisons:0,unique:false,candidates:vec![],warnings:vec!["Exhaustive concrete antiparallel overlap model under the recorded bounds; every compared overlapping base must pair. Stability, temperature, concentration, hairpins, self-association, mismatches and reaction kinetics are not modeled. Input chemistry is an explicit caller assumption, not inferred from imported metadata.".into()]};
    let n = top.len() as i64;
    let m = bottom.len() as i64;
    if options.min_paired_bases > top.len().min(bottom.len()) {
        report
            .warnings
            .push("Neither supplied oligo pair can meet the minimum paired length.".into());
        return Ok(report);
    }
    let min = options.min_paired_bases as i64;
    let first = min - m;
    let last = n - min;
    for start in first..=last {
        let paired_start = start.max(0);
        let paired_end = (start + m).min(n);
        let mut complementary = true;
        for coordinate in paired_start..paired_end {
            if report.comparisons == options.max_comparisons {
                report.complete = false;
                report.warnings.push("Comparison budget exhausted; the current offset and remaining offsets were not accepted as complete evidence.".into());
                complementary = false;
                break;
            }
            report.comparisons += 1;
            let bottom_index = (m - 1 - (coordinate - start)) as usize;
            if top.as_str().as_bytes()[coordinate as usize]
                != collagene_model::complement(
                    bottom.as_str().as_bytes()[bottom_index],
                    Alphabet::Dna,
                )
            {
                complementary = false;
                break;
            }
        }
        if !report.complete {
            break;
        }
        if !complementary {
            continue;
        }
        let duplex = Duplex {
            top: top.clone(),
            bottom: bottom.clone(),
            bottom_start: start,
            topology: Topology::Linear,
            left: Terminus {
                five_prime: options.top_chemistry.five_prime,
                three_prime: options.bottom_chemistry.three_prime,
            },
            right: Terminus {
                five_prime: options.bottom_chemistry.five_prime,
                three_prime: options.top_chemistry.three_prime,
            },
            nicks: vec![],
        };
        duplex.validate()?;
        let identifier = serde_json::to_vec(&(
            ANNEAL_ENGINE_VERSION,
            &report.top_sequence_blake3,
            &report.bottom_sequence_blake3,
            options,
            start,
        ))
        .map_err(|e| Error::Format(e.to_string()))?;
        let candidate = AnnealCandidate {
            id: format!("anneal-{}", blake3::hash(&identifier).to_hex()),
            duplex,
            paired_start: paired_start as usize,
            paired_end: paired_end as usize,
            paired_bases: (paired_end - paired_start) as usize,
            mappings: vec![
                identity_mapping(Strand::Top, top.len()),
                identity_mapping(Strand::Bottom, bottom.len()),
            ],
        };
        report.candidates.push(candidate);
        if report.candidates.len() == options.max_candidates && start != last {
            report.complete = false;
            report.warnings.push("Candidate limit reached; later offsets remain unsearched. No uniqueness claim is available.".into());
            break;
        }
    }
    report.unique = report.complete && report.candidates.len() == 1;
    if report.candidates.len() > 1 {
        report.warnings.push("Multiple geometries satisfy the minimum overlap. No longest or preferred offset was selected; choosing a candidate does not establish uniqueness.".into());
    }
    if report.candidates.is_empty() && report.complete {
        report
            .warnings
            .push("No fully complementary offset satisfies the requested minimum overlap.".into());
    }
    bounded_chemistry_json(&report)?;
    Ok(report)
}
