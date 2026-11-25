//! Editable sequence reuse only. Anchors describe an explicit selected reference,
//! never physical oligo identity, current binding, or successful amplification.
use collagene_core::{
    Alphabet, Molecule, Sequence, Topology,
    analysis::{Primer, base_mask},
    primers::MAX_OLIGO,
    run_lineage::EventBinding,
};
use serde::{Deserialize, Serialize};

pub const MAX_NAME_BYTES: usize = 256;
pub const MAX_RAW_BASE_BYTES: usize = MAX_OLIGO * 8;

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SelectionAnchor {
    pub source: EventBinding,
    /// Source-forward normalized start, strictly less than source length.
    pub start: usize,
    /// Unrolled exclusive boundary; may exceed source length on a circle.
    pub end: usize,
    pub reverse: bool,
    /// Original selected DNA in this oligo's 5′→3′ orientation, excluding tail.
    pub reference: String,
}
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SelectionSeed {
    pub source: EventBinding,
    pub start: usize,
    pub end: usize,
    /// Original source-forward selected DNA; end may pass a circular origin once.
    pub reference: String,
}
#[derive(Clone, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct OligoDraft {
    pub name: String,
    pub annealing: String,
    pub tail: String,
    pub anchor: Option<SelectionAnchor>,
}
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PairDraft {
    pub forward: OligoDraft,
    pub reverse: OligoDraft,
}
impl Default for PairDraft {
    fn default() -> Self {
        Self {
            forward: OligoDraft {
                name: "Forward".into(),
                ..OligoDraft::default()
            },
            reverse: OligoDraft {
                name: "Reverse".into(),
                ..OligoDraft::default()
            },
        }
    }
}
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Comparison {
    pub reference: String,
    pub domain: String,
    /// Disjoint IUPAC call sets, zero-based in the oligo's 5′→3′ domain.
    pub mismatches: Vec<usize>,
    /// Compatible call sets with at least one ambiguous call, including N/N.
    pub ambiguous: Vec<usize>,
}
fn hex(s: &str) -> bool {
    s.len() == 64
        && s.bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}
fn binding_guard(b: &EventBinding) -> Result<(), String> {
    if b.document_id.is_empty()
        || b.document_id.len() > 128
        || !b
            .document_id
            .bytes()
            .all(|v| v.is_ascii_alphanumeric() || matches!(v, b'_' | b'-' | b'.' | b':'))
        || !hex(&b.revision)
        || !b.event_id.strip_prefix("event-").is_some_and(hex)
    {
        return Err("Selection requires an exact native document/event/revision binding.".into());
    }
    Ok(())
}
pub(crate) fn selected_guard(
    source: &EventBinding,
    start: usize,
    end: usize,
    reference: &str,
) -> Result<(), String> {
    binding_guard(source)?;
    if start >= end || end - start > MAX_OLIGO || reference.len() != end - start {
        return Err("Selection must be one exact contiguous domain of 1–200 bases.".into());
    }
    let sequence = Sequence::dna(reference).map_err(|e| e.to_string())?;
    if sequence.as_str() != reference {
        return Err("Selection evidence must retain its canonical original DNA.".into());
    }
    Ok(())
}

/// Retain actual selected DNA in oligo 5′→3′ orientation, without a binding scan.
/// Start is normalized within the source; end is an unrolled boundary. The
/// nonempty domain is at most one circular turn and at most 200 bases. The host
/// resolves exact event membership; this verifies its shape and actual revision.
pub fn anchor_reference(
    molecule: &Molecule,
    source: &EventBinding,
    start: usize,
    end: usize,
    reverse: bool,
) -> Result<String, String> {
    binding_guard(source)?;
    let n = molecule.sequence.len();
    if start >= n
        || end <= start
        || end - start > n
        || end - start > MAX_OLIGO
        || molecule.sequence.alphabet() != Alphabet::Dna
        || molecule.topology == Topology::Linear && end > n
    {
        return Err("Choose one nonempty DNA domain of at most 200 bases and one circular turn, with its start inside the source.".into());
    }
    collagene_core::model::bounds::model_preflight(&[molecule]).map_err(|e| e.to_string())?;
    molecule.validate().map_err(|e| e.to_string())?;
    if molecule.revision() != source.revision {
        return Err("The selected source changed; resolve its exact retained event again.".into());
    }
    let calls = molecule.sequence.as_str().as_bytes();
    let selected: String = (0..end - start)
        .map(|i| calls[(start + i) % n] as char)
        .collect();
    if reverse {
        Sequence::dna(&selected)
            .and_then(|sequence| sequence.reverse_complement())
            .map(|sequence| sequence.as_str().to_owned())
            .map_err(|e| e.to_string())
    } else {
        Ok(selected)
    }
}
impl SelectionSeed {
    /// The host supplies the actual selected event binding. This checks its shape
    /// and actual molecular revision; only the workspace can verify event ancestry.
    pub fn capture(
        molecule: &Molecule,
        source: &EventBinding,
        start: usize,
        end: usize,
    ) -> Result<Self, String> {
        let reference = anchor_reference(molecule, source, start, end, false)?;
        Ok(Self {
            source: source.clone(),
            start,
            end,
            reference,
        })
    }
    pub fn oligo(&self, name: &str, reverse: bool, tail: &str) -> Result<OligoDraft, String> {
        selected_guard(&self.source, self.start, self.end, &self.reference)?;
        // Guard caller text before copying it into the editable draft.
        text_guard(name, &self.reference, tail)?;
        let reference = if reverse {
            Sequence::dna(&self.reference)
                .and_then(|s| s.reverse_complement())
                .map_err(|e| e.to_string())?
                .as_str()
                .to_owned()
        } else {
            self.reference.clone()
        };
        let mut draft = OligoDraft {
            name: name.into(),
            annealing: reference.clone(),
            tail: tail.into(),
            anchor: Some(SelectionAnchor {
                source: self.source.clone(),
                start: self.start,
                end: self.end,
                reverse,
                reference,
            }),
        };
        let primer = draft.primer()?;
        draft.annealing = primer.annealing.as_str().into();
        draft.tail = primer.tail;
        Ok(draft)
    }
}
fn text_guard(name: &str, domain: &str, tail: &str) -> Result<(), String> {
    if name.is_empty()
        || name.trim() != name
        || name.len() > MAX_NAME_BYTES
        || name.chars().any(char::is_control)
    {
        return Err(
            "Give the primer a trimmed name of at most 256 bytes without control characters."
                .into(),
        );
    }
    if domain
        .len()
        .checked_add(tail.len())
        .is_none_or(|n| n > MAX_RAW_BASE_BYTES)
        || !domain.is_ascii()
        || !tail.is_ascii()
    {
        return Err("Primer text exceeds the bounded DNA input profile.".into());
    }
    Ok(())
}
impl OligoDraft {
    pub fn primer(&self) -> Result<Primer, String> {
        text_guard(&self.name, &self.annealing, &self.tail)?;
        let p = Primer::new(self.name.clone(), &self.annealing, &self.tail)
            .map_err(|e| e.to_string())?;
        if p.annealing
            .len()
            .checked_add(p.tail.len())
            .is_none_or(|n| n > MAX_OLIGO)
        {
            return Err(
                "The complete ordered primer, including its 5′ tail, must be at most 200 bases."
                    .into(),
            );
        }
        Ok(p)
    }
    /// Literal comparison to retained selection evidence, not a new binding search.
    /// A length edit leaves placement unresolved rather than guessing an alignment.
    pub fn comparison(&self) -> Result<Option<Comparison>, String> {
        let Some(anchor) = &self.anchor else {
            return Ok(None);
        };
        selected_guard(&anchor.source, anchor.start, anchor.end, &anchor.reference)?;
        let p = self.primer()?;
        if p.annealing.len() != anchor.reference.len() {
            return Err("Annealing length changed; select the intended domain again. Placement is unresolved.".into());
        }
        let mut mismatches = Vec::new();
        let mut ambiguous = Vec::new();
        for (i, (r, d)) in anchor
            .reference
            .bytes()
            .zip(p.annealing.as_str().bytes())
            .enumerate()
        {
            let a = base_mask(r);
            let b = base_mask(d);
            if a & b == 0 {
                mismatches.push(i);
            } else if a.count_ones() != 1 || b.count_ones() != 1 {
                ambiguous.push(i);
            }
        }
        Ok(Some(Comparison {
            reference: anchor.reference.clone(),
            domain: p.annealing.as_str().into(),
            mismatches,
            ambiguous,
        }))
    }
}
impl PairDraft {
    pub fn primers(&self) -> Result<(Primer, Primer), String> {
        Ok((self.forward.primer()?, self.reverse.primer()?))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use collagene_core::{Topology, workspace::Workspace};
    fn source(dna: &str) -> (Molecule, EventBinding) {
        source_topology(dna, Topology::Linear)
    }
    fn source_topology(dna: &str, topology: Topology) -> (Molecule, EventBinding) {
        let m = Molecule::dna("source", dna, topology).unwrap();
        let w = Workspace::from_documents("primer draft", vec![m.clone()], 0).unwrap();
        let d = &w.documents[0];
        (
            m,
            EventBinding {
                document_id: d.id.clone(),
                event_id: d.head_event.clone().unwrap(),
                revision: d.head.clone(),
            },
        )
    }
    #[test]
    fn forward_and_reverse_preserve_explicit_geometry_and_full_order() {
        let (m, b) = source("GGATGCRYTT");
        let seed = SelectionSeed::capture(&m, &b, 2, 8).unwrap();
        let f = seed.oligo("F", false, "TT").unwrap();
        let r = seed.oligo("R", true, "AA").unwrap();
        assert_eq!(f.primer().unwrap().ordered_sequence(), "TTATGCRY");
        assert_eq!(r.primer().unwrap().ordered_sequence(), "AARYGCAT");
        assert_eq!(r.anchor.as_ref().unwrap().reference, "RYGCAT");
        assert!(r.anchor.as_ref().unwrap().reverse);
        assert_eq!(r.anchor.as_ref().unwrap().start, 2);
        assert_eq!(r.anchor.as_ref().unwrap().end, 8);
        assert_eq!(seed.reference, "ATGCRY");
        assert_eq!(m.sequence.as_str(), "GGATGCRYTT");
    }
    #[test]
    fn comparison_distinguishes_disjoint_and_compatible_uncertainty_in_oligo_coordinates() {
        let (m, b) = source("ARNT");
        let seed = SelectionSeed::capture(&m, &b, 0, 4).unwrap();
        let mut f = seed.oligo("F", false, "").unwrap();
        f.annealing = "AGNC".into();
        let c = f.comparison().unwrap().unwrap();
        assert_eq!(c.mismatches, vec![3]);
        assert_eq!(c.ambiguous, vec![1, 2]);
        let mut r = seed.oligo("R", true, "CC").unwrap();
        assert_eq!(r.annealing, "ANYT");
        r.annealing = "TNGT".into();
        let c = r.comparison().unwrap().unwrap();
        assert_eq!(c.reference, "ANYT");
        assert_eq!(c.mismatches, vec![0, 2]);
        assert_eq!(c.ambiguous, vec![1]);
        r.annealing = "TNG".into();
        assert!(r.comparison().unwrap_err().contains("unresolved"));
        r.anchor = None;
        assert_eq!(r.comparison().unwrap(), None);
    }
    #[test]
    fn stale_malformed_and_noncontiguous_selections_are_refused() {
        let (m, mut b) = source("ACGT");
        b.revision = "0".repeat(64);
        assert!(SelectionSeed::capture(&m, &b, 0, 4).is_err());
        let (m, b) = source("ACGT");
        for (a, z) in [(0, 0), (3, 2), (0, 5), (usize::MAX, usize::MAX)] {
            assert!(SelectionSeed::capture(&m, &b, a, z).is_err());
        }
        let mut bad = b.clone();
        bad.event_id = "0".repeat(64);
        assert!(SelectionSeed::capture(&m, &bad, 0, 4).is_err());
        let mut bad = b;
        bad.document_id = "not a native identifier".into();
        assert!(SelectionSeed::capture(&m, &bad, 0, 4).is_err());
        let rna = Molecule {
            sequence: Sequence::new(Alphabet::Rna, "ACGU").unwrap(),
            ..m
        };
        assert!(SelectionSeed::capture(&rna, &bad, 0, 4).is_err());
    }
    #[test]
    fn full_order_and_raw_text_are_bounded_before_draft_copy_and_parse() {
        let (m, b) = source(&"A".repeat(MAX_OLIGO));
        let seed = SelectionSeed::capture(&m, &b, 0, MAX_OLIGO).unwrap();
        assert!(seed.oligo("F", false, "").is_ok());
        assert!(seed.oligo("F", false, "A").is_err());
        assert!(
            seed.oligo("F", false, &" ".repeat(MAX_RAW_BASE_BYTES))
                .is_err()
        );
        assert!(seed.oligo(&"x".repeat(257), false, "").is_err());
        let mut pair = PairDraft::default();
        assert_eq!(pair.forward.name, "Forward");
        assert_eq!(pair.reverse.name, "Reverse");
        assert!(pair.primers().is_err());
        pair.forward = SelectionSeed::capture(&m, &b, 0, 8)
            .unwrap()
            .oligo("Forward", false, "")
            .unwrap();
        pair.reverse = pair.forward.clone();
        assert!(pair.primers().is_ok());
        pair.forward.annealing = "A".repeat(MAX_RAW_BASE_BYTES + 1);
        assert!(pair.primers().is_err());
    }

    #[test]
    fn circular_anchor_uses_normalized_start_unrolled_end_and_oriented_calls() {
        let (m, b) = source_topology("ACGTAAGC", Topology::Circular);
        let seed = SelectionSeed::capture(&m, &b, 6, 11).unwrap();
        assert_eq!(seed.reference, "GCACG");
        let f = seed.oligo("forward", false, "TT").unwrap();
        let r = seed.oligo("reverse", true, "GG").unwrap();
        assert_eq!(f.primer().unwrap().ordered_sequence(), "TTGCACG");
        assert_eq!(r.primer().unwrap().ordered_sequence(), "GGCGTGC");
        let anchor = r.anchor.unwrap();
        assert_eq!((anchor.start, anchor.end, anchor.reverse), (6, 11, true));
        assert_eq!(anchor_reference(&m, &b, 6, 11, true).unwrap(), "CGTGC");
        assert_eq!(anchor_reference(&m, &b, 3, 11, false).unwrap(), "TAAGCACG");
    }

    #[test]
    fn one_turn_and_oligo_limits_reject_unnormalized_or_multiturn_coordinates() {
        let (m, b) = source_topology("ACGTAAGC", Topology::Circular);
        for (start, end) in [(8, 9), (6, 6), (7, 16), (usize::MAX - 1, usize::MAX)] {
            assert!(anchor_reference(&m, &b, start, end, false).is_err());
        }
        let (linear, lb) = source("ACGTAAGC");
        assert!(anchor_reference(&linear, &lb, 6, 11, false).is_err());
        let (large, large_binding) = source_topology(&"A".repeat(201), Topology::Circular);
        assert_eq!(
            anchor_reference(&large, &large_binding, 200, 400, false)
                .unwrap()
                .len(),
            200
        );
        assert!(anchor_reference(&large, &large_binding, 200, 401, false).is_err());
        let mut stale = b;
        stale.revision = "0".repeat(64);
        assert!(anchor_reference(&m, &stale, 6, 11, true).is_err());
    }
}
