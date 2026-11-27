//! Editable sequence reuse only. Anchors describe an explicit selected reference,
//! never physical oligo identity, current binding, or successful amplification.
use collagene_core::{
    Alphabet, Location, Molecule, Sequence, Topology,
    analysis::{Primer, base_mask},
    primers::{BindingHit, MAX_OLIGO},
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
    #[serde(deserialize_with = "required_anchor")]
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

    /// Capture comparison context for an existing oligo without replacing any
    /// editable text. The host must resolve the seed's actual retained event.
    /// This is neither an alignment guess nor a binding-site search: the domain
    /// and explicitly selected reference must have equal canonical lengths.
    pub fn compare_existing(
        &self,
        current: &OligoDraft,
        reverse: bool,
    ) -> Result<OligoDraft, String> {
        selected_guard(&self.source, self.start, self.end, &self.reference)?;
        let primer = current.primer()?;
        if primer.annealing.len() != self.reference.len() {
            return Err("Select exactly as many reference bases as this oligo's annealing domain. The existing oligo was not changed.".into());
        }
        let reference = if reverse {
            Sequence::dna(&self.reference)
                .and_then(|s| s.reverse_complement())
                .map_err(|e| e.to_string())?
                .as_str()
                .to_owned()
        } else {
            self.reference.clone()
        };
        Ok(OligoDraft {
            name: current.name.clone(),
            annealing: current.annealing.clone(),
            tail: current.tail.clone(),
            anchor: Some(SelectionAnchor {
                source: self.source.clone(),
                start: self.start,
                end: self.end,
                reverse,
                reference,
            }),
        })
    }
}

/// Retain one explicitly chosen binding placement as editable reference context.
/// The host guards the original report, options, epoch and exact source-event
/// membership. This checks the chosen hit's actual geometry and per-base facts
/// without rescanning the reference or retaining a claim of complete/unique
/// binding. Subsequent draft edits still require explicit Apply and analysis.
pub fn capture_binding_hit(
    molecule: &Molecule,
    source: &EventBinding,
    primer: &Primer,
    hit: &BindingHit,
) -> Result<OligoDraft, String> {
    text_guard(&primer.name, primer.annealing.as_str(), &primer.tail)?;
    if primer.annealing.alphabet() != Alphabet::Dna
        || primer.annealing.len() > MAX_OLIGO
        || hit.length != primer.annealing.len()
        || hit.mismatches.len() > MAX_OLIGO
        || hit.ambiguous_positions.len() > MAX_OLIGO
    {
        return Err("The selected binding hit does not match a bounded DNA primer.".into());
    }
    let current = OligoDraft {
        name: primer.name.clone(),
        annealing: primer.annealing.as_str().into(),
        tail: primer.tail.clone(),
        anchor: None,
    };
    current.primer()?;
    let end = hit
        .start
        .checked_add(hit.length)
        .ok_or("The selected binding interval exceeds the coordinate bounds.")?;
    let seed = SelectionSeed::capture(molecule, source, hit.start, end)?;
    let n = molecule.sequence.len();
    let crosses_origin = end > n;
    let location = if crosses_origin {
        Location::Join(vec![
            Location::range(hit.start, n),
            Location::range(0, end - n),
        ])
    } else {
        Location::range(hit.start, end)
    };
    let location = if hit.reverse {
        Location::Complement(Box::new(location))
    } else {
        location
    };
    if hit.crosses_origin != crosses_origin || hit.location != location {
        return Err(
            "The selected binding hit's location or circular crossing is inconsistent.".into(),
        );
    }
    let compared = seed.compare_existing(&current, hit.reverse)?;
    let comparison = compared
        .comparison()?
        .ok_or("The selected binding reference was not captured.")?;
    let three_prime_exact = (0..hit.length)
        .rev()
        .take_while(|index| {
            !comparison.mismatches.contains(index) && !comparison.ambiguous.contains(index)
        })
        .count();
    if hit.mismatches != comparison.mismatches
        || hit.ambiguous_positions != comparison.ambiguous
        || hit.three_prime_exact != three_prime_exact
        || hit.certain != (comparison.mismatches.is_empty() && comparison.ambiguous.is_empty())
    {
        return Err("The selected binding hit's mismatch, ambiguity or 3′ evidence does not match its actual reference.".into());
    }
    Ok(compared)
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

    #[test]
    fn compare_existing_preserves_raw_oligo_and_uses_explicit_direction() {
        let (m, b) = source("GGAGNTCC");
        let seed = SelectionSeed::capture(&m, &b, 2, 6).unwrap();
        for (reverse, calls, reference, mismatch, ambiguous) in [
            (false, "t g n t", "AGNT", vec![0], vec![2]),
            (true, "t n c t", "ANCT", vec![0], vec![1]),
        ] {
            let original = OligoDraft {
                name: "Intentional edits".into(),
                annealing: calls.into(),
                tail: "g g".into(),
                anchor: None,
            };
            let compared = seed.compare_existing(&original, reverse).unwrap();
            assert_eq!(compared.name, original.name);
            assert_eq!(compared.annealing, original.annealing);
            assert_eq!(compared.tail, original.tail);
            assert!(original.anchor.is_none());
            let comparison = compared.comparison().unwrap().unwrap();
            assert_eq!(comparison.reference, reference);
            assert_eq!(comparison.mismatches, mismatch);
            assert_eq!(comparison.ambiguous, ambiguous);
            let anchor = compared.anchor.unwrap();
            assert_eq!((anchor.start, anchor.end, anchor.reverse), (2, 6, reverse));
            assert_eq!(anchor.source, b);
        }
    }

    #[test]
    fn compare_existing_replaces_only_anchor_and_handles_wrapped_reverse_context() {
        let (m, b) = source_topology("ACGTAAGC", Topology::Circular);
        let seed = SelectionSeed::capture(&m, &b, 6, 11).unwrap();
        let old = SelectionSeed::capture(&m, &b, 0, 5).unwrap();
        let mut original = old.oligo("Retained mismatch", false, "TT").unwrap();
        original.annealing = "c g t N c".into();
        let compared = seed.compare_existing(&original, true).unwrap();
        assert_eq!(compared.annealing, "c g t N c");
        assert_eq!(compared.primer().unwrap().ordered_sequence(), "TTCGTNC");
        let comparison = compared.comparison().unwrap().unwrap();
        assert_eq!(comparison.reference, "CGTGC");
        assert_eq!(comparison.mismatches, Vec::<usize>::new());
        assert_eq!(comparison.ambiguous, vec![3]);
        assert_eq!(original.anchor.unwrap().start, 0);
        let anchor = compared.anchor.unwrap();
        assert_eq!((anchor.start, anchor.end, anchor.reverse), (6, 11, true));
    }

    #[test]
    fn compare_existing_refuses_length_invalid_text_and_malformed_seeds_without_repair() {
        let (m, b) = source("AGNT");
        let seed = SelectionSeed::capture(&m, &b, 0, 4).unwrap();
        let original = OligoDraft {
            name: "Existing".into(),
            annealing: "t g n t".into(),
            tail: "TT".into(),
            anchor: None,
        };
        for domain in ["AGT", "", "AG?T"] {
            let mut changed = original.clone();
            changed.annealing = domain.into();
            assert!(seed.compare_existing(&changed, false).is_err());
            assert_eq!(changed.annealing, domain);
            assert!(changed.anchor.is_none());
        }
        let mut oversized = original.clone();
        oversized.tail = "A".repeat(MAX_OLIGO);
        assert!(seed.compare_existing(&oversized, false).is_err());
        let mut unnamed = original.clone();
        unnamed.name = " ".into();
        assert!(seed.compare_existing(&unnamed, false).is_err());
        let mut malformed = seed.clone();
        malformed.reference = "agnt".into();
        assert!(malformed.compare_existing(&original, false).is_err());
        malformed = seed.clone();
        malformed.source.event_id = "unbound".into();
        assert!(malformed.compare_existing(&original, false).is_err());
        assert_eq!(original.annealing, "t g n t");
        assert!(original.anchor.is_none());
    }

    #[test]
    fn chosen_binding_hits_capture_actual_forward_and_wrapped_reverse_reference() {
        let (m, b) = source("GGAGNTCC");
        let primer = Primer::new("Forward hit", "TGNT", "GG").unwrap();
        let hit = BindingHit {
            start: 2,
            length: 4,
            reverse: false,
            crosses_origin: false,
            certain: false,
            mismatches: vec![0],
            ambiguous_positions: vec![2],
            three_prime_exact: 1,
            location: Location::range(2, 6),
        };
        let captured = capture_binding_hit(&m, &b, &primer, &hit).unwrap();
        assert_eq!(captured.primer().unwrap().ordered_sequence(), "GGTGNT");
        assert_eq!(captured.anchor.unwrap().reference, "AGNT");

        let (circle, cb) = source_topology("ACGTAAGC", Topology::Circular);
        let primer = Primer::new("Wrapped reverse", "CATNC", "GG").unwrap();
        let hit = BindingHit {
            start: 6,
            length: 5,
            reverse: true,
            crosses_origin: true,
            certain: false,
            mismatches: vec![1],
            ambiguous_positions: vec![3],
            three_prime_exact: 1,
            location: Location::Complement(Box::new(Location::Join(vec![
                Location::range(6, 8),
                Location::range(0, 3),
            ]))),
        };
        let captured = capture_binding_hit(&circle, &cb, &primer, &hit).unwrap();
        assert_eq!(captured.primer().unwrap().ordered_sequence(), "GGCATNC");
        let anchor = captured.anchor.unwrap();
        assert_eq!((anchor.start, anchor.end, anchor.reverse), (6, 11, true));
        assert_eq!(anchor.reference, "CGTGC");
    }

    #[test]
    fn chosen_binding_hit_refuses_stale_source_and_inconsistent_placement_facts() {
        let (m, b) = source("GGAGNTCC");
        let primer = Primer::new("Chosen", "TGNT", "").unwrap();
        let hit = BindingHit {
            start: 2,
            length: 4,
            reverse: false,
            crosses_origin: false,
            certain: false,
            mismatches: vec![0],
            ambiguous_positions: vec![2],
            three_prime_exact: 1,
            location: Location::range(2, 6),
        };
        let variants = [
            BindingHit {
                start: usize::MAX,
                ..hit.clone()
            },
            BindingHit {
                start: 7,
                ..hit.clone()
            },
            BindingHit {
                length: 3,
                ..hit.clone()
            },
            BindingHit {
                reverse: true,
                ..hit.clone()
            },
            BindingHit {
                crosses_origin: true,
                ..hit.clone()
            },
            BindingHit {
                certain: true,
                ..hit.clone()
            },
            BindingHit {
                mismatches: vec![1],
                ..hit.clone()
            },
            BindingHit {
                ambiguous_positions: vec![],
                ..hit.clone()
            },
            BindingHit {
                three_prime_exact: 2,
                ..hit.clone()
            },
            BindingHit {
                location: Location::range(1, 5),
                ..hit.clone()
            },
            BindingHit {
                location: Location::Order(vec![Location::range(2, 6)]),
                ..hit.clone()
            },
        ];
        for changed in variants {
            assert!(capture_binding_hit(&m, &b, &primer, &changed).is_err());
        }
        let mut stale = b.clone();
        stale.revision = "0".repeat(64);
        assert!(capture_binding_hit(&m, &stale, &primer, &hit).is_err());
        let mut rna = primer.clone();
        rna.annealing = Sequence::new(Alphabet::Rna, "UGNU").unwrap();
        assert!(capture_binding_hit(&m, &b, &rna, &hit).is_err());
        let mut oversized = primer;
        oversized.tail = "A".repeat(MAX_OLIGO);
        assert!(capture_binding_hit(&m, &b, &oversized, &hit).is_err());
    }
}

fn required_anchor<'de, D: serde::Deserializer<'de>>(
    deserializer: D,
) -> Result<Option<SelectionAnchor>, D::Error> {
    Option::deserialize(deserializer)
}
