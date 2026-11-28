//! Ephemeral, exact-source display evidence. No sequence annotations are installed.
use crate::biology::antibody::antibody_model::Capture;
use collagene_core::antibody_numbering::domain_tracks as n;
use std::sync::Arc;
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Overlay {
    pub capture: Capture,
    pub track: Arc<n::Track>,
}
#[derive(Clone, Debug, PartialEq)]
pub struct Label {
    pub source_index: usize,
    pub text: String,
    pub left: f64,
    pub width: f64,
}
impl Overlay {
    pub fn new(capture: &Capture, report: &n::Report, candidate: usize) -> Result<Self, String> {
        capture.validate_report(report)?;
        let track = report
            .tracks
            .get(candidate)
            .ok_or("Choose a retained antibody candidate.")?;
        Ok(Self {
            capture: capture.clone(),
            track: Arc::new(track.clone()),
        })
    }
    pub fn identity(&self) -> String {
        sequence_identity(self.capture.epoch, &self.capture.source)
    }
    pub fn matches(&self, identity: &str, revision: &str, length: usize) -> bool {
        self.identity() == identity
            && self.capture.source.revision == revision
            && self.capture.source_length == length
            && self.track.source_end <= length
    }
    pub fn at(&self, index: usize) -> Option<&n::Residue> {
        self.track
            .residues
            .get(index.checked_sub(self.track.first_supported_source_index)?)
    }
    /// Sparse labels in source order, never floating-point position order. All
    /// labels remain available in the readout/inspector. Return nonoverlapping
    /// pixel intervals within this logical row; no biology is inferred here.
    pub fn labels(&self, start: usize, end: usize, cell: f64) -> Vec<Label> {
        if start >= end || !cell.is_finite() || cell <= 0. || end > self.capture.source_length {
            return vec![];
        }
        let available = (end - start) as f64 * cell;
        let mut right = f64::NEG_INFINITY;
        let mut out = vec![];
        for index in
            start.max(self.track.first_supported_source_index)..end.min(self.track.source_end)
        {
            if index % 5 != 0 && index != self.track.first_supported_source_index {
                continue;
            }
            let Some(row) = self.at(index) else {
                continue;
            };
            let text = row
                .assignment
                .as_ref()
                .map(|a| a.position.to_string())
                .unwrap_or_else(|| "?".into());
            let width = text.len() as f64 * 6.;
            if width > available {
                continue;
            }
            let left =
                (((index - start) as f64 + 0.5) * cell - width / 2.).clamp(0., available - width);
            if left < right + 6. {
                continue;
            }
            out.push(Label {
                source_index: index,
                text,
                left,
                width,
            });
            right = left + width;
        }
        out
    }
    pub fn description(&self, index: usize) -> String {
        let Some(row) = self.at(index) else {
            return format!(
                "Source {}: no supported assignment in this candidate; not proof of an outside-domain residue.",
                index + 1
            );
        };
        let region = row
            .region
            .as_ref()
            .map(crate::biology::antibody::antibody_overlay::region)
            .unwrap_or("unresolved region");
        match &row.assignment {
            Some(a) => format!(
                "Source {} {} · conditional {} {} · {} {} · {}",
                index + 1,
                row.residue,
                scheme_label(&self.track.numbering_scheme),
                a.position,
                definition_label(&self.track.region_definition),
                region,
                origin(&a.evidence)
            ),
            None => format!(
                "Source {} {} · unresolved {} position · {} {}",
                index + 1,
                row.residue,
                scheme_label(&self.track.numbering_scheme),
                definition_label(&self.track.region_definition),
                region
            ),
        }
    }
}
/// Unambiguous exact-event identity for the shared sequence viewer. IDs remain
/// arbitrary strings; delimiter characters cannot create identity collisions.
pub fn sequence_identity(
    epoch: u64,
    source: &collagene_core::material::lineage::EventBinding,
) -> String {
    serde_json::to_string(&(epoch, source)).expect("integer/string event identity")
}
pub fn region(r: &n::Region) -> &'static str {
    match r {
        n::Region::Fr1 => "FR1",
        n::Region::Cdr1 => "CDR1",
        n::Region::Fr2 => "FR2",
        n::Region::Cdr2 => "CDR2",
        n::Region::Fr3 => "FR3",
        n::Region::Cdr3 => "CDR3",
        n::Region::Fr4 => "FR4",
    }
}
pub fn origin(e: &n::AssignmentEvidence) -> &'static str {
    match e {
        n::AssignmentEvidence::InvariantVReference => "invariant V-reference correspondence",
        n::AssignmentEvidence::ConditionalLoopRule => "conditional loop rule",
        n::AssignmentEvidence::JFr4Reference => "matched J/FR4 reference",
        n::AssignmentEvidence::PublishedFrameworkCorrespondence => {
            "published invariant framework correspondence"
        }
        n::AssignmentEvidence::PublishedWindowRule => "published conditional window rule",
    }
}
pub fn color(r: &n::Region) -> Option<&'static str> {
    match r {
        n::Region::Cdr1 => Some("rgba(52,134,184,0.22)"),
        n::Region::Cdr2 => Some("rgba(195,144,38,0.25)"),
        n::Region::Cdr3 => Some("rgba(140,97,175,0.22)"),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use collagene_core::{antibody_numbering::domain_tracks, material::lineage::EventBinding};
    use collagene_model::{Alphabet, Molecule, Sequence, Topology};
    const V: &str = "EVQLLESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYC";
    const J: &str = "WGQGTLVTVSS";
    fn fixture(sequence: &str) -> (Capture, domain_tracks::Report) {
        let m = Molecule {
            name: "Conditional overlay fixture".into(),
            sequence: Sequence::new(Alphabet::Protein, sequence).unwrap(),
            topology: Topology::Linear,
            features: vec![],
            metadata: Default::default(),
            original_genbank: None,
        };
        let (o, l) = crate::biology::antibody::antibody_model::preset(
            collagene_core::antibody_numbering::v_mapping::ReferenceSelection::human_all(),
            300,
            1,
            60,
            32,
        );
        let c = Capture::prepare(
            7,
            EventBinding {
                document_id: "document-a".into(),
                event_id: "event-a".into(),
                revision: m.revision(),
            },
            &m,
            o.clone(),
            l.clone(),
            test_selection(),
        )
        .unwrap();
        (
            c,
            domain_tracks::analyze(&m, &o, &l, &test_selection(), &mut || false).unwrap(),
        )
    }
    #[test]
    fn publication_binds_candidate_document_event_epoch_revision_and_length() {
        let (c, r) = fixture(&format!("{V}ARDYGMDV{J}"));
        let o = Overlay::new(&c, &r, 0).unwrap();
        assert!(o.matches(&o.identity(), &c.source.revision, c.source_length));
        for (epoch, document, event) in [
            (7, "document-b", "event-a"),
            (7, "document-a", "event-b"),
            (8, "document-a", "event-a"),
        ] {
            let mut source = c.source.clone();
            source.document_id = document.into();
            source.event_id = event.into();
            assert!(!o.matches(
                &sequence_identity(epoch, &source),
                &c.source.revision,
                c.source_length
            ));
        }
        let mut a = c.source.clone();
        a.document_id = "a:b".into();
        a.event_id = "c".into();
        let mut b = c.source.clone();
        b.document_id = "a".into();
        b.event_id = "b:c".into();
        assert_ne!(sequence_identity(7, &a), sequence_identity(7, &b));
        assert!(!o.matches(&o.identity(), "stale", c.source_length));
        assert!(!o.matches(&o.identity(), &c.source.revision, c.source_length + 1));
        assert!(Overlay::new(&c, &r, usize::MAX).is_err());
        let mut forged = r;
        forged.tracks[0].residues[0].assignment = None;
        assert!(Overlay::new(&c, &forged, 0).is_err());
    }
    #[test]
    fn sparse_labels_fit_rows_without_sorting_or_losing_exact_insertion_syntax() {
        let (c, r) = fixture(&format!("{V}{}{J}", "A".repeat(40)));
        let o = Overlay::new(&c, &r, 0).unwrap();
        for columns in (10..=100).step_by(10) {
            for start in (0..c.source_length).step_by(columns) {
                let end = (start + columns).min(c.source_length);
                let labels = o.labels(start, end, 8.5);
                assert!(
                    labels
                        .windows(2)
                        .all(|w| w[0].source_index < w[1].source_index
                            && w[0].left + w[0].width + 6. <= w[1].left)
                );
                for l in labels {
                    assert!(l.source_index >= start && l.source_index < end);
                    assert!(l.left >= 0. && l.left + l.width <= (end - start) as f64 * 8.5);
                    assert_eq!(
                        l.text,
                        o.at(l.source_index)
                            .unwrap()
                            .assignment
                            .as_ref()
                            .map(|a| a.position.to_string())
                            .unwrap_or_else(|| "?".into())
                    );
                }
            }
        }
        let left = o
            .track
            .residues
            .iter()
            .find(|r| {
                r.assignment
                    .as_ref()
                    .is_some_and(|a| a.position.to_string() == "111.10")
            })
            .unwrap();
        assert!(o.description(left.source_index).contains("111.10"));
        let right = o
            .track
            .residues
            .iter()
            .filter_map(|r| r.assignment.as_ref().map(|a| a.position.to_string()))
            .filter(|p| p.starts_with("112."))
            .collect::<Vec<_>>();
        assert!(right.iter().position(|p| p == "112.10") < right.iter().position(|p| p == "112.1"));
        assert!(o.labels(0, 10, f64::NAN).is_empty());
        assert!(o.labels(0, 10, 0.).is_empty());
        assert!(o.labels(10, 10, 8.5).is_empty());
        assert!(o.labels(0, c.source_length + 1, 8.5).is_empty());
    }
    #[test]
    fn unresolved_observed_rows_never_acquire_a_number_or_cdr_color() {
        let panel: serde_json::Value = serde_json::from_str(include_str!(
            "../../../../antibody-numbering/data/rcsb-assessment/panel.json"
        ))
        .unwrap();
        let case = panel["cases"]
            .as_array()
            .unwrap()
            .iter()
            .find(|c| c["id"] == "1HZH:2")
            .unwrap();
        let (c, r) = fixture(case["sequence"].as_str().unwrap());
        let o = Overlay::new(&c, &r, 0).unwrap();
        assert_eq!(o.track.unassigned_residues, 8);
        for row in o.track.residues.iter().filter(|r| r.assignment.is_none()) {
            let text = o.description(row.source_index);
            assert!(text.contains("unresolved") && !text.contains("conditional IMGT"));
            assert!(
                o.at(row.source_index)
                    .unwrap()
                    .region
                    .as_ref()
                    .and_then(color)
                    .is_none()
            );
        }
        assert!(o.at(c.source_length).is_none());
        assert!(o.description(c.source_length).contains("not proof"));
    }
    #[test]
    fn tandem_candidates_and_truncated_prefixes_keep_absolute_source_ownership() {
        let block = format!("{V}ARDYGMDV{J}");
        let offset = block.len() + 10;
        let (c, r) = fixture(&format!("{block}GGGGSGGGGS{block}"));
        let index = r
            .tracks
            .iter()
            .position(|t| t.first_supported_source_index >= offset)
            .unwrap();
        let o = Overlay::new(&c, &r, index).unwrap();
        assert!(o.at(offset - 1).is_none());
        assert_eq!(o.at(offset).unwrap().residue, 'E');
        let (c, r) = fixture(&format!("{}ARDYGMDV{J}", &V[8..]));
        let o = Overlay::new(&c, &r, 0).unwrap();
        assert_eq!(o.track.invariant_position_one, None);
        assert!(
            o.track
                .residues
                .iter()
                .filter_map(|r| r.assignment.as_ref())
                .all(|a| a.position.number != 1)
        );
    }
}

pub fn scheme_label(s: &n::NumberingScheme) -> &'static str {
    match s {
        n::NumberingScheme::Imgt => "IMGT",
        n::NumberingScheme::Kabat1991 => "Kabat (1991)",
        n::NumberingScheme::Chothia1997 => "Chothia (1997 profiles)",
    }
}
pub fn definition_label(d: &n::RegionDefinition) -> &'static str {
    match d {
        n::RegionDefinition::Imgt => "IMGT",
        n::RegionDefinition::Kabat1991 => "Kabat (1991)",
        n::RegionDefinition::ChothiaCanonical => "Chothia (canonical charts)",
        n::RegionDefinition::North2011 => "North (2011)",
    }
}
pub fn issue_text(i: &collagene_core::antibody_numbering::alternate_tracks::Issue) -> String {
    use collagene_core::antibody_numbering::{
        alternate_tracks::WindowIssue as W, numbering_windows::Window,
    };
    let name = match i.window {
        Window::CysToTrp => "first-loop window",
        Window::MiddleLoop => "middle-loop window",
        Window::CysToJ => "junction window",
    };
    let reason = match &i.issue {
        W::MissingInvariantPositions { positions } => format!(
            "missing invariant reference positions {}",
            positions
                .iter()
                .map(ToString::to_string)
                .collect::<Vec<_>>()
                .join(", ")
        ),
        W::NonContiguousFlanks => "framework flank spacing is unresolved".into(),
        W::UnsupportedLength { residues } => {
            format!("no supported published profile for {residues} interior residues")
        }
        W::NonCanonicalAnchors => "anchor residues differ from this profile".into(),
    };
    format!("{} {name}: {reason}", scheme_label(&i.numbering_scheme))
}
pub fn diagnostics(report: &n::Report, index: usize) -> Vec<String> {
    use collagene_core::antibody_numbering::domain_candidates as d;
    let mut out = report
        .tracks
        .get(index)
        .map(|t| t.issues.iter().map(issue_text).collect::<Vec<_>>())
        .unwrap_or_default();
    out.push(crate::biology::antibody::antibody_model::reference_label(
        &report.evidence.options.references,
    ));
    out.push(format!(
        "Reference corpus identity: {}",
        report.evidence.v_evidence.corpus_revision
    ));
    for mapping in &report.evidence.v_evidence.mappings {
        out.push(format!(
            "{:?} model identity: {}",
            mapping.locus, mapping.model_revision
        ));
        out.push(format!(
            "{:?} reference admission: {} training alleles; {} omitted alleles.",
            mapping.locus,
            mapping.training_alleles.len(),
            mapping.omitted_alleles.len()
        ));
        if !mapping.omitted_alleles.is_empty() {
            out.push(format!(
                "Omitted allele identities: {}",
                mapping.omitted_alleles.join(", ")
            ));
        }
    }
    if let Some(c) = report.evidence.candidates.get(index) {
        for (name, o) in [("CDR1", &c.cdr1), ("CDR2", &c.cdr2)] {
            if let d::LoopAssignment::Unresolved { issue, .. } = o {
                let reason = match issue {
                    d::LoopIssue::MissingInvariantBoundaries { positions } => format!(
                        "missing invariant reference boundaries {}",
                        positions
                            .iter()
                            .map(ToString::to_string)
                            .collect::<Vec<_>>()
                            .join(", ")
                    ),
                    d::LoopIssue::UnsupportedLength { residues } => {
                        format!("unsupported length {residues}")
                    }
                };
                out.push(format!("IMGT {name} source context: {reason}"));
            }
        }
    }
    out
}
pub fn cdr_summaries(report: &n::Report, index: usize) -> Vec<String> {
    let Some(t) = report.tracks.get(index) else {
        return vec![];
    };
    t.cdrs
        .iter()
        .map(|c| {
            let value = match c.status {
                n::CdrStatus::Unresolved => {
                    "partial/unresolved; inspect mapped residues and diagnostics".into()
                }
                n::CdrStatus::Mapped { start, end } if start == end => "(empty)".into(),
                n::CdrStatus::Mapped { start, end } => t
                    .residues
                    .iter()
                    .filter(|r| r.source_index >= start && r.source_index < end)
                    .map(|r| r.residue)
                    .collect::<String>(),
            };
            format!(
                "{} {}: {}",
                definition_label(&t.region_definition),
                region(&c.kind),
                value
            )
        })
        .collect()
}

#[cfg(test)]
fn test_selection() -> n::Selection {
    n::Selection {
        numbering_scheme: n::NumberingScheme::Imgt,
        region_definition: n::RegionDefinition::Imgt,
    }
}
