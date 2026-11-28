//! Exact selected locations retain biological order and strand independently of highlights.
use collagene_core::{cloning::restriction::selection::footprint, material::lineage::EventBinding};
use collagene_model::{Location, Molecule};

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DirectedSelection {
    pub source: EventBinding,
    pub location: Location,
    pub spans: Vec<(usize, usize)>,
}

impl DirectedSelection {
    pub fn new(
        source: EventBinding,
        molecule: &Molecule,
        location: Location,
    ) -> Result<Self, String> {
        let spans = footprint(&location, molecule.sequence.len()).map_err(|e| e.to_string())?;
        if source.revision != molecule.revision() {
            return Err("The selected ORF source changed.".into());
        }
        Ok(Self {
            source,
            location,
            spans,
        })
    }

    pub fn contiguous(&self) -> Option<(usize, usize)> {
        if self.spans.len() == 1 {
            self.spans.first().copied()
        } else {
            None
        }
    }
}

pub fn capture(
    molecule: &Molecule,
    source: &EventBinding,
    selected_feature: Option<usize>,
    selection: Option<(usize, usize)>,
    regions: Option<&(String, Vec<(usize, usize)>)>,
    directed: Option<&DirectedSelection>,
) -> Result<Location, String> {
    if let Some(index) = selected_feature {
        let feature = molecule
            .features
            .get(index)
            .ok_or("The selected annotation changed.")?;
        footprint(&feature.location, molecule.sequence.len()).map_err(|e| e.to_string())?;
        return Ok(feature.location.clone());
    }
    if let Some(saved) = directed {
        if saved.source == *source
            && selection == saved.contiguous()
            && regions.is_some_and(|(revision, spans)| {
                *revision == source.revision && *spans == saved.spans
            })
        {
            return Ok(saved.location.clone());
        }
    }
    let spans = match regions
        .filter(|(revision, spans)| *revision == source.revision && !spans.is_empty())
    {
        Some((_, spans)) => spans.clone(),
        None => vec![selection.ok_or("Select reference bases or an exact annotation first.")?],
    };
    if spans.len() > collagene_core::cloning::restriction::selection::MAX_REGION_SPANS {
        return Err("Too many selected reference parts.".into());
    }
    let mut parts = spans
        .into_iter()
        .map(|(start, end)| Location::range(start, end))
        .collect::<Vec<_>>();
    let location = if parts.len() == 1 {
        parts.remove(0)
    } else {
        Location::Join(parts)
    };
    footprint(&location, molecule.sequence.len()).map_err(|e| e.to_string())?;
    Ok(location)
}

/// A declared origin is navigable only when its exact event is retained in the
/// protein's ancestry. This does not assert that later protein edits still match.
pub fn origin<'a>(
    frame: &'a collagene_core::documents::workspace::Workspace,
    protein: &EventBinding,
    requested: &EventBinding,
) -> Result<&'a Molecule, String> {
    frame.validate().map_err(|e| e.to_string())?;
    let events = frame
        .history
        .iter()
        .map(|e| (e.id.as_str(), e))
        .collect::<std::collections::BTreeMap<_, _>>();
    let event = events
        .get(protein.event_id.as_str())
        .ok_or("The protein event is unavailable.")?;
    if event.document_id != protein.document_id || event.revision != protein.revision {
        return Err("The protein event differs from its source link.".into());
    }
    let molecule = frame
        .revisions
        .get(&protein.revision)
        .ok_or("The protein revision is unavailable.")?;
    if molecule.sequence.alphabet() != collagene_model::Alphabet::Protein
        || molecule.metadata.get("translation:source_document") != Some(&requested.document_id)
        || molecule.metadata.get("translation:source_event") != Some(&requested.event_id)
        || molecule.metadata.get("translation:source_revision") != Some(&requested.revision)
    {
        return Err("This document has no matching declared translation origin.".into());
    }
    let mut pending = vec![protein.event_id.as_str()];
    let mut visited = std::collections::BTreeSet::new();
    while let Some(id) = pending.pop() {
        if !visited.insert(id) {
            continue;
        }
        let event = events
            .get(id)
            .ok_or("A translation ancestor is unavailable.")?;
        if id == requested.event_id {
            if event.document_id != requested.document_id || event.revision != requested.revision {
                return Err(
                    "The translation origin event has different ownership or content.".into(),
                );
            }
            let source = frame
                .revisions
                .get(&requested.revision)
                .ok_or("The translation origin revision is unavailable.")?;
            if source.sequence.alphabet() == collagene_model::Alphabet::Protein {
                return Err("A translation origin must be a nucleotide reference.".into());
            }
            return Ok(source);
        }
        pending.extend(event.parents.iter().map(String::as_str));
    }
    Err("The declared translation source is not retained in this document's ancestry.".into())
}

#[cfg(test)]
mod tests {
    use super::*;
    use collagene_model::{Feature, Topology};
    fn source(m: &Molecule) -> EventBinding {
        EventBinding {
            document_id: "doc".into(),
            event_id: "event".into(),
            revision: m.revision(),
        }
    }
    #[test]
    fn circular_reverse_orf_highlights_all_parts_without_losing_biological_direction() {
        let m = Molecule::dna("circle", "ATGCCCAAAT", Topology::Circular).unwrap();
        let event = source(&m);
        let location = Location::Complement(Box::new(Location::Join(vec![
            Location::range(8, 10),
            Location::range(0, 4),
        ])));
        let directed = DirectedSelection::new(event.clone(), &m, location.clone()).unwrap();
        assert_eq!(directed.spans, vec![(8, 10), (0, 4)]);
        let regions = (event.revision.clone(), directed.spans.clone());
        assert_eq!(
            capture(&m, &event, None, None, Some(&regions), Some(&directed)).unwrap(),
            location
        );
        let changed = (event.revision.clone(), vec![(1, 4)]);
        assert_eq!(
            capture(
                &m,
                &event,
                None,
                Some((1, 4)),
                Some(&changed),
                Some(&directed)
            )
            .unwrap(),
            Location::range(1, 4)
        );
        let mut alias = event.clone();
        alias.document_id = "other".into();
        assert_eq!(
            capture(&m, &alias, None, None, Some(&regions), Some(&directed)).unwrap(),
            Location::Join(vec![Location::range(8, 10), Location::range(0, 4)])
        );
    }
    #[test]
    fn exact_feature_orientation_wins_and_carets_or_uncertain_locations_do_not_invent_bases() {
        let mut m = Molecule::dna("source", "ATGAAATAG", Topology::Linear).unwrap();
        m.features.push(Feature {
            kind: "gene".into(),
            location: Location::Complement(Box::new(Location::range(0, 9))),
            qualifiers: vec![],
            needs_review: false,
        });
        let event = source(&m);
        assert_eq!(
            capture(&m, &event, Some(0), Some((0, 9)), None, None).unwrap(),
            m.features[0].location
        );
        assert!(capture(&m, &event, None, Some((4, 4)), None, None).is_err());
        m.features[0].location = Location::Range {
            start: 0,
            end: 9,
            before: true,
            after: false,
        };
        assert!(capture(&m, &source(&m), Some(0), None, None, None).is_err());
    }
    #[test]
    fn translation_origin_keeps_the_historical_source_and_refuses_metadata_only_aliases() {
        use collagene_core::documents::{history::DocumentHistory, workspace::Workspace};
        let dna = Molecule::dna("template", "ATGGCTTAA", Topology::Linear).unwrap();
        let mut source = DocumentHistory::new("dna", dna.clone()).unwrap();
        let original = EventBinding {
            document_id: source.id.clone(),
            event_id: source.head_event().into(),
            revision: source.head_revision().into(),
        };
        let mut protein = Molecule::dna("protein", "ATGGCT", Topology::Linear).unwrap();
        protein.sequence =
            collagene_model::Sequence::new(collagene_model::Alphabet::Protein, "MA").unwrap();
        protein.metadata.insert(
            "translation:source_document".into(),
            original.document_id.clone(),
        );
        protein
            .metadata
            .insert("translation:source_event".into(), original.event_id.clone());
        protein.metadata.insert(
            "translation:source_revision".into(),
            original.revision.clone(),
        );
        let derived = source
            .derive_event(
                &original.event_id,
                "protein",
                protein.clone(),
                "Computed translation",
            )
            .unwrap();
        let protein_binding = EventBinding {
            document_id: derived.id.clone(),
            event_id: derived.head_event().into(),
            revision: derived.head_revision().into(),
        };
        let mut next = dna.clone();
        next.name = "later template name".into();
        source.commit(next, "Later edit").unwrap();
        let mut frame = Workspace::from_documents("translation", vec![dna.clone()], 0).unwrap();
        frame.documents.clear();
        frame.revisions.clear();
        frame.history.clear();
        frame.active_document = "protein".into();
        source.contribute(&mut frame).unwrap();
        derived.contribute(&mut frame).unwrap();
        assert_eq!(origin(&frame, &protein_binding, &original).unwrap(), &dna);
        let imported = DocumentHistory::new("imported", protein).unwrap();
        let imported_binding = EventBinding {
            document_id: imported.id.clone(),
            event_id: imported.head_event().into(),
            revision: imported.head_revision().into(),
        };
        imported.contribute(&mut frame).unwrap();
        assert!(origin(&frame, &imported_binding, &original).is_err());
        let mut alias = original;
        alias.document_id = "other".into();
        assert!(origin(&frame, &protein_binding, &alias).is_err());
    }
}
