//! Exact event guards around reference-only collection similarity results.
use collagene_core::{run_lineage::EventBinding, sequence_similarity as native};
use collagene_model::Molecule;
#[derive(Clone, Debug, PartialEq, Eq, serde::Serialize)]
pub struct Pin {
    pub index: usize,
    pub binding: EventBinding,
    pub length: usize,
    pub alphabet: collagene_model::Alphabet,
    pub topology: collagene_model::Topology,
}
#[derive(Clone, Debug, PartialEq, Eq, serde::Serialize)]
pub struct Capture {
    pub epoch: u64,
    pub query: EventBinding,
    pub query_length: usize,
    pub query_alphabet: collagene_model::Alphabet,
    pub subjects: Vec<Pin>,
    pub options: native::Options,
    pub context_revision: String,
}
pub struct Prepared {
    pub capture: Capture,
    pub query: native::Source,
    pub subjects: Vec<native::Source>,
}
pub fn prepare(
    epoch: u64,
    query: EventBinding,
    molecule: &Molecule,
    subjects: &[(usize, EventBinding, &Molecule)],
    options: native::Options,
) -> Result<Prepared, String> {
    let refs = subjects
        .iter()
        .map(|(_, b, m)| (b.document_id.as_str(), *m))
        .collect::<Vec<_>>();
    native::preflight_borrowed((&query.document_id, molecule), &refs, &options)
        .map_err(|e| e.to_string())?;
    if molecule.revision() != query.revision
        || subjects.iter().any(|(_, b, m)| m.revision() != b.revision)
    {
        return Err("Similarity source events changed.".into());
    }
    let pins = subjects
        .iter()
        .map(|(index, b, m)| Pin {
            index: *index,
            binding: b.clone(),
            length: m.sequence.len(),
            alphabet: m.sequence.alphabet(),
            topology: m.topology,
        })
        .collect();
    let q = native::Source {
        key: query.document_id.clone(),
        molecule: molecule.clone(),
    };
    let sources = subjects
        .iter()
        .map(|(_, b, m)| native::Source {
            key: b.document_id.clone(),
            molecule: (*m).clone(),
        })
        .collect::<Vec<_>>();
    let context_revision =
        native::context_revision(&q, &sources, &options).map_err(|e| e.to_string())?;
    Ok(Prepared {
        capture: Capture {
            epoch,
            query,
            query_length: molecule.sequence.len(),
            query_alphabet: molecule.sequence.alphabet(),
            subjects: pins,
            options,
            context_revision,
        },
        query: q,
        subjects: sources,
    })
}
impl Capture {
    pub fn validate_current(
        &self,
        epoch: u64,
        query: &EventBinding,
        current: &[(usize, EventBinding)],
    ) -> Result<(), String> {
        if self.epoch != epoch
            || self.query != *query
            || current.len() != self.subjects.len()
            || current
                .iter()
                .zip(&self.subjects)
                .any(|((index, b), p)| *index != p.index || *b != p.binding)
        {
            return Err("Similarity sources or workspace changed. Run the search again.".into());
        }
        Ok(())
    }
    pub fn validate_report(&self, report: &native::Report) -> Result<(), String> {
        report
            .validate_binding(&self.context_revision)
            .map_err(|e| e.to_string())?;
        if report.query_key != self.query.document_id
            || report.query_revision != self.query.revision
            || report.query_length != self.query_length
            || report.query_alphabet != self.query_alphabet
            || report.options != self.options
            || report.subjects.len() != self.subjects.len()
            || report.subjects.iter().zip(&self.subjects).any(|(row, p)| {
                row.key != p.binding.document_id
                    || row.revision != p.binding.revision
                    || row.length != p.length
                    || row.alphabet != p.alphabet
                    || row.topology != p.topology
            })
        {
            return Err("Similarity reply does not match its exact source collection.".into());
        }
        Ok(())
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    use collagene_core::workspace::Workspace;
    use collagene_model::Topology;
    fn binding(w: &Workspace, i: usize) -> EventBinding {
        let d = &w.documents[i];
        EventBinding {
            document_id: d.id.clone(),
            event_id: d.head_event.clone().unwrap(),
            revision: d.head.clone(),
        }
    }
    #[test]
    fn equal_content_event_and_epoch_changes_refuse_navigation_and_reply_retargeting() {
        let molecule = Molecule::dna("Same content", "ACGT", Topology::Circular).unwrap();
        let w = Workspace::from_documents("Search", vec![molecule.clone(), molecule], 0).unwrap();
        let q = binding(&w, 0);
        let s = binding(&w, 1);
        let prepared = prepare(
            7,
            q.clone(),
            w.head(&q.document_id).unwrap(),
            &[(1, s.clone(), w.head(&s.document_id).unwrap())],
            native::Options::default(),
        )
        .unwrap();
        let report = native::search(
            &prepared.query,
            &prepared.subjects,
            &prepared.capture.options,
            &mut || false,
        )
        .unwrap();
        prepared.capture.validate_report(&report).unwrap();
        prepared
            .capture
            .validate_current(7, &q, &[(1, s.clone())])
            .unwrap();
        let mut same = s.clone();
        same.event_id = "another-retained-event".into();
        assert!(
            prepared
                .capture
                .validate_current(7, &q, &[(1, same)])
                .is_err()
        );
        assert!(
            prepared
                .capture
                .validate_current(8, &q, &[(1, s.clone())])
                .is_err()
        );
        assert!(
            prepared
                .capture
                .validate_current(7, &s, &[(1, s.clone())])
                .is_err()
        );
        let mut altered = report;
        altered.subjects[0].key = "other-document".into();
        assert!(prepared.capture.validate_report(&altered).is_err());
    }
    #[test]
    fn borrowed_preflight_and_actual_revision_checks_precede_owned_job_construction() {
        let m = Molecule::dna("m", "ACG", Topology::Linear).unwrap();
        let q = EventBinding {
            document_id: "query".into(),
            event_id: "event-1".into(),
            revision: m.revision(),
        };
        let mut stale = q.clone();
        stale.revision = "a".repeat(64);
        assert!(
            prepare(
                0,
                stale,
                &m,
                &[(0, q.clone(), &m)],
                native::Options::default()
            )
            .is_err()
        );
        let p = prepare(0, q.clone(), &m, &[(0, q, &m)], native::Options::default()).unwrap();
        assert_eq!(p.capture.subjects.len(), 1);
    }
}
