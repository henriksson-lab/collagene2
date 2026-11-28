//! Ephemeral source/event guards around the single native options/report format.
use collagene_core::{protein_liabilities as native, run_lineage::EventBinding};
use collagene_model::Molecule;
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ProvidedContext {
    pub source: EventBinding,
    pub context: native::ProvidedCdrs,
}
impl ProvidedContext {
    pub fn options(
        &self,
        source: &EventBinding,
        max_hits: usize,
    ) -> Result<native::Options, String> {
        if self.source != *source || self.context.source_revision != source.revision {
            return Err(
                "Supplied CDRs belong to another source version. Reset the context explicitly."
                    .into(),
            );
        }
        Ok(native::Options {
            max_hits,
            provided_cdrs: Some(self.context.clone()),
        })
    }
}
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Capture {
    pub epoch: u64,
    pub source: EventBinding,
    pub options: native::Options,
    pub source_length: usize,
    pub context_revision: String,
}
impl Capture {
    pub fn prepare(
        epoch: u64,
        source: EventBinding,
        m: &Molecule,
        options: native::Options,
    ) -> Result<Self, String> {
        let context_revision = native::request_revision(m, &options).map_err(|e| e.to_string())?;
        if native::source_revision(m).map_err(|e| e.to_string())? != source.revision {
            return Err("Protein source revision differs.".into());
        }
        Ok(Self {
            epoch,
            source,
            options,
            source_length: m.sequence.len(),
            context_revision,
        })
    }
    pub fn guard(
        &self,
        epoch: u64,
        source: &EventBinding,
        options: &native::Options,
    ) -> Result<(), String> {
        if self.epoch != epoch || self.source != *source || self.options != *options {
            return Err("Protein source, options or workspace changed; screen again.".into());
        }
        Ok(())
    }
    pub fn validate_report(&self, r: &native::Report) -> Result<(), String> {
        r.validate_binding(&self.context_revision)
            .map_err(|e| e.to_string())?;
        if r.source_revision != self.source.revision
            || r.source_length != self.source_length
            || r.options != self.options
        {
            return Err("Liability reply belongs to another source/options scope.".into());
        }
        Ok(())
    }
}
pub fn interval(start: &str, end: &str, length: usize) -> Result<(usize, usize), String> {
    let a = start.trim().parse::<usize>().ok();
    let b = end.trim().parse::<usize>().ok();
    match (a, b) {
        (Some(a), Some(b)) if a > 0 && a <= b && b <= length => Ok((a - 1, b)),
        _ => Err(format!(
            "Enter a one-based inclusive range within 1–{length}."
        )),
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    use collagene_core::workspace::Workspace;
    use collagene_model::{Alphabet, Sequence, Topology};
    fn molecule() -> Molecule {
        Molecule {
            name: "Protein".into(),
            sequence: Sequence::new(Alphabet::Protein, "NGNAT").unwrap(),
            topology: Topology::Linear,
            features: vec![],
            metadata: Default::default(),
            original_genbank: None,
        }
    }
    fn binding(w: &Workspace, i: usize) -> EventBinding {
        let d = &w.documents[i];
        EventBinding {
            document_id: d.id.clone(),
            event_id: d.head_event.clone().unwrap(),
            revision: d.head.clone(),
        }
    }
    #[test]
    fn equal_content_foreign_events_and_epochs_do_not_retarget_context_or_replies() {
        let m = molecule();
        let w =
            Workspace::from_documents("Protein sources", vec![m.clone(), m.clone()], 0).unwrap();
        let a = binding(&w, 0);
        let b = binding(&w, 1);
        assert_eq!(a.revision, b.revision);
        let context = ProvidedContext {
            source: a.clone(),
            context: native::ProvidedCdrs {
                source_revision: a.revision.clone(),
                definition_label: "Provided CDR only".into(),
                provenance: "Authored fixture".into(),
                regions: vec![native::CdrRegion {
                    domain_id: "domain-a".into(),
                    kind: native::CdrKind::Cdr1,
                    start: 0,
                    end: 1,
                }],
            },
        };
        assert!(context.options(&b, 100).is_err());
        let o = context.options(&a, 100).unwrap();
        let capture = Capture::prepare(3, a.clone(), &m, o.clone()).unwrap();
        let r = native::screen(&m, &o, &mut || false).unwrap();
        capture.validate_report(&r).unwrap();
        assert!(capture.guard(3, &b, &o).is_err());
        assert!(capture.guard(4, &a, &o).is_err());
        let mut changed = o;
        changed.max_hits = 0;
        assert!(capture.guard(3, &a, &changed).is_err());
    }
    #[test]
    fn inclusive_form_ranges_preserve_single_residues_and_refuse_invalid_bounds() {
        assert_eq!(interval("1", "1", 5).unwrap(), (0, 1));
        assert_eq!(interval("2", "5", 5).unwrap(), (1, 5));
        for (a, b) in [
            ("0", "1"),
            ("3", "2"),
            ("1", "6"),
            ("-1", "2"),
            ("1.5", "3"),
        ] {
            assert!(interval(a, b, 5).is_err());
        }
    }
}
