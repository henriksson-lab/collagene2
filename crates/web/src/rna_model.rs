//! RNA form edits and exact-source response guards; no alternate saved encoding.
use collagene_core::{Molecule, Topology, rna_chemistry as rna, run_lineage::EventBinding};
pub fn blank(revision: String, topology: Topology) -> rna::Declaration {
    rna::Declaration {
        version: rna::DECLARATION_VERSION,
        source_revision: revision,
        coverage: rna::Coverage::Partial,
        ends: if topology == Topology::Circular {
            rna::Ends::Circular
        } else {
            rna::Ends::Linear {
                five: rna::FivePrime::Unknown {
                    label: "Unspecified 5' terminus".into(),
                },
                three: rna::ThreePrime::Unknown {
                    label: "Unspecified 3' terminus".into(),
                },
            }
        },
        sites: vec![],
        links: vec![],
        annotation_reviews: vec![],
        provenance: rna::Provenance {
            statement: String::new(),
            reference: None,
        },
    }
}
pub fn position(text: &str, len: usize) -> Result<usize, String> {
    text.trim()
        .parse::<usize>()
        .ok()
        .filter(|&i| i > 0 && i <= len)
        .map(|i| i - 1)
        .ok_or_else(|| format!("Enter a position from 1 to {len}."))
}
pub fn upsert_site(d: &mut rna::Declaration, site: rna::Site) -> Result<(), String> {
    if let Some(s) = d.sites.iter_mut().find(|s| s.index == site.index) {
        *s = site;
    } else {
        if d.sites.len() >= rna::MAX_SITES {
            return Err("Residue declaration limit reached.".into());
        }
        d.sites.push(site);
    }
    Ok(())
}
pub fn remove_site(d: &mut rna::Declaration, entry: usize) {
    if entry >= d.sites.len() {
        return;
    }
    d.sites.remove(entry);
    d.annotation_reviews.retain(|r| r.site_entry != entry);
    for r in &mut d.annotation_reviews {
        if r.site_entry > entry {
            r.site_entry -= 1;
        }
    }
}
pub fn accept_annotation(
    d: &mut rna::Declaration,
    m: &Molecule,
    feature_index: usize,
) -> Result<(), String> {
    if d.annotation_reviews
        .iter()
        .any(|r| r.feature_index == feature_index)
    {
        return Err("This feature already has an explicit review.".into());
    }
    let c = rna::annotation_candidates(m)
        .map_err(|e| e.to_string())?
        .into_iter()
        .find(|c| c.feature_index == feature_index)
        .ok_or("Modification feature is unavailable.")?;
    let suggested = c.suggested_site.ok_or_else(|| {
        c.unsupported_reason
            .unwrap_or("Unsupported annotation.".into())
    })?;
    let code = c.raw_code.as_ref().unwrap().trim().to_ascii_lowercase();
    let profile = rna::annotation_profiles()
        .into_iter()
        .find(|p| p.code == code)
        .ok_or("Unsupported annotation code.")?;
    let mut trial = d.clone();
    let entry = if let Some(i) = trial.sites.iter().position(|s| s.index == suggested.index) {
        if let Some(n) = profile.nucleoside {
            trial.sites[i].nucleoside = n;
        }
        if let Some(s) = profile.sugar {
            trial.sites[i].sugar = s;
        }
        i
    } else {
        let i = trial.sites.len();
        upsert_site(&mut trial, suggested)?;
        i
    };
    trial.annotation_reviews.push(rna::FeatureReview {
        feature_index,
        expected_feature_revision: c.feature_revision,
        site_entry: entry,
    });
    rna::declaration_revision(m, &trial).map_err(|e| e.to_string())?;
    *d = trial;
    Ok(())
}
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Capture {
    pub epoch: u64,
    pub source: EventBinding,
    pub declaration: rna::Declaration,
    pub context_revision: String,
}
impl Capture {
    pub fn prepare(
        epoch: u64,
        source: EventBinding,
        m: &Molecule,
        declaration: rna::Declaration,
    ) -> Result<Self, String> {
        let context_revision = rna::request_revision(m, &declaration).map_err(|e| e.to_string())?;
        if declaration.source_revision != source.revision {
            return Err("RNA source event revision differs.".into());
        }
        Ok(Self {
            epoch,
            source,
            declaration,
            context_revision,
        })
    }
    pub fn guard(
        &self,
        epoch: u64,
        source: &EventBinding,
        d: &rna::Declaration,
    ) -> Result<(), String> {
        if self.epoch != epoch || self.source != *source || self.declaration != *d {
            return Err("RNA source, declaration or workspace changed; calculate again.".into());
        }
        Ok(())
    }
    pub fn validate_report(&self, r: &rna::Report) -> Result<(), String> {
        r.validate_binding(&self.context_revision)
            .map_err(|e| e.to_string())?;
        if r.declaration != self.declaration || r.source_revision != self.source.revision {
            return Err("RNA report belongs to another declaration.".into());
        }
        Ok(())
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    use collagene_core::{Alphabet, Feature, Location, Sequence};
    fn source() -> Molecule {
        Molecule {
            name: "RNA".into(),
            sequence: Sequence::new(Alphabet::Rna, "CC").unwrap(),
            topology: Topology::Linear,
            features: vec![],
            metadata: Default::default(),
            original_genbank: None,
        }
    }
    #[test]
    fn coordinate_and_row_edits_preserve_exact_review_entry_identity() {
        assert_eq!(position("1", 2).unwrap(), 0);
        assert_eq!(position("2", 2).unwrap(), 1);
        for s in ["0", "-1", "3", "1.5", ""] {
            assert!(position(s, 2).is_err());
        }
        let m = source();
        let mut d = blank(m.revision(), m.topology);
        d.provenance.statement = "Explicit model".into();
        let provenance = d.provenance.clone();
        for i in 0..2 {
            upsert_site(
                &mut d,
                rna::Site {
                    index: i,
                    nucleoside: rna::Nucleoside::Canonical,
                    sugar: rna::Sugar::Ribose,
                    provenance: provenance.clone(),
                },
            )
            .unwrap();
        }
        d.annotation_reviews = vec![
            rna::FeatureReview {
                feature_index: 7,
                expected_feature_revision: "a".repeat(64),
                site_entry: 1,
            },
            rna::FeatureReview {
                feature_index: 8,
                expected_feature_revision: "b".repeat(64),
                site_entry: 0,
            },
        ];
        remove_site(&mut d, 0);
        assert_eq!(d.sites[0].index, 1);
        assert_eq!(d.annotation_reviews.len(), 1);
        assert_eq!(d.annotation_reviews[0].feature_index, 7);
        assert_eq!(d.annotation_reviews[0].site_entry, 0);
    }
    #[test]
    fn combined_annotation_components_merge_and_stale_reply_guards_are_exact() {
        let mut m = source();
        for code in ["m5c", "cm"] {
            m.features.push(Feature {
                kind: "modified_base".into(),
                location: Location::range(0, 1),
                qualifiers: vec![("mod_base".into(), Some(code.into()))],
                needs_review: false,
            });
        }
        let mut d = blank(m.revision(), m.topology);
        d.provenance.statement = "Explicit hypothesis".into();
        accept_annotation(&mut d, &m, 0).unwrap();
        accept_annotation(&mut d, &m, 1).unwrap();
        assert_eq!(d.sites.len(), 1);
        assert_eq!(d.sites[0].nucleoside, rna::Nucleoside::FiveMethylcytidine);
        assert_eq!(d.sites[0].sugar, rna::Sugar::TwoPrimeOMethyl);
        assert!(accept_annotation(&mut d, &m, 0).is_err());
        let source = EventBinding {
            document_id: "doc-a".into(),
            event_id: "event-a".into(),
            revision: m.revision(),
        };
        let c = Capture::prepare(5, source.clone(), &m, d.clone()).unwrap();
        let r = rna::analyze(&m, &d).unwrap();
        c.validate_report(&r).unwrap();
        c.guard(5, &source, &d).unwrap();
        assert!(c.guard(6, &source, &d).is_err());
        let mut other = source.clone();
        other.event_id = "event-b".into();
        assert!(c.guard(5, &other, &d).is_err());
        other = source.clone();
        other.document_id = "doc-b".into();
        assert!(c.guard(5, &other, &d).is_err());
        d.coverage = rna::Coverage::ClosedModel;
        assert!(c.guard(5, &source, &d).is_err());
    }
}
