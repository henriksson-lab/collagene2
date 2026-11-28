//! Exact event/options guards for ephemeral conditional antibody analysis.
use collagene_core::{
    antibody_numbering::{self as n, domain_candidates as d, domain_tracks},
    material::lineage::EventBinding,
};
use collagene_model::Molecule;
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Capture {
    pub epoch: u64,
    pub source: EventBinding,
    pub source_length: usize,
    pub options: d::Options,
    pub limits: d::Limits,
    pub selection: domain_tracks::Selection,
    pub context_revision: String,
}
pub fn reference_label(references: &n::v_mapping::ReferenceSelection) -> String {
    let species = match references.species {
        n::references::Species::Human => "Human",
        n::references::Species::MouseC57Bl6 => "Mouse C57BL/6",
    };
    let loci = references
        .loci
        .iter()
        .map(|l| match l {
            n::v_mapping::Locus::Heavy => "IGH",
            n::v_mapping::Locus::Kappa => "IGK",
            n::v_mapping::Locus::Lambda => "IGL",
        })
        .collect::<Vec<_>>()
        .join(", ");
    format!(
        "{species} reference scope: {loci}. This is an explicit reference choice, not a species or chain prediction."
    )
}
pub fn validate_display_scope(
    references: &n::v_mapping::ReferenceSelection,
    selection: &domain_tracks::Selection,
) -> Result<(), String> {
    selection.validate_for_references(references).map_err(|e| {
        format!("{e}. Choose a supported explicit reference, numbering and region combination; existing choices are not converted.")
    })
}
pub fn numbering_scheme(value: &str) -> Result<domain_tracks::NumberingScheme, String> {
    match value {
        "Imgt" => Ok(domain_tracks::NumberingScheme::Imgt),
        "Kabat1991" => Ok(domain_tracks::NumberingScheme::Kabat1991),
        "Chothia1997" => Ok(domain_tracks::NumberingScheme::Chothia1997),
        _ => Err("Unsupported antibody numbering scheme.".into()),
    }
}
pub fn region_definition(value: &str) -> Result<domain_tracks::RegionDefinition, String> {
    match value {
        "Imgt" => Ok(domain_tracks::RegionDefinition::Imgt),
        "Kabat1991" => Ok(domain_tracks::RegionDefinition::Kabat1991),
        "ChothiaCanonical" => Ok(domain_tracks::RegionDefinition::ChothiaCanonical),
        "North2011" => Ok(domain_tracks::RegionDefinition::North2011),
        _ => Err("Unsupported antibody region definition.".into()),
    }
}
pub fn preset(
    references: n::v_mapping::ReferenceSelection,
    minimum_score: i32,
    max_fr4_mismatches: usize,
    max_cdr3_residues: usize,
    max_candidates: usize,
) -> (d::Options, d::Limits) {
    (
        d::Options {
            references,
            v_scoring: n::v_mapping::Scoring {
                observed_reward: 5,
                frequency_bonus: 8,
                unseen_penalty: -6,
                gap_open: -12,
                gap_extend: -2,
                minimum_score,
            },
            min_cdr3_residues: 2,
            max_cdr3_residues,
            max_fr4_mismatches,
        },
        d::Limits {
            v_mapping: n::v_mapping::Limits {
                max_cells: n::profile::MAX_CELLS,
                max_endpoints_per_locus: n::profile::MAX_ENDPOINTS,
                max_trace_states: n::profile::MAX_TRACE_STATES,
            },
            max_scan_windows: d::MAX_SCAN_WINDOWS,
            max_candidates,
        },
    )
}
impl Capture {
    pub fn prepare(
        epoch: u64,
        source: EventBinding,
        m: &Molecule,
        options: d::Options,
        limits: d::Limits,
        selection: domain_tracks::Selection,
    ) -> Result<Self, String> {
        let context_revision = domain_tracks::request_revision(m, &options, &limits, &selection)
            .map_err(|e| e.to_string())?;
        if n::source_revision(m).map_err(|e| e.to_string())? != source.revision {
            return Err("Protein source revision differs.".into());
        }
        Ok(Self {
            epoch,
            source,
            source_length: m.sequence.len(),
            options,
            limits,
            selection,
            context_revision,
        })
    }
    pub fn guard(
        &self,
        epoch: u64,
        source: &EventBinding,
        options: &d::Options,
        limits: &d::Limits,
        selection: &domain_tracks::Selection,
    ) -> Result<(), String> {
        if self.epoch != epoch
            || self.source != *source
            || self.options != *options
            || self.limits != *limits
            || self.selection != *selection
        {
            return Err("Source, analysis settings or workspace changed; analyze again.".into());
        }
        Ok(())
    }
    pub fn validate_report(&self, r: &domain_tracks::Report) -> Result<(), String> {
        r.validate_binding(&self.context_revision)
            .map_err(|e| e.to_string())?;
        if r.source_revision != self.source.revision
            || r.evidence.options != self.options
            || r.evidence.limits != self.limits
            || r.selection != self.selection
            || r.tracks.iter().any(|t| t.source_end > self.source_length)
        {
            return Err("Antibody reply belongs to another source or analysis scope.".into());
        }
        Ok(())
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    use collagene_core::documents::workspace::Workspace;
    use collagene_model::{Alphabet, Sequence, Topology};
    #[test]
    fn identical_documents_different_events_epochs_and_settings_do_not_retarget_replies() {
        let m = Molecule {
            name: "Protein".into(),
            sequence: Sequence::new(Alphabet::Protein, "NATMW").unwrap(),
            topology: Topology::Linear,
            features: vec![],
            metadata: Default::default(),
            original_genbank: None,
        };
        let w = Workspace::from_documents("Sources", vec![m.clone(), m.clone()], 0).unwrap();
        let binding = |i: usize| {
            let d = &w.documents[i];
            EventBinding {
                document_id: d.id.clone(),
                event_id: d.head_event.clone().unwrap(),
                revision: d.head.clone(),
            }
        };
        let (o, l) = preset(
            n::v_mapping::ReferenceSelection::human_all(),
            300,
            1,
            60,
            32,
        );
        let a = binding(0);
        let b = binding(1);
        assert_eq!(a.revision, b.revision);
        let c = Capture::prepare(3, a.clone(), &m, o.clone(), l.clone(), test_selection()).unwrap();
        let r = domain_tracks::analyze(&m, &o, &l, &test_selection(), &mut || false).unwrap();
        c.validate_report(&r).unwrap();
        let mut selection = test_selection();
        selection.numbering_scheme = domain_tracks::NumberingScheme::Kabat1991;
        assert!(c.guard(3, &a, &o, &l, &selection).is_err());
        let mut definition = test_selection();
        definition.region_definition = domain_tracks::RegionDefinition::ChothiaCanonical;
        assert!(c.guard(3, &a, &o, &l, &definition).is_err());
        c.guard(3, &a, &o, &l, &test_selection()).unwrap();
        let mut reference_change = o.clone();
        reference_change.references = n::v_mapping::ReferenceSelection {
            species: n::references::Species::MouseC57Bl6,
            loci: vec![n::v_mapping::Locus::Heavy],
        };
        assert!(
            c.guard(3, &a, &reference_change, &l, &test_selection())
                .is_err()
        );
        let mut locus_change = o.clone();
        locus_change.references.loci = vec![n::v_mapping::Locus::Heavy];
        assert!(
            c.guard(3, &a, &locus_change, &l, &test_selection())
                .is_err()
        );
        assert!(c.guard(3, &b, &o, &l, &test_selection()).is_err());
        assert!(c.guard(4, &a, &o, &l, &test_selection()).is_err());
        let mut changed = l;
        changed.max_candidates = 0;
        assert!(c.guard(3, &a, &o, &changed, &test_selection()).is_err());
        let mut changed = o;
        changed.max_fr4_mismatches = 0;
        assert!(
            c.guard(3, &a, &changed, &c.limits, &test_selection())
                .is_err()
        );
    }
    #[test]
    fn mouse_display_does_not_silently_convert_unsupported_human_rules() {
        let mouse = n::v_mapping::ReferenceSelection {
            species: n::references::Species::MouseC57Bl6,
            loci: vec![n::v_mapping::Locus::Heavy],
        };
        validate_display_scope(&mouse, &test_selection()).unwrap();
        let mut alternate = test_selection();
        alternate.numbering_scheme = domain_tracks::NumberingScheme::Kabat1991;
        assert!(validate_display_scope(&mouse, &alternate).is_err());
        alternate = test_selection();
        alternate.region_definition = domain_tracks::RegionDefinition::North2011;
        assert!(validate_display_scope(&mouse, &alternate).is_err());
        assert!(numbering_scheme("unknown").is_err());
        assert!(region_definition("unknown").is_err());
    }
}

#[cfg(test)]
fn test_selection() -> domain_tracks::Selection {
    domain_tracks::Selection {
        numbering_scheme: collagene_core::antibody_numbering::domain_tracks::NumberingScheme::Imgt,
        region_definition:
            collagene_core::antibody_numbering::domain_tracks::RegionDefinition::Imgt,
    }
}
