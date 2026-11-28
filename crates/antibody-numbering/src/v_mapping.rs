//! Source-bound correspondence to explicitly selected species/strain V models.
//! Integer count scores are declared heuristics, not calibrated domain/chain/
//! species classification or full rearranged antibody numbering.
use super::*;
use crate::{profile, references};
use collagene_model::PROTEIN_RESIDUES;
use std::collections::BTreeSet;
pub const MAPPING_ENGINE: &str = "antibody-v-reference-mapping-5";
pub const MAX_MODEL_ALLELES: usize = 354;
data! {
    pub enum Locus {
        Heavy,
        Kappa,
        Lambda
    }
}
impl Locus {
    pub(crate) fn source_index(&self) -> usize {
        match self {
            Self::Heavy => 0,
            Self::Kappa => 1,
            Self::Lambda => 2,
        }
    }
}
data! {
    pub struct ReferenceSelection {
        pub species: references::Species,
        #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, Locus, 3>")]
        pub loci: Vec<Locus>
    }
}
impl ReferenceSelection {
    pub fn human_all() -> Self {
        Self {
            species: references::Species::Human,
            loci: vec![Locus::Heavy, Locus::Kappa, Locus::Lambda],
        }
    }
    pub fn validate(&self) -> Result<()> {
        if self.loci.is_empty()
            || self.loci.len() > 3
            || self
                .loci
                .iter()
                .enumerate()
                .any(|(i, locus)| self.loci[..i].contains(locus))
        {
            return Err(Error::Invalid(
                "explicit distinct nonempty reference loci required",
            ));
        }
        if self.species == references::Species::MouseC57Bl6 && self.loci != [Locus::Heavy] {
            return Err(Error::Invalid(
                "mouse C57BL/6 reference scope is heavy-chain IGH only",
            ));
        }
        Ok(())
    }
}
data! {
    pub struct Scoring {
        pub observed_reward: i16,
        pub frequency_bonus: i16,
        pub unseen_penalty: i16,
        pub gap_open: i16,
        pub gap_extend: i16,
        pub minimum_score: i32,
    }
}
impl Scoring {
    pub(crate) fn validate(&self) -> Result<()> {
        if !(1..=500).contains(&self.observed_reward)
            || !(0..=500).contains(&self.frequency_bonus)
            || !(-1000..=0).contains(&self.unseen_penalty)
            || !(-1000..=0).contains(&self.gap_open)
            || !(-1000..=0).contains(&self.gap_extend)
        {
            return Err(Error::Invalid("V-reference scoring weights"));
        }
        Ok(())
    }
}
data! {
    pub struct ColumnCounts {
        pub residues: [usize; 22],
        pub gaps: usize
    }
}
data! {
    pub struct Model {
        pub engine: String,
        pub corpus_revision: String,
        pub species: references::Species,
        pub locus: Locus,
        pub scoring: Scoring,
        pub training_alleles: Vec<String>,
        pub omitted_alleles: Vec<String>,
        pub distinct_proteins: usize,
        /// Model states retain corpus-declared column labels and their source policy;
        /// all-gap positions have
        /// no emitting state and stay explicit in omitted_positions instead.
        pub positions: Vec<Position>,
        pub omitted_positions: Vec<Position>,
        pub counts: Vec<ColumnCounts>,
        pub profile: profile::Profile,
    }
}
fn admitted(entry: &references::VReference, species: &references::Species) -> bool {
    entry.non_cysteine_anchors.is_empty()
        && (*species == references::Species::Human || entry.aligned_protein.as_bytes()[40] == b'W')
}
fn model_revision(model: &Model) -> Result<String> {
    hash(
        b"antibody-v-reference-model-5\0",
        &(profile::PROFILE_ENGINE, model),
        1 << 20,
    )
}
/// Reconstruct only explicitly selected models from the pinned corpus. No caller can
/// silently substitute a scoring matrix under an unchanged reference revision.
pub fn models(
    selection: &ReferenceSelection,
    scoring: &Scoring,
    cancel: &mut impl FnMut() -> bool,
) -> Result<Vec<Model>> {
    scoring.validate()?;
    selection.validate()?;
    if cancel() {
        return Err(Error::Cancelled);
    }
    let corpus = references::load(&selection.species)?;
    let corpus_revision = corpus.revision()?;
    let mut models = vec![];
    for locus in &selection.loci {
        let source = locus.source_index();
        let mut training_alleles = vec![];
        let mut omitted_alleles = vec![];
        let mut proteins = BTreeSet::new();
        for entry in corpus.v_entries().iter().filter(|e| e.source == source) {
            if cancel() {
                return Err(Error::Cancelled);
            }
            if admitted(entry, &selection.species) {
                training_alleles.push(entry.allele_ref.clone());
                proteins.insert(entry.aligned_protein.as_str());
            } else {
                omitted_alleles.push(entry.allele_ref.clone());
            }
        }
        let mut positions = vec![];
        let mut omitted_positions = vec![];
        let mut counts = vec![];
        let mut states = vec![];
        for column in 0..104 {
            let mut count = ColumnCounts {
                residues: [0; 22],
                gaps: 0,
            };
            for protein in &proteins {
                let b = protein.as_bytes()[column];
                if b == b'-' {
                    count.gaps += 1;
                } else {
                    let index = PROTEIN_RESIDUES
                        .iter()
                        .position(|&a| a == b)
                        .ok_or(Error::Invalid("reference protein alphabet"))?;
                    count.residues[index] += 1;
                }
            }
            let occupied = proteins.len() - count.gaps;
            let label = position(column as u16 + 1);
            if occupied == 0 {
                omitted_positions.push(label);
                continue;
            }
            let required = matches!(column, 22 | 103);
            let mut emissions = [None; 22];
            for (index, &n) in count.residues.iter().enumerate() {
                if required && PROTEIN_RESIDUES[index] != b'C' {
                    continue;
                }
                emissions[index] = Some(if n == 0 {
                    scoring.unseen_penalty
                } else {
                    scoring.observed_reward
                        + (usize::try_from(scoring.frequency_bonus).unwrap() * n / occupied) as i16
                });
            }
            positions.push(label);
            counts.push(count);
            states.push(profile::MatchState {
                emissions,
                required,
                allow_insert_after: column < 103,
            });
        }
        let admission = if selection.species == references::Species::Human {
            "Human admission: provider-functional V with canonical C23/C104; atypical reference facts remain in corpus"
        } else {
            "Mouse C57BL/6 admission: provider-functional IGH V with canonical C23/W41/C104; three atypical functional reference facts remain in corpus"
        };
        let profile = profile::Profile {
            version: profile::PROFILE_VERSION,
            name: format!(
                "{} {} V-prefix IMGT reference scores",
                if selection.species == references::Species::Human {
                    "Human"
                } else {
                    "Mouse C57BL/6"
                },
                corpus.sources()[source].locus
            ),
            provenance: format!(
                "{}; {}; {}; distinct protein count scores; not a calibrated classifier",
                MAPPING_ENGINE, corpus_revision, admission
            ),
            states,
            gap_open: scoring.gap_open,
            gap_extend: scoring.gap_extend,
            minimum_score: scoring.minimum_score,
        };
        profile.revision()?;
        models.push(Model {
            engine: MAPPING_ENGINE.into(),
            corpus_revision: corpus_revision.clone(),
            species: selection.species.clone(),
            locus: locus.clone(),
            scoring: scoring.clone(),
            training_alleles,
            omitted_alleles,
            distinct_proteins: proteins.len(),
            positions,
            omitted_positions,
            counts,
            profile,
        });
    }
    if cancel() {
        return Err(Error::Cancelled);
    }
    Ok(models)
}
impl Model {
    pub fn validate(&self) -> Result<()> {
        if self.engine != MAPPING_ENGINE
            || models(
                &ReferenceSelection {
                    species: self.species.clone(),
                    loci: vec![self.locus.clone()],
                },
                &self.scoring,
                &mut || false,
            )?
            .iter()
            .find(|m| m.locus == self.locus)
                != Some(self)
        {
            return Err(Error::Binding);
        }
        Ok(())
    }
}
data! {
    pub struct Limits {
        /// Aggregate grid cells across all explicitly selected models, with each matrix reused for
        /// lower/upper scoring bounds. Not a separate budget for each locus.
        pub max_cells: usize,
        pub max_endpoints_per_locus: usize,
        pub max_trace_states: usize,
    }
}
data! {
    pub enum ProjectionIssue {
        AmbiguousSource,
        TiedOptimalPaths
    }
}
data! {
    pub struct Correspondence {
        pub source_index: usize,
        pub residue: char,
        pub reference_position: Position,
    }
}
data! {
    pub struct Projection {
        /// Candidate reference correspondence, not biological scheme certification.
        pub matches: Vec<Correspondence>,
        pub deleted_positions: Vec<Position>,
        /// Query insertions deliberately receive no invented scheme labels.
        pub unassigned_residues: Vec<usize>,
    }
}
data! {
    pub struct EndpointProjection {
        pub end: usize,
        pub projection: Option<Projection>,
        pub projection_issue: Option<ProjectionIssue>,
        /// Positions shared by every optimal state path at this fixed end. None for
        /// ambiguous sources; tied paths may still share framework/loop boundaries.
        pub invariant_correspondences: Option<Vec<Correspondence>>,
    }
}
data! {
    pub struct Mapping {
        pub locus: Locus,
        pub model_revision: String,
        pub alignment: profile::EndpointReport,
        /// Exact original reference identities admitted to/omitted from this model;
        /// identical proteins may share scoring weight without losing allele provenance.
        #[serde(
            deserialize_with = "collagene_model::decode::bounded_vec::<_, String, MAX_MODEL_ALLELES>"
        )]
        pub training_alleles: Vec<String>,
        #[serde(
            deserialize_with = "collagene_model::decode::bounded_vec::<_, String, MAX_MODEL_ALLELES>"
        )]
        pub omitted_alleles: Vec<String>,
        pub projections: Vec<EndpointProjection>,
    }
}
data! {
    pub struct Report {
        pub engine: String,
        pub source_revision: String,
        pub corpus_revision: String,
        pub context_revision: String,
        pub selection: ReferenceSelection,
        pub scoring: Scoring,
        pub limits: Limits,
        pub cells: usize,
        pub invariant_trace_states: usize,
        pub mappings: Vec<Mapping>,
    }
}
/// Return every locus independently; uncalibrated scores of differently sized
/// profiles must not choose a winning chain type or imply that a domain exists.
pub fn map_query(
    source: &Molecule,
    selection: &ReferenceSelection,
    scoring: &Scoring,
    limits: &Limits,
    cancel: &mut impl FnMut() -> bool,
) -> Result<Report> {
    if cancel() {
        return Err(Error::Cancelled);
    }
    let source_revision = source_revision(source)?;
    scoring.validate()?;
    if limits.max_cells == 0
        || limits.max_cells > profile::MAX_CELLS
        || limits.max_endpoints_per_locus > profile::MAX_ENDPOINTS
        || limits.max_trace_states > profile::MAX_TRACE_STATES
    {
        return Err(Error::Limit("V-reference aggregate cells"));
    }
    let models = models(selection, scoring, cancel)?;
    let cells = models
        .iter()
        .try_fold(0_usize, |n, m| {
            (source.sequence.len() + 1)
                .checked_mul(m.profile.states.len() + 1)
                .and_then(|v| v.checked_add(n))
        })
        .filter(|&n| n <= limits.max_cells)
        .ok_or(Error::Limit("V-reference aggregate cells"))?;
    let corpus_revision = models[0].corpus_revision.clone();
    let context_revision = hash(
        b"antibody-v-reference-context-5\0",
        &(
            MAPPING_ENGINE,
            profile::PROFILE_ENGINE,
            profile::ENDPOINT_ENGINE,
            &source_revision,
            &corpus_revision,
            selection,
            scoring,
            limits,
            "all-optimal-path invariant correspondences with aggregate trace work; per-endpoint optima; exact endpoint counts and explicit truncation; distinct protein weighting; C23/C104 required; all-gap positions omitted; no classification; no inferred insertion labels; ambiguous/tied projections withheld",
        ),
        65536,
    )?;
    let mut remaining_trace = limits.max_trace_states;
    let mut mappings = vec![];
    for model in models {
        let model_revision = model_revision(&model)?;
        if model
            .training_alleles
            .len()
            .checked_add(model.omitted_alleles.len())
            .is_none_or(|n| n > MAX_MODEL_ALLELES)
        {
            return Err(Error::Limit("model admission identities"));
        }
        // Move exact selected-model facts before alignment; no second model or
        // result-sized identity clone is constructed for the report.
        let training_alleles = model.training_alleles;
        let omitted_alleles = model.omitted_alleles;
        let individual_cells = (source.sequence.len() + 1) * (model.profile.states.len() + 1);
        let alignment = profile::align_endpoints(
            source,
            &model.profile,
            &profile::EndpointLimits {
                max_cells: individual_cells,
                max_endpoints: limits.max_endpoints_per_locus,
                max_trace_states: remaining_trace,
            },
            cancel,
        )?;
        remaining_trace -= alignment.invariant_trace_states;
        let last = model.profile.states.len() - 1;
        let state = &model.profile.states[last];
        if model.positions[last] != position(104)
            || !state.required
            || state.allow_insert_after
            || state
                .emissions
                .iter()
                .enumerate()
                .any(|(i, e)| e.is_some() && PROTEIN_RESIDUES[i] != b'C')
        {
            return Err(Error::Invalid("V model invariant terminal anchor"));
        }
        let mut projections = vec![];
        for endpoint in &alignment.endpoints {
            let invariant_correspondences = endpoint.invariant_matches.as_ref().map(|matches| {
                matches
                    .iter()
                    .map(|m| Correspondence {
                        source_index: m.source_index,
                        residue: source.sequence.as_str().as_bytes()[m.source_index] as char,
                        reference_position: model.positions[m.model_index].clone(),
                    })
                    .collect::<Vec<_>>()
            });
            if alignment.ambiguous_residues == 0
                && !invariant_correspondences.as_ref().is_some_and(|matches| {
                    matches.iter().any(|m| {
                        m.reference_position == position(104)
                            && m.source_index + 1 == endpoint.end
                            && m.residue == 'C'
                    })
                })
            {
                return Err(Error::Binding);
            }
            let projection_issue = if alignment.ambiguous_residues > 0 {
                Some(ProjectionIssue::AmbiguousSource)
            } else if endpoint.upper.alternative_optimal_paths {
                Some(ProjectionIssue::TiedOptimalPaths)
            } else {
                None
            };
            let projection = if projection_issue.is_none() {
                let mut matches = vec![];
                let mut deleted_positions = vec![];
                let mut unassigned_residues = vec![];
                for step in &endpoint.upper.steps {
                    match *step {
                        profile::Step::Match {
                            source_index,
                            model_index,
                        } => matches.push(Correspondence {
                            source_index,
                            residue: source.sequence.as_str().as_bytes()[source_index] as char,
                            reference_position: model.positions[model_index].clone(),
                        }),
                        profile::Step::Delete { model_index } => {
                            deleted_positions.push(model.positions[model_index].clone())
                        }
                        profile::Step::Insert { source_index, .. } => {
                            unassigned_residues.push(source_index)
                        }
                    }
                }
                Some(Projection {
                    matches,
                    deleted_positions,
                    unassigned_residues,
                })
            } else {
                None
            };
            projections.push(EndpointProjection {
                end: endpoint.end,
                projection,
                projection_issue,
                invariant_correspondences,
            });
        }
        mappings.push(Mapping {
            locus: model.locus,
            model_revision,
            training_alleles,
            omitted_alleles,
            alignment,
            projections,
        });
    }
    if cancel() {
        return Err(Error::Cancelled);
    }
    Ok(Report {
        engine: MAPPING_ENGINE.into(),
        source_revision,
        corpus_revision,
        context_revision,
        selection: selection.clone(),
        scoring: scoring.clone(),
        limits: limits.clone(),
        cells,
        invariant_trace_states: limits.max_trace_states - remaining_trace,
        mappings,
    })
}
impl Mapping {
    /// Cheap current-corpus admission facts only; no model scoring/alignment work.
    pub(crate) fn validate_admission(&self, species: &references::Species) -> Result<()> {
        if self
            .training_alleles
            .len()
            .checked_add(self.omitted_alleles.len())
            .is_none_or(|n| n > MAX_MODEL_ALLELES)
        {
            return Err(Error::Limit("model admission identities"));
        }
        let corpus = references::load(species)?;
        let entries = || {
            corpus
                .v_entries()
                .iter()
                .filter(|entry| entry.source == self.locus.source_index())
        };
        if !self
            .training_alleles
            .iter()
            .map(String::as_str)
            .eq(entries()
                .filter(|entry| admitted(entry, species))
                .map(|entry| entry.allele_ref.as_str()))
            || !self.omitted_alleles.iter().map(String::as_str).eq(entries()
                .filter(|entry| !admitted(entry, species))
                .map(|entry| entry.allele_ref.as_str()))
        {
            return Err(Error::Binding);
        }
        Ok(())
    }
}
impl Report {
    pub fn validate_for(
        &self,
        source: &Molecule,
        selection: &ReferenceSelection,
        scoring: &Scoring,
        limits: &Limits,
    ) -> Result<()> {
        for mapping in &self.mappings {
            mapping.validate_admission(&selection.species)?;
        }
        if self.engine != MAPPING_ENGINE
            || self.selection != *selection
            || self.mappings.len() != selection.loci.len()
            || self.mappings.iter().any(|m| {
                m.projections.len() > profile::MAX_ENDPOINTS
                    || m.projections.iter().any(|e| {
                        e.projection.as_ref().is_some_and(|p| {
                            p.matches.len() > MAX_SOURCE_RESIDUES
                                || p.deleted_positions.len() > 104
                                || p.unassigned_residues.len() > MAX_SOURCE_RESIDUES
                        })
                    })
            })
            || *self != map_query(source, selection, scoring, limits, &mut || false)?
        {
            return Err(Error::Binding);
        }
        Ok(())
    }
}
