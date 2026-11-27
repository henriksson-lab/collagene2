//! Source-bound correspondence to human IMGT V-prefix reference models.
//! Integer count scores are declared heuristics, not calibrated domain/chain/
//! species classification or full rearranged antibody numbering.
use super::*;
use crate::{profile, references};
use collagene_model::PROTEIN_RESIDUES;
use std::collections::BTreeSet;
pub const MAPPING_ENGINE: &str = "antibody-v-reference-mapping-1";
data! {pub enum Locus {Heavy,Kappa,Lambda}}
data! {pub struct Scoring {
    pub observed_reward:i16,pub frequency_bonus:i16,pub unseen_penalty:i16,
    pub gap_open:i16,pub gap_extend:i16,pub minimum_score:i32,
}}
impl Scoring {
    fn validate(&self) -> Result<()> {
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
data! {pub struct ColumnCounts {pub residues:[usize;22],pub gaps:usize}}
data! {pub struct Model {
    pub engine:String,pub corpus_revision:String,pub locus:Locus,pub scoring:Scoring,
    pub training_alleles:Vec<String>,pub omitted_alleles:Vec<String>,pub distinct_proteins:usize,
    /// Model states retain the explicit provider labels; all-gap positions have
    /// no emitting state and stay explicit in omitted_positions instead.
    pub positions:Vec<Position>,pub omitted_positions:Vec<Position>,
    pub counts:Vec<ColumnCounts>,pub profile:profile::Profile,
}}
fn model_revision(model: &Model) -> Result<String> {
    hash(
        b"antibody-v-reference-model-1\0",
        &(profile::PROFILE_ENGINE, model),
        1 << 20,
    )
}
/// Reconstruct all three models from the current embedded corpus. No caller can
/// silently substitute a scoring matrix under an unchanged reference revision.
pub fn models(scoring: &Scoring, cancel: &mut impl FnMut() -> bool) -> Result<Vec<Model>> {
    scoring.validate()?;
    if cancel() {
        return Err(Error::Cancelled);
    }
    let corpus = references::human_v_references()?;
    let corpus_revision = corpus.revision()?;
    let mut models = vec![];
    for (source, locus) in [Locus::Heavy, Locus::Kappa, Locus::Lambda]
        .into_iter()
        .enumerate()
    {
        let mut training_alleles = vec![];
        let mut omitted_alleles = vec![];
        let mut proteins = BTreeSet::new();
        for entry in corpus.entries().iter().filter(|e| e.source == source) {
            if cancel() {
                return Err(Error::Cancelled);
            }
            if entry.non_cysteine_anchors.is_empty() {
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
        let profile = profile::Profile {
            version: profile::PROFILE_VERSION,
            name: format!(
                "Human {} V-prefix IMGT reference scores",
                corpus.sources()[source].locus
            ),
            provenance: format!(
                "{}; {}; distinct protein count scores; not a calibrated classifier",
                MAPPING_ENGINE, corpus_revision
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
            locus,
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
            || models(&self.scoring, &mut || false)?
                .iter()
                .find(|m| m.locus == self.locus)
                != Some(self)
        {
            return Err(Error::Binding);
        }
        Ok(())
    }
}
data! {pub struct Limits {
    /// Aggregate grid cells across all three models, with each matrix reused for
    /// lower/upper scoring bounds. Not a separate budget for each locus.
    pub max_cells:usize,
}}
data! {pub enum ProjectionIssue {NoPath,BelowScoreThreshold,AmbiguousSource,TiedOptimalPaths}}
data! {pub struct Correspondence {
    pub source_index:usize,pub residue:char,pub reference_position:Position,
}}
data! {pub struct Projection {
    /// Candidate reference correspondence, not biological scheme certification.
    pub matches:Vec<Correspondence>,pub deleted_positions:Vec<Position>,
    /// Query insertions deliberately receive no invented scheme labels.
    pub unassigned_residues:Vec<usize>,
}}
data! {pub struct Mapping {
    pub locus:Locus,pub model_revision:String,pub alignment:profile::Report,
    pub projection:Option<Projection>,pub projection_issue:Option<ProjectionIssue>,
}}
data! {pub struct Report {
    pub engine:String,pub source_revision:String,pub corpus_revision:String,pub context_revision:String,
    pub scoring:Scoring,pub limits:Limits,pub cells:usize,pub mappings:Vec<Mapping>,
}}
/// Return every locus independently; uncalibrated scores of differently sized
/// profiles must not choose a winning chain type or imply that a domain exists.
pub fn map_query(
    source: &Molecule,
    scoring: &Scoring,
    limits: &Limits,
    cancel: &mut impl FnMut() -> bool,
) -> Result<Report> {
    if cancel() {
        return Err(Error::Cancelled);
    }
    let source_revision = source_revision(source)?;
    scoring.validate()?;
    if limits.max_cells == 0 || limits.max_cells > profile::MAX_CELLS {
        return Err(Error::Limit("V-reference aggregate cells"));
    }
    let models = models(scoring, cancel)?;
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
        b"antibody-v-reference-context-1\0",
        &(
            MAPPING_ENGINE,
            profile::PROFILE_ENGINE,
            &source_revision,
            &corpus_revision,
            scoring,
            limits,
            "distinct protein weighting; C23/C104 required; all-gap positions omitted; no classification; no inferred insertion labels; ambiguous/tied projections withheld",
        ),
        65536,
    )?;
    let mut mappings = vec![];
    for model in models {
        let model_revision = model_revision(&model)?;
        let individual_cells = (source.sequence.len() + 1) * (model.profile.states.len() + 1);
        let alignment = profile::align(
            source,
            &model.profile,
            &profile::Limits {
                max_cells: individual_cells,
            },
            cancel,
        )?;
        let projection_issue = if alignment.upper.is_none() {
            Some(ProjectionIssue::NoPath)
        } else if alignment.threshold_status == profile::ThresholdStatus::BelowThreshold {
            Some(ProjectionIssue::BelowScoreThreshold)
        } else if alignment.ambiguous_residues > 0 {
            Some(ProjectionIssue::AmbiguousSource)
        } else if alignment.upper.as_ref().unwrap().alternative_optimal_paths {
            Some(ProjectionIssue::TiedOptimalPaths)
        } else {
            None
        };
        let projection = if projection_issue.is_none() {
            let mut matches = vec![];
            let mut deleted_positions = vec![];
            let mut unassigned_residues = vec![];
            for step in &alignment.upper.as_ref().unwrap().steps {
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
        mappings.push(Mapping {
            locus: model.locus,
            model_revision,
            alignment,
            projection,
            projection_issue,
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
        scoring: scoring.clone(),
        limits: limits.clone(),
        cells,
        mappings,
    })
}
impl Report {
    pub fn validate_for(
        &self,
        source: &Molecule,
        scoring: &Scoring,
        limits: &Limits,
    ) -> Result<()> {
        if self.engine != MAPPING_ENGINE
            || self.mappings.len() != 3
            || self.mappings.iter().any(|m| {
                m.projection.as_ref().is_some_and(|p| {
                    p.matches.len() > MAX_SOURCE_RESIDUES
                        || p.deleted_positions.len() > 104
                        || p.unassigned_residues.len() > MAX_SOURCE_RESIDUES
                })
            })
            || *self != map_query(source, scoring, limits, &mut || false)?
        {
            return Err(Error::Binding);
        }
        Ok(())
    }
}
