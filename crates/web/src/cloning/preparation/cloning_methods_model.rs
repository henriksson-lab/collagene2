//! Editable preparation fields and stable source identities, independent of UI.
use crate::assembly::draft::assembly_role_ids::Role;
use crate::cloning::preparation::golden_gate_pool_form::PoolForm;
use collagene_core::{
    chemistry::{AnnealOptions, EndRepairMethod, EndTreatment, FivePrime, Terminus, ThreePrime},
    cloning::{
        ligation::{
            tailing::TailingBase,
            tailing_recipe::{Action, Profile},
        },
        overlap::profiles::{OverlapMethod, ProfileOptions},
    },
    sequence::analysis::Primer,
};
use serde::{Deserialize, Serialize};

#[derive(Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum Mode {
    Tailing,
    Overlap,
    Ligation,
    GoldenGate,
    MultiSiteReference,
}
#[derive(Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum Treatment {
    Keep,
    Phosphorylate,
    Dephosphorylate,
}
impl Treatment {
    pub(crate) fn operation(self) -> Option<EndTreatment> {
        match self {
            Self::Keep => None,
            Self::Phosphorylate => Some(EndTreatment::Phosphorylate),
            Self::Dephosphorylate => Some(EndTreatment::Dephosphorylate),
        }
    }
}
#[derive(Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct TailInput {
    pub action: Action,
    pub treatment: Treatment,
}
impl TailInput {
    pub(crate) fn new(action: Action) -> Self {
        Self {
            action,
            treatment: Treatment::Keep,
        }
    }
}
#[derive(Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum OverlapPreparation {
    WholeLinear,
    DeclaredDuplex,
    Digest,
    CatalogDigest,
    ExactPcr,
    PcrCatalogDigest,
    AnnealedOligos,
    MismatchPcr,
    MismatchPcrCatalogDigest,
    OverlapExtensionPcr,
    OverlapExtensionPcrCatalogDigest,
}
/// Editable oligos remain strings until the operator requests reconstruction.
/// Incomplete fields are retained by draft persistence, without scientific status.
#[derive(Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ExtensionPrimerInput {
    pub name: String,
    pub annealing: String,
    pub tail: String,
}
#[derive(Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ExtensionPrimerPairInput {
    pub template_strand: collagene_core::chemistry::Strand,
    pub forward: ExtensionPrimerInput,
    pub reverse: ExtensionPrimerInput,
}
#[derive(Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct OverlapExtensionPcrInput {
    pub left: ExtensionPrimerPairInput,
    pub right: ExtensionPrimerPairInput,
    pub outer: ExtensionPrimerPairInput,
    pub overlap: String,
    pub min_paired_bases: String,
    pub max_comparisons: String,
    pub max_candidates: String,
}
impl ExtensionPrimerPairInput {
    fn empty(label: &str) -> Self {
        let primer = |direction: &str| ExtensionPrimerInput {
            name: format!("{label} {direction}"),
            annealing: String::new(),
            tail: String::new(),
        };
        Self {
            template_strand: collagene_core::chemistry::Strand::Top,
            forward: primer("forward"),
            reverse: primer("reverse"),
        }
    }
    fn typed(
        &self,
    ) -> Result<collagene_core::cloning::overlap::extension::pcr::PrimerPair, String> {
        Ok(
            collagene_core::cloning::overlap::extension::pcr::PrimerPair {
                template_strand: self.template_strand,
                forward: Primer::new(
                    &self.forward.name,
                    &self.forward.annealing,
                    &self.forward.tail,
                )
                .map_err(|e| e.to_string())?,
                reverse: Primer::new(
                    &self.reverse.name,
                    &self.reverse.annealing,
                    &self.reverse.tail,
                )
                .map_err(|e| e.to_string())?,
            },
        )
    }
}
impl Default for OverlapExtensionPcrInput {
    fn default() -> Self {
        let options = collagene_core::cloning::overlap::extension::Options::default();
        Self {
            left: ExtensionPrimerPairInput::empty("Left PCR"),
            right: ExtensionPrimerPairInput::empty("Right PCR"),
            outer: ExtensionPrimerPairInput::empty("Outer PCR"),
            overlap: options.overlap.to_string(),
            min_paired_bases: options.min_paired_bases.to_string(),
            max_comparisons: options.max_comparisons.to_string(),
            max_candidates: options.max_candidates.to_string(),
        }
    }
}
impl OverlapExtensionPcrInput {
    /// Bounds only: typing an incomplete number or oligo must remain saveable.
    pub(crate) fn validate_raw(&self) -> Result<(), String> {
        for pair in [&self.left, &self.right, &self.outer] {
            for primer in [&pair.forward, &pair.reverse] {
                if primer.name.len() > 256
                    || primer.annealing.len() > 1000
                    || primer.tail.len() > 1000
                {
                    return Err(
                        "Overlap-extension primer fields exceed their supported bounds.".into(),
                    );
                }
            }
        }
        if [
            &self.overlap,
            &self.min_paired_bases,
            &self.max_comparisons,
            &self.max_candidates,
        ]
        .iter()
        .any(|s| s.len() > 20)
        {
            return Err("Overlap-extension numeric fields exceed20 characters.".into());
        }
        Ok(())
    }
    pub(crate) fn request(
        &self,
    ) -> Result<collagene_core::cloning::overlap::extension::pcr::Request, String> {
        self.validate_raw()?;
        let number = |text: &str, label: &str| {
            if text.is_empty() || !text.bytes().all(|byte| byte.is_ascii_digit()) {
                return Err(format!("Enter a positive whole number for {label}."));
            }
            text.parse::<usize>()
                .map_err(|_| format!("{label} exceeds the supported range."))
        };
        let request = collagene_core::cloning::overlap::extension::pcr::Request {
            left: self.left.typed()?,
            right: self.right.typed()?,
            outer: self.outer.typed()?,
            extension: collagene_core::cloning::overlap::extension::Options {
                overlap: number(&self.overlap, "selected terminal homology")?,
                min_paired_bases: number(&self.min_paired_bases, "minimum terminal pairing")?,
                max_comparisons: number(&self.max_comparisons, "terminal comparison limit")?,
                max_candidates: number(&self.max_candidates, "terminal candidate limit")?,
            },
        };
        request.validate_shape().map_err(|e| e.to_string())?;
        Ok(request)
    }
}
/// Editable source declarations; invalid unfinished text remains draft data only.
#[derive(Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SourceGenerationInput {
    pub dam: String,
    pub dcm: String,
    pub cpg: String,
}
impl Default for SourceGenerationInput {
    fn default() -> Self {
        Self {
            dam: "unknown".into(),
            dcm: "unknown".into(),
            cpg: "unknown".into(),
        }
    }
}
impl SourceGenerationInput {
    pub(crate) fn validate_raw(&self) -> Result<(), String> {
        if [&self.dam, &self.dcm, &self.cpg]
            .iter()
            .any(|s| s.len() > 32)
        {
            return Err("Source generation fields exceed32 bytes.".into());
        }
        Ok(())
    }
    pub(crate) fn typed(
        &self,
    ) -> Result<collagene_core::chemistry::tracked_duplex::GenerationAssumptions, String> {
        use collagene_core::chemistry::{
            methylation::MethylationState as S, tracked_duplex::GenerationAssumptions,
        };
        self.validate_raw()?;
        let parse = |v: &str| {
            match v {"unknown"=>Ok(S::Unknown),"unmodified"=>Ok(S::Unmodified),"fully"=>Ok(S::FullyModified),"mixed"=>Ok(S::HemiOrMixed),_=>Err("Complete the explicit Dam/Dcm/CpG source-generation choices before preparing material.".to_string())}
        };
        Ok(GenerationAssumptions {
            dam: parse(&self.dam)?,
            dcm: parse(&self.dcm)?,
            cpg: parse(&self.cpg)?,
        })
    }
}

#[derive(Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct OverlapInput {
    #[serde(deserialize_with = "required_option")]
    pub run_origin: Option<crate::history::run_origin_model::PinnedOrigin>,
    #[serde(deserialize_with = "required_option")]
    pub source_generation: Option<SourceGenerationInput>,
    pub document: usize,
    pub interval_start: String,
    pub interval_end: String,
    pub reverse: bool,
    pub treatment: Treatment,
    pub preparation: OverlapPreparation,
    pub template_strand: collagene_core::chemistry::Strand,
    #[serde(deserialize_with = "required_option")]
    pub overlap_extension_pcr: Option<OverlapExtensionPcrInput>,
    #[serde(deserialize_with = "required_option")]
    pub catalog_selection: Option<collagene_core::cloning::restriction::catalog::CatalogSelection>,
    pub left_terminus: Terminus,
    pub right_terminus: Terminus,
    pub enzymes: String,
    pub fragment: String,
    #[serde(deserialize_with = "required_option")]
    pub repair: Option<EndRepairMethod>,
    pub forward_name: String,
    pub reverse_name: String,
    pub forward: String,
    pub reverse_primer: String,
    pub forward_tail: String,
    pub reverse_tail: String,
    pub mismatch_max_mismatches: String,
    pub mismatch_three_prime_exact: String,
    pub mismatch_max_hits: String,
    pub mismatch_max_work: String,
    pub mismatch_max_pairs: String,
    pub mismatch_max_products: String,
    pub mismatch_max_product_bases: String,
    pub right_fusion: String,
    pub anneal_options: AnnealOptions,
    #[serde(deserialize_with = "required_option")]
    pub anneal_candidate: Option<usize>,
    #[serde(deserialize_with = "required_option")]
    pub anneal_bottom_document: Option<usize>,
    /// Documentary selection context only; physical authority requires replay.
    #[serde(deserialize_with = "required_option")]
    pub anneal_source_context: Option<String>,
}
fn default_mismatch_max() -> String {
    collagene_core::cloning::pcr::mutagenesis::PcrOptions::default()
        .binding
        .max_mismatches
        .to_string()
}
fn default_mismatch_exact() -> String {
    collagene_core::cloning::pcr::mutagenesis::PcrOptions::default()
        .binding
        .min_three_prime_exact
        .to_string()
}

fn default_mismatch_hits() -> String {
    collagene_core::cloning::pcr::mutagenesis::PcrOptions::default()
        .binding
        .max_hits
        .to_string()
}

fn default_mismatch_work() -> String {
    collagene_core::cloning::pcr::mutagenesis::PcrOptions::default()
        .binding
        .max_work
        .to_string()
}

fn default_mismatch_pairs() -> String {
    collagene_core::cloning::pcr::mutagenesis::PcrOptions::default()
        .max_pair_checks
        .to_string()
}

fn default_mismatch_products() -> String {
    collagene_core::cloning::pcr::mutagenesis::PcrOptions::default()
        .max_products
        .to_string()
}

fn default_mismatch_bases() -> String {
    collagene_core::cloning::pcr::mutagenesis::PcrOptions::default()
        .max_product_bases
        .to_string()
}

impl Default for OverlapInput {
    fn default() -> Self {
        Self {
            template_strand: collagene_core::chemistry::Strand::Top,
            run_origin: None,
            source_generation: None,
            document: 0,
            interval_start: "1".into(),
            interval_end: String::new(),
            reverse: false,
            treatment: Treatment::Keep,
            preparation: OverlapPreparation::WholeLinear,
            overlap_extension_pcr: None,
            catalog_selection: None,
            left_terminus: Terminus::UNKNOWN,
            right_terminus: Terminus::UNKNOWN,
            enzymes: "EcoRV".into(),
            fragment: "1".into(),
            repair: None,
            forward_name: "Forward".into(),
            reverse_name: "Reverse".into(),
            forward: String::new(),
            reverse_primer: String::new(),
            forward_tail: String::new(),
            reverse_tail: String::new(),
            mismatch_max_mismatches: default_mismatch_max(),
            mismatch_three_prime_exact: default_mismatch_exact(),
            mismatch_max_hits: default_mismatch_hits(),
            mismatch_max_work: default_mismatch_work(),
            mismatch_max_pairs: default_mismatch_pairs(),
            mismatch_max_products: default_mismatch_products(),
            mismatch_max_product_bases: default_mismatch_bases(),
            right_fusion: String::new(),
            anneal_options: AnnealOptions::default(),
            anneal_candidate: None,
            anneal_bottom_document: None,
            anneal_source_context: None,
        }
    }
}
impl OverlapInput {
    /// Original declaration settings remain inactive while an exact saved origin
    /// supplies its own freshly reconstructed positional state.
    pub(crate) fn generation(
        &self,
    ) -> Result<Option<collagene_core::chemistry::tracked_duplex::GenerationAssumptions>, String>
    {
        if let Some(raw) = &self.source_generation {
            raw.validate_raw()?;
        }
        if self.run_origin.is_some() {
            Ok(None)
        } else {
            self.source_generation
                .as_ref()
                .map(SourceGenerationInput::typed)
                .transpose()
        }
    }

    pub(crate) fn is_pcr(&self) -> bool {
        matches!(
            self.preparation,
            OverlapPreparation::ExactPcr
                | OverlapPreparation::PcrCatalogDigest
                | OverlapPreparation::MismatchPcr
                | OverlapPreparation::MismatchPcrCatalogDigest
        )
    }
    pub(crate) fn is_overlap_extension_pcr(&self) -> bool {
        matches!(
            self.preparation,
            OverlapPreparation::OverlapExtensionPcr
                | OverlapPreparation::OverlapExtensionPcrCatalogDigest
        )
    }
    pub(crate) fn has_pcr_amplicon(&self) -> bool {
        self.is_pcr() || self.is_overlap_extension_pcr()
    }
    /// Select the complete PCR product without replacing the declared PCR policy.
    pub(crate) fn use_complete_pcr_amplicon(&mut self) -> Result<(), String> {
        if !self.has_pcr_amplicon() {
            return Err(
                "Choose an explicit PCR preparation before using its complete amplicon.".into(),
            );
        }
        if self.preparation == OverlapPreparation::PcrCatalogDigest {
            self.preparation = OverlapPreparation::ExactPcr;
        } else if self.preparation == OverlapPreparation::MismatchPcrCatalogDigest {
            self.preparation = OverlapPreparation::MismatchPcr;
        } else if self.preparation == OverlapPreparation::OverlapExtensionPcrCatalogDigest {
            self.preparation = OverlapPreparation::OverlapExtensionPcr;
        }
        self.fragment = "1".into();
        Ok(())
    }
    pub(crate) fn mismatch_options(
        &self,
    ) -> Result<collagene_core::cloning::pcr::mutagenesis::PcrOptions, String> {
        fn number(value: &str, name: &str) -> Result<usize, String> {
            if value.is_empty()
                || value.len() > 20
                || !value.bytes().all(|byte| byte.is_ascii_digit())
            {
                return Err(format!("Enter a nonnegative whole number for {name}."));
            }
            value
                .parse::<usize>()
                .map_err(|_| format!("{name} exceeds the supported range."))
        }
        let max_mismatches = number(&self.mismatch_max_mismatches, "allowed mismatches")?;
        let min_three_prime_exact = number(&self.mismatch_three_prime_exact, "exact 3′ bases")?;
        let domain = self.forward.len().min(self.reverse_primer.len());
        if domain == 0 || domain > 1000 || max_mismatches > domain || min_three_prime_exact > domain
        {
            return Err("Mismatch and exact 3′ counts must fit both explicit primer annealing domains (1–1000 bases).".into());
        }
        let mut options = collagene_core::cloning::pcr::mutagenesis::PcrOptions::default();
        options.binding.max_mismatches = max_mismatches;
        options.binding.min_three_prime_exact = min_three_prime_exact;
        options.binding.max_hits = number(&self.mismatch_max_hits, "maximum primer hits")?;
        options.binding.max_work = number(&self.mismatch_max_work, "binding work limit")?;
        options.max_pair_checks = number(&self.mismatch_max_pairs, "primer pair check limit")?;
        options.max_products = number(&self.mismatch_max_products, "PCR product limit")?;
        options.max_product_bases = number(
            &self.mismatch_max_product_bases,
            "total PCR product base limit",
        )?;
        if options.binding.max_hits == 0
            || options.binding.max_hits > 10_000
            || options.binding.max_work == 0
            || options.binding.max_work > 200_000_000
            || options.max_pair_checks == 0
            || options.max_pair_checks > 1_000_000
            || options.max_products == 0
            || options.max_products > 10_000
            || options.max_product_bases == 0
            || options.max_product_bases > 32_000_000
        {
            return Err(
                "PCR search limits must be positive and inside the recorded native profile.".into(),
            );
        }
        Ok(options)
    }
    pub(crate) fn summary(&self) -> String {
        let preparation = match self.preparation {
            OverlapPreparation::WholeLinear => {
                "Whole linear source (circular DNA requires preparation)".into()
            }
            OverlapPreparation::DeclaredDuplex => format!(
                "Declared intact linear duplex assumption; left {}; right {}",
                declared_terminus(self.left_terminus),
                declared_terminus(self.right_terminus)
            ),
            OverlapPreparation::CatalogDigest
            | OverlapPreparation::PcrCatalogDigest
            | OverlapPreparation::MismatchPcrCatalogDigest => format!(
                "{}: {}; fragment {}",
                if self.preparation == OverlapPreparation::MismatchPcrCatalogDigest {
                    "Intentional-mismatch PCR followed by pinned catalog digest"
                } else if self.preparation == OverlapPreparation::PcrCatalogDigest {
                    "Exact PCR followed by pinned catalog digest"
                } else {
                    "Pinned catalog digest"
                },
                self.catalog_selection
                    .as_ref()
                    .map(|s| s
                        .enzymes
                        .iter()
                        .map(|e| e.enzyme.name.as_str())
                        .collect::<Vec<_>>()
                        .join(", "))
                    .unwrap_or_else(|| "select a virtual digest fragment".into()),
                self.fragment
            ),
            OverlapPreparation::Digest => format!(
                "Restriction digest: {}; fragment {} (1-based); {}",
                self.enzymes,
                self.fragment,
                match self.repair {
                    None => "cohesive ends retained",
                    Some(EndRepairMethod::FillFivePrimeOverhangs) => "explicit 5′ overhang fill-in",
                    Some(EndRepairMethod::RemoveThreePrimeOverhangs) =>
                        "explicit 3′ overhang removal",
                    Some(EndRepairMethod::Blunt) => "explicit blunting",
                }
            ),
            OverlapPreparation::MismatchPcr => format!(
                "Intentional-mismatch PCR: annealing {} / {} bp; 5′ extensions {} / {} bp; at most {} mismatches, {} exact 3′ bases per primer",
                self.forward.len(),
                self.reverse_primer.len(),
                self.forward_tail.len(),
                self.reverse_tail.len(),
                self.mismatch_max_mismatches,
                self.mismatch_three_prime_exact,
            ),
            OverlapPreparation::ExactPcr => format!(
                "Exact PCR: annealing {} / {} bp; 5′ extensions {} / {} bp",
                self.forward.len(),
                self.reverse_primer.len(),
                self.forward_tail.len(),
                self.reverse_tail.len()
            ),
            OverlapPreparation::AnnealedOligos => format!(
                "Anneal original top/bottom oligos · {} paired bases minimum · {}",
                self.anneal_options.min_paired_bases,
                self.anneal_candidate
                    .map(|index| format!("explicit candidate {}", index + 1))
                    .unwrap_or_else(|| "choose an annealing candidate".into())
            ),
            OverlapPreparation::OverlapExtensionPcr
            | OverlapPreparation::OverlapExtensionPcrCatalogDigest => format!(
                "Two source PCRs → reciprocal extension ({} bp homology) → outer-primer PCR{}",
                self.overlap_extension_pcr
                    .as_ref()
                    .map(|request| request.overlap.as_str())
                    .unwrap_or("enter"),
                if self.preparation == OverlapPreparation::OverlapExtensionPcrCatalogDigest {
                    format!(
                        " → pinned catalog digest: {}; fragment {}",
                        self.catalog_selection
                            .as_ref()
                            .map(|selection| selection
                                .enzymes
                                .iter()
                                .map(|entry| entry.enzyme.name.as_str())
                                .collect::<Vec<_>>()
                                .join(", "))
                            .unwrap_or_else(|| "choose enzymes".into()),
                        self.fragment
                    )
                } else {
                    String::new()
                }
            ),
        };
        format!(
            "{} · {}",
            preparation,
            if self.reverse {
                "reverse complement after preparation"
            } else {
                "forward after preparation"
            }
        )
    }
}
/// Editable generic 5′ resection policy, never evidence of a completed reaction.
#[derive(Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct HomologyForm {
    pub enabled: bool,
    pub exact_homology: bool,
    pub cuts: String,
    pub all_accessible_ends: bool,
    pub outer_left: String,
    pub outer_right: String,
    pub min_paired_bases: String,
    pub max_comparisons: String,
    pub max_candidates: String,
}
impl Default for HomologyForm {
    fn default() -> Self {
        Self {
            enabled: true,
            exact_homology: true,
            cuts: String::new(),
            all_accessible_ends: false,
            outer_left: "0".into(),
            outer_right: "0".into(),
            min_paired_bases: "12".into(),
            max_comparisons: "20000000".into(),
            max_candidates: "100".into(),
        }
    }
}
impl HomologyForm {
    pub(crate) fn validate_raw(&self) -> Result<(), String> {
        if self.cuts.len() > 1024
            || [
                &self.outer_left,
                &self.outer_right,
                &self.min_paired_bases,
                &self.max_comparisons,
                &self.max_candidates,
            ]
            .iter()
            .any(|s| s.len() > 20)
        {
            return Err(
                "Strand reconstruction settings exceed their editable field bounds.".into(),
            );
        }
        Ok(())
    }
    pub(crate) fn policy(
        &self,
        overlaps: Vec<usize>,
        circular: bool,
        rows: usize,
    ) -> Result<collagene_core::cloning::overlap::assembly::Options, String> {
        use collagene_core::cloning::overlap::assembly::{
            EndCuts, Options, OuterEndPolicy, Resection,
        };
        self.validate_raw()?;
        let number = |s: &str| -> Result<usize, String> {
            if s.is_empty() || !s.bytes().all(|b| b.is_ascii_digit()) {
                return Err(
                    "Complete the strand reconstruction settings with nonnegative whole numbers."
                        .into(),
                );
            }
            s.parse()
                .map_err(|_| "A strand reconstruction number is too large.".into())
        };
        let resection = if self.exact_homology {
            Resection::ExactHomology
        } else {
            let cuts = self
                .cuts
                .split(',')
                .map(|pair| {
                    let (left, right) = pair.trim().split_once(':').ok_or(
                        "Supply one left:right 5′ cut pair per fragment, separated by commas.",
                    )?;
                    Ok(EndCuts {
                        left: number(left.trim())?,
                        right: number(right.trim())?,
                    })
                })
                .collect::<Result<Vec<_>, String>>()?;
            if cuts.len() != rows {
                return Err(
                    "Supply exactly one 5′ cut pair per oriented prepared fragment.".into(),
                );
            }
            Resection::Explicit { cuts }
        };
        let outer_end_policy = if !circular && self.all_accessible_ends {
            OuterEndPolicy::AllAccessibleEnds {
                left: number(&self.outer_left)?,
                right: number(&self.outer_right)?,
            }
        } else {
            OuterEndPolicy::JoiningEndsOnly
        };
        let policy = Options {
            overlaps,
            topology: if circular {
                collagene_model::Topology::Circular
            } else {
                collagene_model::Topology::Linear
            },
            resection,
            outer_end_policy,
            min_paired_bases: number(&self.min_paired_bases)?,
            max_comparisons: number(&self.max_comparisons)?,
            max_candidates: number(&self.max_candidates)?,
        };
        policy.validate_shape(rows).map_err(|e| e.to_string())?;
        Ok(policy)
    }
}
/// Editable explicit 3′ cuts. No gap-fill, nick-seal or host-rescue assumption.
#[derive(Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct AnnealingForm {
    pub enabled: bool,
    pub planned_ecoli_rescue: bool,
    pub cuts: String,
    pub min_paired_bases: String,
    pub max_comparisons: String,
    pub max_candidates: String,
}
impl Default for AnnealingForm {
    fn default() -> Self {
        Self {
            enabled: true,
            planned_ecoli_rescue: false,
            cuts: String::new(),
            min_paired_bases: "12".into(),
            max_comparisons: "20000000".into(),
            max_candidates: "100".into(),
        }
    }
}
impl AnnealingForm {
    pub(crate) fn validate_raw(&self) -> Result<(), String> {
        if self.cuts.len() > 1024
            || [
                &self.min_paired_bases,
                &self.max_comparisons,
                &self.max_candidates,
            ]
            .iter()
            .any(|s| s.len() > 20)
        {
            return Err("3′ annealing settings exceed their editable field bounds.".into());
        }
        Ok(())
    }
    pub(crate) fn policy(
        &self,
        overlaps: Vec<usize>,
        circular: bool,
        rows: usize,
    ) -> Result<collagene_core::cloning::overlap::annealing::Options, String> {
        use collagene_core::cloning::overlap::{annealing::Options, resection::Cuts};
        use collagene_model::Topology;
        self.validate_raw()?;
        if rows == 1 && !circular {
            return Err("A single linear source has no annealing junction. Use an explicit 3′ resection step for standalone source processing.".into());
        }
        let number = |s: &str| -> Result<usize, String> {
            if s.is_empty() || !s.bytes().all(|b| b.is_ascii_digit()) {
                return Err(
                    "Complete 3′ cuts and annealing limits with nonnegative whole numbers.".into(),
                );
            }
            s.parse()
                .map_err(|_| "A 3′ annealing number is too large.".into())
        };
        let cuts =
            self.cuts
                .split(',')
                .map(|pair| {
                    let (left, right) = pair.trim().split_once(':').ok_or(
                        "Supply one left:right 3′ cut pair per oriented prepared fragment.",
                    )?;
                    Ok(Cuts {
                        left: number(left.trim())?,
                        right: number(right.trim())?,
                    })
                })
                .collect::<Result<Vec<_>, String>>()?;
        if cuts.len() != rows {
            return Err("Supply exactly one 3′ cut pair per oriented prepared fragment.".into());
        }
        let policy = Options {
            planned_rescue: self.planned_ecoli_rescue.then_some(
                collagene_core::cloning::overlap::annealing::PlannedRescue::EcoliTransformation,
            ),
            overlaps,
            topology: if circular {
                Topology::Circular
            } else {
                Topology::Linear
            },
            cuts,
            min_paired_bases: number(&self.min_paired_bases)?,
            max_comparisons: number(&self.max_comparisons)?,
            max_candidates: number(&self.max_candidates)?,
        };
        policy.validate_shape(rows).map_err(|e| e.to_string())?;
        Ok(policy)
    }
}
#[derive(Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct MethodsState {
    pub mode: Mode,
    pub vector: TailInput,
    pub insert: TailInput,
    pub tailing_profile: Profile,
    pub circular: bool,
    pub allow_nicked: bool,
    #[serde(deserialize_with = "deserialize_method_rows")]
    pub rows: Vec<OverlapInput>,
    pub overlaps: String,
    pub method: OverlapMethod,
    pub scan_repeats: bool,
    pub comparisons: String,
    pub occurrences: String,
    pub golden_gate: GoldenGateForm,
    /// Editable pool/reaction settings and explicit candidate intent only;
    /// inspections, plans and scientific authority remain transient.
    #[serde(deserialize_with = "required_option")]
    pub golden_gate_pool: Option<PoolForm>,
    #[serde(deserialize_with = "required_option")]
    pub multisite: Option<crate::cloning::methods::multisite_form::Form>,
    /// Stable source-choice and bin identities; explicit identities are never repinned.
    #[serde(deserialize_with = "deserialize_source_roles")]
    pub source_roles: Option<Vec<Role>>,
    #[serde(deserialize_with = "required_option")]
    pub homology: Option<HomologyForm>,
    #[serde(deserialize_with = "required_option")]
    pub annealing: Option<AnnealingForm>,
}
#[derive(Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct GoldenGateForm {
    pub enzyme: collagene_core::cloning::golden_gate::GoldenGateEnzyme,
    pub clamp: String,
    pub spacer: String,
    pub forward_length: String,
    pub reverse_length: String,
}
impl Default for GoldenGateForm {
    fn default() -> Self {
        Self {
            enzyme: collagene_core::cloning::golden_gate::GoldenGateEnzyme::BsaI,
            clamp: "ATGC".into(),
            spacer: "A".into(),
            forward_length: "20".into(),
            reverse_length: "20".into(),
        }
    }
}
fn deserialize_method_rows<'de, D: serde::Deserializer<'de>>(
    deserializer: D,
) -> Result<Vec<OverlapInput>, D::Error> {
    struct Visitor;
    impl<'de> serde::de::Visitor<'de> for Visitor {
        type Value = Vec<OverlapInput>;
        fn expecting(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
            f.write_str("at most12 assembly preparation rows")
        }
        fn visit_seq<A: serde::de::SeqAccess<'de>>(
            self,
            mut seq: A,
        ) -> Result<Self::Value, A::Error> {
            if seq.size_hint().is_some_and(|n| n > 12) {
                return Err(serde::de::Error::custom(
                    "at most12 assembly preparation rows",
                ));
            }
            let mut rows = Vec::with_capacity(2);
            while rows.len() < 12 {
                match seq.next_element()? {
                    Some(row) => rows.push(row),
                    None => return Ok(rows),
                }
            }
            if seq.next_element::<serde::de::IgnoredAny>()?.is_some() {
                return Err(serde::de::Error::custom(
                    "at most12 assembly preparation rows",
                ));
            }
            Ok(rows)
        }
    }
    deserializer.deserialize_seq(Visitor)
}
fn deserialize_source_roles<'de, D: serde::Deserializer<'de>>(
    deserializer: D,
) -> Result<Option<Vec<Role>>, D::Error> {
    struct Optional;
    struct Roles;
    impl<'de> serde::de::Visitor<'de> for Roles {
        type Value = Vec<Role>;
        fn expecting(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
            f.write_str("at most12 stable assembly identities")
        }
        fn visit_seq<A: serde::de::SeqAccess<'de>>(
            self,
            mut seq: A,
        ) -> Result<Self::Value, A::Error> {
            if seq.size_hint().is_some_and(|n| n > 12) {
                return Err(serde::de::Error::custom(
                    "at most12 stable assembly identities",
                ));
            }
            let mut roles = Vec::with_capacity(2);
            while roles.len() < 12 {
                match seq.next_element()? {
                    Some(role) => roles.push(role),
                    None => return Ok(roles),
                }
            }
            if seq.next_element::<serde::de::IgnoredAny>()?.is_some() {
                return Err(serde::de::Error::custom(
                    "at most12 stable assembly identities",
                ));
            }
            Ok(roles)
        }
    }
    impl<'de> serde::de::Visitor<'de> for Optional {
        type Value = Option<Vec<Role>>;
        fn expecting(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
            f.write_str("optional stable assembly identities")
        }
        fn visit_none<E: serde::de::Error>(self) -> Result<Self::Value, E> {
            Ok(None)
        }
        fn visit_unit<E: serde::de::Error>(self) -> Result<Self::Value, E> {
            Ok(None)
        }
        fn visit_some<D: serde::Deserializer<'de>>(
            self,
            deserializer: D,
        ) -> Result<Self::Value, D::Error> {
            deserializer.deserialize_seq(Roles).map(Some)
        }
    }
    deserializer.deserialize_option(Optional)
}
impl Default for MethodsState {
    fn default() -> Self {
        Self {
            mode: Mode::Overlap,
            vector: TailInput::new(Action::AddTail {
                base: TailingBase::T,
            }),
            insert: TailInput::new(Action::AddTail {
                base: TailingBase::A,
            }),
            tailing_profile: Profile::Ta,
            circular: true,
            allow_nicked: false,
            rows: vec![OverlapInput::default(), OverlapInput::default()],
            overlaps: "20".into(),
            method: OverlapMethod::HiFi,
            scan_repeats: true,
            comparisons: "5000000".into(),
            occurrences: "1000".into(),
            golden_gate: GoldenGateForm::default(),
            golden_gate_pool: None,
            multisite: None,
            source_roles: None,
            homology: None,
            annealing: None,
        }
    }
}
impl MethodsState {
    pub(crate) fn cleavage_policy(
        &self,
    ) -> collagene_core::cloning::assembly::recipe::CleavagePolicy {
        if self.mode == Mode::GoldenGate {
            collagene_core::cloning::assembly::recipe::CleavagePolicy::PositionalChecked
        } else {
            collagene_core::cloning::assembly::recipe::CleavagePolicy::AssumedComplete
        }
    }
    pub(crate) fn validate_processing_selection(&self) -> Result<(), String> {
        if self.homology.as_ref().is_some_and(|f| f.enabled)
            && self.annealing.as_ref().is_some_and(|f| f.enabled)
        {
            return Err(
                "Choose one explicit processing route: 5′ reconstruction or 3′ annealing.".into(),
            );
        }
        Ok(())
    }
    pub(crate) fn bottom_input_id(&self, row: usize) -> Result<String, String> {
        if self
            .rows
            .get(row)
            .is_none_or(|input| input.preparation != OverlapPreparation::AnnealedOligos)
        {
            return Err("Choose an active paired-oligo preparation row.".into());
        }
        let roles = self.role_ids(self.rows.len())?;
        collagene_core::cloning::assembly::recipe::anneal_bottom_input_id(&roles[row].input_id)
            .map_err(|error| error.to_string())
    }
    pub(crate) fn all_input_ids(&self) -> Result<Vec<String>, String> {
        let mut ids = self
            .role_ids(self.rows.len())?
            .into_iter()
            .map(|role| role.input_id)
            .collect::<Vec<_>>();
        for (index, row) in self.rows.iter().enumerate() {
            if row.preparation == OverlapPreparation::AnnealedOligos {
                ids.push(self.bottom_input_id(index)?);
            }
        }
        if ids.iter().collect::<std::collections::BTreeSet<_>>().len() != ids.len() {
            return Err("Primary and paired-oligo input identities collide.".into());
        }
        Ok(ids)
    }
    pub(crate) fn source_indices(&self) -> Result<Vec<usize>, String> {
        self.all_input_ids()?;
        let mut indices = self.rows.iter().map(|row| row.document).collect::<Vec<_>>();
        for row in &self.rows {
            if row.preparation == OverlapPreparation::AnnealedOligos {
                indices.push(row.anneal_bottom_document.ok_or(
                    "Capture the original bottom oligo before reconstructing its preparation.",
                )?);
            }
        }
        Ok(indices)
    }
    pub(crate) fn role_ids(&self, count: usize) -> Result<Vec<Role>, String> {
        if self.mode == Mode::MultiSiteReference && self.source_roles.is_none() {
            return crate::cloning::methods::multisite_form::ordinary_roles(count);
        }
        crate::assembly::draft::assembly_role_ids::resolve(&self.source_roles, count)
    }
    pub(crate) fn preparation_context(&self) -> Self {
        let mut context = self.clone();
        context.golden_gate_pool = None;
        context
    }
    pub(crate) fn pool_form(&self) -> PoolForm {
        self.golden_gate_pool.clone().unwrap_or_else(|| PoolForm {
            circular: self.circular,
            choices: vec![0; self.rows.len()],
            ..PoolForm::default()
        })
    }
    pub(crate) fn swap_rows(&mut self, a: usize, b: usize) {
        if a >= self.rows.len()
            || b >= self.rows.len()
            || (self.mode == Mode::Tailing && self.rows.len() != 2)
            || self.role_ids(self.rows.len()).is_err()
            || self.golden_gate_pool.as_ref().is_some_and(|form| {
                !form.choices.is_empty() && form.choices.len() != self.rows.len()
            })
        {
            return;
        }
        self.rows.swap(a, b);
        if self.mode == Mode::Tailing && a != b {
            // These are the only tail controls used by the source-slot recipe.
            // Preparation/roles travel with rows; sync_rows rebinds the
            // document/reverse fields to the new slot positions afterwards.
            std::mem::swap(&mut self.vector.action, &mut self.insert.action);
            std::mem::swap(&mut self.vector.treatment, &mut self.insert.treatment);
        }
        if let Some(roles) = &mut self.source_roles {
            roles.swap(a, b);
        }
        if let Some(form) = &mut self.golden_gate_pool {
            form.choices.resize(self.rows.len(), 0);
            form.choices.swap(a, b);
            if a != b {
                form.candidate = None;
            }
        }
    }
    pub(crate) fn remove_row(&mut self, index: usize) {
        if index >= self.rows.len()
            || self.rows.len() <= 1
            || self.role_ids(self.rows.len()).is_err()
            || self.golden_gate_pool.as_ref().is_some_and(|form| {
                !form.choices.is_empty() && form.choices.len() != self.rows.len()
            })
        {
            return;
        }
        self.rows.remove(index);
        if let Some(roles) = &mut self.source_roles {
            roles.remove(index);
        }
        if let Some(form) = &mut self.golden_gate_pool {
            form.choices.resize(self.rows.len() + 1, 0);
            form.choices.remove(index);
            form.candidate = None;
        }
    }
    pub(crate) fn push_row(&mut self, row: OverlapInput) -> Result<(), String> {
        self.role_ids(self.rows.len())?;
        if self.source_roles.is_some() {
            return Err(
                "Add an explicit source choice before adding a stable assembly role.".into(),
            );
        }
        if self.rows.len() >= 12
            || self.golden_gate_pool.as_ref().is_some_and(|form| {
                !form.choices.is_empty() && form.choices.len() != self.rows.len()
            })
        {
            return Err("Assembly rows or reaction choices exceed their supported bounds.".into());
        }
        self.rows.push(row);
        if let Some(form) = &mut self.golden_gate_pool {
            form.choices.resize(self.rows.len(), 0);
            form.candidate = None;
        }
        Ok(())
    }
    pub(crate) fn input(&self, vector: bool) -> &TailInput {
        if vector { &self.vector } else { &self.insert }
    }
    pub(crate) fn input_mut(&mut self, vector: bool) -> &mut TailInput {
        if vector {
            &mut self.vector
        } else {
            &mut self.insert
        }
    }
    pub(crate) fn overlap_finishing(
        &self,
        overlaps: Vec<usize>,
        rows: usize,
    ) -> Result<collagene_core::cloning::assembly::recipe::OverlapFinishing, String> {
        use collagene_core::cloning::assembly::recipe::OverlapFinishing;
        let finishing = match (
            self.annealing.as_ref().filter(|form| form.enabled),
            self.homology.as_ref().filter(|form| form.enabled),
        ) {
            (Some(_), Some(_)) => {
                return Err("Choose one strand-processing route for this construct.".into());
            }
            (Some(form), None) => OverlapFinishing::HomologyAnnealing {
                policy: form.policy(overlaps.clone(), self.circular, rows)?,
            },
            (None, Some(form)) => OverlapFinishing::HomologyAssembly {
                policy: form.policy(overlaps.clone(), self.circular, rows)?,
            },
            (None, None) => OverlapFinishing::ReferenceOnly,
        };
        finishing
            .validate_geometry(rows, &overlaps, self.circular, &self.profile_options()?)
            .map_err(|e| e.to_string())?;
        Ok(finishing)
    }
    pub(crate) fn profile(&self) -> Result<(Vec<usize>, ProfileOptions), String> {
        // A linear singleton has no joining boundary. An empty declaration is
        // meaningful only in the explicitly selected physical processing route.
        if self.overlaps.trim().is_empty()
            && self.rows.len() == 1
            && !self.circular
            && self.homology.as_ref().is_some_and(|form| form.enabled)
        {
            return Ok((Vec::new(), self.profile_options()?));
        }
        let overlaps = self
            .overlaps
            .split(',')
            .map(|s| {
                s.trim()
                    .parse::<usize>()
                    .map_err(|_| "Supply comma-separated overlap lengths".to_string())
            })
            .collect::<Result<Vec<_>, _>>()?;
        Ok((overlaps, self.profile_options()?))
    }
    pub(crate) fn profile_options(&self) -> Result<ProfileOptions, String> {
        let options = ProfileOptions {
            method: self.method,
            scan_repeated_overlaps: self.scan_repeats,
            max_comparisons: self
                .comparisons
                .parse()
                .map_err(|_| "Invalid comparison budget")?,
            max_occurrences: self
                .occurrences
                .parse()
                .map_err(|_| "Invalid occurrence limit")?,
        };
        options.validate().map_err(|e| e.to_string())?;
        Ok(options)
    }
}

pub(crate) fn declared_terminus(end: Terminus) -> String {
    format!(
        "5′ {}; 3′ {}",
        match end.five_prime {
            FivePrime::Phosphate => "phosphate",
            FivePrime::Hydroxyl => "hydroxyl",
            FivePrime::Unknown => "unknown",
        },
        match end.three_prime {
            ThreePrime::Hydroxyl => "hydroxyl",
            ThreePrime::Blocked => "blocked",
            ThreePrime::Unknown => "unknown",
        }
    )
}

#[cfg(test)]
mod mismatch_tests {
    use super::*;
    #[test]
    fn current_row_json_retains_explicit_defaults_and_mismatch_policy_bounds() {
        let ordinary = OverlapInput::default();
        let json = serde_json::to_value(&ordinary).unwrap();
        assert!(json.get("mismatch_max_mismatches").is_some());
        assert!(json.get("mismatch_three_prime_exact").is_some());
        let mut row: OverlapInput = serde_json::from_value(json).unwrap();
        row.forward = "ACGT".into();
        row.reverse_primer = "TGCA".into();
        let native = collagene_core::cloning::pcr::mutagenesis::PcrOptions::default();
        let options = row.mismatch_options().unwrap();
        assert_eq!(
            options.binding.max_mismatches,
            native.binding.max_mismatches
        );
        assert_eq!(
            options.binding.min_three_prime_exact,
            native.binding.min_three_prime_exact
        );
        for value in ["-1", "1.5", "", "5", "999999999999999999999"] {
            row.mismatch_max_mismatches = value.into();
            assert!(row.mismatch_options().is_err());
        }
        row.mismatch_max_mismatches = "2".into();
        row.mismatch_three_prime_exact = "3".into();
        assert_eq!(row.mismatch_options().unwrap().binding.max_mismatches, 2);
        row.mismatch_three_prime_exact = "5".into();
        assert!(row.mismatch_options().is_err());
    }
    #[test]
    fn complete_amplicon_keeps_intentional_mismatch_policy_and_rejects_non_pcr_atomically() {
        let mut row = OverlapInput {
            preparation: OverlapPreparation::MismatchPcr,
            mismatch_max_mismatches: "2".into(),
            mismatch_three_prime_exact: "3".into(),
            fragment: "9".into(),
            ..Default::default()
        };
        let before = row.clone();
        row.use_complete_pcr_amplicon().unwrap();
        let mut expected = before;
        expected.fragment = "1".into();
        assert!(row == expected);
        row.preparation = OverlapPreparation::PcrCatalogDigest;
        row.use_complete_pcr_amplicon().unwrap();
        assert!(row.preparation == OverlapPreparation::ExactPcr);
        row.preparation = OverlapPreparation::Digest;
        let before = row.clone();
        assert!(row.use_complete_pcr_amplicon().is_err());
        assert!(row == before);
    }
}

#[cfg(test)]
mod overlap_extension_tests {
    use super::*;
    #[test]
    fn partial_six_oligo_input_roundtrips_current_row_shape() {
        let current = serde_json::to_value(OverlapInput::default()).unwrap();
        assert!(current["overlap_extension_pcr"].is_null());
        let mut row: OverlapInput = serde_json::from_value(current).unwrap();
        let mut raw = OverlapExtensionPcrInput::default();
        for (index, pair) in [&mut raw.left, &mut raw.right, &mut raw.outer]
            .into_iter()
            .enumerate()
        {
            pair.forward.name = format!("Pair {index} forward");
            pair.reverse.name = format!("Pair {index} reverse");
            pair.forward.annealing = format!("incomplete forward {index}");
            pair.reverse.tail = format!("incomplete reverse tail {index}");
        }
        raw.overlap = "2e".into();
        raw.max_comparisons = String::new();
        raw.validate_raw().unwrap();
        assert!(raw.request().is_err());
        row.preparation = OverlapPreparation::OverlapExtensionPcr;
        row.overlap_extension_pcr = Some(raw.clone());
        let decoded: OverlapInput =
            serde_json::from_slice(&serde_json::to_vec(&row).unwrap()).unwrap();
        assert!(decoded == row);
        assert!(decoded.overlap_extension_pcr.as_ref().unwrap() == &raw);
        assert!(!decoded.is_pcr());
        let before = row.clone();
        row.use_complete_pcr_amplicon().unwrap();
        assert!(row == before);
        assert!(row.has_pcr_amplicon());
    }
    #[test]
    fn complete_oe_amplicon_keeps_all_six_oligos_and_catalog_settings() {
        let mut row = OverlapInput::default();
        row.preparation = OverlapPreparation::OverlapExtensionPcrCatalogDigest;
        let mut raw = OverlapExtensionPcrInput::default();
        raw.left.forward.annealing = "unfinished".into();
        raw.outer.forward.tail = "CCGAATTC".into();
        raw.outer.reverse.tail = "GGCTCGAG".into();
        row.overlap_extension_pcr = Some(raw);
        row.fragment = "8".into();
        row.repair = Some(EndRepairMethod::Blunt);
        let catalog = collagene_core::cloning::restriction::catalog::builtin_catalog().unwrap();
        let id = catalog
            .edition()
            .enzymes
            .iter()
            .find(|entry| entry.name == "EcoRI")
            .unwrap()
            .id
            .clone();
        row.catalog_selection =
            Some(collagene_core::cloning::restriction::catalog::select(&catalog, &[id]).unwrap());
        let mut expected = row.clone();
        expected.preparation = OverlapPreparation::OverlapExtensionPcr;
        expected.fragment = "1".into();
        assert!(!row.is_pcr() && row.has_pcr_amplicon());
        row.use_complete_pcr_amplicon().unwrap();
        assert!(row == expected);
    }
    #[test]
    fn raw_bounds_are_separate_from_scientific_oligo_and_search_validity() {
        let mut raw = OverlapExtensionPcrInput::default();
        raw.left.forward.name = String::new();
        raw.overlap = "-".into();
        raw.validate_raw().unwrap();
        assert!(raw.request().is_err());
        for which in 0..3 {
            let mut large = raw.clone();
            let pair = match which {
                0 => &mut large.left,
                1 => &mut large.right,
                _ => &mut large.outer,
            };
            pair.reverse.tail = "A".repeat(1001);
            assert!(large.validate_raw().is_err());
        }
        raw.max_comparisons = "9".repeat(21);
        assert!(raw.validate_raw().is_err());
    }
}

#[cfg(test)]
mod homology_tests {
    use super::*;
    #[test]
    fn physical_linear_singleton_has_no_phantom_junction_and_reference_csv_stays_strict() {
        let mut state = MethodsState::default();
        state.circular = false;
        state.rows.truncate(1);
        state.overlaps.clear();
        assert!(state.profile().is_err());
        state.homology = Some(HomologyForm::default());
        let (overlaps, _) = state.profile().unwrap();
        assert!(overlaps.is_empty());
        let policy = state
            .homology
            .as_ref()
            .unwrap()
            .policy(overlaps, false, 1)
            .unwrap();
        assert!(policy.overlaps.is_empty());
        state.circular = true;
        assert!(state.profile().is_err());
        state.circular = false;
        state.rows.push(OverlapInput::default());
        assert!(state.profile().is_err());
        state.rows.truncate(1);
        state.homology.as_mut().unwrap().enabled = false;
        assert!(state.profile().is_err());
    }
    #[test]
    fn raw_policy_maps_oriented_cuts_outer_tails_and_bounds_without_promoting_partials() {
        use collagene_core::cloning::overlap::assembly::{OuterEndPolicy, Resection};
        let mut form = HomologyForm::default();
        let policy = form.policy(vec![16], false, 2).unwrap();
        assert_eq!(policy.overlaps, vec![16]);
        assert_eq!(policy.topology, collagene_model::Topology::Linear);
        assert!(matches!(policy.resection, Resection::ExactHomology));
        assert!(matches!(
            policy.outer_end_policy,
            OuterEndPolicy::JoiningEndsOnly
        ));
        form.exact_homology = false;
        form.cuts = "2:19, 17:3".into();
        form.all_accessible_ends = true;
        form.outer_left = "2".into();
        form.outer_right = "3".into();
        let policy = form.policy(vec![16], false, 2).unwrap();
        assert!(
            matches!(policy.resection,Resection::Explicit{cuts} if cuts.len()==2&&cuts[0].left==2&&cuts[0].right==19&&cuts[1].left==17&&cuts[1].right==3)
        );
        assert!(matches!(
            policy.outer_end_policy,
            OuterEndPolicy::AllAccessibleEnds { left: 2, right: 3 }
        ));
        form.cuts = "2:unfinished,17:3".into();
        form.validate_raw().unwrap();
        assert!(form.policy(vec![16], false, 2).is_err());
        form.exact_homology = true;
        form.all_accessible_ends = false;
        form.outer_left = "inactive".into();
        assert!(form.policy(vec![16, 16], true, 2).is_ok());
        let current = serde_json::to_value(MethodsState::default()).unwrap();
        assert!(current["homology"].is_null());
    }
}

#[cfg(test)]
mod annealing_tests {
    use super::*;
    #[test]
    fn explicit_three_prime_policy_keeps_cut_orientation_and_never_invents_cuts() {
        let mut form = AnnealingForm::default();
        assert!(form.policy(vec![16], false, 2).is_err());
        form.cuts = "2:19, 17:3".into();
        let policy = form.policy(vec![16], false, 2).unwrap();
        assert_eq!(policy.topology, collagene_model::Topology::Linear);
        assert_eq!(policy.overlaps, vec![16]);
        assert_eq!((policy.cuts[0].left, policy.cuts[0].right), (2, 19));
        assert_eq!((policy.cuts[1].left, policy.cuts[1].right), (17, 3));
        assert_eq!(policy.min_paired_bases, 12);
        assert!(form.policy(vec![16, 16], true, 2).is_ok());
        assert!(form.policy(vec![16], true, 2).is_err());
        assert!(form.policy(vec![], false, 1).is_err());
        form.cuts = "16:16".into();
        assert!(form.policy(vec![16], true, 1).is_ok());
        form.cuts = "2:unfinished,17:3".into();
        assert!(form.validate_raw().is_ok());
        assert!(form.policy(vec![16], false, 2).is_err());
    }
    #[test]
    fn three_prime_work_limits_and_saved_conflicting_modes_fail_closed() {
        let mut form = AnnealingForm::default();
        form.cuts = "0:16,16:0".into();
        for limit in ["0", "-1", "1.5", "100000001", "999999999999999999999"] {
            form.max_comparisons = limit.into();
            assert!(form.policy(vec![16], false, 2).is_err());
        }
        form.max_comparisons = "20000000".into();
        form.max_candidates = "0".into();
        assert!(form.policy(vec![16], false, 2).is_err());
        form.max_candidates = "100".into();
        form.min_paired_bases = "17".into();
        assert!(form.policy(vec![16], false, 2).is_err());
        let mut state = MethodsState::default();
        assert!(serde_json::to_value(&state).unwrap()["annealing"].is_null());
        state.annealing = Some(AnnealingForm::default());
        assert!(state.validate_processing_selection().is_ok());
        state.homology = Some(HomologyForm::default());
        assert!(state.validate_processing_selection().is_err());
        state.annealing.as_mut().unwrap().enabled = false;
        assert!(state.validate_processing_selection().is_ok());
    }
}

fn required_option<'de, D: serde::Deserializer<'de>, T: Deserialize<'de>>(
    deserializer: D,
) -> Result<Option<T>, D::Error> {
    Option::deserialize(deserializer)
}

#[cfg(test)]
mod source_generation_tests {
    use super::*;
    use collagene_core::chemistry::methylation::MethylationState;
    #[test]
    fn defaults_are_unknown_and_partial_text_is_never_a_scientific_declaration() {
        let mut row = OverlapInput::default();
        assert!(row.generation().unwrap().is_none());
        row.source_generation = Some(SourceGenerationInput::default());
        let unknown = row.generation().unwrap().unwrap();
        assert!(
            unknown.dam == MethylationState::Unknown
                && unknown.dcm == MethylationState::Unknown
                && unknown.cpg == MethylationState::Unknown
        );
        row.source_generation.as_mut().unwrap().dam = "unfinished".into();
        row.source_generation
            .as_ref()
            .unwrap()
            .validate_raw()
            .unwrap();
        assert!(row.generation().is_err());
        let restored: OverlapInput =
            serde_json::from_value(serde_json::to_value(&row).unwrap()).unwrap();
        assert!(restored == row);
        row.source_generation.as_mut().unwrap().dam = "fully".into();
        row.source_generation.as_mut().unwrap().dcm = "unmodified".into();
        row.source_generation.as_mut().unwrap().cpg = "mixed".into();
        let declared = row.generation().unwrap().unwrap();
        assert!(
            declared.dam == MethylationState::FullyModified
                && declared.dcm == MethylationState::Unmodified
                && declared.cpg == MethylationState::HemiOrMixed
        );
        row.source_generation.as_mut().unwrap().dam = "x".repeat(33);
        assert!(
            row.source_generation
                .as_ref()
                .unwrap()
                .validate_raw()
                .is_err()
        );
    }
    #[test]
    fn explicit_generation_follows_reordered_sources_and_changes_context() {
        let mut state = MethodsState::default();
        state.rows[1].document = 1;
        assert!(matches!(
            state.cleavage_policy(),
            collagene_core::cloning::assembly::recipe::CleavagePolicy::AssumedComplete
        ));
        state.mode = Mode::GoldenGate;
        assert!(matches!(
            state.cleavage_policy(),
            collagene_core::cloning::assembly::recipe::CleavagePolicy::PositionalChecked
        ));
        let before = state.preparation_context();
        state.rows[0].source_generation = Some(SourceGenerationInput {
            dam: "fully".into(),
            dcm: "unmodified".into(),
            cpg: "unknown".into(),
        });
        assert!(state.preparation_context() != before);
        let declaration = state.rows[0].source_generation.clone();
        state.swap_rows(0, 1);
        assert!(state.rows[0].source_generation.is_none());
        assert!(state.rows[1].source_generation == declaration);
        assert_eq!(state.rows[1].document, 0);
        state.swap_rows(1, 0);
        assert!(state.rows[0].source_generation == declaration);
    }
}
