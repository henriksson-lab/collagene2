//! Editable preparation fields and stable source identities, independent of UI.
use crate::{assembly_role_ids::Role, golden_gate_pool_form::PoolForm};
use collagene_core::{
    analysis::Primer,
    chemistry::{AnnealOptions, EndRepairMethod, EndTreatment, FivePrime, Terminus, ThreePrime},
    overlap_profiles::{OverlapMethod, ProfileOptions},
    tailing::TailingBase,
    verification::TailSubstrate,
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
pub struct TailInput {
    pub document: usize,
    pub pcr: bool,
    pub enzymes: String,
    pub fragment: String,
    pub blunt: bool,
    pub forward_name: String,
    pub reverse_name: String,
    pub forward: String,
    pub reverse: String,
    pub forward_tail: String,
    pub reverse_tail: String,
    pub base: TailingBase,
    pub treatment: Treatment,
    pub reversed: bool,
}
impl TailInput {
    pub(crate) fn new(pcr: bool, base: TailingBase) -> Self {
        Self {
            document: 0,
            pcr,
            enzymes: "EcoRV".into(),
            fragment: "1".into(),
            blunt: false,
            forward_name: "Forward".into(),
            reverse_name: "Reverse".into(),
            forward: String::new(),
            reverse: String::new(),
            forward_tail: String::new(),
            reverse_tail: String::new(),
            base,
            treatment: Treatment::Keep,
            reversed: false,
        }
    }
    pub(crate) fn recipe(&self) -> Result<TailSubstrate, String> {
        if self.pcr {
            Ok(TailSubstrate::ExactPcr {
                forward: Primer::new(&self.forward_name, &self.forward, &self.forward_tail)
                    .map_err(|e| e.to_string())?,
                reverse: Primer::new(&self.reverse_name, &self.reverse, &self.reverse_tail)
                    .map_err(|e| e.to_string())?,
            })
        } else {
            let enzymes = self
                .enzymes
                .split(',')
                .map(str::trim)
                .filter(|s| !s.is_empty())
                .map(str::to_owned)
                .collect::<Vec<_>>();
            if enzymes.is_empty() {
                return Err("Choose explicit restriction enzymes".into());
            }
            let fragment_index = self
                .fragment
                .parse::<usize>()
                .ok()
                .and_then(|n| n.checked_sub(1))
                .ok_or("Fragment numbers start at 1")?;
            Ok(TailSubstrate::Digest {
                enzymes,
                fragment_index,
                blunt_before: self.blunt,
            })
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
            forward: primer("forward"),
            reverse: primer("reverse"),
        }
    }
    fn typed(&self) -> Result<collagene_core::overlap_extension_pcr::PrimerPair, String> {
        Ok(collagene_core::overlap_extension_pcr::PrimerPair {
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
        })
    }
}
impl Default for OverlapExtensionPcrInput {
    fn default() -> Self {
        let options = collagene_core::overlap_extension::Options::default();
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
    pub(crate) fn request(&self) -> Result<collagene_core::overlap_extension_pcr::Request, String> {
        self.validate_raw()?;
        let number = |text: &str, label: &str| {
            if text.is_empty() || !text.bytes().all(|byte| byte.is_ascii_digit()) {
                return Err(format!("Enter a positive whole number for {label}."));
            }
            text.parse::<usize>()
                .map_err(|_| format!("{label} exceeds the supported range."))
        };
        let request = collagene_core::overlap_extension_pcr::Request {
            left: self.left.typed()?,
            right: self.right.typed()?,
            outer: self.outer.typed()?,
            extension: collagene_core::overlap_extension::Options {
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
#[derive(Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct OverlapInput {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub run_origin: Option<crate::run_origin_model::PinnedOrigin>,
    pub document: usize,
    pub interval_start: String,
    pub interval_end: String,
    pub reverse: bool,
    pub treatment: Treatment,
    pub preparation: OverlapPreparation,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub overlap_extension_pcr: Option<OverlapExtensionPcrInput>,
    #[serde(default)]
    pub catalog_selection: Option<collagene_core::catalog_digest::CatalogSelection>,
    pub left_terminus: Terminus,
    pub right_terminus: Terminus,
    pub enzymes: String,
    pub fragment: String,
    pub repair: Option<EndRepairMethod>,
    pub forward_name: String,
    pub reverse_name: String,
    pub forward: String,
    pub reverse_primer: String,
    pub forward_tail: String,
    pub reverse_tail: String,
    #[serde(
        default = "default_mismatch_max",
        skip_serializing_if = "is_default_mismatch_max"
    )]
    pub mismatch_max_mismatches: String,
    #[serde(
        default = "default_mismatch_exact",
        skip_serializing_if = "is_default_mismatch_exact"
    )]
    pub mismatch_three_prime_exact: String,
    #[serde(
        default = "default_mismatch_hits",
        skip_serializing_if = "is_default_mismatch_hits"
    )]
    pub mismatch_max_hits: String,
    #[serde(
        default = "default_mismatch_work",
        skip_serializing_if = "is_default_mismatch_work"
    )]
    pub mismatch_max_work: String,
    #[serde(
        default = "default_mismatch_pairs",
        skip_serializing_if = "is_default_mismatch_pairs"
    )]
    pub mismatch_max_pairs: String,
    #[serde(
        default = "default_mismatch_products",
        skip_serializing_if = "is_default_mismatch_products"
    )]
    pub mismatch_max_products: String,
    #[serde(
        default = "default_mismatch_bases",
        skip_serializing_if = "is_default_mismatch_bases"
    )]
    pub mismatch_max_product_bases: String,
    #[serde(default)]
    pub right_fusion: String,
    #[serde(default, skip_serializing_if = "default_anneal_options")]
    pub anneal_options: AnnealOptions,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub anneal_candidate: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub anneal_bottom_document: Option<usize>,
    /// Documentary selection context only; physical authority requires replay.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub anneal_source_context: Option<String>,
}
fn default_mismatch_max() -> String {
    collagene_core::mutagenesis::PcrOptions::default()
        .binding
        .max_mismatches
        .to_string()
}
fn default_mismatch_exact() -> String {
    collagene_core::mutagenesis::PcrOptions::default()
        .binding
        .min_three_prime_exact
        .to_string()
}
fn is_default_mismatch_max(value: &String) -> bool {
    *value == default_mismatch_max()
}
fn is_default_mismatch_exact(value: &String) -> bool {
    *value == default_mismatch_exact()
}
fn default_mismatch_hits() -> String {
    collagene_core::mutagenesis::PcrOptions::default()
        .binding
        .max_hits
        .to_string()
}
fn is_default_mismatch_hits(value: &String) -> bool {
    *value == default_mismatch_hits()
}
fn default_mismatch_work() -> String {
    collagene_core::mutagenesis::PcrOptions::default()
        .binding
        .max_work
        .to_string()
}
fn is_default_mismatch_work(value: &String) -> bool {
    *value == default_mismatch_work()
}
fn default_mismatch_pairs() -> String {
    collagene_core::mutagenesis::PcrOptions::default()
        .max_pair_checks
        .to_string()
}
fn is_default_mismatch_pairs(value: &String) -> bool {
    *value == default_mismatch_pairs()
}
fn default_mismatch_products() -> String {
    collagene_core::mutagenesis::PcrOptions::default()
        .max_products
        .to_string()
}
fn is_default_mismatch_products(value: &String) -> bool {
    *value == default_mismatch_products()
}
fn default_mismatch_bases() -> String {
    collagene_core::mutagenesis::PcrOptions::default()
        .max_product_bases
        .to_string()
}
fn is_default_mismatch_bases(value: &String) -> bool {
    *value == default_mismatch_bases()
}
fn default_anneal_options(options: &AnnealOptions) -> bool {
    *options == AnnealOptions::default()
}
impl Default for OverlapInput {
    fn default() -> Self {
        Self {
            run_origin: None,
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
    pub(crate) fn is_pcr(&self) -> bool {
        matches!(
            self.preparation,
            OverlapPreparation::ExactPcr
                | OverlapPreparation::PcrCatalogDigest
                | OverlapPreparation::MismatchPcr
                | OverlapPreparation::MismatchPcrCatalogDigest
        )
    }
    /// Select the complete PCR product without replacing the declared PCR policy.
    pub(crate) fn use_complete_pcr_amplicon(&mut self) -> Result<(), String> {
        if !self.is_pcr() {
            return Err(
                "Choose an explicit PCR preparation before using its complete amplicon.".into(),
            );
        }
        if self.preparation == OverlapPreparation::PcrCatalogDigest {
            self.preparation = OverlapPreparation::ExactPcr;
        } else if self.preparation == OverlapPreparation::MismatchPcrCatalogDigest {
            self.preparation = OverlapPreparation::MismatchPcr;
        }
        self.fragment = "1".into();
        Ok(())
    }
    pub(crate) fn mismatch_options(
        &self,
    ) -> Result<collagene_core::mutagenesis::PcrOptions, String> {
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
        let mut options = collagene_core::mutagenesis::PcrOptions::default();
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
            OverlapPreparation::OverlapExtensionPcr => format!(
                "Two source PCRs → reciprocal extension ({} bp homology) → outer-primer PCR",
                self.overlap_extension_pcr
                    .as_ref()
                    .map(|request| request.overlap.as_str())
                    .unwrap_or("enter")
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
#[derive(Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct MethodsState {
    pub mode: Mode,
    pub vector: TailInput,
    pub insert: TailInput,
    pub circular: bool,
    pub allow_nicked: bool,
    #[serde(deserialize_with = "deserialize_method_rows")]
    pub rows: Vec<OverlapInput>,
    pub overlaps: String,
    pub method: OverlapMethod,
    pub scan_repeats: bool,
    pub comparisons: String,
    pub occurrences: String,
    #[serde(default)]
    pub golden_gate: GoldenGateForm,
    /// Editable pool/reaction settings and explicit candidate intent only;
    /// inspections, plans and scientific authority remain transient.
    /// Missing in older drafts: initialize using their preparation topology and row count.
    #[serde(default)]
    pub golden_gate_pool: Option<PoolForm>,
    #[serde(default)]
    pub multisite: Option<crate::multisite_form::Form>,
    /// Stable source-choice and bin identities. Older ordinary drafts keep the
    /// existing ordinal identity contract; explicit identities are never repinned.
    #[serde(default, deserialize_with = "deserialize_source_roles")]
    pub source_roles: Option<Vec<Role>>,
}
#[derive(Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct GoldenGateForm {
    pub enzyme: collagene_core::golden_gate::GoldenGateEnzyme,
    pub clamp: String,
    pub spacer: String,
    pub forward_length: String,
    pub reverse_length: String,
}
impl Default for GoldenGateForm {
    fn default() -> Self {
        Self {
            enzyme: collagene_core::golden_gate::GoldenGateEnzyme::BsaI,
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
            vector: TailInput::new(false, TailingBase::T),
            insert: TailInput::new(true, TailingBase::A),
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
        }
    }
}
impl MethodsState {
    pub(crate) fn bottom_input_id(&self, row: usize) -> Result<String, String> {
        if self
            .rows
            .get(row)
            .is_none_or(|input| input.preparation != OverlapPreparation::AnnealedOligos)
        {
            return Err("Choose an active paired-oligo preparation row.".into());
        }
        let roles = self.role_ids(self.rows.len())?;
        Ok(format!(
            "anneal-bottom-{}",
            blake3::hash(roles[row].input_id.as_bytes()).to_hex()
        ))
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
        if self.mode == Mode::Tailing {
            return Ok(vec![self.vector.document, self.insert.document]);
        }
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
            return crate::multisite_form::ordinary_roles(count);
        }
        crate::assembly_role_ids::resolve(&self.source_roles, count)
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
            // Preparation/roles travel with rows; sync_rows rebinds the legacy
            // document/reverse fields to the new slot positions afterwards.
            std::mem::swap(&mut self.vector.base, &mut self.insert.base);
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
    pub(crate) fn profile(&self) -> Result<(Vec<usize>, ProfileOptions), String> {
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
    fn legacy_row_json_defaults_and_mismatch_policy_bounds_are_explicit() {
        let ordinary = OverlapInput::default();
        let json = serde_json::to_value(&ordinary).unwrap();
        assert!(json.get("mismatch_max_mismatches").is_none());
        assert!(json.get("mismatch_three_prime_exact").is_none());
        let mut row: OverlapInput = serde_json::from_value(json).unwrap();
        row.forward = "ACGT".into();
        row.reverse_primer = "TGCA".into();
        let native = collagene_core::mutagenesis::PcrOptions::default();
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
    fn partial_six_oligo_input_roundtrips_without_changing_legacy_row_json() {
        let legacy = serde_json::to_value(OverlapInput::default()).unwrap();
        assert!(legacy.get("overlap_extension_pcr").is_none());
        let mut row: OverlapInput = serde_json::from_value(legacy).unwrap();
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
        assert!(row.use_complete_pcr_amplicon().is_err());
        assert!(row == before);
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
