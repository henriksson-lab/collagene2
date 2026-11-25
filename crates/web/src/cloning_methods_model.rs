//! Editable preparation fields and stable source identities, independent of UI.
use crate::{assembly_role_ids::Role, golden_gate_pool_form::PoolForm};
use collagene_core::{
    analysis::Primer,
    chemistry::{EndRepairMethod, EndTreatment, FivePrime, Terminus, ThreePrime},
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
}
#[derive(Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct OverlapInput {
    pub document: usize,
    pub interval_start: String,
    pub interval_end: String,
    pub reverse: bool,
    pub treatment: Treatment,
    pub preparation: OverlapPreparation,
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
    #[serde(default)]
    pub right_fusion: String,
}
impl Default for OverlapInput {
    fn default() -> Self {
        Self {
            document: 0,
            interval_start: "1".into(),
            interval_end: String::new(),
            reverse: false,
            treatment: Treatment::Keep,
            preparation: OverlapPreparation::WholeLinear,
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
            right_fusion: String::new(),
        }
    }
}
impl OverlapInput {
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
            OverlapPreparation::CatalogDigest | OverlapPreparation::PcrCatalogDigest => format!(
                "{}: {}; fragment {}",
                if self.preparation == OverlapPreparation::PcrCatalogDigest {
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
            OverlapPreparation::ExactPcr => format!(
                "Exact PCR: annealing {} / {} bp; 5′ extensions {} / {} bp",
                self.forward.len(),
                self.reverse_primer.len(),
                self.forward_tail.len(),
                self.reverse_tail.len()
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
            || self.role_ids(self.rows.len()).is_err()
            || self.golden_gate_pool.as_ref().is_some_and(|form| {
                !form.choices.is_empty() && form.choices.len() != self.rows.len()
            })
        {
            return;
        }
        self.rows.swap(a, b);
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
