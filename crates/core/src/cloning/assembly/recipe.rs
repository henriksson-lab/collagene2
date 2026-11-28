//! Compile explicit preparation rows into existing typed operations. Fresh
//! reconstruction remains the verifier's job; this module never slices a source
//! document and calls that physical preparation, or searches toward a target.
use crate::{
    cloning::overlap::profiles::{self as overlap_profiles, ProfileOptions, ProfileReport},
    encoding::{bounded_json::hash, required_option},
    planning::plan::{DesignPlan, Operation, OrientedSource, PlanInput, Source},
    sequence::analysis::Primer,
    verification::{self, evidence::VerificationReport},
};
use collagene_model::{Alphabet, Error, Molecule, Topology};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};
pub const ENGINE_VERSION: &str = "overlap-recipe-3";
pub const MAX_ROWS: usize = 128;
const MAX_INPUT_BASES: usize = 8 << 20;
pub(crate) const MAX_ROLE_BASES: usize = 16 << 20;
const MAX_PREPARATION_WORK: usize = 100_000_000;
pub(crate) const MAX_STEPS_PER_ROW: usize = 7;
macro_rules! dto{($($t:tt)*)=>{#[derive(Clone,Debug,Serialize,Deserialize)]#[cfg_attr(feature="json-schema",derive(schemars::JsonSchema))]#[serde(deny_unknown_fields)]$($t)*};}
// Explicit finishing semantics shared by ordinary recipes and source libraries.
// A planned host rescue remains documentary; it does not fill or seal strands.
dto! {pub enum OverlapFinishing {
    ReferenceOnly,
    HomologyAssembly {policy:crate::cloning::overlap::assembly::Options},
    HomologyAnnealing {policy:crate::cloning::overlap::annealing::Options},
}}
dto! {pub enum Preparation {WholeLinear,Digest{enzymes:Vec<String>,fragment_index:usize,repair:Option<crate::chemistry::EndRepairMethod>},ExactPcr{forward:Primer,reverse:Primer,template_strand:crate::chemistry::Strand},DeclaredDuplex{left:crate::chemistry::Terminus,right:crate::chemistry::Terminus},CatalogDigest{selection:crate::cloning::restriction::catalog::CatalogSelection,#[serde(deserialize_with="crate::encoding::required_option")]#[cfg_attr(feature="json-schema",schemars(required))]conditions:Option<crate::cloning::restriction::conditions::ReactionConditions>,fragment_index:usize,repair:Option<crate::chemistry::EndRepairMethod>},PcrCatalogDigest{forward:Primer,reverse:Primer,template_strand:crate::chemistry::Strand,selection:crate::cloning::restriction::catalog::CatalogSelection,#[serde(deserialize_with="crate::encoding::required_option")]#[cfg_attr(feature="json-schema",schemars(required))]conditions:Option<crate::cloning::restriction::conditions::ReactionConditions>,fragment_index:usize,repair:Option<crate::chemistry::EndRepairMethod>},AnnealedOligos{bottom_input_id:String,options:crate::chemistry::AnnealOptions,candidate_index:usize},FromRun{origin:crate::material::capsule::Capsule,preparation:RunPreparation},MismatchPcr{forward:Primer,reverse:Primer,options:crate::cloning::pcr::mutagenesis::PcrOptions},MismatchPcrCatalogDigest{forward:Primer,reverse:Primer,options:crate::cloning::pcr::mutagenesis::PcrOptions,selection:crate::cloning::restriction::catalog::CatalogSelection,#[serde(deserialize_with="crate::encoding::required_option")]#[cfg_attr(feature="json-schema",schemars(required))]conditions:Option<crate::cloning::restriction::conditions::ReactionConditions>,fragment_index:usize,repair:Option<crate::chemistry::EndRepairMethod>},OverlapExtensionPcr{request:crate::cloning::overlap::extension::pcr::Request},OverlapExtensionPcrCatalogDigest{request:crate::cloning::overlap::extension::pcr::Request,selection:crate::cloning::restriction::catalog::CatalogSelection,#[serde(deserialize_with="crate::encoding::required_option")]#[cfg_attr(feature="json-schema",schemars(required))]conditions:Option<crate::cloning::restriction::conditions::ReactionConditions>,fragment_index:usize,repair:Option<crate::chemistry::EndRepairMethod>},FromDerived{origin:crate::material::derivation::Capsule,preparation:RunPreparation}}}
impl Preparation {
    /// Literal catalog pin used by a preparation, including a flat saved-origin
    /// child. This is structural lookup only, never replay or cleavage authority.
    pub fn catalog_selection(
        &self,
    ) -> Option<&crate::cloning::restriction::catalog::CatalogSelection> {
        match self {
            Preparation::CatalogDigest { selection, .. }
            | Preparation::PcrCatalogDigest { selection, .. }
            | Preparation::MismatchPcrCatalogDigest { selection, .. }
            | Preparation::OverlapExtensionPcrCatalogDigest { selection, .. } => Some(selection),
            Preparation::FromRun { preparation, .. }
            | Preparation::FromDerived { preparation, .. } => match preparation {
                RunPreparation::CatalogDigest { selection, .. }
                | RunPreparation::PcrCatalogDigest { selection, .. }
                | RunPreparation::MismatchPcrCatalogDigest { selection, .. }
                | RunPreparation::OverlapExtensionPcrCatalogDigest { selection, .. } => {
                    Some(selection)
                }
                _ => None,
            },
            _ => None,
        }
    }
    /// Explicit documentary choices; absence does not imply supplier compatibility.
    pub fn reaction_conditions(
        &self,
    ) -> Option<&crate::cloning::restriction::conditions::ReactionConditions> {
        match self {
            Self::CatalogDigest { conditions, .. }
            | Self::PcrCatalogDigest { conditions, .. }
            | Self::MismatchPcrCatalogDigest { conditions, .. }
            | Self::OverlapExtensionPcrCatalogDigest { conditions, .. } => conditions.as_ref(),
            Self::FromRun { preparation, .. } | Self::FromDerived { preparation, .. } => {
                match preparation {
                    RunPreparation::CatalogDigest { conditions, .. }
                    | RunPreparation::PcrCatalogDigest { conditions, .. }
                    | RunPreparation::MismatchPcrCatalogDigest { conditions, .. }
                    | RunPreparation::OverlapExtensionPcrCatalogDigest { conditions, .. } => {
                        conditions.as_ref()
                    }
                    _ => None,
                }
            }
            _ => None,
        }
    }
    /// Exact flat operation count, for host/preflight limits before any model copy.
    /// It counts the declared graph only and does not establish reaction completion.
    pub(crate) fn operation_count(&self) -> usize {
        let repair =
            |repair: &Option<crate::chemistry::EndRepairMethod>| usize::from(repair.is_some());
        match self {
            Preparation::WholeLinear => 0,
            Preparation::ExactPcr { .. }
            | Preparation::MismatchPcr { .. }
            | Preparation::DeclaredDuplex { .. }
            | Preparation::AnnealedOligos { .. } => 1,
            Preparation::Digest { repair: r, .. }
            | Preparation::CatalogDigest { repair: r, .. } => 1 + repair(r),
            Preparation::PcrCatalogDigest { repair: r, .. }
            | Preparation::MismatchPcrCatalogDigest { repair: r, .. } => 2 + repair(r),
            Preparation::OverlapExtensionPcr { .. } => 4,
            Preparation::OverlapExtensionPcrCatalogDigest { repair: r, .. } => 5 + repair(r),
            Preparation::FromRun { preparation, .. }
            | Preparation::FromDerived { preparation, .. } => {
                1 + match preparation {
                    RunPreparation::WholeLinear => 0,
                    RunPreparation::ExactPcr { .. } | RunPreparation::MismatchPcr { .. } => 1,
                    RunPreparation::Digest { repair: r, .. }
                    | RunPreparation::CatalogDigest { repair: r, .. } => 1 + repair(r),
                    RunPreparation::PcrCatalogDigest { repair: r, .. }
                    | RunPreparation::MismatchPcrCatalogDigest { repair: r, .. } => 2 + repair(r),
                    RunPreparation::OverlapExtensionPcr { .. } => 4,
                    RunPreparation::OverlapExtensionPcrCatalogDigest { repair: r, .. } => {
                        5 + repair(r)
                    }
                }
            }
        }
    }
}
dto! {pub enum RunPreparation {WholeLinear,Digest{enzymes:Vec<String>,fragment_index:usize,repair:Option<crate::chemistry::EndRepairMethod>},CatalogDigest{selection:crate::cloning::restriction::catalog::CatalogSelection,#[serde(deserialize_with="crate::encoding::required_option")]#[cfg_attr(feature="json-schema",schemars(required))]conditions:Option<crate::cloning::restriction::conditions::ReactionConditions>,fragment_index:usize,repair:Option<crate::chemistry::EndRepairMethod>},ExactPcr{forward:Primer,reverse:Primer,template_strand:crate::chemistry::Strand},PcrCatalogDigest{forward:Primer,reverse:Primer,template_strand:crate::chemistry::Strand,selection:crate::cloning::restriction::catalog::CatalogSelection,#[serde(deserialize_with="crate::encoding::required_option")]#[cfg_attr(feature="json-schema",schemars(required))]conditions:Option<crate::cloning::restriction::conditions::ReactionConditions>,fragment_index:usize,repair:Option<crate::chemistry::EndRepairMethod>},MismatchPcr{forward:Primer,reverse:Primer,options:crate::cloning::pcr::mutagenesis::PcrOptions},MismatchPcrCatalogDigest{forward:Primer,reverse:Primer,options:crate::cloning::pcr::mutagenesis::PcrOptions,selection:crate::cloning::restriction::catalog::CatalogSelection,#[serde(deserialize_with="crate::encoding::required_option")]#[cfg_attr(feature="json-schema",schemars(required))]conditions:Option<crate::cloning::restriction::conditions::ReactionConditions>,fragment_index:usize,repair:Option<crate::chemistry::EndRepairMethod>},OverlapExtensionPcr{request:crate::cloning::overlap::extension::pcr::Request},OverlapExtensionPcrCatalogDigest{request:crate::cloning::overlap::extension::pcr::Request,selection:crate::cloning::restriction::catalog::CatalogSelection,#[serde(deserialize_with="crate::encoding::required_option")]#[cfg_attr(feature="json-schema",schemars(required))]conditions:Option<crate::cloning::restriction::conditions::ReactionConditions>,fragment_index:usize,repair:Option<crate::chemistry::EndRepairMethod>}}}
/// Actual cleavage model selected by the compiling workflow, never inferred from metadata.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum CleavagePolicy {
    AssumedComplete,
    PositionalChecked,
}
fn digest_operation(
    template: Source,
    selection: crate::cloning::restriction::catalog::CatalogSelection,
    fragment_index: usize,
    policy: CleavagePolicy,
    conditions: Option<crate::cloning::restriction::conditions::ReactionConditions>,
) -> Operation {
    match policy {
        CleavagePolicy::AssumedComplete => Operation::CatalogDigest {
            template,
            selection,
            conditions,
            fragment_index,
        },
        CleavagePolicy::PositionalChecked => Operation::CheckedDigest {
            template,
            selection,
            conditions,
            options: crate::chemistry::tracked_duplex::checked_inspection_options(),
            fragment_index,
        },
    }
}
fn starter_selection(
    enzymes: &[String],
) -> Result<crate::cloning::restriction::catalog::CatalogSelection, Error> {
    let catalog = crate::cloning::restriction::catalog::builtin_catalog()?;
    let ids = enzymes
        .iter()
        .map(|name| {
            catalog
                .edition()
                .enzymes
                .iter()
                .find(|e| e.name.eq_ignore_ascii_case(name))
                .map(|e| e.id.clone())
                .ok_or_else(|| bad("starter digest enzyme has no exact current catalog identity"))
        })
        .collect::<Result<Vec<_>, _>>()?;
    crate::cloning::restriction::catalog::select(&catalog, &ids)
}

impl RunPreparation {
    /// Bounded conversion to the existing preparation vocabulary; no recursion,
    /// chemistry restoration, or scientific computation.
    pub(crate) fn preparation(&self) -> Result<Preparation, Error> {
        match self {
            Self::CatalogDigest {
                selection,
                conditions: Some(conditions),
                ..
            }
            | Self::PcrCatalogDigest {
                selection,
                conditions: Some(conditions),
                ..
            }
            | Self::MismatchPcrCatalogDigest {
                selection,
                conditions: Some(conditions),
                ..
            }
            | Self::OverlapExtensionPcrCatalogDigest {
                selection,
                conditions: Some(conditions),
                ..
            } => conditions.validate(selection)?,
            _ => {}
        }
        match self {
            Self::OverlapExtensionPcrCatalogDigest {
                request,
                selection,
                conditions,
                fragment_index,
                repair,
            } => {
                overlap_extension_catalog_shape(request, selection, *fragment_index)?;
                Ok(Preparation::OverlapExtensionPcrCatalogDigest {
                    request: request.clone(),
                    selection: selection.clone(),
                    conditions: conditions.clone(),
                    fragment_index: *fragment_index,
                    repair: *repair,
                })
            }
            Self::OverlapExtensionPcr { request } => {
                request.validate_shape()?;
                Ok(Preparation::OverlapExtensionPcr {
                    request: request.clone(),
                })
            }
            Self::WholeLinear => Ok(Preparation::WholeLinear),
            Self::Digest {
                enzymes,
                fragment_index,
                repair,
            } => {
                if enzymes.is_empty()
                    || enzymes.len() > 32
                    || enzymes.iter().any(|e| e.is_empty() || e.len() > 128)
                    || *fragment_index > 100_000
                {
                    return Err(Error::TooLarge);
                }
                Ok(Preparation::Digest {
                    enzymes: enzymes.clone(),
                    fragment_index: *fragment_index,
                    repair: *repair,
                })
            }
            Self::CatalogDigest {
                selection,
                conditions,
                fragment_index,
                repair,
            } => {
                selection.validate()?;
                if *fragment_index > 100_000 {
                    return Err(Error::TooLarge);
                }
                Ok(Preparation::CatalogDigest {
                    selection: selection.clone(),
                    conditions: conditions.clone(),
                    fragment_index: *fragment_index,
                    repair: *repair,
                })
            }
            Self::MismatchPcrCatalogDigest {
                forward,
                reverse,
                options,
                selection,
                conditions,
                fragment_index,
                repair,
            } => {
                mismatch_pcr_catalog_shape(forward, reverse, options, selection, *fragment_index)?;
                Ok(Preparation::MismatchPcrCatalogDigest {
                    forward: forward.clone(),
                    reverse: reverse.clone(),
                    options: options.clone(),
                    selection: selection.clone(),
                    conditions: conditions.clone(),
                    fragment_index: *fragment_index,
                    repair: *repair,
                })
            }
            Self::MismatchPcr {
                forward,
                reverse,
                options,
            } => {
                mismatch_pcr_shape(forward, reverse, options)?;
                Ok(Preparation::MismatchPcr {
                    forward: forward.clone(),
                    reverse: reverse.clone(),
                    options: options.clone(),
                })
            }
            Self::ExactPcr {
                forward,
                reverse,
                template_strand,
            } => {
                primer(forward)?;
                primer(reverse)?;
                Ok(Preparation::ExactPcr {
                    forward: forward.clone(),
                    reverse: reverse.clone(),
                    template_strand: *template_strand,
                })
            }
            Self::PcrCatalogDigest {
                forward,
                reverse,
                template_strand,
                selection,
                conditions,
                fragment_index,
                repair,
            } => {
                primer(forward)?;
                primer(reverse)?;
                selection.validate()?;
                if *fragment_index > 100_000 {
                    return Err(Error::TooLarge);
                }
                Ok(Preparation::PcrCatalogDigest {
                    forward: forward.clone(),
                    reverse: reverse.clone(),
                    template_strand: *template_strand,
                    selection: selection.clone(),
                    conditions: conditions.clone(),
                    fragment_index: *fragment_index,
                    repair: *repair,
                })
            }
        }
    }
}
pub(crate) fn from_run_shape(
    origin: &crate::material::capsule::Capsule,
    preparation: &RunPreparation,
) -> Result<(), Error> {
    origin.validate_documentary()?;
    preparation.preparation()?;
    Ok(())
}

pub(crate) fn from_derived_shape(
    origin: &crate::material::derivation::Capsule,
    preparation: &RunPreparation,
) -> Result<(), Error> {
    origin.validate_documentary()?;
    preparation.preparation()?;
    Ok(())
}
dto! {pub struct FragmentSpec {pub role:String,pub input_id:String,pub preparation:Preparation,pub reverse:bool,#[serde(deserialize_with="required_option")]#[cfg_attr(feature="json-schema",schemars(required))]pub generation:Option<crate::chemistry::tracked_duplex::GenerationAssumptions>}}
dto! {
    pub struct Recipe {
    pub name:String,
    pub fragments:Vec<FragmentSpec>,
    pub overlaps:Vec<usize>,
    pub circular:bool,
    pub options:ProfileOptions,
    #[serde(deserialize_with="crate::encoding::required_option")]
    #[cfg_attr(feature="json-schema",schemars(required))]
    pub standard:Option<crate::cloning::assembly::standards::declaration::Declaration>,
    }
}
dto! {pub struct PreparedRow {pub row_index:usize,pub role:String,pub input_id:String,pub input_revision:String,pub preparation_step:Option<usize>,pub repair_step:Option<usize>,pub prepared_source:Source,pub reverse:bool}}
dto! {pub struct Compiled {pub plan:DesignPlan,pub rows:Vec<PreparedRow>,pub assembly_step:usize}}
dto! {pub struct Preview {pub context_revision:String,pub compiled:Compiled,pub verification:VerificationReport,
/// Original-row-indexed prepared documents, before optional orientation.
/// Failed or blocked preparations retain explicit empty slots.
pub prepared:Vec<Option<Molecule>>,pub profile:Option<ProfileReport>,pub preparation_complete:bool}}
fn bad(s: &str) -> Error {
    Error::Format(format!("overlap recipe: {s}"))
}
fn hex(s: &str) -> bool {
    s.len() == 64
        && s.bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}
fn token(s: &str, cap: usize) -> Result<(), Error> {
    if s.is_empty() || s.len() > cap || s.chars().any(char::is_control) {
        Err(bad("invalid bounded source/role/name"))
    } else {
        Ok(())
    }
}
/// Stable secondary-input identity shared by ordinary and library paired-oligo recipes.
pub fn anneal_bottom_input_id(primary: &str) -> Result<String, Error> {
    if primary.is_empty()
        || primary.len() > 1000
        || primary.trim() != primary
        || primary.chars().any(char::is_control)
    {
        return Err(bad("invalid primary input identity for paired oligos"));
    }
    Ok(format!(
        "anneal-bottom-{}",
        blake3::hash(primary.as_bytes()).to_hex()
    ))
}
/// Borrowed bounds only: candidate availability is reconstructed by the verifier.
pub(crate) fn anneal_shape(
    bottom: &str,
    options: &crate::chemistry::AnnealOptions,
    candidate: usize,
) -> Result<(), Error> {
    token(bottom, 128)?;
    if options.min_paired_bases == 0
        || options.min_paired_bases > crate::chemistry::MAX_ANNEAL_OLIGO_BASES
        || options.max_comparisons == 0
        || options.max_comparisons > crate::chemistry::MAX_ANNEAL_COMPARISONS
        || options.max_candidates == 0
        || options.max_candidates > crate::chemistry::MAX_ANNEAL_CANDIDATES
        || candidate >= options.max_candidates
    {
        return Err(bad("invalid annealing limits or selected candidate"));
    }
    Ok(())
}
pub(crate) fn anneal_sources(top: &Molecule, bottom: &Molecule) -> Result<(), Error> {
    for source in [top, bottom] {
        if source.sequence.is_empty()
            || source.sequence.len() > crate::chemistry::MAX_ANNEAL_OLIGO_BASES
        {
            return Err(Error::TooLarge);
        }
        verification::terminal::declared_oligo(source)?;
        if !source
            .sequence
            .as_str()
            .bytes()
            .all(|b| b"ACGT".contains(&b))
        {
            return Err(bad("annealing requires concrete DNA oligos"));
        }
    }
    Ok(())
}
/// Check every original binding, including the secondary strand of a paired preparation.
/// This is structural validation; it never enumerates annealing placements.
pub(crate) fn validate_prepared_inputs(
    inputs: &[PlanInput],
    steps: &[Operation],
) -> Result<(), Error> {
    if inputs.is_empty() || inputs.len() > MAX_ROWS {
        return Err(Error::TooLarge);
    }
    crate::material::capsule::material_capsules_preflight(steps.iter().filter_map(
        |step| match step {
            Operation::ReconstructRunProduct { origin, .. } => {
                Some(crate::material::capsule::MaterialOriginRef::Run(origin))
            }
            Operation::ReconstructMaterial { origin, .. } => {
                Some(crate::material::capsule::MaterialOriginRef::Derived(origin))
            }
            _ => None,
        },
    ))?;
    let mut ids = BTreeSet::new();
    for input in inputs {
        token(&input.id, 128)?;
        collagene_model::bounds::model_preflight(&[&input.molecule])?;
        if !ids.insert(input.id.as_str()) || input.revision != input.molecule.revision() {
            return Err(bad("original input binding changed"));
        }
    }
    for step in steps {
        if let Operation::Pcr {
            forward, reverse, ..
        } = step
        {
            primer(forward)?;
            primer(reverse)?;
        }
        if let Operation::OverlapExtension { options, .. } = step {
            crate::cloning::overlap::extension::pcr::extension_shape(options)?;
        }
        if let Operation::PcrWithMismatches {
            forward,
            reverse,
            options,
            ..
        } = step
        {
            mismatch_pcr_shape(forward, reverse, options)?;
        }
        if let Operation::ReconstructRunProduct { template, origin } = step {
            origin.validate_documentary()?;
            let Source::Input(id) = template else {
                return Err(bad(
                    "material reconstruction must consume the exact saved original input",
                ));
            };
            let input = inputs
                .iter()
                .find(|i| &i.id == id)
                .ok_or_else(|| bad("missing saved material source"))?;
            if input.revision != origin.selection.event.revision {
                return Err(bad("saved material source revision changed"));
            }
        }

        if let Operation::ReconstructMaterial { template, origin } = step {
            origin.validate_documentary()?;
            let Source::Input(id) = template else {
                return Err(bad(
                    "material derivation reconstruction must consume its exact original input",
                ));
            };
            let input = inputs
                .iter()
                .find(|i| &i.id == id)
                .ok_or_else(|| bad("missing saved material derivation source"))?;
            if input.revision != origin.product.revision {
                return Err(bad("saved material derivation source revision changed"));
            }
        }
        if let Operation::AnnealOligos {
            top,
            bottom,
            options,
            candidate_index,
        } = step
        {
            let (Source::Input(top), Source::Input(bottom)) = (top, bottom) else {
                return Err(bad("annealing requires original input oligos"));
            };
            anneal_shape(bottom, options, *candidate_index)?;
            if top == bottom {
                return Err(bad(
                    "annealing requires two distinct original input identities",
                ));
            }
            let resolve = |id: &str| {
                inputs
                    .iter()
                    .find(|i| i.id == id)
                    .map(|i| &i.molecule)
                    .ok_or_else(|| bad("missing annealing original"))
            };
            anneal_sources(resolve(top)?, resolve(bottom)?)?;
        }
    }
    Ok(())
}
pub(crate) fn validate_annealed_row(row: &PreparedRow, steps: &[Operation]) -> Result<(), Error> {
    if let Some(i) = row.preparation_step {
        if let Some(Operation::AnnealOligos { top, .. }) = steps.get(i) {
            if !matches!(top, Source::Input(id) if id == &row.input_id)
                || row.repair_step.is_some()
                || !matches!(row.prepared_source, Source::Step(j) if j == i)
            {
                return Err(bad(
                    "annealed row lost its primary original or prepared step",
                ));
            }
        }
    }
    Ok(())
}
fn preparation_template(op: &Operation) -> Option<&Source> {
    match op {
        Operation::Pcr { template, .. }
        | Operation::PcrWithMismatches { template, .. }
        | Operation::Digest { template, .. }
        | Operation::CatalogDigest { template, .. }
        | Operation::CheckedDigest { template, .. }
        | Operation::DeclareSourceModifications { template, .. }
        | Operation::DeclareSourceTermini { template, .. }
        | Operation::ReconstructRunProduct { template, .. }
        | Operation::ReconstructMaterial { template, .. } => Some(template),
        _ => None,
    }
}
fn same_source(a: &Source, b: &Source) -> bool {
    match (a, b) {
        (Source::Input(a), Source::Input(b)) => a == b,
        (Source::Step(a), Source::Step(b)) => a == b,
        _ => false,
    }
}
pub(crate) fn prepared_chain(
    plan: &DesignPlan,
    row: &PreparedRow,
    index: usize,
    cursor: &mut usize,
    assembly: usize,
    owners: &mut [Option<usize>],
) -> Result<Vec<usize>, Error> {
    if row.row_index != index
        || row.role.trim().is_empty()
        || row.role.len() > 128
        || row.input_id.is_empty()
        || row.input_id.len() > 128
        || !plan
            .inputs
            .iter()
            .any(|i| i.id == row.input_id && i.revision == row.input_revision)
    {
        return Err(bad("row identity/source binding changed"));
    }
    let input = Source::Input(row.input_id.clone());
    let mut steps = Vec::new();
    if let Some(first) = row.preparation_step {
        if first != *cursor || first >= assembly {
            return Err(bad("preparation chain is not in compiler order"));
        }
        let op = &plan.steps[first];
        let initial_bound = match op {
            Operation::AnnealOligos {
                top,
                bottom: Source::Input(bottom),
                ..
            } => {
                same_source(top, &input)
                    && bottom != &row.input_id
                    && plan.inputs.iter().any(|i| &i.id == bottom)
            }
            Operation::ReconstructRunProduct { .. }
            | Operation::ReconstructMaterial { .. }
            | Operation::Pcr { .. }
            | Operation::PcrWithMismatches { .. }
            | Operation::Digest { .. }
            | Operation::CatalogDigest { .. }
            | Operation::CheckedDigest { .. }
            | Operation::DeclareSourceModifications { .. }
            | Operation::DeclareSourceTermini { .. } => {
                preparation_template(op).is_some_and(|s| same_source(s, &input))
            }
            _ => false,
        };
        if !initial_bound {
            return Err(bad(
                "initial preparation must consume its own original input",
            ));
        }
        steps.push(first);
        *cursor += 1;
        if matches!(
            op,
            Operation::ReconstructRunProduct { .. }
                | Operation::ReconstructMaterial { .. }
                | Operation::DeclareSourceModifications { .. }
        ) && *cursor < assembly
            && matches!(plan.steps.get(*cursor),Some(Operation::Pcr{template:Source::Step(n),..}|Operation::PcrWithMismatches{template:Source::Step(n),..}|Operation::Digest{template:Source::Step(n),..}|Operation::CatalogDigest{template:Source::Step(n),..}|Operation::CheckedDigest{template:Source::Step(n),..}) if *n==first)
        {
            steps.push(*cursor);
            *cursor += 1;
        }
        let preceding = *steps.last().ok_or_else(|| bad("empty preparation chain"))?;
        if matches!(plan.steps.get(preceding), Some(Operation::Pcr { .. }))
            && matches!(
                plan.steps.get(preceding + 2),
                Some(Operation::OverlapExtension { .. })
            )
        {
            let second = preceding + 1;
            let extension = preceding + 2;
            let outer = preceding + 3;
            let first_source = preparation_template(&plan.steps[preceding])
                .ok_or_else(|| bad("overlap extension lost its first PCR source"))?;
            if outer >= assembly
                || *cursor != second
                || !matches!(plan.steps.get(second), Some(Operation::Pcr { template, .. }) if same_source(template, first_source))
                || !matches!(plan.steps.get(extension), Some(Operation::OverlapExtension { left, right, .. })
                    if !left.reverse && !right.reverse
                    && same_source(&left.source, &Source::Step(preceding))
                    && same_source(&right.source, &Source::Step(second)))
                || !matches!(plan.steps.get(outer), Some(Operation::Pcr { template, .. }) if same_source(template, &Source::Step(extension)))
            {
                return Err(bad(
                    "overlap-extension PCR lost its original-source branches, selected pair or outer PCR",
                ));
            }
            steps.extend([second, extension, outer]);
            *cursor = outer + 1;
        }
        let preceding = *steps.last().ok_or_else(|| bad("empty preparation chain"))?;
        if matches!(
            plan.steps[preceding],
            Operation::Pcr { .. } | Operation::PcrWithMismatches { .. }
        ) && *cursor < assembly
            && matches!(plan.steps.get(*cursor),Some(Operation::CatalogDigest{template:Source::Step(n),..}|Operation::CheckedDigest{template:Source::Step(n),..})if *n==preceding)
        {
            steps.push(*cursor);
            *cursor += 1;
        }
        if let Some(repair) = row.repair_step {
            let last = *steps
                .last()
                .ok_or_else(|| bad("repair has no preparation"))?;
            if repair != *cursor
                || repair >= assembly
                || !matches!(
                    plan.steps[last],
                    Operation::Digest { .. }
                        | Operation::CatalogDigest { .. }
                        | Operation::CheckedDigest { .. }
                )
                || !matches!(plan.steps.get(repair),Some(Operation::RepairEnds{template:Source::Step(n),..})if *n==last)
            {
                return Err(bad("repair lost its exact digest source"));
            }
            steps.push(repair);
            *cursor += 1;
        }
        if !same_source(
            &row.prepared_source,
            &Source::Step(*steps.last().ok_or_else(|| bad("empty preparation chain"))?),
        ) {
            return Err(bad("prepared source is not the chain's final step"));
        }
    } else if row.repair_step.is_some() || !same_source(&row.prepared_source, &input) {
        return Err(bad("unprepared row must retain its own original reference"));
    }
    if steps.len() > MAX_STEPS_PER_ROW {
        return Err(Error::TooLarge);
    }
    for &i in &steps {
        if owners
            .get_mut(i)
            .ok_or_else(|| bad("preparation exceeds assembly boundary"))?
            .replace(index)
            .is_some()
        {
            return Err(bad("operation assigned to multiple rows"));
        }
    }
    Ok(steps)
}
pub(crate) fn validate_preparation_rows(
    plan: &DesignPlan,
    rows: &[PreparedRow],
    boundary: usize,
) -> Result<(), Error> {
    if boundary > plan.steps.len() || boundary > MAX_STEPS_PER_ROW * MAX_ROWS {
        return Err(Error::TooLarge);
    }
    let mut owners = vec![None; boundary];
    let mut cursor = 0;
    for (i, row) in rows.iter().enumerate() {
        prepared_chain(plan, row, i, &mut cursor, boundary, &mut owners)?;
    }
    if cursor != boundary || owners.iter().any(Option::is_none) {
        return Err(bad("unassigned or skipped preparation operation"));
    }
    Ok(())
}
/// Same explicit bounded policy as the native mismatch-PCR model; no binding
/// search is performed by this structural guard.
pub(crate) fn mismatch_pcr_shape(
    forward: &Primer,
    reverse: &Primer,
    options: &crate::cloning::pcr::mutagenesis::PcrOptions,
) -> Result<(), Error> {
    for p in [forward, reverse] {
        primer(p)?;
        if p.annealing
            .len()
            .checked_add(p.tail.len())
            .is_none_or(|n| n > 200)
            || options.binding.max_mismatches > p.annealing.len()
            || options.binding.min_three_prime_exact > p.annealing.len()
        {
            return Err(bad(
                "mismatch PCR requires bounded oligos and valid mismatch/3-prime policy",
            ));
        }
    }
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
        return Err(bad("invalid mismatch PCR search/materialization bounds"));
    }
    Ok(())
}
pub(crate) fn overlap_extension_catalog_shape(
    request: &crate::cloning::overlap::extension::pcr::Request,
    selection: &crate::cloning::restriction::catalog::CatalogSelection,
    fragment_index: usize,
) -> Result<(), Error> {
    request.validate_shape()?;
    selection.validate()?;
    if fragment_index > 100_000 {
        return Err(Error::TooLarge);
    }
    Ok(())
}
pub(crate) fn mismatch_pcr_catalog_shape(
    forward: &Primer,
    reverse: &Primer,
    options: &crate::cloning::pcr::mutagenesis::PcrOptions,
    selection: &crate::cloning::restriction::catalog::CatalogSelection,
    fragment_index: usize,
) -> Result<(), Error> {
    mismatch_pcr_shape(forward, reverse, options)?;
    selection.validate()?;
    if fragment_index > 100_000 {
        return Err(Error::TooLarge);
    }
    Ok(())
}
pub(crate) fn mismatch_pcr_work(
    n: usize,
    forward: &Primer,
    reverse: &Primer,
    options: &crate::cloning::pcr::mutagenesis::PcrOptions,
) -> Result<usize, Error> {
    mismatch_pcr_shape(forward, reverse, options)?;
    let binding = n
        .checked_mul(forward.annealing.len() + reverse.annealing.len())
        .and_then(|n| n.checked_mul(4))
        .ok_or(Error::TooLarge)?
        .min(
            options
                .binding
                .max_work
                .checked_mul(2)
                .ok_or(Error::TooLarge)?,
        );
    let pairs = n
        .checked_mul(n)
        .and_then(|n| n.checked_mul(4))
        .ok_or(Error::TooLarge)?
        .min(options.max_pair_checks);
    let products = n
        .checked_add(forward.tail.len())
        .and_then(|n| n.checked_add(reverse.tail.len()))
        .and_then(|n| n.checked_mul(options.max_products))
        .ok_or(Error::TooLarge)?
        .min(options.max_product_bases);
    binding
        .checked_add(pairs)
        .and_then(|n| n.checked_add(products))
        .ok_or(Error::TooLarge)
}

fn primer(p: &Primer) -> Result<(), Error> {
    token(&p.name, 256)?;
    if p.annealing.alphabet() != Alphabet::Dna
        || p.annealing.is_empty()
        || p.annealing.len() > 1000
        || p.tail.len() > 1000
        || !p
            .annealing
            .as_str()
            .bytes()
            .chain(p.tail.bytes())
            .all(|b| b"ACGT".contains(&b))
    {
        return Err(bad(
            "exact-PCR primers require concrete DNA, annealing1..1000 and tail0..1000 bases",
        ));
    }
    Ok(())
}
/// Borrowed guards for hosts before cloning source documents or hashing revisions.
/// Input IDs identify originals; distinct row roles may reuse the same original.
/// Document/event identity is captured by the host and is never inferred here.
pub fn source_preflight(
    inputs: &[(&str, &Molecule)],
    target: &Molecule,
    recipe: &Recipe,
) -> Result<(), Error> {
    if let Some(standard) = &recipe.standard {
        standard.validate()?;
    }

    recipe.options.validate()?;
    if recipe.fragments.is_empty() || recipe.fragments.len() > MAX_ROWS {
        return Err(Error::TooLarge);
    }
    if recipe.overlaps.len() != recipe.fragments.len() - 1 + usize::from(recipe.circular)
        || recipe.overlaps.iter().any(|n| *n == 0 || *n > 10_000)
    {
        return Err(bad(
            "one explicit positive overlap1..10000 is required per junction",
        ));
    }
    let total_homology = recipe
        .overlaps
        .iter()
        .try_fold(0usize, |a, b| a.checked_add(*b).ok_or(Error::TooLarge))?;
    if total_homology > 1_000_000 || total_homology > recipe.options.max_comparisons {
        return Err(Error::TooLarge);
    }
    preparation_source_preflight(inputs, target, &recipe.name, &recipe.fragments)?;
    hash(recipe)?;
    Ok(())
}
/// Shared pure preparation guard; independent of overlap/ligation chemistry.
pub fn preparation_source_preflight(
    inputs: &[(&str, &Molecule)],
    target: &Molecule,
    name: &str,
    fragments: &[FragmentSpec],
) -> Result<(), Error> {
    preparation_preflight_inner(inputs, Some(target), name, fragments)
}
/// Source/preparation-only guard for callers that independently reserve their
/// real intended target before constructing it. No dummy sequence is required.
pub(crate) fn preparation_rows_preflight(
    inputs: &[(&str, &Molecule)],
    name: &str,
    fragments: &[FragmentSpec],
) -> Result<(), Error> {
    preparation_preflight_inner(inputs, None, name, fragments)
}
fn preparation_preflight_inner(
    inputs: &[(&str, &Molecule)],
    target: Option<&Molecule>,
    name: &str,
    fragments: &[FragmentSpec],
) -> Result<(), Error> {
    if inputs.is_empty()
        || inputs.len() > MAX_ROWS
        || fragments.is_empty()
        || fragments.len() > MAX_ROWS
    {
        return Err(Error::TooLarge);
    }
    token(name, 1024)?;
    crate::material::capsule::material_capsules_preflight(fragments.iter().filter_map(|row| {
        match &row.preparation {
            Preparation::FromRun { origin, .. } => {
                Some(crate::material::capsule::MaterialOriginRef::Run(origin))
            }
            Preparation::FromDerived { origin, .. } => {
                Some(crate::material::capsule::MaterialOriginRef::Derived(origin))
            }
            _ => None,
        }
    }))?;
    let mut models = inputs.iter().map(|(_, m)| *m).collect::<Vec<_>>();
    if let Some(target) = target {
        models.push(target);
    }
    collagene_model::bounds::model_preflight(&models)?;
    if let Some(target) = target {
        if target.sequence.alphabet() != Alphabet::Dna {
            return Err(bad("independent target must be DNA"));
        }
        target.validate()?;
    }
    let mut lookup = BTreeMap::new();
    let mut bases = target.map_or(0, |m| m.sequence.len());
    for &(id, m) in inputs {
        token(id, 128)?;
        if lookup.insert(id, m).is_some() {
            return Err(bad("duplicate original input ID"));
        }
        if m.sequence.alphabet() != Alphabet::Dna {
            return Err(bad("source preparation requires DNA"));
        }
        m.validate()?;
        bases = bases
            .checked_add(m.sequence.len())
            .filter(|n| *n <= MAX_INPUT_BASES)
            .ok_or(Error::TooLarge)?;
    }
    let mut roles = BTreeSet::new();
    let mut role_bases = 0usize;
    let mut work = 0usize;
    let catalog = crate::sequence::analysis::starter_enzymes();
    let pinned_catalog = if fragments
        .iter()
        .any(|f| f.preparation.catalog_selection().is_some())
    {
        Some(crate::cloning::restriction::catalog::builtin_catalog()?)
    } else {
        None
    };
    for row in fragments {
        token(&row.role, 128)?;
        token(&row.input_id, 128)?;
        if !roles.insert(row.role.as_str()) {
            return Err(bad(
                "duplicate row role ID; reused originals require distinct role identities",
            ));
        }
        let m = lookup
            .get(row.input_id.as_str())
            .copied()
            .ok_or_else(|| bad("row refers to an unknown original input"))?;
        role_bases = role_bases
            .checked_add(m.sequence.len())
            .filter(|n| *n <= MAX_ROLE_BASES)
            .ok_or(Error::TooLarge)?;
        let effective;
        let preparation = if let Preparation::FromRun {
            origin,
            preparation,
        } = &row.preparation
        {
            from_run_shape(origin, preparation)?;
            if origin.selection.event.revision != m.revision() {
                return Err(bad("saved material origin source revision changed"));
            }
            effective = preparation.preparation()?;
            &effective
        } else if let Preparation::FromDerived {
            origin,
            preparation,
        } = &row.preparation
        {
            from_derived_shape(origin, preparation)?;
            if origin.product.revision != m.revision() {
                return Err(bad("saved derivation origin source revision changed"));
            }
            effective = preparation.preparation()?;
            &effective
        } else {
            &row.preparation
        };
        if row.generation.is_some()
            && matches!(
                row.preparation,
                Preparation::FromRun { .. }
                    | Preparation::FromDerived { .. }
                    | Preparation::DeclaredDuplex { .. }
                    | Preparation::AnnealedOligos { .. }
            )
        {
            return Err(bad(
                "Source generation may only be declared on an original input and cannot reset saved, declared-terminal or paired-oligo state.",
            ));
        }
        if row.generation.is_some() && crate::chemistry::has_unmapped_physical_evidence(m) {
            return Err(bad(
                "Imported physical evidence cannot become original generation authority.",
            ));
        }
        let (units, extra) = preparation_cost(
            m,
            preparation,
            &row.input_id,
            &lookup,
            &catalog,
            pinned_catalog.as_ref(),
            None,
        )?;
        role_bases = role_bases
            .checked_add(extra)
            .filter(|n| *n <= MAX_ROLE_BASES)
            .ok_or(Error::TooLarge)?;
        work = work
            .checked_add(units)
            .filter(|n| *n <= MAX_PREPARATION_WORK)
            .ok_or(Error::TooLarge)?;
    }
    hash(&(inputs, target, name, fragments))?;
    Ok(())
}
/// Existing native declared-preparation work, also used by aggregate batches.
/// This performs no search or replay and copies no closed origin. The flat child
/// conversion copies only the bounded current preparation declaration.
pub(crate) fn preparation_work(
    source: &Molecule,
    preparation: &Preparation,
) -> Result<usize, Error> {
    let child;
    let effective = match preparation {
        Preparation::FromRun { preparation, .. } | Preparation::FromDerived { preparation, .. } => {
            child = preparation.preparation()?;
            &child
        }
        _ => preparation,
    };
    let catalog = crate::sequence::analysis::starter_enzymes();
    let pinned = if effective.catalog_selection().is_some() {
        Some(crate::cloning::restriction::catalog::builtin_catalog()?)
    } else {
        None
    };
    preparation_cost(
        source,
        effective,
        "",
        &BTreeMap::new(),
        &catalog,
        pinned.as_ref(),
        None,
    )
    .map(|(work, _)| work)
}
/// Account for both ordered oligos and the actual bounded annealing search.
/// Uses the same admission/cost kernel as recipe compilation, without enumeration.
pub(crate) fn preparation_work_for_annealing(
    top: &Molecule,
    bottom: &Molecule,
    primary: &str,
    preparation: &Preparation,
) -> Result<usize, Error> {
    let Preparation::AnnealedOligos {
        bottom_input_id, ..
    } = preparation
    else {
        return Err(bad("paired work requires an annealed-oligo preparation"));
    };
    let lookup = BTreeMap::from([(primary, top), (bottom_input_id.as_str(), bottom)]);
    let (work, secondary_bases) =
        preparation_cost(top, preparation, primary, &lookup, &[], None, None)?;
    work.checked_add(secondary_bases).ok_or(Error::TooLarge)
}
/// Actual already-resolved source accounting. No origin replay or model clone;
/// callers which own a native context must not budget Bottom work as top length.
pub(crate) fn preparation_work_from_frame(
    source: crate::cloning::pcr::source::Template<'_>,
    preparation: &Preparation,
) -> Result<usize, Error> {
    let child;
    let effective = match preparation {
        Preparation::FromRun { preparation, .. } | Preparation::FromDerived { preparation, .. } => {
            child = preparation.preparation()?;
            &child
        }
        _ => preparation,
    };
    let strand = match effective {
        Preparation::ExactPcr {
            template_strand, ..
        }
        | Preparation::PcrCatalogDigest {
            template_strand, ..
        } => *template_strand,
        Preparation::OverlapExtensionPcr { request }
        | Preparation::OverlapExtensionPcrCatalogDigest { request, .. } => {
            let left = source.with_strand(request.left.template_strand)?;
            let right = source.with_strand(request.right.template_strand)?;
            if left.sequence.len() >= right.sequence.len() {
                left.strand
            } else {
                right.strand
            }
        }
        _ => crate::chemistry::Strand::Top,
    };
    let frame = source.with_strand(strand)?;
    let catalog = crate::sequence::analysis::starter_enzymes();
    let pinned = if effective.catalog_selection().is_some() {
        Some(crate::cloning::restriction::catalog::builtin_catalog()?)
    } else {
        None
    };
    preparation_cost(
        source.reference,
        effective,
        "",
        &BTreeMap::new(),
        &catalog,
        pinned.as_ref(),
        Some(frame),
    )
    .map(|(units, _)| units)
}
fn preparation_cost(
    m: &Molecule,
    preparation: &Preparation,
    input_id: &str,
    lookup: &BTreeMap<&str, &Molecule>,
    catalog: &[crate::sequence::analysis::Enzyme],
    pinned_catalog: Option<&crate::cloning::restriction::catalog::catalog::Catalog>,
    frame: Option<crate::cloning::pcr::source::Template<'_>>,
) -> Result<(usize, usize), Error> {
    if let Some(conditions) = preparation.reaction_conditions() {
        conditions.validate_for(
            pinned_catalog.ok_or(Error::TooLarge)?,
            preparation.catalog_selection().ok_or(Error::TooLarge)?,
        )?;
    }
    let source_bases = frame.map_or(m.sequence.len(), |v| v.sequence.len());
    let comparison_bases = frame
        .map(|v| v.comparison_bases())
        .transpose()?
        .unwrap_or(source_bases);
    let mut extra_bases = 0usize;
    let units = match preparation {
        Preparation::OverlapExtensionPcrCatalogDigest {
            request,
            selection,
            fragment_index,
            repair,
            ..
        } => {
            overlap_extension_catalog_shape(request, selection, *fragment_index)?;
            let (pcr_work, extra) = request.work(source_bases)?;
            let pcr_work = pcr_work
                .checked_mul(comparison_bases.div_ceil(source_bases))
                .ok_or(Error::TooLarge)?;
            let enzymes = selection.resolved_enzymes(pinned_catalog.ok_or(Error::TooLarge)?)?;
            let width = enzymes.iter().try_fold(0usize, |sum, e| {
                sum.checked_add(e.recognition.len() * 2)
                    .ok_or(Error::TooLarge)
            })?;
            let amplicon = source_bases
                .checked_mul(2)
                .and_then(|n| n.checked_add(request.left.forward.tail.len()))
                .and_then(|n| n.checked_add(request.left.reverse.tail.len()))
                .and_then(|n| n.checked_add(request.right.forward.tail.len()))
                .and_then(|n| n.checked_add(request.right.reverse.tail.len()))
                .and_then(|n| n.checked_sub(request.extension.overlap))
                .and_then(|n| n.checked_add(request.outer.forward.tail.len()))
                .and_then(|n| n.checked_add(request.outer.reverse.tail.len()))
                .ok_or(Error::TooLarge)?;
            let digest_material = amplicon
                .checked_mul(1 + usize::from(repair.is_some()))
                .ok_or(Error::TooLarge)?;
            extra_bases = extra_bases
                .checked_add(extra)
                .and_then(|n| n.checked_add(digest_material))
                .filter(|n| *n <= MAX_ROLE_BASES)
                .ok_or(Error::TooLarge)?;
            amplicon
                .checked_mul(width)
                .and_then(|n| n.checked_add(pcr_work))
                .ok_or(Error::TooLarge)?
        }
        Preparation::OverlapExtensionPcr { request } => {
            let (units, extra) = request.work(source_bases)?;
            let units = units
                .checked_mul(comparison_bases.div_ceil(source_bases))
                .ok_or(Error::TooLarge)?;
            extra_bases = extra_bases
                .checked_add(extra)
                .filter(|n| *n <= MAX_ROLE_BASES)
                .ok_or(Error::TooLarge)?;
            units
        }
        Preparation::FromRun { .. } | Preparation::FromDerived { .. } => {
            return Err(bad("nested material preparation is unsupported"));
        }
        Preparation::AnnealedOligos {
            bottom_input_id,
            options,
            candidate_index,
        } => {
            anneal_shape(bottom_input_id, options, *candidate_index)?;
            if bottom_input_id == input_id {
                return Err(bad(
                    "annealing requires two distinct original input identities",
                ));
            }
            let bottom = lookup
                .get(bottom_input_id.as_str())
                .copied()
                .ok_or_else(|| bad("unknown bottom oligo original input"))?;
            anneal_sources(m, bottom)?;
            extra_bases = extra_bases
                .checked_add(bottom.sequence.len())
                .filter(|n| *n <= MAX_ROLE_BASES)
                .ok_or(Error::TooLarge)?;
            // Reserve the actual exhaustive offset upper bound, capped by the recorded search budget.
            m.sequence
                .len()
                .checked_add(bottom.sequence.len())
                .and_then(|n| n.checked_mul(m.sequence.len().min(bottom.sequence.len())))
                .ok_or(Error::TooLarge)?
                .min(options.max_comparisons)
        }
        Preparation::DeclaredDuplex { .. } => {
            verification::sources::declared_source_termini_preflight(m)?;
            0
        }
        Preparation::WholeLinear => {
            if m.topology != Topology::Linear {
                return Err(bad(
                    "WholeLinear cannot materialize a circular source; specify actual PCR or digest preparation",
                ));
            }
            0
        }
        Preparation::MismatchPcr {
            forward,
            reverse,
            options,
        } => {
            mismatch_pcr_shape(forward, reverse, options)?;
            mismatch_pcr_work(m.sequence.len(), forward, reverse, options)?
        }
        Preparation::ExactPcr {
            forward, reverse, ..
        } => {
            primer(forward)?;
            primer(reverse)?;
            comparison_bases
                .checked_mul(forward.annealing.len() + reverse.annealing.len())
                .and_then(|n| n.checked_mul(4))
                .ok_or(Error::TooLarge)?
        }
        Preparation::PcrCatalogDigest {
            forward,
            reverse,
            selection,
            fragment_index,
            ..
        }
        | Preparation::MismatchPcrCatalogDigest {
            forward,
            reverse,
            selection,
            fragment_index,
            ..
        } => {
            if let Preparation::MismatchPcrCatalogDigest { options, .. } = preparation {
                mismatch_pcr_catalog_shape(forward, reverse, options, selection, *fragment_index)?;
            } else {
                primer(forward)?;
                primer(reverse)?;
            }
            if *fragment_index > 100_000 {
                return Err(Error::TooLarge);
            }
            let enzymes = selection.resolved_enzymes(pinned_catalog.ok_or(Error::TooLarge)?)?;
            let width = enzymes.iter().try_fold(0usize, |sum, e| {
                sum.checked_add(e.recognition.len() * 2)
                    .ok_or(Error::TooLarge)
            })?;
            let amplicon = source_bases
                .checked_add(forward.tail.len())
                .and_then(|n| n.checked_add(reverse.tail.len()))
                .ok_or(Error::TooLarge)?;
            extra_bases = extra_bases
                .checked_add(forward.tail.len())
                .and_then(|n| n.checked_add(reverse.tail.len()))
                .filter(|n| *n <= MAX_ROLE_BASES)
                .ok_or(Error::TooLarge)?;
            let pcr = if let Preparation::MismatchPcrCatalogDigest { options, .. } = preparation {
                mismatch_pcr_work(m.sequence.len(), forward, reverse, options)?
            } else {
                comparison_bases
                    .checked_mul(forward.annealing.len() + reverse.annealing.len())
                    .and_then(|n| n.checked_mul(4))
                    .ok_or(Error::TooLarge)?
            };
            amplicon
                .checked_mul(width)
                .and_then(|n| n.checked_add(pcr))
                .ok_or(Error::TooLarge)?
        }
        Preparation::CatalogDigest {
            selection,
            fragment_index,
            ..
        } => {
            if *fragment_index > 100_000 {
                return Err(Error::TooLarge);
            }
            let enzymes = selection.resolved_enzymes(pinned_catalog.ok_or(Error::TooLarge)?)?;
            let width = enzymes.iter().try_fold(0usize, |sum, e| {
                sum.checked_add(e.recognition.len() * 2)
                    .ok_or(Error::TooLarge)
            })?;
            m.sequence.len().checked_mul(width).ok_or(Error::TooLarge)?
        }
        Preparation::Digest {
            enzymes,
            fragment_index,
            ..
        } => {
            if enzymes.is_empty() || enzymes.len() > 32 || *fragment_index > 100_000 {
                return Err(Error::TooLarge);
            }
            let mut names = BTreeSet::new();
            let mut width = 0usize;
            for name in enzymes {
                token(name, 128)?;
                if !names.insert(name) {
                    return Err(bad("duplicate digest enzyme"));
                }
                let enzyme = catalog
                    .iter()
                    .find(|e| e.name.eq_ignore_ascii_case(name))
                    .ok_or_else(|| {
                        bad("digest enzyme is not supported by the existing typed Digest operation")
                    })?;
                width = width
                    .checked_add(enzyme.recognition.len() * 2)
                    .ok_or(Error::TooLarge)?;
            }
            m.sequence.len().checked_mul(width).ok_or(Error::TooLarge)?
        }
    };
    Ok((units, extra_bases))
}
pub fn preflight(inputs: &[PlanInput], target: &Molecule, recipe: &Recipe) -> Result<(), Error> {
    if inputs.is_empty() || inputs.len() > MAX_ROWS {
        return Err(Error::TooLarge);
    }
    source_preflight(
        &inputs
            .iter()
            .map(|i| (i.id.as_str(), &i.molecule))
            .collect::<Vec<_>>(),
        target,
        recipe,
    )?;
    if inputs.iter().any(|i| !hex(&i.revision)) {
        return Err(bad("source revisions must be explicit canonical hashes"));
    }
    Ok(())
}
pub(crate) fn bound_revisions(inputs: &[PlanInput], target: &Molecule) -> Result<String, Error> {
    let target_revision = hash(target)?;
    for input in inputs {
        let actual = hash(&input.molecule)?;
        if actual != input.revision {
            return Err(bad("original input revision changed"));
        }
    }
    Ok(target_revision)
}
pub fn preview_revision(
    inputs: &[PlanInput],
    target: &Molecule,
    recipe: &Recipe,
) -> Result<String, Error> {
    preflight(inputs, target, recipe)?;
    bound_revisions(inputs, target)?;
    hash(&(ENGINE_VERSION, "preview", inputs, target, recipe))
}
/// Pure compilation: a structurally valid plan is not a verified assembly.
/// Original immutable inputs and actual preparation operations, before assembly.
pub struct PreparedCompilation {
    pub inputs: Vec<PlanInput>,
    pub target_revision: String,
    pub steps: Vec<Operation>,
    pub rows: Vec<PreparedRow>,
}
pub fn compile_preparations(
    inputs: &[PlanInput],
    target: &Molecule,
    name: &str,
    fragments: &[FragmentSpec],
    policy: CleavagePolicy,
) -> Result<PreparedCompilation, Error> {
    if inputs.is_empty() || inputs.len() > MAX_ROWS {
        return Err(Error::TooLarge);
    }
    preparation_source_preflight(
        &inputs
            .iter()
            .map(|i| (i.id.as_str(), &i.molecule))
            .collect::<Vec<_>>(),
        target,
        name,
        fragments,
    )?;
    if inputs.iter().any(|i| !hex(&i.revision)) {
        return Err(bad("source revisions must be explicit canonical hashes"));
    }
    let target_revision = bound_revisions(inputs, target)?;
    let (steps, rows) = preparation_steps(inputs, fragments, policy)?;
    Ok(PreparedCompilation {
        inputs: inputs.to_vec(),
        target_revision,
        steps,
        rows,
    })
}
/// Pure compiler shared by readonly recipe correspondence checks. Callers must
/// perform borrowed source/recipe bounds first; no molecular operation runs.
pub(crate) fn preparation_steps(
    inputs: &[PlanInput],
    fragments: &[FragmentSpec],
    policy: CleavagePolicy,
) -> Result<(Vec<Operation>, Vec<PreparedRow>), Error> {
    let mut steps = Vec::new();
    let mut rows = Vec::with_capacity(fragments.len());
    for (index, row) in fragments.iter().enumerate() {
        let original = inputs
            .iter()
            .find(|i| i.id == row.input_id)
            .ok_or_else(|| bad("missing original"))?;
        let original_input = Source::Input(row.input_id.clone());
        let effective;
        let (input, first_step, preparation) = if let Preparation::FromRun {
            origin,
            preparation,
        } = &row.preparation
        {
            let first = steps.len();
            steps.push(Operation::ReconstructRunProduct {
                template: original_input,
                origin: origin.clone(),
            });
            effective = preparation.preparation()?;
            (Source::Step(first), Some(first), &effective)
        } else if let Preparation::FromDerived {
            origin,
            preparation,
        } = &row.preparation
        {
            let first = steps.len();
            steps.push(Operation::ReconstructMaterial {
                template: original_input,
                origin: origin.clone(),
            });
            effective = preparation.preparation()?;
            (Source::Step(first), Some(first), &effective)
        } else if let Some(assumptions) = &row.generation {
            let first = steps.len();
            steps.push(Operation::DeclareSourceModifications {
                template: original_input,
                assumptions: assumptions.clone(),
            });
            (Source::Step(first), Some(first), &row.preparation)
        } else {
            (original_input, None, &row.preparation)
        };
        if let Preparation::OverlapExtensionPcr { request }
        | Preparation::OverlapExtensionPcrCatalogDigest { request, .. } = preparation
        {
            let preparation_step = first_step.unwrap_or(steps.len());
            let mut prepared_source = request.append(input, &mut steps)?;
            let mut repair_step = None;
            if let Preparation::OverlapExtensionPcrCatalogDigest {
                selection,
                conditions,
                fragment_index,
                repair,
                ..
            } = preparation
            {
                let index = steps.len();
                steps.push(digest_operation(
                    prepared_source,
                    selection.clone(),
                    *fragment_index,
                    policy,
                    conditions.clone(),
                ));
                prepared_source = Source::Step(index);
                if let Some(method) = repair {
                    let index = steps.len();
                    steps.push(Operation::RepairEnds {
                        template: prepared_source,
                        method: *method,
                    });
                    prepared_source = Source::Step(index);
                    repair_step = Some(index);
                }
            }
            rows.push(PreparedRow {
                row_index: index,
                role: row.role.clone(),
                input_id: row.input_id.clone(),
                input_revision: original.revision.clone(),
                preparation_step: Some(preparation_step),
                repair_step,
                prepared_source,
                reverse: row.reverse,
            });
            continue;
        }
        let operation = match preparation {
            Preparation::OverlapExtensionPcr { .. }
            | Preparation::OverlapExtensionPcrCatalogDigest { .. } => {
                unreachable!("handled explicit fork above")
            }
            Preparation::FromRun { .. } | Preparation::FromDerived { .. } => {
                return Err(bad("nested material preparation is unsupported"));
            }
            Preparation::WholeLinear => None,
            Preparation::AnnealedOligos {
                bottom_input_id,
                options,
                candidate_index,
            } => Some(Operation::AnnealOligos {
                top: input.clone(),
                bottom: Source::Input(bottom_input_id.clone()),
                options: options.clone(),
                candidate_index: *candidate_index,
            }),
            Preparation::DeclaredDuplex { left, right } => Some(Operation::DeclareSourceTermini {
                template: input.clone(),
                left: *left,
                right: *right,
            }),
            Preparation::CatalogDigest {
                selection,
                conditions,
                fragment_index,
                ..
            } => Some(digest_operation(
                input.clone(),
                selection.clone(),
                *fragment_index,
                policy,
                conditions.clone(),
            )),
            Preparation::Digest {
                enzymes,
                fragment_index,
                ..
            } => Some(match policy {
                CleavagePolicy::AssumedComplete => Operation::Digest {
                    template: input.clone(),
                    enzymes: enzymes.clone(),
                    fragment_index: *fragment_index,
                },
                CleavagePolicy::PositionalChecked => digest_operation(
                    input.clone(),
                    starter_selection(enzymes)?,
                    *fragment_index,
                    policy,
                    None,
                ),
            }),
            Preparation::MismatchPcr {
                forward,
                reverse,
                options,
            }
            | Preparation::MismatchPcrCatalogDigest {
                forward,
                reverse,
                options,
                ..
            } => Some(Operation::PcrWithMismatches {
                template: input.clone(),
                forward: forward.clone(),
                reverse: reverse.clone(),
                options: options.clone(),
            }),
            Preparation::ExactPcr {
                forward,
                reverse,
                template_strand,
            }
            | Preparation::PcrCatalogDigest {
                forward,
                reverse,
                template_strand,
                ..
            } => Some(Operation::Pcr {
                template: input.clone(),
                forward: forward.clone(),
                reverse: reverse.clone(),
                template_strand: *template_strand,
            }),
        };
        let preparation_step = first_step.or_else(|| operation.as_ref().map(|_| steps.len()));
        let mut prepared_source = if let Some(operation) = operation {
            let step = steps.len();
            steps.push(operation);
            Source::Step(step)
        } else {
            input
        };
        if let Preparation::PcrCatalogDigest {
            selection,
            conditions,
            fragment_index,
            ..
        }
        | Preparation::MismatchPcrCatalogDigest {
            selection,
            conditions,
            fragment_index,
            ..
        } = preparation
        {
            let index = steps.len();
            steps.push(digest_operation(
                prepared_source,
                selection.clone(),
                *fragment_index,
                policy,
                conditions.clone(),
            ));
            prepared_source = Source::Step(index);
        }
        let repair_step = if let Preparation::Digest {
            repair: Some(method),
            ..
        }
        | Preparation::CatalogDigest {
            repair: Some(method),
            ..
        }
        | Preparation::PcrCatalogDigest {
            repair: Some(method),
            ..
        }
        | Preparation::MismatchPcrCatalogDigest {
            repair: Some(method),
            ..
        } = preparation
        {
            let index = steps.len();
            steps.push(Operation::RepairEnds {
                template: prepared_source,
                method: *method,
            });
            prepared_source = Source::Step(index);
            Some(index)
        } else {
            None
        };
        rows.push(PreparedRow {
            row_index: index,
            role: row.role.clone(),
            input_id: row.input_id.clone(),
            input_revision: original.revision.clone(),
            preparation_step,
            repair_step,
            prepared_source: prepared_source.clone(),
            reverse: row.reverse,
        });
    }
    Ok((steps, rows))
}

pub fn compile(
    inputs: &[PlanInput],
    target: &Molecule,
    recipe: &Recipe,
) -> Result<Compiled, Error> {
    preflight(inputs, target, recipe)?;
    let prepared = compile_preparations(
        inputs,
        target,
        &recipe.name,
        &recipe.fragments,
        CleavagePolicy::AssumedComplete,
    )?;
    let mut steps = prepared.steps;
    let rows = prepared.rows;
    let fragments = rows
        .iter()
        .map(|row| OrientedSource {
            source: row.prepared_source.clone(),
            reverse: row.reverse,
        })
        .collect();
    let assembly_step = steps.len();
    steps.push(Operation::ProfiledOverlap {
        fragments,
        overlaps: recipe.overlaps.clone(),
        circular: recipe.circular,
        options: recipe.options.clone(),
    });
    let compiled = Compiled {
        plan: DesignPlan {
            standard: recipe.standard.clone(),
            name: recipe.name.clone(),
            target_revision: prepared.target_revision,
            inputs: prepared.inputs,
            steps,
        },
        rows,
        assembly_step,
    };
    hash(&compiled)?;
    Ok(compiled)
}
pub fn preview(inputs: &[PlanInput], target: &Molecule, recipe: &Recipe) -> Result<Preview, Error> {
    let context_revision = preview_revision(inputs, target, recipe)?;
    let compiled = compile(inputs, target, recipe)?;
    finish_preview(context_revision, compiled, target, recipe)
}

fn finish_preview(
    context_revision: String,
    compiled: Compiled,
    target: &Molecule,
    recipe: &Recipe,
) -> Result<Preview, Error> {
    // Exactly the existing verifier prepares PCR/digest products and projects
    // annotations. Neither the target nor a preview chooses a fragment or geometry.
    let mut verification = verification::engine::verify(&compiled.plan, target);
    let available = compiled
        .rows
        .iter()
        .map(|row| match &row.prepared_source {
            Source::Input(id) => compiled
                .plan
                .inputs
                .iter()
                .find(|input| &input.id == id)
                .map(|input| &input.molecule),
            Source::Step(index) => verification.product(*index),
        })
        .collect::<Vec<_>>();
    let preparation_complete = available.iter().all(Option::is_some);
    let models = available
        .iter()
        .filter_map(|model| *model)
        .collect::<Vec<_>>();
    collagene_model::bounds::model_preflight(&models)?;
    let mut bases = 0usize;
    for model in &models {
        bases = bases
            .checked_add(model.sequence.len())
            .filter(|n| *n <= MAX_ROLE_BASES)
            .ok_or(Error::TooLarge)?;
    }
    hash(&models)?;
    let prepared = available
        .into_iter()
        .map(|model| model.cloned())
        .collect::<Vec<_>>();
    // Complete-only collection preserves original row order and moves, rather
    // than clones, the already admitted documents. Sparse rows stay sparse.
    let (prepared, profile) = if preparation_complete {
        let complete = prepared
            .into_iter()
            .collect::<Option<Vec<_>>>()
            .ok_or_else(|| bad("complete preparation lost its original-row product"))?;
        let profile = match overlap_profiles::check_oriented_fragments(
            &complete,
            &compiled
                .rows
                .iter()
                .map(|row| row.reverse)
                .collect::<Vec<_>>(),
            &recipe.overlaps,
            recipe.circular,
            &recipe.options,
        ) {
            Ok(report) => Some(report),
            Err(error) => {
                use crate::verification::checks::{CheckFinding, CheckScope, CheckVerdict, RuleId};
                verification.findings.push(CheckFinding::new(
                    RuleId::HomologyPolicy,
                    CheckVerdict::Indeterminate,
                    CheckScope::Step {
                        step_index: compiled.assembly_step,
                    },
                    format!("Recipe advisory junction profile withheld: {error}"),
                ));
                if verification.status
                    == crate::verification::evidence::Status::ConsistentUnderAssumptions
                {
                    verification.status = crate::verification::evidence::Status::Incomplete;
                }
                None
            }
        };
        (complete.into_iter().map(Some).collect(), profile)
    } else {
        (prepared, None)
    };
    let result = Preview {
        context_revision,
        compiled,
        verification,
        prepared,
        profile,
        preparation_complete,
    };
    result.validate()?;
    Ok(result)
}
/// Physical homology requests retain the original preparation recipe and require
/// its geometry to agree exactly with the separate versioned material policy.
pub fn homology_source_preflight(
    inputs: &[(&str, &Molecule)],
    target: &Molecule,
    recipe: &Recipe,
    policy: &crate::cloning::overlap::assembly::Options,
) -> Result<(), Error> {
    policy.validate_shape(recipe.fragments.len())?;
    if recipe.options.method == overlap_profiles::OverlapMethod::InFusion {
        return Err(bad(
            "In-Fusion requires its distinct 3′ processing/host-rescue model. Keep reference-only interpretation; this declared 5′ resection/fill/ligase model cannot establish an In-Fusion physical product.",
        ));
    }
    if policy.overlaps != recipe.overlaps
        || policy.topology
            != if recipe.circular {
                Topology::Circular
            } else {
                Topology::Linear
            }
    {
        return Err(bad(
            "physical policy geometry differs from the declared recipe",
        ));
    }
    source_preflight(inputs, target, recipe)
}
pub fn homology_preflight(
    inputs: &[PlanInput],
    target: &Molecule,
    recipe: &Recipe,
    policy: &crate::cloning::overlap::assembly::Options,
) -> Result<(), Error> {
    homology_source_preflight(
        &inputs
            .iter()
            .map(|i| (i.id.as_str(), &i.molecule))
            .collect::<Vec<_>>(),
        target,
        recipe,
        policy,
    )?;
    preflight(inputs, target, recipe)
}
pub fn homology_preview_revision(
    inputs: &[PlanInput],
    target: &Molecule,
    recipe: &Recipe,
    policy: &crate::cloning::overlap::assembly::Options,
) -> Result<String, Error> {
    homology_preflight(inputs, target, recipe, policy)?;
    bound_revisions(inputs, target)?;
    hash(&(
        ENGINE_VERSION,
        crate::cloning::overlap::assembly::ENGINE_VERSION,
        "physical-preview",
        inputs,
        target,
        recipe,
        policy,
    ))
}
pub fn compile_homology(
    inputs: &[PlanInput],
    target: &Molecule,
    recipe: &Recipe,
    policy: &crate::cloning::overlap::assembly::Options,
) -> Result<Compiled, Error> {
    homology_preflight(inputs, target, recipe, policy)?;
    let mut compiled = compile(inputs, target, recipe)?;
    let last = compiled
        .plan
        .steps
        .last_mut()
        .ok_or_else(|| bad("missing assembly operation"))?;
    let Operation::ProfiledOverlap {
        fragments, options, ..
    } = last
    else {
        return Err(bad("unexpected final preparation operation"));
    };
    *last = Operation::HomologyAssembly {
        fragments: fragments.clone(),
        options: options.clone(),
        policy: policy.clone(),
    };
    hash(&compiled)?;
    Ok(compiled)
}
pub fn preview_homology(
    inputs: &[PlanInput],
    target: &Molecule,
    recipe: &Recipe,
    policy: &crate::cloning::overlap::assembly::Options,
) -> Result<Preview, Error> {
    let context = homology_preview_revision(inputs, target, recipe, policy)?;
    let compiled = compile_homology(inputs, target, recipe, policy)?;
    finish_preview(context, compiled, target, recipe)
}

/// Separate generic3′ processing/annealing; no in-vitro sealing or cellular rescue.
pub fn annealing_source_preflight(
    inputs: &[(&str, &Molecule)],
    target: &Molecule,
    recipe: &Recipe,
    policy: &crate::cloning::overlap::annealing::Options,
) -> Result<(), Error> {
    policy.validate_shape(recipe.fragments.len())?;
    policy.validate_profile(recipe.options.method)?;
    if recipe.fragments.len() == 1 && policy.topology != Topology::Circular {
        return Err(bad(
            "Single linear3′ resection uses a standalone ResectThreePrime operation; an annealing recipe must join2+ sources or self-circularize one source",
        ));
    }
    if policy.overlaps != recipe.overlaps
        || policy.topology
            != if recipe.circular {
                Topology::Circular
            } else {
                Topology::Linear
            }
    {
        return Err(bad(
            "physical policy geometry differs from the declared recipe",
        ));
    }
    source_preflight(inputs, target, recipe)
}
pub fn annealing_preflight(
    inputs: &[PlanInput],
    target: &Molecule,
    recipe: &Recipe,
    policy: &crate::cloning::overlap::annealing::Options,
) -> Result<(), Error> {
    annealing_source_preflight(
        &inputs
            .iter()
            .map(|i| (i.id.as_str(), &i.molecule))
            .collect::<Vec<_>>(),
        target,
        recipe,
        policy,
    )?;
    preflight(inputs, target, recipe)
}
pub fn annealing_preview_revision(
    inputs: &[PlanInput],
    target: &Molecule,
    recipe: &Recipe,
    policy: &crate::cloning::overlap::annealing::Options,
) -> Result<String, Error> {
    annealing_preflight(inputs, target, recipe, policy)?;
    bound_revisions(inputs, target)?;
    hash(&(
        ENGINE_VERSION,
        crate::cloning::overlap::annealing::ENGINE_VERSION,
        "three-prime-annealing-preview",
        inputs,
        target,
        recipe,
        policy,
    ))
}
pub fn compile_annealing(
    inputs: &[PlanInput],
    target: &Molecule,
    recipe: &Recipe,
    policy: &crate::cloning::overlap::annealing::Options,
) -> Result<Compiled, Error> {
    annealing_preflight(inputs, target, recipe, policy)?;
    let mut compiled = compile(inputs, target, recipe)?;
    let last = compiled
        .plan
        .steps
        .last_mut()
        .ok_or_else(|| bad("missing assembly operation"))?;
    let Operation::ProfiledOverlap {
        fragments, options, ..
    } = last
    else {
        return Err(bad("unexpected final preparation operation"));
    };
    *last = Operation::HomologyAnnealing {
        fragments: fragments.clone(),
        options: options.clone(),
        policy: policy.clone(),
    };
    hash(&compiled)?;
    Ok(compiled)
}
pub fn preview_annealing(
    inputs: &[PlanInput],
    target: &Molecule,
    recipe: &Recipe,
    policy: &crate::cloning::overlap::annealing::Options,
) -> Result<Preview, Error> {
    let context = annealing_preview_revision(inputs, target, recipe, policy)?;
    let compiled = compile_annealing(inputs, target, recipe, policy)?;
    finish_preview(context, compiled, target, recipe)
}

impl Preview {
    /// Readonly shape/binding checks only. No reconstruction or authority restoration.
    pub fn validate(&self) -> Result<(), Error> {
        self.verification.validate_documentary()?;
        self.verification
            .validate_standard_for_plan(&self.compiled.plan)?;
        if !hex(&self.context_revision)
            || self.compiled.rows.is_empty()
            || self.compiled.rows.len() > MAX_ROWS
            || self.prepared.len() != self.compiled.rows.len()
            || self.compiled.plan.inputs.len() > MAX_ROWS
            || self.verification.steps.len() > MAX_STEPS_PER_ROW * MAX_ROWS + 1
            || self.verification.planned_steps != self.compiled.plan.steps.len()
            || self.verification.step_reports.len() > MAX_STEPS_PER_ROW * MAX_ROWS + 1
            || self.compiled.plan.steps.len() > MAX_STEPS_PER_ROW * MAX_ROWS + 1
            || self.compiled.assembly_step.checked_add(1) != Some(self.compiled.plan.steps.len())
        {
            return Err(bad("invalid preview shape"));
        }
        if self.preparation_complete != self.prepared.iter().all(Option::is_some)
            || self.profile.is_some() && !self.preparation_complete
            || self.verification.target_revision != self.compiled.plan.target_revision
        {
            return Err(bad("preview plan/evidence binding mismatch"));
        }
        let models = self
            .compiled
            .plan
            .inputs
            .iter()
            .map(|i| &i.molecule)
            .chain(self.prepared.iter().filter_map(Option::as_ref))
            .chain(
                self.verification
                    .steps
                    .iter()
                    .filter_map(|step| step.product.as_ref()),
            )
            .collect::<Vec<_>>();
        for model in models {
            collagene_model::bounds::model_preflight(&[model])?;
        }
        for report in &self.verification.step_reports {
            report.validate_documentary()?;
            if report.step_index >= self.compiled.plan.steps.len() {
                return Err(Error::TooLarge);
            }
        }
        validate_prepared_inputs(&self.compiled.plan.inputs, &self.compiled.plan.steps)?;
        if self.verification.plan_revision != hash(&self.compiled.plan)? {
            return Err(bad("preview plan revision changed"));
        }
        for (index, row) in self.compiled.rows.iter().enumerate() {
            validate_annealed_row(row, &self.compiled.plan.steps)?;
            token(&row.role, 128)?;
            token(&row.input_id, 128)?;
            if row.row_index != index || !hex(&row.input_revision) {
                return Err(bad("invalid preview row identity"));
            }
            let input = self
                .compiled
                .plan
                .inputs
                .iter()
                .find(|i| i.id == row.input_id)
                .ok_or_else(|| bad("preview row lost its source"))?;
            if input.revision != row.input_revision {
                return Err(bad("preview row revision mismatch"));
            }
        }
        validate_preparation_rows(
            &self.compiled.plan,
            &self.compiled.rows,
            self.compiled.assembly_step,
        )?;
        let (fragments, overlaps, circular, options) = match self.compiled.plan.steps.last() {
            Some(Operation::ProfiledOverlap {
                fragments,
                overlaps,
                circular,
                options,
            }) => (fragments, overlaps, *circular, options),
            Some(Operation::HomologyAssembly {
                fragments,
                options,
                policy,
            }) => {
                policy.validate_shape(fragments.len())?;
                (
                    fragments,
                    &policy.overlaps,
                    policy.topology == Topology::Circular,
                    options,
                )
            }
            Some(Operation::HomologyAnnealing {
                fragments,
                options,
                policy,
            }) => {
                policy.validate_shape(fragments.len())?;
                if fragments.len() == 1 && policy.topology != Topology::Circular {
                    return Err(bad("single-linear processing is not an annealing assembly"));
                }
                if self.verification.final_product().is_some() {
                    let actual = self
                        .verification
                        .final_product()
                        .ok_or_else(|| bad("missing annealing diagnostic reference"))?;
                    let evidence=self.verification.physical_processing_reports.iter().find_map(|step|if step.step_index==self.compiled.assembly_step{match &step.evidence{crate::verification::evidence::PhysicalProcessingEvidence::HomologyAnnealing(report)=>Some(report),_=>None}}else{None}).ok_or_else(||bad("missing exact annealing strand evidence"))?;
                    evidence.validate_documentary()?;
                    if hash(policy)? != hash(&evidence.options)?
                        || evidence.preparations.len() != fragments.len()
                        || evidence.reference != actual.sequence
                        || policy.topology != actual.topology
                        || self.verification.status
                            == crate::verification::evidence::Status::ConsistentUnderAssumptions
                    {
                        return Err(bad("annealing diagnostic/policy/status binding differs"));
                    }
                }
                (
                    fragments,
                    &policy.overlaps,
                    policy.topology == Topology::Circular,
                    options,
                )
            }
            _ => {
                return Err(bad(
                    "preview has no final reference or physical homology operation",
                ));
            }
        };
        if fragments.len() != self.compiled.rows.len()
            || self.compiled.plan.steps[..self.compiled.assembly_step]
                .iter()
                .any(|step| {
                    !matches!(
                        step,
                        Operation::Pcr { .. }
                            | Operation::PcrWithMismatches { .. }
                            | Operation::OverlapExtension { .. }
                            | Operation::ReconstructRunProduct { .. }
                            | Operation::ReconstructMaterial { .. }
                            | Operation::AnnealOligos { .. }
                            | Operation::Digest { .. }
                            | Operation::CatalogDigest { .. }
                            | Operation::RepairEnds { .. }
                            | Operation::DeclareSourceTermini { .. }
                    )
                })
        {
            return Err(bad("preview operation shape changed"));
        }
        for (fragment, row) in fragments.iter().zip(&self.compiled.rows) {
            let same = match (&fragment.source, &row.prepared_source) {
                (Source::Input(a), Source::Input(b)) => a == b,
                (Source::Step(a), Source::Step(b)) => a == b,
                _ => false,
            };
            if !same || fragment.reverse != row.reverse {
                return Err(bad("final overlap source mapping changed"));
            }
        }
        if let Some(profile) = &self.profile {
            profile.validate()?;
            if &profile.overlaps != overlaps
                || profile.circular != circular
                || hash(&profile.options)? != hash(options)?
            {
                return Err(bad("profile recipe options changed"));
            }
            if profile.source_revisions
                != self
                    .prepared
                    .iter()
                    .map(|model| {
                        model
                            .as_ref()
                            .ok_or_else(|| bad("profile requires every original-row product"))
                            .and_then(hash)
                    })
                    .collect::<Result<Vec<_>, _>>()?
                || profile.orientations
                    != self
                        .compiled
                        .rows
                        .iter()
                        .map(|r| r.reverse)
                        .collect::<Vec<_>>()
            {
                return Err(bad("profile row bindings changed"));
            }
        }
        hash(self)?;
        Ok(())
    }
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        self.validate()?;
        if self.context_revision != expected {
            return Err(bad(
                "preview belongs to a different source/target/recipe context",
            ));
        }
        Ok(())
    }
}

impl OverlapFinishing {
    pub fn validate_geometry(
        &self,
        rows: usize,
        overlaps: &[usize],
        circular: bool,
        options: &overlap_profiles::ProfileOptions,
    ) -> Result<(), Error> {
        options.validate()?;
        let topology = if circular {
            Topology::Circular
        } else {
            Topology::Linear
        };
        let (actual_overlaps, actual_topology) = match self {
            Self::ReferenceOnly => return Ok(()),
            Self::HomologyAssembly { policy } => {
                policy.validate_shape(rows)?;
                if options.method == overlap_profiles::OverlapMethod::InFusion {
                    return Err(bad(
                        "In-Fusion requires its distinct 3′ annealing and planned host-rescue interpretation, not 5′ fill/seal processing",
                    ));
                }
                (&policy.overlaps, policy.topology)
            }
            Self::HomologyAnnealing { policy } => {
                policy.validate_shape(rows)?;
                policy.validate_profile(options.method)?;
                if rows == 1 && !circular {
                    return Err(bad(
                        "a single linear source requires standalone 3′ resection rather than an annealing assembly",
                    ));
                }
                (&policy.overlaps, policy.topology)
            }
        };
        if actual_overlaps != overlaps || actual_topology != topology {
            return Err(bad(
                "finishing policy differs from the independent declared junctions or topology",
            ));
        }
        Ok(())
    }
    pub fn preflight(
        &self,
        inputs: &[PlanInput],
        target: &Molecule,
        recipe: &Recipe,
    ) -> Result<(), Error> {
        self.validate_geometry(
            recipe.fragments.len(),
            &recipe.overlaps,
            recipe.circular,
            &recipe.options,
        )?;
        match self {
            Self::ReferenceOnly => preflight(inputs, target, recipe),
            Self::HomologyAssembly { policy } => homology_preflight(inputs, target, recipe, policy),
            Self::HomologyAnnealing { policy } => {
                annealing_preflight(inputs, target, recipe, policy)
            }
        }
    }
    pub fn compile(
        &self,
        inputs: &[PlanInput],
        target: &Molecule,
        recipe: &Recipe,
    ) -> Result<Compiled, Error> {
        self.preflight(inputs, target, recipe)?;
        match self {
            Self::ReferenceOnly => compile(inputs, target, recipe),
            Self::HomologyAssembly { policy } => compile_homology(inputs, target, recipe, policy),
            Self::HomologyAnnealing { policy } => compile_annealing(inputs, target, recipe, policy),
        }
    }
    pub fn preview_revision(
        &self,
        inputs: &[PlanInput],
        target: &Molecule,
        recipe: &Recipe,
    ) -> Result<String, Error> {
        self.preflight(inputs, target, recipe)?;
        match self {
            Self::ReferenceOnly => preview_revision(inputs, target, recipe),
            Self::HomologyAssembly { policy } => {
                homology_preview_revision(inputs, target, recipe, policy)
            }
            Self::HomologyAnnealing { policy } => {
                annealing_preview_revision(inputs, target, recipe, policy)
            }
        }
    }
    pub fn preview(
        &self,
        inputs: &[PlanInput],
        target: &Molecule,
        recipe: &Recipe,
    ) -> Result<Preview, Error> {
        self.preflight(inputs, target, recipe)?;
        match self {
            Self::ReferenceOnly => preview(inputs, target, recipe),
            Self::HomologyAssembly { policy } => preview_homology(inputs, target, recipe, policy),
            Self::HomologyAnnealing { policy } => preview_annealing(inputs, target, recipe, policy),
        }
    }
}
