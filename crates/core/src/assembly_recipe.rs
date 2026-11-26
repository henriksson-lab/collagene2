//! Compile explicit preparation rows into existing typed operations. Fresh
//! reconstruction remains the verifier's job; this module never slices a source
//! document and calls that physical preparation, or searches toward a target.
use crate::{
    Alphabet, Error, Molecule, Topology,
    analysis::Primer,
    overlap_profiles::{self, ProfileOptions, ProfileReport},
    verification::{
        self, DesignPlan, Operation, OrientedSource, PlanInput, Source, VerificationReport,
    },
};
use serde::{Deserialize, Serialize};
use std::{
    collections::{BTreeMap, BTreeSet},
    io::Write,
};
pub const ENGINE_VERSION: &str = "overlap-recipe-1";
pub const MAX_ROWS: usize = 128;
pub(crate) const MAX_BYTES: usize = 32 << 20;
const MAX_INPUT_BASES: usize = 8 << 20;
pub(crate) const MAX_ROLE_BASES: usize = 16 << 20;
const MAX_PREPARATION_WORK: usize = 100_000_000;
pub(crate) const MAX_STEPS_PER_ROW: usize = 7;
macro_rules! dto{($($t:tt)*)=>{#[derive(Clone,Debug,Serialize,Deserialize)]#[cfg_attr(feature="json-schema",derive(schemars::JsonSchema))]#[serde(deny_unknown_fields)]$($t)*};}
dto! {pub enum Preparation {WholeLinear,Digest{enzymes:Vec<String>,fragment_index:usize,repair:Option<crate::chemistry::EndRepairMethod>},ExactPcr{forward:Primer,reverse:Primer},DeclaredDuplex{left:crate::chemistry::Terminus,right:crate::chemistry::Terminus},CatalogDigest{selection:crate::catalog_digest::CatalogSelection,fragment_index:usize,repair:Option<crate::chemistry::EndRepairMethod>},PcrCatalogDigest{forward:Primer,reverse:Primer,selection:crate::catalog_digest::CatalogSelection,fragment_index:usize,repair:Option<crate::chemistry::EndRepairMethod>},AnnealedOligos{bottom_input_id:String,options:crate::chemistry::AnnealOptions,candidate_index:usize},FromRun{origin:crate::run_material::Capsule,preparation:RunPreparation},MismatchPcr{forward:Primer,reverse:Primer,options:crate::mutagenesis::PcrOptions},MismatchPcrCatalogDigest{forward:Primer,reverse:Primer,options:crate::mutagenesis::PcrOptions,selection:crate::catalog_digest::CatalogSelection,fragment_index:usize,repair:Option<crate::chemistry::EndRepairMethod>},OverlapExtensionPcr{request:crate::overlap_extension_pcr::Request},OverlapExtensionPcrCatalogDigest{request:crate::overlap_extension_pcr::Request,selection:crate::catalog_digest::CatalogSelection,fragment_index:usize,repair:Option<crate::chemistry::EndRepairMethod>},FromDigest{origin:crate::run_digest_history::Capsule,preparation:RunPreparation}}}
dto! {pub enum RunPreparation {WholeLinear,Digest{enzymes:Vec<String>,fragment_index:usize,repair:Option<crate::chemistry::EndRepairMethod>},CatalogDigest{selection:crate::catalog_digest::CatalogSelection,fragment_index:usize,repair:Option<crate::chemistry::EndRepairMethod>},ExactPcr{forward:Primer,reverse:Primer},PcrCatalogDigest{forward:Primer,reverse:Primer,selection:crate::catalog_digest::CatalogSelection,fragment_index:usize,repair:Option<crate::chemistry::EndRepairMethod>},MismatchPcr{forward:Primer,reverse:Primer,options:crate::mutagenesis::PcrOptions},MismatchPcrCatalogDigest{forward:Primer,reverse:Primer,options:crate::mutagenesis::PcrOptions,selection:crate::catalog_digest::CatalogSelection,fragment_index:usize,repair:Option<crate::chemistry::EndRepairMethod>},OverlapExtensionPcr{request:crate::overlap_extension_pcr::Request},OverlapExtensionPcrCatalogDigest{request:crate::overlap_extension_pcr::Request,selection:crate::catalog_digest::CatalogSelection,fragment_index:usize,repair:Option<crate::chemistry::EndRepairMethod>}}}
/// Literal catalog pin used by a preparation, including a flat saved-origin
/// child. This is structural lookup only, never replay or cleavage authority.
pub fn preparation_catalog_selection(
    preparation: &Preparation,
) -> Option<&crate::catalog_digest::CatalogSelection> {
    match preparation {
        Preparation::CatalogDigest { selection, .. }
        | Preparation::PcrCatalogDigest { selection, .. }
        | Preparation::MismatchPcrCatalogDigest { selection, .. }
        | Preparation::OverlapExtensionPcrCatalogDigest { selection, .. } => Some(selection),
        Preparation::FromRun { preparation, .. } | Preparation::FromDigest { preparation, .. } => {
            match preparation {
                RunPreparation::CatalogDigest { selection, .. }
                | RunPreparation::PcrCatalogDigest { selection, .. }
                | RunPreparation::MismatchPcrCatalogDigest { selection, .. }
                | RunPreparation::OverlapExtensionPcrCatalogDigest { selection, .. } => {
                    Some(selection)
                }
                _ => None,
            }
        }
        _ => None,
    }
}
/// Exact flat operation count, for host/preflight limits before any model copy.
/// It counts the declared graph only and does not establish reaction completion.
pub(crate) fn preparation_operation_count(preparation: &Preparation) -> usize {
    let repair = |repair: &Option<crate::chemistry::EndRepairMethod>| usize::from(repair.is_some());
    match preparation {
        Preparation::WholeLinear => 0,
        Preparation::ExactPcr { .. }
        | Preparation::MismatchPcr { .. }
        | Preparation::DeclaredDuplex { .. }
        | Preparation::AnnealedOligos { .. } => 1,
        Preparation::Digest { repair: r, .. } | Preparation::CatalogDigest { repair: r, .. } => {
            1 + repair(r)
        }
        Preparation::PcrCatalogDigest { repair: r, .. }
        | Preparation::MismatchPcrCatalogDigest { repair: r, .. } => 2 + repair(r),
        Preparation::OverlapExtensionPcr { .. } => 4,
        Preparation::OverlapExtensionPcrCatalogDigest { repair: r, .. } => 5 + repair(r),
        Preparation::FromRun { preparation, .. } | Preparation::FromDigest { preparation, .. } => {
            1 + match preparation {
                RunPreparation::WholeLinear => 0,
                RunPreparation::ExactPcr { .. } | RunPreparation::MismatchPcr { .. } => 1,
                RunPreparation::Digest { repair: r, .. }
                | RunPreparation::CatalogDigest { repair: r, .. } => 1 + repair(r),
                RunPreparation::PcrCatalogDigest { repair: r, .. }
                | RunPreparation::MismatchPcrCatalogDigest { repair: r, .. } => 2 + repair(r),
                RunPreparation::OverlapExtensionPcr { .. } => 4,
                RunPreparation::OverlapExtensionPcrCatalogDigest { repair: r, .. } => 5 + repair(r),
            }
        }
    }
}
impl RunPreparation {
    /// Bounded conversion to the existing preparation vocabulary; no recursion,
    /// chemistry restoration, or scientific computation.
    pub(crate) fn preparation(&self) -> Result<Preparation, Error> {
        match self {
            Self::OverlapExtensionPcrCatalogDigest {
                request,
                selection,
                fragment_index,
                repair,
            } => {
                overlap_extension_catalog_shape(request, selection, *fragment_index)?;
                Ok(Preparation::OverlapExtensionPcrCatalogDigest {
                    request: request.clone(),
                    selection: selection.clone(),
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
                fragment_index,
                repair,
            } => {
                selection.validate()?;
                if *fragment_index > 100_000 {
                    return Err(Error::TooLarge);
                }
                Ok(Preparation::CatalogDigest {
                    selection: selection.clone(),
                    fragment_index: *fragment_index,
                    repair: *repair,
                })
            }
            Self::MismatchPcrCatalogDigest {
                forward,
                reverse,
                options,
                selection,
                fragment_index,
                repair,
            } => {
                mismatch_pcr_catalog_shape(forward, reverse, options, selection, *fragment_index)?;
                Ok(Preparation::MismatchPcrCatalogDigest {
                    forward: forward.clone(),
                    reverse: reverse.clone(),
                    options: options.clone(),
                    selection: selection.clone(),
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
            Self::ExactPcr { forward, reverse } => {
                primer(forward)?;
                primer(reverse)?;
                Ok(Preparation::ExactPcr {
                    forward: forward.clone(),
                    reverse: reverse.clone(),
                })
            }
            Self::PcrCatalogDigest {
                forward,
                reverse,
                selection,
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
                    selection: selection.clone(),
                    fragment_index: *fragment_index,
                    repair: *repair,
                })
            }
        }
    }
}
pub(crate) fn from_run_shape(
    origin: &crate::run_material::Capsule,
    preparation: &RunPreparation,
) -> Result<(), Error> {
    origin.validate_documentary()?;
    preparation.preparation()?;
    Ok(())
}
pub(crate) fn from_digest_shape(
    origin: &crate::run_digest_history::Capsule,
    preparation: &RunPreparation,
) -> Result<(), Error> {
    origin.validate_documentary()?;
    preparation.preparation()?;
    Ok(())
}
dto! {pub struct FragmentSpec {pub role:String,pub input_id:String,pub preparation:Preparation,pub reverse:bool}}
dto! {pub struct Recipe {pub name:String,pub fragments:Vec<FragmentSpec>,pub overlaps:Vec<usize>,pub circular:bool,pub options:ProfileOptions}}
dto! {pub struct PreparedRow {pub row_index:usize,pub role:String,pub input_id:String,pub input_revision:String,pub preparation_step:Option<usize>,pub repair_step:Option<usize>,pub prepared_source:Source,pub reverse:bool}}
dto! {pub struct Compiled {pub plan:DesignPlan,pub rows:Vec<PreparedRow>,pub assembly_step:usize}}
dto! {pub struct Preview {pub context_revision:String,pub compiled:Compiled,pub verification:VerificationReport,
/// Fresh prepared documents in recipe order, before optional orientation.
/// On a failed preparation only the available leading row prefix is returned.
pub prepared:Vec<Molecule>,pub profile:Option<ProfileReport>,pub preparation_complete:bool}}
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
struct Hash {
    hash: blake3::Hasher,
    bytes: usize,
}
impl Write for Hash {
    fn write(&mut self, b: &[u8]) -> std::io::Result<usize> {
        self.bytes = self
            .bytes
            .checked_add(b.len())
            .filter(|n| *n <= MAX_BYTES)
            .ok_or_else(|| std::io::Error::other("bounded assembly recipe"))?;
        self.hash.update(b);
        Ok(b.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}
pub(crate) fn hash(v: &impl Serialize) -> Result<String, Error> {
    let mut w = Hash {
        hash: blake3::Hasher::new(),
        bytes: 0,
    };
    serde_json::to_writer(&mut w, v).map_err(|_| Error::TooLarge)?;
    Ok(w.hash.finalize().to_hex().to_string())
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
fn anneal_sources(top: &Molecule, bottom: &Molecule) -> Result<(), Error> {
    for source in [top, bottom] {
        if source.sequence.is_empty()
            || source.sequence.len() > crate::chemistry::MAX_ANNEAL_OLIGO_BASES
        {
            return Err(Error::TooLarge);
        }
        verification::declared_oligo(source)?;
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
    crate::run_material::material_capsules_preflight(steps.iter().filter_map(|step| match step {
        Operation::ReconstructRunProduct { origin, .. } => {
            Some(crate::run_material::MaterialOriginRef::Run(origin))
        }
        Operation::ReconstructDigestFragment { origin, .. } => {
            Some(crate::run_material::MaterialOriginRef::Digest(origin))
        }
        _ => None,
    }))?;
    let mut ids = BTreeSet::new();
    for input in inputs {
        token(&input.id, 128)?;
        crate::model::bounds::model_preflight(&[&input.molecule])?;
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
            crate::overlap_extension_pcr::extension_shape(options)?;
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
        if let Operation::ReconstructDigestFragment { template, origin } = step {
            origin.validate_documentary()?;
            let Source::Input(id) = template else {
                return Err(bad(
                    "digest reconstruction must consume its exact original input",
                ));
            };
            let input = inputs
                .iter()
                .find(|i| &i.id == id)
                .ok_or_else(|| bad("missing saved digest source"))?;
            if input.revision != origin.record.product.revision {
                return Err(bad("saved digest source revision changed"));
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
        | Operation::DeclareSourceTermini { template, .. }
        | Operation::ReconstructRunProduct { template, .. }
        | Operation::ReconstructDigestFragment { template, .. } => Some(template),
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
            | Operation::ReconstructDigestFragment { .. }
            | Operation::Pcr { .. }
            | Operation::PcrWithMismatches { .. }
            | Operation::Digest { .. }
            | Operation::CatalogDigest { .. }
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
            Operation::ReconstructRunProduct { .. } | Operation::ReconstructDigestFragment { .. }
        ) && *cursor < assembly
            && matches!(plan.steps.get(*cursor),Some(Operation::Pcr{template:Source::Step(n),..}|Operation::PcrWithMismatches{template:Source::Step(n),..}|Operation::Digest{template:Source::Step(n),..}|Operation::CatalogDigest{template:Source::Step(n),..}) if *n==first)
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
            && matches!(plan.steps.get(*cursor),Some(Operation::CatalogDigest{template:Source::Step(n),..})if *n==preceding)
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
                    Operation::Digest { .. } | Operation::CatalogDigest { .. }
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
    options: &crate::mutagenesis::PcrOptions,
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
    request: &crate::overlap_extension_pcr::Request,
    selection: &crate::catalog_digest::CatalogSelection,
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
    options: &crate::mutagenesis::PcrOptions,
    selection: &crate::catalog_digest::CatalogSelection,
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
    options: &crate::mutagenesis::PcrOptions,
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
    crate::run_material::material_capsules_preflight(fragments.iter().filter_map(|row| {
        match &row.preparation {
            Preparation::FromRun { origin, .. } => {
                Some(crate::run_material::MaterialOriginRef::Run(origin))
            }
            Preparation::FromDigest { origin, .. } => {
                Some(crate::run_material::MaterialOriginRef::Digest(origin))
            }
            _ => None,
        }
    }))?;
    let mut models = inputs.iter().map(|(_, m)| *m).collect::<Vec<_>>();
    if let Some(target) = target {
        models.push(target);
    }
    crate::model::bounds::model_preflight(&models)?;
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
    let catalog = crate::analysis::starter_enzymes();
    let pinned_catalog = if fragments
        .iter()
        .any(|f| preparation_catalog_selection(&f.preparation).is_some())
    {
        Some(crate::catalog_digest::builtin_catalog()?)
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
        } else if let Preparation::FromDigest {
            origin,
            preparation,
        } = &row.preparation
        {
            from_digest_shape(origin, preparation)?;
            if origin.record.product.revision != m.revision() {
                return Err(bad("saved digest origin source revision changed"));
            }
            effective = preparation.preparation()?;
            &effective
        } else {
            &row.preparation
        };
        let units = match preparation {
            Preparation::OverlapExtensionPcrCatalogDigest {
                request,
                selection,
                fragment_index,
                repair,
            } => {
                overlap_extension_catalog_shape(request, selection, *fragment_index)?;
                let (pcr_work, extra) = request.work(m.sequence.len())?;
                let enzymes =
                    selection.resolved_enzymes(pinned_catalog.as_ref().ok_or(Error::TooLarge)?)?;
                let width = enzymes.iter().try_fold(0usize, |sum, e| {
                    sum.checked_add(e.recognition.len() * 2)
                        .ok_or(Error::TooLarge)
                })?;
                let amplicon = m
                    .sequence
                    .len()
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
                role_bases = role_bases
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
                let (units, extra) = request.work(m.sequence.len())?;
                role_bases = role_bases
                    .checked_add(extra)
                    .filter(|n| *n <= MAX_ROLE_BASES)
                    .ok_or(Error::TooLarge)?;
                units
            }
            Preparation::FromRun { .. } | Preparation::FromDigest { .. } => {
                return Err(bad("nested material preparation is unsupported"));
            }
            Preparation::AnnealedOligos {
                bottom_input_id,
                options,
                candidate_index,
            } => {
                anneal_shape(bottom_input_id, options, *candidate_index)?;
                if bottom_input_id == &row.input_id {
                    return Err(bad(
                        "annealing requires two distinct original input identities",
                    ));
                }
                let bottom = lookup
                    .get(bottom_input_id.as_str())
                    .copied()
                    .ok_or_else(|| bad("unknown bottom oligo original input"))?;
                anneal_sources(m, bottom)?;
                role_bases = role_bases
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
                verification::declared_source_termini_preflight(m)?;
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
            Preparation::ExactPcr { forward, reverse } => {
                primer(forward)?;
                primer(reverse)?;
                m.sequence
                    .len()
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
                    mismatch_pcr_catalog_shape(
                        forward,
                        reverse,
                        options,
                        selection,
                        *fragment_index,
                    )?;
                } else {
                    primer(forward)?;
                    primer(reverse)?;
                }
                if *fragment_index > 100_000 {
                    return Err(Error::TooLarge);
                }
                let enzymes =
                    selection.resolved_enzymes(pinned_catalog.as_ref().ok_or(Error::TooLarge)?)?;
                let width = enzymes.iter().try_fold(0usize, |sum, e| {
                    sum.checked_add(e.recognition.len() * 2)
                        .ok_or(Error::TooLarge)
                })?;
                let amplicon = m
                    .sequence
                    .len()
                    .checked_add(forward.tail.len())
                    .and_then(|n| n.checked_add(reverse.tail.len()))
                    .ok_or(Error::TooLarge)?;
                role_bases = role_bases
                    .checked_add(forward.tail.len())
                    .and_then(|n| n.checked_add(reverse.tail.len()))
                    .filter(|n| *n <= MAX_ROLE_BASES)
                    .ok_or(Error::TooLarge)?;
                let pcr = if let Preparation::MismatchPcrCatalogDigest { options, .. } = preparation
                {
                    mismatch_pcr_work(m.sequence.len(), forward, reverse, options)?
                } else {
                    m.sequence
                        .len()
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
                let enzymes =
                    selection.resolved_enzymes(pinned_catalog.as_ref().ok_or(Error::TooLarge)?)?;
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
                    let enzyme=catalog.iter().find(|e|e.name.eq_ignore_ascii_case(name)).ok_or_else(||bad("digest enzyme is not supported by the existing typed Digest operation"))?;
                    width = width
                        .checked_add(enzyme.recognition.len() * 2)
                        .ok_or(Error::TooLarge)?;
                }
                m.sequence.len().checked_mul(width).ok_or(Error::TooLarge)?
            }
        };
        work = work
            .checked_add(units)
            .filter(|n| *n <= MAX_PREPARATION_WORK)
            .ok_or(Error::TooLarge)?;
    }
    hash(&(inputs, target, name, fragments))?;
    Ok(())
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
    let (steps, rows) = preparation_steps(inputs, fragments)?;
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
        } else if let Preparation::FromDigest {
            origin,
            preparation,
        } = &row.preparation
        {
            let first = steps.len();
            steps.push(Operation::ReconstructDigestFragment {
                template: original_input,
                origin: origin.clone(),
            });
            effective = preparation.preparation()?;
            (Source::Step(first), Some(first), &effective)
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
                fragment_index,
                repair,
                ..
            } = preparation
            {
                let index = steps.len();
                steps.push(Operation::CatalogDigest {
                    template: prepared_source,
                    selection: selection.clone(),
                    fragment_index: *fragment_index,
                });
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
            Preparation::FromRun { .. } | Preparation::FromDigest { .. } => {
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
                fragment_index,
                ..
            } => Some(Operation::CatalogDigest {
                template: input.clone(),
                selection: selection.clone(),
                fragment_index: *fragment_index,
            }),
            Preparation::Digest {
                enzymes,
                fragment_index,
                ..
            } => Some(Operation::Digest {
                template: input.clone(),
                enzymes: enzymes.clone(),
                fragment_index: *fragment_index,
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
            Preparation::ExactPcr { forward, reverse }
            | Preparation::PcrCatalogDigest {
                forward, reverse, ..
            } => Some(Operation::Pcr {
                template: input.clone(),
                forward: forward.clone(),
                reverse: reverse.clone(),
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
            fragment_index,
            ..
        }
        | Preparation::MismatchPcrCatalogDigest {
            selection,
            fragment_index,
            ..
        } = preparation
        {
            let index = steps.len();
            steps.push(Operation::CatalogDigest {
                template: prepared_source,
                selection: selection.clone(),
                fragment_index: *fragment_index,
            });
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
    let prepared = compile_preparations(inputs, target, &recipe.name, &recipe.fragments)?;
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
    let mut verification = verification::verify(&compiled.plan, target);
    let mut available = Vec::with_capacity(compiled.rows.len());
    for row in &compiled.rows {
        let m = match &row.prepared_source {
            Source::Input(id) => compiled
                .plan
                .inputs
                .iter()
                .find(|i| &i.id == id)
                .map(|i| &i.molecule),
            Source::Step(i) => verification.products.get(*i),
        };
        if let Some(m) = m {
            available.push(m)
        } else {
            break;
        }
    }
    let preparation_complete = available.len() == compiled.rows.len();
    crate::model::bounds::model_preflight(&available)?;
    let mut bases = 0usize;
    for m in &available {
        bases = bases
            .checked_add(m.sequence.len())
            .filter(|n| *n <= MAX_ROLE_BASES)
            .ok_or(Error::TooLarge)?;
    }
    hash(&available)?;
    let prepared = available.into_iter().cloned().collect::<Vec<_>>();
    // An advisory junction report remains useful after a target mismatch. Physical
    // validity/blocked/staggered substrates remain the verifier's separate result.
    let profile = if preparation_complete {
        match overlap_profiles::check_oriented_fragments(
            &prepared,
            &compiled.rows.iter().map(|r| r.reverse).collect::<Vec<_>>(),
            &recipe.overlaps,
            recipe.circular,
            &recipe.options,
        ) {
            Ok(report) => Some(report),
            Err(error) => {
                verification.findings.push(format!(
                    "Recipe advisory junction profile withheld: {error}"
                ));
                None
            }
        }
    } else {
        None
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
    policy: &crate::homology_assembly::Options,
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
    policy: &crate::homology_assembly::Options,
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
    policy: &crate::homology_assembly::Options,
) -> Result<String, Error> {
    homology_preflight(inputs, target, recipe, policy)?;
    bound_revisions(inputs, target)?;
    hash(&(
        ENGINE_VERSION,
        crate::homology_assembly::ENGINE_VERSION,
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
    policy: &crate::homology_assembly::Options,
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
    policy: &crate::homology_assembly::Options,
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
    policy: &crate::homology_annealing::Options,
) -> Result<(), Error> {
    policy.validate_shape(recipe.fragments.len())?;
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
    policy: &crate::homology_annealing::Options,
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
    policy: &crate::homology_annealing::Options,
) -> Result<String, Error> {
    annealing_preflight(inputs, target, recipe, policy)?;
    bound_revisions(inputs, target)?;
    hash(&(
        ENGINE_VERSION,
        crate::homology_annealing::ENGINE_VERSION,
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
    policy: &crate::homology_annealing::Options,
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
    policy: &crate::homology_annealing::Options,
) -> Result<Preview, Error> {
    let context = annealing_preview_revision(inputs, target, recipe, policy)?;
    let compiled = compile_annealing(inputs, target, recipe, policy)?;
    finish_preview(context, compiled, target, recipe)
}

impl Preview {
    /// Readonly shape/binding checks only. No reconstruction or authority restoration.
    pub fn validate(&self) -> Result<(), Error> {
        if !hex(&self.context_revision)
            || self.compiled.rows.is_empty()
            || self.compiled.rows.len() > MAX_ROWS
            || self.prepared.len() > self.compiled.rows.len()
            || self.compiled.plan.inputs.len() > MAX_ROWS
            || self.verification.products.len() > MAX_STEPS_PER_ROW * MAX_ROWS + 1
            || self.verification.step_reports.len() > MAX_STEPS_PER_ROW * MAX_ROWS + 1
            || self.verification.annotation_reports.len() > MAX_STEPS_PER_ROW * MAX_ROWS + 1
            || self.compiled.plan.steps.len() > MAX_STEPS_PER_ROW * MAX_ROWS + 1
            || self.compiled.assembly_step.checked_add(1) != Some(self.compiled.plan.steps.len())
        {
            return Err(bad("invalid preview shape"));
        }
        if self.preparation_complete != (self.prepared.len() == self.compiled.rows.len())
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
            .chain(self.prepared.iter())
            .chain(self.verification.products.iter())
            .collect::<Vec<_>>();
        for model in models {
            crate::model::bounds::model_preflight(&[model])?;
        }
        for report in &self.verification.annotation_reports {
            if report.step_index >= self.compiled.plan.steps.len()
                || report.evidence.sources.len() > MAX_ROWS
            {
                return Err(Error::TooLarge);
            }
            let mut budget = crate::model::bounds::ModelBudget::default();
            budget.features(&report.evidence.features)?;
            for role in &report.evidence.sources {
                if role.unmapped.len() > 4096 {
                    return Err(Error::TooLarge);
                }
                for unmapped in &role.unmapped {
                    budget.features(std::slice::from_ref(&unmapped.feature))?;
                }
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
                if self.verification.products.len() == self.compiled.plan.steps.len() {
                    let actual = self
                        .verification
                        .products
                        .last()
                        .ok_or_else(|| bad("missing annealing diagnostic reference"))?;
                    let evidence=self.verification.physical_processing_reports.iter().find_map(|step|if step.step_index==self.compiled.assembly_step{match &step.evidence{crate::verification::PhysicalProcessingEvidence::HomologyAnnealing(report)=>Some(report),_=>None}}else{None}).ok_or_else(||bad("missing exact annealing strand evidence"))?;
                    evidence.validate_documentary()?;
                    if hash(policy)? != hash(&evidence.options)?
                        || evidence.preparations.len() != fragments.len()
                        || evidence.reference != actual.sequence
                        || policy.topology != actual.topology
                        || self.verification.status
                            == crate::verification::Status::ConsistentUnderAssumptions
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
                            | Operation::ReconstructDigestFragment { .. }
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
                    .map(hash)
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
