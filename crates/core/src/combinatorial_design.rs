//! Explicit source alternatives over existing preparation compilers.
//! This is a design-choice layer, independent of reaction-product enumeration.
use crate::{
    assembly_edit_intent, assembly_intent, assembly_junction_edit,
    assembly_recipe::{self, OverlapFinishing},
    assembly_substitution_intent,
    chemistry::EndTreatment,
    golden_gate::{GoldenGateEnzyme, GoldenGateOptions},
    ligation_recipe,
    overlap_profiles::ProfileOptions,
    run_lineage::{EventBinding, InputEventBinding},
    tailing_recipe,
    verification::PlanInput,
    workspace::Workspace,
};
use collagene_model::{Error, Molecule, Topology};
use serde::{Deserialize, Serialize};
use std::{collections::BTreeSet, io};

const VERSION: &str = "combinatorial-source-design-3";
pub const MAX_BINS: usize = 12;
pub const MAX_CHOICES: usize = 256;
pub const MAX_SELECTED: usize = 256;
pub const MAX_PAGE: usize = 64;
const MAX_BYTES: usize = 32 * 1024 * 1024;
macro_rules! dto {($($t:tt)*)=>{#[derive(Clone,Debug,Serialize,Deserialize)]#[cfg_attr(feature="json-schema",derive(schemars::JsonSchema))]#[serde(deny_unknown_fields)]$($t)*};}
dto! {pub struct Choice {
    pub id:String,pub label:String,pub source:EventBinding,
    /// Explicit independent contribution domain. Stored native strand offsets
    /// are distinct from original-reference coordinates and physical orientation.
    pub selector:assembly_intent::Selector,pub preparation_reverse:bool,
    /// Unfinished alternatives can be retained; an explicitly chosen one needs preparation.
    pub preparation:Option<assembly_recipe::Preparation>,
    pub treatment:Option<EndTreatment>,
    #[serde(deserialize_with="assembly_recipe::required_option")]#[cfg_attr(feature="json-schema",schemars(required))]pub generation:Option<crate::modified_duplex::GenerationAssumptions>,
}}
dto! {pub struct Bin {pub id:String,pub label:String,pub choices:Vec<Choice>}}
dto! {pub struct TailingAction {pub bin_id:String,pub action:tailing_recipe::Action}}
dto! {pub enum Method {
    Tailing{profile:tailing_recipe::Profile,actions:Vec<TailingAction>,allow_nicked:bool},
    Overlap{overlaps:Vec<usize>,options:ProfileOptions,finishing:OverlapFinishing},
    Ligation{allow_nicked:bool},
    /// Same existing pinned-site geometry requirement; no kinetic-mixture inference.
    OrderedTypeIis{enzyme:GoldenGateEnzyme,allow_nicked:bool},
    /// Physical reaction substrates are independent of original target-contribution
    /// coordinates. An absent pin is an editable inspection draft, never candidate zero.
    GoldenGatePool{options:GoldenGateOptions,substrates:Vec<PoolSubstrate>,selection:Option<PoolSelection>},
    /// Whole original reference parents; contribution intervals only define independent intent.
    MultiSiteReference{request:crate::multisite_recombination::Request,selection:Option<MultiSiteSelection>},
}}
dto! {pub enum PoolSubstrateKind {WholeSource,Prepared}}
dto! {pub struct PoolSubstrate {pub bin_id:String,pub kind:PoolSubstrateKind}}
// Documentary user-selection intent, never imported scientific authority.
// The outer revision also pins equal-DNA distinct workspace events and geometry.
dto! {pub struct PoolSelection {pub selection:crate::prepared_golden_gate::Selection,pub source_revision:String}}
dto! {pub struct MultiSiteSelection {pub selection:crate::recombination::Selection,pub source_revision:String}}
dto! {pub struct MultiSiteInspection {pub context_revision:String,pub source_revision:String,pub batch:crate::multisite_recombination::Batch}}
dto! {pub struct LibraryMaterialBinding {pub event:EventBinding,pub origin:crate::run_material::MaterialOrigin}}
dto! {pub struct Design {pub id:String,pub name:String,pub circular:bool,pub bins:Vec<Bin>,pub method:Method,pub material_bindings:Vec<LibraryMaterialBinding>}}
dto! {pub struct ChoiceSelection {pub bin_id:String,pub choice_id:String}}
dto! {pub struct Combination {pub id:String,pub choices:Vec<ChoiceSelection>}}
dto! {pub struct MatrixPage {pub total:u128,pub first:u128,pub combinations:Vec<Combination>}}
dto! {pub enum TargetDeclaration {
    /// An explicitly selected independent retained reference event.
    Existing(EventBinding),
    /// Literal original contribution coordinates + optional declared junction edits.
    /// Uses existing independent intent builder, never the reconstructed product.
    OriginalRegions{junction_edits:Vec<assembly_junction_edit::Edit>,substitutions:Vec<assembly_substitution_intent::Substitution>,edits:Vec<assembly_edit_intent::Edit>},
}}
dto! {pub struct PreparationOverride {pub bin_id:String,pub preparation_reverse:bool,pub preparation:Option<assembly_recipe::Preparation>,pub treatment:Option<EndTreatment>,#[serde(deserialize_with="assembly_recipe::required_option")]#[cfg_attr(feature="json-schema",schemars(required))]pub generation:Option<crate::modified_duplex::GenerationAssumptions>}}
dto! {pub struct ConstructOverrides {pub preparations:Vec<PreparationOverride>,pub method:Option<Method>,pub circular:Option<bool>}}
dto! {pub struct SelectedConstruct {
    pub id:String,pub name:String,pub choices:Vec<ChoiceSelection>,pub target:TargetDeclaration,
    #[serde(default)] pub overrides:Option<ConstructOverrides>,
}}
dto! {pub struct SelectedMatrix {pub design:Design,pub constructs:Vec<SelectedConstruct>}}
#[derive(Clone, Debug)]
pub enum IndependentIntent {
    Regions(assembly_intent::Report),
    Junctions(assembly_junction_edit::Intent),
    Substitutions(assembly_substitution_intent::Report),
    Edits(assembly_edit_intent::Report),
}
#[derive(Clone, Debug)]
pub enum Recipe {
    Tailing {
        recipe: tailing_recipe::Recipe,
    },
    Overlap {
        recipe: assembly_recipe::Recipe,
        finishing: OverlapFinishing,
    },
    Ligation {
        recipe: ligation_recipe::Recipe,
        cleavage_policy: assembly_recipe::CleavagePolicy,
    },
    GoldenGatePool {
        request: crate::golden_gate_reaction::Request,
        selection: Option<PoolSelection>,
        source_revision: String,
    },
    MultiSiteReference {
        request: crate::multisite_recombination::Request,
        selection: Option<MultiSiteSelection>,
        source_revision: String,
        sources: Vec<crate::verification::TopoSource>,
    },
}
pub struct Resolved {
    pub construct_id: String,
    pub construct_name: String,
    pub combination_id: String,
    pub declaration_revision: String,
    pub inputs: Vec<PlanInput>,
    pub bindings: Vec<InputEventBinding>,
    pub target: Molecule,
    pub target_binding: Option<EventBinding>,
    pub intent: Option<IndependentIntent>,
    pub recipe: Recipe,
    preparation_work: usize,
}
impl Resolved {
    /// Actual declared preparation/search cost resolved in the same native
    /// source context as independent intent. This scalar grants no chemistry.
    pub fn declared_work(&self) -> usize {
        self.preparation_work
    }
}
pub struct Compiled {
    pub resolved: Resolved,
    pub plan: crate::verification::DesignPlan,
}
pub enum Preview {
    Tailing(tailing_recipe::Preview),
    Overlap(assembly_recipe::Preview),
    Ligation(ligation_recipe::Preview),
    GoldenGatePool(PoolPreview),
    MultiSiteReference(PoolPreview),
}
pub struct PoolPreview {
    pub plan: crate::verification::DesignPlan,
    pub verification: crate::verification::VerificationReport,
    pub preparation_complete: bool,
}
pub struct Evaluation {
    pub resolved: Resolved,
    pub preview: Preview,
}
fn bad(message: &str) -> Error {
    Error::Format(format!("Combinatorial design: {message}"))
}
fn identifier(value: &str) -> bool {
    !value.is_empty()
        && value.len() <= 48
        && value
            .bytes()
            .all(|b| b.is_ascii_alphanumeric() || matches!(b, b'-' | b'_' | b'.'))
}
fn label(value: &str, limit: usize) -> bool {
    !value.trim().is_empty() && value.len() <= limit && !value.chars().any(char::is_control)
}
fn digest(value: &str) -> bool {
    value.len() == 64
        && value
            .bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}
fn binding(value: &EventBinding) -> Result<(), Error> {
    if value.document_id.is_empty()
        || value.document_id.len() > 128
        || !value
            .document_id
            .bytes()
            .all(|b| b.is_ascii_alphanumeric() || matches!(b, b'_' | b'-' | b'.' | b':'))
        || !digest(&value.revision)
        || !value.event_id.strip_prefix("event-").is_some_and(digest)
    {
        return Err(bad("invalid exact source binding"));
    }
    Ok(())
}
struct Counter(usize);
impl io::Write for Counter {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        self.0 = self
            .0
            .checked_add(bytes.len())
            .filter(|n| *n <= MAX_BYTES)
            .ok_or_else(|| io::Error::other("combinatorial design byte limit"))?;
        Ok(bytes.len())
    }
    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}
fn count(value: &impl Serialize) -> Result<(), Error> {
    serde_json::to_writer(&mut Counter(0), value).map_err(|_| Error::TooLarge)
}
/// Stable operation IDs preserve design/construct/bin/choice identity even when
/// labels, DNA contents or display indices are equal. Events remain separate pins.
pub fn source_role_ids(
    design_id: &str,
    construct_id: &str,
    bin_id: &str,
    choice_id: &str,
) -> Result<(String, String), Error> {
    if [design_id, construct_id, bin_id, choice_id]
        .iter()
        .any(|id| !identifier(id))
    {
        return Err(bad("invalid stable source/role identity"));
    }
    let hash = assembly_recipe::hash(&(
        VERSION,
        "source-role",
        design_id,
        construct_id,
        bin_id,
        choice_id,
    ))?;
    Ok((format!("source-{hash}"), format!("fragment-{bin_id}")))
}
/// Mixed-radix dimensions shared by typed scientific designs and unfinished UI
/// drafts. Counting/paging itself grants no source/preparation validity.
pub fn cartesian_size(counts: &[usize]) -> Result<u128, Error> {
    if counts.is_empty() || counts.len() > MAX_BINS || counts.iter().any(|n| *n == 0 || *n > 64) {
        return Err(bad("invalid alternative-bin dimensions"));
    }
    let total_choices = counts.iter().try_fold(0usize, |n, count| {
        n.checked_add(*count).ok_or(Error::TooLarge)
    })?;
    if total_choices > MAX_CHOICES {
        return Err(Error::TooLarge);
    }
    counts.iter().try_fold(1u128, |total, count| {
        total.checked_mul(*count as u128).ok_or(Error::TooLarge)
    })
}
pub fn cartesian_page(
    counts: &[usize],
    first: u128,
    limit: usize,
) -> Result<Vec<Vec<usize>>, Error> {
    let total = cartesian_size(counts)?;
    if first > total || limit == 0 || limit > MAX_PAGE {
        return Err(bad("invalid bounded matrix page"));
    }
    let end = first
        .checked_add(limit as u128)
        .unwrap_or(u128::MAX)
        .min(total);
    let mut rows = Vec::with_capacity((end - first) as usize);
    for rank in first..end {
        let mut remainder = rank;
        let mut indices = vec![0usize; counts.len()];
        for (i, count) in counts.iter().enumerate().rev() {
            indices[i] = (remainder % *count as u128) as usize;
            remainder /= *count as u128;
        }
        rows.push(indices);
    }
    Ok(rows)
}
fn preparation_origin(
    preparation: &assembly_recipe::Preparation,
) -> Option<crate::run_material::NativeOriginRef<'_>> {
    match preparation {
        assembly_recipe::Preparation::FromRun { origin, .. } => {
            Some(crate::run_material::NativeOriginRef::Run(origin))
        }
        assembly_recipe::Preparation::FromDerived { origin, .. } => {
            Some(crate::run_material::NativeOriginRef::Derived(origin))
        }
        _ => None,
    }
}
/// Shallow allocation limits only; source/literal/conflict semantics are assessed
/// per selected row so one invalid row does not suppress the rest of a batch.
pub(crate) fn declaration_bounds(target: &TargetDeclaration) -> Result<(), Error> {
    if let TargetDeclaration::OriginalRegions {
        junction_edits,
        substitutions,
        edits,
    } = target
    {
        if junction_edits.len() > MAX_BINS
            || substitutions.len() > 256
            || edits.len() > 256
            || junction_edits.iter().any(|e| e.inserted.len() > 200)
            || substitutions
                .iter()
                .any(|e| e.role.len() > 128 || e.expected.len() > 1 || e.replacement.len() > 1)
            || edits.iter().any(|e| e.role.len() > 128)
        {
            return Err(Error::TooLarge);
        }
        let mut bytes = 0usize;
        for e in edits {
            bytes = bytes
                .checked_add(e.expected.len())
                .and_then(|n| n.checked_add(e.replacement.len()))
                .filter(|n| *n <= 100_000)
                .ok_or(Error::TooLarge)?;
        }
    }
    Ok(())
}
fn material_binding<'a>(d: &'a Design, event: &EventBinding) -> Option<&'a LibraryMaterialBinding> {
    d.material_bindings.iter().find(|b| &b.event == event)
}
pub(crate) fn preparation_binding(
    d: &Design,
    event: &EventBinding,
    p: &assembly_recipe::Preparation,
) -> Result<(), Error> {
    let claimed = match p {
        assembly_recipe::Preparation::FromRun { origin, .. } => {
            Some((&origin.selection.event, assembly_recipe::hash(origin)?))
        }
        assembly_recipe::Preparation::FromDerived { origin, .. } => {
            Some((&origin.product, assembly_recipe::hash(origin)?))
        }
        _ => None,
    };
    if let Some((pin, hash)) = claimed {
        let actual = material_binding(d, event).ok_or_else(|| {
            bad("saved preparation requires its explicit material source binding")
        })?;
        let actual_hash = match &actual.origin {
            crate::run_material::MaterialOrigin::Run(o) => assembly_recipe::hash(o)?,
            crate::run_material::MaterialOrigin::Derived(o) => assembly_recipe::hash(o)?,
        };
        if pin != event || hash != actual_hash {
            return Err(bad(
                "saved preparation differs from declared exact material origin",
            ));
        }
    }
    Ok(())
}
/// Borrowed exact documentary ownership; no material replay or target construction.
/// Material bounds/state are freshly checked by the shared intent resolver on use.
pub fn source_preflight(design: &Design, workspace: &Workspace) -> Result<(), Error> {
    validate_design(design)?;
    workspace_guard(workspace)?;
    crate::run_material::validate_native_origins(
        workspace,
        design.material_bindings.iter().map(|b| match &b.origin {
            crate::run_material::MaterialOrigin::Run(o) => {
                crate::run_material::NativeOriginRef::Run(o)
            }
            crate::run_material::MaterialOrigin::Derived(o) => {
                crate::run_material::NativeOriginRef::Derived(o)
            }
        }),
    )?;
    for c in design.bins.iter().flat_map(|b| &b.choices) {
        let m = resolve_event(workspace, &c.source)?;
        if m.sequence.alphabet() != collagene_model::Alphabet::Dna {
            return Err(bad("DNA source required"));
        }
        if let assembly_intent::Selector::ReferenceRange { start, end, .. } = &c.selector {
            let n = m.sequence.len();
            if *start >= n
                || end.checked_sub(*start).is_none_or(|k| k == 0 || k > n)
                || (m.topology == Topology::Linear && *end > n)
            {
                return Err(bad("original reference contribution outside source"));
            }
        }
    }
    Ok(())
}
fn preparation_shape(preparation: &assembly_recipe::Preparation) -> Result<(), Error> {
    match preparation {
        assembly_recipe::Preparation::FromRun {
            origin,
            preparation,
        } => assembly_recipe::from_run_shape(origin, preparation)?,
        assembly_recipe::Preparation::FromDerived {
            origin,
            preparation,
        } => assembly_recipe::from_derived_shape(origin, preparation)?,
        assembly_recipe::Preparation::AnnealedOligos { .. } => {
            return Err(bad(
                "paired oligo preparation requires separately bound original inputs and is not supported by library projection",
            ));
        }
        assembly_recipe::Preparation::ExactPcr {
            forward, reverse, ..
        }
        | assembly_recipe::Preparation::PcrCatalogDigest {
            forward, reverse, ..
        } => {
            for primer in [forward, reverse] {
                if primer.name.len() > 256
                    || primer.annealing.len() > 200
                    || primer.tail.len() > 200
                    || primer.annealing.len().saturating_add(primer.tail.len())
                        > crate::primers::MAX_OLIGO
                {
                    return Err(Error::TooLarge);
                }
                if !label(&primer.name, 256)
                    || primer.name.trim() != primer.name
                    || primer.annealing.alphabet() != collagene_model::Alphabet::Dna
                {
                    return Err(bad("invalid explicit primer"));
                }
                let parsed = crate::analysis::Primer::new(
                    &primer.name,
                    primer.annealing.as_str(),
                    &primer.tail,
                )?;
                if parsed.tail != primer.tail || parsed.annealing != primer.annealing {
                    return Err(bad("typed preparation must retain canonical primer DNA"));
                }
            }
        }
        assembly_recipe::Preparation::MismatchPcr {
            forward,
            reverse,
            options,
        } => assembly_recipe::mismatch_pcr_shape(forward, reverse, options)?,
        assembly_recipe::Preparation::MismatchPcrCatalogDigest {
            forward,
            reverse,
            options,
            selection,
            fragment_index,
            ..
        } => assembly_recipe::mismatch_pcr_catalog_shape(
            forward,
            reverse,
            options,
            selection,
            *fragment_index,
        )?,
        assembly_recipe::Preparation::Digest { enzymes, .. } => {
            if enzymes.len() > 64 || enzymes.iter().any(|e| e.len() > 256) {
                return Err(Error::TooLarge);
            }
        }
        _ => {}
    }
    if let Some(selection) = assembly_recipe::preparation_catalog_selection(preparation) {
        selection.validate()?;
    }

    Ok(())
}
fn method_shape(method: &Method, bins: usize, circular: bool) -> Result<(), Error> {
    if let Method::Tailing {
        profile, actions, ..
    } = method
    {
        if actions.len() != bins || actions.len() > MAX_BINS {
            return Err(bad("declare exactly one explicit tail action per bin"));
        }
        let mut ids = BTreeSet::new();
        for row in actions {
            if !identifier(&row.bin_id) || !ids.insert(&row.bin_id) {
                return Err(bad("invalid or duplicate tail-action bin"));
            }
            if let tailing_recipe::Action::AddTail { base } = row.action {
                use crate::tailing::TailingBase;
                let allowed = match profile {
                    tailing_recipe::Profile::Ta => matches!(base, TailingBase::A | TailingBase::T),
                    tailing_recipe::Profile::Gc => matches!(base, TailingBase::G | TailingBase::C),
                };
                if !allowed {
                    return Err(bad("added base differs from the declared TA/GC profile"));
                }
            }
        }
    }
    if let Method::Overlap {
        overlaps,
        options,
        finishing,
    } = method
    {
        if overlaps.len() != bins.saturating_sub(1) + usize::from(circular)
            || overlaps.iter().any(|n| *n == 0 || *n > 10_000)
        {
            return Err(bad(
                "declare one positive overlap length per effective junction",
            ));
        }
        finishing.validate_geometry(bins, overlaps, circular, options)?;
    }
    if let Method::MultiSiteReference { request, selection } = method {
        if !(2..=5).contains(&bins) {
            return Err(bad(
                "multi-site reference exchange requires two to five original parent bins",
            ));
        }
        crate::multisite_recombination::request_preflight(request)?;
        if let Some(pin) = selection {
            if !digest(&pin.source_revision)
                || !digest(&pin.selection.batch_revision)
                || pin.selection.candidate_id >= 256
                || pin.selection.component_id >= 32
            {
                return Err(bad("invalid explicit multi-site candidate/component pin"));
            }
        }
    }
    if let Method::GoldenGatePool {
        options,
        substrates,
        selection,
    } = method
    {
        if substrates.len() != bins
            || options.circular != circular
            || options.min_fragments == 0
            || options.min_fragments > options.max_fragments
            || options.max_fragments > 32
            || options.max_products == 0
            || options.max_products > 1000
            || options.max_states == 0
            || options.max_states > 1_000_000
        {
            return Err(bad(
                "declare one substrate per bin and bounded reaction options matching construct topology",
            ));
        }
        let mut ids = BTreeSet::new();
        if substrates
            .iter()
            .any(|row| !identifier(&row.bin_id) || !ids.insert(&row.bin_id))
        {
            return Err(bad("invalid or duplicate reaction substrate bin"));
        }
        if let Some(pin) = selection {
            if !digest(&pin.source_revision)
                || !digest(&pin.selection.context_revision)
                || !digest(&pin.selection.directed_sequence_key)
                || !digest(&pin.selection.physical_sequence_key)
                || pin.selection.candidate_index >= options.max_products
            {
                return Err(bad("invalid bounded explicit reaction candidate pin"));
            }
        }
    }
    Ok(())
}
fn pool_bins(method: &Method, bins: &[Bin]) -> Result<(), Error> {
    if let Method::Tailing { actions, .. } = method {
        if actions
            .iter()
            .any(|row| !bins.iter().any(|bin| bin.id == row.bin_id))
        {
            return Err(bad("tail action references an unknown bin"));
        }
    }
    if let Method::GoldenGatePool { substrates, .. } = method {
        if substrates
            .iter()
            .any(|row| !bins.iter().any(|bin| bin.id == row.bin_id))
        {
            return Err(bad("reaction substrate references an unknown bin"));
        }
    }
    Ok(())
}
fn effective<'a>(design: &'a Design, construct: &'a SelectedConstruct) -> (&'a Method, bool) {
    let method = construct
        .overrides
        .as_ref()
        .and_then(|o| o.method.as_ref())
        .unwrap_or(&design.method);
    let circular = construct
        .overrides
        .as_ref()
        .and_then(|o| o.circular)
        .unwrap_or(design.circular);
    (method, circular)
}
fn preparation_override<'a>(
    construct: &'a SelectedConstruct,
    bin_id: &str,
) -> Option<&'a PreparationOverride> {
    construct
        .overrides
        .as_ref()?
        .preparations
        .iter()
        .find(|p| p.bin_id == bin_id)
}
pub fn validate_design(design: &Design) -> Result<(), Error> {
    if !identifier(&design.id)
        || !label(&design.name, 256)
        || design.bins.is_empty()
        || design.bins.len() > MAX_BINS
        || design.material_bindings.len() > MAX_CHOICES
    {
        return Err(bad("invalid design identity or bin count"));
    }
    let mut total = 0usize;
    for bin in &design.bins {
        if bin.choices.len() > 64 {
            return Err(Error::TooLarge);
        }
        total = total
            .checked_add(bin.choices.len())
            .filter(|n| *n <= MAX_CHOICES)
            .ok_or(Error::TooLarge)?;
    }
    crate::run_material::native_origins_preflight(
        design
            .bins
            .iter()
            .flat_map(|b| &b.choices)
            .filter_map(|c| c.preparation.as_ref())
            .filter_map(preparation_origin)
            .chain(design.material_bindings.iter().map(|b| match &b.origin {
                crate::run_material::MaterialOrigin::Run(o) => {
                    crate::run_material::NativeOriginRef::Run(o)
                }
                crate::run_material::MaterialOrigin::Derived(o) => {
                    crate::run_material::NativeOriginRef::Derived(o)
                }
            })),
    )?;
    let mut sources = BTreeSet::new();
    for b in &design.material_bindings {
        binding(&b.event)?;
        if !sources.insert((&b.event.document_id, &b.event.event_id, &b.event.revision))
            || b.origin.event() != &b.event
            || !design
                .bins
                .iter()
                .flat_map(|b| &b.choices)
                .any(|c| c.source == b.event)
        {
            return Err(bad(
                "duplicate, unused or mismatched exact native material binding",
            ));
        }
    }
    let (mut bin_ids, mut choice_ids) = (BTreeSet::new(), BTreeSet::new());
    let mut count_choices = 0usize;
    for bin in &design.bins {
        if !identifier(&bin.id)
            || !bin_ids.insert(&bin.id)
            || !label(&bin.label, 80)
            || bin.choices.is_empty()
            || bin.choices.len() > 64
        {
            return Err(bad("invalid/duplicate bin identity or choice count"));
        }
        count_choices = count_choices
            .checked_add(bin.choices.len())
            .filter(|n| *n <= MAX_CHOICES)
            .ok_or(Error::TooLarge)?;
        for choice in &bin.choices {
            if !identifier(&choice.id)
                || !choice_ids.insert(&choice.id)
                || !label(&choice.label, 256)
                || choice.selector.length()? > collagene_model::MAX_SEQUENCE_LEN
            {
                return Err(bad(
                    "invalid/duplicate choice identity or original contribution",
                ));
            }
            binding(&choice.source)?;
            if let assembly_intent::Selector::MaterialStrandRange {
                expected_state_revision,
                ..
            } = &choice.selector
            {
                if !digest(expected_state_revision)
                    || material_binding(design, &choice.source).is_none()
                {
                    return Err(bad(
                        "actual material contribution requires exact state and native source binding",
                    ));
                }
            }
            if material_binding(design, &choice.source).is_some() && choice.generation.is_some() {
                return Err(bad(
                    "native source binding cannot receive an original-generation reset",
                ));
            }
            if let Some(p) = &choice.preparation {
                preparation_binding(design, &choice.source, p)?;
            }

            if let Some(preparation) = &choice.preparation {
                preparation_shape(preparation)?;
            }
        }
    }
    method_shape(&design.method, design.bins.len(), design.circular)?;
    pool_bins(&design.method, &design.bins)?;
    count(design)
}
fn selected<'a>(
    design: &'a Design,
    choices: &[ChoiceSelection],
) -> Result<Vec<(&'a Bin, &'a Choice)>, Error> {
    if choices.len() != design.bins.len() {
        return Err(bad("choose exactly one alternative per ordered bin"));
    }
    let mut seen = BTreeSet::new();
    for choice in choices {
        if !seen.insert(&choice.bin_id) {
            return Err(bad("duplicate bin in selected construct"));
        }
    }
    design
        .bins
        .iter()
        .map(|bin| {
            let selection = choices
                .iter()
                .find(|selection| selection.bin_id == bin.id)
                .ok_or_else(|| bad("selected construct is missing a bin"))?;
            let choice = bin
                .choices
                .iter()
                .find(|choice| choice.id == selection.choice_id)
                .ok_or_else(|| bad("selected choice does not belong to its bin"))?;
            Ok((bin, choice))
        })
        .collect()
}
pub(crate) fn combination_key(
    design: &Design,
    choices: &[ChoiceSelection],
) -> Result<String, Error> {
    let canonical = selected(design, choices)?
        .iter()
        .map(|(bin, choice)| (&bin.id, &choice.id))
        .collect::<Vec<_>>();
    Ok(format!(
        "combination-{}",
        assembly_recipe::hash(&(VERSION, &design.id, canonical))?
    ))
}
pub fn combination_id(design: &Design, choices: &[ChoiceSelection]) -> Result<String, Error> {
    validate_design(design)?;
    combination_key(design, choices)
}
fn construct_shape(design: &Design, construct: &SelectedConstruct) -> Result<(), Error> {
    if !identifier(&construct.id) || !label(&construct.name, 256) {
        return Err(bad("invalid selected construct identity"));
    }
    selected(design, &construct.choices)?;
    let (method, circular) = effective(design, construct);
    method_shape(method, design.bins.len(), circular)?;
    pool_bins(method, &design.bins)?;
    if matches!(
        method,
        Method::GoldenGatePool { .. } | Method::MultiSiteReference { .. }
    ) && selected(design, &construct.choices)?.iter().any(|(b, c)| {
        preparation_override(construct, &b.id)
            .map(|o| o.preparation_reverse)
            .unwrap_or(c.preparation_reverse)
    }) {
        return Err(bad(
            "reaction-pool orientation is selected by its explicit native path, not an ordered preparation reversal",
        ));
    }
    if let Some(overrides) = &construct.overrides {
        if overrides.preparations.len() > MAX_BINS {
            return Err(Error::TooLarge);
        }
        crate::run_material::native_origins_preflight(
            overrides
                .preparations
                .iter()
                .filter_map(|p| p.preparation.as_ref())
                .filter_map(preparation_origin),
        )?;
        let mut bins = BTreeSet::new();
        for p in &overrides.preparations {
            if !bins.insert(&p.bin_id) || !design.bins.iter().any(|b| b.id == p.bin_id) {
                return Err(bad("invalid or duplicate per-construct preparation bin"));
            }
            if let Some(preparation) = &p.preparation {
                preparation_shape(preparation)?;
                let (_, choice) = selected(design, &construct.choices)?
                    .into_iter()
                    .find(|(b, _)| b.id == p.bin_id)
                    .ok_or_else(|| bad("missing override source"))?;
                preparation_binding(design, &choice.source, preparation)?;
                if material_binding(design, &choice.source).is_some() && p.generation.is_some() {
                    return Err(bad(
                        "native source cannot receive original-generation reset",
                    ));
                }
            }
        }
    }
    if matches!(method, Method::MultiSiteReference { .. }) {
        for (bin, choice) in selected(design, &construct.choices)? {
            let preparation = preparation_override(construct, &bin.id)
                .map(|p| &p.preparation)
                .unwrap_or(&choice.preparation);
            let treatment = preparation_override(construct, &bin.id)
                .map(|p| p.treatment)
                .unwrap_or(choice.treatment);
            if preparation.is_some() || treatment.is_some() {
                return Err(bad(
                    "multi-site reference parents must use no preparation or terminal treatment; clear those declarations explicitly instead of treating a reference as prepared physical material",
                ));
            }
        }
    }
    declaration_bounds(&construct.target)?;
    match &construct.target {
        TargetDeclaration::Existing(target) => binding(target)?,
        TargetDeclaration::OriginalRegions {
            junction_edits,
            substitutions: _,
            edits: _,
        } => {
            let count = design.bins.len().saturating_sub(1) + usize::from(circular);
            let mut junctions = BTreeSet::new();
            if junction_edits.len() > MAX_BINS
                || junction_edits.iter().any(|edit| {
                    edit.junction >= count
                        || !junctions.insert(edit.junction)
                        || edit.inserted.len() > 200
                })
            {
                return Err(bad("invalid literal junction edits"));
            }
        }
    }
    count(construct)
}
pub fn validate_construct(design: &Design, construct: &SelectedConstruct) -> Result<(), Error> {
    validate_design(design)?;
    construct_shape(design, construct)
}
fn resolved_work(
    design: &Design,
    construct: &SelectedConstruct,
    sources: &assembly_intent::Resolved<'_>,
    identities: &[(String, String)],
) -> Result<usize, Error> {
    construct_shape(design, construct)?;
    let (method, _) = effective(design, construct);
    let mut work = 0usize;
    for ((bin, choice), (id, _)) in selected(design, &construct.choices)?
        .into_iter()
        .zip(identities)
    {
        let source = &sources
            .inputs
            .iter()
            .find(|i| i.id == *id)
            .ok_or_else(|| bad("resolved source input is absent"))?
            .molecule;
        let n = source.sequence.len();
        let preparation = preparation_override(construct, &bin.id)
            .map(|p| &p.preparation)
            .unwrap_or(&choice.preparation)
            .as_ref();
        let preparation = if matches!(method, Method::GoldenGatePool { substrates, .. }
            if substrates.iter().any(|row| row.bin_id == bin.id && matches!(row.kind, PoolSubstrateKind::WholeSource)))
        {
            None
        } else {
            preparation
        };
        // Use the very same native preparation cost as recipe preflight. Saved
        // Run/Derived wrappers retain their child recipe costs; no origin state
        // or primer search is replayed merely to reserve batch work.
        let preparation_units = if let Some(preparation) = preparation {
            use assembly_recipe::{Preparation as P, RunPreparation as R};
            let framed = matches!(
                preparation,
                P::ExactPcr { .. }
                    | P::PcrCatalogDigest { .. }
                    | P::OverlapExtensionPcr { .. }
                    | P::OverlapExtensionPcrCatalogDigest { .. }
                    | P::FromRun {
                        preparation: R::ExactPcr { .. }
                            | R::PcrCatalogDigest { .. }
                            | R::OverlapExtensionPcr { .. }
                            | R::OverlapExtensionPcrCatalogDigest { .. },
                        ..
                    }
                    | P::FromDerived {
                        preparation: R::ExactPcr { .. }
                            | R::PcrCatalogDigest { .. }
                            | R::OverlapExtensionPcr { .. }
                            | R::OverlapExtensionPcrCatalogDigest { .. },
                        ..
                    }
            );
            if let Some(fresh) = sources.material(id).filter(|_| framed) {
                let frame = crate::pcr_source::Template::fresh(
                    source,
                    fresh.tracked(),
                    crate::chemistry::Strand::Top,
                )?;
                assembly_recipe::preparation_work_from_frame(frame, preparation)?
            } else {
                assembly_recipe::preparation_work(source, preparation)?
            }
        } else {
            0
        };
        work = work
            .checked_add(n)
            .and_then(|n| n.checked_add(preparation_units))
            .ok_or(Error::TooLarge)?;
        if let Method::GoldenGatePool { options, .. } = method {
            work = work
                .checked_add(
                    n.checked_add(400)
                        .and_then(|n| n.checked_mul(options.enzyme.profile().recognition.len()))
                        .ok_or(Error::TooLarge)?,
                )
                .ok_or(Error::TooLarge)?;
        }
    }
    if let Method::Overlap {
        options, finishing, ..
    } = method
    {
        let physical_search = match finishing {
            OverlapFinishing::ReferenceOnly => 0,
            OverlapFinishing::HomologyAssembly { policy } => policy.max_comparisons,
            OverlapFinishing::HomologyAnnealing { policy } => policy.max_comparisons,
        };
        work = work.checked_add(physical_search).ok_or(Error::TooLarge)?;
        if options.scan_repeated_overlaps {
            work = work
                .checked_add(options.max_comparisons)
                .ok_or(Error::TooLarge)?;
        }
    }
    if let Method::GoldenGatePool { options, .. } = method {
        work = work
            .checked_add(
                options
                    .max_states
                    .checked_mul(options.max_fragments)
                    .ok_or(Error::TooLarge)?,
            )
            .ok_or(Error::TooLarge)?;
    }
    if let Method::MultiSiteReference { request, .. } = method {
        let limit = &request.reaction.limits;
        work = work
            .checked_add(usize::try_from(limit.max_comparisons).map_err(|_| Error::TooLarge)?)
            .and_then(|work| {
                work.checked_add(
                    limit
                        .max_combinations
                        .checked_mul(request.reaction.specificities.len())?,
                )
            })
            .and_then(|work| work.checked_add(limit.max_product_bases))
            .ok_or(Error::TooLarge)?;
    }
    Ok(work.max(1))
}
/// Count only. Source alternatives are never reaction substrates/products.
pub fn matrix_size(design: &Design) -> Result<u128, Error> {
    validate_design(design)?;
    cartesian_size(
        &design
            .bins
            .iter()
            .map(|bin| bin.choices.len())
            .collect::<Vec<_>>(),
    )
}
/// Mixed-radix paging, last bin varying fastest. Never allocate the entire product.
pub fn matrix_page(design: &Design, first: u128, limit: usize) -> Result<MatrixPage, Error> {
    let total = matrix_size(design)?;
    let pages = cartesian_page(
        &design
            .bins
            .iter()
            .map(|bin| bin.choices.len())
            .collect::<Vec<_>>(),
        first,
        limit,
    )?;
    let mut combinations = Vec::with_capacity(pages.len());
    for indices in pages {
        let choices = design
            .bins
            .iter()
            .zip(indices)
            .map(|(bin, index)| ChoiceSelection {
                bin_id: bin.id.clone(),
                choice_id: bin.choices[index].id.clone(),
            })
            .collect::<Vec<_>>();
        let canonical = design
            .bins
            .iter()
            .zip(&choices)
            .map(|(bin, choice)| (&bin.id, &choice.choice_id))
            .collect::<Vec<_>>();
        let id = format!(
            "combination-{}",
            assembly_recipe::hash(&(VERSION, &design.id, canonical))?
        );
        combinations.push(Combination { id, choices });
    }
    Ok(MatrixPage {
        total,
        first,
        combinations,
    })
}
pub fn validate_matrix(matrix: &SelectedMatrix) -> Result<(), Error> {
    validate_design(&matrix.design)?;
    if matrix.constructs.is_empty() || matrix.constructs.len() > MAX_SELECTED {
        return Err(bad("select 1–256 named constructs"));
    }
    let mut ids = BTreeSet::new();
    for construct in &matrix.constructs {
        construct_shape(&matrix.design, construct)?;
        if !identifier(&construct.id) || !ids.insert(&construct.id) || !label(&construct.name, 256)
        {
            return Err(bad("invalid/duplicate selected construct identity"));
        }
    }
    count(matrix)
}
/// Workspace validation establishes actual event ownership, including historical
/// inputs, rather than accepting matching DNA/revision assertions as membership.
pub fn resolve_event<'a>(
    workspace: &'a Workspace,
    event: &EventBinding,
) -> Result<&'a Molecule, Error> {
    binding(event)?;
    if !workspace
        .documents
        .iter()
        .any(|d| d.id == event.document_id)
    {
        return Err(bad("source document is not retained"));
    }
    crate::read_checks::resolve(workspace, event)
}
pub(crate) fn workspace_guard(workspace: &Workspace) -> Result<(), Error> {
    if workspace.documents.len() > crate::workspace::MAX_DOCUMENTS
        || workspace.revisions.len() > crate::workspace::MAX_REVISIONS
        || workspace.history.len() > crate::workspace::MAX_HISTORY_EVENTS
        || workspace.attachments.len() > crate::workspace::MAX_ATTACHMENTS
        || workspace.candidates.len() > crate::workspace::MAX_CANDIDATES
    {
        return Err(Error::TooLarge);
    }
    // Native workspace validation permits unrelated supported rich documents.
    // Scientific ModelBudget applies only to this construct's actual sources.
    // Guard native annotation depth/counts before workspace JSON sizing. This
    // permits native rich unrelated documents; scientific ModelBudget is applied
    // only to the selected inputs and independent target below.
    for molecule in workspace.revisions.values() {
        molecule.validate()?;
    }
    workspace.validate()?;
    Ok(())
}
pub fn resolve(
    design: &Design,
    construct: &SelectedConstruct,
    workspace: &Workspace,
) -> Result<Resolved, Error> {
    validate_design(design)?;
    construct_shape(design, construct)?;
    workspace_guard(workspace)?;
    let resolved = resolve_validated(design, construct, workspace)?;
    preflight_resolved(&resolved)?;
    Ok(resolved)
}
pub(crate) fn resolve_validated(
    design: &Design,
    construct: &SelectedConstruct,
    workspace: &Workspace,
) -> Result<Resolved, Error> {
    construct_shape(design, construct)?;
    let (method, circular) = effective(design, construct);
    let choices = selected(design, &construct.choices)?;
    let mut originals = Vec::with_capacity(choices.len());
    for (_, choice) in &choices {
        let molecule = resolve_event(workspace, &choice.source)?;
        let n = molecule.sequence.len();
        if molecule.sequence.alphabet() != collagene_model::Alphabet::Dna
            || matches!(&choice.selector,
            assembly_intent::Selector::ReferenceRange{start,end,..} if *start>=n || end.checked_sub(*start).is_none_or(|k|k==0 || k>n) || (molecule.topology==Topology::Linear && *end>n))
        {
            return Err(bad(
                "selected contribution is outside its original source or exceeds one turn",
            ));
        }
        originals.push(molecule);
    }
    let preparations = choices
        .iter()
        .map(|(bin, choice)| {
            preparation_override(construct, &bin.id)
                .map(|p| p.preparation.as_ref())
                .unwrap_or(choice.preparation.as_ref())
        })
        .collect::<Vec<_>>();
    for (((_, choice), molecule), preparation) in choices.iter().zip(&originals).zip(&preparations)
    {
        if let Some(preparation) = preparation {
            let event = match preparation {
                assembly_recipe::Preparation::FromRun { origin, .. } => {
                    Some(&origin.selection.event)
                }
                assembly_recipe::Preparation::FromDerived { origin, .. } => Some(&origin.product),
                _ => None,
            };
            if let Some(event) = event {
                if event != &choice.source || event.revision != molecule.revision() {
                    return Err(bad(
                        "saved preparation origin differs from exact selected source event",
                    ));
                }
            }
        }
    }
    crate::run_material::validate_native_origins(
        workspace,
        preparations
            .iter()
            .filter_map(|p| *p)
            .filter_map(preparation_origin)
            .chain(
                choices
                    .iter()
                    .filter_map(|(_, c)| material_binding(design, &c.source))
                    .map(|b| match &b.origin {
                        crate::run_material::MaterialOrigin::Run(o) => {
                            crate::run_material::NativeOriginRef::Run(o)
                        }
                        crate::run_material::MaterialOrigin::Derived(o) => {
                            crate::run_material::NativeOriginRef::Derived(o)
                        }
                    }),
            ),
    )?;
    if let Method::MultiSiteReference { request, .. } = method {
        crate::multisite_recombination::source_preflight(&originals, request)?;
    }
    let existing = match &construct.target {
        TargetDeclaration::Existing(target) => Some(resolve_event(workspace, target)?),
        _ => None,
    };
    if matches!(method, Method::MultiSiteReference { .. }) {
        if let Some(target) = existing {
            if target.sequence.alphabet() != collagene_model::Alphabet::Dna
                || (target.topology == Topology::Circular) != circular
            {
                return Err(bad(
                    "multi-site independent target alphabet/topology differs from the construct declaration",
                ));
            }
        }
    }
    let mut borrowed = originals.clone();
    if let Some(target) = existing {
        borrowed.push(target);
    }
    collagene_model::bounds::model_preflight(&borrowed)?;
    let identities = choices
        .iter()
        .map(|(bin, choice)| source_role_ids(&design.id, &construct.id, &bin.id, &choice.id))
        .collect::<Result<Vec<_>, _>>()?;
    let ids = identities
        .iter()
        .map(|(id, _)| id.clone())
        .collect::<Vec<_>>();
    let fragments = choices
        .iter()
        .zip(&identities)
        .map(|((_, choice), (id, role))| assembly_intent::Contribution {
            role: role.clone(),
            input_id: id.clone(),
            selector: choice.selector.clone(),
        })
        .collect::<Vec<_>>();
    let request = assembly_intent::Request {
        name: construct.name.clone(),
        fragments,
        circular,
        material_bindings: choices
            .iter()
            .zip(&ids)
            .filter_map(|((_, choice), id)| {
                material_binding(design, &choice.source).map(|b| assembly_intent::MaterialBinding {
                    input_id: id.clone(),
                    event: b.event.clone(),
                    origin: b.origin.clone(),
                })
            })
            .collect(),
    };
    let native = ids
        .iter()
        .zip(&originals)
        .map(|(id, m)| (id.as_str(), *m))
        .collect::<Vec<_>>();
    let edited_request = match &construct.target {
        TargetDeclaration::OriginalRegions {
            junction_edits,
            substitutions,
            edits,
        } => Some(assembly_edit_intent::Request {
            base: assembly_substitution_intent::Request {
                base: request.clone(),
                junction_edits: junction_edits.clone(),
                substitutions: substitutions.clone(),
            },
            edits: edits.clone(),
        }),
        _ => None,
    };
    if let Some(edited) = &edited_request {
        assembly_edit_intent::source_preflight(&native, edited, Some(workspace))?;
    } else {
        assembly_intent::preflight(&native, &request, Some(workspace))?;
    }
    let inputs = ids
        .iter()
        .zip(&originals)
        .map(|(id, m)| PlanInput {
            id: id.clone(),
            revision: m.revision(),
            molecule: (**m).clone(),
        })
        .collect::<Vec<_>>();
    let bindings = ids
        .iter()
        .zip(&choices)
        .map(|(id, (_, choice))| InputEventBinding {
            input_id: id.clone(),
            event: choice.source.clone(),
        })
        .collect::<Vec<_>>();
    let source_context = assembly_intent::resolve(&inputs, &request, Some(workspace))?;
    let preparation_work = resolved_work(design, construct, &source_context, &identities)?;
    let (target, target_binding, intent) = match &construct.target {
        TargetDeclaration::Existing(event) => (
            existing.expect("resolved existing target").clone(),
            Some(event.clone()),
            None,
        ),
        TargetDeclaration::OriginalRegions {
            junction_edits,
            substitutions,
            edits,
        } => {
            let edited = edited_request
                .as_ref()
                .expect("typed independent declaration");
            if !edits.is_empty() {
                let intent = assembly_edit_intent::build_resolved(&source_context, edited)?;
                (
                    intent.target.clone(),
                    None,
                    Some(IndependentIntent::Edits(intent)),
                )
            } else if !substitutions.is_empty() {
                let intent =
                    assembly_substitution_intent::report_resolved(&source_context, &edited.base)?;
                (
                    intent.target.clone(),
                    None,
                    Some(IndependentIntent::Substitutions(intent)),
                )
            } else if junction_edits.is_empty() {
                let layout = source_context
                    .rows
                    .iter()
                    .map(|r| (0, r.original_sequence.len()))
                    .collect::<Vec<_>>();
                let (target, fragments, warnings) = assembly_intent::build_layout(
                    &source_context,
                    &request,
                    &layout,
                    &source_context.context,
                )?;
                let intent = assembly_intent::Report {
                    context_revision: source_context.context.clone(),
                    request: request.clone(),
                    target,
                    fragments,
                    warnings,
                };
                intent.validate()?;
                (
                    intent.target.clone(),
                    None,
                    Some(IndependentIntent::Regions(intent)),
                )
            } else {
                let intent = assembly_junction_edit::report_resolved(
                    &source_context,
                    &request,
                    junction_edits,
                )?;
                (
                    intent.target.clone(),
                    None,
                    Some(IndependentIntent::Junctions(intent)),
                )
            }
        }
    };
    drop(source_context);
    let recipe = if let Method::MultiSiteReference { request, selection } = method {
        let sources = identities
            .iter()
            .map(|(id, role)| crate::verification::TopoSource {
                role: role.clone(),
                source: crate::verification::Source::Input(id.clone()),
            })
            .collect::<Vec<_>>();
        // No component selection is included: applying an explicit pin must not
        // invalidate its own source/intent context. Every exact event and literal
        // contribution remains bound even when source DNA or names are equal.
        let source_revision = assembly_recipe::hash(&(
            "combinatorial-multisite-source-intent-1",
            &design.id,
            &construct.id,
            &construct.name,
            &construct.choices,
            &bindings,
            &construct.target,
            circular,
            choices
                .iter()
                .map(|(bin, choice)| {
                    (
                        &bin.id,
                        &choice.id,
                        &choice.selector,
                        choice.preparation_reverse,
                    )
                })
                .collect::<Vec<_>>(),
            target.revision(),
            request,
            &sources,
        ))?;
        Recipe::MultiSiteReference {
            request: request.clone(),
            selection: selection.clone(),
            source_revision,
            sources,
        }
    } else if let Method::GoldenGatePool {
        options,
        substrates,
        selection,
    } = method
    {
        let rows = choices
            .iter()
            .zip(&identities)
            .map(|((bin, choice), (id, role))| {
                let treatment = preparation_override(construct, &bin.id)
                    .map(|o| o.treatment)
                    .unwrap_or(choice.treatment);
                let kind = &substrates
                    .iter()
                    .find(|row| row.bin_id == bin.id)
                    .ok_or_else(|| bad("reaction substrate bin is missing"))?
                    .kind;
                let substrate = match kind {
                    PoolSubstrateKind::WholeSource => {
                        if treatment.is_some() {
                            return Err(bad("whole-source reaction substrate cannot silently apply a terminal treatment; choose an explicit prepared route"));
                        }
                        let origin=material_binding(design,&choice.source).map(|b|b.origin.clone());
                        crate::golden_gate_reaction::Substrate::WholeSource{origin}
                    }
                    PoolSubstrateKind::Prepared => {
                        let preparation = preparation_override(construct, &bin.id)
                            .map(|o| &o.preparation)
                            .unwrap_or(&choice.preparation)
                            .clone()
                            .ok_or_else(|| bad("prepared reaction substrate needs an explicit preparation"))?;
                        crate::golden_gate_reaction::Substrate::Prepared { preparation, treatment }
                    }
                };
                Ok(crate::golden_gate_reaction::Row {
                    role: role.clone(),
                    input_id: id.clone(),
                    substrate,
                    generation:preparation_override(construct,&bin.id).map(|o|o.generation.clone()).unwrap_or_else(||choice.generation.clone()),
                })
            })
            .collect::<Result<Vec<_>, Error>>()?;
        let request = crate::golden_gate_reaction::Request {
            name: construct.name.clone(),
            rows,
            options: options.clone(),
        };
        let source_revision = assembly_recipe::hash(&(
            "combinatorial-pool-source-intent-1",
            &design.id,
            &construct.id,
            &construct.choices,
            &bindings,
            &construct.target,
            choices
                .iter()
                .map(|(bin, choice)| {
                    (
                        &bin.id,
                        &choice.id,
                        &choice.selector,
                        choice.preparation_reverse,
                    )
                })
                .collect::<Vec<_>>(),
            target.revision(),
            &request,
        ))?;
        Recipe::GoldenGatePool {
            request,
            selection: selection.clone(),
            source_revision,
        }
    } else {
        let mut rows = Vec::with_capacity(choices.len());
        for ((bin, choice), (id, role)) in choices.iter().zip(&identities) {
            let preparation = preparation_override(construct, &bin.id)
                .map(|o| &o.preparation)
                .unwrap_or(&choice.preparation)
                .clone()
                .ok_or_else(|| bad("chosen alternative has no explicit preparation"))?;
            rows.push(assembly_recipe::FragmentSpec {
                generation: preparation_override(construct, &bin.id)
                    .map(|o| o.generation.clone())
                    .unwrap_or_else(|| choice.generation.clone()),
                role: role.clone(),
                input_id: id.clone(),
                preparation,
                reverse: preparation_override(construct, &bin.id)
                    .map(|o| o.preparation_reverse)
                    .unwrap_or(choice.preparation_reverse),
            });
        }
        match method {
            Method::Tailing {
                profile,
                actions,
                allow_nicked,
            } => {
                let fragments = rows
                    .into_iter()
                    .zip(&choices)
                    .map(|(source, (bin, choice))| {
                        let action = actions
                            .iter()
                            .find(|row| row.bin_id == bin.id)
                            .ok_or_else(|| bad("chosen bin has no explicit tail action"))?
                            .action;
                        Ok(tailing_recipe::Fragment {
                            source,
                            action,
                            treatment: preparation_override(construct, &bin.id)
                                .map(|o| o.treatment)
                                .unwrap_or(choice.treatment),
                        })
                    })
                    .collect::<Result<Vec<_>, Error>>()?;
                Recipe::Tailing {
                    recipe: tailing_recipe::Recipe {
                        name: construct.name.clone(),
                        profile: *profile,
                        fragments,
                        circular,
                        allow_nicked: *allow_nicked,
                    },
                }
            }
            Method::Overlap {
                overlaps,
                options,
                finishing,
            } => {
                if choices.iter().any(|(bin, choice)| {
                    preparation_override(construct, &bin.id)
                        .map(|o| o.treatment)
                        .unwrap_or(choice.treatment)
                        .is_some()
                }) {
                    return Err(bad(
                        "terminal treatments require an explicit ligation route",
                    ));
                }
                let recipe = assembly_recipe::Recipe {
                    name: construct.name.clone(),
                    fragments: rows,
                    overlaps: overlaps.clone(),
                    circular,
                    options: options.clone(),
                };
                Recipe::Overlap {
                    recipe,
                    finishing: finishing.clone(),
                }
            }
            Method::Ligation { allow_nicked } | Method::OrderedTypeIis { allow_nicked, .. } => {
                if let Method::OrderedTypeIis { enzyme, .. } = method {
                    let profile = enzyme.profile();
                    for row in &rows {
                        let selection=assembly_recipe::preparation_catalog_selection(&row.preparation).ok_or_else(||bad("Golden Gate needs pinned existing-site or explicit PCR/catalog-digest preparation"))?;
                        if selection.enzymes.len() != 1
                            || !selection.enzymes.iter().all(|entry| {
                                entry.enzyme.recognition == profile.recognition
                                    && entry.enzyme.top_cut == profile.top_cut
                                    && entry.enzyme.bottom_cut == profile.bottom_cut
                            })
                        {
                            return Err(bad(
                                "chosen Type IIS pin differs from the design enzyme geometry",
                            ));
                        }
                    }
                }
                let recipe = ligation_recipe::Recipe {
                    name: construct.name.clone(),
                    fragments: rows
                        .into_iter()
                        .zip(&choices)
                        .map(|(source, (bin, choice))| ligation_recipe::Fragment {
                            source,
                            treatment: preparation_override(construct, &bin.id)
                                .map(|o| o.treatment)
                                .unwrap_or(choice.treatment),
                        })
                        .collect(),
                    circular,
                    allow_nicked: *allow_nicked,
                };
                Recipe::Ligation {
                    recipe,
                    cleavage_policy: if matches!(method, Method::OrderedTypeIis { .. }) {
                        assembly_recipe::CleavagePolicy::PositionalChecked
                    } else {
                        assembly_recipe::CleavagePolicy::AssumedComplete
                    },
                }
            }
            Method::GoldenGatePool { .. } | Method::MultiSiteReference { .. } => {
                unreachable!("whole-parent reaction compiled above")
            }
        }
    };
    let combination_id = combination_key(design, &construct.choices)?;
    let declaration_revision =
        assembly_recipe::hash(&(VERSION, design, construct, &bindings, target.revision()))?;
    Ok(Resolved {
        construct_id: construct.id.clone(),
        construct_name: construct.name.clone(),
        combination_id,
        declaration_revision,
        inputs,
        bindings,
        target,
        target_binding,
        intent,
        recipe,
        preparation_work,
    })
}

fn preflight_resolved(resolved: &Resolved) -> Result<(), Error> {
    let borrowed = resolved
        .inputs
        .iter()
        .map(|i| (i.id.as_str(), &i.molecule))
        .collect::<Vec<_>>();
    match &resolved.recipe {
        Recipe::Tailing { recipe } => {
            tailing_recipe::source_preflight(&borrowed, &resolved.target, recipe)
        }
        Recipe::Overlap { recipe, finishing } => {
            finishing.preflight(&resolved.inputs, &resolved.target, recipe)
        }
        Recipe::Ligation { recipe, .. } => {
            ligation_recipe::source_preflight(&borrowed, &resolved.target, recipe)
        }
        Recipe::GoldenGatePool { request, .. } => {
            crate::golden_gate_reaction::preflight(&borrowed, &resolved.target, request)
        }
        Recipe::MultiSiteReference { .. } => Ok(()),
    }
}
fn multisite_inputs(
    resolved: &Resolved,
) -> Result<Vec<crate::recombination::reference::Input>, Error> {
    let Recipe::MultiSiteReference {
        request, sources, ..
    } = &resolved.recipe
    else {
        return Err(bad(
            "resolved construct is not a multi-site reference design",
        ));
    };
    if sources.len() != resolved.inputs.len() {
        return Err(bad("multi-site source role/input counts differ"));
    }
    let mut borrowed = Vec::with_capacity(sources.len());
    let mut roles = BTreeSet::new();
    for source in sources {
        if !roles.insert(&source.role) {
            return Err(bad("multi-site parent roles must be distinct"));
        }
        let crate::verification::Source::Input(id) = &source.source else {
            return Err(bad(
                "multi-site library parents must be whole original inputs",
            ));
        };
        let input = resolved
            .inputs
            .iter()
            .find(|input| input.id == *id)
            .ok_or_else(|| bad("multi-site parent input is absent"))?;
        borrowed.push((&source.role, &input.molecule));
    }
    crate::multisite_recombination::source_preflight(
        &borrowed
            .iter()
            .map(|(_, molecule)| *molecule)
            .collect::<Vec<_>>(),
        request,
    )?;
    Ok(borrowed
        .into_iter()
        .map(|(role, molecule)| crate::recombination::reference::Input {
            role: role.clone(),
            molecule: molecule.clone(),
        })
        .collect())
}
/// Resolve whole original parents and freeze independent target intent. A missing
/// component selection remains an inspection draft, not a default candidate.
pub fn resolve_multisite(
    design: &Design,
    construct: &SelectedConstruct,
    workspace: &Workspace,
) -> Result<Resolved, Error> {
    let resolved = resolve(design, construct, workspace)?;
    if !matches!(&resolved.recipe, Recipe::MultiSiteReference { .. }) {
        return Err(bad(
            "selected construct is not a multi-site reference design",
        ));
    }
    Ok(resolved)
}
/// Borrowed inspection context before source cloning, target construction or native
/// enumeration. The workspace fingerprint retains equal-DNA distinct events.
pub fn inspect_multisite_revision(
    design: &Design,
    construct: &SelectedConstruct,
    workspace: &Workspace,
) -> Result<String, Error> {
    validate_construct(design, construct)?;
    workspace_guard(workspace)?;
    let (method, circular) = effective(design, construct);
    let Method::MultiSiteReference { request, .. } = method else {
        return Err(bad(
            "selected construct is not a multi-site reference design",
        ));
    };
    let selected = selected(design, &construct.choices)?;
    let mut parents = Vec::with_capacity(selected.len());
    for (_, choice) in selected {
        let molecule = resolve_event(workspace, &choice.source)?;
        let n = molecule.sequence.len();
        if matches!(&choice.selector,assembly_intent::Selector::ReferenceRange{start,end,..} if *start>=n || end.checked_sub(*start).is_none_or(|k|k==0 || k>n) || (molecule.topology==Topology::Linear && *end>n))
        {
            return Err(bad(
                "independent contribution is outside its original parent",
            ));
        }
        parents.push(molecule);
    }
    crate::multisite_recombination::source_preflight(&parents, request)?;
    if let TargetDeclaration::Existing(binding) = &construct.target {
        let target = resolve_event(workspace, binding)?;
        collagene_model::bounds::model_preflight(&[target])?;
        if target.sequence.alphabet() != collagene_model::Alphabet::Dna
            || (target.topology == Topology::Circular) != circular
        {
            return Err(bad(
                "independent target alphabet/topology differs from this construct",
            ));
        }
    }
    assembly_recipe::hash(&(
        "combinatorial-multisite-inspection-1",
        crate::run_lineage::workspace_fingerprint(workspace)?,
        design,
        construct,
    ))
}
/// Native reference discovery runs only after the independent intent was frozen.
/// Neither that target nor a stored pin selects the enumeration's pairings.
pub fn inspect_multisite(
    design: &Design,
    construct: &SelectedConstruct,
    workspace: &Workspace,
) -> Result<MultiSiteInspection, Error> {
    let context_revision = inspect_multisite_revision(design, construct, workspace)?;
    let resolved = resolve_multisite(design, construct, workspace)?;
    let Recipe::MultiSiteReference {
        request,
        source_revision,
        ..
    } = &resolved.recipe
    else {
        unreachable!()
    };
    let inputs = multisite_inputs(&resolved)?;
    let batch = crate::multisite_recombination::simulate(&inputs, request)?;
    let result = MultiSiteInspection {
        context_revision: context_revision.clone(),
        source_revision: source_revision.clone(),
        batch,
    };
    result.validate_binding(&context_revision)?;
    Ok(result)
}
impl MultiSiteInspection {
    /// Readonly current response/source-selection context, not a scientific replay.
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        if !digest(expected) || self.context_revision != expected || !digest(&self.source_revision)
        {
            return Err(bad("multi-site inspection context differs"));
        }
        self.batch.validate_binding(&self.batch.batch_revision)?;
        count(self)
    }
}

/// Pure plan compilation. The returned plan has not been scientifically verified.
pub fn compile(
    design: &Design,
    construct: &SelectedConstruct,
    workspace: &Workspace,
) -> Result<Compiled, Error> {
    let resolved = resolve(design, construct, workspace)?;
    let plan = compile_resolved(&resolved)?;
    Ok(Compiled { resolved, plan })
}
/// Compile an already source-resolved declaration without running scientific
/// enumeration. Callers must not treat the resulting plan as verification evidence.
pub fn compile_resolved(resolved: &Resolved) -> Result<crate::verification::DesignPlan, Error> {
    match &resolved.recipe {
        Recipe::Tailing { recipe } => {
            Ok(tailing_recipe::compile(&resolved.inputs, &resolved.target, recipe)?.plan)
        }
        Recipe::Overlap { recipe, finishing } => Ok(finishing
            .compile(&resolved.inputs, &resolved.target, recipe)?
            .plan),
        Recipe::Ligation {
            recipe,
            cleavage_policy,
        } => Ok(ligation_recipe::compile(
            &resolved.inputs,
            &resolved.target,
            recipe,
            *cleavage_policy,
        )?
        .plan),
        Recipe::MultiSiteReference {
            request,
            selection,
            source_revision,
            sources,
        } => {
            let pin=selection.as_ref().ok_or_else(||bad("inspect the multi-site network and explicitly select a candidate and reciprocal component before compilation"))?;
            if &pin.source_revision != source_revision {
                return Err(bad(
                    "multi-site selection belongs to different exact source events, contribution geometry, or independent target intent",
                ));
            }
            let native = multisite_inputs(resolved)?;
            crate::multisite_recombination::product_revision(&native, request, &pin.selection)?;
            Ok(crate::verification::DesignPlan {
                name: resolved.construct_name.clone(),
                target_revision: resolved.target.revision(),
                inputs: resolved.inputs.clone(),
                steps: vec![
                    crate::verification::Operation::MultiSiteRecombineReference {
                        sources: sources.clone(),
                        request: request.clone(),
                        selection: pin.selection.clone(),
                    },
                ],
            })
        }
        Recipe::GoldenGatePool {
            request,
            selection,
            source_revision,
        } => {
            let pin = selection.as_ref().ok_or_else(|| {
                bad("inspect the reaction and explicitly select a candidate before compilation")
            })?;
            if &pin.source_revision != source_revision {
                return Err(bad(
                    "reaction selection belongs to different exact source events, contribution geometry, or independent target intent",
                ));
            }
            crate::golden_gate_reaction::compile_selected(
                &resolved.inputs,
                &resolved.target,
                request,
                &pin.selection,
            )
        }
    }
}
/// Resolve exact original sources and the independent target for pool discovery.
/// A missing candidate selection is intentional and does not select candidate zero.
pub fn resolve_pool(
    design: &Design,
    construct: &SelectedConstruct,
    workspace: &Workspace,
) -> Result<Resolved, Error> {
    let resolved = resolve(design, construct, workspace)?;
    if !matches!(&resolved.recipe, Recipe::GoldenGatePool { .. }) {
        return Err(bad(
            "selected construct is not an explicit Golden Gate reaction pool",
        ));
    }
    Ok(resolved)
}
/// Fresh bounded discovery only. A report or target match grants no output authority.
pub fn inspect_pool(
    design: &Design,
    construct: &SelectedConstruct,
    workspace: &Workspace,
) -> Result<crate::golden_gate_reaction::Report, Error> {
    let resolved = resolve_pool(design, construct, workspace)?;
    let Recipe::GoldenGatePool { request, .. } = &resolved.recipe else {
        unreachable!()
    };
    crate::golden_gate_reaction::inspect(&resolved.inputs, &resolved.target, request)
}
/// One explicitly selected construct, reconstructed by the existing native
/// preparation compiler/verifier. This is not batch completion or mixture search.
pub fn evaluate(
    design: &Design,
    construct: &SelectedConstruct,
    workspace: &Workspace,
) -> Result<Evaluation, Error> {
    let resolved = resolve(design, construct, workspace)?;
    evaluate_resolved(resolved)
}
pub(crate) fn evaluate_resolved(resolved: Resolved) -> Result<Evaluation, Error> {
    let preview = match &resolved.recipe {
        Recipe::Tailing { recipe } => Preview::Tailing(tailing_recipe::preview(
            &resolved.inputs,
            &resolved.target,
            recipe,
        )?),
        Recipe::Overlap { recipe, finishing } => {
            Preview::Overlap(finishing.preview(&resolved.inputs, &resolved.target, recipe)?)
        }
        Recipe::Ligation {
            recipe,
            cleavage_policy,
        } => Preview::Ligation(ligation_recipe::preview(
            &resolved.inputs,
            &resolved.target,
            recipe,
            *cleavage_policy,
        )?),
        Recipe::MultiSiteReference { .. } => {
            let plan = compile_resolved(&resolved)?;
            let verification = crate::verification::verify(&plan, &resolved.target);
            let preparation_complete = verification.products.len() == plan.steps.len();
            Preview::MultiSiteReference(PoolPreview {
                plan,
                verification,
                preparation_complete,
            })
        }
        Recipe::GoldenGatePool { .. } => {
            let plan = compile_resolved(&resolved)?;
            let verification = crate::verification::verify(&plan, &resolved.target);
            let preparation_complete =
                verification.products.len() >= plan.steps.len().saturating_sub(1);
            Preview::GoldenGatePool(PoolPreview {
                plan,
                verification,
                preparation_complete,
            })
        }
    };
    Ok(Evaluation { resolved, preview })
}
impl Evaluation {
    pub fn verification(&self) -> &crate::verification::VerificationReport {
        match &self.preview {
            Preview::Tailing(preview) => &preview.verification,
            Preview::Overlap(preview) => &preview.verification,
            Preview::Ligation(preview) => &preview.verification,
            Preview::GoldenGatePool(preview) | Preview::MultiSiteReference(preview) => {
                &preview.verification
            }
        }
    }
    pub fn plan(&self) -> &crate::verification::DesignPlan {
        match &self.preview {
            Preview::Tailing(preview) => &preview.compiled.plan,
            Preview::Overlap(preview) => &preview.compiled.plan,
            Preview::Ligation(preview) => &preview.compiled.plan,
            Preview::GoldenGatePool(preview) | Preview::MultiSiteReference(preview) => {
                &preview.plan
            }
        }
    }
    pub fn preparation_complete(&self) -> bool {
        match &self.preview {
            Preview::Tailing(preview) => preview.preparation_complete,
            Preview::Overlap(preview) => preview.preparation_complete,
            Preview::Ligation(preview) => preview.preparation_complete,
            Preview::GoldenGatePool(preview) | Preview::MultiSiteReference(preview) => {
                preview.preparation_complete
            }
        }
    }
}
