//! Explicit source alternatives over existing preparation compilers.
//! This is a design-choice layer, independent of reaction-product enumeration.
use crate::{
    Error, Molecule, Topology, assembly_intent, assembly_junction_edit, assembly_recipe,
    chemistry::EndTreatment,
    golden_gate::GoldenGateEnzyme,
    ligation_recipe,
    overlap_profiles::ProfileOptions,
    run_lineage::{EventBinding, InputEventBinding},
    verification::PlanInput,
    workspace::Workspace,
};
use serde::{Deserialize, Serialize};
use std::{collections::BTreeSet, io};

const VERSION: &str = "combinatorial-source-design-1";
pub const MAX_BINS: usize = 12;
pub const MAX_CHOICES: usize = 256;
pub const MAX_SELECTED: usize = 256;
pub const MAX_PAGE: usize = 64;
const MAX_BYTES: usize = 32 * 1024 * 1024;
macro_rules! dto {($($t:tt)*)=>{#[derive(Clone,Debug,Serialize,Deserialize)]#[cfg_attr(feature="json-schema",derive(schemars::JsonSchema))]#[serde(deny_unknown_fields)]$($t)*};}
dto! {pub struct Choice {
    pub id:String,pub label:String,pub source:EventBinding,
    /// Literal original-source contribution, not a prepared/PCR output.
    /// Start normalized, end unrolled at most one circular turn.
    pub start:usize,pub end:usize,pub reverse:bool,
    /// Unfinished alternatives can be retained; an explicitly chosen one needs preparation.
    pub preparation:Option<assembly_recipe::Preparation>,
    pub treatment:Option<EndTreatment>,
}}
dto! {pub struct Bin {pub id:String,pub label:String,pub choices:Vec<Choice>}}
dto! {pub enum Method {
    Overlap{overlaps:Vec<usize>,options:ProfileOptions},
    Ligation{allow_nicked:bool},
    /// Same existing pinned-site geometry requirement; no kinetic-mixture inference.
    GoldenGate{enzyme:GoldenGateEnzyme,allow_nicked:bool},
}}
dto! {pub struct Design {pub id:String,pub name:String,pub circular:bool,pub bins:Vec<Bin>,pub method:Method}}
dto! {pub struct ChoiceSelection {pub bin_id:String,pub choice_id:String}}
dto! {pub struct Combination {pub id:String,pub choices:Vec<ChoiceSelection>}}
dto! {pub struct MatrixPage {pub total:u128,pub first:u128,pub combinations:Vec<Combination>}}
dto! {pub enum TargetDeclaration {
    /// An explicitly selected independent retained reference event.
    Existing(EventBinding),
    /// Literal original contribution coordinates + optional declared junction edits.
    /// Uses existing independent intent builder, never the reconstructed product.
    OriginalRegions{junction_edits:Vec<assembly_junction_edit::Edit>},
}}
dto! {pub struct SelectedConstruct {
    pub id:String,pub name:String,pub choices:Vec<ChoiceSelection>,pub target:TargetDeclaration,
}}
dto! {pub struct SelectedMatrix {pub design:Design,pub constructs:Vec<SelectedConstruct>}}
#[derive(Clone, Debug)]
pub enum IndependentIntent {
    Regions(assembly_intent::Report),
    Junctions(assembly_junction_edit::Intent),
}
#[derive(Clone, Debug)]
pub enum Recipe {
    Overlap(assembly_recipe::Recipe),
    Ligation(ligation_recipe::Recipe),
}
pub struct Resolved {
    pub construct_id: String,
    pub combination_id: String,
    pub declaration_revision: String,
    pub inputs: Vec<PlanInput>,
    pub bindings: Vec<InputEventBinding>,
    pub target: Molecule,
    pub target_binding: Option<EventBinding>,
    pub intent: Option<IndependentIntent>,
    pub recipe: Recipe,
}
pub struct Compiled {
    pub resolved: Resolved,
    pub plan: crate::verification::DesignPlan,
}
pub enum Preview {
    Overlap(assembly_recipe::Preview),
    Ligation(ligation_recipe::Preview),
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
pub fn validate_design(design: &Design) -> Result<(), Error> {
    if !identifier(&design.id)
        || !label(&design.name, 256)
        || design.bins.is_empty()
        || design.bins.len() > MAX_BINS
    {
        return Err(bad("invalid design identity or bin count"));
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
                || choice.start >= choice.end
                || choice.end - choice.start > crate::model::MAX_SEQUENCE_LEN
            {
                return Err(bad(
                    "invalid/duplicate choice identity or original contribution",
                ));
            }
            binding(&choice.source)?;
            if let Some(preparation) = &choice.preparation {
                match preparation {
                    assembly_recipe::Preparation::ExactPcr { forward, reverse }
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
                                || primer.annealing.alphabet() != crate::Alphabet::Dna
                            {
                                return Err(bad("invalid explicit primer"));
                            }
                            let parsed = crate::analysis::Primer::new(
                                &primer.name,
                                primer.annealing.as_str(),
                                &primer.tail,
                            )?;
                            if parsed.tail != primer.tail || parsed.annealing != primer.annealing {
                                return Err(bad(
                                    "typed preparation must retain canonical primer DNA",
                                ));
                            }
                        }
                    }
                    assembly_recipe::Preparation::Digest { enzymes, .. } => {
                        if enzymes.len() > 64 || enzymes.iter().any(|e| e.len() > 256) {
                            return Err(Error::TooLarge);
                        }
                    }
                    _ => {}
                }
                match preparation {
                    assembly_recipe::Preparation::CatalogDigest { selection, .. }
                    | assembly_recipe::Preparation::PcrCatalogDigest { selection, .. } => {
                        selection.validate()?
                    }
                    _ => {}
                }
            }
        }
    }
    if let Method::Overlap { overlaps, options } = &design.method {
        if overlaps.len() != design.bins.len().saturating_sub(1) + usize::from(design.circular)
            || overlaps
                .iter()
                .any(|length| *length == 0 || *length > 10_000)
        {
            return Err(bad("declare one overlap length per junction"));
        }
        options.validate()?;
    }
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
fn combination_key(design: &Design, choices: &[ChoiceSelection]) -> Result<String, Error> {
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
    match &construct.target {
        TargetDeclaration::Existing(target) => binding(target)?,
        TargetDeclaration::OriginalRegions { junction_edits } => {
            let count = design.bins.len().saturating_sub(1) + usize::from(design.circular);
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
    let mut combinations = BTreeSet::new();
    for construct in &matrix.constructs {
        construct_shape(&matrix.design, construct)?;
        if !identifier(&construct.id) || !ids.insert(&construct.id) || !label(&construct.name, 256)
        {
            return Err(bad("invalid/duplicate selected construct identity"));
        }
        let combination = combination_key(&matrix.design, &construct.choices)?;
        if !combinations.insert(combination) {
            return Err(bad("the same combination was selected twice"));
        }
    }
    count(matrix)
}
/// Workspace validation establishes actual event ownership, including historical
/// inputs, rather than accepting matching DNA/revision assertions as membership.
fn resolve_event<'a>(
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
pub fn resolve(
    design: &Design,
    construct: &SelectedConstruct,
    workspace: &Workspace,
) -> Result<Resolved, Error> {
    validate_design(design)?;
    construct_shape(design, construct)?;
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
    let choices = selected(design, &construct.choices)?;
    let mut originals = Vec::with_capacity(choices.len());
    for (_, choice) in &choices {
        let molecule = resolve_event(workspace, &choice.source)?;
        let n = molecule.sequence.len();
        if molecule.sequence.alphabet() != crate::Alphabet::Dna
            || choice.start >= n
            || choice.end - choice.start > n
            || (molecule.topology == Topology::Linear && choice.end > n)
        {
            return Err(bad(
                "selected contribution is outside its original source or exceeds one turn",
            ));
        }
        originals.push(molecule);
    }
    let existing = match &construct.target {
        TargetDeclaration::Existing(target) => Some(resolve_event(workspace, target)?),
        _ => None,
    };
    let mut borrowed = originals.clone();
    if let Some(target) = existing {
        borrowed.push(target);
    }
    crate::model::bounds::model_preflight(&borrowed)?;
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
        .map(
            |((_, choice), (id, role))| crate::overlap_primer_design::Fragment {
                role: role.clone(),
                input_id: id.clone(),
                start: choice.start,
                end: choice.end,
                reverse: choice.reverse,
            },
        )
        .collect::<Vec<_>>();
    let request = assembly_intent::Request {
        name: construct.name.clone(),
        fragments,
        circular: design.circular,
    };
    let native = ids
        .iter()
        .zip(&originals)
        .map(|(id, m)| (id.as_str(), *m))
        .collect::<Vec<_>>();
    assembly_intent::preflight(&native, &request)?;
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
    let (target, target_binding, intent) = match &construct.target {
        TargetDeclaration::Existing(event) => (
            existing.expect("resolved existing target").clone(),
            Some(event.clone()),
            None,
        ),
        TargetDeclaration::OriginalRegions { junction_edits } => {
            if junction_edits.is_empty() {
                let intent = assembly_intent::build(&inputs, &request)?;
                (
                    intent.target.clone(),
                    None,
                    Some(IndependentIntent::Regions(intent)),
                )
            } else {
                let intent =
                    assembly_junction_edit::build_intent(&inputs, &request, junction_edits)?;
                (
                    intent.target.clone(),
                    None,
                    Some(IndependentIntent::Junctions(intent)),
                )
            }
        }
    };
    let mut rows = Vec::with_capacity(choices.len());
    for ((_, choice), (id, role)) in choices.iter().zip(&identities) {
        let preparation = choice
            .preparation
            .clone()
            .ok_or_else(|| bad("chosen alternative has no explicit preparation"))?;
        rows.push(assembly_recipe::FragmentSpec {
            role: role.clone(),
            input_id: id.clone(),
            preparation,
            reverse: choice.reverse,
        });
    }
    let recipe = match &design.method {
        Method::Overlap { overlaps, options } => {
            if choices.iter().any(|(_, choice)| choice.treatment.is_some()) {
                return Err(bad(
                    "terminal treatments require an explicit ligation route",
                ));
            }
            let recipe = assembly_recipe::Recipe {
                name: construct.name.clone(),
                fragments: rows,
                overlaps: overlaps.clone(),
                circular: design.circular,
                options: options.clone(),
            };
            assembly_recipe::preflight(&inputs, &target, &recipe)?;
            Recipe::Overlap(recipe)
        }
        Method::Ligation { allow_nicked } | Method::GoldenGate { allow_nicked, .. } => {
            if let Method::GoldenGate { enzyme, .. } = &design.method {
                let profile = enzyme.profile();
                for row in &rows {
                    let selection = match &row.preparation {
                        assembly_recipe::Preparation::CatalogDigest { selection, .. }
                        | assembly_recipe::Preparation::PcrCatalogDigest { selection, .. } => {
                            selection
                        }
                        _ => {
                            return Err(bad(
                                "Golden Gate needs pinned existing-site or PCR/catalog-digest preparation",
                            ));
                        }
                    };
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
                    .map(|(source, (_, choice))| ligation_recipe::Fragment {
                        source,
                        treatment: choice.treatment,
                    })
                    .collect(),
                circular: design.circular,
                allow_nicked: *allow_nicked,
            };
            ligation_recipe::source_preflight(&native, &target, &recipe)?;
            Recipe::Ligation(recipe)
        }
    };
    let combination_id = combination_key(design, &construct.choices)?;
    let declaration_revision =
        assembly_recipe::hash(&(VERSION, design, construct, &bindings, target.revision()))?;
    Ok(Resolved {
        construct_id: construct.id.clone(),
        combination_id,
        declaration_revision,
        inputs,
        bindings,
        target,
        target_binding,
        intent,
        recipe,
    })
}

/// Pure plan compilation. The returned plan has not been scientifically verified.
pub fn compile(
    design: &Design,
    construct: &SelectedConstruct,
    workspace: &Workspace,
) -> Result<Compiled, Error> {
    let resolved = resolve(design, construct, workspace)?;
    let plan = match &resolved.recipe {
        Recipe::Overlap(recipe) => {
            assembly_recipe::compile(&resolved.inputs, &resolved.target, recipe)?.plan
        }
        Recipe::Ligation(recipe) => {
            ligation_recipe::compile(&resolved.inputs, &resolved.target, recipe)?.plan
        }
    };
    Ok(Compiled { resolved, plan })
}
/// One explicitly selected construct, reconstructed by the existing native
/// preparation compiler/verifier. This is not batch completion or mixture search.
pub fn evaluate(
    design: &Design,
    construct: &SelectedConstruct,
    workspace: &Workspace,
) -> Result<Evaluation, Error> {
    let resolved = resolve(design, construct, workspace)?;
    let preview = match &resolved.recipe {
        Recipe::Overlap(recipe) => Preview::Overlap(assembly_recipe::preview(
            &resolved.inputs,
            &resolved.target,
            recipe,
        )?),
        Recipe::Ligation(recipe) => Preview::Ligation(ligation_recipe::preview(
            &resolved.inputs,
            &resolved.target,
            recipe,
        )?),
    };
    Ok(Evaluation { resolved, preview })
}
impl Evaluation {
    pub fn verification(&self) -> &crate::verification::VerificationReport {
        match &self.preview {
            Preview::Overlap(preview) => &preview.verification,
            Preview::Ligation(preview) => &preview.verification,
        }
    }
    pub fn plan(&self) -> &crate::verification::DesignPlan {
        match &self.preview {
            Preview::Overlap(preview) => &preview.compiled.plan,
            Preview::Ligation(preview) => &preview.compiled.plan,
        }
    }
}
