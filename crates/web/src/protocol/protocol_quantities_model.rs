//! Documentary numeric declarations with exact native event identity.
use collagene_core::{
    documents::workspace::Workspace,
    material::lineage::{EventBinding, InputEventBinding},
    output::quantities::{self as native, units},
    planning::plan::Source,
};
use serde::Serialize;
#[derive(Clone, Serialize)]
pub struct Binding {
    pub plan_source: native::PlanSource,
    pub plan_revision: String,
    pub target: EventBinding,
    pub inputs: Vec<InputEventBinding>,
    pub step_count: usize,
    pub component_sources: Vec<Vec<Source>>,
    #[serde(skip)]
    pub operation_kinds: Vec<std::mem::Discriminant<collagene_core::planning::plan::Operation>>,
}
fn key(value: &impl Serialize) -> Result<Vec<u8>, String> {
    struct Bytes(Vec<u8>);
    impl std::io::Write for Bytes {
        fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
            if bytes.len() > native::MAX_BYTES.saturating_sub(self.0.len()) {
                return Err(std::io::Error::other(
                    "Quantity context exceeds its native bound",
                ));
            }
            self.0.extend_from_slice(bytes);
            Ok(bytes.len())
        }
        fn flush(&mut self) -> std::io::Result<()> {
            Ok(())
        }
    }
    let mut out = Bytes(Vec::new());
    serde_json::to_writer(&mut out, value).map_err(|e| e.to_string())?;
    Ok(out.0)
}
#[derive(Clone)]
pub struct Scope {
    pub epoch: u64,
    pub fingerprint: String,
    binding: Vec<u8>,
}
impl Scope {
    pub fn capture(w: &Workspace, binding: &Binding, epoch: u64) -> Result<Self, String> {
        w.validate().map_err(|e| e.to_string())?;
        Ok(Self {
            epoch,
            fingerprint: collagene_core::material::lineage::workspace_fingerprint(w)
                .map_err(|e| e.to_string())?,
            binding: key(binding)?,
        })
    }
    pub fn permit(&self, w: &Workspace, binding: &Binding, epoch: u64) -> Result<(), String> {
        if epoch != self.epoch
            || key(binding)? != self.binding
            || collagene_core::material::lineage::workspace_fingerprint(w)
                .map_err(|e| e.to_string())?
                != self.fingerprint
        {
            return Err("The exact plan, target/source events or workspace changed; review quantities again.".into());
        }
        Ok(())
    }
}
#[derive(Clone, PartialEq, Eq)]
pub struct ComponentForm {
    pub id: String,
    pub label: String,
    pub source: String,
    pub mode: String,
    pub volume: String,
    pub dilution_volume: String,
    pub dilution_unit: String,
    pub unit: String,
    pub stock: String,
    pub target: String,
    pub concentration_unit: String,
    pub rounding: String,
    pub increment: String,
}
impl ComponentForm {
    pub fn blank(id: String) -> Self {
        Self {
            id,
            label: String::new(),
            source: String::new(),
            mode: "direct".into(),
            volume: String::new(),
            dilution_volume: String::new(),
            dilution_unit: "uL".into(),
            unit: "uL".into(),
            stock: String::new(),
            target: String::new(),
            concentration_unit: "nM".into(),
            rounding: String::new(),
            increment: String::new(),
        }
    }
}
#[derive(Clone, PartialEq, Eq)]
pub struct ReactionForm {
    pub id: String,
    pub step: String,
    pub total: String,
    pub unit: String,
    pub fill: bool,
    pub fill_id: String,
    pub fill_label: String,
    pub components: Vec<ComponentForm>,
}
#[derive(Clone, PartialEq, Eq)]
pub struct Form {
    pub name: String,
    pub reactions: Vec<ReactionForm>,
    pub next_id: u64,
}
impl Default for Form {
    fn default() -> Self {
        Self {
            name: "Reaction quantities".into(),
            reactions: vec![],
            next_id: 1,
        }
    }
}
fn optional_number(text: &str) -> Result<Option<units::Rational>, String> {
    if text.len() > 64 {
        return Err("Decimal text exceeds64bytes.".into());
    }
    if text.trim().is_empty() {
        Ok(None)
    } else {
        let text = text.trim();
        let value = if let Some((n, d)) = text.split_once('/') {
            let parse = |v: &str| {
                v.parse::<u64>().map_err(|_| {
                    "Use an unsigned integer numerator and nonzero denominator.".to_string()
                })
            };
            units::Rational::new(parse(n)?, parse(d)?).map_err(|e| e.to_string())?
        } else {
            units::Rational::parse_decimal(text).map_err(|e| e.to_string())?
        };
        Ok(Some(value))
    }
}
fn volume(text: &str, unit: &str) -> Result<Option<units::Volume>, String> {
    let unit = units::VolumeUnit::parse(unit).map_err(|e| e.to_string())?;
    optional_number(text)?
        .map(|v| units::Volume::from_unit(v, unit).map_err(|e| e.to_string()))
        .transpose()
}
impl Form {
    pub fn add_reaction(&mut self) -> Result<(), String> {
        if self.reactions.len() >= 128 {
            return Err("At most128 reaction declarations.".into());
        }
        let id = self.allocate_id("reaction")?;
        self.reactions.push(ReactionForm {
            id,
            step: String::new(),
            total: String::new(),
            unit: "uL".into(),
            fill: false,
            fill_id: String::new(),
            fill_label: String::new(),
            components: vec![],
        });
        Ok(())
    }
    pub fn add_component(&mut self, index: usize) -> Result<(), String> {
        if self
            .reactions
            .iter()
            .map(|r| r.components.len())
            .sum::<usize>()
            >= 4096
        {
            return Err("At most4096 declared components.".into());
        }
        let row = self.reactions.get(index).ok_or("Reaction unavailable")?;
        if row.components.len() >= 128 {
            return Err("At most128 components per reaction.".into());
        }
        let id = self.allocate_id("component")?;
        self.reactions[index]
            .components
            .push(ComponentForm::blank(id));
        Ok(())
    }
    pub fn request(&self, binding: &Binding) -> Result<native::Request, String> {
        if self.reactions.len() > 128
            || self
                .reactions
                .iter()
                .map(|r| r.components.len())
                .sum::<usize>()
                > 4096
            || self.name.len() > 256
        {
            return Err("Quantity declaration bounds exceeded.".into());
        }
        for row in &self.reactions {
            if row.components.len() > 128
                || [&row.total, &row.step, &row.unit]
                    .iter()
                    .any(|v| v.len() > 64)
                || row.fill_id.len() > 128
                || row.fill_label.len() > 256
            {
                return Err("Reaction fields exceed native bounds.".into());
            }
            for c in &row.components {
                if c.label.len() > 256
                    || c.source.len() > 1024
                    || [
                        &c.stock,
                        &c.target,
                        &c.volume,
                        &c.dilution_volume,
                        &c.dilution_unit,
                        &c.increment,
                        &c.mode,
                        &c.unit,
                        &c.rounding,
                        &c.concentration_unit,
                    ]
                    .iter()
                    .any(|v| v.len() > 64)
                {
                    return Err("Component fields exceed native bounds.".into());
                }
            }
        }
        let mut reactions = Vec::new();
        for row in &self.reactions {
            let step_index = row
                .step
                .parse::<usize>()
                .map_err(|_| "Choose an explicit operation for every reaction.")?;
            if step_index >= binding.step_count {
                return Err("Reaction operation no longer exists.".into());
            }
            let final_volume = volume(&row.total, &row.unit)?;
            let mut components = Vec::new();
            for c in &row.components {
                let source = if c.source.is_empty() {
                    None
                } else if let Some(id) = c.source.strip_prefix("input:") {
                    if !binding.inputs.iter().any(|i| i.input_id == id) {
                        return Err("Planned input identity is unavailable.".into());
                    }
                    Some(Source::Input(id.into()))
                } else if let Some(index) = c.source.strip_prefix("step:") {
                    let index = index
                        .parse::<usize>()
                        .map_err(|_| "Choose a prior operation reference")?;
                    if index >= step_index {
                        return Err("Component source must precede the reaction operation.".into());
                    }
                    Some(Source::Step(index))
                } else {
                    return Err("Unknown planned material reference.".into());
                };
                if let Some(selected) = &source {
                    if !binding
                        .component_sources
                        .get(step_index)
                        .is_some_and(|choices| choices.iter().any(|s| same_source(s, selected)))
                    {
                        return Err(
                            "Component material is not an operand of this operation.".into()
                        );
                    }
                }
                let quantization = if c.mode == "direct" || c.rounding.is_empty() {
                    None
                } else {
                    let increment = volume(&c.increment, &c.unit)?
                        .ok_or("Declare the rounding increment explicitly.")?;
                    let mode = match c.rounding.as_str() {
                        "exact" => units::RoundingMode::ExactOnly,
                        "floor" => units::RoundingMode::Floor,
                        "ceiling" => units::RoundingMode::Ceiling,
                        "nearest" => units::RoundingMode::NearestTiesUp,
                        _ => return Err("Choose an explicit rounding mode.".into()),
                    };
                    Some(native::Quantization { increment, mode })
                };
                let recipe_volume = if c.mode == "direct" {
                    None
                } else {
                    volume(&c.dilution_volume, &c.dilution_unit)?
                };
                let spec = match c.mode.as_str() {
                    "direct" => native::VolumeSpec::Direct {
                        volume: volume(&c.volume, &c.unit)?,
                    },
                    "mass" => {
                        let parse = |v: &str| {
                            optional_number(v)?
                                .map(|x| {
                                    units::MassConcentration::from_units(
                                        x,
                                        units::MassUnit::Nanogram,
                                        match c.concentration_unit.as_str() {
                                            "ng/uL" => units::VolumeUnit::Microliter,
                                            "ng/nL" => units::VolumeUnit::Nanoliter,
                                            _ => return Err(
                                                "Choose ng/µL or ng/nL mass concentration units."
                                                    .into(),
                                            ),
                                        },
                                    )
                                    .map_err(|e| e.to_string())
                                })
                                .transpose()
                        };
                        native::VolumeSpec::MassDilution {
                            recipe: units::DilutionRecipe {
                                schema_version: units::SCHEMA_VERSION,
                                stock: parse(&c.stock)?,
                                target: parse(&c.target)?,
                                final_volume: recipe_volume,
                            },
                            quantization,
                        }
                    }
                    "molar" => {
                        let parse = |v: &str| {
                            optional_number(v)?
                                .map(|x| {
                                    let value = match c.concentration_unit.as_str() {
                                        "nM" => units::MolarConcentration::from_unit(
                                            x,
                                            units::MolarUnit::Nanomolar,
                                        ),
                                        "uM" => units::MolarConcentration::from_unit(
                                            x,
                                            units::MolarUnit::Micromolar,
                                        ),
                                        "pmol/nL" => units::MolarConcentration::from_amount_volume(
                                            units::Amount::from_unit(
                                                x,
                                                units::AmountUnit::Picomole,
                                            )
                                            .map_err(|e| e.to_string())?,
                                            units::Volume::from_unit(
                                                units::Rational::ONE,
                                                units::VolumeUnit::Nanoliter,
                                            )
                                            .map_err(|e| e.to_string())?,
                                        ),
                                        _ => {
                                            return Err(
                                                "Choose nM, µM or pmol/nL concentration units."
                                                    .into(),
                                            );
                                        }
                                    };
                                    value.map_err(|e| e.to_string())
                                })
                                .transpose()
                        };
                        native::VolumeSpec::MolarDilution {
                            recipe: units::DilutionRecipe {
                                schema_version: units::SCHEMA_VERSION,
                                stock: parse(&c.stock)?,
                                target: parse(&c.target)?,
                                final_volume: recipe_volume,
                            },
                            quantization,
                        }
                    }
                    _ => {
                        return Err("Choose direct volume, mass dilution or molar dilution.".into());
                    }
                };
                components.push(native::Component {
                    id: c.id.clone(),
                    label: c.label.clone(),
                    source,
                    volume: spec,
                });
            }
            reactions.push(native::Reaction {
                id: row.id.clone(),
                step_index,
                final_volume,
                fill: if row.fill {
                    native::FillPolicy::ToFinal {
                        component_id: row.fill_id.clone(),
                        label: row.fill_label.clone(),
                    }
                } else {
                    native::FillPolicy::Unallocated
                },
                components,
            });
        }
        let request = native::Request {
            name: self.name.clone(),
            plan_source: binding.plan_source.clone(),
            plan_revision: binding.plan_revision.clone(),
            target: binding.target.clone(),
            inputs: binding.inputs.clone(),
            reactions,
        };
        key(&request)?;
        Ok(request)
    }
}
impl Form {
    fn allocate_id(&mut self, prefix: &str) -> Result<String, String> {
        loop {
            let number = self.next_id;
            self.next_id = number
                .checked_add(1)
                .ok_or("Stable identifier bound exceeded")?;
            let id = format!("{prefix}-{number}");
            if !self
                .reactions
                .iter()
                .any(|r| r.id == id || r.fill_id == id || r.components.iter().any(|c| c.id == id))
            {
                return Ok(id);
            }
        }
    }
    /// Restore declarations exactly, including missing values and diagnostic
    /// dilution volumes. Canonical-unit fractions avoid rounding or overflow.
    pub fn from_request(request: &native::Request, binding: &Binding) -> Result<Self, String> {
        if !request_matches(binding, request)? {
            return Err(
                "Saved quantities belong to another exact plan/source/target binding.".into(),
            );
        }
        let volume_text =
            |v: Option<units::Volume>| v.map(|v| v.canonical().to_string()).unwrap_or_default();
        let mut result = Self {
            name: request.name.clone(),
            reactions: vec![],
            next_id: 1,
        };
        for reaction in &request.reactions {
            let (fill, fill_id, fill_label) = match &reaction.fill {
                native::FillPolicy::Unallocated => (false, String::new(), String::new()),
                native::FillPolicy::ToFinal {
                    component_id,
                    label,
                } => (true, component_id.clone(), label.clone()),
            };
            let mut components = Vec::new();
            for component in &reaction.components {
                let mut c = ComponentForm::blank(component.id.clone());
                c.label = component.label.clone();
                c.source = component
                    .source
                    .as_ref()
                    .map(source_text)
                    .unwrap_or_default();
                c.unit = "nL".into();
                c.dilution_unit = "nL".into();
                let quantization = match &component.volume {
                    native::VolumeSpec::Direct { volume } => {
                        c.volume = volume_text(*volume);
                        None
                    }
                    native::VolumeSpec::MassDilution {
                        recipe,
                        quantization,
                    } => {
                        c.mode = "mass".into();
                        c.concentration_unit = "ng/nL".into();
                        c.stock = recipe
                            .stock
                            .map(|v| v.canonical().to_string())
                            .unwrap_or_default();
                        c.target = recipe
                            .target
                            .map(|v| v.canonical().to_string())
                            .unwrap_or_default();
                        c.dilution_volume = volume_text(recipe.final_volume);
                        quantization.as_ref()
                    }
                    native::VolumeSpec::MolarDilution {
                        recipe,
                        quantization,
                    } => {
                        c.mode = "molar".into();
                        c.concentration_unit = "pmol/nL".into();
                        c.stock = recipe
                            .stock
                            .map(|v| v.canonical().to_string())
                            .unwrap_or_default();
                        c.target = recipe
                            .target
                            .map(|v| v.canonical().to_string())
                            .unwrap_or_default();
                        c.dilution_volume = volume_text(recipe.final_volume);
                        quantization.as_ref()
                    }
                };
                if let Some(q) = quantization {
                    c.increment = q.increment.canonical().to_string();
                    c.rounding = match q.mode {
                        units::RoundingMode::ExactOnly => "exact",
                        units::RoundingMode::Floor => "floor",
                        units::RoundingMode::Ceiling => "ceiling",
                        units::RoundingMode::NearestTiesUp => "nearest",
                    }
                    .into();
                }
                components.push(c);
            }
            result.reactions.push(ReactionForm {
                id: reaction.id.clone(),
                step: reaction.step_index.to_string(),
                total: volume_text(reaction.final_volume),
                unit: "nL".into(),
                fill,
                fill_id,
                fill_label,
                components,
            });
        }
        if key(&result.request(binding)?)? != key(request)? {
            return Err(
                "Saved quantity declarations cannot be represented exactly by this form.".into(),
            );
        }
        Ok(result)
    }
    /// Copy only explicitly selected rows. Amounts remain declarations; operands
    /// map by their native operation positions, never by sequence equality.
    pub fn apply_to_reactions(
        &mut self,
        binding: &Binding,
        source_id: &str,
        destinations: &[String],
    ) -> Result<(), String> {
        if destinations.is_empty() {
            return Err("Select at least one destination reaction.".into());
        }
        let source = self
            .reactions
            .iter()
            .find(|r| r.id == source_id)
            .ok_or("Choose an existing source reaction.")?
            .clone();
        let step = source
            .step
            .parse::<usize>()
            .map_err(|_| "Choose the source reaction operation.")?;
        let kind = binding
            .operation_kinds
            .get(step)
            .ok_or("Source operation is unavailable.")?;
        let operands = binding
            .component_sources
            .get(step)
            .ok_or("Source operation operands are unavailable.")?;
        let mut seen = std::collections::BTreeSet::new();
        let mut seen_steps = std::collections::BTreeSet::new();
        let mut next = self.clone();
        for id in destinations {
            if id == source_id || !seen.insert(id) {
                return Err(
                    "Choose distinct destination reactions separate from the source.".into(),
                );
            }
            let index = next
                .reactions
                .iter()
                .position(|r| &r.id == id)
                .ok_or("A selected destination reaction disappeared.")?;
            let target_step = next.reactions[index]
                .step
                .parse::<usize>()
                .map_err(|_| "Choose every destination operation explicitly.")?;
            if step == target_step
                || !seen_steps.insert(target_step)
                || binding.operation_kinds.get(target_step) != Some(kind)
            {
                return Err(
                    "Selected reactions must reference distinct operations of the same type."
                        .into(),
                );
            }
            let target_operands = binding
                .component_sources
                .get(target_step)
                .filter(|v| v.len() == operands.len())
                .ok_or("The destination operation has a different operand layout.")?;
            let mut components = source.components.clone();
            for c in &mut components {
                if !c.source.is_empty() {
                    let mut mapped = std::collections::BTreeSet::new();
                    for (position, operand) in operands.iter().enumerate() {
                        if source_text(operand) == c.source {
                            mapped.insert(source_text(&target_operands[position]));
                        }
                    }
                    if mapped.len() != 1 {
                        return Err("A component source has no unique positional match in the selected destination.".into());
                    }
                    c.source = mapped.into_iter().next().unwrap();
                }
                c.id = next.allocate_id("component")?;
            }
            let mut copy = source.clone();
            copy.id = next.reactions[index].id.clone();
            copy.step = target_step.to_string();
            copy.components = components;
            if copy.fill {
                copy.fill_id = next.allocate_id("fill")?;
            }
            next.reactions[index] = copy;
        }
        if next
            .reactions
            .iter()
            .map(|r| r.components.len())
            .sum::<usize>()
            > 4096
        {
            return Err("At most4096 declared components.".into());
        }
        *self = next;
        Ok(())
    }
}
fn source_text(source: &Source) -> String {
    match source {
        Source::Input(id) => format!("input:{id}"),
        Source::Step(i) => format!("step:{i}"),
    }
}

pub fn amount(volume: units::Volume) -> String {
    volume
        .in_unit(units::VolumeUnit::Microliter)
        .map(|r| format!("{r} µL (exact)"))
        .unwrap_or_else(|e| e.to_string())
}

pub fn summary_matches(binding: &Binding, summary: &native::Summary) -> Result<bool, String> {
    Ok(key(&(
        &summary.plan_source,
        &summary.plan_revision,
        &summary.target,
        &summary.inputs,
    ))? == key(&(
        &binding.plan_source,
        &binding.plan_revision,
        &binding.target,
        &binding.inputs,
    ))?)
}

pub fn same_source(a: &Source, b: &Source) -> bool {
    match (a, b) {
        (Source::Input(a), Source::Input(b)) => a == b,
        (Source::Step(a), Source::Step(b)) => a == b,
        _ => false,
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    use collagene_core::{
        documents::{history::DocumentHistory, transfer::Candidate},
        material::lineage as run_lineage,
        planning::plan::{DesignPlan, Operation, PlanInput},
        sequence::analysis::Primer,
        verification,
    };
    use collagene_model::{Molecule, Topology};
    fn fixture() -> (Workspace, DesignPlan, Binding) {
        let source = Molecule::dna("Stock", "AGTCCGATGCTACGTA", Topology::Linear).unwrap();
        let target =
            Molecule::dna("Independent target", "AGTCCGATGCTACGTA", Topology::Linear).unwrap();
        let mut w =
            Workspace::from_documents("Quantity fixture", vec![source.clone(), target.clone()], 1)
                .unwrap();
        let pin = |i: usize| {
            let d = &w.documents[i];
            EventBinding {
                document_id: d.id.clone(),
                event_id: d.head_event.clone().unwrap(),
                revision: d.head.clone(),
            }
        };
        let source_pin = pin(0);
        let target_pin = pin(1);
        let plan = DesignPlan {
            name: "Actual PCR".into(),
            target_revision: target.revision(),
            inputs: vec![PlanInput {
                id: "template".into(),
                revision: source.revision(),
                molecule: source,
            }],
            steps: vec![Operation::Pcr {
                template_strand: collagene_core::chemistry::Strand::Top,
                template: Source::Input("template".into()),
                forward: Primer::new("F", "AGTCC", "").unwrap(),
                reverse: Primer::new("R", "TACGT", "").unwrap(),
            }],
        };
        let candidate = Candidate::new(plan.clone());
        let revision = candidate.plan_revision.clone();
        w.candidates.push(candidate);
        w.validate().unwrap();
        let binding = Binding {
            plan_source: native::PlanSource::Candidate,
            plan_revision: revision,
            target: target_pin,
            inputs: vec![InputEventBinding {
                input_id: "template".into(),
                event: source_pin,
            }],
            step_count: plan.steps.len(),
            operation_kinds: plan.steps.iter().map(std::mem::discriminant).collect(),
            component_sources: plan
                .steps
                .iter()
                .map(|o| {
                    collagene_core::planning::plan::Operation::sources(o)
                        .into_iter()
                        .cloned()
                        .collect()
                })
                .collect(),
        };
        (w, plan, binding)
    }
    fn form() -> Form {
        let mut f = Form::default();
        f.add_reaction().unwrap();
        f.reactions[0].step = "0".into();
        f.add_component(0).unwrap();
        f.reactions[0].components[0].label = "Declared DNA addition".into();
        f.reactions[0].components[0].source = "input:template".into();
        f
    }
    #[test]
    fn blank_amounts_are_missing_and_explicit_mixture_remainder_is_exact() {
        let (w, _, binding) = fixture();
        let mut f = form();
        let blank = f.request(&binding).unwrap();
        assert!(matches!(
            blank.reactions[0].components[0].volume,
            native::VolumeSpec::Direct { volume: None }
        ));
        let preview = native::preview(&w, &blank).unwrap();
        assert!(!preview.complete);
        assert!(preview.reactions[0].components[0].volume.is_none());
        f.reactions[0].total = "20".into();
        f.reactions[0].components[0].volume = "2".into();
        f.reactions[0].fill = true;
        f.reactions[0].fill_id = "declared-fill".into();
        f.reactions[0].fill_label = "Explicitly chosen fill".into();
        let resolved = native::preview(&w, &f.request(&binding).unwrap()).unwrap();
        assert_eq!(
            resolved.reactions[0]
                .fill
                .as_ref()
                .unwrap()
                .volume
                .in_unit(units::VolumeUnit::Microliter)
                .unwrap(),
            units::Rational::integer(18)
        );
        f.reactions[0].components[0].volume = "21".into();
        let over = native::preview(&w, &f.request(&binding).unwrap()).unwrap();
        assert!(!over.complete);
        assert!(!over.reactions[0].findings.is_empty());
        for invalid in ["NaN", "1e3", "-1"] {
            f.reactions[0].components[0].volume = invalid.into();
            assert!(f.request(&binding).is_err());
        }
        f.reactions[0].components[0].volume = "1".into();
        f.reactions[0].components[0].source = "step:0".into();
        assert!(f.request(&binding).is_err());
    }
    #[test]
    fn explicit_quantization_keeps_achieved_mass_concentration_and_inactive_raw_text() {
        let (w, _, binding) = fixture();
        let mut f = form();
        let row = &mut f.reactions[0];
        row.total = "20".into();
        let c = &mut row.components[0];
        c.mode = "mass".into();
        c.concentration_unit = "ng/uL".into();
        c.dilution_volume = "20".into();
        c.stock = "7".into();
        c.target = "3".into();
        c.rounding = "floor".into();
        c.increment = "1".into();
        c.volume = "unfinished".into();
        let p = native::preview(&w, &f.request(&binding).unwrap()).unwrap();
        let Some(native::DilutionResult::QuantizedMass(result)) =
            &p.reactions[0].components[0].calculation
        else {
            panic!("missing exact mass dilution")
        };
        assert_eq!(
            result
                .stock
                .rounded
                .in_unit(units::VolumeUnit::Microliter)
                .unwrap(),
            units::Rational::integer(8)
        );
        assert!(!result.target_reached);
        assert_eq!(
            result
                .achieved
                .in_units(units::MassUnit::Nanogram, units::VolumeUnit::Microliter)
                .unwrap(),
            units::Rational::new(14, 5).unwrap()
        );
        assert_eq!(f.reactions[0].components[0].volume, "unfinished");
        f.reactions[0].components[0].mode = "direct".into();
        assert!(f.request(&binding).is_err());
    }
    #[test]
    fn saved_fractional_declarations_round_trip_without_filling_missing_or_diagnostic_values() {
        let (_, _, binding) = fixture();
        let mut request = form().request(&binding).unwrap();
        let r = &mut request.reactions[0];
        let vol = |n, d| {
            units::Volume::from_unit(
                units::Rational::new(n, d).unwrap(),
                units::VolumeUnit::Nanoliter,
            )
            .unwrap()
        };
        r.final_volume = Some(vol(1, 3));
        r.fill = native::FillPolicy::ToFinal {
            component_id: "retained-fill".into(),
            label: "A separately declared fill".into(),
        };
        let mass = |n, d| {
            units::MassConcentration::from_units(
                units::Rational::new(n, d).unwrap(),
                units::MassUnit::Nanogram,
                units::VolumeUnit::Nanoliter,
            )
            .unwrap()
        };
        r.components = vec![
            native::Component {
                id: "component-2".into(),
                label: "Long label retained exactly ".repeat(4),
                source: Some(Source::Input("template".into())),
                volume: native::VolumeSpec::MassDilution {
                    recipe: units::DilutionRecipe {
                        schema_version: units::SCHEMA_VERSION,
                        stock: Some(mass(u64::MAX, 1)),
                        target: Some(mass(1, 7)),
                        final_volume: None,
                    },
                    quantization: Some(native::Quantization {
                        increment: vol(1, 11),
                        mode: units::RoundingMode::Floor,
                    }),
                },
            },
            native::Component {
                id: "component-3".into(),
                label: "Missing direct addition".into(),
                source: None,
                volume: native::VolumeSpec::Direct { volume: None },
            },
            native::Component {
                id: "molar".into(),
                label: "Canonical molar declaration".into(),
                source: None,
                volume: native::VolumeSpec::MolarDilution {
                    recipe: units::DilutionRecipe {
                        schema_version: units::SCHEMA_VERSION,
                        stock: Some(
                            units::MolarConcentration::from_amount_volume(
                                units::Amount::from_unit(
                                    units::Rational::integer(u64::MAX),
                                    units::AmountUnit::Picomole,
                                )
                                .unwrap(),
                                vol(1, 1),
                            )
                            .unwrap(),
                        ),
                        target: None,
                        final_volume: Some(vol(2, 3)),
                    },
                    quantization: None,
                },
            },
        ];
        let restored = Form::from_request(&request, &binding).unwrap();
        assert_eq!(
            key(&restored.request(&binding).unwrap()).unwrap(),
            key(&request).unwrap()
        );
        assert_eq!(restored.reactions[0].total, "1/3");
        assert!(
            restored.reactions[0].components[0]
                .dilution_volume
                .is_empty()
        );
        assert_eq!(restored.reactions[0].components[2].dilution_volume, "2/3");
        let mut f = restored;
        f.add_component(0).unwrap();
        let ids = f.reactions[0]
            .components
            .iter()
            .map(|c| c.id.as_str())
            .collect::<std::collections::BTreeSet<_>>();
        assert_eq!(ids.len(), 4);
        for text in ["1/0", "1/2/3", "-1/2", "NaN", "1e3"] {
            assert!(optional_number(text).is_err());
        }
    }
    #[test]
    fn revising_retained_quantities_creates_a_new_report_and_preserves_original_attachment_and_heads()
     {
        let (w, _, binding) = fixture();
        let mut initial = form();
        initial.reactions[0].total = "20".into();
        initial.reactions[0].components[0].volume = "2".into();
        let first = native::capture(&w, &initial.request(&binding).unwrap()).unwrap();
        let frame = collagene_core::documents::workspace::load_workspace(
            &collagene_core::documents::workspace::save_workspace(&first.workspace).unwrap(),
        )
        .unwrap();
        let original_bytes = frame.attachments[&first.attachment_hash].bytes.clone();
        let heads = key(&frame.documents).unwrap();
        let opened = native::open(&frame, &first.attachment_hash).unwrap();
        let mut revised = Form::from_request(&opened.preview.request, &binding).unwrap();
        revised.reactions[0].components[0].volume = "3000".into();
        let second = native::capture(&frame, &revised.request(&binding).unwrap()).unwrap();
        assert_ne!(first.attachment_hash, second.attachment_hash);
        assert_eq!(
            second.workspace.attachments[&first.attachment_hash].bytes,
            original_bytes
        );
        assert_eq!(key(&second.workspace.documents).unwrap(), heads);
        let old = native::open(&second.workspace, &first.attachment_hash).unwrap();
        let new = native::open(&second.workspace, &second.attachment_hash).unwrap();
        assert_eq!(
            old.preview.reactions[0].components[0]
                .volume
                .unwrap()
                .canonical(),
            units::Rational::integer(2000)
        );
        assert_eq!(
            new.preview.reactions[0].components[0]
                .volume
                .unwrap()
                .canonical(),
            units::Rational::integer(3000)
        );
        let mut wrong = binding;
        wrong.target.event_id = "different-event".into();
        assert!(Form::from_request(&opened.preview.request, &wrong).is_err());
    }
    #[test]
    fn selected_reaction_copy_maps_native_operands_and_rejects_incompatible_destinations_atomically()
     {
        let (_, mut plan, mut binding) = fixture();
        let mut second_input = plan.inputs[0].clone();
        second_input.id = "second".into();
        plan.inputs.push(second_input);
        let mut second = plan.steps[0].clone();
        let Operation::Pcr { template, .. } = &mut second else {
            unreachable!()
        };
        *template = Source::Input("second".into());
        plan.steps.push(second);
        plan.steps.push(Operation::TreatEnds {
            template: Source::Step(1),
            treatment: collagene_core::chemistry::EndTreatment::Phosphorylate,
        });
        binding.plan_revision = collagene_core::planning::plan::DesignPlan::revision(&plan);
        binding.step_count = 3;
        binding.component_sources = plan
            .steps
            .iter()
            .map(|op| {
                collagene_core::planning::plan::Operation::sources(op)
                    .into_iter()
                    .cloned()
                    .collect()
            })
            .collect();
        binding.operation_kinds = plan.steps.iter().map(std::mem::discriminant).collect();
        let mut pin = binding.inputs[0].clone();
        pin.input_id = "second".into();
        binding.inputs.push(pin);
        let mut f = form();
        f.reactions[0].total = "20".into();
        f.reactions[0].components[0].volume = "2".into();
        f.reactions[0].fill = true;
        f.reactions[0].fill_id = "source-fill".into();
        f.reactions[0].fill_label = "Declared fill".into();
        let source = f.reactions[0].id.clone();
        f.add_reaction().unwrap();
        f.reactions[1].step = "1".into();
        let destination = f.reactions[1].id.clone();
        f.add_reaction().unwrap();
        f.reactions[2].step = "2".into();
        let incompatible = f.reactions[2].id.clone();
        let before = f.clone();
        assert!(
            f.apply_to_reactions(&binding, &source, &[destination.clone(), incompatible])
                .is_err()
        );
        assert!(f == before);
        f.apply_to_reactions(&binding, &source, &[destination.clone()])
            .unwrap();
        assert_eq!(f.reactions[1].id, destination);
        assert_eq!(f.reactions[1].components[0].source, "input:second");
        assert_eq!(f.reactions[1].components[0].volume, "2");
        assert_ne!(
            f.reactions[1].components[0].id,
            f.reactions[0].components[0].id
        );
        assert_ne!(f.reactions[1].fill_id, f.reactions[0].fill_id);
        assert!(f.reactions[0] == before.reactions[0]);
        f.request(&binding).unwrap();
        let before = f.clone();
        let mut ambiguous = binding;
        ambiguous.component_sources[0] = vec![
            Source::Input("template".into()),
            Source::Input("template".into()),
        ];
        ambiguous.component_sources[1] = vec![
            Source::Input("template".into()),
            Source::Input("second".into()),
        ];
        assert!(
            f.apply_to_reactions(&ambiguous, &source, &[destination])
                .is_err()
        );
        assert!(f == before);
    }
    #[test]
    fn saved_run_quantities_keep_exact_retained_events_and_native_archive_without_moving_heads() {
        let (w, plan, mut binding) = fixture();
        let kept = run_lineage::record_run(
            &w,
            &plan,
            binding.target.clone(),
            binding.inputs.clone(),
            false,
        )
        .unwrap();
        binding.plan_source = native::PlanSource::SavedRun {
            run_id: kept.receipt.run_id,
        };
        let mut frame = kept.workspace;
        let mut docs = DocumentHistory::restore_all(&frame).unwrap();
        let alias = docs[0]
            .fork_event(&binding.inputs[0].event.event_id, "equal-snapshot-alias")
            .unwrap();

        let alias_pin = EventBinding {
            document_id: alias.id.clone(),
            event_id: alias.head_event().into(),
            revision: alias.head_revision().into(),
        };
        for d in docs.iter_mut().take(2) {
            let mut private = d.molecule.clone();
            private.name.push_str(" private working head");
            d.commit(private, "Unpublished working branch").unwrap();
        }
        frame.documents.clear();
        frame.history.clear();
        for document in &docs {
            document.contribute(&mut frame).unwrap();
        }
        alias.contribute(&mut frame).unwrap();
        frame.validate().unwrap();
        let scope = Scope::capture(&frame, &binding, 5).unwrap();
        let form = form();
        let request = form.request(&binding).unwrap();
        let historical = native::preview(&frame, &request).unwrap();
        historical
            .validate_binding(&native::preview_revision(&frame, &request).unwrap())
            .unwrap();
        let mut wrong = binding.clone();
        wrong.inputs[0].event = alias_pin;
        assert!(scope.permit(&frame, &wrong, 5).is_err());
        assert!(scope.permit(&frame, &binding, 6).is_err());
        assert!(native::preflight(&frame, &form.request(&wrong).unwrap()).is_err());
        let before = frame.documents.clone();
        let captured = native::capture(&frame, &request).unwrap();
        captured
            .validate_binding(&historical.context_revision)
            .unwrap();
        let reopened = collagene_core::documents::workspace::load_workspace(
            &collagene_core::documents::workspace::save_workspace(&captured.workspace).unwrap(),
        )
        .unwrap();
        assert_eq!(
            serde_json::to_vec(&reopened.documents).unwrap(),
            serde_json::to_vec(&before).unwrap()
        );
        let opened = native::open(&reopened, &captured.attachment_hash).unwrap();
        opened
            .validate_binding(&native::open_revision(&reopened, &captured.attachment_hash).unwrap())
            .unwrap();
        assert_eq!(
            opened.preview.request.target.event_id,
            binding.target.event_id
        );
        assert!(!opened.preview.complete);
        let listed = native::list(&reopened).unwrap();
        assert!(
            listed
                .entries
                .iter()
                .any(|s| s.attachment_hash == captured.attachment_hash
                    && summary_matches(&binding, s).unwrap())
        );
    }
}

pub fn request_matches(binding: &Binding, request: &native::Request) -> Result<bool, String> {
    Ok(key(&(
        &request.plan_source,
        &request.plan_revision,
        &request.target,
        &request.inputs,
    ))? == key(&(
        &binding.plan_source,
        &binding.plan_revision,
        &binding.target,
        &binding.inputs,
    ))?)
}

pub fn same_binding(a: &Binding, b: &Binding) -> Result<bool, String> {
    Ok(key(a)? == key(b)?)
}
