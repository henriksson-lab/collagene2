//! Documentary numeric declarations with exact native event identity.
use collagene_core::{
    protocol_quantities::{self as native, units},
    run_lineage::{EventBinding, InputEventBinding},
    verification::Source,
    workspace::Workspace,
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
            fingerprint: collagene_core::run_lineage::workspace_fingerprint(w)
                .map_err(|e| e.to_string())?,
            binding: key(binding)?,
        })
    }
    pub fn permit(&self, w: &Workspace, binding: &Binding, epoch: u64) -> Result<(), String> {
        if epoch != self.epoch
            || key(binding)? != self.binding
            || collagene_core::run_lineage::workspace_fingerprint(w).map_err(|e| e.to_string())?
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
        units::Rational::parse_decimal(text.trim())
            .map(Some)
            .map_err(|e| e.to_string())
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
        let id = self.next_id;
        self.next_id = self
            .next_id
            .checked_add(1)
            .ok_or("Stable identifier bound exceeded")?;
        self.reactions.push(ReactionForm {
            id: format!("reaction-{id}"),
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
        let row = self
            .reactions
            .get_mut(index)
            .ok_or("Reaction unavailable")?;
        if row.components.len() >= 128 {
            return Err("At most128 components per reaction.".into());
        }
        let id = self.next_id;
        self.next_id = self
            .next_id
            .checked_add(1)
            .ok_or("Stable identifier bound exceeded")?;
        row.components
            .push(ComponentForm::blank(format!("component-{id}")));
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
            || self.name.chars().count() > 64
        {
            return Err("Quantity declaration bounds exceeded.".into());
        }
        for row in &self.reactions {
            if row.components.len() > 128
                || [&row.total, &row.step, &row.unit]
                    .iter()
                    .any(|v| v.len() > 64)
                || row.fill_id.len() > 64
                || row.fill_label.chars().count() > 64
            {
                return Err("Reaction fields exceed native bounds.".into());
            }
            for c in &row.components {
                if c.label.chars().count() > 64
                    || c.source.len() > 256
                    || [
                        &c.stock,
                        &c.target,
                        &c.volume,
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
                let recipe_volume = final_volume;
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
                                        units::VolumeUnit::Microliter,
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
                        let unit = match c.concentration_unit.as_str() {
                            "nM" => units::MolarUnit::Nanomolar,
                            "uM" => units::MolarUnit::Micromolar,
                            _ => return Err("Choose nM or µM concentration units.".into()),
                        };
                        let parse = |v: &str| {
                            optional_number(v)?
                                .map(|x| {
                                    units::MolarConcentration::from_unit(x, unit)
                                        .map_err(|e| e.to_string())
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
        Molecule, Topology,
        analysis::Primer,
        history::DocumentHistory,
        run_lineage,
        transfer::Candidate,
        verification::{self, DesignPlan, Operation, PlanInput},
    };
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
            component_sources: plan
                .steps
                .iter()
                .map(|o| {
                    verification::operation_sources(o)
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
        let reopened = collagene_core::workspace::load_workspace(
            &collagene_core::workspace::save_workspace(&captured.workspace).unwrap(),
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
