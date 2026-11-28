//! Exact declared reaction quantities bound to retained plans and native events.
//! Arithmetic documents confer no chemistry authority. Combined protocol exports
//! freshly verify their retained graph; no measured stock or execution is inferred.
use crate::{
    documents::workspace::{Workspace, WorkspaceAttachment},
    material::lineage::{self as run_lineage, EventBinding, InputEventBinding, RunRecord},
    planning::plan::{DesignPlan, Source},
    verification,
};
use collagene_model::Error;
pub use collagene_protocol_units as units;
use serde::{Deserialize, Serialize};
use std::{
    collections::BTreeSet,
    io::{self, Write},
};
pub const KIND: &str = "protocol-quantities-v1";
pub const ENGINE_VERSION: &str = "protocol-quantities-1";
pub const MAX_BYTES: usize = 4 * 1024 * 1024;
pub const MAX_EXPORT_BYTES: usize = 16 * 1024 * 1024;
const MAX_EXPORT_ENVELOPE_BYTES: usize = 32 * 1024 * 1024 - 64 * 1024;
pub const MAX_RECORDS: usize = 128;
pub const MAX_REACTIONS: usize = 128;
pub const MAX_COMPONENTS: usize = 4096;
macro_rules! dto{($($t:tt)*)=>{#[derive(Debug,Clone,Serialize,Deserialize)]#[cfg_attr(feature="json-schema",derive(schemars::JsonSchema))]#[serde(deny_unknown_fields)]$($t)*};}
dto! {#[derive(PartialEq,Eq)]pub enum PlanSource{Candidate,SavedRun{run_id:String}}}
dto! {pub struct Quantization{#[cfg_attr(feature="json-schema",schemars(with="serde_json::Value"))]pub increment:units::Volume,#[cfg_attr(feature="json-schema",schemars(with="serde_json::Value"))]pub mode:units::RoundingMode}}
dto! {pub enum VolumeSpec{
 Direct{#[serde(deserialize_with="crate::cloning::assembly::recipe::required_option")]#[cfg_attr(feature="json-schema",schemars(with="serde_json::Value",required))]volume:Option<units::Volume>},
 MassDilution{#[cfg_attr(feature="json-schema",schemars(with="serde_json::Value"))]recipe:units::DilutionRecipe<units::MassConcentration>,#[serde(deserialize_with="crate::cloning::assembly::recipe::required_option")]#[cfg_attr(feature="json-schema",schemars(required))]quantization:Option<Quantization>},
 MolarDilution{#[cfg_attr(feature="json-schema",schemars(with="serde_json::Value"))]recipe:units::DilutionRecipe<units::MolarConcentration>,#[serde(deserialize_with="crate::cloning::assembly::recipe::required_option")]#[cfg_attr(feature="json-schema",schemars(required))]quantization:Option<Quantization>},
}}
dto! {pub struct Component{pub id:String,pub label:String,#[serde(deserialize_with="crate::cloning::assembly::recipe::required_option")]#[cfg_attr(feature="json-schema",schemars(required))]pub source:Option<Source>,pub volume:VolumeSpec}}
dto! {pub enum FillPolicy{Unallocated,ToFinal{component_id:String,label:String}}}
dto! {pub struct Reaction{pub id:String,pub step_index:usize,#[serde(deserialize_with="crate::cloning::assembly::recipe::required_option")]#[cfg_attr(feature="json-schema",schemars(with="serde_json::Value",required))]pub final_volume:Option<units::Volume>,pub fill:FillPolicy,#[serde(deserialize_with="components")]pub components:Vec<Component>}}
dto! {pub struct Request{pub name:String,pub plan_source:PlanSource,pub plan_revision:String,pub target:EventBinding,#[serde(deserialize_with="inputs")]pub inputs:Vec<InputEventBinding>,#[serde(deserialize_with="reactions")]pub reactions:Vec<Reaction>}}
dto! {#[derive(PartialEq,Eq)]pub struct Finding{pub path:String,pub message:String}}
dto! {pub enum DilutionResult{
 Mass(#[cfg_attr(feature="json-schema",schemars(with="serde_json::Value"))]units::DilutionPlan<units::MassConcentration>),
 Molar(#[cfg_attr(feature="json-schema",schemars(with="serde_json::Value"))]units::DilutionPlan<units::MolarConcentration>),
 QuantizedMass(#[cfg_attr(feature="json-schema",schemars(with="serde_json::Value"))]units::QuantizedDilution<units::MassConcentration>),
 QuantizedMolar(#[cfg_attr(feature="json-schema",schemars(with="serde_json::Value"))]units::QuantizedDilution<units::MolarConcentration>),
}}
dto! {pub struct ComponentReport{pub component_id:String,#[cfg_attr(feature="json-schema",schemars(with="serde_json::Value"))]pub volume:Option<units::Volume>,pub calculation:Option<DilutionResult>,pub findings:Vec<Finding>}}
dto! {pub struct FillVolume{pub component_id:String,pub label:String,#[cfg_attr(feature="json-schema",schemars(with="serde_json::Value"))]pub volume:units::Volume}}
dto! {pub struct ReactionReport{pub id:String,pub step_index:usize,pub components:Vec<ComponentReport>,#[cfg_attr(feature="json-schema",schemars(with="serde_json::Value"))]pub mixture:Option<units::MixturePlan>,pub fill:Option<FillVolume>,pub complete:bool,pub findings:Vec<Finding>}}
dto! {pub struct Preview{pub schema_version:u16,pub engine_version:String,pub units_engine:String,pub context_revision:String,pub source_workspace_fingerprint:String,pub request:Request,pub plan_name:String,pub complete:bool,pub findings:Vec<Finding>,pub warnings:Vec<String>,pub reactions:Vec<ReactionReport>}}
dto! {pub struct Capture{#[cfg_attr(feature="json-schema",schemars(with="serde_json::Value"))]pub workspace:Workspace,pub context_revision:String,pub source_workspace_fingerprint:String,pub attachment_hash:String,pub preview:Preview}}
dto! {pub struct Opened{pub context_revision:String,pub source_workspace_fingerprint:String,pub quantity_hash:String,pub preview:Preview}}
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum ExportFormat {
    Csv,
    Html,
    Pdf,
    ProtocolHtml,
    ProtocolPdf,
}
dto! {pub struct MolecularVerification{pub engine_version:String,pub plan_revision:String,pub target_revision:String,pub status:verification::evidence::Status,pub report_revision:String,pub findings_total:usize,pub exhaustive_alternative_search:bool}}
dto! {pub struct Exported{pub context_revision:String,#[serde(deserialize_with="crate::cloning::assembly::recipe::required_option")]pub source_workspace_fingerprint:Option<String>,pub quantity_hash:String,pub plan_revision:String,pub target:EventBinding,pub inputs:Vec<InputEventBinding>,pub arithmetic_complete:bool,#[serde(deserialize_with="crate::cloning::assembly::recipe::required_option")]pub molecular:Option<MolecularVerification>,pub format:ExportFormat,pub filename:String,pub media_type:String,pub bytes:Vec<u8>,pub content_hash:String,pub warnings:Vec<String>}}
dto! {pub struct Summary{pub attachment_hash:String,pub name:String,pub plan_source:PlanSource,pub plan_revision:String,pub target:EventBinding,pub inputs:Vec<InputEventBinding>,pub complete:bool,pub reaction_count:usize,pub available:bool,pub findings:Vec<Finding>}}
dto! {pub struct Listed{pub context_revision:String,pub source_workspace_fingerprint:String,pub entries:Vec<Summary>}}
#[derive(Debug, Clone)]
pub struct Dependencies {
    pub events: Vec<EventBinding>,
    pub runs: Vec<String>,
    pub candidates: Vec<String>,
}
fn bad(s: &str) -> Error {
    Error::Format(format!("reaction quantities: {s}"))
}
fn digest(s: &str) -> bool {
    s.len() == 64
        && s.bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}
fn text(s: &str, n: usize) -> bool {
    !s.trim().is_empty() && s.len() <= n && !s.chars().any(char::is_control)
}
struct Bytes(Vec<u8>);
impl Write for Bytes {
    fn write(&mut self, b: &[u8]) -> io::Result<usize> {
        if self
            .0
            .len()
            .checked_add(b.len())
            .is_none_or(|n| n > MAX_BYTES)
        {
            return Err(io::Error::other("quantity document bound"));
        }
        self.0.extend_from_slice(b);
        Ok(b.len())
    }
    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}
fn encode<T: Serialize>(v: &T) -> Result<Vec<u8>, Error> {
    let mut b = Bytes(vec![]);
    serde_json::to_writer(&mut b, v).map_err(|_| Error::TooLarge)?;
    Ok(b.0)
}
fn hash<T: Serialize>(v: &T) -> Result<String, Error> {
    Ok(blake3::hash(&encode(v)?).to_hex().to_string())
}
fn bounded<'de, D: serde::Deserializer<'de>, T: Deserialize<'de>, const N: usize>(
    d: D,
) -> Result<Vec<T>, D::Error> {
    struct V<T, const N: usize>(std::marker::PhantomData<T>);
    impl<'de, T: Deserialize<'de>, const N: usize> serde::de::Visitor<'de> for V<T, N> {
        type Value = Vec<T>;
        fn expecting(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
            write!(f, "at most {N} quantity entries")
        }
        fn visit_seq<A: serde::de::SeqAccess<'de>>(
            self,
            mut a: A,
        ) -> Result<Self::Value, A::Error> {
            if a.size_hint().is_some_and(|n| n > N) {
                return Err(serde::de::Error::custom("quantity entry bound"));
            }
            let mut v = vec![];
            while v.len() < N {
                match a.next_element()? {
                    Some(x) => v.push(x),
                    None => return Ok(v),
                }
            }
            if a.next_element::<serde::de::IgnoredAny>()?.is_some() {
                return Err(serde::de::Error::custom("quantity entry bound"));
            }
            Ok(v)
        }
    }
    d.deserialize_seq(V::<T, N>(std::marker::PhantomData))
}
fn inputs<'de, D: serde::Deserializer<'de>>(d: D) -> Result<Vec<InputEventBinding>, D::Error> {
    bounded::<D, _, 128>(d)
}
fn reactions<'de, D: serde::Deserializer<'de>>(d: D) -> Result<Vec<Reaction>, D::Error> {
    bounded::<D, _, MAX_REACTIONS>(d)
}
fn components<'de, D: serde::Deserializer<'de>>(d: D) -> Result<Vec<Component>, D::Error> {
    bounded::<D, _, 128>(d)
}
fn shape(r: &Request) -> Result<(), Error> {
    if !text(&r.name, 256)
        || !digest(&r.plan_revision)
        || r.inputs.is_empty()
        || r.inputs.len() > 128
        || r.reactions.is_empty()
        || r.reactions.len() > MAX_REACTIONS
    {
        return Err(bad(
            "choose a named retained plan, exact events and1–128 explicit reactions",
        ));
    }
    if let PlanSource::SavedRun { run_id } = &r.plan_source {
        if !digest(run_id) {
            return Err(bad("invalid saved run identity"));
        }
    }
    let mut count = 0usize;
    let mut ids = BTreeSet::new();
    let mut steps = BTreeSet::new();
    for reaction in &r.reactions {
        if !text(&reaction.id, 128)
            || !ids.insert(&reaction.id)
            || reaction.step_index >= 1000
            || !steps.insert(reaction.step_index)
            || reaction.components.len() > 128
        {
            return Err(bad("invalid/duplicate reaction or operation selection"));
        }
        count = count
            .checked_add(reaction.components.len())
            .filter(|n| *n <= MAX_COMPONENTS)
            .ok_or(Error::TooLarge)?;
        let mut components = BTreeSet::new();
        for c in &reaction.components {
            if !text(&c.id, 128) || !text(&c.label, 256) || !components.insert(c.id.as_str()) {
                return Err(bad("invalid/duplicate component identifier or label"));
            }
            match &c.volume {
                VolumeSpec::MassDilution { recipe, .. }
                    if recipe.schema_version != units::SCHEMA_VERSION =>
                {
                    return Err(Error::Unsupported(
                        "quantity dilution recipe version".into(),
                    ));
                }
                VolumeSpec::MolarDilution { recipe, .. }
                    if recipe.schema_version != units::SCHEMA_VERSION =>
                {
                    return Err(Error::Unsupported(
                        "quantity dilution recipe version".into(),
                    ));
                }
                _ => {}
            }
        }
        if let FillPolicy::ToFinal {
            component_id,
            label,
        } = &reaction.fill
        {
            if !text(component_id, 128)
                || !text(label, 256)
                || components.contains(component_id.as_str())
            {
                return Err(bad(
                    "explicit fill must have a distinct identifier and label",
                ));
            }
        }
    }
    encode(r)?;
    Ok(())
}
const MAX_REFERENCED_BYTES: usize = 32 * 1024 * 1024;
struct Scope<'a> {
    w: &'a Workspace,
    index: std::collections::BTreeMap<String, &'a WorkspaceAttachment>,
    runs: std::collections::BTreeMap<String, RunRecord>,
    bytes: usize,
}
impl<'a> Scope<'a> {
    fn new(w: &'a Workspace) -> Result<Self, Error> {
        crate::planning::combinatorial::workspace_guard(w)?;
        Self::validated(w)
    }
    fn validated(w: &'a Workspace) -> Result<Self, Error> {
        let mut index = std::collections::BTreeMap::new();
        // Build discovery indexes through the current native reader, after the
        // caller has validated the complete archive and event cross references.
        for a in w
            .attachments
            .values()
            .filter(|a| a.kind == run_lineage::WORKSPACE_ATTACHMENT_KIND)
        {
            let record = run_lineage::load_run_record(&a.bytes)?;
            if index.insert(record.run_id, a).is_some() {
                return Err(bad("duplicate native run"));
            }
        }
        Ok(Self {
            w,
            index,
            runs: std::collections::BTreeMap::new(),
            bytes: 0,
        })
    }
    fn resolve(&mut self, r: &Request) -> Result<&DesignPlan, Error> {
        shape(r)?;
        let (p, record) = match &r.plan_source {
            PlanSource::Candidate => {
                let c = self
                    .w
                    .candidates
                    .iter()
                    .find(|c| c.plan_revision == r.plan_revision)
                    .ok_or_else(|| bad("the exact candidate is not retained in this workspace"))?;
                (&c.plan, None)
            }
            PlanSource::SavedRun { run_id } => {
                if !self.runs.contains_key(run_id) {
                    let asset = self
                        .index
                        .get(run_id)
                        .ok_or_else(|| bad("native source run is unavailable"))?;
                    self.bytes = self
                        .bytes
                        .checked_add(asset.bytes.len())
                        .filter(|n| *n <= MAX_REFERENCED_BYTES)
                        .ok_or(Error::TooLarge)?;
                    let record = run_lineage::load_run_record(&asset.bytes)?;
                    if record.run_id != *run_id {
                        return Err(bad("native run identity differs"));
                    }
                    self.runs.insert(run_id.clone(), record);
                }
                let record = self
                    .runs
                    .get(run_id)
                    .ok_or_else(|| bad("native source run unavailable"))?;
                (&record.plan, Some(record))
            }
        };
        if p.steps.is_empty()
            || p.inputs.len() != r.inputs.len()
            || p.inputs.len() > 128
            || crate::planning::plan::DesignPlan::revision(p) != r.plan_revision
            || p.target_revision != r.target.revision
        {
            return Err(bad("plan, target or original input identities changed"));
        }
        let target = crate::planning::combinatorial::resolve_event(self.w, &r.target)?;
        if target.revision() != p.target_revision {
            return Err(bad("exact independent target differs from retained plan"));
        }
        if let Some(record) = record {
            if record.target != r.target || encode(&record.inputs)? != encode(&r.inputs)? {
                return Err(bad(
                    "saved run quantities require its exact original source and target events",
                ));
            }
        }
        let mut models = vec![target];
        for (input, b) in p.inputs.iter().zip(&r.inputs) {
            let m = crate::planning::combinatorial::resolve_event(self.w, &b.event)?;
            if input.id != b.input_id || input.revision != b.event.revision || &input.molecule != m
            {
                return Err(bad(
                    "exact input event differs from the retained plan snapshot",
                ));
            }
            models.push(m);
        }
        collagene_model::bounds::model_preflight(&models)?;
        // Candidate/run structural graph validation was performed once by the native
        // workspace validator; quantity references must name actual selected inputs.
        for reaction in &r.reactions {
            let op = p
                .steps
                .get(reaction.step_index)
                .ok_or_else(|| bad("reaction operation is not retained by this plan"))?;
            let actual = crate::planning::plan::Operation::sources(op);
            for c in &reaction.components {
                if let Some(source) = &c.source {
                    if !actual.iter().any(|s| match (source, *s) {
                        (Source::Input(a), Source::Input(b)) => a == b,
                        (Source::Step(a), Source::Step(b)) => a == b,
                        _ => false,
                    }) {
                        return Err(bad(
                            "component material is not an actual input of the selected operation",
                        ));
                    }
                }
            }
        }
        Ok(p)
    }
}
/// Borrowed bounds and native plan/event membership, without chemical replay.
pub fn preflight(w: &Workspace, r: &Request) -> Result<(), Error> {
    shape(r)?;
    Scope::new(w)?.resolve(r).map(|_| ())
}
fn context(fp: &str, r: &Request, name: &str) -> Result<String, Error> {
    hash(&(ENGINE_VERSION, units::ENGINE_VERSION, fp, r, name))
}
pub fn preview_revision(w: &Workspace, r: &Request) -> Result<String, Error> {
    shape(r)?;
    let mut scope = Scope::new(w)?;
    let p = scope.resolve(r)?;
    context(&run_lineage::workspace_fingerprint(w)?, r, &p.name)
}
fn finding(path: &str, message: impl ToString) -> Finding {
    Finding {
        path: path.into(),
        message: message.to_string(),
    }
}
fn component(reaction: &Reaction, c: &Component) -> ComponentReport {
    let path = format!("reactions.{}.components.{}", reaction.id, c.id);
    let mut findings = vec![];
    let calculated: Result<(units::Volume, Option<DilutionResult>), units::Error> = match &c.volume
    {
        VolumeSpec::Direct { volume } => {
            volume.map(|v| (v, None)).ok_or(units::Error::MissingInput {
                field: units::InputField::ComponentVolume,
                index: None,
            })
        }
        VolumeSpec::MassDilution {
            recipe,
            quantization,
        } => {
            if recipe.final_volume.is_some() && recipe.final_volume != reaction.final_volume {
                findings.push(finding(&path,"Dilution final volume must equal this reaction’s explicitly declared final volume."));
            }
            match quantization {
                None => {
                    units::dilute(recipe).map(|p| (p.stock_volume, Some(DilutionResult::Mass(p))))
                }
                Some(q) => units::quantized_dilution(recipe, q.increment, q.mode)
                    .map(|p| (p.stock.rounded, Some(DilutionResult::QuantizedMass(p)))),
            }
        }
        VolumeSpec::MolarDilution {
            recipe,
            quantization,
        } => {
            if recipe.final_volume.is_some() && recipe.final_volume != reaction.final_volume {
                findings.push(finding(&path,"Dilution final volume must equal this reaction’s explicitly declared final volume."));
            }
            match quantization {
                None => {
                    units::dilute(recipe).map(|p| (p.stock_volume, Some(DilutionResult::Molar(p))))
                }
                Some(q) => units::quantized_dilution(recipe, q.increment, q.mode)
                    .map(|p| (p.stock.rounded, Some(DilutionResult::QuantizedMolar(p)))),
            }
        }
    };
    let (volume, calculation) = match calculated {
        Ok((v, p)) if findings.is_empty() => (Some(v), p),
        Ok(_) => (None, None),
        Err(e) => {
            findings.push(finding(&path, e));
            (None, None)
        }
    };
    ComponentReport {
        component_id: c.id.clone(),
        volume,
        calculation,
        findings,
    }
}
fn calculations(r: &Request) -> (Vec<ReactionReport>, Vec<Finding>, bool) {
    let mut reports = vec![];
    let mut all = vec![];
    for reaction in &r.reactions {
        let path = format!("reactions.{}", reaction.id);
        let components = reaction
            .components
            .iter()
            .map(|c| component(reaction, c))
            .collect::<Vec<_>>();
        let mut findings = components
            .iter()
            .flat_map(|c| c.findings.clone())
            .collect::<Vec<_>>();
        let mut mixture = None;
        let mut fill = None;
        match reaction.final_volume {
            None => findings.push(finding(
                &path,
                "Final reaction volume is missing; no fill volume was inferred.",
            )),
            Some(final_volume) => {
                let recipe = units::MixtureRecipe {
                    schema_version: units::SCHEMA_VERSION,
                    final_volume,
                    components: components
                        .iter()
                        .map(|c| units::MixtureComponent {
                            role: c.component_id.clone(),
                            volume: c.volume,
                        })
                        .collect(),
                };
                match units::balance_mixture(&recipe) {
                    Err(e) => findings.push(finding(&path, e)),
                    Ok(plan) => {
                        match &reaction.fill{
      FillPolicy::ToFinal{component_id,label}=>fill=Some(FillVolume{component_id:component_id.clone(),label:label.clone(),volume:plan.remaining_diluent}),
      FillPolicy::Unallocated if !plan.remaining_diluent.is_zero()=>findings.push(finding(&path,"Unallocated remaining volume: choose an explicit fill component or complete the specified mixture. No water/reagent was inferred.")),
      FillPolicy::Unallocated=>{},
     }
                        mixture = Some(plan);
                    }
                }
            }
        }
        let complete = findings.is_empty();
        all.extend(findings.clone());
        reports.push(ReactionReport {
            id: reaction.id.clone(),
            step_index: reaction.step_index,
            components,
            mixture,
            fill,
            complete,
            findings,
        });
    }
    let complete = all.is_empty();
    (reports, all, complete)
}
fn warnings() -> Vec<String> {
    vec!["Declared concentrations, volumes and dilution arithmetic are user inputs, not measured stock availability, inventory or chemical compatibility.".into(),"No molecular verifier, chemical replay, laboratory execution, material isolation, supplier conditions or reaction success is established by this document.".into(), "Arithmetic completeness covers only declared reactions and components. Unlisted operations and omitted reagents remain unspecified, even if an explicit fill alone balances the total.".into(), "Candidate event bindings are caller declarations checked against full retained snapshots, not proof of historical candidate ancestry. SavedRun bindings must match its native original target and input events exactly.".into(),"Dilution stock volumes use C1×V1=C2×V2 with additive volumes and zero-concentration diluent. Fill is included only when explicitly named. Quantization reports achieved concentration and deviation, not instrument accuracy. Explicit rounding can leave the requested concentration unmet even when arithmetic is complete; inspect achieved and target_reached.".into()]
}
pub fn preview(w: &Workspace, r: &Request) -> Result<Preview, Error> {
    shape(r)?;
    let mut scope = Scope::new(w)?;
    let p = scope.resolve(r)?;
    let fp = run_lineage::workspace_fingerprint(w)?;
    let name = p.name.clone();
    let (reactions, findings, complete) = calculations(r);
    let result = Preview {
        schema_version: 1,
        engine_version: ENGINE_VERSION.into(),
        units_engine: units::ENGINE_VERSION.into(),
        context_revision: context(&fp, r, &name)?,
        source_workspace_fingerprint: fp,
        request: r.clone(),
        plan_name: name,
        complete,
        findings,
        warnings: warnings(),
        reactions,
    };
    result.validate_binding(&result.context_revision)?;
    Ok(result)
}
impl Preview {
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        shape(&self.request)?;
        if self.schema_version != 1
            || self.engine_version != ENGINE_VERSION
            || self.units_engine != units::ENGINE_VERSION
        {
            return Err(Error::Unsupported(
                "reaction quantity document version or arithmetic engine".into(),
            ));
        }
        if !digest(&self.source_workspace_fingerprint)
            || self.plan_name.trim().is_empty()
            || self.plan_name.len() > 1024
            || context(
                &self.source_workspace_fingerprint,
                &self.request,
                &self.plan_name,
            )? != expected
            || self.context_revision != expected
        {
            return Err(bad("quantity document context, engine or version changed"));
        }
        let (reactions, findings, complete) = calculations(&self.request);
        if self.complete != complete
            || encode(&self.reactions)? != encode(&reactions)?
            || self.findings != findings
            || self.warnings != warnings()
        {
            return Err(bad(
                "quantity results differ from exact declared arithmetic",
            ));
        }
        encode(self)?;
        Ok(())
    }
    pub fn validate_for_workspace(&self, w: &Workspace) -> Result<(), Error> {
        self.validate_binding(&self.context_revision)?;
        let mut scope = Scope::new(w)?;
        let p = scope.resolve(&self.request)?;
        if p.name != self.plan_name {
            return Err(bad("retained plan identity changed"));
        }
        Ok(())
    }
}
fn records(w: &Workspace) -> Result<(), Error> {
    if w.attachments.len() > crate::documents::workspace::MAX_ATTACHMENTS
        || w.attachments.values().filter(|a| a.kind == KIND).count() > MAX_RECORDS
    {
        return Err(Error::TooLarge);
    }
    Ok(())
}
pub fn capture(w: &Workspace, r: &Request) -> Result<Capture, Error> {
    records(w)?;
    let preview = preview(w, r)?;
    let bytes = encode(&preview)?;
    let attachment_hash = blake3::hash(&bytes).to_hex().to_string();
    if !w.attachments.contains_key(&attachment_hash)
        && (w.attachments.len() >= crate::documents::workspace::MAX_ATTACHMENTS
            || w.attachments.values().filter(|a| a.kind == KIND).count() >= MAX_RECORDS)
    {
        return Err(Error::TooLarge);
    }
    let mut workspace = w.clone();
    let attachment = WorkspaceAttachment {
        kind: KIND.into(),
        name: r.name.clone(),
        media_type: "application/json".into(),
        bytes,
    };
    if let Some(old) = workspace.attachments.get(&attachment_hash) {
        if old.kind != attachment.kind
            || old.name != attachment.name
            || old.media_type != attachment.media_type
            || old.bytes != attachment.bytes
        {
            return Err(bad("quantity hash has another attachment representation"));
        }
    } else {
        workspace
            .attachments
            .insert(attachment_hash.clone(), attachment);
    }
    workspace.validate()?;
    let result = Capture {
        workspace,
        context_revision: preview.context_revision.clone(),
        source_workspace_fingerprint: preview.source_workspace_fingerprint.clone(),
        attachment_hash,
        preview,
    };
    result.validate_binding(&result.context_revision)?;
    Ok(result)
}
impl Capture {
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        self.preview.validate_binding(expected)?;
        if self.context_revision != expected
            || self.source_workspace_fingerprint != self.preview.source_workspace_fingerprint
        {
            return Err(bad("quantity capture request changed"));
        }
        let saved = load(&self.workspace, &self.attachment_hash)?;
        if encode(&saved)? != encode(&self.preview)? {
            return Err(bad(
                "capture attachment differs from the exact quantity preview",
            ));
        }
        let mut original = self.workspace.clone();
        original.attachments.remove(&self.attachment_hash);
        if run_lineage::workspace_fingerprint(&original)? != self.source_workspace_fingerprint {
            return Err(bad(
                "capture changed documents, runs, candidate acceptance or unrelated attachments",
            ));
        }
        Ok(())
    }
}
pub fn load_record(bytes: &[u8]) -> Result<Preview, Error> {
    if bytes.len() > MAX_BYTES {
        return Err(Error::TooLarge);
    }
    let p: Preview = serde_json::from_slice(bytes).map_err(|e| bad(&e.to_string()))?;
    p.validate_binding(&p.context_revision)?;
    if encode(&p)?.as_slice() != bytes {
        return Err(bad(
            "quantity record bytes are not the current canonical representation",
        ));
    }
    Ok(p)
}
pub fn load(w: &Workspace, quantity_hash: &str) -> Result<Preview, Error> {
    records(w)?;
    if !digest(quantity_hash) {
        return Err(bad("invalid quantity attachment hash"));
    }
    let a = w
        .attachments
        .get(quantity_hash)
        .ok_or_else(|| bad("quantity document is not retained"))?;
    if a.kind != KIND
        || a.media_type != "application/json"
        || blake3::hash(&a.bytes).to_hex().as_str() != quantity_hash
    {
        return Err(bad("quantity attachment kind, media type or hash differs"));
    }
    let p = load_record(&a.bytes)?;
    if a.name != p.request.name {
        return Err(bad("quantity document name differs"));
    }
    p.validate_for_workspace(w)?;
    Ok(p)
}
fn list_hashes(w: &Workspace) -> Result<Vec<String>, Error> {
    records(w)?;
    let mut total = 0usize;
    w.attachments
        .iter()
        .filter(|(_, a)| a.kind == KIND)
        .map(|(h, a)| {
            if !digest(h)
                || a.media_type != "application/json"
                || a.bytes.len() > MAX_BYTES
                || blake3::hash(&a.bytes).to_hex().as_str() != h
            {
                return Err(bad("quantity listing attachment identity differs"));
            }
            total = total
                .checked_add(a.bytes.len())
                .filter(|n| *n <= MAX_REFERENCED_BYTES)
                .ok_or(Error::TooLarge)?;
            Ok(h.clone())
        })
        .collect()
}
fn list_context(fp: &str, hashes: &[String]) -> Result<String, Error> {
    hash(&(ENGINE_VERSION, "list", fp, hashes))
}
pub fn list_revision(w: &Workspace) -> Result<String, Error> {
    let hashes = list_hashes(w)?;
    list_context(&run_lineage::workspace_fingerprint(w)?, &hashes)
}
/// Discovery-only summaries. Open/load validates the complete record before use.
pub fn list(w: &Workspace) -> Result<Listed, Error> {
    let hashes = list_hashes(w)?;
    let fp = run_lineage::workspace_fingerprint(w)?;
    let mut scope = Scope::validated(w)?;
    let mut entries = vec![];
    for h in &hashes {
        let a = w
            .attachments
            .get(h)
            .ok_or_else(|| bad("quantity attachment absent"))?;
        let p = load_record(&a.bytes)?;
        if a.name != p.request.name {
            return Err(bad("quantity attachment name changed"));
        }
        let binding = scope.resolve(&p.request).and_then(|plan| {
            if plan.name == p.plan_name {
                Ok(())
            } else {
                Err(bad("retained plan identity differs"))
            }
        });
        let available = binding.is_ok();
        let findings = binding
            .err()
            .map(|e| vec![finding("binding", e)])
            .unwrap_or_default();
        entries.push(Summary {
            attachment_hash: h.clone(),
            name: p.request.name,
            plan_source: p.request.plan_source,
            plan_revision: p.request.plan_revision,
            target: p.request.target,
            inputs: p.request.inputs,
            complete: p.complete,
            reaction_count: p.request.reactions.len(),
            available,
            findings,
        });
    }
    let result = Listed {
        context_revision: list_context(&fp, &hashes)?,
        source_workspace_fingerprint: fp,
        entries,
    };
    result.validate_binding(&result.context_revision)?;
    Ok(result)
}
impl Listed {
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        if !digest(&self.source_workspace_fingerprint) || self.entries.len() > MAX_RECORDS {
            return Err(bad("quantity listing bounds/context changed"));
        }
        let mut hashes = vec![];
        let mut previous = None;
        for e in &self.entries {
            if !digest(&e.attachment_hash)
                || !digest(&e.plan_revision)
                || !text(&e.name, 256)
                || e.inputs.is_empty()
                || e.inputs.len() > 128
                || e.reaction_count == 0
                || e.reaction_count > MAX_REACTIONS
                || e.available != e.findings.is_empty()
                || e.findings.len() > 64
                || e.findings
                    .iter()
                    .any(|f| f.path.len() > 256 || f.message.len() > 4096)
                || previous.is_some_and(|p: &str| p >= e.attachment_hash.as_str())
            {
                return Err(bad("quantity listing summary malformed or unordered"));
            }
            crate::material::run::event_shape(&e.target)?;
            for i in &e.inputs {
                if !text(&i.input_id, 128) {
                    return Err(bad("quantity listing input ID malformed"));
                }
                crate::material::run::event_shape(&i.event)?;
            }
            if let PlanSource::SavedRun { run_id } = &e.plan_source {
                if !digest(run_id) {
                    return Err(bad("quantity listing run malformed"));
                }
            }
            previous = Some(e.attachment_hash.as_str());
            hashes.push(e.attachment_hash.clone());
        }
        if self.context_revision != expected
            || list_context(&self.source_workspace_fingerprint, &hashes)? != expected
        {
            return Err(bad(
                "quantity listing current workspace/attachment set changed",
            ));
        }
        encode(self)?;
        Ok(())
    }
}
fn open_context(fp: &str, h: &str) -> Result<String, Error> {
    hash(&(ENGINE_VERSION, "open", fp, h))
}
pub fn open_revision(w: &Workspace, h: &str) -> Result<String, Error> {
    load(w, h)?;
    open_context(&run_lineage::workspace_fingerprint(w)?, h)
}
pub fn open(w: &Workspace, h: &str) -> Result<Opened, Error> {
    let preview = load(w, h)?;
    let fp = run_lineage::workspace_fingerprint(w)?;
    let result = Opened {
        context_revision: open_context(&fp, h)?,
        source_workspace_fingerprint: fp,
        quantity_hash: h.into(),
        preview,
    };
    result.validate_binding(&result.context_revision)?;
    Ok(result)
}
impl Opened {
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        self.preview
            .validate_binding(&self.preview.context_revision)?;
        if !digest(&self.quantity_hash)
            || open_context(&self.source_workspace_fingerprint, &self.quantity_hash)? != expected
            || self.context_revision != expected
            || blake3::hash(&encode(&self.preview)?).to_hex().as_str() != self.quantity_hash
        {
            return Err(bad(
                "opened quantity scope, hash or current workspace changed",
            ));
        }
        Ok(())
    }
}
pub fn source_bindings(w: &Workspace, p: &Preview) -> Result<Dependencies, Error> {
    p.validate_for_workspace(w)?;
    let mut events = vec![p.request.target.clone()];
    events.extend(p.request.inputs.iter().map(|i| i.event.clone()));
    let (runs, candidates) = match &p.request.plan_source {
        PlanSource::Candidate => (vec![], vec![p.request.plan_revision.clone()]),
        PlanSource::SavedRun { run_id } => (vec![run_id.clone()], vec![]),
    };
    Ok(Dependencies {
        events,
        runs,
        candidates,
    })
}
fn export_context(
    h: &str,
    f: ExportFormat,
    fp: Option<&str>,
    plan: &str,
    target: &EventBinding,
    inputs: &[InputEventBinding],
) -> Result<String, Error> {
    if matches!(f, ExportFormat::ProtocolHtml | ExportFormat::ProtocolPdf) {
        let fp = fp
            .filter(|fp| digest(fp))
            .ok_or_else(|| bad("combined export workspace attestation missing"))?;
        hash(&(
            ENGINE_VERSION,
            units::ENGINE_VERSION,
            verification::ENGINE_VERSION,
            "export",
            h,
            f,
            fp,
            plan,
            target,
            inputs,
        ))
    } else {
        if fp.is_some() {
            return Err(bad("standalone quantity export has no workspace receipt"));
        }
        hash(&(ENGINE_VERSION, units::ENGINE_VERSION, "export", h, f))
    }
}
pub fn export_revision(w: &Workspace, h: &str, f: ExportFormat) -> Result<String, Error> {
    let p = load(w, h)?;
    let fp = if matches!(f, ExportFormat::ProtocolHtml | ExportFormat::ProtocolPdf) {
        Some(run_lineage::workspace_fingerprint(w)?)
    } else {
        None
    };
    export_context(
        h,
        f,
        fp.as_deref(),
        &p.request.plan_revision,
        &p.request.target,
        &p.request.inputs,
    )
}
// Borrowed serialization cap before native verification or large renderer allocations.
// Hashing the fresh report never clones/restores its products or scientific evidence.
struct CountHash {
    count: usize,
    limit: usize,
    hash: blake3::Hasher,
}
impl Write for CountHash {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        self.count = self
            .count
            .checked_add(bytes.len())
            .filter(|n| *n <= self.limit)
            .ok_or_else(|| io::Error::other("quantity protocol export bound"))?;
        self.hash.update(bytes);
        Ok(bytes.len())
    }
    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}
fn bounded_hash<T: Serialize>(value: &T, limit: usize) -> Result<String, Error> {
    let mut sink = CountHash {
        count: 0,
        limit,
        hash: blake3::Hasher::new(),
    };
    serde_json::to_writer(&mut sink, value).map_err(|_| Error::TooLarge)?;
    Ok(sink.hash.finalize().to_hex().to_string())
}
fn appendix(p: &Preview, h: &str) -> Result<verification::output::ProtocolAppendix, Error> {
    let mut tables = Vec::with_capacity(p.reactions.len());
    for (reaction, report) in p.request.reactions.iter().zip(&p.reactions) {
        let mut rows = Vec::with_capacity(reaction.components.len() + 1);
        for (component, resolved) in reaction.components.iter().zip(&report.components) {
            rows.push(vec![
                component.id.clone(),
                component.label.clone(),
                String::from_utf8(encode(&component.source)?).map_err(|_| Error::TooLarge)?,
                volume(resolved.volume),
            ]);
        }
        if let Some(fill) = &report.fill {
            rows.push(vec![
                fill.component_id.clone(),
                fill.label.clone(),
                "Explicit diluent fill".into(),
                volume(Some(fill.volume)),
            ]);
        }
        tables.push(verification::output::ProtocolTable {
            heading: format!(
                "Reaction {} / operation {} / final {}",
                reaction.id,
                reaction.step_index + 1,
                volume(reaction.final_volume)
            ),
            columns: vec![
                "Component ID".into(),
                "Declared label".into(),
                "Planned source".into(),
                "Applied volume".into(),
            ],
            rows,
        });
    }
    let appendix = verification::output::ProtocolAppendix {
        paragraphs: paragraph_lines(p, h)?,
        tables,
    };
    appendix.validate()?;
    Ok(appendix)
}
fn composed(
    w: &Workspace,
    p: &Preview,
    h: &str,
    format: ExportFormat,
) -> Result<(Vec<u8>, MolecularVerification, Vec<String>), Error> {
    let mut scope = Scope::new(w)?;
    let plan = scope.resolve(&p.request)?;
    let target = crate::planning::combinatorial::resolve_event(w, &p.request.target)?;
    // Scope resolution enforces exact native models/pins. This additionally bounds
    // borrowed source+plan bytes before fresh reconstruction and rendering.
    crate::output::pdf::preflight(plan, target)?;
    bounded_hash(&(plan, target, p), MAX_REFERENCED_BYTES)?;
    let appendix = appendix(p, h)?;
    let checked =
        verification::output::CheckedProtocol::new(plan, target).with_appendix(&appendix)?;
    let report = checked.report();
    let receipt = MolecularVerification {
        engine_version: report.engine_version.clone(),
        plan_revision: report.plan_revision.clone(),
        target_revision: report.target_revision.clone(),
        status: report.status,
        report_revision: bounded_hash(report, MAX_REFERENCED_BYTES)?,
        findings_total: report.findings.len(),
        exhaustive_alternative_search: report.exhaustive_alternative_search,
    };
    let (bytes, warnings) = match format {
        ExportFormat::ProtocolHtml => (
            verification::output::CheckedProtocol::html(&checked).into_bytes(),
            vec![],
        ),
        ExportFormat::ProtocolPdf => {
            let pdf = crate::output::pdf::export_protocol_pdf_checked(&checked)?;
            (pdf.bytes, pdf.warnings)
        }
        _ => return Err(bad("combined export format required")),
    };
    Ok((bytes, receipt, warnings))
}
fn volume(v: Option<units::Volume>) -> String {
    v.map(|v| {
        format!(
            "{}/{} nL",
            v.canonical().numerator(),
            v.canonical().denominator()
        )
    })
    .unwrap_or_else(|| "missing".into())
}
fn escape(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
        .replace('\'', "&#39;")
}
fn paragraph_lines(p: &Preview, h: &str) -> Result<Vec<String>, Error> {
    let mut lines = vec![
        format!(
            "Quantity document {}. Context {}. Plan {} / {}.",
            h, p.context_revision, p.request.plan_revision, p.plan_name
        ),
        format!(
            "Independent target: {}",
            String::from_utf8(encode(&p.request.target)?).map_err(|_| Error::TooLarge)?
        ),
        format!(
            "Exact original input events: {}",
            String::from_utf8(encode(&p.request.inputs)?).map_err(|_| Error::TooLarge)?
        ),
        format!(
            "Arithmetic completeness for declared reactions only: {}. Unlisted operations and omitted reagents remain unspecified. This is independent of molecular verification.",
            p.complete
        ),
    ];
    lines.extend(p.warnings.clone());
    for (reaction, report) in p.request.reactions.iter().zip(&p.reactions) {
        lines.push(format!(
            "Reaction {} / operation {}: final {}; complete {}",
            reaction.id,
            reaction.step_index + 1,
            volume(reaction.final_volume),
            report.complete
        ));
        for (component, resolved) in reaction.components.iter().zip(&report.components) {
            lines.push(format!(
                "{} / {}. Planned source {}. Applied stock/component volume {}.",
                component.id,
                component.label,
                String::from_utf8(encode(&component.source)?).map_err(|_| Error::TooLarge)?,
                volume(resolved.volume)
            ));
            if let Some(c) = &resolved.calculation {
                lines.push(format!("Explicit dilution calculation (achieved concentration and rounding retained): {}",String::from_utf8(encode(c)?).map_err(|_|Error::TooLarge)?));
            }
        }
        if let Some(m) = &report.mixture {
            lines.push(format!(
                "Specified {}; arithmetic remainder {}.",
                volume(Some(m.specified_volume)),
                volume(Some(m.remaining_diluent))
            ));
        }
        if let Some(fill) = &report.fill {
            lines.push(format!(
                "Explicit fill {} / {}: {}.",
                fill.component_id,
                fill.label,
                volume(Some(fill.volume))
            ));
        }
    }
    lines.extend(
        p.findings
            .iter()
            .map(|f| format!("Unresolved {}: {}", f.path, f.message)),
    );
    Ok(lines)
}
pub fn export(w: &Workspace, h: &str, f: ExportFormat) -> Result<Exported, Error> {
    let p = load(w, h)?;
    let source_workspace_fingerprint =
        if matches!(f, ExportFormat::ProtocolHtml | ExportFormat::ProtocolPdf) {
            Some(run_lineage::workspace_fingerprint(w)?)
        } else {
            None
        };
    let context_revision = export_context(
        h,
        f,
        source_workspace_fingerprint.as_deref(),
        &p.request.plan_revision,
        &p.request.target,
        &p.request.inputs,
    )?;
    let mut molecular = None;
    let mut extra = vec![];
    let (filename, media_type) = match f {
        ExportFormat::Csv => ("reaction-quantities.csv", "text/csv"),
        ExportFormat::Html => ("reaction-quantities.html", "text/html"),
        ExportFormat::Pdf => ("reaction-quantities.pdf", "application/pdf"),
        ExportFormat::ProtocolHtml => ("cloning-protocol-quantities.html", "text/html"),
        ExportFormat::ProtocolPdf => ("cloning-protocol-quantities.pdf", "application/pdf"),
    };
    let bytes = match f {
        ExportFormat::Csv => {
            let mut writer = csv::Writer::from_writer(Bytes(vec![]));
            writer
                .write_record([
                    "record",
                    "quantity_hash",
                    "quantity_name",
                    "quantity_context",
                    "plan_revision",
                    "target_event_json",
                    "reaction",
                    "step_index",
                    "component",
                    "label",
                    "planned_source_json",
                    "volume_nL",
                    "quantity_complete",
                    "detail_json",
                ])
                .map_err(|_| Error::TooLarge)?;
            let target =
                String::from_utf8(encode(&p.request.target)?).map_err(|_| Error::TooLarge)?;
            for (r, report) in p.request.reactions.iter().zip(&p.reactions) {
                for (c, v) in r.components.iter().zip(&report.components) {
                    writer
                        .write_record([
                            "component".into(),
                            h.into(),
                            p.request.name.clone(),
                            p.context_revision.clone(),
                            p.request.plan_revision.clone(),
                            target.clone(),
                            r.id.clone(),
                            r.step_index.to_string(),
                            c.id.clone(),
                            c.label.clone(),
                            String::from_utf8(encode(&c.source)?).map_err(|_| Error::TooLarge)?,
                            volume(v.volume),
                            report.complete.to_string(),
                            String::from_utf8(encode(v)?).map_err(|_| Error::TooLarge)?,
                        ])
                        .map_err(|_| Error::TooLarge)?;
                }
                writer
                    .write_record([
                        "reaction".into(),
                        h.into(),
                        p.request.name.clone(),
                        p.context_revision.clone(),
                        p.request.plan_revision.clone(),
                        target.clone(),
                        r.id.clone(),
                        r.step_index.to_string(),
                        String::new(),
                        String::new(),
                        String::new(),
                        volume(r.final_volume),
                        report.complete.to_string(),
                        String::from_utf8(encode(report)?).map_err(|_| Error::TooLarge)?,
                    ])
                    .map_err(|_| Error::TooLarge)?;
            }
            writer
                .write_record([
                    "scope".into(),
                    h.into(),
                    p.request.name.clone(),
                    p.context_revision.clone(),
                    p.request.plan_revision.clone(),
                    target,
                    String::new(),
                    String::new(),
                    String::new(),
                    String::new(),
                    String::new(),
                    String::new(),
                    p.complete.to_string(),
                    String::from_utf8(encode(&(&p.request.inputs, &p.findings, &p.warnings))?)
                        .map_err(|_| Error::TooLarge)?,
                ])
                .map_err(|_| Error::TooLarge)?;
            writer.into_inner().map_err(|_| Error::TooLarge)?.0
        }
        ExportFormat::Html => {
            let mut out = Bytes(vec![]);
            write!(out,"<!doctype html><html><head><meta charset=\"utf-8\"><title>{}</title><style>body{{font-family:system-ui;max-width:70em;margin:2em auto;padding:1em}}p{{white-space:pre-wrap;overflow-wrap:anywhere}}@media print{{body{{margin:0}}}}</style></head><body><h1>{}</h1>",escape(&p.request.name),escape(&p.request.name)).map_err(|_|Error::TooLarge)?;
            for line in paragraph_lines(&p, h)? {
                write!(out, "<p>{}</p>", escape(&line)).map_err(|_| Error::TooLarge)?;
            }
            write!(out, "</body></html>").map_err(|_| Error::TooLarge)?;
            out.0
        }
        ExportFormat::Pdf => {
            let blocks = paragraph_lines(&p, h)?
                .into_iter()
                .map(crate::output::pdf::ReviewBlock::Paragraph)
                .collect::<Vec<_>>();
            let (bytes, _, warnings) = crate::output::pdf::export_documentary_pdf(
                crate::output::pdf::DocumentaryKind::ReactionQuantities,
                &p.request.name,
                "Declared quantities; no chemical verification",
                h,
                &p.context_revision,
                &blocks,
            )?;
            extra = warnings;
            bytes
        }
        ExportFormat::ProtocolHtml | ExportFormat::ProtocolPdf => {
            let (bytes, receipt, warnings) = composed(w, &p, h, f)?;
            molecular = Some(receipt);
            extra = warnings;
            bytes
        }
    };
    if bytes.len() > MAX_EXPORT_BYTES
        || (!matches!(f, ExportFormat::ProtocolHtml | ExportFormat::ProtocolPdf)
            && bytes.len() > MAX_BYTES)
    {
        return Err(Error::TooLarge);
    }
    let content_hash = blake3::hash(&bytes).to_hex().to_string();
    let result = Exported {
        context_revision,
        quantity_hash: h.into(),
        source_workspace_fingerprint,
        plan_revision: p.request.plan_revision.clone(),
        target: p.request.target.clone(),
        inputs: p.request.inputs.clone(),
        arithmetic_complete: p.complete,
        molecular,
        format: f,
        filename: filename.into(),
        media_type: media_type.into(),
        bytes,
        content_hash,
        warnings: extra,
    };
    result.validate_for_preview(&p)?;
    Ok(result)
}
impl Exported {
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        let (name, mime) = match self.format {
            ExportFormat::Csv => ("reaction-quantities.csv", "text/csv"),
            ExportFormat::Html => ("reaction-quantities.html", "text/html"),
            ExportFormat::Pdf => ("reaction-quantities.pdf", "application/pdf"),
            ExportFormat::ProtocolHtml => ("cloning-protocol-quantities.html", "text/html"),
            ExportFormat::ProtocolPdf => ("cloning-protocol-quantities.pdf", "application/pdf"),
        };
        if !digest(&self.quantity_hash)
            || self.context_revision != expected
            || !digest(&self.plan_revision)
            || crate::material::run::event_shape(&self.target).is_err()
            || self.inputs.len() > 128
            || self.inputs.iter().any(|i| {
                !text(&i.input_id, 128) || crate::material::run::event_shape(&i.event).is_err()
            })
            || export_context(
                &self.quantity_hash,
                self.format,
                self.source_workspace_fingerprint.as_deref(),
                &self.plan_revision,
                &self.target,
                &self.inputs,
            )? != expected
            || self.bytes.len() > MAX_EXPORT_BYTES
            || (!matches!(
                self.format,
                ExportFormat::ProtocolHtml | ExportFormat::ProtocolPdf
            ) && self.bytes.len() > MAX_BYTES)
            || self.bytes.is_empty()
            || blake3::hash(&self.bytes).to_hex().as_str() != self.content_hash
            || self.filename != name
            || self.media_type != mime
            || self.warnings.len() > 128
            || self.warnings.iter().any(|w| w.len() > 4096)
        {
            return Err(bad(
                "quantity export source, content, format or context changed",
            ));
        }
        bounded_hash(self, MAX_EXPORT_ENVELOPE_BYTES)?;
        let combined = matches!(
            self.format,
            ExportFormat::ProtocolHtml | ExportFormat::ProtocolPdf
        );
        if combined != self.molecular.is_some() {
            return Err(bad("molecular verification receipt/format mismatch"));
        }
        if let Some(m) = &self.molecular {
            m.validate_binding(&self.plan_revision, &self.target.revision)?;
        }
        if matches!(self.format, ExportFormat::Pdf | ExportFormat::ProtocolPdf)
            && !self.bytes.starts_with(b"%PDF-")
        {
            return Err(bad("quantity PDF bytes malformed"));
        }
        Ok(())
    }
}

impl MolecularVerification {
    /// A documentary receipt of fresh native verification, not reusable material.
    pub fn validate_binding(
        &self,
        plan_revision: &str,
        target_revision: &str,
    ) -> Result<(), Error> {
        if self.engine_version != verification::ENGINE_VERSION
            || self.plan_revision != plan_revision
            || self.target_revision != target_revision
            || !digest(&self.plan_revision)
            || !digest(&self.target_revision)
            || !digest(&self.report_revision)
            || self.findings_total > 1_000_000
        {
            return Err(bad("molecular export receipt binding changed"));
        }
        Ok(())
    }
}
impl Exported {
    pub fn validate_for_preview(&self, p: &Preview) -> Result<(), Error> {
        self.validate_binding(&self.context_revision)?;
        p.validate_binding(&p.context_revision)?;
        if blake3::hash(&encode(p)?).to_hex().as_str() != self.quantity_hash
            || self.plan_revision != p.request.plan_revision
            || self.target != p.request.target
            || encode(&self.inputs)? != encode(&p.request.inputs)?
            || self.arithmetic_complete != p.complete
        {
            return Err(bad("quantity export differs from exact captured preview"));
        }
        Ok(())
    }
}
