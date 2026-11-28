//! Explicit source-generation declarations for checked restriction planning.
use collagene_core::chemistry::{
    methylation::MethylationState, tracked_duplex::GenerationAssumptions,
};
use collagene_model::Molecule;
use leptos::prelude::*;

#[derive(Clone, Copy)]
pub struct AssumptionSignals {
    pub host: RwSignal<Option<collagene_core::chemistry::source_host::SourceHost>>,
    pub dam: RwSignal<MethylationState>,
    pub dcm: RwSignal<MethylationState>,
    pub cpg: RwSignal<MethylationState>,
}
impl AssumptionSignals {
    pub fn new() -> Self {
        Self {
            host: RwSignal::new(None),
            dam: RwSignal::new(MethylationState::Unknown),
            dcm: RwSignal::new(MethylationState::Unknown),
            cpg: RwSignal::new(MethylationState::Unknown),
        }
    }
    pub fn track(self) {
        let _ = (
            self.host.get(),
            self.dam.get(),
            self.dcm.get(),
            self.cpg.get(),
        );
    }
    pub fn current(self) -> GenerationAssumptions {
        GenerationAssumptions {
            host: self.host.get_untracked(),
            dam: self.dam.get_untracked(),
            dcm: self.dcm.get_untracked(),
            cpg: self.cpg.get_untracked(),
        }
    }
}
pub fn state_label(state: MethylationState) -> &'static str {
    match state {
        MethylationState::Unknown => "unknown",
        MethylationState::Unmodified => "assume unmodified",
        MethylationState::FullyModified => "assume fully modified",
        MethylationState::HemiOrMixed => "hemi / mixed (unsupported)",
    }
}
fn state_key(state: MethylationState) -> &'static str {
    match state {
        MethylationState::Unknown => "unknown",
        MethylationState::Unmodified => "unmodified",
        MethylationState::FullyModified => "fully",
        MethylationState::HemiOrMixed => "mixed",
    }
}
fn parse_state(value: &str) -> MethylationState {
    match value {
        "unmodified" => MethylationState::Unmodified,
        "fully" => MethylationState::FullyModified,
        "mixed" => MethylationState::HemiOrMixed,
        _ => MethylationState::Unknown,
    }
}
#[component]
pub fn SourceAssumptions(side: &'static str, signals: AssumptionSignals) -> impl IntoView {
    view! {<fieldset class="checked-source-assumptions"><legend>{format!("{side} source generation")}</legend><HostChoice selected=Signal::derive(move || signals.host.get()) on_change=Callback::new(move |host| signals.host.set(host))/>{[("Dam",signals.dam),("Dcm",signals.dcm),("CpG",signals.cpg)].into_iter().map(|(category,value)|view!{<label>{format!("{category} source state")}<select aria-label=format!("{side} checked {category} assumption") prop:value=move||state_key(value.get()) on:change=move|e|value.set(parse_state(&event_target_value(&e)))>{[MethylationState::Unknown,MethylationState::Unmodified,MethylationState::FullyModified,MethylationState::HemiOrMixed].into_iter().map(|state|view!{<option value=state_key(state)>{state_label(state)}</option>}).collect_view()}</select></label>}).collect_view()}</fieldset>}
}
/// A producer summary is display evidence only; core verification reconstructs the
/// private positional state from declaration and operation history every time.
pub fn modification_summary(molecule: &Molecule) -> Option<String> {
    let value = molecule.metadata.get("tracked_modifications")?;
    let summary: collagene_core::chemistry::tracked_duplex::ModificationSummary =
        serde_json::from_str(value).ok()?;
    Some(format!(
        "Declared-source positional model · top {} nt / bottom {} nt · {} modified, {} unknown eligible and {} unmodified chemical calls · {} source declarations. Newly synthesized bases are unmodified under this closed three-system model; inherited positions are propagated without rescanning descendant motifs. State {}. This producer summary is not reusable state or evidence of a physical reagent.",
        summary.top_length,
        summary.bottom_length,
        summary.known_modified,
        summary.unknown_eligible,
        summary.known_unmodified,
        summary.generation.len(),
        summary.state_revision
    ))
}

pub fn prepare_report_export(
    report: &collagene_core::verification::evidence::VerificationReport,
) -> Result<crate::io::download::ExportArtifact, String> {
    struct Json(Vec<u8>);
    impl std::io::Write for Json {
        fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
            if bytes.len() > 32 * 1024 * 1024 - self.0.len() {
                return Err(std::io::Error::other(
                    "Reconstruction JSON exceeds the 32 MiB browser output bound",
                ));
            }
            self.0.extend_from_slice(bytes);
            Ok(bytes.len())
        }
        fn flush(&mut self) -> std::io::Result<()> {
            Ok(())
        }
    }
    let mut bytes = Json(vec![]);
    serde_json::to_writer_pretty(&mut bytes, report).map_err(|e| e.to_string())?;
    crate::io::download::ExportArtifact::prepare(
        "cloning-reconstruction.json",
        "application/json",
        &bytes.0,
    )
}

/// A documentary choice never changes the independent chemical-state controls.
#[component]
pub(crate) fn HostChoice(
    selected: Signal<Option<collagene_core::chemistry::source_host::SourceHost>>,
    on_change: Callback<Option<collagene_core::chemistry::source_host::SourceHost>>,
    #[prop(default = Signal::derive(|| false))] disabled: Signal<bool>,
) -> impl IntoView {
    use collagene_core::cloning::restriction::catalog::catalog::Catalog;
    let choices = StoredValue::new(Catalog::builtin().map(|catalog| catalog.hosts().to_vec()));
    let error = RwSignal::new(String::new());
    view! {
        <details class="source-host-choice">
            <summary>"Source host (documentary)"</summary>
            <label>"Host strain / supplier product"
                <select aria-label="Source host preset" disabled=move || disabled.get()
                    prop:value=move || selected.get().map(|host| host.preset.id).unwrap_or_default()
                    on:change=move |event| {
                        let id = event_target_value(&event);
                        if id.is_empty() {
                            on_change.run(None);
                            error.set(String::new());
                            return;
                        }
                        let result = collagene_core::chemistry::source_host::SourceHost::from_current(&id)
                            .map_err(|error| error.to_string());
                        match result {
                            Ok(host) => { on_change.run(Some(host)); error.set(String::new()); },
                            Err(message) => error.set(message),
                        }
                    }>
                    <option value="">"No host declaration"</option>
                    {choices.with_value(|choices| choices.as_ref().ok().into_iter().flatten().map(|host| {
                        view! {<option value=host.id.clone()>{format!("{} · {} {}", host.name, host.vendor, host.product_code)}</option>}
                    }).collect_view())}
                </select>
            </label>
            <p class="helper">"Host expectations are documentary. Choose Dam, Dcm and CpG chemical assumptions separately; strain identity does not establish measured methylation or a complete modification system."</p>
            {move || selected.get().map(|host| view! {
                <div>
                    <p class="helper">{format!("{} · {} · Dam {:?}, Dcm {:?}, CpG {:?}. {}", host.preset.species, host.preset.host_scope, host.preset.dam, host.preset.dcm, host.preset.cpg, host.preset.qualifiers.join(" "))}</p>
                    <p class="helper">{format!("Catalogue {} · {}", host.catalog_edition, host.catalog_content_hash)}</p>
                    <details><summary>"Host evidence sources"</summary>
                        {host.sources.iter().map(|source| view! {
                            <p class="helper"><a href=source.url.clone() target="_blank" rel="noopener noreferrer">{source.publisher.clone()}</a>
                                {format!(" · {} · retrieved {} · facts {}", source.edition, source.retrieved_on, source.facts_blake3)}</p>
                        }).collect_view()}
                    </details>
                </div>
            })}
            <p class="helper" role="status">{move || error.get()}</p>
        </details>
    }
}
