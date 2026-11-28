//! Revision-bound reference gel predictions and diagnostic digest ranking.
use crate::platform::worker;
use collagene_core::{
    cloning::restriction::setups::Recipe,
    material::lineage::EventBinding,
    output::gel::{
        self, Calibration, CalibrationPoint, DiagnosticOptions, DiagnosticReport, GelReport,
        LaneInput, Preparation,
    },
    protocol::{
        request::Request,
        response::Response,
        schema::{AnalysisJob, AnalysisResult},
    },
    sequence::analysis::Primer,
};
use collagene_model::Molecule;
use leptos::prelude::*;
use std::sync::{
    Arc,
    atomic::{AtomicBool, Ordering},
};
use wasm_bindgen_futures::spawn_local;

/// One saved, typed digest setup queued by its workspace host.
#[derive(Clone)]
pub struct DigestLaneSeed {
    pub name: String,
    pub recipe: Recipe,
    pub epoch: u64,
}
#[derive(Clone, PartialEq, Eq)]
enum LaneSource {
    Open(usize),
    Historical { binding: EventBinding, epoch: u64 },
}
#[derive(Clone)]
struct LaneChoice {
    key: u32,
    source: LaneSource,
    label: String,
    setup_name: Option<String>,
    preparation: Preparation,
    mass_ng: Option<f64>,
}
#[derive(Clone, serde::Serialize)]
struct PinnedDigestLane {
    /// Zero-based lane position in the returned gel report.
    lane_index: usize,
    label: String,
    setup_name: Option<String>,
    mass_ng: Option<f64>,
    recipe: Recipe,
}
#[derive(Clone, serde::Serialize)]
struct GelEvidence {
    schema_version: u16,
    report: GelReport,
    /// Captured from the successful request, never from mutable current forms.
    pinned_digest_lanes: Vec<PinnedDigestLane>,
}
fn resolve_lane(
    lane: &LaneChoice,
    epoch: u64,
    get_query: Callback<usize, Option<Molecule>>,
    resolve_event: Callback<EventBinding, Result<Molecule, String>>,
) -> Result<Molecule, String> {
    match &lane.source {
        LaneSource::Open(index) => get_query
            .run(*index)
            .ok_or_else(|| "A lane source is no longer open".into()),
        LaneSource::Historical {
            binding,
            epoch: saved_epoch,
        } => {
            if *saved_epoch != epoch {
                return Err("A historical digest lane belongs to a previous workspace. Remove it and add its saved setup again.".into());
            }
            let molecule = resolve_event.run(binding.clone())?;
            if molecule.revision() != binding.revision {
                return Err("The retained digest source does not match its exact revision.".into());
            }
            Ok(molecule)
        }
    }
}
fn bounded_lane_label(mut label: String) -> String {
    if label.len() > 256 {
        let mut boundary = 256;
        while !label.is_char_boundary(boundary) {
            boundary -= 1;
        }
        label.truncate(boundary);
    }
    label
}
fn preparation_label(preparation: &Preparation) -> String {
    match preparation {
        Preparation::Uncut => "Uncut reference".into(),
        Preparation::Digest { enzymes } => {
            format!("Complete restriction digest · {}", enzymes.join(" + "))
        }
        Preparation::Pcr { forward, reverse } => format!(
            "Possible exact-PCR products · {} / {}",
            forward.name, reverse.name
        ),
        Preparation::CatalogDigest { recipe } => format!(
            "Historical digest setup · {} · pinned catalog {} · {} · assumed complete cleavage; methylation unassessed",
            recipe
                .selection
                .enzymes
                .iter()
                .map(|e| e.enzyme.name.as_str())
                .collect::<Vec<_>>()
                .join(" + "),
            recipe.selection.catalog_edition,
            match &recipe.origin {
                collagene_core::material::derivation::SourceOrigin::OriginalReference =>
                    "original-reference source",
                collagene_core::material::derivation::SourceOrigin::Saved(_) =>
                    "exact saved-material origin",
            },
        ),
    }
}
fn num<T: std::str::FromStr>(text: &str, label: &str) -> Result<T, String> {
    text.trim()
        .parse()
        .map_err(|_| format!("{label} must be a valid number"))
}
fn names(text: &str) -> Vec<String> {
    text.split(',')
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .map(str::to_owned)
        .collect()
}
fn calibration(text: &str, resolution: &str) -> Result<Calibration, String> {
    if text.len() > 8192 {
        return Err("Ladder exceeds its text limit".into());
    }
    let mut points = vec![];
    for line in text.lines().filter(|l| !l.trim().is_empty()) {
        let words = line
            .split(|c: char| c == ',' || c.is_ascii_whitespace())
            .filter(|s| !s.is_empty())
            .collect::<Vec<_>>();
        if words.len() != 2 {
            return Err("Each ladder line needs size in bp, then normalized migration".into());
        }
        points.push(CalibrationPoint {
            size_bp: num(words[0], "Ladder size")?,
            migration: num(words[1], "Ladder migration")?,
        });
    }
    let result = Calibration {
        label: "User-supplied illustrative/calibration points; no run conditions inferred".into(),
        points,
        resolution: num(resolution, "Migration resolution")?,
    };
    result.validate().map_err(|e| e.to_string())?;
    Ok(result)
}
fn short(s: &str) -> String {
    s.chars().take(12).collect()
}
#[component]
pub fn GelPanel(
    molecule: Memo<Molecule>,
    query_options: Memo<Vec<(usize, String)>>,
    get_query: Callback<usize, Option<Molecule>>,
    resolve_event: Callback<EventBinding, Result<Molecule, String>>,
    seeds: RwSignal<Vec<DigestLaneSeed>>,
    epoch: RwSignal<u64>,
) -> impl IntoView {
    let source = RwSignal::new(0usize);
    let prep = RwSignal::new("digest".to_string());
    let enzymes = RwSignal::new("EcoRI".to_string());
    let mass = RwSignal::new(String::new());
    let f = RwSignal::new(String::new());
    let r = RwSignal::new(String::new());
    let ft = RwSignal::new(String::new());
    let rt = RwSignal::new(String::new());
    let ladder = RwSignal::new("100, 0.9\n1000, 0.5\n10000, 0.1".to_string());
    let resolution = RwSignal::new("0.02".to_string());
    let lanes = RwSignal::new(Vec::<LaneChoice>::new());
    let next_key = RwSignal::new(0u32);
    let diagnostic_indices = RwSignal::new(Vec::<usize>::new());
    let diagnostic_enzymes =
        RwSignal::new("EcoRI,BamHI,HindIII,EcoRV,SmaI,NotI,XbaI,BsaI,PstI".to_string());
    let double = RwSignal::new(true);
    let report = RwSignal::new(None::<GelEvidence>);
    let diagnostics = RwSignal::new(None::<DiagnosticReport>);
    let status = RwSignal::new(String::new());
    let busy = RwSignal::new(false);
    let generation = RwSignal::new(0u32);
    let alive = Arc::new(AtomicBool::new(true));
    {
        let alive = alive.clone();
        on_cleanup(move || alive.store(false, Ordering::Relaxed));
    }
    let previous_epoch = StoredValue::new(epoch.get_untracked());
    Effect::new(move |_| {
        let current = epoch.get();
        if previous_epoch.get_value() == current {
            return;
        }
        previous_epoch.set_value(current);
        // Open-document indices and historical event lanes belong to the old
        // workspace alike. Keep calibration/preferences, never retarget sources.
        batch(|| {
            generation.update(|g| *g = g.wrapping_add(1));
            lanes.set(Vec::new());
            diagnostic_indices.set(Vec::new());
            report.set(None);
            diagnostics.set(None);
            source.set(0);
            f.set(String::new());
            r.set(String::new());
            ft.set(String::new());
            rt.set(String::new());
            mass.set(String::new());
            seeds.update(|queued| queued.retain(|seed| seed.epoch == current));
            busy.set(false);
        });
        status.set("Workspace changed; gel sources and predictions were cleared. Add lanes from this workspace.".into());
    });
    let snapshot = move || {
        let workspace_epoch = epoch.get();
        let active = molecule.get().revision();
        let choices = lanes.get();
        let candidates = diagnostic_indices.get();
        let revisions = choices
            .iter()
            .map(|l| {
                (
                    l.key,
                    l.source.clone(),
                    resolve_lane(l, workspace_epoch, get_query, resolve_event)
                        .map(|m| m.revision()),
                )
            })
            .collect::<Vec<_>>();
        let candidate_revisions = candidates
            .iter()
            .map(|i| (*i, get_query.run(*i).map(|m| m.revision())))
            .collect::<Vec<_>>();
        // Include editable calibration/options directly so a just-completed
        // asynchronous request cannot win before the invalidation effect runs.
        (
            workspace_epoch,
            active,
            revisions,
            candidate_revisions,
            ladder.get(),
            resolution.get(),
            diagnostic_enzymes.get(),
            double.get(),
        )
    };
    let binding = Memo::new(move |_| snapshot());
    Effect::new(move |_| {
        let _ = (
            binding.get(),
            lanes.get(),
            ladder.get(),
            resolution.get(),
            diagnostic_enzymes.get(),
            double.get(),
        );
        generation.update(|g| *g = g.wrapping_add(1));
        report.set(None);
        diagnostics.set(None);
        status.set(String::new());
        busy.set(false);
    });
    Effect::new(move |_| {
        let queued_count = seeds.with(|queued| queued.len());
        if queued_count == 0 {
            return;
        }
        if queued_count > 32
            || lanes
                .with_untracked(|current| current.len())
                .saturating_add(queued_count)
                > 32
        {
            seeds.set(Vec::new());
            status.set("Queued digest setups exceed the 32-lane limit. Existing lanes were preserved; remove lanes and add the setups again.".into());
            return;
        }
        let queued = seeds.get_untracked();
        let current_epoch = epoch.get_untracked();
        let mut additions = Vec::with_capacity(queued.len());
        let mut key = next_key.get_untracked();
        for seed in queued {
            let checked = (|| -> Result<LaneChoice, String> {
                if seed.epoch != current_epoch {
                    return Err(
                        "The queued digest setups belong to a previous workspace; add them again."
                            .into(),
                    );
                }
                if seed.name.trim().is_empty()
                    || seed.name.len() > 256
                    || seed.name.chars().any(char::is_control)
                {
                    return Err("The saved digest setup name is invalid.".into());
                }
                seed.recipe.validate().map_err(|e| e.to_string())?;
                let source = seed.recipe.source.clone();
                let original = resolve_event.run(source.clone())?;
                if original.revision() != source.revision {
                    return Err(
                        "The retained digest source revision differs from the saved setup.".into(),
                    );
                }
                Ok(LaneChoice {
                    key,
                    source: LaneSource::Historical {
                        binding: source,
                        epoch: current_epoch,
                    },
                    label: bounded_lane_label(format!(
                        "{} · {}",
                        seed.name,
                        original.name.chars().take(64).collect::<String>()
                    )),
                    setup_name: Some(seed.name.clone()),
                    preparation: Preparation::CatalogDigest {
                        recipe: seed.recipe,
                    },
                    mass_ng: None,
                })
            })();
            match checked {
                Ok(lane) => {
                    additions.push(lane);
                    key = key.wrapping_add(1);
                }
                Err(error) => {
                    seeds.set(Vec::new());
                    status.set(format!("No queued lanes were added: {error}"));
                    return;
                }
            }
        }
        let count = additions.len();
        batch(|| {
            seeds.set(Vec::new());
            lanes.update(|current| current.extend(additions));
            next_key.set(key);
        });
        status.set(format!("Added {count} historical digest setup lane(s). Predict gel lanes freshly inspects the exact retained source events."));
    });
    let run = {
        let alive = alive.clone();
        move |job: AnalysisJob, expected: Vec<String>| {
            let revision = molecule.get_untracked().revision();
            let bound = untrack(snapshot);
            let submitted_choices = lanes.get_untracked();
            let (pinned_digest_lanes, expected_calibration, expected_labels, expect_pinned) =
                match &job {
                    AnalysisJob::PredictGel {
                        lanes, calibration, ..
                    }
                    | AnalysisJob::PredictPinnedGel {
                        lanes, calibration, ..
                    } => {
                        let pinned = lanes
                            .iter()
                            .enumerate()
                            .filter_map(|(lane_index, lane)| match &lane.preparation {
                                Preparation::CatalogDigest { recipe } => Some(PinnedDigestLane {
                                    lane_index,
                                    label: lane.label.clone(),
                                    setup_name: submitted_choices
                                        .get(lane_index)
                                        .and_then(|choice| choice.setup_name.clone()),
                                    mass_ng: lane.mass_ng,
                                    recipe: recipe.clone(),
                                }),
                                _ => None,
                            })
                            .collect::<Vec<_>>();
                        let is_pinned = matches!(&job, AnalysisJob::PredictPinnedGel { .. });
                        (
                            pinned,
                            Some(calibration.clone()),
                            lanes
                                .iter()
                                .map(|lane| lane.label.clone())
                                .collect::<Vec<_>>(),
                            is_pinned,
                        )
                    }
                    _ => (Vec::new(), None, Vec::new(), false),
                };
            let token = generation.get_untracked().wrapping_add(1);
            generation.set(token);
            busy.set(true);
            status.set("Computing reference predictions in background…".into());
            let alive = alive.clone();
            spawn_local(async move {
                let result = worker::analyze(Request {
                    id: token,
                    revision: revision.clone(),
                    job,
                })
                .await;
                if !alive.load(Ordering::Relaxed)
                    || token != generation.get_untracked()
                    || untrack(snapshot) != bound
                {
                    return;
                }
                busy.set(false);
                match result {
                    Ok(Response {
                        id,
                        revision: returned,
                        result,
                    }) if id == token && returned == revision => match match result {
                        AnalysisResult::PinnedGel(result) if expect_pinned => {
                            AnalysisResult::Gel(result)
                        }
                        AnalysisResult::Gel(_) if expect_pinned => {
                            AnalysisResult::Error("Expected a pinned gel worker response".into())
                        }
                        AnalysisResult::PinnedGel(_) => {
                            AnalysisResult::Error("Unexpected pinned gel worker response".into())
                        }
                        result => result,
                    } {
                        AnalysisResult::Gel(result)
                            if result
                                .lanes
                                .iter()
                                .map(|l| l.source_revision.clone())
                                .collect::<Vec<_>>()
                                == expected
                                && expected_calibration.as_ref() == Some(&result.calibration)
                                && result
                                    .lanes
                                    .iter()
                                    .map(|lane| lane.label.clone())
                                    .collect::<Vec<_>>()
                                    == expected_labels =>
                        {
                            report.set(Some(GelEvidence {
                                schema_version: 2,
                                report: result,
                                pinned_digest_lanes,
                            }));
                            status.set("Gel prediction complete under the recorded calibration and cleavage assumptions.".into());
                        }
                        AnalysisResult::DiagnosticDigest(result)
                            if result.source_revisions == expected =>
                        {
                            diagnostics.set(Some(result));
                            status.set("Diagnostic digest ranking complete; inspect unknown pairs and search completeness.".into());
                        }
                        AnalysisResult::Error(error) => status.set(error),
                        _ => status.set("Worker source bindings did not match the request".into()),
                    },
                    Err(error) => status.set(error),
                    _ => status.set("Unexpected gel worker response".into()),
                }
            });
        }
    };
    let add = move |_| {
        let Some(m) = get_query.run(source.get_untracked()) else {
            status.set("Choose an open source document".into());
            return;
        };
        if lanes.with_untracked(|l| l.len() >= 32) {
            status.set("At most 32 gel lanes".into());
            return;
        }
        let mass_ng = if mass.get_untracked().trim().is_empty() {
            None
        } else {
            match num::<f64>(&mass.get_untracked(), "Loaded mass") {
                Ok(v) if v.is_finite() && v > 0. && v <= 1e9 => Some(v),
                _ => {
                    status.set(
                        "Loaded mass must be positive and finite, or blank for unknown".into(),
                    );
                    return;
                }
            }
        };
        let preparation = match prep.get_untracked().as_str() {
            "uncut" => Preparation::Uncut,
            "pcr" => match Primer::new("F", &f.get_untracked(), &ft.get_untracked()).and_then(|f| {
                Primer::new("R", &r.get_untracked(), &rt.get_untracked()).map(|r| (f, r))
            }) {
                Ok((forward, reverse)) => Preparation::Pcr { forward, reverse },
                Err(error) => {
                    status.set(error.to_string());
                    return;
                }
            },
            _ => Preparation::Digest {
                enzymes: names(&enzymes.get_untracked()),
            },
        };
        let key = next_key.get_untracked();
        next_key.set(key.wrapping_add(1));
        let label = m.name.chars().take(48).collect::<String>();
        lanes.update(|lanes| {
            lanes.push(LaneChoice {
                key,
                source: LaneSource::Open(source.get_untracked()),
                label,
                setup_name: None,
                preparation,
                mass_ng,
            })
        });
    };
    let predict = {
        let run = run.clone();
        move |_| {
            let calibration =
                match calibration(&ladder.get_untracked(), &resolution.get_untracked()) {
                    Ok(c) => c,
                    Err(e) => {
                        status.set(e);
                        return;
                    }
                };
            let inputs = lanes
                .get_untracked()
                .into_iter()
                .map(|l| {
                    resolve_lane(&l, epoch.get_untracked(), get_query, resolve_event).map(|m| {
                        LaneInput {
                            label: l.label,
                            molecule: m,
                            preparation: l.preparation,
                            mass_ng: l.mass_ng,
                        }
                    })
                })
                .collect::<Result<Vec<_>, _>>();
            match inputs {
                Ok(lanes) => {
                    let revisions = lanes.iter().map(|l| l.molecule.revision()).collect();
                    let pinned = lanes
                        .iter()
                        .any(|lane| matches!(&lane.preparation, Preparation::CatalogDigest { .. }));
                    let job = if pinned {
                        AnalysisJob::PredictPinnedGel {
                            molecule: molecule.get_untracked(),
                            lanes,
                            calibration,
                        }
                    } else {
                        AnalysisJob::PredictGel {
                            molecule: molecule.get_untracked(),
                            lanes,
                            calibration,
                        }
                    };
                    run(job, revisions);
                }
                Err(error) => status.set(error),
            }
        }
    };
    let diagnose = move |_| {
        let calibration = match calibration(&ladder.get_untracked(), &resolution.get_untracked()) {
            Ok(c) => c,
            Err(e) => {
                status.set(e);
                return;
            }
        };
        let candidates = diagnostic_indices
            .get_untracked()
            .into_iter()
            .map(|i| {
                get_query
                    .run(i)
                    .ok_or_else(|| "A diagnostic source is no longer open".to_string())
            })
            .collect::<Result<Vec<_>, _>>();
        match candidates {
            Ok(candidates) => {
                let revisions = candidates.iter().map(Molecule::revision).collect();
                run(
                    AnalysisJob::DiagnosticDigests {
                        molecule: molecule.get_untracked(),
                        candidates,
                        options: DiagnosticOptions {
                            enzymes: names(&diagnostic_enzymes.get_untracked()),
                            include_double_digests: double.get_untracked(),
                            calibration,
                            ..Default::default()
                        },
                    },
                    revisions,
                );
            }
            Err(error) => status.set(error),
        }
    };
    view! {
        <style>".gel-panel textarea{font:11px/1.6 monospace}.gel-panel .gel-lane{padding:8px 0;border-bottom:1px solid #e1e7e0;font-size:12px;overflow-wrap:anywhere}.gel-panel .gel-check{display:flex;align-items:center;gap:8px}.gel-panel .gel-check input{width:auto}.gel-panel .gel-source-picker{max-height:180px;overflow:auto}.gel-panel .gel-svg{overflow:auto;max-height:440px;border:1px solid #dce5de}.gel-panel .gel-table-wrap{overflow:auto}.gel-panel table{border-collapse:collapse;font-size:11px;white-space:nowrap}.gel-panel th,.gel-panel td{padding:5px;border-bottom:1px solid #e1e7e0;text-align:left}.gel-panel .gel-candidate{padding:10px 0;border-bottom:1px solid #e1e7e0;font-size:12px;line-height:1.6}.gel-panel ul{padding-left:18px}.gel-status{font-size:12px;line-height:1.6;overflow-wrap:anywhere}"</style>
        <div class="inspector-content gel-panel"><h2>"Predicted gel"</h2><p class="helper">"Reference predictions only. Circular conformation and out-of-ladder fragments remain unplaced. No observed-gel comparison."</p>
            <details><summary>"Ladder and resolution"</summary><label>"Size bp, migration fraction (one point per line)"<textarea aria-label="Gel ladder calibration" rows="4" prop:value=move||ladder.get() on:input=move|e|ladder.set(event_target_value(&e))/></label><label>"Migration resolution (0–1 scale)"<input aria-label="Gel migration resolution" type="number" step="0.001" prop:value=move||resolution.get() on:input=move|e|resolution.set(event_target_value(&e))/></label><p class="helper">"Default ladder is illustrative, not an empirical gel calibration. Positions interpolate versus log size only between supplied points."</p></details>
            <h3>"Add a source lane"</h3><label>"Source document"<select aria-label="Gel source document" prop:value=move||source.get().to_string() on:change=move|e|source.set(event_target_value(&e).parse().unwrap_or(0))>{move||query_options.get().into_iter().map(|(i,name)|view!{<option value=i.to_string()>{name}</option>}).collect_view()}</select></label>
            <label>"Preparation"<select aria-label="Gel preparation" prop:value=move||prep.get() on:change=move|e|prep.set(event_target_value(&e))><option value="digest">"Complete restriction digest"</option><option value="uncut">"Uncut reference"</option><option value="pcr">"Possible exact-PCR products"</option></select></label>
            <Show when=move||prep.get()=="digest"><label>"Enzymes (comma separated)"<input aria-label="Gel digest enzymes" prop:value=move||enzymes.get() on:input=move|e|enzymes.set(event_target_value(&e))/></label></Show>
            <Show when=move||prep.get()=="pcr">{[("Forward annealing","Gel forward primer",f),("Forward tail","Gel forward tail",ft),("Reverse annealing","Gel reverse primer",r),("Reverse tail","Gel reverse tail",rt)].into_iter().map(|(name,aria,value)|view!{<label>{name}<input aria-label=aria prop:value=move||value.get() on:input=move|e|value.set(event_target_value(&e))/></label>}).collect_view()}</Show>
            <label>"Loaded DNA mass (ng; blank means unknown)"<input aria-label="Gel loaded mass" type="number" prop:value=move||mass.get() on:input=move|e|mass.set(event_target_value(&e))/></label><button class="button secondary full" on:click=add>"Add gel lane"</button>
            <div class="gel-lanes"><For each=move||lanes.get() key=|l|l.key children=move|lane|{let key=lane.key;view!{<div class="gel-lane"><strong>{lane.label}</strong>{match &lane.source {LaneSource::Historical{binding,..}=>Some(view!{<p class="helper">{format!("Exact retained event {} · revision {}…",binding.event_id,short(&binding.revision))}</p>}),LaneSource::Open(_)=>None}}<p>{format!("{} · mass {}",preparation_label(&lane.preparation),lane.mass_ng.map(|m|format!("{m} ng")).unwrap_or_else(||"unknown".into()))}</p><button class="text-tool" on:click=move|_|lanes.update(|l|l.retain(|l|l.key!=key))>"Remove lane"</button></div>}}/></div>
            <button class="button primary full" disabled=move||busy.get()||lanes.with(|l|l.is_empty()) on:click=predict>"Predict gel lanes"</button><p class="gel-status" role="status">{move||status.get()}</p>
            {move||report.get().map(|evidence|{let downloadable=evidence.clone();let result=evidence.report;let svg=gel::render_svg(&result);let for_svg=result.clone();view!{<section class="gel-result"><h3>"Predicted lanes"</h3>{match svg{Ok(svg)=>Some(view!{<div class="gel-svg" inner_html=svg/>}),Err(_)=>None}}
                <button class="button secondary full" on:click=move|_|match gel::render_svg(&for_svg){Ok(svg)=>crate::io::download::download_with_completion("predicted-gel.svg","image/svg+xml",svg.as_bytes(),move|result|{if status.try_get_untracked().is_none(){return;}if let Err(e)=result {status.set(e)}}),Err(e)=>status.set(e.to_string())}>"Download gel SVG"</button><button class="button secondary full" on:click=move|_|match if downloadable.pinned_digest_lanes.is_empty(){serde_json::to_vec_pretty(&downloadable.report)}else{serde_json::to_vec_pretty(&downloadable)}{Ok(bytes)=>crate::io::download::download_with_completion("predicted-gel.json","application/json",&bytes,move|result|{if status.try_get_untracked().is_none(){return;}if let Err(e)=result {status.set(e)}}),Err(e)=>status.set(e.to_string())}>"Download gel report"</button>
                {result.lanes.into_iter().map(|lane|view!{<details class="gel-lane-report"><summary>{format!("{} · {} bands · {} unplaced",lane.label,lane.bands.len(),lane.unpositioned_fragments.len())}</summary><p>{format!("{} · source {}…",lane.preparation,short(&lane.source_revision))}</p>{lane.error.map(|e|view!{<p role="alert">{e}</p>})}<ul>{lane.warnings.into_iter().map(|w|view!{<li>{w}</li>}).collect_view()}</ul><div class="gel-table-wrap"><table><thead><tr><th>"Fragment"</th><th>"Mean bp"</th><th>"Strands nt"</th><th>"Migration"</th><th>"Mass ng"</th></tr></thead><tbody>{lane.fragments.iter().take(200).enumerate().map(|(i,f)|view!{<tr><td>{i+1}</td><td>{format!("{:.2}",f.size_bp)}</td><td>{format!("{} / {}",f.top_nt,f.bottom_nt)}</td><td>{f.migration.map(|p|format!("{p:.4}")).unwrap_or_else(||format!("Unplaced {:?}",f.state))}</td><td>{f.amount_ng.map(|m|format!("{m:.3}")).unwrap_or_else(||"unknown".into())}</td></tr>}).collect_view()}</tbody></table></div><p class="helper">{format!("Showing {} of {} fragments; JSON retains the full returned lane.",lane.fragments.len().min(200),lane.fragments.len())}</p></details>}).collect_view()}
                <details><summary>"Prediction assumptions"</summary><ul>{result.assumptions.into_iter().map(|a|view!{<li>{a}</li>}).collect_view()}</ul></details></section>}})}
            <h3>"Find diagnostic digests"</h3><p class="helper">"Choose 2–16 candidate documents. Unknown circular/unsupported/out-of-ladder profiles never count as identification evidence."</p>
            <div class="gel-source-picker">{move||query_options.get().into_iter().map(|(i,name)|view!{<label class="gel-check"><input type="checkbox" aria-label=format!("Diagnostic candidate {i}") prop:checked=move||diagnostic_indices.with(|v|v.contains(&i)) on:change=move|e|{let checked=event_target::<web_sys::HtmlInputElement>(&e).checked();diagnostic_indices.update(|v|{if checked&&!v.contains(&i){v.push(i);v.sort();}else if !checked{v.retain(|n|*n!=i);}});}/>{name}</label>}).collect_view()}</div>
            <label>"Candidate enzymes"<input aria-label="Diagnostic digest enzymes" prop:value=move||diagnostic_enzymes.get() on:input=move|e|diagnostic_enzymes.set(event_target_value(&e))/></label><label class="gel-check"><input aria-label="Diagnostic double digests" type="checkbox" prop:checked=move||double.get() on:change=move|e|double.set(event_target::<web_sys::HtmlInputElement>(&e).checked())/>"Include double digests"</label>
            <button class="button primary full" disabled=move||busy.get()||diagnostic_indices.with(|v|v.len()<2||v.len()>16) on:click=diagnose>"Rank diagnostic digests"</button>
            {move||diagnostics.get().map(|result|{let downloaded=result.clone();let source_indices=diagnostic_indices.get_untracked();let source_labels=result.source_revisions.iter().enumerate().map(|(i,revision)|format!("{}: {} · {}…",i,source_indices.get(i).and_then(|i|get_query.run(*i)).map(|m|m.name).unwrap_or_else(||"missing source".into()),short(revision))).collect::<Vec<_>>();view!{<section class="gel-diagnostics"><h3>{format!("{} retained digest candidates",result.candidates.len())}</h3><ol class="gel-diagnostic-bindings" start="0">{source_labels.into_iter().map(|label|view!{<li>{label}</li>}).collect_view()}</ol><p>{format!("{} enzyme sets evaluated · exhaustive requested search: {}",result.evaluated_digests,result.exhaustive_search)}</p><ul>{result.warnings.into_iter().map(|w|view!{<li>{w}</li>}).collect_view()}</ul>{result.candidates.into_iter().map(|candidate|{let names=candidate.enzymes.clone();let indices=source_indices.clone();view!{<article class="gel-candidate"><strong>{candidate.enzymes.join(" + ")}</strong><p>{format!("Distinguished {} pairs · indistinguishable {} · unknown {}",candidate.distinguished_pairs.len(),candidate.indistinguishable_pairs.len(),candidate.unknown_pairs.len())}</p><p>{format!("Candidate-index pairs (zero-based): {:?}",candidate.distinguished_pairs)}</p><button class="button secondary" on:click=move|_|{let mut choices=vec![];for index in &indices{if let Some(m)=get_query.run(*index){let key=next_key.get_untracked();next_key.set(key.wrapping_add(1));choices.push(LaneChoice{key,source:LaneSource::Open(*index),label:m.name.chars().take(48).collect(),setup_name:None,preparation:Preparation::Digest{enzymes:names.clone()},mass_ng:None});}}lanes.set(choices);}>"Use digest as gel lanes"</button></article>}}).collect_view()}<button class="button secondary full" on:click=move|_|match serde_json::to_vec_pretty(&downloaded){Ok(bytes)=>crate::io::download::download_with_completion("diagnostic-digests.json","application/json",&bytes,move|result|{if status.try_get_untracked().is_none(){return;}if let Err(e)=result {status.set(e)}}),Err(e)=>status.set(e.to_string())}>"Download digest ranking"</button></section>}})}
        </div>
    }
}
