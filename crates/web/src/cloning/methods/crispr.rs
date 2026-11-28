//! Revision-bound guide design, selected-record search and literal donor review.
mod search_view;
use crate::io::download::ExportArtifact;
use crate::platform::worker::CachedWorker;
use crate::reads::contigs::ReadBinding;
use crate::reads::contigs::ReadNavigation;
use collagene_core::{
    protocol::{
        request::Request,
        schema::{AnalysisJob, AnalysisResult},
    },
    sequence::{
        crispr::{
            self, BulgePolicy, GuideCandidate, GuideOptions, GuideQuery, GuideReport,
            IntendedTarget, NucleaseProfile, OffTargetOptions, OffTargetReport, PamPolicy, PamSide,
            ReferenceInput,
        },
        template::{
            self as edit_template, EditTemplateCheck, EditTemplateOptions, EditTemplateReport,
            GuideBinding,
        },
    },
};
use collagene_model::{Alphabet, Molecule, Topology};
use leptos::prelude::*;
use serde::Serialize;
use std::{
    collections::BTreeMap,
    rc::Rc,
    sync::{
        Arc,
        atomic::{AtomicBool, Ordering},
    },
};
use wasm_bindgen_futures::spawn_local;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DocumentInfo {
    pub binding: ReadBinding,
    pub name: String,
    pub length: usize,
    pub alphabet: Alphabet,
    pub topology: Topology,
}
#[derive(Clone, Debug, PartialEq, Eq)]
struct ReferenceBinding {
    document: DocumentInfo,
    build: Option<String>,
    accession: Option<String>,
}
#[derive(Clone)]
struct BoundGuides {
    source: DocumentInfo,
    report: Arc<GuideReport>,
    epoch: u64,
}
#[derive(Clone)]
struct BoundSearch {
    source: DocumentInfo,
    guide: GuideCandidate,
    guide_options: GuideOptions,
    references: Vec<ReferenceBinding>,
    report: Arc<OffTargetReport>,
    epoch: u64,
}
#[derive(Clone)]
pub struct DonorRequest {
    pub source: DocumentInfo,
    pub target: DocumentInfo,
    pub report: Arc<EditTemplateReport>,
    pub checked: EditTemplateCheck,
    pub edited: bool,
    pub epoch: u64,
}
#[derive(Clone, Copy)]
struct JobSlot {
    token: RwSignal<u32>,
    busy: RwSignal<bool>,
    worker: StoredValue<Option<Rc<CachedWorker>>, LocalStorage>,
}
impl JobSlot {
    fn new() -> Self {
        Self {
            token: RwSignal::new(0),
            busy: RwSignal::new(false),
            worker: StoredValue::new_local(None),
        }
    }
    fn stop(self) {
        self.token.update(|t| *t = t.wrapping_add(1));
        if let Some(w) = self.worker.get_value() {
            w.close();
        }
        self.worker.set_value(None);
        self.busy.set(false);
    }
    fn start(self) -> Result<(Rc<CachedWorker>, u32), String> {
        self.stop();
        let w = Rc::new(CachedWorker::new()?);
        self.worker.set_value(Some(w.clone()));
        self.busy.set(true);
        Ok((w, self.token.get_untracked()))
    }
    fn current(self, token: u32) -> bool {
        self.token.get_untracked() == token
    }
}
fn whole<T: std::str::FromStr>(value: &str, label: &str) -> Result<T, String> {
    value
        .trim()
        .parse()
        .map_err(|_| format!("Enter a whole number for {label}."))
}
fn percent(value: &str, label: &str) -> Result<u16, String> {
    let n = value
        .trim()
        .parse::<f64>()
        .map_err(|_| format!("Enter a percentage for {label}."))?;
    if !n.is_finite()
        || !(0.0..=100.0).contains(&n)
        || (n * 100.0 - (n * 100.0).round()).abs() > 0.00001
    {
        return Err(format!(
            "{label} must be 0–100%, with at most two decimal places."
        ));
    }
    Ok((n * 100.0).round() as u16)
}
fn optional(value: &str) -> Option<String> {
    let s = value.trim();
    (!s.is_empty()).then(|| s.into())
}
fn positions(values: &[usize]) -> String {
    let mut s = values
        .iter()
        .take(40)
        .map(|v| (v + 1).to_string())
        .collect::<Vec<_>>()
        .join(", ");
    if values.len() > 40 {
        s.push_str(&format!(" … ({} source positions)", values.len()));
    }
    s
}
fn boundary(value: Option<usize>) -> String {
    value
        .map(|v| format!("{v} bases from source origin"))
        .unwrap_or_else(|| "outside available source".into())
}
fn preview(value: &str) -> String {
    if value.len() > 1000 {
        format!(
            "{}… ({} nt; full sequence in export)",
            &value[..1000],
            value.len()
        )
    } else {
        value.into()
    }
}
fn response(
    result: Result<collagene_core::protocol::response::Response, String>,
    id: u32,
    revision: &str,
) -> Result<AnalysisResult, String> {
    let value = result?;
    if value.id != id || value.revision != revision {
        return Err("The worker response does not match this source-bound request.".into());
    }
    match value.result {
        AnalysisResult::Error(e) => Err(e),
        v => Ok(v),
    }
}
fn document(
    get: Callback<String, Option<Molecule>>,
    info: &DocumentInfo,
) -> Result<Molecule, String> {
    get.run(info.binding.document_id.clone())
        .ok_or_else(|| "A bound document is no longer open.".into())
}
struct Json(Vec<u8>);
impl std::io::Write for Json {
    fn write(&mut self, b: &[u8]) -> std::io::Result<usize> {
        if b.len() > 32 * 1024 * 1024 - self.0.len() {
            return Err(std::io::Error::other(
                "CRISPR JSON exceeds the 32 MiB browser output bound",
            ));
        }
        self.0.extend_from_slice(b);
        Ok(b.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}
fn json_bytes(value: &impl Serialize) -> Result<Vec<u8>, String> {
    let mut bytes = Json(vec![]);
    serde_json::to_writer_pretty(&mut bytes, value).map_err(|e| e.to_string())?;
    Ok(bytes.0)
}
#[derive(Serialize)]
struct Binding<'a> {
    document_id: &'a str,
    event_id: &'a str,
    revision: &'a str,
}
fn binding(info: &DocumentInfo) -> Binding<'_> {
    Binding {
        document_id: &info.binding.document_id,
        event_id: &info.binding.event_id,
        revision: &info.binding.revision,
    }
}
#[derive(Serialize)]
struct DonorEvidence<'a> {
    kind: &'static str,
    authority: &'static str,
    source: Binding<'a>,
    independent_target: Binding<'a>,
    report: &'a EditTemplateReport,
    independent_check: &'a EditTemplateCheck,
}
#[derive(Serialize)]
struct GuideEvidence<'a> {
    kind: &'static str,
    source: Binding<'a>,
    report: &'a GuideReport,
}
#[derive(Serialize)]
struct SearchEvidence<'a> {
    kind: &'static str,
    source: Binding<'a>,
    reference_events: Vec<Binding<'a>>,
    selected_guide: &'a GuideCandidate,
    guide_options: &'a GuideOptions,
    report: &'a OffTargetReport,
}
fn donor_json(request: &DonorRequest) -> Result<Vec<u8>, String> {
    json_bytes(&DonorEvidence {
        kind: "crispr-donor-review-v1",
        authority: "Historical producer evidence; native restoration does not restore current verification or establish cellular editing",
        source: binding(&request.source),
        independent_target: binding(&request.target),
        report: &request.report,
        independent_check: &request.checked,
    })
}

#[component]
pub fn CrisprPanel(
    source: Memo<DocumentInfo>,
    documents: Memo<Vec<DocumentInfo>>,
    get_document: Callback<String, Option<Molecule>>,
    epoch: RwSignal<u64>,
    on_navigate: Callback<ReadNavigation>,
    on_record: Callback<DonorRequest, Result<(), String>>,
) -> impl IntoView {
    let preset = RwSignal::new("spcas9".to_string());
    let custom_id = RwSignal::new("Custom nuclease rule".to_string());
    let custom_version = RwSignal::new("user-declared-v1".to_string());
    let custom_pam = RwSignal::new("NGG".to_string());
    let custom_side = RwSignal::new("three".to_string());
    let custom_length = RwSignal::new("20".to_string());
    let custom_top = RwSignal::new("17".to_string());
    let custom_bottom = RwSignal::new("17".to_string());
    let min_gc = RwSignal::new("0".to_string());
    let max_gc = RwSignal::new("100".to_string());
    let homopolymer = RwSignal::new("64".to_string());
    let guide_limit = RwSignal::new("4096".to_string());
    let guide_budget = RwSignal::new("20000000".to_string());
    let options = Memo::new(move |_| -> Result<GuideOptions, String> {
        let profile = match preset.get().as_str() {
            "custom" => NucleaseProfile {
                id: custom_id.get(),
                version: custom_version.get(),
                guide_length: whole(&custom_length.get(), "custom guide length")?,
                pam: custom_pam.get().trim().to_uppercase(),
                pam_side: if custom_side.get() == "five" {
                    PamSide::FivePrime
                } else {
                    PamSide::ThreePrime
                },
                nominal_top_offset: whole(&custom_top.get(), "nominal top offset")?,
                nominal_bottom_offset: whole(&custom_bottom.get(), "nominal bottom offset")?,
                custom: true,
                sources: vec![],
            },
            id => crispr::catalogue()
                .into_iter()
                .find(|profile| profile.id == id)
                .ok_or("Choose an installed nuclease profile or declare a custom rule.")?,
        };
        Ok(GuideOptions {
            profile,
            min_gc_bp: percent(&min_gc.get(), "minimum GC")?,
            max_gc_bp: percent(&max_gc.get(), "maximum GC")?,
            max_homopolymer: whole(&homopolymer.get(), "homopolymer limit")?,
            max_candidates: whole(&guide_limit.get(), "guide count limit")?,
            max_comparisons: whole(&guide_budget.get(), "guide comparison budget")?,
        })
    });
    let guides = RwSignal::new(None::<Arc<BoundGuides>>);
    let selected = RwSignal::new(None::<String>);
    let guide_page = RwSignal::new(0usize);
    let rejected = RwSignal::new(false);
    let expanded = RwSignal::new(false);
    let guide_status = RwSignal::new(String::new());
    let guide_job = JobSlot::new();
    let search_job = JobSlot::new();
    let donor_job = JobSlot::new();
    let alive = Arc::new(AtomicBool::new(true));
    let cleanup = alive.clone();
    on_cleanup(move || {
        cleanup.store(false, Ordering::Relaxed);
        guide_job.stop();
        search_job.stop();
        donor_job.stop();
    });
    let selected_guide = Memo::new(move |_| {
        guides.with(|data| {
            data.as_ref().and_then(|b| {
                selected
                    .get()
                    .and_then(|id| b.report.candidates.iter().find(|c| c.id == id).cloned())
            })
        })
    });
    let passing_indices = Memo::new(move |_| {
        guides.with(|data| {
            data.as_ref()
                .map(|b| {
                    b.report
                        .candidates
                        .iter()
                        .enumerate()
                        .filter(|(_, c)| rejected.get() || c.passes_filters)
                        .map(|(i, _)| i)
                        .collect::<Vec<_>>()
                })
                .unwrap_or_default()
        })
    });
    let selected_refs = RwSignal::new(Vec::<String>::new());
    let reference_page = RwSignal::new(0usize);
    let builds = RwSignal::new(BTreeMap::<String, String>::new());
    let accessions = RwSignal::new(BTreeMap::<String, String>::new());
    let scope =
        RwSignal::new("Explicitly selected records; no whole-genome coverage claim".to_string());
    let intended = RwSignal::new(false);
    let mismatches = RwSignal::new("3".to_string());
    let dna_bulges = RwSignal::new("0".to_string());
    let rna_bulges = RwSignal::new("0".to_string());
    let total_bulges = RwSignal::new("0".to_string());
    let alignment_limit = RwSignal::new("128".to_string());
    let column_limit = RwSignal::new("262144".to_string());
    let separate_pam = RwSignal::new(false);
    let screening_pam = RwSignal::new(String::new());
    let screening_rationale = RwSignal::new(String::new());
    let screening_sources = RwSignal::new(String::new());
    let search_budget = RwSignal::new("100000000".to_string());
    let index_limit = RwSignal::new("1000000".to_string());
    let hit_limit = RwSignal::new("4096".to_string());
    let search = RwSignal::new(None::<Arc<BoundSearch>>);
    let search_status = RwSignal::new(String::new());
    let hit_page = RwSignal::new(0usize);
    let selected_hit = RwSignal::new(None::<String>);
    let reference_documents = Memo::new(move |_| -> Result<Vec<DocumentInfo>, String> {
        documents.with(|docs| {
            selected_refs
                .get()
                .iter()
                .map(|id| {
                    let info = docs
                        .iter()
                        .find(|d| d.binding.document_id == *id)
                        .ok_or("A selected reference is no longer open.")?;
                    if info.alphabet != Alphabet::Dna {
                        return Err("Off-target references must be DNA.".into());
                    }
                    Ok(info.clone())
                })
                .collect()
        })
    });
    let references = Memo::new(move |_| -> Result<Vec<ReferenceBinding>, String> {
        Ok(reference_documents
            .get()?
            .into_iter()
            .map(|info| {
                let id = &info.binding.document_id;
                let build = builds.with(|v| v.get(id).and_then(|v| optional(v)));
                let accession = accessions.with(|v| v.get(id).and_then(|v| optional(v)));
                ReferenceBinding {
                    document: info,
                    build,
                    accession,
                }
            })
            .collect())
    });
    let search_options = Memo::new(move |_| -> Result<OffTargetOptions, String> {
        let profile = options.get()?.profile;
        let refs = references.get()?;
        let current = source.get();
        let locus = if intended.get() {
            let index=refs.iter().position(|r|r.document.binding==current.binding).ok_or("Select the exact source event as a reference before marking its intended locus.")?;
            let guide = selected_guide.get().ok_or("Select a guide first.")?;
            Some(IntendedTarget {
                reference_id: format!("reference_{}", index + 1),
                reference_revision: current.binding.revision.clone(),
                guide_positions: guide.guide_positions,
                pam_positions: guide.pam_positions,
                reverse: guide.reverse,
            })
        } else {
            None
        };
        Ok(OffTargetOptions {
            profile,
            bulges: BulgePolicy {
                max_dna_bulge_bases: whole(&dna_bulges.get(), "DNA bulge bases")?,
                max_rna_bulge_bases: whole(&rna_bulges.get(), "RNA bulge bases")?,
                max_total_bulge_bases: whole(&total_bulges.get(), "total bulge bases")?,
                max_alignments_per_hit: whole(&alignment_limit.get(), "alignments per hit")?,
                max_retained_alignment_columns: whole(
                    &column_limit.get(),
                    "retained alignment columns",
                )?,
            },
            pam_policy: if separate_pam.get() {
                PamPolicy::DeclaredScreening {
                    pam: screening_pam.get().trim().to_uppercase(),
                    rationale: screening_rationale.get(),
                    sources: screening_sources
                        .get()
                        .lines()
                        .map(str::trim)
                        .filter(|line| !line.is_empty())
                        .map(str::to_string)
                        .collect(),
                }
            } else {
                PamPolicy::DesignProfile
            },
            max_mismatches: whole(&mismatches.get(), "maximum mismatches")?,
            intended: locus,
            reference_scope: scope.get(),
            max_comparisons: whole(&search_budget.get(), "search comparison budget")?,
            max_index_entries: whole(&index_limit.get(), "index entry limit")?,
            max_hits: whole(&hit_limit.get(), "hit limit")?,
        })
    });
    let target_id = RwSignal::new(String::new());
    let edit_start = RwSignal::new("1".to_string());
    let delete_length = RwSignal::new("0".to_string());
    let replacement = RwSignal::new("A".to_string());
    let left_arm = RwSignal::new("30".to_string());
    let right_arm = RwSignal::new("30".to_string());
    let reverse_donor = RwSignal::new(false);
    let cut_distance = RwSignal::new("100".to_string());
    let donor = RwSignal::new(None::<Arc<DonorRequest>>);
    let donor_status = RwSignal::new(String::new());
    let target = Memo::new(move |_| {
        documents.with(|docs| {
            docs.iter()
                .find(|d| {
                    d.binding.document_id == target_id.get()
                        && d.binding.document_id != source.get().binding.document_id
                        && d.alphabet == Alphabet::Dna
                })
                .cloned()
        })
    });
    let donor_options = Memo::new(move |_| -> Result<EditTemplateOptions, String> {
        Ok(EditTemplateOptions {
            edit_start: whole::<usize>(&edit_start.get(), "1-based edit start")?
                .checked_sub(1)
                .ok_or("Edit start must be at least 1.")?,
            delete_length: whole(&delete_length.get(), "deleted length")?,
            replacement: replacement
                .get()
                .split_whitespace()
                .collect::<String>()
                .to_uppercase(),
            left_arm: whole(&left_arm.get(), "left arm length")?,
            right_arm: whole(&right_arm.get(), "right arm length")?,
            reverse_donor: reverse_donor.get(),
            max_cut_distance: whole(&cut_distance.get(), "maximum nominal cut distance")?,
        })
    });
    Effect::new(move |_| {
        let _ = (source.get(), options.get(), epoch.get());
        guide_job.stop();
        guides.set(None);
        selected.set(None);
        guide_page.set(0);
        guide_status.set(String::new());
    });
    Effect::new(move |_| {
        let _ = (
            selected_guide.get(),
            references.get(),
            search_options.get(),
            epoch.get(),
        );
        search_job.stop();
        search.set(None);
        hit_page.set(0);
        selected_hit.set(None);
        search_status.set(String::new());
    });
    Effect::new(move |_| {
        let _ = (
            selected_guide.get(),
            source.get(),
            target.get(),
            donor_options.get(),
            options.get(),
            epoch.get(),
        );
        donor_job.stop();
        donor.set(None);
        donor_status.set(String::new());
    });
    Effect::new(move |_| {
        let _ = rejected.get();
        guide_page.set(0);
    });
    let guide_guard = Callback::new(move |data: Arc<BoundGuides>| -> Result<(), String> {
        if data.source != source.get_untracked()
            || data.epoch != epoch.get_untracked()
            || options.get_untracked() != Ok(data.report.options.clone())
        {
            return Err("The source event, profile or guide policy changed. Design again.".into());
        }
        Ok(())
    });
    let donor_guard = Callback::new(move |data: Arc<DonorRequest>| -> Result<(), String> {
        let guide = selected_guide
            .get_untracked()
            .ok_or("Select a current guide.")?;
        if data.source != source.get_untracked()
            || target.get_untracked() != Some(data.target.clone())
            || data.epoch != epoch.get_untracked()
            || donor_options.get_untracked() != Ok(data.report.options.clone())
            || data.report.guide.as_ref().is_none_or(|g| {
                g.candidate != guide || options.get_untracked() != Ok(g.options.clone())
            })
        {
            return Err("Source, independent target, guide or donor options changed. Design and check again.".into());
        }
        Ok(())
    });
    let design = {
        let alive = alive.clone();
        move |_| {
            let result = (|| -> Result<_, String> {
                let opts = options.get_untracked()?;
                let info = source.get_untracked();
                if info.alphabet != Alphabet::Dna {
                    return Err("Guide design requires DNA.".into());
                }
                let molecule = document(get_document, &info)?;
                let revision =
                    crispr::guide_revision(&molecule, &opts).map_err(|e| e.to_string())?;
                if molecule.revision() != info.binding.revision {
                    return Err("Source event and molecule do not match.".into());
                }
                let (w, token) = guide_job.start()?;
                Ok((
                    opts,
                    info,
                    molecule,
                    revision,
                    w,
                    token,
                    epoch.get_untracked(),
                ))
            })();
            let (opts, info, molecule, revision, w, token, expected_epoch) = match result {
                Ok(v) => v,
                Err(e) => {
                    guide_status.set(e);
                    return;
                }
            };
            guides.set(None);
            guide_status.set("Scanning both strands for guide/PAM geometry…".into());
            let alive = alive.clone();
            spawn_local(async move {
                let result = w
                    .ask(Request {
                        id: token,
                        revision: revision.clone(),
                        job: AnalysisJob::DesignCrisprGuides {
                            molecule: molecule.clone(),
                            options: opts.clone(),
                        },
                    })
                    .await;
                if !alive.load(Ordering::Relaxed)
                    || !guide_job.current(token)
                    || source.get_untracked() != info
                    || options.get_untracked() != Ok(opts.clone())
                    || epoch.get_untracked() != expected_epoch
                {
                    return;
                }
                guide_job.busy.set(false);
                let result = (|| -> Result<GuideReport, String> {
                    let AnalysisResult::CrisprGuides(report) = response(result, token, &revision)?
                    else {
                        return Err("Unexpected guide worker result.".into());
                    };
                    report
                        .validate_bindings_for(&molecule, &opts)
                        .map_err(|e| e.to_string())?;
                    Ok(report)
                })();
                match result {
                    Ok(report) => {
                        selected.set(
                            report
                                .candidates
                                .iter()
                                .find(|c| c.passes_filters)
                                .map(|c| c.id.clone()),
                        );
                        guide_status.set(format!(
                            "{} candidates · {} scan · no empirical activity score",
                            report.candidates.len(),
                            if report.complete {
                                "complete"
                            } else {
                                "incomplete"
                            }
                        ));
                        guides.set(Some(Arc::new(BoundGuides {
                            source: info,
                            report: Arc::new(report),
                            epoch: expected_epoch,
                        })));
                    }
                    Err(e) => guide_status.set(e),
                }
            });
        }
    };
    let search_run = {
        let alive = alive.clone();
        move |_| {
            let prepared = (|| -> Result<_, String> {
                let base = guides.get_untracked().ok_or("Design guides first.")?;
                let guide_options = base.report.options.clone();
                guide_guard.run(base)?;
                let guide = selected_guide
                    .get_untracked()
                    .ok_or("Select a guide first.")?;
                if !guide.concrete || !guide.pam_certain {
                    return Err(
                        "Off-target querying requires a concrete selected guide and certain PAM."
                            .into(),
                    );
                }
                let refs = references.get_untracked()?;
                if refs.is_empty() || refs.len() > 64 {
                    return Err("Explicitly select 1–64 reference records.".into());
                }
                refs.iter().try_fold(0usize, |total, reference| {
                    total
                        .checked_add(reference.document.length)
                        .filter(|total| *total <= 20_000_000)
                        .ok_or("Selected reference sequences exceed 20 million bases.")
                })?;
                let opts = search_options.get_untracked()?;
                opts.validate().map_err(|error| error.to_string())?;
                let inputs = refs
                    .iter()
                    .enumerate()
                    .map(|(i, r)| {
                        let molecule = document(get_document, &r.document)?;
                        if molecule.revision() != r.document.binding.revision {
                            return Err(
                                "A reference no longer matches its exact bound event.".into()
                            );
                        }
                        Ok(ReferenceInput {
                            id: format!("reference_{}", i + 1),
                            molecule,
                            build: r.build.clone(),
                            accession: r.accession.clone(),
                        })
                    })
                    .collect::<Result<Vec<_>, String>>()?;
                let query = GuideQuery {
                    id: guide.id.clone(),
                    sequence: guide.guide_sequence.clone(),
                };
                let revision = crispr::off_target_revision(&query, &inputs, &opts)
                    .map_err(|e| e.to_string())?;
                let (w, token) = search_job.start()?;
                Ok((
                    source.get_untracked(),
                    guide,
                    guide_options,
                    refs,
                    inputs,
                    query,
                    opts,
                    revision,
                    w,
                    token,
                    epoch.get_untracked(),
                ))
            })();
            let (
                info,
                guide,
                guide_options,
                refs,
                inputs,
                query,
                opts,
                revision,
                w,
                token,
                expected_epoch,
            ) = match prepared {
                Ok(v) => v,
                Err(e) => {
                    search_status.set(e);
                    return;
                }
            };
            search.set(None);
            search_status.set("Searching only the explicitly selected records…".into());
            let alive = alive.clone();
            spawn_local(async move {
                let result = w
                    .ask(Request {
                        id: token,
                        revision: revision.clone(),
                        job: AnalysisJob::SearchCrisprOffTargets {
                            query: query.clone(),
                            references: inputs.clone(),
                            options: opts.clone(),
                        },
                    })
                    .await;
                if !alive.load(Ordering::Relaxed)
                    || !search_job.current(token)
                    || source.get_untracked() != info
                    || selected_guide.get_untracked() != Some(guide.clone())
                    || options.get_untracked() != Ok(guide_options.clone())
                    || references.get_untracked() != Ok(refs.clone())
                    || search_options.get_untracked() != Ok(opts.clone())
                    || epoch.get_untracked() != expected_epoch
                {
                    return;
                }
                search_job.busy.set(false);
                let result = (|| -> Result<OffTargetReport, String> {
                    let AnalysisResult::CrisprOffTargets(report) =
                        response(result, token, &revision)?
                    else {
                        return Err("Unexpected off-target worker result.".into());
                    };
                    report
                        .validate_bindings_for(&query, &inputs, &opts)
                        .map_err(|e| e.to_string())?;
                    Ok(report)
                })();
                match result {
                    Ok(report) => {
                        search_status.set(format!(
                            "{} retained matches · {} search of selected records only",
                            report.hits.len(),
                            if report.complete {
                                "complete"
                            } else {
                                "incomplete"
                            }
                        ));
                        search.set(Some(Arc::new(BoundSearch {
                            source: info,
                            guide,
                            guide_options,
                            references: refs,
                            report: Arc::new(report),
                            epoch: expected_epoch,
                        })));
                    }
                    Err(e) => search_status.set(e),
                }
            });
        }
    };
    let donor_design = {
        let alive = alive.clone();
        move |_| {
            let prepared = (|| -> Result<_, String> {
                let base = guides.get_untracked().ok_or("Design guides first.")?;
                guide_guard.run(base.clone())?;
                let candidate = selected_guide
                    .get_untracked()
                    .filter(|g| g.passes_filters)
                    .ok_or("Select an eligible guide first.")?;
                let target = target
                    .get_untracked()
                    .ok_or("Choose a separate open DNA document as the intended edited target.")?;
                let original = source.get_untracked();
                let source_molecule = document(get_document, &original)?;
                let target_molecule = document(get_document, &target)?;
                let opts = donor_options.get_untracked()?;
                let guide = Some(GuideBinding {
                    options: base.report.options.clone(),
                    candidate,
                });
                let revision = edit_template::batch_revision(
                    &source_molecule,
                    &target_molecule,
                    &opts,
                    guide.as_ref(),
                )
                .map_err(|e| e.to_string())?;
                if source_molecule.revision() != original.binding.revision
                    || target_molecule.revision() != target.binding.revision
                {
                    return Err("A donor input does not match its bound source event.".into());
                }
                let (w, token) = donor_job.start()?;
                Ok((
                    original,
                    target,
                    source_molecule,
                    target_molecule,
                    opts,
                    guide,
                    revision,
                    w,
                    token,
                    epoch.get_untracked(),
                ))
            })();
            let (
                original,
                target_info,
                source_molecule,
                target_molecule,
                opts,
                guide,
                revision,
                w,
                token,
                expected_epoch,
            ) = match prepared {
                Ok(v) => v,
                Err(e) => {
                    donor_status.set(e);
                    return;
                }
            };
            donor.set(None);
            donor_status.set("Constructing literal homology arms and independently checking the intended target…".into());
            let alive = alive.clone();
            spawn_local(async move {
                let result = w
                    .ask(Request {
                        id: token,
                        revision: revision.clone(),
                        job: AnalysisJob::DesignEditTemplate {
                            source: source_molecule.clone(),
                            target: target_molecule.clone(),
                            options: opts.clone(),
                            guide: guide.clone(),
                        },
                    })
                    .await;
                if !alive.load(Ordering::Relaxed)
                    || !donor_job.current(token)
                    || source.get_untracked() != original
                    || target.get_untracked() != Some(target_info.clone())
                    || donor_options.get_untracked() != Ok(opts.clone())
                    || epoch.get_untracked() != expected_epoch
                    || selected_guide.get_untracked() != guide.as_ref().map(|g| g.candidate.clone())
                    || guide
                        .as_ref()
                        .is_none_or(|g| options.get_untracked() != Ok(g.options.clone()))
                {
                    return;
                }
                donor_job.busy.set(false);
                let result = (|| -> Result<EditTemplateReport, String> {
                    let AnalysisResult::EditTemplate(report) = response(result, token, &revision)?
                    else {
                        return Err("Unexpected donor design result.".into());
                    };
                    report
                        .validate_bindings_for(
                            &source_molecule,
                            &target_molecule,
                            &opts,
                            guide.as_ref(),
                        )
                        .map_err(|e| e.to_string())?;
                    Ok(report)
                })();
                match result {
                    Ok(report) => {
                        donor_status.set(if report.check.consistent{"Literal donor and independent target are consistent under this sequence model. Cellular editing is not predicted.".into()}else{"The independent sequence check did not pass. Inspect target, arms and nominal cut distance; release remains disabled.".into()});
                        let checked = report.check.clone();
                        donor.set(Some(Arc::new(DonorRequest {
                            source: original,
                            target: target_info,
                            report: Arc::new(report),
                            checked,
                            edited: false,
                            epoch: expected_epoch,
                        })));
                    }
                    Err(e) => donor_status.set(e),
                }
            });
        }
    };
    let release = {
        let alive = alive.clone();
        Callback::new(move |kind: &'static str| {
            let Some(old) = donor.get_untracked() else {
                return;
            };
            if let Err(e) = donor_guard.run(old.clone()) {
                donor_status.set(e);
                return;
            }
            if !old.checked.consistent {
                donor_status
                    .set("A failed independent donor check cannot release a product.".into());
                return;
            }
            let source_molecule = match document(get_document, &old.source) {
                Ok(v) => v,
                Err(e) => {
                    donor_status.set(e);
                    return;
                }
            };
            let target_molecule = match document(get_document, &old.target) {
                Ok(v) => v,
                Err(e) => {
                    donor_status.set(e);
                    return;
                }
            };
            let (w, token) = match donor_job.start() {
                Ok(v) => v,
                Err(e) => {
                    donor_status.set(e);
                    return;
                }
            };
            donor_status.set(
                "Rebuilding the donor and running a fresh independent check before release…".into(),
            );
            let alive = alive.clone();
            spawn_local(async move {
                let fresh = replay_donor(&source_molecule, &target_molecule, &old, &w).await;
                if !alive.load(Ordering::Relaxed) || !donor_job.current(token) {
                    return;
                }
                donor_job.busy.set(false);
                if let Err(e) = donor_guard.run(old.clone()) {
                    donor_status.set(e);
                    return;
                }
                let mut fresh = match fresh {
                    Ok(v) => v,
                    Err(e) => {
                        donor_status.set(e);
                        return;
                    }
                };
                match kind {
                    "open-donor" | "open-edited" => {
                        fresh.edited = kind == "open-edited";
                        if let Err(e) = on_record.run(fresh) {
                            donor_status.set(e);
                        }
                    }
                    "fasta" | "genbank" => {
                        let result = (|| -> Result<crate::io::download::ExportArtifact, String> {
                            let product = producer_product(&fresh, false)?;
                            if kind == "fasta" {
                                ExportArtifact::prepare(
                                    "crispr-donor.fasta",
                                    "text/plain",
                                    collagene_core::formats::sequence::export_fasta(&product)
                                        .as_bytes(),
                                )
                            } else {
                                ExportArtifact::prepare(
                                    "crispr-donor.gb",
                                    "text/plain",
                                    collagene_core::formats::sequence::export_genbank(&product)
                                        .map_err(|e| e.to_string())?
                                        .as_bytes(),
                                )
                            }
                        })();
                        let dispatch_old = old.clone();
                        crate::io::download::complete_download_guarded(
                            result,
                            move || donor_guard.run(dispatch_old.clone()),
                            move |result| {
                                if donor_status.try_get_untracked().is_none()
                                    || donor_guard.run(old.clone()).is_err()
                                {
                                    return;
                                }
                                match result{Ok(())=>donor_status.set("Exported after fresh source, guide and independent-target reconstruction/check. Keep JSON/native evidence for provenance.".into()),Err(e)=>donor_status.set(e)}
                            },
                        );
                    }
                    _ => {}
                }
            });
        })
    };
    let export_guides = move |_| {
        if let Some(data) = guides.get_untracked() {
            let result = guide_guard
                .run(data.clone())
                .and_then(|_| {
                    json_bytes(&GuideEvidence {
                        kind: "crispr-guide-review-v1",
                        source: binding(&data.source),
                        report: data.report.as_ref(),
                    })
                })
                .and_then(|bytes| {
                    ExportArtifact::prepare("crispr-guides.json", "application/json", &bytes)
                });
            crate::io::download::complete_download(result, move |result| {
                if guide_status.try_get_untracked().is_some() {
                    if let Err(e) = result {
                        guide_status.set(e);
                    }
                }
            });
        }
    };
    let export_search = Callback::new(move |(): ()| {
        if let Some(data) = search.get_untracked() {
            let valid = data.source == source.get_untracked()
                && data.epoch == epoch.get_untracked()
                && selected_guide.get_untracked() == Some(data.guide.clone())
                && options.get_untracked() == Ok(data.guide_options.clone())
                && references.get_untracked() == Ok(data.references.clone())
                && search_options.get_untracked() == Ok(data.report.options.clone());
            if !valid {
                search_status
                    .set("Reference roles, labels, source or query changed. Search again.".into());
                return;
            }
            let refs = data
                .references
                .iter()
                .map(|r| binding(&r.document))
                .collect::<Vec<_>>();
            let result = json_bytes(&SearchEvidence {
                kind: "crispr-off-target-review-v1",
                source: binding(&data.source),
                reference_events: refs,
                selected_guide: &data.guide,
                guide_options: &data.guide_options,
                report: data.report.as_ref(),
            })
            .and_then(|bytes| {
                ExportArtifact::prepare("crispr-off-targets.json", "application/json", &bytes)
            });
            crate::io::download::complete_download(result, move |result| {
                if search_status.try_get_untracked().is_some() {
                    if let Err(e) = result {
                        search_status.set(e);
                    }
                }
            });
        }
    });
    let export_donor = move |_| {
        if let Some(data) = donor.get_untracked() {
            let result = donor_guard
                .run(data.clone())
                .and_then(|_| donor_json(&data))
                .and_then(|bytes| {
                    ExportArtifact::prepare("crispr-donor-review.json", "application/json", &bytes)
                });
            crate::io::download::complete_download(result, move |result| {
                if donor_status.try_get_untracked().is_some() {
                    if let Err(e) = result {
                        donor_status.set(e);
                    }
                }
            });
        }
    };
    let navigate_guide = Callback::new(move |at: usize| {
        if let Some(data) = guides.get_untracked() {
            if guide_guard.run(data.clone()).is_ok() && at < data.source.length {
                on_navigate.run(ReadNavigation {
                    binding: data.source.binding.clone(),
                    call_index: at,
                    epoch: data.epoch,
                });
            }
        }
    });
    let navigate_hit = Callback::new(move |id: String| {
        if let Some(data) = search.get_untracked() {
            if data.epoch != epoch.get_untracked()
                || data.source != source.get_untracked()
                || selected_guide.get_untracked() != Some(data.guide.clone())
                || options.get_untracked() != Ok(data.guide_options.clone())
                || references.get_untracked() != Ok(data.references.clone())
                || search_options.get_untracked() != Ok(data.report.options.clone())
            {
                return;
            }
            if let Some(hit) = data.report.hits.iter().find(|h| h.id == id) {
                if let (Some(reference), Some(at)) = (
                    data.references.get(hit.reference_index),
                    hit.guide_positions.first(),
                ) {
                    on_navigate.run(ReadNavigation {
                        binding: reference.document.binding.clone(),
                        call_index: *at,
                        epoch: data.epoch,
                    });
                }
            }
        }
    });
    let can_release = Memo::new(move |_| {
        !donor_job.busy.get() && donor.with(|d| d.as_ref().is_some_and(|d| d.checked.consistent))
    });
    view! {<section class="crispr-panel" class:crispr-expanded=move||expanded.get()><div class="crispr-heading"><h2>"CRISPR design"</h2><button class="button secondary" on:click=move|_|expanded.update(|v|*v=!*v)>{move||if expanded.get(){"Close expanded CRISPR"}else{"Expand CRISPR"}}</button></div><p class="helper">"Design guide/PAM geometry, search explicitly selected reference records and review literal donor templates. Nominal cut coordinates are advisory; no cellular editing, repair pathway or empirical efficacy is predicted."</p><p class="crispr-source">{move||source.with(|s|format!("Source: {} · {} nt · {}",s.name,s.length,if s.topology==Topology::Circular{"circular"}else{"linear"}))}</p>
        <label>"Nuclease rule"
            <select aria-label="CRISPR nuclease" prop:value=move || preset.get()
                on:change=move |e| preset.set(event_target_value(&e))>
                {crispr::catalogue().into_iter().map(|profile| {
                    let label = format!("{} · nominal {} / {} nt", profile.id, profile.pam, profile.guide_length);
                    view! {<option value=profile.id>{label}</option>}
                }).collect_view()}
                <option value="custom">"Explicit custom PAM and geometry"</option>
            </select>
        </label>
        <Show when=move||preset.get()=="custom"><fieldset><legend>"Custom rule (user declared)"</legend><label>"Identity"<input aria-label="CRISPR custom identity" prop:value=move||custom_id.get() on:input=move|e|custom_id.set(event_target_value(&e))/></label><label>"Version"<input aria-label="CRISPR custom version" prop:value=move||custom_version.get() on:input=move|e|custom_version.set(event_target_value(&e))/></label><label>"PAM (IUPAC)"<input aria-label="CRISPR custom PAM" prop:value=move||custom_pam.get() on:input=move|e|custom_pam.set(event_target_value(&e))/></label><label>"PAM side"<select aria-label="CRISPR custom PAM side" prop:value=move||custom_side.get() on:change=move|e|custom_side.set(event_target_value(&e))><option value="three">"3′ of guide"</option><option value="five">"5′ of guide"</option></select></label><label>"Guide length"<input aria-label="CRISPR custom length" prop:value=move||custom_length.get() on:input=move|e|custom_length.set(event_target_value(&e))/></label><label>"Nominal top boundary offset"<input aria-label="CRISPR custom top offset" prop:value=move||custom_top.get() on:input=move|e|custom_top.set(event_target_value(&e))/></label><label>"Nominal bottom boundary offset"<input aria-label="CRISPR custom bottom offset" prop:value=move||custom_bottom.get() on:input=move|e|custom_bottom.set(event_target_value(&e))/></label><p class="helper">"Offsets count bases from the oriented DNA spacer’s 5′ boundary; reverse orientation maps and swaps the source-strand labels. Custom spacers are 12–64 nt and PAMs are 1–16 nt. These nominal rules have no implied experimental calibration."</p></fieldset></Show>
        <details><summary>"Guide filters and scan limits"</summary><div class="crispr-options"><label>"Minimum GC%"<input aria-label="CRISPR minimum GC" prop:value=move||min_gc.get() on:input=move|e|min_gc.set(event_target_value(&e))/></label><label>"Maximum GC%"<input aria-label="CRISPR maximum GC" prop:value=move||max_gc.get() on:input=move|e|max_gc.set(event_target_value(&e))/></label><label>"Longest homopolymer"<input aria-label="CRISPR homopolymer limit" prop:value=move||homopolymer.get() on:input=move|e|homopolymer.set(event_target_value(&e))/></label><label>"Candidate limit"<input aria-label="CRISPR candidate limit" prop:value=move||guide_limit.get() on:input=move|e|guide_limit.set(event_target_value(&e))/></label><label>"Comparison budget"<input aria-label="CRISPR guide budget" prop:value=move||guide_budget.get() on:input=move|e|guide_budget.set(event_target_value(&e))/></label></div></details>
        <div class="crispr-actions"><button class="button primary" disabled=move||guide_job.busy.get()||source.with(|s|s.alphabet!=Alphabet::Dna) on:click=design>"Design CRISPR guides"</button><button class="button secondary" disabled=move||!guide_job.busy.get() on:click=move|_|{guide_job.stop();guide_status.set("Guide design cancelled.".into());}>"Cancel guide design"</button></div><p class="crispr-guide-status" role="status">{move||guide_status.get()}</p>
        {move||guides.get().map(|data|view!{<div class="crispr-guides"><p>{format!("Profile {} · {} · engine {} · {} candidates · {} comparisons",data.report.options.profile.id,data.report.options.profile.version,data.report.engine_version,data.report.candidates.len(),data.report.comparisons)}</p><details><summary>"Profile and source binding"</summary><p>{format!("Source document {} · event {} · revision {}",data.source.binding.document_id,data.source.binding.event_id,data.source.binding.revision)}</p><p>{format!("PAM {} on {} side · {} nt · nominal top/bottom offsets {}/{}",data.report.options.profile.pam,if data.report.options.profile.pam_side==PamSide::ThreePrime{"3′"}else{"5′"},data.report.options.profile.guide_length,data.report.options.profile.nominal_top_offset,data.report.options.profile.nominal_bottom_offset)}</p>{data.report.options.profile.sources.iter().cloned().map(|url|view!{<a href=url.clone() target="_blank" rel="noopener noreferrer">{url.clone()}</a>}).collect_view()}</details><button class="button secondary" on:click=export_guides>"Download guide report"</button><label class="checkbox-label"><input type="checkbox" aria-label="CRISPR show rejected guides" prop:checked=move||rejected.get() on:change=move|e|rejected.set(event_target_checked(&e))/>"Include rejected and ambiguous guide diagnostics"</label>
            {move||guides.get().map(|b|passing_indices.get().into_iter().skip(guide_page.get()*25).take(25).filter_map(|i|b.report.candidates.get(i).map(|g|{let id=g.id.clone();let attr_id=id.clone();let selected_id=id.clone();view!{<button class="button secondary crispr-guide-entry" data-guide-id=attr_id class:active=move||selected.get()==Some(selected_id.clone()) on:click=move|_|selected.set(Some(id.clone()))>{format!("{} · PAM {} · {} · {} · GC {}",g.guide_sequence,g.pam_sequence,if g.reverse{"reverse"}else{"forward"},if g.passes_filters{"passes declared filters"}else{"diagnostic only"},g.gc_bp.map(|v|format!("{:.2}%",v as f64/100.)).unwrap_or_else(||"uncertain".into()))}</button>}})).collect_view())}
            <div class="crispr-actions"><button class="button secondary" disabled=move||guide_page.get()==0 on:click=move|_|guide_page.update(|p|*p=p.saturating_sub(1))>"Previous guides"</button><button class="button secondary" disabled=move||{(guide_page.get()+1)*25>=passing_indices.with(Vec::len)} on:click=move|_|guide_page.update(|p|*p+=1)>"Next guides"</button></div>{data.report.limits.iter().map(|v|view!{<p class="crispr-limit">{v.clone()}</p>}).collect_view()}{data.report.warnings.iter().map(|v|view!{<p class="helper">{v.clone()}</p>}).collect_view()}</div>})}
        {move||selected_guide.get().map(|g|{let start=g.guide_positions.first().copied();view!{<article class="crispr-guide-detail"><h3>"Selected guide"</h3><p><strong>{g.guide_sequence.clone()}</strong>" · DNA spacer (5′→3′), not a complete guide RNA/scaffold"</p><p>{format!("PAM {} · {} · {} · homopolymer {}",g.pam_sequence,if g.reverse{"reverse strand"}else{"forward strand"},if g.pam_certain{"certain PAM"}else{"possible PAM"},g.longest_homopolymer.map(|n|n.to_string()).unwrap_or_else(||"uncertain".into()))}</p><p>{format!("Guide source bases {} · PAM source bases {} (1-based, ordered)",positions(&g.guide_positions),positions(&g.pam_positions))}</p><p>{format!("Nominal top boundary: {}; nominal bottom boundary: {}. These advisory coordinates do not establish physical cut products.",boundary(g.nominal_top_cut),boundary(g.nominal_bottom_cut))}</p>{g.rejections.into_iter().map(|r|view!{<p class="helper">{r}</p>}).collect_view()}<button class="button secondary" disabled=start.is_none() on:click=move|_|if let Some(at)=start{navigate_guide.run(at)}>"Open guide 5′ position"</button></article>}})}
        <h3>"Search selected reference records"</h3><p class="helper">"Choose named DNA records and optional build/accession labels. Limits are 64 records, 20 million aggregate bases and 16 MiB of serialized reference data. The declared global alignment and screening-PAM policy define the searched sites, including bounded DNA/RNA gaps when enabled. Chromatin, expression, whole-genome completeness and cellular activity are not assessed. Intended-site tagging keeps all retained paths and additional loci."</p><details class="crispr-reference-picker"><summary>"Reference records and labels"</summary><button class="button secondary" on:click=move|_|selected_refs.set(vec![source.get_untracked().binding.document_id])>"Use source as sole reference"</button><button class="button secondary" on:click=move|_|selected_refs.set(vec![])>"Clear CRISPR references"</button>{move||documents.with(|docs|docs.iter().skip(reference_page.get()*25).take(25).map(|d|{let id=d.binding.document_id.clone();let check_id=id.clone();let changed=id.clone();view!{<label class="crispr-reference-choice"><input type="checkbox" aria-label=format!("CRISPR reference {id}") disabled=d.alphabet!=Alphabet::Dna prop:checked=move||selected_refs.with(|v|v.contains(&check_id)) on:change=move|e|{let yes=event_target_checked(&e);selected_refs.update(|v|if yes{if !v.contains(&changed){v.push(changed.clone());}}else{v.retain(|id|id!=&changed);});}/><span><strong>{d.name.clone()}</strong><small>{format!("{} nt · {}",d.length,d.binding.document_id)}</small></span></label>}}).collect_view())}<button class="button secondary" disabled=move||reference_page.get()==0 on:click=move|_|reference_page.update(|p|*p=p.saturating_sub(1))>"Previous references"</button><button class="button secondary" disabled=move||{(reference_page.get()+1)*25>=documents.with(Vec::len)} on:click=move|_|reference_page.update(|p|*p+=1)>"Next references"</button>
        {move||reference_documents.get().ok().map(|refs|refs.into_iter().enumerate().map(|(i,r)|{let build_id=r.binding.document_id.clone();let set_build=build_id.clone();let accession_id=build_id.clone();let set_accession=build_id.clone();let build_label=format!("CRISPR build {build_id}");let accession_label=format!("CRISPR accession {accession_id}");view!{<article class="crispr-reference-label"><strong>{format!("reference_{} · {}",i+1,r.name)}</strong><p>{format!("Document {} · event {} · revision {}",r.binding.document_id,r.binding.event_id,r.binding.revision)}</p><label>"Build label (optional)"<input aria-label=build_label prop:value=move||builds.with(|v|v.get(&build_id).cloned().unwrap_or_default()) on:input=move|e|builds.update(|v|{v.insert(set_build.clone(),event_target_value(&e));})/></label><label>"Accession (optional)"<input aria-label=accession_label prop:value=move||accessions.with(|v|v.get(&accession_id).cloned().unwrap_or_default()) on:input=move|e|accessions.update(|v|{v.insert(set_accession.clone(),event_target_value(&e));})/></label></article>}}).collect_view())}</details>
        <label>"Declared reference scope"<input aria-label="CRISPR reference scope" prop:value=move||scope.get() on:input=move|e|scope.set(event_target_value(&e))/></label><label class="checkbox-label"><input type="checkbox" aria-label="CRISPR mark intended locus" prop:checked=move||intended.get() on:change=move|e|intended.set(event_target_checked(&e))/>"Tag only this selected guide locus in the exact source reference role as intended"</label>
        <details>
            <summary>"Mismatch, bulge and search limits"</summary>
            <p class="helper">"Zero bulge limits restrict the model to ungapped alignments. Each kind and the total cap apply together. Bulged searches include terminal and mixed gaps under a full-input global alignment model, not a physical heteroduplex or cleavage prediction."</p>
            <div class="crispr-options">
                <search_view::WholeField label="Maximum mismatches (0–6)" value=mismatches/>
                <search_view::WholeField label="Extra target DNA bases (0–4)" value=dna_bulges/>
                <search_view::WholeField label="Extra guide RNA bases (0–4)" value=rna_bulges/>
                <search_view::WholeField label="Total bulge bases (0–4)" value=total_bulges/>
                <search_view::WholeField label="Retained alignments per hit" value=alignment_limit/>
                <search_view::WholeField label="Retained alignment columns" value=column_limit/>
                <search_view::WholeField label="Comparison budget" value=search_budget/>
                <search_view::WholeField label="Index entries (ungapped search only)" value=index_limit/>
                <search_view::WholeField label="Retained hits" value=hit_limit/>
            </div>
        </details>
        <details>
            <summary>"Screening PAM policy"</summary>
            <label class="checkbox-label">
                <input type="checkbox" prop:checked=move || separate_pam.get()
                    on:change=move |e| separate_pam.set(event_target_checked(&e))/>
                "Declare a separate screening PAM"
            </label>
            <p class="helper">"The design profile keeps its own PAM. A screening mask may deliberately include additional PAM sequences; its length must match the design PAM. These declarations do not establish activity for the extra sites."</p>
            <Show when=move || separate_pam.get()>
                <label>"Screening PAM (IUPAC)"<input maxlength="16"
                    prop:value=move || screening_pam.get()
                    on:input=move |e| screening_pam.set(event_target_value(&e))/></label>
                <label>"Reason for this screening policy"<input maxlength="1024"
                    prop:value=move || screening_rationale.get()
                    on:input=move |e| screening_rationale.set(event_target_value(&e))/></label>
                <label>"Source URLs (one per line)"<textarea maxlength="8192"
                    prop:value=move || screening_sources.get()
                    on:input=move |e| screening_sources.set(event_target_value(&e))/></label>
            </Show>
        </details>
        <div class="crispr-actions">
            <button class="button primary" disabled=move || search_job.busy.get()
                || selected_guide.get().is_none_or(|g| !g.concrete || !g.pam_certain)
                on:click=search_run>"Search CRISPR references"</button>
            <button class="button secondary" disabled=move || !search_job.busy.get()
                on:click=move |_| { search_job.stop(); search_status.set("Reference search cancelled.".into()); }>
                "Cancel reference search"
            </button>
        </div>
        <p class="crispr-search-status" role="status">{move || search_status.get()}</p>
        <search_view::SearchResults
            search page=hit_page selected=selected_hit navigate=navigate_hit export=export_search
        />
        <h3>"Literal homology donor template"</h3><p class="helper">"Provide an intended edited DNA target as a separate document. The donor uses original-source arms and your literal replacement; the target is a noncausal comparison only. Donors are bounded to 100,000 nt, source to 1 million nt and target to 1.1 million nt, with 8 MiB rich-input bounds. This does not model cellular HDR, editing efficiency, repair mechanisms or physical synthesis."</p><label>"Independent intended edited target"<select aria-label="CRISPR independent target" prop:value=move||target_id.get() on:change=move|e|target_id.set(event_target_value(&e))><option value="">"Choose a separate target document"</option>{move||documents.with(|docs|docs.iter().filter(|d|d.alphabet==Alphabet::Dna&&d.binding.document_id!=source.get().binding.document_id).map(|d|view!{<option value=d.binding.document_id.clone()>{format!("{} · {}",d.name,d.binding.document_id)}</option>}).collect_view())}</select></label><div class="crispr-options"><label>"Edit starts before source base (1-based; linear length+1 appends)"<input aria-label="CRISPR edit start" prop:value=move||edit_start.get() on:input=move|e|edit_start.set(event_target_value(&e))/></label><label>"Deleted source bases"<input aria-label="CRISPR delete length" prop:value=move||delete_length.get() on:input=move|e|delete_length.set(event_target_value(&e))/></label><label>"Left arm length"<input aria-label="CRISPR left arm" prop:value=move||left_arm.get() on:input=move|e|left_arm.set(event_target_value(&e))/></label><label>"Right arm length"<input aria-label="CRISPR right arm" prop:value=move||right_arm.get() on:input=move|e|right_arm.set(event_target_value(&e))/></label><label>"Maximum nominal cut-to-edit distance"<input aria-label="CRISPR cut distance" prop:value=move||cut_distance.get() on:input=move|e|cut_distance.set(event_target_value(&e))/></label></div><label>"Literal replacement (empty permits deletion)"<textarea aria-label="CRISPR replacement" prop:value=move||replacement.get() on:input=move|e|replacement.set(event_target_value(&e))></textarea></label><label class="checkbox-label"><input type="checkbox" aria-label="CRISPR reverse donor" prop:checked=move||reverse_donor.get() on:change=move|e|reverse_donor.set(event_target_checked(&e))/>"Write donor as reverse complement"</label><div class="crispr-actions"><button class="button primary" disabled=move||donor_job.busy.get()||selected_guide.get().is_none_or(|g|!g.passes_filters)||target.get().is_none() on:click=donor_design>"Design and check donor"</button><button class="button secondary" disabled=move||!donor_job.busy.get() on:click=move|_|{donor_job.stop();donor_status.set("Donor computation cancelled.".into());}>"Cancel donor computation"</button></div><p class="crispr-donor-status" role="status">{move||donor_status.get()}</p>
        {move||donor.get().map(|data|view!{<article class="crispr-donor-result"><h4>{if data.checked.consistent{"Independent sequence check passed"}else{"Independent sequence check failed"}}</h4><p>{format!("Left arm: {} · replacement: {} · right arm: {} · target: {} · sequence changed: {}",if data.checked.left_matches{"matches"}else{"does not match"},if data.checked.replacement_matches{"matches"}else{"does not match"},if data.checked.right_matches{"matches"}else{"does not match"},if data.checked.target_matches{"matches"}else{"does not match"},if data.checked.changes_sequence{"yes"}else{"no"})}</p><p>{format!("{} donor bases · nominal cut distance {} · {}",data.checked.donor_length,data.checked.guide_cut_distance.map(|v|v.to_string()).unwrap_or_else(||"unavailable".into()),match data.checked.guide_within_distance{Some(true)=>"within declared distance",Some(false)=>"outside declared distance",None=>"no guide-distance assertion"})}</p><p>{format!("Left original source bases {} · right original source bases {} · deleted source bases {}",positions(&data.report.left.source_positions),positions(&data.report.right.source_positions),positions(&data.report.deleted_positions))}</p><pre>{preview(data.report.donor.sequence.as_str())}</pre><details><summary>"Source, guide and noncausal target bindings"</summary><p>{format!("Source document {} · event {} · revision {}",data.source.binding.document_id,data.source.binding.event_id,data.source.binding.revision)}</p><p>{format!("Target requirement document {} · event {} · revision {}",data.target.binding.document_id,data.target.binding.event_id,data.target.binding.revision)}</p><p>{format!("Source-only product ancestry. Engine {} · options {} · batch {}",data.report.engine_version,data.report.options_revision,data.report.batch_revision)}</p></details><p class="helper">"FASTA and GenBank omit full source-event and per-call design evidence. Keep the JSON review/native workspace. Restoring producer provenance never restores current verification."</p><button class="button secondary" on:click=export_donor>"Download donor review JSON"</button><div class="crispr-actions"><button class="button primary" disabled=move||!can_release.get() on:click=move|_|release.run("open-donor")>"Open checked donor template"</button><button class="button secondary" disabled=move||!can_release.get() on:click=move|_|release.run("open-edited")>"Open checked edited reference"</button><button class="button secondary" disabled=move||!can_release.get() on:click=move|_|release.run("fasta")>"Export checked donor FASTA"</button><button class="button secondary" disabled=move||!can_release.get() on:click=move|_|release.run("genbank")>"Export checked donor GenBank"</button></div>{data.report.warnings.iter().chain(data.checked.warnings.iter()).map(|v|view!{<p class="helper">{v.clone()}</p>}).collect_view()}</article>})}
    </section>}
}

async fn replay_donor(
    source: &Molecule,
    target: &Molecule,
    old: &DonorRequest,
    worker: &CachedWorker,
) -> Result<DonorRequest, String> {
    let opts = &old.report.options;
    let guide = old.report.guide.as_ref();
    let revision =
        edit_template::batch_revision(source, target, opts, guide).map_err(|e| e.to_string())?;
    if source.revision() != old.source.binding.revision
        || target.revision() != old.target.binding.revision
    {
        return Err("Source and target no longer match their exact event bindings.".into());
    }
    let result = worker
        .ask(Request {
            id: u32::MAX - 50,
            revision: revision.clone(),
            job: AnalysisJob::DesignEditTemplate {
                source: source.clone(),
                target: target.clone(),
                options: opts.clone(),
                guide: guide.cloned(),
            },
        })
        .await;
    let AnalysisResult::EditTemplate(report) = response(result, u32::MAX - 50, &revision)? else {
        return Err("Unexpected fresh donor reconstruction.".into());
    };
    report
        .validate_bindings_for(source, target, opts, guide)
        .map_err(|e| e.to_string())?;
    if report.donor.revision() != old.report.donor.revision()
        || report.edited_reference.revision() != old.report.edited_reference.revision()
        || report.options != old.report.options
        || report.guide != old.report.guide
    {
        return Err("Fresh donor reconstruction differs from the displayed design. Inspect and check the current design again.".into());
    }
    let donor = report.donor.sequence.as_str().to_string();
    let revision = edit_template::check_batch_revision(source, target, opts, guide, &donor)
        .map_err(|e| e.to_string())?;
    let result = worker
        .ask(Request {
            id: u32::MAX - 51,
            revision: revision.clone(),
            job: AnalysisJob::CheckEditTemplate {
                source: source.clone(),
                target: target.clone(),
                options: opts.clone(),
                guide: guide.cloned(),
                donor: donor.clone(),
            },
        })
        .await;
    let AnalysisResult::EditTemplateCheck(checked) = response(result, u32::MAX - 51, &revision)?
    else {
        return Err("Unexpected independent donor-check result.".into());
    };
    checked
        .validate_bindings_for(source, target, opts, guide, &donor)
        .map_err(|e| e.to_string())?;
    if !checked.consistent || !report.check.consistent {
        return Err(
            "Fresh independent donor/target check failed. No product can be released.".into(),
        );
    }
    Ok(DonorRequest {
        source: old.source.clone(),
        target: old.target.clone(),
        report: Arc::new(report),
        checked,
        edited: old.edited,
        epoch: old.epoch,
    })
}
fn producer_product(request: &DonorRequest, edited: bool) -> Result<Molecule, String> {
    let mut product = if edited {
        request.report.edited_reference.clone()
    } else {
        request.report.donor.clone()
    };
    let scientific = product.revision();
    let bytes = donor_json(request)?;
    let provenance = String::from_utf8(bytes).map_err(|_| "Cannot encode donor provenance.")?;
    product
        .metadata
        .insert("crispr:scientific_product_revision".into(), scientific);
    product
        .metadata
        .insert("crispr:producer_json".into(), provenance);
    product.metadata.insert("crispr:authority".into(),"Unverified producer evidence on restoration; fresh exact source/guide/independent-target checks are required".into());
    product.metadata.insert(
        "crispr:product_kind".into(),
        if edited {
            "edited reference design"
        } else {
            "donor template design"
        }
        .into(),
    );
    product.metadata.insert(
        "crispr:source_document".into(),
        request.source.binding.document_id.clone(),
    );
    product.metadata.insert(
        "crispr:source_event".into(),
        request.source.binding.event_id.clone(),
    );
    product.metadata.insert(
        "crispr:target_document".into(),
        request.target.binding.document_id.clone(),
    );
    product.metadata.insert(
        "crispr:target_event".into(),
        request.target.binding.event_id.clone(),
    );
    product.metadata.insert(
        "crispr:target_revision".into(),
        request.target.binding.revision.clone(),
    );
    Ok(product)
}

#[component]
pub fn DonorOpenPanel(
    request: DonorRequest,
    snapshot: Callback<(), Result<collagene_core::documents::workspace::Workspace, String>>,
    validate: Callback<collagene_core::documents::workspace::Workspace, Result<(), String>>,
    on_restore: Callback<collagene_core::documents::workspace::Workspace>,
    db: StoredValue<Option<crate::platform::storage::Database>, LocalStorage>,
    busy: RwSignal<bool>,
    loaded: RwSignal<bool>,
    save_generation: RwSignal<u64>,
    epoch: RwSignal<u64>,
    message: RwSignal<String>,
) -> impl IntoView {
    let request = StoredValue::new(request);
    let status = RwSignal::new(String::new());
    let alive = Arc::new(AtomicBool::new(true));
    let worker = StoredValue::new_local(None::<Rc<CachedWorker>>);
    let cleanup = alive.clone();
    on_cleanup(move || {
        cleanup.store(false, Ordering::Relaxed);
        if let Some(w) = worker.get_value() {
            w.close();
        }
    });
    let record = {
        let alive = alive.clone();
        move |_| {
            if busy.get_untracked() {
                return;
            }
            let Some(database) = db.get_value() else {
                status.set("Browser storage is unavailable; no product was opened.".into());
                return;
            };
            let input = request.get_value();
            if input.epoch != epoch.get_untracked() {
                status.set(
                    "The workspace changed. Reopen the current source-bound donor design.".into(),
                );
                return;
            }
            busy.set(true);
            save_generation.update(|v| *v = v.wrapping_add(1));
            let prepared = (|| -> Result<_, String> {
                let frame = snapshot.run(())?;
                let fingerprint = collagene_core::material::lineage::workspace_fingerprint(&frame)
                    .map_err(|e| e.to_string())?;
                let histories =
                    collagene_core::documents::history::DocumentHistory::restore_all(&frame)
                        .map_err(|e| e.to_string())?;
                let exact = |info: &DocumentInfo| -> Result<Molecule, String> {
                    let history = histories
                        .iter()
                        .find(|h| {
                            h.id == info.binding.document_id
                                && h.head_event() == info.binding.event_id
                        })
                        .ok_or(
                            "A source or independent target event changed. Design and check again.",
                        )?;
                    history
                        .revision_molecule(&info.binding.event_id)
                        .map_err(|e| e.to_string())
                };
                let source = exact(&input.source)?;
                let target = exact(&input.target)?;
                let w = Rc::new(CachedWorker::new()?);
                worker.set_value(Some(w.clone()));
                Ok((frame, fingerprint, source, target, w))
            })();
            let (frame, fingerprint, source, target, w) = match prepared {
                Ok(v) => v,
                Err(e) => {
                    busy.set(false);
                    status.set(e);
                    return;
                }
            };
            let alive = alive.clone();
            status.set("Reconstructing and independently checking the exact source/target events before an atomic native update…".into());
            spawn_local(async move {
                let result=async{let fresh=replay_donor(&source,&target,&input,&w).await?;let guard=||->Result<(),String>{if !alive.load(Ordering::Relaxed)||epoch.get_untracked()!=input.epoch{return Err("Dialog closed or workspace changed. Nothing was opened.".into());}let current=snapshot.run(())?;if collagene_core::material::lineage::workspace_fingerprint(&current).map_err(|e|e.to_string())?!=fingerprint{return Err("The workspace changed during donor reconstruction. Nothing was opened.".into());}Ok(())};guard()?;let product=producer_product(&fresh,input.edited)?;let mut next=frame.clone();let histories=collagene_core::documents::history::DocumentHistory::restore_all(&next).map_err(|e|e.to_string())?;let source_history=histories.iter().find(|h|h.id==input.source.binding.document_id).ok_or("Original source history is unavailable.")?;let prefix=format!("crispr-{}-{}",if input.edited{"edited"}else{"donor"},&product.revision()[..24]);let mut id=prefix.clone();let mut suffix=1usize;while next.documents.iter().any(|d|d.id==id){id=format!("{prefix}-{suffix}");suffix+=1;}let derived=collagene_core::documents::history::DocumentHistory::derive_from_events(&[(source_history,input.source.binding.event_id.as_str())],id.clone(),product,format!("Literal CRISPR {} design; source event {}; independently compared target requirement {} / {}; no cellular editing or experimental ancestry is established",if input.edited{"edited reference"}else{"donor template"},input.source.binding.event_id,input.target.binding.document_id,input.target.binding.event_id)).map_err(|e|e.to_string())?;derived.contribute(&mut next).map_err(|e|e.to_string())?;next.active_document=id;next.validate().map_err(|e|e.to_string())?;validate.run(next.clone())?;guard()?;let current=crate::workspace::session::encode_workspace(next.clone(),u32::MAX-52).await?;let previous=crate::workspace::session::encode_workspace(frame,u32::MAX-53).await?;guard()?;crate::platform::storage::replace_workspace(&database,&current,&previous).await?;Ok::<_,String>(next)}.await;
                match result {
                    Ok(frame) => {
                        loaded.set(true);
                        message.set("Opened a freshly checked sequence design with source-only ancestry and noncausal target provenance. Native restoration does not restore current verification.".into());
                        on_restore.run(frame);
                    }
                    Err(e) => {
                        if alive.load(Ordering::Relaxed) {
                            status.set(format!("Donor product was not opened. {e}"));
                        }
                    }
                }
                busy.set(false);
            });
        }
    };
    view! {<section class="crispr-open-panel"><p>{request.with_value(|r|format!("Open checked {}",if r.edited{"edited reference design"}else{"donor template design"}))}</p><p class="helper">"The original source event is the sole design parent. The independent target remains a comparison requirement in producer JSON, not a causal or experimental parent. A fresh worker reconstruction and separate independent donor check run before the product and native recovery slots are updated atomically."</p><button class="button primary" disabled=move||busy.get() on:click=record>"Recheck and open CRISPR product"</button><p role="status">{move||status.get()}</p></section>}
}
