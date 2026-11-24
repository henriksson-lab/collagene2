//! Virtual sequence canvas with revision-bound, cached-worker viewport tracks.
use crate::{
    annotation_layout::AnnotationLayout,
    renderer::{ResizeListener, css_color},
    worker::CachedWorker,
};
use collagene_core::{
    model::complement,
    protocol::{AnalysisJob, AnalysisResult, Request},
    view_tracks::*,
    *,
};
use leptos::{html, prelude::*};
use std::{
    rc::Rc,
    sync::{
        Arc,
        atomic::{AtomicBool, Ordering},
    },
};
use wasm_bindgen::{JsCast, closure::Closure};
use wasm_bindgen_futures::spawn_local;
use web_sys::CanvasRenderingContext2d as Ctx;

type Regions = RwSignal<Option<(String, Vec<(usize, usize)>)>>;
struct Timer {
    id: i32,
    _callback: Closure<dyn FnMut()>,
}
impl Drop for Timer {
    fn drop(&mut self) {
        if let Some(w) = web_sys::window() {
            w.clear_timeout_with_handle(self.id);
        }
    }
}
struct DragTimer {
    id: i32,
    _callback: Closure<dyn FnMut()>,
}
impl Drop for DragTimer {
    fn drop(&mut self) {
        if let Some(window) = web_sys::window() {
            window.clear_interval_with_handle(self.id);
        }
    }
}
#[derive(Clone)]
struct HitBox {
    x: f64,
    y: f64,
    width: f64,
    height: f64,
    item: TrackItem,
}
// The shared selection is half-open. A zero-width range is an insertion caret.
// Pointer gestures retain their original base; keyboard extensions retain their
// original boundary, so reversing direction never moves the anchor.
#[derive(Clone)]
struct EditCursor {
    identity: String,
    revision: String,
    anchor: usize,
    caret: usize,
    range: (usize, usize),
    pointer_anchor: Option<usize>,
}
#[derive(Clone)]
struct BaseDrag {
    identity: String,
    revision: String,
    pointer_id: i32,
    anchor: usize,
}
fn contiguous_location(mut location: &Location) -> bool {
    loop {
        match location {
            Location::Range { .. } => return true,
            Location::Complement(inner) => location = inner,
            _ => return false,
        }
    }
}
fn base_at(
    x: f64,
    y: f64,
    top: f64,
    rh: f64,
    cols: usize,
    len: usize,
    dragging: bool,
) -> Option<usize> {
    if len == 0 || cols == 0 || rh <= 0. {
        return None;
    }
    let absolute = (y + top).max(0.);
    let within = absolute % rh;
    let gx = x - 64.;
    if !dragging && (gx < 0. || gx >= cols as f64 * 8.5 || !(12. ..=54.).contains(&within)) {
        return None;
    }
    let row = (absolute / rh).floor() as usize;
    let col = ((gx / 8.5).floor().max(0.) as usize).min(cols - 1);
    let base = row.saturating_mul(cols).saturating_add(col);
    if dragging {
        Some(base.min(len - 1))
    } else {
        (base < len).then_some(base)
    }
}
fn lane(item: &TrackItem) -> usize {
    match item.kind {
        TrackKind::Feature => 0,
        TrackKind::Primer => 1,
        TrackKind::Enzyme => 2,
        TrackKind::Translation => 3,
        TrackKind::Cds => 4,
        TrackKind::Orf => 5,
    }
}
fn lanes(o: &TrackOptions, layout: Option<&AnnotationLayout>) -> [usize; 6] {
    [
        if o.features {
            layout.map(|l| l.features.lane_count).unwrap_or(0)
        } else {
            0
        },
        if o.primers {
            layout.map(|l| l.primers.lane_count).unwrap_or(0)
        } else {
            0
        },
        if o.enzyme_recognition || o.enzyme_cuts {
            2
        } else {
            0
        },
        if o.translations { 6 } else { 0 },
        if o.known_cds == KnownCdsPolicy::Validated {
            2
        } else {
            0
        },
        if o.orfs { 2 } else { 0 },
    ]
}
fn offsets(counts: [usize; 6]) -> [usize; 6] {
    let mut total = 0;
    std::array::from_fn(|i| {
        let at = total;
        total += counts[i];
        at
    })
}
fn slot(
    item: &TrackItem,
    ordinal: usize,
    counts: [usize; 6],
    layout: Option<&AnnotationLayout>,
) -> Option<usize> {
    let group = lane(item);
    if counts[group] == 0 {
        return None;
    }
    if group == 0 || group == 1 {
        let layout = layout?;
        let feature = item.source_feature?;
        if group == 0 {
            layout.features.placements.get(&feature).copied()
        } else {
            layout.primers.placements.get(&feature).copied()
        }
    } else if group == 2 && item.id.starts_with("restriction:") {
        Some(
            item.id
                .trim_start_matches("restriction:")
                .parse::<usize>()
                .unwrap_or(ordinal)
                % counts[group],
        )
    } else if group == 3 {
        Some(usize::from(item.frame.unwrap_or(0)) + if item.reverse == Some(true) { 3 } else { 0 })
    } else {
        Some(ordinal % counts[group])
    }
}
fn track_color(item: &TrackItem, m: &Molecule) -> String {
    if let Some(i) = item.source_feature {
        if let Some(f) = m.features.get(i) {
            return css_color(f);
        }
    }
    match item.kind {
        TrackKind::Enzyme => "#a66c3c",
        TrackKind::Translation => "#537a91",
        TrackKind::Cds => "#7d66a4",
        TrackKind::Orf => "#458666",
        _ => "#6e8d7e",
    }
    .into()
}

#[component]
pub fn SequenceView(
    molecule: Memo<Molecule>,
    selection: RwSignal<Option<(usize, usize)>>,
    selected_feature: RwSignal<Option<usize>>,
    highlighted_spans: Memo<Vec<(usize, usize)>>,
    search_regions: Regions,
    on_show_features: Callback<()>,
    identity: Memo<String>,
    locked: Memo<bool>,
    on_edit: Callback<(String, usize, usize, String), Result<(), String>>,
    on_history: Callback<bool, Result<(), String>>,
    on_message: Callback<String>,
    restriction_sites: Memo<Option<Arc<Vec<restriction_inspection::InspectionSite>>>>,
    restriction_visible: RwSignal<bool>,
    on_show_restrictions: Callback<()>,
) -> impl IntoView {
    let canvas = NodeRef::<html::Canvas>::new();
    let container = NodeRef::<html::Div>::new();
    let scroll = RwSignal::new(0f64);
    let columns = RwSignal::new(60usize);
    let height = RwSignal::new(500f64);
    let settings = RwSignal::new(TrackOptions::default());
    let show_controls = RwSignal::new(false);
    let feature_rows = RwSignal::new(3usize);
    let primer_rows = RwSignal::new(2usize);
    let annotation_layout =
        Memo::new(move |_| {
            Arc::new(molecule.with(|m| {
                crate::annotation_layout::prepare(m, feature_rows.get(), primer_rows.get())
            }))
        });
    let report = RwSignal::new(None::<TrackReport>);
    let status = RwSignal::new("Loading sequence tracks…".to_string());
    let selected_item = RwSignal::new(None::<TrackItem>);
    let page = RwSignal::new(0usize);
    let alive = Arc::new(AtomicBool::new(true));
    let prepare_alive = alive.clone();
    let query_alive = alive.clone();
    let source = Memo::new(move |_| molecule.with(|m| m.revision()));
    let edit_cursor = StoredValue::new(None::<EditCursor>);
    let base_drag = StoredValue::new(None::<BaseDrag>);
    let drag_active = RwSignal::new(false);
    // Coordinates remain in CSS client pixels; each tick resolves them against
    // the current canvas box, so parent page scrolling cannot move the anchor.
    let drag_pointer = StoredValue::new(None::<(f64, f64)>);
    let drag_timer = StoredValue::new_local(None::<DragTimer>);
    let stop_drag = move || {
        drag_active.set(false);
        drag_timer.set_value(None);
        drag_pointer.set_value(None);
        let previous = base_drag.get_value();
        base_drag.set_value(None);
        if let (Some(previous), Some(canvas)) = (previous, canvas.get_untracked()) {
            let _ = canvas.release_pointer_capture(previous.pointer_id);
        }
    };
    Effect::new(move |_| {
        let id = identity.get();
        let revision = source.get();
        let len = molecule.with_untracked(|m| m.sequence.len());
        if selection
            .get_untracked()
            .is_some_and(|(a, b)| a > b || b > len)
        {
            batch(|| {
                selection.set(None);
                selected_feature.set(None);
                search_regions.set(None);
                selected_item.set(None);
            });
        }
        if base_drag.with_value(|d| {
            d.as_ref()
                .is_some_and(|d| d.identity != id || d.revision != revision)
        }) {
            stop_drag();
        }
        if edit_cursor.with_value(|c| {
            c.as_ref()
                .is_some_and(|c| c.identity != id || c.revision != revision)
        }) {
            edit_cursor.set_value(None);
        }
    });
    let set_cursor = Callback::new(
        move |(anchor, caret, pointer_anchor): (usize, usize, Option<usize>)| {
            let range = (anchor.min(caret), anchor.max(caret));
            edit_cursor.set_value(Some(EditCursor {
                identity: identity.get_untracked(),
                revision: source.get_untracked(),
                anchor,
                caret,
                range,
                pointer_anchor,
            }));
            batch(|| {
                search_regions.set(None);
                selected_item.set(None);
                selected_feature.set(None);
                selection.set(Some(range));
            });
        },
    );
    let current_cursor = move || {
        edit_cursor.get_value().filter(|c| {
            c.identity == identity.get_untracked()
                && c.revision == source.get_untracked()
                && selection.get_untracked() == Some(c.range)
                && selected_feature.get_untracked().is_none()
                && search_regions.get_untracked().is_none()
        })
    };

    let options = Memo::new(move |_| {
        let mut o = settings.get();
        // Restriction geometry belongs to the parent's fresh, source-bound scan.
        // Never run a second starter-enzyme scan through the viewport worker.
        o.enzyme_recognition = false;
        o.enzyme_cuts = false;
        molecule.with(|m| {
            if m.sequence.alphabet() != Alphabet::Dna {
                o.enzyme_recognition = false;
                o.enzyme_cuts = false;
            }
            if m.sequence.alphabet() == Alphabet::Protein {
                o.translations = false;
                o.orfs = false;
                o.known_cds = KnownCdsPolicy::Off;
            }
        });
        o
    });
    let lane_counts = Memo::new(move |_| {
        annotation_layout.with(|layout| {
            let mut counts = lanes(&options.get(), layout.as_ref().as_ref().ok());
            counts[2] = if restriction_visible.get()
                && molecule.with(|m| m.sequence.alphabet() == Alphabet::Dna)
            {
                2
            } else {
                0
            };
            counts
        })
    });
    let row_height = Memo::new(move |_| 66. + lane_counts.get().iter().sum::<usize>() as f64 * 18.);
    let drag_alive = alive.clone();
    Effect::new(move |_| {
        let active = drag_active.get();
        let id = identity.get();
        let revision = source.get();
        let is_locked = locked.get();
        if !active
            || is_locked
            || base_drag.with_value(|d| {
                d.as_ref()
                    .is_none_or(|d| d.identity != id || d.revision != revision)
            })
        {
            if base_drag.with_value(Option::is_some) || drag_timer.with_value(Option::is_some) {
                stop_drag();
            }
            return;
        }
        if drag_timer.with_value(Option::is_some) {
            return;
        }
        let tick_alive = drag_alive.clone();
        let callback = Closure::<dyn FnMut()>::new(move || {
            if !tick_alive.load(Ordering::Relaxed) {
                return;
            }
            let Some(drag) = base_drag.get_value() else {
                drag_active.set(false);
                return;
            };
            if !drag_active.get_untracked()
                || locked.get_untracked()
                || drag.identity != identity.get_untracked()
                || drag.revision != source.get_untracked()
            {
                drag_active.set(false);
                return;
            }
            let (Some(canvas), Some(container), Some((client_x, client_y))) = (
                canvas.get_untracked(),
                container.get_untracked(),
                drag_pointer.get_value(),
            ) else {
                drag_active.set(false);
                return;
            };
            let rect = canvas.get_bounding_client_rect();
            if rect.width() < 1.
                || rect.height() < 1.
                || container.client_height() < 1
                || !canvas.has_pointer_capture(drag.pointer_id)
            {
                drag_active.set(false);
                return;
            }
            let x = client_x - rect.left();
            let y = client_y - rect.top();
            // At most 32 CSS pixels per 32 ms tick (~1000 px/s); no frame loop,
            // source cloning, DNA scanning, or scroll-dependent annotation layout.
            let edge = 24.;
            let speed = if y < edge {
                -((edge - y) / edge).clamp(0., 1.)
            } else if y > rect.height() - edge {
                ((y - (rect.height() - edge)) / edge).clamp(0., 1.)
            } else {
                0.
            };
            if speed == 0. {
                return;
            }
            let before = container.scroll_top();
            container.set_scroll_top(
                (before as f64 + speed * 32.)
                    .max(0.)
                    .round()
                    .min(i32::MAX as f64) as i32,
            );
            let top = container.scroll_top();
            if top == before {
                return;
            }
            scroll.set(top as f64);
            if let Some(base) = base_at(
                x,
                y.clamp(0., (rect.height() - 1.).max(0.)),
                top as f64,
                row_height.get_untracked(),
                columns.get_untracked(),
                molecule.with_untracked(|m| m.sequence.len()),
                true,
            ) {
                let (a, p) = if base >= drag.anchor {
                    (drag.anchor, base + 1)
                } else {
                    (drag.anchor + 1, base)
                };
                set_cursor.run((a, p, Some(drag.anchor)));
            }
        });
        if let Some(window) = web_sys::window() {
            match window.set_interval_with_callback_and_timeout_and_arguments_0(
                callback.as_ref().unchecked_ref(),
                32,
            ) {
                Ok(id) => drag_timer.set_value(Some(DragTimer {
                    id,
                    _callback: callback,
                })),
                Err(_) => stop_drag(),
            }
        } else {
            stop_drag();
        }
    });
    let layout_overflow = Memo::new(move |_| {
        annotation_layout.with(|layout| match layout.as_ref() {
            Ok(layout) => {
                let o = options.get();
                usize::from(o.features) * layout.features.overflow.len()
                    + usize::from(o.primers) * layout.primers.overflow.len()
            }
            Err(_) => 0,
        })
    });
    // Display reflow preserves the visible reference position for this document
    // and unchanged reference symbols. Annotation edits can change row height.
    let anchor_identity = Memo::new(move |_| {
        (
            identity.get(),
            molecule.with(|m| {
                blake3::hash(m.sequence.as_str().as_bytes())
                    .to_hex()
                    .to_string()
            }),
        )
    });
    let geometry = StoredValue::new_local(None::<((String, String), usize, f64)>);
    Effect::new(move |_| {
        let key = anchor_identity.get();
        let cols = columns.get();
        let rh = row_height.get();
        let previous = geometry.get_value();
        geometry.set_value(Some((key.clone(), cols, rh)));
        let Some((old_key, old_cols, old_rh)) = previous else {
            return;
        };
        if old_key != key || (old_cols == cols && old_rh == rh) {
            return;
        }
        let old_top = scroll.get_untracked();
        let base = (old_top / old_rh).floor() as usize * old_cols;
        let fraction = (old_top % old_rh) / old_rh;
        let new_top = ((base / cols) as f64 + fraction) * rh;
        // Render bindings update spacer geometry before this post-render effect.
        // Apply immediately so consecutive reflows use the adjusted viewport.
        if let Some(element) = container.get_untracked() {
            element.set_scroll_top(new_top.round().min(i32::MAX as f64) as i32);
            scroll.set(element.scroll_top() as f64);
        }
    });
    // Quantized overlapping tiles avoid preparing a chromosome or transferring it per scroll tick.
    let window = Memo::new(move |_| {
        let cols = columns.get();
        let row = (scroll.get() / row_height.get()).floor() as usize;
        let start = row.saturating_sub(6) / 12 * 12 * cols;
        let rows = (height.get() / row_height.get()).ceil() as usize + 24;
        let n = molecule.with(|m| m.sequence.len());
        TrackWindow {
            start: start.min(n.saturating_sub(1)),
            end: (start + rows * cols).min(n),
        }
    });
    let shared_restrictions = Memo::new(move |_| {
        if !restriction_visible.get() || molecule.with(|m| m.sequence.alphabet() != Alphabet::Dna) {
            return Arc::new(Vec::<TrackItem>::new());
        }
        let window = window.get();
        let Some(sites) = restriction_sites.get() else {
            return Arc::new(Vec::new());
        };
        let mut items = Vec::new();
        for site in sites.iter().filter(|site| {
            site.spans
                .iter()
                .any(|span| span.start < window.end && span.end > window.start)
                || [site.top_cut, site.bottom_cut]
                    .into_iter()
                    .flatten()
                    .any(|cut| cut >= window.start && cut <= window.end)
        }) {
            let label = if site.diagnostics.is_empty() {
                site.enzyme.clone()
            } else {
                format!(
                    "{} · raw cuts top {} / bottom {} · {:?}",
                    site.enzyme, site.top_boundary, site.bottom_boundary, site.diagnostics
                )
            };
            items.push(TrackItem {
                id: format!("restriction:{}", site.index),
                kind: TrackKind::Enzyme,
                label,
                parts: site
                    .spans
                    .iter()
                    .map(|span| ReferencePart {
                        start: span.start,
                        end: span.end,
                    })
                    .collect(),
                reverse: Some(site.reverse),
                source_feature: None,
                enzyme: Some(site.enzyme.clone()),
                certain: site.certain,
                cuts: Some(CutPositions {
                    top: site.top_cut,
                    bottom: site.bottom_cut,
                }),
                frame: None,
                protein: None,
                residues: vec![],
            });
        }
        Arc::new(items)
    });
    let client = StoredValue::new_local(None::<Rc<CachedWorker>>);
    let source_generation = RwSignal::new(0u32);
    let retry = RwSignal::new(0u32);
    let prepared = RwSignal::new(false);
    let in_flight = RwSignal::new(false);
    let attempted = StoredValue::new_local(None::<(u32, TrackWindow, TrackOptions)>);
    let timer = StoredValue::new_local(None::<Timer>);
    let resize = StoredValue::new_local(None::<ResizeListener>);
    let hit_boxes = StoredValue::new_local(Vec::<HitBox>::new());
    Effect::new(move |_| {
        let rev = source.get();
        let _ = retry.get();
        let id = source_generation.get_untracked().wrapping_add(1);
        source_generation.set(id);
        prepared.set(false);
        in_flight.set(false);
        report.set(None);
        selected_item.set(None);
        page.set(0);
        timer.set_value(None);
        status.set("Loading sequence tracks…".into());
        if let Some(previous) = client.get_value() {
            previous.close();
        }
        let worker = match CachedWorker::new() {
            Ok(w) => Rc::new(w),
            Err(e) => {
                client.set_value(None);
                status.set(e);
                return;
            }
        };
        client.set_value(Some(worker.clone()));
        let mol = molecule.get_untracked();
        let alive = prepare_alive.clone();
        spawn_local(async move {
            let result = worker
                .ask(Request {
                    id,
                    revision: rev.clone(),
                    job: AnalysisJob::PrepareTracks { molecule: mol },
                })
                .await;
            if !alive.load(Ordering::Relaxed)
                || source_generation.get_untracked() != id
                || source.get_untracked() != rev
            {
                return;
            }
            match result {
                Ok(r) => match r.result {
                    AnalysisResult::TracksPrepared(value) if value == rev => prepared.set(true),
                    AnalysisResult::Error(e) => status.set(e),
                    _ => status.set("Unexpected track preparation response".into()),
                },
                Err(e) => status.set(e),
            }
        });
    });
    Effect::new(move |_| {
        let ready = prepared.get();
        let busy = in_flight.get();
        let desired = window.get();
        let opts = options.get();
        let rev = source.get();
        let generation = source_generation.get();
        // A prior tile remains valid while its replacement loads. Settings/source
        // changes invalidate scientific meaning immediately, even during a job.
        if report.with_untracked(|r| {
            r.as_ref()
                .is_some_and(|r| r.source_revision != rev || r.options != opts)
        }) {
            report.set(None);
            selected_item.set(None);
        }
        timer.set_value(None);
        if !ready || busy {
            return;
        }
        if attempted.with_value(|last| {
            last.as_ref()
                .is_some_and(|(g, w, o)| *g == generation && *w == desired && *o == opts)
        }) {
            return;
        }
        if report.with_untracked(|r| {
            r.as_ref().is_some_and(|r| {
                r.source_revision == rev && r.window == desired && r.options == opts
            })
        }) {
            return;
        }
        status.set("Updating visible tracks…".into());
        page.set(0);
        let alive = query_alive.clone();
        let callback = Closure::<dyn FnMut()>::new(move || {
            let Some(worker) = client.get_value() else {
                return;
            };
            let opts = opts.clone();
            let rev = rev.clone();
            let alive = alive.clone();
            attempted.set_value(Some((generation, desired, opts.clone())));
            in_flight.set(true);
            spawn_local(async move {
                let result = worker
                    .ask(Request {
                        id: generation,
                        revision: rev.clone(),
                        job: AnalysisJob::CachedViewTracks {
                            window: desired,
                            options: opts.clone(),
                        },
                    })
                    .await;
                if !alive.load(Ordering::Relaxed)
                    || source_generation.get_untracked() != generation
                    || source.get_untracked() != rev
                {
                    return;
                }
                if window.get_untracked() == desired && options.get_untracked() == opts {
                    match result {
                        Ok(response) => match response.result {
                            AnalysisResult::ViewTracks(value)
                                if value.source_revision == rev
                                    && value.window == desired
                                    && value.options == opts =>
                            {
                                status.set(if value.complete {
                                    "Visible tracks complete".into()
                                } else {
                                    "Visible tracks incomplete; inspect warnings".into()
                                });
                                report.set(Some(value));
                            }
                            AnalysisResult::Error(e) => status.set(e),
                            _ => status.set("Unexpected viewport track response".into()),
                        },
                        Err(e) => status.set(e),
                    }
                }
                in_flight.set(false);
            });
        });
        if let Some(w) = web_sys::window() {
            if let Ok(id) = w.set_timeout_with_callback_and_timeout_and_arguments_0(
                callback.as_ref().unchecked_ref(),
                140,
            ) {
                timer.set_value(Some(Timer {
                    id,
                    _callback: callback,
                }));
            }
        }
    });
    // Stop the retained worker as soon as this canvas is unmounted.
    on_cleanup(move || {
        stop_drag();
        alive.store(false, Ordering::Relaxed);
        timer.set_value(None);
        if let Some(previous) = client.get_value() {
            previous.close();
        }
        client.set_value(None);
    });
    Effect::new(move |_| {
        if let Some(element) = container.get() {
            if resize.with_value(|v| v.is_none()) {
                if let Ok(listener) = ResizeListener::new(element.unchecked_ref(), move || {
                    if let Some(v) = scroll.try_get_untracked() {
                        let _ = scroll.try_set(v);
                    }
                }) {
                    resize.set_value(Some(listener));
                }
            }
        }
    });
    let choose = Callback::new(move |item: TrackItem| {
        let parts = item
            .parts
            .iter()
            .map(|p| (p.start, p.end))
            .collect::<Vec<_>>();
        batch(|| {
            selected_feature.set(item.source_feature);
            selection.set(if parts.len() == 1 || item.source_feature.is_some() {
                parts.first().copied()
            } else {
                None
            });
            search_regions.set(Some((source.get_untracked(), parts)));
            selected_item.set(Some(item));
        });
    });
    let press_base = Callback::new(move |(x, y, shift): (f64, f64, bool)| {
        let len = molecule.with_untracked(|m| m.sequence.len());
        let Some(base) = base_at(
            x,
            y,
            scroll.get_untracked(),
            row_height.get_untracked(),
            columns.get_untracked(),
            len,
            false,
        ) else {
            return;
        };
        let anchor = if shift {
            current_cursor()
                .and_then(|c| {
                    c.pointer_anchor
                        .or(Some(c.anchor.min(len.saturating_sub(1))))
                })
                .or_else(|| {
                    selection
                        .get_untracked()
                        .map(|(a, _)| a.min(len.saturating_sub(1)))
                })
                .unwrap_or(base)
        } else {
            base
        };
        let (a, p) = if base >= anchor {
            (anchor, base + 1)
        } else {
            (anchor + 1, base)
        };
        set_cursor.run((a, p, Some(anchor)));
    });
    let keydown = move |e: web_sys::KeyboardEvent| {
        if e.is_composing() || e.alt_key() {
            return;
        }
        let key = e.key();
        if e.ctrl_key() || e.meta_key() {
            let shortcut = key.to_ascii_lowercase();
            let select_all = shortcut == "a" && !e.shift_key();
            let history = if shortcut == "z" {
                Some(e.shift_key())
            } else if shortcut == "y" && e.ctrl_key() && !e.shift_key() {
                Some(true)
            } else {
                None
            };
            if !select_all && history.is_none() {
                return;
            } // Leave C/X/V and other shortcuts untouched.
            e.prevent_default();
            if locked.get_untracked() {
                return;
            }
            stop_drag();
            if select_all {
                let len = molecule.with_untracked(|m| m.sequence.len());
                set_cursor.run((0, len, None));
            } else if let Some(redo) = history {
                match on_history.run(redo) {
                    Ok(()) => edit_cursor.set_value(None),
                    Err(error) => on_message.run(error),
                }
            }
            return;
        }
        if key == "Escape" {
            e.prevent_default();
            stop_drag();
            edit_cursor.set_value(None);
            batch(|| {
                selection.set(None);
                selected_feature.set(None);
                selected_item.set(None);
                search_regions.set(None);
            });
            return;
        }
        let len = molecule.with_untracked(|m| m.sequence.len());
        let range = selection
            .get_untracked()
            .filter(|(a, b)| a <= b && *b <= len);
        let cursor = current_cursor();
        let (mut anchor, mut caret) = cursor
            .as_ref()
            .map(|c| (c.anchor, c.caret))
            .unwrap_or_else(|| range.unwrap_or((0, 0)));
        if e.shift_key() {
            if let Some(c) = cursor
                .as_ref()
                .filter(|c| c.pointer_anchor.is_some() && c.range.1 - c.range.0 == 1)
            {
                let base = c.pointer_anchor.unwrap();
                if matches!(key.as_str(), "ArrowLeft" | "ArrowUp" | "Home") {
                    anchor = base + 1;
                    caret = base;
                } else if matches!(key.as_str(), "ArrowRight" | "ArrowDown" | "End") {
                    anchor = base;
                    caret = base + 1;
                }
            }
        }
        let next = match key.as_str() {
            "ArrowRight" => Some(if !e.shift_key() && range.is_some_and(|(a, b)| a < b) {
                range.unwrap().1
            } else {
                caret.saturating_add(1).min(len)
            }),
            "ArrowLeft" => Some(if !e.shift_key() && range.is_some_and(|(a, b)| a < b) {
                range.unwrap().0
            } else {
                caret.saturating_sub(1)
            }),
            "ArrowDown" => Some(caret.saturating_add(columns.get_untracked()).min(len)),
            "ArrowUp" => Some(caret.saturating_sub(columns.get_untracked())),
            "Home" => Some(0),
            "End" => Some(len),
            _ => None,
        };
        if let Some(next) = next {
            e.prevent_default();
            stop_drag();
            set_cursor.run((if e.shift_key() { anchor } else { next }, next, None));
            return;
        }
        let replacement = if key == "Backspace" || key == "Delete" {
            Some(String::new())
        } else if key.chars().count() == 1 {
            Some(key.clone())
        } else {
            None
        };
        let Some(replacement) = replacement else {
            return;
        };
        e.prevent_default();
        if locked.get_untracked() {
            on_message.run("The workspace is busy. Wait before editing the sequence.".into());
            return;
        }
        if key == " " {
            on_message.run("Type sequence symbols without spaces.".into());
            return;
        }
        let compound = selected_feature.get_untracked().is_some_and(|i| {
            molecule.with_untracked(|m| {
                m.features
                    .get(i)
                    .is_some_and(|f| !contiguous_location(&f.location))
            })
        }) || search_regions
            .with_untracked(|r| r.as_ref().is_some_and(|(_, parts)| parts.len() != 1));
        if compound {
            on_message.run("This annotation has multiple or unresolved parts. Select contiguous reference bases before typing.".into());
            return;
        }
        let Some((mut a, mut b)) = range else {
            on_message.run("Click a base or place the caret before editing.".into());
            return;
        };
        if a == b {
            if key == "Backspace" {
                if a == 0 {
                    return;
                }
                a -= 1;
            } else if key == "Delete" {
                if b == len {
                    return;
                }
                b += 1;
            }
        }
        let revision = source.get_untracked();
        match on_edit.run((revision, a, b, replacement.clone())) {
            Ok(()) => {
                stop_drag();
                set_cursor.run((a + replacement.len(), a + replacement.len(), None));
                on_message.run(String::new());
            }
            Err(error) => on_message.run(error),
        }
    };
    Effect::new(move |_| {
        molecule.with(|mol| {
            let selected = highlighted_spans.get();
            let caret = selection.get().filter(|(a, b)| a == b).map(|(a, _)| a);
            let top = scroll.get();
            let rh = row_height.get();
            let opts = options.get();
            let revision = source.get();
            let current = report
                .get()
                .filter(|r| r.source_revision == revision && r.options == opts);
            let packed = annotation_layout.get();
            let layout = packed.as_ref().as_ref().ok();
            let Some(canvas) = canvas.get() else {
                return;
            };
            let Some(container) = container.get() else {
                return;
            };
            // DOM client sizes are rounded. A rounded-up canvas can create a
            // scrollbar, changing the very viewport being measured. Preserve
            // the fractional usable CSS box and reserve whole backing pixels.
            let rect = container.get_bounding_client_rect();
            let width = (rect.width()
                - f64::from(container.offset_width() - container.client_width()))
            .max(1.);
            let h = (rect.height()
                - f64::from(container.offset_height() - container.client_height()))
            .max(1.);
            if height.get_untracked() != h {
                height.set(h);
            }
            let cols =
                ((((container.client_width().max(200) as f64) - 75.) / 8.5).floor() as usize / 10
                    * 10)
                    .clamp(10, 100);
            if columns.get_untracked() != cols {
                columns.set(cols);
            }
            let dpr = web_sys::window()
                .map(|w| w.device_pixel_ratio())
                .unwrap_or(1.)
                .min(2.);
            let pw = (width * dpr).ceil() as u32;
            let ph = (h * dpr).ceil() as u32;
            if canvas.width() != pw {
                canvas.set_width(pw);
            }
            if canvas.height() != ph {
                canvas.set_height(ph);
            }
            let style = canvas.unchecked_ref::<web_sys::HtmlElement>().style();
            let _ = style.set_property("height", &format!("{h}px"));
            let _ = style.set_property("width", &format!("{width}px"));
            let Some(ctx) = canvas
                .get_context("2d")
                .ok()
                .flatten()
                .and_then(|c| c.dyn_into::<Ctx>().ok())
            else {
                return;
            };
            let _ = ctx.set_transform(f64::from(pw) / width, 0., 0., f64::from(ph) / h, 0., 0.);
            ctx.set_fill_style_str("#ffffff");
            ctx.fill_rect(0., 0., width, h);
            let first = (top / rh).floor() as usize;
            let count = (h / rh).ceil() as usize + 1;
            let start_x = 64.;
            let cell = 8.5;
            let counts = lane_counts.get();
            let offset = offsets(counts);
            let shared_items = shared_restrictions.get();
            let mut boxes = vec![];
            let mut painted = 0usize;
            for row in first..first + count {
                let start = row * cols;
                if start >= mol.sequence.len() {
                    break;
                }
                let end = (start + cols).min(mol.sequence.len());
                let y = row as f64 * rh - top + 30.;
                ctx.set_font("14px monospace");
                ctx.set_fill_style_str("#89938e");
                ctx.set_text_align("right");
                let _ = ctx.fill_text(&(start + 1).to_string(), 48., y);
                ctx.set_text_align("left");
                for &(a, b) in &selected {
                    let left = a.max(start);
                    let right = b.min(end);
                    if left < right {
                        ctx.set_fill_style_str("#d5eee5");
                        ctx.fill_rect(
                            start_x + (left - start) as f64 * cell,
                            y - 15.,
                            (right - left) as f64 * cell,
                            39.,
                        );
                    }
                }
                for (i, b) in mol.sequence.as_str().as_bytes()[start..end]
                    .iter()
                    .enumerate()
                {
                    ctx.set_fill_style_str(match b {
                        b'A' => "#337c6d",
                        b'C' => "#4476a3",
                        b'G' => "#b18b38",
                        b'T' | b'U' => "#ad6b6e",
                        _ => "#5d6a64",
                    });
                    let _ = ctx.fill_text(&(*b as char).to_string(), start_x + i as f64 * cell, y);
                    if mol.sequence.alphabet() != Alphabet::Protein {
                        ctx.set_fill_style_str("#98a39d");
                        let _ = ctx.fill_text(
                            &(complement(*b, mol.sequence.alphabet()) as char).to_string(),
                            start_x + i as f64 * cell,
                            y + 19.,
                        );
                    }
                }
                if let Some(p) = caret.filter(|p| *p <= mol.sequence.len()) {
                    let caret_row = if p > 0 && p == mol.sequence.len() && p % cols == 0 {
                        (p - 1) / cols
                    } else {
                        p / cols
                    };
                    if caret_row == row {
                        let x = start_x + (p - row * cols) as f64 * cell;
                        ctx.set_stroke_style_str("#195d4f");
                        ctx.set_line_width(2.);
                        ctx.begin_path();
                        ctx.move_to(x, y - 15.);
                        ctx.line_to(x, y + 23.);
                        ctx.stroke();
                    }
                }
                {
                    let mut ordinals = [0usize; 6];
                    // Keep the scientific viewport report unchanged. External RE
                    // evidence is painted from the shared parent scan separately.
                    for item in current
                        .as_ref()
                        .into_iter()
                        .flat_map(|r| r.items.iter())
                        .chain(shared_items.iter())
                    {
                        let group = lane(item);
                        let ordinal = ordinals[group];
                        ordinals[group] += 1;
                        let Some(item_slot) = slot(item, ordinal, counts, layout) else {
                            continue;
                        };
                        let lane_y = y + 36. + (offset[group] + item_slot) as f64 * 18.;
                        let color = track_color(item, mol);
                        for part in &item.parts {
                            let left = part.start.max(start);
                            let right = part.end.min(end);
                            if left >= right {
                                continue;
                            }
                            let x = start_x + (left - start) as f64 * cell;
                            let size = (right - left) as f64 * cell;
                            ctx.set_fill_style_str(&color);
                            ctx.fill_rect(x, lane_y, size, 3.);
                            if item.reverse == Some(true) {
                                ctx.begin_path();
                                ctx.move_to(x + 4., lane_y - 3.);
                                ctx.line_to(x, lane_y + 1.);
                                ctx.line_to(x + 4., lane_y + 5.);
                                ctx.stroke();
                            }
                            if item.residues.is_empty() {
                                ctx.set_font("10px system-ui");
                                if size > 25. {
                                    let text = format!(
                                        "{}{}",
                                        if item.certain { "" } else { "? " },
                                        item.label
                                    );
                                    let label =
                                        text.chars().take((size / 6.) as usize).collect::<String>();
                                    let _ = ctx.fill_text(&label, x, lane_y + 13.);
                                }
                            }
                            boxes.push(HitBox {
                                x,
                                y: lane_y - 5.,
                                width: size,
                                height: 18.,
                                item: item.clone(),
                            });
                            painted += 1;
                        }
                        ctx.set_font("12px monospace");
                        for residue in &item.residues {
                            if let Some(anchor) = residue.codon_parts.first() {
                                let pos = if item.reverse == Some(true) {
                                    anchor.end - 1
                                } else {
                                    anchor.start
                                };
                                if pos >= start && pos < end {
                                    ctx.set_fill_style_str(&color);
                                    let _ = ctx.fill_text(
                                        &residue.amino_acid,
                                        start_x + (pos - start) as f64 * cell,
                                        lane_y + 13.,
                                    );
                                    let codon_item = TrackItem {
                                        id: format!("{}:{}", item.id, residue.residue_index),
                                        kind: item.kind,
                                        label: format!(
                                            "{} · residue {} ({})",
                                            item.label,
                                            residue.residue_index + 1,
                                            residue.amino_acid
                                        ),
                                        parts: residue.codon_parts.clone(),
                                        reverse: item.reverse,
                                        source_feature: None,
                                        enzyme: None,
                                        certain: item.certain,
                                        cuts: None,
                                        frame: item.frame,
                                        protein: Some(residue.amino_acid.clone()),
                                        residues: vec![residue.clone()],
                                    };
                                    boxes.push(HitBox {
                                        x: start_x + (pos - start) as f64 * cell,
                                        y: lane_y + 3.,
                                        width: cell
                                            * if item.reverse == Some(true) { 1. } else { 3. },
                                        height: 14.,
                                        item: codon_item,
                                    });
                                }
                            }
                        }
                        if opts.enzyme_cuts || item.id.starts_with("restriction:") {
                            if let Some(cuts) = &item.cuts {
                                if item.id.starts_with("restriction:") {
                                    // Show staggered cleavage against the actual two DNA rows,
                                    // separately from the labeled recognition track below them.
                                    if let (Some(top), Some(bottom)) = (cuts.top, cuts.bottom) {
                                        if top >= start
                                            && top <= end
                                            && bottom >= start
                                            && bottom <= end
                                        {
                                            ctx.set_stroke_style_str("#8b5ca4");
                                            ctx.set_line_width(0.8);
                                            ctx.begin_path();
                                            ctx.move_to(
                                                start_x + (top - start) as f64 * cell,
                                                y + 5.,
                                            );
                                            ctx.line_to(
                                                start_x + (bottom - start) as f64 * cell,
                                                y + 5.,
                                            );
                                            ctx.stroke();
                                        }
                                    }
                                }
                                for (pos, bottom) in [(cuts.top, false), (cuts.bottom, true)] {
                                    if let Some(pos) = pos {
                                        if pos >= start && pos <= end {
                                            let x = start_x + (pos - start) as f64 * cell;
                                            ctx.set_stroke_style_str(if bottom {
                                                "#b45c19"
                                            } else {
                                                "#2563a6"
                                            });
                                            ctx.set_line_width(1.2);
                                            if item.id.starts_with("restriction:") {
                                                ctx.begin_path();
                                                ctx.move_to(
                                                    x,
                                                    if bottom { y + 5. } else { y - 15. },
                                                );
                                                ctx.line_to(
                                                    x,
                                                    if bottom { y + 23. } else { y + 5. },
                                                );
                                                ctx.stroke();
                                            }
                                            ctx.begin_path();
                                            ctx.move_to(x, lane_y - 5.);
                                            ctx.line_to(x, lane_y + if bottom { 14. } else { 7. });
                                            ctx.stroke();
                                            if item.id.starts_with("restriction:") {
                                                boxes.push(HitBox {
                                                    x: x - 4.,
                                                    y: lane_y - 5.,
                                                    width: 8.,
                                                    height: 20.,
                                                    item: item.clone(),
                                                });
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            hit_boxes.set_value(boxes);
            let _ = canvas.set_attribute("data-row-height", &rh.to_string());
            let _ = canvas.set_attribute("data-feature-lanes", &counts[0].to_string());
            let _ = canvas.set_attribute("data-primer-lanes", &counts[1].to_string());
            let _ =
                canvas.set_attribute("data-layout-overflow", &layout_overflow.get().to_string());
            let _ = canvas.set_attribute("data-layout-source", &revision);

            let _ = canvas.set_attribute("data-painted-tracks", &painted.to_string());
            let _ = canvas.set_attribute(
                "data-shared-restriction-count",
                &shared_items.len().to_string(),
            );
            let _ = canvas.set_attribute(
                "data-shared-restriction-ready",
                if restriction_sites.with(|sites| sites.is_some()) {
                    "true"
                } else {
                    "false"
                },
            );
            let _ = canvas.set_attribute(
                "data-track-item-count",
                &current
                    .as_ref()
                    .map(|r| r.items.len())
                    .unwrap_or(0)
                    .to_string(),
            );
            let _ = canvas.set_attribute(
                "data-track-complete",
                if current.as_ref().is_some_and(|r| r.complete) {
                    "true"
                } else {
                    "false"
                },
            );
        });
    });
    Effect::new(move |_| {
        let ranges = highlighted_spans.get();
        let _ = selection.get();
        if base_drag.get_value().is_some() {
            return;
        }
        let active = current_cursor().map(|c| {
            if c.caret > c.anchor {
                c.caret.saturating_sub(1)
            } else {
                c.caret
            }
        });
        if let Some(a) = active.or_else(|| ranges.first().map(|(a, _)| *a)) {
            if let Some(c) = container.get() {
                let rh = row_height.get_untracked();
                let y = (a / columns.get_untracked()) as f64 * rh;
                let top = c.scroll_top() as f64;
                if y < top || y + rh > top + c.client_height() as f64 {
                    c.set_scroll_top((y - rh).max(0.) as i32);
                    scroll.set(c.scroll_top() as f64);
                }
            }
        }
    });
    view! {
        <div class="sequence-editor-hint">"Drag bases to select · type reference symbols to edit · Shift + arrows to extend"</div>
        <div class="track-toolbar"><button class="text-tool" aria-expanded=move||show_controls.get().to_string() on:click=move|_|show_controls.update(|v|*v=!*v)>"Sequence tracks"</button><span class="track-status" role="status">{move||status.get()}</span><button class="text-tool" title="Reload sequence tracks" on:click=move|_|retry.update(|n|*n=n.wrapping_add(1))>"↻"</button></div>
        <Show when=move||restriction_visible.get()&&molecule.with(|m|m.sequence.alphabet()==Alphabet::Dna)><div class="shared-restriction-status" role="status">{move||restriction_sites.with(|sites|sites.as_ref().map(|sites|format!("Restriction tracks: {} selected sites from the shared source scan.",sites.len())).unwrap_or_else(||"Restriction tracks: shared source scan pending or unavailable.".into()))}</div></Show>
        <Show when=move||{layout_overflow.get()>0 || (annotation_layout.with(|p|p.as_ref().is_err()) && options.with(|o|o.features||o.primers))}>
            <div class="annotation-layout-note" role="status"><span>{move||if annotation_layout.with(|p|p.as_ref().is_err()) {"Annotation display is limited. Inspect source annotations in Features.".to_string()}else{format!("{} additional annotations are available in Features.",layout_overflow.get())}}</span><button class="text-tool" on:click=move|_|on_show_features.run(())>"Inspect annotations"</button></div>
        </Show>
        <Show when=move||show_controls.get()><div class="track-controls" aria-label="Sequence track settings">
            <div class="annotation-row-controls"><label>"Annotation rows"<select aria-label="Annotation rows" prop:value=move||feature_rows.get().to_string() on:change=move|e|feature_rows.set(event_target_value(&e).parse::<usize>().ok().filter(|v|[1,3,6,8].contains(v)).unwrap_or(3))>{[1,3,6,8].map(|n|view!{<option value=n.to_string()>{n.to_string()}</option>})}</select></label><label>"Primer rows"<select aria-label="Primer rows" prop:value=move||primer_rows.get().to_string() on:change=move|e|primer_rows.set(event_target_value(&e).parse::<usize>().ok().filter(|v|[1,2,4,8].contains(v)).unwrap_or(2))>{[1,2,4,8].map(|n|view!{<option value=n.to_string()>{n.to_string()}</option>})}</select></label></div>
            <div class="track-toggles">
                <label><input type="checkbox" aria-label="Feature tracks" prop:checked=move||settings.get().features on:change=move|e|settings.update(|o|o.features=event_target_checked(&e))/><span>"Features"</span></label>
                <label><input type="checkbox" aria-label="Primer tracks" prop:checked=move||settings.get().primers on:change=move|e|settings.update(|o|o.primers=event_target_checked(&e))/><span>"Primers"</span></label>
                <label><input type="checkbox" aria-label="Restriction site tracks" prop:checked=move||restriction_visible.get() disabled=move||molecule.with(|m|m.sequence.alphabet()!=Alphabet::Dna) on:change=move|e|restriction_visible.set(event_target_checked(&e))/><span>"Restriction sites and cuts"</span></label>
                <label><input type="checkbox" aria-label="Translation tracks" prop:checked=move||options.get().translations disabled=move||molecule.with(|m|m.sequence.alphabet()==Alphabet::Protein) on:change=move|e|settings.update(|o|o.translations=event_target_checked(&e))/><span>"Six-frame translation"</span></label>
                <label><input type="checkbox" aria-label="CDS tracks" prop:checked=move||options.get().known_cds==KnownCdsPolicy::Validated disabled=move||molecule.with(|m|m.sequence.alphabet()==Alphabet::Protein) on:change=move|e|settings.update(|o|o.known_cds=if event_target_checked(&e){KnownCdsPolicy::Validated}else{KnownCdsPolicy::Off})/><span>"Reviewed CDS"</span></label>
                <label><input type="checkbox" aria-label="ORF tracks" prop:checked=move||options.get().orfs disabled=move||molecule.with(|m|m.sequence.alphabet()==Alphabet::Protein) on:change=move|e|settings.update(|o|o.orfs=event_target_checked(&e))/><span>"Predicted ORFs"</span></label>
            </div>
            <button class="text-tool" on:click=move|_|on_show_restrictions.run(())>"Choose restriction enzymes"</button>
            <p class="helper">"Restriction sites and cut boundaries come from the shared source scan. Viewport report downloads contain the sequence worker tracks only."</p>
            <Show when=move||restriction_visible.get()&&restriction_sites.with(|sites|sites.is_none())><p class="helper" role="status">"Shared restriction sites are not available yet. Choose restriction enzymes to inspect the source scan."</p></Show>
            <label>"Track direction"<select aria-label="Track direction" prop:value=move||format!("{:?}",settings.get().orientation) on:change=move|e|settings.update(|o|o.orientation=match event_target_value(&e).as_str(){"Forward"=>TrackOrientation::Forward,"Reverse"=>TrackOrientation::Reverse,_=>TrackOrientation::Both})><option value="Both">"Both strands"</option><option value="Forward">"Forward"</option><option value="Reverse">"Reverse"</option></select></label>
            <div class="track-toggles">{(0u8..3).map(|frame|view!{<label><input type="checkbox" aria-label=format!("Track frame {}",frame+1) prop:checked=move||settings.with(|o|o.frames.contains(&frame)) on:change=move|e|settings.update(|o|{o.frames.retain(|f|*f!=frame);if event_target_checked(&e){o.frames.push(frame);o.frames.sort_unstable();}})/><span>{format!("Frame {}",frame+1)}</span></label>}).collect_view()}</div>
            <label>"Minimum ORF length (aa)"<input type="number" min="0" max="100000" aria-label="Track minimum ORF length" prop:value=move||settings.get().min_orf_amino_acids.to_string() on:change=move|e|settings.update(|o|o.min_orf_amino_acids=event_target_value(&e).parse::<usize>().unwrap_or(20).min(100000))/></label>
            <label>"Track genetic code"<select aria-label="Track genetic code" on:change=move|e|settings.update(|o|o.genetic_code=event_target_value(&e).parse().unwrap_or(1))>{collagene_core::translation::genetic_codes().iter().map(|c|view!{<option value=c.id.to_string()>{format!("{} · {}",c.id,c.name)}</option>}).collect_view()}</select></label>
            <p class="helper">"Nonoverlapping annotations share rows. Joined annotations stay together; additional annotations remain in Features and the viewport list. Predicted ORFs and ambiguous sites are provisional. Blue/red ticks are reference top/bottom cut boundaries. Tracks are bounded and refresh after scrolling."</p>
            <details class="viewport-tracks"><summary>"Viewport items and warnings"</summary>{move||report.get().map(|r|{let count=r.items.len();view!{<p class="helper">{format!("{} items · {} work units · {}",count,r.work.used_units,if r.complete{"complete"}else{"incomplete"})}</p>{r.warnings.into_iter().map(|w|view!{<p class="helper">{w}</p>}).collect_view()}{r.items.into_iter().skip(page.get()*50).take(50).map(|item|{let text=format!("{:?} · {} · {}{}",item.kind,item.label,item.parts.iter().map(|p|format!("{}–{}",p.start+1,p.end)).collect::<Vec<_>>().join(", "),if item.certain{""}else{" · possible/review"});view!{<button class="result-entry viewport-item" data-kind=format!("{:?}",item.kind) on:click=move|_|choose.run(item.clone())>{text}</button>}}).collect_view()}<Show when=move||{count>50}><button class="text-tool" disabled=move||page.get()==0 on:click=move|_|page.update(|p|*p=p.saturating_sub(1))>"Previous track items"</button><button class="text-tool" disabled=move||{(page.get()+1)*50>=count} on:click=move|_|page.update(|p|*p+=1)>"Next track items"</button></Show>}})}<button class="text-tool" disabled=move||report.get().is_none() on:click=move|_|if let Some(r)=report.get_untracked(){match serde_json::to_vec_pretty(&r){Ok(bytes)=>{let _=crate::app::download("viewport-tracks.json","application/json",&bytes);},Err(e)=>status.set(e.to_string())}}>"Download viewport report"</button></details>
        </div></Show>
        {move||selected_item.get().map(|item|view!{<div class="track-selection"><strong>{item.label}</strong><span>{format!(" · {:?}{}",item.kind,if item.certain{""}else{" · possible/review"})}</span><span>{item.cuts.as_ref().map(|c|format!(" · top {}, bottom {}",crate::search_view::boundary_label(c.top),crate::search_view::boundary_label(c.bottom)))}</span></div>})}
        <div class="sequence-scroll" node_ref=container on:scroll=move|_|if let Some(c)=container.get(){scroll.set(c.scroll_top()as f64);}>
            <div class="sequence-spacer" style:height=move||format!("{}px",molecule.with(|m|m.sequence.len()).div_ceil(columns.get())as f64*row_height.get())>
                <canvas node_ref=canvas class="sequence-canvas" role="textbox" aria-multiline="true" aria-readonly=move||locked.get().to_string() tabindex="0" aria-label="Sequence editor. Drag reference or complement bases to select. Shift extends selection; type to replace reference symbols. Arrow keys move the caret; Shift and arrows extend. Home and End move to document boundaries." on:keydown=keydown
                    on:pointerdown=move|e:web_sys::PointerEvent|{
                        if e.button()!=0 || locked.get_untracked(){return;}
                        if let Some(c)=canvas.get(){
                            let _=c.focus();let rect=c.get_bounding_client_rect();let x=e.client_x()as f64-rect.left();let y=e.client_y()as f64-rect.top();
                            let hit=hit_boxes.with_value(|boxes|boxes.iter().rev().find(|b|x>=b.x&&x<b.x+b.width&&y>=b.y&&y<b.y+b.height).map(|b|b.item.clone()));
                            if let Some(hit)=hit{e.prevent_default();stop_drag();edit_cursor.set_value(None);choose.run(hit);return;}
                            if base_at(x,y,scroll.get_untracked(),row_height.get_untracked(),columns.get_untracked(),molecule.with_untracked(|m|m.sequence.len()),false).is_some(){
                                e.prevent_default();press_base.run((x,y,e.shift_key()));
                                if let Some(cursor)=current_cursor(){base_drag.set_value(Some(BaseDrag{identity:identity.get_untracked(),revision:source.get_untracked(),pointer_id:e.pointer_id(),anchor:cursor.pointer_anchor.unwrap_or(cursor.anchor)}));drag_pointer.set_value(Some((e.client_x()as f64,e.client_y()as f64)));if c.set_pointer_capture(e.pointer_id()).is_ok(){drag_active.set(true);}else{stop_drag();}}
                            }
                        }
                    }
                    on:pointermove=move|e:web_sys::PointerEvent|{
                        let Some(drag)=base_drag.get_value()else{return;};
                        if drag.pointer_id!=e.pointer_id(){return;}if e.buttons()&1==0{stop_drag();return;}
                        if drag.identity!=identity.get_untracked()||drag.revision!=source.get_untracked()||locked.get_untracked(){stop_drag();return;}
                        if let Some(c)=canvas.get(){
                            e.prevent_default();let rect=c.get_bounding_client_rect();let x=e.client_x()as f64-rect.left();let y=e.client_y()as f64-rect.top();
                            drag_pointer.set_value(Some((e.client_x()as f64,e.client_y()as f64)));
                            if rect.width()<1.||rect.height()<1.{stop_drag();return;}
                            if let Some(base)=base_at(x,y.clamp(0.,(rect.height()-1.).max(0.)),scroll.get_untracked(),row_height.get_untracked(),columns.get_untracked(),molecule.with_untracked(|m|m.sequence.len()),true){let(a,p)=if base>=drag.anchor{(drag.anchor,base+1)}else{(drag.anchor+1,base)};set_cursor.run((a,p,Some(drag.anchor)));}
                        }
                    }
                    on:pointerup=move|e:web_sys::PointerEvent|{
                        let Some(drag)=base_drag.get_value()else{return;};if drag.pointer_id!=e.pointer_id(){return;}
                        if let Some(c)=canvas.get(){
                            if drag.identity==identity.get_untracked()&&drag.revision==source.get_untracked()&&!locked.get_untracked(){let rect=c.get_bounding_client_rect();if let Some(base)=base_at(e.client_x()as f64-rect.left(),(e.client_y()as f64-rect.top()).clamp(0.,(rect.height()-1.).max(0.)),scroll.get_untracked(),row_height.get_untracked(),columns.get_untracked(),molecule.with_untracked(|m|m.sequence.len()),true){let(a,p)=if base>=drag.anchor{(drag.anchor,base+1)}else{(drag.anchor+1,base)};set_cursor.run((a,p,Some(drag.anchor)));}}
                            stop_drag();
                        }else{stop_drag();}
                    }
                    on:pointercancel=move|_|stop_drag()
                    on:lostpointercapture=move|_|stop_drag()
                    on:click=move|e:web_sys::MouseEvent|{
                        // Preserve programmatic/accessibility clicks; a real pointer click
                        // was already handled on down/up and must not collapse its drag.
                        if e.detail()!=0 || locked.get_untracked(){return;}
                        if let Some(c)=canvas.get(){let rect=c.get_bounding_client_rect();let x=e.client_x()as f64-rect.left();let y=e.client_y()as f64-rect.top();let hit=hit_boxes.with_value(|boxes|boxes.iter().rev().find(|b|x>=b.x&&x<b.x+b.width&&y>=b.y&&y<b.y+b.height).map(|b|b.item.clone()));if let Some(hit)=hit{edit_cursor.set_value(None);choose.run(hit);}else{press_base.run((x,y,e.shift_key()));}}
                    }></canvas>
            </div>
            <div class="sr-only" aria-live="polite">{move||{if let Some(feature)=selected_feature.get().and_then(|i|molecule.with(|m|m.features.get(i).cloned())){let size=feature.location.spans().iter().try_fold(0usize,|n,(a,b,_)|n.checked_add(b-a)).unwrap_or(usize::MAX);let sequence=if size<=10_000{molecule.with(|m|feature.location.extract(&m.sequence).unwrap_or_else(|e|e.to_string()))}else{"Sequence display limit exceeded".into()};return format!("Selected annotation {}: {}",feature.label(),sequence);}selection.get().map(|(a,b)|molecule.with(|m|{if a==b{return format!("Insertion caret at boundary {} of {}. Type to insert reference symbols.",a,m.sequence.len());}if b.saturating_sub(a)>10_000 {format!("Selected {} {} to {}. Sequence display limit exceeded; use a sequence export for the complete selection.",m.sequence.alphabet().selected_symbols(),a+1,b)}else{format!("Selected {} {} to {}: {}",m.sequence.alphabet().selected_symbols(),a+1,b,m.sequence.as_str().get(a..b).unwrap_or("selection changed"))}})).unwrap_or_else(||format!("{} reference parts highlighted",highlighted_spans.get().len()))}}</div>
        </div>
    }
}
