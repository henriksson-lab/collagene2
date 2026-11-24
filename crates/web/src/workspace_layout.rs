//! Desktop display preferences and pane controls. No molecule/history/worker
//! state is read or written. Register only with the matching preference helpers.
use crate::workspace_preferences::{self, Preferences, ViewMode};
use leptos::prelude::*;
use std::{cell::RefCell, rc::Rc};
use wasm_bindgen::{JsCast, closure::Closure};
use wasm_bindgen_futures::spawn_local;
use web_sys::{Element, HtmlElement, IdbDatabase};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Pane {
    Documents,
    Sequence,
}

struct Drag {
    pane: Pane,
    pointer: i32,
    element: Element,
    start_x: f64,
    start_actual: f64,
    available: f64,
    pending: f64,
}
struct State {
    alive: bool,
    prefs: Preferences,
    db: Option<IdbDatabase>,
    pending_save: Option<Vec<u8>>,
    writing: bool,
    saved: Option<Vec<u8>>,
    drag: Option<Drag>,
    frame_pending: bool,
}

#[derive(Clone, Copy)]
pub struct Layout {
    /// Root must also disable its existing mode/library buttons until ready.
    pub ready: RwSignal<bool>,
    /// Empty on success. Read/write failures leave view controls usable locally.
    pub status: RwSignal<String>,
    view: RwSignal<&'static str>,
    library_visible: RwSignal<bool>,
    locked: Memo<bool>,
    desktop: RwSignal<bool>,
    document_now: RwSignal<u16>,
    document_max: RwSignal<u16>,
    sequence_now: RwSignal<u16>,
    sequence_min: RwSignal<u16>,
    sequence_max: RwSignal<u16>,
    sequence_available: RwSignal<f64>,
    sequence_usable: RwSignal<bool>,
    state: StoredValue<Rc<RefCell<State>>, LocalStorage>,
}

struct Listeners {
    observer: web_sys::ResizeObserver,
    _observer_callback: Closure<dyn FnMut(js_sys::Array, web_sys::ResizeObserver)>,
    window: web_sys::Window,
    resize: Closure<dyn FnMut(web_sys::Event)>,
}
impl Drop for Listeners {
    fn drop(&mut self) {
        self.observer.disconnect();
        let _ = self
            .window
            .remove_event_listener_with_callback("resize", self.resize.as_ref().unchecked_ref());
    }
}

fn release(drag: Drag) {
    let _ = drag.element.release_pointer_capture(drag.pointer);
}

impl Layout {
    fn enabled(self, pane: Pane) -> bool {
        self.ready.try_get().unwrap_or(false)
            && !self.locked.try_get().unwrap_or(true)
            && self.desktop.try_get().unwrap_or(false)
            && match pane {
                Pane::Documents => self.library_visible.try_get().unwrap_or(false),
                Pane::Sequence => {
                    self.view.try_get() == Some("split")
                        && self.sequence_usable.try_get().unwrap_or(false)
                }
            }
    }
    fn bounds(self, pane: Pane) -> (u16, u16, u16) {
        match pane {
            Pane::Documents => (
                180,
                self.document_max.get_untracked(),
                self.document_now.get_untracked(),
            ),
            Pane::Sequence => (
                self.sequence_min.get_untracked(),
                self.sequence_max.get_untracked(),
                self.sequence_now.get_untracked(),
            ),
        }
    }
    fn tracked_bounds(self, pane: Pane) -> (u16, u16, u16) {
        match pane {
            Pane::Documents => (180, self.document_max.get(), self.document_now.get()),
            Pane::Sequence => (
                self.sequence_min.get(),
                self.sequence_max.get(),
                self.sequence_now.get(),
            ),
        }
    }
    fn cancel(self) {
        let state = self.state.get_value();
        let drag = state.borrow_mut().drag.take();
        if let Some(drag) = drag {
            release(drag);
        }
        refresh(self);
    }
    fn commit(self, pane: Pane, value: u16) {
        let state = self.state.get_value();
        let changed = {
            let mut state = state.borrow_mut();
            let target = match pane {
                Pane::Documents => &mut state.prefs.library_width,
                Pane::Sequence => &mut state.prefs.split_per_mille,
            };
            let value = match pane {
                Pane::Documents => value.clamp(180, 360),
                Pane::Sequence => value.clamp(100, 900),
            };
            if *target == value {
                false
            } else {
                *target = value;
                true
            }
        };
        refresh(self);
        if changed {
            queue_save(self, state);
        }
    }
}

fn queue_frame(layout: Layout) {
    let state = layout.state.get_value();
    {
        let mut value = state.borrow_mut();
        if !value.alive || value.frame_pending {
            return;
        }
        value.frame_pending = true;
    }
    let action = Closure::once_into_js(move || {
        {
            let mut value = state.borrow_mut();
            if !value.alive {
                return;
            }
            value.frame_pending = false;
        }
        refresh(layout);
    });
    if let Some(window) = web_sys::window() {
        if window
            .request_animation_frame(action.unchecked_ref())
            .is_err()
        {
            layout.state.get_value().borrow_mut().frame_pending = false;
        }
    }
}

fn rect_width(document: &web_sys::Document, selector: &str) -> f64 {
    document
        .query_selector(selector)
        .ok()
        .flatten()
        .map(|e| e.get_bounding_client_rect().width())
        .unwrap_or(0.)
}

fn refresh(layout: Layout) {
    let Some(window) = web_sys::window() else {
        return;
    };
    let Some(document) = window.document() else {
        return;
    };
    let state = layout.state.get_value();
    let value = state.borrow();
    if !value.alive {
        return;
    }
    let width = window
        .inner_width()
        .ok()
        .and_then(|w| w.as_f64())
        .unwrap_or(0.);
    let document_max = ((width - 828.).floor() as i32).clamp(180, 360) as u16;
    let mut document_now = value.prefs.library_width.clamp(180, document_max);
    // Sum the actual two grid panes, avoiding assumptions about inspector width,
    // border rounding or the separating control's width.
    let available = rect_width(&document, ".editor-grid .sequence-panel")
        + rect_width(&document, ".editor-grid .map-panel");
    let minimum = if available > 0. {
        (210. / available * 1000.).ceil() as u16
    } else {
        900
    };
    let maximum = if available > 0. {
        ((available - 220.) / available * 1000.).floor().max(0.) as u16
    } else {
        100
    };
    let minimum = minimum.clamp(100, 900);
    let maximum = maximum.clamp(100, 900);
    let usable = available >= 430. && minimum <= maximum;
    let (minimum, maximum) = if usable {
        (minimum, maximum)
    } else {
        (100, 900)
    };
    let mut sequence_now = if usable {
        value.prefs.split_per_mille.clamp(minimum, maximum)
    } else {
        value.prefs.split_per_mille
    };
    if let Some(drag) = &value.drag {
        match drag.pane {
            Pane::Documents => {
                document_now = (drag.pending.round() as u16).clamp(180, document_max)
            }
            Pane::Sequence if usable => {
                sequence_now = (drag.pending.round() as u16).clamp(minimum, maximum)
            }
            Pane::Sequence => {}
        }
    }
    drop(value);
    if let Ok(Some(shell)) = document.query_selector(".workspace-shell") {
        if let Some(shell) = shell.dyn_ref::<HtmlElement>() {
            let style = shell.style();
            let _ = style.set_property("--document-width", &format!("{document_now}px"));
            let _ = style.set_property("--sequence-fr", &format!("{sequence_now}fr"));
            let _ = style.set_property("--map-fr", &format!("{}fr", 1000 - sequence_now));
        }
    }
    batch(|| {
        if layout.desktop.get_untracked() != (width >= 1100.) {
            layout.desktop.set(width >= 1100.);
        }
        if layout.document_max.get_untracked() != document_max {
            layout.document_max.set(document_max);
        }
        if layout.document_now.get_untracked() != document_now {
            layout.document_now.set(document_now);
        }
        if layout.sequence_available.get_untracked() != available {
            layout.sequence_available.set(available);
        }
        if layout.sequence_usable.get_untracked() != usable {
            layout.sequence_usable.set(usable);
        }
        if layout.sequence_min.get_untracked() != minimum {
            layout.sequence_min.set(minimum);
        }
        if layout.sequence_max.get_untracked() != maximum {
            layout.sequence_max.set(maximum);
        }
        if layout.sequence_now.get_untracked() != sequence_now {
            layout.sequence_now.set(sequence_now);
        }
    });
}

fn queue_save(layout: Layout, state: Rc<RefCell<State>>) {
    let encoded = workspace_preferences::encode(&state.borrow().prefs);
    let bytes = match encoded {
        Ok(bytes) => bytes,
        Err(e) => {
            layout.status.set(e);
            return;
        }
    };
    {
        let mut value = state.borrow_mut();
        if !value.alive {
            return;
        }
        if value.db.is_none() {
            layout
                .status
                .set("View changes are temporary: browser storage is unavailable".into());
            return;
        }
        if !value.writing && value.saved.as_ref() == Some(&bytes) {
            value.pending_save = None;
            return;
        }
        // Keep newest committed request even if it equals the last saved bytes:
        // an already in-flight DIFFERENT request may still overwrite that value.
        value.pending_save = Some(bytes);
        if value.writing {
            return;
        }
        value.writing = true;
    }
    spawn_local(async move {
        loop {
            let next = {
                let mut value = state.borrow_mut();
                if !value.alive {
                    value.writing = false;
                    return;
                }
                match value.pending_save.take() {
                    Some(bytes) => value.db.clone().map(|db| (db, bytes)),
                    None => {
                        value.writing = false;
                        None
                    }
                }
            };
            let Some((db, bytes)) = next else {
                return;
            };
            match crate::storage::save_view_preferences(&db, &bytes).await {
                Ok(()) => {
                    let mut value = state.borrow_mut();
                    value.saved = Some(bytes);
                    if value.alive {
                        layout.status.try_set(String::new());
                    }
                }
                Err(error) => {
                    let mut value = state.borrow_mut();
                    value.writing = false;
                    if value.alive {
                        layout
                            .status
                            .try_set(format!("View changes are temporary: {error}"));
                    }
                    return;
                }
            }
        }
    });
}

/// Call once at App scope; no child context is required. Preferences live only in
/// the separate browser view-preferences record, never a workspace/archive event.
/// Root must disable mode/library controls until `ready` to prevent input racing
/// the initial asynchronous read. Missing preferences are NOT written on mount.
pub fn use_layout(
    view: RwSignal<&'static str>,
    library_visible: RwSignal<bool>,
    locked: Memo<bool>,
) -> Layout {
    let state = Rc::new(RefCell::new(State {
        alive: true,
        prefs: Preferences::default(),
        db: None,
        pending_save: None,
        writing: false,
        saved: None,
        drag: None,
        frame_pending: false,
    }));
    let layout = Layout {
        ready: RwSignal::new(false),
        status: RwSignal::new("Loading view preferences…".into()),
        view,
        library_visible,
        locked,
        desktop: RwSignal::new(false),
        document_now: RwSignal::new(244),
        document_max: RwSignal::new(360),
        sequence_now: RwSignal::new(524),
        sequence_min: RwSignal::new(100),
        sequence_max: RwSignal::new(900),
        sequence_available: RwSignal::new(0.),
        sequence_usable: RwSignal::new(false),
        state: StoredValue::new_local(state.clone()),
    };
    let listeners = StoredValue::new_local(None::<Listeners>);
    let read_state = state.clone();
    spawn_local(async move {
        let loaded = async {
            let db = crate::storage::open().await?;
            let bytes = match crate::storage::load_view_preferences(&db).await {
                Ok(bytes) => bytes,
                Err(error) => {
                    db.close();
                    return Err(error);
                }
            };
            let preferences = match bytes
                .as_ref()
                .map(|b| workspace_preferences::decode(b))
                .transpose()
            {
                Ok(preferences) => preferences,
                Err(error) => {
                    db.close();
                    return Err(error);
                }
            };
            Ok::<_, String>((db, bytes, preferences))
        }
        .await;
        if !read_state.borrow().alive {
            if let Ok((db, _, _)) = loaded {
                db.close();
            }
            return;
        }
        let (preferences, message) = match loaded {
            Ok((db, bytes, preferences)) => {
                let preferences = preferences.unwrap_or_default();
                let mut value = read_state.borrow_mut();
                value.db = Some(db);
                value.saved = bytes;
                value.prefs = preferences.clone();
                (preferences, String::new())
            }
            Err(error) => (
                Preferences::default(),
                format!("View preferences unavailable; using temporary defaults: {error}"),
            ),
        };
        batch(|| {
            view.set(preferences.view.as_str());
            library_visible.set(preferences.library_visible);
            layout.status.set(message);
            layout.ready.set(true);
        });
        queue_frame(layout);
    });
    // View changes are committed controls, not pointer previews. Initial loaded
    // values already equal state.prefs, so this effect never writes on mount.
    Effect::new(move |_| {
        let mode = view.get();
        let visible = library_visible.get();
        let ready = layout.ready.get();
        let locked = locked.get();
        let desktop = layout.desktop.get();
        if !ready {
            return;
        }
        let state = layout.state.get_value();
        let cancel = state.borrow().drag.as_ref().is_some_and(|d| {
            locked
                || !desktop
                || (d.pane == Pane::Documents && !visible)
                || (d.pane == Pane::Sequence && mode != "split")
        });
        if cancel {
            layout.cancel();
        }
        if locked {
            return;
        }
        let Some(mode) = ViewMode::from_str(mode) else {
            return;
        };
        let changed = {
            let mut value = state.borrow_mut();
            let changed = value.prefs.view != mode || value.prefs.library_visible != visible;
            value.prefs.view = mode;
            value.prefs.library_visible = visible;
            changed
        };
        if changed {
            queue_save(layout, state);
        }
        queue_frame(layout);
    });
    let observer_state = state.clone();
    Effect::new(move |_| {
        let _ = layout.ready.get();
        if listeners.with_value(|v| v.is_some()) {
            return;
        }
        let Some(window) = web_sys::window() else {
            return;
        };
        let Some(document) = window.document() else {
            return;
        };
        let Ok(Some(shell)) = document.query_selector(".workspace-shell") else {
            return;
        };
        let retained = observer_state.clone();
        let callback =
            Closure::<dyn FnMut(js_sys::Array, web_sys::ResizeObserver)>::new(move |_, _| {
                let alive = retained.borrow().alive;
                if alive {
                    queue_frame(layout);
                }
            });
        let Ok(observer) = web_sys::ResizeObserver::new(callback.as_ref().unchecked_ref()) else {
            return;
        };
        observer.observe(&shell);
        if let Ok(Some(grid)) = document.query_selector(".editor-grid") {
            observer.observe(&grid);
        }
        let retained = observer_state.clone();
        let resize = Closure::<dyn FnMut(web_sys::Event)>::new(move |_| {
            let alive = retained.borrow().alive;
            if alive {
                queue_frame(layout);
            }
        });
        let _ = window.add_event_listener_with_callback("resize", resize.as_ref().unchecked_ref());
        listeners.set_value(Some(Listeners {
            observer,
            _observer_callback: callback,
            window,
            resize,
        }));
        queue_frame(layout);
    });
    on_cleanup(move || {
        listeners.update_value(|value| *value = None);
        let state = layout.state.get_value();
        let (drag, db) = {
            let mut value = state.borrow_mut();
            value.alive = false;
            value.pending_save = None;
            (value.drag.take(), value.db.take())
        };
        if let Some(drag) = drag {
            release(drag);
        }
        if let Some(db) = db {
            db.close();
        }
    });
    layout
}

#[component]
pub fn PaneResize(layout: Layout, pane: Pane) -> impl IntoView {
    let cleanup_state = layout.state;
    on_cleanup(move || {
        let state = cleanup_state.get_value();
        let drag = {
            let mut value = state.borrow_mut();
            if value.drag.as_ref().is_some_and(|d| d.pane == pane) {
                value.drag.take()
            } else {
                None
            }
        };
        if let Some(drag) = drag {
            release(drag);
            let alive = state.borrow().alive;
            if alive {
                queue_frame(layout);
            }
        }
    });
    let pointer_down = move |event: web_sys::PointerEvent| {
        if !layout.enabled(pane) || !event.is_primary() || event.button() != 0 {
            return;
        }
        let state = layout.state.get_value();
        if state.borrow().drag.is_some() {
            return;
        }
        let Some(element) = event
            .current_target()
            .and_then(|e| e.dyn_into::<Element>().ok())
        else {
            return;
        };
        if element.set_pointer_capture(event.pointer_id()).is_err() {
            return;
        }
        let (_, _, actual) = layout.bounds(pane);
        state.borrow_mut().drag = Some(Drag {
            pane,
            pointer: event.pointer_id(),
            element,
            start_x: event.client_x() as f64,
            start_actual: actual as f64,
            available: layout.sequence_available.get_untracked(),
            pending: actual as f64,
        });
        event.prevent_default();
    };
    let pointer_move = move |event: web_sys::PointerEvent| {
        let state = layout.state.get_value();
        let mut value = state.borrow_mut();
        let Some(drag) = value
            .drag
            .as_mut()
            .filter(|d| d.pane == pane && d.pointer == event.pointer_id())
        else {
            return;
        };
        let delta = event.client_x() as f64 - drag.start_x;
        drag.pending = drag.start_actual
            + match pane {
                Pane::Documents => delta,
                Pane::Sequence => {
                    if drag.available > 0. {
                        delta / drag.available * 1000.
                    } else {
                        0.
                    }
                }
            };
        drop(value);
        event.prevent_default();
        queue_frame(layout);
    };
    let pointer_up = move |event: web_sys::PointerEvent| {
        let state = layout.state.get_value();
        if !state
            .borrow()
            .drag
            .as_ref()
            .is_some_and(|d| d.pane == pane && d.pointer == event.pointer_id())
        {
            return;
        }
        let drag = state
            .borrow_mut()
            .drag
            .take()
            .expect("matching retained drag");
        let (min, max, _) = layout.bounds(pane);
        let value = (drag.pending.round() as u16).clamp(min, max);
        let changed = value != drag.start_actual.round() as u16;
        release(drag);
        if changed && layout.enabled(pane) {
            layout.commit(pane, value);
        } else {
            refresh(layout);
        }
        event.prevent_default();
    };
    let cancel = move |event: web_sys::PointerEvent| {
        let state = layout.state.get_value();
        if state
            .borrow()
            .drag
            .as_ref()
            .is_some_and(|d| d.pane == pane && d.pointer == event.pointer_id())
        {
            layout.cancel();
        }
    };
    view! {
        <button type="button" role="separator" aria-orientation="vertical"
            class=match pane {Pane::Documents=>"pane-resize pane-resize-documents",Pane::Sequence=>"pane-resize pane-resize-sequence"}
            aria-label=match pane {Pane::Documents=>"Resize document panel",Pane::Sequence=>"Resize sequence and map"}
            title="Drag or use Left/Right arrows; Shift moves farther; Home/End reach bounds; double-click resets"
            disabled=move||!layout.enabled(pane)
            style:display=move||if layout.desktop.get() && match pane {Pane::Documents=>layout.library_visible.get(),Pane::Sequence=>layout.view.get()=="split" && layout.sequence_usable.get()} {"block"}else{"none"}
            aria-valuemin=move||{let(min,_,_)=layout.tracked_bounds(pane);if pane==Pane::Documents{min.to_string()}else{format!("{:.1}",min as f64/10.)}}
            aria-valuemax=move||{let(_,max,_)=layout.tracked_bounds(pane);if pane==Pane::Documents{max.to_string()}else{format!("{:.1}",max as f64/10.)}}
            aria-valuenow=move||{let(_,_,now)=layout.tracked_bounds(pane);if pane==Pane::Documents{now.to_string()}else{format!("{:.1}",now as f64/10.)}}
            aria-valuetext=move||{let(_,_,now)=layout.tracked_bounds(pane);match pane {Pane::Documents=>format!("Document panel {now} pixels"),Pane::Sequence=>format!("Sequence {:.1} percent, map {:.1} percent",now as f64/10.,(1000-now) as f64/10.)}}
            on:pointerdown=pointer_down on:pointermove=pointer_move on:pointerup=pointer_up
            on:pointercancel=cancel on:lostpointercapture=cancel
            on:keydown=move|event| {
                if !layout.enabled(pane){return;}
                if event.key()=="Escape" {layout.cancel();event.prevent_default();return;}
                let (min,max,now)=layout.bounds(pane);
                let step=match(pane,event.shift_key()){(Pane::Documents,false)=>10,(Pane::Documents,true)=>50,(Pane::Sequence,false)=>20,(Pane::Sequence,true)=>100};
                let value=match event.key().as_str(){"ArrowLeft"=>Some(now.saturating_sub(step).max(min)),"ArrowRight"=>Some(now.saturating_add(step).min(max)),"Home"=>Some(min),"End"=>Some(max),_=>None};
                if let Some(value)=value {event.prevent_default();if value!=now{layout.cancel();layout.commit(pane,value);}}
            }
            on:dblclick=move|_|if layout.enabled(pane){layout.cancel();layout.commit(pane,match pane{Pane::Documents=>244,Pane::Sequence=>524});}
        ></button>
    }
}
