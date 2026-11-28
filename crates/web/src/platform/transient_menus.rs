//! Shared lifecycle for transient menus; ordinary disclosure panels are unaffected.
use leptos::prelude::*;
use wasm_bindgen::{JsCast, closure::Closure};
use web_sys::{Document, Element, Event, HtmlElement, KeyboardEvent};

const MENU: &str = "details[data-transient-menu]";
const OPEN_MENU: &str = "details[data-transient-menu][open]";

fn target(event: &Event) -> Option<Element> {
    let target = event.target()?;
    target.dyn_ref::<Element>().cloned().or_else(|| {
        target
            .dyn_ref::<web_sys::Node>()
            .and_then(|node| node.parent_element())
    })
}
fn closest(element: &Element, selector: &str) -> Option<Element> {
    element.closest(selector).ok().flatten()
}
fn open_menus(document: &Document) -> Vec<Element> {
    let Ok(nodes) = document.query_selector_all(OPEN_MENU) else {
        return Vec::new();
    };
    (0..nodes.length())
        .filter_map(|i| nodes.item(i)?.dyn_into::<Element>().ok())
        .collect()
}
fn close_branch(document: &Document, menu: &Element) {
    for open in open_menus(document) {
        if menu.contains(Some(&open)) {
            let _ = open.remove_attribute("open");
        }
    }
}
fn close_others(document: &Document, retained: Option<&Element>) {
    for open in open_menus(document) {
        // Retain the current menu and its ancestors, allowing nested submenus.
        if retained.is_none_or(|element| !open.contains(Some(element))) {
            let _ = open.remove_attribute("open");
        }
    }
}
fn summary(menu: &Element) -> Option<Element> {
    let mut child = menu.first_element_child();
    while let Some(element) = child {
        if element.tag_name().eq_ignore_ascii_case("summary") {
            return Some(element);
        }
        child = element.next_element_sibling();
    }
    None
}
fn direct_summary(element: &Element) -> Option<(Element, Element)> {
    let trigger = closest(element, "summary")?;
    let menu = trigger.parent_element()?;
    menu.matches(MENU).ok().filter(|matches| *matches)?;
    // Only the first direct summary is the native details trigger.
    (summary(&menu).as_ref() == Some(&trigger)).then_some((menu, trigger))
}

/// Install once under the App owner. All handlers are removed with that owner.
pub fn install() {
    let Some(document) = web_sys::window().and_then(|window| window.document()) else {
        return;
    };
    let mut listeners: Vec<(&'static str, Closure<dyn FnMut(Event)>)> = Vec::new();
    for name in ["pointerdown", "click", "keydown", "focusin", "toggle"] {
        let document = document.clone();
        let handler = Closure::<dyn FnMut(Event)>::new(move |event: Event| {
            let element = target(&event);
            match name {
                "pointerdown" | "focusin" => {
                    let retained = element
                        .as_ref()
                        .and_then(|target| closest(target, OPEN_MENU));
                    close_others(&document, retained.as_ref());
                }
                "toggle" => {
                    if let Some(menu) =
                        element.filter(|element| element.matches(OPEN_MENU).unwrap_or(false))
                    {
                        close_others(&document, Some(&menu));
                    }
                }
                "click" => {
                    let Some(element) = element else { return };
                    if let Some((menu, _)) = direct_summary(&element) {
                        // Cancel the native delayed toggle so coordination is synchronous.
                        event.prevent_default();
                        if menu.has_attribute("open") {
                            close_branch(&document, &menu);
                        } else {
                            close_others(&document, Some(&menu));
                            let _ = menu.set_attribute("open", "");
                        }
                        return;
                    }
                    if closest(&element, OPEN_MENU).is_none()
                        || closest(&element, "[data-menu-keep-open]").is_some()
                        || closest(&element, "input,label,select,textarea").is_some()
                    {
                        return;
                    }
                    let Some(action) = closest(&element, "button,a[href],[data-menu-dismiss]")
                    else {
                        return;
                    };
                    if action
                        .matches(":disabled,[aria-disabled='true']")
                        .unwrap_or(true)
                    {
                        return;
                    }
                    // Do not move focus or unmount anything: the action still receives
                    // this click after the capture handler closes its containing menu.
                    close_others(&document, None);
                }
                "keydown" => {
                    let Some(key) = event.dyn_ref::<KeyboardEvent>() else {
                        return;
                    };
                    if key.key() != "Escape" {
                        return;
                    }
                    let menu = element
                        .as_ref()
                        .and_then(|target| closest(target, OPEN_MENU))
                        .or_else(|| open_menus(&document).into_iter().last());
                    if let Some(menu) = menu {
                        event.prevent_default();
                        event.stop_propagation();
                        close_branch(&document, &menu);
                        if let Some(trigger) =
                            summary(&menu).and_then(|s| s.dyn_into::<HtmlElement>().ok())
                        {
                            let _ = trigger.focus();
                        }
                    }
                }
                _ => {}
            }
        });
        listeners.push((name, handler));
    }
    for (name, handler) in &listeners {
        let _ = document.add_event_listener_with_callback_and_bool(
            name,
            handler.as_ref().unchecked_ref(),
            true,
        );
    }
    let installed = StoredValue::new_local(Some((document, listeners)));
    on_cleanup(move || {
        installed.update_value(|installed| {
            if let Some((document, listeners)) = installed.take() {
                for (name, handler) in listeners {
                    let _ = document.remove_event_listener_with_callback_and_bool(
                        name,
                        handler.as_ref().unchecked_ref(),
                        true,
                    );
                }
            }
        });
    });
}
