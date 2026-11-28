use leptos::prelude::*;
use wasm_bindgen::JsCast;

pub(super) fn dismiss_viewer_menus(event: web_sys::MouseEvent) {
    let Some(target) = event
        .target()
        .and_then(|t| t.dyn_into::<web_sys::Element>().ok())
    else {
        return;
    };
    let menus =
        ".selection-action-menu,.toolbar-action-menu,.document-file-menu,.viewer-tools-menu";
    let inside = target.closest(menus).ok().flatten().is_some();
    if inside
        && (target.closest("summary").ok().flatten().is_some()
            || target.closest("button").ok().flatten().is_none())
    {
        return;
    }
    if let Some(document) = web_sys::window().and_then(|w| w.document()) {
        if let Ok(open)=document.query_selector_all(".selection-action-menu[open],.toolbar-action-menu[open],.document-file-menu[open],.viewer-tools-menu[open]") {
            for i in 0..open.length(){if let Some(node)=open.item(i).and_then(|n|n.dyn_into::<web_sys::Element>().ok()){let _=node.remove_attribute("open");}}
        }
    }
}
