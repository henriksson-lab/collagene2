//! Bounded, virtualized sequence window with explicit global contig coordinates.
use crate::viewer::map::renderer::ResizeListener;
use collagene_core::tracks::genome::{GenomicInterval, ReferenceLink};
use leptos::{html, prelude::*};
use std::sync::Arc;
use wasm_bindgen::JsCast;
#[derive(Clone)]
pub struct WindowDisplay {
    pub epoch: u64,
    pub link: ReferenceLink,
    pub cache_id: String,
    pub interval: GenomicInterval,
    pub contig_length: u64,
    pub sequence: String,
}
#[component]
pub fn GenomeWindowView(
    data: Arc<WindowDisplay>,
    on_position: Callback<(String, String, u64)>,
) -> impl IntoView {
    let canvas = NodeRef::<html::Canvas>::new();
    let scroll = NodeRef::<html::Div>::new();
    let top = RwSignal::new(0usize);
    let columns = RwSignal::new(80usize);
    let resized = RwSignal::new(0usize);
    let resize = StoredValue::new_local(None::<ResizeListener>);
    let draw = data.clone();
    Effect::new(move |_| {
        if let Some(element) = scroll.get() {
            if resize.with_value(|r| r.is_none()) {
                resize.set_value(
                    ResizeListener::new(element.unchecked_ref(), move || {
                        let _ = resized.try_update(|v| *v = v.wrapping_add(1));
                    })
                    .ok(),
                );
            }
        }
    });
    Effect::new(move |_| {
        let _ = resized.get();
        let Some(c) = canvas.get() else {
            return;
        };
        let Some(container) = scroll.get() else {
            return;
        };
        let offset = top.get();
        let width = container.client_width().max(200) as f64;
        let height = container.client_height().max(250) as f64;
        let per = ((width - 106.) / 9.).floor().clamp(8., 120.) as usize;
        columns.set(per);
        let ratio = web_sys::window()
            .map(|w| w.device_pixel_ratio())
            .unwrap_or(1.)
            .clamp(1., 2.);
        c.set_width((width * ratio) as u32);
        c.set_height((height * ratio) as u32);
        let _ = web_sys::HtmlElement::style(c.unchecked_ref::<web_sys::HtmlElement>())
            .set_property("height", &format!("{height}px"));
        let Ok(Some(ctx)) = c.get_context("2d") else {
            return;
        };
        let Ok(ctx) = ctx.dyn_into::<web_sys::CanvasRenderingContext2d>() else {
            return;
        };
        let _ = ctx.set_transform(ratio, 0., 0., ratio, 0., 0.);
        ctx.set_fill_style_str("#ffffff");
        ctx.fill_rect(0., 0., width, height);
        ctx.set_font("12px monospace");
        ctx.set_text_baseline("middle");
        let first = offset / 22;
        let y_offset = offset % 22;
        let visible = (height / 22.).ceil() as usize + 1;
        let bytes = draw.sequence.as_bytes();
        for row in first..first + visible {
            let start = row * per;
            if start >= bytes.len() {
                break;
            }
            let y = (row - first) as f64 * 22. + 11. - y_offset as f64;
            let global = (draw.interval.start + start as u64) % draw.contig_length;
            ctx.set_fill_style_str("#718076");
            let _ = ctx.fill_text(&(global + 1).to_string(), 4., y);
            for (i, b) in bytes.iter().skip(start).take(per).enumerate() {
                ctx.set_fill_style_str(match b {
                    b'A' => "#347d5e",
                    b'C' => "#3d6b9e",
                    b'G' => "#9a722b",
                    b'T' => "#aa5650",
                    _ => "#9a7c91",
                });
                let _ = ctx.fill_text(&(*b as char).to_string(), 106. + i as f64 * 9., y);
            }
        }
    });
    let len = data.sequence.len();
    let role = data.link.role.clone();
    let cache = data.cache_id.clone();
    let origin = data.interval.start;
    let contig_length = data.contig_length;
    let picked = data.clone();
    view! {<section class="genome-window"><p>{format!("Local window only · reference {} ({}) · {} · global source bases {}–{} (1-based) · {} of {} nt. Circular coordinates wrap at the source origin.",data.link.role,data.link.name,data.interval.contig_id,data.interval.start+1,(data.interval.end-1)%data.contig_length+1,data.sequence.len(),data.contig_length)}</p><p class="helper">"This bounded window is not a chromosome Molecule/document. Scrolling paints only visible rows and does not load the full contig or rerun guide/search science."</p><div class="genome-sequence-scroll" node_ref=scroll on:scroll=move|_|{if let Some(e)=scroll.get(){top.set(e.scroll_top().max(0)as usize);}}><div class="genome-sequence-spacer" style:height=move||format!("{}px",len.div_ceil(columns.get())*22)><canvas node_ref=canvas class="genome-sequence-canvas" tabindex="0" aria-label="Verified local genome window with global original coordinates" on:click=move|event|{let Some(c)=canvas.get()else{return;};let r=c.get_bounding_client_rect();let x=event.client_x()as f64-r.left();let y=event.client_y()as f64-r.top();if x<106.||y<0.{return;}let col=((x-106.)/9.).floor()as usize;let row=(top.get_untracked()+y as usize)/22;let index=row*columns.get_untracked()+col;if col<columns.get_untracked()&&index<picked.sequence.len(){on_position.run((role.clone(),cache.clone(),(origin+index as u64)%contig_length));}}/></div></div></section>}
}
