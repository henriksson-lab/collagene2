//! Visible-row drawing only; exact source guards live in the shared model/view.
use crate::biology::antibody::antibody_overlay;
use crate::biology::antibody::antibody_overlay::Overlay;
use web_sys::CanvasRenderingContext2d as Ctx;
pub fn row(ctx: &Ctx, overlay: &Overlay, start: usize, end: usize, x: f64, y: f64, cell: f64) {
    for index in
        start.max(overlay.track.first_supported_source_index)..end.min(overlay.track.source_end)
    {
        let Some(residue) = overlay.at(index) else {
            continue;
        };
        let left = x + (index - start) as f64 * cell;
        if let Some(color) = residue.region.as_ref().and_then(antibody_overlay::color) {
            ctx.set_fill_style_str(color);
            ctx.fill_rect(left, y - 13., cell, 17.);
        }
        if residue.assignment.is_none() || residue.region.is_none() {
            ctx.set_fill_style_str("#808080");
            ctx.fill_rect(left, y + 3., cell - 1., 2.);
        }
    }
    ctx.set_font("10px monospace");
    ctx.set_fill_style_str("#52635d");
    ctx.set_text_align("left");
    for label in overlay.labels(start, end, cell) {
        let _ = ctx.fill_text(&label.text, x + label.left, y + 19.);
    }
}
