# Sequence glyph rendering

The sequence editor now has a separate transparent WebGL2 glyph surface above
its existing2D selection, caret, cut-marker and annotation canvas. Both surfaces
share one sticky viewport box, so browser compositor scrolling moves them
together. Existing row, character, hit-test and keyboard geometry is retained.

A Canvas-measured raster atlas uses the existing14px monospace font and actual
alphabetic baseline bounds. Visible reference and complement glyphs are rendered
with instanced quads and their original colors. The atlas uploads on device-scale
changes rather than per character. Premultiplied alpha and matching blend state
keep the transparent overlay composited over the existing editor graphics.

The GPU path has explicit atlas, framebuffer, device-scale and instance bounds.
Initialization, upload or context-loss errors hide the overlay and draw the same
visible glyphs on2D in that paint pass. Context restoration permits reinitializing
the GPU resources. Listeners and GL objects are disposed on unmount; drawing
requests no additional scientific worker analysis or full-sequence copy.

Native/WASM compiler checks passed. Browser rendering, GPU hardware and frame-time
validation remain deferred; this implementation makes no measured performance
claim. Existing historical scrolling observations do not validate this new path.
