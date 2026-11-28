# Documentary image evidence

Pure Rust JPEG, PNG and TIFF/BigTIFF admission and preview generation for immutable
workspace attachments. MIT. Native and `wasm32-unknown-unknown`; no external process,
filesystem, network, system font, browser, or platform image-codec dependency.

`Evidence::prepare(original_bytes, Limits::default())` preserves the exact original
and BLAKE3 hash, records encoded dimensions and EXIF orientation, and creates a
small orientation-correct RGBA8 PNG. `original_bytes()` / `into_original()` re-export
the exact source. The owner cannot be constructed through unchecked deserialization.

The default input limit is 16 MiB, raster limit 4,194,304 pixels / 16,384 pixels per
axis, decoded-output admission 64 MiB, and conservative working admission 256 MiB.
The working estimate includes padded source/intermediate space and small preview;
thus not every raster below the pixel ceiling is admitted. The codec additionally
receives allocation and dimension limits. Codec allocation limits are best effort,
not a process-wide allocator guarantee. Callers may lower all settings, never raise
them beyond these ceilings. Preview limits are 1,000,000 pixels, 1,600 per axis and
5 MiB encoded. Preview sampling is nearest-neighbour, has no scientific quantitative
meaning, and does not apply ICC profiles, physical-unit or acquisition calibration.
Original high-bit samples remain available only in the original file.

Container admission happens before codec construction: bounded chunk/marker scans,
PNG dimensions, JPEG sampling/scan count, and classic/BigTIFF IFD offsets, typed value
sizes/counts, graph cycles, strip/tile ranges and tile intermediate sizes. Embedded
TIFF JPEG block headers are admitted independently; their dimensions cannot bypass
the outer image admission. Metadata IFD traversal is at most 32 directories, 4,096
entries, 65,536 values per entry and 2 MiB referenced values. Overlapping block
references and shared JPEG table replays have a bounded aggregate compressed extent. Unsupported pixel encodings
produce explicit errors; no external fallback decoder is invoked.

Only the first/default image is previewed. Additional TIFF pages/subimages and PNG
animation images generate a warning; other images are retained as original bytes,
not claimed to be decoded or verified. Unusual valid files can be rejected by these
resource bounds. Metadata stays in the original; it is intentionally not transferred
to the PNG preview. Metadata and uncalibrated-display warnings always accompany the
preview; precision and downsampling warnings appear when applicable. The preview
is documentary only: it supplies no band identification, quantification, sequence
verification, cloning success, sample identity, or physical-material authority.

## Host integration

This crate is registered and wired through core/archive/worker, MCP and the web
File → Images dialog. Native/WASM checks and the302-case connected batch pass.
The importing host should run `prepare` in its existing private compute worker and
publish a preview result only after context/fingerprint validation. Keep original
and preview as immutable hash-addressed workspace attachments with distinct kind /
media type. A bounded, current-version typed manifest should bind source and preview
hashes, format, dimensions, orientation, preview engine version and warnings, and an
explicit user-authored association to a document/run when requested. Validate both
attachment hashes and sizes on load; recompute preview or validate the manifest in
the trusted worker rather than deserializing an `Evidence` owner from remote JSON.

Attachment addition must leave molecule sequence, features, material origins and
scientific verification status unchanged. Existing workspace attachment/export and
MCP controllers can expose add/list/download operations with their current context
and authorization rules. Do not insert base64 blobs into molecule metadata or infer
verification from the image. PNG preview is suitable for a blob URL or escaped image
source; never render original metadata as HTML. No OpenGel comparison is included.

Pin: `image = 0.25.10`, default features disabled, only `jpeg`, `png`, `tiff` (MIT or
Apache-2.0; Rust 1.88). Source hash and errors use existing workspace dependencies.
The shared native/WASM build and `evidence` integration harness pass. The
repository integration record retains the exact executed cases and source scope. Hand-built PNG (stored DEFLATE), grayscale JPEG (DC/EOB Huffman), and
classic/BigTIFF little/big endian fixtures are independent of image's encoders.
