# Desktop MVP: sequence viewing and virtual digest

The current delivery priority is an ordinary desktop workflow: open GenBank, inspect the sequence and map, select enzymes, inspect a virtual digest and explicitly open/export a fragment. Features and Digest are primary tools. The metadata registry in `desktop_shell.rs` keeps advanced modules behind More tools without removing their algorithms or changing their stable identities. Mobile scope and acceptance remain deferred until desktop completion and user review.

## Workspace layout in development

The early preview was not accepted as a finished GUI. The desktop frame now
separates a narrow functional navigation rail, a searchable document explorer,
document switching tabs and the sequence editor. There is one simple workspace;
no project hierarchy or excluded inventory/notebook modules are implied.

The explorer searches names, stable document IDs, recorded sources and molecule
types, and sorts by workspace order, name, length or source. Its bounded pages
contain at most 100 rows. Tab order follows the workspace, independently of the
explorer's current search and sort. Arrow keys, Home and End switch and focus
tabs. Collapsing the explorer changes only layout. Rows and tabs address stable
document IDs and react to current metadata; they do not create, close or delete
documents, alter science or replace source history. Imported source labels are
shown only when recorded in the document.

`scripts/workspace_navigation_browser_qa.py` exercises the actual imported-file,
pointer and keyboard workflow, including native history/source preservation
across navigation and reload, and an intentional same-ID rename. Passing it is
a functional gate for this frame, not acceptance of the overall GUI design.

## Annotation display layout

The sequence view prepares a whole-document annotation plan with the independent
[`collagene-view-layout`](../crates/view-layout/README.md) library. Nonoverlapping
annotations share rows; every part of a joined or circular-split annotation stays
in one row. Primer annotations have a separate row budget. Counts reflect rows
actually used, so empty annotation and primer rows collapse. Source annotations
remain available in Features; the `source` annotation kind is omitted from the
sequence tracks, matching the existing scientific viewport index.

**Sequence tracks → Annotation rows / Primer rows** sets display maxima. Additional
annotations are identified explicitly and remain inspectable in Features and in
the bounded viewport report. Those reports retain their own completeness and
warning semantics. Unknown local positions remain unplaced evidence. A display
row never implies a biological relationship or higher annotation confidence.

The plan is cached by the document and display settings, independently of incoming
viewport results. Reflow preserves the visible source-base anchor for an unchanged
document/reference, with normal clamping at the document end. Original DNA,
annotations, source records and revision history retain their semantics. The
layout library has bounded sorting, ordered interval indexes and preflight work
and storage limits; the adapter checks borrowed feature/location bounds before
constructing interval vectors. Both support the existing 1 Mb/10,000-feature case.

`scripts/annotation_layout_browser_qa.py` uses actual GenBank imports, delayed real
worker messages, exact source-feature identities, visible annotation/base clicks,
middle/end reflows and a 1 Mb/10,000-annotation record. It compares native source
and history before/after navigation and reload.

## Boundaries retained for later features

- `Molecule`, `Sequence` and `Location` remain the common scientific model. GenBank adapters retain original records, topology, annotations and unsupported evidence; the renderer never owns the biological coordinate model.
- `restriction_inspection` supplies bounded multi-enzyme recognition and lazy physical-product reconstruction. The previous strict analysis and chemistry APIs remain separate and compatible. Enzyme definitions are data, including exact recognition and strand-cut geometry.
- Worker schema 20 supplied typed inspection, fragment materialization and source-event installation. Schema 21 appends catalog-pinned variants and retains those earlier layouts. Normal scientific envelopes remain bounded at 32 MiB; only the explicit installation variants join the 128 MiB workspace-message whitelist.
- `restriction_history` uses the existing `Workspace`, revision events and document histories. It reconstructs a selected fragment afresh and retains its exact source event as its sole computational parent. History labels and producer records cannot restore physical or verified-plan authority.
- The digest panel controls requests and evidence pages; `digest_product` handles guarded atomic native installation, while `storage` supplies the reusable cancellable current/previous transaction. Later workflows can use the same state and persistence seams.

## Scientific interpretation

Recognition enumeration and determined physical products are separate. A possible IUPAC recognition site, an out-of-substrate Type IIS cut, unsupported palindromic strand geometry, a work/output limit or unavailable physical substrate can leave useful site evidence while withholding products. A complete scan is not a statement that a laboratory digest will work. Blocking methylation, supplier buffers, reaction conditions and experimental success are not inferred by the assumed-complete model.

Reference lengths are counts along the declared top reference strand. Product summaries separately retain top/bottom lengths, original cut boundaries and end geometry. A recognition site is selected using its full motif and original circular spans, not an assumed six-base window. Multi-enzyme fragments retain deterministic source-cut ordering.

Product construction runs fresh chemistry and the established annotation projection. Partial/unresolved annotations retain review findings; derived products do not reuse source GenBank headers or source coordinate records as their own. Native products retain a bounded versioned recipe with exact source, enzyme geometry, options and fragment identity. Annotated exports retain unverified recipe/end claims; importing those claims does not reconstruct chemical state. FASTA is a sequence-only exchange and its losses remain explicit.

## Current-state and persistence rules

A request binds the exact document event, sequence revision, enzyme definitions and options. Switching documents, advancing an equal-content event, changing options or restoring a workspace invalidates prior results. Source/epoch/envelope guards reject obsolete responses. Fragment opening operates on an immutable snapshot, reconstructs in the worker and prepares a native archive before committing it. Current and previous recovery slots commit together; cancellation can abort an uncommitted replacement. Validation/capacity/quota errors preserve the prior workspace.

The viewer and format adapters retain their existing limits and extension points. This MVP is a delivery gate within the full project plan; it does not replace the remaining cloning, catalog, collaboration, scale or Tauri work.

## Validation checkpoint

The release WASM/web build, native workspace checks (1,040 tests plus 14 doctests), real stdio MCP smoke and all 30 desktop browser workflow suites pass. The dedicated desktop MVP browser suite exercises actual GenBank import, map rotation and sequence navigation, independently specified multi-enzyme products, non-six-base and circular-origin sites, fresh fragment construction/export, exact source history, ambiguous/limited/uncut cases, paged evidence, stale-result cancellation and atomic storage rollback. Desktop widths of 1600, 1280 and 1200 pixels are covered. These are local Chrome/software-renderer observations; they do not establish Firefox/Safari, hardware GPU or installed desktop parity.

The navigation suite additionally covers 1366×768 and 1280×800 laptop windows, a 101-record import across explorer pages, same-ID rename updates, keyboard tab focus and native source/history preservation. The canvas regression reproduces scrollbar-only client viewport changes without synthetic resize/scroll events, verifies backing/CSS dimensions and exact base selection, and checks natural short/long document switches. Fractional CSS changes with unchanged rounded client dimensions must also repaint. Backing dimensions round outward to whole pixels, while drawing and hit testing retain CSS-coordinate geometry.

Run the desktop suites against the isolated QA browser and release preview with `python3 scripts/desktop_qa.py`. See the README for prerequisites. Mobile and machine-specific performance acceptance remain separate from this gate.

## Catalog identity and release assets

The desktop digest now selects from the forty-identity built-in catalog. A selection pins the edition, full catalog content hash, ordered enzyme IDs and exact strand geometry. Catalog wrappers reuse the original inspection/physical-product algorithms; raw legacy inspection artifacts and recipe version 1 retain their layouts. Supplier observations have separate documentary contexts, qualifiers and citations and do not decide cleavage or compatibility. Fresh worker installation adds the catalog pin before deriving the exact-parent native history event. GenBank/EMBL carry bounded producer claims; restored claims remain unverified. [Desktop catalog workflow](CATALOG_DIGEST_UI.md).

Both isolated analysis jobs and retained viewport workers read one content-addressed bundle URL from the generated HTML. The build packages the loader, JavaScript glue and WASM as one immutable directory; no fixed worker URL fallback is accepted. The post-build step runs against staged assets before publication ([Trunk hook lifecycle](https://github.com/trunk-rs/trunk/blob/main/Trunk.toml)). Deployment should revalidate HTML and retain previously published hashed assets for already-open tabs. The release regression keeps the browser HTTP cache enabled and tests both worker factories with unusable legacy fixed URLs.
