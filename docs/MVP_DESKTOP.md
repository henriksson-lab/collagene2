# Desktop MVP: sequence viewing and virtual digest

The current delivery priority is an ordinary desktop workflow: open GenBank, inspect the sequence and map, select enzymes, inspect a virtual digest and explicitly open/export a fragment. Features and Digest are primary tools. The metadata registry in `desktop_shell.rs` keeps advanced modules behind More tools without removing their algorithms or changing their stable identities. Mobile scope and acceptance remain deferred until desktop completion and user review.

## Boundaries retained for later features

- `Molecule`, `Sequence` and `Location` remain the common scientific model. GenBank adapters retain original records, topology, annotations and unsupported evidence; the renderer never owns the biological coordinate model.
- `restriction_inspection` supplies bounded multi-enzyme recognition and lazy physical-product reconstruction. The previous strict analysis and chemistry APIs remain separate and compatible. Enzyme definitions are data, including exact recognition and strand-cut geometry.
- Worker schema 20 adds typed inspection, fragment materialization and source-event installation. Normal scientific envelopes remain bounded at 32 MiB; the installation job alone joins the explicit 128 MiB workspace-message whitelist.
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

The release WASM/web build, native workspace checks (800 tests), real stdio MCP smoke and all 24 desktop browser workflow suites pass. The dedicated desktop MVP browser suite exercises actual GenBank import, map rotation and sequence navigation, independently specified multi-enzyme products, non-six-base and circular-origin sites, fresh fragment construction/export, exact source history, ambiguous/limited/uncut cases, paged evidence, stale-result cancellation and atomic storage rollback. Desktop widths of 1600, 1280 and 1200 pixels are covered. These are local Chrome/software-renderer observations; they do not establish Firefox/Safari, hardware GPU or installed desktop parity.

Run the desktop suites against the isolated QA browser and release preview with `python3 scripts/desktop_qa.py`. See the README for prerequisites. Mobile and machine-specific performance acceptance remain separate from this gate.
