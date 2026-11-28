# Web code organization audit

Date: 2026-11-28. Scope: `crates/web`, current working tree, including unfinished SBOL work. This is a static organization audit, not verification of biological correctness or of the uncommitted changes. No source edits, builds, tests, browser checks, or application requests were performed.

## Findings and priorities

Implementation update: draft production code now lives in the canonical
`assembly_draft_model/mod.rs` owner (1,594 formatted lines), with declaration
identity/comparison methods beside `Draft`. The workspace reexport facade and
model-as-view alias are removed. Nine extracted test modules preserve all 33
cases and their filters; all 33 passed in the native harness. The remaining
findings below describe the audit baseline and still guide the App/source-tree
refactor. WASM validation is tracked by the integration checkpoint.

1. **The directory limit is exceeded by more than an order of magnitude.** All **228 Rust files** are directly in `crates/web/src`; there are no source subdirectories. `main.rs` is effectively a growing flat module catalog. Four additional physical siblings are declared through `#[path]`: `alignment_index.rs`, `dotplot_projection.rs`, `history_graph.rs`, and `storage_browser.rs`. Include those in relocation, rather than counting only `mod` declarations.
2. **`app.rs` violates the absolute file limit:** **7,273 lines**, with `App` starting at line 1094 and its final shell view at line 6848. The file contains workspace restoration, document mutation, import classification, export infrastructure, history, PCR, assembly review, library operations, read review, rendering setup, and modal bodies. Static lexical counts are **124 `RwSignal::new`**, **18 `Effect::new`**, and **174 `Callback::new`** occurrences. These are indicators of orchestration concentration, not 124 independent state owners.
3. **`assembly_draft_model.rs` crosses the warning threshold:** **3,174 lines**. It combines `Draft`, `Slot`, capture/geometry types, bounded serde visitors and writers, validation, draft transformations, request construction, labels, and nine inline test modules. Some production functions appear between/after test modules, making ownership harder to follow. Extract coherent codec and test modules; do not split `impl Draft` into arbitrary line ranges.
4. **Several substantial files deserve monitoring but are below 3,000:** `cloning_methods.rs` 2475; `sequence_view.rs` 2231; `assembly_review.rs` 2083; `assembly_workspace.rs` 2031; `renderer.rs` 1752. The current refactor should resolve their responsibility boundaries without gratuitously splitting every file.

### Owners and dependencies that need correction

- **Assembly draft:** `Draft` is declared at `assembly_draft_model.rs:79`, its default at 219, and its main methods at 711; a second inherent implementation containing `identity` lives at `assembly_workspace.rs:33`. That component module reexports `Draft`, capture/geometry types, and a mutation helper at line 3. Consumers can therefore reach the same owner through two UI module paths. Move declaration identity and comparison to the draft owner, keep frozen review state private to the controller, and import the owner directly everywhere. Remove the facade when moving callers; this is not a reason to add another facade.
- **Preparation configuration:** `cloning_methods_model.rs` owns `MethodsState` and preparation types, while `cloning_methods.rs:28` glob reexports the entire model and mixes panel state, suggestion evidence, presentation, and async orchestration. Keep the model and its methods together; use explicit imports from that owner. Do not preserve the old flattened API with `pub use` aliases.
- **Workspace session:** `RestoredWorkspace` is at `app.rs:157`, restoration starts at 173, snapshot generation at 98, and application of restored state at 2205. Live workspace signals start at 1106; IndexedDB/save state at 1662. There is no cohesive live session owner today. Introduce a `WorkspaceSession` owning document histories, active document, workspace epoch/restoration state, attachments, retained candidates, and guarded mutation/snapshot operations. `RestoredWorkspace` should be the validated snapshot used by this owner, rather than state installed by many unrelated callbacks. Keep durable persistence mechanics in the existing `workspace_commit` module and storage adapter.
- **Document/view state:** tab reconciliation (`app.rs:1212` onward), sequence editing (1279), history editing (1326), restriction navigation (1502), explorer activation (1616), and tab closure (1629) share document identity and active-document transitions. Give a `DocumentController` the transitions and a `ViewerSelection` the selection/feature/search overlays. Presentation state must remain distinct from scientific source bindings and retained histories.
- **Assembly session:** `App` holds draft, previous draft, pending finalization, selected release, pending material, and multiple review callbacks (`app.rs:1151`, 4216–4410). `assembly_workspace.rs` also owns `IntentIdentity`, `Frozen`, and selection/review messages. Put draft replacement and guarded review/finalize coordination in an assembly controller; keep immutable scientific declarations in the draft model and verified staging in the review owner. Preserve exact source-event and workspace-epoch checks.
- **Export service:** `ExportArtifact` and its implementation already sit together (`app.rs:437`, 483), which is good. Their placement in `app` is wrong: feature modules call `crate::app::download_with_completion` (for example `alignment_view.rs:523`, `dotplot.rs:818`, and `cloning_methods.rs:1211`). Move the artifact, reservation, browser download, and completion APIs together into `io/download.rs`; update all consumers directly. Consider an `ExportService` owner for notice state instead of the thread-local signal at `app.rs:435`. Preserve bounded byte copies and native reservation-before-copy behavior.
- **Shared source identity:** `source_identity.rs` declares `SourceInfo`; `digest.rs:2` reexports it, and assembly/preparation modules then import it through `digest`. Source identity belongs to a shared source owner, not a digest panel. Move its sole definition to `sources/identity.rs` and update imports. Do not create parallel identities for PCR, assembly, and restriction.
- **Renderer state:** `MapRenderer`, private `State`, GPU cleanup, and interaction helpers are in one file. Keep that ownership. Helpers `advance_gesture(&mut State, ...)` at `renderer.rs:1251` and `draw(&mut State)` at 1494 have a natural `State` receiver; geometrical helpers such as arc tessellation remain free functions. Do not scatter one renderer across separate files per drawing primitive just to satisfy a count target.

### Concrete method opportunities

Review the resulting ownership before changing signatures; a common parameter alone is not sufficient evidence.

| Current item | Proposed owner/API | Reason |
| --- | --- | --- |
| `assembly_workspace::same_declaration` (83), `Draft::identity` (33) | Draft owner: `same_declaration(&self, other)`, declaration identity method | Both inspect the same draft declaration; active-slot UI state is deliberately excluded from declaration comparison. |
| `assembly_draft_model::set_digest_fragment` (2355), `orient_fragment` (3092), `set_material_origin` (3129) | `Draft` mutation/transformation methods | Validate and change that draft's slots; keep all guards and atomic replacement semantics. |
| `assembly_draft_model::inputs_request` (2862), `edit_request` (2992) | Draft request methods, or a single explicit draft-to-intent builder if independent context dominates | Draft owns the declaration being converted; do not disperse its rules among views. |
| `assembly_workspace::context` (133), `warnings` (141) | `IntendedReference` methods alongside `target` and `fragments` in `assembly_review.rs:19` | Removes repeated enum knowledge from the view/controller. |
| `app::workspace_snapshot`, `restore_workspace_state`, `commit_document`, `export_context_guard` | Session snapshot/restore/mutation methods and captured source guard | The same workspace/document state and revision checks recur in callbacks. Keep guard snapshots explicit across awaits. |
| `renderer::advance_gesture`, `draw` | Private `State` methods | Their primary state parameter already identifies the owner. |
| `app::prepare_download` | `ExportArtifact::prepare`, or `ExportService::prepare` when it owns notice state | Construction belongs with the bounded export artifact/reservation lifecycle. |

Pure geometry, symmetric comparisons without an owner, format encoding services, and UI event adapters may remain free functions. Separate model/view files are justified when they separate portable declarations from Leptos/browser lifetimes; they are not intrinsically a problem.

## Complete physical allocation

The table assigns **all 228 existing `.rs` files exactly once**. Names are current file stems. Counts exclude a new `mod.rs`; add one per non-root row. The largest proposed directory starts with 15 existing files and becomes **16**, leaving room below the absolute limit of 20. `assembly`, `cloning`, and `viewer` additionally need declaration-only parent `mod.rs` files. Prefer natural shorter names inside their new namespaces, updating imports directly.

| Destination under `src` | Existing files | File stems |
| --- | ---: | --- |
| `assembly/draft` | 10 | assembly_alternatives, assembly_alternatives_view, assembly_draft_actions, assembly_draft_model, assembly_multisite_view, assembly_pool_view, assembly_profile_view, assembly_role_ids, assembly_strategy_view, assembly_workspace |
| `assembly/edit` | 9 | assembly_edit_model, assembly_edit_view, assembly_fusion_design, assembly_fusion_view, assembly_junction_navigation_model, assembly_junctions_view, assembly_substitution_model, assembly_substitution_view, literal_intent_view |
| `assembly/primer_design` | 8 | assembly_edited_primer_commit, assembly_edited_primer_model, assembly_edited_primer_view, assembly_extension_primer_commit, assembly_extension_primer_model, assembly_extension_primer_view, assembly_inverse_pcr_model, assembly_inverse_pcr_view |
| `assembly/library` | 12 | assembly_batch_projection, assembly_batch_view, assembly_comparison_model, assembly_comparison_view, assembly_library_add_view, assembly_library_commit, assembly_library_history_view, assembly_library_topology, assembly_library_view, assembly_saved_designs_commit, assembly_saved_designs_model, assembly_saved_designs_view |
| `assembly/review` | 10 | assembly_archive_view, assembly_finalize, assembly_output_model, assembly_outputs, assembly_product_export_model, assembly_product_exports, assembly_product_model, assembly_product_view, assembly_release_view, assembly_review |
| `cloning/preparation` | 13 | anneal_preparation_view, cloning_methods, cloning_methods_model, cloning_preparations, end_geometry, end_repair, end_view, golden_gate_candidate_view, golden_gate_pool_form, ligation_outcomes_view, overlap_extension_preparation_view, pcr_digest_model, pcr_digest_view |
| `cloning/methods` | 15 | checked_cloning, crispr, domestication_form, domestication_review, domestication_view, multisite, multisite_form, mutagenesis, mutation_review, recombination, topo_actions_view, topo_graph_model, topo_graph_view, topo_preparation, topo_view |
| `restriction` | 9 | digest, digest_fragment_model, digest_product, gel, material_derivation_model, restriction_label_groups, restriction_map, restriction_view, reusable_material_view |
| `primers` | 15 | mismatch_pcr_outputs, pcr_outputs, pcr_overlay, pcr_overlay_draw, pcr_source_view, primer_binding_view, primer_domain_edit, primer_domain_edit_view, primer_orders, primer_pair_editor, primer_pair_model, primer_pair_storage, primer_picker, primer_reference_capture_view, primers |
| `viewer/sequence` | 11 | annotation_geometry, annotation_layout, annotation_table, annotation_visibility, editors, feature_tools, search_view, selection_toolbar, sequence_edit_transition, sequence_geometry, sequence_view |
| `viewer/map` | 3 | glyph_renderer, map_labels, renderer |
| `viewer/alignment` | 13 | alignment_index, alignment_view, compare, dotplot, dotplot_projection, pairwise_context_view, pairwise_export_view, pairwise_pages, pairwise_result_view, pairwise_search_view, pairwise_source_view, sequence_similarity_model, sequence_similarity_view |
| `biology` | 11 | biology, biology_selection, biology_source_view, codon_design, methylation, protein_export, protein_liability_model, protein_liability_view, protein_properties, rna_model, rna_view |
| `biology/antibody` | 6 | antibody_model, antibody_overlay, antibody_overlay_draw, antibody_view, alignment_antibody_model, alignment_antibody_view |
| `biology/structure` | 3 | structure_camera, structure_renderer, structure_view |
| `genomes` | 5 | genome_store, genome_view, genomes, variant_source_library, variants_view |
| `reads` | 11 | contigs, read_alignments_view, read_checks_view, read_review_capture, read_review_history, read_review_view, read_source_library, scf_export_model, scf_export_view, trace_validation, trace_view |
| `protocol` | 8 | plan_editing_host, plan_editing_sources, plan_editing_view, plan_processing_model, protocol_layout_model, protocol_layout_view, protocol_quantities_model, protocol_quantities_view |
| `history` | 9 | history_graph, history_view, run_capture, run_event_model, run_event_view, run_origin_model, run_origin_view, saved_run_review_model, saved_run_review_view |
| `library` | 13 | catalog_backup_model, catalog_backup_view, library_bindings, library_capture, library_client, library_discussion_client, library_discussion_model, library_discussion_view, library_documents, library_view, remote_library_view, remote_removal_model, remote_removal_view |
| `io` | 9 | map_export, map_export_model, oligo_import, sbol_exchange, sequence_compat_export, sequence_compat_import, sequence_export_model, source_gateway_client, source_import |
| `workspace` | 9 | clipboard, clipboard_model, document_explorer, document_tabs_model, reference_transfer_view, workspace_commit, workspace_layout, workspace_preferences, workspace_view |
| `platform` | 10 | commit_lease, desktop_rpc, desktop_shell, desktop_sync_client, native_publication_model, native_publication_view, storage, storage_browser, transient_menus, worker |
| `ai` | 3 | web_mcp, web_mcp_assembly, web_mcp_design |
| `sources` | 1 | source_identity |
| root | 2 | app, main |

This is a physical starting allocation, not a demand to retain every current file boundary. In particular:

- Replace root `app.rs` with an `app/` module: shell component, composition/context wiring, menu/modal router, and demonstration data. Target 4–6 files including `mod.rs`; keep the shell comfortably below 1,000 lines.
- Put session snapshot/restoration/guard/mutation ownership in `workspace/session.rs` or a small `workspace/session/` subtree if needed, rather than moving thousands of root callbacks to another god module. Put feature-specific controllers beside their feature, not all under `app/`.
- Within `assembly/draft`, colocate `Draft`, its transformations, and declaration identity; bounded serialization helpers may be a private codec module. Place its nine cohesive test modules under `assembly/draft/tests/`. This leaves the draft directory below 20 even with a codec and controller.
- Add `io/download.rs` for exports. Move `protein_export` there if its only responsibility remains file export; its initial table assignment groups existing files by present feature responsibility, so this is an explicit refinement, not a second copy.
- Do not keep physical `#[path]` siblings at the old root. Declare their actual logical children at their new paths.
- Module declarations contain real children only. No root facade reexports, aliases, old-path wrappers, duplicate definitions, or `include!` fragments.

## Ordered implementation batches

1. **Establish owners before relocation.** Review the uncommitted SBOL changes with their author, inventory every `crate::app` consumer, settle session/document/selection ownership, and move the export service. Convert the reviewed owner methods above. Keep serialized draft fields, versions, source bindings, resource bounds, and scientific behavior unchanged for this refactor. Introduce no historical decoder.
2. **Extract `App` by feature responsibility.** Session restore/save and document mutation first; assembly and PCR orchestration next; library/read/protocol controllers next; import dispatch and export dialogs next; shell/menu/modal composition last. Controllers own their feature signals and guarded transitions; Leptos components receive typed owners or narrow callbacks. Do not pass a new universal context containing all 124 signals to every child. Do not weaken guards to simplify signatures.
3. **Complete the physical module move as one shared integration batch.** Workers can own disjoint assembly, viewer/biology, and workspace/platform/library trees. A single integrator owns `main.rs`, cross-tree imports, cfg gates, and fixtures. Workers do not individually compile. Update direct consumers, test paths, worker/binary configuration, and build includes in the same batch; delete superseded paths.
4. **Finish draft cohesion.** Move draft identity back to its owner, remove view/model facade imports, extract bounded codec helpers and named tests, and consolidate tiny single-purpose fragments where they only service one owner. Keep genuine UI/scientific evidence boundaries and explicit revision-guarded commit modules.
5. **Verify the completed batch once, collect failures together, repair, and recheck only affected selections.** Run the native web harness and WASM target check using the shared existing cache and consistent package/features. Build selected integration harnesses and the web binary together if both are needed. Run focused existing tests for draft current-format rejection/bounds, exact-event selection, stale proposal commits, tab reconciliation, guarded persistence, and source/export cancellation. No browser/application checks during these edits.
6. **Add a cheap structural guard and review boundaries.** Enforce the absolute limits after meaningful owner splits. Review crossing visibility, module cycles, and remaining owner fragmentation; avoid turning every child into `pub` merely to make relocation compile.

## Refactoring helper assessment

Read `/home/mahogny/github/claude/rust_refactor/README.md`, `docs/MODULE_REFACTOR_WORKFLOW.md`, and `docs/MODULE_SIZE_CHECK.md`. Useful optional discovery: `struct-impl-scatter-stats`, `struct-module-stats`, `module-coupling-stats`, and `to-oop-stats`. UI discovery needs `--target wasm32-unknown-unknown` because many owners are cfg gated. No discovery commands were executed during this audit; the findings above come from direct source inspection.

The helper workflow explicitly says its verification uses the release profile; its automated writing commands may check/rollback after individual moves and describe temporary compatibility imports/delegators. Those defaults conflict with this repository's instructions. Use read-only reports/dry-run plans selectively; do not adopt delegators, compatibility imports, per-file release builds, or an enforced green intermediate step. `apply-batch` is only appropriate if its verification can match the repository's approved shared-cache native/WASM workflow; otherwise apply reviewed edits directly in a coherent batch.

The existing `module-size-check --max-lines 5000` does not implement the whole user policy: it excludes integration-test trees and some generated files, and does not enforce files per directory or the 3000-line warning. Add a small repository-wide structural check covering every Rust file and counting `mod.rs`/tests/generated sources normally. Report 3000–5000 lines as warnings; reject >5000 and >20 `.rs` files directly in a directory. Do not invent exemptions.

## Acceptance

- Every directory contains at most 20 Rust files; the proposed web layout starts at at most 16, including `mod.rs`.
- No Rust source exceeds 5000 lines. Every 3000+ file is visible in the structural report and has a reviewed cohesion justification or follow-up; the two web offenders are resolved in this batch.
- No production code is inserted between inline test blocks. `Draft`, `MethodsState`, source identity, export artifacts, session transitions, and renderer state each have an identifiable owner and a single canonical internal path.
- No imports from `crate::app` outside application composition; feature services do not depend on shell components. No new old-path facades, aliases, fallback internal decoders, or duplicated plan representations.
- All original physical sources are accounted for by a move, deliberate consolidation, or deletion; the four `#[path]` modules and target-gated declarations still have their intended target coverage.
- Native focused harnesses and WASM compilation pass for the completed batch. Source/event/epoch guards, intended-versus-alternative products, serialization/resource bounds, cancellation, and persistent workspace integrity remain enforced by meaningful existing tests.
- SBOL remains subject to its own pending validation; this audit does not claim it is verified.
