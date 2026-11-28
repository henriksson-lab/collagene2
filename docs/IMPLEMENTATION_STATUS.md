# Collagene 2 implementation status

The implementation goal remains active. The [project plan](PROJECT_PLAN.md) is the
current scope and the [delivery queue](DELIVERY_AUDIT.md#current-implementation-queue)
records remaining feature work. The user's final product review is the primary
acceptance step. Historical checkpoints retain their evidence in Git and the
linked audit records; they are not additional unfinished plans.

Application/scientific code is independently implemented and MIT licensed. Runtime
scientific operations use Rust libraries, without external scientific programs.

## Current source

| Item | Current version |
| --- | --- |
| Worker protocol | 127 |
| Candidate / workspace plan | 31 / 29 |
| Assembly draft / batch / release | 20 / 3 / 2 |
| Scientific / search engine | science-35 / sequence-search-2 |
| Documentary image engine | documentary-image-1 |
| Desktop sync | 5 |
| Native MCP tools | 212 |

Cloning methods, supported external format profiles and platforms are summarized
in the [README](../README.md). Current workflows include target verification,
intended and alternative products, PCR-added bases and internal edits, reusable
materials, saved native history, HTML/PDF protocols, shared plasmid folders and
offline desktop storage/sync. Their detailed contracts are linked from the plan.

The current image batch adds **File → Images…**, bounded JPEG/PNG/TIFF/BigTIFF
original retention, private worker previews, historical document associations,
native workspace save/load and MCP attach/list/export. Images remain documentary
attachments and leave sequence, features, material origins and scientific verdicts
unchanged. Archive encoding is checked before UI attachment installation, with a
fresh workspace/document/lifetime check afterward. Standalone sequence exports
omit workspace attachments. [Image workflow](IMAGE_EVIDENCE_WEB.md).

## Verification and CI

Native and web/worker WASM checks pass. The connected batch executes **302 unique
focused cases**, including the four oversized-source cases found by the preceding
scientific CI, image codec/admission/orientation, native archive/transport, actual
in-memory MCP and UI state models. All 1,449 prior names in the 19 compared harnesses
remain; expanded selected discovery has 1,484 cases in 21 harnesses. Discovery is
not execution. [Current evidence](audits/image-evidence-current-integration.json).

The preceding PNG/catalogue/hosted-MCP/security batch passed 490 focused cases;
its [evidence](audits/features-security-current-integration.json) retains that
scope. Counts from different batches overlap and must not be added together.
All 27 extant NCBI genetic codes are implemented, with explicit context-dependent
stop uncertainty. Exact search/PCR has bounded work; model/plan/archive admission
precedes amplified construction; hosted MCP pins member/lab/workspace lifecycles.
[Security findings](DELIVERY_AUDIT.md#security-findings).

At `5f72160`, desktop CI produced the **Windows portable executable, macOS
installer/app and Linux Debian package**, including packaged-executable stdio MCP
smoke checks. These are unsigned development artifacts, not installation or GUI
acceptance. [Package evidence](audits/desktop-packages-current.json) ·
[Package instructions](DESKTOP_PACKAGING.md).

Scientific CI at that checkpoint found four resource-error classification
regressions: oversized models were rejected early as `Format` instead of
`TooLarge`. Typed scope outcomes now preserve resource status, malformed input
remains `Format`, and the unchanged strict cases pass in the current batch.
Exact-source scientific/frontend/package CI follows publication. The final
File-picker ordering adjustment is explicitly recorded in the current evidence.

No browser, CDP, screenshot, application HTTP or live-provider acceptance was run
in these batches. Such checks remain deferred until in-scope code is complete.

## Public preview

Port **8082** remains pinned to CI frontend `48a82bc` (worker 121). It does not
represent the newer source or image work. Its static server exposes no application
filesystem API. Do not infer current product acceptance from this older preview.

## Remaining feature work

- Integrate a curated host codon-frequency profile. Pinned MG1655 CDS/annotation
  sources are captured; an unregistered bounded derivation owner and seven tests
  are authored. Fresh derivation, a compact runtime profile and UI/MCP selection
  remain. This provides observed preferences, not calibrated expression claims.
- Complete the selected-reference CRISPR bulge/variant extensions and additional
  supported nuclease/scoring profiles recorded in the delivery queue.
- Broaden independent antibody reference/species coverage and biological
  calibration. Conditional human-reference numbering and North CDR regions are
  implemented; unsupported mappings remain explicit.
- Finish connected ownership follow-up with the next relevant feature changes,
  then current-source graphical/performance/platform review. The published image
  layout inventory covers 1,273 Rust files in 223 directories, at most 17 direct
  files per directory, with no 3,000-line warnings or hard-limit violations.
  App is 2,454 lines. [Refactor queue](audits/REFACTOR_PLAN.md).

ELN/inventory/freezers/stocks/barcodes, instrument/automation interfaces,
scheduling/reporting, OpenGel comparison/formats, external Primer3/MSA engine
ports and mobile expansion retain their exclusions or deferrals. Further BAM/BAI
work remains behind cloning features. Broader acceptance is not an invitation to
reopen those exclusions or repeat unchanged passing checks.
