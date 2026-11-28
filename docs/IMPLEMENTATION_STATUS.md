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
| Worker protocol | 131 |
| Candidate / workspace plan | 32 / 30 |
| Assembly draft / batch / release | 21 / 4 / 2 |
| Scientific / search engine | science-36 / sequence-search-2 |
| Construct QC options / engine | 2 / construct-qc-2 |
| Material derivation / digest export | 3 / material-digest-export-3 |
| Documentary image engine | documentary-image-1 |
| Desktop sync | 5 |
| Native MCP tools | 219 |

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

The pinned MG1655 host-preference catalogue is integrated with exclusive uniform,
explicit or reference-host selection, retained count/weight provenance and web/MCP
controls. Fresh Rust derivation matches an independent source tally of all 192
count values. Native/WASM checks and **61 focused cases pass**. Observed codon
preferences do not predict expression. [Host integration evidence](audits/codon-host-current-integration.json).

The user confirmed antibody capabilities remain in scope where supported by
Benchling and Geneious. Parallel primary-source GUI/workflow research is mapping
concrete gaps beyond the existing conditional human-reference tracks.

The restriction-condition batch adds explicit formulations, contextual unit-assay
assessment/common-condition suggestions, guarded digest controls and paired
historical setup evidence. Catalogue schema 2 fixes typed Postcard admission;
named digest setups use version 3. Native/WASM checks and **103 focused cases
pass**. This does not establish cleavage, kinetics, star risk or experimental
success. Broader condition coverage and full cloning-step protocol integration
remain. [Integration evidence](audits/restriction-conditions-current-integration.json).

The F36–F38 batch adds source-bound structural QC, the scoped RFC10 part-format
profile and collection-wide scored oligo footprints through the worker, Analyze
panels and immutable-document MCP tools. Final native/WASM checks and all 50
focused cases pass, including the budget-admission correction. [Integration
evidence](audits/design-checks-current-integration.json).
These are partial implementations: broader QC predicates, other assembly profiles,
wizard/protocol integration and an applicable mismatch thermodynamic model remain.
[Controls and supported scope](DESIGN_CHECKS.md).

The next connected source batch retains explicit reaction choices through cloning
preparations, current verification, native material ancestry and HTML/PDF protocols.
It adds declared element-chain geometry and literal reference-copy requirements,
with shared bounded searches and indeterminate uncertainty. The assembly-standard
catalogue now includes the scoped RFC1000 Level-0 documented-body profile. Final
native/WASM compilation and all 86 focused cases pass. [Integration evidence](audits/reaction-elements-current-integration.json) · [Reaction choices](CLONING_REACTION_CONDITIONS.md).

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

At `2bfcc6e`, desktop CI produced the **Windows portable executable, macOS
installer/app and Linux Debian package**, including packaged-executable stdio MCP
smoke checks. These are unsigned development artifacts, not installation or GUI
acceptance. [Package evidence](audits/desktop-packages-current.json) ·
[Package instructions](DESKTOP_PACKAGING.md).

The resource-error and captured-byte fixes pass scientific CI on all three native
platforms plus web at `9d3c76e`; desktop packages also succeed there. At the host
checkpoint `476127e`, scientific cases and web/package jobs pass, while the final
MCP smoke step still used the removed codon field. Its caller now uses current
Explicit preferences, and its catalogue assertion uses schema-2 evidence. Full
local application/HTTP smoke remains deferred. The corrected caller and catalogue
assertion pass CI at `095fe39`. The first connected F36–F38 paths also pass all
native platforms and web at `2bfcc6e`; current reaction-retention/QC-element/RFC1000
source is outside that CI evidence until publication.
[Cross-platform evidence](audits/cross-platform-fixture-fix-current.json).

No browser, CDP, screenshot, application HTTP or live-provider acceptance was run
in these batches. Such checks remain deferred until in-scope code is complete.

## Public preview

Port **8082** remains pinned to CI frontend `48a82bc` (worker 121). It does not
represent the newer source or image work. Its static server exposes no application
filesystem API. Do not infer current product acceptance from this older preview.

## Remaining feature work

- Complete the selected-reference CRISPR bulge/variant extensions and additional
  supported nuclease/scoring profiles recorded in the delivery queue.
- Extend antibody analysis from the Benchling/Geneious GUI and feature comparison,
  including independently supported scheme/species/chain coverage and calibration.
- Complete restriction condition data/presets. Cloning protocol-step and native
  ancestry retention pass the current verification batch; target-free PCR/digest
  inspection remains explicitly limited to geometry.
- Complete the reviewed F36–F41 additions and queue extensions. The first
  assembly-standard, oligo-binding and construct-QC paths are integrated; broader
  predicates/profiles and wizard/protocol coverage remain. Barcode scope is
  awaiting clarification. [Current partition](PROJECT_PLAN.md#implementation-partition-of-the-expanded-feature-queue).
- Retain the documented advanced contig-layout and real corpus gaps alongside
  current implemented read assembly; do not reopen deferred BAM/BAI work.
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
