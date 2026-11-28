# Current delivery gap audit

This is a source-based triage of current delivery gaps, not a completion
certificate. The full feature matrix and M0–M6 gates remain in
[the current plan](PROJECT_PLAN.md). Browser/application HTTP checks remain deferred
until in-scope implementation is complete.

| Area | Current source evidence | Remaining delivery work |
| --- | --- | --- |
| Saved-material digest continuation | Implemented: origin-aware recipe v2, fresh tracked-strand fragment export, named setup reopening, predicted gel and exact-event sharing closure. Current worker schema 126; closure engine 5. | Native/WASM compilation and 823 affected regression cases pass; cross-platform scientific CI and desktop packaging passed at `83eb6c4`. Current-source GUI/platform acceptance remains deferred; GenBank is a detached reference with explicit physical/history losses. |
| Paired-oligo source libraries | Both ordered source events survive choice projection, immutable batch capture, release, native run records, protocol output and complete plate publication. Current draft 20, batch 3, release 2, worker 126 and sharing closure 5. | The base passed focused checks/529 cases and scientific/desktop CI at `ded0067`; current-source CI and deferred graphical/platform acceptance remain. |
| Protected codon ranges | `codon_design::ConstrainedPolicy` includes protected ranges and allowed codons; its independent checker reports violations | Current platform/workflow acceptance; the older blanket “protected ranges missing” entry is obsolete |
| Plate layouts | `protocol_layout` plus web/MCP capture and CSV/HTML/PDF export | Broader workflow acceptance; this is documentary cloning layout, not inventory or instruments |
| Lab comments/review | Native collaboration store, individual member registry, discussion SDK and shared browser/Tauri remote controls | Exact-source CI, graphical acceptance and deliberate deployment cutover; no automatic offline comment queue |
| F33 collection similarity | `sequence_similarity` plus current binary worker schema 126, event-bound Compare controls and immutable-reference MCP discovery rank explicitly selected DNA/RNA sources | Current-source graphical/platform acceptance; server-wide/indexed similarity is outside this implemented profile |
| F34 specialized biology | Protein properties model an unmodified free chain; structure import/viewing has worker, workspace and MCP paths; `collagene-rna-chemistry` adds source-bound base/sugar/backbone/end declarations, stoichiometry, explicit imported-annotation review, event-bound saved declarations and Biology controls | RNA graphical/platform acceptance, broader independently validated antibody numbering/region coverage, graphical/performance acceptance remain. [Structures](STRUCTURE_IMPORT.md) · [RNA chemistry](RNA_MODIFICATIONS.md) |
| M4 source interoperability | URL adapters, retained-source parsers and native MCP binary sequence/trace import plus XDNA/SCF export exist | Real provider/interoperability evidence remains separate from mocked fixtures |
| M6 installed desktop | Native CI packaging exists for Windows/macOS/Linux | Actual installation, offline graphical use and chosen-remote workflow acceptance on the supported platforms |
| Viewer/performance/accessibility | Existing GPU/layout sources and historical measurements exist | Current-source browser/hardware evidence and the user's visual review; historical screenshots/tests do not prove current acceptance |

This audit found no reason to reopen deferred BAM/BAI work or external Primer3 /
multiple-alignment engine ports. OpenGel, inventory/ELN, instruments, scheduling,
organizational reporting and mobile expansion retain their exclusions/deferrals.

## Previous verified trace and recovery batch

Canonical mutation-intent annotation admission and domestication feature handling,
desktop folder-clone recovery (sync schema 5), and worker trace preparation with
document-bound SCF export (worker 125) pass native/WASM compilation and 111
focused scientific/transport/UI-model cases. Test discovery preserves all 2,540
prior cases, with 39 intended trace-test namespace changes; expanded discovery
contains 2,617 cases in 52 harnesses. Discovery includes a deliberately ignored
live provider test; that test was not executed. The captured-record offline case
passes. Exact-source cross-platform CI and deferred product acceptance remain.
[Detailed integration evidence](audits/refactor-current-integration.json).

## Current connected delivery batch

Portable map/gel SVG/PNG export (worker 126), all 27 extant NCBI genetic codes,
authenticated hosted MCP and the confirmed decoding/search/concurrency corrections
pass native/WASM checks and **490 unique focused cases**. Discovery retains all
2,617 preceding names and exposes 2,675 cases in 54 harnesses; discovery includes
ignored HTTP/provider cases and is not an execution count. The initial execution
caught and corrected dual-stop report validation and four new fixture diagnostic
mistakes. No browser, application HTTP or live-provider acceptance ran.
[Detailed current integration evidence](audits/features-security-current-integration.json).
At `5f72160`, desktop packaging succeeded on all three platforms and produced
the requested portable executable, macOS installer/app and Debian package.
[Exact package evidence](audits/desktop-packages-current.json). Scientific CI
found four resource-error classification regressions in existing contig/CRISPR/
methylation/codon transport cases: admission correctly rejects oversized models,
but native decode must preserve `TooLarge` instead of reporting `Format`. The
current connected batch fixes typed scope propagation; the tests remain strict.
Final product review remains.

## Image and resource-error follow-up

The next connected batch implements documentary JPEG/PNG/TIFF/BigTIFF workspace
attachments, worker127 preparation/reopening, exact-original download and three
MCP tools (212total). Native/WASM checks and **302 unique focused cases pass**.
Compared discovery preserves all1,449 preceding names in19 existing harnesses;
expanded selected discovery has1,484 cases in21 harnesses, including35 new cases.
This is selected discovery, not a claim that all preceding54 harnesses ran.
[Current image integration evidence](audits/image-evidence-current-integration.json).
The four scientific CI failures described above now pass unchanged, with typed
resource status preserved through native decoding. Exact-source CI follows the
publication checkpoint; browser/installed-GUI review remains deferred.

## Security findings

Parallel source audits cover web/HTML boundaries, native authority/storage and
scientific/resource limits. Confirmed findings are fixed in source or included
in the current implementation batch; this is not an all-clear certificate.

| Confirmed issue | Current action |
| --- | --- |
| Stale hosted-MCP DELETE could retire a replacement workspace under the same label. | Exact captured-owner comparison under the pool mutex is implemented; deterministic interleaving regression passes. [Finding](audits/security-hosted-delete-current.md). |
| An initialize request naming the modern MCP version still creates a legacy SDK session, bypassing version-based quota admission. | Every initialize now shares the same handshake mutex and session cap. Pure quota regression passes; application HTTP fixtures remain deliberately ignored. [Native audit](audits/security-native-auth-current.md). |
| Native model, plan and snapshot collections could allocate amplified owned data before validation. | Model-owned seeded decoding, frame-wide aggregate budgets and preallocation collection limits are implemented; matching save/readback admission and bounded streaming ZIP workspace JSON decoding are implemented; focused decoder, archive and save/reopen regressions pass. [Input audit](audits/security-web-inputs-current.md). |
| Exact search/PCR could perform quadratic work on a long concrete near-match or ambiguous input. | One shared KMP/bounded-IUPAC kernel now serves exact and typed searches; PCR shares a finite work budget and never returns a truncated success. Focused kernel and actual native MCP PCR regressions pass. [Scientific audit](audits/security-scientific-current.md). |

These fixes pass the connected native/WASM checks and focused regressions above. This does not establish exhaustive security or deployment acceptance.
A suspected normal desktop folder-generation race was refuted by the existing
atomic repository checkpoint CAS; it does not justify a mutation-schema rewrite.
The dedicated [HTML/DOM source review](audits/security-dom-html-current.md) found no confirmed injection path in the inspected code; this is a bounded review, not a general security guarantee.

## Current implementation queue

This is the feature queue. The user's final product review is the primary
acceptance step. Existing code awaiting graphical/live-provider/platform review
is not an additional missing feature. Keep compiler, scientific correctness and
data-integrity checks proportional to the implementation changes.


The current parallel cloning batch passes native/WASM checks and **349 unique
focused cases**. All 86 previous focused names remain in the compiled inventory;
3,023 discovered names are not all executed cases. Two strict catalogue-pin
failures were corrected without weakening admission. The final MCP wording-only
change is recorded separately from compilation. [Current evidence](audits/cloning-extensions-current-integration.json).

| Priority | Concrete source gap | Required outcome |
| --- | --- | --- |
| Complete | F26 host codon data | Pinned MG1655 runtime profile, independent all-192-count agreement, exclusive preference sources, retained evidence and UI/MCP selection are implemented. Native/WASM checks and 61 focused cases pass; observed preferences do not imply calibrated expression. |
| Complete within declared scope | F19 non-destructive annotation merging | Explicit append/skip/attribution, exact full-record donor transfer, accumulated source contributions and native history/worker/MCP/web review pass the [500-case effective batch](audits/annotation-merge-current-integration.json). Actual library rescans remain idempotent despite target-dependent proposal IDs; dataset reimports preserve their source-derived identities. Partial homology and generic qualifier reconciliation remain unsupported. |
| 2 | F25 bounded CRISPR extensions | The current source integrates bounded full-input DNA/RNA bulges, five primary-source nuclease presets, explicit design/screening PAM policies and retained path columns through core/worker/web/MCP; the shared verification below passes. Supplied variants and appropriately licensed empirical scoring remain; deferred large-reference/BAM work does not gain priority. [Research](audits/crispr-extensions-current.md). |
| 1 | F21 restriction reaction conditions | Exact choices, joint unit-assay assessment/suggestions, worker/MCP/UI and named setup retention are implemented. Current source adds plan/protocol and native-ancestry retention with fresh source/material binding; native/WASM checks and 86 focused cases pass. The current parallel batch adds eight cited host presets and four further formulations, in the checked batch. Hosts retain exact catalogue/source pins and leave chemical calls explicit. Broader formulations and appropriate condition-specific star evidence remain. Target-free PCR/digest inspection states its geometry-only scope. [Reaction workflow](CLONING_REACTION_CONDITIONS.md) · [Preset/condition research](audits/restriction-source-presets-current.md). |
| 3 | F34 antibody competitor coverage | User-confirmed scope: map Benchling and Geneious Biologics GUI/workflows and actual supported analysis, schemes, species/chains, liabilities and exports. Current source adds pinned mouse C57BL/6 IGH-only IMGT/IMGT references, explicit species/locus controls, simultaneous conditional interval tracks and training/omission diagnostics; the shared verification below passes. Broader chain/species/scheme coverage and independent residue-level calibration remain. |
| 4 | F36 construct QC | CDS/reference/register paths are integrated. Current source adds structural element chains and literal copy requirements through the same UI/MCP owner. The current source adds a pinned literal initiation-context comparison and bounded canonical ORF enumeration with explicit annotation coverage; the shared verification below passes. Native edit/assembly ancestry, causal junction checks and scoped sourced vector-class checks are integrated in checkpoint `190ae08`, with native/WASM checks and **2,934 distinct** effective selected cases. Preserve individual evidence and indeterminate results; broader biological calibration remains outside that evidence. |
| 4 | F37 assembly standards | Scoped RFC10, RFC1000 Level-0 and original MoClo 2011 profiles use one independent grammar engine. Current source connects exact-role declarations, multi-part review/hints and retained HTML/PDF advice to the wizard. The combined native/WASM checks and 2,959 effective selected cases pass; broader licensed toolkit data and genuine prepared-part coordinate projection remain. Advice does not establish full-standard or physical compliance. |
| 4 | F20 screening strategy | Native named intended/decoy mixtures, exact digest/PCR hypotheses, ideal size witnesses and bounded panel search are connected through binary worker, MCP and web controls. Native/WASM checks and the 23 screening-specific cases pass within the [shared execution record](audits/connected-design-current-integration.json). No success probability, measured gel resolution or intensity claim. |
| 4 | F28/F17 design versus supplied consensus | Current source replaces rectangular affine storage with one sparse row-grid recurrence and typed adaptive/global search evidence. Native/WASM checks and independent literal/oracle fixtures pass within the shared execution record. Explicit design/consensus roles, retained supplied-evidence semantics and conditional CDS consequences pass the [675-case semantic batch](audits/semantic-consensus-current-integration.json). Exact local/URL FASTA retention passes a separate [700-case batch](audits/fasta-retention-current-integration.json); no raw-read pipeline or BAM/BAI expansion. |
| 4 | F38 collection binding | Bounded ungapped suffix scoring, both strands, circular coordinates, per-oligo completion and supported perfect-duplex Tm are integrated. The current parallel source adds independently implemented, licensed single internal DNA-mismatch parameter replacement and the shared sodium/equal-concentration solver; the shared verification below passes. Ambiguous, multiple, near-end or self-associating unsupported footprints retain unavailable Tm. No PCR-success inference. |
| 5 | F39–F41 and reviewed queue extensions | The current verifier/IVT/base-editor/MoClo batch passes native/WASM checks and 1,860 unique selected cases ([execution evidence](audits/verifier-and-design-current-integration.json)); it implements typed cited findings, indexed outcomes, independent branches and shared exact-PCR work; see the [concrete contract batch](audits/verifier-next-contract-current.md). Then deliver the admitted IVT/editor rules and ranked extensions in the current plan. Preserve existing exclusions and the pending genetic-barcode scope decision. [Implementation partition](PROJECT_PLAN.md#implementation-partition-of-the-expanded-feature-queue). |

The antibody comparison now has inspected primary-source GUI evidence from
[Benchling](audits/benchling-antibody-gui-current.md) and
[Geneious Biologics](audits/geneious-antibody-gui-current.md). The concrete table in
[the current plan](PROJECT_PLAN.md#antibody-gui-and-competitor-coverage-requirements)
retains Martin/AHo, species/chain/domain coverage, compact controls, simultaneous
alignment tracks, liability presentation and region/detail/export work. Existing
human numbering, North regions and focused-row mapping remain implemented.

Hosted MCP, portable SVG/PNG export, the27-table genetic-code catalogue and
image attachments are implemented with their verification scopes above. Their
remaining graphical/live/platform acceptance is not another missing feature.

The captured NCBI/DDBJ J02482.1 records provide real public-format input for the
new offline import/provenance/native-save/GenBank-export and independent HaeIII
reference-length fixture. ENA retrieval timed out; no ENA payload is claimed.
This is not live application SourceClient, CORS or authenticated Addgene evidence.
The explicitly ignored live provider test remains deferred with application HTTP
acceptance. Unsigned installer-package generation is a separate build-artifact
gate; executable CI is not proof of installer creation or installation.

## Ownership follow-up

The current [ownership batch](audits/organization-followup-2026-11-28.md)
implements typed worker execution, complete material-preview owners, narrow
restoration state owners, natural planning/record/library methods and shared
bounded encoding. The paired-oligo base `ded0067` passed scientific CI on all
three native platforms and desktop packaging before these moves.

Focused execution passed 1,745 unique cases across affected scientific, worker,
MCP, sharing and web-model paths. All 1,380 baseline case names match across the
same 14 harnesses; expanded discovery exposes 2,430 cases in 27 harnesses. Final
native/WASM checks pass after removal of obsolete imports and the unused recipe-side
hash-limit constant; the shared encoding owner retains its checked 32 MiB stream
bound. Unchanged passing execution evidence is retained. No
scientific, binary selector or current version changes are intended. Independent
static reviews cover both worker execution and material authority; they do not
replace runtime or graphical acceptance. The ownership checkpoint `ba76722` passed scientific CI on Windows/macOS/Linux and desktop packaging; later source changes retain their own CI gate.

The proposal/pool, natural Design/Draft/TreeJob methods and import-facade batch
already passes native/WASM checks and 495 cases; the
[implementation record](audits/ownership-implementation-2026-11-28.md) retains its
scope. These are completed findings, not a queue to repeat.

The connected [controller/owner batch](audits/controller-ownership-2026-11-28.md)
is now implemented in source: App's verification and persistence callbacks use
existing domain controllers, Preparation/TrackedDuplex own their operations, and
renderer State owns its gesture/draw methods. App is 2,413 lines. Each verification
now advances generation and rejects superseded/disposed completions; complete
PlanState invalidation also clears input-event bindings. Native/WASM checks and
107 focused cases pass, preserving all 2,509 case names / 33 harnesses. Exact-source
CI passed at `7f6951b`; the current North batch retains its own exact-source gate.
This adds no scientific representation,
schema, hash or engine version. The [current refactor queue](audits/REFACTOR_PLAN.md)
retains only connected remaining work. Graphical/platform acceptance stays open.

## Completed connected cloning batches

Saved-material digest export, named setups and gel handoff are implemented in
one current recipe representation. `MaterialOrigin` retains exact Run/Derived
ancestry; fresh reconstruction supplies products and gel lengths. Current named
setup schema 2 and worker schema 125 reject superseded internal encodings.
GenBank exports state lost physical/history evidence in both comments and export
warnings. Exact source, origin, options, workspace epoch and generation guards
remain in place. A reference fork cannot inherit a saved physical origin.

Recognition-only reference inspections retain sites and diagnostics when actual
strand material cannot be established. Their required `source_material` field is
explicitly null; products and fragment chemistry are withheld. Saved-origin
replay failures propagate and cannot become reference hypotheses. Historical
setups remain documentary and can reopen without installed catalog membership.

Native compilation and web/worker WASM checks pass. Focused execution covers 823
unique affected cases: 224 material, 118 restriction, 402 web, 60 MCP and 19
sharing cases. One export assertion was corrected for recognition-only reimports;
all four export cases were rerun, while 220 unchanged material cases retained
their passing result against unchanged production code. The prior gel withholding
assertion passes unchanged. Formatting and layout checks pass; browser/application
HTTP and installed-platform acceptance remain deferred. The refactor checkpoint
separately passed all native platforms and desktop packaging at `416b26c`.

Paired-oligo combinatorial choices now retain an explicitly present nullable
bottom-source binding and one canonical ordered input vector. Twelve target
contributions can use twenty-four physical inputs; bottom operands are not
additional target fragments. Equal DNA and repeated events retain distinct roles.
Capture, release, native run records, plate sharing and protocol output use those
same bindings. Web candidate contexts pin both ordered events and options while
allowing construct-local role namespaces to change without losing a still-valid
explicit choice. Current schemas advance together without fallback decoders.
See the [paired-library contract](PAIRED_OLIGO_LIBRARIES.md).

Integrated native harness compilation and current web/worker WASM checks pass.
Focused execution passed 529 unique cases: 407 web state, 66 assembly, 14 history,
eight material, nine protocol-layout format, one core primer, five qualifier
owner-unit, ten public MCP and nine sharing cases. The revised capture case
preserves row-local preview errors while rejecting an immutable capture with
unavailable declared sources; the new version-refusal boundary retains explicit
native errors. Unchanged passing cases were retained while the revised case and
remaining harnesses completed. Layout, formatting, whitespace and current document
links pass. Current Rust layout remains 1,213 files/210 directories, no violations
or 3,000-line warnings. Exact-source CI remains required; browser/application HTTP
and installed-platform acceptance stay deferred.

Coordinate-bearing qualifier projection now reuses the complete location codec
and a bounded, fully consumed qualifier-wrapper parser. Exact local
`transl_except` and `anticodon` positions follow biological parent order and strand
through extraction, rotation, reverse complement and physical product projection.
Translation phase and suffix-only one-/two-base `TERM` exceptions are checked;
clipped or interrupted claims become bounded source evidence with review required.
Unknown, remote, fuzzy, ordered and unresolved positions retain explicit refusal
boundaries. Identity projections preserve source text. This preserves declarations;
it does not interpret recoding biology or extend the CDS translator.

Manual existing-feature location edits now refuse active coordinate qualifiers
on either old or incoming annotations, including removal-and-move in one request.
Other location edits deactivate translation/phase fields case-insensitively and
retain their source evidence. An explicit review/remove step at the unchanged
location permits fresh reannotation. Separate native replace, destination paste
and annotation-only strand-flip refusal boundaries remain in place.

The native MCP registry now contains 209 tools, including
`edit_reference_document`: nine explicit immutable reference operations use the
same core guards, exact source identity and atomic result admission. Receipts
state reference-only authority; no physical material, cloning verification or
native history is invented. Thin WebMCP permissions do not expand.

Scientific engine `science-33` and affected TOPO/recombination/reference-transfer
engine stamps invalidate stale annotation contexts. No binary layout or plan
schema changes are required. Final native harness compilation and web/worker WASM
checks pass. Focused execution passes 2,152 unique scientific, SDK, sharing and
web-model cases, including all 21 new core/SDK fixtures. After advancing standalone
engine stamps and correcting one obsolete clipped-feature expectation, 366 affected
cases were selected for rerun while unchanged passing results were retained.
Discovery exposes 2,451 cases in 27 harnesses. The 4,096-feature replacement case
also verifies that admission does not double-count the old annotation. Exact-source
CI remains required for this annotation batch. See the
[coordinate qualifier profile](COORDINATE_QUALIFIERS.md) for bounds and unsupported
cases. Browser/application HTTP and installed-platform acceptance remain deferred.

The current GenBank integrity batch implements independent location parsing and
formatting, feature-table framing, complete qualifier parsing and safe feature
output. The rich dependency reader receives no feature data and only validated
CONTIG placeholders. This removes its prefix-consumption, signed-arithmetic and
quote/width hazards without fabricating source positions. Native export validates
field tokens and line widths; original records remain exact. Native/WASM compilation,
formatting and layout checks pass. Focused verification covers 450 unique cases:
166 format, 255 sequence, 20 MCP sequence and nine owner-unit cases. The 19 failed
cases were corrected and rerun, alongside the changed owner checks, without
repeating unrelated passing cases. All 18 new public integrity cases are
discoverable. Exact-commit web and desktop packaging CI passed at `3553cdb`;
scientific CI exposed a shared qualifier-key regression affecting saved primer
`5prime_tail` fields on all three platforms. The current fix admits safe numeric
custom keys while retaining delimiter/injection rejection. The two original
core/MCP primer regression cases and five qualifier-owner cases now pass;
scientific CI and desktop packaging passed for the fix at `ded0067`. Current
ownership changes retain their own exact-source CI gate.
[Profile](GENBANK_FORMAT.md).

## Connected native annotation editing

The [native editing batch](NATIVE_ANNOTATION_EDITING.md) is implemented in source.
Projection, qualifier proofs and INSDC codec now have a dependency-safe model
owner. All callers use canonical paths; old core forwarding entry points are
removed. Native replace preserves endpoint envelopes and supported biological
claims, while disrupted assertions become bounded source evidence. Paste uses
that policy without index-pairing remaining annotations. Annotation-only strand
flips deactivate claims; arbitrary manual relocation remains guarded.

Native editing retains the source-format scale of 100,000 features and 32 MiB of
original GenBank under a separate bounded edit profile; ordinary shared-request
admission is unchanged. Positive fixtures cover 4,097 rich features and an original
record larger than 5 MiB. Invalid Unicode whitespace cannot become silent deletion.
MCP adds two typed operations to its existing reference-edit tool, retaining 209
tools and immutable reference-only publication.

Independent comparison accounts for 54 unchanged moved function bodies and source
strings, with new behavior covered by 19 core and two real SDK fixtures. Existing
native-model expectations now check literal geometry, source-evidence notes and
citation replay instead of the removed blanket refusal. Scientific engine
`science-34` and directly affected standalone edit-engine stamps reject stale
contexts; binary layouts remain unchanged. Native harness compilation and the final web/worker WASM check pass. Focused
execution passes 1,349 distinct cases: 852 core, 79 MCP SDK, 11 model and 407 web.
The late continuity correction replays 622 affected cases, retaining 727 unchanged
passing results. Discovery preserves all 2,463 baseline cases after two intentional
semantic test renames and exposes 2,495 cases in 32 harnesses. Source comparison
and execution evidence is retained in the
[owner record](audits/native-annotation-owner-comparison.json). Layout passes
1,245 Rust files / 215 directories with no violations or size warnings. The
preceding `4b87ac0` passed scientific and desktop CI; exact-source CI remains
required for this batch. Installed-platform and graphical acceptance remain open.

## Connected citation and edit-log batch

The current source supports legal GenBank `(sites)`, bare and prose references,
repeated `bases` after semicolons, and fully deleted citation coverage. The
[contract](GENBANK_CITATIONS.md) describes the supported grammar and strict proof
requirements. `model::genbank_edit` owns the serialized ledger, bounded reference
admission and runtime verified log; root forwarding APIs are removed. The current
seven-field encoding, order, key and version remain unchanged.

Citation mapping tracks exact surviving source bases through current-coordinate
edits. Inserted/replacement DNA receives no source mapping. References with no
surviving coverage, or exact coverage beyond the supported output representation,
retain their bibliography and numbering as explicit source-only evidence. Numeric
input, framing, stale binding, version, replay and resource failures remain errors.
The native feature endpoint envelope is deliberately a different operation.

Scientific engine `science-35` and affected mutation-intent, assembly-edit,
junction-edit, domestication and reference-transfer stamps invalidate stale
contexts. No native binary schema changes are introduced. The final exact-map
revision passes native/WASM compilation and **427 selected cases**, including all
**14 new model/core/actual MCP SDK fixtures**. Discovery preserves the preceding
2,495 cases with one deliberate semantic rename and adds 14 cases, yielding
**2,509 cases / 33 harnesses**. The two initial fixture failures concerned wrapped
whitespace and writer quoting; assertions now compare normalized remarks and
parsed citation identities, including the large portable reference number. Only
those failed cases were rerun after the fixture correction. The
[integration record](audits/citation-integration.json) separates discovery from
execution. The complete ledger/proof owner is 733 lines; Molecule remains with
its implementation in the 527-line model root. The prior annotation
checkpoint `1ca0b26` passed scientific CI and desktop builds. This audit does not
establish installed-platform or graphical acceptance.

## Antibody analysis acceptance scope

Benchling's current [antibody annotation reference](https://help.benchling.com/hc/en-us/articles/9684246077325-Annotate-antibody-specific-information-on-protein-sequences)
distinguishes variable-domain detection, alternative residue numbering, CDR/framework
region definitions and candidate liability tracks on individual proteins and
alignments. Numbering schemes include Kabat, Chothia and IMGT; CDR definitions
also include North. These are separate choices. Basic unmodified-chain protein
properties do not fulfill this feature requirement.

Independent native implementation still needs validated domain/numbering references,
explicit supported chain/species profiles and independent numbering fixtures.
Motif-only guesses must not be presented as validated domain numbering. Candidate
PTM patterns are sequence liabilities, not measured chemical state. The source discrepancy is resolved for the new screening core using the canonical
dipeptide definitions in Lu et al.; the exact profile and limitations are documented
in [protein liabilities](PROTEIN_LIABILITIES.md). Worker, Biology and native MCP integration are implemented. Conditional human
V/J analysis now derives IMGT/Kabat1991/Chothia1997 tracks with independently
selected IMGT/Kabat1991/canonical-chart Chothia regions on proteins and exact imported alignment rows.
Unsupported rule profiles, ambiguous boundaries and incomplete regions remain
explicit. Natural residue-level validation, broader chain/species calibration,
North2011 regions now have primary-source boundaries and independent fixtures;
native/WASM checks and 48 focused cases pass, including 10 new cases. See the
[integration evidence](audits/north-integration.json). Current-source CI remains
required. Broader biology and full application
acceptance remain open.
Existing external-tool, MIT-implementation and provenance requirements remain unchanged.


The [scientific contract](ANTIBODY_NUMBERING.md) documents the single current
source-bound report, all-optimal invariant correspondences, finite published
window profiles, separate numbering/region membership and complete versus
unresolved CDR intervals. Current implementation evidence does not establish
universal automatic domain classification or independent biological numbering
accuracy. The fixed external panel establishes deposited-chain identity only;
it is not residue-position ground truth.


## Current restriction-condition checkpoint

The integrated source passes native and web/worker WASM checks and **103 unique
focused cases**. All707 prior names in14 compared harnesses remain; expanded
selected discovery has1,416 names in22 harnesses. Discovery is not execution.
The F21 owner uses current catalogue schema2, worker129, named setup3 and material
digest export2. Full reaction-condition dataset coverage and cloning-step protocol
integration remain in the queue. The newly authored parallel crates/QC and mouse
reference foundation are outside this proof.
[Current evidence](audits/restriction-conditions-current-integration.json).


## Published checkpoint CI and current parallel source

Published checkpoint `190ae08e7842322adbfae9d887d76dfcaf319bae` passes
[Check run 46178350510](https://github.com/henriksson-lab/collagene2/actions/runs/46178350510),
including scientific jobs on Windows/macOS/Linux and the web job. The two stale
PDF prose assertions from `39d70f4` were corrected while preserving typed target
checks. Desktop packaging's fresh offline frontend and all three platform package jobs
pass at that same head. The [current package inventory](audits/desktop-packages-current.json)
records artifact IDs, byte sizes and SHA256 digests for the Windows portable
executable, macOS installer/app archive and Linux Debian package. This is CI
artifact evidence, not installed or graphical acceptance.
These CI results do not verify the uncommitted F20/F28/F37 integration source.

Vector-class QC, IVT linearization selection and native causal/edit ancestry are integrated. Native/WASM checks and **2,934 distinct selected cases** pass across initial and affected corrective executions. [Exact scope and corrections](audits/causal-history-current-integration.json) retain original failures and distinguish source stages; this is not a claim that every case ran on the final binary. GUI/platform/provider acceptance and the broader feature goal remain open.


The connected screening/affine-alignment/standard-workflow batch passes selected
native/WASM checks and **2,959 distinct effective cases**, including the actual
workspace-codec declaration fix and affected SDK round trips.
[Case/source stages and original failures](audits/connected-design-current-integration.json)
distinguish earlier passing binaries from targeted correction runs. The next
design/consensus semantic owner and adapters were excluded from that checkpoint;
their separate connected verification is recorded below.


The connected design/consensus semantic batch passes native/WASM checks and
**675 focused cases**, with no native rerun after its WASM-only syntax fix.
Exact code pins, explicit provenance-preview completeness, one retained history
record, HTML/PDF evidence and aligned-text loss warnings are connected.
[Source and execution evidence](audits/semantic-consensus-current-integration.json).
Graphical/platform/provider acceptance and the broader project goal remain open.
