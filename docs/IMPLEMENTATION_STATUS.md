# Collagene 2 implementation status

The implementation goal is active. The working application covers substantial parts of M0–M3, initial M4 cloning, gel and export workflows, an initial F26 codon-design workflow and F24 dot plots/reference-free read assembly; the complete [project plan](PROJECT_PLAN.md) is not finished. Application/scientific code is independently implemented and MIT licensed. No installed scientific program, downloader or sidecar is invoked at runtime.

## Current source and public preview

Published6af40c7 uses worker109, candidate30, native plan28, science31, draft18 and
189 MCP tools. There is one current internal representation; unsupported internal
schemas are rejected without compatibility readers. Native TOPO graph records
use schema2 and a single ordered transfer/flap/sealing transcript, freshly
replayed from original substrates. One physical bridge replaces the separate
phase executors. Full reviewed-reaction graph capture retains processing steps,
all components, effective cut-end chemistry and consumed activation through restart.

The published batch also connects generic Run/Derived mismatch-PCR and reciprocal-
extension digest inspectors with exact workspace ownership, stages literal
WebMCP design proposals for visible apply/review, and fixes saved-source inverse-
PCR's missing workspace attestation. Nine browser tools share existing Rust
commands. The integrated batch passes234 distinct focused native cases across53 target/
filter groups, including the52-call post-flap/sealed continuation, root ownership
and endpoint/backbone/nick tampering refusals. The web/worker WASM check and static
WebMCP JavaScript syntax check pass. Exact-source web CI passes. Native CI found two stale tag expectations, corrected
in the new working batch without changing scientific assertions. All three
Windows/macOS/Linux desktop build jobs pass for6af40c7.

Publishedc06cdeb uses worker110 with candidate30, plan28, science31, draft18 and
189 MCP tools. It adds exact Run/Derived domestication preparation and TA/GC
library method/row evidence. Current batch records use schema2; domestication
preparation/review use engine/evidence3. Saved sources retain exact owning-
workspace guards, and tailing-family assessment stays distinct from target/
ligation consistency. All193 selected native/SDK/model cases and the integrated
web/worker WASM check pass; exact-source web CI passes and all three desktop builds pass. Native CI found
a missed sharing-test caller and outdated unsupported-version assertion; those
are corrected and checked in the working batch.

Publishedf65723f uses worker111 for target-free reusable OE PCR material.
The existing native material-derivation install/replay path gains a source-bound
OE preview and selected outer-PCR product. Core execution is shared with OE digest
inspection, and explicit web/MCP review retains all six newly performed primer
uses. All97 distinct selected native/SDK/web-model cases and the integrated web/worker
WASM check pass. Native fixture repairs were checked directly. Exact-source
web CI passes. Linux native CI found one outdated nested-method
enum tag; its scientific assertions are unchanged and the corrected case passes
in the terminal-processing batch. [Contract](OE_REUSABLE_MATERIAL.md).

Published3f5ad10 uses worker112 to retain target-free terminal treatments and
single-base 3′ tails through the existing material derivation lifecycle. Native
chemistry, explicit source/origin guards, modification mappings and inherited
limitations remain authoritative. All48 distinct selected native/SDK/worker/web-model cases and the integrated
web/worker WASM check pass. These include real PCR → tails → phosphorylation
→ native restart → ligation with retained alternatives and HTML instructions;
source markers, positional modifications, internal nick behavior and inherited
limitations are checked independently. Exact-source scientific CI passes on
Windows, macOS and Linux, and web CI passes for3f5ad10.
[Contract](TERMINAL_PROCESSING_MATERIAL.md).

Publishedd5e824c uses worker113, candidate31, native plan29, editor draft19
and science32 for reusable annealed-oligo material. Native material records and
capsules use schema2 with complete role-bound sources/snapshots and plural
receipts. Template methods retain one role; annealing retains both oligos,
including two legitimate uses of one self-complementary event. Source history,
the mixed ancestry ledger, selected sharing, MCP and native archive installation
are updated directly. All239 distinct selected native/SDK/worker/web-model
cases and the integrated web/worker WASM check pass. Fixtures cover literal sticky
adapters, same-event role reuse, reverse-strand annotation mapping, archive replay,
public/private source closure, later independent-target ligation and retained
alternatives/incompleteness. Only failed or initially uncompiled cases were
rerun after fixture repairs. Exact-source cross-platform CI remains open.
[Contract](ANNEALED_OLIGO_MATERIAL.md).

Working source connects exact saved-quantity editing and selected reaction
copying to existing native review/capture. Native schemas and tool count are
unchanged. New model fixtures cover fractional/canonical-unit round trips,
missing/diagnostic declarations, immutable revision and atomic positional source
mapping. All six focused model cases and the integrated web/worker WASM check
pass. Exact-source CI for this batch remains open.
[Workflow](REACTION_QUANTITY_EDITING.md).


README documentation was published separately as9da73ad, with two current
screenshots, compact feature/method/format/platform tables and detailed operating
instructions moved to APPLICATION_GUIDE.md. Browser use was restricted to the
user-requested screenshots; general browser/application acceptance remains deferred.
Public8082 was restored using the successful3f5ad10 CI frontend (worker112).


Published19a1a1c (worker108/candidate29/plan27/science30/draft18) passed118 focused
native cases and the integrated WASM check. Exact-source web CI and all three
desktop build jobs subsequently pass. Full native Windows/macOS/Linux CI found
the same eight failing fixture groups: outdated library request shapes, changed
physical-source diagnostics, and the saved-source inverse-PCR attestation gap.
The current source addresses these without relaxed scientific acceptance; focused
rechecks pass in the current batch. Full graphical and installed-desktop acceptance, and the
complete project goal, remain open.

The prior443f51b native CI on Windows, macOS and Linux found the same five stale
binary selector/size-lane fixture expectations, corrected in f66deff without
relaxing bounds. Its web and all three desktop build jobs passed. Full graphical
and installed-desktop acceptance remain open.

Port8082 serves the immutable successful CI frontend for `3f5ad10` (worker112),
including reusable OE PCR, terminal treatment/tailing, saved-material
domestication, TA/GC libraries, ordered TOPO continuation and staged WebMCP changes. Its14 assets,
SRI references, WASM headers, worker namespace and listening socket were checked
without executing the application. Private library assets/storage are untouched.
Browser and application HTTP acceptance remains deferred until in-scope code is
complete. The worker113 plural-source/annealed-material increment is not yet published there.

Saved-run review now explicitly opens a separate candidate in the existing DAG
editor, preserving original source and target events after document heads advance.
It clears cached assessments and verifies edits afresh. Output/protocol/capture
paths use the same retained independent target; recording creates a new run.
Eighteen focused native cases pass, including saved-PCR editing, phosphorylation,
HTML/PDF, archive restart and physical replay. Editing preferences are session-only;
after restart prepare the saved run again. [Contract](RUN_REVIEW.md).

Documentary96/384-well layouts are integrated in the saved-batch panel, including
explicit planned/released identities, replicates, immutable capture and CSV/HTML/PDF.
MCP supports explicit whole-layout dependency review/projection; ordinary plasmid
sharing omits these sidecars. Exact reaction quantities and combined HTML/PDF companions are now implemented.
Saved-report revision and selected compatible reaction copying are in the
current validated batch; six model cases and the integrated WASM check pass. [Model](PROTOCOL_LAYOUT.md) · [Sharing](PROTOCOL_LAYOUT_TRANSPORT.md).

Saved-template inverse PCR now reuses the existing compiler, worker and MCP tool.
It preserves exact origin events, independently frozen targets and fresh native
reconstruction. Newly synthesized PCR strands are unmodified under the existing
model. Intact circular Run/Derived templates are supported; nicked/unsupported
circular contours remain diagnostic. Inherited incomplete origins retain their
findings. [Native workflow](INVERSE_PCR_MATERIAL.md).

The checkpoint evidence below records earlier implementation stages. Current
schema and public-preview claims are the ones above; historical claims of missing
features do not reopen capabilities implemented in later checkpoints.

Parallel format work now adds selected XDNA v0 reference export, SCF regeneration
from retained real trace signals, and direct MCP annotated text exports. These
pass25 focused native/SDK format cases. The two binary writers have reviewed
worker-backed controls in the existing File/trace menus. Fourteen export
core/protocol cases (including six repeated SCF cases) and seven native UI-model
cases pass, alongside the final integrated web/worker WASM check;
the [format matrix](FORMAT_DELIVERY_MATRIX.md) records existing support and gaps.

Native saved-catalog backup controls now create/list/verify and restore an
exact selected archive to a new managed copy through a separate service lane.
They preserve durable retry identities and leave the working editor state in
its existing recovery scope. Reviewed database disconnection resolves pending
work first, forgets its scoped credentials, then atomically retires only that
remote's bookkeeping while retaining local designs/history. Controlled active
catalog cutover, repeated adoption/rollback and restart-bound activation are
implemented on main. Native fixtures cover exact review, editor-only recovery,
scoped cancellation and preserved catalog branches. Installed/offline graphical
acceptance and broader collaboration remain unfinished.
[Backup controls](NATIVE_CATALOG_BACKUP_CONTROLS.md) ·
[Removal lifecycle](REMOTE_REMOVAL.md).

Sequence editing now preserves exact local caret continuity through typing,
cut and normalized paste, while clearing stale external/undo/restore selections.
GenBank RNA T/U handling and citation replay are corrected; explicit FASTA
alphabet review is connected to file/paste and MCP imports. The committed
target-free digest foundation passes38 focused cases and web/worker WASM checking.
Saved digest opening and mixed Run/Digest preparations are committed. New source
adds target-free exact-PCR and repeated digest derivations, native archive/replay,
atomic primer collection installation, Derived-origin preparation and selected
sharing closure. Ordinary PCR controls explicitly request a reusable native
product; ambiguous/multiple-product and mismatch-PCR previews remain
reference-only. Assembly still checks a separately declared independent target. [Editor](SEQUENCE_EDIT_TRANSITIONS.md) · [RNA](GENBANK_RNA.md) ·
[FASTA](FASTA_ALPHABETS.md) · [Digest foundation](SAVED_MATERIAL_DIGESTS.md).

The resumed format batch adds fresh reviewed SnapGene DNA/protein reference
export, explicit-alphabet GCG/Wisconsin single-sequence exchange and opt-in
lower-strand-only original Collagene XML migration. The File menu, guarded import
reviews and three MCP tools expose these selected profiles. Native originals,
uncertain geometry, unsupported semantics and actual read/write integration are
tracked separately in the [delivery matrix](FORMAT_DELIVERY_MATRIX.md) and
[controls contract](SEQUENCE_COMPATIBILITY_CONTROLS.md). All64 focused format,
MCP, protocol, native-service and UI-model cases pass, alongside native SDK and
web/worker WASM compiler checks. New cross-platform CI is pending; browser
acceptance remains deferred.

The parallel format batch adds annotated GenPept/UniProt, PHYLIP/MSF/NEXUS and
PIR/GDE/SELEX alignment exchange, plus an independent SBOL3.1 structural JSON-LD
crate and component/sequence review. Source data, exact labels/coordinates,
explicit losses and native history remain separate from physical authority.
Known profiles have focused fixtures and native/worker/WASM checks. SBOL sharing
now supports an explicitly reviewed canonical structural subset, a new hash and
an omission receipt while retaining the exact original locally. Native identities
and source hashes remain unchanged. The receiver independently validates the
minimal graph and claim coverage; missing-original receipts remain documentary.
Reachable alternative sequences require explicit choices. Native sharing fixtures
include real run/capsule/PCR reuse. [Subset contract](SBOL_DOCUMENTARY_SUBSETS.md).

Explicit 3′ shortening and unsealed homology annealing now connect to typed
verifier operations, source-bound assembly recipes, worker/MCP previews, draft11
persistence and physical strand/gap/nick review. Standalone resection retains its
actual shortened duplex; joined annealing remains incomplete and cannot supply
reusable assembled material. Focused scientific, sharing and MCP fixtures,
integrated UI/worker WASM checking and the native SDK check pass. All57 focused cases pass, including draft11
roundtrip/migration. The full product goal remains active.

Public8082 now serves the exact CI frontend for `6f3a8c7`, with worker88 hash
`819025a440572207fab6a0f6429e3450754841f03f2f06b1aa0c21027de1ae06`.
Its16 immutable assets total149,678,894 bytes, retaining the immediately previous
worker86 for existing tabs. Static references/integrity and exact process/socket
metadata were verified without browser or application HTTP checks. The preview
exposes no host filesystem API. The separately authenticated loopback8084 preview
is unchanged. Current work on private local saves and saved-template TOPO is ahead
of the published source. Visual acceptance remains deferred until all in-scope
implementation is complete. Earlier sections below describe historical checkpoints.

## Current priority

Finish the cloning application first: editor/maps/digests, PCR and primer handling,
restriction/ligation and overlap assembly, product/junction inspection, deterministic
manual/MCP plan verification, retained history and HTML/PDF instructions, then
remaining cloning methods and practical multi-part/batch workflows. The detailed
[execution queue](PROJECT_PLAN.md#current-execution-order--cloning-first-24-november-2026)
overrides the older milestone/format expansion order. Further BAM/BAI/index work
is postponed until the rest of the software is implemented. Existing modules are
retained; no more alignment/index features or integration are scheduled now.
Browser checks remain deferred and compiler checks run alongside coding.

Three parallel follow-up reviews of official Benchling documentation confirm the
existing core workflow and narrow the remaining UI connections: literal-preserving
primer comparison/selected-hit handoff, local recovery, computed pairwise
search/context/exports and genuine combinatorial source bins. Existing read checks
and review proposals are implemented, despite older gap tables. See
[workflow decisions](BENCHLING_WORKFLOW_DECISIONS.md). This documentation-only
checkpoint adds no runtime/build/browser checks.

The latest two parallel official-documentation reviews found presentation seams,
not missing PCR/export engines. The assembly dock now shows method and actual
per-slot preparation before construct review. A full call-site check confirms that selection-to-primer-editor shortcuts and
selected-reference GC already exist in the separate toolbar; the compact exact-event
provenance overview is now implemented over the existing exact-event viewers.
[Research and decisions](BENCHLING_WORKFLOW_DECISIONS.md).

## Manual cloning-plan composition

New manual plans and explicit candidate editing now connect original inputs and
earlier native products through typed PCR, digest, annealing, end treatment/repair,
tailing and ordered ligation controls. Exact source event choices and the separate
target are preserved; edits create an unverified candidate, clear old result state
and reuse the existing verification/output/protocol route. The core remaps producer
identities on valid graph edits and refuses unsafe ordering/removal atomically.
[Workflow, independent adaptor fixture and limits](MANUAL_PLAN_EDITING.md).

A current source audit found the practical gap: the annealing candidate could not
add a separately digested backbone. All three engines already existed. A parallel
[Gateway audit](GATEWAY_MULTISITE_SCOPE.md) records the broader remaining MultiSite
and physical work without reopening the existing BP/LR reference workflow.
Validation passes: 15 focused native cases (eight graph edits, three complete
adaptor/backbone/protocol/history cases and four exact-event catalog cases), plus
the full application/worker WASM build and formatting checks. Port8082 now serves
an immutable worker67 snapshot: 21 assets /150,564,558 bytes, retaining65/61/59
compatibility workers. The static preview exposes no host filesystem API. No
application browser or HTTP smoke checks ran. A [loadable example](../fixtures/manual-adaptor-demo.md)
provides the separate original inputs and frozen intended target. Older deployment
descriptions below are historical checkpoints.

The preceding pushed checkpoint passed all native macOS/Windows/Linux CI jobs and
WASM ([run](https://github.com/henriksson-lab/collagene2/actions/runs/45689734905)).
Browser checks remain deferred.

## Upstream PCR synthesis and compact assembly provenance

A prior native product can now supply a template to an explicitly new ExactPcr
preparation inside TopoActions. Actual tracked complete aligned unnicked material
is required; a role cannot also feed Prepared correspondence. Template modifications
remain on its own retained product, and new strands inherit no atoms/modification
history. Imported summaries and public standalone source APIs remain strict.
[Contract and remaining material/UI work](TOPO_UPSTREAM_PCR.md).

Historical evidence binds each upstream role to its exact prior retained product
event/revision; fresh replay reconstructs the complete plan. Worker67/candidate15/
plan13/science18 add the capability without changing report or operation fields.
The run panel now shows paged exact sources → typed operation → retained product,
with reference-only/historical labels and guarded exact-version navigation.
66 focused native cases pass, including modified-template PCR composition through
real MCP SDK calls, exact archive/replay and legacy version gates. The application/
worker WASM compiler check, core/MCP/web fixture compiler check and formatting
checks pass.
The preceding pushed checkpoint's CI found a duplicate field in the legacy test
fixture, which prevented that test from compiling; the redundant initializer is
removed here and all seven legacy tests pass. This was a test-source compile
error, not a scientific acceptance threshold. [Prior run](https://github.com/henriksson-lab/collagene2/actions/runs/45686799125).
The previous selection-shortcut gap was corrected after checking the separate
existing toolbar; no duplicate controls were added. Browser checks remain deferred,
and public8082 remains on the earlier worker65 snapshot.

## Ordered source-conserving TOPO processing

Explicit anchored-flap hydrolysis and nick-sealing actions now share one native
reconstruction pipeline. The actual excised fragment remains in the graph; a
per-call partition audit preserves all original source material. Fresh contour and
tracked-duplex bridging feed operation20, later digest, independent target checking,
HTML/PDF and retained run history. Desktop review exposes explicit ordered actions
and actual components before a declaration-only candidate handoff. MCP165 exposes
the same inspection/projection and generic verification evidence.
[Full contract and remaining scope](TOPO_PROCESSING_ACTIONS.md).

Worker66/candidate14/plan12/science17 add a generic typed physical-evidence slot.
All nine positional legacy paths preserve old nonempty processed reports through65.
Historical generic/mixed-run downloads bundle exact bound JSON attachments;
history never restores current chemical authority. 90 focused native tests and
the application/worker WASM check pass; browser checks remain deferred.
Public8082 still serves worker65. The preceding pushed checkpoint passed Linux,
Windows and WASM CI; macOS could not acquire a hosted runner after repeated attempts and was cancelled
before any steps ran, leaving that platform unverified for that checkpoint ([run](https://github.com/henriksson-lab/collagene2/actions/runs/45681960208)).

## Processed TOPO in retained plans and downstream operations

The explicit processed-contour action now opens a declaration-only candidate
using exact original source/target events. Fresh verification reconstructs its
full typed graphs and supplies actual private duplex state to downstream digest
operations. Unknown modifications, nicks, displaced material and other components
remain explicit; imported sequence metadata cannot restore that state. Withheld
contours retain typed evidence and return Incomplete. Original StepReport remains
unchanged; dedicated processed evidence is separate.
[Full contract](TOPO_PROCESSED_PLAN.md).

Retained processed runs atomically add an immutable historical graph sidecar bound
to exact run/plan/engine/source/target/product events. Historical inspection checks
documentary bindings, not current chemical rules; old inner model labels stay
inspectable. Fresh replay supplies current state. GUI worker download returns
exact bytes with source/run guards; MCP163 exposes historical inspection and full
processed verification artifacts. Readable HTML/PDF steps retain actual chemical
declarations and source-coordinate selections, with complete graphs in matching
JSON evidence. No raw graph transcript overwhelms the default PDF.

Worker65/candidate13/plan11/science16 introduce the interfaces. All nine old
report-bearing wire paths have explicit positional converters, including the
shorter schema24–26 nested inverse-PCR/replay reports. Native version gates reject
new semantics under old headers; actual old science15 layouts retain their old
minimum versions. Current-writer assertions were updated without relaxing those
historical gates.

Validation passed:62 unique focused native checks (58 core,4 MCP SDK), application/
worker WASM checking and a final complete Trunk build. Tests include actual sealed
TA→EcoRI material/ends against frozen independent DNA, retained nicks/flaps,
wrong targets, byte-exact historical download, native archive/replay, old-layout
migration and related Golden Gate version/biology regressions. The preceding
c288c98 and db8c125 checkpoints passed every configured native Windows/macOS/Linux
and web CI job ([latest prior run](https://github.com/henriksson-lab/collagene2/actions/runs/45675283403)).

Port8082 serves the new immutable worker65 snapshot:18 assets /125,813,134 bytes,
retaining both previously served worker bundles. Exact process/root/socket and
local asset identities were verified; no browser/HTTP checks ran. The server
still exposes no host filesystem API. Local publication evidence is
`/tmp/collagene2-processed-plan-public-evidence.json`.

Accounted directional flap processing, upstream tracked-source chemistry,
durable graph-only review without material output, broader cloning methods and
remaining delivery/GUI acceptance stay active. The full goal is not complete.

## Fresh post-treatment contour and native duplex material

Processed TOPO review now reconstructs original transfer and explicit ligase
treatment once, then proves a source-seeded contour over the actual treated graph.
The native bridge preserves representable top/bottom, topology, free-end and nick
state; unsupported chemistry withholds material and internal modification calls
remain Unknown. Source-coordinate annotations and complete original/treated
component graphs remain retained. The standalone product document carries
restrictive physical metadata, not a reusable proof of material state.
[Contract and remaining integration](TOPO_PROCESSED_PRODUCT.md).

Worker64/MCP162 append the reconstruction interfaces. Existing plan/candidate/
science layouts are unchanged. The published afe2acf checkpoint passed every
configured CI job ([run](https://github.com/henriksson-lab/collagene2/actions/runs/45671141448)).
The preceding c288c98 checkpoint has passed Linux and web CI; macOS/Windows were
still running at the last observation. Current validation passed40 focused
native/SDK cases (processed/source/annotation/binding fixtures, original projection
regressions and actual MCP routes), plus the corrected application/worker WASM
compiler check. The unrelated four-call byproduct fixture selects its component
by original source seed rather than assuming component ordering. No browser/HTTP
check or full preview rebuild ran.

Typed-plan consumption, retained historical graph/run attachments, HTML/PDF
handoff, accounted directional flap processing and remaining delivery/acceptance
work remain active. The next operation/schema migration inventory is explicit;
no old TOPO report is fabricated and no imported sequence grants material state.

## Explicit downstream TOPO nick sealing

An inspected TOPO orientation now opens a downstream ligase review with actual
source-call nicks, explicit selection and a required ligase/conditions assumption.
Fresh reconstruction performs original preparation/transfers once, then joins
only supported exposed3′OH/5′phosphate nicks across a continuous antiparallel
template. The original and treated graphs retain every source nucleotide,
pairing, released activation, displaced flap and coupled component. Blocked,
unknown and unsupported steps preserve their incoming graph.

The shared physical engine, source adapter, worker63 and MCP161 use the same
request and scientific checks. Original wire/plan/candidate/science layouts are
unchanged. A treated graph cannot automatically replace the original projection,
become a tracked duplex or establish independent-target acceptance.
[Contract and remaining F15 work](TOPO_NICK_SEALING.md).

Five source-backed native checks and three real MCP SDK cases passed, including
binary version/downgrade guards, exact source contexts, conditional two-nick TA
closure and retained directional flap. The corrected application/worker WASM
compiler check passed. Seven independent nick-sealing fixtures and29 existing
physical/projection regressions passed (44 total focused native/SDK cases).
No browser/HTTP check or full preview build ran. Port8082 remains the prior
immutable snapshot. Accounted flap processing, post-treatment contour/duplex
composition, durable graph-only handoff and remaining delivery work stay active.

## Published ligation profiles and directional TOPO preparation

Five Pryor2020 matrices now have exact source/model/formulation bindings and
separate CC-BY attribution. Rust imports the original data, ranks retained fusion
sets and exposes detailed pair evidence through workers, the cloning dock and MCP.
An experiment is explicitly selected; geometry alone never substitutes a physical
enzyme formulation. Linear free ends compete in their actual orientation, without
invented complements or extra intended-join factors. The old supplied-data report
layout remains readable. [Contract](PUBLISHED_LIGATION_PROFILES.md).

Directional CACC and TA share native source-backed preparation compilation.
CACC primer additions require review and Apply. Actual PCR and explicitly declared
activated vector strands feed the existing orientation/physical graph review;
unaccounted flap repair cannot produce an accepted paired contour. F15 graph
composition and broader repair chemistry remain unfinished.
[Preparation contract](TOPO_WEB_WORKFLOW.md).

Two focused parallel reviews of current official Benchling documentation confirm
connected preparation/primer/product review and row-specific errors/warnings.
Dedicated public directional TOPO or BP/LR procedures were not located; this is an
evidence limit, not proof of absence. Recommendations and documented behavior are
separate in the [assembly follow-up](BENCHLING_ASSEMBLY_WORKFLOW_REVIEW.md) and
[method follow-up](research/BENCHLING_CLONING_WORKFLOWS.md).

The preceding c8cf748 checkpoint passed every configured CI job, including native
Windows/macOS/Linux and web ([run](https://github.com/henriksson-lab/collagene2/actions/runs/45665755108)).
Current source uses worker62/MCP160, with candidate12/plan10/science15 unchanged.
Focused scientific and binary checks passed:47 unique native core cases and four
actual MCP SDK cases, plus eight MCP internal checks. These cover original-source
preparation, literal end competitors, independent published numerical examples,
hostile data files, stale bindings and append-only protocol/version gates. The
application/worker WASM compiler check and three native TOPO form checks passed.
No full application bundle was rebuilt this checkpoint.
No browser/HTTP checks run; port8082 still serves the preceding immutable snapshot.
Remaining cloning/delivery work and final GUI acceptance stay active.

## Source/CDS domestication and explicit physical preparation

**Remove internal sites** now connects source inspection, retained/protected sites,
selected coding features, a synonymous proposal, independent full-source checking,
reference review and explicit original-template inverse-PCR preparation. The
sequence-selection action opens the same workflow. The shared constrained codon
search supports positional domains, protected bases and scoped motifs while
preserving the old codon API. Exact search/checker work is charged across genes;
unused successful reservations are refunded.

The checker preserves feature coordinates and qualifiers, phase/strand/code and
initiation/stop semantics, rescans both strands and circular origins, and blocks
unsupported uncertain coverage. Real primers must introduce every declared
substitution. Opening reference intent and opening a physical candidate retain
exact original event ancestry through separate guarded actions; neither restores
positive reaction verification. Existing fresh plan review, actual products,
oligo sheets and HTML/PDF exports use the same preparation.
[Contracts and remaining work](DOMESTICATION.md).

Validation passed:187 unique focused native checks (core59, web126, MCP SDK2),
the WASM compiler check, and a final full Trunk application/worker build after the
aggregate-budget fix. The last build includes worker61/MCP156 with unchanged
candidate12/plan10/science15. Source/feature/catalog/candidate drift, phase/reverse/
compound/origin cases, uncertain CDSs, two-gene budget reuse, annotation retention,
physical primer limits, exact event aliases and protocol/oligo export are covered.
No browser/HTTP interaction check ran; GUI/performance acceptance stays deferred.

Port8082 now serves the new immutable snapshot:15 assets /99,419,439 bytes, with
the previously active worker retained for open tabs. The server still has no host
filesystem API. Local asset identity and the exact listening process were checked.
Earlier preview descriptions below are historical checkpoints.

Previous commit4d21bd5 passed all native Windows/macOS/Linux CI jobs. Its web job
failed at formatting before WASM checking; this checkpoint fixes that formatting
and passes `cargo fmt --all -- --check`. New-commit CI remains separate from these
local results.

Joint overlapping-frame optimization, broader method/translation profiles,
multi-round and other physical preparations, curated empirical profiles/ranking,
remaining cloning work and portable desktop delivery stay active.

## Signed Golden Gate fusion placement and supplied-count scoring

The junction editor now exposes signed placement windows, explicit candidate
materialization and Apply through the existing source/preparation transaction.
Original independent DNA, source slots and earlier literal edits stay unchanged;
physical primers explicitly introduce any transferred prefix/suffix bases. Fixed
manual preparations remain literal. Both boundary directions, reverse sources and
circular closing shifts use the actual PCR/digest engines and directed target
policy. Search limits remain visible and do not establish an optimal or unique
mixture outcome.

Optional supplied 3/4-base directed count matrices retain conditions, provenance,
license and data revisions. Missing observations stay unknown. Exact wide totals
use decimal JSON strings and binary 128-bit counts; the declared conditional model
does not claim commercial optimizer equivalence or experimental efficiency.
Worker60 appends advisory search/materialization/scoring jobs; MCP152 exposes the
same paths. Candidate12, plan10 and science15 remain unchanged.
[Contracts and unfinished F13 work](GOLDEN_GATE_FUSION_DESIGN.md).

Validation passed:153 unique focused native checks (core31, web119, MCP SDK3),
including signed/reverse/circular actual preparation and ligation, retained
manual settings, original intent, bounds, stale selections, worker downgrade
guards, partial count data and exact wide JSON/binary count roundtrips. The final
application/worker WASM compiler check passed. No browser/HTTP checks or new full
Trunk preview build ran for this checkpoint.

The complete goal remains active: broader biological constraints and curated
empirical profiles/ranking, internal-site domestication and physical mutagenic
preparation, richer methods/mixtures and remaining delivery/acceptance still need
work. Port8082 still serves the earlier reaction-matrix snapshot774db29.

## Reviewed library bin editing and exact recovery

The construct library now supports explicit add/remove/reorder proposals with
per-construct primer, junction, preparation and candidate impacts. Changes preserve
stable source identities and raw manual settings. Recovery is staged with the new
active draft and installed in one guarded current/previous transaction; cancellation
uses the existing worker and owned commit-lease pattern.

Exact raw snapshots preserve unsaved working fields separately from saved row
settings. Restoration is another reviewed edit with source-event attestation and
monotonic allocation across current, recovered and retained history. Source DNA
heads and unrelated workspace data remain unchanged. History is readable as data
through two MCP additions (149 tools), without restoring scientific authority.
[Contracts, limits and remaining work](ASSEMBLY_LIBRARY_EDITS.md).

Validation passed:135 focused native checks (recovery7, workspace12, web114,
MCP SDK2), followed by the final application/worker WASM check. Add-bin selection
uses active working settings, and root-owned restoration callbacks survive a dock
closing after commit. Browser checks remain deferred. Port8082 retains the checked
reaction-matrix snapshot774db29; this checkpoint was not republished there.

## Whole-reaction construct matrices and exact candidate identity

Source-matrix constructs can now discover finite Golden Gate pools from prepared
material, whole originals and complete exact-PCR amplicons, then save an explicit
candidate declaration. Exact source events, contribution geometry and independent
intent are pinned separately from DNA equality. The new native operation checks
index and both sequence/physical keys in one fresh enumeration; a target-matching
alternative is never substituted.

The web provides paged raw fragment/path/recut/chemistry evidence, a read-only
reference canvas, complete guarded JSON export and direct shared product/junction/
output review. Ordered-method controls cannot silently replace the selected pool;
returning to ordered assembly is explicit. Draft3 migrates earlier drafts without restoring proof. Worker59,
candidate12, plan10 and science15 preserve typed compatibility; MCP147 includes
source-bound pool discovery before existing batch/release/export operations.
[Contract and remaining boundaries](COMBINATORIAL_POOLS.md).

Validation passed:71 focused core tests, six real MCP SDK fixture groups, all101
native web tests and the final application/worker WASM check. The final full
Trunk build also passed, including direct selected-reaction review/output routing.
The public preview at http://130.239.192.39:8082/ now serves worker59 and includes
the previous source-bin, batch and release checkpoints. The immutable snapshot
contains51 assets /250,201,860 bytes; the server has no host filesystem API. Listener
and local asset identity checks passed; browser/HTTP interaction checks were not
run. Final application acceptance remains deferred until implementation completes.

## Fresh atomic construct outputs and bulk ZIP export

Saved batch declarations now support selected per-row named outputs through a
fresh private native finalizer. Actual staged product/target review precedes an
explicit all-or-nothing workspace publication. Any failed row discards the stage;
historical positive statuses never authorize output creation. Immutable manifests
retain exact source, target, run, product, named construct and primer links.

The web uses guarded cancellable same-worker persistence and owned commit leases.
Read-only export supports selected GenBank/FASTA/primer CSV/TSV and HTML/PDF files
in a bounded ZIP, with complete paged loss/provenance review before download.
Sequence-only exchange remains historical; protocol options freshly check current
independent declarations and share one native report between HTML/PDF. Worker58
adds private release/inspection and archive jobs. Four MCP additions bring the
registry to146 operations. [Contracts and remaining scope](COMBINATORIAL_RELEASE.md).

Validation passed:54 focused core tests, four real MCP SDK fixture groups, all98
native web tests and the final WASM check for both application and worker. The
protocol checks distinguish a rejected out-of-order request from an actual failed
reconstruction: the former preserves an untouched stage, the latter discards it.
No browser checks or public preview refresh were run; port8082 retains schema56.

## Selected batch evaluation, historical capture and MCP

Selected constructs now run sequentially in a private worker with progress,
cancellation and row-local draft/reconstruction/evaluation outcomes. Per-row
manual preparations, method, topology, treatment and independent junction intent
are retained; bad editable rows do not hide valid ones. Historical native records
retain exact source/target/plan/product identities and current-engine evidence
without releasing products or restoring fresh authority on load.

The web saves evidence through a guarded current/previous recovery transaction,
with cancellable same-worker encoding and token-owned commit leases. The dock's
editing controls lock during capture; cancelled old completions cannot unlock
later writes. Saved matrices reopen read-only and replay explicitly. Six MCP tools
provide bounded design paging, inspection, batch capture/list/open and row replay;
the registry is142 tools. Worker57 appends private session begin/advance/capture
while preserving existing tags and native schemas.
[Detailed contracts, limits and remaining output work](COMBINATORIAL_BATCHES.md).

Validation:19 focused core tests, two real MCP SDK groups and all98 native web
tests passed, followed by the final WASM check. This includes the cancellation
lease regression and all prior web fixtures. No browser checks or public preview
refresh were run; port8082 retains the schema56 comparison/primer/digest snapshot.

## Alternative source bins and per-construct edits

The existing assembly dock now captures immutable source alternatives, pages the
Cartesian matrix and opens explicitly selected constructs in the same preparation
and verification route. Each row preserves its raw primer/tail/method/junction
settings. Stable design/bin/choice/construct identities flow into original-source
inputs and proposal/candidate guards. Pending source edits persist but cannot
inherit old choice identity. Draft version2 migrates legacy data, isolates old
readers and attests every inactive source on restoration.

UI-independent Rust models support native state/persistence checks; the typed core
resolves actual historical events and dispatches existing compilers against an
independently declared target. Bulk evaluation/output records, MCP design tools,
release and library topology editing remain delivery work.
[Behavior, limits and remaining scope](COMBINATORIAL_SOURCE_DESIGN.md).

Eight core design fixtures and all93 native web tests passed, including eight new
web state/identity groups. The WASM compile passed. No browser checks or public
preview refresh were run at this checkpoint; port8082 retains the preceding
schema56 comparison/primer/digest snapshot.

## Computed comparison search, source context and portable review

The event-bound pairwise viewer now searches displayed rows with literal or
compatible-IUPAC policies, reports all returned hits through pages and navigates
their original source ranges. Annotation and full-source CDS context project onto
the same600-column page with circular/reverse/gap-aware mappings and explicit
unsupported/budget notes. Core preparation and search run in cancellable workers;
page changes reuse the context cache.

Direct aligned FASTA/Clustal and complete HTML/PDF review reuse existing exchange
and PDF primitives, preserving gaps/orientation and disclosing format losses.
PDF base columns have explicit fixed cell positions; rich review distinguishes
source events, original coordinates and historical/fresh status from experimental
verification. Worker56 appends three jobs/results with the ordinary32 MiB cap;
native archive/candidate/plan/scientific-engine schemas stay unchanged.
[Presentation, limits and guards](PAIRWISE_PRESENTATION.md).

Validation passed: 40 focused native tests (eight search/context, five export,
seven binary protocol, twelve existing alignment-exchange and eight existing PDF
cases) and the combined Trunk app/worker build. No browser or HTTP smoke checks
were run. Port8082 now serves the immutable schema56 snapshot:48 assets,
228,876,923 bytes, retaining older versioned workers and exposing no filesystem
API or writes. Index SHA-256:
`429aec162ab5d2c3457394fb2c7d0b6b65342e81d37a329bc081a9cb2910aa1a`.

## Existing-oligo comparison and selected binding handoff

The shared pair editor can preview and accept a new directional reference
comparison while preserving exact oligo text, including intentional mismatches.
Single-oligo results page every returned hit and transfer an explicitly chosen
placement into either pair slot with actual source/binding-fact checks. Navigation
highlights every circular-origin part. Local edit-domain, review-policy and rerun
controls reuse the existing form without automatic threshold changes. Epoch,
event, form and calculation guards protect actions; pair Apply stays separate.
No new parser, scientific engine or storage/worker schema is introduced.
[Semantics and limits](PRIMER_REFERENCE_RECOVERY.md).

All 85 native web tests and the WASM compile passed, including five new native
capture-contract fixture groups.
Browser acceptance and preview refresh remain deferred at this checkpoint.

## Primer endpoint and tail-boundary previews

The shared pair editor now previews and explicitly accepts endpoint or tail-boundary
changes against an exact retained source event. Literal edits survive coordinate
overlap; new reference context, removed calls and ambiguity remain visible.
Tail repartition preserves the complete ordered oligo. Circular-origin domains
retain normalized/unrolled geometry through capture, storage and overlays, while
wrapped direct base/mutation editing stays disabled. Source/pair/epoch/control
guards plus fresh reconstruction protect acceptance; applying the pair and running
PCR remain separate. Pair attachment version2 reads/migrates legacy version1 data.
[Editing and migration semantics](PRIMER_DOMAIN_EDITS.md).

Validation: all 80 native web tests and the WASM compile passed. No browser checks
or preview publication were run at this checkpoint.

## Named digest setups and predicted gels

Digest setups now retain exact source events, full catalog pins/options and
historical inspection evidence in native workspaces. Current/saved setup actions
send typed recipes to predicted gels without enzyme retyping. Historical review
never restores fresh evidence; explicit settings reuse and source forks retain
strict identity and require a new inspection. Gel lanes freshly resolve/inspect
their retained source, preserve existing limits and withhold incomplete products.
Multiple lanes survive inspector navigation, while workspace restoration clears
old source selections. Pinned JSON retains exact submitted recipes and lane data.
[Semantics, compatibility and limits](NAMED_DIGESTS.md).

Validation passed: 35 focused core tests, 65 native web tests, two real MCP gel
tool tests and the web WASM compile. Coverage includes unrelated 10,000-feature
documents, retained historical events, nonstarter/staggered digestion, legacy
compatibility and binary limits. No browser checks or preview refresh were run
for this checkpoint.

## Synchronized PCR source overlays

Applied captured domains and explicitly selected exact/mismatch PCR products now
share transient source overlays in sequence and circular/linear maps. Directional
primer placements, amplified intervals and diagnostic source bases are separate
from synthetic tails. Captured-reference comparisons remain visibly weaker than
computed binding evidence; incomplete searches and ambiguous role assignments stay
explicit. Cross-document source inspection opens a read-only snapshot without
retargeting an independent target. Owner/event/epoch guards prevent stale overlays,
and a separate context retains the editor's cursor/scroll identity.
[Projection semantics and navigation](PCR_OVERLAYS.md).

All 65 native web tests passed, including seven focused PCR projection fixtures.
The combined Trunk web/worker build passed and port 8082 was refreshed with an immutable
schema-54 snapshot at that checkpoint, including annotation, translation, enzyme-selection
and historical-source work. The server loaded 45 assets (210,781,733 bytes) into RAM,
retains old versioned workers, and exposes no filesystem API or writes. No browser
or HTTP smoke checks were run. Index SHA-256:
`982729a5963514a86eab2886a051f024afcfee97094a70b17a70ea736f25c206`.

## Selected translation and linked protein outputs

The biology inspector now translates an exact selected location or whole reference,
with explicit frame/orientation/initiation policy, and retains the separate
qualifier-directed CDS route. Complete protein review is paged and exports FASTA.
Fresh source reconstruction precedes named protein creation or a review-required
CDS proposal. Protein documents retain exact nucleotide-event ancestry, selection,
codon options and warnings; the origin action opens the retained historical source
without changing its head. ORF paging and selection retain all compound spans and
biological orientation. [Semantics and source linkage](TRANSLATION_REFERENCE.md).

Focused native validation passed: 11 translation/output fixtures, seven binary
protocol/version/budget fixtures and 58 native web tests, including compound
selection and historical-origin ownership. The final web WASM compile passed.
Browser checks and preview publication were not performed for this checkpoint.

## Annotation management

Annotations now have a searchable paged table, source-bound multi-selection and
one-commit type/strand/color edits. Native batch operations preserve unfamiliar
qualifiers and compound/fuzzy locations; strand changes retain stale coding
assertions as source evidence and withhold unsupported coordinate reinterpretation.

Independent show/hide preferences filter map/sequence geometry, labels, hits and
lane occupancy without altering annotations or scientific reports. Exact-event
preferences persist in native workspaces and survive supported batch edits.
Literal RGB colors, including safe imported ApE strand colors, are shared across
renderers. [Controls, persistence and current resource bounds](ANNOTATION_MANAGEMENT.md).

Native checks passed: ten core scientific fixtures and 55 native web tests,
including exact-event display/archive and layout cases. The final web WASM compile
passed. Browser checks and preview publication remain deferred at this checkpoint.

## Retained Golden Gate reaction settings

The assembly draft now retains prepared-pool/reaction mode, each role's substrate
choice, topology/orientation policy and finite search limits. Choices follow roles
when reordered/removed; new roles default to their current prepared fragment.
Older drafts initialize the existing defaults using their topology and row count.
Receipts remain transient, and reopened settings require explicit fresh inspection.
Options edit directly in shared state without remounting their panel; stale worker
results and candidate selections remain guarded.
[Persistence and invalidation](GOLDEN_GATE_POOL_SETTINGS.md).

Focused native validation covers editable/runnable state, hostile or oversized
decoding, evidence-field rejection and agreement with existing native reaction
limits. A read-only review checked reactive mounts, source/receipt guards and
role remapping before compilation. All 51 native web tests and the final WASM
compile passed. No browser checks or preview update were run.

## Parallel workflow review follow-up

Three further official-document/source reviews cover
[assembly](BENCHLING_ASSEMBLY_WORKFLOW_CHECK.md),
[PCR/primers](BENCHLING_PCR_WORKFLOW_CHECK.md) and
[editor/supporting work](BENCHLING_EDITOR_FOLLOWUP_CHECK.md).
They distinguish implemented capabilities from actual remaining gaps, and documented
Benchling interaction from recommendations. The plan now prioritizes reaction-setting
persistence, connected source/preparation/junction review, synchronized primer/product
overlays, annotation operations and genuine combinatorial source bins within the
existing cloning-first queue. No local browser checks were performed.

## Restriction selection and named sets

The shared restriction scan now supports pinned inside/outside/crossing filters
using either recognition bases or both strand-cut boundaries. Sequence, map and
picker use the same scoped evidence. Site navigation preserves the pinned region;
source changes invalidate it. Possible recognition, incomplete scans and unavailable
cut geometry cannot become exact zero/one-site claims after filtering.

Named enzyme sets persist in native workspaces with ordered choices, stable catalog
IDs and pinned geometry. Missing/changed definitions prevent loading; harmless
catalog changes remain visible. Saving/deleting sets changes attachments only.
[Semantics and persistence](RESTRICTION_SELECTION.md).

Validation: 20 focused native region/preset fixtures and 43 native web tests passed;
the web WASM compile passed in the shared cache. No browser checks or preview
publication were performed for this checkpoint.

## Annotated reference transfer

The editor's internal reference clipboard copies exact selected regions or
biological join/complement selections with annotation projection. It works without
system clipboard permission and survives source navigation while original events
remain retained. Paste previews all product bases/annotations and warnings in
pages, then freshly reconstructs before applying one source-bound undoable edit.
Both source and previous destination events are retained as computational parents;
plain OS clipboard actions remain separate.

Projection uses the original source directly against final destination topology.
Fixtures distinguish a linear preview's broken-origin gene from a restored circular
gene, and from a genuinely interrupted gene in a larger destination. Source citation
numbers, raw annotation/header evidence and unmapped features remain scoped evidence.
Current physical/trace/FASTQ descriptors do not become invented prepared material.
Native reference editing also handles translation/unsupported coordinate-qualifier
case variants consistently, preserving untouched assertions.

Worker53 adds copy120/paste121 and results112/113 under ordinary 32 MiB limits.
Candidate11, plan9, science-14, archive1 and MCP136 remain unchanged. Sixty-four
native cases passed: six transfer groups, two binary groups, four multi-parent
history cases, seven existing history cases, two qualifier-case fixtures and 43
native web cases. The final web/worker WASM compiler check passed. Browser checks
and full preview rebuilding were not run; port8082 still serves the last immutable
primer/mismatch publication. The full product objective remains unfinished.
[Transfer contract](REFERENCE_TRANSFER.md).

## Retained construct comparison

The comparison panel now keeps a reference pinned while navigating reference/query
differences. Stored origin/orientation supports global, local and query-in-reference
alignment; explicit whole-construct policy reuses native bounded circular-origin
and reverse-query search. Candidate limits, best-score scope, uncertain calls,
coverage and unaligned flanks remain visible. No alignment is a sample identity
or physical cloning verification claim.

All alignment columns, differences and CIGAR text are paged, with discrepancy
stepping and number/column jumps. Native-coordinate navigation supports insertion
and deletion carets, reverse queries and two-part circular-origin spans. Current
one-span events navigate the main document without retargeting the comparison;
historical/origin-spanning positions use the shared read-only reference canvas.

Comparison records use bounded `pairwise-review-v1` native attachments with exact
reference/query events and options. Saving stages only the new attachment after
source/path consistency checks; it does not run a second dynamic program or change
heads. Historical reopening restores no current/fresh flag; explicit Compute
reanalyzes the pinned versions. Native archive loading checks path letters,
statistics, differences, CIGAR and path score, without treating an imported
optimality assertion as current evidence.

Worker52 adds job119/result111 under the ordinary 32 MiB limit. Candidate11, plan9,
science-14, archive1 and MCP136 remain unchanged. Fifty native cases passed: four
core review groups, three archive groups and 43 native web cases, including two
gap/late-page display groups and a supported unrelated-candidate history regression.
The final web/worker WASM compiler check passed. Browser checks
and full preview rebuilding remain deferred for this checkpoint; port8082 retains
the last successful primer/mismatch publication. Integrated trace review and the
broader cloning/editor/native delivery scope remain unfinished.
[Comparison contract](PAIRWISE_REVIEW.md).

## Exact versions from saved cloning runs

Saved runs now offer exact target/source/product-version navigation. The host
re-resolves the requested run role and document/event/revision binding, then opens
a frozen read-only sequence/map view using the shared reference canvas. Opening
and returning to the run do not rewind heads, create documents or restore a
positive assessment. The current-head notice is explicitly the state at opening.

Fork editable copy is separate: the host rechecks epoch, displayed request and run
membership, calls the existing exact-event fork operation and validates capacity
before publishing only the new document/active selection. Original heads, recorded
assessments, candidate data and unfinished primer/assembly drafts are preserved.
Document-ID syntax is not imposed on otherwise valid Unicode/spaced plan input IDs.

Forty native web cases passed, including three new real recorded-run fixture
groups for all roles, later edit/revert events with equal DNA, alias/missing-event
rejection and explicit fork ancestry/rollback. The final web/analysis-worker WASM
check passed. Worker51, candidate11, plan9, science-14, archive1 and MCP136 remain
unchanged. No browser checks or full preview rebuild ran for this connection;
port8082 retains the last successful primer/mismatch snapshot below. Ordinary
pairwise evidence/navigation, integrated trace review, remaining cloning/editor
work and portable delivery remain active. [Run UI contract](RUN_LINEAGE_UI.md).

## Primer-led reference outputs

Mismatch-PCR now retains the submitted named pair and exact source/viewed-reference
events with the current form and worker evidence. Product, DNA, substitution and
finding pages replace the former first-50-only preview. Users select and name one
geometry, then choose reference and/or oligo outputs; fresh worker reconstruction
checks the original preview revision before atomic installation with source-event
ancestry and a compact search/geometry receipt. Capped, ambiguous and unsupported
search states remain explicit. This route asserts no independent-target agreement
or physical PCR success.

Single analyzed oligos and individual ranked pairs now export CSV/TSV or save an
oligo collection. Source/event, literal form, epoch and generation guards prevent
stale actions. The shared generic order constructor retains one row per declared
oligo, full named tails/domains and existing evidence-only/size policies; exact-PCR
pair tables delegate to it. Pair order actions are collapsed in their result row.

Worker51 adds job118/result110 under the ordinary 32 MiB limit. Candidate11,
plan9, science-14, archive1 and MCP136 remain unchanged. Sixty-two native cases
passed: four selected-mismatch-output groups, two generic order groups, four
existing selected-PCR groups, 15 oligo-format regressions and 37 native web cases.
The full Trunk web/worker build passed and the immutable preview was published
on port8082. [Output contract](PCR_OUTPUTS.md).

### Current public preview checkpoint

The successful primer/mismatch build replaces the earlier diagnostics/pair-draft
snapshot at `http://130.239.192.39:8082/`. It includes the intervening whole-reaction
Golden Gate, exact/inverse/mismatch-PCR and primer-output connections. The server
loads only the staged compiled assets into RAM and exposes no filesystem API,
uploads or writes. Old versioned worker bundles remain available for open tabs.

- Immutable root: `/tmp/collagene2-primer-mismatch-public-dist`
- Assets: 42 files, 192,291,826 bytes including retained older workers
- Index SHA256: `ebbbde7ff68c0f74b49a9f7ab4671ebc40ff6bc41fb308ac5ef7b382b4c993b7`
- Active worker bundle: `fcfbaa00f4db4cf8278aac24f2dae380eb9971dcdc5d0ba49da159d936ef7d41`
- Publication evidence: `/tmp/collagene2-primer-mismatch-public-evidence.json`
- Native log: `/tmp/collagene2-primer-mismatch-native.log`
- Full build log: `/tmp/collagene2-primer-mismatch-build.log`

The replacement listener was confirmed on `0.0.0.0:8082`. No browser, HTTP smoke,
screenshot or GUI acceptance checks were run. This proves build/publication, not
finished GUI readiness or completion of the full product objective. Earlier
checkpoint statements about the old preview record their publication state then.

## Explicit ordinary and inverse PCR outputs

Follow-up: inverse-PCR now forwards the original source event, independent target
event and epoch into candidate review. The main plan panel reuses the existing
actual-product view and finalizer, including named constructs, optional primers
and prepared-reference visibility, CSV/TSV and HTML/PDF. Target-event receipts
guard later actions; explicit historical source events remain usable after head
advancement without substituting matching DNA in another document. Failed-final
plans retain diagnostic protocol export, without opening a preceding preparation
as a final product. Existing-target staging creates no synthetic assembly intent.

Forty-five native cases passed: six existing cloning-output cases, two end-to-end
inverse-PCR output fixtures, and 37 native web cases including four new exact-event
staging cases. No new scientific engine, runtime dependency, worker schema or MCP
tool is introduced by this connection. The web application and analysis worker
passed the final WASM compiler check; no browser checks were run.
[Output contract](CLONING_OUTPUTS.md).

Ordinary exact-PCR no longer opens every returned product immediately. The review
pages every geometry and selected DNA, shows actual oligo roles, retains names and
tails, and requires a choice when several geometries exist. Users choose a named
reference preview, optional oligo collection, or export the named pair as CSV/TSV.
Output creation freshly reconstructs the selected preview in the worker and checks
its original revision; the host checks the source event, literal draft and epoch
before staging a complete validated workspace. Output history and a compact
reference-provenance attachment retain the original event and selected geometry.
No independently verified target or physical chemistry is asserted by this route.

Worker50 adds job117/result109; existing candidate/plan/science/archive versions
and the MCP registry are unchanged. Twenty-four focused native cases passed,
covering selected alternatives, annotations, stale inputs, complete named orders,
optional outputs, IUPAC/long tails and old-schema binary rejection, alongside
existing PCR-document and oligo-format regression cases. The web application and
analysis worker also passed the WASM compiler check. Browser acceptance remains
deferred. The wider product goal remains unfinished; saved-run event navigation
and the remaining cloning/editor/delivery work stay active.
[Output contract](PCR_OUTPUTS.md).

## Coherent assembly junction edits

The persistent assembly workspace can trim the two oriented source cores around a
junction and insert literal DNA. Intended DNA is rebuilt independently from those
declarations. The proposal analyzes original-source PCR primers with the complete
insert-plus-homology tails present before ranking, including reverse sources and
circular closure. Explicit Apply updates the declared edits, expected reference,
all named PCR primers and overlap vector together. Original source slots remain
unchanged; fresh method reconstruction still checks the recipe against the target.

This first route deliberately proposes every row as PCR-prepared overlap assembly;
it does not silently alter a mixed digest/ligation recipe. The review shows all
oligos before application. Applied edits survive native draft reload; analysis does
not. Undo restores the prior draft only while no later declaration or primer change
would be overwritten, and never restores a positive result. Source navigation uses
actual adjusted original coordinates; synthetic contributions are labelled.

Worker44 adds independent edited-intent and primer-design jobs109/110 with
results101/102. MCP has130 tools. Candidate11/plan9/science-14/archive1 remain
unchanged. Sixteen focused core/binary/MCP/staging cases passed, including the four
existing region-staging regressions. The full 130-tool stdio smoke and final
web/worker WASM compiler check also passed. Browser checks remain deferred. Selective
primer preservation across mixed preparations, Golden Gate junction editing and
broader bulk-cloning flows remain pending. [Edit model](ASSEMBLY_JUNCTION_EDIT.md) ·
[Exact-event staging](ASSEMBLY_JUNCTION_STAGING.md).

## Explicit sequencing review and reference proposals

Saved construct checks now open a review workspace with a required fresh replay.
The operator can retain a reference call, leave evidence unresolved, or nominate an
exact eligible observed base, deletion or insertion. The choice retains source
read/call/boundary/allele identity; optional notes default to factual choice
provenance. Changed placement or allele meaning is refused, and unselected
positions remain unreviewed. Conflicting observations stay in the full report.

Preview shows literal original-coordinate edits and an optional revised reference.
Save repeats fresh reconstruction and commits decisions plus the optional derived
document in one recovery transaction. Original reference/read heads and raw quality
are unchanged. New references retain all original source events as computational
parents, remain restricted proposals, and require a new read check. Retain/unresolved
only saves decisions without manufacturing a new reference. Edited earlier outputs
are preserved; retries reuse only the matching proposal and historical record.

Historical decision records reopen as editable choices without restored evidence.
Worker43 adds preview/commit jobs107/108 and results99/100; MCP has128 tools.
Candidate11/plan9/science-14/archive1 remain unchanged. Nine focused core/binary/SDK
cases pass, along with the full 128-tool stdio smoke and final web/worker WASM
compiler check. Browser checks remain deferred. Integrated assembly junction/oligo
transactions and unfinished design-draft persistence remain pending.
[Review contract](READ_REVIEW.md).

## Construct-owned sequencing checks

A construct can open Read validation with its reference pinned independently of
the currently viewed read. Exact reference/read events, options and workspace epoch
bind fresh reports and all actions; equal-DNA document/event switches cannot reuse
an old result. Discrepancies, uncovered ranges, insertion observations/calls and
read columns are paginated, with navigation to exact current/original chromatogram
call/sample coordinates. Reference insertion boundaries use a caret, not a guessed
base. Historical references/reads open as retained forks without rewinding heads.

Checks can be saved on the construct after fresh reconstruction and one durable
workspace transaction. Opening a saved record is historical; explicit replay
reconstructs the original events with the current engine, even after source heads
change. Navigating a read does not replace the pinned construct or invalidate an
unchanged saved check. Partial coverage remains incomplete with uncovered consensus
calls N. No sample identity, inferred correction or experimental success is asserted.

Worker42 adds capture/replay jobs105/106 and results97/98; only workspace-bearing
request/capture-result tags use128MiB, while ordinary analysis retains32MiB. MCP has
125 tools, including immutable native-workspace record/list/replay. Candidate11,
plan9, science-14 and archive1 remain unchanged. Nine focused core/protocol/real-SDK
cases pass, along with the full 125-tool stdio smoke and development Trunk build
of the main application, worker and packaged assets. Port8082 serves the refreshed
immutable read-only snapshot, retaining older versioned workers for open tabs.
Browser checks remain deferred. Explicit discrepancy-review decisions
and derived-reference revision proposals remain pending. [Saved check contract](READ_CHECKS.md).

## Selection-centered mutation intent

The sequence toolbar can send an exact DNA range or insertion caret to Mutagenesis.
The operator previews an explicit replacement/insertion/deletion, then opens the
intended reference before designing inverse PCR against it. The intended DNA is
built only from the original source and requested edit, independently of reaction
reconstruction. Circular origin-crossing edits disclose their origin shift; an
insertion at the displayed end caret preserves the original origin.

Opening retains the exact original event as a parent and saves source-bound intent
evidence in a native attachment. Repeat creation reuses only a matching target with
the same source event. Equal-DNA aliases keep separate ancestry; forged intended
DNA and stale events fail before publication. Source heads, manual primers and the
mutation form remain intact. The reference is visibly labelled intended and remains
restricted against physical-source promotion.

Worker41 adds source-edit job104/result96; MCP now exposes122 tools, including
atomic `build_mutation_intent` target/evidence retention. Candidate11, plan9,
science-14 and archive1 remain unchanged. Twelve focused new core/binary/MCP/history
staging cases pass; the existing native web unit target, full 122-tool stdio smoke
and final web/worker WASM compiler check also passed. Browser checks remain deferred. Shared unfinished-draft persistence and sequence-centered assembly
junction editing remain pending. [Intent model](MUTATION_INTENT.md) ·
[MCP flow](MCP_MUTATION_INTENT_FLOW.md).

## Golden Gate PCR preparation

The assembly dock now offers an ordered Golden Gate preparation/ligation method.
Existing pinned Type IIS digest fragments can be mixed with actual original-source
PCR→catalog digest rows. Site-bearing primer proposals show clamp, motif, spacer,
fusion and annealing domains separately; explicit Apply preserves other rows.
Fresh method evidence uses the same product/junction view and output finalizer.

The shared compiler emits PCR, exact catalog digest and optional end repair as
separate retained operations. The new source-bound planner handles BsaI, BsmBI,
BbsI and three-base SapI ends, including reverse-row bottom-strand geometry and
circular intervals. It refuses alternate binding or extra internal cuts under its
conservative profile. Right fusion DNA and source regions are explicit; the bound
independent target never supplies oligo sequence. Ordered reconstruction does not
establish ligation-mixture uniqueness or experimental yield.

Worker 40 adds planner job103/result95 and the new preparation tag5 in recipe
requests. Candidate11, plan9, science-14 and archive1 remain unchanged. MCP now has
121 tools. Twenty-one focused native preparation/planner/protocol/MCP and existing
recipe regression cases pass, as do the full 121-tool stdio smoke and final
web/worker WASM compiler check. Browser checks remain deferred. [Preparation](PCR_DIGEST_PREPARATION.md) ·
[Planner and strand convention](GOLDEN_GATE_PRIMER_DESIGN.md).

## Shared directional primer pair editing

Primers and Mutagenesis now share editable named forward/reverse oligos, separate
annealing domains and tails, explicit saved-primer reuse and explicit Apply.
Selection actions offer forward or reverse-complement DNA with exact event/range
anchors. The editor displays substitutions and compatible ambiguity in oligo
orientation; changing domain length leaves placement unresolved. Pair drafts survive
panel/document switches within the session; restoring a workspace resets them.

Mutation results are invalidated by draft changes and bound to exact source/target
events, including equal-DNA switches. Ordinary primer evidence retains its original
event for annotation/navigation. Four focused native model cases and the final
web/worker WASM compiler check pass. Browser checks remain deferred. Integrated junction editing and native persistence of
unfinished pair drafts remain pending. [Draft model](PRIMER_PAIR_DRAFT.md).

## Finalizing reviewed assemblies

The assembly dock now offers a named construct, optional primer collection and
optional preparation-document visibility. Finalization reconstructs the bound plan
against its independent intended reference and commits the complete replayable run;
hiding preparation documents never deletes their events. Existing edited outputs
are preserved. The chooser stages reference/candidate/output changes privately and
publishes them together only after a durable workspace transaction. Cancellation
or failure leaves the original workspace unchanged.

Worker 39 adds finalization job 102/result 94; MCP now exposes 120 tools and reserves
both result artifacts before publishing a changed workspace. Candidate 11, plan 9,
science-14 and native archive 1 remain unchanged. Fifteen focused native cases and
the full 120-tool stdio smoke passed. The full development Trunk build of the
web application, scientific worker and packaged assets also passed. Port 8082
serves the new immutable read-only snapshot, retaining older versioned workers
for open tabs. Browser checks remain deferred.
[Output model](CLONING_OUTPUTS.md).

## Direct cut choices and reviewed design outputs

Sequence/map restriction-site clicks can hand exact site identities to Virtual
digest while the assembly dock is open. Two-boundary controls select a directed
circular path or its complement, the full circle from a single cut, or explicit
linear source ends. Fresh selection must match one actual both-strand fragment;
internal cuts cannot be bypassed by slicing a desired reference interval.

The reconstructed-design view now exports primer CSV/TSV and HTML/PDF steps from
the reviewed plan and target. Order tables preserve per-use domain/tail and
mismatch details, use exact order-DNA/modification deduplication and spreadsheet
escaping, and expose unavailable order data. The web action blocks an incomplete
all-required sheet. Chemistry assumptions and downstream end treatments do not
become invented supplier modification requests.

Worker 38 adds three output/selection jobs; MCP `export_plan_oligos` brings the tool
count to 119. Thirteen focused native core/transport/MCP cases passed. The full 119-tool stdio smoke and final web/worker WASM compiler check also
passed; browser checks remain deferred.
[Cut selection](DIGEST_CUT_SELECTION.md) · [Oligo output](PLAN_OLIGOS.md).

## Digest fragments and actual product review

The digest fragment table can now place a freshly determined linear fragment into
the active assembly slot, retaining the original source event and exact catalog
selection. The shared preparation compiler emits an explicit catalog-pinned digest
and optional end repair. Sticky-ended reversal uses the original bottom-strand
source interval; it does not substitute a reverse-complemented top interval.

Fresh assembly evidence now includes a read-only actual-product map/sequence view,
coordinate and annotation navigation and opt-in fresh junction inspection. A
mismatched intended target can remain visible alongside actual reconstructed DNA;
a failed final step never substitutes a leading preparation as the final product.
Pointer/keyboard selection is available without permitting edits to the snapshot.

Worker 37, candidate 11, plan 9 and science-14 version the new operation/preparation.
Seventeen focused core/MCP/preparation/history/compatibility cases passed.
The full 118-tool stdio smoke and integrated web/worker WASM compiler check also
passed. Browser checks remain deferred.
[Catalog preparation](CATALOG_DIGEST_PREPARATION.md) · [Product review](ASSEMBLY_PRODUCT_REVIEW.md).

## Persistent assembly draft and independently declared intent

The first dedicated assembly dock now sits below the sequence/map
workspace. Its named backbone/insert slots retain source events, regions and
orientations across document changes. A separate pure intent builder concatenates
those declared regions into an expected reference before PCR/digest/assembly
simulation. Existing overlap-primer suggestions and preparation/verification
remain separate steps, using the same original sources.

Draft declarations and manual preparation settings persist in a bounded native
workspace attachment. Reopening checks the original source events and restores no
analysis badge. Source navigation can open the retained event as a separate
document, and candidate staging creates the intended-reference target with exact
source parents. No GUI/browser acceptance check has been run.

MCP now has `build_assembly_intent` (118 tools), producing an explicitly unverified
intended reference usable as the target of subsequent planning. Worker 36 adds
only its request/result tags. Five core intent fixtures, one binary transport case
and one actual MCP intent→primer→PCR/assembly case passed. Three native staging/history cases and the full 118-tool stdio smoke also passed.
The final integrated web/worker WASM compiler check passed. [Intent model](ASSEMBLY_INTENT.md).

## Sequence-centred cloning workflow and primer suggestions

The reviewed Benchling reference is method selection, backbone/insert slots filled
from sequence selections, preparation/primer review, a central junction/product
preview, then output creation. The current inspector is groundwork; a persistent
assembly draft and full source-view/product-preview connections remain required.
[Interaction design and references](CLONING_INTERACTION_DESIGN.md). Three parallel official-document/walkthrough reviews now map restriction/Golden Gate, primers/mutagenesis and verification/history to specific missing UI connections; their findings guide the existing cloning-first queue.

Deterministic PCR overlap-primer suggestions now connect to the existing
prepared-overlap compiler, web controls and MCP. Explicit source intervals and
orientations supply all sequence; a separately supplied target stays independent.
Suggestions require an explicit apply and fresh assembly check. Worker 35 adds
append-only planner request/result tags, without changing plan, candidate or
scientific-operation schemas. Twelve focused native core/transport/MCP cases and
the complete 117-tool stdio smoke and final web/worker WASM compiler check passed.
No browser checks were run.
[Scientific model](OVERLAP_PRIMER_DESIGN.md).

## Explicit starting-material termini

Original intact linear duplex references can now carry explicit operator terminal
hypotheses as a typed preparation step. Both ends initially remain unknown;
choices distinguish 5′ P/OH/unknown and 3′ OH/blocked/unknown. Ordinary WholeLinear
preparation is unchanged. This enables declared starting material to feed ordered
ligation without pretending that file metadata proves phosphorylation or 3′-OH.
Source DNA/features/raw citations remain unchanged, and the declared assumptions
are retained in plan history and HTML/PDF instructions. Imported physical or
reference-only intermediates cannot be cleared through this declaration.

Worker 34, candidate 10, workspace-plan 8 and science-13 explicitly version the
new operation/preparation variants. Historical layouts remain supported under
their original semantics, including readonly old reports; nested declarations
cannot hide in old worker/candidate/archive headers. Positional modifications
remain unknown rather than receiving an unmodified-DNA claim. The existing 116
MCP tools expose the new typed options without adding tools.

Twenty-five focused core/preparation/history/compatibility cases and two MCP
cases, the full 116-tool stdio smoke and the WASM compiler check passed. Browser
checks remain deferred. See [model and bounds](DECLARED_SOURCE_TERMINI.md) and
[declared starting-material example](../fixtures/declared-ends-demo.md).

## Junction DNA and reading-frame review

The reconstructed-plan panel now offers explicit junction inspection for ordered
ligation and overlap assembly. A fresh worker pass reconstructs original sources,
then reports actual final-product boundaries/homology, bounded DNA windows,
requested window-relative forward/reverse translations and available local CDS
coding evidence. Failed final reconstruction never substitutes a preparation as
the product; wrong independent targets retain their invalid assessment alongside
actual reconstructed DNA. Ordinary gene-length CDS evidence is supported with
bounded local output and explicit review/phase limitations.

Known sequence ends and short circular products clip the requested display flank
without making a complete inspection incomplete. Missing/unsupported scientific
evidence and omitted requested junctions remain explicit. Detailed findings and
source hashes are expandable, rather than filling the main cloning workflow.

Six focused core/worker cases and two MCP cases passed. Native web test targets
and the web/worker WASM compiler check passed. The full actual stdio smoke check
passed with 116 tools and the current candidate-schema expectation. Browser checks remain
deferred. Worker schema 33 adds one readonly analysis; MCP adds
`inspect_assembly_junctions` (116 tools), preserving existing candidate/plan and
scientific-engine layouts. See [review UI](ASSEMBLY_JUNCTION_REVIEW.md),
[analysis contract](ASSEMBLY_JUNCTIONS.md) and
[genetic-code example](../fixtures/junction-frame-demo.md).

The isolated native-build correction `99c4c02` cleared the compiler failure in
CI. Its WASM job passed; Linux native tests/build then reached stdio smoke, where
a stale candidate-schema-7 assertion rejected the current schema-9 writer.
The script now reads the current writer constant while retaining magic, source
binding and scientific assertions. The following pushed integration `be78a8a` passed the complete native Windows/macOS/Linux and WASM CI matrix ([run](https://github.com/henriksson-lab/collagene2/actions/runs/45507077269)). Subsequent terminal-declaration work `c301597` also passed all four CI jobs: native Windows/macOS/Linux and WASM ([run](https://github.com/henriksson-lab/collagene2/actions/runs/45509508425)).

## Ordered multi-fragment ligation

The shared ligation recipe compiler uses the existing source preparation compiler
for actual PCR/digest/repair steps, adds explicit terminal treatments, then emits
ordered oriented ligation under declared closure/nick policy. Fresh end/junction
and product-nick evidence comes from private material state in the same verifier
pass; imported metadata cannot supply that state. The UI and MCP share this
compiler, preserve original source roles, and retain independent-target checking
and existing candidate/history/HTML/PDF workflows. Returning to Assemble preserves
the selected workflow instead of resetting ligation to overlap.

Worker schema 32 adds one preview job/result; MCP adds
`propose_ligation_assembly` (115 tools). Existing plan/candidate/scientific schemas
are unchanged. Twenty-one focused core/worker/regression cases and two MCP
ligation/protocol cases and the corresponding WASM compiler check passed;
browser checks remain deferred. See [workflow](LIGATION_WORKFLOW.md),
[core contract](LIGATION_RECIPES.md) and
[three-piece example](../fixtures/ligation-demo.md).

## Saved primer reuse and PCR template history

The shared `core::primer_reuse` resolver validates explicit saved primer sequence
fields, rejects duplicate/conflicting qualifiers, preserves names and 5′ tails,
and uses exact current geometry only when required fields are absent. Stored
mismatch primers remain explicit sequence data with warnings. Actual imported
oligo annotations are supported; bare DNA does not invent a primer partition.

PCR and both overlap/TA–GC PCR preparation rows now offer an explicit saved-primer
picker across open documents. Selection is bound to the exact document/head event
and actual revision, with bounded searchable lists and full domain/tail review.
Copied fields remain editable and downstream binding/reconstruction is fresh.
MCP `resolve_primer` provides the same immutable resolution; the registry now has
114 tools. No worker/candidate/plan schema change is needed for literal primers.

Direct PCR previews now retain their actual template event as a parent, stage all
products before workspace validation/publication, preserve named ordered primers,
and reject asynchronous results after template identity, event, workspace epoch
or primer fields change. This closes the previous source-parent omission in this
UI path. Four focused core fixtures and one MCP resolution/PCR fixture passed.
The corresponding WASM web/worker compiler check passed; browser acceptance remains deferred. See [reuse contract](PRIMER_REUSE.md) and
[explicit-tail example](../fixtures/primer-reuse-demo.md).

## Prepared overlap cloning workflow

The shared `core::assembly_recipe` compiler now connects original source roles to
explicit whole-linear, exact-PCR or selected-digest preparations, optional typed
end repair, ordered orientation and profiled overlap assembly. The worker freshly
reconstructs preparations and checks the independently selected target. The UI
shows prepared fragments and actual junction words, and the retained candidate
uses original sources and preparation steps rather than cached products. PCR,
Assemble and Plan now sit directly beside the editor/digest inspector controls.

MCP `propose_overlap_assembly` uses the same compiler and verifier, atomically
retaining a candidate and complete evidence artifact for existing verification
and HTML/PDF exports. Worker schema 31 appends
only preview job/result tags, while candidate/plan/scientific versions stay
unchanged. Six focused core/worker fixtures and the MCP preparation/protocol fixture passed.
The following primer-reuse integration has separate native/WASM check results. Browser checks remain deferred. See [workflow](ASSEMBLY_WORKFLOW.md), [core contract](ASSEMBLY_RECIPES.md)
and [literal PCR/digest example](../fixtures/prepared-overlap-demo.md).

## Current editor implementation

Direct sequence typing/deletion, pointer and keyboard range selection, insertion carets and keyboard undo/redo are now implemented in the sequence canvas. A shared source-bound restriction scan drives the enzyme picker, site-count filters, sequence recognition/cut tracks and circular/linear map markers. The scan-only worker job appends schema 23 without reconstructing digest fragments. Browser verification is deliberately deferred until the in-scope code is complete; historical acceptance results below do not validate these new controls. [Controls and remaining boundaries](EDITOR_CONTROLS.md).

## Existing read-alignment integration — further work deferred

SAM/BAM collections now have an immutable source-backed core adapter, explicit actual-DNA reference checking, four MCP tools (101 total) and schema-24 worker sessions. The browser viewer and raw-source workspace retention/reopening are implemented. These sources are separate from editable molecules; saved reports grant no restored mapping authority. Pane sizes/view choices also have bounded, separately persisted display preferences. Background native core/MCP/asset and WASM web compiler checks passed at this checkpoint; new fixtures have been authored but not executed, and browser checks remain deferred. [Workflow and boundaries](READ_ALIGNMENT_WORKFLOW.md).

## Variant browser and MCP integration

VCF/BCF sources now have retained worker29 sessions, bounded record inspection, explicit DNA/ALT mapping and fresh reference-only previews in the web app, plus seven MCP tools (112 total). Opening a browser preview retains exact original variant bytes, exact source-reference event ancestry and a bounded unverified record/ALT recipe. Native/MCP replay reconstructs against those actual inputs rather than trusting saved results. Workspace reopening reparses files and restores no scientific authority. Standard gzip VCF and strict BGZF inputs use bounded decoding and integrity checks. Native type checks and the development WASM build passed; port8082 serves the updated immutable snapshot. Browser checks remain deferred. [Workflow and profiles](VARIANT_FORMATS.md).

## GPU sequence text and variant exchange

The sequence editor now uses a bounded WebGL2 glyph atlas/instancing overlay for visible reference/complement letters, synchronized with the existing2D interaction canvas and immediate2D fallback. Browser/GPU performance validation remains deferred. [Rendering boundaries](GPU_SEQUENCE_TEXT.md).

An immutable VCF/BCF exchange crate and explicit actual-DNA reference-allele checker are implemented as native/WASM APIs. Original bytes, typed field/genotype evidence, bounded conversion and allocation/framing guards remain separate from allele truth. Concrete edit previews require a fresh matching REF and explicit alternative selection. GUI/MCP import integration now uses private retained collections, explicit mappings and freshly computed reference-only previews. Native core/MCP/web fixture type checks and WASM compiler checks passed; new fixtures are authored but not executed locally. [Format and checking scope](VARIANT_FORMATS.md).

## Final physical TOPO acceptance

Reviewed web candidates now carry final `CheckedTopo` operations with all six explicit policy choices. The verifier reconstructs once against the independently bound target, retains graph/annotation evidence before rejection or withholding and can report consistency under the declared graph assumptions for fully specified accepted DNA. Native capture shares the same rotation/reverse comparison. Exported reference views retain their chemistry restriction; no nick repair or downstream intact-duplex authority is inferred. Background native core/MCP/web fixture type checks and the development WASM build passed. Port8082 serves the new immutable compiled snapshot. New scientific fixtures were authored but not executed locally; browser checks remain deferred.

## TOPO annotation integration

Annotation correspondence now follows actual selected top-call provenance through explicit reverse/circular preparation and PCR template intervals. Primer tails and synthesized additions supply no invented source coordinates; PCR-derived features require review. Repeated/unsupported mappings retain original source features and reasons. Additive annotated wrappers preserve the older physical DTOs and feed reconstruction evidence, MCP projection documents, web inspection and HTML/PDF summaries. Background native core/MCP/web fixture type checks and WASM web compiler checks passed. New scientific fixtures were authored but not executed locally; browser checks remain deferred. [Mapping boundaries](TOPO_ADAPTER.md#source-coordinate-annotations).

## TOPO web preparation

The web integration adds an explicit TA preparation form for declared activated vector strands and exact-PCR insert preparation, worker-backed orientation/projection/target inspection, six separate acceptance-policy allowances and source-event-bound candidate review. The form constructs typed declarations; it performs no scientific reconstruction or target-guided path selection on the main thread. A literal [synthetic TA example](../fixtures/topo-ta-demo.md) supplies independent source/target records and retains two nicks. Background native core/MCP/web fixture type checking and WASM web compiler checks passed. New fixtures were authored but not executed locally; no browser checks have run. [Workflow](TOPO_WEB_WORKFLOW.md).

## TOPO plan integration

The explicit source-backed TOPO adapter now has a typed plan operation with full physical step evidence, versioned candidate/workspace/worker transport and four MCP entry points. HTML/PDF export records actual preparation assumptions, PCR oligos, transfers and retained physical state. Browser candidate review exposes the physical evidence separately from reference DNA. A unique paired contour does not seal nicks or remove flaps; the initial plan integration remains Incomplete without explicit physical target acceptance and downstream graph composition. Source/event history remains replayable. Background native core/MCP type checking (including fixture sources) and WASM web compiler checks passed. New scientific fixtures were authored but not executed locally; browser checks remain deferred. [Adapter and remaining boundaries](TOPO_ADAPTER.md).

## Working scientific foundation

The independent `collagene-recombination` workspace crate now connects to source annotation projection, typed plans and native history, worker schema 22, MCP and a desktop panel. Explicitly choose two sources, the pinned manufacturer-context profile and a reciprocal component; the separately specified target never selects or repairs a product. Fresh reconstruction, exact source events, reviewed/unmapped annotation evidence and all reciprocal components remain available. Matching results stay **Incomplete**: duplex chemistry, biological repair, efficacy and broader Gateway coverage remain unfinished. Native and GenBank/EMBL reference-only guards prevent later operations from silently establishing chemistry. [Web workflow/example](RECOMBINATION_UI.md) · [Adapter and bounds](RECOMBINATION_ADAPTER.md) · [Scientific model](../crates/recombination/README.md) · [MCP](RECOMBINATION_MCP.md).

The published reference-recombination checkpoint `3ea00bc` passed **939 native tests and eight doctests**, including 25 engine tests, 13 adapter adversarial/annotation tests, 11 typed-plan/history/export tests, two worker boundary tests and the imported GenBank demo. Native SDK tests and the actual 96-tool stdio server cover source inspection → typed proposal → independent verification → source-only history replay. The new browser workflow passes reciprocal source oracles, explicit selection, wrong-target rejection, invalidation and HTML/PDF draft output; all 27 desktop suites pass against the final staged release, including the existing HTTP-cache worker-chain proof. WASM/schema checks and a full staged release build pass.

The shared Rust core validates DNA/RNA/protein alphabets, topology, content revisions and half-open compound annotation locations. Exact extraction rejects unresolved/fuzzy/order semantics. Edits remap supported coordinates, flag affected annotations and remove stale translations. Known unremappable GenBank reference/qualifier evidence blocks unsafe edits/exports.

Unified bounded search covers DNA/RNA/protein sequences, both nucleotide strands/circular origins, annotations/qualifiers, saved primer annotations and starter-catalog enzyme names/recognition sites. Stable hits retain compound coordinates, source revisions, actual cut boundaries and uncertainty; hit/row/work exhaustion is explicit. [Search profile](SEARCH.md). The starter catalog now has twelve enzymes, including a 3′-overhang example; it is not a supplier/methylation database. Reference-strand digest lengths remain distinct from the new two-strand fragment geometry. Exact PCR enumerates competing oligo orientations and single-oligo products, including introduced tails, within its documented binding model.

The duplex model tracks both actual 5′→3′ strands, their relative extents, blunt/5′/3′ ends, terminal phosphorylation/3′ chemistry and strand-specific nick boundaries. Complete digestion handles reverse Type IIS sites and circular-origin cuts; fragment indexes follow ascending top-reference cut order. Ordered ligation checks complementarity and junction bonds, supports explicit orientation and terminal phosphorylation/dephosphorylation, and reports retained nicks only when allowed. No nick repair or exhaustive mixture-product uniqueness is assumed. Uncut circles stay circular. Unsupported crossing cuts/single-stranded intermediates, unknown terminal states and unsupported physical templates remain explicit errors/incomplete results.

End repair models explicit recessed 3′ fill-in and protruding 3′ removal, with actual source/product strands, base ancestry, end chemistry and retained nicks. Affected 3′ ends must have accessible OH; no phosphorylation or nick sealing is invented. Exhaustive bounded antiparallel oligo placement reports all concrete compatible geometries, with uniqueness only for complete one-candidate searches. The verifier replays selected physical operations against the independent target; preview metadata has no physical authority. HTML/PDF protocols include actual ordered oligos and terminal declarations. [End-repair model](END_REPAIR.md).

Explicit single-nucleotide 3′ tailing now transforms established blunt duplexes with accessible hydroxyls, preserves 5′ chemistry/internal nicks and gives new residues no invented template ancestry. The Cloning methods panel reconstructs digest/optional blunting or unique exact PCR before TA/GC tailing, applies declared phosphorylation/dephosphorylation and checks ordered ligation against the independent target. Both added junction bases matter; nicked and sealed products remain distinct. Generic endpoint states are modeled, not factory-kit activation or polymerase completion. [Tailing model](TAILING.md).

Oriented Gibson/HiFi/In-Fusion profiles now check exact end homology and original-reference occurrence coordinates, with versioned advisory length guidance and bounded repeated-word searches. Those sequence reports do not establish reagent compatibility or unique mixture products. The verifier independently reconstructs the candidate and marks incomplete searches; known blocked 3′ physical substrates cannot silently enter either legacy or profiled overlap assembly. Strict operation/orientation/primer/sequence JSON parsing rejects unsupported chemistry fields. [Overlap profiles](OVERLAP_PROFILES.md).

Typed PCR/overlap/digest/end-treatment/ligation plans bind immutable inputs and an independently specified target. The verifier reconstructs products, compares final sequence/topology (circular comparison permits rotation), and distinguishes invalid, incomplete and consistent-under-assumptions results. Failed plans retain already reconstructed intermediates. Imported physical metadata cannot grant end chemistry or erase nicks. Known single/mixed-stranded GenBank declarations are retained and rejected by the intact-duplex model; PCR on known staggered/nicked intermediates is incomplete until modeled explicitly.

Cloning-plan products and standalone PCR reference previews project source annotations through actual PCR tail geometry; plans also map digestion, overlap and oriented ligation. Intact biological extraction order survives circular-origin joins; clipped/interrupted features are reviewed, and unprojectable or homology-collapsed features remain source evidence. Equal PCR sequence bytes with different template origins remain separate preview geometries and ambiguous verification inputs. Standalone previews bind their parent revision and oligos without asserting target verification or physical chemistry. Ancestral opposite-strand-only annotations and additional coordinate-bearing qualifiers remain pending. GenBank/EMBL preserve Collagene review-required flags and validated wrapped library hash/ID qualifiers; text exports explicitly report omitted quality, library and annotation-source evidence. [Annotation mapping rules](ANNOTATIONS.md).

Embedded translation supports NCBI codes 1, 2 and 11, ambiguity-aware codons, explicit initiation policies and compound CDS phase/code qualifiers. Bounded linear/circular six-frame ORF predictions retain exact/partial locations and possible-stop evidence. Independently implemented affine pairwise alignment supports global, local and query-in-reference modes with certain versus ambiguous-compatible matches, CIGAR, coverage and differences. Native bounded circular/reverse comparison exposes search-completeness flags; the initial web comparison uses the supplied origins/orientation. These are pairwise tools, not the deferred multiple-alignment engines. [Translation](TRANSLATION.md) · [Pairwise alignment](PAIRWISE_ALIGNMENT.md).

Basic primer tools calculate GC bounds and condition-bound perfect-DNA nearest-neighbor Tm, screen contiguous hairpin/dimer runs, and find both-strand/circular binding with explicit mismatches, ambiguity and certain 3′ runs. Anchored endpoint selection ranks bounded primer pairs with separate tails and recorded filters; candidates require independent binding/PCR review. Document-based primer annotations can be saved/reused and passed into PCR controls. No Primer3 or advanced folding engine is used. [Primer models and limits](PRIMERS.md).

Golden Gate supports four factual Type IIS profiles, complete physical cleavage, oriented cohesive-end assembly and junction-spanning recut checks. A bounded one-copy fragment-pool enumerator retains carryover, alternative paths/source histories, stable sequence groups and explicit truncation; it never claims exhaustive kinetic mixture search. The selected typed operation is reconstructed against a separate target, including trusted prior-step duplex states and source annotation mappings. [Golden Gate model](GOLDEN_GATE.md).

Intentional-mismatch PCR reconstructs actual oligo-encoded bases with recorded mismatch/3′ rules, competing roles and changed-feature review. Ambiguous/overlapping products remain explicit previews. Inverse-PCR substitution/insertion/deletion design creates a separately target-bound PCR → phosphorylation → sealed ligation candidate. Complete selected products that conclusively miss the target remain invalid even if other paths were unsearched; hard partial-step failures are not compared as completed products. [Mutagenesis models](MUTAGENESIS.md).

Reference gels derive bands from actual digest/PCR products under an explicitly supplied log-size ladder and resolution. Circular conformation, unsupported lanes and out-of-ladder fragments remain unpositioned; co-migration and optional mass assumptions stay visible. Bounded diagnostic single/double digests rank distinguished, indistinguishable and unknown candidate pairs. SVG/JSON exports retain calibration and source revisions. The default ladder is illustrative, not empirical; observed-gel comparison and OpenGel remain deferred. [Gel prediction](GEL_PREDICTION.md).

Quality-aware read validation uses explicit source-score thresholds and trimming, both orientations, affine query-in-reference alignment and independently counted optimal paths. Repeated/tied/unproven placements are excluded from support; duplicate original source files cannot inflate it. Circular origins retain source call/sample coordinates. Consensus leaves uncovered or insufficiently supported sites as N, with explicit gaps, indels and budget limits. FASTQ Phred encoding must be supplied; SCF/caller scores are not relabeled as calibrated Phred. [Read validation model](TRACE_VALIDATION.md).

Versioned user-owned feature libraries extract exact biological templates from selected source annotations. Bounded whole-template ungapped searches cover both strands/circular origins, allowed mismatches and IUPAC identity bounds. Palindrome/repeated placements remain distinct and limits explicit. Selected proposals are rederived against the same source/library/options before application; transferred annotations retain provenance and require review, with sequence-specific phase/translation/identifier claims kept as source evidence. Web library definitions persist in workspace attachments. GFF3/BED imports require explicit reference mapping, preserve source/rows/relationships/phases and apply only selected proposals. Original sources remain retrievable after edits/reopening; stale source claims cannot be exported as current coordinates. [Feature libraries](FEATURE_LIBRARY.md) · [Annotation formats](ANNOTATION_FORMATS.md).

## File formats and evidence

| Format | Current import/export | Evidence and limits |
| --- | --- | --- |
| GenBank/ApE | Nucleotide import/export; protein import, no GenPept writer | Standalone `gb-io` syntax adapter, compound locations/duplicate qualifiers, original record retained once; explicit strand type preserved; complete interoperability corpus still pending |
| GFF3 / GTF / BED annotations | Reference-bound import; supported-profile export | Explicit reference IDs, coordinate/strand/CDS-phase and Parent evidence, selected application, original sources and stale-source export guards; GTF supported-profile exchange is implemented |
| FASTA/plain DNA | FASTA import/export; plain DNA import | Explicit non-DNA alphabet tags; FASTA is a sequence-only export with loss reporting |
| EMBL/ENA nucleotide | Import/export | Modern/legacy IDs, DNA/RNA, shared INSDC feature grammar, original headers/evidence and unchanged record terminators; conservative guards for edited references/unknown coordinate evidence; no Swiss-Prot/submission validation |
| SnapGene `.dna` | Known uncompressed dsDNA import; original-byte download | Topology, compound/reverse/origin features, duplicate qualifiers and supported colors; original bytes/hash retained; unsupported encodings/flags rejected and unmapped semantics warned. No edited `.dna` writer or real vendor-file corpus validation yet |
| ABIF/AB1 and SCF | ABIF major 1 and SCF v2/v3 import; original-byte download | Canonically labeled signed signal channels, peaks, calls, raw caller/SCF accuracy evidence and exact original attachment; no invented Phred calibration or edited trace writer |
| FASTQ | Explicit Phred+33/64 import/export | Original records, scores and encoding retained; unsupported modifications never become plain DNA; generic call edits invalidate qualities; explicit trim/reversal retains mapped evidence |
| Oligo CSV/TSV | Explicit column mapping, selected DNA documents, original/normalized order exports | Full source table retained once; modified/unsupported oligos stay evidence; units and amounts remain recorded source values |
| Original Collagene XML / XDNA v0 | Supported-profile import; exact original download | Audited original serializer, sequences/aligned strands/notes/features/colors/primer evidence; binary header/feature coordinates and raw fields retained; unsupported physical states quarantined. No producer-generated biological corpus or edited legacy writer yet |
| Native `.cg2` | Experimental single-document save/load | Checksummed bounded binary payload; retains supported source attachments/metadata |
| Native `.cg2workspace` | Full workspace archive save/load | Strict bounded ZIP with hashed deduplicated blobs, document identities, revision DAG, exact undo/redo events, candidate snapshots and attachments; no imported verification authority |
| `.cg2design` / `.cg2plan` | Planning snapshot / immutable candidate exchange | Exact input/target revisions; the receiver recomputes status. Initial MCP snapshot exchange is limited to 8 MiB |

Workspace structure and limits are in [Workspace archive](WORKSPACE_ARCHIVE.md), table/quality rules in [Oligo and FASTQ formats](OLIGO_FORMATS.md), and read verification in [Trace validation](TRACE_VALIDATION.md). Format-specific coverage, sources and omissions are in [SnapGene](SNAPGENE_FORMAT.md), [EMBL](EMBL_FORMAT.md) and [trace](TRACE_FORMATS.md) profiles. Text conversions explicitly report unsupported source evidence. Annotated GenBank/EMBL exports include bounded retrieval provenance comments; FASTA still omits rich metadata/attachments. Original-source downloads are distinct from exporting a current edited view.

Codon design now supports reverse translation and synonymous optimization from canonical protein or complete concrete linear coding DNA, preserving the independently source-derived protein and a single terminal stop under declared ordinary-codon semantics. Explicit integer codon preferences, both-strand IUPAC motifs, fixed contexts, exact global/full-window GC and homopolymer constraints are checked by a separate full-scan algorithm. Suffix/backpointer search retains bounded work, frontier and prefix-arena data, distinguishes exhaustive unsatisfiability from incomplete search, and never claims global optimality after pruning. The web uses separate design/check worker jobs before Open/GenBank, exposes an independent supplied-DNA check, cancels source/options/event changes and records an exact source-event parent when opening a product. Bounded GenBank/EMBL comments preserve producer evidence without imported verification authority. Initial native debug checks produced checked 3 kb genes from 1,000-residue sources in 1.041 s and 0.460 s; these are local debug observations, not cross-platform performance acceptance. [Model and limits](CODON_DESIGN.md) · [Web workflow](CODON_DESIGN_UI.md).

Source-bound dot plots now compare fixed windows across DNA, RNA or Protein references, retaining both nucleotide orientations, original query coordinates and origin-crossing spans. Exact concrete-word indexing retains repeated positions and full collision checks; lower-threshold or ambiguous inputs use an explicit bounded window grid. Compatible and certain identities remain distinct, and sampling or any work/index/hit/output limit prevents an exhaustive absence claim. The web uses a retained WebGL2 anchor buffer with camera-relative integer/fraction coordinates for deep zoom, exact accessible hit pages and source-span navigation; a disclosed density fallback retains every scientific hit. JSON and standalone vector SVG preserve the model/options/source bindings. A sampled 1 Mb native debug case completed in 472 ms with 1,000 hits and 60,000 work units; its declared grid was complete, its whole domain was explicitly sampled. Reference-free read contigs now have a separate bounded embedded engine and web workflow. [Scientific model](DOTPLOT.md) · [Interactive view](DOTPLOT_UI.md).

## Current delivery priority

Prioritize the desktop viewer, GenBank import and virtual restriction digest as a coherent MVP. Provider expansion and additional advanced interface work are paused as delivery priorities while the existing core, worker/render boundaries and everyday UI are audited and completed. Independent post-MVP libraries may progress in parallel without changing the MVP source or its acceptance gate. Existing implemented capabilities stay in the full product scope and are preserved as extension modules; the immediate interface should concentrate on the MVP workflow. This is a delivery-order change, not a reduced completion definition. Mobile remains deferred until desktop completion and user review/edits.

## Web application and interfaces

Desktop browser workflows are the current delivery and UI acceptance target. Mobile design, control selection and acceptance testing are deferred until the desktop product is finished and the user has reviewed and edited it. Earlier exploratory mobile-width checks are historical observations; they do not establish a mobile product scope or gate desktop delivery.

The Leptos workspace opens documents directly, without organization/project setup. It provides selection/editing/undo, explicit clipboard copy/reverse/cut/replace/insert with asynchronous source/selection guards, contextual digest/PCR/plan controls, direct format exports and IndexedDB recovery. Annotation controls create/edit/delete compound locations and duplicate/flag qualifiers, with paginated access beyond 500 annotations; compound selections highlight every span and show exact biological sequence where defined. Reference controls rename, change topology, rotate origin and reverse-complement with coordinate/source-evidence checks. Biology and Compare panels run revision-bound worker jobs with protein/report downloads. The initial workspace allows 1,000 documents and an estimated 64 MiB of current document data, including retained metadata. Undo navigation retains up to 100 states while older revision branches remain retained. The history controller bounds new retained content to 64 MiB and 10,000 events per document, rejecting new edits at capacity rather than silently dropping evidence. Native workspace autosave preserves exact undo/redo events, branches, candidate snapshots and original attachments. Workspace replacement validates first and atomically writes the new archive plus a previous-workspace recovery copy; stale autosave generations cannot overwrite newer work.

WebGL2 supplies circular/linear map geometry, rotation and zoom; sequence text is virtualized Canvas2D. The viewport borrows the document without per-scroll cloning. A retained worker validates/indexes one immutable source and serves small binary viewport requests for configurable feature/primer, enzyme recognition/cut, six-frame translation, reviewed CDS and provisional ORF tracks. Reports retain exact codon mappings and disclose unsupported evidence/work/output limits. The web canvas and accessible paginated viewport list synchronize selections; frame/direction/code/enzyme controls and JSON report download are available. [Track profile](VIEW_TRACKS.md). Interactive circular-map and SVG terminal arrows share a biological-order rule; unordered, unresolved, fuzzy and protein locations cannot acquire a directional arrow. Map geometry uses an explicit 12 MiB budget with a simplification notice, rather than an arbitrary first-feature cutoff. Before the configurable worker tracks were introduced, a local 1 Mb/10,000-feature full-reference scroll measurement improved p95 frame intervals from 33.332 to 16.666 ms with all 10,000 map features included; this was short headless Chrome/SwiftShader testing, not cross-platform acceptance. [Measurements and open gates](PERFORMANCE.md). Trace documents have a virtualized Canvas2D chromatogram with signed amplitudes, overview envelopes, pan/zoom, shared call selection, keyboard navigation and bounded accessible evidence. Editing calls makes the immutable trace view stale rather than silently realigning it; original evidence still downloads. Native attachments and candidate proposals survive recovery; verification is recomputed after reload.

The Plan panel can retain a complete reconstructed graph through **Keep reconstructed products**. Exact source events are selected explicitly when identical content is ambiguous; operation inputs become product parents and the independent target stays separate and active. Capture re-verifies before an atomic current/previous workspace installation. History lists saved reconstructions as unverified historical evidence, with explicit current-engine replay, product/status/engine comparison and JSON download. Source heads and undo/redo remain unchanged. This covers retained full plans; standalone preview exports still need equivalent source-event retention. [Run lineage](RUN_LINEAGE.md).

Rust/WASM workers use Postcard and transferable binary buffers, bounded job limits and obsolete-result rejection. Current worker100/candidate26/plan24/science27 encodings are updated together; unsupported internal versions are rejected without fallback readers or migrations. Exact workspace-bearing jobs permit bounded128MiB payloads, while ordinary scientific messages retain32MiB limits. Recorded report engines remain documentary provenance; scientific reuse requires fresh current-engine replay. Sequence tracks debounce viewport jobs, reject obsolete replies and terminate retained workers on source replacement/unmount. Large-reference work uses private retained sessions and bounded owner-block advances with cancellation. A general worker pool and chunked plasmid-document caches remain pending. [Binary protocol](WORKER_PROTOCOL.md).

The session-only stdio MCP server uses the official Rust SDK, typed tool schemas, immutable resources and bounded outputs. Browser snapshots can pass through the actual MCP server and return as candidates for UI review, independent current-revision binding and WASM verification. Its190 tools expose restriction/PCR/duplex analysis, URL import and Addgene sequence choices, translation/CDS/ORFs, primer analysis/design, PDF export, Golden Gate/mismatch-PCR simulation, inverse-PCR candidates, gel/diagnostic artifacts, pairwise/circular comparison, read validation, explicit-encoding FASTQ, oligo collections/selections and full native workspace exchange, versioned feature libraries and GFF3/BED annotation datasets, plus end-repair/oligo-annealing previews, source-bound unified search/viewport tracks, immutable alignment collections/exchange/column inspection/row extraction, protein-property models, explicit tailing previews, oriented overlap-profile checks and immutable retained-run capture/list/replay plus source-bound codon design, independent candidate checking and paginated dot-plot analysis with JSON/SVG artifacts and source-bound read contigs with immutable report/consensus resources, plus paginated methylation assessments, CRISPR guide/off-target reports and independently checked edit-template resources, plus managed streamed reference import/export, global interval reads/guides and private cooperative genome searches. Native workspace upload uses bounded offset/hash-checked chunks and atomic installation, including a tested archive above 8 MiB and a 1,000-document archive with retained history, candidates and evidence. It has no implicit browser/project-write authority. Live synchronization, persistent authorized writes, hosted HTTP and cooperative cancellation for older one-shot scientific analyses remain pending. See [MCP usage](MCP.md).

The History inspector now renders recorded event DAG branches/merges and exact parent links with bounded SVG/DOM pages, keyboard inspection, navigation-stack badges and external-parent details. Equal-content events remain distinct. Opening an old event as a new document retains its exact source-event parent through native archive/reload. Graph content/candidate matches carry no scientific verification authority. [History graph](HISTORY_GRAPH.md).

The Biology inspector provides explicit protein mass, theoretical charge/pI and A280 models, declared cysteine/stop/modification assumptions, ambiguity intervals and unavailable reasons. Changing inputs cancels work and clears stale reports; no folding/activity/empirical spectrum is inferred. Independent numeric tests include primary residue data, high-precision alternative root solving, ambiguous-resolution enumeration and disulfide/modification checks. [Protein properties](PROTEIN_PROPERTIES.md).

Reference-free read assembly uses the independent `contig-1` affine overlap/layout/consensus model. It retains native raw-score semantics, original calls/samples, explicit duplicate groups, gap eligibility, branch/repeat alternatives and incomplete resource states; unresolved evidence yields `N` rather than reference-filled bases. The bounded Canvas viewer exposes exact evidence and source navigation. Native records bind original source events; restoring them grants no scientific authority, and consensus use requires fresh worker replay. Atomic evidence/product transactions retain every contributor as a history parent. [Scientific model](CONTIGS.md) · [Web workflow](CONTIGS_UI.md) · [Native evidence](CONTIG_HISTORY.md).

F21 now has a source-bound methylation assessment engine and web planning review. Its versioned supplier profile covers all twelve starter identities, records Dam/Dcm/CpG assumptions separately, inspects both-strand/circular motif context and preserves exact potential modified-base coordinates. Unknown or unsupported hemi/mixed states, possible-only contexts, unresolved some-overlap rules and conflicting formulation evidence stay explicit. Generic BsmBI is not silently replaced by BsmBI-v2, and the conflicting generic BsaI Dcm descriptions remain unsupported. Reports do not measure modifications or change existing assumed-cleavage products. A separate tracked-strand foundation and checked plan operations now preserve positional state through physical reconstruction. Supplier conditions, broader catalogs, more positional susceptibility rules and modification chemistries remain required F21 work. [Scientific profile](METHYLATION.md) · [Web review](METHYLATION_UI.md).

Current plan reconstruction carries a sparse `TrackedDuplex` alongside each established physical product. Original-input `DeclareSourceModifications` binds closed Dam/Dcm/CpG assumptions; `CheckedDigest` pins the exact catalog selection and reports the current positional profile, retains blocked sites and withholds the product set for unknown/impaired/incomplete sites. Reverse, assumed-complete/checked digestion, ligation including bottom-origin realignment, end treatment, repair, tailing and Golden Gate strand mappings preserve inherited calls. Ideal PCR and supported standard unmodified oligo annealing use explicit unmodified-base models; sequence-only overlap still provides no established physical state. New motifs are not methylated automatically. Assumed-cleavage operations on explicitly tracked modified/uncertain material remain incomplete. Bounded producer summaries never restore chemistry from imported metadata; native source/operation graphs and fresh replay remain necessary. [Positional model](MODIFIED_DUPLEX.md) · [Checked planning UI](CHECKED_RESTRICTION_UI.md).

F25 now has a declared-profile guide/PAM engine and a bounded pigeonhole-seed off-target index with exhaustive ambiguous-reference fallback. SpCas9 NGG/20-nt, declared Cas12a TTTV/23-nt and custom geometry retain both strands, circular-origin coordinates, nominal cuts, exact rational GC filters and explicit model limits. Selected references retain named roles, revisions, build/accession labels and min/max IUPAC mismatches; only an exact intended tuple is tagged, and equal DNA elsewhere stays visible. Coverage is selected-record scope, not a verified whole build or biological specificity. Literal homology donors derive from the source and replacement while a separate checker reconstructs the edit and compares an independent target. Fresh worker design/check, source-only design parents, atomic native installation and unverified producer provenance prevent restored claims from becoming current checks. [Guide/search profile](CRISPR.md) · [Donor model](EDIT_TEMPLATE.md) · [Web workflow](CRISPR_UI.md).

F25/F28 now has an initial immutable large-reference engine and native managed-file adapter. Streaming FASTA ingestion emits packed IUPAC blocks; `.cgref` framing separates exact layout/provenance identity from canonical DNA collection identity. Global `u64` interval reads/designs load bounded verified blocks, while cooperative whole-manifest search validates every owner and canonical contig hash using a private cursor. Exhaustive imported-contig coverage, source integrity and completeness of a declared expected build stay separate. Limits, cancellation, tampering and missing components prevent complete coverage claims. Native independent checks include a 33 MiB contig/full search and a 70 MiB FASTA to a native store exceeding 32 MiB, followed by streamed reimport/export; ordinary plasmid limits stay unchanged. The release WASM build and desktop browser suite pass with bounded Blob slices, separate IndexedDB blocks, semantic readback before publication, scoped rollback, global-coordinate navigation and exact-hash reattachment. Four reference links occupy a 13,000-byte workspace archive with only 10,932 uncompressed entry bytes; bulk DNA stays separate. [Desktop workflow](GENOME_UI.md). [Reference format](GENOME_FORMAT.md) · [Genomic search model](GENOME_CRISPR.md).

## URL import and exports

F35 now explicitly requires Addgene/NCBI URL import in the plan. The new `collagene-sources` crate keeps networking outside the scientific core. NCBI record URLs/accession.version resolve through documented EFetch, and direct public HTTPS sequence files use the same parsers as local import. Native retrieval uses embedded Rust HTTP/TLS, bounded streaming, rate limits and validated redirects; the browser uses cancellable bounded fetch where CORS permits. Canonical source, actual accession/version, retrieval timestamp and evidence hash are retained without credentials/signed query parameters. The MCP `import_url` tool uses the native adapter; `export_document` transfers the exact fetched native document into the browser when CORS prevents direct retrieval.

The authenticated Addgene catalog adapter uses its published OpenAPI endpoint and token scope. The browser dialog discovers deposited/Addgene and full/partial choices, requires an explicit sequence ID, and keeps tokens only in dialog memory. Native/MCP retrieval handles authorized download redirects; the browser needs direct CORS-permitted responses or an explicit native/MCP file transfer. Catalog/download bases and lengths are cross-checked when evidence exists. Selected identity/categories, response hash and original GenBank evidence survive save/export; signed links and tokens do not become provenance. **Live authenticated Addgene interoperability remains unverified**; the public schema and synthetic provider fixtures were used. Caching and live browser-native bridges remain pending. [URL adapter details](URL_IMPORT.md).

Standalone escaped HTML protocols re-verify the exact bound plan and include inputs, intermediates, chemistry/nicks, assumptions, missing laboratory conditions, provenance and an embedded vector target map. Embedded deterministic PDF export now does the same through a background worker or MCP artifact, with A4 pagination, repeated draft/status/revision headers, an embedded Unicode font and vector maps from the SVG scene. Unsupported text shaping is visibly escaped with warnings; this is not PDF/A or tagged PDF. Browser printing remains available. SVG reference maps also download directly. The interactive WebGL renderer still has its own layout implementation, so a shared interactive/export scene is not claimed. [PDF profile](PDF_EXPORT.md) · [Map export profile](MAP_EXPORT.md).

## Shared model and source snapshots

`collagene-model` now owns the existing sequence/location/document types and validation; `collagene-core::model` and its top-level exports retain their public paths. The extracted implementation is unchanged, and literal pre-extraction JSON, Postcard and revision goldens preserve data compatibility. New scientific engine crates can depend on the shared model without introducing a core dependency cycle. [Model contract](../crates/model/README.md).

Addgene catalogs now keep their source identity, response hash and choices private with read-only getters; the serializable summary retains its four-field JSON layout. A detached choice or JSON summary cannot change the original catalog or its download/import provenance. Five external compile-fail examples enforce that public mutation is unavailable. A live unauthenticated provider audit established sequence login requirements and rejected localhost authorization preflight; authenticated native interoperability still requires approved API access. [Live access evidence](URL_IMPORT.md).

## Independent post-MVP libraries

The reusable foundation checkpoint `7c239c0` passes **1,040 native tests and 14 doctests**, plus the SAM/BAM interoperability configuration and native/WASM library checks. The folder-policy package and MCP checks pass 80 focused cases; the actual 96-tool stdio executable passes with default and configured storage, leaving the configured parent empty after normal disposal. Local-storage cross-compilation checks pass for Windows and macOS. The current staged GUI passes all 30 desktop suites, including compact annotation rows, explicit overflow, stable source anchors, fractional canvas viewports, laptop dimensions, document pagination, exact base selection and retained source/history. Overall design polish remains unfinished. The foundation checkpoint passed all three native OS jobs and WASM ([CI run](https://github.com/henriksson-lab/collagene2/actions/runs/45450438682)). The subsequent compact-row and fractional-viewport checkpoint `4b031da` also passes all four jobs ([CI run](https://github.com/henriksson-lab/collagene2/actions/runs/45457750569)).

`collagene-map-label-layout` supplies dependency-free measured-box placement with explicit overflow and camera-independent budget selection. Nine native tests and a doctest cover boundaries, literal collisions, invalid inputs and 10,000-item cases. The web adapter caches source names and glyph metrics, projects at most 128 candidates per camera frame, retains original feature identities and exposes omitted/abbreviated counts. The final browser check passes actual circular/linear GenBank maps, identical-span feature selection, source-switch gesture cancellation, whole-cache name limits, zero rotation remeasurement, the 1 Mb/10,000-annotation case and exact native history/reload preservation. Split/Sequence/Map controls keep panels mounted. The integration passed all 30 existing desktop regression suites; final source-transition/CSS corrections are additionally checked with the map-label, navigation, viewport-resize and desktop MVP workflows. Native workspace checks now pass **1,054 tests and 15 doctests**. [Interactive label architecture and remaining work](MAP_LABELS.md).

`collagene-topo` models declared activated substrates and actual elementary backbone transfers, preserving every source strand, nick, flap/free displaced strand and enzyme release. Independent review distinguished the standard attached-flap diagram from complete four-base displacement; both retain their actual source calls. Twenty-nine mechanism/projection fixtures, nine source-adapter fixtures and a doctest pass, with native/WASM checks. The source adapter reconstructs declared material recipes or fresh exact PCR with explicit tailing, reserves aggregate primer-search work, retains the full physical graph and compares a separately supplied target under an explicit nick/end policy. A proven reference contour is a tainted coordinate view; it does not replace the physical graph. An independent exhaustive oracle checked 9,624 small backbone graphs for contour proof consistency. Typed-plan/history installation, worker/MCP/UI integration, kit preparation profiles, annotation projection and broader chemistry remain unfinished. [Source adapter and boundaries](TOPO_ADAPTER.md). [API, sources and scope](../crates/topo/README.md).

`collagene-view-layout` provides dependency-free deterministic annotation lane packing with bounded mergesorts, ordered disjoint interval indexes, whole-item joined placement and explicit overflow. The original feature identities and all intervals remain accounted for. Twelve native tests, a doctest, the exhaustive small-domain oracle and two 10,000-item scale cases pass; the pure viewer adapter adds seven source-preservation/resource fixtures. The viewer integration passes the real 1 Mb/10,000-annotation browser case, delayed viewport response and middle/end reflow tests, exact annotation/base selection and native source/history comparisons; it is included in the 30-suite desktop gate. [Layout API and bounds](../crates/view-layout/README.md).

`collagene-read-alignments` provides bounded streaming SAM/BAM/BGZF parsing/writing, sealed source evidence, CIGAR runs, reverse/hard-clip coordinates and loss-explicit conversion. Source integrity/EOF receipts remain separate from reported alignment assertions. Thirty-five native fixtures, four doctests, an optional independent noodles interoperability case and a WASM release build pass. Host/UI integration, reference-bound validation, coverage/variant evidence, indexed queries and CRAM/VCF/BCF remain future work. [API and format contract](../crates/read-alignments/README.md).

`collagene-local-storage` confines managed native data to directory capabilities. MCP accepts an existing `--storage-root`, creates private session storage and never accepts a tool-controlled host path. Replacement-aware cleanup preserves nonowned files. A public preview helper serves only an immutable read-only compiled-asset snapshot, with no filesystem API. The actual public-address viewer/import/digest workflow and traversal/write/symlink rejection checks pass. [Folder setup and precise boundary](LOCAL_STORAGE.md).

`collagene-protocol-units` provides checked exact nonnegative rational quantities, dimension-specific units, dilution and named mixture balances, explicit quantization with achieved concentrations, and DNA molecule ratios under a caller-declared mass model. It has no default dependencies or implicit molecular-weight coefficient. Optional serialization preserves full u64 precision as decimal strings and rejects dimensional mismatches. Numerical calculations do not establish reaction compatibility or experimental success. [API and assumptions](../crates/protocol-units/README.md).

`collagene-enzyme-catalog` supplies a sealed versioned edition with forty enzyme identities, four distinct supplier formulations and thirty-two documentary condition observations. Exact recognition/cut descriptors, reverse-complement equivalence, distinct isoschizomer/neoschizomer identities and explicit ambiguous lookup are separate from sequence matching. Context-specific resolution preserves unknown/conflict states; curated-fact checksums do not authenticate vendor sources. The SmaI chart-temperature observation retains its supplier-documented extended-digestion exception. The desktop digest now binds the forty-identity catalog through `catalog_digest`, reusing the established inspection, physical-product and exact-parent history logic. Catalog selection is retained as unverified producer evidence in native and GenBank/EMBL products. Reaction-condition interpretation remains separate work. [Catalog API, source citations and limits](../crates/enzyme-catalog/README.md).

The integrated libraries pass forty-one independent native tests plus three doctests, the zero-default-dependency units configuration and native/WASM compilation. They do not change the previously accepted desktop MVP renderer or digest pipeline. MCP exposes `calculate_dilution`, `balance_mixture` and `calculate_dna_molar_ratios` using these same algorithms, with decimal-string inputs, closed unit choices, explicit model assumptions and reported achieved concentrations after requested rounding. All sixty-six MCP tests and a freshly built ninety-three-tool executable stdio smoke pass. The server also exposes versioned catalog discovery/evidence, preserves strict legacy inspection artifacts alongside complete catalog-pinned artifacts, and rejects stale pins.

## Annotation and index format additions

GTF exchange now uses the existing reference-bound annotation pipeline in the core, browser and MCP. It retains row order, repeated quoted attributes, phase and strand without guessing transcript assemblies. Original input remains separately downloadable and application requires selected proposals against the same reference revision. Worker schema 30 adds only the nested GTF selector, retaining historical GFF3/BED binary tags. See [GTF profile](GTF.md).

The dependency-free `collagene-genomic-index` crate supplies bounded BAI and decoded-CSI parsing with strict BGZF transport in the core adapter and conservative merged BGZF candidate chunks. `core::indexed_bam` binds actual private BAM source bytes and header order, filters reported mapped spans, and independently scans all retained records to detect omitted matches. It rejects missing-index matches rather than reporting false absence. Fresh bounded BAI construction uses actual record offsets/spans and checks real coordinate sorting. This initial bounded adapter does not yet provide streamed indexed browsing, CRAM or variant-index hosts. Index metadata and candidate chunks are not scientific verification. Focused GTF and compatibility fixtures (six), BAM selection/BGZF CSI fixtures (three), index parser/query fixtures (five) and variant exchange/compression fixtures (nine) passed locally. Fresh BAI construction also passed its four focused selection/transport/writer cases. The final WASM and MCP GTF checks stopped on disk exhaustion; neither is recorded as passed. Redundant local build artifacts are being removed. No browser checks ran. Further BAM/BAI work remains deferred. See [index profile](GENOMIC_INDEXES.md).

## Validation and remaining delivery work

The shared-model/Addgene checkpoint `43e1913` passed **882 tests: seven shared-model, 721 core, sixty-six MCP unit/integration, twenty-five source-adapter, eleven native web-view, twenty-two exact-unit, nineteen catalog and eleven build-time worker-packager tests**, plus eight doctests. These include independently hand-specified cleavage/junction cases, exhaustive small alignment/ORF oracles, all-origin/reverse checks, ambiguity and coordinate failures, binary corruption/migration, parser fixtures, evidence preservation and independent physical/annotation/primer-model audits. Primer audit oracles independently check mass-action thermodynamics, 46,656 ambiguous dimer pairs, 6,561 hairpins, 18,750 oriented binding cases and ranked circular designs. New independent review checks cover six-frame circular codon mappings, ambiguity certainty, palindromic recognition, bounded labels and empty qualifier accounting. Alignment exchange has exhaustive gap/column mapping and annotated round-trip oracles. Worker jobs have independent legacy-layout fixtures; scientific MCP tools and binary PDF artifacts are exercised over real SDK JSON-RPC with revision bindings, pagination and uncertainty/limit checks. Independent Poppler tools also parsed, extracted Unicode text and rendered generated linear/circular PDFs. Codon audit oracles independently enumerate every canonical amino-acid pair and all 216 KLS-stop codon combinations, with literal preferences, exact GC inequalities, both-strand/context motifs and full-window checks; two default-budget kiloresidue sources yield independently checked 3 kb designs. Dot-plot oracles enumerate short ambiguous pairs in all source topologies, both orientations, every tested threshold and exact original spans; all protein symbols, a complete indexed 20 kb comparison and sampled 1 Mb comparison have independent checks. Native camera cases compare f32 projection arithmetic with f64 at one-base zoom near 10 Mb. Contig audits independently enumerate affine boundary paths for 900 small read pairs and recount original source calls, eligible gaps, duplicate groups, strands and consensus; native history tests preserve all contributing event parents and historical unverified evidence. Methylation audits enumerate all 256 concrete XbaI flank combinations across linear/reverse/circular origins, independently check both-strand Dam/Dcm/CpG target coordinates and asymmetric Type IIS cut boundaries, and exercise unknown contexts, limits, formulation gaps and fresh replay. Tracked-state audits independently recount original strand chemical positions across all origins/orientations, verify staggered fragment and circular bottom-coordinate mappings, retain modified bases after generation motifs break, distinguish copied versus newly synthesized repair bases, and reject forged summaries, stale profiles and descendant re-generation. CRISPR audits independently enumerate full Hamming/IUPAC windows for both PAM sides, orientations, circular origins and every supported mismatch threshold; they check rational GC boundaries, exact intended roles and every resource cap. Donor audits independently check literal replacements, source annotations, target independence, every supplied-base mutation, and all edits/orientations of short binary sources. No other cloning application is a correctness oracle. The release web/worker build, integrated WASM checks and formatting check pass. One live native Rust read of `NC_001422.1` returned 5,386 bases with matching accession/version; the offline suites do not require provider availability.

Twenty-six desktop Chrome workflow suites pass against the final release: the focused viewer/GenBank/virtual-digest MVP, catalog digest and worker-release cache acceptance, plus twenty-three existing regression suites. The cache acceptance preserves the existing HTTP cache, arms unusable legacy worker URLs, verifies actual loader/JavaScript/WASM requests share the generated content hash, and completes isolated analysis and retained viewport jobs before and after reload. The focused suite independently checks physical multi-enzyme products, fresh GenBank headers/annotation projection, exact source-event history, stale/cancelled work, uncertainty and resource bounds, two-slot quota/cancellation rollback, retained unverified provenance and 1600/1280/1200-pixel desktop layouts. The final regression run completes all twenty-three suites without a fresh runtime error. Mobile acceptance is deferred; earlier exploratory mobile-width observations below are historical, not current delivery gates. They cover real map painting, selection/edit/undo/recovery, annotation creation/compound extraction/duplicate qualifiers/deletion, reference topology/origin/reversal, worker translation/CDS/ORFs and review gates, pairwise coordinates/coverage/JSON export and option invalidation, worker PCR/digests, correct/incorrect overlap and physical ligation plans, explicit nicks, actual stdio binary exchange, HTML/SVG/browser PDF, SnapGene/EMBL/SCF file opening, chromatogram selection/stale-call recovery/original evidence, mocked NCBI fetch/provenance, Addgene explicit catalog choices/auth routing/in-flight selection cancellation/partial-record evidence/secret-free exports, primer save/reuse/design and tailed-PCR annotation retention, embedded PDF worker download with independent text extraction, Golden Gate candidate review/reverification/caps, mismatch changes and reviewed gene annotations, inverse-PCR physical closure, gel co-migration/mass/unknown circles/SVG, diagnostic ranking, and mobile width. The new suites also exercise native workspace archive download/reopen, exact undo/redo events and retained branches, candidate active/closed state, correct-CRC/stale-hash tamper rejection with both IndexedDB slots unchanged, previous-workspace recovery after corrupt autosave, explicit FASTQ encoding/original export, threshold-invalidated support-only consensus, selected oligo documents and exact retained table recovery with modified rows and normalized quantities. Browser helpers wait for storage restoration and completed autosave before starting imports or reloading. Annotation/library browser checks cover explicit GFF3 reference mapping/coordinate extraction, BED12 reverse compounds, original raw input downloads after reopening, both-orientation library proposals and identity bounds, selected reviewed application, undo/redo, persisted library definitions, duplicate prevention, worker exports, invalid-UTF-8 rejection and mobile width. Physical-operation checks additionally verify actual 5′ fill-in/3′ removal products, candidate v5, chemistry-sensitive phosphorylation/closure, independent wrong-target rejection, complete/incomplete annealing evidence and actual oligo strings in independently extracted PDF text. Clipboard/migration checks exercise real clipboard permissions and keyboard actions, cut/insert/replace/undo, stale-selection and identical-content/different-document cancellation, failed-write Cut rollback, exact legacy XML/XDNA original bytes, aligned-source displays and mixed-strand guards. Search/track checks additionally exercise literal EcoRI boundaries, CDS/ORF/primer evidence, real canvas codon selection, sub-512-byte cached viewport requests after source installation, option/scroll invalidation, both paginated result lists, compound-origin navigation, protein search and mobile controls. The preceding cached-track checkpoint retained 16.666 ms p95 rotation/scroll intervals for the 20 kb and 1 Mb fixtures on the same local software renderer; it does not establish per-frame scientific track completeness. Alignment UI checks additionally cover explicit supplied alphabets, literal gap/residue coordinates and imported annotations, original/canonical exports, all-gap atomic rejection, exact source/native archive recovery, cancelled stale imports, row pagination and navigation to column 4,000,000 with bounded painted cells. History checks exercise real equal-content events, branch edges, navigation stacks, source-bound forks/reload and more than 100 retained events. Protein checks run literal mass/extinction/disulfide cases in WASM, ambiguity/unsupported-state bounds, JSON artifacts and cancellation. Method UI checks exercise literal T/A and C/G junction products, free-OH versus phosphorylated PCR products, sealed/nicked ligation, wrong independent targets, actual HTML/PDF oligos, candidate 6/plan 4, reverse-source homology coordinates, advisory length failures, bounded incomplete scans, source/option cancellation and pointer access to the final document in a scrollable list. Retained-run checks additionally verify explicit equal-content source choices, exact six-step TA product edges, independent target exclusion, unchanged source heads/navigation, unmount cancellation, synchronous transaction rollback, recovery copies, fresh replay/JSON and historical authority reset after archive/reload, including mobile width. Codon UI checks additionally exercise literal weighted codons, independently matching required/candidate proteins, both-strand/context motifs, GC/homopolymer/full-window constraints, exhaustive failures versus incomplete search, ordinary code-11 and code-2 assignments, separate supplied-DNA checking, equal-content source cancellation, exact derived parent/recovery evidence, GenBank comments and a default-budget 1,000-residue to 3,000-base WASM design, including mobile width. Dot-plot UI checks additionally compare literal full-window/reverse/circular/ambiguous/protein/RNA cases, original spans and previews, both origin-crossing reference/query navigation paths, sampling/domain flags and all declared resource caps, exact-source cancellation, disclosed Canvas aggregation, JSON/SVG artifacts, alphabet-correct labels, visible WebGL pixels and a complete 30 kb self-plot with 29,986 retained windows. That local SwiftShader run measured 16.666 ms p95 zoom frames and no vertex-buffer uploads during camera movement; it is not a general hardware/browser performance gate. A later final-release shared-host run measured 49.998 ms for the same local 30 kb dot-plot zoom case; this checkpoint does not establish a fixed 60 fps guarantee. Read-assembly browser checks additionally reconstruct a literal three-read gapped/reverse/quality fixture with an independent source-call oracle, exercise trimming/support/budget policies and repeat/branch diagnostics, preserve all contributor parents, reject injected atomic storage failures, replay archived historical records, fork advanced source events, invalidate equal-content/different-source work and inspect a bounded 1,000-column layout at desktop/mobile widths. Methylation browser checks additionally verify explicit source states, XbaI Dam and NotI/SmaI CpG blocks, ambiguous and unresolved supplier rules, independent IUPAC/reverse/circular coordinates, optional Type IIS cut substrate, original target navigation, exact equal-revision/different-document guards, option/cancel/unmount invalidation, bounded paging/resource evidence, unchanged source history, readable boundary/sensitivity labels, responsive checkbox rows and mobile width. Checked-plan browser checks additionally verify unchanged blocked circular digest products, blocked self-ligation and missing-product recording gates, unknown/hemi states, explicit legacy behavior, exact equal-content source/target invalidation, preferred event bindings and causal parents, archive/current-engine replay, literal repaired 14-nt and Dam-retaining 18-nt products, conditional HTML/PDF/candidate-6 exports, wrong independent targets, imported-summary nonauthority and desktop/mobile evidence. CRISPR browser checks additionally verify both built-in PAM layouts and declared custom geometry, exact circular/reverse coordinates, an independent selected-record Hamming/IUPAC oracle, intended-only tagging with additional perfect matches retained, focused reference-label edits, every search limit and paging, literal independent-target and reverse/circular donor checks, fresh release gates, atomic IndexedDB rollback, source-only causal parents, native unverified producer provenance, cancellation/equal-content identity and desktop/mobile layouts. Large-reference browser checks additionally stream more than 32 MiB of canonical DNA, verify every owner and full canonical contig hashes, distinguish exhaustive imported-contig coverage from a missing declared chromosome, compare small-reference hits with an independent Hamming/IUPAC oracle, exercise both PAM sides and circular/reverse seams, actual guide/hit paging and navigation, both work/hit caps, cancellation, quota and corrupt/missing-block failures, native round trips, archive/reload without restored authority and exact-hash reattachment. Desktop window scrolling sends no scientific worker requests; corrected pagination labels and expanded/compact panel widths pass inspection. Provider mocking is explicit; it is not a live Addgene/API test. The published desktop MVP commit `33df2d8` has observed successful GitHub native tests, MCP builds and real stdio smoke on Windows/macOS/Linux, plus the separate WASM check ([run](https://github.com/henriksson-lab/collagene2/actions/runs/45384873995)). The post-MVP libraries and ninety-one-tool MCP checkpoint `77fcd64` also pass native tests/build/stdio on all three operating systems and the separate WASM check ([run](https://github.com/henriksson-lab/collagene2/actions/runs/45387324520)). The catalog/cache release `68980d9` also passed all three native OS jobs and WASM ([run](https://github.com/henriksson-lab/collagene2/actions/runs/45398885253)). Native desktop installations remain unverified. Firefox/Safari, hardware GPU execution and the remaining stated performance gates remain unverified.

Remaining work preserves the full plan:

- M0/M1: larger real-world interoperability corpora, real legacy producer-file validation after the source/schema audit, broader performance measurements including main-thread preflight/hash/serialization of large rich documents, richer graphical annotation editing, broader alphabet-aware editing labels, more configurable track layout/accessibility and browser coverage. GPU sequence glyphs are implemented with a same-frame 2D fallback; final browser/performance verification and a shared interactive/export scene remain pending.
- M2: fuller ancestral strand/qualifier annotation propagation, recorded source-event/operation links for standalone generated previews beyond retained full plans, codon designs and explicit historical forks, richer overlapping-primer/heteroduplex and overlap-assembly models, repeated-copy reaction alternatives, additional physical/end-repair states beyond ideal endpoint transforms, additional positional susceptibility rules, modified-oligo/other chemical models and supplier-condition data beyond the initial tracked profile, method-specific rules and bounded alternative-product exploration. Current ordered candidates do not establish uniqueness or empirical success.
- M3: broader MCP tools, persistent scoped transactional application, live authenticated bridges/hosted transports and cancellable progress/status.
- M4–M6: live Addgene interoperability, Gateway/TOPO chemistry, richer overlap/recombination and tailing-mixture models, curated shared feature catalogs, broader primer libraries/order formats, broader reference-free contig repeat/layout handling and real caller calibration datasets, broader trace/comparison UI, empirical ladder/catalog datasets and richer gel models, richer PDF/script shaping, protein processing/fusion models beyond the declared basic free-chain property profile, provider-derived reference build manifests, streamed remote reference import, persistent genomic indexes and real build/performance corpora, global edit-patch validation, variant/bulge search and calibrated activity/specificity datasets plus broader editing mechanisms, richer F26 host/codon-pair/RNA objectives, protected ranges and repeat/structure/supplier constraints, batch/library work, collaboration and Tauri packaging/testing on all three operating systems.

Aligned FASTA, CLUSTAL and Stockholm exchange is independently implemented in the core with strict bounded JSON retention, exact original bytes, row/column mappings, metadata and explicit export losses. The web now imports explicit-alphabet collections, displays a two-axis virtualized column/gap viewer with readable mapped coordinates, reviews export loss, opens selected ungapped rows and retains collections through full workspace exchange. Schema-10 workers and MCP expose the same typed exchange; recognized alignment attachments validate before archive/storage installation. [Alignment exchange](ALIGNMENT_FORMATS.md).

Primer3, multiple alignment engines/ports and OpenGel retain their deferrals. ELN/inventory/freezers/sample/barcode/stock management, laboratory automation/scheduling and organizational reporting remain outside scope.


## Selective junction preparation and public workflow research

The junction editor supports selected PCR regeneration while retaining unaffected
manual PCR and catalog digest preparations exactly. Variable overlaps, explicit
Apply and protected Undo remain source-bound; retained choices can still fail the
fresh complete-recipe verifier. Worker45/MCP131 provide the same advisory contract.
Official Benchling assembly and primer/mutation/validation workflow reviews are
recorded in the project plan. They guide interaction work, do not establish GUI
acceptance, and introduce no new browser-test gate. Further Golden Gate junction
editing, draft persistence and remaining cloning/bulk workflows remain active.

Validation for this checkpoint: eight focused native tests passed (five independent
core fixtures, one append-only worker/binding fixture and two real MCP SDK groups).
The WASM web target passed `cargo check`. Checks ran in the existing shared cache
while implementation/research continued. No browser checks or preview refresh ran.


## Local preparation recovery and persistent primer pair draft

The cloning dock shows typed per-row preparation progress, first failed operation,
blocked later steps and actual profile junction evidence. Direct Edit preparation
focuses the existing row without changing its contents; source navigation retains
exact context and pinned original regions. The same read-only locations are added
to existing MCP overlap/ligation responses and full evidence artifacts. Complete
reconstruction with an invalid independent target remains distinct from a failed
preparation, and an incomplete result is never promoted by these locations.

Unfinished shared primer pairs are saved as bounded native attachments with their
names, raw domain/tail text and exact historical selection anchors. Native export
and browser autosave use the existing workspace path. Restore verifies actual
original oriented bases and source events, preserves intentional edits and does
not restore prior PCR/mutation results. No new native/worker/science version or
external runtime dependency is introduced.

Validation: the final unified native batch passed 40 tests (five core diagnostic
fixtures, two actual MCP SDK groups and 33 native web cases, including the three
new pair-persistence cases). The full development Trunk build passed for web, worker and asset packaging.
No separate repeated WASM compiler pass or browser check ran. Port 8082 now serves
`/tmp/collagene2-diagnostics-public-dist`: 39 immutable assets / 175,726,276 bytes,
retaining previously versioned workers for existing tabs. Index SHA256:
`9e1828f66d6e8d1cc938ba931aaaeb7b4f61e31267dc0058443efaf1a2ad4309`; current worker bundle:
`5338dcc50e30f7b717c2ad6620b9cf5b961c5be0986666a03902f9be88b4ce58`. The server exposes no writes or filesystem API. Detailed local evidence:
`/tmp/collagene2-diagnostics-public-evidence.json`. This is build and publication
evidence, not GUI acceptance or completion of the product.


## Selected Golden Gate junction preparations

The assembly dock's existing junction editor now supports Golden Gate alongside
overlap designs. Source-core trims and literal inserts define independent intended
DNA first. Selected originals produce explicit PCR→pinned Type IIS digest
proposals; tail segmentation distinguishes actual own-5′ fusion and inserted
bases. The core reconstructs real PCR products, eligible cuts and oriented strand
ends before reporting new preparations. Retained rows preserve named oligos,
catalog pins, digest choices, repair and end treatment byte-for-byte. Their
compatibility is decided by the later fresh complete-recipe verifier.

Apply installs selected replacements and edited target declarations together,
reuses existing actual-product/diagnostic/output workflows and keeps protected
Undo. Unedited Golden Gate suggestions cannot overwrite an applied edited intent.
Last-linear selected PCR rows require an explicit terminal fusion; circles use
the first adjusted core's prefix. No terminal choice is invented for a retained
last row. The shared intent supports up to 200 normalized inserted bases per
junction (400 raw bytes), while complete new ordered oligos must fit 200 bases.

Worker46/job112/result104 and MCP132 provide the same advisory route. The older
Golden Gate DTOs and public empty-insert behavior remain unchanged. Candidate11,
plan9, science14 and archive1 remain unchanged. No external scientific program,
fusion optimizer, internal-site domestication or cycling-mixture solver is added.

The focused native batch passed 20 tests: five new scientific junction fixtures,
one new worker binding/selector fixture, two real MCP SDK groups, five existing
Golden Gate designer cases, two existing designer protocol cases and five existing
junction-intent cases. They include independent 42-base targets, true reversed
strand preparations, BsaI/SapI geometry, exact retained-row identity and full
fresh target mismatch checks. No browser checks ran.

The web WASM target also passed its background compiler check. Native/WASM checks
used the existing shared cache; no full Trunk rebuild, browser check or preview
refresh was run for this checkpoint. Port8082 retains the previous successful
immutable diagnostic/pair-draft snapshot until the next suitable publication.


## Prepared Golden Gate alternatives

The Golden Gate method evidence can now inspect a finite pool of explicitly
selected prepared roles, using raw physical strands retained privately by the
same fresh preparation verifier pass. The pool engine redigests those substrates
and retains self/subset backgrounds, alternate paths, recut sites, target
correspondence and explicit search caps. Missing material withholds the pool.
Discarded digest pieces, unrepresented original templates/amplicons and their
carryover are outside this declared scope. This is not a complete original one-pot
reaction or experimental-yield model.

Candidates and recut sites are paged without hiding reported alternatives. Exact
context/index/directed+physical keys select an unverified native plan retaining
all original inputs and preparation/treatment steps. Only the final operation
changes from ordered Ligate to existing GoldenGate. Review and output creation
use the existing fresh-verification path; no pool duplex becomes an imported
physical substrate or automatically installed DNA document. The exact selected
GoldenGate path's junction-window adapter remains pending; existing actual-product
sequence/map review and final verification remain available.

Worker47 adds inspection113/result105 and selection114/result106. MCP134 adds
inspection plus atomic selected-plan/binary publication using the same original
source/target bindings. Existing candidate11, plan9, science14 and archive1 formats
are unchanged. Ligation's public preview still uses the same single native pass;
its internal helper now makes the actual prepared raw strands available privately.

The focused native batch passed 15 cases: five new scientific prepared-pool fixtures,
one two-way worker selector/binding fixture, two real MCP SDK groups and seven
existing ligation/worker regressions. Independent literal targets distinguish the
full desired construct from self-products; fresh selection verification confirms
the difference. Recutting, missing preparation, metadata-only physical claims,
truncated search, raw orientation and stale/key tampering remain explicit.

The web WASM compiler check also passed. Both checks used the existing shared
cache while documentation work continued. No browser check or preview rebuild
ran; port8082 retains the last published immutable snapshot. Further one-pot
reaction coverage and selected-path junction windows remain active work.


## Actual Golden Gate path junctions

The existing junction inspector now supports the actual final GoldenGate path.
A private same-pass verifier receipt supplies the selected fragments, true source
cuts, actual opposite-strand orientation and final physical state. The inspector
checks product coverage and reuses existing DNA-window, frame and CDS calculations.
The UI labels selected path parts explicitly and offers a source/output mapping
view; MCP includes the same typed path and retains full evidence. Search caps and
overall target status remain distinct from junction-window coverage.

Worker48 appends nested method3 while preserving job96/result88, all old method
layouts and candidate11/plan9/science14/archive1. No scientific operation or second
search pass is added.

The native batch passed 23 cases: four new scientific/protocol fixtures and one
real SDK group (forward, actual reverse-bottom and wrong-target cases), plus 18
existing junction/worker/Golden Gate regression cases. Independent literals check
cuts 11/15, product spans, circle closure, boundary-crossing codons, capped-search
status, target-only annotation exclusion and failed-final-product withholding.

The background web WASM compiler check passed as well. No browser check, full
preview build or publication ran; port8082 retains the prior immutable snapshot.
Whole-reaction Golden Gate preparation/cycling and the remaining cloning product
requirements remain active.


## Whole-reaction Golden Gate composition

The shared alternatives panel now supports explicitly declared whole original
circles/linear DNA, complete PCR amplicons and current isolated preparations.
It is source-bound and available before ordered-method evidence. Whole-input
inspection retains all cut pieces and reports uncut circular carryover separately.
Native preparation plus actual enumeration execute once; index0 diagnostics retain
the pool even when that probe is unavailable or recuttable. Selection creates no
positive authority and changes only the final native candidate index.

Worker49/MCP136 expose the new typed request/report while old prepared-pool formats,
candidate11, plan9, science14 and archive1 remain unchanged. Complete cleavage and
one-copy limits remain assumptions, with broader kinetic/partial/repeated-copy
reaction modeling, optimization/domestication, bulk design and unfinished form
persistence still active.

The final native batch passed 25 cases: three new scientific/protocol fixtures,
one actual MCP SDK group and 21 existing Golden Gate, prepared-pool and junction
regressions. The literal 40-base target is declared independently. Tests retain
six input digest/carryover records, distinguish self-products and recutting,
preserve complete amplicons and verify selected plans freshly. The empty capped
case explicitly uses minimum two fragments, since a one-state search can already
find a one-fragment self-product.
The web WASM check passed. No browser check or preview publication ran; the public
8082 snapshot remains the last successful immutable build.


## Connected private-save and saved-template TOPO workflows

Current source separates private native saves from public sharing. Local catalog
snapshots retain exact original attachments and the complete bounded workspace;
opening restores only the selected native dependency graph plus its required
private documentary sources. Browser transport rejects private workspace mutation
payloads before serialization. Public selected-document and explicit SBOL subset
review paths remain separate. The native local capture admission check now accepts
its intended whole-workspace scope; public admission still requires document
closure. [Contract](PRIVATE_LOCAL_LIBRARY.md).

The candidate editor now inserts and edits standalone3′ resection with explicit
bounded counts while preserving its selected source. Saved-run ancestry labels
identify incomplete annealing references instead of implying reusable material.

Saved-native-template TOPO inspection replays the originating run before new PCR,
keeps vectors as explicit original preparations, and uses a separate worker89
wrapper binding the full source/capsule/task context. New PCR oligos are retained
for ordinary primer collection/export; inherited source-run oligos are not
reordered. Origin incompleteness and bounded finding-excerpt counts remain visible
in the UI and tool181 `inspect_saved_topo`. Candidate output is still freshly
reconstructed and checked against its independent target. Thirty focused native/static cases and integrated UI/worker WASM checking pass.
New cross-platform CI is pending the commit. [Workflow](TOPO_SAVED_TEMPLATE.md) ·
[MCP](MCP_SAVED_TOPO.md).


## Native publication reviews and operator backups

Native Databases controls can review/queue a new local folder, a single rename or
move, or deletion of a known locally removed bound empty folder. Parent bindings,
remote configuration, local generations and checkpoint are rechecked. Root and
selected-subtree mounts retain their acknowledged remote placement unless an
actual relocation is declared. Deleted-folder choices have checkpoint-bound
pagination, including access beyond item200. Queueing remains separate from
network resume and durable receipts preserve retry identity.

Saved native SBOL archives can now use structured reachable-sequence choices and
an explicit omission-count/hash review before queueing the public subset. Private
originals and native identities remain unchanged. Approval objects stay native,
source-bound and non-deserializable; bounded review IDs expire on session reset,
explicit discard, successful queue or replacement by newer reviews. Cancellation
cannot retain an undisclosed prepared token. [Metadata](REPOSITORY_METADATA_PUBLICATION.md) ·
[Native SBOL review](REPOSITORY_PUBLICATION_REVIEW.md).

The shared-library operator CLI now creates immutable WAL-consistent backups,
verifies database/schema/BLOB/history/receipt/checkpoint integrity and exact lab
scope, and restores only into a new private directory. It invokes no external
SQLite program and never overwrites a live database. Existing destination
replacement, automatic cutover and desktop backup controls remain later work.
[Backup and recovery](LAB_LIBRARY_BACKUP.md).

Thirty-seven focused cases pass across repository, raw host protocol, review/UI
models, backup/restore and related compatibility gates. Integrated UI/worker WASM
and native SDK checks pass. No browser/application HTTP checks or live-user data
backup/restore operations were performed.

## Ligation outcomes and Tauri recovery integration

The current source has one `Ligate` operation with explicit outcome options.
TA/GC, blunt and ordinary restriction ligation retain the intended order and
orientation while exposing bounded one-copy alternatives. Enumeration never
chooses the target. Fresh material remains private; complete outcomes are bound
to recorded products and history. HTML exposes JSON/strand FASTA; PDF includes
printable outcomes and embedded JSON/FASTA. Other copy numbers, vector-only
self-ligation and partial assemblies remain outside the declared scope.

Obsolete development codecs, compatibility conversions and version ladders have
been removed. Current internal schemas are worker95, candidate25, native plan23,
DTQ/DTR5 and DSQ/DSR3; scientific engine26 is unchanged. Unsupported internal
versions are rejected. External biological formats remain intentional features.
Current source passes417 focused cases, the affected native caller checks and
integrated web/worker WASM compilation. Cross-platform CI for this increment is
pending.

Desktop has one Bootstrap API with an optional genuine catalog context. Catalog
failure preserves the working editor. Exact pending-switch cancellation appends a
terminal journal generation and requires restart. Corrupt journals block catalog
access/cancellation while preserving editor startup. Repositories and backup
services retain their native session lease. Historical or refused cancellation
requests do not themselves pause a healthy catalog.

Saved-run fresh review, canonical SVG map export and their MCP tools are being
integrated. These preserve recorded evidence and require fresh verification before
producing new verified protocols; they do not rewrite historical assessments or
install candidate plans into the current editor.

The requested [SciLifeLab assessment](SCILIFELAB_LLM_ASSESSMENT.md) is committed
as `7eb3e22`. Native OpenAI-compatible HTTP planning is feasible; current model
capabilities, data terms and desktop/public-web entitlement remain unconfirmed.
No inference endpoint was contacted and no credentials or sequence data were sent.

Earlier CI exposed fixtures relying on obsolete development encodings. Those
positive compatibility fixtures are removed; retained tests exercise current
scientific behavior, integrity bounds and explicit unsupported-version refusal.


### Revision-bound lab discussion and human review

Native storage, individual server credentials, bounded binary discussion messages
and web revision-history controls are implemented. The operator registry supplies
trusted author/role; readers cannot mutate shared documents or submit reliable
sync writes. Threads, replies, resolution and review decisions remain immutable
and bound to exact stored revisions. Partial/stale pages cannot present a current
approval summary or enable posting. Human review never replaces deterministic
scientific verification or propagates to changed revisions.

All 29 native storage/backup cases and 18 binary SDK/credential/permission cases
pass. Both focused web-model cases and the final integrated WASM check also pass
(20 distinct native/model cases in this service/UI batch). HTTP/router fixtures are compiled but not executed. Graphical
acceptance and controlled service cutover remain open. SQLite schema 2 rejects
older schemas, and `serve` now requires `--members-file`; existing running storage
and services remain untouched. Offline/remote plasmid synchronization does not
silently transfer discussion or actor identity. See [lab collaboration](LAB_COLLABORATION.md).


### Desktop chosen-remote human review

The Tauri remote revision history now reuses the browser discussion/review controls
through native credentials and the exact saved remote/lab binding. Appends pin the
loaded actor ID, so a changed authenticated person cannot silently author the
previous person's review. Discussion remains explicit online work, independent of
local plasmids and their durable push queue. Discussion schema 2 and desktop sync
IPC schema 4 reject previous internal versions directly. All 25 focused native/model
cases pass, including identity and scope refusal, uncertain replies, full-page
bounds and native cancellation. The integrated WASM check also passes;
exact-source cross-platform CI remains open. No running service/database or application HTTP/browser check
was used for this batch. [Controls and scope](LAB_COLLABORATION.md).


### Collection sequence similarity in the worker, GUI and MCP

The shared affine collection-search profile now runs through worker114 (job191,
result181), Compare's expandable workspace search and MCP's immutable-reference
`search_sequence_library` tool. The current engine/report is
`sequence-similarity-2`; native tool inventory is190. Results retain original
reverse/circular coordinates, uncertain bases, unsearched subjects, origin/work
bounds and hit truncation. Context/shape guards check source revisions, dimensions,
score evidence, statistics, ranked rows and completeness consistency.

Web controls pin a query explicitly, select documents by stable IDs, preflight
borrowed models before cloning, cancel isolated workers on changes and navigate
exact captured events. Equal-content event changes are refused. Guarded JSON
export includes source captures; no document, candidate, material or run is
installed. All13 focused scientific/worker/SDK/model cases and the integrated
WASM check pass. Two fixture issues were corrected directly; passing scientific
cases were not rerun. Current-source browser acceptance and cross-platform CI remain open.
The implemented search covers explicitly loaded workspace/session subjects.
[Profile and controls](SEQUENCE_SIMILARITY.md) · [Remaining delivery triage](DELIVERY_AUDIT.md).

### Desktop catalog recovery CI regression

Cross-platform scientific CI exposed a production mismatch: native catalog schema
2 was still rejected by schema-1 checks in desktop startup/backup and recovery
review. One shared database-version constant now drives the store and protocol;
native read-only probes reuse the store's application/version guard. Unsupported
versions remain refused, with no migration or fallback. Working-state databases
and activation journals retain their separate formats. Running services and user
databases remain untouched.

All 53 focused native/model cases pass: 24 desktop unit cases, five real temporary
startup/recovery fixtures, 19 storage/backup/discussion cases, one binary receipt
case and four web recovery-model cases. The receipt fixture also rejects schema
0, superseded schema 1 and the next unknown schema. The integrated WASM check
also passes. Exact-source cross-platform CI remains a delivery gate.
[Recovery contract](DESKTOP_CATALOG_BACKUP.md).

### Imported structure coordinate and scene core

The new MIT `collagene-structures` crate implements bounded PDB and mmCIF
coordinate import with exact retained bytes/hash, separate models and identity
namespaces, alternate/insertion qualifiers and declared connections. The scene
preparer keeps source atom indices, normalizes f64 geometry into finite f32 GPU
positions and refuses ambiguous, conflicting or unsupported connection resolution.
No distance-based bonds, folding, symmetry expansion or sequence correspondence
are inferred. Deserialized projections must reparse from the exact current source.

Ten distinct authored fixture groups pass across the integration batch, including
fixed columns, quoting, nonpolymer identities, unknown/not-applicable/literal
selectors, category admission, original-byte retention, geometric expectations,
extreme finite coordinates and tampering. Only affected fixture groups were
rechecked after parser/selector corrections. The final crate WASM check also
passes. Native three-platform and WASM CI now include this crate. Worker, GUI,
workspace and MCP integration remain to be implemented; this core does not
complete F34 or establish graphical acceptance. [Profile](STRUCTURE_IMPORT.md).

### macOS collaboration fixture isolation

The published catalog fix passed Linux scientific CI and the web job. macOS
reached the discussion backup fixture but failed creating its temporary directory:
two concurrent tests received the same timestamp-derived name. The helper now
includes an atomic process-local counter, as the other native disk fixtures do.
No production behavior, recovery assertion or scientific passing criterion is
relaxed. Both affected disk-backed fixture groups pass locally; cross-platform CI
must confirm the correction on macOS.
