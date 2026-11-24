# Source-backed physical TOPO adapter

The `core::topo` adapter reconstructs prepared strands from exact immutable
`Input { role, molecule }` records and an explicit versioned `GenerationRecipe`.
It does not infer activation from a product name, vendor, annotation, provenance
comment or imported physical summary. It does not assert commercial kit parity.

A `Prepared` material declares named actual 5′→3′ strands, each assembled from
bounded source intervals, actual pairing spans, and explicit terminal chemistry.
Intervals may reverse-complement a reference region or cross one circular origin.
The declaration includes its own ID and human-readable preparation assumption.
All strands are derived from the source DNA; an arbitrary unbound product string
is not a preparation input. References are generation templates: repeated copying
is permitted as an explicit hypothesis, and source correspondence is not a claim
of atomic inheritance. Unrepresented reference intervals are reported separately.
They are not invented manufacturing waste or simulated byproducts.

`ExactPcr` is an alternative material recipe. It freshly runs the existing exact
PCR geometry enumerator on an original source with explicit primers and requires
exactly one geometry. It reconstructs the current ideal blunt OH-ended duplex,
then optionally invokes the existing explicit single-nucleotide 3′ tail transform.
The generated strand evidence records template interval, primer-tail lengths and
actual strand direction; added terminal bases are separately identified as new
synthesis. No inherited strand modification or polymerase terminal addition is
inferred. Imported graph-only intermediates cannot be used as ordinary templates.

`simulate` retains all prepared-source maps, unrepresented-reference evidence and
the entire physical TOPO engine report. Transfers use actual generated-strand
indices in declared material order. `inspect_orientations` accepts an empty
transfer request and explicit vector/insert strand pairs, computes both declared
orientations with the existing engine, and returns both full graphs plus their
exact candidate-specific requests. It selects no winner and claims no exhaustive
biochemical mixture. Pairing indices inside each prepared material
are local to that material. Reverse-oriented preparations require explicit reverse
source intervals and pairing; names do not assign a productive orientation. All
unreacted chains, reciprocal physical components and displaced calls remain.

`product` freshly reconstructs the recipe and selected physical component. Its
`Selection` pins adapter batch revision, component and an actual source-strand
call as reference seed. The independent paired-contour proof is required before
creating a reference molecule. This proof does not seal nicks or remove flaps.
A directional flap therefore remains a full graph with a withheld contour view.
A separate free strand remains in the graph even when another selected component
has a proven paired contour. Circular reference origin and sense follow the seed,
not an intended target. A linear origin follows the physical endpoint.

The optional `reference` molecule is a coordinate view with projected top DNA,
not a replacement physical substrate. It carries the existing generic
`source:legacy:chemistry_unmapped=true` guard. Full graph, actual bottom DNA,
source-call maps, activations, nick chemistry and end chemistry remain in the
product envelope. No parent annotations are silently copied into this initial
adapter; its source correspondence and immutable originals remain evidence.

`check_target` accepts a separately supplied target and explicit acceptance policy.
DNA/topology comparison is independent of component reconstruction, orientation
and origin. Rotation and reverse complement are opt-in. Nicked products, unknown
chemistry, unapplied transfers and other components each require separate opt-in.
A match can be accepted only under this declared graph policy; it does not prove
reaction success, repair, efficiency, stable short DNA, or a proprietary kit state.
Without a unique paired contour, comparison is withheld rather than extracting a
target-guided path through a directional flap.

Adapter, source, recipe, engine, proof and policy domains are pinned separately.
Changing source metadata changes immutable source identity as well. Imported
reports do not restore chemical authority. `replay_product` is a heavy fresh native
audit, not a main-thread validator or a pathway from serialized claims to a physical
input. A new source-based reconstruction is required to materialize a product.

Borrowed aggregate model guards run before rich validation/hash/copy. The initial
adapter bounds source DNA and prepared DNA to 100,000 bases, 32 input roles,
16 materials, 32 generated strands, 2,048 reference segments, 1,024 pair spans,
16 transfers and 4,096 unrepresented intervals. Recipe serialization is capped at
1 MiB; complete reports at 32 MiB. Existing engine/proof work, memory, component
and report budgets remain additional constraints. PCR enumeration is conservatively
preflighted by possible primer-pair count times maximum product length before
calling its existing enumerator. Limits withhold/error rather than claim absence.

Before any PCR scan, preparation also reserves an aggregate search allowance of
`4 × template length × (forward annealing length + reverse annealing length)`
per PCR material. This covers both orientations in the candidate-count pass and
again in the geometry enumerator. The recipe-wide sum must fit the requested
`engine_limits.max_work_cells`; the physical engine independently reserves its
own graph work under that cap. This is a conservative comparison-work bound,
not a wall-clock timeout or a combined preparation-and-graph work counter.

The adapter now connects to typed `Operation::Topo` plans, full reconstruction
JSON evidence, HTML/PDF protocol export and four MCP tools. Each role binds an
explicit immutable source. Plans freshly reconstruct the graph and retain it in
`VerificationReport.step_reports`; a paired reference view keeps the chemistry
boundary sticky through subsequent operations. A missing contour preserves the
graph while withholding DNA. The initial plan operation deliberately has no
implicit physical target-acceptance policy and remains Incomplete; standalone
`check_topo_target` requires a separately supplied target and explicit policy.

The browser can review and verify imported TOPO candidates and download the full
physical reconstruction evidence. Existing source-event run capture/replay
retains source recipes and immutable inputs, not a restored proof. Candidate
schema8/workspace plan schema6/worker schema25 and science10 identify the new
semantics; older worker verification layouts migrate without invented graphs.

Kit-specific preparation wizards, full annotation projection and downstream
physical graph composition remain implementation work. No contour preview can
be used as an intact-duplex substrate by silently dropping its graph. The
underlying engine and its bounded elementary-transfer model are documented by
`collagene-topo`. Compiler validation is separate from scientific/runtime
validation; browser checks remain deferred until in-scope code is complete.

The source-bound [web preparation workflow](TOPO_WEB_WORKFLOW.md) adds an explicit TA form and worker-backed orientation/projection/target inspection. Its recipe builder does not run science or choose a winner on the main thread.
