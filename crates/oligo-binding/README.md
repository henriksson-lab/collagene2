# Oligo binding foundation

MIT Rust library for bounded DNA oligo-collection sweeps. Native and WASM use the
same pure-Rust algorithm. Core adds source-document identity and supported
thermodynamics; the worker, web Analyze panel and MCP share that wrapper.

`search(source, topology, collection, options, &mut budget)` returns explicit
ungapped, physical-3′-anchored annealing suffixes on both strands. Each endpoint is
extended incrementally in O(oligo length), including unfavorable intermediate
extensions. Concrete inputs retain **every admissible maximum-score suffix tie**
at that endpoint/orientation. Scores rank sequence alignments; they are not
binding efficiency, probability, PCR yield or extension predictions.

IUPAC inputs produce minimum/maximum mismatches and score bounds. An ambiguous
candidate is retained unless its upper score is below a *certainly admissible*
competitor's lower score. This conservative envelope includes possible concrete
optima, but can contain candidates that never win under any consistent concrete
realization. It is not an exact ambiguity-resolved optimizer. Physical 3′ pairing
is reported as certain/possible runs. The default required run is zero, so terminal
mismatches are visible. An explicit caller policy can require possible pairing.

A separately declared 5′ tail is never aligned. Additional unpaired 5′ sequence is
reported separately. Aligned columns cover the whole annealing footprint in oligo
5′→3′ order, with original source coordinates, source strand bases and actual
antiparallel partner bases. Circular hits traverse at most one source circumference.
No indels, bulges or gap-aware alignment are implemented.

A reusable `Budget` accounts input hashing, complete endpoint extensions, retained
hits and conservative output allocation charges. Bounds: 1,024 oligos, 200 bases
per oligo, 200,000 collection bases, 200 million work units, 10,000 retained hits,
16 MiB allocation charge. All ties for an endpoint are admitted atomically before
columns/IDs are cloned. Metadata/input allowance failure is an explicit error.
Other exhaustion gives a limited current oligo and unsearched later oligos.
Only a completed oligo with zero hits supports a zero-site result **within this
explicit model and policy**. Shared budgets never reset implicitly. Work units
are documented algorithmic charges, not wall-clock time or exact CPU instructions.
`Budget::admission_cost(source, collection)` exposes the minimum work and output
charge needed to bind inputs before examining endpoints. A smaller allowance is
refused; an allowance equal to that cost can report an admitted but unsearched
endpoint set, with an explicit work limit rather than a negative binding result.

`Report::validate_binding` checks current engine/model, source/collection/options
hashes, per-oligo accounting and recorded columns. `Hit::validate_columns` checks
individual alignment facts independently of the extension recurrence. Neither
method authenticates an untrusted completeness/optimality claim; recompute in the
trusted worker for current claims. Source hashes bind DNA plus topology, not
molecule revision, feature metadata, physical material or experimental authority.
Core integration must additionally bind the current document/event, collection
snapshot and worker lifecycle generation before accepting results.

Public DTOs use Postcard-compatible external enums. Optional `json-schema` enables
MCP schemas. Report decoding bounds individual hits and aggregate retained columns
(65,536), including native typed resource admission through the existing model
decoder. Callers must still impose an outer transport-frame byte
limit, validate inputs and validate current report bindings. General serde formats
can allocate strings before their bounded string visitor receives ownership.

Thermodynamics remain in the existing core owner. Its integrated wrapper reports
whole-oligo and actual-annealed-stretch calculations as distinct, context-bound
results. Mismatched/ambiguous duplexes unsupported by the current thermodynamic
model must report unavailable; do not apply a perfect-duplex Tm as mismatch Tm.
No external software, third-party alignment runtime or experimental claims.

Authored tests compare the production recurrence to an independent exhaustive tiny
concrete footprint oracle and ambiguity-realization envelope oracle. They also
cover ties, terminal mismatches, circular/reverse coordinates, tails, work/output
exhaustion, context tampering and transport admission. The integration owner runs
these together with native/WASM checks in the connected integration batch.
