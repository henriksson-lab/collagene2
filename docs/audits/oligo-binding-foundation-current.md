# F38 oligo collection binding foundation — 2026-11-28

This independent crate addresses the collection binding gap in
[BEYOND_COMPETITION](../research/BEYOND_COMPETITION.md), with no borrowed cloning
application code. Existing `core/src/cloning/pcr/primers.rs` owns fixed-domain
primer binding and the perfect-duplex nearest-neighbor model; mutagenesis owns
physical product reconstruction. Neither owner is replaced by this foundation.

The new explicit model is ungapped, physical-3′-anchored suffix extension, not a
gap-aware aligner. An incremental recurrence considers every eligible suffix
length once per endpoint; independent brute-force fixtures rescore full suffixes
from scratch. Every best-score concrete tie is retained or the endpoint is
explicitly omitted under a limited result. For ambiguity, only a certainly
admissible competitor can eliminate a candidate using score bounds. Fixtures
expand tiny IUPAC inputs and require every realization's concrete optimum in the
retained conservative envelope. Scores are alignment ranking only.

The primary thermodynamic reference [SantaLucia 1998,
“A unified view of polymer, dumbbell, and oligonucleotide DNA nearest-neighbor
thermodynamics”](https://pmc.ncbi.nlm.nih.gov/articles/PMC19045/) motivates keeping
sequence-alignment ranking distinct from context-dependent thermodynamics. This
crate does not implement mismatch thermodynamics or call its score an energy,
Tm, binding probability, extension efficiency or PCR success prediction. Existing
perfect-duplex calculations must be unavailable for unsupported mismatch contexts.

Integration inputs: `Collection`, `Options`, reusable `Budget`, `search`, current
`Report::validate_binding`, individual `Hit::validate_columns`. The report binds
source DNA/topology and ordered collection/options hashes. Core/worker integration
must add current source document/revision, collection snapshot and lifecycle
acceptance. UI can display mismatch ranges, exact footprint, tail vs extra unpaired
bases, 3′ pairing evidence, per-oligo completion and explicit limited/unsearched
states. It must not interpret a limited or unsearched zero-hit result as absence.

The crate is now registered with its core wrapper, worker schema 130, Analyze menu
panel and immutable-document MCP tool. The wrapper uses the existing perfect-duplex
thermodynamic owner only for supported concrete pairing. Whole-oligo Tm is labelled
as a hypothetical perfect complement, including the ordered tail; mismatched or
uncertain actual footprints report unavailable Tm. The existing model's typed
decode-resource admission hook is public for bounded decoding.

Shared final native/WASM compilation passes. The worker ran no compiler, test,
browser or application HTTP checks. Fourteen independent-library fixtures, core
transport cases, six web input-model cases and actual in-memory MCP SDK cases are
all pass. The budget owner exposes minimum source/collection admission cost;
insufficient admission is refused, while admitted scans stopped before examining
endpoints report an explicit work limit. This distinction is tested at native,
transport and actual MCP SDK boundaries. [Execution evidence](design-checks-current-integration.json).
