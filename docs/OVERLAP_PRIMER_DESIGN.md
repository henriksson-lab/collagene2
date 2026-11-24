# Deterministic PCR overlap-primer planning

`core::overlap_primer_design` plans one anchored PCR primer pair for each ordered
original-source interval, then returns a practical `assembly_recipe::Recipe`.
Selection uses the existing primer engine and homology profile. It is advisory:
no field claims PCR uniqueness, experimental compatibility or target verification.
The original source molecules, annotations, citations and event references are
never replaced by extracted fragments or mutated during planning.

`Request` specifies a name, ordered `Fragment` rows, a common overlap length,
linear/circular product topology, existing `primers::DesignOptions` and existing
`overlap_profiles::ProfileOptions`. Each fragment has a distinct role ID, original
input ID, zero-based `start`, unrolled exclusive `end` and assembly `reverse`
orientation. `start < source_length`; an interval contains at most one turn.
Only an actually circular source may have an unrolled end beyond source length.
The requested core must be concrete DNA, exceed the homology length and permit
both nonoverlapping endpoint annealing domains. Reusing an original under distinct
roles is supported; this does not establish reaction copy number or uniqueness.

## Tail strategy and origin

Let each **oriented** requested core be `C[i]`, and let `H` be the declared overlap.
Every upstream core with a successor receives `prefix(C[i+1], H)` on its right:

- Forward row: its forward primer tail is empty; its reverse primer 5′ tail is the
  reverse complement of that next oriented prefix.
- Reversed row: its original-source forward primer 5′ tail is that reverse
  complement; its reverse primer tail is empty. The assembly reverse flag acts
  after PCR, moving the synthetic extension to the oriented core's right end.
- The last linear row has no synthetic tail. The last circular row extends toward
  the first oriented prefix. Trimming the explicit closing homology preserves the
  first requested core's origin.

Every returned row remains an `ExactPcr` operation on the **original** input.
Each junction in the recipe uses `H`; sequence is not borrowed from the target.
For example, a next prefix `GCTTACCG` gives a primer tail `CGGTAAGC` in the relevant
original-source primer. The downstream core supplies the matching overlap once;
assembly removes the duplicate arm and retains the concatenated requested cores.

Existing `design_primers` ranks the requested annealing-length range under explicit
Tm, GC, homopolymer and contiguous-complementarity filters. It receives the actual
source interval, calculated tail positions and `max_pairs = 1`. The request's
preexisting forward/reverse tails must be exactly empty; supplied sequence is
never silently discarded. Request-level target endpoints and output count are
explicitly overridden by each row and one pair; other settings are preserved.
Returned primer names and analyses use `role-forward` and `role-reverse`.
No passing pair is an explicit error rather than an empty or partially usable plan.

## Independent target and evidence

`design(inputs, independent_target, request)` returns:

- `context_revision`, binding the engine, ordered input ID/revision tuples, actual
  target revision and complete request.
- `source_revisions`, aligned with original input order.
- The original request, generated recipe, one `DesignReport` per requested row
  and explicit scope warnings.

The independently supplied target is never used to select endpoints, rank primers,
choose orientations or synthesize tails. A target mismatch therefore does not
change the proposed oligos. It changes the context binding and must be detected
by fresh verification of the compiled recipe. A best-ranked pair can have repeated
or alternative PCR geometries; the existing verifier still enumerates those and
withholds a unique physical preparation when appropriate. Primer ranking does not
replace that independent scientific check.

`preflight` accepts borrowed `(input_id, &Molecule)` pairs before host cloning.
`design_revision` validates actual source revisions and computes the bounded
context binding. `Report::validate`/`validate_binding` check returned row/order,
source revision membership, names, domain/tail equality, requested filter settings
and recipe consistency without rerunning primer selection or PCR. A decoded
report restores no private chemical state and authorizes no document edit.

## Bounds and unsupported input

There are at most 32 inputs/rows and homology is 1–80 bases. Existing borrowed
native-model limits apply, with at most 8 MiB total original-source-plus-target
DNA, the inherited generated-recipe bounds and a 32 MiB canonical request/output
ceiling. Source metadata and actual primer endpoints are checked before copying
original documents. Physical/reference-only/variant producer evidence and
single/mixed-strand intermediates cannot become new PCR substrates through this
planner. Ambiguous requested core calls are refused; binding elsewhere in a source
remains the verifier's separate question.

The existing primer filter profile allows annealing lengths 8–80, with a bounded
range of candidates. Candidate structure work is conservatively charged as
`8*K*(A²+B²)` and pair alignment work as `K²*(A+B)*max(A,B)` per row, where `K` is
the requested count of lengths and `A/B` include the larger possible synthetic
tail. The aggregate ceiling is 100 million units. The existing per-row engine's
20-million pair reservation also applies. Later exact-PCR scan reservation is
precharged as `8*N*max_annealing_length` per row, at most 100 million aggregate
comparisons; existing PCR alternative/product bounds remain authoritative.
No external executable, alignment engine, Primer3 or target-guided DNA search is
called.

Focused source fixtures use independently written linear, circular, reversed and
origin-crossing targets; compare originals before/after; reject tails/bad ranges,
stale contexts and failed filters; and retain repeated-PCR ambiguity for fresh
verification. Worker and host fixtures separately bind the source/target/request
through the ordinary typed transport.
