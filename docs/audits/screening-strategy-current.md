# Named-product screening foundation

Current integration source, 2026-11-28. `core::cloning::screening` is registered
with the binary worker, MCP and the Analyze menu's screening panel. Native/WASM checks and 23 screening-specific core, worker, SDK and web-model
cases pass in the [shared execution record](connected-design-current-integration.json).
No application or graphical acceptance check ran.

## Scope and inputs

The request contains exact-revision molecule components, named subjects, explicit
subject component lists and explicit candidate tests. Exactly one subject is
intended; every other subject has an explicitly declared decoy role. Labels such
as empty vector, wrong orientation, tandem insert or parent carryover **never
construct sequence geometry**. A carryover subject can contain intended material
plus its separately supplied parent. Identical DNA in different named outcomes
remains separately represented and may be nondiscriminable.

Existing `material::derivation::SourceOrigin` resolves each component through the
canonical original-reference guard or fresh native saved-origin replay. Original
modifications remain unknown unless an explicit `GenerationAssumptions` is
supplied for an original reference. A saved source cannot be relabeled by such a
declaration. Unresolved native source replay prevents a complete signature.
No imported physical metadata is removed or promoted into native authority.

## One current scientific implementation

* Digest candidates use the installed exact catalogue selection and
  `tracked_duplex::checked_digest_with_budget`, including native recognition,
  signed cut geometry, positional modification rules and fragment projection.
  Unknown cleavage is retained as unresolved, and every component/candidate
  shares the same actual `CleavageBudget`. A digest candidate may explicitly
  combine up to four selected enzymes; candidates are separate hypotheses.
* PCR candidates use `pcr::source::Template::fresh` and
  `geometries_with_budget` with one shared `ExactSearchBudget`. Native exact
  binding, both primer roles, circular coordinates, actual template-strand
  coverage, declared tails and all returned geometries stay in their existing
  owners. Finite enumeration completeness is separate from product uniqueness
  or the experimental amplification of any/all products. Native unsupported
  ambiguity or resource refusal never becomes a complete empty signal.
* Native catalogue provenance/primary numerical facts remain owned by
  `cloning::restriction::catalog` and `chemistry::tracked_duplex`; this owner
  imports no new parameter table, biological calibration or external algorithm.
  Existing `output::gel` migration/calibration/ranking is not used as experimental
  resolution evidence. Existing ligation outcome enumeration can supply actual
  explicitly selected molecules through host integration; no outcome is inferred
  by comparing DNA, names or lengths.

## Explicit ideal readout

Each complete component yields a sorted set of distinct **ideal modeled sizes**.
Digest sizes are integer twice-mean-strand-length values (top length + bottom
length); PCR sizes are twice product length. A subject mixture takes their union.
Equal-size fragments collapse: multiplicity has no band-count/intensity meaning.
No concentrations, intensity, detection limits, migration, yield or PCR efficiency
are supplied. An uncut circular digest is unsupported for this linear-size model;
it is not an empty signal. A complete zero-product PCR search is an empty modeled
signal. The latter may be nondiscriminating if another mixture component supplies
the intended signal.

A separation witness is one size in either complete signature having **no**
opposite-signature size within the user-declared closed tolerance. Different set
counts do not prove separation: `{100,104}` versus `{102}` at tolerance 2 has no
witness. This tolerance is a declared comparison policy, not calibrated gel
resolution. Unsearched, ambiguous, unsupported or limited components prevent
mixture absence/distinction claims. Partial data may remain visible as such.

Test ranking sorts proven-separated decoy count, complete comparison count and
original request order. Optional panel search enumerates combinations of one to
three tests, unioning only proven per-test distinctions. Actual visited comparison
work is charged, with complete-search/best-among-searched disclosed independently.
Equal coverage retains the first result in size-then-lexicographic order, favoring
fewer tests. No biological winner or full optimality claim accompanies a limited
panel search.

## Admission and identities

At most 16 components, 16 named subjects and 64 tests; aggregate source length is
one million bases, each digest has at most four exact selected enzymes, each
retained native cell has at most 64 sites/products. Source requests are at most
16 MiB and reports at most 8 MiB, with a smaller explicit report ceiling admitted.
Native shared cleavage work is at most 100 million and PCR work at most 200
million. Native PCR independently bounds pair enumeration and materialized DNA;
the wrapper's retained-product ceiling does not change that native algorithm.
At most 100,000 panel comparison-work units are admitted. Retained cells use a
streamed aggregate byte ledger; final report overhead and identity are checked
before publication. Refusal returns an explicit resource error rather than a
silently truncated report. Every requested cell remains represented when only
scientific search budgets exhaust. Saved-origin replay and direct tests share
the same exact-PCR and cleavage ledgers. Other algorithms in saved-origin replay
retain their separately documented native caps.

The context hashes the exact current request, including all source snapshots,
origins, labels, mixture lists, tests, scope, tolerance and limits. The report
retains source/origin pins, explicit test definitions, input-indexed cells, native
cleavage evidence/PCR geometry coordinates, named signatures, witnesses, ranking
and panel completeness/work. `validate_binding` checks bounded documentary
identity, original slot structure and recomputes only signature/rank/panel
aggregation; fresh `validate_for` replays native assessment. Documentary content
alone does not attest an actual workspace origin or experimental completion.

## Independent authored fixtures and next integration

Ten core source fixtures cover literal EcoRI stagger-derived sizes and parent
mixtures, separately named identical DNA, complete empty PCR versus uncut circle,
finite multiple PCR products versus retained-output refusal, shared exhausted
work with unsearched rows, exhaustive tiny-set tolerance witnesses, independent
small-subset panel optimum/completeness, and stale/tampered/bounded requests.
The small-set oracle computes closed-interval presence directly, while the panel
oracle enumerates bitmask subsets rather than the production combination loops.
Additional fixtures cover integer-boundary tolerance and self-resealed evidence
that exceeds work limits or contradicts actual PCR geometry. Three binary
fixtures cover worker/session dispatch, exact context and source pins, current
schema refusal, altered report bodies and limited-search uncertainty. Four MCP
fixtures and six web model fixtures cover retained source selection and adapters.
All 23 registered screening cases pass. The SDK helper was corrected to send
an immutable workspace reference; its four cases reran, including a rejection
case whose former pass could have had the wrong cause.

The web panel captures exact retained events and native material origins, with
request/settings/workspace-lifetime guards before accepting responses. Named
subjects and mixtures remain explicit. It does not turn ranking into an automatic
choice or reuse documentary results as native material authority. The owner has
no separate plan/protocol encoding and writes no molecules or history. Final user acceptance of the finished UI remains.
