# F37 assembly-standard foundation and first supported profile

Integration checkpoint, 2026-11-28. The new
[assembly-standards crate](../../crates/assembly-standards/README.md) is registered
with its core wrapper, worker schema 130, Analyze menu panel and MCP catalogue and
assessment tools. Shared final native/WASM compiler checks pass. No compiler,
tests, browser or application checks were run by the source-authoring worker.

The first installed-data definition is **BioBrick RFC10 part-sequence format
(2007)**, not a generic grammar alone and not full RFC10 compliance. Primary
evidence is [Knight's MIT-deposited draft](https://dspace.mit.edu/entities/publication/764237ec-294e-4c93-b224-2703800bf358)
and [the author's OpenWetWare standard](https://openwetware.org/wiki/The_BioBricks_Foundation:BBFRFC10).
The original 2007 text includes mandatory NotI absence, distinct coding/noncoding
prefixes and advisory sites/termination rules; these are preserved rather than
replaced by an unverified modern four-site summary. Its exact scope/rights and
revision appear in [the data notice](../../crates/assembly-standards/data/rfc10-2007/NOTICE.md).
An original raw-text snapshot hash is not claimed after a DSpace HTTP429 response.

Implemented source:

- Data-driven slots, declared local successors, start/end roles, occurrence/total
  cardinality, enzyme IDs/exact geometry, left/right additions and overhang
  signatures, body constraints and default/extra avoided sites.
- Immutable validated `Profile` bound to catalog edition/content and a derived
  profile identity; source/revision/license metadata remain present. Reserved
  built-in identity accepts only value-exact installed facts. Modified rules must
  use a declared `Custom` profile. No caller `verified` flag exists.
- Bounded sequence-format assessment, both-orientation motif scanning with IUPAC
  certainty, exact formatted flanks and separate advisory recommendations.
  Exhausted work/findings/output budgets remain explicit; partial evidence never
  establishes conformance. Reports retain exact declared scope, profile kind,
  source evidence, sequence/order/body/limit binding and fresh replay validation.
- Postcard-compatible external enums plus optional JSON-schema derives on public
  transport DTOs. Passing fixtures include native binary roundtrips.

No functional promoter/CDS/terminator order is invented. No duplicate physical
assembly simulator, cleavage prediction, ligation success model or external tool
is introduced. Existing core target reconstruction/verification must remain the
physical authority when this crate is integrated. Grammar hints are local,
not exhaustive future-path reachability or chemical compatibility proofs.

Seven source-authored integration cases cover literal required/advisory facts,
coding prefix/start/stop distinctions, reverse and ambiguous sites, custom
cardinality/order, built-in/catalog tamper refusal, bounded incomplete outputs,
fresh binding and Postcard roundtrips. All seven pass with the connected native/
WASM checks, core transport and actual MCP SDK cases. The core wrapper binds original
source intervals, orientation and formatted-body coordinates; findings can select
the corresponding original circular/reverse spans without inventing assembly
history or chemical compatibility. [Execution evidence](design-checks-current-integration.json).

Other toolkit tables remain unadmitted. The existing research document's advice
to consult MIT repositories is a source lead, not evidence that every toolkit
table/sequence is licensed or complete. Primary/data-specific review is still
required for MoClo, YTK, EMMA, PhytoBricks, SEVA and RFC1000. RFC10 original prose
is OpenWetWare CC BY-SA/GFDL; this independent implementation/factual selection
does not copy that prose or any competitor's code.
