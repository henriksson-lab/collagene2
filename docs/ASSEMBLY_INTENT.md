# Explicit frozen assembly intent

`core::assembly_intent` builds an intended reference from explicitly selected
original-source regions, order and orientation. It is independent of every PCR,
digest, primer-ranking and assembly simulator. There is no target argument,
preparation operation, predicted physical product or experimental status.

The workspace first captures this intended reference privately. Later primer or
assembly planning receives it as a separately frozen target; those planners cannot
choose its calls to make their own output pass. Host code owns exact source
document/event bindings, stale-selection guards and any explicit publication.

## API

`Request {name, fragments: Vec<overlap_primer_design::Fragment>, circular}` uses
existing region rows: `role`, `input_id`, zero-based `start`, unrolled exclusive
`end`, and `reverse`. Roles are distinct stable contribution identities; the same
original may appear more than once under different roles.

- `preflight(&[(&str,&Molecule)], &Request)` guards borrowed originals before host
  copies, canonical hashing or coordinate extraction.
- `intent_revision(&[PlanInput], &Request)` verifies actual original revisions,
  then binds the engine, ordered input ID/revision tuples and full request.
- `build(&[PlanInput], &Request)` concatenates exactly those declared reference
  intervals, applying reverse complement where requested.
- `Report::validate` / `validate_binding(expected)` check readonly scope, request,
  contribution placements, sequence receipt and context without replaying a
  scientific assembly. Decoding a report restores no material or edit authority.

`Report` returns `context_revision`, the request, intended `target`, ordered
`PlacedFragment` records and warnings. Each placement retains role, original ID
and revision, requested start/end/reverse, plus half-open `output_start/output_end`.
Those boundaries cover the intended reference exactly and support a main-pane
source contribution preview.

## Coordinates and independence

`start < original_length`; the interval is nonempty and contains at most one
source turn. Linear ends stay inside the actual source. Circular origins may use
an unrolled end, for example source `AACCGGTT`, interval `[5,10)`, reverse, gives
`TTAAC`. Source orientation/topology is never changed to obtain this output.
Output topology is exactly the declared circular/linear intent. IUPAC notation
is preserved with an explicit warning; it is not resolved into concrete bases.
RNA/protein-to-DNA conversion is not inferred.

Readable reference sources remain eligible even when their metadata carries
physical or reference-only provenance. Reading coordinates for design intent is
not a claim that the source can serve as a PCR/ligation substrate. Later physical
preparation guards remain unchanged.

The target carries `assembly:intent:reference_only=true` and an honest textual
`duplex_model` descriptor: `Reference-only assembly design intent; physical material
state unavailable`. This invokes existing restrictive metadata guards without
fabricating JSON chemical state or labeling intent as a legacy/recombination
product. Other metadata records the intent engine/context, contribution mapping,
sequence integrity receipt and scope. A sequence receipt detects accidental
inconsistency; it does not authenticate a report or prove source ancestry.

## Source annotations and raw evidence

Existing explicit-segment projection maps source coordinates into the selected
output contributions, including reversed and origin-crossing intervals. Fully
retained safe features preserve their biological extracted sequence. Clipped,
interrupted, fuzzy or unresolved annotations remain reviewed/unmapped under the
existing projection policy; translation/phase claims are never fabricated.

Unsafe projection or output annotation-budget exhaustion withholds that role's
annotations with a notice, preserving the intended DNA. Counts identify reviewed,
unprojectable and outside-contribution features; original IDs/revisions identify
where full features and raw evidence remain. The original input molecules,
metadata, annotations and citations are untouched. The new concatenated target
receives no active original GenBank/citation record pretending to describe the
whole design. Native workspace/source-event retention is required for full
original evidence; a standalone intended-reference report does not embed all
original files or prove physical molecular ancestry.

## Bounds

At most 32 originals and 32 contribution rows are accepted. Original DNA and
output DNA each have an 8 MiB aggregate limit. Standard borrowed native-model
limits and a 32 MiB canonical request/output ceiling apply. Intervals are checked
before extraction/allocation. The target annotation count is at most 4,096 and
aggregate native feature/qualifier bounds apply before copying projected features.
A conservative 30 MiB annotation-copy serialization budget reserves space for
request/provenance/warnings; excess annotations stay with the originals and
produce an explicit notice. Warnings are bounded and abbreviation is identified.

No primer design, digestion, overlap matching, physical materialization,
experimental verification, external executable or network call occurs. Literal
fixtures cover linear/reversed/circular-origin contributions, safe projection,
clipped/unresolved annotations, readable physical/IUPAC source references,
immutable originals, bad intervals, stale revisions and changed result receipts.
