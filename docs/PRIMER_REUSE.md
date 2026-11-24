# Reusing saved primer sequences

`core::primer_reuse` is the shared sequence-reuse boundary for primer panels,
overlap preparation rows and MCP. It accepts an actual DNA `Molecule` and an
explicit original `feature_index`; only `primer_bind` features are eligible.
Bare DNA and a document named “oligo” do not define an annealing domain or tail.
The explicit feature emitted by `oligo_formats` is supported directly.

APIs:

- `preflight(&Molecule)` applies borrowed aggregate native-model bounds before
  cloning, hashing or extraction.
- `resolve(&Molecule, feature_index)` returns `Resolved` with the actual canonical
  source revision, feature index, named `Primer`, sequence origin and warnings.
- `resolve_at_revision(&Molecule, feature_index, expected_revision)` additionally
  rejects a stale actual-document revision before resolving the feature.
- `list(&Molecule, limit)` returns a bounded `List` with one actual source revision,
  total/omitted counts and `Summary` entries. Original indices are preserved;
  unusable entries carry a diagnostic instead of disappearing. Maximum limit is
  1,000. Summary fields include name, reusable flag, annealing/tail/ordered lengths
  when available and the source feature's review flag. Listing hashes the source
  once, never once per feature.

## Sequence policy

Canonical fields are `annealing_sequence`, `5prime_tail` and `oligo_sequence`.
Names use `label` or explicit `primer_name`, otherwise the existing feature label.
Case-insensitive duplicate relevant keys are rejected, including equal duplicate
values; conflicting explicit names are rejected. Missing qualifier values are
errors. Present DNA strings normalize through the existing DNA/IUPAC parser.
An absent tail means empty, never an invented prefix. A present empty tail is
an explicit empty declaration.

A stored annealing domain takes precedence over reference extraction. If it is
missing and both ordered sequence and a declared tail exist, the domain is the
nonempty remainder after that exact prefix. Otherwise missing domains require
exact current coordinate extraction. A supplied ordered sequence must always
equal tail plus annealing; disagreement or an undeclared tail split is rejected.
The complete ordered primer is limited to 200 bases. IUPAC ambiguity is retained
with a notice, without claiming a concrete binding site or synthesized mixture.
Exact-PCR consumers retain their own stricter ambiguity policy.

Geometry fallback requires exact intervals on one strand, contiguous in reference
order, possibly crossing one explicitly circular origin. Outer complements
preserve reverse-primer orientation. Gaps, overlapping/repeated calls, fuzzy or
unresolved coordinates, `Order` locations, mixed strands and locations longer
than one source turn or 200 bases are refused. A feature marked for review cannot
supply missing sequence fields through coordinates, since it may be clipped or
otherwise biologically uncertain.

Explicit stored sequences remain reusable when their recorded location differs
from current template calls, since mismatch primers can be intentional. The
resolver reports that disagreement instead of silently replacing the primer
with reference sequence. Unknown/reviewed locations also generate a notice when
stored fields are complete; they do not become current binding evidence.
A historical `source_revision` qualifier is checked for unambiguous canonical
syntax but never treated as the actual current document revision or a fresh
binding receipt. Hosts pin document/event identity independently.

## Scope and bounds

The operation neither mutates the source nor creates physical oligo material.
Native metadata, raw attachments, old binding reports, modifiers and primer names
do not establish terminal chemistry, annealing, specificity or amplification.
Fresh downstream analyses remain necessary; a decoded `Resolved` object is a
readonly sequence description, not edit or verification authority.

Borrowed native-model bounds guard all metadata, features and location trees
before serialization. Source revision hashing streams canonical JSON with a
32 MiB ceiling. Relevant raw qualifier values are at most 8,192 bytes, primer
names at most 256 bytes, and normalized ordered sequences at most 200 bases.
Coordinate extraction is bounded before allocating a sequence. No external
software, network request or full-reference search is performed.

Literal fixtures cover forward/reverse/origin geometry, independent stored
mismatch and tail fields, duplicate/conflicting qualifiers, IUPAC preservation,
real imported oligo features, missing bare-DNA domains, actual list indices and
stale document pins.

A [saved-primer example](../fixtures/primer-reuse-demo.md) uses explicit domain/tail columns and independently specifies the 22-base PCR reference preview.
