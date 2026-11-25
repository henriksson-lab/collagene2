# Source-bound Type IIS primer proposals

`golden_gate_primer_design` proposes actual ordered primers for explicit source
intervals, then checks exact PCR and pinned-catalog digestion. This is a preparation
planner, not a target-directed sequence generator or a Golden Gate mixture/yield
prediction. The independent target is included in the context binding and never
chooses a core, endpoint, fusion or primer base. Returned preparations can be passed
to the existing ordered ligation compiler and fresh verifier.

## API

- `Request { name, fragments, enzyme, clamp, spacer }`.
- Each `Fragment { role, input_id, start, end, reverse, forward_length,
  reverse_length, right_overhang }` names an actual original input, not a copied
  fragment posing as a template. Coordinates are zero-based, end-exclusive.
  Circular sources may cross the origin once using an unwrapped end.
- `source_preflight(&[(&str, &Molecule)], target, request)` guards borrowed models
  and reserves aggregate work before source copies, hashing or scientific scans.
- `preview_revision(&[PlanInput], target, request)` binds exact ordered originals,
  canonical revisions, independent target and every request field.
- `design(...) -> Report { context_revision, source_revisions, request,
  fragments, rows, warnings }` returns existing
  `assembly_recipe::Preparation::PcrCatalogDigest` rows, with the full embedded
  catalog edition/content/identity/geometry pin and no end repair.
- `Report::validate_binding(expected)` is a bounded readonly consistency check.
  It does not replay PCR, restore imported physical authority or authorize a
  saved positive conclusion.

`GoldenGateEnzyme` is the existing BsaI/BsmBI/BbsI/SapI enum. Actual geometry is
checked against the selected embedded catalog record. Clamp bases are an explicit
operator choice of 1–12 uppercase A/C/G/T bases, not a manufacturer efficiency
recommendation. Spacer length must equal top cut offset minus motif length:
BsaI/BsmBI/SapI require one base; BbsI requires two. Fusion length is four for
BsaI/BsmBI/BbsI and three for SapI. No kit conditions, incubation, methylation
clearance, star activity or ligation-efficiency values are inferred.

Each returned `Row` carries role/input/revision, named forward/reverse `Primer`,
`TailSegments { clamp, motif, spacer, fusion }` for each oligo's own 5′→3′ order,
the catalog selection, exact digest fragment index, actual oriented retained top
and both physical `End` descriptions. Empty fusion segments mean no synthetic
fusion in that primer. The actual ordered tail is the segment concatenation.

## Strand convention and reverse rows

Let `S` be the selected original top interval, `k` the profile stagger, `F` the
explicit right fusion in the desired downstream assembly-top direction, and `T`
the common clamp + forward recognition motif + spacer.

For a forward row:

- Forward oligo tail is `T`; its annealing domain starts at the source interval.
- Reverse oligo tail is `T + RC(F)`; its annealing domain reverse-complements the
  selected original interval's end.
- The retained digest top is exactly `S`. Its left protruding top is `S[..k]`;
  its right protruding bottom, in that strand's own 5′→3′ order, is `RC(F)`.

For a reverse row, PCR binds only the original interval `S[..len-k]`:

- Forward tail is `T + RC(F)`.
- Reverse tail is `T + RC(S[len-k..])`.
- Before orientation, the digest top is `RC(F) + S[..len-k]`, with the original
  last `k` source bases retained on the opposite-strand flap.
- Reversing this actual duplex gives retained top **`RC(S)`**, left protrusion
  `RC(S)[..k]` and right protruding bottom `RC(F)`.

Thus reverse orientation is not implemented by reverse-complementing a full
selected top and pretending its ends remained unchanged. The planner compares
these independently specified source formulas to fresh exact-PCR sequence and
actual two-strand digest geometry. Both primer domains must fit without overlap;
reverse rows need the two domain lengths to fit within `len(S)-k`.

The left fusion comes from the selected oriented source itself. No neighboring
row need be present in this request: callers may propose only their PCR rows and
preserve other existing-site catalog preparations. Compatibility across all
rows, circular closure and the independent final target remain the subsequent
ligation verifier's responsibility. Linear endpoints still require explicit
right fusion choices; the planner does not invent terminal chemistry or discard
end evidence.

## Conservative applicability and bounds

The current defined profile requires a concrete, physically eligible original
DNA template. Detached reference-only, variant, tracked/physical intermediate and
imported restrictive state claims cannot be promoted to a new substrate. Original
sources and their annotations/raw evidence are never changed. No source sequence
is silently domesticated to remove internal recognition sites.

Before PCR materialization, both oligos are scanned on the actual original:
each must have its single intended exact orientation/site, with no alternative
or swapped binding. Extra binding sites, ambiguous source calls or unsupported
geometry withhold the proposal rather than select a favorable PCR product.
Actual amplicon digestion must have exactly the two eligible added cuts, isolate
the intended middle fragment, and reproduce the independently specified retained
top, stagger and cleaved 5′-phosphate/3′-OH ends. This complete-digest geometry is
not a claim about experimental enzyme performance.

There are at most 32 source inputs/rows, annealing domains are 8–80 bases, aggregate
source plus target DNA is at most 8 MiB, proposed row output is capped at 2 MiB,
and aggregate reserved source-scan work is at most 100 million comparison units.
Shared borrowed model shape limits precede validation/hashing, and serialized
report hashing is capped by the existing 32 MiB recipe budget. Historical report
validation bounds vectors/strings before serialization and checks returned recipe,
named oligos, segmented tails, source revision references, catalog profile,
orientation and end consistency. A host must still compare the live request
revision and freshly reconstruct before treating any preparation as current.

## Source-authored fixtures

`core/tests/golden_gate_primer_design.rs` covers literal forward BsaI/BsmBI/BbsI/
SapI tails and strand ends, reverse BsaI/SapI stagger handling, circular-origin
source intervals, independent wrong targets leaving proposed oligos unchanged,
internal motifs, binding alternatives, stale source revisions, invalid fusions/
spacers/intervals, row capacity and historical-row tampering. These five fixture groups passed, together with four preparation cases, two
binary protocol cases, one real MCP SDK pipeline group and nine existing recipe
regression cases (21 total). The full 121-tool stdio smoke and final web/worker
WASM compiler check also passed. Browser acceptance remains deferred; compiler
and native checks do not establish GUI acceptance.
