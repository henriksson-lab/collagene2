# Source-bound internal-site domestication

Open **Cloning → Remove internal sites**, or use the sequence-selection action.
The active DNA document is the original source. Choose catalog enzymes and a
source scope, select coding features permitted for synonymous search, and declare
protected bases or recognition sites to retain. Coordinates are 1-based inclusive
in the form; circular ranges may cross the origin. Inspect the actual sites and
feature diagnostics before proposing changes.

## Reference constraints and explicit review

The native adapter pins the complete source revision, exact feature definitions,
catalog geometry and every option. Hosts additionally retain the actual document
event; equal DNA in a different event is not interchangeable ancestry.

Scoped unwanted motifs, retained recognition occurrences and protected bases are
distinct constraints. Search may change complete codons intersecting an unwanted
site, including codon bases outside the motif when necessary. It does not silently
remove wanted sites elsewhere in the same CDS. The shared constrained codon
engine applies positional domains and locks while searching, with changed-base
count first, declared preference second and deterministic lexical ties. Global
codon weights cannot substitute for this positional objective.

Choose a proposal explicitly and run **Check substitutions**. The independent
checker starts from original source letters and literal single-base patches,
rescans selected enzymes across the whole physical DNA on both strands and at the
circular origin, and checks every affected coding feature. Phase, biological
join/complement order, genetic code, initiation and terminal-stop semantics remain
those of the original features. Partial codon bases stay protected. New selected-
enzyme sites and changes to retained/unscoped recognition sites are rejected.

Feature locations and qualifiers retain their exact original coordinates; touched
features are marked for review. A spanning replacement is not used to remap
annotations between distant substitutions. Unsupported translation exceptions
remain explicit. Fuzzy or partly unresolved CDS coverage cannot be assumed absent
at another edited position; uncertain coverage blocks unsupported claims.

The checked reference is frozen before molecular reconstruction. It can be opened
as a derived reference with its original event parent and full intent evidence.
This does not prepare DNA or restore positive reaction verification. Opening a
different equal-content source produces a separately parented reference.

## Physical preparation and outputs

Physical preparation is a separate action. One route is original-template
inverse PCR with an explicit replacement interval, forward/reverse annealing
lengths and tail split. The interval may wrap around a circular source. Every
selected substitution must lie inside it and must be encoded in the actual ordered
primers. Unchanged intervening bases are copied literally. Both complete ordered
primers must fit the 200-base limit; retained annealing domains must be nonoverlapping.

The shared current inverse-PCR compiler builds actual primers and the PCR →
phosphorylation → circular ligation plan. Native materialization verifies this
canonical plan once against the already frozen reference.

The overlap-extension route prepares a whole original linear or circular source
using two actual branch PCRs, reciprocal extension and outer PCR. Circular sources
then receive explicit terminal phosphorylation and intended circular ligation.
Native anchored primer ranking receives only original DNA and explicit literal
patches. The independently checked target is frozen before planning and never
selects a pair or product. All six primer uses, native call/filter/search evidence,
and the one canonical physical plan remain inspectable. No recoded reference is
invented as pre-existing template material.

Both unchanged flanks must support the requested nonoverlapping domains; long
bridges must fit actual oligos. Edits near an original linearization boundary,
unsupported distant edits or saved-material templates require another explicit
route or reviewed rounds. This slice does not silently rotate the original,
recover source chemistry from metadata or add external primer tails.

Prepared.proposal_ready describes only satisfaction of requested primer-route
constraints. A ranked diagnostic fallback is visibly nonready even if some final
sequence could match. The actual VerificationReport is retained unchanged;
readiness does not assert reaction acceptance, globally optimal primers,
experimental isolation or absence of bounded terminal alternatives.

Review actual oligos, findings and the final reconstructed product before opening
the plan. An invalid preparation does not acquire authority because its reference
constraints passed. Plan review retains the exact original input and independent
target events. Existing fresh verification, product/junction inspection, output
selection, primer CSV/TSV and HTML/PDF instructions use that same typed plan.
Prepared-product annotations retain the existing conservative source projection
and review flags. The checked reference's annotations are not substituted for the
actual reconstructed product or treated as physical source ancestry.

## Bounds, interfaces and remaining work

The native source adapter permits up to 256 CDS assessments, 4096 recognition sites
and 256 explicit substitutions, with source/mapping/serialization and aggregate
work limits. Search and independent checker reservations are both charged; pruning
and exhausted work stay visible. Conservative per-feature search locks overlapping
or unselected CDSs. Its failed search does not prove that no joint synonymous
solution exists, and its candidate does not prove a global biological optimum.
Explicit literal candidates can still be independently checked against all
supported affected frames.

The current worker exposes inspection, proposal, checked-reference and preparation jobs.
The MCP surface includes `inspect_domestication`, `propose_domestication`,
`check_domestication` and `prepare_domestication`. Inspection exposes exact feature
pins for subsequent requests. Prepared MCP candidates retain the independent
target, plan and complete evidence atomically, then use existing `verify_plan`,
`export_protocol`, `export_protocol_pdf` and `export_plan_oligos`. Imported receipts
are descriptive data; only fresh checking creates the private native checked intent.

Joint overlapping-frame optimization, broader initiation/translation profiles,
more complete compound-gap search contexts, multi-round/oligo preparation
and curated empirical method guidance remain active delivery work. Actual GUI and
performance acceptance remains deferred until all in-scope code is implemented.
No browser, external scientific program, kit outcome or template-clearance claim
is introduced by this workflow.
