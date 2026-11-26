# Assembly output finalization

`cloning_outputs::finalize` is a private-workspace transaction. It first asks the existing run recorder to freshly reconstruct the complete plan against its exact independent target and historical source-event bindings. Invalid or mismatching products are refused; an incomplete reconstruction additionally requires explicit `allow_incomplete`. No stored receipt supplies scientific authority.

The returned workspace retains every original source, preparation product, recorded product event and source parent. The final document receives a separate naming event. Its recorded original event remains replayable. If a user has subsequently edited that document beyond its name, finalization opens the recorded product event as a new named descendant instead of replacing the edited head. Repeating the same operation reuses an exact matching named descendant. The returned construct binding always identifies its current document head.

`show_prepared` controls only an attachment of kind `cloning-output-visibility-v1`. The attachment stores version 1 and a unique list of existing hidden document IDs; `load_hidden` checks these references. The default hides earlier products of the selected run, keeping its original final product and selected construct visible. Other runs' visibility is retained. Applications should keep the active document accessible and offer advanced access to the full retained history. Hiding never deletes history or changes molecular metadata.

Optional primer saving uses the shared declared-plan oligo collector, with stable first-use order and deduplication. Only complete, nonempty tables are saved as existing `oligo-table-v1` JSON attachments. No quantities, stocks, modifications or synthesis success are invented. Identical tables reuse their content-addressed attachment.

`finalize_revision` binds the incoming workspace fingerprint, plan, exact source/target bindings and options. `Commit::validate_binding` checks a decoded result's workspace, run attachment, recorded assessment/products, construct content and event ancestry, visibility partition and optional primer table. This is a readonly consistency check, not fresh scientific verification. Hosts must compare the request token and current workspace before replacing live state. Any failure leaves the caller's workspace unchanged.

Names must be trimmed, nonempty, free of control characters and at most 256 bytes. Existing workspace, model, plan and output byte limits remain in force. Output selection does not weaken physical-state rules or experiment uncertainty.

## Existing-target plans and inverse PCR

General candidate review now uses the same actual-product viewer, primer CSV/TSV,
HTML/PDF steps and output chooser as the assembly dock. Inverse-PCR forwards its
exact original source event, independent target event and epoch with the candidate.
The handoff checks the displayed design, literal form, generation and current
source/target bindings before accepting it. It does not discard those preferences
or locate another source merely because its sequence matches.

`web::assembly_review::stage_existing` validates an already specified active target
and ordered input-event bindings, then stages the candidate without creating an
assembly intent or changing document heads. Explicit historical source events
remain usable after their source head advances; native reconstruction reads that
exact original event. The target must still be the exact active head event. An
edit/revert or equal-content target alias is a different event and requires review.
Imported candidates without exact source preferences are not silently given them.

The main review retains the target event with its fresh verifier receipt. Actual
product review, output creation and exports check that context. Shared exports
also recheck their host capture after the worker returns, before downloading.
This catches event changes even when target letters remain equal. The finalizer
freshly reconstructs again and retains all intermediate/source events in one
durable commit, using its existing invalid/incomplete assessment policy.

Typical inverse-PCR flow: declare/open the intended reference → choose the
original source and explicit preparation → review candidate → Verify candidate →
inspect the actual circular product → create named construct, optional primers
and prepared references. Order exports retain the named PCR oligos and tails;
terminal phosphorylation is a reaction treatment, not a supplier modification.
The output path adds no scientific engine or external executable.

A failed final operation still permits exporting its diagnostic HTML/PDF protocol
from a correctly bound fresh receipt. It does not permit opening a preceding
preparation as the final product or enabling construct finalization. Product
availability and receipt identity are separate checks.

## Reusing a saved product as material

`run_material::originating_products` discovers documentary run/step choices for an
exact retained event. It follows an unambiguous chain whose reference content is
unchanged except for its name; equal letters in an unrelated document do not
establish ancestry. An intervening sequence, annotation or metadata edit stops
that chain even if a later edit restores the original content. A discovered origin
is a reconstruction request, not a statement that reusable material exists.

`run_material::reconstruct` resolves the original run, exact original input events
and independent target, then reruns the verifier. The entire recorded product set
must match the current reconstruction. The selected output must have actual
private `TrackedDuplex` state. Reference-only products, withheld graph bridges and
other outputs without established strands remain unavailable as material.

The returned `FreshProduct` has private fields and no serialization/deserialization
implementation. Its actual strands, nicks, termini, positional modifications and
generation records come from that same verifier pass. The complete report remains
available so assumptions, incomplete findings and unselected components are not
erased. Its selected reference retains a checked presentation name. Returning a
summary over MCP or the worker protocol cannot restore this in-process receipt.

Ordinary source selection in the assembly dock is still to be connected. The
durable native handoff below now persists executable origin declarations, retains
source-event bindings and independently verifies the downstream target. Simply
concatenating plans is insufficient: changing a plan
changes intermediate reference revisions, invalidating revision-pinned TOPO and
reference-network selections. Preserve their originating reconstruction boundary.
Never turn a saved duplex summary into an original-input declaration, reapply
source-generation assumptions to descendants, or infer nick repair from a matching
reference sequence. Existing PCR restrictions for staggered/nicked templates still
apply. No external scientific executable is involved.

`run_material::capsule` packages the unchanged origin using the existing typed
candidate encoding. `Operation::ReconstructRunProduct` consumes the selected saved
input reference and freshly evaluates that closed origin against its own target.
Existing downstream PCR, digest, repair and ligation operations consume its actual
tracked step material. Selected origin targets remain acceptance evidence, not
nucleotide parents. Incomplete inner findings remain incomplete in the enclosing
verdict. The complete typed originating report is retained as bounded documentary
`PhysicalProcessingEvidence::RunProduct`, including graph evidence and byproducts.

Iterative preflight checks at most eight origin generations, 1,000 aggregate steps,
8 MiB of embedded candidate bytes and 8 MiB of aggregate input/target bases before
hashing or science. A retained origin report is limited to 32 MiB. Nested candidate
bytes are raw in binary communications and base64 in JSON archives/MCP, avoiding
JSON integer-array expansion. Saved reports have the same compact byte handling.

Candidate verification reconstructs declared science. Native workspace adoption
separately attests every retained origin run, input, target and product event, and
the outer source binding must equal the explicitly selected saved event. An
unrelated equal-content alias is rejected. No serialized graph or tracked summary
is used to supply material. Replayed origin primers belong to already-prepared
material and are excluded from new primer orders; explicit downstream PCR primers
remain orderable. HTML/PDF protocols identify the origin and retain its limitations.
Manual plan review exposes the reconstruction dependency and immutable origin.
The ordinary preparation editor now also provides **Reuse a saved cloning
product…** beneath the row's preparation selector. Find retained origins, choose
one explicitly, and let the worker reconstruct it. The selection changes only
that row's material origin; primer domains, enzyme pins, repair/treatment,
orientation and independent intended target remain explicit settings. A changed
source event clears the pin. Original-only duplex declarations and ssDNA annealing
cannot consume reconstructed duplex material. Removing the origin is explicit.

`Preparation::FromRun` is a flat preparation wrapper around the same origin and
one supported WholeLinear, Digest, CatalogDigest, ExactPcr or PcrCatalogDigest
preparation. Actual ends and modifications feed the existing operation chain.
PCR/digest fragment inspection has an origin-bound worker route, so saved sources
use fresh native reconstruction without deleting physical-marker metadata. Draft
v5 persists the cheap shared documentary pin; opening a draft restores no verdict.
Older ordinary/paired drafts retain their existing v1–v4 layouts. Library-choice
projection of material origins remains deferred until per-choice bindings exist.

MCP clients use `list_run_material_origins` and `prepare_run_material` with an
explicit immutable workspace reference and exact source event. The resulting
artifact contains a recipe, not material authority; normal recipe proposals and
verification reconstruct it again. Neither tool publishes documents or a target.

Validation includes a named 40 bp circle assembled from two PCR fragments, with
phosphorylation explicitly applied to one. Fresh reconstruction retains the top
nick at 20 and bottom nick at 0; it never upgrades the circular reference to a
covalently closed molecule. Circular outer-end fields remain non-terminal
placeholders, not evidence of nick chemistry.

Reusing material does not extend the downstream reaction model. In particular,
PCR of a nicked/staggered template and digest propagation of existing nicks remain
withheld by the current engines. The replayed material and its nicks remain in
the reconstructed prefix; an unsupported downstream step cannot produce or
finalize a silently repaired product.

## Direct virtual digestion of saved constructs

VirtualDigest now offers the same explicit saved-origin picker as preparation
rows. Choose a retained origin, then run the digest. The worker freshly replays
that origin and applies the existing tracked digest; it never reconstructs
material from a saved end summary. The ordinary inspector remains unchanged.

Both the fragment-list and cut-boundary assembly buttons retain the exact saved
source event, original source molecule, origin recipe, catalog/options and actual
two-strand geometry. The draft emits `FromRun::CatalogDigest` and reconstructs
again during method verification. Cut-defined selection uses the source's actual
bottom-strand free endpoints; a staggered saved source need not start at zero or
end at the top reference length. Forward intent can retain such fragments. If a
reverse orientation requires bases outside the saved top reference, the current
draft projection refuses it rather than inventing sequence; an explicit intended
sequence/strand-aware intent extension remains necessary.

The report shows original ends/nicks and source/per-fragment modification
summaries. Complete cleavage is still an explicit conditional assumption, not a
methylation-permission assessment. Known and unknown calls and generation records
remain intact. Existing nicks preserve recognition evidence but withhold physical
products under the current digest engine.

The saved-origin route currently supports inspection, JSON evidence and assembly
handoff. Ordinary named digest setups, gel handoff, GenBank product export and
reference installation remain disabled for that route until they can retain the
same origin scope. They do not silently fall back to an ordinary source. Those
actions remain available for ordinary-reference digests.


## Intentional-mismatch PCR fragment preparation

Fragment preparation can use explicit forward/reverse primer annealing domains,
separate 5′ tails and a recorded mismatch/3′ binding policy. The shared recipe
compiler emits the existing mismatch-PCR operation; overlap assembly, ordered
ligation and TA/GC preparation retain the same source event and policy. Saved-run
sources are freshly reconstructed before amplification. The complete-amplicon
reaction option retains this policy rather than converting it to exact PCR.

The intended reference remains independent. Primer changes never rewrite the
selected target to make a result pass. The standalone method inspector can use
an independently edited/imported target. The ordinary dock declares its reference from chosen source intervals, explicit
junction edits and independently entered internal substitutions. An internal
primer mutation absent from that declaration must fail target comparison.

Unique concrete binding geometry, search completeness and strand compatibility
are checked by the existing verifier. Truncated searches cannot establish a
complete result. Reconstructed PCR products have newly synthesized unphosphorylated
5′ ends; ordered ligation still requires an explicit terminal treatment where
necessary. Binding compatibility does not establish amplification efficiency.


## Internal substitutions in ordinary assembly intent

Open Internal substitutions in the assembly dock. Select one source base in the
sequence pane, or enter a one-based position within the active contribution in
construct orientation, then choose its replacement. Each declaration stores its
stable fragment role, zero-based original contribution offset, expected base and
replacement. Forward/reverse and origin-crossing contributions use the original
source coordinates; edits remain with their source when fragments are reordered.
Replacing a source, reversing a contribution or swapping oligo roles clears the
affected contribution's edits. Unrelated substitutions are retained.

The native builder first applies explicit junction trims/insertions, then maps
the substitutions into retained original cores. A trimmed-away, duplicate,
unknown-role, ambiguous or stale-expected-base edit fails. Original documents,
primers and physical preparations remain untouched. Affected projected features
need review and lose stale translations; full original annotations remain in
source documents. Unsupported affected coordinate-bearing qualifiers fail with
an explicit diagnostic. This interface supports up to256 single concrete-base
substitutions; general indels and per-library-row internal edit intent remain
separate required integrations.

Review intended construct freezes the complete declaration. Finalization rebuilds
it from the candidate's original inputs and checks all exact workspace source
events and the complete target revision. Editing primers cannot rewrite that
target. Exact-PCR primer proposals currently exclude internal substitutions;
enter intentional-mismatch primers explicitly while retaining literal junction
editing. Draftv6 persists declarations, while legacyv1–5 drafts keep their prior
meaning. No saved draft or reference report restores scientific verification.
