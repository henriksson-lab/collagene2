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
