# Reviewed SBOL documentary subsets

Selected plasmid publication can explicitly replace an original multi-component
SBOL attachment with its selected structural graph and a documentary receipt.
Ordinary publication never silently takes this alternative. The local workspace
and its exact original bytes remain untouched. Native molecule/revision/event/
capsule identifiers and the original `sbol:source_hash` remain identical.

The two generated attachment kinds are `sbol-selected-subset-v1`
(`application/ld+json`) and `sbol-selected-source-receipt-v1`
(`application/json`). The receipt states the original hash, NEW subset hash,
source profile, versioned subset policy, required component/sequence root pairs,
explicit selected alternatives, omitted object/sequence-edge COUNTS,
`original_bytes_included:false`, and `derivation:"unverified-documentary"`.
It includes no omitted identifiers, names, values or original filenames. Retained
annotations are retained verbatim: this is graph scoping, not free-text redaction.

The canonical profile sorts JSON object keys and graph nodes by IRI, preserves
declared array order, and serializes expanded JSON-LD. It is not RDF URDNA2015.
Component/Sequence/Feature/Location/Constraint structural edges define closure.
Selected additional classes or unknown properties are unsupported. Vocabulary
references that name omitted objects in the original graph are rejected too.
Nested components with alternative sequences require explicit choices; no first
sequence is chosen. Locations or source locations cannot force an unselected
alternative: the canonical result is re-imported with ownership validation.
A shared Sequence retained by another explicitly selected root is permissible.

Native frontend API:

1. `sbol_subset::choice_requirements(&workspace,&selection)` lists reachable
   alternatives for local review. These local candidates are not public receipts.
2. `prepare_documentary_subset(&workspace,&selection,&PublicationChoices)` derives
   native closure and returns `PublicationReview`. `entries()` and `inventory()`
   are read-only. Choices map original hash to component-to-sequence-set.
3. An explicit UI acknowledgement calls `review.approve()`, producing opaque
   `ReviewedPublication` (Clone, no deserialization constructor).
4. `project_document_closure_reviewed` re-derives the exact source/selection/
   choices/counts context and rejects stale acknowledgement. The returned NEW
   `Projection` can use the existing archive worker. `document_closure_reviewed`
   is the synchronous encoded convenience.

The acknowledgement is local UI intent, not remotely supplied authorization.
Normal trusted-host authentication and publication scope checks still apply.
The backend and actual-versus-permitted archive inventory validation independently
check missing-original receipt policy, hash pins, native claim pairs, canonical
minimal graph and complete partition of required claim pairs. Independent public
parts of one original may cover different pairs; overlapping receipts fail as
ambiguous. Extra supplied assets are not made permitted by a receipt. A received
public part can be republished through ordinary closure without recovering the
private original. Exact-original source bytes, when present, retain the strict
original rule; a receipt cannot bypass it.

Without original bytes the receiver CANNOT establish derivation or omission
counts. They are bounded documentary declarations, not a scientific proof,
identity reassignment, executable native operation, material authority or claim
that edited DNA still equals its historical SBOL source. An unknown receipt
policy or unsupported selected graph is rejected. Equal original/subset hashes
use the ordinary exact-original path instead of falsely declaring original
absence. Preparing a new review currently requires the exact original locally;
existing original-absent public parts can be republished unchanged.

Bounds reuse the design-exchange 16 MiB/10,000 object/200,000 value limits,
workspace archive bounds and one-million native publication reference budget.
Each receipt is at most 1 MiB. Tests in design-exchange `subset` and
library-sharing `sbol_subsets` use independent literal graphs, private sentinels,
exact native pins, stale review and tampered receipt cases. The six structural exchange and ten new/legacy native sharing fixture groups
pass in parent validation. The web sharing controls now expose explicit reachable
sequence choices and a separate approval step bound to source/head/workspace,
login session and destination. Integrated UI/worker WASM checking passes;
the authored WASM-only capture fixture has not been executed. No browser or
application HTTP checks were run.
