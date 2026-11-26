# Target-free material derivation publication closure

`material_derivation_closure` is an independent child adapter for the selected-document `Graph`. It is registered alongside the corresponding native schema/validator and follows the ordinary public/private selected-closure paths. Documentary dependency retention alone does not grant physical material reuse.

The native `material-derivation-v1` attachment contains the canonical fixed-header record package, with MIME `application/vnd.collagene.material-derivation`. Its content hash covers the whole package. `Index::new` invokes native documentary workspace validation and maps exact product events to their single typed owning record. It also checks package bounds/hash, MIME, exact producer operation (`Target-free native material derivation v1`), source-parent pin, product document/revision/event and duplicate ownership. Run and saved-digest ownership collisions are separately checked independent of attachment iteration order.

No arbitrary metadata value, substring, matching sequence, current document head or source filename establishes ownership. Historical scientific-engine labels remain readable documentary evidence according to native validation; current chemistry is not replayed or granted by the sharing index. Unsupported structural schemas/packages remain explicit errors.

Selecting a derived product retains its exact record and follows the native helpers `source_bindings`, `required_run_ids` and `required_derivation_products`. An original reference adds its exact source event; a saved run/digest/derived origin adds its typed closed dependencies. Nested opaque packages are interpreted only by the native bounded walker, not independently recursed through JSON in the sharing layer. Every dependency enters the ordinary Graph's global edge counter and deduplicated event/run queues. A referenced nested derivation must have a retained exact typed owner; a substituted or missing record is refused.

A reconstruction capsule must equal the retained canonical attachment byte for byte and pass native workspace attestation. Its exact product enters the same event queue. No capsule creates a new ownership claim merely by being present in a plan. A retained event with the exact reserved producer operation but no typed owner is an orphan, not an ordinary source reference.

Index construction alone retains no record assets. Thus unrelated derivations and their private sources/targets stay outside a selected publication. Existing Graph inventory code chooses only required dependency heads and explicitly preserves the selected document's namespace/head/event. It does not move dependency documents to their unrelated private current heads. Existing presentation/draft exclusions remain in effect; this adapter introduces no broad unknown-attachment escape.

## Graph integration

1. The child module and native core schema are registered together. `Graph::new` explicitly recognizes the new attachment kind while retaining rejection of other unknown kinds.
2. At the start of `Graph::resolve`, construct a local `Index::new(self.source)` and call `validate_existing_owners(&self.owners, &self.digest_owners)`. A local index avoids borrowing an index field while mutating the graph.
3. In the retained-event loop, after existing run/digest owner handling, call `follow_product(self, &id)`. If no owner handles the event, call `reject_orphan(event)` before treating it as an ordinary reference.
4. For `ReconstructMaterial`, `follow_capsule` handles its explicitly typed derived capsule. Do not scan candidate bytes or metadata for event-looking strings. Existing typed original-run and saved-digest hooks remain unchanged.
5. `tests/material_derivation_closure.rs` constructs actual target-free PCR→digest→new PCR/physical homology assembly saved as a native run. Public/private projections reload and replay that final product while excluding private current source heads, unrelated documents/candidates, draft-v13 state and the optional historical primer-order table. A selected digest alone retains its own ancestors without retaining a downstream run. Missing owners, wrong attachment hashes, exact-event aliases and a well-formed orphan producer event are refused. Native historical documentary readability/current-engine authority remains the core validator's separate contract.

No source/archive mutation, external process, SQL, filesystem operation, desired-target selection or new scientific authority occurs in this adapter.

Worker91 uses a shared32MiB serialized-origin ledger across Run/Digest/Derived roots and descendants, while legacy CGP packages and individual derivation records remain bounded8MiB. Depth, steps, source bases and model text are counted across the same tree; a nested opaque byte package does not reset them. No scientific status is upgraded by that resource-policy expansion.
