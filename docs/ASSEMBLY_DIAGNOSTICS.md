# Local assembly diagnostics from existing evidence

`assembly_diagnostics::overlap(&assembly_recipe::Preview)` and `ligation(&ligation_recipe::Preview)` derive typed row and junction locations from an already bound preview. They do not run PCR, digestion, repair, overlap, ligation or target simulation again, and never parse human finding strings to discover a location. Full original findings remain available unchanged beside concise local recovery guidance.

The adapter first calls the preview's readonly binding validator, checks the actual canonical plan/target/source bindings and supported current engine, and validates product-prefix identity. Each native reconstructed product is appended exactly once at the end of its successful operation, after materialization limits; reconstruction stops at its first error. Thus, after the dense prefix and compiler mapping have been checked, `products.len()` identifies the first operation without a retained product. This is a reconstruction location, not evidence that a laboratory operation failed.

Native shared preparation chains are checked explicitly: original input → PCR, digest, catalog digest or source declaration; optionally PCR → catalog digest; then an explicit repair. Steps must be in compiler order, own exactly one row and terminate at that row's actual prepared source. Ligation end treatments are appended **after all preparation chains**, not interleaved with each row. An earlier PCR may therefore have produced a reference while its later treatment is blocked by another row's failed PCR. `completed_steps` and `total_steps` expose this partial work rather than falsely attributing the failure to the earlier row.

Rows have four structural states:

- `InputReference`: the row directly uses its existing original reference and has no generated preparation steps. This does not establish material state.
- `Reconstructed`: all of the row's declared preparation/treatment steps retained products. This is not a positive chemical or experimental verdict.
- `Failed`: the row owns the first operation that retained no product.
- `Blocked`: one or more owned operations were not reached after reconstruction stopped elsewhere.

Prepared references must equal the actual bound input or retained step product, in the compiler's available leading-row order. Whole-input references can remain available even if the preview's leading prepared-row display stopped earlier. A full set of retained operation products with an invalid assessment identifies the final independent-target comparison, not a failed preparation. Full retained products with `Incomplete` status remain incomplete without inventing a missing step or inferring a local cause from prose. Target agreement never grants authority to a row or molecule.

For overlap recipes, local exact/mismatched/unavailable homology states come only from the bound typed profile. Retained profile strings are checked against their actual oriented prepared-source spans; repeat-search and physical-state limitations remain separate. A typed local mismatch identifies that row pair, but does not prove it is the only reason the molecular operation could not proceed. Generic assembly failure therefore does not arbitrarily select one pair as the cause.

For ligation recipes, junction states say only whether both retained oriented endpoint summaries are available. They do not independently evaluate compatibility or identify one failed pair. Chemical groups, blocked ends, nicks, unsupported physical state and the original full reconstruction findings remain the proper evidence to review. No repair or phosphorylation is inferred from a displayed endpoint.

The returned report preserves the preview context, plan and target revisions, original overall status, row locations, optional first missing-operation/target-comparison location, junction observations and full findings. Messages and recovery text ask the user to review actual declarations and rerun the original method; they never silently change domains, enzyme choices, overlap length, orientation or nick policy.

This adapter is a bounded readonly consistency check, not authentication of arbitrary imported reports or a new source of scientific authority. Hosts must first bind a fresh worker response to the current expected request. Unsupported/old engine semantics, stale input/product bindings, malformed contexts and inconsistent row chains fail closed. The adapter neither upgrades archived evidence nor introduces a new verification engine, scientific DTO version or scientific operation.

Five focused fixtures independently cover middle-PCR failure versus blocked later work, all-prepared assembly failure and typed local mismatch, target disagreement versus complete/incomplete evidence, the actual treatment scheduling and partial progress of ligation recipes, endpoint availability without a guessed cause, and malformed/stale/forged row or receipt bindings. No local browser checks or extra scientific reconstruction are required to render these locations.
