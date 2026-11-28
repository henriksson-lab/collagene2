# Selective junction primer proposals through MCP

Build an independent reference first with `build_edited_assembly_intent`: original
immutable source bindings, declared region order/orientation and explicit junction
trims/inserts determine the expected DNA. Then call
`design_selective_assembly_primers` with those same originals and a request containing
`base`, `edits`, ordered `rows`, per-junction `overlaps`, `primer_options` and
`profile_options`.

Each ordered row is either `DesignPcr { role, input_id, reverse }` or
`Retain(FragmentSpec)`. A retained row is an exact existing preparation: its names,
annealing domains, tails, catalog pin, fragment index, end repair and orientation
must not be silently replaced. The returned advisory recipe can therefore combine
selected new primer pairs, reviewed manual oligos and pinned existing-site digests.
Only selected rows appear in the new ordered-oligo proposal list; retained indices
and row-aligned optional design reports expose what was preserved.

Compare the returned `intent_context_revision` and `target_revision` to the
previously built intent. The design tool does not install a target or prove a
successful assembly. Pass the entire returned recipe and that independent target
to `propose_overlap_assembly`, then inspect its fresh verifier and saved plan.
Individual primer rankings do not establish exact amplification, cleavage,
compatible chemistry, experimental yield or a unique biological outcome.

The new `crates/mcp/tests/assembly/selective_junction_tools.rs` fixture declares a 68-base
three-source target before designing anything. It replans only row A after a `GC`
junction insertion, preserves named manual PCR row B exactly, and preserves an
EcoRV catalog-pinned row C with two actual flanking cuts. Junction homologies are
15 and 18 bases. The SDK source covers the complete preparation/verifier flow,
wrong-target verification, immutable original exports, artifact deduplication,
stale document references, row role/order mismatch and invalid retained catalog
pins. These are authored fixtures; no executed result is claimed here.
