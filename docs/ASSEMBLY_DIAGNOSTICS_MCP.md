# Typed assembly diagnostics in existing MCP proposals

`propose_overlap_assembly` and `propose_ligation_assembly` now return an additive
`diagnostics` field. They keep their existing arguments, tool count, candidate
status, saved plan identity and export/verification routes. Their full JSON
artifacts retain all existing top-level preview fields and add the same
`diagnostics` object; clients reading `compiled`, `prepared` or `verification`
continue to find those fields at the original paths.

Diagnostics are computed from the checked preview and its typed compiled rows,
step ownership and dense execution prefix **before** any plan/artifact publication.
The module validates the source preview; the host also compares its context to
the exact expected preview revision. It does not interpret finding text or
reconstruct chemistry from names/metadata. Existing atomic capacity accounting
covers the expanded artifact and candidate together.

The report contains exact context/plan/target revisions, original verifier status,
rows, optional failure, junction diagnostics, original findings and advisory
warnings. Each row identifies its index, role, input, orientation, preparation/
repair steps and completed/total step counts. Its state distinguishes input
reference, reconstructed preparation, the actual failed row and later blocked
rows. A terminal-treatment step belongs to its row; row numbers and operation
indices are not interchangeable.

If every step produced its result but the independent target disagrees, the
failure is a target comparison. Successfully prepared rows remain reconstructed;
they are not relabeled failed PCR. Incomplete status remains incomplete even if
no specific failed operation is identified. Overlap evidence distinguishes exact,
mismatching and unavailable homology. Ligation endpoint availability is reported
without independently declaring compatible/sealed ends or experimental success.
Full verifier findings and physical evidence remain authoritative within their
explicit computational scope.

`crates/mcp/tests/assembly_diagnostics_tools.rs` authors two actual SDK fixture
groups for both existing tools: a failed middle PCR maps to the correct operation
index. Ligation emits all preparation steps before end treatments, so the first
row can have its PCR complete while its treatment remains blocked; later rows
remain blocked, and
wrong targets preserve completed preparations/junction evidence. It also checks
response/artifact diagnostic equality, existing artifact field paths, immutable
candidate/artifact deduplication and unchanged dynamic tool count. Validation results are recorded in [implementation status](IMPLEMENTATION_STATUS.md).
