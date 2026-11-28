# Selective Golden Gate junction edits through MCP

`build_edited_assembly_intent` first fixes the expected sequence from immutable
original source regions, their order/orientation and explicit junction trims or
insertions. It does not digest or amplify DNA. Keep its returned target DocumentRef
separate from any later preparation simulator.

`design_golden_gate_junction_primers` takes original source bindings and a typed
request for that same intent, a selected Type IIS enzyme, explicit clamp/spacer,
and ordered retained-versus-designed preparation rows. It returns an advisory
recipe, context/intent/target revisions, retained-row indices, selected new ordered
oligos and a complete evidence artifact. It does not install a target or change
source documents. Retained manual names, primer domains/tails, catalog pins,
fragment choices, orientation and terminal treatment remain exact.

Review the new primer-tail segments in their own 5′→3′ direction, including
recognition site, spacer, fusion and any synthetic junction bases. Then compare
the designer's intent/target revisions with the independently built reference.
Pass the returned **whole** recipe and that target to `propose_ligation_assembly`
for fresh PCR, pinned digestion, strand-end and final-target checking. The ordered
ligation check is not an exhaustive Golden Gate reaction-mixture/yield model.
No experimental efficiency, supplier condition, methylation clearance, repair or
winner is inferred from a reference sequence or advisory primer design.

The new SDK fixture uses two 20-base originals, a circular expected construct and
an explicit `GC` insertion after the first core. It regenerates only the first
row and retains the second row's named manual PCR→BsaI preparation. The expected
42-base target is asserted before primer design; exact retained recipe fields and
new ordered tail are checked before fresh ligation/target verification. Negative
cases cover wrong independent target, stale source references, changed retained
role/order and invalid catalog pins. Repeated requests must reuse immutable
artifacts, and original exports must remain unchanged.

`crates/mcp/tests/methods/golden_gate/golden_gate_junction_tools.rs` contains the SDK fixtures.
Validation results are recorded in [implementation status](IMPLEMENTATION_STATUS.md).
