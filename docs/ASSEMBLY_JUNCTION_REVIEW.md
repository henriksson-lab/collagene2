# Assembly junction review

A reviewed cloning plan can be inspected explicitly against its independently selected target. The panel reconstructs the original source-bound plan in a private worker and presents actual junction DNA, homology and available reading-frame/translation evidence. Opening the panel does not automatically rerun scientific calculations.

Inspection options are explicit, including flank size and genetic code. Displayed peptide frames are sequence translations, not a claim that a gene is expressed or functional. Annotation overlaps, disrupted annotations, ambiguous calls and unsupported junction types remain visible limitations. Retained nicks, chemical state and experimental success remain separate from sequence/frame correspondence.

Changing the plan, active document/event, workspace epoch or inspection options clears current evidence and cancels the private worker. Cancellation and unmount also close it. The host capture callback validates the original plan/target and current active document before cloning source models. A response must match the exact requested context before it is displayed or exported. Full JSON export retains the inspection report; importing old report text cannot restore current verification authority.

This is a desktop implementation. Focused native core/MCP fixtures and native/WASM compiler checks cover the integration; browser acceptance remains deferred until the full in-scope implementation is complete. Mobile work remains deferred.

A [literal stop-codon example](../fixtures/junction-frame-demo.md) independently specifies the same window as `K*F` under code 1 and `KWF` under code 2.
