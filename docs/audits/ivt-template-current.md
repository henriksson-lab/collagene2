# F40 IVT-template foundation — 2026-11-28

The module in `crates/rna-chemistry/src/ivt/` is integrated with binary worker,
MCP and web controls. Native/WASM checks pass, including nine core fixtures,
three binary fixtures, one actual MCP SDK case and three web-model cases.
[Combined execution evidence](verifier-and-design-current-integration.json).
Graphical acceptance remains open. Existing
covalent RNA chemistry remains a separate owner; no verifier representation,
restriction scanner, RNA simulator or physical-material authority is introduced.

The admitted profile checks one declared consensus T7 promoter and original
CleanCap AG N-7113 initiator. [Exact source/rights ledger](../../crates/rna-chemistry/data/ivt/NOTICE.md).

| Requirement | Exact admission | Deliberate boundary |
| --- | --- | --- |
| T7 consensus | Coding-strand `TAATACGACTCACTATA` at −17..−1 relative to declared +1, from the [NEB primary FAQ](https://www.neb.com/en/faqs/do-i-need-to-change-my-promoter-sequence-to-use-cleancap-reagent-ag) | Other promoter variants/polymerases are not evaluated; profile mismatch is not biological inactivity. Missing linear-source flanks and compatible IUPAC ambiguity are indeterminate. A circular transcript overlapping its separate promoter is an explicit geometry conflict. |
| Original AG cap | Exact +1/+2 `AG`, from [TriLink N-7113](https://shop.trilinkbiotech.com/shopping/products/ivt-and-capping-reagents/rna-cap-analogs/cleancap-reagent-ag/) | No cap-occupancy, incorporation or efficiency claim. N-7413/M6/AU are separate products with no admitted rules. |
| Additional start requirement | Caller explicitly chooses none, G, GG or GGG | No implicit GGG default. AG plus a G-start declaration is a conflict even when one individual sequence test passes. [Conrad et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC7429497/) studied particular T7 conditions and distinguishes SP6; those measurements do not calibrate this checker. |
| Runoff end | Supplied coding/template strand-cut offsets; compare template cut with declared end; classify blunt/5′/3′ overhang | [NEB guidance](https://www.neb.com/en-gb/faqs/how-to-prepare-plasmid-templates) supports a 3′-end advisory. No cleavage, enzyme recognition, retained physical promoter, isolation or actual end homogeneity is established. An absent cut declaration stays indeterminate, including linear sequence records. |
| Encoded poly(A) | Declared interval, exact A/other-concrete/ambiguous counts, longest exact A run, all run counts and bounded run preview; terminal-position check | No stability, expression or half-life inference. [Trepotec et al.](https://doi.org/10.1261/rna.069286.118) is an experiment, not a universal host-risk rule; no quantitative data admitted under its NC article licence. |

Coordinates always identify the original DNA record. `Transcript.start` is the
first coding-strand base, not a leftmost coordinate; reverse walks towards
decreasing original indices. Its complement is the polymerase-read template.
Transcript and poly(A) ranges map to at most two original half-open spans in
traversal order. Linearization offsets are bonds measured along that same coding
axis from +1; original bond coordinates are retained separately for both strand
roles. Reverse-complement and every circular-origin fixture use literal independent
coordinate oracles. No source or intended transcript is edited.

The public API exposes installed profiles, immutable profile references, a
source-revision-bound request, assessment/revision and report binding/fresh replay
validation. Findings have named rules, scoped sequence locations and source IDs;
there is no aggregate biological pass. Report hashes establish documentary
integrity only. Scientific consumption uses `validate_for(source, request)` to
recompute; a caller-supplied ID, source URL or recalculated hash cannot admit an
altered built-in profile. Declaration provenance is retained for both transcript
intent and supplied cut geometry.

Resource admission precedes cloning/hashing: one million template bases, 16KiB
request, 32MiB source model, one MiB report, 16 findings and eight assumptions.
The promoter retains at most17 original positions, initiation context at most8.
Poly(A) scanning is linear; at most1,024 runs and64 DNA preview calls per run are
retained. Full counts remain exact when segmentation is incomplete, and the
incompleteness is explicit. Nested transport vectors have deserialize bounds;
required nullable options do not silently disappear from the current format.

Next integration: exact document capture and worker/MCP/UI controls, then reuse
the existing restriction catalogue/material owners to propose and independently
verify downstream cleavage candidates. No new restriction algorithms, external
executables, external scientific parameters or dependencies are needed here.
