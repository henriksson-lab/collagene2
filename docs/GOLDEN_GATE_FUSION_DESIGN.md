# Golden Gate fusion placement and empirical counts

This workflow extends the existing junction editor and original-template PCR →
pinned Type IIS digest preparation. It does not introduce another cloning wizard
or replace the whole-reaction candidate verifier.

## Independent DNA and physical preparation

The source intervals, orientation and existing trim/insert declarations define the
intended construct before placement search. A signed junction window asks where
to place a fusion without changing those intended letters. Bases transferred
between contributions are recorded as primer-added prefixes or suffixes, not
invented original-template material. Original source slots and intended edits stay
unchanged when actual oligos are applied.

Search and materialization are separate. The search checks declared palindrome,
reverse-complement uniqueness and Hamming-distance constraints. It reports bounded
work and completeness; a capped search cannot prove impossibility or optimality.
The user explicitly chooses a candidate. Materialization reconstructs real PCR
and catalog-pinned cuts from the original inputs, rejecting ambiguous preparation
or added/internal restriction sites through the existing scientific engines.
Manual preparations remain fixed. A candidate satisfying cheap end constraints
may still fail preparation and remains a proposal until that check succeeds.

Native requests permit up to 32 contributions, signed shifts of at most 128 bases,
100,000 assignments and 256 retained candidates, with an additional 20 million
comparison-work budget. The desktop draft retains its existing 12-bin limit.
Unspecified junction windows are fixed at zero. Fixed manual neighbors safely
restrict their shared placement to zero; contradictory windows can therefore be
exhaustively empty without executing PCR. Other domains visit zero and nearby
shifts first. These are bounded search rules, not a proof of a globally best
biological design.

The junction editor applies actual named oligos, tails, catalog selection and
fragment index through its existing guarded transaction. The frozen independent
intent remains separate from the physical recipe. Subsequent reaction discovery,
explicit product selection and fresh target verification remain necessary.
Circular origin handling is directed; visual rotation cannot substitute a
different target or reverse-complement policy. Synthetic copied bases do not
automatically acquire the original feature annotations or chemical provenance.

## Explicit empirical data

`ligation_fidelity` accepts supplied 3- or 4-base directed count matrices. No
experimental data is implied by this supplied-data API. Separately attributed
published profiles are now available through [their explicit API](PUBLISHED_LIGATION_PROFILES.md).
No third-party optimizer code is bundled. Each supplied dataset binds
its conditions, source citation, attribution, license, documentary source hash,
normalization description and canonical content revision. A documentary hash is
not proof that Collagene fetched or authenticated the cited observations.

The independently written `directed-count-pooled-physical-pairs-1` model pools the
two observed rows for each nonpalindromic reverse-complement pair. For end `x`,
reverse complement `r`, and all declared physical end sequences `U`, its conditional
estimate is:

```
(C[x,r] + C[r,x]) / (sum(C[x,y], y in U) + sum(C[r,y], y in U))
```

Every required directed cell, including self-pair cells, must be present. Missing
observations are unknown, rather than zero. Zero totals and palindromic ends also
withhold an estimate. Repeated declarations retain their request positions but
contribute once per physical reverse-complement pair. The reported product uses
each such pair once; a log product preserves comparisons when a positive floating
product underflows. Per-orientation counts and missing partners remain visible.
Aggregate counts use exact 128-bit arithmetic; JSON receipts encode those totals
as decimal strings, while Postcard retains binary 128-bit values. This avoids
rounding large observation totals through JSON numbers.

This is a declared conditional count model, not a claim of equivalence to a
commercial optimizer, a calibrated probability of a successful assembly, or
experimental yield. It does not model abundance, reaction kinetics, partial
cleavage or recutting. Exact end compatibility, empirical ranking and verified
agreement with the intended construct remain distinct evidence.

## Interfaces and remaining scope

Binary worker schema 60 appends fusion search, selected materialization and count
scoring jobs. Historical envelopes cannot carry these new selectors. MCP exposes
`search_golden_gate_fusions`, `materialize_golden_gate_fusions` and
`score_ligation_fidelity` as advisory tools; full bounded evidence is returned in
artifacts. Original immutable source references are retained with that evidence.
No source is edited and no product is released by these proposal tools.

The desktop editor provides signed windows, explicit candidate/materialization
review and optional supplied-matrix scoring. Source, draft, form and workspace
changes invalidate pending proposals. No browser acceptance or scrolling/performance
claim follows from native fixtures or compilation.

Published-profile ranking and an initial original-template domestication route
now exist; see [published profiles](PUBLISHED_LIGATION_PROFILES.md) and
[domestication](DOMESTICATION.md). The UI uses the separate competition report
when linear outer ends are present, avoiding invented complementary competitors
or extra intended joins. The original supplied-data DTO remains readable.

Complete F13 still includes joint biological constraints, broader physical
domestication routes, calibrated experimental profiles and richer cycling-mixture
behavior. These increments do not redefine that completion boundary.
[Research and licensing decisions](GOLDEN_GATE_OPTIMIZATION_DESIGN.md).
