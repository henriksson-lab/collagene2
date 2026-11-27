# Protein sequence liability screening

`collagene-protein-liabilities` is a new MIT Rust/WASM core for bounded,
source-bound sequence screening, with worker, Biology-panel and native MCP paths.
Automatic antibody domain detection, residue numbering
and CDR/framework derivation remain required F34 work, not inferred by this screen.

## Current profile

Engine `protein-liabilities-1` defines the following literal candidate patterns:

| Rule | Match |
| --- | --- |
| Canonical Asn deamidation motif | `N[GSNTH]` |
| Canonical Asp isomerization motif | `D[GSDTH]` |
| N-glycosylation sequon | `N[^P][ST]` |
| Methionine oxidation review | `M` |
| Tryptophan oxidation review | `W` |
| Supplied-CDR Asn review | `N`, anchored inside an explicitly provided CDR |

The dipeptide sets follow the definitions in
[Lu et al., 2018](https://pmc.ncbi.nlm.nih.gov/articles/PMC6343770/).
NH belongs to that paper's canonical deamidation set; ND does not. The study also
observed modifications outside canonical motifs and structural dependence, so
matches cannot establish chemical modification or its absence. The profile does
not silently adopt conflicting statements in the product-reference article.

Sequon and oxidation screening are consistent with the scientific discussion in
[the LAP paper](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1011881).
Broad supplied-CDR Asn review corresponds to the separate review behavior described
by [Benchling](https://help.benchling.com/hc/en-us/articles/9684246077325-Annotate-antibody-specific-information-on-protein-sequences).
No source code, antibody datasets, trained models or third-party sequences are
bundled. The algorithm and fixtures are independently authored.

## Interpretation

Results are sequence-pattern matches, not measured PTMs, folding, accessibility,
activity, rates, reaction yield or clinical predictions. They do not modify a
sequence, alter protein mass, prove an antibody identity or install cloning material.
A zero count means no match to this particular profile, not biological safety or
absence of modifications. There are no probabilities or high/medium/low risk scores.

The source must be explicitly protein. B expands to D/N, Z to E/Q, J to I/L and
X to the 22 supported amino acids (20 canonical plus U/O). `CertainPattern` means
every assignment matches the pattern; `PossiblePattern` means only some do. These
labels describe matching certainty, not chemical state or resolved residue identity.
Stop symbols block motif windows; original source coordinates remain unchanged.

Optional provided CDR intervals are zero-based, half-open and bound to the exact
full molecule revision, including annotations. Each has a domain ID and CDR1/2/3
label under a descriptive definition/provenance statement. Overlaps, repeated
(domain, CDR) pairs, invalid bounds and stop-crossing intervals refuse. These are
user-supplied regions, not validated Kabat/Chothia/IMGT numbering. Inside/outside
labels refer only to the supplied regions. A motif anchored inside a provided CDR
may use adjacent source residues beyond that region's end. Missing/empty context
marks the conditional CDR rule unassessed; nonempty context does not prove complete
CDR coverage or absence of other domains.

## Binding, limits and evidence

Full source, options, rule masks, citations and ambiguity/boundary interpretation
bind to the report context. Source hashing uses bounded streaming serialization
and matches the shared molecule revision. `Report::validate_for` independently
reconstructs the complete result from its source/options, including counts and
retained hit order; a saved hash alone grants no authority.

Limits are one million source positions, 256 supplied CDR intervals, 100,000 retained
hits and 512 KiB serialized options, plus the shared model/text/location budgets
and 32 MiB streamed source bound. `max_hits=0` allows counts only. Every position
is scanned even after the hit-retention cap; exact per-rule certain/possible counts
and `hits_truncated` distinguish that cap from a negative screen. Cancellation
returns an explicit error, never an apparently complete partial report.

Seven native fixture groups cover literal published motifs and exclusions, authored
sequons, B/Z/J/X/U/O handling, supplied-CDR boundaries and stale context, retention
caps, cancellation, resource bounds, report tampering, response binding and
independently supplied match validation. An independent concrete-
string oracle checks unions/intersections over all 258 words of lengths 1–3 from
a six-symbol ambiguous alphabet. Native and optional-schema WASM checks pass.
These checks do not establish numbering, integrated application behavior or
experimental/platform acceptance.

## Application controls

With an explicit protein open, expand **Biology → Sequence liability screening**.
Choose the retained-hit limit (zero for counts only) and run the screen. Counts
always cover the complete source, while retained hits can be filtered by rule or
possible matching and are paged at 50 rows. Selection navigates to exact source
coordinates. **Annotate candidate** independently rechecks the match against the
current source and adds a review-required `misc_feature`; it never declares a
measured chemical modification or changes protein mass calculations.

Optional supplied-CDR controls accept domain labels, CDR1/2/3 labels, inclusive
one-based start/end coordinates, definition text and provenance. Current selection
can fill the range. A context is pinned to its exact document history event;
changing source, even to an equal-content document, requires explicit reset/clear.
Native validation rejects overlapping/repeated, stop-crossing and invalid regions.
Source/options/workspace changes cancel or invalidate earlier worker results.

The full report download includes native options and any provided regions. Form
settings and reports are not saved as workspace analysis records; explicit candidate
annotations persist with ordinary source history. No automatic numbering or full
CDR/domain coverage is inferred from supplied intervals or zero motif counts.

## Transport and MCP

Current worker 117 adds job 194 `ScreenProteinLiabilities` and result 184
`ProteinLiabilities`, under the ordinary binary message limits. Stateless and
retained-session dispatch bind the exact native request context. Response decoding
checks catalog/source/options bindings, counts, retention/truncation, hit ordering,
pattern certainty, supplied-CDR membership and consistent observed overlapping
residues. Without the original source these are transport checks; accepting an
imported report as reconstructed evidence still requires `validate_for` with its
independently retained source/options.

The native MCP registry exposes `list_protein_liability_profiles` and
`screen_protein_liabilities`. Screening requires an immutable protein document
reference and typed native options. Source/options are bounded before cloning;
responses preview at most 100 retained hits and store the full JSON report in a
session artifact. Counts-only results remain explicit. These tools do not edit a
source, save a declaration, install a candidate/material or expand thin WebMCP.

The integration batch passes 19 distinct native cases: seven scientific/binding
fixtures, two binary-worker/annotation fixtures, one real MCP SDK round trip,
two web-model source/form guards and seven affected RNA/structure/similarity
worker cases. The final integrated web/worker WASM check and repository formatting
check also pass. No browser/application HTTP checks were run. These results do
not establish graphical/performance, installed-platform or antibody numbering
acceptance.
