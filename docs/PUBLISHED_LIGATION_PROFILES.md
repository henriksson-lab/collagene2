# Published ligation profiles and advisory fusion ranking

Five versioned Pryor et al. (2020) supplementary matrices are bundled for
BsaI-HFv2, BsmBI-v2, Esp3I, BbsI-HF and SapI. Software remains independently
implemented MIT code. The original XLSX data are separately attributed CC-BY-4.0;
see [data attribution](../data/pryor2020/ATTRIBUTION.md) and
[the primary publication](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0238592).
No third-party optimization code is incorporated.

`core::ligation_profiles` pins source bytes, conversion/model version, exact
formulation, experimental conditions and data revision. Bounded ZIP/XML parsing
checks literal row/column labels, complete cells, symmetry and unsupported parts;
rows must not be interpreted as lexicographic sequence order. Explicit zero
weights remain zero. These are published symmetric pair frequencies, not newly
reconstructed independent directional read counts. The existing supplied directed
count model retains its original API and binary layout.

Catalog listing and source/scope revision checks do not parse XLSX on the GUI
thread. Worker scoring loads and prepares a pinned matrix once. The UI requires
an explicit experiment choice and shows its formulation, conditions, provenance
and separate data license. Matching enzyme geometry alone does not establish that
the physical preparation used the same formulation or reaction protocol.

## Intended joins and competing physical ends

An intended junction contributes its two actual reverse-complementary 5′ ends.
For a linear assembly, the outer left end and the reverse complement of the
stored assembly-top right terminal key additionally compete as literal singleton
physical ends. Their absent complements are never invented. Outer ends contribute
to denominators, never to the product of intended successful-join estimates.
A one-fragment linear design with no intended joins has no combined estimate.
The model uses distinct physical partner sequences, not concentration-weighted
stoichiometry; repeated declarations do not supply abundance.

`fusion::end_scope` centralizes that orientation rule. Both bundled profiles and
the new supplied-data `CompetitionReport` use the same numerical kernel. Existing
`ligation_fidelity::Report` and its original endpoint keep their historical
meaning. Missing cells or unsupported intended palindromes remain unknown;
unknown and measured zero remain distinct. Combined values assume independent
intended-pair outcomes and do not predict measured reaction yield.

## Review and application

`golden_gate_fusion_rank` searches the bounded mathematical placement space once,
then scores retained candidates with one prepared profile. Positive scores sort
by descending log score, followed by measured zero and unknown; deterministic
candidate identity breaks ties. Search truncation and scientific assumptions stay
visible. Compact ranking rows retain exact integer weights; detailed scoring
shows physical end orientations and denominators.

Ranking never changes intended DNA, picks a physical reaction winner, applies
primers or establishes preparation feasibility. The user selects a placement,
materializes original-source PCR/digest preparations, and explicitly applies it.
The existing full-recipe reconstruction against independent intended DNA follows.
Source/form/profile changes invalidate pending work and displayed evidence.

Worker62 adds jobs143–146/results135–138 without changing candidate12/plan10 or
science15. MCP adds catalog, published score and ranking tools. The supplied-data
score tool accepts optional `physical_competitors`; omission preserves its old
report model. HTML/PDF cloning instructions continue to describe the actual
reviewed physical plan, not an automatically substituted ranking recommendation.

Remaining F13 work includes joint biological constraints and optimization,
concentration-aware or richer cycling-mixture models, additional calibrated
experimental profiles and remaining method/delivery acceptance. This checkpoint
is neither a commercial optimizer equivalence claim nor GUI acceptance.
