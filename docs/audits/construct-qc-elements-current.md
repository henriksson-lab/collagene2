# Declared element chains and required reference intervals

Source batch reviewed 28 November 2026. No compiler, tests or browser checks were
run by this owner. This extends the earlier construct-QC foundation; root owns
current plan, transport and presentation updates.

The single current construct-QC options version is 2 and its engine identity is
`construct-qc-2`. There are no fallback decoders or old-version adapters.

## Structural predicates

An element chain explicitly declares a CDS, an upstream promoter role, a downstream
terminator/polyadenylation role and optional distance limits. Roles are assumptions,
not discoveries. Both annotations must have current source/feature-pinned role
declarations, need no review and define supported exact simple intervals. A
captured concrete reference must match the oriented element sequence exactly.
Only then are same-strand placement and directed spacing checked. A mismatching
reference makes the role premise indeterminate; it does not establish absence of
a promoter or termination element.

Linear checks use the declared CDS direction. Circular checks use directed distances
and require the complete flanking chain to fit into one turn without reusing
overlapping flanks. A user distance bound is an explicit structural requirement,
never a biological activity threshold. Compound/fuzzy paths remain indeterminate.
Missing chain declarations produce indeterminate findings for the annotated CDS.

Required-part declarations name a literal reference, an assumed role, orientation
and copy bounds. The existing canonical sequence-search owner examines both
strands and circular origins. The result counts distinct literal intervals,
including overlapping occurrences; it does not establish functional copy number.
A complete zero-hit search can establish absence of the **specified literal
reference**, never absence of the biological role or of divergent homologues.

Ambiguous source matches retain certain/possible count bounds. Ambiguous reference
calls are not admitted to a cardinality proof. Palindromic references count once
per interval, while their direction remains unresolved. Any ambiguous relevant
match, unresolved required direction or incomplete search yields indeterminate.

## Bounds and binding

- At most 64 chains and 64 required-part declarations. Identity comparisons share
  the existing one-million extracted-base budget with CDS/reference claims.
- Captured reference identities and concreteness are computed once per job.
  Reusing a long reference in many short claims cannot repeatedly hash/scan it.
- Canonical search admits references up to 4,096 bases. One reference search is
  reused for multiple roles. The job shares at most 50 million canonical search
  work units, 1,024 retained matches and 256 matches per reference. Defaults are
  20 million / 512 / 128. No search is started after a shared budget is exhausted.
- Report search rows retain coordinates once, with index references from role
  assessments. Finding selections copy at most 16 ranges each. Existing source,
  annotation, location-tree and two-MiB report bounds remain in force.
- Search work is the canonical engine's accounting, not elapsed CPU cycles.
  Its ordinary validation/source-identity overhead is independently bounded by
  64 searches against the admitted source envelope.
- Context identity includes the current canonical search engine as well as the
  exact source, annotations, references, genetic codes and declared options.
  Fresh replay reconstructs every new count, distance and verdict. Restored
  checksums are historical content integrity, not scientific authority.

Ten source-only tests cover both strands, zero distances, spacing/orientation
failures, missing/unreviewed/compound premises, circular order, overlapping and
origin-spanning matches, certain/possible ambiguity, palindromes, shared query
reuse, work/match limits, pin/context tampering and aggregate extraction admission.
All ten pass alongside the earlier nine QC fixtures, current MCP calls and final
native/WASM checks. [Integration evidence](reaction-elements-current-integration.json).

No element-function inference, initiation-strength model, vector-class threshold,
expression prediction or overall biological pass has been added.
