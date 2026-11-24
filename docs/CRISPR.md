# Guide geometry and selected-reference off-target search

`core::crispr` is an independently implemented embedded Rust engine
(`crispr-1`). It provides bounded guide/PAM enumeration and ungapped Hamming
search against explicit immutable reference records. It does not call an
external aligner, borrow another cloning application's implementation, compute
an empirical efficiency score, or create physical nuclease-cleavage products.
Editing-template construction and independent donor checks are documented in
[EDIT_TEMPLATE.md](EDIT_TEMPLATE.md).

## Declared nuclease profiles

The SpCas9 preset declares a 20-base protospacer followed by `NGG`, with nominal
strand-cut offsets 17/17 measured as boundaries from the protospacer's oriented
5′ start. The primary mechanistic reference is [Jinek et al. (2012)](https://pmc.ncbi.nlm.nih.gov/articles/PMC6286148/).

The Cas12a preset explicitly selects a 23-base protospacer preceded by `TTTV`,
with nominal offsets 18/23. This is a **declared planning geometry**, not a claim
that every Cas12a ortholog uses exactly this PAM, guide length or pair of cuts.
[Zetsche et al. (2015)](https://pmc.ncbi.nlm.nih.gov/articles/PMC4638220/)
characterize T-rich upstream PAMs and differences among Cpf1/Cas12a enzymes.
[Cofsky et al. (2020)](https://pmc.ncbi.nlm.nih.gov/articles/PMC7286691/)
show cleavage distributions, trimming and strand gaps. The two nominal
boundaries neither encode that heterogeneity nor establish a particular
physical overhang or repair outcome. Other PAMs outside the selected profile
are not searched.

A custom profile declares its own identity, version, PAM, PAM side, guide
length, nominal offsets and optional HTTPS sources. Lengths are bounded to
12–64 bases and PAMs to 1–16 IUPAC DNA calls. Custom profiles require a distinct
identity and an explicit custom flag. An altered built-in cannot silently keep
its built-in identity. All reports retain the complete profile and its hash.
Custom parameters are caller declarations, not curated empirical activity data.

## Guide enumeration and original coordinates

Both source orientations are searched, including circular-origin crossings.
A complete guide plus PAM occupies at most one source turn. Guide and PAM
position vectors contain exact original-source base indices in their oriented
5′→3′ order. Reverse-strand vectors therefore decrease, with modulo wrapping
for circles. The guide string is the oriented DNA protospacer representation;
scaffold RNA and chemical synthesis modifications are not modeled.

Nominal cut fields are original-reference **boundaries**, not base indices.
Reverse orientation swaps the profile's strand offsets before mapping them
back to the original reference. Circular boundaries are modulo source length;
unavailable linear boundaries are absent. A nominal boundary alone supplies no
physical cleavage or product authority.

Possible IUPAC PAM matches are retained. `pam_certain` means every represented
base possibility satisfies the declared PAM. Ambiguity at an unconstrained `N`
position can still be certain; ambiguity that might violate a constrained PAM
position cannot. Ambiguous protospacers remain diagnostic candidates with no
GC/homopolymer measurement and are not selectable. Selectable guides must have
concrete A/C/G/T protospacers, certain PAMs and passing explicit filters.

GC bounds are evaluated by exact integer cross-products against the guide
length. `gc_bp` is a floored basis-point display value; it is not used to round
a marginal candidate into the permitted range. The homopolymer filter counts
the longest identical-base run. These are explicit deterministic filters, not
predictions of efficiency, expression or laboratory performance.

Candidate IDs bind the exact source revision, complete profile and original
coordinates/orientation. Changing GC or homopolymer thresholds can change
eligibility without changing a candidate's ID. `validate_candidate` freshly
reconstructs the local candidate and requires exact equality and eligibility;
it does not trust an imported guide report or require an entire new search.

## Named immutable reference collection

Each selected reference has a unique stable role ID, exact molecule snapshot,
and optional build/accession labels. All are included in request/report
bindings. Two roles containing identical DNA remain two records. Labels and a
free-text collection scope never establish chromosome completeness, assembly
coverage or genome-wide absence. The engine reports completeness **only for
the exact selected records under the declared Hamming/PAM model**.

The query is a named concrete DNA guide of the declared profile length. The
model examines both orientations and circular windows, with at most six
mismatches. It excludes indels, bulges, chromatin, unrepresented haplotypes,
structural variation and PAM alternatives outside the profile. A zero-hit
complete selected-record result does not establish biological specificity.

Intended-target tagging requires an exact tuple: reference role and revision,
guide position vector, PAM position vector and orientation. It removes no
hits. Equal DNA at another position or in another selected role remains an
ordinary additional hit. Stale intended revisions and malformed geometry are
errors rather than sequence-wide exemptions.

For ambiguous reference calls, each hit retains minimum and maximum possible
Hamming mismatch counts. A query-compatible ambiguous call contributes zero
to the minimum and one to the maximum unless it is the exact concrete query
base. A hit is possible when its minimum is within the threshold and its PAM
is possible. It is definite only when its maximum is within the threshold and
its PAM is certain. No guessed base or probability replaces the ambiguity.

## Embedded index and completeness

For concrete references, split the query into `d+1` disjoint nonempty blocks
for mismatch bound `d`. Any length-preserving match with at most `d` mismatches
must match at least one block exactly. The implementation uses rolling two-bit
keys to find those query-specific seed occurrences on both orientations,
including seed windows crossing circular origins. Seed offsets determine full
candidate windows; a bounded ordered set deduplicates them before exact full
Hamming/PAM checking. This is an algorithmic pigeonhole index, not a biological
PAM-proximal seed score. Query seed lengths are part of the declared search
algorithm, not empirical mismatch weights.

A reference containing IUPAC ambiguity uses a full window scan instead. This
preserves possible matches that a concrete-only seed index could miss. A fully
completed index plus verification/fallback is exhaustive for the supported
ungapped model on the selected snapshots. Index, comparison or hit limits
produce explicit incomplete reports. Missing candidates after a limit never
establish absence or uniqueness.

## Bounds and authority

Source shape and byte preflight runs before hashing. It bounds rich-source
nodes and location depth, enzyme/profile strings, role/provenance labels and
reference count. There are at most 64 selected records, a 20-million aggregate
base ceiling, and a separate **16 MiB aggregate serialized-source ceiling**.
The byte ceiling can reject a plain sequence collection before the base ceiling
is reached; this is not a promise that 20 million bases always fit. Each
sequence also retains the molecule model's ten-million-base limit.

Guide search has at most 4,096 retained candidates and a hard 200-million
comparison budget. Off-target search has at most 4,096 hits, two million seed
entries and 500 million comparisons; options normally choose smaller limits.
Seed-entry accounting bounds considered seed-occurrence/offset entries before
window deduplication, so allocated candidate count cannot exceed it. Reference
classification, rolling index scans, PAM checks and Hamming verification share
the comparison budget. Candidate/output structures have bounded lengths and
16 MiB serialized reports. Index memory is bounded separately from source and
output bytes. No whole-human-genome scale is claimed by this selected-record
profile.

`GuideReport::validate` and `OffTargetReport::validate` check bounded historical
structure and cross-field consistency. Current binding helpers additionally
require the exact current engine, source/query, options, profile and reference
provenance, without rerunning the search. Imported reports remain historical
claims. `replay_guides` and `replay_off_targets` freshly recompute and require
exact report equality. Donor construction must use fresh local candidate
validation, and browser/MCP flows must retain the explicit selected-record and
nominal-geometry scope.
