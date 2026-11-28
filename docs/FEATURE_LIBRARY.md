# Feature and part libraries

The independently implemented `feature_library` module supplies user-created, versioned DNA feature templates and deterministic automatic-annotation proposals. It includes no proprietary or scraped annotation database. Application code is MIT licensed, and matching invokes no external software.

## Library records and source evidence

A schema-1 `FeatureLibrary` has a name, a user-supplied version, and immutable `LibraryEntry` records. `add_feature` extracts an existing document annotation in biological order: exact reverse/complement, joined and origin-spanning locations are supported. The DNA template, original feature kind/qualifiers, source name/content revision/feature index/location, and source review flag stay together. Re-adding the identical feature is idempotent. Source documents are not copied into every entry.

Unresolved, fuzzy and `order` locations are rejected. Joined source parts become a concatenated template; the matching model does not infer introns, new exon positions or partial genes. RNA/protein templates and references are outside this initial DNA profile.

Entry revisions are BLAKE3 hashes of a domain-separated typed record, and IDs contain the full revision hash. Library revisions bind the schema, name, version, ordered entries and source claims. Editing any bound entry field without updating its content identity fails validation. These hashes detect changed content; **imported source names/revisions are provenance claims, not authenticated proof that the original document exists or its annotation is correct**. Source material must be reviewed independently when that matters.

`save_library` writes bounded canonical typed JSON. `load_library` validates the schema, fields, hashes, locations, exact template lengths and total sizes before use. Unknown record and nested model fields, duplicate known fields, duplicate entry IDs and noncanonical template strings are rejected rather than silently discarded. Normal whitespace or JSON object-key ordering does not alter biological semantics or content revisions. Libraries fit ordinary native workspace attachments; the web application retains immutable snapshots and selects an explicit library revision.

## Matching and interpretation

`scan` visits every permitted fixed-length placement in library-entry, forward-reference-start, then forward/reverse order. Linear references require the entire template to fit. Circular references include origin-spanning placements but permit at most one turn; templates longer than the reference cannot fit. Reverse proposals use an explicit complemented location. Palindromic templates retain both orientation proposals, with a warning that sequence alone cannot resolve their orientation.

Matching uses literal IUPAC allele-set compatibility:

- A concrete equal A/C/G/T pair contributes a certain match.
- Intersecting ambiguity sets contribute an uncertain compatible column.
- Disjoint allele sets contribute a definite mismatch.

The explicit `max_mismatches` bound counts definite mismatches. It does not turn ambiguity into known identity or bound every possible allele combination. Identity is reported as `[concrete matches / length, (concrete + ambiguous-compatible columns) / length]`. Coverage is exactly 1 because this model searches only whole-template ungapped placements. Partial/homology alignments, indels, sequence-specific functional prediction and empirical validation are not performed.

All retained locations and orientations stay separate, including repeats and palindrome alternatives. Entry/proposal IDs bind exact source, library, options, coordinates and engine version. Report completeness describes only the declared whole-template search. Reaching comparison, hit or serialized-output limits marks it incomplete; unfinished alignments and unsearched positions cannot become annotations or evidence of absence. Longer-than-reference templates have no permitted placement under this model and produce an explanatory warning.

## Reviewing and applying proposals

All inferred annotations have `needs_review = true`, including exact, concrete matches. A sequence match alone does not validate the source annotation's biological function. Descriptive/display qualifiers are copied conservatively; source translations, CDS phase/code, database/protein identifiers, coordinate-bearing qualifiers and unknown semantics remain in the library evidence rather than becoming current target assertions. Every proposal records the omitted qualifier names.

`apply_selected` requires the exact target and library revisions and reconstructs the complete bounded scan using the recorded options/engine. Any altered report/proposal, stale revision, unknown selected ID, duplicate selection or existing equivalent same-kind/name/location/orientation annotation causes an error. Existing annotations are never overwritten, and rejected selections leave the input document intact. The entire selected addition is atomic, including cross-entry duplicates. The sequence and original source evidence are unchanged. Added features carry entry/library/source/match provenance, and the resulting document records source/library/report/proposal identities. No imported verification status grants scientific authority.

## Limits

| Item | Bound |
| --- | --- |
| Serialized library JSON | 8 MiB |
| Entries | 1,000 |
| Per-template DNA | 1,000,000 bases |
| Total template DNA | 10,000,000 bases; the JSON bound also applies |
| Original qualifiers per entry | 500 |
| Per-original-qualifier value | 1 MiB; aggregate library size still applies |
| Comparison budget | 20 million default; 200 million maximum |
| Retained proposals | 1,000 default; 10,000 maximum |
| Definite mismatches permitted | 0 default; 1,000 maximum |
| Serialized report/current document | 32 MiB |

The budget counts compared bases, including failed placements. Setup, source validation, reference-feature indexing and bounded evidence serialization also consume time. Existing annotation geometries are indexed before matching, avoiding a full annotation-list scan per candidate. Serialized evidence accounting is linear in retained/omitted evidence, rather than repeatedly serializing every report prefix. These are bounded-work/storage rules, not a claim that process RSS is limited to the serialized byte count. Options and resulting completeness must remain visible to callers.

## Independent validation

`cargo test -p collagene-core --test sequence annotations::feature_library:: --features json-schema --locked` includes 4,096 independent small allele-set/orientation/topology/mismatch cases plus all 225 IUPAC symbol pairs and their reverse complements under both exact/mismatch budgets (450 additional cases). The oracle uses explicit allele lists, manually defined complements and forward slices, rather than the scanner's bit masks or coordinate builder. Separate fixtures cover every circular origin/orientation, joined source extraction, palindrome alternatives, repeated placements, definite versus possible identity, partial-work and hit limits, large evidence/output caps, idempotent source extraction, strict JSON/tamper rejection, nonterminating identity fractions after JSON reopening, long original CDS translations, conservative qualifiers, stale/forged reports, atomic application and existing/cross-entry duplicates. No other cloning application supplies the correctness oracle.

The web Feature libraries panel retains up to 128 immutable library revisions in workspace attachments, supports individual or whole-document feature extraction, paginates entries/proposals and applies only explicitly selected matches in a Rust worker. Native `.cg2workspace` archives retain library definitions; a single-document `.cg2` or INSDC export carries references rather than the complete library. GenBank/EMBL exchange normalizes wrapped machine identifiers only when they validate as their declared hash/ID types, and carries the `collagene_review_required=true` marker for inferred/review-required annotations. Text export reports omitted library/source evidence; it does not confer verification.

There is no bundled common-feature catalog yet. Users create/import their own versioned libraries; curated public-source catalogs and broader indel/spliced homology searches remain planned work.
