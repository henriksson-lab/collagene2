# Exact overlap profiles and advisory method guidance

`core::overlap_profiles` independently checks one ordered, explicitly oriented exact-homology candidate. It does not implement or invoke an external assembly tool, proprietary reagent chemistry, thermodynamic folding, reaction kinetics, product yield, cleanup, transformation, or empirical kit validation. The existing deterministic overlap reconstruction still constructs the actual candidate product and compares it with an independently bound target.

## API and identity

`check_fragments(molecules, overlaps, circular, options)` is an all-forward convenience. `check_oriented_fragments(molecules, reverse_flags, overlaps, circular, options)` supports explicit orientation without constructing new reverse-complement molecule copies. Both return `ProfileReport`. Source revisions refer to the **original** input documents, and the report retains the full orientation vector, overlaps, options, model and engine version.

A junction connects an oriented upstream suffix to the next oriented downstream prefix; a circular candidate also contains the final-to-first junction. Its two ordered 5′→3′ strings must be identical. `upstream_span` and `downstream_span` are zero-based half-open **original source-reference** coordinates with their orientation, including reverse inputs. GC fraction is a direct homology composition count, not a binding probability or a melting-temperature model. The checker does not guess alternate orientations or shorten an overlap to obtain a match.

Inputs must be linear, unambiguous DNA. Other alphabets, ambiguous bases and intersecting/exhaustive overlap arms are explicit unsupported states, matching the existing reconstruction profile. Zero/out-of-range overlaps or mismatching terminal strings produce `geometry_valid: false` with diagnostics. Physical accessibility, nicked/staggered intermediates and imported physical-chemistry evidence are assessed by the originating verifier graph; this sequence-only report never establishes those facts.

The report's `validate()` checks bounded structure, model identity, coordinate derivations, homology comparisons, GC values, guidance flags, occurrence fields and completeness. It does not establish that arbitrary serialized strings match the named source documents. Scientific consumers must recompute the checker from bound sources/options. No imported report becomes verification authority.

## Versioned length guidance

The factual table was reviewed on 2026-11-24 and is deliberately identified by a versioned model string. Guidance is **advisory** and separate from exact sequence geometry. A candidate can have exact compatible homology while being outside a published length recommendation. `recommendation_compliant` and `within_recommended` refer only to the implemented length screen, not every condition for a branded reaction.

| Profile | Declared fragments | Screened homology length |
| --- | ---: | --- |
| HiFi | 2–3 | 15–20 nt |
| HiFi | 4–6 | 20–30 nt |
| Gibson | 2–3 | 15–25 nt |
| Gibson | 4–6 | 20–80 nt |
| In-Fusion | 2 | At least the 15 nt recommendation |
| In-Fusion | Multiple inserts, 3+ fragments | At least the 20 nt recommendation |

The HiFi/Gibson values follow [NEB's optimization guidelines](https://www.neb.com/en-us/tools-and-resources/usage-guidelines/optimization-tips-for-nebuilder-hifi-dna-assembly-and-neb-gibson-assembly). Other fragment counts return no inferred length recommendation for those profiles. That is not a claim that other counts cannot assemble. DNA amounts, ratios, purity, reaction timing and other vendor conditions remain outside this sequence check.

The In-Fusion thresholds follow [Takara's homologous-overlap guidance](https://help.takarabio.com/36168/kb/article/109845/what-are-the-requirements-for-a-homologous-overlap-that-will-facilitate-a-successful-infusion-cloning-reaction). The checker does not infer an upper cutoff or guarantee performance for longer overlaps or arbitrary fragment counts. Applying the multiple-insert length threshold does not claim that a kit supports every fragment count allowed by the bounded sequence model. A one-fragment self-circularization has no inferred profile recommendation.

These facts are not an implementation copied from vendor or cloning software. Names identify the source of advisory compatibility guidance; no affiliation, reagent endorsement or complete proprietary reaction model is implied.

## Repeated-homology screen

When enabled, each declared upstream homology word is searched exactly against every supplied fragment on both strands. Occurrences retain original fragment indices, original reference intervals, strand and whether they correspond to the declared junction ends. Palindromic words are reported once per coordinate, canonically forward; that canonical representation cannot establish a unique physical strand orientation. Identical source revisions used twice remain separate ordered fragment instances.

Extra occurrences expose additional concrete sequence contexts. They do **not** enumerate alternative reaction products or establish that an extra site is kinetically accessible. Shorter microhomologies, mismatches, secondary structures, fragment-copy multiplicity, orders/orientations and reaction mixtures are not exhaustively explored. `exhaustive_alternative_search` is always false, including when all requested word searches complete. Disabling repeat screening records an explicit warning rather than making a uniqueness claim.

## Bounds

The checker accepts up to 1,000 ordered fragment instances, at most 16 MiB total source bases and 32 MiB source evidence with conservative per-entry overhead. It bounds source entries to 500,000 and location nesting to 64 before validation/hashing. Ordinary preflight validation and content hashing have a separate bounded source cost; comparison counters do not claim elapsed-time or whole-program CPU bounds.

Each declared overlap is at most 10,000 bases, with at most 1,000,000 declared homology bases across junctions. The comparison budget is 1–50,000,000 (default 5,000,000), counting actual exact-junction and motif-letter comparisons. A budget too small even for the specified exact junction lengths rejects the request rather than returning an unestablished geometry as valid. The shared occurrence budget is 1–10,000 (default 1,000). Exhausted valid repeat-scan budgets produce `complete: false` and explicit warnings. Omitted sites never establish absence or uniqueness.

Reference spans, source bindings and homology output are bounded; the report contains no complete source molecule copies. The implementation creates only bounded overlap/reverse-word strings rather than full reverse-complement reference copies. Report occurrence validation stops at the configured count before traversing excess entries.

## Independent evidence

`tests/overlap_profiles.rs` checks literal 15/25 nt overlaps and GC counts, all four source-orientation choices with original-reference coordinates, three-fragment circular closure, advisory range differences, invalid/unsupported states and declared work limits. An independent literal-string oracle enumerates every four-base DNA reference and all four input orientation choices (1,024 cases), checking both-strand occurrence coordinates, including palindromes. Serialized reports are structurally checked without being treated as source evidence.

## Application workflow

The **Cloning methods** inspector supports up to 12 ordered source rows, explicit reverse-complement orientation, one homology length per linear junction plus a circular closure when selected, profile choice and repeat-search limits. The current interface preserves settings across candidate review but invalidates evidence when source/options/workspace/target changes. It displays exact geometry, requested-search completeness and advisory length guidance separately; the full JSON retains all reported occurrences beyond the bounded on-screen preview. Unknown guidance is not displayed as a failed experiment.

**Review method candidate** creates a typed `ProfiledOverlap` operation. Plan verification replays source orientation, physical guards and exact reconstruction against the independent target, rather than trusting the displayed profile report. A short exact overlap may reconstruct the intended sequence while remaining outside guidance; an incomplete requested scan remains incomplete. HTML/PDF retain the method, orientation, overlaps and advisory/experimental limits. Source-derived known blocked 3′ ends require an explicitly modeled unblocking/resection mechanism in both legacy and profiled overlap paths.
