# F41: declared base-editor reference states

Integrated checkpoint: 2026-11-28. Native and WASM checks pass, with all 13
base-editor core fixtures, two binary fixtures and two actual MCP SDK fixtures
passing. [Combined source/execution evidence](verifier-and-design-current-integration.json).
Graphical acceptance remains open.

The existing CRISPR owner now exposes `sequence::crispr::base_editing`. It reuses
fresh canonical guide/PAM reconstruction and the existing translation/CDS owners;
there is no new guide matcher, codon decoder, empirical score, executable or plan
representation. The owner contains the current profile catalogue and its revision.

The builtin is deliberately named **BE1**, not BE3 or a universal CBE: the primary
2016 study describes a typical window at protospacer positions 4–8, counting from
the PAM-distal end, and describes cytidine deamination and eventual C→T/G→A
substitutions. This implementation declares conditional resolved reference states;
it does not simulate the uracil intermediate or repair. [Komor et al., Nature
533:420–424 (2016), Results and Figures 1–2](https://doi.org/10.1038/nature17946).
No sequence-dependent rate table, source implementation or empirical training data
was copied. The new Rust implementation remains under the repository MIT licence.

`EditorProfile` captures an exact compatible guide/PAM geometry revision, nominal
window, conversion, edition, rationale and citations. Altered builtin profiles are
refused; custom declarations need a distinct id and explicit content pin. A
citation does not authenticate a custom window or confer calibrated activity.
Three-prime-PAM geometry is the current supported coordinate profile.

`Request` binds the exact current molecule, freshly selected `GuideCandidate`,
`GuideOptions`, editor profile, intended guide positions, output bound and optional
CDS assumption. Guide positions are **one-based in biological guide orientation**;
source positions are **zero-based forward-reference calls**. The canonical guide
owner supplies actual positions across reverse strands and circular origins.
Consequently oriented C→T on a reverse guide is forward-reference G→A. Nominal cut
marks are never used to infer editing positions or actual cleavage.

The report lists each concrete editable site once, then enumerates all subsets in
deterministic bit-mask order, including the unedited reference. It retains intended
satisfaction and explicit bystander positions; an empty intent declaration produces
nullable satisfaction and does not silently invent an intended path. An intended
position outside the declared editable set fails **model reachability**, not the
possibility of biological editing beyond that model.

`finite_enumeration_complete` means only that all declared subsets were retained.
`biological_outcomes_exhaustive` is always false. No independence, equal frequency,
probability, efficiency, isolation, repair chemistry or functional-expression claim
follows from the subset enumeration. Outside-window edits, other substitutions,
indels and off-target outcomes are outside this finite set.

Optional `CdsRequest` reuses the current explicit `AnnotationAssumption`. Its exact
source/annotation pins, review flag, declared code, codon start and biological path
are checked through existing owners. The current consequence coordinate profile
supports an exact range or its complement; joined/fuzzy/partial/exception paths,
stale or review-required annotations, mismatched codes and oversized extraction
remain explicitly unresolved while window states remain available. Source labels
are not proof of annotation truth.

Affected ordinary codons are bundled into one canonical translation per state.
The table's actual 64-codon interpretation and descriptor are pinned in the CDS
evidence. Output distinguishes synonymous/missense/stop-gained/stop-lost/retained
stop/unresolved calls, and preserves three original reference coordinates in
biological codon order. First-codon flags deliberately describe ordinary table
calls; initiation-M recoding, actual start use, protein function and expression
are not inferred. [NCBI genetic-code tables](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi).

Admission/output bounds:

- Current canonical CRISPR source/model/guide guards; one source and one guide.
- Request 64 KiB; window width at most 12; at most 12 intended positions.
- Exact finite count at most 4,096; requested retained count 1–4,096.
- Report 2 MiB, with summary/citation/hash headroom reserved before retaining states.
- One optional CDS extraction at most 1M bases; at most 12 affected codons per state;
  no full genome or full CDS copy per reference state.
- Restored content/shape validation is documentary. `validate_for` freshly replays
  source, guide, declarations and consequences; imported reports never grant edit
  or chemical authority.

Thirteen source fixtures cover independent literal subsets/bystanders; reverse
G→A and circular-origin coordinates; explicit output and byte truncation;
unreachable/no-declared intent; source/profile/guide/report tampering; forward and
reverse ordinary stop-gain codons; stale/review-required/joined/wrong-code CDS;
required-null/current limits; declared A→G geometry; and ambiguous-guide refusal
through the existing selection owner. Worker/UI/MCP integration belongs to the
same parent source batch; no separate build was initiated for these files.
