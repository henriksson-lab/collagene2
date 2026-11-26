# Declared physical homology assembly

This increment adds an independently implemented material model for ordered
homology assembly. It preserves the old sequence-only `ProfiledOverlap` operation
and its historical inability to supply downstream physical material.

The new operation consumes actual, oriented paired linear duplexes. It applies
explicit 5′ strand cuts, pairs the declared exact homologies, fills only gaps
supported by actual opposite template calls and available 3′OH groups, and seals
only eligible gapless 3′OH/5′ phosphate nicks. Both product strands are tiled
independently from surviving source calls and synthesized complements. No target
sequence is supplied to the kernel, and no PCR constructor is applied to assembled
letters. Surviving modifications are mapped; synthesized bases are unmodified.
Positive supported backbone cleavage exposes a new 5′ phosphate; zero cuts
preserve the original terminal group.

`ExactHomology` cuts each joining 5′ end by its declared overlap length. `Explicit`
accepts a left/right cut pair for every oriented input and rejects crossed or
exhausted source domains. Single-fragment circularization and single-linear-input
processing are supported. Circular origin stays at the first source's declared
reference origin. Linear `JoiningEndsOnly` explicitly exempts outer ends as a
model abstraction; `AllAccessibleEnds` retains the actual outer 3′ overhangs.
Neither silently creates a blunt outer product. Source annotations are projected
only after physical strands exist and their deterministic reference interval is
checked. Clipped features retain partial/fuzzy boundaries and require review.

The evidence retains source geometry revisions, removed source intervals,
copied/synthesized mappings, each joining gap, exact nick boundaries, terminal
chemical groups, seal decisions and bounded alternative-homology findings.
Ordered-path selection does not establish a unique reaction mixture, isolation,
yield or proprietary reagent behavior. Gibson/HiFi sequence guidance remains
advisory. In-Fusion's different 3′ processing/host-repair route cannot be selected
as this 5′ processing model; its reference interpretation remains available.
SLIC and broader host-repair chemistry remain separate unfinished work.

The current private nick representation cannot faithfully retain blocked or
unknown internal nick chemistry. Such paired products remain explicitly
incomplete, reference-only diagnostics with their full nick/chemical evidence;
the bridge withholds reusable material. Unsupported gap filling withholds the
product rather than inventing DNA. A sealed, faithfully represented product can
be finalized with its exact run, archived, replayed through `FromRun`, and used
for downstream PCR/digestion. Detached summaries and producer metadata cannot
restore that authority. Search caps remain sticky during replay.

The normal cloning dock exposes these assumptions in a collapsed strand-processing
control. Old drafts retain reference-only behavior by default. Draft10 preserves
unfinished and disabled processing settings. New worker86/job167 carries a
physical recipe preview; operation24 and physical evidence3 require candidate19
and native-plan17. Science22 semantics for existing operations remain unchanged,
so valid existing science22 material does not require an artificial retag.
MCP tool179 uses the same source/target-bound compiler and verifier. HTML/PDF
instructions retain the actual processing steps and unresolved findings.

Validation covers seven strand-kernel groups, four native integration groups,
two binary groups, one real MCP SDK group, eight library-history cases, twenty
historical binary cases, five existing PCR-ancestry cases and five web
persistence/mapping cases: 52 distinct passing cases. Integration includes an
independently written 40-base circle from 5′OH PCR products, native output
finalization/archive, exact saved-origin replay, PCR/EcoRI digestion/fill to an
independently written 50-base product, and fresh HTML/PDF exports. Native and WASM
compilation and formatting pass. No browser checks were run, and these new
controls await the next batched public build.

[Primary mechanism references](HOMOLOGY_MATERIAL_SOURCES.md) ·
[Original material-reuse audit](HOMOLOGY_MATERIAL_SCOPE.md)
