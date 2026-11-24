# Portable SVG sequence maps

`crates/core/src/svg.rs` provides `layout_map(&Molecule) -> Result<MapScene, Error>` and `render_map(&Molecule) -> Result<String, Error>`. These are portable Rust functions with no new dependencies, filesystem, network, subprocess, browser or graphics-runtime requirement. The standalone SVG can be downloaded directly or embedded into an HTML protocol; rasterization is a separate export boundary.

## Scientific coordinate scene

`MapScene` contains bounded typed `MapFeature`, `FeatureSpan`, `MapLabel` and legend data. The scene retains the molecule revision, alphabet, topology, feature source index, zero-based half-open ranges, per-span strand and fuzzy-end flags, ordered/unresolved/partial flags and annotation-review status. A renderer can reuse these data without re-parsing feature locations or extracting a different sequence.

Compound spans remain separate, including origin-spanning joins. Complement toggles strand and reverses compound-part concatenation order, matching the biological order of reverse-complement extraction. Exact nucleotide features receive one terminal arrow on the biological 3′ part, rather than an arrow implying that every exon is a separate feature. Fuzzy, order, unresolved or export-truncated features have no direction arrow; known indicative spans are dashed. Unknown portions are not guessed onto the backbone. Protein maps retain positional annotations without nucleotide strand arrows.

Circular origin is the boundary before base 1 at twelve o'clock; coordinates increase clockwise. Full-molecule features use two semicircular SVG arcs so they cannot disappear as a zero-length start/end arc. Linear feature positions derive directly from `start / length` and `end / length`; one-based residue/base ticks are placed at base centers. Circular visible ticks use the first base after each selected boundary. SVG `data-start`, `data-end` and `data-terminal` attributes expose the underlying exact reference boundaries. Drawing coordinates are rounded to three decimal places; this is a bounded figure, not a base-resolution editor for an arbitrarily long molecule.

The export uses the same feature-category palette as the current web renderer. Greedy tracks keep adjacent half-open intervals on the same lane while separating overlaps where possible. Label anchors use a real component of a compound feature rather than averaging an origin-spanning join onto the wrong side of the molecule. Circular labels use separated left/right columns; linear labels use collision-checked rows and compact vertical layout. Annotation-review flags add red outlines and a label asterisk.

The current WebGL renderer still owns its interactive mesh and viewport layout. This new portable scene supplies the export foundation and can be adopted by that renderer; the export does not claim that the two rendering backends already share all drawing code. Backbone topology describes the reference sequence. Duplex end chemistry, nicks, methylation, secondary structure and experimental feasibility are explicitly not visualized or established by a sequence map. Protocol verification carries that separate evidence.

## Bounds and visible fidelity notes

- At most 512 features and 4,096 mapped spans are drawn. Additional features/spans are explicitly counted or flagged in the figure. Affected partial features do not receive arrows.
- At most 32 feature labels are placed, with at most 16 on either circular side. Labels, titles and feature-type strings are abbreviated, and unusually wide text is clipped to its bounded label viewport. The document's original annotations/qualifiers are not changed.
- Four geometry tracks are available. Unavoidable overlaps and missing map labels are reported visibly, along with unresolved/fuzzy/order/review counts. Source features are included rather than silently discarded.
- SVG output is capped at 2 MiB. Exceeding that limit returns an error instead of silently returning a more severely truncated artifact.
- Input model validation and its location-depth/range limits apply before layout. Empty annotation sets and single-base/residue sequences are valid. Circular protein documents remain unsupported by the core model.

Each figure exposes the complete molecule content revision in `data-revision` and a visible provenance line. It contains no sequence-data mutation, asserted verification result, hidden external attachment or link to a network source. The native document remains the source of complete annotations and evidence; an SVG figure is not a lossless project archive.

## Safe markup and validation

Only fixed SVG drawing primitives and styles are emitted. Names, feature types and labels are escaped as XML text; invalid XML 1.0 control characters become U+FFFD. Qualifiers cannot create markup, event handlers, scripts, style sheets, external resources or permissions. Label clipping uses local nested SVG viewports and needs no potentially colliding global IDs when multiple maps appear in the same report.

`cargo test -p collagene-core --test svg --locked` covers independent scientific expectations for origin-spanning reverse joins, biological terminal arrows, half-open overlap/track behavior, full-circle arcs, uncertain/unresolved geometry, DNA/RNA/protein semantics, visible limits, malformed ranges, XML/control-character escaping and absence of executable/external SVG constructs. Tests inspect typed coordinates and parsed XML rather than treating another cloning application or a pixel screenshot as a correctness oracle.

Development previews were rendered from an independently authored synthetic 4,800-base feature fixture through CairoSVG to inspect circular and linear label layouts. CairoSVG is only a development inspection tool and is not an application/runtime dependency. The implementation itself uses the repository's model and standard Rust math/string operations.
