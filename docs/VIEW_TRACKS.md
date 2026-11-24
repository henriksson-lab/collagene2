# Bounded reference browser tracks

`view_tracks` implements F04's reference-coordinate track calculation in portable Rust. It uses the embedded codon model and factual restriction starter catalog. It does not run external tools. Drawing, WebGL, animation, worker scheduling, and rejecting results belonging to an obsolete UI document remain client responsibilities.

## Source lifetime and API

Create `PreparedTracks::new(Arc<Molecule>)` when the document changes. Preparation validates the molecule, hashes its revision, and builds a bounded interval index (at most one million spans). It owns an immutable shared source, so a cached index cannot silently refer to edited coordinates. `report(window, options)` neither clones nor hashes the complete source. `report_for_revision` also rejects an explicitly stale requested revision. `source_revision()` identifies the prepared source. `view_tracks(&Molecule, window, &options)` is a one-shot convenience that clones and prepares; it should not be used for scrolling a cached document.

All serializable configuration/report types derive JSON Schema under the core `json-schema` feature. `TrackReport::validate()` checks decoded structural invariants and budgets. Structural validation does not establish scientific agreement with a source: callers must match the report revision to their document. The report preserves effective options, engine/version, original window, indexed source size, charged work, completeness, and warnings.

`TrackWindow { start, end }` is nonempty, zero-based, half-open. Start is within the source. Linear windows end within the source. Circular windows may end beyond source length to describe one origin-crossing window, with width at most one source turn. Report parts always use physical reference coordinates within the source. A renderer can place a low-coordinate part in the next visual turn using the window. Cut positions are strand boundaries and may equal source length on a linear molecule.

## Options and bounds

Features, primer-binding annotations, enzyme recognition regions, actual enzyme cut positions, ordinary frame translations, validated known-CDS translations, and ORF predictions have separate switches. Primers are existing `primer_bind` annotations, not invented oligonucleotide designs. Source annotations are omitted from the display index. Feature indices and deterministic item IDs identify original evidence within the reported revision. Display labels are cached at preparation and capped to 512 Unicode characters. Only visible annotation parts are returned; total item coordinate parts are capped at 10,000 and warnings at 1,024. These caps are disclosed as incomplete.

Only names from `analysis::starter_enzymes()` are accepted. Those twelve factual profiles do not provide methylation, supplier conditions, star activity, or experimental cleavage evidence. Recognition matching preserves compatible ambiguity versus definite recognition. Type IIS cut boundaries are calculated from the profile offsets and orientation, independently of recognition region endpoints. Cuts outside a linear source are absent and warned about; circular cuts are normalized modulo source length. A cut-only viewport can include a site whose recognition region is outside that viewport. Recognition evidence is retained in the item's parts even if only cut drawing is enabled.

Frames are unique offsets 0, 1, and 2. Orientation selects forward, reverse, or both for frame translation and ORF candidates (and directed recognition hits); annotations preserve their recorded orientation. Genetic codes are the embedded supported codes, currently 1, 2, and 11. Frame translation always uses ordinary internal codons; it never infers an initiator from the start of a viewport. Known CDS uses its own recorded `/transl_table` and `/codon_start` through the existing CDS algorithm rather than overriding annotation evidence with the frame-track controls.

| Budget | Default | Hard maximum |
| --- | ---: | ---: |
| Track items | 1,000 | 1,000 |
| Emitted residue mappings | 4,000 | 10,000 |
| Window bases | 10,000 | 100,000 |
| Work units | 100,000 | 1,000,000 |

Invalid budgets/configuration are errors. Exhausted valid budgets produce `complete: false` and explicit warnings. An omitted or capped result never means that no feature/site/ORF exists. Work units charge index spans examined, frame-coordinate and codon-base visits, recognition comparisons, bounded full-CDS extraction/translation and qualifier metadata scanning, and conservative whole-source ORF scanning/materialization. Preparation validation/hash/index cost is disclosed separately as source bases/indexed spans. These are algorithmic work units rather than elapsed-time guarantees; ambiguity decoding has an additional bounded codon-table expansion factor. No full source-sized allocation is made by ordinary viewport frame translation. Known-CDS/ORF full-source work must fit the configured budget before allocation.

## Translation coordinates and scientific limits

Forward frame offsets start at reference zero. Reverse frame offsets start at the last reference base and proceed on its complement. Codons are selected when any of their three reference bases intersects the viewport, including codons beginning before a line boundary. Each visible residue records its index within the full oriented frame, amino acid, and three-base reference coordinate parts in biological order. Reverse parts retain ascending reference endpoints; their biological traversal direction is the item's reverse flag. Compound CDS codons may span multiple joined parts.

Circular ordinary frames are anchored to the reference origin: accept codon starts `frame + 3*k` within the source, and allow their final codon to cross the origin. For source lengths not divisible by three this is an origin-anchored finite display convention, not a claim of a unique continuous reading frame around an unlimited number of turns. The same codon retains the same mapping in windows on either side of the origin. At least three source bases are required for frame display. Origin-crossing biological CDS/ORF interpretations use their recorded location or bounded ORF start, independently of this display convention.

`KnownCdsPolicy::Validated` calls the existing reviewed CDS translation model. Reviewed-required, unordered, fuzzy, and mixed-strand locations do not receive invented coding frames. Unsupported qualifiers, stored-translation disagreements, partial translations, and other evidence limitations remain warnings from the existing algorithm. A known CDS's visible protein excludes the terminal stop following the established CDS report convention; ordinary frame tracks include `*`.

ORFs use existing `find_orfs` logic with ATG-only starts, definite stop boundaries, complete candidates, configured minimum amino-acid length, and the selected genetic code. They are candidates rather than evidence of expression. They are searched on the whole source so a cropped viewport does not invent starts or stops. Whole-source scan and worst-case bounded protein materialization are prepaid; chromosomes exceeding that budget have an explicitly omitted ORF track. A prediction cap can be reached before viewport/orientation/frame filtering, which is disclosed as incomplete. This does not silently reinterpret a cropped sequence as an entire chromosome. Definite terminal ORF stops are present in location footprints but absent from residue mappings/protein, matching the original ORF report.

`tests/view_tracks.rs` uses literal independent phase, reverse-complement, joined-CDS, circular-origin, and Type IIS coordinate oracles. It also checks toggles, stale revisions, malformed options/source/report, work/item/residue caps, and a chromosome-sized source with small-window frame computation.

## Web viewport

**Sequence tracks** exposes separate switches, enzyme names, strand direction, frames, genetic code and ORF minimum length. The canvas shares reference selection with the map; clicking an amino-acid symbol selects its mapped codon, and the accessible viewport list exposes returned items with paging. Compound items highlight all parts without constructing a contiguous editable interval. Overlapping annotations share a fixed number of display lanes; the list and downloadable JSON retain every returned item.

The web transport installs one immutable source into a retained worker and sends only revision/window/options for subsequent queries. Viewport tiles are quantized and requests debounced; source/settings/window changes reject obsolete replies. Source replacement and unmount terminate the prior worker. The reload control explicitly retries a failed session. This is a dedicated sequence-view cache, not a general worker pool or incremental genome store. Whole-source ORF work can exceed the budget even when ordinary viewport frame translation succeeds.
