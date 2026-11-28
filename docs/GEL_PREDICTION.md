# Reference gel prediction and diagnostic digests

`collagene_core::output::gel` independently implements reference-only lane preparation, bounded gel diagrams and diagnostic digest ranking. It invokes no external scientific tool, image processor or observation service. OpenGel formats/integration and predicted-versus-observed comparisons remain deferred.

## Calibration and physical assumptions

The illustrated migration coordinate is a piecewise linear interpolation **against the logarithm of linear dsDNA size** between user-supplied ladder points. Each point gives size in bp and normalized distance from the well in the range 0–1. Sizes must increase while distance strictly decreases. No extrapolation is made outside the supplied size interval. The default 100/1,000/10,000 bp ladder at 0.9/0.5/0.1 is deliberately synthetic: it is an illustrative scale, not a measured agarose gel or a global physical mobility law.

Restriction-fragment sizing by reference electrophoresis has an experimental foundation in [Helling, Goodman and Boyer 1974](https://pubmed.ncbi.nlm.nih.gov/4372397/). More general size/mobility relationships can be nonlinear, as investigated by [Willis, Willis and Holmquist 1988](https://pubmed.ncbi.nlm.nih.gov/3154955/). This module makes only the stated local interpolation assumption; it does not infer agarose percentage, field strength, run time, buffer, temperature or staining chemistry from a size list. The calibration label and all coordinates remain in exported reports.

Uncut circular DNA is not positioned on the linear ladder. A circular reference does not establish its supercoiled, relaxed or nicked/open-circular state. The original conformation-dependent work, including [Aaij and Borst 1972](https://www.sciencedirect.com/science/article/pii/0005278772904261) and [Serwer and Allen's conformation/field study](https://pubmed.ncbi.nlm.nih.gov/6370305/), motivates that explicit limitation rather than an invented fixed mobility multiplier. Known unmodeled physical-intermediate or single-strand evidence also prevents silent treatment as ordinary intact duplex DNA.

## Actual source-derived fragments

`LaneInput` binds a molecule and one preparation:

- `Uncut`: a linear intact DNA reference can be positioned; circular conformation stays unplaced.
- `Digest { enzymes }`: the existing two-strand complete-digest engine reconstructs real fragments from the embedded enzyme definitions. Methylation, incomplete digestion, star activity and supplier conditions are not inferred.
- `Pcr { forward, reverse }`: the existing exact-PCR enumerator reconstructs possible products, including tails and competing oligo orientations. Multiple alternatives are not presented as known accumulated yields or a uniquely amplified product.

Both strand lengths are retained. For staggered-end fragments, `(top nt + bottom nt)/2` is used as an explicitly approximate mass-equivalent dsDNA size; short-overhang mobility effects are not calculated. Blunt fragments retain their exact bp counts. Failed/unsupported lanes carry an error, not a fabricated empty negative result. Fragments outside the ladder and circular fragments remain in the full lane table with no drawn position.

An optional loaded mass in ng allows a complete digest's total mass to be allocated in proportion to nucleotide count, assuming one copy of each fragment per source molecule. Equal-length fragments keep separate table entries even when they form one visible band. A single PCR product may receive its declared loaded mass; masses across alternative PCR products remain unknown because yields are not predicted. Uniform SVG strokes do not claim intensity or exposure calibration. There is no invented detection threshold.

## Co-migration and diagnostic ranking

`Calibration.resolution` is a declared normalized migration span, from 0.001 to 0.2. Sorted fragments are grouped only while a band's entire min–max span fits that threshold. This complete-link rule does not merge a chain whose adjacent gaps are small but whose endpoints are well separated. The illustrated band location is the group's migration midrange, and every member fragment remains in `Band.fragment_indexes`.

`diagnostic_digests` accepts 2–16 immutable candidate documents and explicit single/double enzyme choices. It evaluates complete-digest lane profiles, compares resolved band counts and sorted band-center separations against the supplied resolution, and records distinguished, indistinguishable and unknown candidate-index pairs separately. Profiles with unsupported chemistry, incomplete model evidence or unplaced fragments cannot establish positive or negative identification evidence. Intensity/copy-count differences do not establish discrimination under this positional model.

Ranking prefers more distinguished pairs, fewer unknown pairs, larger minimum reported diagnostic separation, then fewer enzymes and lexical enzyme names. A different visible band count receives a count-separation sentinel of 1.0; this is a deterministic rank convention, not a physical separation distance. Results bind source revisions, calibration, enzymes and budgets. They do not establish whole-plasmid sequence identity, experimental sensitivity or a globally optimal laboratory assay.

The requested enzyme-set search has explicit work/evaluation limits. A truncated search marks `exhaustive_search=false` and warns that the retained digests are best only among evaluated sets. Limiting retained results does not itself mean an incomplete search: all allowed sets may have been scored before retaining the best bounded subset.

## Limits and web workflow

The initial core accepts 1–32 lanes, 2–64 calibration points, at most 1,024 materialized fragments per lane, and 32 MiB of source bases and streamed serialized source evidence per gel call. Input label/mass/calibration validation runs before prediction. Reports/SVGs are bounded to 16 MiB; the binary worker envelope can impose a tighter limit when the active context and referenced documents are repeated in the same request. Diagnostic searches allow 1–16 embedded enzyme names, at most 256 evaluations, 20 retained candidates, and up to 200 million charged source-base/enzyme visits; defaults use 128 evaluations, 10 retained candidates and 50 million visits. These bounds are application limits, not recommendations for experimental gels.

The web **Gel** panel configures source lanes, uncut/digest/PCR preparation, optional loaded mass, ladder points and resolution. It runs predictions and ranking in Rust/WASM workers and checks the active document plus every referenced lane/candidate revision before accepting results. Changing a bound source, calibration, lane configuration or diagnostic options invalidates obsolete reports. SVG export is a standalone vector diagram with embedded calibration/model/source-revision metadata; JSON retains all returned fragment details, assumptions and errors. Accessible lane tables show up to 200 fragments per lane with an explicit preview count; JSON keeps the full bounded lane. Diagnostic candidate choices and source-index pairs are explicit, and a ranked enzyme set can populate the gel's lanes for inspection.

## Validation

`cargo test -p collagene-core --test restriction catalog::gel:: --locked` uses hand-derived blunt/staggered fragment lengths, independent log-coordinate midpoint calculations, mass accounting, circular versus cleaved references, actual tailed-PCR reconstruction, non-transitive co-migration fixtures, matched versus unknown diagnostic profiles, intensity-independent indistinguishability, budget limits and hostile-looking SVG labels. No other cloning application supplies expected results.
