# CRISPR on large immutable reference stores

The `genome` store and `genome_crispr` algorithms keep chromosome data outside ordinary plasmid documents. Streaming FASTA ingestion writes packed IUPAC owner blocks to a `.cgref` file or browser storage. An interval operation loads only the required blocks. A cooperative search retains one bounded owner window, its right halo, its query-specific candidate index and a bounded hit list. It never constructs a whole chromosome as a `Molecule`.

All public genome positions and lengths use `u64`. Intervals are zero-based and end-exclusive. Guide intervals own **forward genomic footprint starts**: the entire guide and PAM can extend to the right of the requested interval. A circular interval may have an unwrapped end past the contig length, with width at most one traversal. Individual returned guide/PAM coordinates are projected to original contig positions in the strand's reading order. Nominal cuts are zero-based boundaries, including the terminal boundary of a linear contig; circular cuts are reduced modulo the contig length. Interval reads return exactly their requested bases. Current interval operations allow at most 200,000 owner starts or read bases.

`GenomeGuideOptions.source_role` and `GenomeSearchOptions.reference_role` distinguish user-visible logical roles even when two imports contain identical DNA. Stable candidate identity includes the manifest, role, contig, declared profile, original coordinate vectors and orientation. Changing the interval or GC filter does not change an existing candidate's locus identity. Expected-component metadata and provenance belong to the exact manifest hash; canonical content identity separately describes the ordered-by-role DNA collection. A name such as a species or genome release is descriptive metadata, not an attestation that a provider's entire build was retrieved.

An interval report verifies the layout and hash of every supplied packed block, but does not read the whole contig to verify its declared canonical DNA hash. Its warnings explicitly state that scope. Imported files receive a separate ordered semantic readback before publication. Reopened browser metadata cannot restore that live verification authority. A whole-manifest search independently verifies the actual canonical DNA of every owner block in order and checks each contig's canonical hash at EOF.

Search assigns every forward footprint start to exactly one owner block and includes the right halo, wrapping circular origins when necessary. Both orientations are analyzed. Bounded linear windows use the independently implemented ungapped Hamming/PAM search, including its query-specific pigeonhole seed index and full compatible fallback for windows containing ambiguous calls. A final short block can require both the end and beginning of a circular contig as halo; repeated block requirements are fetched once. Hits in a neighboring owner's halo are discarded by ownership, not by sequence equality. Whole-contig resident indexes and persisted FM indexes are later performance work.

`GenomeSearchSession` is deliberately neither serializable nor clonable. Its private cursor and canonical hash state establish actual traversal. Each advance supplies exactly the next owner's requested blocks in the requested order. Invalid, missing, reordered or tampered blocks poison the session. Serialized progress, coverage flags, receipts or historical reports cannot resume it. Cancellation, work/index/output limits and failed semantic checks prevent complete coverage claims. Pages contain only bounded retained hits; reading all retained pages does not remove a truncation limit.

The progress fields separate three facts:

- `source_integrity_verified`: all actual contig DNA has passed canonical semantic validation in this live search.
- `exhaustive_for_manifest`: every manifest owner was searched under the declared Hamming/PAM model without limits or cancellation, and source integrity passed.
- `complete`: exhaustive manifest search plus all declared expected components being present. An unenumerated expected build does not establish whole-build completeness.

The current `genome-crispr-2` engine explicitly uses zero bulges and the selected
design PAM. The separate selected-record bulge and screening-PAM controls do not
apply to this API; genome requests have no such options.

None of these fields predicts biological specificity, cleavage activity or repair. The current model excludes bulges, population variation, chromatin and alternative PAMs outside the selected profile. Intended-site tagging requires the exact reference role, manifest and canonical content revisions, contig, guide/PAM vectors and orientation. It never removes the hit or suppresses identical DNA elsewhere.

Report `validate` and `validate_bindings_for` methods check shape and request binding only; they do not replay analysis or grant positive scientific authority. Native/browser job registries must keep their private sessions and bind every result to the current reference role and request.

Independent tests exercise linear owner seams, short final blocks, both strands, small circular genomes, circular origins, Cas12a upstream PAMs and unequal nominal boundaries, compatible ambiguous calls, exact intended-role tagging, missing/tampered blocks, forged canonical source hashes, cancellation and resource limits. A streamed FASTA over 33 MiB contains a contig over 10 million bases and verifies the final global interval and every owner during full search, without raising ordinary document/file limits.
