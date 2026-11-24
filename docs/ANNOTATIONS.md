# Annotation coordinate projection

`collagene_core::annotations` is independently implemented Rust. It uses explicit coordinate maps and invokes no external tool. A `MappingSegment` maps a half-open source interval into a destination interval of equal length, either unchanged or reverse-complemented. Source and destination intervals may not overlap within one call. To place repeated copies of a template, project each copy separately and combine the resulting features deliberately.

`project_features(source, segments, destination_length)` assumes a linear destination. `project_features_with_topology(..., destination_topology)` also supports a circular destination. The caller constructs the destination bases independently; coordinate projection is not evidence that a claimed product sequence was reconstructed correctly. Check that each mapped sequence interval actually matches the chosen orientation.

## Biological semantics

Join/complement structure follows biological extraction order rather than sorting by destination display coordinates. A full retained exact feature can preserve its extracted sequence and original CDS qualifiers through reverse complementation or circular-origin changes. An origin crossing in a circle is distinct from a break across the ends of a linearized molecule. Previously coherent exon order/orientation cannot be silently replaced with a new permutation.

Partial, interrupted, fuzzy and ordered locations carry review flags. Clipped or interrupted ranges are marked partial so the model's exact-extraction function cannot silently assert a complete feature. `/translation` and `/codon_start` assertions are removed from active qualifiers for these reviewed projections and retained explicitly as source-evidence notes. No new CDS phase or translation is guessed. Retained exact features do not have their source translation assertions independently validated by this module.

Unresolved locations cannot be projected into new exact positions. Nonidentity maps return their original features in `Projection.unmapped_features`, with warnings. They are source evidence, not product annotations. A caller building a cloning graph must preserve those originals in the immutable input graph or metadata; silently discarding the vector would violate the API's fidelity contract. `dropped_features` counts known features wholly outside the retained source intervals; those are deliberately absent from the product.

The initial adapter does not parse/remap coordinates embedded in `transl_except`, `anticodon`, `rpt_unit_range` or `coded_by`. A retained feature carrying one of these rejects a nonidentity projection explicitly. A wholly omitted feature does not block projection. Mapping a full identity partition preserves these unchanged local coordinates. More complete qualifier remapping remains pending.

Location operators, partial positions and CDS phase semantics follow the [INSDC feature table](https://www.insdc.org/submitting-standards/feature-table/) and [NCBI partial-feature guidance](https://www.ncbi.nlm.nih.gov/genbank/feature_table/). The transformation policies above are application rules, not a claim that coordinates establish biological function or gene-expression feasibility.

## Document operations and source evidence

The module also provides `reverse_complement`, `rotate_origin`, `extract_region` and `change_reference_topology`. Circular extraction uses unwrapped half-open endpoints: on a ten-base circle, `start=8, end=13` retains five bases across the origin. An extraction is always a linear sequence document. Setting a circular origin is a reference-coordinate operation. Changing declared reference topology is not physical ligation, nick repair or evidence of experimental circularity.

Derived documents move exporter-active original GenBank/EMBL headers into explicit `annotation_mapping:*` source-evidence metadata. This prevents stale reference ranges from leaking into a new annotated export even when palindromic or periodic sequence bytes happen to be unchanged. Native files retain the original records and any unprojectable source features. Original headers are not automatically remapped current database annotations. Existing provenance URLs remain source provenance, accompanied by the derived-operation/source-revision fields. The metadata describes the latest operation; it is not the future complete replayable project-history archive.

Sequence-only transforms reject physical-end/nick evidence that requires the originating operation graph. Full trace reverse complementation validates the old attachment and switches its supported orientation, preserving original bytes and transforming its signal view through the trace API. Region extraction retains the full original attachment as evidence; a partial called sequence is not presented as a re-aligned chromatogram. Active trace topology changes are unsupported.

Limits are 4,096 mapping segments, one million traversal/intersection work units, the model's existing 100,000-feature/32-level location bounds, and bounded warning output. Document-operation results must fit the existing native serialization limit with all evidence retained. Real-world qualifier/schema coverage and larger annotation interoperability corpora remain future work.

## Validation

`cargo test -p collagene-core --test annotations --locked` covers independently specified biological extraction/orientation, circular-origin rotations versus linear cuts, partial exon deletion and new exon permutations, fuzzy/order/unresolved semantics, coordinate-qualifier failures, source-header preservation, native evidence round trips and trace orientation. These tests do not use another cloning application's output as an oracle.

## Reconstructed cloning products

The verifier maps PCR annotations after the actual forward synthetic tail and retains distinct PCR geometries even when their final sequence bytes coincide. A unique product sequence does not establish a unique annotation origin. Digest maps preserve circular source order while exposing new free reference ends. Reversed ligation maps old top coordinates using the actual bottom strand extent, then clips against the final linear reference or wraps around the final circular reference. Neighboring ligated strands can restore pairing for an old top overhang; clipping to the old paired intersection would incorrectly discard it.

A one-fragment overlap closure identifies two source homology copies. Its internal projection path tests injectivity separately for each feature: intact prefix/suffix annotations can remain exact, including origin-spanning features, while a feature covering both collapsed copies remains reviewed source evidence. It cannot become an asserted shorter exact gene. The public segment-map API remains one-to-one.

Product metadata retains unprojectable features by source revision rather than copying arbitrary input attachments into every intermediate. The immutable plan inputs/intermediates carry source evidence. An annotation absent from a digest fragment's top reference but present on its opposite-strand protrusion is currently not restored automatically when that protrusion later becomes top reference after reversal; the original graph and omission warning remain available. Restoring it requires explicit ancestral strand provenance.

`cargo test -p collagene-core --test annotation_verification --locked` covers manually specified PCR tail origins and competing geometries, circular digest clipping, reversed staggered ligation, linear/circular overlap closure, self-circular homology collapse, unresolved source evidence and bounded metadata duplication.
