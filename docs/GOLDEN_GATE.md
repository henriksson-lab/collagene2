# Golden Gate / Type IIS reaction geometry

`collagene_core::golden_gate` independently implements complete Type IIS digestion followed by bounded exact assembly enumeration. It invokes no scientific executable, copies no cloning application code and assigns no experimental success probability. Recognition/cut facts are represented directly; the algorithm reuses the independently tested physical duplex digestion/ligation primitives.

The method follows the restriction-site elimination principle described by [Engler, Kandzia and Marillonnet (2008)](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0003647). The factual starter profiles are BsaI `GGTCTC(1/5)`, BsmBI `CGTCTC(1/5)`, BbsI `GAAGAC(2/6)` and SapI `GCTCTTC(1/4)`, from [NEB's recognition table](https://www.neb.com/en-us/tools-and-resources/selection-charts/alphabetized-list-of-recognition-specificities) and [Type IIS table](https://www.neb.com/en/tools-and-resources/selection-charts/type-iis-restriction-enzymes). Parenthesized offsets are measured after the motif; internal `Enzyme` offsets are measured from its first base. These facts do not supply supplier buffers, methylation rules, temperature schedules or current commercial formulation behavior.

## API and immutable source evidence

`simulate(&[Molecule], &GoldenGateOptions)` accepts intact, unambiguous double-stranded reference inputs. Imported metadata cannot establish phosphorylation, repair nicks or authorize flattening a physical intermediate. `simulate_reconstructed(inputs, substrates, options)` is the explicit operation-graph boundary: callers supply already reconstructed `Duplex` values, with top sequence/topology checked against their corresponding source snapshots. It must not be used to interpret a JSON chemistry assertion in imported metadata as a reconstructed state. Existing duplex validation and digest rules reject unsupported nicks, crossing cuts or insufficient double-stranded substrates.

`GoldenGateReport` records source indexes/revisions, all complete-digest fragments, original digest indexes, top and unwrapped bottom source boundaries, actual 5′→3′ strands and end structures. Uncut circular sources remain explicit `carryover_sources`; they are not quietly counted as newly assembled products. The fragment pool includes carrier/backbone and terminal fragments, not only fragments a design wizard expects to use.

Candidate paths contain `OrientedFragment { fragment, reverse }`, referencing the complete pool. This preserves the data needed to independently reconstruct strand geometry and project source annotations. The module returns physical duplex products, not an annotated document or an asserted target match. The parent graph retains original annotations/evidence and must project them through the recorded digest/orientation maps. Unknown ancestry cannot be inferred from matching sequence text.

## Finite reaction model and completeness

The initial model uses **at most one copy of each complete-digest pool fragment**. It explores unordered subsets, all compatible orders and, by default, both physical fragment orientations. Minimum fragment count defaults to one so self-ligation/background candidates are visible. Circular paths are anchored at the lowest pool index to remove cyclic path duplicates; global reverse-complement orientations remain separate directed candidates so target and annotation mappings can be inspected.

Every junction must be exactly complementary, cohesive and fully sealable under known 5′-phosphate/3′-OH states. Products do not silently acquire nick repair. Circular outputs require both strands to close; linear outputs retain their actual terminal chemistry and strand extents. Original free ends of linear input molecules have unknown chemistry, so those boundaries cannot be asserted to seal.

`enumeration_complete` concerns this finite pool and the requested orientation scope. State/output/materialization caps or an extendable path stopped by the fragment-depth limit make it false, with explicit reasons. It never establishes real-mixture completeness: `exhaustive_mixture_search` is always false. Repeated-copy concatemers, additional molecular copies, partial digestion, star activity, ligation mismatch fidelity, kinetics, concentrations, selection/host effects and experimental outcomes remain outside this model.

`directed_sequence_key` groups equal circular reference rotations while retaining reverse orientation. `physical_sequence_key` additionally groups global reverse complements for circles; linear grouping compares actual full duplex geometry and terminal states under whole-molecule reversal. Equal sequence keys do not establish a unique original source or annotation history. All distinct source paths remain returned, even when their sequence bytes coincide.

## Recut checks

Each actual ligation product is scanned on both orientations, including circular-origin joins. This catches recognition motifs retained in fragments and motifs newly completed across junctions. `RecutSite` retains the recognition start/orientation, actual signed linear or modulo circular cut boundaries, whether the motif is fully double-stranded, and whether cleavage has a complete substrate under the physical model.

`restriction_site_free` requires no complete recognition motif in the product's top reference. It is deliberately conservative for linear overhangs/near-end motifs: a retained motif remains visible even when there is insufficient paired sequence for cleavage. It does not infer enzyme inactivity or empirical stability from absent cleavage geometry. Raw recuttable candidates remain available for inspection; they are not relabeled as site-free constructions.

`stable_physical_product_count` is the number of physical sequence groups among **returned recognition-site-free candidates only**. The field name reflects site-elimination under this model, not empirical reaction stability or a guarantee of unique construction. It does not include uncut carryover or omitted candidates, and must be interpreted with truncation/orientation flags.

## Bounds and validation

Limits are 32 input documents, 128 complete-digest pool elements, 32 fragments per path, 1,000 returned candidates and one million visited states. Defaults are 16 fragments/path, 128 candidates and 50,000 states. Stored pool/product strand sequence is capped at 16 MiB; reconstruction work is capped at 100 million strand bases. Unsupported pool/input bounds reject the operation; an output/work cap retains current evidence and marks enumeration incomplete. Additional shared model/binary/worker limits still apply.

`cargo test -p collagene-core --test assembly recipes::golden_gate:: --locked` covers hand-derived nonpalindromic BsaI strand diagrams, an exact two-part circular product, all donor circular origins and reverse complementation, junction-created recognition, self-ligation/palindromic ends, an independently enumerated 34-path subset/permutation oracle, output/depth/state truncation, unknown chemistry, uncut carryover, near-end noncleavable motifs, all four enzyme profiles, physical PCR-to-Golden-Gate graph inputs, a 32-donor/96-fragment assembly, the addressable 128/129 pool boundary, and independent rotation-equivalence classes over 192 directed candidates. Expected products are manually specified or checked through small independent rotation/permutation enumeration; another cloning application's output is not a correctness oracle.

## Web workflow

The Golden Gate panel selects explicit open templates, enzyme, output topology and finite search limits. It shows complete-digest pool origins/ends, directed candidate paths, retained recognition sites, carryover and truncation evidence. Downloaded JSON retains the full returned report; the panel pages through 50 candidates at a time.

**Review Golden Gate candidate** creates a typed plan bound to immutable input snapshots and the independently active target revision. It grants no verification status. The Plan panel must reconstruct that explicit candidate again, check current bindings/target and carry its method assumptions into HTML/PDF exports. Changing a selected source, target or calculation option invalidates the prior panel result; in-flight jobs also check request identity, generation and component lifetime.
