# Ordered TOPO source graphs

An independent MIT Rust library for **physical strand/bond models**, using `collagene-model::Sequence` for validated DNA. It does not depend on `collagene-core`, other cloning software, an external scientific executable, a network service or a kit database.

This implements elementary transfer from a declared Vaccinia type-IB topoisomerase/Tyr274–DNA3′-phosphotyrosyl complex to an actual free DNA5′OH. A successful step adds one3′–5′ backbone bond and records the released enzyme-activation identity. The opposite strand is not automatically sealed. Source activation is a required chemical input, not an interpretation of a motif, plasmid name, GenBank annotation or supplier metadata.

## Supported physical geometry

All source strands are contiguous chemical chains written **5′→3′**. `PairSpan` intervals are antiparallel: `a[start_a+k]` pairs with `b[start_b+length-1-k]`. Existing nicks are represented by separate actual strands paired to a continuous opposite template. Free ends have explicit hydroxyl, phosphate, blocked, unknown or phosphotyrosyl states.

- `NickReseal`: actual paired donor3′ and acceptor5′ calls are adjacent across a continuous antiparallel template. This can seal the specified nick or close a single activated strand. It does not repair another nick.
- `Flush`: a terminal donor3′ activation and free acceptor5′OH meet at two flush duplex faces. The model joins the specified strands only.
- `ThreePrimeTA`: exactly one unpaired donor3′T faces an acceptor3′A overhang. The new T:A pairing is explicit, as is the donor→acceptor backbone bond.
- `DirectionalCacc`: an incoming blunt strand begins5′CACC. A donor duplex has the complementary four-base5′ overhang on its opposite strand. In the manufacturer's drawing this bottom sequence reads `GTGG` **3′→5′**; its actual stored5′→3′ sequence is `GGTG`. Four original incoming base pairs are replaced. Their fate is explicit `Displacement` evidence with actual5′→3′ source calls and their remaining5′ continuation. A duplex-anchored continuation gives the **attached3′ flap** drawn in the standard directional diagram. An exact four-base complementary chain can instead be completely displaced as a **free four-mer**; an unpaired continuation is classified separately as an unpaired-tail extension. These complete-displacement/tail-extension cases are declared physical graph branches, not claims of primary-diagram kit geometry or efficiency. There is no covalent excision, fill-in or nucleotide loss.

Both donor pentamers `CCCTT` and `TCCTT` are recognized by this declared Vaccinia profile. This does not produce an activation from an uncharged substrate. Uncharacterized enzymes, noncanonical donor contexts and unsupported geometries remain explicit unsupported outcomes. Unknown terminal chemistry or ambiguous DNA remains unknown. Incompatible declared base pairs are invalid inputs. Reported experimental mismatch transfer does not justify assigning rates or declaring unsupported mismatches unable to react.

## Current ordered API

`Step::Transfer(Transfer { donor: Call, acceptor: Call, geometry })` and
`Step::Process(processing::Action)` form one chronological transcript. Calls name
actual immutable source nucleotides. All motif, end and pairing checks use their
current private backbone and effective chemistry. Original source ends are not
substituted for newly exposed endpoints.

- `preflight(&Substrates, &[Step], &Limits)` bounds the full native execution.
- `batch_revision(...)` binds source hashes, ordered requests and limits.
- `simulate(...)` returns one canonical `Report` with `request`, typed `steps`,
  complete `product` and `all_requested_steps_applied`.
- `simulate_with_checkpoints(...)` constructs a private `Executed` capability;
  checkpoints use `after_steps`, exact whole-source scope and `step_map`.
- `projection::project` / `project_fresh` prove a component/original-call contour
  over the same current graph; only a native capability enters fresh projection.
- `replay(...)` freshly reconstructs every requested step. Readonly binding and
  source-partition validation do not manufacture material.

The private cursor holds backbone, pairing, displacement provenance, per-call
terminal overrides and consumed/released activation identities. Explicit ideal
flap hydrolysis creates3′OH/5′P; nick sealing joins only a supported actual
continuous-template3′OH/5′P nick. Neither step supplies activation,5′OH, missing
phosphate, polymerase extension or host repair. All source calls partition exactly
once, including the excised four-call component. Unsupported/unknown/blocked steps
retain the incoming state. Later steps describe that declared branch, without
probabilities or a mixture-success claim.

Source, steps, work, memory, output and report bounds are explicit. The work ledger
is shared across private orientation branches, source checkpoints and subsequent
fresh projection. A literal directional continuation fixture preserves52 calls:
48 in two closed24-call complementary contours and four in a separate component.
Captured or projected selection does not imply isolation. No supplied graph,
report, reference metadata or serialized cursor establishes source chemistry.

There is one current transcript and projection representation. Separate sealing
and processing report executors/readers are removed. Native application graph
capture, continuation and faithful material conversion are described in
[graph history](../../docs/TOPO_GRAPH_HISTORY.md),
[ordered processing](../../docs/TOPO_PROCESSING_ACTIONS.md) and
[the material bridge](../../docs/TOPO_PROCESSED_PRODUCT.md).

## Primary evidence

- [Shuman1994, JBC269:32678–32684](https://pubmed.ncbi.nlm.nih.gov/7798275/): prepared covalent complexes and acceptor-dependent strand transfer.
- [Cheng and Shuman2000, NAR28:1893–1898](https://doi.org/10.1093/nar/28.9.1893), [full article](https://pmc.ncbi.nlm.nih.gov/articles/PMC103307/):3′T donor/3′A acceptor geometry, a nicked complementary strand, alternative insert orientations and condition-dependent mismatch transfer. This engine does not convert those observations into a kinetic prediction.
- [Manufacturer directional mechanism, printed p2](https://documents.thermofisher.com/TFS-Assets/LSG/manuals/pcdna3_1dv5histopo_man.pdf): Tyr274/3′-phosphotyrosyl chemistry,5′OH attack, CACC invasion and retained displaced strand in the depicted product.
- [Manufacturer enzyme mechanism](https://www.thermofisher.com/us/en/home/life-science/cloning/topo/topo-resources/the-technology-behind-topo-cloning.html): Vaccinia(Y)CCTT specificity and prepared activated vectors.
- [Manufacturer PCR-cloning technical support](https://www.thermofisher.com/uk/en/home/technical-resources/technical-reference-library/cloning-technical-support-center/pcr-cloning/pcr-cloning-getting-started.html):5′-phosphate inhibition and vector terminal-chemistry distinction.

No implementation code from these sources or other cloning applications was consulted.
