# Explicit TOPO strand transfer

An independent MIT Rust library for **physical strand/bond models**, using `collagene-model::Sequence` for validated DNA. It does not depend on `collagene-core`, other cloning software, an external scientific executable, a network service or a kit database.

This implements elementary transfer from a declared Vaccinia type-IB topoisomerase/Tyr274–DNA3′-phosphotyrosyl complex to an actual free DNA5′OH. A successful step adds one3′–5′ backbone bond and records the released enzyme-activation identity. The opposite strand is not automatically sealed. Source activation is a required chemical input, not an interpretation of a motif, plasmid name, GenBank annotation or supplier metadata.

## Supported physical geometry

All source strands are contiguous chemical chains written **5′→3′**. `PairSpan` intervals are antiparallel: `a[start_a+k]` pairs with `b[start_b+length-1-k]`. Existing nicks are represented by separate actual strands paired to a continuous opposite template. Free ends have explicit hydroxyl, phosphate, blocked, unknown or phosphotyrosyl states.

- `NickReseal`: actual paired donor3′ and acceptor5′ calls are adjacent across a continuous antiparallel template. This can seal the specified nick or close a single activated strand. It does not repair another nick.
- `Flush`: a terminal donor3′ activation and free acceptor5′OH meet at two flush duplex faces. The model joins the specified strands only.
- `ThreePrimeTA`: exactly one unpaired donor3′T faces an acceptor3′A overhang. The new T:A pairing is explicit, as is the donor→acceptor backbone bond.
- `DirectionalCacc`: an incoming blunt strand begins5′CACC. A donor duplex has the complementary four-base5′ overhang on its opposite strand. In the manufacturer's drawing this bottom sequence reads `GTGG` **3′→5′**; its actual stored5′→3′ sequence is `GGTG`. Four original incoming base pairs are replaced. Their fate is explicit `Displacement` evidence with actual5′→3′ source calls and their remaining5′ continuation. A duplex-anchored continuation gives the **attached3′ flap** drawn in the standard directional diagram. An exact four-base complementary chain can instead be completely displaced as a **free four-mer**; an unpaired continuation is classified separately as an unpaired-tail extension. These complete-displacement/tail-extension cases are declared physical graph branches, not claims of primary-diagram kit geometry or efficiency. There is no covalent excision, fill-in or nucleotide loss.

Both donor pentamers `CCCTT` and `TCCTT` are recognized by this declared Vaccinia profile. This does not produce an activation from an uncharged substrate. Uncharacterized enzymes, noncanonical donor contexts and unsupported geometries remain explicit unsupported outcomes. Unknown terminal chemistry or ambiguous DNA remains unknown. Incompatible declared base pairs are invalid inputs. Reported experimental mismatch transfer does not justify assigning rates or declaring unsupported mismatches unable to react.

## API

```rust
use collagene_model::Sequence;
use collagene_topo::*;
let source = Substrates {
    strands: vec![
        Strand {role:"donor".into(),id:"top".into(),sequence:Sequence::dna("GACGCCCTT")?,
            five_prime:FivePrime::Phosphate,
            three_prime:ThreePrime::Phosphotyrosyl{enzyme:Enzyme::VacciniaTypeIBTyr274,activation_id:"declared-enzyme-1".into()}},
        Strand {role:"donor".into(),id:"opposite".into(),sequence:Sequence::dna("AGGGCGTC")?,
            five_prime:FivePrime::Phosphate,three_prime:ThreePrime::Hydroxyl},
        Strand {role:"insert".into(),id:"top".into(),sequence:Sequence::dna("CGTACGGA")?,
            five_prime:FivePrime::Hydroxyl,three_prime:ThreePrime::Hydroxyl},
        Strand {role:"insert".into(),id:"opposite".into(),sequence:Sequence::dna("CCGTACGA")?,
            five_prime:FivePrime::Hydroxyl,three_prime:ThreePrime::Hydroxyl},
    ],
    pairing:vec![
        PairSpan{strand_a:0,start_a:0,strand_b:1,start_b:0,length:8},
        PairSpan{strand_a:2,start_a:0,strand_b:3,start_b:0,length:7},
    ],
};
let steps = [Transfer{donor:0,acceptor:2,geometry:Geometry::ThreePrimeTA}];
let report = simulate(&source,&steps,&Limits::default())?;
assert_eq!(report.steps[0].outcome,Outcome::Applied);
assert_eq!(report.product.chains[0].sequence.as_str(),"GACGCCCTTCGTACGGA");
assert_eq!(report.product.nicks.len(),1);
audit_sources(&source,&report)?;
replay(&source,&steps,&Limits::default(),&report)?;
# Ok::<(),Box<dyn std::error::Error>>(())
```

`simulate` applies a bounded ordered list of specified elementary transfers. Every blocked/unknown/unsupported step remains evidence; any previously applied bonds remain explicit in the partial physical state. An unknown/unsupported step retains its incoming state, so later transfers describe that retained-state branch, not a prediction of the experimental distribution. `orientations` evaluates both source-strand orientations for one explicitly supplied vector/insert quartet. The caller supplies actual strands; the engine does not synthesize reverse complements to invent missing DNA. Other source strands/components are retained. This two-orientation search is never an exhaustive chemical reaction-mixture or kinetic search.

`ProductGraph` contains actual covalent `Chain`s, source placements, connected molecular components, open paired-backbone junctions with explicit nick chemistry, unpaired runs and released activation identities. `Displacement.kind` distinguishes duplex-anchored3′ flap, fully displaced chain and unpaired-tail continuation. Its continuation pairing is explicit; primary-diagram geometry requires that pairing to follow the incoming CACC duplex. `source_pairing`, `Nick` calls and step pair changes use **immutable source-strand coordinates**; `Chain.sources` translates them into output coordinates. Each source strand is present exactly once, without orientation reversal or nucleotide loss. `Component.paired_nucleotides` counts strand nucleotides, twice the number of base pairs. A component is covalently closed only when all of its chemical chains are circular; this does not assert supercoiling or biological success.

The whole-vector synthetic fixtures distinguish:

| Declared substrate geometry | Actual output | Unrepaired structures |
|---|---|---|
| Two charged TA vector ends plus uncharged3′A PCR ends | Two20-nt actual strands; all40 source nucleotides retained | Two nicks |
| Charged directional vector, blunt5′CACC insert | Actual strands24/28nt; all52 source nucleotides retained | One nick and attached4-nt3′ flap |
| Exact four-base CACC duplex (declared graph extension) | Joined donor/acceptor duplex plus separate unpaired GGTG strand | No excision; all four displaced nucleotides retained |
| One activated nick in a paired circular contour | One actual strand becomes circular | The opposite nick remains |

These are independent literal **synthetic mechanism fixtures**, not claimed sequences or activation states of a commercial vector. No host repair is assumed. `data/primary_sources.json` records the mechanism sources and hash of the official directional diagram consulted during development; no copyrighted PDF is bundled.

## Bounds and evidence

Borrowed source/string/count/span checks run before hashing, source-sized allocation or copying. Defaults permit100,000 total actual-strand nucleotides,32 source strands,1,024 pair spans and16 specified transfers. All limits apply together. Separate limits cover a conservative graph-work charge, graph-work allocation charge, cumulative orientation output bases and serialized report bytes. Hard caps are1,000,000 source nucleotides,64 strands,4,096 pair spans,64 steps and32MiB serialized reports. Workspace charges include per-base arrays and per-strand/step overhead; these are model resource budgets, not a measurement of process RSS or caller-owned input storage. A failed budget returns a typed `Error::Limit`; it cannot become negative evidence or uniqueness.

Hashes bind immutable declared content, including exact DNA, terminal states, identities, pairings, model version, steps and limits. They do not authenticate the person declaring activation. Public serialized reports are data. `audit_sources` independently verifies exact source identity and every output nucleotide/placement; it does not certify biochemical applicability. Only `replay` compares an imported report against fresh eligibility, bond, pairing and component reconstruction.

The engine does not model activation preparation/cleavage of a circular reference plasmid, enzyme rebinding/recycling, secondary cleavage, ligase, polymerase, covalent repair/excision of displaced flaps, mismatched transfer rates, transformation, host repair, selection or universal proprietary kit states. Those remain real integration/model work; reference-only substitution is not provided as a replacement. A future adapter must retain this full physical graph where a simpler duplex/end model cannot represent the flap.

## Primary evidence

- [Shuman1994, JBC269:32678–32684](https://pubmed.ncbi.nlm.nih.gov/7798275/): prepared covalent complexes and acceptor-dependent strand transfer.
- [Cheng and Shuman2000, NAR28:1893–1898](https://doi.org/10.1093/nar/28.9.1893), [full article](https://pmc.ncbi.nlm.nih.gov/articles/PMC103307/):3′T donor/3′A acceptor geometry, a nicked complementary strand, alternative insert orientations and condition-dependent mismatch transfer. This engine does not convert those observations into a kinetic prediction.
- [Manufacturer directional mechanism, printed p2](https://documents.thermofisher.com/TFS-Assets/LSG/manuals/pcdna3_1dv5histopo_man.pdf): Tyr274/3′-phosphotyrosyl chemistry,5′OH attack, CACC invasion and retained displaced strand in the depicted product.
- [Manufacturer enzyme mechanism](https://www.thermofisher.com/us/en/home/life-science/cloning/topo/topo-resources/the-technology-behind-topo-cloning.html): Vaccinia(Y)CCTT specificity and prepared activated vectors.
- [Manufacturer PCR-cloning technical support](https://www.thermofisher.com/uk/en/home/technical-resources/technical-reference-library/cloning-technical-support-center/pcr-cloning/pcr-cloning-getting-started.html):5′-phosphate inhibition and vector terminal-chemistry distinction.

No implementation code from these sources or other cloning applications was consulted.

## Explicit downstream ligase treatment

`nick_sealing` composes a fresh source-based TOPO reconstruction with an explicit
list of ideal templated DNA ligase nick seals. It requires actual 3′OH/5′phosphate
ends and continuous antiparallel pairing. This is separate from topoisomerase
`NickReseal`, which consumes a declared phosphotyrosyl activation and a 5′OH.

Every source call, pairing and displaced strand remains in the resulting graph.
The original transfer report is retained, and blocked/unknown/unsupported selected
treatments preserve their incoming state. No phosphorylation, flap processing,
kinetics, host repair or automatic target acceptance is inferred. Public entry
points reconstruct from original substrates and transfers; serialized product
reports are not accepted as physical starting state.

The application adapter, binary worker, MCP and downstream review controls are
described in [the source-bound contract](../../docs/TOPO_NICK_SEALING.md). Ordinary
paired-contour projection and tracked-duplex composition remain separate stages.

`processed_projection` now proves a source-seeded paired contour directly from
that freshly treated graph, reusing the private contour algorithm. It retains both
ancestral graphs and all other components. Strict requested-operation policy checks
original transfers and nick treatment; directional unpaired flaps remain withheld.
[Application bridge and remaining plan handoff](../../docs/TOPO_PROCESSED_PRODUCT.md).
