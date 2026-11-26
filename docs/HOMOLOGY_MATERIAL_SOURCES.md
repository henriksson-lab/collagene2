# Primary references for a physical homology-assembly model

Research baseline: 26 November 2026. These sources describe chemical activities;
they do not establish deterministic resection lengths, experimental success or
an implementation of a proprietary formulation. No implementation code was used.

| Source | Supported interpretation | Boundary |
| --- | --- | --- |
| [Gibson et al., 2009, DOI](https://doi.org/10.1038/nmeth.1318), [paper PDF](https://static.igem.org/mediawiki/2013/d/d0/3.1-_Enzymatic_Assembly_of_DNA_Molecules_up_to_Several_Hundred_Kilobases.pdf), p.343 Fig.1 | T5-mediated removal from 5′ ends, complementary-tail pairing, polymerase gap filling and nick sealing by ligase | p.344/Fig.2 also describes linear products as further exonuclease substrates; unchanged linear outer ends cannot be assumed |
| [Pickering et al., NAR 27:730–735, 1999](https://academic.oup.com/nar/article/27/3/730/2902074), Results after Fig.2 | Direct mass-spectrometric evidence of T5 cleavage yielding 3′OH and 5′ phosphate | A defined overhanging hairpin establishes cleavage chemistry, not the distribution of assembly-reaction cut positions |
| [NEB Taq DNA Ligase](https://www.neb.com/en/products/m0208-taq-dna-ligase), Product Information | Sealing requires adjacent, correctly paired 3′OH and 5′ phosphate without a gap | No authority to seal blocked/unknown ends, 5′OH, gaps or an unsupported collision |
| [NEB Gibson mutagenesis note, Fall 2012](https://www.neb.com/en/-/media/nebus/files/expressions/expressions_fall12.pdf?hash=8CB799C5EC605C3D5F83CA7414396153&rev=5af05b5ad9eb41e48f5b87302b581dc4), [T5 technical reference](https://media.neb.com/m/bd28765a734dba9/original/us_neb_catalog_technicalreference_2025_26.pdf) | Phosphorylated primers are unnecessary for ordinary Gibson; T5 accepts phosphorylated and dephosphorylated substrates | Positive backbone cleavage exposes a new terminus; it does not phosphorylate an untouched original PCR end |
| [NEBuilder HiFi manual v6.0](https://www.neb.com/-/media/nebus/files/manuals/manuale2621_e5520.pdf?hash=02CF464149634283F5633CA1DE4D2FDF&rev=67d09a762b2c47ed8df3a0bd1cb59a1b), p.2/Fig.1 | Broad resection/pairing/filling/sealing model and use of products in subsequent PCR/RCA | Proprietary mismatch handling, concentrations and kinetics are not disclosed sufficiently for reproduction |
| [Takara In-Fusion planning guidance](https://help.takarabio.com/36168/kb/article/109886/how-should-i-plan-my-in-fusion-experiment) | A distinct 3′ exonuclease route creates 5′ overhangs; cellular rescue completes the construct | The guidance explicitly distinguishes this from covalent in-vitro assembly of linear molecules; do not inherit a sealed Gibson product |
| [Li and Elledge, 2007](https://www.nature.com/articles/nmeth1010), paper and supplementary protocols | SLIC uses exonuclease-generated intermediates and recombination/annealing, with distinct RecA/no-RecA protocols | This does not establish in-vitro Taq-ligase sealing; use a separate model |

A first independently implemented model should use explicit bounded resection
lengths on intact paired inputs, keep surviving strand mappings and chemistry,
assign a new 5′ phosphate only at an actual supported backbone cleavage, fill
only template-supported gaps from available 3′OH ends, and seal only eligible
nicks. Zero resection must preserve the original terminal group. An exact-H
resection fixture can isolate gapless homology pairing; asymmetric extra
resection needs independently calculated gap/fill fixtures.

Circular assemblies are the clearest initial material workflow because every
end participates in a declared junction. Linear workflows need an explicit
external-end policy and must retain the actual resulting overhangs rather than
silently restoring blunt outer ends. This initial order does not remove linear
assembly or broader method coverage from the product goal.

Required fixtures include ordinary 5′OH PCR circles, zero-cut chemistry, asymmetric
gaps, retained/removed/synthesized modifications, reversed strands, blocked ends,
excessive or conflicting resection, competing homologies, linear outer ends and
fresh material replay into PCR/digestion. Target letters must never select a
physical path or authorize unsupported repair. See [implementation scope](HOMOLOGY_MATERIAL_SCOPE.md).
