# External chain-identity assessment panel

This is a development assessment dataset from the public wwPDB archive through
RCSB PDB, retrieved 27 November 2026. The 30 distinct included protein sequences
comprise 12 antibody heavy/light chains and 18 non-antibody controls. Deposited
molecule identities and associated structure publications define these **chain-
level** labels. They do not provide reference residue numbering, exact CDR
boundaries or statistical population-level calibration.

Archive data is dedicated under CC0-1.0; Collagene's generator and assessment code
are MIT licensed. Integrated IMGT, SAbDab, Pfam, SCOP and other external annotation
records are not copied, used as numbering truth or assumed to share archive terms.
[RCSB archive usage/citation policy](https://www.rcsb.org/pages/policies).

The source manifest pins 48 entry/entity API responses by exact SHA-256 and URL.
The derived panel retains canonical deposited sequence, molecule description,
entity identity and source provenance. All canonical sequences are used whole,
including constant domains and expression tags; no candidate-guided crop is used.
Identical sequences would be grouped while retaining every source. Four source
entities are excluded explicitly: three peptides shorter than 50 residues and a
GFP record whose canonical sequence contains unresolved X. Excluded sequences and
reasons remain in the panel instead of inflating the negative-control count.

| Entries | Assessment roles |
| --- | --- |
| [1HZH](https://www.rcsb.org/structure/1HZH) | Intact b12 heavy/light chains |
| [1N8Z](https://www.rcsb.org/structure/1N8Z) | Herceptin Fab chains; HER2 control |
| [1BJ1](https://www.rcsb.org/structure/1BJ1) | Humanized Fab chains; VEGF control |
| [5GGQ](https://www.rcsb.org/structure/5GGQ), [5GGS](https://www.rcsb.org/structure/5GGS), [5GGT](https://www.rcsb.org/structure/5GGT) | Therapeutic Fab chains; PD-1/PD-L1 controls |
| [1AO7](https://www.rcsb.org/structure/1AO7), [1D9K](https://www.rcsb.org/structure/1D9K) | TCR alpha/beta and MHC controls |
| [1HSA](https://www.rcsb.org/structure/1HSA) | MHC and beta-2 microglobulin controls |
| [1CD8](https://www.rcsb.org/structure/1CD8), [1CDH](https://www.rcsb.org/structure/1CDH) | CD8/CD4 immunoglobulin-fold controls |
| [1TEN](https://www.rcsb.org/structure/1TEN), [1LYZ](https://www.rcsb.org/structure/1LYZ) | Tenascin/lysozyme controls |
| [1EMA](https://www.rcsb.org/structure/1EMA) | Explicitly excluded unresolved GFP sequence |

Related or homologous controls are not independent statistical samples. Antibody
chains include human and humanized constructs; archived source-organism labels
are not reinterpreted as a species-classification target. Truth labels were fixed
before inspecting algorithm output. This small panel cannot establish broad
specificity, sensitivity, species/chain assignment accuracy or numbering accuracy.

Rebuild from the exact original local responses named in `sources.json`:

```sh
python3 scripts/build_antibody_assessment.py ORIGINAL_RESPONSES_DIRECTORY \
  --output crates/antibody-numbering/data/rcsb-assessment/panel.json
```

The generator verifies source hashes; it performs no network or scientific-tool
calls. The native `assess_rcsb` example emits reproducible summaries, exact
engine/options/limit identities, dataset revision, every missed positive and
non-antibody candidate, and source coordinates. It is a development tool, not
an installed external runtime dependency or a production classification profile.


Run the embedded panel from the repository root with the normal repository Cargo
cache/profile settings:

```sh
cargo run -p collagene-antibody-numbering --example assess_rcsb --offline --locked > assessment.json
```

The sequence identities and truth labels stay fixed across the three explicit
FR4 mismatch budgets. This panel supplies no lambda-antibody positive, no complete
scheme assignment and no residue-level numbering reference; those remain required
independent validation work.
