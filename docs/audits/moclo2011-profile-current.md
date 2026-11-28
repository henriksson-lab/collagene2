# F37 original MoClo Level-0 profile

2026-11-28: integrated profile in `crates/assembly-standards`. Native/WASM
checks and all six new profile fixtures pass, alongside the existing RFC10 and
RFC1000 cases. [Combined source/execution evidence](verifier-and-design-current-integration.json).
No browser or application HTTP checks run.

The primary [Weber et al. 2011 article](https://doi.org/10.1371/journal.pone.0016765)
and Figure1B/C establish these two paths:

| Role | Left fusion | Right fusion |
| --- | --- | --- |
| promoter | GGAG | TACT |
| 5′ UTR | TACT | AATG |
| signal peptide | AATG | AGGT |
| coding module after signal peptide | AGGT | GCTT |
| cytosolic coding module, replacing the preceding two | AATG | GCTT |
| terminator | GCTT | CGCT |

The article's permissions grant Creative Commons Attribution reuse; its captured
XML omits a license version. Source attribution is retained in the profile and
NOTICE. Publisher XML (75,479 bytes) SHA256:
`cbc7be2584380029a96188834bff83a13eec15c0d35e83dc830101ad25ade0ad`.
Figure1 PNG (104,623 bytes), actually visually inspected, SHA256:
`1b6b0313fe229271d33b0484fdcf19cec8ee3562bb93e4bc5f86df5e707240c4`.
Exact URLs and capture facts are in `data/moclo-l0-2011/source-provenance.json`;
raw assets remain outside the repository. Our factual curation/code are MIT;
the upstream work is not relabeled MIT.

BpiI entry and BsaI release are documented in the article's module-generation
section. The [official recognition/isoschizomer table](https://www.neb.com/en-de/tools-and-resources/selection-charts/isoschizomers)
crosschecks BpiI/BbsI `GAAGAC(2/6)` and Esp3I/BsmBI `CGTCTC(1/5)`.
The current catalogue additionally binds BsaI `GGTCTC(1/5)`. BsaI/BpiI motifs
are required absent under the selected entry/release body policy; Esp3I is
advisory for later hierarchy. The profile does not claim an Esp3I-containing
body necessarily prevents an isolated BsaI-only assembly.

`MOCLO_2011_L0_ID` extends the existing `Profile` registry directly. The existing
finite grammar represents both branches with four-to-five parts and each role
used at most once. No alternate representation, new schema, adapter or algorithm
is introduced. All expected bodies have empty additions. AATG/AGGT overlap facts
do not manufacture bases, coding frames, signal processing or function. The
result means conformance to declared module syntax and selected body-site rules,
not vector or complete transcription-unit validation.

Six independently authored fixtures cover the literal table, both complete
paths and successor hints, wrong/cross-branch/omitted order, literal forward and
reverse motifs with required/advisory severity, ambiguity/work caps/unclaimed
flanks, and reserved identity/catalogue/retained-report tampering. The two current
RFC profiles retain their original definitions and distinct assessment identities.

Remaining grammar gaps are explicit: no combined promoter/UTR modules, Level1/2
destination/linker variants, physical-end admission, translation/frame policy or
backbone conformity. Those need separate source-curated scope. The upstream YTK
p4a→p4b/CGAC branch is not approximated or grafted onto a presumed GCTG rejoin.
