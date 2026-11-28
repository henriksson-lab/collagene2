# Mouse C57BL/6 IGH curation checkpoint

Source-only checkpoint: 2026-11-28. This module is intentionally **unregistered**;
no existing Cargo/module/schema/UI/enum is changed and no compiler or Rust test
has run. Root owns the next connected species-aware registration and `.gitattributes`
update preserving raw data line endings.

The owner is [references/mouse_c57bl6](../../crates/antibody-numbering/src/references/mouse_c57bl6/mod.rs).
`Prepared::new(json, gapped)` accepts only the bounded exact pinned OGRDB v5
byte pair, then verifies unique full-coding-DNA correspondence, source identities,
functional admission, whole V-prefix triplet gaps, source CDR intervals and J
anchor/extent facts. `dataset()`/`into_dataset()` expose serializable reference
facts; `identity()` binds bounded derived serialization and raw BLAKE3 provenance.
This is offline curation, not another numbering/detection engine or dropdown-time
runtime parser. A small prepared runtime dataset can be emitted once during the
integration batch.

Original inputs, independently captured hashes/timestamps/license proof and
notices are in [the dedicated data directory](../../crates/antibody-numbering/data/ogrdb-mouse-c57bl6-igh-v5/NOTICE.md).
The raw JSON is 254,707 bytes; provider-gapped FASTA is 35,235 bytes. Both exact
SHA-256 values were checked again after copying. BLAKE3 and final derived-identity
literals await execution of the connected owner; the source computes them without
new dependencies. No hash of curated assertions is presented as an upstream hash.

The independent worker oracle decoded the same exact sources and NCBI's pinned
standard genetic-code entry without using this Rust parser/table. Its full report
and script remain in `/tmp/collagene2-antibody-reference-research-2026-11-28/`.
A compact expected-translation fixture is retained with its script/source hashes.
Observed source facts:

- 114 records, all one-to-one unique coding-DNA joins; differing FASTA/JSON names
  stay separate provenance fields.
- 95 functional V, seven nonfunctional V, eight functional D and four functional J.
  D and nonfunctional V are excluded from V/J model facts, never discarded from
  the source/exclusion ledger.
- All 95 functional V have usable first-104 whole-triplet gapped columns and
  matching JSON CDR nucleotide intervals. Only 92 have canonical C23/W41/C104.
  IGHV-ZPQB has S23; IGHV-IQE5 C41; IGHV-OSZZ R41/S104. All provider-functional
  references remain retained with an explicit four-observation atypical ledger.
- Full gapped V lengths need not be divisible by three; unlabeled germline tails
  are outside the translated 104-column prefix.
- All four J references yield 11 complete FR4 codons and trailing G when anchored
  to the **actual pinned** `j_cdr3_end` field. Their phase disagrees with literal
  AIRR prose; raw frame, observed phase and documented phase remain distinct.
  No J numbering labels are manufactured.
- All retained V-prefix protein strings and J-FR4 strings are distinct in this
  capture. The owner nevertheless retains allele identity and records duplicate
  protein groups without collapsing them when present; a synthetic tie fixture
  exercises that behavior.

Eight authored nested tests cover inventory, independent full translations,
literal anchors, J phase discrepancy, raw bytes/line endings/bounds, ambiguous
DNA joins and duplicate source use, duplicate labels, partial-codon gaps, altered
source coordinates/frame and protein-tie retention. These are new unregistered
tests, not moved Cargo coverage. Their execution remains pending the integration
batch; formatting and repository layout/whitespace checks are separate static checks.

Public `VReference` facts are reused where their invariants fit. Their mouse
labels explicitly mean ordinals of **provider-declared IMGT-gapped columns**,
not a nonexistent JSON label array. A cohesive `JSegment` holds unnumbered
provider/FR4 facts because the human `JReference` shape requires authored 118+
labels. Future integration must explicitly choose and validate that general
numbering rule, not pretend the provider supplied it. The existing human corpus
and its scientific behavior are unchanged.

The broader [primary-source expansion audit](antibody-reference-expansion-current.md)
records mouse/light-chain, alpaca, TCR and further species paths. IMGT's current
CC BY4.0 policy conflicts with older download-specific notices; those raw files
remain research-only in `/tmp`. This checkpoint redistributes only the separately
CC0 OGRDB inputs and independently derived factual fixture, with no provider
implementation code or external scientific runtime.
