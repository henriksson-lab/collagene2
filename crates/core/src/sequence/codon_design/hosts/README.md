# Reference-CDS host preferences

Registered bounded derivation and compact pinned runtime catalogue for codon
design. Uniform, explicit and reference-host preferences are mutually exclusive.

`HostProfile::derive_mg1655()` freshly derives the captured reference dataset.
`derive_gzip(Dataset, Compressed)` and `derive_expanded(Dataset, Expanded)` accept
explicit bounded source sets for curation and fixtures. The sources are CDS FASTA,
genomic GFF3, assembly report and annotation-hash metadata. Parsing binds every
FASTA entry to one grouped GFF chromosome/locus/protein identity and exact strand/
coordinates, checks parent gene identity and containment, and requires every GFF
CDS group to be represented once. Distinct paralog loci remain distinct even when
protein identifiers or complete sequences are identical. This does not independently
check CDS bases against a genomic FASTA, which is not an input.

Counts exclude initiation and terminal termination. Start and stop observations
are separate. Admission requires complete single-part, nonpseudogene,
nonexceptional code-11 CDS in phase zero, a recognized start, a definite terminal
stop, no internal stop and only concrete DNA. Biological exclusions are explicit
and mutually exclusive; malformed syntax, mismatched locus bindings and exhausted
resources reject the whole derivation rather than returning a partial table.
Recognized initiation and codon assignments come from the canonical translator;
no table copy or visibility change is needed.

Combined compressed input is capped at 20 MiB and combined expansion at 32 MiB,
with 20,000 CDS/gene identities, 100,000 GFF rows, 4 KiB lines, 64 pieces per grouped
CDS and 300,000 bases per CDS. The parser has no external process, network or
filesystem operation. Curation uses this pinned Rust backend:

```toml
flate2 = { version = "=1.1.11", default-features = false, features = ["zlib-rs"] }
```

`usage()` returns observed counts, integer millionths overall and within a
synonymous family, and explicit unknown-family `None`. `preferences()` produces
64 integer weights for the current additive codon-design objective. Observed
within-family maxima get 1,000,000; zero observations and unobserved families get
1; terminal-stop preferences stay uniform. These are arithmetic preferences, not
CAI, expression efficiency, host fitness or calibrated synthesis-success claims.

Profile identity binds dataset/assembly metadata, raw and expanded BLAKE3 hashes,
exact annotation timestamp/hash text, derivation engine/policy, counts/exclusions
and sorted admitted-locus identity hash. Structural validation does not authenticate
an arbitrary imported profile's claimed organism or count completeness. The curated catalogue pins the freshly derived identity and compact result, and
binds selected identity/provenance to expanded optimizer weights and retained
reports.

The runtime catalogue loads only the compact `profile.json`; it never derives the
full genome on UI selection. The maintainer `derive_codon_host` example freshly
emits the profile. Compressed originals remain curation/test evidence. The separate
stdlib Python oracle in `data/mg1655/independent_oracle.py` tallies source counts
without the Rust parser or translator; all 192 internal/start/terminal values match.
It is a maintainer verification tool, never an application runtime dependency.

Private tests (`sequence::codon_design::hosts::tests`) cover independent tiny count
oracles and start/stop policy, distinct paralog loci, exclusion accounting,
locus/coordinate/code mismatches, duplicate/missing records, parser/gzip bounds,
identity conservation and fresh captured-MG1655 derivation. The integrated native/WASM checks and all nine host cases pass. The connected
61-case scope also covers design/check, retained evidence, MCP and UI state.
