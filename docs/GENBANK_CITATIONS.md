# GenBank citations after sequence edits

A literature reference can remain useful source evidence without describing a
precise interval on the edited sequence. GenBank permits annotation-site references,
an omitted parenthetical scope and prose within the reference description.
[NCBI format change](https://www.ncbi.nlm.nih.gov/genbank/-release/123/),
[current reference format](https://www.ncbi.nlm.nih.gov/genbank/release/current/).
These bibliography forms do not block an otherwise supported native DNA/RNA edit.

## One edit log and a fresh proof

`model::genbank_edit` owns the current serialized `EditLog`, its `SequenceEdit`
values, framing/resource admission and runtime `VerifiedLog`. Molecule and its
meaningful edit implementation remain together. Old root exports and forwarding
functions are removed; callers use the canonical owner directly.

The edit log retains its current seven fields, field order, metadata key and
version 1. It is untrusted data, not a citation authority. Unsupported versions
are rejected. A runtime verified log cannot be deserialized from saved claims.
Export independently parses the selected original record and normalizes its
DNA/RNA before checking the original bytes, source index, baseline sequence,
current sequence and complete edit replay. The retained source is a single exact
record, so its selected index is zero. A changed index cannot gain authority by
changing the ledger, metadata and supplied index together.

Malformed or stale logs, missing/changed originals, wrong baselines, wrong current
DNA, invalid intervals and exceeded replay/input limits remain errors. Neither a
metadata boolean nor a generated remark can waive these checks. Edited records
receive borrowed framing/reference admission before rich-header allocation.

## Citation outcomes

| Source description or resulting coverage | Exported outcome |
| --- | --- |
| Fully consumed numeric `bases` grammar with surviving original bases | Mapped scope containing only retained source positions |
| `n (sites)` | Source-only bibliography, bare `REFERENCE n` |
| Bare `REFERENCE n` | Source-only bibliography without an invented whole-record scope |
| Other bounded prose or an unrecognized complete description | Source-only bibliography; no parsed-prefix mapping |
| All originally cited bases replaced/deleted, even if new DNA remains | Source-only bibliography; later inserted DNA cannot revive the scope |
| Valid retained coverage exceeds the supported output descriptor | Source-only bibliography with an explicit representation-limit reason |
| Invalid numeric bounds, reference identity/framing, input limits or replay | Error, not a source-only fallback |

Known numeric grammar is recognized completely before numeric interpretation.
Repeated `bases` after semicolons is supported. Citation numbers use `u64`, so their
identity does not depend on native versus WASM pointer width. Sequence positions
are converted to native indices only after checking the sequence length.

Coverage follows a monotone map of original source bases through the edit log.
Inserted and replacement bases have no source mapping. Insertions can split a
cited interval; removing the inserted gap can join its surviving coverage again.
This differs from the annotation endpoint envelopes retained by native editing.
A geometrical envelope containing only replacement DNA is not citation coverage.

Source-only references retain their numbers, authors, title, journal and database
identifiers. The original description, reason and explicit lack of destination
coverage appear in `REMARK`; `/citation=[n]` links remain stable ancestral links.
Mapped references also state that citation support for edited DNA/RNA is
unverified. Source records remain exact in native storage. Repeated exports and
native reopening do not multiply generated review prose. Remarks are presentation
and provenance, never an instruction to restore coverage authority.

## Limits and evidence

Input limits retain 1,024 references, 4,096 description bytes, 256 numeric spans
per reference, 4,096 total numeric spans and the bounded scope-work policy. The
edit log retains 1,024 steps, 512 KiB of inserted sequence, a 1 MiB serialized limit
and 64 million replay units. Original records remain bounded to 32 MiB. Output
representation limits withhold the current scope rather than assigning new DNA
coverage; they do not hide malformed-input or replay failures.

No native binary layout changes or second ledger encoding are introduced. The
scientific engine and affected native-edit engine stamps advance to invalidate
stale contexts. Native/GenBank/RNA round trips, source-only and mixed references,
replacement/deletion, complete-consumption parsing, portable identifiers, strict
proof failures and bounds have independent fixtures. Current integration evidence
is recorded in the [delivery audit](DELIVERY_AUDIT.md). Installed-platform and
visual acceptance remain separate requirements.
