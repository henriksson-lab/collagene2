# Pinned mouse C57BL/6 IGH curation inputs

Source: Open Germline Receptor Database (OGRDB), set `OGRDB:G00004`,
`OGRDB:C57BL/6 IGH.5`, release 5 (2024-10-15), author Katherine JL Jackson.
Species is *Mus musculus*, NCBITAXON:10090, **C57BL/6 strain**, IGH only.

OGRDB's [content dedication](https://ogrdb.airr-community.org/render_page/licensing_statement.html)
is CC0-1.0, separate from this repository's MIT software license and from
OGRDB's EUPL website implementation. No OGRDB implementation code is included.
The original license-page capture is retained in the task research directory;
its byte count, SHA-256 and timestamp appear in `capture.json`.

Raw files are copied **byte-for-byte without newline normalization**:

- `source.json`: [exact API v5 JSON](https://ogrdb.airr-community.org/api_v2/germline/set/10090.C57BL.sl.6%20IGH.C57BL.sl.6/5),
  254,707 bytes, SHA-256 `a557dc47e4946b015fcb5414c29c275783c464848fd13b8f4a0386c4d6e686e8`.
- `gapped.fasta`: [same release's provider-gapped representation](https://ogrdb.airr-community.org/api_v2/germline/set/10090.C57BL.sl.6%20IGH.C57BL.sl.6/5/gapped),
  35,235 bytes, SHA-256 `27403baa82141e88a065ddedbc9fc6b4bc4594fc9ecf64d8c4056e98faf4d160`.

The [provider API guide](https://github.com/airr-community/ogrdb/blob/master/schema/ogrdb_api_v2_guide.md)
describes `gapped` as IMGT-gapped FASTA. The JSON alone does **not** contain
aligned V strings or position-label arrays in this release. Curation joins
the paired files by exact full ungapped coding DNA, not their different labels.
All 114 source records have a unique match. Individual provider allele identities
remain distinct even if later protein models tie.

The maintainer curation source retains 95 provider-functional V-prefix facts
and four functional J-segment facts. Eight D records and seven provider
nonfunctional V records remain in the original source and explicit exclusion
ledger. Provider functionality does not certify a canonical anchor or numbering
model: four noncanonical C23/W41/C104 observations occur across three functional
V references. No residue is repaired and no reference is silently reclassified.
Only the first 104 gapped triplet columns are translated for V; full V strings
can end with incomplete codons and are not required to be divisible by three.

J anchor/frame observations are release-specific: raw `j_cdr3_end-1` points to
W with phase `j_codon_frame-1`. The literal AIRR prose interpretation
`(4-frame)%3` differs for these four records. Both values, original fields,
coding DNA, 11 complete WGxG-containing FR4 codons and the trailing G are retained;
this source does not supply J position labels. The curation owner does not invent
IMGT118 labels or normalize that discrepancy.
[AIRR germline field documentation](https://docs.airr-community.org/en/stable/datarep/germline.html).

`independent-oracle.json` retains expected reference translations derived by an
independent standard-library parser of these inputs and the authoritative NCBI
`gc.prt` standard-code entry (SHA-256
`2aecbdd09ecc5e0d233c105090a593737298ef93d88f6bc0c1bef26df6051a04`).
It does not call this Rust curation or copy its codon table. The original oracle
script/report remain under `/tmp/collagene2-antibody-reference-research-2026-11-28/`;
the fixture pins the script hash and primary format URLs. The source-authored Rust comparison tests are part of the connected species-aware batch.

`Prepared::new` requires this exact pinned byte pair and emits raw BLAKE3 values
and a bounded derived dataset identity using the existing Rust dependency.
Runtime compact facts and maintainer provenance are described below; raw BLAKE3
values are computed by Rust rather than fabricated by the Python emitter.
The paired raw sources are development curation inputs, not data to parse on
every GUI dropdown change. No classifier, species probability, experimental
fitness, VHH subtype, TCR or non-IMGT numbering claim is made.

## Current runtime admission

`references.json` is a compact current Corpus-v3 projection: 95 functional V facts
and four J records. `runtime-provenance.json` pins its SHA-256, emitter SHA-256 and
both source hashes. The maintainer-only stdlib `emit_runtime.py` verified those
raw hashes and the independent 95-V/four-J translation fixture before emitting.
It is never invoked by the application or by Rust tests as an authority. The
Rust `Prepared::runtime_corpus` test independently reconstructs the same complete
Corpus projection from the raw pair. Software remains MIT; these biological data
remain CC0. Corpus/model identity is computed with the existing Rust BLAKE3 owner.

Application loading caches this compact corpus once. The paired raw parsers are
compiled only for tests or explicit `reference-curation` maintainer builds.
`Species::MouseC57Bl6` permits `Locus::Heavy` only. The mouse model omits all three
functional references with atypical C23/W41/C104 residues, retaining 92 distinct
V-prefix proteins; it does not repair or remove their source facts. Human model
admission remains its existing C23/C104 policy and is separately identified.

The runtime J record adds **authored conditional** sequential IMGT118–128 labels
for the 11 observed complete FR4 codons. Their origin is the existing
[IMGT V-domain numbering chart](https://www.imgt.org/IMGTScientificChart/Numbering/IMGTIGVLsuperfamily.html),
not a provider J-label field or a generic interpretation of `j_codon_frame`.
Original frame/anchor fields, the contrary AIRR prose phase and the raw source
remain available. The curation `JSegment` continues to have no numbering labels.
Only IMGT numbering with IMGT regions is admitted for the mouse reference model;
Kabat/Chothia/North, light chains, other strains/species, confidence calibration
and experimental performance remain outside this admission.
