# Mouse C57BL/6 IGH reference admission — 2026-11-28

The existing antibody-numbering path now selects `references::Species` and an
explicit `v_mapping::ReferenceSelection { species, loci }`. The human reference
and scientific scoring policies are retained. Mouse C57BL/6 accepts heavy IGH
only, with IMGT numbering and IMGT regions; unsupported combinations refuse
instead of being relabeled or converted silently. No species classification,
accuracy calibration or other-species/strain claim is introduced.

The paired [OGRDB release-5 JSON](https://ogrdb.airr-community.org/api_v2/germline/set/10090.C57BL.sl.6%20IGH.C57BL.sl.6/5)
and [provider IMGT-gapped FASTA](https://ogrdb.airr-community.org/api_v2/germline/set/10090.C57BL.sl.6%20IGH.C57BL.sl.6/5/gapped)
are preserved under their [CC0 content dedication](https://ogrdb.airr-community.org/render_page/licensing_statement.html).
See [complete attribution, hashes and phase caveat](../../crates/antibody-numbering/data/ogrdb-mouse-c57bl6-igh-v5/NOTICE.md).
The compact corpus retains 95 functional V facts, four J records and 15 excluded
source records. Four atypical C23/W41/C104 observations across three functional V
alleles remain explicit. Mouse model training admits 92 canonical distinct
V-prefixes; human model admission continues its existing C23/C104 rule.

The only authorized maintainer emitter run produced a 123,204-byte runtime file,
SHA-256 `ea94e910b54fdf3318470c05ef0f382ad83ab7e189646bf14c0028b8285954cf`.
It used Python's standard library, checked both raw SHA-256 pins and all 95 V/four
J translations against the existing independent oracle. It did not call Rust,
external software or a biological runtime. Script and data hashes are in
`runtime-provenance.json`. Rust tests compare every compact Corpus fact with the
existing `Prepared` raw-pair projection; they never run the emitter as an oracle.

The compact corpus is cached once per species. Raw curation parsing is test or
explicit maintainer-feature code, not a dropdown operation. Corpus-v3, model,
query context and candidate reports bind species, chosen chain models, source
pins and policies. Current-format caller changes are direct; old corpus versions
and missing explicit selections are rejected.

J118–128 positions are an explicitly authored conditional IMGT complete-codon
rule, not provider labels. The original J frame/anchor values and observed
frame-1 versus AIRR-prose `(4-frame)%3` disagreement remain documented. The raw
curation's `JSegment` still has no labels. Runtime use does not normalize or hide
that difference.

New source fixtures cover exact Prepared/runtime equality, retained ledgers,
92/3 model admission, unchanged human subset models, unsupported combinations,
synthetic heavy V+ARDY+J routing, current bindings and Unicode/cross-species
corpus tampering. No compiler, test, browser or application HTTP checks were run
by this worker; parent owns the combined native/WASM verification batch.


Parent integration checkpoint: the shared native/WASM batch and 349 unique focused
cases pass. See [source and execution evidence](cloning-extensions-current-integration.json)
for the exact scope, initial catalogue-pin corrections and final metadata-only
wording edit. No graphical or application HTTP acceptance is claimed.
