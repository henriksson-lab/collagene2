# GenBank RNA encoding and exact source evidence

GenBank import normalizes T to native U only when the record declares an RNA
molecule type, including mRNA. Other alphabets keep their own validation rules.
Features, compound/reverse coordinates, qualifiers and source record boundaries
remain unchanged. The original record is retained byte for byte.

This follows the INSDC RNA flat-file representation documented by
[DDBJ](https://www.ddbj.nig.ac.jp/ddbj/sequence-e.html) and
[NCBI](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/data-processing/taxonomy-processing/genetic-codes/).
The application keeps its RNA alphabet as U; GenBank export writes T in ORIGIN.
Legacy U-encoded RNA records remain readable and export in that canonical form.
Import/export reports disclose the encoding projection.

Source comparison and citation edit replay use the original record's declared
alphabet and normalized biological sequence. Encoding normalization alone leaves
unchanged accession/version and citation scope intact. Actual edits still require
the existing source-bound replay recipe, remap citation positions and clear stale
active database identity. Changing the original bytes or supplying an unrelated
baseline cannot bypass that replay. RNA citation notices identify edited RNA.

Independent literal fixtures cover T-encoded mRNA, reverse compound locations,
synthetic native U export, unchanged source identity, citation remapping after
insertion, changed original bytes, unlogged edits, per-record RNA/DNA distinction
and older U-encoded records. They use no external scientific program or network
request. Browser acceptance remains deferred.
