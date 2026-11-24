# Synthetic declared-allele example

Open `variants-demo.gb`, then open `variants-demo.vcf` separately. Select the
first record, explicitly choose `Reference_DNA`, enter chromosome mapping
`chr_demo`, and choose ALT index **1**. The genotype column `0|1` is file evidence
and never selects an edit automatically.

The literal REF `CC` occupies reference positions3–4 (zero-based interval2..4).
Fresh checking must compare those actual DNA calls before enabling a preview.
The selected ALT `CTC` produces `AACTCGTT`; the original `AACCGTT` reference stays
unchanged. The independently supplied second GenBank record contains the same
literal desired design. Review the intersected source landmark annotation.

Opening the preview retains the original variant bytes and exact reference-event
parent, and marks the new document as a chemistry-unmapped reference design.
A matching file allele does not establish that the designed DNA exists as cloning
material. Use the independent design target in a separately reconstructed cloning
plan rather than silently treating it as a prepared substrate.

The second record has two alternatives and a missing/declared genotype column;
choose an ALT explicitly to inspect it. The symbolic `<DEL>` row remains source
evidence and must not manufacture a deletion from its END field. Original and
normalized exports remain separate. The example is synthetic coordinate data,
not a measured genotype. Browser validation remains deferred.

The design retains its explicit record/ALT recipe as unverified metadata. Fresh replay requires the original variant bytes and exact reference revision; a changed design is compared with a newly reconstructed result rather than trusting its saved check hash.
