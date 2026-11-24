# GTF annotation example

Open `variants-demo.gb`, then open `annotations-demo.gtf` and choose **GTF** in
annotation import. Map its sequence name explicitly to `chr_demo`, review the
rows and select the desired proposals before applying them.

The seven-base reference is `AACCGTT`. The reverse-strand exon is the zero-based
half-open interval `[1,5)` (`ACCG` on the forward reference; `CGGT` in the reported
feature orientation). The separately reported reverse CDS is `[1,4)` and retains
frame 0. These synthetic rows do not establish a functional CDS or transcript.
Repeated `tag` attributes stay separate and the numeric `exon_number` is retained
as source text. No transcript joining is inferred.

GTF export normalizes quoting while preserving the retained row meaning. Download
the original source for exact lexical bytes; save the native workspace to retain
source evidence alongside edits. Editing coordinates or qualifiers requires
review before source-backed normalized export can be trusted.
