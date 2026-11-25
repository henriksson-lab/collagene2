# Selected translation and linked protein outputs

Translation now has an explicit whole-reference or current-selection scope.
An exact selected feature retains its biological join/complement order; an
origin-spanning search/ORF selection retains every part. Frame skip occurs at the
oriented 5′ end. Global reverse orientation is separate from the selected location's
own biological orientation. Exact repeated parts are translated in their declared
order rather than silently merged. Fuzzy, ordered-uncertain and unresolved region
selections are rejected instead of inventing exact codon boundaries.

The selected-CDS action uses the existing native CDS translator and the actual
feature's code, phase and qualifiers. It does not use a stored protein string as
the answer. Its warnings and stored-translation comparison remain visible. Literal
frame translation keeps every stop symbol; native CDS translation omits one
definite terminal stop. Neither establishes expression or experimental success.

Protein review pages the complete result in 600-symbol pages and exports full
protein FASTA. Named output actions ask a worker to reconstruct from the original
source and request again. Source document/event, workspace epoch, selected request
and name must still match. The host independently validates source reconstruction
before installing a document or annotation.

Open protein creates a clean linear protein document, with the exact retained
nucleotide event as its computational parent. It inherits no DNA chemistry,
sequencing descriptors or original GenBank header authority. Its metadata retains
source identity, original selection, code/frame/initiation options, output name,
context hashes and warnings. The origin panel shows the codon policy and warnings.
Its source action checks ownership, revision and actual ancestry, then opens the
historical nucleotide version read-only without moving its current head. A
metadata-only imported link cannot stand in for retained ancestry. This origin
link does not verify later protein edits.

Retain CDS creates a proposal at the effective selected location with genetic
code and codon-start qualifiers. It deliberately requires review and carries no
invented partial flag or stored translation assertion. Notes retain the original
initiation policy: biological CDS initiation and terminal-stop interpretation can
differ from an arbitrary literal frame translation. The existing feature editor
can record human review. Adding the proposal is one undoable document edit and
invalidates the previous source-bound result.

ORF review now pages all returned candidates in groups of 100. Select ORF
highlights its complete compound footprint, retains its original strand/order for
subsequent selected translation, and never reduces a circular ORF to its first
span. Different source or selection geometry invalidates that directed location.
Source/epoch guards also cover ORF annotation and FASTA export. Search truncation
and the native finite model remain visible; paging does not imply exhaustive ORFs.

The native `translation_reference` layer composes existing codon algorithms.
Worker schema 54 appends translation/reference-output jobs and responses under the
ordinary 32 MiB binary message limit. Native fixtures independently cover frame,
compound/reverse/circular order, RNA, ambiguous calls, stops/trailing bases, CDS
policy, clean outputs, stale bindings, tampering, historical ancestry and binary
version/budget gates. Local browser/performance acceptance remains deferred until
the full in-scope implementation is complete.
