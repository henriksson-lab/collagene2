# Inspecting expected assembly junctions

`core::assembly_junctions::inspect_plan(plan, independent_target, options)` runs
fresh verification once and uses the physical trace from that same reconstruction.
It never accepts an imported report or product metadata as chemical authority.
The original plan and independent target remain distinct bindings. The returned
context hash covers the complete original plan, actual target and options.

Initial supported final operations are ordered `Ligate`, `Overlap` and
`ProfiledOverlap`. Other final operations return their fresh verification and,
when available, exact final product receipt, with junction inspection explicitly
incomplete. A final product exists only if every plan step reconstructed; the
last successful preparation is never substituted for a failed final assembly.
A wrong target can still leave an actual reconstructed final product available
for inspection while verification remains invalid.

## Coordinates and actual calls

Ranges are zero-based half-open coordinates of the freshly reconstructed product.
Ligation boundaries sum the **oriented actual duplex top lengths**. Reversing a
staggered fragment can change its top length, so raw reference lengths cannot
replace this calculation. Adjacent end summaries and retained final-product
nicks come from the fresh physical trace. Endpoint descriptions are not an
independent claim of experimental bond formation.

For overlap assembly, each shared homology interval is derived from the actual
resolved prepared source lengths and explicit overlap. The displayed boundary is
its **end**, where downstream unique sequence begins. For example, `AAAACCCC`
followed by `CCCCGGGG`, overlap four, gives homology `[4,8)` and boundary 8.
Circular closure uses boundary 0 and the retained prefix homology `[0,k)`.
Homology exceeding one final-product turn is withheld instead of repeatedly
cycling calls or inventing a target-guided placement. Overlap windows carry DNA
correspondence only; missing physical state does not become material authority.

Windows contain at most the requested number of actual calls on either side of
the boundary. Linear ends clip the available flank. Circular windows split into
explicit product ranges across the origin; short products contribute each call
at most once. A short circular product is never repeated to fill a display.
Clipping at these known boundaries is a display condition, not missing scientific
information: it sets `Window.clipped` without making `Report.complete` false.
Omitted junctions, unsupported reconstruction and withheld requested CDS evidence
remain distinct reasons for an incomplete inspection.
`Window.boundary_offset` locates the boundary in its forward-reference DNA string.
Full homology ranges can exceed the bounded DNA window; the window is not claimed
to contain the entire homology in that case.

## Frame and codon interpretation

**Frame choices are window-relative, not product-origin frames.** Each
`Frame {reverse, frame}` applies the existing translation engine to the forward
window or its reverse complement, omitting 0, 1 or 2 calls at that oriented
window's beginning. The boundary offset is reversed with the window. An optional
`crossing_codon` is the complete translated codon strictly spanning that boundary;
`window_start` is its start offset in the frame's oriented window. A boundary
exactly between codons does not produce a falsely crossing codon.

The genetic code is explicit (existing audited codes 1, 2 or 11). Translation uses
ordinary codons, retains stops/ambiguity, and never infers initiation, reading-frame
intent or fusion function from labels. Changing the flank length can change local
frame alignment; the report states this coordinate system instead of silently
calling those translations product-origin frames.

## CDS evidence

Current projected CDS annotations intersecting a flank are reported by actual
feature index. `spans_boundary` initially describes coverage of both adjacent
reference calls, without assuming an exact biological connection. Exact uniform
strand/join geometry must additionally place those calls consecutively in the
CDS biological coordinate order before any coding offset is used.

The existing feature-directed CDS engine validates phase, genetic code,
exceptions and review state. Ordinary genes up to 30,000 coding bases are
supported, with an aggregate one-million-base translation budget. Only a bounded
local codon DNA/peptide window is returned, not the full gene/protein. Its
`first_codon`, `junction_coding_offset`, `codon_start`, strand and actual genetic
code identify the coordinate interpretation. CDS initiation/terminal-stop rules
remain those of the existing engine. The frame options for generic flank windows
do not override a CDS's declared translation table.

Review-marked, fuzzy, ordered, mixed-strand, repeated-coordinate, unsupported-phase
or exception-bearing annotations retain explicit warnings rather than a guessed
fusion frame. A nearby CDS that does not span the boundary is distinguished from
one that does. Oversized genes, exhausted aggregate work, unplaced annotations or
omitted candidates make requested inspection incomplete. A complete local reading
frame still does not establish fusion function or an experimental outcome.

## Bounds and readonly binding

`Options` defaults to 30-base flanks, six explicit window-relative strand/frame
choices, code 1, at most 128 junctions and requested CDS evidence. Flanks are
1–150 bases; frames are unique; no more than 128 CDS candidates are returned over
the whole report. Nick evidence copies are precharged at an aggregate 100,000
limit. Original input-plus-target DNA is capped at 8 MiB; existing borrowed model,
verifier and scientific-engine bounds remain in force. Canonical request/output
is capped at 32 MiB. `source_preflight(plan,target)` guards originals before host
copies; `inspection_revision` also validates options.

`Report::validate`/`validate_binding` perform readonly shape, interval, DNA-source
and context checks without rerunning verification or translation. A decoded
report grants no authority to edit a document, recover chemistry or promote a
candidate. Existing workspace and plan schemas remain unchanged; only the new
worker analysis envelope carries this report.
