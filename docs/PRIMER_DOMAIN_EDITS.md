# Primer endpoints and tail-boundary editing

The shared pair editor now offers an explicit preview for source-forward endpoint
changes or the boundary between a synthetic 5′ tail and annealing domain. This
works in both ordinary PCR and mutagenesis editing. It changes editable oligo
text/reference comparison only; it does not search binding sites, infer synthesis,
validate specificity or claim amplification.

Endpoint inputs use one-based inclusive display coordinates. Internally the start
is normalized within the source and the end is an unrolled half-open boundary.
A circular-origin domain can end beyond the reference length; it must remain at
most one source turn and within the existing 200-base complete-oligo profile.
The preview shows both actual source spans, orientation and the 3′ base.

Endpoint changes map overlap by actual source coordinate in the original oligo's
orientation. Every retained literal call survives, including compatible ambiguous
edits, rather than replacing edits with template bases. Newly included positions
come from the actual oriented reference. Names and tails stay unchanged. Removed
edited calls are listed explicitly, so trimming cannot silently discard an intended
change. The operation never guesses an alignment after an unrelated length edit.

Tail-boundary changes preserve the complete canonical ordered 5′→3′ oligo. Moving
annealing-prefix calls into the tail reduces selected 5′ template context; moving
tail-suffix calls into annealing captures fresh 5′ context. Promoted calls remain
literal and may mismatch that context. Forward and reverse extensions use their
own biological 5′ end. Empty domains, extensions beyond a linear source and more
than one circular turn are rejected. Calls leaving annealing may remain in the
tail; this is distinguished from sequence removed by endpoint trimming.

Before/after review includes full ordered oligos, separate tails/domains, new
reference comparisons, mismatch/ambiguity positions, removed edits and warnings.
Accept resolves the exact source event again, checks the entire original pair,
workspace epoch and unchanged control values, then recomputes the proposal before
replacing only that slot. The existing pair Apply action remains separate, followed
by explicit scientific analysis. Invalid/incomplete editable text is preserved on
failure; preview does not repair or normalize it into an unintended primer.

Retained historical source versions remain usable without moving heads or silently
rebasing to the current source. Their exact document/event/revision is displayed.
The comparison is against that version, not evidence of binding to another source.
A missing event, changed original anchor, changed pair or stale preview prevents
acceptance.

Circular support is shared by selection capture, endpoint planning, saved anchors,
the toolbar and captured PCR overlays. Exactly two ordered origin pieces can seed
one contiguous primer domain. Gaps, repeated/multiturn selections and fuzzy or
unresolved locations do not become a fabricated contiguous template. Wrapped
primer eligibility does not enable contiguous sequence or mutation editing.

Native pair attachments now write `primer-pair-draft-v2` with JSON envelope version2.
The reader explicitly supports versions1/2 and both attachment kinds, validating
original event bases. Version1 remains nonwrapping. Autosave replaces either known
kind atomically with the version2 representation; default reset removes both.
Multiple active known pair drafts are rejected rather than silently choosing one.
Unfinished names/domain/tail text retains case/spacing and invalid input. The outer
native archive and binary worker schemas are unchanged.

Implementation: `web::primer_domain_edit`, its shared review component, existing
pair model/storage, selection toolbar and captured overlays. Fifteen new native
cases cover coordinate/order/edit preservation, version compatibility, source
identity, stale/tampered proposals and wrapping eligibility. Browser interaction
acceptance remains deferred until full in-scope implementation is complete.
