# Explicit primer pair draft

`web::primer_pair_model` is a pure Rust view model shared by manual and mutation-PCR editors. It does not run a binding search or mutate a sequence. A pair defaults to named Forward and Reverse slots; each slot retains its editable name, annealing domain, 5′ tail and optional selection anchor.

`SelectionSeed::capture` takes the actual molecule, an explicit native document/event/revision binding, and one zero-based half-open, nonempty contiguous DNA interval of at most 200 bases. It checks bounded model content before hashing and requires the actual molecular revision to match. The host remains responsible for establishing that the supplied event belongs to that document; an event-shaped string alone is not authority. Circular wrap, multipart and caret selections are not silently flattened.

`seed.oligo(name, reverse, tail)` copies the selected reference bases or their reverse complement into a 5′→3′ oligo. The anchor retains the same original source coordinates, explicit direction and original selected DNA in oligo orientation. Introduced tails are separate fields. Reverse means reverse complement, not reversed letters. Neither selection nor a saved primer proves current binding or a physical oligo stock.

`OligoDraft::primer` validates a trimmed, nonempty, control-free name of at most 256 bytes and bounds raw ASCII domain/tail input to 1600 bytes before parsing. DNA normalization accepts existing IUPAC calls and ASCII whitespace; the normalized complete ordered oligo, including its tail, must fit the existing 200-base profile. Conversion does not validate intended binding or uniqueness.

`comparison` compares only an explicit retained anchor with an equal-length edited domain. Offsets are zero-based in that oligo's 5′→3′ annealing domain, excluding the tail. Disjoint IUPAC sets are mismatches; compatible sets with any ambiguous call are uncertain comparisons, including N against N. Equal concrete calls are exact. A length change returns an unresolved-placement error; it never guesses a realignment or moves source coordinates. A manual or saved oligo without an anchor returns no comparison.

Changing sources or options must invalidate prior scientific results in the containing controller. Draft editing and comparison alone do not generate a PCR product, infer synthesis, or modify a candidate plan. Four inline native fixtures cover independently specified forward/reverse IUPAC examples, literal mismatch offsets, stale/malformed bindings and bounded input/order conversion.

## Connected controls

The app owns a session-local pair draft shared by Primers and Mutagenesis. Switching
panels or source documents retains editable fields; replacing/restoring a workspace
resets the pair. Persisting unfinished pair drafts across reload is still pending.
The selection toolbar offers forward and reverse-complement creation; it and the
editor use the same exact contiguous-selection eligibility rules. Existing source
annotations and ordered assembly primer tables remain independently saved data.

Either slot can copy a forward or reverse selected domain while preserving its name
and tail, or explicitly reuse a named saved primer through the existing strict
source-event resolver. Selection comparisons retain the original event, coordinates
and oriented bases; manual/saved copies do not invent a selection anchor. Local
comparisons inspect at most the existing 200-base oligo profile. An index of retained
event identities avoids cloning/hashing a full source molecule on every keystroke.

Applying the shared pair from Primers fills the existing ordinary PCR controls.
Applying it in Mutagenesis retains the exact draft and source/target event context;
simulation checks that binding directly and results compare the captured draft,
events and workspace epoch before installation. Editing the draft invalidates the
applied pair and prior mutation evidence. Standalone primer evidence also retains
its exact originating event for navigation and annotation; identical DNA in another
document or historical event does not authorize applying that evidence there.

Four focused native model cases and the final web/worker WASM compiler check pass.
Full GUI behavior remains subject to the
user-deferred browser acceptance review. No new scientific engine, protocol tag,
archive version, external runtime program or positive experimental claim is added.
