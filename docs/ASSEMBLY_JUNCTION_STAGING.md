# Staging edited assembly intent

The web host carries `assembly_review::IntendedReference` as either
`Regions(assembly_intent::Report)` or
`Junctions(assembly_junction_edit::Intent)`. Both expose `target()` and
`fragments()`; edited contribution output spans include outgoing synthetic
junction bases, while their source coordinates remain original-source intervals.
Synthetic bases are explicit design intent, not observed source DNA.

`stage_any(frame, intent, source_events, candidate)` dispatches the old region
path unchanged or rebuilds the edited intent from the candidate's original
`PlanInput`s, `base_request` and explicit junction edits. The edited path performs
borrowed source guards before model copies, then resolves every pinned event and
compares actual content/revisions with the preparation inputs. A fresh core intent
must match the report context, actual target and candidate's independent target
revision. Simulator products never become the expected reference.

The target and candidate are staged together in a private workspace. A target is
reused only when its complete canonical reference identity and exact source-parent
event set agree; identical DNA in a different source document/event is not enough.
Otherwise staging derives a new target with all source parents and the core's
restrictive intended-reference markers. Original source heads remain unchanged.
No scientific status, physical material state or experimental outcome is restored
by this operation; later preparation reconstruction and target verification remain
separate.

Inline native fixtures cover literal trim/insert intent (`ACG + GG + CGT`), exact
source ancestry, immutable original workspace, retry reuse, equal-DNA source-event
aliases, forged target/context, stale bindings and unchanged region dispatch.
All seven staging cases passed, including the four existing region-only cases;
new core/binary/MCP fixtures also passed. Browser checks remain deferred.
