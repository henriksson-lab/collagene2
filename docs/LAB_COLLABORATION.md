# Revision-bound lab discussion and human review

The native storage layer in `collagene-lab-library::collaboration` now stores
comments, thread resolution and human review against an exact **library document
ID and library revision ID**. It is not yet exposed by the server or application
UI. Individual authentication, binary service messages and review controls remain
required before this is an available shared-library feature.

A library revision pins the complete immutable stored payload and its selected
native document/event. Renaming or moving its catalog entry preserves discussion.
Saving another revision starts with no discussion or human approvals; the previous
revision retains its own history. Discussion does not change molecular data,
catalog generations, native verification reports or scientific authority.

## Identity and permissions

The operator's `LabAccess` capability manages storage. Authoring discussion instead
requires a non-serializable `MemberAccess`, constructed by trusted host code with
an authenticated stable actor ID, current display name and current role. Actor
IDs are scoped by the lab. Historical events retain the posting identity/name and
role. Display names and text embedded in a comment cannot grant permission.

- Reader: read discussion only.
- Contributor: start threads, reply to open threads, submit or withdraw their own
  revision-specific review; resolve or reopen threads they started.
- Maintainer: the same actions, plus resolve or reopen any thread in the lab.

The host must check current membership for every request and issue a fresh
capability. The existing shared lab key does **not** establish individual authors;
it must not be presented as doing so. Role restrictions here govern discussion,
not the separate operator-only native document APIs. The service must enforce its
member policy on all exposed writes when individual authentication is added.

## Events, concurrency and bounds

Events are append-only: thread opening, reply, resolution/reopening and human
review (`Approved`, `ChangesRequested`, `Withdrawn`). There is no destructive edit
or deletion API. A thread's opening sequence number identifies it within its exact
revision; replies cannot target another reply or a thread from another revision.
A resolved thread must be reopened before replying.

Each append requires the exact discussion generation. SQLite performs permission,
target, generation and quota checks in one immediate transaction. Two simultaneous
writers cannot both append at the same generation. After a lost acknowledgment,
resubmitting that generation conflicts instead of duplicating the event. The
client must refresh and review before choosing to resubmit. Paged reads bind the
discussion generation; an intervening append invalidates further pages.

A review decision applies only to its author and exact revision. Later decisions
from that author supersede their earlier decision when presenting current review;
all events remain available. A human approval never proves that a cloning design
is valid and must not substitute for deterministic reconstruction or verification.

Bodies are plain text, limited to 8,192 UTF-8 bytes. Render them as escaped text,
never executable HTML/Markdown. Discussion is bounded to 4,096 events per revision,
100,000 events and 16 MiB of charged event data per lab. Quotas include action,
author ID, display name and role bytes, independently of sequence archive quotas.
Pages contain at most 200 events. Failed writes leave events and counters unchanged.

## Storage and recovery

The current library SQLite schema is **2**. Schema 1 is explicitly rejected;
there is no automatic migration, fallback reader or compatibility mode. Do not
replace a running service's binary against an older database and assume it will
open. Service cutover requires an explicitly chosen current-schema storage area;
existing live databases have not been opened or altered by this change.

Whole-database operator backups retain discussion along with documents, revisions
and private repository records. Verification replays the event grammar and checks
thread ownership, resolution transitions, per-revision sequence continuity and
lab quota accounting. Restoring a backup does not authenticate its historical
actors or turn their reviews into molecular evidence. Public document export and
chosen-remote sync do not yet transfer discussion; no identity or discussion
copying is implied by ordinary plasmid sharing.

The native tests cover revision isolation, role enforcement, exact authorship,
concurrent independent SQLite writers, lost acknowledgments, snapshot paging,
text bounds, immutable events, quota refusal, recovery and unsupported-schema
rejection. They do not establish completion of the service or graphical workflow.
