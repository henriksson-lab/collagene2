# Revision-bound lab discussion and human review

The native storage layer in `collagene-lab-library::collaboration` now stores
comments, thread resolution and human review against an exact **library document
ID and library revision ID**. The server exposes bounded binary discussion messages
using operator-configured individual credentials. The web library's revision
history provides discussion/review controls; graphical acceptance remains open.

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
capability. The server now requires individual member credentials, resolved through a trusted
operator registry; it accepts no client author/role authority. Every append pins
the actor ID seen during review; a changed authenticated actor refuses the write
rather than silently authoring it as another person. Reader credentials also
refuse ordinary document writes and reliable sync mutations. Role restrictions
in the storage module govern discussion; separate native operator APIs retain
their trusted administrative capability. [Server configuration](SERVER_LIBRARY.md).

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
rejection. Additional native binary SDK and credential fixtures check the service contract,
author binding, role restrictions and registry revocation. They do not establish
HTTP or graphical acceptance.


## Web workflow

Connect using your individual member key, choose a plasmid's **Revision history**,
then expand **Comments and human review** on the exact saved revision. **Load /
refresh discussion** fetches its current generation and trusted member identity.
Use **Load more events** until all current events are available before posting or
interpreting current review decisions. Partial/stale pages do not present a
current approval summary.

Post a new comment, choose **Reply** on an open thread, or resolve/reopen an
eligible thread. Review controls append your decision for that exact revision:
**Approve revision**, **Request changes**, or **Withdraw my review**. The optional
note uses the same plain-text draft. Resolved threads retain their messages;
withdrawal retains earlier review history. Latest decisions are shown per author,
without turning them into a collective molecular-verification status.

A generation conflict or uncertain transport result invalidates posting controls
and retains an unacknowledged draft. Refresh and review before resubmitting; the
application never retries a write automatically. Once a write is acknowledged,
the draft clears and the first page refreshes. Failure of that refresh is reported
as a successful save with unavailable refresh, not a failed write. Changing the
member session or closing the revision view discards pending display updates.

In Tauri, choose a saved remote database, open **Remote revision history**, then
expand the same **Comments and human review** controls on an exact remote revision.
Requests use the native credential vault and the complete configured remote/lab
binding; browser code receives neither the remote bearer nor a filesystem path.
Changing remotes or replacing/forgetting a credential closes prior revision views.
The identity precondition also refuses a different person after another window or
operator changes a credential.

Offline documents and chosen-remote plasmid synchronization remain independent.
Loading/posting discussion requires an explicit online action. No discussion or
actor identity is silently copied to another database or attached to local
scientific history. There is no offline comment queue or automatic write retry.
Lost acknowledgments require refreshing the exact remote discussion before choosing
to resubmit. These are service discussion records, not reagent/ELN/stock records.

The discussion wire schema is **2**, and desktop sync IPC schema is **4**. Older
internal versions are rejected directly without fallback decoding. Ordinary
library document and reliable-mutation message formats remain their independent
current schemas. Native discussion paging permits 2 MiB responses; desktop IPC
permits 3 MiB frames to carry a full bounded discussion page plus trusted remote
context. The sequence/workspace archive limits are unchanged.
