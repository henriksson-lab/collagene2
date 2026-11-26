# Snapshot subtree clone and mount-aware folder synchronization

`repository::folders` provides an explicit, persistent remote-tree clone job.
It extends the real offline SQLite repository; it does not replace the injected
authenticated transport, upload private repository checkpoints, or run external
`git`. The user selects either a lab's virtual root or one folder with its known
parent, plus a local destination. Beginning a job persists that scope before any
network access.

The module exposes `begin_tree_clone`, `resume_tree_clone`, `tree_clone_job(s)`,
`restart_tree_scan`, `retry_tree_clone`, `accept_tree_parent_change` and
`queue_tree_folder_metadata`. Hosts can run bounded steps off the UI thread and
display the returned progress. Creating a job is not successful completion.

## One bounded catalog snapshot

Folder selection is resolved by paging its actual parent; root selection browses
the virtual lab root. The first response pins a catalog generation. Every later
discovery page and every page of every child folder requires **that same**
generation, including each child's first page. The breadth-first scan completes
before any local node is created.

Duplicate IDs, repeated ancestors/cycles, parent mismatches and excess depth/node/
page counts are rejected. Validated wire DTOs bind each reply to the requested
operation, destination folder, page and authenticated lab. Explicit hard bounds
are 2,000 nodes, depth 128, 10,000 pages and 200 entries/page; defaults are 1,000
nodes, depth 64, 2,000 pages and 100 entries/page. Existing repository state,
payload, storage and checkpoint limits still apply.

A catalog change invalidates the scan. Missing/deleted scan parents also prevent
claiming a complete snapshot. No local data was created by such a scan. The user
may explicitly restart it before imports begin; the helper never silently mixes
generations or carries partial old pages into a new snapshot.

## Exact immutable document checkout

After the scan, folders are created parent-first and documents are imported into
their mapped local parents. Each document is opened with
`revision: Some(scanned_document.head.id)`. The returned immutable revision
descriptor must equal the complete scanned descriptor, and its payload hash,
kind, selected native pins and sharing descriptor are checked by existing wire
validation. The remote's current head/name/parent may have changed; those current
fields are **not substituted** for the scanned snapshot.

A completed scan can therefore continue fetching its immutable revisions after
the remote catalog advances. Bindings retain the scanned remote metadata and
generation, so a later offline metadata push detects intervening changes.

Native documents are loaded normally. Native workspace payloads must contain the
independently derived exact selected-document dependency closure, including its
original source/event/run/capsule ancestry. Extra native documents or private UI
records are rejected, rather than silently imported under one plasmid label.
Pinned archive bytes are retained unchanged. A real saved PCR product remains
available to the ordinary fresh `run_material::reconstruct` algorithm afterward;
clone metadata itself grants no material authority.

## Persistent progress and local changes

Jobs and their bounded metadata index live in private repository checkpoints.
For each created local node, its native storage mutation, remote binding and job
mapping/progress commit in the **same SQLite transaction**, using both repository
state and job checkpoint CAS. A crash cannot commit a node while forgetting its
mapping. Restart reads the existing job and skips completed nodes, so an
interrupted fetch/response does not create duplicates.

Completed local documents are never reopened or overwritten by resume. Private
edits made after their import remain their local heads. A remote node already
bound by another checkout is an explicit blocked condition, not permission to
replace or move that existing local item. Name conflicts also block visibly.

The destination and created local folder generations are recorded. A locally
renamed/moved import parent blocks further descendant creation until the user
reviews its current generation with `accept_tree_parent_change`, then explicitly
retries. That acknowledgment keeps the private metadata; it does not restore old
names or alter remote records. Concurrent repository-state changes roll back the
new node transaction instead of overwriting another handle's checkpoint.

Progress reports distinguish `Discovering`, `Scanning`, `Importing`, `Complete`,
`Invalidated` and `Blocked`, with counts and a bounded issue. Transport failure
leaves the current phase incomplete. Explicit retry of a blocked completed scan
keeps all applied nodes and immutable remote pins. Restart cannot discard a
partially imported snapshot. Successfully importing some documents never causes
a partially completed job to be labeled Complete.

## Local mounts are not remote moves

A selected subtree can be mounted beneath an arbitrary local destination. For
its unchanged mounted root, `queue_tree_folder_metadata` preserves the original
external remote parent; it does not infer a move from where the clone was placed.
A virtual lab root maps to the explicitly selected local destination, not a
fabricated catalog folder ID. No new unbound local parent is guessed to mean the
remote root.

Actual offline rename/move intent uses the acknowledged remote generation through
the existing durable mutation outbox and receipt recovery. A remote member's
intervening rename/move yields the existing typed, durable conflict; the private
local folder remains unchanged. Renaming a mounted root can proceed without
having cloned its external remote parent. Moving it elsewhere requires an
explicit known parent mapping. Callers should use this mount-aware method for
folders from a subtree job rather than interpret the generic `queue_folder`
placement comparison as mount intent.

This increment supplies point-in-time recursive checkout and protected folder
metadata publication. It does not automatically reconcile deleted remote items,
merge divergent document heads, recursively publish an entire local tree or
silently refresh a completed clone to today's remote snapshot. Those actions need
their own explicit comparison/selection workflows. Local unpublished branches
and existing per-remote bindings remain independent of the job's publication view.

## Focused source fixtures

The fixtures use a real SQLite lab and actual native host/sync dispatch as the
injected binary transport. They cover an interrupted two-level checkout followed
by repository restart, a newer remote head after the scan, preservation of an
already imported private local edit, exact old archive bytes and fresh PCR-origin
replay, catalog invalidation before imports, root mounting and offline rename CAS
conflicts, cyclic/duplicate-ID replies, and a changed local destination requiring
explicit generation review. No live HTTP/browser checks or builds were run by
this worker; the parent coordinates the targeted compilation/test batch.
