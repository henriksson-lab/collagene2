# Adopt a checked recovered desktop library

The desktop saved-library backup controls now distinguish three actions:

1. Verify a known backup and restore it to a new managed copy.
2. Select that copy from its receipt or the paged managed-recovery list, then
   review its native hashes, lab scope, target contents and the current catalog generation.
3. Explicitly choose **Use recovered library after restart**. A native receipt
   confirms durable selection for the next desktop start, not live activation.

Until Collagene is closed and reopened, the current window states **pending
restart**. Saved-library writes and synchronization are paused. Current editor
contents, working-editor autosave/recovery and ordinary file export remain
available. On the next successful native startup, the host independently checks
and opens the selected managed database; only its native `ActiveAfterRestart`
status permits the UI to say the recovered library is active.

The private storage root, working.sqlite and its current/previous recovery frames
remain in place. The original repository.sqlite is preserved as a rollback source;
no database bytes are overwritten or merged. Credential storage stays scoped to
the existing private root and trusted native identities. Restart clears in-memory
sessions; restoring a saved catalog does not copy secrets into a backup/archive.
An append-only native activation journal supports repeated checked selections.
The controls offer explicit rollback to the preserved original saved library and
selection of prior managed recovery copies; no old catalog is deleted. Each
review pins the current activation generation, source state, and a fresh
WAL-aware SQLite snapshot of the target. Previously active copies are reviewed
using their current contents, rather than pretending their original backup
receipt describes later edits.

The shell uses a shared catalog admission gate for local library requests, native
synchronization and catalog management. Adoption requires an exclusive lease;
it refuses while a catalog job is running, including remote I/O. After native
acceptance, even a lost response leaves the shell's catalog gate paused. Small
activation-status queries remain available. Working-editor commands, controls,
source import and confined file exports do not take this catalog gate.

The UI initially pauses saved-library controls until it has checked activation
status. A pending result invalidates the local library client epoch and clears
its selected destination/context; old callbacks cannot select the recovered
catalog. Native window capabilities, exact review generation/checkpoint and
host-side connection checks remain the authority. A UI boolean cannot adopt a
library or grant access to another path.

Source fixtures cover catalog-gate exclusivity/pending behavior and exact
copy/scope/generation/receipt binding. Native backend fixtures independently cover
actual restore, durable adoption and startup. Parent native/WASM validation is
recorded separately from the source contract. Packaging and GUI acceptance remain
unfinished; these controls do not restart the application automatically.

Rollback follows the same review, native acceptance and restart flow. The UI
shows both current and selected locations and never treats rollback as a merge
or rewrites the working editor to match the older catalog. Old local/remote
revision bindings may require explicit reconciliation after a rollback; ordinary
immutable-head and checkpoint CAS continue to reject silent overwrites.

Every writable desktop catalog connection retains the same native process lease
inside its SQLite owner. An escaped synchronization repository or backup service
therefore continues to exclude another desktop startup after the main host is
dropped. Opening ordinary secondary services shares this lease rather than
acquiring another OS lock; a standalone native backup service acquires its own.
The one-shot native keepalive cannot replace a previously retained owner.
Read-only MCP connections continue to coexist and neither initialize the catalog
nor apply pending selections. Generic repository/storage constructors remain
trusted native APIs; they do not discover or arbitrate desktop roots themselves.

The activation journal accepts at most 4096 generations. Once full, a new
adoption is refused before creating review snapshots or journal records, while
exact retries of previously accepted operations still return their original
receipt. Source fixtures cover the last legal generation and overflow refusal,
and escaped writable-handle lifetime with concurrent read-only MCP access.
These new fixtures have been added but have not yet been run.

A failed startup adoption no longer blocks the working editor. A changed source,
corrupt target or unreadable selector starts a native editor-only session:
working-file metadata, current/previous autosave, preferences, source import and
confined exports remain usable. Saved-library and synchronization services stay
paused; the interface never invents a lab identity or silently switches catalogs.
The new version-4 recovery bootstrap returns real working state with an optional
catalog context. The original bootstrap format remains available for healthy
catalog sessions.

When the journal still identifies an exact pending selection, the recovery panel
shows it and offers **Cancel pending switch; keep previous library**. Cancellation
is pinned to the accepted operation, review digest and activation generation.
Native acceptance appends a terminal cancellation record and consumes that
generation without deleting its acceptance or changing the selected catalog.
Both database branches and working-editor files are preserved. Cancellation does
not rebind live connections: restart before using the previous saved library.
A lost response leaves controls paused and triggers a fresh native status query.
A corrupt journal with no authenticated pending identity offers no cancellation;
its diagnostic remains visible while editor autosave and export remain available.

Pure source fixtures pin cancellation to the exact accepted review and generation,
and ensure editor-only admission permits status/cancellation while rejecting
catalog work. These additions await the parent-managed compiler/test batch; no
GUI or browser check was performed.

After an adoption or cancellation attempt, shell admission rechecks the durable
process-generation status. A rejected request or exact historical receipt replay
without a journal transition cannot pause a healthy current session. A new pending
selection, consumed cancellation generation or unreadable journal remains paused.
