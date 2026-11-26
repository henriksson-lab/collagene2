# Reviewed offline folder and plasmid metadata publication

`Repository::prepare_metadata_publication` produces an immutable review for a
chosen local item and configured remote. Targets are Folder, DocumentMetadata
and DeletedFolder. The returned description records trusted local lab identity,
exact configured remote, repository checkpoint, current local generation/name/
parent, acknowledged remote parent binding and ONE proposed metadata action.
`None` remote parent means the remote lab root, explicitly visible for review.
A nested local parent must already have a binding to the chosen remote. Preparing
never creates parents, queues changes, deletes data or contacts a network.

The dedicated `MetadataAction` type admits only CreateFolder, RenameFolder,
MoveFolder, RenameDocument, MoveDocument and DeleteEmptyFolder. It has no biological
payload, SQL, credential or arbitrary mutation extension. Folder/document names
use the existing bounded Unicode display-name rules. Renaming precedes moving
when both differ; the second operation requires a new review and the first
acknowledged generation. Cloned mount-root placement is derived from actual persisted tree-job membership
before ordinary local-parent bindings. A virtual lab root mounted inside a local
folder stays the remote lab root, even if that mount folder later receives its
own remote binding. A selected subtree root retains its acknowledged external
remote parent while its local placement still equals the recorded destination.
Direct plasmids of a mounted lab root follow the same rule. Harmless rename is
therefore possible under unbound mounts; unchanged mount placement never becomes
a proposed remote move. Explicit relocation outside the recorded mount uses the
new local parent's acknowledged binding. Conflicting recorded mount meanings
refuse instead of choosing one.

An explicit UI action queues the reviewed description and its versioned digest
through `queue_metadata_publication_checked`. The digest checks representation
and context, not user authorization. Host/window approval and configured native
repository authority are separate. The native method reloads trusted state and
checks the current local lab and full remote configuration. If a matching outbox
intent already exists for the same item/remote, it returns the same immutable job
and mutation nonce, including after restart or a lost reply. An incompatible
pending content/metadata intent returns Pending. This is a replay of an already
queued operation, not fresh permission to act on current data.

Otherwise the method re-derives the complete current description and requires
its digest and captured repository checkpoint to match. It queues the EXACT
reviewed metadata operation, not a reread-dependent replacement. Existing outbox
checkpoint CAS commits make a concurrent supported repository mutation invalidate
the review/queue instead of overwriting state. All local writers must retain the
existing Repository transaction/checkpoint discipline; bypassing it with direct
SQL is not a supported concurrent application path.

Queuing performs no network request. Existing `resume_push` first looks up the
same durable mutation receipt, then sends only an unknown intent. Expected remote
metadata generations, folder-cycle checks, cross-lab scope and empty-folder
checks remain enforced by the remote transaction. A local folder must already
have been explicitly deleted (and previously bound) before a DeleteEmptyFolder
review can be prepared. A nonempty remote folder is rejected without deleting
its other members' data. No document delete or recursive deletion is introduced.

The new repository `folder_publication` fixture target uses real SQLite and
in-process native binary dispatch, covering root creation, nested-parent binding,
review staleness across two local connections, descriptor tampering, lost reply/
restart/exact nonce reuse, separate document rename+move, and remote nonempty
folder protection. It performs no HTTP, GUI, vault access or external process.
Implementation-worker formatting/diff checks pass; parent runs native compilation
and fixtures. Native DSQ/UI integration is parent-owned follow-up, not claimed by
this storage module alone. Automatic folder reconciliation remains separate work.

Deleted-folder selection requires no pasted native ID. The read-only
`deleted_folder_publications(remote)` discovery returns
`DeletedFolderPublication { local_folder, remote_folder }` records only for
that configured remote's retained bindings whose local folder lookup is actually
NotFound. Present local folders and never-bound folders are excluded; other
lookup errors fail instead of being mistaken for deletion. The acknowledged
remote folder snapshot supplies the display name. Results are a stable local-ID
ordered complete discovery of locally deleted bindings for that remote, capped
at 10,000 matches (explicit error above the cap). It is not a remote deletion
scan. The native host pages these records with a repository checkpoint pin so
entry 201 can be inspected without deleting earlier entries. Discovery creates no mutation, outbox or remote request;
selected items still require a fresh ordinary descriptor review and CAS-protected
remote empty-folder operation.

Two additional real tree-clone fixtures cover virtual lab-root folder/direct
plasmid and selected subtree/external-parent mounts. They verify unbound and
later-bound destination behavior, harmless separately reviewed renames, absence
of phantom moves, and explicitly requested relocation into another bound parent.
