# Private native catalog saves and public publication

A Tauri Local library save stores the complete bounded native workspace archive,
with the exact selected document/revision/head-event pin. Original SBOL exchange
bytes, other local documents, native undo/redo branches, candidates and cloning
history stay in the private local repository. Saving locally needs neither a
server login nor a public omission review. This matches the repository's
`private_workspace` contract rather than selected-document publication.

The App chooses `library_documents::capture_private_local` only from the trusted
native IPC capability (`desktop_rpc::is_desktop`), never from a URL or hostname.
It rejects public review choices/approval and mismatched optional source pins,
validates the complete workspace, and uses the existing bounded archive worker.
The App rechecks its source epoch after encoding. The payload uses existing
`ExplicitWholeWorkspace` scope; the embedded DesktopSession sends it only through
the approved local DTQ command. The native host validates exact archive/native
selection and commits through durable intent-bound local mutation/CAS receipts.
No native identity or scientific status is fabricated by capture.

The browser's shared-library Session rejects non-DocumentClosure workspace
mutations before encoding, request admission or HTTP fetch. It cannot upload a
private catalog snapshot accidentally. Public capture retains strict native
selected-document closure and the separate explicitly reviewed SBOL subset path.
Generic server whole-workspace scope policy is unchanged; this browser pane is
not an interface for granting it. Remote pushes independently derive public scope
from saved local data; native reviewed-SBOL push remains separate future work.

Pure native tests in web `library_capture` preserve exact original attachments
and full native graph through local capture, check selected pins and refuse both
private import/save requests at the public transport gate. The desktop-host
`private_catalog` fixture exercises actual offline IPC import, host restart,
receipt replay without duplicate revisions, open and immutable history bytes.
These are in-process fixtures, with no browser, HTTP, vault or GUI calls.

Reopening a native private catalog record uses an opaque `OpenGuard` origin bit
set only by the embedded `DesktopSession` capability. It is not a wire field or
hostname assertion. For that trusted local origin with `ExplicitWholeWorkspace`,
`project_private_document_closure` restores only the selected document and native
run/dependency closure, together with exact original SBOL attachments required by
those immutable revisions. It does not restore unrelated snapshot documents,
candidates or UI/draft state. The shared private merge entry point retains the
same immutable collision, head/CAS, ancestry and displaced-bookmark checks as
public import. The current working workspace and its unrelated local documents
remain protected. Historical local `DocumentClosure` payloads continue through
the original strict path, as do ALL browser/server-origin payloads regardless of
their claimed sharing enum. A private scope inventory never authorizes upload.

The sharing `private_open` fixture target exercises exact-original private open
and merge, exclusion of unrelated snapshot documents, unchanged native identities,
public-path refusal, stale/dirty local head protection and saved PCR/capsule reuse.
Private documentary traversal has its own one-million-reference budget in
addition to the bounded native dependency traversal; each SBOL document retains
the 16 MiB/10,000-object parser bounds. Missing-original public subsets keep their
existing exact canonical/receipt checks instead of fabricating an original.
