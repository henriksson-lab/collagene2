# Embedded desktop MCP

`collagene-desktop-mcp` prepares a memory-only MCP session from explicitly
selected **saved local catalog revisions**. The desktop executable can serve that
session over stdio in the same process; it neither launches the GUI nor runs an
external scientific executable. The existing 179-tool MCP surface is unchanged.

The shell/operator chooses selections as `doc-32hex@revision64hex`. The revision
is an immutable local library revision ID, not a molecule content revision or a
request to follow the latest head. At most 32 unique pairs are accepted. Missing,
malformed or corrupt pins fail session creation; no partially prepared session
is returned or served.

## Native API

- `NativeSelection::parse(text)` / `FromStr` parse exact catalog pins.
- `Session::empty()` makes a bare session without filesystem access.
- `create(repository_root, selections)` opens only the fixed private native
  repository when at least one selection exists.
- `create_with_storage(repository_root, selections, Some(managed_storage_root))`
  uses the existing `CollageneServer::with_storage_root` policy. Managed genome
  scratch remains under the operator-configured root, rather than falling back
  to the OS temporary directory. This does not broaden managed-genome tools.
- `Session::server()` supports an embedded in-process transport; `into_server()`
  transfers ownership. `serve(session).await` serves stdio directly.

The root arguments are trusted shell configuration. MCP tools cannot provide
paths, select further native database records, read credentials or write the
local catalog. `open_local_repository` checks the native private-root policy and
opens `repository.sqlite` only. No MCP path opens, reads or writes
`working.sqlite`, current/previous recovery slots or view preferences. Bare
`--mcp` requires neither login nor a repository.

## Exact source scope and material evidence

For each workspace payload, the factory loads and validates the exact saved
archive and its selected document/revision/event pin. It derives the selected
history/dependency closure with `collagene-library-sharing`, including original
run sources, independent targets, product events and supported origin-capsule
dependencies. It preloads **that closure**, never the whole private archive.
Unrelated private documents, unpublished branches and presentation sidecars do
not enter the MCP session. Repository bindings, outbox, remote configuration and
credentials are outside biological payloads and are never preloaded.

A bare `.cg2` record loads through the existing native document decoder as
literal molecule data. It receives no fabricated event, run history or physical
material authority. Imported history is documentary evidence: ordinary MCP
replay and `prepare_run_material` independently reconstruct retained exact
sources with the existing deterministic algorithms. Neither a closure inventory
nor a historical successful report conveys current private tracked material.

Preload input is capped at 64 MiB aggregate in addition to the existing MCP
session memory/item limits and native archive bounds. A failure discards the
unserved session. Native helper `preload_archive_bytes` avoids a base64 roundtrip
but uses the existing validated archive import path; `preload_document_bytes`
uses the existing `.cg2` loader and immutable document insertion. These helpers
are native embedding APIs, not additional tools.

All AI-created documents, plans, artifacts and finalized outputs remain in the
MCP session. Exporting bytes to the MCP caller does not save or sync them into
the local catalog; later user-directed import is a separate operation. stdout
is exclusively MCP protocol traffic. The shell owns CLI routing, stderr
messages and choosing whether to start MCP instead of the GUI.

## Desktop command

Configure the AI client's stdio server with the installed desktop executable
and `--mcp`. This empty session can import user-provided sequence data and run
the embedded algorithms fully offline.

To expose saved local plasmids, use the local library document menu's **Copy MCP
reference** action, then configure arguments as:

```text
--mcp --storage-root /operator/chosen/private-root --mcp-document ID@REVISION
```

Replace `ID@REVISION` with the copied reference and pass each argument separately
in the client's configuration; the directory must match the desktop repository.
Repeat `--mcp-document` for additional explicitly chosen revisions. The copied
reference contains no filesystem path, remote credentials or live-head alias.
Saving a newer revision does not change a running session's selected revision.

## Focused fixtures

The in-process SDK fixtures use a real native SQLite repository and duplex MCP
transport. They cover strict pin parsing, empty-session zero filesystem access,
a saved PCR product with exact FromRun reconstruction and independent target,
exclusion of a private unrelated plasmid, unchanged local catalog/checkpoint,
absence of the working database, literal `.cg2` without invented workspace
history, and atomic failure when a later selected revision is missing. They use
no browser, GUI, HTTP server, network service or operating-system vault.
