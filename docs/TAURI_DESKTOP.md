# Collagene desktop shell

`crates/desktop` is the native Tauri 2 shell for Windows, macOS and Linux. It embeds
the existing Leptos/Rust/WASM interface and calls the native desktop host through
one binary command. Local analysis and repository storage work without a network
connection. Mobile is outside this shell's scope.

The shell introduces no TypeScript application, external scientific executable,
local HTTP listener, database service, shell plugin, filesystem plugin or SQL
plugin. The native repository retains plasmids, folders, immutable revisions and
working-state recovery. Explicit synchronization with configured remotes remains
a native repository concern; remote credentials never belong in plasmid archives.

## Assets and startup

Prepare a fresh frontend build in the repository's `dist` directory before
building the shell. `tauri.conf.json` defaults `frontendDist` to `../../dist`, and
the default `custom-protocol` Cargo feature embeds those bytes for offline use.
There is no development-server URL or automatic external build command. The build
script rejects symlinks, special files, excessive asset size/depth, a missing
index or multiple historical worker versions. Do not point it at the public
server's archive of old workers. Keep the current hashed worker paths intact.

Packaging also follows the main `collagene-web-…_bg.wasm` URL actually referenced
by Trunk's module bootstrap and matching WASM preload. That exact main WASM must
contain the retained `collagene_rpc` native invoke literal. A worker or stale
unreferenced sibling cannot satisfy it. Old browser-only assets are rejected:
otherwise an ephemeral desktop window could silently fall back to temporary
IndexedDB and lose edits on restart. This is an asset compatibility check, not
an authentication or scientific-authority mechanism; native IPC still validates
its protocol and window capability normally.

The WebView runs in ephemeral/incognito mode. Native working-state storage and
the local repository provide plasmid persistence; WebView IndexedDB is not a
second durable plasmid database outside the configured root. The separate,
deferred genome/dataset IndexedDB cache is temporary on desktop until an explicit
native dataset-storage adapter is added.

The desktop GUI command accepts:

```text
collagene-desktop [--storage-root DIRECTORY] [--allow-loopback-http]
```

The same executable also accepts `--mcp` for headless stdio operation, with
optional exact saved-document selections. It does not open a GUI, initialize
remote credentials or modify the editor's working/recovery slots. See
[embedded MCP](DESKTOP_MCP.md) for arguments and selection scope. `--help` and
`--version` also exit before GUI initialization.

Remote databases normally require HTTPS. The explicit per-launch
`--allow-loopback-http` flag additionally permits literal loopback addresses for
local development or SSH forwarding; DNS names and ordinary LAN HTTP still fail
transport validation. This flag does not relax expected-lab or credential scoping.

Without an override it uses the persistent configured root, or initializes a
repository under the application's private app-data directory. The native host
uses fixed database filenames inside that root. The selected canonical root is
saved in the application's configuration directory as versioned `desktop.json`;
future launches use it again. A new explicit CLI override changes the next
startup selection after the host successfully opens it. No frontend message can
read arbitrary paths or change this restriction. A writable root overlapping the
frontend build assets is rejected so a later package cannot accidentally embed
local databases. There is no database or credential file under bundled assets.

Settings writes use exclusive temporary files and atomic replacement; Unix files
are owner-only (`0600`). Native host directory checks constrain the repository
files. On Windows, protect the application data/configuration locations with the
user's normal private ACLs. Configuration contains only a storage path and schema
version, not passwords or remote access keys.

File exports use the separate native [export service](DESKTOP_EXPORTS.md) and
the fixed `exports` child of this storage area. The UI waits for a hash/length
receipt and shows the saved filename. Export I/O holds neither the working-state
mutex nor the remote-sync lane. All export formats use this same boundary;
native mode has no browser-download fallback.

## IPC and window permissions

The embedded main window is created by native code. Remote navigation and new
windows are denied. The command checks the main label and exact platform local
app origin/document before accepting a request. Windows uses the local HTTPS
custom-scheme origin; macOS/Linux use `tauri://localhost`. There is no localhost
development-server exception. CSP forbids frames, embedded objects, arbitrary
connections, base-URL changes and form navigation; local WASM, styles, images and
workers remain available. The filesystem asset protocol is disabled.

`build.rs` declares `collagene_rpc` through
`tauri_build::AppManifest::commands`. The only capability targets the local main
window and grants only `allow-collagene-rpc`. It has no remote URL grant, default
core/window grants or plugins. This prevents newly added commands from silently
becoming available to the frontend. Native path configuration and authorization
handles are not deserialized from frontend claims.

The frontend uses the public `window.__TAURI__.core.invoke` API enabled by
`withGlobalTauri`. It supplies a `Uint8Array`/`ArrayBuffer` request, and the native
command returns `tauri::ipc::Response` containing raw response bytes. Native
bounded/versioned protocol decoding, exact operation IDs and repository CAS
remain authoritative; the shell does not reconstruct scientific state from
presentation JSON. JSON IPC bodies are rejected.

Blocking native work runs off the UI thread. One admitted native job bounds
pending repository work and request copies (maximum 129 MiB for a paired working
archive replacement, maximum 65 MiB per response), while small native
control/cancellation messages can be handled separately without waiting for the repository mutex. Dropping a frontend
promise cannot cancel a transaction that has already committed: retain the same
operation ID and query/retry through the durable protocol. Native host fixtures
exercise actual headless repository and cancellation behavior without launching
a WebView.

## Build and validation status

The new shell is registered by the root agent only after its host/protocol graph
is ready. Build checks may compile cached Tauri/GTK/WebKit dependencies; they do
not launch the GUI or install SDKs. Packaging/signing for all three operating
systems is distinct from a successful native type check and remains explicit.
No browser, GUI, listener or application HTTP check is required by this shell's
source integration. The startup ACL/CSP and configuration are reviewable in the
new crate; actual native host tests live with that host.

Tauri's official [raw request documentation](https://v2.tauri.app/develop/calling-rust/#accessing-raw-request)
describes the binary invoke boundary; its [capabilities documentation](https://v2.tauri.app/security/capabilities/)
describes local-window grants and why remote/iframe access must be constrained.

The asset compatibility parser also has three std-only regression fixtures. They
can be compiled directly without the Tauri/native graph:

```text
rustc --edition 2024 --test crates/desktop/tests/assets_guard.rs -o /tmp/collagene-desktop-assets-guard
/tmp/collagene-desktop-assets-guard
```

They exercise the actually referenced native main (including a streaming-buffer
boundary), rejection of a marker found only in a worker/stale sibling, and exact
Trunk bootstrap/preload agreement. They do not instantiate WASM or launch a GUI.


## Native remote controls and independent saves

The shell routes both DTQ working/catalog operations and DSQ remote controls
through `collagene_rpc`, preserving the one approved main-window capability.
Each namespace has its own bounded admission counter and frontend request lane.
DSQ uses a separate native repository connection; HTTPS never holds the
working-state host mutex or a working SQLite transaction. Intervening local
catalog edits are checked by immutable-head and repository-checkpoint CAS.
Cancellation/status still share the branded window operation registry.

The desktop interface exposes a `Databases` accordion beside the local library.
Remote keys may remain in session memory or be explicitly remembered in the
OS vault; no key belongs in a plasmid, exported workspace or repository setting.
Clone and pull update the local catalog. Users reopen those records to review
them; sync never silently replaces the active working editor.

Tauri 2.11.5 is pinned with its matching runtime 2.11.3/runtime-wry 2.11.4,
macros/codegen 2.6.3 and utils 2.9.3 family. Newer runtime/macro releases selected
by their broad transitive constraints currently use APIs absent from that exact
facade version. Upgrade this cohort together and retain the lockfile.
