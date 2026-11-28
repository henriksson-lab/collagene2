# Restrict managed native storage

The MCP server can keep all its managed file-backed genome imports/datasets beneath a folder you choose:

```sh
mkdir -p "$HOME/CollageneData"
chmod 700 "$HOME/CollageneData"
collagene-mcp --storage-root "$HOME/CollageneData"
```

On Windows, create a folder protected for your account, then use its absolute path:

```powershell
New-Item -ItemType Directory -Force "$env:LOCALAPPDATA\CollageneData"
collagene-mcp.exe --storage-root "$env:LOCALAPPDATA\CollageneData"
```

For an MCP client configuration:

```json
{
  "mcpServers": {
    "collagene": {
      "command": "/absolute/path/collagene-mcp",
      "args": ["--storage-root", "/absolute/path/CollageneData"]
    }
  }
}
```

The folder must already exist. A missing path, regular file, final symlink or Windows junction is rejected before MCP starts; there is no fallback to another folder when an explicit root fails. The server opens a canonical directory capability, creates a unique session child and uses generated exclusive filenames and retained file handles. It never opens existing user files in that folder. Ordinary documents, candidate plans and artifacts remain memory-only. Without this option the existing private OS-temp session behavior is preserved.

Only trusted native startup/configuration accepts the root path. No MCP tool or browser endpoint accepts an arbitrary local read/write path. The public GUI preview serves a separate fixed static asset root and does not grant access to this native storage folder. Public sequence URL imports remain bounded network downloads; selecting a storage root does not turn a URL into local filesystem access.

`--help` and `--version` write to stderr. Argument/configuration errors exit2 with a JSON error category/message on stderr. Stdout remains MCP-only. These commands do not start the MCP session:

```sh
collagene-mcp --help
collagene-mcp --version
```

## Boundary and cleanup

The reusable [collagene-local-storage](../crates/local-storage/README.md) crate has no science, MCP, HTTP or async dependency. It provides `StorageRoot`, private session creation and opaque owned files; a future Tauri native host can reuse it. Later I/O and cleanup are capability-relative, not reopened absolute paths. File imports/exports still use bounded chunks and opaque session handles, with the existing reference reservation/capacity rules.

Unix session/file creation modes are0700/0600, subject to a stricter umask. Windows inherits the protected parent folder's ACL. Managed Windows files request DELETE access, deny FILE_SHARE_DELETE and use DELETE_ON_CLOSE, so cleanup targets the retained owned handle. Unix cleanup checks the retained file's device/inode against its session-relative entry before unlinking while open; observed replacement files are left untouched. Session cleanup only removes an owned empty child directory, never recursively removes user content or the configured parent.

This is an application-managed I/O boundary, **not an OS sandbox** against other processes or code running as the same user. Such code can change permissions, rename directories, retain file handles or race namespace operations. Portable Windows empty-directory removal requires closing its pinned handle before removing its name; no complete same-user race-isolation claim is made. Capability-relative operations still cannot be redirected by a client-supplied path or symlink to arbitrary host files.

Normal disposal cleans owned files and empty session directories. Crashes, forced process termination, external modifications or foreign handles can leave files behind; old sessions and unrelated files are never automatically swept. The configured folder should be private and have adequate disk space.

## Verification

Focused native checks:

```sh
cargo test --locked -p collagene-local-storage
cargo test --locked -p collagene-mcp --test workspace storage_root::
cargo test --locked -p collagene-mcp
```

The development stdio smoke accepts the same optional root and exercises the real executable:

```sh
python3 scripts/mcp_smoke.py /absolute/path/collagene-mcp --storage-root /absolute/path/CollageneData
```

Tests cover root validity, canonicalization, session isolation, platform permissions, exclusive files, symlink/replacement sentinels, cleanup ownership, CLI stdout isolation and a complete configured-root SDK reference import/export/disposal lifecycle. Windows-specific delete-sharing/cleanup tests remain enabled on Windows. Scientific core, web and worker dependency graphs do not acquire this native storage crate.
