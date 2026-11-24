# Managed native local storage

This MIT crate provides the local I/O policy boundary for native Collagene hosts. It has no molecular model, MCP, HTTP, async-runtime or scientific dependency. A future Tauri host can reuse the same policy instead of giving web content ambient filesystem access.

`StorageRoot::existing(path)` accepts a trusted native configuration path to an **already existing directory**. It rejects a file, missing directory and a final symlink/junction. Parent aliases resolve once to a canonical root. The root directory is opened as a capability; its final component is opened without following symlinks, and Unix identity is checked against the initial metadata. `temporary()` preserves the default OS-temp parent behavior.

`StorageRoot::session()` creates a unique nonrecursive session child; `ManagedSession::create_file()` generates an opaque filename and uses exclusive creation. There is no read-by-path, caller-selected filename, recursive deletion or public filesystem server API. Actual reads/writes use the retained file handle. Creation and cleanup use retained directory capabilities rather than reopening an absolute path. Session/root path strings are diagnostic only.

On Unix, session/file creation modes are0700/0600, subject to a more restrictive process umask. Windows uses the configured parent's inherited ACL: choose a folder already protected for the intended account. Managed Windows file handles exclude delete sharing and use delete-on-close. Normal Unix file cleanup verifies the retained file identity against its capability-relative name while the handle is still open; observed replacements are left untouched. Empty session directories are removed nonrecursively after owned files close. External additions/renames or process termination can leave data behind; there is no automatic cleanup of unrelated files or previous sessions.

This confines **managed application I/O**. It is not an OS sandbox and cannot stop other code or processes running with the same account from reading files, moving directories, changing permissions or racing namespace changes. In particular, portable empty-directory removal on Windows must close its pinned handle before removing the name. The capability boundary still prevents an MCP/browser-supplied path or symlink from granting access to another host path. No claim of complete same-user race isolation is made.

Dependencies are pinned to cap-std4.0.3 and cap-fs-ext4.0.3. Their directory-capability primitives support Linux, macOS and Windows; their permissive dependency licenses are retained through Cargo metadata. This native-only crate is not part of the scientific-core/web/WASM dependency graph.

Validation: `cargo test --locked -p collagene-local-storage`. Tests cover existing-root validation/canonicalization, session isolation/lifecycle, Unix modes, configured/root/session/file symlink replacements, preservation of outside sentinels and regular-file replacement ownership. Windows tests cover delete-sharing denial and handle-based cleanup.

The native MCP entry point exposes this policy as `collagene-mcp --storage-root /path/to/existing/CollageneData` (Windows: `collagene-mcp.exe --storage-root "C:\Users\you\AppData\Local\CollageneData"`). See [folder setup and client configuration](../../docs/LOCAL_STORAGE.md) for complete examples and the policy boundary.
