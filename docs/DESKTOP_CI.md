# Desktop build CI

`.github/workflows/desktop.yml` automatically builds downloadable desktop
artifacts for relevant source changes: a portable Windows `.exe`, macOS `.dmg`
and `.app.tar.gz`, and Linux `.deb`. The [packaging contract](DESKTOP_PACKAGING.md)
describes file contents, prerequisites and manual development/release selection.
Scientific CI remains independent, and a platform failure does not cancel the
other package builds.

One Linux frontend job builds the current Leptos application and scientific
worker. The reusable packaging workflow receives that exact same-run artifact,
then performs one native build per platform. It skips its own frontend job when
called this way; manual package dispatch builds a fresh frontend itself. Relevant
Rust/frontend/Cargo/Trunk, packaging-script and workflow changes trigger automatic
CI. There is no extra native-only compilation after creating the packages.

Frontend construction starts with an empty `dist`, permits exactly one current
content-addressed worker bundle, and bounds the snapshot to 512 files, 128 MiB
per file, 256 MiB total and depth 16. Symlinks and special files are rejected.
Desktop compilation independently checks those bounds and the `collagene_rpc`
bridge in the actual referenced main WASM. A browser-only frontend cannot silently
replace durable native storage. The default `custom-protocol` feature embeds the
frontend, including for the portable Windows executable.

Trunk 0.21.14, wasm-bindgen-cli 0.2.128 and cargo-tauri 2.11.3 are pinned upstream
prebuilt tools, with checksums and no source-build fallback. Native dependencies
remain locked. Shared frontend/package caches are reused; Cargo uses two jobs,
no incremental compilation and no development/test debug information. Linux
packages build on Ubuntu 22.04 using `scripts/install_linux_desktop_deps.sh`, which
retains Ubuntu TLS/signature verification and avoids unrelated third-party feeds.
Windows explicitly targets x86_64 MSVC with static CRT settings for portability.

`scripts/stage_desktop_packages.py` inspects the actual platform executable,
extracting Debian packages without installing them and using the executable
inside the macOS app. It runs `desktop_mcp_smoke.py` with `--mcp`, checks expected
artifact counts and size limits, and records exact source/frontend/profile/hash
metadata. The Windows portable executable additionally passes a native-import
check. Development and manual release builds use this same headless path.

The stdio check launches no GUI, browser, listener, repository or credential
vault and invokes no external scientific program. It verifies MCP initialization,
cloning-tool availability, synthetic FASTA import and a virtual EcoRI digest.
Windows release-mode checks retain explicit inherited pipes; no console attachment
replaces them. [Microsoft's inherited standard-handle documentation](https://learn.microsoft.com/en-us/windows/console/getstdhandle)
describes that distinction.

Artifacts expire after seven days. These are unsigned prerelease downloads, not
a signed release or updater pipeline. The user's review of the finished product
is the primary usability acceptance step. No local desktop build, installer run,
GUI or browser check was added while changing this CI configuration.
