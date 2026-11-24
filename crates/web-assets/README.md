# collagene-web-assets

MIT, build-time pure Rust packaging of a complete immutable Trunk analysis-worker bundle. It launches no subprocesses and makes no network requests. The application runtime does not depend on this crate.

```rust,no_run
use collagene_web_assets::package_worker;
let manifest = package_worker(std::path::Path::new("dist/.stage"), "/collagene2/")?;
println!("{}", manifest.loader_url);
# Ok::<(), collagene_web_assets::Error>(())
```

`package_worker(staging: &Path, public_url: &str) -> Result<Manifest, Error>` consumes:

- `analysis-worker_loader.js`, with only `importScripts('./analysis-worker.js'); wasm_bindgen('./analysis-worker_bg.wasm');` (either quote style, optional whitespace and `./` accepted).
- `analysis-worker.js`, a UTF-8 classic `let wasm_bindgen` glue file containing the WebAssembly binding.
- `analysis-worker_bg.wasm`, with the WebAssembly version 1 header.
- `index.html`, containing one `meta[name=collagene-worker]` with the **entire quoted content value** `__COLLAGENE_WORKER_URL__`. A sentinel in a script/comment is not accepted. Duplicate sentinels or worker meta elements fail.

The packager domain-separates BLAKE3 with `collagene-worker-bundle-v1\0`, then hashes each filename's u64 little-endian byte length/name, its u64 byte length, and exact bytes in the fixed loader/glue/WASM order. This is a content-integrity address, not an authenticity or JavaScript safety signature. Glue and WASM checks diagnose unsupported shapes, not full JavaScript parsing or WebAssembly validation; Trunk/rustc own those build semantics.

All three unchanged byte streams are written into `workers/<64lowerhex>/`. The loader's relative imports consequently stay inside that same immutable directory. The HTML meta receives the public-prefix versioned loader URL, escaping `&`, both quotes, `<` and `>`. The entire raw URL is limited to **4096 bytes**, matching the browser worker factory. Directory prefixes may be root-relative, relative or HTTP(S); query/fragment, credentials, backslash and parent-directory traversal are rejected. The public prefix is not fetched or authenticated.

Root worker copies remain unchanged for compatibility, but the application must use the versioned metadata URL exclusively. Older immutable worker directories are retained; pruning is a separate deployment/retention policy and is not implemented here.

A temporary complete directory is renamed into place before the staged index is atomically replaced. Existing hash directories must contain exactly the three expected regular files with identical contents; corruption is an error, never repaired in place. Failed HTML validation leaves the input unchanged. An index-write failure can leave an unused complete immutable bundle; it cannot publish a partial three-file directory or update only part of the index. Concurrent cooperative packagers for identical input can reuse a verified bundle.

The packager appends one `collagene-worker-bundle:v1` HTML comment. This permits idempotent repeated packaging and rebuilding an already packaged stage after a worker file changes, while a fresh index missing its sentinel still fails. It is an idempotence marker, not a cryptographic assertion about arbitrary HTML.

Byte limits are 8 MiB HTML, 1 MiB loader, 16 MiB glue, 64 MiB WASM and 80 MiB total worker contents. Reads also check file length and modification metadata before/after. The explicitly chosen staging directory is canonicalized; its own final component and every managed input/output directory/file reject symlinks and Windows reparse points. Ancestors outside the explicitly chosen staging boundary may be canonicalized (including normal platform temporary-directory aliases). Run after staging finishes in a build-owned tree without hostile concurrent filesystem mutation: portable `std` path checks do not supply openat-style isolation against an actively racing filesystem attacker.

The CLI requires `TRUNK_STAGING_DIR`, takes `TRUNK_PUBLIC_URL` (default `/`), packages the files, and prints the resulting loader URL. Errors are nonzero and avoid printing input bytes or credentials. With Trunk, run it in a `post_build` hook after both main and worker staging outputs exist.

The library and CLI have only a BLAKE3 dependency. Tests use synthetic valid relative-loader/glue/WASM fixtures, check every file's hash influence and exact relative chain, both HTML quote styles and escaping, idempotence, changed inputs, corrupt existing bundles, missing/invalid input, metadata/size limits, CLI environment behavior and symlink rejection. No browser or external scientific program is used by the library.
