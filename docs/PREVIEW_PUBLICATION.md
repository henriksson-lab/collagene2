# Publishing the read-only development preview

The preview serves an immutable compiled asset snapshot. All served content comes
from its explicit `--root` folder, read at startup; requests use the in-memory
snapshot. There is no
host filesystem API, upload endpoint, directory listing or application disk-write
route. Browser-selected sequence files stay in that visitor's browser.

Build with the existing project toolchain. On beagle, reuse the established lean
Cargo target and profile settings; do not create another target/cache merely for
publication. Keep the currently served directory while producing the next build.

Stage into a **new** directory and preserve the current snapshot's versioned
analysis workers. Example after a successful build:

```sh
python3 scripts/stage_preview.py \
  --build /tmp/collagene2-preview-build \
  --previous /tmp/collagene2-viewer-assembly-public-dist \
  --output /tmp/collagene2-preview-next \
  --verify-assets \
  --manifest /tmp/collagene2-preview-next-validation.json
```

Choose the actual current snapshot for `--previous` and unused output/manifest
paths. Additional `--previous` arguments are supported. Staging drops obsolete
main bundles but retains old worker bytes so open tabs can finish their work.
Conflicting bytes at the same versioned worker path are rejected.

`--verify-assets` checks the index's local script/style/preload references and
all declared SHA256/384/512 integrity digests, one active worker meta reference,
complete active/retained versioned workers, and WebAssembly binary headers.
Validation runs before creating the destination. It does not execute JavaScript,
start a browser/server or make HTTP requests. Header/reference checks do not
replace a successful compiler/build result or establish interactive acceptance.
The manifest records the actual asset count, byte count, index digest and worker
identities; it makes no claim that an arbitrary directory was built from a given
source revision. Record the build's source commit separately.

For this development preview, check the build transcript, source revision and
static evidence before replacing the listener. Verify the old process's exact
command and root before stopping it. Start the replacement with:

```sh
python3 -u scripts/serve_preview.py \
  --root /tmp/collagene2-preview-next \
  --bind 0.0.0.0 \
  --port 8082
```

Use a compiled distribution directory as `--root`. Keep the previous served
snapshot as rollback until the new listener is confirmed. Socket/process identity
can be checked without an application HTTP/browser smoke test. Browser checks
remain deferred under the user's instruction.

The staging tests are entirely offline:

```sh
PYTHONPATH=scripts python3 -m unittest scripts/test_stage_preview.py
```

They cover fresh assets and exact old-worker preservation, invalid references/
integrity/WASM/worker bundles before destination creation, worker identity
conflicts and refusal to overwrite an existing snapshot. This command does not
run the separate HTTP server tests.
