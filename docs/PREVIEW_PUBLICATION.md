# Publishing the read-only development preview

The preview serves an immutable compiled asset snapshot from its explicit
`--root` folder, read at startup. Requests use the in-memory snapshot. There is no
host filesystem API, upload endpoint, directory listing or application disk-write
route. Browser-selected sequence files stay in that visitor's browser.

Use a successful CI frontend artifact and record its exact source commit and
artifact identity. Do not start a local release/Trunk build on beagle for each
preview update. Keep the currently served directory while preparing the next
snapshot. Extract only bounded regular artifact files inside a fresh directory;
reject absolute/traversal paths, links and special files.

Stage the current compiled distribution into a **new** directory:

```sh
python3 scripts/stage_preview.py \
  --build /tmp/collagene2-preview-build \
  --output /tmp/collagene2-preview-next \
  --manifest /tmp/collagene2-preview-next-validation.json
```

Choose unused output and manifest paths. A snapshot contains exactly one current
versioned analysis-worker bundle. Staging rejects additional worker versions;
there is no merge or older-worker retention mode. Existing tabs may need reload
after publication. Keep the previous snapshot outside the active serving root
for manual rollback.

Staging always validates local script/style/preload references, all declared
SHA256/384/512 integrity digests, one active worker meta reference, its complete
versioned worker bundle and WebAssembly binary headers. Validation runs before
creating the destination. It does not execute JavaScript, start a browser/server
or make HTTP requests. These checks do not replace a successful build or establish
interactive acceptance. The manifest records asset count, byte count, index digest
and worker identity; build source provenance comes from the recorded CI artifact.

Check the CI build result and static evidence before replacing the listener.
Verify the old process's exact command and root before stopping it. Start the
replacement with:

```sh
python3 -u scripts/serve_preview.py \
  --root /tmp/collagene2-preview-next \
  --bind 0.0.0.0 \
  --port 8082
```

Use only the staged compiled distribution as `--root`. Confirm the new listener
through socket/process identity. Keep private services and their storage untouched.
Browser and application HTTP checks remain deferred until all in-scope code is
complete, as requested by the user.

The staging tests are entirely offline:

```sh
PYTHONPATH=scripts python3 -m unittest scripts/test_stage_preview.py
```

They cover a fresh current bundle, invalid references/integrity/WASM/worker assets,
additional worker-version refusal and refusal to overwrite an existing snapshot.
This command does not run the separate HTTP server tests.
