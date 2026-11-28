# Frontend bundle size correction

2026-11-28, source-only follow-up to desktop run 46165940229 at 39d70f4.
The captured log `/tmp/collagene2-ci-46165940229-clean.log` records a successful
Trunk `dev` unoptimized WASM build followed by the native packaging hook's
`worker build input exceeds its byte limit`. That old diagnostic did not retain
the affected input or observed size. It does **not** establish which asset failed
or whether the individual-file or complete-worker limit was exceeded.

The packager already streams compiler WASM inputs up to 1 GiB when stripping
optional names. Its final bounds remain 8 MiB HTML, 1 MiB loader, 16 MiB JS,
64 MiB worker WASM and 80 MiB complete worker content. A size refusal now retains
the exact managed input name, observed bytes and admitted limit. Staging-entry
refusals separately report entry counts. Existing adversarial sparse-file cases
assert all four individual input boundaries and the actual diagnostic facts.

Cargo's existing `profile.dev.package."*"` optimization override applies to
non-workspace dependencies; workspace core/UI/worker packages remained at default
development optimization 0. This behavior is explicit in the
[Cargo profile documentation](https://doc.rust-lang.org/cargo/reference/profiles.html#overrides).

The dedicated `frontend` profile inherits `dev`, with optimization 1 and debug
information disabled. Trunk selects it for both Rust assets. The pinned
[Trunk 0.21.14 build configuration](https://github.com/trunk-rs/trunk/blob/v0.21.14/src/config/models/build.rs)
supports `cargo_profile`; its
[Rust pipeline](https://github.com/trunk-rs/trunk/blob/v0.21.14/src/pipelines/rust/mod.rs)
passes the chosen profile into Cargo. Both desktop frontend workflows explicitly
select the same profile. Ordinary native/WASM Cargo development, check and test
profiles remain unchanged; the native post-build helper stays in `dev`. Existing
Cargo target/cache locations remain unchanged. There is no extra release build,
LTO requirement or wasm-opt invocation.

This correction optimizes the actual delivered development frontend rather than
raising the application/packaging bounds. The optimized full frontend size and
CI success still need the next coordinated delivery build; no local Cargo build,
test, release verification, browser or application HTTP check was run here.

Native packager compilation and its focused fixtures pass in the connected causal/history batch. The optimized real frontend and requested desktop packages still need exact-head CI; no bundle size or installer acceptance is inferred from compiler checks.
