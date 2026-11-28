# Model, worker, draft and integration-suite refactoring verification

This checkpoint covers the current integration batch. The full product goal and
remaining organization work are still open.

## Native compilation and test discovery

The combined native harness build passed:

```sh
CARGO_TARGET_DIR=/tmp/collagene2-lean-target CARGO_INCREMENTAL=0 \
CARGO_PROFILE_DEV_DEBUG=0 CARGO_PROFILE_TEST_DEBUG=0 \
cargo test --no-run -p collagene-core -p collagene-mcp -p collagene-web \
  -p collagene-antibody-numbering --features collagene-core/json-schema \
  --offline --locked -j2 --message-format=json
```

Cargo metadata registers nine core and eight MCP suites. Actual executable
`--list --format terse` output matches every name in the
[current inventory](integration-test-inventory.json): 1,755 core cases and
180 MCP cases. All 1,931 pre-refactor cases remain discoverable; four new cases
exercise native SBOL exchange. This is discovery evidence, not a claim that all
1,935 cases were executed during this batch.

Static review found no process-wide environment or working-directory mutations
introduced by grouping tests. Existing storage fixtures keep their unique
temporary directories; optional PDF fixture writes were already in one target.

## Focused execution

Executed the built harnesses directly with two test threads, retaining one
Cargo package/feature selection and avoiding additional build cycles.

| Harness | Filter | Passed |
| --- | --- | ---: |
| MCP `protocol` | All | 9 |
| MCP `sequence` | `sbol_tools::` | 4 |
| MCP `sequence` | `sequence_exchange_tools::` | 3 |
| Web unit harness | `assembly_draft_model::` | 33 |
| Core `sequence` | `protocol` | 47 |
| Core `verification` | `protocol` | 24 |
| Core `cloning` | `protocol` | 30 |
| Core `formats` | `transport` | 8 |
| Total | | **158** |

The SBOL SDK cases cover retained original/expanded graphs, stale and foreign
source roles, bounded paging, explicit omission review, native workspace/archive
roundtrip, refusal of unsafe single-stranded projection, valid native export,
invalid identifiers, physical-evidence refusal, and source immutability.

All 203 previously published MCP tool declarations remain present; exactly five
SBOL declarations are added. Existing SDK registry assertions observe 208 tools.
The two worker enums retain their exact declarations and variant order after
normalizing whitespace and canonical model paths. Worker/candidate versions
remain 121/31. Current-format/header/selector and source-binding cases pass.

## WASM and remaining gates

The web/worker WASM check passes after fixing nine private-field reads
introduced by moving draft identity into its owner. Those callers now use
identity methods for labels and edit/material decisions. Encoding is unchanged.

```sh
CARGO_TARGET_DIR=/tmp/collagene2-lean-target CARGO_INCREMENTAL=0 \
CARGO_PROFILE_DEV_DEBUG=0 CARGO_PROFILE_TEST_DEBUG=0 \
cargo check -p collagene-web -p analysis-worker \
  --features collagene-core/json-schema --target wasm32-unknown-unknown \
  --offline --locked -j2 --message-format=json
```

Final formatting and whitespace checks pass.

The repository layout checker intentionally still fails on four violations:
core source directory (154 files), web source directory (227), MCP root (9,005
lines), and web App (7,273). Verifier (4,528) retains its review warning. These
are current work, not exemptions. Draft production code is now below 1,600 lines;
every worker protocol source is below 3,000. The
[current refactoring plan](CODE_ORGANIZATION.md) owns the remaining work.

No release build, browser, screenshot, CDP or application HTTP check was run.
