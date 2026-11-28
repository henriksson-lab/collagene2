# MCP saved-run review and map export

The native MCP registry has 209 tools. Saved-run review and map exports use the same
deterministic core owners as the current web worker (worker 126); image presentation
does not change scientific, candidate or native-plan versions. No external executables or network services
are involved.

## prepare_saved_run_review

Request fields are exactly `workspace: { id, revision }` and `run_id`. Unknown
fields and stale/unknown workspace/run IDs refuse. Both workspace IDs identify an
immutable session snapshot, not a filesystem path, desktop catalog or remote.

The tool calls `run_review::request_revision`, `prepare` and the readonly receipt
validator. Preparation resolves the retained independent target and input events
from the native run record. Current document heads may differ: they never replace
those original bindings. Preparation does not replay physics, migrate historical
material models, adopt recorded status or manufacture a target from a product.

Results include:

- `state: "candidate"`, `verified: false`, the exact plan revision and target event,
  input event pins (first100 inline with total/truncation), source-workspace and
  receipt identities, historical assessment and engine/limit warnings.
- `artifact_id`: full `saved-run-review.json`, including every input binding,
  original target, candidate and historical assessment. Inline historical finding/
  assumption previews cap20 rows of1024 characters and expose counts/truncation;
  full text remains in this artifact.
- `candidate_artifact_id`: binary `saved-run-review.cg2plan`, read through ordinary
  base64 artifact chunks. Decode concatenated chunks to recover exact bytes.
- An unverified session plan bound to the retained independent target, so the
  returned `plan_revision` can be passed directly to `verify_plan` for fresh native
  checks. Preparing a review never automatically calls that verifier.

Full JSON, candidate bytes and plan storage are staged and charged atomically
against existing item/session limits before publishing any of them. Exact repeats
reuse immutable artifacts. A conflicting representation, different exact target/
source data, serialization error or capacity failure publishes no partial result.
The uploaded workspace, its active document, heads, private branches, run records
and candidates are unchanged; the new plan lives only in the MCP session.

Historical engine agreement is metadata, not current verification. In particular,
old saved-material recipes can still be refused by a later current verifier.
Inspection/preparation does not bypass that boundary or create reusable stock.

## export_map

Request fields are exactly `document: {id, revision}` and `format`: `"Svg"` or
`{"Png":{"width":1600}}`. The tool performs native model preflight, calls
`core::output::map_export::Request::render`, validates source identity/image bounds, and
returns a complete artifact: UTF-8 `image/svg+xml` named `map.svg`, or base64
binary `image/png` named `map.png`. The source remains immutable. All bounded
renderer warnings and `restriction_overlay: "none_requested"` are included.

This tool exports canonical native feature geometry only. It accepts no raw SVG,
site lists, paths, arbitrary URLs, hidden viewport state or invented enzyme calls.
Viewport rotation/zoom/selection, feature-visibility filters, search/primer
overlays and selected-site highlights are not reproduced. Native warnings state
layout/label limitations. Rendering is documentary display; it does not prove
reaction activity, intended target success, isolation or reusable material.

PNG uses the same canonical SVG geometry, native portable Rust rasterization and
bundled Liberation Sans font. Missing glyphs are visibly escaped and reported.
Width is 16–4096 px, height preserves aspect ratio, each side is at most 4096 px and
total pixels at most 4 million. Encoded bytes are capped at 24 MiB and carry the
source revision in PNG text metadata. PNG retains pixels rather than full SVG
coordinate metadata; keep the native source for complete annotation evidence.

Source-bound bytes/result metadata are prepared before artifact reservation and
publication. Existing representation and capacity checks stay under one store
lock. SVG retains its 2 MiB ceiling; binary PNG uses the same bounded artifact
pagination mechanism. Neither format invokes an external program or grants
filesystem/network authority.

## Authored SDK coverage

`crates/mcp/tests/protocol/review_map_tools.rs` uses a real in-process MCP client/server transport
and a literal saved PCR. It advances both source and independent-target heads,
checks that review still binds their exact older events, validates downloaded
full review/candidate artifacts, and explicitly calls fresh verification separately.
It compares exported SVG bytes to the native renderer, checks XML escaping and
source/warning identity, refuses stale/unknown/extra-field requests, and proves the
original workspace archive and session document list stay unchanged. The parent
runs the fixture; this source worker ran no compiler, browser, HTTP request, remote
model, external process, live-data access, commit or push.
