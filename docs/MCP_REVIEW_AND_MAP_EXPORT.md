# MCP saved-run review and map export

The native MCP registry has189 tools. Two additions use the same deterministic
Core adapters as the web worker (worker93); they do not require new scientific,
candidate or native-plan versions. No external executables or network services
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

## export_map_svg

Request fields are exactly `id` and `revision` of an immutable session document.
The tool performs native model preflight, calls `svg_export::render`, validates the
returned source identity and bounds, and returns a complete UTF-8 `image/svg+xml`
artifact named `map.svg`. The response includes document/source pins, byte count,
all renderer warnings (natively bounded32×4096), and
`restriction_overlay: "none_requested"`.

This initial tool profile exports canonical native feature geometry only. It does
not accept raw SVG, raw site lists, paths, arbitrary URLs, hidden viewport state or
invented enzyme calls. It does not reproduce viewport rotation/zoom/selection,
feature-visibility filters, search/primer overlays or selected-site highlights.
Native renderer warnings state layout/feature-label limitations. Map output is a
reference display, not evidence of reaction activity, target success, isolation or
physical-material authority. The full native source remains unchanged.

The source-bound SVG and result metadata are prepared before the artifact is
reserved/published. Existing artifact representation and capacity are checked
under the same store lock. SVG uses the native2MiB ceiling and standard bounded
`read_artifact` chunks; there is no external rendering or filesystem write.

## Authored SDK coverage

`mcp/tests/review_map_tools.rs` uses a real in-process MCP client/server transport
and a literal saved PCR. It advances both source and independent-target heads,
checks that review still binds their exact older events, validates downloaded
full review/candidate artifacts, and explicitly calls fresh verification separately.
It compares exported SVG bytes to the native renderer, checks XML escaping and
source/warning identity, refuses stale/unknown/extra-field requests, and proves the
original workspace archive and session document list stay unchanged. The parent
runs the fixture; this source worker ran no compiler, browser, HTTP request, remote
model, external process, live-data access, commit or push.
