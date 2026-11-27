# Per-construct internal contribution edits

Each selected library construct declares its own required `junction_edits`,
`substitutions`, and `edits` vectors in `TargetDeclaration::OriginalRegions`.
Empty vectors explicitly declare no change. Literal expected and replacement DNA
uses the original **oriented contribution** coordinates, before junction trims,
substitutions, or length-changing edits. Native intent validation checks exact
expected letters, conflicts and retained intervals; it does not reinterpret
coordinates against an already edited target.

The target is frozen from original sources and these declarations before physical
preparation. The resulting `IndependentIntent` retains the actual native region,
junction, substitution or general-edit report. Editing a row does not alter
another construct that selected the same source, and does not modify upstream
biological documents. Existing retained targets remain a distinct explicit target
choice and never receive hidden edit declarations.

Preparations stay explicit. Keeping the old PCR or ligation recipe may produce a
target mismatch. A literal edit is neither proof that a primer search can meet its
constraints nor permission to adopt a simulated product as the intended target.
The existing edited-primer planner can propose ordinary PCR and reciprocal
extension-PCR preparations for a selected overlap construct, followed by native
whole-recipe verification. Applying a proposal updates that construct's existing
preparation overrides; it is not a library-wide automatic operation.

The planner uses the intent's single explicit per-input material-binding roster,
with exact events and `MaterialOrigin::{Run,Derived}`. Saved-template proposals
require the exact retained workspace; no second role-origin list is retained.
Borrowed input, request and aggregate mixed-origin limits precede source copying;
native ownership attestation checks the whole retained molecule/event before
primer search. Each source preparation is emitted as `FromRun` or `FromDerived`,
so the final recipe verifier reconstructs actual native strands and modifications.
Native source views feed one nonverifying primer-compilation phase. Per-row
ranked calls and search notes are diagnostics; the composed recipe has one fresh
scientific assessment. No local original-reference status stands in for actual
material verification.
`Report::validate_sources` rechecks exact source/workspace correspondence without
repeating primer ranking or PCR. Read-only `validate_binding` checks documentary
report correspondence and does not grant native material authority.

Library resolution likewise validates selected saved-origin events against each
choice's exact event and attests all selected origins through the shared native
ledger before cloning physical inputs. An equal-DNA original stock, another event,
or an unretained capsule cannot substitute for that choice. Library previews,
batch capture/replay and release use the same resolver and canonical recipes.
A stale literal in one bounded batch row becomes that row's error while other
rows continue; malformed allocation limits reject the request before evaluation.

Limits remain 12 bins, 256 choices/selected rows, at most 256 substitutions and
256 edits per construct, 100,000 total expected/replacement edit bytes, and the
existing bounded row/request/native-origin envelopes. Physical recipe constraints
and all source-event, target, intended-path and search-limit findings remain in
force. Internal declaration identity follows the current source-selector and batch contracts;
there is one current representation, with no old-shape reader or conversion.

Aggregate batch reservations use the same declared preparation-cost calculation
as native recipe preflight, including overlap-extension and saved Run/Derived
children. Reserving work neither clones a closed origin nor runs PCR. Exhausting
the bounded batch budget leaves later rows explicitly uncomputed; it cannot
silently omit their expensive preparations or report scientific acceptance.

The independent `library_internal_edits` fixtures exercise a real target-free
Derived PCR source through row-specific extension-PCR proposal, batch restart,
release and fresh saved-run replay, plus exact ownership, stale literals,
substitution/indel conflicts and allocation refusal. These are authored native
fixtures; validation is performed by the parent integration batch.

Whole-primer planning retains the selected finishing policy in its request and
single composed final preview. Reference-only,5′ fill/seal and3′ unsealed annealing
are explicit alternatives sharing the existing recipe compiler. Planned host
rescue is documentary and cannot promote the final reference to physical material.
The physical terminal comparison allowance is reserved before primer search in
the existing whole-request work budget.
