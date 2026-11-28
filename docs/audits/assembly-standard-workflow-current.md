# Assembly-standard workflow: current implementation

Current status, 2026-11-28. The registered F37 workflow passes the selected native
and WASM checks in the combined F20/F28/F37 batch, including corrected native
archive declaration retention and real SDK round trips.
[Execution scope and source stages](connected-design-current-integration.json).
Browser/application and final user acceptance remain unperformed.

## What is connected

The assembly wizard now has a collapsed **Declared part standard (advisory)**
section. The operator chooses an installed profile and explicitly binds it to the
current stable source roles and exact retained events. Each part needs a declared
grammar slot and explicit zero-based body start/exclusive end within its selected
oriented source interval. Binding a profile leaves those fields blank; gene
names, feature labels, primers, tails and clamps do not supply their values.

The current native contract is `DesignPlan.standard: Option<Declaration>`, a
required-current nullable field on the single canonical plan. The declaration
contains an installed profile ID/identity, ordered role/input/revision/interval/
orientation/slot/body descriptors, and native assessment limits. It contains no
new molecule copies, executable DNA operation or alternate plan. The native
standards owner resolves the plan's actual original input snapshots, checks exact
revisions, and calls the existing `standards::Request`/`assess` implementation.

`VerificationReport.standard_review` is a separate required-current nullable
advisory field. Its closed state is `Assessed(existing standards Report)` or
`Unavailable(typed reasons)`. Findings remain outside the status-changing cloning
check vector. A source-format violation does not remove a product, change
chemical authority, alter target agreement, add a cloning step or change overall
cloning status. A scoped format pass is not full standard/vector compliance.

## Installed profile scope

| Installed profile | Actual checked scope | Limits of the claim |
|---|---|---|
| RFC10 2007 part sequence | Explicit formatted prefixes/suffixes, body constraints, required/advisory avoided sites, installed coding rules | No biological order, complete vector/registry compliance or expression/physical success claim |
| RFC1000 Level-0 documented body | Nonempty contiguous subpaths of four declared roles, empty additions, documentary fusion syntax, installed CDS-body prefix rule | Not a required complete transcription unit, physical termini, complete frame or vector compliance |
| MoClo 2011 Level-0 | Installed complete four-module cytosolic or five-module signal-peptide path; explicit bodies, empty additions, original scoped entry/release site rules | No inferred targeting/function/translation/frame, automatic ATG/GGT addition, destination vector or higher-level compliance |

The installed definitions retain exact catalogue pins, data revisions, source
citations, rights and available capture hashes. Their primary sources are
[Knight's RFC10 text](https://openwetware.org/mediawiki/index.php?title=The_BioBricks_Foundation:BBFRFC10&oldid=297083),
[iGEM assembly documentation](https://technology.igem.org/technologies/assembly/),
and [Weber et al., 2011, Figure 1B/C](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0016765).
This integration adds no biological parameter table, copied application code,
external runtime executable or fresh upstream-source verification claim.

RFC10 formatted flanks are actual selected sequence. RFC1000/MoClo body profiles
require their explicit empty additions; PCR recognition sites/clamps/spacers/tails
are neither silently stripped nor inserted to obtain a favorable result.

## Source ownership, ordering and editable state

`MethodsState.standard` stores a `StandardDraft` inside the existing current
assembly draft. It retains profile identity and per-role exact source event,
selector, preparation kind/orientation, slot and body text. Source-role IDs are
made persistent when the operator binds the profile. Actual reorder operations
retain those identities and rebuild declared order from the current role order;
the existing native grammar checker evaluates that order.

Changing the exact source event, selector, orientation or preparation makes the
binding stale. Editable stale settings remain stored and visible, rather than
being cleared or reinterpreted. Execution/review guards require explicit rebind
or clearing of that declaration. Rebinding is an explicit action that clears
prior slot/body values. Scientific nonconformity in a fresh native review remains
advisory; an unfinished/stale editable declaration is a separate workflow issue.

Current support is original-reference input intervals only. Saved/prepared strand
coordinates are not converted by length, DNA equality or labels. The UI refuses
to bind a material-origin selector as an original part, and native physical-marker
inputs produce typed unavailable advice without dropping physical products.
Future prepared-part support needs a genuine native reference-coordinate receipt.
The independent source-format review does not infer the final prepared fragment's
ends, inserted PCR sequence or biological identity from its declared role.

Native overlap, ligation and tailing recipes and Golden Gate reaction requests
carry the optional declaration into the plan before its final hash/verification.
Combinatorial declarations belong to each `SelectedConstruct`, not the shared
source-choice catalogue: different selected alternatives cannot silently inherit
another choice's body coordinates or source pins. Proposal-only request types
that have no declaration explicitly construct `None`; a wizard review guard
refuses a candidate that omits its current declared standard and requires a fresh
method preview.

## Persistence and protocol output

The editable draft payload is version 24. Named designs retain the current typed
settings and exact standard event pins; library recovery inventories also retain
working and inactive construct declarations. Source pin attestation checks actual
retained events, not equal DNA or current heads. Restoring settings does not
restore a scientific review or ready state. Library impact presentation explains
retained/stale standard bindings without rewriting their body offsets.

The saved-design wrapper stays `assembly-design-v1`/schema 1 and library recovery
stays schema 3 because their outer layouts did not change; embedded draft admission
is strict current-version admission. There is no draft migration or fallback
reader. The current cloning RunRecord is **schema 2**, because its recorded
assessment now retains the typed advisory review. Its attachment category and
separate run-capsule wrapper remain version 1; those are different contracts.

HTML and PDF render the same fresh `CheckedProtocol` verification review through
`output::standard_review`; they do not run a second checker per renderer. Output
includes profile/data/catalogue/declaration/request/report pins, explicit part/body
bindings, native outcomes/completeness/budgets, findings, citations and scope.
MCP compact verification exposes an explicit nullable advisory summary and the
full artifact. Web method/source views paginate native findings and use the
existing native `SourceBinding::map_interval` for displayed original coordinates.
Plan-aware response consumers invoke `validate_standard_for_plan`, preventing a
report with absent advice from being accepted for a declaration-bearing plan.

## Bounds and refusal semantics

The declaration is limited to 128 descriptors/128 KiB. The wizard has at most 12
rows; its stored field text and event identifiers are bounded. Native request
admission is 16 MiB, enclosing plan admission 32 MiB, and the native review wrapper
reserves the existing 2 MiB report envelope plus 64 KiB for review metadata. Native
assessment limits also apply: at most 128 parts, two million selected bases,
64 million sequence comparisons, 4,096 findings and the library's report limit.

`DeclarationBudget` is shared in native replay scopes. `charged()` includes
explicit reservations, whereas `observed_comparisons()` contains measured counts
from successfully retained, validated native reports. The effective native limit
and remaining-before value are explicit. A native error that returns no comparison
counter retains a `ReservedOnNativeRefusal` charge; it is not reported as measured
work. Settlement happens after final review association/hash/size checks.
Admitted evaluation/retention failures become compact typed unavailable advice,
with reservations retained. Malformed/unrepresentable declaration admission still
refuses explicitly; it is not replaced by an absent declaration or physical pass.

The existing aggregate draft, named-design, workspace, recovery and source-union
limits remain in force. Bounded combinations of many inactive/stale declarations
may be refused by those aggregate limits; there is no silent pin truncation.

## Authored fixtures and remaining limits

Registered native declaration fixtures cover real RFC10 additions versus body-only
profiles, RFC1000 subpaths, both complete MoClo paths/wrong order, literal reverse
circular BsaI/source mapping, profile/source/body/missing-role refusal, physical
marker refusal, ambiguity/shared exhaustion, real late native report-size refusal
and current documentary integrity. Registered wizard fixtures cover declaration
in the hashed plan, advisory failure with unchanged actual products/status,
stable-role reorder, genuine same-DNA/new-event staleness, selector/preparation
changes, explicit blank body fields, current named-design roundtrip and HTML/PDF
pins. MCP fixtures cover current transport/retention/admission paths. These registered cases pass within the shared execution record. Historical
source audits alone were insufficient: the SDK round trip caught an omitted
custom workspace-codec field, which is now fixed and independently covered.

No full standards-compliance gate, prepared/material coordinate conversion,
automatic biological slot classifier, frame/function inference, destination-vector
validation, tail/clamp synthesis or experimental-success model was added. The
standalone single-part inspector remains useful for individual source format;
whole ordered wizard review is a distinct use of the same native checker.
