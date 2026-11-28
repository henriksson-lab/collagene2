# Typed history admission review

2026-11-28. Current source fixes and adversarial fixtures pass selected native/WASM and connected execution checks. [Combined compiler and execution evidence](causal-history-current-integration.json).

| Confirmed issue | Current correction | Evidence to execute |
| --- | --- | --- |
| Supplemental source events were deduplicated after applying edits, before the parent-count cap. Programmatic input could cause quadratic copying outside the declared bound. | Parent count and exact bounded event-ID shape are admitted before action execution or deduplication. The shared analysis commit also checks them before cloning parent IDs. | A request with 129 parents and an independently invalid edit must return `TooLarge` first and retain the original head and event count. |
| Typed source ownership validation scanned the complete history for every reference/fork event and paste donor. | Validation reuses the canonical event index already built by workspace admission. Source and donor lookups use that index; attachment membership uses its existing map. | Existing source-ownership and rehashed-tamper cases retain their rejection semantics. This is an algorithmic source finding, not a measured latency claim. |
| An action array could fit its own 32 MiB cap while its required transition envelope exceeded the native archive entry limit. History retention could succeed before later export refused it. | `RevisionEvent::check_fields` admits the complete typed transition under the same 32 MiB archive-entry cap after iterative action admission and before event hashing. New, fork, derive, commit and load use this canonical check. | An exact-limit action array with an oversized transition envelope must be refused by event construction without retaining an event. |

The native operation owner separately admits action count, borrowed model shape,
aggregate literal/report bytes and transaction copy work. Native workspace
admission bounds the event collection, descriptors, expanded transition payloads,
source models and attachments. These checks do not grant scientific authority:
causal reconstruction still needs exact event bindings and fresh supported native
recipes. No obsolete internal decoder or alternate commit implementation was
introduced to resolve these findings.
