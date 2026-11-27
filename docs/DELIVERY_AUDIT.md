# Current delivery gap audit

This is a source-based triage of the old remaining-work list, not a completion
certificate. The full feature matrix and M0–M6 gates remain in
[the current plan](PROJECT_PLAN.md). Browser/application HTTP checks remain deferred
until in-scope implementation is complete.

| Area | Current source evidence | Remaining delivery work |
| --- | --- | --- |
| Protected codon ranges | `codon_design::ConstrainedPolicy` includes protected ranges and allowed codons; its independent checker reports violations | Current platform/workflow acceptance; the older blanket “protected ranges missing” entry is obsolete |
| Plate layouts | `protocol_layout` plus web/MCP capture and CSV/HTML/PDF export | Broader workflow acceptance; this is documentary cloning layout, not inventory or instruments |
| Lab comments/review | Native collaboration store, individual member registry, discussion SDK and shared browser/Tauri remote controls | Exact-source CI, graphical acceptance and deliberate deployment cutover; no automatic offline comment queue |
| F33 collection similarity | `sequence_similarity` plus binary worker114, event-bound Compare controls and immutable-reference MCP discovery rank explicitly selected DNA/RNA sources | Current-source graphical/platform acceptance; server-wide/indexed similarity is outside this implemented profile |
| F34 specialized biology | The existing protein profile explicitly models an unmodified free chain and reports unavailable modified quantities | RNA modification semantics, antibody-specific analysis and imported structure visualization are not established by that profile; keep them in the full plan after cloning work |
| M4 source interoperability | URL adapters and retained-source parsers exist | Real provider/interoperability evidence remains separate from mocked fixtures |
| M6 installed desktop | Native CI packaging exists for Windows/macOS/Linux | Actual installation, offline graphical use and chosen-remote workflow acceptance on the supported platforms |
| Viewer/performance/accessibility | Existing GPU/layout sources and historical measurements exist | Current-source browser/hardware evidence and the user's visual review; historical screenshots/tests do not prove current acceptance |

This audit found no reason to reopen deferred BAM/BAI work or external Primer3 /
multiple-alignment engine ports. OpenGel, inventory/ELN, instruments, scheduling,
organizational reporting and mobile expansion retain their exclusions/deferrals.
