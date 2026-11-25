# Retained Golden Gate pool settings

Assembly drafts now include the editable prepared-pool/reaction settings alongside
ordinary preparation rows. These include the reaction mode, per-role substrate
choices, product topology, orientation policy and finite enumeration limits.
The existing assembly draft JSON and native-workspace attachment path carry these
settings; no new scientific receipt or history format is introduced.

Each substrate choice follows its source role when rows are reordered or removed.
A newly added role defaults to its current prepared fragment. Changing the source
or preparation keeps that explicit role choice; an unavailable whole-amplicon
choice fails inspection rather than silently choosing different material.

Older drafts omit the new optional field. They initialize with the existing
prepared-pool defaults, the draft's topology and one prepared choice per role.
Saving an untouched old draft remains valid. Saved settings must match the bounded
role count; unknown substrate values, excessive choices and oversized fields are
rejected. Empty or invalid numeric text may remain unfinished draft data; it cannot
be used for inspection until corrected. The native source-bound scientific request
still performs its own preflight.

The panel reads and edits the shared settings directly. A memo containing only
preparation/source context keeps an options edit from remounting the controls.
Receipts, candidate choices, current pages and worker handles remain transient.
Editing options invalidates the current receipt and cancels pending work. Source,
target or workspace changes invalidate the context. Navigation/reopen restores
editable settings and requires explicit inspection; it cannot restore a current
success flag or bypass fresh candidate selection.

Implementation is in `web::golden_gate_pool_form`, `web::cloning_methods` and
`web::assembly_workspace`. The pure portable form has focused native tests for
editable/runnable distinction, finite-limit parity with scientific preflight,
bounded decoding and rejection of injected evidence fields. Local browser checks
remain deferred until full in-scope implementation is complete.
