# Construct-owned saved read checks

`read_checks` captures a fresh existing trace/FASTQ validation against an exact retained construct event and an ordered list of exact retained read events. It does not infer document identity from equal sequence contents. Switching the active viewer to a read does not change the saved construct reference.

`capture_revision` binds the complete incoming workspace fingerprint, reference/read document-event-revision bindings and validation options. `capture` resolves these events, runs `trace_validation::validate_reads` afresh, and returns a private workspace copy with an immutable `read-check-v1` JSON attachment. Existing source documents, heads, active document, annotations and history are unchanged. Invalid bindings or resource failures leave the original workspace intact.

The saved check has schema version 1 and retains the exact reference, ordered reads and full native report/options. Its own context revision identifies those event bindings and recorded engine/model/options, independently of the installation token. Thus repeating a capture still performs fresh science but identical evidence reuses its content-addressed attachment, even though the incoming workspace installation token changed. Equal-content aliases with different native events remain different saved checks. Duplicate original read files are still diagnosed by the native validator and never become independent support merely because they have different document names.

`list` returns historical summaries, including the recorded engine and recorded status. `load` checks bounded bytes, content hash, immutable event references and report consistency. Neither operation establishes a current verdict. Imported reports, including records naming another engine, are historical assertions; serializing and decoding them cannot restore fresh analysis authority.

`replay_revision` additionally binds the current workspace fingerprint, saved attachment hash and current engine/model. `replay` resolves the saved original events, independently of later source-head edits, and invokes the current validator. It returns both the original saved check and a newly computed report. `Checked::validate_binding` checks readonly source/options identity and the current engine/model pair; it does not rerun read alignment or basecalling. Hosts must compare the request token and current workspace before installing a response. The reference remains the saved construct even if the user is navigating its trace document.

Native quality scales, trimming, placement ambiguity, support requirements, uncovered calls and discrepancies retain the existing engine's meaning. A partial read remains incomplete with uncovered consensus calls `N`; reference bases are not copied into those gaps. A saved check is not physical sample identity, experimental success, an inferred correction, or permission to edit a construct.

Limits are 128 saved checks per workspace, 1–32 reads per check, at most 17 MiB saved JSON, and existing bounded models, work/consensus/column limits and workspace budgets. Options and borrowed source models are guarded before copying read models or hashing rich workspace content. Readonly checks bound report vectors, alignment runs/differences, consensus/coverage, call/reference positions, identities, engine labels and evidence-derived status consistency. They reject malformed or internally contradictory receipts without turning a historical claim into proof.

Five independent fixtures cover literal partial FASTQ coverage, repeated capture deduplication, same-content native aliases, original-file duplicate diagnostics, historical replay after user edits/navigation, archived engine assertions, current replay, tampered status/positions, wrong event binding, option limits and attachment corruption. No scientific algorithm or verification engine version changes.

## Host integration

The main viewer pins the checked construct independently of its active document.
Viewing a supporting read cannot silently turn that read into the reference.
Changing the construct head or selected read event invalidates live evidence,
including equal-content aliases. A saved replay remains bound to its historical
reference/read events. Its action guards check source-event and immutable-attachment
availability, so ordinary view navigation does not discard unchanged evidence.

Live validation pages all sites, uncovered ranges, insertions/calls and per-read
columns. Reference navigation retains exact base ranges or insertion carets; read
navigation opens current call coordinates with the stored/original mappings kept
visible. Historical sources fork without rewinding heads. Saving recomputes in a
private worker, validates the response token and live workspace, then commits current
and previous recovery states in one cancellable IndexedDB transaction. A failed or
cancelled capture installs no partial check. Positive reports are not restored by
loading native files, JSON or UI state.

Five core fixture groups, two binary worker cases and two real MCP SDK groups passed,
including literal6/12 coverage, archived-engine replay and unchanged edited source
heads. The full 125-tool stdio smoke and integrated development application/worker
build also passed; the read-only8082 preview was refreshed. Browser checks remain
deferred; source/native checks do not establish GUI
acceptance. Human review decisions and derived-reference proposals are still pending.
