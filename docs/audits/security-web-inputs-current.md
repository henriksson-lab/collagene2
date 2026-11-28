# Web and untrusted-input security audit

The image/resource-error follow-up passes302 focused cases and native/WASM checks.
`with_native_frame_budget` now returns an explicit `NativeFrameOutcome`: typed
resource admission is carried separately from Serde diagnostics. Native decoding
restores `TooLarge` for quota failures and preserves `Format` for malformed data,
without matching error strings. Existing contig/CRISPR/methylation/codon assertions
remain unchanged and pass. [Current evidence](image-evidence-current-integration.json).

Integrated verification: native/WASM compilation and the connected **490-case**
focused batch pass, including the relevant admission/owner/search regressions.
[Exact case list and scope](features-security-current-integration.json).
Application HTTP, browser and live-provider acceptance remain unexecuted.

28 November 2026. Source review of the working tree based on published
`f3d39cf`, including the connected PNG, hosted-MCP and model-admission changes.
This is a bounded security review, not a claim that the application has no other
vulnerabilities. No browser, CDP, screenshot, application HTTP request, network
attack, Cargo build or test execution was performed. Production edits below were
explicitly assigned by the integrating parent after the finding was reported.

## SW-01 — native model collection amplification before validation

**Confirmed availability/resource-admission defect. Source correction complete;
compiler checks and focused regressions pass in the integrated batch.** A malicious native file
can allocate large owned collections before the existing semantic resource checks
run. This requires importing attacker-controlled data (or an authenticated MCP
client submitting it). No unauthenticated network exploit or code execution is
claimed.

The relevant paths before correction were:

- `core::protocol::load_document` called generic Postcard decode, constructing a
  full `Molecule`, and only then `Molecule::validate`.
- `Molecule.features`, `Feature.qualifiers`, `Molecule.metadata`, and `Location`
  `Join`/`Order` children used unbounded derived collection decoding. Location
  already guarded recursive depth at 32; that did not limit collection width.
- `Molecule::validate` rejected more than 100,000 features, but only after their
  construction. It did not apply the native qualifier/map/text aggregate budget.
- `web/src/io/controller.rs` called native `.cg2` decoding in the browser's main
  task after reading the file. `.cg2design` snapshots and `.cg2candidate` plans
  similarly reached native decoding. Native MCP `sequence_exchange_tools` uses
  the same document loader; trusted native desktop embedding also exposes it.
- The 32 MiB frame bound and correct BLAKE3 checksum limit encoded bytes, not
  allocated object counts. A malicious author can compute the public checksum.
  Postcard's cautious initial collection allocation does not cap the subsequent
  element loop.

A minimal feature occupies **8 Postcard bytes**: empty kind length, a Range
variant with tiny coordinates/flags, empty qualifier count and a false review
flag. A qualifier `("", None)` occupies **2 bytes**. Large complete collections
therefore fit within the encoded frame cap while allocating many larger owned
Rust elements. Qualifier-heavy models could even pass the old shallow
`Molecule::validate`. Exact memory consumption depends on target layout; no OOM
measurement or browser crash was executed during this audit.

### Current correction

`model/src/decode.rs` now owns bounded current-format decoding, using explicit
per-molecule seeds for features, locations, metadata and text. It preserves the
current Serde field/variant order; there is no old decoder, migration or alternate
model. Standalone Feature and Location decoding is bounded too. Known oversized
collection hints fail before reading elements; unknown hints stop before the
first element beyond the cap. No untrusted length is used for a Vec reservation.

The existing annotation-edit budget remains the authority: 100,000 features,
1,000,000 total qualifier tuples (16,384 per feature), 1,000,000 location nodes,
16,384 Join/Order parts, 32 levels, 4,096 metadata entries, 8 MiB annotation text
and a separate 32 MiB retained-original allowance. Counts and text are charged
before constructing their corresponding owned fields.

A synchronous scoped frame guard shares an explicit envelope-lifecycle budget:
1,000,000 each of features, nodes, qualifier tuples and metadata entries, with
annotation plus retained-source text sharing the admitted envelope maximum
(32 MiB for ordinary frames; up to 128 MiB for the current workspace profile).
Per-molecule limits still apply independently. This closes cross-molecule
amplification without imposing one molecule's 4,096 metadata/32 MiB originals
allowance on a complete 1,000-document workspace. The private non-Send guard
never crosses an await; nested frames share the outer counter, thread-local state
is isolated between threads, and RAII restores the scope on error or panic.
No mutable RefCell borrow is held across recursive deserialization.

Exact retained-evidence metadata keys use the separate original-record budget:
`trace_original_base64`, `source:snapgene:original_base64`,
`legacy:original_base64`, `fastq:original_base64`, `fastq:quality_base64`, and
`annotation_import:<64 lowercase hex digits>:original_base64`. Arbitrary
`*_base64` keys do not gain an exception. This preserves intentional current
ABIF/SCF, SnapGene, original Collagene, FASTQ and annotation-source attachments.
The borrowed annotation-edit guard uses the same classification. Ordinary
scientific computation keeps its existing smaller admission profile.

Eleven active unit fixtures pass in the parent integration: claimed/unknown
collection lengths, remaining aggregate qualifier budget, 100,000-feature current
round trip, separately budgeted original text, fanout/depth, standalone fields,
cross-molecule aggregate rejection, 1,000-document metadata support, multiple
retained sources and combined frame-text limits, large exact-key attachments and
non-exempt arbitrary keys, scope error/panic/nesting/thread isolation, and literal
current Postcard records independent of the serializer. Attack-sized allocation
is avoided by precharging private counters where appropriate.

This scope wraps native envelopes, not every unrelated JSON container in the
program. Standalone Molecule/Feature/Location deserialization retains local limits;
an arbitrary outer non-envelope collection must still own its count, frame and
aggregate admission policy. The parent is reviewing current plan/worker DTOs for
those enclosing collection boundaries. Text checks occur before copying borrowed
Postcard strings; a generic deserializer may construct its own temporary string
before calling a visitor, so this does not claim to cap arbitrary parser internals.

### Save/load integrity symmetry

`protocol::save_document` now applies borrowed annotation-edit admission before
validation/serialization, so molecularly valid but oversized qualifiers, location
fanout, metadata or text cannot be emitted as a CG2 file the current decoder
cannot reopen. `bounds::native_frame_preflight` applies each model's independent
native budget and the enclosing frame budget without cloning or serializing.
PlanningWorkspace validation uses this before revision hashing; Candidate creation
and validation apply it to their input models before computing plan revisions.
Capsule origin preflight remains authoritative for its stricter shared nested
model budget. Its no-origin fast path skips ordinary inputs, so it could not
replace the new direct-input native admission.

Two additional integration fixtures were authored in the existing registered
workspace/transfer test module: illegal in-memory model rejection across CG2/CGD/
CGP saves, and successful native 100,000-feature plus retained-evidence, 1,000-doc
workspace and ordinary candidate round trips. They remain unexecuted until the
integrating parent's focused batch. No save-time roundtrip decoding is used.

### Typed resource rejection through Serde

The native admission scope now returns `NativeFrameOutcome`, retaining a typed
resource-exhaustion flag alongside the closure result. Model-budget TooLarge,
sequence length admission, collection hints and unknown-length overflow,
location fanout/depth, and bounded binary byte admission mark this scoped status.
The envelope maps a failed decode with that status to `Error::TooLarge`; malformed
non-resource input remains `Error::Format`. This preserves scientific source
bounds diagnostics when the bounded decoder rejects before later source checks.
There is no matching of Serde error-message text. Archive loading consumes the
same outcome while preserving its other existing typed errors. Nested scopes
share admission state; errors/panics restore it and independent threads isolate it.

Two focused model cases and one envelope case were added for typed classification,
including oversized versus invalid sequences, source metadata versus malformed
input, bulk bytes, unknown lengths and nested/panic/thread reset. They are source
complete and await the integrating parent's existing combined recheck; no separate
Cargo/build/test invocation was run by this worker.

## SW-02 — native candidate/snapshot list allocation before plan admission

**Confirmed same-class availability defect. Parent owns the source correction and
its verification.** `transfer::load_candidate` decoded a complete Candidate before
its plan's 1,000-input/1,000-step limits were checked. A minimal
`Operation::Overlap { fragments: [], overlaps: [], circular: false }` occupies
4 encoded bytes. Roughly one million such steps fit in a 4 MiB candidate body,
while constructing one million much larger Operation enums before rejection.
The attack does not need a scientifically valid plan because rejection was late.

`PlanningWorkspace.documents` was likewise only checked after full decode. Its
models could also split an aggregate model budget across individually legal
members. The parent applied the reusable `decode::bounded_vec` resource
visitor to the current DesignPlan, snapshot and operation source collections,
using the existing canonical limits. No encoding version change is needed for
these unchanged valid structures. A fixture can use a valid frame prefix with an
oversized declared collection and assert early resource rejection without building
an attack-sized collection.

## Reviewed guards and boundaries retained

- The web's raw `inner_html` sink found in `restriction/gel.rs` receives only
  `gel::render_svg`, not arbitrary imported SVG. Its renderer escapes lane names,
  preparations, revisions and serialized metadata, validates numerical fields and
  limits lanes/fragments. The typed map SVG renderer also escapes text and uses
  fixed styles/colors. The new raster owner rejects external resource references;
  its implementation and independent fixtures belong to the PNG worker.
- Protocol HTML uses `escape_html`, including the source-label closure; alignment,
  layout and quantity HTML exporters likewise escape text. No confirmed XSS path
  was found in these reviewed sinks. This is not exhaustive coverage of every UI
  callback or future renderer.
- Reviewed profile/citation href sources require `https://` before use. Leptos
  text and attribute bindings are distinct from raw markup. Browser download hrefs
  are generated Blob URLs.
- Browser library credentials are session-memory state. Requests use same-origin
  mode, omitted ambient credentials and redirect rejection. Disconnect aborts
  outstanding requests and advances an epoch checked before response publication.
- Worker clients admit a borrowed envelope prefix before copying complete reply
  bytes; current schema, length, revision, checksum and request ID are checked.
  Component controllers must still guard source/event/generation and disposal.
  This review did not establish a new exploitable stale-result path there.
- Workspace ZIP loading validates the directory/entry structure and declared
  aggregate expansion, restricts entry kinds/compression and uses bounded reads.
  XML importers reviewed use explicit node/depth limits and disallow DTD/entity
  resolution. No ZIP extraction or ambient filesystem authority was introduced.

Parent integration should update the single current feature/security plan with
these findings and their source/check status. No extra generic acceptance project
is proposed; the user's finished-product review remains the main usability gate.
