# Sanger trace import

`collagene_core::traces::import_trace(&[u8])` recognizes ABIF and SCF by their signatures. This independently written MIT implementation uses embedded Rust only. It does not run a basecaller, instrument interface or external converter.

The typed `Trace` retains called symbols and normalized DNA, zero-based called-peak sample positions, canonically labeled A/C/G/T signal vectors, typed score evidence, format metadata, warnings and unchanged original bytes with a BLAKE3 hash. Canonical labels do not imply amplitude normalization: ABIF signed short intensities remain signed `i32`; SCF unsigned intensities remain their original magnitudes.

## ABIF profile

The layout and tag meanings come from Applied Biosystems' July 2006 [ABIF File Format Specification and Sample File Schema](https://archive.gfjc.fiu.edu/workshops/resources/literature/Genetic_Analysis_and_Interpretation/06_ABIF_File_Format.pdf), a vendor document hosted by FIU.

The reader supports major version 1, bounded big-endian directories and inline values. As required for vendor-compatible readers, root directory count/offset are checked while other root hints are ignored; older libraries may reserve extra directory space. Synthetic fixtures nevertheless write the canonical root fields. Analyzed `DATA` 9–12 channels are labeled using the required `FWO_` permutation. `PBAS`/`PLOC` edition 2 is preferred; edition 1 is used only if edition 2 calls are absent. Peak positions and optional `PCON` scores must match that edition. Editions are never mixed. Integer widths and payload dimensions are checked; duplicate tags, invalid offsets, malformed channel orders and unsupported types are errors.

`PCON` bytes are retained as **raw basecaller scores**, without asserting a Phred scale. Optional sample names are read from supported ASCII string types. Uninterpreted tags remain in the original attachment. `RevC != 0` is explicitly unsupported pending an audited call/peak/trace orientation adapter. Raw-only, uncalled fragment-analysis files and alternative instrument/software schemas are not claimed as supported.

## SCF profile

The authoritative description is maintained by the [Staden project](https://staden.sourceforge.net/formats.html): [header](https://staden.sourceforge.net/manual/formats_unix_3.html), [sample representation](https://staden.sourceforge.net/manual/formats_unix_4.html), [sequence/accuracy data](https://staden.sourceforge.net/manual/formats_unix_5.html), and [section layout](https://staden.sourceforge.net/manual/formats_unix_8.html).

Versions 2.00, 2.01, 2.02 and 3.00 are supported. Version 2 samples are interleaved; version 3 samples are channel-major and recovered by two modulo-width prefix sums. Both 8-bit and 16-bit precision are supported where the version declares them; 2.00 uses its documented 8-bit/default-code profile. Sections may occur in any order, but must be bounded and non-overlapping.

Four stored accuracy arrays are preserved without conversion to Phred or a guessed called-base quality. Implemented call code sets are default 0, IUPAC 2 and ABI 4. In code set 0, `-` remains in `called_symbols` and is explicitly projected to `N` in DNA with a warning, retaining its position. Other encodings are rejected. Obsolete clip fields are retained as raw metadata and never applied automatically. Comments/private data remain in the original attachment; a supported `NAME` comment can label the view.

## Views, original evidence and limits

`Trace::reverse_complement()` explicitly complements/reverses called DNA, A↔T/C↔G channels, call positions and score arrays together. Original bytes, hash and raw source metadata remain unchanged. Reversing twice restores the view.

`Trace::molecule()` is a **called-sequence projection**, with a clear metadata note that signals are absent. Use `Trace::attached_molecule()` for native saving: it reparses original evidence, checks the derived view, stores the unchanged attachment as base64 and records orientation. The entire serialized `.cg2` document must fit 32 MiB, including base64 and native overhead. This is a bridge for the experimental single-document native format, not a complete attachment/archive design.

`original_trace_bytes()` decodes and hash-checks the attachment independently of current calls, allowing original evidence downloads after sequence edits. `trace_from_molecule()` checks the attachment hash, reparses it and reapplies the explicit view orientation. It rejects sequence/topology changes that make the current document disagree with the immutable trace. Editing a called-sequence document does not rewrite its original chromatogram. A future edited-call/derived-record model must record those changes separately.

Import limits are 32 MiB of source data, 2,000,000 samples per channel, 100,000 calls and 4,096 ABIF directory entries. Peak indexes must lie within the signal arrays and be non-descending; coincident positions are retained with a warning. These are data limits, not a claim about total process RSS. Compressed wrappers, MacBinary wrappers, unsupported uncertainty alphabets, re-basecalling, quality calibration and instrument automation are outside this parser implementation. Optional raw-score trimming and reference-free assembly are implemented separately in [the read-contig engine](CONTIGS.md), using retained call/quality evidence.

## Validation and remaining fidelity work

`cargo test -p collagene-core --test traces --locked` checks hand-authored ABIF directories/inline values, signed signals, edition isolation, four-channel labeling, SCF v2/v3 layouts, independently calculated 8-/16-bit second differences, score preservation, reverse complements, corruption/bounds failures, native attachment round trips, stale evidence rejection and attachment-size accounting. No other cloning application's reader or output is a test oracle.

The fixtures are synthetic, not a representative instrument migration corpus. Real-file interoperability, additional ABIF variants, `RevC` handling and chromatogram viewing/selection still need validation. Correct parser fields alone do not establish complete trace-workflow parity.

## Initial web chromatogram component

`crates/web/src/trace_view.rs` provides a Rust Canvas2D view for validated attached traces. The canvas backing store stays viewport-sized; native horizontal scrolling supplies pan, buttons zoom/fit, and precomputed min/max bins retain signal extrema at overview scales. All channels share one amplitude range, including negative ABIF values. Call labels are thinned at overview zoom; an accessible table shows up to 200 visible calls and their original score evidence. Click/shift-click and keyboard navigation update the workspace base selection, and selected calls scroll into view.

The component reparses original evidence and refuses stale call mappings after sequence edits. Its WASM compilation is checked. Browser painting, native scroll behavior, accessibility and performance budgets must be verified by the integrating application's browser tests; compilation alone is not a rendered or performance pass. Initial parsing and min/max preparation are bounded but currently run on the browser thread; worker-side preparation is a later improvement.
