# Read evidence against a reference

`core::trace_validation::validate_reads` compares 1–32 original-evidence read documents with a separately selected DNA reference. The implementation is Rust and runs in the browser worker or natively without external scientific programs. It is pairwise read placement and evidence aggregation, not a new basecaller, multiple-sequence aligner, or proof of sample identity.

## Quality and trimming

ABIF PCON values remain raw caller scores. SCF contributes the stored accuracy value for the called concrete base; an ambiguous SCF call has no invented confidence. Neither trace format is silently converted to calibrated Phred. FASTQ requires an explicitly selected Phred+33 or Phred+64 encoding and preserves that declaration in the report. Declared FASTQ scores are distinct from the raw trace scales; the software does not independently calibrate any caller.

The SCF specification describes separate per-base accuracy arrays, called symbols and sample peak positions; see the [Staden base information specification](https://staden.sourceforge.net/manual/formats_unix_5.html). Original attachments, their hashes and current called bases are rechecked. Editing calls without a valid mapped quality transformation makes the evidence unusable. A document containing competing trace and FASTQ sources is rejected rather than choosing silently.

Quality-window trimming retains the interval from the first passing window to the end of the last passing window. A window passes when the configured number of scores meets the recorded threshold in that source's units. Short windows clamp both width and required passing count. Internal low or unknown confidence calls remain visible and are individually excluded from consensus. All-calls mode disables end trimming, not confidence filtering. Reads below the minimum retained length produce no usable alignment. Missing scores default to rejection; explicit keep-unscored mode allows alignment only.

Defaults are threshold 20, width 5, four passing calls, and at least 20 retained bases. These are recorded policy choices, not universal biological cutoffs. Set a shorter retained length explicitly when inspecting short fixtures.

## Placement and coordinate proof

Every retained read is placed using the embedded affine-gap query-in-reference algorithm described in [PAIRWISE_ALIGNMENT.md](PAIRWISE_ALIGNMENT.md). Both orientations are searched by default. Concrete exact matches use a bounded exhaustive scan; a complete single mapping reaches the maximum possible substitution score. Imperfect reads also run an independently implemented rolling dynamic program which counts all optimal paths, saturated at two. A tied alignment, repeated placement, orientation tie, interrupted search or skipped orientation cannot contribute consensus evidence. Palindromic sequence orientations are distinct original-call mappings and remain ambiguous unless orientation is constrained explicitly.

Circular placement uses a doubled reference, canonical starts in its first turn, and at most one reference turn per accepted read. Queries longer than one turn and an unconstrained best path that spans more than one turn are conservatively rejected. The implementation does not claim a unique valid circular alternative when its proof cannot establish one. It does not claim whole-plasmid verification from a short origin-crossing read.

Report coordinates are zero-based. Reference positions identify actual source bases; insertion boundaries can equal the linear reference length and wrap to zero on circular references. Trim bounds are half-open in the current read. Every column retains current called-base index and original stored called-base index. Trace columns additionally retain current and original sample positions; gaps and FASTQ reads have no fabricated sample coordinates.

For an original FASTQ forward-record retained interval `[a,b)`, current call `i` maps to `a+i` in forward display or `b-1-i` after reverse complementation. The original record length and retained bounds remain explicit in each read report. Additional analysis trimming and alignment reversal compose with that mapping. A reverse-complemented display does not create an opposite biological read: strand-support accounting uses the original stored read relative to the reference.

## Coverage, discrepancies and consensus

Paired alignment coverage and eligible evidence coverage are reported separately. Only uniquely placed, quality-passing concrete calls vote. A deletion needs passing calls on both sides. Insertions and deletions obey the same minimum file-support and optional both-strand policies as substitutions. An insertion is retained as an event at a reference boundary; it is not silently incorporated into the reference-coordinate consensus.

The consensus never fills gaps from the reference. Uncovered, unknown-confidence, low-confidence, insufficiently supported and deletion sites remain `N`. Conflicting eligible concrete calls produce the corresponding IUPAC ambiguity code; a majority vote does not hide a minority call. Per-site counts, read indices, strand support and low/ambiguous/deletion observations remain available. Reference ambiguities do not become confirmed matches merely because they overlap a read call.

Original-file hashes prevent renamed copies or reversed views of the same source from counting as independent support. Different files still do not establish independent molecules or biological samples.

`Incomplete` means the required evidence has not established full agreement. `Discrepancies` means at least one sufficiently supported substitution, conflict or indel exists, even if other regions remain uncovered. `ConsistentUnderRecordedThresholds` requires every reference base to have sufficient unanimous matching evidence and every submitted read to have an error-free, complete, unique placement. None of these statuses constitutes experimental confirmation, and no imported or AI-provided verdict substitutes for computation.

The report binds reference revision, read revisions, original evidence hashes, options, engine version and model identifier. Opening a consensus creates a separate document containing the report and warning metadata; it copies no reference annotations. Keep the original read documents and reference alongside that derived document.

## Browser use and bounds

Select the independent reference, choose open trace or explicitly encoded FASTQ read documents, set trimming/support policies and run **Validate selected reads**. Results are calculated through the binary worker protocol. Changing a bound source or option invalidates the result. Clicking an ordinary reference interval selects it; a circular crossing remains represented as explicit mapped positions. JSON evidence, support-only FASTA and a separate consensus document are available. Imported reads are never automatically marked as verified.

Inputs are bounded to 32 MiB of serialized reference/read evidence, at most 32 reads and 100,000 current called bases per read. Default per-alignment work is 64 million cells and aggregate work 128 million cells; the aggregate hard maximum is 256 million. Alignment and independent path proof are both charged. Exact scanning defaults to 50 million comparisons, with a hard maximum of 200 million. Exact placement previews cap at 64. Reports cap at 200,000 columns and 16 MiB serialized output; consensus evidence sites default to 100,000, with a maximum of 200,000. Work exhaustion remains explicit and never becomes a negative experimental finding.

## Verification

`crates/core/tests/analysis/reads/trace_validation.rs` constructs independent ABIF/SCF fixtures and explicit FASTQ records. Cases cover raw confidence trimming, missing scores, partial coverage, quality-gated substitutions, conflicts, duplicates, indels, circular origin crossing, palindromes, reverse/current/original coordinates, stale evidence, work limits, and strand/support policies. Independent exhaustive substring/global-path enumeration checks scores and uniqueness claims across 288 small linear and 288 small circular cases. These tests validate the stated bounded algorithm; they do not calibrate real sequencing instruments.
