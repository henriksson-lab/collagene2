# Construct QC foundation audit

Source inspection and primary-source review, 28 November 2026. The foundation is
registered in the shared native/WASM integration batch. Final compiler checks and
the focused cases pass. [Execution evidence](design-checks-current-integration.json).

## Literature correction

[Bai et al., NAR 2025, doi:10.1093/nar/gkaf697](https://academic.oup.com/nar/article/53/14/gkaf697/8222435)
reports design errors in **384 of 2,521 plasmids**. The **852** denominator belongs
to restriction-digestion validation, with **128** inconsistent patterns. The direct
sequencing cohort was **259**, with **91** variant plasmids. These cohorts must not
be collapsed into one denominator. The paper does not establish 384 independent
error instances or validate a software detector. Its taxonomy includes contextual
promoter, silencing and instability claims; describing every category as a
machine-provable sequence error is unjustified. The cohort came from a commercial
provider's clients, rather than random sampling of all laboratories.

These qualifications are reflected in the current plan and research summary.

## Owner and scientific scope

The new `core/src/construct_qc` owner has one request/options model
and one report model. It calls the existing CDS and literal translation owners;
there is no second codon decoder, plan representation or verification engine.
[NCBI's genetic-code definitions](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi)
provide codon and initiation semantics. Each used code is pinned to the current
descriptor and the existing translator's complete 64-codon ordinary output.
Context-dependent termination stays indeterminate.

Checks cover CDS length/frame, recognized start, definite terminal stop, internal
stops and stored-translation disagreement. Explicit contained fusion and junction
claims check declared register at the current coordinates. Literal reference
comparison proves only matching captured sequence; it does not authenticate a
reference's declared provenance or establish an element's function.

Annotation assumptions name exact source and feature identities and retain the
caller's rationale. They are declarations, not experimental verification. Every
feature-dependent finding labels this assumption. Imported annotations without a
current declaration, stale declarations, review-required features and unsupported
biology produce dependent indeterminate findings. Unrelated checks continue.
There is no overall pass or expression-readiness verdict.

The report preserves the required no-guarantee statement. A restored report is
documentary content; only fresh `validate_for` recomputation establishes agreement
with the current source and declared options. Hashes are integrity and binding
evidence, never biological authority.

## Resource profile

Admission uses borrowed model bounds and streaming byte measurement before source,
qualifier, reference or coordinate copies. Limits are 2 million source bases,
256 features, 4,096 location nodes, 256 declarations of each kind, 1 million total
CDS/reference-claim extracted bases, and 1 million reference bases. Repeated join
parts and repeated reference claims share the extraction budget. Translation is
performed once per eligible CDS; fusion/junction inspection reuses its cached
frame and amino-acid calls. Findings share a streaming budget before coordinate
copies. The final report is bounded to 2 MiB and 4,096 findings.

## Remaining integration and product work

- Worker schema 130, the Analyze menu panel and immutable-document MCP requests
  now route through this owner. Native/WASM compilation and focused scientific,
  transport and SDK checks pass; graphical acceptance remains deferred.
- Adapt the established annotation-review flow to explicit assumptions; a producer
  boolean or feature label must not become an authenticity receipt.
- Derive causal edit/assembly junction evidence through the existing fresh
  `assembly::junctions::inspect_plan` owner. A supplied boundary/history identifier
  cannot prove that an edit caused a frameshift. Compound-path fusion geometry
  currently remains indeterminate.
- Add declared element-chain, initiation-context, required-parts and unannotated-ORF
  checks using existing search/ORF owners and pinned reference data.
- Add evidence-bearing generic structural profiles after reviewing their sources.
  No vector-class thresholds, empirical predictions or specialized biological
  capabilities are included in this foundation.

Nine authored literal tests cover independent CDS defects, ambiguity, code changes,
assumption dependencies, translation disagreement, forward/reverse register,
origin-spanning joins, reference identity, restored-report tampering and repeated-location amplification.
All nine pass, alongside the connected transport and actual MCP SDK cases.
