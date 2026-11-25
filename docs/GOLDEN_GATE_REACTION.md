# Explicit Golden Gate reaction substrates

`golden_gate_reaction` compiles one declared material per role. `WholeSource`
uses the entire original DNA document, including a circular backbone. `Prepared`
uses the existing preparation enum and optional terminal treatment. In particular,
`Prepared { ExactPcr, treatment: None }` supplies the complete amplicon before
Type IIS digestion; `CatalogDigest`/`PcrCatalogDigest` deliberately supply the
selected isolated fragment. Mixed reactions can state exactly which material
enters the tube, without silently discarding fragments from whole inputs.

The shared preparation compiler emits original-source operations. The final
existing GoldenGate operation digests all declared substrates and enumerates the
bounded one-copy pool, including every digest piece and uncut circular carryover.
The actual native enumeration is captured in the same verifier pass, even when
probe index0 is unavailable or recuttable. No second search, detached physical
summary, target-derived substrate or implicit circular linearization is used.
The target is independently supplied and never updated from a product.

`inspect` returns original request/target, material-to-source mapping, a probe
plan, native diagnostics, all pool records and per-candidate target/role assessment.
The index0 probe exists to execute preparation and collect the actual enumeration;
it is not nominated, saved or installed as a design. Search absence under caps is
not complete absence. Uncut circles remain carryover rather than newly assembled
candidates. Recuttables remain inspectable but cannot be selected as site-free.

`select` requires exact context, index and directed/physical keys. It preserves
originals and every explicit preparation/treatment step and changes only the final
candidate index. The result is an unverified native plan. Fresh verification must
reconstruct it, including the exact reaction pool, against the independent target.
Existing actual-product/junction review, output creation and HTML/PDF steps use
that same native plan. Selection never restores imported physical authority.

The GUI shares its pool view across prepared-fragment and reaction modes. Users
explicitly choose whole source, current complete PCR amplicon or current prepared
fragment for each role. Whole-amplicon choice preserves primer strings but omits
operations that were specified after digestion; the choice is shown before
inspection and never edits the original primer/preparation draft. Source-event,
settings, target, epoch and displayed-report guards reject stale actions. Reaction
inspection is available before ordered-method computation, so an isolated-fragment
recipe is not a prerequisite for whole-source analysis. Per-view scope/options are
currently reset on source/configuration changes; broader unfinished-form persistence
remains required. Selected native plans preserve all submitted operations/options.

The model still assumes complete cleavage and one copy per declared role under
explicit fragment/orientation/state/product limits. It does not predict thermal
cycling kinetics, partial cleavage, repeated-copy stoichiometry, methylation
clearance, unlisted contaminants, experimental yield or exhaustive real-mixture
outcomes. Original modification/chemical limitations remain in native findings.
Existing source/model/serialization budgets apply before material copies.

Worker49 adds inspect115/result107 and select116/result108. Old prepared-pool
request/report layouts remain unchanged. Candidate11, plan9, science14 and native
archive1 are unchanged; no new scientific operation or external runtime is added.
