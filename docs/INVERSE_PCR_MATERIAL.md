# Inverse-PCR templates from saved native material

The existing inverse-PCR compiler now accepts either an original circular DNA reference or an explicitly selected generic native material origin. There is one compiler and one scientific operation graph, not a separate saved-template method.

```
mutagenesis::inverse_pcr_revision(source, target, origin, options)
mutagenesis::compile_inverse_pcr(source, target, origin, options)
mutagenesis::design_inverse_pcr(source, target, origin, options)
```

`origin` is `Option<&run_material::MaterialOrigin>` before `options`. `None` retains the original-reference evidence guard. `Some(Run)` emits `ReconstructRunProduct`; `Some(Derived)` emits `ReconstructMaterial`. The next three operations are actual exact PCR, explicit terminal phosphorylation and the current intended circular ligation with bounded alternatives. Each operation consumes the preceding actual private tracked state. Compilation does not replay a source or claim molecular success; fresh ordinary verification performs that work.

The independently supplied circular target never chooses PCR geometry. Oligos come from the declared source interval and literal insertion, with the explicit insertion-tail split. Newly synthesized PCR strands are unmodified under the native ideal-PCR model. This does not remove template methylation in the saved source, implement a proprietary kit, destroy the original template, predict yields or establish isolation.

Both preparation and design results include `context_revision` and `validate_binding(expected)`. The context binds current preparation/scientific engines, exact source revision, full native origin, independent target revision and all declared options. Returned options retain the exact request, including insertion case/whitespace; oligo DNA is normalized separately. Thus readonly validation never mistakes normalized output for a different request. It checks declared primer/step correspondence and report pins without PCR search or native replay. Fresh verification and exact workspace attestation still determine acceptance.

Before hashing or cloning origin/source models, combined source/target models and borrowed request serialization are bounded. A 64KiB envelope reserve leaves source+target+origin+options below the ordinary 32MiB request ceiling; complete results also remain below32MiB. The existing shared native ancestry ledger bounds all nested Run/Derived capsules. No metadata is stripped to evade the source guard.

An origin's revision alone does not prove which historical event supplied an equal snapshot. Hosts must attest the selected origin against the exact captured document/event/revision in the retained native workspace. A matching molecule in another document is insufficient. Web and MCP adapters must freeze the origin and independent target along with the source event, and reject changed source/session/epoch/context on acceptance. The returned worker assessment remains bound to its original plan if a UI rebinds an input role; that rebound candidate requires fresh verification.

Intact, paired, covalent circular templates use the current native PCR model. Nicked circles or unsupported circular strand geometry remain diagnostic; there is no automatic repair or linearization. Origin incompleteness and original assumptions remain sticky even when all new operations produce the intended sequence. A complete new path is not an assertion of mixture uniqueness or experimental confirmation.

Native fixtures cover actual two-fragment circular ligation Run material, cut-free EcoRI digestion retained as Derived circular material, archive restart, literal insertion/deletion/replacement, inherited annotation and modification states, original-versus-derived source prefixes, a new saved run with HTML/PDF/oligo outputs, exact-event alias refusal, wrong targets, tampered receipt/steps/options, nicked-circle rejection and capped-origin incompleteness. Original source heads/history remain unchanged. No new target-free derivation kind or scientific operation is required.
