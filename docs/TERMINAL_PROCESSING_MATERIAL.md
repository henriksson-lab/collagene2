# Source-backed terminal processing material

The **DNA ends and annealing** panel includes **Tail or treat source material
directly**. Its action selector begins blank: explicitly choose 5′ phosphorylation,
5′ dephosphorylation, or a single A/T/G/C nucleotide added to each strand’s 3′ end.
Preview shows the actual source and product geometry and retained limitations.
Explicit installation reuses the existing native material review and atomic
archive transaction. This operation creates no primers or primer-order attachment.

The core adds one recipe variant to the existing native material representation:

```rust
Derivation::TerminalProcessing {
    action: TerminalAction::TreatEnds { treatment },
    // or TerminalAction::Tail3Prime { base }
}
```

`TreatEnds` accepts only the existing typed phosphorylation/dephosphorylation
operations. There is no “none,” implicit action, declared replacement end state,
imported stock state or intended assembly target. The ordinary surrounding
`material_derivation::Request` retains exact source event, explicit original or
saved origin and output name; `save_primers` must be false.

`terminal_processing_preview_revision(workspace, source_event, origin, action)`
and `preview_terminal_processing(...)` require the owning native workspace. Before
source cloning or scientific work they attest exact event/revision and all native
saved Run/Derived dependencies. The shared source admission also serves reusable
OE previews. Equal DNA under another event is not the selected source. Installation
additionally requires the captured source head and whole-workspace transaction to
remain current. Preview uses the original bound event, not substituted live heads.

An original reference is an explicit intact DNA hypothesis with unknown terminal
chemistry/modifications. Processing its 5′ ends does not establish 3′ hydroxyls.
Consequently raw original-reference DNA can be phosphorylated/dephosphorylated,
but it cannot gain tailing eligibility through that treatment alone. Physical or
reference-only markers require their exact saved native origins, including when
only the `duplex_model` marker was removed from a detached preview. Terminal
preview provenance uses the existing restrictive `source:material_derivation:`
metadata namespace.

The single native execution path calls `TrackedDuplex::treat` or `.tail` on freshly
reconstructed private source state. Treatment changes the two external 5′ groups
on intact linear strands; the existing kernel rejects circular or nicked sources.
It does not seal a nick, infer repair or change strand letters/modifications.

Tailing requires an actual fully paired blunt linear duplex and known accessible
3′OH on both affected ends. Exactly one declared base is appended to each stored
5′→3′ strand. The top-reference coordinates stay fixed, and the bottom start becomes
−1. Original calls, 5′ groups, internal nicks and positional modification mappings
remain; new bases are unmodified under the existing standard-nucleotide model.
Tailing does not establish kit/enzyme activity, homogeneous tail length, yield or
isolation. It does not remove an internal nick or make a second tailing action on
already staggered material eligible. Repair is a separate explicit operation.

Preview context binds exact source event, native-origin hash, source sequence hash
and explicit action. `TerminalProcessingPreview::validate_binding` checks bounded
documentary relationships without replay or endpoint science: unchanged top letters
for treatment or exact one-base append for tailing, strand lengths/offset, end
groups, retained nicks, modification totals, generation declarations and sticky
incompleteness. Its summaries cannot restore private physical state or authorize a
product. Installation/restart freshly replay the native source recipe and compare
the recorded output snapshot/state/incompleteness.

A standalone PCR can therefore be tailed and phosphorylated, saved and reopened,
and chosen as Derived material for later TA/GC ligation. The later cloning plan
still declares an independent target and intended fragment orientations; native
ligation retains its alternative products for review/protocol output. Processing
never upgrades upstream incomplete evidence or claims that alternative products
cannot occur. Source features keep their unchanged top-reference coordinates;
appended terminal bases are outside those original features.

Native fixtures are in `crates/core/tests/materials/derivation/material_terminal_processing.rs`:
all four literal tails, 5′ treatment/dephosphorylation, archive replay, TA ligation
with retained alternatives and protocol instructions, original unknown-chemistry
and ownership refusals, native cut-created hydroxyls with retained Dam calls,
actual internal nicks, and conditional upstream PCR. Root integration supplies
worker/MCP/UI coverage and records the integrated validation results.
