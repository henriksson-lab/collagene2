# TOPO preparation and inspection

The Cloning → TOPO panel is a source-bound inspection workflow. It begins with
explicit preparation declarations, not a commercial kit name or an inferred
activation state.

For either preparation mode, open separate reference documents for the actual
prepared vector top and bottom strands (both written 5′→3′), the PCR template,
and the independent intended target. Set the target as the active document.
Choose each source explicitly; a shared display name does not bind a source.
The same immutable source may fill distinct declared roles, but the target
must be a separate document.

The TA form declares equal-length vector strands paired over all but one base,
with one unpaired 3′ nucleotide on each strand. It requires explicit confirmation
of each vector strand's 5′ phosphate and 3′ phosphotyrosyl chemistry and a written
preparation assumption. These declarations are hypotheses; a filename, vendor
or annotation supplies no chemical authority. Actual DNA and pair compatibility
are checked by the physical engine.

Enter complete forward and reverse primers in their ordered 5′→3′ direction,
including any 5′ tails, and declare each 3′ annealing-suffix length. PCR geometry
must reconstruct uniquely from the selected original template. Explicitly opt in
to adding one A to each freshly generated PCR strand's 3′ end. No terminal
addition, inherited modification or polymerase activity is inferred from a name.

Inspect both orientations before selecting an exact candidate-specific transfer
request. Then select a physical component and source-call seed for projection.
A missing unique paired contour withholds reference DNA; it does not choose a
path using the target. The full graph retains other components, unpaired calls,
nicks and transfer outcomes. Resource limits and incomplete evidence remain
visible.

Check the selected component against the separate target. Each allowance in the
physical acceptance policy is explicit: retained nicks, unknown terminal
chemistry, unapplied transfers, other components, circular rotation and reverse
complement comparison. An accepted result means agreement under that declared
graph policy; it does not establish reaction efficiency or experimental success.
Download source-bound JSON evidence or, after explicit review of an accepted fully specified graph, open a `CheckedTopo` candidate for ordered-plan verification and HTML/PDF protocol export. The plan freshly checks the same explicit policy; a match can be consistent under its declared graph assumptions. Source/event changes revoke pending results and
require fresh inspection. Imported reports cannot reinstall physical authority.

The panel now also offers **Directional CACC — actual blunt PCR**. This is a
second declared preparation profile over the same native general recipe and
physical inspection; it adds no kit-activation or flap-repair model.

For the directional profile, supply the complete actual vector top and bottom
strand reference sequences in their own 5′→3′ orientations. The bottom source is
four bases longer: its initial **GGTG** is an unpaired 5′ extension, and its paired
core starts at bottom offset four. `GTGG` read from a 3′→5′ drawing is not the stored
sequence. Both strand records, terminal declarations and distinct activation IDs
remain explicit; neither reverse strand nor activation is synthesized from a
plasmid annotation.

Directional insertion uses unique exact PCR with **no 3′ terminal addition**. If a
previous TA configuration selected A addition, explicitly uncheck it; switching
modes does not silently rewrite that choice. The actual incoming PCR top must
start with CACC. This may already be in the template or be supplied by a literal
forward-primer 5′ tail. The adjacent CACC-prefix action only proposes a changed
complete ordered forward oligo. Review it and explicitly Apply; the annealing
suffix length, reverse oligo, original template and independent target are
unchanged. Existing manual primer fields remain the editing surface. The proposal
is bound to the complete current source/event/epoch/primer context and disappears
when it changes. No CACC, A base or phosphorylation is silently introduced.

A synthetic mechanism example uses vector top `AAGGGGCATCCCTT` (14 bases),
vector bottom `GGTGAAGGGATGCCCCTT` (18 bases), and PCR template `CACCATGACG`.
The complete ordered primers `CACCA` and `CGTCA`, each with five annealing bases,
reconstruct the blunt ten-base insert. These short primers and declared activation
are a bounded mechanism fixture, not a recommended experiment or commercial stock.
The separately declared reference `AAGGGGCATCCCTTCACCATGACG` does **not** turn the
result into closed DNA: the native graph retains the four displaced calls as an
attached flap and withholds the paired contour even under permissive target
policy. All actual strands, nicks, source calls and reciprocal components remain
available through the existing review and source-bound JSON download.

## Shared native declaration builder

`core::topo_preparation::PreparationRequest` contains an explicit mode, prepared
vector declaration, literal exact-PCR declaration, engine limits and projection
options. The vector supplies its two strand interval lists, actual termini and
pairing span. `compile` validates bounded declarations and returns the ordinary
`topo::Request`, strand indices and directional flag. `source_preflight` and
`compile_preparation` additionally apply the existing native borrowed-original
model guards. These functions do not simulate activation or PCR and confer no
scientific authority. Both web and MCP adapters can use this same DTO; heavy
orientation work continues through the existing native worker tools.

The builder preserves an empty transfer list for discovery. Its returned
orientation-specific requests must be used for subsequent projection and target
checking. The web TA compatibility helper retains its old role names, material
IDs, pairing and request bytes; the directional profile has distinct declaration
IDs and no implicit terminal tail. Complete ordered primers in these preparation
profiles are bounded at 200 bases. More general source/graph recipes remain
available through the existing adapter rather than being relabeled as this profile.

## Remaining physical scope

This increment does not complete F15. Further transfer composition, source-call
accounted directional flap processing, generic nick sealing and a proven
fresh-graph-to-duplex bridge remain separate work. A target-policy allowance must
never delete an adduct, nick, flap, unpaired call or coupled component, or promote
a reference-only contour to digest/PCR/ligation material. No current action claims
host repair, transformation, template clearance or kit efficacy. A withheld graph
can still be inspected and downloaded; opening it as closed reference DNA is not
a workaround.

Browser/visual/performance acceptance remains deferred until all in-scope code is
complete. Source review, native fixtures and compiler results do not establish GUI
usability or experimental performance.

Try the supplied [synthetic TA example](../fixtures/topo-ta-demo.md) for literal source strands, primers and an independent target with two retained nicks.
