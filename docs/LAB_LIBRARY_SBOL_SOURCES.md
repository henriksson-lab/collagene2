# Selected publication of SBOL documentary sources

`collagene-library-sharing` recognizes exactly the attachment kind
`sbol-design-exchange-v1`. Its body is the exact imported structural SBOL 3.1
JSON-LD file, with media type `application/ld+json`; it is not a native cloning
run, capsule, target or material receipt. The SBOL file can contain multiple
components, private activities and alternative sequence representations, so its
mere presence is not permission to publish the entire file.

## Exact-original publication contract

The adapter runs after the ordinary native selected-document/event/run closure
has been derived and validated. It reads only native molecule revisions retained
by that closure. The source projection's existing four documentary metadata keys
are interpreted together:

- `sbol:component_iri`
- `sbol:sequence_iri`
- `sbol:source_hash`
- `sbol:exchange_profile`

The profile must equal `collagene-sbol-3.1-structural-jsonld-1`. The source hash
must identify the unchanged raw SBOL attachment. The adapter parses that graph
through the bounded design-exchange crate and confirms both selected IRIs and
their declared Component→Sequence relationship. Incomplete claims, absent source
bytes, wrong attachment contracts and unknown profiles are explicit errors.
Descriptions, roles and arbitrary strings are never scanned for native event
IDs or presumed dependencies.

The adapter independently follows the supported structural reference fields:

- Component: sequences, features and constraints.
- Feature: locations; SubComponent additionally references its instance Component
  and source locations.
- Location: its referenced sequence.
- Constraint: subject and object features.
- Sequence: no further structural object references.

Referenced objects must have the expected supported class, exactly one class and
no retained unknown properties. Design-exchange's structural parser checks its
existing ownership/reference rules. No remote reference is fetched. Namespaces,
ontology/type/role IRIs and literals are not treated as instructions or as generic
foreign keys. Unsupported graph-bearing property contracts are refused, rather
than guessed, stripped or traversed through arbitrary JSON text.

For a Component with multiple sequence representations, every alternative must
be explicitly selected by a molecule revision in the retained native closure.
A second local document existing elsewhere in the workspace is insufficient.
A real retained cloning run using both declared source representations can make
both necessary; an unrelated local/private sequence remains excluded.

The original file is retained **only if its complete object set equals this
independently derived structural closure**. A sibling Component, Collection,
unrelated Activity or orphan graph object makes automatic publication fail.
The whole raw graph is never an automatically allowed opaque leaf. Recognized
SBOL attachments unreferenced by retained native revisions are omitted, without
blocking publication of an ordinary unrelated plasmid. Other unknown attachment
kinds still use the existing explicit unsupported-contract error.

## Native history and documentary meaning

The adapter adds only the proven necessary original attachment hash to the
existing permitted asset inventory. It does not rewrite molecule metadata, DNA,
features, native revisions, source events, plans, run records, nested capsules or
physical evidence. It does not replace source bytes under an old hash. The
original local workspace and its unpublished branches remain untouched.

A retained source hash records an original documentary source, not a claim that
an edited native molecule still equals the original Sequence. Edited DNA is not
matched against SBOL letters to select a source or decide which product is
physical. Ordinary native reconstruction remains responsible for all physical
material and target checks. The adapter never creates a native run or bootstraps
physical state from the SBOL graph.

Both browser/native selected export and an independently validating server use
the same closure rule. A supplied archive containing extra SBOL assets still
fails the host's actual-versus-permitted inventory comparison. No additional
whole-workspace sharing permission is inferred from a source attachment.

## Explicit reviewed selected-subset publication

The exact-original path remains strict. An explicit alternative now exists in
`library_sharing::sbol_subset`; see [SBOL documentary subsets](SBOL_DOCUMENTARY_SUBSETS.md).
It creates a new selected archive containing canonical structural SBOL bytes and
a separate omission receipt. It never rewrites a native molecule, revision,
event, capsule or original `sbol:source_hash`. The original file remains local.
Original-absent receiving checks establish exact scope and bytes, not proof of
derivation from the absent original or scientific correctness. Unsupported
selected properties/classes and forced unselected alternatives remain errors.

## Verification fixtures and bounds

The `sbol_publication` target contains six independent fixture groups. One
projects literal SBOL DNA through the real converter, records an actual native
PCR, publishes the saved product, reopens the archive, reconstructs a fresh
capsule and performs a downstream PCR. It checks unchanged original bytes,
source revision, native product sequence/status and omitted unrelated native/SBOL
records. Other fixtures cover private sibling components; alternative sequence
selection by actual retained PCR inputs; unknown property/class contracts and
missing/wrong pins; unrelated recognized SBOL files versus unknown kinds; and an
edited native molecule whose original source remains documentary.

The existing native workspace bounds and design-exchange bounds remain in force.
SBOL traversal consumes the same aggregate million-edge publication budget as
native event/run dependencies, uses a visited object set and performs no network,
filesystem or scientific subprocess operation. Tests are intended for the parent
integration batch; this source work did not run compiler or browser checks.
