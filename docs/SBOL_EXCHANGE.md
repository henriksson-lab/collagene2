# SBOL 3.1 structural JSON-LD exchange

Native MCP exposes `import_sbol`, `inspect_sbol`, `project_sbol`,
`export_sbol_graph`, and `export_sbol`. Import retains an immutable original
graph; projection requires explicit component/sequence selection and omission
review. Sources are either imported evidence handles or an exact workspace's
original-graph attachment. Original/expanded graph export creates no cloning
authority. Native reference export uses the same physical-evidence refusal guard
as the web dialog. The current integration batch has dedicated SDK cases;
validation is recorded in the implementation status.

`collagene-design-exchange` is an independently implemented MIT Rust crate using
only the shared `collagene-model`, serde/JSON, hashing and error types. It is
portable to WASM and native targets. There is no external executable, RDF server,
JSON-LD download, subprocess or network dependency.

The selected standard is **SBOL 3.1.0**, with implementation profile
`collagene-sbol-3.1-structural-jsonld-1`. This is deliberately a documented
structural subset, **not a claim of full SBOL structural-subset compliance or a
complete SBOL validator**. No third-party implementation code was used.

## Primary references and examples

The official [SBOL data-model release page](https://sbolstandard.org/datamodel-specification/)
identifies 3.1.0 as the current SBOL 3 release. The normative
[3.1.0 specification](https://sbolstandard.org/docs/SBOL3.1.0.pdf) supplies object
semantics: §§6.1–6.4 cover identities, sequences, components, hierarchy, locations
and constraints; §6.9 covers collections; §§8–9 distinguish RDF serialization and
compliance scope. Tables1–3 specify supported encoding/category/topology terms;
tables5–6 specify orientation terms. Its pLac serialization example was reviewed
as an elementary relationship fixture, not copied as application code.

The official [SBOLTestSuite SBOL3 collection](https://github.com/SynBioDex/SBOLTestSuite/tree/master/SBOL3)
and [toggle-switch N-Triples fixture](https://raw.githubusercontent.com/SynBioDex/SBOLTestSuite/master/SBOL3/toggle_switch/toggle_switch.nt)
were inspected to establish that real designs carry components, child instances,
constraints, component references and interactions rather than just sequences.
The latter two object families are retained as unsupported objects in this
profile. These external examples are **research inputs, not an executed
interoperability/conformance corpus**. Repository fixtures are newly authored
literal designs with independently specified expected relationships and native
coordinates. Research baseline:26 November2026.

## Accepted JSON-LD profile

Input is a flat expanded JSON-LD node array, or a wrapper containing `@graph`
and optionally `@context`. Every node has an absolute stable `@id` and one or
more `@type` IRIs. Relative identities, blank nodes, duplicate JSON properties,
duplicate object identities and aliased duplicate predicates are refused.
Objects with identical IDs are not silently merged.

The only compact context accepted is an inline object of exact, pinned prefix
mappings:

| Prefix | Exact expansion |
| --- | --- |
| `sbol` | `http://sbols.org/v3#` |
| `prov` | `http://www.w3.org/ns/prov#` |
| `xsd` | `http://www.w3.org/2001/XMLSchema#` |
| `rdf` | `http://www.w3.org/1999/02/22-rdf-syntax-ns#` |

Any subset is permitted. Prefix replacement, remote/array/scoped contexts,
`@vocab`, keyword remapping, nested object definitions, blank nodes, RDF lists,
`@set`, `@reverse`, JSON literals and implicit reference coercions are refused.
An IRI-valued property must use `{"@id":"..."}`; a plain string remains a literal.
Scalar JSON literals and expanded scalar `@value` objects are supported. Supported
text fields require unlocalized string literals; integer fields accept JSON
integers or explicit `xsd:integer` lexical strings. No context URL is fetched.

Unsupported **typed graph objects and predicates** remain present in the immutable
retained graph and original bytes; their presence produces bounded diagnostics.
Unsupported JSON-LD syntax is a refusal, not an invitation to guess its RDF
meaning. This distinction prevents an unknown context from silently changing
which sequence, coordinate, authority or hierarchy the program believes it read.

## Typed design graph

The typed graph implements:

- Component types, roles, sequence links, owned features and owned constraints.
- Sequence elements/encoding, including sequence-less documentary designs.
- SequenceFeature and SubComponent, with instanceOf, sourceLocation,
  roleIntegration, roles, orientation and explicit locations.
- Range, Cut and EntireSequence, each linked to its sequence, with optional
  orientation and signed integer order.
- Constraint subject/object/restriction and Collection membership.
- Common identity/name/displayId/description and required TopLevel namespaces.

Namespaces are explicit HTTP(S) URLs. URL-valued TopLevel identities must have
the declared namespace prefix. Other stable absolute IRIs, including URNs, are
retained without remapping. SubComponent references are checked for cycles.
Supported feature/constraint/location ownership is unique, constraint endpoints
belong to the owning component, and normal/source locations refer respectively
to the parent/child component's declared sequences. Present ranges and cuts must
fit the referenced explicit elements. Unresolved references remain typed with
visible diagnostics; nothing is fetched or assembled to fill them.

Constraints are **declarations**, not proof of biological validity or a native
cloning plan. The graph does not evaluate every positional/functional constraint
or reconstruct interactions, combinatorial derivations, component references,
implementations or computational models. Those objects/properties remain in
exchange data; their semantics are not silently flattened.

## Selected component → native sequence reference

`Document::project(&Selection, ProjectionPolicy)` requires an explicit component
IRI and one of its linked sequence IRIs. The conversion uses that sequence's
actual elements; it does not assemble child sequences or synthesize sequence
from hierarchy/constraints. Native DNA/RNA/protein is selected from the explicit
Component SBO category, with EDAM1207/1208 encoding correspondence checked. It
never guesses DNA versus RNA from A/C/G-only letters.

Missing DNA topology requires an explicit caller choice. A caller choice cannot
contradict declared topology. RNA/protein use the applicable linear default;
unsupported circular proteins fail native validation. Native sequence case is
normalized with a visible warning; whitespace/non-ASCII IUPAC elements are
refused so coordinates cannot shift silently. Unsupported extra physical types,
including single-stranded DNA, are refused **even with AllowOmissions**: their
restrictions cannot safely disappear into an ordinary duplex-capable DNA model.

Ranges convert from1-based inclusive SBOL coordinates to0-based half-open native
coordinates. Each location needs explicit supported orientation, either locally
or at its feature; conflicting orientations are refused/omitted explicitly.
Multipart locations need unique signed `order` values and retain that biological
join order. Circular origin crossing is represented by multiple explicit ranges,
not a reversed start/end or inferred sorting. Reverse-complement protein
locations, Cut→base-range conversion, missing coordinates and ambiguous order
remain explicit refusals or omissions. EntireSequence maps only when its linked
sequence elements exist.

Located SubComponents can appear as native reference annotations with an exact
instance-IRI qualifier, but child hierarchy/source locations/role integration
are explicitly reported as unrepresented semantics; child DNA is not flattened.
Constraints, additional sequence representations, foreign predicates/classes and
other unsupported selected design semantics also produce named conversion
losses. Imported features are marked for review.

`Strict` refuses any representation loss. `AllowOmissions` returns the native
reference plus the full visible loss list; it never bypasses missing sequence,
conflicting alphabet/topology or unsafe physical-type checks. The original graph
remains available independently. Native metadata contains only documentary SBOL
identities, source hash, profile and selected descriptions—not foreign duplex
models, verified status, run receipts or physical authority.

## Export and public API

- `import(bytes) -> Document` validates this bounded profile.
- `Document::{objects,issues,revision}` provide typed graph/diagnostics/evidence
  identity, without mutable raw fields.
- `Document::original_bytes()` returns exact original evidence.
- `Document::expanded_bytes()` emits the entire retained normalized expanded
  graph, including unsupported predicates/objects. It is deterministic for the
  retained graph, **not RDF dataset canonicalization**.
- `Document::project(selection, policy) -> Projection` returns a native molecule,
  explicit selected IDs, warnings and losses.
- `export_molecule(molecule, ExportIds, policy) -> Export` emits a real Component,
  Sequence, owned SequenceFeatures and per-part Ranges. IDs/namespace are explicit
  caller choices; generated child IDs are deterministic and identity collisions
  are rejected. Reverse multipart biological extraction order is preserved.

A new native-molecule export is intentionally distinct from retained-graph
export. Native metadata, original files, source/chemistry evidence, fuzzy or
unordered/unresolved locations and unsupported feature fields cannot silently
become SBOL design semantics. `Strict` refuses these losses; `AllowOmissions`
reports them. Exported bytes are re-imported under the profile before return.
Full graph exchange preserves design relationships; a native sequence-reference
projection cannot promise full design-graph fidelity.

No graph object or imported annotation restores Collagene material, current
verification, host repair, read validation or native historical event authority.
URLs are identifiers/data, never executable operations or automatically opened
resources. Browser/worker/MCP integration and native attachment validation are
separate, coordinated follow-up work; this crate alone registers no file handler.

## Bounds and fixture coverage

Inputs and encoded outputs are capped at16MiB, graph objects at10,000, JSON
values at200,000, nesting at32, absolute IRIs at2,048bytes, sequence elements
aggregate at8MiB and diagnostics at1,000. Duplicate keys and value/depth bounds are
checked during JSON decoding. Hierarchy cycle checking is iterative. Evidence
output uses a bounded writer, not an unbounded serializer followed by truncation.

Nine fixture groups cover real child/sequence/constraint relationships,
source-location ownership, circular multipart order, exact retained graph
exchange, DNA/RNA/protein and encoding separation, unsafe type refusal, context
and duplicate-key attacks, cyclic/cross-parent graphs, missing orientation/order,
stable-ID native export and reverse joins, explicit export losses, input/depth
bounds, missing external sequences, and Cut/EntireSequence retention. No compiler,
browser, external validator or network runtime checks were run by the author;
these fixtures await the parent's coordinated native/WASM verification batch.

## Web exchange review and native export boundary

The File-menu SBOL dialog uses the same Rust profile with bounded, filterable
Component/Sequence choices. Both selections start empty; loading a graph never
chooses or concatenates a construct automatically. Topology and loss policy are
explicit. A separate review shows the projected alphabet, sequence size, feature
count, normalization warnings and every paged loss before the user opens a
native reference. The hierarchy/constraint summary is paged; identifiers are
rendered as data, never automatically visited.

The App callback receives `(Molecule, Option<WorkspaceAttachment>)` and commits
the new reference and exact original input together. The attachment contract is
`kind = sbol-design-exchange-v1`, `media_type = application/ld+json`, with the
original JSON-LD bytes directly (no envelope). Its BLAKE3 key equals the molecule's
`sbol:source_hash`. That original may contain unrelated components or unsupported
objects. **Local retention is not permission to publish the entire source**:
selected lab-library publication needs a separately validated graph-closure
contract, or an explicit refusal, before those bytes can cross that boundary.
No unknown attachment may be treated as a blanket publishable opaque leaf.

Current-molecule export uses editable explicit namespace/component/sequence IDs.
Placeholder defaults are deterministic from the captured native document identity
and exact revision. New export is prepared and reviewed before save. Original
source, retained expanded graph and selected native exports use the common
native-aware export helper: admission/epoch/source/selection guards are rechecked
before the write, and saved status waits for the real native receipt. Browser
export remains the common browser dispatch; no native failure falls back to a
Blob download. Disposed/stale dialogs cannot report completion.

The lightweight exchange crate cannot itself inspect private Core graph state.
Its `AllowOmissions` model-level export therefore **is not an application-level
permission to erase source restrictions**. Every Core-facing adapter must gate
known restrictive evidence before calling it. The web wrapper refuses documents
recognized by `chemistry::has_unmapped_physical_evidence`, restrictive physical/
reference-only metadata and exact native producer comments. It refuses even with
loss allowance, offering native document/workspace export instead. This profile
has no safe producer marker that would preserve such a restriction on reimport.
Future MCP/file adapters must apply the same boundary; ordinary reference DNA,
RNA and protein files remain supported.

## Explicit SO identifier aliases

The 3.1.0 specification's Table3 writes topology/strandedness IRIs with HTTP.
The official [pySBOL3 constants](https://raw.githubusercontent.com/SynBioDex/pySBOL3/main/sbol3/constants.py)
define an HTTPS SO namespace and consequently publish HTTPS variants for
`SO_LINEAR`, `SO_CIRCULAR`, `SO_SINGLE_STRANDED` and `SO_DOUBLE_STRANDED`.
Only that format vocabulary was consulted; no parser or validation code was
borrowed. The importer now recognizes precisely the HTTP/HTTPS pairs for
SO0000987, SO0000988, SO0000984 and SO0000985 as semantic aliases.

This lookup does not rewrite identities, change original bytes, normalize the
retained graph, follow an HTTP redirect or equate arbitrary HTTP/HTTPS URLs.
Retained graph exports keep the original term IRIs. New molecule exports retain
the specification's canonical HTTP topology spelling. Mixed-alias linear versus
circular declarations, and mixed-alias single versus double strandedness, are
rejected. Single-stranded DNA remains a hard native-conversion refusal under
both spellings; loss allowance cannot remove it.

Feature/location orientation is compared after decoding supported semantic
orientation, so the legacy `sbol:inline`/`sbol:reverseComplement` terms do not
conflict with their equivalent SO0001030/SO0001031 terms. A genuinely opposite
orientation still refuses strict projection. The ninth independent literal
fixture combines HTTPS topology/strand terms with legacy and SO orientations,
checks HTTP/HTTPS DNA/RNA behavior, retains source IRIs exactly, rejects mixed
conflicts/ssDNA, and confirms canonical new export. This adds no remote context
or general RDF/JSON-LD processing capability.
