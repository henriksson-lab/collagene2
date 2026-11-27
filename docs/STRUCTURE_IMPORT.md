# Imported structure coordinates

`collagene-structures` is an independent MIT Rust crate for bounded PDB/mmCIF
coordinate import and source-indexed scene preparation, available on native and
WASM targets. This core batch does not yet expose a graphical viewer, worker job,
workspace attachment control or MCP tool. Those integrations remain part of F34.
RNA modification and antibody-specific analysis remain separate planned work.

## Source and meaning

`import(bytes, Format)` retains the exact input bytes, BLAKE3 hash, format and
current parser identity. Parsed atoms keep original f64 Angstrom coordinates,
record kind, model, author and label identities, entity, alternate location,
insertion code, deposited element, occupancy, B factor and formal charge when
present. Missing author identifiers are never replaced with label identifiers.
An absent element remains unknown; an atom named `CA` does not establish calcium.
PDB has author identities and no invented mmCIF label namespace.

`Structure::validate()` reparses the retained bytes and compares the entire
projection. Unsupported internal parser identities and edited projections refuse.
This binds an import to its bytes; it does not prove experimental validity,
sequence correspondence, chemical completeness or suitability for cloning.
Imported coordinate files can themselves contain predicted models. No provenance
claim about experimental origin is inferred from a file extension or syntax.

Original bytes remain the lossless representation. Unprojected PDB record types
and CIF categories are listed in import notes. Unknown versus not-applicable
connection values and quoted literal identifiers remain distinct. Other unknown
fields, missing-value spellings, cell/symmetry data, metadata and annotations stay
in the original source, without being promoted to supported computation.

## Profiles and bounds

| Profile | Implemented | Explicit limits |
| --- | --- | --- |
| PDB | Fixed-column ATOM/HETATM, MODEL/ENDMDL and CONECT; implicit single model | Decimal serial/model IDs; no hybrid-36; no mixture of implicit and explicit models |
| mmCIF coordinate data | CIF1 tokenization, scalar/loop `_atom_site`, deposited `_struct_conn` rows | One data block, one table per projected category; no CIF2, dictionary/save frames, mixed-category loops or uncertainty-bearing numeric fields |
| Retained connection declarations | PDB serial references and complete bounded `_struct_conn` row values | No distance/name-based bond inference, hydrogen addition, bond-order inference or symmetry expansion |

Sources are bounded to 16 MiB, 100,000 atoms, 400,000 connection declarations and
1,024 populated models. Identifier fields are at most 256 bytes. CIF processing
also bounds tokens, columns, tags, category counts and text-field size; PDB records
are at most 1,024 bytes. Non-ASCII input and nonfinite numeric coordinates refuse.
The supported CIF profile does not claim complete dictionary validation. Unknown
deposited quantities are not replaced by estimated occupancy or displacement.

## Scene contract

`scene::build` validates the import and selects one model, all chains or an exact
author/label chain, all alternates or an explicit alternate label, and whether to
show declared H/D atoms. Named-alternate filtering includes unlabelled atoms; a
global label filter is a display choice, not evidence for a complete conformer.
Models are never overlaid implicitly. Points preserve original atom indices.

The camera normalization first rescales f64 coordinates, then centers and fits
them before producing f32 GPU positions. Original coordinates stay unchanged;
even extreme finite input cannot overflow the fit calculation. Element colors
are display conventions. Missing elements receive an unknown color.

Edges come only from deposited declarations. Repeated PDB declarations retain
their source multiplicity but draw one edge, without interpreting bond order.
mmCIF edges require a supported covalent/disulfide relation, explicit identity
symmetry at both ends and unique endpoints satisfying both provided identity
namespaces and alternate/insertion/model qualifiers. Unknown qualifiers,
conflicting identities, ambiguous endpoints and incompatible alternate labels
withhold edges. Other relations/symmetry operations are counted separately.
Filtering can also leave declarations unresolved; their counts remain visible
in the scene result. Endpoint resolution has an aggregate work bound and uses
chain/atom indexes, rather than an unbounded all-pairs scan.

No folds, ribbons, secondary structure, symmetry mates, biological assemblies,
distance-inferred bonds or structure-to-sequence correspondence are generated.
Such additions require explicit models and independent checks.

## Specifications and verification

The coordinate profiles follow the published
[wwPDB PDB coordinate specification](https://www.wwpdb.org/documentation/file-format-content/format33/sect9.html)
and [PDBx coordinate guide](https://mmcif.wwpdb.org/docs/user-guide/resources/coordinate_section.html).
The streaming lexer implements the selected data-file syntax from the
[IUCr CIF1 specification](https://www.iucr.org/what-we-do/digital-standards/cif/cif1/file-syntax),
including quoted values and line-delimited text fields.

Authored fixtures cover fixed columns, separate models/namespaces, nonpolymer
identifiers, alternates, original bytes, connection ambiguity, missing versus
literal identifiers, malformed syntax, nonfinite values, bounded input, exact
camera geometry and tampered projections. Native fixtures and a WASM compiler
check validate this core; they do not establish GUI or vendor-corpus acceptance.
