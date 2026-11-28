# collagene-enzyme-catalog

MIT Rust library containing a versioned restriction enzyme catalog.
This crate has no dependency on Collagene's core, browser, native process APIs,
or matching algorithms. It performs no network requests, digestion simulation,
site matching, methylation propagation, or reaction validation.

The independently curated `2026-11-28.2` edition contains **40 enzyme identities**
(including the original twelve starter geometries) and **eight explicitly scoped
NEB formulations**, with 64 documentary condition observations and eight
product-scoped documentary source-host presets. Additional
geometries do not imply that vendor conditions have been curated for them.
Product formulations are separate records referencing an enzyme identity;
EcoRI-HF is a formulation of the EcoRI recognition/cut geometry, not another
biochemical identity counted toward the forty.

```rust
use collagene_enzyme_catalog::{Catalog, Condition, Lookup, Resolution};
let catalog = Catalog::builtin()?;
let enzyme = match catalog.lookup_enzyme("ecori") {
    Lookup::Unique(enzyme) => enzyme,
    _ => return Err("missing or ambiguous enzyme".into()),
};
assert_eq!(enzyme.geometry.recognition, "GAATTC");
let condition = Condition::BufferActivity {
    buffer: "rCutSmart".into(), temperature_milli_c: 37_000,
};
let evidence = catalog.resolve("R3101", &condition, None)?;
assert!(matches!(evidence, Resolution::Known { .. }));
let pinned_content = catalog.content_hash();
# Ok::<(), Box<dyn std::error::Error>>(())
```

## API and scientific semantics

`Catalog` owns a private validated snapshot. `edition()` returns shared access;
callers may clone a proposal, but cannot mutate an installed catalog.
`content_hash()` is BLAKE3 of canonical compact JSON, including provenance.
Pin both the edition name and content hash when retaining a selection: an
edition string alone is not an authenticated uniqueness guarantee.

`Geometry` has an uppercase IUPAC recognition sequence and signed top/bottom
cut **boundaries**, using the forward recognition start as zero on both strands.
For example, BsaI `GGTCTC(1/5)` becomes `(7, 11)`; BbsI is `(8, 12)` and SapI
is `(8, 11)`. Bottom is not an index counted from the opposite strand's 5′ end.
Geometry relation comparison canonicalizes reverse-complement descriptors:
`(motif, top, bottom)` equals `(reverse_complement(motif), n-bottom, n-top)`.
This transforms descriptors only and does not search a DNA sequence.

`relation` first requires both supplied records to exactly equal records installed
in this sealed catalog, including citations and aliases, and returns `Result`.
Value-identical clones are accepted; foreign or altered records are rejected
before ID comparison or coordinate arithmetic. Records unchanged across editions
are value-identical; edition provenance itself must be pinned with the catalog
content hash. It then distinguishes identical catalog IDs, same-recognition/same-cut
isoschizomers, same-recognition/different-cut neoschizomers, and different
recognition sequences. HpaII/MspI are separate same-cut identities;
SmaI/XmaI and Acc65I/KpnI are neoschizomers. They are never aliases. Some vendor
charts use “isoschizomer” as an umbrella term including neoschizomers; this API
uses the more specific mutually exclusive results. Same-cut geometry does not
imply equal methylation sensitivity, conditions, or guaranteed ligation.

Lookup is ASCII case-insensitive with surrounding whitespace ignored. Canonical
IDs and names cannot collide. Aliases need citations. Unacknowledged alias
collisions are rejected; shared aliases require every owner to opt in and
produce `Lookup::Ambiguous`, sorted by ID. No arbitrary first-match selection.
The builtin enzyme alias lists are intentionally empty; formulation product
codes are documented aliases. `search` returns ID-ordered substring matches.

`Formulation` separates product code, vendor, and explicit lot/formulation scope
from geometry. Evidence keys and values are typed: incubation temperature,
buffer activity at a stated temperature, and Dam/Dcm/CpG sensitivity categories.
A key also has a named documentary assay context and source citations.
`resolve` returns `Known`, `Unknown`, `Conflict`, or `ContextRequired`:

- No evidence is unknown. It is not compatible, incompatible, or zero activity.
- Missing formulations return an error; facts are never inherited from another
  product, enzyme family, or older formulation.
- Different contexts require an explicit choice instead of treating historical
  or assay differences as an equivalent measurement.
- Conflicting known values within one context remain inspectable as conflict.
- Unknown observations and their reasons remain in the returned evidence even
  when another observation supplies a known value.

Activity percentages are approximate vendor unit-assay measurements, not
predicted digestion efficiency. The EcoRI r2.1/rCutSmart star-activity qualifier
is retained alongside percentages. The current SmaI chart reports **37 °C**. Its qualifier also retains NEB’s
recommendation to use 25 °C for digestions longer than one hour; a chart
temperature does not replace duration and context limitations. Methylation
categories are documentary evidence, including overlapping-combination
uncertainty. There is no inference about a particular molecule, motif context,
hemimethylation, incubation time, mixed buffers, or reaction outcome.

## Documentary source hosts

`Catalog::hosts()` exposes eight supplier/product-scoped E. coli presets.
`Catalog::host(id)` resolves an exact current ID; family names and supplier
variants are not automatic aliases. `Catalog::validate_host(record)` requires
full current value membership after bounds checks, including scope and citations.
The cohesive `hosts` module owns `HostPreset` and
`MethylaseExpectation::{Active, Deficient, Unknown}`. These expectations concern
host methylase documentation, not measured site occupancy in the submitted DNA.
All current CpG expectations are unknown; TOP10 Dcm remains unknown in the reviewed
source set. Host selection cannot overwrite existing positional chemistry, imply
propagation or erase incoming modifications. An explicit molecular generation
assumption is a separate scientific input owned by the application core.

The current schema requires `Edition.hosts`, with at most 128 records and bounded
nested qualifications/citations. Host assertions participate in catalogue/source
integrity hashes and current JSON/Postcard encoding. No migration or fallback
accepts superseded internal schemas.

## Provenance and independent curation

Every assertion has source IDs and specific row/section locators. Each source
records publisher, exact URL, review date, source edition, and a BLAKE3 digest
of the curated assertions citing it. These hashes bind this independently
curated subset; they **do not authenticate the upstream vendor document**.

For the buffer chart and EcoRI-HF product page, source records additionally
contain SHA-256 and byte length of the actual retrieved HTML snapshot. Other
pages were reviewed through rendered official documentation while raw retrieval
returned HTTP 403; their upstream snapshot hash is explicitly absent, with a
reason. No hash is fabricated and no full HTML/PDF is redistributed. An import
checks curated-fact hashes and the schema; it cannot independently establish
that a source assertion was accurate or that the URL still has those contents.

Primary documentation reviewed 2026-11-24:

- [NEB recognition chart](https://www.neb.com/tools-and-resources/selection-charts/alphabetized-list-of-recognition-specificities)
- [NEB isoschizomer chart](https://www.neb.com/en/tools-and-resources/selection-charts/isoschizomers)
- [NEB 2019 catalog](https://www.neb-online.de/wp-content/uploads/2019/05/NEB_Catalog_RGB_150dpi.pdf), for explicitly historical base-enzyme geometries whose current chart emphasizes engineered products
- [NEB buffer performance chart](https://www.neb.com/en/tools-and-resources/usage-guidelines/nebuffer-performance-chart-with-restriction-enzymes)
- [NEB EcoRI R0101](https://www.neb.com/en/products/r0101-ecori)
- [NEB EcoRI-HF R3101](https://www.neb.com/en-us/products/r3101-ecori-hf)
- [NEB XmaI R0180](https://www.neb.com/en/products/r0180-xmai)
- [NEB SmaI temperature-change FAQ](https://www.neb.com/en/faqs/which-incubation-temperature-should-i-use)
- [NEB Dam/Dcm/CpG guide](https://www.neb.com/en/tools-and-resources/selection-charts/dam-dcm-and-cpg-methylation)

Code and independent fact selection are MIT. Vendor names and trademarks identify
products and imply no endorsement. No cloning-app code or REBASE bulk database
was reused. Future additions require source review and a new pinned edition.
The initial geometry schema intentionally excludes nickases, multiple paired
cut sites, and methylation-dependent recognition such as DpnI rather than
misrepresenting them as ordinary single-pair-cut IUPAC enzymes.

## Import, limits, and development

`Catalog::from_json` refuses input beyond 2 MiB before deserialization. Serde's
recursion limit and unknown-field rejection remain enabled. Schema version must
be 3. Caps are 128 sources, 1,000 enzyme identities, 2,000 formulations, 16 aliases
per record, 64 observations per formulation, 16 citations per assertion, and
bounded strings/recognition lengths. Cuts must be within -1,024..=1,024. Invalid
IUPAC, duplicate identities, unresolved references, mismatched condition/value
kinds, invalid unknown-condition keys, malformed hashes, and unsupported schema
versions fail closed. Dates use the proleptic Gregorian calendar, years
0001–9999, including actual month lengths and century leap-year rules. Native
DTO proposals also pass a capped, nonallocating serialized-size preflight before
canonicalization, hashing or edition cloning. Export uses the same validated
bounded compact JSON.

`seal(Edition)` is an authoring API that validates a proposal and computes its
curated-fact checksums. It must not be substituted for validating received JSON:
sealing deliberately creates new checksums. The offline authoring utility is
`cargo run --example seal -- proposal.json`. It does not fetch or curate facts.
Canonical order is entity ID order, alias-name order, and serialized evidence
order; documentary qualifier order is retained. Import has no mutable global
registry or filesystem requirement. Only the optional authoring example uses
filesystem I/O, reading at most the byte limit plus one byte before rejecting an oversized proposal.

The crate is a member of the Collagene workspace. Catalog-backed application
selection exposes all 40 geometries with separate formulation evidence. The core
reaction-condition owner compares exact selected formulations and assay contexts;
it does not turn documentary activity into complete-digest success. Current schema
3 uses externally tagged JSON enums and the same Rust variants in Postcard. The
superseded schemas 1 and 2 are rejected. Optional `json-schema` exposes current typed DTOs.
Verification:

```sh
cargo test -p collagene-enzyme-catalog --locked
cargo check -p collagene-enzyme-catalog --locked --target wasm32-unknown-unknown
```

Tests use literal independently written geometry/condition expectations and
synthetic invalid/conflicting/ambiguous proposals, not another cloning program
as an oracle. Windows/macOS/Linux portability follows from Rust + serde/blake3
and no platform APIs in the library; those native platforms still require CI
execution. WASM compilation is checked locally, without an external scientific
runtime.

The 2026-11-28.2 admission adds BamHI, BamHI-HF, BsaI-HFv2 and BsmBI-v2
formulation facts. The BsmBI-v2 four-buffer matrix belongs to its 55 C unit-assay
context; a separate 37 C activity observation is not extrapolated into an
uncurated buffer/temperature pair. Formulation qualifications retain narrowly
specified star-assay observations without predicting star risk or declaring
current-lot safe-dose thresholds. See [data notices](data/NOTICE.md) and
[the current admission audit](../../docs/audits/restriction-next-data-current.md).
