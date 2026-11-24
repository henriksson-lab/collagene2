# Catalog-backed desktop virtual digest

The existing Digest inspector uses the embedded, versioned enzyme catalog through
worker schema 21 jobs 70–72. The original twelve enzyme controls, labels, ordering,
and default EcoRI selection remain. The same scrollable picker now includes all
forty distinct identities and searches their names and IUPAC recognition sequences.
No extra top-level tool or project hierarchy was introduced. Mobile review is deferred.

Each explicit run pins the catalog edition and canonical content hash, ordered
canonical enzyme IDs, and their actual recognition/strand-cut geometry. The order
is the original twelve followed by the remaining identities sorted by name; the
checkboxes retain independent identities even when their recognition and cuts
match. Formulation product codes do not implicitly select enzymes. Aliases are not
expanded by the digest bridge.

The collapsed **Catalog and formulation evidence** disclosure shows that pin and
selected geometry citations. Vendor formulations are nested separately. Their
incubation temperatures, buffer activities, and methylation categories retain the
exact documentary context, scope, qualifiers, review dates, and official source
links. Unknown data stays unknown, multiple contexts remain separate, and conflicting
observations have no selected value. Documentary activity percentages do not alter
predicted cuts or assert reaction compatibility. In particular, the SmaI extended
incubation exception (25 C for longer than one hour) accompanies the chart's 37 C
observation. EcoRI and EcoRI-HF buffer measurements are separate formulation facts.

The existing full-recognition-span navigation, source/event identity, global epoch,
ordered options, request generation, worker envelope checks, cancellation,
comparison/site/fragment/report limits, and bounded result paging remain. The main
thread performs bounded catalog resolution and lightweight result binding; matching,
physical fragment reconstruction, annotation projection, and fresh installation
remain in the worker.

Inspection JSON retains its underlying `report` field for existing consumers and
adds `catalog_pin`, `catalog_engine_version`, and `catalog_batch_revision`. Lazy
GenBank fragment export starts a fresh catalog product request and validates its
catalog/source/options/index binding before exporting. Opening a fragment uses the
existing single `DigestOpenPanel` installation pipeline, including a fresh worker
request, decoded archive/options/index binding, source-history parent, workspace
fingerprint/epoch checks, and atomic current/previous IndexedDB transaction.

Saved native and imported GenBank/EMBL catalog selections are explicitly unverified
producer claims. They can be downloaded as bounded JSON from the provenance
disclosure. They do not select current enzymes, restore a current result, validate
physical end state, or assert laboratory cleavage. The legacy digest recipe download
and descriptive end information remain available.

`scripts/catalog_digest_browser_qa.py` exercises real desktop controls and the fresh
WASM worker with independently authored GenBank fixtures and literal cut oracles:
all40/original12 controls, new-enzyme search, exact pin/geometry evidence, same-cut
identities, different-cut neoschizomers, five-base IUPAC circular-origin spans,
formulation conditions/qualifiers/citations/unknowns, annotated GenBank exports,
native exact-parent installation, forced second-put failure and pending transaction
abort, native/text/reload unverified evidence, limits/cancel/identical-content
source guards, and 1600/1280/1200 desktop layout. Actual browser results are recorded
in the milestone status after a fresh release build; the script is not a substitute
for an executed acceptance run.

## Worker release cache integrity

The main app and the complete analysis-worker loader/JavaScript/WASM chain must
belong to one build. Content-hashed main assets alone cannot establish that when
workers use fixed resource URLs. Release packaging therefore stages the complete
worker bundle under one content-hashed directory and declares its loader in the
index metadata; both isolated jobs and retained worker sessions use that loader.
This is a release-integrity pin, separate from scientific source/catalog bindings.

`scripts/worker_release_browser_qa.py` preserves the existing profile and HTTP
cache, arms invalid old fixed-URL responses, and checks fresh startup/reload with
actual constructor URLs, child-resource requests, schema21 envelopes, isolated
recovery saves, retained tracks, and a literal catalog digest. All fresh startup
exceptions and console/Worker errors remain checked. Missing/legacy loader metadata
must fail closed rather than falling back to a cached fixed-URL chain.
