# Open source molecular cloning software comparison

Research baseline: 23 November 2026. Scope: cloning and sequence design applications and engines relevant to Collagene 2, with Windows/macOS/Linux, embedded processing, Rust/WASM, minimal TypeScript and MCP as comparison criteria.

Project decision: build Collagene 2 independently under the MIT license. None of these projects is selected as a codebase, UI foundation or source of cloning/verification algorithms. Their scientific correctness has not been established by this survey. Consultation is limited to a careful look at file-format handling when needed to understand undocumented formats; implementations and validation remain our own. The entries below retain the research findings, not a reuse recommendation.

There are active projects with substantial scope overlap. GENtle Rust is the closest architectural comparator for a shared deterministic engine and MCP. OpenCloning is a strong workflow and provenance comparator. PlasCAD supplies practical Rust prior art, while SeqNexus demonstrates a broad browser-local design. None of the inspected sources establishes a ready-made match for all Collagene 2 requirements.

Two parallel research agents covered applications and cloning engines, with additional browser-project research. This comparison uses primary documentation, dated release/package records and accessible commit histories. It is a research review, not a build, performance benchmark or scientific code audit. Dates below identify the observed evidence; cached histories may omit newer activity. A large feature list is not evidence of experimental validation.

## Closest application projects

| Project | Scope and implementation | License | Maintenance evidence | Fit for Collagene 2 |
| --- | --- | --- | --- | --- |
| [GENtle Rust](https://github.com/magnusmanske/gentle_rs) | Rust cloning workbench with a shared engine, GUI/CLI/MCP, workflows, PCR, Gibson planning, lineage and graphical explanations | GPL-2.0-or-later as stated in project documentation | [smoe fork history](https://github.com/smoe/gentle_rs/commits/main/) shows commits on 9 July 2026; [agent interface](https://github.com/magnusmanske/gentle_rs/blob/main/docs/agent_interface.md) dated 21 June 2026 | Closest AI/engine architecture. Desktop orientation and external-tool paths need review; Leptos/WASM suitability is unestablished |
| [OpenCloning](https://github.com/OpenCloning/OpenCloning) | Python/pydna backend and React/JavaScript frontend; cloning plans, common assemblies, primers, genome engineering and retained history | MIT for application repositories | [Cloning API package](https://pypi.org/project/opencloning/) 2.1.0 uploaded 6 November 2026; [organization](https://github.com/OpenCloning) lists frontend/backend updates on 2 November | Strong workflow/data-model reference. Python/JS implementation does not directly supply the requested embedded Rust stack |
| [PlasCAD](https://github.com/David-OConnor/plascad) | Rust/egui plasmid and primer editor; PCR, digest/ligation, maps, primer QC and AB1 viewing | MIT | Latest entry in the retrieved [commit history](https://github.com/David-OConnor/plascad/commits/main/) is 21 December 2025; [crate](https://docs.rs/crate/plascad/latest) 0.7.9 published 7 June 2025 | Useful Rust components and parsers. Windows/Linux downloads and macOS source builds documented; late-2026 maintenance not established |
| [SeqNexus](https://github.com/cacat-ayy/seqnexus) | Browser-local JavaScript application with workers/WASM; editing, assemblies, primers, traces, read analysis and exports | AGPL-3.0 | Public repository and [live feature documentation](https://seqnexus.app/features.html) available; exact latest commit/release date not verified | Strong browser UX/coverage comparator. New project; uses bundled WASM scientific tools and remote BLAST, so assess individual paths against our scope |
| [TeselaGen Open Vector Editor](https://github.com/TeselaGen/tg-oss) | React/JS/TS sequence editor and associated parser packages; circular/linear views, annotations, alignments and chromatograms | MIT | Current [OVE changelog](https://github.com/TeselaGen/tg-oss/blob/master/packages/ove/CHANGELOG.md) includes 0.8.42 on 16 April 2026 | Useful interaction, format and fixture reference. Direct UI reuse would introduce the JavaScript stack we want to minimize |
| [UGENE](https://github.com/ugeneunipro/ugene) | Broad C++/Qt bioinformatics desktop suite, including cloning, restriction/PCR analysis and sequence visualization | GPL-2.0 | Official [news](https://ugene.net/news.html) reports version 53.1 released 3 March 2026 | Established Windows/macOS/Linux comparator. Broad scope and external-tool integrations make it a poor direct architecture match |

### GENtle Rust architecture findings

Its [agent interface](https://github.com/magnusmanske/gentle_rs/blob/main/docs/agent_interface.md) documents `capabilities`, `state_summary`, `op`, `workflow`, `agent_preflight`, `agent_plan` and `agent_execute_plan` over stdio MCP. This establishes architectural overlap with AI-controlled cloning workflows, without establishing scientific correctness.

The [README](https://github.com/magnusmanske/gentle_rs) distinguishes implementation maturity from biological validation, with stronger support for core PCR/single-insert Gibson and caveats for other workflows. It links release artifacts to the smoe fork. The documented operation catalog does not establish complete method parity.

### PlasCAD implementation findings

The [project](https://github.com/David-OConnor/plascad) documents a narrower Rust application and uses an ecosystem including `na_seq`, `bio_files` and `gb-io`. Its [official site](https://www.athanorlab.com/plascad) describes partial SnapGene support, reinforcing the need for explicit import/export fidelity tests.

### Browser implementation findings

[SeqNexus](https://seqnexus.app/) describes local processing, worker jobs, canvas rendering and IndexedDB persistence. Its [alignment documentation](https://seqnexus.app/sequence-alignment.html) identifies MAFFT/MUSCLE/Kalign WASM builds; these are relevant prior art, not additions to our deferred alignment scope. A browser application can embed compiled scientific tools without spawning OS executables, but that does not establish that those tools satisfy the desired Rust implementation policy.

OVE's current development is in the [monorepo](https://github.com/TeselaGen/tg-oss). The old `TeselaGen/openVectorEditor` repository is archived; that alone does not mean its successor is inactive.

## Deterministic cloning engines

| Project | Documented scope | License and maintenance evidence | Runtime implication |
| --- | --- | --- | --- |
| [pydna](https://github.com/pydna-group/pydna) | Double-stranded DNA, ends, circularity, PCR and assembly; executable cloning documentation and verification | Current [pyproject](https://github.com/pydna-group/pydna/blob/master/pyproject.toml) declares BSD-3-Clause; [5.5.16](https://pypi.org/project/pydna/) released 9 June 2026 | Not selected for implementation or runtime use |
| [DNA Cauldron](https://github.com/Edinburgh-Genome-Foundry/DnaCauldron) | Assembly graphs, alternate products, hierarchical/combinatorial assembly, errors and reports | MIT; [2.0.12](https://pypi.org/project/dnacauldron/) published 12 May 2025; newer release not verified | Not selected for implementation or runtime use |

[pydna's documentation](https://pydna.readthedocs.io/latest/getting_started.html) connects provenance with the OpenCloning data model. Its overview labels some Golden Gate work as in progress, so the broad project description does not establish uniform method coverage.

## Additional leads and exclusions

- [GreenPeas](https://github.com/ryuabiru/GreenPeas) is a particularly relevant Tauri/Rust lead. The [author's announcement](https://www.reddit.com/r/bioinformaticstools/comments/1xgenfz/i_built_greenpeas_an_ultralightweight_opensource/) describes beta 0.1.0, canvas maps, digests, ligation, traces and Windows/Apple Silicon downloads. Repository contents/license and current activity could not be verified through the available fetches; Linux support was not established. Keep it on the follow-up list rather than classifying it as a verified foundation.
- [ApE](https://jorgensen.biology.utah.edu/wayned/ape/) remains a relevant active behavior/interoperability comparator; its site lists 3.1.11 on 19 November 2026. The [author's paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC9580900/) explains that redistribution of modified source is not granted. Treat it as source available, not a reusable open-source codebase.
- [AssemblyTron](https://github.com/PlantSynBioLab/AssemblyTron) centers on robotic execution and external design workflows, which are outside the agreed product scope.

## Implementation decision

Implement Collagene's UI, molecule model, cloning algorithms and verifier independently. Establish correctness using published scientific methods, manually derived cases, property tests and native/WASM parity tests. Other cloning software is not the authority for expected scientific results.

Where an undocumented file format blocks interoperability, a narrow parser inspection may establish layout and compatibility details. Record those observations and implement our own adapter with independent fixtures. Preserve the explicit deferrals for Primer3, multiple alignment engines and OpenGel.
