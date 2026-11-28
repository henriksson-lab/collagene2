# Reaction choices in cloning plans

In **Virtual digest → Reaction conditions**, opt in to an explicit setup. Choose
the exact formulation for each enzyme, buffer, temperature, duration and required
recorded activity threshold. Assay contexts distinguish records made under different
conditions. Missing formulation selection means unknown; an HF formulation is not
substituted automatically.

The same editor is available within a catalog-pinned digest preparation in the
assembly workflow, including PCR, mismatch PCR and overlap-extension PCR followed
by digestion. Incomplete fields remain editable and can be saved in the current
draft. They must parse and match the exact enzyme selection before compiling the
plan. An enabled setup on another preparation type must be explicitly disabled or
the preparation changed; it is not silently discarded.

**Use in active assembly slot** carries the current setup along with the exact
source and chosen fragment. **Open fragment** retains it in the native derivation.
Saved-origin preparations retain that ancestry and separately declare any new
digest setup. Editing a setup while an open-fragment operation is pending invalidates
the operation's guard. Invalid setup fields prevent transfer rather than becoming
an omitted setup.

Fresh verification assesses the declared choices against the pinned catalogue and
the actual source entering the digest. For a PCR/digest chain, that source is the
newly reconstructed PCR product. Source and physical-material identities are retained
with the documentary assessment. Replay recomputes the evidence from retained
inputs; imported reports alone do not establish a current assessment.
Ancestry follows actual typed dependencies of the selected product. A disconnected
earlier digest is excluded, and a shared ancestor is retained once.

Read sequence reconstruction and condition evidence separately. Unknown, conflicting,
below-threshold or otherwise unresolved declared conditions leave the corresponding
constraint incomplete while retaining available product geometry and independent
target comparison. Satisfying the exact documentary request does not prove cleavage,
reaction compatibility or experimental success. Recorded unit-assay activity is not
percent cleavage. Star activity and incubation-duration suitability remain unassessed.

HTML and PDF protocols retain choices, source/formulation citations, contexts and
the documentary limitations. Originating assessments are identified as provenance;
they do not automatically describe a later reaction. Plans with no setup say so.
Uncatalogued digest and Golden Gate mixture steps explicitly state that these
supplier-condition choices are unassessed by those operation profiles. Target-free
PCR/digest inspection remains geometry inspection; declare the setup on its plan
preparation for retained condition assessment and protocol output.

The catalogue currently covers a limited set of exact formulations and assay facts.
Broader curated coverage, source-host presets and appropriately supported star facts
remain on the [current plan](PROJECT_PLAN.md). No external scientific executable is
used. See [design checks](DESIGN_CHECKS.md) for the separate construct-QC, assembly
format and oligo-collection tools.


## Documentary source hosts

Original-source generation controls may retain one exact catalogue host preset.
The declaration includes the catalogue edition/content hash, full preset and
cited source records; editing any pin or source record rejects current admission.
Host expectations are strain/product-specific documentary evidence. Selecting a
host does not change the independently declared Dam, Dcm or CpG chemical state.
CpG is unknown for all eight current E. coli presets. Unknown and deficient
activity are distinct. Other modification systems remain outside the closed
source-generation model. Native history and HTML/PDF protocols retain the
declaration; imported metadata does not establish actual methylation.
