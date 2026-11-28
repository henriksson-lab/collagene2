# Reusable overlap-extension PCR material

The ordinary Mutagenesis panel includes a collapsed **Overlap-extension PCR
material** section. Select the source and declare all six oligos: left source
PCR, right source PCR and outer PCR. Exact saved Run/Derived origins can be
selected with the existing source picker; bottom source strands require that
saved origin. Preview reconstructs the actual preparation. A separate explicit
review installs the reusable product through the existing native material
transaction. Assembly primer proposals continue to enter declarations into the
assembly draft and require ordinary complete-assembly review.

The native representation adds one `material_derivation::Derivation` variant:

```rust
OverlapExtensionPcr {
    preparation: overlap_extension_pcr::Request,
    selection: OverlapExtensionSelection {
        expected_preview_revision: String,
        expected_product_revision: String,
    },
}
```

The surrounding `material_derivation::Request` retains its mandatory exact
`source: EventBinding`, explicit `origin: SourceOrigin`, product name and primer
order choice. There is no target input, expected product sequence, status override
or imported physical-state object. The selected preview revision binds the exact
source document/event/revision, complete native origin hash, all three primer
pairs, selected source strands and extension/search options. The product revision
binds the actual singleton outer-PCR reference. Changing declarations requires a
new preview and explicit selection.

`overlap_extension_preview_revision(workspace, source_event, origin, preparation)`
and `preview_overlap_extension(...)` require the owning native workspace even for
original references. The core validates the exact source event and saved ancestry
before source cloning or scientific execution. Original-reference inputs reject
physical/reference-only markers; matching sequence content cannot replace an
origin or selected event. Historical preview does not waive installation's
current source-head and whole-workspace transaction guards.

One internal executor supplies source PCRs, tracked reciprocal extension,
annotation projection and outer PCR to both the existing OE digest inspector and
reusable material installation/replay. The digest inspector adds cleavage after
that shared execution. A verifier target is never fabricated to obtain material.
Neither serialized reports nor saved metadata can restore physical state.

Every source PCR and the outer PCR must expose exactly one supported geometry
before reusable material is available. Competing PCR products remain visible
stage diagnostics and withhold material. The declared terminal pairing is
retained separately from enumeration completeness/uniqueness; source or terminal
uncertainty stays incomplete through selection, installation and archive restart.
Selecting material is not evidence that alternative reactions cannot happen.

The final ideal-PCR product has newly synthesized unmodified blunt strands with
5′OH/3′OH. Source positional modifications are not copied into those new strands;
unsupported source nicks/coverage remain diagnostic. This model establishes no
amplification efficiency, cleanup, experimental isolation or intended-target
success. Source annotations follow the actual branch and outer-PCR mappings.

The native record keeps the complete recipe and exact source ancestry. Archive
reopen freshly reconstructs the source and all preparation steps, and compares
the saved product/state/incompleteness. Installed products can be chosen as
Derived origins for digest, PCR and assembly. Detached OE outer-PCR previews
carry restrictive physical-model markers and are not stock.

Primer export uses the existing native `OligoTable`, retaining six ordered rows
(left F/R, right F/R, outer F/R), full sequences and domain/tail decomposition.
Repeated names or sequences still retain their declared rows. Historical source
primers remain ancestry and are not silently added to this new order.

The record/capsule representation remains the current native material schema1:
the existing JSON recipe enum gains this operation directly, without another
record shape or compatibility reader. The application binary worker schema
advances with the new typed preview job/result. Root integration provides worker,
MCP and UI coverage; core fixtures live in
`crates/core/tests/materials/overlap_extension/overlap_extension_material.rs`. Validation results belong in the
implementation status after the integrated batch is checked.
