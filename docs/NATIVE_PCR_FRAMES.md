# Exact PCR from real native strands

Current exact-PCR declarations explicitly select `template_strand:Top|Bottom`. This frame determines stored5′→3′ coordinates and product orientation; it does not crop the physical source, choose a product by an intended target, repair a template or establish material isolation.

Ordinary reference PCR uses explicit Top. Saved Run/Derived sources must freshly reconstruct their actual private tracked strands before a Bottom frame can be used. The native borrowed `pcr_source::Template` lends the existing top/bottom Sequence and exact source identity, not a substituted Molecule or new stock. The exact enumerator and primer ranking consume the same borrowed view. Coordinate receipts (`PcrProduct`, `DesignReport`) carry the selected strand explicitly.

For staggered sources, a diagnostic union search detects additional geometries before the real selected intact strand is used. The union is never a polymerase seed or input document. A geometry extending outside the selected covalent strand remains Unsupported. Selecting an actual bottom17 frame can enumerate both13- and17-base geometries; neither an intended target nor a selected contribution interval may suppress one. Nicked sources remain unsupported by this exact-PCR propagation model. Mismatch-PCR staggered-source physics is unchanged.

PCR products are synthesized as standard-unmodified blunt duplexes with the existing5′OH/3′OH model. Original source groups/modifications and all source calls remain unchanged in their native history. Source incompleteness/findings remain sticky. Annotation projection maps only real original top coordinates; Bottom maps reverse, and bottom-only extensions receive no fabricated top coordinates. Synthetic tails add only their explicitly declared output bases.

OE left/right first-PCR pairs carry their actual source frame; outer PCR explicitly consumes Top coordinates of the actual newly synthesized extension. Nonverifying OE compilation returns native preparation/search/oligo declarations and independently edited expected DNA, without a fake executable plan missing its native origin. The enclosing workflow supplies actual reconstruction operations and freshly verifies its complete plan. Ordinary standalone reference design builds and verifies its real original-source plan at that boundary.

One shared planner ledger charges actual source coverage scans and all materialized diagnostic/returned PCR products, including discarded union-search products. Origin replay, annotation construction and source-copy bounds remain their separately declared native limits.

The whole-construct planner uses the intent's single explicit material-binding
roster. Its shared resolved context retains each freshly reconstructed source
once for source-strand selection, independent target construction and primer
planning. Saved ReferenceRange contributions may also carry a native binding:
reading reference DNA does not turn it into physical stock, and preparation
still reconstructs that declared source. Target DNA is frozen before selecting
primer pairs or producing any preparation.

Source-declaration identity binds exact input revisions, events, closed origin
instructions and the scientific model. Ownership still requires the current
native workspace, but unrelated document heads and newly saved attachments do
not change the intended target's identity. Host request/adoption guards bind
their workspace scope separately. Capturing a batch and replaying it after
restart must therefore retain its original declaration and target revisions.

`RowDesign` retains the actual template strand, ranked calls and search notes.
There are no per-row original-reference verification badges. The composed
recipe's one fresh preview supplies the scientific assessment. Bounded
`Report.source_frames` provide original native sequence/frame/state evidence for
read-only correspondence and display; they cannot be deserialized into a
`FreshMaterial` or supplied to native PCR in place of origin reconstruction.
Rendering does not replay chemistry. Material-source tail calculations run in
the native worker, while ordinary reference-only tails retain a pure early check.

Library resolution uses the same cached source context, charges actual selected
frame/coverage work before preparing a row, and evaluates that exact resolved
declaration. Physical preparation reversal is independent of the source selector.
Primer order records retain the declared template strand and native step source;
exporting supplied oligos adds no new simulation or inferred top coordinate.

The whole planner's `ProfiledOverlap` final product is a reference projection.
Saving and freshly replaying that run verifies its declared sequence workflow;
it does not supply a reusable duplex. Its actual Bottom-PCR intermediate can
supply reusable material. The separately implemented physical `HomologyAssembly`
route retains its own processing assumptions and material-completion checks.

The `native_pcr_frames` fixture target records independently literal top13 `AGTCCGATGCTAC` and stored bottom17 `GTAGGTAGCATCGGACT`. It checks actual Bottom synthesis, source-top feature mapping to output4..17 reversed, native archive/reusable-product replay and standard product ends; explicitly competing13/17 products; literal tailed output; and a source whose21-base coverage union is not spanned by either intact17-base strand. Tests were authored/formatted without compiler or runtime execution during implementation.
