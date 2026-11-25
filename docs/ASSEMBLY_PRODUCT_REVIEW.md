# Reconstructed assembly product review

`AssemblyProductView` connects fresh overlap/ligation evidence to the existing
WebGL map and virtual sequence renderer. Its caller supplies the actual typed
plan, independent target and freshly computed verification report. This is not an
imported-report or historical-positive-status entry point.

Before selecting any product, it checks the current engine and exact plan/target
revisions, and requires one reported product for every plan step. A failed final
step cannot display a leading preparation as the assembled result. If the final
reconstruction exists but disagrees with the target, its actual DNA remains
viewable with the invalid assessment. No expected sequence is copied from that
product, and display checks do not perform another scientific reconstruction.

The view offers sequence, map and split presentation; map rotation and zoom;
reference-coordinate navigation; pointer/keyboard range selection; annotations
with bounded paging; and explicit opening of the product reference. It reuses
the existing virtual renderer rather than expanding a full DNA string into the
DOM. Product selection is local and does not edit the source documents. The new
`SequenceView.read_only` option disables editing/history operations while keeping
selection available; the existing workspace lock still disables interaction when
required. Shared restriction controls are omitted because this snapshot has no
shared source scan to control.

Source/preparation summaries, findings and actual physical receipts are expandable.
A displayed duplex summary comes only from this fresh report and is never used to
restore material state or authorize a later reaction. The sequence/complement
tracks are reference-coordinate views; the separate physical receipt describes
available actual strand ends/nicks.

Junction inspection is opt-in and uses the existing private-worker reconstruction
of the exact plan and target. Its context binds both revisions and workspace epoch;
changing context cancels obsolete analysis. The inspection provides actual product
coordinates, homology/end evidence and requested local translations under the
existing model. It is not based on intended-source ribbon positions when actual
reconstruction differs.

The containing methods panel clears evidence when its source/settings/target
context changes. A workspace epoch change invalidates actions on an old snapshot.
Opening a product is explicit; retaining a complete producing run and its original
source-event graph uses the existing candidate verification/run workflow. Merely
opening a reference snapshot does not retain a physical graph or claim an
experimental product.

Browser interaction acceptance remains deferred until the in-scope implementation
is complete. Source/compiler checks do not establish that acceptance.
