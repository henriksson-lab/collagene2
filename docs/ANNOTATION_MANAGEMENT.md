# Annotation management

The Annotations inspector now has a searchable table with 100-row pages,
checkbox selection across pages, check-page/check-matching actions, Inspect/Edit
and New controls. Rows retain original source feature indices, exact location
text, type, review status and qualifier count. Batch selection is pinned to the
document/event context and revision; changing source clears it. Inspect still
shows all compound parts and the existing biological-sequence/qualifier details.

Apply changes the optional type, explicitly flips each checked feature's biological
strand, and/or applies a literal RGB color in one source-bound native operation.
Reset colors affects only the recognized color qualifiers. The host stages one
history commit and the resulting workspace before installing anything. Undo/redo
use the existing document history. A no-op creates no extra event.

Strand flip toggles the outer complement of the original location, retaining
compound/fuzzy/unknown geometry. It is not a command to make every mixed-strand
part point forward. Translation and codon-start assertions become source-evidence
notes, and the changed feature requires review. Unsupported coordinate-bearing
qualifiers block strand changes rather than being silently reinterpreted. Protein
features do not gain invented nucleotide strands. Unselected features, unfamiliar
qualifiers and their multiplicity, sequence calls and original evidence remain
unchanged. The current native annotation-edit model budget permits 4,096 features,
16,384 aggregate location nodes/qualifiers and the existing text limits; large
read-only display limits are separate and do not imply unrestricted scientific
edits.

Only literal `#RRGGBB` colors are interpreted. Canonical color and imported ApE
forward/reverse colors feed the map, label and sequence renderers. Conflicting or
unsupported color text falls back safely to type colors. Per-segment rendering
can honor distinct strand colors. Changing type or color preserves other qualifiers.

Show/hide actions are display preferences, not sequence edits or deletions. Hidden
features retain original IDs and consume no annotation lanes; their geometry,
labels and hit targets disappear consistently from map and sequence. Scientific
worker reports and exports retain the annotations. Selection can still highlight
reference bases independently of an annotation being hidden.

Display preferences are content-hashed native-workspace attachments bound to the
exact document, retained event, revision and feature count. Reopening validates
ownership and index bounds. Undo can recover a prior event's display preference.
A batch type/strand/color edit preserves indices, so its preferences are explicitly
carried to the new event during the same staged installation. Other source edits
start with all annotations visible; an old index must not accidentally hide a new
annotation after insertion/deletion or unrelated import. When no indices remain hidden, that
event's preference attachment is removed. Native attachment limits apply.

Implementation: `core::annotation_batch`, `web::annotation_table`,
`web::annotation_visibility`, shared map/sequence rendering and layout. Focused
native fixtures cover atomic edits and qualifier evidence, exact-event visibility
and native archive roundtrip, lane release and preserved original indices. Local
browser and performance acceptance remain deferred until full in-scope code is
complete.
