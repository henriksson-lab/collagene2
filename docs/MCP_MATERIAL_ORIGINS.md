# MCP saved material origin inspection

`inspect_material_origin` accepts an immutable session `WorkspaceRef`, an exact
native event binding, and an optional typed `MaterialSelection`. Omission lists
origins without choosing one. An explicit `Run` or `Digest` selection must name
the same event; native replay and exact workspace attestation precede any returned
closed origin. The complete bounded inspection is available through its immutable
JSON artifact; the response summarizes choices, findings and availability.

This read-only tool creates no document, stock, candidate or independent target.
A returned serialized origin is documentary. Later `FromRun` or `FromDigest`
preparations replay from original retained inputs before acquiring private
material. Unknown fields, stale workspace revisions and crossed source selections
are rejected. Existing run-only discovery/preparation tools remain compatible.

The SDK fixture imports an actual native saved-PCR/digest archive, discovers its
exact digest event, explicitly reconstructs the selected origin, checks the
literal product, and verifies malformed/stale requests leave the session unchanged.
