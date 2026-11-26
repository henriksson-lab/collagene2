# Annotated sequence exports through MCP

`export_annotated_sequence` exports one exact immutable `DocumentRef`. It uses the
existing native writer profiles; it does not add a conversion engine or restore
cloning authority. Registry tool count is 186.

| `format` JSON | Supported alphabet | Explicit options |
|---|---|---|
| `"GenBank"` | DNA, RNA | None |
| `"Embl"` | DNA, RNA | None |
| `{"GenPept":{"locus_date":"26-NOV-2026","entry_name":"LITERAL"}}` | Linear Protein | Calendar export date and stable entry ID; absent values use representable retained source metadata |
| `{"UniProt":{"entry_name":"LITERAL","accession":"P69905"}}` | Linear Protein | Actual caller-supplied or retained stable accession and entry ID |

No alphabet is inferred from residue letters. RNA GenBank exchange uses INSDC T
encoding while native RNA retains U. UniProt exports are marked Unreviewed;
accession assignment, reviewed status and unsupported molecular mass are never
invented. Stopped/unsupported protein states, unsupported annotations, invalid
calendar dates and irrelevant/unknown options remain errors.

`accept_omissions` defaults to false. When the native writer reports warnings, the
result is `state: "review_required"`, includes those warnings, and has no artifact.
Repeat the exact pinned request with `accept_omissions: true` after review to
obtain `state: "prepared"` and its immutable UTF-8 artifact. Acknowledgement does
not bypass a writer error. Output is limited to 32 MiB and total session artifact
quotas still apply; `read_artifact` pages the text. Native `export_document` keeps
its existing `.cg2` behavior. Native workspace archives remain necessary for exact
history, attachments and retained source evidence.

Focused SDK target: `annotated_sequence_exports`. It covers joined/circular and
reverse primer annotations, RNA T/U exchange, explicit synthetic protein metadata,
warning review, immutable sources, wrong alphabets, stale references and unknown
options. Compiler/test execution is delegated to the parent batch; no browser or
provider calls are part of these fixtures.
