# Explicit FASTA alphabets and stable export headers

FASTA itself does not require Collagene's private alphabet marker. For example,
[NCBI's FASTA documentation](https://www.ncbi.nlm.nih.gov/genbank/fastaformat/)
uses sequence IDs/descriptions/source modifiers rather than our protein marker.
Letters such as A/C/G/T are valid protein calls too, so sequence composition
cannot safely establish an alphabet. Filename hints are not parser authority.

`formats::import_fasta` retains its established default-DNA behavior for unmarked
records and recognizes the exact compatibility tokens `[molecule=DNA]`,
`[molecule=RNA]` and `[molecule=protein]`. A header declaring different recognized
alphabets now refuses instead of giving protein accidental priority. Repeated
matching declarations do not change the type. The full header remains the native
name, preserving current callers' marker/name expectations.

`import_fasta_with_alphabet(text, Alphabet)` is the explicit profile for unmarked
RNA/protein or DNA inputs. Recognized header markers must agree with that choice.
A multi-record file containing a conflict fails as a whole; no partial import is
returned. Native sequence validation and existing source/record bounds remain.
No translation, transcription, T/U conversion or alphabet guess is performed by
this API. The ordinary default importer remains usable for existing DNA callers;
a caller opening external protein FASTA must choose the explicit profile.

Export always represents the actual native sequence alphabet. If a literal
native name contains a contradictory known marker, ONLY that token is escaped
to `[molecule=TYPE; label]`; a canonical actual-type marker is added or retained.
Correct headers remain unchanged, and ordinary unmarked DNA output stays
unchanged. Reimport/export is stable: escaped label tokens are not semantic
markers and further rounds do not accumulate markers or escaping. The source
native molecule/name/revision remains untouched. `export_with_report(...,Fasta)`
reports the name-header transformation explicitly as well as existing sequence-
only annotation/history omissions. Plain `export_fasta` remains its legacy
string convenience; it cannot return omission notices.

Five independent fixture groups cover DNA names containing a protein marker,
RNA with inherited contradictory labels, protein sequences also valid as DNA,
correct canonical legacy names, stable repeated round trips, atomic conflict
rejection, explicit unmarked profiles and existing malformed/empty-record rules.
Opening FASTA through the file picker shows a profile review and a bounded record
summary with alphabet and length. `.faa`/`.frn` supply editable filename hints;
reviewed declarations and explicit user choice determine the actual alphabet.
Record-declaration mode supports mixed explicitly marked records and shows that
unmarked records use DNA. A forced choice checks every record for agreement.
Multi-record FASTA is supported; files needing profile review are opened
separately before atomic import. The paste dialog also offers **Choose FASTA
alphabet**. No changes are committed until the reviewed records are accepted
against the same workspace/source context.

MCP `import_compat_sequence` accepts `Fasta { alphabet }`: `null` retains declared
per-record types/default DNA, while a chosen native `Alphabet` selects the strict
explicit profile. The existing GenBank text tool also accepts standard T-encoded
RNA through the corrected shared adapter. Browser checks remain deferred.
