# GenBank nucleotide exchange

Collagene parses locations and feature qualifiers independently. The rich header
and sequence reader receives no feature-table data and only validated CONTIG
placeholders. Original records remain unchanged in native storage. GenPept has a
separate protein adapter.

The [INSDC feature-table specification](https://www.insdc.org/submitting-standards/feature-table/)
is the syntax reference. Exact local points/ranges, fuzzy range boundaries,
complement, join and order retain native location trees. Valid remote, between-base,
uncertain-point, gap, bond and one-of descriptors remain complete unresolved source
expressions. They cannot supply exact extracted sequence. Remote positions are not
compared against the local sequence length; known local positions are checked even
inside an unresolved mixed expression.

Parsing consumes the complete descriptor and uses checked positive coordinates.
Limits include 64 KiB per location, 32 levels, 4,096 parser work nodes, 1,000 records,
100,000 features, 16,384 qualifiers per feature, and matching file-wide location/qualifier
work bounds for reading and writing. Malformed later
records reject the import; no partial document collection is returned. CONTIG-only
records still require reference resolution and are unsupported as sequence input.

Qualifiers preserve order, repeated keys, flags, doubled quotes, literal Unicode
and backslashes. Quoted prose retains continuation newlines; translation and bare
values follow their existing continuation rules. Native machine fields such as
source hashes and primer sequences retain their explicit normalization rules.
Unknown biological declarations are source annotations, not validated claims.

Feature-table output uses an independent bounded writer with Unicode-safe wrapping
and atomic quote escaping. Its physical lines are at most 79 bytes. Feature keys
must fit the 15-column key field; qualifier keys must fit the declared output
prefix. Invalid fields, oversized output and unsafe unquoted continuations are
rejected. Position qualifiers require ASCII descriptors. Retain the native
workspace when a text profile cannot represent an annotation.

GenBank text does not preserve physical duplex material or native source ancestry.
Source comments remain unverified evidence. Native workspaces retain exact source
bytes, immutable events and the recipes needed for fresh material reconstruction.

Qualifier keys preserve explicit custom fields, including positional primer keys
such as `5prime_tail`. ASCII digits do not introduce a field delimiter; whitespace,
slashes, quotes and statement delimiters remain rejected.
