# Captured public INSDC records

These are unmodified public GenBank-format responses for accession **J02482.1**
from NCBI and DDBJ. `provenance.json` records source URLs, retrieval times, exact
byte counts and SHA-256 hashes. Original ORIGIN-line trailing spaces and the
NCBI response's final blank line are intentional captured bytes, not formatting
defects to normalize. ENA retrieval timed out; there is no ENA capture.

The active offline case in `provider_interoperability.rs` checks provider and
version identity, annotated native save/reopen, GenBank sequence/topology/feature
round trips, and matching normalized DNA across these two providers. It compares
reference HaeIII fragment lengths against the primary NEB digest listing linked
in the manifest.

NCBI labels the record ss-DNA. Mathematical double-strand restriction geometry
on that reference does not establish physical sample chemistry or a laboratory
digest. No physical-material authority is claimed.

This evidence exercises captured-record import, not live SourceClient transport,
CORS or provider availability. The separate live test is deliberately ignored
and remains deferred under the application HTTP acceptance restriction.
