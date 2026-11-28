# Importing sequence records by URL

The `collagene-sources` crate keeps source resolution and networking outside the scientific core. A downloaded sequence passes the same bounded format parsers as a local file. The native client embeds Rust HTTP/TLS libraries; it invokes no downloader, browser, shell or scientific executable.

## Supported routes

- NCBI `nuccore`/`nucleotide` record URLs, accession.version identifiers, and explicit `ncbi:` identifiers resolve to documented EFetch requests for `nuccore`, `gbwithparts`, text, and a single bioseq. Numeric IDs are accepted only in explicit NCBI URLs/prefixes; bare positive numbers mean Addgene plasmid IDs. NCBI search, multiple-ID and subsequence/reverse-strand URLs require an explicit record selection instead of silently ignoring those choices.
- Direct public HTTPS sequence-file URLs are fetched and content-sniffed. Currently recognized data includes GenBank, FASTA, EMBL, SnapGene DNA, ABIF/SCF and native Collagene documents. Existing parser limitations still apply. HTML/login/JSON error pages are not passed to a DNA parser. Compressed downloads and arbitrary plain-text DNA responses are not accepted by this URL adapter.
- Addgene plasmid page URLs/IDs resolve through the official authenticated `catalog/plasmid-with-sequences/{id}/` endpoint. The catalog lists all four documented deposited/Addgene and full/partial sequence categories. An explicit `sequence_id` selection is required, even when only one record is available. The selected GenBank file uses the provider's `genbank_api_url` or documented deprecated `genbank_url` alias; its bases/length are cross-checked against catalog sequence evidence when supplied. No endpoint, response schema or automatic sequence choice is invented. Official temporary links can also be imported as direct HTTPS files, but direct-file import cannot infer catalog identity/categories.

The Addgene adapter is implemented from the [official developer documentation](https://docs.developers.addgene.org/docs/) and its [linked OpenAPI schema](https://docs.developers.addgene.org/docs/schema/), retrieved 2026-11-23 (schema build `b2c720c11f`). Native code sends `Authorization: Token KEY` only to `https://api.developers.addgene.org`; catalog discovery requires the approved `catalog:retrieve-with-sequences` scope. Sequence downloads may redirect to a signed public file URL, where the token is omitted. Catalog redirects cannot leave the API origin. Missing/expired access, unavailable sequences, changed schema and mismatched sequence data remain explicit errors. No authenticated live account has been used to validate provider behavior; routine fixtures are independently authored against the published schema.

## Native API and configuration

```rust,no_run
use collagene_sources::{Config, SourceClient};

# async fn example() -> Result<(), Box<dyn std::error::Error>> {
let config = Config::from_env()?;
let client = SourceClient::new(config)?;
let report = client.import_url("NC_001422.1").await?;
# Ok(()) }
```

`Config` supports private constructor/builder credentials and optional `COLLAGENE_NCBI_API_KEY` / `COLLAGENE_NCBI_EMAIL` environment variables. It is neither serializable nor debug-printable. Configure the user's contact identity according to NCBI guidance; no contact email is invented. Configured keys/contact are sent only to the fixed NCBI EFetch origin, never to direct-file providers or another redirect origin. They are absent from molecule metadata, public portable fetch URLs and error messages.

Addgene uses `Config::with_addgene_api_token(...)` or `COLLAGENE_ADDGENE_API_TOKEN`. Discover and select separately:

```rust,no_run
# use collagene_sources::{Config, SourceClient};
# async fn example() -> Result<(), Box<dyn std::error::Error>> {
let client = SourceClient::new(Config::from_env()?)?;
let catalog = client.addgene_catalog("https://www.addgene.org/12345/").await?;
// Present catalog.choices() to the user. Do not choose by array order or first full record.
let selected_sequence_id = 42; // replace with the user's explicit choice
let report = client.import_addgene_selection(&catalog, selected_sequence_id).await?;
# let _ = report;
# Ok(()) }
```

`import_url_selected(input, Some(sequence_id))` refreshes the catalog before retrieving that selected record. The existing `import_url(input)` returns `AddgeneSelectionRequired` with the count, or `AddgeneNoSequences`, after an authorized catalog discovery. `AddgeneCatalog` is an immutable parsed snapshot: `plasmid_id()`, `name()`, `choices()` and `response_blake3()` provide read-only access, with no mutable accessor. Choice records may be cloned for presentation; modifying an owned clone does not relabel the retained import binding. Serialization keeps the same JSON fields (`plasmid_id`, `name`, `choices`, `response_blake3`) and debug output contains only public summaries and response hashes; transient signed URLs and raw bases remain private. `download_url(sequence_id)` deliberately returns a transient URL for a browser adapter, never provenance.

The native client uses a process-wide 350 ms NCBI request interval, conservatively below three requests/second even with a configured API key. It surfaces `Retry-After` and applies NCBI 429/503 cooldowns. There are no automatic HTTP retries. Downloads are streamed with a 32 MiB limit, 30-second request timeout, 60-second whole-chain timeout and at most five redirects. HTTP status failures, partial bodies, missing records, returned identifier mismatches and ambiguous NCBI responses are distinct typed errors.

Addgene API calls share a conservative process-wide one-second interval and respect 429/503 `Retry-After` cooldowns; this is a client policy, not a claim about a provider contract. Catalog and selected-file retrieval each have the same bounded chain timeout. Catalogs are capped at 1,000 choices and 32 MiB; duplicate IDs, unknown sequence categories and invalid lengths/locations/URLs reject the response rather than hiding choices. Nullable lengths and empty catalog bases remain explicit missing evidence. Partial records cannot be imported as circular complete references.

HTTPS certificate verification remains enabled. The native transport disables environment proxies, cookies and automatic redirects. Each redirect is revalidated; DNS resolution rejects non-public or mixed public/private address sets before connecting. IP literals, transition ranges, local hostnames, credentials in URL authority and nonstandard ports are rejected. This conservative policy also excludes some globally reachable special-purpose ranges. It is intended for public sequence sources, not intranet/LIMS access. It does not claim an authenticated hosted endpoint or project authorization; those remain separate application responsibilities.

## Browser API

`parse_source`, `Source::download_url`, `Source::canonical_url`, `detect_format`, `annotate_import` and `import_download` compile to WASM. A browser can fetch a supported public resource where its CORS policy permits, parse it with `formats::import_bytes`, then call:

```rust,no_run
# use collagene_sources::{parse_source, annotate_import};
# use collagene_core::formats::sequence as formats;
# fn example(raw: &[u8]) -> Result<(), Box<dyn std::error::Error>> {
let source = parse_source("NC_001422.1")?;
let transient_fetch_url = source.download_url()?; // never save/log signed query parameters
let mut report = formats::import_bytes(raw)?;
annotate_import(&mut report, &source, raw, "2026-11-23T18:00:00Z")?;
# let _ = transient_fetch_url;
# Ok(()) }
```

Native DNS/IP enforcement cannot be reproduced by ordinary browser fetch. A browser CORS/network failure is not resolved by changing the URL parser or exposing credentials; protected/non-CORS providers need an authenticated native/local or hosted adapter. The source helpers do not create that bridge, issue Addgene tokens, or bypass authentication.

The portable Addgene API comprises `addgene_catalog_url`, `parse_addgene_catalog`, the safe catalog summary, and `catalog.import_download`. The web **Import URL** dialog accepts an optional password-type token kept only in component memory. **Import source** first discovers the catalog, then displays an unselected sequence list; choose a record and click again to retrieve it. Editing source/token aborts the request, discards pending responses and clears prior choices. Authorization-bearing browser requests reject redirects because browser fetch cannot inspect cross-origin redirect locations safely; use native/MCP retrieval if the provider returns its authenticated redirect route or lacks CORS. Signed external download URLs carry no Authorization header. No login cookies are sent. A token alone does not resolve browser CORS: the live localhost preflight described below failed, so the current observed route requires authorized native/MCP retrieval or a file personally downloaded through the provider account.

## Provenance and evidence

Each imported molecule records meaningful `source:provider`, `source:canonical_url`, `source:identifier`, `source:accession_version`, `source:retrieved_at`, `source:content_blake3` and `source:format` entries. Retrieval times are validated RFC3339 and normalized to UTC. Canonical direct URLs omit the entire query and fragment. Content hashes cover fetched bytes; a record index distinguishes multi-record direct files. A direct file without a declared database accession uses the content hash as its source identifier.

NCBI responses must contain exactly one GenBank record with an actual `VERSION` accession. A requested version must match; an unversioned accession is resolved to its returned version and canonical versioned URL. GenBank record terminators and declared sequence lengths are checked in addition to the shared parser. Provenance does not assert experimental function or complete-plasmid identity.

Selected Addgene imports additionally retain `source:sequence_id`, `source:addgene_origin`, `source:addgene_completeness`, the sequence description, the fetched catalog response hash and a bounded selected public-fact JSON snapshot. Full/partial/deposited/Addgene categories are provider assertions; they are not inferred from sequence length or treated as independent experimental verification. Original GenBank evidence remains intact. API tokens and presigned link strings are absent from the portable summary, saved provenance and diagnostics. The catalog response is hashed but not copied verbatim, because it contains expiring signed links.

Original GenBank/EMBL/SnapGene/trace evidence remains under existing format-parser ownership. FASTA imports additionally retain each exact original record as base64 and a record hash in native metadata. All resulting documents must still fit the native file payload limit. FASTA has no declared biological length: HTTP framing can detect a partial transfer when a length is supplied, but cannot prove a biologically complete plasmid. Multiple direct-file records are returned explicitly; they are not silently reduced to one. Re-import should add a revision/document through application transactions, never overwrite an existing design implicitly.

Native `.cg2` persistence retains full provenance and evidence. Derived GenBank/EMBL/FASTA export preservation is implemented by the core/application format adapters; these URL helpers do not reconstruct an original remote document from edited bases. A canonical URL may not be directly fetchable after a signed link expires; the content hash and retained evidence still identify what was imported.

## Validation and remaining work

```sh
cargo test -p collagene-sources --locked
cargo check -p collagene-sources --target wasm32-unknown-unknown --locked
cargo run -p collagene-sources --example import_url -- NC_001422.1
```

Offline tests use independently authored GenBank/FASTA fixtures and a mock HTTP transport. They cover source resolution, actual accession pinning, source/feature fidelity, query-secret redaction, original FASTA evidence, multiple/partial responses, HTML and expired links, redirects, public-DNS filtering and request spacing. The development example performs an optional live public read; routine tests do not require a live provider. A local live run against `NC_001422.1` succeeded and returned 5,386 bases, circular topology, and the matching canonical versioned NCBI URL. Native tests and the portable WASM check also pass. An isolated mocked-browser suite verifies explicit catalog selection, origin-scoped Authorization, aborting an old in-flight choice, partial-record labels, retained annotated exports and absent tokens/signed links in saved evidence. This does not replace a live authenticated interoperability check.

Persistent content-addressed response caching, conditional HTTP retrieval, live authenticated Addgene interoperability checks, per-project source authorization and browser-native bridging remain future work. The current native client fetches fresh data for each import; its evidence hash is not a freshness cache. Offline Addgene tests cover the published four-category schema, explicit IDs, original annotated GenBank evidence, partial records, nullable/empty evidence, schema changes, duplicate/mismatched records, origin-scoped token forwarding and expired/rate-limited/redirected responses.

The implementation uses [NCBI's E-utilities reference](https://www.ncbi.nlm.nih.gov/sites/books/NBK25499/) for request and identity parameters. [Addgene's developer portal](https://developers.addgene.org/access-options/) documents approved access options, and its [official MCP help](https://help.addgene.org/hc/en-us/articles/48534756585869-How-do-I-use-the-Addgene-MCP-Connector) describes the provider-managed connector. Public-destination policy follows a conservative subset of the [IANA IPv4](https://www.iana.org/assignments/iana-ipv4-special-registry/) and [IPv6](https://www.iana.org/assignments/iana-ipv6-special-registry/) special-purpose registries. No cloning application code was copied or used as a scientific correctness oracle.

## Application integration

The web toolbar **Import URL** uses bounded browser fetch and the portable resolver/provenance adapter. The MCP `list_source_sequences` tool lists bounded Addgene choices plus a full safe-summary artifact; `import_url` accepts an explicit `sequence_id`, uses `SourceClient` and returns immutable document references. If browser CORS prevents retrieval, call MCP `export_document` with the returned `{id, revision}`, read/concatenate/decode its base64 artifact and open the `.cg2` file in the browser. This preserves the exact retrieved revision and source evidence; separately re-fetching a record creates a distinct provenance revision. Browser imports time out after 45 seconds; native client limits are described above.

## Subsequent provider adapters

ENA/DDBJ adapter foundations and independent offline fixtures are now present in the source crate. Dedicated provider expansion, live interoperability and further source UI work are deferred while the desktop viewer/GenBank/virtual-digest MVP is prioritized. The routes below distinguish those foundations from the remaining SEVA/iGEM candidates. Research on 2026-11-24 found:

| Provider | Official route and implementation boundary |
|---|---|
| ENA | The [Browser API](https://ena-docs.readthedocs.io/en/latest/retrieval/programmatic-access/browser-api.html) documents `/ena/browser/api/embl/{accession}` and `/fasta/{accession}`. Prefer annotated EMBL for a single sequence accession/version. Project/study/assembly pages may identify collections and require a separate selection or reference-store workflow. |
| DDBJ | [getentry](https://www.ddbj.nig.ac.jp/services/getentry-e.html) documents `https://getentry.ddbj.nig.ac.jp/getentry/na/{accession.version}?format=flatfile&filetype=text`, with FASTA as an alternative. Validate the annotated flatfile against parser fixtures, preserve actual versions and reject ranges initially; versions in accession ranges are ignored by the provider. |
| SEVA | Resolve exact plasmid IDs through the [official catalog](https://www.seva-plasmids.com/v2/table-all.php), then follow its actual NCBI accession or GenBank link. Do not construct filenames from plasmid IDs. No public API was established by this inspection; direct GBK payload interoperability remains untested. |
| iGEM | The [rebuilt Registry](https://registry.igem.org/about/whats-new) advertises annotated GenBank, FASTA and SBOL exports. Current parts use pages such as `https://registry.igem.org/parts/BBa_B0015`; an [API documentation shell](https://api.registry.igem.org/docs) exists. Exact export endpoints and anonymous API access still need verification. Do not assume the legacy XML interface describes the rebuilt service. |

Browser CORS remains unverified for these routes; public pages and server-side retrieval do not establish cross-origin browser access. Preserve provider identity/version, retrieval evidence and attribution, keep multiple/partial records explicit, and disclose annotation loss for sequence-only fallbacks. Large remote chromosomes must stream into the separate reference store rather than bypass ordinary document limits. No other cloning application's source code was consulted for these adapters.

## Addgene anonymous-access audit (2026-11-24)

Small, unauthenticated read-only requests to official Addgene hosts established a
current access boundary. No account, token, cookies, protected file URL guessing,
or authenticated session was used. This is an access/interoperability audit, not
a successful live sequence import.

| Official resource | Observed response | Publicly visible evidence |
|---|---|---|
| `https://www.addgene.org/52961/` | 200, HTML, 237,587 bytes | Canonical plasmid identity/name, map and sequence-list navigation; sequence access login notice. |
| `https://www.addgene.org/52961/sequences/` | 200, HTML, 214,958 bytes | Ten analyzer links: Addgene full 1 and partial 9. Full sequence ID 244694 is publicly listed, but no GenBank or DNA download anchor is offered anonymously. |
| `https://www.addgene.org/125752/sequences/` | 200, HTML, 206,467 bytes | Addgene full 1 (ID 249038) and depositor full 1 (ID 240779), with the same login requirement and no anonymous sequence-file link. |
| `https://api.developers.addgene.org/catalog/plasmid-with-sequences/52961/` | 401, JSON, 58 bytes; `WWW-Authenticate: Token` | No anonymous catalog/sequence response. |

The public sequence-list sections distinguish Addgene/depositor and full/partial
categories, but analyzer links and map images are not sequence-file authority.
Both sampled lists and the linked analyzer pages explicitly require login for
sequence access. Addgene's [official access help](https://help.addgene.org/hc/en-us/articles/44210206209549-Are-there-any-requirements-for-viewing-plasmid-sequence-data-on-the-Addgene-website)
confirms this boundary; its [August 2026 explanation](https://blog.addgene.org/protecting-the-science-you-share-why-sequence-data-requires-a-login)
describes the change from historical anonymous downloads. Its
[developer access options](https://developers.addgene.org/access-options/) require
approved scopes and data-access license acceptance for programmatic sequence
access. Consequently, constructing alternative backing-file URLs from public
map filenames is not an implemented or validated anonymous import strategy.

The native responses above had no `Access-Control-Allow-Origin` header when sent
an `Origin: http://localhost:8080` header. A separate API authorization preflight
(`OPTIONS`, origin `http://127.0.0.1:8080`, requested method `GET`, requested header
`authorization`) returned 401 with `WWW-Authenticate: Token`, allowed methods
`GET, HEAD, OPTIONS`, and no CORS permission header. This preflight fails for that
browser origin. It does not establish how an approved authenticated response or
another explicitly permitted origin behaves. The native Rust client is unaffected
by browser CORS; the browser dialog's mocked success is not evidence of live
cross-origin API interoperability.

Audit retrieval times were 2026-11-24T13:04:59Z–13:06:28Z. SHA-256 identifiers of
the exact observed response bodies (not redistributions or authenticity signatures):

- Plasmid 52961 HTML: `25d2a65c1f0a1731e9cbe8ddf59846d218d73f25ba90a368013b748112b3c506`.
- Plasmid 52961 sequence-list HTML: `be6c81347b7cc9d3f8d5c9aa4bcce8b6f4645a0b3603fbc613bb7906075f4a35`.
- Plasmid 125752 sequence-list HTML: `9b379cbb2db621fa7a1ee9478cbf18f36698a3892a4dd4c94af714fe8695471b`.
- Unauthenticated API 401 body: `dbff0849bc109eee7f9b7b1fe12e40fc8add21823b2f1fb600b37f620cbd7862`.

A safe future public-page adapter could expose bounded, hashed **metadata only**:
canonical plasmid identity, public sequence IDs, provider category, analyzer URL,
and an explicit login-required capability state. It must remain separate from
`AddgeneCatalog`'s authenticated sequence/download evidence and cannot manufacture
a download URL, base/length agreement, or verified file identity. Public-page
schema changes, absent identity evidence and unavailable CORS must fail explicitly.
There is no automatic fallback from a rejected protected API call to a different
public sequence or another sequence category.

The current safe sequence-import choices remain the approved native API with an
explicit catalog sequence selection, or a GenBank/SnapGene file personally
downloaded through the provider's normal account flow and opened locally. A
user-supplied official temporary file link remains a bounded direct-file import;
its filename alone cannot establish plasmid/category identity. Live authenticated
API/download behavior still requires an authorized interoperability check.
