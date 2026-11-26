# Source-import gateway prototype

This workspace MIT crate supplies a source-only native HTTP capability for a
Leptos frontend or Tauri host. It does not serve files, open workspace paths,
execute scientific programs, keep accounts, or change the static Python preview.
No server or live provider request has been started while writing it.

## Host API

- `Gateway::native(Config::default())`: loopback `127.0.0.1:8084`, exact matching
  frontend origin, two concurrent imports, four accepted HTTP connections.
- `Config::loopback_origins(...)`: explicit frontend origins for an owning local
  host. The POST route must be mounted **at the same origin as the frontend**.
  Opening another gateway port without a same-origin route does not fix CORS.
- `Config::remote(bind, origins, capability)`: explicit remote deployment policy,
  exact HTTPS origins and a 32–256 character bearer capability. No wildcard origin
  or credentials in URLs. The owning host/proxy supplies TLS and preserves Host;
  this crate never trusts forwarded-origin headers. Do not expose its plain HTTP
  listener as a public credential-bearing frontend.
- `handle_http(Request<Bytes>, Cancellation)`: reusable same-origin route handler.
- `serve(gateway, shutdown)`: optional bounded HTTP/1 listener for native hosting.
- `execute(Request, Cancellation)`: trusted native/Tauri host entry point; HTTP
  callers must use `handle_http` so host/origin/access checks cannot be skipped.
- `with_backend`: native construction-time fake/provider injection for in-process
  fixtures. Browser request data cannot replace a backend.

The route is POST `/api/source-import/v1` with content type
`application/vnd.collagene.source+postcard`. No CORS allowance, cookie, cacheable
response, or request/provider-body logging is added. Browser origins must match
both the configured origin and the current Host. Non-browser clients may omit
Origin but still require the correct Host, media type and configured capability.

## Shared binary API

`wire` compiles independently of the native listener on WASM. Request/reply
payloads use postcard inside a 42-byte integrity envelope: `CGSQ`/`CGSR` magic,
little-endian schema version and body length, and BLAKE3 checksum. Version1 rejects
future versions, changed hashes, trailing bytes and oversized payloads.
Serialization is bounded while allocating, not after an unbounded intermediate.
The checksum establishes integrity, not identity; use the owning same-origin TLS
host for remote operation.

The frontend must still check its current request ID, dialog generation and source
selection before installing an imported reference. Neither an import response nor
provider provenance creates a cloning verdict or prepared physical material.

A request has an ID and a fresh cryptographically random 32-byte client-session
nonce. Retain that nonce for one import dialog/session; do not reuse a public
constant. Responses retain the request ID. Operations:

1. `Import`: complete NCBI/ENA/DDBJ identifiers or their recognized official URLs.
2. `Catalog`: Addgene plasmid ID/page URL and an optional approved API token.
3. `Select`: returned opaque catalog handle, exact response hash, explicit sequence
   ID and the current optional API token. Signed download URLs are never accepted
   from the client or included in discovery summaries.
4. `Cancel`: cancels one active request within that client-session nonce.
5. `DisposeCatalog`: releases a reviewed catalog within that same client scope.

Addgene catalogs remain immutable native `AddgeneCatalog` values. The selected
record is downloaded using `SourceClient::import_addgene_selection` from the
reviewed catalog; no refetch silently changes the reviewed response or selects the
first sequence. The source client still performs its original download/identity/
sequence agreement checks. Expired or disposed handles require a new discovery
and explicit choice. Dispose old catalogs when a dialog closes or its source
changes; their deadline also removes idle cached catalogs through one weakly-owned
expiry task per gateway. No API credentials are kept in catalog state or request
logs. Request Debug omits input/token/client nonce; catalog Debug omits its handle.

## Resource and source boundaries

Requests are at most16KiB; replies at most64MiB with no partial document fallback.
Native application model guards also apply before returning an imported document.
The reused source client limits downloads to32MiB, uses 30-second HTTP requests
and a60-second redirect/rate-wait chain, and retains its DNS/redirect/credential
and parser checks. Native provider errors are sanitized and their previews bounded.

Gateway defaults: two concurrent native jobs, no waiting job queue, two retained
catalogs and a120-second catalog lifetime, with a70-second job deadline. Connection
headers/body have separate10/15-second deadlines. Explicit host limits can select
up to8 native jobs,16 connections,4 catalogs,90-second jobs and300-second catalog
lifetimes. Each native catalog is bounded by the source client's32MiB input and
its field/count limits; catalog memory is bounded by count, not inaccurately
estimated from the small public summary. Slow reply buffering is additionally
bounded by the connection count and64MiB reply cap. Choose host limits for its RAM.

Cancellation drops network/rate-wait futures and aborts the owned native task.
A synchronous parser already running cannot be forcibly interrupted; its bounded
native phase finishes with its semaphore permit still held. This prevents early
cancellation from admitting unlimited abandoned work. Dropping an execution future
aborts its child job. Native/Tauri hosts can provide explicit Cancellation tokens;
HTTP clients also have the owner-scoped Cancel operation.

Generic `Source::Direct` URLs are deliberately rejected. Initial-host filtering
would be insufficient because the current native Direct client may redirect to
other validated public hosts. Supporting arbitrary hosted GenBank files here
requires adding a complete redirect-origin policy to `collagene-sources`, not a
second fetch implementation or a misleading first-hop allowlist. Addgene's signed
links are different: they originate from its authenticated, reviewed native
catalog and still pass the existing public-address/credential-scope checks.

## Validation and integration still required

The package is registered in workspace members. Offline Cargo metadata resolved
the existing dependencies and added `httpdate1.0.3` for Hyper's server feature;
no dependency versions were refreshed.

Seven fixture groups are written for binary envelope bounds/redaction, explicit
remote configuration, scoped immutable catalog selection, cancellation/admission,
same-origin document exchange, rejected generic/private URLs and authorization.
They use only fake backends and in-process HTTP values. They have **not** been
compiled or run. Frontend/controller integration, serving-host integration,
actual native/WASM checks and deployment remain parent work.
