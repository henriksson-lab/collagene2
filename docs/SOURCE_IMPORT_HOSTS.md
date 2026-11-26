# Source import in native, hosted and static applications

The existing SourceImport dialog keeps its original ImportReport callback. App
captures its workspace epoch when opening the dialog and checks that epoch again
before committing an imported document. A source response cannot bypass that
workspace guard or grant physical material authority merely by containing DNA.

## Native application

A native window always uses `source_gateway_client::Client`, including when its
IPC initialization is broken. It never falls back to browser fetch or requires a
lab login. Root provides `desktop_rpc::source_raw(bytes, before_dispatch)` over the
single `collagene_rpc` IPC command. CGSQ requests and CGSR responses retain their
existing 16 KiB / 64 MiB total frame bounds. Source traffic has its own native
lane so a provider request cannot block saving the working file or preferences.
Cancel and DisposeCatalog controls bypass that lane. For ordinary requests the
raw adapter checks the provided AbortSignal after taking the lane lock and before
invoking native work; this prevents a cancelled queued operation from starting.
It leaves late replies available for the client to dispose a returned catalog.

The native source service scopes the dialog nonce to native window authority out
of band. The client nonce is not a lab identity, login or filesystem capability.
Provider credentials are request-only dialog memory. No credential is added to
working archives, plasmid metadata, native credential vaults or remote bindings.

## Authenticated hosted application

A server enables hosted source retrieval by serving the marker:

```html
<meta name="collagene-source-gateway" content="enabled">
```

The dialog posts binary CGSQ to `/api/source-import/v1` on the same origin. This is
not a provider proxy accepting arbitrary HTTP requests: the native gateway owns
its existing provider/source allowlists, DNS/redirect policy, content and identity
checks, download limits and cancellation/catalog scope. The owning lab server
must authenticate and check Host, Origin and expected lab before reading any
request body, because an Addgene request may contain its provider API token.
Source responses use `application/vnd.collagene.source+postcard`.

`library_client::browser()` and the source dialog use one app-scoped
`Rc<Session>`. `Session::source_request` requires both its private lab bearer and
an established LabContext. Users connect to the Lab library and load/review its
context first. Source import never receives or exposes that bearer string; the
session adds Authorization and `X-Collagene-Expected-Lab`. Requests omit browser
cookies and reject redirects. A source dialog closing does not disconnect the
shared lab session. Library sign-out aborts its pending source fetches and late
responses from the old login are rejected. A lab identity change requires a new
login/context review; source requests never bootstrap a new context silently.

The response body is streamed with a 64 MiB cap rather than buffered through an
unbounded arrayBuffer. Requests have a 90-second browser timeout, leaving room
for the native gateway's default 70-second work bound and bounded response
transfer. The native gateway remains responsible for provider/network deadlines.
The client rejects unexpected response content types, failed HTTP statuses,
changed login epochs, mismatched response IDs and result variants. Provider/raw
HTTP failure strings are replaced with typed safe messages so a signed download
URL or provider credential cannot reach the error UI.

## Catalog and dialog lifetime

Each dialog creates a cryptographically random 32-byte client nonce and fresh
monotonically increasing request IDs. Addgene catalog discovery returns only a
handle, response hash, plasmid identity, expiry and bounded sequence choices.
Signed download URLs never cross this frontend gateway boundary. Choices retain
full/partial and depositor/Addgene classifications; no sequence is auto-selected.
Select sends the exact handle, response hash and chosen sequence ID for native
validation against the retained original catalog. The Addgene API token remains
in this dialog's password signal only; it is never logged, persisted or reused
for non-Addgene imports. Approved API access is still required. Downloading an
authorized file and opening it locally remains an alternative.

Input or token changes abort the active fetch/native queued operation, send a
scoped Cancel control and dispose retained catalogs. Changing a sequence choice
cancels an active selection while keeping the reviewed catalog. Closing the
dialog marks it closed, cancels work and disposes catalogs without touching the
lab login. If a catalog response arrives after cancellation, its handle is
explicitly disposed. Server TTL remains the backstop when transport/login loss
prevents the best-effort control from arriving. A fresh dialog uses a fresh nonce;
it does not replay prior requests or provider credentials.

The dialog's own generation guard rejects late UI results before invoking App.
App's existing captured workspace epoch remains the final commit guard. Catalog
expiry, provider rejection and transport failure remain explicit errors. There
is no switch to a different retrieval authority after a hosted/native failure.

## Static preview fallback

Only an ordinary browser page without the gateway marker keeps the original
DirectSourceImport implementation. It uses bounded CORS fetch directly to the
recognized provider/direct URL with no website cookies. Its existing 32 MiB and
45-second limits remain. Addgene credential-bearing requests are restricted to
the official API origin and reject redirects; browser CORS can still prevent
these requests. This fallback is for static preview, not a hidden workaround for
missing hosted authorization or broken native IPC.

The frontend files contain no HTTP/browser execution checks or provider traffic.
Host route tests, native service tests and the next unified WASM build establish
integration; enabling the marker alone does not create the authenticated route.
