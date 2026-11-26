# Native desktop credentials

`crates/desktop-credentials` implements the existing
`repository::http::CredentialProvider` boundary using real native credential
services. It stores no access token in SQLite, configuration, scientific archives,
remote bindings, outbox records or synchronization frames. It is a native-only
service component; this change does not add a shell command or alter IPC schemas.

## Trusted scope and public API

Construct `Credentials::new(local_repository_lab)` using the lab ID derived from
the actual trusted local repository. A caller-provided identity is not authority.
Each `RemoteConfig` is resolved by native repository state before credential use.
The account name is a domain-separated BLAKE3 hash binding, in order:

1. the local repository lab ID;
2. the stable `rem-` remote configuration ID;
3. the canonical HTTP(S) origin, including a nondefault port;
4. the expected remote lab ID.

Each field is length-prefixed. Origin spelling/case and an explicit default port
normalize consistently; a renamed display label does not change the scope.
Changing any security pin changes the account. Fixed application/version service
names and opaque account hashes avoid placing raw origins, lab/member names or
credentials in vault labels. There is no wildcard lookup or arbitrary service
name API.

Origin parsing rejects user information, paths beyond `/`, query strings,
fragments and non-HTTP(S) schemes. Computing an HTTP account is not permission to
send it over HTTP: `HttpRemoteTransport` independently enforces trusted HTTPS,
or its explicitly enabled literal-loopback exception. The credential adapter
does not make hostname guesses or change transport policy.

The API provides:

- `Secret::new(Vec<u8>)`, consuming and validating a token;
- `Credentials::store(remote, secret, Persistence::Session | Remember)`;
- token-free `status(remote)` returning `Session`, `Remembered` or `Missing`;
- exact-scoped `delete(remote)`, `forget_session(remote)` and `clear_session()`;
- `bearer(remote)` or the existing `CredentialProvider::bearer_for` method;
- nonsecret `local_repository_lab()` for a native service scope assertion.

`Secret` has no `Debug`, `Display`, `Clone` or serde implementation. It zeroizes
its owned buffer on drop, including rejected input. A backend explicitly borrows
its bytes during an OS call. OS, D-Bus and HTTP libraries necessarily manage
additional transient copies; this is not a claim that every native allocation
can be erased by the application. Tokens contain 1–2048 printable ASCII bytes,
with no whitespace/control/header injection. The common bound also fits Windows
Generic Credential blob limits. Current 64-hex-character lab keys fit naturally.

`HttpRemoteTransport::new(Credentials::new(local_lab), options)` provides the
production transport. Its `credentials_mut()` method supports native status,
store, clear and delete operations. `Backend`/`with_backend` permit controlled
headless fixtures and service injection without using a real vault.

## Session versus remembered credentials

Construction is lazy and performs no OS access. Session mode never contacts the
vault, so it works when Secret Service/Keychain/Credential Manager is unavailable.
Memory is scoped by the same account identity and is erased when cleared or the
provider is dropped. No file fallback exists.

Remembering is an explicit operation. A locked, denied, ambiguous or unavailable
vault produces a typed, sanitized error. Failure removes an active session entry
for that scope and does not install the new secret in memory as a hidden
fallback. A user can deliberately choose Session mode afterward. Retrieval uses
an explicitly active session entry first, otherwise the exact OS entry; it never
looks up another remote or repository as a substitute.

Delete removes the matching active session entry before attempting OS deletion.
Missing OS entries count as already deleted. A locked/unavailable vault means
persistent deletion was **not confirmed**; the error remains explicit. Clearing
session memory alone does not delete remembered credentials. Status reports no
token bytes and does not turn an OS failure into a false `Missing` result.

Backend errors are reduced to `Missing`, `Locked`, `Denied`, `Unavailable`,
`Ambiguous`, `InvalidScope` or `InvalidSecret`, without native error descriptions,
labels, paths, secret text or backend response bodies. The transport trait reduces
these further to its generic redacted `CredentialError`.

## Platform backends

Windows uses `CredReadW`, `CredWriteW`, `CredDeleteW` and matching `CredFree` with
`CRED_TYPE_GENERIC`. Persistence is local-machine/current-user rather than
enterprise roaming. Exact service/account targets are NUL-terminated UTF-16;
read size/type are validated before copying. Returned native blob data are
cleared before release. These APIs use the logged-on user's credential set and
do not invoke an external program. See Microsoft's
[CredWriteW reference](https://learn.microsoft.com/en-us/windows/win32/api/wincred/nf-wincred-credwritew).

macOS uses Security framework generic-password `SecItem` calls for the exact
service/account, explicitly non-synchronizable. `kSecUseAuthenticationUIFail`
forbids an implicit unlock/authentication dialog. A locked or interaction-required
operation fails instead of prompting or storing elsewhere. Core Foundation
objects are owned for the complete native call and returned data are checked
before copying. Apple documents the
[Keychain services API](https://developer.apple.com/documentation/security/keychain-services)
and [no-UI policy](https://developer.apple.com/documentation/security/ksecuseauthenticationuifail).

Linux uses the existing `org.freedesktop.secrets` service over an existing local
Unix D-Bus session socket. It rejects TCP, `unixexec`, autolaunch and other bus
transports before connecting; there is no external executable or daemon-launch
fallback. Each operation captures the existing service's validated unique bus
owner, then addresses every session/item/property call to that same owner with
`NoAutoStart`; owner replacement cannot silently receive a pending secret.
Every proxy has `CacheProperties::No`, and property reads use explicit
`Properties.Get` with `NoAutoStart` rather than implicit `GetAll`/`get_property`.
Calls are
bounded by a five-second method timeout and a small message queue. The standard
`plain` session algorithm transfers bytes only over this authenticated local bus;
it is not application-layer encryption against the bus administrator. The OS
Secret Service owns protected at-rest storage. The
[Secret Service specification](https://specifications.freedesktop.org/secret-service/latest/)
defines the session, collection, item and prompt protocol.

Remembering uses an unlocked default collection; retrieval/deletion searches only
the exact application/account attributes. Locked matches or multiple unlocked
matches fail explicitly. The adapter never invokes `Unlock`, `Prompt`, collection
creation or a graphical password dialog. A returned prompt requirement is an
error, not permission to prompt. Sessions are explicitly closed on drop.

## Dependencies and validation

The native workspace registers this crate with the selected native versions: `zbus 5.19.0`
(MIT), `security-framework 3.7.0`/`security-framework-sys 2.17.0`
(MIT or Apache-2.0), `core-foundation 0.10.1` (MIT or Apache-2.0),
`windows-sys 0.61.2` (MIT or Apache-2.0) and `zeroize 1.9.0`
(MIT or Apache-2.0), plus existing URL/BLAKE3/error dependencies. Platform crates
are target-specific. There is no keyring executable, credential database format
or plaintext file storage implementation.

Five injected-backend fixture groups exercise scope pinning/canonicalization,
Session mode without OS access, remembered-provider recreation, exact deletion,
locked failures without fallback, redacted errors and rejected token/URL input.
Linux-only pure address/owner-policy fixtures verify that executable/network bus
transports cannot enter the backend and only unique bus peers can become the
captured destination. These tests do not connect to a user bus,
read/write a real vault or produce a prompt. This source increment has not run
builds or platform vault acceptance tests. Live Windows/macOS/Linux acceptance
with a temporary explicitly approved vault record remains a separate step;
compilation/injected fixtures do not establish that every desktop keyring is
available, unlocked or configured correctly.
