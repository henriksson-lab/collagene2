# RFC1000 Level-0 profile admission — 2026-11-28

Admitted a second built-in profile in the existing `assembly-standards::Profile`
owner: `igem-rfc1000-level0-documented-body-2026`. No separate grammar
representation, standard-name alias or compatibility reader was introduced.

Primary authority is [iGEM Foundation's DNA Assembly technical specification](https://technology.igem.org/technologies/assembly/),
Assembly Compatibility and Basic Parts (Level 0s), Fusion Sites. Four oriented
role pairs are captured in [the source/rights notice](../../crates/assembly-standards/data/rfc1000-l0-2026/NOTICE.md).
The current live source was reviewed through indexed primary-page retrieval;
direct opening returned a JavaScript shell. No immutable upstream revision or
raw-byte digest is invented. The exact curated facts, citations and policy are
covered by the installed profile's existing BLAKE3 identity and catalog pin.

Documented bodies have empty external additions. Declared CDS `ATG` overlaps the
last three positions of its `AATG` fusion label; this profile checks the body
prefix and does not concatenate that fusion label to the body. Right fusion
labels stay in assembled-reference 5′→3′ orientation. Required internal absence
rules are BsaI/SapI only, scanned on both strands. Unsupported roles, levels,
vector flanks and suggested extra restriction sites are refused or omitted.

Our admission grammar accepts contiguous nonempty subpaths of the four supplied
roles. This supports individual-part linting and positive adjacency hints without
mislabeling a partial chain a full transcription unit. All roles may begin/end a
subpath; only the three adjacent transitions are legal. No biological-role
classifier, protein validation, cleavage simulator or experimental-success claim
is added. Existing assessment limits and replay remain authoritative.

Only independently authored factual curation and implementation are MIT. No
license for upstream prose/images was verified, and none is relabeled or
redistributed. The notice distinguishes those rights and retains exact source
attribution. No cloning application implementation was consulted.

Independent source-only fixtures use literal four-pair and restriction-motif
expectations, body-only intervals, ATG overlap, correct/wrong role ordering,
ambiguous sites, extra-site nonrequirements, reserved-ID tampering and
cross-profile replay rejection. Parent owns combined compilation/native/WASM
checks; this worker ran none and made no core/web/MCP/schema changes.


Parent integration now exposes both profiles through the existing worker, web
selector and MCP catalogue. Final native/WASM checks pass, as do all six RFC1000
fixtures and seven RFC10 regressions. [Integration evidence](reaction-elements-current-integration.json).
