# collagene-assembly-standards

Independent MIT Rust grammar and part-format checker. Dependencies are limited
to the shared molecule model, enzyme catalog, serde and BLAKE3; there is no core,
browser, network, filesystem or external scientific-tool dependency.

This crate is registered with source-bound core, worker, web and MCP paths. Its
seven integration cases and connected native/WASM checks pass. Public DTOs
have an optional `json-schema` feature and externally represented enums suitable
for Postcard; the authored tests include binary roundtrips.

The first supported built-in is **BioBrick RFC10 part-sequence format (2007)**,
sections 1–3 of Tom Knight's 3 May 2007 draft. It checks the exact non-coding and
coding prefixes, suffix, mandatory coding-start pattern, five required internal
site absences including NotI, and separately marked recommendations. It does
not claim full vector, host, sequencing-primer or registry compliance.
[Profile data and source/rights notice](data/rfc10-2007/NOTICE.md).

```rust
use collagene_assembly_standards::profile::{Profile, RFC10_ID};
use collagene_enzyme_catalog::Catalog;
let catalog = Catalog::builtin()?;
let profile = Profile::builtin(RFC10_ID, &catalog)?;
let hints = profile.legal_successors(&[])?;
assert!(!hints.may_finish);
assert_eq!(hints.allowed_slots.len(), 2);
# Ok::<(), Box<dyn std::error::Error>>(())
```

`Profile` owns its validated data definition and exact catalog edition/content
pin. Each enzyme ID must have the same recognition and strand-cut geometry as
the pinned catalog. Sources, profile revision, slot names, order, cardinality,
declared overhangs, additions and extra avoided sites are data. Custom grammars
are explicitly `Custom`; changing built-in facts while retaining its reserved ID
is rejected. A user-provided source URI, asserted name or hash is not authority.

`Part` supplies an oriented formatted DNA string and an explicit half-open body
interval. Expected sites in the format additions are excluded from the internal
body scan. Forward/reverse IUPAC comparisons distinguish certain violations from
possible sites; advisory motifs do not invalidate conformance. Resource limits
bound source count/bases, comparisons, findings and serialized output. Partial
checks cannot establish conformance. The immutable assessment binds profile,
catalog, source sequence hashes, order, body intervals and limits; scientific use
of a retained producer report requires `Assessment::validate_for` fresh replay.

Legal-successor hints establish one-step declared grammar consistency only.
RFC10 specifies no promoter/CDS/terminator functional order, so this profile
allows its two declared part classes in either order and permits repetition.
It does not predict ligation, biochemical activity, reading-frame preservation,
function or expression. The existing core remains responsible for actual
physical assembly, intended-target checking and protocol execution.

Default bounds are 64 parts, one million aggregate DNA bases, eight million
residue comparisons, 1,024 retained findings and 512 KiB report JSON. Hard bounds
are 128 parts, two million bases, 64 million comparisons, 4,096 findings and
one MiB report. Profile JSON is limited to 128 KiB before decoding. Typed native
definitions receive the same bounded preflight before hashing/installation.

MoClo, YTK, EMMA, PhytoBricks, SEVA and RFC1000 profiles are not installed by this
crate. They require their own exact source, data-license and complete grammar
curation; lookup refuses unsupported built-in IDs rather than guessing tables.
