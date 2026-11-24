# Methylation assessment

`core::methylation` evaluates sequence context under explicitly recorded
**source-generation assumptions**. It does not measure methylation, infer it
from a host name or imported metadata, propagate physical modifications, change
existing digest operations, or predict experimental success. No external
executable or network request is used at runtime.

## Versioned profile and factual provenance

The engine is `methylation-1`; the initial factual profile is
`neb-starter-2026-11-24-v1`. Reports retain the complete selected profiles,
formulation identity, source URLs, profile hash, exact source revision and
options. The following compact factual classification was independently
transcribed from [NEB's Dam/Dcm/CpG guide](https://www.neb.com/en/tools-and-resources/selection-charts/dam-dcm-and-cpg-methylation),
checked on 2026-11-24. It is a scoped supplier profile, not a universal enzyme
property or a complete supplier catalog.

| Starter identity | Dam | Dcm | CpG |
|---|---|---|---|
| BamHI, HindIII, PstI, BbsI, SapI | Insensitive | Insensitive | Insensitive |
| XbaI | Overlap blocking | Insensitive | Insensitive |
| NotI, SmaI | Insensitive | Insensitive | Blocking |
| EcoRI | Insensitive | Insensitive | Some overlaps block |
| EcoRV | Insensitive | Insensitive | Some overlaps impair |
| BsaI R0535 | Insensitive | Conflicting published descriptions | Some overlaps block |
| Generic BsmBI | Unsupported formulation | Unsupported formulation | Unsupported formulation |

[XbaI's manufacturer FAQ](https://www.neb.com/faqs/is-xbai-affected-by-methylation)
identifies preceding `GA` and following `TC` contexts; the implementation scans
actual overlapping `GATC` motifs on the original reference coordinates. Merely
classifying an enzyme as sensitive never blocks all of its sites.

The exact [generic BsaI R0535 product page](https://www.neb.com/en-gb/products/r0535-bsai)
is retained, rather than silently substituting BsaI-HFv2. Its Dcm sensitivity
section describes some-overlap impairment while a product note describes
blocking. This conflict is explicitly unsupported. The CpG some-overlap rule
and Dam insensitivity come from that same generic product page. The current
NEB table identifies **BsmBI-v2**, which does not establish a formulation-specific
profile for the generic BsmBI starter identity. Its three sensitivity rules
therefore remain unsupported. These unavailable/conflicting facts are useful
unknown results rather than borrowed formulation claims.

The [EcoRV product page](https://www.neb.com/en-us/products/r0195-ecorv)
also describes impairment for only some CpG overlap combinations. This profile
has no independently verified enumeration of the exact affected combinations
for EcoRI, EcoRV or BsaI. A relevant overlap under those rules returns `Unknown`,
not an assertion that every overlap blocks or impairs. `Impaired` is an explicit
outcome type for a supported impairment rule; the initial supported rules do
not invent a determinate impairment outcome from a some-combination label.

## Scientific scope and outcomes

The three independent source-generation categories are:

* Dam: N6-methyladenine in `GATC`.
* Dcm: 5-methylcytosine at the second C in `CCWGG`.
* CpG: 5-methylcytosine in `CG`.

Each category is `Unmodified`, `FullyModified`, `Unknown`, or `HemiOrMixed`.
Fully modified means both strands at every matching motif **on this exact
source** under the caller's declared assumption. Hemi/mixed states are not
modeled and return unknown; no binary complete-cleavage rule or invented
probability is applied. Other modifications and methylases remain outside this
profile even if they could influence an actual digest.

Per-category outcomes are `PermittedUnderProfile`, `Blocked`, `Impaired`, or
`Unknown`. Insensitivity or an explicitly unmodified category can permit that
category. For supported sensitive rules, absence of any possible overlapping
motif permits that category under the current-source model. Unknown states,
possible-only motif overlaps and unsupported/context-incomplete rules remain
unknown. A definite blocking effect dominates other category outcomes; otherwise
unknown dominates impairment. A possible restriction recognition site always
has aggregate outcome unknown, regardless of its category evidence.

Recognition certainty, methylation assessment and available cut substrate are
separate facts. Sites whose Type IIS cuts lie outside a linear source are
retained with absent cut coordinates, not discarded or converted into complete
digest products. Conditions, substrate accessibility, kinetics, star activity,
reaction mixtures and partial-digest products are not evaluated.

## Source coordinates and ambiguity

Both recognition orientations and circular-origin crossings are searched.
Palindromic recognition patterns are represented once. Each motif overlap
contains its exact ordered original-reference positions; circular coordinates
wrap modulo the source length. `certain` distinguishes a guaranteed motif
under IUPAC possibilities from a possible one. Pattern constraints are checked
base by base; ambiguous flanks never establish a certain blocking context.

`ModifiedBase.position` is an original top-reference base index.
`reverse=false` identifies the top strand, `reverse=true` the antiparallel
bottom strand; `base` is the physical modified-strand A or C, not necessarily
the top-reference character. Each motif records both potential modification
positions. These are inspected motif targets under the declared state, not
observed methyl groups, including when the category is unmodified or unknown.

Do not infer inherited modifications by rescanning an isolated product: a
modified base can survive cutting after its original methylase motif is broken.
Physical propagation is provided separately by the [tracked-strand model](MODIFIED_DUPLEX.md), using explicit source assumptions and actual strand mappings.
This assessment neither asserts nor restores such chemistry.

## Bounds and authority

Source preflight bounds sequence length (one million bases), serialized rich
source size (8 MiB), annotation/metadata structure and location depth before
validation and source hashing. Options limit work comparisons (default 20
million, hard 100 million) and retained sites (default/hard 4,096). Recognition
and context base comparisons share the work counter. A bounded output ledger
caps serialized site evidence within the 8 MiB report limit. Preparation errors
are errors; exhausted work, site or output bounds produce explicit incomplete
reports. Omitted evidence never establishes absence of a site or motif.

`assess` computes fresh evidence. `validate` checks bounded historical report
structure and cross-field consistency without scanning a source sequence.
`validate_bindings_for` additionally binds the current engine/profile, exact
source and options; it does not grant imported evidence scientific authority.
`replay` recomputes fresh evidence and requires exact report equality. Browser
and MCP interfaces must disclose planning-only scope and incomplete searches;
no assessment creates a molecule, edits source metadata or silently upgrades
an existing physical digest.
