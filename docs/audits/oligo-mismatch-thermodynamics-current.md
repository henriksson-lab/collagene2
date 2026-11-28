# F38 isolated internal DNA mismatch thermodynamics

Status, 2026-11-28: parameter owner and source fixtures authored in
`crates/oligo-binding/src/thermodynamics/`. Parent integration into the existing
core thermodynamics calculation is a separate part of this coherent batch.
No compiler, test harness, browser or application HTTP check was run by this
worker. This is a scoped equilibrium model, not PCR efficiency or calibration.

## Data admission and primary sources

The original nearest-neighbor measurements support sequence-context-specific
internal mismatches rather than a generic temperature penalty per mismatch:

| Mismatch family | Primary source |
| --- | --- |
| G.T | [Allawi and SantaLucia 1997](https://pubs.acs.org/doi/10.1021/bi962590c) |
| A.C, pH dependence | [Allawi and SantaLucia 1998](https://pubs.acs.org/doi/10.1021/bi9803729) |
| G.A | [Allawi and SantaLucia 1998](https://pubs.acs.org/doi/10.1021/bi9724873) |
| C.T | [Allawi and SantaLucia 1998](https://doi.org/10.1093/nar/26.11.2694) |
| A.A, C.C, G.G, T.T | [Peyret et al. 1999](https://pubs.acs.org/doi/10.1021/bi9825091) |

The numerical capture is the general-purpose Biopython maintainer's
[pinned 1.86 table](https://raw.githubusercontent.com/biopython/biopython/biopython-186/Bio/SeqUtils/MeltingTemp.py),
whose individual header explicitly offers Biopython License Agreement or BSD
3-Clause. The permissive Biopython alternative was selected, with original
copyright, permission and disclaimer retained under
`crates/oligo-binding/data/mismatch-nn/`. Our MIT implementation is independent;
no cloning application or upstream numerical algorithm was copied. The source
SHA256 is `4f3cbc585ffd68d54386d40f5edb113a7c2eae4aa4388c6b496bceecdec06500`;
the [original license file](https://raw.githubusercontent.com/biopython/biopython/biopython-186/LICENSE.rst)
SHA256 is `d4bbb0c4e5345e6f2508c3b93fb247fcb2db95783fb8a584c23cd88e3c03df6c`.
The full capture is 43,151 bytes; only 48 numeric facts are incorporated.

The author-hosted [Moore et al. supplement, Table 2](https://www.maranasgroup.com/pub/2001/sup_moore-etal01.pdf)
also attributes the compilation to SantaLucia and supplies independent literal
parameter context. PDF text extraction has font corruption; this worker did not
claim visual inspection of the PDF. No publisher prose/PDF/image is incorporated
or assigned a software license. The table's inosine and three tandem-GT records
are excluded, leaving all 48 physical single-mismatch step classes.

## Conventions and supported footprint

`DuplexContext::new(source, topology)` checks bounded source admission and hashes
the source once. Its `assess(oligo, options, hit, assumptions)` independently
validates every retained hit column against the original source before using
the at-most-200-base footprint. Top sequence runs oligo 5′→3′; bottom is the actual
antiparallel partner 3′→5′. Forward hits use the complement of the reference;
reverse hits use the reference base. Original coordinates and circular wrapping
are retained. Declared tails and additional unpaired prefixes remain outside
the calculation.

For actual stack `XY/ZW`, physical 180-degree rotation gives `WZ/YX`.
This is neither top-only reverse complement nor top-only reversal. Each empirical
step is represented exactly in tenths of kcal/mol and cal/(mol K). For one
mismatch, its two flanking steps replace the corresponding canonical steps
exactly once. They are absolute ΔH/ΔS parameters, not penalties added on top of
canonical values. No second canonical table or Tm solver is added to this crate.

The admitted mismatch profile requires exactly one internal A/C/G/T mismatch,
at least two canonical pairs at each end, and distinct physical strands. The
two-pair margin is an explicit conservative software scope policy, not a newly
measured universal boundary. Identical physical strands with a central mismatch
are explicitly `SelfAssociation`, since the distinct-strand concentration rule
does not apply. A self-complementary top sequence can still have a distinct
mismatched partner; the caller must remove its hypothetical perfect-baseline
symmetry correction before solving that heteroduplex. Fully canonical duplexes
return `PerfectDuplex` and keep the existing canonical calculation unchanged.

`Assumptions` explicitly records pH7, sodium-only entropy correction and equal
concentrations of two distinct strands. It validates sodium 0.01–1 M and total
strand concentration 1 pM–10 mM, matching the current core owner; these are
implementation guardrails, not evidence of calibration throughout that range.
The underlying studies used controlled solution conditions. A.C protonation is
pH sensitive; this profile does not infer the actual pH of a PCR buffer.

Parent integration uses the existing
[SantaLucia 1998 canonical model](https://doi.org/10.1073/pnas.95.4.1460):
replace the two canonical stack terms; keep canonical end initiation; remove any
hypothetical top-strand self-complementarity entropy term; apply
`ΔS_Na = ΔS_1M + 0.368*(n-1)*ln[Na+]`; then solve
`Tm(K) = 1000*ΔH / (ΔS_Na + R*ln(C_total/4))` for distinct equal-concentration
strands. Units are mol/L and R is cal/(mol K). This is not a primer-excess PCR
mass-action model. The 2004 review's updated canonical parameter set must not
silently replace the current 1998 baseline.

## Explicit refusal and evidence boundaries

Typed unsupported outcomes cover IUPAC ambiguity, short duplexes, multiple
mismatches, near-terminal mismatch, self-association and missing parameters.
Invalid inputs/forged coordinate bindings remain errors. There is no silent
fallback to perfect-duplex Tm. Gap/bulge, dangling/terminal mismatch, modified
bases/inosine, competing hairpins/dimers, Mg, dNTPs, additives, RNA hybrids and
other pH models require their own admitted parameter sets and conventions.
In particular [Owczarzy et al. 2008](https://doi.org/10.1021/bi702363u) treats
mixed Mg/monovalent effects separately; sodium equivalence is not inserted here.
[SantaLucia and Hicks 2004](https://doi.org/10.1146/annurev.biophys.32.110601.141800)
distinguishes structural motifs requiring additional datasets. None of this
predicts extension, kinetics, PCR yield, biological success or specificity.

`Profile::current()` retains the source URL and exact source/license hashes,
parameter count, model scope and BLAKE3 profile identity. Each assessment binds
that profile, source/topology, oligo identity and ordered sequence/tail, options,
orientation/coordinates and concentration/salt/pH assumptions. Typed facts are
Postcard-compatible with bounded strings/vectors and optional JSON Schema.
`Assessment::validate_for` recomputes the bounded replacement evidence rather
than trusting a retained claimed result. Source lifecycle/current worker checks
and the final numerical Tm remain the core integration owner's responsibility.

## Independent source fixtures

Authored tests include literal AG/TT and TT/AG terms and a two-stack arithmetic
fixture, all 96 oriented single-mismatch steps forming 48 rotational classes,
family representatives and explicit tandem exclusion. They construct forward,
reverse, circular and tail-bearing hit columns without calling the production
search recurrence, verify source/conditions/tamper bindings, and exercise every
admission refusal including identical-strand central mismatch and a
self-complementary top with a distinct mismatched partner. These source fixtures
are meaningful arithmetic/geometry checks, not experimental accuracy evidence.
Parent will run focused native/WASM integration checks in the shared batch.


Parent integration checkpoint: the shared native/WASM batch and 349 unique focused
cases pass. See [source and execution evidence](cloning-extensions-current-integration.json)
for the exact scope, initial catalogue-pin corrections and final metadata-only
wording edit. No graphical or application HTTP acceptance is claimed.
