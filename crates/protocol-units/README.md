# collagene-protocol-units

MIT Rust primitives for exact cloning-protocol arithmetic. No sequence model, inventory, reagent-stock database, scheduling, automation, network or external-program runtime is included. Stock concentrations here are explicitly supplied numerical inputs, not assertions about stored reagents.

The default library has **zero dependencies**. Optional `serde` adds typed serialization. This library is a member of the Collagene workspace; its arithmetic does not depend on the scientific molecule model or UI.

```rust
use collagene_protocol_units::*;

let recipe = DilutionRecipe {
    schema_version: SCHEMA_VERSION,
    stock: Some(MolarConcentration::from_unit(
        Rational::integer(10), MolarUnit::Micromolar)?),
    target: Some(MolarConcentration::from_unit(
        Rational::integer(100), MolarUnit::Nanomolar)?),
    final_volume: Some(Volume::from_unit(
        Rational::integer(20), VolumeUnit::Microliter)?),
};
let plan = dilute(&recipe)?;
assert_eq!(plan.stock_volume.in_unit(VolumeUnit::Microliter)?,
           Rational::new(1, 5)?); // exactly 0.2 µL
assert_eq!(plan.diluent_volume.in_unit(VolumeUnit::Microliter)?,
           Rational::new(99, 5)?); // exactly 19.8 µL
# Ok::<(), Error>(())
```

## Quantities and bounds

`Rational` stores a reduced nonnegative `u64/u64` fraction with nonzero denominator. Construction, decimal parsing, addition, subtraction, multiplication and division are checked. Comparisons use exact `u128` cross-products. Compound calculations cross-cancel the **whole** product/quotient before range checks; at most 16 factors are supported. Results outside the reduced `u64/u64` representation return `Overflow`. Mixture sums use reduced checked `u128` intermediates, so an intermediate prefix may exceed the public numerator range and still produce a valid final result. An intermediate denominator/product beyond `u128` is conservatively rejected even if later terms could permit reduction; arbitrary-precision arithmetic is not claimed. Nothing saturates, wraps or silently becomes floating point. The API deliberately supplies checked methods instead of panic-prone arithmetic operators.

Decimal parsing accepts ordinary nonnegative finite decimal notation, up to 64 input bytes and 19 fractional digits. Exponent notation, NaN and infinity are rejected. `Rational::new` supports exact fractions beyond this decimal input grammar.

The distinct types are `Volume` (canonical nL), `Mass` (ng), `Amount` (pmol), `MassConcentration` (ng/nL), and `MolarConcentration` (pmol/nL). Unit enums support nL/µL/mL, ng/µg, pmol/nmol and nM/µM. ASCII `u` and both Unicode micro spellings are accepted by the unit parsers. Units remain case-sensitive: `nm` is not `nM`. Concentration dimensions cannot be mixed in a dilution recipe.

## Dilutions, mixtures and rounding

- `dilute` implements C1×V1 = C2×V2 for a single stock and a zero-concentration diluent, under additive-volume assumptions. It returns exact stock and diluent volumes. Missing fields, zero stock concentration, zero final volume and a target above the stock have distinct errors. A zero target with a positive stock produces an all-diluent recipe.
- `balance_mixture` sums at most 128 uniquely named components with nonblank role identifiers and returns the remaining fill volume. Missing component volumes are reported with their index. Overfilled mixtures fail; they never produce a negative fill volume. Role identifiers are checked without normalizing their identity: nonblank surrounding whitespace is retained, and duplicates are compared exactly. The fill diluent is additional to the supplied components, so do not supply the same fill twice.
- `round_volume` requires an explicit increment and mode: `ExactOnly`, `Floor`, `Ceiling` or `NearestTiesUp`. It returns the requested and rounded volumes, direction, and exact absolute deviation. The increment is a declared resolution, not a claim about an instrument's accuracy.
- `quantized_dilution` rounds the stock volume, refills to the original final volume and reports the **achieved** concentration. It does not preserve a requested target label after rounding changes that concentration. Rounding that would exceed the final volume fails.

These are arithmetic recipes, not measured concentrations or recommendations about pipettes, solvents, buffers, salt, ligation chemistry, volume contraction or reaction success. Inputs and outputs stay exact under the declared arithmetic assumptions; physical measurements and the DNA mass approximation do not acquire extra certainty.

## Explicit duplex DNA mass model

There is **no default DNA mass coefficient**. Supply a `DuplexDnaMassModel` with a positive `grams_per_mole_per_base_pair`, a nonblank model identity, nonblank assumptions, and optional source attribution. Model identifiers are limited to 128 bytes; assumptions/source attribution to 2048 bytes each. Length must be positive. This model counts duplex DNA **molecules**, not molar ends; it does not model ssDNA, sequence composition, modifications, counterions or asymmetric overhangs.

For a declared coefficient `k` in g/mol/bp:

```
amount_pmol = mass_ng × 1000 / (length_bp × k)
mass_ng     = amount_pmol × length_bp × k / 1000
```

For example, **if the caller explicitly declares 660**, 50 ng of 1000 bp duplex DNA gives exactly `5/66` pmol **under that approximation**. A declared 650 gives `1/13` pmol instead. The crate intentionally retains the chosen coefficient and assumptions rather than silently selecting one.

`plan_dna_ratios` takes one weighed anchor and up to 127 target components. For one shared mass model, target mass is anchor mass × requested molecule ratio × target length / anchor length. Thus 50 ng of a 3000 bp vector and a declared 3:1 molecule ratio require 50 ng of a 1000 bp insert, not 150 ng. Optional supplied mass concentrations provide required volumes. Set `require_stock_volumes=true` to reject missing concentrations; otherwise the corresponding volume is explicitly `None`. The ratio is a caller choice; no optimal assembly ratio or efficiency is inferred.

Primary references consulted for the dimensional model only:

- [NEB nucleic acid data](https://www.neb.com/tools-and-resources/usage-guidelines/nucleic-acid-data) presents length-based duplex molar-mass estimates and rounded mass/pmol tables. Its page contains differing numerical conventions; this crate does not treat them as one exact universal coefficient.
- [NEB E3340 manual](https://www.neb.com/en/-/media/nebus/files/manuals/manuale3340.pdf) explicitly illustrates a 660 g/mol × base-pair-length concentration conversion.

No other cloning application's source code was consulted or copied.

## Recipes, serialization and replay

Recipes/results carry schema version 1 and result engine `protocol-units-1`. Results retain the exact recipe. `validate()` recomputes the corresponding calculation and rejects altered results or unsupported versions; serialized output is not an independent authority. Deserializing a recipe is followed by its calculation, and deserializing a result by `validate()`.

With `features=["serde"]`, rationals serialize as **decimal strings**, e.g. `{"numerator":"1","denominator":"5"}`, avoiding JavaScript number precision loss. DNA base-pair lengths also serialize as exact decimal strings. Quantities carry a strict canonical unit tag, e.g. `{"unit":"nL","value":{"numerator":"200","denominator":"1"}}`. A volume payload cannot deserialize into mass. Unknown recipe fields, invalid fractions and wrong canonical unit tags are rejected. Consumers must bound incoming serialized bytes before deserializing; Serde does not supply a transport or streaming byte budget. Collection/text limits are checked before calculation or cloning input recipes.

```sh
cargo test -p collagene-protocol-units --locked
CARGO_TARGET_DIR=/tmp/collagene2-next-targets/units cargo test --offline --locked --features serde -j2
```

Tests use literal independent arithmetic fixtures, exhaustive small rational/concentration/rounding grids, representation boundaries, dimensional serialization guards, and altered-result rejection. They require no external scientific software or network access.
