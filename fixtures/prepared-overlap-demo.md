# PCR + digest → ordered overlap example

These short synthetic sequences demonstrate coordinate and junction arithmetic,
not a laboratory-ready protocol. Open `prepared-overlap-demo.fa`; choose
`independent_target` as the active target, then **Assemble** in the inspector.
Select HiFi guidance and a **linear** product.

For fragment 1, choose `left_template` and **Exact PCR**:

- Forward annealing, 5′→3′: `TTGCAATCGT`
- Reverse annealing, 5′→3′: `GCTAGCATCG`
- Both 5′ extensions: empty

For fragment 2, choose `right_digest_source` and **Restriction digest**:

- Enzyme: `EcoRV`
- Fragment number: `2` in the UI (zero-based index `1` in the typed core)
- End repair: none; EcoRV produces blunt ends

Leave both orientations forward and enter one overlap length, `15`.
Compute method evidence and inspect the exact junction:
`ATCGTCGATGCTAGC`.

The complete right source is `GATATCGTCGATGCTAGCGGTAC`. Cutting `GAT|ATC`
leaves the selected right fragment `ATCGTCGATGCTAGCGGTAC`. Exact PCR yields
`TTGCAATCGTCGATGCTAGC`. The 15-base overlap appears once in the 25-base
independently supplied target `TTGCAATCGTCGATGCTAGCGGTAC`.

Review the candidate, independently verify it, and use the existing HTML/PDF
protocol exports. The retained plan contains the original template and complete
digest source, with PCR and digestion as explicit steps; cached fragments do not
replace original inputs. Switching to `wrong_target` must invalidate old evidence
and a fresh check must reject the changed final base. Reversing fragment 2 must
also fail the intended ordered junction.

Method guidance, complete cleavage and idealized preparation assumptions remain
visible. Sequence consistency does not establish experimental success, reaction
mixture uniqueness or empirical primer performance.
