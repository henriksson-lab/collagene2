# Remaining homology-assembly material workflow

The current ordinary overlap recipe ends in `Operation::ProfiledOverlap`
(`assembly_recipe.rs`, `compile`). The verifier checks oriented sequence,
homology and the declared profile, then produces a reference molecule. It does
not establish private tracked duplex state for that final step.

This is deliberate. `run_material` and `run_material_capsule::replay` require
fresh tracked state before a selected saved product can supply PCR or digestion.
The negative fixtures in `physical_pcr_ancestry.rs` retain this boundary. The
Gibson/HiFi/In-Fusion profile labels currently give advisory overlap guidance;
they do not establish reagent chemistry or physical product ancestry. SLIC
coverage also remains outstanding.

Consequently, the practical “assemble a plasmid, save it, and use that product in
another cloning plan” workflow remains incomplete for overlap assembly. Primer
proposal and independent target agreement do not close that gap.

The next implementation should add a distinct physical homology-assembly
operation, preserving historical `ProfiledOverlap` behavior. Its inputs must be
actual tracked, oriented substrates. An explicit versioned model must account
for resection, strand pairing, any newly synthesized DNA, external end chemistry,
modification propagation and the conditions for sealing each nick. Unsupported
sealing must retain nicks or withhold reusable material. Start with bounded,
ordered intact linear substrates and explicit linear/circular product intent.

The desired target may be compared only after material reconstruction. Do not
manufacture tracked state from assembled letters, apply `exact_pcr` to a desired
product, erase producer markers or reuse the denaturation-based OE operation as
an assumed vendor assembly mechanism. Once actual tracked products exist, the
existing exact `FromRun` event/origin checks can retain their ancestry.

This is a source audit and implementation requirement, not a new chemical model
or a claim of completed method support. Primary mechanism references and
independent strand-level fixtures are needed before registering any new model.
