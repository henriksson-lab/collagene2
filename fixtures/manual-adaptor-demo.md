# Synthetic adaptor-into-backbone example

Open [manual-adaptor-demo.gb](manual-adaptor-demo.gb). Select **Independent_target**
as the active document, then choose Cloning plan → **New manual plan**.

Add three inputs with these explicit IDs and source document choices:

| Input ID | Source document |
| --- | --- |
| `adaptor-top` | Adaptor_top |
| `adaptor-bottom` | Adaptor_bottom |
| `original-backbone` | Original_backbone |

Append these declarations in order:

1. **Oligo annealing**: top `adaptor-top`, bottom `adaptor-bottom`, minimum paired
   bases **6**, candidate index **0**. Explicitly declare **5′ phosphate / 3′OH**
   for both oligos. This is an input hypothesis, not a measurement of a reagent.
2. **Restriction digest**: template `original-backbone`, enzyme **EcoRI**,
   zero-based fragment index **0**.
3. **Ordered ligation**: first **Step 2 product**, then **Step 1 product**, both
   forward. Select circular product and leave retained-nick allowance disabled.

Use **Insert step declaration** for each edited operation, then **Apply edited
candidate**. Verify the candidate, inspect the actual intermediate products and
junctions, then use the existing named-output, history and HTML/PDF actions.

The independently frozen circle is `AATTCTTACGGGAATTGCTACG` (22 bp). The annealed
insert has six paired bases and two AATT cohesive ends; the circular backbone is
cut at its actual EcoRI site. These strings describe a small geometry fixture,
not functional DNA, a kit protocol, reagent quantities or laboratory success.
Native tests independently exercise the full preparation/ligation/history path.
Application browser acceptance is deferred.
