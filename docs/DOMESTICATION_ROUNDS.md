# Sequential domestication preparations

`domestication_preparation::Route::InversePcrRounds` declares two to32 inverse-PCR
rounds for scattered substitutions in one original circular DNA source. The normal
source/CDS/restriction checker freezes the complete intended target before primer
construction. Each selected exact original-base substitution belongs to exactly
one declared round, and every round must carry a selected change. Round order is
part of the request identity; the software does not choose a successful order
from a reconstructed product.

Each round names an original-source half-open replacement interval, forward and
reverse annealing lengths, and a replacement-prefix length supplied by the reverse
primer tail. Origin-crossing intervals may unroll at most one original turn.
The web form uses original 1-based inclusive `start-end:F:R:reverse_tail` entries,
one per line. Raw values survive route switching. Source and independent intended
reference review precede the physical preparation action.

The single canonical plan contains one original input and three operations per
round: exact PCR, explicit phosphorylation, and intended closed blunt circular
ligation. Each later PCR references the preceding ligation step. No recoded
reference or intermediate intended genotype becomes an additional plan input,
saved material capsule or template stock.

Primer declarations use an explicitly intended intermediate sequence frame,
represented only as `Sequence`. They incorporate previously declared substitutions
and retain actual inverse-PCR origin shifts. For an original interval starting at
`s` and reverse-tail prefix length `r`, the intended next PCR reference origin is
`(s + r) mod n`: forward5′ tail suffix, retained template, reverse-tail complement.
The next round translates original coordinates through that origin. The molecular
verifier independently reconstructs every amplification and ligation from actual
preceding material; ambiguity or an unsuccessful preceding step cannot be repaired
from the intended sequence frame.

Admission bounds apply before frame copies or primer extraction: at most32 rounds,
8–80-base nonoverlapping annealing domains, complete oligos at most200 bases, one
proper replaced interval per round, aggregate100M declared comparison reservation
and64MiB material reservation. Existing native operation, source, alternatives and
report limits remain separate. These are model bounds, not a laboratory efficiency
or absolute end-to-end runtime claim.

Readonly proposal validation regenerates the exact declared oligos, source-frame
origins and canonical step references without ranking primers or replaying PCR.
Fresh scientific verification remains required at use/export. Normal candidate,
run history, oligo order tables and HTML/PDF paths consume the same canonical plan,
including all rounds and retained alternative-product evidence. No template
clearance, purification, transformation or experimental success is inferred.

Independent native fixtures cover two distant BsaI changes, a real preceding-step
versus original-stock negative, reordered requests, missing/duplicate patch
coverage, partial reverse-tail splits and an original origin-crossing second round.
The expected final121-base genotype is a literal oracle, not a simulated product.
