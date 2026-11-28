# Document design checks

Open a DNA document and choose a check from **Analyze**. These checks inspect the
current source and explicit declarations. Editing the document or changing the
declarations discards the result; run again before using its findings. Coordinate
buttons select the corresponding sequence interval, including split intervals at
a circular origin. None of these tools changes the sequence.

## Construct QC

Choose **Construct QC** and a default genetic code. A CDS `/transl_table` qualifier
takes precedence. Under **Declared annotation roles**, select features whose roles
and coordinates you explicitly assume for the check. Imported labels alone do not
establish their biological identity; features requiring review remain indeterminate.

Run **Run structural checks** to inspect CDS frame, start and termination, internal
stops and disagreement with a stored translation. Optional declarations add:

- Literal feature-to-reference comparison, retaining reference sequence, source URL,
  retrieval date and licence. Matching sequence establishes identity to the supplied
  reference, not the authenticity of its provenance or biological function.
- Fusion boundary register against an explicitly selected CDS and component.
- Current junction register at a declared boundary. This does not prove that an
  earlier assembly or edit caused a frame change.
- Structural element chains around a CDS. Declare upstream/downstream annotations,
  their roles and pinned references, and optional distance limits. Current role
  assumptions and exact reference identities are required before assessing strand,
  order and spacing. A distance limit is your requirement, not a biological activity
  threshold. Circular chains must fit within one turn without overlapping flanks.
- Required literal reference intervals with minimum/optional maximum counts and
  orientation. Capture references without selecting an annotation to use them for
  these searches. Counts include overlapping literal occurrences; they do not
  establish functional copy number. Palindromic intervals count once, with unresolved
  direction. Ambiguity, unsupported references and incomplete searches stay visible
  as indeterminate results. Search limits are shared across requirements.

Read each finding's verdict and explanation independently. Missing or unsupported
evidence yields indeterminate findings. There is no combined biological pass or
expression prediction. Initiation-context, unannotated-ORF, causal assembly-junction
and vector-class profiles remain planned. Remove dependent declarations before
deleting a pinned reference; removing one comparison does not delete its reference.

## Assembly standard format

Choose **Assembly standards…**, an installed profile and a declared slot.
Two deliberately scoped profiles are available:

| Profile | Checked scope |
| --- | --- |
| BioBrick RFC10 part-sequence format (2007) | Distinct coding/noncoding flanks, mandatory forbidden sites and advisory rules; no whole-vector, strain or registry certification |
| iGEM RFC1000 Level-0 documented body | Four declared role/fusion pairs, CDS `ATG` overlap, BsaI/SapI absence and adjacent role-subpath hints; no full transcription-unit or vector certification |

RFC1000 bodies exclude external flanks. Its CDS `ATG` overlaps the declared fusion
label; do not prepend an extra fusion label to the body. The profile does not guess
vector clamps or spacers.

Enter a one-based source start and part length. Circular intervals may cross the
origin once. Reverse orientation applies to the extracted part. Enter body start
and inclusive body end relative to that oriented part, or explicitly use the
flank-based body helper. The displayed required flanks are checked rather than
silently added to the sequence. Choose **Inspect declared part** and review each
finding, its mandatory/advisory status and any search limits.

This panel inspects one formatted source part. Multi-part slot/order assessment is
available through MCP; grammar hints are not chemical compatibility or successful
assembly predictions. Other toolkit profiles and assembly-wizard integration
remain planned.

## Design versus supplied consensus

In the comparison panel, select **Design versus supplied consensus** under
**Workflow**. Use **Use active document as reference** to pin the intended design,
then choose the supplied consensus explicitly. Source navigation keeps the design
pinned; names or equal DNA do not merge the two roles.

Choose **Whole construct** and **Global** explicitly in the shared comparison
settings. Strand/origin search and compute limits are the same controls used for
ordinary pairwise review. Stored topology is the default consensus-axis basis.
Declare a single circular turn only when that is your explicit interpretation of
the supplied sequence, and provide its rationale. It does not establish physical
circularity or permit an arbitrary repeated consensus to count as one turn.

For a conditional coding summary, enable design-CDS projection, load the exact
design annotations and select the intended CDS. The menu offers the first 100
CDS annotations; the annotation-number field can select another exact annotation.
Declare the genetic-code number and explain the role assumption. Results retain
the resolved code descriptor, ordinary codon assignments and exact code pin.
Unsupported or stale annotations, ambiguous calls and incomplete projection stay
explicit. No consensus annotation is manufactured or accepted as biological truth.

Choose **Compare design and supplied consensus**. The worker computes one
pairwise path, then reports coding register and codon consequences conditional
on that path. Review origin/strand completeness, paired input fractions and
selected-path proof separately: even a proven best score need not identify a
unique correspondence. Coordinate buttons use native original-source spans;
split codons can open a multi-interval source view. They do not infer locations
from matching text.

**Save event-bound comparison** retains the path, semantic roles, source events,
axis/CDS/code declarations, provenance and coding evidence in the existing
comparison record. Saved comparisons reopen as historical evidence, with the
same semantic details; opening does not recompute them. A historical semantic
record is not silently recomputed as an ordinary pairwise comparison. Choose
the design/consensus workflow and explicitly declare the intended inputs and
assumptions for a new computation.

Complete JSON retains all bounded evidence. HTML/PDF use the same retained report
and include its semantic assumptions, provenance, coding results and limitations.
Aligned FASTA/Clustal omit that semantic evidence and disclose the omission.
Changed sources, workspace or settings refuse stale responses, downloads and
navigation; run again. No source sequence, history or material is edited by the
comparison itself. Saving the comparison adds its evidence attachment.

A supplied consensus is a source observation, not experimental verification,
read-depth/quality evidence or proof of sample identity, expression or function.
The current source implementation is awaiting the shared compiler/execution
checks; no graphical acceptance is claimed.

## Oligo collection binding

Choose **Oligo collection…** and paste FASTA or tab-separated rows. FASTA
headers may end in `tail=N`; tab-separated rows contain name, ordered DNA and an
optional declared 5′ tail length. For example:

```text
>forward tail=4
GGGGACGTCAGTTCGA
>reverse
TCGAACTGACGT
```

Set the minimum annealing length, allowed extra unpaired prefix, mismatch policy
and required 3′ pairing. Advanced controls expose alignment scoring, collection-wide
work/output limits and sodium/strand-concentration conditions. Choose **Search
collection**. Declared tails are excluded from the footprint; additional unpaired
5′ bases are reported separately. Both strands and circular origins are searched.

Review each oligo's completion state before interpreting zero hits. A limited or
unsearched oligo does not establish absence. Budgets below the minimum needed to
bind inputs are refused with the required admission allowance. After admission,
an exhausted search reports the limit explicitly. Scores rank ungapped suffix
alignments; they do not predict PCR extension or amplification.

Each retained hit shows source coordinates, strand, aligned columns, mismatch and
3′ pairing bounds, and two distinct thermodynamic contexts. Annealed-stretch Tm is
available only for a concrete perfect duplex supported by the sodium-only model.
Mismatch and ambiguous-duplex Tm are unavailable. Whole-oligo Tm assumes a hypothetical
perfect complementary partner for the entire ordered oligo, including its tail;
it is not its Tm against this template or a recommended PCR temperature. Magnesium,
additives and modified bases are not modelled. Downloaded JSON retains the complete
bounded evidence; the screen paginates hits.

The corresponding [MCP tools](MCP.md) use immutable document references and the
same scientific owners. Saved report hashes bind content and declared context;
fresh recomputation establishes agreement with the current inputs.
