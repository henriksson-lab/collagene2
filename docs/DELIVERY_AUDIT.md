# Current delivery gap audit

This is a source-based triage of current delivery gaps, not a completion
certificate. The full feature matrix and M0–M6 gates remain in
[the current plan](PROJECT_PLAN.md). Browser/application HTTP checks remain deferred
until in-scope implementation is complete.

| Area | Current source evidence | Remaining delivery work |
| --- | --- | --- |
| Protected codon ranges | `codon_design::ConstrainedPolicy` includes protected ranges and allowed codons; its independent checker reports violations | Current platform/workflow acceptance; the older blanket “protected ranges missing” entry is obsolete |
| Plate layouts | `protocol_layout` plus web/MCP capture and CSV/HTML/PDF export | Broader workflow acceptance; this is documentary cloning layout, not inventory or instruments |
| Lab comments/review | Native collaboration store, individual member registry, discussion SDK and shared browser/Tauri remote controls | Exact-source CI, graphical acceptance and deliberate deployment cutover; no automatic offline comment queue |
| F33 collection similarity | `sequence_similarity` plus current binary worker117, event-bound Compare controls and immutable-reference MCP discovery rank explicitly selected DNA/RNA sources | Current-source graphical/platform acceptance; server-wide/indexed similarity is outside this implemented profile |
| F34 specialized biology | Protein properties model an unmodified free chain; structure import/viewing has worker, workspace and MCP paths; `collagene-rna-chemistry` adds source-bound base/sugar/backbone/end declarations, stoichiometry, explicit imported-annotation review, event-bound saved declarations and Biology controls | RNA graphical/platform acceptance, broader independently validated antibody numbering/region coverage, graphical/performance acceptance remain. [Structures](STRUCTURE_IMPORT.md) · [RNA chemistry](RNA_MODIFICATIONS.md) |
| M4 source interoperability | URL adapters and retained-source parsers exist | Real provider/interoperability evidence remains separate from mocked fixtures |
| M6 installed desktop | Native CI packaging exists for Windows/macOS/Linux | Actual installation, offline graphical use and chosen-remote workflow acceptance on the supported platforms |
| Viewer/performance/accessibility | Existing GPU/layout sources and historical measurements exist | Current-source browser/hardware evidence and the user's visual review; historical screenshots/tests do not prove current acceptance |

This audit found no reason to reopen deferred BAM/BAI work or external Primer3 /
multiple-alignment engine ports. OpenGel, inventory/ELN, instruments, scheduling,
organizational reporting and mobile expansion retain their exclusions/deferrals.

## Antibody analysis acceptance scope

Benchling's current [antibody annotation reference](https://help.benchling.com/hc/en-us/articles/9684246077325-Annotate-antibody-specific-information-on-protein-sequences)
distinguishes variable-domain detection, alternative residue numbering, CDR/framework
region definitions and candidate liability tracks on individual proteins and
alignments. Numbering schemes include Kabat, Chothia and IMGT; CDR definitions
also include North. These are separate choices. Basic unmodified-chain protein
properties do not fulfill this feature requirement.

Independent native implementation still needs validated domain/numbering references,
explicit supported chain/species profiles and independent numbering fixtures.
Motif-only guesses must not be presented as validated domain numbering. Candidate
PTM patterns are sequence liabilities, not measured chemical state. The source discrepancy is resolved for the new screening core using the canonical
dipeptide definitions in Lu et al.; the exact profile and limitations are documented
in [protein liabilities](PROTEIN_LIABILITIES.md). Worker, Biology and native MCP integration are implemented. Conditional human
V/J analysis now derives IMGT/Kabat1991/Chothia1997 tracks with independently
selected IMGT/Kabat1991/canonical-chart Chothia regions on proteins and exact imported alignment rows.
Unsupported rule profiles, ambiguous boundaries and incomplete regions remain
explicit. Natural residue-level validation, broader chain/species calibration,
North region definitions and full application acceptance remain open.
Existing external-tool, MIT-implementation and provenance requirements remain unchanged.


The [scientific contract](ANTIBODY_NUMBERING.md) documents the single current
source-bound report, all-optimal invariant correspondences, finite published
window profiles, separate numbering/region membership and complete versus
unresolved CDR intervals. Current implementation evidence does not establish
universal automatic domain classification or independent biological numbering
accuracy. The fixed external panel establishes deposited-chain identity only;
it is not residue-position ground truth.
