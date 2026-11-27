# Published numbering windows

`antibody-numbering-windows-1` implements conditional source-order window rules
for `Kabat1991` and `Chothia1997`. These are native scientific components, not
a standalone domain classifier or CDR definition. The current
[domain-track assembly](ANTIBODY_NUMBERING.md#alternate-frameworks-and-independent-region-definitions)
combines them with invariant reference frameworks and offers explicit application
numbering choices, retaining unsupported mappings as unresolved. No third-party code, trained models or numbering database is
bundled. The implementation encodes published numerical rule facts independently.

## Scope and sources

Three windows use independently supplied counterpart IMGT boundary assumptions:

| Window | Provided counterpart boundaries | Meaning |
| --- | --- | --- |
| `CysToTrp` | 23 / 41 | First loop plus three preceding and two following framework residues |
| `MiddleLoop` | 55 / 66 | The source residues between these framework boundaries |
| `CysToJ` | 104 / 118 | Junction interior, including residues outside the Kabat heavy CDR3 definition |

The first window's supported length is the stated IMGT-loop length plus five
flank residues. Other window lengths count only the actual interior residues.
Anchors are excluded from the numbered residue list and retain their separate
published boundary labels. These windows must not be installed as CDR annotations.

| Convention / locus | First-loop lengths | Middle lengths | Junction lengths |
| --- | --- | --- | --- |
| Kabat 1991 / heavy | 8–12 | 6–10 | 9–21 |
| Kabat 1991 / kappa | 5–12 | 3 | 7–15 |
| Kabat 1991 / lambda | 6–9 | 3 | 9–15 |
| Chothia 1997 / heavy | 8–10 | 6–10 | 11, 13, 19, 21 |
| Chothia 1997 / kappa | 5, 6, 7, 10, 11, 12 | 3 | 9 |
| Chothia 1997 / lambda | 7–9 | 3 | 9, 11 |

This finite coverage follows the explicitly inspected correspondence examples,
including chart rows added for completeness. It is not a biological length limit
or evidence that other loops cannot be numbered. Unsupported inputs return an
error, pending separately verified rules. Insertion letters use the module's
canonical uppercase spelling; historical diagrams also use lowercase for Chothia.

[Heavy-chain correspondence chart](https://www.imgt.org/IMGTScientificChart/Numbering/IMGTnumberingCDR_VH.html)
provides the inspected first/middle windows, Kabat junction cases and common
historical Chothia junction numbering. Heavy junction numbering here retains the
shared labels; historical CDR delimitations are not inferred from those labels.
[First-loop comparison](https://www.imgt.org/IMGTScientificChart/Numbering/CDR1-IMGTgaps.html)
and [junction comparison](https://www.imgt.org/IMGTScientificChart/Numbering/CDR3-IMGTgaps.html)
distinguish numbering from the separately bracketed CDR definitions.

[Kappa chart](https://www.imgt.org/IMGTScientificChart/Numbering/IMGTnumberingCDR_VK.html)
and [lambda chart](https://www.imgt.org/IMGTScientificChart/Numbering/IMGTnumberingCDR_VL.html)
specify locus-dependent first-window patterns and middle/junction examples.
The 1997 lambda junction convention uses the later insertion placement, distinct
from the earlier Chothia illustration. The kappa first-window example of length
11 and the lambda example of length 7 require their literal distinct patterns;
a generic sequential insertion counter would change those published labels.

Original references are Kabat et al., *Sequences of Proteins of Immunological
Interest*, NIH 91-3242 (1991); Chothia and Lesk, *J. Mol. Biol.* 196:901–917 (1987);
and [Al-Lazikani et al., *J. Mol. Biol.* 273:927–948 (1997)](https://pubmed.ncbi.nlm.nih.gov/9367782/).
Charts were inspected on 28 November 2026. Published rule facts are used here;
the original diagrams, page text and antibody databases are not redistributed.

## Binding and verification

`layout` returns exact integer/letter labels in source order for a declared
convention, locus, window and length. It allocates at most 21 interior labels.
`number_anchored_window` additionally validates an immutable explicit protein
source, exact source revision, provided boundary indices and bounded nonblank
provenance. It checks stop-free interiors/boundaries, literal Cys/Trp first-window
anchors, and locus-specific Cys/J plus G-X-G context for junction windows. Middle
anchors have no fixed amino-acid identity. These checks cannot establish that the
provided positions actually have the stated biological meaning or flank lengths.

The report binds source, assumptions, convention/locus/window, versioned rule
catalog and results. `validate_for` reconstructs from independently supplied source
and parameters. Source or annotation changes, altered labels/anchors, stops,
invalid coordinates and changed provenance refuse. No label-only offset mapping,
region inference or automatic selection of a convention occurs.

Five native fixture groups and optional-schema WASM checks pass. Literal published
examples cover distinct insertion positions, short-loop gaps and locus/revision
exceptions. A source-bound sweep exercises all 79 supported configurations for
length, label uniqueness, anchor exclusion, source ownership and reconstruction.
That sweep does not independently certify every biological assignment. Full-domain
assembly, region-definition rules, worker/UI/MCP selection and broader independent
residue-level validation remain required.
