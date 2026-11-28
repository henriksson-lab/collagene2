# IVT-template factual curation

Independent implementation and curation: MIT, as for this crate. No supplier
prose, images, parameter tables, product inserts or scientific software are
redistributed. Small sequence/coordinate facts are independently encoded in
`src/ivt/profiles.rs`; that source defines their exact installed identity.

Reviewed 2026-11-28:

| Source | Admitted fact and limit |
| --- | --- |
| [NEB T7/CleanCap AG FAQ](https://www.neb.com/en/faqs/do-i-need-to-change-my-promoter-sequence-to-use-cleancap-reagent-ag) | Consensus `TAATACGACTCACTATA` precedes +1; original AG capping requires +1/+2 `AG`. An altered promoter is a conflict with this selected profile, not proof of no transcription. |
| [TriLink original CleanCap AG N-7113](https://shop.trilinkbiotech.com/shopping/products/ivt-and-capping-reagents/rna-cap-analogs/cleancap-reagent-ag/) | Original N-7113 requires `AG` at +1/+2. Supplier page version 2026.11.25/build5011 was reviewed. No rule is transferred to AG3′OMe, M6, AU or another analogue. |
| [NEB plasmid-template FAQ](https://www.neb.com/en-gb/faqs/how-to-prepare-plasmid-templates) | Downstream linearization and preference for blunt/5′ overhangs. A 3′ overhang produces a documentary advisory; no yield or cleavage calculation is derived. |
| [Conrad et al. 2020](https://doi.org/10.1038/s42003-020-01167-x), [primary full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC7429497/) | Study-specific initiation context and end-heterogeneity observations motivate explicit constraints and uncertainty. No experimental scores, ranks, frequencies or universal GGG requirement are installed. |

NEB and TriLink retain rights in their content. Their pages are live documents;
an immutable upstream version/hash is not asserted. The profile identity binds
our exact factual curation and reviewed source metadata, not future web content
or source authenticity. Product names identify the documented product; no patent,
product-use or technology licence is implied by this MIT implementation.

Conrad et al., *Communications Biology*3:439 (2020), is © the authors and
[CC BY4.0](https://creativecommons.org/licenses/by/4.0/). Our qualitative summary
is independently worded; no numerical ranking table is copied.

[Trepotec et al. 2019](https://doi.org/10.1261/rna.069286.118), *RNA*25:507–518,
describes segmented poly(A) experiments. The publisher/author-posted article
identifies CC BY-NC4.0 after its initial exclusive period. No article, figure,
measured table, plasmid sequence or host-specific recombination model is copied
into this MIT module. Poly(A) output is arithmetic on the caller's exact declared
sequence only. A separately licensed, context-bound evidence admission is needed
before providing measured propagation guidance.
