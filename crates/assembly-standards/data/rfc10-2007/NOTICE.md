# RFC10 part-sequence format: curation and boundary

Primary definition: Tom Knight, **Draft Standard for Biobrick Biological Parts**,
3 May 2007, BioBricks Foundation RFC10.

- [MIT DSpace publication and MD5 metadata](https://dspace.mit.edu/entities/publication/764237ec-294e-4c93-b224-2703800bf358), persistent handle `1721.1/45138`.
- [Author's exact OpenWetWare revision 297083](https://openwetware.org/mediawiki/index.php?title=The_BioBricks_Foundation:BBFRFC10&oldid=297083).
- [Readable standard text](https://openwetware.org/wiki/The_BioBricks_Foundation:BBFRFC10).
- [OpenWetWare content license](https://openwetware.org/wiki/OpenWetWare:Copyrights):
  CC BY-SA 3.0 or later, alternatively GFDL 1.2 or later. These upstream text
  terms are distinct from this crate's independently authored MIT code/fact
  selection; no standard prose or upstream implementation is copied into it.

Reviewed 2026-11-28. MIT's original raw-text download returned HTTP429. The official
text was read through the documentation tool, so `upstream_sha256` is absent with
an explicit reason. DSpace's published MD5 is not represented as a locally verified
SHA-256. The profile content identity binds our exact curated definition and the
current pinned enzyme catalog, not the authenticity of a remote page.

`profile.json` independently selects the concrete sequence requirements of
sections 1–3: five mandatory internal sites, five advisory site identities
(SapI's reverse complement is searched too), exact format additions, required
ATG start for the declared coding class and preferred in-frame TAATAA ending.
The normal prefix has 22 nt; the coding prefix has 20 nt before the part's ATG;
the suffix has 21 nt. The suffix's NotI site overlaps the PstI-leading C.
XbaI and SpeI's declared CTAG overhang signatures follow the exact cut geometry
already curated in the separately pinned enzyme catalog. All original standard
sequence strings are in 5′→3′ reference orientation; no primer is generated here.

This is deliberately the **2007 part-sequence profile**, not a claim of all
RFC10 compliance or modern iGEM submission eligibility. The standard's vector,
primer and host sections are outside this profile's scope, not dummy future
implementations. Its registry/host requirements do not enter a sequence-only
conformance result. No biological part ordering is prescribed by sections 1–3;
our grammar encodes either declared class followed by either class, with at least
one part. Operational part-count limits do not invent a biological maximum.

Newer summaries sometimes list only four mandatory sites. This source explicitly
includes NotI, so its absence is required under this exact profile. Optional
restriction sites and the preferred stop pattern remain advisory. A later
standard revision needs a new source-bound profile, not a silent alteration or
fallback reader.

No BioBrick/SEVA/toolkit sequence collection or third-party application code is
redistributed. Names identify standards without endorsement. Further profile
tables require independent data-specific licensing and source review; an MIT
software repository does not automatically license every incorporated toolkit
sequence or borrowed table.
