# Synthetic TA geometry example

Open `topo-ta-demo.gb` and select `TA_independent_target` as the active document.
In Cloning → TOPO choose `TA_vector_top`, `TA_vector_bottom` and `TA_PCR_template`
for their corresponding roles. The vector records are actual declared5′→3′
strand references, not an intact-duplex or supplier-state assertion.

Enter ordered forward primer **CGA**, reverse primer **CAT**, and annealing
length **3** for each. Explicitly declare the vector5′ phosphate/3′
phosphotyrosyl assumption, describe this synthetic preparation, and opt in to
one A addition at both PCR3′ ends. The primer pair reconstructs the literal
six-base template; the tailed insert strands are `CGAATGA` and `CATTCGA`.

Inspect both orientations. The forward hypothesis uses prepared vector strands
0/1 and insert strands2/3. Select its actual returned request, component0 and
reference seed strand0/index0; do not use the target to choose a contour.
The separately supplied target is the literal circular sequence
`AGGGGCATCCCTTCGAATGA`. The physical example retains **two unsealed nicks**;
the strict target policy must reject those. Explicitly allowing nicks records a
different policy, rather than repairing them. Review all other physical evidence
and retained components before accepting any graph comparison.

This is a small synthetic coordinate/chemistry example. It establishes no
commercial kit identity, cloning efficiency, cellular repair, functional DNA or
experimental recommendation. The ordered plan currently remains Incomplete
because downstream physical graph composition and plan-level physical target
acceptance are still separate implementation work. Browser validation is deferred.
