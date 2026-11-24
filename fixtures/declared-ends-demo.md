# Ligation from explicitly declared starting-material ends

Open `declared-ends-demo.fa`, make `declared-target` active, choose **Assemble**
and the ordered-ligation workflow. Use a linear product and two rows in order:
`declared-left`, then `declared-right`, both forward.

For each row, choose **Declare intact linear duplex**. Explicitly set each source
end to **5′ phosphate / 3′ hydroxyl**. Leave the separate terminal-treatment
control at none: this example declares starting material rather than describing
a phosphorylation reaction. Require sealed junctions and compute evidence.

The independently expected product is `AAAACCCC` (8 bases), with one sealed blunt
junction and no retained nicks. The original source DNA and files are unchanged.
The candidate and HTML/PDF protocol retain the terminal declarations as explicit
operator hypotheses, alongside source revisions, orientation and target checking.

Selecting ordinary **Whole linear DNA** instead leaves free-end chemistry
unknown and cannot establish the same sealed ligation. Likewise, leaving a
required declared group unknown keeps its uncertainty; a blocked 3′ group is not
silently converted to OH. Nick acceptance does not invent phosphorylation or a
bond where neither strand can seal.

This is a short coordinate/chemistry example, not a laboratory protocol. A
terminal declaration does not measure chemistry, synthesize DNA, repair nicks or
flaps, clear a reference-only warning, or establish enzyme activity. It accepts
only an eligible original intact linear DNA reference. Imported single-stranded,
circular or detached physical-intermediate sources require their actual supported
preparation graph; metadata is not a route to restore physical state.
