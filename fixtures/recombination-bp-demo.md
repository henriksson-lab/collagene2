# BP reference-exchange example

Import [recombination-bp-demo.gb](recombination-bp-demo.gb) to open four documents. `BP_insert` and `BP_donor` are the two sources. `BP_target` is a separately specified circular target; `BP_reciprocal` is the other expected, linear component. The target must never be supplied as a source or used to choose the component automatically. This fixture exercises the reference-sequence model and does not establish experimental success or end chemistry.

The library supports this example now. GUI, plan and MCP integration is being implemented. Use the opt-in `manufacturer-vector-contexts-att1-att2` profile, edition `2026-11-24.1`, reaction BP, preserving source order and orientation.

The donor is the manufacturer's 4,761-base [pDONR221 sequence text](https://documents.thermofisher.com/TFS-Assets/vectors/pdonr221_seq.txt), raw-file SHA256 `9f6fd695e2eaa95d7153dfde0786fb8d1a0ad8471ccde73a9689d4679af78424`. Its [map](https://documents.thermofisher.com/TFS-Assets/vectors/pdonr221_pdonrzeo_map.pdf) places attP1 at 570–801 and reverse attP2 at 2753–2984 (one-based inclusive). The synthetic insert uses the 25-base attB1/attB2 primer strings from printed page 13 of [MAN0000282 Rev1.0](https://assets.thermofisher.com/TFS-Assets/LSG/manuals/gatewayman.pdf), surrounding the explicitly chosen sequence `ATGCCGTAA`. These are factual sequence inputs, with no source code from another cloning application.

Let `d` be the donor sequence and `b` be `attB1 + ATGCCGTAA + reverse_complement(attB2)`. With zero-based half-open slices, an independent reference-switch calculation gives:

```text
target     = d[:651] + b[8:51] + d[2902:]
reciprocal = b[:8] + d[651:2902] + b[51:]
```

The target has 2,553 bases, canonical sequence SHA256 `8195a6c13ff3c9171355b6199cec0a8d8d5087a1216e564f17e661fde8bcca90`. The reciprocal has 2,267 bases, SHA256 `0c5532f1c8ecb3ad26a95f588fc384ee331cf11b317b1b8fa25677f1e8610667`. Their combined 4,820 bases equal the two source lengths. Circular products can have another origin while retaining the same biological reference sequence. The short insert annotation and a donor backbone landmark permit inspection of annotation mapping. No host selection, survival, expression or full functional-site generality is inferred.
