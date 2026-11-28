# Independent documentary image fixtures

These small files were constructed directly from format fields, not exported by
our image decoder/preview encoder or another cloning application.

- `known_2x3.png`: RGB pixels with red values10,20,30,40,50,60 in row order;
  green/blue zero. PNG uses one independently constructed stored DEFLATE block.
- `other_2x3.png`: same geometry/format, independent red values90,80,70,60,50,40.
  Used to ensure fresh preview reopening checks pixels rather than dimensions.
- `known_2x3_orientation6.tif`: the same RGB pixels, classic little-endian TIFF,
  one uncompressed strip, explicit orientation6. Documentary preview is3×2;
  the exact original remains2×3 with its original orientation tag.
- `known_8x8.jpg`: baseline grayscale JPEG with a unit quantization table,
  single-symbol DC0/EOB Huffman tables and one all-zero DCT block (pixels128).

The MCP fixtures check exact source bytes, dimensions/orientation, immutable
workspace publication and historical association. They confer no experimental,
sequence, physical-material, gel-comparison or verification authority.
