# Isolated internal DNA mismatch parameters

Only numerical parameter facts from `DNA_IMM1` in Biopython **1.86**
`Bio/SeqUtils/MeltingTemp.py` are retained in `src/thermodynamics/parameters.rs`.
No upstream thermodynamics, alignment, or solver implementation is copied.
Collagene's implementation is MIT; the incorporated data retain the upstream
permission and attribution below. This does not relabel upstream material MIT.

Source: https://raw.githubusercontent.com/biopython/biopython/biopython-186/Bio/SeqUtils/MeltingTemp.py

Captured 2026-11-28, 43,151 bytes, SHA256
`4f3cbc585ffd68d54386d40f5edb113a7c2eae4aa4388c6b496bceecdec06500`.

Copyright 2004-2008 by Sebastian Bassi.
Copyright 2013-2018 by Markus Piotrowski.
All rights reserved.

The file explicitly permits choosing the Biopython License Agreement or BSD
3-Clause License. We select the Biopython License Agreement. The full original
license file, including both alternatives, is retained as
[BIOPYTHON-LICENSE.txt](BIOPYTHON-LICENSE.txt), captured from
https://raw.githubusercontent.com/biopython/biopython/biopython-186/LICENSE.rst
(3,267 bytes, SHA256
`d4bbb0c4e5345e6f2508c3b93fb247fcb2db95783fb8a584c23cd88e3c03df6c`).

Permission to use, copy, modify, and distribute this software and its
documentation with or without modifications and for any purpose and
without fee is hereby granted, provided that any copyright notices
appear in all copies and that both those copyright notices and this
permission notice appear in supporting documentation, and that the
names of the contributors or copyright holders not be used in
advertising or publicity pertaining to distribution of the software
without specific prior permission.
THE CONTRIBUTORS AND COPYRIGHT HOLDERS OF THIS SOFTWARE DISCLAIM ALL
WARRANTIES WITH REGARD TO THIS SOFTWARE, INCLUDING ALL IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS, IN NO EVENT SHALL THE
CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY SPECIAL, INDIRECT
OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS
OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE
OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE
OR PERFORMANCE OF THIS SOFTWARE.

Selection: retain exactly 48 A/C/G/T dinucleotide records with one canonical pair
and one mismatch; discard every inosine record and the three tandem-mismatch
records `GG/TT`, `GT/TG`, `TG/GT`. Values are represented exactly as integer tenths
of kcal/mol (enthalpy) and cal/(mol K) (entropy). The upstream dictionary was read
as a literal data object during curation, never executed or required at runtime.
The model admits one internal mismatch and two canonical pairs at each end;
its admission policy is narrower than the table's full possible uses.

The original experimental papers are cited in
[the research audit](../../../../docs/audits/oligo-mismatch-thermodynamics-current.md).
Publisher article/PDF assets and prose are not redistributed here. No claim of
MIT licensing for those articles or their supplemental assets is made.
