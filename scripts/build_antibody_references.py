#!/usr/bin/env python3
"""Build the pinned CC0 human IMGT V-prefix/J-FR4 corpus from local OGRDB AIRR JSON.

Development tooling only: no network, scientific executable or Python runtime is
used by Collagene. Inputs are the three original API response files; obtain their
exact URLs from the output source manifest. Unsupported input fails explicitly.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path

RELEASES = [
    ('G00060', 10, 'IGH', 'b2784200e1942b63e4f8423bc46a7db328eff0842345302ccebf440b8ef16c26'),
    ('G00061', 5, 'IGK', 'fd5fd5d696ae20b3dc203c61da8083102ee282e4d5255a274f3f882806879453'),
    ('G00062', 4, 'IGL', '62b27be9656d051a11845695c08fc73bae22c3e298e0abb0fab8a2cf31f5b64f'),
]
# Standard code, grouped by amino acid rather than copied from a runtime table.
CODONS = {
    'F': 'TTT TTC', 'L': 'TTA TTG CTT CTC CTA CTG', 'I': 'ATT ATC ATA',
    'M': 'ATG', 'V': 'GTT GTC GTA GTG', 'S': 'TCT TCC TCA TCG AGT AGC',
    'P': 'CCT CCC CCA CCG', 'T': 'ACT ACC ACA ACG', 'A': 'GCT GCC GCA GCG',
    'Y': 'TAT TAC', '*': 'TAA TAG TGA', 'H': 'CAT CAC', 'Q': 'CAA CAG',
    'N': 'AAT AAC', 'K': 'AAA AAG', 'D': 'GAT GAC', 'E': 'GAA GAG',
    'C': 'TGT TGC', 'W': 'TGG', 'R': 'CGT CGC CGA CGG AGA AGG',
    'G': 'GGT GGC GGA GGG',
}
TRANSLATE = {c: aa for aa, codons in CODONS.items() for c in codons.split()}
assert set(TRANSLATE) == {''.join(c) for c in itertools.product('ACGT', repeat=3)}


def require(ok, message):
    if not ok:
        raise ValueError(message)


def build(paths):
    sources, v_entries, j_entries = [], [], []
    for index, ((gid, version, locus, digest), path) in enumerate(zip(RELEASES, paths, strict=True)):
        raw = path.read_bytes()
        require(len(raw) <= 16 << 20, 'source byte bound')
        require(hashlib.sha256(raw).hexdigest() == digest, f'{gid}: source hash changed')
        sets = json.loads(raw)['GermlineSet']
        require(len(sets) == 1, 'one pinned set required')
        s = sets[0]
        require((s['germline_set_id'], s['release_version'], s['locus'], s['species']['id']) ==
                ('OGRDB:' + gid, version, locus, 'NCBITAXON:9606'), 'source identity')
        alleles = s['allele_descriptions']
        excluded = []
        kept = 0
        for a in alleles:
            if a['sequence_type'] not in ['V', 'J'] or a['functional'] is not True:
                excluded.append({'allele_ref': a['allele_description_ref'],
                                 'reason': 'not-functional-V-or-J', 'sequence_type': a['sequence_type']})
                continue
            require(a['locus'] == locus and a['species']['id'] == 'NCBITAXON:9606', 'allele identity')
            if a['sequence_type'] == 'J':
                coding = a['coding_sequence']
                require(a['sequence'][a['gene_start']-1:a['gene_end']] == coding, 'J gene/coding interval')
                anchor = a['j_cdr3_end']-1
                frame = a['j_codon_frame']
                # The pinned OGRDB release values use a coding-sequence anchor
                # and a 1-based complete-codon start. Preserve the raw fields;
                # this is not a generic AIRR coordinate reinterpretation.
                require(frame in [1, 2, 3] and anchor >= 0 and anchor % 3 == frame-1,
                        'J provider anchor/frame agreement')
                require(set(coding) <= set('ACGT') and anchor < len(coding), 'J coding alphabet/anchor')
                complete = (len(coding)-anchor)//3*3
                dna = coding[anchor:anchor+complete]
                trailing = coding[anchor+complete:]
                protein = ''.join(TRANSLATE[dna[n:n+3]] for n in range(0, len(dna), 3))
                require(len(protein) == (11 if locus == 'IGH' else 10) and len(trailing) == 1,
                        'unsupported J terminal extent')
                require(protein[0] == ('W' if locus == 'IGH' else 'F') and protein[1] == 'G'
                        and protein[3] == 'G' and '*' not in protein, 'J canonical motif/translation')
                j_entries.append({
                    'source': index, 'allele_ref': a['allele_description_ref'],
                    'allele_id': a['allele_description_id'], 'label': a['label'],
                    'release_version': a['release_version'],
                    'original_sequence': a['sequence'], 'coding_nucleotides': coding,
                    'gene_start': a['gene_start'], 'gene_end': a['gene_end'],
                    'provider_codon_frame': frame, 'provider_j_cdr3_end': a['j_cdr3_end'],
                    'anchor_offset': anchor, 'fr4_nucleotides': dna, 'fr4_protein': protein,
                    'fr4_positions': [str(n) for n in range(118, 118+len(protein))],
                    'trailing_nucleotides': trailing,
                })
                continue
            delineations = [d for d in a['v_gene_delineations'] if d['delineation_scheme'] == 'IMGT']
            require(len(delineations) == 1, 'one explicit IMGT delineation required')
            d = delineations[0]
            require(d['alignment_labels'] == [str(n) for n in range(1, 105)], 'unsupported IMGT label set')
            aligned = d['aligned_sequence']
            require(len(aligned) >= 312 and set(aligned) <= set('ACGT.'), 'aligned DNA alphabet/length')
            # Labels cover V positions 1..104 only. Germline CDR3 tails have no
            # labels here and are deliberately excluded; never invent their map.
            aligned = aligned[:312]
            prefix = aligned.replace('.', '')
            require(d['unaligned_sequence'].startswith(prefix), 'gapped/ungapped prefix mismatch')
            require(a['coding_sequence'].startswith(prefix), 'coding/delineation prefix mismatch')
            codons = [aligned[n:n+3] for n in range(0, 312, 3)]
            require(all(c == '...' or '.' not in c for c in codons), 'partial-codon gap')
            protein = ''.join('-' if c == '...' else TRANSLATE[c] for c in codons)
            require('*' not in protein, 'stop in functional V prefix')
            ranges = []
            for kind, low, high in [(1, 27, 38), (2, 56, 65)]:
                start, end = d[f'cdr{kind}_start'], d[f'cdr{kind}_end']
                require((start-1) % 3 == 0 and end % 3 == 0, 'CDR codon boundary')
                start, end = (start-1)//3, end//3
                observed = [i for i, n in enumerate(
                    [n for n, aa in enumerate(protein, 1) if aa != '-']) if low <= n <= high]
                require(observed and (start, end) == (observed[0], observed[-1]+1), 'CDR/labels disagree')
                ranges.append({'kind': kind, 'start': start, 'end': end})
            require(d['fwr3_end'] == len(prefix), 'V-prefix terminus/delineation mismatch')
            v_entries.append({
                'source': index, 'allele_ref': a['allele_description_ref'],
                'allele_id': a['allele_description_id'], 'label': a['label'],
                'release_version': a['release_version'],
                'aligned_nucleotides': aligned, 'aligned_protein': protein,
                'alignment_labels': d['alignment_labels'], 'cdrs': ranges,
                'non_cysteine_anchors': [n for n in [23, 104] if protein[n-1] != 'C'],
            })
            kept += 1
        sources.append({
            'set_id': s['germline_set_id'], 'set_ref': s['germline_set_ref'],
            'release_version': version, 'release_date': s['release_date'],
            'species': s['species']['label'], 'taxonomy_id': s['species']['id'],
            'locus': locus, 'author': s['author'],
            'url': f'https://ogrdb.airr-community.org/api/germline/set/{gid}/{version}',
            'sha256': digest, 'original_records': len(alleles), 'v_records': kept,
            'j_records': sum(e['source'] == index for e in j_entries),
            'excluded': sorted(excluded, key=lambda x: x['allele_ref']),
        })
    v_entries.sort(key=lambda x: (x['source'], x['allele_ref']))
    j_entries.sort(key=lambda x: (x['source'], x['allele_ref']))
    return {'version': 2, 'license': 'CC0-1.0',
            'license_url': 'https://ogrdb.airr-community.org/render_page/licensing_statement.html',
            'scope': 'human-functional-V-prefix-and-J-FR4',
            'j_label_rule': 'IMGT-118-sequential-complete-J-codons',
            'sources': sources, 'v_entries': v_entries, 'j_entries': j_entries}


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('sources', nargs=3, type=Path, help='G00060 v10, G00061 v5, G00062 v4 original AIRR JSON')
    ap.add_argument('--output', required=True, type=Path)
    args = ap.parse_args()
    result = build(args.sources)
    args.output.write_text(json.dumps(result, separators=(',', ':'), ensure_ascii=True)+'\n')
    print(f"{len(result['v_entries'])} V and {len(result['j_entries'])} J references; {args.output.stat().st_size} bytes")
