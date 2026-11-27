#!/usr/bin/env python3
"""Build a fixed CC0 RCSB chain-level assessment panel from pinned local responses.

Only deposited sequence/molecule identity facts are retained. Integrated IMGT,
SAbDab, Pfam and other third-party annotations are not copied or used as truth.
No network requests or scientific executables are used by this generator.
"""
import argparse
import hashlib
import json
from pathlib import Path

# Reviewed chain identities from archival entities and the corresponding entry
# publications. These labels are set before inspecting algorithm results.
TRUTH = {
    '1AO7': ['Other', 'Other', 'Other', 'Other', 'Other'],
    '1BJ1': ['IgLight', 'IgHeavy', 'Other'],
    '1CD8': ['Other'], '1CDH': ['Other'],
    '1D9K': ['Other', 'Other', 'Other', 'Other', 'Other'],
    '1EMA': ['Other'], '1HSA': ['Other', 'Other', 'Other'],
    '1HZH': ['IgHeavy', 'IgLight'], '1LYZ': ['Other'],
    '1N8Z': ['IgLight', 'IgHeavy', 'Other'], '1TEN': ['Other'],
    '5GGQ': ['IgLight', 'IgHeavy'],
    '5GGS': ['IgHeavy', 'IgLight', 'Other'],
    '5GGT': ['Other', 'IgHeavy', 'IgLight'],
}
ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT/'crates/antibody-numbering/data/rcsb-assessment/sources.json'


def build(directory):
    manifest = json.loads(MANIFEST.read_text())
    if manifest['version'] != 1:
        raise ValueError('unsupported source manifest')
    records = {}
    for name, expected in manifest['inputs'].items():
        with (directory/name).open('rb') as source:
            raw = source.read((4 << 20)+1)
        if len(raw) > 4 << 20 or hashlib.sha256(raw).hexdigest() != expected['sha256']:
            raise ValueError(f'source hash/size: {name}')
        records[name] = json.loads(raw)
    cases, excluded, unique = [], [], {}
    for pdb, labels in sorted(TRUTH.items()):
        entry = records[pdb+'-entry.json']
        if entry['rcsb_entry_container_identifiers']['polymer_entity_ids'] != [str(n) for n in range(1, len(labels)+1)]:
            raise ValueError('unexpected entity list')
        for entity_id, truth in enumerate(labels, 1):
            key = f'{pdb}-{entity_id}-entity.json'
            entity = records[key]
            seq = entity['entity_poly']['pdbx_seq_one_letter_code_can']
            if entity['entity_poly']['type'] != 'polypeptide(L)' or not seq.isascii():
                raise ValueError('unsupported polymer')
            source = {'pdb': pdb, 'entity': str(entity_id),
                      'description': entity['rcsb_polymer_entity']['pdbx_description'],
                      'sequence_source': key, 'entry_source': pdb+'-entry.json'}
            row = {'id': f'{pdb}:{entity_id}', 'truth': truth,
                   'sequence': seq, 'sequence_sha256': hashlib.sha256(seq.encode()).hexdigest(),
                   'sources': [source]}
            reason = 'shorter-than-50-residues' if len(seq) < 50 else 'ambiguous-canonical-sequence' if set(seq)-set('ACDEFGHIKLMNPQRSTVWY') else None
            if reason:
                row['reason'] = reason
                excluded.append(row)
            elif seq in unique:
                previous = unique[seq]
                if previous['truth'] != truth:
                    raise ValueError('duplicate sequence has conflicting truth')
                previous['sources'].append(source)
            else:
                cases.append(row)
                unique[seq] = row
    return {'version': 1, 'license': 'CC0-1.0',
            'license_url': 'https://www.rcsb.org/pages/policies',
            'retrieved': '2026-11-27',
            'truth_scope': 'deposited-chain-identity-only-not-residue-numbering',
            'inputs': manifest['inputs'], 'cases': cases, 'excluded': excluded}


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('input_directory', type=Path)
    ap.add_argument('--output', required=True, type=Path)
    args = ap.parse_args()
    data = build(args.input_directory)
    args.output.write_text(json.dumps(data, indent=2)+'\n', encoding='utf-8', newline='\n')
    print(f"{len(data['cases'])} distinct cases; {len(data['excluded'])} explicit exclusions")
