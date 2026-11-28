#!/usr/bin/env python3
"""Independent maintainer count oracle, Python standard library only.

Not imported/called by Collagene. Uses the pinned source's direct CDS exclusion
sets and simple single-range coordinate binding, not the Rust generic parser,
translation engine, skip-priority implementation or profile serializer.
"""
import argparse
import collections
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import re

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--data', type=Path, default=Path(__file__).resolve().parent)
parser.add_argument('--profile', type=Path)
parser.add_argument('--output', type=Path)
parser.add_argument('--write-fixture', type=Path)
args = parser.parse_args()
SOURCE = args.data
RUST = args.profile or SOURCE / 'profile.json'
SOURCE_FILES = {
    'GCF_000005845.2_ASM584v2_cds_from_genomic.fna.gz': 'cds.fna.gz',
    'GCF_000005845.2_ASM584v2_genomic.gff.gz': 'genomic.gff.gz',
    'GCF_000005845.2_ASM584v2_assembly_report.txt': 'assembly_report.txt',
    'annotation_hashes.txt': 'annotation_hashes.txt',
}
PINNED = {
    'GCF_000005845.2_ASM584v2_cds_from_genomic.fna.gz': (
        '8b7b48b7b905313082141db4b77807b7aafebca7a5c40bbd42832baacd1cfd9b',
        '8800d3da1d7f2e7bb2e4a51ab1246701',
        'a5efdd03ccf26b27d8319b4e8bcfdde1da7df1b2a4baaf3edd1431c21128a0cd'),
    'GCF_000005845.2_ASM584v2_genomic.gff.gz': (
        'afdf03dc1d06e423d874ee29d9e0df14f5d32baf893f4da9f93aa721eb5f495e',
        '0f52ffc94af5ddf544ff89cc6f546b0c',
        '68d25735ccde6d37af68613dca57a588ae69423f0ecee99b0f29d596eff7accb'),
    'GCF_000005845.2_ASM584v2_assembly_report.txt': (
        'ca161bcfe8d1842b2160ca3fa20b480956d8242895115c53f5f819a11e895362', None, None),
    'annotation_hashes.txt': (
        'f4bf0e53a7b4ea673c0bad04a677893648b297764c74feed2c0e9521861e5f3e',
        '3809c06d06018f13ad258cbea5f8ea0c', None),
}
source_text = {}
for name, (sha256, md5, expanded_sha256) in PINNED.items():
    raw = (SOURCE / SOURCE_FILES[name]).read_bytes()
    assert len(raw) < 20 * 1024 * 1024
    assert hashlib.sha256(raw).hexdigest() == sha256, name
    if md5:
        assert hashlib.md5(raw).hexdigest() == md5, name
    expanded = gzip.decompress(raw) if name.endswith('.gz') else raw
    assert len(expanded) < 20 * 1024 * 1024
    if expanded_sha256:
        assert hashlib.sha256(expanded).hexdigest() == expanded_sha256, name
    source_text[name] = expanded.decode('ascii')

# Build direct locus/product exclusion sets from the captured GFF, independently
# of the FASTA text. A gene can encode both ordinary and exceptional products:
# keys must include protein identity, rather than discarding an entire locus.
groups = collections.defaultdict(list)
for line in source_text['GCF_000005845.2_ASM584v2_genomic.gff.gz'].splitlines():
    if line.startswith('#'):
        continue
    fields = line.split('\t')
    if fields[2] != 'CDS':
        continue
    attributes = dict(item.split('=', 1) for item in fields[8].split(';'))
    assert fields[0] == 'NC_000913.3'
    assert attributes['transl_table'] == '11'
    key = (attributes['locus_tag'], attributes.get('protein_id', ''))
    groups[key].append((fields, attributes))
assert len(groups) == 4318
pseudogenes = {key for key, rows in groups.items()
              if any(attrs.get('pseudo') == 'true' for _, attrs in rows)}
exceptions = {key for key, rows in groups.items()
              if any('exception' in attrs or 'transl_except' in attrs for _, attrs in rows)}
assert len(pseudogenes) == 18 and len(exceptions) == 6
assert not pseudogenes & exceptions

# Code-11 starts and ordinary stops are explicit assertions for this one pinned
# dataset; no copied codon-to-amino-acid translation table or bio library is used.
STARTS = {'TTG', 'CTG', 'ATT', 'ATC', 'ATA', 'ATG', 'GTG'}
STOPS = {'TAA', 'TAG', 'TGA'}
counts = collections.Counter()
starts = collections.Counter()
stops = collections.Counter()
seen = set()
admitted = 0
fasta = source_text['GCF_000005845.2_ASM584v2_cds_from_genomic.fna.gz']
for record in re.split(r'(?m)^>', fasta)[1:]:
    lines = record.splitlines()
    header, sequence = lines[0], ''.join(lines[1:])
    def attribute(name, default=None):
        match = re.search(r'\[' + re.escape(name) + r'=([^\]]*)\]', header)
        return match.group(1) if match else default
    key = (attribute('locus_tag'), attribute('protein_id', ''))
    assert key in groups and key not in seen, (key, header)
    seen.add(key)
    if key in pseudogenes or key in exceptions:
        continue
    assert attribute('pseudo') is None
    assert attribute('exception') is None and attribute('transl_except') is None
    assert attribute('gbkey') == 'CDS' and key[1]
    rows = groups[key]
    assert len(rows) == 1
    fields, attrs = rows[0]
    assert fields[7] == '0'
    assert not any(k in attrs for k in ('partial', 'start_range', 'end_range'))
    start, end = int(fields[3]), int(fields[4])
    expected_location = f'{start}..{end}'
    if fields[6] == '-':
        expected_location = f'complement({expected_location})'
    else:
        assert fields[6] == '+'
    assert attribute('location') == expected_location
    assert len(sequence) == end - start + 1
    assert len(sequence) >= 6 and len(sequence) % 3 == 0
    assert set(sequence) <= set('ACGT')
    assert sequence[:3] in STARTS
    assert sequence[-3:] in STOPS
    internal = [sequence[i:i+3] for i in range(3, len(sequence)-3, 3)]
    assert not STOPS & set(internal)
    counts.update(internal)
    starts[sequence[:3]] += 1
    stops[sequence[-3:]] += 1
    admitted += 1
assert seen == set(groups)
assert admitted == 4294 and sum(counts.values()) == 1322130
assert sum(starts.values()) == admitted == sum(stops.values())
order = [''.join(codon) for codon in itertools.product('ACGT', repeat=3)]
arrays = {
    'counts': [counts[codon] for codon in order],
    'initiation_counts': [starts[codon] for codon in order],
    'terminal_stop_counts': [stops[codon] for codon in order],
}
rust = json.loads(RUST.read_text())
assert rust['identity'] == '52d52265b4f997e5f95c2d4178bcfc575e5002bef278a7a327a9efe6f8d61f43'
for field, values in arrays.items():
    mismatches = [(codon, values[i], rust[field][i]) for i, codon in enumerate(order)
                  if values[i] != rust[field][i]]
    assert not mismatches, (field, mismatches)
assert rust['observed_loci'] == len(seen)
assert rust['admitted_loci'] == admitted
assert rust['internal_codons'] == sum(counts.values())

result = {
    'method': 'Independent Python stdlib: direct GFF locus-product exclusion sets; per-CDS codon tallies; explicit code-11 boundary assertions',
    'compared_values': 192,
    'all_values_match': True,
    'assembly': 'GCF_000005845.2',
    'rust_profile_identity': rust['identity'],
    'source_sha256': {name: values[0] for name, values in PINNED.items()},
    'codon_order': ' '.join(order),
    'observed_loci': len(seen), 'admitted_loci': admitted,
    'internal_codons': sum(counts.values()), 'pseudogenes': len(pseudogenes),
    'translation_exceptions': len(exceptions), **arrays,
}
if args.output:
    args.output.write_text(json.dumps(result, indent=2) + '\n')
# Compact numeric oracle fixture, with the independent method and raw source pins.
fixture = json.dumps(result, indent=2)
fixture = re.sub(r'\[\n\s+([0-9,\s]+)\n\s+\]',
                 lambda m: '[' + re.sub(r'\s+', '', m.group(1)) + ']', fixture)
if args.write_fixture:
    args.write_fixture.write_text(fixture + '\n')
print(json.dumps({'matched_values': 192, 'admitted': admitted,
                  'internal_codons': sum(counts.values()),
                  'starts': dict(sorted(starts.items())),
                  'stops': dict(sorted(stops.items())),
                  'fixture': str(args.write_fixture) if args.write_fixture else None}, indent=2))
