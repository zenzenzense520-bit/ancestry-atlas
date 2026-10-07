"""Deterministic array-export QC; no medical or ancestry inference."""
from pathlib import Path
from collections import Counter, defaultdict
import csv
import hashlib
import json
import re

DNA = set('ACGT')
INDEL = set('ID')

def kind(gt):
    if gt == '--':
        return 'no_call'
    if len(gt) == 2 and set(gt) <= DNA:
        return 'acgt'
    if len(gt) == 2 and set(gt) <= INDEL:
        return 'indel_code'
    return 'invalid'

def normalize(gt):
    return ''.join(sorted(gt)) if kind(gt) in ('acgt', 'indel_code') else gt

def parse(path):
    groups = defaultdict(list)
    rows = []
    with Path(path).open(newline='', encoding='utf-8-sig') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        if reader.fieldnames != ['gxid', 'chromosome', 'position', 'genotype']:
            raise ValueError(f'Unsupported columns: {reader.fieldnames}')
        for line, row in enumerate(reader, 2):
            if None in row or any(v is None for v in row.values()):
                raise ValueError(f'Wrong column count at line {line}')
            row['position'] = int(row['position'])
            if row['position'] <= 0:
                raise ValueError(f'Invalid position at line {line}')
            row['source_line'] = line
            row['call_kind'] = kind(row['genotype'])
            row['normalized_genotype'] = normalize(row['genotype'])
            rows.append(row)
            groups[row['chromosome'], row['position']].append(row)
    return rows, groups

def consensus(group):
    dna = {x['normalized_genotype'] for x in group if x['call_kind'] == 'acgt'}
    indels = {x['normalized_genotype'] for x in group if x['call_kind'] == 'indel_code'}
    if any(x['call_kind'] == 'invalid' for x in group):
        return 'invalid_group', None
    if dna and indels:
        return 'mixed_variant_encoding', None
    if len(dna) > 1:
        return 'acgt_discordant', None
    if len(indels) > 1:
        return 'indel_code_discordant', None
    if dna:
        return 'acgt_consensus', next(iter(dna))
    if indels:
        return 'indel_only', next(iter(indels))
    return 'no_call', None

def analyze(path, outdir):
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    rows, groups = parse(path)
    states = Counter()
    normalized = []
    issues = []
    duplicate_rows = []
    chroms = []
    order = list(map(str, range(1, 23))) + ['X', 'Y', 'MT']
    for (chrom, position), group in groups.items():
        state, gt = consensus(group)
        states[state] += 1
        record = {'chromosome': chrom, 'position': position,
                  'ids': ';'.join(x['gxid'] for x in group),
                  'source_lines': ';'.join(str(x['source_line']) for x in group),
                  'source_genotypes': ';'.join(x['genotype'] for x in group),
                  'state': state, 'normalized_genotype': gt or ''}
        if state == 'acgt_consensus':
            normalized.append(record)
        if state in ('invalid_group', 'mixed_variant_encoding', 'acgt_discordant', 'indel_code_discordant'):
            issues.append(record)
        if len(group) > 1:
            duplicate_rows.append(record)
    for chrom in order:
        rr = [r for r in rows if r['chromosome'] == chrom]
        if not rr and chrom != 'Y':
            continue
        unique = [(k, v) for k, v in groups.items() if k[0] == chrom]
        clean = [x for x in normalized if x['chromosome'] == chrom]
        raw_no_call = sum(r['call_kind'] == 'no_call' for r in rr)
        raw_acgt = [r for r in rr if r['call_kind'] == 'acgt']
        chroms.append({'chromosome': chrom, 'raw_rows': len(rr), 'unique_coordinates': len(unique),
                      'raw_no_call': raw_no_call,
                      'raw_no_call_rate': raw_no_call / len(rr) if rr else None,
                      'raw_acgt_calls': len(raw_acgt), 'clean_acgt_coordinates': len(clean),
                      'raw_acgt_heterozygotes': sum(len(set(r['genotype'])) == 2 for r in raw_acgt),
                      'clean_heterozygotes': sum(len(set(r['normalized_genotype'])) == 2 for r in clean),
                      'clean_heterozygosity': sum(len(set(r['normalized_genotype'])) == 2 for r in clean) / len(clean) if clean else None,
                      'position_min': min((r['position'] for r in rr), default=None),
                      'position_max': max((r['position'] for r in rr), default=None)})
    autosomal = [x for x in normalized if x['chromosome'] in set(map(str, range(1, 23)))]
    raw_counts = Counter(r['call_kind'] for r in rows)
    qc = {'sample_alias': 'Sample01', 'version': '0.5',
          'input_sha256': hashlib.sha256(Path(path).read_bytes()).hexdigest(),
          'input_bytes': Path(path).stat().st_size, 'columns': ['gxid', 'chromosome', 'position', 'genotype'],
          'raw_rows': len(rows), 'unique_coordinates': len(groups),
          'coordinate_duplicate_groups': len(duplicate_rows),
          'redundant_coordinate_rows': len(rows) - len(groups),
          'exact_duplicate_ids': sum(v - 1 for v in Counter(r['gxid'] for r in rows).values() if v > 1),
          'suffixed_id_rows': sum(bool(re.fullmatch(r'rs\d+\.\d+', r['gxid'])) for r in rows),
          'non_rs_id_rows': sum(not r['gxid'].startswith('rs') for r in rows),
          'raw_call_kinds': dict(raw_counts),
          'raw_no_call_rate': raw_counts['no_call'] / len(rows),
          'raw_reported_nonmissing_rate': 1 - raw_counts['no_call'] / len(rows),
          'coordinate_states': dict(states),
          'clean_acgt_coordinates': len(normalized),
          'autosomal_clean_acgt_coordinates': len(autosomal),
          'autosomal_heterozygosity': sum(len(set(x['normalized_genotype'])) == 2 for x in autosomal) / len(autosomal),
          'duplicate_encoding_issues': issues, 'chromosomes': chroms,
          'y_status': 'NO_Y_ROWS_IN_EXPORT', 'mt_status': 'TWO_SITES_INSUFFICIENT_FOR_HAPLOGROUP',
          'provider': None, 'chip_version': None, 'reported_build': None,
          'reported_strand': None, 'no_call_uses': '--',
          'limitations': ['No Y rows describe export coverage, not a biological absence of the Y chromosome.',
                          'The no-call rate applies only to rows in this export, not the entire original chip.',
                          'Coordinate consensus groups repeated rows; mixed indel/base encodings are excluded.',
                          'I/D encodings have no explicit inserted/deleted sequence and are not converted to SNP alleles.',
                          'Heterozygosity is descriptive for this selected array subset, not an inbreeding or contamination estimate.']}
    (out / 'sample_qc.json').write_text(json.dumps(qc, indent=2, ensure_ascii=False))
    def write(name, data, fields):
        with (out / name).open('w', newline='', encoding='utf-8') as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t')
            writer.writeheader()
            writer.writerows(data)
    fields = ['chromosome', 'position', 'ids', 'source_lines', 'source_genotypes', 'state', 'normalized_genotype']
    write('normalized_acgt.tsv', normalized, fields)
    write('duplicate_coordinate_audit.tsv', duplicate_rows, fields)
    write('complex_coordinate_issues.tsv', issues, fields)
    write('chromosome_qc.tsv', chroms, list(chroms[0]))
    return qc, normalized

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('input')
    parser.add_argument('--out', default='results')
    args = parser.parse_args()
    qc, _ = analyze(args.input, args.out)
    print(json.dumps({k: qc[k] for k in ['raw_rows', 'unique_coordinates', 'clean_acgt_coordinates',
                     'autosomal_clean_acgt_coordinates', 'raw_no_call_rate', 'coordinate_states']}, indent=2))
