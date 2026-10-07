"""Coordinate/allele harmonization and traceable public reference extraction."""
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path
import numpy as np
import pgenlib

ROOT = Path(__file__).resolve().parents[1]
COMPLEMENT = str.maketrans('ACGT', 'TGCA')

def write_tsv(path, rows):
    with path.open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)

def main():
    out, ref = ROOT / 'results', ROOT / 'references'
    psam = list(csv.DictReader((ref / 'phase3_orig.psam').open(), delimiter='\t'))
    panel = {r['sample']: r for r in csv.DictReader(
        (ref / 'integrated_call_samples_v3.20130502.ALL.panel').open(), delimiter='\t')}
    if len(psam) != 2504 or len({r['#IID'] for r in psam}) != 2504 or set(panel) != {r['#IID'] for r in psam}:
        raise ValueError('Reference IDs disagree with official Phase 3 panel')
    for row in psam:
        official = panel[row['#IID']]
        if row['Population'] != official['pop'] or row['SuperPop'] != official['super_pop']:
            raise ValueError('Reference labels disagree with official panel')
    samples = [{'reference_row_index': i, 'sample_id': r['#IID'], 'group': r['SuperPop'],
                'population': r['Population']} for i, r in enumerate(psam)]
    write_tsv(out / 'reference_samples.tsv', samples)
    clean = [r for r in csv.DictReader((out / 'normalized_acgt.tsv').open(), delimiter='\t')
             if r['chromosome'] in set(map(str, range(1, 23)))]
    target = {(r['chromosome'], r['position']): r for r in clean}
    candidates = defaultdict(list)
    headers, variant_count = [], 0
    with (ref / 'all_phase3.pvar').open() as f:
        for line in f:
            if line.startswith('#'):
                headers.append(line.rstrip())
                continue
            cols = line.rstrip('\n').split('\t', 5)
            key = cols[0], cols[1]
            if key in target:
                candidates[key].append({'reference_variant_index': variant_count,
                    'reference_id': cols[2], 'ref': cols[3], 'alt': cols[4]})
            variant_count += 1
    if not any('b37' in h for h in headers):
        raise ValueError('PVAR header does not identify the expected b37 reference')
    audit = []
    for row in clean:
        key = row['chromosome'], row['position']
        matches = candidates.get(key, [])
        snps = [r for r in matches if len(r['ref']) == len(r['alt']) == 1
                and set(r['ref'] + r['alt']) <= set('ACGT') and r['ref'] != r['alt']]
        gt = row['normalized_genotype']
        record = {'chromosome': key[0], 'position': int(key[1]), 'source_ids': row['ids'],
                  'source_lines': row['source_lines'], 'observed_genotype': gt,
                  'reference_variant_index': '', 'reference_id': '', 'ref': '', 'alt': '',
                  'alt_dosage': '', 'state': 'not_in_reference'}
        if matches:
            record['state'] = 'no_unique_biallelic_snp'
        if len(snps) == 1:
            record.update(snps[0])
            alleles = {record['ref'], record['alt']}
            flipped = gt.translate(COMPLEMENT)
            if alleles in (set('AT'), set('CG')):
                record['state'] = 'palindromic_excluded'
            elif set(gt) <= alleles and not set(flipped) <= alleles:
                record['state'] = 'direct_compatible'
                record['alt_dosage'] = gt.count(record['alt'])
            elif set(flipped) <= alleles and not set(gt) <= alleles:
                record['state'] = 'unique_complement'
                record['alt_dosage'] = flipped.count(record['alt'])
            else:
                record['state'] = 'incompatible'
        audit.append(record)
    eligible = sorted([r for r in audit if r['state'] in ('direct_compatible', 'unique_complement')],
                      key=lambda r: r['reference_variant_index'])
    indices = np.array([r['reference_variant_index'] for r in eligible], dtype=np.uint32)
    dosage = np.array([r['alt_dosage'] for r in eligible], dtype=np.int8)
    write_tsv(out / 'harmonization_full_audit.tsv', audit)
    write_tsv(out / 'reference_overlap_sites.tsv', eligible)
    np.savez_compressed(out / 'reference_overlap_index.npz', reference_variant_indices=indices,
                        target_dosages=dosage)
    # Multiallelic metadata is required by PGEN even though selected variants are SNPs.
    pvar = pgenlib.PvarReader(bytes(str(ref / 'all_phase3.pvar'), 'utf-8'))
    try:
        reader = pgenlib.PgenReader(bytes(str(ref / 'all_phase3.pgen'), 'utf-8'),
                                    raw_sample_ct=len(samples), pvar=pvar)
        try:
            if reader.get_variant_ct() != variant_count or reader.get_raw_sample_ct() != len(samples):
                raise ValueError('PGEN dimensions do not match PVAR/PSAM')
            matrix = np.empty((len(indices), len(samples)), dtype=np.int8)
            for start in range(0, len(indices), 4096):
                reader.read_list(indices[start:start + 4096], matrix[start:start + 4096])
        finally:
            reader.close()
    finally:
        pvar.close()
    np.savez_compressed(out / 'reference_overlap.npz', genotypes=matrix,
                        target_dosages=dosage, reference_variant_indices=indices)
    qc = {'reference_n': len(samples), 'reference_variant_count': variant_count,
          'target_autosomal_coordinates': len(clean), 'coordinate_hits': len(candidates),
          'harmonization': dict(Counter(r['state'] for r in audit)),
          'usable_before_maf_ld': len(indices), 'reference_missing_cells': int((matrix < 0).sum()),
          'reference_population_counts': dict(Counter(r['population'] for r in samples)),
          'official_panel_labels_verified': True, 'assembly': 'GRCh37',
          'download_manifest': json.loads((ref / 'download_manifest.json').read_text())}
    (out / 'expanded_reference_qc.json').write_text(json.dumps(qc, indent=2))
    print(json.dumps(qc, indent=2), flush=True)

if __name__ == '__main__':
    main()
