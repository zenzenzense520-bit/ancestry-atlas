"""Exploratory, reference-only PCA. Uses peddy's public 1000G genotype panel.

This is an independently implemented projection, not a peddy CLI run,
not a supervised admixture model, and not an ethnic-identity classifier.
"""
import ast
import csv
import gzip
import hashlib
import json
import platform
from collections import Counter
from pathlib import Path
import sys

import numpy as np
import sklearn
from sklearn.decomposition import PCA
from sklearn.model_selection import train_test_split

from qc import analyze

GROUPS = ['AFR', 'AMR', 'EAS', 'EUR', 'SAS']
SEED = 20261005

def load_reference(directory):
    root = Path(directory)
    sites = [line.strip().split(':') for line in (root / 'GRCH37.sites').read_text().splitlines()]
    data = np.frombuffer(gzip.decompress((root / 'GRCH37.sites.bin.gz').read_bytes()), dtype=np.uint8)
    if len(data) != len(sites) * 2504:
        raise ValueError('Reference genotype/site dimensions disagree')
    matrix = data.reshape(len(sites), 2504).T
    # Read label constants as data; never execute downloaded package code.
    tree = ast.parse((root / 'pca.py').read_text())
    text = next(ast.literal_eval(n.value) for n in tree.body if isinstance(n, ast.Assign)
                and any(isinstance(t, ast.Name) and t.id == '_str' for t in n.targets))
    labels = np.array([int(v) for v in text.split('|')])
    if len(labels) != 2504 or not set(np.unique(matrix)) <= {0, 1, 2}:
        raise ValueError('Unexpected reference genotype encoding or label dimension')
    return sites, matrix, labels

def harmonize(clean, sites):
    by_coord = {(x['chromosome'], int(x['position'])): x for x in clean}
    translation = str.maketrans('ACGT', 'TGCA')
    idx, target, audit = [], [], []
    for i, (chrom, pos, ref, alt) in enumerate(sites):
        row = by_coord.get((chrom, int(pos)))
        if row is None or chrom not in set(map(str, range(1, 23))):
            continue
        gt = row['normalized_genotype']
        flipped = gt.translate(translation)
        state = 'incompatible'
        if {ref, alt} in (set('AT'), set('CG')):
            state = 'palindromic_excluded'
        elif set(gt) <= {ref, alt} and not set(flipped) <= {ref, alt}:
            state = 'direct_compatible'
        elif set(flipped) <= {ref, alt} and not set(gt) <= {ref, alt}:
            state = 'unique_complement'
            gt = flipped
        record = {'reference_site_index': i, 'chromosome': chrom, 'position': int(pos),
                  'ref': ref, 'alt': alt, 'source_ids': row['ids'],
                  'source_lines': row['source_lines'], 'observed_genotype': row['normalized_genotype'],
                  'state': state, 'alt_dosage': gt.count(alt) if state in ('direct_compatible', 'unique_complement') else ''}
        audit.append(record)
        if state in ('direct_compatible', 'unique_complement'):
            idx.append(i)
            target.append(gt.count(alt))
    return np.array(idx), np.array(target, dtype=float), audit

def select_sites(genotypes, site_rows, maf_min=.01, window_bp=500_000, r2_max=.2):
    frequency = genotypes.mean(axis=0) / 2
    variable = np.where((frequency >= maf_min) & (frequency <= 1 - maf_min))[0]
    centered = genotypes[:, variable].astype(float) - genotypes[:, variable].mean(axis=0)
    norm = np.sqrt((centered * centered).sum(axis=0))
    correlation_columns = centered / norm
    kept_local = []
    for local, original in enumerate(variable):
        chrom, pos = site_rows[original][0], int(site_rows[original][1])
        neighbors = [j for j in kept_local if site_rows[variable[j]][0] == chrom
                     and abs(int(site_rows[variable[j]][1]) - pos) <= window_bp]
        if neighbors:
            r2 = (correlation_columns[:, neighbors].T @ correlation_columns[:, local]) ** 2
            if np.any(r2 > r2_max):
                continue
        kept_local.append(local)
    return variable[kept_local], len(variable)

def fit_reference(genotypes, target, labels):
    frequency = genotypes.mean(axis=0) / 2
    mean = 2 * frequency
    scale = np.sqrt(2 * frequency * (1 - frequency))
    z = (genotypes.astype(float) - mean) / scale
    model = PCA(n_components=4, svd_solver='full')
    scores = model.fit_transform(z)
    projection = model.transform(((target - mean) / scale).reshape(1, -1))[0]
    pc_scale = np.sqrt(model.explained_variance_)
    centers = np.stack([scores[labels == g].mean(axis=0) for g in range(5)])
    distances = np.sqrt((((projection - centers) / pc_scale) ** 2).sum(axis=1))
    return model, scores, projection, mean, scale, centers, distances

def run(input_path, ref_path, out_path):
    out = Path(out_path)
    out.mkdir(parents=True, exist_ok=True)
    qc, clean = analyze(input_path, out)
    sites, reference, labels = load_reference(ref_path)
    idx, target, audit = harmonize(clean, sites)
    if len(idx) < 100:
        raise ValueError('Insufficient overlap even for the exploratory workflow (<100)')
    g = reference[:, idx]
    rows = [sites[i] for i in idx]
    selected, after_maf = select_sites(g, rows)
    gp, t = g[:, selected], target[selected]
    if len(selected) < 100:
        raise ValueError('Insufficient markers after pruning for exploratory workflow (<100)')
    model, scores, projection, mean, scale, centers, distances = fit_reference(gp, t, labels)
    nearest = GROUPS[int(distances.argmin())]
    # Honest diagnostic: feature selection is re-fit on training individuals.
    train, test = train_test_split(np.arange(2504), test_size=.2, random_state=SEED, stratify=labels)
    held_select, _ = select_sites(g[train], rows)
    held_model, train_scores, _, held_mean, held_scale, held_centers, _ = fit_reference(
        g[train][:, held_select], target[held_select], labels[train])
    test_scores = held_model.transform((g[test][:, held_select] - held_mean) / held_scale)
    d = (((test_scores[:, None, :] - held_centers[None, :, :]) /
          np.sqrt(held_model.explained_variance_)[None, None, :]) ** 2).sum(axis=2)
    predicted = d.argmin(axis=1)
    heldout = {'train_n': len(train), 'test_n': len(test), 'selected_sites': len(held_select),
               'centroid_agreement': float((predicted == labels[test]).mean()),
               'per_group': {GROUPS[k]: {'n': int((labels[test] == k).sum()),
                                        'agreement': float((predicted[labels[test] == k] == k).mean())} for k in range(5)},
               'interpretation': 'Internal reference holdout diagnostic, not an individual accuracy probability or an admixture estimate.'}
    # Leave-one-chromosome-out refits measure only this sparse panel's sensitivity.
    sensitivity = []
    pruned_rows = [rows[i] for i in selected]
    for chrom in sorted({r[0] for r in pruned_rows}, key=int):
        keep = np.array([r[0] != chrom for r in pruned_rows])
        _, _, _, _, _, _, dd = fit_reference(gp[:, keep], t[keep], labels)
        sensitivity.append({'omitted_chromosome': chrom, 'sites': int(keep.sum()),
                            'nearest_center': GROUPS[int(dd.argmin())],
                            'distances': {GROUPS[i]: float(dd[i]) for i in range(5)}})
    ref_hashes = {name: hashlib.sha256((Path(ref_path) / name).read_bytes()).hexdigest()
                  for name in ['GRCH37.sites', 'GRCH37.sites.bin.gz', 'pca.py']}
    result = {'sample_alias': 'Sample01', 'scope': 'exploratory_low_density_reference_PCA',
              'inferred_build': 'GRCh37', 'build_status': 'database_anchors_and_coordinate_overlap_supported',
              'reference': {'source': 'peddy 0.4.8 public 1000G reference assets',
                            'reference_n': 2504, 'panel_sites': len(sites),
                            'group_counts': {GROUPS[k]: int((labels == k).sum()) for k in range(5)},
                            'hashes': ref_hashes,
                            'identity_note': 'Only package row indices and coarse labels are available; actual individual IDs are not reconstructed.'},
              'overlap_clean': len(audit), 'harmonization': dict(Counter(x['state'] for x in audit)),
              'usable_before_maf_ld': len(idx), 'after_maf': after_maf, 'used_after_ld': len(selected),
              'used_chromosomes': sorted({r[0] for r in pruned_rows}, key=int),
              'method': {'fit': 'reference individuals only; target is projected and never changes the fit',
                         'standardization': '(dosage - 2p) / sqrt(2p(1-p)); p estimated only in fitting reference',
                         'maf_min': .01, 'ld_window_bp': 500_000, 'ld_r2_max': .2,
                         'ld_algorithm': 'greedy full-reference Pearson r2 within physical window',
                         'pca_components': 4, 'svd_solver': 'full', 'seed': SEED,
                         'distance': 'Euclidean distance of PC1..PC4 after division by reference PC standard deviations',
                         'comparison': 'reference centroids; no SVM and no ancestry percentages'},
              'target_pcs': projection.tolist(), 'variance_explained': model.explained_variance_ratio_.tolist(),
              'centers': {GROUPS[k]: centers[k].tolist() for k in range(5)},
              'center_distances': {GROUPS[k]: float(distances[k]) for k in range(5)},
              'nearest_reference_center': nearest, 'heldout_diagnostic': heldout,
              'leave_one_chromosome_out': sensitivity,
              'loco_nearest_counts': dict(Counter(x['nearest_center'] for x in sensitivity)),
              'software': {'python': platform.python_version(), 'numpy': np.__version__, 'sklearn': sklearn.__version__},
              'limitations': ['Very sparse overlap; this is not a genome-wide ancestry estimate.',
                              'Public package labels are coarse groups; no 26-population fine analysis is claimed.',
                              'No supervised ADMIXTURE, NNLS proportions, ancient affinity model, or haplogroup was run.',
                              'Nearest reference center is a coordinate-space comparison, not nationality or ethnicity.',
                              'Selected-site biases and reference coverage can change the projected position.',
                              'Leave-one-chromosome-out stability is not a confidence interval or a correct-probability estimate.']}
    (out / 'reference_pca.json').write_text(json.dumps(result, indent=2, ensure_ascii=False))
    with (out / 'reference_pcs.tsv').open('w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['reference_row_index', 'group', 'PC1', 'PC2', 'PC3', 'PC4'])
        for i, row in enumerate(scores):
            writer.writerow([i, GROUPS[labels[i]], *row])
    selected_indices = set(idx[selected].tolist())
    with (out / 'harmonization_audit.tsv').open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(audit[0]) + ['used_after_ld'], delimiter='\t')
        writer.writeheader()
        for row in audit:
            writer.writerow({**row, 'used_after_ld': int(row['reference_site_index'] in selected_indices)})
    np.savez_compressed(out / 'pca_model.npz', mean=mean, scale=scale, pca_mean=model.mean_, components=model.components_,
                        eigenvalues=model.explained_variance_, reference_site_indices=idx[selected],
                        target_dosages=t, target_projection=projection)
    print(json.dumps({k: result[k] for k in ['usable_before_maf_ld', 'used_after_ld',
                     'nearest_reference_center', 'center_distances', 'heldout_diagnostic', 'loco_nearest_counts']}, indent=2), flush=True)
    return result

if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument('input')
    p.add_argument('--reference', required=True)
    p.add_argument('--out', default='results')
    args = p.parse_args()
    run(args.input, args.reference, args.out)
