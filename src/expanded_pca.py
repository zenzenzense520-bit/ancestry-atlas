"""Reference-only PCA with separate global and East Asian feature selection."""
import argparse
import csv
import json
import platform
from collections import Counter, deque
from pathlib import Path
import numpy as np
import sklearn
from sklearn.decomposition import PCA
from sklearn.model_selection import train_test_split

ROOT = Path(__file__).resolve().parents[1]
SEED = 20261006
GROUPS = {'global': ['AFR', 'AMR', 'EAS', 'EUR', 'SAS'],
          'east_asian': ['CHB', 'CHS', 'JPT', 'CDX', 'KHV']}

def select_sites(genotypes, sites, maf_min=.01, window_bp=500000, r2_max=.2):
    # genotype matrix: reference individuals x overlap variants
    complete = np.all(genotypes >= 0, axis=0)
    p = genotypes.mean(axis=0, dtype=np.float32) / 2
    variable = np.flatnonzero(complete & (p >= maf_min) & (p <= 1 - maf_min))
    centered = genotypes[:, variable].T.astype(np.float32)
    centered -= centered.mean(axis=1, keepdims=True)
    centered /= np.sqrt(np.sum(centered * centered, axis=1))[:, None]
    kept, neighbors = [], deque()
    last_chrom = None
    for local, original in enumerate(variable):
        chrom, pos = sites[original]['chromosome'], int(sites[original]['position'])
        if chrom != last_chrom:
            neighbors.clear()
            last_chrom = chrom
        while neighbors and pos - int(sites[variable[neighbors[0]]]['position']) > window_bp:
            neighbors.popleft()
        if neighbors:
            r2 = (centered[list(neighbors)] @ centered[local]) ** 2
            if np.any(r2 > r2_max):
                continue
        kept.append(local)
        neighbors.append(local)
    return variable[kept], len(variable)

def fit_reference(genotypes, target, labels, groups):
    p = genotypes.mean(axis=0, dtype=np.float32) / 2
    mean, scale = 2 * p, np.sqrt(2 * p * (1 - p))
    z = (genotypes.astype(np.float32) - mean) / scale
    solver = 'full' if len(genotypes) <= 600 else 'randomized'
    model = PCA(n_components=4, svd_solver=solver, random_state=SEED,
                n_oversamples=20, iterated_power=7)
    scores = model.fit_transform(z)
    projection = model.transform(((target - mean) / scale).reshape(1, -1))[0]
    centers = np.stack([scores[labels == group].mean(axis=0) for group in groups])
    distances = np.sqrt(np.sum(((projection - centers) / np.sqrt(model.explained_variance_)) ** 2, axis=1))
    return model, scores, projection, mean, scale, centers, distances

def run(scope):
    out = ROOT / 'results'
    data = np.load(out / 'reference_overlap.npz')
    samples = list(csv.DictReader((out / 'reference_samples.tsv').open(), delimiter='\t'))
    sites = list(csv.DictReader((out / 'reference_overlap_sites.tsv').open(), delimiter='\t'))
    groups = GROUPS[scope]
    sample_indices = np.array([i for i, row in enumerate(samples) if scope == 'global' or row['group'] == 'EAS'])
    labels = np.array([samples[i]['group' if scope == 'global' else 'population'] for i in sample_indices])
    g, target = data['genotypes'][:, sample_indices].T, data['target_dosages']
    maf = .01 if scope == 'global' else .05
    selected, after_maf = select_sites(g, sites, maf_min=maf)
    if len(selected) < 100:
        raise ValueError('Fewer than 100 markers after filtering')
    model, scores, projection, mean, scale, centers, distances = fit_reference(g[:, selected], target[selected], labels, groups)
    print(f'{scope}: {len(selected)} selected sites, {len(g)} reference individuals; main fit complete', flush=True)
    train, test = train_test_split(np.arange(len(g)), test_size=.2, random_state=SEED, stratify=labels)
    hold_selected, _ = select_sites(g[train], sites, maf_min=maf)
    hmodel, hscores, htarget, hmean, hscale, hcenters, hd = fit_reference(g[train][:, hold_selected], target[hold_selected], labels[train], groups)
    test_scores = hmodel.transform((g[test][:, hold_selected].astype(np.float32) - hmean) / hscale)
    dd = np.sum(((test_scores[:, None, :] - hcenters[None, :, :]) / np.sqrt(hmodel.explained_variance_)) ** 2, axis=2)
    predicted = np.array(groups)[dd.argmin(axis=1)]
    holdout = {'train_n': len(train), 'test_n': len(test), 'selected_sites': len(hold_selected),
               'centroid_agreement': float(np.mean(predicted == labels[test])),
               'target_nearest_center': groups[int(hd.argmin())],
               'per_group': {group: {'n': int(np.sum(labels[test] == group)),
                                    'agreement': float(np.mean(predicted[labels[test] == group] == group))} for group in groups},
               'interpretation': 'Internal reference agreement; not a personal accuracy probability or an ancestry proportion.'}
    sensitivity = []
    chromosomes = np.array([r['chromosome'] for r in sites])[selected]
    for chrom in sorted(set(chromosomes), key=int):
        keep = chromosomes != chrom
        *_, distances_loco = fit_reference(g[:, selected[keep]], target[selected[keep]], labels, groups)
        sensitivity.append({'omitted_chromosome': str(chrom), 'sites': int(keep.sum()),
                            'nearest_center': groups[int(distances_loco.argmin())],
                            'distances': dict(zip(groups, map(float, distances_loco)))})
        print(f'{scope}: LOCO {chrom}/22 complete', flush=True)
    pc2 = np.sqrt(np.sum(((projection[:2] - centers[:, :2]) / np.sqrt(model.explained_variance_[:2])) ** 2, axis=1))
    result = {'sample_alias': 'Sample01', 'scope': scope, 'version': '0.5',
              'reference': {'source': 'PLINK public 1000 Genomes Phase 3 GRCh37', 'reference_n': len(g),
                            'group_counts': dict(Counter(labels)), 'actual_ids_and_official_labels_verified': True},
              'usable_before_maf_ld': len(sites), 'after_maf': after_maf, 'used_after_ld': len(selected),
              'used_chromosomes': sorted(set(chromosomes), key=int),
              'method': {'fit': 'Reference only; target never affects selection, allele frequencies, LD or PCA',
                         'standardization': '(g - 2p) / sqrt(2p(1-p))', 'maf_min': maf,
                         'ld_window_bp': 500000, 'ld_r2_max': .2,
                         'ld_algorithm': 'Greedy reference Pearson r2, ordered by PVAR index',
                         'pca_components': 4, 'svd_solver': model.svd_solver, 'seed': SEED,
                         'missing_data': 'Only complete-call reference variants; no imputation',
                         'loco': 'Refit after removing one chromosome from the main selected set; selection not repeated',
                         'holdout': 'Stratified 80/20; feature selection refitted on training references only'},
              'target_pcs': projection.tolist(), 'variance_explained': model.explained_variance_ratio_.tolist(),
              'centers': dict(zip(groups, centers.tolist())), 'center_distances': dict(zip(groups, map(float, distances))),
              'nearest_reference_center': groups[int(distances.argmin())],
              'two_pc_nearest_center': groups[int(pc2.argmin())],
              'two_pc_center_distances': dict(zip(groups, map(float, pc2))),
              'heldout_diagnostic': holdout, 'leave_one_chromosome_out': sensitivity,
              'loco_nearest_counts': dict(Counter(r['nearest_center'] for r in sensitivity)),
              'software': {'python': platform.python_version(), 'numpy': np.__version__, 'sklearn': sklearn.__version__},
              'limitations': ['Distances are not ancestry percentages or ethnicity/nationality labels.',
                              'Reference ascertainment and chip coverage can affect position.',
                              'Small PCs and adjacent population centers can be sensitive to model choices.',
                              'LOCO stability and internal agreement are not confidence intervals or personal accuracy.',
                              'No ADMIXTURE, ancient affinity, Y or mitochondrial haplogroup inference was run.']}
    # Convert numpy str keys to JSON-native strings.
    result['reference']['group_counts'] = {str(k): int(v) for k, v in result['reference']['group_counts'].items()}
    result['used_chromosomes'] = [str(c) for c in result['used_chromosomes']]
    (out / f'{scope}_pca.json').write_text(json.dumps(result, indent=2, ensure_ascii=False))
    with (out / f'{scope}_reference_pcs.tsv').open('w', newline='') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['reference_row_index', 'sample_id', 'group', 'PC1', 'PC2', 'PC3', 'PC4'])
        for local, row in enumerate(scores):
            w.writerow([sample_indices[local], samples[sample_indices[local]]['sample_id'], labels[local], *row])
    np.savez_compressed(out / f'{scope}_model.npz', mean=mean, scale=scale, pca_mean=model.mean_,
        components=model.components_, eigenvalues=model.explained_variance_, selected_overlap_indices=selected,
        target_dosages=target[selected], target_projection=projection, reference_sample_indices=sample_indices)
    print(json.dumps({k: result[k] for k in ['used_after_ld', 'nearest_reference_center', 'center_distances',
                                           'heldout_diagnostic', 'loco_nearest_counts']}, indent=2), flush=True)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--scope', choices=['all', *GROUPS], default='all')
    args = parser.parse_args()
    for name in GROUPS if args.scope == 'all' else [args.scope]:
        run(name)
