"""Compare exact and randomized East Asian PCs on the identical selected matrix."""
import json
import csv
from pathlib import Path
import numpy as np
from sklearn.decomposition import PCA
ROOT=Path(__file__).resolve().parents[1]
def main():
    with np.load(ROOT/'results/reference_overlap.npz') as data:
        genotypes=data['genotypes'];target=data['target_dosages']
    with np.load(ROOT/'results/east_asian_model.npz') as saved:
        model={k:saved[k] for k in saved.files}
    with (ROOT/'results/reference_samples.tsv').open() as f:samples=list(csv.DictReader(f,delimiter='\t'))
    indices=model['selected_overlap_indices'];refs=model['reference_sample_indices']
    labels=np.array([samples[i]['population'] for i in refs]);groups=['CHB','CHS','JPT','CDX','KHV']
    g=genotypes[indices][:,refs].T.astype(np.float32)
    z=(g-model['mean'])/model['scale'];t=(target[indices]-model['mean'])/model['scale']
    exact=json.loads((ROOT/'results/east_asian_pca.json').read_text());pca=PCA(n_components=4,svd_solver='randomized',random_state=20261006,n_oversamples=20,iterated_power=7)
    scores=pca.fit_transform(z);projection=pca.transform(t.reshape(1,-1))[0]
    centers=np.stack([scores[labels==g].mean(axis=0) for g in groups]);distances=np.linalg.norm((projection-centers)/np.sqrt(pca.explained_variance_),axis=1)
    result={'scope':'east_asian','sites':len(indices),'same_selected_variants':True,'reference_only':True,
        'full':{'nearest':exact['nearest_reference_center'],'distances':exact['center_distances']},
        'randomized':{'nearest':groups[int(distances.argmin())],'distances':dict(zip(groups,map(float,distances))),
            'target_pcs':projection.tolist(),'seed':20261006,'n_oversamples':20,'iterated_power':7},
        'interpretation':'Small PC solver sensitivity limits adjacent-centroid identity interpretations.'}
    (ROOT/'qa/solver_sensitivity.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
if __name__=='__main__':main()
