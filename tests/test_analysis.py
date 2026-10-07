"""Independent numerical checks of persisted models and QC edge cases."""
import csv
import json
import sys
import unittest
from collections import deque
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from qc import consensus, normalize, kind

def read_tsv(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))

def rows(*calls):
    return [{'normalized_genotype':normalize(call),'call_kind':kind(call)} for call in calls]

class QCEdgeCases(unittest.TestCase):
    def test_unordered_calls_agree(self):
        self.assertEqual(consensus(rows('TC','CT')),('acgt_consensus','CT'))
    def test_missing_does_not_replace_observed_call(self):
        self.assertEqual(consensus(rows('--','GG')),('acgt_consensus','GG'))
    def test_mixed_indel_base_group_is_not_forced_to_snp(self):
        self.assertEqual(consensus(rows('II','AA')),('mixed_variant_encoding',None))
    def test_incompatible_calls_excluded(self):
        self.assertEqual(consensus(rows('AA','AG')),('acgt_discordant',None))
    def test_indel_order_normalized_but_different_calls_excluded(self):
        self.assertEqual(consensus(rows('ID','DI')),('indel_only','DI'))
        self.assertEqual(consensus(rows('ID','II')),('indel_code_discordant',None))

@unittest.skipUnless((ROOT/'results/global_model.npz').exists() and (ROOT/'results/east_asian_model.npz').exists(), 'Private model artifacts are required for numerical checks')
class ArtifactChecks(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        with np.load(ROOT/'results/reference_overlap.npz') as archive:
            cls.data={key:archive[key] for key in archive.files}
        cls.samples=read_tsv(ROOT/'results/reference_samples.tsv')
        cls.sites=read_tsv(ROOT/'results/reference_overlap_sites.tsv')
        cls.scopes=['global','east_asian'];cls.metrics={}
    @classmethod
    def tearDownClass(cls):
        (ROOT/'qa/analysis_checks.json').write_text(json.dumps(cls.metrics,indent=2))
    def artifacts(self,scope):
        return np.load(ROOT/f'results/{scope}_model.npz'),json.loads((ROOT/f'results/{scope}_pca.json').read_text())
    def test_raw_coordinate_and_alignment_accounting(self):
        q=json.loads((ROOT/'results/sample_qc.json').read_text());r=json.loads((ROOT/'results/expanded_reference_qc.json').read_text())
        self.assertEqual(sum(q['raw_call_kinds'].values()),q['raw_rows'])
        self.assertEqual(sum(q['coordinate_states'].values()),q['unique_coordinates'])
        self.assertEqual(q['raw_rows']-q['unique_coordinates'],q['redundant_coordinate_rows'])
        self.assertEqual(sum(r['harmonization'].values()),q['autosomal_clean_acgt_coordinates'])
        self.assertEqual(len(self.sites),r['usable_before_maf_ld'])
        self.assertEqual(int((self.data['genotypes']<0).sum()),r['reference_missing_cells'])
        b=json.loads((ROOT/'results/build_evidence.json').read_text())
        self.assertEqual(b['anchors_matching_grch37'],len(b['anchors']))
        self.assertEqual(b['anchors_matching_grch38'],0)
    def test_official_ids_and_labels_in_saved_reference_order(self):
        panel={r['sample']:r for r in read_tsv(ROOT/'references/integrated_call_samples_v3.20130502.ALL.panel')}
        self.assertEqual(len(set(r['sample_id'] for r in self.samples)),len(panel))
        for i,row in enumerate(self.samples):
            self.assertEqual(int(row['reference_row_index']),i)
            self.assertEqual(row['group'],panel[row['sample_id']]['super_pop'])
            self.assertEqual(row['population'],panel[row['sample_id']]['pop'])
        for scope in self.scopes:
            m,p=self.artifacts(scope)
            saved=read_tsv(ROOT/f'results/{scope}_reference_pcs.tsv')
            for r,i in zip(saved,m['reference_sample_indices']):self.assertEqual(r['sample_id'],self.samples[i]['sample_id'])
            self.assertEqual(len(saved),p['reference']['reference_n'])
    def test_all_target_alleles_trace_to_source(self):
        clean={(r['chromosome'],r['position']):r for r in read_tsv(ROOT/'results/normalized_acgt.tsv')}
        translation=str.maketrans('ACGT','TGCA')
        for i,r in enumerate(self.sites):
            call=clean[r['chromosome'],r['position']]['normalized_genotype']
            self.assertNotIn({r['ref'],r['alt']},[set('AT'),set('CG')])
            if r['state']=='unique_complement':call=call.translate(translation)
            self.assertTrue(set(call)<={r['ref'],r['alt']})
            self.assertEqual(call.count(r['alt']),self.data['target_dosages'][i])
            self.assertEqual(int(r['reference_variant_index']),self.data['reference_variant_indices'][i])
    def test_projection_scores_centers_and_distances_independently_reconstructed(self):
        for scope in self.scopes:
            m,p=self.artifacts(scope)
            projection=((m['target_dosages'].astype(float)-m['mean'])/m['scale']-m['pca_mean'])@m['components'].T
            error=float(np.max(abs(projection-m['target_projection'])))
            np.testing.assert_allclose(projection,m['target_projection'],atol=2e-3)
            np.testing.assert_allclose(projection,p['target_pcs'],atol=2e-3)
            g=self.data['genotypes'][m['selected_overlap_indices']][:,m['reference_sample_indices']].T.astype(float)
            z=(g-m['mean'])/m['scale']-m['pca_mean'];scores=z@m['components'].T
            saved=read_tsv(ROOT/f'results/{scope}_reference_pcs.tsv')
            expected=np.array([[float(r[f'PC{i}']) for i in range(1,5)] for r in saved]);np.testing.assert_allclose(scores,expected,atol=2e-3)
            labels=np.array([r['group'] for r in saved])
            for group,center in p['centers'].items():
                reconstructed=scores[labels==group].mean(axis=0);np.testing.assert_allclose(reconstructed,center,atol=2e-3)
                distance=np.linalg.norm((projection-reconstructed)/np.sqrt(m['eigenvalues']))
                self.assertAlmostEqual(distance,p['center_distances'][group],places=4)
            self.metrics[scope]={'max_float64_projection_error':error}
    def test_reference_frequencies_pruning_and_eigenvectors(self):
        for scope in self.scopes:
            m,p=self.artifacts(scope);idx=m['selected_overlap_indices']
            g=self.data['genotypes'][idx][:,m['reference_sample_indices']].astype(float)
            freq=g.mean(axis=1)/2;self.assertTrue(np.all((freq>=p['method']['maf_min'])&(freq<=1-p['method']['maf_min'])))
            np.testing.assert_allclose(m['mean'],2*freq,atol=3*np.finfo(m['mean'].dtype).eps)
            np.testing.assert_allclose(m['scale'],np.sqrt(2*freq*(1-freq)),atol=3*np.finfo(m['scale'].dtype).eps)
            normalized=g-g.mean(axis=1,keepdims=True);normalized/=np.linalg.norm(normalized,axis=1,keepdims=True)
            active=deque();last=None;max_r2=0.
            for local,original in enumerate(idx):
                site=self.sites[original];chrom,pos=site['chromosome'],int(site['position'])
                if chrom!=last:active.clear();last=chrom
                while active and pos-int(self.sites[idx[active[0]]]['position'])>p['method']['ld_window_bp']:active.popleft()
                if active:max_r2=max(max_r2,float(np.max((normalized[list(active)]@normalized[local])**2)))
                active.append(local)
            self.assertLessEqual(max_r2,.2+1e-6)
            z=(g.T-m['mean'])/m['scale']-m['pca_mean'];u=z@m['components'].T
            residual=np.linalg.norm(z.T@u/(len(u)-1)-m['components'].T*m['eigenvalues'],axis=0)/m['eigenvalues']
            self.assertTrue(np.all(residual<1e-3))
            self.metrics.setdefault(scope,{}).update({'max_retained_window_r2':max_r2,'eigenvector_relative_residual':residual.tolist()})
    def test_loco_and_holdout_diagnostic_accounting(self):
        for scope in self.scopes:
            _,p=self.artifacts(scope)
            self.assertEqual({r['omitted_chromosome'] for r in p['leave_one_chromosome_out']},set(map(str,range(1,23))))
            self.assertEqual(sum(p['loco_nearest_counts'].values()),22)
            h=p['heldout_diagnostic'];self.assertEqual(h['train_n']+h['test_n'],p['reference']['reference_n'])
            self.assertEqual(sum(r['n'] for r in h['per_group'].values()),h['test_n'])
            self.assertAlmostEqual(sum(r['n']*r['agreement'] for r in h['per_group'].values())/h['test_n'],h['centroid_agreement'])

if __name__=='__main__':unittest.main()
