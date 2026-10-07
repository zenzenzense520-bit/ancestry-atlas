"""Independent PLINK CLI check against the Python-extracted reference matrix."""
import argparse
import csv
import json
import subprocess
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
def main(binary):
    with (ROOT/'results/reference_overlap_sites.tsv').open() as f:sites=list(csv.DictReader(f,delimiter='\t'))
    data=np.load(ROOT/'results/reference_overlap.npz')
    rng=np.random.default_rng(20261006);indices=set(rng.choice(len(sites),100,replace=False).tolist())
    for chrom in map(str,range(1,23)):
        idx=[i for i,r in enumerate(sites) if r['chromosome']==chrom]
        if idx:indices.update([idx[0],idx[-1]])
    selected=sorted(indices);qa=ROOT/'qa';selection=qa/'plink_selection.txt'
    selection.write_text('\n'.join(sites[i]['reference_id'] for i in selected)+'\n')
    ref=ROOT/'references';output=qa/'plink_crosscheck'
    subprocess.run([str(Path(binary).resolve()),'--pgen',str(ref/'all_phase3.pgen'),'--pvar',str(ref/'all_phase3.pvar'),
        '--psam',str(ref/'phase3_orig.psam'),'--extract',str(selection),'--export','A-transpose',
        '--threads','2','--memory','6900','--out',str(output)],check=True,stdout=subprocess.DEVNULL)
    lookup={sites[i]['reference_id']:i for i in selected};checked=0
    with output.with_suffix('.traw').open() as f:
        reader=csv.reader(f,delimiter='\t');header=next(reader)
        with (ROOT/'results/reference_samples.tsv').open() as sf:samples=list(csv.DictReader(sf,delimiter='\t'))
        ids=[r['sample_id'] for r in samples]
        # PLINK sample headers may include a zero FID separated by an underscore.
        if [s.removeprefix('0_') for s in header[6:]]!=ids:raise ValueError('PLINK sample order mismatch')
        for row in reader:
            i=lookup[row[1]];dosage=np.array(row[6:],dtype=np.int8)
            if row[4]==sites[i]['ref']:dosage=2-dosage
            elif row[4]!=sites[i]['alt']:raise ValueError('Counted allele mismatch')
            np.testing.assert_array_equal(dosage,data['genotypes'][i]);checked+=len(dosage)
    result={'selected_variants':len(selected),'reference_n':len(ids),'checked_genotype_cells':checked,
            'all_match':True,'selection':'100 deterministic random variants plus first/last eligible SNP of each autosome',
            'plink_version':subprocess.check_output([str(Path(binary).resolve()),'--version'],text=True).strip()}
    (qa/'plink_genotype_crosscheck.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--binary',required=True);main(parser.parse_args().binary)
