"""Package private results and a separate source-only GitHub deliverable."""
import hashlib
import json
import shutil
from pathlib import Path
import zipfile
ROOT=Path(__file__).resolve().parents[1]

def tree(directory):
    return [p for p in directory.rglob('*') if p.is_file() and '__pycache__' not in p.parts]
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def main():
    out=ROOT.parent/'deliverables';out.mkdir(exist_ok=True)
    source=[]
    for d in ['src','tests','docs','tools','.github']:source+=tree(ROOT/d)
    source += [ROOT/name for name in ['README.md','LICENSE','requirements.txt','.gitignore','Publish-GitHub.cmd']]
    publication=ROOT.parent/'ancestry_atlas_github';publication.mkdir(exist_ok=True)
    for p in source:
        destination=publication/p.relative_to(ROOT);destination.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,destination)
    source_zip=out/'AncestryAtlas_source_v05.zip'
    with zipfile.ZipFile(source_zip,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as archive:
        for p in sorted(source):archive.write(p,Path('ancestry-atlas')/p.relative_to(ROOT))
    files=source.copy()
    files += [p for p in tree(ROOT/'results') if p.name not in ['reference_pca.json','reference_pcs.tsv','harmonization_audit.tsv','pca_model.npz']]
    files += list((ROOT/'figures').glob('global_*'))+list((ROOT/'figures').glob('east_asian_*'))
    files += [ROOT/'report/Sample01_reference_report_v05.html',ROOT/'report/Sample01_reference_findings_v05.md']
    files += [p for p in tree(ROOT/'qa') if p.name in ['analysis_checks.json','browser_checks.json','design_review.md','solver_sensitivity.json','plink_genotype_crosscheck.json','plink_binary.json','plink_crosscheck.log','plink_crosscheck.traw','plink_selection.txt'] or p.name.startswith(('hero-','east-','report-'))]
    files += list((ROOT/'references').glob('*.json'))
    files += [ROOT/'references'/n for n in ['phase3_orig.psam','integrated_call_samples_v3.20130502.ALL.panel']]
    panel=ROOT/'references/peddy-0.4.8';files += [panel/'LICENSE',panel/'peddy/GRCH37.sites',panel/'peddy/GRCH38.sites']
    font=ROOT/'references/fontsource/package';files += [font/'LICENSE',font/'package.json']+list((font/'files').glob('*.woff2'))
    files=sorted(set(files));manifest={'version':'0.5','sample_alias':'Sample01','contains_personal_genotypes':True,'original_input_included':False,
       'large_public_reference_files_included':False,'files':{str(p.relative_to(ROOT)):{'bytes':p.stat().st_size,'sha256':sha(p)} for p in files}}
    (ROOT/'MANIFEST.json').write_text(json.dumps(manifest,indent=2));files.append(ROOT/'MANIFEST.json')
    private_zip=out/'Sample01_research_v05.zip'
    with zipfile.ZipFile(private_zip,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as archive:
        for p in files:archive.write(p,Path('Sample01_research_v05')/p.relative_to(ROOT))
    for package in [source_zip,private_zip]:
        with zipfile.ZipFile(package) as archive:
            if archive.testzip() is not None:raise ValueError('ZIP integrity error')
    html=out/'Sample01_reference_report_v05.html';shutil.copyfile(ROOT/'report/Sample01_reference_report_v05.html',html)
    print(json.dumps({'html_bytes':html.stat().st_size,'research_zip_bytes':private_zip.stat().st_size,'source_zip_bytes':source_zip.stat().st_size,'private_files':len(files),'source_files':len(source)},indent=2))
if __name__=='__main__':main()
