"""Reject accidental publication of private inputs, outputs and runtime files."""
import subprocess
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
DENIED={'input','upload','results','report','figures','qa','references','runtime','deliverables','analysis-runtime','.venv','venv','__pycache__'}
EXTENSIONS={'.npz','.npy','.pgen','.pvar','.zip','.zst'}
def allowed(path):
    p=Path(path)
    return not (set(p.parts)&DENIED or p.suffix in EXTENSIONS or any(x=='.env' or x.startswith('.env.') for x in p.parts))
def main():
    proc=subprocess.run(['git','ls-files','-z'],cwd=ROOT,capture_output=True)
    files=proc.stdout.decode().strip('\0').split('\0') if proc.returncode==0 else [str(p.relative_to(ROOT)) for d in ['src','tests','docs','tools'] for p in (ROOT/d).rglob('*') if p.is_file() and '__pycache__' not in p.parts]
    invalid=[p for p in files if p and not allowed(p)]
    if invalid:raise SystemExit('Private files in publication scope: '+', '.join(invalid))
    print(f'Publication scope checked: {len(files)} tracked/source files')
if __name__=='__main__':main()
