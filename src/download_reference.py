"""Download pinned public GRCh37 Phase 3 resources with byte/hash validation."""
import concurrent.futures
import hashlib
import json
import time
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RESOURCES = [
    ('all_phase3.pgen.zst', 'https://www.dropbox.com/s/y6ytfoybz48dc0u/all_phase3.pgen.zst?dl=1', 2414715690, '1feb503994db9d7e4f85bdcd2dd7bdd58fccb08b15bdaa2d6b4cf2256773f71e'),
    ('all_phase3.pvar.zst', 'https://www.dropbox.com/s/c95n8quqwqww4s0/all_phase3_noannot.pvar.zst?dl=1', 643457044, '82b3462007ae44e0293500d0b71fa24c4d8cc4231139993af5bbb82037a5f8cc'),
    ('phase3_orig.psam', 'https://www.dropbox.com/s/c8z9ws4xwx046kv/phase3_orig.psam?dl=1', 45101, '3384a707515ed9aa4c005f186daf7c04e517dba1f6450219d881f8e327394161'),
    ('integrated_call_samples_v3.20130502.ALL.panel', 'https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel', 55156, 'b4023dc6ee2d62ee89c8d4d347db4d348e65518d66d346574cdae7a4bbd76858'),
]

def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda: f.read(8 << 20), b''):
            h.update(chunk)
    return h.hexdigest()

def download(resource):
    name, url, size, sha = resource
    out = ROOT / 'references' / name
    out.parent.mkdir(exist_ok=True)
    part = out.with_suffix(out.suffix + '.part')
    if not (out.exists() and out.stat().st_size == size and digest(out) == sha):
        for attempt in range(6):
            offset = part.stat().st_size if part.exists() else 0
            try:
                headers = {'User-Agent': 'AncestryAtlas/0.5 public-reference-download'}
                if offset:
                    headers['Range'] = f'bytes={offset}-'
                with urllib.request.urlopen(urllib.request.Request(url, headers=headers), timeout=120) as response:
                    if offset and response.status != 206:
                        offset = 0
                    with part.open('ab' if offset else 'wb') as f:
                        for chunk in iter(lambda: response.read(4 << 20), b''):
                            f.write(chunk)
                if part.stat().st_size != size:
                    raise ValueError(f'{name}: size {part.stat().st_size} != {size}')
                if digest(part) != sha:
                    raise ValueError(f'{name}: SHA256 mismatch')
                part.replace(out)
                break
            except Exception as error:
                print(f'{name}: attempt {attempt + 1}: {error}', flush=True)
                if attempt == 5:
                    raise
                time.sleep(2)
    meta = {'filename': name, 'url': url, 'bytes': size, 'sha256': sha,
            'source': 'PLINK 2 resources / IGSR public Phase 3 GRCh37'}
    out.with_suffix(out.suffix + '.download.json').write_text(json.dumps(meta, indent=2))
    print(f'Verified {name}: {size:,} bytes', flush=True)
    return meta

def main():
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
        manifest = list(pool.map(download, RESOURCES))
    (ROOT / 'references/download_manifest.json').write_text(json.dumps(manifest, indent=2))

if __name__ == '__main__':
    main()
