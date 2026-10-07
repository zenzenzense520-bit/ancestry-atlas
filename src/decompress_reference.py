"""Full stream decompression; never accept a truncated output as complete."""
import json
from pathlib import Path
import zstandard
from download_reference import digest

ROOT = Path(__file__).resolve().parents[1]
SIZES = {'all_phase3.pgen': 6696871489, 'all_phase3.pvar': 2318137353}

def main():
    manifest = []
    for name, expected in SIZES.items():
        target = ROOT / 'references' / name
        if not (target.exists() and target.stat().st_size == expected):
            temp = target.with_suffix(target.suffix + '.decompressing')
            count = 0
            with target.with_suffix(target.suffix + '.zst').open('rb') as source:
                with zstandard.ZstdDecompressor().stream_reader(source, read_across_frames=True) as stream:
                    with temp.open('wb') as output:
                        for block in iter(lambda: stream.read(8 << 20), b''):
                            count += len(block)
                            output.write(block)
            if count != expected:
                raise ValueError(f'{name}: decompressed {count}, expected {expected}')
            temp.replace(target)
        manifest.append({'filename': name, 'bytes': expected, 'sha256': digest(target)})
        print(f'Decompressed and verified {name}', flush=True)
    (ROOT / 'references/decompression_manifest.json').write_text(json.dumps(manifest, indent=2))

if __name__ == '__main__':
    main()
