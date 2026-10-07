"""Restore fixed-sample archived assembly evidence and the low-density baseline."""
import json
import hashlib
import urllib.request
import tarfile
import io
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def main():
    ref=ROOT/'references';ref.mkdir(exist_ok=True)
    # Public response snapshots distributed as data, never executed.
    for p in (ROOT/'src/public_evidence').glob('*.json'):
        (ref/p.name).write_bytes(p.read_bytes())
    request=urllib.request.Request('https://pypi.org/pypi/peddy/0.4.8/json')
    metadata=json.load(urllib.request.urlopen(request,timeout=60))
    release=next(r for r in metadata['urls'] if r['packagetype']=='sdist')
    archive=urllib.request.urlopen(release['url'],timeout=120).read()
    if hashlib.sha256(archive).hexdigest()!=release['digests']['sha256']:
        raise ValueError('peddy release archive hash mismatch')
    with tarfile.open(fileobj=io.BytesIO(archive),mode='r:gz') as tf:
        names=['peddy-0.4.8/LICENSE','peddy-0.4.8/peddy/GRCH37.sites','peddy-0.4.8/peddy/GRCH38.sites']
        for name in names:
            destination=ref/name;destination.parent.mkdir(parents=True,exist_ok=True)
            destination.write_bytes(tf.extractfile(name).read())
    print('Restored public assembly evidence and baseline coordinate panels')
if __name__=='__main__':main()
