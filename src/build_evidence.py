"""Rebuild assembly evidence from archived public responses and reference sites."""
import csv
import json
from pathlib import Path
from qc import parse

ROOT = Path(__file__).resolve().parents[1]

def main(input_path):
    ref = ROOT / 'references'
    anchors = json.loads((ref / 'build_anchors.json').read_text())
    # Ensembl coordinates are one based. NCBI SPDI positions are zero based.
    anchors.insert(0, {'id': 'rs2710890', 'chromosome': '1', 'observed_position': 958905,
                      'grch37': 958905, 'grch38': 1023525, 'retrieved': True})
    ncbi = json.loads((ref / 'ncbi_rs3748595.json').read_text())
    positions = {}
    for place in ncbi['primary_snapshot_data']['placements_with_allele']:
        for assembly in place.get('placement_annot', {}).get('seq_id_traits_by_assembly', []):
            name = assembly['assembly_name']
            if place['seq_id'].startswith('NC_000001.') and name.startswith(('GRCh37', 'GRCh38')):
                allele = place['alleles'][0]['allele']['spdi']
                positions[name.split('.')[0].lower()] = allele['position'] + 1
    assert positions == {'grch37': 887560, 'grch38': 952180}, positions
    anchors.insert(0, {'id': 'rs3748595', 'chromosome': '1', 'observed_position': 887560,
                      **positions, 'retrieved': True})
    raw, raw_groups = parse(input_path)
    for anchor in anchors:
        observed = {r['position'] for r in raw if r['gxid'].split('.')[0] == anchor['id'] and r['chromosome'] == anchor['chromosome']}
        assert observed == {anchor['observed_position']}, (anchor['id'], observed)
        if anchor['id'] != 'rs3748595':
            for assembly, key in [('grch37','grch37'),('grch38','grch38')]:
                response = json.loads((ref / f'{assembly}_{anchor["id"]}.json').read_text())
                matched = [m for m in response['mappings'] if m['seq_region_name'] == anchor['chromosome']
                           and m['assembly_name'].lower() == assembly]
                assert any(m['start'] == anchor[key] for m in matched), (anchor['id'], assembly)
    clean = {(r['chromosome'], int(r['position'])) for r in csv.DictReader(
        (ROOT / 'results/normalized_acgt.tsv').open(), delimiter='\t')}
    overlap = {}
    raw_overlap = {}
    for build in ['GRCH37', 'GRCH38']:
        coords = {(r[0], int(r[1])) for r in [s.split(':') for s in
                  (ref / 'peddy-0.4.8/peddy' / f'{build}.sites').read_text().splitlines()]}
        overlap[build] = len(clean & coords)
        raw_overlap[build] = len(set(raw_groups) & coords)
    evidence = {'inferred_build': 'GRCh37', 'reported_build': None,
                'status': 'strongly_supported_inference_not_provider_confirmation',
                'anchors': anchors, 'anchors_matching_grch37': sum(a['observed_position'] == a['grch37'] for a in anchors),
                'anchors_matching_grch38': sum(a['observed_position'] == a['grch38'] for a in anchors),
                'clean_coordinate_overlap': overlap, 'raw_coordinate_overlap': raw_overlap, 'retrieved_date': '2026-10-05',
                'coordinate_note': 'One-based export/Ensembl; NCBI SPDI positions converted from zero based.',
                'limitations': ['Anchor agreement supports assembly, not chip provenance or strand metadata.',
                                'No coordinate lift-over was applied.']}
    (ROOT / 'results/build_evidence.json').write_text(json.dumps(evidence, ensure_ascii=False, indent=2))
    print(json.dumps({k: evidence[k] for k in ['anchors_matching_grch37','anchors_matching_grch38','clean_coordinate_overlap']}))

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('input')
    main(parser.parse_args().input)
