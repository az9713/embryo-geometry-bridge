"""Verify the16 publicly linked supplementary assets against official metadata."""
from pathlib import Path
from urllib.parse import urlsplit, parse_qs
import hashlib
import json
import re
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT/'evidence/public_supplement_asset_manifest.json'


def main():
    if OUT.exists():
        raise RuntimeError('Preserve asset manifest')
    xml_file = ROOT/'data/raw/primary.xml'
    metadata_file = ROOT/'data/raw/PMC12361397.1.json'
    listing_file = ROOT/'data/raw/oa_object_listing.xml'
    metadata = json.loads(metadata_file.read_text(encoding='utf-8'))
    if metadata['pmcid'] != 'PMC12361397' or metadata['doi'] != '10.1038/s41467-025-62483-6':
        raise RuntimeError('Source identity differs')
    urls = {urlsplit(u).path.rsplit('/', 1)[-1]: u for u in metadata['media_urls']}
    objects = {}
    for e in ET.parse(listing_file).getroot():
        d = {c.tag.rsplit('}', 1)[-1]: c.text for c in e}
        if 'Key' in d:
            objects[d['Key'].rsplit('/', 1)[-1]] = d
    assets = []
    for e in ET.parse(xml_file).findall('.//supplementary-material'):
        if not re.fullmatch(r'MOESM\d+', e.get('id', '')):
            continue
        media = e.find('media')
        name = media.get('{http://www.w3.org/1999/xlink}href')
        p = ROOT/'data/raw'/name
        data = p.read_bytes()
        expected = parse_qs(urlsplit(urls[name]).query)['md5'][0]
        if len(data) != int(objects[name]['Size']) or hashlib.md5(data).hexdigest() != expected:
            raise RuntimeError('Official size/MD5 mismatch: '+name)
        ext = p.suffix
        magic = (ext == '.pdf' and data.startswith(b'%PDF-')) or (ext == '.xlsx' and data.startswith(b'PK\x03\x04')) or (ext == '.avi' and data[:4] == b'RIFF' and data[8:12] == b'AVI ')
        if not magic:
            raise RuntimeError('File signature mismatch: '+name)
        assets.append({'source_id': e.get('id'), 'filename': name,
            'description': ' '.join(''.join(media.itertext()).split()),
            'file': p.relative_to(ROOT).as_posix(), 'bytes': len(data),
            'sha256': hashlib.sha256(data).hexdigest(), 'official_MD5': expected,
            'official_download_URL': 'https://pmc-oa-opendata.s3.amazonaws.com/PMC12361397.1/'+name,
            'content_signature_verified': True})
    if len(assets) != 16 or len({a['filename'] for a in assets}) != 16:
        raise RuntimeError('Exactly16 unique linked supplement files required')
    types = {suffix: sum(a['filename'].endswith(suffix) for a in assets) for suffix in ('.pdf', '.xlsx', '.avi')}
    if types != {'.pdf': 4, '.xlsx': 2, '.avi': 10}:
        raise RuntimeError('Supplement coverage differs')
    result = {'status': 'PASS_ALL16_LINKED_SUPPLEMENT_ASSETS',
        'checks': {'exact16_unique_primary_XML_links': True, 'all16_official_sizes_and_MD5s_match': True,
                   'all16_file_signatures_match': True, 'fourPDF_twoXLSX_tenAVI': True,
                   'no_scientific_usability_inferred_from_download': True},
        'assets': assets, 'total_bytes': sum(a['bytes'] for a in assets), 'formats': types,
        'input_hashes': {p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest() for p in (xml_file, metadata_file, listing_file)},
        'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'scope': 'Complete bytes/checksum coverage for the16 supplementary assets linked in the inspected article version. This does not imply complete primary microscopy, a global repository search, specimen-level geometry, accurate segmentation, or a usable cohort observation contract. Scientific inspection is recorded separately.'}
    OUT.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    print(json.dumps({'status': result['status'], 'files': len(assets), 'bytes': result['total_bytes'], 'sha256': hashlib.sha256(OUT.read_bytes()).hexdigest()}))


if __name__ == '__main__':
    main()
