"""Inspect public workbook package structure beyond already audited cell values.

No extraction to arbitrary paths, macros, formula execution, or source changes.
Opaque binary metadata remains opaque; its package relationship is identified.
"""
from pathlib import Path
from zipfile import ZipFile
from urllib.parse import urlsplit, parse_qs
import hashlib
import json
import posixpath
import re
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[2]
RUN = ROOT/'runs/source_contract_supplements_v1/workbook_package'
OUT = ROOT/'evidence/source_contract_workbook_package_audit.json'
sha = lambda b: hashlib.sha256(b).hexdigest()
S = '{http://schemas.openxmlformats.org/spreadsheetml/2006/main}'


def main():
    if OUT.exists():
        raise RuntimeError('Preserve completed package audit')
    metadata_file = ROOT/'data/raw/PMC12361397.1.json'
    metadata = json.loads(metadata_file.read_text(encoding='utf-8'))
    listing_file = ROOT/'data/raw/oa_object_listing.xml'
    listing = ET.parse(listing_file).getroot()
    objects = {}
    for obj in listing:
        fields = {c.tag.rsplit('}', 1)[-1]: c.text for c in obj}
        if 'Key' in fields:
            objects[fields['Key'].rsplit('/', 1)[-1]] = fields
    urls = {urlsplit(u).path.rsplit('/', 1)[-1]: u for u in metadata['media_urls']}
    reports = []
    formulas = []
    for number in (4, 16):
        name = f'41467_2025_62483_MOESM{number}_ESM.xlsx'
        source = ROOT/'data/raw'/name
        before = source.read_bytes()
        expected_md5 = parse_qs(urlsplit(urls[name]).query)['md5'][0]
        if hashlib.md5(before).hexdigest() != expected_md5 or len(before) != int(objects[name]['Size']):
            raise RuntimeError('Official source checksum/size mismatch')
        with ZipFile(source) as z:
            if z.testzip() is not None:
                raise RuntimeError('ZIP CRC failure')
            names = z.namelist()
            if len(set(names)) != len(names):
                raise RuntimeError('Ambiguous duplicate ZIP member')
            entries = []
            parsed = {}
            for entry in z.infolist():
                data = z.read(entry.filename)
                entries.append({'part': entry.filename, 'bytes': len(data), 'sha256': sha(data)})
                if entry.filename.endswith(('.xml', '.rels')):
                    parsed[entry.filename] = ET.fromstring(data)
            wb = parsed['xl/workbook.xml']
            sheet_map = {e.attrib.get('{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id'): e.attrib for e in wb.findall(S+'sheets/'+S+'sheet')}
            relationships = []
            external = []
            for part, tree in parsed.items():
                if not part.endswith('.rels'):
                    continue
                # Rel parts live one directory below their owning part.
                base = '' if part == '_rels/.rels' else posixpath.dirname(posixpath.dirname(part))
                for e in tree:
                    record = {'part': part, **e.attrib}
                    relationships.append(record)
                    if e.get('TargetMode') == 'External':
                        external.append(record)
                    elif e.get('Target'):
                        target = e.get('Target')
                        resolved = target.lstrip('/') if target.startswith('/') else posixpath.normpath(posixpath.join(base, target))
                        record['resolved_internal_part'] = resolved
                        record['internal_part_present'] = resolved in names
            hidden_rows, hidden_columns, hyperlinks, hidden_sheets, defined_names = [], [], [], [], []
            for e in wb.findall(S+'sheets/'+S+'sheet'):
                if e.get('state', 'visible') != 'visible':
                    hidden_sheets.append(dict(e.attrib))
            for e in wb.findall(S+'definedNames/'+S+'definedName'):
                defined_names.append({'attributes': dict(e.attrib), 'text': e.text})
            worksheet_parts = [n for n in parsed if n.startswith('xl/worksheets/') and n.endswith('.xml')]
            formula_count = 0
            for part in worksheet_parts:
                tree = parsed[part]
                for e in tree.iter():
                    tag = e.tag.rsplit('}', 1)[-1]
                    if tag == 'row' and e.get('hidden') in ('1', 'true'):
                        hidden_rows.append({'part': part, **e.attrib})
                    elif tag == 'col' and e.get('hidden') in ('1', 'true'):
                        hidden_columns.append({'part': part, **e.attrib})
                    elif tag == 'hyperlink':
                        hyperlinks.append({'part': part, **e.attrib})
                    elif tag == 'c':
                        f = e.find(S+'f')
                        if f is not None:
                            formula_count += 1
                            v = e.find(S+'v')
                            formulas.append({'workbook': name, 'part': part, 'cell': e.get('r'), 'formula': f.text, 'formula_attributes': dict(f.attrib), 'cached_value': v.text if v is not None else None})
            binaries = []
            for entry in entries:
                n = entry['part']
                if n.endswith(('.xml', '.rels')):
                    continue
                data = z.read(n)
                rels = [r for r in relationships if r.get('resolved_internal_part') == n]
                binaries.append({**entry, 'declared_relationship_types': [r.get('Type') for r in rels],
                    'schema_decoded': False,
                    'printable_URLs': [x.decode('ascii') for x in re.findall(rb'https?://[^\x00-\x20<>"\x7f-\xff]+', data)],
                    'interpretation': 'Opaque bytes retained in original workbook; relationship identity does not decode scientific meaning.'})
            report = {'source_file': source.relative_to(ROOT).as_posix(), 'source_sha256': sha(before), 'official_md5': expected_md5,
                'source_bytes': len(before), 'ZIP_member_count': len(entries), 'ZIP_uncompressed_bytes': sum(e['bytes'] for e in entries),
                'entries': entries, 'all_XML_and_relationship_parts_parsed': len(parsed), 'worksheet_count': len(worksheet_parts),
                'sheet_definitions': list(sheet_map.values()), 'hidden_worksheets': hidden_sheets, 'hidden_rows': hidden_rows,
                'hidden_column_ranges': hidden_columns, 'defined_names': defined_names, 'external_relationships': external,
                'worksheet_hyperlinks': hyperlinks, 'relationships': relationships, 'formula_cells': formula_count,
                'comment_parts': [n for n in names if 'comment' in n.lower()], 'chart_parts': [n for n in names if '/charts/' in n],
                'embedded_files_or_customXML_parts': [n for n in names if '/embeddings/' in n or n.startswith('customXml/')],
                'opaque_non_XML_parts': binaries,
                'source_unchanged_after_audit': source.read_bytes() == before}
            reports.append(report)
    checks = {
        'both_workbooks_official_size_and_MD5_verified': len(reports) == 2,
        'all_ZIP_CRCs_and_XML_members_parsed': True,
        'all_sources_unchanged': all(r['source_unchanged_after_audit'] for r in reports),
        'hidden_and_external_parts_examined_not_assumed_absent': True,
        'no_formula_or_macro_execution': True,
        'opaque_metadata_not_falsely_decoded': all(not b['schema_decoded'] for r in reports for b in r['opaque_non_XML_parts']),
        'no_new_anatomical_error_or_cohort_identity_claim': True,
    }
    RUN.mkdir(parents=True, exist_ok=True)
    fp = RUN/'formula_inventory.json'
    if fp.exists():
        raise RuntimeError('Preserve formula inventory')
    fp.write_text(json.dumps(formulas, indent=2)+'\n', encoding='utf-8')
    result = {'status': 'PASS_PACKAGE_INVENTORY_ONLY', 'checks': checks, 'workbooks': reports,
        'formula_inventory_file': fp.relative_to(ROOT).as_posix(), 'formula_inventory_sha256': sha(fp.read_bytes()),
        'source_sha256': sha(Path(__file__).read_bytes()),
        'official_metadata_hashes': {metadata_file.relative_to(ROOT).as_posix(): sha(metadata_file.read_bytes()), listing_file.relative_to(ROOT).as_posix(): sha(listing_file.read_bytes())},
        'scope': 'Package structures beyond the prior cell-value audit. Numeric formula caches are inventoried, not recomputed. Public binary workbook metadata remains semantically opaque. No new geometry/time/ROI/acquisition crosswalk or scientific uncertainty follows from this structural audit.'}
    OUT.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    print(json.dumps({'status': result['status'], 'workbooks': len(reports), 'formula_cells': len(formulas), 'opaque_parts': sum(len(r['opaque_non_XML_parts']) for r in reports), 'output_sha256': sha(OUT.read_bytes())}))


if __name__ == '__main__':
    main()
