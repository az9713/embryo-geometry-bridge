"""Refresh publication metadata while guarding original numerical/source bytes."""
from pathlib import Path
from hashlib import sha256
import json

ROOT=Path(__file__).resolve().parents[1]
EXCLUDED={'PUBLICATION_MANIFEST.json','EXPORT_RECEIPT.json'}
def digest(p):return sha256(p.read_bytes()).hexdigest()
def main():
    path=ROOT/'PUBLICATION_MANIFEST.json';manifest=json.loads(path.read_bytes())
    old={r['path']:r for r in manifest['files']};rows=[];variants=[]
    for p in sorted(ROOT.rglob('*')):
        if not p.is_file() or '.git' in p.relative_to(ROOT).parts:continue
        relative=p.relative_to(ROOT).as_posix()
        if relative in EXCLUDED:continue
        h=digest(p);row=old.get(relative)
        if row:
            row=row.copy()
            if row['snapshot_sha256']!=h and row['kind']!='generated':
                if p.suffix.lower() not in ('.md','.html'):
                    raise ValueError('Protected original scientific/source bytes changed: '+relative)
                row.update(kind='derived_document',source_sha256=row.get('source_sha256',row['snapshot_sha256']),
                           previous_snapshot_sha256=row['snapshot_sha256'],
                           reason='Public-site presentation, live links and static equation typesetting; numerical scientific files unchanged. Legacy unit typography corrections are logged separately.')
                variants.append(relative)
        else:
            row={'path':relative,'kind':'generated'}
            if relative=='CHARTER.md':row.update(kind='original',source_sha256=h)
            elif relative=='docs/ARCHIVED-AVIAN-FINDINGS.md':
                source=json.loads((ROOT/'docs/JOURNEY-EVIDENCE.json').read_bytes())
                origin=next(x for x in source['sources'] if x['record']=='earlier_avian_findings')
                row.update(kind='derived_document',source_sha256=origin['source_sha256'],
                           reason='Archived historical findings; unavailable local references point to explicit omission notices.')
        row.update(bytes=p.stat().st_size,snapshot_sha256=h);rows.append(row)
    missing=set(old)-{r['path'] for r in rows}
    if missing:raise ValueError('Unexpected removed snapshot files: '+str(sorted(missing)))
    manifest.update(files=rows,included_data_files=len(rows),included_data_bytes=sum(x['bytes'] for x in rows),
                    scope='Curated scientific source, development journey, static web pages and selected receipts. Full historical raw integrations are omitted.',
                    five_gate_goal_complete=False)
    path.write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf8',newline='\n')
    print(json.dumps({'status':'PASS_MANIFEST_REFRESH_WITH_ORIGINAL_DATA_GUARD','files':len(rows),
                      'bytes':manifest['included_data_bytes'],'document_variants_this_refresh':len(variants),'manifest_sha256':digest(path)}))
if __name__=='__main__':main()
