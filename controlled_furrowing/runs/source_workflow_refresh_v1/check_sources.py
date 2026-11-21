"""Bounded primary availability refresh: metadata/text only, no experimental files."""
from pathlib import Path
import datetime, hashlib, json, re, urllib.request, xml.etree.ElementTree as ET
from html.parser import HTMLParser

OUT = Path(__file__).resolve().parent
class Links(HTMLParser):
    def __init__(self): super().__init__(); self.links=[]
    def handle_starttag(self, tag, attrs):
        a=dict(attrs)
        if tag=='a' and a.get('href'): self.links.append(a['href'])
def sha(b): return hashlib.sha256(b).hexdigest()
urls=[
('publisher','https://www.nature.com/articles/s41467-025-62483-6'),
('europepmc','https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12361397/fullTextXML'),
('lab_publications','https://www.kaszalab.com/publications'),
('lab_home','https://www.kaszalab.com/'),
('coauthor_publications','https://herreralab.org/publications/'),
('crossmark','https://crossmark.crossref.org/dialog/?doi=10.1038%2Fs41467-025-62483-6')]
records=[]
for label,url in urls:
    rec={'label':label,'url':url,'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
    try:
        req=urllib.request.Request(url,headers={'User-Agent':'Mozilla/5.0 research-availability-check'})
        with urllib.request.urlopen(req,timeout=30) as f:
            b=f.read(8_000_001);rec.update(http_status=f.status,final_url=f.url,content_type=f.headers.get('Content-Type'))
        if len(b)>8_000_000: raise ValueError('metadata budget exceeded')
        p=OUT/(label+'.txt');p.write_bytes(b);rec.update(file=str(p.relative_to(OUT)),bytes=len(b),sha256=sha(b))
        text=b.decode('utf-8',errors='replace')
        if label=='europepmc':
            root=ET.fromstring(b)
            rec['availability_sections']=[{'kind':s.attrib.get('sec-type'),'text':' '.join(s.itertext())} for s in root.iter() if s.tag in ('sec','notes') and ('availability' in s.attrib.get('sec-type','') or s.find('title') is not None and 'availability' in ''.join(s.find('title').itertext()).lower())]
            rec['external_links']=[x.attrib.get('{http://www.w3.org/1999/xlink}href') for x in root.findall('.//ext-link')]
            rec['supplementary_material_links']=[x.attrib.get('{http://www.w3.org/1999/xlink}href') for x in root.findall('.//supplementary-material')]
            rec['article_title']=''.join(root.find('.//article-title').itertext())
        else:
            p=Links();p.feed(text);rec['links']=sorted(set(p.links))
        rec['repository_links']=[x for x in rec.get('links',rec.get('external_links',[])) if x and any(s in x.lower() for s in ('github','zenodo','figshare','bioimage','biostudies','osf.io','dryad'))]
        rec['status']='ACCESSED'
    except Exception as e:
        rec.update(status='ACCESS_FAILED',error=repr(e))
    records.append(rec)
receipt={'recorded_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'records':records,'scope':'Direct public primary text/link inspection only; no primary microscopy acquisition or workbook re-decode.'}
(OUT/'direct_access.json').write_text(json.dumps(receipt,indent=2,ensure_ascii=True)+'\n',encoding='utf-8')
print(json.dumps([{'label':x['label'],'status':x['status'],'bytes':x.get('bytes'),'repositories':x.get('repository_links'),'availability':x.get('availability_sections'),'error':x.get('error')} for x in records],ensure_ascii=True))
