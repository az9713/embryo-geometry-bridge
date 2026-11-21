"""Supplement browser checks with a static backslash-letter macro scan."""
from html.parser import HTMLParser
from pathlib import Path
import hashlib,json,re

ROOT=Path(__file__).resolve().parents[1]

class Prose(HTMLParser):
    def __init__(self):super().__init__();self.stack=[];self.parts=[]
    def handle_starttag(self,tag,attrs):
        if tag in ('br','img','hr','meta','link','input','area','base','embed','param','source','track','wbr'):return
        hidden=tag in ('script','style','pre','code') or 'math' in dict(attrs).get('class','').split()
        self.stack.append((tag,hidden))
    def handle_endtag(self,tag):
        for i in range(len(self.stack)-1,-1,-1):
            if self.stack[i][0]==tag:
                self.stack=self.stack[:i];break
    def handle_data(self,value):
        if not any(hidden for _,hidden in self.stack):self.parts.append(value)

def main():
    rows=[]
    for path in sorted((ROOT/'docs').glob('*.html')):
        parser=Prose();parser.feed(path.read_text(encoding='utf-8'))
        rows.append({'page':path.name,'raw_TeX_commands_in_prose':sorted(set(re.findall(r'\\[A-Za-z]+',' '.join(parser.parts)))),
            'HTML_sha256':hashlib.sha256(path.read_bytes()).hexdigest()})
    receipt={'status':'PASS' if all(not r['raw_TeX_commands_in_prose'] for r in rows) else 'FAIL',
        'html_pages':len(rows),'rows':rows,'verifier_source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'scope':'Static prose scan outside declared math/code/scripts/styles. Detects backslash-letter macros, not every possible plain-text mathematical expression. Actual browser parsing/rendering is checked separately.'}
    (ROOT/'evidence/raw_prose_math_scan.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'status':receipt['status'],'pages':len(rows),'findings':[r for r in rows if r['raw_TeX_commands_in_prose']]}))
    if receipt['status']!='PASS':raise SystemExit(1)

if __name__=='__main__':main()
