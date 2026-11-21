"""Check every live HTML response against the local publication bytes."""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from hashlib import sha256
from urllib.parse import quote
from urllib.request import Request,urlopen
import argparse
import json
import subprocess

ROOT=Path(__file__).resolve().parents[1]
def main():
    cli=argparse.ArgumentParser();cli.add_argument('--base-url',required=True);cli.add_argument('--output',type=Path,required=True)
    args=cli.parse_args();head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    paths=sorted(p for p in ROOT.rglob('*.html') if '.git' not in p.relative_to(ROOT).parts)
    def check(p):
        relative=p.relative_to(ROOT).as_posix();address=args.base_url.rstrip('/')+'/'+quote(relative,safe='/')+'?revision='+head
        with urlopen(Request(address,headers={'User-Agent':'EmbryoGeometryVerification/1.0'}),timeout=45) as response:
            data=response.read();status=response.status
        expected=sha256(p.read_bytes()).hexdigest();actual=sha256(data).hexdigest()
        return {'path':relative,'http_status':status,'expected_sha256':expected,'served_sha256':actual,'pass':status==200 and expected==actual}
    with ThreadPoolExecutor(max_workers=4) as pool:results=list(pool.map(check,paths))
    record={'status':'PASS_ALL_LIVE_HTML_BYTES' if all(r['pass'] for r in results) else 'FAIL_LIVE_HTML_BYTES',
            'commit':head,'base_url':args.base_url,'html_files':len(paths),'results':results,
            'scope':'HTTP availability and exact served HTML bytes; mathematical and browser evidence are separate.'}
    args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(record,indent=2)+'\n',encoding='utf8')
    print(json.dumps({'status':record['status'],'commit':head,'html_files':len(paths)}))
    if not all(r['pass'] for r in results):raise SystemExit(1)
if __name__=='__main__':main()
