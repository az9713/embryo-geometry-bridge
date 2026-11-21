"""Render every published HTML file at desktop and mobile widths."""
from pathlib import Path
from functools import partial
from hashlib import sha256
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from html.parser import HTMLParser
import argparse
import json
import threading
from urllib.parse import urlsplit,quote,urljoin
import re
from playwright.sync_api import sync_playwright

ROOT=Path(__file__).resolve().parents[1]
class Quiet(SimpleHTTPRequestHandler):
    def log_message(self,*args): pass
def digest(p): return sha256(p.read_bytes()).hexdigest()
class Refresh(HTMLParser):
    def __init__(self):super().__init__();self.target=None
    def handle_starttag(self,tag,attrs):
        d=dict(attrs)
        if tag=='meta' and d.get('http-equiv','').lower()=='refresh':
            m=re.search(r'url\s*=\s*(.*)',d.get('content',''),re.I)
            if m:self.target=m[1].strip().strip('"\'')

def main():
    cli=argparse.ArgumentParser();cli.add_argument('--output',type=Path,required=True)
    cli.add_argument('--base-url');cli.add_argument('--capture-previews',action='store_true')
    cli.add_argument('--javascript',action='store_true')
    args=cli.parse_args();args.output.mkdir(parents=True,exist_ok=False)
    server=None;thread=None
    if args.base_url:
        base=args.base_url.rstrip('/')+'/'
    else:
        server=ThreadingHTTPServer(('127.0.0.1',0),partial(Quiet,directory=str(ROOT)))
        thread=threading.Thread(target=server.serve_forever,daemon=True);thread.start()
        base='http://127.0.0.1:'+str(server.server_port)+'/'
    paths=sorted(p for p in ROOT.rglob('*.html') if '.git' not in p.relative_to(ROOT).parts)
    results=[];failed=[]
    try:
        with sync_playwright() as pw:
            browser=pw.chromium.launch(channel='msedge',headless=True);version=browser.version
            for width in (1440,390):
                context=browser.new_context(viewport={'width':width,'height':1000},java_script_enabled=args.javascript)
                page=context.new_page()
                for p in paths:
                    relative=p.relative_to(ROOT).as_posix();errors=[]
                    def error(e): errors.append(str(e))
                    page.on('pageerror',error)
                    reply=page.goto(base+quote(relative,safe='/'),wait_until='load',timeout=30000)
                    redirect=Refresh();redirect.feed(p.read_text(encoding='utf8'))
                    if redirect.target:
                        page.wait_for_url(urljoin(base+quote(relative,safe='/'),redirect.target),wait_until='load')
                    page.evaluate('document.fonts.ready')
                    info=page.evaluate('''() => ({width:document.documentElement.clientWidth,scroll:document.documentElement.scrollWidth,
                      math:document.querySelectorAll('.math').length,rendered:document.querySelectorAll('.katex').length,
                      accessible:document.querySelectorAll('math').length,failures:document.querySelectorAll('.katex-error').length,
                      missingImages:Array.from(document.images).filter(i=>!i.complete||i.naturalWidth===0).map(i=>i.getAttribute('src'))})''')
                    ok=reply is not None and reply.status==200 and info['scroll']<=width+1 and info['math']==info['rendered']==info['accessible'] and info['failures']==0 and not info['missingImages'] and not errors
                    row={'path':relative,'width':width,'status':'PASS' if ok else 'FAIL','http_status':reply.status if reply else None,
                         'source_sha256':digest(p),'geometry_and_math':info,'page_errors':errors}
                    results.append(row)
                    if not ok: failed.append(row)
                    page.remove_listener('pageerror',error)
                if args.capture_previews and width==1440:
                    for relative,name in [('docs/index.html','dossier-preview.png'),('docs/DEVELOPMENT-JOURNEY.html','journey-preview.png')]:
                        page.goto(base+relative,wait_until='networkidle');page.evaluate('document.fonts.ready')
                        page.screenshot(path=str(ROOT/'docs/assets'/name),full_page=False)
                context.close()
            browser.close()
    finally:
        if server:
            server.shutdown();server.server_close();thread.join(timeout=5)
    receipt={'status':'PASS_ALL_HTML_BROWSER' if not failed else 'FAIL_ALL_HTML_BROWSER','scope':'Every HTML file at two widths, with the declared JavaScript mode.',
             'source_sha256':digest(Path(__file__)),'browser':'Microsoft Edge headless','version':version,'javascript_disabled':not args.javascript,'html_files':len(paths),'view_count':len(results),
             'math_expressions_per_width':sum(x['geometry_and_math']['math'] for x in results if x['width']==1440),
             'base_url':base if args.base_url else 'local static server','views':results,'failures':failed}
    (args.output/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf8')
    print(json.dumps({'status':receipt['status'],'html_files':len(paths),'views':len(results),'math_per_width':receipt['math_expressions_per_width'],
                      'failures':[{'path':x['path'],'width':x['width'],'geometry':x['geometry_and_math']} for x in failed]}))
    if failed: raise SystemExit(1)

if __name__=='__main__':main()
