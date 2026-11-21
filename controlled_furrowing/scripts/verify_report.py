"""Exercise every study HTML page in Edge at desktop and mobile widths."""
from contextlib import contextmanager
from functools import partial
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
import json
from pathlib import Path
import threading
from urllib.parse import unquote,urlsplit
from html.parser import HTMLParser
from playwright.sync_api import sync_playwright

ROOT=Path(__file__).resolve().parents[1]
REPO=ROOT.parent


class Links(HTMLParser):
    def __init__(self):
        super().__init__();self.links=[]
    def handle_starttag(self,tag,attrs):
        for key,value in attrs:
            if key in ('href','src') and value:
                self.links.append(value)


class Quiet(SimpleHTTPRequestHandler):
    def log_message(self,*args):
        pass


def main():
    output=ROOT/"evidence/browser"
    output.mkdir(parents=True,exist_ok=True)
    server=ThreadingHTTPServer(('127.0.0.1',0),partial(Quiet,directory=str(REPO)))
    thread=threading.Thread(target=server.serve_forever,daemon=True)
    thread.start()
    checks=[]
    try:
        with sync_playwright() as p:
            browser=p.chromium.launch(channel='msedge',headless=True)
            for width in [1440,390]:
                page=browser.new_page(viewport={'width':width,'height':1000})
                for source in sorted((ROOT/'docs').glob('*.html')):
                    errors=[]
                    callback=lambda error:errors.append(str(error))
                    page.on('pageerror',callback)
                    address=f'http://127.0.0.1:{server.server_port}/'+source.relative_to(REPO).as_posix()
                    response=page.goto(address,wait_until='networkidle')
                    result=page.evaluate(r'''() => {
                        const prose=document.body.cloneNode(true);
                        prose.querySelectorAll('.math,script,style,pre,code').forEach(x=>x.remove());
                        return {errors:document.getElementById('math-errors')?.textContent||'',
                            math:document.querySelectorAll('[data-tex]').length,
                            rendered:document.querySelectorAll('.katex').length,
                            rawProseTeX:[...new Set(prose.textContent.match(/\\[A-Za-z]+/g)||[])],
                            images:[...document.images].every(i=>i.complete&&i.naturalWidth>0),
                            overflow:document.documentElement.scrollWidth>innerWidth+1};
                    }''')
                    parser=Links();parser.feed(source.read_text(encoding='utf-8'))
                    missing=[]
                    for link in parser.links:
                        parsed=urlsplit(link)
                        if parsed.scheme or parsed.netloc or not parsed.path:
                            continue
                        target=(source.parent/unquote(parsed.path)).resolve()
                        if not target.exists():missing.append(link)
                    passed=response.status==200 and not errors and not result['errors'] and result['math']==result['rendered'] and not result['rawProseTeX'] and result['images'] and not result['overflow'] and not missing
                    checks.append({'page':source.name,'width':width,'pass':passed,'browser_errors':errors,'missing_links':missing,**result})
                    if source.stem in ['index','LOCAL-INFERENCE']:
                        page.screenshot(path=str(output/f'{source.stem}_{width}.png'),full_page=True)
                    page.remove_listener('pageerror',callback)
                page.close()
            browser.close()
    finally:
        server.shutdown();server.server_close()
    result={'status':'PASS' if all(c['pass'] for c in checks) else 'FAIL','browser':'Microsoft Edge, Playwright headless','checks':checks,'scope':'HTML/math/assets/local-links/layout; not scientific proof or public deployment'}
    (ROOT/'evidence/report_browser_verification.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'status':result['status'],'page_viewport_checks':len(checks),'failures':[c for c in checks if not c['pass']]}))
    if result['status']!='PASS':raise SystemExit(1)


if __name__=='__main__':
    main()
