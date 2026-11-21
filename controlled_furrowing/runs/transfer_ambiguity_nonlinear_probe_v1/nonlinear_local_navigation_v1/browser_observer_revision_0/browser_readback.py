"""Check only the new local navigation section at both dossier landings."""
from functools import partial
from hashlib import sha256
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
import json
from pathlib import Path
from threading import Thread
from urllib.parse import urljoin
from playwright.sync_api import sync_playwright

HERE=Path(__file__).resolve().parent
SITE=HERE.parents[3]
OUT=HERE/'browser_readback'
PAGES=('docs/index.html','controlled_furrowing/docs/index.html')


def sha(p):
    return sha256(p.read_bytes()).hexdigest()


class Quiet(SimpleHTTPRequestHandler):
    def log_message(self,*a):
        pass


def main():
    OUT.mkdir(exist_ok=False)
    nav=json.loads((HERE/'receipt.json').read_bytes())
    paths=list(nav['before_after_navigation_files'])
    before={x:sha(SITE/x) for x in paths}
    assert all(before[x]==v['after'] for x,v in nav['before_after_navigation_files'].items())
    server=ThreadingHTTPServer(('127.0.0.1',0),partial(Quiet,directory=str(SITE)))
    Thread(target=server.serve_forever,daemon=True).start()
    checks=[]
    try:
        with sync_playwright() as p:
            browser=p.chromium.launch(channel='msedge',headless=True)
            for name in PAGES:
                for width in (1440,390):
                    errors=[]; consoles=[]; missing=[]; failed=[]
                    page=browser.new_page(viewport={'width':width,'height':1000})
                    page.on('pageerror',lambda e:errors.append(str(e)))
                    page.on('console',lambda e:consoles.append(e.text) if e.type=='error' else None)
                    page.on('response',lambda r:missing.append(r.url) if r.status>=400 else None)
                    page.on('requestfailed',lambda r:failed.append(r.url))
                    url=f'http://127.0.0.1:{server.server_port}/'+name
                    result=page.goto(url,wait_until='networkidle')
                    data=page.evaluate('''() => {
                       const h=[...document.querySelectorAll('h2')].find(e=>e.textContent==='Nonlinear transfer depends on the observer');
                       if(!h) return {found:false};
                       let n=h.nextElementSibling,links=[],text='';
                       while(n&&n.tagName!=='H2'){text+=n.textContent;links.push(...[...n.querySelectorAll('a[href]')].map(a=>a.getAttribute('href')));n=n.nextElementSibling;}
                       h.scrollIntoView();return {found:true,text,links,width:innerWidth,document_width:document.documentElement.scrollWidth,
                       mathErrors:[...document.querySelectorAll('#math-errors,#latex-error')].map(e=>e.textContent).join('')};
                    }''')
                    links=[{'href':x,'status':page.request.get(urljoin(url,x)).status} for x in data.pop('links',[])]
                    shot=OUT/(name.replace('/','_')+f'_{width}.png')
                    page.screenshot(path=str(shot),full_page=False)
                    known_favicon=all(x.endswith('/favicon.ico') for x in missing)
                    ok=result.status==200 and data.get('found') and data['document_width']<=width+1 and not data['mathErrors']
                    ok=ok and not errors and not failed and (not consoles or bool(missing) and known_favicon)
                    ok=ok and len(links)>0 and all(x['status']==200 for x in links)
                    ok=ok and all(x in data['text'] for x in ('0.154755','5.3898'))
                    checks.append({'status':'PASS' if ok else 'FAIL','page':name,'width':width,**data,
                                   'new_section_links':links,'page_errors':errors,'failed_requests':failed,
                                   'console_errors':consoles,'http_error_resources':missing,
                                   'favicon_scope':'Unchanged implicit favicon omission is recorded, not asserted repaired.',
                                   'screenshot':shot.name,'screenshot_sha256':sha(shot)})
                    page.close()
            version=browser.version
            browser.close()
    finally:
        server.shutdown();server.server_close()
    after={x:sha(SITE/x) for x in paths}
    ok=before==after and all(x['status']=='PASS' for x in checks)
    receipt={'status':'PASS_NEW_LOCAL_NAVIGATION_BROWSER' if ok else 'FAIL_NAVIGATION_BROWSER',
             'browser':'Microsoft Edge headless Playwright','version':version,'checks':checks,
             'navigation_static_receipt_sha256':sha(HERE/'receipt.json'),'protected_hashes':after,
             'files_unchanged_during_browser':before==after,
             'scope':'New links, section text and page-width delta on existing landings; not a fresh whole-dossier scientific or resource audit.'}
    dest=OUT/'receipt.json'
    with dest.open('x',encoding='utf-8') as f:
        json.dump(receipt,f,indent=2);f.write('\n')
    print(json.dumps({'status':receipt['status'],'receipt_sha256':sha(dest),'views':len(checks)}))
    assert ok


if __name__=='__main__':
    main()
