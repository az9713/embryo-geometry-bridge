"""Read actual local HTML in Edge at desktop/mobile widths; no page edits."""
from functools import partial
from hashlib import sha256
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
import json
from pathlib import Path
from threading import Thread
from urllib.parse import urljoin, urlparse

from playwright.sync_api import sync_playwright

HERE = Path(__file__).resolve().parent
STUDY = HERE.parents[2]
SITE = STUDY.parent
OUT = HERE/'browser_readback'


def sha(p):
    return sha256(p.read_bytes()).hexdigest()


class Quiet(SimpleHTTPRequestHandler):
    def log_message(self,*args):
        pass


def main():
    OUT.mkdir(exist_ok=False)
    protected = [HERE/'REPORT.html', HERE/'REPORT.md', HERE/'report_build_receipt.json',
                 *HERE.glob('*.png'), STUDY/'docs/index.html', SITE/'docs/index.html']
    before = {p.relative_to(SITE).as_posix():sha(p) for p in protected}
    server = ThreadingHTTPServer(('127.0.0.1',0), partial(Quiet,directory=str(SITE)))
    Thread(target=server.serve_forever,daemon=True).start()
    checks=[]
    try:
        with sync_playwright() as p:
            browser=p.chromium.launch(channel='msedge',headless=True)
            for width in (1440,390):
                errors=[]; failed=[]
                page=browser.new_page(viewport={'width':width,'height':1000})
                page.on('pageerror',lambda e:errors.append(str(e)))
                page.on('requestfailed',lambda e:failed.append({'url':e.url,'error':e.failure}))
                page.on('console',lambda e:errors.append(e.text) if e.type=='error' else None)
                url=f'http://127.0.0.1:{server.server_port}/'+(HERE/'REPORT.html').relative_to(SITE).as_posix()
                response=page.goto(url,wait_until='networkidle')
                data=page.evaluate(r'''() => {
                    const d=document.body.cloneNode(true);
                    d.querySelectorAll('.math,script,style').forEach(e=>e.remove());
                    return {math:document.querySelectorAll('[data-tex]').length,
                      rendered:document.querySelectorAll('.katex').length,
                      math_errors:document.getElementById('math-errors').textContent,
                      raw_tex:[...new Set(d.textContent.match(/\\[A-Za-z]+/g)||[])],
                      images:[...document.images].map(e=>({src:e.getAttribute('src'),complete:e.complete,width:e.naturalWidth})),
                      width:innerWidth,document_width:document.documentElement.scrollWidth,
                      hrefs:[...document.querySelectorAll('a[href]')].map(e=>e.getAttribute('href'))};
                }''')
                local=[]; external=[]
                for href in data.pop('hrefs'):
                    target=urljoin(url,href)
                    if urlparse(target).netloc==urlparse(url).netloc:
                        local.append({'href':href,'status':page.request.get(target).status})
                    else:
                        assert urlparse(target).scheme=='https'
                        external.append(target)
                shots=[]
                for full in (False,True):
                    shot=OUT/f'REPORT_{width}_{"full" if full else "top"}.png'
                    page.screenshot(path=str(shot),full_page=full)
                    shots.append({'file':shot.name,'sha256':sha(shot)})
                ok=response.status==200 and not errors and not failed and data['math']==data['rendered']==34
                ok=ok and not data['math_errors'] and not data['raw_tex'] and data['document_width']<=width+1
                ok=ok and len(data['images'])==3 and all(x['complete'] and x['width']>0 for x in data['images'])
                ok=ok and all(x['status']==200 for x in local)
                checks.append({'status':'PASS' if ok else 'FAIL','page_errors':errors,'failed_requests':failed,
                               **data,'local_links':local,'external_https_links_not_refetched':external,'screenshots':shots})
                page.close()
            version=browser.version
            browser.close()
    finally:
        server.shutdown(); server.server_close()
    after={p.relative_to(SITE).as_posix():sha(p) for p in protected}
    ok=before==after and all(x['status']=='PASS' for x in checks)
    receipt={'status':'PASS_BROWSER_READBACK' if ok else 'FAIL_BROWSER_READBACK',
             'browser':'Microsoft Edge headless Playwright','browser_version':version,'checks':checks,
             'protected_documents_unchanged':before==after,'document_hashes':after,
             'scope':'Local HTML, math, images and local links. External scientific sources were audited separately; this is not a scientific proof.'}
    path=OUT/'receipt.json'
    with path.open('x',encoding='utf-8') as f:
        json.dump(receipt,f,indent=2); f.write('\n')
    print(json.dumps({'status':receipt['status'],'receipt_sha256':sha(path),
                      'views':[{'width':x['width'],'math':x['math'],'links':len(x['local_links']),
                                'errors':x['page_errors']} for x in checks]}))
    assert ok


if __name__=='__main__':
    main()
