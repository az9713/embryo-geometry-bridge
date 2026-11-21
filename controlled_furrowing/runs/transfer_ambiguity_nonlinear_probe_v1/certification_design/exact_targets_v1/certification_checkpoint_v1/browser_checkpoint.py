"""Read the checkpoint in local headless Edge with no JavaScript execution."""

from functools import partial
from hashlib import sha256
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
import json
from pathlib import Path
from threading import Thread
from urllib.parse import urljoin, urlparse

from playwright.sync_api import sync_playwright


HERE = Path(__file__).resolve().parent
SITE = HERE.parents[4].parent
RUN = HERE.parents[2]
OUT = HERE/'browser_readback'
PAGE = HERE/'CHECKPOINT.html'
BUILD = HERE/'checkpoint_build_receipt.json'


def sha(path):
    return sha256(path.read_bytes()).hexdigest()


class Quiet(SimpleHTTPRequestHandler):
    def log_message(self,*args):
        pass


VISUALS = r'''() => {
    const clone=document.body.cloneNode(true);
    clone.querySelectorAll('.math,script,style').forEach(el=>el.remove());
    const scrolls=[...document.querySelectorAll('table,.math.display')].map(el=>({
        kind:el.tagName, client:el.clientWidth, scroll:el.scrollWidth, view:el.getBoundingClientRect().right
    }));
    return {
      width:innerWidth, document_width:document.documentElement.scrollWidth,
      rendering:{
        wrappers:document.querySelectorAll('.math').length,
        katex:document.querySelectorAll('.katex').length,
        accessible_mathml:document.querySelectorAll('.katex-mathml math').length,
        display:document.querySelectorAll('.katex-display').length,
        errors:document.querySelectorAll('.katex-error').length,
        unresolved:document.querySelectorAll('[data-tex]').length,
      },
      math_fonts_loaded:document.fonts.check('16px KaTeX_Main'),
      raw_tex:[...new Set(clone.textContent.match(/\\[A-Za-z]+/g)||[])],
      images:[...document.images].map(e=>({src:e.getAttribute('src'),alt:e.alt,loaded:e.complete,width:e.naturalWidth})),
      hrefs:[...document.querySelectorAll('a[href],link[href]')].map(el=>el.getAttribute('href')),
      title:document.querySelector('h1').textContent, scrolls,
    };
}'''


def main():
    OUT.mkdir(exist_ok=False)
    build=json.loads(BUILD.read_bytes())
    assert all(sha(HERE/name)==h for name,h in build['file_sha256'].items())
    protected=(PAGE,HERE/'CHECKPOINT.md',BUILD,HERE/'build_checkpoint.py',HERE/'receipt.json',
               HERE/'cached_derivative_headroom.png',RUN/'report/REPORT.html',SITE/'docs/index.html',
               SITE/'controlled_furrowing/docs/index.html')
    before={path.relative_to(SITE).as_posix():sha(path) for path in protected}
    server=ThreadingHTTPServer(('127.0.0.1',0),partial(Quiet,directory=str(SITE)))
    Thread(target=server.serve_forever,daemon=True).start()
    url=f'http://127.0.0.1:{server.server_port}/'+PAGE.relative_to(SITE).as_posix()
    checks=[]
    try:
        with sync_playwright() as p:
            browser=p.chromium.launch(channel='msedge',headless=True)
            version=browser.version
            for width in (1440,390):
                errors,failed,resources=[],[],[]
                page=browser.new_page(viewport={'width':width,'height':1000},java_script_enabled=False)
                page.on('pageerror',lambda e:errors.append(str(e)))
                page.on('console',lambda e:errors.append(e.text) if e.type=='error' else None)
                page.on('requestfailed',lambda r:failed.append({'url':r.url,'error':r.failure}))
                page.on('response',lambda r:resources.append({'url':r.url,'status':r.status}))
                response=page.goto(url,wait_until='networkidle')
                page.evaluate('document.fonts.ready.then(() => true)')
                data=page.evaluate(VISUALS)
                statuses=[]
                for href in dict.fromkeys(data.pop('hrefs')):
                    target=urljoin(url,href)
                    assert urlparse(target).netloc==urlparse(url).netloc
                    statuses.append({'href':href,'status':page.request.get(target).status})
                shots=[]
                for full in (False,True):
                    shot=OUT/f'CHECKPOINT_{width}_{"full" if full else "top"}.png'
                    page.screenshot(path=str(shot),full_page=full)
                    shots.append({'file':shot.name,'sha256':sha(shot)})
                math=data['rendering']
                ok=(response.status==200 and not errors and not failed and not data['raw_tex']
                    and data['document_width']<=width and all(x['view']<=width for x in data['scrolls'])
                    and sum(x['scroll']>x['client']+1 for x in data['scrolls'])==0
                    and math['wrappers']==math['katex']==math['accessible_mathml']==build['compiled_math_blocks']
                    and math['display']==3 and math['errors']==math['unresolved']==0 and data['math_fonts_loaded']
                    and len(data['images'])==1 and all(i['loaded'] and i['width']>0 and i['alt'] for i in data['images'])
                    and all(x['status']==200 for x in statuses)
                    and all(x['status']==200 and urlparse(x['url']).netloc==urlparse(url).netloc for x in resources))
                checks.append({'status':'PASS' if ok else 'FAIL',**data,'page_errors':errors,'failed_requests':failed,
                               'actual_local_resources':resources,'local_links':statuses,'screenshots':shots})
                page.close()
            browser.close()
    finally:
        server.shutdown()
        server.server_close()
    after={path.relative_to(SITE).as_posix():sha(path) for path in protected}
    ok=before==after and all(x['status']=='PASS' for x in checks)
    result={'status':'PASS_EDGE_CHECKPOINT_DESKTOP_MOBILE' if ok else 'FAIL_EDGE_CHECKPOINT',
            'browser':'Microsoft Edge headless Playwright','browser_version':version,'javascript_disabled':True,
            'browser_source_sha256':sha(Path(__file__)),'build_receipt_sha256':sha(BUILD),'protected_hashes':after,
            'protected_documents_unchanged':before==after,'viewports':checks,
            'scope':'Local compilation, appearance and navigation checks, with static KaTeX and local fonts.'}
    receipt=OUT/'receipt.json'
    with receipt.open('x',encoding='utf-8') as f:
        f.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':result['status'],'receipt_sha256':sha(receipt),
                      'views':[{k:x[k] for k in ('status','width','document_width','rendering')} for x in checks]}))
    assert ok


if __name__=='__main__':
    main()
