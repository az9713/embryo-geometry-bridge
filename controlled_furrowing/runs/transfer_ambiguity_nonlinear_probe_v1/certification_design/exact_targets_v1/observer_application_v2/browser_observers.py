"""Verify static accessible math and local evidence links in headless Edge."""

from functools import partial
from hashlib import sha256
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
import json
from pathlib import Path
from threading import Thread
from urllib.parse import urljoin, urlsplit

from playwright.sync_api import sync_playwright


HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[4]
SITE=ROOT.parent
PAGE=HERE/'OBSERVERS.html'
BUILD=HERE/'observer_report_build_receipt.json'
OUT=HERE/'browser_readback'


class Quiet(SimpleHTTPRequestHandler):
    def log_message(self,*args):
        pass


def sha(path):
    return sha256(path.read_bytes()).hexdigest()


DOM=r'''() => {
    const clone=document.body.cloneNode(true);
    clone.querySelectorAll('.math').forEach(e=>e.remove());
    const extent=e=>({width:e.clientWidth,scroll_width:e.scrollWidth,right:e.getBoundingClientRect().right});
    const count=s=>document.querySelectorAll(s).length;
    return {width:innerWidth, document_width:document.documentElement.scrollWidth,
        title:document.querySelector('h1').textContent, table_rows:count('tbody tr'),
        math: {wrappers:count('.math'), katex:count('.katex'), mathml:count('.katex-mathml math'),
               errors:count('.katex-error'), unresolved:count('[data-tex]'), display:count('.katex-display')},
        fonts:document.fonts.check('16px KaTeX_Main'),
        images:[...document.images].map(e=>({loaded:e.complete && e.naturalWidth>0, alt:e.alt,src:e.getAttribute('src')})),
        extent:[...document.querySelectorAll('table,.math.display')].map(extent),
        raw_tex:[...new Set(clone.textContent.match(/\\[A-Za-z]+/g)||[])],
        links:[...document.querySelectorAll('[href]')].map(e=>e.getAttribute('href'))};
}'''


def inspect(browser,url,width,build):
    page=browser.new_page(viewport={'width':width,'height':960},java_script_enabled=False)
    errors,failed,pages=[],[],[]
    page.on('pageerror',lambda e:errors.append(str(e)))
    page.on('console',lambda e:errors.append(e.text) if e.type=='error' else None)
    page.on('requestfailed',lambda r:failed.append({'url':r.url,'error':r.failure}))
    page.on('response',lambda r:pages.append({'url':r.url,'status':r.status}))
    response=page.goto(url,wait_until='networkidle')
    page.evaluate('document.fonts.ready.then(() => true)')
    observed=page.evaluate(DOM)
    links=[]
    for href in dict.fromkeys(observed.pop('links')):
        target=urljoin(url,href)
        assert urlsplit(target).netloc==urlsplit(url).netloc
        links.append({'href':href,'status':page.request.get(target).status})
    screenshots=[]
    for full in (False,True):
        img=OUT/f'OBSERVERS_{width}_{"full" if full else "top"}.png'
        page.screenshot(path=str(img),full_page=full)
        screenshots.append({'path':img.name,'sha256':sha(img)})
    page.close()
    math=observed['math']
    good=(response.status==200 and observed['title']=='Observer bounds at four exact-spatial roots'
          and not errors and not failed and not observed['raw_tex'] and observed['table_rows']==7
          and observed['document_width']==width and all(e['right']<=width and e['scroll_width']<=e['width']+1
                                                        for e in observed['extent'])
          and math['wrappers']==math['katex']==math['mathml']==build['rendered_math_count']
          and math['display']==3 and math['errors']==math['unresolved']==0 and observed['fonts']
          and len(observed['images'])==1 and all(x['loaded'] and x['alt'] for x in observed['images'])
          and set(x['href'] for x in links)|{observed['images'][0]['src']}==set(build['local_links'])
          and all(x['status']==200 for x in links) and all(x['status']==200 and urlsplit(x['url']).netloc==urlsplit(url).netloc
                                                        for x in pages))
    return {'status':'PASS' if good else 'FAIL',**observed,'errors':errors,'request_failures':failed,
            'resources':pages,'linked_local_responses':links,'screenshots':screenshots}


def main():
    OUT.mkdir(exist_ok=False)
    build=json.loads(BUILD.read_bytes())
    assert all(sha(HERE/name)==h for name,h in build['file_sha256'].items())
    inputs={name:h for name,h in build['input_sha256'].items()}
    assert all(sha(SITE/name)==h for name,h in inputs.items())
    server=ThreadingHTTPServer(('127.0.0.1',0),partial(Quiet,directory=str(SITE)))
    Thread(target=server.serve_forever,daemon=True).start()
    url=f'http://127.0.0.1:{server.server_port}/'+PAGE.relative_to(SITE).as_posix()
    try:
        with sync_playwright() as playwright:
            browser=playwright.chromium.launch(channel='msedge',headless=True)
            version=browser.version
            views=[inspect(browser,url,width,build) for width in (1440,390)]
            browser.close()
    finally:
        server.shutdown()
        server.server_close()
    same=all(sha(SITE/name)==h for name,h in inputs.items())
    good=same and all(view['status']=='PASS' for view in views)
    report=OUT/'receipt.json'
    with report.open('x',encoding='utf-8') as f:
        json.dump({'status':'PASS_EDGE_EXACT_SPATIAL_OBSERVERS_DESKTOP_MOBILE' if good else 'FAIL_EDGE_EXACT_SPATIAL_OBSERVERS',
                   'browser':'Microsoft Edge headless','browser_version':version,'javascript_disabled':True,
                   'source_sha256':sha(Path(__file__)),'build_receipt_sha256':sha(BUILD),
                   'input_evidence_and_old_reports_remained_byte_identical':same,'views':views},f,indent=2)
        f.write('\n')
    print(json.dumps({'views':[(v['width'],v['document_width'],v['math'],v['status']) for v in views],
                      'browser_receipt_sha256':sha(report)}))
    assert good


if __name__=='__main__':
    main()
