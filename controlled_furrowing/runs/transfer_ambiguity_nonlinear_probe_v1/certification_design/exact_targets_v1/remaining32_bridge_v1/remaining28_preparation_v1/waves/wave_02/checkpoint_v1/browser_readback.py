"""Read the frozen wave 02 checkpoint in local headless Edge."""

from functools import partial
from hashlib import sha256
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
import json
from pathlib import Path
from threading import Thread
from urllib.parse import urljoin, urlsplit

from playwright.sync_api import sync_playwright


HERE=Path(__file__).resolve().parent
SITE=next(p for p in HERE.parents if (p/'scripts/build_report.py').is_file()).parent
OUT=HERE/'browser_readback'
BUILD=HERE/'build_receipt.json'
PAGE=HERE/'REPORT.html'


class Quiet(SimpleHTTPRequestHandler):
    def log_message(self,*_):
        pass


def sha(path):
    return sha256(path.read_bytes()).hexdigest()


DOM=r"""() => {
    const n=s=>document.querySelectorAll(s).length;
    const clone=document.body.cloneNode(true);
    clone.querySelectorAll('.math').forEach(e=>e.remove());
    return {width:innerWidth, document_width:document.documentElement.scrollWidth,
        title:document.querySelector('h1').textContent, rows:n('tbody tr'),
        fonts:document.fonts.check('16px KaTeX_Main'),
        raw_tex:(clone.textContent.match(/\\[A-Za-z]+/g)||[]),
        math:[n('.math'),n('.katex'),n('.katex-mathml math'),n('.math.display'),n('.katex-error'),n('[data-tex]')],
        images:[...document.images].map(e=>({src:e.getAttribute('src'),loaded:e.complete&&e.naturalWidth>0,alt:e.alt})),
        extents:[...document.querySelectorAll('table,.math.display')].map(e=>({
            client:e.clientWidth, scroll:e.scrollWidth, right:e.getBoundingClientRect().right
        })),
        links:[...new Set([...document.querySelectorAll('[href]')].map(e=>e.getAttribute('href')))],
        includes_summary: ['12 of 36','40 primary IN','111 alternative IN / 9 OUT'].every(s=>document.body.innerText.includes(s))
    };
}"""


def view(browser,url,width,build):
    page=browser.new_page(viewport={'width':width,'height':960},java_script_enabled=False)
    errors,failures,net=[],[],[]
    page.on('pageerror',lambda e:errors.append(str(e)))
    page.on('console',lambda e:errors.append(e.text) if e.type=='error' else None)
    page.on('requestfailed',lambda e:failures.append({'url':e.url,'failure':e.failure}))
    page.on('response',lambda e:net.append({'url':e.url,'status':e.status}))
    main=page.goto(url,wait_until='networkidle').status
    page.evaluate('document.fonts.ready.then(() => true)')
    data=page.evaluate(DOM)
    links={}
    for v in data.pop('links'):
        target=urljoin(url,v)
        assert urlsplit(target).netloc==urlsplit(url).netloc
        links[v]=page.request.get(target).status
    shots={}
    for full in (False,True):
        path=OUT/f'wave02_{width}_{"full" if full else "top"}.png'
        page.screenshot(path=str(path),full_page=full)
        shots[path.name]=sha(path)
    page.close()
    images=data['images']
    fit=(main==200 and data['title']=='Wave 02 exact-spatial roots and observer bounds'
         and data['width']==data['document_width']==width and data['rows']==8 and data['includes_summary']
         and data['fonts'] and not data['raw_tex'] and data['math']==[build['rendered_math_count']]*3+[3,0,0]
         and all(r['client']>0 and r['scroll']<=r['client']+1 and r['right']<=width for r in data['extents'])
         and {r['src'] for r in images}=={'selfmap_wave02.png','observer_failures_wave02.png'}
         and all(r['loaded'] and r['alt'] for r in images)
         and set(links)|{r['src'] for r in images}==set(build['local_links'])
         and all(v==200 for v in links.values()) and not errors and not failures
         and all(x['status']==200 and urlsplit(x['url']).netloc==urlsplit(url).netloc for x in net))
    return {'status':'PASS' if fit else 'FAIL',**data,
            'request_failures':failures,'errors':errors,'local_links':links,'resources':net,'screenshots':shots}


def main():
    OUT.mkdir(exist_ok=False)
    build=json.loads(BUILD.read_bytes())
    unchanged=lambda:all(sha(SITE/k)==v for k,v in build['input_sha256'].items())
    assert unchanged() and all(sha(HERE/k)==v for k,v in build['file_sha256'].items())
    server=ThreadingHTTPServer(('127.0.0.1',0),partial(Quiet,directory=str(SITE)))
    Thread(target=server.serve_forever,daemon=True).start()
    url=f'http://127.0.0.1:{server.server_port}/'+PAGE.relative_to(SITE).as_posix()
    try:
        with sync_playwright() as p:
            browser=p.chromium.launch(channel='msedge',headless=True)
            version=browser.version
            views=[view(browser,url,w,build) for w in (1440,390)]
    finally:
        server.shutdown()
        server.server_close()
    same=unchanged()
    passed=same and all(x['status']=='PASS' for x in views)
    path=OUT/'receipt.json'
    with path.open('x',encoding='utf-8',newline='\n') as f:
        json.dump({'status':'PASS_EDGE_WAVE02_DESKTOP_MOBILE' if passed else 'FAIL_EDGE_WAVE02',
            'browser':'Microsoft Edge headless','version':version,'javascript_disabled':True,
            'source_sha256':sha(Path(__file__)),'build_receipt_sha256':sha(BUILD),
            'frozen_inputs_unchanged':same,'views':views},f,indent=2)
        f.write('\n')
    print(json.dumps({'receipt_sha256':sha(path), 'views':[(x['width'],x['document_width'],x['status']) for x in views]}))
    if not passed:
        raise ValueError('browser readback failed')


if __name__=='__main__':
    main()
