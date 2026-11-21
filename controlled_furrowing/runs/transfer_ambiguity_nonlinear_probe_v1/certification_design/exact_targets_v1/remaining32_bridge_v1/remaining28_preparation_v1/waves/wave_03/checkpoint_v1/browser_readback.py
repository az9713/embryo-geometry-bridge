"""Verify the frozen Wave 03 report in Edge with local HTTP and JavaScript disabled."""

from functools import partial
from hashlib import sha256
from http.server import SimpleHTTPRequestHandler,ThreadingHTTPServer
import json
from pathlib import Path
from threading import Thread
from urllib.parse import urljoin,urlsplit

from playwright.sync_api import sync_playwright


HERE=Path(__file__).resolve().parent
SITE=next(p for p in HERE.parents if (p/'scripts/build_report.py').is_file()).parent
PAGE=HERE/'REPORT.html'
BUILD=HERE/'build_receipt.json'
OUT=HERE/'browser_readback'


class Quiet(SimpleHTTPRequestHandler):
    def log_message(self,*_):
        pass


def sha(path):
    return sha256(path.read_bytes()).hexdigest()


DOM=r"""() => {
    const n=s=>document.querySelectorAll(s).length;
    const clone=document.body.cloneNode(true);
    clone.querySelectorAll('.math').forEach(x=>x.remove());
    return {width:innerWidth, document_width:document.documentElement.scrollWidth,
        title:document.querySelector('h1').textContent, rows:n('tbody tr'),
        fonts:document.fonts.check('16px KaTeX_Main'),
        raw_tex:(clone.textContent.match(/\\[A-Za-z]+/g)||[]),
        math:[n('.math'),n('.katex'),n('.katex-mathml math'),n('.katex-error'),n('[data-tex]')],
        images:[...document.images].map(e=>({src:e.getAttribute('src'),loaded:e.complete&&e.naturalWidth>0,alt:e.alt})),
        extents:[...document.querySelectorAll('table,.math.display')].map(e=>({
            client:e.clientWidth,scroll:e.scrollWidth,right:e.getBoundingClientRect().right
        })),
        links:[...new Set([...document.querySelectorAll('[href]')].map(e=>e.getAttribute('href')))],
        includes_summary:['16 of 36','58 strictly negative','56 primary IN','155 alternative IN / 13 OUT',
            'arbitrarily short wavelengths','below cell-scale validity'].every(s=>document.body.innerText.includes(s))
    };
}"""


def check(browser,url,width,build):
    view=browser.new_page(viewport={'width':width,'height':960},java_script_enabled=False)
    errors,failures,network=[],[],[]
    view.on('console',lambda e:errors.append(e.text) if e.type=='error' else None)
    view.on('pageerror',lambda e:errors.append(str(e)))
    view.on('requestfailed',lambda e:failures.append({'url':e.url,'failure':e.failure}))
    view.on('response',lambda e:network.append({'url':e.url,'status':e.status}))
    status=view.goto(url,wait_until='networkidle').status
    view.evaluate('document.fonts.ready.then(() => true)')
    d=view.evaluate(DOM)
    link_status={}
    externals=[]
    for link in d.pop('links'):
        target=urljoin(url,link)
        if urlsplit(target).netloc!=urlsplit(url).netloc:
            externals.append(target)
        else:
            link_status[link]=view.request.get(target).status
    images={x['src'] for x in d['images']}
    screenshots={}
    for full in (False,True):
        path=OUT/f'wave03_{width}_{"full" if full else "top"}.png'
        view.screenshot(path=str(path),full_page=full)
        screenshots[path.name]=sha(path)
    view.close()
    html_ok=(status==200 and d['width']==d['document_width']==width
         and d['title']=='Wave 03 exact-spatial roots and transverse stability' and d['rows']==10
         and d['includes_summary'] and d['fonts'] and not d['raw_tex']
         and d['math']==[build['rendered_math_count']]*3+[0,0]
         and all(r['client']>0 and r['scroll']<=r['client']+1 and r['right']<=width for r in d['extents']))
    resources_ok=(not errors and not failures and all(x['status']==200 and urlsplit(x['url']).netloc==urlsplit(url).netloc
          for x in network) and all(x==200 for x in link_status.values())
          and images=={'selfmap.png','transverse_windows.png','slab_guards.png'}
          and all(x['loaded'] and x['alt'] for x in d['images'])
          and set(build['local_links'])==images|link_status.keys()
          and externals==[build['external_primary_source']])
    return {'status':'PASS' if html_ok and resources_ok else 'FAIL','render_ok':html_ok,'resources_ok':resources_ok,
            **d,'local_links_http_status':link_status,'external_source_from_independent_web_check':externals,
            'resources':network,'errors':errors,'request_failures':failures,'screenshots_sha256':screenshots}


def main():
    build=json.loads(BUILD.read_bytes())

    def unchanged():
        return all(sha(SITE/p)==h for p,h in build['input_sha256'].items()) and all(
            sha(HERE/p)==h for p,h in build['file_sha256'].items())

    assert unchanged()
    OUT.mkdir(exist_ok=False)
    server=ThreadingHTTPServer(('127.0.0.1',0),partial(Quiet,directory=str(SITE)))
    Thread(target=server.serve_forever,daemon=True).start()
    url=f'http://127.0.0.1:{server.server_port}/'+PAGE.relative_to(SITE).as_posix()
    try:
        with sync_playwright() as playwright:
            browser=playwright.chromium.launch(channel='msedge',headless=True)
            version=browser.version
            results=[check(browser,url,w,build) for w in (1440,390)]
    finally:
        server.shutdown()
        server.server_close()
    retained=unchanged()
    passed=retained and all(v['status']=='PASS' for v in results)
    p=OUT/'receipt.json'
    with p.open('x',encoding='utf-8',newline='\n') as file:
        json.dump({'status':'PASS_EDGE_WAVE03_DESKTOP_MOBILE' if passed else 'FAIL_EDGE_WAVE03',
           'browser':'Microsoft Edge headless','version':version,'javascript_disabled':True,
           'source_sha256':sha(Path(__file__)),'build_receipt_sha256':sha(BUILD),
           'frozen_inputs_and_outputs_unchanged':retained,'views':results},file,indent=2,ensure_ascii=False)
        file.write('\n')
    print(json.dumps({'receipt_sha256':sha(p),'views':[(r['width'],r['status'],r['render_ok'],r['resources_ok']) for r in results]}))
    if not passed:
        raise ValueError('browser verification failed')


if __name__=='__main__':
    main()
