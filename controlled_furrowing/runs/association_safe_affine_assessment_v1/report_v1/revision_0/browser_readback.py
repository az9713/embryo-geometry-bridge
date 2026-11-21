"""Read the local checkpoint in headless Edge and retain browser evidence."""

from functools import partial
from hashlib import sha256
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
import json
from pathlib import Path
from threading import Thread
from urllib.parse import urljoin, urlsplit

from playwright.sync_api import sync_playwright


HERE = Path(__file__).resolve().parent
SITE = HERE.parents[3]
OUT = HERE / 'browser_readback'
BUILD = HERE / 'build_receipt.json'
PAGE = HERE / 'ASSESSMENT.html'

DOM = r"""() => {
  const n = s => document.querySelectorAll(s).length;
  const clone = document.body.cloneNode(true);
  clone.querySelectorAll('.math').forEach(e => e.remove());
  return {
    width: innerWidth,
    document_width: document.documentElement.scrollWidth,
    title: document.querySelector('h1').textContent,
    rows: [...document.querySelectorAll('tbody tr')].map(e => [...e.children].map(c => c.textContent.trim())),
    math: [n('.math'), n('.katex'), n('.katex-mathml math'), n('.math.display'), n('.katex-error'), n('[data-tex]')],
    hidden_raw_tex: (clone.textContent.match(/\\[A-Za-z]+/g) || []),
    fonts: document.fonts.check('16px KaTeX_Main'),
    regions: [...document.querySelectorAll('table,.math.display')].map(e => ({
      empty: !e.innerText.trim(), visible: e.clientWidth > 0, client_width: e.clientWidth,
      scroll_width: e.scrollWidth, right: e.getBoundingClientRect().right
    })),
    images: [...document.images].map(e => ({source:e.getAttribute('src'), alt:e.alt, loaded:e.complete&&e.naturalWidth>0})),
    hrefs: [...new Set([...document.querySelectorAll('[href]')].map(e=>e.getAttribute('href')))]
  };
}"""


class Quiet(SimpleHTTPRequestHandler):
    def log_message(self, *_):
        pass


def sha(p):
    return sha256(p.read_bytes()).hexdigest()


def screen(browser, url, build, width):
    page = browser.new_page(viewport={'width':width, 'height':960}, java_script_enabled=False)
    errors, failures, resources = [], [], []
    page.on('pageerror', lambda e: errors.append(str(e)))
    page.on('console', lambda e: errors.append(e.text) if e.type == 'error' else None)
    page.on('requestfailed', lambda req: failures.append({'url':req.url, 'reason':req.failure}))
    page.on('response', lambda res: resources.append({'url':res.url, 'status':res.status}))
    status = page.goto(url, wait_until='networkidle').status
    page.evaluate('document.fonts.ready.then(() => true)')
    shown = page.evaluate(DOM)
    paths = shown.pop('hrefs')
    checks = {}
    for path in paths:
        address = urljoin(url, path)
        assert urlsplit(address).netloc == urlsplit(url).netloc
        checks[path] = page.request.get(address).status
    images = {entry['source'] for entry in shown['images']}
    views = []
    for full in (False, True):
        image = OUT / f'assessment_{width}_{"full" if full else "top"}.png'
        page.screenshot(path=str(image), full_page=full)
        views.append({'path':image.name, 'sha256':sha(image)})
    page.close()
    whole = (status == 200 and shown['title'] == 'How source linkage changes the affine transfer range'
             and shown['width'] == shown['document_width'] == width and len(shown['rows']) == 4
             and shown['rows'][0][1:] == ['0.06042515','0.08470768','40.19%']
             and shown['rows'][2][1:] == ['0.06258378','0.09363464','49.61%']
             and shown['math'] == [build['rendered_math_count']]*3+[3,0,0]
             and shown['fonts'] and not shown['hidden_raw_tex']
             and all(r['visible'] and not r['empty'] and r['right'] <= width
                     and r['scroll_width'] <= r['client_width'] + 1 for r in shown['regions'])
             and images == {'affine_widths.png','source_associations.png'}
             and all(x['loaded'] and x['alt'] for x in shown['images'])
             and images | set(paths) == set(build['local_links'])
             and all(v == 200 for v in checks.values()) and not errors and not failures
             and all(r['status'] == 200 and urlsplit(r['url']).netloc == urlsplit(url).netloc for r in resources))
    return {'status':'PASS' if whole else 'FAIL', **shown, 'request_failures':failures, 'errors':errors,
            'screenshots':views, 'resources':resources, 'local_href_status':checks}


def main():
    OUT.mkdir(exist_ok=False)
    build = json.loads(BUILD.read_bytes())
    assert all(sha(SITE/name) == h for name,h in build['input_sha256'].items())
    assert all(sha(HERE/name) == h for name,h in build['file_sha256'].items())
    server = ThreadingHTTPServer(('127.0.0.1', 0), partial(Quiet, directory=str(SITE)))
    Thread(target=server.serve_forever, daemon=True).start()
    url = f'http://127.0.0.1:{server.server_port}/'+PAGE.relative_to(SITE).as_posix()
    try:
        with sync_playwright() as p:
            browser = p.chromium.launch(channel='msedge', headless=True)
            version = browser.version
            views = [screen(browser, url, build, w) for w in (1440,390)]
            browser.close()
    finally:
        server.shutdown()
        server.server_close()
    same = all(sha(SITE/name) == h for name,h in build['input_sha256'].items())
    passed = same and all(v['status'] == 'PASS' for v in views)
    receipt = OUT / 'receipt.json'
    with receipt.open('x', encoding='utf-8', newline='\n') as f:
        json.dump({'status':'PASS_EDGE_ASSOCIATION_SAFE_AFFINE_DESKTOP_MOBILE' if passed else 'FAIL_EDGE_ASSOCIATION_SAFE_AFFINE',
                   'browser':'Microsoft Edge headless', 'version':version, 'javascript_disabled':True,
                   'source_sha256':sha(Path(__file__)), 'build_receipt_sha256':sha(BUILD),
                   'frozen_inputs_unchanged':same, 'views':views}, f, indent=2)
        f.write('\n')
    print(json.dumps({'receipt_sha256':sha(receipt), 'views':[(v['width'],v['document_width'],v['status']) for v in views]}))
    if not passed:
        raise RuntimeError('browser readback failed')


if __name__ == '__main__':
    main()
