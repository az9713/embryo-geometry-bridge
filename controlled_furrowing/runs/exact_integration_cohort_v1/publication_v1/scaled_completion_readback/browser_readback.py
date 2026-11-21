"""Targeted live Edge readback of this publication's new pages and landings."""
from functools import partial
from hashlib import sha256
from http.server import ThreadingHTTPServer, SimpleHTTPRequestHandler
import json
from pathlib import Path
import threading

from playwright.sync_api import sync_playwright

ROOT = Path(__file__).resolve().parents[5]
OUT = Path(__file__).resolve().parent
PAGES = ['controlled_furrowing/docs/COHORT-INTEGRATION-L.html']


class Quiet(SimpleHTTPRequestHandler):
    def log_message(self, *args):
        pass


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def main():
    files = [ROOT / page for page in PAGES]
    before = {page: digest(ROOT / page) for page in PAGES}
    assert all(path.exists() for path in files)
    server = ThreadingHTTPServer(('127.0.0.1', 0), partial(Quiet, directory=str(ROOT)))
    threading.Thread(target=server.serve_forever, daemon=True).start()
    rows = []
    try:
        with sync_playwright() as p:
            browser = p.chromium.launch(channel='msedge', headless=True)
            for width in (390, 1440):
                page = browser.new_page(viewport={'width': width, 'height': 1000})
                for relative in PAGES:
                    errors = []
                    callback = lambda value: errors.append(str(value))
                    page.on('pageerror', callback)
                    response = page.goto(f'http://127.0.0.1:{server.server_port}/{relative}', wait_until='networkidle')
                    result = page.evaluate(r'''() => {
                      const prose=document.body.cloneNode(true);
                      prose.querySelectorAll('.math,script,style,pre,code').forEach(e=>e.remove());
                      return {math:document.querySelectorAll('[data-tex]').length,
                        rendered:document.querySelectorAll('.katex').length,
                        renderErrors:document.getElementById('math-errors')?.textContent||'',
                        rawTeX:[...new Set(prose.textContent.match(/\\[A-Za-z]+/g)||[])],
                        images:[...document.images].every(e=>e.complete&&e.naturalWidth>0),
                        overflow:document.documentElement.scrollWidth>innerWidth+1};
                    }''')
                    passed = response.status == 200 and not errors and not result['renderErrors'] \
                        and result['math'] == result['rendered'] and not result['rawTeX'] \
                        and result['images'] and not result['overflow']
                    name = relative.replace('/', '_').replace('.html', '') + f'_{width}.png'
                    image = OUT / name
                    page.screenshot(path=str(image), full_page=False)
                    rows.append({'page': relative, 'width': width, 'pass': passed,
                                 'page_errors': errors, **result, 'screenshot': name,
                                 'screenshot_sha256': digest(image)})
                    page.remove_listener('pageerror', callback)
                page.close()
            browser.close()
    finally:
        server.shutdown()
        server.server_close()
    after = {page: digest(ROOT / page) for page in PAGES}
    receipt = {'status': 'PASS' if all(x['pass'] for x in rows) and before == after else 'FAIL',
               'browser': 'Microsoft Edge via headless Playwright', 'checks': rows,
               'documents_unchanged': before == after, 'document_hashes': after,
               'scope': 'Actual HTML, KaTeX, images and desktop/mobile layout, not scientific correctness.'}
    target = OUT / 'browser_readback.json'
    with target.open('x', encoding='utf8') as stream:
        json.dump(receipt, stream, indent=2)
        stream.write('\n')
    print(json.dumps({'status': receipt['status'], 'views': len(rows),
                      'receipt_sha256': digest(target), 'failures': [r for r in rows if not r['pass']]}))
    if receipt['status'] != 'PASS':
        raise SystemExit(1)


if __name__ == '__main__':
    main()
