"""Check the canonical guide's static rendering and retained file identities."""
from pathlib import Path
from functools import partial
from http.server import ThreadingHTTPServer, SimpleHTTPRequestHandler
from html.parser import HTMLParser
from urllib.parse import urlparse, unquote
from hashlib import sha256
import json
import re
import threading
from playwright.sync_api import sync_playwright

ROOT = Path(__file__).resolve().parents[1]
DOC = ROOT/'docs'


class Quiet(SimpleHTTPRequestHandler):
    def log_message(self, *args):
        pass


class Links(HTMLParser):
    def __init__(self):
        super().__init__()
        self.links = []

    def handle_starttag(self, tag, attrs):
        d = dict(attrs)
        for key in ['href', 'src']:
            if key in d:
                self.links.append(d[key])


def main():
    target = DOC/'SHAPE-DESCRIPTION.html'
    build = json.loads((DOC/'SHAPE-DESCRIPTION-BUILD.json').read_text(encoding='utf8'))
    for path, expected in build['sha256'].items():
        assert sha256((ROOT/path).read_bytes()).hexdigest() == expected, path
    parser = Links()
    parser.feed(target.read_text(encoding='utf8'))
    local = []
    for address in parser.links:
        u = urlparse(address)
        if not u.scheme and u.path:
            assert (DOC/unquote(u.path)).resolve().is_file(), address
            local.append(address)
    expected_images = len(re.findall(r'!\[', (DOC/'SHAPE-DESCRIPTION.md').read_text(encoding='utf8')))
    server = ThreadingHTTPServer(('127.0.0.1', 0), partial(Quiet, directory=str(ROOT)))
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    rows = []
    try:
        with sync_playwright() as pw:
            browser = pw.chromium.launch(channel='msedge', headless=True)
            version = browser.version
            for width in [1440, 390]:
                ctx = browser.new_context(viewport={'width': width, 'height': 1000}, java_script_enabled=False)
                page = ctx.new_page()
                errors = []
                page.on('pageerror', lambda e: errors.append(str(e)))
                reply = page.goto('http://127.0.0.1:'+str(server.server_port)+'/docs/SHAPE-DESCRIPTION.html', wait_until='load')
                page.evaluate('document.fonts.ready')
                info = page.evaluate('''()=>({math:document.querySelectorAll('.math').length,katex:document.querySelectorAll('.katex').length,accessible:document.querySelectorAll('math').length,scroll:document.documentElement.scrollWidth,width:document.documentElement.clientWidth,images:document.images.length,missingImages:Array.from(document.images).filter(x=>!x.complete||!x.naturalWidth).length,errors:document.querySelectorAll('.katex-error').length})''')
                assert reply.status == 200 and info['math'] == info['katex'] == info['accessible'] == build['math_expressions'], info
                assert info['scroll'] <= width and info['images'] == expected_images and info['missingImages'] == 0 and info['errors'] == 0 and not errors, info
                rows.append({'viewport_width': width, 'status': 'PASS', 'rendering': info, 'page_errors': errors})
                ctx.close()
            browser.close()
    finally:
        server.shutdown()
        server.server_close()
        thread.join(timeout=5)
    receipt = {'status': 'PASS_SHAPE_GUIDE_BROWSER', 'browser': 'Microsoft Edge headless',
               'version': version, 'javascript_disabled': True,
               'html_sha256': sha256(target.read_bytes()).hexdigest(),
               'verifier_sha256': sha256(Path(__file__).read_bytes()).hexdigest(),
               'retained_hashes_checked': len(build['sha256']),
               'local_links_checked': len(local), 'views': rows,
               'scope': 'Local explanatory-document rendering and file identity. No new research theorem or live-site deployment.'}
    (DOC/'SHAPE-DESCRIPTION-BROWSER.json').write_text(json.dumps(receipt, indent=2)+'\n', encoding='utf8')
    print(json.dumps({'status': receipt['status'], 'math_per_view': build['math_expressions'],
                      'views': len(rows), 'images_per_view': expected_images,
                      'local_links_checked': len(local)}))


if __name__ == '__main__':
    main()
