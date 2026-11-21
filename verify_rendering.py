"""Read-only browser checks of the generated portable landing page and report."""
import argparse
from functools import partial
from hashlib import sha256
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
import json
from pathlib import Path
import threading
from playwright.sync_api import sync_playwright

REPORT = 'controlled_furrowing/runs/transfer_ambiguity_nonlinear_probe_v1/certification_design/exact_targets_v1/remaining32_bridge_v1/remaining28_preparation_v1/waves/wave_03/checkpoint_v1/REPORT.html'

class Quiet(SimpleHTTPRequestHandler):
    def log_message(self, *args): pass

def digest(p): return sha256(p.read_bytes()).hexdigest()

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('snapshot',type=Path)
    parser.add_argument('output',type=Path)
    args=parser.parse_args()
    root=args.snapshot.resolve()
    args.output.mkdir(parents=True,exist_ok=False)
    pages=['docs/index.html',REPORT]
    hashes={p:digest(root/p) for p in pages}
    server=ThreadingHTTPServer(('127.0.0.1',0),partial(Quiet,directory=str(root)))
    thread=threading.Thread(target=server.serve_forever,daemon=True);thread.start()
    result=[]
    try:
        with sync_playwright() as pw:
            browser=pw.chromium.launch(channel='msedge',headless=True)
            version=browser.version
            for path in pages:
                for width in (1440,390):
                    context=browser.new_context(viewport={'width':width,'height':1000},java_script_enabled=False)
                    page=context.new_page();errors=[]
                    page.on('pageerror',lambda error:errors.append(str(error)))
                    response=page.goto('http://127.0.0.1:'+str(server.server_port)+'/'+path,wait_until='networkidle')
                    assert response.status==200
                    page.evaluate('document.fonts.ready')
                    info=page.evaluate('''() => ({width:document.documentElement.clientWidth,
                      scroll:document.documentElement.scrollWidth,
                      math:document.querySelectorAll('.math').length,
                      rendered:document.querySelectorAll('.katex').length,
                      accessible:document.querySelectorAll('math').length,
                      failed:document.querySelectorAll('.katex-error').length,
                      images:Array.from(document.images).map(i=>({src:i.getAttribute('src'),ok:i.complete&&i.naturalWidth>0}))})''')
                    expected=24 if path==REPORT else 0
                    assert info['math']==info['rendered']==info['accessible']==expected
                    assert info['failed']==0 and info['scroll']<=width and all(i['ok'] for i in info['images']) and not errors
                    links=page.locator('a[href]').evaluate_all('(xs)=>xs.map(x=>({raw:x.getAttribute("href"),url:x.href}))')
                    checked=[]
                    for link in links:
                        if link['url'].startswith('http://127.0.0.1:'):
                            reply=context.request.get(link['url'])
                            assert reply.status==200,link['raw']
                            checked.append({'path':link['raw'],'status':reply.status})
                    name=('report' if path==REPORT else 'landing')+'_'+str(width)+'.png'
                    page.screenshot(path=str(args.output/name),full_page=True)
                    result.append({'page':path,'width':width,'status':'PASS','geometry':info,'local_links':checked,
                                   'screenshot':name,'screenshot_sha256':digest(args.output/name)})
                    context.close()
            browser.close()
    finally:
        server.shutdown();server.server_close();thread.join(timeout=5)
    assert all(digest(root/p)==h for p,h in hashes.items())
    receipt={'status':'PASS_PORTABLE_SITE_BROWSER','browser':'Microsoft Edge headless','version':version,
             'javascript_disabled':True,'page_sha256':hashes,'views':result,
             'scope':'Portable generated landing and byte-identical latest report at two widths; other pages statically checked.'}
    (args.output/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf8')
    print(json.dumps({'status':receipt['status'],'receipt_sha256':digest(args.output/'receipt.json'),
                      'views':[(r['page'].split('/')[-1],r['width'],r['status']) for r in result]}))

if __name__=='__main__': main()
