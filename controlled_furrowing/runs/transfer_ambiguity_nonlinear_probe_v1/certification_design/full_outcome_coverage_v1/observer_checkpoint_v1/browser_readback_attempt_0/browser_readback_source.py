"""Read the unchanged standalone report at mobile and desktop widths."""
from functools import partial
from hashlib import sha256
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
import json
from pathlib import Path
from threading import Thread
from urllib.parse import quote, urljoin, urlsplit

from playwright.sync_api import sync_playwright

HERE=Path(__file__).resolve().parent
ROOT=next(p for p in HERE.parents if (p/"PLAN.md").is_file())
PAGE=HERE/"OBSERVERS.html"
OUT=HERE/"browser_readback"
EXPECTED="Observer comparisons at 36 cached-quadrature equilibria"


class Quiet(SimpleHTTPRequestHandler):
    def log_message(self,*args):
        pass


def sha(p):
    return sha256(p.read_bytes()).hexdigest()


def main():
    OUT.mkdir(exist_ok=False)
    build=json.loads((HERE/"build_receipt.json").read_bytes())
    assert all(sha(HERE/item)==h for item,h in build["artifacts"].items())
    server=ThreadingHTTPServer(("127.0.0.1",0),partial(Quiet,directory=str(ROOT)))
    worker=Thread(target=server.serve_forever,daemon=True)
    worker.start()
    href=quote(PAGE.relative_to(ROOT).as_posix())
    address=f"http://127.0.0.1:{server.server_port}/{href}"
    views=[]
    try:
        with sync_playwright() as runtime:
            browser=runtime.chromium.launch(channel="msedge",headless=True)
            for width in (390,1440):
                p=browser.new_page(viewport={"width":width,"height":900},java_script_enabled=False)
                console,failures=[],[]
                p.on("pageerror",lambda err:console.append(str(err)))
                p.on("console",lambda m:console.append(m.text) if m.type=="error" else None)
                p.on("requestfailed",lambda req:failures.append(req.url))
                resp=p.goto(address,wait_until="networkidle")
                content=p.evaluate("""() => ({
                  headline:document.querySelector('h1').textContent,
                  names:[...document.querySelectorAll('h2')].map(x=>x.textContent),
                  tableRows:[...document.querySelectorAll('table')].map(x=>x.tBodies[0].rows.length),
                  text:document.querySelector('main').innerText,
                  width:document.documentElement.scrollWidth,
                  overflow:[...document.querySelectorAll('main,table,figure,img')].map(x=>({
                    left:x.getBoundingClientRect().left,right:x.getBoundingClientRect().right,
                    width:x.clientWidth,scroll:x.scrollWidth})),
                  images:[...document.images].map(x=>({load:x.complete&&x.naturalWidth>0, alt:x.alt,
                                                       w:x.naturalWidth, h:x.naturalHeight})),
                  links:[...document.querySelectorAll('a[href]')].map(x=>x.getAttribute('href'))
                })""")
                link_checks={}
                for h in content.pop("links"):
                    url=urljoin(address,h)
                    assert urlsplit(url).netloc==urlsplit(address).netloc
                    link_checks[h]=p.request.get(url).status
                found=content.pop("text")
                assert all(k in found for k in ("68 / 68","0.15475515644","positive E and negative K",
                                                 "cached-quadrature","exceeded the same frozen allowance"))
                image=content["images"][0]
                assert resp.status==200 and content["headline"]==EXPECTED and content["width"]==width
                assert content["tableRows"]==[2,3] and len(content["names"])==3
                assert len(content["images"])==1 and image["load"] and image["alt"] and image["w"]>800
                assert len(link_checks)==4 and all(x==200 for x in link_checks.values())
                assert all(x["right"]<=width and x["left"]>=0 and x["scroll"]<=x["width"]+1 for x in content["overflow"])
                assert not console and not failures
                png=OUT/f"observers_{width}.png"
                p.screenshot(path=str(png),full_page=True)
                views.append({"width":width,"HTTP_status":resp.status,"DOM":content,"links":link_checks,
                              "console_errors":console,"failed_requests":failures,
                              "screen_file":png.name,"screen_sha256":sha(png)})
                p.close()
            name=browser.version
            browser.close()
    finally:
        server.shutdown()
        server.server_close()
        worker.join(timeout=5)
    out={
        "status":"PASS_ACTUAL_TWO_VIEWPORTS", "browser":"Microsoft Edge "+name,
        "javascript_enabled":False,"browser_script_sha256":sha(Path(__file__)),
        "build_receipt_sha256":sha(HERE/"build_receipt.json"),"HTML_sha256":sha(PAGE),"views":views,
    }
    file=OUT/"receipt.json"
    with file.open("x",encoding="utf-8",newline="\n") as f:
        f.write(json.dumps(out,indent=2)+"\n")
    print(json.dumps({"receipt_sha256":sha(file),"widths":[v["width"] for v in views]}))


if __name__=="__main__":
    main()

