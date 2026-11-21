"""Verify the frozen local report and its 72 plotted values in two actual browser viewports."""
from fractions import Fraction
from functools import partial
from hashlib import sha256
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
import csv
import json
from pathlib import Path
from threading import Thread
from urllib.parse import quote, urljoin, urlsplit

from playwright.sync_api import sync_playwright

HERE=Path(__file__).resolve().parent
ROOT=next(p for p in HERE.parents if (p/"PLAN.md").is_file())
OUT=HERE/"browser_readback"
EXPECTED="36 strict minima of the cached-quadrature energies"


class Quiet(SimpleHTTPRequestHandler):
    def log_message(self,*args):
        pass

    def do_GET(self):
        if self.path=="/favicon.ico":
            self.send_response(204)
            self.end_headers()
        else:
            super().do_GET()


def sha(p):
    return sha256(p.read_bytes()).hexdigest()


def main():
    build=json.loads((HERE/"build_receipt.json").read_bytes())
    assert all(sha(HERE/item)==h for item,h in build["artifacts"].items())
    proof=HERE.parent/"independent_readback"/"receipt.json"
    assert sha(proof)==build["independent_minima_sha256"]
    proof_values={r["pilot_id"]:r for r in json.loads(proof.read_bytes())["certificate_cases"]}
    with (HERE/"bounds.csv").open(encoding="utf-8",newline="") as f:
        csv_values=list(csv.DictReader(f))
    fields=("cached_H_eigenvalue_lower_exact","strict_cube_interior_margin_exact")
    assert len(csv_values)==len(proof_values)==36 and len({r["pilot_id"] for r in csv_values})==36
    assert all(Fraction(r[k])==Fraction(proof_values[r["pilot_id"]][k]) for r in csv_values for k in fields)
    expected_table=[[r["pilot_id"],*[f"{float(Fraction(r[k])):.8e}" for k in fields]] for r in csv_values]
    OUT.mkdir(exist_ok=False)
    server=ThreadingHTTPServer(("127.0.0.1",0),partial(Quiet,directory=str(ROOT)))
    thread=Thread(target=server.serve_forever,daemon=True)
    thread.start()
    address=f"http://127.0.0.1:{server.server_port}/"+quote((HERE/"MINIMA.html").relative_to(ROOT).as_posix())
    views=[]
    try:
        with sync_playwright() as P:
            browser=P.chromium.launch(channel="msedge",headless=True)
            for width in (390,1440):
                p=browser.new_page(viewport={"width":width,"height":900},java_script_enabled=False)
                errors,failures=[],[]
                p.on("pageerror",lambda e:errors.append(str(e)))
                p.on("console",lambda e:errors.append(e.text) if e.type=="error" else None)
                p.on("requestfailed",lambda req:failures.append(req.url))
                http=p.goto(address,wait_until="networkidle")
                d=p.evaluate("""() => {
                      const root=document.documentElement;
                      const main=document.querySelector('main');
                      const box=(s)=>[...main.querySelectorAll(s)].map(n=>[n.getBoundingClientRect().left,n.getBoundingClientRect().right]);
                      return {title:document.querySelector('h1').innerText,text:main.innerText,width:root.scrollWidth,
                        heading:[...main.querySelectorAll('h2')].map(n=>n.innerText), box:box('table,figure,img,.metrics'),
                        math:[...main.querySelectorAll('p sub,p sup')].map(n=>n.textContent),
                        rows:[...main.querySelectorAll('tbody>tr')].map(t=>[t.dataset.case,t.cells[1].innerText,t.cells[2].innerText]),
                        images:[...main.querySelectorAll('img')].map(i=>[i.complete,i.naturalWidth,i.naturalHeight,i.alt]),
                        links:[...main.querySelectorAll('a')].map(n=>n.getAttribute('href')) };
                    }""")
                assert http.status==200 and d["title"]==EXPECTED and d["width"]==width
                assert d["rows"]==expected_table and d["math"]==["Q","Q","∞","*","*","2","2"]
                assert all(k in d["text"] for k in ("original native posterior","Exact-spatial","2⁻³⁴","6.59641","4.09887","1.06061"))
                assert len(d["images"])==1 and d["images"][0][0] and d["images"][0][1]>1000 and d["images"][0][-1]
                assert all(left>=0 and right<=width for left,right in d["box"])
                checks={}
                for href in d["links"]:
                    url=urljoin(address,href)
                    assert urlsplit(url).netloc==urlsplit(address).netloc
                    checks[href]=p.request.get(url).status
                assert len(checks)==4 and all(status==200 for status in checks.values()) and not errors and not failures
                screen=OUT/f"cached_minima_{width}.png"
                p.screenshot(path=str(screen),full_page=True)
                views.append({"width":width,"table_rows":len(d["rows"]),"math_tokens":d["math"],
                              "images":d["images"],"scroll_width":d["width"],"links":checks,"console_errors":errors,
                              "failed_requests":failures,"screenshot_sha256":sha(screen),"screenshot":screen.name})
                p.close()
            version=browser.version
            browser.close()
    finally:
        server.shutdown()
        server.server_close()
        thread.join(timeout=5)
    rec={"status":"PASS_ACTUAL_DESKTOP_AND_MOBILE_CACHED_MINIMUM_REPORT","browser":"Microsoft Edge "+version,
         "javascript_enabled":False,"source_sha256":sha(Path(__file__)),"build_receipt_sha256":sha(HERE/"build_receipt.json"),
         "verified_exact_csv_scalar_values":72,"views":views,"raw_physics_jobs":0}
    path=OUT/"receipt.json"
    with path.open("x",encoding="utf-8",newline="\n") as f:
        json.dump(rec,f,indent=2)
        f.write("\n")
    print(sha(path))


if __name__=="__main__":
    main()
