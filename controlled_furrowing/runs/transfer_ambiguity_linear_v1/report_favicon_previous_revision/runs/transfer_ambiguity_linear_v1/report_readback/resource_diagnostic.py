"""Locate the desktop browser resource failure using unchanged report bytes."""
from functools import partial
from http.server import SimpleHTTPRequestHandler,ThreadingHTTPServer
import json
from threading import Thread

from playwright.sync_api import sync_playwright
from readback import OUT,REPORT,SITE,sha,write


def examine():
    original="04d3d381a547a195fb959be0a4915dc37e745899e9913e984f8c6e31d77ebb82"
    if sha(REPORT/"REPORT.html")!=original:
        raise ValueError("wrong frozen HTML")
    wire=[]
    class Track(SimpleHTTPRequestHandler):
        def log_message(self,format,*args):
            wire.append({"path":self.path,"message":format%args})
    server=ThreadingHTTPServer(("127.0.0.1",0),partial(Track,directory=str(SITE)))
    Thread(target=server.serve_forever,daemon=True).start()
    output=[]
    try:
        with sync_playwright() as p:
            browser=p.chromium.launch(channel="msedge",headless=True)
            for width in (1440,390):
                page=browser.new_page(viewport={"width":width,"height":1000})
                logs=[];errors=[];response=[]
                page.on("console",lambda c:logs.append({"type":c.type,"location":c.location,"text":c.text}) if c.type=="error" else None)
                page.on("pageerror",lambda c:errors.append(str(c)))
                page.on("response",lambda c:response.append({"url":c.url,"status":c.status}) if c.status>=400 else None)
                rel=(REPORT/"REPORT.html").relative_to(SITE).as_posix()
                page.goto(f"http://127.0.0.1:{server.server_port}/{rel}",wait_until="networkidle")
                page.wait_for_timeout(450)
                output.append({"width":width,"console_errors":logs,"page_errors":errors,"response_errors":response})
                page.close()
            browser.close()
    finally:
        server.shutdown()
        server.server_close()
    errors=[x for x in wire if "404" in x["message"]]
    status="FAVICON_ONLY_404" if errors and all(x["path"]=="/favicon.ico" for x in errors) else "OTHER_RESPONSE_PATTERN"
    data={"status":status,"report_sha256":original,"browser":output,"local_http_404s":errors,"report_unchanged":sha(REPORT/"REPORT.html")==original}
    digest=write(OUT/"resource_diagnostic.json",data)
    print(json.dumps({**data,"receipt_sha256":digest}))


if __name__=="__main__":
    examine()

