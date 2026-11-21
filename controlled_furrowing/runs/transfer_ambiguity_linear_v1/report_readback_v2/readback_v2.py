"""Prove the favicon-only change and repeat the fixed desktop/mobile browser audit."""
import ast
from pathlib import Path
import json
import sys
from urllib.parse import unquote
from xml.etree import ElementTree

HERE=Path(__file__).resolve().parent
RUN=HERE.parent
STUDY=HERE.parents[2]
SITE=STUDY.parent
OLD=RUN/"report_readback"
sys.path.insert(0,str(OLD))
import readback as R

ARCHIVE_SHA="17ac74f8467190db6c3f6044ec6988ff9d99510e12b9c399b8e308b8c9ac54d9"
CURRENT_SHA="79ae626a72f45dda6a7b1a038b5900f27adf4e6ce96884101ebb55d4783bec8b"
ORIGINAL_REVIEW_SHA="32c0e5563dce8832fadc58677c341fb2ea15f8c28d39eb9e5bc11e783152d1e2"
ORIGINAL_BROWSER_SHA="8c1174a5b2e1c522bfcc662cd9d9e6c6ae76388ba8bd0688ecbc2aebaa61ce2c"


def review_bytes():
    current=R.read(RUN/"report_build_receipt.json",CURRENT_SHA)
    archive=R.read(STUDY/current["preservation_manifest"]["file"],ARCHIVE_SHA)
    R.check(archive["status"]=="ARCHIVED_BEFORE_FAVICON_ONLY_EDIT" and current["schema_version"]==2
            and current["preservation_manifest"]["sha256"]==ARCHIVE_SHA and current["generator_rerun"] is False,
            "preserved complete historical build before the edit")
    unchanged={}
    for relative,h in archive["unmodified_input_baseline"].items():
        R.check(R.sha(SITE/relative)==h,f"protected byte change {relative}")
        unchanged[relative]=h
    archived=archive["archived_files"]
    R.check(len(archived)==10 and all(R.sha(STUDY/v["archived_file"])==v["sha256"] for v in archived.values()),
            "all original source, HTML, receipts, review scripts and screenshots preserved")
    for key,v in archived.items():
        if "/report_readback/" in key:
            R.check(R.sha(STUDY/key)==v["sha256"],"unchanged original scientific and browser readback")
    source=STUDY/"scripts/report_transfer_ambiguity_linear.py"
    ast.parse(source.read_bytes())
    html=RUN/"REPORT.html"
    add=current["embedded_favicon_markup"].encode("utf8")
    R.check(add.startswith(b'<link rel="icon" type="image/svg+xml" href="data:image/svg+xml,') and add.endswith(b'">'),
            "embedded local SVG favicon")
    uri=add.decode().split('href="',1)[1][:-2]
    svg=ElementTree.fromstring(unquote(uri.split(",",1)[1]))
    R.check(svg.tag=="{http://www.w3.org/2000/svg}svg" and len(svg)==1,"valid embedded image without external fetch")
    for p in (source,html):
        old=archived[p.relative_to(STUDY).as_posix()]
        original=(STUDY/old["archived_file"]).read_bytes()
        new=p.read_bytes()
        R.check(new.count(add)==1 and new.replace(add,b"",1)==original and
                new.index(add)==new.index(b"</title>")+len(b"</title>"),
                "only one identical favicon insertion in generator and current HTML")
    R.check(R.sha(source)==current["source_sha256"] and
            all(R.sha(STUDY/k)==v for k,v in current["files"].items()),"current v2 report binding")
    science=R.read(OLD/"scientific_readback.json",ORIGINAL_REVIEW_SHA)
    browser=R.read(OLD/"browser_readback.json",ORIGINAL_BROWSER_SHA)
    R.check(science["status"]=="PASS_SCIENTIFIC_REPORT_READBACK"
            and browser["status"]=="FAIL_BROWSER_READBACK","previous science pass and favicon-only browser failure retained")
    old_build,= [v for k,v in archived.items() if k.endswith("/report_build_receipt.json")]
    original=R.read(STUDY/old_build["archived_file"],old_build["sha256"])
    R.check(old_build["sha256"]==current["previous_build_receipt_sha256"] and
            original["results_sha256"]==current["results_sha256"]==R.sha(RUN/"actual/results.json") and
            original["review_sha256"]==current["review_sha256"]==R.sha(RUN/"actual_independent/receipt.json"),
            "same original results and independent scientific evidence")
    return {"status":"PASS_FAVICON_ONLY_SCIENTIFIC_AND_SOURCE_READBACK",
            "current_build_receipt_sha256":CURRENT_SHA,"archive_sha256":ARCHIVE_SHA,
            "original_scientific_readback_sha256":ORIGINAL_REVIEW_SHA,"protected_input_hashes":unchanged,
            "unchanged_science_figures_markdown_KaTeX_and_dossier":True,
            "current_source_sha256":current["source_sha256"],"current_html_sha256":R.sha(html)}


def main():
    science=review_bytes()
    science_pin=R.write(HERE/"favicon_only_readback.json",science)
    R.OUT=HERE
    browser=R.browser()
    R.check(review_bytes()==science,"unchanged scientific, report and old dossier inputs throughout browser test")
    browser["current_build_receipt_sha256"]=CURRENT_SHA
    browser["favicon_only_readback_sha256"]=science_pin
    pin=R.write(HERE/"browser_readback.json",browser)
    print(json.dumps({"status":browser["status"],"browser_receipt_sha256":pin,"favicon_only_receipt_sha256":science_pin,
                      "views":[{k:b[k] for k in ("viewport_width_px","math","rendered","page_errors","request_failures",
                                                 "mathErrors","documentWidth","images")} for b in browser["checks"]]}))
    R.check(browser["status"]=="PASS_BROWSER_READBACK","current report browser audit failed")


if __name__=="__main__":
    main()

