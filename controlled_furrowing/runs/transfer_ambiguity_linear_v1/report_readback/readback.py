"""Local independent science and Edge readback of the transfer-ambiguity report."""
import ast
from fractions import Fraction as Q
from functools import partial
from hashlib import sha256
from http.server import SimpleHTTPRequestHandler,ThreadingHTTPServer
import json
from pathlib import Path
import re
from threading import Thread
from urllib.parse import urljoin,urlparse

import numpy as np
from playwright.sync_api import sync_playwright

OUT=Path(__file__).resolve().parent
REPORT=OUT.parent
STUDY=OUT.parents[2]
SITE=STUDY.parent
SOURCE=STUDY/"scripts/report_transfer_ambiguity_linear.py"
SOURCE_SHA="9895b0225d3585bea2bc74cac5b9c50ddd867cd6ff4ac1ff8de14da2589294b9"
RESULT_SHA="68ddc210523cb40001c99babb251916ac22c62eecd836d976fac1ab78dc0fe3f"
REVIEW_SHA="8fc50b42694b3db513e6d2166c6d539fbe67c66c2b3c6895c90df12f3856d072"
BUILD_SHA="0398fa96f7cc09a2c4ab94ad55376a110094d09a203a633534eaf5220e9cf1ff"


def check(ok,msg):
    if not ok:
        raise ValueError(msg)


def sha(p):
    return sha256(Path(p).read_bytes()).hexdigest()


def read(p,pin=None):
    raw=Path(p).read_bytes()
    check(pin is None or pin==sha256(raw).hexdigest(),f"changed input {p}")
    return json.loads(raw)


def write(path,d):
    with path.open("x",encoding="utf8",newline="\n") as w:
        json.dump(d,w,indent=2,allow_nan=False)
        w.write("\n")
    return sha(path)


def exact_rank(a):
    a=[r.copy() for r in a]
    piv=0
    for j in range(len(a[0])):
        ix=next((i for i in range(piv,len(a)) if a[i][j]),None)
        if ix is None:
            continue
        a[piv],a[ix]=a[ix],a[piv]
        for i in range(piv+1,len(a)):
            if a[i][j]:
                factor=a[i][j]/a[piv][j]
                a[i]=[u-factor*v for u,v in zip(a[i],a[piv])]
        piv+=1
    return piv


def science():
    built=read(REPORT/"report_build_receipt.json",BUILD_SHA)
    check(sha(SOURCE)==built["source_sha256"]==SOURCE_SHA and built["results_sha256"]==RESULT_SHA
          and built["review_sha256"]==REVIEW_SHA,"frozen builder and inputs")
    before={Path(k):v for k,v in built["files"].items()}
    check({p.name for p in before}=={"REPORT.md","REPORT.html","local_affine_frontier.png"}
          and all(sha(STUDY/p)==v for p,v in before.items()),"three frozen report outputs")
    actual,review=read(REPORT/"actual/results.json",RESULT_SHA),read(REPORT/"actual_independent/receipt.json",REVIEW_SHA)
    check(review["status"]=="PASS_INDEPENDENT_ACTUAL_LOCAL_AFFINE_DIAGNOSTIC" and all(review["checks"].values())
          and all(review[flag] is False for flag in ("source_compatible_witness","validated_equilibrium_or_PDE_claim","novelty_established",
                                                    "larger_response_files_read","new_force_or_equilibrium_computation"))
          and actual["heldout_response_tables_read"] is False and actual["new_equilibrium_solves"] is False,
          "only rounded affine diagnostic")
    old_path="runs/material_cohort_descent_v1/results.json"
    small_path="runs/material_cohort_descent_v1/protocol.json"
    pilot_path="runs/transfer_pair_pilot_v1/preparation/protocol.json"
    parent=read(STUDY/old_path,review["bindings"][old_path])
    small=read(STUDY/small_path,review["bindings"][small_path])
    pilot=read(STUDY/pilot_path,review["bindings"][pilot_path])
    pp_path="runs/transfer_ambiguity_linear_v1/preparation/protocol.json"
    pp=read(STUDY/pp_path,review["bindings"][pp_path])
    large=read(STUDY/pp["large_field"],review["bindings"][pp["large_field"]])
    sm=next(x for x in parent["candidates"] if x["variant"]=="current_disk__spatial_max")["affine_LP"]
    r=list(map(Q,sm["affine_residual_intercepts_exact"]))
    J=[list(map(Q,row)) for row in sm["affine_derivative_coefficients_exact"]]
    check(len(r)==len(J)==68 and all(len(row)==5 for row in J),"68 selected measured-observable rows")
    rk=exact_rank(J)
    singular=np.linalg.svd(np.array(J,dtype=float),compute_uv=False)
    check(rk==5 and singular[-1]>0 and singular[-1]/singular[0]<0.01,"weak direction and exact rank of rounded matrix")
    c=list(map(Q,actual["rounded_affine_contrast_gradient_exact"]))
    md=(REPORT/"REPORT.md").read_text(encoding="utf8")
    check(all(f"{float(v):.8f}" in md for v in c) and pp["columns"]==actual["columns"]
          and pp["columns"]==["log_apical_area","log_basal_area","apical_tension_additive","log_basal_tension","log_common_gain"],
          "reported signed five-component shape gradient and coordinate order")
    trust=[tuple(map(Q,x)) for x in pp["trust_box_exact"]]
    check([(float(a),float(b)) for a,b in trust]==[(-.1,.1),(-.1,.1),(0,.02),(-.1,.1),(-.1,.1)],"trust range and one-sided tension")
    selection=pilot["input_selection"]
    large_diameter,dose,shoulder=Q(selection["large_reported_diameter_um_exact"]),Q(selection["matched_dose_exact"]),Q(pp["shoulder_um_exact"])
    check(large_diameter==2*shoulder and large_diameter==Q("49.8939753") and shoulder==Q("24.94698765")
          and dose==Q("0.89767374") and Q(small["baseline_gain_exact"])==Q(pilot["calibration_candidates"]["baseline"]["gain_exact"])
          and {k:v for k,v in small["parameters"].items() if k!="patch_radius"}
              =={k:v for k,v in large["parameters"].items() if k!="patch_radius"}
          and small["parameters"]["patch_radius"]!=large["parameters"]["patch_radius"],"same material, gain and geometry parameters apart from activation radius")
    g=actual["shoulder"]
    check(0<Q(g["t_exact"][0])<=Q(g["t_exact"][1])<1 and min(map(Q,g["radial_derivative_exact"]))>0
          and Q(g["endpoint_margin_um_exact"])>0,"stored quadratic graph shoulder strictly inside unique edge and positive slope denominator")
    check(r"\delta Q=\delta z_a(0)-\sum_jN_j(t)\,\delta Z_j" in md
          and r"+\frac{Z'(t)}{R'(t)}\sum_jN_j(t)\,\delta R_j" in md,
          "implicit d t = -sum(N_j dR_j)/R' produces correct plus term in dQ")
    r0=max(map(abs,r))
    check(r0==Q(actual["source_baseline_residual_um_exact"]) and f"{float(r0):.8f}"=="5.37983567",
          "baseline source residual, not noise")
    witnesses={(x["mode"],Q(x["budget_exact"]),x["objective"]):x for x in review["LP_endpoints"]}
    worst=Q(0)
    all_series={}
    width_row={}
    for row in actual["rows"]:
        mode,budget=row["mode"],Q(row["budget_exact"])
        allowable=Q(row["allowance_um_exact"])
        check(allowable==(r0+budget if mode=="source_residual_allowance" else budget),"exact displayed axes allowance")
        dirs=[]
        for key,label,sign in (("max_contrast_change","maximize_Q",1),("max_negative_contrast_change","minimize_Q",-1)):
            proof=row[key]
            qlo,qhi=(Q(proof[k]) for k in ("primal_lower_exact","dual_upper_exact"))
            reference=witnesses[(mode,budget,label)]
            step=list(map(Q,proof["step_exact"]))
            check([qlo,qhi]==list(map(Q,reference["signed_objective_exact"]))
                  and step==list(map(Q,reference["step_exact"])) and all(a<=d<=b for d,(a,b) in zip(step,trust))
                  and qlo==sign*sum(a*b for a,b in zip(c,step))<=qhi
                  and qhi-qlo==Q(proof["duality_gap_exact"]) and proof["exact_primal_feasible"] and proof["exact_weak_duality"],
                  "same 20 exact certified endpoints, parameter units and signed objective")
            fx=[sum(a*b for a,b in zip(j,step)) for j in J]
            local=max(abs(v+(r[j] if mode=="source_residual_allowance" else 0)) for j,v in enumerate(fx))
            check(local<=allowable,"exact appropriate tolerance applied per row")
            worst=max(worst,qhi-qlo)
            dirs.append((qlo,qhi,fx))
        wl,wh=[dirs[0][i]+dirs[1][i] for i in range(2)]
        check((wl,wh)==tuple(map(Q,row["contrast_range_width_um_exact"]))
              and [float(wl),float(wh)]==row["contrast_range_width_um_display"],"exact contrast range from both signs and correct float projection")
        all_series.setdefault(mode,[]).append({"x_allowance_um":float(allowable),"y_upper_endpoint_display_um":float(wh),
                                              "x_exact":str(allowable),"y_endpoints_exact":[str(wl),str(wh)]})
        if mode=="small_prediction_change":
            line=f"| {float(allowable)} | {float(wl):.9f} | {float(wh):.9f} |"
            check(line in md,"published baseline near-invisibility table from exact endpoints rounded to nine decimals")
        if mode=="small_prediction_change" and budget==Q(1,1000):
            diff=max(abs(a-b) for a,b in zip(dirs[0][2],dirs[1][2]))
            check(diff<=Q(1,500) and f"{float(wh):.8f}"=="0.06042515"
                  and "0.001 µm of the baseline prediction" in md and "0.002 µm" in md,
                  "0.001 per endpoint relative to baseline and 0.002 maximum possible pair change")
            width_row["near_invisibility"]={"range_approx_um":float(wh),"pair_max_change_um":float(diff),
                                            "exact_pair_max_change_um":str(diff)}
        if mode=="source_residual_allowance" and budget==0:
            check(f"{float(wh):.8f}"=="0.30181702" and "not a measured noise level" in md,
                  "source-relative residual correctly distinguished from experimental error")
            width_row["source_relative"]={"range_approx_um":float(wh),"allowance_um":float(r0)}
    check(len(witnesses)==20 and float(worst)==review["max_exact_duality_gap_um_display"] and
          f"{float(worst):.2e}"=="1.12e-15","twenty exact optimization endpoint gaps")
    ast_tree=ast.parse(SOURCE.read_bytes())
    calls=[x for x in ast.walk(ast_tree) if isinstance(x,ast.Call) and isinstance(x.func,ast.Attribute)
           and x.func.attr in ("plot","loglog")]
    check(len(calls)==2 and {x.func.attr for x in calls}=={"plot","loglog"}
          and all("float(Q(x['allowance_um_exact']))" in ast.unparse(x.args[0])
                  and "x['contrast_range_width_um_display'][1]" in ast.unparse(x.args[1]) for x in calls),
          "plot shows rounded upper endpoints at exact sampled allowances, no interval width depiction")
    check(all(v in md for v in ("rounded affine optimization","floating diagnostics at the saved state",
                                "not a measured noise level","does not establish two source-compatible nonlinear equilibria",
                                "central differences of maxima","one-sided physical directions")),
          "scientific scope and sensitivity kink qualification")
    os=read(REPORT/"actual/os_exit_observation.json")
    term=read(REPORT/"actual/terminal.json",os["worker_terminal_sha256"])
    check(os["observed_exit_code"]==term["exit_code"]==0 and term["pid"]==os["worker_terminal_pid"]
          and os["results_sha256"]==term["results_sha256"]==RESULT_SHA,"true displayed execution evidence")
    return {"status":"PASS_SCIENTIFIC_REPORT_READBACK","shape_values":width_row,
            "reported_rounded_matrix_exact_rank":rk,"numerical_singular_values":[float(v) for v in singular],
            "max_primal_dual_gap_um":float(worst),"verified_LP_optimizations":20,"verified_small_rows":68,
            "figure_points":all_series,
            "figure_semantics":"Float projection of upper exact optimizer width bounds at the ten sampled allowances. Connecting segments are visual guides, and no physical or numerical uncertainty envelope is shown.",
            "sensitivity_scope":"Mechanically floating at saved field. Moving-radius formula differentiates the stored P2 graph only. The exact rank and primal/dual certificates apply to the rounded affine matrix."}


class Quiet(SimpleHTTPRequestHandler):
    def log_message(self,*args):
        pass


def browser():
    files=[REPORT/name for name in ("REPORT.md","REPORT.html","local_affine_frontier.png","report_build_receipt.json")]
    old=[STUDY/"docs/index.html",STUDY/"docs/TRANSFER-CLAIM-AUDIT.html",SITE/"docs/index.html",
         SITE/"docs/assets/katex/katex.min.js",SITE/"docs/assets/katex/katex.min.css"]
    before={str(p.relative_to(SITE)):sha(p) for p in files+old}
    checked=[]
    server=ThreadingHTTPServer(("127.0.0.1",0),partial(Quiet,directory=str(SITE)))
    Thread(target=server.serve_forever,daemon=True).start()
    try:
        with sync_playwright() as p:
            br=p.chromium.launch(channel="msedge",headless=True)
            for width in (1440,390):
                page=br.new_page(viewport={"width":width,"height":1000})
                errors=[];requests=[]
                page.on("pageerror",lambda e:errors.append(str(e)))
                page.on("requestfailed",lambda v:requests.append((v.url,v.failure)))
                page.on("console",lambda e:errors.append(f"console: {e.text}") if e.type=="error" else None)
                url=f"http://127.0.0.1:{server.server_port}/"+(REPORT/"REPORT.html").relative_to(SITE).as_posix()
                reply=page.goto(url,wait_until="networkidle")
                data=page.evaluate(r"""() => {
                    const doc=document.body.cloneNode(true);
                    doc.querySelectorAll('.math,script,style').forEach(e=>e.remove());
                    return {math:document.querySelectorAll('[data-tex]').length,
                      rendered:document.querySelectorAll('.katex').length,
                      mathErrors:document.getElementById('math-errors').textContent,
                      strayTeX:[...new Set(doc.textContent.match(/\\[A-Za-z]+/g)||[])],
                      images:[...document.images].map(e=>({complete:e.complete,width:e.naturalWidth,height:e.naturalHeight,alt:e.alt})),
                      documentWidth:document.documentElement.scrollWidth,viewport:innerWidth,
                      localScrollable:[...document.querySelectorAll('.math.display,table')].filter(e=>e.scrollWidth>e.clientWidth+1).length,
                      links:[...document.querySelectorAll('a[href]')].map(e=>e.getAttribute('href'))};
                }""")
                links=[]
                for href in data.pop("links"):
                    target=urljoin(url,href)
                    check(urlparse(target).netloc==urlparse(url).netloc,"unexpected remote link")
                    v=page.request.get(target)
                    links.append({"href":href,"status":v.status})
                screenshot=OUT/f"REPORT_{width}.png"
                check(not screenshot.exists(),"preserve screenshot readback")
                page.screenshot(path=str(screenshot),full_page=True)
                passed=(reply.status==200 and not errors and not requests and not data["mathErrors"] and not data["strayTeX"]
                        and data["math"]==data["rendered"]>0 and len(data["images"])==1 and data["images"][0]["complete"]
                        and data["images"][0]["width"]>0 and data["documentWidth"]<=width+1
                        and len(links)==10 and all(l["status"]==200 for l in links))
                checked.append({"status":"PASS" if passed else "FAIL","viewport_width_px":width,"http_status":reply.status,
                                "page_errors":errors,"request_failures":requests,**data,"links":links,
                                "screenshot":screenshot.name,"screenshot_sha256":sha(screenshot)})
                page.close()
            version=br.version
            br.close()
    finally:
        server.shutdown()
        server.server_close()
    after={str(p.relative_to(SITE)):sha(p) for p in files+old}
    ok=before==after and all(x["status"]=="PASS" for x in checked)
    return {"status":"PASS_BROWSER_READBACK" if ok else "FAIL_BROWSER_READBACK",
            "browser":"Microsoft Edge via headless Playwright","browser_version":version,
            "checks":checked,"source_documents_unchanged":before==after,"document_hashes":after}


if __name__=="__main__":
    s=science()
    b=browser()
    s_sha=write(OUT/"scientific_readback.json",s)
    b_sha=write(OUT/"browser_readback.json",b)
    print(json.dumps({"science":s["status"],"science_sha256":s_sha,"browser":b["status"],"browser_sha256":b_sha,
          "browser_details":[{k:r[k] for k in ("status","viewport_width_px","math","rendered","documentWidth","page_errors","request_failures")} for r in b["checks"]],
          "singular_values":s["numerical_singular_values"]}))
    check(b["status"]=="PASS_BROWSER_READBACK","browser readback failed")

