"""Verify the portable Wave 03 report snapshot with the Python standard library."""
from collections import Counter
from fractions import Fraction as Q
from hashlib import sha256
from html.parser import HTMLParser
import json
from pathlib import Path
import re
from urllib.parse import unquote, urlsplit

ROOT=Path(__file__).resolve().parent
STUDY=ROOT/"controlled_furrowing"
EX=STUDY/"runs/transfer_ambiguity_nonlinear_probe_v1/certification_design/exact_targets_v1"
W=EX/"remaining32_bridge_v1/remaining28_preparation_v1/waves/wave_03"
T=STUDY/"runs/nonaxisymmetric_surface_sign_preparation_v1"
G=STUDY/"runs/exact_loading_tube_pilot_geometry_actual_readback_v1/root_segment_actual_readback"
MANIFEST="PUBLICATION_MANIFEST.json"
SELF_METADATA={MANIFEST,"EXPORT_RECEIPT.json"}
SCOPE="Manifest hashes, static local links, and selected scalar summaries. No raw integrations or historical process reproduction."


def need(ok,message):
    if not ok:
        raise ValueError(message)


def digest(p):
    return sha256(p.read_bytes()).hexdigest()


def read(p):
    return json.loads(p.read_bytes())


class Page(HTMLParser):
    def __init__(self):
        super().__init__()
        self.urls=[]
        self.ids=set()

    def handle_starttag(self,tag,attrs):
        data=dict(attrs)
        self.ids.update(data[k] for k in ("id","name") if data.get(k))
        for k in ("href","src","poster"):
            if data.get(k):
                self.urls.append((tag,k,data[k]))


def local(path,href,root=ROOT):
    u=urlsplit(href)
    # Every scheme is external for offline traversal, except file escapes.
    need(u.scheme.lower()!="file","local host file URL")
    if u.scheme or u.netloc:
        return None
    base=unquote(u.path)
    candidate=((root/base.lstrip("/")) if base.startswith("/") else (path.parent/base)).resolve() if base else path
    need(candidate.is_relative_to(root),"escaped publication tree")
    if candidate.is_dir():
        candidate=candidate/"index.html"
    return candidate,unquote(u.fragment)


def anchor_exists(p,fragment,parsed):
    if not fragment:
        return True
    suffix=p.suffix.lower()
    if suffix in (".html",".htm"):
        return fragment in parsed[p].ids
    if re.fullmatch(r"L[1-9]\d*(?:-L[1-9]\d*)?",fragment):
        return max(map(int,re.findall(r"\d+",fragment)))<=len(p.read_bytes().splitlines())
    if suffix==".md":
        headings=[t.strip().lower() for t in re.findall(r"^#{1,6}\s+(.*?)\s*#*\s*$",p.read_text(encoding="utf8"),re.M)]
        return fragment.lower() in {re.sub(r"[^\w -]","",t).replace(" ","-") for t in headings}
    if suffix==".svg":
        return fragment in re.findall(r"""\bid=["']([^"']+)["']""",p.read_text(encoding="utf8"))
    return False


def verify_links(files):
    parsed={}
    for path in files:
        if path.suffix.lower() in (".html",".htm"):
            page=Page()
            page.feed(path.read_text(encoding="utf8"))
            parsed[path]=page
    count=0
    for path,page in parsed.items():
        for tag,attr,url in page.urls:
            result=local(path,url)
            if result is None:
                continue
            target,anchor=result
            need(target in files and anchor_exists(target,anchor,parsed),
                 "unresolved static link in "+path.relative_to(ROOT).as_posix())
            count+=1
    for path in files:
        if path.suffix.lower()!=".css":
            continue
        for url in re.findall(r"""url\(\s*['"]?([^'")\s]+)""",path.read_text(encoding="utf8")):
            dest=local(path,url)
            if dest is not None:
                need(dest[0] in files,"unresolved CSS asset")
                count+=1
    return {"html_pages":len(parsed),"local_links_and_assets":count}


def pair(a):
    lo,hi=map(Q,a)
    need(lo<=hi,"reversed scalar interval")
    return lo,hi


def replay():
    root=W/"global_reduction_v1/independent"
    b=read(root/"results.json")
    o=read(root/"root_observed_exit.json")
    need(type(o["observed_exit_code"]) is int and o["observed_exit_code"]==0
         and digest(root/"results.json")==o["results_sha256"],"published observed global summary does not match")
    need(len(b["cases"])==4 and len({x["case_id"] for x in b["cases"]})==4,"four Wave 03 case summaries")
    validated={}
    for record in b["cases"]:
        case=record["case_id"]
        file=root/case/"results.json"
        need(digest(file)==record["case_results_sha256"],"matched final per-case root result")
        r=read(file)
        a=r["arithmetic"]
        eta,z,radius=(Q(v) for v in (a["eta_upper_exact"],a["z_upper_exact"],r["radius_exact"]))
        lhs=eta+z*radius
        e=eta/(1-z)
        need(z==Q(a["ZQ_upper_exact"])+Q(a["b_H_upper_exact"]) and 0<=z<1 and radius==Q(1,2**34)
            and lhs==Q(a["selfmap_lhs_upper_exact"])<radius and eta>=0
            and a["strict_selfmap"] is True and r["root_asserted"] is True
            and record["exact_spatial_root_certified"] is True
            and e==Q(a["posterior_error_upper_exact"])<radius,"strict root self-map and posterior")
        transfer=r["observer_and_optional_minimum"]
        t=r["target"]; sigma,xi,goal=map(Q,(t["archived_center_shift_exact"],t["original_addition_gap_exact"],t["observer_goal_exact"]))
        need(sigma+e<goal and Q(transfer["physical_nodal_error_um_upper_exact"])==Q(67,2)*(sigma+e+xi)
             and transfer["strict_minimum_proved"] is True and Q(transfer["normalized_H_euclidean_lower_exact"])>0
             and record["exact_spatial_minimum_certified"] is True,"published same-root observer and minimum scalar premises")
        validated[case]=(e,sigma,xi,record["case_results_sha256"])
    o=read(W/"observer_application_v1/independent_actual_readback/receipt.json")
    allowance=Q(o["allowance_um_exact"])
    need(set(validated)=={a["case_id"] for a in o["cases"]} and len(o["cases"])==4,"same four archived observer cases")
    obs=Counter()
    for c in o["cases"]:
        e,s,xi,h=validated[c["case_id"]]
        eps=Q(67,2)*(e+s+xi)
        need(c["root_proof_sha256"]==h and (e,s,xi,eps)==tuple(Q(c[k]) for k in
             ("posterior_exact","archived_sigma_exact","original_xi_exact","physical_nodal_error_um_exact")),
             "same frozen archived field and center error")
        roles=("apical","basal")
        errs={}
        for role in roles:
            x=c["graph_premises"][role]
            radial,slope,delta=map(Q,(x["stored_radial_min_um_exact"],x["slope_cap_exact"],x["graph_error_um_exact"]))
            need(radial>8*eps and radial-8*eps==Q(x["strict_radial_margin_um_exact"]) and
                 delta==Q(5,4)*(1+slope)*eps,"strict radial and graph error")
            errs[role]=delta
        need(len(c["comparisons"])==len({(r["variant"],r["observable"]) for r in c["comparisons"]})==16,"complete detailed observer rows")
        for r in c["comparisons"]:
            lo,hi=pair(r["stored_interval_um_exact"])
            point=Q(r["observed_um_exact"])
            err=Q(r["extra_observer_error_um_exact"])
            name=r["observable"]
            expected=errs["apical"] if name.startswith("apical_") else errs["basal"] if name.startswith("basal_") else sum(errs.values())
            a,b=lo-err-point,hi+err-point
            low,high=max(0,a,-b),max(abs(a),abs(b))
            status="IN" if high<=allowance else "OUT" if low>allowance else "UNRESOLVED"
            need(err==expected and pair(r["exact_spatial_root_observer_interval_um_exact"])==(lo-err,hi+err)
                and pair(r["absolute_residual_um_exact"])==(low,high) and r["standing"]==status,
                "observer residual absolute intervals and allowance")
            obs[(r["variant"]=="current_disk__spatial_max",status)]+=1
    need(obs=={(True,"IN"):16,(False,"IN"):44,(False,"OUT"):4},"Wave 03 observer summary counts")
    t=read(T/"independent_actual_readback.json")
    tw=Counter()
    for c in t["case_readback"]:
        ids=set()
        signs=Counter()
        for k in c["all_windows"]:
            a,b=pair(k["window_t_exact"])
            rm,rp=pair(k["reference_R_exact"])
            need(0<a<b<1 and rm>0 and pair(k["independent_r_exact"])[0]>0
                 and pair(k["independent_r_R_exact"])[0]>0,"strict annular material and current trace")
            ids.add((k["face"],k["edge"],a,b))
            d=k["independent_D_exact"]
            lo,hi=(None,None) if d is None else pair(d)
            sign="NONNEGATIVE_ON_WINDOW" if d is None or lo>=0 else "STRICT_NEGATIVE_PATCH" if hi<0 else "UNRESOLVED"
            need(sign==k["independent_standing"],"exact independently bounded stress sign")
            if sign=="STRICT_NEGATIVE_PATCH":
                need(Q(k["independent_tau_negative_upper_exact"])<0
                    and Q(k["independent_minus_tau_r_R_squared_over_Js_lower_exact"])>0,"strict same-subedge bound")
            signs[sign]+=1
        need(len(ids)==len(c["all_windows"])==336 and signs==c["independent_counts"],"same complete stress edge windows")
        tw.update(signs)
    need(tw=={"STRICT_NEGATIVE_PATCH":116,"NONNEGATIVE_ON_WINDOW":552,"UNRESOLVED":4},"all 672 stress windows")
    gp=read(G/"receipt.json")
    levels={}
    for level in gp["levels"]:
        nn=level["slabs"]
        flags=[]
        for slot,row in enumerate(level["rows"]):
            body=G/f"body_{nn}_{slot}.json"
            need(digest(body)==row["body_sha256"] and row["slot"]==slot,"same exact body scalar receipt")
            v=read(body)
            good=all(v[k] is True for k in ("all_element_J_positive","uniform_strong_monotonicity",
                                            "axis_preserved_for_entire_box","global_injectivity_for_entire_box"))
            if good:
                need(Q(v["strong_monotonicity_lower_exact"])>0 and Q(row["max_log_rho_exact"])<1
                     and Q(v["J_bound"]["lower_exact"])>Q(1,10**6),"whole accepted slab kernel geometry")
            need(good is row["passes"],"preserved slab disposition")
            flags.append(good)
        need(len(flags)==nn and all(flags) is level["all_pass"],"all declared dyadic slab slots")
        levels[nn]=flags
    need(set(levels)=={2,4} and not all(levels[2]) and all(levels[4])
         and gp["surface_cache_force_H_or_loading_root_proved"] is False,"geometry only, complete four-slab pass")
    return {"Wave03_selfmaps":4,"Wave03_observer_rows":64,"stress_windows":sum(tw.values()),"whole_box_geometry_slots":6,
            "passed_four_slab_geometry_only":True}


def main():
    manifest=read(ROOT/MANIFEST)
    rows=manifest["files"]
    files=set()
    for row in rows:
        v=row["path"]; p=(ROOT/v).resolve()
        need(p.is_relative_to(ROOT) and p.relative_to(ROOT).as_posix()==v and p not in files and p.is_file(),
             "invalid publication path")
        need(digest(p)==row["snapshot_sha256"] and p.stat().st_size==row["bytes"],"snapshot byte drift "+v)
        if row["kind"]=="original":
            need(row["source_sha256"]==row["snapshot_sha256"],"original proof bytes changed")
        elif row["kind"]=="derived_document":
            need(p.suffix.lower() in (".md",".html") and bool(row["reason"]) and bool(row["source_sha256"]),"explicit document derivative")
        else:
            need(row["kind"]=="generated","known manifest role")
        files.add(p)
    actual={p.resolve() for p in ROOT.rglob("*") if p.is_file()
            and ".git" not in p.relative_to(ROOT).parts and p.relative_to(ROOT).as_posix() not in SELF_METADATA}
    need(actual==files,"unmanifested or missing publication files")
    report={"status":"PASS_PUBLICATION_SNAPSHOT_SCOPE","scope":SCOPE,"included_hashes":len(rows),
            **verify_links(files|{ROOT/MANIFEST}),"selected_scalar_replay":replay()}
    print(json.dumps(report,indent=2))


if __name__=="__main__":
    main()
