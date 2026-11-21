"""Cold and actual dyadic body/log-domain screening of the same complete K–L trial."""
import argparse
from fractions import Fraction as Q
from hashlib import sha256
from itertools import chain
import json
import os
from pathlib import Path

import readback as R

HERE=Path(__file__).resolve().parent
SELF=Path(__file__).resolve()
READBACK_SHA="ba5450b546371c730389a41f93d1e3d88b37c08c19f613c096f051cf8c0eb403"
COLD_SHA="3fa79d7c14865978c4b891f4b4c44ccb388f23c0ea061004c9a807c714e5584e"
INITIAL_OS="d33e301a76decc888c803227a88058fab7aa1714c373536f49acd36750c340a0"
SPLITS=(2,4,8,16)


def secant_box(box,y0,y1,radius,a,b):
    """Return hull of every original-radius affine fiber in rational normalized-step [a,b]."""
    free=set(box["free"])
    ref=[Q(v) for row in box["reference_nodes_exact"] for v in row]
    center=[]; err=[]
    for i,(u,v,x) in enumerate(zip(y0,y1,ref,strict=True)):
        p,q=(u+t*(v-u) for t in (a,b))
        if i not in free:
            R.require(p==q==x,"fixed coordinates frozen")
        center.append((p+q)/2)
        err.append(abs(q-p)/2+(radius if i in free else 0))
    def paired(values):
        return [[values[i],values[i+1]] for i in range(0,len(values),2)]
    return {**box,"current_centers_exact":paired(center),"coordinate_radii_exact":paired(err)}


def anchors(plan,links):
    bb=[R.load(R.ROOT/path,h,links) for path,h in zip(plan["anchor_box_files"],plan["anchor_box_sha256"],strict=True)]
    return tuple([Q(v) for row in b["current_centers_exact"] for v in row] for b in bb)


def initial_phase(links):
    p=R.load(HERE/"protocol.json",COLD_SHA,links)
    R.require(R.sha(HERE/"readback.py")==READBACK_SHA==p["source_sha256"],"same prior independent body replay")
    o=R.load(HERE/"actual/os_exit_observation.json",INITIAL_OS,links)
    r=R.load(HERE/"actual/receipt.json",o["receipt_sha256"],links)
    t=R.load(HERE/"actual/terminal.json",o["terminal_sha256"],links)
    st=R.load(HERE/"actual/started.json",o["started_sha256"],links)
    R.require(o["status"]=="PASS_OBSERVED_INDEPENDENT_KL_GEOMETRY_ACTUAL_READBACK_EXIT"
              and type(o["pid"]) is type(t["pid"]) is type(st["pid"]) is int and o["pid"]==t["pid"]==st["pid"]==r["pid"]>0
              and type(o["observed_exit_code"]) is type(t["exit_code"]) is int and o["observed_exit_code"]==t["exit_code"]==0
              and t["status"]=="FINISHED" and t["receipt_sha256"]==o["receipt_sha256"]
              and r["status"]=="VERIFIED_UNRESOLVED_KL_WHOLE_ENVELOPE_BODY"
              and r["source_sha256"]==READBACK_SHA and r["protocol_sha256"]==COLD_SHA
              and r["root_or_loading_branch_proved"] is False,"actual observed original-scope unresolved interval bound")


def fixture():
    toy={"free":[1],"reference_nodes_exact":[[0,0],[1,1]]}
    r=Q(1,8); y0=[Q(x) for x in (0,0,1,1)]; y1=[Q(x) for x in (0,2,1,1)]
    b=secant_box(toy,y0,y1,r,Q(1,4),Q(1,2))
    R.require(b["current_centers_exact"]==[[0,Q(3,4)],[1,1]]
              and b["coordinate_radii_exact"]==[[0,Q(3,8)],[0,0]],"exact model-independent slab and unchanged constraints")
    for n in SPLITS:
        links=[Q(i,n) for i in range(n+1)]
        R.require(links[0]==0 and links[-1]==1 and all(a<b for a,b in zip(links[:-1],links[1:],strict=True)),
                  "complete exact dyadic cover")
    return {"toy_closed_subhulls_exact":True,"all_dyadic_protocol_levels_cover_closed_original_interval":True}


def review_context():
    bindings,plan,box,_,_,B,links=R.pinned_sources()
    initial_phase(links)
    y0,y1=anchors(plan,links)
    r=Q(plan["fiber_radius_exact"])
    R.require(r==Q(1,4096) and secant_box(box,y0,y1,r,Q(0),Q(1))["coordinate_radii_exact"]
              == [[Q(v) for v in row] for row in box["coordinate_radii_exact"]],
              "unchanged full-span affine predictor and corrector radius")
    return bindings,plan,box,B,links,y0,y1,r


def policy(bindings,links):
    return {"status":"FROZEN_SAME_KL_SECANT_DYADIC_BODY_LOG_DIAGNOSTIC",
            "source_sha256":R.sha(SELF),"original_bound_phase_sha256":bindings,
            "source_and_observed_unresolved_sha256":links,"checks":fixture(),
            "dyadic_levels":list(SPLITS),"retain_first_level_with_every_slab_body_and_analytic_bulk_rho":True,
            "retains_original_dose_range_linear_predictor_free_constraints_and_radius":True,
            "no_surface_or_cached_analytic_domain_test":True,
            "no_new_force_Hessian_moment_AD_solve_evaluations":True,"actual_slabs_evaluated":0}


def prepare():
    b,_,_,_,links,_,_,_=review_context()
    print(R.save(HERE/"segment_protocol.json",policy(b,links)),flush=True)


def box_screen(B,box):
    ref,nodes=R.axis_hulls(box)
    j=[];m=[]; failed_j=[];failed_spd=[];failed_rho=[];rho=[]
    for i,t in enumerate(box["triangles"]):
        v=R.tri_replay(B,[ref[k] for k in t[:3]],[nodes[k] for k in t])
        if not v["positive_J_for_entire_box"]:
            failed_j.append(i)
        if not v["strong_monotonicity_for_entire_box"]:
            failed_spd.append(i)
        else:
            m.append(Q(v["strong_monotonicity_lower_exact"]))
        if v["J_bound"] is None:
            failed_rho.append(i)
        else:
            pair=R.interval(v["J_bound"])
            j.append(pair)
            _,beta=R.exact_Rho(pair)
            if pair[0]<=Q(1,10**6) or beta>=1:
                failed_rho.append(i)
            else:
                rho.append(beta)
    radial=min(nodes[i//2][0][0]-ref[i//2][0]-Q(float(-(1-1e-7)*float(ref[i//2][0])))
               for i in box["free"] if i%2==0)
    return {"passes_body_J_monotonicity_log_rho_and_radial_bound":not (failed_j or failed_spd or failed_rho) and radial>0,
            "J_enclosure_exact":R.hull(j) if len(j)==144 else None,
            "minimum_monotonicity_exact":min(m) if len(m)==144 else None,
            "maximum_bulk_log_rho_exact":max(rho) if len(rho)==144 else None,
            "radial_nodal_margin_exact":radial,"insufficient_J_slots":failed_j,
            "insufficient_monotonicity_slots":failed_spd,"insufficient_bulk_log_rho_slots":failed_rho}


def run(ph):
    bindings,plan,box,B,links,y0,y1,r=review_context()
    p=R.load(HERE/"segment_protocol.json",ph,{})
    R.require(p==policy(bindings,links),"same retained closed-box protocol")
    out=HERE/"segment_actual"
    out.mkdir(exist_ok=False)
    stamp={"source_sha256":R.sha(SELF),"protocol_sha256":ph,"pid":os.getpid(),
           "prior_full_scope_geometry_readback_os_sha256":INITIAL_OS}
    R.save(out/"started.json",{"status":"RUNNING",**stamp})
    terminal={"status":"FAILED",**stamp,"exit_code":1}
    try:
        history=[]
        picked=None
        i0,i1=map(Q,plan["dose_interval_exact"])
        for n in SPLITS:
            layer=[]
            for j in range(n):
                a,b=Q(j,n),Q(j+1,n)
                subbox=secant_box(box,y0,y1,r,a,b)
                h=box_screen(B,subbox)
                layer.append({"slot":j,"exact_dose_subinterval":[i0+a*(i1-i0),i0+b*(i1-i0)],**h})
            history.append({"slabs":n,"checks":layer})
            if all(x["passes_body_J_monotonicity_log_rho_and_radial_bound"] for x in layer):
                picked=n
                break
        h=R.save(out/"receipt.json",{"status":"CANDIDATE_DYADIC_KL_BODY_LOG_COVERAGE_AWAITS_READBACK",
            **stamp,"exact_original_closed_interval":[i0,i1],"constant_corrector_radius_exact":r,
            "first_complete_tested_dyadic_level":picked,"levels":history,
            "tested_domain_scope":"body monotonicity determinant axis and bulk exact-spatial log expansion only",
            "surface_and_cached_domains_certified":False,"root_or_loading_branch_proved":False,
            "force_H_AD_solve_moment_or_spatial_integration_jobs":0})
        terminal.update(status="FINISHED",exit_code=0,receipt_sha256=h)
        print(json.dumps({"receipt_sha256":h,"first_complete_tested_dyadic_level":picked}),flush=True)
    finally:
        R.save(out/"terminal.json",terminal)


if __name__=="__main__":
    c=argparse.ArgumentParser(description=__doc__)
    s=c.add_subparsers(dest="mode",required=True)
    s.add_parser("prepare")
    r=s.add_parser("run")
    r.add_argument("--protocol-sha256",required=True)
    a=c.parse_args()
    prepare() if a.mode=="prepare" else run(a.protocol_sha256)

