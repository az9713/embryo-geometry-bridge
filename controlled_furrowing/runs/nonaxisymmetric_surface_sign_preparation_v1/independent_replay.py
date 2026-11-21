"""Independent exact center-P2 and perturbation stress-window enclosure. No producer imports."""
import argparse
import base64
from collections import Counter
from fractions import Fraction as Q
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from io import BytesIO
import json
from math import comb
from pathlib import Path

import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
PROD=ROOT/"runs/nonaxisymmetric_adequacy_audit_v1"
EX=ROOT/"runs/transfer_ambiguity_nonlinear_probe_v1/certification_design/exact_targets_v1"
OLD=EX/"negative_KL_active_energy_v1"
INPUT_HELPER=OLD/"independent_actual/replay.py"
INPUT_PIN="932515cc6725ecc48fe29dd58e0d4b181cc9dd58fdb391a96da49e9d8c049353"
PHASE={
"source_sha256":"f33aaf277eea10fd364a529f95d369551c99b73747f7856cee3af4318b60b3a3",
"protocol_sha256":"a09fb767d9990fce5816203d45579a5b400a1271248788b685a91aebda355e32",
"preflight_sha256":"7eb3c81037d1ecac47f5e5ca45122f8abd0a316f5022314c14dd6e1470220206",
"review_sha256":"64464f4f1825a78ec478490b804772b87eafeb0b8671846b5becd5cbcd5fbb3b"}
COLD={
PROD/"surface_stress.py":PHASE["source_sha256"],PROD/"protocol.json":PHASE["protocol_sha256"],
PROD/"preflight.json":PHASE["preflight_sha256"],PROD/"independent_prelaunch_review.json":PHASE["review_sha256"],
INPUT_HELPER:INPUT_PIN,OLD/"protocol.json":"06dd35434358ad7dcb993349281ed6bf8a84fd0a34f7ca2d4932ec0792842361"}
OS_STATUS="PASS_ROOT_OBSERVED_TRANSVERSE_SURFACE_STRESS_EXIT"
WINDOWS=tuple((Q(i,16),Q(i+1,16)) for i in range(1,15))
P2=((Q(1),Q(-3),Q(2)),(Q(0),Q(4),Q(-4)),(Q(0),Q(-1),Q(2)))
MAX_SUBDEPTH=4
FACES=("basal","apical")


def require(ok,why):
    if not ok:
        raise ValueError(why)


def sha(b):
    return sha256(b).hexdigest()


def rel(p):
    return p.resolve().relative_to(ROOT).as_posix()


def load(path,want,checked,parse=True):
    data=path.read_bytes()
    name=rel(path)
    require(sha(data)==want==checked.setdefault(name,want),"closed immutable bytes "+name)
    return json.loads(data) if parse else data


def save(path,doc):
    with path.open("x",encoding="utf-8",newline="\n") as f:
        json.dump(doc,f,indent=2)
        f.write("\n")
    return sha(path.read_bytes())


def exact(x):
    require(type(x) is str,"exact rational string required")
    return Q(x)


def pair(v):
    require(type(v) is list and len(v)==2,"exact ordered pair")
    a,b=map(exact,v)
    require(a<=b,"interval in correct order")
    return a,b


def plus(*v):
    n=max(map(len,v))
    return [sum(r[i] for r in v if i<len(r)) for i in range(n)]


def product(a,b):
    entries=[Q(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):
            entries[i+j]+=x*y
    return entries


def scale(a,k):
    return [k*x for x in a]


def diff(a):
    return [k*x for k,x in enumerate(a) if k]


def center_bound(a,lo,hi):
    # Translate the exact center polynomial into x=t-mid. Enclose signed
    # powers directly, retaining nonnegative even powers. No Bernstein basis.
    mid=(lo+hi)/2
    h=(hi-lo)/2
    q=[sum(a[j]*comb(j,i)*mid**(j-i) for j in range(i,len(a))) for i in range(len(a))]
    lb=ub=q[0]
    for i,v in enumerate(q[1:],1):
        w=v*h**i
        lb+=-abs(w) if i%2 else min(0,w)
        ub+=abs(w) if i%2 else max(0,w)
    return lb,ub


def maxabs(p,lo,hi):
    return max(map(abs,center_bound(p,lo,hi)))


def d_at(p,t):
    return p[0]+p[1]*t


def trace(center,ref,free,e,K,T):
    l=ref[1]-ref[0]
    require(l>0 and ref[0]>=0 and e>=0 and K>=0 and T>=0,"ordered flat material trace and nonnegative law")
    radii=[ref[0],l]
    coeff=[plus(*(scale(p,pt[i]) for p,pt in zip(P2,center,strict=True))) for i in (0,1)]
    r,zz=coeff
    rR,zR=(scale(diff(p),1/l) for p in coeff)
    alpha=1-T/K if K else Q(0)
    Dn=None
    if alpha>0:
        Dn=plus(product(product(r,r),plus(product(rR,rR),product(zR,zR))),
                    scale(product(radii,radii),-alpha**2))
    slopes=[scale(diff(p),1/l) for p in P2]
    def d_error(axis,lo,hi):
        return e*max(sum(abs(d_at(p,t)) for p,f in zip(slopes,free,strict=True) if f[axis]) for t in (lo,hi))
    def compute(lo,hi):
        rr=(radii[0]+l*lo,radii[0]+l*hi)
        # Sum |N_j| <= 5/4 on [0,1] (the single negative endpoint basis
        # on each half is bounded below by -1/8).
        ar=Q(5,4)*e if any(f[0] for f in free) else Q(0)
        dr,dz=d_error(0,lo,hi),d_error(1,lo,hi)
        a,b=center_bound(r,lo,hi)
        x,y=center_bound(rR,lo,hi)
        ab=(a-ar,b+ar)
        xy=(x-dr,y+dr)
        if Dn is None:
            return rr,ab,xy,None
        # For S=r_R^2+z_R^2, use |delta(r^2 S)| <=
        # |delta r^2| max(S0) + (max(r0^2)+|delta r^2|)|delta S|.
        mr=max(abs(a),abs(b))
        ms=max(abs(x),abs(y))
        mz=maxabs(zR,lo,hi)
        rv=2*mr*ar+ar**2
        ss=2*ms*dr+dr**2+2*mz*dz+dz**2
        offset=rv*(ms**2+mz**2)+(mr**2+rv)*ss
        lo0,hi0=center_bound(Dn,lo,hi)
        return rr,ab,xy,(lo0-offset,hi0+offset)
    return alpha,compute


def merge_parts(a,b):
    return tuple(None if left is None else (min(left[0],right[0]),max(left[1],right[1]))
                 for left,right in zip(a,b,strict=True))


def sign(D):
    return "NONNEGATIVE_ON_WINDOW" if D is None or D[0]>=0 else "STRICT_NEGATIVE_PATCH" if D[1]<0 else "UNRESOLVED"


def refined(compute,lo,hi,depth=0):
    out=compute(lo,hi)
    if out[0][0]>0 and out[1][0]>0 and out[2][0]>0 and sign(out[3])!="UNRESOLVED":
        return out,1
    if depth==MAX_SUBDEPTH:
        return out,1
    m=(lo+hi)/2
    a,n=refined(compute,lo,m,depth+1)
    b,k=refined(compute,m,hi,depth+1)
    return merge_parts(a,b),n+k


def bound_row(trace_data,K,win):
    alpha,calc=trace_data
    (R,r,s,D),parts=refined(calc,*win)
    require(R[0]>0 and r[0]>0 and s[0]>0,"annular geometry must be strictly positive throughout each checked window")
    choose=sign(D)
    tau=K*D[1]/(2*alpha*R[1]**2) if choose=="STRICT_NEGATIVE_PATCH" else None
    desc=-tau*s[0]**2/alpha if tau is not None else None
    return R,r,s,D,choose,tau,desc,parts


def old_result(row,R,r,s,D,alpha,K):
    require(pair(row["ref_radial_interval_exact"])==R,"same exact annular material window")
    pr,ps=(exact(row[k]) for k in ("deformed_r_lower_exact","deformed_r_R_lower_exact"))
    require(0<pr<=r[1] and 0<ps<=s[1],"producer and independent radius/slope artifacts do not contradict")
    pd=None if D is None else pair(row["D_R2_times_Js2_minus_alpha2_exact"])
    require((pd is None)==(D is None) and (D is not None or row["D_R2_times_Js2_minus_alpha2_exact"] is None),
            "same exact alpha branch")
    overlap=pd is None or max(D[0],pd[0])<=min(D[1],pd[1])
    stand=row["standing"]
    require(stand==sign(pd) and overlap,"producer exact sign arithmetic or actual scalar bounds conflict")
    if stand=="STRICT_NEGATIVE_PATCH":
        tau=K*pd[1]/(2*alpha*R[1]**2)
        c=-tau*ps**2/alpha
        require(exact(row["tau_negative_upper_exact"])==tau<0 and exact(row["minus_tau_r_R_squared_over_Js_lower_exact"])==c>0,
                "producer strict tau and c exact inequalities")
    else:
        require(row["tau_negative_upper_exact"] is None and row["minus_tau_r_R_squared_over_Js_lower_exact"] is None,
                "nonstrict producer has no witness")
    return stand,pd,overlap


def cold(checked):
    docs={p:load(p,h,checked,p.suffix==".json") for p,h in COLD.items()}
    prot=docs[PROD/"protocol.json"]
    pf=docs[PROD/"preflight.json"]
    review=docs[PROD/"independent_prelaunch_review.json"]
    require(prot["source_sha256"]==pf["source_sha256"]==PHASE["source_sha256"]
            and review["status"]=="PASS_INDEPENDENT_COLD_TRANSVERSE_SURFACE_STRESS_SCREEN"
            and all(review[k]==v for k,v in PHASE.items() if k!="review_sha256")
            and pf["status"]=="PASS" and pf["protocol_sha256"]==PHASE["protocol_sha256"]
            and prot["actual_stress_evaluations"]==pf["actual_stress_evaluations"]==0
            and len(prot["premise_sha256"])==26 and prot["total_windows"]==672
            and tuple(map(pair,prot["fixed_windows_on_each_ordered_edge"]))==WINDOWS,"current signed screen cold source and window list")
    load(PROD/"READBACK.md",prot["derivation_sha256"],checked,False)
    return prot,docs[OLD/"protocol.json"]


def prepare():
    checked={}
    p,k=cold(checked)
    points=[(Q(R,2),Q(3)*R/8) for R in (Q(1),Q(3,2),Q(2))]
    tr=trace(points,(Q(1),Q(2)),[(True,True)]*3,Q(0),Q(8),Q(2))
    _,r,s,D,j,tau,c,b=bound_row(tr,Q(8),(Q(1,4),Q(3,4)))
    require(j=="STRICT_NEGATIVE_PATCH" and Q(-7,2)<=tau<0 and 0<c<=Q(14,5),"independent tilted exact point test")
    cross=[(R*R/2,Q(0)) for R in (Q(1),Q(3,2),Q(2))]
    tr=trace(cross,(Q(1),Q(2)),[(True,True)]*3,Q(0),Q(4),Q(0))
    require(bound_row(tr,Q(4),(Q(1,4),Q(1,2)))[4]=="UNRESOLVED","sign crossing cannot pass")
    flat=[(R,Q(0)) for R in (Q(1),Q(3,2),Q(2))]
    tr=trace(flat,(Q(1),Q(2)),[(True,True)]*3,Q(0),Q(4),Q(0))
    require(bound_row(tr,Q(4),(Q(1,4),Q(3,4)))[3:5]==((0,0),"NONNEGATIVE_ON_WINDOW"),"equality not strict")
    tr=trace(points,(Q(1),Q(2)),[(True,True)]*3,Q(1,4096),Q(8),Q(2))
    require(bound_row(tr,Q(8),(Q(1,4),Q(3,4)))[4]=="STRICT_NEGATIVE_PATCH","proof with nonzero independent nodal uncertainty")
    d={"status":"PASS","source_sha256":sha(Path(__file__).read_bytes()),"method_sha256":sha((HERE/"METHOD.md").read_bytes()),
       "checks":{"cold_current_sources_and_windows":True,"exact_centered_power_intact":True,
                 "strict_negative_crossing_equality_and_positive_radius":True,"nonzero_posterior_trace_bound":True},
       "cold_actual_stress_windows_read_or_evaluated":0,"cold_bindings":checked}
    print(save(HERE/"independent_preflight.json",d),flush=True)


def authority(expected, checked):
    o=load(PROD/"actual/os_exit_observation.json",expected,checked)
    # Only root's actual observed record (not an internal terminal) authorizes
    # access to the three producer result files.
    tag=o.get("unified_exec_chunk_id")
    sess=o.get("unified_exec_session")
    require(o["status"]==OS_STATUS and type(o["pid"]) is int and o["pid"]>0
            and type(o["observed_exit_code"]) is int and o["observed_exit_code"]==0
            and (type(sess) is int and sess>0 or type(tag) is str and bool(tag)),
            "actual root-observed producer process exit required")
    start,term,r=(load(PROD/"actual"/(n+".json"),o[n+"_sha256"],checked) for n in ("started","terminal","results"))
    phase={**PHASE,"pid":o["pid"]}
    require(all(type(x["pid"]) is int and all(x[k]==v for k,v in phase.items()) for x in (start,term,r))
            and start["status"]=="RUNNING" and term["status"]=="FINISHED" and type(term["exit_code"]) is int
            and term["exit_code"]==0 and term["results_sha256"]==o["results_sha256"]
            and r["status"]=="CANDIDATE_TRANSVERSE_SURFACE_STRESS_AWAITS_INDEPENDENT_READBACK" and r["evaluated_windows"]==672
            and r["continuum_equilibrium_or_physical_wavelength_claimed"] is False
            and r["new_3D_solve_or_AD_or_spatial_integration"]==0,"same actual producer source and exactly 672 completed windows")
    return r


def pure_input_helper():
    require(sha(INPUT_HELPER.read_bytes())==INPUT_PIN,"exact already independently reviewed K/L proof input provenance")
    spec=spec_from_file_location("_independently_checked_K_L_root_state",INPUT_HELPER)
    m=module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def cases_input(p,k,checked):
    docs={h:load(ROOT/n,h,checked,n.endswith(".json")) for n,h in p["premise_sha256"].items()}
    h=pure_input_helper()
    cases=[]
    for c,q in zip(k["cases"],p["two_exact_spatial_roots"],strict=True):
        for key,v in c.items():
            require(q[key]==v,"same already independently checked K/L exact root target")
        u=h.actual_inputs(k,docs,checked,c)
        g=docs[q["geometry_proof_sha256"]]
        require(g["pilot_id"]==c["case_id"] and g["box_sha256"]==c["box_sha256"] and exact(g["radius_exact"])==h.RADIUS
                and g["all_144_body_and_24_trace_domains"] is True
                and exact(g["strong_monotonicity_lower_exact"])>0 and exact(g["min_trace_hoop_exact"])>0
                and exact(g["body_J_bounds_exact"][0])>Q(1,10**6),"same whole-cube geometry and derivative laws")
        indata=docs[c["cache_sha256"]]
        b=base64.b64decode(indata["frozen_binary_cache"]["data"],validate=True)
        with np.load(BytesIO(b),allow_pickle=False) as n:
            surfaces=[n[f"surface_{i}_edges"].tolist() for i in (0,1)]
        u["surfaces"]=surfaces
        cases.append(u)
    require(cases[0]["operator"]==cases[1]["operator"] and cases[0]["material"]==cases[1]["material"]
             and cases[0]["scale"]==cases[1]["scale"]>0 and cases[0]["gain"]==cases[1]["gain"]>0,"same exact fixed model")
    return cases


def compare_case(c,a,producer):
    require(all(producer[k]==v for k,v in a.items()),"same externally certified exact-spatial root metadata")
    result=[]
    n=0
    for face,edges in zip(FACES,c["surfaces"],strict=True):
        require(len(edges)==12 and all(x[2]==y[0] for x,y in zip(edges,edges[1:])),"complete ordered fixed surface")
        K,T=(pair(c["material"][face+name]) for name in ("_area_modulus","_surface_tension"))
        require(K[0]==K[1] and T[0]==T[1],"singleton nonnegative surface material")
        for j,(left,middle,right) in enumerate(edges):
            idx=(left,middle,right)
            ref=(c["reference"][left][0],c["reference"][right][0])
            data=trace([c["nodes"][i] for i in idx],ref,[c["free"][i] for i in idx],c["posterior"],K[0],T[0])
            for win in WINDOWS:
                raw=producer["rows"][n]
                R,r,s,D,own,tau,desc,used=bound_row(data,K[0],win)
                require(raw["face"]==face and type(raw["edge"]) is int and raw["edge"]==j and pair(raw["window_t_exact"])==win,
                        "all exact same-case ordered native windows")
                old,z,overlap=old_result(raw,R,r,s,D,data[0],K[0])
                result.append({"face":face,"edge":j,"window_t_exact":list(map(str,win)),"reference_R_exact":list(map(str,R)),
                    "independent_r_exact":list(map(str,r)),"independent_r_R_exact":list(map(str,s)),
                    "independent_D_exact":None if D is None else list(map(str,D)),
                    "producer_D_exact":None if z is None else list(map(str,z)),
                    "independent_tau_negative_upper_exact":None if tau is None else str(tau),
                    "independent_minus_tau_r_R_squared_over_Js_lower_exact":None if desc is None else str(desc),
                    "producer_standing":old,"independent_standing":own,"D_bounds_nondisjoint":overlap,
                    "independent_dyadic_subwindows":used})
                n+=1
    require(n==len(producer["rows"])==336,"24 surfaces x14 windows per case")
    return {"case_id":a["case_id"],"all_windows":result,"independent_counts":dict(Counter(x["independent_standing"] for x in result)),
         "producer_counts":dict(Counter(x["producer_standing"] for x in result)),
         "producer_strict_negative_not_independently_confirmed":sum(x["producer_standing"]=="STRICT_NEGATIVE_PATCH"
                      and x["independent_standing"]!="STRICT_NEGATIVE_PATCH" for x in result),
         "angular_instability_direction_exists_in_stated_enriched_natural_3D_model":
                      any(x["independent_standing"]=="STRICT_NEGATIVE_PATCH" for x in result)}


def main(args):
    checked={}
    pf=load(HERE/"independent_preflight.json",args.preflight_sha256,checked)
    require(pf["status"]=="PASS" and pf["source_sha256"]==sha(Path(__file__).read_bytes())
            and pf["method_sha256"]==sha((HERE/"METHOD.md").read_bytes()) and pf["cold_actual_stress_windows_read_or_evaluated"]==0,
            "frozen independent exact source, method and cold fixtures")
    p,k=cold(checked)
    require(all(checked[x]==h for x,h in pf["cold_bindings"].items()),"matched own cold source closure")
    stored=authority(args.observed_sha256,checked)
    require(stored["input_sha256"]==p["premise_sha256"],"same entire root and geometry proof premises")
    inputs=cases_input(p,k,checked)
    rows=[compare_case(v,a,r) for v,a,r in zip(inputs,p["two_exact_spatial_roots"],stored["cases"],strict=True)]
    require(all(stored["negative_windows_per_case"][q["case_id"]]==q["producer_counts"].get("STRICT_NEGATIVE_PATCH",0) for q in rows),
            "actual producer classifications summed exactly")
    h=save(HERE/"independent_actual_readback.json",{"status":"PASS_INDEPENDENT_ACTUAL_ALL672_TRANSVERSE_STRESS_WINDOWS",
         "independent_Preflight_sha256":args.preflight_sha256,"independent_source_sha256":sha(Path(__file__).read_bytes()),
         "verified_source_and_actual_sha256":checked,"genuine_producer_observed_exit_sha256":args.observed_sha256,
         "checked_windows":672,"case_readback":rows,"all_scalar_enclosures_nondisjoint":True,
         "no_nesting_or_intersection_used_as_proof":True,"physical_wavelength_or_continuum_stability_claimed":False,
         "actual_AD_matrix_solve_or_integration_jobs":0})
    print(json.dumps({"receipt_sha256":h,"independent_counts":{r["case_id"]:r["independent_counts"] for r in rows},
                     "confirmed_negative":{r["case_id"]:r["angular_instability_direction_exists_in_stated_enriched_natural_3D_model"] for r in rows}}),flush=True)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    sub=parser.add_subparsers(dest="phase",required=True)
    sub.add_parser("prepare")
    a=sub.add_parser("run")
    a.add_argument("--observed-sha256",required=True)
    a.add_argument("--preflight-sha256",required=True)
    args=parser.parse_args()
    prepare() if args.phase=="prepare" else main(args)
