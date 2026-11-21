"""Independent monomial replay for two sealed active-energy scalar endpoints."""
import argparse
import base64
from fractions import Fraction as Q
from hashlib import sha256
from io import BytesIO
import json
from math import factorial
from pathlib import Path
from zipfile import ZipFile

import numpy as np

HERE=Path(__file__).resolve().parent
RUN=HERE.parent
EX=RUN.parent
ROOT=next(p for p in EX.parents if (p/"PLAN.md").is_file())
PINS={
    "active_energy.py":"072ed83a30565be34bd5fe7dde559526e089a070df0ad95f0e61611db9e9e76d",
    "protocol.json":"06dd35434358ad7dcb993349281ed6bf8a84fd0a34f7ca2d4932ec0792842361",
    "preflight.json":"ae11e52676f1e7ecacc5b925de47a553024dec6604d4302ade9ddd2bad43738c",
    "independent_prelaunch_review.json":"7cb39b176a8221029dd101b6c30437c54ad238dd7b24928c2609a718d9b2e344",
    "actual/os_exit_observation.json":"0087a7f44dd09188062bed3b621dc1460b72df64ce4640a9f36e682903362f0f",
    "actual/results.json":"18a94a0fed94471b00a32fc0617581fe3b6b603581f63d442e5b584b96f8ef48",
}
PHASE={k+"_sha256":PINS[p] for k,p in (("source","active_energy.py"),("protocol","protocol.json"),
    ("preflight","preflight.json"),("review","independent_prelaunch_review.json"))}
CASES=("max_negative_contrast_change_small_16","max_negative_contrast_change_small_17")
RADIUS=Q(1,2**34)
VALUE="unit_drive_energy_over_normalization_exact"
OUT_DIFF="V_L_minus_V_K_over_shared_normalization_exact"
EDGES=((0,1),(1,2),(2,0))
EXPONENTS=((0,0),(1,0),(0,1),(2,0),(1,1),(0,2))
P2=tuple({p:Q(a) for p,a in zip(EXPONENTS,row,strict=True) if a} for row in
    ((1,-3,-3,2,4,2),(0,-1,0,2,0,0),(0,0,-1,0,0,2),(0,4,0,-4,-4,0),(0,0,0,0,4,0),(0,0,4,0,-4,-4)))
VERTEX=((0,0),(1,0),(0,1))


def require(ok,msg):
    if not ok:
        raise ValueError(msg)


def sha(b):
    return sha256(b).hexdigest()


def load(p,h,closed,parse=True):
    p=p.resolve()
    name=p.relative_to(ROOT).as_posix()
    b=p.read_bytes()
    require(sha(b)==h==closed.setdefault(name,h),"invalid closed authority "+name)
    return json.loads(b) if parse else b


def write(p,d):
    with p.open("x",encoding="utf-8",newline="\n") as f:
        json.dump(d,f,indent=2)
        f.write("\n")
    return sha(p.read_bytes())


def mul(p,q):
    terms={}
    for (i,j),a in p.items():
        for (k,l),b in q.items():
            key=(i+k,j+l)
            terms[key]=terms.get(key,0)+a*b
    return terms


def add(p,q,c=Q(1)):
    terms=dict(p)
    for k,v in q.items():
        terms[k]=terms.get(k,0)+c*v
    return terms


def derivative(p,g0,g1):
    out={}
    for (s,t),x in p.items():
        if s:
            key=(s-1,t)
            out[key]=out.get(key,0)+s*x*g0
        if t:
            key=(s,t-1)
            out[key]=out.get(key,0)+t*x*g1
    return out


def monomial_integral(p):
    return sum((a*Q(factorial(i)*factorial(j),factorial(i+j+2)) for (i,j),a in p.items()),Q(0))


def at(p,x):
    s,t=x
    return sum(v*s**i*t**j for (i,j),v in p.items())


def active_element(ref,y,move,e):
    r0,z0=ref[0]
    drs,dzs=(b-a for a,b in zip(ref[0],ref[1],strict=True))
    drt,dzt=(b-a for a,b in zip(ref[0],ref[2],strict=True))
    area=drs*dzt-drt*dzs
    require(area!=0 and min(r for r,z in ref)>=0 and max(r for r,z in ref)>0,"positive material radius and nonzero element map")
    radii={(0,0):r0,(1,0):drs,(0,1):drt}
    dz=[derivative(p,-drt/area,drs/area) for p in P2]
    energy=Q(0)
    perturbation=Q(0)
    max_center=[]
    max_change=[]
    for j in (0,1):
        q={}
        for ni,x in zip(dz,y,strict=True):
            q=add(q,ni,x[j])
        center=max(abs(at(q,p)) for p in VERTEX)
        delta=e*max(sum(abs(at(ni,p)) for ni,f in zip(dz,move,strict=True) if f[j]) for p in VERTEX)
        energy+=abs(area)*monomial_integral(mul(radii,mul(q,q)))/2
        perturbation+=abs(area)*monomial_integral(radii)*(center*delta+delta**2/2)
        max_center.append(center)
        max_change.append(delta)
    require(energy>=0 and perturbation>=0,"energy nonnegative sum of per-axis squares")
    return energy,perturbation,(max(Q(0),energy-perturbation),energy+perturbation),max_center,max_change


def pair(seq):
    require(type(seq) is list and len(seq)==2 and all(type(x) is str for x in seq),"explicit ordered exact rational endpoints")
    lo,hi=map(Q,seq)
    require(lo<=hi,"reversed interval")
    return lo,hi


def overlap(a,b):
    return max(a[0],b[0])<=min(a[1],b[1])


def interval_compare(a,b):
    return {"independent_interval_exact":list(map(str,a)),"producer_interval_exact":list(map(str,b)),
            "independent_contains_producer":a[0]<=b[0] and b[1]<=a[1],
            "producer_contains_independent":b[0]<=a[0] and a[1]<=b[1],"nonempty_intersection":overlap(a,b)}


def toy():
    ref=[(Q(0),Q(0)),(Q(2),Q(0)),(Q(0),Q(3))]
    pts=ref+[tuple((p+q)/2 for p,q in zip(ref[i],ref[j],strict=True)) for i,j in EDGES]
    data=[(z*z,r*z) for r,z in pts]
    dirs=[(j!=1,j!=0) for j in range(6)]
    e=Q(1,37)
    c,d,interval,*_=active_element(ref,data,dirs,e)
    require(c==Q(24,5) and d>0,"center polynomial oracle")
    for sign in (-1,1):
        varied=[tuple(v+sign*e if dirs[j][k] else v for k,v in enumerate(row)) for j,row in enumerate(data)]
        proof=active_element(ref,varied,dirs,Q(0))[0]
        require(interval[0]<=proof<=interval[1],"perturbation witness")
    no_move=[(False,False)]*6
    require(active_element(ref,data,no_move,e)[1]==0,"fixed coordinates receive no error")
    p=write(HERE/"preflight.json",{"status":"PASS","source_sha256":sha(Path(__file__).read_bytes()),
        "method_sha256":sha((HERE/"METHOD.md").read_bytes()),"checks":{"exact_center_known_axis_polynomial":True,
        "exact_fixed_and_moving_coordinate_error_bounds":True},"actual_energy_evaluations":0})
    print(p,flush=True)


def require_scope(obs,phase,start,end,res):
    require(obs["status"]=="PASS_ROOT_OBSERVED_NEGATIVE_KL_ACTIVE_ENERGY_EXIT" and obs["pid"]==1996
            and obs["unified_exec_chunk_id"]=="a76a43" and type(obs["observed_exit_code"]) is int
            and obs["observed_exit_code"]==0 and obs["direct_root_tool_observation"] is True
            and obs["independent_actual_arithmetic_readback_required"] is True,"observed genuine producer exit required")
    require(start["status"]=="RUNNING" and end["status"]=="FINISHED" and type(end["exit_code"]) is int and end["exit_code"]==0
            and end["results_sha256"]==obs["results_sha256"]==PINS["actual/results.json"]
            and all(type(x["pid"]) is int and x["pid"]==obs["pid"] and all(x[k]==v for k,v in phase.items()) for x in (obs,start,end,res))
            and res["status"]=="CANDIDATE_ACTIVE_ENERGY_ORDER_AWAITS_INDEPENDENT_READBACK"
            and res["actual_energy_evaluations"]==2 and res["certified_loading_path"] is False
            and res["source_fit_or_biological_dynamics_claimed"] is False,"complete same frozen producer candidate")


def actual_inputs(prot,docs,closed,row):
    contract=docs[row["contract_sha256"]]
    packet=docs[row["cache_sha256"]]
    box=docs[row["box_sha256"]]
    native=docs[row["exact_spatial_proof_sha256"]]
    case=row["case_id"]
    match=case==CASES[0]
    target=native["target"] if match else native["target_identity"]
    a=native["arithmetic"] if match else native["independent_bounds"]
    eta,z,e=map(Q,(a["eta_upper_exact"],a["z_upper_exact"],a["posterior_error_upper_exact"]))
    r=Q(native["radius_exact"])
    require(r==RADIUS and a["strict_selfmap"] is True and 0<=z<1 and eta>=0
            and e==eta/(1-z)==Q(row["posterior_exact"]) and e<r and eta+z*r==Q(a["selfmap_lhs_upper_exact"])<r,
            "own exact posterior from current native stationarity proof")
    if match:
        h=native["observer_and_optional_minimum"]["normalized_H_euclidean_lower_exact"]
        require(native["status"]=="PASS_INDEPENDENT_FIRST_WAVE_EXACT_SPATIAL_ROOT"
                and target["archived_center_shift_exact"]!="0" and target["trial_box_sha256"]==row["box_sha256"],"K genuine recentered proof")
    else:
        h=native["whole_box_normalized_H_euclidean_lower_exact"]
        require(native["status"]=="PASS_INDEPENDENT_EXACT_SPATIAL_ROOT_AND_STRICT_LOCAL_MINIMUM","L genuine original proof")
        need=[name for name,digest in prot["inputs_sha256"].items() if digest==row["box_sha256"]]
        require(len(need)==1 and native["source_and_actual_sha256"][need[0]]==row["box_sha256"],"exact L proof cube pin")
    require(Q(h)>0 and contract["pilot_id"]==case and
            target["model_contract_sha256"]==row["contract_sha256"]==box["model_observation_contract_sha256"]
            and target["record_id"]==contract["record_id"]==box["record_id"]
            and target.get("point_cache_sha256",target.get("cache_packet_sha256"))==row["cache_sha256"]
            ==contract["point_cache_sha256"]==box["point_cache_sha256"],"same finite root and proof operator identification")
    gain,dose,normal,drive=(Q(row[k]) for k in ("gain_exact","dose_exact","normalization_exact","exact_drive_exact"))
    require(drive==gain*dose==Q(target.get("exact_drive",target.get("drive_exact")))
            ==Q(packet["activation_strength_exact"]) and dose==Q(box["dose_exact"])>0 and gain>0 and
            Q(contract["energy_scale_binary_exact"])==Q(packet["normalization_binary_exact"])==normal>0
            and Q(contract["material_parameters_binary_exact"]["reference_curvature"])==0,"fixed flat law and already shared normalization")
    data=base64.b64decode(packet["frozen_binary_cache"]["data"],validate=True)
    require(sha(data)==packet["frozen_binary_cache"]["sha256"]==contract["frozen_binary_cache_sha256"],"exact binary state archive")
    with ZipFile(BytesIO(data)) as archive:
        operator={n:sha(archive.read(n)) for n in archive.namelist() if n not in ("u.npy","preconditioner.npy")}
    with np.load(BytesIO(data),allow_pickle=False) as n:
        names=("reference","u","free","triangles","active","energy_scale","surface_0_edges","surface_1_edges")
        u={k:n[k] for k in names}
    for k,h in prot["same_model_array_bytes_sha256"].items():
        require(sha(u[k].tobytes())==h,"same fixed material operator array "+k)
    require(len(operator)==21 and u["reference"].shape==(325,2) and u["u"].shape==(611,)
            and u["triangles"].shape==(144,6) and u["free"].shape==(611,) and u["active"].shape==(144,1)
            and Q(float(u["energy_scale"][0]))==normal,"case exact target dimensions")
    X=[list(map(Q,x)) for x in u["reference"].tolist()]
    flat=[x for row in X for x in row]
    rim=Q(contract["material_parameters_binary_exact"]["outer_radius"])
    free=[2*i+k for i,(R,Z) in enumerate(X) for k in (0,1) if R!=rim and (R!=0 or k==1)]
    require(free==u["free"].tolist()==box["free"] and len(free)==611
            and [[Q(x) for x in row] for row in box["reference_nodes_exact"]]==X
            and box["triangles"]==u["triangles"].tolist(),"same exact axis/rim map and reference coordinates")
    shift={j:Q(v) for j,v in zip(free,u["u"].tolist(),strict=True)}
    Y=[x+shift.get(j,0) for j,x in enumerate(flat)]
    for j,(y,w,p) in enumerate(zip((v for c in box["current_centers_exact"] for v in c),
             (v for c in box["coordinate_radii_exact"] for v in c),(v for c in box["coordinate_hulls_exact"] for v in c),strict=True)):
        radius=r if j in shift else Q(0)
        require(Y[j]==Q(y) and Q(w)==radius and pair(p)==(Y[j]-radius,Y[j]+radius),"same exact proof center and literal free/fixed coordinate hull")
    nodes=list(zip(Y[::2],Y[1::2],strict=True))
    moving=[(2*i in shift,2*i+1 in shift) for i in range(len(X))]
    mask=[int(m[0]) for m in u["active"].tolist()]
    require(set(mask)=={0,1} and sum(mask)==24,"24 active original reference elements")
    return {"nodes":nodes,"free":moving,"reference":X,"triangles":u["triangles"].tolist(),"mask":mask,
            "posterior":e,"dose":dose,"gain":gain,"material":contract["parameter_cell"],"scale":normal,"operator":operator,
            "center_u_sha256":sha(u["u"].tobytes())}


def main(preflight_sha):
    closed={}
    pre=load(HERE/"preflight.json",preflight_sha,closed)
    require(pre["status"]=="PASS" and pre["source_sha256"]==sha(Path(__file__).read_bytes())
            and pre["method_sha256"]==sha((HERE/"METHOD.md").read_bytes()) and pre["actual_energy_evaluations"]==0,
            "own frozen replay source and method")
    files={p:load(RUN/p,h,closed,p.endswith(".json")) for p,h in PINS.items()}
    p=files["protocol.json"]
    o=files["actual/os_exit_observation.json"]
    results=files["actual/results.json"]
    s,t=(load(RUN/"actual"/(x+".json"),o[x+"_sha256"],closed) for x in ("started","terminal"))
    require_scope(o,PHASE,s,t,results)
    require(len(p["inputs_sha256"])==22 and p["cold_actual_energy_evaluations"]==0
            and p["same_effective_binary_material_and_exact_gain"] is True,"original 22 frozen evidence pins")
    docs={h:load(ROOT/k,h,closed,k.endswith(".json")) for k,h in p["inputs_sha256"].items()}
    require(tuple(x["case_id"] for x in p["cases"])==CASES,"exact ordered K–L")
    pk,pl=(actual_inputs(p,docs,closed,x) for x in p["cases"])
    require(pk["operator"]==pl["operator"] and pk["material"]==pl["material"] and pk["scale"]==pl["scale"] and
            pk["gain"]==pl["gain"] and pl["dose"]-pk["dose"]==Q(22159027,10**8),"same physical operator increasing dose")
    values=[]
    cases=[]
    for frozen,actual,case in zip(p["cases"],results["two_certified_same_operator_cases"],(pk,pl),strict=True):
        require(all(actual[k]==v for k,v in frozen.items()),"same both producer and independent root inputs")
        records=[]
        total_center,total_error=Q(0),Q(0)
        counts=[]
        for i,(flag,tri) in enumerate(zip(case["mask"],case["triangles"],strict=True)):
            if not flag:
                continue
            require(len(set(tri))==6 and all(0<=j<325 for j in tri),"valid exact P2 element connectivity")
            c,d,intv,M,E=active_element([case["reference"][j] for j in tri[:3]],[case["nodes"][j] for j in tri],
                                       [case["free"][j] for j in tri],case["posterior"])
            slot=actual["active_elements"][len(records)]
            require(slot["element"]==i,"same active element index")
            corr=interval_compare(intv,pair(slot[VALUE]))
            require(corr["nonempty_intersection"],"valid same element scalar enclosures disjoint")
            counts.append(corr)
            records.append({"element":i,"center_exact":str(c),"perturbation_exact":str(d),
                "derivative_abs_center_max_exact":list(map(str,M)),"derivative_perturbation_max_exact":list(map(str,E)),**corr})
            total_center+=c
            total_error+=d
        require(len(records)==len(actual["active_elements"])==24,"same complete active support")
        ours=(max(Q(0),total_center-total_error),total_center+total_error)
        other=pair(actual[VALUE])
        require(other==tuple(sum(pair(x[VALUE])[k] for x in actual["active_elements"]) for k in (0,1)),
                "stored producer total agrees with its sealed exact local packets")
        corr=interval_compare(ours,other)
        require(corr["nonempty_intersection"],"whole energy comparison disjoint")
        cases.append({**frozen,"independent_center_energy_exact":str(total_center),"independent_perturbation_error_exact":str(total_error),
          "center_binary_sha256":case["center_u_sha256"],"comparison":corr,"all24_active_elements":records,
          "disjoint_elements":0,"nonnested_elements":sum(not x["independent_contains_producer"] and not x["producer_contains_independent"] for x in counts)})
        values.append((ours,other))
    ownK,prK=values[0]
    ownL,prL=values[1]
    ours=(ownL[0]-ownK[1],ownL[1]-ownK[0])
    old=(prL[0]-prK[1],prL[1]-prK[0])
    out=interval_compare(ours,old)
    require(ours[1]<0 and old==pair(results[OUT_DIFF]) and out["nonempty_intersection"]
            and results["classification"]=="NECESSARY_ENERGY_ORDER_VERIFIED_ONLY","independent strict negative signed difference")
    h=write(HERE/"receipt.json",{"status":"PASS_INDEPENDENT_ACTUAL_NEGATIVE_KL_NECESSARY_ENERGY_ORDER",
       "independent_replay_source_sha256":sha(Path(__file__).read_bytes()),"preflight_sha256":preflight_sha,
       "verified_input_sha256":closed,"two_endpoint_proofs":cases,"difference":out,
       "independent_normalized_difference_divided_by_common_factor_display":list(map(float,ours)),
       "case_count":2,"active_elements_each":24,"disjointness_is_diagnostic_only":True,
       "positive_common_factor_2pi_over_scale":True,"successful_necessary_order_only":True,
       "certified_loading_path":False,"source_fit_or_biological_dynamics_claimed":False,
       "new_global_matrix_AD_or_equilibrium_solves":0})
    print(json.dumps({"receipt_sha256":h,"own_difference":list(map(float,ours)),
          "own_endpoints":[list(map(float,v[0])) for v in values]}),flush=True)


if __name__=="__main__":
    p=argparse.ArgumentParser()
    sub=p.add_subparsers(dest="mode",required=True)
    sub.add_parser("prepare")
    a=sub.add_parser("run")
    a.add_argument("--preflight-sha256",required=True)
    args=p.parse_args()
    toy() if args.mode=="prepare" else main(args.preflight_sha256)
