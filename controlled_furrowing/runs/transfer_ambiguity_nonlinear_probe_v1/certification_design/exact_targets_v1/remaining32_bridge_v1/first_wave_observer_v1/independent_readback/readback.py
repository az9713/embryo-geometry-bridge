"""Independent rational archived-geometry readback of first-wave spatial observers."""
import argparse
from collections import Counter
from fractions import Fraction as Q
from hashlib import sha256
import json
import os
from pathlib import Path
import traceback

HERE=Path(__file__).resolve().parent
OB=HERE.parent
BR=OB.parent
EX=BR.parent
DESIGN=EX.parent
PROBE=DESIGN.parent
ROOT=PROBE.parents[1]
SRC=Path(__file__).resolve()
NATIVE=BR/"first_wave_reduction_v1/independent_actual_preparation/actual_replay"
GEOM=EX/"remaining32_geometry_v1/preparation/protocol.json"
BUD=DESIGN/"budgets_v1"
OLD=DESIGN/"full_outcome_coverage_v1/observer_application_v1"
SEC=PROBE/"observer_robustness"
PRIMARY=PROBE/"independent_actual_preparation/actual_replay/receipt.json"
R=Q(1,2**34)
L=Q(67,2)
CASES=("max_contrast_change_small_01","max_contrast_change_small_13",
       "max_negative_contrast_change_small_01","max_negative_contrast_change_small_16")
VAR=("current_disk__spatial_max","current_disk__regional_mean",
     "current_radial_line__spatial_max","current_radial_line__regional_mean")
OBS=("maximum_shortening_um","region_averaged_shortening_um","apical_max_shortening_um","basal_max_shortening_um")
PASS_ROOT="PASS_INDEPENDENT_FIRST_WAVE_EXACT_SPATIAL_ROOT"
PINS={
 OB/"apply.py":"6db0999d7a11c6bcc7fbd4030a890c7ea80dcca0a41264f40329fb7070264b37",
 OB/"protocol.json":"a8852b364005c87b81f54b466515ebae460f4c9ebac66bd1063aa62764b36b9f",
 OB/"preflight.json":"e9ba72a948ae421e9b3923049bae4e9e64300be3272760398f90ff5dcbb6368d",
 OB/"independent_source_review/receipt.json":"9876af4b080d59c78045c0d1b815b6dc7374438efc751d8352991cd40d7125cc",
 OB/"actual/results.json":"4213df2fc7af939b0d44e27e7982829c1f7f28ec0a521bb7a10117516a72d761",
 OB/"actual/terminal.json":"76bc88f692a7a62a2887b626321f9f25eaafdf41ec8d4b81a468c4b869d6b34d",
 OB/"actual/started.json":"de1720519c5963e614a987386131cbe06e1919c15a5802b69dce005877c3f7fe",
 NATIVE/"os_exit_observation.json":"aa3c90101c64ccefe1403af2c296bcb22cafad3b2e35d73539aa7e3c5227f281",
 NATIVE/"receipt.json":"e73ea05a4703256235a46ce7b914298d3f8258e9fef029f4dbe952b04bb9798a",
 GEOM:"754da64bb96cd8a19707181315e6c3e3f014e8d37a898b83c1693342c82c9850",
 PRIMARY:"c795ae78ba0b53a01c3d65f466d5fa2f1d233c2a68c4275c04260c3cb5120ebd",
 BUD/"actual/results.json":"689948766db950ee4ddf399e6731aea08890cc7a1f84c5615d19e7ded68b0872",
 BUD/"independent_actual/receipt.json":"978a2a4095bb8735103a5ee2c11ad266a6b0e052e22af2b7a36ef29f1f38265c",
 OLD/"results.json":"1df2bc372edd12687fad17edf63013509994e6bbe41c02150cefb25952eab9f6",
 OLD/"independent_readback/receipt.json":"a06549f048abba4eb2b4d1c1fb398f06e8f524196c92a673bd0493c527c8f40d",
 SEC/"actual/results.json":"53a17296ddae08b4d1fe6b7e9b20ad0a377bf8dc962142bbfb9eb9877a381647",
 SEC/"independent_actual/receipt.json":"703b574bc51a083b10c7328d9622862e939ddabe791b18ce9e1817f07718f69f",
 ROOT/"docs/NODAL-OBSERVER-ERROR-TRANSFER.md":"d0d5f711a6a8f2ebe46e8dcb76c70234753bda714836dbcf7d4be888ff7ba171",
}
# Producer and independently observed phase details are confirmed outside its evaluator.
ROOT_OS_STATUS="PASS_ROOT_OBSERVED_FIRST_WAVE_EXACT_OBSERVER_APPLICATION_EXIT"


def must(ok,message):
    if not ok:
        raise ValueError(message)


def rel(path):
    return path.resolve().relative_to(ROOT).as_posix()


def sha(path):
    return sha256(path.read_bytes()).hexdigest()


def doc(path,h=None,links=None):
    path=(ROOT/path).resolve()
    h=PINS[path] if h is None else h
    b=path.read_bytes()
    must(path.is_relative_to(ROOT) and sha256(b).hexdigest()==h,"same immutable file "+str(path))
    if links is not None:
        must(links.setdefault(rel(path),h)==h,"conflicting phase file hash")
    return json.loads(b)


def save(path,data):
    with path.open("x",encoding="utf8",newline="\n") as f:
        json.dump(data,f,indent=2,allow_nan=False)
        f.write("\n")
    return sha(path)


def q(x):
    must(type(x) in (int,str,Q),"literal exact rational")
    return Q(x)


def b(x):
    must(type(x) in (int,float),"exact binary coordinate")
    return Q(x)


def pair(x):
    must(type(x) in (list,tuple) and len(x)==2,"complete exact pair")
    lo,hi=map(q,x)
    must(lo<=hi,"ordered interval")
    return lo,hi


def pick(rows,**keys):
    matched=[r for r in rows if all(r[k]==v for k,v in keys.items())]
    must(len(matched)==1,"unique matched original field")
    return matched[0]


def native_rows(links):
    obs=doc(NATIVE/"os_exit_observation.json",links=links)
    r=doc(NATIVE/"receipt.json",obs["receipt_sha256"],links)
    term=doc(NATIVE/"terminal.json",obs["terminal_sha256"],links)
    must(obs["status"]=="PASS_ROOT_OBSERVED_FIRST_WAVE_INDEPENDENT_GLOBAL_EXIT" and obs["direct_root_tool_observation"] is True
         and type(obs["pid"]) is type(term["pid"]) is type(r["pid"]) is int and obs["pid"]==term["pid"]==r["pid"]>0
         and type(obs["observed_exit_code"]) is type(term["exit_code"]) is int and obs["observed_exit_code"]==term["exit_code"]==0
         and term["status"]=="FINISHED" and term["receipt_sha256"]==obs["receipt_sha256"]
         and r["status"]=="COMPLETE_INDEPENDENT_FIRST_WAVE_SIGNED_ENDPOINT_REPLAY"
         and r["exact_spatial_roots"]==4 and tuple(x["case_id"] for x in r["case_results"])==CASES,
         "genuine four independently completed root purposes")
    must(obs["native_cases_sha256"]=={r["native_file"]:r["native_sha256"] for r in r["case_results"]},
         "native archived root leaf SHA set")
    return r["case_results"],term


def original_sources(links):
    table=doc(BUD/"actual/results.json",links=links)
    audit=doc(BUD/"independent_actual/receipt.json",links=links)
    old=doc(OLD/"results.json",links=links)
    old_a=doc(OLD/"independent_readback/receipt.json",links=links)
    secondary=doc(SEC/"actual/results.json",links=links)
    secondary_a=doc(SEC/"independent_actual/receipt.json",links=links)
    primary=doc(PRIMARY,links=links)
    geometry=doc(GEOM,links=links)
    allow=q(table["conditional_source_allowance_um_exact"])
    must(allow==q(old["frozen_allowance_um_exact"])==q(old_a["frozen_model_discrepancy_um_exact"])
         ==q(secondary["frozen_primary_source_allowance_um_exact"])==q(secondary_a["source_configured_delta_um_exact"])>0,
         "one retained allowance, no rescore")
    must(audit["status"]=="PASS_INDEPENDENT_ACTUAL_CONDITIONAL_OBSERVER_BUDGETS"
         and old_a["status"]=="PASS_INDEPENDENT_ACTUAL_36_CACHED_ROOT_OBSERVERS"
         and secondary_a["status"]=="PASS_INDEPENDENT_SECONDARY_STORED_OBSERVER_SWEEP"
         and audit["source_and_input_hashes"][rel(BUD/"actual/results.json")]==PINS[BUD/"actual/results.json"]
         and old_a["current_input_sha256"][rel(OLD/"results.json")]==PINS[OLD/"results.json"]
         and secondary_a["input_and_output_sha256"][rel(SEC/"actual/results.json")]==PINS[SEC/"actual/results.json"],
         "same separately reviewed immutable budget and stored observer proofs")
    return table,audit,old,old_a,secondary,secondary_a,primary,geometry,allow


def independent_root(leaf):
    a=leaf["arithmetic"]
    t=leaf["target"]
    h=q(a["b_H_upper_exact"])
    z=q(a["z_upper_exact"])
    eta=q(a["eta_upper_exact"])
    rf=list(map(pair,a["signed_R_F_I_rows_exact"]))
    rh=list(map(q,a["RH_absolute_row_bounds_exact"]))
    must(len(rf)==len(rh)==611 and min(rh)>=0 and max(rh)==h
         and eta==max(max(map(abs,e)) for e in rf) and z==q(a["ZQ_upper_exact"])+h
         and 0<=z<1 and eta+z*R==q(a["selfmap_lhs_upper_exact"])<R
         and a["strict_selfmap"] is True and a["strict_cube_exclusion"] is False,
         "native signed scalar root bounds reproduced")
    e=eta/(1-z)
    shift,xi,goal=(q(t[k]) for k in ("archived_center_shift_exact","original_addition_gap_exact","observer_goal_exact"))
    v=leaf["observer_and_optional_minimum"]
    eps=L*(e+shift+xi)
    must(0<=e==q(a["posterior_error_upper_exact"])<R and xi>=0 and 0<=e+shift<goal
         and e+shift==q(v["archived_shift_plus_posterior"]) and v["archived_state_budget_proved"] is True
         and eps==q(v["physical_nodal_error_um_upper_exact"]),"original observer budget exact displacement")
    return e,shift,xi,eps

def archived_geometry(state,box,budget,prior,saved,eps,shift,xi):
    ref=[tuple(map(b,node)) for node in state["reference_nodes"]]
    raw=[tuple(map(b,node)) for node in state["current_nodes"]]
    flat=[v for node in ref for v in node]
    saved_flat=[v for node in raw for v in node]
    rim=b(state["parameters"]["outer_radius"])
    height=b(state["parameters"]["height"])
    free=[2*i+j for i,(r,_) in enumerate(ref) for j in (0,1) if r!=rim and (r!=0 or j==1)]
    fixed=set(range(650))-set(free)
    model=budget["center_and_model"]
    must(len(ref)==len(raw)==325 and len(free)==len(state["free_displacement"])==611
         and free==box["free"]==model["free_map"] and state["triangles"]==box["triangles"]
         and len(box["triangles"])==144 and state["parameters"]==model["material_parameters"],
         "independently rebuilt ordered fixed boundary elimination")
    exact_arch=flat.copy()
    for k,v in zip(free,state["free_displacement"],strict=True):
        exact_arch[k]+=b(v)
    center=[q(x) for y in box["current_centers_exact"] for x in y]
    radii=[q(x) for y in box["coordinate_radii_exact"] for x in y]
    ends=[pair(x) for y in box["coordinate_hulls_exact"] for x in y]
    must(len(center)==len(radii)==len(ends)==650 and
         [q(x) for y in box["reference_nodes_exact"] for x in y]==flat
         and max(abs(p-s) for p,s in zip(exact_arch,saved_flat,strict=True))==xi
         and all(Q.from_float(float(a))==s for a,s in zip(exact_arch,saved_flat,strict=True))
         and max(abs(p-g) for p,g in zip(exact_arch,center,strict=True))==shift
         and all(center[i]==exact_arch[i]==saved_flat[i] for i in fixed)
         and all(radius==(R if i in free else 0) and ep==(coord-radius,coord+radius)
                 for i,(coord,radius,ep) in enumerate(zip(center,radii,ends,strict=True))),
         "exact archived addition gap, proof center shift and 650 trial endpoint identities")
    packed=saved["observer_contract"]
    def rational_object(v):
        return Q(v["numerator"],v["denominator"])
    roi=(Q(0),L*b(state["parameters"]["patch_radius"]))
    must(roi==pair(budget["fixed_roi_um_exact"])==pair(prior["ROI_um_exact"])
         ==tuple(map(rational_object,packed["current_radius_ROI_um"]))
         and 0<roi[1]<=L*rim and rational_object(packed["reference_apical_um"])==q(model["reference_height_um_exact"])==L*height==L
         and rational_object(packed["reference_basal_um"])==0,"fixed nonempty shared dimensional ROI and reference")
    graph={}
    for face,side,z in (("top","apical",height),("bottom","basal",Q(0))):
        unique=set()
        # Recover physical boundary edge orientation from the mesh itself.
        for t in state["triangles"]:
            for i,j,k in ((0,3,1),(1,4,2),(2,5,0)):
                nodes=[t[m] for m in (i,j,k)]
                if all(ref[n][1]==z for n in nodes):
                    unique.add(tuple(sorted(nodes,key=lambda n:ref[n][0])))
        edges=sorted(unique,key=lambda t:ref[t[0]][0])
        base=model["fixed_trace_endpoints"][face]
        must(len(edges)==12 and list(map(list,edges))==base["trace_nodes"]
             and all(e[-1]==f[0] for e,f in zip(edges[:-1],edges[1:],strict=True))
             and [2*edges[0][0],2*edges[-1][-1]]==base["radial_dofs"]
             and all(k in fixed for k in base["radial_dofs"])
             and tuple(L*center[k] for k in base["radial_dofs"])==(0,rim*L)
             ==pair(base["proved_fixed_radii_um_exact"]),"12 conforming edges with pinned radial axis and rim")
        controls=[[tuple(L*x for x in raw[i]) for i in edge] for edge in edges]
        must(controls==[[tuple(L*b(x) for x in row) for row in edge]
                     for edge in state["boundary_"+face+"_quadratic_controls"]],"same archived P2 nodal trace")
        lower=[]
        slopes=[]
        for points in controls:
            # From P(t)=A*t^2+B*t+C using the 0,1/2,1 interpolation equations.
            aa=[2*(m+q_-2*n) for m,n,q_ in zip(*points,strict=True)]
            bb=[q_-m-a for m,n,q_,a in zip(*points,aa,strict=True)]
            endpoints=((bb[0],bb[1]),(2*aa[0]+bb[0],2*aa[1]+bb[1]))
            rmin=min(x[0] for x in endpoints)
            must(rmin>8*eps,"strict radial derivative lower bound including root perturbation")
            lower.append(rmin)
            slopes.append(max(abs(x[1]) for x in endpoints)/rmin)
        m,zcap=min(lower),max(slopes)
        old_eps=q(prior["root"]["epsilon_um"])
        old_side=prior["apical_and_basal_graph_error_bounds"][side]
        expected=budget["stored_face_derivatives"][face]
        must(m==q(expected["radial_min_um_exact"])==q(old_side["radial_derivative_margin_exact"])+8*old_eps
             and zcap==q(expected["global_stored_slope_cap_exact"])==q(old_side["base_graph_Lipschitz_upper_exact"])
             and q(old_side["Eulerian_vertical_error_um_exact"])==Q(5,4)*(1+zcap)*old_eps,
             "independently recreated archived per-side radial and slope constants")
        graph[side]={"stored_radial_min_um_exact":str(m),"slope_cap_exact":str(zcap),
                     "strict_radial_margin_um_exact":str(m-8*eps),"graph_error_um_exact":str(Q(5,4)*eps*(1+zcap))}
    return roi,graph


def interval_check(stored,sample,err,allow):
    s0,s1=pair(stored)
    must(err>=0 and allow>=0,"positive error and applicable fixed allowance")
    z0,z1=s0-err,s1+err
    a,b=z0-sample,z1-sample
    lo=Q(0) if a<=0<=b else min(abs(a),abs(b))
    hi=max(abs(a),abs(b))
    where="IN" if hi<=allow else "OUT" if lo>allow else "UNRESOLVED"
    return (z0,z1),(lo,hi),where


def read_case(spec,phase,sources,reported,links):
    name=spec["case_id"]
    endpoint,number=name.rsplit("_small_",1)
    checkpoint="small_"+number
    table,audit,old,old_a,secondary,secondary_a,primary,geom,allow=sources
    leaf=doc(spec["native_file"],spec["native_sha256"],links)
    tgt=leaf["target"]
    trial=pick(geom["trials"],case_id=name)
    ctr=doc(trial["model_contract"],trial["model_contract_sha256"],links)
    bdoc=pick(table["field_budgets"],endpoint_id=endpoint,checkpoint=checkpoint)
    bd=pick(audit["fields"],endpoint_id=endpoint,checkpoint=checkpoint)
    prior=pick(old["cases"],case_id=name)
    file=rel(SEC/"actual"/endpoint/(checkpoint+".json"))
    sp=pick(secondary["field_observer_receipts"],path=file)["sha256"]
    saved=doc(file,sp,links)
    must(spec["status"]==leaf["status"]==PASS_ROOT and spec["observer_state_budget_proved"] is True
         and all(leaf[k]==phase[k] for k in ("source_sha256","protocol_sha256","preflight_sha256","pid",
                  "producer_root_terminal_admission_sha256","producer_reduction_input_admission_sha256"))
         and leaf["closed_cube_invertible_R_from_contraction"] is True and leaf["original_observer_field_replaced"] is False
         and leaf["case_id"]==tgt["case_id"]==name and tgt["record_id"]==trial["record_id"]==ctr["record_id"]
         and all(tgt[k]==trial[k] for k in ("trial_box_sha256","model_contract_sha256","point_cache_sha256",
                          "R_binary_sha256","center_binary_sha256","archived_center_shift_exact","exact_drive","energy_scale_binary_exact"))
         and q(leaf["radius_exact"])==q(trial["trial_radius_exact"])==R,
         "same native independently observed root, proof center, material and R")
    box=doc(trial["trial_box"],trial["trial_box_sha256"],links)
    must(box["record_id"]==tgt["record_id"] and box["point_cache_sha256"]==tgt["point_cache_sha256"]
         and box["model_observation_contract_sha256"]==tgt["model_contract_sha256"] and q(box["corrector_radius_exact"])==R
         and bdoc["record_id"]==prior["record_id"]==bd["record_id"]==ctr["physical_source_record_id"]
         and saved["step"]==bdoc["step"] and saved["role"]==bdoc["role"]=="small_record"
         and secondary_a["input_and_output_sha256"][file]==sp
         and saved["endpoint_id"]==endpoint and saved["checkpoint"]==checkpoint
         and saved["field_sha256"]==prior["field_sha256"]==bd["field_sha256"]==ctr["archived_field_sha256"]
                                  ==tgt["original_observer_field_sha256"]
         and saved["field"]==bd["field"]==ctr["archived_field"]
         and saved["verification"]==ctr["archived_verification"] and saved["verification_sha256"]==ctr["archived_verification_sha256"]
         and all(bdoc[k]==saved[k] for k in ("field","field_sha256","verification","verification_sha256","record_id"))
         and q(tgt["observer_goal_exact"])==q(bd["e_goal_dimensionless_exact"])==q(bdoc["e_goal_dimensionless_exact"])
               ==q(ctr["observer_goal_exact"])
         and q(tgt["original_addition_gap_exact"])==q(bd["binary_gap_exact"])
               ==q(bdoc["center_and_model"]["node_addition_gap_dimensionless_exact"])==q(ctr["archived_addition_gap_exact"])
         and q(tgt["archived_center_shift_exact"])==q(ctr["archived_shift_exact"])==q(prior["root"]["sigma"])
         and ctr["archived_observer_variant"]==VAR[0] and q(ctr["archived_observer_allowance_um_exact"])==allow,
         "same source-linked original archive, verification, original observer budget and recenter provenance")
    requested=ctr["exact_requested_drive"]
    must(q(tgt["exact_drive"])==q(requested["activation_strength_exact"])==q(bdoc["center_and_model"]["exact_target_strength"])
         and q(requested["common_gain_exact"])==q(bdoc["center_and_model"]["gain_exact"])
         and q(requested["optical_dose_exact"])==q(box["dose_exact"])==q(saved["step"]["dose_exact"])
         and q(requested["activation_strength_exact"])==q(requested["common_gain_exact"])*q(requested["optical_dose_exact"]),
         "exact gain-dose product, original fixed input")
    for k in ("field","verification"):
        path=saved[k]
        expected=saved[k+"_sha256"]
        must(primary["all_hashes"][path]==old_a["current_input_sha256"][path]
            ==secondary_a["input_and_output_sha256"][path]==expected,"all three independent recorded archived input digests agree")
    field=doc(saved["field"],saved["field_sha256"],links)
    verify=doc(saved["verification"],saved["verification_sha256"],links)
    e,sigma,xi,eps=independent_root(leaf)
    roi,graph=archived_geometry(field,box,bdoc,prior,saved,eps,sigma,xi)
    ea,eb=(q(graph[k]["graph_error_um_exact"]) for k in ("apical","basal"))
    same={"case_id":name,"root_proof_sha256":spec["native_sha256"],"record_id":bdoc["record_id"],
          "original_field_sha256":saved["field_sha256"],"posterior_exact":str(e),"archived_sigma_exact":str(sigma),
          "original_xi_exact":str(xi),"physical_nodal_error_um_exact":str(eps),"ROI_um_exact":list(map(str,roi)),
          "graph_premises":graph}
    must(all(reported[k]==v for k,v in same.items()) and len(reported["comparisons"])==16
         and verify["variant"]==VAR[0] and set(saved["variant_observables_um_exact"])==set(VAR),
         "independent root, source archive and geometry exactly equal")
    comps=[]
    for var in VAR:
        stored=saved["variant_observables_um_exact"][var]
        must(set(stored)==set(OBS),"same archived four observable labels")
        for obs,error in zip(OBS,(ea+eb,ea+eb,ea,eb),strict=True):
            actual=pick(reported["comparisons"],variant=var,observable=obs)
            sample=q(saved["step"]["observed_exact"][obs])
            if var==VAR[0]:
                prev=pick(prior["comparisons"],observable=obs)
                must(pair(stored[obs])==pair(verify["observer_intervals_um_exact"][obs])==pair(prev["stored_interval_um_exact"])
                     and sample==q(prev["observed_um_exact"]),"same archived primary baseline interval")
            roots,res,standing=interval_check(stored[obs],sample,error,allow)
            cell={"variant":var,"observable":obs,"stored_interval_um_exact":list(map(str,pair(stored[obs]))),
                 "observed_um_exact":str(sample),"extra_observer_error_um_exact":str(error),
                 "exact_spatial_root_observer_interval_um_exact":list(map(str,roots)),
                 "absolute_residual_um_exact":list(map(str,res)),"standing":standing}
            must(cell==actual,"independent exact observer or residual result mismatches producer")
            comps.append(cell)
    return {**same,"comparisons":comps}

def analytical_fixtures():
    must(interval_check(("0","1"),Q(2),Q(1),Q(1))[2]=="UNRESOLVED","interval crossing zero has exact zero lower")
    must(interval_check(("2","2"),Q(0),Q(1),Q(1))[2]=="UNRESOLVED","lower cutoff equality is not OUT")
    must(interval_check(("1","1"),Q(0),Q(0),Q(1))[2]=="IN","inclusive upper cutoff")
    must(interval_check(("-3","-2"),Q(0),Q(0),Q(1))[2]=="OUT","both signs exact")
    # P2 |N_0|+|N_1|+|N_2| has maximum 5/4 on the half-cell
    # from the quadratic 1+2t-4t^2 at t=1/4, with symmetry across t=1/2.
    t=Q(1,4)
    basis=((1-t)*(1-2*t),4*t*(1-t),t*(2*t-1))
    d0=(-3,4,-1)
    d1=(1,-4,3)
    must(sum(map(abs,basis))==Q(5,4) and sum(map(abs,d0))==sum(map(abs,d1))==8,
         "exact Lagrange trace bound constants")
    return {"closed_allowance_boundaries_and_both_signed_residuals":True,"exact_basis_value_and_derivative_constants":True}


def prepare():
    for path,h in PINS.items():
        must(sha(path)==h,"same immutable source bytes")
    p={"status":"FROZEN_FIRST_WAVE_OBSERVER_INDEPENDENT_RATIONAL_READBACK",
       "source_sha256":sha(SRC),"upstream_sha256":{rel(p):h for p,h in PINS.items()},
       "case_order":list(CASES),"variant_order":list(VAR),"length_um_exact":str(L),"radius_exact":str(R),
       "checks":analytical_fixtures(),"observer_producer_os_observation_required":True,
       "no_new_FE_AD_spatial_observer_kernel_jobs":True}
    ph=save(HERE/"protocol.json",p)
    print({"source_sha256":sha(SRC),"protocol_sha256":ph},flush=True)


def producer_os(h,links):
    osr=doc(OB/"actual/os_exit_observation.json",h,links)
    actual=doc(OB/"actual/results.json",links=links)
    term=doc(OB/"actual/terminal.json",links=links)
    started=doc(OB/"actual/started.json",links=links)
    must(osr["status"]==ROOT_OS_STATUS and osr["direct_root_tool_observation"] is True
         and osr["unified_exec_chunk_id"]=="b1f783" and osr["started_sha256"]==PINS[OB/"actual/started.json"]
         and osr["independent_actual_scalar_readback_required"] is True
         and all(osr[k]==actual[k] for k in ("source_sha256","protocol_sha256","preflight_sha256","source_review_sha256","native_os_sha256"))
         and type(osr["pid"]) is int and type(osr["observed_exit_code"]) is int and osr["observed_exit_code"]==0
         and osr["results_sha256"]==term["results_sha256"]==PINS[OB/"actual/results.json"]
         and osr["terminal_sha256"]==PINS[OB/"actual/terminal.json"]
         and term["status"]=="FINISHED" and type(term["exit_code"]) is int and term["exit_code"]==0
         and started["status"]=="RUNNING"
         and actual["status"]=="FIRST_WAVE_EXACT_SPATIAL_OBSERVERS_AWAIT_INDEPENDENT_READBACK"
         and actual["full_calibration_or_source_fit_claimed"] is False
         and all(type(x["pid"]) is int and x["pid"]==osr["pid"] and
           all(x[k]==actual[k] for k in ("source_sha256","protocol_sha256","preflight_sha256","source_review_sha256","native_os_sha256"))
                 for x in (term,started,actual)),"same genuine observed scalar application phase")
    return actual,osr


def run(args):
    p=doc(HERE/"protocol.json",args.protocol_sha256)
    must(p["source_sha256"]==sha(SRC) and p["upstream_sha256"]=={rel(k):h for k,h in PINS.items()}
         and p["checks"]==analytical_fixtures(),"same exact frozen independent source and protocol")
    links={}
    for path,h in PINS.items():
        must(sha(path)==h,"changed static original bytes")
        links[rel(path)]=h
    a,osr=producer_os(args.producer_os_sha256,links)
    phase,native_term=native_rows(links)
    sources=original_sources(links)
    allow=sources[-1]
    must(q(a["allowance_um_exact"])==allow and tuple(x["case_id"] for x in a["cases"])==CASES
         and type(a["new_force_AD_or_integral_jobs"]) is int and a["new_force_AD_or_integral_jobs"]==0,
         "same delivered firstwave fixed observed allowance, no mechanics")
    out=HERE/"actual"
    out.mkdir(exist_ok=False)
    mark={"source_sha256":sha(SRC),"protocol_sha256":args.protocol_sha256,"pid":os.getpid(),
          "producer_observed_exit_sha256":args.producer_os_sha256}
    save(out/"started.json",{"status":"RUNNING",**mark})
    terminal={"status":"FAILED","exit_code":1,**mark}
    try:
        cases=[read_case(s,native_term,sources,r,links) for s,r in zip(phase,a["cases"],strict=True)]
        counts={v:dict(Counter(x["standing"] for d in cases for x in d["comparisons"] if x["variant"]==v)) for v in VAR}
        must(counts==a["counts"]=={VAR[0]:{"IN":16},VAR[1]:{"IN":15,"OUT":1},VAR[2]:{"IN":16},VAR[3]:{"IN":16}}
             and a["total_primary_comparisons"]==16 and a["total_alternative_comparisons"]==48,
             "all exact64 outcome rows and four independently rebuilt root observer envelopes")
        outs=[{"case_id":d["case_id"],"variant":c["variant"],"observable":c["observable"],
               "strict_OUT_margin_um_exact":str(pair(c["absolute_residual_um_exact"])[0]-allow),
               "absolute_residual_um_exact":c["absolute_residual_um_exact"]}
              for d in cases for c in d["comparisons"] if c["standing"]=="OUT"]
        must(len(outs)==1 and q(outs[0]["strict_OUT_margin_um_exact"])>0,"the single exact signed OUT has positive separation")
        h=save(out/"receipt.json",{"status":"PASS_INDEPENDENT_ACTUAL_FIRST_WAVE_SPATIAL_OBSERVER_READBACK",
               **mark,"source_and_data_sha256":links,"scope":"Four exact spatial small fields only, fixed archived observers and allowance.",
               "producer_observed_pid":osr["pid"],"allowance_um_exact":str(allow),"cases":cases,"counts":counts,"strict_out":outs,
               "all_64_exact_comparisons_recreated":True,"complete_136_or_source_fit_claimed":False,
               "new_force_FE_AD_spatial_integral_or_solve_jobs":0,"producer_apply_imported":False})
        terminal.update(status="FINISHED",exit_code=0,receipt_sha256=h)
        print({"receipt_sha256":h,"counts":counts,"strict_out":outs,**mark},flush=True)
    except BaseException as err:
        save(out/"failure.json",{"error":str(err),"traceback":traceback.format_exc(),**mark})
        raise
    finally:
        save(out/"terminal.json",terminal)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    cmd=parser.add_subparsers(dest="phase",required=True)
    cmd.add_parser("prepare")
    actual=cmd.add_parser("run")
    actual.add_argument("--protocol-sha256",required=True)
    actual.add_argument("--producer-os-sha256",required=True)
    args=parser.parse_args()
    prepare() if args.phase=="prepare" else run(args)
