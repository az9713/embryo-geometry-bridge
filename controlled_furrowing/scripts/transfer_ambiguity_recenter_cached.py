"""Prepare two explicit cached-quadrature recenter trials and gate all actual work."""
import argparse
import base64
from fractions import Fraction as Q
from hashlib import sha256
import io
import json
from math import isfinite, nextafter, inf
import os
from pathlib import Path
import struct
import subprocess
import sys
import traceback
from unittest.mock import patch

import numpy as np
import review_algebraic_root as ad

SOURCE=Path(__file__).resolve()
ROOT=SOURCE.parents[1]
DESIGN=ROOT/"runs/transfer_ambiguity_nonlinear_probe_v1/certification_design"
OLD=DESIGN/"full_outcome_coverage_v1/remaining_cached_30_v1"
NATIVE=OLD/"independent_actual_preparation"
OUT=OLD/"recenter_v1"
PREP=OUT/"preparation"
ACTUAL=OUT/"actual"
CASES={
    "max_negative_contrast_change_small_16":("wave_03","optoGEF2_small_col_K"),
    "max_contrast_change_small_13":("wave_04","optoGEF2_small_col_E"),
}
OLD_PINS={
    OLD/"negative_K_recenter_v1/receipt.json":"d9621eef3cc7fec129b2f597dd4ab584393a6ed0b610051498734fdd2cad1d04",
    NATIVE/"replay.py":"19bc3760b36f10d9fd4dddea576d4af232af6b0667fd34e1f1580bbffac3529c",
    NATIVE/"protocol.json":"ab70239db02c62a7de672ef7b1f258779d7ecc06d0d431e4c4a4648a134a1349",
    NATIVE/"preflight.json":"8ae4f1b36cf4ad29ba0e200b0181876127f001fd11aeafd65ade81f05007d9c8",
    ROOT/"scripts/review_algebraic_root.py":"6af2311d0667a28af263b6fb6237063481aa8e65209ae9722c2b0246a8cff27a",
    DESIGN/"budgets_v1/actual/results.json":"689948766db950ee4ddf399e6731aea08890cc7a1f84c5615d19e7ded68b0872",
    DESIGN/"budgets_v1/independent_actual/receipt.json":"978a2a4095bb8735103a5ee2c11ad266a6b0e052e22af2b7a36ef29f1f38265c",
    NATIVE/"actual_replay/wave_03/os_exit_observation.json":"5b7636f8ec7bedd6fd7a3502c17f9e32cc723a84dbf8e273650d7049d1df35bf",
    NATIVE/"actual_replay/wave_04/os_exit_observation.json":"4f7e0968f8462570e505b98135192d4dbb2bb52c7208538fec81c32c237439e8",
}
RADIUS=Q(1,2**36)
LENGTH=Q(67,2)
PRECISION=256
NFREE=611
COUNTS={"bulk":9216,"surface":288}
MARGINS={"bulk_det_min","bulk_hoop_min","bulk_J_min","normalized_weight_min",
         "surface_radius_min","surface_tangent_squared_min","surface_J_min"}
SCOPE=("One binary recenter per approved unresolved archived state with the same exact-real frozen Q8/S12 "
       "potential, exact requested strength and unchanged R. The old saved field remains the scoring and "
       "observer authority. No spatial integral, whole-body embedding, minimum, source or continuum proof.")
REVIEW="PASS_INDEPENDENT_COLD_CACHED_RECENTER_PRELAUNCH"


def require(ok,message):
    if not ok:
        raise ValueError(message)


def sha(data):
    return sha256(data).hexdigest()


def filename(path):
    return Path(path).resolve().relative_to(ROOT).as_posix()


def path_for(name):
    require(type(name) is str,"canonical source path required")
    p=(ROOT/name).resolve()
    require(p.is_relative_to(ROOT) and filename(p)==name,"outside or noncanonical source authority")
    return p


def file_sha(path):
    return sha(Path(path).read_bytes())


def link(links,name,digest):
    path_for(name)
    require(type(digest) is str and len(digest)==64 and all(x in "0123456789abcdef" for x in digest)
            and links.setdefault(name,digest)==digest,"missing or conflicting exact artifact fingerprint")


def read(path,digest):
    data=Path(path).read_bytes()
    require(sha(data)==digest,f"changed frozen file: {path}")
    return json.loads(data)


def load(links,name,digest):
    link(links,name,digest)
    return read(path_for(name),digest)


def save(path,value):
    with path.open("x",encoding="utf8",newline="\n") as f:
        json.dump(value,f,indent=2,allow_nan=False)
        f.write("\n")
    return file_sha(path)


def yes(d):
    return type(d) is dict and bool(d) and all(x is True for x in d.values())


def zero(code):
    return type(code) is int and code==0


def positive_pid(x):
    return type(x) is int and x>0


def select(rows,**fields):
    matching=[x for x in rows if all(x[k]==v for k,v in fields.items())]
    require(len(matching)==1,"expected exact single retained case")
    return matching[0]


def pair(d):
    l,u=Q(d["lower_exact"]),Q(d["upper_exact"])
    require(l<=u,"ordered native signed force interval required")
    return l,u


def model_identity(case):
    d=case["drive"]
    dose,gain,drive,rounded=(Q(d[k]) for k in ("optical_dose_exact","common_gain_exact",
                                                "activation_strength_exact","activation_strength_binary_exact"))
    cell=case["closed_parameter_cell_exact"]
    material=case["material"]
    names=tuple(case["parameter_order"])
    require(set(names)==set(cell) and len(names)==5 and all(len(v)==2 and Q(v[0])==Q(v[1]) for v in cell.values())
            and dose>0 and gain>0 and drive==dose*gain>0 and drive!=rounded==Q.from_float(float(drive))
            and gain==Q(cell["common_gain"][0])==Q(material["gain_exact"])
            and dose==Q(case["step"]["dose_exact"]) and drive==Q(case["step"]["requested_strength_exact"]),
            "same fixed material cell and unrounded exact drive required")
    for k in set(names)-{"common_gain"}:
        require(Q(cell[k][0])==Q(float(material["parameters"][k]))==Q(material["parameter_storage_exact"][k]),
                "fixed material differs from its exact binary cell")
    return dose,gain,drive


def old_authority():
    links={}
    for path,h in OLD_PINS.items():
        name=filename(path)
        link(links,name,h)
        require(file_sha(path)==h,"frozen prerequisite bytes changed")
    design=read(OLD/"negative_K_recenter_v1/receipt.json",OLD_PINS[OLD/"negative_K_recenter_v1/receipt.json"])
    require(design["status"]=="READ_ONLY_NEGATIVE_K_CACHED_RECENTER_DESIGN_COMPLETE"
            and design["actual_recentered_seed_constructed"] is False and design["new_model_cache_AD_or_solve_jobs"]==0
            and len(design["source_and_input_sha256"])==18,"read-only 18-link case design required")
    # Cache/field data are still unopened. Copy published digests only; other
    # design links are verified as metadata/provenance in this cold phase.
    heavy={"cache_packet.json","field.json","verification.json"}
    for name,h in design["source_and_input_sha256"].items():
        link(links,name,h)
        if Path(name).name not in heavy:
            require(file_sha(path_for(name))==h,"design dependency changed")
    bp=DESIGN/"budgets_v1/actual/results.json"
    bi=DESIGN/"budgets_v1/independent_actual/receipt.json"
    budgets=read(bp,OLD_PINS[bp])
    review=read(bi,OLD_PINS[bi])
    require(review["status"]=="PASS_INDEPENDENT_ACTUAL_CONDITIONAL_OBSERVER_BUDGETS"
            and yes(review["all_literal_checks"]) and review["root_error_bounds_proved"]==0
            and positive_pid(review["observed_budget_terminal"]["pid"]) and zero(review["observed_budget_terminal"]["exit_code"]),
            "same independently observed frozen science budgets required")
    npf=read(NATIVE/"preflight.json",OLD_PINS[NATIVE/"preflight.json"])
    require(npf["status"]=="PASS" and npf["protocol_sha256"]==OLD_PINS[NATIVE/"protocol.json"]
            and npf["source_sha256"]==OLD_PINS[NATIVE/"replay.py"],"native phase changed")
    native_protocol=read(NATIVE/"protocol.json",OLD_PINS[NATIVE/"protocol.json"])
    for name,h in native_protocol["immutable_source_and_input_hashes"].items():
        if name.endswith(".py"):
            link(links,name,h)
            require(file_sha(path_for(name))==h,"native kernel or dependency changed")
    ids=[]
    for ident,(wave,record) in CASES.items():
        rw=NATIVE/"actual_replay"/wave
        a=read(rw/"os_exit_observation.json",OLD_PINS[rw/"os_exit_observation.json"])
        require(a["status"]=="PASS_OBSERVED_INDEPENDENT_CACHED_WAVE_TRUE_OS_EXIT" and a["wave_id"]==wave
                and positive_pid(a["pid"]) and positive_pid(a["unified_exec_session"]) and zero(a["observed_exit_code"]),
                "native observed OS zero exit required")
        nr=load(links,filename(rw/"receipt.json"),a["receipt_sha256"])
        nt=load(links,filename(rw/"terminal.json"),a["terminal_sha256"])
        require(nr["status"]=="INDEPENDENT_CACHED_WAVE_REPLAY_WITH_UNRESOLVED_CASES"
                and nr["wave_id"]==nt["wave_id"]==wave and nr["pid"]==nt["pid"]==a["pid"]
                and nt["status"]=="FINISHED" and zero(nt["exit_code"]) and nt["receipt_sha256"]==a["receipt_sha256"]
                and nr["root_terminal_admission_sha256"]==nt["root_terminal_admission_sha256"],
                "complete independently replayed original wave required")
        for k,n in (("source_sha256","replay.py"),("protocol_sha256","protocol.json"),("preflight_sha256","preflight.json")):
            require(nr[k]==nt[k]==OLD_PINS[NATIVE/n],"native old phase hash changed")
        row=select(nr["cases"],case_id=ident,record_id=record)
        require(row["status"]=="UNRESOLVED_INDEPENDENT_CACHED_ROOT" and row["native_file"]==filename(rw/(ident+".json")),
                "only designated previously unresolved native case allowed")
        proof=load(links,row["native_file"],row["native_sha256"])
        require(all(row[k]==proof[k] for k in ("status","case_id","record_id","field_sha256","cache_packet_sha256"))
                and proof["rounded_strength_used_in_proof"] is False and proof["site_counts"]==COUNTS
                and proof["chosen_radius_exact"] is None and proof["posterior_error_upper_exact"] is None,
                "old unresolved native anchor only")
        aroot=OLD/"actual"/wave
        ao=load(links,filename(aroot/"root_terminal_admission.json"),nr["root_terminal_admission_sha256"])
        require(ao["status"]=="PASS_ROOT_OBSERVED_REMAINING30_CACHED_WAVE_TERMINALS"
                and ao["wave_id"]==wave and zero(ao["observed_exit_code"])
                and positive_pid(ao["supervisor_pid"]) and positive_pid(ao["unified_exec_session"]),
                "original producer true supervisor exit required")
        trpath=filename(aroot/"terminal.json")
        end=load(links,trpath,ao["sha256"][trpath])
        require(end["status"]=="FINISHED" and end["pid"]==ao["supervisor_pid"] and end["wave_id"]==wave,
                "same original producer terminal")
        w=select(end["workers"],case_id=ident)
        child=aroot/ident
        os_path=filename(child/"observed_os_exit.json")
        old_os=load(links,os_path,ao["sha256"][os_path])
        child_path=filename(child/"terminal.json")
        child_end=load(links,child_path,ao["sha256"][child_path])
        require(w["observed_exit_code"]==0 and zero(old_os["observed_exit_code"]) and positive_pid(w["pid"])
                and old_os["pid"]==child_end["pid"]==w["pid"] and w["supervisor_pid"]==end["pid"]
                and w["terminal_sha256"]==ao["sha256"][child_path] and w["os_exit_receipt_sha256"]==ao["sha256"][os_path]
                and child_end["status"]=="FINISHED","original same-case child completion evidence")
        old_result=filename(child/"results.json")
        old=load(links,old_result,ao["sha256"][old_result])
        require(child_end["results_sha256"]==ao["sha256"][old_result]
                and old["status"]=="UNRESOLVED_CACHED_CONTRACTION_OR_BUDGET" and old["accepted_radius"] is None
                and old["cache_packet_sha256"]==proof["cache_packet_sha256"]
                and old["embedded_binary_cache_sha256"]==proof["embedded_cache_sha256"]
                and old["cache_packet"]==filename(child/"cache_packet.json"),"same unresolved producer cache evidence")
        c=old["case"]
        require(c["case_id"]==ident and c["record_id"]==record and c["field"]==proof["field"]
                and c["field_sha256"]==proof["field_sha256"],"wrong physical or prediction identity")
        _,_,drive=model_identity(c)
        goal=Q(c["e_goal_dimensionless_exact"])
        xi=Q(c["xi_dimensionless_exact"])
        b=select(budgets["field_budgets"],field=c["field"],field_sha256=c["field_sha256"])
        v=select(review["fields"],field=c["field"],field_sha256=c["field_sha256"],record_id=record)
        require(drive==Q(proof["exact_requested_strength"]) and 0<RADIUS<goal==Q(proof["observer_radius_budget_exact"])
                ==Q(v["e_goal_dimensionless_exact"])<=Q(b["e_goal_dimensionless_exact"])
                and 0<=xi==Q(v["binary_gap_exact"])==Q(b["center_and_model"]["node_addition_gap_dimensionless_exact"])
                and c["observer_variant"]=="current_disk__spatial_max"
                and Q(c["frozen_conditional_allowance_um_exact"])==Q(budgets["conditional_source_allowance_um_exact"]),
                "same original frozen observer allowance, goal and addition-rounding gap")
        attempt=select(proof["own_radius_attempts"],radius_exact=str(RADIUS))
        require(attempt["status"]=="UNRESOLVED_INDEPENDENT_CACHED_CONTRACTION"
                and attempt["whole_cube_AD_recomputed_for_this_radius"] is True
                and attempt["posterior_error_upper_exact"] is None,
                "only original independently enclosed signed center anchor authorized")
        require(len(attempt["signed_center_force_intervals"])==NFREE
                and max(max(map(abs,pair(x))) for x in attempt["signed_center_force_intervals"])
                   <=Q(attempt["eta_upper_exact"]),"complete original signed RF rows required")
        for name,h in ((c["field"],c["field_sha256"]),(c["verification"],c["verification_sha256"]),
                       (old["cache_packet"],old["cache_packet_sha256"])):
            require(nr["producer_and_payload_sha256"].get(name)==h,"native replay lacks exact original operand")
            link(links,name,h)
        for r in old["radius_attempts"]:
            name=filename(child/r["path"])
            require(nr["producer_and_payload_sha256"].get(name)==r["sha256"],"original radius attempt changed")
            link(links,name,r["sha256"])
        if wave=="wave_03":
            require(design["case_id"]==ident and design["exact_prescribed_drive"]==c["drive"]
                    and design["R_binary_sha256"]==proof["per_tensor_sha256"]["preconditioner"]
                    and design["original_per_tensor_sha256_from_independent_replay"]==proof["per_tensor_sha256"],
                    "negative K explicit design identity changed")
        ids.append({"case":c,"native_proof":row["native_file"],"native_proof_sha256":row["native_sha256"],
                    "native_os_sha256":OLD_PINS[rw/"os_exit_observation.json"],
                    "native_attempt":str(RADIUS),"cache_packet":old["cache_packet"],
                    "cache_packet_sha256":old["cache_packet_sha256"],"old_result":old_result,
                    "old_result_sha256":ao["sha256"][old_result],"original_failed_radii":old["radius_attempts"],
                    "cache_payload_sha256":old["embedded_binary_cache_sha256"],
                    "per_tensor_sha256":proof["per_tensor_sha256"]})
    return dict(sorted(links.items())),ids


def array_sha(a):
    return sha(a.tobytes(order="C"))


def decode(packet):
    d=packet["frozen_binary_cache"]
    data=base64.b64decode(d["data"],validate=True)
    require(len(data)==d["bytes"] and sha(data)==d["sha256"],"exact embedded cache bytes changed")
    with np.load(io.BytesIO(data),allow_pickle=False) as f:
        a={k:f[k].copy() for k in f.files}
    require({k:list(x.shape) for k,x in a.items()}==d["array_shapes"] and all(
        v.dtype.kind in ("i","u","f","b") and np.isfinite(v).all() for v in a.values()),"checked finite cache tensors")
    return a


def verify_tensors(before,after,wanted,seed):
    require(before.keys()==after.keys()==wanted.keys(),"changed tensor list")
    for k,a in before.items():
        b=after[k]
        require(a.shape==b.shape and a.dtype==b.dtype and array_sha(a)==wanted[k],"wrong original tensor/type/shape")
        if k!="u":
            require(array_sha(a)==array_sha(b),"changed frozen map, coefficient, weight or preconditioner")
    require(array_sha(after["u"])==array_sha(seed) and array_sha(before["u"])!=array_sha(seed),
            "new exact center must be the single proposed binary vector")


def round_binary(x):
    y=float(x)
    require(isfinite(y),"finite dyadic proposal required")
    prev,nxt=nextafter(y,-inf),nextafter(y,inf)
    require(isfinite(prev) and isfinite(nxt),"finite rounding-neighbor proof required")
    ym=Q.from_float(y)
    lo=(Q.from_float(prev)+ym)/2
    hi=(ym+Q.from_float(nxt))/2
    even=struct.unpack("<Q",struct.pack("<d",y))[0]%2==0
    require(lo<=x<=hi and (lo<x<hi or even),"nearest ties-to-even witness failed")
    return y


def new_seed(u,rows):
    require(u.ndim==1 and u.dtype.str=="<f8" and len(u)==len(rows),"complete exact old vector and signed rows")
    exact=[Q.from_float(float(x)) for x in u]
    mid=[(l+h)/2 for l,h in map(pair,rows)]
    target=[x-d for x,d in zip(exact,mid,strict=True)]
    seed=np.array(list(map(round_binary,target)),dtype="<f8")
    witness=[]
    for x,d,t,y in zip(exact,mid,target,seed,strict=True):
        q=Q.from_float(float(y))
        witness.append([str(x),str(d),str(t),str(q),str(q-x)])
    sigma=max(abs(Q(row[4])) for row in witness)
    return seed,witness,sigma


def classify(eta,z0,dz,sigma,goal,xi,r=RADIUS):
    z=z0+dz
    require(0<=min(eta,z0,dz,sigma,xi) and 0<r<goal,"positive fixed observer goal and radius")
    lhs=eta+z*r
    selfmap=z<1 and lhs<r
    post=eta/(1-z) if selfmap else None
    passes=selfmap and sigma+post<goal
    return {"eta_upper_exact":str(eta),"center_Z0_upper_exact":str(z0),"variation_upper_exact":str(dz),
            "ZQ_upper_exact":str(z),"strict_selfmap_lhs_upper_exact":str(lhs),"strict_contraction_selfmap":selfmap,
            "posterior_error_upper_exact":str(post) if selfmap else None,
            "shift_plus_posterior_exact":str(sigma+post) if selfmap else None,
            "strict_old_observer_goal_pass":passes,
            "physical_nodal_coordinate_error_um_exact":str(LENGTH*(sigma+post+xi)) if passes else None}


def prove_after_recenter(a,parameters,strength):
    """Fresh signed RF, point derivative and whole-cube variation of the same cache."""
    ad.ctx.prec=PRECISION
    ad.A.cache_clear()
    g,H,sites,count=ad.assemble(ad.terms_from_arrays(a,parameters,strength),NFREE)
    require(count==COUNTS and len(sites)==sum(COUNTS.values()),"fresh complete AD center coverage")
    rr=[Q.from_float(x) for row in a["preconditioner"].tolist() for x in row]
    R=ad.arb_mat(NFREE,NFREE,list(map(ad.A,rr)))
    AR=ad.arb_mat(NFREE,NFREE,[ad.A(abs(x)) for x in rr])
    I=ad.arb_mat(NFREE,NFREE)
    for i in range(NFREE):
        I[i,i]=1
    def bound(M):
        return max(Q(str(sum((M[i,j].abs_upper() for j in range(M.ncols())),ad.arb(0)).upper().fmpq()))
                   for i in range(M.nrows()))
    rf=R*g
    eta=bound(rf)
    z0=bound(I-R*H)
    spread,domain=ad.variation(sites,NFREE,ad.A(RADIUS))
    require(domain["counts"]==COUNTS and set(domain["strict_lower_margins"])==MARGINS,
            "same fresh complete analytic cube")
    mins={k:pair(v)[0] for k,v in domain["strict_lower_margins"].items()}
    require(min(mins.values())>0 and mins["bulk_J_min"]>Q(1,10**6),"fresh analytic positive-domain guard")
    return eta,z0,bound(AR*spread),domain,[ad.pack(rf[i,0]) for i in range(NFREE)]


def extract_actual(item):
    c=item["case"]
    proof=read(path_for(item["native_proof"]),item["native_proof_sha256"])
    old=read(path_for(item["cache_packet"]),item["cache_packet_sha256"])
    state=read(path_for(c["field"]),c["field_sha256"])
    read(path_for(c["verification"]),c["verification_sha256"])
    dose,gain,drive=model_identity(c)
    require(old["status"]=="FROZEN_FRESH_EXACT_BINARY_QUADRATURE_CACHE" and old["case_id"]==c["case_id"]
            and old["record_id"]==c["record_id"] and old["input_state"]==c["field"]
            and old["input_state_sha256"]==c["field_sha256"] and old["field_verification_sha256"]==c["verification_sha256"]
            and old["parameters"]==state["parameters"]==c["material"]["parameters"]
            and old["material_parameter_storage_exact"]==c["material"]["parameter_storage_exact"]
            and old["closed_parameter_cell_exact"]==c["closed_parameter_cell_exact"]
            and old["parameter_order"]==c["parameter_order"] and all(old[k]==v for k,v in c["drive"].items())
            and old["free_dof"]==NFREE and old["frozen_binary_cache"]["sha256"]==item["cache_payload_sha256"],
            "exact admitted source, physical cache and archived observer target required")
    arr=decode(old)
    require({k:array_sha(a) for k,a in arr.items()}==item["per_tensor_sha256"]
            and arr["u"].shape==arr["free"].shape==(NFREE,) and arr["preconditioner"].shape==(NFREE,NFREE)
            and arr["u"].dtype.str==arr["preconditioner"].dtype.str=="<f8"
            and list(map(int,arr["free"]))==c["free_map"] and len(set(c["free_map"]))==NFREE
            and np.array_equal(arr["u"],state["free_displacement"])
            and np.array_equal(arr["reference"],state["reference_nodes"])
            and np.array_equal(arr["triangles"],state["triangles"])
            and Q(float(arr["energy_scale"][0]))==Q(old["normalization_binary_exact"])>0,
            "current actual old center, full tensors and exact normalization required")
    rr=select(proof["own_radius_attempts"],radius_exact=str(RADIUS))
    return old,arr,drive,rr["signed_center_force_intervals"]


def case_trial(item,folder,bindings):
    c=item["case"]
    old,a,drive,rows=extract_actual(item)
    u,witness,sigma=new_seed(a["u"],rows)
    require(sigma>0,"single recenter must change the exact binary center")
    digest=array_sha(u)
    step=save(folder/"exact_binary_step.json",{"case_id":c["case_id"],"ordered_rows":{
            "columns":["old_u","native_signed_RF_midpoint","exact_unrounded_u1","rounded_binary_u1","exact_binary_shift"],
            "values":witness},"sigma_exact":str(sigma),"u0_sha256":item["per_tensor_sha256"]["u"],
            "u1_sha256":digest,"native_proof_sha256":item["native_proof_sha256"],"round_once_ties_to_even":True})
    changed={**a,"u":u}
    blob=io.BytesIO()
    np.savez_compressed(blob,**changed)
    data=blob.getvalue()
    old_payload=old["frozen_binary_cache"]
    checkpoint={**old,"status":"FROZEN_SAME_CACHED_OPERATOR_AT_NEW_PROOF_CENTER",
        "frozen_binary_cache":{**old_payload,"data":base64.b64encode(data).decode(),
                               "sha256":sha(data),"bytes":len(data)},
        "recenter":{"original_point_cache_sha256":item["cache_packet_sha256"],"original_payload_sha256":item["cache_payload_sha256"],
                    "original_field_retains_scoring_authority":True,"updated_tensor_keys":["u"],"u1_sha256":digest,
                    "old_field_is_not_new_proof_center":True,"exact_step_sha256":step,"bindings":bindings}}
    fresh=decode(checkpoint)
    verify_tensors(a,fresh,item["per_tensor_sha256"],u)
    require(all(v==checkpoint[k] for k,v in old.items() if k not in ("status","frozen_binary_cache")),
            "other original header or model identity changed")
    ch=save(folder/"recentered_cache_packet.json",checkpoint)
    goal,xi=Q(c["e_goal_dimensionless_exact"]),Q(c["xi_dimensionless_exact"])
    base={"case_id":c["case_id"],"record_id":c["record_id"],"scope":SCOPE,"bindings":bindings,
          "archived_field_sha256":c["field_sha256"],"same_observer_variant":c["observer_variant"],
          "same_frozen_allowance_um_exact":c["frozen_conditional_allowance_um_exact"],
          "old_failed_result_sha256":item["old_result_sha256"],"old_failed_radius_attempts":item["original_failed_radii"],
          "old_native_anchor_sha256":item["native_proof_sha256"],"original_cache_sha256":item["cache_packet_sha256"],
          "recentered_cache_packet_sha256":ch,"exact_step_sha256":step,"u1_binary_sha256":digest,
          "new_per_tensor_sha256":{k:array_sha(v) for k,v in fresh.items()},
          "original_R_binary_sha256":item["per_tensor_sha256"]["preconditioner"],
          "trial_radius_exact":str(RADIUS),"requested_strength_exact":str(drive),"actual_shift_sigma_exact":str(sigma),
          "original_e_goal_exact":str(goal),"original_binary_addition_gap_exact":str(xi),
          "original_observer_field_replaced":False,"original_scores_or_allowance_changed":False,
          "new_equilibrium_solves_refits_or_preconditioners":0,"exact_spatial_or_minimum_claimed":False,
          "independent_actual_recentered_AD_replay_required":True,"root_certified":False}
    try:
        eta,z0,dz,domain,force=prove_after_recenter(fresh,old["parameters"],drive)
    except ValueError as e:
        return {**base,"status":"UNRESOLVED_FRESH_CACHED_ANALYTIC_DOMAIN","error":str(e)}
    dec=classify(eta,z0,dz,sigma,goal,xi)
    return {**base,**dec,"status":"CANDIDATE_RECENTERED_CACHED_ROOT_AWAITING_INDEPENDENT_REVIEW" if dec[
        "strict_old_observer_goal_pass"] else "UNRESOLVED_RECENTERED_CACHED_ROOT_OR_SHIFT_BUDGET",
        "fresh_point_and_whole_cube_AD":True,"fresh_derivative_center_u1_sha256":digest,
        "signed_center_force_intervals":force,"analytic_domain":domain,"AD_precision_bits":PRECISION}


def protocol():
    links,cases=old_authority()
    return {"status":"FROZEN_TWO_OLD_CACHED_OPERATOR_RECENTER_TRIALS","source_sha256":file_sha(SOURCE),
            "scope":SCOPE,"source_and_input_sha256":links,"case_order":list(CASES),"cases":cases,
            "radius_exact":str(RADIUS),"precision_bits":PRECISION,"length_um_exact":str(LENGTH),
            "seed":"round64(exact archived u minus exact midpoint of admitted signed RFQ row), once",
            "observer":"sigma + eta/(1-z) < original e_goal; physical 33.5*(sigma + eta/(1-z) + old xi)",
            "actual_seed_cache_model_AD_proposal_or_solve_jobs":0}


def gate(path,digest):
    binding={"source_sha256":file_sha(SOURCE),"protocol_sha256":file_sha(PREP/"protocol.json"),
             "preflight_sha256":file_sha(PREP/"preflight.json"),"recenter_design_sha256":OLD_PINS[OLD/"negative_K_recenter_v1/receipt.json"]}
    rv=read(path,digest)
    pf=read(PREP/"preflight.json",binding["preflight_sha256"])
    p=read(PREP/"protocol.json",binding["protocol_sha256"])
    require(rv["status"]==REVIEW and yes(rv["checks"]) and all(rv[k]==v for k,v in binding.items())
            and pf["status"]=="PASS" and yes(pf["checks"]) and pf["actual_production_jobs"]==0
            and all(pf[k]==binding[k] for k in ("source_sha256","protocol_sha256"))
            and file_sha(PREP/"transfer_ambiguity_recenter_cached_frozen.py")==binding["source_sha256"]
            and p==protocol(),"independent current cold launch approval required")
    return p,{**binding,"prelaunch_review_sha256":digest}


def worker(args):
    p,links=gate(args.review,args.review_sha256)
    parent=read(ACTUAL/"started.json",file_sha(ACTUAL/"started.json"))
    require(parent["pid"]==os.getppid()==args.supervisor_pid and parent["bindings"]==links
            and parent["case_order"]==list(CASES),"authorized live serial supervisor")
    folder=ACTUAL/args.case_id
    mark={"case_id":args.case_id,"pid":os.getpid(),"supervisor_pid":args.supervisor_pid,"bindings":links}
    save(folder/"started.json",{"status":"RUNNING",**mark})
    end={"status":"FAILED","exit_code":1,**mark}
    try:
        item=select(p["cases"],native_proof=filename(NATIVE/"actual_replay"/CASES[args.case_id][0]/(args.case_id+".json")))
        result=case_trial(item,folder,links)
        rh=save(folder/"results.json",{**result,**mark})
        end.update(status="FINISHED",exit_code=0,results_sha256=rh)
    except BaseException as e:
        save(folder/"failure.json",{"error":str(e),"traceback":traceback.format_exc(),**mark})
        raise
    finally:
        save(folder/"terminal.json",end)


def child_admit(x,terminal,identity,h):
    require(zero(x) and zero(terminal["exit_code"]) and terminal["status"]=="FINISHED"
            and terminal["results_sha256"]==h and all(terminal[k]==v for k,v in identity.items()),
            "actual child OS and same current result required")


def run(args):
    _,links=gate(args.review,args.review_sha256)
    ACTUAL.mkdir(exist_ok=False)
    pid=os.getpid()
    save(ACTUAL/"started.json",{"status":"RUNNING","pid":pid,"bindings":links,"case_order":list(CASES)})
    end={"status":"FAILED","exit_code":1,"pid":pid,"bindings":links,"workers":[]}
    try:
        rows=[]
        for case in CASES:
            folder=ACTUAL/case
            folder.mkdir(exist_ok=False)
            cmd=[sys.executable,"-B",str(SOURCE),"worker","--case-id",case,
                 "--supervisor-pid",str(pid),"--review",str(args.review),"--review-sha256",args.review_sha256]
            with (folder/"stdout.log").open("x",encoding="utf8") as log:
                child=subprocess.Popen(cmd,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT)
                save(folder/"supervisor_child_started.json",{"pid":child.pid,"supervisor_pid":pid,"case_id":case,"command":cmd})
                code=child.wait()
            t=folder/"terminal.json"
            o={"case_id":case,"pid":child.pid,"supervisor_pid":pid,"observed_exit_code":code,
               "terminal_sha256":file_sha(t) if t.is_file() else None}
            oh=save(folder/"observed_os_exit.json",o)
            end["workers"].append({**o,"os_exit_sha256":oh})
            require(o["terminal_sha256"] is not None,"missing worker terminal")
            term=read(t,o["terminal_sha256"])
            rfile=folder/"results.json"
            require(rfile.is_file(),"missing worker output")
            h=file_sha(rfile)
            ident={"case_id":case,"pid":child.pid,"supervisor_pid":pid,"bindings":links}
            child_admit(code,term,ident,h)
            rec=read(rfile,h)
            require(all(rec[k]==v for k,v in ident.items()) and rec["root_certified"] is False,"same-case candidate only")
            rows.append({**o,"results_sha256":h,"status":rec["status"]})
        rh=save(ACTUAL/"results.json",{"status":"COMPLETE_TWO_RECENTER_ATTEMPTS_AWAITING_INDEPENDENT_REVIEW",
                 "scope":SCOPE,"bindings":links,"cases":rows,"root_certified":False})
        end.update(status="FINISHED",exit_code=0,results_sha256=rh)
        print(json.dumps({"results_sha256":rh,"case_statuses":[r["status"] for r in rows]}))
    except BaseException as e:
        end["error"]=str(e)
        raise
    finally:
        save(ACTUAL/"terminal.json",end)


def reject(f):
    try:
        f()
    except ValueError:
        return True
    return False


def fixtures():
    x=np.array([1.,0.5],dtype="<f8")
    rows=[{"lower_exact":"1/8","upper_exact":"1/8"},{"lower_exact":"-1/16","upper_exact":"-1/16"}]
    seed,w,sigma=new_seed(x,rows)
    require(w[0]==["1","1/8","7/8","7/8","-1/8"] and sigma==Q(1,8),"exact binary step oracle")
    require(round_binary(1+Q(1,2**53))==1,"nearest ties-to-even boundary")
    old={"u":x,"preconditioner":np.eye(2),"weights":np.ones(2),"free":np.array([0,1])}
    new={**old,"u":seed}
    h={k:array_sha(v) for k,v in old.items()}
    verify_tensors(old,new,h,seed)
    good=classify(Q(1,10),Q(1,8),Q(1,8),Q(1,5),Q(2),Q(1,100),Q(1))
    tight=classify(Q(1,10),Q(1,8),Q(1,8),Q(9,10),Q(31,30),Q(0),Q(1))
    tie=classify(Q(3,4),Q(1,8),Q(1,8),Q(0),Q(2),Q(0),Q(1))
    actual={"pid":2,"case_id":"a","supervisor_pid":1}
    end={**actual,"status":"FINISHED","exit_code":0,"results_sha256":"a"}
    child_admit(0,end,actual,"a")
    c={"drive":{"optical_dose_exact":"1/3","common_gain_exact":"3/2","activation_strength_exact":"1/2",
                   "activation_strength_binary_exact":"1/2"},"parameter_order":["a","b","t","s","common_gain"],
       "closed_parameter_cell_exact":{**{k:["1","1"] for k in "abts"},"common_gain":["3/2","3/2"]},
       "material":{"gain_exact":"3/2","parameters":dict.fromkeys("abts",1.),
                   "parameter_storage_exact":dict.fromkeys("abts","1")},
       "step":{"dose_exact":"1/3","requested_strength_exact":"1/2"}}
    # Use a nondyadic exact drive, as in both authorized real cases.
    c["drive"].update(optical_dose_exact="1/5",activation_strength_exact="3/10",
                      activation_strength_binary_exact=str(Q.from_float(.3)))
    c["step"].update(dose_exact="1/5",requested_strength_exact="3/10")
    require(model_identity(c)[2]==Q(3,10),"unchanged exact drive")
    return {"exact_once_rounded_binary_step_and_tie_fixture":w[0][3]=="7/8",
        "reject_changed_R_weight_free_map_or_u1_bytes":all(reject(lambda:verify_tensors(old,d,h,seed)) for d in (
            {**new,"preconditioner":2*np.eye(2)}, {**new,"weights":2*np.ones(2)},
            {**new,"free":np.array([1,0])}, {**new,"u":x})),
        "reject_changed_drive_or_mixed_case":reject(lambda:model_identity(
            {**c,"drive":{**c["drive"],"activation_strength_exact":"2/5"}}))
            and reject(lambda:select([{"case_id":"a"}],case_id="b")),
        "strict_selfmap_and_shift_inclusive_observer_bound":good["strict_old_observer_goal_pass"] is True
            and tight["strict_contraction_selfmap"] is True and tight["strict_old_observer_goal_pass"] is False
            and tie["strict_contraction_selfmap"] is False,
        "reject_nonliteral_or_mismatched_process_exit":reject(lambda:child_admit(False,end,actual,"a"))
            and reject(lambda:child_admit(0,{**end,"pid":3},actual,"a"))}


def prepare():
    PREP.mkdir(parents=True,exist_ok=False)
    with (PREP/"transfer_ambiguity_recenter_cached_frozen.py").open("xb") as f:
        f.write(SOURCE.read_bytes())
    try:
        blocked=AssertionError("actual seed, cache or AD entered cold phase")
        with patch.object(sys.modules[__name__],"extract_actual",side_effect=blocked), \
             patch.object(sys.modules[__name__],"prove_after_recenter",side_effect=blocked), \
             patch.object(sys.modules[__name__],"decode",side_effect=blocked), \
             patch.object(ad,"assemble",side_effect=blocked),patch.object(ad,"variation",side_effect=blocked):
            p=protocol()
            checks=fixtures()
        checks["two_explicit_old_cases_and_no_production_artifacts"]=len(p["cases"])==2 and not ACTUAL.exists()
        require(yes(checks),"cold recenter fixture failed")
        ph=save(PREP/"protocol.json",p)
        pf=save(PREP/"preflight.json",{"status":"PASS","source_sha256":file_sha(SOURCE),"protocol_sha256":ph,
             "checks":checks,"actual_production_jobs":0,"old_saved_fields_or_observer_scores_changed":False})
        print(json.dumps({"source_sha256":file_sha(SOURCE),"protocol_sha256":ph,"preflight_sha256":pf}))
    except BaseException as e:
        save(PREP/"failure.json",{"error":str(e),"traceback":traceback.format_exc()})
        raise


if __name__=="__main__":
    cli=argparse.ArgumentParser()
    sub=cli.add_subparsers(dest="stage",required=True)
    sub.add_parser("prepare")
    for n in ("run","worker"):
        s=sub.add_parser(n)
        s.add_argument("--review",required=True,type=lambda p:Path(p).resolve())
        s.add_argument("--review-sha256",required=True)
        if n=="worker":
            s.add_argument("--case-id",choices=CASES,required=True)
            s.add_argument("--supervisor-pid",type=int,required=True)
    a=cli.parse_args()
    {"prepare":lambda _:prepare(),"run":run,"worker":worker}[a.stage](a)

