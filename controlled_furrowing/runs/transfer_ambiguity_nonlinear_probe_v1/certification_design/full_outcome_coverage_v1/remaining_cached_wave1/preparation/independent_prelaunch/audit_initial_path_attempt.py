"""Independent cold check of the two fixed M-record cached root targets."""
import ast
import importlib.metadata
from copy import deepcopy
from fractions import Fraction as Q
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys
from unittest.mock import patch

HERE=Path(__file__).resolve().parent
PRE=HERE.parent
ROOT=PRE.parents[6]
PRIMARY=ROOT/"runs/transfer_ambiguity_nonlinear_probe_v1"
BASE=PRIMARY/"certification_design"
SOURCE=ROOT/"scripts/transfer_ambiguity_remaining_cached.py"
PROT=PRE/"protocol.json"
PREF=PRE/"preflight.json"
LEDGER=BASE/"full_outcome_coverage_v1/receipt.json"
MANIFEST=PRIMARY/"field_manifest.json"
BUDGET=BASE/"budgets_v1/independent_actual/receipt.json"
PINS={
    SOURCE:"a0ec457505c6a2e82aab91bd985d6f4383864124df4a1e295e3119e1f8694480",
    PROT:"c1543cbbd2520664313b56114b66a885faaf4d75869f56dbb1aa83a43a78fe3a",
    PREF:"cf328f574402555363ee4f74e2ca1afbb9cf507ee871d1a6501ca1bace9eb0dc",
    LEDGER:"235e8336c735b1bf1a75889018ae6a18586fd0905a65d38cb683c5ed236df19a",
    MANIFEST:"308ca4c05e476ba6a14cbc7b8d34018983b8f11263e2c224488672b35eb38ac0",
    BUDGET:"978a2a4095bb8735103a5ee2c11ad266a6b0e052e22af2b7a36ef29f1f38265c",
}
CASES=["max_contrast_change_small_09","max_negative_contrast_change_small_09"]
RID="optoGEF2_small_col_M"
DOSE=Q(78208913,100000000)
ORDER=("apical_area_modulus","basal_area_modulus","apical_surface_tension","basal_surface_tension","common_gain")
RADII=[Q(1,2**k) for k in (20,24,28,32,36)]
OBS=("maximum_shortening_um","region_averaged_shortening_um","apical_max_shortening_um","basal_max_shortening_um")


def sha(f):
    return sha256(f.read_bytes()).hexdigest()


def read(f):
    return json.loads(f.read_bytes())


def save(path,doc):
    with path.open("x",encoding="utf-8",newline="\n") as f:
        json.dump(doc,f,indent=2,allow_nan=False)
        f.write("\n")
    return sha(path)


def closure():
    seen=set()
    stack=[SOURCE]
    while stack:
        f=stack.pop()
        if f in seen:
            continue
        seen.add(f)
        for node in ast.walk(ast.parse(f.read_text(encoding="utf-8-sig"))):
            names=[n.name for n in node.names] if isinstance(node,ast.Import) else (
                [node.module] if isinstance(node,ast.ImportFrom) else [])
            for n in names:
                if n and (p:=ROOT/"scripts"/(n.split(".")[0]+".py")).exists():
                    stack.append(p)
    return {f.relative_to(ROOT).as_posix() for f in seen if f!=SOURCE}


def reject(fn):
    try:
        fn()
    except ValueError:
        return True
    return False


def review_case(c,ep,old_primary,b0,b1,manifest,ledger):
    assert c["case_id"]==f"{ep['endpoint_id']}_small_09"
    assert c["record_id"]==ledger["record_id"]==b0["record_id"]==b1["record_id"]==RID
    assert all(c[k]==manifest[k]==ledger[k]==b0[k] for k in ("field","field_sha256","verification","verification_sha256"))
    assert old_primary["all_hashes"][c["field"]]==c["field_sha256"]
    assert old_primary["all_hashes"][c["verification"]]==c["verification_sha256"]
    assert ledger["pilot_cached_and_geometry_coverage"] is None
    assert ledger["fresh_cache_and_preconditioner_needed"] is True
    assert c["step"]==manifest["step"]==b0["step"]==ep["small_records"][8]
    f=read(ROOT/c["field"])
    verification=read(ROOT/c["verification"])
    data=c["material"]
    m=f["parameters"]
    assert data==ep["materials"]["small"] and m==data["parameters"]
    assert {k:Q(v) for k,v in m.items()}=={k:Q(v) for k,v in data["parameter_storage_exact"].items()}
    assert f["converged"] is f["reached_target"] is True
    assert len(f["triangles"])==144 and len(f["reference_nodes"])==len(f["current_nodes"])==325
    ref=[tuple(map(Q,x)) for x in f["reference_nodes"]]
    rim=Q(m["outer_radius"])
    free=[2*i+j for i,(r,z) in enumerate(ref) for j in (0,1) if r!=rim and (r!=0 or j!=0)]
    assert c["free_map"]==free and len(free)==len(f["free_displacement"])==611
    assert c["active_elements"]==[i for i,v in enumerate(f["activation_elements"]) if v]
    assert len(c["active_elements"])==24
    offsets=dict(zip(free,map(Q,f["free_displacement"]),strict=True))
    gaps=[]
    for i,(a,b) in enumerate(zip(ref,f["current_nodes"],strict=True)):
        for j,(x,y) in enumerate(zip(a,b,strict=True)):
            k=2*i+j
            y=Q(y)
            true=x+offsets.get(k,0)
            assert Q(float(true))==y
            g=abs(y-true)
            assert k in offsets or g==0
            gaps.append(g)
    xi=max(gaps)
    assert xi==Q(c["xi_dimensionless_exact"])==Q(ledger["node_addition_gap_exact"])==Q(b1["binary_gap_exact"])
    gain=Q(data["gain_exact"])
    d=Q(c["step"]["dose_exact"])
    exact=gain*d
    rounded=Q(float(exact))
    assert d==DOSE and 0<exact!=rounded==Q(f["activation_strength"])
    assert {k:Q(v) for k,v in c["drive"].items()}=={
        "optical_dose_exact":d,"common_gain_exact":gain,"activation_strength_exact":exact,
        "activation_strength_binary_exact":rounded,"binary_input_minus_exact_strength":rounded-exact}
    vals=dict(zip(ORDER,[Q(m[k]) for k in ORDER[:-1]]+[gain],strict=True))
    assert c["parameter_order"]==list(ORDER)
    assert {k:list(map(Q,v)) for k,v in c["closed_parameter_cell_exact"].items()}=={k:[v,v] for k,v in vals.items()}
    assert c["observer_variant"]==verification["variant"]=="current_disk__spatial_max"
    delta=Q(c["frozen_conditional_allowance_um_exact"])
    ca,cb=(Q(b0["stored_face_derivatives"][name]["C_um_per_dimensionless_coordinate_exact"])
           for name in ("top","bottom"))
    A=(ca+cb,ca+cb,ca,cb)
    caps=[]
    for k,a in zip(OBS,A,strict=True):
        p=list(map(Q,verification["observer_intervals_um_exact"][k]))
        assert p[0]<=p[1]
        error=max(abs(y-Q(c["step"]["observed_exact"][k])) for y in p)
        margin=delta-error
        assert margin==Q(b0["stored_residuals"][k]["margin_um_exact"])==Q(ledger["stored_primary_observer_margins_um_exact"][k])>0
        caps.append(margin/(2*a))
    radial=min(Q(b0["stored_face_derivatives"][x]["radial_min_um_exact"]) for x in ("top","bottom"))/16/Q(67,2)
    goal=min(*caps,radial)-xi
    assert goal==Q(c["e_goal_dimensionless_exact"])==Q(b1["e_goal_dimensionless_exact"])==Q(ledger["observer_goal_exact"])
    assert 0<RADII[-1]<RADII[0]<goal/2
    return {"case_id":c["case_id"],"field_sha256":c["field_sha256"],"verification_sha256":c["verification_sha256"],
            "independent_radius_goal_exact":str(goal),"goal_display":float(goal),
            "binary_coordinate_gap_exact":str(xi),"exact_drive":str(exact),
            "rounded_reserve_not_used_as_force_input":True}


def main():
    assert all(sha(f)==h for f,h in PINS.items())
    p,pre,l,f,b1=map(read,(PROT,PREF,LEDGER,MANIFEST,BUDGET))
    assert sha(PRE/"transfer_ambiguity_remaining_cached_frozen.py")==PINS[SOURCE]
    assert p["source_sha256"]==pre["source_sha256"]==PINS[SOURCE] and pre["protocol_sha256"]==PINS[PROT]
    assert pre["status"]=="PASS" and len(pre["checks"])==9 and all(x is True for x in pre["checks"].values())
    assert all(type(pre[k]) is int and pre[k]==0 for k in (
        "actual_model_instances_constructed","actual_preconditioners_recovered","actual_binary_caches_frozen",
        "actual_cached_AD_sites_assembled","new_equilibrium_solves"))
    assert p["ordered_case_ids"]==CASES and p["old_pilot_case_reused"] is False
    assert p["site_counts"]=={"bulk":9216,"surface":288} and p["precision_bits"]==256 and p["torch_cpu_threads"]==1
    assert list(map(Q,p["finite_radius_ladder_exact"]))==RADII
    assert all(sha(ROOT/k)==v for k,v in (p["source_closure_sha256"]|p["input_sha256"]).items())
    assert closure()==set(p["source_closure_sha256"])
    assert all(importlib.metadata.version(k)==v for k,v in p["runtime_packages"].items())
    ps=read(PRIMARY/"preparation/protocol.json")
    nb=read(PRIMARY/"independent_actual_preparation/actual_replay/receipt.json")
    b0=read(BASE/"budgets_v1/actual/results.json")
    assert f["independent_replay_observed_os_exit_code"]==0 and nb["checkpoints"]==len(f["fields"])==44
    def key(r):
        return r["endpoint_id"],r["checkpoint"]
    tab=[{key(row):row for row in table} for table in (b0["field_budgets"],b1["fields"],f["fields"],l["outcome_states"])]
    result=[]
    for ep,c in zip(ps["endpoints"],p["cases"],strict=True):
        rows=[d[key(c)] for d in tab]
        result.append(review_case(c,ep,nb,*rows))
    assert [x["case_id"] for x in result]==CASES and len({x["field_sha256"] for x in result})==2
    assert not (PRE.parent/"actual").exists()

    sys.path.insert(0,str(ROOT/"scripts"))
    import transfer_ambiguity_remaining_cached as driver
    bound=dict(zip(("source_sha256","protocol_sha256","preflight_sha256","coverage_ledger_sha256",
        "whole_44_field_manifest_sha256","budget_independent_actual_review_sha256"),PINS.values(),strict=True))
    def forbidden(*a,**k):
        raise AssertionError("no model, cache or worker during readback")
    with (patch.object(driver,"actual_case",side_effect=forbidden),
          patch.object(driver,"worker",side_effect=forbidden),patch.object(driver,"launch",side_effect=forbidden),
          patch.object(subprocess,"Popen",side_effect=forbidden)):
        assert driver.prepare_protocol()==p
        assert reject(lambda:driver.unique([p["cases"][0]]*2,lambda x:True))
        assert reject(lambda:driver.case_inputs("max_contrast_change_small_17",driver.authority_inputs()))
        expected={"pid":45,"supervisor_pid":43,"case_id":CASES[0],"bindings":bound}
        terminal=expected|{"status":"FINISHED"}
        driver.admit_cached_worker(0,terminal,expected,expected)
        for code in (False,0.0,-1):
            assert reject(lambda:driver.admit_cached_worker(code,terminal,expected,expected))
        assert reject(lambda:driver.admit_cached_worker(0,terminal|{"pid":1},expected,expected))
        g=HERE/"synthetic_bad_gate.json"
        gh=save(g,{"status":"PASS_INDEPENDENT_REMAINING_M_CACHED_PRELAUNCH","checks":{"nonliteral":1},
                   "ordered_case_ids":CASES,**bound})
        assert reject(lambda:driver.prelaunch(g,gh))
        for i,v in enumerate(({"ordered_case_ids":list(reversed(CASES))},{"coverage_ledger_sha256":"3"*64})):
            w=HERE/f"synthetic_mixed_gate_{i}.json"
            h=save(w,{"status":"PASS_INDEPENDENT_REMAINING_M_CACHED_PRELAUNCH","checks":{"test":True},
                      "ordered_case_ids":CASES,**bound,**v})
            assert reject(lambda:driver.prelaunch(w,h))

    assert not (PRE.parent/"actual").exists()
    checks={"independently_reproduced_both_M_exact_coordinates_drives_and_observer_budgets":True,
            "same_full_local_Python_source_closure_and_exclusive_two_point_scope":True,
            "fresh_point_cache_R_and_high_precision_exact_force_driver_source_review":True,
            "actual_radii_strictly_below_both_observer_goals":True,
            "independent_nonliteral_mixed_duplicate_and_wrong_exit_fixtures":True,
            "serial_true_OS_child_observation_scope_and_no_actual_launch":True}
    support=HERE/"support.json"
    sh=save(support,{"status":"PASS","checks":checks,"scope":"Two fresh fixed-point cached-quadrature attempts only.",
            "new_model_cache_force_hessian_or_solve_evaluations":0,"cases":result,
            "reused_four_point_driver_scope":"Four-pilot source supplies explicit-argument numerical primitives only; no earlier pilot identity, cache or R loaded.",
            "dependencies":p["source_closure_sha256"]|p["input_sha256"],"review_source_sha256":sha(Path(__file__))})
    gate=HERE/"review.json"
    gh=save(gate,{"status":"PASS_INDEPENDENT_REMAINING_M_CACHED_PRELAUNCH","checks":checks,
         "ordered_case_ids":CASES,**bound,"support_file":support.relative_to(ROOT).as_posix(),"support_sha256":sh})
    print(gh,sh,flush=True)


if __name__=="__main__":
    main()
