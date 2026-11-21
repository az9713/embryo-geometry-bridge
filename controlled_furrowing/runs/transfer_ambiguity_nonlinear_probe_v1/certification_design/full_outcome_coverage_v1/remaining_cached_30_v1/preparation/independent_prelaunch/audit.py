"""Common cold source and case audit with one gate per immutable wave."""
import ast
from collections import Counter
from fractions import Fraction as Q
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys
from unittest.mock import patch

HERE=Path(__file__).resolve().parent
PRE=HERE.parent
ROOT=PRE.parents[5]
PRIMARY=ROOT/"runs/transfer_ambiguity_nonlinear_probe_v1"
DESIGN=PRIMARY/"certification_design"
COVER=DESIGN/"full_outcome_coverage_v1"
SRC=ROOT/"scripts/transfer_ambiguity_remaining30_cached.py"
PROTO=PRE/"protocol.json"
PRE_FLIGHT=PRE/"preflight.json"
INDEX=PRE/"wave_index.json"
LEDGER=COVER/"receipt.json"
M=COVER/"remaining_cached_wave1/VERIFIED-CACHED-COVERAGE.json"
FOUR=DESIGN/"cached_pilots_v1/independent_actual_preparation/actual_replay/receipt.json"
PINS={
    SRC:"23c6e4033bcb9f718a31b75a1c6238afe611fab0f99bf5cd5b16c17298abe702",
    PROTO:"e31870b0c6ac61e95a7b2657eab73d8231657af86c65f7a1d1172ff0b45a1e4c",
    PRE_FLIGHT:"1cad3ad4f7a3b8d486b16d511b0748d215830cab2fb6e74e16143bd3f066a11d",
    INDEX:"6aa242a734bfd6a2524aa8fae6f899285a0ba35305731abb16961f3b7c4d057c",
    LEDGER:"235e8336c735b1bf1a75889018ae6a18586fd0905a65d38cb683c5ed236df19a",
    M:"d5591fe4cfa66ef6e1d25eeaafde6df1000352c94225dd680584dab3018901b1",
    FOUR:"4e891afd79dd2f73508d1f878ddc1113ef107bdff8d9e0f673894f4d2f10ae17",
}
EP=("max_contrast_change","max_negative_contrast_change")
WAVES=tuple(f"wave_{i:02d}" for i in range(1,9))
RADII=[Q(1,2**k) for k in (20,24,28,32,34,36)]
OBS=("maximum_shortening_um","region_averaged_shortening_um","apical_max_shortening_um","basal_max_shortening_um")
MAT=("apical_area_modulus","basal_area_modulus","apical_surface_tension","basal_surface_tension")


def sha(path):
    return sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_bytes())


def save(p,x):
    with p.open("x",encoding="utf-8",newline="\n") as h:
        json.dump(x,h,indent=2,allow_nan=False)
        h.write("\n")
    return sha(p)


def name(p):
    return p.relative_to(ROOT).as_posix()


def imports():
    todo=[SRC]
    seen=set()
    while todo:
        p=todo.pop()
        if p in seen:
            continue
        seen.add(p)
        for n in ast.walk(ast.parse(p.read_text(encoding="utf-8-sig"))):
            mods=[i.name for i in n.names] if isinstance(n,ast.Import) else [n.module] if isinstance(n,ast.ImportFrom) else []
            for token in mods:
                if token and (f:=ROOT/"scripts"/(token.split(".")[0]+".py")).exists():
                    todo.append(f)
    return {name(f) for f in seen if f!=SRC}


def case(c,ledger,manifest,b0,b1,native,ep,delta):
    ident=c["case_id"]
    sid=c["checkpoint"]
    assert ident==ledger["case_id"]==f"{ep['endpoint_id']}_{sid}"
    step=ep["small_records"][int(sid[-2:])-1]
    assert c["step"]==step==manifest["step"]==b0["step"]
    assert c["record_id"]==step["record_id"]==ledger["record_id"]==b1["record_id"]
    assert ledger["fresh_cache_and_preconditioner_needed"] is True and ledger["pilot_cached_and_geometry_coverage"] is None
    for k in ("field","verification"):
        v=c[k+"_sha256"]
        assert c[k]==manifest[k]==ledger[k]==b0[k]
        assert v==ledger[k+"_sha256"]==manifest[k+"_sha256"]==b0[k+"_sha256"]==native["all_hashes"][c[k]]==sha(ROOT/c[k])
    state=read(ROOT/c["field"])
    obs=read(ROOT/c["verification"])
    material=c["material"]
    p=state["parameters"]
    assert material==ep["materials"]["small"] and p==material["parameters"]
    assert state["reached_target"] is state["converged"] is True
    assert c["parameter_order"]==[*MAT,"common_gain"]
    assert {k:Q(v) for k,v in p.items()}=={k:Q(v) for k,v in material["parameter_storage_exact"].items()}
    assert Q(p["mu"])==1 and Q(p["bulk"])==5000
    gain=Q(material["gain_exact"])
    d=Q(step["dose_exact"])
    exact=gain*d
    rounded=Q(float(exact))
    assert gain>0 and d>0 and d==Q(ledger["dose_exact"])
    assert exact==Q(ledger["requested_drive_exact"])==Q(step["requested_strength_exact"])
    assert rounded==Q(state["activation_strength"])==Q(step["solver_strength_binary64_exact"])!=exact
    assert {k:Q(v) for k,v in c["drive"].items()}=={"optical_dose_exact":d,"common_gain_exact":gain,
        "activation_strength_exact":exact,"activation_strength_binary_exact":rounded,"binary_input_minus_exact_strength":rounded-exact}
    closed={k:[Q(p[k]),Q(p[k])] for k in MAT}|{"common_gain":[gain,gain]}
    assert {k:list(map(Q,v)) for k,v in c["closed_parameter_cell_exact"].items()}==closed
    ref=[list(map(Q,v)) for v in state["reference_nodes"]]
    free=[2*i+k for i,(r,z) in enumerate(ref) for k in (0,1)
          if r!=Q(p["outer_radius"]) and (r!=0 or k==1)]
    assert free==c["free_map"] and len(free)==len(state["free_displacement"])==611
    u=dict(zip(free,map(Q,state["free_displacement"]),strict=True))
    xi=Q(0)
    for i,(a,b) in enumerate(zip(ref,state["current_nodes"],strict=True)):
        for k,(x,y) in enumerate(zip(a,b,strict=True)):
            v=x+u.get(2*i+k,0)
            assert Q(float(v))==Q(y)
            xi=max(xi,abs(v-Q(y)))
    assert xi==Q(c["xi_dimensionless_exact"])==Q(b1["binary_gap_exact"])==Q(ledger["node_addition_gap_exact"])
    assert len(ref)==325 and len(state["triangles"])==144
    assert [i for i,r in enumerate(state["activation_elements"]) if r]==c["active_elements"]
    assert len(c["active_elements"])==24
    assert c["observer_variant"]==obs["variant"]=="current_disk__spatial_max"
    assert Q(c["frozen_conditional_allowance_um_exact"])==delta
    D=b0["stored_face_derivatives"]
    a,b=(Q(D[k]["C_um_per_dimensionless_coordinate_exact"]) for k in ("top","bottom"))
    coeffs=[a+b,a+b,a,b]
    used={}
    for k,A in zip(OBS,coeffs,strict=True):
        y=Q(step["observed_exact"][k])
        l,r=map(Q,obs["observer_intervals_um_exact"][k])
        assert l<=r
        margin=delta-max(abs(y-l),abs(y-r))
        assert margin==Q(ledger["stored_primary_observer_margins_um_exact"][k])==Q(b0["stored_residuals"][k]["margin_um_exact"])>0
        used[k]=margin/(2*A)
    radial=min(Q(D[k]["radial_min_um_exact"]) for k in ("top","bottom"))/(16*Q(67,2))
    g=min(*used.values(),radial)-xi
    assert g==Q(c["e_goal_dimensionless_exact"])==Q(b1["e_goal_dimensionless_exact"])==Q(ledger["observer_goal_exact"])
    assert 0<RADII[-1]<RADII[0]<g
    return {"record_id":c["record_id"],"case_id":ident,"field_sha256":c["field_sha256"],
            "exact_requested_drive":str(exact),"observer_goal_exact":str(g),"observer_goal_display":float(g),
            "checked_primary_inequalities":len(used)}


def rejects(fn):
    try:
        fn()
    except ValueError:
        return True
    return False


def main():
    assert all(sha(f)==v for f,v in PINS.items())
    p,pf,idx,ld,m,old=map(read,(PROTO,PRE_FLIGHT,INDEX,LEDGER,M,FOUR))
    assert pf["status"]=="PASS" and len(pf["checks"])==10 and all(x is True for x in pf["checks"].values())
    assert pf["source_sha256"]==p["source_sha256"]==PINS[SRC]==sha(PRE/"transfer_ambiguity_remaining30_cached_frozen.py")
    assert pf["protocol_sha256"]==PINS[PROTO] and pf["wave_index_sha256"]==PINS[INDEX]
    assert all(type(pf[k]) is int and pf[k]==0 for k in (
        "models_constructed","preconditioners_recovered","binary_caches_frozen","actual_AD_sites_assembled","new_equilibrium_solves"))
    assert not (PRE.parent/"actual").exists()
    hashes=p["source_closure_sha256"]|p["input_sha256"]|ld["readback_sha256"]
    assert len(ld["readback_sha256"])==131
    assert imports()<=set(p["source_closure_sha256"])
    assert all(sha(ROOT/f)==h for f,h in hashes.items())
    assert list(map(Q,p["finite_radius_ladder_exact"]))==RADII
    assert p["site_counts"]=={"bulk":9216,"surface":288}
    assert p["precision_bits"]==256 and p["torch_cpu_threads"]==1
    assert p["prior_case_cache_or_R_reused"] is False
    assert p["requires_independent_actual_review_for_every_current_point"] is True
    assert old["root_certificates_cached_quadrature_only"]==4 and m["new_cached_local_root_count"]==2
    coverage={r["pilot_id"]:r["receipt_sha256"] for r in old["field_results"]}|{
              r["case_id"]:r["native_proof_sha256"] for r in m["case_results"]}
    assert p["covered_cached_cases_excluded"]==coverage and len(coverage)==6
    keys={r["case_id"]:r for r in ld["outcome_states"]}
    assert len(keys)==36 and set(keys)-set(coverage)==set(p["cases"])
    expected={}
    for wid,g in zip(WAVES,ld["recommended_batch_waves"],strict=True):
        expected[wid]=[r for r in g["case_ids"] if r not in coverage]
    assert expected==p["waves"]
    assert sorted(r for v in expected.values() for r in v)==sorted(p["cases"])
    assert [len(v) for v in expected.values()]==[2]+[4]*7
    assert [p["cases"][k]["record_id"] for k in expected["wave_01"]]==["optoGEF2_small_col_Q"]*2
    lens=Counter(c["record_id"] for c in p["cases"].values())
    assert len(lens)==15 and set(lens.values())=={2}
    assert p["numerical_coverage"]["associated_primary_small_inequalities"]==120
    s=read(PRIMARY/"preparation/protocol.json")
    native=read(PRIMARY/"independent_actual_preparation/actual_replay/receipt.json")
    mf=read(PRIMARY/"field_manifest.json")
    bs=read(DESIGN/"budgets_v1/actual/results.json")
    br=read(DESIGN/"budgets_v1/independent_actual/receipt.json")
    assert len(mf["fields"])==native["checkpoints"]==44 and mf["independent_replay_observed_os_exit_code"]==0
    delta=Q(s["selection"]["conditional_source_allowance_um_exact"])
    key=lambda r:f"{r['endpoint_id']}_{r['checkpoint']}"
    tables=[{key(x):x for x in a} for a in (mf["fields"],bs["field_budgets"],br["fields"])]
    eps={x["endpoint_id"]:x for x in s["endpoints"]}

    docs=[]
    available=[]
    assert idx["source_sha256"]==PINS[SRC] and idx["protocol_sha256"]==PINS[PROTO] and idx["total_unique_cases"]==30
    for w,entry in zip(WAVES,idx["waves"],strict=True):
        assert w==entry["wave_id"] and entry["file"]==w+".json" and entry["ordered_case_ids"]==expected[w]
        f=PRE/entry["file"]
        assert sha(f)==entry["sha256"]
        z=read(f)
        assert z["ordered_case_ids"]==expected[w] and z["wave_id"]==w and z["source_sha256"]==PINS[SRC]
        assert z["protocol_sha256"]==PINS[PROTO] and type(z["producer_actual_work_in_preparation"]) is int and z["producer_actual_work_in_preparation"]==0
        for ident,anchor in zip(z["ordered_case_ids"],z["case_anchors"],strict=True):
            c=p["cases"][ident]
            assert all(c[k]==v for k,v in anchor.items())
            available.append(case(c,keys[ident],*[d[ident] for d in tables],native,eps[c["endpoint_id"]],delta))
        for pair in range(0,len(z["ordered_case_ids"]),2):
            a,b=(p["cases"][k] for k in z["ordered_case_ids"][pair:pair+2])
            assert (a["endpoint_id"],b["endpoint_id"])==EP and a["record_id"]==b["record_id"]
            assert a["drive"]["optical_dose_exact"]==b["drive"]["optical_dose_exact"]
        docs.append((f,z,entry))
    assert len(available)==len({x["field_sha256"] for x in available})==30
    assert sum(r["checked_primary_inequalities"] for r in available)==120
    sys.path.insert(0,str(ROOT/"scripts"))
    import transfer_ambiguity_remaining30_cached as driver
    def stop(*_a,**_b):
        raise AssertionError("actual work forbidden")
    with (patch.object(driver,"actual_case",side_effect=stop),patch.object(driver,"worker",side_effect=stop),
          patch.object(driver,"launch",side_effect=stop),patch.object(subprocess,"Popen",side_effect=stop)):
        assert driver.prepare_protocol()==p
        for f,doc,entry in docs:
            assert doc==driver.wave_document(p,entry["wave_id"],PINS[PROTO])
        def mock(_,r):
            return {"radius_requested_exact":str(r),"strict_selfmap_contraction_and_observer_radius":r==RADII[3]}
        hist,ok,back=driver.run_radius_ladder(RADII[0],mock)
        assert hist[0]["status"]=="SKIPPED_OUTSIDE_OBSERVER_GOAL"
        assert Q(ok["radius_requested_exact"])==RADII[3]
        assert Q(back["radius_requested_exact"])==RADII[4] and back["strict_selfmap_contraction_and_observer_radius"] is False
        full,absent,other=driver.run_radius_ladder(Q(1),lambda i,r:mock(i,r)|{
               "strict_selfmap_contraction_and_observer_radius":False})
        assert len(full)==len(RADII) and absent is other is None
        wrong={"case_id":"test","pid":10}
        assert rejects(lambda:driver.admit_cached_worker(False,wrong|{"status":"FINISHED"},wrong,wrong))
        assert rejects(lambda:driver.admit_cached_worker(0,wrong|{"status":"FINISHED","pid":11},wrong,wrong))
        assert rejects(lambda:driver.admit_wave_review({"status":driver.REVIEW_STATUS,"checks":{"x":1}},
               "wave_01",{},expected["wave_01"]))

    assert not (PRE.parent/"actual").exists()
    checks={"verified_131_ledger_pins_and_transitive_local_source_dependencies":True,
       "exact_independent_30_maps_material_cells_drives_and_120_observation_inequalities":True,
       "original_four_and_two_M_excluded_once_and_eight_disjoint_waves_verified":True,
       "fresh_binary_R_and_cache_explicit_exact_force_drive_source_review":True,
       "radius_ladder_preserves_unresolved_smaller_backup_without_promotion":True,
       "serial_os_exit_and_exclusive_per_wave_gate_source_review":True,
       "no_actual_models_cache_AD_solve_or_spatial_work":True}
    h=save(HERE/"support.json",{"status":"PASS","checks":checks,"source_sha256":sha(Path(__file__)),
       "source_and_input_hashes":hashes|{name(f):s for f,s in PINS.items()},
       "rows":available,"eight_wave_ids":list(WAVES),"new_model_cache_AD_or_solve_work":0,
       "scope":"Prelaunch review of independent fixed-parameter cached quadrature attempts only."})
    output={}
    for f,doc,entry in docs:
        wid=entry["wave_id"]
        g=HERE/f"{wid}_review.json"
        digest=save(g,{"status":"PASS_INDEPENDENT_REMAINING30_CACHED_WAVE_PRELAUNCH","checks":checks,
           "wave_id":wid,"ordered_case_ids":doc["ordered_case_ids"],"source_sha256":PINS[SRC],
           "protocol_sha256":PINS[PROTO],"preflight_sha256":PINS[PRE_FLIGHT],"wave_index_sha256":PINS[INDEX],
           "wave_manifest_sha256":entry["sha256"],"coverage_ledger_sha256":PINS[LEDGER],
           "verified_two_M_coverage_sha256":PINS[M],"first_four_cached_coverage_sha256":PINS[FOUR],
           "shared_review_file":name(HERE/"support.json"),"shared_review_sha256":h})
        output[wid]={"gate_file":name(g),"sha256":digest,"manifest_sha256":entry["sha256"]}
    index=save(HERE/"review_index.json",{"status":"PASS_ALL_EIGHT_DISJOINT_WAVE_GATES",
            "all_30_cases_verified_cold":True,"all_120_small_observer_inequalities_checked":True,
            "support_sha256":h,"waves":output})
    print(index,h,json.dumps(output),flush=True)


if __name__=="__main__":
    main()
