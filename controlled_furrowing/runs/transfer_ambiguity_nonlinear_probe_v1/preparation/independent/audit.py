"""Independent prelaunch checks for the frozen two-endpoint replay."""
import csv
from fractions import Fraction as Q
from hashlib import sha256
from itertools import product
import json
import math
from pathlib import Path
from statistics import median
import sys
from unittest.mock import patch

HERE=Path(__file__).resolve().parent
PRE=HERE.parent
ROOT=HERE.parents[3]
SOURCE=ROOT/"scripts/transfer_ambiguity_nonlinear_probe.py"
LIN=ROOT/"runs/transfer_ambiguity_linear_v1"
SMALL=ROOT/"data/derived/calibration_small_18um.csv"
LARGE=ROOT/"data/derived/spatial_inputs_only.json"
PINS={
 SOURCE:"87ca09b46f9569903276923d937bef63d00bba9ad1b99d78b7e1708168f7e631",
 PRE/"protocol.json":"9d5d38bb5f9785d40bf843dd225eb93b1d5ecbe96f20942aec0023765d119d94",
 PRE/"preflight.json":"f225616c4274a856576f587d40a93c89c64f88a9e4b3f1522d9f026434477211",
 LIN/"actual_independent/receipt.json":"8fc50b42694b3db513e6d2166c6d539fbe67c66c2b3c6895c90df12f3856d072",
}
MATERIAL=["apical_area_modulus","basal_area_modulus","apical_surface_tension","basal_surface_tension"]
OBS=["maximum_shortening_um","region_averaged_shortening_um","apical_max_shortening_um","basal_max_shortening_um"]


def sha(p):
    return sha256(p.read_bytes()).hexdigest()


def read(p):
    return json.loads(p.read_bytes())


def write(path,doc):
    with path.open("x",encoding="utf-8") as stream:
        json.dump(doc,stream,indent=2)
        stream.write("\n")


def bf(x):
    return Q.from_float(float(x))


def closed_dose(record,dose,gain):
    assert Q(record["dose_exact"])==dose
    exact=Q(record["requested_strength_exact"])
    rounded=Q(record["solver_strength_binary64_exact"])
    assert exact==dose*gain and rounded==bf(exact)
    assert Q(record["binary64_rounding_difference_exact"])==rounded-exact


def transformed(base,step,target,names):
    changes=dict(zip(MATERIAL,(step[0],step[1],step[2],step[3]),strict=True))
    assert list(map(Q,target["coordinate_binary64_exact"]))==list(map(bf,step))
    assert list(map(Q,target["coordinate_rounding_error_exact"]))==[bf(x)-x for x in step]
    expq={names[i]:bf(math.exp(float(step[i]))) for i in (0,1,3,4)}
    assert {k:Q(v) for k,v in target["exp_of_binary_coordinate_binary64_exact"].items()}==expq
    expected={k:bf(value) for k,value in base["parameters"].items()}
    before={}
    for i,k in enumerate(MATERIAL):
        before[k]=expected[k]+bf(changes[k]) if i==2 else expected[k]*expq[names[i]]
        expected[k]=bf(before[k])
    assert {k:Q(v) for k,v in target["material_before_final_rounding_exact"].items()}==before
    assert {k:Q(v) for k,v in target["material_final_rounding_error_exact"].items()}=={k:expected[k]-v for k,v in before.items()}
    assert {k:bf(v) for k,v in target["parameters"].items()}==expected
    assert {k:Q(v) for k,v in target["parameter_storage_exact"].items()}==expected
    assert target["parameters"]["apical_surface_tension"]==0
    assert all(expected[k]>0 for k in ("apical_area_modulus","basal_area_modulus","basal_surface_tension"))
    gain=Q(base["gain_exact"])*expq[names[-1]]
    assert gain>0 and Q(target["gain_exact"])==gain
    return gain


def main():
    assert all(sha(p)==h for p,h in PINS.items())
    p=read(PRE/"protocol.json")
    pre=read(PRE/"preflight.json")
    assert pre["status"]=="PASS" and pre["checks"] and all(v is True for v in pre["checks"].values())
    assert pre["source_sha256"]==p["source_sha256"]==PINS[SOURCE]
    assert pre["protocol_sha256"]==PINS[PRE/"protocol.json"]
    indep=PINS[LIN/"actual_independent/receipt.json"]
    assert pre["linear_actual_review_sha256"]==p["authority"]["independent_linear_actual_review_sha256"]==indep
    assert pre["nonlinear_jobs_launched"]==0 and pre["heldout_response_tables_read"] is False
    assert not (PRE/"failure.json").exists() and not (PRE.parent/"actual").exists()
    assert sha(PRE/"transfer_ambiguity_nonlinear_probe_frozen.py")==PINS[SOURCE]
    pins={**p["input_hashes"],**p["source_hashes"]}
    assert all(sha(ROOT/k)==v for k,v in pins.items())
    expected_data={
        LIN/"actual/results.json",LIN/"actual_independent/receipt.json",LIN/"actual/terminal.json",
        LIN/"actual/started.json",LIN/"actual/os_exit_observation.json",SMALL,LARGE,
        ROOT/"runs/transfer_pair_pilot_v1/preparation/protocol.json",
        ROOT/"runs/material_cohort_direct_v1/protocol.json",
    }
    assert set(p["input_hashes"])=={x.relative_to(ROOT).as_posix() for x in expected_data}
    assert all(p[k] is False for k in ("heldout_response_tables_read","retuning_or_fitting","source_witness_certified",
                                      "exact_integral_FE_root_certified","continuum_accuracy_certified"))
    result=read(LIN/"actual/results.json")
    review=read(LIN/"actual_independent/receipt.json")
    tr,st,os=[read(LIN/"actual"/name) for name in ("terminal.json","started.json","os_exit_observation.json")]
    assert tr["pid"]==st["pid"]==os["worker_terminal_pid"] and tr["status"]=="FINISHED"
    assert tr["exit_code"]==os["observed_exit_code"]==0 and os["worker_terminal_sha256"]==sha(LIN/"actual/terminal.json")
    assert tr["results_sha256"]==os["results_sha256"]==sha(LIN/"actual/results.json")
    assert all(sha(ROOT/f)==h for f,h in review["bindings"].items())
    assert review["status"]=="PASS_INDEPENDENT_ACTUAL_LOCAL_AFFINE_DIAGNOSTIC"
    assert review["source_compatible_witness"] is False
    assert result["heldout_response_tables_read"] is False
    (row,)=[r for r in result["rows"] if r["mode"]=="source_residual_allowance" and Q(r["budget_exact"])==Q(1,100)]
    chosen=[r for r in review["LP_endpoints"] if r["mode"]=="source_residual_allowance" and Q(r["budget_exact"])==Q(1,100)]
    assert len(chosen)==2
    assert p["coordinates"]==result["columns"] and p["coordinates"][2]=="apical_tension_additive"
    c=list(map(Q,result["rounded_affine_contrast_gradient_exact"]))
    d0=Q(result["source_baseline_residual_um_exact"])
    selection=p["selection"]
    assert selection["mode"]=="source_residual_allowance" and Q(selection["increment_um_exact"])==Q(1,100)
    assert Q(selection["scale_each_frozen_step_exact"])==Q(1,2)
    assert Q(selection["conditional_source_allowance_um_exact"])==Q(row["allowance_um_exact"])==d0+Q(1,100)
    assert Q(selection["affine_baseline_source_residual_um_exact"])==d0
    with SMALL.open(encoding="utf-8-sig",newline="") as stream:
        small=list(csv.DictReader(stream))
    small.sort(key=lambda x:(Q(x["activation_magnitude"]),x["record_id"]))
    large=json.loads(LARGE.read_bytes(),parse_float=str)
    assert len(small)==17 and len(large)==31
    assert len({r["record_id"] for r in small})==17 and len({r["record_id"] for r in large})==31
    assert all(set(x)=={"record_id","reported_diameter_um","activation_magnitude"} for x in large)
    largest=max(Q(x["reported_diameter_um"]) for x in large)
    m=median([Q(r["activation_magnitude"]) for r in small])
    selected=min((r for r in large if Q(r["reported_diameter_um"])==largest),
                 key=lambda x:(abs(Q(x["activation_magnitude"])-m),x["record_id"]))
    dose=Q(selected["activation_magnitude"])
    assert p["matched_input"]["large_record_id"]==selected["record_id"]=="fig4f_col_I"
    assert Q(p["matched_input"]["matched_dose_exact"])==dose
    assert Q(p["large_observable"]["shoulder_radius_um_exact"])*2==largest
    assert p["observable_variant"]=="current_disk__spatial_max"
    assert Q(p["large_observable"]["length_unit_um_exact"])==Q("33.5")
    assert p["loading"]["serial_torch_threads"]==1
    bas=p["baseline"]
    assert {k for k in bas["small"]["parameters"] if bas["small"]["parameters"][k]!=bas["large"]["parameters"][k]}=={"patch_radius"}
    pilot=read(ROOT/"runs/transfer_pair_pilot_v1/preparation/protocol.json")
    for level in ("small","large"):
        assert bas[level]==next(x for x in pilot["cases"] if x["case_id"]=="baseline_"+level)
    sh=[]
    for candidate,independent in zip(p["endpoints"],chosen,strict=True):
        k=candidate["endpoint_id"]
        assert independent["objective"]=={"max_contrast_change":"maximize_Q","max_negative_contrast_change":"minimize_Q"}[k]
        orig=list(map(Q,row[k]["step_exact"]))
        assert orig==list(map(Q,independent["step_exact"]))==list(map(Q,candidate["full_linear_step_exact"]))
        step=[x/2 for x in orig]
        assert step==list(map(Q,candidate["half_step_exact"])) and step[2]==0
        score=sum(a*b for a,b in zip(c,step,strict=True))
        assert score==Q(candidate["affine_large_Q_change_um_exact"])
        sh.append(score)
        gains=[transformed(bas[level],step,candidate["materials"][level],p["coordinates"]) for level in ("small","large")]
        assert gains[0]==gains[1]
        sp=candidate["materials"]["small"]["parameters"]
        lg=candidate["materials"]["large"]["parameters"]
        assert {k for k in sp if sp[k]!=lg[k]}=={"patch_radius"}
        assert sp["mu"]==1 and sp["bulk"]==5000 and sp["height"]==1
        gain=gains[0]
        closed_dose(candidate["small_zero"],Q(0),gain)
        for rr,ss in zip(candidate["small_records"],small,strict=True):
            assert rr["record_id"]==ss["record_id"] and rr["source_csv_dose_literal"]==ss["activation_magnitude"]
            assert set(rr["observed_exact"])==set(OBS) and all(Q(rr["observed_exact"][k])==Q(ss[k]) for k in OBS)
            assert Q(rr["reported_diameter_um_exact"])==Q(ss["reported_diameter_um"])==Q(bas["small"]["diameter_um_exact"])
            closed_dose(rr,Q(ss["activation_magnitude"]),gain)
        for i,cp in enumerate(candidate["large_checkpoints"]):
            closed_dose(cp,i*dose/3,gain)
        assert len(candidate["large_checkpoints"])==4
    assert sh[0]>0>sh[1] and Q(selection["half_pair_affine_large_Q_range_um_exact"])==sh[0]-sh[1]
    assert p["numerical_policy"]==pilot["numerical_policy"]
    for key, v,surface in (("solve_quadrature",8,12),("cross_force_quadrature",12,16)):
        policy=p["numerical_policy"][key]
        assert (policy["volume_duffy_gauss_order"],policy["surface_gauss_order"])==(v,surface)
    refs=pre["unsolved_zero_reference_receipts"]
    assert len(refs)==4
    e0,e1=(x["endpoint_id"] for x in p["endpoints"])
    for level in ("small","large"):
        a,b=(refs[f"{name}_{level}"] for name in (e0,e1))
        for k in ("reference_coordinate_sha256","connectivity_sha256","fixed_mask_sha256","activation_mask_sha256"):
            assert a[k]==b[k]
        for x in (a,b):
            assert x["reference_node_count"]==325 and x["triangle_count"]==144 and x["free_dof_count"]==611
            assert x["activation_element_count"]==24 and Q(x["exact_J_lower"])>0 and x["Q_um_exact"]==["0","0"]
    assert refs[f"{e0}_small"]["reference_coordinate_sha256"]!=refs[f"{e0}_large"]["reference_coordinate_sha256"]

    # No solves or child launch are admitted for this replay of metadata and guards.
    sys.path.insert(0,str(ROOT/"scripts"))
    import transfer_ambiguity_nonlinear_probe as m
    m.torch.set_num_threads(1)  # Match the frozen CLI runtime before comparison.
    def forbid(*args,**kwargs):
        raise AssertionError("nonlinear operations disabled in independent prelaunch")
    with (patch.object(m.pilot.AxisymmetricModel,"_attempt",side_effect=forbid),
          patch.object(m.pilot.AxisymmetricModel,"advance",side_effect=forbid),
          patch.object(m.subprocess,"Popen",side_effect=forbid)):
        assert m.prepare_protocol()==p
        for c in p["endpoints"]:
            plan=list(m.checkpoint_plan(c))
            assert [r[0] for r in plan]==["small_zero"]+[f"small_{i:02d}" for i in range(1,18)]+[f"large_{i:02d}" for i in range(4)]
            assert [role for _,_,role in plan]==["small_zero"]+["small_record"]*17+["large_zero","large_checkpoint","large_checkpoint","large_final"]
        c=p["endpoints"][0]
        predictions={r["record_id"]:{k:[r["observed_exact"][k],r["observed_exact"][k]] for k in OBS} for r in c["small_records"]}
        r=c["small_records"][-1]
        one=r["observed_exact"][OBS[0]]
        predictions[r["record_id"]][OBS[0]]=[str(Q(one)+1),str(Q(one)+1)]
        accepted=m.residual_summary(c["small_records"],predictions,Q(1))
        excessive=m.residual_summary(c["small_records"],predictions,Q(1,2))
        predictions[r["record_id"]][OBS[0]]=[str(Q(one)),str(Q(one)+2)]
        unresolved=m.residual_summary(c["small_records"],predictions,Q(1))
        assert [x["status"] for x in (accepted,excessive,unresolved)]==["ALLOWANCE_SATISFIED","ALLOWANCE_EXCEEDED","UNRESOLVED"]
        assert all(x["all68_residuals"]["constraints"]==68 for x in (accepted,excessive,unresolved))
        binds={"source_sha256":PINS[SOURCE],"protocol_sha256":PINS[PRE/"protocol.json"],
               "preflight_sha256":PINS[PRE/"preflight.json"],"linear_actual_review_sha256":indep}
        assert m.review_bindings({"status":"PASS","checks":{"fixture":True},**binds},binds)
        for k in binds:
            assert not m.review_bindings({"status":"PASS","checks":{"fixture":True},**binds,k:"x"},binds)
        for v in ({},{"fixture":"PASS"},{"fixture":1},{"fixture":False}):
            assert not m.review_bindings({"status":"PASS","checks":v,**binds},binds)
    assert not (PRE.parent/"actual").exists() and not (PRE/"failure.json").exists()
    assert all(sha(f)==h for f,h in PINS.items())
    artifact=HERE/"support.json"
    support={
        "status":"PASS_INDEPENDENT_PRELAUNCH_NONLINEAR_PAIR",
        "checks":{
            "source_protocol_preflight_linear_review_and_os_exit_frozen":True,
            "exact_independently_reviewed_row_half_steps_unchanged":True,
            "binary_exp_material_gain_and_dose_rounding_reconstructed":True,
            "same_frozen_nonradius_geometry_and_matching_reference_receipts":True,
            "all_17_small_doses_68_values_and_input_only_large_selection_reconstructed":True,
            "four_fresh_zeros_and_44_serial_checkpoint_contract_audited":True,
            "separate_Q8_Q12_stored_geometry_exact_graph_and_retained_failures_audited":True,
            "conditional_residual_outcomes_do_not_stop_field_collection":True,
            "exclusive_serial_worker_OS_pid_terminal_and_hash_binding_audited":True,
            "rejects_mixed_failed_empty_or_nonnative_review_checks":True,
        },
        "checks_scope":"Frozen-input and source audit with graph/residual fixture checks. The four zero-reference certificates are existing producer preparation checks; independent review ran no optimizer, continuation, or actual workers.",
        "primary_pair_scope":"Half steps on the 0.01 µm baseline-residual-allowance row, frozen before any actual outcomes. The source tolerance is a conditional diagnostic threshold.",
        "source_history_note":"Small historical archive informs retrospective parameters and ascending small-dose order; the new runs start from freshly solved zero for each endpoint and geometry. Finite checkpoints do not prove continuous biological or mechanically stable history.",
        "terminal_note":"Each worker is gated, must match the live parent PID, and is waited to true process exit before advancing to the next endpoint. A pre-worker admission failure is retained by the supervisor and blocks aggregation. Residual allowance is computed after all 22 states per endpoint.",
        "input_sha256":pins,
        "audit_source_sha256":sha(Path(__file__)),
        "actual_workers_started":0,"outcome_files_read":False,"source_witness_claimed":False,
    }
    write(artifact,support)
    path=PRE/"prelaunch_independent_review.json"
    write(path,{"status":"PASS","checks":support["checks"],**binds,
                "support_file":artifact.relative_to(ROOT).as_posix(),"support_sha256":sha(artifact),
                "scope":support["checks_scope"]})
    print(sha(path),sha(artifact))


if __name__=="__main__":
    main()

