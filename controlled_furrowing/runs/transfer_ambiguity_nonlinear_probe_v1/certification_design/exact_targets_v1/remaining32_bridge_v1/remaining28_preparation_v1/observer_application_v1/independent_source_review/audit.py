"""Cold independent readback for later-wave observer transfer.

Reads archived fields and frozen proof schemas. Uses only inert proof and
process fixtures for later-wave replay. Never opens a later-wave native result.
"""
from collections import Counter
from copy import deepcopy
from fractions import Fraction as Q
from hashlib import sha256
import importlib.util
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

HERE = Path(__file__).resolve().parent
APP = HERE.parent
PINS = {
    "apply.py": "df0ae901ad4c5bc9c976002c4c62e11e46a5a41b33c747779dfca4fa22871e6b",
    "protocol.json": "d20799dc6c69d9f8ab4f36e4d9bd085b2f48b5abf68d821296c3a9ddca74231e",
    "preflight.json": "e4978572f081a662adfbd77789171fe66da40365df47c66054f7d90dd45ff357",
}
L = Q(67, 2)


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_bytes())


def load_source():
    for name, h in PINS.items():
        assert digest(APP / name) == h
    spec = importlib.util.spec_from_file_location("_cold_later_obs", APP / "apply.py")
    A = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(A)
    return A


def qpin(pair):
    return tuple(Q(z["numerator"], z["denominator"]) for z in pair)


def reject(fn):
    try:
        fn()
    except (KeyError, ValueError):
        return
    raise AssertionError("bad fixture accepted")


def cold_rebuild(A, M, case, target, p, docs, allowance, source_links):
    # Reconstruct from explicit archive paths rather than trusting the adapter binding.
    endpoint, slot = case.rsplit("_small_", 1)
    checkpoint = "small_" + slot
    item, = [t for t in docs[M.GEOMETRY]["trials"] if t["case_id"] == case]
    bud, = [b for b in docs[M.BUDGET/"actual/results.json"]["field_budgets"]
            if (b["endpoint_id"], b["checkpoint"]) == (endpoint, checkpoint)]
    prior, = [t for t in docs[M.CACHED/"results.json"]["cases"] if t["case_id"] == case]
    sweep = M.SECONDARY/"actual"/endpoint/(checkpoint+".json")
    links = p["source_and_archive_sha256"]
    contract, box, sec = [read(M.ROOT/f) for f in (item["model_contract"], item["trial_box"], M.rel(sweep))]
    state, verify = read(M.ROOT/sec["field"]), read(M.ROOT/sec["verification"])
    assert all(links[f] == digest(M.ROOT/f) for f in
               (item["model_contract"], item["trial_box"], M.rel(sweep), sec["field"], sec["verification"]))
    assert (contract["pilot_id"], contract["record_id"], contract["physical_source_record_id"]) == (
        case, target["record_id"], sec["record_id"])
    assert sec["step"] == bud["step"] and bud["record_id"] == sec["record_id"] == prior["record_id"]
    assert sec["field_sha256"] == bud["field_sha256"] == target["original_observer_field_sha256"]
    assert sec["verification_sha256"] == bud["verification_sha256"] == contract["archived_verification_sha256"]
    drive = contract["exact_requested_drive"]
    assert Q(drive["common_gain_exact"])*Q(drive["optical_dose_exact"]) == Q(drive["activation_strength_exact"]) == Q(target["exact_drive"])
    assert Q(drive["optical_dose_exact"]) == Q(sec["step"]["dose_exact"]) == Q(box["dose_exact"])
    x = [Q(z) for node in state["reference_nodes"] for z in node]
    stored_nodes = [Q(z) for node in state["current_nodes"] for z in node]
    rim = Q(state["parameters"]["outer_radius"])
    fixed = {i for i in range(650) if x[i//2*2] == rim or (i%2 == 0 and x[i] == 0)}
    free = [i for i in range(650) if i not in fixed]
    assert free == box["free"] and len(free) == 611
    assert len(state["triangles"]) == 144 and state["triangles"] == box["triangles"]
    exact_y = x.copy()
    for i, u in zip(free, state["free_displacement"], strict=True):
        exact_y[i] += Q(u)
    center = [Q(z) for node in box["current_centers_exact"] for z in node]
    binary_gap = max(abs(u-v) for u,v in zip(stored_nodes, exact_y, strict=True))
    shift = max(abs(u-v) for u,v in zip(center, exact_y, strict=True))
    assert all(Q(float(y)) == u for y,u in zip(exact_y, stored_nodes, strict=True))
    assert all(center[i] == exact_y[i] == stored_nodes[i] for i in fixed)
    assert binary_gap == Q(target["original_addition_gap_exact"])
    assert shift == Q(target["archived_center_shift_exact"])
    roi = (Q(0), L*Q(state["parameters"]["patch_radius"]))
    packed = sec["observer_contract"]
    assert qpin(packed["current_radius_ROI_um"]) == roi == tuple(map(Q,bud["fixed_roi_um_exact"])) == tuple(map(Q,prior["ROI_um_exact"]))
    assert qpin([packed["reference_apical_um"],packed["reference_basal_um"]]) == (L, Q(0))
    assert roi[1] <= L*rim and roi[1] > 0
    r=A.RADIUS
    e, z = r/9, Q(3, 7)
    synthetic_eta = e*(1-z)
    goal=Q(target["observer_goal_exact"])
    eps=L*(e+shift+binary_gap)
    face_bounds=[]
    for face in ("top","bottom"):
        edge_ids=bud["center_and_model"]["fixed_trace_endpoints"][face]["trace_nodes"]
        assert len(edge_ids)==12 and all(a[-1]==b[0] for a,b in zip(edge_ids,edge_ids[1:]))
        assert stored_nodes[2*edge_ids[0][0]]==0 and stored_nodes[2*edge_ids[-1][-1]]==rim
        a_rad=[]; a_slope=[]
        for a,b,c in edge_ids:
            for k in (0,1):
                f0,fm,f1 = (L*stored_nodes[2*i+k] for i in (a,b,c))
                # Polynomial in t: f0 + (4*fm-3*f0-f1)*t + 2*(f0-2*fm+f1)*t**2.
                d0 = 4*fm-3*f0-f1
                d1 = d0+4*(f0-2*fm+f1)
                if k==0: radial=min(d0,d1)
                else: slope=max(abs(d0),abs(d1))/radial
            a_rad.append(radial); a_slope.append(slope)
        radial, slope = min(a_rad),max(a_slope)
        assert radial == Q(bud["stored_face_derivatives"][face]["radial_min_um_exact"]) > 8*eps
        assert slope == Q(bud["stored_face_derivatives"][face]["global_stored_slope_cap_exact"])
        face_bounds.append(Q(5,4)*(1+slope)*eps)
    stand=shift+e<goal
    leaf={"target":target,"status":A.ROOT_PASS,"root_asserted":True,"radius_exact":str(r),
          "arithmetic":{"eta_upper_exact":str(synthetic_eta),"z_upper_exact":str(z),"RF_lower_exact":"0",
                        "strict_selfmap":True,"strict_cube_exclusion":False,
                        "selfmap_lhs_upper_exact":str(synthetic_eta+z*r),"posterior_error_upper_exact":str(e)},
          "observer_and_optional_minimum":{"archived_state_budget_proved":stand,"archived_shift_plus_posterior":str(shift+e),
                         "physical_nodal_error_um_upper_exact":str(eps)}}
    fixture={"status":A.ROOT_PASS,"exact_spatial_root_certified":True,"shift_inclusive_observer_budget_proved":stand,
             "case_results_sha256":"a"*64}
    frozen=A.bind_archive(M,docs,source_links,case,target)
    out=A.apply_case(M,case,target,fixture,leaf,frozen,allowance)
    assert Q(out["physical_nodal_error_um_exact"])==eps and tuple(map(Q,out["ROI_um_exact"]))==roi
    compare_count=0
    for item,(variant,obs) in zip(out["comparisons"],[(v,o) for v in M.VARIANTS for o in M.OBS],strict=True):
        assert (item["variant"], item["observable"]) == (variant,obs)
        clo,chi=map(Q, sec["variant_observables_um_exact"][variant][obs])
        sto=Q(bud["step"]["observed_exact"][obs])
        if variant==M.VARIANTS[0]:
            assert (clo,chi)==tuple(map(Q,verify["observer_intervals_um_exact"][obs]))
        a,b=face_bounds
        err=a if obs=="apical_max_shortening_um" else b if obs=="basal_max_shortening_um" else a+b
        residual_lower=max(clo-err-sto,sto-chi-err,0)
        residual_upper=max(abs(clo-err-sto),abs(chi+err-sto))
        want="IN" if residual_upper<=allowance else "OUT" if residual_lower>allowance else "UNRESOLVED"
        assert item["standing"] == want
        assert tuple(map(Q,item["exact_spatial_root_observer_interval_um_exact"])) == (clo-err,chi+err)
        assert tuple(map(Q,item["absolute_residual_um_exact"])) == (residual_lower,residual_upper)
        compare_count+=1
    bad_archive = (deepcopy(frozen[0]),*frozen[1:])
    bad_archive[0]["stored_face_derivatives"]["top"]["radial_min_um_exact"]=str(8*eps)
    for rr,ll,record,reason in (
        (fixture, leaf, bad_archive, "STRICT_RADIAL_GRAPH_GUARD_NOT_ESTABLISHED"),
        ({**fixture,"status":A.UNRESOLVED,"exact_spatial_root_certified":False,"shift_inclusive_observer_budget_proved":False},
         {**leaf,"status":A.UNRESOLVED,"root_asserted":False,
          "arithmetic":{**leaf["arithmetic"],"eta_upper_exact":str((1-z)*r),"selfmap_lhs_upper_exact":str(r),
                        "strict_selfmap":False,"posterior_error_upper_exact":None}},
         frozen,"ROOT_NOT_CERTIFIED_IN_THIS_CUBE"),
    ):
        blocked=A.apply_case(M,case,target,rr,ll,record,allowance)
        assert len(blocked["comparisons"])==16 and blocked["graph_premises"] is None
        assert all(v["standing"]==reason and v["exact_spatial_root_observer_interval_um_exact"] is None
                   and v["absolute_residual_um_exact"] is None for v in blocked["comparisons"])
    reject(lambda:A.locate(M,fixture,{**leaf,"target":{**target,"original_addition_gap_exact":str(binary_gap+r)}}))
    return {"case_id":case,"archived_record_id":sec["record_id"],"field_sha256":sec["field_sha256"],
            "archived_shift_exact":str(shift),"original_addition_gap_exact":str(binary_gap),
            "test_only_reconstructed_rows":compare_count}




def schema(A, M, wave, profile):
    """Independent 11-file fixture with real digest I/O and both metadata phases."""
    G=A.G
    def write(path,obj):
        path.parent.mkdir(parents=True,exist_ok=True)
        path.write_text(json.dumps(obj),encoding="utf8")
        return digest(path)
    with TemporaryDirectory(prefix="inert_",dir=HERE) as scratch:
        wave_dir=Path(scratch)
        red=wave_dir/"global_reduction_v1"
        ps=red/"producer"
        nv=red/"independent"
        order=profile["case_order"]
        phase={**A.GLOBAL_PHASE,"wave_id":wave,"pid":321,
               "hull_admission_sha256":"a"*64,"raw_wave_admission_sha256":"b"*64,
               "raw_root_observed_exit_sha256":"c"*64}
        def native(root, stamp, status, independent, bump=None):
            case_index=[]
            pins={}
            def retain(path,obj):
                h=write(path,obj)
                pins[M.rel(path)]=h
                return h
            for case in order:
                p=root/case
                op=retain(p/"operands.json",{"inert":case})
                item={**stamp,"case_id":case,"target":profile["raw_worker_targets"][case],
                      "operands_sha256":op,"status":A.ROOT_PASS if independent else "INERT_PRODUCER",
                      "radius_exact":str(A.RADIUS),"full_five_gate_complete":False,
                      "archived_observer_field_replaced":False}
                if bump and case==order[-1]:
                    item.update(bump)
                dat=retain(p/"results.json",item)
                case_index.append({"case_id":case,"status":item["status"],"operands_sha256":op,"case_results_sha256":dat})
            s=retain(root/"started.json",{"status":"RUNNING",**stamp})
            q=retain(root/"results.json",{"status":status,**stamp,"cases":case_index,
                      "signed_hull_components":672,"new_force_AD_or_spatial_integrals":0,"full_five_gate_complete":False})
            e=retain(root/"terminal.json",{"status":"FINISHED",**stamp,"exit_code":0,"results_sha256":q})
            assert len(pins)==11
            return {**stamp,"observed_exit_code":0,"unified_exec_chunk_id":"fixture",
                    "started_sha256":s,"terminal_sha256":e,"results_sha256":q,"sha256":pins}
        p=native(ps,{**phase,"mode":"producer"},G.PROD_STATUS,False)
        ph=write(ps/"root_observed_exit.json",{**p,"status":G.PROD_OS})
        stamp={**phase,"pid":654,"mode":"independent","producer_observed_os_sha256":ph}
        row=native(nv,stamp,G.INDEP_STATUS,True)
        def run(row):
            f=nv/"root_observed_exit.json"
            h=write(f,{"status":A.OBSERVED,**row})
            with patch.object(A.X,"wave_folder",return_value=wave_dir):
                return A.native_gate(M,wave,profile,h,{})
        assert tuple(run(row))==tuple(order)
        reject(lambda:run({**row,"wave_id":"foreign"}))
        reject(lambda:run({**row,"status":"PASS_ROOT_OBSERVED_FIRST_WAVE_INDEPENDENT_GLOBAL_EXIT"}))
        reject(lambda:run({**row,"pid":654.0}))
        reject(lambda:run({**row,"observed_exit_code":False}))
        reject(lambda:run({**row,"unified_exec_chunk_id":None}))
        reject(lambda:run({**row,"sha256":{**row["sha256"],"foreign":"9"*64}}))
        reject(lambda:run(native(nv,stamp,G.INDEP_STATUS,True,{"pid":654.0})))
        reject(lambda:run(native(nv,stamp,G.INDEP_STATUS,True,{"target":profile["raw_worker_targets"][order[0]]})))
        reject(lambda:run(native(nv,stamp,G.INDEP_STATUS,True,{"source_sha256":"7"*64})))
        row=native(nv,stamp,G.INDEP_STATUS,True)
        write(nv/order[-1]/"operands.json",{"mutated":True})
        reject(lambda:run(row))
        row=native(nv,stamp,G.INDEP_STATUS,True)
        write(ps/order[0]/"operands.json",{"mutated":True})
        reject(lambda:run(row))


def main():
    A=load_source()
    M,profiles,docs,allowance,links=A.source_inputs()
    p=read(APP/"protocol.json")
    f=read(APP/"preflight.json")
    assert p["global_phase"]==A.GLOBAL_PHASE and p["source_sha256"]==f["source_sha256"]==PINS["apply.py"]
    assert f["protocol_sha256"]==PINS["protocol.json"] and all(x is True for x in f["checks"].values())
    assert list(profiles)==list(A.X.WAVES)==[f"wave_{i:02}" for i in range(2,9)]
    assert p["ordered_cases"]=={w:d["case_order"] for w,d in profiles.items()}
    assert allowance==Q(p["allowance_um_exact"])>0
    evidence=[]
    for wave, profile in profiles.items():
        assert len(profile["case_order"])==4
        evidence.extend(cold_rebuild(A,M,case,profile["raw_worker_targets"][case],
                                     p,docs,allowance,links) for case in profile["case_order"])
        schema(A,M,wave,profile)
    assert len({r["case_id"] for r in evidence})==28
    assert sum(r["test_only_reconstructed_rows"] for r in evidence)==448
    assert Counter(bool(Q(r["archived_shift_exact"])) for r in evidence)=={False:28}
    assert M.compare((Q(0),Q(1)),Q(2),Q(0),Q(2))[2]=="IN"
    assert M.compare((Q(0),Q(1)),Q(2),Q(0),Q(1))[2]=="UNRESOLVED"
    assert M.compare((Q(0),Q(1)),Q(2),Q(0),Q(9,10))[2]=="OUT"
    assert links==p["source_and_archive_sha256"]
    assert all(digest(M.ROOT/path)==h for path,h in links.items())
    out={
        "status":"PASS_SYNTHETIC_AND_FROZEN_ARCHIVE_REPLAY",
        "source_sha256":PINS["apply.py"],"protocol_sha256":PINS["protocol.json"],"preflight_sha256":PINS["preflight.json"],
        "later_waves_release_sha256":A.G.RELEASE,"profiles_sha256":A.X.PROFILE_SHA,
        "global_source_sha256":A.GLOBAL_PHASE["source_sha256"],"global_source_review_sha256":A.GLOBAL_PHASE["source_review_sha256"],
        "approved_waves":list(A.X.WAVES),"replayed_archive_cases":evidence,"archive_pins_rehashed":len(links),
        "fixture_groups":{"native_metadata_waves":7,"synthetic_valid_comparisons":448,"graph_boundary_blocked_rows":448,
                          "selfmap_boundary_blocked_rows":448,"materialized_actual_observer_comparisons":0},
        "checks":{
            "all_seven_wave_profiles_and_28_frozen_fields_with_current_pins":True,
            "genuine_external_OS_gate_requiring_all_eleven_current_native_files_and_separate_producer":True,
            "fixture_rejects_wrong_wave_phase_types_operand_and_target":True,
            "pure_proven_P2_helpers_use_same_fixed_free_map_trace_span_ROI_and_references":True,
            "exact_shift_posterior_original_binary_gap_and_synthetic_64_rows_per_wave":True,
            "uncertified_root_and_graph_guard_equality_retain_only_archived_values":True,
            "allowance_equality_IN_and_OUT_strict_signed_residuals":True,
            "cold_only_no_native_later_output_and_no_new_scientific_jobs":True,
        },
        "opened_actual_later_global_or_observer_payloads":0,"new_spatial_integrations_authorized":False,
    }
    with (HERE/"audit_result.json").open("x",encoding="utf8") as f:
        json.dump(out,f,indent=2); f.write("\n")
    print(json.dumps({"status":out["status"],"audit_result_sha256":digest(HERE/"audit_result.json"),
                      "archive_pins":len(links),"profile_cases":28,"synthetic_comparisons":448}))


if __name__=="__main__":
    main()

