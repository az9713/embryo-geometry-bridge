"""Join complete observed geometry and derivative premise bundles without numerical kernels."""
import argparse
from base64 import b64decode
from fractions import Fraction as Q
import importlib.util
import io
import json
from pathlib import Path

import numpy as np

HERE=Path(__file__).resolve().parent
BR=HERE.parent
COLD=BR/"independent_prelaunch/audit.py"
COLD_SOURCE_SHA="e77519527a7fe405b0ed868f54df823a0328839765a449af01cce9ce3cbaba8c"
COLD_SHA="89f9b5cacf89b4055fc404c00d3e3c42cf10fed20337fe455e996c2944ef73b2"
GEOMETRY_PHASE={
    "source_sha256":"6630ad02df027293b248498fefcd6c8532bc95b9919d2669a1ef096389028895",
    "protocol_sha256":"673d83339d13d65e5464d035986e1da2cbb6df4977c5187b36448391c6ac0a3c",
    "preflight_sha256":"5e5903f504dbea6a29bc6a1b0667b0a7b7f92524516df919237ecfdef5215a3d",
}
GEOMETRY_PRELAUNCH="0554a8b74023aa1823b0355c062187976665adcf13abe9f736d7e175eead800b"
RAW_GEOM_SHA="50ce32f50b5e330b9c78fc1863527c49eb6c5e4394bdb1c4e4313cea8628efe0"


def load_cold():
    from hashlib import sha256
    if sha256(COLD.read_bytes()).hexdigest()!=COLD_SOURCE_SHA:
        raise ValueError("independent cold audit changed")
    spec=importlib.util.spec_from_file_location("_independent_cold_bridge_audit",COLD)
    a=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(a)
    a.read(COLD,COLD_SOURCE_SHA,False)
    return a


def dimensions(a,x,ordered,digest,scope):
    pin=(a.G if scope=="geometry" else a.D)/"actual"
    raw=a.read(pin/"root_terminal_admission.json",digest)
    assert raw["stage"]==scope and a.zero(raw["observed_exit_code"])
    assert type(raw["pid"]) is type(raw["supervisor_pid"]) is type(raw["unified_exec_session"]) is int
    assert raw["pid"]==raw["supervisor_pid"]>0 and raw["unified_exec_session"]>0
    grp=a.sha(json.dumps(list(ordered),sort_keys=True,separators=(",",":")).encode())
    assert raw["launch_group"]==grp
    folder=pin/"launches"/grp
    group_files={k:folder/(k+".json") for k in ("started","terminal","manifest")}
    leaf_names=("launched","started","terminal","results","observed_os_exit")+(("body_geometry",) if scope=="geometry" else ())
    known={a.relative(p) for p in group_files.values()}|{
        a.relative(pin/"cases"/c/(leaf+".json")) for c in ordered for leaf in leaf_names}
    assert set(raw["sha256"])==known and len(known)==(195 if scope=="geometry" else 53)
    assert all(x[k]==h for k,h in raw["sha256"].items())
    group_docs={k:a.read(p,raw["sha256"][a.relative(p)]) for k,p in group_files.items()}
    start,term,manifest=(group_docs[k] for k in ("started","terminal","manifest"))
    assert start["status"]=="RUNNING" and term["status"]=="FINISHED" and a.zero(term["exit_code"])
    assert manifest["status"]=="COMPLETE_ACTUAL_WORK_AWAITS_INDEPENDENT_REPLAY"
    assert manifest["root_asserted"] is False and term["manifest_sha256"]==raw["sha256"][a.relative(group_files["manifest"])]
    assert all(b["pid"]==raw["pid"] and b["launch_group"]==grp and b["stage"]==scope
               and tuple(b["case_order"])==tuple(ordered) and b["bindings"]==start["bindings"] for b in (start,term,manifest))
    assert isinstance(term["workers"],list) and len({x["pid"] for x in term["workers"]})==len(ordered)
    assert term["workers"]==manifest["workers"] and [q["case_id"] for q in term["workers"]]==list(ordered)
    out={}
    for w in term["workers"]:
        cid=w["case_id"]
        p=pin/"cases"/cid
        sealed={k:raw["sha256"][a.relative(p/(k+".json"))] for k in leaf_names}
        metadata={k:a.read(p/(k+".json"),sealed[k]) for k in ("launched","started","terminal","observed_os_exit")}
        obs,end,st=(metadata[k] for k in ("observed_os_exit","terminal","started"))
        assert type(w["pid"]) is int and w["pid"]>0 and a.zero(w["observed_exit_code"]) and a.zero(obs["observed_exit_code"])
        common={"pid":w["pid"],"supervisor_pid":raw["pid"],"case_id":cid,"bindings":term["bindings"],"stage":scope,"launch_group":grp}
        assert all(all(j[k]==v for k,v in common.items()) for j in (*metadata.values(),w))
        assert end["status"]=="FINISHED" and a.zero(end["exit_code"]) and st["status"]=="RUNNING"
        assert end["results_sha256"]==w["results_sha256"]==sealed["results"]
        assert obs["terminal_sha256"]==w["terminal_sha256"]==sealed["terminal"] and w["os_exit_sha256"]==sealed["observed_os_exit"]
        out[cid]=(p,sealed,w)
    return raw,out


def replay_receipt(a,folder,observed_hash,kind):
    o=a.read(folder/"os_exit_observation.json",observed_hash)
    e=a.read(folder/"terminal.json",o["terminal_sha256"])
    n=a.read(folder/"receipt.json",o["receipt_sha256"])
    assert type(o["pid"]) is type(e["pid"]) is type(o["unified_exec_session"]) is int
    assert o["pid"]==e["pid"]>0 and o["unified_exec_session"]>0
    assert a.zero(o["observed_exit_code"]) and a.zero(e["exit_code"]) and e["status"]=="FINISHED"
    assert e["receipt_sha256"]==o["receipt_sha256"]
    if kind=="geometry":
        phase=n["source_and_phase"]
        assert o["status"]=="PASS_ROOT_OBSERVED_REMAINING32_GEOMETRY_REPLAY_EXIT"
        assert n["status"]=="PASS_INDEPENDENT_ACTUAL_REMAINING32_SINGLETON_GEOMETRY"
        assert n["independent_all_32_geometry"] is True and n["successfully_replayed_geometry_cases"]==n["total_geometry_cases"]==32
        assert n["root_asserted"] is n["exact_spatial_bridge_authorized"] is False
        assert a.zero(n["new_spatial_integral_force_Hessian_AD_or_solve_jobs"])
        assert phase["producer_root_admission_sha256"]==RAW_GEOM_SHA
        assert phase["source_review_sha256"]==GEOMETRY_PRELAUNCH
        assert all(phase[k]==v==e[k] for k,v in GEOMETRY_PHASE.items())
        for k,name in (("source_sha256","replay_geometry.py"),("protocol_sha256","protocol.json"),
                       ("preflight_sha256","preflight.json")):
            a.read(folder.parent/name,GEOMETRY_PHASE[k],False)
        p=a.read(folder.parent/"root_source_review.json",GEOMETRY_PRELAUNCH)
        assert p["status"]=="PASS_ROOT_REMAINING32_GEOMETRY_REPLAY_SOURCE" and a.truth(p["checks"])
        assert all(p[k]==v for k,v in GEOMETRY_PHASE.items())
    else:
        keys=("pid","source_sha256","protocol_sha256","preflight_sha256",
              "producer_admission_sha256","root_source_review_sha256")
        phase={k:n[k] for k in keys}
        assert observed_hash==a.PINS[folder/"os_exit_observation.json"]
        assert o["status"]=="PASS_OBSERVED_INDEPENDENT_REMAINING32_VARIATION_TRUE_OS_EXIT"
        assert n["status"]=="PASS_INDEPENDENT_ACTUAL_TEN_FRESH_CACHED_DERIVATIVES" and n["fresh_valid_derivative_count"]==10
        assert a.zero(n["spatial_derivative_root_or_integration_count"])
        assert a.PINS[a.D/"actual/root_terminal_admission.json"]==phase["producer_admission_sha256"]
    assert all(e[k]==v for k,v in phase.items()) and phase["pid"]==o["pid"]
    return n,e,o,phase


def checked_geo(a,e,native_row,geophase,leafmap,folder):
    case=e["case_id"]
    tr=e["trial"]
    p,seals,worker=leafmap[case]
    leaf_path=folder/(case+".json")
    assert native_row["receipt"]==a.relative(leaf_path)
    leaf=a.read(leaf_path,native_row["sha256"])
    assert leaf["status"]==native_row["status"]=="PASS_INDEPENDENT_REMAINING32_WHOLE_BOX_GEOMETRY_ONLY"
    assert all(leaf[k]==native_row[k] for k in ("pilot_id","record_id","radius_exact","box_sha256","body_geometry_sha256"))
    keys=("record_id","model_contract_sha256","point_cache_sha256","R_binary_sha256",
          "center_binary_sha256","exact_drive","energy_scale_binary_exact","derivative_policy")
    expected={k:tr[k] for k in keys}
    expected.update(pilot_id=case,radius_exact=str(a.R),box_sha256=tr["trial_box_sha256"],
                    archived_field_sha256=e["original_observer_field_sha256"],
                    archived_addition_gap_exact=e["original_addition_gap_exact"],
                    archived_observer_goal_exact=e["observer_goal_exact"],
                    archived_shift_exact=tr["archived_center_shift_exact"])
    assert all(leaf[k]==v for k,v in expected.items()) and leaf["replay_phase"]==geophase
    assert leaf["producer_case_results_file"]==a.relative(p/"results.json")
    assert leaf["producer_case_body_file"]==a.relative(p/"body_geometry.json")
    assert leaf["body_geometry_sha256"]==seals["body_geometry"]
    assert leaf["producer_case_results_sha256"]==seals["results"]
    assert leaf["producer_worker_terminal_sha256"]==seals["terminal"]
    assert leaf["producer_worker_observed_os_exit_sha256"]==seals["observed_os_exit"]
    assert leaf["all_144_body_and_24_trace_domains"] is True and leaf["root_asserted"] is False
    assert type(leaf["closed_full_nodal_coordinates_rebuilt"]) is int and leaf["closed_full_nodal_coordinates_rebuilt"]==650
    assert type(leaf["raw_surface_cache_geometry_lines"]) is int and leaf["raw_surface_cache_geometry_lines"]==288
    jlo,jhi=map(a.rational,leaf["body_J_bounds_exact"])
    assert Q(1,10**6)<jlo<=jhi and Q(0)<=a.rational(leaf["max_bulk_rho_outward_exact"])<1
    assert Q(0)<=a.rational(leaf["max_surface_rho_raw_exact"])<1
    assert all(a.rational(leaf[k])>0 for k in ("strong_monotonicity_lower_exact","min_surface_G_exact",
                                               "min_trace_hoop_exact","min_site_tangent_squared_exact"))
    body=a.read(p/"body_geometry.json",seals["body_geometry"])
    results=a.read(p/"results.json",seals["results"])
    assert results["trial"]==tr and results["whole_body_sha256"]==seals["body_geometry"]
    assert worker["candidate_status"]==results["status"]=="PASS_CANDIDATE_WHOLE_BOX_GEOMETRY_AWAITS_INDEPENDENT_REVIEW"
    assert body["status"]=="CERTIFIED_ENTIRE_NODAL_BOX_EMBEDDING"
    assert len(results["bulk_geometry_only"])==len(body["triangles"])==144 and len(results["trace_geometry_only"])==24
    assert all(len(x["cache_site_geometry"])==12 for x in results["trace_geometry_only"])
    return leaf


def direct_R(a,trial):
    contract=a.read(trial["model_contract"],trial["model_contract_sha256"])
    packet=a.read(contract["point_cache"],contract["point_cache_sha256"])
    blob=b64decode(packet["frozen_binary_cache"]["data"],validate=True)
    assert a.sha(blob)==contract["frozen_binary_cache_sha256"]==packet["frozen_binary_cache"]["sha256"]
    with np.load(io.BytesIO(blob),allow_pickle=False) as z:
        R=z["preconditioner"]
        u=z["u"]
    assert R.dtype.str==u.dtype.str=="<f8" and R.shape==(611,611) and u.shape==(611,)
    assert a.sha(R.tobytes(order="C"))==contract["binary_R_sha256"]==trial["R_binary_sha256"]
    assert a.sha(u.tobytes(order="C"))==contract["center_binary_sha256"]==trial["center_binary_sha256"]
    def row_sum(row):
        rs=[abs(v).as_integer_ratio() for v in row]
        denominator=max(d for _,d in rs)
        return Q(sum(n*(denominator//d) for n,d in rs),denominator)
    return max(map(row_sum,R.tolist()))


def run(observed_geometry):
    a=load_cold()
    cold=a.read(BR/"independent_prelaunch/review.json",COLD_SHA)
    assert cold["status"]=="PASS_INDEPENDENT_REMAINING32_EXACT_BRIDGE_COLD" and a.truth(cold["checks"])
    assert all(cold[k]==v for k,v in a.PHASE.items()) and cold["cold_review_source_sha256"]==COLD_SOURCE_SHA
    assert cold["selector_only_authorized_in_this_gate"] is True and cold["spatial_workers_authorized_in_this_gate"] is False
    p=a.read(a.P/"protocol.json",a.PHASE["protocol_sha256"])
    m=a.read(a.P/"manifest.json",a.PHASE["manifest_sha256"])
    a.read(a.SOURCE,a.PHASE["source_sha256"],False)
    a.read(a.P/"preflight.json",a.PHASE["preflight_sha256"])
    for path,h in cold["input_sha256"].items():
        a.read(path,h,False)
    gfolder=a.G/"independent_actual_preparation/actual_replay"
    draw=a.PINS[a.DR/"os_exit_observation.json"]
    geom,gt,go,gphase=replay_receipt(a,gfolder,observed_geometry,"geometry")
    der,dt,do,dphase=replay_receipt(a,a.DR,draw,"variation")
    rawg,gclose=dimensions(a,geom["checked_input_sha256"],a.CASES,RAW_GEOM_SHA,"geometry")
    vorder=[e["case_id"] for e in m["cases"] if e["cached_input_kind"]=="FRESH_R34"]
    rawd,dclose=dimensions(a,der["input_sha256"],vorder,a.PINS[a.D/"actual/root_terminal_admission.json"],"variation")
    assert [x["pilot_id"] for x in geom["geometry_records"]]==list(a.CASES)
    assert [x["case_id"] for x in der["case_results"]]==vorder
    assert geom["observed_producer_os"]["pid"]==rawg["supervisor_pid"] and geom["observed_producer_os"]["session"]==rawg["unified_exec_session"]
    assert der["root_observed_producer_pid"]==rawd["supervisor_pid"] and der["root_observed_producer_exec"]==rawd["unified_exec_session"]
    assert list(x["case_id"] for x in m["cases"])==list(a.CASES) and len(vorder)==10
    proof=[]
    for e,gr in zip(m["cases"],geom["geometry_records"],strict=True):
        case=e["case_id"]
        gl=checked_geo(a,e,gr,gphase,gclose,gfolder)
        cat=a.catalogue_match(e)
        trial=e["trial"]
        z,norm,signed=a.independently_read_derivative(e,der,rawd)
        if case in vorder:
            dr=a.single(der["case_results"],"case_id",case)
            leaf=a.read(dr["native_file"],dr["native_sha256"])
            assert leaf["phase"]==dphase and dr["producer_result_sha256"]==dclose[case][1]["results"]
            dh=dr["native_sha256"]
        else:
            dh=e["native_case_sha256"]
        assert direct_R(a,trial)==norm
        f,h,caps=a.exact_budget(cat,z,norm)
        assert len(signed)==611
        previous=a.single(cold["case_rows"],"case_id",case)
        assert previous["ZQ_upper_exact"]==str(z) and previous["R_inf_exact"]==str(norm)
        assert previous["per_component_force_and_hessian_tail_caps_exact"]==list(map(str,caps))
        proof.append({"case_id":case,"record_id":e["record_id"],"trial_box_sha256":trial["trial_box_sha256"],
                       "model_contract_sha256":trial["model_contract_sha256"],"cache_packet_sha256":trial["point_cache_sha256"],
                       "geometry_replay_case_sha256":gr["sha256"],"cached_derivative_case_sha256":dh,
                       "cached_input_kind":e["cached_input_kind"],"geometry_J_lower_exact":gl["body_J_bounds_exact"][0],
                       "same_signed_cache_anchor_rows":611,"R_inf_exact":str(norm),"reviewed_ZQ_upper_exact":str(z),
                       "localized_force_and_H_tail_caps_exact":list(map(str,caps))})
        print(case,flush=True)
    assert len(proof)==32
    out={"status":"PASS_INDEPENDENT_REMAINING32_SAME_TARGET_BRIDGE_INPUTS",
         **a.PHASE,"cold_review_sha256":COLD_SHA,
         "geometry_observed_exit_sha256":observed_geometry,"variation_observed_exit_sha256":draw,
         "geometry_replay_phase":GEOMETRY_PHASE,"geometry_receipt_sha256":go["receipt_sha256"],
         "variation_receipt_sha256":do["receipt_sha256"],"geometry_external_producer_admission_sha256":RAW_GEOM_SHA,
         "variation_external_producer_admission_sha256":a.PINS[a.D/"actual/root_terminal_admission.json"],
         "same_full_32_geometry_10_fresh_derivative_required":True,"spatial_workers_authorized":False,
         "checks":{"both_complete_new_native_replays_observed_OS0":True,
            "direct_canonical_raw_producer_supervisor_and_child_chains":True,
            "all_32_closed_body_and_trace_domain_receipts_same_trial":True,
            "same_22_native_and_10_fresh_derivatives_signed_centers_and_R":True,
            "all_32_matrix_infinity_norms_recomputed_from_exact_binary_bytes":True,
            "all_168_routes_per_case_and_strict_tail_sizes_match_independent_cold_review":True},
         "cases":proof,"source_and_input_sha256":a.LINKS,"join_source_sha256":a.sha(Path(__file__).read_bytes()),
         "no_new_selector_geometry_derivative_integral_or_solve":True,
         "scope":"The joined proof supplies fixed-cube inputs for separately observed and reviewed selectors. No numerical component, exact-spatial root or minimum is asserted."}
    path=HERE/"bridge_input_admission.json"
    with path.open("x",encoding="utf-8",newline="\n") as f:
        json.dump(out,f,indent=2)
        f.write("\n")
    print(json.dumps({"joined_sha256":a.sha(path.read_bytes()),"files_rehashed":len(a.LINKS)}))


if __name__=="__main__":
    cli=argparse.ArgumentParser()
    cli.add_argument("--geometry-observed-os-sha256",required=True)
    args=cli.parse_args()
    if not __debug__:
        raise RuntimeError("exact premise assertions required")
    run(args.geometry_observed_os_sha256)

