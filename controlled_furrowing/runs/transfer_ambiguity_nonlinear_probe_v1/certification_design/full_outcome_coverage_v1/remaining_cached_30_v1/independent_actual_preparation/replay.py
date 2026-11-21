"""One independent cache replay for eight separately observed cached wave terminals."""
import argparse
from fractions import Fraction as Q
from hashlib import sha256
import importlib.metadata
import importlib.util
import json
import os
from pathlib import Path
import re
import sys
import time
import traceback

SOURCE=Path(__file__).resolve()
HERE=SOURCE.parent
WAVESET=HERE.parent
DESIGN=WAVESET.parents[1]
ROOT=DESIGN.parents[2]
PRE=WAVESET/"preparation"
ACTUAL=WAVESET/"actual"
REPLAY=HERE/"actual_replay"
PRIOR=DESIGN/"cached_pilots_v1/independent_actual_preparation/replay.py"
AD=ROOT/"scripts/review_algebraic_root.py"
DRIVER=ROOT/"scripts/transfer_ambiguity_remaining30_cached.py"
GATES=PRE/"independent_prelaunch"
PH=PRE/"protocol.json"
PF=PRE/"preflight.json"
WI=PRE/"wave_index.json"
RI=GATES/"review_index.json"
CHILD=("child_started","worker_started","observed_os_exit","terminal","results")
COUNTS={"bulk":9216,"surface":288}
N=611
LADDER=tuple(Q(1,2**k) for k in (20,24,28,32,34,36))
PINS={
    DRIVER:"23c6e4033bcb9f718a31b75a1c6238afe611fab0f99bf5cd5b16c17298abe702",
    PH:"e31870b0c6ac61e95a7b2657eab73d8231657af86c65f7a1d1172ff0b45a1e4c",
    PF:"1cad3ad4f7a3b8d486b16d511b0748d215830cab2fb6e74e16143bd3f066a11d",
    WI:"6aa242a734bfd6a2524aa8fae6f899285a0ba35305731abb16961f3b7c4d057c",
    RI:"d91b2e899813d29a0d2816d487d5adccb1efd718a386ebe8ca298e1793704ab1",
    GATES/"support.json":"79dddd0eed61dcf67e3264e3b75b7047f6a3d6aa4f4ceabe7f5e8b63c44fb4c8",
    PRIOR:"3f71e5cae5ee314e260e3bef6c478a923fb7c28f446c184c53d912d6ae296a94",
    AD:"6af2311d0667a28af263b6fb6237063481aa8e65209ae9722c2b0246a8cff27a",
}
SCOPE=("Independent local root tests for the prescribed exact-real binary Q8/S12 cached potentials. "
       "No exact spatial integral, stability, whole-body geometry, continuum, observer transfer or "
       "experimental test. A failed sufficient test is unresolved.")


def name(p):
    return p.relative_to(ROOT).as_posix()


def sha(p):
    return sha256(p.read_bytes()).hexdigest()


def read(p,d):
    assert type(d) is str and re.fullmatch("[0-9a-f]{64}",d)
    raw=p.read_bytes()
    assert sha256(raw).hexdigest()==d,p
    return json.loads(raw)


def save(p,out):
    with p.open("x",encoding="utf-8",newline="\n") as h:
        json.dump(out,h,indent=2,allow_nan=False)
        h.write("\n")
    return sha(p)


def yes(d):
    return type(d) is dict and bool(d) and all(v is True for v in d.values())


def authority():
    assert all(sha(f)==h for f,h in PINS.items())
    p,producer_index,review_index=[read(f,PINS[f]) for f in (PH,WI,RI)]
    assert review_index["status"]=="PASS_ALL_EIGHT_DISJOINT_WAVE_GATES"
    assert producer_index["status"]=="FROZEN_EIGHT_DISJOINT_CACHED_WAVES"
    assert review_index["support_sha256"]==PINS[GATES/"support.json"]
    assert tuple(map(Q,p["finite_radius_ladder_exact"]))==LADDER
    assert p["site_counts"]==COUNTS and p["precision_bits"]==256
    assert p["new_equilibrium_solves"]==0 and p["prior_case_cache_or_R_reused"] is False
    hashes={name(f):h for f,h in PINS.items()}|p["source_closure_sha256"]|p["input_sha256"]
    waves={}
    for entry in producer_index["waves"]:
        w=entry["wave_id"]
        g=review_index["waves"][w]
        gate_file,gate_sha=g["gate_file"],g["sha256"]
        assert gate_file==name(GATES/f"{w}_review.json")
        assert entry["file"]==f"{w}.json" and entry["sha256"]==g["manifest_sha256"]
        mpath=PRE/entry["file"]
        wave=read(mpath,entry["sha256"])
        gate=read(ROOT/gate_file,gate_sha)
        ids=wave["ordered_case_ids"]
        assert ids==p["waves"][w]==entry["ordered_case_ids"]==gate["ordered_case_ids"]
        assert gate["status"]=="PASS_INDEPENDENT_REMAINING30_CACHED_WAVE_PRELAUNCH" and yes(gate["checks"])
        assert gate["wave_id"]==w==wave["wave_id"] and gate["wave_manifest_sha256"]==entry["sha256"]
        assert wave["source_sha256"]==gate["source_sha256"]==PINS[DRIVER]
        assert wave["protocol_sha256"]==gate["protocol_sha256"]==PINS[PH]
        assert gate["preflight_sha256"]==PINS[PF] and gate["wave_index_sha256"]==PINS[WI]
        assert gate["shared_review_sha256"]==PINS[GATES/"support.json"]
        assert wave["producer_actual_work_in_preparation"]==0
        paths=[name(ACTUAL/w/f"{k}.json") for k in ("started","results","terminal")]
        paths+=[name(ACTUAL/w/ident/f"{k}.json") for ident in ids for k in CHILD]
        assert len(paths)==len(set(paths))==3+5*len(ids)
        hashes[name(mpath)]=entry["sha256"]
        hashes[gate_file]=gate_sha
        waves[w]={"ordered_case_ids":ids,"producer_wave_manifest_sha256":entry["sha256"],
                  "producer_wave_gate_file":gate_file,"producer_wave_gate_sha256":gate_sha,
                  "required_root_terminal_pins":paths}
    assert list(waves)==list(p["waves"])==list(review_index["waves"])
    assert sorted(ident for w in waves.values() for ident in w["ordered_case_ids"])==sorted(p["cases"])
    assert len(p["cases"])==30 and len(waves)==8
    return p,hashes,waves


def protocol():
    _,hashes,waves=authority()
    return {"status":"FROZEN_ALL_EIGHT_INDEPENDENT_CACHED_WAVE_REPLAYS","source_sha256":sha(SOURCE),
            "immutable_source_and_input_hashes":hashes,"waves":waves,"site_counts":COUNTS,"precision_bits":256,
            "radius_ladder_exact":list(map(str,LADDER)),"scope":SCOPE,
            "reuse_limits":"Reuse independent interval AD and binary cache regeneration primitives only. "
               "No prior pilot target, point cache, R, convergence claim or physical input authority is loaded. "
               "Fresh own whole-cube variation for each actually evaluated radius."}


def external_shape(d,plan,wave):
    assert d["status"]=="PASS_ROOT_OBSERVED_REMAINING30_CACHED_WAVE_TERMINALS"
    assert d["wave_id"]==wave
    assert type(d["unified_exec_session"]) is int and d["unified_exec_session"]>0
    assert type(d["supervisor_pid"]) is int and d["supervisor_pid"]>0
    assert type(d["observed_exit_code"]) is int and d["observed_exit_code"]==0
    assert set(d["sha256"])==set(plan["waves"][wave]["required_root_terminal_pins"])
    assert all(type(h) is str and re.fullmatch("[0-9a-f]{64}",h) for h in d["sha256"].values())
    return d["sha256"]


def radius_candidates(attempt_ids,primary,backup):
    if primary is None:
        # If the producer retained no strict pass, the final evaluated radius
        # can still receive one independent sufficient test.
        assert backup is None
        return [attempt_ids[-1]] if attempt_ids else []
    assert primary in attempt_ids
    if backup is None:
        assert primary==attempt_ids[-1]
        return [primary]
    assert attempt_ids[-2:]==[primary,backup] and backup<primary
    assert LADDER.index(backup)==LADDER.index(primary)+1
    return [primary,backup]


def prepare():
    plan=protocol()
    assert not REPLAY.exists()
    for w in plan["waves"]:
        d={"status":"PASS_ROOT_OBSERVED_REMAINING30_CACHED_WAVE_TERMINALS","wave_id":w,
           "unified_exec_session":5,"supervisor_pid":6,"observed_exit_code":0,
           "sha256":dict.fromkeys(plan["waves"][w]["required_root_terminal_pins"],"a"*64)}
        assert external_shape(d,plan,w)==d["sha256"]
        for bad in (d|{"observed_exit_code":False},d|{"sha256":{}},d|{"wave_id":"invalid"}):
            try:
                external_shape(bad,plan,w)
            except AssertionError:
                continue
            raise AssertionError("bad OS gate accepted")
    assert radius_candidates(list(LADDER[:5]),LADDER[3],LADDER[4])==list(LADDER[3:5])
    assert radius_candidates(list(LADDER),None,None)==[LADDER[-1]]
    assert radius_candidates([],None,None)==[]
    h=save(HERE/"protocol.json",plan)
    f=save(HERE/"preflight.json",{"status":"PASS","source_sha256":sha(SOURCE),"protocol_sha256":h,
        "synthetic_all_eight_terminal_shapes_and_radius_handoff":True,
        "actual_producer_payloads_opened":0,"actual_site_AD":0,"actual_models_or_preconditioners":0})
    print(h,f,sha(SOURCE),flush=True)


def admit(plan,wave,path,digest):
    a=read(path,digest)
    pins=external_shape(a,plan,wave)
    docs={f:read(ROOT/f,h) for f,h in pins.items()}
    get=lambda p:docs[name(p)]
    root=ACTUAL/wave
    started,end,all_results=[get(root/f"{k}.json") for k in ("started","terminal","results")]
    p,_,waves=authority()
    won=waves[wave]
    gate=read(ROOT/won["producer_wave_gate_file"],won["producer_wave_gate_sha256"])
    bind={k:gate[k] for k in ("source_sha256","protocol_sha256","preflight_sha256","wave_index_sha256",
        "wave_manifest_sha256","coverage_ledger_sha256","verified_two_M_coverage_sha256","first_four_cached_coverage_sha256")}
    bind["independent_wave_review_sha256"]=won["producer_wave_gate_sha256"]
    ids=won["ordered_case_ids"]
    assert started["pid"]==end["pid"]==a["supervisor_pid"] and end["status"]=="FINISHED"
    assert end["results_sha256"]==pins[name(root/"results.json")]
    assert all(r["wave_id"]==wave and r["bindings"]==bind for r in (started,end,all_results))
    assert started["ordered_case_ids"]==end["ordered_case_ids"]==ids
    assert all_results["status"]=="COMPLETE_CACHED_WAVE_ATTEMPTS_UNREVIEWED"
    assert all_results["independent_actual_review_required"] is True and all_results["certified_cached_roots_by_this_supervisor"]==0
    assert all_results["wave_manifest_sha256"]==won["producer_wave_manifest_sha256"]
    assert [v["case_id"] for v in end["workers"]]==[v["case_id"] for v in all_results["cases"]]==ids
    assert len({v["pid"] for v in end["workers"]})==len(ids)
    cases=[]
    for ident,r,w in zip(ids,all_results["cases"],end["workers"],strict=True):
        x=root/ident
        cs,ws,o,t,result=[get(x/f"{k}.json") for k in CHILD]
        assert w==o|{"os_exit_receipt_sha256":pins[name(x/"observed_os_exit.json")]}
        assert type(o["observed_exit_code"]) is int and o["observed_exit_code"]==0
        assert o["terminal_sha256"]==pins[name(x/"terminal.json")]
        assert t["status"]=="FINISHED" and t["results_sha256"]==r["results_sha256"]==pins[name(x/"results.json")]
        assert ws=={"pid":o["pid"],"supervisor_pid":a["supervisor_pid"],"wave_id":wave,"case_id":ident,"bindings":bind}
        assert all(all(doc[k]==v for k,v in ws.items()) for doc in (cs,o,t))
        assert result["case"]==p["cases"][ident] and result["bindings"]==bind and result["status"]==r["status"]
        assert result["status"] in ("CANDIDATE_CACHED_ROOT_AWAITING_INDEPENDENT_REVIEW","UNRESOLVED_CACHED_CONTRACTION_OR_BUDGET")
        assert result["cache_packet"]==r["cache_packet"]==name(x/"cache_packet.json")
        assert result["cache_packet_sha256"]==r["cache_packet_sha256"]
        assert result["independent_actual_review_required"] is True and result["new_nonlinear_equilibrium_solves"]==0
        assert result["exact_spatial_integration_established"] is result["whole_body_geometry_established"] is result["local_minimum_established"] is False
        cases.append((x,result))
    return a,p,pins,cases


def verify_radii(result,case,d,links):
    rows=result["radius_attempts"]
    goal=Q(case["e_goal_dimensionless_exact"])
    assert 0<len(rows)<=len(LADDER)
    primary=result["accepted_radius"]
    backup=result["narrower_predeclared_replay_backup"]
    evaluated=[]
    for k,(row,r) in enumerate(zip(rows,LADDER,strict=False)):
        assert Q(row["radius_requested_exact"])==r
        passed=row["strict_selfmap_contraction_and_observer_radius"]
        assert type(passed) is bool
        if r>=goal:
            assert row["status"]=="SKIPPED_OUTSIDE_OBSERVER_GOAL" and passed is False and "path" not in row
            continue
        path=d/f"radius_{k:02d}.json"
        evaluated.append(r)
        assert row["path"]==path.name
        links[name(path)]=row["sha256"]
        stored=read(path,row["sha256"])
        assert stored=={k:v for k,v in row.items() if k not in ("path","sha256")}|{
              "case_id":case["case_id"],"cache_packet_sha256":result["cache_packet_sha256"]}
        if row["status"]=="UNRESOLVED_QUADRATURE_ANALYTIC_DOMAIN":
            assert passed is False
            continue
        eta,z=map(Q,(row["eta_upper_exact"],row["z_upper_exact"]))
        assert 0<=eta==Q(result["eta_upper_exact"]) and 0<=z
        assert z==Q(result["center_Z0_upper_exact"])+Q(row["preconditioned_variation_upper_exact"])
        assert Q(row["eta_plus_zr_upper_exact"])==eta+z*r
        assert Q(row["strict_selfmap_gap_lower_exact"])==r-(eta+z*r)
        expected=eta/(1-z) if z<1 else None
        assert row["posterior_error_upper_exact"]==(str(expected) if expected is not None else None)
        assert passed is (z<1 and eta+z*r<r)
    primary_record=next((row for row in rows if row["strict_selfmap_contraction_and_observer_radius"]),None)
    assert primary_record==primary
    if primary is None:
        assert len(rows)==len(LADDER) and result["status"]=="UNRESOLVED_CACHED_CONTRACTION_OR_BUDGET"
    else:
        assert result["status"]=="CANDIDATE_CACHED_ROOT_AWAITING_INDEPENDENT_REVIEW"
        first=rows.index(primary)
        if first+1<len(LADDER):
            assert len(rows)==first+2 and backup==rows[first+1]
        else:
            assert backup is None and len(rows)==first+1
    return radius_candidates(evaluated,Q(primary["radius_requested_exact"]) if primary else None,
                              Q(backup["radius_requested_exact"]) if backup else None)


def verified_case(helper,ad,p,d,result,links):
    t0=time.perf_counter()
    c=result["case"]
    ident=c["case_id"]
    radii=verify_radii(result,c,d,links)
    cache_path=d/"cache_packet.json"
    center=d/"point_enclosure.json"
    links[name(cache_path)]=result["cache_packet_sha256"]
    links[name(center)]=result["point_enclosure_sha256"]
    packet=read(cache_path,links[name(cache_path)])
    point=read(center,links[name(center)])
    assert packet["case_id"]==point["case_id"]==ident
    assert packet["bindings"]==point["bindings"]==result["bindings"]
    assert packet["status"]=="FROZEN_FRESH_EXACT_BINARY_QUADRATURE_CACHE"
    assert packet["closed_parameter_cell_exact"]==c["closed_parameter_cell_exact"]
    assert packet["parameter_order"]==c["parameter_order"]
    assert point["expected_and_rebuilt_site_counts"]==COUNTS and point["cache_packet_sha256"]==result["cache_packet_sha256"]
    assert packet["frozen_binary_cache"]["sha256"]==result["embedded_binary_cache_sha256"]
    assert Q(point["eta_upper_exact"])==Q(result["eta_upper_exact"])
    assert Q(point["center_Z0_upper_exact"])==Q(result["center_Z0_upper_exact"])
    assert all(c["drive"][k]==packet[k] for k in c["drive"])
    # Spell the case identity for the explicit old mathematical helper only.
    explicit=c|{"pilot_id":ident,"exact_drive":c["drive"]}
    arrays,strength,_=helper.cache_extract(p,packet|{"pilot_id":ident},explicit,d,links)
    want=Q(c["material"]["gain_exact"])*Q(c["step"]["dose_exact"])
    assert strength==want==Q(c["drive"]["activation_strength_exact"])
    assert strength!=Q(c["drive"]["activation_strength_binary_exact"])
    ad.A.cache_clear()
    try:
        g,H,cache,counts=ad.assemble(ad.terms_from_arrays(arrays,packet["parameters"],strength),N)
    except ValueError as e:
        return {"case_id":ident,"record_id":c["record_id"],"scope":SCOPE,
                "status":"UNRESOLVED_INDEPENDENT_CACHED_POINT_DOMAIN","detail":str(e)}
    assert counts==COUNTS and len(cache)==sum(COUNTS.values())
    own_attempts=[]
    chosen=None
    goal=Q(c["e_goal_dimensionless_exact"])
    for attempt,r in enumerate(radii):
        assert 0<r<goal
        try:
            rg,eta,z0,dz,Rnorm,domain=helper.replay_bound(ad,g,H,cache,arrays["preconditioner"],r,COUNTS)
        except ValueError as e:
            own_attempts.append({"status":"UNRESOLVED_INDEPENDENT_WHOLE_CUBE_DOMAIN","radius_exact":str(r),
                                 "reason":str(e),"new_own_variation_evaluation_started":True})
            continue
        assert domain["counts"]==COUNTS and set(domain["strict_lower_margins"])==helper.MARGINS
        assert all(Q(v["lower_exact"])>0 for v in domain["strict_lower_margins"].values())
        assert Q(domain["strict_lower_margins"]["bulk_J_min"]["lower_exact"])>Q(1,10**6)
        assert Rnorm<=Q(point["preconditioner_infinity_norm_upper_exact"])
        signed=[ad.pack(rg[i,0]) for i in range(N)]
        for v,w in zip(signed,point["preconditioned_center_gradient_by_row"],strict=True):
            assert Q(v["lower_exact"])<=Q(w["upper_exact"]) and Q(w["lower_exact"])<=Q(v["upper_exact"])
        z=z0+dz
        lhs=eta+z*r
        ok=0<=z<1 and lhs<r
        post=eta/(1-z) if ok else None
        if ok:
            assert post<r<goal
        row={"status":"PASS_INDEPENDENT_CACHED_LOCAL_ROOT" if ok else "UNRESOLVED_INDEPENDENT_CACHED_CONTRACTION",
            "radius_exact":str(r),"is_independent_backup_attempt":attempt>0,
            "whole_cube_AD_recomputed_for_this_radius":True,"analytic_domain":domain,
            "eta_upper_exact":str(eta),"center_Z0_upper_exact":str(z0),"variation_upper_exact":str(dz),
            "ZQ_upper_exact":str(z),"strict_selfmap_lhs_upper_exact":str(lhs),
            "preconditioner_infinity_norm_exact":str(Rnorm),
            "signed_center_force_intervals":signed,
            "posterior_error_upper_exact":str(post) if ok else None,
            "display":{"z":float(z),"eta":float(eta),"lhs":float(lhs),
                       "error":float(post) if ok else None}}
        own_attempts.append(row)
        if ok:
            chosen=row
            break
    return {"status":"PASS_INDEPENDENT_CACHED_LOCAL_ROOT" if chosen else "UNRESOLVED_INDEPENDENT_CACHED_ROOT",
            "case_id":ident,"record_id":c["record_id"],"field":c["field"],"field_sha256":c["field_sha256"],
            "scope":SCOPE,"site_counts":COUNTS,"exact_requested_strength":str(strength),
            "observer_radius_budget_exact":str(goal),"rounded_strength_used_in_proof":False,
            "cache_packet_sha256":result["cache_packet_sha256"],"embedded_cache_sha256":packet["frozen_binary_cache"]["sha256"],
            "per_tensor_sha256":{k:sha256(x.tobytes(order="C")).hexdigest() for k,x in arrays.items()},
            "approved_radii_for_independent_replay_exact":list(map(str,radii)),
            "own_radius_attempts":own_attempts,"chosen_radius_exact":chosen["radius_exact"] if chosen else None,
            "posterior_error_upper_exact":chosen["posterior_error_upper_exact"] if chosen else None,
            "new_R_or_solver_evaluations":0,"seconds":time.perf_counter()-t0}


def run(args):
    ph,fh=args.protocol_sha256,args.preflight_sha256
    plan=read(HERE/"protocol.json",ph)
    pf=read(HERE/"preflight.json",fh)
    assert plan==protocol() and pf["status"]=="PASS" and pf["source_sha256"]==sha(SOURCE)
    assert pf["protocol_sha256"]==ph and pf["actual_producer_payloads_opened"]==pf["actual_site_AD"]==0
    for f,h in plan["immutable_source_and_input_hashes"].items():
        assert sha(ROOT/f)==h,f
    a,p,links,case_inputs=admit(plan,args.wave_id,args.admission.resolve(),args.admission_sha256)
    spec=importlib.util.spec_from_file_location("prior_independent_cache_primitives",PRIOR)
    helper=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(helper)
    assert helper.ROOT==ROOT and helper.PINS[helper.AD_SOURCE]==PINS[AD]
    ad=helper.load_ad()
    import torch
    torch.set_num_threads(1)
    assert torch.get_num_threads()==1 and ad.ctx.prec==256
    assert all(importlib.metadata.version(k)==v for k,v in p["runtime_packages"].items())
    folder=REPLAY/args.wave_id
    folder.mkdir(parents=True,exist_ok=False)
    tag={"pid":os.getpid(),"wave_id":args.wave_id,"source_sha256":sha(SOURCE),
         "protocol_sha256":ph,"preflight_sha256":fh,"root_terminal_admission_sha256":args.admission_sha256}
    save(folder/"started.json",tag|{"status":"RUNNING"})
    t=tag|{"status":"FAILED","exit_code":1}
    t0=time.perf_counter()
    try:
        checked=[]
        for d,r in case_inputs:
            result=verified_case(helper,ad,p,d,r,links)
            f=folder/f"{result['case_id']}.json"
            h=save(f,result)
            row={k:result.get(k) for k in ("status","case_id","record_id","field_sha256","cache_packet_sha256",
                 "chosen_radius_exact","posterior_error_upper_exact")}
            checked.append(row|{"native_file":name(f),"native_sha256":h})
            print(json.dumps(checked[-1]),flush=True)
        assert [r["case_id"] for r in checked]==plan["waves"][args.wave_id]["ordered_case_ids"]
        assert all(sha(ROOT/f)==h for f,h in links.items())
        n=sum(r["status"]=="PASS_INDEPENDENT_CACHED_LOCAL_ROOT" for r in checked)
        sha_r=save(folder/"receipt.json",{"status":"PASS_INDEPENDENT_CACHED_WAVE_LOCAL_ROOTS" if n==len(checked) else
                  "INDEPENDENT_CACHED_WAVE_REPLAY_WITH_UNRESOLVED_CASES","scope":SCOPE,**tag,
                  "producer_root_observed_exit":{k:a[k] for k in ("supervisor_pid","unified_exec_session","observed_exit_code")},
                  "producer_and_payload_sha256":links,
                  "producer_wave_manifest_sha256":plan["waves"][args.wave_id]["producer_wave_manifest_sha256"],
                  "producer_wave_gate_sha256":plan["waves"][args.wave_id]["producer_wave_gate_sha256"],
                  "checks":{"all_wave_cases_read_at_same_exact_fixed_inputs":True,
                            "accepted_radii_independently_enclosed":True,
                            "all_wave_cases_passed_independent_cached_root_test":n==len(checked)},
                  "cached_local_root_count":n,"exact_spatial_root_count":0,"new_equilibrium_solves_or_preconditioners":0,
                  "cases":checked,"seconds":time.perf_counter()-t0})
        t|={"status":"FINISHED","exit_code":0,"receipt_sha256":sha_r}
    except BaseException as error:
        save(folder/"failure.json",{"message":str(error),"traceback":traceback.format_exc()})
        raise
    finally:
        save(folder/"terminal.json",t)
        print(json.dumps(t),flush=True)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    sub=parser.add_subparsers(dest="stage",required=True)
    sub.add_parser("prepare")
    x=sub.add_parser("run")
    x.add_argument("--wave-id",required=True)
    x.add_argument("--admission",required=True,type=Path)
    x.add_argument("--admission-sha256",required=True)
    x.add_argument("--protocol-sha256",required=True)
    x.add_argument("--preflight-sha256",required=True)
    args=parser.parse_args()
    prepare() if args.stage=="prepare" else run(args)
