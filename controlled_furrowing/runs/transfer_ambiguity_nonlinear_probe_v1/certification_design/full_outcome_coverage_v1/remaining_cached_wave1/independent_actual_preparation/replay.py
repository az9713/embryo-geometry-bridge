"""Independent replay of the two true-terminal fresh M cached potentials."""
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
WAVE=HERE.parent
DESIGN=WAVE.parents[1]
ROOT=DESIGN.parents[2]
PRE=WAVE/"preparation"
ACTUAL=WAVE/"actual"
REPLAY=HERE/"actual_replay"
PRIOR=DESIGN/"cached_pilots_v1/independent_actual_preparation/replay.py"
AD=ROOT/"scripts/review_algebraic_root.py"
DRIVER=ROOT/"scripts/transfer_ambiguity_remaining_cached.py"
IDS=("max_contrast_change_small_09","max_negative_contrast_change_small_09")
CHILD=("child_started","worker_started","observed_os_exit","terminal","results")
PINS={
    DRIVER:"a0ec457505c6a2e82aab91bd985d6f4383864124df4a1e295e3119e1f8694480",
    PRE/"protocol.json":"c1543cbbd2520664313b56114b66a885faaf4d75869f56dbb1aa83a43a78fe3a",
    PRE/"preflight.json":"cf328f574402555363ee4f74e2ca1afbb9cf507ee871d1a6501ca1bace9eb0dc",
    PRE/"independent_prelaunch/review.json":"3be5762af861524a7084e2a09a38389d69130cf987fb0dd530618a0e4f7981b5",
    PRIOR:"3f71e5cae5ee314e260e3bef6c478a923fb7c28f446c184c53d912d6ae296a94",
    AD:"6af2311d0667a28af263b6fb6237063481aa8e65209ae9722c2b0246a8cff27a",
}
COUNTS={"bulk":9216,"surface":288}
N=611
SCOPE=("Independent local roots of two exact-real binary Q8/S12 quadrature potentials at the "
       "exact prescribed M drives. This supplies no exact spatial integration, whole-body "
       "geometry, stability, continuum, source or experimental compatibility claim.")


def name(p):
    return p.relative_to(ROOT).as_posix()


def sha(p):
    return sha256(p.read_bytes()).hexdigest()


def pinned(p,d):
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


def finite_plan():
    assert all(sha(f)==h for f,h in PINS.items())
    prot=pinned(PRE/"protocol.json",PINS[PRE/"protocol.json"])
    gate=pinned(PRE/"independent_prelaunch/review.json",PINS[PRE/"independent_prelaunch/review.json"])
    assert gate["status"]=="PASS_INDEPENDENT_REMAINING_M_CACHED_PRELAUNCH" and yes(gate["checks"])
    assert tuple(prot["ordered_case_ids"])==tuple(gate["ordered_case_ids"])==IDS
    assert prot["site_counts"]==COUNTS and prot["precision_bits"]==256
    assert prot["new_nonlinear_equilibrium_solves"]==0 and prot["old_pilot_case_reused"] is False
    paths=[name(ACTUAL/f"{k}.json") for k in ("started","results","terminal")]
    paths+=[name(ACTUAL/p/f"{k}.json") for p in IDS for k in CHILD]
    assert len(paths)==len(set(paths))==13
    return {"status":"FROZEN_TWO_M_INDEPENDENT_CACHED_REPLAY","source_sha256":sha(SOURCE),
            "input_hashes":{name(f):h for f,h in PINS.items()}|prot["source_closure_sha256"]|prot["input_sha256"],
            "expected_external_paths":paths,"expected_producer_session":22696,
            "precision_bits":256,"site_counts":COUNTS,"ordered_case_ids":list(IDS),
            "scope":SCOPE,"reuse":"Independent interval AD derivatives, binary model regeneration, "
              "and matrix kernels from the prior independently reviewed source. No prior pilot target, "
              "preconditioner, field, cache or proof input is reused."}


def external_shape(d,p):
    assert d["status"]=="PASS_ROOT_OBSERVED_REMAINING_M_CACHED_TERMINALS"
    assert d["unified_exec_session"]==p["expected_producer_session"]
    assert type(d["supervisor_pid"]) is int and d["supervisor_pid"]>0
    assert type(d["observed_exit_code"]) is int and d["observed_exit_code"]==0
    assert set(d["sha256"])==set(p["expected_external_paths"])
    assert all(type(h) is str and re.fullmatch("[0-9a-f]{64}",h) for h in d["sha256"].values())
    return d["sha256"]


def prepare():
    p=finite_plan()
    doc={"status":"PASS_ROOT_OBSERVED_REMAINING_M_CACHED_TERMINALS","unified_exec_session":22696,
         "supervisor_pid":1,"observed_exit_code":0,"sha256":dict.fromkeys(p["expected_external_paths"],"a"*64)}
    assert external_shape(doc,p)==doc["sha256"]
    for wrong in (doc|{"observed_exit_code":False},doc|{"sha256":{}},doc|{"unified_exec_session":-1}):
        try:
            external_shape(wrong,p)
        except AssertionError:
            continue
        raise AssertionError("failed or partial producer admitted")
    q=Q(1,10)
    assert q!=Q(float(q))
    r,e,z=Q(1,2**36),Q(1,2**40),Q(1,3)
    assert e+z*r<r and e/(1-z)<r
    h=save(HERE/"protocol.json",p)
    f=save(HERE/"preflight.json",{"status":"PASS","source_sha256":sha(SOURCE),"protocol_sha256":h,
         "local_fixtures_check_exact_fraction_radius_and_exit":True,"actual_producer_payloads_opened":0,
         "actual_site_ad_replays":0,"new_point_state_or_model_built":0})
    print(h,f,sha(SOURCE),flush=True)


def admit(p,file,digest):
    a=pinned(file,digest)
    links=external_shape(a,p)
    docs={path:pinned(ROOT/path,h) for path,h in links.items()}
    get=lambda f:docs[name(f)]
    start,end,all_results=[get(ACTUAL/f"{k}.json") for k in ("started","terminal","results")]
    prot=pinned(PRE/"protocol.json",PINS[PRE/"protocol.json"])
    gate=pinned(PRE/"independent_prelaunch/review.json",PINS[PRE/"independent_prelaunch/review.json"])
    b={k:gate[k] for k in ("source_sha256","protocol_sha256","preflight_sha256","coverage_ledger_sha256",
                           "whole_44_field_manifest_sha256","budget_independent_actual_review_sha256")}
    b["independent_prelaunch_review_sha256"]=PINS[PRE/"independent_prelaunch/review.json"]
    assert end["status"]=="FINISHED" and start["pid"]==end["pid"]==a["supervisor_pid"]
    assert end["results_sha256"]==links[name(ACTUAL/"results.json")]
    assert start["bindings"]==end["bindings"]==all_results["bindings"]==b
    assert all_results["status"]=="COMPLETE_TWO_M_CACHED_ATTEMPTS_UNREVIEWED"
    assert all_results["independent_actual_review_required"] is True
    assert tuple(start["ordered_case_ids"])==IDS
    assert tuple(x["case_id"] for x in end["workers"])==tuple(x["case_id"] for x in all_results["cases"])==IDS
    assert len({x["pid"] for x in end["workers"]})==2
    output=[]
    for c,r,w in zip(prot["cases"],all_results["cases"],end["workers"],strict=True):
        ident=c["case_id"]
        d=ACTUAL/ident
        assert ident==r["case_id"]==w["case_id"]
        cs,ws,exit,t,result=[get(d/f"{k}.json") for k in CHILD]
        assert w==exit|{"os_exit_receipt_sha256":links[name(d/"observed_os_exit.json")]}
        assert type(exit["observed_exit_code"]) is int and exit["observed_exit_code"]==0
        assert exit["terminal_sha256"]==links[name(d/"terminal.json")]
        assert t["status"]=="FINISHED"
        assert t["results_sha256"]==r["results_sha256"]==links[name(d/"results.json")]
        assert ws=={"pid":exit["pid"],"supervisor_pid":a["supervisor_pid"],"case_id":ident,"bindings":b}
        assert all(all(x[k]==v for k,v in ws.items()) for x in (cs,exit,t))
        assert result["case"]==c and result["bindings"]==b and result["status"]==r["status"]
        assert result["status"] in ("CANDIDATE_CACHED_ROOT_AWAITING_INDEPENDENT_REVIEW",
                                    "UNRESOLVED_CACHED_CONTRACTION_OR_BUDGET")
        assert result["cache_packet"]==r["cache_packet"]==name(d/"cache_packet.json")
        assert result["cache_packet_sha256"]==r["cache_packet_sha256"]
        assert result["independent_actual_review_required"] is True
        assert result["new_nonlinear_equilibrium_solves"]==0
        assert result["exact_spatial_integration_established"] is result["whole_body_geometry_established"] is result["local_minimum_established"] is False
        output.append((d,result))
    return a,prot,links,output


def attempts(p,root,case,result,links):
    rows=result["radius_attempts"]
    radii=list(map(Q,p["finite_radius_ladder_exact"]))
    assert 0<len(rows)<=len(radii)
    assert all(0<r<Q(case["e_goal_dimensionless_exact"]) for r in radii)
    pass_last=result["accepted_radius"] is not None
    for k,row in enumerate(rows):
        f=root/f"radius_{k:02d}.json"
        assert row["path"]==f.name and Q(row["radius_requested_exact"])==radii[k]
        links[name(f)]=row["sha256"]
        exact=pinned(f,row["sha256"])
        assert exact=={a:b for a,b in row.items() if a not in ("path","sha256")}|{
             "case_id":case["case_id"],"cache_packet_sha256":result["cache_packet_sha256"]}
        assert row["strict_selfmap_contraction_and_observer_radius"] is (pass_last and k==len(rows)-1)
        if row["status"]=="UNRESOLVED_QUADRATURE_ANALYTIC_DOMAIN":
            continue
        e,z=Q(row["eta_upper_exact"]),Q(row["z_upper_exact"])
        assert e==Q(result["eta_upper_exact"])
        assert z==Q(result["center_Z0_upper_exact"])+Q(row["preconditioned_variation_upper_exact"])
        assert Q(row["eta_plus_zr_upper_exact"])==e+z*radii[k]
        assert Q(row["strict_selfmap_gap_lower_exact"])==radii[k]-e-z*radii[k]
        if z<1:
            assert Q(row["posterior_error_upper_exact"])==e/(1-z)
        else:
            assert row["posterior_error_upper_exact"] is None
    if pass_last:
        assert result["status"]=="CANDIDATE_CACHED_ROOT_AWAITING_INDEPENDENT_REVIEW"
        assert result["accepted_radius"]=={a:b for a,b in rows[-1].items() if a not in ("path","sha256")}
        assert rows[-1]["status"]=="STRICT_CACHED_QUADRATURE_ROOT_CANDIDATE"
    else:
        assert len(rows)==len(radii)
    return radii[len(rows)-1],pass_last


def cached_case(helper,ad,prot,folder,result,links):
    t0=time.perf_counter()
    c=result["case"]
    ident=c["case_id"]
    assert c["record_id"]=="optoGEF2_small_col_M"
    assert c["checkpoint"]=="small_09"
    r,pass_last=attempts(prot,folder,c,result,links)
    f=folder/"cache_packet.json"
    center=folder/"point_enclosure.json"
    links[name(f)]=result["cache_packet_sha256"]
    links[name(center)]=result["point_enclosure_sha256"]
    packet=pinned(f,links[name(f)])
    point=pinned(center,links[name(center)])
    assert packet["case_id"]==point["case_id"]==ident
    assert packet["bindings"]==point["bindings"]==result["bindings"]
    assert packet["status"]=="FROZEN_FRESH_EXACT_BINARY_QUADRATURE_CACHE"
    assert packet["closed_parameter_cell_exact"]==c["closed_parameter_cell_exact"]
    assert packet["parameter_order"]==c["parameter_order"]
    assert point["expected_and_rebuilt_site_counts"]==COUNTS
    assert point["cache_packet_sha256"]==result["cache_packet_sha256"]
    assert packet["frozen_binary_cache"]["sha256"]==result["embedded_binary_cache_sha256"]
    assert Q(point["eta_upper_exact"])==Q(result["eta_upper_exact"])
    assert Q(point["center_Z0_upper_exact"])==Q(result["center_Z0_upper_exact"])
    assert all(pt in c["drive"] and c["drive"][pt]==packet[pt] for pt in c["drive"])
    # Translate only the explicit identity-field spelling for the proven
    # helper. No path, hash, model, coordinates or drive is substituted.
    old_spelling=c|{"pilot_id":ident,"exact_drive":c["drive"]}
    arrays,strength,_=helper.cache_extract(prot,packet|{"pilot_id":ident},old_spelling,folder,links)
    assert strength==Q(c["drive"]["common_gain_exact"])*Q(78208913,100000000)
    assert strength!=Q(c["drive"]["activation_strength_binary_exact"])
    ad.A.cache_clear()
    g,H,cache,count=ad.assemble(ad.terms_from_arrays(arrays,packet["parameters"],strength),N)
    assert count==COUNTS and len(cache)==sum(COUNTS.values())
    rg,eta,z0,dz,Rnorm,domain=helper.replay_bound(ad,g,H,cache,arrays["preconditioner"],r,COUNTS)
    assert set(domain["strict_lower_margins"])==helper.MARGINS
    assert all(Q(v["lower_exact"])>0 for v in domain["strict_lower_margins"].values())
    assert Q(domain["strict_lower_margins"]["bulk_J_min"]["lower_exact"])>Q(1,10**6)
    assert Rnorm<=Q(point["preconditioner_infinity_norm_upper_exact"])
    arrays_info={k:sha256(v.tobytes(order="C")).hexdigest() for k,v in arrays.items()}
    own=[ad.pack(rg[i,0]) for i in range(N)]
    for v,w in zip(own,point["preconditioned_center_gradient_by_row"],strict=True):
        assert Q(v["lower_exact"])<=Q(w["upper_exact"]) and Q(w["lower_exact"])<=Q(v["upper_exact"])
    z=z0+dz
    lhs=eta+z*r
    goal=Q(c["e_goal_dimensionless_exact"])
    ok=z<1 and lhs<r<goal
    bound=eta/(1-z) if ok else None
    if ok:
        assert bound<r
    out={"status":"PASS_INDEPENDENT_M_CACHED_LOCAL_ROOT" if ok else "UNRESOLVED_INDEPENDENT_M_CACHED_ROOT",
         "case_id":ident,"record_id":c["record_id"],"field":c["field"],"field_sha256":c["field_sha256"],
         "cache_packet_sha256":result["cache_packet_sha256"],"cache_array_contents_sha256":packet["frozen_binary_cache"]["sha256"],
         "cache_tensor_content_sha256":arrays_info,"exact_strength":str(strength),"rounded_strength_used_in_proof":False,
         "site_counts":count,"analytic_domain":domain,"signed_center_force_intervals":own,
         "producer_candidate_at_same_radius":pass_last,"radius_exact":str(r),"observer_goal_exact":str(goal),
         "independent_eta_upper_exact":str(eta),"independent_Z0_upper_exact":str(z0),
         "independent_variation_upper_exact":str(dz),"independent_ZQ_upper_exact":str(z),
         "preconditioner_infinity_norm_exact":str(Rnorm),"strict_selfmap_lhs_upper_exact":str(lhs),
         "independent_posterior_error_upper_exact":str(bound) if bound is not None else None,
         "display":{"eta":float(eta),"z":float(z),"lhs":float(lhs),"root_displacement_dimensionless":float(bound) if ok else None},
         "scope":SCOPE,"seconds":time.perf_counter()-t0}
    return out


def run(args):
    p=pinned(HERE/"protocol.json",args.protocol_sha256)
    pf=pinned(HERE/"preflight.json",args.preflight_sha256)
    assert p==finite_plan() and pf["status"]=="PASS" and pf["source_sha256"]==sha(SOURCE)
    assert pf["protocol_sha256"]==args.protocol_sha256
    assert pf["actual_producer_payloads_opened"]==pf["actual_site_ad_replays"]==0
    for f,d in p["input_hashes"].items():
        assert sha(ROOT/f)==d
    a,prot,links,cases=admit(p,args.admission.resolve(),args.admission_sha256)
    spec=importlib.util.spec_from_file_location("prior_independent_cached_helpers",PRIOR)
    helper=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(helper)
    assert helper.ROOT==ROOT and helper.PINS[helper.AD_SOURCE]==PINS[AD]
    ad=helper.load_ad()
    import torch
    torch.set_num_threads(1)
    assert torch.get_num_threads()==1
    assert all(importlib.metadata.version(k)==v for k,v in prot["runtime_packages"].items())
    t0=time.perf_counter()
    REPLAY.mkdir(exist_ok=False)
    tag={"pid":os.getpid(),"source_sha256":sha(SOURCE),"protocol_sha256":args.protocol_sha256,
         "preflight_sha256":args.preflight_sha256,"producer_admission_sha256":args.admission_sha256}
    save(REPLAY/"started.json",tag|{"status":"RUNNING"})
    terminal=tag|{"status":"FAILED","exit_code":1}
    try:
        results=[]
        for d,r in cases:
            out=cached_case(helper,ad,prot,d,r,links)
            f=REPLAY/f"{out['case_id']}.json"
            h=save(f,out)
            results.append({k:out[k] for k in ("case_id","status","record_id","display","field_sha256",
                                    "cache_packet_sha256","radius_exact","independent_posterior_error_upper_exact")}
                           |{"native_file":name(f),"native_sha256":h})
            print(json.dumps(results[-1]),flush=True)
        assert len(results)==2 and all(sha(ROOT/f)==h for f,h in links.items())
        n=sum(row["status"]=="PASS_INDEPENDENT_M_CACHED_LOCAL_ROOT" for row in results)
        receipt=save(REPLAY/"receipt.json",{"status":"PASS_INDEPENDENT_TWO_M_CACHED_LOCAL_ROOTS" if n==2 else
                "INDEPENDENT_TWO_M_CACHED_REPLAY_WITH_UNRESOLVED_CASE","scope":SCOPE,"source_and_phase":tag,
                "producer_result_bytes_and_payloads":links,
                "producer_supervisor_pid":a["supervisor_pid"],"producer_supervisor_exec_session":a["unified_exec_session"],
                "producer_observed_exit_code":a["observed_exit_code"],"results":results,"cached_local_roots":n,
                "spatial_root_count":0,"new_equilibrium_solves_or_preconditioners":0,
                "checks":{"current_two_authoritative_M_cases":True,"same_exact_drive_and_fresh_binary_R":True,
                          "independent_ad_for_every_site":True,"own_exact_contraction_on_both_cubes":n==2},
                "seconds":time.perf_counter()-t0})
        terminal|={"status":"FINISHED","exit_code":0,"receipt_sha256":receipt}
    except BaseException as error:
        save(REPLAY/"failure.json",{"message":str(error),"traceback":traceback.format_exc()})
        raise
    finally:
        save(REPLAY/"terminal.json",terminal)
        print(terminal,flush=True)


if __name__=="__main__":
    a=argparse.ArgumentParser()
    q=a.add_subparsers(dest="action",required=True)
    q.add_parser("prepare")
    r=q.add_parser("run")
    r.add_argument("--admission",type=Path,required=True)
    r.add_argument("--admission-sha256",required=True)
    r.add_argument("--protocol-sha256",required=True)
    r.add_argument("--preflight-sha256",required=True)
    opts=a.parse_args()
    prepare() if opts.action=="prepare" else run(opts)
