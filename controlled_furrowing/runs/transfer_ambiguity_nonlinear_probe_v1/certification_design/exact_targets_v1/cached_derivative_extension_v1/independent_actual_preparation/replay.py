"""Independently extend the cached-quadrature derivative to the larger cube."""
import argparse
from fractions import Fraction as Q
from hashlib import sha256
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
BASE=HERE.parent
EXACT=BASE.parent
DESIGN=EXACT.parent
ROOT=DESIGN.parents[2]
PRE=BASE/"preparation"
ACTUAL=BASE/"actual"
RUN=HERE/"actual_replay"
PRIOR=DESIGN/"cached_pilots_v1/independent_actual_preparation/replay.py"
AD=ROOT/"scripts/review_algebraic_root.py"
DX=ROOT/"scripts/transfer_ambiguity_cached_derivative_extension.py"
RADIUS=Q(1,2**34)
COUNTS={"bulk":9216,"surface":288}
N=611
PINS={
    DX:"7ad0b7d4d42efe634123989c6092ccf768cedb92a79244cc7aaa5c10a62f370c",
    PRE/"protocol.json":"4af874d6e45836d11b608ae2929ff5120088fdbd549fd3615398c0068b29ebd7",
    PRE/"preflight.json":"899d9c9486aa0412934bf8866d24c68842c0d2085c380089e1852d906f343f4f",
    PRE/"prelaunch_independent_review.json":"a831fa55d0b30196360c27679f3cfffd18d4fdbdcd3ae2143e6837a9e823a594",
    PRIOR:"3f71e5cae5ee314e260e3bef6c478a923fb7c28f446c184c53d912d6ae296a94",
    AD:"6af2311d0667a28af263b6fb6237063481aa8e65209ae9722c2b0246a8cff27a",
    PRIOR.parent/"actual_replay/receipt.json":"4e891afd79dd2f73508d1f878ddc1113ef107bdff8d9e0f673894f4d2f10ae17",
    EXACT/"independent_actual/replay/receipt.json":"884e0bc76a27465c86b62a202cab491d764aedf030623592238e09adee86d77d",
}
OUTPUT={
    ACTUAL/"started.json":"83cdf9d16dcef98676c80b50e09d317cae70ef8c3410a35dfd1562bad1eff5b0",
    ACTUAL/"results.json":"d9a5c18d6117ba64d0040f371e8e749363e8ceff070fe30fe88b2f95b1348275",
    ACTUAL/"terminal.json":"4cee960506c06a997371b10e5e23eb816f2411189efce7d907acd0dca2d46c5c",
}
SCOPE=("Uniform derivative enclosure of the four exact-real binary Q8/S12 cache potentials on their "
       "own closed r=2^-34 free-coordinate cubes. Exact rational drives; prior independently "
       "verified point eta/Z0 and frozen binary R. No spatial force, exact spatial derivative, "
       "spatial integration, energy-stability, continuum or source-data claim.")


def name(f):
    return f.relative_to(ROOT).as_posix()


def sha(f):
    return sha256(f.read_bytes()).hexdigest()


def pinned(f,h):
    raw=f.read_bytes()
    assert sha256(raw).hexdigest()==h,f
    return json.loads(raw)


def save(f,doc):
    with f.open("x",encoding="utf-8",newline="\n") as w:
        json.dump(doc,w,indent=2,allow_nan=False)
        w.write("\n")
    return sha(f)


def yes(d):
    return bool(d) and all(v is True for v in d.values())


def plan():
    assert all(sha(f)==h for f,h in PINS.items())
    producer=pinned(PRE/"protocol.json",PINS[PRE/"protocol.json"])
    assert len(producer["input_sha256"])==194
    assert len(producer["ordered_cases"])==4
    assert Q(producer["larger_radius_exact"])==RADIUS
    assert producer["whole_box_site_counts"]==COUNTS
    return {"status":"FROZEN_INDEPENDENT_LARGER_RADIUS_DERIVATIVE_REPLAY",
            "source_sha256":sha(SOURCE),"prior_pins":{name(f):h for f,h in PINS.items()},
            "fixed_root_observed_producer":{"supervisor_pid":28720,"unified_exec_session":51655,"observed_os_exit_code":0},
            "producer_admission_pins":{name(f):h for f,h in OUTPUT.items()},
            "pilot_order":[x["pilot_id"] for x in producer["ordered_cases"]],
            "site_counts":COUNTS,"radius_exact":str(RADIUS),"precision_bits":256,
            "scope":SCOPE}


def prepare():
    p=plan()
    assert not RUN.exists()
    # Independent finite-dimensional row majorant and strict radius reserve fixture.
    rows=[[Q(-3),Q(1,2)],[Q(1,4),Q(-2)]]
    V=[Q(1,8),Q(1,32)]
    propagated=[sum(abs(a)*b for a,b in zip(row,V,strict=True)) for row in rows]
    assert propagated==[Q(25,64),Q(3,32)]
    z=Q(1,3)
    h=1-z
    f=h*RADIUS/4
    assert f+(z+h/2)*RADIUS==RADIUS-f<RADIUS
    g=save(HERE/"protocol.json",p)
    h=save(HERE/"preflight.json",{"status":"PASS","source_sha256":sha(SOURCE),"protocol_sha256":g,
                                 "actual_payloads_opened":0,"actual_site_AD":0,"synthetic_exact_fixture":True})
    print(g,h,sha(SOURCE),flush=True)


def admit(p):
    producer=pinned(PRE/"protocol.json",PINS[PRE/"protocol.json"])
    gate=pinned(PRE/"prelaunch_independent_review.json",PINS[PRE/"prelaunch_independent_review.json"])
    assert gate["status"]=="PASS" and yes(gate["checks"])
    bind={k:v for k,v in gate.items() if k in (
          "source_sha256","protocol_sha256","preflight_sha256","actual_geometry_review_sha256","cached_root_independent_review_sha256")}
    assert len(bind)==5
    bind["prelaunch_review_sha256"]=PINS[PRE/"prelaunch_independent_review.json"]
    start,results,terminal=[pinned(ACTUAL/f"{k}.json", OUTPUT[ACTUAL/f"{k}.json"]) for k in ("started","results","terminal")]
    obs=p["fixed_root_observed_producer"]
    assert obs["observed_os_exit_code"]==0 and start["pid"]==terminal["pid"]==obs["supervisor_pid"]
    assert start["bindings"]==terminal["bindings"]==results["bindings"]==bind
    assert terminal["status"]=="FINISHED" and terminal["results_sha256"]==OUTPUT[ACTUAL/"results.json"]
    assert results["status"]=="COMPLETE_FOUR_CACHED_DERIVATIVE_ATTEMPTS_UNREVIEWED"
    assert results["exact_spatial_bridge_admitted"] is False
    order=p["pilot_order"]
    assert start["pilot_order"]==order
    assert [x["pilot_id"] for x in terminal["workers"]]==[x["pilot_id"] for x in results["case_results"]]==order
    assert len({x["pid"] for x in terminal["workers"]})==4
    files={name(f):h for f,h in OUTPUT.items()} | producer["input_sha256"] | p["prior_pins"]
    assert all(sha(ROOT/f)==h for f,h in files.items())
    rows=[]
    for case,top,w in zip(producer["ordered_cases"],results["case_results"],terminal["workers"],strict=True):
        pid=case["pilot_id"]
        assert pid==top["pilot_id"]==w["pilot_id"] and re.fullmatch("[a-z0-9_]+",pid)
        d=ACTUAL/pid
        e=pinned(d/"observed_os_exit.json",w["os_exit_receipt_sha256"])
        assert w==e|{"os_exit_receipt_sha256":w["os_exit_receipt_sha256"]}
        assert e["supervisor_pid"]==obs["supervisor_pid"] and type(e["observed_exit_code"]) is int and e["observed_exit_code"]==0
        t=pinned(d/"terminal.json",e["terminal_sha256"])
        assert t["status"]=="FINISHED" and t["results_sha256"]==top["results_sha256"]
        key={k:e[k] for k in ("pid","pilot_id","supervisor_pid","bindings")}
        assert all(t[k]==v for k,v in key.items()) and e["bindings"]==bind
        for label in ("started","supervisor_child_started"):
            f=d/f"{label}.json"
            files[name(f)]=sha(f)
            authored=pinned(f,files[name(f)])
            assert all(authored[k]==v for k,v in key.items())
        res=pinned(d/"results.json",top["results_sha256"])
        assert res["case"]==case and res["bindings"]==bind
        assert res["status"]==top["status"]=="FRESH_LARGER_RADIUS_CACHED_DERIVATIVE_UNREVIEWED"
        assert res["whole_box_controls_only_cached_quadrature"] is True
        assert res["exact_spatial_derivatives_computed"] is res["exact_spatial_root_certified"] is False
        assert res["independent_derivative_review_required"] is True
        assert Q(res["whole_box_cache_radius_exact"])==RADIUS
        assert res["ZQ_upper_exact"]==top["ZQ_upper_exact"]
        center=d/"center_and_analytic_cache.json"
        files[name(center)]=sha(center)
        assert pinned(center,files[name(center)])==res["center"]
        files.update({name(d/"observed_os_exit.json"):w["os_exit_receipt_sha256"],
                      name(d/"terminal.json"):e["terminal_sha256"],name(d/"results.json"):top["results_sha256"]})
        rows.append((case,res,top["results_sha256"]))
    return producer, files, rows


def prior_module():
    spec=importlib.util.spec_from_file_location("prior_cached_pilot_replay",PRIOR)
    helper=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(helper)
    return helper


def actual_one(helper,ad,ctx,case,result,ph,files):
    started=time.perf_counter()
    assert ctx.drive==Q(case["actual_exact_drive"])!=Q(case["rounded_solver_strength"])
    assert ctx.pilot_id==case["pilot_id"] and ctx.record_id==case["record_id"]
    assert sha256(ctx.R.tobytes(order="C")).hexdigest()==case["preconditioner_binary_sha256"]
    assert ctx.target.packet["frozen_binary_cache"]["sha256"]==case["embedded_cache_sha256"]
    point=pinned(ROOT/case["verified_point"],case["verified_point_sha256"])
    assert point["status"]=="PASS_INDEPENDENT_CACHED_LOCAL_ROOT"
    assert Q(point["exact_strength"])==ctx.drive and point["rounded_strength_used_in_proof"] is False
    assert point["cache_packet_sha256"]==case["cache_packet_sha256"] and point["field_sha256"]==case["field_sha256"]
    assert point["cache_array_contents_sha256"]==case["embedded_cache_sha256"]
    eta=Q(point["independent_eta_upper_exact"])
    z0=Q(point["independent_Z0_upper_exact"])
    assert eta==Q(case["point_eta_upper_exact"]) and z0==Q(case["center_Z0_upper_exact"])
    assert Q(point["radius_exact"])==Q(1,2**36)<RADIUS
    box=pinned(ROOT/case["box"],case["box_sha256"])
    assert Q(box["corrector_radius_exact"])==RADIUS and box["root_membership_assumed"] is False
    a=ctx.target.a
    n=len(a["free"])
    assert n==N and n==len(a["u"]) and a["reference"].shape==(325,2) and a["triangles"].shape==(144,6)
    assert a["shape"].shape==(64,6) and a["shape_grad"].shape==(144,64,6,2)
    assert all(not row.flags.writeable for row in a.values())
    ad.ctx.prec=256
    ad.A.cache_clear()
    # Only site Hessians at the center are required. The prior independently
    # verified point Rg and Z0 are valid at this unchanged center.
    cache=[]
    counts=dict.fromkeys(COUNTS,0)
    for term in ad.terms_from_arrays(a,ctx.target.packet["parameters"],ctx.drive):
        assert term["weight"]>0 and term["indices"]
        fn=ad.bulk_AD if term["kind"]=="bulk" else ad.surface_AD
        _,H,_,_=fn(term["q"],*term["parameters"])
        counts[term["kind"]]+=1
        term["Hcenter"]=H
        term["absB"]=[[ad.UP(x) for x in row] for row in term["B"]]
        term["rowsumB"]=[sum(row,ad.arb(0)).upper() for row in term["absB"]]
        cache.append(term)
    assert counts==COUNTS and len(cache)==9504
    vr,domain=ad.variation(cache,N,ad.A(RADIUS))
    assert domain["counts"]==COUNTS
    assert set(domain["strict_lower_margins"])==helper.MARGINS
    assert all(Q(b["lower_exact"])>0 for b in domain["strict_lower_margins"].values())
    assert Q(domain["strict_lower_margins"]["bulk_J_min"]["lower_exact"])>Q(1,10**6)
    v=[helper.qhi(vr[i,0]) for i in range(N)]
    assert min(v)>=0
    abs_rows=[[abs(Q(x)) for x in row] for row in ctx.R.tolist()]
    norm=max(sum(row) for row in abs_rows)
    assert norm==Q(point["preconditioner_infinity_norm_exact"])==Q(case["R_inf_exact"])>0
    weighted=[sum(x*y for x,y in zip(row,v,strict=True)) for row in abs_rows]
    dz=max(weighted)
    z=z0+dz
    used=Q(result["ZQ_upper_exact"])
    v0=list(map(Q,result["whole_row_variation_upper_exact"]))
    d0=list(map(Q,result["preconditioned_row_variation_upper_exact"]))
    assert len(v0)==len(d0)==N and min(v0)>=0 and min(d0)>=0
    assert all(sum(x*y for x,y in zip(row,v0,strict=True))<=upper for row,upper in zip(abs_rows,d0,strict=True))
    assert used==z0+max(d0)==z0+Q(result["preconditioned_variation_upper_exact"])
    assert max(z,used)<1
    assert eta+z*RADIUS<RADIUS<Q(case["observer_goal_exact"])
    h=1-max(used,z)
    b=h/2
    f=h*RADIUS/4
    assert f+(max(used,z)+b)*RADIUS==RADIUS-f<RADIUS
    detail={"status":"PASS_INDEPENDENT_LARGER_RADIUS_CACHED_DERIVATIVE","pilot_id":case["pilot_id"],
            "producer_result_sha256":ph,"box_sha256":case["box_sha256"],"cache_packet_sha256":case["cache_packet_sha256"],
            "radius_exact":str(RADIUS),"proof_drive_exact":str(ctx.drive),"precision_bits":256,
            "certified_site_counts":counts,"prior_point_eta_upper_exact":str(eta),"prior_center_Z0_upper_exact":str(z0),
            "preconditioner_infinity_norm_exact":str(norm),"whole_box_domains":domain,
            "own_whole_row_bounds_exact":list(map(str,v)),"own_mapped_row_bounds_exact":list(map(str,weighted)),
            "own_variation_upper_exact":str(dz),"reviewed_ZQ_upper_exact":str(z),
            "producer_ZQ_upper_exact":str(used),"z_display":float(z),
            "conservative_headroom_exact":str(h),"planning_b_cap_exact":str(b),
            "planning_total_signed_force_cap_exact":str(f),
            "own_cached_selfmap_lhs_upper_exact":str(eta+z*RADIUS),
            "own_cached_displacement_upper_exact":str(eta/(1-z)),
            "scope":SCOPE,"seconds":time.perf_counter()-started}
    return detail


def run(p_hash,pf_hash):
    p=pinned(HERE/"protocol.json",p_hash)
    assert p==plan()
    pf=pinned(HERE/"preflight.json",pf_hash)
    assert pf["status"]=="PASS" and pf["source_sha256"]==sha(SOURCE) and pf["protocol_sha256"]==p_hash
    assert pf["actual_payloads_opened"]==pf["actual_site_AD"]==0
    producer,files,cases=admit(p)
    helper=prior_module()
    assert helper.ROOT==ROOT
    ad=helper.load_ad()
    import transfer_ambiguity_exact_targets as exact
    RUN.mkdir(exist_ok=False)
    mark={"pid":os.getpid(),"source_sha256":sha(SOURCE),"protocol_sha256":p_hash,"preflight_sha256":pf_hash}
    save(RUN/"started.json",mark|{"status":"RUNNING"})
    terminal=mark|{"status":"FAILED","exit_code":1}
    try:
        _,man,data=exact.prepared()
        out=[]
        with exact.no_integrals():
            for c,res,h in cases:
                ctx,_,entry=exact.open_context(c["pilot_id"],man,data)
                assert c["model_contract_sha256"]==entry["model_contract_sha256"]
                case=actual_one(helper,ad,ctx,c,res,h,files)
                file=RUN/f"{c['pilot_id']}.json"
                h=save(file,case)
                out.append({k:case[k] for k in ("pilot_id","reviewed_ZQ_upper_exact","producer_result_sha256",
                                                "box_sha256","cache_packet_sha256","radius_exact")}
                           | {"native_case_file":name(file),"native_case_sha256":h})
                print(json.dumps({"pilot_id":c["pilot_id"],"ZQ":case["z_display"],
                                  "conservative_h":float(Q(case["conservative_headroom_exact"]))}),flush=True)
        assert len(out)==4
        assert all(sha(ROOT/f)==h for f,h in files.items())
        digest=save(RUN/"receipt.json",{"status":"PASS_INDEPENDENT_FOUR_LARGER_CACHED_DERIVATIVES",
                 "source_and_phase":mark,"scope":SCOPE,"frame":"Whole closed free-coordinate infinity-norm cube",
                 "producer_exit":p["fixed_root_observed_producer"],"producer_and_input_sha256":files,
                 "checks":{"all_four_same_verified_centers_binary_R_and_exact_drives":True,
                           "all_38016_independent_local_cache_sites_replayed":True,
                           "all_four_cached_analytic_domains_positive":True,
                           "all_four_independent_ZQ_below_one_and_cached_selfmap_strict":True,
                           "all_producer_terminal_and_case_pins_verified":True},
                 "reviewed_cases":out,"spatial_force_hessian_or_solve_jobs":0})
        terminal.update(status="FINISHED",exit_code=0,receipt_sha256=digest)
    except BaseException as error:
        save(RUN/"failure.json",{"error":str(error),"traceback":traceback.format_exc()})
        raise
    finally:
        save(RUN/"terminal.json",terminal)
        print(terminal,flush=True)


def bridge(exit_session, native_sha, terminal_sha):
    native=pinned(RUN/"receipt.json",native_sha)
    t=pinned(RUN/"terminal.json",terminal_sha)
    assert native["status"]=="PASS_INDEPENDENT_FOUR_LARGER_CACHED_DERIVATIVES" and yes(native["checks"])
    assert t["status"]=="FINISHED" and type(t["exit_code"]) is int and t["exit_code"]==0 and t["receipt_sha256"]==native_sha
    assert t["pid"]==native["source_and_phase"]["pid"]
    assert exit_session>0
    for c in native["reviewed_cases"]:
        p=pinned(ROOT/c["native_case_file"],c["native_case_sha256"])
        assert all(p[k]==v for k,v in c.items() if k not in ("native_case_file","native_case_sha256"))
        assert max(Q(p["reviewed_ZQ_upper_exact"]), Q(p["producer_ZQ_upper_exact"]))<1
    o=save(RUN/"os_exit_observation.json",{"status":"PASS_OBSERVED_INDEPENDENT_DERIVATIVE_REPLAY_TRUE_OS_EXIT",
           "pid":t["pid"],"unified_exec_session":exit_session,"observed_exit_code":0,
           "terminal_sha256":terminal_sha,"receipt_sha256":native_sha})
    digest=save(HERE/"bridge_admission.json",{
        "status":"PASS_REVIEWED_LARGER_CACHED_DERIVATIVE_BRIDGE_INPUTS",
        "checks":{"native_independent_numeric_review_and_observed_exit":True,"all_four_cached_ZQ_below_one":True,
                  "matching_box_exact_drive_and_fixed_point_cache":True},
        "radius_exact":str(RADIUS),"native_independent_review_file":name(RUN/"receipt.json"),
        "native_independent_review_sha256":native_sha,
        "native_independent_terminal_sha256":terminal_sha,"native_independent_os_exit_sha256":o,
        "derivative_results_sha256":OUTPUT[ACTUAL/"results.json"],"derivative_terminal_sha256":OUTPUT[ACTUAL/"terminal.json"],
        "reviewed_cases":native["reviewed_cases"],"scope":SCOPE})
    print(digest,o,flush=True)


if __name__=="__main__":
    x=argparse.ArgumentParser()
    sub=x.add_subparsers(dest="stage",required=True)
    sub.add_parser("prepare")
    r=sub.add_parser("run")
    r.add_argument("--protocol-sha256",required=True)
    r.add_argument("--preflight-sha256",required=True)
    s=sub.add_parser("seal")
    s.add_argument("--observed-replay-session",type=int,required=True)
    s.add_argument("--native-sha256",required=True)
    s.add_argument("--terminal-sha256",required=True)
    a=x.parse_args()
    if a.stage=="prepare":
        prepare()
    elif a.stage=="run":
        run(a.protocol_sha256,a.preflight_sha256)
    else:
        bridge(a.observed_replay_session,a.native_sha256,a.terminal_sha256)
