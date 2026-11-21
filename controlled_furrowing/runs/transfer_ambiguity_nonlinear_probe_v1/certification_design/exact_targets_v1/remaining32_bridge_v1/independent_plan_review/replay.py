"""Independently replay a sealed first-wave selector plan without moments or AD."""
import argparse
from collections import Counter
from contextlib import ExitStack
from fractions import Fraction as Q
from hashlib import sha256
import importlib.util
import json
import os
from pathlib import Path
import sys
import time
import traceback
from types import SimpleNamespace
from unittest.mock import patch

HERE=Path(__file__).resolve().parent
BR=HERE.parent
ROOT=next(r for r in BR.parents if (r/"PLAN.md").is_file())
EX=BR.parent
OLD=EX/"bridge_v1/independent_plan_plus_L/audit.py"
SRC=ROOT/"scripts/transfer_ambiguity_remaining32_exact_bridge.py"
JOIN=BR/"admission/bridge_input_admission.json"
COLD=BR/"independent_prelaunch/review.json"
P=BR/"preparation"
BASIS={
    SRC:"4d52e8fb8dbc01ba58c38052bbff537443693234499a9a3b9b3e138c82d9b51b",
    P/"protocol.json":"ade6158b96bec63665880beaae1447057b36f162863211acf1f8d24f95519664",
    P/"manifest.json":"d42b7260283ba02ccc34b9e9ba69ce76d58507115c22615c90d28fa38020cb41",
    P/"preflight.json":"a6741fb642a4006b126ac71db4fdb964761f5245e2563cd70a0175db6ee18bda",
    COLD:"89f9b5cacf89b4055fc404c00d3e3c42cf10fed20337fe455e996c2944ef73b2",
    JOIN:"2e2e388e798ed509bb52bc72a6169b8b6ee5db55173b3529ac73c3812b23a6cb",
    OLD:"e74c2f2076aaef11e1abcd24d3895a3aa5a1a26f611b7e8b87b057a919fb2d18",
}
FIRST=("max_contrast_change_small_01","max_contrast_change_small_13",
       "max_negative_contrast_change_small_01","max_negative_contrast_change_small_16")
ORDER=tuple((kind,i) for kind,n in (("bulk",144),("surface",24)) for i in range(n))
R=Q(1,2**34)
POINT="EXACT_CENTER_POINT_ONLY"
BOX="ENTIRE_CLOSED_FIXED_PARAMETER_CUBE"


def demand(ok,message):
    if not ok:
        raise ValueError(message)


def sha(raw):
    return sha256(raw).hexdigest()


def rel(p):
    return p.resolve().relative_to(ROOT).as_posix()


def load(path,want,links,parse=True):
    p=(ROOT/path).resolve()
    name=rel(p)
    demand(not isinstance(path,str) or name==path,"invalid selector proof path")
    raw=p.read_bytes()
    got=sha(raw)
    demand(got==want==links.get(name,got),f"changed checked bytes: {name}")
    links[name]=got
    return json.loads(raw) if parse else raw


def save(path,result):
    with path.open("x",encoding="utf-8",newline="\n") as f:
        json.dump(result,f,indent=2)
        f.write("\n")
    return sha(path.read_bytes())


def yes(d):
    return type(d) is dict and bool(d) and all(v is True for v in d.values())


def zero(x):
    return type(x) is int and x==0


def pid(x):
    return type(x) is int and x>0


def verify_targets(doc,links):
    m=doc[P/"manifest.json"]
    joint=doc[JOIN]
    cold=doc[COLD]
    demand(cold["status"]=="PASS_INDEPENDENT_REMAINING32_EXACT_BRIDGE_COLD"
           and yes(cold["checks"]) and cold["selector_only_authorized_in_this_gate"] is True
           and cold["spatial_workers_authorized_in_this_gate"] is False
           and joint["status"]=="PASS_INDEPENDENT_REMAINING32_SAME_TARGET_BRIDGE_INPUTS"
           and yes(joint["checks"]) and joint["spatial_workers_authorized"] is False,
           "admitted work is context and selectors only")
    indexed={e["case_id"]:e for e in m["cases"]}
    joined={e["case_id"]:e for e in joint["cases"]}
    is_chosen=[indexed[k] for k in FIRST]
    demand(len(indexed)==len(joined)==32
           and tuple(k for k in indexed if k in FIRST)==FIRST
           and all(indexed[k]["cached_input_kind"]==joined[k]["cached_input_kind"]=="FRESH_R34" for k in FIRST)
           and all(indexed[k]["trial"]["trial_box_sha256"]==joined[k]["trial_box_sha256"] for k in FIRST),
           "wrong signed first-wave case/box/target")
    for e in is_chosen:
        demand(tuple((c["kind"],c["slot"]) for c in e["component_catalogue"])==ORDER
               and (e["point_force_scope"],e["whole_H_scope"])==(POINT,BOX),
               "not a complete signed catalogue")
    return indexed,joined


def inputs(links):
    doc={p:load(p,h,links,p.suffix==".json") for p,h in BASIS.items()}
    i,j=verify_targets(doc,links)
    return i,j


def independently_size(catalogue,z,norm):
    f=Counter(i for r in catalogue for i in r["global_indices"])
    h=Counter()
    for r in catalogue:
        idx=r["global_indices"]
        h.update({i:len(idx) for i in idx})
    demand(sorted(h)==sorted(f)==list(range(611)) and 0<=z<1 and norm>0,
           "independent all-free projection or cached derivative range failed")
    f0,h0=max(f.values()),max(h.values())
    cap_f=(1-z)*R/4
    cap_h=(1-z)/2
    return {"cF":f0,"cH":h0,"capF":cap_f,"capH":cap_h,
            "local":(cap_f/(10000*norm*f0),cap_h/(10000*norm*h0))}


def validate_external(case,observed,native,plan,start,observed_manifest):
    expected={"case_id":case,"pid":observed["pid"]}
    demand(observed["status"]=="PASS_OBSERVED_REMAINING32_SELECTOR_TRUE_OS_EXIT"
           and zero(observed["observed_exit_code"]) and pid(observed["unified_exec_session"]) and pid(observed["pid"])
           and native["status"]=="FINISHED" and zero(native["exit_code"]) and start["status"]=="RUNNING_SELECTOR_PLAN"
           and plan["status"]=="COMPLETE_REMAINING32_CASE_SELECTOR_PLAN_AWAITS_INDEPENDENT_REVIEW"
           and plan["spatial_workers_authorized"] is False and zero(plan["new_spatial_integrations"])
           and all(all(d[k]==v for k,v in expected.items()) and pid(d["pid"]) for d in (observed,native,plan,start))
           and start["bindings"]==native["bindings"]==plan["bindings"]
           and native["plan_sha256"]==observed["plan_sha256"]==observed_manifest,
           "missing real clean selected-plan terminal or wrong process phase")
    demand(tuple((t["kind"],t["slot"]) for t in plan["components"])==ORDER
           and all(type(t["slot"]) is int and t["file"]==f"{t['kind']}_{t['slot']:03d}.json" for t in plan["components"])
           and plan["whole_box_radius_exact"]==str(R)
           and (plan["point_force_scope"],plan["whole_H_scope"])==(POINT,BOX),
           "wrong typed component scope and requested radius")


def fixtures(document,joined):
    n=Counter()
    for cid in FIRST:
        base=document[cid]
        y=joined[cid]
        exact=independently_size(base["component_catalogue"],Q(y["reviewed_ZQ_upper_exact"]),Q(y["R_inf_exact"]))
        assert exact["cF"]==6 and exact["cH"]==72 and all(q>0 for q in exact["local"])
        assert list(map(str,exact["local"]))==y["localized_force_and_H_tail_caps_exact"]
        assert exact["capF"]+(Q(y["reviewed_ZQ_upper_exact"])+exact["capH"])*R==R-exact["capF"]
        n[base["cached_input_kind"]]+=1
    return {"four_explicit_target_contexts_and_all672_routes":n=={"FRESH_R34":4},
            "independent_exact_strict_budget_scaling":True,
            "prior_independent_order_tail_math_and_guard_source_pinned":True,
            "no_live_plans_opened_during_cold_preparation":True}


def prepare():
    links={}
    a,b=inputs(links)
    c=fixtures(a,b)
    demand(yes(c),"cold selector review fixtures failed")
    src=sha(Path(__file__).read_bytes())
    ph=save(HERE/"protocol.json",{"status":"FROZEN_INDEPENDENT_FIRST4_SELECTOR_PLAN_REPLAY",
           "source_sha256":src,"input_sha256":links,"first_wave":list(FIRST),"component_count_each":168,
           "radius_exact":str(R),"prior_math_reuse":"Pinned +L independent bulk and trace replay functions",
           "actual_plan_payloads_opened":0,"force_H_AD_moment_or_solve_jobs":0})
    pf=save(HERE/"preflight.json",{"status":"PASS","checks":c,"source_sha256":src,"protocol_sha256":ph,
                                  "actual_plan_payloads_opened":0,"force_H_AD_moment_or_solve_jobs":0})
    print(json.dumps({"source_sha256":src,"protocol_sha256":ph,"preflight_sha256":pf}))


def load_math():
    spec=importlib.util.spec_from_file_location("_prior_independent_selector_math",OLD)
    mod=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def run(args):
    links={}
    ph0=load(HERE/"protocol.json",args.protocol_sha256,links)
    pf=load(HERE/"preflight.json",args.preflight_sha256,links)
    demand(ph0["status"]=="FROZEN_INDEPENDENT_FIRST4_SELECTOR_PLAN_REPLAY"
           and ph0["source_sha256"]==pf["source_sha256"]==sha(Path(__file__).read_bytes())
           and pf["status"]=="PASS" and yes(pf["checks"]) and pf["protocol_sha256"]==args.protocol_sha256
           and ph0["first_wave"]==list(FIRST),"wrong independent selector review phase")
    current={}
    catalogue,joined=inputs(current)
    demand(ph0["input_sha256"]==current,"frozen cold selector proof inputs changed")
    links.update(current)
    case=args.case_id
    directory=BR/"plans"/case
    obs=load(directory/"root_os_exit_observation.json",args.observed_os_sha256,links)
    demand(obs["status"]=="PASS_OBSERVED_REMAINING32_SELECTOR_TRUE_OS_EXIT"
           and obs["case_id"]==case and zero(obs["observed_exit_code"])
           and pid(obs["pid"]) and pid(obs["unified_exec_session"]),"full external plan OS0 admission required")
    term=load(directory/"terminal.json",obs["terminal_sha256"],links)
    plan=load(directory/"plan.json",obs["plan_sha256"],links)
    start_path=directory/"started.json"
    starth=obs.get("started_sha256",sha(start_path.read_bytes()))
    start=load(start_path,starth,links)
    validate_external(case,obs,term,plan,start,obs["plan_sha256"])
    for item in plan["components"]:
        load(directory/item["file"],item["sha256"],links,False)
    sys.path.insert(0,str(ROOT/"scripts"))
    import transfer_ambiguity_remaining32_exact_bridge as api
    demand(Path(api.__file__).resolve()==SRC,"reviewer runtime bridge source mismatch")
    math=load_math()
    folder=HERE/case
    folder.mkdir(exist_ok=False)
    stamp={"case_id":case,"pid":os.getpid(),"observed_os_exit_sha256":args.observed_os_sha256,
           "source_sha256":sha(Path(__file__).read_bytes()),"protocol_sha256":args.protocol_sha256,
           "preflight_sha256":args.preflight_sha256}
    save(folder/"started.json",{"status":"RUNNING",**stamp})
    end={"status":"FAILED","exit_code":1,**stamp}
    start_time=time.perf_counter()
    try:
        bound=SimpleNamespace(case_id=case,protocol_sha256=BASIS[P/"protocol.json"],
             manifest_sha256=BASIS[P/"manifest.json"],preflight_sha256=BASIS[P/"preflight.json"],
             review=COLD,review_sha256=BASIS[COLD],admission=JOIN,admission_sha256=BASIS[JOIN])
        d,sizing,signed,bindings=api.admitted_domain(bound)
        domain=joined[case]
        demand(d["entry"]==catalogue[case] and d["entry"]["trial"]["trial_box_sha256"]==plan["trial_box_sha256"]==domain["trial_box_sha256"]
               and plan["point_box_sha256"]==math.digest(d["point_box"]) and plan["signed_cached_RF_center_rows_exact"]==signed
               and plan["bindings"]==bindings and plan["sizing"]==sizing,
               "same independent joined context/raw plan/cached anchor cannot be reconstructed")
        z,norm=map(Q,(domain["reviewed_ZQ_upper_exact"],domain["R_inf_exact"]))
        caps=independently_size(d["entry"]["component_catalogue"],z,norm)
        tail=caps["local"]
        demand(list(map(Q,plan["sizing"]["caps_exact"]))==list(tail) and plan["sizing"]["strictly_tail_only"] is True,
               "different local independent force or Hessian budgets")
        force,scaled=api.kernels()
        summaries=[]
        with ExitStack() as guard:
            guard.enter_context(math.guard(force,scaled))
            for name in ("plan_case","evaluate_reviewed_component"):
                guard.enter_context(patch.object(api,name,side_effect=AssertionError("producer selector/evaluation path in audit")))
            for slot_info,(k,i) in zip(plan["components"],ORDER,strict=True):
                proof=load(directory/slot_info["file"],slot_info["sha256"],links)
                if k=="bulk":
                    n,jmin=math.bulk(d,i,proof,tail,scaled)
                    ords={side:proof[side]["selector"]["degree"] for side in (POINT,BOX)}
                else:
                    n=math.trace(d,i,proof,tail)
                    jmin=None
                    ords={side:proof[side]["order"] for side in (POINT,BOX)}
                summaries.append({"kind":k,"slot":i,"checked_leaf_sha256":slot_info["sha256"],
                                  "point_and_box_selected_orders":ords,"coefficient_hulls_checked":n,
                                  "whole_same_Y_cache_J_lower_exact":str(jmin) if jmin is not None else None})
        checks={"genuine_observed_complete_plan_root_and_typed_168_leaf_digests":True,
                "point_and_full_cube_polynomial_domains_exact_rule_and_free_map":True,
                "all_bulk_exact_tail_orders_scaled_or_axis_routes_and_same_Y_site_guards":True,
                "all_trace_binomial_orders_and_point_box_288_each_cache_line_domains":True,
                "fresh_independent_local_tail_caps_and_common_proof_center":True,
                "no_spatial_moment_force_H_AD_or_solver_jobs_by_guard":True}
        out={"status":"PASS_INDEPENDENT_REMAINING32_ACTUAL_SELECTOR_PLAN",**stamp,
             "record_id":d["entry"]["record_id"],"plan_sha256":obs["plan_sha256"],"terminal_sha256":obs["terminal_sha256"],
             "selector_started_sha256":starth,"bindings":bindings,"checks":checks,
             "independently_verified_component_count":len(summaries),"case_component_review":summaries,
             "observed_selector_pid":obs["pid"],"observed_selector_unified_exec_session":obs["unified_exec_session"],
             "checked_input_sha256":links,"force_H_AD_moment_or_solve_jobs":0,
             "spatial_worker_authorized_by_this_receipt":False,"seconds":time.perf_counter()-start_time,
             "shared_primitives":"Independent +L selector review math, shared reviewed local polynomial, norm, map and analytic tail-bound kernels; producer selectors and expensive kernels guarded."}
        h=save(folder/"review.json",out)
        end.update(status="FINISHED",exit_code=0,review_sha256=h)
        print(json.dumps({"case_id":case,"review_sha256":h,"components":len(summaries)}))
    except BaseException as err:
        save(folder/"failure.json",{"error":str(err),"traceback":traceback.format_exc(),**stamp})
        raise
    finally:
        save(folder/"terminal.json",end)


if __name__=="__main__":
    if not __debug__:
        raise RuntimeError("independent exact selector premise assertions must remain enabled")
    cli=argparse.ArgumentParser()
    sub=cli.add_subparsers(dest="phase",required=True)
    sub.add_parser("prepare")
    actual=sub.add_parser("run")
    actual.add_argument("--case-id",required=True,choices=FIRST)
    for name in ("protocol-sha256","preflight-sha256","observed-os-sha256"):
        actual.add_argument("--"+name,required=True)
    op=cli.parse_args()
    prepare() if op.phase=="prepare" else run(op)

