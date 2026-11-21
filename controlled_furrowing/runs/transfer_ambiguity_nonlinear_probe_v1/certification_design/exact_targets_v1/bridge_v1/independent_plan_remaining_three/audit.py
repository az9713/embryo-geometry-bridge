"""Cold preparation and explicit per-target review of three singleton plans."""
import argparse
import ast
from collections import Counter
from contextlib import ExitStack
from fractions import Fraction as Q
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import sys
from unittest.mock import patch

SOURCE = Path(__file__).resolve()
OUT = SOURCE.parent
ROOT = SOURCE.parents[6]
sys.path.insert(0,str(ROOT/"scripts"))
import transfer_ambiguity_exact_bridge as J

SESSIONS = {
    "max_contrast_change_large_03": 5202,
    "max_negative_contrast_change_small_17": 6260,
    "max_negative_contrast_change_large_03": 13271,
}
SOURCE_PIN = "58ae92f849414f185c4e28bbb2ec3b9d05361391c4d6e3fa71f519162e30e5e7"
GATE = "e9d4f429e467f10cb0d1c3d736fed5b6cd88a195c1873138a8b2d8d7730d416c"
WRAPPER = "269229312793ec101cc301086d634dd4de89c62a5050ccfc60f185cf761dbc7d"
HELPER = J.OUT/"independent_plan_plus_L/audit.py"
HELPER_SHA = "e74c2f2076aaef11e1abcd24d3895a3aa5a1a26f611b7e8b87b057a919fb2d18"
HELPER_NAMES = {"digest","guard","canonical","inside","bulk_orders","site_map_points","bulk","trace_orders","trace"}
ORDER = tuple((kind,k) for kind,n in (("bulk",144),("surface",24)) for k in range(n))
SCOPES = (J.POINT,J.BOX)
E,S,T,C = J.E,J.SUR,J.T,J.C


def sha(path):
    return sha256(path.read_bytes()).hexdigest()


def arithmetic_helpers():
    raw = HELPER.read_bytes()
    assert sha256(raw).hexdigest() == HELPER_SHA
    parsed = ast.parse(raw)
    funcs = [n for n in parsed.body if isinstance(n,ast.FunctionDef) and n.name in HELPER_NAMES]
    assert {n.name for n in funcs} == HELPER_NAMES
    # Reuse pure definitions only. No +L pilot paths, session, review or results.
    names = ("J","Q","json","sha256","E","S","T","C","SCOPES","product","patch","ExitStack")
    scope = {k:globals()[k] for k in names}
    exec(compile(ast.Module(body=funcs,type_ignores=[]),str(HELPER),"exec"),scope)
    return scope


def state_and_derivatives():
    source_gate = J.PREP/"prelaunch_independent_review.json"
    assert sha(J.SOURCE)==SOURCE_PIN and sha(source_gate)==GATE
    state,b = J.plan_gate(source_gate,GATE)
    assert J.DERIV_ADMISSION_SHA == WRAPPER
    native = J.derivative_inputs(J.DERIV_ADMISSION,WRAPPER,state[1])
    return state,b,native


def pilot_binding(pid, base, native):
    d = native[pid]
    chain = {**base,"derivative_admission_sha256":WRAPPER,
             "independent_native_derivative_review_sha256":d["independent_native_sha256"]}
    tails = tuple(Q(d["planning"][k]) for k in ("per_component_force_entry_tail_cap_exact",
                                                "per_component_H_entry_tail_cap_exact"))
    assert min(tails)>0 and Q(d["ZQ_upper_exact"])<1
    return d,chain,tails


def protocol(base, native):
    return {"status":"FROZEN_INDEPENDENT_THREE_SINGLETON_PLAN_READBACKS",
            "source_sha256":sha(SOURCE),"shared_pure_helper_source_sha256":HELPER_SHA,
            "explicit_pilot_and_producer_sessions":SESSIONS,
            "plan_authority":{k:{"bindings":b,"tail_caps_exact":list(map(str,t))}
                              for k in SESSIONS for _,b,t in [pilot_binding(k,base,native)]},
            "requires_root_os_observation_hash_before_plan_or_leaf_payloads":True,
            "actual_plan_payloads_opened":False,
            "scope":"Serial one-pilot pure arithmetic readbacks. Does not open any other pilot or prior +L plan."}


def prepare():
    _,base,native = state_and_derivatives()
    assert len(ORDER)==len(set(ORDER))==168
    hs=arithmetic_helpers()
    a={"q":[{(0,0):Q(1)}],"B":[[{}]]}
    b={"q":[{(0,0):Q(1)}],"B":[[{(0,0):C.A(0).union(C.A(1))}]]}
    jet={"value":{},"gradient":[{}],"hessian_upper":{(0,0):{}}}
    assert hs["inside"](("q","B"),(a,jet),(b,jet))==2
    assert all(not (OUT/k).exists() for k in SESSIONS)
    p=J.save(OUT/"protocol.json",protocol(base,native))
    k=J.save(OUT/"preflight.json",{"status":"PASS_PREPARATION_ONLY","protocol_sha256":p,
        "source_sha256":sha(SOURCE),"checks":{"frozen_source_and_explicit_distinct_pilot_sessions":True,
            "current_cold_gate_and_larger_derivative_rebound":True,"reviewed_pure_definition_reuse_only":True,
            "coefficient_containment_synthetic":True,"no_actual_plan_payloads_or_integrations":True}})
    print(json.dumps({"source_sha256":sha(SOURCE),"protocol_sha256":p,"preflight_sha256":k}))


def admitted_plan(pid, obs_sha, links, derivative, tails):
    folder=J.PLAN/pid
    obs = J.pinned(folder/"root_observed_os_exit.json",obs_sha)
    assert obs["status"]=="PASS_ROOT_OBSERVED_SINGLETON_BRIDGE_PLAN_OS_EXIT"
    assert obs["pilot_id"]==pid and obs["unified_exec_session"]==SESSIONS[pid]
    assert type(obs["observed_os_exit_code"]) is int and obs["observed_os_exit_code"]==0
    assert type(obs["pid"]) is int and obs["pid"]>0
    term=J.pinned(folder/"terminal.json",obs["terminal_sha256"])
    started=J.pinned(folder/"started.json",obs["started_sha256"])
    assert type(term["exit_code"]) is int and term["exit_code"]==0 and term["status"]=="FINISHED"
    assert term["pilot_id"]==started["pilot_id"]==pid and term["pid"]==started["pid"]==obs["pid"]
    assert started["bindings"]==links and term["plan_sha256"]==obs["plan_sha256"]
    plan=J.pinned(folder/"plan.json",obs["plan_sha256"])
    assert plan["status"]=="COMPLETE_POINT_FORCE_WHOLE_H_SELECTOR_PLAN_UNREVIEWED"
    assert plan["pilot_id"]==pid and plan["bindings"]==links and plan["derivative"]==derivative
    assert plan["tail_caps_exact"]==list(map(str,tails))
    assert plan["point_membership"]==J.POINT and plan["H_membership"]==J.BOX
    assert Q(plan["whole_box_radius_exact"])==J.RADIUS and plan["spatial_component_launch_authorized"] is False
    assert tuple((x["kind"],x["slot"]) for x in plan["components"])==ORDER
    for item in plan["components"]:
        assert type(item["slot"]) is int and item["file"]==f"{item['kind']}_{item['slot']:03d}.json"
        assert sha(folder/item["file"])==item["sha256"]
    return folder,obs,plan


def review(args):
    pid=args.pilot_id
    frozen=J.pinned(OUT/"protocol.json",args.protocol_sha)
    pre=J.pinned(OUT/"preflight.json",args.preflight_sha)
    assert pre["status"]=="PASS_PREPARATION_ONLY" and J.sweep.all_true(pre["checks"])
    assert pre["source_sha256"]==frozen["source_sha256"]==sha(SOURCE)
    assert pre["protocol_sha256"]==args.protocol_sha
    state,base,native=state_and_derivatives()
    assert frozen==protocol(base,native)
    derivative,links,tails=pilot_binding(pid,base,native)
    folder,obs,plan=admitted_plan(pid,args.root_observation_sha,links,derivative,tails)
    helpers=arithmetic_helpers()
    dest=OUT/pid
    dest.mkdir(exist_ok=False)
    J.save(dest/"started.json",{"pilot_id":pid,"root_plan_os_observation_sha256":args.root_observation_sha,
                               "source_sha256":sha(SOURCE),"protocol_sha256":args.protocol_sha,
                               "preflight_sha256":args.preflight_sha,"plan_sha256":obs["plan_sha256"]})
    force,scaled=J.kernels({})
    lengths=Counter()
    orders={f"{kind}/{mode}":Counter() for kind in ("bulk","surface") for mode in ("point","box")}
    summary=[]
    min_j=[]
    with helpers["guard"](force,scaled):
        d=J.domain(state,pid)
        assert plan["record_id"]==d["ctx"].record_id and C.ctx.prec==640
        assert plan["context_box_sha256"]==d["entry"]["box_1_sha256"]
        assert plan["point_box_sha256"]==helpers["digest"](d["point_box"])
        for item in plan["components"]:
            k,i=item["kind"],item["slot"]
            component=J.pinned(folder/item["file"],item["sha256"])
            if k=="bulk":
                count,minbox=helpers["bulk"](d,i,component,tails,scaled)
                min_j.append(minbox)
                s=component[J.BOX]["scaled_or_axis_proof"]
                lengths[s["integration_chart"]]+=1
                lengths["off_axis_exact_identities"]+=s["scaled_polynomial_identity_count"]
                for side,mode in zip(("point","box"),SCOPES,strict=True):
                    orders[f"{k}/{side}"][component[mode]["selector"]["degree"]]+=1
            else:
                count=helpers["trace"](d,i,component,tails)
                lengths["surface"]+=1
                for side,mode in zip(("point","box"),SCOPES,strict=True):
                    orders[f"{k}/{side}"][component[mode]["order"]]+=1
            lengths["coefficient_inclusion_checks"]+=count
            summary.append({"kind":k,"slot":i,"leaf_sha256":item["sha256"],"coefficient_checks":count})
            if len(summary) in (36,72,108,144,168):
                print(json.dumps({"pilot_id":pid,"plan_components_recomputed":len(summary),"total":168}),flush=True)
    assert lengths["scaled_positive_R_affine_material"]==132
    assert lengths["axis_vertex"]==lengths["axis_edge_physical"]==6 and lengths["surface"]==24
    assert lengths["off_axis_exact_identities"]==19824 and lengths["coefficient_inclusion_checks"]>0
    assert len(summary)==168 and len(min_j)==144 and min(min_j)>Q(1,10**6)
    assert not (J.ACTUAL/pid).exists()
    result={"status":"PASS_INDEPENDENT_ACTUAL_SINGLETON_BRIDGE_PLAN","pilot_id":pid,"source_sha256":SOURCE_PIN,
            "plan_sha256":obs["plan_sha256"],"plan_terminal_sha256":obs["terminal_sha256"],
            "plan_observed_pid":obs["pid"],"plan_observed_os_exit_code":obs["observed_os_exit_code"],
            "plan_observed_session":obs["unified_exec_session"],
            "plan_observed_started_sha256":obs["started_sha256"],
            "plan_root_observation_sha256":args.root_observation_sha,
            "derivative_admission_sha256":WRAPPER,"audit_source_sha256":sha(SOURCE),
            "audit_protocol_sha256":args.protocol_sha,"audit_preflight_sha256":args.preflight_sha,
            "checks":{"true_root_plan_exit_and_all_168_leaf_digests":True,
                      "all_exact_orders_rhos_and_two_mass_tail_bounds_recomputed":True,
                      "point_and_box_polynomial_containments_and_digests_recomputed":True,
                      "all_9216_point_and_box_bulk_and_576_trace_sites_checked":True,
                      "whole_box_scaled_and_axis_structure_replayed":True,
                      "no_spatial_integration_AD_or_root_claim":True},
            "route_counts_and_identity_checks":dict(lengths),
            "selected_orders":{k:dict(v) for k,v in orders.items()},
            "bulk_same_Y_point_and_box_sites_each":144*64,
            "trace_same_Y_point_and_box_sites_each":24*12,
            "whole_box_bulk_min_line_J_exact":str(min(min_j)),
            "component_readback":summary,
            "scope":"Pure selector readback. Reuses pinned reviewed P2/jet/chart/preflight primitives and explicit pilot authorities. No force, H, moment, cache AD, solve or root inference."}
    digest=J.save(dest/"review.json",result)
    print(json.dumps({"review_sha256":digest,"pilot_id":pid,"route_counts":dict(lengths),
                      "orders":result["selected_orders"]}))


if __name__=="__main__":
    cli=argparse.ArgumentParser()
    sub=cli.add_subparsers(dest="stage",required=True)
    sub.add_parser("prepare")
    act=sub.add_parser("review")
    act.add_argument("--pilot-id",choices=SESSIONS,required=True)
    for name in ("protocol-sha","preflight-sha","root-observation-sha"):
        act.add_argument("--"+name,required=True)
    args=cli.parse_args()
    prepare() if args.stage=="prepare" else review(args)
