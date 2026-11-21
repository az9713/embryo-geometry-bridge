"""Later-wave independent selector replay. No integration or mechanics."""
import argparse
from contextlib import ExitStack
from fractions import Fraction
import importlib.util
import json
import os
from pathlib import Path
import sys
from time import perf_counter
import traceback
from types import SimpleNamespace
from unittest.mock import patch

HERE=Path(__file__).resolve().parent
ROOT=next(p for p in HERE.parents if (p/"PLAN.md").is_file())
BR=HERE.parents[1]
SYS=HERE.parent
sys.path.insert(0,str(SYS))
import wave_execution as E

OLD=BR/"independent_plan_review/replay.py"
APPROVAL=HERE/"root_source_review.json"
STATUS="PASS_ROOT_REMAINING28_SELECTOR_REPLAY_SOURCE"
ADAPTER_REVIEW="bf1af4575dffa36a7ddade7df0137156ec7a9ae610265e4f4d17386b9956149a"
RELEASE="002b471dceed8666aa0f70e1e9673744f49d0bf3a70da65a194c3faf34b0adbe"
PINS={
    OLD:"c4b096ba3d93ffa17abe5eb5a0b20bec2bd41a0c97ddb0fe0b13ed9c99616e9b",
    E.SOURCE:"44c4fea65168a77926f93a08c3817edc0645c94ec7ef625fa34b4656beade8b7",
    E.PREP/"protocol.json":"5057fecf76205bd0eee84cde3aff730b0011299d87c36fdebaad1be53cf10c72",
    E.PREP/"preflight.json":"aa03be5fd4d2341326a830bea1081e421b3739711d51b50d9398dbe6c764e214",
    E.REVIEW_PATH:ADAPTER_REVIEW,
    E.RELEASE_PATH:RELEASE,
    OLD.with_name("root_source_review.json"): "da4939cbad4c64c7ad9a0b3d7dd909dae2bb87e9f76e53d665122bb33194dfc5",
}
require,save,sha,read=E.require,E.save,E.sha,E.read
Q=Fraction


def old_helpers():
    require(sha(OLD)==PINS[OLD],"reviewed first-wave proof primitives")
    spec=importlib.util.spec_from_file_location("_old_independent_selectors",OLD)
    helpers=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(helpers)
    return helpers


def adapter_args(wave):
    return SimpleNamespace(wave=wave,protocol_sha256=PINS[E.PREP/"protocol.json"],
         preflight_sha256=PINS[E.PREP/"preflight.json"],review=E.REVIEW_PATH,review_sha256=ADAPTER_REVIEW,
         release=E.RELEASE_PATH,release_sha256=RELEASE)


def route(profiles,wave,case):
    require(wave in E.WAVES and case in profiles[wave]["case_order"] and case not in E.W.FIRST,
            "selector replay case must be explicitly assigned to this frozen later wave")
    return profiles[wave]["raw_worker_targets"][case]


def inputs():
    pins={E.rel(p):h for p,h in PINS.items() if sha(p)==h}
    require(len(pins)==len(PINS),"reviewed execution and prior independent proof source closure")
    M=old_helpers()
    catalogue,joined=M.inputs(pins)
    profiles,bridge=E.cold()
    for wave,p in profiles.items():
        actual,b=E.release_gate(adapter_args(wave))
        require(actual==p and b["later_waves_release_sha256"]==RELEASE,"same actual first-wave release")
        for c in p["case_order"]:
            require(route(profiles,wave,c)==E.W.case_payload(catalogue[c]),"same signed target")
    return M,catalogue,joined,profiles,pins


def admitted(M,case,outer_pin,links):
    # External observation is loaded first. The producer plan components
    # stay unopened until this complete case-specific process chain is checked.
    folder=E.S.PLAN/case
    obs=M.load(folder/"root_os_exit_observation.json",outer_pin,links)
    require(obs["status"]=="PASS_OBSERVED_REMAINING32_SELECTOR_TRUE_OS_EXIT" and obs["case_id"]==case
            and M.zero(obs["observed_exit_code"]) and M.pid(obs["pid"]) and M.pid(obs["unified_exec_session"]),
            "literal observed same-case producer exit0 required")
    terminal=M.load(folder/"terminal.json",obs["terminal_sha256"],links)
    started=M.load(folder/"started.json",obs["started_sha256"],links)
    plan=M.load(folder/"plan.json",obs["plan_sha256"],links)
    M.validate_external(case,obs,terminal,plan,started,obs["plan_sha256"])
    for item in plan["components"]:
        M.load(folder/item["file"],item["sha256"],links,False)
    return folder,obs,plan


def prepare():
    M,catalogue,joined,profiles,pins=inputs()
    modes={}
    for w,p in profiles.items():
        for case in p["case_order"]:
            frozen=catalogue[case]
            verified=joined[case]
            caps=M.independently_size(frozen["component_catalogue"],Q(verified["reviewed_ZQ_upper_exact"]),Q(verified["R_inf_exact"]))
            require(caps["cF"]==6 and caps["cH"]==72
                    and [str(x) for x in caps["local"]]==verified["localized_force_and_H_tail_caps_exact"]
                    and tuple((v["kind"],v["slot"]) for v in frozen["component_catalogue"])==M.ORDER,
                    "each of 28 closed-cube target proofs has exact independent incidence/tail budgets")
            modes[verified["cached_input_kind"]]=modes.get(verified["cached_input_kind"],0)+1
    require(sum(modes.values())==28,"all 28 covered once")
    first=profiles[E.WAVES[0]]["case_order"][0]
    require(E.V.A.reject(lambda:route(profiles,E.WAVES[1],first)) and
            E.V.A.reject(lambda:route(profiles,E.WAVES[0],E.W.FIRST[0])),"cross-wave and first-wave selectors rejected")
    require(not any("/plans/" in k for k in pins),"no actual selector plan packets opened during cold setup")
    h=sha(Path(__file__))
    pw={w:p["case_order"] for w,p in profiles.items()}
    ph=save(HERE/"protocol.json",{"status":"FROZEN_INDEPENDENT_REMAINING28_ACTUAL_SELECTOR_REPLAY",
         "source_sha256":h,"later_release_sha256":RELEASE,"profiles_sha256":E.PROFILE_SHA,
         "input_sha256":pins,"ordered_case_ids_by_wave":pw,"component_count_each":168,"radius_exact":str(M.R),
         "actual_plan_payloads_opened":0,"force_H_AD_moment_or_solve_jobs":0})
    f=save(HERE/"preflight.json",{"status":"PASS","source_sha256":h,"protocol_sha256":ph,
            "checks":{"same_all_28_immutable_targets_and_component_catalogues":True,
                      "cross_wave_and_first_wave_case_rejected":True,
                      "independent_all28_exact_common_tail_budget_scaling":True,
                      "source_and_original_replay_primitives_pinned_without_rebinding":True,
                      "no_actual_plan_payloads_or_numerical_kernels_opened":True},
            "exact_cached_authority_kind_counts":modes,"actual_plan_payloads_opened":0,
            "force_H_AD_moment_or_solve_jobs":0})
    print(json.dumps({"source_sha256":h,"protocol_sha256":ph,"preflight_sha256":f}),flush=True)


def run(args):
    M,catalogue,joined,profiles,pins=inputs()
    target=route(profiles,args.wave,args.case_id)
    p=M.load(HERE/"protocol.json",args.protocol_sha256,pins)
    f=M.load(HERE/"preflight.json",args.preflight_sha256,pins)
    g=M.load(APPROVAL,args.source_review_sha256,pins)
    bind={"source_sha256":sha(Path(__file__)),"protocol_sha256":args.protocol_sha256,"preflight_sha256":args.preflight_sha256}
    expected={**bind,"later_release_sha256":RELEASE,"profiles_sha256":E.PROFILE_SHA}
    require(p["status"]=="FROZEN_INDEPENDENT_REMAINING28_ACTUAL_SELECTOR_REPLAY"
            and g["status"]==STATUS and f["status"]=="PASS" and M.yes(f["checks"]) and M.yes(g["checks"])
            and all(g[k]==v for k,v in expected.items()) and all(p[k]==g[k] for k in ("source_sha256","later_release_sha256","profiles_sha256"))
            and f["source_sha256"]==bind["source_sha256"] and f["protocol_sha256"]==args.protocol_sha256
            and p["input_sha256"]=={k:v for k,v in pins.items() if not k.startswith(E.rel(HERE)+"/")}
            and p["ordered_case_ids_by_wave"]=={w:u["case_order"] for w,u in profiles.items()}
            and g["new_spatial_integrations_authorized"] is False,"approved current independent later-wave selector replay")
    case=args.case_id
    directory,obs,plan=admitted(M,case,args.observed_os_sha256,pins)
    folder=OLD.parent/case
    folder.mkdir(exist_ok=False)
    stamp={**bind,"case_id":case,"pid":os.getpid(),"wave_id":args.wave,
           "observed_os_exit_sha256":args.observed_os_sha256,"root_source_review_sha256":args.source_review_sha256,
           "later_release_sha256":RELEASE}
    save(folder/"started.json",{"status":"RUNNING",**stamp})
    end={"status":"FAILED","exit_code":1,**stamp}
    t0=perf_counter()
    try:
        api=E.S
        math=M.load_math()
        d,sizing,signed,bindings=api.admitted_domain(E.selector_args(case))
        domain=joined[case]
        require(d["entry"]==catalogue[case] and E.W.case_payload(d["entry"])==target
               and target["trial_box_sha256"]==plan["trial_box_sha256"]==domain["trial_box_sha256"]
               and plan["point_box_sha256"]==math.digest(d["point_box"]) and plan["signed_cached_RF_center_rows_exact"]==signed
               and plan["bindings"]==bindings and plan["sizing"]==sizing,"same exact target, maps and parameter center")
        z,norm=map(Q,(domain["reviewed_ZQ_upper_exact"],domain["R_inf_exact"]))
        tail=M.independently_size(d["entry"]["component_catalogue"],z,norm)["local"]
        require(list(map(Q,plan["sizing"]["caps_exact"]))==list(tail) and
                [str(x) for x in tail]==domain["localized_force_and_H_tail_caps_exact"]
                and plan["sizing"]["strictly_tail_only"] is True,"same independent analytic budgets")
        force,scaled=api.kernels()
        rows=[]
        with ExitStack() as stack:
            stack.enter_context(math.guard(force,scaled))
            for name in ("plan_case","evaluate_reviewed_component"):
                stack.enter_context(patch.object(api,name,side_effect=AssertionError("producer selector/evaluation is excluded")))
            for item,(kind,i) in zip(plan["components"],M.ORDER,strict=True):
                raw=M.load(directory/item["file"],item["sha256"],pins)
                if kind=="bulk":
                    n,lower=math.bulk(d,i,raw,tail,scaled)
                    orders={scope:raw[scope]["selector"]["degree"] for scope in (M.POINT,M.BOX)}
                else:
                    n=math.trace(d,i,raw,tail)
                    lower=None
                    orders={scope:raw[scope]["order"] for scope in (M.POINT,M.BOX)}
                rows.append({"kind":kind,"slot":i,"checked_leaf_sha256":item["sha256"],
                    "point_and_box_selected_orders":orders,"coefficient_hulls_checked":n,
                    "whole_same_Y_cache_J_lower_exact":None if lower is None else str(lower)})
        require(len(rows)==168,"all same-target plan leaves required")
        proof={"status":"PASS_INDEPENDENT_REMAINING32_ACTUAL_SELECTOR_PLAN",**stamp,
            "record_id":target["record_id"],"plan_sha256":obs["plan_sha256"],"terminal_sha256":obs["terminal_sha256"],
            "selector_started_sha256":obs["started_sha256"],"bindings":bindings,
            "checks":{"observed_native_exit_and_exact_168_selector_packet_digests":True,
                      "whole_box_and_point_coefficient_geometry_maps_orders_tails":True,
                      "independently_checked_analytic_budget_at_exact_target_and_cache":True,
                      "guard_no_force_H_AD_moment_or_solve_jobs":True},
            "independently_verified_component_count":len(rows),"case_component_review":rows,
            "observed_selector_pid":obs["pid"],"observed_selector_unified_exec_session":obs["unified_exec_session"],
            "checked_input_sha256":pins,"force_H_AD_moment_or_solve_jobs":0,
            "spatial_worker_authorized_by_this_receipt":False,"seconds":perf_counter()-t0}
        h=save(folder/"review.json",proof)
        end.update(status="FINISHED",exit_code=0,review_sha256=h)
        print(json.dumps({"case_id":case,"review_sha256":h,"pid":stamp["pid"],"components":len(rows)}),flush=True)
    except BaseException as error:
        save(folder/"failure.json",{"error":str(error),"traceback":traceback.format_exc(),**stamp})
        raise
    finally:
        save(folder/"terminal.json",end)


if __name__=="__main__":
    if not __debug__:
        raise RuntimeError("reviewed math requires enabled Python assertions")
    cli=argparse.ArgumentParser()
    sub=cli.add_subparsers(dest="phase",required=True)
    sub.add_parser("prepare")
    actual=sub.add_parser("run")
    actual.add_argument("--wave",required=True,choices=E.WAVES)
    actual.add_argument("--case-id",required=True)
    for n in ("protocol-sha256","preflight-sha256","source-review-sha256","observed-os-sha256"):
        actual.add_argument("--"+n,required=True)
    args=cli.parse_args()
    prepare() if args.phase=="prepare" else run(args)
