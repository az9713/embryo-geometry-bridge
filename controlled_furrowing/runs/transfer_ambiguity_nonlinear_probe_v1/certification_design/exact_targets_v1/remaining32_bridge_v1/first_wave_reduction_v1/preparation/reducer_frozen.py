"""First-wave exact signed reduction of independently valid local hulls."""
import argparse
from collections import Counter
from fractions import Fraction
from itertools import combinations_with_replacement
import os
from pathlib import Path
import traceback
from types import SimpleNamespace

import transfer_ambiguity_exact_bridge_reduction as K
import transfer_ambiguity_remaining32_spatial_worker as W

SOURCE=Path(__file__).resolve()
S=W.S
V=W.V
ROOT=W.ROOT
OUT=S.OUT/"first_wave_reduction_v1"
PREP=OUT/"preparation"
ACTUAL=OUT/"actual"
REVIEW_FILE=PREP/"independent_prelaunch/review.json"
ADMIT_FILE=OUT/"admission/actual_input_admission.json"
HULL=S.OUT/"independent_component_replay_preparation"
HULL_ACTUAL=HULL/"actual_replay"
SPD=S.OUT.parent/"remaining32_preconditioner_spd_v1/independent_actual_preparation"
RAW=W.ACTUAL/"root_terminal_admission.json"
HULL_APPROVAL=HULL/"root_source_review.json"
W_GATE="5463da7c9c36a30007513dce4788f7114878978bb7274d2e484b470c6916e85b"
W_ADMIT="8468da6cd555797298df182d6213a93a4d054029bdfca9b332aea80ca39f6855"
PINS={
 Path(K.__file__).resolve():"d31ca58afe2e1b935d723513fb9fa18f2f978885f8f6ec32e48766f8b15bf0e8",
 W.SOURCE:"61b0082a68a297f3ae287232704a9b72ed4fbf53e96090db990ba0d8a8458f2d",
 W.PREP/"protocol.json":"434c3786447d92c75bd13880d57ad090cbc36c2aff5b795515ceed526d110e42",
 W.PREP/"preflight.json":"2f9c6d5e1f3c88f50389774e6afc67b1e3ed18030c6b3f12af7804cb4f41d918",
 W.INDEPENDENT:W_GATE,W.WAVE_GATE:W_ADMIT,
 HULL/"replay.py":"283044f9d213bdfa82d1333ae9be1e86b8d2dc94f5ceb7d8186e9c6e1a1c0720",
 HULL/"hull_math.py":"f8b4e1e8b4fb9069d235303f3fd0a45ad629ccd9d0bf37bf30ceca509c96b38b",
 HULL/"protocol.json":"7d5f357d87b7bc1d0370b7cab5a1473c4a8a9174523070e461c267e2c9d4423c",
 HULL/"preflight.json":"f461429d67e6214c9dfdc9d80f7539eb1538105c7baa3af7ae930b95a4597ec5",
 S.OUT/"admission/root_hull_design_review.json":"6442253147ebc1a73dff1358f8a44ebba8425ddb3e0af18ba01fb3b6a8feFC03".lower(),
}
REVIEW="PASS_INDEPENDENT_FIRST_WAVE_SPATIAL_REDUCER_COLD"
ADMISSION="PASS_INDEPENDENT_FIRST_WAVE_WHOLE_HULLS_REDUCTION_INPUTS"
HULL_PASS="PASS_INDEPENDENT_ACTUAL_REMAINING32_COMPONENT_HULLS"
LEAF_PASS="PASS_INDEPENDENT_COMPONENT_WITH_EXACT_HULLS"
OS_HULL="PASS_OBSERVED_INDEPENDENT_REMAINING32_HULL_REPLAY_TRUE_OS_EXIT"
SPD_PASS="PASS_INDEPENDENT_REMAINING32_EXACT_SYMMETRIC_PART_SPD"
CHANNELS=("point_force","whole_box_H","point_same_Y_cache_correction")
SCOPE="Conservative exact signed union hulls for the same four fixed spatial targets at r34. Raw producer endpoints or overlap alone supply no proof."
ZERO={"actual_reduction_jobs":0,"new_spatial_integrals":0,"new_AD_H_moments_or_solves":0}
Q=K.rational
require,save,read,sha=S.require,S.save,S.read,S.sha
rel,canonical=S.rel,S.canonical


def worker_args():
    return SimpleNamespace(protocol_sha256=PINS[W.PREP/"protocol.json"],preflight_sha256=PINS[W.PREP/"preflight.json"],
                           review=W.INDEPENDENT,review_sha256=W_GATE,admission=W.WAVE_GATE,admission_sha256=W_ADMIT)


def static_phase():
    require(all(sha(p)==h for p,h in PINS.items()),"unchanged reviewed premises and reducer kernel")
    a,b=W.gate(worker_args())
    p=read(HULL/"protocol.json",PINS[HULL/"protocol.json"])
    f=read(HULL/"preflight.json",PINS[HULL/"preflight.json"])
    h=read(S.OUT/"admission/root_hull_design_review.json",PINS[S.OUT/"admission/root_hull_design_review.json"])
    require(p["first_wave"]==list(W.FIRST) and f["status"]=="PASS" and V.yes(f["checks"])
            and h["status"]=="PASS_ROOT_CONDITIONAL_INDEPENDENT_ENCLOSURE_HULL_DESIGN" and V.yes(h["checks"])
            and h["actual_component_jobs_authorized"] is False and tuple(x["case_id"] for x in a["plans"])==W.FIRST
            and K.RADIUS==S.RADIUS and K.POINT==S.POINT and K.BOX==S.BOX,
            "raw, independent and reducer scope all match no-later-wave authority")
    return a,b


def gate(args):
    a,worker=static_phase()
    p=read(PREP/"protocol.json",args.protocol_sha256)
    f=read(PREP/"preflight.json",args.preflight_sha256)
    require(args.review.resolve()==REVIEW_FILE.resolve() and args.admission.resolve()==ADMIT_FILE.resolve(),
            "canonical independent reducer review and actual input gate")
    r=read(args.review,args.review_sha256)
    joined=read(args.admission,args.admission_sha256)
    bind={"source_sha256":sha(SOURCE),"protocol_sha256":args.protocol_sha256,"preflight_sha256":args.preflight_sha256,
          "first_wave_worker_admission_sha256":W_ADMIT,"joined_bridge_input_admission_sha256":a["joined_bridge_input_admission_sha256"]}
    require(sha(PREP/"reducer_frozen.py")==bind["source_sha256"]==p["source_sha256"]==f["source_sha256"]
            and p["status"]=="FROZEN_FIRST_WAVE_HULL_ONLY_SIGNED_REDUCER" and f["status"]=="PASS" and V.yes(f["checks"])
            and all(V.zero(f[k]) for k in ZERO) and p["worker_runtime_binding"]==worker
            and tuple(p["first_wave"])==W.FIRST and f["protocol_sha256"]==args.protocol_sha256,
            "current exact cold proof wrapper")
    check_gates(r,joined,bind,args.review_sha256)
    return a,worker,joined,{**bind,"reducer_review_sha256":args.review_sha256,"actual_input_admission_sha256":args.admission_sha256}


def check_gates(review,joined,bind,reviewsha):
    require(review["status"]==REVIEW and V.yes(review["checks"])
            and all(review[k]==v for k,v in bind.items()) and review["case_order"]==list(W.FIRST)
            and review["new_spatial_integrals_authorized"] is False,
            "independent cold current global reducer source review")
    require(joined["status"]==ADMISSION and V.yes(joined["checks"])
            and all(joined[k]==v for k,v in bind.items()) and joined["reducer_review_sha256"]==reviewsha
            and tuple(x["case_id"] for x in joined["hull_case_observations"])==W.FIRST
            and joined["all_four_cases_required_before_arithmetic"] is True
            and joined["new_spatial_integrals_authorized"] is False,
            "full first-wave externally verified hull closure only")


def observed_raw(sha_raw,a,worker):
    d=read(RAW,sha_raw)
    require(d["status"]=="PASS_ROOT_OBSERVED_REMAINING32_FIRST_WAVE_SPATIAL_TERMINALS"
            and V.positive_int(d["pid"]) and V.positive_int(d["unified_exec_session"]) and V.zero(d["observed_exit_code"]),
            "full genuine observed raw producer exit required")
    links=d["sha256"]
    tests={rel(W.ACTUAL/(n+".json")) for n in ("started","manifest","terminal")}
    for case,part in W.JOBS:
        folder=W.path_for(case,part)
        tests|={rel(folder/(n+".json")) for n in ("started","launched","observed_os_exit","terminal","manifest")}
        tests|={rel(folder/f"{part}_{i:03d}.json") for i in range(W.COUNTS[part])}
    require(set(links)==tests and len(links)==715,"all raw actual eight-part pins required")
    def doc(path):
        return read(path,links[rel(path)])
    start,end,raw=(doc(W.ACTUAL/(n+".json")) for n in ("started","terminal","manifest"))
    require(all(V.positive_int(x["pid"]) and x["pid"]==d["pid"] and x["stage"]==W.STAGE
                and x["bindings"]==worker and x["jobs"]==[list(job) for job in W.JOBS] for x in (start,end,raw))
            and start["status"]=="RUNNING" and end["status"]=="FINISHED" and V.zero(end["exit_code"])
            and end["manifest_sha256"]==links[rel(W.ACTUAL/"manifest.json")]
            and raw["status"]=="COMPLETE_FOUR_RAW_SPATIAL_BRIDGES_AWAITING_INDEPENDENT_PACKET_REPLAY"
            and raw["raw_components"]==672 and raw["root_asserted"] is False,
            "same complete raw supervisor and true OS chain")
    require(tuple(row["case_id"] for row in raw["cases"])==W.FIRST,"exact raw wave order")
    parts={}
    pids=set()
    for rec in raw["cases"]:
        case=rec["case_id"]
        selected=V.A.pick(a["plans"],case)
        require(rec["target"]==selected["target"] and tuple(p["part"] for p in rec["parts"])==tuple(W.COUNTS),
                "foreign raw case/part")
        for part in rec["parts"]:
            kind=part["part"]
            base=W.path_for(case,kind)
            begun,launch,native,obs,pm=(doc(base/(n+".json")) for n in
                ("started","launched","terminal","observed_os_exit","manifest"))
            pid=part["pid"]
            common={"pid":pid,"supervisor_pid":d["pid"],"stage":W.STAGE,"case_id":case,"part":kind,"bindings":worker}
            require(V.positive_int(pid) and pid not in pids and V.zero(part["observed_exit_code"])
                    and all(V.positive_int(x["pid"]) and V.positive_int(x["supervisor_pid"])
                            and all(x[k]==v for k,v in common.items()) for x in (part,begun,launch,native,obs,pm))
                    and V.zero(obs["observed_exit_code"]) and begun["status"]=="RUNNING"
                    and part["observed_os_exit_sha256"]==links[rel(base/"observed_os_exit.json")]
                    and part["terminal_sha256"]==obs["terminal_sha256"]==links[rel(base/"terminal.json")]
                    and part["manifest_sha256"]==native["results_sha256"]==links[rel(base/"manifest.json")],
                    "literal typed native child identities/true exits")
            V.child_admit(obs["observed_exit_code"],pm,native,common,part["manifest_sha256"])
            pids.add(pid)
            require(pm["status"]==W.STATUS and pm["target"]==selected["target"] and pm["plan"]==selected
                    and pm["packet_independent_replay_required"] is True and pm["radius_exact"]==str(S.RADIUS)
                    and (pm["point_force_scope"],pm["whole_H_scope"])==(S.POINT,S.BOX)
                    and len(pm["components"])==W.COUNTS[kind],"complete raw part and unchanged plan")
            for i,c in enumerate(pm["components"]):
                require(type(c["slot"]) is int and c["slot"]==i and c["kind"]==kind and c["file"]==f"{kind}_{i:03d}.json"
                        and c["sha256"]==links[rel(base/c["file"])],"exact raw leaf provenance")
            parts[(case,kind)]=(base,pm)
    require(len(parts)==len(pids)==8,"all eight separate raw successful native parts")
    return d,parts

def hull_static(h):
    rev=read(HULL_APPROVAL,h)
    static={"source_sha256":PINS[HULL/"replay.py"],
            "local_bound_and_hull_source_sha256":PINS[HULL/"hull_math.py"],
            "protocol_sha256":PINS[HULL/"protocol.json"],"preflight_sha256":PINS[HULL/"preflight.json"]}
    require(rev["status"]=="PASS" and V.yes(rev["checks"]) and all(rev[k]==v for k,v in static.items()),
            "approved independent own-interval mathematical implementation")
    return static


def open_all_hull_metadata(a,admission,raw):
    source_review=admission["hull_source_review_sha256"]
    raw_h=admission["root_producer_admission_sha256"]
    static=hull_static(source_review)
    native={}
    pids=set()
    for spec in admission["hull_case_observations"]:
        case=spec["case_id"]
        opts=V.A.pick(a["plans"],case)
        target=opts["target"]
        folder=HULL_ACTUAL/case
        osr=read(folder/"os_exit_observation.json",spec["observed_os_exit_sha256"])
        require(osr["status"]==OS_HULL and osr["case_id"]==case and V.positive_int(osr["pid"])
                and osr["pid"] not in pids and V.positive_int(osr["unified_exec_session"])
                and V.zero(osr["observed_exit_code"]) and osr["root_producer_admission_sha256"]==raw_h,
                "true typed successful independent actual whole-case replay only")
        pids.add(osr["pid"])
        end=read(folder/"terminal.json",osr["terminal_sha256"])
        r=read(folder/"receipt.json",osr["receipt_sha256"])
        phase={**static,"pid":osr["pid"],"case_id":case,"root_source_review_sha256":source_review,
               "root_producer_admission_sha256":raw_h}
        require(all(V.positive_int(x["pid"]) and all(x[k]==v for k,v in phase.items()) for x in (r,end))
                and end["status"]=="FINISHED" and V.zero(end["exit_code"]) and end["receipt_sha256"]==osr["receipt_sha256"]
                and r["status"]==HULL_PASS and r["full_independent_hulls_computed"] is True
                and r["admit_hull_only_global_reduction_after_separate_review"] is True
                and r["overlap_or_producer_containment_used_as_proof"] is False and r["root_asserted"] is False
                and set(r["disjoint_counts"])==set(CHANNELS) and all(V.zero(n) for n in r["disjoint_counts"].values()),
                "same full native mathematical phase without any disjoint HOLD")
        identity={**{k:target[k] for k in ("record_id","trial_box_sha256","model_contract_sha256",
                                           "exact_drive","R_binary_sha256","center_binary_sha256")},
                  "original_archived_shift_exact":target["archived_center_shift_exact"],
                  "cache_packet_sha256":target["point_cache_sha256"],"radius_exact":str(S.RADIUS),
                  "plan_sha256":opts["plan_sha256"],"plan_review_sha256":opts["independent_plan_review_sha256"],
                  "producer_observed_pid":raw["pid"],"producer_observed_unified_exec_session":raw["unified_exec_session"]}
        require(all(r[k]==v for k,v in identity.items()) and len(r["case_components"])==168
                and tuple((x["kind"],x["slot"]) for x in r["case_components"])==S.OPS
                and all(type(x["slot"]) is int and x["status"]==LEAF_PASS for x in r["case_components"])
                and all(r["source_and_actual_sha256"].get(name)==h for name,h in raw["sha256"].items()),
                "full 168 same-case proof leaves and all 715 observed raw pins")
        native[case]=(folder,r,phase,osr["receipt_sha256"])
    require(tuple(native)==W.FIRST,"all four independent receipts required ahead of any hull leaf")
    return native


def checked_interval(x):
    require(type(x) is list and len(x)==2 and all(type(v) is str for v in x),"explicit string signed interval pair")
    return K.interval(x)


def hull_channel(ch,raw_indices,producer):
    require(type(ch["rows"]) is list and len(ch["rows"])==len(raw_indices)==len(producer),
            "complete independent signed channel")
    results=[]
    misses=[]
    for row,idx,raw in zip(ch["rows"],raw_indices,producer,strict=True):
        i=row["index"]
        valid_type=(type(i) is int if type(idx) is int else
                    type(i) is list and len(i)==2 and all(type(s) is int for s in i))
        require(valid_type and i==idx,"strict same typed local row index")
        a,b=checked_interval(row["independent_exact"])
        c,d=checked_interval(row["producer_exact"])
        enclosure=checked_interval(row["hull_exact"])
        nested=c<=a<=b<=d
        disjoint=max(a,c)>min(b,d)
        require((c,d)==checked_interval(raw) and row["producer_contains_independent"] is nested
                and row["disjoint"] is disjoint and not disjoint and enclosure==(min(a,c),max(b,d)),
                "only exact full hull of valid independent and matching producer intervals; any disjoint row holds")
        if not nested:
            misses.append(idx)
        results.append(list(map(str,enclosure)))
    require(ch["disjoint_indices"]==[] and ch["producer_containment_misses"]==misses,
            "retain exact diagnostics and reject all true disjoint")
    return results,len(misses)


def packet_operands(leaf,meta,prod,par,proof,phase,tail,target):
    kind,slot=meta["kind"],meta["slot"]
    ids=proof["global_indices"]
    require(leaf["status"]==LEAF_PASS and leaf["phase"]==phase and leaf["packet_sha256"]==meta["packet_sha256"]
            and leaf["selector_leaf_content_sha256"]==W.B.digest(proof)
            and (leaf["kind"],leaf["slot"],leaf["global_indices"])==(kind,slot,ids)
            and ids==meta["global_indices"] and leaf["record_id"]==target["record_id"]
            and type(leaf["slot"]) is int and leaf["only_independently_valid_enclosures_used_for_hull_truth"] is True
            and leaf["overlap_not_used_as_proof"] is True and leaf["root_asserted"] is False
            and leaf["checked_cache_sites"]==(64 if kind=="bulk" else 12)
            and (leaf["point_force_scope"],leaf["whole_H_scope"])==(S.POINT,S.BOX)
            and Q(leaf["radius_exact"])==S.RADIUS,
            "admit only native proved own-interval leaf with same cache, chart and selector")
    W.validate_packet(prod,proof,target,slot,kind,tail)
    require(prod["trial_box_sha256"]==target["trial_box_sha256"] and prod["selector_plan_sha256"]==par["plan"]["plan_sha256"]
            and prod["worker_identity"]=={k:par[k] for k in ("pid","supervisor_pid","stage","case_id","part","bindings")},
            "same original observed raw leaf and worker identity")
    F=prod["point_force"]
    H=prod["whole_box_hessian"]
    n=len(ids)
    ij=list(combinations_with_replacement(range(n),2))
    art=leaf["independent_point_and_whole_tails_exact"]
    require(len(art["F"])==n and all(0<=Q(x)<=tail[0] for x in art["F"])
            and leaf["producer_force_tails_exact"]==F["force_tail_upper_exact"],"separate positive F tail arithmetic retained")
    htail=art["H"]
    wh=leaf["independent_same_Y_whole_H_correction_exact"]
    require(type(htail) is list and type(wh) is list
            and all(type(x) is list and len(x)==3 and type(x[0]) is type(x[1]) is int for x in htail)
            and all(type(x) is list and len(x)==4 and type(x[0]) is type(x[1]) is int for x in wh)
            and [tuple(x[:2]) for x in htail]==[tuple(x[:2]) for x in wh]==ij
            and all(0<=Q(x[2])<=tail[1] for x in htail),"separate complete positive H tails and finite cache correction")
    for x in wh:
        checked_interval(x[2:])
    ch=leaf["signed_channels"]
    require(set(ch)==set(CHANNELS),"three independent signed channels required")
    f,nf=hull_channel(ch["point_force"],list(range(n)),F["F_I_minus_F_cache_intervals_exact"])
    h,nh=hull_channel(ch["whole_box_H"],list(map(list,ij)),[x[2:] for x in H["H_I_minus_H_cache_upper_exact"]])
    _,nc=hull_channel(ch["point_same_Y_cache_correction"],list(range(n)),F["Qstar_minus_cache_correction_exact"])
    route={"kind":kind,"slot":slot,"global_indices":ids}
    return route,{**route,"point_force":{"scope":S.POINT,"F_I_minus_F_cache_intervals_exact":f},
        "whole_box_hessian":{"scope":S.BOX,"radius_exact":str(S.RADIUS),
                            "H_I_minus_H_cache_upper_exact":[[i,j,*ep] for (i,j),ep in zip(ij,h,strict=True)]}},(nf,nh,nc)


def case_hulls(a,case,phases,raw_parts):
    folder,r,phase,receipt_h=phases[case]
    opts,d,tail,proofs=W.reviewed_context(a,case)
    plan=read(S.PLAN/case/"plan.json",opts["plan_sha256"])
    require(len(plan["signed_cached_RF_center_rows_exact"])==611 and tuple(map(Q,plan["sizing"]["caps_exact"]))==tail
            and Q(plan["whole_box_radius_exact"])==S.RADIUS and plan["trial_box_sha256"]==opts["target"]["trial_box_sha256"],
            "current independently bound same R cached anchor and exact radius")
    total=Counter()
    routes=[]
    packets=[]
    pins={}
    for meta in r["case_components"]:
        k,i=meta["kind"],meta["slot"]
        rawbase,part=raw_parts[(case,k)]
        raw_meta=part["components"][i]
        fname=f"{k}_{i:03d}.json"
        rawpath=rawbase/fname
        leafpath=folder/fname
        require(meta["packet_file"]==rel(rawpath) and meta["native_file"]==rel(leafpath)
                and meta["packet_sha256"]==raw_meta["sha256"],"strict current canonical proof leaf file and digest")
        leaf=read(leafpath,meta["native_sha256"])
        prod=read(rawpath,meta["packet_sha256"])
        route,p,n=packet_operands(leaf,meta,prod,part,proofs[(k,i)],phase,tail,opts["target"])
        require(route["global_indices"]==next(x["global_indices"] for x in d["entry"]["component_catalogue"]
                if (x["kind"],x["slot"])==(k,i)),"same independently frozen global projection")
        routes.append(route)
        packets.append(p)
        for key,num in zip(CHANNELS,n,strict=True):
            total[key]+=num
        pins[rel(leafpath)]=meta["native_sha256"]
    require(dict(total)==r["producer_containment_miss_counts"],"all nonnested signed scalar diagnostics accounted")
    return {"opts":opts,"d":d,"routes":routes,"packets":packets,"native_leaf_sha256":pins,
            "native_hull_receipt_sha256":receipt_h,"containment_misses":dict(total),"plan":plan}

def optional_spd(spec):
    if spec is None:
        return None
    folder=SPD/"actual_replay"
    ext=read(folder/"os_exit_observation.json",spec["observed_os_exit_sha256"])
    require(ext["status"]=="PASS_ROOT_OBSERVED_REMAINING32_SPD_REPLAY_EXIT" and V.zero(ext["observed_exit_code"])
            and V.positive_int(ext["pid"]) and V.positive_int(ext["unified_exec_session"]),
            "optional SPD minimum requires separate actual observed native exit")
    terminal=read(folder/"terminal.json",ext["terminal_sha256"])
    receipt=read(folder/"receipt.json",ext["receipt_sha256"])
    approved=spec["replay_phase"]
    paths={"source_sha256":SPD/"replay_spd.py","protocol_sha256":SPD/"protocol.json","preflight_sha256":SPD/"preflight.json"}
    require(type(approved) is dict and set(approved)==set(paths) and all(sha(p)==approved[k] for k,p in paths.items()),
            "explicit three-pin reviewed same-R SPD source phase")
    phase=receipt["source_and_phase"]
    require(terminal["status"]=="FINISHED" and V.zero(terminal["exit_code"]) and V.positive_int(terminal["pid"])
            and V.positive_int(phase["pid"]) and phase["pid"]==terminal["pid"]==ext["pid"]
            and terminal["receipt_sha256"]==ext["receipt_sha256"]
            and all(terminal[k]==phase[k]==v for k,v in approved.items())
            and terminal["producer_root_admission_sha256"]==phase["producer_root_admission_sha256"]
                ==ext["producer_root_admission_sha256"]
            and receipt["status"] in ("PASS_INDEPENDENT_ACTUAL_REMAINING32_EXACT_SYMMETRIC_PARTS",
                                       "REMAINING32_SPD_REPLAY_RETAINS_UNRESOLVED_CASES")
            and receipt["root_or_minimum_claimed"] is False and receipt["preconditioner_or_archived_field_changed"] is False,
            "optional independently verified finite same-R congruence premise only")
    rdoc=read(SPD/"root_source_review.json",phase["source_review_sha256"])
    require(rdoc["status"]=="PASS_ROOT_REMAINING32_SPD_REPLAY_SOURCE" and V.yes(rdoc["checks"])
            and all(rdoc[k]==v for k,v in approved.items()),"matching reviewed optional native SPD derivation")
    rows=receipt["case_results"]
    require(tuple(x["case_id"] for x in rows)==V.A.IDS
            and sum(x["status"]==SPD_PASS for x in rows)==receipt["independently_verified_symmetric_parts_SPD"],
            "complete independently reviewed 32 matching SPD cases")
    return {r["case_id"]:r for r in rows},phase,ext["receipt_sha256"]


def min_premise(spd,obj,rnorm):
    if spd is None:
        return None
    case=obj["opts"]["case_id"]
    rows,phase,receipt_h=spd
    row=rows[case]
    path=SPD/"actual_replay"/(case+".json")
    require(row["receipt"]==rel(path),"same independently reviewed optional SPD case")
    leaf=read(path,row["sha256"])
    if leaf["status"]!=SPD_PASS:
        require(row["status"]==leaf["status"],"unresolved SPD remains unresolved")
        return None
    t=obj["opts"]["target"]
    d=obj["d"]
    ctx=d["ctx"]
    cell=ctx.row["parameter_cell_sha256"]
    ident=leaf["target_identity"]
    expect={"pilot_id":case,"record_id":t["record_id"],"field_sha256":t["original_observer_field_sha256"],
            "cache_packet_sha256":t["point_cache_sha256"],"embedded_cache_sha256":ctx.row["frozen_binary_cache_sha256"],
            "R_binary_sha256":t["R_binary_sha256"],"center_binary_sha256":t["center_binary_sha256"],
            "model_contract_sha256":t["model_contract_sha256"],"parameter_cell_sha256":cell,
            "dose_exact":str(ctx.dose),"gain_exact":str(ctx.base["common_gain"]),
            "drive_exact":t["exact_drive"],"energy_scale_binary_exact":t["energy_scale_binary_exact"]}
    require(leaf["replay_phase"]==phase and leaf["strict_congruence_bound"] is True and row["status"]==SPD_PASS
            and leaf["root_or_minimum_claimed"] is False and ident==row["target_identity"]==expect
            and leaf["trial_box_sha256"]==row["trial_box_sha256"]==t["trial_box_sha256"]
            and leaf["R_file_sha256"]==t["point_cache_sha256"] and leaf["R_file"]==ctx.row["point_cache_file"]
            and leaf["native_cached_derivative_case_sha256"]==obj["plan"]["bindings"]["cached_derivative_case_sha256"]
            and leaf["native_geometry_case_sha256"]==obj["plan"]["bindings"]["geometry_replay_case_sha256"]
            and Q(leaf["archived_shift_exact"])==Q(t["archived_center_shift_exact"])
            and Q(leaf["archived_addition_gap_exact"])==Q(t["original_addition_gap_exact"])
            and rnorm==Q(leaf["exact_R_infinity_norm"])==Q(row["exact_R_infinity_norm"]),
            "optional proof is for the identical live binary R and exact target")
    bounds=list(map(Q,leaf["independent_row_absolute_sum_upper_exact"]))
    require(len(bounds)==611 and min(bounds)>=0 and max(bounds)==Q(leaf["independent_q_infinity_upper_exact"])<1,
            "strict positive exact Sym(R) follows only from admitted congruence bound")
    return {"same_Sym_R_SPD_native_case_sha256":row["sha256"],"same_Sym_R_SPD_native_receipt_sha256":receipt_h}


def transfer(arithmetic,target,premise,rnorm):
    sigma,xi,goal=map(Q,(target["archived_center_shift_exact"],target["original_addition_gap_exact"],target["observer_goal_exact"]))
    require(0<=sigma<goal and xi>=0 and rnorm>0,"fixed archived observer and error accounting")
    passed=arithmetic["strict_contraction_and_selfmap"] is True
    posterior=Q(arithmetic["posterior_error_upper_exact"]) if passed else None
    if passed:
        require(0<=posterior<S.RADIUS,"posterior root error inside same certified trial")
    z=Q(arithmetic["z_upper_exact"])
    eligible=passed and sigma+posterior<goal
    # Sym(R)>0 and ||I-RH||<1 imply H>0 for symmetric H:
    # for A=RH and P=H^{-1}, AP+PA^T=2 Sym(R). The stable
    # Lyapunov integral is positive. The inverse infinity bound
    # then supplies lambda_min(H) >= (1-z)/||R||_infinity.
    minimum=passed and premise is not None and 0<=z<1
    return {"archived_center_shift_exact":str(sigma),"archived_addition_gap_exact":str(xi),
            "archived_observer_goal_exact":str(goal),"strict_shift_inclusive_observer_budget_candidate":eligible,
            "shift_plus_root_posterior_exact":str(sigma+posterior) if passed else None,
            "physical_nodal_error_um_upper_exact":str(Fraction(67,2)*(sigma+posterior+xi)) if passed else None,
            "independent_observer_transfer_still_required":True,
            "same_target_SPD_premise":premise,
            "candidate_strict_minimum_requires_independent_final_review":minimum,
            "candidate_whole_box_H_euclidean_lower_exact":str((1-z)/rnorm) if minimum else None}


def reduce_case(obj,spd,folder,mark):
    ctx=obj["d"]["ctx"]
    p=obj["plan"]
    t=obj["opts"]["target"]
    K.ctx.prec=640
    mat=[[V.A.T.rational(x) for x in row] for row in ctx.R]
    norm=max(sum(abs(x) for x in row) for row in mat)
    z=Q(p["sizing"]["ZQ_exact"])
    require(norm==Q(p["sizing"]["R_inf_exact"])>0 and V.A.a_sha(ctx.R)==t["R_binary_sha256"]
            and ctx.row["point_cache_sha256"]==t["point_cache_sha256"]
            and V.A.a_sha(ctx.u)==t["center_binary_sha256"] and ctx.drive==Q(t["exact_drive"]),"one exact R and signed anchor")
    anchor=[checked_interval(x) for x in p["signed_cached_RF_center_rows_exact"]]
    maps=obj["routes"]
    force,H=K.assemble(obj["packets"],maps,611)
    prem=min_premise(spd,obj,norm)
    fields={"status":"ADMITTED_ALL_INDEPENDENT_SIGNED_HULLS","target":t,"radius_exact":str(S.RADIUS),
        "exact_material_singleton":ctx.row["parameter_cell"],"dose_exact":str(ctx.dose),"exact_R_infinity_norm":str(norm),
        "ZQ_whole_cube_upper_exact":str(z),"cached_anchor_already_preconditioned":True,
        "signed_native_R_F_Q_exact":[list(map(str,x)) for x in anchor],
        "F_I_minus_F_Q_hull_global_exact":[list(map(str,x)) for x in force],
        "H_I_minus_H_Q_hull_global_upper_exact":[[i,j,*map(str,x)] for (i,j),x in sorted(H.items())],
        "maps":maps,"whole_case_hull_receipt_sha256":obj["native_hull_receipt_sha256"],
        "native_component_sha256":obj["native_leaf_sha256"],"producer_containment_miss_counts":obj["containment_misses"],
        "selector_and_cached_geometry_bindings":p["bindings"],"root_asserted":False,**mark}
    oh=save(folder/"assembled_operands.json",fields)
    arithmetic=K.contract(mat,force,H,anchor,z,radius=S.RADIUS)
    assess=transfer(arithmetic,t,prem,norm)
    result={"status":"FIRST_WAVE_SPATIAL_REDUCTION_AWAITS_INDEPENDENT_GLOBAL_REVIEW",
            "scope":SCOPE,"target":t,"assembled_operands_sha256":oh,"arithmetic":arithmetic,
            "observer_and_optional_minimum":assess,"independent_global_endpoint_replay_required":True,
            "exact_spatial_root_certified":False,"exact_spatial_minimum_certified":False,
            "new_spatial_integrals":0,"new_AD_H_moments_or_solves":0,"root_asserted":False,**mark}
    h=save(folder/"results.json",result)
    return {"case_id":t["case_id"],"results_sha256":h,"assembled_operands_sha256":oh,
            "classification":arithmetic["classification"],"observer_budget_candidate":assess["strict_shift_inclusive_observer_budget_candidate"],
            "minimum_candidate":assess["candidate_strict_minimum_requires_independent_final_review"]}


def run(args):
    a,w,j,b=gate(args)
    raw,parts=observed_raw(j["root_producer_admission_sha256"],a,w)
    phases=open_all_hull_metadata(a,j,raw)
    spd=optional_spd(j["optional_spd"])
    ACTUAL.mkdir(parents=True,exist_ok=False)
    mark={"pid":os.getpid(),"case_order":list(W.FIRST),"bindings":b}
    save(ACTUAL/"started.json",{"status":"RUNNING",**mark})
    result={"status":"FAILED","exit_code":1,**mark}
    rows=[]
    try:
        inputs=[case_hulls(a,case,phases,parts) for case in W.FIRST]
        # Every individual endpoint/provenance check precedes any matrix reduction.
        for case,obj in zip(W.FIRST,inputs,strict=True):
            folder=ACTUAL/case
            folder.mkdir(exist_ok=False)
            row=reduce_case(obj,spd,folder,mark)
            rows.append(row)
            print(row,flush=True)
        h=save(ACTUAL/"results.json",{"status":"COMPLETE_FIRST_WAVE_SIGNED_REDUCTION_AWAITS_INDEPENDENT_GLOBAL_REVIEW",
          **mark,"cases":rows,"root_asserted":False,"full_five_gate_complete":False})
        result.update(status="FINISHED",exit_code=0,results_sha256=h)
    except BaseException as e:
        save(ACTUAL/"failure.json",{**mark,"error":str(e),"traceback":traceback.format_exc(),"closed_candidate_cases":rows})
        raise
    finally:
        save(ACTUAL/"terminal.json",result)


def fixtures():
    c={"rows":[{"index":0,"independent_exact":["-1","1"],"producer_exact":["0","2"],
                "hull_exact":["-1","2"],"producer_contains_independent":False,"disjoint":False}],
       "producer_containment_misses":[0],"disjoint_indices":[]}
    require(hull_channel(c,[0],[["0","2"]])==([["-1","2"]],1),"exact nonnested safe hull")
    variants=[{"hull_exact":["0","1"]},{"hull_exact":["0","2"]},{"producer_contains_independent":True},
              {"index":False},{"independent_exact":["-3","-2"],"disjoint":True,"hull_exact":["-3","2"]}]
    fail=all(V.A.reject(lambda:hull_channel({**c,"rows":[{**c["rows"][0],**v}]},[0],[["0","2"]])) for v in variants)
    K.ctx.prec=640
    raw=K.contract([[Fraction(2)]],[(Fraction(-1,8),Fraction(-1,8))],{(0,0):(Fraction(1,8),Fraction(1,8))},
       [(Fraction(1,4),Fraction(1,4))],Fraction(1,4),S.RADIUS)
    require(Q(raw["eta_I_upper_exact"])==0 and Q(raw["z_upper_exact"])==Fraction(1,2)
            and raw["strict_contraction_and_selfmap"] is True,"cached preconditioned RF anchor added once")
    target={"archived_center_shift_exact":"1/8","original_addition_gap_exact":"1/16","observer_goal_exact":"1/4"}
    tr=transfer(raw,target,{"synthetic_Sym_R":True},Fraction(2))
    require(tr["strict_shift_inclusive_observer_budget_candidate"] is True
            and Q(tr["physical_nodal_error_um_upper_exact"])==Fraction(201,32)
            and Q(tr["candidate_whole_box_H_euclidean_lower_exact"])==Fraction(1,4)
            and not transfer(raw,target,None,Fraction(2))["candidate_strict_minimum_requires_independent_final_review"],
            "exact archived origin and optional SPD classification")
    return {"non_nested_full_union_admitted_and_shrunk_intersection_disjoint_fake_index_rejected":fail,
            "already_R_weighted_anchor_added_once_before_norm":True,
            "strict_sigma_budget_exact_original_xi_and_optional_SPD_only":True,
            "cold_no_live_component_or_reduction_payloads_opened":True}


def prepare():
    a,w=static_phase()
    checks=fixtures()
    require(V.yes(checks) and not ACTUAL.exists(),"cold fixtures and actual output absence required")
    PREP.mkdir(parents=True,exist_ok=False)
    source=sha(SOURCE)
    p={"status":"FROZEN_FIRST_WAVE_HULL_ONLY_SIGNED_REDUCER","scope":SCOPE,"source_sha256":source,
       **ZERO,"first_wave":list(W.FIRST),"worker_runtime_binding":w,"worker_target_authority":a["targets"],
       "fixed_source_and_input_sha256":{rel(k):v for k,v in PINS.items()},
       "radius_exact":str(S.RADIUS),"analytic_point_force_scope":S.POINT,"analytic_whole_H_scope":S.BOX,
       "input_admission_path":rel(ADMIT_FILE),"independent_reducer_source_review_required":True,
       "independent_final_global_endpoint_replay_required":True,"no_later_waves_authorized":True,
       "optional_SPD_proof_producer_scope":"Exact Sym(R) only; minimum conditional on strict contraction and final replay.",
       "full_five_gate_complete":False}
    ph=save(PREP/"protocol.json",p)
    fh=save(PREP/"preflight.json",{"status":"PASS","source_sha256":source,"protocol_sha256":ph,**ZERO,"checks":checks})
    with (PREP/"reducer_frozen.py").open("xb") as f:
        f.write(SOURCE.read_bytes())
    print({"source_sha256":source,"protocol_sha256":ph,"preflight_sha256":fh},flush=True)


if __name__=="__main__":
    cli=argparse.ArgumentParser(description=__doc__)
    sub=cli.add_subparsers(dest="phase",required=True)
    sub.add_parser("prepare")
    actual=sub.add_parser("run")
    for key in ("protocol","preflight","review","admission"):
        actual.add_argument("--"+key+"-sha256",required=True)
        if key in ("review","admission"):
            actual.add_argument("--"+key,type=lambda x:Path(x).resolve(),required=True)
    args=cli.parse_args()
    prepare() if args.phase=="prepare" else run(args)
