"""Reduce one fixed singleton only after verified raw v1 plus endpoint companion v3."""
import argparse
from copy import deepcopy
from fractions import Fraction as Q
from hashlib import sha256
import json
import os
from pathlib import Path
import traceback

import transfer_ambiguity_exact_bridge_reduce_driver as base

SOURCE=Path(__file__).resolve()
ROOT=base.ROOT
OUT=base.OUT
PREP=OUT/"caller_preparation_v3"
ACTUAL=OUT/"actual_v3"
COMPANION=base.REPLAY/"exact_endpoint_companion"
V3=COMPANION/"actual_v3"
RAW=base.REPLAY/"actual_replay"
SOURCE_PINS={
    base.SOURCE:"6eafff943ae116d4f9550f480f5a3c6e75648e5428291aa263238a762555a998",
    base.PREP/"protocol.json":"bdb7d8576aad46407a2a6caef8dfc40b79b54034feb0c01803d77e470df94a40",
    base.PREP/"preflight.json":"7b61454b090aa18cfb4fcf6188a0a3b79047b8ea323e3b678d7dd8c67cce8213",
    base.PREP/"independent_prelaunch/review.json":"02df786806f0950cbab78a8303d9d6a3acd4236cf34ecbbd9293d667bc4c9aa2",
    COMPANION/"tighten_v3.py":"d9ca0122c62a6bd5f94e8bdba078c715adca9056d8b8822e693e2d4a866419de",
    COMPANION/"protocol_v3.json":"a825830947d53daee36a9f75904a85d2097375fadcf070fcd3ee30d1ffb7c9d7",
    COMPANION/"preflight_v3.json":"a62668937b49f9e76e65706d70f2ba015820923150183396815ff26c82a3c01d",
    COMPANION/"independent_review_v3/receipt.json":"3e91cfe40c26688466b0a3e05cca7716e87f53016afce37d06e17f304b0847c9",
}
PHASE={
    "source_sha256":SOURCE_PINS[COMPANION/"tighten_v3.py"],
    "protocol_sha256":SOURCE_PINS[COMPANION/"protocol_v3.json"],
    "preflight_sha256":SOURCE_PINS[COMPANION/"preflight_v3.json"],
}
REVIEW_STATUS="PASS_INDEPENDENT_TWO_PROOF_SINGLETON_REDUCTION_CALLER_PRELAUNCH"
RAW_OS="OBSERVED_RAW_BRIDGE_REPLAY_SHELL_TERMINAL"
FRESH_OS="PASS_OBSERVED_SURFACE_H_ENDPOINT_COMPANION_TRUE_OS_EXIT"
V3_CHECKS=base.NATIVE_CHECKS|{"exact_tail_inflation_erosion_replayed_and_old_result_preserved",
                             "fresh_all_24_surface_H_validated_with_signed_containment"}
FKEY="F_I_minus_F_cache_intervals_exact"
HKEY="H_I_minus_H_cache_upper_exact"
SCOPE=("One fixed exact-spatial finite-element target and radius-2^-34 cube, using explicit full raw v1 "
       "bulk-H and direct-cache proofs plus v3 refined F and fresh 24 surface H proofs. "
       "Signed matrix contraction is an unreviewed arithmetic candidate. No minimum, "
       "continuum, loading, source or cross-target claim.")
require=base.require
sha,pinned,relative,save=base.sha,base.pinned,base.rel,base.save
rat=base.rational


def protocol():
    for path,digest in SOURCE_PINS.items():
        require(sha(path)==digest,f"changed two-proof authority: {relative(path)}")
    source=base.authority_protocol()
    old=pinned(base.PREP/"protocol.json",SOURCE_PINS[base.PREP/"protocol.json"])
    require(source==old,"original reviewed caller prerequisites must be current")
    old_review=base.PREP/"independent_prelaunch/review.json"
    base.review_gate(old_review,SOURCE_PINS[old_review])
    peer=COMPANION/"independent_review_v3/receipt.json"
    p=pinned(peer,SOURCE_PINS[peer])
    require(p["status"]=="PASS_INDEPENDENT_COLD_SURFACE_H_ENDPOINT_COMPANION" and base.yes(p["checks"])
            and all(p[k]==v for k,v in PHASE.items()),"current v3 cold peer proof required")
    frozen=pinned(COMPANION/"protocol_v3.json",PHASE["protocol_sha256"])
    pre=pinned(COMPANION/"preflight_v3.json",PHASE["preflight_sha256"])
    require(frozen["status"]=="FROZEN_EXACT_ENDPOINT_AND_FRESH_SURFACE_H_REVIEW"
            and frozen["source_sha256"]==PHASE["source_sha256"] and frozen["precision_bits"]==base.PRECISION
            and rat(frozen["radius_exact"])==base.RADIUS and pre["status"]=="PASS"
            and pre["source_sha256"]==PHASE["source_sha256"] and pre["protocol_sha256"]==PHASE["protocol_sha256"]
            and pre["actual_raw_component_packets_opened"]==pre["force_H_cache_or_moment_work"]==0,
            "frozen unchanged v3 source, preflight and radius required")
    return {"status":"FROZEN_TWO_PROOF_SINGLETON_BRIDGE_REDUCTION_CALLER",
            "source_sha256":sha(SOURCE),"scope":SCOPE,"base_authority":source,
            "v3_companion_phase":PHASE,"source_pins":{relative(p):h for p,h in SOURCE_PINS.items()},
            "pilot_choices":list(base.PILOTS),"radius_exact":str(base.RADIUS),"precision_bits":base.PRECISION,
            "packet_counts":{"bulk":144,"surface":24},"fresh_surface_H_count":24,
            "raw_OS_receipt_status":RAW_OS,"v3_OS_receipt_status":FRESH_OS,
            "required_external_admission_sha256s":["caller_prelaunch_review","original_root_producer_admission",
                "observed_raw_v1_shell_exit","observed_companion_v3_exit"],
            "actual_packet_payloads_opened":0,"actual_matrix_work":0}


def gate(path,digest):
    current={"source_sha256":sha(SOURCE),"protocol_sha256":sha(PREP/"protocol.json"),
             "preflight_sha256":sha(PREP/"preflight.json"),
             "old_caller_sha256":SOURCE_PINS[base.SOURCE],
             "v3_cold_peer_review_sha256":SOURCE_PINS[COMPANION/"independent_review_v3/receipt.json"],
             "v3_source_sha256":PHASE["source_sha256"],"v3_protocol_sha256":PHASE["protocol_sha256"],
             "v3_preflight_sha256":PHASE["preflight_sha256"]}
    r=pinned(path,digest)
    require(r["status"]==REVIEW_STATUS and base.yes(r["checks"]) and all(r[k]==v for k,v in current.items()),
            "current independent two-proof caller review required")
    p=pinned(PREP/"protocol.json",current["protocol_sha256"])
    pf=pinned(PREP/"preflight.json",current["preflight_sha256"])
    require(p==protocol() and pf["status"]=="PASS" and base.yes(pf["checks"])
            and all(pf[k]==current[k] for k in ("source_sha256","protocol_sha256"))
            and base.true_zero(pf["actual_packet_payloads_opened"]) and base.true_zero(pf["actual_matrix_work"])
            and sha(PREP/"transfer_ambiguity_exact_bridge_reduce_driver_v3_frozen.py")==current["source_sha256"],
            "unmodified cold two-proof phase required")
    return p,{**current,"prelaunch_review_sha256":digest}


def raw_status(obs,native,end,start,pilot,root_sha):
    facts=(end["status"],native["status"],end["exit_code"],obs["observed_exit_code"])
    require(obs["status"]==RAW_OS and base.positive_int(obs["pid"]) and base.positive_int(obs["unified_exec_session"])
            and type(end["exit_code"]) is type(obs["observed_exit_code"]) is type(obs["native_documented_exit_code"]) is int
            and obs["native_documented_exit_code"]==end["exit_code"] and obs["terminal_status"]==end["status"]
            and facts in {("FINISHED",base.NATIVE_STATUS,0,0),
                 ("FINISHED_WITH_HOLD","HOLD_INDEPENDENT_COMPONENT_CONTAINMENT",2,1),
                 ("FINISHED_WITH_HOLD","HOLD_INDEPENDENT_COMPONENT_CONTAINMENT",2,2)},
            "native declared exit and true enclosing-host observation must remain separate")
    phase={"pid":obs["pid"],"root_observation_sha256":root_sha,
           "source_sha256":base.PINS[base.REPLAY_SOURCE],"protocol_sha256":base.PINS[base.REPLAY_PROTOCOL],
           "preflight_sha256":base.PINS[base.REPLAY_PREFLIGHT]}
    require(all(all(x[k]==v for k,v in phase.items()) for x in (native,end,start))
            and native["pilot_id"]==pilot and end["receipt_sha256"]==obs["receipt_sha256"] and start["status"]=="RUNNING"
            and set(native["checks"])==base.NATIVE_CHECKS
            and all(native["checks"][k] is True for k in (
                "true_os_complete_producer_and_both_child_terminals","same_exact_point_and_whole_box_all_168_bindings"))
            and all(native["checks"][k] is (facts[-1]==0) for k in (
                "all_point_cache_force_reconstructed_and_contained","all_point_force_and_whole_H_contained"))
            and rat(native["radius_exact"])==base.RADIUS,
            "raw replay closure must match the exact signed native phase and target")


def v3_status(obs,native,end,start,pilot,raw_obs_sha,raw_obs):
    require(obs["status"]==FRESH_OS and obs["pilot_id"]==native["pilot_id"]==pilot
            and base.positive_int(obs["pid"]) and base.positive_int(obs["unified_exec_session"])
            and base.true_zero(obs["observed_exit_code"]) and end["status"]=="FINISHED"
            and base.true_zero(end["exit_code"]) and end["receipt_sha256"]==obs["receipt_sha256"]
            and base.digest_text(obs["launch_command_sha256"])
            and obs["peer_review_sha256"]==SOURCE_PINS[COMPANION/"independent_review_v3/receipt.json"],
            "genuine zero exit for this v3 companion required")
    peer=COMPANION/"independent_review_v3/receipt.json"
    phase={"pid":obs["pid"],**PHASE,"raw_component_os_exit_sha256":raw_obs_sha,
           "raw_component_receipt_sha256":raw_obs["receipt_sha256"],
           "raw_component_terminal_sha256":raw_obs["terminal_sha256"],
           "peer_review_file":relative(peer),"peer_review_sha256":SOURCE_PINS[peer]}
    require(all(all(x[k]==v for k,v in phase.items()) for x in (native,end,start))
            and start["status"]=="RUNNING" and native["status"]==base.NATIVE_STATUS
            and set(native["checks"])==V3_CHECKS and base.yes(native["checks"])
            and native["remaining_noncontained"]==[]
            and type(native["fresh_surface_H_jobs"]) is int and native["fresh_surface_H_jobs"]==24
            and base.true_zero(native["new_bulk_force_or_Hessian_moment_jobs"])
            and native["original_raw_component_receipt_file"]==relative(RAW/pilot/"receipt.json")
            and rat(native["radius_exact"])==base.RADIUS,
            "same complete approved v3 phase with 24 freshly replayed surface Hessians required")


def phases(args):
    pilot=args.pilot_id
    raw=RAW/pilot
    fresh=V3/pilot
    require(args.raw_observation==raw/"os_exit_observation.json",
            "canonical sealed raw OS observation path required")
    o0=pinned(args.raw_observation,args.raw_observation_sha256)
    require(args.companion_observation.is_relative_to(ROOT),"companion observations must be inside this study")
    o1=pinned(args.companion_observation,args.companion_observation_sha256)
    r0=pinned(raw/"receipt.json",o0["receipt_sha256"])
    t0=pinned(raw/"terminal.json",o0["terminal_sha256"])
    s0h=sha(raw/"started.json")
    s0=pinned(raw/"started.json",s0h)
    raw_status(o0,r0,t0,s0,pilot,args.admission_sha256)
    r1=pinned(fresh/"receipt.json",o1["receipt_sha256"])
    t1=pinned(fresh/"terminal.json",o1["terminal_sha256"])
    s1h=sha(fresh/"started.json")
    s1=pinned(fresh/"started.json",s1h)
    v3_status(o1,r1,t1,s1,pilot,args.raw_observation_sha256,o0)
    same=("record_id","producer_source_sha256","producer_manifest_sha256","producer_terminal_sha256",
          "producer_observed_os_exit_code","producer_observed_pid","producer_observed_exec_session",
          "plan_sha256","plan_review_sha256","derivative_admission_sha256","box_sha256","parameter_cell_sha256",
          "cache_packet_sha256","binary_preconditioner_sha256","signed_cached_anchor_file","signed_cached_anchor_sha256",
          "reviewed_cached_ZQ_upper_exact","radius_exact","checked_packet_sha256","plan_packet_sha256",
          "root_observation_sha256")
    require(all(r0[k]==r1[k] for k in same),"v3 and original replay must attest the same exact physical proof")
    return r0,r1,o0,o1,s0h,s1h


def adopt_catalogue(pilot,r0,r1,plan,manifest,raw_packets,proofs):
    old=base.catalogues(pilot,r0,plan,manifest,raw_packets,proofs)
    base.check_order(r1["components"])
    new=[]
    for (raw,prod,proof,prior),row in zip(old,r1["components"],strict=True):
        path=V3/pilot/prior.name
        require(all(row[k]==v for k,v in raw.items() if k not in ("replay_file","replay_sha256"))
                and row["replay_file"]==relative(path) and base.digest_text(row["replay_sha256"]),
                "noncanonical v3 path or changed raw/producer/route binding")
        base.typed_map(row["global_indices"])
        new.append((raw,row,prod,proof,prior,path))
    require(all(sha(old_path)==raw["replay_sha256"] and sha(new_path)==new["replay_sha256"]
                for raw,new,_,_,old_path,new_path in new),
            "168 native plus 168 companion leaves must be authenticated before parsing any")
    return new


def signed_misses(outer,inner):
    require(len(outer)==len(inner),"matching complete signed arrays required")
    return [i for i, (a,b) in enumerate(zip(outer,inner,strict=True)) if not (a[0]<=b[0]<=b[1]<=a[1])]


def refined_endpoints(old,s,t0,t,limit,given,kind):
    l,u=old
    require(s.denominator.bit_count()==1 and 0<=t<=t0<=s
            and t0<=given<=limit and l+s<=u-s,
            "actual inserted old radius witness and nonnegative valid analytic tail required")
    if kind=="bulk":
        require(t==t0==given,"bulk uses exactly the old mathematically justified tail")
    return l+s-t,u-s+t


def leaf_refinement(raw_row,new_row,raw,fresh,producer,tails):
    common=("kind","slot","global_indices","point_force_scope","whole_box_H_scope","plan_file_sha256",
            "packet_file","packet_sha256","checked_force_rows","checked_upper_H_entries","checked_cache_sites")
    require(all(raw[k]==fresh[k] for k in common)
            and all(fresh[k]==raw_row[k]==new_row[k] for k in common[:8])
            and type(raw["slot"]) is type(fresh["slot"]) is int
            and fresh["raw_component_replay_file"]==raw_row["replay_file"]
            and fresh["raw_component_replay_sha256"]==raw_row["replay_sha256"],
            "companion leaf must directly bind its same typed raw proof")
    ids=base.typed_map(raw["global_indices"])
    m=len(ids)
    count=m*(m+1)//2
    sites=64 if raw["kind"]=="bulk" else 12
    for k,n in (("checked_force_rows",m),("checked_upper_H_entries",count),("checked_cache_sites",sites)):
        require(type(raw[k]) is type(fresh[k]) is int and raw[k]==fresh[k]==n,"wrong local count")
    bad=raw["containment_failures"]
    repaired=fresh["containment_failures"]
    require(set(bad)==set(repaired)==base.FAILURE_KEYS and bad["direct_point_cache"]==[]
            and all(v==[] for v in repaired.values()),"original direct cache and new all-channel strict containment required")
    bF=[base.interval(x) for x in raw["independent_point_F_exact"]]
    nF=[base.interval(x) for x in fresh["independent_point_F_exact"]]
    pF=[base.interval(x) for x in producer["point_force"][FKEY]]
    qtail=list(map(rat,producer["point_force"]["force_tail_upper_exact"]))
    require(len(bF)==len(nF)==len(pF)==len(qtail)==m and all(0<=x<=tails[0] for x in qtail),
            "complete point rows and declared producer tail range required")
    ids_bad=signed_misses(pF,bF)
    require(type(bad["point_force"]) is list and all(type(i) is int for i in bad["point_force"])
            and bad["point_force"]==ids_bad,"raw F miss indices must be exact")
    witness=fresh["point_original_excess_rounded_tails"]
    require(type(witness) is list and all(type(row) is list and len(row)==5 and type(row[0]) is int and row[4] is True
                                          for row in witness)
            and [k for k,*_ in witness]==ids_bad,"one exact refinement witness for each and only each raw point miss")
    corrected=bF.copy()
    for j,s,t0,t,_ in witness:
        corrected[j]=refined_endpoints(bF[j],rat(s),rat(t0),rat(t),tails[0],qtail[j],raw["kind"])
    require(corrected==nF and not signed_misses(pF,nF),
            "no fabricated narrowing, and no overlap substituted for point containment")
    oldH=base.upper_intervals(raw["independent_whole_box_H_exact"],m)
    newH=base.upper_intervals(fresh["independent_whole_box_H_exact"],m)
    prodH=base.upper_intervals(producer["whole_box_hessian"][HKEY],m)
    order=[row[:2] for row in fresh["independent_whole_box_H_exact"]]
    old_bad=[order[i] for i in signed_misses(prodH,oldH)]
    require(type(bad["whole_box_H"]) is list and all(type(i) is list and len(i)==2
            and all(type(k) is int for k in i) for i in bad["whole_box_H"])
            and bad["whole_box_H"]==old_bad and not signed_misses(prodH,newH),
            "old H miss witness and new exact H containment required")
    new_tail=fresh["fresh_surface_H_tails"]
    if raw["kind"]=="bulk":
        require(old_bad==new_tail==[] and raw["independent_whole_box_H_exact"]==fresh["independent_whole_box_H_exact"],
                "bulk H must be the immutable already-contained raw proof")
    else:
        require(type(new_tail) is list and all(type(row) is list and len(row)==3 and type(row[0]) is type(row[1]) is int
                                               and 0<=rat(row[2])<=tails[1] for row in new_tail)
                and [row[:2] for row in new_tail]==order,"all exact fresh surface Hessian tail entries required")
    return len(witness),len(old_bad)


def admitted(args,bindings):
    r0,r1,o0,o1,s0h,s1h=phases(args)
    review=base.loader(base.REPLAY_SOURCE,"_singleton_bridge_v1_original_admission")
    a,plan,links,m,packets,proofs=review.admitted(args.admission,args.admission_sha256,args.pilot_id)
    run=base.BRIDGE/"actual"/args.pilot_id
    require(a["producer_protocol_file"]==relative(base.BRIDGE_PROTOCOL)
            and a["producer_preflight_file"]==relative(base.BRIDGE_PREFLIGHT)
            and a["producer_frozen_source_file"]==relative(base.BRIDGE/"preparation_v2/transfer_ambiguity_exact_bridge_frozen.py")
            and links["source_sha256"]==r0["producer_source_sha256"]==base.PINS[base.BRIDGE_SOURCE]
            and links["protocol_sha256"]==base.PINS[base.BRIDGE_PROTOCOL]
            and links["preflight_sha256"]==base.PINS[base.BRIDGE_PREFLIGHT]
            and links["derivative_admission_sha256"]==r0["derivative_admission_sha256"]==base.PINS[base.DERIV_ADMISSION]
            and links["plan_sha256"]==r0["plan_sha256"] and links["plan_review_sha256"]==r0["plan_review_sha256"]
            and base.positive_int(a["supervisor_pid"]) and base.positive_int(r0["producer_observed_pid"])
            and base.positive_int(r0["producer_observed_exec_session"])
            and a["supervisor_pid"]==r0["producer_observed_pid"] and a["unified_exec_session"]==r0["producer_observed_exec_session"]
            and base.true_zero(r0["producer_observed_os_exit_code"])
            and a["sha256"][relative(run/"manifest.json")]==r0["producer_manifest_sha256"]
            and a["sha256"][relative(run/"terminal.json")]==r0["producer_terminal_sha256"],
            "same admitted producer, current plan and genuine worker process exits required")
    end=pinned(run/"terminal.json",r0["producer_terminal_sha256"])
    require(base.positive_int(end["pid"]) and end["pid"]==a["supervisor_pid"]
            and all(base.positive_int(w["pid"]) and base.positive_int(w["supervisor_pid"])
                and w["supervisor_pid"]==end["pid"] and base.true_zero(w["observed_exit_code"]) for w in end["workers"]),
            "literal producer supervisor and child identities/exits")
    entries=adopt_catalogue(args.pilot_id,r0,r1,plan,m,packets,proofs)
    actual_links={**bindings,"producer_admission_file":relative(args.admission),"producer_admission_sha256":args.admission_sha256,
        "raw_os_file":relative(args.raw_observation),"raw_os_sha256":args.raw_observation_sha256,
        "raw_native_receipt_sha256":o0["receipt_sha256"],"raw_native_terminal_sha256":o0["terminal_sha256"],
        "raw_native_started_sha256":s0h,"raw_host_observed_code":o0["observed_exit_code"],
        "raw_native_declared_exit_code":o0["native_documented_exit_code"],
        "companion_os_file":relative(args.companion_observation),"companion_os_sha256":args.companion_observation_sha256,
        "companion_receipt_sha256":o1["receipt_sha256"],"companion_terminal_sha256":o1["terminal_sha256"],
        "companion_started_sha256":s1h,"companion_host_observed_code":o1["observed_exit_code"],
        "plan_sha256":r0["plan_sha256"],"producer_manifest_sha256":r0["producer_manifest_sha256"]}
    return r0,r1,plan,links,entries,actual_links


def run(args):
    prot,links=gate(args.review,args.review_sha256)
    r0,r1,plan,producer_links,entries,links=admitted(args,links)
    folder=ACTUAL/args.pilot_id
    folder.mkdir(parents=True,exist_ok=False)
    mark={"pid":os.getpid(),"pilot_id":args.pilot_id,"bindings":links}
    save(folder/"started.json",{"status":"RUNNING",**mark})
    end={"status":"FAILED","exit_code":1,**mark}
    try:
        ctx,case,entry,routes,anchor,z,digest,eta_q=base.context(args.pilot_id,r0,plan,producer_links,prot["base_authority"])
        tails=tuple(map(rat,plan["tail_caps_exact"]))
        require(len(tails)==2 and min(tails)>0,"matched exact positive tail caps")
        maps,packets=[],[]
        fixed,h_misses=0,0
        for raw,new,prod,proof,old_path,new_path in entries:
            p=pinned(prod,raw["packet_sha256"])
            sel=pinned(proof,raw["plan_file_sha256"])
            old=pinned(old_path,raw["replay_sha256"])
            fresh=pinned(new_path,new["replay_sha256"])
            nf,nh=leaf_refinement(raw,new,old,fresh,p,tails)
            fixed+=nf
            h_misses+=nh
            k,i=raw["kind"],raw["slot"]
            route,packet=base.component_operands(raw,fresh,p,sel,routes[k][i],ctx.record_id,producer_links,digest(sel))
            maps.append(route)
            packets.append(packet)
        require(r1["exact_endpoint_refinement_counts"]=={"point_force":fixed,"surface_H_recomputations":24},
                "all refined point rows and fresh surface H recertifications accounted for")
        kernel=base.loader(base.KERNEL_SOURCE,"_same_reviewed_spatial_signed_reduction")
        kernel.ctx.prec=base.PRECISION
        require(kernel.RADIUS==base.RADIUS and kernel.POINT==base.POINT and kernel.BOX==base.BOX,"kernel scope")
        force,H=kernel.assemble(packets,maps,base.N)
        R=[list(map(base.binary_fraction,row)) for row in ctx.R.tolist()]
        payload={"status":"ADMITTED_TWO_PROOF_SIGNED_SPATIAL_OPERANDS","pilot_id":args.pilot_id,
                 "record_id":ctx.record_id,"scope":SCOPE,"bindings":links,"radius_exact":str(base.RADIUS),
                 "ZQ_upper_exact":str(z),"closed_cell_exact":ctx.row["parameter_cell"],
                 "parameter_cell_sha256":ctx.row["parameter_cell_sha256"],"same_exact_drive":str(ctx.drive),
                 "dose_exact":str(ctx.dose),"binary_R_sha256":case["preconditioner_binary_sha256"],
                 "binary_R_source":case["cache_packet"],"cache_packet_sha256":case["cache_packet_sha256"],
                 "model_contract_sha256":entry["model_contract_sha256"],"geometry_box_sha256":case["box_sha256"],
                 "cached_signed_anchor_sha256":case["verified_point_sha256"],"cached_eta_upper_exact":str(eta_q),
                 "signed_R_F_cache_exact":[list(map(str,row)) for row in anchor],
                 "signed_F_defect_exact":[list(map(str,row)) for row in force],
                 "whole_box_signed_H_defect_upper_exact":[[i,j,*map(str,v)] for (i,j),v in sorted(H.items())],
                 "maps":maps,"packet_sha256":r0["checked_packet_sha256"],
                 "raw_native_component_sha256":{e["replay_file"]:e["replay_sha256"] for e,*_ in entries},
                 "companion_component_sha256":{e["replay_file"]:e["replay_sha256"] for _,e,*_ in entries},
                 "old_noncontained_F_rows_refined":fixed,"old_surface_H_misses_replaced":h_misses}
        oh=save(folder/"assembled_operands.json",payload)
        arithmetic=kernel.contract(R,force,H,anchor,z,radius=base.RADIUS)
        goal=rat(case["observer_goal_exact"])
        posterior=arithmetic["posterior_error_upper_exact"]
        rh=save(folder/"results.json",{"status":"SPATIAL_REDUCTION_CANDIDATE_AWAITING_INDEPENDENT_REVIEW",
             "scope":SCOPE,**mark,"assembled_operands_sha256":oh,"arithmetic":arithmetic,
             "observer_design_goal_exact":str(goal),"posterior_below_observer_design_goal":posterior is not None and rat(posterior)<goal,
             "independent_actual_reduction_review_required":True,"exact_spatial_root_certified":False,
             "force_H_AD_or_equilibrium_jobs_run":0})
        end.update(status="FINISHED",exit_code=0,results_sha256=rh,assembled_operands_sha256=oh)
        print(json.dumps({"result_sha256":rh,"classification":arithmetic["classification"],"pilot_id":args.pilot_id}))
    except BaseException as error:
        save(folder/"failure.json",{"error":str(error),"traceback":traceback.format_exc()})
        raise
    finally:
        save(folder/"terminal.json",end)


def fixtures():
    h="a"*64
    p=base.PILOTS[0]
    base_phase={"pid":7,"root_observation_sha256":h,"source_sha256":base.PINS[base.REPLAY_SOURCE],
                "protocol_sha256":base.PINS[base.REPLAY_PROTOCOL],"preflight_sha256":base.PINS[base.REPLAY_PREFLIGHT]}
    o={"status":RAW_OS,"pid":7,"unified_exec_session":23,"observed_exit_code":1,"native_documented_exit_code":2,
       "terminal_status":"FINISHED_WITH_HOLD","receipt_sha256":h,"terminal_sha256":h}
    n={**base_phase,"pilot_id":p,"status":"HOLD_INDEPENDENT_COMPONENT_CONTAINMENT","radius_exact":str(base.RADIUS),
       "checks":{k:k in {"true_os_complete_producer_and_both_child_terminals","same_exact_point_and_whole_box_all_168_bindings"}
                 for k in base.NATIVE_CHECKS}}
    t={**base_phase,"status":"FINISHED_WITH_HOLD","exit_code":2,"receipt_sha256":h}
    s={**base_phase,"status":"RUNNING"}
    raw_status(o,n,t,s,p,h)
    peer=COMPANION/"independent_review_v3/receipt.json"
    phase={**PHASE,"pid":8,"raw_component_receipt_sha256":h,"raw_component_terminal_sha256":h,
           "raw_component_os_exit_sha256":h,"peer_review_file":relative(peer),"peer_review_sha256":SOURCE_PINS[peer]}
    newo={"status":FRESH_OS,"pid":8,"pilot_id":p,"unified_exec_session":25,"observed_exit_code":0,
          "peer_review_sha256":SOURCE_PINS[peer],"launch_command_sha256":h,"receipt_sha256":h}
    newn={**phase,"pilot_id":p,"status":base.NATIVE_STATUS,"radius_exact":str(base.RADIUS),"checks":dict.fromkeys(V3_CHECKS,True),
          "remaining_noncontained":[],"fresh_surface_H_jobs":24,"new_bulk_force_or_Hessian_moment_jobs":0,
          "original_raw_component_receipt_file":relative(RAW/p/"receipt.json")}
    newt={**phase,"status":"FINISHED","exit_code":0,"receipt_sha256":h}
    news={**phase,"status":"RUNNING"}
    v3_status(newo,newn,newt,news,p,h,o)
    row={"kind":"surface","slot":2,"global_indices":[17],"point_force_scope":base.POINT,"whole_box_H_scope":base.BOX,
         "plan_file_sha256":h,"packet_file":"fixture","packet_sha256":h,"replay_file":"raw","replay_sha256":h}
    old={**row,"checked_force_rows":1,"checked_upper_H_entries":1,"checked_cache_sites":12,
         "containment_failures":{"direct_point_cache":[],"point_force":[0],"whole_box_H":[[0,0]]},
         "independent_point_F_exact":[["-2","2"]],"independent_whole_box_H_exact":[[0,0,"-4","4"]]}
    new={**old,"raw_component_replay_file":"raw","raw_component_replay_sha256":h,
         "containment_failures":dict.fromkeys(base.FAILURE_KEYS,[]),
         "independent_point_F_exact":[["-1","1"]],"independent_whole_box_H_exact":[[0,0,"-1","1"]],
         "point_original_excess_rounded_tails":[[0,"3/2","1","1/2",True]],
         "fresh_surface_H_tails":[[0,0,"1/2"]]}
    actual={"point_force":{FKEY:[["-1","1"]],"force_tail_upper_exact":["1"]},
            "whole_box_hessian":{HKEY:[[0,0,"-2","2"]]}}
    def leaf(r=old,n=new):
        return leaf_refinement(row,row,r,n,actual,(Q(1),Q(1)))
    require(leaf()==(1,1),"valid exact Minkowski truncation must pass")
    # For bulk no raw H miss can be bypassed by a companion H or tightened scalar.
    bulk_row={**row,"kind":"bulk"}
    bulk_old={**old,"kind":"bulk","checked_cache_sites":64}
    bulk_new={**new,"kind":"bulk","checked_cache_sites":64,"point_original_excess_rounded_tails":[[0,"2","1","1",True]]}
    return {"raw_host1_requires_explicit_declared_HOLD2_and_full_source_phase":(
                base.rejects(lambda:raw_status({**o,"observed_exit_code":0},n,t,s,p,h))
                and base.rejects(lambda:raw_status({**o,"native_documented_exit_code":1},n,t,s,p,h))),
            "reject_v3_failed_nonliteral_or_mixed_scope_exit":(
                base.rejects(lambda:v3_status({**newo,"observed_exit_code":False},newn,newt,news,p,h,o))
                and base.rejects(lambda:v3_status(newo,{**newn,"radius_exact":str(2*base.RADIUS)},newt,news,p,h,o))
                and base.rejects(lambda:v3_status(newo,{**newn,"fresh_surface_H_jobs":23},newt,news,p,h,o))),
            "exact_point_tail_erosion_and_strict_containment_rechecked":leaf()==(1,1),
            "reject_fake_exact_endpoints_or_overlapping_H":(
                base.rejects(lambda:leaf(n={**new,"independent_point_F_exact":[["-1/2","1/2"]]}))
                and base.rejects(lambda:leaf(n={**new,"independent_whole_box_H_exact":[[0,0,"-3","1"]]}))),
            "reject_direct_cache_miss_and_changed_projected_coordinates":(
                base.rejects(lambda:leaf(r={**old,"containment_failures":{**old["containment_failures"],"direct_point_cache":[0]}}))
                and base.rejects(lambda:leaf(n={**new,"global_indices":[18]}))),
            "reject_new_bulk_H_or_any_original_bulk_H_miss":base.rejects(
                lambda:leaf_refinement(bulk_row,bulk_row,bulk_old,bulk_new,actual,(Q(1),Q(1)))),
            "reject_whole_H_tail_typing_or_cap_violations":(
                base.rejects(lambda:leaf(n={**new,"fresh_surface_H_tails":[[False,0,"1/2"]]}))
                and base.rejects(lambda:leaf(n={**new,"fresh_surface_H_tails":[[0,0,"2"]]}))),
            "same_old_pure_kernel_and_all168_order_guards":all(base.fixtures().values())}


def prepare():
    PREP.mkdir(exist_ok=False)
    with (PREP/"transfer_ambiguity_exact_bridge_reduce_driver_v3_frozen.py").open("xb") as f:
        f.write(SOURCE.read_bytes())
    try:
        p=protocol()
        checks=fixtures()
        require(base.yes(checks) and not ACTUAL.exists(),"cold two-proof preflight failed or actual already exists")
        ph=save(PREP/"protocol.json",p)
        pf=save(PREP/"preflight.json",{"status":"PASS","source_sha256":sha(SOURCE),"protocol_sha256":ph,
              "checks":checks,"actual_packet_payloads_opened":0,"actual_matrix_work":0,
              "actual_root_claimed":False})
        print(json.dumps({"source_sha256":sha(SOURCE),"protocol_sha256":ph,"preflight_sha256":pf}))
    except BaseException as error:
        save(PREP/"failure.json",{"error":str(error),"traceback":traceback.format_exc()})
        raise


if __name__=="__main__":
    cli=argparse.ArgumentParser()
    sub=cli.add_subparsers(dest="stage",required=True)
    sub.add_parser("prepare")
    actual=sub.add_parser("run")
    actual.add_argument("--pilot-id",choices=base.PILOTS,required=True)
    for n in ("review","admission","raw-observation","companion-observation"):
        actual.add_argument("--"+n,required=True,type=lambda s:Path(s).resolve())
        actual.add_argument("--"+n+"-sha256",required=True)
    args=cli.parse_args()
    if args.stage=="prepare":
        prepare()
    else:
        run(args)

